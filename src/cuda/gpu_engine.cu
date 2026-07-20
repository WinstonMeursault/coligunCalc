#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "coilgun/simulation/cuda/gpu_mutual_pipeline.hpp"
#include "coilgun/simulation/cuda/gpu_state_kernels.hpp"
#include "coilgun/simulation/cuda/gpu_thermal.hpp"

#include <cuda_runtime_api.h>

#include <chrono>
#include <cmath>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <numeric>
#include <atomic>
#include <thread>
#include <sstream>

namespace coilgun::simulation::cuda {

namespace {

class ScopedCudaDevice {
public:
    explicit ScopedCudaDevice(int target) : target_(target) {
        const auto status = cudaGetDevice(&previous_);
        if (status != cudaSuccess)
            throw std::runtime_error(std::string("cudaGetDevice: ") + cudaGetErrorString(status));
        if (previous_ != target_) {
            const auto set_status = cudaSetDevice(target_);
            if (set_status != cudaSuccess)
                throw std::runtime_error(std::string("cudaSetDevice: ") + cudaGetErrorString(set_status));
            switched_ = true;
        }
    }

    ~ScopedCudaDevice() noexcept {
        if (switched_) (void)cudaSetDevice(previous_);
    }

    ScopedCudaDevice(const ScopedCudaDevice&) = delete;
    ScopedCudaDevice& operator=(const ScopedCudaDevice&) = delete;

private:
    int target_ = -1;
    int previous_ = -1;
    bool switched_ = false;
};

} // namespace

GpuEngine::Resources::~Resources() {
    int previous = -1;
    if (device_id >= 0 && cudaGetDevice(&previous) == cudaSuccess)
        (void)cudaSetDevice(device_id);
    thermal_workspace.reset();
    for (void* buffer : buffers) cudaFree(buffer);
    if (previous >= 0 && previous != device_id) (void)cudaSetDevice(previous);
}

void GpuEngine::shutdown_runtime() noexcept {
    release_runtime_resources();
}

void GpuEngine::release_runtime_resources() noexcept {
    int previous_device = -1;
    const bool restore_device = config_.device_id >= 0 &&
        cudaGetDevice(&previous_device) == cudaSuccess &&
        cudaSetDevice(config_.device_id) == cudaSuccess;
    if (context_) {
        try { context_->synchronize(); } catch (...) { }
    }
    if (persistent_active_) {
        free_persistent_buffers(persistent_buffers_);
        persistent_active_ = false;
    }
    persistent_adaptor_ = GpuAdaptor{};
    graph_cache_.reset();
    // GpuSolver may retain a non-owning context pointer, so destroy it before
    // the context itself. The replacement solver is created below on CPU.
    solver_.reset();
    context_.reset();
    if (resources_) {
        resources_->thermal_workspace.reset();
        for (void* buffer : resources_->buffers) (void)cudaFree(buffer);
        resources_->buffers.clear();
        resources_->device_id = -1;
    }
    if (restore_device && previous_device != config_.device_id)
        (void)cudaSetDevice(previous_device);
}

std::vector<std::uintptr_t> GpuEngine::device_buffer_addresses() const {
    if (!resources_) return {};
    std::vector<std::uintptr_t> result;
    result.reserve(resources_->buffers.size());
    for (void* buffer : resources_->buffers)
        result.push_back(reinterpret_cast<std::uintptr_t>(buffer));
    if (resources_->thermal_workspace) {
        const auto thermal = resources_->thermal_workspace->device_addresses();
        result.insert(result.end(), thermal.begin(), thermal.end());
    }
    return result;
}

bool cuda_device_available() noexcept {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

std::unique_ptr<GpuExecutionContext> make_gpu_execution_context() {
    return std::make_unique<GpuExecutionContext>();
}

void GpuEngine::initialize_runtime() {
    config_.validate();
#if defined(COILGUN_CUDA_AVAILABLE)
    // Auto is a host-only conservative fallback until the runtime confirms a
    // usable device. Once confirmed, use the direct physical pipeline; the
    // wrapper's explicit Graph request remains the graph path.
    if (config_.backend == BackendMode::Auto) {
        policy_.backend = BackendMode::Direct;
        policy_.backend_fallback_reason = FallbackReason::None;
        report_.backend = BackendMode::Direct;
        report_.static_fallback_reason = FallbackReason::None;
        report_.fallback_reason.clear();
    }
#endif
    // Runtime ownership follows the resolved policy. A requested Graph that
    // failed capability planning is already a CPU contract and must not even
    // create a CUDA context.
    const bool needs_cuda = policy_.backend != BackendMode::Fallback;
    if (!needs_cuda) {
        policy_.solver = SolverMode::Eigen;
        report_.solver = SolverMode::Eigen;
        if (policy_.thermal == ThermalMode::Gpu) {
            policy_.thermal = ThermalMode::Cpu;
            report_.thermal = ThermalMode::Cpu;
        }
        solver_ = std::make_unique<GpuSolver>(SolverMode::Eigen,
                                               SolverBatchLayout{layout_.B, layout_.D});
        const auto status = solver_->initialize_workspace();
        if (!status.ok) throw std::runtime_error(status.message);
        if (config_.backend == BackendMode::Fallback && report_.fallback_reason.empty())
            report_.fallback_reason = "CPU fallback explicitly requested";
        select_graph_variant_at_boundary();
        return;
    }
    if (fault_injection_.fail_device_initialization) {
        record_runtime_failure(SolverStatus::failure_status(
            SolverFailure::InvalidArgument, "injected CUDA device initialization failure"));
        return;
    }
    int device_count = 0;
    const auto count_status = cudaGetDeviceCount(&device_count);
    const bool unavailable = count_status == cudaErrorNoDevice ||
                            count_status == cudaErrorInsufficientDriver;
    if (count_status != cudaSuccess && !unavailable) {
        record_runtime_failure(SolverStatus::failure_status(
            SolverFailure::InvalidArgument,
            std::string("CUDA device enumeration failed: ") + cudaGetErrorString(count_status)));
        return;
    }
    if (unavailable || device_count == 0) {
        report_.runtime_fallback_reason = FallbackReason::CapabilityUnavailable;
        report_.fallback_reason = "CUDA device unavailable; using CPU fallback";
        ++report_.fallback_count;
        fallback_locked_ = true;
        policy_.backend = BackendMode::Fallback;
        policy_.solver = SolverMode::Eigen;
        normalize_cpu_thermal();
        report_.backend = BackendMode::Fallback;
        report_.solver = SolverMode::Eigen;
        solver_ = std::make_unique<GpuSolver>(SolverMode::Eigen,
                                               SolverBatchLayout{layout_.B, layout_.D});
        const auto status = solver_->initialize_workspace();
        if (!status.ok) throw std::runtime_error(status.message);
        return;
    }
    if (config_.device_id >= device_count) {
        std::ostringstream message;
        message << "CUDA device_id " << config_.device_id << " is unavailable (" << device_count
                << " device(s) detected)";
        throw std::invalid_argument(message.str());
    }
    ScopedCudaDevice device(config_.device_id);
    int selected_device = -1;
    const auto selected_status = cudaGetDevice(&selected_device);
    if (selected_status != cudaSuccess || selected_device != config_.device_id) {
        record_runtime_failure(SolverStatus::failure_status(
            SolverFailure::InvalidArgument,
            std::string("configured CUDA device was not selected: ") +
                cudaGetErrorString(selected_status)));
        return;
    }
    resources_->device_id = config_.device_id;
    if (policy_.backend == BackendMode::Persistent) {
        const auto flags = cudaSetDeviceFlags(cudaDeviceMapHost);
        if (flags != cudaSuccess && flags != cudaErrorSetOnActiveProcess) {
            record_runtime_failure(SolverStatus::failure_status(
                SolverFailure::InvalidArgument, "CUDA host mapping is unavailable"));
            return;
        }
    }
    if (needs_cuda && cuda_device_available()) {
        try {
            context_ = std::make_unique<GpuExecutionContext>(GpuExecutionContextConfig{
                config_.device_id, cudaStreamNonBlocking, 0, config_.enable_profiling});
        } catch (const std::exception& error) {
            record_runtime_failure(SolverStatus::failure_status(
                SolverFailure::InvalidArgument,
                std::string("execution context initialization failed: ") + error.what()));
        }
    } else if (needs_cuda) {
        policy_.backend = BackendMode::Fallback;
        policy_.solver = SolverMode::Eigen;
        report_.backend = BackendMode::Fallback;
        report_.solver = SolverMode::Eigen;
        normalize_cpu_thermal();
        report_.runtime_fallback_reason = FallbackReason::CapabilityUnavailable;
        report_.fallback_reason = "CUDA device unavailable; using CPU fallback";
        ++report_.fallback_count;
        fallback_locked_ = true;
    }

    if (policy_.solver == SolverMode::Batched && !context_available()) {
        report_.solver = SolverMode::Eigen;
        report_.runtime_fallback_reason = report_.runtime_fallback_reason == FallbackReason::None
            ? FallbackReason::CapabilityUnavailable : report_.runtime_fallback_reason;
        if (report_.fallback_reason.empty()) report_.fallback_reason = "batched solver requires a CUDA execution context";
        ++report_.fallback_count;
        fallback_locked_ = true;
        solver_ = std::make_unique<GpuSolver>(SolverMode::Eigen,
                                              SolverBatchLayout{layout_.B, layout_.D});
    } else if (context_available()) {
        solver_ = std::make_unique<GpuSolver>(*context_, policy_,
                                              SolverBatchLayout{layout_.B, layout_.D});
    } else {
        solver_ = std::make_unique<GpuSolver>(SolverMode::Eigen,
                                              SolverBatchLayout{layout_.B, layout_.D});
    }

    if (context_available()) {
        const auto allocate = [&](std::size_t bytes) {
            void* pointer = nullptr;
            if (fault_injection_.fail_allocation && device_allocation_count_ >= 3)
                throw std::runtime_error("injected CUDA allocation failure");
            if (cudaMalloc(&pointer, bytes) != cudaSuccess)
                throw std::runtime_error("persistent pipeline CUDA allocation failed");
            try {
                resources_->buffers.push_back(pointer);
            } catch (...) {
                (void)cudaFree(pointer);
                throw;
            }
            ++device_allocation_count_;
            return pointer;
        };
        const auto B = layout_.B;
        const auto S = layout_.S;
        const auto F = layout_.F;
        try {
            allocate(sizeof(CoilGeo) * S);
            allocate(sizeof(FilGeo) * F);
            allocate(sizeof(double) * B * S * F);
            allocate(sizeof(std::uint8_t) * B * S * F);
            allocate(sizeof(std::uint8_t) * B);
            allocate(sizeof(std::uint8_t) * B * S);
            allocate(sizeof(double) * B * S * F);
            allocate(sizeof(double) * B * S * F);
            allocate(sizeof(double) * B * (S + F));
            allocate(sizeof(double) * B * (S + F));
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * B);
            if (policy_.thermal != ThermalMode::Disabled) {
                material_tables_ = generate_material_tables();
            }
            if (policy_.thermal == ThermalMode::Gpu) {
                resources_->thermal_workspace = std::make_unique<ThermalWorkspace>(
                    material_tables_, B * F);
            }

            if (policy_.backend == BackendMode::Graph) {
                graph_cache_ = std::make_unique<GpuGraphCache>();
            } else if (policy_.backend == BackendMode::Persistent) {
                // The current engine's synchronous CUDA stream is not a safe
                // host for a resident polling kernel. Keep the protocol object
                // available for its dedicated tests, but resolve this engine
                // path to fallback until a separate control stream is owned.
                throw std::runtime_error("persistent protocol requires a dedicated control stream");
            }
        } catch (const std::exception& error) {
            record_runtime_failure(SolverStatus::failure_status(
                SolverFailure::InvalidArgument,
                std::string("CUDA runtime resource initialization failed: ") + error.what()));
        }
    }

    const SolverStatus workspace = solver_->initialize_workspace();
    if (!workspace.ok) {
        record_runtime_failure(workspace);
    }
    if (config_.enable_calibration && !fallback_locked_) {
        calibrate_solver();
    }
    select_graph_variant_at_boundary();
}

void GpuEngine::record_runtime_failure(const SolverStatus& status) {
    fallback_locked_ = true;
    policy_.backend = BackendMode::Fallback;
    policy_.solver = SolverMode::Eigen;
    report_.backend = BackendMode::Fallback;
    report_.solver = SolverMode::Eigen;
    report_.runtime_fallback_reason = FallbackReason::RuntimeFailure;
    report_.fallback_reason = status.message;
    ++report_.fallback_count;
    if (policy_.thermal == ThermalMode::Gpu) {
        policy_.thermal = ThermalMode::Cpu;
        report_.thermal = ThermalMode::Cpu;
    }
    release_runtime_resources();
    solver_ = std::make_unique<GpuSolver>(SolverMode::Eigen,
                                           SolverBatchLayout{layout_.B, layout_.D});
    const SolverStatus fallback = solver_->initialize_workspace();
    if (!fallback.ok) {
        report_.fallback_reason += "; CPU fallback initialization failed: " + fallback.message;
        throw std::runtime_error(report_.fallback_reason);
    }
    report_.solver = SolverMode::Eigen;
    report_.backend = BackendMode::Fallback;
    select_graph_variant_at_boundary();
}

void GpuEngine::execute_solver_step() {
    if (!solver_) return;

    const auto B = layout_.B;
    const auto D = layout_.D;
    if (fallback_locked_ && solver_ && solver_->resolved_mode() != SolverMode::Eigen) {
        solver_ = std::make_unique<GpuSolver>(SolverMode::Eigen, SolverBatchLayout{B, D});
        (void)solver_->initialize_workspace();
        report_.solver = SolverMode::Eigen;
    }

    const auto start = std::chrono::steady_clock::now();
    const SolverStatus status = solver_->solve_batch(matrices_.data(), rhs_.data(), solution_.data());
    const auto stop = std::chrono::steady_clock::now();
    report_.solver_time_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    current_derivatives_ = solution_;
    state_.current_derivatives = solution_;
    if (!status.ok) {
        record_runtime_failure(status);
        throw std::runtime_error("solver failed; current step was not committed: " + status.message);
    }
}

void GpuEngine::execute_physical_pipeline() {
    if (fallback_locked_ || !context_available()) {
        execute_cpu_physical_pipeline();
        return;
    }
    // Host wall time for a successful CUDA-backed physical pipeline, including
    // host orchestration and synchronous transfers; not device-only time.
    const auto gpu_start = std::chrono::steady_clock::now();
    ScopedCudaDevice device(config_.device_id);
    pipeline_order_.clear();
    pipeline_order_.push_back(PipelineStage::Mutual);

    // The host representation is deliberately kept canonical.  The CUDA path
    // mirrors it into transient launch buffers and copies only the physical
    // results back at the synchronous boundary.
    const std::size_t B = layout_.B;
    const std::size_t S = layout_.S;
    const std::size_t F = layout_.F;
    const auto state_snapshot = state_;
    bool solver_done = false;
    std::vector<double> separations(B * S * F, 0.0);
    std::vector<std::uint8_t> active(B * S * F, 0);
    std::vector<CoilGeo> coils(S);
    std::vector<FilGeo> filaments(F);
    for (std::size_t s = 0; s < S; ++s) {
        coils[s] = CoilGeo{geometry_.stage_inner_radii[s], geometry_.stage_outer_radii[s],
                           geometry_.stage_lengths[s], geometry_.stage_positions[s],
                           geometry_.stage_turns[s]};
    }
    for (std::size_t f = 0; f < F; ++f) {
        filaments[f] = FilGeo{geometry_.filament_inner_radii[f], geometry_.filament_outer_radii[f],
                              geometry_.filament_lengths[f]};
    }
    for (std::size_t b = 0; b < B; ++b) {
        for (std::size_t s = 0; s < S; ++s) {
            for (std::size_t f = 0; f < F; ++f) {
                const auto index = mutual_pipeline_index(b, s, f, S, F);
                active[index] = state_.active_mask[b] != 0 && mutual_stage_mask_[b * S + s] != 0;
                separations[index] = geometry_.stage_positions[s] -
                    (geometry_.filament_positions[f] - state_.position[b]);
            }
        }
    }

#if defined(COILGUN_CUDA_AVAILABLE)
    if (context_available()) {
        auto* d_coils = static_cast<CoilGeo*>(resources_->buffers[0]);
        auto* d_filaments = static_cast<FilGeo*>(resources_->buffers[1]);
        auto* d_separations = static_cast<double*>(resources_->buffers[2]);
        auto* d_active = static_cast<std::uint8_t*>(resources_->buffers[3]);
        auto* d_batch_active = static_cast<std::uint8_t*>(resources_->buffers[4]);
        auto* d_trigger = static_cast<std::uint8_t*>(resources_->buffers[5]);
        auto* d_mutual = static_cast<double*>(resources_->buffers[6]);
        auto* d_gradient = static_cast<double*>(resources_->buffers[7]);
        auto* d_currents = static_cast<double*>(resources_->buffers[8]);
        auto* d_derivative = static_cast<double*>(resources_->buffers[9]);
        auto* d_acceleration = static_cast<double*>(resources_->buffers[10]);
        auto* d_velocity = static_cast<double*>(resources_->buffers[11]);
        auto* d_position = static_cast<double*>(resources_->buffers[12]);
        auto* d_force = static_cast<double*>(resources_->buffers[13]);
        try {
            const auto copy = [&](void* destination, const void* source, std::size_t bytes,
                                  cudaMemcpyKind kind, const char* operation) {
                const auto transfer_start = std::chrono::steady_clock::now();
                const auto status = cudaMemcpy(destination, source, bytes, kind);
                if (status != cudaSuccess)
                    throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
                const auto transfer_stop = std::chrono::steady_clock::now();
                report_.transfer_time_ms +=
                    std::chrono::duration<double, std::milli>(transfer_stop - transfer_start).count();
            };
            std::vector<std::uint8_t> effective_trigger = state_.trigger_mask;
            for (std::size_t b = 0; b < B; ++b) {
                for (std::size_t s = 0; s < S; ++s) {
                    if (state_.active_mask[b] == 0 || mutual_stage_mask_[b * S + s] == 0)
                        effective_trigger[b * S + s] = 0;
                }
            }
            copy(d_coils, coils.data(), sizeof(CoilGeo) * S, cudaMemcpyHostToDevice, "coil upload");
            copy(d_filaments, filaments.data(), sizeof(FilGeo) * F, cudaMemcpyHostToDevice, "filament upload");
            copy(d_separations, separations.data(), sizeof(double) * separations.size(), cudaMemcpyHostToDevice, "separation upload");
            copy(d_active, active.data(), sizeof(std::uint8_t) * active.size(), cudaMemcpyHostToDevice, "active upload");
            copy(d_batch_active, state_.active_mask.data(), sizeof(std::uint8_t) * B, cudaMemcpyHostToDevice, "batch active upload");
            copy(d_trigger, effective_trigger.data(), sizeof(std::uint8_t) * state_.trigger_mask.size(), cudaMemcpyHostToDevice, "trigger upload");
            copy(d_currents, state_.currents.data(), sizeof(double) * state_.currents.size(), cudaMemcpyHostToDevice, "current upload");
            copy(d_velocity, state_.velocity.data(), sizeof(double) * B, cudaMemcpyHostToDevice, "velocity upload");
            copy(d_position, state_.position.data(), sizeof(double) * B, cudaMemcpyHostToDevice, "position upload");
            const auto launch_mutual = [&](cudaStream_t stream) {
                launch_mutual_pipeline(MutualPipelineView{d_coils, d_filaments, d_separations, d_active,
                                                           d_mutual, d_gradient, B, S, F, 9},
                                         policy_.precision == PrecisionMode::Aggressive
                                             ? GpuOptLevel::Aggressive
                                             : (policy_.precision == PrecisionMode::Full
                                                    ? GpuOptLevel::Full : GpuOptLevel::Standard),
                                          config_.threads_per_block, stream);
                return GraphCaptureStatus::success();
             };
              if (policy_.backend == BackendMode::Graph && graph_cache_) {
                  if (fault_injection_.fail_graph_capture) {
                      fault_injection_.fail_graph_capture = false;
                      throw std::runtime_error("injected CUDA graph capture failure");
                  }
                 const auto key = graph_key();
                 if (!graph_cache_->has_current() || !(graph_cache_->current_key() == key)) {
                     // Upload the constant quadrature tables before capture. CUDA
                     // graph capture does not reliably retain symbol copies.
                     (void)launch_mutual(context_->stream());
                     context_->synchronize();
                 }
                  const auto capture_count = graph_cache_->capture_count();
                  auto selected = graph_cache_->capture_and_select(key, context_->stream(), launch_mutual);
                 if (!selected.ok) {
                    record_runtime_failure(SolverStatus::failure_status(
                        SolverFailure::InvalidArgument,
                        "CUDA graph capture failed: " + selected.failure.message));
                     throw std::runtime_error(report_.fallback_reason);
                 }
                 if (graph_cache_->capture_count() > capture_count)
                     ++report_.graph_rebuild_count;
                auto replayed = graph_cache_->replay(context_->stream());
                if (!replayed.ok) {
                    record_runtime_failure(SolverStatus::failure_status(
                        SolverFailure::InvalidArgument,
                        "CUDA graph replay failed: " + replayed.failure.message));
                    throw std::runtime_error(report_.fallback_reason);
                }
            } else if (persistent_active_) {
                execute_persistent_mutual();
            } else {
                (void)launch_mutual(context_->stream());
            }
            if (!persistent_active_)
                context_->synchronize();
             copy(state_.m1.data(), d_mutual, sizeof(double) * state_.m1.size(), cudaMemcpyDeviceToHost, "mutual download");
             copy(state_.dm1.data(), d_gradient, sizeof(double) * state_.dm1.size(), cudaMemcpyDeviceToHost, "gradient download");
              if (fault_injection_.fail_after_mutual) {
                  fault_injection_.fail_after_mutual = false;
                 throw std::runtime_error("injected GPU failure after mutual segment");
             }
             assemble_physical_system(true);
            pipeline_order_.push_back(PipelineStage::Matrix);
            execute_solver_step();
            pipeline_order_.push_back(PipelineStage::Solver);
            solver_done = true;
            copy(d_derivative, solution_.data(), sizeof(double) * solution_.size(), cudaMemcpyHostToDevice, "derivative upload");
            pipeline_order_.push_back(PipelineStage::Force);
            if (launch_state_update_masked(
                    B, S, F, d_currents, d_derivative, d_gradient,
                    d_trigger, d_batch_active, state_.mass, state_.dt,
                    d_acceleration, d_velocity, d_position, d_force,
                    StateKernelConfig{config_.deterministic,
                                      static_cast<unsigned int>(config_.threads_per_block)},
                    context_->stream()) != cudaSuccess) {
                throw std::runtime_error("state kernel launch failed");
            }
            context_->synchronize();
             copy(state_.currents.data(), d_currents, sizeof(double) * state_.currents.size(), cudaMemcpyDeviceToHost, "current download");
             copy(state_.velocity.data(), d_velocity, sizeof(double) * B, cudaMemcpyDeviceToHost, "velocity download");
             copy(state_.position.data(), d_position, sizeof(double) * B, cudaMemcpyDeviceToHost, "position download");
             pipeline_order_.push_back(PipelineStage::State);
             if (policy_.thermal != ThermalMode::Disabled && !state_.temperatures.empty()) {
                 const auto thermal_start = std::chrono::steady_clock::now();
                 std::vector<double> filament_currents(B * F, 0.0);
                for (std::size_t b = 0; b < B; ++b)
                    if (state_.active_mask[b] != 0)
                        std::copy_n(state_.currents.data() + b * layout_.D + S, F,
                                    filament_currents.data() + b * F);
                state_.resistivities.resize(B * F);
                state_.resistances.resize(B * F);
                state_.joule_energy.resize(B * F);
                const auto precision = policy_.precision == PrecisionMode::Aggressive
                    ? ThermalPrecision::Aggressive
                    : (policy_.precision == PrecisionMode::Full ? ThermalPrecision::Full : ThermalPrecision::Standard);
                if (policy_.thermal == ThermalMode::Cpu) {
                    update_thermal_batch_cpu(material_tables_, precision, B, F, filament_currents.data(),
                        state_.filament_masses.data(), state_.reference_resistances.data(), state_.filament_materials.data(),
                        state_.dt, state_.temperatures.data(), state_.resistivities.data(), state_.resistances.data(), state_.joule_energy.data());
                 } else {
                     resources_->thermal_workspace->update(material_tables_, precision, B, F, filament_currents.data(),
                         state_.filament_masses.data(), state_.reference_resistances.data(), state_.filament_materials.data(),
                         state_.dt, state_.temperatures.data(), state_.resistivities.data(), state_.resistances.data(), state_.joule_energy.data(), context_->stream());
                  }
                  context_->synchronize();
                 const auto thermal_stop = std::chrono::steady_clock::now();
                 report_.thermal_time_ms +=
                     std::chrono::duration<double, std::milli>(thermal_stop - thermal_start).count();
                pipeline_order_.push_back(PipelineStage::Thermal);
             }
             restore_inactive_state(state_snapshot);
         } catch (...) { throw; }
    }
#endif

    if (!solver_done) {
        pipeline_order_.push_back(PipelineStage::Matrix);
        execute_solver_step();
        pipeline_order_.push_back(PipelineStage::Solver);
    }
    if (!solver_done) {
        pipeline_order_.push_back(PipelineStage::Force);
    }
    if (policy_.thermal != ThermalMode::Disabled && !solver_done) {
        pipeline_order_.push_back(PipelineStage::Thermal);
        if (!state_.temperatures.empty()) {
            std::vector<double> filament_currents(B * F, 0.0);
            for (std::size_t b = 0; b < B; ++b)
                std::copy_n(state_.currents.data() + b * layout_.D + S,
                            F, filament_currents.data() + b * F);
            state_.resistivities.resize(B * F);
            state_.resistances.resize(B * F);
            state_.joule_energy.resize(B * F);
            const auto& tables = material_tables_;
            const auto precision = policy_.precision == PrecisionMode::Aggressive
                ? ThermalPrecision::Aggressive
                : (policy_.precision == PrecisionMode::Full ? ThermalPrecision::Full
                                                              : ThermalPrecision::Standard);
            if (policy_.thermal == ThermalMode::Cpu || !context_available()) {
                update_thermal_batch_cpu(tables, precision, B, F, filament_currents.data(),
                                 state_.filament_masses.data(), state_.reference_resistances.data(),
                                 state_.filament_materials.data(), state_.dt,
                                 state_.temperatures.data(), state_.resistivities.data(),
                                 state_.resistances.data(), state_.joule_energy.data());
             } else resources_->thermal_workspace->update(tables, precision, B, F, filament_currents.data(),
                                 state_.filament_masses.data(), state_.reference_resistances.data(),
                                 state_.filament_materials.data(), state_.dt,
                                 state_.temperatures.data(), state_.resistivities.data(),
                                 state_.resistances.data(), state_.joule_energy.data(), context_->stream());
        }
    }
     pipeline_order_.push_back(PipelineStage::State);
     const auto gpu_stop = std::chrono::steady_clock::now();
     report_.gpu_time_ms +=
         std::chrono::duration<double, std::milli>(gpu_stop - gpu_start).count();
}

GpuGraphVariantKey GpuEngine::graph_key() const {
    GpuGraphVariantKey key;
    key.batch_capacity = layout_.B;
    key.layout_signature = static_cast<std::uint64_t>(layout_.D) * 1099511628211ull ^ layout_.S ^ (layout_.F << 16);
    key.precision = policy_.precision;
    key.thermal = policy_.thermal;
    key.solver = policy_.solver;
    for (const auto bit : stage_mask_)
        key.stage_signature = key.stage_signature * 1315423911u + bit + 1;
    for (const auto bit : mutual_stage_mask_)
        key.stage_signature = key.stage_signature * 1315423911u + bit + 1;
    return key;
}

void GpuEngine::execute_persistent_mutual() {
    if (layout_.B != 1)
        throw std::logic_error("persistent mutual protocol requires a single-batch engine");
    const std::size_t pairs = layout_.S * layout_.F;
    auto* generation = reinterpret_cast<volatile int*>(persistent_buffers_.generation);
    auto* doorbell = reinterpret_cast<volatile int*>(persistent_buffers_.doorbell);
    auto* out_m = reinterpret_cast<volatile double*>(persistent_buffers_.out_M);
    auto* out_dm = reinterpret_cast<volatile double*>(persistent_buffers_.out_dM);
    for (std::size_t b = 0; b < layout_.B; ++b) {
        if (b != 0) break;
        for (std::size_t s = 0; s < layout_.S; ++s) {
            for (std::size_t f = 0; f < layout_.F; ++f) {
                const auto index = s * layout_.F + f;
                persistent_buffers_.seps[index] = geometry_.stage_positions[s] -
                    (geometry_.filament_positions[f] - state_.position[b]);
                doorbell[index] = 1;
            }
        }
        std::atomic_thread_fence(std::memory_order_seq_cst);
        ++*generation;
        std::atomic_thread_fence(std::memory_order_seq_cst);
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(250);
        for (std::size_t index = 0; index < pairs; ++index) {
            while (doorbell[index] != 0) {
                if (std::chrono::steady_clock::now() >= deadline)
                    throw std::runtime_error("persistent protocol timed out");
                if (cudaPeekAtLastError() != cudaSuccess)
                    throw std::runtime_error("persistent protocol execution failed");
                std::atomic_thread_fence(std::memory_order_seq_cst);
                std::this_thread::yield();
            }
        }
        for (std::size_t s = 0; s < layout_.S; ++s) {
            for (std::size_t f = 0; f < layout_.F; ++f) {
                const auto pair = s * layout_.F + f;
                state_.m1[pair] = out_m[pair];
                state_.dm1[pair] = out_dm[pair];
            }
        }
    }
    (void)pairs;
}

void GpuEngine::calibrate_solver() {
    const auto B = layout_.B;
    const auto D = layout_.D;
    std::fill(matrices_.begin(), matrices_.end(), 0.0);
    std::fill(rhs_.begin(), rhs_.end(), 0.0);
    std::fill(solution_.begin(), solution_.end(), 0.0);
    for (std::size_t batch = 0; batch < B; ++batch) {
        for (std::size_t index = 0; index < D; ++index) {
            matrices_[batch * D * D + index * D + index] = 1.0;
            rhs_[batch * D + index] = state_.currents[batch * D + index];
        }
    }
    const auto start = std::chrono::steady_clock::now();
    const SolverStatus status = solver_->solve_batch(matrices_.data(), rhs_.data(), solution_.data());
    const auto stop = std::chrono::steady_clock::now();
    report_.solver_time_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    calibration_count_ = 1;
    if (!status.ok) {
        record_runtime_failure(status);
        return;
    }
    report_.calibrated = true;
    report_.max_condition_estimate = 1.0;
}

} // namespace coilgun::simulation::cuda
