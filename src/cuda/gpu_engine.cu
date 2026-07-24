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

enum DeviceBufferIndex : std::size_t {
    Coils = 0,
    Filaments,
    Separations,
    PairActive,
    BatchActive,
    TriggerMask,
    Mutual,
    Gradient,
    Currents,
    PreStepCurrents,
    Derivative,
    Acceleration,
    Velocity,
    Position,
    Force,
    StageInductance,
    StageResistance,
    StageMutual,
    FilamentInductance,
    FilamentReferenceResistance,
    FilamentMutual,
    StageMask,
    MutualStageMask,
    StageVoltage,
    DynamicResistance,
    Matrix,
    Rhs,
    Solution,
    Residual,
    CompactStatus,
    TriggerMode,
    TriggerValue,
    ExcitationFinished,
    StageCompleted,
    TriggerTime,
    TriggerPosition,
    PositionOffset,
    ControlTime,
    StagePosition,
    FilamentPosition,
    BufferCount,
};

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

void GpuEngine::sync_runtime_state_after_reset() {
    if (!context_available() || resources_->buffers.size() < BufferCount) return;
    ScopedCudaDevice device(config_.device_id);
    const auto B = layout_.B;
    const auto D = layout_.D;
    const auto copy = [](void* destination, const void* source, std::size_t bytes,
                         const char* operation) {
        const auto status = cudaMemcpy(destination, source, bytes, cudaMemcpyHostToDevice);
        if (status != cudaSuccess)
            throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
    };
    copy(resources_->buffers[Currents], state_.currents.data(), B * D * sizeof(double),
         "reset current upload");
    copy(resources_->buffers[Velocity], state_.velocity.data(), B * sizeof(double),
         "reset velocity upload");
    copy(resources_->buffers[Position], state_.position.data(), B * sizeof(double),
         "reset position upload");
    if (policy_.thermal == ThermalMode::Gpu && resources_->thermal_workspace) {
        resources_->thermal_workspace->initialize_device_state(
            material_tables_, B, layout_.F, state_.filament_masses.data(),
            state_.reference_resistances.data(), state_.filament_materials.data(),
            state_.temperatures.data(), state_.resistances.data(), context_->stream());
        context_->synchronize();
    }
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

GpuAssemblySnapshot GpuEngine::assemble_reference_for_test() {
    ensure_running();
    assemble_physical_system(false);
    return {matrices_, rhs_};
}

void GpuEngine::complete_stage(std::size_t batch, std::size_t stage) {
    ensure_running();
    if (batch >= layout_.B || stage >= layout_.S)
        throw std::out_of_range("completed stage index is outside GPU layout");
    state_.currents[batch * layout_.D + stage] = 0.0;
    state_.stage_mask[batch * layout_.S + stage] = 0;
    state_.mutual_stage_mask[batch * layout_.S + stage] = 0;
    stage_mask_[batch * layout_.S + stage] = 0;
    mutual_stage_mask_[batch * layout_.S + stage] = 0;
    select_graph_variant_at_boundary();
    if (!context_available() || resources_->buffers.size() < BufferCount) return;
    ScopedCudaDevice device(config_.device_id);
    const auto error = cudaMemsetAsync(
        static_cast<double*>(resources_->buffers[Currents]) + batch * layout_.D + stage,
        0, sizeof(double), context_->stream());
    if (error != cudaSuccess)
        throw std::runtime_error(std::string("completed stage current clear: ") +
                                 cudaGetErrorString(error));
    context_->synchronize();
}

GpuAssemblySnapshot GpuEngine::assemble_device_for_test() {
    ensure_running();
    if (!context_available())
        throw std::logic_error("device assembly requires an initialized CUDA context");
    ScopedCudaDevice device(config_.device_id);
    const auto B = layout_.B;
    const auto S = layout_.S;
    const auto F = layout_.F;
    const auto D = layout_.D;
    auto copy = [](void* destination, const void* source, std::size_t bytes,
                   cudaMemcpyKind kind, const char* operation) {
        const auto status = cudaMemcpy(destination, source, bytes, kind);
        if (status != cudaSuccess)
            throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
    };
    std::vector<double> voltages = state_.stage_voltages;
    if (voltages.empty()) voltages.assign(B * S, 0.0);
    copy(resources_->buffers[Currents], state_.currents.data(), B * D * sizeof(double),
         cudaMemcpyHostToDevice, "test current upload");
    copy(resources_->buffers[Velocity], state_.velocity.data(), B * sizeof(double),
         cudaMemcpyHostToDevice, "test velocity upload");
    copy(resources_->buffers[BatchActive], state_.active_mask.data(), B * sizeof(std::uint8_t),
         cudaMemcpyHostToDevice, "test active upload");
    copy(resources_->buffers[TriggerMask], state_.trigger_mask.data(), B * S * sizeof(std::uint8_t),
         cudaMemcpyHostToDevice, "test trigger upload");
    copy(resources_->buffers[StageMask], stage_mask_.data(), B * S * sizeof(std::uint8_t),
         cudaMemcpyHostToDevice, "test stage mask upload");
    copy(resources_->buffers[MutualStageMask], mutual_stage_mask_.data(), B * S * sizeof(std::uint8_t),
         cudaMemcpyHostToDevice, "test mutual mask upload");
    copy(resources_->buffers[StageVoltage], voltages.data(), B * S * sizeof(double),
         cudaMemcpyHostToDevice, "test voltage upload");
    copy(resources_->buffers[Mutual], state_.m1.data(), B * S * F * sizeof(double),
         cudaMemcpyHostToDevice, "test mutual upload");
    copy(resources_->buffers[Gradient], state_.dm1.data(), B * S * F * sizeof(double),
         cudaMemcpyHostToDevice, "test gradient upload");
    if (state_.resistances.size() == B * F)
        copy(resources_->buffers[DynamicResistance], state_.resistances.data(),
             B * F * sizeof(double), cudaMemcpyHostToDevice,
             "test dynamic resistance upload");
    const auto launch = launch_device_assembly(DeviceAssemblyView{
        B, S, F,
        static_cast<double*>(resources_->buffers[StageInductance]),
        static_cast<double*>(resources_->buffers[StageResistance]),
        static_cast<double*>(resources_->buffers[StageMutual]),
        static_cast<double*>(resources_->buffers[FilamentInductance]),
        static_cast<double*>(resources_->buffers[FilamentReferenceResistance]),
        static_cast<double*>(resources_->buffers[FilamentMutual]),
        state_.resistances.size() == B * F
            ? static_cast<double*>(resources_->buffers[DynamicResistance]) : nullptr,
        static_cast<double*>(resources_->buffers[Mutual]),
        static_cast<double*>(resources_->buffers[Gradient]),
        static_cast<double*>(resources_->buffers[Currents]),
        static_cast<double*>(resources_->buffers[Velocity]),
        static_cast<double*>(resources_->buffers[StageVoltage]),
        static_cast<std::uint8_t*>(resources_->buffers[BatchActive]),
        static_cast<std::uint8_t*>(resources_->buffers[TriggerMask]),
        static_cast<std::uint8_t*>(resources_->buffers[StageMask]),
        static_cast<std::uint8_t*>(resources_->buffers[MutualStageMask]),
        static_cast<double*>(resources_->buffers[Matrix]),
        static_cast<double*>(resources_->buffers[Rhs])}, context_->stream());
    if (launch != cudaSuccess)
        throw std::runtime_error(std::string("test device assembly: ") + cudaGetErrorString(launch));
    context_->synchronize();
    GpuAssemblySnapshot result;
    result.matrix.resize(B * D * D);
    result.rhs.resize(B * D);
    copy(result.matrix.data(), resources_->buffers[Matrix], result.matrix.size() * sizeof(double),
         cudaMemcpyDeviceToHost, "test matrix download");
    copy(result.rhs.data(), resources_->buffers[Rhs], result.rhs.size() * sizeof(double),
         cudaMemcpyDeviceToHost, "test rhs download");
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
            context_->ensure_quadrature9_loaded();
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
        if (policy_.backend == BackendMode::Graph && policy_.solver != SolverMode::Batched) {
            policy_.solver = SolverMode::Batched;
            report_.solver = SolverMode::Batched;
        }
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
            // Keep this allocation order synchronized with DeviceBufferIndex.
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
            allocate(sizeof(double) * B * (S + F));
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * B);
            allocate(sizeof(double) * S);
            allocate(sizeof(double) * S);
            allocate(sizeof(double) * S * S);
            allocate(sizeof(double) * F);
            allocate(sizeof(double) * F);
            allocate(sizeof(double) * F * F);
            allocate(sizeof(std::uint8_t) * B * S);
            allocate(sizeof(std::uint8_t) * B * S);
            allocate(sizeof(double) * B * S);
            allocate(sizeof(double) * B * F);
            allocate(sizeof(double) * B * (S + F) * (S + F));
            allocate(sizeof(double) * B * (S + F));
            allocate(sizeof(double) * B * (S + F));
            allocate(sizeof(double) * B);
            allocate(sizeof(DeviceStepStatus) * B);
            allocate(sizeof(std::uint8_t) * B * S);
            allocate(sizeof(double) * B * S);
            allocate(sizeof(std::uint8_t) * B * S);
            allocate(sizeof(std::uint8_t) * B * S);
            allocate(sizeof(double) * B * S);
            allocate(sizeof(double) * B * S);
            allocate(sizeof(double) * B);
            allocate(sizeof(double));
            allocate(sizeof(double) * S);
            allocate(sizeof(double) * F);

            auto copy_to_device = [](void* destination, const void* source,
                                     std::size_t bytes, const char* operation) {
                const auto status = cudaMemcpy(destination, source, bytes,
                                               cudaMemcpyHostToDevice);
                if (status != cudaSuccess)
                    throw std::runtime_error(std::string(operation) + ": " +
                                             cudaGetErrorString(status));
            };
            std::vector<CoilGeo> host_coils(S);
            std::vector<FilGeo> host_filaments(F);
            std::vector<double> stage_inductances(S);
            std::vector<double> stage_resistances(S);
            std::vector<double> filament_inductances(F);
            std::vector<double> filament_resistances(F);
            for (std::size_t stage = 0; stage < S; ++stage) {
                host_coils[stage] = {geometry_.stage_inner_radii[stage],
                    geometry_.stage_outer_radii[stage], geometry_.stage_lengths[stage],
                    geometry_.stage_positions[stage], geometry_.stage_turns[stage]};
                stage_inductances[stage] = stage_inductance(stage);
                stage_resistances[stage] = stage_resistance(stage);
            }
            for (std::size_t filament = 0; filament < F; ++filament) {
                host_filaments[filament] = {geometry_.filament_inner_radii[filament],
                    geometry_.filament_outer_radii[filament],
                    geometry_.filament_lengths[filament]};
                filament_inductances[filament] = filament_inductance(filament);
                filament_resistances[filament] = filament_resistance(filament);
            }
            copy_to_device(resources_->buffers[Coils], host_coils.data(),
                           sizeof(CoilGeo) * S, "immutable coil upload");
            copy_to_device(resources_->buffers[Filaments], host_filaments.data(),
                           sizeof(FilGeo) * F, "immutable filament upload");
            copy_to_device(resources_->buffers[StageInductance], stage_inductances.data(),
                           sizeof(double) * S, "stage inductance upload");
            copy_to_device(resources_->buffers[StageResistance], stage_resistances.data(),
                           sizeof(double) * S, "stage resistance upload");
            copy_to_device(resources_->buffers[StageMutual],
                           geometry_.stage_mutual_inductances.data(),
                           sizeof(double) * S * S, "stage mutual upload");
            copy_to_device(resources_->buffers[FilamentInductance],
                           filament_inductances.data(), sizeof(double) * F,
                           "filament inductance upload");
            copy_to_device(resources_->buffers[FilamentReferenceResistance],
                           filament_resistances.data(), sizeof(double) * F,
                           "filament resistance upload");
            copy_to_device(resources_->buffers[FilamentMutual],
                           geometry_.filament_mutual_inductances.data(),
                           sizeof(double) * F * F, "filament mutual upload");
            copy_to_device(resources_->buffers[StagePosition],
                           geometry_.stage_positions.data(), sizeof(double) * S,
                           "stage position upload");
            copy_to_device(resources_->buffers[FilamentPosition],
                           geometry_.filament_positions.data(), sizeof(double) * F,
                           "filament position upload");
            copy_to_device(resources_->buffers[Currents], state_.currents.data(),
                           sizeof(double) * B * (S + F), "initial current upload");
            copy_to_device(resources_->buffers[Velocity], state_.velocity.data(),
                           sizeof(double) * B, "initial velocity upload");
            copy_to_device(resources_->buffers[Position], state_.position.data(),
                           sizeof(double) * B, "initial position upload");
            if (policy_.thermal != ThermalMode::Disabled) {
                material_tables_ = generate_material_tables();
            }
            if (policy_.thermal == ThermalMode::Gpu) {
                resources_->thermal_workspace = std::make_unique<ThermalWorkspace>(
                    material_tables_, B * F);
                resources_->thermal_workspace->initialize_device_state(
                    material_tables_, B, F, state_.filament_masses.data(),
                    state_.reference_resistances.data(), state_.filament_materials.data(),
                    state_.temperatures.data(), state_.resistances.data(), context_->stream());
                context_->synchronize();
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
    const auto& state_snapshot = step_workspace_.state_snapshot;
    bool solver_done = false;
    auto& active = step_workspace_.active_pairs;
    for (std::size_t b = 0; b < B; ++b) {
        for (std::size_t s = 0; s < S; ++s) {
            for (std::size_t f = 0; f < F; ++f) {
                const auto index = mutual_pipeline_index(b, s, f, S, F);
                active[index] = state_.active_mask[b] != 0 && mutual_stage_mask_[b * S + s] != 0;
            }
        }
    }

#if defined(COILGUN_CUDA_AVAILABLE)
    if (context_available()) {
        auto* d_coils = static_cast<CoilGeo*>(resources_->buffers[Coils]);
        auto* d_filaments = static_cast<FilGeo*>(resources_->buffers[Filaments]);
        auto* d_separations = static_cast<double*>(resources_->buffers[Separations]);
        auto* d_active = static_cast<std::uint8_t*>(resources_->buffers[PairActive]);
        auto* d_batch_active = static_cast<std::uint8_t*>(resources_->buffers[BatchActive]);
        auto* d_trigger = static_cast<std::uint8_t*>(resources_->buffers[TriggerMask]);
        auto* d_mutual = static_cast<double*>(resources_->buffers[Mutual]);
        auto* d_gradient = static_cast<double*>(resources_->buffers[Gradient]);
        auto* d_currents = static_cast<double*>(resources_->buffers[Currents]);
        auto* d_pre_step_currents = static_cast<double*>(resources_->buffers[PreStepCurrents]);
        auto* d_derivative = static_cast<double*>(resources_->buffers[Derivative]);
        auto* d_acceleration = static_cast<double*>(resources_->buffers[Acceleration]);
        auto* d_velocity = static_cast<double*>(resources_->buffers[Velocity]);
        auto* d_position = static_cast<double*>(resources_->buffers[Position]);
        auto* d_force = static_cast<double*>(resources_->buffers[Force]);
        auto* d_stage_mask = static_cast<std::uint8_t*>(resources_->buffers[StageMask]);
        auto* d_mutual_stage_mask = static_cast<std::uint8_t*>(resources_->buffers[MutualStageMask]);
        auto* d_stage_voltage = static_cast<double*>(resources_->buffers[StageVoltage]);
        auto* d_dynamic_resistance = static_cast<double*>(resources_->buffers[DynamicResistance]);
        auto* d_matrix = static_cast<double*>(resources_->buffers[Matrix]);
        auto* d_rhs = static_cast<double*>(resources_->buffers[Rhs]);
        auto* d_solution = static_cast<double*>(resources_->buffers[Solution]);
        auto* d_residual = static_cast<double*>(resources_->buffers[Residual]);
        auto* d_status = static_cast<DeviceStepStatus*>(resources_->buffers[CompactStatus]);
        auto* d_trigger_modes = static_cast<std::uint8_t*>(resources_->buffers[TriggerMode]);
        auto* d_trigger_values = static_cast<double*>(resources_->buffers[TriggerValue]);
        auto* d_excitation_finished = static_cast<std::uint8_t*>(resources_->buffers[ExcitationFinished]);
        auto* d_stage_completed = static_cast<std::uint8_t*>(resources_->buffers[StageCompleted]);
        auto* d_trigger_times = static_cast<double*>(resources_->buffers[TriggerTime]);
        auto* d_trigger_positions = static_cast<double*>(resources_->buffers[TriggerPosition]);
        auto* d_position_offsets = static_cast<double*>(resources_->buffers[PositionOffset]);
        auto* d_control_time = static_cast<double*>(resources_->buffers[ControlTime]);
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
            copy(d_active, active.data(), sizeof(std::uint8_t) * active.size(), cudaMemcpyHostToDevice, "active upload");
            copy(d_batch_active, state_.active_mask.data(), sizeof(std::uint8_t) * B, cudaMemcpyHostToDevice, "batch active upload");
            copy(d_trigger, state_.trigger_mask.data(), sizeof(std::uint8_t) * state_.trigger_mask.size(), cudaMemcpyHostToDevice, "trigger upload");
            copy(d_stage_mask, stage_mask_.data(), sizeof(std::uint8_t) * B * S,
                 cudaMemcpyHostToDevice, "stage mask upload");
            copy(d_mutual_stage_mask, mutual_stage_mask_.data(), sizeof(std::uint8_t) * B * S,
                 cudaMemcpyHostToDevice, "mutual stage mask upload");
            auto& stage_voltages = step_workspace_.stage_voltages;
            if (state_.stage_voltages.empty())
                std::fill(stage_voltages.begin(), stage_voltages.end(), 0.0);
            else
                stage_voltages = state_.stage_voltages;
            copy(d_stage_voltage, stage_voltages.data(), sizeof(double) * B * S,
                 cudaMemcpyHostToDevice, "stage voltage upload");
            if (policy_.thermal != ThermalMode::Gpu && state_.resistances.size() == B * F)
                copy(d_dynamic_resistance, state_.resistances.data(), sizeof(double) * B * F,
                     cudaMemcpyHostToDevice, "dynamic resistance upload");
            if (device_control_enabled()) {
                copy(d_trigger_modes, state_.trigger_modes.data(), B * S * sizeof(std::uint8_t),
                     cudaMemcpyHostToDevice, "trigger mode upload");
                copy(d_trigger_values, state_.trigger_values.data(), B * S * sizeof(double),
                     cudaMemcpyHostToDevice, "trigger value upload");
                copy(d_excitation_finished, state_.excitation_finished.data(),
                     B * S * sizeof(std::uint8_t), cudaMemcpyHostToDevice,
                     "excitation completion upload");
                copy(d_stage_completed, state_.stage_completed.data(),
                     B * S * sizeof(std::uint8_t), cudaMemcpyHostToDevice,
                     "stage completion upload");
                copy(d_trigger_times, state_.trigger_times.data(), B * S * sizeof(double),
                     cudaMemcpyHostToDevice, "trigger time upload");
                copy(d_trigger_positions, state_.trigger_positions.data(), B * S * sizeof(double),
                     cudaMemcpyHostToDevice, "trigger position upload");
                copy(d_position_offsets, state_.position_offsets.data(), B * sizeof(double),
                     cudaMemcpyHostToDevice, "position offset upload");
                const double control_time = (result_.completed_steps + 1) * state_.dt;
                copy(d_control_time, &control_time, sizeof(double), cudaMemcpyHostToDevice,
                     "control time upload");
            }
            const auto opt_level = policy_.precision == PrecisionMode::Aggressive
                ? GpuOptLevel::Aggressive
                : (policy_.precision == PrecisionMode::Full
                       ? GpuOptLevel::Full : GpuOptLevel::Standard);
            const auto thermal_precision = policy_.precision == PrecisionMode::Aggressive
                ? ThermalPrecision::Aggressive
                : (policy_.precision == PrecisionMode::Full
                       ? ThermalPrecision::Full : ThermalPrecision::Standard);
            const double* assembly_resistances = policy_.thermal == ThermalMode::Gpu
                ? resources_->thermal_workspace->device_resistances()
                : (state_.resistances.size() == B * F ? d_dynamic_resistance : nullptr);
            const auto assembly_view = DeviceAssemblyView{
                B, S, F,
                static_cast<double*>(resources_->buffers[StageInductance]),
                static_cast<double*>(resources_->buffers[StageResistance]),
                static_cast<double*>(resources_->buffers[StageMutual]),
                static_cast<double*>(resources_->buffers[FilamentInductance]),
                static_cast<double*>(resources_->buffers[FilamentReferenceResistance]),
                static_cast<double*>(resources_->buffers[FilamentMutual]),
                assembly_resistances, d_mutual, d_gradient, d_currents, d_velocity,
                d_stage_voltage, d_batch_active, d_trigger, d_stage_mask,
                d_mutual_stage_mask, d_matrix, d_rhs};

            const auto launch_device_step = [&](cudaStream_t stream) {
                if (cudaMemcpyAsync(d_pre_step_currents, d_currents,
                                    B * (S + F) * sizeof(double),
                                    cudaMemcpyDeviceToDevice, stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "pre-step current snapshot failed");
                if (launch_separation_update(B, S, F,
                        static_cast<double*>(resources_->buffers[StagePosition]),
                        static_cast<double*>(resources_->buffers[FilamentPosition]),
                        d_position, d_separations, stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "separation update failed");
                launch_mutual_pipeline(MutualPipelineView{
                        d_coils, d_filaments, d_separations, d_active,
                        d_mutual, d_gradient, B, S, F, 9},
                    opt_level, config_.threads_per_block, stream);
                if (cudaGetLastError() != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "mutual pipeline failed");
                if (launch_device_assembly(assembly_view, stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "device assembly failed");
                const auto solve = solver_->solve_device(
                    DeviceMatrixView{d_matrix, B, S + F},
                    DeviceVectorView{d_rhs, B, S + F},
                    DeviceVectorView{d_solution, B, S + F},
                    DeviceResidualView{d_residual, B});
                if (!solve.ok)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody, 0,
                                                       solve.message);
                if (launch_state_update_masked(
                        B, S, F, d_currents, d_solution, d_gradient,
                        d_trigger, d_batch_active, state_.mass, state_.dt,
                        d_acceleration, d_velocity, d_position, d_force,
                        StateKernelConfig{config_.deterministic,
                            static_cast<unsigned int>(config_.threads_per_block)},
                        stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "state update failed");
                if (policy_.thermal == ThermalMode::Gpu &&
                    resources_->thermal_workspace->launch_device(
                        thermal_precision, B, S, F, d_pre_step_currents,
                        d_batch_active, state_.dt, stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "thermal update failed");
                if (device_control_enabled() && launch_device_control(DeviceControlView{
                        B, S, F, S + F, 1.0e-6, d_control_time,
                        d_currents, d_position, d_position_offsets, d_trigger_modes,
                        d_trigger_values, d_excitation_finished, d_batch_active, d_trigger,
                        d_stage_mask, d_mutual_stage_mask, d_stage_completed, d_active,
                        d_trigger_times, d_trigger_positions}, stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "device control failed");
                if (launch_compact_status(B, S + F, d_currents, d_velocity,
                        d_position, d_residual, d_batch_active, d_status, stream) != cudaSuccess)
                    return GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                        static_cast<int>(cudaGetLastError()), "compact status failed");
                return GraphCaptureStatus::success(
                    reinterpret_cast<std::uintptr_t>(d_status),
                    B * sizeof(DeviceStepStatus));
            };

            bool device_step_complete = false;
            if (policy_.backend == BackendMode::Graph && graph_cache_) {
                if (fault_injection_.fail_graph_capture) {
                    fault_injection_.fail_graph_capture = false;
                    throw std::runtime_error("injected CUDA graph capture failure");
                }
                GpuGraphBoundaryState boundary;
                boundary.topology = graph_key();
                boundary.runtime_masks.stage_mask = stage_mask_;
                boundary.runtime_masks.mutual_stage_mask = mutual_stage_mask_;
                const auto capture_count = graph_cache_->capture_count();
                const auto selected = graph_cache_->capture_and_select(
                    boundary, context_->stream(), launch_device_step,
                    GraphWorkspace{reinterpret_cast<std::uintptr_t>(d_status),
                                   B * sizeof(DeviceStepStatus)});
                if (!selected.ok)
                    throw std::runtime_error("CUDA graph capture failed: " +
                                             selected.failure.message);
                if (graph_cache_->capture_count() > capture_count)
                    ++report_.graph_rebuild_count;
                const auto replayed = graph_cache_->replay(context_->stream());
                if (!replayed.ok)
                    throw std::runtime_error("CUDA graph replay failed: " +
                                             replayed.failure.message);
                device_step_complete = true;
            } else if (solver_->resolved_mode() == SolverMode::Batched) {
                const auto launched = launch_device_step(context_->stream());
                if (!launched.ok) throw std::runtime_error(launched.failure.message);
                device_step_complete = true;
            }

            if (!device_step_complete) {
                if (launch_separation_update(B, S, F,
                        static_cast<double*>(resources_->buffers[StagePosition]),
                        static_cast<double*>(resources_->buffers[FilamentPosition]),
                        d_position, d_separations, context_->stream()) != cudaSuccess)
                    throw std::runtime_error("separation update kernel launch failed");
                launch_mutual_pipeline(MutualPipelineView{
                        d_coils, d_filaments, d_separations, d_active,
                        d_mutual, d_gradient, B, S, F, 9},
                    opt_level, config_.threads_per_block, context_->stream());
                context_->synchronize();
                if (fault_injection_.fail_after_mutual) {
                    fault_injection_.fail_after_mutual = false;
                    throw std::runtime_error("injected GPU failure after mutual segment");
                }
                if (launch_device_assembly(assembly_view, context_->stream()) != cudaSuccess)
                    throw std::runtime_error("device matrix/RHS assembly failed");
                context_->synchronize();
                copy(matrices_.data(), d_matrix, sizeof(double) * matrices_.size(),
                     cudaMemcpyDeviceToHost, "matrix fallback download");
                copy(rhs_.data(), d_rhs, sizeof(double) * rhs_.size(),
                     cudaMemcpyDeviceToHost, "rhs fallback download");
                execute_solver_step();
                copy(d_derivative, solution_.data(), sizeof(double) * solution_.size(),
                     cudaMemcpyHostToDevice, "derivative fallback upload");
                if (cudaMemcpyAsync(d_pre_step_currents, d_currents,
                                    B * (S + F) * sizeof(double),
                                    cudaMemcpyDeviceToDevice, context_->stream()) != cudaSuccess)
                    throw std::runtime_error("pre-step current snapshot failed");
                if (launch_state_update_masked(
                        B, S, F, d_currents, d_derivative, d_gradient,
                        d_trigger, d_batch_active, state_.mass, state_.dt,
                        d_acceleration, d_velocity, d_position, d_force,
                        StateKernelConfig{config_.deterministic,
                            static_cast<unsigned int>(config_.threads_per_block)},
                        context_->stream()) != cudaSuccess)
                    throw std::runtime_error("state kernel launch failed");
                if (policy_.thermal == ThermalMode::Gpu &&
                    resources_->thermal_workspace->launch_device(
                        thermal_precision, B, S, F, d_pre_step_currents,
                        d_batch_active, state_.dt, context_->stream()) != cudaSuccess)
                    throw std::runtime_error("thermal update failed");
                if (device_control_enabled() && launch_device_control(DeviceControlView{
                        B, S, F, S + F, 1.0e-6, d_control_time,
                        d_currents, d_position, d_position_offsets, d_trigger_modes,
                        d_trigger_values, d_excitation_finished, d_batch_active, d_trigger,
                        d_stage_mask, d_mutual_stage_mask, d_stage_completed, d_active,
                        d_trigger_times, d_trigger_positions}, context_->stream()) != cudaSuccess)
                    throw std::runtime_error("device control failed");
                if (launch_compact_status(B, S + F, d_currents, d_velocity,
                        d_position, nullptr, d_batch_active, d_status,
                        context_->stream()) != cudaSuccess)
                    throw std::runtime_error("compact status failed");
            }

            context_->synchronize();
            if (solver_->resolved_mode() == SolverMode::Batched) {
                const auto solver_status = solver_->validate_device_result(
                    DeviceResidualView{d_residual, B});
                if (!solver_status.ok) throw std::runtime_error(solver_status.message);
            }
            auto& compact_status = step_workspace_.compact_status;
            compact_status.resize(B);
            copy(compact_status.data(), d_status, B * sizeof(DeviceStepStatus),
                 cudaMemcpyDeviceToHost, "compact status download");
            for (std::size_t batch = 0; batch < B; ++batch) {
                if (compact_status[batch].finite == 0 ||
                    compact_status[batch].solver_ok == 0)
                    throw std::runtime_error("device compact status rejected the step");
            }

            copy(state_.m1.data(), d_mutual, sizeof(double) * state_.m1.size(),
                 cudaMemcpyDeviceToHost, "mutual download");
            copy(state_.dm1.data(), d_gradient, sizeof(double) * state_.dm1.size(),
                 cudaMemcpyDeviceToHost, "gradient download");
            copy(state_.currents.data(), d_currents, sizeof(double) * state_.currents.size(),
                 cudaMemcpyDeviceToHost, "current download");
            copy(state_.velocity.data(), d_velocity, sizeof(double) * B,
                 cudaMemcpyDeviceToHost, "velocity download");
            copy(state_.position.data(), d_position, sizeof(double) * B,
                 cudaMemcpyDeviceToHost, "position download");
            state_.current_derivatives.resize(B * (S + F));
            copy(state_.current_derivatives.data(),
                 solver_->resolved_mode() == SolverMode::Batched ? d_solution : d_derivative,
                 sizeof(double) * state_.current_derivatives.size(),
                 cudaMemcpyDeviceToHost, "derivative download");
            current_derivatives_ = state_.current_derivatives;
            if (device_control_enabled()) {
                copy(state_.active_mask.data(), d_batch_active, B * sizeof(std::uint8_t),
                     cudaMemcpyDeviceToHost, "active status download");
                copy(state_.trigger_mask.data(), d_trigger, B * S * sizeof(std::uint8_t),
                     cudaMemcpyDeviceToHost, "trigger status download");
                copy(state_.stage_mask.data(), d_stage_mask, B * S * sizeof(std::uint8_t),
                     cudaMemcpyDeviceToHost, "stage status download");
                copy(state_.mutual_stage_mask.data(), d_mutual_stage_mask,
                     B * S * sizeof(std::uint8_t), cudaMemcpyDeviceToHost,
                     "mutual status download");
                copy(state_.stage_completed.data(), d_stage_completed,
                     B * S * sizeof(std::uint8_t), cudaMemcpyDeviceToHost,
                     "completion status download");
                copy(state_.trigger_times.data(), d_trigger_times, B * S * sizeof(double),
                     cudaMemcpyDeviceToHost, "trigger time download");
                copy(state_.trigger_positions.data(), d_trigger_positions,
                     B * S * sizeof(double), cudaMemcpyDeviceToHost,
                     "trigger position download");
                stage_mask_ = state_.stage_mask;
                mutual_stage_mask_ = state_.mutual_stage_mask;
                select_graph_variant_at_boundary();
            }
            pipeline_order_.push_back(PipelineStage::Matrix);
            pipeline_order_.push_back(PipelineStage::Solver);
            pipeline_order_.push_back(PipelineStage::Force);
            solver_done = true;
            pipeline_order_.push_back(PipelineStage::State);

            if (policy_.thermal != ThermalMode::Disabled && !state_.temperatures.empty()) {
                const auto thermal_start = std::chrono::steady_clock::now();
                state_.resistivities.resize(B * F);
                state_.resistances.resize(B * F);
                state_.joule_energy.resize(B * F);
                if (policy_.thermal == ThermalMode::Cpu) {
                    auto& filament_currents = step_workspace_.thermal_filament_currents;
                    std::fill(filament_currents.begin(), filament_currents.end(), 0.0);
                    for (std::size_t b = 0; b < B; ++b)
                        if (state_.active_mask[b] != 0)
                            std::copy_n(state_snapshot.currents.data() + b * layout_.D + S,
                                        F, filament_currents.data() + b * F);
                    update_thermal_batch_cpu(
                        material_tables_, thermal_precision, B, F,
                        filament_currents.data(), state_.filament_masses.data(),
                        state_.reference_resistances.data(), state_.filament_materials.data(),
                        state_.dt, state_.temperatures.data(), state_.resistivities.data(),
                        state_.resistances.data(), state_.joule_energy.data());
                } else if (resources_->thermal_workspace->download_device_state(
                               state_.temperatures.data(), state_.resistivities.data(),
                               state_.resistances.data(), state_.joule_energy.data(),
                               context_->stream()) != cudaSuccess) {
                    throw std::runtime_error("resident thermal observation failed");
                } else {
                    context_->synchronize();
                }
                report_.thermal_time_ms += std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - thermal_start).count();
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
            auto& filament_currents = step_workspace_.thermal_filament_currents;
            std::fill(filament_currents.begin(), filament_currents.end(), 0.0);
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

GpuGraphTopologyKey GpuEngine::graph_key() const {
    GpuGraphTopologyKey key;
    key.batch_capacity = layout_.B;
    key.layout_signature = static_cast<std::uint64_t>(layout_.D) * 1099511628211ull ^ layout_.S ^ (layout_.F << 16);
    key.stage_signature = static_cast<std::uint64_t>(layout_.S) ^
        (static_cast<std::uint64_t>(layout_.F) << 32);
    key.precision = policy_.precision;
    key.thermal = policy_.thermal;
    key.solver = policy_.solver;
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
