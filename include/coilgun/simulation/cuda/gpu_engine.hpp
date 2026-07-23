/**
 * @file gpu_engine.hpp
 * @brief Synchronous contract for the unified GPU execution engine.
 *
 * This first contract is host-callable and deliberately kernel-neutral. The
 * implementation keeps ownership and step-boundary semantics stable while the
 * CUDA pipeline performs the physical mutual-inductance work.
 */

#pragma once

#include "coilgun/simulation/cuda/gpu_execution_config.hpp"
#include "coilgun/simulation/cuda/gpu_execution_report.hpp"
#include "coilgun/simulation/cuda/gpu_state_layout.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/self_inductance.hpp"
#include <Eigen/Dense>

#if defined(COILGUN_CUDA_AVAILABLE)
#include "coilgun/simulation/cuda/gpu_execution_context.hpp"
#include "coilgun/simulation/cuda/gpu_graph.hpp"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/simulation/cuda/gpu_solver.hpp"
#include "coilgun/simulation/cuda/gpu_mutual_pipeline.hpp"
#include "coilgun/simulation/cuda/gpu_thermal.hpp"
#include "coilgun/simulation/cuda/persistent_kernel.cuh"
#endif

#include <cstddef>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace coilgun::simulation::cuda {

namespace detail {

/** Test-only engine seam; never part of ordinary execution configuration. */
struct GpuEngineFaultInjection {
    bool fail_after_mutual = false;
    bool fail_allocation = false;
    bool fail_device_initialization = false;
    bool fail_graph_capture = false;
};

} // namespace detail

enum class PipelineStage { Mutual, Matrix, Solver, Force, Thermal, State };

#if defined(COILGUN_CUDA_AVAILABLE)
bool cuda_device_available() noexcept;
std::unique_ptr<GpuExecutionContext> make_gpu_execution_context();
#endif

struct GpuGeometryInput {
    std::size_t n_stages = 0;
    std::size_t n_filaments = 0;
    std::vector<double> stage_geometry;
    std::vector<double> filament_geometry;
    bool thermal_enabled = false;
    std::vector<double> stage_inner_radii;
    std::vector<double> stage_outer_radii;
    std::vector<double> stage_lengths;
    std::vector<int> stage_turns;
    std::vector<double> stage_positions;
    std::vector<double> filament_inner_radii;
    std::vector<double> filament_outer_radii;
    std::vector<double> filament_lengths;
    std::vector<double> filament_positions;
    std::vector<double> stage_resistances;
    std::vector<double> stage_inductances;
    std::vector<double> filament_resistances;
    std::vector<double> filament_inductances;
    std::vector<double> stage_mutual_inductances;
    std::vector<double> filament_mutual_inductances;

    void validate() const {
        if (n_stages == 0 || n_filaments == 0) {
            throw std::invalid_argument("GPU geometry dimensions must be positive");
        }
        if (stage_geometry.size() != n_stages || filament_geometry.size() != n_filaments) {
            throw std::invalid_argument("GPU geometry buffers do not match dimensions");
        }
        for (const double value : stage_geometry) {
            if (!std::isfinite(value)) throw std::invalid_argument("non-finite stage geometry");
        }
        for (const double value : filament_geometry) {
            if (!std::isfinite(value)) throw std::invalid_argument("non-finite filament geometry");
        }
        const auto check_size = [](std::size_t actual, std::size_t expected, const char* name) {
            if (actual != expected) throw std::invalid_argument(name);
        };
        check_size(stage_inner_radii.size(), n_stages, "stage inner radii are required");
        check_size(stage_outer_radii.size(), n_stages, "stage outer radii are required");
        check_size(stage_lengths.size(), n_stages, "stage lengths are required");
        check_size(stage_turns.size(), n_stages, "stage turns are required");
        check_size(stage_positions.size(), n_stages, "stage positions are required");
        check_size(filament_inner_radii.size(), n_filaments, "filament inner radii are required");
        check_size(filament_outer_radii.size(), n_filaments, "filament outer radii are required");
        check_size(filament_lengths.size(), n_filaments, "filament lengths are required");
        check_size(filament_positions.size(), n_filaments, "filament positions are required");
        if (!stage_resistances.empty() && stage_resistances.size() != n_stages)
            throw std::invalid_argument("stage resistances do not match dimensions");
        if (!stage_inductances.empty() && stage_inductances.size() != n_stages)
            throw std::invalid_argument("stage inductances do not match dimensions");
        if (!filament_resistances.empty() && filament_resistances.size() != n_filaments)
            throw std::invalid_argument("filament resistances do not match dimensions");
        if (!filament_inductances.empty() && filament_inductances.size() != n_filaments)
            throw std::invalid_argument("filament inductances do not match dimensions");
        if (!stage_mutual_inductances.empty() && stage_mutual_inductances.size() != n_stages * n_stages)
            throw std::invalid_argument("stage mutual inductances do not match dimensions");
        if (!filament_mutual_inductances.empty() && filament_mutual_inductances.size() != n_filaments * n_filaments)
            throw std::invalid_argument("filament mutual inductances do not match dimensions");
        const auto finite_positive = [](double value, const char* name) {
            if (!std::isfinite(value) || value <= 0.0) throw std::invalid_argument(name);
        };
        for (std::size_t i = 0; i < n_stages; ++i) {
            finite_positive(stage_inner_radii[i], "invalid stage inner radius");
            finite_positive(stage_outer_radii[i], "invalid stage outer radius");
            finite_positive(stage_lengths[i], "invalid stage length");
            if (stage_outer_radii[i] <= stage_inner_radii[i] || stage_turns[i] < 2 ||
                !std::isfinite(stage_positions[i])) throw std::invalid_argument("invalid stage geometry");
        }
        for (std::size_t i = 0; i < n_filaments; ++i) {
            finite_positive(filament_inner_radii[i], "invalid filament inner radius");
            finite_positive(filament_outer_radii[i], "invalid filament outer radius");
            finite_positive(filament_lengths[i], "invalid filament length");
            if (filament_outer_radii[i] <= filament_inner_radii[i] ||
                !std::isfinite(filament_positions[i])) throw std::invalid_argument("invalid filament geometry");
        }
    }
};

struct GpuEngineState {
    std::vector<double> currents;
    std::vector<double> m1;
    std::vector<double> dm1;
    std::vector<double> temperatures;
    std::vector<std::uint8_t> active_mask;
    std::vector<std::uint8_t> trigger_mask;
    std::vector<std::uint8_t> stage_mask;
    std::vector<std::uint8_t> mutual_stage_mask;
    std::vector<double> velocity;
    std::vector<double> position;
    std::vector<double> filament_masses;
    std::vector<double> reference_resistances;
    std::vector<int> filament_materials;
    std::vector<double> resistivities;
    std::vector<double> resistances;
    std::vector<double> joule_energy;
    std::vector<double> current_derivatives;
    std::vector<double> stage_voltages;
    std::vector<std::uint8_t> trigger_modes;
    std::vector<double> trigger_values;
    std::vector<std::uint8_t> excitation_finished;
    std::vector<std::uint8_t> stage_completed;
    std::vector<double> trigger_times;
    std::vector<double> trigger_positions;
    std::vector<double> position_offsets;
    double dt = 0.0;
    double mass = 0.0;
    double reference_temperature = 293.0;
    double material_density = 0.0;
};

struct GpuEngineResult {
    std::size_t completed_steps = 0;
    bool finished = false;
};

struct GpuAssemblySnapshot {
    std::vector<double> matrix;
    std::vector<double> rhs;
};

struct GpuRunBoundary {
    std::size_t max_steps = std::numeric_limits<std::size_t>::max();
    bool stop_when_inactive = true;
};

struct GpuGraphVariant {
    std::vector<std::uint8_t> stage_mask;
    std::vector<std::uint8_t> mutual_stage_mask;
    std::size_t batch_size = 0;
    std::size_t batch_capacity = 0;
    std::size_t layout_dimension = 0;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Disabled;
    SolverMode solver = SolverMode::Eigen;

    bool operator==(const GpuGraphVariant& other) const noexcept {
        return batch_size == other.batch_size && batch_capacity == other.batch_capacity &&
               layout_dimension == other.layout_dimension && precision == other.precision &&
               thermal == other.thermal && solver == other.solver && stage_mask == other.stage_mask &&
               mutual_stage_mask == other.mutual_stage_mask;
    }
};

class GpuEngine {
public:
    GpuEngine(GpuGeometryInput geometry,
              GpuEngineState initial_state,
              GpuExecutionConfig config = {},
              GpuCapability capability = {},
              detail::GpuEngineFaultInjection fault_injection = {})
        : resources_(std::make_unique<Resources>()),
          geometry_(std::move(geometry)),
          layout_(validate_and_layout(geometry_, initial_state, config)),
          initial_state_(std::move(initial_state)),
          state_(initial_state_),
          config_(config),
           policy_(GpuExecutionPlanner::plan(geometry_.n_stages,
                                             geometry_.n_filaments,
                                             layout_.B,
                                            geometry_.thermal_enabled ||
                                                config.thermal == ThermalMode::Cpu ||
                                                config.thermal == ThermalMode::Gpu,
                                                 capability,
                                                  config_)) {
          fault_injection_ = fault_injection;
          if (geometry_.stage_mutual_inductances.empty()) {
              geometry_.stage_mutual_inductances.assign(
                  geometry_.n_stages * geometry_.n_stages, 0.0);
              for (std::size_t first = 0; first < geometry_.n_stages; ++first) {
                  for (std::size_t second = first + 1; second < geometry_.n_stages; ++second) {
                      const double mutual = physics::mutual_inductance_coil(
                          geometry_.stage_inner_radii[first], geometry_.stage_outer_radii[first],
                          geometry_.stage_lengths[first], geometry_.stage_turns[first],
                          geometry_.stage_inner_radii[second], geometry_.stage_outer_radii[second],
                          geometry_.stage_lengths[second], geometry_.stage_turns[second],
                          std::abs(geometry_.stage_positions[first] - geometry_.stage_positions[second]),
                          9, true);
                      geometry_.stage_mutual_inductances[first * geometry_.n_stages + second] = mutual;
                      geometry_.stage_mutual_inductances[second * geometry_.n_stages + first] = mutual;
                  }
              }
          }
          if (geometry_.filament_mutual_inductances.empty()) {
              geometry_.filament_mutual_inductances.assign(
                  geometry_.n_filaments * geometry_.n_filaments, 0.0);
              const auto radius = [](const GpuGeometryInput& geometry, std::size_t filament) {
                  return 0.5 * (geometry.filament_inner_radii[filament] +
                                geometry.filament_outer_radii[filament]);
              };
              for (std::size_t first = 0; first < geometry_.n_filaments; ++first) {
                  for (std::size_t second = first + 1; second < geometry_.n_filaments; ++second) {
                      const double mutual = physics::mutual_inductance_filament(
                          radius(geometry_, first), radius(geometry_, second),
                          std::abs(geometry_.filament_positions[first] -
                                   geometry_.filament_positions[second]), true);
                      geometry_.filament_mutual_inductances[first * geometry_.n_filaments + second] = mutual;
                      geometry_.filament_mutual_inductances[second * geometry_.n_filaments + first] = mutual;
                  }
              }
          }
          bool has_active_simulation = false;
         for (const auto active : initial_state_.active_mask) {
             if (active != 0) {
                 has_active_simulation = true;
                 break;
             }
         }
  #if defined(COILGUN_CUDA_AVAILABLE)
           if (config_.backend == BackendMode::Graph && capability.supports_graph &&
               has_active_simulation) {
               policy_.backend = BackendMode::Graph;
               policy_.backend_fallback_reason = FallbackReason::None;
           } else if (config_.backend == BackendMode::Persistent && layout_.B == 1 &&
                      capability.supports_persistent &&
                       (!config_.deterministic || capability.persistent_is_deterministic)) {
                policy_.backend = BackendMode::Persistent;
                policy_.backend_fallback_reason = FallbackReason::None;
            } else if (config_.backend == BackendMode::Graph && !capability.supports_graph) {
                policy_.backend_fallback_reason = FallbackReason::CapabilityUnavailable;
            } else if (config_.backend == BackendMode::Persistent &&
                       (layout_.B != 1 || !capability.supports_persistent)) {
                policy_.backend_fallback_reason = FallbackReason::CapabilityUnavailable;
            }
 #endif
          if (state_.velocity.empty()) state_.velocity.assign(layout_.B, 0.0);
          if (state_.position.empty()) state_.position.assign(layout_.B, 0.0);
        if (state_.velocity.size() != layout_.B || state_.position.size() != layout_.B) {
              throw std::invalid_argument("GPU position and velocity buffers do not match batch size");
          }
          if (!std::isfinite(state_.dt) || state_.dt <= 0.0 || !std::isfinite(state_.mass) || state_.mass <= 0.0) {
              throw std::invalid_argument("GPU state requires positive finite dt and mass");
          }
           if (geometry_.thermal_enabled || config_.thermal == ThermalMode::Cpu || config_.thermal == ThermalMode::Gpu) {
              if (!std::isfinite(state_.reference_temperature) || state_.reference_temperature <= 0.0 ||
                  !std::isfinite(state_.material_density) || state_.material_density <= 0.0 ||
                  state_.filament_masses.size() != layout_.B * layout_.F ||
                  state_.reference_resistances.size() != layout_.B * layout_.F ||
                  state_.filament_materials.size() != layout_.B * layout_.F) {
                   throw std::invalid_argument("thermal-enabled GPU state requires material tables and references");
               }
               if (state_.resistances.empty()) {
                   state_.resistances = state_.reference_resistances;
                   initial_state_.resistances = state_.reference_resistances;
               }
           }
         matrices_.assign(layout_.system_matrix_size(), 0.0);
         rhs_.assign(layout_.rhs_size(), 0.0);
         solution_.assign(layout_.rhs_size(), 0.0);
        report_.requested_backend = policy_.requested_backend;
        report_.requested_solver = policy_.requested_solver;
        report_.requested_precision = policy_.requested_precision;
        report_.requested_thermal = policy_.requested_thermal;
        report_.backend = policy_.backend;
        report_.solver = policy_.solver;
        report_.precision = policy_.precision;
         report_.thermal = policy_.thermal;
         report_.device_id = config_.device_id;
         report_.threads_per_block = config_.threads_per_block;
         report_.profiling_enabled = config_.enable_profiling;
         report_.static_fallback_reason =
             policy_.backend_fallback_reason != FallbackReason::None
                 ? policy_.backend_fallback_reason
                 : (policy_.solver_fallback_reason != FallbackReason::None
                        ? policy_.solver_fallback_reason : policy_.thermal_fallback_reason);
          if (policy_.backend_fallback_reason != FallbackReason::None) {
              report_.fallback_reason = policy_.backend_fallback_reason == FallbackReason::CapabilityUnavailable
                  ? "requested backend capability is unavailable; using CPU fallback"
                  : "requested backend is unavailable in the current engine; using CPU fallback";
               ++report_.fallback_count;
           }
#if !defined(COILGUN_CUDA_AVAILABLE)
            if (config_.backend != BackendMode::Fallback) {
               policy_.backend = BackendMode::Fallback;
               policy_.solver = SolverMode::Eigen;
              if (policy_.thermal == ThermalMode::Gpu) policy_.thermal = ThermalMode::Cpu;
              report_.backend = BackendMode::Fallback;
              report_.solver = SolverMode::Eigen;
              report_.thermal = policy_.thermal;
              report_.static_fallback_reason = FallbackReason::CapabilityUnavailable;
               report_.fallback_reason = "CUDA support is not compiled; using CPU fallback";
               ++report_.fallback_count;
           }
           if (config_.backend == BackendMode::Fallback) {
               policy_.solver = SolverMode::Eigen;
               if (policy_.thermal == ThermalMode::Gpu) policy_.thermal = ThermalMode::Cpu;
               report_.solver = SolverMode::Eigen;
               report_.thermal = policy_.thermal;
               if (report_.fallback_reason.empty())
                   report_.fallback_reason = "CPU fallback explicitly requested";
           }
#endif
         normalize_boundary_masks(initial_state_);
         normalize_boundary_masks(state_);
         stage_mask_ = state_.stage_mask;
         mutual_stage_mask_ = state_.mutual_stage_mask;
        select_graph_variant_at_boundary();
#if defined(COILGUN_CUDA_AVAILABLE)
        initialize_runtime();
#endif
    }

    ~GpuEngine() { shutdown(); }

    GpuEngine(const GpuEngine&) = delete;
    GpuEngine& operator=(const GpuEngine&) = delete;
    GpuEngine(GpuEngine&&) = delete;
    GpuEngine& operator=(GpuEngine&&) = delete;

    /** Execute exactly one synchronous step through the resolved solver. */
    void step() {
        ensure_running();
#if defined(COILGUN_CUDA_AVAILABLE)
        const auto state_snapshot = state_;
        const auto result_snapshot = result_;
        const auto stage_mask_snapshot = stage_mask_;
        const auto mutual_stage_mask_snapshot = mutual_stage_mask_;
        const auto variant_snapshot = variant_;
        const auto selected_variant_snapshot = selected_variant_;
        const auto pipeline_snapshot = pipeline_order_;
        const auto report_snapshot = report_;
        const auto policy_snapshot = policy_;
        const auto matrices_snapshot = matrices_;
        const auto rhs_snapshot = rhs_;
        const auto solution_snapshot = solution_;
        const auto derivatives_snapshot = current_derivatives_;
        const bool variant_valid_snapshot = variant_valid_;
        const bool fallback_locked_snapshot = fallback_locked_;
        try {
            execute_physical_pipeline();
            report_.gpu_executed = report_.gpu_executed ||
                                   (context_available() && !fallback_locked_);
        } catch (const std::exception& error) {
            state_ = state_snapshot;
            result_ = result_snapshot;
            stage_mask_ = stage_mask_snapshot;
            mutual_stage_mask_ = mutual_stage_mask_snapshot;
            variant_ = variant_snapshot;
            selected_variant_ = selected_variant_snapshot;
            pipeline_order_ = pipeline_snapshot;
            report_ = report_snapshot;
            policy_ = policy_snapshot;
            matrices_ = matrices_snapshot;
            rhs_ = rhs_snapshot;
            solution_ = solution_snapshot;
            current_derivatives_ = derivatives_snapshot;
            variant_valid_ = variant_valid_snapshot;
            fallback_locked_ = fallback_locked_snapshot;
            if (!fallback_locked_snapshot) {
                SolverStatus status = SolverStatus::failure_status(
                    SolverFailure::InvalidArgument, error.what());
                record_runtime_failure(status);
            }
            try {
                execute_cpu_physical_pipeline();
            } catch (...) {
                state_ = state_snapshot;
                result_ = result_snapshot;
                stage_mask_ = stage_mask_snapshot;
                mutual_stage_mask_ = mutual_stage_mask_snapshot;
                variant_ = variant_snapshot;
                selected_variant_ = selected_variant_snapshot;
                pipeline_order_ = pipeline_snapshot;
                matrices_ = matrices_snapshot;
                rhs_ = rhs_snapshot;
                solution_ = solution_snapshot;
                current_derivatives_ = derivatives_snapshot;
                variant_valid_ = variant_valid_snapshot;
                fallback_locked_ = true;
                if (fallback_locked_snapshot) {
                    policy_ = policy_snapshot;
                    report_ = report_snapshot;
                }
                throw;
            }
        }
#else
        execute_cpu_physical_pipeline();
#endif
        ++result_.completed_steps;
    }

    /** Execute a fixed number of synchronous steps. */
    void run(std::size_t steps) {
        GpuRunBoundary boundary;
        boundary.max_steps = steps;
        run(boundary);
    }

    /** Execute until the requested run boundary is reached. */
    void run(const GpuRunBoundary& boundary) {
        ensure_running();
        if (boundary.max_steps == std::numeric_limits<std::size_t>::max() && any_active()) {
            throw std::logic_error("unbounded GPU run is unavailable for the stub backend");
        }
        const auto start = result_.completed_steps;
        while (result_.completed_steps - start < boundary.max_steps) {
            if (boundary.stop_when_inactive && !any_active()) {
                result_.finished = true;
                break;
            }
            step();
        }
        if (result_.completed_steps - start >= boundary.max_steps) {
            result_.finished = boundary.max_steps != std::numeric_limits<std::size_t>::max();
        }
    }

    void reset() {
        ensure_running();
        state_ = initial_state_;
        result_ = {};
        stage_mask_ = state_.stage_mask;
        mutual_stage_mask_ = state_.mutual_stage_mask;
        select_graph_variant_at_boundary();
#if defined(COILGUN_CUDA_AVAILABLE)
        sync_runtime_state_after_reset();
#endif
        // Execution diagnostics are cumulative audit data. Reset clears the
        // simulation state/history and boundary variant, but retains timings,
        // fallback count, calibration, and graph rebuild history.
    }

    /** Change stage participation only at a step boundary. */
    void set_stage_mask(std::vector<std::uint8_t> mask) {
        ensure_running();
        validate_mask(mask, layout_.B * layout_.S, "stage mask size does not match layout");
        state_.stage_mask = mask;
        stage_mask_ = std::move(mask);
        select_graph_variant_at_boundary();
    }

    /** Change mutual-inductance and force participation at a step boundary. */
    void set_mutual_stage_mask(std::vector<std::uint8_t> mask) {
        ensure_running();
        validate_mask(mask, layout_.B * layout_.S, "mutual stage mask size does not match layout");
        state_.mutual_stage_mask = mask;
        mutual_stage_mask_ = std::move(mask);
        select_graph_variant_at_boundary();
    }

    /** Replace all step-boundary controls atomically before the next step. */
    void set_step_boundary_state(std::vector<std::uint8_t> active_mask,
                                std::vector<std::uint8_t> trigger_mask,
                                std::vector<std::uint8_t> stage_mask,
                                std::vector<std::uint8_t> mutual_stage_mask,
                                std::vector<double> stage_voltages) {
        ensure_running();
        validate_mask(active_mask, layout_.B, "active mask size does not match layout");
        validate_mask(trigger_mask, layout_.B * layout_.S,
                     "trigger mask size does not match layout");
        validate_mask(stage_mask, layout_.B * layout_.S,
                     "stage mask size does not match layout");
        validate_mask(mutual_stage_mask, layout_.B * layout_.S,
                     "mutual stage mask size does not match layout");
        validate_stage_voltages(stage_voltages, layout_.B * layout_.S);

        state_.active_mask = std::move(active_mask);
        state_.trigger_mask = std::move(trigger_mask);
        state_.stage_mask = stage_mask;
        state_.mutual_stage_mask = mutual_stage_mask;
        state_.stage_voltages = std::move(stage_voltages);
        stage_mask_ = std::move(stage_mask);
        mutual_stage_mask_ = std::move(mutual_stage_mask);
        select_graph_variant_at_boundary();
    }

    void set_control_boundary_state(
        std::vector<std::uint8_t> trigger_modes,
        std::vector<double> trigger_values,
        std::vector<std::uint8_t> excitation_finished,
        std::vector<std::uint8_t> stage_completed,
        std::vector<double> trigger_times,
        std::vector<double> trigger_positions,
        std::vector<double> position_offsets) {
        ensure_running();
        const auto stage_values = layout_.B * layout_.S;
        validate_mask(trigger_modes, stage_values, "trigger mode size does not match layout");
        validate_mask(excitation_finished, stage_values,
                      "excitation finished size does not match layout");
        validate_mask(stage_completed, stage_values,
                      "stage completed size does not match layout");
        if (trigger_values.size() != stage_values || trigger_times.size() != stage_values ||
            trigger_positions.size() != stage_values || position_offsets.size() != layout_.B)
            throw std::invalid_argument("GPU control buffers do not match layout");
        state_.trigger_modes = std::move(trigger_modes);
        state_.trigger_values = std::move(trigger_values);
        state_.excitation_finished = std::move(excitation_finished);
        state_.stage_completed = std::move(stage_completed);
        state_.trigger_times = std::move(trigger_times);
        state_.trigger_positions = std::move(trigger_positions);
        state_.position_offsets = std::move(position_offsets);
    }

    /** Select or capture a graph variant at the current step boundary. */
    void select_graph_variant_at_boundary() {
        ensure_running();
        variant_.stage_mask = stage_mask_;
        variant_.mutual_stage_mask = mutual_stage_mask_;
        variant_.batch_size = layout_.B;
        variant_.batch_capacity = layout_.B;
        variant_.layout_dimension = layout_.D;
        variant_.precision = policy_.precision;
        variant_.thermal = policy_.thermal;
        variant_.solver = policy_.solver;
        if (!variant_valid_ || !(variant_ == selected_variant_)) {
            selected_variant_ = variant_;
            variant_valid_ = true;
        }
    }

    void shutdown() noexcept {
#if defined(COILGUN_CUDA_AVAILABLE)
        shutdown_runtime();
#endif
        if (resources_) {
            resources_->shutdown = true;
        }
    }

    bool is_shutdown() const noexcept {
        return !resources_ || resources_->shutdown;
    }

    const GpuStateLayout& layout() const noexcept { return layout_; }
    const GpuEngineState& state() const noexcept { return state_; }
    const GpuEngineResult& result() const noexcept { return result_; }
    const ExecutionReport& report() const noexcept { return report_; }
    const GpuExecutionPolicy& policy() const noexcept { return policy_; }
    const std::vector<PipelineStage>& pipeline_order() const noexcept { return pipeline_order_; }
    const GpuGraphVariant& graph_variant() const noexcept { return variant_; }
    void set_stage_voltage(std::size_t stage, double voltage) {
        ensure_running();
        if (stage >= layout_.S || !std::isfinite(voltage))
            throw std::invalid_argument("invalid stage voltage");
        if (layout_.B != 1)
            throw std::invalid_argument("set_stage_voltage requires a single-batch engine");
        if (state_.stage_voltages.empty())
            state_.stage_voltages.assign(layout_.B * layout_.S, 0.0);
        state_.stage_voltages[stage] = voltage;
    }
    void complete_stage(std::size_t batch, std::size_t stage);
    std::size_t calibration_count() const noexcept { return calibration_count_; }
    GpuAssemblySnapshot assemble_reference_for_test();
#if defined(COILGUN_CUDA_AVAILABLE)
    GpuAssemblySnapshot assemble_device_for_test();
#endif
#if defined(COILGUN_CUDA_AVAILABLE)
    bool context_available() const noexcept { return context_ != nullptr && context_->valid(); }
    bool solver_workspace_initialized() const noexcept {
        return solver_ != nullptr && solver_->workspace().initialized;
    }
    std::vector<std::uintptr_t> device_buffer_addresses() const;
    std::size_t device_allocation_count() const noexcept { return device_allocation_count_; }
#else
    bool context_available() const noexcept { return false; }
    bool solver_workspace_initialized() const noexcept { return false; }
#endif

private:
    struct Resources {
        bool shutdown = false;
#if defined(COILGUN_CUDA_AVAILABLE)
        ~Resources();
#else
        ~Resources() = default;
#endif
#if defined(COILGUN_CUDA_AVAILABLE)
        std::vector<void*> buffers;
        std::unique_ptr<ThermalWorkspace> thermal_workspace;
        int device_id = -1;
#endif
    };

    static GpuStateLayout validate_and_layout(const GpuGeometryInput& geometry,
                                              const GpuEngineState& state,
                                              const GpuExecutionConfig& config) {
        config.validate();
        const auto thermal_mode = config.thermal;
        geometry.validate();
        const std::size_t batch_size = state.active_mask.size();
        if (batch_size == 0) {
            throw std::invalid_argument("GPU state batch size must be positive");
        }
        const GpuStateLayout layout(batch_size, geometry.n_stages, geometry.n_filaments);
        if (state.currents.size() != layout.currents_size() ||
            state.m1.size() != layout.state_size() ||
            state.dm1.size() != layout.state_size() ||
            state.trigger_mask.size() != layout.trigger_mask_size()) {
            throw std::invalid_argument("GPU state buffers do not match layout");
        }
        validate_optional_mask(state.stage_mask, layout.B * layout.S,
                               "stage mask size does not match layout");
        validate_optional_mask(state.mutual_stage_mask, layout.B * layout.S,
                               "mutual stage mask size does not match layout");
        validate_stage_voltages(state.stage_voltages, layout.B * layout.S);
        const auto control_size = layout.B * layout.S;
        const bool control_enabled = !state.trigger_modes.empty() ||
            !state.trigger_values.empty() || !state.excitation_finished.empty() ||
            !state.stage_completed.empty() || !state.trigger_times.empty() ||
            !state.trigger_positions.empty() || !state.position_offsets.empty();
        if (control_enabled &&
            (state.trigger_modes.size() != control_size ||
             state.trigger_values.size() != control_size ||
             state.excitation_finished.size() != control_size ||
             state.stage_completed.size() != control_size ||
             state.trigger_times.size() != control_size ||
             state.trigger_positions.size() != control_size ||
             state.position_offsets.size() != layout.B))
            throw std::invalid_argument("GPU device control buffers do not match layout");
        const bool thermal_required = geometry.thermal_enabled ||
            thermal_mode == ThermalMode::Cpu || thermal_mode == ThermalMode::Gpu;
        if (thermal_required && state.temperatures.size() != layout.temperatures_size()) {
            throw std::invalid_argument("thermal-enabled GPU state requires temperatures");
        }
        if (!thermal_required && !state.temperatures.empty() &&
            state.temperatures.size() != layout.temperatures_size()) {
            throw std::invalid_argument("GPU temperature buffer does not match layout");
        }
        const auto finite_values = [](const std::vector<double>& values, const char* name) {
            for (const double value : values)
                if (!std::isfinite(value)) throw std::invalid_argument(name);
        };
        finite_values(state.currents, "non-finite current state");
        finite_values(state.m1, "non-finite mutual state");
        finite_values(state.dm1, "non-finite mutual gradient state");
        finite_values(state.velocity, "non-finite velocity state");
        finite_values(state.position, "non-finite position state");
        if (thermal_required) {
            finite_values(state.temperatures, "non-finite temperature state");
            finite_values(state.filament_masses, "invalid filament masses");
            finite_values(state.reference_resistances, "invalid reference resistances");
            for (std::size_t i = 0; i < state.filament_masses.size(); ++i)
                if (state.filament_masses[i] <= 0.0 || state.reference_resistances[i] <= 0.0 ||
                    (state.filament_materials[i] != 0 && state.filament_materials[i] != 1))
                    throw std::invalid_argument("invalid thermal material state");
        }
        return layout;
    }

    static void validate_mask(const std::vector<std::uint8_t>& mask,
                              std::size_t expected,
                              const char* message) {
        if (mask.size() != expected) throw std::invalid_argument(message);
    }

    static void validate_optional_mask(const std::vector<std::uint8_t>& mask,
                                       std::size_t expected,
                                       const char* message) {
        if (!mask.empty()) validate_mask(mask, expected, message);
    }

    static void validate_stage_voltages(const std::vector<double>& voltages,
                                        std::size_t expected) {
        if (!voltages.empty() && voltages.size() != expected)
            throw std::invalid_argument("GPU stage voltage buffer does not match batch-stage layout");
        for (const double voltage : voltages)
            if (!std::isfinite(voltage)) throw std::invalid_argument("non-finite GPU stage voltage");
    }

    void normalize_boundary_masks(GpuEngineState& state) const {
        const auto expected = layout_.B * layout_.S;
        if (state.stage_mask.empty()) state.stage_mask.assign(expected, 1);
        if (state.mutual_stage_mask.empty()) state.mutual_stage_mask.assign(expected, 1);
    }

    bool device_control_enabled() const noexcept {
        return state_.trigger_modes.size() == layout_.B * layout_.S &&
            state_.trigger_values.size() == layout_.B * layout_.S &&
            state_.excitation_finished.size() == layout_.B * layout_.S &&
            state_.stage_completed.size() == layout_.B * layout_.S &&
            state_.trigger_times.size() == layout_.B * layout_.S &&
            state_.trigger_positions.size() == layout_.B * layout_.S &&
            state_.position_offsets.size() == layout_.B;
    }

    void restore_inactive_state(const GpuEngineState& snapshot) {
        const auto B = layout_.B;
        const auto S = layout_.S;
        const auto F = layout_.F;
        const auto D = layout_.D;
        for (std::size_t b = 0; b < B; ++b) {
            if (state_.active_mask[b] != 0) continue;
            std::copy_n(snapshot.currents.data() + b * D, D, state_.currents.data() + b * D);
            std::copy_n(snapshot.m1.data() + b * S * F, S * F, state_.m1.data() + b * S * F);
            std::copy_n(snapshot.dm1.data() + b * S * F, S * F, state_.dm1.data() + b * S * F);
            if (snapshot.current_derivatives.size() == B * D &&
                state_.current_derivatives.size() == B * D)
                std::copy_n(snapshot.current_derivatives.data() + b * D, D,
                            state_.current_derivatives.data() + b * D);
            if (current_derivatives_.size() == B * D && snapshot.current_derivatives.size() == B * D)
                std::copy_n(snapshot.current_derivatives.data() + b * D, D,
                            current_derivatives_.data() + b * D);
            if (state_.temperatures.size() == B * F && snapshot.temperatures.size() == B * F)
                std::copy_n(snapshot.temperatures.data() + b * F, F, state_.temperatures.data() + b * F);
            if (state_.resistivities.size() == B * F && snapshot.resistivities.size() == B * F)
                std::copy_n(snapshot.resistivities.data() + b * F, F, state_.resistivities.data() + b * F);
            if (state_.resistances.size() == B * F && snapshot.resistances.size() == B * F)
                std::copy_n(snapshot.resistances.data() + b * F, F, state_.resistances.data() + b * F);
            if (state_.joule_energy.size() == B * F && snapshot.joule_energy.size() == B * F)
                std::copy_n(snapshot.joule_energy.data() + b * F, F, state_.joule_energy.data() + b * F);
            state_.velocity[b] = snapshot.velocity[b];
            state_.position[b] = snapshot.position[b];
        }
    }

    void ensure_running() const {
        if (is_shutdown()) {
            throw std::logic_error("GPU engine has been shut down");
        }
    }

    bool any_active() const noexcept {
        for (const auto active : state_.active_mask) {
            if (active != 0) {
                return true;
            }
        }
        return false;
    }

    double stage_resistance(std::size_t stage) const {
        return geometry_.stage_resistances.empty() ? 0.0 : geometry_.stage_resistances[stage];
    }

    double stage_inductance(std::size_t stage) const {
        if (!geometry_.stage_inductances.empty()) return geometry_.stage_inductances[stage];
        return physics::self_inductance(
            geometry_.stage_inner_radii[stage], geometry_.stage_outer_radii[stage],
            geometry_.stage_lengths[stage],
            static_cast<double>(geometry_.stage_turns[stage]) /
                ((geometry_.stage_outer_radii[stage] - geometry_.stage_inner_radii[stage]) *
                 geometry_.stage_lengths[stage]));
    }

    double filament_resistance(std::size_t filament) const {
        return geometry_.filament_resistances.empty() ? 0.0 : geometry_.filament_resistances[filament];
    }

    double filament_inductance(std::size_t filament) const {
        if (!geometry_.filament_inductances.empty()) return geometry_.filament_inductances[filament];
        return physics::self_inductance(
            geometry_.filament_inner_radii[filament], geometry_.filament_outer_radii[filament],
            geometry_.filament_lengths[filament],
            1.0 / ((geometry_.filament_outer_radii[filament] - geometry_.filament_inner_radii[filament]) *
                   geometry_.filament_lengths[filament]));
    }

    double stage_mutual(std::size_t first, std::size_t second) const {
        if (!geometry_.stage_mutual_inductances.empty())
            return geometry_.stage_mutual_inductances[first * layout_.S + second];
        return physics::mutual_inductance_coil(
            geometry_.stage_inner_radii[first], geometry_.stage_outer_radii[first],
            geometry_.stage_lengths[first], geometry_.stage_turns[first],
            geometry_.stage_inner_radii[second], geometry_.stage_outer_radii[second],
            geometry_.stage_lengths[second], geometry_.stage_turns[second],
            std::abs(geometry_.stage_positions[first] - geometry_.stage_positions[second]), 9, true);
    }

    double filament_mutual(std::size_t first, std::size_t second) const {
        if (!geometry_.filament_mutual_inductances.empty())
            return geometry_.filament_mutual_inductances[first * layout_.F + second];
        const auto radius = [](const GpuGeometryInput& g, std::size_t f) {
            return 0.5 * (g.filament_inner_radii[f] + g.filament_outer_radii[f]);
        };
        return physics::mutual_inductance_filament(
            radius(geometry_, first), radius(geometry_, second),
            std::abs(geometry_.filament_positions[first] - geometry_.filament_positions[second]), true);
    }

    void assemble_physical_system(bool use_cached_mutual = false) {
        const auto B = layout_.B;
        const auto S = layout_.S;
        const auto F = layout_.F;
        const auto D = layout_.D;
        std::fill(matrices_.begin(), matrices_.end(), 0.0);
        std::fill(rhs_.begin(), rhs_.end(), 0.0);
        for (std::size_t b = 0; b < B; ++b) {
            if (state_.active_mask[b] == 0) {
                for (std::size_t i = 0; i < D; ++i)
                    matrices_[b * D * D + i * D + i] = 1.0;
                continue;
            }
            const double velocity = state_.velocity[b];
            for (std::size_t s = 0; s < S; ++s) {
                if (stage_mask_[b * S + s] == 0 || state_.trigger_mask[b * S + s] == 0) {
                    matrices_[b * D * D + s * D + s] = 1.0;
                    continue;
                }
                const auto row = b * D * D + s * D;
                matrices_[row + s] = stage_inductance(s);
                rhs_[b * D + s] = (state_.stage_voltages.empty() ? 0.0 : state_.stage_voltages[b * S + s]) -
                    stage_resistance(s) * state_.currents[b * D + s];
                for (std::size_t other = 0; other < S; ++other) {
                    if (other != s && stage_mask_[b * S + s] != 0 &&
                        stage_mask_[b * S + other] != 0 &&
                        state_.trigger_mask[b * S + other] != 0)
                        matrices_[row + other] = stage_mutual(s, other);
                }
            }
            for (std::size_t f = 0; f < F; ++f) {
                const auto row = b * D * D + (S + f) * D;
                matrices_[row + S + f] = filament_inductance(f);
                const auto resistance_index = b * F + f;
                const double resistance = state_.resistances.size() == B * F &&
                    std::isfinite(state_.resistances[resistance_index]) &&
                    state_.resistances[resistance_index] > 0.0
                    ? state_.resistances[resistance_index] : filament_resistance(f);
                rhs_[b * D + S + f] = -resistance * state_.currents[b * D + S + f];
                for (std::size_t other = 0; other < F; ++other) {
                    if (other != f) matrices_[row + S + other] = filament_mutual(f, other);
                }
            }
            for (std::size_t s = 0; s < S; ++s) {
                for (std::size_t f = 0; f < F; ++f) {
                    const auto pair = (b * S + s) * F + f;
                    const double separation = geometry_.stage_positions[s] -
                        (geometry_.filament_positions[f] - state_.position[b]);
                    double m = 0.0;
                    double dm = 0.0;
                    const bool pair_active = mutual_stage_mask_[b * S + s] != 0 &&
                        state_.trigger_mask[b * S + s] != 0;
                    if (pair_active && use_cached_mutual) {
                        m = state_.m1[pair];
                        dm = state_.dm1[pair];
                    } else if (pair_active && !use_cached_mutual) {
                        m = physics::mutual_inductance_coil(
                            geometry_.stage_inner_radii[s], geometry_.stage_outer_radii[s],
                            geometry_.stage_lengths[s], geometry_.stage_turns[s],
                            geometry_.filament_inner_radii[f], geometry_.filament_outer_radii[f],
                            geometry_.filament_lengths[f], 1, separation, 9, true);
                        dm = physics::mutual_inductance_gradient_coil(
                            geometry_.stage_inner_radii[s], geometry_.stage_outer_radii[s],
                            geometry_.stage_lengths[s], geometry_.stage_turns[s],
                            geometry_.filament_inner_radii[f], geometry_.filament_outer_radii[f],
                            geometry_.filament_lengths[f], 1, separation, 9, true);
                    }
                    state_.m1[pair] = m;
                    state_.dm1[pair] = dm;
                    if (mutual_stage_mask_[b * S + s] != 0 && state_.trigger_mask[b * S + s] != 0) {
                        matrices_[b * D * D + s * D + S + f] = m;
                        matrices_[b * D * D + (S + f) * D + s] = m;
                        rhs_[b * D + s] -= velocity * dm * state_.currents[b * D + S + f];
                        rhs_[b * D + S + f] -= velocity * dm * state_.currents[b * D + s];
                    }
                }
            }
        }
    }

    void normalize_cpu_thermal() {
        if (policy_.thermal == ThermalMode::Gpu) {
            policy_.thermal = ThermalMode::Cpu;
            report_.thermal = ThermalMode::Cpu;
        }
    }

#if defined(COILGUN_CUDA_AVAILABLE)
    void initialize_runtime();
    void execute_solver_step();
    void execute_physical_pipeline();
    void execute_persistent_mutual();
    GpuGraphVariantKey graph_key() const;
    void calibrate_solver();
    void record_runtime_failure(const SolverStatus& status);
    void shutdown_runtime() noexcept;
    void release_runtime_resources() noexcept;
    void sync_runtime_state_after_reset();
#endif
    void execute_cpu_physical_pipeline() {
        const auto state_snapshot = state_;
        const auto matrices_snapshot = matrices_;
        const auto rhs_snapshot = rhs_;
        const auto solution_snapshot = solution_;
        const auto derivatives_snapshot = current_derivatives_;
        const auto pipeline_snapshot = pipeline_order_;
        const auto report_snapshot = report_;
        const auto policy_snapshot = policy_;
        try {
            normalize_cpu_thermal();
            pipeline_order_ = {PipelineStage::Mutual, PipelineStage::Matrix, PipelineStage::Solver,
                               PipelineStage::Force};
            const auto d = layout_.D;
            assemble_physical_system();
            std::vector<double> candidate_solution(solution_);
            double solver_time_ms = 0.0;
            for (std::size_t b = 0; b < layout_.B; ++b) {
                if (state_.active_mask[b] == 0) continue;
                Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix(
                    matrices_.data() + b * d * d, static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d));
                Eigen::Map<const Eigen::VectorXd> rhs(rhs_.data() + b * d, static_cast<Eigen::Index>(d));
                Eigen::Map<Eigen::VectorXd> solution(candidate_solution.data() + b * d,
                                                     static_cast<Eigen::Index>(d));
                const auto solver_start = std::chrono::steady_clock::now();
                Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> ldlt(matrix);
                if (ldlt.info() != Eigen::Success) {
                    throw std::runtime_error("Eigen LDLT factorization failed");
                }
                const auto diagonal = ldlt.vectorD();
                const double diagonal_scale = std::max(1.0, diagonal.cwiseAbs().maxCoeff());
                if (!diagonal.allFinite() ||
                    (diagonal.cwiseAbs().array() <=
                     std::numeric_limits<double>::epsilon() * diagonal_scale).any()) {
                    throw std::runtime_error("Eigen LDLT factorization is singular");
                }
                solution = ldlt.solve(rhs);
                const auto solver_stop = std::chrono::steady_clock::now();
                solver_time_ms +=
                    std::chrono::duration<double, std::milli>(solver_stop - solver_start).count();
                if (!solution.allFinite()) {
                    throw std::runtime_error("Eigen LDLT produced a non-finite solution");
                }
                double residual = 0.0;
                double scale = 1.0;
                for (Eigen::Index row = 0; row < static_cast<Eigen::Index>(d); ++row) {
                    double reconstructed = 0.0;
                    for (Eigen::Index col = 0; col < static_cast<Eigen::Index>(d); ++col)
                        reconstructed += matrix(row, col) * solution(col);
                    residual = std::max(residual, std::abs(reconstructed - rhs(row)));
                    scale = std::max(scale, std::abs(rhs(row)));
                }
                if (!std::isfinite(residual) || residual > 1.0e-10 * scale) {
                    throw std::runtime_error("Eigen LDLT residual exceeds sanity tolerance");
                }
            }
            solution_ = std::move(candidate_solution);
            report_.solver_time_ms += solver_time_ms;
            current_derivatives_ = solution_;
            state_.current_derivatives = solution_;
        for (std::size_t b = 0; b < layout_.B; ++b) {
            if (state_.active_mask[b] == 0) continue;
            double force = 0.0;
            for (std::size_t s = 0; s < layout_.S; ++s) {
                if (mutual_stage_mask_[b * layout_.S + s] == 0 || state_.trigger_mask[b * layout_.S + s] == 0) continue;
                for (std::size_t f = 0; f < layout_.F; ++f) {
                    force += state_.currents[b * d + s] * state_.currents[b * d + layout_.S + f] *
                             state_.dm1[(b * layout_.S + s) * layout_.F + f];
                }
            }
            state_.velocity[b] += state_.dt * force / state_.mass;
            state_.position[b] += state_.dt * (state_.velocity[b] - state_.dt * force / state_.mass);
            for (std::size_t i = 0; i < d; ++i)
                state_.currents[b * d + i] += state_.dt * current_derivatives_[b * d + i];
        }
        pipeline_order_.push_back(PipelineStage::State);
        if (policy_.thermal != ThermalMode::Disabled && !state_.temperatures.empty()) {
            const auto thermal_start = std::chrono::steady_clock::now();
            state_.resistivities.resize(layout_.B * layout_.F);
            state_.resistances.resize(layout_.B * layout_.F);
            state_.joule_energy.resize(layout_.B * layout_.F);
            for (std::size_t b = 0; b < layout_.B; ++b) {
                if (state_.active_mask[b] == 0) continue;
                for (std::size_t f = 0; f < layout_.F; ++f) {
                    const auto index = b * layout_.F + f;
                    const double current = state_.currents[b * d + layout_.S + f];
                    const double previous_resistance = state_.resistances[index] > 0.0
                        ? state_.resistances[index] : state_.reference_resistances[index];
                    const double q = current * current * previous_resistance * state_.dt;
                    const auto material = state_.filament_materials[index] == 1
                        ? physics::ArmatureMaterial::Copper : physics::ArmatureMaterial::Aluminum;
                    state_.temperatures[index] += q /
                        (state_.filament_masses[index] * physics::material_cp(material, state_.temperatures[index]));
                    const double beta = state_.filament_materials[index] == 1
                        ? physics::COPPER.temp_coefficient : physics::ALUMINUM.temp_coefficient;
                    state_.resistivities[index] = state_.filament_materials[index] == 1
                        ? physics::resistivity_copper(state_.temperatures[index])
                        : physics::resistivity_aluminum(state_.temperatures[index]);
                    state_.resistances[index] = state_.reference_resistances[index] *
                        (1.0 + beta * (state_.temperatures[index] - state_.reference_temperature));
                    state_.joule_energy[index] = q;
                }
            }
            const auto thermal_stop = std::chrono::steady_clock::now();
            report_.thermal_time_ms +=
                std::chrono::duration<double, std::milli>(thermal_stop - thermal_start).count();
            pipeline_order_.push_back(PipelineStage::Thermal);
        }
        restore_inactive_state(state_snapshot);
        } catch (...) {
            state_ = state_snapshot;
            matrices_ = matrices_snapshot;
            rhs_ = rhs_snapshot;
            solution_ = solution_snapshot;
            current_derivatives_ = derivatives_snapshot;
            pipeline_order_ = pipeline_snapshot;
            report_ = report_snapshot;
            policy_ = policy_snapshot;
            throw;
        }
    }

    std::unique_ptr<Resources> resources_;
    GpuGeometryInput geometry_;
    GpuStateLayout layout_;
    GpuEngineState initial_state_;
    GpuEngineState state_;
    GpuExecutionConfig config_;
    detail::GpuEngineFaultInjection fault_injection_;
    GpuExecutionPolicy policy_;
    ExecutionReport report_;
    GpuEngineResult result_;
    std::vector<std::uint8_t> stage_mask_;
    std::vector<std::uint8_t> mutual_stage_mask_;
    GpuGraphVariant variant_;
    GpuGraphVariant selected_variant_;
    bool variant_valid_ = false;
    std::size_t calibration_count_ = 0;
    std::size_t device_allocation_count_ = 0;
    std::vector<double> matrices_;
    std::vector<double> rhs_;
    std::vector<double> solution_;
    std::vector<double> current_derivatives_;
    std::vector<PipelineStage> pipeline_order_;
#if defined(COILGUN_CUDA_AVAILABLE)
    std::unique_ptr<GpuExecutionContext> context_;
    std::unique_ptr<GpuSolver> solver_;
    std::unique_ptr<GpuGraphCache> graph_cache_;
    GpuAdaptor persistent_adaptor_;
    PersistentBuffers persistent_buffers_;
    bool persistent_active_ = false;
    bool fallback_locked_ = false;
    MaterialTables material_tables_;
#endif
};

} // namespace coilgun::simulation::cuda
