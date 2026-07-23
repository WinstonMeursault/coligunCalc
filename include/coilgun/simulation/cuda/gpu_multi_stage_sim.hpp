/**
 * @file gpu_multi_stage_sim.hpp
 * @brief GPU-accelerated multi-stage coilgun simulation.
 * @author Winston Meursault
 *
 * Compatible core stepping API to MultiStageSim. Physical state advancement
 * is delegated to GpuEngine using its normalized B=1, S-stage layout. The
 * optional explicit backend mode disambiguates migrated backend selection
 * from the legacy use_persistent flag.
 */

#pragma once

#include "coilgun/simulation/multi_stage_result.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/time_stepper.hpp"
#include "coilgun/simulation/trigger_config.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <vector>

#include "coilgun/simulation/cuda/persistent_kernel.cuh"

namespace coilgun::simulation::cuda {

inline GpuBackend multi_stage_default_backend() {
    GpuBackend backend;
    backend.use_persistent = false;
    backend.backend = BackendMode::Direct;
    return backend;
}

/**
 * @brief GPU-accelerated multi-stage coilgun simulator.
 *
 * The core stepping surface follows MultiStageSim. Trigger, excitation,
 * result, and termination policy state remains in this compatibility wrapper while
 * GpuEngine owns the normalized physical state and execution policy.
 *
 * @tparam SP Time-stepping policy (EulerStepper or RK4Stepper).
 *
 * @see MultiStageSim
 */
template<typename SP = EulerStepper>
class GpuMultiStageSim {
public:
    /// @brief Maximum number of stages that can be configured (safety cap).
    static constexpr int kMaxStages = 50;

    /**
     * @brief Construct and initialise a GPU-accelerated multi-stage simulation.
     *
     * @param coils Driving coil geometries (one per stage, copied).
     * @param armature Armature geometry and discretisation (copied).
     * @param excitations Excitation sources (one per stage, moved-in).
     * @param trigger_configs Trigger configs for stages 1..N-1 (stage 0 auto-triggers at t=0).
     * @param dt Fixed time step, s.
     * @param enable_thermal Enable adiabatic filament heating.
     * @param opt_level GPU optimisation level (default: Full).
      * @param backend GPU backend configuration. The omitted/default backend
      *         selects Direct (or its documented CUDA fallback). The legacy
      *         `use_persistent` flag is consulted only when `backend.backend`
      *         is `Auto`; an explicit backend mode takes precedence.
     * @param explicit_backend Explicit backend selection. When not Auto it
     *         takes precedence over `use_persistent` and `backend.backend`.
     *
     * @throws std::invalid_argument if the stage/excitation/trigger counts
     *         are inconsistent, stage count exceeds kMaxStages, dt is not
     *         positive and finite, an excitation is null or has a non-finite
     *         voltage, or the backend configuration is invalid.
     */
    GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,
        components::Armature                      armature,
        std::vector<std::unique_ptr<Excitation>>  excitations,
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal = false,
         GpuOptLevel                               opt_level = GpuOptLevel::Full,
          const GpuBackend&                         backend = multi_stage_default_backend(),
         BackendMode                               explicit_backend = BackendMode::Auto);

    ~GpuMultiStageSim();

    GpuMultiStageSim(const GpuMultiStageSim&) = delete;
    GpuMultiStageSim& operator=(const GpuMultiStageSim&) = delete;

    /**
     * @brief Advance one time step.
     * @return Reference to the newly recorded MultiStageStep.
     */
    const MultiStageStep& step();

    /**
     * @brief Run to completion using the default termination policy.
     * @return Simulation result (history + summary).
     */
    const MultiStageResult& run();

    /**
     * @brief Run to completion with a custom termination policy.
     * @param policy Termination criteria.
     * @return Simulation result.
     */
    const MultiStageResult& run(const TerminationPolicy& policy);

    /**
     * @brief Reset to initial conditions.
     *
     * Restores all currents to zero, armature position/velocity, resets
     * all excitations, clears the result history.
     */
    void reset();

    /// @brief Result after the last run().
    const MultiStageResult& result()  const { return result_; }
    /// @brief Current internal state.
    const MultiStageState&  state()   const { return state_; }
    /// @brief Fixed time step, s.
    double  dt()          const { return dt_; }
    /// @brief Step count since construction or last reset().
    int     step_count()  const { return step_count_; }
    /// @brief Number of configured stages.
    int     num_stages()  const { return n_stages_; }
     /// @brief Resolved backend and execution diagnostics.
     ///
    /// For Graph, the report describes graph-assisted execution of the
    /// mutual-inductance segment only. Matrix assembly, solve, force/state
    /// orchestration, and thermal updates remain outside the captured graph.
    const ExecutionReport& execution_report() const { return engine_->report(); }
    /// @brief True only after a complete CUDA-backed step used the partial Graph path.
    bool graph_assisted() const {
        const auto& report = execution_report();
        return report.backend == BackendMode::Graph && report.gpu_executed;
    }
    /// @brief Current per-filament resistances in thermal mode.
     std::vector<double> filament_resistances() const { return engine_->state().resistances; }
     /// @brief Snapshot of stage/filament mutual inductances, row-major.
     std::vector<double> mutual_inductances() const { return engine_->state().m1; }
     /// @brief Snapshot of stage/filament dM/dx values in engine coordinates.
     std::vector<double> mutual_gradients() const { return engine_->state().dm1; }

private:
    void sync_state_from_engine();
    void configure_engine_boundary();
    double compute_force(const std::vector<double>& pre_step_currents) const;
    std::vector<double> compute_stage_forces(const std::vector<double>& pre_step_currents) const;
    std::vector<double> compute_pre_step_gradients() const;
    std::vector<double> compute_recorded_stage_forces(
        const std::vector<double>& gradients) const;
    void check_triggers(double pre_position, double post_position, double next_time);
    void extinguish_quiet_stages();
    void record_step();
    bool check_all_finished() const;
    bool check_termination(const TerminationPolicy& policy);
    void prepare_summary();
    bool is_stage_within_range(int stage_idx) const;

    int n_stages_;
    std::vector<components::DrivingCoil> coils_;
    components::Armature armature_;
    std::vector<std::unique_ptr<Excitation>> excitations_;
    std::vector<TriggerConfig> trigger_configs_;
    double dt_;
    bool enable_thermal_;
    GpuOptLevel opt_level_;
    GpuBackend backend_;
    BackendMode explicit_backend_;
    int N_fil_;

    std::vector<bool> triggered_, excitation_finished_, stage_completed_;
    std::vector<double> trigger_times_;
    std::vector<double> trigger_positions_;
    std::vector<double> initial_stage_energies_;
    std::vector<std::uint8_t> active_stage_mask_;
    std::vector<std::uint8_t> mutual_stage_mask_;
    std::vector<std::vector<double>> recorded_stage_force_history_;

    std::unique_ptr<GpuEngine> engine_;
    MultiStageState state_;
    MultiStageResult result_;
    int step_count_ = 0;
    double applied_force_ = 0.0;
    double recorded_force_ = 0.0;
    SP stepper_;
};

} // namespace coilgun::simulation::cuda
