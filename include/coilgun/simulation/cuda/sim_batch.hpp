/**
 * @file sim_batch.hpp
 * @brief SimBatch — batch simulation container for GPU-accelerated parameter sweeps.
 * @author Winston Meursault
 *
 * All simulations in a batch share the same driving coil and armature
 * geometry. Only excitation parameters (voltage, capacitance) and
 * trigger configurations may differ between simulations.
 */

#pragma once

#include "coilgun/simulation/multi_stage_result.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/trigger_config.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

#include "coilgun/simulation/cuda/persistent_kernel.cuh"

namespace coilgun::simulation::cuda {

/**
 * @brief Batch container for GPU-accelerated parameter sweeps.
 *
 * All simulations share identical coil and armature geometry.
 * Use set_excitations() to configure per-simulation excitation
 * and trigger parameters. Batch GPU kernel launches amortise
 * kernel launch overhead across all simulations. The wrapper owns one
 * GpuEngine with a fixed row-major [B][S][F] state layout; inactive rows are
 * frozen in place rather than compacted.
 *
 * GpuEngine-backed SimBatch supports EulerStepper only. RK4Stepper is
 * explicitly rejected: construction and run() throw std::logic_error rather
 * than silently executing Euler.
 *
 * @tparam StepperPolicy Time-stepping policy (EulerStepper or RK4Stepper).
 */
template<typename StepperPolicy = EulerStepper>
class SimBatch {
public:
    static constexpr int kMaxStages = 50;

    /**
     * @brief Construct a batch of geometrically identical simulations.
     * @param coils Shared driving coil geometry (copied).
     * @param armature Shared armature geometry (copied).
     * @param num_sims Number of concurrent simulations.
     * @param dt Fixed time step, s.
     * @param backend GPU backend configuration.
     */
    SimBatch(std::vector<components::DrivingCoil> coils,
             components::Armature                  armature,
             int                                   num_sims,
             double                                dt,
             const GpuBackend&                     backend = {},
             BackendMode                           explicit_backend = BackendMode::Auto);

    ~SimBatch();
    SimBatch(const SimBatch&) = delete;
    SimBatch& operator=(const SimBatch&) = delete;

    /**
     * @brief Configure excitation and triggers for one simulation.
     * @param sim_id Simulation index (0..num_sims-1).
     * @param excitations Per-stage excitation sources (moved-in).
     * @param trigger_configs Trigger configs (size = n_stages-1).
     */
    void set_excitations(int sim_id,
        std::vector<std::unique_ptr<Excitation>> excitations,
        std::vector<TriggerConfig>               trigger_configs);

     /// Run all simulations to their respective termination criteria.
     /// RK4Stepper instantiations throw std::logic_error.
    void run();

     /// Run all simulations with a custom termination policy.
     /// RK4Stepper instantiations throw std::logic_error.
    void run(const TerminationPolicy& policy);

    /// @brief Result for a single simulation.
    const MultiStageResult& result(int sim_id) const;

    /// @brief Number of simulations in the batch.
    int num_sims() const { return num_sims_; }

    /// @brief Resolved execution diagnostics for the shared engine.
    ///
    /// The report retains requested/resolved backend and solver modes,
    /// fallback reasons, and whether CUDA actually executed a step.
    const ExecutionReport& execution_report() const { return engine_->report(); }

    /// @brief Whether the mutual-inductance graph path has executed.
    bool graph_assisted() const {
        const auto& report = execution_report();
        return report.backend == BackendMode::Graph && report.gpu_executed;
    }

    /// @brief Snapshot of row-major [B][S][F] mutual-inductance gradients.
    std::vector<double> mutual_gradients() const { return engine_->state().dm1; }

private:
    struct SimInstance {
        std::vector<std::unique_ptr<Excitation>> excitations;
        std::vector<TriggerConfig>               trigger_configs;
        std::vector<bool>                        triggered;
        std::vector<bool>                        excitation_finished;
        std::vector<bool>                        stage_completed;
        std::vector<double>                      trigger_times;
        std::vector<double>                      trigger_positions;
        std::vector<double>                      initial_stage_energies;
        bool                                     configured = false;
        MultiStageState                          state;
        MultiStageResult                         result;
        std::vector<std::vector<double>>         stage_force_history;
        int                                      step_count = 0;
        bool                                     active = true;
    };

    void configure_engine_boundary();
    void sync_states_from_engine();
    void check_triggers(SimInstance& sim);
    void extinguish_quiet_stages(SimInstance& sim);
    std::vector<double> compute_stage_forces(std::size_t sim_id,
                                             const std::vector<double>& currents,
                                             const std::vector<double>& gradients) const;
    void record_step(std::size_t sim_id,
                     const std::vector<double>& stage_forces);
    bool check_termination(SimInstance& sim, const TerminationPolicy& policy) const;
    void prepare_summary(SimInstance& sim);
    bool terminally_ineligible(const SimInstance& sim, int stage) const;
    bool is_stage_within_range(int stage, const SimInstance& sim) const;

    int n_stages_;
    int N_fil_;
    int num_sims_;
    double dt_;
    GpuBackend backend_;
    BackendMode explicit_backend_;
    std::vector<components::DrivingCoil> coils_;
    components::Armature armature_;

    std::vector<SimInstance> sims_;
    std::unique_ptr<GpuEngine> engine_;
};

} // namespace coilgun::simulation::cuda
