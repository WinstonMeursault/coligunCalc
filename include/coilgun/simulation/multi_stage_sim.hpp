/**
 * @file multi_stage_sim.hpp
 * @brief Multi-stage synchronous induction coilgun simulation engine.
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.3.3, Sec.3.4, Sec.5, Sec.8.
 */

#pragma once

#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/multi_stage_result.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/time_stepper.hpp"
#include "coilgun/simulation/trigger_config.hpp"
#include "coilgun/simulation/integration_state.hpp"
#include "coilgun/simulation/derivative_workspace.hpp"
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace coilgun::simulation {

/**
 * @brief Computational optimisation levels for the multi-stage simulator.
 *
 * Levels are cumulative: each level includes all optimisations of lower levels.
 */
enum class OptimizationLevel {
    Reference   = 0,   ///< Fixed 9-point mutual-inductance quadrature, no distance cutoff.
    LookupTable = 1,   ///< Reserved for the component-level T(q,p) lookup selection; same runtime path as Reference.
    Full        = 2    ///< Distance cutoff plus adaptive 9/4-point mutual-inductance quadrature.
};

/**
 * @brief Multi-stage synchronous induction coilgun simulation engine.
 *
 * @tparam StepperPolicy Time-integration policy (EulerStepper or RK4Stepper).
 *
 * Solves the coupled circuit ODE (NumericalModel Eq.3.12) for @em n driving
 * coils + one discretised armature.  Supports capacitor, crowbar-diode, and
 * arbitrary-waveform excitation per stage.  Stages are triggered sequentially
 * by user-specified position or time-delay criteria.
 *
 * The system matrix has dimension @f$ S + F @f$ where @f$ S @f$ is the number
 * of stages and @f$ F @f$ the number of armature filaments.  Inactive stages
 * (not yet triggered or extinguished) have their rows/columns set to identity
 * to keep the matrix non-singular.
 *
 * @see NumericalModel Sec.3.3, Sec.5, Sec.8.
 */
template<typename StepperPolicy = EulerStepper>
class MultiStageSim {
public:
    /// @brief Maximum number of stages that can be configured (safety cap).
    static constexpr int kMaxStages = 50;

    /**
     * @brief Construct and initialise a multi-stage simulation.
     *
     * @param coils Driving coil geometries (one per stage, copied).
     * @param armature Armature geometry and discretisation (copied).
     * @param excitations Excitation sources (one per stage, moved-in).
     * @param trigger_configs Trigger configs for stages 1..N-1 (stage 0 auto-triggers at t=0).
     * @param dt Fixed time step, s.
     * @param enable_thermal Enable adiabatic filament heating.
     * @param opt_level Optimisation level (default: Reference).
     *
     * @throws std::invalid_argument if coils.size() != excitations.size()
     *         or trigger_configs.size() != coils.size()-1 or coils.size() > kMaxStages.
     */
    MultiStageSim(
        std::vector<components::DrivingCoil>        coils,
        components::Armature                         armature,
        std::vector<std::unique_ptr<Excitation>>     excitations,
        std::vector<TriggerConfig>                   trigger_configs,
        double                                       dt,
        bool                                         enable_thermal = false,
        OptimizationLevel                            opt_level = OptimizationLevel::Reference);

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

    /// @brief Result after the last run().
    const MultiStageResult& result()  const { return result_; }
    /// @brief Current internal state.
    const MultiStageState& state() const { return integration_state_.physical; }
    const StageRuntimeState& stage_state(std::size_t stage) const;
    bool circuit_active(std::size_t stage) const noexcept;
    /// @brief Fixed time step, s.
    double  dt()          const { return dt_; }
    /// @brief Step count since construction or last reset().
    int     step_count()  const { return step_count_; }
    /// @brief Number of configured stages.
    int     num_stages()  const { return n_stages_; }

    /**
     * @brief Reset to initial conditions.
     *
     * Restores all currents to zero, armature position/velocity, resets
     * all excitations, clears the result history.
     */
    void    reset();

private:
    DerivativeResult evaluate_derivatives(const IntegrationState& state,
                                          DerivativeWorkspace& workspace) const;
    IntegrationState advance_euler(const IntegrationState& pre, double dt);
    IntegrationState advance_rk4_segment(const IntegrationState& pre, double dt);
    IntegrationState advance_rk4_event_aware(const IntegrationState& pre, double dt);
    IntegrationState make_trial(const IntegrationState& initial,
                                const DerivativeResult& derivative,
                                double scale) const;
    Eigen::VectorXd derive_resistance(const MultiStageState& state) const;
    void apply_boundary_events(IntegrationState& state, const IntegrationState& pre,
                               double segment_time, double absolute_time) const;
    void check_triggers(IntegrationState& state, const IntegrationState& pre,
                        double absolute_time) const;
    void complete_quiet_stages(IntegrationState& state) const;
    void     build_filament_M_matrix();
    void     precompute_M_cc();
    bool     is_stage_within_range(int stage_idx) const;
    double   compute_force(const MultiStageState& state,
                           const Eigen::MatrixXd& mutual_gradient,
                           const std::vector<StageRuntimeState>& stages) const;
    void     record_step(double post_time);
    bool     check_all_finished() const;
    bool     check_termination(const TerminationPolicy& policy);
    void     prepare_summary();

    int n_stages_;

    std::vector<components::DrivingCoil> coils_;
    components::Armature armature_;

    std::vector<std::unique_ptr<Excitation>> excitations_;
    std::vector<TriggerConfig> trigger_configs_;
    std::vector<double> initial_stage_energies_;

    double dt_;
    bool   enable_thermal_;
    OptimizationLevel opt_level_;

    int N_fil_;

    Eigen::VectorXd R_diag_;       // [n_stages]
    Eigen::VectorXd L_diag_;       // [n_stages]
    Eigen::MatrixXd M_cc_;         // [n_stages x n_stages], inter-coil mutual

    Eigen::VectorXd L_fil_;        // [N_fil]
    Eigen::VectorXd R_fil_ref_;    // [N_fil]
    Eigen::VectorXd mass_fil_;     // [N_fil]
    Eigen::MatrixXd M_mat_;        // [N_fil x N_fil], inter-filament mutual

    IntegrationState integration_state_;
    IntegrationState initial_integration_state_;
    DerivativeWorkspace workspace_;
    Eigen::VectorXd pre_step_mutual_gradient_;
    MultiStageResult result_;
    int              step_count_ = 0;
    StepperPolicy    stepper_;
};

} // namespace coilgun::simulation
