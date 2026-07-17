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
    Reference   = 0,   ///< All reference-grade: exact self-inductance, no distance cutoff, 9-pt GL.
    LookupTable = 1,   ///< T(q,p) lookup table for self-inductance only; other calcs reference-grade.
    Full        = 2    ///< All optimisations: table lookup + distance cutoff + adaptive GL order.
};

/**
 * @brief Full state vector for the multi-stage ODE system.
 *
 * Layout: @p currents[0..S-1] = coil currents @f$ I_{di} @f$ (S = n_stages),
 * @p currents[S..S+F-1] = filament currents @f$ I_{ij} @f$ (F = N_fil).
 */
struct MultiStageState {
    Eigen::VectorXd currents;              ///< [n_stages + N_fil] current vector, A.
    double arm_position = 0.0;             ///< Armature centre position, m.
    double arm_velocity = 0.0;             ///< Armature velocity, m/s.
    Eigen::VectorXd filament_temperatures; ///< [N_fil] filament temperatures, K.

    MultiStageState& operator+=(const MultiStageState& rhs);
    MultiStageState& operator*=(double scalar);
};

/// @related MultiStageState
MultiStageState operator+(MultiStageState lhs, const MultiStageState& rhs);
/// @related MultiStageState
MultiStageState operator*(double scalar, MultiStageState s);

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
    const MultiStageState&  state()   const { return state_; }
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
    MultiStageState compute_derivatives(const MultiStageState& s);
    void     build_filament_M_matrix();
    void     precompute_M_cc();
    bool     is_stage_within_range(int stage_idx) const;
    void     extinguish_quiet_stages();
    double   compute_force(const MultiStageState& s);
    void     update_temperatures(MultiStageState& s, double dt_sub);
    void     check_triggers();
    void     record_step();
    bool     check_all_finished() const;
    bool     check_termination(const TerminationPolicy& policy);
    void     prepare_summary();
    void     build_system_matrix(const MultiStageState& s);

    int n_stages_;

    std::vector<components::DrivingCoil> coils_;
    components::Armature armature_;

    std::vector<std::unique_ptr<Excitation>> excitations_;
    std::vector<TriggerConfig> trigger_configs_;
    std::vector<bool> triggered_;
    std::vector<bool> finished_;
    std::vector<double> trigger_times_;

    double dt_;
    bool   enable_thermal_;
    OptimizationLevel opt_level_;

    int N_fil_;

    Eigen::VectorXd R_diag_;       // [n_stages]
    Eigen::VectorXd L_diag_;       // [n_stages]
    Eigen::MatrixXd M_cc_;         // [n_stages x n_stages], inter-coil mutual

    Eigen::VectorXd L_fil_;        // [N_fil]
    Eigen::VectorXd R_fil_ref_;    // [N_fil]
    Eigen::VectorXd R_fil_;        // [N_fil]
    Eigen::VectorXd mass_fil_;     // [N_fil]
    Eigen::MatrixXd M_mat_;        // [N_fil x N_fil], inter-filament mutual

    Eigen::MatrixXd M1_mat_;       // [n_stages x N_fil], coil<->filament M
    Eigen::MatrixXd dM1_mat_;      // [n_stages x N_fil], coil<->filament dM/dx

    Eigen::MatrixXd L_total_;      // [(n_stages+N_fil) x (n_stages+N_fil)]
    Eigen::VectorXd RHS_;          // [n_stages + N_fil]

    MultiStageState  state_;
    MultiStageResult result_;
    int              step_count_ = 0;
    StepperPolicy    stepper_;
};

} // namespace coilgun::simulation
