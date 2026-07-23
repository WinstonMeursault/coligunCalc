/**
 * @file single_stage_sim.hpp
 * @brief Single-stage synchronous induction coilgun simulation engine.
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.3.2, Sec.5, Sec.8.
 */

#pragma once

#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/sim_result.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/time_stepper.hpp"
#include "coilgun/simulation/integration_state.hpp"
#include "coilgun/simulation/derivative_workspace.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include <Eigen/Dense>
#include <memory>

namespace coilgun::simulation {

/**
 * @brief Single-stage synchronous induction coilgun simulation engine.
 *
 * @tparam StepperPolicy Time-integration policy (EulerStepper or RK4Stepper).
 *
 * Solves the coupled circuit ODE (§3.2) and motion equations (§5.3)
 * for one driving coil + one discretised armature. Supports capacitor,
 * crowbar-diode, and arbitrary-waveform excitation.
 *
 * @see NumericalModel Sec.3.2, Sec.5, Sec.8.
 */
template<typename StepperPolicy = EulerStepper>
class SingleStageSim {
public:
    /**
     * @brief Construct and initialise a single-stage simulation.
     * @param coil Driving coil geometry (copied).
     * @param armature Armature geometry and discretisation (copied).
     * @param excitation Excitation source (moved-in).
     * @param dt Fixed time step, s.
     * @param enable_thermal Enable adiabatic filament heating.
     */
    SingleStageSim(
        components::DrivingCoil coil,
        components::Armature    armature,
        std::unique_ptr<Excitation> excitation,
        double                     dt,
        bool                       enable_thermal = false);

    /**
     * @brief Advance one time step.
     * @return Reference to the newly recorded SimStep.
     */
    const SimStep&   step();

    /**
     * @brief Run to completion using the default termination policy.
     * @return Simulation result (history + summary).
     */
    const SimResult& run();

    /**
     * @brief Run to completion with a custom termination policy.
     * @param policy Termination criteria.
     * @return Simulation result.
     */
    const SimResult& run(const TerminationPolicy& policy);

    /// @brief Result after the last run().
    const SimResult& result()  const { return result_; }
    /// @brief Current internal state.
    const SimState& state() const { return integration_state_.physical; }
    const StageRuntimeState& stage_state() const { return integration_state_.stages.front(); }
    /// @brief Fixed time step, s.
    double  dt()          const { return dt_; }
    /// @brief Step count since construction or last reset().
    int     step_count()  const { return step_count_; }

    /**
     * @brief Reset to initial conditions.
     *
     * Restores all currents to zero, armature position/velocity to their
     * initial values, resets the excitation, and clears the result history.
     */
    void    reset();

private:
    DerivativeResult evaluate_derivatives(const SimState& state,
                                          const ExcitationSnapshot& excitation,
                                          DerivativeWorkspace& workspace,
                                          bool circuit_active) const;
    IntegrationState advance_euler(const IntegrationState& pre, double dt);
    IntegrationState advance_rk4_segment(const IntegrationState& pre, double dt);
    IntegrationState advance_rk4_event_aware(const IntegrationState& pre, double dt);
    IntegrationState make_trial(const IntegrationState& initial,
                                const DerivativeResult& derivative,
                                double scale) const;
    void apply_excitation_events(IntegrationState& state, double coil_current) const;
    double excitation_event_value(const ExcitationSnapshot& snapshot) const;
    double compute_force(const SimState& state,
                         const Eigen::VectorXd& mutual_gradient) const;
    Eigen::VectorXd derive_resistance(const SimState& state) const;
    void     build_filament_M_matrix();
    void     record_step(double post_time);
    bool     check_termination(const TerminationPolicy& policy);
    void     prepare_summary();

    components::DrivingCoil coil_;
    components::Armature    armature_;
    std::unique_ptr<Excitation>    excitation_;
    double dt_;
    bool   enable_thermal_;

    int    N_fil_;
    double R_d_, L_d_;

    Eigen::MatrixXd M_mat_;
    Eigen::VectorXd L_fil_;
    Eigen::VectorXd R_fil_ref_;
    Eigen::VectorXd mass_fil_;

    IntegrationState integration_state_;
    IntegrationState initial_integration_state_;
    DerivativeWorkspace workspace_;
    Eigen::VectorXd pre_step_mutual_gradient_;
    SimResult result_;
    int       step_count_ = 0;
    StepperPolicy stepper_;
};

} // namespace coilgun::simulation
