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
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include <Eigen/Dense>
#include <memory>

namespace coilgun::simulation {

/**
 * @brief Internal state vector for the ODE system.
 *
 * Layout: @p currents[0] = driving coil current @f$ I_d @f$,
 * @p currents[1..N_fil] = filament currents @f$ I_{ij} @f$.
 */
struct SimState {
    Eigen::VectorXd currents;              ///< [N_fil+1] current vector, A.
    double arm_position = 0.0;             ///< Armature centre position, m.
    double arm_velocity = 0.0;             ///< Armature velocity, m/s.
    Eigen::VectorXd filament_temperatures; ///< [N_fil] filament temperatures, K.

    SimState& operator+=(const SimState& rhs);
    SimState& operator*=(double scalar);
};

/// @related SimState
SimState operator+(SimState lhs, const SimState& rhs);
/// @related SimState
SimState operator*(double scalar, SimState s);

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
    const SimState&  state()   const { return state_; }
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
    SimState compute_derivatives(const SimState& s);
    void     build_filament_M_matrix();
    double   compute_force(const SimState& s);
    void     update_temperatures(SimState& s, double dt_sub);
    void     record_step(double cap_voltage);
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
    Eigen::VectorXd R_fil_;
    Eigen::VectorXd mass_fil_;

    Eigen::VectorXd M1_;
    Eigen::VectorXd dM1_;

    Eigen::MatrixXd L_total_;
    Eigen::VectorXd RHS_;

    SimState  state_;
    SimResult result_;
    int       step_count_ = 0;
    StepperPolicy stepper_;
};

} // namespace coilgun::simulation
