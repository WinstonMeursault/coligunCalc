/**
 * @file gpu_single_stage_sim.hpp
 * @brief GPU-accelerated single-stage coilgun simulation.
 * @author Winston Meursault
 *
 * Identical public API to SingleStageSim. Internally uses CUDA kernels
 * for the 4D mutual inductance integration and Eigen LDLT for the linear
 * system solve (identical to CPU path).
 */

#pragma once

#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/time_stepper.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include <Eigen/Dense>
#include <memory>

namespace coilgun::simulation::cuda {

/**
 * @brief GPU-accelerated single-stage coilgun simulator.
 *
 * Public API is identical to SingleStageSim. Internally, the 4D
 * Gauss-Legendre mutual inductance integration is offloaded to CUDA.
 * The linear system solve (Eigen LDLT) and kinematic/thermal updates
 * remain on CPU.
 *
 * @tparam StepperPolicy Time-stepping policy (EulerStepper or RK4Stepper).
 *
 * @see SingleStageSim
 */
template<typename StepperPolicy = EulerStepper>
class GpuSingleStageSim {
public:
    /**
     * @brief Construct and initialise a GPU-accelerated single-stage simulation.
     * @param coil Driving coil geometry (copied).
     * @param armature Armature geometry and discretisation (copied).
     * @param excitation Excitation source (moved-in).
     * @param dt Fixed time step, s.
     * @param enable_thermal Enable adiabatic filament heating.
     * @param opt_level GPU optimisation level (default: Full).
     * @param backend GPU backend configuration.
     */
    GpuSingleStageSim(
        components::DrivingCoil                coil,
        components::Armature                   armature,
        std::unique_ptr<Excitation>           excitation,
        double                                 dt,
        bool                                   enable_thermal = false,
        GpuOptLevel                            opt_level = GpuOptLevel::Full,
        const GpuBackend&                      backend = {});

    ~GpuSingleStageSim();

    GpuSingleStageSim(const GpuSingleStageSim&) = delete;
    GpuSingleStageSim& operator=(const GpuSingleStageSim&) = delete;

    /**
     * @brief Advance one time step.
     * @return Reference to the newly recorded SimStep.
     */
    const SimStep& step();

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

    /**
     * @brief Reset to initial conditions.
     *
     * Restores all currents to zero, armature position/velocity to their
     * initial values, resets the excitation, and clears the result history.
     */
    void reset();

    /// @brief Result after the last run().
    const SimResult& result()  const { return result_; }
    /// @brief Current internal state.
    const SimState&  state()   const { return state_; }
    /// @brief Fixed time step, s.
    double  dt()          const { return dt_; }
    /// @brief Step count since construction or last reset().
    int     step_count()  const { return step_count_; }

private:
    /// @brief Compute M1 and dM1 matrices using GPU kernel.
    void compute_M1_dM1();

    SimState compute_derivatives(const SimState& s);
    void     build_filament_M_matrix();
    double   compute_force(const SimState& s);
    void     update_temperatures(SimState& s, double dt_sub);
    void     record_step(double cap_voltage);
    bool     check_termination(const TerminationPolicy& policy);
    void     prepare_summary();

    components::DrivingCoil coil_;
    components::Armature    armature_;
    std::unique_ptr<Excitation> excitation_;
    double dt_;
    bool   enable_thermal_;
    GpuOptLevel opt_level_;
    GpuBackend  backend_;

    int N_fil_;
    double R_d_, L_d_;

    Eigen::MatrixXd M_mat_;
    Eigen::VectorXd L_fil_;
    Eigen::VectorXd R_fil_ref_;
    Eigen::VectorXd R_fil_;
    Eigen::VectorXd mass_fil_;

    GpuAdaptor adaptor_;

    Eigen::MatrixXd M1_mat_;
    Eigen::MatrixXd dM1_mat_;

    Eigen::MatrixXd L_total_;
    Eigen::VectorXd RHS_;

    SimState  state_;
    SimResult result_;
    int       step_count_ = 0;
    StepperPolicy stepper_;
};

} // namespace coilgun::simulation::cuda
