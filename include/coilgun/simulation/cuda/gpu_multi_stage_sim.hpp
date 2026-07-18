/**
 * @file gpu_multi_stage_sim.hpp
 * @brief GPU-accelerated multi-stage coilgun simulation.
 * @author Winston Meursault
 *
 * Identical public API to MultiStageSim. Internally uses CUDA kernels
 * for the 4D mutual inductance integration and Eigen LDLT for the linear
 * system solve (identical to CPU path).
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
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace coilgun::simulation::cuda {

/**
 * @brief GPU-accelerated multi-stage coilgun simulator.
 *
 * Public API is identical to MultiStageSim. Internally, the 4D
 * Gauss-Legendre mutual inductance integration is offloaded to CUDA.
 * The linear system solve (Eigen LDLT) and kinematic/thermal updates
 * remain on CPU. Inter-coil mutual inductance is omitted for
 * simplicity.
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
     * @param backend GPU backend configuration.
     *
     * @throws std::invalid_argument if coils.size() != excitations.size()
     *         or trigger_configs.size() != coils.size()-1 or coils.size() > kMaxStages.
     */
    GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,
        components::Armature                      armature,
        std::vector<std::unique_ptr<Excitation>>  excitations,
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal = false,
        GpuOptLevel                               opt_level = GpuOptLevel::Full,
        const GpuBackend&                         backend = {});

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

private:
    void compute_M1_dM1();
    void build_system_matrix(const MultiStageState& s);
    void build_filament_M_matrix();
    MultiStageState compute_derivatives(const MultiStageState& s);
    double compute_force(const MultiStageState& s);
    void update_temperatures(MultiStageState& s, double dt_sub);
    void check_triggers();
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

    int N_fil_;
    Eigen::VectorXd R_diag_, L_diag_;
    Eigen::VectorXd R_fil_ref_, R_fil_, L_fil_, mass_fil_;
    Eigen::MatrixXd M_mat_;
    Eigen::MatrixXd M_cc_;

    std::vector<bool> triggered_, finished_;
    std::vector<double> trigger_times_;

    GpuAdaptor adaptor_;

    Eigen::MatrixXd M1_mat_, dM1_mat_;
    Eigen::MatrixXd L_total_;
    Eigen::VectorXd RHS_;

    MultiStageState state_;
    MultiStageResult result_;
    int step_count_ = 0;
    SP stepper_;
};

} // namespace coilgun::simulation::cuda
