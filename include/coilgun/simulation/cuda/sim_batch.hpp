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
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace coilgun::simulation::cuda {

/**
 * @brief Batch container for GPU-accelerated parameter sweeps.
 *
 * All simulations share identical coil and armature geometry.
 * Use set_excitations() to configure per-simulation excitation
 * and trigger parameters. Batch GPU kernel launches amortise
 * kernel launch overhead across all simulations.
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
             const GpuBackend&                     backend = {});

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
    void run();

    /// Run all simulations with a custom termination policy.
    void run(const TerminationPolicy& policy);

    /// @brief Result for a single simulation.
    const MultiStageResult& result(int sim_id) const;

    /// @brief Number of simulations in the batch.
    int num_sims() const { return num_sims_; }

private:
    struct SimInstance {
        std::vector<std::unique_ptr<Excitation>> excitations;
        std::vector<TriggerConfig>               trigger_configs;
        std::vector<bool>                        triggered;
        std::vector<bool>                        finished;
        std::vector<double>                      trigger_times;
        MultiStageState                          state;
        Eigen::VectorXd                          R_fil;
        Eigen::VectorXd                          R_fil_ref;
        MultiStageResult                         result;
        int                                      step_count = 0;
        bool                                     all_finished = false;
    };

    void step_all();
    void compute_all_M1_dM1();
    void build_system_matrix(const SimInstance& sim,
                              Eigen::MatrixXd& L,
                              const Eigen::MatrixXd& M1) const;
    void solve_and_update(SimInstance& sim);
    double compute_force(const SimInstance& sim,
                         const Eigen::MatrixXd& dM) const;
    void check_triggers(SimInstance& sim);
    void extinguish_quiet_stages(SimInstance& sim);
    void record_step(SimInstance& sim);
    bool check_termination(SimInstance& sim, const TerminationPolicy& policy) const;
    void prepare_summary(SimInstance& sim);

    int n_stages_;
    int N_fil_;
    int num_sims_;
    double dt_;
    GpuBackend backend_;
    GpuAdaptor adaptor_;

    std::vector<components::DrivingCoil> coils_;
    components::Armature armature_;

    Eigen::VectorXd R_diag_;
    Eigen::VectorXd L_diag_;
    Eigen::VectorXd L_fil_;
    Eigen::VectorXd mass_fil_;
    Eigen::MatrixXd M_mat_;
    Eigen::MatrixXd M_cc_;

    std::vector<SimInstance> sims_;
    Eigen::MatrixXd batch_M1_;
    Eigen::MatrixXd batch_dM1_;
    Eigen::MatrixXd L_total_;
    Eigen::VectorXd RHS_;
    StepperPolicy stepper_;
};

} // namespace coilgun::simulation::cuda
