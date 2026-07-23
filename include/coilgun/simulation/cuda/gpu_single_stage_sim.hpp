/**
 * @file gpu_single_stage_sim.hpp
 * @brief GPU-accelerated single-stage coilgun simulation.
 */

#pragma once

#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"

#include <memory>

namespace coilgun::simulation::cuda {

inline GpuBackend single_stage_default_backend() {
    GpuBackend backend;
    backend.use_persistent = false;
    backend.backend = BackendMode::Graph;
    return backend;
}

template<typename StepperPolicy = EulerStepper>
class GpuSingleStageSim {
public:
    // GpuEngine currently exposes a committed-step API, not the four staged
    // evaluations required by RK4. RK4 instantiations therefore fail loudly
    // in step() instead of silently changing the requested integration method.
    GpuSingleStageSim(
        components::DrivingCoil coil,
        components::Armature armature,
        std::unique_ptr<Excitation> excitation,
        double dt,
         bool enable_thermal = false,
         GpuOptLevel opt_level = GpuOptLevel::Full,
         const GpuBackend& backend = single_stage_default_backend(),
         BackendMode explicit_backend = BackendMode::Auto);

    ~GpuSingleStageSim() = default;

    GpuSingleStageSim(const GpuSingleStageSim&) = delete;
    GpuSingleStageSim& operator=(const GpuSingleStageSim&) = delete;

    const SimStep& step();
    const SimResult& run();
    const SimResult& run(const TerminationPolicy& policy);
    void reset();

    const SimResult& result() const { return result_; }
    const SimState& state() const { return state_; }
    double dt() const { return dt_; }
    int step_count() const { return step_count_; }
    // The wrapper owns the engine and its CUDA/CPU resources through RAII.
    // The moved-in excitation remains owned by this wrapper; no raw device
    // pointer is exposed by this interface.
    const ExecutionReport& execution_report() const { return engine_->report(); }
    std::vector<double> filament_resistances() const { return engine_->state().resistances; }

private:
    void sync_state_from_engine();
    double compute_force() const;
    double compute_force_at(double position, const Eigen::VectorXd& currents) const;
    void record_step(double force);
    bool check_termination(const TerminationPolicy& policy) const;
    void prepare_summary();

    components::DrivingCoil coil_;
    components::Armature armature_;
    std::unique_ptr<Excitation> excitation_;
    double dt_;
    bool enable_thermal_;
    GpuOptLevel opt_level_;
    GpuBackend backend_;
    BackendMode explicit_backend_;
    std::unique_ptr<GpuEngine> engine_;
    SimState state_;
    SimResult result_;
    int step_count_ = 0;
};

} // namespace coilgun::simulation::cuda
