/**
 * @file gpu_single_stage_sim.cu
 * @brief GPU single-stage compatibility wrapper over GpuEngine.
 */

#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"

#include "coilgun/physics/constants.hpp"

#include <cuda_runtime_api.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <stdexcept>

namespace coilgun::simulation::cuda {

namespace {

GpuGeometryInput make_geometry(const components::DrivingCoil& coil,
                               const components::Armature& armature,
                               bool thermal) {
    GpuGeometryInput geometry;
    const auto filament_count = static_cast<std::size_t>(armature.total_filaments());
    geometry.n_stages = 1;
    geometry.n_filaments = filament_count;
    geometry.thermal_enabled = thermal;
    geometry.stage_geometry = {coil.mean_radius()};
    geometry.stage_inner_radii = {coil.inner_radius()};
    geometry.stage_outer_radii = {coil.outer_radius()};
    geometry.stage_lengths = {coil.length()};
    geometry.stage_turns = {coil.turns()};
    geometry.stage_positions = {coil.position()};
    geometry.filament_geometry.resize(filament_count);
    geometry.filament_inner_radii.resize(filament_count);
    geometry.filament_outer_radii.resize(filament_count);
    geometry.filament_lengths.resize(filament_count);
    geometry.filament_positions.resize(filament_count);
    for (std::size_t k = 0; k < filament_count; ++k) {
        const int radial = static_cast<int>(k % armature.radial_filaments()) + 1;
        const int axial = static_cast<int>(k / armature.radial_filaments()) + 1;
        geometry.filament_geometry[k] = armature.filament_mean_radius(radial);
        geometry.filament_inner_radii[k] = armature.filament_inner_radius(radial);
        geometry.filament_outer_radii[k] = armature.filament_outer_radius(radial);
        geometry.filament_lengths[k] = armature.length() / armature.axial_filaments();
        geometry.filament_positions[k] = armature.filament_axial_position(axial);
    }
    geometry.stage_resistances = {coil.resistance()};
    geometry.stage_inductances = {coil.self_inductance()};
    geometry.filament_resistances = armature.resistances();
    geometry.filament_inductances = armature.inductances();
    return geometry;
}

GpuEngineState make_state(const components::Armature& armature, double dt, bool thermal) {
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    GpuEngineState state;
    state.currents.assign(filaments + 1, 0.0);
    state.m1.assign(filaments, 0.0);
    state.dm1.assign(filaments, 0.0);
    state.active_mask = {1};
    state.trigger_mask = {1};
    state.velocity = {-armature.velocity()};
    state.position = {0.0};
    state.stage_voltages = {0.0};
    state.dt = dt;
    state.mass = armature.mass();
    if (thermal) {
        state.temperatures.assign(filaments, physics::T_REFERENCE);
        state.filament_masses = armature.masses();
        state.reference_resistances = armature.resistances();
        state.filament_materials.assign(filaments,
            armature.material() == physics::ArmatureMaterial::Copper ? 1 : 0);
        state.material_density = physics::ALUMINUM.density;
        if (armature.material() == physics::ArmatureMaterial::Copper)
            state.material_density = physics::COPPER.density;
    }
    return state;
}

BackendMode effective_backend(const GpuBackend& backend, BackendMode explicit_backend) {
    if (explicit_backend != BackendMode::Auto) return explicit_backend;
    if (backend.backend != BackendMode::Auto) return backend.backend;
    return backend.use_persistent ? BackendMode::Persistent : BackendMode::Direct;
}

GpuExecutionConfig make_config(GpuOptLevel opt_level, const GpuBackend& backend,
                               bool thermal, BackendMode explicit_backend) {
    switch (opt_level) {
    case GpuOptLevel::Standard:
    case GpuOptLevel::Full:
    case GpuOptLevel::Aggressive:
        break;
    default:
        throw std::invalid_argument("GPU optimization level is invalid");
    }
    GpuExecutionConfig config;
    config.precision = static_cast<PrecisionMode>(static_cast<int>(opt_level));
    config.thermal = thermal ? ThermalMode::Cpu : ThermalMode::Disabled;
    config.backend = effective_backend(backend, explicit_backend);
    config.solver = SolverMode::Eigen;
    config.deterministic = config.backend != BackendMode::Persistent;
    config.device_id = backend.device_id;
    config.threads_per_block = backend.threads_per_block;
    config.enable_profiling = backend.enable_profiling;
    return config;
}

} // namespace

template<typename SP>
GpuSingleStageSim<SP>::GpuSingleStageSim(
    components::DrivingCoil coil,
    components::Armature armature,
    std::unique_ptr<Excitation> excitation,
    double dt,
    bool enable_thermal,
    GpuOptLevel opt_level,
    const GpuBackend& backend,
    BackendMode explicit_backend)
    : coil_(std::move(coil)), armature_(std::move(armature)),
      excitation_(std::move(excitation)), dt_(dt), enable_thermal_(enable_thermal),
      opt_level_(opt_level), backend_(backend), explicit_backend_(explicit_backend) {
    if (!excitation_) throw std::invalid_argument("GPU excitation must not be null");
    if (!std::isfinite(dt_) || dt_ <= 0.0) {
        throw std::invalid_argument("GPU dt must be positive and finite");
    }
    backend_.validate();
    switch (explicit_backend_) {
    case BackendMode::Auto:
    case BackendMode::Graph:
    case BackendMode::Persistent:
    case BackendMode::Fallback:
    case BackendMode::Direct:
        break;
    default:
        throw std::invalid_argument("GPU single-stage explicit backend mode is invalid");
    }
    const double voltage = excitation_->voltage();
    if (!std::isfinite(voltage)) throw std::invalid_argument("GPU excitation voltage must be finite");
    auto geometry = make_geometry(coil_, armature_, enable_thermal_);
    auto initial_state = make_state(armature_, dt_, enable_thermal_);
    engine_ = std::make_unique<GpuEngine>(
        std::move(geometry), std::move(initial_state),
        make_config(opt_level_, backend_, enable_thermal_, explicit_backend_));
    sync_state_from_engine();
}

template<typename SP>
void GpuSingleStageSim<SP>::sync_state_from_engine() {
    const auto& source = engine_->state();
    state_.currents = Eigen::Map<const Eigen::VectorXd>(source.currents.data(),
                                                        static_cast<Eigen::Index>(source.currents.size()));
    state_.arm_position = armature_.position() - source.position[0];
    state_.arm_velocity = -source.velocity[0];
    if (enable_thermal_)
        state_.filament_temperatures = Eigen::Map<const Eigen::VectorXd>(
            source.temperatures.data(), static_cast<Eigen::Index>(source.temperatures.size()));
    else
        state_.filament_temperatures.resize(0);
}

template<typename SP>
double GpuSingleStageSim<SP>::compute_force() const {
    const auto& source = engine_->state();
    double force = 0.0;
    for (std::size_t k = 0; k < engine_->layout().F; ++k)
        force -= source.currents[k + 1] * source.currents[0] * source.dm1[k];
    return force;
}

template<typename SP>
void GpuSingleStageSim<SP>::record_step() {
    SimStep entry;
    entry.time = step_count_ * dt_;
    entry.cap_voltage = excitation_->voltage();
    entry.coil_current = state_.currents(0);
    entry.arm_position = state_.arm_position;
    entry.arm_velocity = state_.arm_velocity;
    entry.force = compute_force();
    entry.filament_currents.resize(engine_->layout().F);
    for (std::size_t k = 0; k < engine_->layout().F; ++k)
        entry.filament_currents[k] = state_.currents(static_cast<Eigen::Index>(k + 1));
    if (enable_thermal_)
        entry.filament_temperatures.assign(state_.filament_temperatures.data(),
                                           state_.filament_temperatures.data() + state_.filament_temperatures.size());
    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool GpuSingleStageSim<SP>::check_termination(const TerminationPolicy& policy) const {
    if (policy.enable_bound_check && state_.arm_position >= policy.barrel_end_position) return true;
    if (step_count_ >= policy.max_steps) return true;
    if (!policy.enable_velocity_check || step_count_ < policy.velocity_decay_steps) return false;
    const auto n = static_cast<int>(result_.history.size());
    bool decaying = true;
    for (int i = 0; i < policy.velocity_decay_steps; ++i) {
        if (n - 2 - i < 0 || result_.history[n - 1 - i].arm_velocity >= result_.history[n - 2 - i].arm_velocity) {
            decaying = false;
            break;
        }
    }
    return decaying && std::abs(compute_force() / armature_.mass()) < policy.accel_threshold;
}

template<typename SP>
const SimStep& GpuSingleStageSim<SP>::step() {
    if constexpr (std::is_same_v<SP, RK4Stepper>) {
        throw std::logic_error("RK4Stepper is not supported by GpuSingleStageSim");
    }
    const bool within_mutual_range = opt_level_ == GpuOptLevel::Standard ||
        std::abs(state_.arm_position - coil_.position()) <= 10.0 * coil_.length();
    engine_->set_mutual_stage_mask({static_cast<std::uint8_t>(within_mutual_range ? 1 : 0)});
    engine_->set_stage_voltage(0, excitation_->voltage());
    engine_->step();
    sync_state_from_engine();
    excitation_->advance(dt_, state_.currents(0));
    record_step();
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const SimResult& GpuSingleStageSim<SP>::run() { return run(TerminationPolicy::defaults()); }

template<typename SP>
const SimResult& GpuSingleStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy) && !excitation_->finished()) step();
    prepare_summary();
    return result_;
}

template<typename SP>
void GpuSingleStageSim<SP>::prepare_summary() {
    auto& summary = result_.summary;
    summary.step_count = step_count_;
    if (result_.history.empty()) return;
    summary.total_time = result_.history.back().time;
    summary.muzzle_velocity = result_.history.back().arm_velocity;
    for (const auto& entry : result_.history) {
        summary.max_force = std::max(summary.max_force, std::abs(entry.force));
        summary.peak_coil_current = std::max(summary.peak_coil_current, entry.coil_current);
    }
    if (const auto* cap = dynamic_cast<const CapacitorExcitation*>(excitation_.get())) {
        const double input = 0.5 * cap->capacitance() * cap->initial_voltage() * cap->initial_voltage();
        const double output = 0.5 * armature_.mass() * summary.muzzle_velocity * summary.muzzle_velocity;
        summary.efficiency = input > 0.0 ? output / input : 0.0;
    }
}

template<typename SP>
void GpuSingleStageSim<SP>::reset() {
    engine_->reset();
    excitation_->reset();
    result_ = {};
    step_count_ = 0;
    sync_state_from_engine();
}

template class GpuSingleStageSim<EulerStepper>;
template class GpuSingleStageSim<RK4Stepper>;

} // namespace coilgun::simulation::cuda
