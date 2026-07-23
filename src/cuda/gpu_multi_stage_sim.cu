/**
 * @file gpu_multi_stage_sim.cu
 * @brief GPU multi-stage compatibility wrapper over GpuEngine.
 */

#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"

#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <cstdint>

namespace coilgun::simulation::cuda {

namespace {

GpuGeometryInput make_geometry(const std::vector<components::DrivingCoil>& coils,
                               const components::Armature& armature,
                               bool thermal) {
    GpuGeometryInput geometry;
    const auto filament_count = static_cast<std::size_t>(armature.total_filaments());
    geometry.n_stages = coils.size();
    geometry.n_filaments = filament_count;
    geometry.thermal_enabled = thermal;
    geometry.stage_geometry.resize(coils.size());
    geometry.stage_inner_radii.resize(coils.size());
    geometry.stage_outer_radii.resize(coils.size());
    geometry.stage_lengths.resize(coils.size());
    geometry.stage_turns.resize(coils.size());
    geometry.stage_positions.resize(coils.size());
    geometry.stage_resistances.resize(coils.size());
    geometry.stage_inductances.resize(coils.size());
    for (std::size_t stage = 0; stage < coils.size(); ++stage) {
        const auto& coil = coils[stage];
        geometry.stage_geometry[stage] = coil.mean_radius();
        geometry.stage_inner_radii[stage] = coil.inner_radius();
        geometry.stage_outer_radii[stage] = coil.outer_radius();
        geometry.stage_lengths[stage] = coil.length();
        geometry.stage_turns[stage] = coil.turns();
        geometry.stage_positions[stage] = coil.position();
        geometry.stage_resistances[stage] = coil.resistance();
        geometry.stage_inductances[stage] = coil.self_inductance();
    }

    geometry.filament_geometry.resize(filament_count);
    geometry.filament_inner_radii.resize(filament_count);
    geometry.filament_outer_radii.resize(filament_count);
    geometry.filament_lengths.resize(filament_count);
    geometry.filament_positions.resize(filament_count);
    const auto filament_length = armature.length() / armature.axial_filaments();
    for (std::size_t k = 0; k < filament_count; ++k) {
        const int radial = static_cast<int>(k % armature.radial_filaments()) + 1;
        const int axial = static_cast<int>(k / armature.radial_filaments()) + 1;
        geometry.filament_geometry[k] = armature.filament_mean_radius(radial);
        geometry.filament_inner_radii[k] = armature.filament_inner_radius(radial);
        geometry.filament_outer_radii[k] = armature.filament_outer_radius(radial);
        geometry.filament_lengths[k] = filament_length;
        geometry.filament_positions[k] = armature.filament_axial_position(axial);
    }
    geometry.filament_resistances = armature.resistances();
    geometry.filament_inductances = armature.inductances();
    return geometry;
}

GpuEngineState make_state(const components::Armature& armature, double dt, bool thermal) {
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    GpuEngineState state;
    state.currents.assign(filaments, 0.0);
    state.m1.assign(filaments, 0.0);
    state.dm1.assign(filaments, 0.0);
    state.active_mask = {1};
    state.trigger_mask.assign(filaments, 1);
    state.velocity = {-armature.velocity()};
    state.position = {0.0};
    state.dt = dt;
    state.mass = armature.mass();
    if (thermal) {
        state.temperatures.assign(filaments, physics::T_REFERENCE);
        state.filament_masses = armature.masses();
        state.reference_resistances = armature.resistances();
        state.filament_materials.assign(filaments,
            armature.material() == physics::ArmatureMaterial::Copper ? 1 : 0);
        state.material_density = armature.material() == physics::ArmatureMaterial::Copper
            ? physics::COPPER.density : physics::ALUMINUM.density;
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
GpuMultiStageSim<SP>::GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,
        components::Armature                      armature,
        std::vector<std::unique_ptr<Excitation>>  excitations,
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal,
        GpuOptLevel                               opt_level,
        const GpuBackend&                         backend,
        BackendMode                               explicit_backend)
    : n_stages_(static_cast<int>(coils.size()))
    , coils_(std::move(coils))
    , armature_(std::move(armature))
    , excitations_(std::move(excitations))
    , trigger_configs_(std::move(trigger_configs))
    , dt_(dt)
    , enable_thermal_(enable_thermal)
    , opt_level_(opt_level)
    , backend_(backend)
    , explicit_backend_(explicit_backend)
    , N_fil_(armature_.total_filaments())
{
    if (n_stages_ <= 0)
        throw std::invalid_argument("GpuMultiStageSim: at least one stage is required");
    if (n_stages_ != static_cast<int>(excitations_.size()))
        throw std::invalid_argument("GpuMultiStageSim: coils.size() != excitations.size()");
    if (n_stages_ > kMaxStages)
        throw std::invalid_argument("GpuMultiStageSim: n_stages exceeds kMaxStages");
    if (static_cast<int>(trigger_configs_.size()) != n_stages_ - 1)
        throw std::invalid_argument(
            "GpuMultiStageSim: trigger_configs.size() must be n_stages-1");
    for (const auto& config : trigger_configs_)
        validate_trigger_config(config);
    if (!std::isfinite(dt_) || dt_ <= 0.0)
        throw std::invalid_argument("GpuMultiStageSim: dt must be positive and finite");
    for (const auto& excitation : excitations_) {
        if (!excitation || !std::isfinite(excitation->voltage()))
            throw std::invalid_argument("GpuMultiStageSim: invalid excitation");
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
        throw std::invalid_argument("GpuMultiStageSim: explicit backend mode is invalid");
    }

    triggered_.assign(n_stages_, false);
    excitation_finished_.assign(n_stages_, false);
    stage_completed_.assign(n_stages_, false);
    active_stage_mask_.assign(n_stages_, 0);
    mutual_stage_mask_.assign(n_stages_, 0);
    trigger_times_.assign(n_stages_, 0.0);
    trigger_positions_.assign(n_stages_, armature_.position());
    initial_stage_energies_.assign(n_stages_, 0.0);
    triggered_[0] = true;
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(excitations_[stage].get()))
            initial_stage_energies_[static_cast<std::size_t>(stage)] =
                0.5 * capacitor->capacitance() * capacitor->initial_voltage() *
                capacitor->initial_voltage();
    }

    auto initial_state = make_state(armature_, dt_, enable_thermal_);
    initial_state.currents.resize(static_cast<std::size_t>(n_stages_) + N_fil_, 0.0);
    initial_state.m1.assign(static_cast<std::size_t>(n_stages_) * N_fil_, 0.0);
    initial_state.dm1.assign(static_cast<std::size_t>(n_stages_) * N_fil_, 0.0);
    initial_state.trigger_mask.assign(n_stages_, 1);
    initial_state.stage_voltages.assign(n_stages_, 0.0);
    auto geometry = make_geometry(coils_, armature_, enable_thermal_);
    engine_ = std::make_unique<GpuEngine>(
        std::move(geometry), std::move(initial_state),
        make_config(opt_level_, backend_, enable_thermal_, explicit_backend_));
    sync_state_from_engine();
}

template<typename SP>
GpuMultiStageSim<SP>::~GpuMultiStageSim() = default;

template<typename SP>
void GpuMultiStageSim<SP>::sync_state_from_engine() {
    const auto& source = engine_->state();
    state_.currents = Eigen::Map<const Eigen::VectorXd>(
        source.currents.data(), static_cast<Eigen::Index>(source.currents.size()));
    state_.arm_position = armature_.position() - source.position[0];
    state_.arm_velocity = -source.velocity[0];
    if (enable_thermal_) {
        state_.filament_temperatures = Eigen::Map<const Eigen::VectorXd>(
            source.temperatures.data(), static_cast<Eigen::Index>(source.temperatures.size()));
    } else {
        state_.filament_temperatures.resize(0);
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::configure_engine_boundary() {
    active_stage_mask_.assign(static_cast<std::size_t>(n_stages_), 0);
    mutual_stage_mask_.assign(static_cast<std::size_t>(n_stages_), 0);
    std::vector<std::uint8_t> trigger_mask(static_cast<std::size_t>(n_stages_), 0);
    std::vector<double> voltages(static_cast<std::size_t>(n_stages_), 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        trigger_mask[static_cast<std::size_t>(stage)] = triggered_[stage] ? 1 : 0;
        if (triggered_[stage] && !stage_completed_[stage]) {
            active_stage_mask_[static_cast<std::size_t>(stage)] = 1;
            if (is_stage_within_range(stage))
                mutual_stage_mask_[static_cast<std::size_t>(stage)] = 1;
        }
        voltages[static_cast<std::size_t>(stage)] =
            triggered_[stage] && !excitation_finished_[stage]
                ? excitations_[stage]->voltage() : 0.0;
    }
    engine_->set_step_boundary_state(
        {1}, std::move(trigger_mask), active_stage_mask_, mutual_stage_mask_,
        std::move(voltages));

}

template<typename SP>
bool GpuMultiStageSim<SP>::is_stage_within_range(int stage_idx) const {
    if (opt_level_ == GpuOptLevel::Standard) return true;
    const auto& coil = coils_[static_cast<std::size_t>(stage_idx)];
    return std::abs(state_.arm_position - coil.position()) <= 10.0 * coil.length();
}

template<typename SP>
double GpuMultiStageSim<SP>::compute_force(
        const std::vector<double>& pre_step_currents) const {
    const auto stage_forces = compute_stage_forces(pre_step_currents);
    double force = 0.0;
    for (const double stage_force : stage_forces) force += stage_force;
    return force;
}

template<typename SP>
std::vector<double> GpuMultiStageSim<SP>::compute_stage_forces(
        const std::vector<double>& pre_step_currents) const {
    const auto& source = engine_->state();
    std::vector<double> stage_forces(static_cast<std::size_t>(n_stages_), 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!triggered_[stage] || stage_completed_[stage] ||
            mutual_stage_mask_[static_cast<std::size_t>(stage)] == 0) continue;
        for (int filament = 0; filament < N_fil_; ++filament) {
            const auto pair = static_cast<std::size_t>(stage) *
                              static_cast<std::size_t>(N_fil_) +
                              static_cast<std::size_t>(filament);
            // GpuEngine uses the opposite armature coordinate direction.
            stage_forces[static_cast<std::size_t>(stage)] -=
                pre_step_currents[static_cast<std::size_t>(stage)] *
                pre_step_currents[static_cast<std::size_t>(n_stages_ + filament)] *
                source.dm1[pair];
        }
    }
    return stage_forces;
}

template<typename SP>
std::vector<double> GpuMultiStageSim<SP>::compute_pre_step_gradients() const {
    std::vector<double> gradients(static_cast<std::size_t>(n_stages_) * N_fil_, 0.0);
    const double initial_center = armature_.position();
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!triggered_[stage] || stage_completed_[stage] ||
            mutual_stage_mask_[static_cast<std::size_t>(stage)] == 0)
            continue;
        for (int filament = 0; filament < N_fil_; ++filament) {
            const int radial = armature_.radial_filaments();
            const int axial = filament / radial + 1;
            const int radial_index = filament % radial + 1;
            const double relative =
                armature_.filament_axial_position(axial) - initial_center;
            const double separation = coils_[static_cast<std::size_t>(stage)].position() -
                (relative + state_.arm_position);
            const auto pair = static_cast<std::size_t>(stage) * N_fil_ +
                              static_cast<std::size_t>(filament);
            gradients[pair] = physics::mutual_inductance_gradient_coil(
                coils_[static_cast<std::size_t>(stage)].inner_radius(),
                coils_[static_cast<std::size_t>(stage)].outer_radius(),
                coils_[static_cast<std::size_t>(stage)].length(),
                coils_[static_cast<std::size_t>(stage)].turns(),
                armature_.filament_inner_radius(radial_index),
                armature_.filament_outer_radius(radial_index),
                armature_.length() / armature_.axial_filaments(), 1,
                separation, 9, false);
        }
    }
    return gradients;
}

template<typename SP>
std::vector<double> GpuMultiStageSim<SP>::compute_recorded_stage_forces(
        const std::vector<double>& gradients) const {
    const auto& source = engine_->state();
    std::vector<double> stage_forces(static_cast<std::size_t>(n_stages_), 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!triggered_[stage] || stage_completed_[stage] ||
            mutual_stage_mask_[static_cast<std::size_t>(stage)] == 0) continue;
        for (int filament = 0; filament < N_fil_; ++filament) {
            const auto pair = static_cast<std::size_t>(stage) *
                              static_cast<std::size_t>(N_fil_) +
                              static_cast<std::size_t>(filament);
            stage_forces[static_cast<std::size_t>(stage)] -=
                source.currents[static_cast<std::size_t>(stage)] *
                source.currents[static_cast<std::size_t>(n_stages_ + filament)] *
                gradients[pair];
        }
    }
    return stage_forces;
}

template<typename SP>
void GpuMultiStageSim<SP>::check_triggers(
        double pre_position, double post_position, double next_time) {
    const double current_time = next_time - dt_;
    for (int stage = 1; stage < n_stages_; ++stage) {
        if (triggered_[stage]) continue;
        if (!triggered_[stage - 1]) continue;
        const auto& config = trigger_configs_[static_cast<std::size_t>(stage - 1)];
        const auto index = static_cast<std::size_t>(stage);
        const bool fire = config.mode == TriggerMode::Position
            ? pre_position >= config.value ||
              (pre_position < config.value && post_position >= config.value)
            : next_time >= trigger_times_[static_cast<std::size_t>(stage - 1)] + config.value;
        if (fire) {
            triggered_[stage] = true;
            const double target = trigger_times_[static_cast<std::size_t>(stage - 1)] + config.value;
            trigger_times_[index] = config.mode == TriggerMode::Position
                ? current_time : std::max(current_time, target);
            trigger_positions_[index] = pre_position;
        }
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::extinguish_quiet_stages() {
    static constexpr double kQuietThreshold = 1e-6;
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!triggered_[stage] || stage_completed_[stage]) continue;
        if (excitations_[stage]->voltage() == 0.0 &&
            std::abs(state_.currents(stage)) < kQuietThreshold) {
            excitation_finished_[stage] = true;
            stage_completed_[stage] = true;
            state_.currents(stage) = 0.0;
            engine_->complete_stage(0, static_cast<std::size_t>(stage));
        }
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::record_step() {
    MultiStageStep entry;
    entry.state.time = step_count_ * dt_;
    entry.state.arm_position = state_.arm_position;
    entry.state.arm_velocity = state_.arm_velocity;
    entry.state.force = recorded_force_;
    entry.state.filament_currents.resize(N_fil_);
    for (int filament = 0; filament < N_fil_; ++filament)
        entry.state.filament_currents[filament] = state_.currents(n_stages_ + filament);
    if (enable_thermal_) {
        entry.state.filament_temperatures.resize(N_fil_);
        for (int filament = 0; filament < N_fil_; ++filament)
            entry.state.filament_temperatures[filament] = state_.filament_temperatures(filament);
    }
    entry.cap_voltages.resize(n_stages_);
    entry.coil_currents.resize(n_stages_);
    entry.stage_forces = recorded_stage_force_history_.back();
    for (int stage = 0; stage < n_stages_; ++stage) {
        entry.cap_voltages[stage] = triggered_[stage] ? excitations_[stage]->voltage() : 0.0;
        entry.coil_currents[stage] = triggered_[stage] ? state_.currents(stage) : 0.0;
    }
    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool GpuMultiStageSim<SP>::check_all_finished() const {
    // A finite trigger policy remains eligible even if earlier stages finished;
    // only completed stages and explicit +infinity policies are terminal.
    auto terminally_ineligible = [this](int stage_idx) {
        for (int stage = stage_idx; stage > 0; --stage) {
            const auto& config = trigger_configs_[static_cast<std::size_t>(stage - 1)];
            if (config.value == INFINITY) return true;
            if (triggered_[stage - 1]) return false;
        }
        return false;
    };
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (stage_completed_[stage]) continue;
        if (!triggered_[stage] && terminally_ineligible(stage)) continue;
        return false;
    }
    return true;
}

template<typename SP>
bool GpuMultiStageSim<SP>::check_termination(const TerminationPolicy& policy) {
    if (policy.enable_bound_check && state_.arm_position >= policy.barrel_end_position)
        return true;
    if (step_count_ >= policy.max_steps) return true;
    if (check_all_finished()) return true;
    for (int stage = 0; stage < n_stages_; ++stage)
        if (triggered_[stage] && !stage_completed_[stage] && !excitation_finished_[stage])
            return false;
    if (policy.enable_velocity_check && step_count_ >= policy.velocity_decay_steps) {
        const double acceleration = recorded_force_ / armature_.mass();
        const auto n = static_cast<int>(result_.history.size());
        bool decaying = true;
        for (int i = 0; i < policy.velocity_decay_steps; ++i) {
            if (n - 2 - i < 0 ||
                result_.history[n - 1 - i].state.arm_velocity >=
                    result_.history[n - 2 - i].state.arm_velocity) {
                decaying = false;
                break;
            }
        }
        if (decaying && std::abs(acceleration) < policy.accel_threshold) return true;
    }
    return false;
}

template<typename SP>
void GpuMultiStageSim<SP>::prepare_summary() {
    auto& summary = result_.summary;
    summary = MultiStageSummary{};
    summary.step_count = step_count_;
    if (result_.history.empty()) return;
    const auto& last = result_.history.back();
    summary.total_time = last.state.time;
    summary.muzzle_velocity = last.state.arm_velocity;

    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!triggered_[stage]) continue;
        PerStageSummary per_stage;
        per_stage.stage_index = stage;
        per_stage.trigger_time = trigger_times_[stage];
        per_stage.trigger_position = trigger_positions_[static_cast<std::size_t>(stage)];
        auto* capacitor = dynamic_cast<CapacitorExcitation*>(excitations_[stage].get());
        bool active = false;
        for (const auto& step : result_.history) {
            per_stage.peak_current = std::max(per_stage.peak_current,
                                              step.coil_currents[stage]);
            const auto history_index = static_cast<std::size_t>(
                &step - result_.history.data());
            if (history_index < recorded_stage_force_history_.size()) {
                per_stage.max_force = std::max(
                    per_stage.max_force,
                    std::abs(recorded_stage_force_history_[history_index][static_cast<std::size_t>(stage)]));
            }
            if (step.coil_currents[stage] > 1e-6) active = true;
            if (active) ++per_stage.step_count_active;
        }
        if (capacitor) {
            const double final_voltage = result_.history.back().cap_voltages[
                static_cast<std::size_t>(stage)];
            const double final_energy = 0.5 * capacitor->capacitance() *
                final_voltage * final_voltage;
            per_stage.energy_depleted =
                initial_stage_energies_[static_cast<std::size_t>(stage)] - final_energy;
        }
        summary.per_stage.push_back(per_stage);
    }

    for (const auto& step : result_.history) {
        summary.max_force = std::max(summary.max_force, std::abs(step.state.force));
        for (const double current : step.coil_currents)
            summary.peak_coil_current = std::max(summary.peak_coil_current, current);
    }
    double input_energy = 0.0;
    for (const auto& excitation : excitations_) {
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(excitation.get()))
            input_energy += 0.5 * capacitor->capacitance() *
                capacitor->initial_voltage() * capacitor->initial_voltage();
    }
    const double kinetic_energy = 0.5 * armature_.mass() *
        summary.muzzle_velocity * summary.muzzle_velocity;
    summary.efficiency = input_energy > 0.0 ? kinetic_energy / input_energy : 0.0;
}

template<typename SP>
const MultiStageStep& GpuMultiStageSim<SP>::step() {
    if constexpr (std::is_same_v<SP, RK4Stepper>)
        throw std::logic_error("RK4Stepper is not supported by GpuMultiStageSim");

    const double pre_position = state_.arm_position;
    check_triggers(pre_position, pre_position, step_count_ * dt_);
    extinguish_quiet_stages();
    const auto triggered_at_step_start = triggered_;
    configure_engine_boundary();
    const auto pre_step_currents = engine_->state().currents;
    const auto pre_step_gradients = compute_pre_step_gradients();
    engine_->step();
    const auto stage_forces = compute_stage_forces(pre_step_currents);
    applied_force_ = 0.0;
    for (const double stage_force : stage_forces) applied_force_ += stage_force;
    sync_state_from_engine();
    check_triggers(pre_position, state_.arm_position, (step_count_ + 1) * dt_);
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (triggered_at_step_start[stage] && !excitation_finished_[stage]) {
            excitations_[stage]->advance(
                dt_, pre_step_currents[static_cast<std::size_t>(stage)]);
            if (excitations_[stage]->finished()) excitation_finished_[stage] = true;
        }
        if (triggered_[stage] && excitation_finished_[stage] &&
            std::abs(state_.currents(stage)) < 1e-6) {
            stage_completed_[stage] = true;
            state_.currents(stage) = 0.0;
            engine_->complete_stage(0, static_cast<std::size_t>(stage));
        }
    }
    const auto recorded_stage_forces = compute_recorded_stage_forces(pre_step_gradients);
    recorded_force_ = 0.0;
    for (const double stage_force : recorded_stage_forces) recorded_force_ += stage_force;
    recorded_stage_force_history_.push_back(recorded_stage_forces);
    record_step();
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const MultiStageResult& GpuMultiStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const MultiStageResult& GpuMultiStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy)) step();
    prepare_summary();
    return result_;
}

template<typename SP>
void GpuMultiStageSim<SP>::reset() {
    engine_->reset();
    result_ = MultiStageResult{};
    step_count_ = 0;
    applied_force_ = 0.0;
    recorded_force_ = 0.0;
    recorded_stage_force_history_.clear();
    for (int stage = 0; stage < n_stages_; ++stage) {
        triggered_[stage] = stage == 0;
        excitation_finished_[stage] = false;
        stage_completed_[stage] = false;
        trigger_times_[stage] = 0.0;
        trigger_positions_[stage] = armature_.position();
        excitations_[stage]->reset();
    }
    initial_stage_energies_.assign(n_stages_, 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(excitations_[stage].get()))
            initial_stage_energies_[static_cast<std::size_t>(stage)] =
                0.5 * capacitor->capacitance() * capacitor->initial_voltage() *
                capacitor->initial_voltage();
    }
    sync_state_from_engine();
}

template class GpuMultiStageSim<EulerStepper>;
template class GpuMultiStageSim<RK4Stepper>;

} // namespace coilgun::simulation::cuda
