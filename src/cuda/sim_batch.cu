/**
 * @file sim_batch.cu
 * @brief Batched compatibility wrapper over the unified GPU engine.
 */

#include "coilgun/simulation/cuda/sim_batch.hpp"

#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/simulation/trigger_config.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace coilgun::simulation::cuda {

namespace {

GpuGeometryInput make_geometry(const std::vector<components::DrivingCoil>& coils,
                               const components::Armature& armature) {
    GpuGeometryInput geometry;
    const auto stages = coils.size();
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    geometry.n_stages = stages;
    geometry.n_filaments = filaments;
    geometry.stage_geometry.resize(stages);
    geometry.stage_inner_radii.resize(stages);
    geometry.stage_outer_radii.resize(stages);
    geometry.stage_lengths.resize(stages);
    geometry.stage_turns.resize(stages);
    geometry.stage_positions.resize(stages);
    geometry.stage_resistances.resize(stages);
    geometry.stage_inductances.resize(stages);
    for (std::size_t stage = 0; stage < stages; ++stage) {
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

    geometry.filament_geometry.resize(filaments);
    geometry.filament_inner_radii.resize(filaments);
    geometry.filament_outer_radii.resize(filaments);
    geometry.filament_lengths.resize(filaments);
    geometry.filament_positions.resize(filaments);
    const auto filament_length = armature.length() / armature.axial_filaments();
    for (std::size_t filament = 0; filament < filaments; ++filament) {
        const int radial = static_cast<int>(filament % armature.radial_filaments()) + 1;
        const int axial = static_cast<int>(filament / armature.radial_filaments()) + 1;
        geometry.filament_geometry[filament] = armature.filament_mean_radius(radial);
        geometry.filament_inner_radii[filament] = armature.filament_inner_radius(radial);
        geometry.filament_outer_radii[filament] = armature.filament_outer_radius(radial);
        geometry.filament_lengths[filament] = filament_length;
        geometry.filament_positions[filament] = armature.filament_axial_position(axial);
    }
    geometry.filament_resistances = armature.resistances();
    geometry.filament_inductances = armature.inductances();

    geometry.stage_mutual_inductances.assign(stages * stages, 0.0);
    for (std::size_t first = 0; first < stages; ++first) {
        for (std::size_t second = first + 1; second < stages; ++second) {
            const double mutual = physics::mutual_inductance_coil(
                geometry.stage_inner_radii[first], geometry.stage_outer_radii[first],
                geometry.stage_lengths[first], geometry.stage_turns[first],
                geometry.stage_inner_radii[second], geometry.stage_outer_radii[second],
                geometry.stage_lengths[second], geometry.stage_turns[second],
                std::abs(geometry.stage_positions[first] - geometry.stage_positions[second]), 9, true);
            geometry.stage_mutual_inductances[first * stages + second] = mutual;
            geometry.stage_mutual_inductances[second * stages + first] = mutual;
        }
    }

    geometry.filament_mutual_inductances.assign(filaments * filaments, 0.0);
    for (std::size_t first = 0; first < filaments; ++first) {
        const double first_radius = geometry.filament_geometry[first];
        for (std::size_t second = first + 1; second < filaments; ++second) {
            const double mutual = physics::mutual_inductance_filament(
                first_radius, geometry.filament_geometry[second],
                std::abs(geometry.filament_positions[first] - geometry.filament_positions[second]), true);
            geometry.filament_mutual_inductances[first * filaments + second] = mutual;
            geometry.filament_mutual_inductances[second * filaments + first] = mutual;
        }
    }
    return geometry;
}

GpuEngineState make_state(const components::Armature& armature, int num_sims,
                          int stages, double dt) {
    const auto filaments = static_cast<std::size_t>(armature.total_filaments());
    const auto batch = static_cast<std::size_t>(num_sims);
    const auto stage_count = static_cast<std::size_t>(stages);
    const auto dimension = stage_count + filaments;
    GpuEngineState state;
    state.currents.assign(batch * dimension, 0.0);
    state.m1.assign(batch * stage_count * filaments, 0.0);
    state.dm1.assign(batch * stage_count * filaments, 0.0);
    state.active_mask.assign(batch, 1);
    state.trigger_mask.assign(batch * stage_count, 0);
    state.stage_mask.assign(batch * stage_count, 0);
    state.mutual_stage_mask.assign(batch * stage_count, 0);
    state.stage_voltages.assign(batch * stage_count, 0.0);
    state.velocity.assign(batch, -armature.velocity());
    state.position.assign(batch, 0.0);
    for (int sim = 0; sim < num_sims; ++sim) {
        state.trigger_mask[static_cast<std::size_t>(sim) * stage_count] = 1;
        state.stage_mask[static_cast<std::size_t>(sim) * stage_count] = 1;
        state.mutual_stage_mask[static_cast<std::size_t>(sim) * stage_count] = 1;
    }
    state.dt = dt;
    state.mass = armature.mass();
    return state;
}

BackendMode effective_backend(const GpuBackend& backend, BackendMode explicit_backend) {
    if (explicit_backend != BackendMode::Auto) return explicit_backend;
    if (backend.backend != BackendMode::Auto) return backend.backend;
    return backend.use_persistent ? BackendMode::Persistent : BackendMode::Direct;
}

GpuExecutionConfig make_config(const GpuBackend& backend, BackendMode explicit_backend) {
    GpuExecutionConfig config;
    const auto resolved_backend = effective_backend(backend, explicit_backend);
    config.backend = resolved_backend;
    config.solver = SolverMode::Auto;
    config.precision = PrecisionMode::Full;
    config.deterministic = resolved_backend != BackendMode::Persistent;
    config.device_id = backend.device_id;
    config.threads_per_block = backend.threads_per_block;
    config.enable_profiling = backend.enable_profiling;
    return config;
}

GpuCapability make_capability(const GpuBackend& backend, int num_sims,
                              BackendMode explicit_backend) {
    GpuCapability capability;
    if (effective_backend(backend, explicit_backend) == BackendMode::Persistent && num_sims > 1)
        capability.supports_persistent = false;
    return capability;
}

} // namespace

template<typename SP>
SimBatch<SP>::SimBatch(std::vector<components::DrivingCoil> coils,
                       components::Armature armature, int num_sims, double dt,
                       const GpuBackend& backend, BackendMode explicit_backend)
    : n_stages_(static_cast<int>(coils.size()))
    , N_fil_(armature.total_filaments())
    , num_sims_(num_sims)
    , dt_(dt)
    , backend_(backend)
    , explicit_backend_(explicit_backend)
    , coils_(std::move(coils))
    , armature_(std::move(armature)) {
    if constexpr (std::is_same_v<SP, RK4Stepper>)
        throw std::logic_error("RK4Stepper is not supported by SimBatch");
    if (n_stages_ <= 0) throw std::invalid_argument("SimBatch: at least one coil stage is required");
    if (num_sims_ <= 0) throw std::invalid_argument("SimBatch: num_sims must be positive");
    if (static_cast<std::size_t>(num_sims_) > backend_.max_batch_sims)
        throw std::invalid_argument("SimBatch: num_sims exceeds backend.max_batch_sims");
    if (n_stages_ > kMaxStages) throw std::invalid_argument("SimBatch: n_stages exceeds kMaxStages");
    if (!std::isfinite(dt_) || dt_ <= 0.0) throw std::invalid_argument("SimBatch: dt must be positive and finite");
    backend_.validate();
    switch (explicit_backend_) {
    case BackendMode::Auto:
    case BackendMode::Graph:
    case BackendMode::Persistent:
    case BackendMode::Fallback:
    case BackendMode::Direct:
        break;
    default:
        throw std::invalid_argument("SimBatch: explicit backend mode is invalid");
    }

    auto state = make_state(armature_, num_sims_, n_stages_, dt_);
    engine_ = std::make_unique<GpuEngine>(make_geometry(coils_, armature_), std::move(state),
                                           make_config(backend_, explicit_backend_),
                                           make_capability(backend_, num_sims_, explicit_backend_));
    sims_.resize(static_cast<std::size_t>(num_sims_));
    for (auto& sim : sims_) {
        sim.triggered.assign(static_cast<std::size_t>(n_stages_), false);
        sim.excitation_finished.assign(static_cast<std::size_t>(n_stages_), false);
        sim.stage_completed.assign(static_cast<std::size_t>(n_stages_), false);
        sim.triggered[0] = true;
        sim.trigger_times.assign(static_cast<std::size_t>(n_stages_), 0.0);
        sim.trigger_positions.assign(static_cast<std::size_t>(n_stages_), armature_.position());
        sim.initial_stage_energies.assign(static_cast<std::size_t>(n_stages_), 0.0);
        sim.active = true;
    }
    sync_states_from_engine();
}

template<typename SP>
SimBatch<SP>::~SimBatch() = default;

template<typename SP>
void SimBatch<SP>::set_excitations(int sim_id,
    std::vector<std::unique_ptr<Excitation>> excitations,
    std::vector<TriggerConfig> trigger_configs) {
    if (sim_id < 0 || sim_id >= num_sims_) throw std::out_of_range("sim_id out of range");
    if (static_cast<int>(excitations.size()) != n_stages_)
        throw std::invalid_argument("SimBatch: excitations size must equal stage count");
    if (static_cast<int>(trigger_configs.size()) != n_stages_ - 1)
        throw std::invalid_argument("SimBatch: trigger_configs size must equal stage count - 1");
    for (const auto& excitation : excitations) {
        if (!excitation || !std::isfinite(excitation->voltage()))
            throw std::invalid_argument("SimBatch: invalid excitation");
    }
    for (const auto& config : trigger_configs) validate_trigger_config(config);

    auto& sim = sims_[static_cast<std::size_t>(sim_id)];
    sim.excitations = std::move(excitations);
    sim.trigger_configs = std::move(trigger_configs);
    sim.triggered.assign(static_cast<std::size_t>(n_stages_), false);
    sim.excitation_finished.assign(static_cast<std::size_t>(n_stages_), false);
    sim.stage_completed.assign(static_cast<std::size_t>(n_stages_), false);
    sim.triggered[0] = true;
    sim.trigger_times.assign(static_cast<std::size_t>(n_stages_), 0.0);
    sim.trigger_positions.assign(static_cast<std::size_t>(n_stages_), armature_.position());
    sim.initial_stage_energies.assign(static_cast<std::size_t>(n_stages_), 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(sim.excitations[stage].get()))
            sim.initial_stage_energies[static_cast<std::size_t>(stage)] =
                0.5 * capacitor->capacitance() * capacitor->initial_voltage() * capacitor->initial_voltage();
    }
    sim.result = MultiStageResult{};
    sim.stage_force_history.clear();
    sim.step_count = 0;
    sim.active = true;
    sim.configured = true;
}

template<typename SP>
void SimBatch<SP>::sync_states_from_engine() {
    const auto& source = engine_->state();
    const auto D = static_cast<std::size_t>(n_stages_ + N_fil_);
    for (int id = 0; id < num_sims_; ++id) {
        auto& state = sims_[static_cast<std::size_t>(id)].state;
        const auto b = static_cast<std::size_t>(id);
        state.currents = Eigen::Map<const Eigen::VectorXd>(source.currents.data() + b * D,
                                                            static_cast<Eigen::Index>(D));
        state.arm_position = armature_.position() - source.position[b];
        state.arm_velocity = -source.velocity[b];
    }
}

template<typename SP>
bool SimBatch<SP>::is_stage_within_range(int stage, const SimInstance& sim) const {
    const auto& coil = coils_[static_cast<std::size_t>(stage)];
    // Keep the same centre-to-stage cutoff as GpuMultiStageSim, regardless of
    // whether the engine ultimately resolves to a CUDA or CPU backend.
    return std::abs(sim.state.arm_position - coil.position()) <= 10.0 * coil.length();
}

template<typename SP>
void SimBatch<SP>::check_triggers(SimInstance& sim) {
    const double time = sim.step_count * dt_;
    for (int stage = 1; stage < n_stages_; ++stage) {
        if (sim.triggered[static_cast<std::size_t>(stage)] ||
            !sim.triggered[static_cast<std::size_t>(stage - 1)]) continue;
        const auto& config = sim.trigger_configs[static_cast<std::size_t>(stage - 1)];
        const bool fire = config.mode == TriggerMode::Position
            ? sim.state.arm_position >= config.value
            : time >= sim.trigger_times[static_cast<std::size_t>(stage - 1)] + config.value;
        if (fire) {
            sim.triggered[static_cast<std::size_t>(stage)] = true;
            sim.trigger_times[static_cast<std::size_t>(stage)] = time;
            sim.trigger_positions[static_cast<std::size_t>(stage)] = sim.state.arm_position;
        }
    }
}

template<typename SP>
void SimBatch<SP>::extinguish_quiet_stages(SimInstance& sim) {
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!sim.triggered[static_cast<std::size_t>(stage)] ||
            sim.stage_completed[static_cast<std::size_t>(stage)]) continue;
        if (sim.excitations[static_cast<std::size_t>(stage)]->voltage() == 0.0 &&
            std::abs(sim.state.currents(stage)) < 1e-6) {
            sim.excitation_finished[static_cast<std::size_t>(stage)] = true;
            sim.stage_completed[static_cast<std::size_t>(stage)] = true;
            sim.state.currents(stage) = 0.0;
            engine_->complete_stage(
                static_cast<std::size_t>(&sim - sims_.data()),
                static_cast<std::size_t>(stage));
        }
    }
}

template<typename SP>
bool SimBatch<SP>::terminally_ineligible(const SimInstance& sim, int stage) const {
    for (int current = stage; current > 0; --current) {
        if (sim.trigger_configs[static_cast<std::size_t>(current - 1)].value == INFINITY) return true;
        if (sim.triggered[static_cast<std::size_t>(current - 1)]) return false;
    }
    return false;
}

template<typename SP>
bool SimBatch<SP>::check_termination(SimInstance& sim, const TerminationPolicy& policy) const {
    if (policy.enable_bound_check && sim.state.arm_position >= policy.barrel_end_position) return true;
    if (sim.step_count >= policy.max_steps) return true;
    bool all_finished = true;
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!sim.stage_completed[static_cast<std::size_t>(stage)] &&
            !terminally_ineligible(sim, stage)) {
            all_finished = false;
            break;
        }
    }
    if (all_finished) return true;
    if (policy.enable_velocity_check && sim.step_count >= policy.velocity_decay_steps) {
        const auto& history = sim.result.history;
        bool decaying = history.size() >= static_cast<std::size_t>(policy.velocity_decay_steps + 1);
        for (int i = 0; decaying && i < policy.velocity_decay_steps; ++i) {
            const auto& newer = history[history.size() - 1 - static_cast<std::size_t>(i)];
            const auto& older = history[history.size() - 2 - static_cast<std::size_t>(i)];
            if (newer.state.arm_velocity >= older.state.arm_velocity) decaying = false;
        }
        if (decaying && !history.empty()) {
            const double force = std::abs(history.back().state.force);
            if (force / armature_.mass() < policy.accel_threshold) return true;
        }
    }
    return false;
}

template<typename SP>
void SimBatch<SP>::configure_engine_boundary() {
    const auto B = static_cast<std::size_t>(num_sims_);
    const auto S = static_cast<std::size_t>(n_stages_);
    std::vector<std::uint8_t> active(B, 0), trigger(B * S, 0), stage(B * S, 0), mutual(B * S, 0);
    std::vector<double> voltages(B * S, 0.0);
    for (std::size_t b = 0; b < B; ++b) {
        auto& sim = sims_[b];
        active[b] = sim.active ? 1 : 0;
        for (std::size_t s = 0; s < S; ++s) {
            const auto index = b * S + s;
            trigger[index] = sim.triggered[s] ? 1 : 0;
            if (!sim.active || !sim.triggered[s] || sim.stage_completed[s]) continue;
            stage[index] = 1;
            mutual[index] = is_stage_within_range(static_cast<int>(s), sim) ? 1 : 0;
            voltages[index] = sim.excitation_finished[s]
                ? 0.0 : sim.excitations[s]->voltage();
        }
    }
    engine_->set_step_boundary_state(std::move(active), std::move(trigger), std::move(stage),
                                     std::move(mutual), std::move(voltages));

    std::vector<std::uint8_t> trigger_modes(B * S, 0);
    std::vector<double> trigger_values(B * S, INFINITY);
    std::vector<std::uint8_t> excitation_finished(B * S, 0);
    std::vector<std::uint8_t> stage_completed(B * S, 0);
    std::vector<double> trigger_times(B * S, 0.0);
    std::vector<double> trigger_positions(B * S, armature_.position());
    std::vector<double> position_offsets(B, armature_.position());
    for (std::size_t b = 0; b < B; ++b) {
        const auto& sim = sims_[b];
        for (std::size_t s = 0; s < S; ++s) {
            const auto index = b * S + s;
            excitation_finished[index] = sim.excitation_finished[s] ? 1 : 0;
            stage_completed[index] = sim.stage_completed[s] ? 1 : 0;
            trigger_times[index] = sim.trigger_times[s];
            trigger_positions[index] = sim.trigger_positions[s];
            if (s == 0) continue;
            const auto& config = sim.trigger_configs[s - 1];
            trigger_modes[index] = config.mode == TriggerMode::Position ? 1 : 2;
            trigger_values[index] = config.value;
        }
    }
    engine_->set_control_boundary_state(
        std::move(trigger_modes), std::move(trigger_values),
        std::move(excitation_finished), std::move(stage_completed),
        std::move(trigger_times), std::move(trigger_positions),
        std::move(position_offsets));
}

template<typename SP>
std::vector<double> SimBatch<SP>::compute_stage_forces(std::size_t sim_id,
    const std::vector<double>& currents, const std::vector<double>& gradients) const {
    const auto S = static_cast<std::size_t>(n_stages_);
    const auto F = static_cast<std::size_t>(N_fil_);
    const auto D = S + F;
    std::vector<double> forces(S, 0.0);
    const auto& sim = sims_[sim_id];
    for (std::size_t s = 0; s < S; ++s) {
        if (!sim.triggered[s] || sim.stage_completed[s]) continue;
        for (std::size_t f = 0; f < F; ++f)
            forces[s] -= currents[sim_id * D + s] * currents[sim_id * D + S + f] *
                         gradients[(sim_id * S + s) * F + f];
    }
    return forces;
}

template<typename SP>
void SimBatch<SP>::record_step(std::size_t sim_id, const std::vector<double>& stage_forces) {
    const auto S = static_cast<std::size_t>(n_stages_);
    const auto F = static_cast<std::size_t>(N_fil_);
    auto& sim = sims_[sim_id];
    MultiStageStep entry;
    entry.state.time = sim.step_count * dt_;
    entry.state.arm_position = sim.state.arm_position;
    entry.state.arm_velocity = sim.state.arm_velocity;
    entry.state.force = 0.0;
    for (const double force : stage_forces) entry.state.force += force;
    entry.state.filament_currents.resize(F);
    for (std::size_t f = 0; f < F; ++f) entry.state.filament_currents[f] = sim.state.currents(static_cast<Eigen::Index>(S + f));
    entry.cap_voltages.resize(S, 0.0);
    entry.coil_currents.resize(S, 0.0);
    for (std::size_t s = 0; s < S; ++s) {
        entry.cap_voltages[s] = sim.triggered[s] ? sim.excitations[s]->voltage() : 0.0;
        entry.coil_currents[s] = sim.triggered[s] ? sim.state.currents(static_cast<Eigen::Index>(s)) : 0.0;
    }
    entry.stage_forces = stage_forces;
    sim.result.history.push_back(std::move(entry));
    sim.stage_force_history.push_back(stage_forces);
    ++sim.step_count;
}

template<typename SP>
void SimBatch<SP>::prepare_summary(SimInstance& sim) {
    auto& summary = sim.result.summary;
    summary = MultiStageSummary{};
    summary.step_count = sim.step_count;
    if (sim.result.history.empty()) return;
    summary.total_time = sim.result.history.back().state.time;
    summary.muzzle_velocity = sim.result.history.back().state.arm_velocity;
    for (std::size_t s = 0; s < static_cast<std::size_t>(n_stages_); ++s) {
        if (!sim.triggered[s]) continue;
        PerStageSummary per;
        per.stage_index = static_cast<int>(s);
        per.trigger_time = sim.trigger_times[s];
        per.trigger_position = sim.trigger_positions[s];
        bool active = false;
        for (std::size_t h = 0; h < sim.result.history.size(); ++h) {
            const auto& entry = sim.result.history[h];
            per.peak_current = std::max(per.peak_current, entry.coil_currents[s]);
            if (h < sim.stage_force_history.size()) per.max_force = std::max(per.max_force, std::abs(sim.stage_force_history[h][s]));
            if (entry.coil_currents[s] > 1e-6) active = true;
            if (active) ++per.step_count_active;
        }
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(sim.excitations[s].get())) {
            per.energy_depleted = sim.initial_stage_energies[s] -
                0.5 * capacitor->capacitance() * capacitor->capacitor_voltage() * capacitor->capacitor_voltage();
        }
        summary.per_stage.push_back(per);
    }
    double input_energy = 0.0;
    for (const auto& excitation : sim.excitations)
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(excitation.get()))
            input_energy += 0.5 * capacitor->capacitance() * capacitor->initial_voltage() * capacitor->initial_voltage();
    for (const auto& entry : sim.result.history) {
        summary.max_force = std::max(summary.max_force, std::abs(entry.state.force));
        for (const double current : entry.coil_currents) summary.peak_coil_current = std::max(summary.peak_coil_current, current);
    }
    const double kinetic = 0.5 * armature_.mass() * summary.muzzle_velocity * summary.muzzle_velocity;
    summary.efficiency = input_energy > 0.0 ? kinetic / input_energy : 0.0;
}

template<typename SP>
void SimBatch<SP>::run() { run(TerminationPolicy::defaults()); }

template<typename SP>
void SimBatch<SP>::run(const TerminationPolicy& policy) {
    if constexpr (std::is_same_v<SP, RK4Stepper>) {
        throw std::logic_error("RK4Stepper is not supported by SimBatch");
    } else {
        for (const auto& sim : sims_)
            if (!sim.configured) throw std::logic_error("SimBatch: set_excitations must configure every simulation before run");

        bool any_active = true;
        while (any_active) {
            any_active = false;
            for (auto& sim : sims_) {
                if (!sim.active) continue;
                if (check_termination(sim, policy)) sim.active = false;
                else any_active = true;
            }
            if (!any_active) break;
            for (auto& sim : sims_) {
                if (!sim.active) continue;
                if (engine_->report().backend == BackendMode::Fallback) check_triggers(sim);
                extinguish_quiet_stages(sim);
            }
            configure_engine_boundary();
            std::vector<std::uint8_t> stepped(sims_.size(), 0);
            for (std::size_t b = 0; b < sims_.size(); ++b)
                stepped[b] = sims_[b].active ? 1 : 0;
            const auto pre_step_currents = engine_->state().currents;
            engine_->step();
            const auto boundary_gradients = engine_->state().dm1;
            sync_states_from_engine();
            if (engine_->report().backend != BackendMode::Fallback) {
                const auto& engine_state = engine_->state();
                for (std::size_t b = 0; b < sims_.size(); ++b) {
                    auto& sim = sims_[b];
                    sim.active = engine_state.active_mask[b] != 0;
                    for (std::size_t s = 0; s < static_cast<std::size_t>(n_stages_); ++s) {
                        const auto index = b * static_cast<std::size_t>(n_stages_) + s;
                        sim.triggered[s] = engine_state.trigger_mask[index] != 0;
                        sim.stage_completed[s] = engine_state.stage_completed[index] != 0;
                        sim.trigger_times[s] = engine_state.trigger_times[index];
                        sim.trigger_positions[s] = engine_state.trigger_positions[index];
                    }
                }
            }
            for (std::size_t b = 0; b < sims_.size(); ++b) {
                auto& sim = sims_[b];
                if (stepped[b] == 0) continue;
                for (std::size_t s = 0; s < static_cast<std::size_t>(n_stages_); ++s) {
                    if (sim.triggered[s] && !sim.excitation_finished[s]) {
                        const auto dimension = static_cast<std::size_t>(n_stages_ + N_fil_);
                        sim.excitations[s]->advance(
                            dt_, pre_step_currents[b * dimension + s]);
                        if (sim.excitations[s]->finished())
                            sim.excitation_finished[s] = true;
                    }
                    if (sim.triggered[s] && sim.excitation_finished[s] &&
                        std::abs(sim.state.currents(static_cast<Eigen::Index>(s))) < 1e-6) {
                        sim.stage_completed[s] = true;
                        sim.state.currents(static_cast<Eigen::Index>(s)) = 0.0;
                        engine_->complete_stage(b, s);
                    }
                }
                // The engine used pre-step currents for physical acceleration.
                // Recorded history uses post-step currents with the gradient
                // cached at the pre-position step boundary, masked after
                // excitation advancement.
                const auto forces = compute_stage_forces(b, engine_->state().currents,
                                                         boundary_gradients);
                record_step(b, forces);
            }
        }
        for (auto& sim : sims_) prepare_summary(sim);
    }
}

template<typename SP>
const MultiStageResult& SimBatch<SP>::result(int sim_id) const {
    if (sim_id < 0 || sim_id >= num_sims_) throw std::out_of_range("sim_id out of range");
    return sims_[static_cast<std::size_t>(sim_id)].result;
}

template class SimBatch<EulerStepper>;
template class SimBatch<RK4Stepper>;

} // namespace coilgun::simulation::cuda
