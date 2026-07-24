#include "coilgun/simulation/multi_stage_sim.hpp"

#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace coilgun::simulation {
namespace {

constexpr double kQuietCurrent = 1e-6;
constexpr double kEventTolerance = 1e-12;
constexpr int kMaxEventSegments = 16;
constexpr int kMaxEventBisections = 64;

bool snapshot_finished(const ExcitationSnapshot& snapshot) {
    if (const auto* value = dynamic_cast<const CapacitorSnapshot*>(&snapshot))
        return value->finished;
    if (const auto* value = dynamic_cast<const CrowbarSnapshot*>(&snapshot))
        return value->finished;
    if (const auto* value = dynamic_cast<const WaveformSnapshot*>(&snapshot))
        return value->finished;
    throw std::invalid_argument("unsupported excitation snapshot");
}

double event_value(const Excitation& source, const ExcitationSnapshot& snapshot) {
    if (const auto* value = dynamic_cast<const CapacitorSnapshot*>(&snapshot))
        return value->finished ? -1.0 : value->capacitor_voltage;
    if (const auto* value = dynamic_cast<const CrowbarSnapshot*>(&snapshot))
        return value->diode_on ? 1.0 : value->capacitor_voltage;
    if (const auto* value = dynamic_cast<const WaveformSnapshot*>(&snapshot)) {
        const auto* waveform = dynamic_cast<const WaveformExcitation*>(&source);
        return waveform ? waveform->end_time() - value->time : 1.0;
    }
    return 1.0;
}

} // namespace

template<typename SP>
MultiStageSim<SP>::MultiStageSim(
    std::vector<components::DrivingCoil> coils,
    components::Armature armature,
    std::vector<std::unique_ptr<Excitation>> excitations,
    std::vector<TriggerConfig> trigger_configs,
    double dt,
    bool enable_thermal,
    OptimizationLevel opt_level)
    : n_stages_(static_cast<int>(coils.size())), coils_(std::move(coils)),
      armature_(std::move(armature)), excitations_(std::move(excitations)),
      trigger_configs_(std::move(trigger_configs)), dt_(dt),
      enable_thermal_(enable_thermal), opt_level_(opt_level),
      N_fil_(armature_.total_filaments()) {
    if (n_stages_ <= 0) throw std::invalid_argument("MultiStageSim requires at least one coil");
    if (n_stages_ > kMaxStages) throw std::invalid_argument("MultiStageSim stage count exceeds kMaxStages");
    if (excitations_.size() != coils_.size())
        throw std::invalid_argument("coils.size() must equal excitations.size()");
    if (trigger_configs_.size() + 1 != coils_.size())
        throw std::invalid_argument("trigger_configs.size() must equal coils.size() - 1");
    if (!std::isfinite(dt_) || dt_ <= 0.0)
        throw std::invalid_argument("dt must be finite and positive");
    for (const auto& excitation : excitations_)
        if (!excitation) throw std::invalid_argument("excitation pointer must not be null");
    for (const auto& trigger : trigger_configs_) validate_trigger_config(trigger);

    R_diag_.resize(n_stages_);
    L_diag_.resize(n_stages_);
    initial_stage_energies_.assign(n_stages_, 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        R_diag_(stage) = coils_[stage].resistance();
        L_diag_(stage) = coils_[stage].self_inductance();
        if (const auto* capacitor =
                dynamic_cast<const CapacitorExcitation*>(excitations_[stage].get())) {
            initial_stage_energies_[stage] = 0.5 * capacitor->capacitance() *
                capacitor->initial_voltage() * capacitor->initial_voltage();
        }
    }

    R_fil_ref_.resize(N_fil_);
    L_fil_.resize(N_fil_);
    mass_fil_.resize(N_fil_);
    filament_inner_radii_.resize(N_fil_);
    filament_outer_radii_.resize(N_fil_);
    filament_mean_radii_.resize(N_fil_);
    filament_relative_axial_positions_.resize(N_fil_);
    filament_lengths_.resize(N_fil_);
    const auto& metadata = armature_.filament_metadata();
    for (int filament = 0; filament < N_fil_; ++filament) {
        R_fil_ref_(filament) = armature_.resistances()[filament];
        L_fil_(filament) = armature_.inductances()[filament];
        mass_fil_(filament) = armature_.masses()[filament];
        filament_inner_radii_[filament] = metadata[filament].inner_radius;
        filament_outer_radii_[filament] = metadata[filament].outer_radius;
        filament_mean_radii_[filament] = metadata[filament].mean_radius;
        filament_relative_axial_positions_[filament] =
            metadata[filament].axial_center - armature_.position();
        filament_lengths_[filament] = metadata[filament].length;
    }
    build_filament_M_matrix();
    precompute_M_cc();
    workspace_.resize(static_cast<std::size_t>(n_stages_),
                      static_cast<std::size_t>(N_fil_));

    auto& physical = integration_state_.physical;
    physical.currents = Eigen::VectorXd::Zero(n_stages_ + N_fil_);
    physical.arm_position = armature_.position();
    physical.arm_velocity = armature_.velocity();
    if (enable_thermal_)
        physical.filament_temperatures =
            Eigen::VectorXd::Constant(N_fil_, physics::T_REFERENCE);

    integration_state_.stages.resize(n_stages_);
    integration_state_.excitations.reserve(n_stages_);
    for (int stage = 0; stage < n_stages_; ++stage) {
        integration_state_.excitations.push_back(excitations_[stage]->snapshot());
        auto& runtime = integration_state_.stages[stage];
        runtime.triggered = stage == 0;
        runtime.circuit_active = stage == 0;
        runtime.trigger_time = 0.0;
        runtime.trigger_position = physical.arm_position;
    }
    initial_integration_state_ = clone_integration_state(integration_state_);
}

template<typename SP>
const StageRuntimeState& MultiStageSim<SP>::stage_state(std::size_t stage) const {
    return coilgun::simulation::stage_state(integration_state_, stage);
}

template<typename SP>
bool MultiStageSim<SP>::circuit_active(std::size_t stage) const noexcept {
    return stage < integration_state_.stages.size() &&
           integration_state_.stages[stage].circuit_active;
}

template<typename SP>
void MultiStageSim<SP>::build_filament_M_matrix() {
    M_mat_ = Eigen::MatrixXd::Zero(N_fil_, N_fil_);
    for (int first = 0; first < N_fil_; ++first) {
        const double first_radius = filament_mean_radii_[first];
        for (int second = first + 1; second < N_fil_; ++second) {
            const double second_radius = filament_mean_radii_[second];
            const double mutual = physics::mutual_inductance_filament(
                first_radius, second_radius,
                std::abs(filament_relative_axial_positions_[first] -
                         filament_relative_axial_positions_[second]), false);
            M_mat_(first, second) = mutual;
            M_mat_(second, first) = mutual;
        }
    }
}

template<typename SP>
void MultiStageSim<SP>::precompute_M_cc() {
    M_cc_ = Eigen::MatrixXd::Zero(n_stages_, n_stages_);
    for (int first = 0; first < n_stages_; ++first) {
        for (int second = first + 1; second < n_stages_; ++second) {
            const double mutual = physics::mutual_inductance_coil(
                coils_[first].inner_radius(), coils_[first].outer_radius(),
                coils_[first].length(), coils_[first].turns(),
                coils_[second].inner_radius(), coils_[second].outer_radius(),
                coils_[second].length(), coils_[second].turns(),
                std::abs(coils_[first].position() - coils_[second].position()),
                9, false);
            M_cc_(first, second) = mutual;
            M_cc_(second, first) = mutual;
        }
    }
}

template<typename SP>
bool MultiStageSim<SP>::is_stage_within_range(int stage) const {
    if (opt_level_ != OptimizationLevel::Full) return true;
    return std::abs(integration_state_.physical.arm_position - coils_[stage].position()) <=
           10.0 * coils_[stage].length();
}

template<typename SP>
Eigen::VectorXd MultiStageSim<SP>::derive_resistance(
    const MultiStageState& state) const {
    Eigen::VectorXd resistance = R_fil_ref_;
    if (!enable_thermal_ || state.filament_temperatures.size() == 0) return resistance;
    const double beta = physics::material_beta(armature_.material());
    for (int filament = 0; filament < N_fil_; ++filament) {
        resistance(filament) = R_fil_ref_(filament) *
            (1.0 + beta * (state.filament_temperatures(filament) - physics::T_REFERENCE));
        if (!std::isfinite(resistance(filament)) || resistance(filament) <= 0.0)
            throw std::runtime_error("derived filament resistance is not positive");
    }
    return resistance;
}

template<typename SP>
double MultiStageSim<SP>::compute_force(
    const MultiStageState& state,
    const Eigen::MatrixXd& gradient,
    const std::vector<StageRuntimeState>& stages) const {
    double force = 0.0;
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!stages[stage].circuit_active || !is_stage_within_range(stage)) continue;
        for (int filament = 0; filament < N_fil_; ++filament)
            force += state.currents(stage) * state.currents(n_stages_ + filament) *
                     gradient(stage, filament);
    }
    return force;
}

template<typename SP>
DerivativeResult MultiStageSim<SP>::evaluate_derivatives(
    const IntegrationState& state,
    DerivativeWorkspace& workspace) const {
#if COILGUN_ENABLE_CPU_PHASE_TIMING
    CpuPhaseTimingScope orchestration_timing(CpuPhase::Orchestration);
#endif
    workspace.resize(static_cast<std::size_t>(n_stages_),
                     static_cast<std::size_t>(N_fil_));
#if COILGUN_ENABLE_CPU_PHASE_TIMING
    if (enable_thermal_) {
        CpuPhaseTimingScope thermal_timing(CpuPhase::Thermal);
        workspace.resistance = derive_resistance(state.physical);
    } else {
        workspace.resistance = derive_resistance(state.physical);
    }
#else
    workspace.resistance = derive_resistance(state.physical);
#endif
    auto& matrix = workspace.system_matrix;
    auto& rhs = workspace.rhs;
    matrix.setZero();
    rhs.setZero();

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mutual(
        workspace.mutual.data(), n_stages_, N_fil_);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> gradient(
        workspace.mutual_gradient.data(), n_stages_, N_fil_);
    const int dimension = n_stages_ + N_fil_;
    mutual.setZero();
    gradient.setZero();

{
#if COILGUN_ENABLE_CPU_PHASE_TIMING
    CpuPhaseTimingScope mutual_timing(CpuPhase::Mutual);
#endif
    for (int stage = 0; stage < n_stages_; ++stage) {
        const auto& runtime = state.stages[stage];
        if (!runtime.circuit_active || !is_stage_within_range(stage)) continue;
#pragma omp parallel for if (N_fil_ >= 8)
        for (int filament = 0; filament < N_fil_; ++filament) {
            const double separation =
                state.physical.arm_position +
                filament_relative_axial_positions_[filament] -
                coils_[stage].position();
            int nodes = 9;
            if (opt_level_ == OptimizationLevel::Full &&
                std::abs(state.physical.arm_position - coils_[stage].position()) >
                    coils_[stage].length())
                nodes = 4;
            const auto pair = physics::mutual_detail::mutual_inductance_coil_pair(
                coils_[stage].inner_radius(), coils_[stage].outer_radius(),
                coils_[stage].length(), coils_[stage].turns(),
                filament_inner_radii_[filament], filament_outer_radii_[filament],
                filament_lengths_[filament], 1,
                separation, nodes, false);
            mutual(stage, filament) = pair.mutual;
            gradient(stage, filament) = pair.gradient;
        }
    }
}

{
#if COILGUN_ENABLE_CPU_PHASE_TIMING
    CpuPhaseTimingScope assembly_timing(CpuPhase::Assembly);
#endif
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!state.stages[stage].triggered || state.stages[stage].stage_completed) {
            matrix(stage, stage) = 1.0;
            continue;
        }
        matrix(stage, stage) = L_diag_(stage);
        for (int other = 0; other < n_stages_; ++other) {
            if (other != stage && state.stages[other].circuit_active)
                matrix(stage, other) = M_cc_(stage, other);
        }
        if (state.stages[stage].circuit_active) {
            for (int filament = 0; filament < N_fil_; ++filament) {
                matrix(stage, n_stages_ + filament) = mutual(stage, filament);
                matrix(n_stages_ + filament, stage) = mutual(stage, filament);
            }
        }
    }
    for (int filament = 0; filament < N_fil_; ++filament) {
        const int row = n_stages_ + filament;
        matrix(row, row) = L_fil_(filament);
        for (int other = 0; other < N_fil_; ++other)
            if (filament != other) matrix(row, n_stages_ + other) = M_mat_(filament, other);
    }

    const double velocity = state.physical.arm_velocity;
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!state.stages[stage].circuit_active) continue;
        double motional = 0.0;
        for (int filament = 0; filament < N_fil_; ++filament)
            motional += gradient(stage, filament) *
                        state.physical.currents(n_stages_ + filament);
        const double source_voltage = state.stages[stage].excitation_finished
            ? 0.0 : excitations_[stage]->voltage(*state.excitations[stage]);
        rhs(stage) = source_voltage - R_diag_(stage) * state.physical.currents(stage) -
                     velocity * motional;
    }
    for (int filament = 0; filament < N_fil_; ++filament) {
        double back_emf = 0.0;
        for (int stage = 0; stage < n_stages_; ++stage) {
            if (state.stages[stage].circuit_active)
                back_emf += gradient(stage, filament) * state.physical.currents(stage);
        }
        rhs(n_stages_ + filament) =
            -workspace.resistance(filament) *
                state.physical.currents(n_stages_ + filament) -
            velocity * back_emf;
    }
}

    DerivativeResult result;
    result.physical_derivative.currents.resize(dimension);
{
#if COILGUN_ENABLE_CPU_PHASE_TIMING
    CpuPhaseTimingScope solve_timing(CpuPhase::Solve);
#endif
    Eigen::LDLT<Eigen::MatrixXd> solver(matrix);
    if (solver.info() == Eigen::Success)
        result.physical_derivative.currents = solver.solve(rhs);
    else
        result.physical_derivative.currents = matrix.colPivHouseholderQr().solve(rhs);
    if (!result.physical_derivative.currents.allFinite())
        throw std::runtime_error("multi-stage circuit solve produced non-finite derivatives");
}

    result.force = compute_force(state.physical, gradient, state.stages);
    result.physical_derivative.arm_position = state.physical.arm_velocity;
    result.physical_derivative.arm_velocity = result.force / armature_.mass();
    if (enable_thermal_) {
#if COILGUN_ENABLE_CPU_PHASE_TIMING
        CpuPhaseTimingScope thermal_timing(CpuPhase::Thermal);
#endif
        result.physical_derivative.filament_temperatures.resize(N_fil_);
        for (int filament = 0; filament < N_fil_; ++filament) {
            const double current = state.physical.currents(n_stages_ + filament);
            const double temperature = state.physical.filament_temperatures(filament);
            result.physical_derivative.filament_temperatures(filament) =
                current * current * workspace.resistance(filament) /
                (mass_fil_(filament) *
                 physics::material_cp(armature_.material(), temperature));
        }
    }

    result.excitation_derivatives.resize(n_stages_);
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (state.stages[stage].triggered && !state.stages[stage].excitation_finished) {
            result.excitation_derivatives[stage] =
                excitations_[stage]->continuous_derivative(
                    *state.excitations[stage], state.physical.currents(stage));
        }
    }
    return result;
}

template<typename SP>
IntegrationState MultiStageSim<SP>::make_trial(
    const IntegrationState& initial,
    const DerivativeResult& derivative,
    double scale) const {
    auto trial = clone_integration_state(initial);
    trial.physical += scale * derivative.physical_derivative;
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (trial.stages[stage].triggered && !trial.stages[stage].excitation_finished) {
            excitations_[stage]->advance_snapshot_derivative(
                *trial.excitations[stage], scale,
                derivative.excitation_derivatives[stage]);
        }
    }
    return trial;
}

template<typename SP>
void MultiStageSim<SP>::check_triggers(
    IntegrationState& state,
    const IntegrationState& pre,
    double absolute_time) const {
    for (int stage = 1; stage < n_stages_; ++stage) {
        auto& runtime = state.stages[stage];
        if (runtime.triggered || !state.stages[stage - 1].triggered) continue;
        const auto& trigger = trigger_configs_[stage - 1];
        if (trigger.value == std::numeric_limits<double>::infinity()) continue;

        bool fire = false;
        double event_time = absolute_time;
        double event_position = pre.physical.arm_position;
        if (trigger.mode == TriggerMode::Position) {
            const double before = pre.physical.arm_position - trigger.value;
            const double after = state.physical.arm_position - trigger.value;
            fire = before >= 0.0 || (before < 0.0 && after >= 0.0);
            if (fire && before < 0.0 && after > before) {
                const double fraction = -before / (after - before);
                event_time = absolute_time + fraction * dt_;
                event_position = pre.physical.arm_position;
            }
        } else {
            const double target = state.stages[stage - 1].trigger_time + trigger.value;
            fire = absolute_time + dt_ >= target;
            if (fire) event_time = std::max(absolute_time, target);
        }
        if (fire) {
            runtime.triggered = true;
            runtime.circuit_active = true;
            runtime.trigger_time = event_time;
            runtime.trigger_position = event_position;
        }
    }
}

template<typename SP>
void MultiStageSim<SP>::complete_quiet_stages(IntegrationState& state) const {
    for (int stage = 0; stage < n_stages_; ++stage) {
        auto& runtime = state.stages[stage];
        if (runtime.circuit_active && runtime.excitation_finished &&
            std::abs(state.physical.currents(stage)) < kQuietCurrent) {
            mark_stage_completed(state, static_cast<std::size_t>(stage));
        }
    }
}

template<typename SP>
void MultiStageSim<SP>::apply_boundary_events(
    IntegrationState& state,
    const IntegrationState& pre,
    double,
    double absolute_time) const {
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!state.stages[stage].triggered || state.stages[stage].excitation_finished) continue;
        auto& snapshot = *state.excitations[stage];
        if (auto* capacitor = dynamic_cast<CapacitorSnapshot*>(&snapshot)) {
            if (capacitor->capacitor_voltage <= 0.0) {
                excitations_[stage]->apply_event(snapshot, ExcitationEvent::CapacitorZero);
                excitations_[stage]->apply_event(snapshot, ExcitationEvent::Finished);
            }
        } else if (auto* crowbar = dynamic_cast<CrowbarSnapshot*>(&snapshot)) {
            if (!crowbar->diode_on && crowbar->capacitor_voltage <= 0.0) {
                excitations_[stage]->apply_event(snapshot, ExcitationEvent::CapacitorZero);
                excitations_[stage]->apply_event(snapshot, ExcitationEvent::CrowbarOn);
            }
            state.stages[stage].crowbar_on = crowbar->diode_on;
            if (crowbar->diode_on && std::abs(state.physical.currents(stage)) < kQuietCurrent)
                excitations_[stage]->apply_event(snapshot, ExcitationEvent::Finished);
        } else if (auto* waveform = dynamic_cast<WaveformSnapshot*>(&snapshot)) {
            const auto* source = dynamic_cast<const WaveformExcitation*>(excitations_[stage].get());
            if (source && waveform->time >= source->end_time())
                excitations_[stage]->apply_event(snapshot, ExcitationEvent::WaveformEnd);
        }
        state.stages[stage].excitation_finished = snapshot_finished(snapshot);
    }
    check_triggers(state, pre, absolute_time);
    complete_quiet_stages(state);
}

template<typename SP>
IntegrationState MultiStageSim<SP>::advance_euler(
    const IntegrationState& original_pre, double dt) {
    auto pre = clone_integration_state(original_pre);
    check_triggers(pre, original_pre, step_count_ * dt_);
    const auto derivative = evaluate_derivatives(pre, workspace_);
    pre_step_mutual_gradient_ = workspace_.mutual_gradient;
    auto post = make_trial(pre, derivative, dt);
    apply_boundary_events(post, pre, dt, step_count_ * dt_);
    return post;
}

template<typename SP>
IntegrationState MultiStageSim<SP>::advance_rk4_segment(
    const IntegrationState& pre, double dt) {
    const auto k1 = evaluate_derivatives(pre, workspace_);
    pre_step_mutual_gradient_ = workspace_.mutual_gradient;
    const auto s2 = make_trial(pre, k1, 0.5 * dt);
    const auto k2 = evaluate_derivatives(s2, workspace_);
    const auto s3 = make_trial(pre, k2, 0.5 * dt);
    const auto k3 = evaluate_derivatives(s3, workspace_);
    const auto s4 = make_trial(pre, k3, dt);
    const auto k4 = evaluate_derivatives(s4, workspace_);

    auto post = clone_integration_state(pre);
    auto physical = k1.physical_derivative;
    physical += 2.0 * k2.physical_derivative;
    physical += 2.0 * k3.physical_derivative;
    physical += k4.physical_derivative;
    post.physical += (dt / 6.0) * physical;
    for (int stage = 0; stage < n_stages_; ++stage) {
        ExcitationDerivative derivative;
        derivative.capacitor_voltage_rate =
            (k1.excitation_derivatives[stage].capacitor_voltage_rate +
             2.0 * k2.excitation_derivatives[stage].capacitor_voltage_rate +
             2.0 * k3.excitation_derivatives[stage].capacitor_voltage_rate +
             k4.excitation_derivatives[stage].capacitor_voltage_rate) / 6.0;
        derivative.waveform_time_rate =
            (k1.excitation_derivatives[stage].waveform_time_rate +
             2.0 * k2.excitation_derivatives[stage].waveform_time_rate +
             2.0 * k3.excitation_derivatives[stage].waveform_time_rate +
             k4.excitation_derivatives[stage].waveform_time_rate) / 6.0;
        if (post.stages[stage].triggered && !post.stages[stage].excitation_finished)
            excitations_[stage]->advance_snapshot_derivative(
                *post.excitations[stage], dt, derivative);
    }
    return post;
}

template<typename SP>
IntegrationState MultiStageSim<SP>::advance_rk4_event_aware(
    const IntegrationState& original_pre, double dt) {
    pre_step_mutual_gradient_.resize(0);
    auto current = clone_integration_state(original_pre);
    check_triggers(current, original_pre, step_count_ * dt_);
    double remaining = dt;
    double elapsed = 0.0;
    for (int segment = 0; segment < kMaxEventSegments && remaining > 0.0; ++segment) {
        auto candidate = advance_rk4_segment(current, remaining);
        int crossing_stage = -1;
        for (int stage = 0; stage < n_stages_; ++stage) {
            if (!current.stages[stage].triggered || current.stages[stage].excitation_finished) continue;
            const double before = event_value(*excitations_[stage], *current.excitations[stage]);
            const double after = event_value(*excitations_[stage], *candidate.excitations[stage]);
            if (before > 0.0 && after <= 0.0) {
                crossing_stage = stage;
                break;
            }
        }
        if (crossing_stage < 0) {
            apply_boundary_events(candidate, current, remaining,
                                  step_count_ * dt_ + elapsed);
            return candidate;
        }

        double lower = 0.0;
        double upper = remaining;
        for (int iteration = 0; iteration < kMaxEventBisections; ++iteration) {
            const double midpoint = 0.5 * (lower + upper);
            const auto trial = advance_rk4_segment(current, midpoint);
            if (event_value(*excitations_[crossing_stage],
                            *trial.excitations[crossing_stage]) > 0.0)
                lower = midpoint;
            else
                upper = midpoint;
            if (upper - lower <= kEventTolerance) break;
        }
        if (upper <= 0.0) throw std::runtime_error("zero-duration repeated RK4 event");
        auto event_state = advance_rk4_segment(current, upper);
        apply_boundary_events(event_state, current, upper,
                              step_count_ * dt_ + elapsed);
        current = std::move(event_state);
        elapsed += upper;
        remaining -= upper;
    }
    if (remaining > kEventTolerance)
        throw std::runtime_error("multi-stage RK4 event segment limit exceeded");
    return current;
}

template<typename SP>
void MultiStageSim<SP>::record_step(double post_time) {
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> gradient(
        pre_step_mutual_gradient_.data(), n_stages_, N_fil_);
    MultiStageStep entry;
    entry.state.time = post_time;
    entry.state.arm_position = integration_state_.physical.arm_position;
    entry.state.arm_velocity = integration_state_.physical.arm_velocity;
    entry.state.filament_currents.resize(N_fil_);
    for (int filament = 0; filament < N_fil_; ++filament)
        entry.state.filament_currents[filament] =
            integration_state_.physical.currents(n_stages_ + filament);
    if (enable_thermal_) {
        entry.state.filament_temperatures.resize(N_fil_);
        for (int filament = 0; filament < N_fil_; ++filament)
            entry.state.filament_temperatures[filament] =
                integration_state_.physical.filament_temperatures(filament);
    }

    entry.cap_voltages.resize(n_stages_);
    entry.coil_currents.resize(n_stages_);
    entry.stage_forces.assign(n_stages_, 0.0);
    for (int stage = 0; stage < n_stages_; ++stage) {
        const auto& runtime = integration_state_.stages[stage];
        if (runtime.triggered)
            entry.cap_voltages[stage] =
                excitations_[stage]->voltage(*integration_state_.excitations[stage]);
        entry.coil_currents[stage] = integration_state_.physical.currents(stage);
        if (!runtime.circuit_active) continue;
        for (int filament = 0; filament < N_fil_; ++filament) {
            entry.stage_forces[stage] += integration_state_.physical.currents(stage) *
                integration_state_.physical.currents(n_stages_ + filament) *
                gradient(stage, filament);
        }
        entry.state.force += entry.stage_forces[stage];
    }
    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool MultiStageSim<SP>::check_all_finished() const {
    for (int stage = 0; stage < n_stages_; ++stage) {
        if (integration_state_.stages[stage].stage_completed) continue;
        if (!integration_state_.stages[stage].triggered && stage > 0 &&
            trigger_configs_[stage - 1].value == std::numeric_limits<double>::infinity())
            continue;
        return false;
    }
    return true;
}

template<typename SP>
bool MultiStageSim<SP>::check_termination(const TerminationPolicy& policy) {
    const auto& state = integration_state_.physical;
    if (policy.enable_bound_check && state.arm_position >= policy.barrel_end_position) return true;
    if (step_count_ >= policy.max_steps || check_all_finished()) return true;
    for (const auto& stage : integration_state_.stages)
        if (stage.triggered && stage.circuit_active && !stage.excitation_finished)
            return false;
    if (policy.enable_velocity_check && step_count_ >= policy.velocity_decay_steps) {
        const double acceleration = result_.history.empty() ? 0.0 :
            result_.history.back().state.force / armature_.mass();
        bool decaying = true;
        const auto count = static_cast<int>(result_.history.size());
        for (int index = 0; index < policy.velocity_decay_steps; ++index) {
            if (count - 2 - index < 0 ||
                result_.history[count - 1 - index].state.arm_velocity >=
                    result_.history[count - 2 - index].state.arm_velocity) {
                decaying = false;
                break;
            }
        }
        if (decaying && std::abs(acceleration) < policy.accel_threshold)
            return true;
    }
    return false;
}

template<typename SP>
void MultiStageSim<SP>::prepare_summary() {
    auto& summary = result_.summary;
    summary = {};
    summary.step_count = step_count_;
    if (result_.history.empty()) return;
    summary.total_time = result_.history.back().state.time;
    summary.muzzle_velocity = result_.history.back().state.arm_velocity;

    for (int stage = 0; stage < n_stages_; ++stage) {
        if (!integration_state_.stages[stage].triggered) continue;
        PerStageSummary stage_summary;
        stage_summary.stage_index = stage;
        stage_summary.trigger_time = integration_state_.stages[stage].trigger_time;
        stage_summary.trigger_position = integration_state_.stages[stage].trigger_position;
        for (const auto& step : result_.history) {
            stage_summary.peak_current = std::max(
                stage_summary.peak_current, std::abs(step.coil_currents[stage]));
            stage_summary.max_force = std::max(
                stage_summary.max_force, std::abs(step.stage_forces[stage]));
            if (std::abs(step.coil_currents[stage]) >= kQuietCurrent)
                ++stage_summary.step_count_active;
        }
        if (const auto* capacitor =
                dynamic_cast<const CapacitorExcitation*>(excitations_[stage].get())) {
            const double voltage = result_.history.back().cap_voltages[stage];
            stage_summary.energy_depleted = initial_stage_energies_[stage] -
                0.5 * capacitor->capacitance() * voltage * voltage;
        }
        summary.per_stage.push_back(stage_summary);
    }
    for (const auto& step : result_.history) {
        summary.max_force = std::max(summary.max_force, std::abs(step.state.force));
        for (const double current : step.coil_currents)
            summary.peak_coil_current =
                std::max(summary.peak_coil_current, std::abs(current));
    }
    double input_energy = 0.0;
    for (const double energy : initial_stage_energies_) input_energy += energy;
    const double kinetic = 0.5 * armature_.mass() *
        summary.muzzle_velocity * summary.muzzle_velocity;
    summary.efficiency = input_energy > 0.0 ? kinetic / input_energy : 0.0;
}

template<typename SP>
const MultiStageStep& MultiStageSim<SP>::step() {
    const auto pre = clone_integration_state(integration_state_);
    if constexpr (std::is_same_v<SP, RK4Stepper>)
        integration_state_ = advance_rk4_event_aware(pre, dt_);
    else
        integration_state_ = advance_euler(pre, dt_);
    for (int stage = 0; stage < n_stages_; ++stage)
        excitations_[stage]->restore(*integration_state_.excitations[stage]);
    record_step((step_count_ + 1) * dt_);
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const MultiStageResult& MultiStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const MultiStageResult& MultiStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy)) step();
    prepare_summary();
    return result_;
}

template<typename SP>
void MultiStageSim<SP>::reset() {
    integration_state_ = clone_integration_state(initial_integration_state_);
    for (int stage = 0; stage < n_stages_; ++stage)
        excitations_[stage]->restore(*integration_state_.excitations[stage]);
    result_ = {};
    step_count_ = 0;
}

template class MultiStageSim<EulerStepper>;
template class MultiStageSim<RK4Stepper>;

} // namespace coilgun::simulation
