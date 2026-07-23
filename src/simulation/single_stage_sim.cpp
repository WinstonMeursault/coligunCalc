#include "coilgun/simulation/single_stage_sim.hpp"

#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace coilgun::simulation {
namespace {

constexpr double kQuietCurrent = 1e-6;
constexpr double kEventTimeTolerance = 1e-12;
constexpr int kMaxEventBisections = 64;
constexpr int kMaxEventSegments = 16;

double filament_axial_sep(const components::Armature& arm, int first, int second) {
    const int radial = arm.radial_filaments();
    const int first_axial = first / radial + 1;
    const int second_axial = second / radial + 1;
    return std::abs(arm.filament_axial_position(first_axial) -
                    arm.filament_axial_position(second_axial));
}

bool snapshot_finished(const ExcitationSnapshot& snapshot) {
    if (const auto* value = dynamic_cast<const CapacitorSnapshot*>(&snapshot))
        return value->finished;
    if (const auto* value = dynamic_cast<const CrowbarSnapshot*>(&snapshot))
        return value->finished;
    if (const auto* value = dynamic_cast<const WaveformSnapshot*>(&snapshot))
        return value->finished;
    throw std::invalid_argument("unsupported excitation snapshot");
}

} // namespace

template<typename SP>
SingleStageSim<SP>::SingleStageSim(
    components::DrivingCoil coil,
    components::Armature armature,
    std::unique_ptr<Excitation> excitation,
    double dt,
    bool enable_thermal)
    : coil_(std::move(coil)), armature_(std::move(armature)),
      excitation_(std::move(excitation)), dt_(dt), enable_thermal_(enable_thermal) {
    if (!excitation_) throw std::invalid_argument("excitation must not be null");
    if (!std::isfinite(dt_) || dt_ <= 0.0)
        throw std::invalid_argument("dt must be finite and positive");

    N_fil_ = armature_.total_filaments();
    R_d_ = coil_.resistance();
    L_d_ = coil_.self_inductance();

    R_fil_ref_.resize(N_fil_);
    L_fil_.resize(N_fil_);
    mass_fil_.resize(N_fil_);
    for (int filament = 0; filament < N_fil_; ++filament) {
        R_fil_ref_(filament) = armature_.resistances()[filament];
        L_fil_(filament) = armature_.inductances()[filament];
        mass_fil_(filament) = armature_.masses()[filament];
    }

    build_filament_M_matrix();
    workspace_.resize(1, static_cast<std::size_t>(N_fil_));

    auto& physical = integration_state_.physical;
    physical.currents = Eigen::VectorXd::Zero(N_fil_ + 1);
    physical.arm_position = armature_.position();
    physical.arm_velocity = armature_.velocity();
    if (enable_thermal_)
        physical.filament_temperatures =
            Eigen::VectorXd::Constant(N_fil_, physics::T_REFERENCE);

    integration_state_.excitations.push_back(excitation_->snapshot());
    StageRuntimeState stage;
    stage.triggered = true;
    stage.circuit_active = true;
    stage.trigger_time = 0.0;
    stage.trigger_position = physical.arm_position;
    integration_state_.stages.push_back(stage);
    initial_integration_state_ = clone_integration_state(integration_state_);
}

template<typename SP>
void SingleStageSim<SP>::build_filament_M_matrix() {
    M_mat_ = Eigen::MatrixXd::Zero(N_fil_, N_fil_);
    const int radial = armature_.radial_filaments();
    for (int first = 0; first < N_fil_; ++first) {
        const double first_radius = armature_.filament_mean_radius(first % radial + 1);
        for (int second = first + 1; second < N_fil_; ++second) {
            const double second_radius = armature_.filament_mean_radius(second % radial + 1);
            const double mutual = physics::mutual_inductance_filament(
                first_radius, second_radius,
                filament_axial_sep(armature_, first, second), false);
            M_mat_(first, second) = mutual;
            M_mat_(second, first) = mutual;
        }
    }
}

template<typename SP>
Eigen::VectorXd SingleStageSim<SP>::derive_resistance(const SimState& state) const {
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
double SingleStageSim<SP>::compute_force(
    const SimState& state, const Eigen::VectorXd& mutual_gradient) const {
    double force = 0.0;
    const double coil_current = state.currents(0);
    for (int filament = 0; filament < N_fil_; ++filament)
        force += coil_current * state.currents(filament + 1) * mutual_gradient(filament);
    return force;
}

template<typename SP>
DerivativeResult SingleStageSim<SP>::evaluate_derivatives(
    const SimState& state,
    const ExcitationSnapshot& excitation,
    DerivativeWorkspace& workspace,
    bool circuit_active) const {
    workspace.resize(1, static_cast<std::size_t>(N_fil_));
    auto& matrix = workspace.system_matrix;
    auto& rhs = workspace.rhs;
    auto& mutual = workspace.mutual;
    auto& gradient = workspace.mutual_gradient;
    workspace.resistance = derive_resistance(state);

    const int radial = armature_.radial_filaments();
    const double initial_center = armature_.position();
#pragma omp parallel for if (N_fil_ >= 8)
    for (int filament = 0; filament < N_fil_; ++filament) {
        const int axial = filament / radial + 1;
        const int radial_index = filament % radial + 1;
        const double relative = armature_.filament_axial_position(axial) - initial_center;
        const double separation = state.arm_position + relative - coil_.position();
        const double filament_length =
            armature_.length() / static_cast<double>(armature_.axial_filaments());
        mutual(filament) = physics::mutual_inductance_coil(
            coil_.inner_radius(), coil_.outer_radius(), coil_.length(), coil_.turns(),
            armature_.filament_inner_radius(radial_index),
            armature_.filament_outer_radius(radial_index), filament_length, 1,
            separation, 9, false);
        gradient(filament) = physics::mutual_inductance_gradient_coil(
            coil_.inner_radius(), coil_.outer_radius(), coil_.length(), coil_.turns(),
            armature_.filament_inner_radius(radial_index),
            armature_.filament_outer_radius(radial_index), filament_length, 1,
            separation, 9, false);
    }

    matrix.setZero();
    matrix(0, 0) = circuit_active ? L_d_ : 1.0;
    for (int filament = 0; filament < N_fil_; ++filament) {
        if (circuit_active) {
            matrix(0, filament + 1) = mutual(filament);
            matrix(filament + 1, 0) = mutual(filament);
        }
        matrix(filament + 1, filament + 1) = L_fil_(filament);
        for (int other = 0; other < N_fil_; ++other)
            if (filament != other) matrix(filament + 1, other + 1) = M_mat_(filament, other);
    }

    const double coil_current = state.currents(0);
    const double velocity = state.arm_velocity;
    double motional_emf = 0.0;
    for (int filament = 0; filament < N_fil_; ++filament)
        motional_emf += gradient(filament) * state.currents(filament + 1);
    rhs(0) = circuit_active
        ? excitation_->voltage(excitation) - R_d_ * coil_current - velocity * motional_emf
        : 0.0;
    for (int filament = 0; filament < N_fil_; ++filament) {
        rhs(filament + 1) = -workspace.resistance(filament) *
            state.currents(filament + 1) -
            velocity * gradient(filament) * coil_current;
    }

    DerivativeResult result;
    result.physical_derivative.currents.resize(N_fil_ + 1);
    Eigen::LDLT<Eigen::MatrixXd> solver(matrix);
    if (solver.info() == Eigen::Success)
        result.physical_derivative.currents = solver.solve(rhs);
    else
        result.physical_derivative.currents = matrix.colPivHouseholderQr().solve(rhs);
    if (!result.physical_derivative.currents.allFinite())
        throw std::runtime_error("single-stage circuit solve produced non-finite derivatives");

    result.force = circuit_active ? compute_force(state, gradient) : 0.0;
    result.mutual_gradient = gradient;
    result.physical_derivative.arm_position = state.arm_velocity;
    result.physical_derivative.arm_velocity = result.force / armature_.mass();
    if (enable_thermal_) {
        result.physical_derivative.filament_temperatures.resize(N_fil_);
        for (int filament = 0; filament < N_fil_; ++filament) {
            const double current = state.currents(filament + 1);
            const double temperature = state.filament_temperatures(filament);
            const double cp = physics::material_cp(armature_.material(), temperature);
            result.physical_derivative.filament_temperatures(filament) =
                current * current * workspace.resistance(filament) /
                (mass_fil_(filament) * cp);
        }
    }
    result.excitation_derivatives.push_back(
        excitation_->continuous_derivative(excitation, coil_current));
    return result;
}

template<typename SP>
IntegrationState SingleStageSim<SP>::make_trial(
    const IntegrationState& initial,
    const DerivativeResult& derivative,
    double scale) const {
    IntegrationState trial = clone_integration_state(initial);
    trial.physical += scale * derivative.physical_derivative;
    excitation_->advance_snapshot_derivative(
        *trial.excitations.front(), scale, derivative.excitation_derivatives.front());
    return trial;
}

template<typename SP>
void SingleStageSim<SP>::apply_excitation_events(
    IntegrationState& state, double coil_current) const {
    auto& snapshot = *state.excitations.front();
    auto& stage = state.stages.front();
    if (auto* capacitor = dynamic_cast<CapacitorSnapshot*>(&snapshot)) {
        if (capacitor->capacitor_voltage <= 0.0) {
            excitation_->apply_event(snapshot, ExcitationEvent::CapacitorZero);
            excitation_->apply_event(snapshot, ExcitationEvent::Finished);
        }
    } else if (auto* crowbar = dynamic_cast<CrowbarSnapshot*>(&snapshot)) {
        if (!crowbar->diode_on && crowbar->capacitor_voltage <= 0.0) {
            excitation_->apply_event(snapshot, ExcitationEvent::CapacitorZero);
            excitation_->apply_event(snapshot, ExcitationEvent::CrowbarOn);
        }
        stage.crowbar_on = crowbar->diode_on;
        if (crowbar->diode_on && std::abs(coil_current) < kQuietCurrent)
            excitation_->apply_event(snapshot, ExcitationEvent::Finished);
    } else if (auto* waveform = dynamic_cast<WaveformSnapshot*>(&snapshot)) {
        const auto* source = dynamic_cast<const WaveformExcitation*>(excitation_.get());
        if (source && waveform->time >= source->end_time())
            excitation_->apply_event(snapshot, ExcitationEvent::WaveformEnd);
    }

    stage.excitation_finished = snapshot_finished(snapshot);
    if (stage.excitation_finished && std::abs(coil_current) < kQuietCurrent)
        mark_stage_completed(state, 0);
}

template<typename SP>
double SingleStageSim<SP>::excitation_event_value(
    const ExcitationSnapshot& snapshot) const {
    if (const auto* capacitor = dynamic_cast<const CapacitorSnapshot*>(&snapshot))
        return capacitor->finished ? -1.0 : capacitor->capacitor_voltage;
    if (const auto* crowbar = dynamic_cast<const CrowbarSnapshot*>(&snapshot))
        return crowbar->diode_on ? 1.0 : crowbar->capacitor_voltage;
    if (const auto* waveform = dynamic_cast<const WaveformSnapshot*>(&snapshot)) {
        const auto* source = dynamic_cast<const WaveformExcitation*>(excitation_.get());
        return source ? source->end_time() - waveform->time : 1.0;
    }
    return 1.0;
}

template<typename SP>
IntegrationState SingleStageSim<SP>::advance_euler(
    const IntegrationState& pre, double dt) {
    const auto derivative = evaluate_derivatives(
        pre.physical, *pre.excitations.front(), workspace_,
        pre.stages.front().circuit_active);
    pre_step_mutual_gradient_ = derivative.mutual_gradient;
    auto post = make_trial(pre, derivative, dt);
    apply_excitation_events(post, post.physical.currents(0));
    return post;
}

template<typename SP>
IntegrationState SingleStageSim<SP>::advance_rk4_segment(
    const IntegrationState& pre, double dt) {
    const auto k1 = evaluate_derivatives(pre.physical, *pre.excitations.front(), workspace_,
                                         pre.stages.front().circuit_active);
    if (pre_step_mutual_gradient_.size() == 0)
        pre_step_mutual_gradient_ = k1.mutual_gradient;
    const auto s2 = make_trial(pre, k1, 0.5 * dt);
    const auto k2 = evaluate_derivatives(s2.physical, *s2.excitations.front(), workspace_,
                                         s2.stages.front().circuit_active);
    const auto s3 = make_trial(pre, k2, 0.5 * dt);
    const auto k3 = evaluate_derivatives(s3.physical, *s3.excitations.front(), workspace_,
                                         s3.stages.front().circuit_active);
    const auto s4 = make_trial(pre, k3, dt);
    const auto k4 = evaluate_derivatives(s4.physical, *s4.excitations.front(), workspace_,
                                         s4.stages.front().circuit_active);

    auto post = clone_integration_state(pre);
    auto weighted = k1.physical_derivative;
    weighted += 2.0 * k2.physical_derivative;
    weighted += 2.0 * k3.physical_derivative;
    weighted += k4.physical_derivative;
    post.physical += (dt / 6.0) * weighted;
    ExcitationDerivative derivative;
    derivative.capacitor_voltage_rate =
        (k1.excitation_derivatives[0].capacitor_voltage_rate +
         2.0 * k2.excitation_derivatives[0].capacitor_voltage_rate +
         2.0 * k3.excitation_derivatives[0].capacitor_voltage_rate +
         k4.excitation_derivatives[0].capacitor_voltage_rate) / 6.0;
    derivative.waveform_time_rate =
        (k1.excitation_derivatives[0].waveform_time_rate +
         2.0 * k2.excitation_derivatives[0].waveform_time_rate +
         2.0 * k3.excitation_derivatives[0].waveform_time_rate +
         k4.excitation_derivatives[0].waveform_time_rate) / 6.0;
    excitation_->advance_snapshot_derivative(*post.excitations.front(), dt, derivative);
    return post;
}

template<typename SP>
IntegrationState SingleStageSim<SP>::advance_rk4_event_aware(
    const IntegrationState& pre, double dt) {
    pre_step_mutual_gradient_.resize(0);
    auto current = clone_integration_state(pre);
    double remaining = dt;
    for (int segment = 0; segment < kMaxEventSegments && remaining > 0.0; ++segment) {
        auto candidate = advance_rk4_segment(current, remaining);
        const double start_value = excitation_event_value(*current.excitations.front());
        const double end_value = excitation_event_value(*candidate.excitations.front());
        if (!(start_value > 0.0 && end_value <= 0.0)) {
            apply_excitation_events(candidate, candidate.physical.currents(0));
            return candidate;
        }

        double lower = 0.0;
        double upper = remaining;
        for (int iteration = 0; iteration < kMaxEventBisections; ++iteration) {
            const double midpoint = 0.5 * (lower + upper);
            const auto trial = advance_rk4_segment(current, midpoint);
            if (excitation_event_value(*trial.excitations.front()) > 0.0)
                lower = midpoint;
            else
                upper = midpoint;
            if (upper - lower <= kEventTimeTolerance) break;
        }
        if (upper <= 0.0)
            throw std::runtime_error("zero-duration repeated excitation event");
        current = advance_rk4_segment(current, upper);
        apply_excitation_events(current, current.physical.currents(0));
        remaining -= upper;
    }
    if (remaining > kEventTimeTolerance)
        throw std::runtime_error("single-stage RK4 event segment limit exceeded");
    return current;
}

template<typename SP>
void SingleStageSim<SP>::record_step(double post_time) {
    SimStep entry;
    entry.time = post_time;
    entry.cap_voltage = excitation_->voltage(*integration_state_.excitations.front());
    entry.coil_current = integration_state_.physical.currents(0);
    entry.arm_position = integration_state_.physical.arm_position;
    entry.arm_velocity = integration_state_.physical.arm_velocity;
    entry.force = 0.0;
    if (integration_state_.stages.front().circuit_active) {
        const double coil_current = integration_state_.physical.currents(0);
        for (int filament = 0; filament < N_fil_; ++filament)
            entry.force += coil_current *
                integration_state_.physical.currents(filament + 1) *
                pre_step_mutual_gradient_(filament);
    }
    entry.filament_currents.resize(N_fil_);
    for (int filament = 0; filament < N_fil_; ++filament)
        entry.filament_currents[filament] = integration_state_.physical.currents(filament + 1);
    if (enable_thermal_) {
        entry.filament_temperatures.resize(N_fil_);
        for (int filament = 0; filament < N_fil_; ++filament)
            entry.filament_temperatures[filament] =
                integration_state_.physical.filament_temperatures(filament);
    }
    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool SingleStageSim<SP>::check_termination(const TerminationPolicy& policy) {
    const auto& state = integration_state_.physical;
    if (policy.enable_bound_check && state.arm_position >= policy.barrel_end_position) return true;
    if (step_count_ >= policy.max_steps) return true;
    if (integration_state_.stages.front().circuit_active) return false;
    if (integration_state_.stages.front().stage_completed) return true;
    if (policy.enable_velocity_check && step_count_ >= policy.velocity_decay_steps) {
        const auto derivative = evaluate_derivatives(
            state, *integration_state_.excitations.front(), workspace_,
            integration_state_.stages.front().circuit_active);
        bool decaying = true;
        const auto count = static_cast<int>(result_.history.size());
        for (int index = 0; index < policy.velocity_decay_steps; ++index) {
            if (count - 2 - index < 0 ||
                result_.history[count - 1 - index].arm_velocity >=
                    result_.history[count - 2 - index].arm_velocity) {
                decaying = false;
                break;
            }
        }
        if (decaying && std::abs(derivative.physical_derivative.arm_velocity) <
                            policy.accel_threshold)
            return true;
    }
    return false;
}

template<typename SP>
void SingleStageSim<SP>::prepare_summary() {
    auto& summary = result_.summary;
    summary = {};
    summary.step_count = step_count_;
    if (result_.history.empty()) return;
    summary.total_time = result_.history.back().time;
    summary.muzzle_velocity = result_.history.back().arm_velocity;
    for (const auto& step : result_.history) {
        summary.max_force = std::max(summary.max_force, std::abs(step.force));
        summary.peak_coil_current =
            std::max(summary.peak_coil_current, std::abs(step.coil_current));
    }
    if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(excitation_.get())) {
        const double input = 0.5 * capacitor->capacitance() *
            capacitor->initial_voltage() * capacitor->initial_voltage();
        const double output = 0.5 * armature_.mass() *
            summary.muzzle_velocity * summary.muzzle_velocity;
        summary.efficiency = input > 0.0 ? output / input : 0.0;
    }
}

template<typename SP>
const SimStep& SingleStageSim<SP>::step() {
    const IntegrationState pre = clone_integration_state(integration_state_);
    if constexpr (std::is_same_v<SP, RK4Stepper>)
        integration_state_ = advance_rk4_event_aware(pre, dt_);
    else
        integration_state_ = advance_euler(pre, dt_);
    excitation_->restore(*integration_state_.excitations.front());
    record_step((step_count_ + 1) * dt_);
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const SimResult& SingleStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const SimResult& SingleStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy)) step();
    prepare_summary();
    return result_;
}

template<typename SP>
void SingleStageSim<SP>::reset() {
    integration_state_ = clone_integration_state(initial_integration_state_);
    excitation_->restore(*integration_state_.excitations.front());
    result_ = {};
    step_count_ = 0;
}

template class SingleStageSim<EulerStepper>;
template class SingleStageSim<RK4Stepper>;

} // namespace coilgun::simulation
