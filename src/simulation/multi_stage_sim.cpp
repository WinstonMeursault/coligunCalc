/**
 * @file multi_stage_sim.cpp
 * @brief Multi-stage coilgun simulation engine implementation.
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.3.3, Sec.3.4, Sec.5, Sec.8.
 */

#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/constants.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace coilgun::simulation {

MultiStageState& MultiStageState::operator+=(const MultiStageState& rhs) {
    currents += rhs.currents;
    arm_position += rhs.arm_position;
    arm_velocity += rhs.arm_velocity;
    if (filament_temperatures.size() > 0 && rhs.filament_temperatures.size() > 0)
        filament_temperatures += rhs.filament_temperatures;
    return *this;
}

MultiStageState& MultiStageState::operator*=(double scalar) {
    currents *= scalar;
    arm_position *= scalar;
    arm_velocity *= scalar;
    filament_temperatures *= scalar;
    return *this;
}

MultiStageState operator+(MultiStageState lhs, const MultiStageState& rhs) { lhs += rhs; return lhs; }
MultiStageState operator*(double scalar, MultiStageState s) { s *= scalar; return s; }

namespace {

double filament_axial_sep(const components::Armature& arm, int a_idx, int b_idx) {
    int nr = arm.radial_filaments();
    int ia = a_idx / nr + 1;
    int ib = b_idx / nr + 1;
    return std::abs(arm.filament_axial_position(ia) - arm.filament_axial_position(ib));
}

} // anonymous namespace

template<typename SP>
MultiStageSim<SP>::MultiStageSim(
        std::vector<components::DrivingCoil>     coils,
        components::Armature                      armature,
        std::vector<std::unique_ptr<Excitation>>  excitations,
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal,
        OptimizationLevel                         opt_level)
    : n_stages_(static_cast<int>(coils.size()))
    , coils_(std::move(coils))
    , armature_(std::move(armature))
    , excitations_(std::move(excitations))
    , trigger_configs_(std::move(trigger_configs))
    , dt_(dt)
    , enable_thermal_(enable_thermal)
    , opt_level_(opt_level)
    , N_fil_(armature_.total_filaments())
{
    // ---- validation ----
    if (n_stages_ <= 0) {
        throw std::invalid_argument("MultiStageSim: at least one stage is required");
    }
    if (n_stages_ != static_cast<int>(excitations_.size())) {
        throw std::invalid_argument(
            "MultiStageSim: coils.size() != excitations.size()");
    }
    if (n_stages_ > kMaxStages) {
        throw std::invalid_argument(
            "MultiStageSim: n_stages exceeds kMaxStages");
    }
    if (static_cast<int>(trigger_configs_.size()) != n_stages_ - 1) {
        throw std::invalid_argument(
            "MultiStageSim: trigger_configs.size() must be n_stages-1");
    }
    for (const auto& config : trigger_configs_)
        validate_trigger_config(config);

    triggered_.resize(n_stages_, false);
    finished_.resize(n_stages_, false);
    trigger_times_.resize(n_stages_, 0.0);
    trigger_positions_.resize(n_stages_, armature_.position());
    initial_stage_energies_.resize(n_stages_, 0.0);

    // Stage 0 auto-triggers at t=0
    triggered_[0] = true;
    for (int i = 0; i < n_stages_; ++i) {
        if (const auto* capacitor = dynamic_cast<const CapacitorExcitation*>(excitations_[i].get()))
            initial_stage_energies_[i] = 0.5 * capacitor->capacitance() *
                capacitor->initial_voltage() * capacitor->initial_voltage();
    }

    // ---- extract coil self-inductances and resistances ----
    R_diag_.resize(n_stages_);
    L_diag_.resize(n_stages_);
    for (int i = 0; i < n_stages_; ++i) {
        R_diag_(i) = coils_[i].resistance();
        L_diag_(i) = coils_[i].self_inductance();
    }

    // ---- extract filament arrays ----
    int N = N_fil_;
    R_fil_ref_.resize(N); R_fil_.resize(N);
    L_fil_.resize(N); mass_fil_.resize(N);
    const auto& R_arm = armature_.resistances();
    const auto& L_arm = armature_.inductances();
    const auto& M_arm = armature_.masses();
    for (int k = 0; k < N; ++k) {
        R_fil_ref_(k) = R_arm[k];
        R_fil_(k)     = R_arm[k];
        L_fil_(k)     = L_arm[k];
        mass_fil_(k)  = M_arm[k];
    }

    // ---- build inter-filament and inter-coil mutual inductance matrices ----
    build_filament_M_matrix();
    precompute_M_cc();

    // ---- allocate working buffers ----
    M1_mat_.resize(n_stages_, N);
    dM1_mat_.resize(n_stages_, N);

    int dim = n_stages_ + N;
    L_total_.resize(dim, dim);
    RHS_.resize(dim);

    state_.currents.resize(dim);
    state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    if (enable_thermal_) {
        state_.filament_temperatures.resize(N);
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    }

}

template<typename SP>
void MultiStageSim<SP>::build_filament_M_matrix() {
    int N = N_fil_;
    M_mat_.resize(N, N);
    M_mat_.setZero();
    int nr = armature_.radial_filaments();
    for (int a = 0; a < N; ++a) {
        int ja = a % nr + 1;
        double ra = armature_.filament_mean_radius(ja);
        for (int b = a + 1; b < N; ++b) {
            int jb = b % nr + 1;
            double rb = armature_.filament_mean_radius(jb);
            double sep = filament_axial_sep(armature_, a, b);
            double m = physics::mutual_inductance_filament(ra, rb, sep, true);
            M_mat_(a, b) = m;
            M_mat_(b, a) = m;
        }
    }
}

template<typename SP>
void MultiStageSim<SP>::precompute_M_cc() {
    M_cc_.resize(n_stages_, n_stages_);
    M_cc_.setZero();
    for (int i = 0; i < n_stages_; ++i) {
        for (int j = i + 1; j < n_stages_; ++j) {
            const auto& a = coils_[i];
            const auto& b = coils_[j];
            double sep = std::abs(a.position() - b.position());
            double m = physics::mutual_inductance_coil(
                a.inner_radius(), a.outer_radius(), a.length(), a.turns(),
                b.inner_radius(), b.outer_radius(), b.length(), b.turns(),
                sep, 9, true);
            M_cc_(i, j) = m;
            M_cc_(j, i) = m;
        }
    }
}

template<typename SP>
bool MultiStageSim<SP>::is_stage_within_range(int stage_idx) const {
    if (opt_level_ != OptimizationLevel::Full) return true;
    const auto& coil = coils_[stage_idx];
    double cutoff = 10.0 * coil.length();
    double dist = std::abs(state_.arm_position - coil.position());
    return dist <= cutoff;
}

template<typename SP>
void MultiStageSim<SP>::build_system_matrix(const MultiStageState& s) {
    int S = n_stages_;
    int F = N_fil_;
    int dim = S + F;
    L_total_.setZero();

    // coil-coil block
    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) {
            L_total_(i, i) = 1.0;
            continue;
        }
        L_total_(i, i) = L_diag_(i);
        for (int j = 0; j < S; ++j) {
            if (i != j && triggered_[j] && !finished_[j])
                L_total_(i, j) = M_cc_(i, j);
        }
    }

    // coil-filament and filament-coil blocks
    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        for (int k = 0; k < F; ++k) {
            L_total_(i, S + k) = M1_mat_(i, k);
            L_total_(S + k, i) = M1_mat_(i, k);
        }
    }

    // filament block: self-inductances on diagonal, inter-filament mutual off-diagonal
    for (int a = 0; a < F; ++a) {
        L_total_(S + a, S + a) = L_fil_(a);
        for (int b = 0; b < F; ++b) {
            if (a != b) L_total_(S + a, S + b) = M_mat_(a, b);
        }
    }
}

template<typename SP>
MultiStageState MultiStageSim<SP>::compute_derivatives(const MultiStageState& s) {
    int S = n_stages_;
    int F = N_fil_;
    int dim = S + F;

    MultiStageState ds;
    ds.currents.resize(dim);
    ds.arm_position = s.arm_velocity;
    if (enable_thermal_ && s.filament_temperatures.size() > 0) {
        ds.filament_temperatures.resize(F);
        ds.filament_temperatures.setZero();
    }

    std::vector<int> active_idx;
    active_idx.reserve(S);
    for (int st = 0; st < S; ++st) {
        if (triggered_[st] && !finished_[st] && is_stage_within_range(st))
            active_idx.push_back(st);
    }
    int n_active = static_cast<int>(active_idx.size());

    int nr = armature_.radial_filaments();
    double arm_center = s.arm_position;
    double arm_center_init = armature_.position();

    M1_mat_.setZero();
    dM1_mat_.setZero();

    if (n_active > 0) {
#pragma omp parallel for
        for (int flat = 0; flat < n_active * F; ++flat) {
            int si = active_idx[flat / F];
            int fi = flat % F;

            const auto& coil = coils_[si];
            int i = fi / nr + 1;
            int j = fi % nr + 1;
            double z_rel = armature_.filament_axial_position(i) - arm_center_init;
            double z_global = arm_center + z_rel;
            double sep = z_global - coil.position();
            double fil_ri = armature_.filament_inner_radius(j);
            double fil_re = armature_.filament_outer_radius(j);
            double fil_l  = armature_.length() / armature_.axial_filaments();

            double dist = std::abs(arm_center - coil.position());
            int n_nodes = 9;
            if (opt_level_ == OptimizationLevel::Full && dist > coil.length())
                n_nodes = 4;

            M1_mat_(si, fi) = physics::mutual_inductance_coil(
                coil.inner_radius(), coil.outer_radius(),
                coil.length(), coil.turns(),
                fil_ri, fil_re, fil_l, 1, sep, n_nodes, false);

            dM1_mat_(si, fi) = physics::mutual_inductance_gradient_coil(
                coil.inner_radius(), coil.outer_radius(),
                coil.length(), coil.turns(),
                fil_ri, fil_re, fil_l, 1, sep, n_nodes, false);
        }
    }

    build_system_matrix(s);

    ds.arm_velocity = compute_force(s) / armature_.mass();

    double v = s.arm_velocity;
    RHS_.setZero();

    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) continue;

        double U   = excitations_[i]->voltage();
        double I_d = s.currents(i);

        double motional_emf = 0.0;
        for (int k = 0; k < F; ++k)
            motional_emf += dM1_mat_(i, k) * s.currents(S + k);

        RHS_(i) = U - R_diag_(i) * I_d - v * motional_emf;
    }

    for (int k = 0; k < F; ++k) {
        double I_f = s.currents(S + k);
        double coil_back_emf = 0.0;
        for (int i = 0; i < S; ++i) {
            if (triggered_[i] && !finished_[i])
                coil_back_emf += dM1_mat_(i, k) * s.currents(i);
        }
        RHS_(S + k) = -R_fil_(k) * I_f - v * coil_back_emf;
    }

    Eigen::LDLT<Eigen::MatrixXd> solver(L_total_);
    if (solver.info() != Eigen::Success)
        ds.currents = L_total_.colPivHouseholderQr().solve(RHS_);
    else
        ds.currents = solver.solve(RHS_);
    return ds;
}

template<typename SP>
double MultiStageSim<SP>::compute_force(const MultiStageState& s) {
    int S = n_stages_;
    int F = N_fil_;
    double F_net = 0.0;
#pragma omp parallel for reduction(+:F_net)
    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        double I_d = s.currents(i);
        for (int k = 0; k < F; ++k)
            F_net += I_d * s.currents(S + k) * dM1_mat_(i, k);
    }
    return F_net;
}

template<typename SP>
void MultiStageSim<SP>::update_temperatures(MultiStageState& s, double dt_sub) {
    int F = N_fil_;
    auto mat = armature_.material();
    double beta = physics::material_beta(mat);
#pragma omp parallel for
    for (int k = 0; k < F; ++k) {
        double I_k = s.currents(n_stages_ + k);
        double m_k = mass_fil_(k);
        double T_k = s.filament_temperatures(k);
        double cp  = physics::material_cp(mat, T_k);
        double dT  = I_k * I_k * R_fil_(k) * dt_sub / (m_k * cp);
        T_k += dT;
        s.filament_temperatures(k) = T_k;
        R_fil_(k) = R_fil_ref_(k) * (1.0 + beta * (T_k - physics::T_REFERENCE));
    }
}

template<typename SP>
void MultiStageSim<SP>::check_triggers() {
    if (n_stages_ <= 1) return;

    double current_time = step_count_ * dt_;

    for (int i = 1; i < n_stages_; ++i) {
        if (triggered_[i]) continue;
        if (!triggered_[i - 1]) continue;

        const auto& cfg = trigger_configs_[i - 1];
        bool fire = false;

        if (cfg.mode == TriggerMode::Position) {
            fire = (state_.arm_position >= cfg.value);
        } else {
            // TimeDelay: relative to previous stage's trigger time
            fire = (current_time >= trigger_times_[i - 1] + cfg.value);
        }

        if (fire) {
            triggered_[i] = true;
            trigger_times_[i] = current_time;
            trigger_positions_[i] = state_.arm_position;
        }
    }
}

template<typename SP>
void MultiStageSim<SP>::extinguish_quiet_stages() {
    static constexpr double kQuietThreshold = 1e-6; // A
    for (int i = 0; i < n_stages_; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        if (excitations_[i]->voltage() == 0.0 &&
            std::abs(state_.currents(i)) < kQuietThreshold) {
            finished_[i] = true;
        }
    }
}

template<typename SP>
void MultiStageSim<SP>::record_step() {
    int S = n_stages_;
    int F = N_fil_;

    MultiStageStep entry;
    entry.state.time         = step_count_ * dt_;
    entry.state.arm_position = state_.arm_position;
    entry.state.arm_velocity = state_.arm_velocity;
    entry.state.force        = compute_force(state_);
    entry.state.filament_currents.resize(F);
    for (int k = 0; k < F; ++k)
        entry.state.filament_currents[k] = state_.currents(S + k);

    if (enable_thermal_ && state_.filament_temperatures.size() > 0) {
        entry.state.filament_temperatures.resize(F);
        for (int k = 0; k < F; ++k)
            entry.state.filament_temperatures[k] = state_.filament_temperatures(k);
    }

    entry.cap_voltages.resize(S);
    entry.coil_currents.resize(S);
    entry.stage_forces.resize(S);
    for (int i = 0; i < S; ++i) {
        entry.cap_voltages[i] = triggered_[i] ? excitations_[i]->voltage() : 0.0;
        entry.coil_currents[i] = triggered_[i] ? state_.currents(i) : 0.0;
        if (!triggered_[i] || finished_[i]) continue;
        for (int k = 0; k < F; ++k)
            entry.stage_forces[i] += state_.currents(i) * state_.currents(S + k) * dM1_mat_(i, k);
    }
    entry.state.force = 0.0;
    for (const double force : entry.stage_forces) entry.state.force += force;

    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool MultiStageSim<SP>::check_all_finished() const {
    // A finite trigger policy remains eligible even if earlier stages finished;
    // only completed stages and explicit +infinity policies are terminal.
    auto terminally_ineligible = [this](int stage_idx) {
        for (int stage = stage_idx; stage > 0; --stage) {
            const auto& cfg = trigger_configs_[static_cast<std::size_t>(stage - 1)];
            if (cfg.value == INFINITY) return true;
            if (triggered_[stage - 1]) return false;
        }
        return false;
    };
    for (int i = 0; i < n_stages_; ++i) {
        if (finished_[i]) continue;
        if (!triggered_[i] && terminally_ineligible(i)) continue;
        return false;
    }
    return true;
}

template<typename SP>
bool MultiStageSim<SP>::check_termination(const TerminationPolicy& policy) {
    if (policy.enable_bound_check && state_.arm_position >= policy.barrel_end_position)
        return true;
    if (step_count_ >= policy.max_steps) return true;
    if (check_all_finished()) return true;

    if (policy.enable_velocity_check && step_count_ >= policy.velocity_decay_steps) {
        double accel = compute_force(state_) / armature_.mass();
        const auto& hist = result_.history;
        bool decaying = true;
        auto n = static_cast<int>(hist.size());
        for (int i = 0; i < policy.velocity_decay_steps; ++i) {
            if (n - 2 - i < 0 || hist[n - 1 - i].state.arm_velocity >= hist[n - 2 - i].state.arm_velocity) {
                decaying = false; break;
            }
        }
        if (decaying && std::abs(accel) < policy.accel_threshold) return true;
    }
    return false;
}

template<typename SP>
void MultiStageSim<SP>::prepare_summary() {
    int S = n_stages_;
    int F = N_fil_;
    auto& s = result_.summary;
    s = MultiStageSummary{};
    s.step_count = step_count_;
    if (result_.history.empty()) return;

    const auto& last = result_.history.back();
    s.total_time      = last.state.time;
    s.muzzle_velocity = last.state.arm_velocity;

    // per-stage summaries
    for (int i = 0; i < S; ++i) {
        if (!triggered_[i]) continue;
        PerStageSummary ps;
        ps.stage_index = i;
        ps.trigger_time = trigger_times_[i];

        ps.trigger_position = trigger_positions_[static_cast<std::size_t>(i)];

        auto* cap = dynamic_cast<CapacitorExcitation*>(excitations_[i].get());
        bool active = false;

        for (const auto& step : result_.history) {
            if (step.coil_currents[i] > ps.peak_current)
                ps.peak_current = step.coil_currents[i];
            if (i < static_cast<int>(step.stage_forces.size()))
                ps.max_force = std::max(ps.max_force, std::abs(step.stage_forces[i]));

            if (triggered_[i] && step.coil_currents[i] > 1e-6)
                active = true;
            if (active) ps.step_count_active++;

        }

        if (cap) {
            double E_end = 0.5 * cap->capacitance() * cap->capacitor_voltage() * cap->capacitor_voltage();
            ps.energy_depleted = initial_stage_energies_[static_cast<std::size_t>(i)] - E_end;
        }

        s.per_stage.push_back(ps);
    }

    // global
    for (const auto& step : result_.history) {
        if (std::abs(step.state.force) > s.max_force)
            s.max_force = std::abs(step.state.force);
    }
    for (int i = 0; i < S; ++i) {
        for (const auto& step : result_.history) {
            if (step.coil_currents[i] > s.peak_coil_current)
                s.peak_coil_current = step.coil_currents[i];
        }
    }

    // efficiency
    double E_total = 0.0;
    for (int i = 0; i < S; ++i) {
        auto* cap = dynamic_cast<CapacitorExcitation*>(excitations_[i].get());
        if (cap)
            E_total += 0.5 * cap->capacitance() * cap->initial_voltage() * cap->initial_voltage();
    }
    double E_kin = 0.5 * armature_.mass() * s.muzzle_velocity * s.muzzle_velocity;
    s.efficiency = (E_total > 0.0) ? E_kin / E_total : 0.0;
}

template<typename SP>
const MultiStageStep& MultiStageSim<SP>::step() {
    check_triggers();
    extinguish_quiet_stages();
    state_ = stepper_.advance(dt_, state_,
        [this](const MultiStageState& s) { return compute_derivatives(s); });

    for (int i = 0; i < n_stages_; ++i) {
        if (triggered_[i] && !finished_[i]) {
            excitations_[i]->advance(dt_, state_.currents(i));
            if (excitations_[i]->finished()) finished_[i] = true;
        }
    }

    if (enable_thermal_) update_temperatures(state_, dt_);

    record_step();
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const MultiStageResult& MultiStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const MultiStageResult& MultiStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy))
        step();
    prepare_summary();
    return result_;
}

template<typename SP>
void MultiStageSim<SP>::reset() {
    state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    if (enable_thermal_ && state_.filament_temperatures.size() > 0)
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);

    R_fil_ = R_fil_ref_;
    result_ = MultiStageResult{};
    step_count_ = 0;

    for (int i = 0; i < n_stages_; ++i) {
        triggered_[i] = (i == 0);
        finished_[i] = false;
        trigger_times_[i] = 0.0;
        trigger_positions_[i] = armature_.position();
        excitations_[i]->reset();
    }
}

template class MultiStageSim<EulerStepper>;
template class MultiStageSim<RK4Stepper>;

} // namespace coilgun::simulation
