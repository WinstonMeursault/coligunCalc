/**
 * @file single_stage_sim.cpp
 * @brief Single-stage coilgun simulation engine implementation.
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.3.2, Sec.5, Sec.8.
 */

#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/constants.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <utility>

namespace coilgun::simulation {

SimState& SimState::operator+=(const SimState& rhs) {
    currents += rhs.currents;
    arm_position += rhs.arm_position;
    arm_velocity += rhs.arm_velocity;
    if (filament_temperatures.size() > 0 && rhs.filament_temperatures.size() > 0)
        filament_temperatures += rhs.filament_temperatures;
    return *this;
}

SimState& SimState::operator*=(double scalar) {
    currents *= scalar;
    arm_position *= scalar;
    arm_velocity *= scalar;
    filament_temperatures *= scalar;
    return *this;
}

SimState operator+(SimState lhs, const SimState& rhs) { lhs += rhs; return lhs; }
SimState operator*(double scalar, SimState s) { s *= scalar; return s; }

namespace {

double filament_axial_sep(const components::Armature& arm, int a_idx, int b_idx) {
    int nr = arm.radial_filaments();
    int ia = a_idx / nr + 1;
    int ib = b_idx / nr + 1;
    return std::abs(arm.filament_axial_position(ia) - arm.filament_axial_position(ib));
}

} // anonymous namespace

template<typename SP>
SingleStageSim<SP>::SingleStageSim(
        components::DrivingCoil coil,
        components::Armature    armature,
        std::unique_ptr<Excitation> excitation,
        double                     dt,
        bool                       enable_thermal)
    : coil_(std::move(coil)), armature_(std::move(armature)),
      excitation_(std::move(excitation)),
      dt_(dt), enable_thermal_(enable_thermal)
{
    N_fil_ = armature_.total_filaments();
    R_d_   = coil_.resistance();
    L_d_   = coil_.self_inductance();

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

    build_filament_M_matrix();

    M1_.resize(N); dM1_.resize(N);
    L_total_.resize(N + 1, N + 1);
    RHS_.resize(N + 1);

    state_.currents.resize(N + 1);
    state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    if (enable_thermal_) {
        state_.filament_temperatures.resize(N);
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    }

}

template<typename SP>
void SingleStageSim<SP>::build_filament_M_matrix() {
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
SimState SingleStageSim<SP>::compute_derivatives(const SimState& s) {
    SimState ds;
    ds.currents.resize(N_fil_ + 1);
    ds.arm_position = s.arm_velocity;
    ds.arm_velocity = compute_force(s) / armature_.mass();
    if (enable_thermal_ && s.filament_temperatures.size() > 0) {
        ds.filament_temperatures.resize(N_fil_);
        ds.filament_temperatures.setZero();
    }

    int nr = armature_.radial_filaments();
    double arm_center = s.arm_position;
    double arm_center_init = armature_.position();

#pragma omp parallel for
    for (int k = 0; k < N_fil_; ++k) {
        int i = k / nr + 1;
        int j = k % nr + 1;
        double z_rel = armature_.filament_axial_position(i) - arm_center_init;
        double z_global = arm_center + z_rel;
        double sep = z_global - coil_.position();
        double fil_ri = armature_.filament_inner_radius(j);
        double fil_re = armature_.filament_outer_radius(j);
        double fil_l  = armature_.length() / armature_.axial_filaments();
        M1_(k) = physics::mutual_inductance_coil(
            coil_.inner_radius(), coil_.outer_radius(),
            coil_.length(), coil_.turns(),
            fil_ri, fil_re, fil_l, 1, sep, 9, false);
        dM1_(k) = physics::mutual_inductance_gradient_coil(
            coil_.inner_radius(), coil_.outer_radius(),
            coil_.length(), coil_.turns(),
            fil_ri, fil_re, fil_l, 1, sep, 9, false);
    }

    L_total_.setZero();
    L_total_(0, 0) = L_d_;
    for (int k = 0; k < N_fil_; ++k) {
        L_total_(0, k + 1) = M1_(k);
        L_total_(k + 1, 0) = M1_(k);
    }
    for (int a = 0; a < N_fil_; ++a) {
        L_total_(a + 1, a + 1) = L_fil_(a);
        for (int b = 0; b < N_fil_; ++b) {
            if (a != b) L_total_(a + 1, b + 1) = M_mat_(a, b);
        }
    }

    double U   = excitation_->voltage();
    double I_d = s.currents(0);
    double v   = s.arm_velocity;

    double motional_emf = 0.0;
    for (int k = 0; k < N_fil_; ++k)
        motional_emf += dM1_(k) * s.currents(k + 1);

    RHS_(0) = U - R_d_ * I_d - v * motional_emf;
    for (int k = 0; k < N_fil_; ++k) {
        RHS_(k + 1) = -R_fil_(k) * s.currents(k + 1) - v * dM1_(k) * I_d;
    }

    Eigen::LDLT<Eigen::MatrixXd> solver(L_total_);
    if (solver.info() != Eigen::Success)
        ds.currents = L_total_.colPivHouseholderQr().solve(RHS_);
    else
        ds.currents = solver.solve(RHS_);
    return ds;
}

template<typename SP>
double SingleStageSim<SP>::compute_force(const SimState& s) {
    double F = 0.0, I_d = s.currents(0);
#pragma omp parallel for reduction(+:F)
    for (int k = 0; k < N_fil_; ++k)
        F += I_d * s.currents(k + 1) * dM1_(k);
    return F;
}

template<typename SP>
void SingleStageSim<SP>::update_temperatures(SimState& s, double dt_sub) {
    auto mat = armature_.material();
    double beta = physics::material_beta(mat);
#pragma omp parallel for
    for (int k = 0; k < N_fil_; ++k) {
        double I_k = s.currents(k + 1);
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
void SingleStageSim<SP>::record_step(double cap_voltage) {
    SimStep entry;
    entry.time         = step_count_ * dt_;
    entry.cap_voltage  = cap_voltage;
    entry.coil_current = state_.currents(0);
    entry.arm_position = state_.arm_position;
    entry.arm_velocity = state_.arm_velocity;
    entry.force        = compute_force(state_);
    entry.filament_currents.resize(N_fil_);
    for (int k = 0; k < N_fil_; ++k)
        entry.filament_currents[k] = state_.currents(k + 1);
    if (enable_thermal_ && state_.filament_temperatures.size() > 0) {
        entry.filament_temperatures.resize(N_fil_);
        for (int k = 0; k < N_fil_; ++k)
            entry.filament_temperatures[k] = state_.filament_temperatures(k);
    }
    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool SingleStageSim<SP>::check_termination(const TerminationPolicy& policy) {
    if (policy.enable_bound_check && state_.arm_position >= policy.barrel_end_position)
        return true;
    if (step_count_ >= policy.max_steps) return true;
    if (policy.enable_velocity_check && step_count_ >= policy.velocity_decay_steps) {
        double accel = compute_force(state_) / armature_.mass();
        const auto& hist = result_.history;
        bool decaying = true;
        auto n = static_cast<int>(hist.size());
        for (int i = 0; i < policy.velocity_decay_steps; ++i) {
            if (n - 2 - i < 0 || hist[n - 1 - i].arm_velocity >= hist[n - 2 - i].arm_velocity) {
                decaying = false; break;
            }
        }
        if (decaying && std::abs(accel) < policy.accel_threshold) return true;
    }
    return false;
}

template<typename SP>
void SingleStageSim<SP>::prepare_summary() {
    auto& s = result_.summary;
    s.step_count = step_count_;
    if (result_.history.empty()) return;
    const auto& last = result_.history.back();
    s.total_time      = last.time;
    s.muzzle_velocity = last.arm_velocity;
    for (const auto& step : result_.history) {
        if (step.force > s.max_force) s.max_force = step.force;
        if (step.coil_current > s.peak_coil_current)
            s.peak_coil_current = step.coil_current;
    }
    auto* cap = dynamic_cast<CapacitorExcitation*>(excitation_.get());
    if (cap) {
        double E_in = 0.5 * cap->capacitance()
                      * cap->initial_voltage()
                      * cap->initial_voltage();
        double E_out = 0.5 * armature_.mass() * s.muzzle_velocity * s.muzzle_velocity;
        s.efficiency = (E_in > 0.0) ? E_out / E_in : 0.0;
    }
}

template<typename SP>
const SimStep& SingleStageSim<SP>::step() {
    state_ = stepper_.advance(dt_, state_,
        [this](const SimState& s) { return compute_derivatives(s); });

    excitation_->advance(dt_, state_.currents(0));

    if (enable_thermal_) update_temperatures(state_, dt_);

    record_step(excitation_->voltage());
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const SimResult& SingleStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const SimResult& SingleStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy) && !excitation_->finished())
        step();
    prepare_summary();
    return result_;
}

template<typename SP>
void SingleStageSim<SP>::reset() {
    state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    if (enable_thermal_ && state_.filament_temperatures.size() > 0)
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    R_fil_ = R_fil_ref_;
    result_ = SimResult{};
    step_count_ = 0;
    excitation_->reset();
}

template class SingleStageSim<EulerStepper>;
template class SingleStageSim<RK4Stepper>;

} // namespace coilgun::simulation
