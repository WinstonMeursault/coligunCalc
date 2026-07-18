/**
 * @file gpu_single_stage_sim.cu
 * @brief GPU-accelerated single-stage coilgun simulation — implementation.
 * @author Winston Meursault
 *
 * Uses CUDA kernels for the 4D Gauss-Legendre mutual inductance
 * integration. Linear system solve (Eigen LDLT) and kinematic/thermal
 * updates remain on the CPU identically to the single_stage_sim.cpp path.
 *
 * @see single_stage_sim.cpp
 */

#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/constants.hpp"
#include <atomic>
#include <cmath>
#include <iostream>
#include <utility>

namespace coilgun::physics {
    void upload_gl_nodes(int n_nodes);
    __global__ void mutual_inductance_coil_pair_kernel(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes,
        double* out_M, double* out_dM);
}

namespace coilgun::simulation::cuda {

namespace {

double filament_axial_sep(const components::Armature& arm, int a_idx, int b_idx) {
    int nr = arm.radial_filaments();
    int ia = a_idx / nr + 1;
    int ib = b_idx / nr + 1;
    return std::abs(arm.filament_axial_position(ia) - arm.filament_axial_position(ib));
}

} // anonymous namespace

template<typename SP>
GpuSingleStageSim<SP>::GpuSingleStageSim(
        components::DrivingCoil coil,
        components::Armature    armature,
        std::unique_ptr<Excitation> excitation,
        double                     dt,
        bool                       enable_thermal,
        GpuOptLevel                opt_level,
        const GpuBackend&          backend)
    : coil_(std::move(coil)), armature_(std::move(armature)),
      excitation_(std::move(excitation)),
      dt_(dt), enable_thermal_(enable_thermal),
      opt_level_(opt_level), backend_(backend)
{
    cudaSetDevice(backend_.device_id);

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

    adaptor_.setup({coil_}, armature_, 9);
    physics::upload_gl_nodes(9);

    init_persistent_mode();

    M1_mat_.resize(1, N);
    dM1_mat_.resize(1, N);
    L_total_.resize(N + 1, N + 1);
    RHS_.resize(N + 1);

    state_.currents.resize(N + 1);
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    if (enable_thermal_) {
        state_.filament_temperatures.resize(N);
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    }
}

template<typename SP>
GpuSingleStageSim<SP>::~GpuSingleStageSim() {
    if (use_persistent_) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers(pbuf_);
    }
}

template<typename SP>
void GpuSingleStageSim<SP>::init_persistent_mode() {
    if (backend_.use_persistent) {
        int N_fil = armature_.total_filaments();
        use_persistent_ = init_persistent_buffers(pbuf_, N_fil, backend_);
        if (use_persistent_) {
            launch_persistent_kernel(pbuf_, adaptor_, N_fil,
                                      backend_.threads_per_block, 9, opt_level_);
        }
    }
}

template<typename SP>
void GpuSingleStageSim<SP>::build_filament_M_matrix() {
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
void GpuSingleStageSim<SP>::compute_M1_dM1_fallback() {
    int N = N_fil_;
    int nr = armature_.radial_filaments();
    double arm_center = state_.arm_position;
    double arm_center_init = armature_.position();

    int n_nodes = 9;

    double rai = coil_.inner_radius();
    double rae = coil_.outer_radius();
    double la  = coil_.length();
    int    na  = coil_.turns();

    // Zero device results buffer before kernel writes (prevents stale entries)
    cudaMemset(adaptor_.d_results_M(),  0, N * sizeof(double));
    cudaMemset(adaptor_.d_results_dM(), 0, N * sizeof(double));

    for (int k = 0; k < N; ++k) {
        int i = k / nr + 1;
        int j = k % nr + 1;
        double z_rel = armature_.filament_axial_position(i) - arm_center_init;
        double z_global = arm_center + z_rel;
        double sep = z_global - coil_.position();

        double fil_ri = armature_.filament_inner_radius(j);
        double fil_re = armature_.filament_outer_radius(j);
        double fil_l  = armature_.length() / armature_.axial_filaments();

        physics::mutual_inductance_coil_pair_kernel<<<1, backend_.threads_per_block>>>(
            rai, rae, la, na, fil_ri, fil_re, fil_l, 1, sep, n_nodes,
            adaptor_.d_results_M() + k,
            adaptor_.d_results_dM() + k);
    }
    cudaDeviceSynchronize();

    std::vector<double> h_M(N), h_dM(N);
    adaptor_.download_results(h_M, h_dM, N);

    for (int k = 0; k < N; ++k) {
        M1_mat_(0, k)  = h_M[k];
        dM1_mat_(0, k) = h_dM[k];
    }
}

template<typename SP>
void GpuSingleStageSim<SP>::compute_M1_dM1_persistent() {
    int F = armature_.total_filaments(), N_b = F;
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();

    // Submit one batch and wait for its mapped results before solving the ODE.
    // The previous double-buffer fields did not advance newly completed data,
    // which made the solver reuse a stale one-step-old matrix.
    double arm_center = state_.arm_position;
    for (int fi = 0; fi < F; ++fi) {
        int i = fi / nr + 1, j = fi % nr + 1;
        double z_rel = armature_.filament_axial_position(i) - arm_center_init;
        pbuf_.seps[fi] = arm_center + z_rel - coil_.position();
    }
    *pbuf_.active_pairs = F;
    __sync_synchronize();
    for (int i = 0; i < N_b; ++i) pbuf_.doorbell[i] = 1;
    __sync_synchronize();
    ++*pbuf_.batch_id;

    for (int i = 0; i < N_b; ++i) {
        while (pbuf_.doorbell[i] != 0) __sync_synchronize();
    }
    __sync_synchronize();
    M1_mat_.resize(1, F); dM1_mat_.resize(1, F);
    for (int fi = 0; fi < F; ++fi) {
        M1_mat_(0, fi)  = pbuf_.out_M[fi];
        dM1_mat_(0, fi) = pbuf_.out_dM[fi];
    }
}

template<typename SP>
void GpuSingleStageSim<SP>::compute_M1_dM1() {
    if (use_persistent_)
        compute_M1_dM1_persistent();
    else
        compute_M1_dM1_fallback();
}

template<typename SP>
SimState GpuSingleStageSim<SP>::compute_derivatives(const SimState& s) {
    compute_M1_dM1();

    SimState ds;
    ds.currents.resize(N_fil_ + 1);
    ds.arm_position = s.arm_velocity;
    ds.arm_velocity = compute_force(s) / armature_.mass();
    if (enable_thermal_ && s.filament_temperatures.size() > 0) {
        ds.filament_temperatures.resize(N_fil_);
        ds.filament_temperatures.setZero();
    }

    L_total_.setZero();
    L_total_(0, 0) = L_d_;
    for (int k = 0; k < N_fil_; ++k) {
        L_total_(0, k + 1) = M1_mat_(0, k);
        L_total_(k + 1, 0) = M1_mat_(0, k);
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
        motional_emf += dM1_mat_(0, k) * s.currents(k + 1);

    RHS_(0) = U - R_d_ * I_d - v * motional_emf;
    for (int k = 0; k < N_fil_; ++k) {
        RHS_(k + 1) = -R_fil_(k) * s.currents(k + 1) - v * dM1_mat_(0, k) * I_d;
    }

    Eigen::LDLT<Eigen::MatrixXd> solver(L_total_);
    if (solver.info() != Eigen::Success)
        ds.currents = L_total_.colPivHouseholderQr().solve(RHS_);
    else
        ds.currents = solver.solve(RHS_);

    return ds;
}

template<typename SP>
double GpuSingleStageSim<SP>::compute_force(const SimState& s) {
    double F = 0.0, I_d = s.currents(0);
    for (int k = 0; k < N_fil_; ++k)
        F += I_d * s.currents(k + 1) * dM1_mat_(0, k);
    return F;
}

template<typename SP>
void GpuSingleStageSim<SP>::update_temperatures(SimState& s, double dt_sub) {
    auto mat = armature_.material();
    double beta = physics::material_beta(mat);
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
void GpuSingleStageSim<SP>::record_step(double cap_voltage) {
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
bool GpuSingleStageSim<SP>::check_termination(const TerminationPolicy& policy) {
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
void GpuSingleStageSim<SP>::prepare_summary() {
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
const SimStep& GpuSingleStageSim<SP>::step() {
    state_ = stepper_.advance(dt_, state_,
        [this](const SimState& s) { return compute_derivatives(s); });

    excitation_->advance(dt_, state_.currents(0));

    if (enable_thermal_) update_temperatures(state_, dt_);

    record_step(excitation_->voltage());
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const SimResult& GpuSingleStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const SimResult& GpuSingleStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy) && !excitation_->finished())
        step();
    prepare_summary();
    return result_;
}

template<typename SP>
void GpuSingleStageSim<SP>::reset() {
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

template class GpuSingleStageSim<EulerStepper>;
template class GpuSingleStageSim<RK4Stepper>;

} // namespace coilgun::simulation::cuda
