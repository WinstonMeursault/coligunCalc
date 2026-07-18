/**
 * @file gpu_multi_stage_sim.cu
 * @brief GPU-accelerated multi-stage coilgun simulation — implementation.
 * @author Winston Meursault
 *
 * Uses CUDA kernels for the 4D Gauss-Legendre mutual inductance
 * integration. Linear system solve (Eigen LDLT) and kinematic/thermal
 * updates remain on the CPU identically to the multi_stage_sim.cpp path.
 * Force computation and temperature update are serial (no OpenMP).
 *
 * @see multi_stage_sim.cpp
 */

#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/constants.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <atomic>
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
GpuMultiStageSim<SP>::GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,
        components::Armature                      armature,
        std::vector<std::unique_ptr<Excitation>>  excitations,
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal,
        GpuOptLevel                               opt_level,
        const GpuBackend&                         backend)
    : n_stages_(static_cast<int>(coils.size()))
    , coils_(std::move(coils))
    , armature_(std::move(armature))
    , excitations_(std::move(excitations))
    , trigger_configs_(std::move(trigger_configs))
    , dt_(dt)
    , enable_thermal_(enable_thermal)
    , opt_level_(opt_level)
    , backend_(backend)
    , N_fil_(armature_.total_filaments())
{
    cudaSetDevice(backend_.device_id);

    if (n_stages_ != static_cast<int>(excitations_.size())) {
        throw std::invalid_argument(
            "GpuMultiStageSim: coils.size() != excitations.size()");
    }
    if (n_stages_ > kMaxStages) {
        throw std::invalid_argument(
            "GpuMultiStageSim: n_stages exceeds kMaxStages");
    }
    if (static_cast<int>(trigger_configs_.size()) != n_stages_ - 1) {
        throw std::invalid_argument(
            "GpuMultiStageSim: trigger_configs.size() must be n_stages-1");
    }

    triggered_.resize(n_stages_, false);
    finished_.resize(n_stages_, false);
    trigger_times_.resize(n_stages_, 0.0);
    triggered_[0] = true;

    R_diag_.resize(n_stages_);
    L_diag_.resize(n_stages_);
    for (int i = 0; i < n_stages_; ++i) {
        R_diag_(i) = coils_[i].resistance();
        L_diag_(i) = coils_[i].self_inductance();
    }

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

    adaptor_.setup(coils_, armature_, 9);
    physics::upload_gl_nodes(9);

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

    init_persistent_mode();
}

template<typename SP>
void GpuMultiStageSim<SP>::init_persistent_mode() {
    if (backend_.use_persistent) {
        int N_b = n_stages_ * N_fil_;
        use_persistent_ = init_persistent_buffers(pbuf_, N_b, backend_);
        if (use_persistent_) {
            launch_persistent_kernel(pbuf_, adaptor_, N_b,
                                          backend_.threads_per_block, 9, opt_level_);
        }
    }
}

template<typename SP>
GpuMultiStageSim<SP>::~GpuMultiStageSim() {
    if (use_persistent_) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers(pbuf_);
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::build_filament_M_matrix() {
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
bool GpuMultiStageSim<SP>::is_stage_within_range(int stage_idx) const {
    if (opt_level_ == GpuOptLevel::Standard) return true;
    const auto& coil = coils_[stage_idx];
    double cutoff = 10.0 * coil.length();
    double dist = std::abs(state_.arm_position - coil.position());
    return dist <= cutoff;
}

template<typename SP>
void GpuMultiStageSim<SP>::compute_M1_dM1_fallback() {
    int S = n_stages_;
    int F = N_fil_;
    int nr = armature_.radial_filaments();
    double arm_center = state_.arm_position;
    double arm_center_init = armature_.position();

    M1_mat_.setZero();
    dM1_mat_.setZero();

    int n_launched = 0;
    for (int si = 0; si < S; ++si) {
        if (!triggered_[si] || finished_[si]) continue;
        if (!is_stage_within_range(si)) continue;

        const auto& coil = coils_[si];

        int n_nodes = 9;

        double rai = coil.inner_radius();
        double rae = coil.outer_radius();
        double la  = coil.length();
        int    na  = coil.turns();

        for (int k = 0; k < F; ++k) {
            int i = k / nr + 1;
            int j = k % nr + 1;
            double z_rel = armature_.filament_axial_position(i) - arm_center_init;
            double z_global = arm_center + z_rel;
            double sep = z_global - coil.position();

            double fil_ri = armature_.filament_inner_radius(j);
            double fil_re = armature_.filament_outer_radius(j);
            double fil_l  = armature_.length() / armature_.axial_filaments();

            int flat_idx = si * F + k;
            physics::mutual_inductance_coil_pair_kernel<<<1, backend_.threads_per_block>>>(
                rai, rae, la, na, fil_ri, fil_re, fil_l, 1, sep, n_nodes,
                adaptor_.d_results_M() + flat_idx,
                adaptor_.d_results_dM() + flat_idx);
            ++n_launched;
        }
    }

    if (n_launched == 0) return;

    cudaDeviceSynchronize();

    int total_pairs = S * F;
    std::vector<double> h_M(total_pairs), h_dM(total_pairs);
    adaptor_.download_results(h_M, h_dM, total_pairs);

    for (int si = 0; si < S; ++si) {
        for (int k = 0; k < F; ++k) {
            int idx = si * F + k;
            M1_mat_(si, k)  = h_M[idx];
            dM1_mat_(si, k) = h_dM[idx];
        }
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::compute_M1_dM1_persistent() {
    int S = n_stages_, F = N_fil_;
    int N_b = S * F;
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();
    double arm_center = state_.arm_position;

    // Keep the persistent-kernel index space identical to the stage/filament
    // matrix.  The kernel derives its geometry from idx, so compacting active
    // pairs would associate a separation with the wrong stage or filament.
    for (int si = 0; si < S; ++si) {
        for (int fi = 0; fi < F; ++fi) {
            int i = fi / nr + 1, j = fi % nr + 1;
            double z_rel = armature_.filament_axial_position(i) - arm_center_init;
            int idx = si * F + fi;
            pbuf_.seps[idx] = arm_center + z_rel - coils_[si].position();
        }
    }

    // Inactive pairs are harmless: their matrix entries are ignored when the
    // corresponding stage is inactive, while processing the full fixed index
    // space preserves the kernel's geometry mapping.
    *pbuf_.active_pairs = N_b;
    __sync_synchronize();

    for (int i = 0; i < N_b; ++i)
        pbuf_.doorbell[i] = 1;
    __sync_synchronize();
    ++*pbuf_.batch_id;

    for (int i = 0; i < N_b; ++i) {
        while (pbuf_.doorbell[i] != 0)
            __sync_synchronize();
    }

    M1_mat_.resize(S, F);
    dM1_mat_.resize(S, F);
    for (int si = 0; si < S; ++si)
        for (int fi = 0; fi < F; ++fi) {
            int idx = si * F + fi;
            M1_mat_(si, fi)  = pbuf_.out_M[idx];
            dM1_mat_(si, fi) = pbuf_.out_dM[idx];
        }
}

template<typename SP>
void GpuMultiStageSim<SP>::compute_M1_dM1() {
    if (use_persistent_)
        compute_M1_dM1_persistent();
    else
        compute_M1_dM1_fallback();
}

template<typename SP>
void GpuMultiStageSim<SP>::build_system_matrix(const MultiStageState& s) {
    int S = n_stages_;
    int F = N_fil_;
    int dim = S + F;
    L_total_.setZero();

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

    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        for (int k = 0; k < F; ++k) {
            L_total_(i, S + k) = M1_mat_(i, k);
            L_total_(S + k, i) = M1_mat_(i, k);
        }
    }

    for (int a = 0; a < F; ++a) {
        L_total_(S + a, S + a) = L_fil_(a);
        for (int b = 0; b < F; ++b) {
            if (a != b) L_total_(S + a, S + b) = M_mat_(a, b);
        }
    }
}

template<typename SP>
MultiStageState GpuMultiStageSim<SP>::compute_derivatives(const MultiStageState& s) {
    int S = n_stages_;
    int F = N_fil_;
    int dim = S + F;

    MultiStageState ds;
    ds.currents.resize(dim);
    ds.arm_position = s.arm_velocity;
    ds.arm_velocity = compute_force(s) / armature_.mass();
    if (enable_thermal_ && s.filament_temperatures.size() > 0) {
        ds.filament_temperatures.resize(F);
        ds.filament_temperatures.setZero();
    }

    compute_M1_dM1();

    build_system_matrix(s);

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
double GpuMultiStageSim<SP>::compute_force(const MultiStageState& s) {
    int S = n_stages_;
    int F = N_fil_;
    double F_net = 0.0;
    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        double I_d = s.currents(i);
        for (int k = 0; k < F; ++k)
            F_net += I_d * s.currents(S + k) * dM1_mat_(i, k);
    }
    return F_net;
}

template<typename SP>
void GpuMultiStageSim<SP>::update_temperatures(MultiStageState& s, double dt_sub) {
    int F = N_fil_;
    auto mat = armature_.material();
    double beta = physics::material_beta(mat);
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
void GpuMultiStageSim<SP>::check_triggers() {
    if (n_stages_ <= 1) return;
    double current_time = step_count_ * dt_;
    for (int i = 1; i < n_stages_; ++i) {
        if (triggered_[i]) continue;
        const auto& cfg = trigger_configs_[i - 1];
        bool fire = false;
        if (cfg.mode == TriggerMode::Position) {
            fire = (state_.arm_position >= cfg.value);
        } else {
            fire = (current_time >= trigger_times_[i - 1] + cfg.value);
        }
        if (fire) {
            triggered_[i] = true;
            trigger_times_[i] = current_time;
        }
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::extinguish_quiet_stages() {
    static constexpr double kQuietThreshold = 1e-6;
    for (int i = 0; i < n_stages_; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        if (excitations_[i]->voltage() == 0.0 &&
            std::abs(state_.currents(i)) < kQuietThreshold) {
            finished_[i] = true;
        }
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::record_step() {
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
    for (int i = 0; i < S; ++i) {
        entry.cap_voltages[i]  = triggered_[i] ? excitations_[i]->voltage() : 0.0;
        entry.coil_currents[i] = triggered_[i] ? state_.currents(i) : 0.0;
    }

    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool GpuMultiStageSim<SP>::check_all_finished() const {
    for (int i = 0; i < n_stages_; ++i) {
        if (triggered_[i] && !finished_[i]) return false;
    }
    return true;
}

template<typename SP>
bool GpuMultiStageSim<SP>::check_termination(const TerminationPolicy& policy) {
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
void GpuMultiStageSim<SP>::prepare_summary() {
    int S = n_stages_;
    auto& s = result_.summary;
    s.step_count = step_count_;
    if (result_.history.empty()) return;

    const auto& last = result_.history.back();
    s.total_time      = last.state.time;
    s.muzzle_velocity = last.state.arm_velocity;

    for (int i = 0; i < S; ++i) {
        if (!triggered_[i]) continue;
        PerStageSummary ps;
        ps.stage_index   = i;
        ps.trigger_time  = trigger_times_[i];

        for (const auto& step : result_.history) {
            if (step.state.time >= ps.trigger_time) {
                ps.trigger_position = step.state.arm_position;
                break;
            }
        }

        auto* cap = dynamic_cast<CapacitorExcitation*>(excitations_[i].get());
        bool active = false;
        double E_init = 0.0;

        for (const auto& step : result_.history) {
            if (step.coil_currents[i] > ps.peak_current)
                ps.peak_current = step.coil_currents[i];
            if (step.state.force > ps.max_force)
                ps.max_force = step.state.force;

            if (triggered_[i] && step.coil_currents[i] > 1e-6)
                active = true;
            if (active) ps.step_count_active++;

            if (cap) {
                if (!active && step.coil_currents[i] > 1e-6 && E_init == 0.0)
                    E_init = 0.5 * cap->capacitance() * step.cap_voltages[i] * step.cap_voltages[i];
            }
        }

        if (cap) {
            double E_end = 0.5 * cap->capacitance() * cap->capacitor_voltage() * cap->capacitor_voltage();
            ps.energy_depleted = E_init - E_end;
        }

        s.per_stage.push_back(ps);
    }

    for (const auto& step : result_.history) {
        if (step.state.force > s.max_force) s.max_force = step.state.force;
    }
    for (int i = 0; i < S; ++i) {
        for (const auto& step : result_.history) {
            if (step.coil_currents[i] > s.peak_coil_current)
                s.peak_coil_current = step.coil_currents[i];
        }
    }

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
const MultiStageStep& GpuMultiStageSim<SP>::step() {
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
const MultiStageResult& GpuMultiStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const MultiStageResult& GpuMultiStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy))
        step();
    prepare_summary();
    return result_;
}

template<typename SP>
void GpuMultiStageSim<SP>::reset() {
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
        excitations_[i]->reset();
    }
}

template class GpuMultiStageSim<EulerStepper>;
template class GpuMultiStageSim<RK4Stepper>;

} // namespace coilgun::simulation::cuda
