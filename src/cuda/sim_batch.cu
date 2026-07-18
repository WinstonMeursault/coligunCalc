/**
 * @file sim_batch.cu
 * @brief SimBatch — batch GPU-accelerated simulation implementation.
 * @author Winston Meursault
 *
 * Uses the proven single-pair kernel mutual_inductance_coil_pair_kernel
 * in a loop over simulations and (stage, filament) pairs. Performance
 * is equivalent to individual GpuMultiStageSim calls.
 *
 * The batch kernel (grid with sim_id dimension) is deferred due to
 * nvcc struct-layout issues in cross-module CUDA device code linking.
 * See docs/bugs/P0-batch-kernel.md for details.
 */

#include "coilgun/simulation/cuda/sim_batch.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/excitation.hpp"
#include <atomic>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace coilgun { namespace physics {
    __global__ void mutual_inductance_coil_pair_kernel(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes,
        double* out_M, double* out_dM);
}}

namespace coilgun::simulation::cuda {

using coilgun::physics::mutual_inductance_filament;
using coilgun::physics::mutual_inductance_coil;

namespace {
double filament_axial_sep(const components::Armature& arm, int a, int b) {
    int nr = arm.radial_filaments();
    int ia = a / nr + 1, ib = b / nr + 1;
    return std::abs(arm.filament_axial_position(ia) - arm.filament_axial_position(ib));
}
}

template<typename SP>
SimBatch<SP>::SimBatch(
        std::vector<components::DrivingCoil> coils,
        components::Armature                   armature,
        int                                    num_sims,
        double                                 dt,
        const GpuBackend&                      backend)
    : n_stages_(static_cast<int>(coils.size()))
    , N_fil_(armature.total_filaments())
    , num_sims_(num_sims)
    , dt_(dt), backend_(backend)
    , coils_(std::move(coils)), armature_(std::move(armature))
{
    cudaSetDevice(backend_.device_id);
    if (n_stages_ <= 0)
        throw std::invalid_argument("SimBatch: at least one coil stage is required");
    if (num_sims_ <= 0)
        throw std::invalid_argument("SimBatch: num_sims must be positive");
    if (static_cast<size_t>(num_sims_) > backend_.max_batch_sims)
        throw std::invalid_argument("SimBatch: num_sims exceeds backend.max_batch_sims");
    if (n_stages_ > kMaxStages)
        throw std::invalid_argument("SimBatch: n_stages exceeds kMaxStages");

    R_diag_.resize(n_stages_); L_diag_.resize(n_stages_);
    for (int i = 0; i < n_stages_; ++i) {
        R_diag_(i) = coils_[i].resistance();
        L_diag_(i) = coils_[i].self_inductance();
    }

    int F = N_fil_;
    L_fil_.resize(F); mass_fil_.resize(F);
    const auto& L_arm = armature_.inductances();
    const auto& M_arm = armature_.masses();
    for (int k = 0; k < F; ++k) { L_fil_(k) = L_arm[k]; mass_fil_(k) = M_arm[k]; }

    M_mat_.resize(F, F); M_mat_.setZero();
    int nr = armature_.radial_filaments();
    for (int a = 0; a < F; ++a) {
        int ja = a % nr + 1; double ra = armature_.filament_mean_radius(ja);
        for (int b = a + 1; b < F; ++b) {
            int jb = b % nr + 1; double rb = armature_.filament_mean_radius(jb);
            double s = filament_axial_sep(armature_, a, b);
            double m = mutual_inductance_filament(ra, rb, s, true);
            M_mat_(a, b) = m; M_mat_(b, a) = m;
        }
    }

    M_cc_.resize(n_stages_, n_stages_); M_cc_.setZero();
    for (int i = 0; i < n_stages_; ++i)
        for (int j = i + 1; j < n_stages_; ++j) {
            const auto& a = coils_[i], &b = coils_[j];
            double s = std::abs(a.position() - b.position());
            double m = mutual_inductance_coil(a.inner_radius(),a.outer_radius(),a.length(),a.turns(),
                b.inner_radius(),b.outer_radius(),b.length(),b.turns(),s,9,true);
            M_cc_(i,j)=m; M_cc_(j,i)=m;
        }

    adaptor_.setup(coils_, armature_, 9);

    sims_.resize(num_sims_);
    for (int s = 0; s < num_sims_; ++s) {
        auto& sim = sims_[s];
        sim.triggered.resize(n_stages_, false); sim.finished.resize(n_stages_, false);
        sim.trigger_times.resize(n_stages_, 0.0); sim.triggered[0] = true;
        const auto& R_arm = armature_.resistances();
        sim.R_fil_ref.resize(F); sim.R_fil.resize(F);
        for (int k = 0; k < F; ++k) { sim.R_fil_ref(k)=R_arm[k]; sim.R_fil(k)=R_arm[k]; }
        int dim = n_stages_ + F;
        sim.state.currents.resize(dim); sim.state.currents.setZero();
        sim.state.arm_position = armature_.position();
        sim.state.arm_velocity = armature_.velocity();
    }

    int pc = n_stages_ * F;
    batch_M1_.resize(num_sims_, pc); batch_M1_.setZero();
    batch_dM1_.resize(num_sims_, pc); batch_dM1_.setZero();
    L_total_.resize(n_stages_+F, n_stages_+F);
    RHS_.resize(n_stages_+F);

    init_persistent_mode();
}

template<typename SP>
SimBatch<SP>::~SimBatch() {
    if (use_persistent_) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers(pbuf_);
    }
}

template<typename SP>
void SimBatch<SP>::set_excitations(int sim_id,
    std::vector<std::unique_ptr<Excitation>> excitations,
    std::vector<TriggerConfig>               trigger_configs)
{
    if (sim_id<0||sim_id>=num_sims_) throw std::out_of_range("sim_id out of range");
    if ((int)excitations.size()!=n_stages_) throw std::invalid_argument("excitations size");
    if ((int)trigger_configs.size()!=n_stages_-1) throw std::invalid_argument("triggers size");
    sims_[sim_id].excitations = std::move(excitations);
    sims_[sim_id].trigger_configs = std::move(trigger_configs);
    sims_[sim_id].configured = true;
}

template<typename SP>
void SimBatch<SP>::compute_all_M1_dM1_fallback() {
    int S = n_stages_, F = N_fil_;
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();

    for (int s = 0; s < num_sims_; ++s) {
        auto& sim = sims_[s];
        sim.all_finished = true;
        for (int i=0;i<S;++i) if(sim.triggered[i]&&!sim.finished[i]){sim.all_finished=false;break;}
        if (sim.all_finished) continue;

        double arm_center = sim.state.arm_position;
        for (int si = 0; si < S; ++si) {
            if (!sim.triggered[si] || sim.finished[si]) continue;
            const auto& coil = coils_[si];
            for (int fi = 0; fi < F; ++fi) {
                int i = fi/nr+1, j = fi%nr+1;
                double z_rel = armature_.filament_axial_position(i) - arm_center_init;
                double sep = arm_center + z_rel - coil.position();
                double fil_ri = armature_.filament_inner_radius(j);
                double fil_re = armature_.filament_outer_radius(j);
                double fil_l  = armature_.length()/armature_.axial_filaments();
                int idx = si*F+fi;
                coilgun::physics::mutual_inductance_coil_pair_kernel<<<1,backend_.threads_per_block>>>(
                    coil.inner_radius(),coil.outer_radius(),coil.length(),coil.turns(),
                    fil_ri,fil_re,fil_l,1,sep,9,
                    adaptor_.d_results_M()+idx, adaptor_.d_results_dM()+idx);
            }
        }
        cudaDeviceSynchronize();
        int pc = S*F;
        std::vector<double> hM(pc),hdM(pc);
        cudaMemcpy(hM.data(),adaptor_.d_results_M(),pc*sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(hdM.data(),adaptor_.d_results_dM(),pc*sizeof(double),cudaMemcpyDeviceToHost);
        for(int si=0;si<S;++si) for(int fi=0;fi<F;++fi) {
            batch_M1_(s,si*F+fi)=hM[si*F+fi]; batch_dM1_(s,si*F+fi)=hdM[si*F+fi];
        }
    }
}

template<typename SP>
void SimBatch<SP>::init_persistent_mode() {
    if (backend_.use_persistent) {
        int N_b = n_stages_ * N_fil_;
        use_persistent_ = init_persistent_buffers(pbuf_, N_b, backend_);
        if (use_persistent_) {
            launch_persistent_kernel(pbuf_, adaptor_, N_b,
                                          backend_.threads_per_block, 9, GpuOptLevel::Full);
        }
    }
}

template<typename SP>
void SimBatch<SP>::compute_all_M1_dM1_persistent() {
    int S = n_stages_, F = N_fil_;
    int N_b = S * F;
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();

    for (int s = 0; s < num_sims_; ++s) {
        auto& sim = sims_[s];
        sim.all_finished = true;
        for (int i = 0; i < S; ++i)
            if (sim.triggered[i] && !sim.finished[i]) {
                sim.all_finished = false; break;
            }
        if (sim.all_finished) continue;

        double arm_center = sim.state.arm_position;
        for (int si = 0; si < S; ++si) {
            for (int fi = 0; fi < F; ++fi) {
                int i = fi / nr + 1, j = fi % nr + 1;
                double z_rel = armature_.filament_axial_position(i) - arm_center_init;
                int idx = si * F + fi;
                pbuf_.seps[idx] = arm_center + z_rel - coils_[si].position();
            }
        }

        // Persistent blocks are indexed as si * F + fi.  Process the full
        // fixed index space so inactive pairs cannot shift later geometry.
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

        for (int si = 0; si < S; ++si)
            for (int fi = 0; fi < F; ++fi) {
                int idx = si * F + fi;
                batch_M1_(s, idx)  = pbuf_.out_M[idx];
                batch_dM1_(s, idx) = pbuf_.out_dM[idx];
            }
    }
}

template<typename SP>
void SimBatch<SP>::compute_all_M1_dM1() {
    if (use_persistent_)
        compute_all_M1_dM1_persistent();
    else
        compute_all_M1_dM1_fallback();
}

template<typename SP>
void SimBatch<SP>::build_system_matrix(const SimInstance& sim, Eigen::MatrixXd& L, const Eigen::MatrixXd& M1) const {
    int S=n_stages_, F=N_fil_; L.setZero();
    for(int i=0;i<S;++i){if(!sim.triggered[i]||sim.finished[i]){L(i,i)=1.0;continue;}L(i,i)=L_diag_(i);
        for(int j=0;j<S;++j)if(i!=j&&sim.triggered[j]&&!sim.finished[j])L(i,j)=M_cc_(i,j);}
    for(int i=0;i<S;++i){if(!sim.triggered[i]||sim.finished[i])continue;
        for(int k=0;k<F;++k){L(i,S+k)=M1(i,k);L(S+k,i)=M1(i,k);}}
    for(int a=0;a<F;++a){L(S+a,S+a)=L_fil_(a);
        for(int b=0;b<F;++b)if(a!=b)L(S+a,S+b)=M_mat_(a,b);}
}

template<typename SP>
double SimBatch<SP>::compute_force(const SimInstance& sim, const Eigen::MatrixXd& dM) const {
    int S=n_stages_,F=N_fil_; double F_net=0.0;
    for(int i=0;i<S;++i){if(!sim.triggered[i]||sim.finished[i])continue;
        double Id=sim.state.currents(i);
        for(int k=0;k<F;++k)F_net+=Id*sim.state.currents(S+k)*dM(i,k);}
    return F_net;
}

template<typename SP>
void SimBatch<SP>::solve_and_update(SimInstance& sim) {
    int S=n_stages_,F=N_fil_,dim=S+F;
    int sid=(int)(&sim-sims_.data());
    Eigen::MatrixXd M1(S,F),dM(S,F);
    for(int si=0;si<S;++si)for(int fi=0;fi<F;++fi){
        M1(si,fi)=batch_M1_(sid,si*F+fi); dM(si,fi)=batch_dM1_(sid,si*F+fi);}
    build_system_matrix(sim,L_total_,M1);
    double v=sim.state.arm_velocity; RHS_.setZero();
    for(int i=0;i<S;++i){if(!sim.triggered[i]||sim.finished[i])continue;
        double U=sim.excitations[i]->voltage(),Id=sim.state.currents(i),emf=0.0;
        for(int k=0;k<F;++k)emf+=dM(i,k)*sim.state.currents(S+k);
        RHS_(i)=U-R_diag_(i)*Id-v*emf;}
    for(int k=0;k<F;++k){double If=sim.state.currents(S+k),back=0.0;
        for(int i=0;i<S;++i)if(sim.triggered[i]&&!sim.finished[i])back+=dM(i,k)*sim.state.currents(i);
        RHS_(S+k)=-sim.R_fil(k)*If-v*back;}
    MultiStageState ds; ds.currents.resize(dim);
    Eigen::LDLT<Eigen::MatrixXd> solver(L_total_);
    if (solver.info() != Eigen::Success)
        ds.currents = L_total_.colPivHouseholderQr().solve(RHS_);
    else
        ds.currents = solver.solve(RHS_);
    ds.arm_position=sim.state.arm_velocity;
    ds.arm_velocity=compute_force(sim,dM)/armature_.mass();
    sim.state=stepper_.advance(dt_,sim.state,[&ds](const MultiStageState&){return ds;});
    for(int i=0;i<S;++i)if(sim.triggered[i]&&!sim.finished[i]){
        sim.excitations[i]->advance(dt_,sim.state.currents(i));
        if(sim.excitations[i]->finished())sim.finished[i]=true;}
}

template<typename SP> void SimBatch<SP>::check_triggers(SimInstance& sim) {
    if(n_stages_<=1)return; double t=sim.step_count*dt_;
    for(int i=1;i<n_stages_;++i){if(sim.triggered[i])continue; bool fire=false;
        fire=(sim.trigger_configs[i-1].mode==TriggerMode::Position)
            ?(sim.state.arm_position>=sim.trigger_configs[i-1].value)
            :(t>=sim.trigger_times[i-1]+sim.trigger_configs[i-1].value);
        if(fire){sim.triggered[i]=true;sim.trigger_times[i]=t;}}
}

template<typename SP> void SimBatch<SP>::extinguish_quiet_stages(SimInstance& sim) {
    for(int i=0;i<n_stages_;++i){if(!sim.triggered[i]||sim.finished[i])continue;
        if(sim.excitations[i]->voltage()==0.0&&std::abs(sim.state.currents(i))<1e-6)
            sim.finished[i]=true;}
}

template<typename SP> void SimBatch<SP>::record_step(SimInstance& sim) {
    int S=n_stages_,F=N_fil_,sid=(int)(&sim-sims_.data());
    Eigen::MatrixXd dM(S,F);
    for(int si=0;si<S;++si)for(int fi=0;fi<F;++fi)dM(si,fi)=batch_dM1_(sid,si*F+fi);
    MultiStageStep e; e.state.time=sim.step_count*dt_;
    e.state.arm_position=sim.state.arm_position; e.state.arm_velocity=sim.state.arm_velocity;
    e.state.force=compute_force(sim,dM); e.state.filament_currents.resize(F);
    for(int k=0;k<F;++k)e.state.filament_currents[k]=sim.state.currents(S+k);
    e.cap_voltages.resize(S);e.coil_currents.resize(S);
    for(int i=0;i<S;++i){e.cap_voltages[i]=sim.triggered[i]?sim.excitations[i]->voltage():0.0;
        e.coil_currents[i]=sim.triggered[i]?sim.state.currents(i):0.0;}
    sim.result.history.push_back(std::move(e)); ++sim.step_count;
}

template<typename SP>
bool SimBatch<SP>::check_termination(SimInstance& sim, const TerminationPolicy& policy) const {
    if(policy.enable_bound_check&&sim.state.arm_position>=policy.barrel_end_position)return true;
    if(sim.step_count>=policy.max_steps)return true;
    bool all_finished = true;
    for(int i=0;i<n_stages_;++i)
        if(sim.triggered[i]&&!sim.finished[i]){all_finished=false;break;}
    if (all_finished) return true;

    if (policy.enable_velocity_check &&
        sim.step_count >= policy.velocity_decay_steps &&
        static_cast<int>(sim.result.history.size()) >= policy.velocity_decay_steps + 1) {
        bool decaying = true;
        const auto& history = sim.result.history;
        const int n = static_cast<int>(history.size());
        for (int i = 0; i < policy.velocity_decay_steps; ++i) {
            if (history[n - 1 - i].state.arm_velocity >=
                history[n - 2 - i].state.arm_velocity) {
                decaying = false;
                break;
            }
        }
        if (decaying) {
            Eigen::MatrixXd dM(n_stages_, N_fil_);
            const int sid = static_cast<int>(&sim - sims_.data());
            for (int stage = 0; stage < n_stages_; ++stage)
                for (int filament = 0; filament < N_fil_; ++filament)
                    dM(stage, filament) = batch_dM1_(sid, stage * N_fil_ + filament);
            const double acceleration = compute_force(sim, dM) / armature_.mass();
            if (std::abs(acceleration) < policy.accel_threshold) return true;
        }
    }
    return false;
}

template<typename SP> void SimBatch<SP>::prepare_summary(SimInstance& sim) {
    auto& s=sim.result.summary; s.step_count=sim.step_count;
    if(sim.result.history.empty())return;
    const auto& last=sim.result.history.back();
    s.total_time=last.state.time; s.muzzle_velocity=last.state.arm_velocity;
    for(const auto& h:sim.result.history){
        if(h.state.force>s.max_force)s.max_force=h.state.force;
        for(int i=0;i<n_stages_;++i)if(h.coil_currents[i]>s.peak_coil_current)s.peak_coil_current=h.coil_currents[i];}
    double Et=0.0;
    for(int i=0;i<n_stages_;++i){auto*cap=dynamic_cast<CapacitorExcitation*>(sim.excitations[i].get());
        if(cap)Et+=0.5*cap->capacitance()*cap->initial_voltage()*cap->initial_voltage();}
    double Ek=0.5*armature_.mass()*s.muzzle_velocity*s.muzzle_velocity;
    s.efficiency=(Et>0.0)?Ek/Et:0.0;
}

template<typename SP> void SimBatch<SP>::step_all() {
    for(int s=0;s<num_sims_;++s){check_triggers(sims_[s]);extinguish_quiet_stages(sims_[s]);}
    compute_all_M1_dM1();
    for(int s=0;s<num_sims_;++s){auto& sim=sims_[s];if(sim.all_finished)continue;solve_and_update(sim);record_step(sim);}
}

template<typename SP> void SimBatch<SP>::run(){run(TerminationPolicy::defaults());}
template<typename SP> void SimBatch<SP>::run(const TerminationPolicy& policy) {
    for (const auto& sim : sims_)
        if (!sim.configured)
            throw std::logic_error("SimBatch: set_excitations must configure every simulation before run");
    while(true){bool any=false;for(int s=0;s<num_sims_;++s)if(!check_termination(sims_[s],policy)){any=true;break;}if(!any)break;step_all();}
    for(int s=0;s<num_sims_;++s)prepare_summary(sims_[s]);
}
template<typename SP> const MultiStageResult& SimBatch<SP>::result(int sid) const {return sims_[sid].result;}

template class SimBatch<EulerStepper>;
template class SimBatch<RK4Stepper>;

} // namespace coilgun::simulation::cuda
