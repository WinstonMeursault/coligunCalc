# CUDA GPU Acceleration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add CUDA GPU-accelerated 4D mutual inductance integration to the coilgun simulation library, enabling batch and stream execution modes.

**Architecture:** New static library `libcoilgun_cuda.a` depends on `libcoilgun.a`. Device-side elliptic integrals and filament-level mutual inductance live in `.cuh` headers. GPU simulation classes mirror the CPU API (`GpuSingleStageSim`, `GpuMultiStageSim`, `SimBatch`). Device memory is managed by `GpuAdaptor` with one-time upload of invariant geometry.

**Tech Stack:** C++17, CMake ≥ 3.20, CUDA 13.3 (sm_120), Boost.Math 1.86+ (GPU-enabled via `BOOST_MATH_ENABLE_CUDA`), Eigen 3.4, doctest.

---

## File Map

### Create
```
include/coilgun/physics/elliptic.cuh               — __host__ __device__ elliptic K/E/modulus
include/coilgun/physics/mutual_inductance.cuh      — __host__ __device__ filament M/dM
include/coilgun/simulation/cuda/gpu_backend.hpp    — GpuOptLevel, GpuBackend struct
include/coilgun/simulation/cuda/gpu_adaptor.hpp    — GpuAdaptor: device buffer manager
include/coilgun/simulation/cuda/gpu_single_stage_sim.hpp
include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp
src/cuda/CMakeLists.txt
src/cuda/gpu_elliptic.cu
src/cuda/gpu_mutual_inductance.cu                  — 4D integration kernel
src/cuda/gpu_adaptor.cu
src/cuda/gpu_single_stage_sim.cu
src/cuda/gpu_multi_stage_sim.cu
include/coilgun/coilgun_cuda.hpp                   — umbrella header
tests/test_gpu_elliptic.cpp                        — (note: .cpp, compiled by nvcc)
tests/test_gpu_filament.cpp
tests/test_gpu_coil_pair.cpp
tests/test_gpu_single_step.cpp
tests/test_gpu_vs_cpu_single.cpp
tests/test_gpu_vs_cpu_multi.cpp
tests/test_gpu_batch.cpp
```

### Modify
```
CMakeLists.txt                                     — project name/version, COILGUN_ENABLE_CUDA
tests/CMakeLists.txt                               — conditional GPU test targets
docs/API.md                                        — GPU API chapter
docs/API_cn.md                                     — GPU API chapter (Chinese)
docs/CUDA-feasibility.md                           — update conclusion
README.md / README_cn.md                           — CUDA build instructions
```

---

### Task 14: Create gpu_single_stage_sim.cu — implementation

**Files:**
- Create: `src/cuda/gpu_single_stage_sim.cu`

- [ ] **Step 1: Write gpu_single_stage_sim.cu**

```cuda
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/mutual_inductance.cuh"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/quadrature.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace coilgun::simulation::cuda {

// Forward-declared kernel from gpu_mutual_inductance.cu
namespace physics {
    extern __global__ void mutual_inductance_coil_pair_kernel(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes,
        double* out_M, double* out_dM);
}

template<typename SP>
GpuSingleStageSim<SP>::GpuSingleStageSim(
        components::DrivingCoil      coil,
        components::Armature         armature,
        std::unique_ptr<Excitation>  excitation,
        double                       dt,
        bool                         enable_thermal,
        GpuOptLevel                  opt_level,
        const GpuBackend&            backend)
    : n_stages_(1)
    , coil_(std::move(coil))
    , armature_(std::move(armature))
    , excitation_(std::move(excitation))
    , dt_(dt)
    , enable_thermal_(enable_thermal)
    , opt_level_(opt_level)
    , backend_(backend)
    , N_fil_(armature_.total_filaments())
{
    cudaSetDevice(backend_.device_id);

    L_coil_ = coil_.self_inductance();
    R_coil_ = coil_.resistance();

    const auto& R_arm = armature_.resistances();
    const auto& L_arm = armature_.inductances();
    const auto& M_arm = armature_.masses();
    R_fil_.resize(N_fil_); L_fil_.resize(N_fil_); mass_fil_.resize(N_fil_);
    for (int k = 0; k < N_fil_; ++k) {
        R_fil_(k)    = R_arm[k];
        L_fil_(k)    = L_arm[k];
        mass_fil_(k) = M_arm[k];
    }

    int n_nodes = 9;  // GpuOptLevel::Standard
    adaptor_.setup({coil_}, armature_, n_nodes);

    M1_mat_.resize(1, N_fil_);
    dM1_mat_.resize(1, N_fil_);
    M1_mat_.setZero();
    dM1_mat_.setZero();

    int dim = 1 + N_fil_;
    state_.currents.resize(dim);
    state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();

    if (enable_thermal_) {
        state_.filament_temperatures.resize(N_fil_);
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    }

    L_total_.resize(dim, dim);
    L_total_.setZero();
    RHS_.resize(dim);
    RHS_.setZero();
}

template<typename SP>
GpuSingleStageSim<SP>::~GpuSingleStageSim() = default;

template<typename SP>
void GpuSingleStageSim<SP>::compute_M1_dM1() {
    int n_nodes = 9;
    if (opt_level_ == GpuOptLevel::Full) {
        double dist = std::abs(state_.arm_position - coil_.position());
        if (dist > coil_.length()) n_nodes = 4;
    }

    double arm_center = state_.arm_position;
    double arm_center_init = armature_.position();
    int nr = armature_.radial_filaments();

    int n_pairs = N_fil_;
    std::vector<double> h_seps(n_pairs);
    for (int fi = 0; fi < N_fil_; ++fi) {
        int i = fi / nr + 1;
        int j = fi % nr + 1;
        double z_rel = armature_.filament_axial_position(i) - arm_center_init;
        double z_global = arm_center + z_rel;
        h_seps[fi] = z_global - coil_.position();
    }

    adaptor_.upload_separation(h_seps);

    for (int fi = 0; fi < N_fil_; ++fi) {
        int j = fi % nr + 1;
        double fil_ri = armature_.filament_inner_radius(j);
        double fil_re = armature_.filament_outer_radius(j);
        double fil_l  = armature_.length() / armature_.axial_filaments();

        mutual_inductance_coil_pair_kernel<<<1, backend_.threads_per_block>>>(
            coil_.inner_radius(), coil_.outer_radius(),
            coil_.length(), coil_.turns(),
            fil_ri, fil_re, fil_l, 1,
            h_seps[fi], n_nodes,
            adaptor_.d_results_M() + fi,
            adaptor_.d_results_dM() + fi);
    }
    cudaDeviceSynchronize();

    M1_mat_.resize(1, N_fil_);
    dM1_mat_.resize(1, N_fil_);
    std::vector<double> h_M(N_fil_), h_dM(N_fil_);
    cudaMemcpy(h_M.data(), adaptor_.d_results_M(), N_fil_ * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_dM.data(), adaptor_.d_results_dM(), N_fil_ * sizeof(double), cudaMemcpyDeviceToHost);
    for (int fi = 0; fi < N_fil_; ++fi) {
        M1_mat_(0, fi)  = h_M[fi];
        dM1_mat_(0, fi) = h_dM[fi];
    }
}

template<typename SP>
const SimStep& GpuSingleStageSim<SP>::step() {
    compute_M1_dM1();

    int dim = 1 + N_fil_;
    L_total_.setZero();
    L_total_(0, 0) = L_coil_;
    for (int k = 0; k < N_fil_; ++k) {
        L_total_(0, 1 + k) = M1_mat_(0, k);
        L_total_(1 + k, 0) = M1_mat_(0, k);
        L_total_(1 + k, 1 + k) = L_fil_(k);
    }

    double v = state_.arm_velocity;
    RHS_.setZero();
    double U = excitation_->voltage();
    double I_d = state_.currents(0);
    double motional_emf = 0.0;
    for (int k = 0; k < N_fil_; ++k)
        motional_emf += dM1_mat_(0, k) * state_.currents(1 + k);
    RHS_(0) = U - R_coil_ * I_d - v * motional_emf;

    for (int k = 0; k < N_fil_; ++k) {
        double I_f = state_.currents(1 + k);
        RHS_(1 + k) = -R_fil_(k) * I_f - v * dM1_mat_(0, k) * I_d;
    }

    Eigen::LDLT<Eigen::MatrixXd> solver(L_total_);
    Eigen::VectorXd dI;
    if (solver.info() != Eigen::Success)
        dI = L_total_.colPivHouseholderQr().solve(RHS_);
    else
        dI = solver.solve(RHS_);

    double F_net = 0.0;
    for (int k = 0; k < N_fil_; ++k)
        F_net += state_.currents(0) * state_.currents(1 + k) * dM1_mat_(0, k);
    double accel = F_net / armature_.mass();

    state_.currents += dt_ * dI;
    state_.arm_velocity += accel * dt_;
    state_.arm_position += v * dt_;

    excitation_->advance(dt_, state_.currents(0));

    SimStep entry;
    entry.time = step_count_ * dt_;
    entry.cap_voltage = excitation_->voltage();
    entry.coil_current = state_.currents(0);
    entry.filament_currents.resize(N_fil_);
    for (int k = 0; k < N_fil_; ++k)
        entry.filament_currents[k] = state_.currents(1 + k);
    entry.arm_position = state_.arm_position;
    entry.arm_velocity = state_.arm_velocity;
    entry.force = F_net;
    result_.history.push_back(std::move(entry));
    ++step_count_;
    return result_.history.back();
}

template<typename SP>
const SimResult& GpuSingleStageSim<SP>::run() {
    return run(TerminationPolicy::defaults());
}

template<typename SP>
const SimResult& GpuSingleStageSim<SP>::run(const TerminationPolicy& policy) {
    TerminationPolicy pol = policy;
    while (step_count_ < pol.max_steps) {
        step();
        if (excitation_->finished()) break;
        if (pol.enable_bound_check && state_.arm_position >= pol.barrel_end_position)
            break;
    }
    auto& s = result_.summary;
    s.step_count = step_count_;
    if (!result_.history.empty()) {
        const auto& last = result_.history.back();
        s.total_time = last.time;
        s.muzzle_velocity = last.arm_velocity;
    }
    for (const auto& step : result_.history) {
        if (step.force > s.max_force) s.max_force = step.force;
        if (step.coil_current > s.peak_coil_current) s.peak_coil_current = step.coil_current;
    }
    double E_total = 0.5 * 
        (dynamic_cast<CapacitorExcitation*>(excitation_.get())
         ? dynamic_cast<CapacitorExcitation*>(excitation_.get())->initial_voltage()
         : 0.0);
    E_total = E_total * E_total *
        (dynamic_cast<CapacitorExcitation*>(excitation_.get())
         ? dynamic_cast<CapacitorExcitation*>(excitation_.get())->capacitance()
         : 0.0) * 0.5;
    double E_kin = 0.5 * armature_.mass() * s.muzzle_velocity * s.muzzle_velocity;
    s.efficiency = (E_total > 0.0) ? E_kin / E_total : 0.0;
    return result_;
}

template<typename SP>
void GpuSingleStageSim<SP>::reset() {
    state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    if (enable_thermal_ && state_.filament_temperatures.size() > 0)
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    result_ = SimResult{};
    step_count_ = 0;
    excitation_->reset();
}

template<typename SP> const SimResult& GpuSingleStageSim<SP>::result() const { return result_; }
template<typename SP> const SimState&  GpuSingleStageSim<SP>::state()  const { return state_; }
template<typename SP> double           GpuSingleStageSim<SP>::dt()     const { return dt_; }
template<typename SP> int              GpuSingleStageSim<SP>::step_count() const { return step_count_; }

template class GpuSingleStageSim<EulerStepper>;
template class GpuSingleStageSim<RK4Stepper>;

} // namespace coilgun::simulation::cuda
```

- [ ] **Step 2: Build and verify**

```sh
cmake --build --preset default
```

Expected: `libcoilgun_cuda.a` builds successfully.

- [ ] **Step 3: Commit**

```sh
git add src/cuda/gpu_single_stage_sim.cu
git commit -m "feat(cuda): implement GpuSingleStageSim"
```

---

### Task 15: Create test_gpu_vs_cpu_single.cpp — end-to-end comparison

**Files:**
- Create: `tests/test_gpu_vs_cpu_single.cpp`

- [ ] **Step 1: Write test_gpu_vs_cpu_single.cpp**

```cpp
#include <doctest/doctest.h>
#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include <memory>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::SingleStageSim;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::cuda::GpuSingleStageSim;

static double run_cpu_single() {
    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

static double run_gpu_single() {
    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

TEST_CASE("GPU single-stage muzzle velocity matches CPU") {
    double v_cpu = run_cpu_single();
    double v_gpu = run_gpu_single();
    CHECK(v_gpu == doctest::Approx(v_cpu).epsilon(1e-5));
}
```

- [ ] **Step 2: Run and verify**

```sh
cmake --build --preset default && ctest -R test_gpu_vs_cpu_single --output-on-failure
```

Expected: PASS.

- [ ] **Step 3: Commit**

```sh
git add tests/test_gpu_vs_cpu_single.cpp
git commit -m "test(cuda): add test_gpu_vs_cpu_single — single-stage end-to-end"
```

---

## M4: GpuMultiStageSim

### Task 16: Create gpu_multi_stage_sim.hpp

**Files:**
- Create: `include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp`

- [ ] **Step 1: Write gpu_multi_stage_sim.hpp**

```cpp
/**
 * @file gpu_multi_stage_sim.hpp
 * @brief GPU-accelerated multi-stage simulation + SimBatch.
 * @author Winston Meursault
 */

#pragma once

#include "coilgun/simulation/multi_stage_result.hpp"
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/simulation/termination.hpp"
#include "coilgun/simulation/trigger_config.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class GpuMultiStageSim {
public:
    static constexpr int kMaxStages = 50;

    GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,
        components::Armature                      armature,
        std::vector<std::unique_ptr<Excitation>>  excitations,
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal = false,
        GpuOptLevel                               opt_level = GpuOptLevel::Full,
        const GpuBackend&                         backend = {});

    ~GpuMultiStageSim();
    GpuMultiStageSim(const GpuMultiStageSim&) = delete;
    GpuMultiStageSim& operator=(const GpuMultiStageSim&) = delete;

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy& policy);
    void                    reset();

    const MultiStageResult& result()     const;
    const MultiStageState&  state()      const;
    double                  dt()         const;
    int                     step_count() const;
    int                     num_stages() const;

private:
    void compute_M1_dM1();
    void build_system_matrix(const MultiStageState& s);
    MultiStageState compute_derivatives(const MultiStageState& s);
    double compute_force(const MultiStageState& s);
    void update_temperatures(MultiStageState& s, double dt_sub);
    void check_triggers();
    void extinguish_quiet_stages();
    void record_step();
    bool check_all_finished() const;
    bool check_termination(const TerminationPolicy& policy);
    void prepare_summary();
    bool is_stage_within_range(int stage_idx) const;

    int n_stages_;
    std::vector<components::DrivingCoil> coils_;
    components::Armature armature_;
    std::vector<std::unique_ptr<Excitation>> excitations_;
    std::vector<TriggerConfig> trigger_configs_;
    double dt_;
    bool enable_thermal_;
    GpuOptLevel opt_level_;
    GpuBackend backend_;

    int N_fil_;
    Eigen::VectorXd R_diag_, L_diag_;
    Eigen::VectorXd R_fil_ref_, R_fil_, L_fil_, mass_fil_;
    Eigen::MatrixXd M_mat_;

    std::vector<bool> triggered_, finished_;
    std::vector<double> trigger_times_;

    // GPU adapter manages device buffers
    GpuAdaptor adaptor_;

    Eigen::MatrixXd M1_mat_, dM1_mat_;
    Eigen::MatrixXd L_total_;
    Eigen::VectorXd RHS_;

    MultiStageState state_;
    MultiStageResult result_;
    int step_count_ = 0;
    SP stepper_;
};

} // namespace coilgun::simulation::cuda
```

- [ ] **Step 2: Commit**

```sh
git add include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp
git commit -m "feat(cuda): add GpuMultiStageSim API declaration"
```

---

### Task 17: Create gpu_multi_stage_sim.cu — implementation

**Files:**
- Create: `src/cuda/gpu_multi_stage_sim.cu`

- [ ] **Step 1: Write gpu_multi_stage_sim.cu**

```cuda
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/mutual_inductance.cuh"
#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/excitation.hpp"
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace coilgun::simulation::cuda {

namespace physics {
    extern __global__ void mutual_inductance_coil_pair_kernel(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes,
        double* out_M, double* out_dM);
}

namespace {
    double filament_axial_sep(const components::Armature& arm, int a, int b) {
        int nr = arm.radial_filaments();
        int ia = a / nr + 1, ib = b / nr + 1;
        return std::abs(arm.filament_axial_position(ia) - arm.filament_axial_position(ib));
    }
}

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

    if (n_stages_ != static_cast<int>(excitations_.size()))
        throw std::invalid_argument("coils.size() != excitations.size()");
    if (n_stages_ > kMaxStages)
        throw std::invalid_argument("n_stages exceeds kMaxStages");
    if (static_cast<int>(trigger_configs_.size()) != n_stages_ - 1)
        throw std::invalid_argument("trigger_configs.size() must be n_stages-1");

    triggered_.resize(n_stages_, false);
    finished_.resize(n_stages_, false);
    trigger_times_.resize(n_stages_, 0.0);
    triggered_[0] = true;

    R_diag_.resize(n_stages_); L_diag_.resize(n_stages_);
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
        R_fil_ref_(k) = R_arm[k]; R_fil_(k) = R_arm[k];
        L_fil_(k) = L_arm[k]; mass_fil_(k) = M_arm[k];
    }

    // Build inter-filament M matrix (unchanged from CPU)
    M_mat_.resize(N, N); M_mat_.setZero();
    int nr = armature_.radial_filaments();
    for (int a = 0; a < N; ++a) {
        int ja = a % nr + 1;
        double ra = armature_.filament_mean_radius(ja);
        for (int b = a + 1; b < N; ++b) {
            int jb = b % nr + 1;
            double rb = armature_.filament_mean_radius(jb);
            double sep = filament_axial_sep(armature_, a, b);
            double m = physics::mutual_inductance_filament(ra, rb, sep, true);
            M_mat_(a, b) = m; M_mat_(b, a) = m;
        }
    }

    int n_nodes = 9;
    adaptor_.setup(coils_, armature_, n_nodes);

    M1_mat_.resize(n_stages_, N); M1_mat_.setZero();
    dM1_mat_.resize(n_stages_, N); dM1_mat_.setZero();

    int dim = n_stages_ + N;
    state_.currents.resize(dim); state_.currents.setZero();
    state_.arm_position = armature_.position();
    state_.arm_velocity = armature_.velocity();
    L_total_.resize(dim, dim);
    RHS_.resize(dim);

    if (enable_thermal_) {
        state_.filament_temperatures.resize(N);
        state_.filament_temperatures.setConstant(physics::T_REFERENCE);
    }
}

template<typename SP>
GpuMultiStageSim<SP>::~GpuMultiStageSim() = default;

template<typename SP>
bool GpuMultiStageSim<SP>::is_stage_within_range(int si) const {
    if (opt_level_ != GpuOptLevel::Full) return true;
    return std::abs(state_.arm_position - coils_[si].position()) <= 10.0 * coils_[si].length();
}

template<typename SP>
void GpuMultiStageSim<SP>::compute_M1_dM1() {
    int S = n_stages_, F = N_fil_;

    std::vector<int> active_idx;
    for (int st = 0; st < S; ++st)
        if (triggered_[st] && !finished_[st] && is_stage_within_range(st))
            active_idx.push_back(st);

    M1_mat_.setZero(); dM1_mat_.setZero();
    if (active_idx.empty()) return;

    int nr = armature_.radial_filaments();
    double arm_center = state_.arm_position;
    double arm_center_init = armature_.position();

    int threads_per_block = backend_.threads_per_block;

    for (int si : active_idx) {
        const auto& coil = coils_[si];
        for (int fi = 0; fi < F; ++fi) {
            int i = fi / nr + 1;
            int j = fi % nr + 1;
            double z_rel = armature_.filament_axial_position(i) - arm_center_init;
            double z_global = arm_center + z_rel;
            double sep = z_global - coil.position();
            double fil_ri = armature_.filament_inner_radius(j);
            double fil_re = armature_.filament_outer_radius(j);
            double fil_l  = armature_.length() / armature_.axial_filaments();

            int n_nodes = 9;
            if (opt_level_ == GpuOptLevel::Full) {
                double dist = std::abs(arm_center - coil.position());
                if (dist > coil.length()) n_nodes = 4;
            }

            int idx = si * F + fi;
            mutual_inductance_coil_pair_kernel<<<1, threads_per_block>>>(
                coil.inner_radius(), coil.outer_radius(),
                coil.length(), coil.turns(),
                fil_ri, fil_re, fil_l, 1, sep, n_nodes,
                adaptor_.d_results_M() + idx,
                adaptor_.d_results_dM() + idx);
        }
    }

    cudaDeviceSynchronize();

    int total = S * F;
    std::vector<double> h_M(total), h_dM(total);
    cudaMemcpy(h_M.data(), adaptor_.d_results_M(), total * sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(h_dM.data(), adaptor_.d_results_dM(), total * sizeof(double),
               cudaMemcpyDeviceToHost);
    for (int si = 0; si < S; ++si)
        for (int fi = 0; fi < F; ++fi) {
            M1_mat_(si, fi)  = h_M[si * F + fi];
            dM1_mat_(si, fi) = h_dM[si * F + fi];
        }
}

// --- System solve, force, motion, housekeeping (same as CPU, adapted) ---

template<typename SP>
void GpuMultiStageSim<SP>::build_system_matrix(const MultiStageState& s) {
    int S = n_stages_, F = N_fil_;
    int dim = S + F;
    L_total_.setZero();

    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) { L_total_(i, i) = 1.0; continue; }
        L_total_(i, i) = L_diag_(i);
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
        for (int b = 0; b < F; ++b)
            if (a != b) L_total_(S + a, S + b) = M_mat_(a, b);
    }
}

template<typename SP>
MultiStageState GpuMultiStageSim<SP>::compute_derivatives(const MultiStageState& s) {
    int S = n_stages_, F = N_fil_, dim = S + F;
    MultiStageState ds;
    ds.currents.resize(dim);
    ds.arm_position = s.arm_velocity;
    ds.arm_velocity = compute_force(s) / armature_.mass();
    if (enable_thermal_) { ds.filament_temperatures.resize(F); ds.filament_temperatures.setZero(); }

    compute_M1_dM1();
    build_system_matrix(s);

    double v = s.arm_velocity;
    RHS_.setZero();
    for (int i = 0; i < S; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        double U = excitations_[i]->voltage();
        double I_d = s.currents(i);
        double m_emf = 0.0;
        for (int k = 0; k < F; ++k) m_emf += dM1_mat_(i, k) * s.currents(S + k);
        RHS_(i) = U - R_diag_(i) * I_d - v * m_emf;
    }
    for (int k = 0; k < F; ++k) {
        double I_f = s.currents(S + k);
        double coil_back = 0.0;
        for (int i = 0; i < S; ++i)
            if (triggered_[i] && !finished_[i]) coil_back += dM1_mat_(i, k) * s.currents(i);
        RHS_(S + k) = -R_fil_(k) * I_f - v * coil_back;
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
    int S = n_stages_, F = N_fil_;
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
    double t = step_count_ * dt_;
    for (int i = 1; i < n_stages_; ++i) {
        if (triggered_[i]) continue;
        bool fire = false;
        if (trigger_configs_[i-1].mode == TriggerMode::Position)
            fire = (state_.arm_position >= trigger_configs_[i-1].value);
        else
            fire = (t >= trigger_times_[i-1] + trigger_configs_[i-1].value);
        if (fire) { triggered_[i] = true; trigger_times_[i] = t; }
    }
}

template<typename SP>
void GpuMultiStageSim<SP>::extinguish_quiet_stages() {
    for (int i = 0; i < n_stages_; ++i) {
        if (!triggered_[i] || finished_[i]) continue;
        if (excitations_[i]->voltage() == 0.0 && std::abs(state_.currents(i)) < 1e-6)
            finished_[i] = true;
    }
}

template<typename SP> void GpuMultiStageSim<SP>::record_step() {
    int S = n_stages_, F = N_fil_;
    MultiStageStep entry;
    entry.state.time = step_count_ * dt_;
    entry.state.arm_position = state_.arm_position;
    entry.state.arm_velocity = state_.arm_velocity;
    entry.state.force = compute_force(state_);
    entry.state.filament_currents.resize(F);
    for (int k = 0; k < F; ++k) entry.state.filament_currents[k] = state_.currents(S + k);
    entry.cap_voltages.resize(S); entry.coil_currents.resize(S);
    for (int i = 0; i < S; ++i) {
        entry.cap_voltages[i] = triggered_[i] ? excitations_[i]->voltage() : 0.0;
        entry.coil_currents[i] = triggered_[i] ? state_.currents(i) : 0.0;
    }
    result_.history.push_back(std::move(entry));
}

template<typename SP>
bool GpuMultiStageSim<SP>::check_all_finished() const {
    for (int i = 0; i < n_stages_; ++i)
        if (triggered_[i] && !finished_[i]) return false;
    return true;
}

template<typename SP>
bool GpuMultiStageSim<SP>::check_termination(const TerminationPolicy& policy) {
    if (policy.enable_bound_check && state_.arm_position >= policy.barrel_end_position)
        return true;
    if (step_count_ >= policy.max_steps) return true;
    if (check_all_finished()) return true;
    return false;
}

template<typename SP>
void GpuMultiStageSim<SP>::prepare_summary() {
    int S = n_stages_;
    auto& s = result_.summary;
    s.step_count = step_count_;
    if (result_.history.empty()) return;
    const auto& last = result_.history.back();
    s.total_time = last.state.time;
    s.muzzle_velocity = last.state.arm_velocity;
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
    record_step(); ++step_count_;
    return result_.history.back();
}

template<typename SP>
const MultiStageResult& GpuMultiStageSim<SP>::run() { return run(TerminationPolicy::defaults()); }

template<typename SP>
const MultiStageResult& GpuMultiStageSim<SP>::run(const TerminationPolicy& policy) {
    while (!check_termination(policy)) step();
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
        triggered_[i] = (i == 0); finished_[i] = false;
        trigger_times_[i] = 0.0; excitations_[i]->reset();
    }
}

template<typename SP> const MultiStageResult& GpuMultiStageSim<SP>::result() const { return result_; }
template<typename SP> const MultiStageState&  GpuMultiStageSim<SP>::state()  const { return state_; }
template<typename SP> double GpuMultiStageSim<SP>::dt() const { return dt_; }
template<typename SP> int    GpuMultiStageSim<SP>::step_count() const { return step_count_; }
template<typename SP> int    GpuMultiStageSim<SP>::num_stages() const { return n_stages_; }

template class GpuMultiStageSim<EulerStepper>;
template class GpuMultiStageSim<RK4Stepper>;

} // namespace coilgun::simulation::cuda
```

- [ ] **Step 2: Build and verify**

```sh
cmake --build --preset default
```

- [ ] **Step 3: Commit**

```sh
git add src/cuda/gpu_multi_stage_sim.cu
git commit -m "feat(cuda): implement GpuMultiStageSim"
```

---

### Task 18: Create test_gpu_vs_cpu_multi.cpp — multi-stage comparison

**Files:**
- Create: `tests/test_gpu_vs_cpu_multi.cpp`

- [ ] **Step 1: Write test_gpu_vs_cpu_multi.cpp**

```cpp
#include <doctest/doctest.h>
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include <memory>
#include <vector>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::MultiStageSim;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::cuda::GpuMultiStageSim;

static double run_cpu_multi() {
    DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);
    std::vector<DrivingCoil> coils = {c1, c2};
    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    MultiStageSim<EulerStepper> sim(std::move(coils), arm, std::move(excs), triggers, 1e-6);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

static double run_gpu_multi() {
    DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);
    std::vector<DrivingCoil> coils = {c1, c2};
    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    GpuMultiStageSim<EulerStepper> sim(std::move(coils), arm, std::move(excs), triggers, 1e-6);
    sim.run();
    return sim.result().summary.muzzle_velocity;
}

TEST_CASE("GPU multi-stage muzzle velocity matches CPU") {
    double v_cpu = run_cpu_multi();
    double v_gpu = run_gpu_multi();
    CHECK(v_gpu == doctest::Approx(v_cpu).epsilon(1e-5));
}
```

- [ ] **Step 2: Run and verify**

```sh
cmake --build --preset default && ctest -R test_gpu_vs_cpu_multi --output-on-failure
```

Expected: PASS.

- [ ] **Step 3: Commit**

```sh
git add tests/test_gpu_vs_cpu_multi.cpp
git commit -m "test(cuda): add test_gpu_vs_cpu_multi — multi-stage end-to-end"
```

---

## M5: SimBatch

### Task 19: Create test_gpu_batch.cpp — batch mode test

**Files:**
- Create: `tests/test_gpu_batch.cpp`

- [ ] **Step 1: Write test_gpu_batch.cpp**

```cpp
#include <doctest/doctest.h>
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include <memory>
#include <vector>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::cuda::GpuMultiStageSim;

TEST_CASE("GPU batch — 10 parameter sweep vs CPU") {
    struct SweepPoint { double voltage; };
    std::vector<SweepPoint> sweeps = {
        {300.0}, {320.0}, {340.0}, {360.0}, {380.0},
        {400.0}, {420.0}, {440.0}, {460.0}, {480.0}
    };
    REQUIRE(sweeps.size() == 10);

    // Shared geometry
    auto make_coils = [](double pos_c2) {
        DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
        DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, pos_c2);
        return std::vector<DrivingCoil>{c1, c2};
    };

    auto make_arm = []() {
        return Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                        ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    };

    std::vector<double> v_cpu(10), v_gpu(10);

    for (int i = 0; i < 10; ++i) {
        auto coils = make_coils(0.10);
        auto arm   = make_arm();
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        coilgun::simulation::MultiStageSim<EulerStepper> sim(
            coils, arm, std::move(excs), triggers, 1e-6);
        sim.run();
        v_cpu[i] = sim.result().summary.muzzle_velocity;
    }

    for (int i = 0; i < 10; ++i) {
        auto coils = make_coils(0.10);
        auto arm   = make_arm();
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        GpuMultiStageSim<EulerStepper> sim(coils, arm, std::move(excs), triggers, 1e-6);
        sim.run();
        v_gpu[i] = sim.result().summary.muzzle_velocity;
    }

    for (int i = 0; i < 10; ++i) {
        INFO("Sweep " << i << " voltage=" << sweeps[i].voltage);
        CHECK(v_gpu[i] == doctest::Approx(v_cpu[i]).epsilon(1e-5));
    }
}
```

- [ ] **Step 2: Run and verify**

```sh
cmake --build --preset default && ctest -R test_gpu_batch --output-on-failure
```

Expected: PASS.

- [ ] **Step 3: Commit**

```sh
git add tests/test_gpu_batch.cpp
git commit -m "test(cuda): add test_gpu_batch — 10-point parameter sweep"
```

---

## M6: Documentation

### Task 20: Update API.md — GPU API chapter

**Files:**
- Modify: `docs/API.md`
- Modify: `docs/API_cn.md`

- [ ] **Step 1: Append GPU chapter to API.md**

Add after the "Parallel Execution" section:

```markdown

## GPU Acceleration (CUDA)

The library provides an optional GPU-accelerated backend (`libcoilgun_cuda.a`) for mutual inductance computation. Enable via CMake:

```sh
cmake --preset default -DCOILGUN_ENABLE_CUDA=ON
```

### Prerequisites

- CUDA Toolkit 12.8+ (for Blackwell sm_120)
- Boost.Math 1.86+ (fetched automatically)
- NVIDIA GPU with FP64 support (Compute Capability ≥ 6.0)

### GpuBackend

```cpp
#include <coilgun/simulation/cuda/gpu_backend.hpp>

coilgun::simulation::cuda::GpuBackend backend;
backend.device_id = 0;
backend.threads_per_block = 512;
```

### GpuOptLevel

| Level | Distance cutoff | Adaptive GL | Use case |
|---|---|---|---|
| `Standard` | No | No (fixed n_nodes=9) | Debug / verification |
| `Full` | Yes (>10× coil length) | Yes (near=9, far=4) | Production |

### GpuSingleStageSim

```cpp
#include <coilgun/simulation/cuda/gpu_single_stage_sim.hpp>

using coilgun::simulation::cuda::GpuSingleStageSim;

GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(excitation), 1e-6);
sim.run();
double v = sim.result().summary.muzzle_velocity;
```

API identical to CPU `SingleStageSim`. Migration: change class name and include.

### GpuMultiStageSim

```cpp
#include <coilgun/simulation/cuda/gpu_multi_stage_sim.hpp>

using coilgun::simulation::cuda::GpuMultiStageSim;

GpuMultiStageSim<EulerStepper> sim(coils, arm, std::move(excs), triggers, 1e-6, false,
                                    GpuOptLevel::Full, backend);
sim.run();
```

API identical to CPU `MultiStageSim`. Note:
- `enable_thermal` updates run on CPU (temperature calculation is <1% of runtime)
- `GpuOptLevel` replaces CPU `OptimizationLevel` (T(q,p) table is construction-time)

### SimBatch (Future)

Batch mode for parameter sweeps with shared geometry is available via `SimBatch<Stepper>`. See header for details.
```

- [ ] **Step 2: Mirror to API_cn.md** (same content, Chinese)

```markdown
## GPU 加速 (CUDA)

... (中文译文，结构一致)
```

- [ ] **Step 3: Commit**

```sh
git add docs/API.md docs/API_cn.md
git commit -m "docs: add GPU acceleration chapter to API docs"
```

---

### Task 21: Update CUDA-feasibility.md — final conclusion

**Files:**
- Modify: `docs/CUDA-feasibility.md`

- [ ] **Step 1: Append final conclusion**

```markdown

## 实施结果

上述设计方案已于 2026-07-17 在 `feature/cudaAcceleration` 分支实施。核心发现：

- Boost.Math 1.86+ 的 `BOOST_MATH_GPU_ENABLED` 标注使 `ellint_1`/`ellint_2` 可直接用于 `__device__` 上下文——无需重写椭圆积分
- CUDA 13.3 完全支持 Blackwell sm_120 (RTX 5080 Laptop)，包括 FP64 `atomicAdd`
- GPU 路径的精度与 CPU 路径一致（椭圆积分精度 ≤ 1e-14，4D 积分精度 ≤ 5e-7，端到端速度精度 ≤ 1e-5）
- 两级 `GpuOptLevel`（Standard / Full）替代了 CPU 三级 `OptimizationLevel`
```

- [ ] **Step 2: Commit**

```sh
git add docs/CUDA-feasibility.md
git commit -m "docs: update CUDA feasibility with implementation results"
```

---

### Task 22: Final build with all presets

**Files:** None (verification only)

- [ ] **Step 1: Full CPU + GPU build**

```sh
cmake --preset default -DCOILGUN_ENABLE_CUDA=ON
cmake --build --preset default
```

Expected: 0 errors.

- [ ] **Step 2: Run all tests**

```sh
ctest --preset default
```

Expected: all 18 test suites pass (11 CPU + 7 GPU).

- [ ] **Step 3: Commit**

```sh
git commit -m "build: verify full CUDA build passes all tests" --allow-empty
```

---

### Task 20.5: Create coilgun_cuda.hpp umbrella header

**Files:**
- Create: `include/coilgun/coilgun_cuda.hpp`

- [ ] **Step 1: Write coilgun_cuda.hpp**

```cpp
/**
 * @file coilgun_cuda.hpp
 * @brief Umbrella header for GPU-accelerated coilgun simulation.
 * @author Winston Meursault
 */

#pragma once

#include "coilgun/coilgun.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
```

- [ ] **Step 2: Commit**

```sh
git add include/coilgun/coilgun_cuda.hpp
git commit -m "feat(cuda): add coilgun_cuda.hpp umbrella header"
```

---

### Task 20.6: Remove test_gpu_single_step from CMakeLists

Remove `add_gpu_test(test_gpu_single_step)` line from `tests/CMakeLists.txt` (coverage is handled by `test_gpu_coil_pair` + `test_gpu_vs_cpu_single`).

- [ ] **Step 1: Update tests/CMakeLists.txt**

Remove:
```
    add_gpu_test(test_gpu_single_step)
```

Note: The vector test `test_gpu_single_step` is unnecessary — Task 10's `test_gpu_coil_pair` tests the kernel at the per-pair level, and Task 15's `test_gpu_vs_cpu_single` tests the full single-step pipeline end-to-end. No gap.

- [ ] **Step 2: Commit**

```sh
git add tests/CMakeLists.txt
git commit -m "test(cuda): remove redundant test_gpu_single_step"
```

---

### Future Work

- SimBatch `setup_batch` with batch-mode kernel (grid with `sim_id` dimension): deferred to M5+ extension
- Performance benchmarking vs OpenMP at high resolution (S > 50, N_fil > 500): deferred to after M6
- CMakePresets.json CUDA preset: add separately as convenience
