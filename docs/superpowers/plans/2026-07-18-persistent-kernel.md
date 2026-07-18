# Persistent Kernel — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace per-sim kernel launch loop in `SimBatch::compute_all_M1_dM1` with a single persistent kernel using mapped memory, achieving 0 kernel launches per step.

**Architecture:** Persistent kernel is defined in `sim_batch.cu` in `namespace coilgun::simulation::cuda`. Mapped memory buffers (`cudaHostAllocMapped`) are managed by `PersistentBuffers` struct in `sim_batch.hpp`. Kernel uses doorbell/batch_id protocol for host-device synchronization. `GpuBackend::use_persistent` flag controls persistent vs. fallback mode.

**Tech Stack:** C++17, CUDA 13.3, Eigen 3.4, doctest.

---

## File Map

### Modify
```
include/coilgun/simulation/cuda/sim_batch.hpp    — PersistentBuffers struct, new private methods
include/coilgun/simulation/cuda/gpu_backend.hpp   — use_persistent flag
src/cuda/sim_batch.cu                             — persistent_batch_kernel + integration
tests/test_gpu_sim_batch.cpp                      — persistent mode test
```

---

### Task 1: Add `use_persistent` to GpuBackend

**Files:**
- Modify: `include/coilgun/simulation/cuda/gpu_backend.hpp`

- [ ] **Step 1: Read current file**

Read `include/coilgun/simulation/cuda/gpu_backend.hpp`.

- [ ] **Step 2: Add use_persistent flag to GpuBackend struct**

```cpp
struct GpuBackend {
    int     device_id         = 0;    ///< cudaSetDevice target.
    int     threads_per_block = 512;  ///< Threads per block for integration kernel.
    size_t  max_batch_sims    = 256;  ///< Pre-allocated buffer size for batch mode.
    bool    enable_profiling  = false; ///< Enable NVTX range annotations.
    bool    use_persistent    = true;  ///< Use persistent kernel (mapped memory). Falls back to per-sim launches if false.
};
```

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/gpu_backend.hpp
git commit -m "feat(cuda): add use_persistent flag to GpuBackend"
```

---

### Task 2: Add PersistentBuffers and declarations to sim_batch.hpp

**Files:**
- Modify: `include/coilgun/simulation/cuda/sim_batch.hpp`

- [ ] **Step 1: Read current file**

Read `include/coilgun/simulation/cuda/sim_batch.hpp`.

- [ ] **Step 2: Add PersistentBuffers struct and private methods**

After the `#include` block, before the `namespace coilgun::simulation::cuda {` line, add:

```cpp
namespace coilgun { namespace simulation { namespace cuda {

/**
 * @brief Mapped memory buffers for persistent kernel execution.
 *
 * Allocated once via cudaHostAllocMapped. Host writes to seps/doorbell;
 * device writes to out_M/out_dM/doorbell. No cudaMemcpy calls needed.
 */
struct PersistentBuffers {
    int*    batch_id      = nullptr;  ///< Monotonic counter — host increments to trigger new batch.
    int*    active_pairs  = nullptr;  ///< Number of valid (stage,filament) pairs in this batch.
    int*    shutdown      = nullptr;  ///< Set to 1 to terminate the persistent kernel.
    int*    doorbell      = nullptr;  ///< [max_pairs] — host writes 1, device writes 0.
    double* seps          = nullptr;  ///< [max_pairs] — separation values set by host.
    double* out_M         = nullptr;  ///< [max_pairs] — mutual inductance results from device.
    double* out_dM        = nullptr;  ///< [max_pairs] — M gradient results from device.

    // Device-side pointer aliases (obtained via cudaHostGetDevicePointer)
    int*    d_batch_id     = nullptr;
    int*    d_active_pairs = nullptr;
    int*    d_shutdown     = nullptr;
    int*    d_doorbell     = nullptr;
    double* d_seps         = nullptr;
    double* d_out_M        = nullptr;
    double* d_out_dM       = nullptr;
};

}}} // namespace coilgun::simulation::cuda
```

Inside the `SimBatch` class private section, add these methods:

```cpp
    void init_persistent_buffers();
    void free_persistent_buffers();
    void launch_persistent_kernel();
```

And add these members:

```cpp
    PersistentBuffers pbuf_;
```

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/sim_batch.hpp
git commit -m "feat(cuda): add PersistentBuffers and persistent kernel declarations"
```

---

### Task 3: Implement persistent kernel buffer management

**Files:**
- Modify: `src/cuda/sim_batch.cu`

- [ ] **Step 1: Read current file**

Read `src/cuda/sim_batch.cu`. Note the current constructor, destructor, and `compute_all_M1_dM1`.

- [ ] **Step 2: Add init/free/launch implementations**

After the anonymous namespace, add these template specializations:

```cpp
template<typename SP>
void SimBatch<SP>::init_persistent_buffers() {
    int N_b = n_stages_ * N_fil_;
    auto ah = [](auto& p, size_t n) {
        cudaError_t e = cudaHostAlloc(&p, n * sizeof(std::remove_pointer_t<decltype(p)>), cudaHostAllocMapped);
        if (e != cudaSuccess)
            throw std::runtime_error(std::string("cudaHostAlloc failed: ") + cudaGetErrorString(e));
    };
    auto gdp = [](auto& d, auto h) {
        cudaHostGetDevicePointer(&d, h, 0);
    };

    ah(pbuf_.batch_id, 1);      gdp(pbuf_.d_batch_id, pbuf_.batch_id);
    ah(pbuf_.active_pairs, 1);  gdp(pbuf_.d_active_pairs, pbuf_.active_pairs);
    ah(pbuf_.shutdown, 1);      gdp(pbuf_.d_shutdown, pbuf_.shutdown);
    ah(pbuf_.doorbell, N_b);    gdp(pbuf_.d_doorbell, pbuf_.doorbell);
    ah(pbuf_.seps, N_b);        gdp(pbuf_.d_seps, pbuf_.seps);
    ah(pbuf_.out_M, N_b);       gdp(pbuf_.d_out_M, pbuf_.out_M);
    ah(pbuf_.out_dM, N_b);      gdp(pbuf_.d_out_dM, pbuf_.out_dM);

    std::memset(pbuf_.doorbell, 0, N_b * sizeof(int));
    std::memset(pbuf_.seps, 0, N_b * sizeof(double));
    *pbuf_.batch_id = 0;
    *pbuf_.active_pairs = 0;
    *pbuf_.shutdown = 0;
}

template<typename SP>
void SimBatch<SP>::free_persistent_buffers() {
    auto fh = [](void*& p) { if (p) { cudaFreeHost(p); p = nullptr; } };
    fh((void*&)pbuf_.batch_id);
    fh((void*&)pbuf_.active_pairs);
    fh((void*&)pbuf_.shutdown);
    fh((void*&)pbuf_.doorbell);
    fh((void*&)pbuf_.seps);
    fh((void*&)pbuf_.out_M);
    fh((void*&)pbuf_.out_dM);
}

template<typename SP>
void SimBatch<SP>::launch_persistent_kernel() {
    int N_b = n_stages_ * N_fil_;
    persistent_batch_kernel<<<N_b, backend_.threads_per_block>>>(
        adaptor_.d_coils(), adaptor_.d_fils(),
        adaptor_.d_nodes(), adaptor_.d_weights(),
        pbuf_.d_batch_id,
        pbuf_.d_active_pairs,
        pbuf_.d_shutdown,
        pbuf_.d_seps,
        pbuf_.d_doorbell,
        pbuf_.d_out_M,
        pbuf_.d_out_dM,
        n_stages_, N_fil_, 9);
}
```

- [ ] **Step 3: Commit**

```sh
git add src/cuda/sim_batch.cu
git commit -m "feat(cuda): implement persistent kernel buffer init/free/launch"
```

---

### Task 4: Implement the persistent_batch_kernel

**Files:**
- Modify: `src/cuda/sim_batch.cu`

- [ ] **Step 1: Add the kernel above the SimBatch class**

Place this kernel definition between the `#include` block and the namespace declarations, in `namespace coilgun::simulation::cuda`:

```cpp
namespace coilgun { namespace simulation { namespace cuda {

/**
 * @brief Persistent batch kernel — processes (stage, filament) pairs via mapped memory.
 *
 * Launched once per SimBatch lifetime. Host triggers work batches via
 * batch_id increment + doorbell writes. Blocks spin-wait between batches.
 * Uses 512-thread shared-memory reduction for 4D GL integration.
 *
 * @param coils, fils    Device arrays of coil/filament geometry (uploaded once by GpuAdaptor).
 * @param gl_nodes       Gauss-Legendre nodes (from adaptor_.d_nodes()).
 * @param gl_weights     Gauss-Legendre weights (from adaptor_.d_weights()).
 * @param batch_id       Host-incremented monotonic counter — triggers new work.
 * @param active_pairs   Number of valid (stage, filament) indices in this batch.
 * @param shutdown       Set to 1 by host to terminate the kernel.
 * @param seps           [N_b] separation values for each (si, fi) pair.
 * @param doorbell       [N_b] — host writes 1 to signal data ready, device writes 0 when done.
 * @param out_M, out_dM  [N_b] — device writes M and dM results.
 * @param n_stages, N_fil, n_nodes  Geometry dimensions.
 */
__global__ void persistent_batch_kernel(
        const CoilGeo* coils, const FilGeo* fils,
        const double* gl_nodes, const double* gl_weights,
        volatile int* batch_id,
        volatile int* active_pairs,
        volatile int* shutdown,
        volatile double* seps,
        volatile int* doorbell,
        volatile double* out_M,
        volatile double* out_dM,
        int n_stages, int N_fil, int n_nodes) {

    int idx = blockIdx.x;              // (stage, filament) pair index
    int N_b = n_stages * N_fil;
    int si = idx / N_fil;              // stage index
    int fi = idx % N_fil;              // filament index
    int my_last_batch = -1;

    const auto& coil = coils[si];
    const auto& fil  = fils[fi];

    // Precompute constant geometry for this block
    const double ra_mid  = 0.5 * (coil.re + coil.ri);
    const double ra_half = 0.5 * (coil.re - coil.ri);
    const double rb_mid  = 0.5 * (fil.re + fil.ri);
    const double rb_half = 0.5 * (fil.re - fil.ri);
    const double la_half = 0.5 * coil.len;
    const double lb_half = 0.5 * fil.len;
    const double prefactor = coil.turns / 16.0;

    int n2 = n_nodes * n_nodes;
    int n4 = n2 * n2;
    int tid = threadIdx.x;
    int total_threads = blockDim.x;

    while (*shutdown == 0) {
        // Wait for next batch
        int cur = *batch_id;
        if (cur == my_last_batch) continue;
        my_last_batch = cur;

        // Wait for host to fill our doorbell
        while (doorbell[idx] == 0) {}

        double M = 0.0, dM = 0.0;

        if (idx < *active_pairs) {
            double separation = seps[idx];

            __shared__ double sM[512];
            __shared__ double sdM[512];
            sM[tid]  = 0.0;
            sdM[tid] = 0.0;

            for (int pt = tid; pt < n4; pt += total_threads) {
                int i1 = pt / (n_nodes * n_nodes * n_nodes);
                int rem1 = pt % (n_nodes * n_nodes * n_nodes);
                int j1 = rem1 / (n_nodes * n_nodes);
                int rem2 = rem1 % (n_nodes * n_nodes);
                int i2 = rem2 / n_nodes;
                int j2 = rem2 % n_nodes;

                double w = gl_weights[i1] * gl_weights[j1]
                         * gl_weights[i2] * gl_weights[j2];

                double ra = ra_mid + ra_half * gl_nodes[i1];
                double rb = rb_mid + rb_half * gl_nodes[i2];
                double za = la_half * gl_nodes[j1];
                double zb = separation + lb_half * gl_nodes[j2];

                double abs_sep = fabs(zb - za);
                sM[tid]  += w * coilgun::physics::mutual_inductance_filament_device(ra, rb, abs_sep);
                sdM[tid] += w * coilgun::physics::mutual_inductance_gradient_filament_device(ra, rb, zb - za);
            }

            __syncthreads();

            for (int stride = total_threads / 2; stride > 0; stride /= 2) {
                if (tid < stride) {
                    sM[tid]  += sM[tid + stride];
                    sdM[tid] += sdM[tid + stride];
                }
                __syncthreads();
            }

            if (tid == 0) {
                M  = prefactor * sM[0];
                dM = prefactor * sdM[0];
            }
        }

        // Write results (tid 0 only, others write 0.0)
        // All threads participate in doorbell sync
        __threadfence();
        if (tid == 0) {
            out_M[idx]  = M;
            out_dM[idx] = dM;
        }
        __threadfence();
        __syncthreads();
        if (tid == 0)
            doorbell[idx] = 0;
    }
}

}}} // namespace coilgun::simulation::cuda
```

- [ ] **Step 2: Build and verify compilation**

```sh
cmake --build build/ninja-debug --target coilgun_cuda
```

Expected: 0 errors.

- [ ] **Step 3: Commit**

```sh
git add src/cuda/sim_batch.cu
git commit -m "feat(cuda): implement persistent_batch_kernel"
```

---

### Task 5: Wire persistent kernel into SimBatch constructor/destructor

**Files:**
- Modify: `src/cuda/sim_batch.cu`

- [ ] **Step 1: Modify constructor to use persistent or fallback**

At the end of the `SimBatch` constructor, after the existing initialization, add:

```cpp
    if (backend_.use_persistent) {
        init_persistent_buffers();
        launch_persistent_kernel();
    }
```

- [ ] **Step 2: Modify destructor to clean up persistent buffers**

Change `SimBatch<SP>::~SimBatch() = default;` to:

```cpp
template<typename SP>
SimBatch<SP>::~SimBatch() {
    if (backend_.use_persistent) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers();
    }
}
```

- [ ] **Step 3: Commit**

```sh
git add src/cuda/sim_batch.cu
git commit -m "feat(cuda): wire persistent kernel into SimBatch lifecycle"
```

---

### Task 6: Implement persistent compute_all_M1_dM1

**Files:**
- Modify: `src/cuda/sim_batch.cu`

- [ ] **Step 1: Replace compute_all_M1_dM1 with persistent version**

When `backend_.use_persistent` is true, use the mapped memory path. When false, use the existing per-sim kernel loop. Replace the entire function body:

```cpp
template<typename SP>
void SimBatch<SP>::compute_all_M1_dM1() {
    int S = n_stages_, F = N_fil_, N_b = S * F;
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();

    if (backend_.use_persistent) {
        // === Persistent kernel path ===
        int active_count = 0;

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
                if (!sim.triggered[si] || sim.finished[si]) continue;
                for (int fi = 0; fi < F; ++fi) {
                    int i = fi / nr + 1, j = fi % nr + 1;
                    double z_rel = armature_.filament_axial_position(i) - arm_center_init;
                    double z_global = arm_center + z_rel;
                    int idx = si * F + fi;
                    pbuf_.seps[idx] = z_global - coils_[si].position();
                    ++active_count;
                }
            }

            // Trigger kernel batch
            *pbuf_.active_pairs = active_count;
            for (int i = 0; i < active_count; ++i)
                pbuf_.doorbell[i] = 1;
            (*pbuf_.batch_id)++;

            // Wait for device completion
            for (int i = 0; i < active_count; ++i)
                while (pbuf_.doorbell[i] != 0) {}
            cudaStreamSynchronize(0);

            // Read results
            for (int si = 0; si < S; ++si)
                for (int fi = 0; fi < F; ++fi) {
                    int idx = si * F + fi;
                    batch_M1_(s, idx)  = pbuf_.out_M[idx];
                    batch_dM1_(s, idx) = pbuf_.out_dM[idx];
                }
        }
    } else {
        // === Fallback: per-sim kernel loop (existing code, unchanged) ===
        for (int s = 0; s < num_sims_; ++s) {
            auto& sim = sims_[s];
            sim.all_finished = true;
            for (int i = 0; i < S; ++i)
                if (sim.triggered[i] && !sim.finished[i]) { sim.all_finished = false; break; }
            if (sim.all_finished) continue;

            double arm_center = sim.state.arm_position;
            for (int si = 0; si < S; ++si) {
                if (!sim.triggered[si] || sim.finished[si]) continue;
                const auto& coil = coils_[si];
                for (int fi = 0; fi < F; ++fi) {
                    int i = fi / nr + 1, j = fi % nr + 1;
                    double z_rel = armature_.filament_axial_position(i) - arm_center_init;
                    double z_global = arm_center + z_rel;
                    double sep = z_global - coil.position();
                    double fil_ri = armature_.filament_inner_radius(j);
                    double fil_re = armature_.filament_outer_radius(j);
                    double fil_l  = armature_.length() / armature_.axial_filaments();
                    int idx = si * F + fi;
                    coilgun::physics::mutual_inductance_coil_pair_kernel<<<1, backend_.threads_per_block>>>(
                        coil.inner_radius(), coil.outer_radius(), coil.length(), coil.turns(),
                        fil_ri, fil_re, fil_l, 1, sep, 9,
                        adaptor_.d_results_M() + idx, adaptor_.d_results_dM() + idx);
                }
            }
            cudaDeviceSynchronize();
            std::vector<double> hM(S * F), hdM(S * F);
            cudaMemcpy(hM.data(),  adaptor_.d_results_M(),  S * F * sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(hdM.data(), adaptor_.d_results_dM(), S * F * sizeof(double), cudaMemcpyDeviceToHost);
            for (int si = 0; si < S; ++si)
                for (int fi = 0; fi < F; ++fi) {
                    int idx = si * F + fi;
                    batch_M1_(s, idx) = hM[idx];
                    batch_dM1_(s, idx) = hdM[idx];
                }
        }
    }
}
```

- [ ] **Step 2: Build and verify**

```sh
cmake --build build/ninja-debug --target coilgun_cuda
```

Expected: 0 errors.

- [ ] **Step 3: Commit**

```sh
git add src/cuda/sim_batch.cu
git commit -m "feat(cuda): implement persistent compute_all_M1_dM1 with fallback"
```

---

### Task 7: Update test_gpu_sim_batch for persistent mode

**Files:**
- Modify: `tests/test_gpu_sim_batch.cpp`

- [ ] **Step 1: Read current test file**

Read `tests/test_gpu_sim_batch.cpp`.

- [ ] **Step 2: Add a test case for persistent mode (default)**

Add after the existing test cases:

```cpp
TEST_CASE("SimBatch — persistent kernel mode (default)") {
    auto coils = make_coils();
    auto arm   = make_arm();
    coilgun::simulation::cuda::GpuBackend backend;
    backend.use_persistent = true;
    backend.threads_per_block = 512;

    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6, backend);

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    batch.set_excitations(0, std::move(excs), triggers);

    auto pol = coilgun::simulation::TerminationPolicy::defaults();
    pol.max_steps = 500;
    batch.run(pol);
    double v = batch.result(0).summary.muzzle_velocity;
    CHECK(v > 0.0);
}

TEST_CASE("SimBatch — fallback kernel mode") {
    auto coils = make_coils();
    auto arm   = make_arm();
    coilgun::simulation::cuda::GpuBackend backend;
    backend.use_persistent = false;  // Force fallback

    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6, backend);

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    batch.set_excitations(0, std::move(excs), triggers);

    auto pol = coilgun::simulation::TerminationPolicy::defaults();
    pol.max_steps = 500;
    batch.run(pol);
    double v = batch.result(0).summary.muzzle_velocity;
    CHECK(v > 0.0);
}
```

- [ ] **Step 3: Build and run tests**

```sh
cmake --build build/ninja-debug --target test_gpu_sim_batch
./build/ninja-debug/tests/test_gpu_sim_batch
```

Expected: all test cases pass (persistent + fallback match).

- [ ] **Step 4: Run with compute-sanitizer**

```sh
compute-sanitizer ./build/ninja-debug/tests/test_gpu_sim_batch
```

Expected: 0 errors.

- [ ] **Step 5: Commit**

```sh
git add tests/test_gpu_sim_batch.cpp
git commit -m "test(cuda): add persistent and fallback mode tests for SimBatch"
```

---

### Task 8: Final verification — full test suite

**Files:** None (verification only)

- [ ] **Step 1: Rebuild everything**

```sh
cmake --build build/ninja-debug
```

Expected: 0 errors.

- [ ] **Step 2: Run all GPU tests**

```sh
ctest --test-dir build/ninja-debug -R "test_gpu"
```

Expected: all pass.

- [ ] **Step 3: Verify callable with correct include**

Build a minimal example:

```cpp
#include <coilgun/coilgun_cuda.hpp>

int main() {
    coilgun::simulation::cuda::GpuBackend be;
    be.use_persistent = true;
    return 0;
}
```

Compile with:
```sh
g++ -std=c++17 -Iinclude -Ibuild/ninja-debug/_deps/eigen-src -o /tmp/quick_test test.cpp
```

Expected: compiles clean.

- [ ] **Step 4: Commit**

```sh
git commit -m "build: verify full test suite passes with persistent kernel" --allow-empty
```
