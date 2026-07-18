# Persistent Kernel — Implementation Plan (v2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace per-step kernel launch loop in all three GPU simulation classes (GpuSingleStageSim, GpuMultiStageSim, SimBatch) with a single persistent kernel using mapped memory, achieving zero kernel launches per step. Expand GpuOptLevel from two to three levels (Standard/Full/Aggressive). Auto-fallback to per-pair kernel launch on cudaHostAllocMapped failure.

**Architecture:** Shared `persistent_kernel.cuh` header defines PersistentBuffers struct and `persistent_batch_kernel`. Included by all three `.cu` files via `#include`. Doorbell/batch_id protocol for host-device sync. `GpuBackend::use_persistent` bool controls persistent vs. fallback mode. `GpuOptLevel::Aggressive` added as third level (FP32 integrand + FP64 reduction).

**Tech Stack:** C++17, CUDA 13.3, Eigen 3.4, Boost.Math, doctest.

**Scope:** Tasks 1-12 cover D1-D4 (persistent kernel for all three classes). D5 (double buffering) and D6 (Aggressive FP32) are separate follow-up plans per the spec.

---

## File Map

### Create
```
include/coilgun/simulation/cuda/persistent_kernel.cuh    — PersistentBuffers + kernel + templates
```

### Modify
```
include/coilgun/simulation/cuda/gpu_backend.hpp           — use_persistent flag + GpuOptLevel extension
include/coilgun/simulation/cuda/gpu_single_stage_sim.hpp  — PersistentBuffers member + declarations
include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp   — PersistentBuffers member + declarations
include/coilgun/simulation/cuda/sim_batch.hpp             — PersistentBuffers member + declarations
src/cuda/gpu_single_stage_sim.cu                          — persistent compute_M1_dM1 + lifecycle
src/cuda/gpu_multi_stage_sim.cu                           — persistent compute_M1_dM1 + lifecycle
src/cuda/sim_batch.cu                                     — persistent compute_all_M1_dM1 + lifecycle
src/cuda/CMakeLists.txt                                   — (no changes needed)
```

### Test Extend
```
tests/test_gpu_sim_batch.cpp                              — persistent/fallback consistency checks
```

### Follow-up Plans (not in this plan)
```
docs/superpowers/plans/YYYY-MM-DD-double-buffering.md     — D5: double-buffered CPU/GPU overlap
docs/superpowers/plans/YYYY-MM-DD-fp32-aggressive.md     — D6: FP32 integrand + AGM elliptic integrals
```

---

### Task 1: Update GpuBackend and GpuOptLevel

**Files:**
- Modify: `include/coilgun/simulation/cuda/gpu_backend.hpp`

- [ ] **Step 1: Read current file**

Read `include/coilgun/simulation/cuda/gpu_backend.hpp`.

- [ ] **Step 2: Rewrite file with use_persistent flag and GpuOptLevel expansion**

Replace the entire file content with:

```cpp
/**
 * @file gpu_backend.hpp
 * @brief GPU backend types and configuration.
 * @author Winston Meursault
 */

#pragma once

#include <cstddef>

namespace coilgun::simulation::cuda {

/**
 * @brief GPU optimization level. Levels are cumulative — each includes
 *        all optimisations of lower levels.
 *
 * | Level      | Precision                | Distance cutoff | GL order | Use case           |
 * |:----------:|:------------------------:|:---------------:|:--------:|--------------------|
 * | Standard   | FP64                     | No              | 9        | Validation/debug   |
 * | Full       | FP64                     | Yes (>10x coil) | 9        | Production default |
 * | Aggressive | FP32 integrand + FP64 red.| Yes (>10x coil) | 9        | Large-scale sweeps |
 */
enum class GpuOptLevel {
    Standard   = 0,   ///< FP64, no distance cutoff, n_nodes=9.
    Full       = 1,   ///< FP64, distance cutoff, n_nodes=9.
    Aggressive = 2,   ///< FP32 integrand, FP64 reduction, distance cutoff, n_nodes=9.
};

/**
 * @brief GPU backend configuration.
 */
struct GpuBackend {
    int     device_id         = 0;    ///< cudaSetDevice target.
    int     threads_per_block = 512;  ///< Threads per block for integration kernel.
    size_t  max_batch_sims    = 256;  ///< Pre-allocated buffer size for batch mode.
    bool    enable_profiling  = false; ///< Enable NVTX range annotations.
    bool    use_persistent    = true;  ///< Use persistent kernel (mapped memory). Falls back to per-pair launches if false or if cudaHostAllocMapped fails.
};

} // namespace coilgun::simulation::cuda
```

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/gpu_backend.hpp
git commit -m "feat(cuda): add use_persistent flag and GpuOptLevel::Aggressive to GpuBackend"
```

---

### Task 2: Create persistent_kernel.cuh — PersistentBuffers

**Files:**
- Create: `include/coilgun/simulation/cuda/persistent_kernel.cuh`

- [ ] **Step 1: Create the file**

Create `include/coilgun/simulation/cuda/persistent_kernel.cuh`:

```cpp
/**
 * @file persistent_kernel.cuh
 * @brief Persistent kernel infrastructure — shared by all three GPU simulation classes.
 * @author Winston Meursault
 *
 * Defines PersistentBuffers (mapped-memory doorbell protocol) and
 * persistent_batch_kernel. Included by gpu_single_stage_sim.cu,
 * gpu_multi_stage_sim.cu, and sim_batch.cu.
 */

#pragma once

#ifdef __CUDACC__

#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/physics/mutual_inductance.cuh"
#include <cuda_runtime.h>
#include <cstring>
#include <stdexcept>
#include <string>
#include <iostream>

namespace coilgun::simulation::cuda {

/**
 * @brief Mapped memory buffers for persistent kernel execution.
 *
 * Allocated once via cudaHostAllocMapped. Host writes seps/doorbell;
 * device writes out_M/out_dM/doorbell. No cudaMemcpy calls needed.
 */
struct PersistentBuffers {
    int*    batch_id      = nullptr;  ///< Monotonic counter — host increments to trigger new batch.
    int*    active_pairs  = nullptr;  ///< Number of valid (stage, filament) indices in this batch.
    int*    shutdown      = nullptr;  ///< Set to 1 to terminate the persistent kernel.
    int*    doorbell      = nullptr;  ///< [max_pairs] — host writes 1, device writes 0.
    double* seps          = nullptr;  ///< [max_pairs] — separation values set by host.
    double* out_M         = nullptr;  ///< [max_pairs] — mutual inductance results from device.
    double* out_dM        = nullptr;  ///< [max_pairs] — M gradient results from device.

    // Device-side pointer aliases (obtained via cudaHostGetDevicePointer)
    int*    d_batch_id      = nullptr;
    int*    d_active_pairs  = nullptr;
    int*    d_shutdown      = nullptr;
    int*    d_doorbell      = nullptr;
    double* d_seps          = nullptr;
    double* d_out_M         = nullptr;
    double* d_out_dM        = nullptr;
};

/**
 * @brief Allocate mapped memory buffers for persistent kernel.
 * @param pbuf PersistentBuffers to initialise.
 * @param N_b Number of (stage, filament) pairs = n_stages × N_fil.
 * @param backend GPU backend config (for fallback on failure).
 * @return true on success, false if cudaHostAllocMapped failed (use_persistent should be set to false).
 */
inline bool init_persistent_buffers(PersistentBuffers& pbuf, int N_b, GpuBackend& backend) {
    auto ah = [&](auto& p, size_t n) -> bool {
        cudaError_t e = cudaHostAlloc(&p, n * sizeof(std::remove_pointer_t<decltype(p)>), cudaHostAllocMapped);
        if (e != cudaSuccess) {
            std::cerr << "[coilgun] cudaHostAllocMapped failed: " << cudaGetErrorString(e)
                      << " — falling back to per-pair kernel launches." << std::endl;
            return false;
        }
        return true;
    };
    auto gdp = [](auto& d, auto h) {
        cudaHostGetDevicePointer(&d, h, 0);
    };

    if (!ah(pbuf.batch_id, 1))      return false;
    gdp(pbuf.d_batch_id, pbuf.batch_id);
    if (!ah(pbuf.active_pairs, 1))  return false;
    gdp(pbuf.d_active_pairs, pbuf.active_pairs);
    if (!ah(pbuf.shutdown, 1))      return false;
    gdp(pbuf.d_shutdown, pbuf.shutdown);
    if (!ah(pbuf.doorbell, N_b))    return false;
    gdp(pbuf.d_doorbell, pbuf.doorbell);
    if (!ah(pbuf.seps, N_b))        return false;
    gdp(pbuf.d_seps, pbuf.seps);
    if (!ah(pbuf.out_M, N_b))       return false;
    gdp(pbuf.d_out_M, pbuf.out_M);
    if (!ah(pbuf.out_dM, N_b))      return false;
    gdp(pbuf.d_out_dM, pbuf.out_dM);

    std::memset(pbuf.doorbell, 0, N_b * sizeof(int));
    std::memset(pbuf.seps, 0, N_b * sizeof(double));
    *pbuf.batch_id = 0;
    *pbuf.active_pairs = 0;
    *pbuf.shutdown = 0;
    return true;
}

/**
 * @brief Free mapped memory buffers.
 */
inline void free_persistent_buffers(PersistentBuffers& pbuf) {
    auto fh = [](void*& p) { if (p) { cudaFreeHost(p); p = nullptr; } };
    fh((void*&)pbuf.batch_id);
    fh((void*&)pbuf.active_pairs);
    fh((void*&)pbuf.shutdown);
    fh((void*&)pbuf.doorbell);
    fh((void*&)pbuf.seps);
    fh((void*&)pbuf.out_M);
    fh((void*&)pbuf.out_dM);
}

} // namespace coilgun::simulation::cuda

#endif // __CUDACC__
```

- [ ] **Step 2: Commit**

```sh
git add include/coilgun/simulation/cuda/persistent_kernel.cuh
git commit -m "feat(cuda): add PersistentBuffers with auto-fallback to persistent_kernel.cuh"
```

---

### Task 3: Add persistent_batch_kernel to persistent_kernel.cuh

**Files:**
- Modify: `include/coilgun/simulation/cuda/persistent_kernel.cuh`

- [ ] **Step 1: Append the kernel after the inline functions**

Append the following code after the `free_persistent_buffers` function, before `#endif // __CUDACC__`:

```cpp
/**
 * @brief Launch the persistent kernel (one grid per class lifetime).
 * @param pbuf Initialised PersistentBuffers.
 * @param adaptor GpuAdaptor with already-uploaded geometry.
 * @param N_b Number of (stage, filament) pairs = n_stages × N_fil.
 * @param threads_per_block Threads per block (must be power of 2, default 512).
 * @param n_nodes Number of GL quadrature nodes per dimension (typ. 9).
 */
inline void launch_persistent_kernel(
        const PersistentBuffers& pbuf,
        const GpuAdaptor& adaptor,
        int N_b, int threads_per_block, int n_nodes);

/**
 * @brief Persistent batch kernel — processes (stage, filament) pairs via mapped memory.
 *
 * Launched once per class lifetime. Host triggers work batches via
 * batch_id increment + doorbell writes. Blocks spin-wait between batches.
 * Uses shared-memory reduction for 4D GL integration.
 *
 * @param coils, fils    Device arrays of coil/filament geometry (uploaded once by GpuAdaptor).
 * @param gl_nodes       Gauss-Legendre nodes (from adaptor.d_nodes()).
 * @param gl_weights     Gauss-Legendre weights (from adaptor.d_weights()).
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

        // Write results (tid 0 only; all threads participate in fence)
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

/**
 * @brief Implementation of launch_persistent_kernel.
 * Must be in a .cu file (not inlined in header) because it uses <<<>>> syntax.
 * Template specializations in each .cu file call this.
 */
inline void launch_persistent_kernel_impl(
        const PersistentBuffers& pbuf,
        const GpuAdaptor& adaptor,
        int N_b, int threads_per_block, int n_nodes) {
    persistent_batch_kernel<<<N_b, threads_per_block>>>(
        adaptor.d_coils(), adaptor.d_fils(),
        adaptor.d_nodes(), adaptor.d_weights(),
        pbuf.d_batch_id,
        pbuf.d_active_pairs,
        pbuf.d_shutdown,
        pbuf.d_seps,
        pbuf.d_doorbell,
        pbuf.d_out_M,
        pbuf.d_out_dM,
        adaptor.n_stages(), N_b / adaptor.n_stages(), n_nodes);
}
```

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/persistent_kernel.cuh
git commit -m "feat(cuda): add persistent_batch_kernel to persistent_kernel.cuh"
```

---

### Task 4: Integrate persistent kernel into GpuSingleStageSim — header

**Files:**
- Modify: `include/coilgun/simulation/cuda/gpu_single_stage_sim.hpp`

- [ ] **Step 1: Read current file**

Read `include/coilgun/simulation/cuda/gpu_single_stage_sim.hpp`.

- [ ] **Step 2: Add PersistentBuffers member and method declarations**

After the `#include` block, add:
```cpp
#include "coilgun/simulation/cuda/persistent_kernel.cuh"
```

In the private section of `GpuSingleStageSim`, add these members:
```cpp
    bool               use_persistent_ = false;
    PersistentBuffers  pbuf_;
```

And add these private method declarations:
```cpp
    void init_persistent_mode();
    void compute_M1_dM1_persistent();
    void compute_M1_dM1_fallback();
```

Note: Read the exact current private section structure of the header to place these correctly. The declarations should be placed after the existing `GpuAdaptor adaptor_;` member and before the closing `};`.

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/gpu_single_stage_sim.hpp
git commit -m "feat(cuda): add persistent kernel declarations to GpuSingleStageSim"
```

---

### Task 5: Integrate persistent kernel into GpuSingleStageSim — implementation

**Files:**
- Modify: `src/cuda/gpu_single_stage_sim.cu`

- [ ] **Step 1: Read current file**

Read `src/cuda/gpu_single_stage_sim.cu`. Note the current constructor, destructor, `compute_M1_dM1`, and `step` methods.

- [ ] **Step 2: Add init_persistent_mode method**

Add this method implementation after the existing constructor definition:

```cpp
template<typename SP>
void GpuSingleStageSim<SP>::init_persistent_mode() {
    if (backend_.use_persistent) {
        int N_fil = armature_.total_filaments();
        // GpuAdaptor already uploaded at construction, but not in batch mode
        // For persistent kernel we use the existing adaptor setup.
        use_persistent_ = init_persistent_buffers(pbuf_, N_fil, backend_);
        if (use_persistent_) {
            launch_persistent_kernel_impl(pbuf_, adaptor_, N_fil,
                                          backend_.threads_per_block, 9);
        }
    }
}
```

- [ ] **Step 3: Wire init_persistent_mode into constructor**

At the end of the constructor body, after `adaptor_.setup(...)`, add:
```cpp
    init_persistent_mode();
```

- [ ] **Step 4: Add destructor cleanup**

Replace the existing destructor (if `= default` or empty) with:
```cpp
template<typename SP>
GpuSingleStageSim<SP>::~GpuSingleStageSim() {
    if (use_persistent_) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers(pbuf_);
    }
}
```

If the existing destructor had other cleanup, preserve that code and add the persistent cleanup block.

- [ ] **Step 5: Add compute_M1_dM1_persistent method**

Add this method that uses the mapped memory path:

```cpp
template<typename SP>
void GpuSingleStageSim<SP>::compute_M1_dM1_persistent() {
    int F = armature_.total_filaments();
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();

    double arm_center = state_.arm_position;
    for (int fi = 0; fi < F; ++fi) {
        int i = fi / nr + 1, j = fi % nr + 1;
        double z_rel = armature_.filament_axial_position(i) - arm_center_init;
        pbuf_.seps[fi] = arm_center + z_rel - coil_.position();
    }

    *pbuf_.active_pairs = F;
    for (int fi = 0; fi < F; ++fi)
        pbuf_.doorbell[fi] = 1;
    (*pbuf_.batch_id)++;

    for (int fi = 0; fi < F; ++fi)
        while (pbuf_.doorbell[fi] != 0) {}
    cudaStreamSynchronize(0);

    M1_mat_.resize(1, F);
    dM1_mat_.resize(1, F);
    for (int fi = 0; fi < F; ++fi) {
        M1_mat_(0, fi)  = pbuf_.out_M[fi];
        dM1_mat_(0, fi) = pbuf_.out_dM[fi];
    }
}
```

- [ ] **Step 6: Add compute_M1_dM1_fallback method**

Extract the existing per-pair kernel launch loop from the current `compute_M1_dM1` into this method. Study the current `compute_M1_dM1` in `gpu_single_stage_sim.cu` and move its body here:

```cpp
template<typename SP>
void GpuSingleStageSim<SP>::compute_M1_dM1_fallback() {
    // [existing compute_M1_dM1 body — per-pair kernel launch loop]
    // Copy the EXACT existing implementation here verbatim.
}
```

- [ ] **Step 7: Rewrite compute_M1_dM1 as dispatcher**

Replace `compute_M1_dM1` with:
```cpp
template<typename SP>
void GpuSingleStageSim<SP>::compute_M1_dM1() {
    if (use_persistent_)
        compute_M1_dM1_persistent();
    else
        compute_M1_dM1_fallback();
}
```

- [ ] **Step 8: Build and verify compilation**

```sh
cmake --build build/ninja-cuda-debug --target coilgun_cuda
```
Expected: 0 errors.

- [ ] **Step 9: Commit**

```sh
git add src/cuda/gpu_single_stage_sim.cu
git commit -m "feat(cuda): implement persistent kernel path in GpuSingleStageSim"
```

---

### Task 6: Test GpuSingleStageSim persistent mode

**Files:**
- Modify: `tests/test_gpu_vs_cpu_single.cpp`

- [ ] **Step 1: Read current test file**

Read `tests/test_gpu_vs_cpu_single.cpp`. Note the existing test structure.

- [ ] **Step 2: Add persistent vs fallback consistency test**

Add after the existing test cases:

```cpp
TEST_CASE("GpuSingleStageSim — persistent vs fallback consistency") {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::GpuSingleStageSim;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;

    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    // Run with persistent kernel (default)
    coilgun::simulation::cuda::GpuBackend be_persistent;
    be_persistent.use_persistent = true;
    auto exc1 = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim_p(coil, arm, std::move(exc1), 1e-6,
                                          false, GpuOptLevel::Full, be_persistent);
    sim_p.run();
    double v_persistent = sim_p.result().summary.muzzle_velocity;

    // Run with fallback kernel
    coilgun::simulation::cuda::GpuBackend be_fallback;
    be_fallback.use_persistent = false;
    auto exc2 = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim_f(coil, arm, std::move(exc2), 1e-6,
                                          false, GpuOptLevel::Full, be_fallback);
    sim_f.run();
    double v_fallback = sim_f.result().summary.muzzle_velocity;

    CHECK(v_persistent > 0.0);
    CHECK(v_fallback > 0.0);
    CHECK(v_persistent == doctest::Approx(v_fallback).epsilon(1e-12));
}
```

- [ ] **Step 3: Build and run**

```sh
cmake --build build/ninja-cuda-debug --target test_gpu_vs_cpu_single
./build/ninja-cuda-debug/tests/test_gpu_vs_cpu_single
```
Expected: all test cases pass, including the new one.

- [ ] **Step 4: Commit**

```sh
git add tests/test_gpu_vs_cpu_single.cpp
git commit -m "test(cuda): add persistent vs fallback consistency test for GpuSingleStageSim"
```

---

### Task 7: Integrate persistent kernel into GpuMultiStageSim — header

**Files:**
- Modify: `include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp`

- [ ] **Step 1: Read current file**

Read `include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp`.

- [ ] **Step 2: Add include and PersistentBuffers member**

After the existing `#include` block, add:
```cpp
#include "coilgun/simulation/cuda/persistent_kernel.cuh"
```

In the private section of `GpuMultiStageSim`, add:
```cpp
    bool               use_persistent_ = false;
    PersistentBuffers  pbuf_;
```

And add these private method declarations:
```cpp
    void init_persistent_mode();
    void compute_M1_dM1_persistent();
    void compute_M1_dM1_fallback();
```

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp
git commit -m "feat(cuda): add persistent kernel declarations to GpuMultiStageSim"
```

---

### Task 8: Integrate persistent kernel into GpuMultiStageSim — implementation

**Files:**
- Modify: `src/cuda/gpu_multi_stage_sim.cu`

- [ ] **Step 1: Read current file**

Read `src/cuda/gpu_multi_stage_sim.cu`. Note: (1) the constructor, (2) the destructor, (3) `compute_M1_dM1`, (4) `is_stage_within_range`, and (5) the existing per-pair kernel launch loop. Pay close attention to the distance cutoff logic in `GpuOptLevel::Full` mode.

- [ ] **Step 2: Add init_persistent_mode method**

```cpp
template<typename SP>
void GpuMultiStageSim<SP>::init_persistent_mode() {
    if (backend_.use_persistent) {
        int N_b = n_stages_ * N_fil_;
        use_persistent_ = init_persistent_buffers(pbuf_, N_b, backend_);
        if (use_persistent_) {
            launch_persistent_kernel_impl(pbuf_, adaptor_, N_b,
                                          backend_.threads_per_block, 9);
        }
    }
}
```

- [ ] **Step 3: Wire into constructor**

At the end of the constructor body, after `adaptor_.setup(...)`, add:
```cpp
    init_persistent_mode();
```

- [ ] **Step 4: Add destructor cleanup**

Replace the existing destructor with:
```cpp
template<typename SP>
GpuMultiStageSim<SP>::~GpuMultiStageSim() {
    if (use_persistent_) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers(pbuf_);
    }
}
```

If the existing destructor had other cleanup code, preserve it and add the persistent cleanup block.

- [ ] **Step 5: Add compute_M1_dM1_persistent**

```cpp
template<typename SP>
void GpuMultiStageSim<SP>::compute_M1_dM1_persistent() {
    int S = n_stages_, F = N_fil_;
    int nr = armature_.radial_filaments();
    double arm_center_init = armature_.position();
    double arm_center = state_.arm_position;

    int active_count = 0;
    for (int si = 0; si < S; ++si) {
        if (!triggered_[si] || finished_[si]) continue;
        if (opt_level_ >= GpuOptLevel::Full && !is_stage_within_range(si)) continue;

        for (int fi = 0; fi < F; ++fi) {
            int i = fi / nr + 1, j = fi % nr + 1;
            double z_rel = armature_.filament_axial_position(i) - arm_center_init;
            int idx = si * F + fi;
            pbuf_.seps[idx] = arm_center + z_rel - coils_[si].position();
            ++active_count;
        }
    }

    *pbuf_.active_pairs = active_count;
    for (int i = 0; i < active_count; ++i)
        pbuf_.doorbell[i] = 1;
    (*pbuf_.batch_id)++;

    for (int i = 0; i < active_count; ++i)
        while (pbuf_.doorbell[i] != 0) {}
    cudaStreamSynchronize(0);

    M1_mat_.resize(S, F);
    dM1_mat_.resize(S, F);
    for (int si = 0; si < S; ++si) {
        for (int fi = 0; fi < F; ++fi) {
            int idx = si * F + fi;
            M1_mat_(si, fi)  = pbuf_.out_M[idx];
            dM1_mat_(si, fi) = pbuf_.out_dM[idx];
        }
    }
}
```

**Important**: The persistent kernel fills ALL N_b slots. Inactive/skipped stages' blocks (idx beyond active_pairs or not in the seps loop) output M=0, dM=0. The host reads these zeros which is correct for downstream computation (force summation with inactive stages yields 0).

The `active_count` in the persistent path is set to the number of pairs the host actually fills separations for. Idle blocks (idx >= active_count) still participate in doorbell sync but skip computation.

- [ ] **Step 6: Extract existing code into compute_M1_dM1_fallback**

Extract the EXISTING per-pair kernel launch loop from the current `compute_M1_dM1` into `compute_M1_dM1_fallback`. Copy the body verbatim with no changes.

- [ ] **Step 7: Rewrite compute_M1_dM1 as dispatcher**

```cpp
template<typename SP>
void GpuMultiStageSim<SP>::compute_M1_dM1() {
    if (use_persistent_)
        compute_M1_dM1_persistent();
    else
        compute_M1_dM1_fallback();
}
```

- [ ] **Step 8: Build and verify**

```sh
cmake --build build/ninja-cuda-debug --target coilgun_cuda
```
Expected: 0 errors.

- [ ] **Step 9: Commit**

```sh
git add src/cuda/gpu_multi_stage_sim.cu
git commit -m "feat(cuda): implement persistent kernel path in GpuMultiStageSim"
```

---

### Task 9: Test GpuMultiStageSim persistent mode

**Files:**
- Modify: `tests/test_gpu_vs_cpu_multi.cpp`

- [ ] **Step 1: Read current test file**

Read `tests/test_gpu_vs_cpu_multi.cpp`. Note the existing test structure.

- [ ] **Step 2: Add persistent vs fallback consistency test**

Add after existing test cases:

```cpp
TEST_CASE("GpuMultiStageSim — persistent vs fallback consistency") {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::GpuMultiStageSim;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;

    DrivingCoil c1(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    auto make_excs = []() -> std::vector<std::unique_ptr<coilgun::simulation::Excitation>> {
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
        return excs;
    };
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};

    // Persistent
    coilgun::simulation::cuda::GpuBackend be_p;
    be_p.use_persistent = true;
    std::vector<DrivingCoil> coils_p = {c1, c2};
    GpuMultiStageSim<EulerStepper> sim_p(
        std::move(coils_p), arm, make_excs(), triggers, 1e-6,
        false, GpuOptLevel::Full, be_p);
    sim_p.run();
    double v_p = sim_p.result().summary.muzzle_velocity;

    // Fallback
    coilgun::simulation::cuda::GpuBackend be_f;
    be_f.use_persistent = false;
    Armature arm_f(0.005, 0.025, 0.08,
                   ALUMINUM.resistivity_ref, ALUMINUM.density,
                   0.0, 0.120, 5, 2, 0.05);
    DrivingCoil c1f(0.01, 0.03, 0.05, 150,
                    COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2f(0.01, 0.03, 0.05, 150,
                    COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    std::vector<DrivingCoil> coils_f = {c1f, c2f};
    GpuMultiStageSim<EulerStepper> sim_f(
        std::move(coils_f), arm_f, make_excs(), triggers, 1e-6,
        false, GpuOptLevel::Full, be_f);
    sim_f.run();
    double v_f = sim_f.result().summary.muzzle_velocity;

    CHECK(v_p > 0.0);
    CHECK(v_f > 0.0);
    CHECK(v_p == doctest::Approx(v_f).epsilon(1e-12));
}
```

- [ ] **Step 3: Build and run**

```sh
cmake --build build/ninja-cuda-debug --target test_gpu_vs_cpu_multi
./build/ninja-cuda-debug/tests/test_gpu_vs_cpu_multi
```
Expected: all test cases pass.

- [ ] **Step 4: Commit**

```sh
git add tests/test_gpu_vs_cpu_multi.cpp
git commit -m "test(cuda): add persistent vs fallback consistency test for GpuMultiStageSim"
```

---

### Task 10: Integrate persistent kernel into SimBatch — header

**Files:**
- Modify: `include/coilgun/simulation/cuda/sim_batch.hpp`

- [ ] **Step 1: Read current file**

Read `include/coilgun/simulation/cuda/sim_batch.hpp`.

- [ ] **Step 2: Add include and PersistentBuffers member**

After the existing `#include` block, add:
```cpp
#include "coilgun/simulation/cuda/persistent_kernel.cuh"
```

In the private section, add these members:
```cpp
    bool               use_persistent_ = false;
    PersistentBuffers  pbuf_;
```

And add these private method declarations:
```cpp
    void init_persistent_mode();
    void compute_all_M1_dM1_persistent();
    void compute_all_M1_dM1_fallback();
```

- [ ] **Step 3: Commit**

```sh
git add include/coilgun/simulation/cuda/sim_batch.hpp
git commit -m "feat(cuda): add persistent kernel declarations to SimBatch"
```

---

### Task 11: Integrate persistent kernel into SimBatch — implementation

**Files:**
- Modify: `src/cuda/sim_batch.cu`

- [ ] **Step 1: Read current file**

Read `src/cuda/sim_batch.cu`. Note the constructor, destructor (`= default`), and `compute_all_M1_dM1`.

- [ ] **Step 2: Add init_persistent_mode method**

```cpp
template<typename SP>
void SimBatch<SP>::init_persistent_mode() {
    if (backend_.use_persistent) {
        int N_b = n_stages_ * N_fil_;
        use_persistent_ = init_persistent_buffers(pbuf_, N_b, backend_);
        if (use_persistent_) {
            launch_persistent_kernel_impl(pbuf_, adaptor_, N_b,
                                          backend_.threads_per_block, 9);
        }
    }
}
```

- [ ] **Step 3: Wire into constructor**

At the end of the constructor body, after `adaptor_.setup(...)`, add:
```cpp
    init_persistent_mode();
```

- [ ] **Step 4: Replace default destructor**

Replace `template<typename SP> SimBatch<SP>::~SimBatch() = default;` with:
```cpp
template<typename SP>
SimBatch<SP>::~SimBatch() {
    if (use_persistent_) {
        *pbuf_.shutdown = 1;
        cudaDeviceSynchronize();
        free_persistent_buffers(pbuf_);
    }
}
```

- [ ] **Step 5: Extract existing compute_all_M1_dM1 into fallback method**

Rename the current `compute_all_M1_dM1` to `compute_all_M1_dM1_fallback`. The body remains EXACTLY as-is.

- [ ] **Step 6: Add persistent compute_all_M1_dM1**

```cpp
template<typename SP>
void SimBatch<SP>::compute_all_M1_dM1_persistent() {
    int S = n_stages_, F = N_fil_;
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
        int active_count = 0;
        for (int si = 0; si < S; ++si) {
            if (!sim.triggered[si] || sim.finished[si]) continue;
            for (int fi = 0; fi < F; ++fi) {
                int i = fi / nr + 1, j = fi % nr + 1;
                double z_rel = armature_.filament_axial_position(i) - arm_center_init;
                int idx = si * F + fi;
                pbuf_.seps[idx] = arm_center + z_rel - coils_[si].position();
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
}
```

- [ ] **Step 7: Add dispatcher compute_all_M1_dM1**

```cpp
template<typename SP>
void SimBatch<SP>::compute_all_M1_dM1() {
    if (use_persistent_)
        compute_all_M1_dM1_persistent();
    else
        compute_all_M1_dM1_fallback();
}
```

- [ ] **Step 8: Build and verify**

```sh
cmake --build build/ninja-cuda-debug --target coilgun_cuda
```
Expected: 0 errors.

- [ ] **Step 9: Commit**

```sh
git add src/cuda/sim_batch.cu
git commit -m "feat(cuda): implement persistent kernel path in SimBatch"
```

---

### Task 12: Test SimBatch persistent vs fallback consistency + full suite

**Files:**
- Modify: `tests/test_gpu_sim_batch.cpp`

- [ ] **Step 1: Read current test file**

Read `tests/test_gpu_sim_batch.cpp`.

- [ ] **Step 2: Add persistent vs fallback consistency test**

Add after the existing test cases:

```cpp
TEST_CASE("SimBatch — persistent vs fallback consistency") {
    auto coils = make_coils();
    auto arm   = make_arm();

    auto make_excs = [](double V) -> std::vector<std::unique_ptr<coilgun::simulation::Excitation>> {
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(V, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(V, 0.001));
        return excs;
    };
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    auto pol = coilgun::simulation::TerminationPolicy::defaults();
    pol.max_steps = 500;

    // Persistent (default)
    coilgun::simulation::cuda::GpuBackend be_p;
    be_p.use_persistent = true;
    auto coils_p = make_coils();
    auto arm_p   = make_arm();
    SimBatch<EulerStepper> batch_p(std::move(coils_p), arm_p, 1, 1e-6, be_p);
    batch_p.set_excitations(0, make_excs(450.0), triggers);
    batch_p.run(pol);
    double v_p = batch_p.result(0).summary.muzzle_velocity;

    // Fallback
    coilgun::simulation::cuda::GpuBackend be_f;
    be_f.use_persistent = false;
    auto coils_f = make_coils();
    auto arm_f   = make_arm();
    SimBatch<EulerStepper> batch_f(std::move(coils_f), arm_f, 1, 1e-6, be_f);
    batch_f.set_excitations(0, make_excs(450.0), triggers);
    batch_f.run(pol);
    double v_f = batch_f.result(0).summary.muzzle_velocity;

    CHECK(v_p > 0.0);
    CHECK(v_f > 0.0);
    CHECK(v_p == doctest::Approx(v_f).epsilon(1e-12));
}
```

- [ ] **Step 3: Build and run SimBatch tests**

```sh
cmake --build build/ninja-cuda-debug --target test_gpu_sim_batch
./build/ninja-cuda-debug/tests/test_gpu_sim_batch
```
Expected: all test cases pass (existing + new).

- [ ] **Step 4: Run with compute-sanitizer**

```sh
compute-sanitizer ./build/ninja-cuda-debug/tests/test_gpu_sim_batch
```
Expected: 0 errors.

- [ ] **Step 5: Commit**

```sh
git add tests/test_gpu_sim_batch.cpp
git commit -m "test(cuda): add persistent vs fallback consistency test for SimBatch"
```

---

### Task 13: Final verification — full test suite

**Files:** None (verification only)

- [ ] **Step 1: Full rebuild**

```sh
cmake --build build/ninja-cuda-debug
```
Expected: 0 errors.

- [ ] **Step 2: Run all GPU tests**

```sh
ctest --test-dir build/ninja-cuda-debug -R "test_gpu"
```
Expected: all pass.

- [ ] **Step 3: Run compute-sanitizer on all GPU tests**

```sh
for t in test_gpu_elliptic test_gpu_filament test_gpu_coil_pair \
         test_gpu_vs_cpu_single test_gpu_vs_cpu_multi \
         test_gpu_batch test_gpu_sim_batch; do
    compute-sanitizer ./build/ninja-cuda-debug/tests/$t || break
done
```
Expected: 0 errors for all.

- [ ] **Step 4: Verify callable with correct include**

Build a minimal example:
```cpp
#include <coilgun/coilgun_cuda.hpp>
int main() {
    coilgun::simulation::cuda::GpuBackend be;
    be.use_persistent = true;
    auto level = coilgun::simulation::cuda::GpuOptLevel::Aggressive;
    (void)level;
    return 0;
}
```

Compile with:
```sh
g++ -std=c++17 -Iinclude -Ibuild/ninja-cuda-debug/_deps/eigen-src \
    -o /tmp/quick_test /tmp/quick_test.cpp
```
Expected: compiles clean.

- [ ] **Step 5: Commit**

```sh
git commit -m "build: verify full GPU test suite passes with persistent kernel" --allow-empty
```

---

## Follow-up Plans

After this plan is complete:
1. **Double Buffering** — D5: CPU/GPU overlap via dual result buffers
2. **GpuOptLevel::Aggressive** — D6: FP32 integrand + hand-written AGM FP32 elliptic integrals

These are tracked as separate implementation plans per the design spec.
