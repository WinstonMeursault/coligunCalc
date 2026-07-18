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
