/**
 * @file gpu_backend.hpp
 * @brief GPU backend types and configuration.
 * @author Winston Meursault
 */

#pragma once

#include <cstddef>

namespace coilgun::simulation::cuda {

/**
 * @brief GPU optimization level.
 *
 * Controls distance cutoff and adaptive quadrature order.
 */
enum class GpuOptLevel {
    Standard = 0,   ///< Fixed n_nodes=9, no distance cutoff.
    Full     = 1,   ///< Distance cutoff (>10× coil length) + fixed n_nodes=9.
};

/**
 * @brief GPU backend configuration.
 */
struct GpuBackend {
    int     device_id         = 0;    ///< cudaSetDevice target.
    int     threads_per_block = 512;  ///< Threads per block for integration kernel.
    size_t  max_batch_sims    = 256;  ///< Pre-allocated buffer size for batch mode.
    bool    enable_profiling  = false; ///< Enable NVTX range annotations.
};

} // namespace coilgun::simulation::cuda
