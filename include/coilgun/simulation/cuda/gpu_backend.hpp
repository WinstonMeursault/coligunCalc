/**
 * @file gpu_backend.hpp
 * @brief GPU backend types and configuration.
 * @author Winston Meursault
 */

#pragma once

#include "coilgun/simulation/cuda/gpu_execution_config.hpp"

#include <cstddef>
#include <stdexcept>

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
    bool    enable_profiling  = false; ///< Retain profiling-request metadata; host-wall timing fields are always collected. No NVTX guarantee.
    // Kept as a request flag for source compatibility. It is consulted only
    // when backend is Auto; explicit backend modes take precedence.
    bool    use_persistent    = true;
    // Graph is the default optimized backend. When selected and supported, it
    // captures the complete fixed-shape resident physical step.
    BackendMode backend = BackendMode::Graph;

    void validate() const {
        switch (backend) {
        case BackendMode::Auto:
        case BackendMode::Graph:
        case BackendMode::Persistent:
        case BackendMode::Fallback:
        case BackendMode::Direct:
            break;
        default:
            throw std::invalid_argument("GPU backend mode is invalid");
        }
        if (device_id < 0) throw std::invalid_argument("GPU device_id must not be negative");
        if (threads_per_block <= 0 || threads_per_block > 512 ||
            (threads_per_block & (threads_per_block - 1)) != 0) {
            throw std::invalid_argument("GPU threads_per_block must be a positive power of two no greater than 512");
        }
        if (max_batch_sims == 0) {
            throw std::invalid_argument("GPU max_batch_sims must be positive");
        }
    }
};

} // namespace coilgun::simulation::cuda
