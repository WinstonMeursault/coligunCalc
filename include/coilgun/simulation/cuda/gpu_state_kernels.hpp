/**
 * @file gpu_state_kernels.hpp
 * @brief Batched force and explicit state-update CUDA kernels.
 */

#pragma once

#include <cuda_runtime_api.h>

#include <cstddef>

namespace coilgun::simulation::cuda {

/** Kernel launch controls. The fixed reduction tree is used when deterministic is true. */
struct StateKernelConfig {
    bool deterministic = true;
    unsigned int threads_per_block = 128;
};

cudaError_t launch_force_reduction(
    std::size_t batch_size, std::size_t stage_count, std::size_t filament_count,
    const double* currents, const double* dm1, const unsigned char* trigger_mask,
    double* force, StateKernelConfig config = {}, cudaStream_t stream = nullptr) noexcept;

cudaError_t launch_acceleration(
    std::size_t batch_size, const double* force, double mass, double* acceleration,
    cudaStream_t stream = nullptr) noexcept;

/**
 * Apply one explicit Euler step. Position uses the old velocity, while
 * velocity uses the current-state acceleration, matching EulerStepper.
 */
cudaError_t launch_state_update(
    std::size_t batch_size, std::size_t stage_count, std::size_t filament_count,
    double* currents, const double* current_derivative, const double* dm1,
    const unsigned char* trigger_mask, double mass, double dt, double* acceleration,
    double* velocity, double* position, double* force, StateKernelConfig config = {},
    cudaStream_t stream = nullptr) noexcept;

cudaError_t launch_state_update_masked(
    std::size_t batch_size, std::size_t stage_count, std::size_t filament_count,
    double* currents, const double* current_derivative, const double* dm1,
    const unsigned char* trigger_mask, const unsigned char* active_mask,
    double mass, double dt, double* acceleration,
    double* velocity, double* position, double* force, StateKernelConfig config = {},
    cudaStream_t stream = nullptr) noexcept;

} // namespace coilgun::simulation::cuda
