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

struct DeviceAssemblyView {
    std::size_t batch_size = 0;
    std::size_t stage_count = 0;
    std::size_t filament_count = 0;
    const double* stage_inductances = nullptr;
    const double* stage_resistances = nullptr;
    const double* stage_mutual = nullptr;
    const double* filament_inductances = nullptr;
    const double* filament_reference_resistances = nullptr;
    const double* filament_mutual = nullptr;
    const double* dynamic_resistances = nullptr;
    const double* mutual = nullptr;
    const double* mutual_gradient = nullptr;
    const double* currents = nullptr;
    const double* velocity = nullptr;
    const double* stage_voltages = nullptr;
    const unsigned char* active_mask = nullptr;
    const unsigned char* trigger_mask = nullptr;
    const unsigned char* stage_mask = nullptr;
    const unsigned char* mutual_stage_mask = nullptr;
    double* matrices = nullptr;
    double* rhs = nullptr;
};

struct DeviceStepStatus {
    unsigned char active = 0;
    unsigned char finite = 0;
    unsigned char solver_ok = 0;
    unsigned char reserved = 0;
};

struct DeviceControlView {
    std::size_t batch_size = 0;
    std::size_t stage_count = 0;
    std::size_t filament_count = 0;
    std::size_t dimension = 0;
    double quiet_current = 1.0e-6;
    const double* current_time = nullptr;
    double* currents = nullptr;
    const double* position = nullptr;
    const double* position_offsets = nullptr;
    const unsigned char* trigger_modes = nullptr;
    const double* trigger_values = nullptr;
    const unsigned char* excitation_finished = nullptr;
    unsigned char* active_mask = nullptr;
    unsigned char* trigger_mask = nullptr;
    unsigned char* stage_mask = nullptr;
    unsigned char* mutual_stage_mask = nullptr;
    unsigned char* stage_completed = nullptr;
    unsigned char* pair_active = nullptr;
    double* trigger_times = nullptr;
    double* trigger_positions = nullptr;
};

cudaError_t launch_device_assembly(const DeviceAssemblyView& view,
                                   cudaStream_t stream = nullptr) noexcept;

cudaError_t launch_separation_update(
    std::size_t batch_size, std::size_t stage_count, std::size_t filament_count,
    const double* stage_positions, const double* filament_positions,
    const double* armature_positions, double* separations,
    cudaStream_t stream = nullptr) noexcept;

cudaError_t launch_compact_status(
    std::size_t batch_size, std::size_t dimension,
    const double* currents, const double* velocity, const double* position,
    const double* residuals, const unsigned char* active_mask,
    DeviceStepStatus* status, cudaStream_t stream = nullptr) noexcept;

cudaError_t launch_device_control(const DeviceControlView& view,
                                  cudaStream_t stream = nullptr) noexcept;

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
