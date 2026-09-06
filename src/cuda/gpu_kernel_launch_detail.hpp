#pragma once

#include "coilgun/simulation/cuda/gpu_mutual_pipeline.hpp"
#include "coilgun/simulation/cuda/gpu_state_kernels.hpp"

namespace coilgun::simulation::cuda::detail {

// These entry points are used only after GpuEngine has allocated and fixed all
// device buffers. Public launch functions retain full pointer and device-limit
// validation for borrowed external views.
cudaError_t launch_device_assembly_unchecked(const DeviceAssemblyView& view,
                                             cudaStream_t stream) noexcept;

cudaError_t launch_mutual_input_update_unchecked(
    std::size_t batch_size, std::size_t stage_count, std::size_t filament_count,
    const double* stage_positions, const double* filament_positions,
    const double* armature_positions, const unsigned char* active_mask,
    const unsigned char* trigger_mask, const unsigned char* mutual_stage_mask,
    unsigned char* pair_active, double* separations,
    cudaStream_t stream) noexcept;

cudaError_t launch_compact_status_unchecked(
    std::size_t batch_size, std::size_t dimension, const double* currents,
    const double* velocity, const double* position, const double* residuals,
    const unsigned char* active_mask, DeviceStepStatus* status,
    cudaStream_t stream) noexcept;

cudaError_t launch_device_control_unchecked(const DeviceControlView& view,
                                            cudaStream_t stream) noexcept;

cudaError_t launch_state_update_masked_unchecked(
    std::size_t batch_size, std::size_t stage_count, std::size_t filament_count,
    double* currents, const double* current_derivative, const double* dm1,
    const unsigned char* trigger_mask, const unsigned char* active_mask,
    double mass, double dt, double* acceleration, double* velocity,
    double* position, double* force, StateKernelConfig config,
    cudaStream_t stream) noexcept;

void launch_mutual_pipeline_unchecked(const MutualPipelineView& view,
                                      GpuOptLevel opt_level,
                                      int threads_per_block,
                                      cudaStream_t stream);

} // namespace coilgun::simulation::cuda::detail
