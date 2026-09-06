#include "coilgun/simulation/cuda/gpu_state_kernels.hpp"
#include "gpu_kernel_launch_detail.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <limits>

namespace coilgun::simulation::cuda {
namespace {

constexpr std::size_t kParallelStateDimension = 32;

template <bool SingleStage>
__global__ void force_reduction_kernel(
    std::size_t S, std::size_t F, const double* currents, const double* dm1,
    const unsigned char* trigger, double* force) {
    extern __shared__ double partial[];
    const std::size_t batch = blockIdx.x;
    const unsigned int thread = threadIdx.x;
    const std::size_t terms = S * F;
    const double* batch_currents = currents + batch * (SingleStage ? 1 + F : S + F);
    const double* batch_dm1 = dm1 + batch * terms;

    double sum = 0.0;
    for (std::size_t term = thread; term < terms; term += blockDim.x) {
        const std::size_t stage = SingleStage ? 0 : term / F;
        const std::size_t filament = SingleStage ? term : term % F;
        if (trigger == nullptr || trigger[batch * S + stage] != 0) {
            sum += batch_currents[stage] * batch_currents[S + filament] * batch_dm1[term];
        }
    }
    partial[thread] = sum;
    __syncthreads();
    for (unsigned int stride = blockDim.x / 2; stride != 0; stride >>= 1) {
        if (thread < stride) partial[thread] += partial[thread + stride];
        __syncthreads();
    }
    if (thread == 0) force[batch] = partial[0];
}

__global__ void acceleration_kernel(
    std::size_t B, const double* force, double inverse_mass, double* acceleration) {
    const std::size_t batch = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch < B) acceleration[batch] = force[batch] * inverse_mass;
}

__global__ void state_update_kernel(
    std::size_t B, std::size_t D, const double* current_derivative, double dt,
    const double* acceleration, double* currents, double* velocity, double* position) {
    const std::size_t batch = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch >= B) return;
    const double old_velocity = velocity[batch];
    position[batch] += dt * old_velocity;
    velocity[batch] = old_velocity + dt * acceleration[batch];
    if (current_derivative != nullptr) {
        double* batch_currents = currents + batch * D;
        const double* batch_derivative = current_derivative + batch * D;
        for (std::size_t i = 0; i < D; ++i) batch_currents[i] += dt * batch_derivative[i];
    }
}

__global__ void masked_state_update_kernel(
    std::size_t B, std::size_t D, const double* current_derivative, double dt,
    const double* acceleration, const unsigned char* active, double* currents,
    double* velocity, double* position) {
    const std::size_t batch = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch >= B || (active != nullptr && active[batch] == 0)) return;
    const double old_velocity = velocity[batch];
    position[batch] += dt * old_velocity;
    velocity[batch] = old_velocity + dt * acceleration[batch];
    if (current_derivative != nullptr) {
        double* batch_currents = currents + batch * D;
        const double* batch_derivative = current_derivative + batch * D;
        for (std::size_t i = 0; i < D; ++i) batch_currents[i] += dt * batch_derivative[i];
    }
}

__global__ void parallel_state_update_kernel(
    std::size_t B, std::size_t D, const double* current_derivative, double dt,
    const double* acceleration, double* currents, double* velocity, double* position) {
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t total = B * D;
    if (index >= total) return;
    const std::size_t batch = index / D;
    const std::size_t dimension = index - batch * D;
    if (dimension == 0) {
        const double old_velocity = velocity[batch];
        position[batch] += dt * old_velocity;
        velocity[batch] = old_velocity + dt * acceleration[batch];
    }
    if (current_derivative != nullptr)
        currents[index] += dt * current_derivative[index];
}

__global__ void parallel_masked_state_update_kernel(
    std::size_t B, std::size_t D, const double* current_derivative, double dt,
    const double* acceleration, const unsigned char* active, double* currents,
    double* velocity, double* position) {
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t total = B * D;
    if (index >= total) return;
    const std::size_t batch = index / D;
    if (active[batch] == 0) return;
    const std::size_t dimension = index - batch * D;
    if (dimension == 0) {
        const double old_velocity = velocity[batch];
        position[batch] += dt * old_velocity;
        velocity[batch] = old_velocity + dt * acceleration[batch];
    }
    if (current_derivative != nullptr)
        currents[index] += dt * current_derivative[index];
}

__global__ void assembly_kernel(DeviceAssemblyView view) {
    const std::size_t batch = blockIdx.z;
    const std::size_t row = blockIdx.y * blockDim.y + threadIdx.y;
    const std::size_t column = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t S = view.stage_count;
    const std::size_t F = view.filament_count;
    const std::size_t D = S + F;
    if (batch >= view.batch_size || row >= D || column >= D) return;

    const std::size_t matrix_index = (batch * D + row) * D + column;
    const bool batch_active = view.active_mask[batch] != 0;
    double value = 0.0;
    if (!batch_active) {
        value = row == column ? 1.0 : 0.0;
    } else if (row < S) {
        const bool row_active = view.trigger_mask[batch * S + row] != 0 &&
            view.stage_mask[batch * S + row] != 0;
        if (!row_active) {
            value = row == column ? 1.0 : 0.0;
        } else if (column < S) {
            if (row == column) value = view.stage_inductances[row];
            else if (view.trigger_mask[batch * S + column] != 0 &&
                     view.stage_mask[batch * S + column] != 0)
                value = view.stage_mutual[row * S + column];
        } else {
            const std::size_t filament = column - S;
            if (view.mutual_stage_mask[batch * S + row] != 0)
                value = view.mutual[(batch * S + row) * F + filament];
        }
    } else {
        const std::size_t filament = row - S;
        if (column < S) {
            if (view.trigger_mask[batch * S + column] != 0 &&
                view.stage_mask[batch * S + column] != 0 &&
                view.mutual_stage_mask[batch * S + column] != 0)
                value = view.mutual[(batch * S + column) * F + filament];
        } else {
            const std::size_t other = column - S;
            value = filament == other
                ? view.filament_inductances[filament]
                : view.filament_mutual[filament * F + other];
        }
    }
    view.matrices[matrix_index] = value;

    if (column != 0) return;
    const std::size_t rhs_index = batch * D + row;
    if (!batch_active) {
        view.rhs[rhs_index] = 0.0;
        return;
    }
    const double speed = view.velocity[batch];
    if (row < S) {
        const bool row_active = view.trigger_mask[batch * S + row] != 0 &&
            view.stage_mask[batch * S + row] != 0;
        if (!row_active) {
            view.rhs[rhs_index] = 0.0;
            return;
        }
        double motional = 0.0;
        if (view.mutual_stage_mask[batch * S + row] != 0) {
            for (std::size_t filament = 0; filament < F; ++filament)
                motional += view.mutual_gradient[(batch * S + row) * F + filament] *
                    view.currents[batch * D + S + filament];
        }
        view.rhs[rhs_index] = view.stage_voltages[batch * S + row] -
            view.stage_resistances[row] * view.currents[batch * D + row] -
            speed * motional;
    } else {
        const std::size_t filament = row - S;
        const double resistance = view.dynamic_resistances != nullptr
            ? view.dynamic_resistances[batch * F + filament]
            : view.filament_reference_resistances[filament];
        double back_emf = 0.0;
        for (std::size_t stage = 0; stage < S; ++stage) {
            if (view.trigger_mask[batch * S + stage] != 0 &&
                view.stage_mask[batch * S + stage] != 0 &&
                view.mutual_stage_mask[batch * S + stage] != 0)
                back_emf += view.mutual_gradient[(batch * S + stage) * F + filament] *
                    view.currents[batch * D + stage];
        }
        view.rhs[rhs_index] = -resistance * view.currents[batch * D + row] -
            speed * back_emf;
    }
}

__global__ void separation_kernel(
    std::size_t B, std::size_t S, std::size_t F,
    const double* stage_positions, const double* filament_positions,
    const double* armature_positions, const unsigned char* active_mask,
    const unsigned char* trigger_mask, const unsigned char* mutual_stage_mask,
    unsigned char* pair_active, double* separations) {
    const std::size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t count = B * S * F;
    if (index >= count) return;
    const std::size_t filament = index % F;
    const std::size_t stage = (index / F) % S;
    const std::size_t batch = index / (S * F);
    const std::size_t stage_index = batch * S + stage;
    pair_active[index] = active_mask[batch] != 0 &&
        trigger_mask[stage_index] != 0 && mutual_stage_mask[stage_index] != 0;
    separations[index] = stage_positions[stage] -
        (filament_positions[filament] - armature_positions[batch]);
}

__global__ void compact_status_kernel(
    std::size_t B, std::size_t D, const double* currents,
    const double* velocity, const double* position, const double* residuals,
    const unsigned char* active_mask, DeviceStepStatus* status) {
    const std::size_t batch = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch >= B) return;
    bool finite = isfinite(velocity[batch]) && isfinite(position[batch]);
    for (std::size_t dimension = 0; finite && dimension < D; ++dimension)
        finite = isfinite(currents[batch * D + dimension]);
    const bool solver_ok = residuals == nullptr ||
        (isfinite(residuals[batch]) && residuals[batch] <= 1.0e-10);
    status[batch] = DeviceStepStatus{
        active_mask[batch], static_cast<unsigned char>(finite),
        static_cast<unsigned char>(solver_ok), 0};
}

__global__ void device_control_kernel(DeviceControlView view) {
    const std::size_t batch = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch >= view.batch_size || view.active_mask[batch] == 0) return;
    const std::size_t S = view.stage_count;
    const std::size_t D = view.dimension;
    const double armature_position = view.position_offsets[batch] - view.position[batch];

    for (std::size_t stage = 1; stage < S; ++stage) {
        const std::size_t index = batch * S + stage;
        if (view.trigger_mask[index] != 0 || view.stage_completed[index] != 0 ||
            view.trigger_mask[index - 1] == 0)
            continue;
        const unsigned char mode = view.trigger_modes[index];
        const double value = view.trigger_values[index];
        const bool triggered = mode == 1
            ? armature_position >= value
            : (mode == 2 && *view.current_time >= view.trigger_times[index - 1] + value);
        if (triggered) {
            view.trigger_mask[index] = 1;
            view.stage_mask[index] = 1;
            view.trigger_times[index] = *view.current_time;
            view.trigger_positions[index] = armature_position;
        }
    }

    bool any_remaining = false;
    for (std::size_t stage = 0; stage < S; ++stage) {
        const std::size_t index = batch * S + stage;
        if (view.stage_completed[index] == 0 && view.trigger_mask[index] != 0 &&
            view.excitation_finished[index] != 0 &&
            fabs(view.currents[batch * D + stage]) < view.quiet_current) {
            view.currents[batch * D + stage] = 0.0;
            view.stage_completed[index] = 1;
            view.stage_mask[index] = 0;
            view.mutual_stage_mask[index] = 0;
        }
        if (view.stage_completed[index] == 0) {
            const bool terminally_ineligible = view.trigger_mask[index] == 0 &&
                isinf(view.trigger_values[index]);
            any_remaining = any_remaining || !terminally_ineligible;
        }
    }
    if (!any_remaining) view.active_mask[batch] = 0;
}

bool valid_threads(unsigned int threads) {
    // force_reduction_kernel owns fixed 512-element shared reduction buffers.
    return threads != 0 && threads <= 512 && (threads & (threads - 1)) == 0;
}

bool checked_product(std::size_t a, std::size_t b, std::size_t& result) {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return false;
    result = a * b;
    return true;
}

enum class LaunchValidation { Checked, Trusted };

bool valid_grid(std::size_t blocks, LaunchValidation validation) {
    if (validation == LaunchValidation::Trusted) return true;
    int device = 0;
    cudaDeviceProp properties{};
    return cudaGetDevice(&device) == cudaSuccess &&
           cudaGetDeviceProperties(&properties, device) == cudaSuccess &&
           blocks <= static_cast<std::size_t>(properties.maxGridSize[0]);
}

bool device_pointer(const void* pointer, LaunchValidation validation) {
    if (!pointer) return false;
    if (validation == LaunchValidation::Trusted) return true;
    cudaPointerAttributes attributes{};
    const auto status = cudaPointerGetAttributes(&attributes, pointer);
#if CUDART_VERSION >= 10000
    return status == cudaSuccess && (attributes.type == cudaMemoryTypeDevice ||
                                     attributes.type == cudaMemoryTypeManaged);
#else
    return status == cudaSuccess && attributes.memoryType == cudaMemoryTypeDevice;
#endif
}

} // namespace

cudaError_t launch_force_reduction_impl(
    std::size_t B, std::size_t S, std::size_t F, const double* currents,
    const double* dm1, const unsigned char* trigger, double* force,
    StateKernelConfig config, cudaStream_t stream,
    LaunchValidation validation) noexcept {
    std::size_t terms = 0, state_size = 0;
    if (B == 0 || S == 0 || F == 0 || !device_pointer(currents, validation) ||
        !device_pointer(dm1, validation) ||
        (trigger != nullptr && !device_pointer(trigger, validation)) ||
        !device_pointer(force, validation) ||
        !valid_threads(config.threads_per_block) || B > std::numeric_limits<unsigned int>::max() ||
        !checked_product(S, F, terms) || !checked_product(B, terms, state_size) ||
        !valid_grid(B, validation)) {
        return cudaErrorInvalidValue;
    }
    const auto launch = [&](auto kernel) {
        kernel<<<static_cast<unsigned int>(B), config.threads_per_block,
                 config.threads_per_block * sizeof(double), stream>>>(
            S, F, currents, dm1, trigger, force);
        return cudaGetLastError();
    };
    if (S == 1) return launch(force_reduction_kernel<true>);
    return launch(force_reduction_kernel<false>);
}

cudaError_t launch_force_reduction(
    std::size_t B, std::size_t S, std::size_t F, const double* currents,
    const double* dm1, const unsigned char* trigger, double* force,
    StateKernelConfig config, cudaStream_t stream) noexcept {
    return launch_force_reduction_impl(B, S, F, currents, dm1, trigger, force,
                                       config, stream, LaunchValidation::Checked);
}

cudaError_t launch_device_assembly_impl(const DeviceAssemblyView& view,
                                        cudaStream_t stream,
                                        LaunchValidation validation) noexcept {
    if (view.batch_size == 0 || view.stage_count == 0 || view.filament_count == 0 ||
        !device_pointer(view.stage_inductances, validation) ||
        !device_pointer(view.stage_resistances, validation) ||
        !device_pointer(view.stage_mutual, validation) ||
        !device_pointer(view.filament_inductances, validation) ||
        !device_pointer(view.filament_reference_resistances, validation) ||
        !device_pointer(view.filament_mutual, validation) ||
        !device_pointer(view.mutual, validation) ||
        !device_pointer(view.mutual_gradient, validation) ||
        !device_pointer(view.currents, validation) ||
        !device_pointer(view.velocity, validation) ||
        !device_pointer(view.stage_voltages, validation) ||
        !device_pointer(view.active_mask, validation) ||
        !device_pointer(view.trigger_mask, validation) ||
        !device_pointer(view.stage_mask, validation) ||
        !device_pointer(view.mutual_stage_mask, validation) ||
        !device_pointer(view.matrices, validation) || !device_pointer(view.rhs, validation) ||
        (view.dynamic_resistances != nullptr &&
         !device_pointer(view.dynamic_resistances, validation)))
        return cudaErrorInvalidValue;
    const std::size_t dimension = view.stage_count + view.filament_count;
    if (dimension > std::numeric_limits<unsigned int>::max() ||
        view.batch_size > std::numeric_limits<unsigned int>::max())
        return cudaErrorInvalidValue;
    const dim3 block(16, 16, 1);
    const dim3 grid(static_cast<unsigned int>((dimension + block.x - 1) / block.x),
                    static_cast<unsigned int>((dimension + block.y - 1) / block.y),
                    static_cast<unsigned int>(view.batch_size));
    assembly_kernel<<<grid, block, 0, stream>>>(view);
    return cudaGetLastError();
}

cudaError_t launch_device_assembly(const DeviceAssemblyView& view,
                                   cudaStream_t stream) noexcept {
    return launch_device_assembly_impl(view, stream, LaunchValidation::Checked);
}

cudaError_t launch_mutual_input_update_impl(
    std::size_t B, std::size_t S, std::size_t F,
    const double* stage_positions, const double* filament_positions,
    const double* armature_positions, const unsigned char* active_mask,
    const unsigned char* trigger_mask, const unsigned char* mutual_stage_mask,
    unsigned char* pair_active, double* separations,
    cudaStream_t stream, LaunchValidation validation) noexcept {
    std::size_t count = 0;
    if (B == 0 || S == 0 || F == 0 || !checked_product(B, S, count) ||
        !checked_product(count, F, count) ||
        !device_pointer(stage_positions, validation) ||
        !device_pointer(filament_positions, validation) ||
        !device_pointer(armature_positions, validation) ||
        !device_pointer(active_mask, validation) ||
        !device_pointer(trigger_mask, validation) ||
        !device_pointer(mutual_stage_mask, validation) ||
        !device_pointer(pair_active, validation) ||
        !device_pointer(separations, validation) ||
        count > std::numeric_limits<std::size_t>::max() - 255 ||
        !valid_grid((count + 255) / 256, validation))
        return cudaErrorInvalidValue;
    separation_kernel<<<static_cast<unsigned int>((count + 255) / 256), 256, 0, stream>>>(
        B, S, F, stage_positions, filament_positions, armature_positions,
        active_mask, trigger_mask, mutual_stage_mask, pair_active, separations);
    return cudaGetLastError();
}

cudaError_t launch_mutual_input_update(
    std::size_t B, std::size_t S, std::size_t F,
    const double* stage_positions, const double* filament_positions,
    const double* armature_positions, const unsigned char* active_mask,
    const unsigned char* trigger_mask, const unsigned char* mutual_stage_mask,
    unsigned char* pair_active, double* separations,
    cudaStream_t stream) noexcept {
    return launch_mutual_input_update_impl(
        B, S, F, stage_positions, filament_positions, armature_positions,
        active_mask, trigger_mask, mutual_stage_mask, pair_active, separations,
        stream, LaunchValidation::Checked);
}

cudaError_t launch_compact_status_impl(
    std::size_t B, std::size_t D, const double* currents,
    const double* velocity, const double* position, const double* residuals,
    const unsigned char* active_mask, DeviceStepStatus* status,
    cudaStream_t stream, LaunchValidation validation) noexcept {
    if (B == 0 || D == 0 || !device_pointer(currents, validation) ||
        !device_pointer(velocity, validation) || !device_pointer(position, validation) ||
        !device_pointer(active_mask, validation) || !device_pointer(status, validation) ||
        (residuals != nullptr && !device_pointer(residuals, validation)) ||
        B > std::numeric_limits<std::size_t>::max() - 127 ||
        !valid_grid((B + 127) / 128, validation))
        return cudaErrorInvalidValue;
    compact_status_kernel<<<static_cast<unsigned int>((B + 127) / 128), 128, 0, stream>>>(
        B, D, currents, velocity, position, residuals, active_mask, status);
    return cudaGetLastError();
}

cudaError_t launch_compact_status(
    std::size_t B, std::size_t D, const double* currents,
    const double* velocity, const double* position, const double* residuals,
    const unsigned char* active_mask, DeviceStepStatus* status,
    cudaStream_t stream) noexcept {
    return launch_compact_status_impl(B, D, currents, velocity, position, residuals,
                                      active_mask, status, stream,
                                      LaunchValidation::Checked);
}

cudaError_t launch_device_control_impl(const DeviceControlView& view,
                                       cudaStream_t stream,
                                       LaunchValidation validation) noexcept {
    if (view.batch_size == 0 || view.stage_count == 0 ||
        view.dimension < view.stage_count || !device_pointer(view.current_time, validation) ||
        !std::isfinite(view.quiet_current) || view.quiet_current < 0.0 ||
        !device_pointer(view.currents, validation) ||
        !device_pointer(view.position, validation) ||
        !device_pointer(view.position_offsets, validation) ||
        !device_pointer(view.trigger_modes, validation) ||
        !device_pointer(view.trigger_values, validation) ||
        !device_pointer(view.excitation_finished, validation) ||
        !device_pointer(view.active_mask, validation) ||
        !device_pointer(view.trigger_mask, validation) ||
        !device_pointer(view.stage_mask, validation) ||
        !device_pointer(view.mutual_stage_mask, validation) ||
        !device_pointer(view.stage_completed, validation) ||
        !device_pointer(view.trigger_times, validation) ||
        !device_pointer(view.trigger_positions, validation) ||
        view.batch_size > std::numeric_limits<std::size_t>::max() - 127 ||
        !valid_grid((view.batch_size + 127) / 128, validation))
        return cudaErrorInvalidValue;
    device_control_kernel<<<static_cast<unsigned int>((view.batch_size + 127) / 128),
                            128, 0, stream>>>(view);
    return cudaGetLastError();
}

cudaError_t launch_device_control(const DeviceControlView& view,
                                  cudaStream_t stream) noexcept {
    return launch_device_control_impl(view, stream, LaunchValidation::Checked);
}

cudaError_t launch_acceleration_impl(
    std::size_t B, const double* force, double mass, double* acceleration,
    cudaStream_t stream, LaunchValidation validation) noexcept {
    if (B == 0 || B > std::numeric_limits<unsigned int>::max() ||
        !device_pointer(force, validation) || !device_pointer(acceleration, validation) ||
        !std::isfinite(mass) || mass <= 0.0 ||
        B > std::numeric_limits<std::size_t>::max() - 255 ||
        !valid_grid((B + 255) / 256, validation)) {
        return cudaErrorInvalidValue;
    }
    acceleration_kernel<<<static_cast<unsigned int>((B + 255) / 256), 256, 0, stream>>>(
        B, force, 1.0 / mass, acceleration);
    return cudaGetLastError();
}

cudaError_t launch_acceleration(
    std::size_t B, const double* force, double mass, double* acceleration,
    cudaStream_t stream) noexcept {
    return launch_acceleration_impl(B, force, mass, acceleration, stream,
                                    LaunchValidation::Checked);
}

cudaError_t launch_state_update(
    std::size_t B, std::size_t S, std::size_t F, double* currents,
    const double* current_derivative, const double* dm1, const unsigned char* trigger,
    double mass, double dt, double* acceleration, double* velocity, double* position,
    double* force, StateKernelConfig config, cudaStream_t stream) noexcept {
    std::size_t dimension = 0, state_count = 0;
    if (B == 0 || S == 0 || F == 0 ||
        !device_pointer(currents, LaunchValidation::Checked) ||
        !device_pointer(dm1, LaunchValidation::Checked) ||
        (trigger != nullptr && !device_pointer(trigger, LaunchValidation::Checked)) ||
        !device_pointer(acceleration, LaunchValidation::Checked) ||
        !device_pointer(velocity, LaunchValidation::Checked) ||
        !device_pointer(position, LaunchValidation::Checked) ||
        !device_pointer(force, LaunchValidation::Checked) ||
        !std::isfinite(dt) || dt <= 0.0 || !std::isfinite(mass) || mass <= 0.0 ||
        B > std::numeric_limits<unsigned int>::max() || S > std::numeric_limits<std::size_t>::max() - F) {
        return cudaErrorInvalidValue;
    }
    dimension = S + F;
    cudaError_t status = launch_force_reduction(B, S, F, currents, dm1, trigger, force, config, stream);
    if (status != cudaSuccess) return status;
    status = launch_acceleration(B, force, mass, acceleration, stream);
    if (status != cudaSuccess) return status;
    if (dimension <= kParallelStateDimension) {
        if (B > std::numeric_limits<std::size_t>::max() - 255 ||
            !valid_grid((B + 255) / 256, LaunchValidation::Checked))
            return cudaErrorInvalidValue;
        state_update_kernel<<<static_cast<unsigned int>((B + 255) / 256), 256, 0, stream>>>(
            B, dimension, current_derivative, dt, acceleration, currents, velocity, position);
    } else {
        if (!checked_product(B, dimension, state_count) ||
            state_count > std::numeric_limits<std::size_t>::max() - 255 ||
            !valid_grid((state_count + 255) / 256, LaunchValidation::Checked))
            return cudaErrorInvalidValue;
        parallel_state_update_kernel<<<static_cast<unsigned int>((state_count + 255) / 256), 256, 0, stream>>>(
            B, dimension, current_derivative, dt, acceleration, currents, velocity, position);
    }
    return cudaGetLastError();
}

cudaError_t launch_state_update_masked_impl(
    std::size_t B, std::size_t S, std::size_t F, double* currents,
    const double* current_derivative, const double* dm1, const unsigned char* trigger,
    const unsigned char* active, double mass, double dt, double* acceleration,
    double* velocity, double* position, double* force, StateKernelConfig config,
    cudaStream_t stream, LaunchValidation validation) noexcept {
    std::size_t terms = 0, dimension = 0, state_count = 0;
    if (!device_pointer(active, validation) || B == 0 || S == 0 || F == 0 ||
        !device_pointer(currents, validation) || !device_pointer(dm1, validation) ||
        (trigger != nullptr && !device_pointer(trigger, validation)) ||
        !device_pointer(acceleration, validation) ||
        !device_pointer(velocity, validation) || !device_pointer(position, validation) ||
        !device_pointer(force, validation) ||
        !std::isfinite(dt) || dt <= 0.0 || !std::isfinite(mass) || mass <= 0.0 ||
        B > std::numeric_limits<unsigned int>::max() || S > std::numeric_limits<std::size_t>::max() - F ||
        !checked_product(S, F, terms) || !checked_product(B, terms, terms) ||
        !valid_grid((B + 255) / 256, validation))
        return cudaErrorInvalidValue;
    dimension = S + F;
    auto status = launch_force_reduction_impl(B, S, F, currents, dm1, trigger, force,
                                              config, stream, validation);
    if (status != cudaSuccess) return status;
    status = launch_acceleration_impl(B, force, mass, acceleration, stream, validation);
    if (status != cudaSuccess) return status;
    if (dimension <= kParallelStateDimension) {
        masked_state_update_kernel<<<static_cast<unsigned int>((B + 255) / 256), 256, 0, stream>>>(
            B, dimension, current_derivative, dt, acceleration, active, currents, velocity, position);
    } else {
        if (!checked_product(B, dimension, state_count) ||
            state_count > std::numeric_limits<std::size_t>::max() - 255 ||
            !valid_grid((state_count + 255) / 256, validation))
            return cudaErrorInvalidValue;
        parallel_masked_state_update_kernel<<<static_cast<unsigned int>((state_count + 255) / 256), 256, 0, stream>>>(
            B, dimension, current_derivative, dt, acceleration, active, currents, velocity, position);
    }
    return cudaGetLastError();
}

cudaError_t launch_state_update_masked(
    std::size_t B, std::size_t S, std::size_t F, double* currents,
    const double* current_derivative, const double* dm1, const unsigned char* trigger,
    const unsigned char* active, double mass, double dt, double* acceleration,
    double* velocity, double* position, double* force, StateKernelConfig config,
    cudaStream_t stream) noexcept {
    return launch_state_update_masked_impl(
        B, S, F, currents, current_derivative, dm1, trigger, active, mass, dt,
        acceleration, velocity, position, force, config, stream,
        LaunchValidation::Checked);
}

namespace detail {

cudaError_t launch_device_assembly_unchecked(const DeviceAssemblyView& view,
                                             cudaStream_t stream) noexcept {
    return launch_device_assembly_impl(view, stream, LaunchValidation::Trusted);
}

cudaError_t launch_mutual_input_update_unchecked(
    std::size_t B, std::size_t S, std::size_t F,
    const double* stage_positions, const double* filament_positions,
    const double* armature_positions, const unsigned char* active_mask,
    const unsigned char* trigger_mask, const unsigned char* mutual_stage_mask,
    unsigned char* pair_active, double* separations,
    cudaStream_t stream) noexcept {
    return launch_mutual_input_update_impl(
        B, S, F, stage_positions, filament_positions, armature_positions,
        active_mask, trigger_mask, mutual_stage_mask, pair_active, separations,
        stream, LaunchValidation::Trusted);
}

cudaError_t launch_compact_status_unchecked(
    std::size_t B, std::size_t D, const double* currents,
    const double* velocity, const double* position, const double* residuals,
    const unsigned char* active_mask, DeviceStepStatus* status,
    cudaStream_t stream) noexcept {
    return launch_compact_status_impl(B, D, currents, velocity, position, residuals,
                                      active_mask, status, stream,
                                      LaunchValidation::Trusted);
}

cudaError_t launch_device_control_unchecked(const DeviceControlView& view,
                                            cudaStream_t stream) noexcept {
    return launch_device_control_impl(view, stream, LaunchValidation::Trusted);
}

cudaError_t launch_state_update_masked_unchecked(
    std::size_t B, std::size_t S, std::size_t F, double* currents,
    const double* current_derivative, const double* dm1,
    const unsigned char* trigger, const unsigned char* active, double mass,
    double dt, double* acceleration, double* velocity, double* position,
    double* force, StateKernelConfig config, cudaStream_t stream) noexcept {
    return launch_state_update_masked_impl(
        B, S, F, currents, current_derivative, dm1, trigger, active, mass, dt,
        acceleration, velocity, position, force, config, stream,
        LaunchValidation::Trusted);
}

} // namespace detail

} // namespace coilgun::simulation::cuda
