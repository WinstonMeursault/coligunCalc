#include "coilgun/simulation/cuda/gpu_state_kernels.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <limits>

namespace coilgun::simulation::cuda {
namespace {

__global__ void force_reduction_kernel(
    std::size_t S, std::size_t F, const double* currents, const double* dm1,
    const unsigned char* trigger, double* force) {
    extern __shared__ double partial[];
    const std::size_t batch = blockIdx.x;
    const unsigned int thread = threadIdx.x;
    const std::size_t terms = S * F;
    const double* batch_currents = currents + batch * (S + F);
    const double* batch_dm1 = dm1 + batch * terms;

    double sum = 0.0;
    for (std::size_t term = thread; term < terms; term += blockDim.x) {
        const std::size_t stage = term / F;
        if (trigger == nullptr || trigger[batch * S + stage] != 0) {
            sum += batch_currents[stage] * batch_currents[S + term % F] * batch_dm1[term];
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

bool valid_threads(unsigned int threads) {
    return threads != 0 && threads <= 1024 && (threads & (threads - 1)) == 0;
}

bool checked_product(std::size_t a, std::size_t b, std::size_t& result) {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return false;
    result = a * b;
    return true;
}

bool valid_grid(std::size_t blocks) {
    int device = 0;
    cudaDeviceProp properties{};
    return cudaGetDevice(&device) == cudaSuccess &&
           cudaGetDeviceProperties(&properties, device) == cudaSuccess &&
           blocks <= static_cast<std::size_t>(properties.maxGridSize[0]);
}

bool device_pointer(const void* pointer) {
    if (!pointer) return false;
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

cudaError_t launch_force_reduction(
    std::size_t B, std::size_t S, std::size_t F, const double* currents,
    const double* dm1, const unsigned char* trigger, double* force,
    StateKernelConfig config, cudaStream_t stream) noexcept {
    std::size_t terms = 0, state_size = 0;
    if (B == 0 || S == 0 || F == 0 || !device_pointer(currents) || !device_pointer(dm1) ||
        (trigger != nullptr && !device_pointer(trigger)) || !device_pointer(force) ||
        !valid_threads(config.threads_per_block) || B > std::numeric_limits<unsigned int>::max() ||
        !checked_product(S, F, terms) || !checked_product(B, terms, state_size) || !valid_grid(B)) {
        return cudaErrorInvalidValue;
    }
    force_reduction_kernel<<<static_cast<unsigned int>(B), config.threads_per_block,
                             config.threads_per_block * sizeof(double), stream>>>(
        S, F, currents, dm1, trigger, force);
    return cudaGetLastError();
}

cudaError_t launch_acceleration(
    std::size_t B, const double* force, double mass, double* acceleration,
    cudaStream_t stream) noexcept {
    if (B == 0 || B > std::numeric_limits<unsigned int>::max() || !device_pointer(force) ||
        !device_pointer(acceleration) || !std::isfinite(mass) || mass <= 0.0 ||
        B > std::numeric_limits<std::size_t>::max() - 255 || !valid_grid((B + 255) / 256)) {
        return cudaErrorInvalidValue;
    }
    acceleration_kernel<<<static_cast<unsigned int>((B + 255) / 256), 256, 0, stream>>>(
        B, force, 1.0 / mass, acceleration);
    return cudaGetLastError();
}

cudaError_t launch_state_update(
    std::size_t B, std::size_t S, std::size_t F, double* currents,
    const double* current_derivative, const double* dm1, const unsigned char* trigger,
    double mass, double dt, double* acceleration, double* velocity, double* position,
    double* force, StateKernelConfig config, cudaStream_t stream) noexcept {
    if (B == 0 || S == 0 || F == 0 || !device_pointer(currents) || !device_pointer(dm1) ||
        (trigger != nullptr && !device_pointer(trigger)) ||
        !device_pointer(acceleration) || !device_pointer(velocity) || !device_pointer(position) || !device_pointer(force) ||
        !std::isfinite(dt) || dt <= 0.0 || !std::isfinite(mass) || mass <= 0.0 ||
        B > std::numeric_limits<unsigned int>::max()) {
        return cudaErrorInvalidValue;
    }
    cudaError_t status = launch_force_reduction(B, S, F, currents, dm1, trigger, force, config, stream);
    if (status != cudaSuccess) return status;
    status = launch_acceleration(B, force, mass, acceleration, stream);
    if (status != cudaSuccess) return status;
    state_update_kernel<<<static_cast<unsigned int>((B + 255) / 256), 256, 0, stream>>>(
        B, S + F, current_derivative, dt, acceleration, currents, velocity, position);
    return cudaGetLastError();
}

cudaError_t launch_state_update_masked(
    std::size_t B, std::size_t S, std::size_t F, double* currents,
    const double* current_derivative, const double* dm1, const unsigned char* trigger,
    const unsigned char* active, double mass, double dt, double* acceleration,
    double* velocity, double* position, double* force, StateKernelConfig config,
    cudaStream_t stream) noexcept {
    std::size_t terms = 0;
    if (!device_pointer(active) || B == 0 || S == 0 || F == 0 || !device_pointer(currents) || !device_pointer(dm1) ||
        (trigger != nullptr && !device_pointer(trigger)) ||
        !device_pointer(acceleration) || !device_pointer(velocity) || !device_pointer(position) || !device_pointer(force) ||
        !std::isfinite(dt) || dt <= 0.0 || !std::isfinite(mass) || mass <= 0.0 ||
        B > std::numeric_limits<unsigned int>::max() ||
        !checked_product(S, F, terms) || !checked_product(B, terms, terms) ||
        !valid_grid((B + 255) / 256))
        return cudaErrorInvalidValue;
    auto status = launch_force_reduction(B, S, F, currents, dm1, trigger, force, config, stream);
    if (status != cudaSuccess) return status;
    status = launch_acceleration(B, force, mass, acceleration, stream);
    if (status != cudaSuccess) return status;
    masked_state_update_kernel<<<static_cast<unsigned int>((B + 255) / 256), 256, 0, stream>>>(
        B, S + F, current_derivative, dt, acceleration, active, currents, velocity, position);
    return cudaGetLastError();
}

} // namespace coilgun::simulation::cuda
