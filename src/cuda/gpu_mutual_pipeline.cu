#include "coilgun/simulation/cuda/gpu_mutual_pipeline.hpp"
#include "coilgun/physics/mutual_inductance.cuh"
#include "coilgun/physics/quadrature.hpp"
#include <cuda_runtime.h>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace {

bool checked_product(std::size_t a, std::size_t b, std::size_t& result) {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return false;
    result = a * b;
    return true;
}

__constant__ double pipeline_gl_nodes[9];
__constant__ double pipeline_gl_weights[9];

__device__ inline std::size_t pipeline_index(int simulation, int stage,
                                              int filament, int stages,
                                              int filaments) {
    return (static_cast<std::size_t>(simulation) * stages + stage) * filaments
           + filament;
}

template <bool Aggressive, bool UseCutoff>
__global__ void mutual_pipeline_kernel(
        const coilgun::simulation::cuda::CoilGeo* coils,
        const coilgun::simulation::cuda::FilGeo* filaments,
        const double* separations, const std::uint8_t* active_mask,
        double* mutual, double* gradient, int batch_size, int stage_count,
        int filament_count, int n_nodes) {
    const int simulation = static_cast<int>(blockIdx.z);
    const int stage = static_cast<int>(blockIdx.y);
    const int filament = static_cast<int>(blockIdx.x);
    if (simulation >= batch_size || stage >= stage_count || filament >= filament_count)
        return;

    const std::size_t out = pipeline_index(simulation, stage, filament,
                                            stage_count, filament_count);
    if (active_mask[out] == 0) {
        if (threadIdx.x == 0) {
            mutual[out] = 0.0;
            gradient[out] = 0.0;
        }
        return;
    }

    const auto coil = coils[stage];
    const auto fil = filaments[filament];
    const double separation = separations[out];
    const double cutoff = 10.0 * (fabs(coil.len) + fabs(fil.len)
                                  + fabs(coil.re - coil.ri)
                                  + fabs(fil.re - fil.ri));
    if constexpr (UseCutoff) {
        if (fabs(separation) > cutoff) {
            if (threadIdx.x == 0) {
                mutual[out] = 0.0;
                gradient[out] = 0.0;
            }
            return;
        }
    }

    const double ra_mid = 0.5 * (coil.re + coil.ri);
    const double ra_half = 0.5 * (coil.re - coil.ri);
    const double rb_mid = 0.5 * (fil.re + fil.ri);
    const double rb_half = 0.5 * (fil.re - fil.ri);
    const double la_half = 0.5 * coil.len;
    const double lb_half = 0.5 * fil.len;
    const double prefactor = coil.turns / 16.0;
    const int n4 = n_nodes * n_nodes * n_nodes * n_nodes;
    const int tid = static_cast<int>(threadIdx.x);
    const int threads = static_cast<int>(blockDim.x);
    __shared__ double sum_m[512];
    __shared__ double sum_gradient[512];
    sum_m[tid] = 0.0;
    sum_gradient[tid] = 0.0;

    for (int point = tid; point < n4; point += threads) {
        const int i1 = point / (n_nodes * n_nodes * n_nodes);
        const int rem1 = point % (n_nodes * n_nodes * n_nodes);
        const int j1 = rem1 / (n_nodes * n_nodes);
        const int rem2 = rem1 % (n_nodes * n_nodes);
        const int i2 = rem2 / n_nodes;
        const int j2 = rem2 % n_nodes;
        const double weight = pipeline_gl_weights[i1] * pipeline_gl_weights[j1]
                            * pipeline_gl_weights[i2] * pipeline_gl_weights[j2];
        const double ra = ra_mid + ra_half * pipeline_gl_nodes[i1];
        const double rb = rb_mid + rb_half * pipeline_gl_nodes[i2];
        const double za = la_half * pipeline_gl_nodes[j1];
        const double zb = separation + lb_half * pipeline_gl_nodes[j2];
        if constexpr (Aggressive) {
            const float fra = static_cast<float>(ra);
            const float frb = static_cast<float>(rb);
            const float fdz = static_cast<float>(zb - za);
            const float fw = static_cast<float>(weight);
            sum_m[tid] += static_cast<double>(fw *
                coilgun::physics::mutual_inductance_filament_f32(
                    fra, frb, fabsf(fdz)));
            sum_gradient[tid] += static_cast<double>(fw *
                coilgun::physics::mutual_inductance_gradient_filament_f32(
                    fra, frb, fdz));
        } else {
            sum_m[tid] += weight * coilgun::physics::mutual_inductance_filament_device(
                ra, rb, fabs(zb - za));
            sum_gradient[tid] += weight *
                coilgun::physics::mutual_inductance_gradient_filament_device(
                    ra, rb, zb - za);
        }
    }
    __syncthreads();
    for (int stride = threads / 2; stride > 0; stride /= 2) {
        if (tid < stride) {
            sum_m[tid] += sum_m[tid + stride];
            sum_gradient[tid] += sum_gradient[tid + stride];
        }
        __syncthreads();
    }
    if (tid == 0) {
        mutual[out] = prefactor * sum_m[0];
        gradient[out] = prefactor * sum_gradient[0];
    }
}

template <typename Kernel>
void launch(Kernel kernel, dim3 grid, int threads, cudaStream_t stream,
            const coilgun::simulation::cuda::MutualPipelineView& view) {
    kernel<<<grid, threads, 0, stream>>>(
        view.coils, view.filaments, view.separations, view.active_mask,
        view.mutual, view.gradient, static_cast<int>(view.batch_size),
        static_cast<int>(view.stage_count), static_cast<int>(view.filament_count),
        view.n_nodes);
    const cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
        throw std::runtime_error(cudaGetErrorString(error));
}

} // namespace

namespace coilgun::simulation::cuda {

void launch_mutual_pipeline(const MutualPipelineView& view, GpuOptLevel opt_level,
                            int threads_per_block, cudaStream_t stream) {
    if (!view.coils || !view.filaments || !view.separations || !view.active_mask
        || !view.mutual || !view.gradient)
        throw std::invalid_argument("mutual pipeline requires non-null buffers");
    const auto device_pointer = [](const void* pointer) {
        cudaPointerAttributes attributes{};
        const auto status = cudaPointerGetAttributes(&attributes, pointer);
#if CUDART_VERSION >= 10000
        return status == cudaSuccess &&
               (attributes.type == cudaMemoryTypeDevice || attributes.type == cudaMemoryTypeManaged);
#else
        return status == cudaSuccess && attributes.memoryType == cudaMemoryTypeDevice;
#endif
    };
    if (!device_pointer(view.coils) || !device_pointer(view.filaments) ||
        !device_pointer(view.separations) || !device_pointer(view.active_mask) ||
        !device_pointer(view.mutual) || !device_pointer(view.gradient))
        throw std::invalid_argument("mutual pipeline buffers must be device-accessible");
    if (view.batch_size == 0 || view.stage_count == 0 || view.filament_count == 0
        || view.batch_size > 0x7fffffff || view.stage_count > 0x7fffffff
        || view.filament_count > 0x7fffffff)
        throw std::invalid_argument("invalid mutual pipeline dimensions");
    if (view.n_nodes != 9)
        throw std::invalid_argument("CUDA mutual pipeline requires n_nodes == 9");
    if (threads_per_block <= 0 || threads_per_block > 512 ||
        (threads_per_block & (threads_per_block - 1)) != 0)
        throw std::invalid_argument("threads_per_block must be a positive power of two no greater than 512");
    std::size_t slots = 0;
    if (!checked_product(view.batch_size, view.stage_count, slots) ||
        !checked_product(slots, view.filament_count, slots))
        throw std::invalid_argument("mutual pipeline dimensions overflow");
    if (view.n_nodes > std::numeric_limits<int>::max() / view.n_nodes ||
        view.n_nodes * view.n_nodes > std::numeric_limits<int>::max() / view.n_nodes ||
        view.n_nodes * view.n_nodes * view.n_nodes > std::numeric_limits<int>::max() / view.n_nodes)
        throw std::invalid_argument("mutual pipeline quadrature dimensions overflow");
    // Geometry buffers are device-resident. Validate their host-side source
    // data before upload; dereferencing them here would be an invalid host read.
    int current_device = -1;
    cudaDeviceProp properties{};
    if (cudaGetDevice(&current_device) != cudaSuccess || current_device < 0 ||
        cudaGetDeviceProperties(&properties, current_device) != cudaSuccess ||
        view.filament_count > static_cast<std::size_t>(properties.maxGridSize[0]) ||
        view.stage_count > static_cast<std::size_t>(properties.maxGridSize[1]) ||
        view.batch_size > static_cast<std::size_t>(properties.maxGridSize[2]))
        throw std::invalid_argument("mutual pipeline grid exceeds device limits");

    const dim3 grid(static_cast<unsigned>(view.filament_count),
                    static_cast<unsigned>(view.stage_count),
                    static_cast<unsigned>(view.batch_size));
    switch (opt_level) {
    case GpuOptLevel::Standard:
        launch(mutual_pipeline_kernel<false, false>, grid, threads_per_block, stream, view);
        break;
    case GpuOptLevel::Full:
        launch(mutual_pipeline_kernel<false, true>, grid, threads_per_block, stream, view);
        break;
    case GpuOptLevel::Aggressive:
        launch(mutual_pipeline_kernel<true, true>, grid, threads_per_block, stream, view);
        break;
    default:
        throw std::invalid_argument("unknown GPU optimization level");
    }
}

void initialize_mutual_pipeline_constants(cudaStream_t stream) {
    cudaStreamCaptureStatus capture_status = cudaStreamCaptureStatusNone;
    const auto capture_query = cudaStreamIsCapturing(stream, &capture_status);
    if (capture_query != cudaSuccess)
        throw std::runtime_error("CUDA stream capture status query failed");
    if (capture_status != cudaStreamCaptureStatusNone)
        throw std::invalid_argument("quadrature constants must be initialized before graph capture");
    const auto& quadrature = physics::gauss_legendre(9);
    const auto nodes = cudaMemcpyToSymbolAsync(
        pipeline_gl_nodes, quadrature.nodes.data(), 9 * sizeof(double), 0,
        cudaMemcpyHostToDevice, stream);
    const auto weights = cudaMemcpyToSymbolAsync(
        pipeline_gl_weights, quadrature.weights.data(), 9 * sizeof(double), 0,
        cudaMemcpyHostToDevice, stream);
    if (nodes != cudaSuccess || weights != cudaSuccess)
        throw std::runtime_error("Gauss-Legendre constant upload failed");
}

std::size_t mutual_pipeline_index(std::size_t simulation, std::size_t stage,
                                  std::size_t filament, std::size_t stage_count,
                                  std::size_t filament_count) {
    return (simulation * stage_count + stage) * filament_count + filament;
}

} // namespace coilgun::simulation::cuda
