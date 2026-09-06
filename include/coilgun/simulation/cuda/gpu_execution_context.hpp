/**
 * @file gpu_execution_context.hpp
 * @brief RAII ownership of the CUDA resources shared by GPU execution paths.
 */

#pragma once

#include <cstddef>
#include <memory>

#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

namespace coilgun::simulation::cuda {

/// Index of the preferred CUDA device: the first discrete (non-integrated)
/// device, or device 0 when only integrated devices exist. Returns -1 when no
/// usable CUDA device is present. The result is resolved once per process.
int preferred_cuda_device() noexcept;

/// True when a usable CUDA device exists; as a side effect the preferred
/// device is made current (cudaSetDevice) so raw runtime calls after this
/// gate land on it.
bool cuda_device_available() noexcept;

struct GpuExecutionContextConfig {
    // -1 requests automatic device placement via preferred_cuda_device().
    int device_id = -1;
    unsigned int stream_flags = cudaStreamNonBlocking;
    std::size_t workspace_bytes = 0;
    // Metadata request only. Host-wall timings are recorded by GpuEngine and
    // this context neither requires nor promises NVTX instrumentation.
    bool enable_profiling = false;
};

class GpuExecutionContext {
public:
    explicit GpuExecutionContext(GpuExecutionContextConfig config = {});
    ~GpuExecutionContext() noexcept;

    GpuExecutionContext(const GpuExecutionContext&) = delete;
    GpuExecutionContext& operator=(const GpuExecutionContext&) = delete;
    GpuExecutionContext(GpuExecutionContext&&) noexcept;
    GpuExecutionContext& operator=(GpuExecutionContext&&) noexcept;

    int device_id() const noexcept;
    cudaStream_t stream() const noexcept;
    cudaEvent_t start_event() const noexcept;
    cudaEvent_t stop_event() const noexcept;
    cublasHandle_t cublas() const noexcept;
    cusolverDnHandle_t cusolver() const noexcept;

    void reserve_workspace(std::size_t bytes);
    void* workspace() const noexcept;
    std::size_t workspace_bytes() const noexcept;
    void synchronize() const;
    void record_start() const;
    void record_stop() const;
    void ensure_quadrature9_loaded();
    bool quadrature9_loaded() const noexcept;
    bool valid() const noexcept { return impl_ != nullptr; }

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace coilgun::simulation::cuda
