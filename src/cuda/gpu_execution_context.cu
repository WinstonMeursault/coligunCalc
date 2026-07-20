/**
 * @file gpu_execution_context.cu
 * @brief Checked RAII implementation of the shared CUDA execution context.
 */

#include "coilgun/simulation/cuda/gpu_execution_context.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <memory>

namespace coilgun::simulation::cuda {
namespace {

[[noreturn]] void throw_cuda(cudaError_t error, const char* operation) {
    throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(error));
}

void check_cuda(cudaError_t error, const char* operation) {
    if (error != cudaSuccess) throw_cuda(error, operation);
}

[[noreturn]] void throw_cublas(cublasStatus_t status, const char* operation) {
    throw std::runtime_error(std::string(operation) + " failed (cuBLAS status " +
                             std::to_string(static_cast<int>(status)) + ")");
}

void check_cublas(cublasStatus_t status, const char* operation) {
    if (status != CUBLAS_STATUS_SUCCESS) throw_cublas(status, operation);
}

[[noreturn]] void throw_cusolver(cusolverStatus_t status, const char* operation) {
    throw std::runtime_error(std::string(operation) + " failed (cuSOLVER status " +
                             std::to_string(static_cast<int>(status)) + ")");
}

void check_cusolver(cusolverStatus_t status, const char* operation) {
    if (status != CUSOLVER_STATUS_SUCCESS) throw_cusolver(status, operation);
}

struct ScopedDevice {
    explicit ScopedDevice(int target) : target(target) {
        check_cuda(cudaGetDevice(&previous), "cudaGetDevice");
        if (previous != target) check_cuda(cudaSetDevice(target), "cudaSetDevice");
    }
    ~ScopedDevice() noexcept {
        if (previous >= 0 && previous != target) cudaSetDevice(previous);
    }
    int target = -1;
    int previous = -1;
};

} // namespace

struct GpuExecutionContext::Impl {
    explicit Impl(GpuExecutionContextConfig config) : device_id(config.device_id) {
        check_cuda(cudaGetDevice(&previous_device), "cudaGetDevice");
        check_cuda(cudaSetDevice(device_id), "cudaSetDevice");

        try {
            check_cuda(cudaStreamCreateWithFlags(&stream, config.stream_flags),
                       "cudaStreamCreateWithFlags");
            check_cuda(cudaEventCreateWithFlags(&start_event, cudaEventDefault),
                       "cudaEventCreateWithFlags(start)");
            check_cuda(cudaEventCreateWithFlags(&stop_event, cudaEventDefault),
                       "cudaEventCreateWithFlags(stop)");
            check_cublas(cublasCreate(&cublas), "cublasCreate");
            check_cusolver(cusolverDnCreate(&cusolver), "cusolverDnCreate");
            check_cublas(cublasSetStream(cublas, stream), "cublasSetStream");
            check_cusolver(cusolverDnSetStream(cusolver, stream), "cusolverDnSetStream");
            if (config.workspace_bytes != 0) reserve_workspace(config.workspace_bytes);
            check_cuda(cudaSetDevice(previous_device), "cudaSetDevice(restore constructor device)");
            previous_device = -1;
        } catch (...) {
            cleanup();
            throw;
        }
    }

    ~Impl() noexcept { cleanup(); }

    void reserve_workspace(std::size_t bytes) {
        if (bytes <= workspace_bytes) return;
        ScopedDevice device(device_id);
        check_cuda(cudaStreamSynchronize(stream), "cudaStreamSynchronize(workspace resize)");
        void* replacement = nullptr;
        check_cuda(cudaMalloc(&replacement, bytes), "cudaMalloc(workspace)");
        if (workspace != nullptr) cudaFree(workspace);
        workspace = replacement;
        workspace_bytes = bytes;
    }

    void cleanup() noexcept {
        std::unique_ptr<ScopedDevice> device;
        if (device_id >= 0) { try { device = std::make_unique<ScopedDevice>(device_id); } catch (...) {} }
        if (stream != nullptr) cudaStreamSynchronize(stream);
        if (cusolver != nullptr) cusolverDnDestroy(cusolver);
        if (cublas != nullptr) cublasDestroy(cublas);
        if (stop_event != nullptr) cudaEventDestroy(stop_event);
        if (start_event != nullptr) cudaEventDestroy(start_event);
        if (stream != nullptr) cudaStreamDestroy(stream);
        if (workspace != nullptr) cudaFree(workspace);
        cusolver = nullptr;
        cublas = nullptr;
        stop_event = nullptr;
        start_event = nullptr;
        stream = nullptr;
        workspace = nullptr;
        workspace_bytes = 0;
        device.reset();
        if (previous_device >= 0) cudaSetDevice(previous_device);
        previous_device = -1;
    }

    int device_id = -1;
    int previous_device = -1;
    cudaStream_t stream = nullptr;
    cudaEvent_t start_event = nullptr;
    cudaEvent_t stop_event = nullptr;
    cublasHandle_t cublas = nullptr;
    cusolverDnHandle_t cusolver = nullptr;
    void* workspace = nullptr;
    std::size_t workspace_bytes = 0;
};

GpuExecutionContext::GpuExecutionContext(GpuExecutionContextConfig config)
    : impl_(std::make_unique<Impl>(config)) {}

GpuExecutionContext::~GpuExecutionContext() noexcept = default;
GpuExecutionContext::GpuExecutionContext(GpuExecutionContext&&) noexcept = default;
GpuExecutionContext& GpuExecutionContext::operator=(GpuExecutionContext&&) noexcept = default;

int GpuExecutionContext::device_id() const noexcept { return impl_ ? impl_->device_id : -1; }
cudaStream_t GpuExecutionContext::stream() const noexcept { return impl_ ? impl_->stream : nullptr; }
cudaEvent_t GpuExecutionContext::start_event() const noexcept { return impl_ ? impl_->start_event : nullptr; }
cudaEvent_t GpuExecutionContext::stop_event() const noexcept { return impl_ ? impl_->stop_event : nullptr; }
cublasHandle_t GpuExecutionContext::cublas() const noexcept { return impl_ ? impl_->cublas : nullptr; }
cusolverDnHandle_t GpuExecutionContext::cusolver() const noexcept { return impl_ ? impl_->cusolver : nullptr; }

void GpuExecutionContext::reserve_workspace(std::size_t bytes) { if (!impl_) throw std::logic_error("moved-from GPU context"); impl_->reserve_workspace(bytes); }
void* GpuExecutionContext::workspace() const noexcept { return impl_ ? impl_->workspace : nullptr; }
std::size_t GpuExecutionContext::workspace_bytes() const noexcept { return impl_ ? impl_->workspace_bytes : 0; }
void GpuExecutionContext::synchronize() const { if (!impl_) throw std::logic_error("moved-from GPU context"); ScopedDevice device(impl_->device_id); check_cuda(cudaStreamSynchronize(impl_->stream), "cudaStreamSynchronize"); }
void GpuExecutionContext::record_start() const { if (!impl_) throw std::logic_error("moved-from GPU context"); ScopedDevice device(impl_->device_id); check_cuda(cudaEventRecord(impl_->start_event, impl_->stream), "cudaEventRecord(start)"); }
void GpuExecutionContext::record_stop() const { if (!impl_) throw std::logic_error("moved-from GPU context"); ScopedDevice device(impl_->device_id); check_cuda(cudaEventRecord(impl_->stop_event, impl_->stream), "cudaEventRecord(stop)"); }

} // namespace coilgun::simulation::cuda

#ifdef COILGUN_GPU_CONTEXT_SMOKE
int main() {
    coilgun::simulation::cuda::GpuExecutionContext context({0, cudaStreamNonBlocking, 4096});
    context.record_start();
    context.record_stop();
    context.synchronize();
    return context.workspace_bytes() == 4096 && context.cublas() != nullptr &&
                   context.cusolver() != nullptr
               ? 0
               : 1;
}
#endif
