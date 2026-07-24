#include "coilgun/simulation/cuda/gpu_solver.hpp"
#include "coilgun/simulation/cuda/gpu_execution_context.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <limits>

namespace coilgun::simulation::cuda {
namespace {
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Vector = Eigen::VectorXd;

bool finite_values(const double* p, std::size_t n) noexcept {
    if (!p) return false;
    for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(p[i])) return false;
    return true;
}
bool checked_size(std::size_t value) noexcept { return value <= static_cast<std::size_t>(std::numeric_limits<int>::max()); }
bool checked_product(std::size_t a, std::size_t b, std::size_t& out) noexcept {
    if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return false;
    out = a * b;
    return true;
}
struct DeviceGuard {
    explicit DeviceGuard(int target) : target(target) {
        if (cudaGetDevice(&previous) != cudaSuccess) throw std::runtime_error("cudaGetDevice failed");
        if (previous != target && cudaSetDevice(target) != cudaSuccess) throw std::runtime_error("cudaSetDevice failed");
    }
    ~DeviceGuard() noexcept { if (previous >= 0 && previous != target) (void)cudaSetDevice(previous); }
    int target; int previous = -1;
};

__global__ void residual_kernel(const double* matrices, const double* rhs,
                                const double* solutions, double* residuals,
                                std::size_t batch_size, std::size_t dimension) {
    const std::size_t batch = blockIdx.x * blockDim.x + threadIdx.x;
    if (batch >= batch_size) return;

    const double* matrix = matrices + batch * dimension * dimension;
    const double* vector = rhs + batch * dimension;
    const double* solution = solutions + batch * dimension;
    double maximum = 0.0;
    double scale = 1.0;
    for (std::size_t row = 0; row < dimension; ++row) {
        double reconstructed = 0.0;
        for (std::size_t column = 0; column < dimension; ++column)
            reconstructed += matrix[row * dimension + column] * solution[column];
        maximum = fmax(maximum, fabs(reconstructed - vector[row]));
        scale = fmax(scale, fabs(vector[row]));
    }
    residuals[batch] = maximum / scale;
}
}

struct GpuSolver::Impl {
    SolverMode requested;
    SolverMode resolved;
    SolverBatchLayout layout;
    SolverWorkspace workspace;
    GpuExecutionContext* context = nullptr;
    Matrix matrix;
    Vector rhs;
    Vector solution;
    Eigen::LDLT<Matrix> ldlt;
    double* d_matrix = nullptr;
    double* d_rhs = nullptr;
    double** d_matrix_ptrs = nullptr;
    double** d_rhs_ptrs = nullptr;
    double** d_device_rhs_ptrs = nullptr;
    int* d_pivots = nullptr;
    int* d_info = nullptr;
    std::vector<double> column_major;
    std::vector<double*> matrices;
    std::vector<double*> rhs_ptrs;
    std::vector<double*> device_rhs_ptrs;
    const double* device_rhs_base = nullptr;
    std::vector<int> info;
    std::vector<int> solve_info;
    std::vector<double> residuals;

    ~Impl() { release(); }
    void release() noexcept {
        if (context) {
            try { context->synchronize(); } catch (...) { }
        }
        if (context) { try { DeviceGuard guard(context->device_id()); (void)cudaFree(d_matrix); (void)cudaFree(d_rhs); (void)cudaFree(d_matrix_ptrs); (void)cudaFree(d_rhs_ptrs); (void)cudaFree(d_device_rhs_ptrs); (void)cudaFree(d_pivots); (void)cudaFree(d_info); } catch (...) {} }
        d_matrix = d_rhs = nullptr;
        d_matrix_ptrs = d_rhs_ptrs = d_device_rhs_ptrs = nullptr; d_pivots = d_info = nullptr;
        device_rhs_base = nullptr;
    }
};

SolverBatchLayout::SolverBatchLayout(std::size_t batch, std::size_t dimension)
    : batch_size(batch), dimension(dimension) {
    std::size_t dd = 0, matrix_count = 0, vector_count = 0;
    if (!batch || !dimension || !checked_size(batch) || !checked_size(dimension) ||
        !checked_product(dimension, dimension, dd) || !checked_product(batch, dd, matrix_count) ||
        !checked_product(batch, dimension, vector_count) ||
        matrix_count > std::numeric_limits<std::size_t>::max() / sizeof(double) ||
        vector_count > std::numeric_limits<std::size_t>::max() / sizeof(double))
        throw std::invalid_argument("solver batch layout overflows host or CUDA dimensions");
}
bool SolverBatchLayout::operator==(const SolverBatchLayout& o) const noexcept {
    return batch_size == o.batch_size && dimension == o.dimension;
}
SolverStatus SolverStatus::success(double residual) { SolverStatus s; s.ok = true; s.max_residual = residual; return s; }
SolverStatus SolverStatus::failure_status(SolverFailure f, std::string m) { SolverStatus s; s.failure = f; s.message = std::move(m); return s; }

GpuSolver::GpuSolver(SolverMode mode, SolverBatchLayout layout)
    : impl_(new Impl{mode, mode == SolverMode::Auto ? SolverMode::Eigen : mode, layout}) {}
GpuSolver::GpuSolver(GpuExecutionContext& c, SolverMode mode, SolverBatchLayout layout)
    : impl_(new Impl{mode, mode == SolverMode::Auto ? SolverMode::Eigen : mode, layout, {}, &c}) {}
GpuSolver::GpuSolver(GpuExecutionContext& c, const GpuExecutionPolicy& p, SolverBatchLayout layout)
    : impl_(new Impl{p.requested_solver, p.solver, layout, {}, &c}) {}
GpuSolver::~GpuSolver() { delete impl_; }
GpuSolver::GpuSolver(GpuSolver&& o) noexcept : impl_(o.impl_) { o.impl_ = nullptr; }
GpuSolver& GpuSolver::operator=(GpuSolver&& o) noexcept { if (this != &o) { delete impl_; impl_ = o.impl_; o.impl_ = nullptr; } return *this; }

SolverMode GpuSolver::requested_mode() const noexcept { return impl_ ? impl_->requested : SolverMode::Auto; }
SolverMode GpuSolver::resolved_mode() const noexcept { return impl_ ? impl_->resolved : SolverMode::Auto; }
const SolverBatchLayout& GpuSolver::layout() const noexcept { static const SolverBatchLayout empty{1, 1}; return impl_ ? impl_->layout : empty; }
const SolverWorkspace& GpuSolver::workspace() const noexcept { static const SolverWorkspace empty{}; return impl_ ? impl_->workspace : empty; }

SolverStatus GpuSolver::initialize_workspace() { return initialize_workspace(layout()); }
SolverStatus GpuSolver::initialize_workspace(const SolverBatchLayout& requested_layout) {
    if (!impl_) return SolverStatus::failure_status(SolverFailure::NotInitialized, "moved-from solver");
    if (!(requested_layout == impl_->layout)) return SolverStatus::failure_status(SolverFailure::LayoutMismatch, "solver workspace layout does not match solver layout");
    if (impl_->workspace.initialized) return SolverStatus::success();
    const auto B = impl_->layout.batch_size, D = impl_->layout.dimension;
    if (!checked_size(B) || !checked_size(D)) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "solver dimensions exceed CUDA int range");
    impl_->matrix.resize(static_cast<Eigen::Index>(D), static_cast<Eigen::Index>(D));
    impl_->rhs.resize(static_cast<Eigen::Index>(D)); impl_->solution.resize(static_cast<Eigen::Index>(D));
    if (impl_->resolved == SolverMode::Batched) {
        if (!impl_->context || !impl_->context->valid()) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "batched solver requires a valid execution context");
        std::size_t matrix_count = 0, vector_count = 0;
        checked_product(B, D * D, matrix_count); checked_product(B, D, vector_count);
        const std::size_t matrix_bytes = matrix_count * sizeof(double), vector_bytes = vector_count * sizeof(double);
        DeviceGuard guard(impl_->context->device_id());
        auto alloc = [](void** p, std::size_t n) { return cudaMalloc(p, n); };
        if (alloc(reinterpret_cast<void**>(&impl_->d_matrix), matrix_bytes) != cudaSuccess ||
            alloc(reinterpret_cast<void**>(&impl_->d_rhs), vector_bytes) != cudaSuccess ||
            alloc(reinterpret_cast<void**>(&impl_->d_matrix_ptrs), B * sizeof(double*)) != cudaSuccess ||
            alloc(reinterpret_cast<void**>(&impl_->d_rhs_ptrs), B * sizeof(double*)) != cudaSuccess ||
            alloc(reinterpret_cast<void**>(&impl_->d_device_rhs_ptrs), B * sizeof(double*)) != cudaSuccess ||
            alloc(reinterpret_cast<void**>(&impl_->d_pivots), B * D * sizeof(int)) != cudaSuccess ||
            alloc(reinterpret_cast<void**>(&impl_->d_info), B * sizeof(int)) != cudaSuccess) {
            impl_->release(); return SolverStatus::failure_status(SolverFailure::InvalidArgument, "CUDA workspace allocation failed");
        }
        impl_->matrices.resize(B); impl_->rhs_ptrs.resize(B); impl_->device_rhs_ptrs.resize(B);
        impl_->info.resize(B); impl_->solve_info.resize(B); impl_->residuals.resize(B);
        impl_->column_major.resize(matrix_count);
        for (std::size_t b = 0; b < B; ++b) { impl_->matrices[b] = impl_->d_matrix + b * D * D; impl_->rhs_ptrs[b] = impl_->d_rhs + b * D; }
        if (cudaMemcpy(impl_->d_matrix_ptrs, impl_->matrices.data(), B * sizeof(double*), cudaMemcpyHostToDevice) != cudaSuccess ||
            cudaMemcpy(impl_->d_rhs_ptrs, impl_->rhs_ptrs.data(), B * sizeof(double*), cudaMemcpyHostToDevice) != cudaSuccess) {
            impl_->release(); return SolverStatus::failure_status(SolverFailure::InvalidArgument, "CUDA pointer table copy failed");
        }
    }
    impl_->workspace.initialized = true; impl_->workspace.batch_capacity = B; impl_->workspace.dimension_capacity = D; ++impl_->workspace.allocation_count;
    return SolverStatus::success();
}

SolverStatus GpuSolver::solve(const double* m, const double* r, double* x) { return solve(layout(), m, r, x); }
SolverStatus GpuSolver::solve(const SolverBatchLayout& l, const double* m, const double* r, double* x) {
    if (!impl_) return SolverStatus::failure_status(SolverFailure::NotInitialized, "moved-from solver");
    if (!(l == impl_->layout) || l.batch_size != 1) return SolverStatus::failure_status(SolverFailure::LayoutMismatch, "single solve requires the configured B=1 layout");
    if (!impl_->workspace.initialized) return SolverStatus::failure_status(SolverFailure::NotInitialized, "solver workspace has not been initialized");
    return solve_at(0, m, r, x);
}

SolverStatus GpuSolver::solve_at(std::size_t b, const double* m, const double* r, double* x) {
    const auto D = layout().dimension;
    if (b >= layout().batch_size) return SolverStatus::failure_status(SolverFailure::LayoutMismatch, "batch index is outside the configured layout");
    if (!finite_values(m, D * D) || !finite_values(r, D)) return SolverStatus::failure_status(SolverFailure::NonFiniteInput, "matrix and rhs must contain finite values");
    if (!x) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "solution pointer must not be null");
    if (impl_->resolved == SolverMode::Batched) return solve_batch(m, r, x);
    Eigen::Map<const Matrix> input(m, D, D); Eigen::Map<const Vector> input_rhs(r, D);
    impl_->matrix = input; impl_->rhs = input_rhs; impl_->ldlt.compute(impl_->matrix);
    if (impl_->ldlt.info() != Eigen::Success) return SolverStatus::failure_status(SolverFailure::FactorizationFailed, "Eigen LDLT factorization failed");
    impl_->solution = impl_->ldlt.solve(impl_->rhs);
    if (!finite_values(impl_->solution.data(), D)) return SolverStatus::failure_status(SolverFailure::NonFiniteOutput, "Eigen LDLT produced a non-finite solution");
    std::copy(impl_->solution.data(), impl_->solution.data() + D, x); return check_residual(m, r, x);
}

SolverStatus GpuSolver::solve_batch(const double* matrices, const double* rhs, double* solutions) {
    if (!impl_) return SolverStatus::failure_status(SolverFailure::NotInitialized, "moved-from solver");
    if (!impl_->workspace.initialized) return SolverStatus::failure_status(SolverFailure::NotInitialized, "solver workspace has not been initialized");
    if (!matrices || !rhs || !solutions) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "batch buffers must not be null");
    const auto B = layout().batch_size, D = layout().dimension;
    if (impl_->resolved != SolverMode::Batched) {
        double max = 0; for (std::size_t b = 0; b < B; ++b) { auto s = solve_at(b, matrices + b * D * D, rhs + b * D, solutions + b * D); if (!s.ok) { s.failed_batch = static_cast<int>(b); return s; } max = std::max(max, s.max_residual); } return SolverStatus::success(max);
    }
     std::fill(impl_->info.begin(), impl_->info.end(), 0);
     std::fill(impl_->solve_info.begin(), impl_->solve_info.end(), 0);
    for (std::size_t b = 0; b < B; ++b) for (std::size_t row = 0; row < D; ++row) for (std::size_t col = 0; col < D; ++col) impl_->column_major[b * D * D + col * D + row] = matrices[b * D * D + row * D + col];
    try { DeviceGuard guard(impl_->context->device_id());
    if (cudaMemcpyAsync(impl_->d_matrix, impl_->column_major.data(), impl_->column_major.size() * sizeof(double), cudaMemcpyHostToDevice, impl_->context->stream()) != cudaSuccess ||
        cudaMemcpyAsync(impl_->d_rhs, rhs, B * D * sizeof(double), cudaMemcpyHostToDevice, impl_->context->stream()) != cudaSuccess) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "CUDA host-to-device copy failed");
    // CUDA 12.8 exposes batched LU/solve through cuBLAS, not cuSOLVER.
    const auto factor = cublasDgetrfBatched(impl_->context->cublas(), static_cast<int>(D), impl_->d_matrix_ptrs, static_cast<int>(D), impl_->d_pivots, impl_->d_info, static_cast<int>(B));
    if (factor != CUBLAS_STATUS_SUCCESS) return SolverStatus::failure_status(SolverFailure::FactorizationFailed, "cublasDgetrfBatched failed");
    if (cudaMemcpyAsync(impl_->info.data(), impl_->d_info, B * sizeof(int), cudaMemcpyDeviceToHost, impl_->context->stream()) != cudaSuccess) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "CUDA info copy failed");
    if (cudaStreamSynchronize(impl_->context->stream()) != cudaSuccess) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "cudaStreamSynchronize factorization");
    for (std::size_t b = 0; b < B; ++b) {
        if (impl_->info[b] != 0) {
            auto s = SolverStatus::failure_status(
                SolverFailure::FactorizationFailed,
                "cuSOLVER failed at batch " + std::to_string(b));
            s.failed_batch = static_cast<int>(b);
            s.backend_info = impl_->info[b];
            return s;
        }
    }
    // cublas<t>getrsBatched reports per-system status to host memory, unlike
    // getrfBatched's device info array.
     const auto solve = cublasDgetrsBatched(impl_->context->cublas(), CUBLAS_OP_N, static_cast<int>(D), 1, const_cast<const double* const*>(impl_->d_matrix_ptrs), static_cast<int>(D), impl_->d_pivots, impl_->d_rhs_ptrs, static_cast<int>(D), impl_->solve_info.data(), static_cast<int>(B));
    if (solve != CUBLAS_STATUS_SUCCESS) return SolverStatus::failure_status(SolverFailure::FactorizationFailed, "cublasDgetrsBatched failed");
     if (cudaMemcpyAsync(solutions, impl_->d_rhs, B * D * sizeof(double), cudaMemcpyDeviceToHost, impl_->context->stream()) != cudaSuccess) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "CUDA solution copy failed");
    if (cudaStreamSynchronize(impl_->context->stream()) != cudaSuccess) return SolverStatus::failure_status(SolverFailure::InvalidArgument, "cudaStreamSynchronize solve");
     for (std::size_t b = 0; b < B; ++b) if (impl_->solve_info[b] != 0) { auto s = SolverStatus::failure_status(SolverFailure::FactorizationFailed, "batched solve failed at batch " + std::to_string(b)); s.failed_batch = static_cast<int>(b); s.backend_info = impl_->solve_info[b]; return s; }
    for (std::size_t b = 0; b < B; ++b) { if (!finite_values(solutions + b * D, D)) { auto s = SolverStatus::failure_status(SolverFailure::NonFiniteOutput, "cuSOLVER produced a non-finite solution"); s.failed_batch = static_cast<int>(b); return s; } auto s = check_residual(matrices + b * D * D, rhs + b * D, solutions + b * D); if (!s.ok) { s.failed_batch = static_cast<int>(b); return s; } }
    } catch (const std::exception& error) { return SolverStatus::failure_status(SolverFailure::InvalidArgument, error.what()); }
    return SolverStatus::success();
}

SolverStatus GpuSolver::solve_device(const DeviceMatrixView& matrix,
                                     const DeviceVectorView& rhs,
                                     DeviceVectorView solution,
                                     DeviceResidualView residual) {
    if (!impl_ || !impl_->workspace.initialized)
        return SolverStatus::failure_status(SolverFailure::NotInitialized,
                                            "solver workspace has not been initialized");
    const auto B = layout().batch_size;
    const auto D = layout().dimension;
    if (impl_->resolved != SolverMode::Batched || !impl_->context)
        return SolverStatus::failure_status(SolverFailure::UnsupportedMode,
                                            "device solve requires batched CUDA solver mode");
    if (!matrix.data || !rhs.data || !solution.data || matrix.batch_size != B ||
        rhs.batch_size != B || solution.batch_size != B || matrix.dimension != D ||
        rhs.dimension != D || solution.dimension != D ||
        (residual.values != nullptr && residual.count != B))
        return SolverStatus::failure_status(SolverFailure::LayoutMismatch,
                                            "device solver views do not match configured layout");
    try {
        DeviceGuard guard(impl_->context->device_id());
        cudaStreamCaptureStatus capture_status = cudaStreamCaptureStatusNone;
        if (cudaStreamIsCapturing(impl_->context->stream(), &capture_status) != cudaSuccess)
            return SolverStatus::failure_status(SolverFailure::InvalidArgument,
                                                "unable to query solver stream capture status");
        const bool capturing = capture_status != cudaStreamCaptureStatusNone;
        double** rhs_ptrs = impl_->d_rhs_ptrs;
        double* rhs_workspace = impl_->d_rhs;
        if (!capturing) {
            // The caller's solution buffer is also the cuBLAS RHS output. This
            // removes the old workspace-RHS-to-solution copy while retaining
            // the input RHS and matrix contracts.
            if (impl_->device_rhs_base != solution.data) {
                for (std::size_t b = 0; b < B; ++b)
                    impl_->device_rhs_ptrs[b] = solution.data + b * D;
                if (cudaMemcpyAsync(impl_->d_device_rhs_ptrs,
                                    impl_->device_rhs_ptrs.data(), B * sizeof(double*),
                                    cudaMemcpyHostToDevice,
                                    impl_->context->stream()) != cudaSuccess)
                    return SolverStatus::failure_status(
                        SolverFailure::InvalidArgument,
                        "CUDA device solution pointer table copy failed");
                impl_->device_rhs_base = solution.data;
            }
            rhs_ptrs = impl_->d_device_rhs_ptrs;
            rhs_workspace = solution.data;
        }
        if (cudaMemcpyAsync(impl_->d_matrix, matrix.data, B * D * D * sizeof(double),
                            cudaMemcpyDeviceToDevice, impl_->context->stream()) != cudaSuccess ||
            cudaMemcpyAsync(rhs_workspace, rhs.data, B * D * sizeof(double),
                            cudaMemcpyDeviceToDevice, impl_->context->stream()) != cudaSuccess)
            return SolverStatus::failure_status(SolverFailure::InvalidArgument,
                                                "device solver pointer setup failed");
        const auto factor = cublasDgetrfBatched(
            impl_->context->cublas(), static_cast<int>(D), impl_->d_matrix_ptrs,
            static_cast<int>(D), impl_->d_pivots, impl_->d_info, static_cast<int>(B));
        if (factor != CUBLAS_STATUS_SUCCESS)
            return SolverStatus::failure_status(SolverFailure::FactorizationFailed,
                                                "cublasDgetrfBatched failed for device views");
        const auto solve = cublasDgetrsBatched(
            impl_->context->cublas(), CUBLAS_OP_N, static_cast<int>(D), 1,
            const_cast<const double* const*>(impl_->d_matrix_ptrs), static_cast<int>(D),
            impl_->d_pivots, rhs_ptrs, static_cast<int>(D),
            impl_->solve_info.data(), static_cast<int>(B));
        if (solve != CUBLAS_STATUS_SUCCESS)
            return SolverStatus::failure_status(SolverFailure::FactorizationFailed,
                                                "cublasDgetrsBatched failed for device views");

        if (capturing &&
            cudaMemcpyAsync(solution.data, impl_->d_rhs, B * D * sizeof(double),
                            cudaMemcpyDeviceToDevice, impl_->context->stream()) != cudaSuccess)
            return SolverStatus::failure_status(SolverFailure::InvalidArgument,
                                                "device solver solution copy failed");

        if (residual.values != nullptr) {
            residual_kernel<<<static_cast<unsigned int>((B + 127) / 128), 128, 0,
                              impl_->context->stream()>>>(
                matrix.data, rhs.data, solution.data, residual.values, B, D);
            if (cudaGetLastError() != cudaSuccess)
                return SolverStatus::failure_status(SolverFailure::InvalidArgument,
                                                    "device residual kernel launch failed");
        }

        if (capture_status != cudaStreamCaptureStatusNone) return SolverStatus::success();
        if (cudaMemcpyAsync(impl_->info.data(), impl_->d_info, B * sizeof(int),
                            cudaMemcpyDeviceToHost, impl_->context->stream()) != cudaSuccess ||
            (residual.values != nullptr &&
             cudaMemcpyAsync(impl_->residuals.data(), residual.values, B * sizeof(double),
                             cudaMemcpyDeviceToHost, impl_->context->stream()) != cudaSuccess) ||
            cudaStreamSynchronize(impl_->context->stream()) != cudaSuccess)
            return SolverStatus::failure_status(SolverFailure::InvalidArgument,
                                                "device solver status synchronization failed");
        double maximum_residual = 0.0;
        for (std::size_t batch = 0; batch < B; ++batch) {
            if (impl_->info[batch] != 0 || impl_->solve_info[batch] != 0) {
                auto status = SolverStatus::failure_status(
                    SolverFailure::FactorizationFailed,
                    "device solver factorization failed at batch " + std::to_string(batch));
                status.failed_batch = static_cast<int>(batch);
                status.backend_info = impl_->info[batch] != 0
                    ? impl_->info[batch] : impl_->solve_info[batch];
                return status;
            }
            if (residual.values != nullptr &&
                (!std::isfinite(impl_->residuals[batch]) || impl_->residuals[batch] > 1.0e-10)) {
                auto status = SolverStatus::failure_status(
                    std::isfinite(impl_->residuals[batch])
                        ? SolverFailure::ResidualTooLarge : SolverFailure::NonFiniteOutput,
                    "device solver residual check failed at batch " + std::to_string(batch));
                status.failed_batch = static_cast<int>(batch);
                status.max_residual = impl_->residuals[batch];
                return status;
            }
            if (residual.values != nullptr)
                maximum_residual = std::max(maximum_residual, impl_->residuals[batch]);
        }
        return SolverStatus::success(maximum_residual);
    } catch (const std::exception& error) {
        return SolverStatus::failure_status(SolverFailure::InvalidArgument, error.what());
    }
}

SolverStatus GpuSolver::validate_device_result(DeviceResidualView residual) {
    if (!impl_ || !impl_->workspace.initialized)
        return SolverStatus::failure_status(SolverFailure::NotInitialized,
                                            "solver workspace has not been initialized");
    const auto B = layout().batch_size;
    if (impl_->resolved != SolverMode::Batched || !impl_->context ||
        !residual.values || residual.count != B)
        return SolverStatus::failure_status(SolverFailure::LayoutMismatch,
                                            "device solver status view does not match layout");
    try {
        DeviceGuard guard(impl_->context->device_id());
        if (cudaMemcpyAsync(impl_->info.data(), impl_->d_info, B * sizeof(int),
                            cudaMemcpyDeviceToHost, impl_->context->stream()) != cudaSuccess ||
            cudaMemcpyAsync(impl_->residuals.data(), residual.values, B * sizeof(double),
                            cudaMemcpyDeviceToHost, impl_->context->stream()) != cudaSuccess ||
            cudaStreamSynchronize(impl_->context->stream()) != cudaSuccess)
            return SolverStatus::failure_status(SolverFailure::InvalidArgument,
                                                "device solver status synchronization failed");
        double maximum_residual = 0.0;
        for (std::size_t batch = 0; batch < B; ++batch) {
            if (impl_->info[batch] != 0 || impl_->solve_info[batch] != 0) {
                auto status = SolverStatus::failure_status(
                    SolverFailure::FactorizationFailed,
                    "device solver factorization failed at batch " + std::to_string(batch));
                status.failed_batch = static_cast<int>(batch);
                status.backend_info = impl_->info[batch] != 0
                    ? impl_->info[batch] : impl_->solve_info[batch];
                return status;
            }
            if (!std::isfinite(impl_->residuals[batch]) || impl_->residuals[batch] > 1.0e-10) {
                auto status = SolverStatus::failure_status(
                    std::isfinite(impl_->residuals[batch])
                        ? SolverFailure::ResidualTooLarge : SolverFailure::NonFiniteOutput,
                    "device solver residual check failed at batch " + std::to_string(batch));
                status.failed_batch = static_cast<int>(batch);
                status.max_residual = impl_->residuals[batch];
                return status;
            }
            maximum_residual = std::max(maximum_residual, impl_->residuals[batch]);
        }
        return SolverStatus::success(maximum_residual);
    } catch (const std::exception& error) {
        return SolverStatus::failure_status(SolverFailure::InvalidArgument, error.what());
    }
}

SolverStatus GpuSolver::check_residual(const double* m, const double* r, const double* x) const {
    const auto D = layout().dimension; if (!finite_values(m, D * D) || !finite_values(r, D)) return SolverStatus::failure_status(SolverFailure::NonFiniteInput, "matrix and rhs must contain finite values"); if (!finite_values(x, D)) return SolverStatus::failure_status(SolverFailure::NonFiniteOutput, "solution must contain finite values");
    double max = 0, scale = 1; for (std::size_t row = 0; row < D; ++row) { double value = 0; for (std::size_t col = 0; col < D; ++col) value += m[row * D + col] * x[col]; max = std::max(max, std::abs(value - r[row])); scale = std::max(scale, std::abs(r[row])); }
    if (max > 1e-10 * scale) { auto s = SolverStatus::failure_status(SolverFailure::ResidualTooLarge, "solver residual exceeds sanity tolerance"); s.max_residual = max; return s; } return SolverStatus::success(max);
}
}
