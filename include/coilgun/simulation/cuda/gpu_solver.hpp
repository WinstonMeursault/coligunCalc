/**
 * @file gpu_solver.hpp
 * @brief Unified dense solver contract shared by CPU and CUDA solver paths.
 */

#pragma once

#include "coilgun/simulation/cuda/gpu_execution_config.hpp"

#include <cstddef>
#include <string>

namespace coilgun::simulation::cuda {

class GpuExecutionContext;

struct SolverBatchLayout {
    std::size_t batch_size;
    std::size_t dimension;

    SolverBatchLayout(std::size_t batch_size, std::size_t dimension);
    bool operator==(const SolverBatchLayout& other) const noexcept;
};

enum class SolverFailure {
    None = 0,
    NotInitialized,
    LayoutMismatch,
    InvalidArgument,
    NonFiniteInput,
    NonFiniteOutput,
    FactorizationFailed,
    ResidualTooLarge,
    UnsupportedMode,
};

struct SolverWorkspace {
    bool initialized = false;
    std::size_t allocation_count = 0;
    std::size_t batch_capacity = 0;
    std::size_t dimension_capacity = 0;
};

struct DeviceMatrixView {
    double* data = nullptr;
    std::size_t batch_size = 0;
    std::size_t dimension = 0;
};

struct DeviceVectorView {
    double* data = nullptr;
    std::size_t batch_size = 0;
    std::size_t dimension = 0;
};

struct DeviceResidualView {
    double* values = nullptr;
    std::size_t count = 0;
};

struct SolverStatus {
    bool ok = false;
    SolverFailure failure = SolverFailure::None;
    std::string message;
    double max_residual = 0.0;
    int failed_batch = -1;
    int backend_info = 0;

    static SolverStatus success(double residual = 0.0);
    static SolverStatus failure_status(SolverFailure failure, std::string message);
};

class GpuSolver {
public:
    GpuSolver(SolverMode requested_mode, SolverBatchLayout layout);
    GpuSolver(GpuExecutionContext& context, SolverMode requested_mode, SolverBatchLayout layout);
    GpuSolver(GpuExecutionContext& context, const GpuExecutionPolicy& policy,
              SolverBatchLayout layout);
    ~GpuSolver();

    GpuSolver(const GpuSolver&) = delete;
    GpuSolver& operator=(const GpuSolver&) = delete;
    GpuSolver(GpuSolver&&) noexcept;
    GpuSolver& operator=(GpuSolver&&) noexcept;

    SolverMode requested_mode() const noexcept;
    SolverMode resolved_mode() const noexcept;
    const SolverBatchLayout& layout() const noexcept;
    const SolverWorkspace& workspace() const noexcept;

    SolverStatus initialize_workspace();
    SolverStatus initialize_workspace(const SolverBatchLayout& layout);

    SolverStatus solve(const double* matrix,
                       const double* rhs,
                       double* solution);
    SolverStatus solve(const SolverBatchLayout& layout,
                       const double* matrix,
                       const double* rhs,
                       double* solution);
    SolverStatus solve_batch(const double* matrices,
                             const double* rhs,
                             double* solutions);
    SolverStatus solve_device(const DeviceMatrixView& matrix,
                              const DeviceVectorView& rhs,
                              DeviceVectorView solution,
                              DeviceResidualView residual = {});
    SolverStatus validate_device_result(DeviceResidualView residual);

    SolverStatus check_residual(const double* matrix,
                                const double* rhs,
                                const double* solution) const;

private:
    SolverStatus solve_at(std::size_t batch,
                          const double* matrix,
                          const double* rhs,
                          double* solution);

    struct Impl;
    Impl* impl_ = nullptr;
};

} // namespace coilgun::simulation::cuda
