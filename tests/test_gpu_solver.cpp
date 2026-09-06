/**
 * @file test_gpu_solver.cpp
 * @brief Contract tests for the unified dense solver interface.
 */

#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_solver.hpp"
#include "coilgun/simulation/cuda/gpu_execution_context.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "gpu_engine_fixture.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include <cstdint>

#include <cuda_runtime_api.h>

using namespace coilgun::simulation::cuda;

namespace {
bool gpu_available() {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

template<typename T>
class DeviceAllocation {
public:
    explicit DeviceAllocation(std::size_t count) {
        if (cudaMalloc(reinterpret_cast<void**>(&pointer_), count * sizeof(T)) != cudaSuccess)
            throw std::runtime_error("CUDA test allocation failed");
    }
    ~DeviceAllocation() { (void)cudaFree(pointer_); }
    DeviceAllocation(const DeviceAllocation&) = delete;
    DeviceAllocation& operator=(const DeviceAllocation&) = delete;
    T* get() const noexcept { return pointer_; }

private:
    T* pointer_ = nullptr;
};
}

TEST_CASE("Eigen solver resolves mode and solves a known FP64 system") {
    const SolverBatchLayout layout{1, 3};
    GpuSolver solver(SolverMode::Eigen, layout);

    CHECK(solver.resolved_mode() == SolverMode::Eigen);
    CHECK(solver.initialize_workspace().ok);

    const std::array<double, 9> matrix{
        4.0, 1.0, 1.0,
        1.0, 3.0, 0.0,
        1.0, 0.0, 2.0,
    };
    const std::array<double, 3> rhs{8.0, 7.0, 5.0};
    std::array<double, 3> solution{};

    const SolverStatus status = solver.solve(matrix.data(), rhs.data(), solution.data());

    REQUIRE(status.ok);
    CHECK(solution[0] == doctest::Approx(1.0));
    CHECK(solution[1] == doctest::Approx(2.0));
    CHECK(solution[2] == doctest::Approx(2.0));
    CHECK(status.max_residual < 1e-12);
}

TEST_CASE("Eigen solver supports a row-major batch layout") {
    const SolverBatchLayout layout{2, 2};
    GpuSolver solver(SolverMode::Eigen, layout);
    REQUIRE(solver.initialize_workspace().ok);

    const std::array<double, 8> matrices{
        2.0, 1.0, 1.0, 2.0,
        3.0, 0.0, 0.0, 4.0,
    };
    const std::array<double, 4> rhs{5.0, 5.0, 6.0, 8.0};
    std::array<double, 4> solutions{};

    const SolverStatus status = solver.solve_batch(
        matrices.data(), rhs.data(), solutions.data());

    REQUIRE(status.ok);
    CHECK(solutions[0] == doctest::Approx(5.0 / 3.0));
    CHECK(solutions[1] == doctest::Approx(5.0 / 3.0));
    CHECK(solutions[2] == doctest::Approx(2.0));
    CHECK(solutions[3] == doctest::Approx(2.0));
    CHECK(status.max_residual < 1e-12);
}

TEST_CASE("Solver rejects a dimension or batch-layout mismatch") {
    GpuSolver solver(SolverMode::Eigen, SolverBatchLayout{2, 2});
    CHECK_FALSE(solver.initialize_workspace(SolverBatchLayout{1, 2}).ok);

    const std::array<double, 4> matrix{2.0, 0.0, 0.0, 2.0};
    const std::array<double, 2> rhs{2.0, 4.0};
    std::array<double, 2> solution{};

    const SolverStatus status = solver.solve(
        SolverBatchLayout{1, 2}, matrix.data(), rhs.data(), solution.data());
    CHECK_FALSE(status.ok);
    CHECK(status.failure == SolverFailure::LayoutMismatch);
}

TEST_CASE("Solver batch layout rejects size_t overflow and CUDA int overflow") {
    CHECK_THROWS_AS(
        SolverBatchLayout(std::numeric_limits<std::size_t>::max(), 2),
        std::invalid_argument);
    CHECK_THROWS_AS(
        SolverBatchLayout(2, std::numeric_limits<std::size_t>::max()),
        std::invalid_argument);
    CHECK_THROWS_AS(
        SolverBatchLayout(static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1, 1),
        std::invalid_argument);
    CHECK_THROWS_AS(
        SolverBatchLayout(1, static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1),
        std::invalid_argument);
}

TEST_CASE("Solver reports non-finite input and output as failures") {
    GpuSolver solver(SolverMode::Eigen, SolverBatchLayout{1, 2});
    REQUIRE(solver.initialize_workspace().ok);

    const std::array<double, 4> matrix{
        2.0, 0.0,
        0.0, 2.0,
    };
    const std::array<double, 2> finite_rhs{2.0, 4.0};
    const std::array<double, 2> nonfinite_rhs{
        std::numeric_limits<double>::quiet_NaN(), 4.0};
    std::array<double, 2> solution{};

    const SolverStatus input_status = solver.solve(
        matrix.data(), nonfinite_rhs.data(), solution.data());
    CHECK_FALSE(input_status.ok);
    CHECK(input_status.failure == SolverFailure::NonFiniteInput);

    const SolverStatus output_status = solver.solve(
        matrix.data(), finite_rhs.data(), solution.data());
    REQUIRE(output_status.ok);
    solution[1] = std::numeric_limits<double>::infinity();
    const SolverStatus sanity = solver.check_residual(
        matrix.data(), finite_rhs.data(), solution.data());
    CHECK_FALSE(sanity.ok);
    CHECK(sanity.failure == SolverFailure::NonFiniteOutput);
}

TEST_CASE("Solver exposes residual and sanity checks") {
    GpuSolver solver(SolverMode::Eigen, SolverBatchLayout{1, 2});
    REQUIRE(solver.initialize_workspace().ok);

    const std::array<double, 4> matrix{2.0, 0.0, 0.0, 4.0};
    const std::array<double, 2> rhs{2.0, 8.0};
    const std::array<double, 2> solution{1.0, 2.0};

    const SolverStatus status = solver.check_residual(
        matrix.data(), rhs.data(), solution.data());
    CHECK(status.ok);
    CHECK(status.max_residual == doctest::Approx(0.0));
}

TEST_CASE("Solver workspace is initialized once and reused across calls") {
    GpuSolver solver(SolverMode::Eigen, SolverBatchLayout{1, 2});
    CHECK(solver.workspace().allocation_count == 0);
    REQUIRE(solver.initialize_workspace().ok);
    const auto allocations = solver.workspace().allocation_count;

    const std::array<double, 4> matrix{2.0, 0.0, 0.0, 2.0};
    const std::array<double, 2> rhs{2.0, 4.0};
    std::array<double, 2> solution{};
    REQUIRE(solver.solve(matrix.data(), rhs.data(), solution.data()).ok);
    REQUIRE(solver.solve(matrix.data(), rhs.data(), solution.data()).ok);

    CHECK(solver.workspace().allocation_count == allocations);
    CHECK(solver.workspace().initialized);
}

TEST_CASE("Moved-from solver is safe to query and use") {
    GpuSolver source(SolverMode::Eigen, SolverBatchLayout{1, 2});
    GpuSolver moved(std::move(source));
    CHECK(source.resolved_mode() == SolverMode::Auto);
    CHECK(source.workspace().initialized == false);
    CHECK_FALSE(source.initialize_workspace().ok);
}

TEST_CASE("CUDA solver performs FP64 column-major batched LU without solve allocations") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping batched LU test");
        return;
    }

    GpuExecutionContext context;
    const GpuExecutionPolicy policy = GpuExecutionPlanner::plan(
        1, 1, 2, false, GpuCapability{}, GpuExecutionConfig{BackendMode::Fallback,
                                                               SolverMode::Batched,
                                                               PrecisionMode::Aggressive,
                                                               ThermalMode::Disabled});
    GpuSolver solver(context, policy, SolverBatchLayout{2, 2});
    REQUIRE(solver.resolved_mode() == SolverMode::Batched);
    REQUIRE(solver.initialize_workspace().ok);
    const auto allocations = solver.workspace().allocation_count;

    const std::array<double, 8> matrices{2.0, 1.0, 1.0, 2.0,
                                         3.0, 0.0, 0.0, 4.0};
    const std::array<double, 4> rhs{5.0, 5.0, 6.0, 8.0};
    std::array<double, 4> solutions{};
    const SolverStatus status = solver.solve_batch(matrices.data(), rhs.data(), solutions.data());

    REQUIRE(status.ok);
    CHECK(solutions[0] == doctest::Approx(5.0 / 3.0));
    CHECK(solutions[1] == doctest::Approx(5.0 / 3.0));
    CHECK(solutions[2] == doctest::Approx(2.0));
    CHECK(solutions[3] == doctest::Approx(2.0));
    CHECK(status.backend_info == 0);
    CHECK(status.failed_batch == -1);
    CHECK(solver.workspace().allocation_count == allocations);
}

TEST_CASE("CUDA solver reports the failed batch and backend info") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping batched LU diagnostics test");
        return;
    }

    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched, SolverBatchLayout{2, 2});
    REQUIRE(solver.initialize_workspace().ok);
    const std::array<double, 8> matrices{2.0, 0.0, 0.0, 2.0,
                                         0.0, 0.0, 0.0, 0.0};
    const std::array<double, 4> rhs{2.0, 4.0, 1.0, 1.0};
    std::array<double, 4> solutions{};
    const SolverStatus status = solver.solve_batch(matrices.data(), rhs.data(), solutions.data());

    CHECK_FALSE(status.ok);
    CHECK(status.failure == SolverFailure::FactorizationFailed);
    CHECK(status.failed_batch == 1);
    CHECK(status.backend_info != 0);
    CHECK(status.message.find("batch 1") != std::string::npos);
}

TEST_CASE("CUDA device solver preserves assembly input and reports residual") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping device solver test");
        return;
    }

    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched, SolverBatchLayout{1, 2});
    REQUIRE(solver.initialize_workspace().ok);
    const std::array<double, 4> matrix{2.0, 1.0, 1.0, 2.0};
    const std::array<double, 2> rhs{5.0, 5.0};
    std::array<double, 4> matrix_after{};
    std::array<double, 2> solution{};
    double residual = -1.0;
    DeviceAllocation<double> device_matrix(matrix.size());
    DeviceAllocation<double> device_rhs(rhs.size());
    DeviceAllocation<double> device_solution(solution.size());
    DeviceAllocation<double> device_residual(1);
    REQUIRE(cudaMemcpy(device_matrix.get(), matrix.data(), sizeof(matrix),
                       cudaMemcpyHostToDevice) == cudaSuccess);
    REQUIRE(cudaMemcpy(device_rhs.get(), rhs.data(), sizeof(rhs),
                       cudaMemcpyHostToDevice) == cudaSuccess);

    const SolverStatus status = solver.solve_device(
        DeviceMatrixView{device_matrix.get(), 1, 2},
        DeviceVectorView{device_rhs.get(), 1, 2},
        DeviceVectorView{device_solution.get(), 1, 2},
        DeviceResidualView{device_residual.get(), 1});

    REQUIRE(status.ok);
    REQUIRE(cudaMemcpy(matrix_after.data(), device_matrix.get(), sizeof(matrix_after),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(solution.data(), device_solution.get(), sizeof(solution),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(&residual, device_residual.get(), sizeof(residual),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    CHECK(matrix_after == matrix);
    CHECK(solution[0] == doctest::Approx(5.0 / 3.0));
    CHECK(solution[1] == doctest::Approx(5.0 / 3.0));
    CHECK(residual < 1.0e-12);
}

TEST_CASE("CUDA device residual preserves large-dimension validation") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping large-dimension residual test");
        return;
    }

    constexpr std::size_t batch_size = 2;
    constexpr std::size_t dimension = 129;
    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched,
                     SolverBatchLayout{batch_size, dimension});
    REQUIRE(solver.initialize_workspace().ok);

    std::vector<double> matrices(batch_size * dimension * dimension, 0.0);
    std::vector<double> rhs(batch_size * dimension, 2.0);
    for (std::size_t batch = 0; batch < batch_size; ++batch)
        for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal)
            matrices[batch * dimension * dimension + diagonal * dimension + diagonal] = 2.0;

    DeviceAllocation<double> device_matrices(matrices.size());
    DeviceAllocation<double> device_rhs(rhs.size());
    DeviceAllocation<double> device_solution(rhs.size());
    DeviceAllocation<double> device_residual(batch_size);
    REQUIRE(cudaMemcpyAsync(device_matrices.get(), matrices.data(),
                            matrices.size() * sizeof(double), cudaMemcpyHostToDevice,
                            context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(device_rhs.get(), rhs.data(),
                            rhs.size() * sizeof(double), cudaMemcpyHostToDevice,
                            context.stream()) == cudaSuccess);

    REQUIRE(solver.solve_device(
        DeviceMatrixView{device_matrices.get(), batch_size, dimension},
        DeviceVectorView{device_rhs.get(), batch_size, dimension},
        DeviceVectorView{device_solution.get(), batch_size, dimension},
        DeviceResidualView{device_residual.get(), batch_size}).ok);
    REQUIRE(solver.validate_device_result(
                DeviceResidualView{device_residual.get(), batch_size}).ok);

    std::vector<double> residual(batch_size, -1.0);
    REQUIRE(cudaMemcpy(residual.data(), device_residual.get(),
                       residual.size() * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess);
    for (const double value : residual) CHECK(value < 1.0e-12);
}

TEST_CASE("CUDA device solver defers singular-system failure to validation") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping deferred device status test");
        return;
    }

    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched, SolverBatchLayout{2, 2});
    REQUIRE(solver.initialize_workspace().ok);
    const std::array<double, 8> matrices{
        2.0, 0.0, 0.0, 2.0,
        0.0, 0.0, 0.0, 0.0};
    const std::array<double, 4> rhs{2.0, 4.0, 1.0, 1.0};
    DeviceAllocation<double> device_matrices(matrices.size());
    DeviceAllocation<double> device_rhs(rhs.size());
    DeviceAllocation<double> device_solution(rhs.size());
    DeviceAllocation<double> device_residual(2);
    REQUIRE(cudaMemcpy(device_matrices.get(), matrices.data(), sizeof(matrices),
                       cudaMemcpyHostToDevice) == cudaSuccess);
    REQUIRE(cudaMemcpy(device_rhs.get(), rhs.data(), sizeof(rhs),
                       cudaMemcpyHostToDevice) == cudaSuccess);

    const SolverStatus enqueue_status = solver.solve_device(
        DeviceMatrixView{device_matrices.get(), 2, 2},
        DeviceVectorView{device_rhs.get(), 2, 2},
        DeviceVectorView{device_solution.get(), 2, 2},
        DeviceResidualView{device_residual.get(), 2});

    REQUIRE(enqueue_status.ok);
    const SolverStatus validation_status = solver.validate_device_result(
        DeviceResidualView{device_residual.get(), 2});
    CHECK_FALSE(validation_status.ok);
    CHECK(validation_status.failure == SolverFailure::FactorizationFailed);
    CHECK(validation_status.failed_batch == 1);
    CHECK(validation_status.backend_info != 0);
}

TEST_CASE("CUDA device solver reuses a stable output view across repeated solves") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping repeated device solver test");
        return;
    }

    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched, SolverBatchLayout{2, 2});
    REQUIRE(solver.initialize_workspace().ok);
    const std::array<double, 8> matrix{
        2.0, 1.0, 1.0, 2.0,
        3.0, 1.0, 1.0, 3.0};
    const std::array<double, 4> rhs{5.0, 5.0, 8.0, 8.0};
    std::array<double, 8> matrix_after{};
    std::array<double, 4> solution{};
    std::array<double, 2> residual{};
    DeviceAllocation<double> device_matrix(matrix.size());
    DeviceAllocation<double> device_rhs(rhs.size());
    DeviceAllocation<double> device_solution(solution.size());
    DeviceAllocation<double> device_residual(residual.size());
    REQUIRE(cudaMemcpy(device_matrix.get(), matrix.data(), sizeof(matrix),
                       cudaMemcpyHostToDevice) == cudaSuccess);
    REQUIRE(cudaMemcpy(device_rhs.get(), rhs.data(), sizeof(rhs),
                       cudaMemcpyHostToDevice) == cudaSuccess);

    for (int repeat = 0; repeat < 3; ++repeat) {
        const SolverStatus status = solver.solve_device(
            DeviceMatrixView{device_matrix.get(), 2, 2},
            DeviceVectorView{device_rhs.get(), 2, 2},
            DeviceVectorView{device_solution.get(), 2, 2},
            DeviceResidualView{device_residual.get(), 2});
        REQUIRE(status.ok);
        REQUIRE(solver.validate_device_result(
                    DeviceResidualView{device_residual.get(), 2}).ok);
    }

    REQUIRE(cudaMemcpy(matrix_after.data(), device_matrix.get(), sizeof(matrix_after),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(solution.data(), device_solution.get(), sizeof(solution),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(residual.data(), device_residual.get(), sizeof(residual),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    CHECK(matrix_after == matrix);
    CHECK(solution[0] == doctest::Approx(5.0 / 3.0));
    CHECK(solution[1] == doctest::Approx(5.0 / 3.0));
    CHECK(solution[2] == doctest::Approx(2.0));
    CHECK(solution[3] == doctest::Approx(2.0));
    CHECK(residual[0] < 1.0e-12);
    CHECK(residual[1] < 1.0e-12);
}

TEST_CASE("CUDA context operations preserve the caller's current device") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping context device guard test");
        return;
    }

    int before = -1;
    REQUIRE(cudaGetDevice(&before) == cudaSuccess);
    {
        GpuExecutionContext context({before, cudaStreamNonBlocking, 4096});
        context.record_start();
        context.record_stop();
        context.synchronize();
        int during = -1;
        REQUIRE(cudaGetDevice(&during) == cudaSuccess);
        CHECK(during == before);
    }
    int after = -1;
    REQUIRE(cudaGetDevice(&after) == cudaSuccess);
    CHECK(after == before);
}

TEST_CASE("CUDA batched solver reuses host staging after initialization") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping staging reuse test");
        return;
    }

    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched, SolverBatchLayout{2, 2});
    REQUIRE(solver.initialize_workspace().ok);
    const auto allocations = solver.workspace().allocation_count;
    const std::array<double, 8> matrices{2.0, 1.0, 1.0, 2.0,
                                         3.0, 0.0, 0.0, 4.0};
    const std::array<double, 4> rhs{5.0, 5.0, 6.0, 8.0};
    std::array<double, 4> solutions{};
    REQUIRE(solver.solve_batch(matrices.data(), rhs.data(), solutions.data()).ok);
    REQUIRE(solver.solve_batch(matrices.data(), rhs.data(), solutions.data()).ok);
    CHECK(solver.workspace().allocation_count == allocations);
}

TEST_CASE("CUDA device solver skips inactive identity rows") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping inactive-row solver test");
        return;
    }

    constexpr std::size_t batch_size = 128;
    constexpr std::size_t dimension = 32;
    constexpr int warmup = 3;
    constexpr int samples = 8;
    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched,
                     SolverBatchLayout{batch_size, dimension});
    REQUIRE(solver.initialize_workspace().ok);

    std::vector<double> active_matrix(batch_size * dimension * dimension, 0.0);
    std::vector<double> inactive_matrix(batch_size * dimension * dimension, 0.0);
    std::vector<double> rhs(batch_size * dimension, 0.0);
    std::vector<std::uint8_t> all_active(batch_size, 1);
    std::vector<std::uint8_t> all_inactive(batch_size, 0);
    for (std::size_t batch = 0; batch < batch_size; ++batch) {
        for (std::size_t diagonal = 0; diagonal < dimension; ++diagonal) {
            active_matrix[batch * dimension * dimension + diagonal * dimension + diagonal] = 2.0;
            inactive_matrix[batch * dimension * dimension + diagonal * dimension + diagonal] = 1.0;
        }
    }
    DeviceAllocation<double> d_active(active_matrix.size());
    DeviceAllocation<double> d_inactive(inactive_matrix.size());
    DeviceAllocation<double> d_rhs(rhs.size());
    DeviceAllocation<double> d_solution(rhs.size());
    DeviceAllocation<double> d_residual(batch_size);
    DeviceAllocation<std::uint8_t> d_all_active(batch_size);
    DeviceAllocation<std::uint8_t> d_all_inactive(batch_size);
    REQUIRE(cudaMemcpyAsync(d_active.get(), active_matrix.data(), sizeof(double) * active_matrix.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_inactive.get(), inactive_matrix.data(), sizeof(double) * inactive_matrix.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_rhs.get(), rhs.data(), sizeof(double) * rhs.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_all_active.get(), all_active.data(), all_active.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_all_inactive.get(), all_inactive.data(), all_inactive.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);

    const auto elapsed = [&](double* matrix, const std::uint8_t* active_mask,
                             std::size_t active_count) {
        for (int repeat = 0; repeat < warmup; ++repeat) {
            REQUIRE(solver.solve_device(
                DeviceMatrixView{matrix, batch_size, dimension, active_mask, active_count},
                DeviceVectorView{d_rhs.get(), batch_size, dimension},
                DeviceVectorView{d_solution.get(), batch_size, dimension},
                DeviceResidualView{d_residual.get(), batch_size}).ok);
            REQUIRE(solver.validate_device_result(
                DeviceResidualView{d_residual.get(), batch_size}).ok);
        }
        const auto start = std::chrono::steady_clock::now();
        for (int repeat = 0; repeat < samples; ++repeat) {
            REQUIRE(solver.solve_device(
                DeviceMatrixView{matrix, batch_size, dimension, active_mask, active_count},
                DeviceVectorView{d_rhs.get(), batch_size, dimension},
                DeviceVectorView{d_solution.get(), batch_size, dimension},
                DeviceResidualView{d_residual.get(), batch_size}).ok);
            REQUIRE(solver.validate_device_result(
                DeviceResidualView{d_residual.get(), batch_size}).ok);
        }
        return std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
    };

    const double active_ms = elapsed(d_active.get(), d_all_active.get(), batch_size);
    const double inactive_ms = elapsed(d_inactive.get(), d_all_inactive.get(), 0);
    CHECK(inactive_ms < active_ms * 0.80);

    std::vector<std::uint8_t> mixed_mask(batch_size, 0);
    std::vector<double> mixed_rhs(rhs.size(), 0.0);
    for (std::size_t batch = 0; batch < batch_size; batch += 2) {
        mixed_mask[batch] = 1;
        for (std::size_t row = 0; row < dimension; ++row)
            mixed_rhs[batch * dimension + row] = 2.0;
    }
    DeviceAllocation<std::uint8_t> d_mixed_mask(batch_size);
    DeviceAllocation<double> d_mixed_rhs(mixed_rhs.size());
    REQUIRE(cudaMemcpyAsync(d_mixed_mask.get(), mixed_mask.data(), mixed_mask.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_mixed_rhs.get(), mixed_rhs.data(),
                            sizeof(double) * mixed_rhs.size(), cudaMemcpyHostToDevice,
                            context.stream()) == cudaSuccess);
    REQUIRE(solver.solve_device(
        DeviceMatrixView{d_active.get(), batch_size, dimension, d_mixed_mask.get(), batch_size / 2},
        DeviceVectorView{d_mixed_rhs.get(), batch_size, dimension},
        DeviceVectorView{d_solution.get(), batch_size, dimension},
        DeviceResidualView{d_residual.get(), batch_size}).ok);
    REQUIRE(solver.validate_device_result(
        DeviceResidualView{d_residual.get(), batch_size}).ok);
    std::vector<double> mixed_solution(mixed_rhs.size(), -1.0);
    REQUIRE(cudaMemcpy(mixed_solution.data(), d_solution.get(),
                       sizeof(double) * mixed_solution.size(), cudaMemcpyDeviceToHost) == cudaSuccess);
    for (std::size_t batch = 0; batch < batch_size; ++batch) {
        for (std::size_t row = 0; row < dimension; ++row) {
            const double expected = mixed_mask[batch] != 0 ? 1.0 : 0.0;
            CHECK(mixed_solution[batch * dimension + row] == doctest::Approx(expected));
        }
    }

    for (const std::size_t ratio_active_count :
         {std::size_t{0}, batch_size / 4, batch_size / 2,
          (batch_size * 3) / 4, batch_size}) {
        std::fill(mixed_mask.begin(), mixed_mask.end(), 0);
        std::fill(mixed_rhs.begin(), mixed_rhs.end(), 0.0);
        for (std::size_t batch = 0; batch < ratio_active_count; ++batch) {
            mixed_mask[batch] = 1;
            for (std::size_t row = 0; row < dimension; ++row)
                mixed_rhs[batch * dimension + row] = 2.0;
        }
        REQUIRE(cudaMemcpyAsync(d_mixed_mask.get(), mixed_mask.data(), mixed_mask.size(),
                                cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
        REQUIRE(cudaMemcpyAsync(d_mixed_rhs.get(), mixed_rhs.data(),
                                sizeof(double) * mixed_rhs.size(), cudaMemcpyHostToDevice,
                                context.stream()) == cudaSuccess);
        REQUIRE(solver.solve_device(
            DeviceMatrixView{d_active.get(), batch_size, dimension, d_mixed_mask.get(),
                             ratio_active_count},
            DeviceVectorView{d_mixed_rhs.get(), batch_size, dimension},
            DeviceVectorView{d_solution.get(), batch_size, dimension},
            DeviceResidualView{d_residual.get(), batch_size}).ok);
        REQUIRE(solver.validate_device_result(
            DeviceResidualView{d_residual.get(), batch_size}).ok);
        REQUIRE(cudaMemcpy(mixed_solution.data(), d_solution.get(),
                           sizeof(double) * mixed_solution.size(), cudaMemcpyDeviceToHost) == cudaSuccess);
        for (std::size_t batch = 0; batch < batch_size; ++batch) {
            const double expected = batch < ratio_active_count ? 1.0 : 0.0;
            for (std::size_t row = 0; row < dimension; ++row)
                CHECK(mixed_solution[batch * dimension + row] == doctest::Approx(expected));
        }
    }
}

TEST_CASE("CUDA device solver ignores residuals for inactive rows") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping inactive residual test");
        return;
    }

    constexpr std::size_t batch_size = 2;
    constexpr std::size_t dimension = 2;
    GpuExecutionContext context;
    GpuSolver solver(context, SolverMode::Batched,
                     SolverBatchLayout{batch_size, dimension});
    REQUIRE(solver.initialize_workspace().ok);

    const std::array<double, 8> matrices{1.0, 0.0, 0.0, 1.0,
                                         1.0, 0.0, 0.0, 1.0};
    const std::array<double, 4> rhs{2.0, 3.0, 11.0, 13.0};
    const std::array<std::uint8_t, 2> active{1, 0};
    DeviceAllocation<double> d_matrix(matrices.size());
    DeviceAllocation<double> d_rhs(rhs.size());
    DeviceAllocation<double> d_solution(rhs.size());
    DeviceAllocation<double> d_residual(batch_size);
    DeviceAllocation<std::uint8_t> d_active(active.size());
    REQUIRE(cudaMemcpyAsync(d_matrix.get(), matrices.data(), sizeof(matrices),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_rhs.get(), rhs.data(), sizeof(rhs),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);
    REQUIRE(cudaMemcpyAsync(d_active.get(), active.data(), active.size(),
                            cudaMemcpyHostToDevice, context.stream()) == cudaSuccess);

    REQUIRE(solver.solve_device(
        DeviceMatrixView{d_matrix.get(), batch_size, dimension, d_active.get(), 1},
        DeviceVectorView{d_rhs.get(), batch_size, dimension},
        DeviceVectorView{d_solution.get(), batch_size, dimension},
        DeviceResidualView{d_residual.get(), batch_size}).ok);
    CHECK(solver.validate_device_result(
        DeviceResidualView{d_residual.get(), batch_size}).ok);
}

TEST_CASE("Engine exposes one resolved policy and one calibration report") {
    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    state.active_mask = {0};
    GpuExecutionConfig config;
    config.solver = SolverMode::Batched;
    config.enable_calibration = true;
    GpuEngine engine(std::move(geometry), std::move(state), config);

    const bool has_gpu = gpu_available();
    CHECK(engine.policy().solver == (has_gpu ? SolverMode::Batched : SolverMode::Eigen));
    CHECK(engine.graph_variant().solver == engine.policy().solver);
    CHECK(engine.report().solver == engine.policy().solver);
    CHECK(engine.report().calibrated == has_gpu);
    CHECK(engine.calibration_count() == (has_gpu ? 1u : 0u));

    engine.step();
    CHECK(engine.report().calibrated == has_gpu);
    CHECK(engine.calibration_count() == (has_gpu ? 1u : 0u));
    engine.step();
    CHECK(engine.calibration_count() == (has_gpu ? 1u : 0u));
}

TEST_CASE("Engine reports fallback when graph execution is not implemented") {
    if (!gpu_available()) {
        MESSAGE("CUDA device unavailable; skipping engine resource test");
        return;
    }

    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    state.active_mask = {0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    config.solver = SolverMode::Eigen;
    GpuCapability capability;
    GpuEngine engine(std::move(geometry), std::move(state), config, capability);

    CHECK_FALSE(engine.context_available());
    CHECK(engine.solver_workspace_initialized());
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().static_fallback_reason == FallbackReason::MetadataConflict);
    CHECK(engine.report().solver == SolverMode::Eigen);
}

TEST_CASE("Engine executes the solver contract for CPU Eigen policy") {
    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    state.currents = {2.0, 3.0};
    state.active_mask = {0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;
    GpuEngine engine(std::move(geometry), std::move(state), config);

    engine.step();
    CHECK(engine.result().completed_steps == 1);
    CHECK(engine.report().solver == SolverMode::Eigen);
    CHECK(engine.report().solver_time_ms >= 0.0);
}
