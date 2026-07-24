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
