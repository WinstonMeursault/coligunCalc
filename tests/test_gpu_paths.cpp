#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_state_kernels.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "coilgun/simulation/cuda/persistent_kernel.cuh"
#include "gpu_engine_fixture.hpp"

#include <cuda_runtime_api.h>

#include <array>
#include <algorithm>
#include <cstdint>
#include <vector>

using namespace coilgun::simulation::cuda;

namespace {

bool require_cuda() {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) {
        MESSAGE("CUDA device unavailable; skipping state-kernel test");
        return false;
    }
    return true;
}

template <typename T>
T* device_copy(const T* host, std::size_t count) {
    T* device = nullptr;
    REQUIRE(cudaMalloc(reinterpret_cast<void**>(&device), count * sizeof(T)) == cudaSuccess);
    REQUIRE(cudaMemcpy(device, host, count * sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess);
    return device;
}

} // namespace

TEST_CASE("engine exposes the unified pipeline order for B=1 and B>1") {
    using namespace coilgun::simulation::cuda;
    for (const std::size_t batch : {std::size_t{1}, std::size_t{2}}) {
        auto geometry = gpu_test::geometry();
        auto state = gpu_test::state(batch, 1, 1);
        state.active_mask.assign(batch, 0);
        GpuExecutionConfig config;
        config.backend = BackendMode::Fallback;
        GpuEngine engine(std::move(geometry), std::move(state), config);
        engine.step();
        CHECK(engine.pipeline_order() == std::vector<PipelineStage>{
            PipelineStage::Mutual, PipelineStage::Matrix, PipelineStage::Solver,
            PipelineStage::Force, PipelineStage::State});
    }
}

TEST_CASE("disabled thermal mode is absent from the engine pipeline") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    state.active_mask = {0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.thermal = ThermalMode::Disabled;
    GpuEngine engine(std::move(geometry), std::move(state), config);
    engine.step();
    CHECK(std::find(engine.pipeline_order().begin(), engine.pipeline_order().end(),
                    PipelineStage::Thermal) == engine.pipeline_order().end());
}

TEST_CASE("engine reuses persistent CUDA pipeline allocations across steps") {
    if (!require_cuda()) return;
    auto geometry = gpu_test::geometry(2, 2, true);
    auto state = gpu_test::state(1, 2, 2, true);
    GpuExecutionConfig config;
    config.backend = BackendMode::Auto;
    config.solver = SolverMode::Batched;
    config.thermal = ThermalMode::Gpu;
    GpuEngine engine(std::move(geometry), std::move(state), config);

    const auto before = engine.device_buffer_addresses();
    const auto allocations = engine.device_allocation_count();
    REQUIRE(allocations > 0);
    engine.step();
    const auto after_first = engine.device_buffer_addresses();
    engine.step();
    const auto after_second = engine.device_buffer_addresses();

    CHECK(after_first == before);
    CHECK(after_second == before);
    CHECK(engine.device_allocation_count() == allocations);
}

TEST_CASE("state kernels preserve CPU Euler ordering for B=1") {
    if (!require_cuda()) return;

    constexpr std::size_t B = 1;
    constexpr std::size_t S = 1;
    constexpr std::size_t F = 3;
    constexpr double dt = 0.25;
    constexpr double mass = 2.0;

    std::array<double, B * (S + F)> currents{2.0, 3.0, -1.0, 0.5};
    const std::array<double, B * S * F> dm1{4.0, -2.0, 1.5};
    const std::array<double, B * (S + F)> current_derivative{1.0, -2.0, 3.0, 4.0};
    std::array<double, B> acceleration{};
    std::array<double, B> velocity{5.0};
    std::array<double, B> position{7.0};
    double* d_currents = device_copy(currents.data(), currents.size());
    double* d_dm1 = device_copy(dm1.data(), dm1.size());
    double* d_current_derivative = device_copy(current_derivative.data(), current_derivative.size());
    double* d_acceleration = device_copy(acceleration.data(), B);
    double* d_velocity = device_copy(velocity.data(), B);
    double* d_position = device_copy(position.data(), B);
    double* d_force = device_copy(acceleration.data(), B);

    REQUIRE(launch_state_update(
                B, S, F, d_currents, d_current_derivative, d_dm1, nullptr,
                mass, dt, d_acceleration, d_velocity, d_position, d_force,
                StateKernelConfig{true, 32}) == cudaSuccess);
    REQUIRE(cudaDeviceSynchronize() == cudaSuccess);
    REQUIRE(cudaMemcpy(acceleration.data(), d_acceleration, sizeof(acceleration), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(velocity.data(), d_velocity, sizeof(velocity), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(position.data(), d_position, sizeof(position), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(currents.data(), d_currents, sizeof(currents), cudaMemcpyDeviceToHost) == cudaSuccess);
    cudaFree(d_currents); cudaFree(d_dm1); cudaFree(d_current_derivative); cudaFree(d_acceleration);
    cudaFree(d_velocity); cudaFree(d_position); cudaFree(d_force);

    CHECK(acceleration[0] == doctest::Approx(29.5 / mass));
    CHECK(velocity[0] == doctest::Approx(5.0 + dt * (29.5 / mass)));
    CHECK(position[0] == doctest::Approx(7.0 + dt * 5.0));
    CHECK(currents[0] == doctest::Approx(2.0 + dt));
    CHECK(currents[1] == doctest::Approx(3.0 - 2.0 * dt));
    CHECK(currents[2] == doctest::Approx(-1.0 + 3.0 * dt));
    CHECK(currents[3] == doctest::Approx(0.5 + 4.0 * dt));
}

TEST_CASE("state kernels reduce independent batch members deterministically") {
    if (!require_cuda()) return;

    constexpr std::size_t B = 2;
    constexpr std::size_t S = 2;
    constexpr std::size_t F = 2;
    constexpr double dt = 0.1;
    constexpr double mass = 4.0;

    std::array<double, B * (S + F)> currents{
        2.0, 0.0, 3.0, -1.0,
        -1.0, 4.0, 2.0, 5.0};
    const std::array<double, B * S * F> dm1{
        1.0, 2.0, -3.0, 4.0,
        0.5, -2.0, 1.0, 3.0};
    const std::array<std::uint8_t, B * S> trigger{1, 0, 1, 1};
    std::array<double, B> acceleration{};
    std::array<double, B> velocity{1.0, -2.0};
    std::array<double, B> position{10.0, 20.0};
    std::array<double, B> force{};
    double* d_currents = device_copy(currents.data(), currents.size());
    double* d_dm1 = device_copy(dm1.data(), dm1.size());
    auto* d_trigger = device_copy(trigger.data(), trigger.size());
    double* d_acceleration = device_copy(acceleration.data(), B);
    double* d_velocity = device_copy(velocity.data(), B);
    double* d_position = device_copy(position.data(), B);
    double* d_force = device_copy(force.data(), B);

    REQUIRE(launch_state_update(
                B, S, F, d_currents, nullptr, d_dm1, d_trigger,
                mass, dt, d_acceleration, d_velocity, d_position, d_force,
                StateKernelConfig{true, 32}) == cudaSuccess);
    REQUIRE(cudaDeviceSynchronize() == cudaSuccess);
    REQUIRE(cudaMemcpy(acceleration.data(), d_acceleration, sizeof(acceleration), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(velocity.data(), d_velocity, sizeof(velocity), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(position.data(), d_position, sizeof(position), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(force.data(), d_force, sizeof(force), cudaMemcpyDeviceToHost) == cudaSuccess);
    cudaFree(d_currents); cudaFree(d_dm1); cudaFree(d_trigger); cudaFree(d_acceleration);
    cudaFree(d_velocity); cudaFree(d_position); cudaFree(d_force);

    CHECK(acceleration[0] == doctest::Approx(2.0 * (3.0 * 1.0 + -1.0 * 2.0) / mass));
    CHECK(acceleration[1] == doctest::Approx(((-1.0 * (2.0 * 0.5 + 5.0 * -2.0)) +
                                               (4.0 * (2.0 * 1.0 + 5.0 * 3.0))) / mass));
    CHECK(velocity[0] == doctest::Approx(1.0 + dt * acceleration[0]));
    CHECK(velocity[1] == doctest::Approx(-2.0 + dt * acceleration[1]));
    CHECK(position[0] == doctest::Approx(10.0 + dt));
    CHECK(position[1] == doctest::Approx(20.0 - 2.0 * dt));
}

TEST_CASE("state kernels parallelize high-dimensional masked current updates") {
    if (!require_cuda()) return;

    constexpr std::size_t B = 2;
    constexpr std::size_t S = 1;
    constexpr std::size_t F = 32;
    constexpr std::size_t D = S + F;
    constexpr double dt = 0.125;
    constexpr double mass = 2.0;
    std::vector<double> currents(B * D), derivative(B * D), dm1(B * S * F, 0.0);
    for (std::size_t i = 0; i < currents.size(); ++i) {
        currents[i] = static_cast<double>(i + 1);
        derivative[i] = 0.25 * static_cast<double>(i + 1);
    }
    const auto initial_currents = currents;
    std::vector<double> acceleration(B, 0.0), velocity{2.0, -3.0}, position{4.0, 5.0};
    const std::array<std::uint8_t, B> active{1, 0};
    double* d_currents = device_copy(currents.data(), currents.size());
    double* d_derivative = device_copy(derivative.data(), derivative.size());
    double* d_dm1 = device_copy(dm1.data(), dm1.size());
    auto* d_active = device_copy(active.data(), active.size());
    double* d_acceleration = device_copy(acceleration.data(), acceleration.size());
    double* d_velocity = device_copy(velocity.data(), velocity.size());
    double* d_position = device_copy(position.data(), position.size());
    double* d_force = device_copy(acceleration.data(), acceleration.size());

    REQUIRE(launch_state_update_masked(
                B, S, F, d_currents, d_derivative, d_dm1, nullptr, d_active,
                mass, dt, d_acceleration, d_velocity, d_position, d_force,
                StateKernelConfig{true, 32}) == cudaSuccess);
    REQUIRE(cudaDeviceSynchronize() == cudaSuccess);
    REQUIRE(cudaMemcpy(currents.data(), d_currents, currents.size() * sizeof(double),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(velocity.data(), d_velocity, velocity.size() * sizeof(double),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(position.data(), d_position, position.size() * sizeof(double),
                       cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(cudaMemcpy(acceleration.data(), d_acceleration,
                       acceleration.size() * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess);

    for (std::size_t i = 0; i < D; ++i) {
        CHECK(currents[i] == doctest::Approx(initial_currents[i] + dt * derivative[i]));
        CHECK(currents[D + i] == doctest::Approx(initial_currents[D + i]));
    }
    CHECK(velocity[0] == doctest::Approx(2.0));
    CHECK(velocity[1] == doctest::Approx(-3.0));
    CHECK(position[0] == doctest::Approx(4.0 + dt * 2.0));
    CHECK(position[1] == doctest::Approx(5.0));
    CHECK(acceleration[0] == doctest::Approx(0.0));
    CHECK(acceleration[1] == doctest::Approx(0.0));

    cudaFree(d_currents); cudaFree(d_derivative); cudaFree(d_dm1); cudaFree(d_active);
    cudaFree(d_acceleration); cudaFree(d_velocity); cudaFree(d_position); cudaFree(d_force);
}

TEST_CASE("persistent buffers expose generation protocol and are safely reusable") {
    PersistentBuffers buffers;

    CHECK(buffers.generation == 0);
    CHECK(buffers.shutdown_ack == 0);
    CHECK(buffers.status == PersistentStatus::Uninitialized);

    CHECK(buffers.batch_id == nullptr);
    CHECK(buffers.d_batch_id == nullptr);
    CHECK(buffers.generation == 0);
    CHECK(buffers.shutdown_ack == 0);
    CHECK(buffers.status == PersistentStatus::Uninitialized);
}

TEST_CASE("requested graph backend is reported honestly and remains selected across steps") {
    if (!require_cuda()) return;

    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    GpuEngine engine(std::move(geometry), std::move(state), config);

    REQUIRE(engine.report().backend == BackendMode::Graph);
    engine.step();
    engine.step();
    CHECK(engine.report().backend == BackendMode::Graph);
    CHECK(engine.report().fallback_count == 0);
}

TEST_CASE("unavailable graph backend locks the engine to reported fallback") {
    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    GpuCapability capability;
    capability.supports_graph = false;
    GpuEngine engine(std::move(geometry), std::move(state), config, capability);

    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().static_fallback_reason == FallbackReason::CapabilityUnavailable);
    engine.step();
    engine.step();
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().fallback_count == 1);
}

TEST_CASE("requested persistent backend reports runtime protocol failure") {
    if (!require_cuda()) return;

    auto geometry = gpu_test::geometry(1, 2);
    auto state = gpu_test::state(1, 1, 2);
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;
    GpuCapability capability;
    capability.supports_persistent_control_stream = true;
    GpuEngine engine(std::move(geometry), std::move(state), config, capability);

    REQUIRE(engine.report().requested_backend == BackendMode::Persistent);
    REQUIRE(engine.report().backend == BackendMode::Fallback);
    REQUIRE_FALSE(engine.report().gpu_executed);
    REQUIRE_FALSE(engine.report().fallback_reason.empty());
    REQUIRE(engine.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    engine.step();
    engine.step();
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK(engine.report().fallback_count == 1);
}
