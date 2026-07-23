#include <doctest/doctest.h>

#include "coilgun/coilgun_cuda.hpp"
#include "gpu_engine_fixture.hpp"

TEST_CASE("GPU assembly contract exposes fixed row-major layout") {
    coilgun::simulation::cuda::GpuStateLayout layout(2, 2, 3);
    CHECK(layout.system_matrix_size() == 2 * 5 * 5);
    CHECK(layout.m1(1, 1, 2) == 11);
}

TEST_CASE("device assembly matches host assembly before solve") {
    using namespace coilgun::simulation::cuda;
    if (!cuda_device_available()) return;
    auto geometry = gpu_test::geometry(2, 2, false);
    geometry.stage_resistances = {0.1, 0.2};
    geometry.stage_inductances = {1.0, 1.5};
    geometry.filament_resistances = {0.3, 0.4};
    geometry.filament_inductances = {0.8, 0.9};
    auto state = gpu_test::state(1, 2, 2, false);
    state.currents = {1.0, 2.0, 3.0, 4.0};
    state.m1 = {0.01, 0.02, 0.03, 0.04};
    state.dm1 = {0.1, 0.2, 0.3, 0.4};
    state.stage_voltages = {10.0, 20.0};
    state.velocity = {2.0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Direct;
    config.solver = SolverMode::Batched;
    GpuEngine engine(std::move(geometry), std::move(state), config);
    const auto host = engine.assemble_reference_for_test();
    const auto device = engine.assemble_device_for_test();
    REQUIRE(host.matrix.size() == device.matrix.size());
    REQUIRE(host.rhs.size() == device.rhs.size());
    for (std::size_t index = 0; index < host.matrix.size(); ++index)
        CHECK(device.matrix[index] == doctest::Approx(host.matrix[index]).epsilon(1e-12));
    for (std::size_t index = 0; index < host.rhs.size(); ++index)
        CHECK(device.rhs[index] == doctest::Approx(host.rhs[index]).epsilon(1e-12));
}
