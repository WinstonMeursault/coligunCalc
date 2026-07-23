#include <doctest/doctest.h>

#include "coilgun/coilgun_cuda.hpp"
#include "gpu_engine_fixture.hpp"

TEST_CASE("resident pipeline layout keeps logical batch slots") {
    coilgun::simulation::cuda::GpuStateLayout layout(4, 2, 3);
    CHECK(layout.currents_size() == 20);
    CHECK(layout.active_mask_size() == 4);
}

TEST_CASE("resident engine reuses allocations across steps and reset") {
    using namespace coilgun::simulation::cuda;
    if (!cuda_device_available()) return;
    auto geometry = gpu_test::geometry(1, 1, false);
    geometry.stage_resistances = {0.1};
    geometry.stage_inductances = {1.0};
    geometry.filament_resistances = {0.2};
    geometry.filament_inductances = {0.8};
    auto state = gpu_test::state(1, 1, 1, false);
    state.stage_voltages = {10.0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Direct;
    config.solver = SolverMode::Batched;
    GpuEngine engine(std::move(geometry), std::move(state), config);
    const auto addresses = engine.device_buffer_addresses();
    const auto allocations = engine.device_allocation_count();
    engine.step();
    engine.step();
    engine.reset();
    CHECK(engine.device_buffer_addresses() == addresses);
    CHECK(engine.device_allocation_count() == allocations);
}
