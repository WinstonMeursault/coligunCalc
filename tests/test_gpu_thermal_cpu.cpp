#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "gpu_engine_fixture.hpp"

using namespace coilgun::simulation::cuda;

TEST_CASE("CPU-only engine thermal contract updates resistance from Euler currents") {
    GpuGeometryInput geometry;
    geometry.n_stages = 1;
    geometry.n_filaments = 1;
    geometry.stage_geometry = {0.05};
    geometry.filament_geometry = {0.03};
    geometry.stage_inner_radii = {0.045};
    geometry.stage_outer_radii = {0.055};
    geometry.stage_lengths = {0.01};
    geometry.stage_turns = {20};
    geometry.stage_positions = {0.0};
    geometry.filament_inner_radii = {0.029};
    geometry.filament_outer_radii = {0.031};
    geometry.filament_lengths = {0.01};
    geometry.filament_positions = {0.04};
    geometry.stage_resistances = {1.0};
    geometry.stage_inductances = {2.0};
    geometry.filament_resistances = {0.5};
    geometry.filament_inductances = {1.5};
    geometry.stage_mutual_inductances = {0.0};
    geometry.filament_mutual_inductances = {0.0};
    geometry.thermal_enabled = true;

    GpuEngineState state;
    state.currents = {0.0, 10.0};
    state.m1 = {0.0};
    state.dm1 = {0.0};
    state.stage_voltages = {0.0};
    state.active_mask = {1};
    state.trigger_mask = {1};
    state.velocity = {0.0};
    state.position = {0.0};
    state.temperatures = {293.0};
    state.filament_masses = {0.01};
    state.reference_resistances = {0.5};
    state.filament_materials = {1};
    state.material_density = 1.0;
    state.dt = 1.0e-3;
    state.mass = 1.0;

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;
    config.thermal = ThermalMode::Cpu;

    GpuEngine engine(std::move(geometry), std::move(state), config);
    engine.step();

    REQUIRE(engine.pipeline_order().back() == PipelineStage::Thermal);
    CHECK(engine.state().temperatures[0] > 293.0);
    CHECK(engine.state().resistances[0] > 0.5);
}

TEST_CASE("CPU fallback normalizes a GPU thermal request to CPU") {
    auto geometry = gpu_test::geometry(1, 1, true);
    auto state = gpu_test::state(1, 1, 1, true);
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.thermal = ThermalMode::Gpu;

    GpuEngine engine(std::move(geometry), std::move(state), config);
    CHECK(engine.report().thermal == ThermalMode::Cpu);
    CHECK(engine.policy().thermal == ThermalMode::Cpu);
    engine.step();
    CHECK(engine.report().thermal == ThermalMode::Cpu);
    CHECK(engine.policy().thermal == ThermalMode::Cpu);
}
