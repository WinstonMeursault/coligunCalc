#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_engine.hpp"

#include <cmath>

using namespace coilgun::simulation::cuda;

namespace {

GpuGeometryInput physics_geometry() {
    GpuGeometryInput geometry;
    geometry.n_stages = 2;
    geometry.n_filaments = 1;
    geometry.stage_geometry = {0.05, 0.05};
    geometry.stage_inner_radii = {0.045, 0.045};
    geometry.stage_outer_radii = {0.055, 0.055};
    geometry.stage_lengths = {0.01, 0.01};
    geometry.stage_turns = {20, 20};
    geometry.stage_positions = {0.0, 0.0};
    geometry.filament_geometry = {0.03};
    geometry.filament_inner_radii = {0.029};
    geometry.filament_outer_radii = {0.031};
    geometry.filament_lengths = {0.01};
    geometry.filament_positions = {0.04};
    geometry.stage_resistances = {1.0, 1.0};
    geometry.stage_inductances = {2.0, 3.0};
    geometry.filament_resistances = {0.5};
    geometry.filament_inductances = {1.5};
    geometry.stage_mutual_inductances = {0.0, 0.25, 0.25, 0.0};
    geometry.filament_mutual_inductances = {0.0};
    return geometry;
}

GpuEngineState physics_state(double position, double velocity) {
    GpuEngineState state;
    state.currents = {0.0, 0.0, 0.0};
    state.m1 = {0.0, 0.0};
    state.dm1 = {0.0, 0.0};
    state.stage_voltages = {10.0, 10.0};
    state.active_mask = {1};
    state.trigger_mask = {1, 1};
    state.velocity = {velocity};
    state.position = {position};
    state.dt = 1.0e-5;
    state.mass = 1.0;
    return state;
}

} // namespace

TEST_CASE("engine assembles coupled circuit derivatives and advances position-dependent mutuals") {
    auto geometry = physics_geometry();
    auto state = physics_state(0.0, 100.0);
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;

    GpuEngine engine(std::move(geometry), std::move(state), config);
    engine.step();

    const auto first_currents = engine.state().currents;
    const double first_position = engine.state().position[0];
    const double first_mutual = engine.state().m1[0];
    REQUIRE(first_currents.size() == 3);
    CHECK(std::abs(first_currents[0]) > 0.0);
    CHECK(std::abs(first_currents[1]) > 0.0);
    CHECK(std::abs(first_currents[2]) > 0.0);
    CHECK(first_position == doctest::Approx(100.0 * engine.state().dt));

    engine.step();
    CHECK(engine.state().position[0] > first_position);
    CHECK(std::abs(engine.state().m1[0] - first_mutual) > 1.0e-9);
}

TEST_CASE("CPU fallback stores fused mutual pair values for every active pair") {
    auto geometry = physics_geometry();
    geometry.n_filaments = 2;
    geometry.filament_geometry = {0.03, 0.032};
    geometry.filament_inner_radii = {0.029, 0.031};
    geometry.filament_outer_radii = {0.031, 0.033};
    geometry.filament_lengths = {0.01, 0.012};
    geometry.filament_positions = {0.04, 0.06};
    geometry.filament_resistances = {0.5, 0.6};
    geometry.filament_inductances = {1.5, 1.6};
    geometry.filament_mutual_inductances = {0.0, 0.1, 0.1, 0.0};

    GpuEngineState state;
    state.currents.assign(4, 0.0);
    state.m1.assign(4, 0.0);
    state.dm1.assign(4, 0.0);
    state.stage_voltages = {10.0, 10.0};
    state.active_mask = {1};
    state.trigger_mask = {1, 1};
    state.velocity = {0.0};
    state.position = {0.0};
    state.dt = 1.0e-5;
    state.mass = 1.0;

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;

    GpuEngine engine(geometry, state, config);
    engine.step();

    for (std::size_t stage = 0; stage < geometry.n_stages; ++stage) {
        for (std::size_t filament = 0; filament < geometry.n_filaments; ++filament) {
            const double separation = geometry.stage_positions[stage] -
                (geometry.filament_positions[filament] - engine.state().position[0]);
            const auto expected = coilgun::physics::mutual_detail::mutual_inductance_coil_pair(
                geometry.stage_inner_radii[stage], geometry.stage_outer_radii[stage],
                geometry.stage_lengths[stage], geometry.stage_turns[stage],
                geometry.filament_inner_radii[filament], geometry.filament_outer_radii[filament],
                geometry.filament_lengths[filament], 1, separation, 9, true);
            const auto pair = stage * geometry.n_filaments + filament;
            CHECK(engine.state().m1[pair] == doctest::Approx(expected.mutual));
            CHECK(engine.state().dm1[pair] == doctest::Approx(expected.gradient));
        }
    }
}

TEST_CASE("engine stage mutual coupling changes the solved derivative") {
    auto coupled_geometry = physics_geometry();
    auto uncoupled_geometry = coupled_geometry;
    uncoupled_geometry.stage_mutual_inductances = {0.0, 0.0, 0.0, 0.0};
    auto coupled_state = physics_state(0.0, 0.0);
    auto uncoupled_state = coupled_state;
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;

    GpuEngine coupled(std::move(coupled_geometry), std::move(coupled_state), config);
    GpuEngine uncoupled(std::move(uncoupled_geometry), std::move(uncoupled_state), config);
    coupled.step();
    uncoupled.step();

    CHECK(std::abs(coupled.state().currents[0] - uncoupled.state().currents[0]) > 1.0e-8);
}

TEST_CASE("engine keeps circuit stages active when their mutual mask is disabled") {
    auto coupled_geometry = physics_geometry();
    auto masked_geometry = coupled_geometry;
    auto coupled_state = physics_state(0.0, 0.0);
    auto masked_state = coupled_state;
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;

    GpuEngine coupled(std::move(coupled_geometry), std::move(coupled_state), config);
    GpuEngine masked(std::move(masked_geometry), std::move(masked_state), config);
    masked.set_mutual_stage_mask({1, 0});
    coupled.step();
    masked.step();

    CHECK(std::abs(masked.state().currents[1]) > 0.0);
    CHECK(std::abs(masked.state().currents[0]) > 0.0);
    CHECK(masked.state().m1[1] == doctest::Approx(0.0));
    CHECK(masked.state().dm1[1] == doctest::Approx(0.0));
}

TEST_CASE("CPU engine applies thermal update after Euler state update") {
    auto geometry = physics_geometry();
    geometry.thermal_enabled = true;
    auto state = physics_state(0.0, 0.0);
    state.currents[2] = 10.0;
    state.temperatures = {293.0};
    state.filament_masses = {0.01};
    state.reference_resistances = {0.5};
    state.filament_materials = {1};
    state.material_density = 1.0;
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;
    config.thermal = ThermalMode::Cpu;

    GpuEngine engine(std::move(geometry), std::move(state), config);
    engine.step();

    const auto& order = engine.pipeline_order();
    REQUIRE(order.size() == 6);
    CHECK(order[4] == PipelineStage::State);
    CHECK(order[5] == PipelineStage::Thermal);
    CHECK(engine.state().temperatures[0] > 293.0);
    CHECK(engine.state().resistances[0] > 0.5);
}
