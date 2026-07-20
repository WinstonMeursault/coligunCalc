/**
 * @file test_gpu_engine_layout.cpp
 * @brief Host-side tests for the fixed GPU state layout.
 */

#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_state_layout.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "gpu_engine_fixture.hpp"

#include <stdexcept>
#include <limits>

using coilgun::simulation::cuda::GpuStateLayout;

TEST_CASE("GPU state layout maps a single-stage state in row-major order") {
    const GpuStateLayout layout(1, 1, 3);

    CHECK(layout.currents(0, 0) == 0);
    CHECK(layout.currents(0, 1) == 1);
    CHECK(layout.m1(0, 0, 0) == 0);
    CHECK(layout.m1(0, 0, 2) == 2);
    CHECK(layout.dm1(0, 0, 2) == 2);
    CHECK(layout.temperatures(0, 2) == 2);
    CHECK(layout.active_mask(0) == 0);
    CHECK(layout.trigger_mask(0, 0) == 0);
}

TEST_CASE("GPU state layout preserves stage and batch positions") {
    const GpuStateLayout layout(2, 3, 2);

    CHECK(layout.currents(1, 3) == 8);
    CHECK(layout.m1(0, 2, 1) == 5);
    CHECK(layout.m1(1, 0, 0) == 6);
    CHECK(layout.dm1(1, 2, 1) == 11);
    CHECK(layout.temperatures(1, 1) == 3);
    CHECK(layout.active_mask(1) == 1);
    CHECK(layout.trigger_mask(1, 2) == 5);
}

TEST_CASE("GPU state layout exposes dimensions and inactive masks do not compact indices") {
    const GpuStateLayout layout(3, 2, 4);

    CHECK(layout.batch_size() == 3);
    CHECK(layout.stage_count() == 2);
    CHECK(layout.filament_count() == 4);
    CHECK(layout.current_dimension() == 6);
    CHECK(layout.m1(2, 1, 3) == 23);
    CHECK(layout.m1(1, 1, 3) == 15);
    CHECK(layout.system_matrix(2, 5, 5) == 107);
    CHECK(layout.rhs(2, 5) == 17);
    CHECK(layout.system_matrix_size() == 3 * 6 * 6);
    CHECK(layout.rhs_size() == 3 * 6);
}

TEST_CASE("GPU state layout checks host indices and dimensions") {
    CHECK_THROWS_AS(GpuStateLayout(0, 1, 1), std::invalid_argument);
    CHECK_THROWS_AS(GpuStateLayout(1, 0, 1), std::invalid_argument);
    CHECK_THROWS_AS(GpuStateLayout(1, 1, 0), std::invalid_argument);

    const GpuStateLayout layout(2, 2, 2);
    CHECK_THROWS_AS(layout.currents(2, 0), std::out_of_range);
    CHECK_THROWS_AS(layout.m1(0, 2, 0), std::out_of_range);
    CHECK_THROWS_AS(layout.dm1(0, 0, 2), std::out_of_range);
    CHECK_THROWS_AS(layout.temperatures(0, 2), std::out_of_range);
    CHECK_THROWS_AS(layout.active_mask(2), std::out_of_range);
    CHECK_THROWS_AS(layout.trigger_mask(2, 0), std::out_of_range);
    CHECK_THROWS_AS(layout.trigger_mask(0, 2), std::out_of_range);
}

TEST_CASE("GPU engine supports normalized B=1 state and boundary reset") {
    using namespace coilgun::simulation::cuda;

    auto geometry = gpu_test::geometry(1, 2);
    auto state = gpu_test::state(1, 1, 2);

    GpuEngine engine(std::move(geometry), std::move(state));
    CHECK(engine.layout().B == 1);
    CHECK(engine.layout().D == 3);
    CHECK(engine.graph_variant().batch_size == 1);
    CHECK(engine.report().graph_rebuild_count == 0);

    engine.step();
    CHECK(engine.result().completed_steps == 1);
    engine.run(GpuRunBoundary{2, false});
    CHECK(engine.result().completed_steps == 3);

    engine.reset();
    CHECK(engine.result().completed_steps == 0);
    CHECK(engine.report().graph_rebuild_count == 0);
    engine.shutdown();
    CHECK_THROWS_AS(engine.step(), std::logic_error);
}

TEST_CASE("GPU engine preserves B>1 physical positions and stage mask boundary") {
    using namespace coilgun::simulation::cuda;

    auto geometry = gpu_test::geometry(2, 2);
    auto state = gpu_test::state(2, 2, 2);
    state.active_mask = {1, 0};

    GpuEngine engine(std::move(geometry), std::move(state));
    engine.set_stage_mask({1, 0, 0, 1});
    CHECK(engine.graph_variant().stage_mask == std::vector<std::uint8_t>{1, 0, 0, 1});
    CHECK(engine.report().graph_rebuild_count == 0);
    engine.set_stage_mask({1, 0, 0, 1});
    CHECK(engine.report().graph_rebuild_count == 0);
    CHECK(engine.layout().m1(1, 1, 1) == 7);
    CHECK(engine.layout().active_mask(1) == 1);
}

TEST_CASE("GPU engine never resolves persistent execution for a multi-batch layout") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(2, 2);
    auto state = gpu_test::state(2, 2, 2);
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;
    config.deterministic = false;

    GpuEngine engine(std::move(geometry), std::move(state), config);
    CHECK(engine.policy().backend != BackendMode::Persistent);
}

TEST_CASE("GPU engine keeps a stage mask for every simulation") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(2, 2);
    auto state = gpu_test::state(2, 2, 2);
    GpuEngine engine(std::move(geometry), std::move(state));
    CHECK_THROWS_AS(engine.set_stage_mask({1, 0}), std::invalid_argument);
    engine.set_stage_mask({1, 0, 0, 1});
    CHECK(engine.graph_variant().stage_mask == std::vector<std::uint8_t>{1, 0, 0, 1});
}

TEST_CASE("GPU engine validates geometry values and layout overflow") {
    using namespace coilgun::simulation::cuda;
    GpuEngineState state;
    state.active_mask = {1};
    state.trigger_mask = {1};
    GpuGeometryInput missing{1, 1, {}, {1.0}};
    CHECK_THROWS_AS(GpuEngine(std::move(missing), state), std::invalid_argument);
    GpuGeometryInput nonfinite{1, 1, {1.0}, {std::numeric_limits<double>::infinity()}};
    CHECK_THROWS_AS(GpuEngine(std::move(nonfinite), state), std::invalid_argument);
    CHECK_THROWS_AS(GpuStateLayout(std::numeric_limits<std::size_t>::max(), 2, 2),
                    std::overflow_error);
}

TEST_CASE("GPU engine rejects placeholder dynamics and incomplete thermal material contracts") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry();
    geometry.stage_inner_radii = {0.045};
    geometry.stage_outer_radii = {0.055};
    geometry.stage_lengths = {0.03};
    geometry.stage_turns = {120};
    geometry.stage_positions = {0.0};
    geometry.filament_inner_radii = {0.008};
    geometry.filament_outer_radii = {0.012};
    geometry.filament_lengths = {0.02};
    geometry.filament_positions = {0.04};

    GpuEngineState state;
    state.currents = {1.0, 2.0};
    state.m1 = {0.01};
    state.dm1 = {0.02};
    state.active_mask = {1};
    state.trigger_mask = {1};
    state.dt = 0.0;
    state.mass = 0.0;
    CHECK_THROWS_AS(GpuEngine(geometry, state), std::invalid_argument);

    state.dt = 2.5e-5;
    state.mass = 0.02;
    state.filament_masses = {0.02};
    state.reference_resistances = {0.003};
    state.filament_materials = {0};
    geometry.thermal_enabled = true;
    state.temperatures = {293.0};
    CHECK_THROWS_AS(GpuEngine(geometry, state), std::invalid_argument);
}

TEST_CASE("GPU engine advances active state and preserves inactive state") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(1, 1);

    GpuEngineState state;
    state.currents = {1.0, 2.0, 9.0, 8.0};
    state.m1 = {0.01, 0.02};
    state.dm1 = {0.0, 0.0};
    state.active_mask = {1, 0};
    state.trigger_mask = {1, 1};
    state.dt = 2.5e-5;
    state.mass = 0.02;
    state.velocity = {0.0, 3.0};
    state.position = {0.1, 0.2};
    state.stage_voltages = {1.0, 0.0};

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Eigen;
    GpuEngine engine(std::move(geometry), std::move(state), config);
    engine.step();

    CHECK(engine.state().currents[0] != doctest::Approx(1.0));
    CHECK(engine.state().currents[2] == doctest::Approx(9.0));
    CHECK(engine.state().currents[3] == doctest::Approx(8.0));
    CHECK(engine.state().position[1] == doctest::Approx(0.2));
    CHECK(engine.state().velocity[1] == doctest::Approx(3.0));
}

TEST_CASE("GPU engine preserves a caller-provided thermal state across first step and reset") {
    using namespace coilgun::simulation::cuda;

    auto geometry = gpu_test::geometry(1, 1, true);
    auto state = gpu_test::state(1, 1, 1, true);
    state.currents = {0.4, -0.2};
    state.m1 = {0.012};
    state.dm1 = {-0.003};
    state.temperatures = {311.0};
    state.active_mask = {1};
    state.trigger_mask = {1};
    state.velocity = {2.5};
    state.position = {0.006};
    state.filament_masses = {0.014};
    state.reference_resistances = {0.004};
    state.filament_materials = {0};
    state.resistivities = {2.1e-8};
    state.resistances = {0.0043};
    state.joule_energy = {0.0007};
    state.current_derivatives = {1.2, -0.8};
    state.stage_voltages = {12.0};
    state.reference_temperature = 293.0;
    state.material_density = 2700.0;

    const auto original_state = state;
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.thermal = ThermalMode::Cpu;
    GpuEngine engine(std::move(geometry), std::move(state), config);

    CHECK(engine.state().currents == original_state.currents);
    CHECK(engine.state().temperatures == original_state.temperatures);
    CHECK(engine.state().resistances == original_state.resistances);
    CHECK(engine.state().current_derivatives == original_state.current_derivatives);

    engine.step();
    const auto first_step_state = engine.state();
    CHECK(first_step_state.currents.size() == 2);
    CHECK(first_step_state.temperatures.size() == 1);
    CHECK(first_step_state.resistances.size() == 1);
    CHECK(first_step_state.joule_energy.size() == 1);
    CHECK(first_step_state.current_derivatives.size() == 2);

    engine.reset();
    const auto reset_state = engine.state();
    CHECK(reset_state.currents == original_state.currents);
    CHECK(reset_state.m1 == original_state.m1);
    CHECK(reset_state.dm1 == original_state.dm1);
    CHECK(reset_state.temperatures == original_state.temperatures);
    CHECK(reset_state.active_mask == original_state.active_mask);
    CHECK(reset_state.trigger_mask == original_state.trigger_mask);
    CHECK(reset_state.velocity == original_state.velocity);
    CHECK(reset_state.position == original_state.position);
    CHECK(reset_state.filament_masses == original_state.filament_masses);
    CHECK(reset_state.reference_resistances == original_state.reference_resistances);
    CHECK(reset_state.filament_materials == original_state.filament_materials);
    CHECK(reset_state.resistivities == original_state.resistivities);
    CHECK(reset_state.resistances == original_state.resistances);
    CHECK(reset_state.joule_energy == original_state.joule_energy);
    CHECK(reset_state.current_derivatives == original_state.current_derivatives);
    CHECK(reset_state.stage_voltages == original_state.stage_voltages);
    CHECK(reset_state.dt == original_state.dt);
    CHECK(reset_state.mass == original_state.mass);
    CHECK(reset_state.reference_temperature == original_state.reference_temperature);
    CHECK(reset_state.material_density == original_state.material_density);

    engine.step();
    const auto repeated_first_step_state = engine.state();
    CHECK(repeated_first_step_state.currents == first_step_state.currents);
    CHECK(repeated_first_step_state.m1 == first_step_state.m1);
    CHECK(repeated_first_step_state.dm1 == first_step_state.dm1);
    CHECK(repeated_first_step_state.temperatures == first_step_state.temperatures);
    CHECK(repeated_first_step_state.resistivities == first_step_state.resistivities);
    CHECK(repeated_first_step_state.resistances == first_step_state.resistances);
    CHECK(repeated_first_step_state.joule_energy == first_step_state.joule_energy);
    CHECK(repeated_first_step_state.current_derivatives == first_step_state.current_derivatives);
    CHECK(repeated_first_step_state.velocity == first_step_state.velocity);
    CHECK(repeated_first_step_state.position == first_step_state.position);
}

TEST_CASE("GPU engine gives bounded semantics for empty and unbounded runs") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry();
    auto state = gpu_test::state(1, 1, 1);
    state.active_mask = {0};
    GpuEngine engine(std::move(geometry), std::move(state));
    engine.run(GpuRunBoundary{0, true});
    CHECK(engine.result().finished);

    auto active_geometry = gpu_test::geometry();
    auto active_state = gpu_test::state(1, 1, 1);
    GpuEngine active(std::move(active_geometry), std::move(active_state));
    CHECK_THROWS_AS(active.run(GpuRunBoundary{}), std::logic_error);
}

TEST_CASE("GPU engine requires temperatures only when thermal mode is enabled") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(1, 1, true);
    auto state = gpu_test::state(1, 1, 1);
    state.active_mask = {1};
    CHECK_THROWS_AS(GpuEngine(std::move(geometry), std::move(state)), std::invalid_argument);
}

TEST_CASE("GPU engine rejects state buffers that do not match normalized layout") {
    using namespace coilgun::simulation::cuda;
    GpuGeometryInput geometry{1, 1, {1.0}, {1.0}};
    GpuEngineState state;
    state.active_mask = {1};
    state.trigger_mask = {1};
    CHECK_THROWS_AS(GpuEngine(std::move(geometry), std::move(state)),
                    std::invalid_argument);
}

TEST_CASE("GPU engine applies independent row-major stage voltages to each batch") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(2, 1);
    geometry.stage_resistances = {0.0, 0.0};
    geometry.stage_inductances = {2.0, 2.0};
    geometry.filament_resistances = {1.0};
    geometry.filament_inductances = {1.0};
    geometry.stage_mutual_inductances = {0.0, 0.0, 0.0, 0.0};
    geometry.filament_mutual_inductances = {0.0};

    auto state = gpu_test::state(2, 2, 1);
    state.stage_voltages = {10.0, 20.0, 30.0, 40.0};
    state.mutual_stage_mask = {0, 0, 0, 0};

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    GpuEngine engine(std::move(geometry), std::move(state), config);
    engine.step();

    CHECK(engine.state().currents[0] == doctest::Approx(10.0 * engine.state().dt / 2.0));
    CHECK(engine.state().currents[1] == doctest::Approx(20.0 * engine.state().dt / 2.0));
    CHECK(engine.state().currents[3] == doctest::Approx(30.0 * engine.state().dt / 2.0));
    CHECK(engine.state().currents[4] == doctest::Approx(40.0 * engine.state().dt / 2.0));
}

TEST_CASE("GPU engine validates step boundary state atomically") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(2, 1);
    auto state = gpu_test::state(2, 2, 1);
    state.stage_mask = {1, 0, 1, 0};
    state.mutual_stage_mask = {0, 1, 0, 1};
    state.stage_voltages = {1.0, 2.0, 3.0, 4.0};
    GpuEngine engine(std::move(geometry), std::move(state));

    const auto before = engine.state();
    const auto variant_before = engine.graph_variant();
    CHECK_THROWS_AS(engine.set_step_boundary_state(
                        {0, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1},
                        {1.0, std::numeric_limits<double>::quiet_NaN(), 3.0, 4.0}),
                    std::invalid_argument);

    CHECK(engine.state().active_mask == before.active_mask);
    CHECK(engine.state().trigger_mask == before.trigger_mask);
    CHECK(engine.state().stage_mask == before.stage_mask);
    CHECK(engine.state().mutual_stage_mask == before.mutual_stage_mask);
    CHECK(engine.state().stage_voltages == before.stage_voltages);
    CHECK(engine.graph_variant() == variant_before);

    CHECK_THROWS_AS(engine.set_step_boundary_state(
                        {1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1},
                        {1.0, 2.0, 3.0, 4.0}),
                    std::invalid_argument);
    CHECK(engine.state().active_mask == before.active_mask);
    CHECK(engine.state().stage_voltages == before.stage_voltages);
}

TEST_CASE("GPU engine preserves every inactive batch state including thermal values") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(1, 1, true);
    auto state = gpu_test::state(2, 1, 1, true);
    state.active_mask = {1, 0};
    state.stage_voltages = {10.0, 20.0};
    state.currents = {0.0, 0.0, 7.0, 8.0};
    state.m1 = {0.1, 0.2};
    state.dm1 = {0.3, 0.4};
    state.temperatures = {301.0, 377.0};
    state.resistivities = {1.0, 2.0};
    state.resistances = {3.0, 4.0};
    state.joule_energy = {5.0, 6.0};
    state.current_derivatives = {9.0, 10.0, 11.0, 12.0};
    state.velocity = {0.0, 13.0};
    state.position = {0.0, 14.0};

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.thermal = ThermalMode::Cpu;
    GpuEngine engine(std::move(geometry), std::move(state), config);
    const auto before = engine.state();
    engine.step();

    CHECK(engine.state().currents[2] == before.currents[2]);
    CHECK(engine.state().currents[3] == before.currents[3]);
    CHECK(engine.state().m1[1] == before.m1[1]);
    CHECK(engine.state().dm1[1] == before.dm1[1]);
    CHECK(engine.state().temperatures[1] == before.temperatures[1]);
    CHECK(engine.state().resistivities[1] == before.resistivities[1]);
    CHECK(engine.state().resistances[1] == before.resistances[1]);
    CHECK(engine.state().joule_energy[1] == before.joule_energy[1]);
    CHECK(engine.state().current_derivatives[2] == before.current_derivatives[2]);
    CHECK(engine.state().current_derivatives[3] == before.current_derivatives[3]);
    CHECK(engine.state().velocity[1] == before.velocity[1]);
    CHECK(engine.state().position[1] == before.position[1]);
}

TEST_CASE("GPU engine keeps fixed row indices when one row becomes inactive") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(1, 1);
    auto state = gpu_test::state(2, 1, 1);
    state.currents = {0.2, 0.4, 0.6, 0.8};
    state.stage_voltages = {10.0, 20.0};

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    GpuEngine batch(std::move(geometry), std::move(state), config);
    batch.step();
    const auto completed_row = batch.state();

    batch.set_step_boundary_state({0, 1}, {1, 1}, {1, 1}, {1, 1}, {0.0, 20.0});
    batch.step();

    CHECK(batch.state().currents[0] == completed_row.currents[0]);
    CHECK(batch.state().currents[1] == completed_row.currents[1]);
    CHECK(batch.state().m1[0] == completed_row.m1[0]);
    CHECK(batch.state().dm1[0] == completed_row.dm1[0]);
    CHECK(batch.state().velocity[0] == completed_row.velocity[0]);
    CHECK(batch.state().position[0] == completed_row.position[0]);

    auto standalone_state = gpu_test::state(1, 1, 1);
    standalone_state.currents = {0.6, 0.8};
    standalone_state.stage_voltages = {20.0};
    GpuEngine standalone(gpu_test::geometry(1, 1), std::move(standalone_state), config);
    standalone.step();
    standalone.set_step_boundary_state({1}, {1}, {1}, {1}, {20.0});
    standalone.step();

    CHECK(batch.state().currents[2] == doctest::Approx(standalone.state().currents[0]));
    CHECK(batch.state().currents[3] == doctest::Approx(standalone.state().currents[1]));
    CHECK(batch.state().m1[1] == doctest::Approx(standalone.state().m1[0]));
    CHECK(batch.state().dm1[1] == doctest::Approx(standalone.state().dm1[0]));
    CHECK(batch.state().velocity[1] == doctest::Approx(standalone.state().velocity[0]));
    CHECK(batch.state().position[1] == doctest::Approx(standalone.state().position[0]));
}

TEST_CASE("GPU engine keeps circuit and mutual masks independent at a boundary") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry(2, 1);
    auto state = gpu_test::state(1, 2, 1);
    state.currents = {1.0, 2.0, 3.0};
    state.stage_voltages = {10.0, 20.0};

    auto masked_state = state;
    masked_state.stage_mask = {1, 0};
    masked_state.mutual_stage_mask = {1, 1};
    auto no_mutual_state = state;
    no_mutual_state.stage_mask = {1, 0};
    no_mutual_state.mutual_stage_mask = {0, 0};

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    GpuEngine masked(std::move(geometry), std::move(masked_state), config);
    GpuEngine no_mutual(gpu_test::geometry(2, 1), std::move(no_mutual_state), config);
    masked.step();
    no_mutual.step();

    CHECK(masked.state().currents[1] == doctest::Approx(2.0));
    CHECK(no_mutual.state().currents[1] == doctest::Approx(2.0));
    CHECK(std::abs(masked.state().dm1[0]) > 0.0);
    CHECK(no_mutual.state().dm1[0] == doctest::Approx(0.0));
    CHECK(masked.graph_variant().stage_mask == std::vector<std::uint8_t>{1, 0});
    CHECK(masked.graph_variant().mutual_stage_mask == std::vector<std::uint8_t>{1, 1});
}
