#include <doctest/doctest.h>

#include "coilgun/coilgun_cuda.hpp"
#include "gpu_engine_fixture.hpp"

#include <algorithm>

TEST_CASE("graph variant key includes execution policy") {
    using namespace coilgun::simulation::cuda;
    GpuGraphVariantKey first;
    GpuGraphVariantKey second;
    second.thermal = ThermalMode::Gpu;
    CHECK_FALSE(first == second);
}

TEST_CASE("graph variant changes for topology but not stage voltage") {
    using namespace coilgun::simulation::cuda;
    auto state = gpu_test::state(1, 2, 1, false);
    state.stage_voltages = {10.0, 20.0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    config.solver = SolverMode::Batched;
    GpuEngine engine(gpu_test::geometry(2, 1, false), std::move(state), config);

    const auto first = engine.graph_variant();
    engine.set_stage_voltage(0, 12.0);
    engine.select_graph_variant_at_boundary();
    CHECK(engine.graph_variant() == first);

    engine.set_mutual_stage_mask({1, 0});
    CHECK_FALSE(engine.graph_variant() == first);
}

TEST_CASE("fixed-shape graph replays complete resident step") {
    using namespace coilgun::simulation::cuda;
    if (!cuda_device_available()) return;
    auto state = gpu_test::state(1, 1, 1, false);
    state.stage_voltages = {25.0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    config.solver = SolverMode::Batched;
    GpuEngine engine(gpu_test::geometry(1, 1, false), std::move(state), config);

    engine.step();
    const auto first_rebuilds = engine.report().graph_rebuild_count;
    engine.set_stage_voltage(0, 30.0);
    engine.step();

    CHECK(engine.report().gpu_executed);
    CHECK(engine.report().backend == BackendMode::Graph);
    CHECK(first_rebuilds == 1);
    CHECK(engine.report().graph_rebuild_count == first_rebuilds);
    CHECK(engine.result().completed_steps == 2);
    const auto& order = engine.pipeline_order();
    CHECK(std::find(order.begin(), order.end(), PipelineStage::Matrix) != order.end());
    CHECK(std::find(order.begin(), order.end(), PipelineStage::Solver) != order.end());
    CHECK(std::find(order.begin(), order.end(), PipelineStage::State) != order.end());
}

TEST_CASE("graph replay updates runtime masks without recapture") {
    using namespace coilgun::simulation::cuda;
    if (!cuda_device_available()) return;
    auto state = gpu_test::state(1, 1, 1, false);
    state.stage_voltages = {25.0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    config.solver = SolverMode::Batched;
    GpuEngine engine(gpu_test::geometry(1, 1, false), std::move(state), config);

    engine.step();
    REQUIRE(engine.report().graph_rebuild_count == 1);
    engine.set_mutual_stage_mask({0});
    engine.step();

    CHECK(engine.report().graph_rebuild_count == 1);
    REQUIRE(engine.state().m1.size() == 1);
    REQUIRE(engine.state().dm1.size() == 1);
    CHECK(engine.state().m1[0] == doctest::Approx(0.0));
    CHECK(engine.state().dm1[0] == doctest::Approx(0.0));
}

TEST_CASE("graph thermal pipeline matches direct resident thermal state") {
    using namespace coilgun::simulation::cuda;
    if (!cuda_device_available()) return;
    auto graph_state = gpu_test::state(1, 1, 1, true);
    graph_state.currents = {1.0, 4.0};
    graph_state.stage_voltages = {20.0};
    auto direct_state = graph_state;

    GpuExecutionConfig graph_config;
    graph_config.backend = BackendMode::Graph;
    graph_config.solver = SolverMode::Batched;
    graph_config.thermal = ThermalMode::Gpu;
    GpuExecutionConfig direct_config = graph_config;
    direct_config.backend = BackendMode::Direct;
    GpuEngine graph(gpu_test::geometry(1, 1, true), std::move(graph_state), graph_config);
    GpuEngine direct(gpu_test::geometry(1, 1, true), std::move(direct_state), direct_config);

    graph.step();
    direct.step();

    CHECK(graph.state().currents == direct.state().currents);
    CHECK(graph.state().position == direct.state().position);
    CHECK(graph.state().velocity == direct.state().velocity);
    REQUIRE(graph.state().temperatures.size() == direct.state().temperatures.size());
    CHECK(graph.state().temperatures[0] ==
          doctest::Approx(direct.state().temperatures[0]).epsilon(1e-12));
    CHECK(graph.state().resistances[0] ==
          doctest::Approx(direct.state().resistances[0]).epsilon(1e-12));
    CHECK(std::find(graph.pipeline_order().begin(), graph.pipeline_order().end(),
                    PipelineStage::Thermal) != graph.pipeline_order().end());
}

TEST_CASE("graph capture failure rolls back and locks CPU fallback") {
    using namespace coilgun::simulation::cuda;
    if (!cuda_device_available()) return;
    auto graph_state = gpu_test::state(1, 1, 1, false);
    graph_state.stage_voltages = {20.0};
    auto fallback_state = graph_state;
    GpuExecutionConfig graph_config;
    graph_config.backend = BackendMode::Graph;
    graph_config.solver = SolverMode::Batched;
    GpuExecutionConfig fallback_config;
    fallback_config.backend = BackendMode::Fallback;
    fallback_config.solver = SolverMode::Eigen;
    detail::GpuEngineFaultInjection fault;
    fault.fail_graph_capture = true;
    GpuEngine graph(gpu_test::geometry(), std::move(graph_state), graph_config, {}, fault);
    GpuEngine fallback(gpu_test::geometry(), std::move(fallback_state), fallback_config);

    graph.step();
    fallback.step();

    CHECK(graph.state().currents == fallback.state().currents);
    CHECK(graph.state().position == fallback.state().position);
    CHECK(graph.state().velocity == fallback.state().velocity);
    CHECK(graph.report().backend == BackendMode::Fallback);
    CHECK(graph.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK_FALSE(graph.report().gpu_executed);
}
