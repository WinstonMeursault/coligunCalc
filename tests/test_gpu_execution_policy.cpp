#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_execution_report.hpp"
#include "coilgun/simulation/cuda/gpu_engine.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "gpu_engine_fixture.hpp"

#include <sstream>
#include <limits>

using namespace coilgun::simulation::cuda;

TEST_CASE("GPU execution report has stable defaults") {
    const ExecutionReport report{};

    CHECK(report.backend == Backend::Auto);
    CHECK(report.solver == Solver::Auto);
    CHECK(report.precision == PrecisionMode::Full);
    CHECK(report.thermal == ThermalMode::Auto);
    CHECK_FALSE(report.calibrated);
    CHECK_FALSE(report.gpu_executed);
    CHECK_FALSE(report.precision_fallback);
    CHECK(report.graph_rebuild_count == 0);
    CHECK(report.fallback_count == 0);
    CHECK(report.gpu_time_ms == doctest::Approx(0.0));
    CHECK(report.solver_time_ms == doctest::Approx(0.0));
    CHECK(report.thermal_time_ms == doctest::Approx(0.0));
    CHECK(report.transfer_time_ms == doctest::Approx(0.0));
    CHECK(report.max_condition_estimate == doctest::Approx(0.0));
    CHECK(report.fallback_reason.empty());
    CHECK(report.static_fallback_reason == FallbackReason::None);
    CHECK(report.runtime_fallback_reason == FallbackReason::None);
}

TEST_CASE("GPU execution report provides stable enum formatting") {
    CHECK(to_string(Backend::Graph) == "graph");
    CHECK(to_string(Solver::CuSolver) == "cusolver");
    CHECK(to_string(PrecisionMode::Aggressive) == "aggressive");
    CHECK(to_string(ThermalMode::Gpu) == "gpu");
    CHECK(to_string(Backend::Direct) == "direct");

    std::ostringstream stream;
    stream << Backend::Persistent << '/' << Solver::Eigen << '/'
           << PrecisionMode::Standard << '/' << ThermalMode::Cpu;
    CHECK(stream.str() == "persistent/eigen/standard/cpu");
}

TEST_CASE("GPU backend defaults to graph and validates launch settings") {
    GpuBackend backend;
    CHECK(backend.backend == BackendMode::Graph);
    CHECK(backend.use_persistent);
    const auto single_stage = single_stage_default_backend();
    CHECK(single_stage.backend == BackendMode::Graph);
    CHECK_FALSE(single_stage.use_persistent);
    CHECK_NOTHROW(backend.validate());

    backend.threads_per_block = 1;
    CHECK_NOTHROW(backend.validate());
    backend.threads_per_block = 128;
    CHECK_NOTHROW(backend.validate());
    backend.threads_per_block = 512;
    CHECK_NOTHROW(backend.validate());

    backend.threads_per_block = 3;
    CHECK_THROWS_AS(backend.validate(), std::invalid_argument);
    backend.threads_per_block = 1024;
    CHECK_THROWS_AS(backend.validate(), std::invalid_argument);
    backend.threads_per_block = 512;
    backend.device_id = -1;
    CHECK_THROWS_AS(backend.validate(), std::invalid_argument);
    backend.device_id = 0;
    backend.max_batch_sims = 0;
    CHECK_THROWS_AS(backend.validate(), std::invalid_argument);
}

TEST_CASE("GPU execution report merge aggregates measurements") {
    ExecutionReport total;
    total.backend = Backend::Graph;
    total.solver = Solver::CuSolver;
    total.precision = PrecisionMode::Full;
    total.thermal = ThermalMode::Gpu;
    total.calibrated = true;
    total.gpu_time_ms = 1.5;
    total.solver_time_ms = 0.5;
    total.max_condition_estimate = 10.0;

    ExecutionReport variant;
    variant.precision_fallback = true;
    variant.graph_rebuild_count = 2;
    variant.fallback_count = 1;
    variant.gpu_time_ms = 2.5;
    variant.solver_time_ms = 1.25;
    variant.thermal_time_ms = 0.75;
    variant.transfer_time_ms = 0.25;
    variant.max_condition_estimate = 25.0;
    variant.fallback_reason = "solver failure";

    total.merge(variant);

    CHECK(total.calibrated);
    CHECK(total.precision_fallback);
    CHECK(total.graph_rebuild_count == 2);
    CHECK(total.fallback_count == 1);
    CHECK(total.gpu_time_ms == doctest::Approx(4.0));
    CHECK(total.solver_time_ms == doctest::Approx(1.75));
    CHECK(total.thermal_time_ms == doctest::Approx(0.75));
    CHECK(total.transfer_time_ms == doctest::Approx(0.25));
    CHECK(total.max_condition_estimate == doctest::Approx(25.0));
    CHECK(total.fallback_reason == "solver failure");
}

TEST_CASE("GPU execution report merge reconciles execution metadata") {
    ExecutionReport total;
    total.device_id = 0;
    total.threads_per_block = 128;
    total.profiling_enabled = false;

    ExecutionReport other = total;
    other.device_id = 1;
    other.threads_per_block = 256;
    other.profiling_enabled = true;
    total.merge(other);

    CHECK(total.metadata_conflict);
    CHECK(total.device_id == 0);
    CHECK(total.threads_per_block == 128);
    CHECK(total.profiling_enabled);
}

TEST_CASE("GPU execution report merge detects conflicting resolved precision") {
    ExecutionReport total;
    ExecutionReport other;
    total.precision = PrecisionMode::Standard;
    other.precision = PrecisionMode::Aggressive;

    total.merge(other);

    CHECK(total.metadata_conflict);
}

TEST_CASE("GPU execution report merge detects conflicting resolved thermal modes") {
    ExecutionReport total;
    ExecutionReport other;
    total.thermal = ThermalMode::Cpu;
    other.thermal = ThermalMode::Gpu;

    total.merge(other);

    CHECK(total.metadata_conflict);
}

TEST_CASE("GPU execution report merge accepts matching resolved modes") {
    ExecutionReport total;
    ExecutionReport other;
    total.precision = PrecisionMode::Standard;
    other.precision = PrecisionMode::Standard;
    total.thermal = ThermalMode::Cpu;
    other.thermal = ThermalMode::Cpu;

    total.merge(other);

    CHECK_FALSE(total.metadata_conflict);
}

TEST_CASE("GPU execution report merge treats resolved Auto thermal as a default") {
    ExecutionReport total;
    ExecutionReport other;
    total.thermal = ThermalMode::Cpu;
    other.thermal = ThermalMode::Auto;

    total.merge(other);

    CHECK_FALSE(total.metadata_conflict);
}

TEST_CASE("GPU execution planner resolves small and batched workloads") {
    const GpuCapability capability{};
    const GpuExecutionConfig config{};

    const auto small = GpuExecutionPlanner::plan(1, 3, 1, false,
                                                 capability, config);
    CHECK(small.backend == BackendMode::Fallback);
    CHECK(small.solver == SolverMode::Eigen);
    CHECK(small.thermal == ThermalMode::Disabled);

    const auto large = GpuExecutionPlanner::plan(4, 128, 8, true,
                                                 capability, config);
    CHECK(large.backend == BackendMode::Fallback);
    CHECK(large.backend_fallback_reason == FallbackReason::MetadataConflict);
    CHECK(large.solver == SolverMode::Batched);
    CHECK(large.thermal == ThermalMode::Gpu);
}

TEST_CASE("GPU execution planner preserves standard precision and explicit choices") {
    GpuExecutionConfig config;
    config.precision = PrecisionMode::Standard;
    config.backend = BackendMode::Persistent;
    config.solver = SolverMode::Batched;
    config.thermal = ThermalMode::Gpu;
    config.deterministic = true;

    GpuCapability capability;
    capability.persistent_is_deterministic = false;
    capability.supports_gpu_thermal = false;

    const auto policy = GpuExecutionPlanner::plan(2, 4, 1, true,
                                                  capability, config);
    CHECK(policy.precision == PrecisionMode::Standard);
    CHECK(policy.backend == BackendMode::Fallback);
    CHECK(policy.solver == SolverMode::Batched);
    CHECK(policy.thermal == ThermalMode::Cpu);
}

TEST_CASE("Deterministic execution falls back from nondeterministic persistent mode") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;
    config.deterministic = true;

    const auto policy = GpuExecutionPlanner::plan(1, 4, 1, false,
                                                  GpuCapability{}, config);
    CHECK(policy.backend == BackendMode::Fallback);
}

TEST_CASE("Persistent execution resolves to fallback until engine runtime selection") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;
    config.deterministic = true;
    GpuCapability capability;
    capability.persistent_is_deterministic = true;
    const auto policy = GpuExecutionPlanner::plan(1, 1, 1, false, capability, config);
    CHECK(policy.backend == BackendMode::Fallback);
    CHECK(policy.backend_fallback_reason == FallbackReason::CapabilityUnavailable);
}

TEST_CASE("Persistent planner requires a dedicated control stream capability") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;
    config.deterministic = false;

    GpuCapability capability;
    capability.supports_persistent = true;
    capability.supports_persistent_control_stream = false;

    const auto policy = GpuExecutionPlanner::plan(1, 1, 1, false, capability, config);
    CHECK(policy.backend == BackendMode::Fallback);
    CHECK(policy.backend_fallback_reason == FallbackReason::CapabilityUnavailable);
}

TEST_CASE("Persistent planner reports determinism rejection after capability checks") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;
    config.deterministic = true;

    GpuCapability capability;
    capability.supports_persistent = true;
    capability.supports_persistent_control_stream = true;
    capability.persistent_is_deterministic = false;

    const auto policy = GpuExecutionPlanner::plan(1, 1, 1, false, capability, config);
    CHECK(policy.backend == BackendMode::Fallback);
    CHECK(policy.backend_fallback_reason == FallbackReason::DeterminismRequired);
}

TEST_CASE("Persistent capability rejection reaches the engine as static fallback") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Persistent;

    GpuEngine engine(gpu_test::geometry(), gpu_test::state(1, 1, 1), config);

    CHECK(engine.policy().backend == BackendMode::Fallback);
    CHECK(engine.policy().backend_fallback_reason == FallbackReason::CapabilityUnavailable);
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().static_fallback_reason == FallbackReason::CapabilityUnavailable);
    CHECK_FALSE(engine.context_available());

    engine.step();
    CHECK_FALSE(engine.report().gpu_executed);
    CHECK(engine.report().runtime_fallback_reason == FallbackReason::None);
}

TEST_CASE("Planner reports requested and resolved modes with stable fallback reasons") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    config.solver = SolverMode::Batched;
    config.thermal = ThermalMode::Gpu;
    GpuCapability capability;
    capability.supports_graph = false;
    capability.supports_batched_solver = false;
    capability.supports_gpu_thermal = false;
    const auto policy = GpuExecutionPlanner::plan(1, 1, 1, true, capability, config);
    CHECK(policy.requested_backend == BackendMode::Graph);
    CHECK(policy.backend == BackendMode::Fallback);
    CHECK(policy.backend_fallback_reason == FallbackReason::CapabilityUnavailable);
    CHECK(policy.requested_solver == SolverMode::Batched);
    CHECK(policy.solver_fallback_reason == FallbackReason::CapabilityUnavailable);
    CHECK(policy.requested_thermal == ThermalMode::Gpu);
    CHECK(policy.thermal_fallback_reason == FallbackReason::CapabilityUnavailable);
}

TEST_CASE("Graph request remains a static fallback until engine runtime selection") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    const auto policy = GpuExecutionPlanner::plan(8, 128, 8, false, {}, config);
    CHECK(policy.backend == BackendMode::Fallback);
    CHECK(policy.backend_fallback_reason == FallbackReason::MetadataConflict);
}

TEST_CASE("Report merge records metadata conflicts and fallback reasons") {
    ExecutionReport total;
    total.backend = Backend::Graph;
    total.solver = Solver::Eigen;
    ExecutionReport other;
    other.backend = Backend::Fallback;
    other.solver = Solver::CuSolver;
    other.runtime_fallback_reason = FallbackReason::RuntimeFailure;
    other.fallback_reason = "runtime failure";
    total.merge(other);
    CHECK(total.metadata_conflict);
    CHECK(total.runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK(total.fallback_reason == "runtime failure");
}

TEST_CASE("Report merge detects conflicting fallback reason enums") {
    ExecutionReport total;
    total.static_fallback_reason = FallbackReason::CapabilityUnavailable;
    ExecutionReport other;
    other.static_fallback_reason = FallbackReason::MetadataConflict;

    total.merge(other);

    CHECK(total.metadata_conflict);
    CHECK(total.static_fallback_reason == FallbackReason::CapabilityUnavailable);
}

TEST_CASE("Thermal policy is disabled independently of missing temperatures") {
    const auto policy = GpuExecutionPlanner::plan(1, 1, 1, false, GpuCapability{}, {});
    CHECK(policy.thermal == ThermalMode::Disabled);
}

TEST_CASE("Direct backend remains direct when explicitly requested") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Direct;
    const auto policy = GpuExecutionPlanner::plan(1, 1, 1, false, GpuCapability{}, config);
    CHECK(policy.backend == BackendMode::Direct);
    CHECK(policy.backend_fallback_reason == FallbackReason::None);
}

TEST_CASE("Unsupported explicit Graph capability is a CPU-only fallback") {
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    GpuCapability capability;
    capability.supports_graph = false;

    GpuEngine engine(gpu_test::geometry(), gpu_test::state(1, 1, 1), config, capability);

    CHECK(engine.policy().backend == BackendMode::Fallback);
    CHECK(engine.policy().solver == SolverMode::Eigen);
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().solver == SolverMode::Eigen);
    CHECK(engine.report().static_fallback_reason == FallbackReason::CapabilityUnavailable);
    CHECK_FALSE(engine.context_available());

    engine.step();
    CHECK_FALSE(engine.report().gpu_executed);
    CHECK(engine.pipeline_order().back() == PipelineStage::State);
}

TEST_CASE("GPU execution configuration rejects invalid public values") {
    GpuExecutionConfig config;
    config.device_id = -1;
    CHECK_THROWS_AS(GpuEngine({}, {}, config), std::invalid_argument);

    config = {};
    config.threads_per_block = 0;
    CHECK_THROWS_AS(GpuEngine({}, {}, config), std::invalid_argument);
}

TEST_CASE("GPU execution report distinguishes actual CUDA execution from fallback") {
    ExecutionReport report;
    CHECK_FALSE(report.gpu_executed);
}

TEST_CASE("GPU engine rejects invalid stage voltages at the host boundary") {
    GpuGeometryInput geometry;
    geometry.n_stages = 1;
    geometry.n_filaments = 1;
    geometry.stage_geometry = {0.05};
    geometry.stage_inner_radii = {0.045};
    geometry.stage_outer_radii = {0.055};
    geometry.stage_lengths = {0.01};
    geometry.stage_turns = {20};
    geometry.stage_positions = {0.0};
    geometry.filament_geometry = {0.03};
    geometry.filament_inner_radii = {0.029};
    geometry.filament_outer_radii = {0.031};
    geometry.filament_lengths = {0.01};
    geometry.filament_positions = {0.04};

    GpuEngineState state;
    state.currents = {0.0, 0.0};
    state.m1 = {0.0};
    state.dm1 = {0.0};
    state.active_mask = {1};
    state.trigger_mask = {1};
    state.velocity = {0.0};
    state.position = {0.0};
    state.stage_voltages = {std::numeric_limits<double>::quiet_NaN()};
    state.dt = 1.0e-6;
    state.mass = 1.0;

    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Batched;
    auto valid_geometry = geometry;
    auto valid_state = state;
    auto batched_geometry = geometry;
    auto batched_state = state;
    valid_state.stage_voltages = {0.0};
    batched_state.stage_voltages = {0.0};
    CHECK_THROWS_AS(GpuEngine(std::move(geometry), std::move(state), config),
                    std::invalid_argument);

    GpuEngine engine(std::move(valid_geometry), std::move(valid_state), config);
    CHECK_THROWS_AS(engine.set_stage_voltage(0, std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
    CHECK_THROWS_AS(engine.set_stage_voltage(1, 1.0), std::invalid_argument);

    GpuEngine explicit_fallback(std::move(batched_geometry), std::move(batched_state), config);
    CHECK(explicit_fallback.policy().backend == BackendMode::Fallback);
    CHECK(explicit_fallback.policy().solver == SolverMode::Eigen);
    CHECK_FALSE(explicit_fallback.context_available());
    explicit_fallback.step();
    CHECK(explicit_fallback.report().backend == BackendMode::Fallback);
    CHECK(explicit_fallback.report().solver == SolverMode::Eigen);
    CHECK_FALSE(explicit_fallback.report().gpu_executed);
}

TEST_CASE("CPU fallback rejects singular systems without committing state") {
    using namespace coilgun::simulation::cuda;
    auto geometry = gpu_test::geometry();
    geometry.stage_inductances = {0.0};
    geometry.filament_inductances = {0.0};
    auto state = gpu_test::state(1, 1, 1);
    state.trigger_mask = {0};
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    GpuEngine engine(std::move(geometry), state, config);
    const auto before = engine.state();

    CHECK_THROWS_AS(engine.step(), std::runtime_error);
    CHECK(engine.state().currents == before.currents);
    CHECK(engine.state().m1 == before.m1);
    CHECK(engine.state().dm1 == before.dm1);
    CHECK(engine.state().velocity == before.velocity);
    CHECK(engine.state().position == before.position);
    CHECK(engine.result().completed_steps == 0);
}
