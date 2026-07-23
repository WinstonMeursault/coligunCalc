#include <doctest/doctest.h>
#include "coilgun/simulation/single_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include "gpu_engine_fixture.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <type_traits>
#include <cuda_runtime_api.h>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::SingleStageSim;
using coilgun::simulation::cuda::GpuSingleStageSim;
using coilgun::simulation::cuda::GpuOptLevel;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::ExecutionReport;
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::GpuCapability;
using coilgun::simulation::cuda::ThermalMode;

static double run_cpu_single() {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false);
    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = 20;
    policy.enable_velocity_check = false;
    sim.run(policy);
    return sim.result().summary.muzzle_velocity;
}

static double run_gpu_single(ExecutionReport* report = nullptr) {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend be;
    be.backend = BackendMode::Direct;
    be.use_persistent = false;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Full, be);
    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = 20;
    policy.enable_velocity_check = false;
    sim.run(policy);
    if (report != nullptr) *report = sim.execution_report();
    return sim.result().summary.muzzle_velocity;
}

TEST_CASE("GPU single-stage muzzle velocity matches CPU") {
    double v_cpu = run_cpu_single();
    ExecutionReport report;
    double v_gpu = run_gpu_single(&report);
    if (report.gpu_executed) {
        CHECK(report.gpu_executed);
        CHECK(report.backend == BackendMode::Direct);
    } else {
        CHECK_FALSE(report.gpu_executed);
        CHECK(report.backend == BackendMode::Fallback);
        CHECK_FALSE(report.fallback_reason.empty());
    }
    CHECK(v_gpu == doctest::Approx(v_cpu).epsilon(5e-3));
}

TEST_CASE("GPU single-stage muzzle velocity direction consistent") {
    double v_gpu = run_gpu_single();
    CHECK(v_gpu > 0.0);
}

TEST_CASE("GPU single-stage exposes the unified engine execution report") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        backend);

    CHECK(sim.execution_report().requested_precision ==
          coilgun::simulation::cuda::PrecisionMode::Standard);
    CHECK(sim.execution_report().requested_backend ==
          coilgun::simulation::cuda::BackendMode::Direct);
    const bool direct_or_fallback =
        sim.execution_report().backend == coilgun::simulation::cuda::BackendMode::Direct ||
        sim.execution_report().backend == coilgun::simulation::cuda::BackendMode::Fallback;
    CHECK(direct_or_fallback);
    if (sim.execution_report().backend == coilgun::simulation::cuda::BackendMode::Fallback)
        CHECK_FALSE(sim.execution_report().fallback_reason.empty());
    sim.step();
    const bool gpu_or_fallback = sim.execution_report().gpu_executed ||
        sim.execution_report().backend == coilgun::simulation::cuda::BackendMode::Fallback;
    CHECK(gpu_or_fallback);
}

TEST_CASE("GPU single-stage omitted backend requests Graph without changing shared legacy defaults") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard);

    CHECK(sim.execution_report().requested_backend == BackendMode::Graph);
    CHECK(GpuBackend{}.use_persistent);
}

TEST_CASE("GPU single-stage explicit backend overrides use_persistent") {
    const auto make_sim = [](BackendMode mode, bool use_persistent) {
        DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
        Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                     ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
        auto excitation = std::make_unique<CrowbarExcitation>(450.0, 0.001);
        GpuBackend backend;
        backend.backend = mode;
        backend.use_persistent = use_persistent;
        return GpuSingleStageSim<EulerStepper>(
            coil, arm, std::move(excitation), 1e-6, false,
            coilgun::simulation::cuda::GpuOptLevel::Standard, backend);
    };

    for (const auto mode : {BackendMode::Direct, BackendMode::Graph,
                            BackendMode::Fallback, BackendMode::Persistent}) {
        auto without_legacy = make_sim(mode, false);
        auto with_legacy = make_sim(mode, true);
        const auto& false_report = without_legacy.execution_report();
        const auto& true_report = with_legacy.execution_report();
        CHECK(false_report.requested_backend == mode);
        CHECK(true_report.requested_backend == mode);
        CHECK(false_report.backend == true_report.backend);
        CHECK(false_report.solver == true_report.solver);
        CHECK(false_report.static_fallback_reason == true_report.static_fallback_reason);
        if (mode == BackendMode::Persistent) {
            CHECK(false_report.backend == BackendMode::Fallback);
            CHECK(true_report.backend == BackendMode::Fallback);
            CHECK_FALSE(false_report.gpu_executed);
            CHECK_FALSE(true_report.gpu_executed);
            CHECK_FALSE(false_report.fallback_reason.empty());
            CHECK_FALSE(true_report.fallback_reason.empty());
        }
        without_legacy.step();
        with_legacy.step();
        CHECK(with_legacy.state().currents.isApprox(without_legacy.state().currents, 1e-10));
        CHECK(with_legacy.state().arm_position == doctest::Approx(without_legacy.state().arm_position));
        CHECK(with_legacy.state().arm_velocity == doctest::Approx(without_legacy.state().arm_velocity));
        CHECK(with_legacy.result().history.front().force ==
              doctest::Approx(without_legacy.result().history.front().force));
    }
}

TEST_CASE("GPU single-stage constructor backend overrides conflicting backend fields") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                 ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto excitation = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = true;
    GpuSingleStageSim<EulerStepper> sim(
        coil, arm, std::move(excitation), 1e-6, false,
        coilgun::simulation::cuda::GpuOptLevel::Standard, backend,
        BackendMode::Graph);
    CHECK(sim.execution_report().requested_backend == BackendMode::Graph);
}

TEST_CASE("GPU single-stage proves CUDA execution when a device is available" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        backend);
    sim.step();
    CHECK(sim.execution_report().gpu_executed);
    CHECK(sim.execution_report().backend == coilgun::simulation::cuda::BackendMode::Direct);
}

TEST_CASE("Direct Standard CUDA step numerically matches CPU when a device is available" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    DrivingCoil cpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature cpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     0.0, 0.120, 5, 2, 0.05);
    auto cpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> cpu(cpu_coil, cpu_arm, std::move(cpu_exc), 1e-6, false);

    DrivingCoil gpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature gpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     0.0, 0.120, 5, 2, 0.05);
    auto gpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> gpu(
        gpu_coil, gpu_arm, std::move(gpu_exc), 1e-6, false,
        coilgun::simulation::cuda::GpuOptLevel::Standard, backend);

    for (int step = 0; step < 3; ++step) {
        const auto& cpu_step = cpu.step();
        const auto& gpu_step = gpu.step();
        REQUIRE(gpu.execution_report().gpu_executed);
        REQUIRE(gpu.execution_report().backend == BackendMode::Direct);
        REQUIRE(gpu.state().currents.size() == cpu.state().currents.size());
        for (Eigen::Index i = 0; i < cpu.state().currents.size(); ++i)
            CHECK(gpu.state().currents(i) == doctest::Approx(cpu.state().currents(i)).epsilon(5e-3));
        CHECK(gpu.state().arm_position == doctest::Approx(cpu.state().arm_position).epsilon(5e-3));
        CHECK(gpu.state().arm_velocity == doctest::Approx(cpu.state().arm_velocity).epsilon(5e-3));
        CHECK(gpu_step.force == doctest::Approx(cpu_step.force).epsilon(5e-3));
    }
}

TEST_CASE("GPU execution report records pipeline timing categories" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.enable_profiling = true;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, true,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        backend);

    sim.step();

    const auto& report = sim.execution_report();
    CHECK(report.profiling_enabled);
    CHECK(report.solver_time_ms > 0.0);
    CHECK(report.thermal_time_ms > 0.0);
    if (report.gpu_executed) {
        CHECK(report.gpu_time_ms > 0.0);
        CHECK(report.transfer_time_ms > 0.0);
        CHECK(report.gpu_time_ms >= report.transfer_time_ms);
    }
}

TEST_CASE("GpuEngine keeps an unsupported Graph request CPU-only" ) {
    using namespace coilgun::simulation::cuda;
    GpuExecutionConfig config;
    config.backend = BackendMode::Graph;
    GpuCapability capability;
    capability.supports_graph = false;

    GpuEngine engine(gpu_test::geometry(), gpu_test::state(1, 1, 1), config, capability);

    CHECK(engine.policy().backend == BackendMode::Fallback);
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().static_fallback_reason == FallbackReason::CapabilityUnavailable);
    CHECK_FALSE(engine.context_available());
    engine.step();
    CHECK_FALSE(engine.report().gpu_executed);
}

TEST_CASE("CUDA allocation failure locks the engine to an Eigen CPU fallback" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    using namespace coilgun::simulation::cuda;
    GpuExecutionConfig config;
    config.backend = BackendMode::Direct;
    config.solver = SolverMode::Eigen;
    detail::GpuEngineFaultInjection fault;
    fault.fail_allocation = true;

    GpuEngine engine(gpu_test::geometry(), gpu_test::state(1, 1, 1), config, {}, fault);

    CHECK(engine.policy().backend == BackendMode::Fallback);
    CHECK(engine.policy().solver == SolverMode::Eigen);
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().solver == SolverMode::Eigen);
    CHECK(engine.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK_FALSE(engine.report().fallback_reason.empty());
    CHECK_FALSE(engine.context_available());
    CHECK(engine.solver_workspace_initialized());
    CHECK(engine.device_allocation_count() >= 3);
    CHECK(engine.device_buffer_addresses().empty());

    CHECK_NOTHROW(engine.step());
    CHECK_FALSE(engine.report().gpu_executed);
    CHECK(engine.result().completed_steps == 1);

    CHECK_NOTHROW(engine.step());
    CHECK_FALSE(engine.report().gpu_executed);
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().solver == SolverMode::Eigen);
    CHECK(engine.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK(engine.result().completed_steps == 2);
}

TEST_CASE("CUDA device initialization runtime failure locks a usable CPU fallback") {
    using namespace coilgun::simulation::cuda;
    auto failed_state = gpu_test::state(1, 1, 1);
    failed_state.stage_voltages = {450.0};
    auto fallback_state = failed_state;
    GpuExecutionConfig failed_config;
    failed_config.backend = BackendMode::Direct;
    detail::GpuEngineFaultInjection fault;
    fault.fail_device_initialization = true;
    GpuExecutionConfig fallback_config;
    fallback_config.backend = BackendMode::Fallback;

    GpuEngine failed(gpu_test::geometry(), std::move(failed_state), failed_config, {}, fault);
    GpuEngine fallback(gpu_test::geometry(), std::move(fallback_state), fallback_config);

    CHECK(failed.report().backend == BackendMode::Fallback);
    CHECK(failed.report().solver == SolverMode::Eigen);
    CHECK(failed.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK_FALSE(failed.report().fallback_reason.empty());
    CHECK_FALSE(failed.context_available());
    CHECK(failed.solver_workspace_initialized());

    failed.step();
    fallback.step();
    CHECK(failed.state().currents == fallback.state().currents);
    CHECK(failed.state().position == fallback.state().position);
    CHECK(failed.state().velocity == fallback.state().velocity);
    CHECK_FALSE(failed.report().gpu_executed);
}

TEST_CASE("GPU single-stage reports fallback when CUDA is explicitly unavailable" *
          doctest::skip(coilgun::simulation::cuda::cuda_device_available())) {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        backend);
    sim.step();
    CHECK_FALSE(sim.execution_report().gpu_executed);
    CHECK(sim.execution_report().backend == coilgun::simulation::cuda::BackendMode::Fallback);
    CHECK(sim.execution_report().fallback_count > 0);
}

TEST_CASE("GPU single-stage preserves backend execution settings") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.device_id = 0;
    backend.threads_per_block = 256;
    backend.enable_profiling = true;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        backend);
    CHECK(sim.execution_report().device_id == backend.device_id);
    CHECK(sim.execution_report().threads_per_block == backend.threads_per_block);
    CHECK(sim.execution_report().profiling_enabled);
}

TEST_CASE("GPU engine restores the caller CUDA device across lifecycle boundaries" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    int count = 0;
    REQUIRE(cudaGetDeviceCount(&count) == cudaSuccess);
    int original = -1;
    REQUIRE(cudaGetDevice(&original) == cudaSuccess);
    const int selected = count >= 2 ? (original == 0 ? 1 : 0) : original;

    using namespace coilgun::simulation::cuda;
    GpuExecutionConfig config;
    config.backend = BackendMode::Direct;
    config.solver = SolverMode::Eigen;
    config.device_id = selected;
    {
        GpuEngine engine(gpu_test::geometry(), gpu_test::state(1, 1, 1), config);
        int after_init = -1;
        REQUIRE(cudaGetDevice(&after_init) == cudaSuccess);
        CHECK(after_init == original);
        engine.step();
        int after_step = -1;
        REQUIRE(cudaGetDevice(&after_step) == cudaSuccess);
        CHECK(after_step == original);
        engine.shutdown();
        int after_shutdown = -1;
        REQUIRE(cudaGetDevice(&after_shutdown) == cudaSuccess);
        CHECK(after_shutdown == original);
    }
    int after_destroy = -1;
    REQUIRE(cudaGetDevice(&after_destroy) == cudaSuccess);
    CHECK(after_destroy == original);
}

TEST_CASE("GPU single-stage graph request uses CUDA when available and matches fallback") {
    DrivingCoil graph_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature graph_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto graph_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend graph_backend;
    graph_backend.backend = BackendMode::Graph;
    graph_backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> graph(graph_coil, graph_arm, std::move(graph_exc), 1e-6,
                                          false, coilgun::simulation::cuda::GpuOptLevel::Standard,
                                          graph_backend);

    DrivingCoil fallback_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature fallback_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto fallback_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend fallback_backend;
    fallback_backend.backend = BackendMode::Fallback;
    fallback_backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> fallback(fallback_coil, fallback_arm, std::move(fallback_exc), 1e-6,
                                             false, coilgun::simulation::cuda::GpuOptLevel::Standard,
                                             fallback_backend);

    DrivingCoil direct_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature direct_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto direct_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend direct_backend;
    direct_backend.backend = BackendMode::Direct;
    direct_backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> direct(direct_coil, direct_arm, std::move(direct_exc), 1e-6,
                                           false, coilgun::simulation::cuda::GpuOptLevel::Standard,
                                           direct_backend);

    graph.step();
    fallback.step();
    direct.step();
    CHECK(graph.state().currents.isApprox(fallback.state().currents, 1e-10));
    CHECK(graph.state().arm_position == doctest::Approx(fallback.state().arm_position));
    CHECK(graph.state().arm_velocity == doctest::Approx(fallback.state().arm_velocity));
    REQUIRE(graph.result().history.size() == 1);
    REQUIRE(fallback.result().history.size() == 1);
    CHECK(graph.result().history.front().force == doctest::Approx(fallback.result().history.front().force));
    CHECK(graph.state().currents.isApprox(direct.state().currents, 1e-10));
    CHECK(graph.result().history.front().force == doctest::Approx(direct.result().history.front().force));
    if (graph.execution_report().gpu_executed) {
        CHECK(graph.execution_report().gpu_executed);
        CHECK(graph.execution_report().backend == BackendMode::Graph);
        CHECK(graph.execution_report().graph_rebuild_count == 1);
    } else {
        CHECK_FALSE(graph.execution_report().gpu_executed);
        CHECK(graph.execution_report().backend == BackendMode::Fallback);
        CHECK_FALSE(graph.execution_report().fallback_reason.empty());
        CHECK(graph.execution_report().graph_rebuild_count == 0);
    }
    CHECK_FALSE(fallback.execution_report().gpu_executed);
    CHECK(fallback.execution_report().backend == BackendMode::Fallback);
}

TEST_CASE("Explicit fallback with batched solver stays CPU-only") {
    using namespace coilgun::simulation::cuda;
    GpuExecutionConfig config;
    config.backend = BackendMode::Fallback;
    config.solver = SolverMode::Batched;

    GpuEngine engine(gpu_test::geometry(), gpu_test::state(1, 1, 1), config);
    CHECK(engine.policy().backend == BackendMode::Fallback);
    CHECK(engine.policy().solver == SolverMode::Eigen);
    CHECK_FALSE(engine.context_available());
    engine.step();
    CHECK(engine.report().backend == BackendMode::Fallback);
    CHECK(engine.report().solver == SolverMode::Eigen);
    CHECK_FALSE(engine.report().gpu_executed);
}

TEST_CASE("Injected GPU failure commits one complete CPU-equivalent step and locks fallback" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    using namespace coilgun::simulation::cuda;
    auto injected_state = gpu_test::state(1, 1, 1);
    injected_state.stage_voltages = {450.0};
    auto fallback_state = injected_state;
    GpuExecutionConfig injected_config;
    injected_config.backend = BackendMode::Direct;
    injected_config.solver = SolverMode::Eigen;
    detail::GpuEngineFaultInjection fault;
    fault.fail_after_mutual = true;
    GpuExecutionConfig fallback_config;
    fallback_config.backend = BackendMode::Fallback;
    fallback_config.solver = SolverMode::Eigen;
    GpuEngine injected(gpu_test::geometry(), std::move(injected_state), injected_config, {}, fault);
    GpuEngine fallback(gpu_test::geometry(), std::move(fallback_state), fallback_config);
    const auto initial = injected.state();

    injected.step();
    fallback.step();
    CHECK(injected.state().currents == fallback.state().currents);
    CHECK(injected.state().m1 == fallback.state().m1);
    CHECK(injected.state().dm1 == fallback.state().dm1);
    CHECK(injected.state().velocity == fallback.state().velocity);
    CHECK(injected.state().position == fallback.state().position);
    CHECK(injected.state().current_derivatives == fallback.state().current_derivatives);
    CHECK(injected.state().currents != initial.currents);
    CHECK(injected.result().completed_steps == 1);
    CHECK_FALSE(injected.report().gpu_executed);
    CHECK(injected.report().backend == BackendMode::Fallback);
    CHECK(injected.report().solver == SolverMode::Eigen);
    CHECK(injected.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK(injected.report().fallback_count > 0);
    CHECK_FALSE(injected.report().fallback_reason.empty());
    const auto fallback_count = injected.report().fallback_count;

    injected.step();
    fallback.step();
    CHECK(injected.state().currents == fallback.state().currents);
    CHECK(injected.state().velocity == fallback.state().velocity);
    CHECK(injected.state().position == fallback.state().position);
    CHECK(injected.result().completed_steps == 2);
    CHECK_FALSE(injected.report().gpu_executed);
    CHECK(injected.report().backend == BackendMode::Fallback);
    CHECK(injected.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK(injected.report().fallback_count == fallback_count);
}

TEST_CASE("Injected Graph capture failure commits one CPU step and retains locked fallback" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    using namespace coilgun::simulation::cuda;
    auto graph_state = gpu_test::state(1, 1, 1);
    graph_state.stage_voltages = {450.0};
    auto fallback_state = graph_state;
    GpuExecutionConfig graph_config;
    graph_config.backend = BackendMode::Graph;
    graph_config.solver = SolverMode::Eigen;
    detail::GpuEngineFaultInjection fault;
    fault.fail_graph_capture = true;
    GpuExecutionConfig fallback_config;
    fallback_config.backend = BackendMode::Fallback;
    fallback_config.solver = SolverMode::Eigen;
    GpuEngine graph(gpu_test::geometry(), std::move(graph_state), graph_config, {}, fault);
    GpuEngine fallback(gpu_test::geometry(), std::move(fallback_state), fallback_config);

    graph.step();
    fallback.step();
    CHECK(graph.state().currents == fallback.state().currents);
    CHECK(graph.state().m1 == fallback.state().m1);
    CHECK(graph.state().dm1 == fallback.state().dm1);
    CHECK(graph.state().position == fallback.state().position);
    CHECK(graph.state().velocity == fallback.state().velocity);
    CHECK(graph.state().current_derivatives == fallback.state().current_derivatives);
    CHECK(graph.report().backend == BackendMode::Fallback);
    CHECK(graph.report().solver == SolverMode::Eigen);
    CHECK(graph.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK(graph.report().fallback_count > 0);
    CHECK_FALSE(graph.report().fallback_reason.empty());
    CHECK_FALSE(graph.report().gpu_executed);
    CHECK(graph.report().graph_rebuild_count == 0);
    const auto fallback_count = graph.report().fallback_count;

    graph.step();
    fallback.step();
    CHECK(graph.state().currents == fallback.state().currents);
    CHECK(graph.state().position == fallback.state().position);
    CHECK(graph.state().velocity == fallback.state().velocity);
    CHECK(graph.report().fallback_count == fallback_count);
    CHECK(graph.report().runtime_fallback_reason == FallbackReason::RuntimeFailure);
    CHECK_FALSE(graph.report().gpu_executed);
}

TEST_CASE("GPU single-stage persistent request safely falls back and matches explicit fallback") {
    DrivingCoil persistent_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature persistent_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto persistent_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend persistent_backend;
    persistent_backend.backend = BackendMode::Persistent;
    persistent_backend.use_persistent = true;
    GpuSingleStageSim<EulerStepper> persistent(persistent_coil, persistent_arm, std::move(persistent_exc), 1e-6,
                                               false, coilgun::simulation::cuda::GpuOptLevel::Standard,
                                               persistent_backend);

    DrivingCoil fallback_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature fallback_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto fallback_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend fallback_backend;
    fallback_backend.backend = BackendMode::Fallback;
    GpuSingleStageSim<EulerStepper> fallback(fallback_coil, fallback_arm, std::move(fallback_exc), 1e-6,
                                             false, coilgun::simulation::cuda::GpuOptLevel::Standard,
                                             fallback_backend);

    persistent.step();
    fallback.step();
    CHECK(persistent.state().currents.isApprox(fallback.state().currents, 1e-10));
    CHECK(persistent.state().arm_position == doctest::Approx(fallback.state().arm_position));
    CHECK(persistent.state().arm_velocity == doctest::Approx(fallback.state().arm_velocity));
    REQUIRE(persistent.result().history.size() == 1);
    REQUIRE(fallback.result().history.size() == 1);
    CHECK(persistent.result().history.front().force == doctest::Approx(fallback.result().history.front().force));
    CHECK(persistent.execution_report().requested_backend == BackendMode::Persistent);
    CHECK_FALSE(persistent.execution_report().gpu_executed);
    CHECK(persistent.execution_report().backend == BackendMode::Fallback);
    CHECK(persistent.execution_report().fallback_count > 0);
    CHECK_FALSE(persistent.execution_report().fallback_reason.empty());
}

TEST_CASE("GPU single-stage Full and Aggressive apply the canonical distant mutual cutoff") {
    constexpr double dt = 1e-6;
    constexpr double distant_coil_position = 1.0;

    auto make_sim = [distant_coil_position](GpuOptLevel opt_level,
                                             BackendMode backend_mode) {
        DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref,
                         1e-6, 0.7, distant_coil_position);
        Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                     ALUMINUM.density, 0.0, 0.120, 5, 2, 0.0);
        auto excitation = std::make_unique<CrowbarExcitation>(450.0, 0.001);
        GpuBackend backend;
        backend.backend = backend_mode;
        backend.use_persistent = false;
        return GpuSingleStageSim<EulerStepper>(
            coil, arm, std::move(excitation), dt, false, opt_level, backend);
    };

    for (const auto opt_level : {GpuOptLevel::Full, GpuOptLevel::Aggressive}) {
        auto optimized = make_sim(opt_level, BackendMode::Direct);
        auto fallback = make_sim(opt_level, BackendMode::Fallback);
        optimized.step();
        fallback.step();

        REQUIRE(optimized.state().currents(0) != doctest::Approx(0.0));
        for (Eigen::Index index = 1; index < optimized.state().currents.size(); ++index)
            CHECK(optimized.state().currents(index) == doctest::Approx(0.0));
        CHECK(optimized.result().history.front().force == doctest::Approx(0.0));
        CHECK(optimized.state().arm_position == doctest::Approx(0.0));
        CHECK(optimized.state().arm_velocity == doctest::Approx(0.0));
        CHECK(optimized.state().currents.isApprox(fallback.state().currents, 1e-10));
        CHECK(optimized.result().history.front().force ==
              doctest::Approx(fallback.result().history.front().force));
        CHECK(optimized.state().arm_position == doctest::Approx(fallback.state().arm_position));
        CHECK(optimized.state().arm_velocity == doctest::Approx(fallback.state().arm_velocity));
    }
}

TEST_CASE("GPU single-stage rejects RK4 instead of silently using Euler") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<coilgun::simulation::RK4Stepper> sim(coil, arm, std::move(exc), 1e-6);
    CHECK_THROWS_WITH(sim.step(), "RK4Stepper is not supported by GpuSingleStageSim");
}

TEST_CASE("GPU single-stage records state and force on every step") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6);
    const auto initial = sim.state();
    sim.step();
    const auto first_velocity = sim.state().arm_velocity;
    REQUIRE(sim.result().history.size() == 1);
    CHECK(sim.state().currents.size() == initial.currents.size());
    DrivingCoil cpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature cpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto cpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> cpu(cpu_coil, cpu_arm, std::move(cpu_exc), 1e-6, false);
    cpu.step();
    CHECK(sim.result().history.front().force == doctest::Approx(cpu.result().history.front().force).epsilon(5e-3));
    CHECK(std::signbit(sim.result().history.front().force) ==
          std::signbit(cpu.result().history.front().force));
    CHECK(std::abs(sim.result().history.front().force) > 0.0);
    CHECK(sim.state().currents(0) != initial.currents(0));
    sim.step();
    REQUIRE(sim.result().history.size() == 2);
    CHECK(sim.state().arm_velocity != first_velocity);
    CHECK(std::isfinite(sim.result().history.back().force));
}

TEST_CASE("GPU single-stage engine advances one step") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        {});
    sim.step();
    CHECK(sim.step_count() == 1);
}

TEST_CASE("GPU single-stage honors a short termination policy") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        {});
    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = 3;
    policy.enable_velocity_check = false;
    sim.run(policy);
    CHECK(sim.step_count() == 3);
}

TEST_CASE("GPU single-stage summary max_force matches CPU peak absolute force") {
    DrivingCoil cpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature cpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     0.0, 0.120, 5, 2, -0.1);
    auto cpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> cpu(cpu_coil, cpu_arm, std::move(cpu_exc), 1e-6, false);

    DrivingCoil gpu_coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature gpu_arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     0.0, 0.120, 5, 2, -0.1);
    auto gpu_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuSingleStageSim<EulerStepper> gpu(
        gpu_coil, gpu_arm, std::move(gpu_exc), 1e-6, false,
        coilgun::simulation::cuda::GpuOptLevel::Standard, backend);

    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = 3;
    policy.enable_velocity_check = false;
    cpu.run(policy);
    gpu.run(policy);

    REQUIRE(cpu.result().history.size() == gpu.result().history.size());
    double expected_max_force = 0.0;
    bool saw_negative_force = false;
    for (std::size_t i = 0; i < cpu.result().history.size(); ++i) {
        const double cpu_force = cpu.result().history[i].force;
        const double gpu_force = gpu.result().history[i].force;
        expected_max_force = std::max(expected_max_force, std::abs(cpu_force));
        saw_negative_force = saw_negative_force || cpu_force < 0.0;
        CHECK(gpu_force == doctest::Approx(cpu_force).epsilon(5e-3));
    }
    REQUIRE(saw_negative_force);
    CHECK(cpu.result().summary.max_force == doctest::Approx(expected_max_force));
    CHECK(gpu.result().summary.max_force == doctest::Approx(expected_max_force).epsilon(5e-3));
    CHECK(gpu.result().summary.max_force ==
          doctest::Approx(cpu.result().summary.max_force).epsilon(5e-3));
}

TEST_CASE("GPU single-stage thermal state and reset remain wrapper-compatible") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, true,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard,
                                        {});
    static_assert(std::is_same_v<decltype(sim.filament_resistances()), std::vector<double>>);
    sim.step();
    REQUIRE(sim.state().filament_temperatures.size() == arm.total_filaments());
    REQUIRE(sim.state().filament_temperatures(0) == doctest::Approx(T_REFERENCE));
    sim.step();
    REQUIRE(sim.state().filament_temperatures(0) > T_REFERENCE);
    REQUIRE(sim.filament_resistances().size() == static_cast<std::size_t>(arm.total_filaments()));
    REQUIRE(sim.filament_resistances()[0] > arm.resistances()[0]);
    CHECK(sim.execution_report().thermal != coilgun::simulation::cuda::ThermalMode::Disabled);
    const auto fallback_count = sim.execution_report().fallback_count;
    const auto graph_rebuild_count = sim.execution_report().graph_rebuild_count;
    sim.reset();
    CHECK(sim.step_count() == 0);
    CHECK(sim.result().history.empty());
    CHECK(sim.state().currents.isZero(0.0));
    CHECK(sim.state().arm_position == doctest::Approx(arm.position()));
    CHECK(sim.state().filament_temperatures(0) == doctest::Approx(T_REFERENCE));
    REQUIRE(sim.filament_resistances().size() == arm.resistances().size());
    for (std::size_t i = 0; i < arm.resistances().size(); ++i)
        CHECK(sim.filament_resistances()[i] == doctest::Approx(arm.resistances()[i]));
    CHECK(sim.execution_report().fallback_count == fallback_count);
    CHECK(sim.execution_report().graph_rebuild_count == graph_rebuild_count);
}

TEST_CASE("GPU single-stage reset reruns the same first state, history, thermal, and diagnostics") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, true,
                                        coilgun::simulation::cuda::GpuOptLevel::Standard, backend);

    sim.step();
    const auto first_state = sim.state();
    const auto first_history = sim.result().history.front();
    sim.step();
    const auto cumulative_report = sim.execution_report();
    sim.reset();

    CHECK(sim.step_count() == 0);
    CHECK(sim.result().history.empty());
    CHECK(sim.execution_report().backend == cumulative_report.backend);
    CHECK(sim.execution_report().solver == cumulative_report.solver);
    CHECK(sim.execution_report().thermal == cumulative_report.thermal);
    CHECK(sim.execution_report().fallback_count == cumulative_report.fallback_count);
    CHECK(sim.execution_report().graph_rebuild_count == cumulative_report.graph_rebuild_count);
    CHECK(sim.execution_report().gpu_executed == cumulative_report.gpu_executed);
    CHECK(sim.execution_report().fallback_reason == cumulative_report.fallback_reason);
    CHECK(sim.execution_report().runtime_fallback_reason == cumulative_report.runtime_fallback_reason);
    CHECK(sim.execution_report().gpu_time_ms == doctest::Approx(cumulative_report.gpu_time_ms));
    CHECK(sim.execution_report().solver_time_ms == doctest::Approx(cumulative_report.solver_time_ms));
    CHECK(sim.execution_report().thermal_time_ms == doctest::Approx(cumulative_report.thermal_time_ms));
    CHECK(sim.execution_report().transfer_time_ms == doctest::Approx(cumulative_report.transfer_time_ms));
    CHECK(sim.execution_report().max_condition_estimate == doctest::Approx(cumulative_report.max_condition_estimate));
    CHECK(sim.execution_report().calibrated == cumulative_report.calibrated);
    CHECK(sim.execution_report().precision_fallback == cumulative_report.precision_fallback);
    CHECK(sim.execution_report().metadata_conflict == cumulative_report.metadata_conflict);

    sim.step();
    CHECK(sim.state().currents.isApprox(first_state.currents, 1e-12));
    CHECK(sim.state().filament_temperatures.isApprox(first_state.filament_temperatures, 1e-12));
    CHECK(sim.state().arm_position == doctest::Approx(first_state.arm_position));
    CHECK(sim.state().arm_velocity == doctest::Approx(first_state.arm_velocity));
    REQUIRE(sim.result().history.size() == 1);
    CHECK(sim.result().history.front().time == doctest::Approx(first_history.time));
    CHECK(sim.result().history.front().force == doctest::Approx(first_history.force));
    CHECK(sim.result().history.front().filament_temperatures == first_history.filament_temperatures);
    CHECK(sim.execution_report().backend == cumulative_report.backend);
    CHECK(sim.execution_report().solver == cumulative_report.solver);
    CHECK(sim.execution_report().thermal == ThermalMode::Cpu);
    CHECK(sim.execution_report().fallback_count == cumulative_report.fallback_count);
    CHECK(sim.execution_report().gpu_executed == cumulative_report.gpu_executed);
}

TEST_CASE("GPU single-stage rejects invalid inputs before CUDA selection") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);

    auto null_excitation = std::unique_ptr<coilgun::simulation::Excitation>{};
    CHECK_THROWS_AS((GpuSingleStageSim<EulerStepper>(coil, arm, std::move(null_excitation), 1e-6)),
                    std::invalid_argument);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    CHECK_THROWS_AS((GpuSingleStageSim<EulerStepper>(coil, arm, std::move(exc), 0.0)),
                    std::invalid_argument);
    auto bad_threads_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend bad_threads;
    bad_threads.threads_per_block = 3;
    CHECK_THROWS_AS((GpuSingleStageSim<EulerStepper>(coil, arm, std::move(bad_threads_exc), 1e-6,
                                                     false, coilgun::simulation::cuda::GpuOptLevel::Full,
                                                     bad_threads)), std::invalid_argument);

    auto valid_threads_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuBackend valid_threads;
    valid_threads.backend = BackendMode::Fallback;
    valid_threads.threads_per_block = 128;
    CHECK_NOTHROW((GpuSingleStageSim<EulerStepper>(coil, arm, std::move(valid_threads_exc), 1e-6,
                                                   false, coilgun::simulation::cuda::GpuOptLevel::Full,
                                                   valid_threads)));

    auto bad_opt_exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    CHECK_THROWS_AS((GpuSingleStageSim<EulerStepper>(
                        coil, arm, std::move(bad_opt_exc), 1e-6, false,
                        static_cast<coilgun::simulation::cuda::GpuOptLevel>(99),
                        valid_threads)), std::invalid_argument);
}
