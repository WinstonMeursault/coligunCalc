/**
 * @file test_gpu_vs_cpu_multi.cpp
 * @brief End-to-end test: GPU vs CPU multi-stage simulation.
 * @author Winston Meursault
 */

#include <doctest/doctest.h>
#include "coilgun/simulation/multi_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/simulation/termination.hpp"
#include "gpu_numerical_tolerances.hpp"
#include <memory>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::MultiStageSim;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::cuda::GpuMultiStageSim;
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::GpuOptLevel;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::ExecutionReport;
using coilgun::simulation::cuda::ThermalMode;

namespace {

using gpu_test::numerically_equal;
using gpu_test::tolerance_for;

struct MultiCase {
    std::vector<DrivingCoil> coils;
    Armature armature;
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excitations;
    std::vector<TriggerConfig> triggers;
};

MultiCase make_multi_case() {
    MultiCase value{
        {
            DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.015),
            DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.085)
        },
        Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                 4.0, 0.120, 5, 2, 0.0),
        {},
        {{TriggerMode::Position, 8e-6}}
    };
    value.excitations.push_back(std::make_unique<CrowbarExcitation>(480.0, 0.001));
    value.excitations.push_back(std::make_unique<CrowbarExcitation>(320.0, 0.0008));
    return value;
}

coilgun::simulation::TerminationPolicy short_policy(int max_steps = 8) {
    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = max_steps;
    policy.enable_velocity_check = false;
    policy.enable_bound_check = false;
    return policy;
}

} // namespace

static double run_cpu_multi() {
    auto value = make_multi_case();
    MultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                    std::move(value.excitations), std::move(value.triggers), 1e-6,
                                    false, coilgun::simulation::OptimizationLevel::Reference);
    sim.run(short_policy());
    return sim.result().summary.muzzle_velocity;
}

static double run_gpu_multi(ExecutionReport* report = nullptr) {
    auto value = make_multi_case();
    GpuBackend be;
    be.backend = BackendMode::Auto;
    be.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                       std::move(value.excitations), std::move(value.triggers),
                                       1e-6, false, GpuOptLevel::Standard, be);
    sim.run(short_policy());
    if (report != nullptr) *report = sim.execution_report();
    return sim.result().summary.muzzle_velocity;
}

TEST_CASE("GPU multi-stage primary muzzle velocity comparison requires CUDA" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    double v_cpu = run_cpu_multi();
    ExecutionReport report;
    double v_gpu = run_gpu_multi(&report);
    REQUIRE(report.gpu_executed);
    REQUIRE(report.backend == BackendMode::Direct);
    CHECK(gpu_test::numerically_equal(v_gpu, v_cpu,
                                      gpu_test::tolerance_for(GpuOptLevel::Standard)));
}

TEST_CASE("GPU multi-stage explicit backend overrides use_persistent") {
    auto make_sim = [](BackendMode backend_mode, bool use_persistent,
                       BackendMode explicit_mode = BackendMode::Auto) {
        auto value = make_multi_case();
        GpuBackend backend;
        backend.backend = backend_mode;
        backend.use_persistent = use_persistent;
        return std::make_unique<GpuMultiStageSim<EulerStepper>>(
            std::move(value.coils), value.armature, std::move(value.excitations),
            std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend,
            explicit_mode);
    };

    for (const auto mode : {BackendMode::Direct, BackendMode::Graph,
                            BackendMode::Fallback, BackendMode::Persistent}) {
        auto without_legacy = make_sim(mode, false);
        auto with_legacy = make_sim(mode, true);
        const auto& false_report = without_legacy->execution_report();
        const auto& true_report = with_legacy->execution_report();
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
        without_legacy->step();
        with_legacy->step();
        CHECK(with_legacy->state().currents.isApprox(without_legacy->state().currents, 1e-10));
        CHECK(with_legacy->state().arm_position == doctest::Approx(without_legacy->state().arm_position));
        CHECK(with_legacy->state().arm_velocity == doctest::Approx(without_legacy->state().arm_velocity));
        CHECK(with_legacy->result().history.front().state.force ==
              doctest::Approx(without_legacy->result().history.front().state.force));
    }

    auto legacy_direct = make_sim(BackendMode::Auto, false);
    auto legacy_persistent = make_sim(BackendMode::Auto, true);
    CHECK(legacy_direct->execution_report().requested_backend == BackendMode::Direct);
    CHECK(legacy_persistent->execution_report().requested_backend == BackendMode::Persistent);

    auto explicit_graph = make_sim(BackendMode::Auto, true, BackendMode::Graph);
    CHECK(explicit_graph->execution_report().requested_backend == BackendMode::Graph);
}

TEST_CASE("GPU multi-stage constructor explicit backend wins over conflicting backend metadata") {
    auto make_sim = [](BackendMode backend_mode, bool use_persistent, BackendMode explicit_mode) {
        auto value = make_multi_case();
        GpuBackend backend;
        backend.backend = backend_mode;
        backend.use_persistent = use_persistent;
        return std::make_unique<GpuMultiStageSim<EulerStepper>>(
            std::move(value.coils), value.armature, std::move(value.excitations),
            std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend,
            explicit_mode);
    };

    auto explicit_graph = make_sim(BackendMode::Direct, true, BackendMode::Graph);
    auto explicit_direct = make_sim(BackendMode::Graph, true, BackendMode::Direct);
    auto explicit_fallback = make_sim(BackendMode::Graph, true, BackendMode::Fallback);

    CHECK(explicit_graph->execution_report().requested_backend == BackendMode::Graph);
    CHECK(explicit_direct->execution_report().requested_backend == BackendMode::Direct);
    CHECK(explicit_fallback->execution_report().requested_backend == BackendMode::Fallback);
}

TEST_CASE("GPU multi-stage Persistent request reports its separate safe fallback") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Auto;
    backend.use_persistent = true;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    sim.step();
    CHECK(sim.execution_report().requested_backend == BackendMode::Persistent);
    CHECK_FALSE(sim.execution_report().gpu_executed);
    CHECK(sim.execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(sim.execution_report().fallback_reason.empty());
}

TEST_CASE("GPU multi-stage omitted backend executes the supported Direct path on CUDA" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    auto value = make_multi_case();
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6);

    CHECK(sim.execution_report().requested_backend == BackendMode::Direct);
    sim.step();
    REQUIRE(sim.execution_report().gpu_executed);
    CHECK(sim.execution_report().backend == BackendMode::Direct);
}

TEST_CASE("GPU multi-stage exposes an honest engine report and preserves stage masks") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                       std::move(value.excitations), std::move(value.triggers),
                                       1e-6, false, GpuOptLevel::Standard, backend);

    CHECK(sim.execution_report().requested_backend == BackendMode::Fallback);
    CHECK(sim.execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(sim.execution_report().gpu_executed);
    sim.step();
    REQUIRE(sim.result().history.size() == 1);
    CHECK(sim.result().history.front().coil_currents.size() == 2);
    CHECK(sim.result().history.front().coil_currents[1] == doctest::Approx(0.0));
    sim.step();
    CHECK(sim.result().history.back().coil_currents[1] == doctest::Approx(0.0));
    sim.step();
    CHECK(sim.result().history.back().coil_currents[1] != 0.0);
}

TEST_CASE("GPU multi-stage position and time triggers activate the expected stages") {
    constexpr double position_trigger = 8e-6;
    auto position_case = make_multi_case();
    const auto position_triggers = position_case.triggers;
    GpuBackend fallback;
    fallback.backend = BackendMode::Fallback;
    fallback.use_persistent = false;
    GpuMultiStageSim<EulerStepper> position_sim(
        std::move(position_case.coils), position_case.armature,
        std::move(position_case.excitations), position_triggers,
        1e-6, false, GpuOptLevel::Standard, fallback);
    position_sim.step();
    REQUIRE(position_sim.result().history.size() == 1);
    REQUIRE(position_sim.result().history.front().coil_currents.size() == 2);
    CHECK(position_sim.result().history.front().coil_currents[1] == doctest::Approx(0.0));
    CHECK(position_sim.state().arm_position < position_trigger);
    position_sim.step();
    CHECK(position_sim.result().history.back().coil_currents[1] == doctest::Approx(0.0));
    position_sim.step();
    CHECK(position_sim.result().history.back().coil_currents[1] != 0.0);

    auto time_case = make_multi_case();
    time_case.triggers = {{TriggerMode::TimeDelay, 1e-6}};
    GpuMultiStageSim<EulerStepper> time_sim(
        std::move(time_case.coils), time_case.armature,
        std::move(time_case.excitations), std::move(time_case.triggers),
        1e-6, false, GpuOptLevel::Standard, fallback);
    time_sim.step();
    REQUIRE(time_sim.result().history.size() == 1);
    CHECK(time_sim.result().history.front().coil_currents[1] == doctest::Approx(0.0));
    time_sim.step();
    REQUIRE(time_sim.result().history.size() == 2);
    CHECK(time_sim.result().history.back().coil_currents[1] != 0.0);
}

TEST_CASE("GPU multi-stage triggers capture the exact pre-step boundary") {
    auto zero_delay_case = make_multi_case();
    zero_delay_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> zero_delay(
        std::move(zero_delay_case.coils), zero_delay_case.armature,
        std::move(zero_delay_case.excitations), std::move(zero_delay_case.triggers),
        1e-6, false, GpuOptLevel::Standard, backend);
    const double zero_delay_position = zero_delay.state().arm_position;
    zero_delay.step();
    REQUIRE(zero_delay.result().history.back().coil_currents[1] != 0.0);
    zero_delay.run(short_policy(0));
    REQUIRE(zero_delay.result().summary.per_stage.size() == 2);
    CHECK(zero_delay.result().summary.per_stage[1].trigger_position ==
          doctest::Approx(zero_delay_position).epsilon(1e-12));

    auto crossed_case = make_multi_case();
    crossed_case.triggers = {{TriggerMode::Position, 8e-6}};
    GpuMultiStageSim<EulerStepper> crossed(
        std::move(crossed_case.coils), crossed_case.armature,
        std::move(crossed_case.excitations), std::move(crossed_case.triggers),
        1e-6, false, GpuOptLevel::Standard, backend);
    double crossed_position = 0.0;
    bool triggered = false;
    for (int i = 0; i < 100; ++i) {
        crossed_position = crossed.state().arm_position;
        crossed.step();
        if (crossed.result().history.back().coil_currents[1] != 0.0) {
            triggered = true;
            break;
        }
    }
    REQUIRE(triggered);
    crossed.run(short_policy(0));
    REQUIRE(crossed.result().summary.per_stage.size() == 2);
    CHECK(crossed.result().summary.per_stage[1].trigger_position ==
          doctest::Approx(crossed_position).epsilon(1e-12));
}

TEST_CASE("GPU multi-stage trigger configuration rejects malformed values") {
    auto construct = [](TriggerConfig config) {
        auto value = make_multi_case();
        GpuBackend backend;
        backend.backend = BackendMode::Fallback;
        backend.use_persistent = false;
        GpuMultiStageSim<EulerStepper> sim(
            std::move(value.coils), value.armature, std::move(value.excitations),
            {config}, 1e-6, false, GpuOptLevel::Standard, backend);
    };
    CHECK_THROWS_AS(construct({TriggerMode::Position,
                               std::numeric_limits<double>::quiet_NaN()}),
                    std::invalid_argument);
    CHECK_THROWS_AS(construct({TriggerMode::TimeDelay,
                               std::numeric_limits<double>::quiet_NaN()}),
                    std::invalid_argument);
    CHECK_THROWS_AS(construct({TriggerMode::TimeDelay, -1.0}),
                    std::invalid_argument);
    CHECK_THROWS_AS(construct({static_cast<TriggerMode>(99), 0.0}),
                    std::invalid_argument);
    CHECK_NOTHROW(construct({TriggerMode::Position,
                             std::numeric_limits<double>::infinity()}));
    CHECK_NOTHROW(construct({TriggerMode::TimeDelay,
                             std::numeric_limits<double>::infinity()}));
}

TEST_CASE("GPU multi-stage completes zero-voltage crowbar stages and records both summaries") {
    auto value = make_multi_case();
    value.excitations.clear();
    value.triggers = {{TriggerMode::TimeDelay, 1e-6}};
    value.excitations.push_back(std::make_unique<CrowbarExcitation>(0.0, 0.001));
    value.excitations.push_back(std::make_unique<CrowbarExcitation>(0.0, 0.001));
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                       std::move(value.excitations), std::move(value.triggers),
                                       1e-6, false, GpuOptLevel::Standard, backend);

    sim.step();
    sim.step();
    REQUIRE(sim.result().history.size() == 2);
    sim.run(short_policy(0));
    REQUIRE(sim.result().summary.per_stage.size() == 2);
    CHECK(sim.result().summary.per_stage[0].stage_index == 0);
    CHECK(sim.result().summary.per_stage[1].stage_index == 1);
}

TEST_CASE("GPU multi-stage thermal state and resistance reset match the CPU contract") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                       std::move(value.excitations), std::move(value.triggers),
                                       1e-6, true, GpuOptLevel::Standard, backend);

    sim.step();
    REQUIRE(sim.state().filament_temperatures.size() == value.armature.total_filaments());
    CHECK(sim.state().filament_temperatures(0) > T_REFERENCE);
    REQUIRE(sim.filament_resistances().size() == value.armature.resistances().size());
    CHECK(sim.filament_resistances()[0] > value.armature.resistances()[0]);
    CHECK(sim.execution_report().thermal == ThermalMode::Cpu);

    sim.reset();
    CHECK(sim.step_count() == 0);
    CHECK(sim.result().history.empty());
    CHECK(sim.state().filament_temperatures(0) == doctest::Approx(T_REFERENCE));
    CHECK(sim.filament_resistances()[0] == doctest::Approx(value.armature.resistances()[0]));
}

TEST_CASE("GPU multi-stage thermal and triggered state match the CPU reference") {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, true,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, true, GpuOptLevel::Standard, backend);

    const double previous_gpu_velocity = gpu.state().arm_velocity;
    const auto& cpu_step = cpu.step();
    const auto& gpu_step = gpu.step();
    REQUIRE(cpu.state().currents.size() == gpu.state().currents.size());
    const auto tolerance = gpu_test::tolerance_for(GpuOptLevel::Standard);
    for (Eigen::Index i = 0; i < cpu.state().currents.size(); ++i)
        CHECK(gpu_test::numerically_equal(gpu.state().currents(i), cpu.state().currents(i), tolerance));
    CHECK(numerically_equal(gpu.state().arm_position, cpu.state().arm_position, tolerance));
    CHECK(numerically_equal(gpu.state().arm_velocity, cpu.state().arm_velocity, tolerance));
    for (Eigen::Index i = 0; i < cpu.state().filament_temperatures.size(); ++i)
        CHECK(numerically_equal(gpu.state().filament_temperatures(i),
                                cpu.state().filament_temperatures(i), tolerance));
    const double applied_gpu_force = gpu_case.armature.mass() *
        (gpu_step.state.arm_velocity - previous_gpu_velocity) / 1e-6;
    CHECK(numerically_equal(gpu_step.state.force, cpu_step.state.force, tolerance));
    REQUIRE(gpu.result().history.front().coil_currents.size() == 2);
    CHECK(gpu.result().history.front().coil_currents[1] == doctest::Approx(0.0));
}

TEST_CASE("GPU multi-stage fallback matches CPU force, position, and velocity on an asymmetric case") {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    cpu_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    gpu_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    const double initial_position = cpu_case.armature.position();
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    const auto tolerance = tolerance_for(GpuOptLevel::Standard);
    double final_cpu_force = 0.0;
    for (int i = 0; i < 6; ++i) {
        const double previous_cpu_velocity = cpu.state().arm_velocity;
        const auto& cpu_step = cpu.step();
        const double previous_gpu_velocity = gpu.state().arm_velocity;
        const auto& gpu_step = gpu.step();
        const double cpu_force = cpu_case.armature.mass() *
            (cpu_step.state.arm_velocity - previous_cpu_velocity) / 1e-6;
        final_cpu_force = cpu_force;
        CHECK(numerically_equal(gpu_step.state.arm_position, cpu_step.state.arm_position,
                                tolerance));
        CHECK(numerically_equal(gpu_step.state.arm_velocity, cpu_step.state.arm_velocity,
                                tolerance));
        const double applied_gpu_force = gpu_case.armature.mass() *
            (gpu_step.state.arm_velocity - previous_gpu_velocity) / 1e-6;
        CHECK(numerically_equal(gpu_step.state.force, cpu_step.state.force, tolerance));
    }
    REQUIRE(std::abs(gpu.result().history.back().state.force) > 1e-18);
    REQUIRE(std::abs(final_cpu_force) > 1e-18);
    CHECK(std::abs(gpu.result().history.back().state.arm_velocity) > 0.0);
    CHECK(gpu.result().history.back().state.arm_position != doctest::Approx(initial_position));
    REQUIRE(gpu.result().history.back().coil_currents[0] != 0.0);
    REQUIRE(gpu.result().history.back().coil_currents[1] != 0.0);
}

TEST_CASE("GPU multi-stage summary captures capacitor energy before excitation depletion") {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    cpu_case.triggers = {{TriggerMode::TimeDelay, 2e-6}};
    gpu_case.triggers = {{TriggerMode::TimeDelay, 2e-6}};
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    cpu.run(short_policy(8));
    gpu.run(short_policy(8));
    REQUIRE(cpu.result().summary.per_stage.size() == 2);
    REQUIRE(gpu.result().summary.per_stage.size() == 2);
    CHECK(cpu.result().summary.per_stage[1].trigger_time == doctest::Approx(2e-6));
    CHECK(gpu.result().summary.per_stage[1].trigger_time == doctest::Approx(2e-6));
    const auto tolerance = tolerance_for(GpuOptLevel::Standard);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        CHECK(cpu.result().summary.per_stage[stage].energy_depleted > 0.0);
        CHECK(gpu.result().summary.per_stage[stage].energy_depleted > 0.0);
        CHECK(numerically_equal(gpu.result().summary.per_stage[stage].energy_depleted,
                                cpu.result().summary.per_stage[stage].energy_depleted,
                                tolerance));
    }
}

TEST_CASE("GPU multi-stage exposes fallback M and dM snapshots without raw pointers") {
    auto value = make_multi_case();
    const auto coils = value.coils;
    const auto armature = value.armature;
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    sim.step();
    const auto mutual = sim.mutual_inductances();
    const auto gradients = sim.mutual_gradients();
    REQUIRE(mutual.size() == 2 * value.armature.total_filaments());
    REQUIRE(gradients.size() == mutual.size());
    const auto& coil = coils[0];
    const double separation = coil.position() - armature.filament_axial_position(1);
    const double expected_mutual = mutual_inductance_coil(
        coil.inner_radius(), coil.outer_radius(), coil.length(), coil.turns(),
        armature.filament_inner_radius(1), armature.filament_outer_radius(1),
        armature.length() / armature.axial_filaments(), 1,
        separation, 9, false);
    const double expected_gradient = mutual_inductance_gradient_coil(
        coil.inner_radius(), coil.outer_radius(), coil.length(), coil.turns(),
        armature.filament_inner_radius(1), armature.filament_outer_radius(1),
        armature.length() / armature.axial_filaments(), 1,
        separation, 9, false);
    const auto tolerance = tolerance_for(GpuOptLevel::Standard);
    CHECK(numerically_equal(mutual[0], expected_mutual, tolerance));
    CHECK(numerically_equal(gradients[0], expected_gradient, tolerance));
}

TEST_CASE("GPU multi-stage completion terminates a nonzero bounded run without advancing again") {
    auto value = make_multi_case();
    value.triggers = {{TriggerMode::TimeDelay, 0.0}};
    value.excitations.clear();
    value.excitations.push_back(std::make_unique<CrowbarExcitation>(0.0, 0.001));
    value.excitations.push_back(std::make_unique<CrowbarExcitation>(0.0, 0.0008));
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend);
    auto policy = short_policy(8);

    sim.run(policy);
    REQUIRE(sim.result().summary.per_stage.size() == 2);
    REQUIRE(sim.step_count() > 0);
    REQUIRE(sim.step_count() < policy.max_steps);
    const auto history_size = sim.result().history.size();
    const auto completed_steps = sim.step_count();
    sim.run(policy);
    CHECK(sim.result().history.size() == history_size);
    CHECK(sim.step_count() == completed_steps);
}

TEST_CASE("GPU Full keeps a distant triggered circuit active until mutual coupling enters range" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    auto make_case = [] {
        MultiCase value{
            {
                DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0),
                DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.6)
            },
            Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     100000.0, 4.0, 5, 2, 0.0),
            {},
            {{TriggerMode::TimeDelay, 0.0}}
        };
        value.excitations.push_back(std::make_unique<CrowbarExcitation>(480.0, 0.001));
        value.excitations.push_back(std::make_unique<CrowbarExcitation>(320.0, 0.0008));
        return value;
    };

    auto cpu_case = make_case();
    auto gpu_case = make_case();
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Full);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Full, backend);

    const auto& cpu_first = cpu.step();
    const auto& gpu_first = gpu.step();
    REQUIRE(gpu.execution_report().gpu_executed);
    REQUIRE(gpu.execution_report().backend == BackendMode::Direct);
    REQUIRE(gpu_first.coil_currents[1] != doctest::Approx(0.0));
    CHECK(gpu_first.coil_currents[1] == doctest::Approx(cpu_first.coil_currents[1]).epsilon(1e-4));
    CHECK(gpu_first.state.force == doctest::Approx(0.0));
    CHECK(gpu.state().arm_position == doctest::Approx(0.1));

    const auto& cpu_second = cpu.step();
    const auto& gpu_second = gpu.step();
    CHECK(numerically_equal(gpu_second.coil_currents[1], cpu_second.coil_currents[1],
                            tolerance_for(GpuOptLevel::Full)));
    CHECK(std::abs(gpu_second.state.force) > 1e-12);
}

TEST_CASE("GPU Aggressive keeps a distant triggered circuit active until mutual coupling enters range" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    auto make_case = [] {
        MultiCase value{
            {
                DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.6),
                DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 1.2)
            },
            Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     110000.0, 4.0, 5, 2, 0.0),
            {},
            {{TriggerMode::TimeDelay, 0.0}}
        };
        value.excitations.push_back(std::make_unique<CrowbarExcitation>(480.0, 0.001));
        value.excitations.push_back(std::make_unique<CrowbarExcitation>(320.0, 0.0008));
        return value;
    };

    auto cpu_case = make_case();
    auto gpu_case = make_case();
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Full);
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Aggressive, backend);

    const auto& cpu_first = cpu.step();
    const auto& gpu_first = gpu.step();
    REQUIRE(gpu.execution_report().gpu_executed);
    REQUIRE(gpu.execution_report().backend == BackendMode::Direct);
    const auto tolerance = tolerance_for(GpuOptLevel::Aggressive);
    CHECK(numerically_equal(gpu_first.coil_currents[1], cpu_first.coil_currents[1], tolerance));
    REQUIRE(gpu_first.stage_forces.size() == 2);
    CHECK(std::abs(gpu_first.stage_forces[1]) <= tolerance.absolute);
    const auto first_gradients = gpu.mutual_gradients();
    const auto filament_count = static_cast<std::size_t>(gpu_case.armature.total_filaments());
    REQUIRE(first_gradients.size() == 2 * filament_count);
    REQUIRE(std::abs(gpu.state().arm_position - 1.2) > 0.5);
    for (std::size_t filament = 0; filament < filament_count; ++filament)
        CHECK(std::abs(first_gradients[filament_count + filament]) <= tolerance.absolute);

    bool crossed_stage_one_cutoff = false;
    for (int step = 1; step < 12; ++step) {
        const double position_before_step = gpu.state().arm_position;
        const bool stage_one_inside_before_step =
            std::abs(position_before_step - 1.2) <= 0.5;
        const auto& cpu_step = cpu.step();
        const auto& gpu_step = gpu.step();
        CHECK(numerically_equal(gpu_step.state.arm_position, cpu_step.state.arm_position,
                                tolerance));
        if (!stage_one_inside_before_step) {
            REQUIRE(gpu_step.stage_forces.size() == 2);
            CHECK(std::abs(gpu_step.stage_forces[1]) <= tolerance.absolute);
            CHECK(numerically_equal(gpu_step.stage_forces[1], cpu_step.stage_forces[1],
                                    tolerance));
            const auto gradients = gpu.mutual_gradients();
            REQUIRE(gradients.size() == 2 * filament_count);
            for (std::size_t filament = 0; filament < filament_count; ++filament)
                CHECK(std::abs(gradients[filament_count + filament]) <= tolerance.absolute);
            if (std::abs(gpu.state().arm_position - 1.2) <= 0.5) {
                crossed_stage_one_cutoff = true;
                CHECK(std::abs(position_before_step - 1.2) > 0.5);
            }
        } else {
            REQUIRE(crossed_stage_one_cutoff);
            REQUIRE(gpu_step.stage_forces.size() == 2);
            CHECK(std::abs(gpu_step.stage_forces[1]) > tolerance.absolute);
            const auto gradients = gpu.mutual_gradients();
            bool stage_one_has_gradient = false;
            for (std::size_t filament = 0; filament < filament_count; ++filament) {
                stage_one_has_gradient = stage_one_has_gradient ||
                    std::abs(gradients[filament_count + filament]) >
                    std::numeric_limits<double>::epsilon();
            }
            CHECK(stage_one_has_gradient);
            break;
        }
    }
    REQUIRE(crossed_stage_one_cutoff);
}

TEST_CASE("GPU multi-stage summaries retain distinct signed stage contributions as peak magnitudes") {
    auto value = make_multi_case();
    value.coils[1] = DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.035);
    value.triggers = {{TriggerMode::TimeDelay, 0.0}};
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    sim.run(short_policy(6));
    REQUIRE(sim.result().summary.per_stage.size() == 2);
    const auto& first = sim.result().summary.per_stage[0];
    const auto& second = sim.result().summary.per_stage[1];
    CHECK(first.max_force >= 0.0);
    CHECK(second.max_force >= 0.0);
    CHECK(first.max_force != doctest::Approx(second.max_force));
    double aggregate_peak_magnitude = 0.0;
    bool saw_negative_force = false;
    for (const auto& step : sim.result().history) {
        aggregate_peak_magnitude = std::max(aggregate_peak_magnitude, std::abs(step.state.force));
        saw_negative_force = saw_negative_force || step.state.force < 0.0;
    }
    CHECK(sim.result().summary.max_force == doctest::Approx(aggregate_peak_magnitude));
    CHECK(saw_negative_force);
}

TEST_CASE("GPU multi-stage summary preparation is idempotent") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    const auto& first_result = sim.run(short_policy(8));
    const auto first_summary = first_result.summary;
    const auto first_history_size = first_result.history.size();
    const auto& second_result = sim.run(short_policy(8));
    CHECK(second_result.history.size() == first_history_size);
    CHECK(second_result.summary.per_stage.size() == first_summary.per_stage.size());
    CHECK(second_result.summary.max_force == doctest::Approx(first_summary.max_force));
    CHECK(second_result.summary.peak_coil_current == doctest::Approx(first_summary.peak_coil_current));
    CHECK(second_result.summary.per_stage[0].max_force ==
          doctest::Approx(first_summary.per_stage[0].max_force));
    CHECK(second_result.summary.per_stage[1].max_force ==
          doctest::Approx(first_summary.per_stage[1].max_force));
}

TEST_CASE("GPU multi-stage filament resistance accessor returns an independent value snapshot") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, true, GpuOptLevel::Standard, backend);

    static_assert(std::is_same_v<decltype(sim.filament_resistances()), std::vector<double>>);
    const auto initial = sim.filament_resistances();
    sim.step();
    const auto heated = sim.filament_resistances();
    REQUIRE(initial.size() == heated.size());
    CHECK(heated[0] > initial[0]);
    sim.reset();
    CHECK(sim.filament_resistances()[0] == doctest::Approx(initial[0]));
    CHECK(heated[0] > initial[0]);
}

TEST_CASE("GPU multi-stage records post-step force while velocity uses applied pre-step force") {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    const auto cpu_armature = cpu_case.armature;
    cpu_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    gpu_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    const double previous_cpu_velocity = cpu.state().arm_velocity;
    const double previous_gpu_velocity = gpu.state().arm_velocity;
    const auto& cpu_step = cpu.step();
    const auto& gpu_step = gpu.step();
    const auto tolerance = tolerance_for(GpuOptLevel::Standard);
    REQUIRE(std::abs(cpu_step.state.force) > 1e-12);
    CHECK(numerically_equal(gpu_step.state.force, cpu_step.state.force, tolerance));

    const double applied_cpu_force = cpu_armature.mass() *
        (cpu_step.state.arm_velocity - previous_cpu_velocity) / 1e-6;
    const double applied_gpu_force = cpu_armature.mass() *
        (gpu_step.state.arm_velocity - previous_gpu_velocity) / 1e-6;
    CHECK(numerically_equal(applied_gpu_force, applied_cpu_force, tolerance));
    CHECK(std::abs(gpu_step.state.force - applied_gpu_force) > 1e-12);

    REQUIRE(gpu_step.stage_forces.size() == cpu_step.stage_forces.size());
    for (std::size_t stage = 0; stage < gpu_step.stage_forces.size(); ++stage)
        CHECK(numerically_equal(gpu_step.stage_forces[stage], cpu_step.stage_forces[stage],
                                tolerance));

    const double previous_cpu_velocity_second = cpu.state().arm_velocity;
    const double previous_gpu_velocity_second = gpu.state().arm_velocity;
    const auto& cpu_step_second = cpu.step();
    const auto& gpu_step_second = gpu.step();
    const double applied_cpu_force_second = cpu_armature.mass() *
        (cpu_step_second.state.arm_velocity - previous_cpu_velocity_second) / 1e-6;
    const double applied_gpu_force_second = cpu_armature.mass() *
        (gpu_step_second.state.arm_velocity - previous_gpu_velocity_second) / 1e-6;
    CHECK(std::abs(applied_gpu_force_second) > 1e-12);
    CHECK(applied_gpu_force_second * applied_cpu_force_second > 0.0);
    CHECK(numerically_equal(gpu_step_second.state.force, cpu_step_second.state.force, tolerance));
    REQUIRE(gpu_step_second.stage_forces.size() == cpu_step_second.stage_forces.size());
    for (std::size_t stage = 0; stage < gpu_step_second.stage_forces.size(); ++stage)
        CHECK(numerically_equal(gpu_step_second.stage_forces[stage],
                                cpu_step_second.stage_forces[stage], tolerance));

    gpu.run(short_policy(0));
    cpu.run(short_policy(0));
    REQUIRE(gpu.result().summary.per_stage.size() == 2);
    REQUIRE(cpu.result().summary.per_stage.size() == 2);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        CHECK(numerically_equal(gpu.result().summary.per_stage[stage].max_force,
                                cpu.result().summary.per_stage[stage].max_force, tolerance));
        CHECK(gpu.result().summary.per_stage[stage].max_force >=
              std::abs(gpu_step_second.stage_forces[stage]));
    }
}

TEST_CASE("GPU multi-stage commits completion-boundary force like the CPU reference") {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    cpu_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    gpu_case.triggers = {{TriggerMode::TimeDelay, 0.0}};
    auto cpu_waveform = std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 480.0; });
    auto gpu_waveform = std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 480.0; });
    cpu_waveform->set_end_time(1e-6);
    gpu_waveform->set_end_time(1e-6);
    cpu_case.excitations[0] = std::move(cpu_waveform);
    gpu_case.excitations[0] = std::move(gpu_waveform);

    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    const auto& cpu_step = cpu.step();
    const auto& gpu_step = gpu.step();
    REQUIRE(std::abs(cpu_step.state.force) > 1e-18);
    REQUIRE(cpu_step.stage_forces.size() == 2);
    REQUIRE(gpu_step.stage_forces.size() == 2);
    CHECK(cpu_step.stage_forces[0] == doctest::Approx(0.0));
    CHECK(gpu_step.stage_forces[0] == doctest::Approx(cpu_step.stage_forces[0]));
    CHECK(numerically_equal(gpu_step.state.force, cpu_step.state.force,
                            tolerance_for(GpuOptLevel::Standard)));
    CHECK(gpu_step.stage_forces[1] == doctest::Approx(cpu_step.stage_forces[1]));
}

TEST_CASE("GPU velocity-decay termination uses the committed completion-boundary force" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    auto make_case = [] {
        MultiCase value{
            {
                DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0),
                DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 1.0)
            },
            Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                     4.0, 0.120, 5, 2, -0.08),
            {},
            {{TriggerMode::TimeDelay, 0.0}}
        };
        auto completing = std::make_unique<coilgun::simulation::WaveformExcitation>(
            [](double) { return 480.0; });
        completing->set_end_time(2e-6);
        value.excitations.push_back(std::move(completing));
        value.excitations.push_back(
            std::make_unique<coilgun::simulation::WaveformExcitation>(
                [](double) { return 1.0; }));
        return value;
    };

    auto cpu_case = make_case();
    auto gpu_case = make_case();
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Full);
    GpuBackend backend;
    backend.backend = BackendMode::Auto;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Full, backend);

    cpu.step();
    gpu.step();
    const double velocity_before_boundary = gpu.state().arm_velocity;
    const auto& cpu_boundary = cpu.step();
    const auto& gpu_boundary = gpu.step();
    REQUIRE(gpu.execution_report().gpu_executed);
    REQUIRE(gpu.execution_report().backend == BackendMode::Direct);
    CAPTURE(velocity_before_boundary);
    CAPTURE(gpu_boundary.state.arm_velocity);
    CAPTURE(gpu.result().history.front().state.force);
    CAPTURE(gpu_boundary.state.force);
    REQUIRE(gpu_boundary.state.arm_velocity < velocity_before_boundary);
    REQUIRE(std::abs(gpu_boundary.state.force) < 1e-18);
    const double applied_acceleration =
        std::abs(gpu_boundary.state.arm_velocity - velocity_before_boundary) / 1e-6;
    REQUIRE(applied_acceleration > 1e-6);
    CHECK(numerically_equal(gpu_boundary.state.force, cpu_boundary.state.force,
                            tolerance_for(GpuOptLevel::Standard)));

    coilgun::simulation::TerminationPolicy policy;
    policy.max_steps = 3;
    policy.velocity_decay_steps = 1;
    policy.accel_threshold = 0.5 * applied_acceleration;
    policy.enable_velocity_check = true;
    policy.enable_bound_check = false;
    cpu.run(policy);
    gpu.run(policy);

    CHECK(cpu.step_count() == 2);
    CHECK(gpu.step_count() == cpu.step_count());
    REQUIRE(gpu.result().history.size() == cpu.result().history.size());
    for (std::size_t step = 0; step < cpu.result().history.size(); ++step) {
        CHECK(gpu.result().history[step].state.time ==
              doctest::Approx(cpu.result().history[step].state.time));
        CHECK(numerically_equal(gpu.result().history[step].state.arm_velocity,
                                cpu.result().history[step].state.arm_velocity,
                                tolerance_for(GpuOptLevel::Standard)));
    }
}

TEST_CASE("GPU multi-stage trigger position matches the CPU crossed-boundary position") {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    cpu_case.triggers = {{TriggerMode::Position, 8e-6}};
    gpu_case.triggers = {{TriggerMode::Position, 8e-6}};
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, false,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, false, GpuOptLevel::Standard, backend);

    bool cpu_triggered = false;
    bool gpu_triggered = false;
    for (int i = 0; i < 100 && (!cpu_triggered || !gpu_triggered); ++i) {
        cpu.step();
        gpu.step();
        cpu_triggered = cpu.result().history.back().coil_currents[1] != 0.0;
        gpu_triggered = gpu.result().history.back().coil_currents[1] != 0.0;
    }
    REQUIRE(cpu_triggered);
    REQUIRE(gpu_triggered);

    cpu.run(short_policy(0));
    gpu.run(short_policy(0));
    REQUIRE(cpu.result().summary.per_stage.size() == 2);
    REQUIRE(gpu.result().summary.per_stage.size() == 2);
    CHECK(numerically_equal(gpu.result().summary.per_stage[1].trigger_position,
                            cpu.result().summary.per_stage[1].trigger_position,
                            tolerance_for(GpuOptLevel::Standard)));
    CHECK(gpu.result().summary.per_stage[1].trigger_position <
          gpu.result().history.back().state.arm_position);
}

TEST_CASE("GPU multi-stage primary comparison requires CUDA execution" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    auto cpu_case = make_multi_case();
    auto gpu_case = make_multi_case();
    MultiStageSim<EulerStepper> cpu(
        std::move(cpu_case.coils), cpu_case.armature, std::move(cpu_case.excitations),
        std::move(cpu_case.triggers), 1e-6, true,
        coilgun::simulation::OptimizationLevel::Reference);
    GpuBackend backend;
    backend.backend = BackendMode::Auto;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> gpu(
        std::move(gpu_case.coils), gpu_case.armature, std::move(gpu_case.excitations),
        std::move(gpu_case.triggers), 1e-6, true, GpuOptLevel::Standard, backend);

    for (int i = 0; i < 4; ++i) {
        const double previous_gpu_velocity = gpu.state().arm_velocity;
        const auto& cpu_step = cpu.step();
        const auto& gpu_step = gpu.step();
        REQUIRE(gpu.execution_report().gpu_executed);
        REQUIRE(gpu.execution_report().backend == BackendMode::Direct);
        REQUIRE(cpu.state().currents.size() == gpu.state().currents.size());
        for (Eigen::Index j = 0; j < cpu.state().currents.size(); ++j)
            CHECK(numerically_equal(gpu.state().currents(j), cpu.state().currents(j), tolerance_for(GpuOptLevel::Standard)));
        CHECK(numerically_equal(gpu.state().arm_position, cpu.state().arm_position,
                                tolerance_for(GpuOptLevel::Standard)));
        CHECK(numerically_equal(gpu.state().arm_velocity, cpu.state().arm_velocity,
                                tolerance_for(GpuOptLevel::Standard)));
        const double cpu_force = cpu_case.armature.mass() *
            (cpu_step.state.arm_velocity - (i == 0 ? 0.0 :
                cpu.result().history[static_cast<std::size_t>(i - 1)].state.arm_velocity)) / 1e-6;
        CHECK(numerically_equal(gpu_step.state.force, cpu_step.state.force,
                                tolerance_for(GpuOptLevel::Standard)));
        REQUIRE(gpu_step.cap_voltages.size() == cpu_step.cap_voltages.size());
        REQUIRE(gpu_step.coil_currents.size() == cpu_step.coil_currents.size());
        for (std::size_t stage = 0; stage < gpu_step.coil_currents.size(); ++stage) {
            CHECK(gpu_step.cap_voltages[stage] == doctest::Approx(cpu_step.cap_voltages[stage]));
            CHECK(numerically_equal(gpu_step.coil_currents[stage], cpu_step.coil_currents[stage],
                                    tolerance_for(GpuOptLevel::Standard)));
        }
        if (i == 3) {
            CHECK(std::abs(gpu_step.state.force) > 1e-12);
            CHECK(std::abs(gpu_step.state.arm_velocity) > 0.0);
            CHECK(gpu_step.state.arm_position > 0.0);
            CHECK(std::abs(gpu_step.coil_currents[0]) > 0.0);
            CHECK(std::abs(gpu_step.coil_currents[1]) > 0.0);
        }
        REQUIRE(gpu.filament_resistances().size() == cpu_case.armature.resistances().size());
        const double beta = material_beta(cpu_case.armature.material());
        for (Eigen::Index filament = 0; filament < cpu.state().filament_temperatures.size(); ++filament) {
            const double expected_resistance = cpu_case.armature.resistances()[static_cast<std::size_t>(filament)] *
                (1.0 + beta * (cpu.state().filament_temperatures(filament) - T_REFERENCE));
            CHECK(numerically_equal(
                gpu.filament_resistances()[static_cast<std::size_t>(filament)],
                expected_resistance, tolerance_for(GpuOptLevel::Standard)));
        }
    }
}

TEST_CASE("GPU multi-stage explicit fallback is a separate CPU-only execution path") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Auto;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend,
        BackendMode::Fallback);
    sim.step();
    CHECK(sim.execution_report().requested_backend == BackendMode::Fallback);
    CHECK_FALSE(sim.execution_report().gpu_executed);
    CHECK(sim.execution_report().backend == BackendMode::Fallback);
}

TEST_CASE("GPU multi-stage Graph executes on CUDA and matches explicit fallback" *
          doctest::skip(!coilgun::simulation::cuda::cuda_device_available())) {
    auto make_sim = [](BackendMode mode, bool persistent, BackendMode explicit_mode) {
        auto value = make_multi_case();
        GpuBackend backend;
        backend.backend = mode;
        backend.use_persistent = persistent;
        return std::make_unique<GpuMultiStageSim<EulerStepper>>(
            std::move(value.coils), value.armature, std::move(value.excitations),
            std::move(value.triggers), 1e-6, false, GpuOptLevel::Standard, backend,
            explicit_mode);
    };

    auto fallback = make_sim(BackendMode::Fallback, false, BackendMode::Fallback);
    auto graph = make_sim(BackendMode::Graph, false, BackendMode::Graph);
    fallback->run(short_policy());
    graph->run(short_policy());

    const auto tolerance = tolerance_for(GpuOptLevel::Standard);
    CHECK(numerically_equal(graph->result().summary.muzzle_velocity,
                            fallback->result().summary.muzzle_velocity, tolerance));
    CHECK(graph->execution_report().requested_backend == BackendMode::Graph);
    CHECK(fallback->execution_report().requested_backend == BackendMode::Fallback);
    REQUIRE(graph->execution_report().gpu_executed);
    CHECK(graph->execution_report().backend == BackendMode::Graph);
    CHECK(graph->execution_report().graph_rebuild_count > 0);
    CHECK(graph->graph_assisted());
    CHECK_FALSE(fallback->execution_report().gpu_executed);
    CHECK(fallback->execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(fallback->graph_assisted());
}

TEST_CASE("GPU multi-stage reset deterministically replays the first step") {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Fallback;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                       std::move(value.excitations), std::move(value.triggers),
                                       1e-6, true, GpuOptLevel::Standard, backend);
    sim.step();
    const auto first = sim.result().history.front();
    sim.step();
    const auto report = sim.execution_report();
    sim.reset();
    sim.step();

    REQUIRE(sim.result().history.size() == 1);
    CHECK(sim.result().history.front().state.time == doctest::Approx(first.state.time));
    CHECK(sim.result().history.front().state.force == doctest::Approx(first.state.force));
    CHECK(sim.state().currents(0) ==
          doctest::Approx(sim.result().history.front().coil_currents.front()));
    CHECK(sim.execution_report().fallback_count == report.fallback_count);
}

TEST_CASE("GPU multi-stage reports safe fallback when no CUDA device is available" *
          doctest::skip(coilgun::simulation::cuda::cuda_device_available())) {
    auto value = make_multi_case();
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = false;
    GpuMultiStageSim<EulerStepper> sim(std::move(value.coils), value.armature,
                                       std::move(value.excitations), std::move(value.triggers),
                                       1e-6, false, GpuOptLevel::Standard, backend);
    sim.step();
    CHECK_FALSE(sim.execution_report().gpu_executed);
    CHECK(sim.execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(sim.execution_report().fallback_reason.empty());
}

TEST_CASE("GPU multi-stage rejects RK4 rather than silently changing integration") {
    auto value = make_multi_case();
    GpuMultiStageSim<coilgun::simulation::RK4Stepper> sim(
        std::move(value.coils), value.armature, std::move(value.excitations),
        std::move(value.triggers), 1e-6);
    CHECK_THROWS_WITH(sim.step(), "RK4Stepper is not supported by GpuMultiStageSim");
}

TEST_CASE("GPU multi-stage validates migrated constructor inputs") {
    auto value = make_multi_case();
    CHECK_THROWS_AS(
        GpuMultiStageSim<EulerStepper>(
            {}, value.armature, std::move(value.excitations), {}, 1e-6),
        std::invalid_argument);

    auto null_case = make_multi_case();
    null_case.excitations[0].reset();
    CHECK_THROWS_AS(
        GpuMultiStageSim<EulerStepper>(
            std::move(null_case.coils), null_case.armature,
            std::move(null_case.excitations), std::move(null_case.triggers), 1e-6),
        std::invalid_argument);

    auto dt_case = make_multi_case();
    CHECK_THROWS_AS(
        GpuMultiStageSim<EulerStepper>(
            std::move(dt_case.coils), dt_case.armature,
            std::move(dt_case.excitations), std::move(dt_case.triggers),
            std::numeric_limits<double>::quiet_NaN()),
        std::invalid_argument);

    auto mode_case = make_multi_case();
    CHECK_THROWS_AS(
        GpuMultiStageSim<EulerStepper>(
            std::move(mode_case.coils), mode_case.armature,
            std::move(mode_case.excitations), std::move(mode_case.triggers),
            1e-6, false, GpuOptLevel::Standard, GpuBackend{},
            static_cast<BackendMode>(99)),
        std::invalid_argument);

    auto bad_opt_multi = make_multi_case();
    CHECK_THROWS_AS(
        GpuMultiStageSim<EulerStepper>(
            std::move(bad_opt_multi.coils), bad_opt_multi.armature,
            std::move(bad_opt_multi.excitations), std::move(bad_opt_multi.triggers),
            1e-6, false, static_cast<GpuOptLevel>(99)),
        std::invalid_argument);
}
