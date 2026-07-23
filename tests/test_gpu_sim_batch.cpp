/**
 * @file test_gpu_sim_batch.cpp
 * @brief Integration test: SimBatch — GPU-accelerated batch parameter sweep.
 * @author Winston Meursault
 *
 * Verifies that SimBatch produces the same results as running the
 * same parameter sweep via individual GpuMultiStageSim objects.
 */

#include <doctest/doctest.h>
#include "coilgun/simulation/cuda/sim_batch.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/components/driving_coil.hpp"
#include "coilgun/components/armature.hpp"
#include "coilgun/simulation/excitation.hpp"
#include "coilgun/physics/constants.hpp"
#include "coilgun/physics/mutual_inductance.hpp"
#include <memory>
#include <limits>
#include <cmath>
#include <vector>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::RK4Stepper;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::cuda::SimBatch;
using coilgun::simulation::cuda::GpuMultiStageSim;
using coilgun::simulation::cuda::GpuBackend;
using coilgun::simulation::cuda::BackendMode;
using coilgun::simulation::cuda::SolverMode;
using coilgun::simulation::cuda::GpuCapability;
using coilgun::simulation::cuda::GpuExecutionConfig;
using coilgun::simulation::cuda::GpuExecutionPlanner;

struct SweepPoint { double voltage; };

static auto make_coils() {
    DrivingCoil c1(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    return std::vector<DrivingCoil>{c1, c2};
}

static auto make_arm() {
    return Armature(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref,
                    ALUMINUM.density, 0.0, 0.120, 5, 2, 0.05);
}

static GpuBackend direct_backend() {
    GpuBackend backend;
    backend.use_persistent = false;
    backend.backend = BackendMode::Direct;
    return backend;
}

TEST_CASE("SimBatch — small sweep vs individual GPU runs") {
    std::vector<SweepPoint> sweeps = {{350.0}, {450.0}};
    REQUIRE(sweeps.size() == 2);

    std::vector<double> v_ref(2);
    for (int i = 0; i < 2; ++i) {
        auto coils = make_coils();
        auto arm   = make_arm();
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        GpuMultiStageSim<EulerStepper> sim(coils, arm, std::move(excs), triggers, 1e-6);
        auto policy = coilgun::simulation::TerminationPolicy::defaults();
        policy.max_steps = 500;
        sim.run(policy);
        v_ref[i] = sim.result().summary.muzzle_velocity;
    }

    auto coils = make_coils();
    auto arm   = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 2, 1e-6, direct_backend());

    for (int i = 0; i < 2; ++i) {
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(sweeps[i].voltage, 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        batch.set_excitations(i, std::move(excs), triggers);
    }

    auto pol = coilgun::simulation::TerminationPolicy::defaults();
    pol.max_steps = 500;
    batch.run(pol);

    for (int i = 0; i < 2; ++i) {
        double v_batch = batch.result(i).summary.muzzle_velocity;
        INFO("Sweep " << i << " voltage=" << sweeps[i].voltage);
        CHECK(v_batch == doctest::Approx(v_ref[i]).epsilon(1e-5));
    }
}

TEST_CASE("SimBatch — representative history and summaries match individual GPU runs") {
    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 8;

    auto individual_excitations = std::vector<std::unique_ptr<coilgun::simulation::Excitation>>{};
    individual_excitations.push_back(std::make_unique<CrowbarExcitation>(425.0, 0.001));
    individual_excitations.push_back(std::make_unique<CrowbarExcitation>(325.0, 0.001));
    GpuMultiStageSim<EulerStepper> individual(
        make_coils(), make_arm(), std::move(individual_excitations),
        {{TriggerMode::Position, 0.09}}, 1e-6);
    individual.run(policy);

    SimBatch<EulerStepper> batch(make_coils(), make_arm(), 1, 1e-6, direct_backend());
    auto batch_excitations = std::vector<std::unique_ptr<coilgun::simulation::Excitation>>{};
    batch_excitations.push_back(std::make_unique<CrowbarExcitation>(425.0, 0.001));
    batch_excitations.push_back(std::make_unique<CrowbarExcitation>(325.0, 0.001));
    batch.set_excitations(0, std::move(batch_excitations), {{TriggerMode::Position, 0.09}});
    batch.run(policy);

    const auto& expected = individual.result();
    const auto& actual = batch.result(0);
    REQUIRE(actual.history.size() == expected.history.size());
    REQUIRE(actual.summary.per_stage.size() == expected.summary.per_stage.size());
    CHECK(actual.summary.step_count == expected.summary.step_count);
    CHECK(actual.summary.muzzle_velocity == doctest::Approx(expected.summary.muzzle_velocity).epsilon(1e-5));
    for (std::size_t step = 0; step < expected.history.size(); ++step) {
        CHECK(actual.history[step].state.arm_position ==
              doctest::Approx(expected.history[step].state.arm_position).epsilon(1e-5));
        CHECK(actual.history[step].state.arm_velocity ==
              doctest::Approx(expected.history[step].state.arm_velocity).epsilon(1e-5));
        CHECK(actual.history[step].state.force ==
              doctest::Approx(expected.history[step].state.force).epsilon(1e-5));
        REQUIRE(actual.history[step].stage_forces.size() == expected.history[step].stage_forces.size());
        for (std::size_t stage = 0; stage < expected.history[step].stage_forces.size(); ++stage)
            CHECK(actual.history[step].stage_forces[stage] ==
                  doctest::Approx(expected.history[step].stage_forces[stage]).epsilon(1e-5));
    }
}

TEST_CASE("SimBatch — recorded force uses post-step currents and pre-position gradient") {
    auto coils = make_coils();
    auto arm = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6, direct_backend());
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excitations;
    excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    batch.set_excitations(0, std::move(excitations), {{TriggerMode::Position, 0.09}});

    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 1;
    batch.run(policy);

    const auto& history = batch.result(0).history;
    REQUIRE(history.size() == 1);
    const auto& entry = history.front();
    const auto gradients = batch.mutual_gradients();
    const auto armature = make_arm();
    const auto coils_for_expected = make_coils();
    double expected_force = 0.0;
    const int radial = armature.radial_filaments();
    for (int filament = 0; filament < armature.total_filaments(); ++filament) {
        const int radial_id = filament % radial + 1;
        const int axial_id = filament / radial + 1;
        const double separation = coils_for_expected[0].position() -
            armature.filament_axial_position(axial_id);
        const double gradient = coilgun::physics::mutual_inductance_gradient_coil(
            coils_for_expected[0].inner_radius(), coils_for_expected[0].outer_radius(),
            coils_for_expected[0].length(), coils_for_expected[0].turns(),
            armature.filament_inner_radius(radial_id), armature.filament_outer_radius(radial_id),
            armature.length() / armature.axial_filaments(), 1, separation, 9, true);
        CHECK(gradients[static_cast<std::size_t>(filament)] ==
              doctest::Approx(gradient).epsilon(1e-10));
        expected_force -= entry.coil_currents[0] * entry.state.filament_currents[filament] * gradient;
    }
    CHECK(entry.stage_forces[0] == doctest::Approx(expected_force).epsilon(1e-8));
}

TEST_CASE("SimBatch — direction consistent") {
    auto coils = make_coils();
    auto arm   = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6, direct_backend());

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
    batch.set_excitations(0, std::move(excs), triggers);

    auto pol = coilgun::simulation::TerminationPolicy::defaults();
    pol.max_steps = 500;
    batch.run(pol);
    double v = batch.result(0).summary.muzzle_velocity;
    CHECK(v > 0.0);
}

TEST_CASE("SimBatch — persistent request reports safe fallback") {
    auto coils = make_coils();
    auto arm = make_arm();
    GpuBackend backend;
    backend.use_persistent = true;
    backend.backend = BackendMode::Persistent;
    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6, backend);
    CHECK(batch.execution_report().requested_backend == BackendMode::Persistent);
    CHECK(batch.execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(batch.execution_report().gpu_executed);
    const auto& report = batch.execution_report();
    CHECK_FALSE(report.fallback_reason.empty());

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
    batch.set_excitations(0, std::move(excs), {{TriggerMode::TimeDelay, 1e-6}});
    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 2;
    batch.run(policy);
    CHECK(batch.result(0).summary.step_count == 2);

    GpuBackend fallback_backend;
    fallback_backend.backend = BackendMode::Fallback;
    SimBatch<EulerStepper> fallback_batch(make_coils(), make_arm(), 1, 1e-6,
                                          fallback_backend);
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> fallback_excs;
    fallback_excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
    fallback_excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
    fallback_batch.set_excitations(0, std::move(fallback_excs),
                                   {{TriggerMode::TimeDelay, 1e-6}});
    fallback_batch.run(policy);
    CHECK(batch.result(0).summary.muzzle_velocity ==
          doctest::Approx(fallback_batch.result(0).summary.muzzle_velocity).epsilon(1e-8));
}

TEST_CASE("SimBatch — persistent batch request keeps metadata and uses bounded fallback") {
    auto coils = make_coils();
    auto arm = make_arm();
    GpuBackend backend;
    backend.use_persistent = true;
    backend.backend = BackendMode::Persistent;
    SimBatch<EulerStepper> batch(std::move(coils), arm, 2, 1e-6, backend);

    CHECK(batch.execution_report().requested_backend == BackendMode::Persistent);
    CHECK(batch.execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(batch.execution_report().gpu_executed);
    const auto& report = batch.execution_report();
    CHECK_FALSE(report.fallback_reason.empty());

    for (int sim_id = 0; sim_id < 2; ++sim_id) {
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
        batch.set_excitations(sim_id, std::move(excs), {{TriggerMode::TimeDelay, 1e-6}});
    }

    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 2;
    batch.run(policy);
    CHECK(batch.result(0).summary.step_count == 2);
    CHECK(batch.result(1).summary.step_count == 2);
}

TEST_CASE("SimBatch — explicit backend metadata takes precedence over legacy persistence") {
    for (const auto mode : {BackendMode::Direct, BackendMode::Graph,
                            BackendMode::Fallback, BackendMode::Persistent}) {
        GpuBackend false_flag;
        false_flag.backend = mode;
        false_flag.use_persistent = false;
        SimBatch<EulerStepper> without_legacy(make_coils(), make_arm(), 1, 1e-6, false_flag);

        GpuBackend true_flag = false_flag;
        true_flag.use_persistent = true;
        SimBatch<EulerStepper> with_legacy(make_coils(), make_arm(), 1, 1e-6, true_flag);

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
    }
}

TEST_CASE("SimBatch — explicit constructor backend overrides conflicting backend fields") {
    GpuBackend backend;
    backend.backend = BackendMode::Direct;
    backend.use_persistent = true;
    SimBatch<EulerStepper> batch(make_coils(), make_arm(), 1, 1e-6, backend,
                                 BackendMode::Graph);
    CHECK(batch.execution_report().requested_backend == BackendMode::Graph);
    CHECK(batch.execution_report().requested_solver == SolverMode::Auto);
}

TEST_CASE("SimBatch — explicit backend paths report bounded equivalent results") {
    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 3;

    struct Run {
        BackendMode requested;
        coilgun::simulation::MultiStageResult result;
        coilgun::simulation::cuda::ExecutionReport report;
    };
    std::vector<Run> runs;
    for (const auto mode : {BackendMode::Fallback, BackendMode::Direct,
                            BackendMode::Graph, BackendMode::Persistent}) {
        GpuBackend backend;
        backend.backend = mode;
        backend.use_persistent = false;
        SimBatch<EulerStepper> batch(make_coils(), make_arm(), 1, 1e-6, backend);
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(425.0, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(325.0, 0.001));
        batch.set_excitations(0, std::move(excs), {{TriggerMode::Position, 0.09}});
        batch.run(policy);

        const auto report = batch.execution_report();
        CHECK(report.requested_backend == mode);
        if (report.gpu_executed) {
            CHECK(report.backend == mode);
        } else {
            CHECK(report.backend == BackendMode::Fallback);
            CHECK_FALSE(report.fallback_reason.empty());
            if (mode == BackendMode::Persistent) {
                CHECK(report.requested_backend == BackendMode::Persistent);
                CHECK_FALSE(report.gpu_executed);
            }
        }
        runs.push_back(Run{mode, batch.result(0), report});
    }

    REQUIRE(runs.size() == 4);
    REQUIRE(runs.front().result.history.size() == runs[1].result.history.size());
    for (std::size_t index = 1; index < runs.size(); ++index) {
        REQUIRE(runs.front().result.history.size() == runs[index].result.history.size());
        for (std::size_t step = 0; step < runs.front().result.history.size(); ++step) {
            CHECK(runs.front().result.history[step].state.arm_position ==
                  doctest::Approx(runs[index].result.history[step].state.arm_position).epsilon(1e-4));
            CHECK(runs.front().result.history[step].state.arm_velocity ==
                  doctest::Approx(runs[index].result.history[step].state.arm_velocity).epsilon(1e-4));
            CHECK(runs.front().result.history[step].state.force ==
                  doctest::Approx(runs[index].result.history[step].state.force).epsilon(1e-4));
        }
    }
    for (std::size_t index = 1; index < runs.size(); ++index) {
        CHECK(runs[index].result.summary.step_count == runs.front().result.summary.step_count);
        CHECK(runs[index].result.summary.muzzle_velocity ==
              doctest::Approx(runs.front().result.summary.muzzle_velocity).epsilon(1e-4));
    }
}

TEST_CASE("SimBatch — Auto uses the legacy persistence request explicitly") {
    auto coils = make_coils();
    auto arm = make_arm();
    GpuBackend backend;
    backend.backend = BackendMode::Auto;
    backend.use_persistent = true;
    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6, backend);
    CHECK(batch.execution_report().requested_backend == BackendMode::Persistent);
    CHECK(batch.execution_report().backend == BackendMode::Fallback);
    CHECK_FALSE(batch.execution_report().gpu_executed);
    const auto& report = batch.execution_report();
    CHECK_FALSE(report.fallback_reason.empty());
}

TEST_CASE("SimBatch — batched solver execution is proven or honestly reported as Eigen fallback") {
    auto coils = make_coils();
    auto arm = make_arm();
    GpuBackend direct;
    direct.backend = BackendMode::Direct;
    direct.use_persistent = true;
    SimBatch<EulerStepper> batch(std::move(coils), arm, 8, 1e-6, direct);
    CHECK(batch.execution_report().requested_solver == coilgun::simulation::cuda::SolverMode::Auto);
    GpuExecutionConfig config;
    config.backend = BackendMode::Direct;
    config.solver = SolverMode::Auto;
    const auto planned = GpuExecutionPlanner::plan(2, 10, 8, false, GpuCapability{}, config);
    CHECK(planned.solver == SolverMode::Batched);
    for (int sim_id = 0; sim_id < 8; ++sim_id) {
        std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(350.0 + sim_id, 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(300.0 + sim_id, 0.001));
        batch.set_excitations(sim_id, std::move(excs), {{TriggerMode::TimeDelay, 1.0}});
    }
    auto solver_policy = coilgun::simulation::TerminationPolicy::defaults();
    solver_policy.max_steps = 1;
    batch.run(solver_policy);
    const auto& report = batch.execution_report();
    if (report.gpu_executed) {
        REQUIRE(report.backend == BackendMode::Direct);
        REQUIRE(report.solver == SolverMode::Batched);
    } else {
        REQUIRE(report.backend == BackendMode::Fallback);
        REQUIRE(report.solver == SolverMode::Eigen);
        REQUIRE_FALSE(report.fallback_reason.empty());
    }

    auto fallback_coils = make_coils();
    auto fallback_arm = make_arm();
    GpuBackend fallback;
    fallback.backend = BackendMode::Fallback;
    SimBatch<EulerStepper> cpu_batch(std::move(fallback_coils), fallback_arm, 8, 1e-6, fallback);
    CHECK(cpu_batch.execution_report().solver == coilgun::simulation::cuda::SolverMode::Eigen);
}

TEST_CASE("SimBatch — heterogeneous rows finish independently without compaction") {
    auto coils = make_coils();
    auto arm = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 2, 1e-6, direct_backend());

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> short_excitations;
    short_excitations.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    short_excitations.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    batch.set_excitations(0, std::move(short_excitations),
                          {{TriggerMode::TimeDelay, std::numeric_limits<double>::infinity()}});

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> long_excitations;
    long_excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    long_excitations.push_back(std::make_unique<CrowbarExcitation>(300.0, 0.001));
    batch.set_excitations(1, std::move(long_excitations),
                          {{TriggerMode::TimeDelay, 1.0}});

    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 4;
    batch.run(policy);

    CHECK(batch.result(0).summary.step_count == 1);
    CHECK(batch.result(1).summary.step_count == 4);
    CHECK(batch.result(0).history.size() == 1);
    CHECK(batch.result(1).history.size() == 4);
    CHECK(batch.result(1).history.front().coil_currents[1] == doctest::Approx(0.0));
    CHECK(batch.result(1).history.back().coil_currents[1] == doctest::Approx(0.0));
    CHECK_FALSE(batch.graph_assisted());
    CHECK(batch.execution_report().requested_backend == BackendMode::Direct);
    const auto resolved_backend = batch.execution_report().backend;
    const bool direct_or_fallback = resolved_backend == BackendMode::Direct ||
        resolved_backend == BackendMode::Fallback;
    CHECK(direct_or_fallback);
}

TEST_CASE("SimBatch — E3 fixed-capacity active_mask rows never compact or reindex") {
    SimBatch<EulerStepper> batch(make_coils(), make_arm(), 2, 1e-6, direct_backend());
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> early;
    early.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    early.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    batch.set_excitations(0, std::move(early), {{TriggerMode::TimeDelay, INFINITY}});

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> continuing;
    continuing.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    continuing.push_back(std::make_unique<CrowbarExcitation>(300.0, 0.001));
    batch.set_excitations(1, std::move(continuing), {{TriggerMode::TimeDelay, 1.0}});

    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 4;
    batch.run(policy);

    const auto& early_result = batch.result(0);
    const auto& continuing_result = batch.result(1);
    // E3 intentionally defers active-list compaction: row 1 remains result(1)
    // after row 0 becomes inactive, rather than being moved to a new index.
    REQUIRE(early_result.summary.step_count == 1);
    REQUIRE(continuing_result.summary.step_count == 4);
    REQUIRE(early_result.history.size() == 1);
    REQUIRE(continuing_result.history.size() == 4);

    SimBatch<EulerStepper> standalone(make_coils(), make_arm(), 1, 1e-6, direct_backend());
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> standalone_excitations;
    standalone_excitations.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    standalone_excitations.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    standalone.set_excitations(0, std::move(standalone_excitations),
                               {{TriggerMode::TimeDelay, INFINITY}});
    standalone.run(policy);

    const auto& independent_early = standalone.result(0);
    REQUIRE(independent_early.history.size() == early_result.history.size());
    CHECK(early_result.summary.step_count == independent_early.summary.step_count);
    CHECK(early_result.summary.muzzle_velocity ==
          doctest::Approx(independent_early.summary.muzzle_velocity).epsilon(1e-8));
    CHECK(early_result.history.front().state.arm_position ==
          doctest::Approx(independent_early.history.front().state.arm_position).epsilon(1e-8));
    CHECK(early_result.history.front().state.arm_velocity ==
          doctest::Approx(independent_early.history.front().state.arm_velocity).epsilon(1e-8));
    CHECK(early_result.history.front().state.force ==
          doctest::Approx(independent_early.history.front().state.force).epsilon(1e-8));
    CHECK(early_result.history.front().coil_currents == independent_early.history.front().coil_currents);
    CHECK(early_result.history.front().cap_voltages == independent_early.history.front().cap_voltages);

    SimBatch<EulerStepper> standalone_continuing(make_coils(), make_arm(), 1, 1e-6, direct_backend());
    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> standalone_continuing_excitations;
    standalone_continuing_excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    standalone_continuing_excitations.push_back(std::make_unique<CrowbarExcitation>(300.0, 0.001));
    standalone_continuing.set_excitations(0, std::move(standalone_continuing_excitations),
                                          {{TriggerMode::TimeDelay, 1.0}});
    standalone_continuing.run(policy);
    CHECK(continuing_result.summary.muzzle_velocity ==
          doctest::Approx(standalone_continuing.result(0).summary.muzzle_velocity).epsilon(1e-8));
}

TEST_CASE("SimBatch — heterogeneous triggers and excitation classes stay per-row") {
    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 4;

    std::vector<DrivingCoil> coils;
    for (int stage = 0; stage < 3; ++stage)
        coils.emplace_back(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref,
                           1e-6, 0.7, 0.10 * stage);
    auto arm = make_arm();
    SimBatch<EulerStepper> batch(coils, arm, 2, 1e-6, direct_backend());

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> early;
    early.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    early.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    early.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 0.0; }));
    batch.set_excitations(0, std::move(early),
                          {{TriggerMode::Position, 0.05},
                           {TriggerMode::TimeDelay, 1e-6}});

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> continuing;
    continuing.push_back(std::make_unique<coilgun::simulation::CapacitorExcitation>(400.0, 2e-3));
    continuing.push_back(std::make_unique<CrowbarExcitation>(300.0, 1e-3));
    continuing.push_back(std::make_unique<coilgun::simulation::WaveformExcitation>(
        [](double) { return 150.0; }));
    batch.set_excitations(1, std::move(continuing),
                          {{TriggerMode::TimeDelay, 2e-6},
                           {TriggerMode::Position, 0.05}});
    batch.run(policy);

    const auto& early_result = batch.result(0);
    const auto& continuing_result = batch.result(1);
    REQUIRE(early_result.summary.per_stage.size() == 3);
    REQUIRE(continuing_result.summary.per_stage.size() == 3);
    CHECK(early_result.summary.step_count < continuing_result.summary.step_count);
    CHECK(early_result.summary.per_stage[1].trigger_time == doctest::Approx(0.0));
    CHECK(early_result.summary.per_stage[1].trigger_position == doctest::Approx(0.05));
    CHECK(early_result.summary.per_stage[2].trigger_time == doctest::Approx(1e-6));
    CHECK(continuing_result.summary.per_stage[1].trigger_time == doctest::Approx(2e-6));
    CHECK(continuing_result.summary.per_stage[1].trigger_position == doctest::Approx(0.05));
    CHECK(continuing_result.summary.per_stage[2].trigger_position == doctest::Approx(0.05));
    CHECK(early_result.history.front().cap_voltages[0] == doctest::Approx(0.0));
    CHECK(continuing_result.history.front().cap_voltages[0] == doctest::Approx(400.0));
    CHECK(std::abs(early_result.summary.muzzle_velocity - continuing_result.summary.muzzle_velocity) > 1e-8);
    CHECK(early_result.summary.per_stage[0].energy_depleted == doctest::Approx(0.0));
    CHECK(continuing_result.summary.per_stage[0].energy_depleted != doctest::Approx(0.0));
}

TEST_CASE("SimBatch — canonical cutoff masks a distant active stage") {
    auto coils = make_coils();
    coils[1] = DrivingCoil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref,
                           1e-6, 0.7, 0.75);
    auto arm = make_arm();
    auto policy = coilgun::simulation::TerminationPolicy::defaults();
    policy.max_steps = 1;

    auto individual_excitations = std::vector<std::unique_ptr<coilgun::simulation::Excitation>>{};
    individual_excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    individual_excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    GpuMultiStageSim<EulerStepper> individual(
        coils, arm, std::move(individual_excitations),
        {{TriggerMode::TimeDelay, 0.0}}, 1e-6);
    individual.run(policy);

    SimBatch<EulerStepper> batch(coils, arm, 1, 1e-6, direct_backend());
    auto batch_excitations = std::vector<std::unique_ptr<coilgun::simulation::Excitation>>{};
    batch_excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    batch_excitations.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    batch.set_excitations(0, std::move(batch_excitations), {{TriggerMode::TimeDelay, 0.0}});
    batch.run(policy);

    REQUIRE(batch.result(0).history.size() == 1);
    REQUIRE(individual.result().history.size() == 1);
    CHECK(batch.result(0).history[0].stage_forces[1] ==
          doctest::Approx(individual.result().history[0].stage_forces[1]).epsilon(1e-5));
}

TEST_CASE("SimBatch — RK4 is explicitly rejected rather than treated as Euler") {
    CHECK_THROWS_AS(SimBatch<RK4Stepper>(make_coils(), make_arm(), 1, 1e-6, direct_backend()),
                    std::logic_error);
}

TEST_CASE("SimBatch — fixed indexing and validation") {
    auto coils = make_coils();
    auto arm = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 2, 1e-6, direct_backend());

    CHECK_THROWS_AS(batch.result(-1), std::out_of_range);
    CHECK_THROWS_AS(batch.result(2), std::out_of_range);
    CHECK_THROWS_AS(batch.run(), std::logic_error);

    std::vector<std::unique_ptr<coilgun::simulation::Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(350.0, 0.001));
    CHECK_THROWS_AS(batch.set_excitations(0, std::move(excs), {}), std::invalid_argument);
}
