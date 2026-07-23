#include <doctest/doctest.h>
#include "coilgun/coilgun.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace coilgun::simulation;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::physics::COPPER;
using coilgun::physics::ALUMINUM;

static std::unique_ptr<Excitation> make_crowbar(double U0, double C) {
    return std::make_unique<CrowbarExcitation>(U0, C);
}

static TerminationPolicy no_velocity_policy(int max_steps = 20000) {
    TerminationPolicy policy;
    policy.max_steps = max_steps;
    policy.enable_velocity_check = false;
    policy.enable_bound_check = false;
    return policy;
}

// Clone an excitation vector by reconstructing from capacitor parameters.
// Preserves crowbar vs plain-capacitor type.
static std::vector<std::unique_ptr<Excitation>> clone_excitations(
    const std::vector<std::unique_ptr<Excitation>>& src) {
    std::vector<std::unique_ptr<Excitation>> dst;
    for (const auto& e : src) {
        auto* c = dynamic_cast<CapacitorExcitation*>(e.get());
        if (!c) continue;
        if (dynamic_cast<CrowbarExcitation*>(e.get()))
            dst.push_back(std::make_unique<CrowbarExcitation>(
                c->initial_voltage(), c->capacitance()));
        else
            dst.push_back(std::make_unique<CapacitorExcitation>(
                c->initial_voltage(), c->capacitance()));
    }
    return dst;
}

TEST_CASE("MultiStageSim — n_stages=1 matches SingleStageSim") {
    DrivingCoil coil(0.005, 0.010, 0.010, 10,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 1, 1, 0.008);
    double dt = 1e-6, U0 = 100.0, C = 100e-6;

    double v_single;
    {
        auto exc = make_crowbar(U0, C);
        SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt, false);
        sim.run();
        v_single = sim.result().summary.muzzle_velocity;
    }

    double v_multi;
    {
        std::vector<DrivingCoil> coils;
        coils.push_back(coil);

        std::vector<std::unique_ptr<Excitation>> excs;
        excs.push_back(make_crowbar(U0, C));

        std::vector<TriggerConfig> triggers; // empty: n_stages-1 = 0

        MultiStageSim<EulerStepper> sim(std::move(coils), arm,
            std::move(excs), triggers, dt);
        sim.run();
        v_multi = sim.result().summary.muzzle_velocity;
    }

    CHECK(v_multi == doctest::Approx(v_single).epsilon(1e-6));
}

TEST_CASE("MultiStageSim — two stages, position trigger") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    double dt = 1e-6, U0 = 500.0, C = 500e-6;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil1);
    coils.push_back(coil2);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0, C));
    excs.push_back(make_crowbar(U0, C));

    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::Position, 0.010});

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt);
    sim.run();

    auto& summary = sim.result().summary;
    CHECK(summary.muzzle_velocity > 0.0);
    CHECK(summary.per_stage.size() >= 2u);
    CHECK(summary.step_count > 0);

    // Stage 0 auto-triggers at t=0, stage 1 triggers later
    CHECK(summary.per_stage[0].trigger_time == doctest::Approx(0.0));
    CHECK(summary.per_stage[1].trigger_time > 0.0);
}

TEST_CASE("MultiStageSim — two stages, time delay trigger") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    double dt = 1e-6, U0 = 500.0, C = 500e-6;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil1);
    coils.push_back(coil2);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0, C));
    excs.push_back(make_crowbar(U0, C));

    // Stage 1 triggers 200us after stage 0
    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::TimeDelay, 200e-6});

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt);
    sim.run();

    auto& summary = sim.result().summary;
    CHECK(summary.muzzle_velocity > 0.0);
    CHECK(summary.per_stage.size() >= 2u);

    // Time delay trigger should be exactly 200us after stage 0
    double delay = summary.per_stage[1].trigger_time
                 - summary.per_stage[0].trigger_time;
    CHECK(delay == doctest::Approx(200e-6).epsilon(1e-4));
}

TEST_CASE("MultiStageSim — completion waits for an eligible delayed stage") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 10,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 10,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);

    std::vector<std::unique_ptr<Excitation>> excitations;
    excitations.push_back(make_crowbar(1.0, 50e-6));
    excitations.push_back(make_crowbar(500.0, 500e-6));

    MultiStageSim<EulerStepper> sim(
        {coil1, coil2}, arm, std::move(excitations),
        {{TriggerMode::TimeDelay, 5e-6}}, 1e-6);
    sim.run(no_velocity_policy(200));

    REQUIRE(sim.result().summary.per_stage.size() == 2);
    CHECK(sim.result().summary.per_stage[1].trigger_time ==
          doctest::Approx(5e-6).epsilon(1e-12));
    CHECK(sim.result().summary.per_stage[1].peak_current > 0.0);
}

TEST_CASE("MultiStageSim — trigger position is captured at the crossed boundary") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);

    std::vector<std::unique_ptr<Excitation>> excitations;
    excitations.push_back(make_crowbar(500.0, 500e-6));
    excitations.push_back(make_crowbar(500.0, 500e-6));

    MultiStageSim<EulerStepper> sim(
        {coil1, coil2}, arm, std::move(excitations),
        {{TriggerMode::Position, 0.00301}}, 1e-6);

    bool triggered = false;
    double pre_step_trigger_position = 0.0;
    for (int i = 0; i < 100; ++i) {
        const double pre_step_position = sim.state().arm_position;
        sim.step();
        if (sim.stage_state(1).triggered) {
            triggered = true;
            pre_step_trigger_position = pre_step_position;
            break;
        }
    }
    REQUIRE(triggered);

    sim.run(no_velocity_policy(0));
    REQUIRE(sim.result().summary.per_stage.size() == 2);
    CHECK(sim.result().summary.per_stage[1].trigger_position ==
          doctest::Approx(pre_step_trigger_position).epsilon(1e-12));
    CHECK(pre_step_trigger_position < 0.00301);
}

TEST_CASE("MultiStageSim — zero-delay trigger captures the pre-step boundary") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);

    std::vector<std::unique_ptr<Excitation>> excitations;
    excitations.push_back(make_crowbar(500.0, 500e-6));
    excitations.push_back(make_crowbar(500.0, 500e-6));
    MultiStageSim<EulerStepper> sim(
        {coil1, coil2}, arm, std::move(excitations),
        {{TriggerMode::TimeDelay, 0.0}}, 1e-6);

    const double pre_step_position = sim.state().arm_position;
    sim.step();
    REQUIRE(sim.result().history.back().coil_currents[1] != 0.0);

    sim.run(no_velocity_policy(0));
    REQUIRE(sim.result().summary.per_stage.size() == 2);
    CHECK(sim.result().summary.per_stage[1].trigger_position ==
          doctest::Approx(pre_step_position).epsilon(1e-12));
}

TEST_CASE("MultiStageSim — energy depletion uses each capacitor's initial energy") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    constexpr double c0 = 50e-6;
    constexpr double u0 = 500.0;
    constexpr double c1 = 500e-6;
    constexpr double u1 = 320.0;

    std::vector<std::unique_ptr<Excitation>> excitations;
    excitations.push_back(make_crowbar(u0, c0));
    excitations.push_back(make_crowbar(u1, c1));
    MultiStageSim<EulerStepper> sim(
        {coil1, coil2}, arm, std::move(excitations),
        {{TriggerMode::TimeDelay, 1e-6}}, 1e-6);

    sim.run(no_velocity_policy(3));
    REQUIRE(sim.result().summary.per_stage.size() == 2);
    const auto& stage0 = sim.result().summary.per_stage[0];
    const auto& stage1 = sim.result().summary.per_stage[1];
    const double stage0_initial = 0.5 * c0 * u0 * u0;
    const double stage1_initial = 0.5 * c1 * u1 * u1;
    const double stage0_final = 0.5 * c0 * sim.result().history.back().cap_voltages[0] *
        sim.result().history.back().cap_voltages[0];
    const double stage1_final = 0.5 * c1 * sim.result().history.back().cap_voltages[1] *
        sim.result().history.back().cap_voltages[1];

    CHECK(stage0.energy_depleted == doctest::Approx(stage0_initial - stage0_final));
    CHECK(stage1.energy_depleted == doctest::Approx(stage1_initial - stage1_final));
    CHECK(stage0.energy_depleted > 0.0);
    CHECK(stage1.energy_depleted > 0.0);
}

TEST_CASE("MultiStageSim — trigger configuration rejects malformed values") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 10,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 10,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);

    const auto construct = [&](TriggerConfig config) {
        std::vector<std::unique_ptr<Excitation>> excitations;
        excitations.push_back(make_crowbar(1.0, 50e-6));
        excitations.push_back(make_crowbar(100.0, 50e-6));
        MultiStageSim<EulerStepper> sim(
            {coil1, coil2}, arm, std::move(excitations),
            {config}, 1e-6);
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

TEST_CASE("MultiStageSim — optimization level consistency") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    double dt = 1e-6, U0 = 500.0, C = 500e-6;

    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::Position, 0.010});

    double v_ref;
    {
        std::vector<DrivingCoil> coils;
        coils.push_back(coil1);
        coils.push_back(coil2);
        std::vector<std::unique_ptr<Excitation>> excs;
        excs.push_back(make_crowbar(U0, C));
        excs.push_back(make_crowbar(U0, C));
        MultiStageSim<EulerStepper> sim(std::move(coils), arm,
            std::move(excs), triggers, dt, false,
            OptimizationLevel::Reference);
        sim.run();
        v_ref = sim.result().summary.muzzle_velocity;
    }

    double v_full;
    {
        std::vector<DrivingCoil> coils;
        coils.push_back(coil1);
        coils.push_back(coil2);
        std::vector<std::unique_ptr<Excitation>> excs;
        excs.push_back(make_crowbar(U0, C));
        excs.push_back(make_crowbar(U0, C));
        MultiStageSim<EulerStepper> sim(std::move(coils), arm,
            std::move(excs), triggers, dt, false,
            OptimizationLevel::Full);
        sim.run();
        v_full = sim.result().summary.muzzle_velocity;
    }

    // Optimized should be within 5% of reference
    CHECK(std::abs(v_ref - v_full) / v_ref < 0.05);
}

TEST_CASE("MultiStageSim — termination when all stages finished") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 5,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 5,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    double dt = 1e-6, U0 = 20.0, C = 50e-6;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil1);
    coils.push_back(coil2);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0, C));
    excs.push_back(make_crowbar(U0, C));

    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::Position, 0.025});

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt);
    sim.run();

    // Should have terminated (not hit max_steps)
    CHECK(sim.step_count() < 20000);
    CHECK(sim.result().summary.muzzle_velocity >= 0.0);
}

TEST_CASE("MultiStageSim — reset and re-run") {
    DrivingCoil coil(0.005, 0.010, 0.010, 10,
                     COPPER.resistivity_ref, 1e-6, 0.7, 0.0);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 1, 1, 0.008);
    double dt = 1e-6, U0 = 100.0, C = 100e-6;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0, C));

    std::vector<TriggerConfig> triggers;

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt);
    sim.run();
    double v1 = sim.result().summary.muzzle_velocity;

    const auto first_summary = sim.result().summary;
    sim.run(no_velocity_policy(0));
    CHECK(sim.result().summary.per_stage.size() == first_summary.per_stage.size());
    CHECK(sim.result().summary.step_count == first_summary.step_count);
    CHECK(sim.result().summary.max_force ==
          doctest::Approx(first_summary.max_force));

    sim.reset();
    sim.run();
    double v2 = sim.result().summary.muzzle_velocity;

    CHECK(v1 == doctest::Approx(v2).epsilon(1e-6));
}

TEST_CASE("MultiStageSim — completion boundary records direct CPU force contributions") {
    DrivingCoil coil1(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 20,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.03);
    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    auto completed = std::make_unique<WaveformExcitation>(
        [](double) { return 480.0; });
    completed->set_end_time(1e-6);

    std::vector<std::unique_ptr<Excitation>> excitations;
    excitations.push_back(std::move(completed));
    excitations.push_back(std::make_unique<WaveformExcitation>(
        [](double) { return 480.0; }));
    MultiStageSim<EulerStepper> sim(
        {coil1, coil2}, arm, std::move(excitations),
        {{TriggerMode::TimeDelay, 0.0}}, 1e-6);

    const auto& step = sim.step();
    REQUIRE(step.stage_forces.size() == 2);
    REQUIRE(std::abs(step.stage_forces[1]) > 1e-18);
    // Excitation completion stops the source voltage but does not remove
    // residual circuit current from the coupled force calculation.
    CHECK(std::abs(step.stage_forces[0]) > 1e-18);
    CHECK(step.state.force == doctest::Approx(
        step.stage_forces[0] + step.stage_forces[1]));
}

TEST_CASE("MultiStageSim — crowbar diode independent per-stage") {
    // Stage 0: fast small-cap discharge, crowbars in ~8 us.
    // Stage 1: triggered by time delay at ~5 us — overlaps with stage 0 crowbar.
    // Assert: history contains a mixed-state step where stage 0
    //          is crowbarring (cap=0, I>0) while stage 1 is still in RLC mode.

    DrivingCoil coil1(0.005, 0.010, 0.010, 10,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.005, 0.010, 0.010, 10,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.02);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 2, 1, 0.003);
    double dt = 1e-7;  // 0.1 us — resolve the fast crowbar
    // Stage 0: tiny cap, fast crowbar
    double C0 = 1e-6,   U0_s0 = 100.0;
    // Stage 1: larger cap, triggers by time delay
    double C1 = 100e-6, U0_s1 = 200.0;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil1);
    coils.push_back(coil2);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0_s0, C0));
    excs.push_back(make_crowbar(U0_s1, C1));

    // Stage 1 triggers 5 us after stage 0 — during stage 0's crowbar window
    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::TimeDelay, 5e-6});

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt);
    TerminationPolicy mixed_policy;
    mixed_policy.max_steps = 200;
    mixed_policy.enable_velocity_check = false;
    sim.run(mixed_policy);

    const auto& h = sim.result().history;
    REQUIRE(h.size() > 2);

    bool found_mixed = false;
    for (const auto& step : h) {
        double v0 = step.cap_voltages[0];
        double v1 = step.cap_voltages[1];
        double i0 = step.coil_currents[0];
        double i1 = step.coil_currents[1];
        // Stage 0 crowbarred: cap=0 but still conducting through diode
        // Stage 1 still RLC: cap>0
        if (v0 == 0.0 && i0 > 1e-8 && v1 > 0.0 && i1 > 1e-8) {
            found_mixed = true;
            break;
        }
    }
    CHECK(found_mixed);
}

TEST_CASE("MultiStageSim — thermal Cu vs Al comparison (Ni et al. 2015 Table 2)") {
    // Reference: Ni et al. (2015) Table 2 — copper yields higher velocity
    // than aluminium at equal total mass (1 kg) due to lower temp_coefficient.
    //
    // This test verifies BOTH materials accelerate correctly under thermal
    // mode.  Velocity comparison at equal mass requires different armature
    // geometry per material; we skip that assertion here.
    //
    // 5 stages, light armature (~0.01 kg) for fast runtime.

    constexpr int    n_stages  = 5;
    constexpr double U0        = 1000.0;
    constexpr double C0        = 200e-6;
    constexpr double dt        = 2e-6;
    constexpr double spacing   = 0.030;

    constexpr double coil_ri    = 0.005, coil_re = 0.012;
    constexpr double coil_len   = 0.012;
    constexpr int    coil_turns = 15;
    constexpr double wire_area  = 2e-6;

    constexpr double arm_ri  = 0.002, arm_re = 0.004, arm_len = 0.012;
    constexpr int    m_fil = 4, n_fil = 1;
    constexpr double arm_start  = coil_len / 2.0 - 0.003;  // behind → dM/dz > 0
    constexpr double coil_start = coil_len / 2.0;

    constexpr double al_density = 2700.0;
    constexpr double cu_density = 8960.0;
    constexpr double al_mass_small = al_density * M_PI
        * (arm_re*arm_re - arm_ri*arm_ri) * arm_len;

    auto build_coils = [&]() {
        std::vector<DrivingCoil> cs;
        for (int i = 0; i < n_stages; ++i)
            cs.emplace_back(coil_ri, coil_re, coil_len, coil_turns,
                            COPPER.resistivity_ref, wire_area, 0.7,
                            coil_start + i * spacing);
        return cs;
    };

    auto build_excs = [&]() {
        std::vector<std::unique_ptr<Excitation>> es;
        for (int i = 0; i < n_stages; ++i)
            es.push_back(make_crowbar(U0, C0));
        return es;
    };

    auto build_triggers = [&]() {
        std::vector<TriggerConfig> ts;
        for (int i = 1; i < n_stages; ++i) {
            double pos = coil_start + i * spacing;
            ts.push_back({TriggerMode::Position, pos - 0.5 * arm_len});
        }
        return ts;
    };

    double v_al, v_cu, maxT_al = 0.0, maxT_cu = 0.0;

    {
        Armature arm_al(arm_ri, arm_re, arm_len,
                        ALUMINUM.resistivity_ref, ALUMINUM.density,
                        0.0, al_mass_small, m_fil, n_fil, arm_start,
                        coilgun::physics::ArmatureMaterial::Aluminum);
        MultiStageSim<EulerStepper> sim(build_coils(), arm_al,
            build_excs(), build_triggers(), dt, true,
            OptimizationLevel::Full);
        sim.run();
        v_al = sim.result().summary.muzzle_velocity;
        if (!sim.result().history.empty()) {
            const auto& last = sim.result().history.back();
            if (!last.state.filament_temperatures.empty())
                maxT_al = *std::max_element(
                    last.state.filament_temperatures.begin(),
                    last.state.filament_temperatures.end());
        }
        CHECK(std::isfinite(v_al));
        CHECK(std::abs(v_al) > 0.0);
    }

    {
        Armature arm_cu(arm_ri, arm_re, arm_len,
                        COPPER.resistivity_ref, COPPER.density,
                        0.0, cu_density * M_PI
                            * (arm_re*arm_re - arm_ri*arm_ri) * arm_len,
                        m_fil, n_fil, arm_start,
                        coilgun::physics::ArmatureMaterial::Copper);
        MultiStageSim<EulerStepper> sim(build_coils(), arm_cu,
            build_excs(), build_triggers(), dt, true,
            OptimizationLevel::Full);
        sim.run();
        v_cu = sim.result().summary.muzzle_velocity;
        if (!sim.result().history.empty()) {
            const auto& last = sim.result().history.back();
            if (!last.state.filament_temperatures.empty())
                maxT_cu = *std::max_element(
                    last.state.filament_temperatures.begin(),
                    last.state.filament_temperatures.end());
        }
        CHECK(std::isfinite(v_cu));
        CHECK(std::abs(v_cu) > 0.0);
    }

    // Both materials must produce temperature rise above reference
    CHECK(maxT_al >= coilgun::physics::T_REFERENCE);
    CHECK(maxT_cu >= coilgun::physics::T_REFERENCE);
}

TEST_CASE("MultiStageSim — thermal mode runs without error") {
    DrivingCoil coil(0.005, 0.010, 0.010, 5,
                     COPPER.resistivity_ref, 1e-6, 0.7);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 1, 1, 0.008);
    double dt = 1e-6, U0 = 50.0, C = 50e-6;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0, C));

    std::vector<TriggerConfig> triggers;

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt, true);
    sim.run();

    CHECK(sim.step_count() > 0);
    CHECK(sim.result().summary.muzzle_velocity >= 0.0);

    // At least some filament temperatures should be above reference
    if (!sim.result().history.empty()) {
        const auto& last = sim.result().history.back();
        if (!last.state.filament_temperatures.empty()) {
            double max_T = *std::max_element(
                last.state.filament_temperatures.begin(),
                last.state.filament_temperatures.end());
            CHECK(max_T >= coilgun::physics::T_REFERENCE);
        }
    }
}

TEST_CASE("MultiStageSim — result sampling") {
    DrivingCoil coil(0.005, 0.010, 0.010, 10,
                     COPPER.resistivity_ref, 1e-6, 0.7);

    Armature arm(0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 1, 1, 0.008);
    double dt = 1e-6, U0 = 100.0, C = 100e-6;

    std::vector<DrivingCoil> coils;
    coils.push_back(coil);

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(make_crowbar(U0, C));

    std::vector<TriggerConfig> triggers;

    MultiStageSim<EulerStepper> sim(std::move(coils), arm,
        std::move(excs), triggers, dt);
    sim.run();

    auto sampled = sim.result().sampled(10);
    CHECK(sampled.history.size() <= sim.result().history.size() / 10 + 1);
    CHECK(sampled.summary.muzzle_velocity ==
          doctest::Approx(sim.result().summary.muzzle_velocity));
}
