#include <doctest/doctest.h>
#include "coilgun/coilgun.hpp"
#include <algorithm>
#include <cmath>
#include <memory>

using namespace coilgun::simulation;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::physics::COPPER;
using coilgun::physics::ALUMINUM;

struct FastConfig {
    DrivingCoil coil{0.005, 0.010, 0.010, 10,
                     COPPER.resistivity_ref, 1e-6, 0.7};
    Armature arm{0.002, 0.008, 0.010,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.005, 1, 1, 0.008};
    double dt  = 1e-6;
    double U0  = 100.0;
    double C   = 100e-6;
};

static double integrate_joules(const SimResult& r,
                               const DrivingCoil& coil,
                               const Armature& arm,
                               double dt) {
    double Q_d = 0.0;
    double Q_a = 0.0;
    double R_d = coil.resistance();
    const auto& R_fil = arm.resistances();
    for (const auto& s : r.history) {
        Q_d += s.coil_current * s.coil_current * R_d * dt;
        for (std::size_t k = 0; k < s.filament_currents.size(); ++k)
            Q_a += s.filament_currents[k] * s.filament_currents[k] * R_fil[k] * dt;
    }
    return Q_d + Q_a;
}

TEST_CASE("SingleStageSim — energy conservation") {
    FastConfig cfg;
    auto exc = std::make_unique<CrowbarExcitation>(cfg.U0, cfg.C);
    SingleStageSim<EulerStepper> sim(cfg.coil, cfg.arm, std::move(exc), cfg.dt, false);
    sim.run();

    double E_cap = 0.5 * cfg.C * cfg.U0 * cfg.U0;
    double Q_joule = integrate_joules(sim.result(), cfg.coil, cfg.arm, cfg.dt);
    double v = sim.result().summary.muzzle_velocity;
    double E_kin = 0.5 * cfg.arm.mass() * v * v;
    double closure = (Q_joule + E_kin) / E_cap;

    CHECK(closure > 0.99);
    CHECK(closure < 1.01);
}

TEST_CASE("SingleStageSim — convergence (Euler)") {
    FastConfig cfg;
    double v_coarse = 0.0, v_fine = 0.0;

    {
        auto exc = std::make_unique<CrowbarExcitation>(cfg.U0, cfg.C);
        SingleStageSim<EulerStepper> sim(cfg.coil, cfg.arm, std::move(exc), 1e-6, false);
        sim.run();
        v_coarse = sim.result().summary.muzzle_velocity;
    }
    {
        auto exc = std::make_unique<CrowbarExcitation>(cfg.U0, cfg.C);
        SingleStageSim<EulerStepper> sim(cfg.coil, cfg.arm, std::move(exc), 5e-7, false);
        sim.run();
        v_fine = sim.result().summary.muzzle_velocity;
    }

    double rel_diff = std::abs(v_coarse - v_fine) / v_fine;
    CHECK(rel_diff < 0.02);
}

TEST_CASE("SingleStageSim — summary max_force is the peak absolute force") {
    DrivingCoil coil(0.01, 0.03, 0.05, 150, COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08, ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, -0.1);
    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false);
    TerminationPolicy policy;
    policy.max_steps = 3;
    policy.enable_velocity_check = false;
    sim.run(policy);

    double expected_max_force = 0.0;
    bool saw_negative_force = false;
    for (const auto& step : sim.result().history) {
        expected_max_force = std::max(expected_max_force, std::abs(step.force));
        saw_negative_force = saw_negative_force || step.force < 0.0;
    }
    CHECK(saw_negative_force);
    CHECK(sim.result().summary.max_force == doctest::Approx(expected_max_force));
}

TEST_CASE("SingleStageSim — inductance benchmark (Paper 4 Ex.1)") {
    double ri = 0.04, re = 0.06, length = 0.20;
    double nc = 250000.0;
    int N = static_cast<int>(nc * (re - ri) * length);
    double L_ref = 34.66e-3;

    DrivingCoil coil(ri, re, length, N, COPPER.resistivity_ref, 1e-6, 0.7);
    CHECK(coil.self_inductance() == doctest::Approx(L_ref).epsilon(5e-4));

    Armature arm(0.035, 0.055, 0.02,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.1, 1, 1, 0.22);

    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6, false);
    sim.run();

    CHECK(sim.step_count() > 0);
    CHECK(sim.result().summary.muzzle_velocity > 0.0);
}

// =========================================================================
// 82 mm single-stage coilgun — validation against Xiang (2015) §9.2
//
// Reference data from Ansoft FEM; our CFM solver is a different numerical
// method, so thresholds are generous (order-of-magnitude validation).
// =========================================================================

TEST_CASE("SingleStageSim — 82mm optimised design (Xiang §9.2.2)") {
    // Final optimised parameters — expected v ≈ 105.9 m/s, η ≈ 12.02%.

    constexpr double coil_ri = 0.043, coil_re = 0.060, coil_len = 0.080;
    constexpr int    coil_turns = 77;
    constexpr double wire_area  = 6e-6;    // 3×2 mm copper core
    constexpr double fill_factor = 0.34;    // 77×6 / (17×80)
    constexpr double coil_pos   = 0.040;    // centre

    DrivingCoil coil(coil_ri, coil_re, coil_len, coil_turns,
                     COPPER.resistivity_ref, wire_area, fill_factor,
                     coil_pos);

    constexpr double arm_ri = 0.026, arm_re = 0.041, arm_len = 0.080;
    constexpr double arm_mass = 3.0;
    constexpr double s = 0.035;             // trigger position

    Armature arm(arm_ri, arm_re, arm_len,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, arm_mass, 15, 3, coil_pos + s);

    constexpr double U0 = 10000.0;
    constexpr double C0 = 2800e-6;
    constexpr double dt = 5e-6;
    constexpr double barrel_end = 0.300;  // stop after armature exits

    auto exc = std::make_unique<CrowbarExcitation>(U0, C0);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt, false);
    TerminationPolicy term;
    term.enable_bound_check = true;
    term.barrel_end_position = barrel_end;
    term.max_steps = 80000;
    sim.run(term);

    const auto& r = sim.result();
    double v   = r.summary.muzzle_velocity;
    double eff = r.summary.efficiency;

    // CFM vs FEM: systematic differences expected
    CHECK(v   > 50.0);    // ref ≈ 105.9
    CHECK(v   < 180.0);
    CHECK(eff > 0.04);    // ref ≈ 12.02 %
    CHECK(eff < 0.20);
    CHECK(r.summary.step_count > 0);
    CHECK(r.summary.step_count < 80000);
}

TEST_CASE("SingleStageSim — 82mm opt round 1 (Xiang Tab.9-5 row 1)") {
    // Table 9-5 row 1: r₂=60, N=70, C=2800, s=35  →  η ≈ 11.84 %
    // Fixed: a=80, b=80, U=10000, r₃=26, r₄=41, r₁=43, m=3

    constexpr double coil_ri = 0.043, coil_re = 0.060, coil_len = 0.080;
    constexpr int    coil_turns = 70;
    constexpr double wire_area  = 6e-6;
    constexpr double fill_factor = 0.34;
    constexpr double coil_pos   = 0.040;

    DrivingCoil coil(coil_ri, coil_re, coil_len, coil_turns,
                     COPPER.resistivity_ref, wire_area, fill_factor,
                     coil_pos);

    constexpr double arm_ri = 0.026, arm_re = 0.041, arm_len = 0.080;
    constexpr double s = 0.035;

    Armature arm(arm_ri, arm_re, arm_len,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 3.0, 15, 3, coil_pos + s);

    constexpr double U0 = 10000.0;
    constexpr double C0 = 2800e-6;
    constexpr double dt = 5e-6;
    constexpr double barrel_end = 0.300;

    auto exc = std::make_unique<CrowbarExcitation>(U0, C0);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt, false);
    TerminationPolicy term;
    term.enable_bound_check = true;
    term.barrel_end_position = barrel_end;
    term.max_steps = 80000;
    sim.run(term);

    const auto& r = sim.result();
    double v   = r.summary.muzzle_velocity;
    double eff = r.summary.efficiency;

    CHECK(v   > 50.0);     // ref: body velocity from η=11.84%
    CHECK(v   < 180.0);
    CHECK(eff > 0.04);     // ref ≈ 11.84 %
    CHECK(eff < 0.20);
    CHECK(r.summary.step_count > 0);
    CHECK(r.summary.step_count < 80000);
}

TEST_CASE("SingleStageSim — 82mm L32 run 3 (Xiang Tab.9-2 row 3)") {
    // Table 9-2 row 3: r₂=73, a=110, N=80, r₃=26, b=80,
    //                  C=2500, U=8000, s=30  →  η ≈ 8.59 %

    constexpr double coil_ri = 0.043, coil_re = 0.073;
    constexpr double coil_len   = 0.110;
    constexpr int    coil_turns = 80;
    constexpr double wire_area  = 6e-6;
    constexpr double fill_factor = 0.34;
    constexpr double coil_pos   = 0.055;        // len / 2

    DrivingCoil coil(coil_ri, coil_re, coil_len, coil_turns,
                     COPPER.resistivity_ref, wire_area, fill_factor,
                     coil_pos);

    constexpr double arm_ri = 0.026, arm_re = 0.041, arm_len = 0.080;
    constexpr double s = 0.030;

    Armature arm(arm_ri, arm_re, arm_len,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 3.0, 15, 3, coil_pos + s);

    constexpr double U0 = 8000.0;
    constexpr double C0 = 2500e-6;
    constexpr double dt = 5e-6;
    constexpr double barrel_end = 0.300;

    auto exc = std::make_unique<CrowbarExcitation>(U0, C0);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt, false);
    TerminationPolicy term;
    term.enable_bound_check = true;
    term.barrel_end_position = barrel_end;
    term.max_steps = 80000;
    sim.run(term);

    const auto& r = sim.result();
    double v   = r.summary.muzzle_velocity;
    double eff = r.summary.efficiency;

    CHECK(v   > 40.0);     // ref η ≈ 8.59 %
    CHECK(v   < 150.0);
    CHECK(eff > 0.03);
    CHECK(eff < 0.15);
    CHECK(r.summary.step_count > 0);
    CHECK(r.summary.step_count < 80000);
}
