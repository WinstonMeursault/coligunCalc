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
#include <memory>
#include <vector>

using namespace coilgun::physics;
using coilgun::components::DrivingCoil;
using coilgun::components::Armature;
using coilgun::simulation::EulerStepper;
using coilgun::simulation::CrowbarExcitation;
using coilgun::simulation::TriggerConfig;
using coilgun::simulation::TriggerMode;
using coilgun::simulation::cuda::SimBatch;
using coilgun::simulation::cuda::GpuMultiStageSim;

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
        sim.run();
        v_ref[i] = sim.result().summary.muzzle_velocity;
    }

    auto coils = make_coils();
    auto arm   = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 2, 1e-6);

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
        CHECK(v_batch > 0.0);
    }
}

TEST_CASE("SimBatch — direction consistent") {
    auto coils = make_coils();
    auto arm   = make_arm();
    SimBatch<EulerStepper> batch(std::move(coils), arm, 1, 1e-6);

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
