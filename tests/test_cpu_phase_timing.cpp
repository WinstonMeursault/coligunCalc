#include <doctest/doctest.h>

#include "coilgun/coilgun.hpp"

using namespace coilgun::simulation;
using namespace coilgun::components;
using namespace coilgun::physics;

TEST_CASE("CPU phase timing scope accounts for nested phases") {
    CpuPhaseTiming timing;
    {
        CpuPhaseTimingCollector collector(timing);
        {
            CpuPhaseTimingScope orchestration(CpuPhase::Orchestration);
            CpuPhaseTimingScope mutual(CpuPhase::Mutual);
        }
    }

#if COILGUN_ENABLE_CPU_PHASE_TIMING
    CHECK(timing.derivative_count == 1);
    CHECK(timing.nanoseconds(CpuPhase::Orchestration) >=
          timing.nanoseconds(CpuPhase::Mutual));
#else
    CHECK(timing.derivative_count == 0);
    CHECK(timing.nanoseconds(CpuPhase::Orchestration) == 0);
    CHECK(timing.nanoseconds(CpuPhase::Mutual) == 0);
#endif
}

TEST_CASE("CPU phase timing records simulation derivative phases when enabled") {
    DrivingCoil coil(
        0.005, 0.010, 0.010, 10, COPPER.resistivity_ref,
        1e-6, 0.7);
    Armature armature(
        0.002, 0.008, 0.010, ALUMINUM.resistivity_ref,
        ALUMINUM.density, 0.0, 0.005, 1, 1, 0.008);
    auto excitation = std::make_unique<WaveformExcitation>(
        [](double) { return 100.0; });
    excitation->set_end_time(1e-3);
    SingleStageSim<EulerStepper> simulation(
        coil, armature, std::move(excitation), 1e-6, true);

    CpuPhaseTiming timing;
    {
        CpuPhaseTimingCollector collector(timing);
        simulation.step();
    }

#if COILGUN_ENABLE_CPU_PHASE_TIMING
    MESSAGE("derivatives=" << timing.derivative_count
            << " orchestration_ns=" << timing.nanoseconds(CpuPhase::Orchestration)
            << " mutual_ns=" << timing.nanoseconds(CpuPhase::Mutual)
            << " assembly_ns=" << timing.nanoseconds(CpuPhase::Assembly)
            << " solve_ns=" << timing.nanoseconds(CpuPhase::Solve)
            << " thermal_ns=" << timing.nanoseconds(CpuPhase::Thermal));
    CHECK(timing.derivative_count > 0);
    CHECK(timing.nanoseconds(CpuPhase::Orchestration) > 0);
    CHECK(timing.nanoseconds(CpuPhase::Mutual) > 0);
    CHECK(timing.nanoseconds(CpuPhase::Assembly) > 0);
    CHECK(timing.nanoseconds(CpuPhase::Solve) > 0);
    CHECK(timing.nanoseconds(CpuPhase::Thermal) > 0);
#else
    CHECK(timing.derivative_count == 0);
    CHECK(timing.nanoseconds(CpuPhase::Orchestration) == 0);
#endif
}
