#include <doctest/doctest.h>

#include "coilgun/coilgun.hpp"

#include <stdexcept>

using namespace coilgun::simulation;

TEST_CASE("sampling rejects non-positive stride") {
    SimResult single;
    MultiStageResult multi;
    CHECK_THROWS_AS(single.sampled(0), std::invalid_argument);
    CHECK_THROWS_AS(single.sampled(-1), std::invalid_argument);
    CHECK_THROWS_AS(multi.sampled(0), std::invalid_argument);
}

TEST_CASE("excitation snapshots restore all runtime state") {
    CapacitorExcitation capacitor(100.0, 1.0);
    auto capacitor_state = capacitor.snapshot();
    capacitor.advance(0.25, 20.0);
    REQUIRE(capacitor.capacitor_voltage() == doctest::Approx(95.0));
    capacitor.restore(*capacitor_state);
    CHECK(capacitor.capacitor_voltage() == doctest::Approx(100.0));

    CrowbarExcitation crowbar(1.0, 1.0);
    crowbar.advance(2.0, 1.0);
    REQUIRE(crowbar.diode_on());
    auto crowbar_state = crowbar.snapshot();
    crowbar.reset();
    crowbar.restore(*crowbar_state);
    CHECK(crowbar.diode_on());

    WaveformExcitation waveform([](double time) { return time; });
    waveform.set_end_time(1.0);
    waveform.advance(0.25, 0.0);
    auto waveform_state = waveform.snapshot();
    waveform.advance(1.0, 0.0);
    waveform.restore(*waveform_state);
    CHECK(waveform.voltage() == doctest::Approx(0.25));
    CHECK_FALSE(waveform.finished());
}

TEST_CASE("stage completion clears current before deactivation") {
    IntegrationState state;
    state.physical.currents = Eigen::VectorXd::Constant(1, 4.0);
    state.stages.resize(1);
    state.stages[0].triggered = true;
    state.stages[0].circuit_active = true;
    mark_stage_completed(state, 0);
    CHECK(state.physical.currents(0) == 0.0);
    CHECK(state.stages[0].stage_completed);
    CHECK_FALSE(state.stages[0].circuit_active);
}

TEST_CASE("derivative workspace owns dimensioned scratch") {
    DerivativeWorkspace workspace;
    workspace.resize(2, 3);
    CHECK(workspace.mutual.size() == 6);
    CHECK(workspace.system_matrix.rows() == 5);
    CHECK(workspace.rhs.size() == 5);
}
