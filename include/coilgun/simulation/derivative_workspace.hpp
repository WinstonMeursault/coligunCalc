#pragma once

#include "coilgun/simulation/cpu_phase_timing.hpp"
#include "coilgun/simulation/excitation_snapshot.hpp"
#include "coilgun/simulation/integration_state.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <vector>

namespace coilgun::simulation {

enum class EventType {
    CapacitorZero,
    CrowbarTransition,
    WaveformEnd,
    ExcitationFinished,
    StageTrigger,
    CurrentDecay,
    StageCompleted,
};

struct EventObservation {
    EventType type = EventType::ExcitationFinished;
    double value = 0.0;
    double normalized_time = 0.0;
    int stage = -1;
};

struct DerivativeWorkspace {
    Eigen::VectorXd mutual;
    Eigen::VectorXd mutual_gradient;
    Eigen::MatrixXd system_matrix;
    Eigen::VectorXd rhs;
    Eigen::VectorXd resistance;
    std::vector<int> active_stages;

    void resize(std::size_t stages, std::size_t filaments);
};

struct DerivativeResult {
    ContinuousState physical_derivative;
    std::vector<ExcitationDerivative> excitation_derivatives;
    Eigen::VectorXd mutual_gradient;
    double force = 0.0;
    std::vector<EventObservation> events;
};

} // namespace coilgun::simulation
