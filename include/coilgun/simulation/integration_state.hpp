#pragma once

#include "coilgun/simulation/excitation_snapshot.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <memory>
#include <vector>

namespace coilgun::simulation {

struct StageRuntimeState {
    bool triggered = false;
    bool excitation_finished = false;
    bool circuit_active = false;
    bool stage_completed = false;
    bool crowbar_on = false;
    double trigger_time = 0.0;
    double trigger_position = 0.0;
};

struct ContinuousState {
    Eigen::VectorXd currents;
    double arm_position = 0.0;
    double arm_velocity = 0.0;
    Eigen::VectorXd filament_temperatures;

    ContinuousState& operator+=(const ContinuousState& rhs);
    ContinuousState& operator*=(double scalar);
};

ContinuousState operator+(ContinuousState lhs, const ContinuousState& rhs);
ContinuousState operator*(double scalar, ContinuousState state);

using SimState = ContinuousState;
using MultiStageState = ContinuousState;

struct IntegrationState {
    ContinuousState physical;
    std::vector<std::unique_ptr<ExcitationSnapshot>> excitations;
    std::vector<StageRuntimeState> stages;
};

IntegrationState clone_integration_state(const IntegrationState& state);
void restore_integration_state(IntegrationState& destination,
                               const IntegrationState& source);
const StageRuntimeState& stage_state(const IntegrationState& state,
                                     std::size_t stage);
void mark_stage_completed(IntegrationState& state, std::size_t stage);

} // namespace coilgun::simulation
