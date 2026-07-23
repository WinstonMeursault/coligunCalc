#include "coilgun/simulation/integration_state.hpp"
#include "coilgun/simulation/derivative_workspace.hpp"

#include <stdexcept>

namespace coilgun::simulation {

ContinuousState& ContinuousState::operator+=(const ContinuousState& rhs) {
    currents += rhs.currents;
    arm_position += rhs.arm_position;
    arm_velocity += rhs.arm_velocity;
    if (filament_temperatures.size() != 0 && rhs.filament_temperatures.size() != 0)
        filament_temperatures += rhs.filament_temperatures;
    return *this;
}

ContinuousState& ContinuousState::operator*=(double scalar) {
    currents *= scalar;
    arm_position *= scalar;
    arm_velocity *= scalar;
    filament_temperatures *= scalar;
    return *this;
}

ContinuousState operator+(ContinuousState lhs, const ContinuousState& rhs) {
    lhs += rhs;
    return lhs;
}

ContinuousState operator*(double scalar, ContinuousState state) {
    state *= scalar;
    return state;
}

IntegrationState clone_integration_state(const IntegrationState& state) {
    IntegrationState result;
    result.physical = state.physical;
    result.stages = state.stages;
    result.excitations.reserve(state.excitations.size());
    for (const auto& excitation : state.excitations) {
        if (!excitation) throw std::invalid_argument("excitation snapshot must not be null");
        result.excitations.push_back(excitation->clone());
    }
    return result;
}

void restore_integration_state(IntegrationState& destination,
                               const IntegrationState& source) {
    destination = clone_integration_state(source);
}

const StageRuntimeState& stage_state(const IntegrationState& state,
                                     std::size_t stage) {
    if (stage >= state.stages.size()) throw std::out_of_range("stage index");
    return state.stages[stage];
}

void mark_stage_completed(IntegrationState& state, std::size_t stage) {
    if (stage >= state.stages.size() || stage >= static_cast<std::size_t>(state.physical.currents.size()))
        throw std::out_of_range("stage index");
    state.physical.currents(static_cast<Eigen::Index>(stage)) = 0.0;
    auto& runtime = state.stages[stage];
    runtime.excitation_finished = true;
    runtime.circuit_active = false;
    runtime.stage_completed = true;
}

void DerivativeWorkspace::resize(std::size_t stages, std::size_t filaments) {
    const auto pairs = stages * filaments;
    const auto dimension = stages + filaments;
    mutual.resize(static_cast<Eigen::Index>(pairs));
    mutual_gradient.resize(static_cast<Eigen::Index>(pairs));
    system_matrix.resize(static_cast<Eigen::Index>(dimension),
                         static_cast<Eigen::Index>(dimension));
    rhs.resize(static_cast<Eigen::Index>(dimension));
    resistance.resize(static_cast<Eigen::Index>(filaments));
    active_stages.clear();
    active_stages.reserve(stages);
}

} // namespace coilgun::simulation
