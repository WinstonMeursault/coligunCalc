/**
 * @file trigger_config.hpp
 * @brief Multi-stage trigger configuration (position or time-delay).
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::simulation {

/**
 * @brief Trigger mode for activating a driving coil stage.
 */
enum class TriggerMode {
    Position,   ///< Trigger when armature centre passes a given position (m).
    TimeDelay   ///< Trigger after a fixed delay from the previous stage (s).
};

/**
 * @brief Trigger configuration for a single driving coil stage.
 *
 * Stage 0 always triggers at t=0 and does not need a TriggerConfig.
 * trigger_configs_[i-1] applies to stage i (i >= 1).
 */
struct TriggerConfig {
    TriggerMode mode;   ///< Trigger mode.
    double value;       ///< Trigger value: position (m) or time delay (s).
};

} // namespace coilgun::simulation
