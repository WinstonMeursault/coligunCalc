/**
 * @file trigger_config.hpp
 * @brief Multi-stage trigger configuration (position or time-delay).
 * @author Winston Meursault
 */

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>

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

/**
 * @brief Validate one stage trigger configuration.
 *
 * NaN is always malformed.  Positive infinity is an explicit terminal
 * trigger: it never fires, so the simulator may terminate once no other
 * stage remains eligible. For Position it represents an unreachable position;
 * for TimeDelay it represents an unreachable delay. Position triggers accept
 * any other finite value; time-delay triggers require a non-negative delay.
 */
inline void validate_trigger_config(const TriggerConfig& config) {
    switch (config.mode) {
    case TriggerMode::Position:
        if (std::isnan(config.value))
            throw std::invalid_argument("position trigger value must not be NaN");
        if (config.value == -std::numeric_limits<double>::infinity())
            throw std::invalid_argument("position trigger value must not be negative infinity");
        break;
    case TriggerMode::TimeDelay:
        if (std::isnan(config.value) || config.value < 0.0)
            throw std::invalid_argument("time-delay trigger value must be non-negative and not NaN");
        if (config.value == -std::numeric_limits<double>::infinity())
            throw std::invalid_argument("time-delay trigger value must not be negative infinity");
        break;
    default:
        throw std::invalid_argument("trigger mode is invalid");
    }
}

} // namespace coilgun::simulation
