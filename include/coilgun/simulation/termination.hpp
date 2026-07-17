/**
 * @file termination.hpp
 * @brief Simulation termination criteria.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::simulation {

/**
 * @brief Configurable termination policy for time-marching simulations.
 *
 * Supports multiple criteria: maximum step count, velocity decay detection
 * (armature has stopped accelerating), and hard position bound (barrel exit).
 */
struct TerminationPolicy {
    int    max_steps             = 20000;  ///< Absolute upper bound on step count.
    int    velocity_decay_steps  = 5;      ///< Consecutive velocity drops needed to trigger decay.
    double accel_threshold       = 0.1;    ///< Acceleration below this is treated as zero, m/s^2.
    bool   enable_velocity_check = true;   ///< Enable velocity-decay termination.
    bool   enable_bound_check    = false;  ///< Enable barrel-end position termination.
    double barrel_end_position   = 1.0;    ///< Barrel exit coordinate, m.

    /// @brief Factory returning the default policy.
    static TerminationPolicy defaults();
};

} // namespace coilgun::simulation
