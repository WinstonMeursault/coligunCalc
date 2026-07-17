/**
 * @file sim_result.hpp
 * @brief Simulation result data structures (step record, summary, result set).
 * @author Winston Meursault
 */

#pragma once

#include <cstddef>
#include <vector>

namespace coilgun::simulation {

/**
 * @brief One time-step snapshot of the simulation state.
 */
struct SimStep {
    double time = 0.0;                      ///< Elapsed time, s.
    double cap_voltage = 0.0;               ///< Capacitor voltage at this step, V.
    double coil_current = 0.0;              ///< Driving coil current, A.
    std::vector<double> filament_currents;  ///< Current in each armature filament, A.
    double arm_position = 0.0;              ///< Armature centre position, m.
    double arm_velocity = 0.0;              ///< Armature velocity, m/s.
    double force = 0.0;                     ///< Net axial Lorentz force, N.
    std::vector<double> filament_temperatures; ///< Filament temperatures, K (thermal mode only).
};

/**
 * @brief Aggregate statistics computed at end-of-run.
 */
struct SimSummary {
    double muzzle_velocity   = 0.0;  ///< Final armature velocity, m/s.
    double total_time        = 0.0;  ///< Elapsed simulation time, s.
    double max_force         = 0.0;  ///< Peak axial force observed, N.
    double peak_coil_current = 0.0;  ///< Peak driving coil current, A.
    double efficiency        = 0.0;  ///< Kinetic energy / initial stored energy, 0..1.
    int    step_count        = 0;    ///< Total number of steps taken.
};

/**
 * @brief Complete simulation result: time history + summary statistics.
 */
struct SimResult {
    std::vector<SimStep> history;  ///< Step-by-step state record.
    SimSummary summary;            ///< End-of-run aggregate statistics.

    /**
     * @brief Down-sample the time history for export/plotting.
     * @param every_n Keep every Nth step.
     * @return A sparse copy of this result.
     */
    SimResult sampled(int every_n) const;
};

} // namespace coilgun::simulation
