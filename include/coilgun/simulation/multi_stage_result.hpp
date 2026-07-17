/**
 * @file multi_stage_result.hpp
 * @brief Multi-stage simulation result data structures.
 * @author Winston Meursault
 */

#pragma once

#include <cstddef>
#include <vector>

namespace coilgun::simulation {

/**
 * @brief Snapshot of physical state shared between single-stage and
 *        multi-stage result records.
 */
struct StepSnapshot {
    double time = 0.0;                      ///< Elapsed time, s.
    double arm_position = 0.0;              ///< Armature centre position, m.
    double arm_velocity = 0.0;              ///< Armature velocity, m/s.
    double force = 0.0;                     ///< Net axial Lorentz force, N.
    std::vector<double> filament_currents;  ///< Current per filament, A.
    std::vector<double> filament_temperatures; ///< Filament temperatures, K (thermal mode only).
};

/**
 * @brief One time-step snapshot of a multi-stage simulation.
 *
 * All per-stage vectors are of size n_stages.  Entries for stages that have
 * not yet triggered are set to zero.
 */
struct MultiStageStep {
    StepSnapshot state;                      ///< Shared physical state.
    std::vector<double> cap_voltages;        ///< [n_stages] capacitor voltages, V.
    std::vector<double> coil_currents;       ///< [n_stages] coil currents, A.
};

/**
 * @brief Per-stage summary statistics.
 */
struct PerStageSummary {
    int    stage_index = 0;         ///< Zero-based stage index.
    double trigger_time = 0.0;      ///< Time at which this stage triggered, s.
    double trigger_position = 0.0;  ///< Armature position at trigger, m.
    double peak_current = 0.0;      ///< Peak coil current in this stage, A.
    double max_force = 0.0;         ///< Peak force contributed by this stage, N.
    double energy_depleted = 0.0;   ///< Capacitor energy lost, J (E_init - E_end).
    int    step_count_active = 0;   ///< Number of steps while this stage was active.
};

/**
 * @brief Aggregate statistics for a multi-stage run.
 */
struct MultiStageSummary {
    double muzzle_velocity = 0.0;           ///< Final armature velocity, m/s.
    double total_time = 0.0;                ///< Elapsed simulation time, s.
    double max_force = 0.0;                 ///< Peak axial force (any stage), N.
    double peak_coil_current = 0.0;         ///< Peak coil current (any stage), A.
    double efficiency = 0.0;                ///< E_kin / Sigma(0.5 * C_i * U0_i^2), 0..1.
    int    step_count = 0;                  ///< Total step count.
    std::vector<PerStageSummary> per_stage; ///< One entry per triggered stage.
};

/**
 * @brief Complete multi-stage simulation result.
 */
struct MultiStageResult {
    std::vector<MultiStageStep> history; ///< Step-by-step state record.
    MultiStageSummary summary;           ///< End-of-run summary.

    /**
     * @brief Down-sample the time history for export/plotting.
     * @param every_n Keep every Nth step.
     * @return A sparse copy of this result.
     */
    MultiStageResult sampled(int every_n) const;
};

} // namespace coilgun::simulation
