/**
 * @file gpu_execution_report.hpp
 * @brief Resolved GPU execution choices and aggregate measurements.
 */

#pragma once

#include <algorithm>
#include <ostream>
#include <string>

#include "coilgun/simulation/cuda/gpu_execution_config.hpp"

namespace coilgun::simulation::cuda {

// The configuration enums are the canonical vocabulary. These aliases keep
// report consumers source-compatible while avoiding duplicate enum types.
using Backend = BackendMode;
using Solver = SolverMode;

inline const char* to_string(Backend value) noexcept {
    switch (value) {
    case Backend::Auto:       return "auto";
    case Backend::Graph:      return "graph";
    case Backend::Persistent: return "persistent";
    case Backend::Fallback:   return "fallback";
    case Backend::Direct:     return "direct";
    }
    return "unknown";
}

inline const char* to_string(Solver value) noexcept {
    switch (value) {
    case Solver::Auto:     return "auto";
    case Solver::Eigen:    return "eigen";
    case Solver::Batched: return "cusolver";
    }
    return "unknown";
}

inline const char* to_string(PrecisionMode value) noexcept {
    switch (value) {
    case PrecisionMode::Standard:   return "standard";
    case PrecisionMode::Full:       return "full";
    case PrecisionMode::Aggressive: return "aggressive";
    }
    return "unknown";
}

inline const char* to_string(ThermalMode value) noexcept {
    switch (value) {
    case ThermalMode::Auto:     return "auto";
    case ThermalMode::Cpu:      return "cpu";
    case ThermalMode::Gpu:      return "gpu";
    case ThermalMode::Disabled: return "disabled";
    }
    return "unknown";
}

inline std::ostream& operator<<(std::ostream& stream, Backend value) {
    return stream << to_string(value);
}

inline std::ostream& operator<<(std::ostream& stream, Solver value) {
    return stream << to_string(value);
}

inline std::ostream& operator<<(std::ostream& stream, PrecisionMode value) {
    return stream << to_string(value);
}

inline std::ostream& operator<<(std::ostream& stream, ThermalMode value) {
    return stream << to_string(value);
}

struct ExecutionReport {
    Backend requested_backend = Backend::Auto;
    Solver requested_solver = Solver::Auto;
    PrecisionMode requested_precision = PrecisionMode::Full;
    ThermalMode requested_thermal = ThermalMode::Auto;
    Backend backend = Backend::Auto;
    Solver solver = Solver::Auto;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Auto;

    bool calibrated = false;
    bool precision_fallback = false;
    int graph_rebuild_count = 0;
    int fallback_count = 0;

    // Cumulative host wall time, in milliseconds, for each successful
    // CUDA-backed physical pipeline. This includes host orchestration and
    // host/device transfers; it is not device-only kernel time.
    double gpu_time_ms = 0.0;
    // Cumulative host wall time for solver and thermal sections on whichever
    // CPU or CUDA-backed path executed.
    double solver_time_ms = 0.0;
    double thermal_time_ms = 0.0;
    // Cumulative host wall time spent in synchronous host/device copies.
    double transfer_time_ms = 0.0;
    double max_condition_estimate = 0.0;
    std::string fallback_reason;
    FallbackReason static_fallback_reason = FallbackReason::None;
    FallbackReason runtime_fallback_reason = FallbackReason::None;
    bool metadata_conflict = false;
    // True after at least one physical CUDA pipeline completes successfully.
    // This is cumulative audit data and is retained across reset().
    bool gpu_executed = false;
    int device_id = 0;
    int threads_per_block = 512;
    // Request metadata only. Timing fields are collected independently, and
    // this field does not imply NVTX annotations.
    bool profiling_enabled = false;

    void merge(const ExecutionReport& other) {
        if (requested_backend != Backend::Auto && other.requested_backend != Backend::Auto &&
            requested_backend != other.requested_backend) metadata_conflict = true;
        if (requested_solver != Solver::Auto && other.requested_solver != Solver::Auto &&
            requested_solver != other.requested_solver) metadata_conflict = true;
        if (requested_precision != other.requested_precision) metadata_conflict = true;
        if (requested_thermal != ThermalMode::Auto && other.requested_thermal != ThermalMode::Auto &&
            requested_thermal != other.requested_thermal) metadata_conflict = true;
        if (backend != Backend::Auto && other.backend != Backend::Auto && backend != other.backend) metadata_conflict = true;
        if (solver != Solver::Auto && other.solver != Solver::Auto && solver != other.solver) metadata_conflict = true;
        if (precision != other.precision) metadata_conflict = true;
        if (thermal != ThermalMode::Auto && other.thermal != ThermalMode::Auto &&
            thermal != other.thermal) metadata_conflict = true;
        if (profiling_enabled != other.profiling_enabled) metadata_conflict = true;
        if (!fallback_reason.empty() && !other.fallback_reason.empty() && fallback_reason != other.fallback_reason) metadata_conflict = true;
        if (device_id != other.device_id) metadata_conflict = true;
        if (threads_per_block != other.threads_per_block) metadata_conflict = true;
        calibrated = calibrated || other.calibrated;
        gpu_executed = gpu_executed || other.gpu_executed;
        precision_fallback = precision_fallback || other.precision_fallback;
        graph_rebuild_count += other.graph_rebuild_count;
        fallback_count += other.fallback_count;
        gpu_time_ms += other.gpu_time_ms;
        solver_time_ms += other.solver_time_ms;
        thermal_time_ms += other.thermal_time_ms;
        transfer_time_ms += other.transfer_time_ms;
        max_condition_estimate = std::max(max_condition_estimate,
                                          other.max_condition_estimate);
        if (fallback_reason.empty() && !other.fallback_reason.empty()) {
            fallback_reason = other.fallback_reason;
        }
        if (static_fallback_reason == FallbackReason::None) static_fallback_reason = other.static_fallback_reason;
        if (runtime_fallback_reason == FallbackReason::None) runtime_fallback_reason = other.runtime_fallback_reason;
        metadata_conflict = metadata_conflict || other.metadata_conflict;
        profiling_enabled = profiling_enabled || other.profiling_enabled;
    }
};

} // namespace coilgun::simulation::cuda
