/**
 * @file gpu_execution_config.hpp
 * @brief Host-only GPU execution configuration and static policy selection.
 */

#pragma once

#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace coilgun::simulation::cuda {

enum class BackendMode {
    Auto = 0,
    Graph = 1,
    Persistent = 2,
    Fallback = 3,
    Direct = 4,
};

enum class SolverMode {
    Auto = 0,
    Eigen = 1,
    Batched = 2,
    CuSolver = Batched,
};

enum class PrecisionMode {
    Standard = 0,
    Full = 1,
    Aggressive = 2,
};

enum class ThermalMode {
    Auto = 0,
    Disabled = 1,
    Cpu = 2,
    Gpu = 3,
};

enum class FallbackReason : int {
    None = 0,
    CapabilityUnavailable = 1,
    DeterminismRequired = 2,
    RuntimeFailure = 3,
    MetadataConflict = 4,
};

struct GpuExecutionConfig {
    BackendMode backend = BackendMode::Auto;
    SolverMode solver = SolverMode::Auto;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Auto;
    bool enable_calibration = false;
    bool deterministic = false;
    int device_id = 0;
    int threads_per_block = 512;
    // Request metadata copied to ExecutionReport. Host-wall timing categories
    // are collected independently; this flag does not enable NVTX.
    bool enable_profiling = false;
    void validate() const {
        if (device_id < 0) throw std::invalid_argument("GPU device_id must not be negative");
        if (threads_per_block <= 0) {
            throw std::invalid_argument("GPU threads_per_block must be positive");
        }
        if (threads_per_block > 512 ||
            (threads_per_block & (threads_per_block - 1)) != 0) {
            throw std::invalid_argument("GPU threads_per_block must be a positive power of two no greater than 512");
        }
    }
};

/** Static capabilities used by the planner; no CUDA headers are required. */
struct GpuCapability {
    bool supports_graph = true;
    bool supports_persistent = true;
    bool supports_batched_solver = true;
    bool supports_gpu_thermal = true;
    bool persistent_is_deterministic = false;
    // The current engine does not own the dedicated control stream required
    // by the resident persistent protocol. Keep this new field last so
    // existing aggregate initializers retain their positional meaning.
    bool supports_persistent_control_stream = false;
};

struct GpuExecutionPolicy {
    BackendMode requested_backend = BackendMode::Auto;
    SolverMode requested_solver = SolverMode::Auto;
    PrecisionMode requested_precision = PrecisionMode::Full;
    ThermalMode requested_thermal = ThermalMode::Auto;
    BackendMode backend = BackendMode::Fallback;
    SolverMode solver = SolverMode::Eigen;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Disabled;
    FallbackReason backend_fallback_reason = FallbackReason::None;
    FallbackReason solver_fallback_reason = FallbackReason::None;
    FallbackReason thermal_fallback_reason = FallbackReason::None;
};

inline constexpr bool is_deterministic_backend(BackendMode backend,
                                                const GpuCapability& capability) noexcept {
    return backend != BackendMode::Persistent || capability.persistent_is_deterministic;
}

class GpuExecutionPlanner {
public:
    /**
     * Resolve one immutable policy from problem dimensions and configuration.
     * This function is pure host code and does not inspect CUDA or timings.
     */
    static constexpr GpuExecutionPolicy plan(
        std::size_t n_stages,
        std::size_t n_filaments,
        std::size_t batch_size,
        bool thermal_enabled,
        const GpuCapability& capability,
        const GpuExecutionConfig& config) noexcept
    {
        const std::size_t dimension = n_stages + n_filaments;
        const bool large_workload = batch_size >= 8 || dimension >= 128;

        GpuExecutionPolicy result;
        result.requested_backend = config.backend;
        result.requested_solver = config.solver;
        result.requested_precision = config.precision;
        result.requested_thermal = config.thermal;
        result.precision = config.precision;

        if (config.backend == BackendMode::Persistent) {
            if (!capability.supports_persistent ||
                !capability.supports_persistent_control_stream) {
                result.backend = BackendMode::Fallback;
                result.backend_fallback_reason = FallbackReason::CapabilityUnavailable;
            } else if (config.deterministic && !capability.persistent_is_deterministic) {
                result.backend = BackendMode::Fallback;
                result.backend_fallback_reason = FallbackReason::DeterminismRequired;
            } else {
                result.backend = BackendMode::Persistent;
            }
        } else if (config.backend == BackendMode::Graph) {
            result.backend = BackendMode::Fallback;
            result.backend_fallback_reason = capability.supports_graph
                ? FallbackReason::MetadataConflict
                : FallbackReason::CapabilityUnavailable;
        } else if (config.backend == BackendMode::Fallback) {
            result.backend = BackendMode::Fallback;
        } else if (config.backend == BackendMode::Direct) {
            result.backend = BackendMode::Direct;
        } else {
            // Auto remains conservative at the host-only planning layer. The
            // CUDA engine may upgrade an explicit Graph request after device
            // capabilities are known.
            result.backend = BackendMode::Fallback;
            if (large_workload) result.backend_fallback_reason = FallbackReason::MetadataConflict;
        }

        if (config.solver == SolverMode::Batched) {
            result.solver = capability.supports_batched_solver ? SolverMode::Batched
                                                                : SolverMode::Eigen;
            if (result.solver != SolverMode::Batched) result.solver_fallback_reason = FallbackReason::CapabilityUnavailable;
        } else if (config.solver == SolverMode::Eigen) {
            result.solver = SolverMode::Eigen;
        } else {
            result.solver = (large_workload && capability.supports_batched_solver)
                ? SolverMode::Batched : SolverMode::Eigen;
            if (large_workload && result.solver != SolverMode::Batched) result.solver_fallback_reason = FallbackReason::CapabilityUnavailable;
        }

        if (!thermal_enabled || config.thermal == ThermalMode::Disabled) {
            result.thermal = ThermalMode::Disabled;
        } else if (config.thermal == ThermalMode::Cpu) {
            result.thermal = ThermalMode::Cpu;
        } else if (config.thermal == ThermalMode::Gpu) {
            result.thermal = capability.supports_gpu_thermal ? ThermalMode::Gpu
                                                              : ThermalMode::Cpu;
            if (result.thermal != ThermalMode::Gpu) result.thermal_fallback_reason = FallbackReason::CapabilityUnavailable;
        } else {
            result.thermal = (large_workload && capability.supports_gpu_thermal)
                ? ThermalMode::Gpu : ThermalMode::Cpu;
            if (large_workload && result.thermal != ThermalMode::Gpu) result.thermal_fallback_reason = FallbackReason::CapabilityUnavailable;
        }

        // An explicit Graph request with no Graph capability is a resolved
        // CPU contract. Do not retain an otherwise independent CUDA solver or
        // thermal selection under that backend.
        if (config.backend == BackendMode::Graph && !capability.supports_graph) {
            result.backend = BackendMode::Fallback;
            result.solver = SolverMode::Eigen;
            result.thermal = thermal_enabled ? ThermalMode::Cpu : ThermalMode::Disabled;
        }

        return result;
    }
};

} // namespace coilgun::simulation::cuda
