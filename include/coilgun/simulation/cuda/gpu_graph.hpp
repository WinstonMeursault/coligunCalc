/**
 * @file gpu_graph.hpp
 * @brief Step-boundary CUDA Graph variant selection and replay cache.
 */

#pragma once

#include "coilgun/simulation/cuda/gpu_execution_config.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#if defined(COILGUN_CUDA_AVAILABLE)
#include <cuda_runtime_api.h>
#endif

namespace coilgun::simulation::cuda {

struct GpuGraphVariantKey {
    std::uint64_t stage_signature = 0;
    std::size_t batch_capacity = 0;
    std::uint64_t layout_signature = 0;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Disabled;
    SolverMode solver = SolverMode::Eigen;

    bool operator==(const GpuGraphVariantKey& other) const noexcept {
        return stage_signature == other.stage_signature &&
               batch_capacity == other.batch_capacity &&
               layout_signature == other.layout_signature &&
               precision == other.precision && thermal == other.thermal &&
               solver == other.solver;
    }
};

struct GpuGraphVariantKeyHash {
    std::size_t operator()(const GpuGraphVariantKey& key) const noexcept {
        std::size_t hash = static_cast<std::size_t>(key.stage_signature);
        hash = hash * 31u + key.batch_capacity;
        hash = hash * 31u + static_cast<std::size_t>(key.layout_signature);
        hash = hash * 31u + static_cast<std::size_t>(key.precision);
        hash = hash * 31u + static_cast<std::size_t>(key.thermal);
        return hash * 31u + static_cast<std::size_t>(key.solver);
    }
};

enum class GraphCapturePhase {
    None,
    BeginCapture,
    CaptureBody,
    EndCapture,
    Instantiate,
    Replay,
};

struct GraphCaptureFailure {
    GraphCapturePhase phase = GraphCapturePhase::None;
    int cuda_error = 0;
    std::string message;
    GpuGraphVariantKey key{};
    bool fallback_locked = false;
};

struct GraphWorkspace {
    std::uintptr_t pointer = 0;
    std::size_t bytes = 0;
};

struct GraphCaptureStatus {
    bool ok = true;
    GraphCaptureFailure failure{};
    GraphWorkspace workspace{};

    static GraphCaptureStatus success(std::uintptr_t pointer = 0,
                                      std::size_t bytes = 0) {
        GraphCaptureStatus result;
        result.workspace = {pointer, bytes};
        return result;
    }

    static GraphCaptureStatus failed(GraphCapturePhase phase, int cuda_error,
                                     std::string message) {
        GraphCaptureStatus result;
        result.ok = false;
        result.failure.phase = phase;
        result.failure.cuda_error = cuda_error;
        result.failure.message = std::move(message);
        return result;
    }

    bool operator==(const GraphCaptureStatus& other) const noexcept {
        return ok == other.ok && workspace.pointer == other.workspace.pointer &&
               workspace.bytes == other.workspace.bytes;
    }
};

/**
 * Cache for immutable graph variants. Selection and capture are intentionally
 * explicit: callers invoke select_or_capture only after a completed step.
 */
class GpuGraphCache {
public:
    using CaptureFn = std::function<GraphCaptureStatus()>;
    using ReplayFn = std::function<GraphCaptureStatus()>;

    ~GpuGraphCache();
    GpuGraphCache() = default;
    GpuGraphCache(const GpuGraphCache&) = delete;
    GpuGraphCache& operator=(const GpuGraphCache&) = delete;
    GpuGraphCache(GpuGraphCache&&) = delete;
    GpuGraphCache& operator=(GpuGraphCache&&) = delete;

    GraphCaptureStatus select_or_capture(const GpuGraphVariantKey& key,
                                          const CaptureFn& capture);
    GraphCaptureStatus replay(const ReplayFn& replay);
#if defined(COILGUN_CUDA_AVAILABLE)
    using CudaCaptureBody = std::function<GraphCaptureStatus(cudaStream_t)>;
    GraphCaptureStatus capture_and_select(const GpuGraphVariantKey& key,
                                          cudaStream_t stream,
                                          const CudaCaptureBody& body,
                                          GraphWorkspace workspace = {});
    GraphCaptureStatus replay(cudaStream_t stream);
#endif

    bool has_current() const noexcept { return current_ != nullptr; }
    bool fallback_locked() const noexcept { return fallback_locked_; }
    std::size_t variant_count() const noexcept { return variants_.size(); }
    std::size_t capture_count() const noexcept { return capture_count_; }
    std::size_t replay_count() const noexcept { return replay_count_; }
    const GpuGraphVariantKey& current_key() const;
    GraphWorkspace current_workspace() const noexcept;

private:
    struct Entry {
        GpuGraphVariantKey key;
        GraphWorkspace workspace;
#if defined(COILGUN_CUDA_AVAILABLE)
        cudaGraphExec_t exec = nullptr;
#endif
    };

    std::unordered_map<GpuGraphVariantKey, Entry, GpuGraphVariantKeyHash> variants_;
    Entry* current_ = nullptr;
    bool fallback_locked_ = false;
    std::size_t capture_count_ = 0;
    std::size_t replay_count_ = 0;

    GraphCaptureStatus lock_failure(const GpuGraphVariantKey& key,
                                    GraphCaptureStatus status);
};

} // namespace coilgun::simulation::cuda
