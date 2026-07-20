#include "coilgun/simulation/cuda/gpu_graph.hpp"

#include <stdexcept>
#include <utility>

namespace coilgun::simulation::cuda {

GpuGraphCache::~GpuGraphCache() {
#if defined(COILGUN_CUDA_AVAILABLE)
    for (auto& item : variants_) {
        if (item.second.exec != nullptr) cudaGraphExecDestroy(item.second.exec);
    }
#endif
}

GraphCaptureStatus GpuGraphCache::lock_failure(const GpuGraphVariantKey& key,
                                               GraphCaptureStatus status) {
    status.failure.key = key;
    status.failure.fallback_locked = true;
    fallback_locked_ = true;
    current_ = nullptr;
    return status;
}

GraphCaptureStatus GpuGraphCache::select_or_capture(const GpuGraphVariantKey& key,
                                                     const CaptureFn& capture) {
    if (fallback_locked_) {
        return lock_failure(key, GraphCaptureStatus::failed(
            GraphCapturePhase::Instantiate, 0,
            "CUDA graph cache is locked to fallback after capture failure"));
    }

    const auto found = variants_.find(key);
    if (found != variants_.end()) {
        current_ = &found->second;
        return GraphCaptureStatus::success(current_->workspace.pointer,
                                            current_->workspace.bytes);
    }

    GraphCaptureStatus status = capture();
    if (!status.ok) return lock_failure(key, std::move(status));

    Entry entry{key, status.workspace};
    const auto inserted = variants_.emplace(key, std::move(entry));
    current_ = &inserted.first->second;
    ++capture_count_;
    return status;
}

GraphCaptureStatus GpuGraphCache::replay(const ReplayFn& replay_fn) {
    if (fallback_locked_) {
        return GraphCaptureStatus::failed(GraphCapturePhase::Replay, 0,
                                           "CUDA graph cache is locked to fallback");
    }
    if (current_ == nullptr) {
        return GraphCaptureStatus::failed(GraphCapturePhase::Replay, 0,
                                           "no graph variant selected at step boundary");
    }

    GraphCaptureStatus status = replay_fn();
    if (!status.ok) {
        status.failure.key = current_->key;
        status.failure.fallback_locked = true;
        fallback_locked_ = true;
        return status;
    }
    ++replay_count_;
    return status;
}

#if defined(COILGUN_CUDA_AVAILABLE)
GraphCaptureStatus GpuGraphCache::capture_and_select(const GpuGraphVariantKey& key,
                                                      cudaStream_t stream,
                                                      const CudaCaptureBody& body,
                                                      GraphWorkspace workspace) {
    if (fallback_locked_) {
        return lock_failure(key, GraphCaptureStatus::failed(
            GraphCapturePhase::BeginCapture, cudaErrorUnknown,
            "CUDA graph cache is locked to fallback after capture failure"));
    }

    const auto found = variants_.find(key);
    if (found != variants_.end()) {
        current_ = &found->second;
        return GraphCaptureStatus::success(current_->workspace.pointer,
                                            current_->workspace.bytes);
    }

    cudaError_t error = cudaStreamBeginCapture(stream, cudaStreamCaptureModeThreadLocal);
    if (error != cudaSuccess) {
        return lock_failure(key, GraphCaptureStatus::failed(
            GraphCapturePhase::BeginCapture, static_cast<int>(error),
            cudaGetErrorString(error)));
    }

    GraphCaptureStatus body_status = body(stream);
    cudaGraph_t graph = nullptr;
    error = cudaStreamEndCapture(stream, &graph);
    if (!body_status.ok) {
        if (graph != nullptr) cudaGraphDestroy(graph);
        return lock_failure(key, std::move(body_status));
    }
    if (error != cudaSuccess) {
        return lock_failure(key, GraphCaptureStatus::failed(
            GraphCapturePhase::EndCapture, static_cast<int>(error),
            cudaGetErrorString(error)));
    }

    cudaGraphExec_t exec = nullptr;
    error = cudaGraphInstantiate(&exec, graph, nullptr, nullptr, 0);
    cudaGraphDestroy(graph);
    if (error != cudaSuccess) {
        return lock_failure(key, GraphCaptureStatus::failed(
            GraphCapturePhase::Instantiate, static_cast<int>(error),
            cudaGetErrorString(error)));
    }

    Entry entry{key, workspace, exec};
    const auto inserted = variants_.emplace(key, std::move(entry));
    current_ = &inserted.first->second;
    ++capture_count_;
    return GraphCaptureStatus::success(workspace.pointer, workspace.bytes);
}

GraphCaptureStatus GpuGraphCache::replay(cudaStream_t stream) {
    if (fallback_locked_) {
        return GraphCaptureStatus::failed(GraphCapturePhase::Replay,
                                          cudaErrorUnknown,
                                          "CUDA graph cache is locked to fallback");
    }
    if (current_ == nullptr || current_->exec == nullptr) {
        return GraphCaptureStatus::failed(GraphCapturePhase::Replay,
                                          cudaErrorInvalidResourceHandle,
                                          "no CUDA graph variant selected at step boundary");
    }
    const cudaError_t error = cudaGraphLaunch(current_->exec, stream);
    if (error != cudaSuccess) {
        GraphCaptureStatus status = GraphCaptureStatus::failed(
            GraphCapturePhase::Replay, static_cast<int>(error),
            cudaGetErrorString(error));
        status.failure.key = current_->key;
        status.failure.fallback_locked = true;
        fallback_locked_ = true;
        return status;
    }
    ++replay_count_;
    return GraphCaptureStatus::success(current_->workspace.pointer,
                                        current_->workspace.bytes);
}
#endif

const GpuGraphVariantKey& GpuGraphCache::current_key() const {
    if (current_ == nullptr) throw std::logic_error("no graph variant selected");
    return current_->key;
}

GraphWorkspace GpuGraphCache::current_workspace() const noexcept {
    return current_ == nullptr ? GraphWorkspace{} : current_->workspace;
}

} // namespace coilgun::simulation::cuda
