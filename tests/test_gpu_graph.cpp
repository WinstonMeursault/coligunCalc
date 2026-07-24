#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_graph.hpp"

#include <cstddef>
#include <cstdint>

#if defined(COILGUN_CUDA_AVAILABLE)
#include <cuda_runtime_api.h>
#endif

using namespace coilgun::simulation::cuda;

TEST_CASE("graph variant key includes stage signature and execution modes") {
    GpuGraphVariantKey first;
    first.stage_signature = 0x11;
    first.batch_capacity = 8;
    first.layout_signature = 0x22;
    first.precision = PrecisionMode::Full;
    first.thermal = ThermalMode::Gpu;
    first.solver = SolverMode::Batched;

    auto changed = first;
    changed.stage_signature++;
    CHECK_FALSE(first == changed);
    changed = first;
    changed.thermal = ThermalMode::Disabled;
    CHECK_FALSE(first == changed);
}

TEST_CASE("graph boundary separates topology from runtime masks") {
    GpuGraphBoundaryState first;
    first.topology.stage_signature = 0x11;
    first.topology.batch_capacity = 8;
    first.topology.layout_signature = 0x22;
    first.runtime_masks.stage_mask = {1, 1, 0};
    first.runtime_masks.mutual_stage_mask = {1, 0, 0};

    auto mask_changed = first;
    mask_changed.runtime_masks.stage_mask[0] = 0;
    CHECK_FALSE(first.runtime_masks == mask_changed.runtime_masks);
    CHECK_FALSE(mask_changed.requires_rebuild_from(first));

    auto topology_changed = mask_changed;
    topology_changed.topology.layout_signature++;
    CHECK(topology_changed.requires_rebuild_from(first));
}

#if defined(COILGUN_CUDA_AVAILABLE)
TEST_CASE("graph cache reuses topology across runtime mask changes") {
    GpuGraphCache cache;
    GpuGraphBoundaryState first;
    first.topology.batch_capacity = 4;
    first.runtime_masks.stage_mask = {1, 1};
    first.runtime_masks.mutual_stage_mask = {1, 1};
    auto mask_changed = first;
    mask_changed.runtime_masks.stage_mask[1] = 0;

    int captures = 0;
    REQUIRE(cache.select_or_capture(first, [&] {
        ++captures;
        return GraphCaptureStatus::success();
    }) == GraphCaptureStatus::success());
    REQUIRE(cache.select_or_capture(mask_changed, [&] {
        ++captures;
        return GraphCaptureStatus::success();
    }) == GraphCaptureStatus::success());
    CHECK(captures == 1);
    CHECK(cache.variant_count() == 1);
    CHECK(cache.capture_count() == 1);

    auto topology_changed = mask_changed;
    topology_changed.topology.layout_signature = 9;
    REQUIRE(cache.select_or_capture(topology_changed, [&] {
        ++captures;
        return GraphCaptureStatus::success();
    }) == GraphCaptureStatus::success());
    CHECK(captures == 2);
    CHECK(cache.variant_count() == 2);
    CHECK(cache.capture_count() == 2);
}

TEST_CASE("graph cache reuses a captured variant at step boundaries") {
    GpuGraphCache cache;
    GpuGraphVariantKey key;
    key.stage_signature = 7;
    key.batch_capacity = 4;
    key.layout_signature = 9;

    int captures = 0;
    int replays = 0;
    REQUIRE(cache.select_or_capture(key, [&] {
        ++captures;
        return GraphCaptureStatus::success();
    }) == GraphCaptureStatus::success());
    CHECK(captures == 1);

    REQUIRE(cache.replay([&] {
        ++replays;
        return GraphCaptureStatus::success();
    }) == GraphCaptureStatus::success());
    REQUIRE(cache.select_or_capture(key, [&] {
        ++captures;
        return GraphCaptureStatus::success();
    }) == GraphCaptureStatus::success());
    CHECK(captures == 1);
    CHECK(replays == 1);
    CHECK(cache.capture_count() == 1);
    CHECK(cache.replay_count() == 1);
}

TEST_CASE("graph cache selects a new variant only at the next boundary") {
    GpuGraphCache cache;
    GpuGraphVariantKey first;
    GpuGraphVariantKey second = first;
    second.stage_signature = 2;

    REQUIRE(cache.select_or_capture(first, [] { return GraphCaptureStatus::success(); })
            == GraphCaptureStatus::success());
    CHECK(cache.current_key() == first);
    REQUIRE(cache.select_or_capture(second, [] { return GraphCaptureStatus::success(); })
            == GraphCaptureStatus::success());
    CHECK(cache.current_key() == second);
    CHECK(cache.variant_count() == 2);
    CHECK(cache.capture_count() == 2);
}

TEST_CASE("graph capture failure is structured and locks fallback") {
    GpuGraphCache cache;
    GpuGraphVariantKey key;

    const auto result = cache.select_or_capture(key, [] {
        return GraphCaptureStatus::failed(GraphCapturePhase::BeginCapture, 719,
                                           "capture rejected");
    });

    CHECK_FALSE(result.ok);
    CHECK(result.failure.phase == GraphCapturePhase::BeginCapture);
    CHECK(result.failure.cuda_error == 719);
    CHECK(result.failure.message == "capture rejected");
    CHECK(result.failure.key == key);
    CHECK(result.failure.fallback_locked);
    CHECK(cache.fallback_locked());
    CHECK(cache.variant_count() == 0);
}

TEST_CASE("graph cache preserves workspace identity across replays") {
    GpuGraphCache cache;
    GpuGraphVariantKey key;
    const std::uintptr_t workspace = 0x1234;

    REQUIRE(cache.select_or_capture(key, [&] {
        return GraphCaptureStatus::success(workspace, 4096);
    }) == GraphCaptureStatus::success(workspace, 4096));
    REQUIRE(cache.replay([] { return GraphCaptureStatus::success(); })
            == GraphCaptureStatus::success());

    CHECK(cache.current_workspace().pointer == workspace);
    CHECK(cache.current_workspace().bytes == 4096);
}
#endif

#if defined(COILGUN_CUDA_AVAILABLE)
TEST_CASE("CUDA graph capture instantiates once and replays repeatedly") {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        MESSAGE("CUDA device unavailable; skipping graph capture test");
        return;
    }

    cudaStream_t stream = nullptr;
    REQUIRE(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) == cudaSuccess);
    int value = 0;
    int* device_value = nullptr;
    REQUIRE(cudaMalloc(reinterpret_cast<void**>(&device_value), sizeof(int)) == cudaSuccess);

    GpuGraphCache cache;
    GpuGraphBoundaryState boundary;
    boundary.topology.batch_capacity = 1;
    boundary.runtime_masks.stage_mask = {1};
    boundary.runtime_masks.mutual_stage_mask = {1};
    const auto captured = cache.capture_and_select(
        boundary, stream,
        [&](cudaStream_t captured_stream) {
            return cudaMemsetAsync(device_value, 0, sizeof(int), captured_stream) == cudaSuccess
                ? GraphCaptureStatus::success()
                : GraphCaptureStatus::failed(GraphCapturePhase::CaptureBody,
                                              static_cast<int>(cudaGetLastError()),
                                              "cudaMemsetAsync failed");
        });
    REQUIRE(captured.ok);
    REQUIRE(cache.capture_count() == 1);
    REQUIRE(cache.replay(stream).ok);
    REQUIRE(cache.replay(stream).ok);
    REQUIRE(cudaStreamSynchronize(stream) == cudaSuccess);
    REQUIRE(cudaMemcpy(&value, device_value, sizeof(value), cudaMemcpyDeviceToHost) == cudaSuccess);

    CHECK(value == 0);
    CHECK(cache.capture_count() == 1);
    CHECK(cache.replay_count() == 2);

    cudaFree(device_value);
    cudaStreamDestroy(stream);
}
#endif
