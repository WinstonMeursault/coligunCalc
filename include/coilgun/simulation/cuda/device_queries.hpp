#pragma once

#include <cuda_runtime_api.h>

#include <atomic>
#include <cstddef>
#include <cstdint>

namespace coilgun::simulation::cuda::detail {

// cudaGetDeviceProperties is a heavyweight driver query and the grid limits
// never change for a given device, so cache them per device id instead of
// querying on every kernel launch. Returns false when the current device
// cannot be queried.
inline bool device_max_grid(std::size_t (&limits)[3]) {
    constexpr int kCachedDevices = 64;
    struct CachedLimits {
        std::atomic<bool> valid;
        std::atomic<std::size_t> x;
        std::atomic<std::size_t> y;
        std::atomic<std::size_t> z;
    };
    static CachedLimits cache[kCachedDevices];
    int device = 0;
    if (cudaGetDevice(&device) != cudaSuccess || device < 0) return false;
    if (device < kCachedDevices) {
        const CachedLimits& entry = cache[device];
        if (entry.valid.load(std::memory_order_relaxed)) {
            limits[0] = entry.x.load(std::memory_order_relaxed);
            limits[1] = entry.y.load(std::memory_order_relaxed);
            limits[2] = entry.z.load(std::memory_order_relaxed);
            return true;
        }
    }
    cudaDeviceProp properties{};
    if (cudaGetDeviceProperties(&properties, device) != cudaSuccess) return false;
    limits[0] = static_cast<std::size_t>(properties.maxGridSize[0]);
    limits[1] = static_cast<std::size_t>(properties.maxGridSize[1]);
    limits[2] = static_cast<std::size_t>(properties.maxGridSize[2]);
    if (device < kCachedDevices) {
        CachedLimits& entry = cache[device];
        entry.x.store(limits[0], std::memory_order_relaxed);
        entry.y.store(limits[1], std::memory_order_relaxed);
        entry.z.store(limits[2], std::memory_order_relaxed);
        entry.valid.store(true, std::memory_order_relaxed);
    }
    return true;
}

// The buffers handed to the per-step launchers are stable cudaMalloc
// allocations, so positive device-pointer verdicts are cached to avoid a
// driver query per checked pointer per step. Negative verdicts are never
// cached: a pointer may become device memory at any time.
inline bool is_device_pointer(const void* pointer) {
    if (!pointer) return false;
    constexpr std::size_t kCacheSlots = 1024;
    static std::atomic<const void*> verified[kCacheSlots];
    const auto slot =
        (reinterpret_cast<std::uintptr_t>(pointer) >> 6) & (kCacheSlots - 1);
    if (verified[slot].load(std::memory_order_relaxed) == pointer) return true;
    cudaPointerAttributes attributes{};
    const auto status = cudaPointerGetAttributes(&attributes, pointer);
#if CUDART_VERSION >= 10000
    const bool result = status == cudaSuccess &&
        (attributes.type == cudaMemoryTypeDevice ||
         attributes.type == cudaMemoryTypeManaged);
#else
    const bool result =
        status == cudaSuccess && attributes.memoryType == cudaMemoryTypeDevice;
#endif
    if (result) verified[slot].store(pointer, std::memory_order_relaxed);
    return result;
}

} // namespace coilgun::simulation::cuda::detail
