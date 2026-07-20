#pragma once

#include "coilgun/simulation/cuda/gpu_adaptor.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include <cstddef>
#include <cstdint>

namespace coilgun::simulation::cuda {

struct MutualPipelineView {
    const CoilGeo* coils;
    const FilGeo* filaments;
    const double* separations;
    const std::uint8_t* active_mask; // [B][S][F], zero disables a physical slot
    double* mutual;
    double* gradient;
    std::size_t batch_size;
    std::size_t stage_count;
    std::size_t filament_count;
    int n_nodes = 9;
};

void launch_mutual_pipeline(const MutualPipelineView& view,
                            GpuOptLevel opt_level = GpuOptLevel::Standard,
                            int threads_per_block = 256,
                            cudaStream_t stream = nullptr);

std::size_t mutual_pipeline_index(std::size_t simulation,
                                  std::size_t stage,
                                  std::size_t filament,
                                  std::size_t stage_count,
                                  std::size_t filament_count);

} // namespace coilgun::simulation::cuda
