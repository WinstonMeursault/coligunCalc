#pragma once

#include "coilgun/simulation/cuda/gpu_engine.hpp"

namespace gpu_test {

inline coilgun::simulation::cuda::GpuGeometryInput geometry(std::size_t stages = 1,
                                                              std::size_t filaments = 1,
                                                              bool thermal = false) {
    using namespace coilgun::simulation::cuda;
    GpuGeometryInput value;
    value.n_stages = stages;
    value.n_filaments = filaments;
    value.thermal_enabled = thermal;
    value.stage_geometry.resize(stages);
    value.stage_inner_radii.resize(stages);
    value.stage_outer_radii.resize(stages);
    value.stage_lengths.resize(stages);
    value.stage_turns.resize(stages, 120);
    value.stage_positions.resize(stages);
    for (std::size_t i = 0; i < stages; ++i) {
        value.stage_geometry[i] = 0.05 + 0.02 * i;
        value.stage_inner_radii[i] = value.stage_geometry[i] - 0.005;
        value.stage_outer_radii[i] = value.stage_geometry[i] + 0.005;
        value.stage_lengths[i] = 0.03;
        value.stage_positions[i] = 0.02 * i;
    }
    value.filament_geometry.resize(filaments);
    value.filament_inner_radii.resize(filaments);
    value.filament_outer_radii.resize(filaments);
    value.filament_lengths.resize(filaments, 0.02);
    value.filament_positions.resize(filaments);
    for (std::size_t i = 0; i < filaments; ++i) {
        value.filament_geometry[i] = 0.01 + 0.002 * i;
        value.filament_inner_radii[i] = value.filament_geometry[i] - 0.001;
        value.filament_outer_radii[i] = value.filament_geometry[i] + 0.001;
        value.filament_positions[i] = 0.04 + 0.002 * i;
    }
    return value;
}

inline coilgun::simulation::cuda::GpuEngineState state(std::size_t batch,
                                                         std::size_t stages,
                                                         std::size_t filaments,
                                                         bool thermal = false) {
    using namespace coilgun::simulation::cuda;
    GpuEngineState value;
    const auto d = stages + filaments;
    value.currents.assign(batch * d, 0.0);
    value.m1.assign(batch * stages * filaments, 0.0);
    value.dm1.assign(batch * stages * filaments, 0.0);
    value.active_mask.assign(batch, 1);
    value.trigger_mask.assign(batch * stages, 1);
    value.velocity.assign(batch, 0.0);
    value.position.assign(batch, 0.0);
    value.dt = 2.5e-5;
    value.mass = 0.02;
    if (thermal) {
        value.temperatures.assign(batch * filaments, 293.0);
        value.filament_masses.assign(batch * filaments, 0.01);
        value.reference_resistances.assign(batch * filaments, 0.003);
        value.filament_materials.assign(batch * filaments, 0);
        value.material_density = 2700.0;
    }
    return value;
}

} // namespace gpu_test
