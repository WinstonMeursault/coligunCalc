/**
 * @file coilgun_cuda.hpp
 * @brief Umbrella header for GPU-accelerated coilgun simulation.
 * @author Winston Meursault
 *
 * Include this single header to access all GPU-accelerated simulation
 * classes. Also transitively includes coilgun.hpp for the full CPU API.
 */

#pragma once

#include "coilgun/coilgun.hpp"
#include "coilgun/simulation/cuda/gpu_backend.hpp"
#include "coilgun/simulation/cuda/gpu_single_stage_sim.hpp"
#include "coilgun/simulation/cuda/gpu_multi_stage_sim.hpp"
#include "coilgun/simulation/cuda/sim_batch.hpp"
