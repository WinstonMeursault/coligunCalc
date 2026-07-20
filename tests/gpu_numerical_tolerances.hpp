#pragma once

#include "coilgun/simulation/cuda/gpu_backend.hpp"

#include <algorithm>
#include <cmath>

namespace gpu_test {

struct NumericalTolerance {
    double relative;
    double absolute;
};

inline NumericalTolerance tolerance_for(coilgun::simulation::cuda::GpuOptLevel mode) {
    using coilgun::simulation::cuda::GpuOptLevel;
    switch (mode) {
    case GpuOptLevel::Standard: return {5e-6, 1e-10};
    case GpuOptLevel::Full: return {1e-4, 1e-9};
    case GpuOptLevel::Aggressive: return {1e-2, 1e-8};
    }
    return {1e-4, 1e-9};
}

inline bool numerically_equal(double actual, double expected, NumericalTolerance tolerance) {
    return std::abs(actual - expected) <= tolerance.absolute +
        tolerance.relative * std::max(std::abs(actual), std::abs(expected));
}

} // namespace gpu_test
