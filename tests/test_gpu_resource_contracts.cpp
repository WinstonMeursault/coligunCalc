#include <doctest/doctest.h>

#include "coilgun/coilgun_cuda.hpp"

#include <stdexcept>
#include <vector>

TEST_CASE("GpuAdaptor rejects unconfigured transfers") {
    coilgun::simulation::cuda::GpuAdaptor adaptor;
    CHECK_THROWS_AS(adaptor.upload_separation({0.0}), std::invalid_argument);
    std::vector<double> mutual;
    std::vector<double> gradient;
    CHECK_THROWS_AS(adaptor.download_results(mutual, gradient, -1),
                    std::invalid_argument);
}
