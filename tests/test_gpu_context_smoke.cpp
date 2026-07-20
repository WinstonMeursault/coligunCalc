#include <doctest/doctest.h>

#include "coilgun/simulation/cuda/gpu_execution_context.hpp"

#include <cuda_runtime_api.h>

TEST_CASE("GPU execution context smoke") {
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) {
        MESSAGE("CUDA device unavailable; skipping context smoke test");
        return;
    }

    int before = -1;
    REQUIRE(cudaGetDevice(&before) == cudaSuccess);
    {
        coilgun::simulation::cuda::GpuExecutionContext context({0, cudaStreamNonBlocking, 4096});
        context.record_start();
        context.record_stop();
        context.synchronize();
        CHECK(context.workspace_bytes() == 4096);
        CHECK(context.cublas() != nullptr);
        CHECK(context.cusolver() != nullptr);
        int during = -1;
        REQUIRE(cudaGetDevice(&during) == cudaSuccess);
        CHECK(during == before);
    }
    int after = -1;
    REQUIRE(cudaGetDevice(&after) == cudaSuccess);
    CHECK(after == before);
}
