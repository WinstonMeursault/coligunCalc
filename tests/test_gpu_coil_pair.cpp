/**
 * @file test_gpu_coil_pair.cpp
 * @brief Unit tests: GPU 4D integration vs CPU for single (coil, filament) pair.
 * @author Winston Meursault
 */

#include <doctest/doctest.h>
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/quadrature.hpp"
#include <cmath>
#include <cuda_runtime.h>
#include <vector>

namespace coilgun::physics {
    void upload_gl_nodes(int n_nodes);
    extern __global__ void mutual_inductance_coil_pair_kernel(
        double rai, double rae, double la, int na,
        double rbi, double rbe, double lb, int nb,
        double separation, int n_nodes,
        double* out_M, double* out_dM);
}

using coilgun::physics::mutual_inductance_coil;
using coilgun::physics::mutual_inductance_gradient_coil;
using coilgun::physics::mutual_inductance_coil_pair_kernel;
using coilgun::physics::upload_gl_nodes;

TEST_CASE("GPU coil pair M matches CPU") {
    upload_gl_nodes(9);

    double rai=0.01, rae=0.03, la=0.05;
    int na = 150;
    double rbi=0.005, rbe=0.01, lb=0.016;
    int nb = 1;

    std::vector<double> separations = {0.0, 0.03, 0.06, 0.10};
    for (double sep : separations) {
        double cpu_M = mutual_inductance_coil(rai, rae, la, na,
            rbi, rbe, lb, nb, sep, 9, false);
        double cpu_dM = mutual_inductance_gradient_coil(rai, rae, la, na,
            rbi, rbe, lb, nb, sep, 9, false);

        double *d_M, *d_dM;
        cudaMalloc(&d_M, sizeof(double));
        cudaMalloc(&d_dM, sizeof(double));
        mutual_inductance_coil_pair_kernel<<<1, 512>>>(
            rai, rae, la, na, rbi, rbe, lb, nb, sep, 9, d_M, d_dM);
        cudaDeviceSynchronize();

        double gpu_M, gpu_dM;
        cudaMemcpy(&gpu_M, d_M, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&gpu_dM, d_dM, sizeof(double), cudaMemcpyDeviceToHost);
        cudaFree(d_M); cudaFree(d_dM);

        CHECK(gpu_M  == doctest::Approx(cpu_M).epsilon(5e-7));
        CHECK(gpu_dM == doctest::Approx(cpu_dM).epsilon(5e-7));
    }
}

TEST_CASE("GPU coil pair dM antisymmetry") {
    upload_gl_nodes(9);

    double rai=0.01, rae=0.03, la=0.05;
    double rbi=0.005, rbe=0.01, lb=0.016;
    double sep = 0.05;

    double *d_M1, *d_dM1, *d_M2, *d_dM2;
    cudaMalloc(&d_M1, sizeof(double)); cudaMalloc(&d_dM1, sizeof(double));
    cudaMalloc(&d_M2, sizeof(double)); cudaMalloc(&d_dM2, sizeof(double));

    mutual_inductance_coil_pair_kernel<<<1, 512>>>(
        rai, rae, la, 150, rbi, rbe, lb, 1, sep, 9, d_M1, d_dM1);
    mutual_inductance_coil_pair_kernel<<<1, 512>>>(
        rai, rae, la, 150, rbi, rbe, lb, 1, -sep, 9, d_M2, d_dM2);
    cudaDeviceSynchronize();

    double gpu_dM1, gpu_dM2;
    cudaMemcpy(&gpu_dM1, d_dM1, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&gpu_dM2, d_dM2, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_M1); cudaFree(d_dM1); cudaFree(d_M2); cudaFree(d_dM2);

    CHECK(gpu_dM1 == doctest::Approx(-gpu_dM2).epsilon(1e-12));
}
