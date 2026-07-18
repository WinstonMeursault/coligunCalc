/**
 * @file test_gpu_filament.cpp
 * @brief Unit tests: GPU vs CPU filament-level mutual inductance.
 * @author Winston Meursault
 */

#include <doctest/doctest.h>
#include "coilgun/physics/mutual_inductance.hpp"
#include "coilgun/physics/mutual_inductance.cuh"
#include <cmath>
#include <cuda_runtime.h>
#include <vector>

using coilgun::physics::mutual_inductance_filament;
using coilgun::physics::mutual_inductance_gradient_filament;
using coilgun::physics::mutual_inductance_filament_device;
using coilgun::physics::mutual_inductance_gradient_filament_device;

__global__ void test_filament_m_kernel(double a, double b, double s, double* out) {
    *out = mutual_inductance_filament_device(a, b, s);
}

__global__ void test_filament_dm_kernel(double a, double b, double s, double* out) {
    *out = mutual_inductance_gradient_filament_device(a, b, s);
}

static double gpu_filament_m(double a, double b, double s) {
    double *d_out;
    cudaMalloc(&d_out, sizeof(double));
    test_filament_m_kernel<<<1, 1>>>(a, b, s, d_out);
    cudaDeviceSynchronize();
    double r;
    cudaMemcpy(&r, d_out, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_out);
    return r;
}

static double gpu_filament_dm(double a, double b, double s) {
    double *d_out;
    cudaMalloc(&d_out, sizeof(double));
    test_filament_dm_kernel<<<1, 1>>>(a, b, s, d_out);
    cudaDeviceSynchronize();
    double r;
    cudaMemcpy(&r, d_out, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_out);
    return r;
}

TEST_CASE("GPU filament M matches CPU — reference values") {
    struct Case { double a, b, s, expected; };
    std::vector<Case> cases = {
        {0.02, 0.03, 0.005, 2.959493865841560e-08},
        {0.0135, 0.0135, 0.001, 4.557715786266084e-08},
        {0.01, 0.01, 0.001, 3.002876303701474e-08},
        {0.01, 0.01, 0.01, 4.940784630798257e-09},
        {0.02, 0.01, 0.01, 6.987324633639440e-09},
    };
    for (auto& c : cases) {
        double cpu = mutual_inductance_filament(c.a, c.b, c.s);
        double gpu = gpu_filament_m(c.a, c.b, c.s);
        CHECK(cpu == doctest::Approx(c.expected).epsilon(1e-10));
        CHECK(gpu == doctest::Approx(cpu).epsilon(1e-13));
    }
}

TEST_CASE("GPU filament dM matches CPU — reference values") {
    struct Case { double a, b, s, expected; };
    std::vector<Case> cases = {
        {0.01, 0.01, 0.001, -1.239937005403076e-05},
        {0.02, 0.01, 0.01, -5.079612386972773e-07},
    };
    for (auto& c : cases) {
        double cpu = mutual_inductance_gradient_filament(c.a, c.b, c.s);
        double gpu = gpu_filament_dm(c.a, c.b, c.s);
        CHECK(cpu == doctest::Approx(c.expected).epsilon(1e-10));
        CHECK(gpu == doctest::Approx(cpu).epsilon(1e-13));
    }
}

TEST_CASE("GPU filament dM — antisymmetry") {
    double cpu = mutual_inductance_gradient_filament(0.01, 0.01, 0.001);
    double gpu = gpu_filament_dm(0.01, 0.01, -0.001);
    CHECK(gpu == doctest::Approx(-cpu).epsilon(1e-13));
}
