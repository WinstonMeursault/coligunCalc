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
using coilgun::physics::mutual_inductance_filament_pair_device;
using coilgun::physics::mutual_inductance_filament_pair_f32;

__global__ void test_filament_m_kernel(double a, double b, double s, double* out) {
    *out = mutual_inductance_filament_device(a, b, s);
}

__global__ void test_filament_dm_kernel(double a, double b, double s, double* out) {
    *out = mutual_inductance_gradient_filament_device(a, b, s);
}

__global__ void test_filament_pair_kernel(double a, double b, double s, double* out) {
    const auto pair = mutual_inductance_filament_pair_device(a, b, s);
    out[0] = pair.mutual;
    out[1] = pair.gradient;
}

__global__ void test_filament_pair_f32_kernel(float a, float b, float s, float* out) {
    const auto pair = mutual_inductance_filament_pair_f32(a, b, s);
    out[0] = pair.mutual;
    out[1] = pair.gradient;
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

static std::vector<double> gpu_filament_pair(double a, double b, double s) {
    double* d_out;
    cudaMalloc(&d_out, 2 * sizeof(double));
    test_filament_pair_kernel<<<1, 1>>>(a, b, s, d_out);
    cudaDeviceSynchronize();
    std::vector<double> result(2);
    cudaMemcpy(result.data(), d_out, 2 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_out);
    return result;
}

static std::vector<float> gpu_filament_pair_f32(float a, float b, float s) {
    float* d_out;
    cudaMalloc(&d_out, 2 * sizeof(float));
    test_filament_pair_f32_kernel<<<1, 1>>>(a, b, s, d_out);
    cudaDeviceSynchronize();
    std::vector<float> result(2);
    cudaMemcpy(result.data(), d_out, 2 * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(d_out);
    return result;
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

TEST_CASE("GPU combined filament path returns M and signed dM") {
    for (double separation : {-0.01, 0.0, 0.01}) {
        const auto gpu = gpu_filament_pair(0.02, 0.01, separation);
        const double cpu_m = mutual_inductance_filament(0.02, 0.01, separation);
        const double cpu_dm = mutual_inductance_gradient_filament(0.02, 0.01, separation);
        CHECK(gpu[0] == doctest::Approx(cpu_m).epsilon(1e-13));
        CHECK(gpu[1] == doctest::Approx(cpu_dm).epsilon(1e-13));
    }
}

TEST_CASE("GPU combined FP32 filament path remains within aggressive tolerance") {
    for (float separation : {-0.01f, 0.0f, 0.01f}) {
        const auto gpu = gpu_filament_pair_f32(0.02f, 0.01f, separation);
        const double cpu_m = mutual_inductance_filament(0.02, 0.01, separation);
        const double cpu_dm = mutual_inductance_gradient_filament(0.02, 0.01, separation);
        CHECK(static_cast<double>(gpu[0]) == doctest::Approx(cpu_m).epsilon(1e-2));
        CHECK(static_cast<double>(gpu[1]) == doctest::Approx(cpu_dm).epsilon(1e-2));
    }
}
