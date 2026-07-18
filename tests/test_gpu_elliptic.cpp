/**
 * @file test_gpu_elliptic.cpp
 * @brief Unit tests: GPU vs CPU elliptic integrals.
 * @author Winston Meursault
 */

#include <doctest/doctest.h>
#include "coilgun/physics/elliptic.hpp"
#include "coilgun/physics/elliptic.cuh"
#include <algorithm>
#include <cmath>
#include <cuda_runtime.h>
#include <vector>

namespace coilgun::physics {
    extern __global__ void elliptic_kernel(double m, double* out_k, double* out_e);
    extern __global__ void elliptic_modulus_kernel(double a, double b, double sep, double* out);
}

using coilgun::physics::elliptic_k;
using coilgun::physics::elliptic_e;
using coilgun::physics::elliptic_modulus;
using coilgun::physics::elliptic_kernel;
using coilgun::physics::elliptic_modulus_kernel;

static double run_elliptic_k_gpu(double m) {
    double *d_k, *d_e;
    cudaMalloc(&d_k, sizeof(double));
    cudaMalloc(&d_e, sizeof(double));
    elliptic_kernel<<<1, 1>>>(m, d_k, d_e);
    cudaDeviceSynchronize();
    double result;
    cudaMemcpy(&result, d_k, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_k);
    cudaFree(d_e);
    return result;
}

static double run_elliptic_e_gpu(double m) {
    double *d_k, *d_e;
    cudaMalloc(&d_k, sizeof(double));
    cudaMalloc(&d_e, sizeof(double));
    elliptic_kernel<<<1, 1>>>(m, d_k, d_e);
    cudaDeviceSynchronize();
    double result;
    cudaMemcpy(&result, d_e, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_k);
    cudaFree(d_e);
    return result;
}

static double run_modulus_gpu(double a, double b, double s) {
    double *d_out;
    cudaMalloc(&d_out, sizeof(double));
    elliptic_modulus_kernel<<<1, 1>>>(a, b, s, d_out);
    cudaDeviceSynchronize();
    double result;
    cudaMemcpy(&result, d_out, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_out);
    return result;
}

TEST_CASE("GPU elliptic_k matches CPU") {
    std::vector<double> test_vals = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99};
    for (double m : test_vals) {
        double cpu = elliptic_k(m);
        double gpu = run_elliptic_k_gpu(m);
        CHECK(doctest::Approx(cpu).epsilon(1e-14) == gpu);
    }
}

TEST_CASE("GPU elliptic_e matches CPU") {
    std::vector<double> test_vals = {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99};
    for (double m : test_vals) {
        double cpu = elliptic_e(m);
        double gpu = run_elliptic_e_gpu(m);
        CHECK(doctest::Approx(cpu).epsilon(1e-14) == gpu);
    }
}

TEST_CASE("GPU elliptic_modulus matches CPU") {
    struct TestCase { double a, b, s; };
    std::vector<TestCase> cases = {
        {0.01, 0.01, 0.001}, {0.02, 0.03, 0.01}, {0.01, 0.05, 0.02}};
    for (auto tc : cases) {
        double cpu = elliptic_modulus(tc.a, tc.b, tc.s);
        double gpu = run_modulus_gpu(tc.a, tc.b, tc.s);
        CHECK(doctest::Approx(cpu).epsilon(1e-14) == gpu);
    }
}

TEST_CASE("GPU elliptic edge cases") {
    CHECK(run_elliptic_k_gpu(0.0) == doctest::Approx(elliptic_k(0.0)).epsilon(1e-14));
    CHECK(run_elliptic_e_gpu(0.0) == doctest::Approx(elliptic_e(0.0)).epsilon(1e-14));
}
