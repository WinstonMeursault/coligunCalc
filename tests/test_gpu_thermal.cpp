#include <doctest/doctest.h>

#include "coilgun/physics/constants.hpp"
#include "coilgun/simulation/cuda/gpu_thermal.hpp"

#include <cuda_runtime.h>
#include <vector>

using namespace coilgun::simulation::cuda;

TEST_CASE("GPU thermal tables interpolate CPU material properties") {
    const auto tables = generate_material_tables(257, 293.0, 1293.0);
    CHECK(interpolate_material_cp(tables, ThermalMaterial::Copper, 673.25, ThermalPrecision::Full) ==
          doctest::Approx(coilgun::physics::specific_heat_capacity_copper(673.25)).epsilon(2e-6));
    CHECK(interpolate_material_resistivity(tables, ThermalMaterial::Aluminum, 673.25,
                                           ThermalPrecision::Full) ==
          doctest::Approx(coilgun::physics::resistivity_aluminum(673.25)).epsilon(1e-10));
}

TEST_CASE("GPU thermal batch returns temperature resistance rho and Joule energy") {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        MESSAGE("CUDA device unavailable; skipping thermal-kernel test");
        return;
    }
    const auto tables = generate_material_tables(1025, 293.0, 2000.0);
    constexpr std::size_t B = 2;
    constexpr std::size_t F = 2;
    const double currents[B * F] = {10.0, -7.0, 3.0, 4.0};
    const double masses[B * F] = {0.01, 0.02, 0.03, 0.04};
    const double r0[B * F] = {1e-3, 2e-3, 3e-3, 4e-3};
    const int material[B * F] = {1, 0, 1, 0};
    std::vector<double> t = {293.0, 300.0, 500.0, 900.0};
    std::vector<double> rho(B * F), resistance(B * F), joule(B * F);

    update_thermal_batch(tables, ThermalPrecision::Standard, B, F, currents, masses, r0,
                         material, 1e-3, t.data(), rho.data(), resistance.data(), joule.data());

    CHECK(t[0] > 293.0);
    CHECK(resistance[0] > r0[0]);
    CHECK(rho[0] > coilgun::physics::resistivity_copper(293.0));
    CHECK(joule[0] == doctest::Approx(currents[0] * currents[0] * r0[0] * 1e-3).epsilon(1e-12));
}

TEST_CASE("GPU thermal precision modes stay within CPU reference thresholds") {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        MESSAGE("CUDA device unavailable; skipping thermal precision test");
        return;
    }
    const auto tables = generate_material_tables(1025, 293.0, 2000.0);
    constexpr std::size_t count = 4;
    const double currents[count] = {31.0, -17.0, 8.0, 2.5};
    const double masses[count] = {0.002, 0.004, 0.01, 0.02};
    const double r0[count] = {8e-4, 1.2e-3, 2e-3, 3e-3};
    const double initial_t[count] = {293.0, 420.0, 800.0, 1200.0};
    const int material[count] = {1, 0, 1, 0};
    const double dt = 2e-3;

    for (const auto mode : {ThermalPrecision::Standard, ThermalPrecision::Full,
                            ThermalPrecision::Aggressive}) {
        std::vector<double> t(initial_t, initial_t + count);
        std::vector<double> rho(count), resistance(count), joule(count);
        update_thermal_batch(tables, mode, 1, count, currents, masses, r0, material, dt,
                             t.data(), rho.data(), resistance.data(), joule.data());
        const double temperature_epsilon = mode == ThermalPrecision::Aggressive ? 2e-4 : 2e-8;
        const double property_epsilon = mode == ThermalPrecision::Aggressive ? 2e-4 : 2e-8;
        for (std::size_t k = 0; k < count; ++k) {
            const auto mat = material[k] ? coilgun::physics::ArmatureMaterial::Copper
                                         : coilgun::physics::ArmatureMaterial::Aluminum;
            const double cp = coilgun::physics::material_cp(mat, initial_t[k]);
            const double expected_t = initial_t[k] + currents[k] * currents[k] * r0[k] * dt /
                                      (masses[k] * cp);
            const double expected_rho = material[k]
                ? coilgun::physics::resistivity_copper(expected_t)
                : coilgun::physics::resistivity_aluminum(expected_t);
            const double beta = material[k] ? coilgun::physics::COPPER.temp_coefficient
                                            : coilgun::physics::ALUMINUM.temp_coefficient;
            const double expected_r = r0[k] * (1.0 + beta * (expected_t - 293.0));
            CHECK(t[k] == doctest::Approx(expected_t).epsilon(temperature_epsilon));
            CHECK(rho[k] == doctest::Approx(expected_rho).epsilon(property_epsilon));
            CHECK(resistance[k] == doctest::Approx(expected_r).epsilon(property_epsilon));
            CHECK(joule[k] == doctest::Approx(currents[k] * currents[k] * r0[k] * dt)
                              .epsilon(mode == ThermalPrecision::Aggressive ? 2e-4 : 2e-8));
        }
    }
}

TEST_CASE("CPU thermal batch uses temperature-dependent specific heat") {
    const auto tables = generate_material_tables(1025, 293.0, 2000.0);
    constexpr double current = 40.0;
    constexpr double mass = 0.01;
    constexpr double reference_resistance = 2e-3;
    constexpr double dt = 0.01;
    double temperature = 900.0;
    double resistivity = 0.0;
    double resistance = 0.0;
    double joule = 0.0;
    const int material = 0;

    update_thermal_batch_cpu(tables, ThermalPrecision::Standard, 1, 1,
                             &current, &mass, &reference_resistance, &material,
                             dt, &temperature, &resistivity, &resistance, &joule);

    const double expected = 900.0 + current * current * reference_resistance * dt /
                            (mass * coilgun::physics::material_cp(
                                coilgun::physics::ArmatureMaterial::Aluminum, 900.0));
    CHECK(temperature == doctest::Approx(expected).epsilon(1e-12));
}
