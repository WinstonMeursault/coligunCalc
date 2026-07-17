#include <doctest/doctest.h>
#include "coilgun/components/driving_coil.hpp"

#include <cmath>

using coilgun::components::DrivingCoil;

namespace {

constexpr double k_ri   = 0.015125;
constexpr double k_re   = 0.029125;
constexpr double k_l    = 0.063;
constexpr int    k_n    = 63;
constexpr double k_rho  = 1.75e-8;
constexpr double k_sw   = 2.8e-5;
constexpr double k_kf   = 0.7;

} // namespace

TEST_CASE("DrivingCoil — turns_density") {
    DrivingCoil dc(k_ri, k_re, k_l, k_n, k_rho, k_sw, k_kf);
    double expected = k_n / ((k_re - k_ri) * k_l);
    CHECK(dc.turns_density() == doctest::Approx(expected).epsilon(1e-10));
}

TEST_CASE("DrivingCoil — resistance") {
    DrivingCoil dc(k_ri, k_re, k_l, k_n, k_rho, k_sw, k_kf);
    double expected = k_rho * k_kf * M_PI * (k_re * k_re - k_ri * k_ri)
                      * k_l / (k_sw * k_sw);
    CHECK(dc.resistance() == doctest::Approx(expected).epsilon(1e-10));
}

TEST_CASE("DrivingCoil — self_inductance") {
    DrivingCoil dc(k_ri, k_re, k_l, k_n, k_rho, k_sw, k_kf);
    double expected = 7.191793975090e-05;
    CHECK(dc.self_inductance() == doctest::Approx(expected).epsilon(1e-10));
}

TEST_CASE("DrivingCoil — getters") {
    DrivingCoil dc(k_ri, k_re, k_l, k_n, k_rho, k_sw, k_kf);
    CHECK(dc.inner_radius() == k_ri);
    CHECK(dc.outer_radius() == k_re);
    CHECK(dc.length() == k_l);
    CHECK(dc.turns() == k_n);
    CHECK(dc.mean_radius() == doctest::Approx((k_ri + k_re) / 2.0));
}

TEST_CASE("DrivingCoil — position") {
    DrivingCoil dc(k_ri, k_re, k_l, k_n, k_rho, k_sw, k_kf);
    CHECK(dc.position() == doctest::Approx(k_l / 2.0));
    dc.set_position(0.1);
    CHECK(dc.position() == doctest::Approx(0.1));
}
