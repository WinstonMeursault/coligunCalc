#include <cmath>
#include <doctest/doctest.h>

#include "coilgun/physics/elliptic.hpp"

using coilgun::physics::elliptic_e;
using coilgun::physics::elliptic_k;
using coilgun::physics::elliptic_modulus;

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double PI_OVER_2 = PI / 2.0;

} // namespace

TEST_CASE("elliptic_k") {
    CHECK(elliptic_k(0.0) == doctest::Approx(PI_OVER_2).epsilon(1e-12));
    CHECK(elliptic_k(0.5) == doctest::Approx(1.8540746773013719).epsilon(1e-12));
}

TEST_CASE("elliptic_e") {
    CHECK(elliptic_e(0.0) == doctest::Approx(PI_OVER_2).epsilon(1e-12));
    CHECK(elliptic_e(0.5) == doctest::Approx(1.3506438810476755).epsilon(1e-12));
}

TEST_CASE("elliptic_modulus") {
    // Zero separation: k = sqrt(4ab / (a+b)^2)
    // a=b => k = sqrt(4a^2 / 4a^2) = 1.0 (clamped to 1 - 1e-10)
    double k_eq = elliptic_modulus(1.0, 1.0, 0.0);
    CHECK(k_eq == doctest::Approx(1.0 - 1e-10).epsilon(1e-12));

    // a=1, b=1, h=1: k = sqrt(4 / (4+1)) = sqrt(0.8)
    double k1 = elliptic_modulus(1.0, 1.0, 1.0);
    CHECK(k1 == doctest::Approx(std::sqrt(0.8)).epsilon(1e-12));

    // a=2, b=1, h=0.5
    double k2 = elliptic_modulus(2.0, 1.0, 0.5);
    double expected_k2 = std::sqrt(4.0 * 2.0 * 1.0 / (9.0 + 0.25));
    CHECK(k2 == doctest::Approx(expected_k2).epsilon(1e-12));

    // Symmetry: elliptic_modulus(a,b,d) == elliptic_modulus(b,a,d)
    CHECK(elliptic_modulus(2.0, 3.0, 0.5)
          == doctest::Approx(elliptic_modulus(3.0, 2.0, 0.5)).epsilon(1e-12));

    // Clamping: very small radii produce near-zero k, clamped to 1e-12
    double k_tiny = elliptic_modulus(1e-20, 1e-20, 0.0);
    CHECK(k_tiny >= 1e-12);
}
