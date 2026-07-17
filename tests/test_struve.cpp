#include <doctest/doctest.h>
#include "coilgun/physics/struve.hpp"

#include <algorithm>
#include <cmath>

using coilgun::physics::struve_h0;
using coilgun::physics::struve_h1;

// Reference values from scipy.special.struve.
// Tolerances: 1e-10 where power series is reliable, 5e-8 where
// asymptotic (Gauss-Laguerre + Y Bessel) must take over.

TEST_CASE("struve_h0 — known values") {
    CHECK(struve_h0(0.0) == doctest::Approx(0.0).epsilon(1e-15));

    // x <= 10 — power series always reliable
    CHECK(struve_h0(0.1) == doctest::Approx(6.359126999493e-02).epsilon(1e-10));
    CHECK(struve_h0(0.5) == doctest::Approx(3.095559145838e-01).epsilon(1e-10));
    CHECK(struve_h0(1.0) == doctest::Approx(5.686566270483e-01).epsilon(1e-10));
    CHECK(struve_h0(2.0) == doctest::Approx(7.908588495081e-01).epsilon(1e-10));
    CHECK(struve_h0(5.0) == doctest::Approx(-1.852168157767e-01).epsilon(1e-10));
    CHECK(struve_h0(10.0)== doctest::Approx( 1.187436836875e-01).epsilon(1e-10));

    // x >= 15 — asymptotic selected after PS crosscheck fails
    CHECK(struve_h0(15.0) == doctest::Approx( 2.477238309812e-01).epsilon(5e-8));
    CHECK(struve_h0(20.0) == doctest::Approx( 9.439369808132e-02).epsilon(5e-8));
    CHECK(struve_h0(25.0) == doctest::Approx(-1.018248201600e-01).epsilon(5e-8));
    CHECK(struve_h0(30.0) == doctest::Approx(-9.609842155416e-02).epsilon(5e-8));
    CHECK(struve_h0(50.0) == doctest::Approx(-8.533767482612e-02).epsilon(5e-8));
}

TEST_CASE("struve_h1 — known values") {
    CHECK(struve_h1(0.0) == doctest::Approx(0.0).epsilon(1e-15));

    CHECK(struve_h1(0.1) == doctest::Approx(2.120651601426e-03).epsilon(1e-10));
    CHECK(struve_h1(0.5) == doctest::Approx(5.217374424234e-02).epsilon(1e-10));
    CHECK(struve_h1(1.0) == doctest::Approx(1.984573362019e-01).epsilon(1e-10));
    CHECK(struve_h1(2.0) == doctest::Approx(6.467637282836e-01).epsilon(1e-10));
    CHECK(struve_h1(5.0) == doctest::Approx(8.078119457941e-01).epsilon(1e-10));
    CHECK(struve_h1(10.0)== doctest::Approx(8.918324920945e-01).epsilon(1e-10));

    CHECK(struve_h1(15.0) == doctest::Approx(6.604872985120e-01).epsilon(5e-8));
    CHECK(struve_h1(20.0) == doctest::Approx(4.726881842910e-01).epsilon(5e-8));
    CHECK(struve_h1(25.0) == doctest::Approx(5.388036213269e-01).epsilon(5e-8));
    CHECK(struve_h1(30.0) == doctest::Approx(7.217503783470e-01).epsilon(5e-8));
    CHECK(struve_h1(50.0) == doctest::Approx(5.800784479454e-01).epsilon(5e-8));
}

TEST_CASE("struve — negative x symmetry") {
    // H0 is odd:  H0(-x) = -H0(x)  (DLMF 11.2.1)
    CHECK(struve_h0(-5.0) == doctest::Approx(-struve_h0(5.0)).epsilon(1e-10));
    CHECK(struve_h0(-10.0) == doctest::Approx(-struve_h0(10.0)).epsilon(1e-8));
    // H1 is even: H1(-x) = H1(x)
    CHECK(struve_h1(-5.0) == doctest::Approx(struve_h1(5.0)).epsilon(1e-10));
    CHECK(struve_h1(-10.0) == doctest::Approx(struve_h1(10.0)).epsilon(1e-8));
}
