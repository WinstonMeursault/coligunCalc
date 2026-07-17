#include <cmath>
#include <stdexcept>

#include "coilgun/physics/quadrature.hpp"
#include "doctest/doctest.h"

using coilgun::physics::gauss_legendre;
using coilgun::physics::gauss_laguerre;

TEST_CASE("gauss_legendre(4) integrates x^6 exactly") {
    auto q = gauss_legendre(4);

    CHECK(q.nodes.size() == 4);
    CHECK(q.weights.size() == 4);

    double integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        double x = q.nodes[i];
        integral += q.weights[i] * std::pow(x, 6);
    }

    // ∫_{-1}^{1} x^6 dx = 2/7 ≈ 0.2857142857142857
    CHECK(integral == doctest::Approx(2.0 / 7.0));
}

TEST_CASE("gauss_legendre(9) integrates x^10 and x^17 exactly") {
    auto q = gauss_legendre(9);
    double integral;

    // ∫_{-1}^{1} x^{10} dx = 2/11 ≈ 0.181818...
    integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        integral += q.weights[i] * std::pow(q.nodes[i], 10);
    }
    CHECK(integral == doctest::Approx(2.0 / 11.0));

    // ∫_{-1}^{1} x^{16} dx = 2/17 ≈ 0.117647... — 9-point exact up to degree 2n-1=17
    integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        integral += q.weights[i] * std::pow(q.nodes[i], 16);
    }
    CHECK(integral == doctest::Approx(2.0 / 17.0));
}

TEST_CASE("gauss_laguerre(15) integrates x^2 e^{-x} → Γ(3)=2") {
    auto q = gauss_laguerre(15);

    CHECK(q.nodes.size() == 15);
    CHECK(q.weights.size() == 15);

    double integral = 0.0;
    for (std::size_t i = 0; i < q.nodes.size(); ++i) {
        double x = q.nodes[i];
        integral += q.weights[i] * x * x; // w_i already absorbs e^{-x}
    }

    CHECK(integral == doctest::Approx(2.0));
}

TEST_CASE("gauss_legendre weight sum equals interval length") {
    for (int n : {4, 9, 16, 32}) {
        auto q = gauss_legendre(n);
        double sum = 0.0;
        for (auto w : q.weights) {
            sum += w;
        }
        CHECK(sum == doctest::Approx(2.0)); // interval [-1, 1] has length 2
    }
}

TEST_CASE("gauss_laguerre weight sum equals 1") {
    for (int n : {15, 30}) {
        auto q = gauss_laguerre(n);
        double sum = 0.0;
        for (auto w : q.weights) {
            sum += w;
        }
        CHECK(sum == doctest::Approx(1.0));
    }
}

TEST_CASE("gauss_legendre rejects unsupported n") {
    CHECK_THROWS_AS(gauss_legendre(3), std::invalid_argument);
    CHECK_THROWS_AS(gauss_legendre(0), std::invalid_argument);
}

TEST_CASE("gauss_laguerre rejects unsupported n") {
    CHECK_THROWS_AS(gauss_laguerre(10), std::invalid_argument);
    CHECK_THROWS_AS(gauss_laguerre(0), std::invalid_argument);
}
