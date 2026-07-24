/**
 * @file quadrature.hpp
 * @brief Gaussian quadrature rule generators (Gauss-Legendre, Gauss-Laguerre).
 * @author Winston Meursault
 */

#pragma once

#include <vector>

namespace coilgun::physics {

/**
 * @brief Result of a Gaussian quadrature rule: nodes and weights.
 */
struct QuadratureNodes {
    std::vector<double> nodes;   ///< Quadrature node positions.
    std::vector<double> weights; ///< Quadrature weights.
};

/**
 * @brief Gauss-Legendre quadrature nodes and weights on [-1, 1].
 * @param n Number of quadrature points.
 * @return A value containing nodes and weights.
 */
QuadratureNodes gauss_legendre(int n);

/**
 * @brief Return the immutable cached Gauss-Legendre rule.
 * @param n Number of quadrature points.
 * @return A process-lifetime cached rule.
 */
const QuadratureNodes& gauss_legendre_cached(int n);

/**
 * @brief Gauss-Laguerre quadrature nodes and weights on [0, +inf).
 * @param n Number of quadrature points.
 * @return A value containing nodes and weights.
 */
QuadratureNodes gauss_laguerre(int n);

/**
 * @brief Return the immutable cached Gauss-Laguerre rule.
 * @param n Number of quadrature points.
 * @return A process-lifetime cached rule.
 */
const QuadratureNodes& gauss_laguerre_cached(int n);

} // namespace coilgun::physics
