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
 * @return Nodes and weights.
 */
QuadratureNodes gauss_legendre(int n);

/**
 * @brief Gauss-Laguerre quadrature nodes and weights on [0, +inf).
 * @param n Number of quadrature points.
 * @return Nodes and weights.
 */
QuadratureNodes gauss_laguerre(int n);

} // namespace coilgun::physics
