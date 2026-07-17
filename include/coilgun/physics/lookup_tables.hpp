/**
 * @file lookup_tables.hpp
 * @brief Precomputed shape-factor lookup for fast self-inductance approximation.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::physics {

/**
 * @brief Shape factor T(q, p) for calculating coil self-inductance.
 * @param q Ratio length / inner_radius (dimensionless), q in [0.05, 4.0].
 * @param p Ratio outer_radius / inner_radius (dimensionless), p in [1.05, 4.0].
 * @return Dimensionless shape factor T(q,p).
 *
 * Precomputed dense 1976×1476 lookup table (dq=dp=0.002) with
 * bilinear interpolation.  Bilinear interpolation error ~ 5×10^{-7},
 * negligible for engineering coilgun simulation.
 *
 * @see inductance_shape_factor_reference() for the reference integration.
 */
double inductance_shape_factor(double q, double p);

} // namespace coilgun::physics
