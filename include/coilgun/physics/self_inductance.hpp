/**
 * @file self_inductance.hpp
 * @brief Self-inductance of hollow cylindrical coils (exact and table-based).
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::physics {

/**
 * @brief Self-inductance of a hollow cylindrical coil (auto-select method).
 * @param inner_radius Coil inner radius, m.
 * @param outer_radius Coil outer radius, m.
 * @param length Coil axial length, m.
 * @param turns_density Turns per unit cross-sectional area, turns/m^2
 *        (= N_turns / ((outer - inner) * length)).
 * @param force_exact If true, always use the reference-grade Bessel/Struve
 *        integration (bypasses the T(q,p) lookup table).
 * @return Self-inductance, H.
 *
 * When @p force_exact is false (default), uses the fast T(q,p) lookup table
 * when (q,p) fall within the table bounds (q in [0.05, 4], p in [1.05, 4]),
 * with bilinear interpolation yielding ~5e-7 relative error in T.
 * Automatically falls back to the exact Bessel/Struve integration for
 * out-of-range geometries.  When @p force_exact is true, always uses the
 * exact integration, corresponding to NumericalModel Sec. 4.1 Eqs. (4.1)-(4.5).
 */
double self_inductance(double inner_radius, double outer_radius,
                       double length, double turns_density,
                       bool force_exact = false);

/**
 * @brief Self-inductance of a hollow cylindrical coil (always reference-grade).
 * @param inner_radius Coil inner radius, m.
 * @param outer_radius Coil outer radius, m.
 * @param length Coil axial length, m.
 * @param turns_density Turns per unit cross-sectional area, turns/m^2
 *        (= N_turns / ((outer - inner) * length)).
 * @return Self-inductance, H.
 *
 * Uses the Bessel/Struve kernel with adaptive numerical integration.
 * Always exact — use when you need reference-grade accuracy
 * regardless of performance. Corresponds to NumericalModel Sec. 4.1
 * Eqs. (4.1)-(4.5).
 */
double self_inductance_exact(double inner_radius, double outer_radius,
                             double length, double turns_density);

/**
 * @brief Inductance shape factor T(q,p) — reference computation.
 * @param q Ratio length / inner_radius (dimensionless).
 * @param p Ratio outer_radius / inner_radius (dimensionless).
 * @return Dimensionless shape factor T(q,p).
 *
 * Uses the full Bessel/Struve kernel with composite Gauss-Legendre
 * integration.  This is the exact reference used to populate lookup
 * tables and should NOT be called in hot simulation loops.
 */
double inductance_shape_factor_reference(double q, double p);

} // namespace coilgun::physics
