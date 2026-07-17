/**
 * @file mutual_inductance.hpp
 * @brief Mutual inductance and spatial gradient at filament and coil level.
 * @author Winston Meursault
 */

#pragma once

namespace coilgun::physics {

/**
 * @brief Mutual inductance between two coaxial circular current filament loops.
 * @param radius_a Radius of the first loop, m.
 * @param radius_b Radius of the second loop, m.
 * @param separation Axial separation between loop planes, m.
 * @return Mutual inductance, H.
 *
 * Closed-form via complete elliptic integrals.
 * NumericalModel Sec. 4.2, Eq. (4.8).
 */
double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation);

double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation, bool use_cache);

/**
 * @brief Spatial gradient dM/dz of mutual inductance between two coaxial loops.
 * @param radius_a Radius of the first loop, m.
 * @param radius_b Radius of the second loop, m.
 * @param separation Axial separation between loop planes, m.
 * @return dM/dz gradient, H/m.
 *
 * Closed-form. NumericalModel Sec. 4.4, Eq. (4.15).
 */
double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation);

double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation, bool use_cache);

/**
 * @brief Mutual inductance between two thick cylindrical coils.
 * @param r_inner_a Inner radius of coil A, m.
 * @param r_outer_a Outer radius of coil A, m.
 * @param length_a Axial length of coil A, m.
 * @param turns_a Number of turns in coil A.
 * @param r_inner_b Inner radius of coil B, m.
 * @param r_outer_b Outer radius of coil B, m.
 * @param length_b Axial length of coil B, m.
 * @param turns_b Number of turns in coil B.
 * @param separation Centre-to-centre axial distance, m.
 * @param n_nodes Number of Gauss-Legendre nodes per dimension (default 9,
 *        4 for fast approximation). Total kernel evaluations = n_nodes^4.
 * @return Mutual inductance, H.
 *
 * 4D Gauss-Legendre tensor-product quadrature.
 * NumericalModel Sec. 4.3, Eq. (4.13).
 */
double mutual_inductance_coil(double r_inner_a, double r_outer_a,
                              double length_a, int turns_a,
                              double r_inner_b, double r_outer_b,
                              double length_b, int turns_b,
                              double separation, int n_nodes = 9);

double mutual_inductance_coil(double r_inner_a, double r_outer_a,
                              double length_a, int turns_a,
                              double r_inner_b, double r_outer_b,
                              double length_b, int turns_b,
                              double separation, int n_nodes, bool use_cache);

/**
 * @brief Spatial gradient dM/dz of mutual inductance at the coil level.
 * @param r_inner_a Inner radius of coil A, m.
 * @param r_outer_a Outer radius of coil A, m.
 * @param length_a Axial length of coil A, m.
 * @param turns_a Number of turns in coil A.
 * @param r_inner_b Inner radius of coil B, m.
 * @param r_outer_b Outer radius of coil B, m.
 * @param length_b Axial length of coil B, m.
 * @param turns_b Number of turns in coil B.
 * @param separation Centre-to-centre axial distance, m.
 * @param n_nodes Number of Gauss-Legendre nodes per dimension (default 9,
 *        4 for fast approximation).
 * @return dM/dz gradient, H/m.
 *
 * NumericalModel Sec. 4.4, Eq. (4.16).
 */
double mutual_inductance_gradient_coil(double r_inner_a, double r_outer_a,
                                       double length_a, int turns_a,
                                       double r_inner_b, double r_outer_b,
                                       double length_b, int turns_b,
                                       double separation, int n_nodes = 9);

double mutual_inductance_gradient_coil(double r_inner_a, double r_outer_a,
                                       double length_a, int turns_a,
                                       double r_inner_b, double r_outer_b,
                                       double length_b, int turns_b,
                                       double separation, int n_nodes,
                                       bool use_cache);

} // namespace coilgun::physics
