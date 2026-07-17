/**
 * @file lookup_tables.cpp
 * @brief T(q,p) inductance shape-factor lookup with bilinear interpolation.
 * @author Winston Meursault
 */

#include "coilgun/physics/lookup_tables.hpp"
#include "coilgun/physics/lookup_table_data.hpp"

#include <algorithm>
#include <cmath>

namespace coilgun::physics {

double inductance_shape_factor(double q, double p) {
    using detail::k_table_nq;
    using detail::k_table_np;
    using detail::k_table_dq;
    using detail::k_table_dp;
    using detail::k_table_q0;
    using detail::k_table_p0;
    using detail::k_t_table_data;

    // Clamp to table bounds
    q = std::clamp(q, k_table_q0, detail::k_table_q_max);
    p = std::clamp(p, k_table_p0, detail::k_table_p_max);

    // Fractional indices
    const double qi = (q - k_table_q0) / k_table_dq;
    const double pi = (p - k_table_p0) / k_table_dp;

    int q0 = static_cast<int>(qi);
    int q1 = q0 + 1;
    int p0 = static_cast<int>(pi);
    int p1 = p0 + 1;

    // Boundary clamp
    if (q0 >= k_table_nq - 1) { q0 = k_table_nq - 1; q1 = k_table_nq - 1; }
    if (q1 >= k_table_nq)     { q1 = k_table_nq - 1; }
    if (p0 >= k_table_np - 1) { p0 = k_table_np - 1; p1 = k_table_np - 1; }
    if (p1 >= k_table_np)     { p1 = k_table_np - 1; }

    // Interpolation weights
    const double t = qi - static_cast<double>(q0);
    const double s = pi - static_cast<double>(p0);

    // Row-major index: idx = row * n_cols + col
    const auto at = [](int r, int c) -> double {
        return k_t_table_data[r * k_table_np + c];
    };

    const double v00 = at(q0, p0);
    const double v10 = at(q1, p0);
    const double v01 = at(q0, p1);
    const double v11 = at(q1, p1);

    // Bilinear interpolation
    return v00 * (1.0 - t) * (1.0 - s)
         + v10 * t * (1.0 - s)
         + v01 * (1.0 - t) * s
         + v11 * t * s;
}

} // namespace coilgun::physics
