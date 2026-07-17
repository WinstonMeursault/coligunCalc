/**
 * @file armature.cpp
 * @brief Solid cylindrical armature discretised into current filaments.
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.2.3, Sec.8.2.
 */

#include "coilgun/components/armature.hpp"
#include "coilgun/physics/self_inductance.hpp"

#include <cmath>

namespace coilgun::components {

Armature::Armature(double inner_radius, double outer_radius, double length,
                   double resistivity, double material_density,
                   double velocity, double mass,
                   int m_axial, int n_radial, double position,
                   physics::ArmatureMaterial material,
                   bool force_exact_self_inductance)
    : ri_(inner_radius)
    , re_(outer_radius)
    , l_(length)
    , rho_(resistivity)
    , m_(m_axial)
    , n_(n_radial)
    , x_(position)
    , v_(velocity)
    , ma_(mass)
    , material_(material)
    , dr_((outer_radius - inner_radius) / n_radial)
    , dl_(length / m_axial)
    , nc_fil_(1.0 / (dr_ * dl_))
{
    int total = m_ * n_;
    R_.reserve(total);
    L_.reserve(total);
    mass_.reserve(total);

    // Precompute per-radial-layer values, then expand in row-major order
    // (i=1,j=1), (i=1,j=2), ..., (i=m,j=n) — matching Python's R * m pattern.
    std::vector<double> R_layer(n_);
    std::vector<double> L_layer(n_);
    std::vector<double> mass_layer(n_);

    for (int j = 1; j <= n_; ++j) {
        double r_inner = filament_inner_radius(j);
        double r_outer = filament_outer_radius(j);

        R_layer[j - 1] = 2.0 * M_PI * rho_
            * (m_ / (2.0 * l_) + static_cast<double>(m_ * n_) * ri_ / (l_ * (re_ - ri_))
               + (j - 1) * m_ / l_);

        L_layer[j - 1] = coilgun::physics::self_inductance(r_inner, r_outer, dl_, nc_fil_,
                                                              force_exact_self_inductance);

        mass_layer[j - 1] = material_density * M_PI * (r_outer * r_outer - r_inner * r_inner) * dl_;
    }

    // Expand to row-major flat array: repeat each entry m times
    for (int i = 0; i < m_; ++i) {
        for (int j = 0; j < n_; ++j) {
            R_.push_back(R_layer[j]);
            L_.push_back(L_layer[j]);
            mass_.push_back(mass_layer[j]);
        }
    }
}

double Armature::inner_radius() const { return ri_; }
double Armature::outer_radius() const { return re_; }
double Armature::length() const { return l_; }
int    Armature::axial_filaments() const { return m_; }
int    Armature::radial_filaments() const { return n_; }
int    Armature::total_filaments() const { return m_ * n_; }

double Armature::filament_inner_radius(int j) const {
    return ri_ + (re_ - ri_) * (j - 1) / n_;
}

double Armature::filament_outer_radius(int j) const {
    return ri_ + (re_ - ri_) * j / n_;
}

double Armature::filament_mean_radius(int j) const {
    return ri_ + dr_ * (j - 0.5);
}

double Armature::filament_axial_position(int i) const {
    return x_ - 0.5 * l_ + (i - 0.5) * dl_;
}

const std::vector<double>& Armature::resistances() const { return R_; }
const std::vector<double>& Armature::inductances() const { return L_; }
const std::vector<double>& Armature::masses() const { return mass_; }

double Armature::position() const { return x_; }
double Armature::velocity() const { return v_; }
double Armature::mass() const { return ma_; }
physics::ArmatureMaterial Armature::material() const { return material_; }

void Armature::update_position(double dx) { x_ += dx; }
void Armature::set_velocity(double v) { v_ = v; }

} // namespace coilgun::components
