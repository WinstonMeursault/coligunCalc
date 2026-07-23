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
#include <limits>
#include <stdexcept>

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
    , dr_(0.0)
    , dl_(0.0)
    , nc_fil_(0.0)
{
    if (!std::isfinite(ri_) || ri_ < 0.0)
        throw std::invalid_argument("Armature inner_radius must be finite and non-negative");
    if (!std::isfinite(re_) || re_ <= ri_)
        throw std::invalid_argument("Armature outer_radius must exceed inner_radius");
    if (!std::isfinite(l_) || l_ <= 0.0)
        throw std::invalid_argument("Armature length must be positive");
    if (!std::isfinite(rho_) || rho_ <= 0.0)
        throw std::invalid_argument("Armature resistivity must be positive");
    if (!std::isfinite(material_density) || material_density <= 0.0)
        throw std::invalid_argument("Armature material_density must be positive");
    if (!std::isfinite(v_)) throw std::invalid_argument("Armature velocity must be finite");
    if (!std::isfinite(ma_) || ma_ <= 0.0)
        throw std::invalid_argument("Armature mass must be positive");
    if (!std::isfinite(x_)) throw std::invalid_argument("Armature position must be finite");
    if (m_ <= 0) throw std::invalid_argument("Armature m_axial must be positive");
    if (n_ <= 0) throw std::invalid_argument("Armature n_radial must be positive");
    if (static_cast<std::size_t>(m_) > std::numeric_limits<std::size_t>::max() /
                                      static_cast<std::size_t>(n_))
        throw std::invalid_argument("Armature filament count overflows");

    dr_ = (re_ - ri_) / static_cast<double>(n_);
    dl_ = l_ / static_cast<double>(m_);
    nc_fil_ = 1.0 / (dr_ * dl_);
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
        if (!std::isfinite(R_layer[j - 1]) || R_layer[j - 1] <= 0.0)
            throw std::invalid_argument("Armature derived filament resistance must be positive");
        if (!std::isfinite(L_layer[j - 1]) || L_layer[j - 1] <= 0.0)
            throw std::invalid_argument("Armature derived filament inductance must be positive");
        if (!std::isfinite(mass_layer[j - 1]) || mass_layer[j - 1] <= 0.0)
            throw std::invalid_argument("Armature derived filament mass must be positive");
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

void Armature::update_position(double dx) {
    if (!std::isfinite(dx)) throw std::invalid_argument("Armature displacement must be finite");
    x_ += dx;
}
void Armature::set_velocity(double v) {
    if (!std::isfinite(v)) throw std::invalid_argument("Armature velocity must be finite");
    v_ = v;
}

} // namespace coilgun::components
