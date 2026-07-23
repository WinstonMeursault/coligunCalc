/**
 * @file driving_coil.cpp
 * @brief Single multi-turn helical driving coil (stator) implementation.
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.8.2.
 */

#include "coilgun/components/driving_coil.hpp"
#include "coilgun/physics/self_inductance.hpp"
#include "coilgun/physics/constants.hpp"

#include <cmath>
#include <stdexcept>

namespace coilgun::components {

DrivingCoil::DrivingCoil(double inner_radius, double outer_radius, double length,
                         int turns, double resistivity, double wire_area,
                         double fill_factor, std::optional<double> position,
                         bool force_exact_self_inductance)
    : ri_(inner_radius)
    , re_(outer_radius)
    , l_(length)
    , n_(turns)
    , rho_(resistivity)
    , wire_area_(wire_area)
    , k_fill_(fill_factor)
    , x_(0.0)
    , nc_(0.0)
    , R_(0.0)
    , L_(0.0)
{
    if (!std::isfinite(ri_) || ri_ < 0.0)
        throw std::invalid_argument("DrivingCoil inner_radius must be finite and non-negative");
    if (!std::isfinite(re_) || re_ <= ri_)
        throw std::invalid_argument("DrivingCoil outer_radius must exceed inner_radius");
    if (!std::isfinite(l_) || l_ <= 0.0)
        throw std::invalid_argument("DrivingCoil length must be positive");
    if (n_ <= 0) throw std::invalid_argument("DrivingCoil turns must be positive");
    if (!std::isfinite(rho_) || rho_ <= 0.0)
        throw std::invalid_argument("DrivingCoil resistivity must be positive");
    if (!std::isfinite(wire_area_) || wire_area_ <= 0.0)
        throw std::invalid_argument("DrivingCoil wire_area must be positive");
    if (!std::isfinite(k_fill_) || k_fill_ <= 0.0 || k_fill_ > 1.0)
        throw std::invalid_argument("DrivingCoil fill_factor must be in (0, 1]");
    if (position && !std::isfinite(*position))
        throw std::invalid_argument("DrivingCoil position must be finite");

    x_ = position.value_or(l_ / 2.0);
    nc_ = static_cast<double>(n_) / ((re_ - ri_) * l_);
    R_ = rho_ * k_fill_ * M_PI * (re_ * re_ - ri_ * ri_) * l_ /
         (wire_area_ * wire_area_);
    L_ = physics::self_inductance(ri_, re_, l_, nc_, force_exact_self_inductance);
    if (!std::isfinite(R_) || R_ <= 0.0)
        throw std::invalid_argument("DrivingCoil derived resistance must be positive");
    if (!std::isfinite(L_) || L_ <= 0.0)
        throw std::invalid_argument("DrivingCoil derived inductance must be positive");
}

double DrivingCoil::inner_radius() const { return ri_; }
double DrivingCoil::outer_radius() const { return re_; }
double DrivingCoil::length() const { return l_; }
double DrivingCoil::mean_radius() const { return (ri_ + re_) / 2.0; }
int    DrivingCoil::turns() const { return n_; }

double DrivingCoil::turns_density() const { return nc_; }
double DrivingCoil::resistance() const { return R_; }
double DrivingCoil::self_inductance() const { return L_; }

double DrivingCoil::position() const { return x_; }
void DrivingCoil::set_position(double x) {
    if (!std::isfinite(x))
        throw std::invalid_argument("DrivingCoil position must be finite");
    x_ = x;
}

} // namespace coilgun::components
