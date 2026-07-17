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

namespace coilgun::components {

DrivingCoil::DrivingCoil(double inner_radius, double outer_radius, double length,
                         int turns, double resistivity, double wire_area,
                         double fill_factor, double position,
                         bool force_exact_self_inductance)
    : ri_(inner_radius)
    , re_(outer_radius)
    , l_(length)
    , n_(turns)
    , rho_(resistivity)
    , wire_area_(wire_area)
    , k_fill_(fill_factor)
    , x_(position >= 0.0 ? position : length / 2.0)
    , nc_(turns / ((outer_radius - inner_radius) * length))
    , R_(rho_ * k_fill_ * M_PI * (re_ * re_ - ri_ * ri_) * l_ / (wire_area_ * wire_area_))
    , L_(physics::self_inductance(ri_, re_, l_, nc_, force_exact_self_inductance))
{ }

double DrivingCoil::inner_radius() const { return ri_; }
double DrivingCoil::outer_radius() const { return re_; }
double DrivingCoil::length() const { return l_; }
double DrivingCoil::mean_radius() const { return (ri_ + re_) / 2.0; }
int    DrivingCoil::turns() const { return n_; }

double DrivingCoil::turns_density() const { return nc_; }
double DrivingCoil::resistance() const { return R_; }
double DrivingCoil::self_inductance() const { return L_; }

double DrivingCoil::position() const { return x_; }
void   DrivingCoil::set_position(double x) { x_ = x; }

} // namespace coilgun::components
