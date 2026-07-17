/**
 * @file constants.cpp
 * @brief Physical constants and material property tables.
 * @author Winston Meursault
 *
 * Implements temperature-dependent resistivity and specific heat capacity
 * for copper and aluminum.
 *
 * @see NumericalModel Sec.6.2, Sec.6.3.
 */

#include "coilgun/physics/constants.hpp"

namespace coilgun::physics {

// ---- material property tables ----

const MaterialProperties COPPER = {
    .resistivity_ref  = 1.75e-8,   // ohm·m  at 293 K
    .temp_coefficient = 4.1e-3,    // K^{-1}
    .density           = 8960.0,   // kg/m^3
};

const MaterialProperties ALUMINUM = {
    .resistivity_ref  = 2.82e-8,   // ohm·m  at 293 K
    .temp_coefficient = 4.2e-3,    // K^{-1}
    .density           = 2700.0,   // kg/m^3
};

// ---- material dispatch ----

double material_cp(ArmatureMaterial material, double T) {
    switch (material) {
    case ArmatureMaterial::Aluminum:
        return specific_heat_capacity_aluminum(T);
    case ArmatureMaterial::Copper:
        return specific_heat_capacity_copper(T);
    }
    throw std::invalid_argument("material_cp: unknown ArmatureMaterial");
}

double material_beta(ArmatureMaterial material) {
    switch (material) {
    case ArmatureMaterial::Aluminum:
        return ALUMINUM.temp_coefficient;
    case ArmatureMaterial::Copper:
        return COPPER.temp_coefficient;
    }
    throw std::invalid_argument("material_beta: unknown ArmatureMaterial");
}

// ---- temperature-dependent specific heat capacity ----

// Eq. (6.3): c_Cu(T) = 0.333 * exp(3.917e-4 * T)   kJ/(kg·K)
// Multiply by 1000 to obtain J/(kg·K).
double specific_heat_capacity_copper(double T) {
    return 333.0 * std::exp(3.917e-4 * T);  // J/(kg·K)
}

// Eq. (6.2): c_Al(T) = 0.819 * exp(3.07e-4 * T)   kJ/(kg·K)
double specific_heat_capacity_aluminum(double T) {
    return 819.0 * std::exp(3.07e-4 * T);  // J/(kg·K)
}

// ---- temperature-dependent resistivity ----

// Eq. (6.9): rho_Cu(T) = -3.5e-9 + 7.2e-11 * T   ohm·m
double resistivity_copper(double T) {
    return -3.5e-9 + 7.2e-11 * T;  // ohm·m
}

// Eq. (6.9): rho_Al(T) = -6.57e-9 + 1.2e-10 * T   ohm·m
double resistivity_aluminum(double T) {
    return -6.57e-9 + 1.2e-10 * T;  // ohm·m
}

} // namespace coilgun::physics
