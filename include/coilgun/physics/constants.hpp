/**
 * @file constants.hpp
 * @brief Physical constants and material properties for coilgun simulation.
 * @author Winston Meursault
 */

#pragma once

#include <cmath>
#include <stdexcept>

namespace coilgun::physics {

/**
 * @brief Vacuum permeability, H/m.
 */
constexpr double MU0 = 4.0 * M_PI * 1e-7;

/**
 * @brief Reference temperature for material properties, K.
 */
constexpr double T_REFERENCE = 293.0;

/**
 * @brief Armature material identifier.
 *
 * Used to select temperature-dependent specific heat capacity and resistivity
 * coefficient at runtime.  To add a new material:
 *   1. Add an entry to this enum.
 *   2. Implement the corresponding @c specific_heat_capacity_X(T) and
 *      @c resistivity_X(T) functions.
 *   3. Add cases in @c material_cp() and @c material_beta().
 */
enum class ArmatureMaterial {
    Aluminum,
    Copper
};

/**
 * @brief Temperature-dependent material properties for a conductor.
 */
struct MaterialProperties {
    double resistivity_ref;   ///< Resistivity at T_REFERENCE, ohm·m.
    double temp_coefficient;  ///< Temperature coefficient beta, K^{-1}.
    double density;           ///< Mass density, kg/m^3.
};

extern const MaterialProperties COPPER;
extern const MaterialProperties ALUMINUM;

/**
 * @brief Specific heat capacity for a given armature material, J/(kg·K).
 * @param material Armature material.
 * @param temperature_kelvin Temperature in kelvin.
 * @return Specific heat capacity at constant pressure.
 */
double material_cp(ArmatureMaterial material, double temperature_kelvin);

/**
 * @brief Temperature coefficient of resistivity for a given armature material.
 * @param material Armature material.
 * @return Resistivity temperature coefficient beta, K^{-1}.
 */
double material_beta(ArmatureMaterial material);

/**
 * @brief Specific heat capacity of copper, J/(kg·K).
 * @param temperature_kelvin Temperature in kelvin.
 * @return Specific heat capacity at constant pressure.
 *
 * Eq. (6.3) in NumericalModel.md — converted from kJ to J by factor 1000.
 */
double specific_heat_capacity_copper(double temperature_kelvin);

/**
 * @brief Specific heat capacity of aluminum, J/(kg·K).
 * @param temperature_kelvin Temperature in kelvin.
 * @return Specific heat capacity at constant pressure.
 */
double specific_heat_capacity_aluminum(double temperature_kelvin);

/**
 * @brief Electrical resistivity of copper, ohm·m.
 * @param temperature_kelvin Temperature in kelvin.
 * @return Resistivity at the given temperature.
 *
 * Eq. (6.9) in NumericalModel.md.
 */
double resistivity_copper(double temperature_kelvin);

/**
 * @brief Electrical resistivity of aluminum, ohm·m.
 * @param temperature_kelvin Temperature in kelvin.
 * @return Resistivity at the given temperature.
 *
 * Eq. (6.9) in NumericalModel.md.
 */
double resistivity_aluminum(double temperature_kelvin);

} // namespace coilgun::physics
