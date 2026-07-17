/**
 * @file armature.hpp
 * @brief Solid cylindrical armature discretised into current filaments.
 * @author Winston Meursault
 */

#pragma once

#include "coilgun/physics/constants.hpp"

#include <vector>

namespace coilgun::components {

/**
 * @brief Solid cylindrical armature discretised into @em m &times; @em n
 *        current filaments.
 *
 * @em m = axial slices, @em n = radial layers.  Total filaments = m * n.
 *
 * @par Filament indexing convention (1-indexed, row-major)
 *      @em i = 1&hellip;m (axial), @em j = 1&hellip;n (radial).\n
 *      Internal storage uses a flat array [0 &hellip; m&middot;n&minus;1]
 *      in row-major order:\n
 *      @f$ \mathrm{index} = (i-1) \cdot n + (j-1) @f$
 *
 * @see NumericalModel Sec.2.3, Sec.8.2 items 4&ndash;5.
 */
class Armature {
public:
    /**
     * @brief Construct an armature with geometry, material, motion state,
     *        and discretisation parameters.
     *
     * @param inner_radius    Inner bore radius (m).
     * @param outer_radius    Outer radius (m).
     * @param length          Axial length (m).
     * @param resistivity     Material resistivity (ohm&middot;m) at initial
     *                        temperature.
     * @param material_density Material mass density (kg/m^3) for filament
     *                         masses.
     * @param velocity        Initial velocity (m/s).
     * @param mass            Total armature mass (kg) &mdash; includes payload.
     * @param m_axial         Number of axial filament divisions (@em m).
     * @param n_radial        Number of radial filament divisions (@em n).
     * @param position        Initial centre position (m).
     * @param material        Armature material (determines cp(T) and beta
     *                        for thermal mode).
     * @param force_exact_self_inductance If true, use reference-grade
     *        filament self-inductance (bypasses T(q,p) lookup table).
     *
     * @see NumericalModel Sec.2.3.
     */
    Armature(double inner_radius, double outer_radius, double length,
             double resistivity, double material_density,
             double velocity, double mass,
             int m_axial, int n_radial, double position,
             physics::ArmatureMaterial material = physics::ArmatureMaterial::Aluminum,
             bool force_exact_self_inductance = false);

    // ---- Geometry queries ----

    /// @brief Inner bore radius (m).
    double inner_radius() const;
    /// @brief Outer radius (m).
    double outer_radius() const;
    /// @brief Axial length (m).
    double length() const;
    /// @brief Number of axial divisions @em m.
    int    axial_filaments() const;
    /// @brief Number of radial divisions @em n.
    int    radial_filaments() const;
    /// @brief Total number of filaments = @em m &times; @em n.
    int    total_filaments() const;

    // ---- Filament geometry queries (1-indexed) ----
    // @see NumericalModel Sec.2.3.

    /**
     * @brief Inner radius of the @em j-th radial filament layer.
     * @param j Radial layer index (1..n).
     * @return Radius at the start of the layer (m).
     */
    double filament_inner_radius(int j) const;

    /**
     * @brief Outer radius of the @em j-th radial filament layer.
     * @param j Radial layer index (1..n).
     * @return Radius at the end of the layer (m).
     */
    double filament_outer_radius(int j) const;

    /**
     * @brief Mean (centre-line) radius of the @em j-th radial layer.
     * @param j Radial layer index (1..n).
     * @return @f$ r_i + (j - 0.5) \cdot \Delta r @f$ (m).
     */
    double filament_mean_radius(int j) const;

    /**
     * @brief Axial position of the @em i-th filament ring centre.
     * @param i Axial layer index (1..m).
     * @return @f$ x - l/2 + (i - 0.5) \cdot \Delta l @f$ (m).
     */
    double filament_axial_position(int i) const;

    // ---- Per-filament data (flat arrays, length m*n) ----
    //
    // Order: filament (i=1,j=1), (i=1,j=2), ..., (i=m,j=n) —
    // i.e. row-major with i as the slowest-varying index.
    //
    // @see NumericalModel Sec.8.2 items 4–5.

    /// @brief Resistance per filament (ohm).
    ///        Row-major order.
    const std::vector<double>& resistances() const;

    /// @brief Self-inductance per filament (H).
    ///        Row-major order.
    const std::vector<double>& inductances() const;

    /// @brief Mass per filament (kg).
    ///        Row-major order.
    const std::vector<double>& masses() const;

    // ---- Motion state ----

    /// @brief Current axial centre position (m).
    double position() const;
    /// @brief Current velocity (m/s).
    double velocity() const;
    /// @brief Total armature mass including payload (kg).
    double mass() const;

    /// @brief Armature material (for thermal dispatch).
    physics::ArmatureMaterial material() const;

    /**
     * @brief Translate the armature by @p dx metres.
     * @param dx Displacement (m).
     */
    void update_position(double dx);

    /**
     * @brief Set the velocity to a new value.
     * @param v New velocity (m/s).
     */
    void set_velocity(double v);

private:
    double ri_, re_, l_;
    double rho_;
    int    m_, n_;
    double x_, v_, ma_;
    physics::ArmatureMaterial material_ = physics::ArmatureMaterial::Aluminum;

    double dr_;
    double dl_;
    double nc_fil_;

    std::vector<double> R_;
    std::vector<double> L_;
    std::vector<double> mass_;
};

} // namespace coilgun::components
