/**
 * @file driving_coil.hpp
 * @brief A single multi-turn helical driving coil (stator).
 * @author Winston Meursault
 */

#pragma once

#include <optional>

namespace coilgun::components {

/**
 * @brief A single multi-turn helical driving coil (stator).
 *
 * Stores geometry and precomputes DC resistance and self-inductance
 * per NumericalModel Sec.8.2.
 */
class DrivingCoil {
public:
    /**
     * @brief Construct a driving coil from geometry and material properties.
     *
     * @param inner_radius  Inner winding radius (m).
     * @param outer_radius  Outer winding radius (m).
     * @param length        Axial length of the winding (m).
     * @param turns         Total number of turns.
     * @param resistivity   Wire material resistivity (ohm·m).
     * @param wire_area     Cross-sectional area of the wire (m^2).
     * @param fill_factor   Winding fill factor (0..1).
     * @param position      Initial centre position (m), default @c length/2.
     * @param force_exact_self_inductance If true, use reference-grade
     *        self-inductance computation (bypasses T(q,p) lookup table).
     */
    DrivingCoil(double inner_radius, double outer_radius, double length,
                int turns, double resistivity, double wire_area,
                 double fill_factor,
                 std::optional<double> position = std::nullopt,
                bool force_exact_self_inductance = false);

    /// @brief Inner winding radius (m).
    double inner_radius() const;
    /// @brief Outer winding radius (m).
    double outer_radius() const;
    /// @brief Axial winding length (m).
    double length() const;
    /// @brief Mean winding radius = (ri + re) / 2 (m).
    double mean_radius() const;
    /// @brief Total number of turns.
    int    turns() const;

    /// @brief Turns density = N / ((re - ri) * l) (turns/m^2).
    double turns_density() const;
    /// @brief DC resistance (ohm), precomputed at construction.
    double resistance() const;
    /// @brief Self-inductance (H), precomputed at construction.
    double self_inductance() const;

    /// @brief Current axial position of the coil centre (m).
    double position() const;
    /// @brief Set the axial position of the coil centre.
    /// @param x New position (m).
    void   set_position(double x);

private:
    double ri_, re_, l_;
    int    n_;
    double rho_, wire_area_, k_fill_;
    double x_;
    double nc_;
    double R_;
    double L_;
};

} // namespace coilgun::components
