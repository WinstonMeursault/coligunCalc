/**
 * @file excitation.hpp
 * @brief Excitation models (capacitor discharge, crowbar diode, arbitrary waveform).
 * @author Winston Meursault
 *
 * @see NumericalModel Sec.3.1, Sec.3.4.
 */

#pragma once

#include "coilgun/simulation/excitation_snapshot.hpp"

#include <functional>
#include <memory>
#include <cmath>

namespace coilgun::simulation {

/**
 * @brief Abstract base class for excitation sources.
 *
 * Represents the voltage source driving a coil stage. Subclasses implement
 * capacitor discharge (with/without crowbar diode) and arbitrary waveform
 * excitation.
 */
class Excitation {
public:
    virtual ~Excitation() = default;

    /// @brief Current terminal voltage, V.
    virtual double voltage() const = 0;

    /**
     * @brief Advance internal state by one time step.
     * @param dt Time step, s.
     * @param coil_current Current drawn from the source, A.
     */
    virtual void advance(double dt, double coil_current) = 0;

    /// @brief Whether the excitation has completed its useful output.
    virtual bool finished() const = 0;

    /// @brief Restore initial state (voltage restored, timer reset).
    virtual void reset() {}

    virtual std::unique_ptr<ExcitationSnapshot> snapshot() const = 0;
    virtual void restore(const ExcitationSnapshot& snapshot) = 0;
    virtual double voltage(const ExcitationSnapshot& snapshot) const = 0;
    virtual ExcitationDerivative continuous_derivative(
        const ExcitationSnapshot& snapshot, double coil_current) const = 0;
    virtual void advance_snapshot(ExcitationSnapshot& snapshot,
                                  double dt, double coil_current) const = 0;
    virtual void advance_snapshot_derivative(
        ExcitationSnapshot& snapshot, double dt,
        const ExcitationDerivative& derivative) const = 0;
    virtual void apply_event(ExcitationSnapshot& snapshot,
                             ExcitationEvent event) const = 0;
};

/**
 * @brief Capacitor discharge excitation (no crowbar diode).
 *
 * Models a simple series RLC discharge. The capacitor voltage decays
 * as @f$ U_C' = -I / C @f$. Finished once @f$ U_C \le 0 @f$.
 */
class CapacitorExcitation : public Excitation {
public:
    /**
     * @brief Construct a capacitor discharge source.
     * @param initial_voltage Initial capacitor voltage, V.
     * @param capacitance Capacitance, F.
     */
    CapacitorExcitation(double initial_voltage, double capacitance);

    double voltage() const override;
    void advance(double dt, double coil_current) override;
    bool finished() const override;
    void reset() override;
    std::unique_ptr<ExcitationSnapshot> snapshot() const override;
    void restore(const ExcitationSnapshot& snapshot) override;
    double voltage(const ExcitationSnapshot& snapshot) const override;
    ExcitationDerivative continuous_derivative(
        const ExcitationSnapshot& snapshot, double coil_current) const override;
    void advance_snapshot(ExcitationSnapshot& snapshot,
                          double dt, double coil_current) const override;
    void advance_snapshot_derivative(
        ExcitationSnapshot& snapshot, double dt,
        const ExcitationDerivative& derivative) const override;
    void apply_event(ExcitationSnapshot& snapshot,
                     ExcitationEvent event) const override;

    /// @brief Capacitance, F.
    double capacitance() const;
    /// @brief Current capacitor voltage, V.
    double capacitor_voltage() const;
    /// @brief Initial capacitor voltage, V.
    double initial_voltage() const;

protected:
    double C_;
    double U_C_;
    double U0_;
    bool finished_ = false;
};

/**
 * @brief Capacitor discharge with crowbar (freewheeling) diode.
 *
 * When the capacitor voltage reaches zero, the diode turns on and
 * the circuit transitions from RLC discharge to RL decay.
 *
 * @see NumericalModel Sec.3.4.
 */
class CrowbarExcitation : public CapacitorExcitation {
public:
    using CapacitorExcitation::CapacitorExcitation;

    double voltage() const override { return diode_on_ ? 0.0 : U_C_; }
    void advance(double dt, double coil_current) override;
    void reset() override;
    std::unique_ptr<ExcitationSnapshot> snapshot() const override;
    void restore(const ExcitationSnapshot& snapshot) override;
    double voltage(const ExcitationSnapshot& snapshot) const override;
    ExcitationDerivative continuous_derivative(
        const ExcitationSnapshot& snapshot, double coil_current) const override;
    void advance_snapshot(ExcitationSnapshot& snapshot,
                          double dt, double coil_current) const override;
    void advance_snapshot_derivative(
        ExcitationSnapshot& snapshot, double dt,
        const ExcitationDerivative& derivative) const override;
    void apply_event(ExcitationSnapshot& snapshot,
                     ExcitationEvent event) const override;

    /// @brief Whether the crowbar diode is currently conducting.
    bool diode_on() const;

private:
    bool diode_on_ = false;
};

/**
 * @brief Arbitrary user-defined voltage waveform @f$ V(t) @f$.
 *
 * Accepts a std::function<double(double)> mapping time to voltage.
 * Optional end-time early termination is supported.
 */
class WaveformExcitation : public Excitation {
public:
    /// @brief Voltage function type: @f$ f: t \mapsto V @f$.
    using Waveform = std::function<double(double)>;

    /**
     * @brief Construct from a user-supplied waveform.
     * @param func Function returning voltage given time, V = func(t).
     */
    explicit WaveformExcitation(Waveform func);

    /**
     * @brief Set an optional end time at which excitation is considered finished.
     * @param t End time, s.
     */
    void set_end_time(double t);

    double voltage() const override;
    void advance(double dt, double coil_current) override;
    bool finished() const override;
    void reset() override;
    std::unique_ptr<ExcitationSnapshot> snapshot() const override;
    void restore(const ExcitationSnapshot& snapshot) override;
    double voltage(const ExcitationSnapshot& snapshot) const override;
    ExcitationDerivative continuous_derivative(
        const ExcitationSnapshot& snapshot, double coil_current) const override;
    void advance_snapshot(ExcitationSnapshot& snapshot,
                          double dt, double coil_current) const override;
    void advance_snapshot_derivative(
        ExcitationSnapshot& snapshot, double dt,
        const ExcitationDerivative& derivative) const override;
    void apply_event(ExcitationSnapshot& snapshot,
                     ExcitationEvent event) const override;

    double time() const { return t_; }
    double end_time() const { return end_time_; }

private:
    Waveform func_;
    double t_ = 0.0;
    double end_time_ = INFINITY;
    bool finished_ = false;
};

} // namespace coilgun::simulation
