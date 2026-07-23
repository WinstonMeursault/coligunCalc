/**
 * @file excitation.cpp
 * @brief Excitation model implementations (capacitor, crowbar, waveform).
 * @author Winston Meursault
 */

#include "coilgun/simulation/excitation.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace coilgun::simulation {

namespace {

template<typename Snapshot>
const Snapshot& checked_snapshot(const ExcitationSnapshot& snapshot, const char* type) {
    const auto* value = dynamic_cast<const Snapshot*>(&snapshot);
    if (!value) throw std::invalid_argument(std::string("snapshot type mismatch for ") + type);
    return *value;
}

template<typename Snapshot>
Snapshot& checked_snapshot(ExcitationSnapshot& snapshot, const char* type) {
    auto* value = dynamic_cast<Snapshot*>(&snapshot);
    if (!value) throw std::invalid_argument(std::string("snapshot type mismatch for ") + type);
    return *value;
}

void validate_step(double dt, double current) {
    if (!std::isfinite(dt) || dt < 0.0) throw std::invalid_argument("dt must be finite and non-negative");
    if (!std::isfinite(current)) throw std::invalid_argument("coil_current must be finite");
}

} // namespace

CapacitorExcitation::CapacitorExcitation(double iv, double cap) : C_(cap), U_C_(iv), U0_(iv) {
    if (!std::isfinite(iv) || iv <= 0.0)
        throw std::invalid_argument("initial_voltage must be finite and positive");
    if (!std::isfinite(cap) || cap <= 0.0)
        throw std::invalid_argument("capacitance must be finite and positive");
}
double CapacitorExcitation::voltage() const { return U_C_; }
void CapacitorExcitation::advance(double dt, double I) {
    auto state = snapshot();
    advance_snapshot(*state, dt, I);
    restore(*state);
}
bool CapacitorExcitation::finished() const { return finished_; }
void CapacitorExcitation::reset() { U_C_ = U0_; finished_ = false; }
double CapacitorExcitation::capacitance() const { return C_; }
double CapacitorExcitation::capacitor_voltage() const { return U_C_; }
double CapacitorExcitation::initial_voltage() const { return U0_; }

std::unique_ptr<ExcitationSnapshot> CapacitorExcitation::snapshot() const {
    return std::make_unique<CapacitorSnapshot>(CapacitorSnapshot{U_C_, finished_});
}
void CapacitorExcitation::restore(const ExcitationSnapshot& state) {
    const auto& value = checked_snapshot<CapacitorSnapshot>(state, "capacitor");
    U_C_ = value.capacitor_voltage;
    finished_ = value.finished;
}
double CapacitorExcitation::voltage(const ExcitationSnapshot& state) const {
    return checked_snapshot<CapacitorSnapshot>(state, "capacitor").capacitor_voltage;
}
ExcitationDerivative CapacitorExcitation::continuous_derivative(
    const ExcitationSnapshot& state, double current) const {
    checked_snapshot<CapacitorSnapshot>(state, "capacitor");
    if (!std::isfinite(current)) throw std::invalid_argument("coil_current must be finite");
    return {-current / C_, 0.0};
}
void CapacitorExcitation::advance_snapshot_derivative(
    ExcitationSnapshot& state, double dt, const ExcitationDerivative& derivative) const {
    if (!std::isfinite(dt) || dt < 0.0)
        throw std::invalid_argument("dt must be finite and non-negative");
    auto& value = checked_snapshot<CapacitorSnapshot>(state, "capacitor");
    value.capacitor_voltage += dt * derivative.capacitor_voltage_rate;
}
void CapacitorExcitation::advance_snapshot(
    ExcitationSnapshot& state, double dt, double current) const {
    validate_step(dt, current);
    advance_snapshot_derivative(state, dt, continuous_derivative(state, current));
    auto& value = checked_snapshot<CapacitorSnapshot>(state, "capacitor");
    if (value.capacitor_voltage <= 0.0) {
        apply_event(state, ExcitationEvent::CapacitorZero);
        apply_event(state, ExcitationEvent::Finished);
    }
}
void CapacitorExcitation::apply_event(ExcitationSnapshot& state, ExcitationEvent event) const {
    auto& value = checked_snapshot<CapacitorSnapshot>(state, "capacitor");
    if (event == ExcitationEvent::CapacitorZero) value.capacitor_voltage = 0.0;
    if (event == ExcitationEvent::Finished) value.finished = true;
}

void CrowbarExcitation::advance(double dt, double I) {
    auto state = snapshot();
    advance_snapshot(*state, dt, I);
    restore(*state);
}
void CrowbarExcitation::reset() {
    CapacitorExcitation::reset();
    diode_on_ = false;
}
bool CrowbarExcitation::diode_on() const { return diode_on_; }

std::unique_ptr<ExcitationSnapshot> CrowbarExcitation::snapshot() const {
    return std::make_unique<CrowbarSnapshot>(CrowbarSnapshot{U_C_, diode_on_, finished_});
}
void CrowbarExcitation::restore(const ExcitationSnapshot& state) {
    const auto& value = checked_snapshot<CrowbarSnapshot>(state, "crowbar");
    U_C_ = value.capacitor_voltage;
    diode_on_ = value.diode_on;
    finished_ = value.finished;
}
double CrowbarExcitation::voltage(const ExcitationSnapshot& state) const {
    const auto& value = checked_snapshot<CrowbarSnapshot>(state, "crowbar");
    return value.diode_on ? 0.0 : value.capacitor_voltage;
}
ExcitationDerivative CrowbarExcitation::continuous_derivative(
    const ExcitationSnapshot& state, double current) const {
    const auto& value = checked_snapshot<CrowbarSnapshot>(state, "crowbar");
    if (!std::isfinite(current)) throw std::invalid_argument("coil_current must be finite");
    return {value.diode_on ? 0.0 : -current / C_, 0.0};
}
void CrowbarExcitation::advance_snapshot_derivative(
    ExcitationSnapshot& state, double dt, const ExcitationDerivative& derivative) const {
    if (!std::isfinite(dt) || dt < 0.0)
        throw std::invalid_argument("dt must be finite and non-negative");
    auto& value = checked_snapshot<CrowbarSnapshot>(state, "crowbar");
    value.capacitor_voltage += dt * derivative.capacitor_voltage_rate;
}
void CrowbarExcitation::advance_snapshot(
    ExcitationSnapshot& state, double dt, double current) const {
    validate_step(dt, current);
    advance_snapshot_derivative(state, dt, continuous_derivative(state, current));
    auto& value = checked_snapshot<CrowbarSnapshot>(state, "crowbar");
    if (!value.diode_on && value.capacitor_voltage <= 0.0) {
        apply_event(state, ExcitationEvent::CapacitorZero);
        apply_event(state, ExcitationEvent::CrowbarOn);
    }
    if (value.diode_on && std::abs(current) < 1e-6)
        apply_event(state, ExcitationEvent::Finished);
}
void CrowbarExcitation::apply_event(ExcitationSnapshot& state, ExcitationEvent event) const {
    auto& value = checked_snapshot<CrowbarSnapshot>(state, "crowbar");
    if (event == ExcitationEvent::CapacitorZero) value.capacitor_voltage = 0.0;
    if (event == ExcitationEvent::CrowbarOn) value.diode_on = true;
    if (event == ExcitationEvent::Finished) value.finished = true;
}

WaveformExcitation::WaveformExcitation(Waveform func) : func_(std::move(func)) {
    if (!func_) throw std::invalid_argument("waveform function must not be empty");
}
void WaveformExcitation::set_end_time(double t) {
    if (std::isnan(t) || t < 0.0 || t == -std::numeric_limits<double>::infinity())
        throw std::invalid_argument("waveform end_time must be non-negative or positive infinity");
    end_time_ = t;
}
double WaveformExcitation::voltage() const { return finished_ ? 0.0 : func_(t_); }
void WaveformExcitation::advance(double dt, double) {
    auto state = snapshot();
    advance_snapshot(*state, dt, 0.0);
    restore(*state);
}
bool WaveformExcitation::finished() const { return finished_; }
void WaveformExcitation::reset() { t_ = 0.0; finished_ = false; }

std::unique_ptr<ExcitationSnapshot> WaveformExcitation::snapshot() const {
    return std::make_unique<WaveformSnapshot>(WaveformSnapshot{t_, finished_});
}
void WaveformExcitation::restore(const ExcitationSnapshot& state) {
    const auto& value = checked_snapshot<WaveformSnapshot>(state, "waveform");
    t_ = value.time;
    finished_ = value.finished;
}
double WaveformExcitation::voltage(const ExcitationSnapshot& state) const {
    const auto& value = checked_snapshot<WaveformSnapshot>(state, "waveform");
    return value.finished ? 0.0 : func_(value.time);
}
ExcitationDerivative WaveformExcitation::continuous_derivative(
    const ExcitationSnapshot& state, double) const {
    const auto& value = checked_snapshot<WaveformSnapshot>(state, "waveform");
    return {0.0, value.finished ? 0.0 : 1.0};
}
void WaveformExcitation::advance_snapshot_derivative(
    ExcitationSnapshot& state, double dt, const ExcitationDerivative& derivative) const {
    if (!std::isfinite(dt) || dt < 0.0)
        throw std::invalid_argument("dt must be finite and non-negative");
    auto& value = checked_snapshot<WaveformSnapshot>(state, "waveform");
    value.time += dt * derivative.waveform_time_rate;
}
void WaveformExcitation::advance_snapshot(
    ExcitationSnapshot& state, double dt, double current) const {
    validate_step(dt, current);
    advance_snapshot_derivative(state, dt, continuous_derivative(state, current));
    auto& value = checked_snapshot<WaveformSnapshot>(state, "waveform");
    if (value.time >= end_time_) apply_event(state, ExcitationEvent::WaveformEnd);
}
void WaveformExcitation::apply_event(ExcitationSnapshot& state, ExcitationEvent event) const {
    auto& value = checked_snapshot<WaveformSnapshot>(state, "waveform");
    if (event == ExcitationEvent::WaveformEnd || event == ExcitationEvent::Finished) {
        if (std::isfinite(end_time_)) value.time = end_time_;
        value.finished = true;
    }
}

} // namespace coilgun::simulation
