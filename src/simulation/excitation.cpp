/**
 * @file excitation.cpp
 * @brief Excitation model implementations (capacitor, crowbar, waveform).
 * @author Winston Meursault
 */

#include "coilgun/simulation/excitation.hpp"
#include <cmath>

namespace coilgun::simulation {

CapacitorExcitation::CapacitorExcitation(double iv, double cap) : C_(cap), U_C_(iv), U0_(iv) {}
double CapacitorExcitation::voltage() const { return U_C_; }
void CapacitorExcitation::advance(double dt, double I) {
    U_C_ -= I * dt / C_;
    if (U_C_ <= 0.0) { U_C_ = 0.0; finished_ = true; }
}
bool CapacitorExcitation::finished() const { return finished_; }
void CapacitorExcitation::reset() { U_C_ = U0_; finished_ = false; }
double CapacitorExcitation::capacitance() const { return C_; }
double CapacitorExcitation::capacitor_voltage() const { return U_C_; }
double CapacitorExcitation::initial_voltage() const { return U0_; }

void CrowbarExcitation::advance(double dt, double I) {
    if (!diode_on_) {
        U_C_ -= I * dt / C_;
        if (U_C_ <= 0.0) { diode_on_ = true; U_C_ = 0.0; }
    }
    if (diode_on_ && std::abs(I) < 1e-6) finished_ = true;
}
void CrowbarExcitation::reset() {
    CapacitorExcitation::reset();
    diode_on_ = false;
}
bool CrowbarExcitation::diode_on() const { return diode_on_; }

WaveformExcitation::WaveformExcitation(Waveform func) : func_(std::move(func)) {}
void WaveformExcitation::set_end_time(double t) { end_time_ = t; }
double WaveformExcitation::voltage() const { return func_(t_); }
void WaveformExcitation::advance(double dt, double) {
    t_ += dt;
    if (t_ >= end_time_) finished_ = true;
}
bool WaveformExcitation::finished() const { return finished_; }
void WaveformExcitation::reset() { t_ = 0.0; finished_ = false; }

} // namespace coilgun::simulation
