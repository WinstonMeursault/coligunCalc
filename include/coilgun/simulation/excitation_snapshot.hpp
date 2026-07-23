#pragma once

#include <memory>

namespace coilgun::simulation {

class ExcitationSnapshot {
public:
    virtual ~ExcitationSnapshot() = default;
    virtual std::unique_ptr<ExcitationSnapshot> clone() const = 0;
};

struct CapacitorSnapshot final : ExcitationSnapshot {
    CapacitorSnapshot() = default;
    CapacitorSnapshot(double voltage, bool is_finished)
        : capacitor_voltage(voltage), finished(is_finished) {}
    double capacitor_voltage = 0.0;
    bool finished = false;

    std::unique_ptr<ExcitationSnapshot> clone() const override {
        return std::make_unique<CapacitorSnapshot>(*this);
    }
};

struct CrowbarSnapshot final : ExcitationSnapshot {
    CrowbarSnapshot() = default;
    CrowbarSnapshot(double voltage, bool is_diode_on, bool is_finished)
        : capacitor_voltage(voltage), diode_on(is_diode_on), finished(is_finished) {}
    double capacitor_voltage = 0.0;
    bool diode_on = false;
    bool finished = false;

    std::unique_ptr<ExcitationSnapshot> clone() const override {
        return std::make_unique<CrowbarSnapshot>(*this);
    }
};

struct WaveformSnapshot final : ExcitationSnapshot {
    WaveformSnapshot() = default;
    WaveformSnapshot(double waveform_time, bool is_finished)
        : time(waveform_time), finished(is_finished) {}
    double time = 0.0;
    bool finished = false;

    std::unique_ptr<ExcitationSnapshot> clone() const override {
        return std::make_unique<WaveformSnapshot>(*this);
    }
};

enum class ExcitationEvent {
    CapacitorZero,
    CrowbarOn,
    WaveformEnd,
    Finished,
};

struct ExcitationDerivative {
    double capacitor_voltage_rate = 0.0;
    double waveform_time_rate = 0.0;
};

} // namespace coilgun::simulation
