#include "coilgun/simulation/cpu_phase_timing.hpp"

#include <chrono>

namespace coilgun::simulation {
namespace {

constexpr std::size_t phase_index(CpuPhase phase) noexcept {
    return static_cast<std::size_t>(phase);
}

#if COILGUN_ENABLE_CPU_PHASE_TIMING
thread_local CpuPhaseTiming* active_timing = nullptr;
#endif

} // namespace

void CpuPhaseTiming::reset() noexcept {
    derivative_count = 0;
    phase_nanoseconds.fill(0);
}

std::uint64_t CpuPhaseTiming::nanoseconds(CpuPhase phase) const noexcept {
    return phase_nanoseconds[phase_index(phase)];
}

double CpuPhaseTiming::milliseconds(CpuPhase phase) const noexcept {
    return static_cast<double>(nanoseconds(phase)) / 1.0e6;
}

#if COILGUN_ENABLE_CPU_PHASE_TIMING

CpuPhaseTimingCollector::CpuPhaseTimingCollector(CpuPhaseTiming& timing) noexcept
    : previous_(active_timing) {
    active_timing = &timing;
}

CpuPhaseTimingCollector::~CpuPhaseTimingCollector() noexcept {
    active_timing = previous_;
}

CpuPhaseTimingScope::CpuPhaseTimingScope(CpuPhase phase) noexcept
    : timing_(active_timing), phase_(phase) {
    if (timing_) {
    start_ticks_ = static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now().time_since_epoch()).count());
    }
}

CpuPhaseTimingScope::~CpuPhaseTimingScope() noexcept {
    if (!timing_) return;
    const auto end = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
    const auto elapsed = static_cast<std::uint64_t>(end) - start_ticks_;
    timing_->phase_nanoseconds[phase_index(phase_)] += elapsed;
    if (phase_ == CpuPhase::Orchestration) ++timing_->derivative_count;
}

#endif

} // namespace coilgun::simulation
