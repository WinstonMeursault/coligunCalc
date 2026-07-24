#pragma once

#include <array>
#include <cstddef>
#include <cstdint>

#ifndef COILGUN_ENABLE_CPU_PHASE_TIMING
#define COILGUN_ENABLE_CPU_PHASE_TIMING 0
#endif

namespace coilgun::simulation {

enum class CpuPhase : std::uint8_t {
    Mutual = 0,
    Assembly,
    Solve,
    Thermal,
    Orchestration,
};

struct CpuPhaseTiming {
    std::uint64_t derivative_count = 0;
    std::array<std::uint64_t, 5> phase_nanoseconds{};

    void reset() noexcept;
    std::uint64_t nanoseconds(CpuPhase phase) const noexcept;
    double milliseconds(CpuPhase phase) const noexcept;
};

#if COILGUN_ENABLE_CPU_PHASE_TIMING

class CpuPhaseTimingCollector {
public:
    explicit CpuPhaseTimingCollector(CpuPhaseTiming& timing) noexcept;
    ~CpuPhaseTimingCollector() noexcept;

    CpuPhaseTimingCollector(const CpuPhaseTimingCollector&) = delete;
    CpuPhaseTimingCollector& operator=(const CpuPhaseTimingCollector&) = delete;

private:
    CpuPhaseTiming* previous_ = nullptr;
};

class CpuPhaseTimingScope {
public:
    explicit CpuPhaseTimingScope(CpuPhase phase) noexcept;
    ~CpuPhaseTimingScope() noexcept;

    CpuPhaseTimingScope(const CpuPhaseTimingScope&) = delete;
    CpuPhaseTimingScope& operator=(const CpuPhaseTimingScope&) = delete;

private:
    CpuPhaseTiming* timing_ = nullptr;
    CpuPhase phase_;
    std::uint64_t start_ticks_ = 0;
};

#else

class CpuPhaseTimingCollector {
public:
    explicit constexpr CpuPhaseTimingCollector(CpuPhaseTiming&) noexcept {}
    constexpr ~CpuPhaseTimingCollector() noexcept = default;

    CpuPhaseTimingCollector(const CpuPhaseTimingCollector&) = delete;
    CpuPhaseTimingCollector& operator=(const CpuPhaseTimingCollector&) = delete;
};

class CpuPhaseTimingScope {
public:
    explicit constexpr CpuPhaseTimingScope(CpuPhase) noexcept {}
    constexpr ~CpuPhaseTimingScope() noexcept = default;

    CpuPhaseTimingScope(const CpuPhaseTimingScope&) = delete;
    CpuPhaseTimingScope& operator=(const CpuPhaseTimingScope&) = delete;
};

#endif

} // namespace coilgun::simulation
