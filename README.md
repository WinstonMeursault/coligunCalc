# coligunCalc

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language](https://img.shields.io/badge/language-C%2B%2B20-blue)](https://isocpp.org/)
[![Status](https://img.shields.io/badge/status-developing-orange)]()

[中文文档](README_cn.md)

---

A multi-stage synchronous induction coilgun simulation library using the **current filament method** (CFM). Models the electromagnetic coupling between driving coils and the armature via numerical integration of complete elliptic integrals, producing Lorentz force, velocity, and efficiency predictions.

Licensed under GPLv3.

---

## Key Features

- **Current filament method** — armature discretised into m × n concentric ring filaments capturing the skin effect
- **Closed-form filament kernel** — mutual inductance via complete elliptic integrals K(k) / E(k) (Maxwell's formula)
- **4D Gauss-Legendre quadrature** — coil-to-coil mutual inductance with configurable node count (n⁴ evaluations)
- **Self-inductance dual-track** — fast T(q,p) lookup table (~µs) with automatic fallback to reference-grade Bessel/Struve integration (~ms)
- **Struve H₀/H₁ implementation** — SciPy-equivalent three-regime strategy (power series / asymptotic / cross-checked)
- **Crowbar (freewheeling) diode** — per-stage autonomous diode state management eliminating braking force
- **Configurable excitation** — capacitor discharge, crowbar diode, or arbitrary waveform V(t)
- **Two time steppers** — forward Euler (fast) and classic RK4 (accurate)
- **Optional thermal coupling** — adiabatic ohmic heating with temperature-dependent cp(T) and ρ(T) for Cu and Al
- **Three optimization levels** — Reference / LookupTable / Full (distance cutoff + adaptive GL order)
- **Stage triggering** — position-based or time-delay, up to 50 stages
- **LRU caching** — 4096-entry filament-level M and dM/dz caches
- **OpenMP multicore support** — the CPU M/dM loop is parallelized when the filament count is at least 8; force and thermal update loops are currently serial
- **GPU acceleration (CUDA)** — synchronous Euler method surfaces align with the CPU API; stepping GPU single-/multi-stage RK4 instantiations currently throws `std::logic_error` and requires the CPU implementation

---

## Project Structure

```
include/coilgun/   — Public headers (core types, physics, components, simulation)
src/               — Library implementation (static library libcoilgun.a)
tests/             — Unit and integration tests (doctest/CTest)
tools/             — T(q,p) lookup table generator
docs/              — Documentation (API reference EN/CN, numerical model, design docs)
.references/       — Reference papers in PDF (gitignored, local-only)
```

## Quick Start

**Requirements**: C++20 compiler with OpenMP support, CMake ≥ 3.20.

Both Boost.Math and Eigen are fetched automatically via CMake `FetchContent` — no manual installation needed.

```sh
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
```

### GPU Acceleration (CUDA)

**Additional requirements**: CUDA Toolkit ≥ 12.8 (for Blackwell sm_120; ≥ 9.0 for earlier architectures), NVIDIA GPU with Compute Capability ≥ 6.0.

```sh
cmake --preset cuda-debug
cmake --build --preset cuda-debug
ctest --preset cuda-debug
ctest --preset cuda-debug -L gpu
```

The measurement-only CUDA benchmark includes CPU Reference rows and records
Direct, Graph, Persistent-request, Fallback, batch-size, solver, and thermal
timings. It is excluded from the default build because results are
machine-dependent:

```sh
cmake --build --preset cuda-debug --target bench_gpu_engine
./build/cuda-debug/src/cuda/bench_gpu_engine
```

See [the benchmark record](docs/benchmarks/2026-07-19-unified-gpu-engine.md)
for the RTX 5080 Laptop measurement and the honest Persistent fallback status.

The synchronous GPU method surface aligns with the CPU API for Euler. Stepping GPU single-/multi-stage RK4 instantiations currently throws `std::logic_error`; use the CPU implementation for RK4.

```cpp
#include <coilgun/coilgun_cuda.hpp>

using coilgun::simulation::cuda::GpuSingleStageSim;
// or: coilgun::simulation::cuda::GpuMultiStageSim

GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt);
sim.run();
double v = sim.result().summary.muzzle_velocity;
```

The CUDA engine keeps the fixed-shape Euler state, matrix/RHS assembly, batched solve, force/motion update, optional thermal update, and compact device control in resident buffers. `Graph` captures and replays that complete fixed-shape device step; `Direct` launches the same resident stages directly. Runtime initialization, allocation, solver, capture, replay, or validation failures restore the pre-step state, record `FallbackReason::RuntimeFailure`, and lock the engine to the CPU/Eigen fallback. `Persistent` remains an explicit safe fallback until a dedicated resident control stream is implemented.

GPU tests may be skipped when no CUDA device is available. A GPU validation run must additionally require `ExecutionReport::gpu_executed == true` and `ExecutionReport::backend != BackendMode::Fallback`; a requested backend or a passing CPU fallback is not proof of GPU execution.

### Minimal Usage Example

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::SingleStageSim;
    using coilgun::simulation::CrowbarExcitation;

    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6);

    sim.run();
    auto& r = sim.result();
    std::cout << "Velocity: " << r.summary.muzzle_velocity << " m/s\n"
              << "Efficiency: " << r.summary.efficiency * 100 << " %\n";
}
```

### Link Against the Library

```cmake
# In your CMakeLists.txt
find_package(OpenMP REQUIRED)
target_link_libraries(your_target PRIVATE coilgun OpenMP::OpenMP_CXX)
```

Or compile directly:

```sh
g++ -std=c++20 -fopenmp -Iinclude your_file.cpp build/src/libcoilgun.a -o your_binary
```

---

## Physics Foundation

The library implements the numerical model described in [NumericalModel.md](docs/NumericalModel.md):

| Module | What it computes | Key formula |
|--------|-----------------|-------------|
| Self-inductance | L of hollow cylindrical coils | L = μ₀ · nc² · ri³ · T(q,p) |
| Mutual inductance (filament) | M between two coaxial loops | Maxwell's elliptic integral formula |
| Mutual inductance (coil) | M between two thick coils | 4D Gauss-Legendre quadrature |
| Circuit ODE | dI/dt for coupled coil+armature system | [dI/dt] = (L-M)^(-1) · (U + v·dM·I - R·I) |
| Force | Axial Lorentz force | F = Σ dM/dx · I_d · I_p |
| Thermal | Adiabatic ohmic heating | ΔT = I²RΔt / (m·c_p(T)) |

---

## Test Suite

CTest discovers the current target set from the selected preset; the documentation does not freeze a target count. Use labels to separate quick CPU checks, numerical validation, and CUDA tests:

```sh
ctest --preset cpu-debug -L fast
ctest --preset cpu-debug -L validation
ctest --preset cuda-debug -L gpu
```

GPU tests use a shared resource lock. Tests that require a physical CUDA device skip when no device is available; actual GPU CI must also assert the execution-report conditions above. The benchmark executable `bench_gpu_engine` is not a CTest target and must be run explicitly.

---

## Documentation

- [API Reference (EN)](docs/API.md) — Complete C++ function and class API
- [API 参考 (中文)](docs/API_cn.md) — Chinese translation
- [Numerical Model](docs/NumericalModel.md) — Detailed physics derivation and algorithm
- [GPU Benchmark Record](docs/benchmarks/2026-07-19-unified-gpu-engine.md) — Machine-specific execution measurements and fallback status

---

## CMake Presets

| Preset | Generator | Build type | Flags | Notes |
|--------|-----------|------------|-------|-------|
| `cpu-debug` | Ninja | Debug | `-march=native`, `-O2 -g` | CPU only, tests ON, compile_commands.json |
| `cpu-release` | Ninja | Release | `-march=native` | CPU only, tests ON |
| `cuda-debug` | Ninja | Debug | `-march=native`, `-O2 -g` | CUDA ON, tests ON, compile_commands.json |
| `cuda-release` | Ninja | Release | `-march=native` | CUDA ON, tests ON |

---

## License

GNU General Public License v3.0. See [LICENSE](LICENSE) for the full text.
