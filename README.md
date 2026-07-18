# coligunCalc

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language](https://img.shields.io/badge/language-C%2B%2B17-blue)](https://isocpp.org/)
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
- **OpenMP multicore parallelism** — parallelized M/dM computation, force summation, and thermal updates per time step
- **GPU acceleration (CUDA)** — mutual inductance 4D integration offloaded to GPU via CUDA 13.3+, with seamless API compatibility

---

## Project Structure

```
include/coilgun/   — Public headers (core types, physics, components, simulation)
src/               — Library implementation (static library libcoilgun.a)
tests/             — Unit and integration tests (doctest, 17 suites)
tools/             — T(q,p) lookup table generator
docs/              — Documentation (API reference EN/CN, numerical model, design docs)
.references/       — Reference papers in PDF (gitignored, local-only)
```

## Quick Start

**Requirements**: C++17 compiler with OpenMP support, CMake ≥ 3.20.

Both Boost.Math and Eigen are fetched automatically via CMake `FetchContent` — no manual installation needed.

```sh
cmake --preset ninja-debug       # configure
cmake --build --preset ninja-debug   # build
ctest --preset debug              # run tests
```

### GPU Acceleration (CUDA)

**Additional requirements**: CUDA Toolkit ≥ 12.8 (for Blackwell sm_120; ≥ 9.0 for earlier architectures), NVIDIA GPU with Compute Capability ≥ 6.0.

```sh
cmake --preset ninja-cuda-debug      # configure with CUDA
cmake --build --preset ninja-cuda-debug  # build CPU + GPU libraries
ctest --preset debug                 # run all tests (CPU + GPU)
```

GPU usage is API-compatible with CPU:

```cpp
#include <coilgun/coilgun_cuda.hpp>

using coilgun::simulation::cuda::GpuSingleStageSim;
// or: coilgun::simulation::cuda::GpuMultiStageSim

GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt);
sim.run();
double v = sim.result().summary.muzzle_velocity;
```

The GPU backend accelerates only the 4D mutual inductance integration (>95% of runtime). Linear system solve, kinematics, and thermal updates remain on CPU.

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
g++ -std=c++17 -fopenmp -Iinclude your_file.cpp build/src/libcoilgun.a -o your_binary
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

## Test Suite (17 suites, all passing)

| Suite | Coverage |
|-------|----------|
| `test_elliptic` | K(k), E(k) against known values |
| `test_struve` | H₀(x), H₁(x) against SciPy reference |
| `test_quadrature` | Gauss-Legendre/Laguerre node generation |
| `test_lookup` | T(q,p) table + bilinear interpolation vs reference |
| `test_self_inductance` | Table vs exact, edge cases, exact/comparison |
| `test_mutual_inductance` | Filament & coil-level M and dM/dz |
| `test_driving_coil` | Geometry, R, L, turns density |
| `test_armature` | Filament geometry, per-filament R/L/mass |
| `test_single_stage_sim` | Capacitor, crowbar, Euler/RK4, thermal |
| `test_multi_stage_sim` | 1-stage equivalence, 2-stage, triggers, optimization levels |
| `test_integration` | End-to-end single-stage scenarios |
| `test_gpu_elliptic` | GPU K(k), E(k) vs CPU reference |
| `test_gpu_filament` | GPU filament-level M, dM/dz vs CPU |
| `test_gpu_coil_pair` | GPU coil-filament M, dM/dz vs CPU |
| `test_gpu_vs_cpu_single` | GPU single-stage end-to-end vs CPU (ε < 0.5%) |
| `test_gpu_vs_cpu_multi` | GPU multi-stage end-to-end vs CPU (ε < 0.5%) |
| `test_gpu_batch` | GPU batch simulation mode |

---

## Documentation

- [API Reference (EN)](docs/API.md) — Complete C++ function and class API
- [API 参考 (中文)](docs/API_cn.md) — Chinese translation
- [Numerical Model](docs/NumericalModel.md) — Detailed physics derivation and algorithm
- [Multi-Stage Design](docs/multi_stage_sim_design.md) — Architecture and implementation plan
- [82mm Test Dataset](docs/test_dataset_82mm_coilgun.md) — Validation dataset from Xiang (2015)

---

## CMake Presets

| Preset | Generator | Build type | Flags | Notes |
|--------|-----------|------------|-------|-------|
| `ninja-debug` | Ninja | Debug | `-march=native` | Tests ON, compile_commands.json |
| `ninja-release` | Ninja | Release | `-march=native -O3` | — |
| `make-debug` | Unix Makefiles | Debug | `-march=native` | Tests ON, compile_commands.json |
| `ninja-cuda-debug` | Ninja | Debug | `-march=native` | CUDA ON, Tests ON |

---

## License

GNU General Public License v3.0. See [LICENSE](LICENSE) for the full text.
