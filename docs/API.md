# coligunCalc C++ API Reference

> Chinese version: [API_cn.md](API_cn.md)

---

## Getting Started

### Include and Link

The entire library lives under `namespace coilgun`. The simplest way to get everything:

```cpp
#include <coilgun/coilgun.hpp>   // all physics + components + simulation

// or cherry-pick what you need:
#include <coilgun/physics/constants.hpp>
#include <coilgun/components/driving_coil.hpp>
#include <coilgun/simulation/single_stage_sim.hpp>
// ...
```

Link against the static library:

```cmake
target_link_libraries(your_target PRIVATE coilgun)
```

For one-off scripts you can also compile directly against the static lib:

```sh
g++ -std=c++17 -Iinclude your_file.cpp build/src/libcoilgun.a -o your_binary
```

### Namespaces

| Namespace | Contents |
|-----------|----------|
| `coilgun::physics` | Physical constants, elliptic integrals, Struve functions, quadrature, self/mutual inductance, LRU cache, lookup tables |
| `coilgun::components` | DrivingCoil and Armature classes |
| `coilgun::simulation` | Simulation engine: time steppers, excitation models, termination, trigger config, SimState/MultiStageState, SingleStageSim, MultiStageSim |
| `coilgun::physics::detail` | Internal helpers (lookup table data) — do not rely on these |

### Minimal Example

Compute the self-inductance and DC resistance of a copper driving coil:

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;

    // 150-turn copper coil: ri=10mm, re=30mm, length=50mm, wire 1mm^2, fill 70%
    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);

    std::cout << "Self-inductance: " << coil.self_inductance() << " H\n";
    std::cout << "DC resistance:   " << coil.resistance() << " ohm\n";
    std::cout << "Turns density:   " << coil.turns_density() << " turns/m^2\n";
    return 0;
}
```

Compute mutual inductance between two coaxial filament loops:

```cpp
#include <coilgun/physics/mutual_inductance.hpp>
#include <iostream>

int main() {
    double M  = coilgun::physics::mutual_inductance_filament(0.02, 0.025, 0.01);
    double dM = coilgun::physics::mutual_inductance_gradient_filament(0.02, 0.025, 0.01);

    std::cout << "M  = " << M  << " H\n";
    std::cout << "dM = " << dM << " H/m\n";
    return 0;
}
```

---

## Architecture

```
include/coilgun/
├── core/types.hpp              — common types, forward declarations
├── physics/
│   ├── constants.hpp           — μ₀, material properties (Cu, Al), cp(T), rho(T)
│   ├── elliptic.hpp            — K(k), E(k) complete elliptic integrals
│   ├── struve.hpp              — H₀(x), H₁(x) Struve functions
│   ├── quadrature.hpp          — Gauss-Legendre & Gauss-Laguerre nodes
│   ├── cache.hpp               — LRU cache template (4096 entries default)
│   ├── lookup_tables.hpp       — T(q,p) inductance shape factor table
│   ├── lookup_table_data.hpp   — precomputed T(q,p) dense table (~2.9M entries)
│   ├── self_inductance.hpp     — self-inductance (exact + fast table)
│   └── mutual_inductance.hpp   — filament & coil-level mutual inductance
├── components/
│   ├── driving_coil.hpp        — DrivingCoil class
│   └── armature.hpp            — Armature class (m×n filament discretisation)
├── simulation/
│   ├── excitation.hpp          — Excitation, CapacitorExcitation, CrowbarExcitation, WaveformExcitation
│   ├── time_stepper.hpp        — EulerStepper, RK4Stepper
│   ├── sim_result.hpp          — SimStep, SimResult, SimSummary
│   ├── termination.hpp         — TerminationPolicy
│   ├── single_stage_sim.hpp    — SimState, SingleStageSim<StepperPolicy>
│   ├── trigger_config.hpp      — TriggerMode, TriggerConfig
│   ├── multi_stage_result.hpp  — StepSnapshot, MultiStageStep, PerStageSummary, MultiStageSummary, MultiStageResult
│   └── multi_stage_sim.hpp     — OptimizationLevel, MultiStageState, MultiStageSim<StepperPolicy>
└── coilgun.hpp                 — convenience umbrella header
```

---

## Physics Constants

```cpp
#include <coilgun/physics/constants.hpp>

namespace coilgun::physics {
    constexpr double MU0 = 4.0 * M_PI * 1e-7;    // vacuum permeability (H/m)
    constexpr double T_REFERENCE = 293.0;          // room temperature (K)

    struct MaterialProperties {
        double resistivity_ref;    // resistivity at T_REFERENCE (Ω·m)
        double temp_coefficient;   // temperature coefficient β (K⁻¹)
        double density;            // mass density (kg/m³)
    };

    extern const MaterialProperties COPPER;        // 1.75e-8 Ω·m, 0.0041 K⁻¹, 8960 kg/m³
    extern const MaterialProperties ALUMINUM;      // 2.82e-8 Ω·m, 0.0042 K⁻¹, 2700 kg/m³

    enum class ArmatureMaterial { Aluminum, Copper };

    double material_cp(ArmatureMaterial m, double T);    // specific_heat_capacity dispatch
    double material_beta(ArmatureMaterial m);             // resistivity temp_coefficient dispatch

    double specific_heat_capacity_copper(double T);       // J/(kg·K), Eq.6.3
    double specific_heat_capacity_aluminum(double T);     // J/(kg·K), Eq.6.2
    double resistivity_copper(double T);                  // Ω·m, Eq.6.9
    double resistivity_aluminum(double T);                // Ω·m, Eq.6.9
}
```

**Usage**:

```cpp
double rho_cu_at_373K = coilgun::physics::resistivity_copper(373.0);
// ~ 2.28e-8 Ω·m  (copper at 100°C)

// Passing material to a component:
DrivingCoil coil(..., COPPER.resistivity_ref, ...);
```

`material_cp` dispatches `specific_heat_capacity_copper` or `specific_heat_capacity_aluminum` based on the enum. `material_beta` dispatches the temperature coefficient.

---

## Elliptic Integrals

```cpp
#include <coilgun/physics/elliptic.hpp>

double coilgun::physics::elliptic_k(double m);    // complete elliptic integral K(m)
double coilgun::physics::elliptic_e(double m);    // complete elliptic integral E(m)
double coilgun::physics::elliptic_modulus(
    double radius_a, double radius_b,
    double separation);                           // modulus k for two coaxial loops
```

Wraps Boost.Math `ellint_1` / `ellint_2`. Uses the **parameter** convention (m = k²), **not** the modulus convention.

`elliptic_modulus` computes the geometric modulus for two coaxial circular loops at a given separation — the building block for all filament-level mutual inductance calculations (Eq. 4.7 in NumericalModel.md):

```
k = sqrt( 4·a·b / ((a + b)² + h²) )
```

where `a`, `b` are the two loop radii and `h` is the axial separation.

**Caution**: internally clamps `k` to avoid singularities at `k → 1` (coincident loops), so `elliptic_modulus(a, a, 0)` produces a finite result slightly less than 1.

**Quick checks** (known values):

```cpp
elliptic_k(0.0);   // M_PI / 2  ≈ 1.570796
elliptic_e(0.0);   // M_PI / 2  ≈ 1.570796
elliptic_k(0.5);   // ≈ 1.854075
elliptic_e(0.5);   // ≈ 1.350644
```

---

## Struve Functions

```cpp
#include <coilgun/physics/struve.hpp>

double coilgun::physics::struve_h0(double x);     // Struve H₀(x)
double coilgun::physics::struve_h1(double x);     // Struve H₁(x)
```

**Parity**: H₀ is odd (H₀(-x) = -H₀(x)), H₁ is even (H₁(-x) = H₁(x)).

SciPy three-regime strategy for |x|:
- **|x| < 8**: power series — fast, machine-precision accurate
- **|x| ≥ 20**: asymptotic expansion via modified Bessel K₀/K₁ integrals + Gauss-Laguerre quadrature
- **8 ≤ |x| < 20**: both methods computed; result cross-checked to 1e-6 relative tolerance

Verified against `scipy.special.struve`. Used internally by `self_inductance_exact`; you rarely need to call these directly.

---

## Quadrature

```cpp
#include <coilgun/physics/quadrature.hpp>

struct coilgun::physics::QuadratureNodes {
    std::vector<double> nodes;    // quadrature node positions
    std::vector<double> weights;  // quadrature weights
};

QuadratureNodes gauss_legendre(int n);    // n-point rule on [-1, 1]
QuadratureNodes gauss_laguerre(int n);    // n-point rule on [0, +∞)
```

Internally uses the Golub-Welsch algorithm (eigenvalues of the Jacobi matrix). Used by `self_inductance_exact`, `mutual_inductance_coil`, and the Struve asymptotic regime. Typical users won't need to call these directly.

---

## Self-Inductance

```cpp
#include <coilgun/physics/self_inductance.hpp>

// Primary API — auto-selects table or exact based on geometry

double coilgun::physics::self_inductance(
    double inner_radius, double outer_radius, double length,
    double turns_density,
    bool force_exact = false);    // (H) — table when in range (default), or skip table if true

// For when you explicitly need reference-grade accuracy

double coilgun::physics::self_inductance_exact(
    double inner_radius, double outer_radius, double length,
    double turns_density);                        // (H) — always Bessel/Struve kernel + GL16

// Shape factor (lower level, dimensionless)

double coilgun::physics::inductance_shape_factor(double q, double p);          // table lookup + bilinear interpolation
double coilgun::physics::inductance_shape_factor_reference(double q, double p); // exact via Bessel/Struve (slow)
```

The self-inductance of a hollow cylindrical coil is:

```
L = 2π · μ₀ · nc² · ri⁵ · T(q, p)
```

where `nc = N / ((outer - inner) × length)` is the turns density, `ri` is the inner radius, and `T(q, p)` is the dimensionless shape factor with `q = length / inner_radius`, `p = outer_radius / inner_radius`.

**`self_inductance` — the recommended default.**

When `force_exact` is false (default) and `(q, p)` fall within the table bounds (`q ∈ [0.05, 4]`, `p ∈ [1.05, 4]`), uses the fast T(q,p) lookup table (~μs, ~5e-7 relative error in T). Otherwise automatically falls back to the exact Bessel/Struve integration (~ms). When `force_exact` is true, the table is bypassed and the reference-grade integration is used unconditionally.

**`self_inductance_exact` — reference-grade guarantee.**

Always uses the full Bessel/Struve kernel with composite GL16 integration. Equivalent to calling `self_inductance` with `force_exact = true`.

| Function | Behaviour | Speed | Use case |
|----------|-----------|-------|----------|
| `self_inductance` | Table when in range, exact fallback | ~μs or ~ms | **Default** — simulation, components |
| `self_inductance_exact` | Always exact | ~ms | Verification, benchmark data |

The `DrivingCoil` and `Armature` classes use `self_inductance` internally.

**Example**:

```cpp
// q = 0.05/0.01 = 5 is outside the table range, so this uses exact fallback:
double L1 = self_inductance(0.01, 0.03, 0.05, 1e5);

// q = 0.0002/0.01 = 0.02 < 0.05 — triggers exact fallback:
double L2 = self_inductance(0.01, 0.03, 0.0002, 1e5);

// Explicitly request reference-grade:
double L3 = self_inductance_exact(0.01, 0.03, 0.05, 1e5);

// Force exact via bool flag (same as self_inductance_exact):
double L4 = self_inductance(0.01, 0.03, 0.05, 1e5, true);
```

---

## Mutual Inductance

```cpp
#include <coilgun/physics/mutual_inductance.hpp>
```

### Filament Level — Closed Form

Two coaxial circular loops. This is the innermost kernel; calls are LRU-cached (4096 entries).

```cpp
double mutual_inductance_filament(
    double radius_a, double radius_b, double separation);     // (H), Eq.4.8

double mutual_inductance_gradient_filament(
    double radius_a, double radius_b, double separation);     // (H/m), Eq.4.15
```

`separation` is the axial distance between the two loop planes. The gradient `dM/dz` is **odd** in separation: swapping sign of separation flips the gradient sign.

**Example**: Two loops of radius 20mm and 25mm, 10mm apart:

```cpp
double M  = mutual_inductance_filament(0.02, 0.025, 0.01);   // ~ nH range
double dM = mutual_inductance_gradient_filament(0.02, 0.025, 0.01);
```

### Coil Level — 4D Gauss-Legendre Quadrature

Full finite-cross-section thick coils. Uses n-point Gauss-Legendre per dimension (n⁴ kernel evaluations), with LRU caching on every inner filament call.

```cpp
double mutual_inductance_coil(
    double r_inner_a, double r_outer_a, double length_a, int turns_a,
    double r_inner_b, double r_outer_b, double length_b, int turns_b,
    double separation,
    int n_nodes = 9);                                       // (H), Eq.4.13

double mutual_inductance_gradient_coil(
    double r_inner_a, double r_outer_a, double length_a, int turns_a,
    double r_inner_b, double r_outer_b, double length_b, int turns_b,
    double separation,
    int n_nodes = 9);                                       // (H/m), Eq.4.16
```

`separation` here is centre-to-centre axial distance (unlike the filament level, where it's plane-to-plane).

**Example** — mutual inductance between a driving coil and an armature segment modelled as a thick coil:

```cpp
double M_coil = mutual_inductance_coil(
    0.01, 0.03, 0.05, 150,    // coil A: ri, re, len, turns
    0.005, 0.025, 0.04, 1,    // coil B: ri, re, len, turns (armature segment)
    0.06);                     // centre-to-centre separation

double dM_coil = mutual_inductance_gradient_coil(
    0.01, 0.03, 0.05, 150,
    0.005, 0.025, 0.04, 1,
    0.06);
```

**Performance note**: coil-level functions are the most expensive calls in the library (~milliseconds each for n_nodes=9, ~tens of microseconds for n_nodes=4). In a time-marching simulation you will typically call these once per stage per step, not in tight inner loops.

### use_cache Overloads

Each filament-level and coil-level function has an overload accepting a trailing `bool use_cache = false` parameter:

```cpp
double mutual_inductance_filament(double radius_a, double radius_b,
                                  double separation, bool use_cache);
double mutual_inductance_gradient_filament(double radius_a, double radius_b,
                                           double separation, bool use_cache);
double mutual_inductance_coil(..., int n_nodes, bool use_cache);
double mutual_inductance_gradient_coil(..., int n_nodes, bool use_cache);
```

When `use_cache` is `true`, the function reads/writes a global 4096-entry LRU cache. This improves performance when the same `(radius_a, radius_b, separation)` triple is queried repeatedly — typical in serial cold-path computations like `[M]` matrix initialisation.

When `use_cache` is `false` (default), the cache is bypassed entirely. Use this on the hot path (per-step `[M_I]` updates) where cache hit rates are low and thread safety is required.

The parameterless overloads (without `use_cache`) are preserved for backward compatibility and behave identically to `use_cache = false`.

---

## LRU Cache

```cpp
#include <coilgun/physics/cache.hpp>

template<typename Key, typename Value,
         std::size_t MaxSize = 4096,
         typename KeyHash = std::hash<Key>>
class coilgun::physics::LRUCache {
public:
    bool        get(const Key& key, Value& value);   // true on hit
    void        put(const Key& key, const Value& value); // insert or update
    std::size_t size() const;                         // current entry count
    void        clear();                              // remove all entries
};
```

Standard `unordered_map + list` LRU. Thread-unsafe — each caller should own its own cache instance.

Used internally by `mutual_inductance_filament` (4096 entries) and `mutual_inductance_gradient_filament` (4096 entries). External users can instantiate their own:

```cpp
LRUCache<std::string, Eigen::Vector3d, 1024> position_cache;
position_cache.put("home", {0.0, 0.0, 0.0});
Eigen::Vector3d pos;
if (position_cache.get("home", pos)) { /* hit */ }
```

---

## DrivingCoil

```cpp
#include <coilgun/components/driving_coil.hpp>

class coilgun::components::DrivingCoil {
public:
    DrivingCoil(double inner_radius, double outer_radius, double length,
                int turns, double resistivity, double wire_area,
                double fill_factor, double position = -1.0,
                bool force_exact_self_inductance = false);

    // Geometry
    double inner_radius() const;     // ri (m)
    double outer_radius() const;     // re (m)
    double length() const;           // axial length (m)
    double mean_radius() const;      // (ri + re) / 2 (m)
    int    turns() const;            // N

    // Precomputed electrical properties
    double turns_density() const;    // nc = N / ((re - ri) × l)  (turns/m²)
    double resistance() const;       // DC resistance (Ω), accounts for fill factor
    double self_inductance() const;  // self-inductance (H) via T(q,p) table or exact

    // Position (for translation during simulation)
    double position() const;         // current centre position (m)
    void   set_position(double x);   // move coil centre to x
};
```

All electrical properties are precomputed at construction. The `position` field exists so the coil can be moved during simulation without recalculating geometry.

**Constructor parameters**:

| Parameter | Description | Typical value |
|-----------|-------------|---------------|
| `inner_radius` | Inner winding radius | 0.01 m |
| `outer_radius` | Outer winding radius | 0.03 m |
| `length` | Axial winding length | 0.05 m |
| `turns` | Total number of turns | 100–500 |
| `resistivity` | Wire resistivity at reference temp | `COPPER.resistivity_ref` |
| `wire_area` | Conductor cross-sectional area | 1e-6 m² (1 mm²) |
| `fill_factor` | Winding fill factor (0–1) | 0.5–0.8 |
| `position` | Initial centre position (m), default = length/2 | 0.0 m |
| `force_exact_self_inductance` | Force exact L calculation (skip T(q,p) table) | false |

The resistance formula accounts for the fill factor: only `fill_factor × cross_section` is conducting.

---

## Armature

```cpp
#include <coilgun/components/armature.hpp>

class coilgun::components::Armature {
public:
    Armature(double inner_radius, double outer_radius, double length,
             double resistivity, double material_density,
             double velocity, double mass,
             int m_axial, int n_radial, double position,
             physics::ArmatureMaterial material = physics::ArmatureMaterial::Aluminum,
             bool force_exact_self_inductance = false);

    // Geometry
    double inner_radius() const;
    double outer_radius() const;
    double length() const;
    int    axial_filaments() const;     // m (axial divisions)
    int    radial_filaments() const;    // n (radial divisions)
    int    total_filaments() const;     // m × n

    // Filament queries — 1-indexed to match NumericalModel notation
    double filament_inner_radius(int j) const;         // j = 1..n
    double filament_outer_radius(int j) const;
    double filament_mean_radius(int j) const;          // centre-line radius
    double filament_axial_position(int i) const;       // i = 1..m

    // Per-filament arrays — length m×n, row-major: (i=1,j=1), (i=1,j=2), ...
    const std::vector<double>& resistances() const;    // Ω per filament
    const std::vector<double>& inductances() const;    // H per filament
    const std::vector<double>& masses() const;         // kg per filament

    // Motion state
    double position() const;            // current centre position (m)
    double velocity() const;            // current velocity (m/s)
    double mass() const;                // total mass including payload (kg)
    void   update_position(double dx);  // translate by dx metres
    void   set_velocity(double v);      // set new velocity (m/s)

    // Material
    physics::ArmatureMaterial material() const;  // for thermal dispatch
};
```

The armature is a thick hollow cylinder discretised into m axial slices × n radial layers of current filaments. Each filament is treated as a ring carrying uniform current; its electrical properties (R, L) and mass are precomputed at construction.

**Constructor parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `inner_radius` | Inner bore radius (m) | — |
| `outer_radius` | Outer radius (m) | — |
| `length` | Axial length (m) | — |
| `resistivity` | Material resistivity (Ω·m), e.g. `ALUMINUM.resistivity_ref` | — |
| `material_density` | Material density (kg/m³) | — |
| `velocity` | Initial velocity (m/s) | — |
| `mass` | Total mass including payload (kg) — the filament masses sum to roughly the armature material mass; the difference is the payload | — |
| `m_axial` | Number of axial divisions (typ. 5–20) | — |
| `n_radial` | Number of radial divisions (typ. 1–5) | — |
| `position` | Initial centre position (m) | — |
| `material` | Armature material: `Aluminum` or `Copper` (determines cp(T) and beta for thermal mode) | `Aluminum` |
| `force_exact_self_inductance` | Force exact filament self-inductance | `false` |

**Filament indexing**: uses 1-based indexing to match the convention in NumericalModel.md. `filament_axial_position(1)` gives the position of the leftmost axial ring, `filament_axial_position(m)` gives the rightmost.

**Data layout**: the three per-filament arrays (`resistances()`, `inductances()`, `masses()`) are flat vectors of length `m × n` in **row-major** order. The i-th axial slice, j-th radial layer is at index `(i - 1) * n + (j - 1)`.

**Example** — a 100g aluminium armature with a 20g payload, 5 axial × 2 radial filaments:

```cpp
Armature arm(
    0.005,                     // inner radius: 5 mm bore
    0.025,                     // outer radius: 25 mm
    0.08,                      // length: 80 mm
    ALUMINUM.resistivity_ref,  // 2.82e-8 Ω·m
    ALUMINUM.density,           // 2700 kg/m³
    0.0,                       // initial velocity
    0.120,                     // total mass = 120 g (100 g armature + 20 g payload)
    5, 2,                      // m=5 axial, n=2 radial = 10 filaments total
    0.0                        // initial position
);

std::cout << "Filaments: " << arm.total_filaments() << "\n";           // 10
std::cout << "Ring (i=3,j=1) radius: " << arm.filament_mean_radius(1) << " m\n";
std::cout << "Ring (i=3,j=1) pos:    " << arm.filament_axial_position(3) << " m\n";

// Access per-filament data by flat index:
int idx = (3 - 1) * 2 + (1 - 1);   // (i=3, j=1) → index 4
std::cout << "Resistance: " << arm.resistances()[idx] << " Ω\n";
std::cout << "Mass:       " << arm.masses()[idx] << " kg\n";
```

---

## Simulation Engine

The `coilgun::simulation` namespace provides a turn-key single-stage coilgun simulator with configurable excitation, time-integration policy, and thermal coupling.

### Excitation Framework

```cpp
#include <coilgun/simulation/excitation.hpp>
```

| Class | Description | Key constructor args |
|-------|-------------|---------------------|
| `Excitation` | Abstract base — `voltage()`, `advance(dt, I_coil)`, `finished()`, `reset()` | — |
| `CapacitorExcitation` | Capacitor discharge (crowbar disabled) | `(initial_voltage, capacitance)` |
| `CrowbarExcitation` | Capacitor discharge **with** crowbar diode | `(initial_voltage, capacitance)` |
| `WaveformExcitation` | Arbitrary voltage source `V(t)` | `(V_of_t function)` |

All excitations provide `voltage()`, `advance(dt, I_coil)`, `finished()`, and `reset()`. `CapacitorExcitation` additionally exposes `capacitance()`, `capacitor_voltage()`, and `initial_voltage()`. `CrowbarExcitation` reports crowbar state via `diode_on()`. `WaveformExcitation` supports optional early termination via `set_end_time(t)`.

```cpp
auto cap  = std::make_unique<CapacitorExcitation>(450.0, 0.001);  // 450 V, 1000 μF
auto crow = std::make_unique<CrowbarExcitation>(450.0, 0.001);    // crowbar variant
auto wfm  = std::make_unique<WaveformExcitation>(
    [](double t) { return 450.0 * std::exp(-t / 0.01); });
wfm->set_end_time(0.05);  // stop waveform at 50 ms

// Query excitation state:
double v_cap  = crow->capacitor_voltage();  // current capacitor voltage
double v_init = crow->initial_voltage();    // initial capacitor voltage (450 V)
bool diode_on = crow->diode_on();           // crowbar diode conducting?
crow->reset();                               // restore initial state
```

### Time Steppers

```cpp
#include <coilgun/simulation/time_stepper.hpp>

// EulerStepper     — forward Euler (1st order)
// RK4Stepper       — classical 4th-order Runge-Kutta
```

These are template policies. Pass them as the template argument to `SingleStageSim` or `MultiStageSim`:

```cpp
SingleStageSim<EulerStepper> sim_euler(...);   // faster, less accurate
SingleStageSim<RK4Stepper>   sim_rk4(...);     // slower, more accurate

MultiStageSim<EulerStepper> ms_euler(...);     // ditto for multi-stage
MultiStageSim<RK4Stepper>   ms_rk4(...);
```

### SimState

```cpp
#include <coilgun/simulation/single_stage_sim.hpp>

struct coilgun::simulation::SimState {
    Eigen::VectorXd currents;              // [N_fil+1], index 0 = coil current, [1..N_fil] = filaments
    double          arm_position = 0.0;    // m
    double          arm_velocity = 0.0;    // m/s
    Eigen::VectorXd filament_temperatures; // [N_fil], K (only populated when thermal enabled)

    SimState& operator+=(const SimState& rhs);
    SimState& operator*=(double scalar);
};

SimState operator+(SimState lhs, const SimState& rhs);
SimState operator*(double scalar, SimState s);
```

`SingleStageSim::state()` returns a reference to the current internal `SimState`. The arithmetic operators support the internal stepper algorithms.

### MultiStageState

```cpp
#include <coilgun/simulation/multi_stage_sim.hpp>

struct coilgun::simulation::MultiStageState {
    Eigen::VectorXd currents;              // [n_stages + N_fil], stages first, then filaments
    double          arm_position = 0.0;    // m
    double          arm_velocity = 0.0;    // m/s
    Eigen::VectorXd filament_temperatures; // [N_fil], K (only populated when thermal enabled)

    MultiStageState& operator+=(const MultiStageState& rhs);
    MultiStageState& operator*=(double scalar);
};

MultiStageState operator+(MultiStageState lhs, const MultiStageState& rhs);
MultiStageState operator*(double scalar, MultiStageState s);
```

`MultiStageSim::state()` returns a reference to the current internal `MultiStageState`.

### Simulation Results

```cpp
#include <coilgun/simulation/sim_result.hpp>

struct SimStep {
    double time;
    double cap_voltage;
    double coil_current;
    std::vector<double> filament_currents;
    double arm_position;
    double arm_velocity;
    double force;               // signed instantaneous net axial force, N
    std::vector<double> filament_temperatures;  // only when thermal enabled
};

struct SimSummary {
    double muzzle_velocity;    // m/s
    double total_time;         // s
    double max_force;          // peak magnitude of signed instantaneous force, N
    double peak_coil_current;  // A
    double efficiency;         // 0–1
    int    step_count;
};

struct SimResult {
    std::vector<SimStep> history;
    SimSummary           summary;
    SimResult            sampled(int every_n) const;  // down-sample for plotting
};

// Use sampled() to reduce data for export:
auto sparse = result.sampled(100);  // keep every 100th step
```

### Termination

```cpp
#include <coilgun/simulation/termination.hpp>

struct TerminationPolicy {
    int    max_steps             = 20000;
    int    velocity_decay_steps  = 5;      // consecutive steps of velocity drop
    double accel_threshold       = 0.1;    // m/s² — treat near-zero acceleration as stalled
    bool   enable_velocity_check = true;
    bool   enable_bound_check    = false;
    double barrel_end_position   = 1.0;

    static TerminationPolicy defaults();
};

// Custom policy:
TerminationPolicy pol = TerminationPolicy::defaults();
pol.max_steps = 50000;
sim.run(pol);
```

### SingleStageSim

```cpp
#include <coilgun/simulation/single_stage_sim.hpp>

template<typename StepperPolicy = EulerStepper>
class SingleStageSim {
public:
    SingleStageSim(
        components::DrivingCoil          coil,          // copied
        components::Armature             armature,      // copied
        std::unique_ptr<Excitation>      excitation,    // moved-in
        double                           dt,
        bool                             enable_thermal = false
    );

    const SimStep&   step();                       // advance one step
    const SimResult& run();                        // run to default termination
    const SimResult& run(const TerminationPolicy& p); // run to custom termination
    void             reset();                      // reset to initial state

    const SimResult& result()     const;   // result after last run()
    const SimState&  state()      const;   // internal state (currents, pos, vel, temps)
    double           dt()         const;   // time step
    int              step_count() const;   // steps since construction or last reset()
};

// Re-run with different parameters:
sim.run();
double v1 = sim.result().summary.muzzle_velocity;
sim.reset();          // restore initial state (including excitation)
sim.run();
double v2 = sim.result().summary.muzzle_velocity;  // identical to v1
```

### OptimizationLevel

```cpp
#include <coilgun/simulation/multi_stage_sim.hpp>

enum class OptimizationLevel {
    Reference   = 0,   // Fixed 9-point mutual-inductance quadrature, no distance cutoff.
    LookupTable = 1,   // Reserved for component-level T(q,p) lookup selection; same runtime path as Reference.
    Full        = 2    // Distance cutoff plus adaptive 9/4-point mutual-inductance quadrature.
};
```

Levels are cumulative: each level includes all optimisations of lower levels.

| Level | T(q,p) table | Distance cutoff | Adaptive GL order |
|:-----:|:------------:|:---------------:|:-----------------:|
| 0 Ref | No | No | No (fixed 9) |
| 1 Table | Yes | No | No (fixed 9) |
| 2 Full | Yes | Yes | Yes (near=9, far=4) |

### Trigger Configuration

```cpp
#include <coilgun/simulation/trigger_config.hpp>

enum class TriggerMode {
    Position,    // Trigger when armature centre passes a given position (m).
    TimeDelay    // Trigger after a fixed delay from the previous stage (s).
};

struct TriggerConfig {
    TriggerMode mode;
    double      value;       // position (m) or time delay (s).
};
```

Stage 0 always triggers at t=0 and does NOT need a TriggerConfig. `trigger_configs[i-1]` applies to stage `i` (i >= 1).

`validate_trigger_config()` rejects invalid modes and NaN for both modes. `Position` accepts every finite value but rejects negative infinity; `TimeDelay` accepts finite non-negative values and rejects negative values and negative infinity. Positive infinity is valid in either mode and is an explicit terminal never-trigger policy: it represents an unreachable position or delay rather than malformed input. A finite untriggered stage remains eligible even after an earlier stage finishes, so it keeps the simulation alive until it triggers or another termination criterion wins. Automatic completion occurs only when every stage is finished or terminally ineligible through a `+infinity` policy. Trigger positions are captured at the pre-step boundary where the stage fires.

### Multi-Stage Simulation Results

```cpp
#include <coilgun/simulation/multi_stage_result.hpp>

struct StepSnapshot {
    double time;
    double arm_position;        // m
    double arm_velocity;        // m/s
    double force;               // signed instantaneous net force, N
    std::vector<double> filament_currents;
    std::vector<double> filament_temperatures;  // K (thermal mode only)
};

struct MultiStageStep {
    StepSnapshot          state;
    std::vector<double>   cap_voltages;     // [n_stages], 0 for inactive stages
    std::vector<double>   coil_currents;    // [n_stages], 0 for inactive stages
    std::vector<double>   stage_forces;     // [n_stages], signed instantaneous post-step contributions, N
};

struct PerStageSummary {
    int    stage_index;
    double trigger_time;         // s
    double trigger_position;     // armature position at the pre-step trigger boundary, m
    double peak_current;         // A
    double max_force;            // Peak magnitude of this stage's signed force contribution, N
    double energy_depleted;      // J (capacitor energy consumed)
    int    step_count_active;
};

struct MultiStageSummary {
    double                        muzzle_velocity;      // m/s
    double                        total_time;           // s
    double                        max_force;            // Peak magnitude of signed net force, N
    double                        peak_coil_current;    // A
    double                        efficiency;           // 0..1, E_kin / Sigma(0.5*C_i*U0_i^2)
    int                           step_count;
    std::vector<PerStageSummary>  per_stage;
};

struct MultiStageResult {
    std::vector<MultiStageStep> history;
    MultiStageSummary           summary;
    MultiStageResult            sampled(int every_n) const;  // down-sample for plotting
};
```

### MultiStageSim

```cpp
#include <coilgun/simulation/multi_stage_sim.hpp>

template<typename StepperPolicy = EulerStepper>
class MultiStageSim {
public:
    static constexpr int kMaxStages = 50;

    MultiStageSim(
        std::vector<components::DrivingCoil>        coils,          // copied
        components::Armature                         armature,      // copied
        std::vector<std::unique_ptr<Excitation>>     excitations,   // moved-in
        std::vector<TriggerConfig>                   trigger_configs,
        double                                       dt,
        bool                                         enable_thermal = false,
        OptimizationLevel                            opt_level = OptimizationLevel::Reference
    );

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy& policy);
    void                    reset();

    const MultiStageResult& result()      const;
    const MultiStageState&  state()       const;
    double  dt()              const;     // time step (s)
    int     step_count()      const;     // step count
    int     num_stages()      const;     // configured stage count
};
```

**Constructor constraints:** `coils.size() == excitations.size()`, `trigger_configs.size() == coils.size() - 1`, `coils.size() <= kMaxStages (50)`.

**System details:** The ODE dimension is `n_stages + N_filaments`. Inactive stages (not yet triggered) have their rows/columns set to identity to keep the matrix non-singular. Each stage's CrowbarExcitation autonomously manages its diode state. Termination occurs when all stages have finished (crowbar diode ON and current decayed).

**Thermal mode:** When `enable_thermal = true`, each armature filament undergoes adiabatic ohmic heating per NumericalModel Sec.6. Resistance updates every step based on temperature-dependent resistivity.

**Optimization:** The `opt_level` parameter controls runtime mutual-inductance work independently of the components' self-inductance calculation mode, which is fixed at component construction. `Reference` and `LookupTable` currently use the same runtime path: fixed 9-point Gauss-Legendre quadrature with no distance cutoff. `LookupTable` is retained for API compatibility and does not retroactively change component self-inductances. `Full` skips distant stages and uses 9-point near / 4-point far quadrature.

---

## Putting It Together

### Single-Stage Example

A complete single-stage simulation with crowbar diode, Euler integration, no thermal coupling:

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

    auto excitation = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(excitation), 1e-6, false);

    sim.run();

    const auto& r = sim.result();
    std::cout << "Muzzle velocity: " << r.summary.muzzle_velocity << " m/s\n";
    std::cout << "Efficiency:      " << r.summary.efficiency * 100 << " %\n";
    std::cout << "Steps:           " << r.summary.step_count << "\n";
    return 0;
}
```

The `SingleStageSim` handles the full circuit ODE (coil + armature filaments + mutual coupling), Lorentz force integration, kinematic update, and automatic termination — all through the chosen stepper policy.

### Multi-Stage Example

A complete two-stage simulation with crowbar diodes and position-based trigger:

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>
#include <memory>
#include <vector>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::MultiStageSim;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;

    // Build two coils for a 2-stage gun
    DrivingCoil coil1(0.01, 0.03, 0.05, 150,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil coil2(0.01, 0.03, 0.05, 150,
                      COPPER.resistivity_ref, 1e-6, 0.7, 0.10);

    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    std::vector<DrivingCoil> coils = {coil1, coil2};

    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));

    // Stage 1 triggers when armature centre passes 0.09 m
    std::vector<TriggerConfig> triggers;
    triggers.push_back({TriggerMode::Position, 0.09});

    MultiStageSim<EulerStepper> sim(
        std::move(coils), arm, std::move(excs), triggers, 1e-6);

    sim.run();

    const auto& r = sim.result();
    std::cout << "Muzzle velocity: " << r.summary.muzzle_velocity << " m/s\n";
    std::cout << "Efficiency:      " << r.summary.efficiency * 100 << " %\n";
    std::cout << "Stages:          " << r.summary.per_stage.size() << "\n";
    for (const auto& ps : r.summary.per_stage) {
        std::cout << "  Stage " << ps.stage_index
                  << ": trigger_t=" << ps.trigger_time
                  << "s, I_peak=" << ps.peak_current << " A\n";
    }
    return 0;
}
```

The `MultiStageSim` handles the full multi-stage circuit ODE (inter-coil mutual coupling, motional back-EMF, autonomous crowbar diodes per stage), Lorentz force integration, kinematic update, automatic stage triggering, and termination — all through the chosen stepper policy.

---

## Parallel Execution

The simulation engine uses OpenMP for shared-memory multicore parallelism. Three computational loops are parallelised per time step:

- Coil-to-filament mutual inductance computation (`M` and `dM/dx`)
- Lorentz force summation (reduction)
- Filament temperature update (when thermal mode is enabled)

### Controlling Thread Count

```cpp
#include <omp.h>
omp_set_num_threads(4);
```

Or via environment variable:

```sh
OMP_NUM_THREADS=4 ./your_simulation
```

Default thread count equals the number of logical CPUs.

### Thread Safety

The `mutual_inductance_filament` and `mutual_inductance_gradient_filament` functions are thread-safe when called with `use_cache = false` (the default). The global LRU cache is only accessed from serial cold-path code (`[M]` matrix initialisation, inter-coil mutual inductance precomputation).

### CMake Integration

```cmake
find_package(OpenMP REQUIRED)
target_link_libraries(your_target PRIVATE coilgun OpenMP::OpenMP_CXX)
```

If linking against the static library `libcoilgun.a`, you must also link OpenMP in your target.

### Compiler Flags

All CMake presets include `-march=native` for SIMD instruction set auto-detection. Eigen leverages this at compile time to enable AVX2/FMA/AVX512 for matrix operations.

---

## GPU Acceleration (CUDA)

The library provides an optional GPU-accelerated backend (`libcoilgun_cuda.a`) that offloads the 4D Gauss-Legendre mutual inductance integration to the GPU. The CPU path is completely unchanged — the GPU classes are separate, drop-in replacements.

```sh
cmake --preset ninja-cuda-debug      # configure with CUDA
cmake --build --preset ninja-cuda-debug  # build
ctest --preset debug                 # run all tests (CPU + GPU)
```

### Prerequisites

| Requirement | Minimum | Notes |
|---|---|---|
| CUDA Toolkit | ≥ 12.8 (Blackwell sm_120); ≥ 9.0 (earlier arch) | `nvcc` must be in PATH |
| Boost.Math | ≥ 1.86 | Fetched automatically via CMake |
| NVIDIA GPU | Compute Capability ≥ 6.0 (Pascal+) | Required for FP64 atomic support |
| Host compiler | g++ ≥ 9 or equivalent | Must match the one nvcc detects |

### Include and Link

```cpp
#include <coilgun/coilgun_cuda.hpp>   // all GPU classes + CPU umbrella
```

```cmake
target_link_libraries(your_target PRIVATE coilgun_cuda coilgun CUDA::cudart)
```

The `coilgun_cuda` target depends on `coilgun` transitively — both libraries are needed.

---

### GpuBackend

```cpp
#include <coilgun/simulation/cuda/gpu_backend.hpp>

namespace coilgun::simulation::cuda {

struct GpuBackend {
    int     device_id         = 0;     ///< cudaSetDevice target.
    int     threads_per_block = 512;   ///< Threads per block for the 4D integration kernel.
    size_t  max_batch_sims    = 256;   ///< Pre-allocation cap for batch simulation buffers.
    bool    enable_profiling  = false; ///< Retain profiling-request metadata; host-wall timing fields are always collected. No NVTX guarantee.
    bool    use_persistent    = true;  ///< Legacy persistent request for multi-stage and batch wrappers.
    BackendMode backend        = BackendMode::Graph; ///< Backend when persistence is not requested.
};

}
```

| Field | Purpose | Default |
|-------|---------|---------|
| `device_id` | Selects which physical GPU to target (relevant on multi-GPU systems). | `0` |
| `threads_per_block` | Requested number of threads per block for the 4D GL integration kernel. Must be a positive power of two no greater than 512; 1, 128, 256 and 512 are valid. | `512` |
| `max_batch_sims` | Maximum number of simulations in a `SimBatch`. `SimBatch` rejects a larger `num_sims` with `std::invalid_argument`; the field itself does not allocate buffers. | `256` |
| `enable_profiling` | When true, retains the request in `ExecutionReport::profiling_enabled`. Host-wall timing categories are collected independently of this flag. This build does not promise NVTX annotations or require an NVTX dependency. | `false` |
| `use_persistent` | Legacy compatibility request. The omitted `GpuMultiStageSim` backend uses `multi_stage_default_backend()` and selects `Direct`; a caller that supplies `use_persistent=true` still intentionally requests `Persistent`. In `SimBatch`, this flag is consulted only when `backend.backend == Auto`; explicit backend modes also keep identical deterministic and capability policy for either flag value. | `true` |
| `backend` | Backend metadata/configuration. In `GpuSingleStageSim`, `Graph`, `Direct`, or `Fallback` may be selected with `use_persistent=false`. In `GpuMultiStageSim`, use the constructor's explicit `BackendMode` parameter to distinguish an intentional `Graph`, `Direct`, or `Fallback` request from the legacy flag. Graph captures the mutual segment only; CUDA unavailability or runtime failure is reported as `Fallback`. | `Graph` |

Fault-injection controls are not part of `GpuExecutionConfig`. Focused tests use the separate internal `detail::GpuEngineFaultInjection` constructor seam; ordinary execution configuration contains only runtime policy and validated launch settings.

**Migration status**: `GpuEngine` is the current execution core. `GpuSingleStageSim`, `GpuMultiStageSim`, and `SimBatch` use the engine contract for backend selection, rollback, reporting, and partial Graph capture. A Graph report means graph-assisted mutual-segment execution only; it does not claim a fully captured multi-stage pipeline.

**Example**:

```cpp
coilgun::simulation::cuda::GpuBackend be;
be.device_id = 0;
be.threads_per_block = 256;   // any positive power of two <= 512
```

---

### GpuOptLevel

```cpp
#include <coilgun/simulation/cuda/gpu_backend.hpp>

namespace coilgun::simulation::cuda {

enum class GpuOptLevel {
    Standard   = 0,   ///< n_nodes=9, no distance cutoff. Use for debugging/verification.
    Full       = 1,   ///< Distance cutoff (>10× coil length) + fixed n_nodes=9. Production.
    Aggressive = 2,   ///< FP32 integrand + FP64 reduction, distance cutoff, n_nodes=9. Large-scale sweeps.
};

}
```

| Level | Distance cutoff | GL order | Behaviour |
|:-----:|:---------------:|:--------:|-----------|
| `Standard` | No | n_nodes = 9 | All (stage, filament) pairs are processed. The kernel uses 9 GL nodes per dimension (6561 integration points). This is the safest setting. |
| `Full` | Yes | n_nodes = 9 | Stages farther than 10× the coil length from the armature are skipped. At typical distances all active stages are within range. |
| `Aggressive` | Yes | n_nodes = 9 | Same as `Full` but uses FP32 arithmetic for the integrand (including FP32 AGM elliptic integrals) with FP64 accumulation. 3-5× faster on consumer GPUs. Slight precision tradeoff. |

**Note**: The GPU backend does NOT support adaptive GL order (n_nodes=4/9 mix). Using n_nodes=4 on the GPU causes non-deterministic floating-point drift in the shared-memory reduction (only 3 active warps out of 8). The CPU `OptimizationLevel::Full` uses both distance cutoff and adaptive quadrature; the GPU `GpuOptLevel::Full` uses distance cutoff only.

---

### Unified Execution Configuration and Planning

These host-only types resolve the requested execution contract before any CUDA resource is created:

```cpp
enum class BackendMode {
    Auto = 0, Graph = 1, Persistent = 2, Fallback = 3, Direct = 4
};
enum class SolverMode {
    Auto = 0, Eigen = 1, Batched = 2, CuSolver = Batched
};
enum class PrecisionMode { Standard = 0, Full = 1, Aggressive = 2 };
enum class ThermalMode { Auto = 0, Disabled = 1, Cpu = 2, Gpu = 3 };
enum class FallbackReason {
    None = 0, CapabilityUnavailable = 1, DeterminismRequired = 2,
    RuntimeFailure = 3, MetadataConflict = 4
};

struct GpuExecutionConfig {
    BackendMode backend = BackendMode::Auto;
    SolverMode solver = SolverMode::Auto;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Auto;
    bool enable_calibration = false;
    bool deterministic = false;
    int device_id = 0;
    int threads_per_block = 512;
    bool enable_profiling = false;
    void validate() const;
};

struct GpuCapability {
    bool supports_graph = true;
    bool supports_persistent = true;
    bool supports_batched_solver = true;
    bool supports_gpu_thermal = true;
    bool persistent_is_deterministic = false;
};

struct GpuExecutionPolicy {
    BackendMode requested_backend = BackendMode::Auto;
    SolverMode requested_solver = SolverMode::Auto;
    PrecisionMode requested_precision = PrecisionMode::Full;
    ThermalMode requested_thermal = ThermalMode::Auto;
    BackendMode backend = BackendMode::Fallback;
    SolverMode solver = SolverMode::Eigen;
    PrecisionMode precision = PrecisionMode::Full;
    ThermalMode thermal = ThermalMode::Disabled;
    FallbackReason backend_fallback_reason = FallbackReason::None;
    FallbackReason solver_fallback_reason = FallbackReason::None;
    FallbackReason thermal_fallback_reason = FallbackReason::None;
};

class GpuExecutionPlanner {
public:
    static constexpr GpuExecutionPolicy plan(
        std::size_t n_stages, std::size_t n_filaments, std::size_t batch_size,
        bool thermal_enabled, const GpuCapability&, const GpuExecutionConfig&) noexcept;
};
```

Constraints and static rules:

- `GpuExecutionPlanner::plan` is pure host code. It uses dimensions, explicit requests, deterministic mode, and the supplied capability snapshot; it does not inspect CUDA or timings.
- An explicit `Fallback` is CPU-only and must not create a CUDA context. The host planner preserves an independently requested solver/thermal choice for policy inspection, while `GpuEngine` normalizes the runtime selection to `Fallback + Eigen` and CPU thermal as needed.
- An explicit `Graph` with `supports_graph == false` resolves to `Fallback` with `CapabilityUnavailable`; the runtime then follows the resolved backend and does not create a context. With graph capability, the current engine may upgrade the request at construction and captures/replays only the mutual segment. Matrix assembly, solve, force/state orchestration, and thermal updates remain outside the graph.
- An explicit `Persistent` is currently resolved to a safe fallback because the synchronous engine does not own the dedicated control stream required by a resident kernel. A deterministic request also rejects a capability marked nondeterministic.
- `Direct` is the direct CUDA mutual pipeline when CUDA is compiled, the selected device is usable, and initialization succeeds. The host-only planner leaves `Auto` conservative at `Fallback`; `GpuEngine` may upgrade `Auto` to `Direct` after runtime device detection. Callers requiring a particular CUDA backend should request `Direct` or `Graph` explicitly.
- `SolverMode::Batched` requires `supports_batched_solver`; otherwise the solver resolves to Eigen with `CapabilityUnavailable`. `SolverMode::Auto` selects Batched only for the planner's large-workload rule, otherwise Eigen.
- `ThermalMode::Gpu` requires `supports_gpu_thermal`; otherwise it resolves to CPU thermal with `CapabilityUnavailable`. `Auto` selects GPU thermal only for a large workload with capability support.
- The current large-workload rule is `batch_size >= 8 || n_stages + n_filaments >= 128`. It is a static rule, not a performance guarantee.
- If `enable_calibration` is true, construction performs one identity-batch solver calibration after workspace initialization. Calibration does not advance physical state and is not repeated by `reset()`; `ExecutionReport::calibrated` is metadata for that one-time operation.
- A capability snapshot is not proof that a device is available. CUDA runtime enumeration, device selection, context creation, and allocation can still fail and are reported as runtime fallback. Callers must inspect the resolved report.

`GpuCapability` and the planner are public for deterministic host testing and policy inspection. They are not a promise that CUDA runtime detection has already occurred.

Enum semantics:

| Type | Values |
|---|---|
| `BackendMode` | `Auto` defers the runtime choice; `Graph` requests graph-assisted mutual-segment capture/replay; `Persistent` requests the resident protocol but currently falls back in `GpuEngine`; `Fallback` is CPU-only; `Direct` launches the CUDA mutual segment directly. |
| `SolverMode` | `Auto` applies the static workload rule; `Eigen` is the host LDLT path; `Batched` requests the CUDA batched solver. `CuSolver` is an alias of `Batched`. |
| `PrecisionMode` | `Standard` is FP64 without the distance cutoff, `Full` is FP64 with the production cutoff, and `Aggressive` uses an FP32 integrand with FP64 reduction. |
| `ThermalMode` | `Auto` applies the workload/capability rule; `Disabled` omits thermal updates; `Cpu` and `Gpu` select the corresponding thermal path when supported. |
| `FallbackReason` | `None` means no fallback; `CapabilityUnavailable` is an unsupported static/runtime capability; `DeterminismRequired` denotes a rejected nondeterministic choice; `RuntimeFailure` is context/allocation/capture/step failure; `MetadataConflict` records a conservative or currently unimplemented request resolution. |

Configuration and policy fields:

| Type | Field group | Meaning |
|---|---|---|
| `GpuExecutionConfig` | `backend`, `solver`, `precision`, `thermal` | Requested execution modes. |
| `GpuExecutionConfig` | `enable_calibration`, `deterministic` | Request one-time construction calibration and deterministic policy constraints. |
| `GpuExecutionConfig` | `device_id`, `threads_per_block` | Runtime device selection and validated launch width. |
| `GpuExecutionConfig` | `enable_profiling` | Metadata copied to the report; timings are collected independently and NVTX is not implied. |
| `detail::GpuEngineFaultInjection` | `fail_after_mutual`, `fail_allocation`, `fail_device_initialization`, `fail_graph_capture` | Internal/test-only constructor seam; not ordinary execution configuration. |
| `GpuCapability` | `supports_graph`, `supports_persistent`, `supports_batched_solver`, `supports_gpu_thermal` | Caller-supplied capability snapshot used only by static planning. |
| `GpuCapability` | `persistent_is_deterministic` | Whether Persistent may satisfy `deterministic=true`. |
| `GpuExecutionPolicy` | `requested_*` | Original requested modes retained for audit. |
| `GpuExecutionPolicy` | `backend`, `solver`, `precision`, `thermal` | Resolved modes that drive resource creation and execution. |
| `GpuExecutionPolicy` | `*_fallback_reason` | Per-dimension static resolution reasons. |

---

### GpuSingleStageSim

```cpp
#include <coilgun/simulation/cuda/gpu_single_stage_sim.hpp>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class GpuSingleStageSim {
public:
    GpuSingleStageSim(
        components::DrivingCoil          coil,          // copied
        components::Armature             armature,      // copied
        std::unique_ptr<Excitation>      excitation,    // moved-in
        double                           dt,
         bool                             enable_thermal = false,
         GpuOptLevel                      opt_level = GpuOptLevel::Full,
         const GpuBackend&                backend = {},
         BackendMode                      explicit_backend = BackendMode::Auto
    );

    const SimStep&   step();
    const SimResult& run();
    const SimResult& run(const TerminationPolicy& p);
    void             reset();

    const SimResult& result()     const;
    const SimState&  state()      const;
    double  dt()         const;
    int     step_count() const;
    const ExecutionReport& execution_report() const;
    std::vector<double> filament_resistances() const;
};

} // namespace
```

**Constructor parameters**:

| Parameter | Description | Typical value |
|-----------|-------------|---------------|
| `coil` | Driving coil geometry (copied) | — |
| `armature` | Armature geometry and filament discretisation (copied) | — |
| `excitation` | Excitation source (moved-in) | `CrowbarExcitation(450.0, 0.001)` |
| `dt` | Fixed time step (s) | `1e-6` |
| `enable_thermal` | Enable adiabatic filament heating (CPU-side) | `false` |
| `opt_level` | GPU optimisation level | `GpuOptLevel::Full` |
| `backend` | GPU backend configuration | `{}` (defaults) |
| `explicit_backend` | Optional constructor-level backend override. `Auto` preserves `backend`; another value overrides `backend.backend` and `use_persistent`. | `BackendMode::Auto` |

The constructor throws `std::invalid_argument` for a null excitation, non-finite or non-positive `dt`, negative `device_id`, invalid `threads_per_block`, `max_batch_sims == 0`, non-finite excitation voltage, invalid geometry/state dimensions, or non-finite stage voltages. CUDA-enabled construction also validates that the configured device exists and remains selected before allocation. `GpuSingleStageSim<RK4Stepper>` remains constructible for source compatibility, but `step()` throws `std::logic_error` because RK4 is unsupported.

**Methods**:

| Method | Behaviour |
|--------|-----------|
| `step()` | Advance one Euler time step. `Direct` launches the existing CUDA mutual segment; `Graph` captures/replays that mutual segment and then runs the existing synchronous matrix, Eigen LDLT, force/state and thermal stages outside the graph. `Fallback` runs the complete step on CPU. In `Full`/`Aggressive`, stage-armature mutual coupling and force are gated off when `abs(armature_position - coil.position()) > 10 * coil.length()`; the driving circuit remains active. Returns the newly recorded `SimStep`. |
| `run()` | Run to default termination. Blocks until the excitation finishes (crowbar diode ON and current decayed). |
| `run(policy)` | Run to custom `TerminationPolicy`. |
| `reset()` | Restore to initial state. Resets all currents to zero, restores armature position/velocity to construction values, resets the excitation, clears result history. |
| `result()` | Result after the last `run()`. Contains `history` (vector of SimStep) and `summary` (SimSummary). |
| `state()` | Current internal state — live reference. |
| `dt()` | Fixed time step. |
| `step_count()` | Steps since construction or last `reset()`. |
| `execution_report()` | Resolved backend/solver/thermal policy and cumulative diagnostics. `gpu_executed` is false on fallback and becomes true only after a complete CUDA physical step commits successfully. GPU/solver/thermal/transfer timings, graph rebuilds, fallback count, calibration, condition estimate and fallback reasons are cumulative across `reset()`. |
| `filament_resistances()` | Returns a value copy of the current filament resistance vector. Thermal construction initializes it from the armature reference resistances; Joule heating updates it, and `reset()` restores those references. The copy remains valid across wrapper mutation and exposes no device pointer. |

`SimStep::force` is the signed force recorded from the committed post-step currents and the mutual-gradient cache computed at the pre-position step boundary. The position update does not trigger a second gradient evaluation; this is the established post-current/pre-position convention shared with the multi-stage and batch wrappers. The force used for the Euler velocity update is the force computed from pre-step currents.

`ExecutionReport` field semantics:

| Fields | Meaning |
|---|---|
| `requested_backend`, `requested_solver`, `requested_precision`, `requested_thermal` | Original requests retained for audit. |
| `backend`, `solver`, `precision`, `thermal` | Current resolved execution modes. |
| `gpu_executed` | Cumulative proof that at least one complete CUDA-backed physical step committed successfully; not a capability/request flag. |
| `calibrated`, `precision_fallback`, `metadata_conflict` | One-time calibration completion and accumulated policy/report diagnostics. |
| `graph_rebuild_count`, `fallback_count` | `graph_rebuild_count` counts only successful captures of a new CUDA Graph variant; host variant selection, Direct/Fallback construction, cache hits, and failed captures do not increment it. `fallback_count` counts fallback events. |
| `gpu_time_ms` | Cumulative host wall time for successful CUDA-backed physical pipelines, including transfers, CPU matrix assembly/Eigen solving, and host orchestration; never device-only kernel time. |
| `solver_time_ms`, `thermal_time_ms` | Cumulative host wall time for those sections on whichever CPU or CUDA-backed path ran. |
| `transfer_time_ms` | Cumulative host wall time spent in synchronous host/device copies. |
| `max_condition_estimate` | Maximum recorded solver condition estimate/calibration diagnostic. |
| `fallback_reason` | Human-readable latest retained fallback message. |
| `static_fallback_reason`, `runtime_fallback_reason` | Machine-readable planning and runtime reasons. |
| `device_id`, `threads_per_block` | Validated execution settings copied from configuration. |
| `profiling_enabled` | Profiling-request metadata only. Timings are collected independently; NVTX is not implied. |

**Ownership and fallback**: The wrapper owns a `GpuEngine` through RAII. The `Excitation` remains the moved-in `std::unique_ptr`; the wrapper owns no raw device pointers. `GpuEngine` owns CUDA context, solver, graph, thermal, and device-buffer lifetimes and frees partially allocated resources during initialization failure. If CUDA runtime enumeration or device selection fails, CUDA is unavailable, context creation/allocation fails, or a runtime pipeline failure occurs, the engine restores the complete pre-step state when needed, initializes an Eigen solver, runs the whole step on CPU, and locks subsequent steps to CPU fallback. Enumeration/selection/context/allocation/pipeline errors retain a nonempty reason and `runtime_fallback_reason=RuntimeFailure`; no-device/insufficient-driver availability uses `CapabilityUnavailable`. A successful CUDA step and its CPU fallback are intended to preserve the same Euler physical pipeline; small CPU/GPU floating-point differences remain possible in mutual-inductance results and can affect long runs.

**Integration contract**: `GpuSingleStageSim<EulerStepper>` is supported. `GpuSingleStageSim<RK4Stepper>` is intentionally unsupported in this migration: `step()` throws `std::logic_error` rather than silently executing Euler. True four-stage RK4 parity is deferred and is not claimed by this API.

**Reset and diagnostics**: `reset()` clears simulation state, excitation state, completed-step count, and result history, then reselects the graph variant at the step boundary. `ExecutionReport` is an execution audit: fallback count, timings, graph rebuild count, calibration status, maximum condition estimate, fallback reasons, and `gpu_executed` are cumulative and are retained across reset. If calibration is enabled, construction performs one identity-batch solver calibration after workspace initialization; it does not advance physical state and is not repeated by `reset()`.

**Migration from HEAD**:

Old declaration and defaults at HEAD:

```cpp
struct GpuBackend {
    int device_id = 0;
    int threads_per_block = 512;
    size_t max_batch_sims = 256;
    bool enable_profiling = false; // NVTX range annotation claim
    bool use_persistent = true;
};
```

Current declaration and defaults:

```cpp
struct GpuBackend {
    int device_id = 0;
    int threads_per_block = 512;
    size_t max_batch_sims = 256;
    bool enable_profiling = false; // metadata and host timings; no NVTX promise
    bool use_persistent = true;
    BackendMode backend = BackendMode::Graph;
};
```

| Contract | HEAD | Current E2 contract |
|---|---|---|
| `GpuBackend` declaration/defaults | `device_id=0`, `threads_per_block=512`, `max_batch_sims=256`, `enable_profiling=false`, `use_persistent=true`; no `backend` field | Shared `GpuBackend` preserves `use_persistent=true` and adds `backend=BackendMode::Graph`. The omitted `GpuMultiStageSim` backend uses `multi_stage_default_backend()` and requests `Direct`; a supplied `use_persistent=true` remains an intentional `Persistent` request. An explicit constructor mode overrides that mapping and may request `Graph`, `Direct`, `Fallback`, or `Persistent`. `enable_profiling` is metadata/timing request, not an NVTX claim. |
| Profiling/timing | `enable_profiling` claimed NVTX range annotations and `gpu_time_ms` was easy to read as device timing | `profiling_enabled` records metadata only. `gpu_time_ms` is host wall time for a successful CUDA-backed physical pipeline, including transfers and host orchestration; it is not device-only time. |
| RK4 | The wrapper executed the selected `StepperPolicy`, including the legacy RK4 template instantiation | `GpuSingleStageSim<RK4Stepper>` remains constructible, but `step()` throws `std::logic_error`; it never silently executes Euler. |
| Report access | No unified wrapper report accessor | `execution_report()` returns the cumulative `ExecutionReport` with requested/resolved modes, fallback reasons, timings, calibration metadata, and `gpu_executed`. |
| Thermal resistance access | No wrapper accessor for updated filament resistance | `filament_resistances()` returns a value copy of the physical resistance vector, initialized/reset to armature references and updated after Joule heating; it exposes no raw device pointer. |
| Validation | Legacy construction could defer or omit several boundary checks | Invalid public configuration, null excitation, invalid dimensions/geometry/state, non-finite voltage, invalid launch settings, and negative/out-of-range configured device IDs throw `std::invalid_argument`. Genuine CUDA enumeration/device-selection failures, no-device/insufficient-driver detection, and context/allocation/pipeline failures instead use an honest locked CPU fallback; runtime failures set `runtime_fallback_reason=RuntimeFailure`, while unavailable runtime capability uses `CapabilityUnavailable`. |
| Ownership | Legacy wrapper owned adaptor/persistent resources directly | The wrapper owns `GpuEngine` by `std::unique_ptr`; `Excitation` remains moved-in; no raw device pointers cross the public wrapper. Engine resources use RAII, including partial allocation cleanup. |
| Fallback | Persistent and device failures could be inferred from the request or legacy path | Inspect resolved `backend`, `solver`, `gpu_executed`, `static_fallback_reason`, `runtime_fallback_reason`, and `fallback_reason`. Explicit Fallback and unsupported Graph are CPU-only and do not create a CUDA context. Runtime availability can still force fallback after a capability snapshot; an invalid configured device is rejected before allocation. |
| CPU/GPU parity | Legacy GPU behavior did not expose the migrated engine contract | The successful GPU path and CPU fallback share the Euler physical pipeline and are compared by focused tests, but CUDA/CPU floating-point differences can remain. Graph currently captures only mutual work; matrix, Eigen solver, force/state, and thermal work remain outside the graph. Tests require `gpu_executed` and the expected resolved backend whenever CUDA is available; fallback is tested separately. |

Migration examples:

```cpp
// HEAD: persistent was the implicit default.
GpuBackend old_style{};

// E1: Graph is the explicit default; inspect the report after construction/step.
GpuBackend graph;
graph.backend = BackendMode::Graph;
graph.use_persistent = false;
GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(excitation), 1e-6,
                                   false, GpuOptLevel::Full, graph);
sim.step();
const auto& report = sim.execution_report();
if (report.gpu_executed) {
    // report.gpu_time_ms is host wall time for the committed CUDA pipeline.
}

// Explicit CPU-only migration.
GpuBackend cpu;
cpu.backend = BackendMode::Fallback;
GpuSingleStageSim<EulerStepper> cpu_sim(coil, arm, std::move(cpu_excitation), 1e-6,
                                        false, GpuOptLevel::Full, cpu);
```

**Complete example**:

```cpp
#include <coilgun/coilgun_cuda.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::GpuSingleStageSim;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;

    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6);

    sim.run();
    std::cout << "Muzzle velocity: " << sim.result().summary.muzzle_velocity << " m/s\n";
    std::cout << "Efficiency:      " << sim.result().summary.efficiency * 100 << " %\n";
    return 0;
}
```

The construction and stepping surface is intentionally close to `SingleStageSim`, but the GPU wrapper adds `execution_report()` and currently supports only Euler stepping. Migrating an RK4 CPU simulation requires an explicit CPU implementation until true GPU RK4 is added; the GPU wrapper does not claim RK4 parity.

---

### GpuMultiStageSim

```cpp
#include <coilgun/simulation/cuda/gpu_multi_stage_sim.hpp>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class GpuMultiStageSim {
public:
    static constexpr int kMaxStages = 50;

    GpuMultiStageSim(
        std::vector<components::DrivingCoil>     coils,          // copied
        components::Armature                      armature,      // copied
        std::vector<std::unique_ptr<Excitation>>  excitations,   // moved-in
        std::vector<TriggerConfig>                trigger_configs,
        double                                    dt,
        bool                                      enable_thermal = false,
        GpuOptLevel                               opt_level = GpuOptLevel::Full,
         const GpuBackend&                         backend = multi_stage_default_backend(),
        BackendMode                               explicit_backend = BackendMode::Auto
    );

    ~GpuMultiStageSim();

    GpuMultiStageSim(const GpuMultiStageSim&) = delete;
    GpuMultiStageSim& operator=(const GpuMultiStageSim&) = delete;

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy& p);
    void                    reset();

    const MultiStageResult& result()      const;
    const MultiStageState&  state()       const;
    double  dt()              const;
    int     step_count()      const;
    int     num_stages()      const;
    const ExecutionReport& execution_report() const;
    bool                    graph_assisted() const;
    std::vector<double> filament_resistances() const;
    std::vector<double> mutual_inductances() const;
    std::vector<double> mutual_gradients() const;
};

} // namespace
```

**Constructor constraints**: `coils.size() == excitations.size()`, `trigger_configs.size() == coils.size() - 1`, `coils.size() <= kMaxStages (50)`, `dt` is positive and finite, every excitation is non-null with a finite voltage, every trigger passes `validate_trigger_config()`, and the backend settings are valid. Violating any throws `std::invalid_argument`, including an invalid explicit backend enum value. Finite positions and finite non-negative delays remain trigger-eligible; `+infinity` is the sole valid terminal never-trigger value.

**Constructor parameters**:

| Parameter | Description | Typical value |
|-----------|-------------|---------------|
| `coils` | Driving coil geometry per stage (copied) | — |
| `armature` | Armature geometry and filament discretisation (copied) | — |
| `excitations` | One `Excitation` per stage (moved-in) | `CrowbarExcitation(450.0, 0.001)` each |
| `trigger_configs` | One `TriggerConfig` per stage except stage 0 | `{Position, 0.09}` |
| `dt` | Fixed time step (s) | `1e-6` |
| `enable_thermal` | Enable adiabatic filament heating (CPU-side) | `false` |
| `opt_level` | GPU optimisation level | `GpuOptLevel::Full` |
| `backend` | GPU backend configuration. The omitted argument uses `multi_stage_default_backend()` (`use_persistent=false`, `backend=Direct`). A supplied `use_persistent=true` remains an intentional Persistent request. | `multi_stage_default_backend()` |
| `explicit_backend` | Optional explicit backend override. `Auto` preserves supplied legacy `use_persistent` behavior; otherwise it overrides `use_persistent` and `backend.backend`. | `BackendMode::Auto` |

**Methods** — identical signatures to CPU `MultiStageSim`:

| Method | Behaviour |
|--------|-----------|
| `step()` | Advance one Euler time step. Triggered unfinished stages always participate in their coil circuit solve. In `Full`/`Aggressive`, the distance cutoff applies only to stage-armature mutual terms and force; those terms are zero until the stage re-enters range. `Direct` launches the CUDA mutual segment; `Graph` captures/replays only that mutual segment and performs matrix assembly, Eigen solve, force/state update, and thermal update outside the graph; `Fallback` performs the complete step on CPU. `GpuMultiStageSim<RK4Stepper>::step()` throws `std::logic_error` and never silently executes Euler. |
| `run()` / `run(policy)` | Run to termination. Automatic completion requires every stage to be completed or terminally ineligible through `+infinity`; finite delayed/position-triggered stages keep the run alive until they trigger. The explicit maximum-step, velocity-decay, and optional position-bound criteria can still terminate an otherwise eligible run. Velocity-decay acceleration uses the latest committed post-step recorded force, matching CPU `compute_force(state_)`, not the pre-step force applied during integration. |
| `reset()` | Restore all state to initial conditions. |
| `result()` / `state()` / `dt()` / `step_count()` / `num_stages()` | Query. |
| `execution_report()` | Read the requested/resolved backend, `gpu_executed`, partial-Graph counters, fallback reasons, timing, and cumulative diagnostics. `gpu_executed` is proof of a complete committed CUDA-backed physical step, not merely a request or capability. |
| `graph_assisted()` | True only after a complete CUDA-backed step resolved to `BackendMode::Graph`; it means mutual-segment graph assistance, not a full captured pipeline. |
| `filament_resistances()` | Returns a value copy of the current filament resistance vector. Thermal construction initializes it from armature references, thermal updates change it, and `reset()` restores it. The copy is independent of later wrapper mutation. |
| `mutual_inductances()` / `mutual_gradients()` | Return value copies of the current stage/filament M1 and dM1 buffers in engine coordinates. These are available for numerical parity checks without exposing raw device pointers. |

`MultiStageStep::state.force` is the signed instantaneous net axial force recorded after the state update. `stage_forces` contains the corresponding signed instantaneous per-stage contributions, using the committed post-step currents and the mutual-gradient cache computed at the pre-position step boundary. If excitation advancement finishes a stage during the step, that stage's committed force contribution is zero, matching the CPU reference. The force used to advance velocity is the separately retained pre-step applied force and is not the public recorded force; velocity-decay termination deliberately uses the committed recorded force. `PerStageSummary::max_force` is the peak magnitude of that stage's signed recorded contribution, while `MultiStageSummary::max_force` is the peak magnitude of the signed net recorded force history. Trigger positions are captured at the actual pre-step boundary where the trigger fires and are never reconstructed from a later timestamped history sample.

**Internals**: Computes M1/dM1 matrices for all active stage-filament pairs through `GpuEngine`. With Graph, only this mutual-inductance segment is captured/replayed. Inter-filament and inter-stage assembly, Eigen LDLT solve, force/state orchestration, position/velocity updates, excitation advancement, trigger handling, and thermal updates remain outside the captured graph. The ODE dimension is `n_stages + N_filaments`; inactive stages use mask/identity semantics. `execution_report().backend` and `gpu_executed` must be inspected to distinguish graph-assisted CUDA execution from CPU fallback.

**E2 migration contract**:

| Concern | Legacy/current behavior | E2 behavior |
|---|---|---|
| `use_persistent` | `GpuBackend` default remains `true` for source compatibility; older multi-stage callers used `false` for the per-pair direct/fallback path. | An omitted `GpuMultiStageSim` backend uses `multi_stage_default_backend()` and requests `Direct`. A supplied `use_persistent=false` requests `Direct` unless `backend=Fallback`; a supplied `use_persistent=true` requests `Persistent` and may honestly fall back. Supplying `explicit_backend=Graph`, `Direct`, `Fallback`, or `Persistent` is an intentional override and is visible in `requested_backend`. |
| Graph | A Graph request could be read as a complete pipeline. | Graph is explicitly graph-assisted: only mutual-inductance capture/replay is in the CUDA Graph. `graph_assisted()` and `execution_report()` identify actual successful partial-Graph execution; CPU matrix/solver/force/state/thermal work remains outside it. A resolved `Fallback` with `gpu_executed=false` is not GPU execution. |
| RK4 | CPU `MultiStageSim<RK4Stepper>` performs true four-stage RK4. | `GpuMultiStageSim<RK4Stepper>` remains constructible for source compatibility, but `step()` throws `std::logic_error`. No pseudo-RK4 or silent Euler substitution is performed. Migrate RK4 callers to CPU `MultiStageSim<RK4Stepper>` until a staged engine API exists. |
| Lifetime | Legacy adaptor/persistent resources were wrapper-specific. | `GpuEngine` owns CUDA context, graph/cache, solver, thermal buffers, and state via RAII. `execution_report()` remains a wrapper-owned reference; `filament_resistances()` returns an independent value copy. Reset preserves cumulative report diagnostics while restoring simulation state. |
| Exceptions | Some invalid inputs were deferred to lower layers. | Constructor throws `std::invalid_argument` for inconsistent stage/excitation/trigger counts, empty or oversized stage sets, non-positive/non-finite `dt`, null/non-finite excitations, invalid backend mode or launch settings, and invalid geometry/state dimensions. Runtime CUDA unavailability/failure is reported as locked CPU fallback rather than confused with validation. |
| Parity | GPU-vs-CPU tests could pass while both paths were CPU fallback. | On a CUDA-capable runtime, the primary non-degenerate comparison requires `gpu_executed=true` and the expected resolved backend, and compares currents, positions, velocities, force, stage outputs, and thermal/resistance state. A separate explicit fallback test asserts CPU-only execution. |

**E2 migration example**:

```cpp
// Legacy multi-stage meaning: false selected the direct/per-pair path.
GpuBackend legacy_direct;
legacy_direct.use_persistent = false;
GpuMultiStageSim<EulerStepper> direct(
    coils, armature, std::move(excitations), triggers, 1e-6,
    false, GpuOptLevel::Standard, legacy_direct); // requested Direct

// Intentional Graph request: explicit mode is required to distinguish it
// from the legacy false flag. It is only graph-assisted.
GpuBackend graph_backend;
graph_backend.use_persistent = false;
GpuMultiStageSim<EulerStepper> graph(
    coils, armature, std::move(graph_excitations), triggers, 1e-6,
    false, GpuOptLevel::Standard, graph_backend, BackendMode::Graph);
graph.step();
const auto& report = graph.execution_report();
if (report.gpu_executed && graph.graph_assisted()) {
    // Only the mutual-inductance segment was captured/replayed.
}

// RK4 migration remains on the CPU.
MultiStageSim<RK4Stepper> cpu_rk4(
    coils, armature, std::move(rk4_excitations), triggers, 1e-6);
```

**Complete example**:

```cpp
#include <coilgun/coilgun_cuda.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::GpuMultiStageSim;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;

    DrivingCoil c1(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    std::vector<DrivingCoil> coils = {c1, c2};
    std::vector<std::unique_ptr<Excitation>> excs;
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    excs.push_back(std::make_unique<CrowbarExcitation>(450.0, 0.001));
    std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};

    GpuMultiStageSim<EulerStepper> sim(
        std::move(coils), arm, std::move(excs), triggers, 1e-6);

    sim.run();
    auto& r = sim.result();
    std::cout << "Muzzle velocity: " << r.summary.muzzle_velocity << " m/s\n";
    std::cout << "Efficiency:      " << r.summary.efficiency * 100 << " %\n";
    for (auto& ps : r.summary.per_stage)
        std::cout << "  Stage " << ps.stage_index
                  << ": I_peak=" << ps.peak_current << " A\n";
    return 0;
}
```

**Differences from CPU `MultiStageSim`**:

| Aspect | CPU | GPU |
|--------|-----|-----|
| M/dM computation | OpenMP-parallelised | CUDA kernel (one block per pair) |
| Optimisation level | `OptimizationLevel` (3 levels) | `GpuOptLevel` (3 levels) |
| Adaptive GL order | n_nodes=9 or 4 (distance-based) | n_nodes=9 only |
| Temperature update | OpenMP-parallelised | Serial CPU loop |
| Force summation | OpenMP reduction | Serial loop |
| T(q,p) lookup table | Used when in table range | Construction-time only (not GPU) |

---

### SimBatch

```cpp
#include <coilgun/simulation/cuda/sim_batch.hpp>

namespace coilgun::simulation::cuda {

template<typename SP = EulerStepper>
class SimBatch {
public:
    static constexpr int kMaxStages = 50;

    SimBatch(
        std::vector<components::DrivingCoil> coils,
        components::Armature                  armature,
         int                                   num_sims,
         double                                dt,
         const GpuBackend&                     backend = {},
         BackendMode                           explicit_backend = BackendMode::Auto);

    ~SimBatch();

    SimBatch(const SimBatch&) = delete;
    SimBatch& operator=(const SimBatch&) = delete;

    void set_excitations(
        int sim_id,
        std::vector<std::unique_ptr<Excitation>> excitations,
        std::vector<TriggerConfig>               trigger_configs);

    void run();
    void run(const TerminationPolicy& policy);

    const MultiStageResult& result(int sim_id) const;
    int num_sims() const;
    const ExecutionReport& execution_report() const;
    bool graph_assisted() const;
};

} // namespace
```

`SimBatch` is a container for **parameter sweeps** — running multiple simulations that share identical coil and armature geometry but differ in excitation parameters (voltage, capacitance) and/or trigger positions.

**Execution contract**: `SimBatch` is a `GpuEngine`-backed wrapper. One engine is constructed with `B = num_sims`; every physical buffer uses a fixed row-major layout, including `[B][S][F]` mutual and gradient arrays and `[B][S+F]` current rows. Active masks freeze completed rows in place. Rows are never compacted or reindexed, so `result(sim_id)` remains stable while shorter simulations finish before longer rows. This fixed-capacity policy is intentional: no active-row compaction is performed between steps. An active-list compaction threshold is explicitly deferred and is not implemented in E3; fixed physical and result indices remain stable.

All simulations share **exactly** the same coil geometry and filament discretisation (`m × n`). Per-simulation excitation sources, trigger configurations, and stage voltages are supplied through `set_excitations()`. Circuit masks select stage participation; mutual masks additionally select stage-armature mutual inductance and force. Distant active stages use the canonical cutoff `abs(armature_position - coil.position()) <= 10 * coil.length()` in every resolved backend, including `Fallback`; outside that range their mutual and recorded force terms are zero.

`SimBatch` supports `EulerStepper` only. `SimBatch<RK4Stepper>` is rejected explicitly: construction and `run()` throw `std::logic_error`. It does not claim RK4 parity or silently substitute Euler.

**Constructor parameters**:

| Parameter | Description |
|-----------|-------------|
| `coils` | Shared driving coil geometry (copied). Same for all simulations. |
| `armature` | Shared armature geometry and filament discretisation (copied). |
| `num_sims` | Number of simulations. Must be positive and ≤ `max_batch_sims` in `GpuBackend`. |
| `dt` | Fixed time step (s), shared across all simulations. |
 | `backend` | GPU backend configuration. An explicit constructor `explicit_backend` value is authoritative first; otherwise an explicit `backend.backend` value is authoritative for `Graph`, `Direct`, `Fallback`, and `Persistent`; the legacy flag is consulted only when both are `Auto`. The resolved mode and fallback reason are available from `execution_report()`. |
 | `explicit_backend` | Optional constructor-level backend override. `Auto` preserves the `GpuBackend` resolution; another value overrides both `backend.backend` and `use_persistent`. |

**Methods**:

| Method | Behaviour |
|--------|-----------|
| `set_excitations(sim_id, excitations, triggers)` | Configure excitation sources and trigger settings for simulation `sim_id`. Must be called for each simulation before `run()`. |
| `run()` / `run(policy)` | Run all simulations simultaneously. Each simulation runs independently; the first to reach its termination condition does not block others. Configure every row before the first call: a `SimBatch` run is one-shot, and a later call does not reset or resume the completed histories. |
| `result(sim_id)` | Retrieve `MultiStageResult` for a specific simulation. |
| `num_sims()` | Number of simulations in this batch. |
| `execution_report()` | Return the shared engine's requested/resolved backend and solver, fallback diagnostics, timing metadata, and `gpu_executed`. |
| `graph_assisted()` | True only after a complete CUDA step resolves to `Graph` and actually executes CUDA. Graph assistance covers the mutual-inductance segment only. |

**Solver and fallback**: The wrapper requests `SolverMode::Auto`, allowing `GpuExecutionPlanner` to select `Batched` for a large batch or dimension. If the CUDA context or batched capability is unavailable, the engine reports `Eigen` and executes the CPU fallback. `SimBatch` does not claim cuSOLVER completion merely because `Batched` was requested; inspect `execution_report().solver`, `backend`, `gpu_executed`, and fallback fields.

**Step and history semantics**: Before each engine step, `SimBatch` checks heterogeneous triggers, updates circuit/mutual masks, and captures pre-step currents and the mutual-gradient cache at the pre-position boundary. The engine uses those values for the physical Euler update. After the step, each excitation advances, finished stages are masked, and recorded per-stage forces are recomputed from post-step currents with that pre-position gradient cache. This matches `GpuMultiStageSim` history semantics. The history stores a row for every executed step; active-list compaction is explicitly deferred and not implemented in E3.

**Backend model**: `Direct` uses the direct mutual pipeline when CUDA is available. `Graph` is graph-assisted only: the mutual-inductance segment may be captured/replayed while matrix assembly, solver, force/state orchestration, and thermal work remain outside the graph. `Persistent` is a request that may resolve to `Fallback` because the current synchronous engine does not provide a resident control stream; multi-row batches also cannot resolve to persistent execution. Explicit `Fallback` is CPU-only and does not create a CUDA context.

**Example — parameter sweep over capacitor voltage**:

```cpp
#include <coilgun/coilgun_cuda.hpp>
#include <iostream>
#include <vector>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::cuda::SimBatch;
    using coilgun::simulation::EulerStepper;
    using coilgun::simulation::CrowbarExcitation;
    using coilgun::simulation::TriggerConfig;
    using coilgun::simulation::TriggerMode;
    using coilgun::simulation::TerminationPolicy;

    // Shared geometry
    DrivingCoil c1(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.0);
    DrivingCoil c2(0.01, 0.03, 0.05, 150,
                   COPPER.resistivity_ref, 1e-6, 0.7, 0.10);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    int N = 5;
    SimBatch<EulerStepper> batch({c1, c2}, arm, N, 1e-6);

    // Configure per-simulation excitations
    double voltages[] = {300.0, 350.0, 400.0, 450.0, 500.0};
    for (int i = 0; i < N; ++i) {
        std::vector<std::unique_ptr<Excitation>> excs;
        excs.push_back(std::make_unique<CrowbarExcitation>(voltages[i], 0.001));
        excs.push_back(std::make_unique<CrowbarExcitation>(voltages[i], 0.001));
        std::vector<TriggerConfig> triggers = {{TriggerMode::Position, 0.09}};
        batch.set_excitations(i, std::move(excs), triggers);
    }

    auto pol = TerminationPolicy::defaults();
    pol.max_steps = 500;  // truncate for quick sweep
    batch.run(pol);

    for (int i = 0; i < N; ++i) {
        auto& r = batch.result(i);
        std::cout << "V=" << voltages[i] << "V  →  v="
                  << r.summary.muzzle_velocity << " m/s\n";
    }
    return 0;
}
```

---

### GpuAdaptor (Legacy/Internal Compatibility)

```cpp
#include <coilgun/simulation/cuda/gpu_adaptor.hpp>

namespace coilgun::simulation::cuda {

struct CoilGeo { double ri, re, len, pos; int turns; };
struct FilGeo  { double ri, re, len; };

class GpuAdaptor {
public:
    GpuAdaptor();
    ~GpuAdaptor();
    GpuAdaptor(const GpuAdaptor&) = delete;
    GpuAdaptor& operator=(const GpuAdaptor&) = delete;
    GpuAdaptor(GpuAdaptor&&) noexcept;
    GpuAdaptor& operator=(GpuAdaptor&&) noexcept;

    void setup(const std::vector<DrivingCoil>& coils, const Armature& arm, int n_nodes);
    void setup_batch(const std::vector<DrivingCoil>& coils, const Armature& arm, int num_sims, int n_nodes);

    void upload_separation(const std::vector<double>& seps);
    void download_results(std::vector<double>& M_out, std::vector<double>& dM_out, int n_pairs);
    void upload_batch_separations(const std::vector<double>& seps);
    void download_batch_results(std::vector<double>& M_out, std::vector<double>& dM_out);

    const CoilGeo*  d_coils()   const;
    const FilGeo*   d_fils()    const;
    const double*   d_nodes()    const;
    const double*   d_weights()  const;
    double*         d_results_M();
    double*         d_results_dM();
    double*         d_seps();
    double*         d_batch_seps();
    double*         d_batch_results_M();
    double*         d_batch_results_dM();

    int n_stages() const; int n_fil() const; int n_nodes() const; int batch_size() const;
};

}
```

`GpuAdaptor` is a legacy, move-only device-memory helper retained for internal compatibility. The current `GpuSingleStageSim` and `GpuMultiStageSim` execute through `GpuEngine`; only `SimBatch` and the engine's private legacy persistent-backend internals still use `GpuAdaptor`. Most users must not interact with it directly.

The `setup()`/upload contract below documents the historical adaptor architecture only. It is not the current `GpuEngine` contract, and it does not promise that the current single-stage engine uploads invariant data once or exposes these device buffers.

**`CoilGeo` / `FilGeo`** are packed POD structs for device transfer. They mirror the geometry fields of `DrivingCoil` and `Armature` respectively, flattened for GPU kernel parameter space.

**Legacy `setup()`** allocates and uploads invariant geometry (coils, filaments, GL nodes/weights) to device memory. In the historical adaptor architecture it was called once before per-step operations; this is not a current `GpuEngine` guarantee.

**`setup_batch()`** extends `setup()` with per-simulation batch buffers for `SimBatch`. Allocates `num_sims`-wide separation and result arrays.

**`upload_separation()`** / **`download_results()`** handle per-step host↔device data movement for a single simulation.

**`upload_batch_separations()`** / **`download_batch_results()`** handle per-step host↔device data movement for a batch of simulations.

The device pointer accessors (`d_*()`) return pointers to allocated device memory — used by CUDA kernels directly. Note that `d_batch_*()` pointers are only valid after `setup_batch()` has been called.

**Note**: `GpuAdaptor` is move-only (deleted copy). The destructor frees all device allocations via `cudaFree()`. This section is retained for legacy/internal compatibility and is not the authoritative single-stage API.

---

### Internal Device Headers

The `.cuh` headers define `__host__ __device__` functions used by GPU kernels. They are not intended for user code.

```cpp
#include <coilgun/physics/elliptic.cuh>          // __host__ __device__ elliptic integrals
#include <coilgun/physics/mutual_inductance.cuh> // __host__ __device__ filament M, dM/dz
```

These are only compilable with `nvcc` (`#ifdef __CUDACC__` guarded). They provide GPU-compatible inline implementations of the elliptic integral and filament mutual inductance functions. The `coilgun_cuda.hpp` umbrella header does **not** include these — they are for internal use only.

---

### Legacy Architecture Reference

The following diagram and transfer tables describe historical `GpuAdaptor`-based paths. They are retained for compatibility documentation only and are not the current `GpuEngine` execution flow. The authoritative current single-stage and multi-stage flows are documented above: `Graph` captures/replays only the mutual-inductance segment, while matrix assembly, Eigen solve, force/state update, and thermal work remain synchronous and direct outside the graph. `Fallback` executes the complete physical step on CPU.

**Compute flow per time step**:

```
┌─────────────────────────────────────────────┐
│ Host (CPU)                                   │
│   check_triggers() → extinguish_quiet()       │
│   fill mapped separations / doorbells         │
│   wait for persistent kernel (or launch pairs)│
│   read mapped results (or cudaMemcpy D→H)     │
│   build_system_matrix [L - M_I]              │
│   Eigen LDLT solve → new currents            │
│   compute_force(F = Σ I_d × I_f × dM)        │
│   update velocity / position                 │
│   update capacitor voltage / temperature     │
├─────────────────────────────────────────────┤
│ Device (GPU)                                  │
│   persistent_batch_kernel or                 │
│   mutual_inductance_coil_pair_kernel          │
│     per block: 512 threads × ~13 loops        │
│     each loop: 1 elliptic integral pair       │
│     shared memory tree reduction              │
│     → 1 double M, 1 double dM per pair       │
└─────────────────────────────────────────────┘
```

**Historical data uploaded once at construction** (via legacy `GpuAdaptor`):

| Buffer | Size | Content |
|--------|------|---------|
| `d_coils_` | `n_stages × sizeof(CoilGeo)` | ri, re, length, position, turns per stage |
| `d_fils_` | `N_fil × sizeof(FilGeo)` | ri, re, length per filament ring |
| `d_nodes_` | `9 × sizeof(double)` | Gauss-Legendre quadrature nodes |
| `d_weights_` | `9 × sizeof(double)` | GL quadrature weights |

**Historical data transferred per step**:

| Direction | Size | Content |
|-----------|------|---------|
| H→D | `n_active × N_fil × sizeof(double)` | Armature position → separation values |
| D→H | `n_stages × N_fil × 2 × sizeof(double)` | M1 and dM1 matrices |

### Performance Characteristics

| Scale (S×F) | CPU (16-core) | GPU (RTX 5080) | GPU advantage |
|---|---|---|---|
| 1×10 | 19 s | 16 s | 1.2× |
| 2×10 | 58 s | 52 s | 1.1× |
| 25×45 (typical) | ~5 min | ~2 min (est.) | ~2.5× |
| 50×200 (high-res) | ~2 h | ~15 min (est.) | ~8× |

GPU advantage increases with problem scale because the 4D integration kernel (6561 elliptic integral evaluations per pair) exposes massive parallelism. At small scales, kernel launch overhead and PCIe transfers dominate. At large scales, GPU compute throughput saturates.

### Thread Safety

The GPU classes are **single-threaded** — they do not support concurrent `step()` calls on the same instance. `SimBatch` serializes host-side simulation updates; the public GPU simulation constructors do not expose a `cudaStream_t` parameter, so independent-stream execution is not part of this API.

The underlying CUDA kernel is thread-safe with respect to the host: mapped memory and kernel launches use the default stream, which serialises operations.

### Known Limitations

| Limitation | Detail |
|---|---|
| Adaptive GL order (n_nodes=4/9) | Removed. Using 4 GL nodes on the GPU causes non-deterministic floating-point drift (B1). |
| Thermal mode | `ThermalMode::Cpu` uses the CPU material-table update; `ThermalMode::Gpu` uses the GPU thermal workspace when supported. The resolved mode and `thermal_time_ms` are reported by `ExecutionReport`. |
| Persistent kernel | `GpuSingleStageSim`, `GpuMultiStageSim`, and `SimBatch` resolve through `GpuEngine`; persistent execution is attempted only when the engine capability and determinism policy allow it, otherwise the report identifies CPU fallback. |
| CUDA Graphs | Implemented for the mutual-inductance segment only. Matrix assembly, Eigen solve, force/state update and thermal work remain outside the captured graph; the report therefore identifies the actual Graph path but does not claim a fully captured pipeline. |
| Double buffering | Not implemented. Persistent results are synchronized before the CPU LDLT solve; no CPU/GPU overlap is promised. |

The measurement target `bench_gpu_engine` includes CPU Reference baselines and
records machine-specific wall, solver, thermal, transfer, and Graph capture
observations. See `docs/benchmarks/2026-07-19-unified-gpu-engine.md`.

---

## Build System Reference

### CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `COILGUN_BUILD_TESTS` | `ON` | Build unit and integration tests |
| `COILGUN_BUILD_GENERATOR` | `OFF` | Build the T(q,p) lookup table generator |
| `COILGUN_ENABLE_CUDA` | `OFF` | Build the GPU-accelerated backend (`libcoilgun_cuda.a`) |

### CMake Presets

| Preset | Generator | Build type | Flags | Notes |
|--------|-----------|------------|-------|-------|
| `ninja-debug` | Ninja | Debug | `-march=native` | Tests enabled, compile_commands.json |
| `ninja-release` | Ninja | Release | `-march=native -O3` | — |
| `make-debug` | Unix Makefiles | Debug | `-march=native` | Tests enabled, compile_commands.json |
| `ninja-cuda-debug` | Ninja | Debug | `-march=native` | CUDA enabled, Tests enabled |

### Test Presets

```sh
ctest --preset debug   # run all 18 configured test suites
```
