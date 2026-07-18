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
    double force;
    std::vector<double> filament_temperatures;  // only when thermal enabled
};

struct SimSummary {
    double muzzle_velocity;    // m/s
    double total_time;         // s
    double max_force;          // N
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

### Multi-Stage Simulation Results

```cpp
#include <coilgun/simulation/multi_stage_result.hpp>

struct StepSnapshot {
    double time;
    double arm_position;        // m
    double arm_velocity;        // m/s
    double force;               // N
    std::vector<double> filament_currents;
    std::vector<double> filament_temperatures;  // K (thermal mode only)
};

struct MultiStageStep {
    StepSnapshot          state;
    std::vector<double>   cap_voltages;     // [n_stages], 0 for inactive stages
    std::vector<double>   coil_currents;    // [n_stages], 0 for inactive stages
};

struct PerStageSummary {
    int    stage_index;
    double trigger_time;         // s
    double trigger_position;     // m
    double peak_current;         // A
    double max_force;            // N
    double energy_depleted;      // J (capacitor energy consumed)
    int    step_count_active;
};

struct MultiStageSummary {
    double                        muzzle_velocity;      // m/s
    double                        total_time;           // s
    double                        max_force;            // N
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
    bool    enable_profiling  = false; ///< Insert NVTX range annotations for nsight profiling.
    bool    use_persistent    = true;  ///< Use persistent kernel with mapped memory. Falls back to per-pair kernel launches if false or if cudaHostAllocMapped fails.
};

}
```

| Field | Purpose | Default |
|-------|---------|---------|
| `device_id` | Selects which physical GPU to target (relevant on multi-GPU systems). | `0` |
| `threads_per_block` | Requested number of threads per block for the 4D GL integration kernel. The persistent launcher clamps values above 512; callers should use a positive power of two no greater than 512. | `512` |
| `max_batch_sims` | Maximum number of simulations in a `SimBatch`. `SimBatch` rejects a larger `num_sims` with `std::invalid_argument`; the field itself does not allocate buffers. | `256` |
| `enable_profiling` | When true, inserts NVTX push/pop ranges around each `step()` for GPU timeline profiling. | `false` |
| `use_persistent` | When true, launches a persistent kernel at construction and uses mapped memory with doorbell protocol for host-device sync (zero kernel launch overhead per step). Falls back to per-pair kernel launches if `false` or if `cudaHostAllocMapped` fails. Requires compute capability ≥ 6.0. | `true` |

**Example**:

```cpp
coilgun::simulation::cuda::GpuBackend be;
be.device_id = 0;
be.threads_per_block = 1024;  // maximise parallelism
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
        const GpuBackend&                backend = {}
    );

    const SimStep&   step();
    const SimResult& run();
    const SimResult& run(const TerminationPolicy& p);
    void             reset();

    const SimResult& result()     const;
    const SimState&  state()      const;
    double  dt()         const;
    int     step_count() const;
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

**Methods**:

| Method | Behaviour |
|--------|-----------|
| `step()` | Advance one time step. Launches CUDA kernels per filament pair, solves the ODE via Eigen LDLT, updates kinematics, records diagnostics. Returns the newly recorded `SimStep`. |
| `run()` | Run to default termination. Blocks until the excitation finishes (crowbar diode ON and current decayed). |
| `run(policy)` | Run to custom `TerminationPolicy`. |
| `reset()` | Restore to initial state. Resets all currents to zero, restores armature position/velocity to construction values, resets the excitation, clears result history. |
| `result()` | Result after the last `run()`. Contains `history` (vector of SimStep) and `summary` (SimSummary). |
| `state()` | Current internal state — live reference. |
| `dt()` | Fixed time step. |
| `step_count()` | Steps since construction or last `reset()`. |

**Internals**: The class internally owns a `GpuAdaptor` that uploads coil geometry, filament discretisation, and GL quadrature nodes to device memory once at construction. Each `step()` call transports the armature position to the device, launches the 4D integration kernel per (coil, filament) pair, and reads back `M1`/`dM1` matrices.

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

The public API is identical to `SingleStageSim`. The only change needed to migrate from CPU to GPU is the class name and include path.

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
        const GpuBackend&                         backend = {}
    );

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy& p);
    void                    reset();

    const MultiStageResult& result()      const;
    const MultiStageState&  state()       const;
    double  dt()              const;
    int     step_count()      const;
    int     num_stages()      const;
};

} // namespace
```

**Constructor constraints**: `coils.size() == excitations.size()`, `trigger_configs.size() == coils.size() - 1`, `coils.size() <= kMaxStages (50)`. Violating any throws `std::invalid_argument`.

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
| `backend` | GPU backend configuration | `{}` (defaults) |

**Methods** — identical signatures to CPU `MultiStageSim`:

| Method | Behaviour |
|--------|-----------|
| `step()` | Advance one time step. Checks triggers, extinguishes quiet stages, launches CUDA kernels for all active (stage, filament) pairs, solves ODE, updates kinematics. |
| `run()` / `run(policy)` | Run to termination. |
| `reset()` | Restore all state to initial conditions. |
| `result()` / `state()` / `dt()` / `step_count()` / `num_stages()` | Query. |

**Internals**: Computes M1/dM1 matrices for all active stage-filament pairs using CUDA kernels. Inter-filament mutual inductance matrix `[M]` is precomputed on CPU at construction (it does not change during simulation). Inter-coil mutual inductance `M_cc` is also precomputed on CPU and populates the off-diagonal coil-coil block of the system matrix. The LDLT linear solve, force computation, and kinematic updates remain on CPU. The ODE dimension is `n_stages + N_filaments`. Inactive stages have identity rows/columns.

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
        const GpuBackend&                     backend = {});

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
};

} // namespace
```

`SimBatch` is a container for **parameter sweeps** — running multiple simulations that share identical coil and armature geometry but differ in excitation parameters (voltage, capacitance) and/or trigger positions.

**Constraints**: All simulations must share **exactly** the same coil geometry and filament discretisation (`m × n`). The `coils` and `armature` passed to the constructor are shared across all simulations. Per-simulation excitation sources and trigger configurations are set via `set_excitations()`.

**Constructor parameters**:

| Parameter | Description |
|-----------|-------------|
| `coils` | Shared driving coil geometry (copied). Same for all simulations. |
| `armature` | Shared armature geometry and filament discretisation (copied). |
| `num_sims` | Number of simulations. Must be positive and ≤ `max_batch_sims` in `GpuBackend`. |
| `dt` | Fixed time step (s), shared across all simulations. |
| `backend` | GPU backend configuration. |

**Methods**:

| Method | Behaviour |
|--------|-----------|
| `set_excitations(sim_id, excitations, triggers)` | Configure excitation sources and trigger settings for simulation `sim_id`. Must be called for each simulation before `run()`. |
| `run()` / `run(policy)` | Run all simulations simultaneously. Each simulation runs independently; the first to reach its termination condition does not block others. |
| `result(sim_id)` | Retrieve `MultiStageResult` for a specific simulation. |
| `num_sims()` | Number of simulations in this batch. |

**Execution model**: `SimBatch` internally manages one `GpuAdaptor` and an array of per-simulation state objects. With `GpuBackend::use_persistent=true` (the default), one persistent CUDA kernel services the fixed `(stage, filament)` pair index space and mapped host buffers carry each simulation's separations and results. The host still checks triggers, solves the ODE, updates kinematics, and services simulations serially. With `use_persistent=false`, it falls back to the per-pair kernel launch path.

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

**Note**: The fallback implementation launches one CUDA kernel per (stage, filament) pair per simulation per step. The default persistent path launches its worker kernel once at construction and uses mapped-memory doorbells thereafter, avoiding per-step kernel launches. Persistent mode is experimental and can be disabled through `GpuBackend::use_persistent`.

---

### GpuAdaptor (Advanced)

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

`GpuAdaptor` manages device memory for GPU-accelerated simulations. It is used internally by `GpuSingleStageSim`, `GpuMultiStageSim`, and `SimBatch`. Most users do not need to interact with it directly.

**`CoilGeo` / `FilGeo`** are packed POD structs for device transfer. They mirror the geometry fields of `DrivingCoil` and `Armature` respectively, flattened for GPU kernel parameter space.

**`setup()`** allocates and uploads invariant geometry (coils, filaments, GL nodes/weights) to device memory. Must be called once before any per-step operations.

**`setup_batch()`** extends `setup()` with per-simulation batch buffers for `SimBatch`. Allocates `num_sims`-wide separation and result arrays.

**`upload_separation()`** / **`download_results()`** handle per-step host↔device data movement for a single simulation.

**`upload_batch_separations()`** / **`download_batch_results()`** handle per-step host↔device data movement for a batch of simulations.

The device pointer accessors (`d_*()`) return pointers to allocated device memory — used by CUDA kernels directly. Note that `d_batch_*()` pointers are only valid after `setup_batch()` has been called.

**Note**: `GpuAdaptor` is move-only (deleted copy). The destructor frees all device allocations via `cudaFree()`.

---

### Internal Device Headers

The `.cuh` headers define `__host__ __device__` functions used by GPU kernels. They are not intended for user code.

```cpp
#include <coilgun/physics/elliptic.cuh>          // __host__ __device__ elliptic integrals
#include <coilgun/physics/mutual_inductance.cuh> // __host__ __device__ filament M, dM/dz
```

These are only compilable with `nvcc` (`#ifdef __CUDACC__` guarded). They provide GPU-compatible inline implementations of the elliptic integral and filament mutual inductance functions. The `coilgun_cuda.hpp` umbrella header does **not** include these — they are for internal use only.

---

### Architecture

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

**Data uploaded once at construction** (via `GpuAdaptor`):

| Buffer | Size | Content |
|--------|------|---------|
| `d_coils_` | `n_stages × sizeof(CoilGeo)` | ri, re, length, position, turns per stage |
| `d_fils_` | `N_fil × sizeof(FilGeo)` | ri, re, length per filament ring |
| `d_nodes_` | `9 × sizeof(double)` | Gauss-Legendre quadrature nodes |
| `d_weights_` | `9 × sizeof(double)` | GL quadrature weights |

**Data transferred per step**:

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
| Thermal mode | Temperature updates run on CPU (single-threaded). The computation is <1% of runtime. |
| Persistent kernel | Experimental. Uses mapped memory with doorbell protocol for host-device sync. It is used by all three GPU simulation classes by default; `use_persistent=false` selects the fallback. Allocation failure falls back to per-pair launches. |
| CUDA Graphs | Not implemented. Would require fixed n_active per step (no dynamic triggering/termination). |
| Double buffering | Not implemented. Persistent results are synchronized before the CPU LDLT solve; no CPU/GPU overlap is promised. |

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
