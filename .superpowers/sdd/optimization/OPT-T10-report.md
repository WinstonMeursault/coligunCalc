# OPT-T10 Report: Add Pareto result selectors

## Scope

Added explicit, non-mutating representative selection to `OptimizationResult`.
The result stores only its existing Pareto front and exposes
`select_representative(const RepresentativeSelector&)`, which returns a copied
candidate in `std::optional`. Empty fronts consistently return `std::nullopt`;
no representative is selected or stored automatically.

The selector API provides:

- `MaxObjective`, honoring each objective's maximize/minimize direction;
- `MinConstraintViolationMargin`, minimizing total normalized hard and soft
  violation while treating unsuccessful candidates as worst;
- `IdealPointDistance`, using per-front min/max normalization;
- `WeightedScore`, using the same normalization and objective directions;
- `LexicographicObjectives`, using explicit objective-id precedence;
- the polymorphic `RepresentativeSelector` interface and a callback adapter for
  application-defined policies.

Equal scores retain Pareto-front order, so selection is deterministic. Invalid
selector configuration or malformed objective data throws
`std::invalid_argument`.

### Pareto normalization invariants

The normalized selectors now treat the first Pareto candidate as the positional
objective schema. Every later candidate must have the same objective count and,
at each position, the same ID and maximize/minimize direction. Reordered
objectives or direction changes throw `std::invalid_argument` instead of being
silently normalized as unrelated quantities.

Normalization scales each objective by the larger absolute endpoint before
subtracting endpoints. This keeps scaled endpoints in `[-1, 1]` and the range
in `[0, 2]`, avoiding the overflow and `NaN` produced by a direct `DBL_MAX -
(-DBL_MAX)` range while retaining neutral `0.5` normalization for zero spans.

## API Integration

`include/coilgun/optimization/types.hpp` receives the member declaration,
`std::optional` include, and selector forward declaration because
`OptimizationResult` already lives there and existing callers include this
header directly. Moving the established result type would create unnecessary
compatibility churn. `result.hpp` is the dedicated result API entry point and
keeps that existing include contract intact. No optimizer behavior or automatic
representative field was added.

## TDD Evidence

### RED

The initial focused build failed at `#include
<coilgun/optimization/result.hpp>` because the result and selector API did not
exist. A subsequent constraint test failed by selecting candidate 3 instead of
candidate 4, demonstrating that soft normalized violations were not yet
included.

### GREEN

The focused selector executable passes 12 test cases and 37 assertions covering
direction handling, deterministic ties, hard and soft constraint violations,
normalization, multi-objective trade-offs, lexicographic ordering, callback and
derived custom selectors, non-mutation, empty fronts, invalid configuration,
the unified single-objective result model, schema mismatch rejection, and
opposite finite extrema.

### Normalization regression fix

The new regression tests first failed: both normalized selectors accepted
reordered IDs and changed directions, while an opposite `DBL_MAX`/
`-DBL_MAX` range overflowed and selected the first candidate as a tie. The
stable scaled normalization and positional schema checks make the same tests
pass for both `IdealPointDistance` and `WeightedScore`.

## Files

- `include/coilgun/optimization/types.hpp`
- `include/coilgun/optimization/result.hpp`
- `include/coilgun/optimization/selectors.hpp`
- `src/optimization/result.cpp`
- `src/optimization/selectors.cpp`
- `src/CMakeLists.txt`
- `tests/test_optimization_selectors.cpp`
- `tests/CMakeLists.txt`
- `.superpowers/sdd/optimization/OPT-T10-report.md`

## Verification

```text
cmake --build --preset cpu-debug --target test_optimization_selectors
./build/cpu-debug/tests/test_optimization_selectors
```

Result: 12 test cases passed, 37 assertions passed.

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug --output-on-failure
```

Result: full CPU build succeeded and 31/31 CPU tests passed.
