# OPT-T2 Report: Implement variable encoding and repair

## Scope

Implemented a physics-independent variable schema and canonical candidate
representation for continuous, integer, and enum variables.

- `VariableSpec` provides validated factory constructors for each variable type.
- `VariableSchema` validates non-empty IDs, finite ordered bounds, enum definitions,
  and duplicate IDs at construction time.
- `repair` returns a new candidate, clips continuous values, rounds and clips
  integer values deterministically, and rejects invalid enum indices.
- `encode` and `decode` use the same canonical positional representation and do
  not mutate their inputs or the schema.

## TDD evidence

### RED

Command:

```text
cmake --preset cpu-debug && cmake --build --preset cpu-debug --target test_optimization_variables -j2
```

Result: failed while compiling `test_optimization_variables.cpp` because
`coilgun/optimization/variables.hpp` did not exist. This was the expected
missing-feature failure before the production implementation.

### GREEN

Command:

```text
cmake --preset cpu-debug && cmake --build --preset cpu-debug --target test_optimization_variables -j2 && ./build/cpu-debug/tests/test_optimization_variables
```

Result: 6 test cases and 19 assertions passed.

## Verification

Commands:

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug -L quick --output-on-failure
```

Result: 22/22 quick tests passed, including the variable and constraints tests.

## Edge-case coverage

Tests cover empty schemas, mixed variable types, clipping at both bounds,
deterministic integer rounding, NaN and signed infinity repair, invalid enum
indices, invalid bounds and enum definitions, duplicate IDs, dimensionality
mismatch, schema immutability, and input immutability.
