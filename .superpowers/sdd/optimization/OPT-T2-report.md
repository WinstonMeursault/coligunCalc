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

## Follow-up validation hardening

The T3 re-review identified that integer bounds were not required to be
integral, unknown `VariableType` values fell through to enum handling, and the
continuous factory did not enforce the validation claimed above. The focused
follow-up adds schema checks for finite integral integer bounds and known
variable types, and validates continuous/integer factory bounds before
constructing a specification.

### RED

Command:

```text
cmake --build --preset cpu-debug --target test_optimization_variables -j2 && ./build/cpu-debug/tests/test_optimization_variables
```

Result: 4 assertions failed in 3 new test cases because fractional integer
bounds, an unknown variable type, and invalid continuous factory bounds were
accepted.

### GREEN

Command:

```text
cmake --build --preset cpu-debug --target test_optimization_variables -j2 && ./build/cpu-debug/tests/test_optimization_variables
```

Result: 9 test cases and 25 assertions passed.

### Follow-up verification

The focused test and the `quick` CTest label were run after the hardening
change; both passed. This follow-up is committed as `Fix optimization variable
validation`.
