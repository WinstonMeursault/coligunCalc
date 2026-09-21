# OPT-T3 Report: Implement objectives and constraints

## Scope

Added physics-independent objective and constraint definitions plus a replaceable feasibility comparator.

- `ObjectiveDefinition` validates finite positive scales, normalizes values, and converts maximize objectives to a minimization-oriented value.
- `ConstraintDefinition` supports equal, less-than/equal, greater-than/equal, and inclusive range relations. Reports contain raw and scale-normalized violations and satisfaction state.
- Hard and soft violations can be aggregated independently; feasibility is determined only by hard constraints.
- `FeasibilityComparator` provides FeasibilityFirst, Penalty, and Lexicographic strategies with deterministic objective tie-breaking.
- Registered the implementation and focused test target in CMake.

## TDD evidence

### RED

Command:

```text
cmake --preset cpu-debug && cmake --build --preset cpu-debug --target test_optimization_constraints -j2
```

The initial focused build was red: the first test draft used braced temporary expressions inside doctest's exception macro, which the compiler rejected. After replacing those with named fixtures, the same target proceeded to the behavioral assertions.

### GREEN

Commands:

```text
cmake --build --preset cpu-debug --target test_optimization_constraints -j2
./build/cpu-debug/tests/test_optimization_constraints
```

Result: 4 test cases and 19 assertions passed.

## Verification

```text
ctest --preset cpu-debug -L quick --output-on-failure
```

All relevant existing and new tests passed. The suite reported one unrelated `Not Run` failure for `test_optimization_variables`, whose executable is not present on this branch (OPT-T2 work is not available in the current checkout).

## Validation behavior

Empty IDs, non-finite or non-positive objective/constraint scales, non-finite bounds/values, and inverted ranges throw `std::invalid_argument`. Comparator penalty weights must be finite and non-negative.

## Review Fix Evidence

The review identified that the original `Lexicographic` comparator collapsed all
hard and soft violations into aggregate sums. Constraint definitions and reports
now carry an explicit integer `priority` (lower values are higher precedence),
and lexicographic comparison walks sorted priority levels independently for hard
then soft constraints. Existing aggregate behavior remains unchanged for the
default priority of zero.

### RED

Command:

```text
cmake --build --preset cpu-debug --target test_optimization_constraints -j2
```

The new priority tests initially failed to compile because the priority field was
not present in `ConstraintDefinition`.

### GREEN

Commands:

```text
cmake --build --preset cpu-debug --target test_optimization_constraints -j2
./build/cpu-debug/tests/test_optimization_constraints
ctest --preset cpu-debug -L quick --output-on-failure
```

Results: 7 focused test cases and 33 assertions passed; the quick suite passed
all 22 tests. Coverage includes NaN/infinite objective values and scales,
NaN/infinite constraint scales, bounds, and values, negative/NaN/infinite
penalty weights, and hard/soft per-priority lexicographic ordering.

## Status Review Fix Evidence

`FeasibilityComparator::compare` now checks `Candidate::evaluation_status` before
constraints or objectives. A `Success` candidate always outranks `Invalid`,
`Failed`, or `Unevaluated` candidates, including under `Penalty`; this prevents
status failures with empty hard-constraint reports from being treated as feasible.

### RED

```text
cmake --build --preset cpu-debug --target test_optimization_constraints -j2
./build/cpu-debug/tests/test_optimization_constraints
```

The status-ordering test failed in all three strategies (6 failed assertions)
before the status gate was added.

### GREEN

```text
cmake --build --preset cpu-debug --target test_optimization_constraints -j2
./build/cpu-debug/tests/test_optimization_constraints
```

Result: 8 test cases and 39 assertions passed.

## Status Semantics Re-review Fix Evidence

Non-success candidates now carry an effective infinite hard violation. A
successful candidate with finite hard violations therefore outranks a failed or
invalid candidate with no reports. When both candidates are non-success and have
the same infinite violation, deterministic status ordering is applied:
`Invalid`, then `Failed`, then `Unevaluated`.

### RED

The mixed status/constraint and invalid-vs-failed tests failed before this change:
the prior status gate treated every non-success candidate as equally worse than a
successful candidate, but did not define ordering among non-success statuses.

### GREEN

```text
cmake --build --preset cpu-debug --target test_optimization_constraints -j2
./build/cpu-debug/tests/test_optimization_constraints
```

Result: 10 test cases and 48 assertions passed.

## Optimization Review Follow-up

Lexicographic comparison now orders `Invalid`, `Failed`, and `Unevaluated`
statuses before comparing hard or soft constraint priority levels when both
candidates are non-successful. This prevents an invalid candidate with a
violation report from being ranked below a failed candidate with no reports.

The focused constraints suite passed with the added invalid-versus-failed
regression test.
