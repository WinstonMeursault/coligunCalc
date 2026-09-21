# A-T3b — Align optimization constraint semantics

## Scope

Aligned `FeasibilityStrategy::Penalty` with the additive design contract. Hard
feasibility and evaluation status remain ordered before objective scoring;
successful hard-infeasible candidates are ordered by aggregate normalized hard
violation. Soft violations are added to normalized/oriented objective values
with the configured weight. NSGA-II uses the same penalized objective view for
dominance and crowding distance.

## TDD RED evidence

Focused tests were added before production edits. The first build exposed the
missing explicit-objective comparator overload:

```text
cmake --build --preset cpu-debug --target test_optimization_constraints test_optimization_nsga2 -j2
error: no matching function for call to FeasibilityComparator::better(..., ObjectiveDefinition)
```

The NSGA-II target was then built and run against the old implementation:

```text
cmake --build --preset cpu-debug --target test_optimization_nsga2 -j2
ctest --preset cpu-debug -R '^test_optimization_nsga2$' --output-on-failure
```

The expected behavior failures were:

```text
NSGA-II penalty applies soft violation to every Pareto objective: 2 fronts, expected 1
NSGA-II crowding uses penalized objective values: 2.0, expected penalized crowding value
```

These failures demonstrate that soft constraints were being compared as a
constraint prefilter and were absent from the objective/crowding view.

## Review-fix evidence

The compatibility regression was added before the production fix. Building the
focused targets against the four-parameter-only declaration failed at the
legacy function-pointer binding:

```text
error: invalid conversion from ‘std::vector<double> (*)(..., const FeasibilityComparator&)’
to ‘LegacyCrowdingFunction’ {aka ‘std::vector<double> (*)(..., const std::vector<ObjectiveDefinition>&)’}
```

`crowding_distances` now has a real three-parameter overload with the original
default definitions argument. It forwards to the distinct four-parameter
comparator-aware overload. The regression binds that legacy function pointer,
invokes it, and checks the resulting crowding values.

The same-status Lexicographic regression covers hard-priority ordering before
soft-priority ordering for failed candidates. Cross-status ordering remains
status-rank-first.

The second re-review found that the compatibility fix still returned early
after same-status non-success constraint handling. That removed the historical
raw-objective tie-break for non-success candidates whose evaluations retained
objective values. A regression covering Invalid, Failed, and Unevaluated under
all three strategies failed before the production change. The comparator now
falls through to the raw oriented objective after same-status handling; Penalty
does not add soft violation for non-success candidates, matching the prior
behavior, while Lexicographic constraint priorities still take precedence.

## Semantic decisions

- Status ordering is deterministic for all strategies: Success, Invalid,
  Failed, Unevaluated. This makes every successful candidate—including a hard-
  infeasible one—better than non-successful candidates.
- For successful candidates, hard feasibility is checked first; when both are
  hard-infeasible, aggregate normalized hard violation is compared before
  objective scoring. Lexicographic hard-priority comparison remains as a
  tie-break after that required aggregate ordering.
- `Penalty` uses `oriented(value / scale) + penalty_weight * aggregate_soft`.
  Objective definitions can be supplied per comparison or stored in the
  comparator constructor; the existing metadata-driven scale-1 behavior stays
  the default.
- Cleared objective arrays perform constraint-only comparisons, so Penalty does
  not turn soft constraints into a prefilter.
- NSGA-II adds the same soft penalty to every normalized/oriented objective and
  uses those values for both Pareto dominance and crowding. The existing
  crowding API remains source-compatible, with an optional comparator added.
- Normalized constraint violations are validated as finite and non-negative;
  penalty weights retain finite/non-negative validation.

## GREEN verification

Focused comparator and NSGA-II tests:

```text
cmake --build --preset cpu-debug --target test_optimization_constraints test_optimization_nsga2 -j2
ctest --preset cpu-debug -R '^(test_optimization_constraints|test_optimization_nsga2)$' --output-on-failure
100% tests passed out of 2
```

All optimization-named CTest tests after a full CPU Debug rebuild:

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug -R optimization --output-on-failure
100% tests passed out of 12
```

Diff hygiene:

```text
git diff --check
clean
```

## Changed files

- `include/coilgun/optimization/comparator.hpp`
- `src/optimization/comparator.cpp`
- `include/coilgun/optimization/nsga2.hpp`
- `src/optimization/nsga2.cpp`
- `tests/test_optimization_constraints.cpp`
- `tests/test_optimization_nsga2.cpp`
- `.superpowers/sdd/optimization-next-phase/A-T3b-report.md`

## Deferred concern

The scalar comparator accepts one explicit `ObjectiveDefinition`; multiobjective
callers should continue to pass definitions to NSGA-II, where the full fixed
objective schema is validated and applied.
