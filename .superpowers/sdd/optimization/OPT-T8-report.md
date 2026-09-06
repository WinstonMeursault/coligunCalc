# OPT-T8 Report: Add automatic strategy routing

## Scope

`GeneticOptimizer` now resolves `SelectionStrategy::Auto` once, immediately
after the first successful objective vector is evaluated:

- one objective routes to the existing single-objective path;
- two or more objectives route to NSGA-II parent/offspring selection;
- explicit `SingleObjective` and `NSGA2` configurations reject incompatible
  objective counts with a configuration error;
- objective count, ids, directions, and finite values are frozen for the run;
  later changes terminate with a configuration error.

The multi-objective path reuses the existing evaluator, population, genetic
operators, feasibility comparator, and `nsga2_select` implementation. It fills
the result Pareto front from the first non-dominated front and does not populate
single-objective `best_by_objective` entries.

## TDD evidence

### RED

Added `tests/test_optimization_routing.cpp` and registered its CMake target.
The new Auto multi-objective test failed before the routing implementation:
the existing optimizer terminated with `EvaluationFailure` instead of routing
to NSGA-II.

### GREEN

```text
cmake --build --preset cpu-debug --target test_optimization_routing -j2
./build/cpu-debug/tests/test_optimization_routing
```

Result: 4 test cases and 25 assertions passed.

The existing single-objective suite was updated to assert the required default
`Auto` to NSGA-II routing and to retain explicit `SingleObjective` mismatch
coverage. The focused optimization suites for types, constraints, variables,
operators, evaluator, single-objective, NSGA-II, and routing all passed.

## Verification

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug --output-on-failure
```

Result: 29/29 CPU tests passed, including physics, single-objective, NSGA-II,
and automatic routing tests.
