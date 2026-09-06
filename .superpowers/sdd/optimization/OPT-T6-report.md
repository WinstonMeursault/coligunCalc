# OPT-T6 Report: Add single objective optimization

## Scope

Implemented a physics-independent single-objective genetic optimizer over the
existing variable schema, batch evaluator, population/operators, and feasibility
comparator.

- `GeneticOptimizer` accepts either a `BatchEvaluator` or an
  `OptimizationProblem` adapter and uses one ordered batch evaluation per generation.
- Exactly one finite objective is required; objective id and direction are frozen
  after the first valid result, and max/min directions are honored.
- Feasibility-first comparator ordering, configurable elite preservation, seeded
  crossover/mutation, and stable candidate ids are used for generation transitions.
- Results consistently populate `pareto_front` with the best feasible candidate and
  `best_by_objective` for the single objective.
- `TerminationConfig` supports maximum generations, maximum evaluations,
  target value, no-improvement generations, and improvement tolerance.
- The pre-existing `TerminationReason` enum has no `MaxEvaluations` member. To keep
  the T1 result model unchanged, evaluation-budget termination is represented as
  `TerminationReason::MaxGenerations` with an explicit message and exposed through
  `genetic_termination_reason()` as `GeneticTerminationReason::MaxEvaluations`.
- Configuration errors, objective schema errors, malformed/non-finite objectives,
  and all-failed batches terminate without producing a misleading best result.

## TDD evidence

### RED

```text
cmake --preset cpu-debug && cmake --build --preset cpu-debug --target test_optimization_single -j2
```

The focused target initially failed to compile because
`coilgun/optimization/genetic_optimizer.hpp` did not exist.

### GREEN

```text
cmake --build --preset cpu-debug --target test_optimization_single -j2
./build/cpu-debug/tests/test_optimization_single
```

Result: 4 test cases and 27 assertions passed.

## Verification

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug -R 'test_optimization_(types|constraints|variables|operators|evaluator|single)' --output-on-failure
```

Result: 6/6 focused optimization tests passed. The complete CPU-debug build also
completed successfully.

Coverage includes max/min direction, deterministic seed behavior, elite
preservation, target/evaluation/no-improvement termination, exact-one-objective
validation, fixed objective schema and direction, all-failed evaluation handling,
and unified single-objective result fields.
