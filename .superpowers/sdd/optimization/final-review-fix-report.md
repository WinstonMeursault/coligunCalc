# Final Optimization Review Fix Report

Date: 2026-09-08

## Review Findings Addressed

- `GeneticOptimizer` now retries a failed direct batch one candidate at a time,
  retaining input order and converting only unsuccessful singleton retries to
  failed results.
- NSGA-II ranking, non-dominated sorting, merge-and-select, and the
  `select_next_generation` convenience API now receive the configured
  `FeasibilityComparator`. Constraint ordering is applied before Pareto
  dominance, with Pareto used when constraint comparison ties.
- `optimize_single_objective(...)` explicitly forces `SingleObjective`; an
  evaluator that returns multiple objectives now reports `ConfigurationError`.
- `BatchEvaluator::statistics_snapshot()` provides optional cumulative
  statistics without wrapper-specific casts. Optimizer results record the
  configured seed and calculate cache hits, actual GPU fallbacks, and elapsed
  time as a per-run snapshot delta.
- `CoilgunOptimizationProblem` counts a GPU fallback only when the injected GPU
  batch callback fails validation or throws and evaluation continues on CPU.
- Candidates skipped by `max_evaluations` stay `Unevaluated`, increment
  `skipped_due_to_budget`, and do not increment failed-evaluation statistics.
- `docs/API.md` and `docs/API_cn.md` now document the matching public
  statistics and fallback contracts.

## RED Evidence

The initial comparator propagation regression test was added before the
`nsga2_rank` signature changed. It failed to compile as intended:

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug \
  --target test_optimization_single test_optimization_nsga2 -j2
```

The compiler reported that `nsga2_rank` had no comparator parameter.

The final convenience-wrapper regression test also failed before its API fix:

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug \
  --target test_optimization_nsga2 -j2
```

The compiler reported too many arguments to `select_next_generation`, proving
that the wrapper did not accept or forward the comparator.

## GREEN Evidence

After the implementation changes, the focused regression suites passed:

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug \
  --target test_optimization_single test_optimization_nsga2 \
  test_optimization_evaluator test_coilgun_optimization -j2
ctest --test-dir /mnt/data/Project/coligunCalc/build/cpu-debug \
  -R '^(test_optimization_single|test_optimization_nsga2|test_optimization_evaluator|test_coilgun_optimization)$' \
  --output-on-failure
```

Result: 4/4 tests passed.

## Final Validation

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug -j2
ctest --test-dir /mnt/data/Project/coligunCalc/build/cpu-debug \
  -R '^(test_optimization_.*|test_coilgun_optimization)$' --output-on-failure
git diff --check
```

Result: the CPU build completed successfully, all 12 optimization-focused CTest
targets passed, and `git diff --check` reported no whitespace errors. The build
emitted pre-existing `nodiscard` warnings in `test_optimization_variables.cpp`.
