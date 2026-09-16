# Optimization Fallback Isolation Fix Report

Date: 2026-09-08

## Findings Fixed

- `CoilgunOptimizationProblem` now treats wrong-size GPU batches, successful
  results with the wrong objective identity/count, empty/non-finite objective or
  constraint fields, `Invalid`, and `Unevaluated` results as malformed. It
  falls back to CPU and increments the fallback counter once per malformed
  batch. Ordinary per-candidate `Failed` results remain valid callback output.
- `StatisticsBatchEvaluator` now uses the same singleton-isolation path as
  `CachedBatchEvaluator`, so nested wrappers retry a thrown batch per candidate
  while preserving order.
- `GeneticOptimizer` retries singleton candidates after both standard and
  unknown batch exceptions. Successful siblings survive; unknown diagnostics
  are emitted only when a singleton retry itself throws an unknown exception.
- NSGA-II validates supplied objective-definition IDs and maximize directions
  against every successful candidate.

## RED Evidence

The four regression tests were added before the fixes and failed as expected:

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug \
  --target test_coilgun_optimization test_optimization_evaluator \
  test_optimization_single test_optimization_nsga2 -j2
ctest --test-dir /mnt/data/Project/coligunCalc/build/cpu-debug \
  -R '^(test_coilgun_optimization|test_optimization_evaluator|test_optimization_single|test_optimization_nsga2)$' \
  --output-on-failure
```

Failures covered missing malformed-GPU fallback, missing nested retry
isolation, all-failed unknown-exception handling, and absent NSGA metadata
validation.

## GREEN Evidence

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug -j2
ctest --test-dir /mnt/data/Project/coligunCalc/build/cpu-debug \
  -R '^(test_optimization_.*|test_coilgun_optimization)$' --output-on-failure
git diff --check
```

Result: CPU build succeeded, all 12 optimization-focused tests passed, and the
whitespace check was clean. The malformed-GPU regression was rerun after the
objective identity validation was tightened and passed.
