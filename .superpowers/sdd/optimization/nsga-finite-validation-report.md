# NSGA-II Finite Validation Fix Report

Date: 2026-09-08

## Finding Fixed

`nsga2_rank` now validates successful candidates before ranking. Objective
values must be finite, and every normalized constraint violation must be finite
and non-negative. Invalid input raises `std::invalid_argument` with a precise
diagnostic. Failed and invalid candidates continue to use the existing
objective-free handling and are not subjected to successful-result validation.

## RED Evidence

The regression tests were added before the validation implementation:

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug \
  --target test_optimization_nsga2 -j2
ctest --test-dir /mnt/data/Project/coligunCalc/build/cpu-debug \
  -R '^test_optimization_nsga2$' --output-on-failure
```

The NaN/Inf objective and NaN/negative normalized-violation assertions failed
because ranking did not throw.

## GREEN Evidence

```sh
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build /mnt/data/Project/coligunCalc/build/cpu-debug -j2
ctest --test-dir /mnt/data/Project/coligunCalc/build/cpu-debug \
  -R '^(test_optimization_.*|test_coilgun_optimization)$' --output-on-failure
git diff --check
```

Result: CPU build succeeded, all 12 optimization-focused tests passed, and the
whitespace check was clean.
