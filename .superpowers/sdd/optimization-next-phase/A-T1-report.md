# A-T1 Implementation Report: Cache identity and isolation

## Summary

Added a stable evaluator-owned cache identity containing a non-empty namespace and a version. Identity-aware cache keys now isolate evaluators and invalidate entries when the version changes, while the existing two-argument `make_cache_key` overload remains available through a documented legacy identity.

`StatisticsBatchEvaluator` and `CachedBatchEvaluator` forward the wrapped evaluator's identity, and `CachedBatchEvaluator` uses that identity for every cache lookup and insertion.

## TDD RED

The regression tests were added before production changes. They covered:

- deterministic keys for an unchanged identity, context, and candidate;
- key changes for namespace and version changes;
- rejection of an empty namespace;
- shared-cache isolation between evaluator namespaces;
- shared-cache invalidation between evaluator versions;
- a cache hit when the same identified evaluator repeats a request.

Command:

```text
cmake --build --preset cpu-debug --target test_optimization_evaluator
```

Expected RED output (exit code 1):

```text
[1/2] Building CXX object tests/CMakeFiles/test_optimization_evaluator.dir/test_optimization_evaluator.cpp.o
FAILED: [code=1] tests/CMakeFiles/test_optimization_evaluator.dir/test_optimization_evaluator.cpp.o
tests/test_optimization_evaluator.cpp:151:53: error: expected ')' before 'identity'
tests/test_optimization_evaluator.cpp:156:19: error: 'EvaluationCacheIdentity' does not name a type
tests/test_optimization_evaluator.cpp:174:5: error: 'EvaluationCacheIdentity' does not name a type
tests/test_optimization_evaluator.cpp:214:11: error: 'EvaluationCacheIdentity' does not name a type
ninja: build stopped: subcommand failed.
```

This was the intended failure: the tests described the identity type, evaluator hook, and identity-aware key overload before those APIs existed.

## Minimal implementation

- Added `EvaluationCacheIdentity` with `namespace_id` and `version` values and constructor validation for an empty namespace.
- Added `make_cache_key(identity, variables, context)` with length-delimited identity components.
- Preserved `make_cache_key(variables, context)` by delegating to the stable legacy identity `coilgun.optimization.legacy` version `1`.
- Added a virtual `BatchEvaluator::cache_identity()` with that same source-compatible default.
- Forwarded identities through statistics and cached evaluator decorators.
- Updated `CachedBatchEvaluator` to take one stable identity snapshot per batch and include it in all keys.

## GREEN and refactor

Initial GREEN command:

```text
cmake --build --preset cpu-debug --target test_optimization_evaluator && ./build/cpu-debug/tests/test_optimization_evaluator
```

Output (exit code 0):

```text
[1/8] Building CXX object src/CMakeFiles/coilgun.dir/optimization/cache.cpp.o
[2/8] Building CXX object src/CMakeFiles/coilgun.dir/coilgun.cpp.o
[3/8] Building CXX object src/CMakeFiles/coilgun.dir/optimization/evaluator.cpp.o
[4/8] Building CXX object src/CMakeFiles/coilgun.dir/optimization/coilgun_problem.cpp.o
[5/8] Building CXX object src/CMakeFiles/coilgun.dir/optimization/genetic_optimizer.cpp.o
[6/8] Building CXX object tests/CMakeFiles/test_optimization_evaluator.dir/test_optimization_evaluator.cpp.o
[7/8] Linking CXX static library src/libcoilgun.a
[8/8] Linking CXX executable tests/test_optimization_evaluator
[doctest] test cases:  15 |  15 passed | 0 failed | 0 skipped
[doctest] assertions: 117 | 117 passed | 0 failed |
[doctest] Status: SUCCESS!
```

The refactor stored the wrapped evaluator identity once per batch instead of requesting it once per candidate. The same full evaluator suite remained green afterward.

Final verification command:

```text
cmake --build --preset cpu-debug && ./build/cpu-debug/tests/test_optimization_evaluator && ctest --preset cpu-debug -R optimization --output-on-failure
```

Final output summary (exit code 0):

```text
[doctest] test cases:  15 |  15 passed | 0 failed | 0 skipped
[doctest] assertions: 118 | 118 passed | 0 failed |
[doctest] Status: SUCCESS!
100% tests passed out of 12
Total Test time (real) = 0.34 sec
```

## Files changed

- `include/coilgun/optimization/cache.hpp`
- `include/coilgun/optimization/evaluator.hpp`
- `src/optimization/cache.cpp`
- `src/optimization/evaluator.cpp`
- `tests/test_optimization_evaluator.cpp`
- `.superpowers/sdd/optimization-next-phase/A-T1-report.md`

## Concerns

None. Existing evaluator subclasses remain source-compatible through the documented legacy identity; evaluators that share the legacy default must opt into distinct identities when their result semantics differ.
