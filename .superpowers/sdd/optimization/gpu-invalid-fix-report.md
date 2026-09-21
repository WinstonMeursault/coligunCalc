# Per-Candidate GPU Invalid Handling Fix

## Finding

An `Invalid` result is a candidate-level evaluation outcome, not a malformed
batch protocol. A mixed GPU batch containing one `Invalid` candidate and valid
successes must preserve all statuses and must not trigger CPU fallback.

## TDD Evidence

The regression test was added before the production change. The initial run
failed because `CoilgunOptimizationProblem` marked `Invalid` as malformed,
reported `last_batch_used_fallback() == true`, and incremented the fallback
counter. The implementation then restricted protocol fallback to wrong-size,
`Unevaluated`, or malformed/non-finite successful records.

## GREEN Verification

```text
TMPDIR=/mnt/data/Project/coligunCalc/build/opt-t11-tmp \
cmake --build build/cpu-debug -j2
ctest --test-dir build/cpu-debug \
  -R '^(test_optimization_.*|test_coilgun_optimization)$' --output-on-failure
git diff --check
```

Result: build succeeded, all 12 optimization-focused tests passed, and the
whitespace check was clean. The mixed invalid/success GPU regression now keeps
the invalid result, preserves successful siblings, and reports zero fallback
events.
