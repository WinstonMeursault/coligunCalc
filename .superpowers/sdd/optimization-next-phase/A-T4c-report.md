# A-T4c — Define optimization statistics ownership

## Outcome

Optimization result statistics are now run-local value snapshots. Every
`GeneticOptimizer::optimize()` call creates a fresh thread-safe
`EvaluationStatisticsCollector`, passes it through `EvaluationContext`, and
finalizes cache-hit, fallback, and wall-clock fields from that collector and
the run's own timer. The optimizer no longer reads or subtracts evaluator
lifetime snapshots.

## RED/GREEN evidence

- RED: before the implementation, the new collector/context test failed to
  compile because `EvaluationStatisticsCollector` and the third
  `EvaluationContext` field did not exist.
- GREEN: the focused evaluator, optimizer, integration, Coilgun, and new
  statistics tests passed after implementation.

## Ownership and counting semantics

- `EvaluationStatisticsCollector` owns one run's additive evaluator metrics;
  all updates and value snapshots are mutex-protected.
- `EvaluationContext{seed, fallback}` remains source-compatible; the optional
  trailing `statistics` sink carries a run collector through every batch call.
- `OptimizationResult::statistics` remains an independent value owned by the
  result. Seed, submitted evaluation count, status counts, budget skips, and
  evaluated generations are maintained by the optimizer. Cache hits and
  backend fallback events are reported by active built-in contributors into
  the run collector. Elapsed time is measured end-to-end by the optimizer and
  is never derived from nested evaluator timers.
- `StatisticsBatchEvaluator` contributes evaluation/status instrumentation and
  the explicit legacy context fallback marker. `CachedBatchEvaluator` reports
  cache hits only to the run sink, and `CoilgunOptimizationProblem` reports
  each actual malformed/failed GPU batch fallback once. This keeps cache and
  backend metrics from being double-counted by wrappers.

## Shared evaluator/cache compatibility

`StatisticsBatchEvaluator::statistics()`,
`CachedBatchEvaluator::statistics()`, and existing
`statistics_snapshot()` surfaces return by-value snapshots. Their aggregate
state is synchronized. `CoilgunOptimizationProblem` synchronizes its retained
fallback aggregate and legacy last-batch flag. `InMemoryEvaluationCache`
synchronizes concurrent `get`/`put` operations.

## Concurrent-run evidence

The deterministic statistics test uses a `std::barrier` to interleave two
optimizer runs sharing one evaluator and verifies distinct collector-specific
cache-hit/fallback metrics. It also runs two optimizers against one cached
evaluator/cache with a preloaded hit and verifies each run receives exactly
one hit while the delegate is not called. A separate barrier test snapshots a
shared lifetime aggregate while concurrent batches are active.

## TSAN evidence

GCC 16.2.1 supported ThreadSanitizer. A dedicated Debug TSAN build ran the
new concurrency test (6/6 cases, 30/30 assertions) and was repeated 20 times
with `TSAN_OPTIONS=halt_on_error=1`, with no reports. TSAN also ran the
evaluator suite (17/17 cases, 125/125 assertions) and Coilgun suite (19/19
cases, 73/73 assertions), with no reports.

## Verification

- All optimization-named CTest tests: 14/14 passed.
- Full CPU Debug CTest: 35/35 passed.
- `git diff --check`: passed.

## Deferred concerns

- GPU execution/transfer/kernel counters remain deferred to B-T3.
- The retained `last_batch_used_fallback()` compatibility flag is synchronized
  but represents the most recently completed batch globally; per-run fallback
  accounting is intentionally provided by the context collector.

## Review follow-up

Two review findings were fixed in the amended task commit:

- Nested `StatisticsBatchEvaluator` instances now share a scoped fallback
  marker in `EvaluationContext`; only the outermost wrapper forwards the
  legacy `fallback` marker to the run collector, while each wrapper retains
  its own lifetime aggregate compatibility behavior. The marker scope resets
  after each call, so reusing a context counts one marker per call.
- `EvaluationStatisticsCollector::add()` remains additive and no longer
  overwrites `seed`; `set_seed()` is the sole seed owner.

RED evidence: the nested-wrapper regression observed collector fallback counts
of 2 then 4 instead of 1 then 2, and the seed regression observed 202 instead
of the run seed 101. GREEN evidence: both regressions pass in the focused
evaluator executable (19/19 cases, 129/129 assertions).

Exact follow-up verification commands and results:

```text
cmake --build --preset cpu-debug --target test_optimization_evaluator -j2
./build/cpu-debug/tests/test_optimization_evaluator
  test cases: 19 | 19 passed | 0 failed
  assertions: 129 | 129 passed | 0 failed
```

The follow-up focused regression command was:

```text
cmake --build --preset cpu-debug --target test_optimization_evaluator test_optimization_statistics test_optimization_single test_optimization_integration test_coilgun_optimization -j2
ctest --preset cpu-debug -R 'test_optimization_evaluator|test_optimization_statistics|test_optimization_single|test_optimization_integration|test_coilgun_optimization' --output-on-failure
git diff --check
```

It passed all five requested
targets (`test_optimization_evaluator`, `test_optimization_statistics`,
`test_optimization_single`, `test_coilgun_optimization`, and
`test_optimization_integration`): 5/5 CTest tests passed. The amended-code
TSAN command was:

```text
cmake --build build/cpu-tsan --target test_optimization_evaluator test_optimization_statistics -j2
TSAN_OPTIONS='halt_on_error=1' ./build/cpu-tsan/tests/test_optimization_evaluator
TSAN_OPTIONS='halt_on_error=1' ./build/cpu-tsan/tests/test_optimization_statistics
```

It rebuilt and ran the evaluator and statistics targets; evaluator
passed 19/19 cases and 129/129 assertions, and statistics passed 6/6 cases and
30/30 assertions, with `TSAN_OPTIONS=halt_on_error=1` and no reports.

### Re-review follow-up: same-context concurrent calls

The shared context-wide fallback scope was found to undercount when two
top-level calls reused one `EvaluationContext` concurrently. A deterministic
barrier regression was added. RED evidence under the shared scope was:

```text
same-context concurrent wrapper calls each forward one fallback marker
CHECK(collector->snapshot().fallbacks == 2)
values: CHECK(1 == 2)
test cases: 20 | 19 passed | 1 failed
```

The fix uses thread-local per-invocation wrapper-chain scope state. Nested
wrappers in one call share the state; concurrent top-level calls, even with
the same context object, use independent state. The context-wide atomic scope
was removed.

Exact GREEN/TSAN commands and results:

```text
cmake --build --preset cpu-debug --target test_optimization_evaluator test_optimization_statistics test_optimization_single test_optimization_integration test_coilgun_optimization -j2
ctest --preset cpu-debug -R 'test_optimization_evaluator|test_optimization_statistics|test_optimization_single|test_optimization_integration|test_coilgun_optimization' --output-on-failure
git diff --check
  5/5 CTest tests passed

cmake --build build/cpu-tsan --target test_optimization_evaluator test_optimization_statistics -j2
TSAN_OPTIONS='halt_on_error=1' ./build/cpu-tsan/tests/test_optimization_evaluator
TSAN_OPTIONS='halt_on_error=1' ./build/cpu-tsan/tests/test_optimization_statistics
  evaluator: 20/20 cases, 130/130 assertions
  statistics: 6/6 cases, 30/30 assertions
  no ThreadSanitizer reports
```
