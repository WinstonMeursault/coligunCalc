# OPT-T12 Report: Validate optimization workflow

## Scope

Added an end-to-end constrained coilgun optimization test, an explicit timing
benchmark executable, and a saved benchmark record. The workload uses seed
`20260908`, optimizes excitation voltage for terminal velocity, enforces a hard
0.0095 m/s velocity floor, verifies deterministic feasible output, and rechecks
the selected candidate through `OptimizationLevel::Reference`.

The integration test also exercises the adapter's actual GPU callback boundary:
a mixed failed/successful callback batch preserves the successful candidates,
while a throwing callback activates CPU fallback. A cached wrapper verifies one
underlying evaluation for duplicate inputs, a later cache hit, failed evaluation
statistics, and fallback-context accounting.

## TDD Evidence

### RED

```text
TMPDIR=/dev/shm cmake --preset cpu-debug
TMPDIR=/dev/shm cmake --build --preset cpu-debug --target test_optimization_integration -j2
```

The initial configure failed as expected because `bench_optimization.cpp`,
registered alongside the new integration target, did not exist yet. After the
benchmark shell was added, the first behavioral run failed the candidate-isolation
assertion: a throwing callback correctly caused full CPU fallback. The test was
then split into a mixed-result callback for per-candidate isolation and a separate
throwing callback for fallback.

### GREEN

```text
TMPDIR=/dev/shm cmake --build --preset cpu-debug \
  --target test_optimization_integration bench_optimization -j2
ctest --preset cpu-debug -R '^test_optimization_integration$' --output-on-failure
```

Result: 1/1 CTest target passed; its 2 doctest cases and 35 assertions passed.
The Reference assertion uses `5e-8 + 1e-6 * abs(reference_velocity)`; measured
absolute error was `2.955109129e-08`. The feasible optimum was 0.009764939958 m/s.

## Benchmark Commands

Each actual project preset was configured, and the benchmark target was built
with compiler temporaries redirected to `/dev/shm` because `/tmp` is quota-limited:

```text
TMPDIR=/dev/shm cmake --preset cpu-debug
TMPDIR=/dev/shm cmake --build --preset cpu-debug --target bench_optimization -j2
./build/cpu-debug/tests/bench_optimization

TMPDIR=/dev/shm cmake --preset cpu-release
TMPDIR=/dev/shm cmake --build --preset cpu-release --target bench_optimization -j2
./build/cpu-release/tests/bench_optimization

TMPDIR=/dev/shm cmake --preset cuda-debug
TMPDIR=/dev/shm cmake --build --preset cuda-debug --target bench_optimization -j2
./build/cuda-debug/tests/bench_optimization

TMPDIR=/dev/shm cmake --preset cuda-release
TMPDIR=/dev/shm cmake --build --preset cuda-release --target bench_optimization -j2
./build/cuda-release/tests/bench_optimization
```

Raw environment and timing data are saved in
`docs/benchmarks/optimization-2026-09-08.md`. CPU-release recorded 0.015047905 s
setup, 0.000378017 s first batch, 0.000285401 s steady-state batch, 7 optimizer
evaluations, 1 cache hit, 1 isolated callback failure, 1 fallback, a feasible
0.009764939958 m/s terminal velocity, and a `2.955109129e-08` Reference recheck error.

## Relevant Regression Verification

```text
TMPDIR=/dev/shm cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug \
  -R '^(test_optimization_.*|test_coilgun_optimization)$' --output-on-failure
```

Result: 12/12 optimization-focused tests passed in 0.24 seconds.

## Performance Decision and Concerns

Decision: **accept as the initial benchmark baseline**. No equivalent historical
measurement exists, so regression classification against an older result is
deferred rather than invented. Numerical tolerances were not relaxed.

The RTX 5080 Laptop GPU is present, but `CoilgunOptimizationProblem` currently
offers only an injectable GPU callback and no concrete CUDA optimizer backend.
CUDA-configured benchmarks therefore measure and explicitly label the callback
failure/fallback boundary; they do not claim GPU throughput or speedup. The
integration target is CPU-only, so it does not claim `RESOURCE_LOCK gpu`.

## Review Fix Evidence

The original centered armature workload produced identical Full and Reference
results because both paths selected 9-point quadrature. A retained regression
assertion was added before changing the workload:

```text
REQUIRE( reference_error > 1e-15 )
values: REQUIRE( 0 > 1e-15 )
validation workload must exercise distinct Full and Reference paths
```

Moving the initial armature position to 0.015 m places it beyond one 0.010 m
coil length while remaining inside the ten-length cutoff. `Full` therefore uses
4-point quadrature while `Reference` uses 9 points. The fixed-seed results are:

```text
Full terminal velocity      = 0.009764939958 m/s
Reference terminal velocity = 0.009764969509 m/s
Absolute error              = 2.955109129e-08 m/s
Tolerance                   = 5e-8 + 1e-6 * abs(reference)
```

`first_step_seconds` was renamed to `first_batch_seconds` because the benchmark
times a four-candidate batch, not one simulation step. CMake now injects both
the actual source revision and dirty/clean worktree state; this run reports
`source_revision=635e1b3` and `worktree_state=dirty`.

The adapter's `last_batch_used_fallback()` is asserted independently after an
actual throwing callback. `CachedBatchEvaluator::statistics().fallbacks` is
explicitly asserted as zero for calls whose `EvaluationContext::fallback` is
false; it is a context counter, not proof that the adapter executed fallback.

Review-fix focused verification:

```text
TMPDIR=/dev/shm cmake --build --preset cpu-debug \
  --target test_optimization_integration bench_optimization -j2
./build/cpu-debug/tests/test_optimization_integration
```

Result: 2/2 cases and 35/35 assertions passed. CPU-debug, CPU-release,
CUDA-debug, and CUDA-release benchmark targets all ran and reported the same
nonzero Full/Reference delta. The CPU-only integration test no longer holds the
GPU CTest resource lock.
