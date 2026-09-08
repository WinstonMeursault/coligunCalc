# OPT-T12 Report: Validate optimization workflow

## Scope

Added an end-to-end constrained coilgun optimization test, an explicit timing
benchmark executable, and a saved benchmark record. The workload uses seed
`20260908`, optimizes excitation voltage for terminal velocity, enforces a hard
0.27 m/s velocity floor, verifies deterministic feasible output, and rechecks
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

Result: 1/1 CTest target passed; its 2 doctest cases and 33 assertions passed.
The Reference assertion uses `1e-9 + 1e-8 * abs(reference_velocity)`; measured
absolute error was zero. The feasible optimum was 0.2951739619 m/s.

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
`docs/benchmarks/optimization-2026-09-08.md`. CPU-release recorded 0.015106393 s
setup, 0.006852812 s first batch, 0.0066235342 s steady-state batch, 7 optimizer
evaluations, 1 cache hit, 1 isolated callback failure, 1 fallback, a feasible
0.2951739619 m/s terminal velocity, and zero Reference recheck error.

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
failure/fallback boundary; they do not claim GPU throughput or speedup. GPU-facing
CTest work is serialized with `RESOURCE_LOCK gpu`.
