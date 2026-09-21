# B-T3 — Integrate GPU optimization metrics

## Status

`DONE`

B-T3 completes the CUDA batch evaluator failure model and exposes run-local GPU
metrics. The default policy is strict. CPU repair is opt-in and distinguishes
per-candidate row repair from whole-batch repair. Local `Invalid` rows remain
in place and are never submitted or repaired; protocol mismatches are
structured failures and are never CPU-repaired.

## RED → GREEN → Refactor

### Recovered RED evidence

The worktree contained a partial interrupted implementation and seven modified
files. The existing focused tests were rebuilt before further edits. The first
rerun failed in two places:

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
build/cuda-debug/tests/test_cuda_batch_evaluator
  test cases: 13 | 11 passed | 2 failed
  backend fallback snapshot: fallbacks 0 != 1
  whole-batch fallback: GPU-eligible rows remained Failed; gpu_fallbacks 0 != 1;
                         cpu_fallback_evaluations 0 != 2
```

The failures showed that an injected non-GPU response with no rows was being
classified as a protocol mismatch before backend fallback handling, and that
the legacy fallback aggregate was no longer updated. The existing strict row,
whole-batch, protocol, invalid-row, and timing tests supplied deterministic
failure seams; no hardware fault was needed.

Two additional RED cycles were recorded while closing the contract:

```text
build/cuda-debug/tests/test_cuda_batch_evaluator --test-case='*strict policy converts*'
  failed: Invalid row status (2) was returned instead of required Failed (3)

build/cuda-debug/tests/test_optimization_statistics --test-case='collector keeps all duration fields*'
  failed: negative/NaN duration deltas were accumulated
```

### GREEN

The implementation now:

- adds `CudaFallbackPolicy` and `CudaFallbackOptions` with strict default,
  per-candidate CPU fallback, and whole-batch CPU fallback;
- injects a deterministic `ExecutionFunction` seam while keeping the default
  path as one `SimBatch<EulerStepper>` invocation;
- compacts locally valid rows, preserves local invalid rows and original order,
  isolates row failures, and attaches structured GPU-failure/CPU-fallback
  metadata and diagnostics;
- treats throw/non-GPU resolution as one failed batch and one backend fallback
  event, while treating count/order mismatch as a non-repairable protocol
  error;
- adds the exact CUDA counters and finite non-negative duration fields to
  `EvaluationStatistics`, the synchronized collector, and
  `OptimizationStatistics`;
- forwards CUDA-owned fields through statistics/cache wrappers without
  duplicate counting and preserves legacy fallback compatibility for older
  contributors;
- documents that `gpu_kernel_seconds` is `ExecutionReport::gpu_time_ms` in
  seconds, i.e. measured CUDA physical-pipeline time rather than a pure kernel
  event timer;
- adds deterministic wrapper/cache/concurrent run-local tests using barriers,
  plus real-device proof of `gpu_executed`, non-Fallback backend, ordered rows,
  and non-negative metrics.

Focused CUDA evaluator result:

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
build/cuda-debug/tests/test_cuda_batch_evaluator
  test cases: 19 | 19 passed | 0 failed | 0 skipped
  assertions: 216 | 216 passed | 0 failed |
```

Focused CPU optimization result:

```text
cmake --build --preset cpu-debug -j2
ctest --preset cpu-debug -R 'optimization|coilgun_optimization' --output-on-failure
  100% tests passed out of 14
```

CUDA optimization and real-device results:

```text
ctest --preset cuda-debug -R 'test_cuda_batch_evaluator|optimization|coilgun_optimization' --output-on-failure
  100% tests passed out of 16
ctest --preset cuda-debug -L gpu --output-on-failure
  100% tests passed out of 19
```

The real-device evaluator test proved `ExecutionReport::gpu_executed == true`,
a resolved backend other than `Fallback`, stable row ordering, and finite
non-negative transfer/kernel/elapsed metric values.

### Refactor and concurrency verification

CPU-side duration sanitization was centralized in the thread-safe collector.
The optimizer finalizer copies only cache/CUDA/backend-owned fields and keeps
generic evaluation/status/generation/wall-time ownership in the optimizer.
The legacy `fallbacks` aggregate is used only as a compatibility fallback when
an older contributor does not populate `gpu_fallbacks`; contributors that
populate both are not double-counted.

TSAN was not a configured preset, so an isolated configuration was created:

```text
cmake -S . -B build/cpu-tsan -G Ninja \
  -DCOILGUN_ENABLE_CUDA=OFF -DCOILGUN_BUILD_TESTS=ON \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS='-fsanitize=thread -fno-omit-frame-pointer -g' \
  -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=thread'
cmake --build build/cpu-tsan --target test_optimization_evaluator test_optimization_statistics -j2
TSAN_OPTIONS='halt_on_error=1' build/cpu-tsan/tests/test_optimization_statistics
  7/7 cases, 38/38 assertions, no reports
TSAN_OPTIONS='halt_on_error=1' build/cpu-tsan/tests/test_optimization_evaluator
  20/20 cases, 130/130 assertions, no reports
```

`test_optimization_statistics` was repeated 20 times under the same TSAN
options; all 20 runs passed with no reports.

`git diff --check` completed with no output (exit status 0).

## Exact metric meanings

| Field | Ownership and meaning |
|---|---|
| `gpu_requested_evaluations` | CUDA-eligible rows submitted after local validation |
| `gpu_executed_evaluations` | Rows with a real GPU-backed result before protocol/CPU repair |
| `gpu_successful_evaluations` | GPU rows accepted as `Success` before CPU fallback |
| `gpu_failed_evaluations` | GPU rows rejected before CPU fallback |
| `cpu_fallback_evaluations` | Candidates actually evaluated by CPU repair |
| `gpu_batches` | CUDA batch invocation attempts |
| `gpu_failed_batches` | Attempts that throw, resolve non-GPU/Fallback, or violate protocol |
| `gpu_fallbacks` | One actual backend fallback event per affected batch |
| `gpu_transfer_seconds` | Actual CUDA report transfer time, converted from milliseconds |
| `gpu_kernel_seconds` | `ExecutionReport::gpu_time_ms`, converted to seconds; measured CUDA physical-pipeline time, not a pure kernel timer |
| `gpu_elapsed_seconds` | Evaluator host elapsed time for CUDA batch attempts |

Cached candidates add cache hits only. Invalid candidates add no GPU request,
execution, batch, or CPU-fallback counts. All duration accumulators reject
non-finite or negative deltas.

## Remaining concerns

- The known B-T1 CUDA `PeakCurrent` discrepancy remains intentionally
  unsupported for production fitness; no physics or tolerance was changed.
- CUDA compilation emits existing Eigen relaxed-constexpr warnings; they do not
  affect the passing build or runtime tests.
- No B-T4 benchmarking was started.

## Review follow-up

The amended task commit addresses two Important findings and one documentation
finding.

### RED

The invalid-only timing regression was added first and reproduced the bug:

```text
build/cuda-debug/tests/test_cuda_batch_evaluator --test-case='*invalid-only batches*'
  failed: gpu_elapsed_seconds 1.323e-05 != 0
          execution_snapshot().host_time_ms 0.01323 != 0
```

The public-injection regression was then added as a compile-time assertion:

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
  failed: static assertion "arbitrary CUDA result injection must not be a public constructor"
```

### GREEN and verification

Invalid-only batches now skip CUDA timing entirely: request, batch, execution,
transfer, kernel, elapsed, and snapshot host timing are all exactly zero. The
arbitrary-result constructor and `ExecutionFunction` alias are private; tests
use only a friend `CudaBatchEvaluatorTestAccess` defined in the test binary.
The installed header therefore exposes only the production constructors.

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
build/cuda-debug/tests/test_cuda_batch_evaluator
  test cases: 21 | 21 passed | 0 failed | 0 skipped
  assertions: 229 | 229 passed | 0 failed |

build/cuda-debug/tests/test_optimization_statistics
  test cases: 7 | 7 passed | 0 failed | 0 skipped
  assertions: 38 | 38 passed | 0 failed |
build/cuda-debug/tests/test_optimization_evaluator
  test cases: 20 | 20 passed | 0 failed | 0 skipped
  assertions: 130 | 130 passed | 0 failed |
```

The privacy follow-up also retained a public production constructor accepting
`CudaFallbackOptions`; a compile attempt before adding it failed with no
matching three-argument constructor. The final public-options regression is
included in the 21-case CUDA run above.

The install/export check was:

```text
bt3_install_dir=$(mktemp -d /tmp/coligun-bt3-install.XXXXXX)
cmake --install build/cuda-debug --prefix "$bt3_install_dir"
  installed libcoilgun_cuda.a, coilgun_cuda.hpp,
  simulation/cuda headers, and optimization/cuda_batch_evaluator.hpp
  exported coilgun::coilgun_cuda and coilgun::coilgun
```

The installed `cuda_batch_evaluator.hpp` inspection showed the three production
constructors under `public:` and the arbitrary-result constructor plus
`ExecutionFunction` under `private:`. The bilingual API docs now describe the
actual CUDA target/header installation and list the CUDA batch adapter in the
umbrella contents. `git diff --check` passed with no output.
