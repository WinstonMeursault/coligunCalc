# B-T2 — CUDA batch optimization evaluator

## Status

`DONE`

The CUDA track now has a production `CudaBatchEvaluator` adapter. It decodes
candidate rows, validates and compacts locally invalid rows, executes all valid
rows in one `SimBatch<EulerStepper>` call, restores the original row order, and
converts GPU results through the same `CoilgunOptimizationProblem` metric and
constraint conversion used by the CPU evaluator.

## RED → GREEN → Refactor

### RED

The focused production-path test was added and registered before the evaluator
header or implementation existed. The first CUDA build failed at the intended
missing-production-API boundary:

```text
cmake --preset cuda-debug
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
fatal error: coilgun/optimization/cuda_batch_evaluator.hpp: No such file or directory
```

The RED test covers supported binding execution/order, invalid-row expansion,
unsupported binding/metric rejection, and explicit backend fallback handling.

### GREEN

The minimal implementation added:

- `CudaBatchEvaluator` under `include/coilgun/optimization/` and its CUDA
  implementation under `src/cuda/`;
- fixed-geometry/shared-armature validation and the voltage, positive
  capacitance, and trigger binding matrix;
- candidate width, finiteness, schema decode, excitation, and trigger checks;
- one compacted `SimBatch<EulerStepper>` invocation for valid rows;
- original-order expansion with `Invalid` results retained in place;
- explicit `gpu_executed`/resolved-backend gating, with fallback diagnostics
  and run-local fallback statistics;
- a synchronized by-value execution snapshot containing `ExecutionReport`
  and host/end-to-end time;
- reuse of the CPU problem's canonical result conversion, including initial
  capacitor voltage in peak-voltage calculation;
- CUDA umbrella and target install/export integration.

Focused GREEN command and result:

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
build/cuda-debug/tests/test_cuda_batch_evaluator
  test cases:  6 |  6 passed | 0 failed | 0 skipped
  assertions: 50 | 50 passed | 0 failed |
```

The real-device cases assert `gpu_executed == true`, a non-fallback resolved
backend, finite objectives, and preserved row identity. The fallback case
requests `BackendMode::Fallback` and verifies failed diagnostics rather than a
silent CPU result.

### Refactor

After the first green run, metric/constraint conversion was factored into the
existing `CoilgunOptimizationProblem` implementation and shared by CPU and
CUDA paths. The CUDA-only public header no longer includes CUDA runtime-heavy
headers, so it remains parseable in CPU-only header consumers while only the
CUDA umbrella exports it. The focused suite was rerun green after the
refactor.

### Review-fix TDD cycle

The review regressions were added before the production edits. The
zero/negative voltage test initially showed all four rows as batch-wide
`Failed` results (including missing result metadata), rather than
`[Success, Invalid, Invalid, Success]`. The wrapper composability regression
initially observed 2 evaluations and 2 successes in the run-local collector
for one row. These RED checks were run independently:

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
build/cuda-debug/tests/test_cuda_batch_evaluator \
  --test-case='CUDA batch evaluator rejects zero and negative voltage rows locally'
  test case failed: valid rows were Failed and zero/negative rows were not Invalid

build/cuda-debug/tests/test_cuda_batch_evaluator \
  --test-case='CUDA batch evaluator composes with statistics wrapper without double counting'
  CHECK(collector->snapshot().evaluations == 1): values 2 == 1
  CHECK(collector->snapshot().successful_evaluations == 1): values 2 == 1
```

The CPU Reference comparison regression was green against the pre-fix
implementation and guards terminal velocity plus a peak-voltage constraint
at the existing B-T1 Full tolerance; it intentionally does not compare or
accept `PeakCurrent` as a production constraint.

The fixes enforce the existing `CapacitorExcitation`/`CrowbarExcitation`
contract (`initial_voltage` finite and strictly positive) for fixed
configuration construction and every decoded row. `CudaBatchEvaluator` now
leaves generic evaluation/success/failure counts to
`StatisticsBatchEvaluator`, retaining only its CUDA-owned fallback count;
the synchronized execution snapshot remains available through wrapper
composition.

Review-fix GREEN command and result:

```text
cmake --build --preset cuda-debug --target test_cuda_batch_evaluator -j2
build/cuda-debug/tests/test_cuda_batch_evaluator
  test cases:  9 |  9 passed | 0 failed | 0 skipped
  assertions: 80 | 80 passed | 0 failed |
```

## Support matrix

| Area | B-T2 behavior |
|---|---|
| Geometry | Fixed shared `DrivingCoil` geometry only; no candidate geometry bindings |
| Armature | One fixed shared armature; position, velocity, and mass bindings rejected |
| Candidate bindings | `ExcitationVoltage`, `ExcitationCapacitance` (finite and positive), `TriggerValue` |
| Excitation voltage | Fixed and decoded values must be finite and strictly positive |
| Stepper | `EulerStepper` through `SimBatch`; RK4 is not exposed |
| Optimization level | `Full` only; other CPU-specific levels rejected |
| Thermal | Disabled only; thermal objectives/constraints rejected |
| Peak current | Rejected as an objective/constraint because of the B-T1 CUDA shortfall |
| Objective | Existing single objective ID/direction and terminal-velocity semantics |
| Constraints | Existing IDs, directions, reports, and metric definitions; supported nonthermal metrics |
| Invalid rows | Per-row `Invalid`; valid rows compacted once and expanded in input order |
| Backend fallback | Per-valid-row `Failed` with `gpu_backend_fallback`; no CPU evaluator substitution |

Peak current remains available only as diagnostic metadata inherited from the
canonical result conversion; it is not accepted as a production objective or
constraint by this evaluator.

When composed with `StatisticsBatchEvaluator`, generic evaluations and
success/failure counts are owned by that wrapper. The CUDA evaluator reports
only CUDA-owned fallback increments in its synchronized statistics snapshot.

## Verification

Exact commands and outcomes:

```text
cmake --preset cpu-debug
cmake --build --preset cpu-debug -j2
  completed successfully (101 build steps)
ctest --preset cpu-debug -R 'optimization|coilgun_optimization' --output-on-failure
  100% tests passed out of 14
```

```text
cmake --preset cuda-debug
cmake --build --preset cuda-debug -j2
  completed successfully (82 build steps)
ctest --preset cuda-debug -R 'test_cuda_batch_evaluator|optimization|coilgun_optimization' --output-on-failure
  100% tests passed out of 16
ctest --preset cuda-debug -L gpu --output-on-failure
  100% tests passed out of 19
```

The full GPU-labelled run includes the B-T1 baseline and the new evaluator;
all real-device tests completed without skips or fallback. The evaluator's
focused device snapshot resolved the explicit Direct request to `direct`,
reported `gpu_executed=true`, and recorded positive GPU time and transfer
time. The B-T1 peak-current discrepancy remains intentionally unclaimed.

The final review-fix focused verification repeated the selected CUDA and CPU
optimization suites:

```text
ctest --preset cuda-debug -R 'test_cuda_batch_evaluator|optimization|coilgun_optimization' --output-on-failure
  100% tests passed out of 16
ctest --preset cpu-debug -R 'optimization|coilgun_optimization' --output-on-failure
  100% tests passed out of 14
ctest --preset cuda-debug -L gpu --output-on-failure
  100% tests passed out of 19
```

The final review-fix run also rebuilt both presets successfully before these
tests (`cmake --build --preset cuda-debug -j2` and
`cmake --build --preset cpu-debug -j2`).

```text
cmake --install build/cuda-debug --prefix <temporary-directory>
  completed successfully; installed libcoilgun_cuda.a,
  coilgun_cuda.hpp, cuda_batch_evaluator.hpp, CUDA headers, and the
  CUDAToolkit dependency hook
git diff --check
  no output; exit status 0
```

## Files

- `include/coilgun/optimization/cuda_batch_evaluator.hpp`
- `src/cuda/cuda_batch_evaluator.cu`
- `include/coilgun/optimization/coilgun_problem.hpp`
- `src/optimization/coilgun_problem.cpp`
- `include/coilgun/coilgun_cuda.hpp`
- `src/cuda/CMakeLists.txt`
- `CMakeLists.txt`
- `cmake/coilgunConfig.cmake.in`
- `tests/test_cuda_batch_evaluator.cpp`
- `tests/CMakeLists.txt`

## Concerns and deferrals

- B-T1's reproducible 5–8% CUDA peak-current shortfall remains a documented
  numerical concern; no tolerance or physics was changed.
- Fine-grained per-candidate device-failure isolation, detailed GPU metric
  aggregation, and expanded fallback schemas remain B-T3 work.
- The evaluator is a standalone `BatchEvaluator`; callers connect it to a
  `CoilgunOptimizationProblem` callback using the existing A-track callback
  contract. No new CPU umbrella dependency was introduced.
