# B-T1 report — CPU/CUDA optimization numerical baseline

## Status

`DONE_WITH_CONCERNS`

The requested real-device baseline is implemented and committed. The focused
probe passes and proves CUDA execution for `SimBatch<EulerStepper>` at batch
sizes 1, 8, 32, and 128. A numerical concern remains: the peak absolute coil
current reported by the CUDA path is consistently about 6.4–6.7% below the
CPU Reference/Full value for the representative rows. This task records the
observation and does not change production simulation code, as required.

## Changes

- Added `tests/test_gpu_optimization_baseline.cpp`, a CUDA-labelled validation
  probe using real `SimBatch<EulerStepper>` execution.
- Registered the probe in `tests/CMakeLists.txt` with `gpu`, `slow`, and
  `integration` labels through the existing CUDA test helpers.
- Expanded the probe to assert all four metrics for CPU Reference ↔ CPU Full
  and CUDA Full ↔ CPU Full. CUDA peak current is an executable known-baseline
  classification: it must miss the normal Full tolerance and remain within a
  named 5–8% relative-shortfall band.
- Added versioned evidence at
  `docs/benchmarks/optimization-cuda-numerical-baseline-2026-09-16.md`.
- No production source or public optimization API was changed.

The deterministic fixed geometry has two coils and a 2 × 2 filament armature.
Rows vary only the two Crowbar excitation voltages. The probe uses `dt=1e-6`
and eight Euler steps, with a fixed 8 µs stage trigger. Representative rows are
0, 7, 31, and 127. It includes the initial capacitor voltage when calculating
the peak capacitor metric. Thermal representative runs use the existing
`GpuMultiStageSim` Full wrapper because `SimBatch` currently has no thermal
enablement parameter; the resolved thermal mode is explicitly recorded as
`cpu`.

## Verification

Commands run:

```text
cmake --preset cuda-debug
cmake --build --preset cuda-debug --target test_gpu_optimization_baseline -j2
build/cuda-debug/tests/test_gpu_optimization_baseline -s
ctest --preset cuda-debug -N -L integration
ctest --preset cuda-debug -L integration --output-on-failure
ctest --preset cuda-debug -L gpu --output-on-failure
git diff --check
```

Initial baseline results were 1 test passed, 86 assertions passed, and 18/18
GPU-labelled tests passed. Follow-up fix validation results are recorded at
the end of this report.

## Follow-up fix validation

Commands run after the review fixes:

```text
$ cmake --preset cuda-debug
-- Configuring done
-- Generating done
$ cmake --build --preset cuda-debug --target test_gpu_optimization_baseline -j2
[1/2] Building CXX object tests/CMakeFiles/test_gpu_optimization_baseline.dir/test_gpu_optimization_baseline.cpp.o
[2/2] Linking CXX executable tests/test_gpu_optimization_baseline
$ build/cuda-debug/tests/test_gpu_optimization_baseline -s
[doctest] test cases:   1 |   1 passed | 0 failed | 0 skipped
[doctest] assertions: 154 | 154 passed | 0 failed |
[doctest] Status: SUCCESS!
$ ctest --preset cuda-debug -N -L integration
Test #40: test_gpu_optimization_baseline
Total Tests: 10
$ ctest --preset cuda-debug -L integration --output-on-failure
100% tests passed out of 10
Total Test time (real) = 35.39 sec
$ ctest --preset cuda-debug -L gpu --output-on-failure
100% tests passed out of 18
Total Test time (real) = 17.98 sec
$ git diff --check
(no output; exit status 0)
```

The focused output includes `BASELINE_PEAK_CURRENT_CLASSIFICATION` records
with `normal_tolerance_match=false` and `baseline_band=0.05..0.08` for every
representative batch and thermal row. The full raw structured baseline output
above remains unchanged; the new classification records are the review-fix
evidence.

Results from the original baseline run:

- Focused probe: 1 test passed, 86 assertions passed, 0 failed.
- Full GPU-labelled suite: 18/18 tests passed.
- Every focused batch run asserted `gpu_executed == true` and
  `backend != BackendMode::Fallback`; all resolved to `direct`.
- Batch solver resolved to `eigen` at B=1 and `cusolver` at B=8/32/128.
- `git diff --check`: passed.

## Device and toolchain evidence

- GPU: NVIDIA GeForce RTX 5080 Laptop GPU, compute capability 12.0,
  16303 MiB.
- Driver: 615.71.09.
- CUDA toolkit: 13.4; `nvcc` 13.4.59.
- Compiler: GCC 16.2.1.
- CMake: 4.4.3.
- Measurement source revision: `7e43c4e`.
- Measurement worktree state: dirty while the uncommitted probe/report files
  were being measured.

## Numerical outcome and concern

CPU Reference and CPU Full agree for all four metrics at the short
deterministic horizon to the existing Full comparison tolerance
(`relative=1e-4`, `absolute=1e-9`). CUDA Full agrees with CPU Full for muzzle
velocity, peak capacitor voltage, and maximum filament temperature. The
peak-current metric is calculated as the maximum absolute recorded coil
current, matching the CPU summary convention. CUDA peak current is explicitly
classified as a known discrepancy: the probe asserts that it is outside the
normal Full tolerance while its relative shortfall remains in the named
5–8% baseline band.

For rows 0, 7, 31, and 127, CUDA-minus-CPU-Full peak-current deltas are,
respectively, -0.071024464 A (-6.68480%), -0.073012205 A (-6.66505%),
-0.079827317 A (-6.60554%), and -0.107087766 A (-6.44851%). The difference
is reproducible in both the batch path and the thermal wrapper, while velocity,
capacitor voltage, and temperature remain aligned. This is retained as a
follow-up concern for fitness metrics; no tolerance was loosened and no
production fix was attempted.

## Timing evidence

| Batch | Host wall (ms) | `gpu_time_ms` | `transfer_time_ms` | Backend | Solver | Precision |
|---:|---:|---:|---:|---|---|---|
| 1 | 6.918381 | 6.842048 | 0.818170 | direct | eigen | full |
| 8 | 16.305367 | 16.085729 | 0.720724 | direct | cusolver | full |
| 32 | 14.484548 | 13.865495 | 0.819694 | direct | cusolver | full |
| 128 | 44.937965 | 43.000193 | 1.000027 | direct | cusolver | full |

These are simulation baseline timings only, not an optimizer speedup claim.
The exact command, complete structured output, numerical deltas, and thermal
report details are in the versioned benchmark evidence document.

## Commit

The requested single logical commit subject is:

```text
Establish optimization numerical baseline
```
