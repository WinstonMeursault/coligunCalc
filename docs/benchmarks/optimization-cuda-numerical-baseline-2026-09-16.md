# CUDA optimization numerical baseline (2026-09-16)

This is the B-T1 simulation baseline for fixed-geometry candidate batches. It
is evidence that the existing `SimBatch<EulerStepper>` path executes on the
real CUDA device and that candidate rows retain their input order. It is not
an optimizer speedup claim: the workload is intentionally small and the
timings include host orchestration and transfers as reported by the library.

## Reproducibility

- Source revision at measurement: `7e43c4e` (`Plan optimization next phase`).
- Worktree state at measurement: dirty, because this probe and this evidence
  document were uncommitted during measurement.
- Preset: `cuda-debug`; C++ compiler: GCC 16.2.1; CMake: 4.4.3.
- CUDA toolkit/compiler: 13.4 (`nvcc` 13.4.59).
- Device: NVIDIA GeForce RTX 5080 Laptop GPU, compute capability 12.0,
  driver 615.71.09, 16303 MiB.
- Workload: two fixed coils, two-stage Euler integration, 2 axial × 2 radial
  armature filaments, `dt=1e-6`, eight steps, time-delay trigger at `8e-6`.
  Candidate row `i` uses stage-0 voltage `280 + 1.25*i` V and stage-1 voltage
  `215 + 0.875*i` V. Geometry and trigger configuration are identical in all
  rows.
- CUDA batch mode requests `BackendMode::Direct`, CUDA optimization precision
  `Full`, and profiling metadata. Thermal batch mode is disabled because the
  current `SimBatch` constructor has no thermal-enable parameter.
- Thermal representative rows use the existing `GpuMultiStageSim` Full wrapper
  with thermal enabled; its resolved thermal mode is CPU, as recorded below.

## Exact commands

```text
cmake --preset cuda-debug
cmake --build --preset cuda-debug --target test_gpu_optimization_baseline -j2
build/cuda-debug/tests/test_gpu_optimization_baseline -s
ctest --preset cuda-debug -L gpu --output-on-failure
git diff --check
```

The focused probe passed (`1` test, `86` assertions). The full GPU-labelled
suite passed (`18/18`).

## Batch timing

The host-wall column is measured around `SimBatch::run`; `gpu_time_ms` and
`transfer_time_ms` are the report's cumulative host timings for the run.

| Batch | Host wall (ms) | GPU report (ms) | Transfer (ms) | Backend | Solver | Precision | GPU executed |
|---:|---:|---:|---:|---|---|---|---|
| 1 | 6.918381 | 6.842048 | 0.818170 | direct | eigen | full | true |
| 8 | 16.305367 | 16.085729 | 0.720724 | direct | cusolver | full | true |
| 32 | 14.484548 | 13.865495 | 0.819694 | direct | cusolver | full | true |
| 128 | 44.937965 | 43.000193 | 1.000027 | direct | cusolver | full | true |

All four runs report `gpu_executed=true` and `backend=direct`; none used the
fallback backend. Row 0, and rows 7 and 31 where present, are identical across
batch sizes to the printed precision. Row IDs and candidate voltages are
checked against the input index in the focused test.

## Numerical deltas

The existing Full tolerance is used (`relative=1e-4`, `absolute=1e-9`) for all
four CPU Reference ↔ CPU Full metrics and for CUDA Full ↔ CPU Full muzzle
velocity, peak-capacitor voltage (including the initial voltage), and thermal
maximum temperature. CPU Reference and CPU Full agree for all four metrics
within that tolerance; they do not necessarily have identical printed
values. CUDA peak current is an explicit known-baseline classification: the
probe asserts it is outside the normal Full tolerance and that its relative
shortfall remains within the named 5–8% band.

The peak-current metric is the maximum absolute current over the recorded
history. The observed CUDA-minus-CPU-Full deltas are deliberately reported,
not hidden by loosening a tolerance:

| Representative row | Muzzle velocity delta | Peak current delta | Peak current relative delta | Peak capacitor voltage delta | Max temperature delta |
|---:|---:|---:|---:|---:|---:|
| 0 | -2.52e-20 m/s | -0.071024464 A | -6.68480% | 0 V | 0 K |
| 7 | -1.47e-19 m/s | -0.073012205 A | -6.66505% | 0 V | 0 K |
| 31 | -2.12e-20 m/s | -0.079827317 A | -6.60554% | 0 V | 0 K |
| 127 | -2.55e-19 m/s | -0.107087766 A | -6.44851% | 0 V | 0 K |

The current discrepancy is a baseline concern for future optimizer fitness
comparisons. It is systematic across candidates and was not fixed in B-T1;
the probe records it as a passing classification only when it remains outside
the normal tolerance but inside the 5–8% band. The thermal wrapper also
confirms that temperature values match exactly, but thermal execution resolves
to `thermal=cpu` rather than GPU thermal in this existing API.

## Follow-up fix validation (2026-09-16)

The review-fix validation was run from the dirty worktree based on
`1ce2ea5` before amending that commit:

```text
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

The focused output emitted `BASELINE_PEAK_CURRENT_CLASSIFICATION` for every
representative batch and thermal row, each with
`normal_tolerance_match=false` and `baseline_band=0.05..0.08`.

## Raw structured probe output

```text
BASELINE_SCHEMA version=1 workload=two-stage-euler-fixed-geometry dt=1e-6 max_steps=8 geometry_rows=2 geometry_axial=2 trigger_mode=time-delay trigger=8e-6 cuda_opt_level=full backend_request=direct thermal_batch=disabled
BASELINE_METRIC path=cpu_reference batch_size=0 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557012734949e-06 peak_coil_current=1.0624789843990838 peak_cap_voltage=280.00164487908188 max_filament_temperature=0
BASELINE_METRIC path=cpu_full batch_size=0 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557012734949e-06 peak_coil_current=1.0624774567950939 peak_cap_voltage=280.00164487908188 max_filament_temperature=0
BASELINE_METRIC path=cpu_reference batch_size=0 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556237175419e-06 peak_coil_current=1.095449673038043 peak_cap_voltage=288.75169628155311 max_filament_temperature=0
BASELINE_METRIC path=cpu_full batch_size=0 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556237175419e-06 peak_coil_current=1.0954481038741668 peak_cap_voltage=288.75169628155311 max_filament_temperature=0
BASELINE_METRIC path=cpu_reference batch_size=0 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921371819077e-06 peak_coil_current=1.2084920340860146 peak_cap_voltage=318.75187251859762 max_filament_temperature=0
BASELINE_METRIC path=cpu_full batch_size=0 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921371819077e-06 peak_coil_current=1.208490322431081 peak_cap_voltage=318.75187251859762 max_filament_temperature=0
BASELINE_METRIC path=cpu_reference batch_size=0 row=127 candidate_id=127 stage0_voltage=438.75 stage1_voltage=326.125 muzzle_velocity=2.8353606734663848e-06 peak_coil_current=1.6606614782839182 peak_cap_voltage=438.75257746677551 max_filament_temperature=0
BASELINE_METRIC path=cpu_full batch_size=0 row=127 candidate_id=127 stage0_voltage=438.75 stage1_voltage=326.125 muzzle_velocity=2.8353606734663848e-06 peak_coil_current=1.6606591966647728 peak_cap_voltage=438.75257746677551 max_filament_temperature=0
BASELINE_BATCH batch_size=1 host_wall_ms=6.9183810000000001 gpu_time_ms=6.8420480000000001 transfer_time_ms=0.81816999999999995 solver_time_ms=0.1883 thermal_time_ms=0 backend=direct solver=eigen precision=full thermal=disabled device_id=0 threads_per_block=512 gpu_executed=true
BASELINE_METRIC path=cuda_batch_full batch_size=1 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557012735053e-06 peak_coil_current=0.99145299275755816 peak_cap_voltage=280.00164487908188 max_filament_temperature=0
BASELINE_BATCH batch_size=8 host_wall_ms=16.305367 gpu_time_ms=16.085728999999997 transfer_time_ms=0.72072400000000003 solver_time_ms=0 thermal_time_ms=0 backend=direct solver=cusolver precision=full thermal=disabled device_id=0 threads_per_block=512 gpu_executed=true
BASELINE_METRIC path=cuda_batch_full batch_size=8 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557012734697e-06 peak_coil_current=0.99145299275754728 peak_cap_voltage=280.00164487908188 max_filament_temperature=0
BASELINE_METRIC path=cuda_batch_full batch_size=8 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556237173797e-06 peak_coil_current=1.0224358987819888 peak_cap_voltage=288.75169628155311 max_filament_temperature=0
BASELINE_BATCH batch_size=32 host_wall_ms=14.484548 gpu_time_ms=13.865495000000001 transfer_time_ms=0.81969399999999992 solver_time_ms=0 thermal_time_ms=0 backend=direct solver=cusolver precision=full thermal=disabled device_id=0 threads_per_block=512 gpu_executed=true
BASELINE_METRIC path=cuda_batch_full batch_size=32 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557012734697e-06 peak_coil_current=0.99145299275754728 peak_cap_voltage=280.00164487908188 max_filament_temperature=0
BASELINE_METRIC path=cuda_batch_full batch_size=32 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556237173797e-06 peak_coil_current=1.0224358987819888 peak_cap_voltage=288.75169628155311 max_filament_temperature=0
BASELINE_METRIC path=cuda_batch_full batch_size=32 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921371818492e-06 peak_coil_current=1.1286630051520412 peak_cap_voltage=318.75187251859762 max_filament_temperature=0
BASELINE_BATCH batch_size=128 host_wall_ms=44.937964999999998 gpu_time_ms=43.000192999999996 transfer_time_ms=1.0000270000000002 solver_time_ms=0 thermal_time_ms=0 backend=direct solver=cusolver precision=full thermal=disabled device_id=0 threads_per_block=512 gpu_executed=true
BASELINE_METRIC path=cuda_batch_full batch_size=128 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557012734697e-06 peak_coil_current=0.99145299275754728 peak_cap_voltage=280.00164487908188 max_filament_temperature=0
BASELINE_METRIC path=cuda_batch_full batch_size=128 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556237173797e-06 peak_coil_current=1.0224358987819888 peak_cap_voltage=288.75169628155311 max_filament_temperature=0
BASELINE_METRIC path=cuda_batch_full batch_size=128 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921371818492e-06 peak_coil_current=1.1286630051520412 peak_cap_voltage=318.75187251859762 max_filament_temperature=0
BASELINE_METRIC path=cuda_batch_full batch_size=128 row=127 candidate_id=127 stage0_voltage=438.75 stage1_voltage=326.125 muzzle_velocity=2.8353606734661633e-06 peak_coil_current=1.5535714306422321 peak_cap_voltage=438.75257746677551 max_filament_temperature=0
BASELINE_THERMAL comparison=cpu_reference,cpu_full,cuda_full execution=GpuMultiStageSim representative_rows=0,7,31,127
BASELINE_THERMAL_REPORT row=0 host_wall_ms=6.2204259999999998 gpu_time_ms=6.1617950000000006 transfer_time_ms=0.44954 solver_time_ms=0.188217 thermal_time_ms=0.0039410000000000001 backend=direct solver=eigen precision=full thermal=cpu device_id=0 gpu_executed=true
BASELINE_METRIC path=cpu_reference_thermal batch_size=0 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557093438099e-06 peak_coil_current=1.0624789980227298 peak_cap_voltage=280.00164487909632 max_filament_temperature=293.00001038620962
BASELINE_METRIC path=cpu_full_thermal batch_size=0 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557093438099e-06 peak_coil_current=1.0624774704187396 peak_cap_voltage=280.00164487909632 max_filament_temperature=293.00001038620962
BASELINE_METRIC path=cuda_full_thermal batch_size=1 row=0 candidate_id=0 stage0_voltage=280 stage1_voltage=215 muzzle_velocity=1.1547557093437801e-06 peak_coil_current=0.99145300634494571 peak_cap_voltage=280.00164487909632 max_filament_temperature=293.00001038620962
BASELINE_THERMAL_REPORT row=7 host_wall_ms=5.9919830000000003 gpu_time_ms=5.9321320000000002 transfer_time_ms=0.41641599999999984 solver_time_ms=0.17680000000000001 thermal_time_ms=0.0036120000000000002 backend=direct solver=eigen precision=full thermal=cpu device_id=0 gpu_executed=true
BASELINE_METRIC path=cpu_reference_thermal batch_size=0 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556328448939e-06 peak_coil_current=1.0954496879792268 peak_cap_voltage=288.75169628156897 max_filament_temperature=293.00001104549045
BASELINE_METRIC path=cpu_full_thermal batch_size=0 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556328448939e-06 peak_coil_current=1.0954481188153471 peak_cap_voltage=288.75169628156897 max_filament_temperature=293.00001104549045
BASELINE_METRIC path=cuda_full_thermal batch_size=1 row=7 candidate_id=7 stage0_voltage=288.75 stage1_voltage=221.125 muzzle_velocity=1.2280556328447472e-06 peak_coil_current=1.0224359136834131 peak_cap_voltage=288.75169628156897 max_filament_temperature=293.00001104549045
BASELINE_THERMAL_REPORT row=31 host_wall_ms=6.1012230000000001 gpu_time_ms=6.0422269999999996 transfer_time_ms=0.57666100000000009 solver_time_ms=0.17530499999999999 thermal_time_ms=0.0035610000000000004 backend=direct solver=eigen precision=full thermal=cpu device_id=0 gpu_executed=true
BASELINE_METRIC path=cpu_reference_thermal batch_size=0 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921507356033e-06 peak_coil_current=1.2084920541847846 peak_cap_voltage=318.75187251861888 max_filament_temperature=293.00001345988693
BASELINE_METRIC path=cpu_full_thermal batch_size=0 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921507356033e-06 peak_coil_current=1.2084903425298565 peak_cap_voltage=318.75187251861888 max_filament_temperature=293.00001345988693
BASELINE_METRIC path=cuda_full_thermal batch_size=1 row=31 candidate_id=31 stage0_voltage=318.75 stage1_voltage=242.125 muzzle_velocity=1.4964921507355821e-06 peak_coil_current=1.1286630251973462 peak_cap_voltage=318.75187251861888 max_filament_temperature=293.00001345988693
BASELINE_THERMAL_REPORT row=127 host_wall_ms=6.2172840000000003 gpu_time_ms=6.1531029999999998 transfer_time_ms=0.50875199999999987 solver_time_ms=0.17233899999999999 thermal_time_ms=0.0035689999999999997 backend=direct solver=eigen precision=full thermal=cpu device_id=0 gpu_executed=true
BASELINE_METRIC path=cpu_reference_thermal batch_size=0 row=127 candidate_id=127 stage0_voltage=438.75 stage1_voltage=326.125 muzzle_velocity=2.8353607221212117e-06 peak_coil_current=1.660661530700718 peak_cap_voltage=438.75257746683093 max_filament_temperature=293.00002550206153
BASELINE_METRIC path=cpu_full_thermal batch_size=0 row=127 candidate_id=127 stage0_voltage=438.75 stage1_voltage=326.125 muzzle_velocity=2.8353607221212117e-06 peak_coil_current=1.6606592490815615 peak_cap_voltage=438.75257746683093 max_filament_temperature=293.00002550206153
BASELINE_METRIC path=cuda_full_thermal batch_size=1 row=127 candidate_id=127 stage0_voltage=438.75 stage1_voltage=326.125 muzzle_velocity=2.8353607221209572e-06 peak_coil_current=1.5535714829195491 peak_cap_voltage=438.75257746683093 max_filament_temperature=293.00002550206153
```

The corresponding `BASELINE_METRIC` records in the focused test output carry
the complete CPU Reference, CPU Full, CUDA batch Full, and thermal
CPU/CUDA-wrapper values for rows `0`, `7`, `31`, and `127`.
