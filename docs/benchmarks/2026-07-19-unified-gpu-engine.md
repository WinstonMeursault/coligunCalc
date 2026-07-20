# Unified GPU Engine Benchmark

This document records measured results for the unified CUDA engine. Timings are
machine-specific observations, not performance guarantees.

## Method

- CPU rows use `MultiStageSim<EulerStepper>` with the Reference optimization level.
- GPU rows use `SimBatch<EulerStepper>` with the same geometry, excitation values,
  time step, and fixed step count.
- CPU/GPU speedup is `CPU wall time / GPU wall time`. For batch rows the CPU
  baseline is multiplied by the batch size because the CPU reference is run as
  independent simulations.
- GPU timing fields are the `ExecutionReport` cumulative host-wall categories:
  `gpu_time_ms`, `solver_time_ms`, `thermal_time_ms`, and `transfer_time_ms`.
- `Graph` describes the current graph-assisted mutual-inductance segment only.
- `Persistent` is expected to resolve to the documented safe fallback until a
  resident control-stream backend is implemented.
- A blank or fallback row is still a valid observation; it must not be treated
  as a successful GPU timing.

## Reproduce

```sh
cmake --build --preset ninja-cuda-debug --target bench_gpu_engine
./build/ninja-cuda-debug/src/cuda/bench_gpu_engine
```

## Measurement Record

The following table is populated from the target machine's benchmark output.
Append a dated subsection for each hardware or build configuration rather than
overwriting prior observations.

### RTX 5080 Laptop, CUDA Compute Capability 12.0

Measured on 2026-07-21 with the CUDA Debug build. Wall time includes wrapper
orchestration; GPU timing fields are the report values for the executed CUDA
path. CPU rows are independently measured Reference simulations.

| Workload | Requested | Batch | Resolved | Solver | Thermal | CPU/GPU speedup | GPU executed | Wall ms | GPU ms | Solver ms | Thermal ms | Transfer ms | Graph rebuilds |
|---|---|---:|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|
| small-single | CPU | 1 | cpu-reference | eigen | cpu | 1.000 | no | 114.944 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| small-single | Direct | 1 | direct | eigen | disabled | 11.658 | yes | 9.860 | 9.772 | 0.355 | 0.000 | 1.037 | 0 |
| small-single | Graph | 1 | graph | eigen | disabled | 12.287 | yes | 9.355 | 9.284 | 0.334 | 0.000 | 1.231 | 1 |
| small-single | Persistent | 1 | fallback | eigen | disabled | 0.161 | no | 714.376 | 0.000 | 0.376 | 0.000 | 0.000 | 0 |
| small-single | Fallback | 1 | fallback | eigen | disabled | 0.160 | no | 719.898 | 0.000 | 0.359 | 0.000 | 0.000 | 0 |
| medium-multi | CPU | 1 | cpu-reference | eigen | cpu | 1.000 | no | 254.073 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| medium-multi | Direct | 1 | direct | eigen | disabled | 26.250 | yes | 9.679 | 9.572 | 1.854 | 0.000 | 1.357 | 0 |
| medium-multi | Graph | 1 | graph | eigen | disabled | 24.863 | yes | 10.219 | 10.124 | 1.803 | 0.000 | 1.295 | 1 |
| large-single | CPU | 1 | cpu-reference | eigen | cpu | 1.000 | no | 231.163 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| large-single | Direct | 1 | direct | cusolver | disabled | 2.105 | yes | 109.811 | 109.637 | 103.602 | 0.000 | 0.517 | 0 |
| large-single | Graph | 1 | graph | cusolver | disabled | 20.626 | yes | 11.208 | 11.140 | 4.082 | 0.000 | 0.511 | 1 |
| batch-medium | Direct | 1 | direct | eigen | disabled | 26.060 | yes | 9.870 | 9.760 | 1.787 | 0.000 | 1.375 | 0 |
| batch-medium | Direct | 8 | direct | cusolver | disabled | 63.514 | yes | 32.396 | 32.033 | 2.706 | 0.000 | 1.117 | 0 |
| batch-medium | Direct | 32 | direct | cusolver | disabled | 71.586 | yes | 114.972 | 113.755 | 5.946 | 0.000 | 1.827 | 0 |
| batch-medium | Direct | 128 | direct | cusolver | disabled | 76.152 | yes | 432.314 | 426.954 | 10.630 | 0.000 | 0.986 | 0 |
| thermal-engine | CPU thermal | 1 | cpu-reference | eigen | cpu | 1.000 | no | 169.713 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| thermal-engine | GPU Full thermal | 1 | direct | eigen | gpu | 12.676 | yes | 13.388 | 13.338 | 1.552 | 0.508 | 0.725 | 0 |
| thermal-engine | GPU Aggressive thermal | 1 | direct | eigen | gpu | 41.896 | yes | 4.051 | 4.011 | 1.511 | 0.463 | 1.114 | 0 |

The run also measured medium-multi Persistent and Fallback requests at
`0.060x` and `0.060x` respectively; both resolved to CPU fallback. Their
fallback reason was `CUDA runtime resource initialization failed: persistent
protocol requires a dedicated control stream` for Persistent and
`CPU fallback explicitly requested` for Fallback. The same behavior was
observed for the large and thermal wrapper cases. These are fallback timings,
not GPU performance results.

The target emits rows for:

- CPU reference versus Direct, Graph, Persistent, and Fallback requests.
- Small single-stage, medium multi-stage, large high-resolution single-stage,
  and thermal single-stage workloads.
- Direct batch sizes 1, 8, 32, and 128.
- CPU thermal, GPU Full thermal, and GPU Aggressive thermal resolution fields
  through the raw engine thermal contract. `SimBatch` currently has no public
  thermal-mode argument and therefore its wrapper benchmark rows are disabled
  thermal unless the raw engine path is used.

## Precision Error Record

The F1 precision fixture uses the design tolerances: Standard `{relative=5e-6,
absolute=1e-10}`, Full `{1e-4, 1e-9}`, and Aggressive `{1e-2, 1e-8}`. On the
same RTX 5080 Laptop run, the maximum relative errors across the 20-step
single-stage CPU comparison were:

| Precision | M/dM fixture | Coil current | Filament current | Position | Velocity | Force |
|---|---|---:|---:|---:|---:|---:|
| Standard | pass; below tolerance | 2.13e-14 | 4.33e-13 | 1.39e-16 | 2.41e-13 | 2.96e-13 |
| Full | pass; below tolerance | 2.13e-14 | 4.33e-13 | 1.39e-16 | 2.41e-13 | 2.96e-13 |
| Aggressive | pass; below tolerance | 1.69e-05 | 3.73e-04 | 2.12e-11 | 4.91e-05 | 1.25e-04 |

These are observations from the test fixture, not universal error bounds. The
test also compares final muzzle velocity, peak force, temperature, and
temperature-derived filament resistance against the CPU reference.

## Solver Decision Record

| Workload shape | `S+F` / batch | Resolved solver | Evidence |
|---|---:|---|---|
| small single-stage | 11 / 1 | Eigen | `small-single` resolved to Eigen; direct wall time 9.860 ms |
| medium multi-stage | 34 / 1 | Eigen | `medium-multi` resolved to Eigen; direct wall time 9.679 ms |
| large single-stage | 129 / 1 | cuSOLVER batched | `large-single` resolved to cuSOLVER; Graph wall time 11.208 ms |
| medium batch | 34 / 8, 32, 128 | cuSOLVER batched | all batch sizes at or above 8 resolved to cuSOLVER |

The static planner rule remains `batch_size >= 8 || S+F >= 128`; no threshold
was tuned from a single machine's timings. Structured/block solver direction 3
(Schur complement, block elimination, and cuSPARSE evaluation) remains future
work and was not implemented in this phase.

## Planner Thresholds

The current static planner thresholds remain unchanged during the first
measurement pass:

- Large workload: `batch_size >= 8` or `stage_count + filament_count >= 128`.
- Auto solver selects batched/cuSOLVER only for a large workload when the
  capability is available.
- Auto GPU thermal selects GPU thermal only for a large workload when the
  capability is available.

These thresholds may be changed only after measured results demonstrate a
stable benefit on representative workloads. Any change must be recorded with
the before/after benchmark rows and must not be described as a universal
guarantee.
