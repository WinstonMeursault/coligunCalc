# GPU Benchmark Schema and Reproduction Contract

This document defines the output contract of
`tests/bench_gpu_engine.cu`. Timings are machine-specific observations. They
are useful for comparing the same workload on the same class of runner, not as
cross-device performance guarantees.

## Reproduce

Run from the repository root:

```sh
cmake --preset cuda-release
cmake --build --preset cuda-release --target bench_gpu_engine
build/cuda-release/src/cuda/bench_gpu_engine > /tmp/coligunCalc-benchmark-output.md
```

Record the commit from `git rev-parse HEAD`, compiler/CMake preset, operating
system, and relevant environment variables beside the output artifact. On
runners with `nvidia-smi`, also record:

```sh
nvidia-smi --query-gpu=name,driver_version,memory.total \
  --format=csv,noheader
```

The benchmark is not a CTest target. It must be invoked explicitly after the
target is built.

## Measurement Protocol

Each GPU repeat constructs a fresh engine. The fixed CUDA Release protocol is:

- 3 independent repeats;
- `setup`: construction and resource initialization, zero simulation steps;
- `first-step/capture-inclusive`: one physical step;
- `replay-only`: the next single step;
- `warm-up`: 5 steps;
- `steady-state`: 10 measured steps.

The current executable measures a two-stage, 32-filament baseline workload and
the same geometry with GPU thermal processing. Baseline requests are Direct /
Eigen, Graph / Batched, Persistent / Eigen, explicit Fallback / Eigen, and
Graph / Batched with alternating runtime masks. The thermal request is Graph /
Batched / GPU thermal.

## Runtime Metadata

The preamble must be retained with every artifact:

| Field | Meaning |
|---|---|
| `CUDA available` | CUDA runtime availability and device count |
| `GPU` | Device name and compute capability when a device is available |
| `CUDA driver API` | Driver API version reported by the CUDA runtime |
| `CUDA runtime` | Runtime version reported by the CUDA runtime |
| `Build contract` | Build type and phase boundaries used by the run |
| `Benchmark repeats` | Number of independent engine constructions |
| `warm-up steps` | Warm-up iterations excluded from steady-state median |
| `measured steady-state steps` | Iterations used for the steady-state sample |
| `Fixed workloads` | Stage, filament, total-filament, and thermal shape |
| `commit` | Repository commit recorded by the reproduction command |

## CPU Reference Rows

The CPU table is emitted before GPU rows:

```text
| Workload | Requested | Thermal | Phase | Iterations | Wall ms | Per-step ms |
```

CPU rows provide the same phase names for speedup comparison. `setup` has zero
iterations and is not a per-step execution comparison. For batch workloads,
the matching CPU per-step value is multiplied by the requested batch size,
because the CPU reference executes independent simulations.

## GPU Row Schema

The GPU table is emitted with this exact header:

```text
| Workload | Requested backend | Requested solver | Requested thermal | Batch | Runtime mask change | Repeat | Phase | Iterations | Resolved backend | Resolved solver | Precision | Resolved thermal | Wall ms | Per-step ms | Steps/s | GPU ms | Solver ms | Thermal ms | Transfer ms | Graph rebuild delta | Graph rebuild total | Fallback delta | Fallback total | CPU/GPU speedup | GPU executed | Finite | Fallback reason |
```

| Field | Type | Meaning |
|---|---|---|
| `Workload` | string | Fixed workload name, currently `baseline` or `thermal` |
| `Requested backend` | enum | Direct, Graph, Persistent, or explicit Fallback request |
| `Requested solver` | enum | Eigen or Batched request |
| `Requested thermal` | enum | Disabled or GPU thermal request |
| `Batch` | integer | Number of independent rows in the engine |
| `Runtime mask change` | boolean | Whether the stage mask is changed during phases |
| `Repeat` | integer | Independent repeat index, starting at zero |
| `Phase` | enum | setup, first-step/capture-inclusive, replay-only, warm-up, steady-state |
| `Iterations` | integer | Physical steps included in the phase |
| `Resolved backend` | enum | Backend actually selected by the engine |
| `Resolved solver` | enum | Solver actually selected by the engine |
| `Precision` | enum | Resolved precision mode |
| `Resolved thermal` | enum | Thermal path actually selected |
| `Wall ms` | number | Host wall time for the entire phase |
| `Per-step ms` | number | `Wall ms / Iterations`, or wall time for zero-step setup |
| `Steps/s` | number | `1000 / Per-step ms` for a nonzero phase, otherwise zero |
| `GPU ms` | number | Cumulative reported GPU execution time |
| `Solver ms` | number | Cumulative solver time |
| `Thermal ms` | number | Cumulative thermal update time |
| `Transfer ms` | number | Cumulative transfer time |
| `Graph rebuild delta` | integer | Rebuilds observed during this phase |
| `Graph rebuild total` | integer | Rebuild count since engine construction |
| `Fallback delta` | integer | Fallback events observed during this phase |
| `Fallback total` | integer | Fallback events since engine construction |
| `CPU/GPU speedup` | number or `n/a` | Comparable CPU per-step time divided by GPU per-step time |
| `GPU executed` | boolean | True only when the report is a non-fallback GPU execution with no new fallback event |
| `Finite` | boolean | All measured wall and report timing fields are finite |
| `Fallback reason` | string | `none` for no fallback, otherwise the report reason |

## Interpretation Rules

1. Count and compare rows by the complete workload/request/batch/mask/repeat/
   phase key. Do not compare setup with steady-state.
2. A speedup is valid only when `GPU executed=yes`, `Finite=yes`, a matching
   CPU phase exists, and the row is not a fallback. `n/a` is required for
   setup, warm-up, missing CPU comparisons, and fallback rows.
3. `gpu_executed=no` with `finite=yes` is a valid fallback observation, not a
   failed benchmark and not a GPU speedup. The fallback reason must be retained.
4. Graph rebuild totals must remain stable when only a fixed-topology runtime
   mask changes. Mask-update cost must be reported separately from topology capture.
5. Use median and dispersion over the three repeats for summaries. Do not claim
   a universal improvement from one machine, one phase, or one microbenchmark.

## Current Schema Smoke

On 2026-07-24, the CUDA Release executable printed the metadata above and the
documented header. The run emitted 90 GPU execution rows: 48 with
`gpu_executed=yes` and 42 valid fallback/no-execution rows. All 90 rows were
finite. Fallback rows included explicit fallback and unavailable/requested
backend cases; none was used as GPU speedup evidence.

The latest detailed repeated-performance record is
`PerformanceOptimizationEffectiveness.md`. It records the retained benchmark
comparisons, median steady-state times, Graph rebuild behavior, and the RTX
5080 Laptop / CUDA 13.3 environment.

## Artifact Checklist

- Keep the complete stdout artifact, not only a summary table.
- Record commit, preset, compiler, GPU, driver, CUDA runtime, and date.
- Verify the table header before parsing rows.
- Verify row counts, `Finite`, `GPU executed`, and fallback reason fields.
- Report median, sample count, and raw values together.
- Never convert fallback latency into a GPU speedup.
