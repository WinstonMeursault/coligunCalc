# GPU Benchmark Schema and Reproduction Contract

This document defines `gpu-benchmark/v2`, the output contract of
`tests/bench_gpu_engine.cu`. Timings are machine-specific observations, not
cross-device performance guarantees.

## Reproduce

Run from the repository root. Every retained artifact must set the commit
explicitly and must pass the C++20 contract checker.

```sh
cmake --preset cuda-release
cmake --build --preset cuda-release --target bench_gpu_engine
c++ -std=c++20 -O2 -Wall -Wextra -Wpedantic \
  tests/check_gpu_benchmark_schema.cpp -o /tmp/check_gpu_benchmark_schema

for run in 1 2 3; do
  COILGUN_BENCHMARK_COMMIT="$(git rev-parse HEAD)" \
    build/cuda-release/src/cuda/bench_gpu_engine \
    > "/tmp/coligunCalc-b0-t1-run-${run}.md"
done

/tmp/check_gpu_benchmark_schema /tmp/coligunCalc-b0-t1-run-{1,2,3}.md
```

The benchmark is not a CTest target. Keep all three raw stdout files. Also
record the operating system and relevant environment variables. On runners
with `nvidia-smi`, retain:

```sh
nvidia-smi --query-gpu=name,driver_version,memory.total \
  --format=csv,noheader
```

## Measurement Protocol

Each executable run contains three independent repeats. Each GPU repeat
constructs a fresh engine and emits five non-overlapping measurement windows:

| `Measurement window` | Steps | Boundary |
|---|---:|---|
| `setup` | 0 | Geometry, state, engine, and resource construction |
| `cold-first-step` | 1 | First physical step after construction |
| `cold-replay` | 1 | Second physical step, before warm-up |
| `warm-up` | 5 | Warm-up steps excluded from steady-state statistics |
| `steady-state` | 10 | Fixed measured steps |

`Capture/replay` is independent of the measurement window and is derived from
the resolved backend, never the requested backend. A successful resolved Graph
window with one or more rebuilds is `capture-inclusive`; a resolved Graph
window with no rebuild is `replay-only`. The capture-inclusive value includes
every physical step in the window and is not a capture-only timer. Resolved
Direct and Fallback rows use `not-applicable`.

The fixed workload matrix is:

| Workload | Stages | Filaments | Required requests |
|---|---:|---:|---|
| `small-single` | 1 | 10 | Direct and Graph, batch 1 |
| `medium-multi` | 2 | 32 | Direct, Graph, Persistent, explicit Fallback, mask-change Graph |
| `large-single` | 1 | 128 | Direct and Graph, batch 1 |
| `medium-multi`, batch 128 | 2 | 32 | Direct throughput |
| `medium-multi-thermal` | 2 | 32 | Graph with GPU thermal |

## Runtime Metadata

The preamble is part of every artifact:

| Field | Meaning |
|---|---|
| `Benchmark schema` | Exact schema identifier, `gpu-benchmark/v2` |
| `Commit` | Value supplied through `COILGUN_BENCHMARK_COMMIT`; `unrecorded` is invalid for retained data |
| `Compiler` | CUDA compiler and host compiler identity |
| `CUDA available` | Runtime availability and device count |
| `GPU` | Device name and compute capability |
| `CUDA driver API` | Driver API version |
| `CUDA runtime` | Runtime version |
| `Build contract` | Preset, windows, and timing relationship |
| `Benchmark repeats` | Independent constructions, warm-up steps, and measured steps |
| `Fixed workloads` | Stage, filament, and thermal shape |

## CPU Reference Rows

The CPU table is emitted before GPU rows:

```text
| Workload | Requested | Thermal | Phase | Iterations | Wall ms | Per-step ms | Finite |
```

CPU phases use the authoritative window names except that there is no CPU
warm-up row. CPU references use the fixed, unmodified stage mask. Batch speedup
multiplies the matching single-simulation CPU per-step value by the requested
GPU batch size. A CPU row with `Finite=no` invalidates its comparisons. A GPU
row with `Runtime mask change=yes` has no identically masked CPU reference in
this contract and therefore must emit `CPU/GPU speedup=n/a`.

## GPU Row Schema

The GPU table contains these columns in this exact order:

```text
| Workload | Requested backend | Requested solver | Requested precision | Requested thermal | Batch | Active ratio | Runtime mask change | Mask updates | Repeat | Phase | Measurement window | Capture/replay | Execution kind | Iterations | Resolved backend | Resolved solver | Precision | Resolved thermal | Timing relation | Setup wall ms | Cold first-step wall ms | Warm-up wall ms | Steady-state wall ms | Capture-inclusive wall ms | Replay wall ms | Fallback wall ms | Wall ms | Per-step ms | Steps/s | Simulations/s | GPU ms | Transfer ms | Mutual ms | Assembly ms | Solver ms | State update ms | Force ms | Thermal ms | Control/status ms | Sync ms | Solver status | Residual | Graph rebuild delta | Graph rebuild total | Fallback delta | Fallback total | CPU/GPU speedup | GPU executed | Finite | Fallback reason |
```

### Identity and State

| Field | Type | Meaning |
|---|---|---|
| `Workload` | string | Fixed workload name |
| `Requested backend/solver/precision/thermal` | enum | Caller request metadata, before capability resolution; never evidence of the executed path |
| `Resolved backend/solver`, `Precision`, `Resolved thermal` | enum | Actual resolved execution policy; resolved backend governs capture/replay classification |
| `Backend selection reason` | report-only enum | `ExecutionReport::backend_selection_reason`; explains Auto/explicit policy resolution but is not emitted as a v2 row column and never proves GPU execution |
| `Batch` | integer | Independent simulations in the engine |
| `Active ratio` | number | Active simulations divided by batch capacity |
| `Runtime mask change` | boolean | Whether stage masks change during measured windows |
| `Mask updates` | integer | Host mask updates in this row |
| `Repeat` | integer | Independent repeat index, starting at zero |
| `Iterations` | integer | Physical steps in this row |
| `Execution kind` | enum | `setup`, `gpu`, `fallback`, or `not-executed` |

### Window and Timing Fields

`Timing relation` is always:

```text
wall-outer;report-timers-nested-not-additive
```

The timing model is:

```text
Wall ms (outer host wall)
  aliases: one matching window wall field
  aliases: Capture-inclusive / Replay / Fallback wall when applicable
  contains or overlaps: ExecutionReport GPU / transfer / solver / thermal timers
```

The word `aliases` is normative: phase-specific wall fields repeat `Wall ms`
under a parseable name and must never be added to it. `GPU ms`, `Transfer ms`,
`Solver ms`, and `Thermal ms` are deltas of cumulative `ExecutionReport`
fields. They are diagnostic nested ranges and must not be summed with `Wall
ms` or with each other. On a CPU fallback, `GPU ms` is zero while solver or
thermal report timers can still be nonzero.

| Field group | Meaning |
|---|---|
| `Measurement window` | `setup`, `cold-first-step`, `cold-replay`, `warm-up`, or `steady-state` |
| `Capture/replay` | `capture-inclusive`, `replay-only`, or `not-applicable` |
| `Setup/Cold first-step/Warm-up/Steady-state wall ms` | Alias of `Wall ms` only for the named window; otherwise `n/a` |
| `Capture-inclusive wall ms` | Alias of a Graph window wall when capture/rebuild occurs |
| `Replay wall ms` | Alias of Graph replay wall for a replay-only window |
| `Fallback wall ms` | Alias of fallback latency; never a GPU timer |
| `Wall ms` | Complete host wall for the row |
| `Per-step ms` | `Wall ms / Iterations`; setup uses its wall value |
| `Steps/s` | `1000 / Per-step ms` for nonzero windows; otherwise zero |
| `Simulations/s` | `Steps/s * Batch` |
| `GPU/Transfer/Solver/Thermal ms` | Current `ExecutionReport` timing deltas |
| `Mutual/Assembly/State update/Force/Control/status/Sync ms` | Reserved stage columns; `n/a` until production instrumentation exists |

`n/a` means not applicable or not instrumented. Numeric zero means an
instrumented field reported zero. Consumers must preserve this distinction.

### Status and Audit Fields

| Field | Meaning |
|---|---|
| `Solver status` | `not-run` for setup or `success` after a completed step; benchmark exceptions terminate the run |
| `Residual` | `n/a` until `ExecutionReport` exposes a residual |
| `Graph rebuild delta/total` | Rebuilds in the row and since engine construction |
| `Fallback delta/total` | Fallback events in the row and since engine construction |
| `CPU/GPU speedup` | Positive finite number only for finite actual GPU execution with a matching CPU phase and identical mask behavior; otherwise `n/a` |
| `GPU executed` | `yes` only for a successful non-fallback GPU execution in this row |
| `Finite` | All available wall and report timing values are finite |
| `Fallback reason` | `none` without fallback; otherwise the retained report reason |

## Compatibility Fields

Version 2 retains the version-1 column names `Phase`, `Wall ms`, `Per-step ms`,
`Steps/s`, `GPU ms`, `Solver ms`, `Thermal ms`, `Transfer ms`, both Graph and
fallback counters, `CPU/GPU speedup`, `GPU executed`, `Finite`, and `Fallback
reason`. Requested/resolved backend, solver, and thermal columns are also
retained.

`Phase` remains for consumers of the old table. Its legacy
`first-step/capture-inclusive` value is not authoritative because Direct and
fallback first steps are not Graph captures. New consumers must use
`Measurement window` plus `Capture/replay`. Requested backend is retained as
metadata; consumers must use `Resolved backend` to identify the executed path.

## Interpretation Rules

1. A v2 artifact contains exactly 165 unique GPU rows: the documented 11
   requests, repeats 0-2, and all five measurement windows. Key rows by the
   complete request, repeat, and measurement window; reject duplicates,
   truncation, missing combinations, and unexpected combinations.
2. Never add a phase-specific wall alias to `Wall ms`.
3. Never derive total time by adding `ExecutionReport` timer columns.
4. `Execution kind=fallback` requires `GPU executed=no`, numeric `Fallback wall
   ms`, a non-`none` reason, and `CPU/GPU speedup=n/a`.
5. `gpu_executed=no` with `Finite=yes` is valid fallback or setup evidence, not
   a failed benchmark and not GPU performance.
6. `Runtime mask change=yes` requires `CPU/GPU speedup=n/a` until an
   identically masked CPU reference is emitted by the contract.
7. Capture/replay classification uses `Resolved backend`; `Requested backend`
   is metadata only. `Backend selection reason` is diagnostic report metadata
   and is not a timing or speedup input.
8. Exclude non-finite rows and CPU comparisons with `Finite=no`.
9. Summaries must retain raw values and report sample count, median, p95, and
   dispersion across independent repeats.

## Contract Check

`tests/check_gpu_benchmark_schema.cpp` rejects artifacts with missing metadata
or columns; malformed numeric, integer, boolean, or enum values; duplicate,
truncated, missing, or unexpected composite rows; capture/replay inconsistent
with the resolved backend; masked work compared to an unmasked CPU reference;
or fallback latency mistaken for GPU speedup. A complete v2 run emits exactly
165 unique GPU rows: 11 requests, 3 independent repeats, and 5 measurement
windows.
