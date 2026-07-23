# Integrated Pipeline Verification Record

**Date:** 2026-07-23
**Configuration:** `cuda-debug`, NVIDIA CUDA device available

## Short Verification

The following commands were run after the resident Graph and device-control changes:

```sh
cmake --build --preset cuda-debug --target test_gpu_vs_cpu_single \
  test_gpu_vs_cpu_multi test_gpu_sim_batch test_gpu_graph_pipeline \
  test_gpu_solver test_gpu_resource_contracts test_gpu_assembly \
  test_gpu_resident_pipeline test_gpu_engine_physics
ctest --preset cuda-debug -R \
  "test_gpu_(solver|resource_contracts|assembly|resident_pipeline|engine_physics|graph_pipeline)$" \
  --output-on-failure
```

Result: build succeeded; 6/6 targeted CTest suites passed.

Additional bounded doctest cases passed:

- Multi-stage Graph execution versus explicit CPU fallback.
- Position and time-delay trigger activation and pre-step trigger boundaries.
- Residual-current and quiet-stage completion lifecycle.
- Heterogeneous SimBatch rows, fixed-capacity non-compaction, and per-row trigger/excitation state.

These cases require `gpu_executed=true` when CUDA is available and separately exercise explicit fallback paths.

## Deferred Evidence

- The complete `test_gpu_sim_batch` executable exceeded the agreed short-test budget and was deferred.
- `compute-sanitizer` memcheck/racecheck was not run in this short verification window.
- Full CPU/CUDA Debug and Release matrix, benchmark executable, and throughput measurements remain to be run.

Fallback timing rows must not be used as GPU speedup evidence. The benchmark must report wall time, solver, thermal, transfer, control/synchronization, Graph rebuild count, steps/s, and simulations/s separately for actual GPU and fallback rows.

## Current Gate Decisions

- Structured solver: **Not approved**, see `2026-07-21-structured-solver-feasibility-review.md`.
- GPU RK4: **Not approved**, see `2026-07-21-gpu-rk4-feasibility-review.md`.
