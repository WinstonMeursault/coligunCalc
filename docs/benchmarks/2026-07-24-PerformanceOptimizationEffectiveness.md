# Performance Optimization Effectiveness

## Scope and verdict

This report evaluates the B0-B7 execution-speed work at branch baseline
`873c72c` plus the B8 report artifacts. The primary metric is execution time;
compile and CTest time are secondary engineering feedback.

Overall verdict: **mixed positive**. The CPU simulation path has a clear,
repeatable short-window gain after the B7 final review wired the fused mutual
and gradient evaluator into both simulation hot loops. CUDA execution remains
healthy and stable, but the available measurements do not prove a broad new
CUDA speedup. Graph runtime-mask changes still cost materially more than normal
Graph replay, and inactive-batch and long-run CPU targets remain workload
specific or unproven.

## CPU fixed-workload result

Normal Release `bench_cpu_sim 3`, median steady-state milliseconds per step:

| Workload | Before | After | Change |
|---|---:|---:|---:|
| single-16 | 1.677208 | 1.246797 | -25.7% |
| multi-32 | 2.349867 | 1.735920 | -26.1% |
| multi-128 | 5.582488 | 3.363284 | -39.8% |
| single-16-thermal | 1.611413 | 1.342622 | -16.7% |
| multi-32-thermal | 2.274128 | 1.683605 | -26.0% |
| multi-128-thermal | 5.423722 | 3.216753 | -40.7% |

The protocol is a four-step warmed window, so these values do not establish
the B0 long physical-duration CPU target. They do establish that the fused
pair integration removes a measured duplicate quadrature traversal in the
representative simulation step.

## CUDA fixed-workload result

CUDA Release `bench_gpu_engine`, median steady-state milliseconds per step:

| Workload | Before | After | Interpretation |
|---|---:|---:|---|
| Direct/Eigen/B=32 | 8.759 | 8.692 | Stable; -0.8% is within run variance |
| Graph/Batched/B=32 | 8.680 | 8.679 | Stable |
| Graph/Batched/B=32, mask change | 13.699 | 13.695 | Stable, but mask update overhead remains |
| Thermal Graph/Batched/B=32 | 16.547 | 8.707 | Inconclusive; same-host runs varied about 8.69-16.55 |

All 90 current GPU phase rows were finite. Graph mask changes retained one
cumulative rebuild and did not recapture. Persistent and explicit fallback
rows correctly reported `gpu_executed=no` and were excluded from speedup.

## Target disposition

| Target | Result |
|---|---|
| CPU local fused mutual/gradient speed | Met for all six aligned short-window workloads |
| CPU representative long run | Not proven; retain and measure with aligned physical workload |
| GPU transfer reduction | Local prior evidence retained; no universal percentage claim |
| Graph topology reuse | Contract and rebuild target pass |
| Graph runtime-mask overhead | Not optimized; B9 candidate |
| Non-Graph GPU regression <=5% | Pass for measured B32 |
| Inactive-batch throughput target | Not proven; focused sample was mixed |
| Fallback correctness | Pass; no fallback speedup claim |
| Numeric behavior | CPU `21/21`, CUDA `37/37`, all benchmark rows finite |

## Environment and verification

- CPU: Linux, 24 logical CPUs, GNU C++ 16.1.1, Ninja, CMake Release,
  `-march=native`, OpenMP 5.2.
- GPU: NVIDIA GeForce RTX 5080 Laptop GPU, compute 12.0, CUDA driver/runtime
  13.3.
- Full CPU Release CTest: `21/21` passed.
- Full CUDA Release CTest: `37/37` passed.
- Focused CPU tests after the fused-pair integration: `5/5` passed.
- No compilation-time improvement is included in the verdict.

This document contains the retained benchmark comparisons, protocol summary,
and task attribution. Temporary per-task benchmark artifacts were removed after
the results were consolidated here.

## Next review baseline

B9 starts from the optimized code and these reports. It is read-only and must
re-examine the whole project through CPU, CUDA, simulation, architecture/debt,
test/feedback, and cross-check lenses. A candidate is counted as new only when
it is not already covered by the historical findings or the deferred items
listed above.
