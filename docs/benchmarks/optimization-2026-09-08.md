# Optimization Workflow Benchmark - 2026-09-08

## Environment

- Validation commit: the commit containing this report (`Validate optimization workflow`)
- Base commit embedded at configure time: `e97f165`
- Compiler: `c++ (GCC) 16.2.1 20260810`
- CPU: Intel Core Ultra 9 275HX, 24 cores / 24 online CPUs
- GPU: NVIDIA GeForce RTX 5080 Laptop GPU, driver 610.57.04
- Presets exercised: `cpu-debug`, `cpu-release`, `cuda-debug`, `cuda-release`
- Fixed seed: `20260908`
- Workload: one-stage coilgun, four voltage candidates per CPU batch, 8 Euler
  steps per evaluation, two optimizer generations, hard terminal-velocity floor
  of 0.27 m/s, Full optimization followed by Reference recheck

## CPU Release Raw Output

```text
base_commit=e97f165
preset=cpu-release
seed=20260908
setup_seconds=0.015106393
first_step_seconds=0.006852812
warmup_runs=2
steady_state_runs=5
steady_state_seconds=0.033117671
steady_state_per_batch_seconds=0.0066235342
evaluations=7
cache_hits=1
optimizer_failed_evaluations=0
failed_evaluations=1
fallback_count=1
terminal_velocity=0.2951739619
feasible=1
reference_terminal_velocity=0.2951739619
reference_recheck_error=0
gpu_batch_backend=unavailable (no concrete CUDA optimizer backend)
gpu_callback_calls=2
gpu_callback_isolated_failures=1
gpu_callback_isolated_successes=3
gpu_callback_fallback_observed=1
gpu_callback_fallback_seconds=0.006269103
gpu_callback_fallback_successes=4
first_batch_successes=4
```

## Preset Comparison

| Preset | Setup (s) | First CPU batch (s) | Steady CPU batch (s) | Callback fallback (s) |
|---|---:|---:|---:|---:|
| `cpu-debug` | 0.014944921 | 0.006941699 | 0.0063305108 | 0.004846527 |
| `cpu-release` | 0.015106393 | 0.006852812 | 0.0066235342 | 0.006269103 |
| `cuda-debug` | 0.015040419 | 0.006953425 | 0.0062473156 | 0.004858903 |
| `cuda-release` | 0.015096251 | 0.006833750 | 0.0062668338 | 0.006077716 |

The CUDA-configured runs intentionally report the same CPU workload. The
adapter exposes a GPU batch callback and CPU fallback boundary, but this branch
does not provide a concrete CUDA optimizer evaluator to inject. Therefore no
GPU throughput or CPU/GPU speedup is claimed. The measured callback tests show
that one failed callback result does not poison three successful candidates,
and a thrown callback falls back to four successful CPU evaluations.

## Decision

**Accept as the initial optimization benchmark baseline.** There is no earlier
equivalent benchmark from which to establish a performance regression. The
measured CPU batch times are stable enough for recorded evidence, all optimized
results are feasible, and Full-to-Reference error is exactly zero for this
fixed workload. Regression classification is deferred until a subsequent run
can compare against this baseline; no numerical tolerance was loosened.
