# Optimization Workflow Benchmark - 2026-09-08

## Environment

- Validation commit: the commit containing this report (`Validate optimization workflow`)
- Source revision embedded at configure time: `677ce3f`
- Worktree state embedded at configure time: `clean`
- Compiler: `c++ (GCC) 16.2.1 20260810`
- CPU: Intel Core Ultra 9 275HX, 24 cores / 24 online CPUs
- GPU: NVIDIA GeForce RTX 5080 Laptop GPU, driver 610.57.04
- Presets exercised: `cpu-debug`, `cpu-release`, `cuda-debug`, `cuda-release`
- Fixed seed: `20260908`
- Workload: one-stage coilgun, four voltage candidates per CPU batch, 8 Euler
  steps per evaluation, two optimizer generations, initial armature position
  0.015 m, hard terminal-velocity floor of 0.0095 m/s, Full optimization
  followed by Reference recheck

## CPU Release Raw Output

```text
source_revision=677ce3f
worktree_state=clean
preset=cpu-release
seed=20260908
setup_seconds=0.015040301
first_batch_seconds=0.000373103
warmup_runs=2
steady_state_runs=5
steady_state_seconds=0.001433422
steady_state_per_batch_seconds=0.0002866844
evaluations=7
cache_hits=1
optimizer_failed_evaluations=0
failed_evaluations=1
fallback_count=1
terminal_velocity=0.009764939958
feasible=1
reference_terminal_velocity=0.009764969509
reference_recheck_error=2.955109129e-08
gpu_batch_backend=unavailable (no concrete CUDA optimizer backend)
gpu_callback_calls=2
gpu_callback_isolated_failures=1
gpu_callback_isolated_successes=3
gpu_callback_fallback_observed=1
gpu_callback_fallback_seconds=0.000313326
gpu_callback_fallback_successes=4
first_batch_successes=4
```

## Preset Comparison

| Preset | Setup (s) | First CPU batch (s) | Steady CPU batch (s) | Callback fallback (s) |
|---|---:|---:|---:|---:|
| `cpu-debug` | 0.014895153 | 0.000360979 | 0.0002855976 | 0.000306574 |
| `cpu-release` | 0.015040301 | 0.000373103 | 0.0002866844 | 0.000313326 |
| `cuda-debug` | 0.014927227 | 0.000386981 | 0.0002815340 | 0.000309433 |
| `cuda-release` | 0.015056662 | 0.000361383 | 0.0002856108 | 0.000311673 |

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
results are feasible, and the nonzero Full-to-Reference error is
`2.955109129e-08`, within the asserted `5e-8 + 1e-6 * |reference|` bound.
Regression classification is deferred until a subsequent run
can compare against this baseline; no numerical tolerance was loosened.
