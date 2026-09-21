# B-T4 — Benchmark GPU optimization workflow

## Status

`DONE`

B-T4 adds the CUDA-only production-path benchmark executable, schema validator,
validator tests, and clean real-device evidence. The fresh measured revision is
`6f6290ad5e2da7dadf0fc40b4a7c0306733e6c88`; it was built and run from a clean
detached worktree on the RTX 5080 Laptop GPU. The final report/ledger commit is
a documentation-only descendant, with the benchmark source/CMake/validator
files unchanged between the measured revision and final descendant.

## TDD evidence

The validator tests were written before the validator implementation. The first
run was RED with a `FileNotFoundError` for the missing validator module. After
the minimal implementation, the exact test command passed:

```text
python3 -m unittest scripts/test_validate_optimization_gpu_workflow.py
Ran 11 tests in 0.002s
OK
```

The tests cover acceptance, missing keys, non-finite/negative/zero timing, false GPU
execution/fallback backend, missing backend/solver/precision and optimizer
fields, derived timing corruption, Reference error/validity corruption,
derived decision corruption, and optimizer count consistency. The validator checks schema version 1, required
keys, finite non-negative timing/summary values, recomputed medians/throughput/
speedup, non-fallback GPU execution, order/numerical flags, full optimizer
statistics, Reference tolerance derivation, and count relationships.

## Benchmark implementation

`tests/bench_gpu_optimization_workflow.cpp` uses one deterministic fixed
geometry (150-turn coil, 2×2 armature), Euler, `OptimizationLevel::Full`,
`dt=1e-6`, 64 max steps, fixed seed `20260917`, and candidate batches
1/8/32/128. It executes two warmups followed by five measured CPU and GPU
iterations per size in one process. The 128 candidate voltages are guaranteed
distinct by the monotone stride `455 + index*(95/128)`. GPU calls use the production
`CudaBatchEvaluator`, explicit Direct backend, and require `gpu_executed`, Full
precision, non-Fallback backend, and zero fallback/failure counters on every
measured iteration. Every returned row is compared by index; representative
first/middle/last rows compare CPU Full and GPU Full on muzzle velocity, peak
voltage, maximum temperature metadata, and all constraint fields. PeakCurrent
is explicitly excluded.

The executable also runs same-seed CPU/GPU four-generation, population-eight
single-objective optimization. It records termination, best candidate/objective,
feasibility, run-local evaluator counters, and CPU Reference rechecks for both
best candidates.

## Clean measurement and results

Exact measurement commands and the complete raw output are versioned in
`docs/benchmarks/optimization-gpu-workflow-2026-09-17.md`. The clean detached
run configured CUDA Release, built `bench_gpu_optimization_workflow`, emitted
schema-versioned JSON, and passed the maintained validator.

Measured median throughput:

| Batch | CPU cand/s | GPU cand/s | Speedup |
|---:|---:|---:|---:|
| 1 | 28.4568 | 19.3077 | 0.6785× |
| 8 | 28.8117 | 153.1679 | 5.3162× |
| 32 | 28.8067 | 262.9542 | 9.1282× |
| 128 | 28.8143 | 355.5852 | 12.3406× |

The B-T4 success criterion passes at batch 32 and 128. Batch 1 latency is
reported separately and is slower on GPU due to launch/transfer overhead.

Every measured batch had `gpu_executed=true`, resolved backend `direct`, zero
fallback events, zero failed batches, and Full numerical/order/constraint checks
passing. Maximum representative supported-metric relative deltas were below
`1.9e-14`. CPU and GPU optimizers both terminated at maximum generations, were
feasible, selected 550 V, and passed Reference recheck using
`5e-8 + 1e-6*abs(reference)` (tolerance `6.294880560104472e-08`); GPU
statistics recorded 32/32 executed/successful evaluations, four batches, and
zero failures, fallbacks, or CPU repairs.

Hardware evidence: NVIDIA GeForce RTX 5080 Laptop GPU, kernel driver
615.71.09, compute capability 12.0, CUDA toolkit/runtime 13.4.59/13.4,
CUDA driver API 13.4, GCC 16.2.1, and CMake 4.4.3.

## Review and concerns

Focused real-device and CPU verification was actually run after the fresh
measurement: CUDA evaluator/baseline tests passed 2/2, CPU optimization,
coilgun, and multi-stage numerical/integration tests passed 3/3, the validator
passed 11/11, and the bilingual API documentation check passed. `git diff --check`
passed on the final tree. The benchmark
does not modify physics, CUDA kernels, tolerances, or the production evaluator.
The known unsupported CUDA PeakCurrent metric remains excluded. Scope remains
fixed-geometry, non-thermal Euler Full; CUDA RK4, thermal metrics, and
arbitrary-geometry candidates are not covered. Existing CUDA compiler warnings
are baseline warnings. This report does not claim final whole-branch
verification.

## Provenance evidence

The measured clean revision is `6f6290ad5e2da7dadf0fc40b4a7c0306733e6c88`.
Before the final amend, the documentation/ledger evidence stage was created as
commit `79d7577f52728e63e90b9cf67e69895a60bdba64`; the exact benchmark-source preservation command
returned exit code 0 with no output:

```sh
git diff --exit-code 6f6290ad5e2da7dadf0fc40b4a7c0306733e6c88 79d7577f52728e63e90b9cf67e69895a60bdba64 -- \
  tests/bench_gpu_optimization_workflow.cpp tests/CMakeLists.txt \
  scripts/validate_optimization_gpu_workflow.py \
  scripts/test_validate_optimization_gpu_workflow.py
```

The measured and post-amend SHA-256 manifests for those four files were
identical:

```text
b8bf6c42c268dbea41c1e64a0c6b8f46e36b404b1895dd1524d68c4470d7737f  tests/bench_gpu_optimization_workflow.cpp
0bffc419373308459b95ed5067228daf7ac99c78f494b4d7ee45d854941b82ed  tests/CMakeLists.txt
83745ae73e6f23fcbe0461de6e56dd344b6b9b4b18c7d1798e3d2aa2f250fb76  scripts/validate_optimization_gpu_workflow.py
9099c98bad825ade76c66534d5bc8c524511aa5c72462c5ddd10d35953b17ff3  scripts/test_validate_optimization_gpu_workflow.py
```

The final commit cannot contain its own SHA. Its immutable identity is the
subject `Benchmark GPU optimization workflow`, parent
`b294888144733f620d0088d77c4a43030d483e83`, and the final `git rev-parse HEAD`
output supplied with this evidence. The measured SHA, evidence-stage SHA,
zero-output diff, and manifest are durable here; no self SHA is invented.
