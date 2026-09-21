# GPU optimization workflow benchmark — 2026-09-17

## Decision

PASS for the B-T4 throughput criterion on the measured workload. GPU
per-candidate throughput exceeded CPU at batch sizes 32 and 128, and all
execution, ordering, numerical, fallback, and optimizer Reference gates passed.
Batch size 1 is reported separately: its latency is slower on GPU because the
fixed CUDA launch/transfer overhead is not amortized.

The measurement was made from clean detached revision `6f6290a` (`git status
--porcelain` was empty before the run). The final report/ledger commit is a
documentation-only descendant; the benchmark source and CMake files are
unchanged between the measured revision and that descendant. This is the
unavoidable self-reference boundary: a Git commit cannot contain its own SHA.
The containing final commit is identified immutably by subject `Benchmark GPU
optimization workflow`, parent `b294888144733f620d0088d77c4a43030d483e83`,
and the final `git rev-parse HEAD` output supplied with this evidence. The
measured and evidence-stage SHAs, source-preservation diff, and manifest remain
durable in this report; the final SHA is deliberately not invented inside it.

## Reproduction

The exact clean-worktree commands were:

```sh
git worktree add --detach /tmp/coligun-bt4-final-measure 6f6290ad5e2da7dadf0fc40b4a7c0306733e6c88
cmake -S /tmp/coligun-bt4-final-measure -B /tmp/coligun-bt4-final-measure/build/cuda-release \
  -G Ninja -DCOILGUN_ENABLE_CUDA=ON -DCOILGUN_BUILD_TESTS=ON \
  -DCOILGUN_BUILD_GENERATOR_TESTS=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS=-march=native -DCMAKE_CUDA_ARCHITECTURES=native
cmake --build /tmp/coligun-bt4-final-measure/build/cuda-release \
  --target bench_gpu_optimization_workflow -j2
/tmp/coligun-bt4-final-measure/build/cuda-release/tests/bench_gpu_optimization_workflow \
  /tmp/optimization-gpu-workflow-2026-09-17-final.json
python3 scripts/validate_optimization_gpu_workflow.py \
  /tmp/optimization-gpu-workflow-2026-09-17-final.json
```

The machine-readable output is schema version 1. It contains source state,
hardware, workload, raw samples, medians, throughput, speedup, backend proof,
numerical deltas, optimizer counters, and the decision booleans.

## Environment and workload

| Field | Measured value |
|---|---|
| GPU | NVIDIA GeForce RTX 5080 Laptop GPU |
| Kernel driver | 615.71.09 (`nvidia-smi`) |
| Compute capability | 12.0 |
| CUDA toolkit/runtime | 13.4.59 / 13.4 |
| CUDA driver API | 13.4 |
| Compiler | GCC 16.2.1 |
| CMake | 4.4.3 |
| Revision | 6f6290a; clean detached worktree |
| Geometry | one fixed 150-turn coil, 2×2 armature filament discretization |
| Precision/path | CPU Full and production `CudaBatchEvaluator` Full; Euler |
| `dt`, termination | `1e-6`; `max_steps=64`, velocity/bound checks disabled |
| Seed | 20260917 |
| Warm-ups/repetitions | 2 / 5 per path and batch size |
| Batch sizes | 1, 8, 32, 128 |
| Excluded metric | `PeakCurrent` (known unsupported CUDA production metric) |

The deterministic candidate generator produces 128 distinct voltage values by
the monotone stride `455 + index*(95/128)` (455.0 through 549.2578125 V), all
inside the schema range 450–550 V. Each CPU and GPU iteration uses the same
candidate vector and fixed geometry. GPU requests use explicit `Direct`; the
resolved backend and `gpu_executed` flag are checked on every measured
iteration.

## Raw timing samples and summaries

Times are host wall milliseconds. Throughput is candidates/second; speedup is
CPU median divided by GPU median (equivalently the per-candidate ratio).

| Batch | CPU raw samples (ms) | GPU raw samples (ms) | CPU median | GPU median | CPU cand/s | GPU cand/s | Speedup |
|---:|---|---|---:|---:|---:|---:|---:|
| 1 | 34.665722, 34.742488, 35.140955, 35.293916, 35.239124 | 51.792762, 58.700427, 51.545049, 58.093642, 51.761262 | 35.140955 | 51.792762 | 28.4568 | 19.3077 | 0.6785× |
| 8 | 277.664509, 277.577191, 277.602420, 279.024654, 277.852371 | 51.780566, 52.909497, 60.230878, 52.000980, 52.230258 | 277.664509 | 52.230258 | 28.8117 | 153.1679 | 5.3162× |
| 32 | 1110.211810, 1111.648846, 1112.064753, 1110.851615, 1110.843584 | 121.694211, 122.334645, 122.301710, 119.654975, 113.376586 | 1110.851615 | 121.694211 | 28.8067 | 262.9542 | 9.1282× |
| 128 | 4440.193372, 4439.247427, 4442.526380, 4442.239252, 4443.064817 | 354.990762, 360.498146, 359.258753, 362.067273, 359.969995 | 4442.239252 | 359.969995 | 28.8143 | 355.5852 | 12.3406× |

Representative rows (first, middle, last) were checked in every batch. Maximum
absolute/relative supported-metric deltas were:

| Batch | Max absolute delta | Max relative delta | Order | Numerical gate |
---:|---:|---:|---|---|
| 1 | 5.684341886080802e-14 | 1.8987268313961294e-14 | pass | pass |
| 8 | 5.684341886080802e-14 | 1.8400033211467645e-14 | pass | pass |
| 32 | 5.684341886080802e-14 | 1.8400033211467645e-14 | pass | pass |
| 128 | 5.684341886080802e-14 | 1.8400033211467645e-14 | pass | pass |

The checked production metrics were muzzle velocity, peak voltage, maximum
temperature metadata, and every constraint row. Constraint count/order,
ID/kind/relation, bounds, values, violation, normalized violation,
`satisfied`, and `priority` were compared with the existing Full tolerance.
`PeakCurrent` was not used for fitness or pass/fail.

Resolved backend was `direct`; batch 1 resolved solver `eigen`, and batches
8/32/128 resolved `cusolver`. Every measured iteration reported
`gpu_executed=true`, `precision=full`, zero fallback events, and zero failed
batches.

## Same-seed end-to-end optimization

Both runs used the same one-variable schema (`voltage`, 450–550 V), geometry,
Full level, seed 20260917, population 8, elite count 1, crossover 0.8,
mutation 0.2, and four generations. CPU used the ordinary CPU problem
evaluator; GPU used the production `CudaBatchEvaluator`.

| Field | CPU | GPU |
|---|---:|---:|
| Termination | maximum generations | maximum generations |
| Evaluations / successes / failures | 32 / 32 / 0 | 32 / 32 / 0 |
| Best voltage | 550 V | 550 V |
| Best Full objective | 0.012948805601044723 | 0.012948805601044600 |
| Feasible | true | true |
| GPU batches / executed / successful | 0 / 0 / 0 | 4 / 32 / 32 |
| GPU failed batches / fallbacks / CPU repairs | 0 / 0 / 0 | 0 / 0 / 0 |
| Reference objective | 0.012948805601044723 | 0.012948805601044723 |
| Reference error / tolerance / valid | 0 / 6.294880560104472e-08 / true | 1.231653667943533e-16 / 6.294880560104472e-08 / true |

The GPU optimizer's run-local counters were
`gpu_transfer_seconds=0.02674878199999997`,
`gpu_kernel_seconds=0.20930406699999998`, and
`gpu_elapsed_seconds=0.215024511`. The CPU run reported all GPU counters as
zero, as expected.

## Validator and verification evidence

TDD validator cycle:

1. Added eleven schema tests and observed RED failures for stale decision flags
   and zero timing samples.
2. Implemented the validator; `python3 -m unittest
   scripts/test_validate_optimization_gpu_workflow.py` passed 11/11.
3. The captured clean-run JSON passed `validate_optimization_gpu_workflow.py`.

The validator rejects missing keys, non-finite, negative, or zero timing samples,
false GPU execution, fallback backend, missing backend/solver/precision fields,
missing optimizer statistics/best/reference fields, non-finite or negative
timing, order/numerical failures, inconsistent counts, and derived median,
throughput, speedup, Reference-tolerance, or derived decision mismatches. It
intentionally allows a measured decision of `false` when timing shows GPU is
slower, so a genuine throughput failure is reportable rather than hidden.

Focused build/run evidence:

```text
cmake --build ... --target bench_gpu_optimization_workflow -j2
  completed successfully (CUDA release target)
python3 -m unittest scripts/test_validate_optimization_gpu_workflow.py
  Ran 11 tests ... OK
python3 scripts/validate_optimization_gpu_workflow.py /tmp/optimization-gpu-workflow-2026-09-17-final.json
  valid B-T4 schema v1
ctest --preset cuda-release -R 'test_cuda_batch_evaluator|test_gpu_optimization_baseline' --output-on-failure
  100% tests passed out of 2
ctest --preset cpu-release -R 'test_optimization_integration|test_coilgun_optimization|test_multi_stage_sim' --output-on-failure
  100% tests passed out of 3
python3 scripts/check_api_optimization_docs.py
  API optimization bilingual check: PASS

The benchmark-affecting manifest (SHA-256, measured clean revision) is:

```text
b8bf6c42c268dbea41c1e64a0c6b8f46e36b404b1895dd1524d68c4470d7737f  tests/bench_gpu_optimization_workflow.cpp
0bffc419373308459b95ed5067228daf7ac99c78f494b4d7ee45d854941b82ed  tests/CMakeLists.txt
83745ae73e6f23fcbe0461de6e56dd344b6b9b4b18c7d1798e3d2aa2f250fb76  scripts/validate_optimization_gpu_workflow.py
9099c98bad825ade76c66534d5bc8c524511aa5c72462c5ddd10d35953b17ff3  scripts/test_validate_optimization_gpu_workflow.py
```

At the evidence stage, before the final amend, the exact source-preservation
check was (recorded evidence-stage commit
`79d7577f52728e63e90b9cf67e69895a60bdba64`):

```sh
git diff --exit-code 6f6290ad5e2da7dadf0fc40b4a7c0306733e6c88 79d7577f52728e63e90b9cf67e69895a60bdba64 -- \
  tests/bench_gpu_optimization_workflow.cpp tests/CMakeLists.txt \
  scripts/validate_optimization_gpu_workflow.py \
  scripts/test_validate_optimization_gpu_workflow.py
```

It returned exit code 0 with no output. The final amend cannot contain its own
SHA; therefore the measured SHA, evidence-stage SHA, zero-output diff, and
post-amend manifest are recorded separately. The post-amend manifest was
reprinted and matched the four lines above exactly.

The focused real-GPU evaluator and CPU optimization/integration tests passed in
the final verification run; B-T4 does not claim unrelated whole-branch
verification. Existing CUDA compiler warnings (Eigen relaxed
constexpr and the volatile generation increment) were unchanged baseline
warnings and did not affect the successful benchmark executable.

## Concerns and scope

- Batch 1 is slower on GPU; it is reported separately and is not used as the
  throughput success gate.
- The production CUDA `PeakCurrent` discrepancy remains excluded as required;
  no physics or tolerance was changed.
- This benchmark validates the fixed-geometry, non-thermal, Euler Full path.
  CUDA RK4, thermal optimization metrics, and arbitrary-geometry candidates
  remain outside the production evaluator contract.
- The benchmark is intentionally versioned evidence, not an ordinary CTest
  wall-clock threshold.
