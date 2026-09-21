# W5 Final Integration Validation — 2026-09-17

## Decision

**PASS for the W5 verification work performed here.** All required fresh
configure/build/full-CTest cycles, real-device CUDA optimization checks,
benchmark/schema gates, documentation checks, install checks, and external
consumer builds passed. The roadmap's final independent whole-branch review is
now complete: the final read-only Luna High review approved
`7e43c4e..8bc0488` with 0 Critical, 0 Important, and 0 Minor findings.

No physics, CUDA kernel, production tolerance, or optimizer search-behavior
changes were made for the validation or the post-review identity/validator
fix wave. The fix wave changes cache identity contracts, callback
synchronization, benchmark provenance validation, documentation, tests, and
ledgers.

## Revision and worktree evidence

- Branch: `feature/optimization-next-phase`
- Starting revision: `084eb66` (`Benchmark GPU optimization workflow`)
- Starting worktree state: clean
- Final validation worktree was kept separate from the primary and task
  worktrees. At inspection time the primary worktree retained its four dirty
  optimization files, and `fix/opt-t6-convergence` retained its dirty
  `tests/test_optimization_single.cpp`; neither was edited or cleaned.
- Other worktrees were not deleted or modified.

The report is being committed as the exact subject `Validate next optimization phase`.

## Fresh required preset verification

Each command was run from this worktree, with output captured in
`/tmp/coligun-w5-{preset}.log`:

```text
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug --output-on-failure
  exit 0; 35/35 passed; 0 skipped; 0 failed

cmake --preset cpu-release
cmake --build --preset cpu-release
ctest --preset cpu-release --output-on-failure
  exit 0; 35/35 passed; 0 skipped; 0 failed

cmake --preset cuda-debug
cmake --build --preset cuda-debug
ctest --preset cuda-debug --output-on-failure
  exit 0; 53/53 passed; 0 skipped; 0 failed

cmake --preset cuda-release
cmake --build --preset cuda-release
ctest --preset cuda-release --output-on-failure
  exit 0; 53/53 passed; 0 skipped; 0 failed
```

The CUDA suites exercised 19 `gpu`-labelled tests in each configuration and
were serialized by the existing CTest resource lock. Totals are discovered
totals from this run, not frozen documentation claims.

## Independent real-GPU optimization evidence

Commands, run independently after the full CUDA Release suite:

```sh
./build/cuda-release/tests/test_gpu_optimization_baseline
./build/cuda-release/tests/test_cuda_batch_evaluator
```

Results:

```text
test_gpu_optimization_baseline: 1 test case, 154/154 assertions passed
test_cuda_batch_evaluator:      21 test cases, 229/229 assertions passed
```

The baseline output proved `gpu_executed=true`, `backend=direct`,
`precision=full` for real `SimBatch` execution at batch sizes 1, 8, 32, and
128, including thermal rows. The evaluator suite directly covered stable row
ordering, local invalid-row isolation, structured protocol failure, strict and
opt-in fallback policies, run-local metrics, and CPU Reference terminal
velocity/constraint parity. No unintended fallback was reported.

Device/toolchain evidence captured during this run:

```text
GPU: NVIDIA GeForce RTX 5080 Laptop GPU
Driver: 615.71.09; compute capability: 12.0
CUDA toolkit/compiler: 13.4.59; runtime/driver API: 13.4
GCC: 16.2.1; CMake: 4.4.3
```

The known B-T1 PeakCurrent discrepancy remained visible in the baseline
classification output (approximately 6.4–6.7% below CPU); it was not hidden,
retuned, or used as a production optimization fitness metric.

## Fresh CUDA Release benchmark and validator

```sh
build/cuda-release/tests/bench_gpu_optimization_workflow \
  /tmp/optimization-gpu-workflow-w5-final.json
python3 scripts/validate_optimization_gpu_workflow.py \
  /tmp/optimization-gpu-workflow-w5-final.json
```

The validator returned `valid B-T4 schema v1`. Fresh median results were:

| Batch | CPU cand/s | GPU cand/s | Speedup | GPU/backend/precision | Order | Numerical | Fallback events |
|---:|---:|---:|---:|---|---|---|---:|
| 1 | 28.6948 | 16.7698 | 0.5844x | true / direct / full | true | true | 0 |
| 8 | 28.8068 | 153.1198 | 5.3156x | true / direct / full | true | true | 0 |
| 32 | 26.5631 | 264.8050 | 9.9688x | true / direct / full | true | true | 0 |
| 128 | 28.8138 | 354.4370 | 12.2854x | true / direct / full | true | true | 0 |

The machine-readable decision was:

```json
{"batch32_gpu_faster":true,"batch128_gpu_faster":true,
 "success_criterion":"gpu per-candidate throughput exceeds CPU at batch 32 and 128"}
```

The GPU optimizer section recorded seed `20260917`, 32/32 successful
evaluations, four GPU batches, 32/32 GPU-executed/successful evaluations,
zero failed batches, zero GPU fallbacks, zero CPU repairs, and a valid CPU
Reference recheck (`1.231653667943533e-16` error versus
`6.294880560104472e-08` tolerance). Batch 1 is reported separately and is
slower on GPU; it is not used for the throughput decision.

## Documentation, validator, and diff checks

```text
python3 scripts/check_api_optimization_docs.py --self-test
  exit 0; negative self-test PASS; bilingual check PASS;
  semantic_contracts=6; common_tokens=30; actual_presets configure=6 build=6 test=12

python3 -m unittest scripts/test_validate_optimization_gpu_workflow.py
  exit 0; Ran 11 tests; OK

git diff --check
  exit 0; no output
```

## Installed package and external consumers

Separate temporary prefixes were used:

```text
CPU:  /tmp/coligun-w5-cpu-install.11O9Dv
CUDA: /tmp/coligun-w5-cuda-install.41C05K
```

```text
cmake --install build/cpu-release --prefix <CPU prefix>
python3 scripts/verify_install_headers.py <CPU prefix>
  exit 0; install header split: PASS (CPU-only)

cmake --install build/cuda-release --prefix <CUDA prefix>
python3 scripts/verify_install_headers.py <CUDA prefix> --cuda
  exit 0; install header split: PASS (CUDA)
```

The external projects were under `/tmp/coligun-w5-cpu-consumer` and
`/tmp/coligun-w5-cuda-consumer`, outside this repository. Both used only
`find_package(coilgun CONFIG REQUIRED)` and the installed prefix.

CPU consumer source:

```cpp
#include <coilgun/coilgun.hpp>
using namespace coilgun::optimization;
VariableSchema schema({VariableSpec::continuous("voltage", 100.0, 500.0)});
ProblemSpec spec(std::move(schema),
                 {ObjectiveDefinition{"muzzle_velocity", true, 1.0}});
CandidateVariables candidate{{250.0}};
const auto repaired = spec.repair(candidate);
```

It linked only `coilgun::coilgun`, built with exit 0, and ran:
`voltage=250`.

CUDA consumer source:

```cpp
#include <coilgun/coilgun_cuda.hpp>
using namespace coilgun::optimization;
static_assert(std::is_class_v<CudaBatchEvaluator>);
const CudaFallbackOptions options{CudaFallbackPolicy::Strict};
```

It linked only `coilgun::coilgun_cuda`, built with exit 0, and ran:
`cuda_batch_evaluator=strict`. Generated Ninja files contained the installed
prefix include path and no repository `include/` or validation build path.

## Requirement audit — roadmap §9

| Requirement | Authoritative evidence | Decision |
|---|---|---|
| §9.1 single objective routes to GA; multi objective to NSGA-II | `tests/test_optimization_routing.cpp` cases “Auto routing is equivalent…” and “Auto routing uses NSGA-II…”; full CPU/CUDA totals | PASS |
| §9.1 default result is Pareto-only; representative is explicit | `tests/test_optimization_selectors.cpp` explicit-selection/front-preservation cases; `test_optimization_integration` | PASS |
| §9.1 CPU behavior and fixed-seed reproducibility | `tests/test_optimization_single.cpp` deterministic case; `tests/test_optimization_integration.cpp`; CPU Debug/Release 35/35 | PASS |
| §9.1 shared-cache evaluator isolation | `tests/test_optimization_evaluator.cpp` namespace/version isolation cases; A-T1 report | PASS |
| §9.1 invalid thermal and NSGA stagnation fail before run | `tests/test_coilgun_optimization.cpp` thermal configuration cases; `tests/test_optimization_problem_spec.cpp` and routing stagnation cases | PASS |
| §9.1 one GPU row failure does not poison successes | `tests/test_cuda_batch_evaluator.cpp` local invalid/row-failure/fallback cases; independent 229/229 run | PASS |
| §9.1 GPU output order equals input order | production evaluator stable-order case; fresh benchmark `order_ok=true` for 1/8/32/128 | PASS |
| §9.1 GPU optimization result CPU Reference recheck | fresh benchmark GPU `reference_valid=true`, error `1.23e-16`; baseline Reference case | PASS |
| §9.1 GPU proof requires executed + non-Fallback backend | baseline and evaluator assertions; benchmark all `gpu_executed=true`, backend `direct` | PASS |
| §9.2 four configure/build/full CTest presets | fresh commands above: 35/35, 35/35, 53/53, 53/53 | PASS |
| §9.2 independent real GPU integration test | independent baseline/evaluator commands: 154/154 and 229/229 | PASS |
| §9.2 device/driver/toolkit/precision/backend recorded | device block, baseline output, and raw benchmark JSON | PASS |
| §9.2 installed-package umbrella consumer | two out-of-tree consumers, header split checks, exit-0 builds/runs | PASS |
| §9.2 bilingual API equivalence | `check_api_optimization_docs.py --self-test`: PASS | PASS |
| §9.2 benchmark schema validation | fresh JSON plus validator: `valid B-T4 schema v1` | PASS |
| §9.2 final independent whole-branch review | Luna High review package `review-7e43c4e..8bc0488.diff`; final verdict Approved, 0 Critical/Important/Minor | PASS |
| §9.3 batch 32/128 GPU per-candidate throughput > CPU | fresh JSON decision booleans both true | PASS |
| §9.3 same round/config/precision for speedups | benchmark source uses one process, same candidate vectors, Full path; raw JSON records workload | PASS |
| §9.3 batch 1 separately reported | fresh JSON has a separate batch-1 row and decision excludes it | PASS |
| §9.3 no tolerance widening for performance | benchmark uses existing Full comparison gates; validator recomputes deltas; no behavior edits | PASS |
| §9.3 FP32 experiment scope/default reliability | no FP32 performance claim; benchmark is Full; existing precision tests passed | PASS / N/A |

## Requirement audit — roadmap §5 A/B deliverables and interface gate

| Deliverable | Evidence | Decision |
|---|---|---|
| A-T1 cache identity/isolation | `A-T1-report.md`; evaluator cache namespace/version tests passed | PASS |
| A-T2 RNG and regression tests | `A-T2a-report.md`, `A-T2b-report.md`; deterministic stream, direction, exception, CoilTurns/TriggerValue tests | PASS |
| A-T3 constraints, penalty, termination | `A-T3a/b/c-report.md`; thermal rejection, additive penalty, structured reason and stagnation fail-fast tests | PASS |
| A-T4 public contracts/selectors/spec/statistics | `A-T4a/b/c-report.md`; selector, ProblemSpec, run-local collector and concurrency tests | PASS |
| A-T5 docs/engineering closure | `A-T5-report.md`; bilingual self-test, real preset names, install split, executable validation | PASS |
| B-T1 numerical baseline/probe | `B-T1-report.md` plus fresh baseline 154/154; real GPU proof and PeakCurrent concern retained | PASS WITH DEFERRED CONCERN |
| B-T2 production CUDA batch evaluator | `B-T2-report.md`, installed CUDA header, fresh evaluator and baseline runs | PASS |
| B-T3 failure isolation/fallback/statistics | `B-T3-report.md`, 229/229 evaluator assertions, benchmark zero-fallback counters | PASS |
| B-T4 CPU/GPU workflow benchmark | `B-T4-report.md`, fresh `/tmp` JSON and validator, batch 32/128 gates | PASS |
| Cross-track evaluator identity/cache key | A-T1 tests and cache implementation | PASS |
| Cross-track failure isolation semantics | B-T3 protocol/row/whole-batch tests | PASS |
| Cross-track run-local statistics | A-T4c/B-T3 tests and fresh GPU stats | PASS |
| Cross-track public API freeze before B-T2 | ledger says gate passed; B-T2 production API and install consumers compile | PASS |

## Requirement audit — roadmap §7 task/commit items

| # | Planned task / subject | Recorded commit evidence | Decision |
|---:|---|---|---|
| 1 | A-T1 — Fix optimization cache identity | `e54e246` | PASS |
| 2 | A-T2a — Advance optimization random stream | `9786e67` | PASS |
| 3 | A-T2b — Strengthen optimization regression tests | `d722fbf` | PASS |
| 4 | A-T3a — Validate thermal optimization constraints | `2ebb1cc` | PASS |
| 5 | A-T3b — Align optimization constraint semantics | `591f504` plus compatibility fix `b98c85c` | PASS |
| 6 | A-T3c — Unify optimization termination reasons | `f06732e` | PASS |
| 7 | A-T4a — Harden optimization result selectors | `f69e81c` | PASS |
| 8 | A-T4b — Stabilize optimization problem contract | `ae70ab9`, `ae7a5aa` | PASS |
| 9 | A-T4c — Define optimization statistics ownership | `f47d7fb`, `620407f` | PASS |
| 10 | B-T1 — Establish optimization numerical baseline | `a521652` | PASS WITH DEFERRED PEAKCURRENT CONCERN |
| 11 | B-T2 — Add CUDA batch optimization evaluator | `76b6438` | PASS |
| 12 | B-T3 — Integrate GPU optimization metrics | `36daa5c` | PASS |
| 13 | A-T5 — Document optimization public API | `b294888` | PASS |
| 14 | B-T4 — Benchmark GPU optimization workflow | `084eb66` | PASS |
| 15 | W5 close — Validate next optimization phase | This report, final correction commits, and independent approval through `8bc0488` | PASS |

## Deferred items and blockers

Deferred, explicitly documented and not blockers for this fixed-geometry Full
Euler phase: the B-T1 PeakCurrent discrepancy (excluded from production fitness
and benchmark gates), arbitrary-geometry CUDA candidates, thermal CUDA
optimization metrics, CUDA RK4, LookupTable implementation, and future
DE/CMA-ES/Bayesian/MOEA-D families. Existing CUDA compiler warnings were
observed during fresh builds and did not produce a failure.

No roadmap gate remains outstanding. The final independent review and the
subsequent user-authorized cleanup are recorded at the end of this report.

## Post-review fix wave — historical intermediate state

The final-review fix wave addressed the release-blocking cache identity finding
and all listed Minor findings without changing physics, tolerances, or optimizer
search behavior:

- `EvaluationCacheIdentity` is immutable and rejects empty namespace or version.
- `CoilgunOptimizationProblem::cache_identity()` uses an unambiguous full
  serialization of physical state, schema/bindings, excitations/triggers,
  constraints/result schema, timestep, thermal/optimization level, and
  termination policy. Injected callback instances and replacements receive
  unique generations under a mutex; evaluation snapshots the callback safely.
- `CudaBatchEvaluator::cache_identity()` includes the problem identity, backend
  configuration, fallback policy, result-schema version, and private injected
  execution seam identity.
- CPU shared-cache tests and the CUDA production evaluator shared-cache test
  prove cross-configuration misses and same-configuration hits. Focused CPU
  tests passed 20/20 (131 assertions) and 21/21 (83 assertions); CUDA evaluator
  tests passed 22/22 (238 assertions).
- Benchmark validation now rejects malformed source revisions and requires
  non-empty toolchain/workload provenance. Its unit suite passed 12/12.
- API.md/API_cn.md now state that optimization statistics use a fresh run-local
  collector and never a lifetime before/after snapshot delta. The obsolete
  immutable-identity deferred item was removed from the ledgers.

At this intermediate point, full four-preset verification and independent
re-review were still pending. Both were completed by the final evidence and
approval recorded below.

## Post-fix fresh evidence (2026-09-17)

After commit `2875457eb554c8688ed25a4e1821c4483d90ebe5`, the required clean
worktree reruns completed with these discovered totals:

```text
cmake --preset cpu-debug && cmake --build --preset cpu-debug -j2 && ctest --preset cpu-debug --output-on-failure
  exit 0; 35/35 passed; 0 skipped; 0 failed
cmake --preset cpu-release && cmake --build --preset cpu-release -j2 && ctest --preset cpu-release --output-on-failure
  exit 0; 35/35 passed; 0 skipped; 0 failed
cmake --preset cuda-debug && cmake --build --preset cuda-debug -j2 && ctest --preset cuda-debug --output-on-failure
  exit 0; 53/53 passed; 0 skipped; 0 failed
cmake --preset cuda-release && cmake --build --preset cuda-release -j2 && ctest --preset cuda-release --output-on-failure
  exit 0; 53/53 passed; 0 skipped; 0 failed
```

Independent CUDA Release production checks passed:

```text
test_gpu_optimization_baseline: 1 test case, 154/154 assertions passed
test_cuda_batch_evaluator:      22 test cases, 238/238 assertions passed
```

The new CUDA cache test exercised two physical problem configurations against
one shared cache, proving both cross-configuration misses and a repeated
same-configuration hit. All measured CUDA rows reported `gpu_executed=true`,
`backend=direct`, `precision=full`, `order_ok=true`, `numerical_ok=true`, and
`fallback_events=0`.

The fresh benchmark was written to
`/tmp/optimization-gpu-workflow-final-fix.json` and validated as
`valid B-T4 schema v1`. It recorded source revision `2875457`, clean state,
toolchain `GNU 16.2.1` / CMake `4.4.3` / Release, and workload provenance
`optimization-gpu-workflow` / `fixed-geometry-euler-full` /
`bench_gpu_optimization_workflow`. Batch results were:

| Batch | CPU cand/s | GPU cand/s | Speedup | GPU/backend | Order | Numerical | Fallback events |
|---:|---:|---:|---:|---|---|---|---:|
| 1 | 28.5909 | 10.2300 | 0.3578x | true/direct | true | true | 0 |
| 8 | 28.5225 | 90.2896 | 3.1656x | true/direct | true | true | 0 |
| 32 | 28.6703 | 187.6973 | 6.5468x | true/direct | true | true | 0 |
| 128 | 28.6464 | 251.4378 | 8.7773x | true/direct | true | true | 0 |

The machine-readable decision kept both throughput gates true. GPU optimizer
statistics recorded 32/32 successful evaluations, 32/32 GPU-executed rows,
four GPU batches, zero failed batches, zero GPU fallbacks, and zero CPU repairs;
the Reference recheck remained valid with error `1.231653667943533e-16` versus
tolerance `6.294880560104472e-08`.

The post-fix documentation and validator gates passed:

```text
check_api_optimization_docs.py --self-test: PASS
test_validate_optimization_gpu_workflow.py: 12 tests, OK
git diff --check: exit 0
```

CPU and CUDA Release installs passed `verify_install_headers.py` (CPU-only and
CUDA modes). Out-of-tree `find_package(coilgun CONFIG REQUIRED)` consumers
built and ran using only the installed prefixes: CPU printed `250`, CUDA
printed `strict`. These were intermediate results; the final review status is
recorded at the end of this report.

## Post-review cache evaluation snapshot fix (2026-09-17)

The remaining Important review finding is fixed in the final-review wave.
`BatchEvaluator::evaluation_snapshot()` is the generic compatibility boundary
for pairing an immutable cache identity with the evaluator callable that
produced the batch. Legacy subclasses retain the default adapter; stateful
evaluators override the hook. `CachedBatchEvaluator` now captures one snapshot
before cache lookup and uses its bound callable for all misses and isolated
retries. `StatisticsBatchEvaluator` forwards and wraps the snapshot, so
statistics and cache wrappers retain their normal forwarding behavior.
`CoilgunOptimizationProblem` copies the callback and generation under its
mutex, builds the canonical identity from that same generation, then releases
the mutex before executing the callback. Callback replacement therefore cannot
pair a new callback with an old cache key, and no long-running evaluation holds
the replacement mutex.

The barrier regression test
`cached evaluation snapshots keep callback identity paired with its result`
forces the old-identity/new-generation interleaving against a shared cache;
the fixed path returns generation 1 and a second evaluator with the old
identity observes the old result rather than a polluted generation-2 result.
The focused evaluator suite passed 21/21 cases and 137/137 assertions; the
Coilgun optimization suite passed 21/21 cases and 83/83 assertions.

Fresh post-fix configure/build/full-CTest evidence:

```text
cpu-debug:    35/35 passed
cpu-release:  35/35 passed
cuda-debug:   53/53 passed
cuda-release: 53/53 passed
```

Real-GPU focused checks passed: `test_cuda_batch_evaluator` 22/22 cases and
238/238 assertions; `test_gpu_optimization_baseline` 1/1 case and 154/154
assertions. The CUDA Release benchmark was rerun from clean source revision
`c367ae0` and passed `validate_optimization_gpu_workflow.py`. It recorded GNU
16.2.1, CMake 4.4.3, Release, workload provenance
`optimization-gpu-workflow` / `fixed-geometry-euler-full` /
`bench_gpu_optimization_workflow`, direct GPU execution and zero fallback
events for all rows. Batch 32 and 128 GPU throughput gates were true. CPU and
CUDA release-library installed-header checks and out-of-tree consumers also
passed. The API checker, validator unit suite (12/12), and `git diff --check`
passed. No TSAN preset or existing sanitizer flow is defined in this project.

## Final whole-branch correction-wave verification (2026-09-21)

The correction wave was rerun through the complete validation gate. The
integration and benchmark callback fixtures now emit the declared
`velocity_floor` report, so strict callback validation does not turn a valid
test workload into an unintended fallback.

```text
cpu-debug:    35/35 passed
cpu-release:  35/35 passed
cuda-debug:   53/53 passed
cuda-release: 53/53 passed
CUDA focused label: 19/19 passed
test_cuda_batch_evaluator: 24/24 cases, 244/244 assertions
test_gpu_optimization_baseline: 1/1 case, 154/154 assertions
```

The CUDA Release benchmark ran on the RTX 5080 Laptop GPU with direct/full
GPU execution, stable ordering, numerical parity, and zero fallback events;
the schema validator accepted the B-T4 JSON. Batch-32 and batch-128
throughput gates passed (9.8805x and 13.0954x GPU per-candidate speedups in
this run), and the GPU optimizer recorded 32/32 successful, GPU-executed
evaluations with four GPU batches and zero repairs. The bilingual API checker,
validator unit suite (12/12), `git diff --check`, CPU/CUDA installed-header
checks, and out-of-tree `find_package(coilgun CONFIG REQUIRED)` CPU/CUDA
consumers all passed.

TDD evidence for this wave includes the RED callback-schema mutations on the
pre-fix implementation (missing/extra/reordered or mismatched constraint
reports were incorrectly accepted), followed by GREEN fallback assertions for
all mutations. The new CUDA source-problem destruction case and the stronger
active-call last-owner barrier are green after the ownership changes.

## Final whole-branch review correction wave (2026-09-21)

The review's Critical CUDA lifetime finding is addressed by making every
`CudaBatchEvaluator` own an immutable `CoilgunOptimizationProblem` copy built
from the source schema and configuration. The evaluator's snapshot therefore
keeps the complete derived evaluator and its CPU-fallback dependency alive
after the source problem is destroyed; cache identity, schema decoding, GPU
execution, and fallback evaluation all use the owned copy.

The callback boundary now accepts a successful GPU result only when it exactly
matches the declared objective and ordered constraint schema: objective count,
ID, direction, and finite value; constraint count/order, ID, kind, relation,
bounds, priority, finite values, and derived violation/satisfaction fields.
The derived normalized violation check enforces the declared scale. Any
malformed callback output is rejected and repaired through the existing CPU
fallback path before caching.

All repository no-op-deleter/aliasing evaluator patterns were removed from the
integration, benchmark, and CUDA tests. Documentation now states that shared
ownership cannot be simulated by a no-op deleter. The active-call regression
now moves the only snapshot copy into the worker, allowing the main thread to
release the last external snapshot and owner while evaluation is blocked.

The review's reported `evaluation_snapshot()` compatibility concern was
verified as a false positive: `git show 7e43c4e:include/coilgun/optimization/evaluator.hpp`
contains neither `BatchEvaluationSnapshot` nor `evaluation_snapshot`; this is
new, unreleased API and the final boundary remains intentionally non-overridable.

## Post-review lifetime-safe snapshot ownership correction (2026-09-21)

The remaining Important finding is fixed by removing the base-destructor drain
protocol. `BatchEvaluator::evaluation_snapshot()` remains the final,
non-overridable lifecycle boundary around the protected
`make_evaluation_snapshot()` hook, but it now first requires genuine
`std::shared_ptr` management and then captures a strong owner in the returned
callable. This keeps the complete derived object alive through saved and
concurrent calls, including the last external-owner reset barrier. A stack or
other unmanaged evaluator is rejected immediately with stable
`std::logic_error` text (`BatchEvaluator::evaluation_snapshot requires shared
ownership`); it never returns a raw-`this` callable. The obsolete
`evaluation_snapshot_expired` post-destruction path and active-call drain were
removed. Cached/statistics/CUDA wrappers and `CoilgunOptimizationProblem` keep
their callback, identity, cache, and statistics forwarding behavior.

TDD regression evidence:

```text
RED (HEAD 930e3e8): unmanaged snapshot test did not throw
GREEN: test_optimization_evaluator 25/25 cases, 148/148 assertions
GREEN: test_coilgun_optimization 24/24 cases, 98/98 assertions
```

The generic suite covers shared snapshots after external reset, an active
snapshot/barrier reset with delayed destruction and no deadlock, explicit
unmanaged rejection, and the shared Cached/Statistics wrapper chain. The
Coilgun suite covers unmanaged rejection, a real callback snapshot after
external reset, and the shared-cache callback-generation barrier. The CUDA
focused suite adds a shared `CudaBatchEvaluator` snapshot regression.

Fresh verification after this correction is green:

```text
cpu-debug:    35/35 passed
cpu-release:  35/35 passed
cuda-debug:   53/53 passed
cuda-release: 53/53 passed
test_cuda_batch_evaluator:      23/23 cases, 241/241 assertions
test_gpu_optimization_baseline: 1/1 case, 154/154 assertions
```

The clean CUDA Release benchmark JSON was accepted by
`validate_optimization_gpu_workflow.py` as `valid B-T4 schema v1`. It recorded
GNU 16.2.1, CMake 4.4.3, Release, direct/full GPU execution, stable ordering,
numerical parity, and zero fallback events:

| Batch | CPU cand/s | GPU cand/s | Speedup |
|---:|---:|---:|---:|
| 1 | 28.9429 | 20.7556 | 0.7171x |
| 8 | 29.1462 | 163.6140 | 5.6136x |
| 32 | 29.0882 | 286.3166 | 9.8430x |
| 128 | 29.0886 | 382.0532 | 13.1341x |

Both batch-32 and batch-128 throughput decisions were true. The GPU optimizer
recorded 32/32 successful evaluations, 32/32 GPU-executed rows, four GPU
batches, zero failed batches/fallbacks/CPU repairs, and a valid CPU Reference
recheck (error `1.231653667943533e-16` versus tolerance
`6.294880560104472e-08`). The API checker and negative self-test passed, the
validator unit suite passed 12/12, and `git diff --check` was clean. CPU-only
and CUDA installed-header checks passed; out-of-tree `find_package(coilgun
CONFIG REQUIRED)` consumers built and ran using only installed prefixes
(outputs `250` and `0`).

## Final independent whole-branch approval and cleanup (2026-09-21)

The final Luna High read-only review covered `7e43c4e..8bc0488` using the
complete review package `review-7e43c4e..8bc0488.diff`, the design and roadmap,
and this validation report. Verdict: **Approved**, with 0 Critical, 0 Important,
and 0 Minor findings.

The reviewer independently confirmed that CUDA snapshots own their complete
problem dependency, shared snapshot ownership has no raw-owner escape or
repository-local no-op-deleter path, callback results are checked against the
complete declared objective and constraint schema before caching, and the
requested Auto-routing, Pareto-only default, explicit representative selection,
CPU/CUDA integration, installation, bilingual documentation, and benchmark
contracts remain satisfied.

After approval, the clean, fully integrated intermediate worktrees and local
branches `feature/cuda-batch-optimization` and
`feature/optimization-quality-closure` were removed. The dirty primary
`feature/optimization-next-phase-roadmap` worktree and dirty
`fix/opt-t6-convergence` worktree were intentionally preserved without edits.
The deliverable remains on `feature/optimization-next-phase` in
`.worktrees/optimization-next-phase-final`.
