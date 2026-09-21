# AGENTS.md — coligunCalc

## What this is
Multi-stage synchronous induction coilgun simulation library using the **current filament method** (CFM). C++20, GPLv3. One static library `coilgun` (CPU/OpenMP) plus an optional CUDA backend (`coilgun_cuda`). Eigen, Boost.Math and doctest are fetched by CMake; OpenMP is required.

## Current state
- Active branch: `feature/optimizationAlgorithm` (2 commits ahead of `origin/feature/optimizationAlgorithm`). `main` / `dev` are the stable lines.
- Physics foundation, CPU simulation engine, CUDA engine, and the optimization framework are implemented. The optimization work is driven task-by-task by the SDD plan below; OPT-T1…T10 are complete, T11 (public API + docs integration) and T12 (workflow validation) remain.
- `docs/PerformanceOptimizationFindings.md` is an untracked, in-progress CPU/GPU performance report (measured on this machine); it is not yet committed.
- The optimization module is **not yet exposed by `coilgun/coilgun.hpp`** and there are no install/export rules yet — that is OPT-T11. Include `coilgun/optimization/*.hpp` directly until then.
- The Python prototype (`basics.py`, `core.py`, `__init__.py`) is gone from the tree but preserved in history (`git show 563c1d1:basics.py`).

## Repository layout
```
include/coilgun/        Public headers
  core/types.hpp        Common types
  physics/              Elliptic/Struve, quadrature, T(q,p) tables, self/mutual inductance, LRU cache
  components/           DrivingCoil, Armature (m×n filament discretisation)
  simulation/           Excitation, steppers, SimState/MultiStageState, results, termination, timing
  simulation/cuda/      GPU backend public surface (gpu_engine, execution policy/report, SimBatch, …)
  optimization/         GA/NSGA-II framework: variables, objectives, constraints, operators,
                        evaluators, routing, selectors, CoilgunOptimizationProblem adapter
  tools/                t_table_workers.hpp (build-tool helper, not part of the runtime API)
  coilgun.hpp           CPU umbrella; coilgun_cuda.hpp adds the CUDA umbrella
src/                    Library implementation (mirrors include/coilgun/); src/cuda/ holds .cu sources
tests/                  doctest suites registered in tests/CMakeLists.txt (CTest)
tools/                  T(q,p) lookup table generator (generate_t_table.cpp)
scripts/                verify_t_table.py
docs/                   API.md / API_cn.md, NumericalModel.md, benchmarks/, superpowers/ (plans, specs)
.references/            Reference PDFs + MinerU OCR output (gitignored, local-only)
.worktrees/             Per-task git worktrees (gitignored); do not edit from the main tree
```

## Build
CPU (default; tests ON):
```sh
cmake --preset cpu-debug && cmake --build --preset cpu-debug && ctest --preset cpu-debug
```
CUDA (needs CUDA Toolkit ≥ 12.8 for Blackwell sm_120, NVIDIA CC ≥ 6.0):
```sh
cmake --preset cuda-debug && cmake --build --preset cuda-debug && ctest --preset cuda-debug -L gpu
```
Presets: `cpu-debug`, `cpu-release`, `cuda-debug`, `cuda-release`, `cpu-release-library`, `cuda-release-library`. All build out-of-source under `build/<preset>/`. Requires CMake ≥ 3.20 and a C++20 compiler with OpenMP.

Other options:
- `COILGUN_ENABLE_CPU_PHASE_TIMING=ON` — opt-in CPU derivative phase instrumentation (off by default, no cost when off).
- `COILGUN_BUILD_GENERATOR=ON` — builds `tools/generate_t_table.cpp`.
- `COILGUN_BUILD_TESTS=ON` (default), `COILGUN_ENABLE_CUDA=OFF` (default).

Test labels: `fast`/`quick`, `validation`, `slow`/`integration`, `gpu`. GPU tests carry `RESOURCE_LOCK gpu` and skip when no device is present. Quick checks: `ctest --preset cpu-debug -L fast`. Do not freeze test counts in docs — CTest discovers them per preset.

## Known limitations and gotchas
- **GPU is Euler-only.** `GpuSingleStageSim` / `GpuMultiStageSim` throw `std::logic_error` for `RK4Stepper` (`src/cuda/gpu_single_stage_sim.cu`, `src/cuda/gpu_multi_stage_sim.cu`); use the CPU simulators for RK4.
- **A requested backend is not proof of GPU execution.** Assert `ExecutionReport::gpu_executed == true` and `backend != BackendMode::Fallback` in any GPU validation.
- `BackendMode::Persistent` is unreachable (`supports_persistent_control_stream == false`) and always falls back; the fallback path is slower than running the CPU engine directly.
- `OptimizationLevel::LookupTable` is reserved: it currently takes the same runtime path as `Reference`.
- Root contains stale in-source CMake artifacts (`CMakeCache.txt`, `CMakeFiles/`, `Makefile`) and stale preset dirs under `build/` (e.g. `ninja-debug`, `cpu-only`) from earlier setups. They are gitignored; ignore them and use the presets above.
- `.worktrees/` and some `/tmp/coligun*` worktrees are per-task branches; several `/tmp` entries are prunable. Never run long builds or edits against the main tree when a task is assigned a worktree.
- `.ckb/`, `.depwire/`, `.reflex/` are local agent-tool caches (gitignored). `.superpowers/sdd/optimization/` holds task briefs/reports and `progress.md` and is intentionally tracked despite the `.superpowers/` ignore rule — keep the ledger updated when executing plan tasks.

## Key documents
- `docs/NumericalModel.md` — canonical physics reference: filament discretization, elliptic-integral kernel, circuit equations, steppers.
- `docs/API.md` / `docs/API_cn.md` — C++ API reference (EN/CN).
- `docs/superpowers/specs/2026-09-06-optimization-algorithm-design.md` — optimization module design.
- `docs/superpowers/plans/2026-09-06-optimization-algorithm-execution-plan.md` — OPT-T1…T12 execution plan and gates.
- `docs/benchmarks/GPU-Benchmark-Schema.md`, `docs/benchmarks/Optional CI Plan.md` — GPU benchmark schema and CI plan.

## Documentation sync rule
**API.md and API_cn.md must always be kept in sync.** Any change to one must be mirrored in the other. Same applies to README.md ↔ README_cn.md. Both language versions must contain equivalent content.

## Commit conventions
History uses capitalized imperative subjects, one logical change per commit: `Add …`, `Fix …`, `Implement …`, `Test …`, `Update …`, `Integrate …`, `Isolate …`. Feature branches use `feature/…`, fixes `fix/…`. Keep commits small and match the subject named in the task brief when one is given.

## MinerU OCR pipeline for .references/
When OCR-ing PDFs in `.references/` via local MinerU API (`http://127.0.0.1:8000`):
- **Backend**: Always `hybrid-engine` (CUDA), not `pipeline`
- **Output**: one folder per PDF under `.references/<basename>/`:
  - `md_content` → `<basename>.md`
  - `middle_json` → `middle_json.json`
  - `content_list` → `content_list.json`
  - OCR images → `images/`
- **No root pollution**: move any auto-created `output/` into the target ref folder and delete from root
