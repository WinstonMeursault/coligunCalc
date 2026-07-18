# Post-Audit Fix Plan — Comprehensive Remediation

> Date: 2026-07-18 | Scope: All audit findings + design spec alignment

## Background

A full codebase audit (see `docs/audit_report.md`) identified 17 issues across documentation, code quality, and tests. Additionally, the GPU optimization design spec (`docs/superpowers/specs/2026-07-18-gpu-optimization-design.md`) defines a target architecture for persistent kernel, GpuOptLevel 3-tier, and auto-fallback — which is partially implemented but has documentation drift and test failures.

## Fix Inventory (17 items → 7 categories)

### Category A: Documentation Sync (API.md / API_cn.md)

| ID | Issue | Action |
|----|-------|--------|
| A1 | CN GpuBackend table missing `use_persistent` row | Add row to CN table |
| A2 | CN Known Limitations missing "Double buffering" row | Add row to CN table |
| A3 | EN/CN Known Limitations row 3 directly contradicts (persistent vs batch kernel status) | Rewrite both to reflect actual state: persistent kernel is experimental+verified on RTX 5080; batch kernel (P0 grid) is blocked |
| A4 | Both files claim GpuOptLevel has "2 levels" (has 3) | Fix "2 levels" → "3 levels" in both EN and CN (line in CPU-vs-GPU diff table) |

### Category B: README Sync (README.md / README_cn.md)

| ID | Issue | Action |
|----|-------|--------|
| R1 | Claims "17 suites, all passing" | Update to "18 suites, 16 passing" with accurate test status |
| R2 | References non-existent docs (multi_stage_sim_design.md, test_dataset_82mm_coilgun.md) | Remove dead references |

### Category C: Code Quality Fixes

| ID | Issue | Action |
|----|-------|--------|
| M1 | `src/coilgun.cpp` empty placeholder stub | Fill with `#include <coilgun/coilgun.hpp>` (ensure symbol export) or delete |
| M3 | GpuMultiStageSim header says M_cc "omitted for simplicity" — WRONG (M_cc IS computed and used) | Fix header comment |

### Category D: Test Failures

| ID | Issue | Action |
|----|-------|--------|
| H1 | `test_gpu_vs_cpu_multi` FAILS (GPU 1.80 vs CPU 1.64 m/s, ~10% error) | Debug root cause; likely persistent kernel doorbell protocol bug in multi-stage path |
| H2 | `test_gpu_batch` TIMEOUT (120s) | Debug; may be same root cause as H1 or test configuration issue |

### Category E: Architecture Verification

| Item | Action |
|------|--------|
| Persistent kernel exists in all 3 GPU classes | Verify all paths compile and function correctly |
| Auto-fallback on cudaHostAllocMapped failure | Verify fallback path works |
| GpuOptLevel 3-tier enum | Already implemented (Standard=0, Full=1, Aggressive=2) |
| GpuBackend::use_persistent | Already implemented (default true) |

---

## Root Cause Analysis for Test Failures

### H1: GPU multi-stage muzzle velocity mismatch (~10%)

**Hypothesis**: The persistent kernel path in `GpuMultiStageSim::compute_M1_dM1_persistent()` has a doorbell synchronization issue.

**Current implementation** (from source analysis):
- Host sets `seps[0..active_count-1]` with separation values
- Host sets `doorbell[0..active_count-1] = 1`
- Host sets `*active_pairs = active_count`
- Host waits for `doorbell[0..active_count-1]` to clear (all 0)
- **Bug**: The kernel waits for `doorbell[idx] == 0 && *shutdown == 0`. On the FIRST call, all doorbells are initialized to 0. The kernel will start processing BEFORE the host sets doorbells to 1 — this is a race condition!
- Host then reads `out_M[0..S*F-1]` and `out_dM[0..S*F-1]` for ALL pairs

**Fix**: Add a `batch_id` synchronization token (as per the design spec §1.4): host increments `batch_id` AFTER filling all data, kernel spins on `batch_id == last_batch_id`. This eliminates the race condition. The batch_id pattern is documented in the superpowers spec but the current implementation uses the older doorbell-only protocol from `persistent_kernel.cu`.

Alternatively, simpler fix: ensure doorbells are set to 0 before host fills data, then set to 1 after all data is ready. The kernel then waits for doorbell[idx]==1. But the key issue is the initial state and ensuring host completes writing before device reads.

**Actual implementation review needed**: Read the exact current `compute_M1_dM1_persistent()` in `gpu_multi_stage_sim.cu` line-by-line.

### H2: GPU batch timeout

Likely the same root cause (persistent kernel sync issue) plus the batch test runs 10 simulations × ~20000 steps = much longer runtime.

---

## Implementation Plan (Ordered by Dependency)

### Step 0: Deep-dive read of failing code paths
Read the exact current implementations of the persistent paths in multi-stage and batch to confirm root cause.

### Step 1: Fix persistent kernel synchronization
- Implement batch_id-based protocol (matching design spec §1.4) OR
- Fix doorbell protocol to eliminate race condition
- Apply fix to all three GPU classes (single-stage already works, but audit for consistency)

### Step 2: Fix documentation (A1–A4, R1–R2)
- API.md: 3 edits
- API_cn.md: 4 edits
- README.md: 2 edits
- README_cn.md: 2 edits

### Step 3: Fix code quality (M1, M3)
- `src/coilgun.cpp` — fill or delete
- `gpu_multi_stage_sim.hpp` — fix line 37 comment

### Step 4: Full build + test verification
- `cmake --build --preset ninja-cuda-debug`
- `ctest --preset debug` — all 18 tests must pass
- `compute-sanitizer` on GPU tests

---

## Success Criteria

1. ✅ All 18 tests pass (no failures, no timeouts)
2. ✅ API.md and API_cn.md are structurally identical for all documented API surfaces
3. ✅ README.md and README_cn.md accurately describe test suite status
4. ✅ No stale/incorrect comments in GPU headers
5. ✅ `src/coilgun.cpp` is either populated or removed
6. ✅ `cmake --build --preset ninja-cuda-debug` compiles with 0 errors
7. ✅ GPU persistent kernel path is thread-safe and deterministic
