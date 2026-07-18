# CUDA Acceleration — Subagent Dispatch Schedule

> Session: 2026-07-17 | Branch: `feature/cudaAcceleration`

## Doxygen 规范

所有新增文件必须使用 Doxygen 风格注释：

```cpp
/**
 * @file filename.hpp
 * @brief 一句话描述
 * @author Winston Meursault
 */

/**
 * @brief 函数/类说明
 * @param x 参数说明
 * @return 返回值说明
 */
```

参考现有文件 `include/coilgun/physics/elliptic.hpp` 的风格。

## Wave 1 — 源码创建（5 parallel）

Wait: 等待全部 agent 完成 → main agent review → 为每个 agent 创建独立 feature commit

| ID | Tasks | Files | Commits |
|---|---|---|---|
| S1 | T1, T4 | CMakeLists.txt, src/cuda/CMakeLists.txt | 2 |
| S2 | T2, T3, T5, T7 | elliptic.cuh, mutual_inductance.cuh, gpu_elliptic.cu, gpu_mutual_inductance.cu | 4 |
| S3 | T8, T11, T12 | gpu_backend.hpp, gpu_adaptor.hpp, gpu_adaptor.cu | 3 |
| S4 | T13, T14, T16, T17 | gpu_single_stage_sim.{hpp,cu}, gpu_multi_stage_sim.{hpp,cu} | 4 |
| S5 | T20.5 | coilgun_cuda.hpp | 1 |

## Wave 2 — 测试创建（4 parallel）

Wait: 等待全部 agent 完成 → main agent review → commit

| ID | Tasks | Files | Commits |
|---|---|---|---|
| S6 | T6, T9, T10, T20.6 | tests/CMakeLists.txt, test_gpu_elliptic.cpp, test_gpu_filament.cpp, test_gpu_coil_pair.cpp | 4 |
| S7 | T15 | test_gpu_vs_cpu_single.cpp | 1 |
| S8 | T18 | test_gpu_vs_cpu_multi.cpp | 1 |
| S9 | T19 | test_gpu_batch.cpp | 1 |

## Wave 3 — 文档（1 serial，不可并行）

Wait: 完成 → main agent review

| ID | Tasks | Files | Commits |
|---|---|---|---|
| S10 | T20, T21 | docs/API.md, docs/API_cn.md, docs/CUDA-feasibility.md | 2 |

## Wave 4 — 验证（1 serial）

Wait: 完成 → main agent review

| ID | Tasks |
|---|---|
| S11 | T22: full build + all tests |

## 进度

| Wave | Status | Agent IDs |
|---|---|---|
| Wave 1 | done | S1-S5 | S1-S5 |
| Wave 2 | done (S9 compile fix applied) | S6-S9 | S6-S9 |
| Wave 3 | done | S10 |
| Wave 4 | done | S11 |
