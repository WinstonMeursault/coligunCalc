# coligunCalc 综合审计报告

> 审计日期：2026-07-18
> 审计范围：docs/ 设计文档 vs 全部源码、API.md vs API_cn.md 一致性、README 准确性

---

## 一、设计文档覆盖度：NumericalModel.md → C++ 实现

C++ 实现**完整覆盖**全部 9 个章节。无核心物理层面的遗漏。

| § | 内容 | 实现文件 | 状态 |
|---|------|----------|:----:|
| 1 | 电感法、洛伦兹力原理 | `mutual_inductance.cpp` + force computation | ✅ |
| 2 | 电流丝法 m×n 离散化 | `components/armature.cpp` Armature 类 | ✅ |
| 3.1 | 单级等效电流环电路方程 | `simulation/single_stage_sim.cpp` | ✅ |
| 3.2 | 电流丝模型电路方程 | `simulation/single_stage_sim.cpp` | ✅ |
| 3.3 | 多级电路方程（Q 截断矩阵） | `simulation/multi_stage_sim.cpp` | ✅ |
| 3.4 | 续流二极管模型 | `simulation/excitation.cpp` CrowbarExcitation | ✅ |
| 4.1 | 自感 T(q,p) 表 + Bessel/Struve 精确 | `physics/self_inductance.cpp`, `lookup_tables.cpp` | ✅ |
| 4.2 | 丝级互感椭圆积分 K(k)/E(k) | `physics/mutual_inductance.cpp` + `elliptic.cpp` | ✅ |
| 4.3 | 线圈级互感 4D Gauss-Legendre | `physics/mutual_inductance.cpp` | ✅ |
| 4.4 | 互感梯度 dM/dx | `physics/mutual_inductance.cpp` | ✅ |
| 4.4.1 | 对称性与奇偶优化 | 丝级函数根据文档实现 | ✅ |
| 4.4.2 | AGM 椭圆积分算法 | 委托给 Boost.Math | ✅ |
| 5.1–5.3 | 力、运动方程、前向 Euler | simulation engine | ✅ |
| 5.3–5.4 | RK4 步进器、触发位置 | `time_stepper.hpp` RK4Stepper、trigger | ✅ |
| 6.1–6.4 | 绝热温升、cp(T)、rho(T)、耦合迭代 | `physics/constants.cpp` + thermal mode | ✅ |
| 7 | 效率公式 η = mv²/(nCU²) | `SimSummary` / `MultiStageSummary` | ✅ |
| 8.1–8.4 | 完整算法流程 | simulation engine | ✅ |

---

## 二、代码层面发现的问题

### 2.1 功能性问题（HIGH）

| # | 问题 | 位置 | 详情 |
|---|------|------|------|
| **H1** | `test_gpu_vs_cpu_multi` **测试失败** | `tests/test_gpu_vs_cpu_multi.cpp:57` | GPU 多级仿真出口速度 1.80 m/s vs CPU 1.64 m/s，误差 ~10%，远超 ε=0.5% 阈值 |
| **H2** | `test_gpu_batch` **测试超时** | 测试目标 #17 | 120 秒内未完成，断言从未运行 |

### 2.2 代码质量问题（MEDIUM）

| # | 问题 | 位置 | 详情 |
|---|------|------|------|
| **M1** | `src/coilgun.cpp` 是空桩文件 | `src/coilgun.cpp` | 仅一行 `// Placeholder — populated in subsequent commits`。应填充或删除 |
| **M2** | GPU 双向缓冲是死代码 | `gpu_single_stage_sim.hpp` | `dbl_M1_`、`dbl_dM1_`、`dbl_buf_active_` 声明但从未激活（Always `false`） |

### 2.3 注释与代码不一致（MEDIUM）

| # | 问题 | 位置 | 详情 |
|---|------|------|------|
| **M3** | GpuMultiStageSim 头文件注释说 `M_cc_`"已省略以简化" | `include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp:37` | 但实际实现 (`src/cuda/gpu_multi_stage_sim.cu:108–120, :325`) **确实计算并使用了** M_cc |

---

## 三、GPU 后端的已知简化与限制

### 3.1 已文档化的简化（非问题）

| 简化项 | 文档位置 | 说明 |
|--------|----------|------|
| GPU 自适应 GL 阶数移除（仅 n_nodes=9） | API.md §Known Limitations | 4 点 GL 导致非确定性浮点漂移 |
| GPU 热模式串行 CPU 执行 | API.md §GpuMultiStageSim | 占比 < 运行时间 1% |
| 批量仿真非单次 grid launch | API.md §SimBatch 注释 | 持久化 kernel 优化计划中 |

### 3.2 未文档化或文档矛盾的问题

| # | 严重度 | 问题 | 详情 |
|---|--------|------|------|
| **D1** | HIGH | 持久化 kernel 状态在 EN/CN 文档中**直接矛盾** | EN 版 API.md: "Experimental. Verified working on RTX 5080 (CC 12.0)"；CN 版 API_cn.md: "持久化 kernel 优化（D1-D2）是计划中的替代方案"——描述的是同一个特性的两个不同实现状态 |
| **D2** | MED | GPU 双向缓冲（D5）在 EN 版列为"Planned"但 CN 版 Known Limitations 完全缺失此行 | EN API.md line 1588 提到了 Double buffering；CN API_cn.md 没有对应行 |

---

## 四、API.md vs API_cn.md 一致性差异

### 4.1 单一来源缺失（MEDIUM）

| # | 缺失项 | 缺失方 | 位置 |
|---|--------|--------|------|
| **A1** | GpuBackend 字段表缺少 `use_persistent` 行 | CN (API_cn.md) | EN line 1064–1070 有5行，CN line 1064–1069 仅4行（struct 定义含 `use_persistent`，仅说明表缺失） |
| **A2** | Known Limitations 表缺少 "Double buffering" 行 | CN (API_cn.md) | EN line 1588 有此项，CN 无对应行 |

### 4.2 直接矛盾（HIGH）

| # | 矛盾内容 | EN 描述 | CN 描述 |
|---|----------|---------|---------|
| **A3** | Known Limitations 表第3行 | `Persistent kernel`: "**Experimental.** Uses mapped memory with doorbell protocol. Verified working on RTX 5080 (CC 12.0)..." | `Batch kernel`: "被 nvcc 跨模块 struct-layout 问题**阻塞**。持久化 kernel 优化（D1-D2）是**计划中的**替代方案" |

### 4.3 两文件共同错误（LOW）

| # | 问题 | 详情 |
|---|------|------|
| **A4** | CPU-vs-GPU 差异表中称 GpuOptLevel 有"2 levels" | 两文件均如此表述，但枚举实际定义 `Standard=0, Full=1, Aggressive=2` 共三级 |

---

## 五、README.md vs README_cn.md

### 5.1 两版一致性

EN 版与 CN 版 README 内容**一一对应**，未发现结构性差异或翻译矛盾。所有章节、预设表、物理公式表、测试套件列表均一致。

### 5.2 两版共同错误

| # | 问题 | 详情 |
|---|------|------|
| **R1** | 声称"17 suites, all passing" | 实际是 **18 个** CTest 测试目标（tests 1–18），其中 `test_gpu_vs_cpu_multi` 失败、`test_gpu_batch` 超时 |
| **R2** | 引用不存在的文档 | `docs/multi_stage_sim_design.md` 和 `docs/test_dataset_82mm_coilgun.md` 在 docs/ 目录下不存在 |

---

## 六、测试结果总览

| 状态 | 数量 | 测试 |
|:----:|:----:|------|
| ✅ 通过 | 16 | test_elliptic, test_lookup, test_struve, test_quadrature, test_self_inductance, test_mutual_inductance, test_driving_coil, test_armature, test_integration, test_single_stage_sim, test_multi_stage_sim, test_gpu_elliptic, test_gpu_filament, test_gpu_coil_pair, test_gpu_vs_cpu_single, test_gpu_sim_batch |
| ❌ 失败 | 1 | **test_gpu_vs_cpu_multi** — GPU 多级出口速度 1.80 vs CPU 1.64 m/s（~10% 误差） |
| ⏱️ 超时 | 1 | **test_gpu_batch** — 120 秒内未完成 |

---

## 七、总结

### 正面

- **物理实现完整**：所有设计文档要求的公式和算法均已落地到代码中，没有理论上的偷工减料
- **CPU 路径质量高**：自感双轨（T表+精确）、Struve 函数 SciPy 等效、LRU 缓存、RK4、多级电路均仔细实现
- **16/18 测试通过**：核心物理层和组件层测试全部通过
- **GPU 后端实际状况优于注释**：M_cc 在代码中已实现（头文件注释有误需修复）

### 需修复（按优先级）

| 优先级 | 问题编号 | 说明 |
|--------|----------|------|
| **P0** | A3, D1 | 持久化 kernel EN/CN 文档矛盾 —— 需统一描述为实际状态 |
| **P0** | H1 | GPU 多级仿真测试失败 —— 误差 10%，需调查根因 |
| **P1** | R1 | README "17 suites, all passing" 不实 —— 更新为准确数字和状态 |
| **P1** | A1, A2, D2 | API_cn.md 缺失的表格行 —— 从 EN 版补充 |
| **P1** | A4 | 两文件 GpuOptLevel "2 levels" → 修正为 "3 levels" |
| **P2** | H2 | GPU batch 测试超时 —— 调查是否为硬件/超时阈值问题 |
| **P2** | R2 | README 引用不存在文档 —— 删除死链或补建文档 |
| **P2** | M3 | GpuMultiStageSim 头文件过时注释 —— 更新为"M_cc included" |
| **P3** | M1 | `src/coilgun.cpp` 空桩 —— 填充或删除 |
| **P3** | M2 | GPU 双向缓冲死代码 —— 启用或清理 |
