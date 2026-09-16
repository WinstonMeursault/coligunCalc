# coligunCalc 性能分析报告（数值模型与仿真计算模块：CPU + GPU）

**日期：** 2026-09-08
**范围：** CPU 数值内核与仿真主循环（`src/physics`、`src/simulation`、`src/optimization`）+ CUDA 后端（`src/cuda`，实测）。
**方法：** 源码走查 + 自建微基准（见 §7）。除标注"推测"外，所有数字均为本机实测。

**测试环境**

| 项目 | 配置 |
|---|---|
| CPU | Intel Core Ultra 9 275HX，8 P-core（0–7，最高 5.4 GHz）+ 16 E-core（8–23，最高 4.7 GHz），24 逻辑核 |
| GPU | NVIDIA GeForce RTX 5080 Laptop（sm_120 / Blackwell），16 GB，驱动 610.57.04 |
| 工具链 | gcc 16.2.1，nvcc 13.3.73（`/opt/cuda`），CMake 4.4.3 |
| 编译 | `-O3 -march=native`，Release，OpenMP 启用 |

---

## 0. 结论摘要

| # | 发现 | 类型 | 实测证据 | 预期收益 | 风险 |
|---|---|---|---|---|---|
| **F1** | 仿真 97% 的时间花在 coil–filament 的 4D 张量积求积上（CPU 与 GPU 同源） | 结构性 | 16 filament：2.19 ms/step，其中 16×140 µs = 2.24 ms 为求积 | — | — |
| **F2** | 求积节点数分配在错误的维度上：filament 自身两个轴用 9 点是巨大浪费 | **性能+精度双赢** | (11,11,3,3)=1089 点 vs 生产 (9,9,9,9)=6561 点，误差反而更小（2.42%→0.76% of peak） | **CPU 4–7×；GPU 5.5–6×；精度同时提升** | 中（需按几何自适应） |
| **F3** | `dM/dx` 的 4D 求积在工作区间**不收敛**，生产 9 点的误差达峰值的 2%–30% | **精度缺陷（阻塞 F2 的朴素版本）** | 与 M 的高阶有限差分对比：n=9 误差 4e-6…1.8e-5（峰 5.6e-5）；n=19 收敛 | 提升精度；是 F2/F4 的前提 | — |
| **F4** | 沿轨迹的 M(x) 未被制表；`OptimizationLevel::LookupTable` 是空实现 | **算法级** | 每径向层一条曲线可服务 m 个 filament；64 点表建表 1.5 ms，之后每步近零成本 | **CPU ~100×（长仿真）** | 中（需回退路径） |
| **F5** | 每次构造 DrivingCoil / Armature 若 (q,p) 出表，走 2048×16 复合 Bessel/Struve 积分 | 设置成本 | 13.5 ms/线圈、16.4 ms/电枢（n=1）；n≥2 时 0.005–0.009 ms | 优化循环每候选省 30–150 ms | 低 |
| **F6** | 并行只做在 filament 内层；候选/种群层完全串行；P/E 核放置敏感 | 并行 | 8 P-core 仅 2.06×，8 E-core 7.48×，16 E-core 10.08× | 2–2.6×（仅调线程放置） | 低 |
| **F7** | 求积规则只有 4/9/16/32，无法表达非对称阶数 | 基础设施 | `gauss_legendre(3)` 抛异常；GPU 侧硬编码 `n_nodes == 9` | 阻塞 F2 | 低 |
| **F8** | 若干次要项：每步堆分配、指标重复全历史扫描、CPU benchmark 未接入 CMake、RK4 事件二分代价、Full 截断用质心而非 filament 位置 | 工程债 | 见 §3.8 | 1–5% | 低 |
| **G1** | GPU 路径默认全程 FP64，而消费级 Blackwell 的 FP64 只有 FP32 的 1/64 | **最大 GPU 杠杆** | 同一内核 FP32 17.4 G 点/s vs FP64 1.45 G 点/s（**12×**）；端到端 batch≥32 快 **9.5×** | **9.5–12×** | 中（精度漂移需按仿真长度验证） |
| **G2** | 内核并行粒度 = 一个 block 负责一对 (sim,stage,filament)，小批量严重欠占用 | 结构 | batch=1、10 对时仅 10 个 block 跑在 84 个 SM 上；点级并行重构后 **0.287→0.081 ms（3.6×）** | 1.4–3.6×（batch=1） | 低 |
| **G3** | 设备端用 Boost `ellint_1`/`ellint_2`（两次独立调用）而非融合 AGM | 内核 | 融合 AGM 快 **1.2–1.4×**，与 Boost 一致到 **4e-12** | 1.2–1.4× | 极低 |
| **G4** | `BackendMode::Persistent` 不可达，且回退路径比直接用 CPU 慢 3.1× | 正确性/可用性 | `supports_persistent_control_stream=false` → 恒 fallback；fallback 9.19 ms/step vs CPU 2.93 ms/step | 避免 3× 性能陷阱 | 低 |
| **G5** | GPU 未接入优化器：`set_gpu_batch_evaluator` 仅在测试中以 mock 出现 | 集成缺口 | 全仓库仅 `tests/test_coilgun_optimization.cpp` 三处 mock 调用 | GA 评估 ~90× | 中 |
| **G6** | batch=1 时 GPU 仅与 CPU 相当（0.6–1.4 ms/step）；batch≥32 才进入吞吐区 | 定位 | 10 filament：0.607 ms/step（batch=1）vs 0.074 ms/sim（batch=128） | 明确分工 | — |

**一句话：** CPU 与 GPU 的瓶颈是同一个 4D 求积内核，而且**该内核既太慢、又在关键区间不准**。最有价值的三项是：把求积节点从 filament 轴挪到线圈轴（F2，双端 4–6× 且更准）、GPU 改用 FP32 积分（G1，9.5–12×）、以及把沿轨迹的重复求积改成一次制表 + 插值求导（F4，CPU ~100×）。三者互不冲突，且 GPU 端叠加后可达 **50–70×**。

---

## 1. 现有数值模型与计算结构

### 1.1 物理模型（`docs/NumericalModel.md`）

- 电流丝法（CFM）：电枢离散为 `m`（轴向）× `n`（径向）个圆环电流丝，状态向量维度 `S + F`（S = 级数，F = m·n）。
- 电路方程（式 3.20）：`[İ] = ([L] − [M_I])⁻¹([U_eff] + v[Ṁ_I][I] − [R(T)][I] − [M][I])`。
- 力（式 5.1）：`F = Σ_s Σ_f I_d,s · I_a,f · dM_{s,f}/dx`。
- 运动（式 5.2/5.3）与绝热温升（式 6.10）显式欧拉耦合；CPU 另有事件感知 RK4（式 5.4/5.5 与 §5.4）。
- 关键量：`M`、`dM/dx` 由 **4D Gauss–Legendre 张量积**（式 4.13/4.16）给出，每对 (线圈, 电流丝) 需 `n_nodes⁴` 次 filament-pair 核函数求值。

### 1.2 CPU 调用链（热路径）

~~~text
MultiStageSim::step()                                  [src/simulation/multi_stage_sim.cpp:728]
  └─ advance_euler / advance_rk4_*                      (Euler 1 次 / RK4 4 次)
       └─ evaluate_derivatives()                        [:232]
            ├─ derive_resistance()                      (热耦合时)
            ├─ #pragma omp parallel for  over filaments [:282]
            │    └─ mutual_detail::mutual_inductance_coil_pair(..., nodes, false)  [:288]
            │         └─ 4D 张量积 9⁴ = 6561 次
            │              └─ compute_filament_pair_with_sqrt_ab()
            │                   └─ elliptic_ke(m)       (AGM，约 7 次 sqrt)
            ├─ 装配 system_matrix / rhs
            ├─ Eigen::LDLT 求解
            └─ compute_force() / 温升导数
  └─ record_step()                                      [:606]
~~~

单阶段路径同构（`src/simulation/single_stage_sim.cpp:183–194`）。

### 1.3 GPU 调用链（热路径）

~~~text
GpuEngine::step()                                      [src/cuda/gpu_engine.cu]
  ├─ 上传运行时掩码/电压等标量                     (每步多次 cudaMemcpyAsync H2D)
  ├─ launch_device_step = lambda:
  │    ├─ launch_mutual_input_update       (separation 更新)
  │    ├─ launch_mutual_pipeline           [src/cuda/gpu_mutual_pipeline.cu:30]
  │    │    grid = (F, S, B)，block = 256，每线程串行 6561/256 ≈ 26 个求积点
  │    │    每个点调用 Boost ellint_1 + ellint_2（两次独立求值，FP64）
  │    │    末尾 shared-memory 树形归约（8 次 __syncthreads）
  │    ├─ launch_device_assembly           (用预上传的 filament-mutual 矩阵)
  │    ├─ solver: Eigen(host) 或 cuBLAS/cuSolver(batched)
  │    ├─ launch_state_update_masked
  │    └─ launch_compact_status
  ├─ cudaGraph capture/replay（BackendMode::Graph）
  ├─ context_->synchronize()
  └─ 每步 D2H 下载状态（mutual/gradient/current/velocity/position/掩码/触发…）
~~~

**关键差异：** CPU 用一次融合 AGM 同时得到 K 与 E（`elliptic_ke`，10.1 ns）；GPU 的 `elliptic.cuh` 直接调用 Boost 的 `ellint_1`/`ellint_2`，**两次独立求值且不共享迭代**。

### 1.4 成本模型

~~~text
每步成本 ≈ (活跃级数) × F × n_nodes⁴ × c_kernel × (Euler:1 / RK4:4)
c_kernel ≈ 21 ns（CPU，含 M 与 dM/dx 融合）
每步成本(GPU) ≈ B × S_active × F × n_nodes⁴ × c_gpu，  c_gpu ≈ 0.69 ns/点(FP64) / 0.057 ns/点(FP32)
~~~

实测：`n_nodes=9` 时 CPU **140 µs / (线圈,电流丝) 对**（6561 × 21 ns）。

---

## 2. 测量基线

### 2.1 CPU 基线

| workload | filament 数 | setup (ms) | ms/step (1 线程) | ms/step (最佳并行) |
|---|---|---|---|---|
| single-16（1 级，16×1） | 16 | 30.0 | 2.19 | 0.29（16 线程 / E-core） |
| multi-32（2 级，16×2） | 32 | 27.2 | 4.50 | — |
| multi-128（8 级，16×8） | 128 | 90.7 | 18.29 | 2.30（16 线程） |

> setup 成本的构成见 F5：single-16 的 30.0 ms = 1 个线圈（13.5 ms）+ 1 个电枢（16.4 ms）；multi-128 的 90.7 ms ≈ 8 个线圈 × 约 11 ms。**两者都出查表范围，全部走了精确 Bessel/Struve 积分。**

单线程下 **每步时间 ≈ F × 140 µs**，与 4D 求积的成本模型完全吻合 → 求积占 **≈97%**。

### 2.2 GPU 基线（`bench_gpu_engine`，全部 141 个 GPU 测试用例通过）

| workload | batch | 后端/精度 | ms/step | 每 sim-step | 同轮 CPU 参考 ms/step | CPU/GPU 加速 |
|---|---|---|---|---|---|---|
| small-single（10 fil） | 1 | graph/full | 0.59 | 0.59 ms | 0.949 | 1.6× |
| medium-multi（32 fil） | 1 | direct/full | 0.67 | 0.67 ms | 2.118 | 3.2× |
| large-single（128 fil） | 1 | direct/full | 1.43 | 1.43 ms | 2.599 | 1.8× |
| medium-multi（32 fil） | 128 | direct/full | 37.4 | 0.292 ms | 2.118 | 7.3× |
| small-single（10 fil） | 128 | graph/full | 9.52 † | 0.074 ms | 0.949 | 12.8× |
| medium-multi（32 fil） | 128 | graph/**aggressive** | **3.94 †** | **0.031 ms** | 2.118 | **~68×** |
| small-single（10 fil） | 128 | graph/**aggressive** | **0.75 †** | **0.0059 ms** | 0.949 | **~161×** |

> † 来自 `gpu_scale_driver`（同几何、300 步）；其余来自 `bench_gpu_engine`。**CPU 参考在同一台机器上不同轮次可相差约 1.9×**（large-single 1.39–2.60 ms/step），这正是 F6 所述的 P/E 核放置敏感问题——比较 CPU/GPU 时必须同轮取数。

- GPU FP64 吞吐随批量饱和于 **约 296 µs / sim-step**（batch ≥ 32）；FP32 饱和于 **约 31 µs / sim-step**。
- **结论：GPU 的价值在批量吞吐（优化循环），而不是单次仿真延迟。** batch=1 时 GPU 与 CPU 相当甚至更慢。

### 2.3 内核微基准（隔离测量）

| 项目 | 实测 |
|---|---|
| CPU 单次 filament 核（M + dM/dz 融合） | 26.8 ns |
| CPU `elliptic_ke`（AGM，K 与 E 一次算） | 10.1 ns |
| CPU 9⁴ 一对 | 140 µs |
| GPU 空 kernel 启动开销 | **1.66 µs** |
| GPU FP64 吞吐（9⁴，大批量） | 1.45 G 点/s |
| GPU FP32 吞吐（9⁴，大批量） | **17.4 G 点/s** |
| GPU 设备端 Boost `ellint_1+ellint_2` vs 融合 AGM | AGM 快 1.2–1.4×，差异 4e-12 |

---

## 3. CPU 候选项

### F1 — 4D 求积主导一切（结构性）

- **位置：** `src/simulation/multi_stage_sim.cpp:273–297`、`src/simulation/single_stage_sim.cpp:183–195`、`src/physics/mutual_inductance.cpp:169–245`。
- **现象：** 每个时间步、每级、每 filament 都重新做一次完整 4D 求积。
- **机制：** 每步 O(S_active · F · n⁴) 次椭圆积分；与时间步数线性相乘。
- **证据：** §2.1。
- **结论：** 任何有意义的优化都必须作用于这一层；其余优化加起来也只有个位数百分比。

### F2 — 节点数分配在错误的轴上（性能与精度双赢，CPU/GPU 通用）

- **位置：** `mutual_inductance_coil_pair(..., int n_nodes, ...)`：同一个 `n_nodes` 同时作用于 **线圈径向、线圈轴向、filament 径向、filament 轴向** 四个轴。GPU 侧 `gpu_mutual_pipeline.cu:157` 更进一步硬编码 `n_nodes == 9`。
- **机制：** filament 的截面（`δr × δl`）比线圈截面小 1–2 个数量级，被积函数在这两个方向上几乎线性，2–4 点 Gauss 已经足够；而误差完全由线圈轴主导。
- **证据**（误差 = |dM/dx − 收敛参考| / 峰值|dM/dx|，参考由高阶 M 的 5 点中心差分给出）：

| 几何 | (9,9,9,9) 6561 点 | (11,11,4,4) 1936 点 | (11,11,3,3) 1089 点 | (13,13,3,3) 1521 点 | (15,15,2,2) 900 点 |
|---|---|---|---|---|---|
| bore_contained（电枢在膛内） | 2.42% / 4.45% | 0.75% / 1.50% | **0.76% / 1.49%** | 0.14% / 0.41% | **0.17% / 0.29%** |
| test_fixture（现有测试夹具） | 8.70% / 26.99% | 9.32% / 31.46% | **3.42% / 16.96%** | 4.82% / 10.77% | 11.99% / 38.43% |
| overlap_bench（径向量叠、非物理） | 15.97% / 32.07% | **6.95% / 16.33%** | 60.9% / 130.9% | 8.13% / 18.18% | 7.97% / 13.14% |

（格式：均值 / 最大值，占峰值 dM/dx 的比例）

- **GPU 侧实测**（隔离内核，`gpu_order_probe`）：

| 形状 | 9⁴ | 11×11×3×3 | 加速 |
|---|---|---|---|
| B128_S2_F32（128×2×32 对） | 43.28 ms | 7.88 ms | **5.5×** |
| B1024_S2_F32 | 357.9 ms | 61.4 ms | **5.8×** |
| B1_S1_F10 | 0.350 ms | 0.066 ms | **5.3×** |

- **结论：** 把 filament 两个轴降到 2–4 点、线圈轴提到 11–13 点，可在 **4–7× 更少的核求值**下获得**相等或更好**的精度；GPU 因吞吐受限，收益与点数近似线性。
- **风险：** 最优组合依赖几何（filament 越"厚"/越贴近线圈，需要更多 filament 节点），应**自适应选择**而非写死；对超范围几何需保守回退。
- **验证：** 用 `opt_probe` 对每种生产几何做"误差-成本"扫描；数值回归用 `test_mutual_inductance`、`test_multi_stage_sim`、`test_coupled_integrator`、`test_gpu_vs_cpu_*`。

### F3 — dM/dx 的 4D 求积在工作区间不收敛（精度缺陷）

- **位置：** 梯度积分路径（`mutual_inductance.cpp:68–74`；GPU 同理）。
- **机制：** 当两个回路半径接近（`ra ≈ rb`）且轴向接近（`za ≈ zb`）时，`dM_loop/dz ~ −2μ₀a/|d|`，被积函数在 4D 域中沿一个 2 维面出现 `1/|d|` 型近奇异。张量积 Gauss 对此收敛很慢（非单调）。**这正是电枢位于线圈内部的整个工作区间。**
- **证据**（`conv_probe2`，真值 = 收敛 M 的 5 点中心差分）：

| 几何 | 分离量 | 真值 dM/dx | n=9 结果 | 绝对误差 | 占峰值 |
|---|---|---|---|---|---|
| bore_contained | −0.006 | 2.161e-6 | 1.764e-6 | 3.97e-7 | 3.9% |
| bore_contained | 0.0008 | −2.686e-7 | −5.278e-7 | 2.59e-7 | 2.5% |
| test_fixture | −0.002 | 2.544e-6 | 3.890e-6 | 1.35e-6 | 27.0% |
| overlap_bench | −0.002 | 3.534e-6 | 1.759e-5 | 1.41e-5 | 25.1% |
| 任意 | ±0.03（远场） | — | — | ~2e-9 | ~0.004% |

- **含义：**
  1. **不能**靠"降低阶数"换性能——在关键区间阶数本来就不够（这也是 F2 必须把节点挪到线圈轴而非简单减少的原因）。
  2. 力的零穿越位置附近，梯度本身接近 0，绝对误差 1e-6 量级会让力的符号/零点位置偏移。
  3. M 本身收敛得比 dM/dx 好得多（n=15 已到 ~1e-4 相对），且 **dM/dx 可由 M 的插值导数精确得到** → 这是 F4 的物理依据。
- **建议：** 属精度问题，需与物理负责人确认期望精度；建议至少建立工作区间梯度收敛性的回归基准（现有测试未覆盖）。

### F4 — 沿轨迹制表（M(x) + 插值导数）——最大 CPU 算法级收益

- **位置：** `OptimizationLevel::LookupTable`（`include/coilgun/simulation/multi_stage_sim.hpp:33`）当前是空实现（"same runtime path as Reference"）。
- **机制：**
  1. 同一径向层 `j` 的所有 `m` 个 filament 几何完全相同（相同的 `r_inner/r_outer/δl`），只有轴向偏移不同 → **一条 M(sep) 曲线可服务 m 个 filament**；
  2. 分离量 `sep = x_arm + rel_axial_i − coil_pos`，只需查表 + 插值；
  3. `dM/dx` 由插值多项式的解析导数给出（见 F3）。
- **成本测算：**
  - 建表（每级 × 每径向层）：`N_table × c_pair`。用 F2 的 (15,15,2,2)=900 点（约 24 µs/点）：64 点表 ≈ **1.5 ms**；
  - 之后每步求积成本 ≈ 0，只剩 `F` 次三阶插值（约 0.05 µs/filament）；
  - single-16：1.5 ms 一次性投入 vs 2.19 ms/步 → **不到 1 步回本**；1000 步仿真约 **100× 加速**。
- **精度证据**（`opt_probe`）：插值本身不是瓶颈——64 点表与 1024 点表误差相同（bore_contained 均为 2.1%），误差完全由建表所用的 M 求积精度决定。**用收敛的 M 规则建表即可。**
- **风险 / 回退：** 表域须覆盖实际位移区间（用 `TerminationPolicy::barrel_end_position` + 裕量），越界或提前终止时回退直接求积；几何随候选变化时表不可复用 → 用**几何为键**的缓存（只变电压/电容/触发时复用）。
- **验证：** 同 workload 下"表路径 vs Reference 路径"的力/速度/温度逐点对比。

### F5 — 构造期的精确 Bessel/Struve 积分（优化循环的主要隐性成本）

- **位置：** `src/physics/self_inductance.cpp:69–94`（`integrate_T_kernel`：2048 个子区间 × 16 点，提前退出阈值 `1e-18·total` 极紧）；触发条件是 `(q,p)` 超出查表范围 `q∈[0.05,4]`、`p∈[1.05,4]`。
- **证据：**

| 对象 | 耗时 |
|---|---|
| `Armature` m=16, n=1（p=5.0 出表） | **16.35 ms** |
| `Armature` m=16, n=2/4/8（在表内） | 0.005–0.009 ms |
| `DrivingCoil`（q=5.0 出表） | **13.53 ms** |
| M 矩阵 128×128（8128 对） | 0.27 ms |

- **机制：** 常见线圈几何（长度/inner_radius ≈ 5）恰好落在表外；而 `CoilgunOptimizationProblem::evaluate_cpu`（`src/optimization/coilgun_problem.cpp:212`）**每个候选都重新构造线圈与电枢**，于是每候选 13.5 ms × 级数 + 16 ms 被反复支付。
- **建议：** (a) 扩展 T(q,p) 表范围或加外推；(b) 按文档 §4.1 的原始建议改用 15 点 Gauss–Laguerre（`gauss_laguerre_cached(15)` 已存在但未被 `self_inductance` 使用）；(c) 在优化层缓存与几何相关的构造结果。

### F6 — 并行策略：内层并行已饱和，候选层完全串行

- **位置：** `multi_stage_sim.cpp:282`（`#pragma omp parallel for if (N_fil_ >= 8)`）；`src/optimization/genetic_optimizer.cpp:18–26`、`coilgun_problem.cpp:255–291`（逐候选串行）。
- **证据（真实椭圆核，16 任务 × 约 110–150 µs，200 个并行区）：**

| 线程放置 | 加速比 |
|---|---|
| 8 线程 @ P-core 0–7 | **2.06×** |
| 8 线程 @ E-core 8–15 | **7.48×** |
| 16 线程 @ 0–15 | 3.09× |
| 16 线程 @ 8–23 | **10.08×** |
| 24 线程（默认放置） | 3.49× |

单任务成本：P-core 106 µs，E-core 152 µs，但**并行下 P-core 集群吞吐严重塌陷**（推测为 AVX2 重负载下的频率/共享 L2 影响）。端到端最佳配置（16 线程 pin 到 E-core）比默认 24 线程快 **2.6×**。

- **结论：**
  1. **零代码改动**的调优项：`OMP_PROC_BIND=spread` / 显式亲和性 / 限制线程数，需以 benchmark 决定。
  2. 中长期应把并行粒度上移到 **候选/种群层**（见 G5：GPU 批量评估是更彻底的方案）。

### F7 — 求积规则种类受限

- **位置：** `src/physics/quadrature.cpp`：Gauss–Legendre 仅支持 n = 4/9/16/32，Gauss–Laguerre 仅 15/30。
- **证据：** `gauss_legendre(3)` 抛 `std::invalid_argument`；`gpu_mutual_pipeline.cu:157` 要求 `n_nodes == 9`。
- **影响：** F2 的"非对称阶数"无法表达；GPU 路径同样被锁死。
- **建议：** 增加通用 Gauss–Legendre 生成器（Newton 迭代 + 进程级缓存，本报告探针已验证可行），保留现有 4/9/16/32 快表以兼容。

### F8 — 次要项（工程债，合计约 1–5%）

| ID | 位置 | 现象 | 建议 |
|---|---|---|---|
| CPU-08 | `tests/bench_cpu_sim.cpp` | **未接入任何 CMakeLists**，正常构建无法产出该可执行文件 | 加入 `tests/CMakeLists.txt`（`EXCLUDE_FROM_ALL` 或独立 benchmark target） |
| CPU-09 | `multi_stage_sim.cpp:606` / `single_stage_sim.cpp:423` | 每步 `push_back` 一个含 2–3 个 `std::vector` 的 step，产生大量堆分配 | 预留容量 / 列式存储 / 可配置采样间隔 |
| CPU-10 | `coilgun_problem.cpp:216–237` | `metric_value` 被调用 6+ 次，每次全历史扫描 | 单次遍历汇总所有指标 |
| CPU-11 | `multi_stage_sim.cpp:562–599` | RK4 事件二分：每次事件最多 64 次迭代 × 4 次导数求值 | 降低最大迭代（dt/tol 只需约 14 次） |
| CPU-12 | `multi_stage_sim.cpp:196–199, 278–281` | `Full` 的截断用**电枢质心**与线圈中心比较，但求积针对每个 filament 的分离量 | 截断判据改为 filament 级 |
| CPU-13 | `mutual_inductance.cpp:57–58` + `elliptic.cpp:62` | 钳位重复且上界不一致（`1−1e-10` vs `1−1e-12`）；`sqrt(m)` 被重复计算 | 统一钳位、把 k 而非 m 传入 AGM（省一次 sqrt，约 14%） |
| CPU-14 | `mutual_inductance.cpp:43–47` | 进程级 LRU 缓存（4096 条）在生产热路径上始终 `use_cache=false` | 要么启用（需线程安全/按线程分片），要么移除以免误导 |

---

## 4. GPU 候选项

### G1 — 默认全程 FP64，而消费级 Blackwell 的 FP64 只有 1/64 速率（最大杠杆）

- **位置：** `include/coilgun/physics/elliptic.cuh`（`double` 版 K/E）、`include/coilgun/physics/mutual_inductance.cuh`（`..._pair_device`）；`gpu_mutual_pipeline.cu:106` 在非 Aggressive 模式下调用它们。
- **机制：** RTX 5080（GB203）的 FP64 FMA 速率约为 FP32 的 1/64，且 `sqrt.f64`/`div.f64` 是软件序列。内核每点约 7–9 次 `sqrt` 加数次除法，因此被 FP64 吞吐钉死。
- **证据（隔离内核，`gpu_prec_probe`）：**

| 形状 | FP64 AGM | FP32 AGM | 加速 |
|---|---|---|---|
| B1_S1_F10 | 0.264 ms | 0.0265 ms | **10.0×** |
| B1_S2_F32 | 0.528 ms | 0.0464 ms | **11.4×** |
| B1_S1_F128 | 0.845 ms | 0.0715 ms | **11.8×** |
| B128_S2_F32 | 37.54 ms | 3.16 ms | **11.9×** |
| B1024_S2_F32 | 296.0 ms | 24.72 ms | **12.0×** |

- **端到端证据（`GpuEngine`，batch=128，300 步）：**

| workload | FP64 ms/step | FP32 ms/step | 加速 | 末速度相对偏差 |
|---|---|---|---|---|
| 1 级 10 filament | 9.52 | 0.75 | **12.7×** | 0.084% |
| 2 级 32 filament | 38.05 | 4.02 | **9.5×** | 0.349% |

- **精度：** FP32 融合 AGM 与 FP64 在 8192 组随机 (ra, rb, d) 上最大相对偏差为 **M: 8.6e-4，dM/dz: 4.8e-4**；项目已为 Aggressive 定义 **1e-2** 的容差（`tests/gpu_numerical_tolerances.hpp:20`），因此**在既有契约内**。但注意：**偏差随步数累积**（50 步 0.02% → 300 步 0.35%），必须按目标仿真长度验证。
- **建议：** 把 FP32 从"可选 Aggressive"提升为**批量评估的默认**（配合每案例的收敛验证），并考虑混合精度（近奇异区用 FP64、远场用 FP32）。当前 `CoilgunOptimizationProblem` 默认 `OptimizationLevel::Full`，GPU 侧默认 `PrecisionMode::Full` → **默认走的正是最慢的 FP64 路径**。

### G2 — 并行粒度 = 一个 block 一对，小批量严重欠占用

- **位置：** `gpu_mutual_pipeline.cu:181–183`：`dim3 grid(filament_count, stage_count, batch_size)`。
- **机制：** batch=1、10 filament、1 级时只有 **10 个 block × 256 线程 = 2560 线程**，而 84 个 SM 可容纳约 17 万线程 → 占用率约 1.5%，实测吞吐 0.19 G 点/s（峰值的 13%）。
- **证据（`gpu_occupancy_probe`，点级并行 = 每对拆成多个 block + atomicAdd）：**

| 形状 | block/对（现状） | 点级并行（最佳） | 加速 |
|---|---|---|---|
| B1_S1_F10 | 0.287 ms（10 block） | 0.081 ms（80 block） | **3.6×** |
| B1_S2_F32 | 0.591 ms（64） | 0.345 ms（512） | **1.7×** |
| B1_S1_F128 | 0.911 ms（128） | 0.663 ms（1024） | **1.4×** |
| B8_S2_F32 | 2.530 ms | 2.511 ms | 1.0×（已达峰值） |
| B128_S2_F32 | 39.89 ms | 39.80 ms | 1.0×（chunk 过大时反而变慢） |

- **建议：** 按 `batch_size × stage × filament` 自适应选择"每对 block 数"：小 workload 拆到 4–8 块/对，大 workload 保持 1 块/对（避免 atomic 竞争）。
- **风险：** 低（结果不变，只是归约方式改变）；需保证确定性（atomicAdd 对 double 的加法顺序不定 → 若要求位级可复现，改用两级归约或 `atomicAdd` + 固定顺序的第二次归约）。

### G3 — 设备端未使用融合 AGM（用 Boost 两次独立调用）

- **位置：** `include/coilgun/physics/elliptic.cuh:28–39`：`elliptic_k`/`elliptic_e` 各调用一次 `boost::math::ellint_1/2`；`mutual_inductance.cuh:33–34` 分别调用。
- **证据：** 融合 AGM（与 CPU `elliptic_ke` 同构）比 Boost 路径快 **1.2–1.4×**（B1_S1_F10 0.512→0.362 ms；B128_S2_F32 53.2→45.9 ms），且两者在 4096 组随机样本上最大相对差异 **3.8e-12**（M）/ **9.5e-13**（dM/dz）。
- **建议：** 把 CPU 的 `elliptic_ke` 逻辑移到 `elliptic.cuh`（纯算术，可直接 `__device__`），设备端复用它。这是**低风险、必做**的一项。

### G4 — `BackendMode::Persistent` 不可达，且回退比直接用 CPU 慢 3.1×

- **位置：** `gpu_execution_config.hpp:96` `supports_persistent_control_stream = false`；`gpu_execution_config.hpp:147–156` → 恒 `Fallback`；`gpu_engine.cu:557–562` 抛异常说明需要专用控制流。
- **证据（`bench_gpu_engine`）：** `persistent` 请求 → `resolved backend=fallback`，9.19 ms/step；而纯 CPU 参考为 2.93 ms/step。**回退路径比不启用 GPU 慢 3.1×**（回退仍走 engine 的装配 + 每步状态快照复制，且 `execute_cpu_physical_pipeline` 每步复制 7 份状态）。
- **风险：** 用户设置 `backend=Auto/Persistent` 时可能悄悄拿到一个比 CPU 慢 3 倍的路径。
- **建议：** (a) 让 `Auto` 在 CUDA 不可用/工作负载过小时**直接选择纯 CPU 仿真器**而非 engine 的 fallback；(b) 修复 `supports_persistent_control_stream`（需要独立控制流）或明确标注为不支持；(c) 让 fallback 路径跳过为 GPU 准备的状态快照与拷贝。

### G5 — GPU 未接入优化器（集成缺口）

- **位置：** `include/coilgun/optimization/coilgun_problem.hpp:81–98`：`GpuBatchEvaluator` 是一个 `std::function` 回调；全仓库仅 `tests/test_coilgun_optimization.cpp` 的三处 mock 使用。
- **现状：** `evaluate_batch` 在无回调时逐候选串行调用 `evaluate_cpu`（`coilgun_problem.cpp:287–290`），即 **GA 的整个种群是串行求值的**，而每次求值内部又只并行到 filament 层（F6）。
- **建议：** 提供一个内置的 GPU 批量适配器：把种群映射为 `GpuEngine` 的 batch 维度（每个候选一个 batch 槽），用 FP32 + 批量路径求值。按 §2.2 的数据，128 个体 × 2000 步：CPU 约 750 s → GPU(FP32) 约 8 s（**~90×**）。
- **风险：** 中（需要处理几何随候选变化、约束指标映射、失败隔离；现有 `CachedBatchEvaluator` 的失败隔离语义要保留）。

### G6 — GPU 的定位：批量吞吐，不是单次延迟

- **证据：** batch=1 时 GPU 0.59–1.43 ms/step，与 CPU 相当（128 filament 时甚至略慢 0.97×）；batch≥32 后才进入 296 µs/sim-step 的吞吐区。
- **建议：** 在文档与 API 中明确：**交互式单次仿真用 CPU；参数扫描/优化用 GPU 批量**。`GpuExecutionPlanner` 已有的 `large_workload = batch_size >= 8 || dimension >= 128` 判据可以复用，但应把"batch=1 时 GPU 无收益"作为显式决策（当前 `Auto` 会为 dimension≥128 选 Graph，即使 batch=1）。

### G7 — 次要项

| ID | 位置 | 现象 | 建议 |
|---|---|---|---|
| GPU-07 | `gpu_engine.cu` 每步 H2D 上传掩码/电压 + D2H 下载状态 | batch=128 时 transfer 0.33 ms/step（约占 0.9%，FP32 下约 8%） | 保留必要项，其余改为设备端常驻；Graph 已把 transfer 降到 0.1 ms/step |
| GPU-08 | `gpu_engine.hpp:262` 构造期用 `mutual_inductance_filament(..., true)` | 该 LRU 是进程级非线程安全缓存；多线程并行构造 engine 存在数据竞争 | 构造期改用 `use_cache=false`，或加锁/分片 |
| GPU-09 | `gpu_mutual_pipeline.cu:78` `__shared__ double sum_m[512]` | 固定 8 KB×2，与 `threads_per_block` 无关 | 按实际 blockDim 分配（占用率可提升） |
| GPU-10 | `tests/bench_gpu_engine.cu` 未覆盖 `PrecisionMode::Aggressive` 与 `GpuOptLevel::Standard` | 最大性能杠杆没有基准覆盖 | 把 Aggressive/Standard 加入固定 request 列表 |

---

## 5. 联合收益与建议实施顺序

按"收益/风险"排序，每阶段独立可回滚：

| 阶段 | 内容 | 适用 | 预期收益 | 风险 |
|---|---|---|---|---|
| 0 | 把 `bench_cpu_sim` 接入构建；GPU benchmark 增加 Aggressive/Standard；固定基线 | 双端 | 可观测性 | 无 |
| 1 | 线程放置与线程数调优（F6-1） | CPU | 2–2.6× | 低 |
| 2 | 设备端融合 AGM（G3）+ 通用 Gauss–Legendre 生成器（F7） | 双端 | 1.2–1.4×（GPU） | 极低 |
| 3 | **按轴独立阶数（F2）**，用误差-成本扫描确定默认与自适应规则 | 双端 | CPU 4–7×；GPU 5.5–6×；**精度同时提升** | 中 |
| 4 | **GPU FP32 批量默认（G1）**，按仿真长度验证漂移 | GPU | 9.5–12× | 中 |
| 5 | 点级并行自适应分块（G2） | GPU batch=1 | 1.4–3.6× | 低 |
| 6 | 轨迹制表 + 解析导数（F4） | CPU | ~100×（长仿真） | 中 |
| 7 | 消除构造期精确积分（F5） | CPU/优化 | 每候选 30–150 ms | 低 |
| 8 | GPU 批量适配器接入优化器（G5）；修复 fallback 语义（G4） | 优化层 | GA ~90× | 中 |
| 9 | F8 / G7 各项清理 | 双端 | 1–5% | 低 |

**组合潜力（推测，基于已测单项）：**
- GPU 批量：FP32（9.5×）× 非对称阶数（5.5×）≈ **52×**，叠加 G3（1.3×）≈ **68×**。
- GPU 单次（batch=1）：点级并行（3.6×）× FP32（10×）× 非对称阶数（5×）≈ **180×** → 单次仿真每步从 0.6 ms 降到约 3 µs。
- CPU：非对称阶数（6×）× 制表（~100×，长仿真）——但两者部分重叠，实际取制表为主。

**不建议**（记录以免重复讨论）：
- 仅靠降低 `n_nodes` 换速度——F3 表明工作区间 9 点尚不足。
- 依赖 `BackendMode::Persistent`——当前恒回退（G4）。
- 在 batch=1 时为了 GPU 而 GPU——G6 显示无收益。

---

## 6. 未验证 / 待确认

- **F3 是否算缺陷：** 需物理/项目负责人确认期望精度；本报告只给出量化误差。
- **G1 的 FP32 精度：** 实测漂移随步数累积（0.35% @ 300 步），必须按目标仿真长度与目标速度量级重新验证；本报告的算例速度极低（0.1 m/s 量级），不代表高能工况。
- **GPU 数值一致性：** 141 个 GPU 测试用例全部通过（`test_gpu_vs_cpu_single` 31、`test_gpu_vs_cpu_multi` 36、`test_gpu_solver` 21、`test_gpu_sim_batch` 17、`test_gpu_paths` 12、`test_gpu_graph` 11、`test_gpu_precision` 5、`test_gpu_thermal` 5、`test_gpu_assembly` 2、`test_gpu_batch` 1），但未覆盖 Aggressive 精度的端到端一致性。
- **CUDA 图捕获：** 首次捕获 6–16 ms（`cold-first-step`），大批量下摊薄后收益有限（Graph 与 Direct 的 steady-state 接近）。
- **AGENTS.md 已过期：** 文档仍写"Python 原型已删除、C++ 物理层完成、18 个测试套件"，但当前树已包含完整 CUDA 后端与优化算法层。

---

## 7. 复现方式

所有探针位于 `build/perf-scratch/`（构建目录，不入库），源码一并保留。

```sh
# ---- CPU ----
cmake --build build/cpu-release --target coilgun -j
EIGEN=build/ninja-debug/_deps/eigen-src
CXX="g++ -O3 -march=native -fopenmp -std=c++20 -DNDEBUG -I include -isystem $EIGEN"
LIB=build/cpu-release/src/libcoilgun.a
$CXX build/perf-scratch/e2e_bench.cpp      $LIB -o build/perf-scratch/e2e_bench      && ./build/perf-scratch/e2e_bench 1 16 1 200
$CXX build/perf-scratch/micro_perf.cpp     $LIB -o build/perf-scratch/micro_perf     && ./build/perf-scratch/micro_perf
$CXX build/perf-scratch/conv_probe2.cpp    $LIB -o build/perf-scratch/conv_probe2    && ./build/perf-scratch/conv_probe2
$CXX build/perf-scratch/opt_probe.cpp      $LIB -o build/perf-scratch/opt_probe      && ./build/perf-scratch/opt_probe
$CXX build/perf-scratch/setup_probe.cpp    $LIB -o build/perf-scratch/setup_probe    && ./build/perf-scratch/setup_probe
taskset -c 8-23 env OMP_NUM_THREADS=16 ./build/perf-scratch/e2e_bench 1 16 1 300   # 线程放置对比

# ---- GPU ----
export PATH=/opt/cuda/bin:$PATH
cmake --build build/cuda-release --target coilgun_cuda bench_gpu_engine -j
./build/cuda-release/src/cuda/bench_gpu_engine > build/perf-scratch/gpu_bench_raw.txt

NVCC="nvcc -O3 -arch=sm_120 -std=c++20 --expt-relaxed-constexpr -Xcompiler=-march=native -I include"
# 隔离内核：Boost vs 融合 AGM、批量扫描
$NVCC build/perf-scratch/gpu_kernel_probe.cu   build/cuda-release/src/cuda/libcoilgun_cuda.a build/cuda-release/src/libcoilgun.a -lcublas -lcusolver -o build/perf-scratch/gpu_kernel_probe   && ./build/perf-scratch/gpu_kernel_probe
# 非对称阶数
$NVCC build/perf-scratch/gpu_order_probe.cu -o build/perf-scratch/gpu_order_probe && ./build/perf-scratch/gpu_order_probe
# FP64 vs FP32 吞吐与精度
$NVCC build/perf-scratch/gpu_prec_probe.cu  -o build/perf-scratch/gpu_prec_probe  && ./build/perf-scratch/gpu_prec_probe
# 点级并行 / 占用率
$NVCC build/perf-scratch/gpu_occupancy_probe.cu -o build/perf-scratch/gpu_occupancy_probe && ./build/perf-scratch/gpu_occupancy_probe
# 引擎端到端（需 -DCOILGUN_CUDA_AVAILABLE=1，否则会静默回退到 CPU）
$NVCC -DCOILGUN_CUDA_AVAILABLE=1 -isystem $EIGEN build/perf-scratch/gpu_scale_driver.cu       build/cuda-release/src/cuda/libcoilgun_cuda.a build/cuda-release/src/libcoilgun.a -lcublas -lcusolver       -o build/perf-scratch/gpu_scale_driver && ./build/perf-scratch/gpu_scale_driver 2 8 4

# GPU 测试套件
cd build/cuda-release/tests && for t in test_gpu_vs_cpu_single test_gpu_vs_cpu_multi test_gpu_precision   test_gpu_solver test_gpu_batch test_gpu_sim_batch test_gpu_paths test_gpu_thermal test_gpu_assembly test_gpu_graph; do ./$t; done
```

> 注意：直接编译的 GPU 驱动若缺少 `-DCOILGUN_CUDA_AVAILABLE=1`，`GpuEngine` 会**静默回退到 CPU**（`fallback_reason = "CUDA support is not compiled"`），且该回退路径比纯 CPU 慢约 3×（见 G4）。
