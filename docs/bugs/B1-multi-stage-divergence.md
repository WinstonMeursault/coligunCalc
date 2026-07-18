# B1: Multi-Stage GPU Simulation — Numerical Divergence Bug

> 日期：2026-07-17 | 状态：已解决 | 严重度：高

## 症状

多级 GPU 仿真出口速度 ~0.001 m/s，CPU 为 ~1.6 m/s（差 3 个数量级）。单级 GPU 仿真正常（0.2% 误差）。

## 定位过程

### 步骤 1：确认 CPU 一致性

单级仿真 GPU vs CPU：v_muzzle 0.947 → 0.949 m/s（0.2% 误差）。✅ kernel 基本正确。

### 步骤 2：找分歧点

50 步诊断输出：

```
step0: F_cpu=0.03127 F_gpu=0.03127 rel_diff=3.3e-14 I_cpu=2.20 I_gpu=2.20  ← 一致
step1: F_cpu=0.12404 F_gpu=0.12404 rel_diff=3.3e-14 I_cpu=4.39 I_gpu=4.39  ← 一致
step2: F_cpu=0.27681 F_gpu=0.27681 rel_diff=3.3e-14 I_cpu=6.55 I_gpu=6.55  ← 一致
step3: F_cpu=0.48808 F_gpu=0.00348 rel_diff=0.993   I_cpu=8.71 I_gpu=7.80  ← 分叉！
```

步骤 0-2 **位级一致**（rel_diff ~1e-14，机器精度极限）。步骤 3 突然分叉，GPU 力量仅为 CPU 的 0.7%。

### 步骤 3：状态完整性验证

步骤 2 结束状态逐字节对比——全部 12 个电流 + 位置 + 速度：

```
CPU step2: pos=0.05 vel=1.29424e-06 I[0]6.55497 ... I[11]-3.02521
GPU step2: pos=0.05 vel=1.29424e-06 I[0]6.55497 ... I[11]-3.02521
```

每一个 double 字段都完全一致。步骤 3 的**输入状态完全相同**。

### 步骤 4：Kernel 确定性验证

独立测试：调用 CPU 版 `mutual_inductance_coil()` 与 GPU kernel 对比，连同参数的独立 Armature 对象计算 vs GPU sim 内部存储：

| 调用路径 | M(0,0) | dM(0,0) |
|---|---|---|
| 独立 CPU 函数 | 7.88e-7 | -3.19e-5 |
| GPU sim 内部 kernel | 7.82e-7 | -3.19e-5 |
| 独立+GPU 参数调用 | 7.88e-7 | -3.19e-5 |

kernel 输出值正确，确定性无漂移。✅

### 步骤 5：Kernel → CPU 函数交换测试

将 GPU `compute_M1_dM1()` 中的 kernel launch 替换为直接调用 CPU `mutual_inductance_coil()`，绕过 GPU 路径：

```
Before (GPU kernel):   force_gpu=0.00348, I0=7.80
After (CPU function):  force_gpu=0.36199, I0=8.68
                           ↑ (改善但仍偏离 CPU 的 0.488)
```

偏差从 99.3% 降到 26%，但依然存在。结论：**kernel 不是根源**。

### 步骤 6：顺序验证

将 GPU 的 `compute_derivatives` 中 `compute_M1_dM1()` 和 `compute_force()` 的顺序调整为与 CPU 一致（先 force 后 M1）：

```cpp
// Before (GPU): compute_M1_dM1() → compute_force()
// After (GPU):  compute_force() → compute_M1_dM1()   ← 与 CPU 一致
```

步骤 3 力量仍然分叉。顺序不是根源。

## 确认排除的假设

| 假设 | 验证方法 | 结论 |
|---|---|---|
| GPU kernel 计算错误 | 独立 CPU 函数对比 | ✅ 排除 |
| 输入状态不同 | 逐字节对比 | ✅ 排除 |
| ODR 违规（Eigen 对齐） | `EIGEN_MAX_ALIGN_BYTES=32` + `-march=native` | ✅ 已修复 |
| CPU/GPU 代码逻辑不一致 | 逐行对比 `build_system_matrix`、RHS、LDLT | ✅ 排除 |
| 设备缓冲区垃圾值 | `cudaMemset` 清零 | ✅ 已修复 |
| 操作顺序 | 调整后无变化 | ✅ 排除 |
| 缓存污染 | `use_cache=false` | ✅ 排除 |

## 当前结论

问题的根源在**更深层**。所有逻辑路径均正确、输入状态一致、kernel 输出正确的条件下，LDLT 求解仍产生不同结果。高度怀疑：

1. **nvcc vs g++ 对 Eigen 模板的编译差异**：`MultiStageState` 的 `operator+=`/`operator*` 在 `libcoilgun.a`（g++）和 lambda 中通过 `EulerStepper::advance` 模板实例化时产生不同的 inlining/SSE 路径
2. **Eigen LDLT 内部矩阵操作的浮点顺序**：相同输入、不同编译器优化导致 LDLT 分解在极远的小数位上分叉，经多步累积放大

## 建议下一步

1. **`compute-sanitizer` + 内存检查**：排除 nvcc 编译的 host 代码中隐藏的内存越界
2. **单步 LDLT 矩阵导出**：在步骤 3 前将 `L_total_` 和 `RHS_` 保存到文件，用 Python/numpy 独立 LDLT 求解对比
3. **最小可复现 case**：用 1 级（n_stages=1）运行 GpuMultiStageSim，与 CPU 版对比——这能排除多级逻辑，聚焦 nvcc 编译的 Eigen 模板问题

## 已修复的相关问题

- **ODR 违规**：`EIGEN_MAX_ALIGN_BYTES=32` + `-Xcompiler=-march=native` 确保 Eigen 对齐在 nvcc/g++ 间一致
- **n_launched=0 保护**：`compute_M1_dM1` 已有早期返回
- **`.cpp` → `.cu` 编译分离**：`add_gpu_test` vs `add_gpu_test_cuda` 避免 nvcc 编译纯 CPU 代码

## 根本原因（2026-07-18 确认）

通过详细的 LD/LT/RHS 诊断输出，发现分歧的直接原因是**N=4 的 GL 积分在 GPU kernel 中数值精度不足**。

### 触发条件

`compute_M1_dM1()` 中的优化逻辑：

```cpp
int n_nodes = (opt_level_ == GpuOptLevel::Full && dist > coil.length()) ? 4 : 9;
```

当电枢中心移动到线圈边缘之外时（`dist > coil.length()`），`n_nodes` 从 9 切换到 4。对于本测试用例：
- 线圈长度 = 0.05 m，线圈中心位置 = 0.0
- 电枢初始位置 = 0.05 m（正好在线圈边缘）
- Step 0-2: `dist = 0.05`，`n_nodes = 9`（`dist` 不大于 `coil.length()`）
- Step 3: `pos = 0.0500000000002606`（移动了 2.6e-13 m），`dist = 0.0500000000002606 > 0.05`，`n_nodes = 4`

### 影响

4 节点 GL 积分（4^4 = 256 个求值点）vs 9 节点（9^4 = 6561 个求值点），在线圈边缘附近（积分核变化剧烈），GPU kernel 的精度不足，导致 M1 值剧烈跳变：

| Step | n_nodes | M1_mat(0,0) |
|------|---------|-------------|
| 0-2  | 9       | 7.88e-7     |
| 3    | 4       | 3.24e-9     |

M1 值下降 240 倍（从 7.88e-7 到 3.24e-9），尽管位置仅移动了 2.6e-13 m。这在物理上是不可能的——ΔM ≈ dM/dz × Δz ≈ (-3e-5) × 2.6e-13 ≈ -7.8e-18，与观测的 7.8e-7 不符。

错误的 M1 值导致 LDLT 求解产生错误的电流，电流又影响力和下一时间步的状态，形成正反馈发散。

### 为什么 CPU 路径不受影响

CPU 路径使用 `physics::mutual_inductance_coil()`（基于 Boost.Math 的椭圆积分），其在低 GL 节点数下的积分精度优于 GPU kernel 的实现（GPU kernel 使用共享内存树归约，累加顺序不同，对边界附近的被积函数变化更敏感）。

### 修复

移除 `n_nodes = 4` 优化，始终使用 `n_nodes = 9`：

```cpp
// src/cuda/gpu_multi_stage_sim.cu:194
int n_nodes = 9;

// src/cuda/gpu_single_stage_sim.cu:125
int n_nodes = 9;
```

### 验证

修复后 CPU 和 GPU 多级仿真结果一致：

```
v_cpu = 1.637651752 m/s
v_gpu = 1.637651752 m/s
rel_err = 3.69e-11  （机器精度）
```

测试 `test_gpu_vs_cpu_multi` 现在完整验证 `v_gpu ≈ v_cpu`，epsilon = 5e-3。
