# 统一 GPU 引擎 Benchmark

本文档记录统一 CUDA 引擎的实测结果。计时数据仅代表特定机器上的观测值，
不构成性能保证。

## 方法

- CPU 行使用 Reference 优化级别的 `MultiStageSim<EulerStepper>`。
- GPU 行使用 `SimBatch<EulerStepper>`，并使用相同的几何参数、激励值、时间步长
  和固定步数。
- CPU/GPU 加速比为 `CPU 墙钟时间 / GPU 墙钟时间`。对于 batch 行，CPU 基线乘以
  batch 大小，因为 CPU Reference 是以独立模拟的方式运行的。
- GPU 计时字段来自 `ExecutionReport` 的累计主机墙钟时间分类：`gpu_time_ms`、
  `solver_time_ms`、`thermal_time_ms` 和 `transfer_time_ms`。
- `Graph` 只表示当前由图辅助的互感计算 segment。
- 在 resident control-stream backend 实现之前，`Persistent` 预期会解析为文档中
  说明的安全 fallback。
- 空行或 fallback 行本身仍是有效观测值，但不得将其视为成功的 GPU 计时。

## 复现

```sh
cmake --build --preset ninja-cuda-debug --target bench_gpu_engine
./build/ninja-cuda-debug/src/cuda/bench_gpu_engine
```

## 测量记录

下表来自目标机器的 benchmark 输出。对于不同硬件或构建配置，应追加带日期的
小节，而不是覆盖既有观测记录。

### RTX 5080 Laptop，CUDA Compute Capability 12.0

测量日期为 2026-07-21，使用 CUDA Debug 构建。墙钟时间包含 wrapper 编排开销；
GPU 计时字段是实际 CUDA 路径的报告值。CPU 行为独立运行的 Reference 模拟。

| 工作负载 | 请求 | Batch | 解析后端 | 求解器 | 热模式 | CPU/GPU 加速比 | GPU 执行 | 墙钟 ms | GPU ms | 求解器 ms | 热路径 ms | 传输 ms | Graph 重建次数 |
|---|---|---:|---|---|---|---:|---|---:|---:|---:|---:|---:|---:|
| small-single | CPU | 1 | cpu-reference | eigen | cpu | 1.000 | no | 114.944 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| small-single | Direct | 1 | direct | eigen | disabled | 11.658 | yes | 9.860 | 9.772 | 0.355 | 0.000 | 1.037 | 0 |
| small-single | Graph | 1 | graph | eigen | disabled | 12.287 | yes | 9.355 | 9.284 | 0.334 | 0.000 | 1.231 | 1 |
| small-single | Persistent | 1 | fallback | eigen | disabled | 0.161 | no | 714.376 | 0.000 | 0.376 | 0.000 | 0.000 | 0 |
| small-single | Fallback | 1 | fallback | eigen | disabled | 0.160 | no | 719.898 | 0.000 | 0.359 | 0.000 | 0.000 | 0 |
| medium-multi | CPU | 1 | cpu-reference | eigen | cpu | 1.000 | no | 254.073 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| medium-multi | Direct | 1 | direct | eigen | disabled | 26.250 | yes | 9.679 | 9.572 | 1.854 | 0.000 | 1.357 | 0 |
| medium-multi | Graph | 1 | graph | eigen | disabled | 24.863 | yes | 10.219 | 10.124 | 1.803 | 0.000 | 1.295 | 1 |
| large-single | CPU | 1 | cpu-reference | eigen | cpu | 1.000 | no | 231.163 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| large-single | Direct | 1 | direct | cusolver | disabled | 2.105 | yes | 109.811 | 109.637 | 103.602 | 0.000 | 0.517 | 0 |
| large-single | Graph | 1 | graph | cusolver | disabled | 20.626 | yes | 11.208 | 11.140 | 4.082 | 0.000 | 0.511 | 1 |
| batch-medium | Direct | 1 | direct | eigen | disabled | 26.060 | yes | 9.870 | 9.760 | 1.787 | 0.000 | 1.375 | 0 |
| batch-medium | Direct | 8 | direct | cusolver | disabled | 63.514 | yes | 32.396 | 32.033 | 2.706 | 0.000 | 1.117 | 0 |
| batch-medium | Direct | 32 | direct | cusolver | disabled | 71.586 | yes | 114.972 | 113.755 | 5.946 | 0.000 | 1.827 | 0 |
| batch-medium | Direct | 128 | direct | cusolver | disabled | 76.152 | yes | 432.314 | 426.954 | 10.630 | 0.000 | 0.986 | 0 |
| thermal-engine | CPU thermal | 1 | cpu-reference | eigen | cpu | 1.000 | no | 169.713 | 0.000 | 0.000 | 0.000 | 0.000 | 0 |
| thermal-engine | GPU Full thermal | 1 | direct | eigen | gpu | 12.676 | yes | 13.388 | 13.338 | 1.552 | 0.508 | 0.725 | 0 |
| thermal-engine | GPU Aggressive thermal | 1 | direct | eigen | gpu | 41.896 | yes | 4.051 | 4.011 | 1.511 | 0.463 | 1.114 | 0 |

本次运行还测量了 `medium-multi` 的 Persistent 和 Fallback 请求，其加速比分别为
`0.060x` 和 `0.060x`；两者都解析为 CPU fallback。Persistent 的 fallback 原因是
`CUDA runtime resource initialization failed: persistent protocol requires a dedicated control stream`，
Fallback 的原因是 `CPU fallback explicitly requested`。大型 workload 和 thermal
wrapper 场景也观察到了相同的行为。这些是 fallback 计时，不是 GPU 性能结果。

该目标会输出以下行：

- CPU Reference 与 Direct、Graph、Persistent 和 Fallback 请求的对比。
- 小型单级、中型多级、大型高分辨率单级以及 thermal 单级 workload。
- Batch 大小为 1、8、32 和 128 的 Direct 路径。
- 通过 raw engine thermal contract 解析 CPU thermal、GPU Full thermal 和 GPU
  Aggressive thermal 的字段。`SimBatch` 当前没有公开的 thermal-mode 参数，因此
  wrapper benchmark 行默认关闭 thermal；只有 raw engine 路径可以启用该模式。

## 精度误差记录

F1 精度 fixture 使用设计容差：Standard `{relative=5e-6, absolute=1e-10}`、
Full `{1e-4, 1e-9}` 和 Aggressive `{1e-2, 1e-8}`。在同一次 RTX 5080 Laptop
运行中，20 步单级 CPU 对比得到的最大相对误差如下：

| 精度 | M/dM fixture | 线圈电流 | 细丝电流 | 位置 | 速度 | 力 |
|---|---|---:|---:|---:|---:|---:|
| Standard | 通过；低于容差 | 2.13e-14 | 4.33e-13 | 1.39e-16 | 2.41e-13 | 2.96e-13 |
| Full | 通过；低于容差 | 2.13e-14 | 4.33e-13 | 1.39e-16 | 2.41e-13 | 2.96e-13 |
| Aggressive | 通过；低于容差 | 1.69e-05 | 3.73e-04 | 2.12e-11 | 4.91e-05 | 1.25e-04 |

这些是测试 fixture 的观测值，不是适用于所有情况的误差上界。测试还会将最终
出口速度、峰值力、温度以及由温度推导出的细丝电阻与 CPU Reference 进行比较。

## 求解器决策记录

| Workload 形状 | `S+F` / Batch | 解析求解器 | 依据 |
|---|---:|---|---|
| 小型单级 | 11 / 1 | Eigen | `small-single` 解析为 Eigen；Direct 墙钟时间 9.860 ms |
| 中型多级 | 34 / 1 | Eigen | `medium-multi` 解析为 Eigen；Direct 墙钟时间 9.679 ms |
| 大型单级 | 129 / 1 | cuSOLVER batched | `large-single` 解析为 cuSOLVER；Graph 墙钟时间 11.208 ms |
| 中型 batch | 34 / 8、32、128 | cuSOLVER batched | 所有大于等于 8 的 batch 大小都解析为 cuSOLVER |

静态 planner 规则保持为 `batch_size >= 8 || S+F >= 128`；本次没有根据单台机器
的计时结果调整阈值。结构化/block solver direction 3，包括 Schur complement、
block elimination 和 cuSPARSE 评估，仍属于后续工作，本阶段未实现。

## Planner 阈值

首次测量期间，当前静态 planner 阈值保持不变：

- 大型 workload：`batch_size >= 8` 或 `stage_count + filament_count >= 128`。
- 在 capability 可用时，Auto solver 仅对大型 workload 选择 batched/cuSOLVER。
- 在 capability 可用时，Auto GPU thermal 仅对大型 workload 选择 GPU thermal。

只有在代表性 workload 的实测结果证明收益稳定后，才能调整这些阈值。任何调整都
必须记录调整前后的 benchmark 行，并且不得将结果描述为通用保证。
