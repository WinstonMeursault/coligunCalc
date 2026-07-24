# coligunCalc

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language](https://img.shields.io/badge/language-C%2B%2B20-blue)](https://isocpp.org/)
[![Status](https://img.shields.io/badge/status-developing-orange)]()

[English Documentation](README.md)

---

基于**电流丝法**（CFM）的多级同步感应线圈炮仿真库。通过对完整椭圆积分的数值积分，对驱动线圈与电枢之间的电磁耦合进行建模，预测洛伦兹力、速度和效率。

采用 GPLv3 许可证。

---

## 主要特性

- **电流丝法** — 电枢离散为 m × n 个同心环形电流丝，捕捉趋肤效应
- **丝级解析核** — 互感通过完全椭圆积分 K(k) / E(k) 的解析表达式（Maxwell 公式）
- **4D Gauss-Legendre 求积** — 线圈对线圈互感，支持可配置求积点数（n⁴ 次核计算）
- **自感双轨计算** — 快速 T(q,p) 查表（~μs），自动回落至参考级 Bessel/Struve 积分（~ms）
- **Struve H₀/H₁ 实现** — SciPy 等效三区间策略（幂级数 / 渐近展开 / 交叉验证）
- **续流（Crowbar）二极管** — 每级自主管理二极管状态，消除制动力
- **可配置激励源** — 电容放电、含续流二极管、或任意波形 V(t)
- **两种时间步进器** — 前向 Euler（快速）和经典 RK4（高精度）
- **可选热耦合** — 绝热焦耳加热，含铜铝的温度相关 cp(T) 和 ρ(T)
- **三级优化等级** — 参考 / 查表 / 全优化（距离截断 + 自适应 GL 阶数）
- **级间触发** — 位置触发或时间延迟触发，最多 50 级
- **LRU 缓存** — 丝级 M 和 dM/dz 均有 4096 条缓存
- **OpenMP 多核支持** — 当丝元数至少为 8 时，CPU 的 M/dM 循环会并行化；力求和和热更新循环当前仍为串行
- **公开状态与诊断 API** — `FilamentMetadata`、激励快照、`IntegrationState`、可复用的 `DerivativeWorkspace`，以及可选的 CPU 阶段计时
- **统一 GPU 执行 API** — 显式后端/求解器/精度/热耦合策略、`ExecutionReport`、行优先 `GpuStateLayout`、Graph 边界掩码和稳定的批次行诊断
- **GPU 加速（CUDA）** — 同步 Euler 方法接口与 CPU API 对齐；GPU 单级/多级 RK4 实例当前在步进时会抛出 `std::logic_error`，需使用 CPU 实现

---

## 项目结构

```
include/coilgun/   — 公开头文件（核心类型、物理层、组件层、仿真器）
src/               — 库实现（静态库 libcoilgun.a）
tests/             — 单元测试与集成测试（doctest/CTest）
tools/             — T(q,p) 查表生成工具
docs/              — 文档（API 参考中/英文、数值模型、review/性能报告、benchmark 记录）
.references/       — 参考论文 PDF（gitignored，仅本地）
```

## 快速开始

**依赖**：C++20 编译器（支持 OpenMP）、CMake ≥ 3.20。

Boost.Math 和 Eigen 均通过 CMake `FetchContent` 自动获取，无需手动安装。

```sh
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug
```

### GPU 加速（CUDA）

**附加依赖**：CUDA Toolkit ≥ 12.8（Blackwell sm_120 需要；早期架构 ≥ 9.0），NVIDIA GPU 计算能力 ≥ 6.0。

```sh
cmake --preset cuda-debug
cmake --build --preset cuda-debug
ctest --preset cuda-debug
ctest --preset cuda-debug -L gpu
```

仅用于测量的 CUDA benchmark 同时包含 CPU Reference 行，并记录 Direct、
Graph、Persistent 请求、Fallback、批量大小、求解器和热路径计时。由于结果
依赖机器，它不加入默认构建：

```sh
cmake --build --preset cuda-debug --target bench_gpu_engine
./build/cuda-debug/src/cuda/bench_gpu_engine
```

RTX 5080 Laptop 的测量结果以及 Persistent 诚实回退状态见
[benchmark 记录](docs/benchmarks/2026-07-19-unified-gpu-engine.md)。

同步 GPU 方法接口在 Euler 下与 CPU API 对齐。GPU 单级/多级 RK4 实例当前在步进时会抛出 `std::logic_error`；RK4 请使用 CPU 实现。

```cpp
#include <coilgun/coilgun_cuda.hpp>

using coilgun::simulation::cuda::GpuSingleStageSim;
// 或: coilgun::simulation::cuda::GpuMultiStageSim

GpuSingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), dt);
sim.run();
double v = sim.result().summary.muzzle_velocity;
```

CUDA 引擎将固定形状 Euler 状态、矩阵/RHS 组装、批量求解、力/运动更新、可选热更新和紧凑设备控制保留在 resident 缓冲区中。`Graph` 捕获并重放完整的固定形状设备步；`Direct` 直接发射同一组 resident 阶段。运行时初始化、分配、求解、捕获、重放或校验失败时，引擎恢复步前状态，记录 `FallbackReason::RuntimeFailure`，并锁定至 CPU/Eigen 回退。`Persistent` 在具备独立 resident 控制 stream 前仍明确使用安全回退。

后端选择可通过 `ExecutionReport` 显式观察。`Fallback` 仅使用 CPU，且不会创建 CUDA context；运行时可用时，`Graph` 和 `Direct` 是支持的 resident CUDA 路径。请求 GPU 路径不等于实际执行 GPU：应检查 `gpu_executed` 和解析后的后端。`SimBatch` 保持稳定的物理行，并提供 `active_row_count()` 和 `mutual_gradients()` 用于批次诊断。CUDA wrapper 仍仅支持 Euler；RK4 由 CPU 仿真器支持。

没有 CUDA 设备时，GPU 测试可以 skip。真实 GPU validation 还必须要求 `ExecutionReport::gpu_executed == true` 且 `ExecutionReport::backend != BackendMode::Fallback`；请求了某个后端或 CPU 回退测试通过都不能证明 GPU 实际执行。

可以在配置时启用可选的 CPU 导数阶段计时：

```sh
cmake --preset cpu-release -DCOILGUN_ENABLE_CPU_PHASE_TIMING=ON
```

计时收集器默认关闭，默认构建不会增加公共运行时开销。

### 最小使用示例

```cpp
#include <coilgun/coilgun.hpp>
#include <iostream>

int main() {
    using namespace coilgun::physics;
    using coilgun::components::DrivingCoil;
    using coilgun::components::Armature;
    using coilgun::simulation::SingleStageSim;
    using coilgun::simulation::CrowbarExcitation;

    DrivingCoil coil(0.01, 0.03, 0.05, 150,
                     COPPER.resistivity_ref, 1e-6, 0.7);
    Armature arm(0.005, 0.025, 0.08,
                 ALUMINUM.resistivity_ref, ALUMINUM.density,
                 0.0, 0.120, 5, 2, 0.05);

    auto exc = std::make_unique<CrowbarExcitation>(450.0, 0.001);
    SingleStageSim<EulerStepper> sim(coil, arm, std::move(exc), 1e-6);

    sim.run();
    auto& r = sim.result();
    std::cout << "速度: " << r.summary.muzzle_velocity << " m/s\n"
              << "效率: " << r.summary.efficiency * 100 << " %\n";
}
```

### 链接库

```cmake
# 在 CMakeLists.txt 中
find_package(OpenMP REQUIRED)
target_link_libraries(your_target PRIVATE coilgun OpenMP::OpenMP_CXX)
```

或直接编译：

```sh
g++ -std=c++20 -fopenmp -Iinclude your_file.cpp build/src/libcoilgun.a -o your_binary
```

---

## 物理基础

库实现了 [NumericalModel.md](docs/NumericalModel.md) 中描述的数值模型：

| 模块 | 计算内容 | 核心公式 |
|--------|-----------------|-------------|
| 自感 | 空心圆柱线圈的 L | L = μ₀ · nc² · ri³ · T(q,p) |
| 互感（丝级） | 两个同轴环之间的 M | Maxwell 椭圆积分公式 |
| 互感（线圈级） | 两个厚线圈之间的 M | 4D Gauss-Legendre 求积 |
| 电路 ODE | 耦合线圈-电枢系统的 dI/dt | [dI/dt] = (L-M)^(-1) · (U + v·dM·I - R·I) |
| 力 | 轴向洛伦兹力 | F = Σ dM/dx · I_d · I_p |
| 热 | 绝热焦耳加热 | ΔT = I²RΔt / (m·c_p(T)) |

---

## 测试套件

CTest 根据所选 preset 发现当前目标集合；文档不固定测试目标数量。使用 label 分离快速 CPU 检查、数值 validation 和 CUDA 测试：

```sh
ctest --preset cpu-debug -L fast
ctest --preset cpu-debug -L validation
ctest --preset cuda-debug -L gpu
```

GPU 测试使用共享资源锁。需要物理 CUDA 设备的测试在无设备时 skip；真实 GPU CI 还必须检查上面的执行报告条件。`bench_gpu_engine` 不是 CTest 目标，必须显式运行。

---

## 文档

- [API 参考（英文）](docs/API.md) — 完整 C++ 函数与类 API
- [API 参考（中文）](docs/API_cn.md) — 中文版
- [数值模型](docs/NumericalModel.md) — 物理推导与算法详解
- [GPU benchmark 记录](docs/benchmarks/2026-07-19-unified-gpu-engine.md) — 机器相关的执行测量和回退状态
- [性能 Review](docs/PerformanceReview.md) — 全项目性能发现与 review 收束记录
- [优化成效](docs/benchmarks/2026-07-24-PerformanceOptimizationEffectiveness.md) — 本轮优化前后 benchmark 对比

---

## CMake 预设

| 预设 | 生成器 | 构建类型 | 标志 | 备注 |
|--------|-----------|------------|-------|-------|
| `cpu-debug` | Ninja | Debug | `-march=native`、`-O2 -g` | 仅 CPU，启用测试，生成 compile_commands.json |
| `cpu-release` | Ninja | Release | `-march=native` | 仅 CPU，启用测试 |
| `cuda-debug` | Ninja | Debug | `-march=native`、`-O2 -g` | 启用 CUDA 和测试，生成 compile_commands.json |
| `cuda-release` | Ninja | Release | `-march=native` | 启用 CUDA 和测试 |
| `cpu-release-library` | Ninja | Release | `-march=native` | 仅构建 CPU library |
| `cuda-release-library` | Ninja | Release | `-march=native` | 仅构建 CUDA library 目标 |

---

## 许可证

GNU General Public License v3.0。全文见 [LICENSE](LICENSE)。
