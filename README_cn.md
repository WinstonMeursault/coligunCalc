# coligunCalc

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language](https://img.shields.io/badge/language-C%2B%2B17-blue)](https://isocpp.org/)
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
- **OpenMP 多核并行** — M/dM 计算、力求和、温度更新均已并行化

---

## 项目结构

```
include/coilgun/   — 公开头文件（核心类型、物理层、组件层、仿真器）
src/               — 库实现（静态库 libcoilgun.a）
tests/             — 单元测试与集成测试（doctest，17 套）
tools/             — T(q,p) 查表生成工具
docs/              — 文档（API 参考中/英文、数值模型、设计文档）
.references/       — 参考论文 PDF（gitignored，仅本地）
```

## 快速开始

**依赖**：C++17 编译器（支持 OpenMP）、CMake ≥ 3.20。

Boost.Math 和 Eigen 均通过 CMake `FetchContent` 自动获取，无需手动安装。

```sh
cmake --preset ninja-debug       # 配置
cmake --build --preset ninja-debug   # 构建
ctest --preset debug              # 运行测试（启用 CUDA 时为 17/17）
```

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
g++ -std=c++17 -fopenmp -Iinclude your_file.cpp build/src/libcoilgun.a -o your_binary
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

## 测试套件（17 套，全部通过）

| 套件 | 覆盖范围 |
|-------|----------|
| `test_elliptic` | K(k)、E(k) 与已知值对比 |
| `test_struve` | H₀(x)、H₁(x) 与 SciPy 参考值对比 |
| `test_quadrature` | Gauss-Legendre/Laguerre 节点生成 |
| `test_lookup` | T(q,p) 表 + 双线性插值 vs 参考值 |
| `test_self_inductance` | 查表 vs 精确、边界情况、精确/比较 |
| `test_mutual_inductance` | 丝级和线圈级 M 及 dM/dz |
| `test_driving_coil` | 几何参数、R、L、匝密度 |
| `test_armature` | 丝几何、逐丝 R/L/质量 |
| `test_single_stage_sim` | 电容、续流、Euler/RK4、热模式 |
| `test_multi_stage_sim` | 单级等效性、两级、触发、优化等级 |
| `test_integration` | 端到端单级场景 |
| `test_gpu_elliptic` | GPU K(k)、E(k) vs CPU 参考 |
| `test_gpu_filament` | GPU 丝级 M、dM/dz vs CPU |
| `test_gpu_coil_pair` | GPU 线圈-丝 M、dM/dz vs CPU |
| `test_gpu_vs_cpu_single` | GPU 单级端到端 vs CPU (ε < 0.5%) |
| `test_gpu_vs_cpu_multi` | GPU 多级端到端 vs CPU (ε < 0.5%) |
| `test_gpu_batch` | GPU 批量仿真模式 |

---

## 文档

- [API 参考（英文）](docs/API.md) — 完整 C++ 函数与类 API
- [API 参考（中文）](docs/API_cn.md) — 中文版
- [数值模型](docs/NumericalModel.md) — 物理推导与算法详解
- [多级仿真设计](docs/multi_stage_sim_design.md) — 架构与实现方案
- [82mm 测试数据集](docs/test_dataset_82mm_coilgun.md) — 向洪 (2015) 验证数据

---

## CMake 预设

| 预设 | 生成器 | 构建类型 | 标志 | 备注 |
|--------|-----------|------------|-------|-------|
| `ninja-debug` | Ninja | Debug | `-march=native` | 启用测试，生成 compile_commands.json |
| `ninja-release` | Ninja | Release | `-march=native -O3` | — |
| `make-debug` | Unix Makefiles | Debug | `-march=native` | 启用测试，生成 compile_commands.json |

---

## 许可证

GNU General Public License v3.0。全文见 [LICENSE](LICENSE)。
