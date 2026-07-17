# CUDA/GPU 加速可行性评估

> 日期: 2026-07-17 | 针对分支 `refactor/cpp`

## 计算热区分析

| 计算阶段 | 每步耗时占比 | 核心操作 | GPU 适配度 |
|---|---|---|---|
| `update_M1_dM1()` (4D GL 积分) | **>95%** | 每对(台级,丝元) n⁴ 次椭圆积分求值 | **极高 (数据并行)** |
| 线性求解 (Eigen LDLT) | <1% | (S+F)×(S+F) 矩阵分解 | 低 (串行为主) |
| 力计算 `compute_force()` | <1% | S×F 点积归约 | 中 |
| 温度更新 `update_temperatures()` | <1% | 逐丝元 cp/ρ 查找 | 低 |

**瓶颈**: 4D Gauss-Legendre 积分的四层嵌套循环 (`integrate_4d`, `src/physics/mutual_inductance.cpp:190-210`)。n_nodes=9 时每个 (台级,丝元) 对需 6561 次独立椭圆积分求值，天然适合 GPU 线程映射。

典型规模: 多级 S=25, N_fil=45 → 1125 对/步 → 约 7.4M 次椭圆积分/步。

## 椭圆积分的 GPU 实现现状

### 现成方案: Boost.Math ≥ 1.86 (2024)

Boost.Math 在 v1.86 版本中由 Matt Borland 引入了完整的 GPU 支持 (`BOOST_MATH_GPU_ENABLED`)。

```cpp
// boost/math/tools/config.hpp
#if defined(__CUDACC__) && defined(BOOST_MATH_ENABLE_CUDA)
#  define BOOST_MATH_GPU_ENABLED __host__ __device__
#  include <cuda/std/type_traits>
// ... CUDA runtime 集成
#endif
```

`ellint_1.hpp` 和 `ellint_2.hpp` 中的所有函数（含多项式逼近优化路径）均已标注 `BOOST_MATH_GPU_ENABLED`，可在 `__device__` 上下文中直接调用。

### 使用方式

编译时传入宏即可:

```cmake
# CMakeLists.txt
target_compile_definitions(coilgun_cuda PRIVATE BOOST_MATH_ENABLE_CUDA)
```

```cpp
// elliptic.cu
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

__host__ __device__ double elliptic_k(double m) {
    return boost::math::ellint_1(std::sqrt(m));     // K(k) — GPU 可用
}
__host__ __device__ double elliptic_e(double m) {
    return boost::math::ellint_2(std::sqrt(m));     // E(k) — GPU 可用
}
```

### 其他选择

| 方案 | 状态 | 说明 |
|---|---|---|
| NVIDIA CUDA Math API | **不包含** | CUDA 内置数学 API 仅有标准函数 (sin/cos/exp...)，无椭圆积分 |
| 自定义 AGM/Carlson 实现 | 可行但冗余 | 5-10 次迭代收敛，可手写 `__device__`，但 Boost 已做此工作 |
| 多项式逼近查表 | 可行 | 项目已有 T(q,p) 查表(23MB)，可类比; 但需额外维护 |

**结论**: 直接使用 Boost.Math 的 CUDA 支持是唯一合理方案。

## 有利因素

1. **4D 积分是典型的 embarrassingly parallel**: n⁴ 个互不依赖的求值点，天然映射到 GPU 线程 grid
2. **多级展平循环**: 项目已将 `(stage, filament)` 双重循环展平为 `for (flat = 0; flat < n_active * N_fil; ++flat)` (见 `multi_stage_sim.cpp`)，正好对应 CUDA 的 1D block 调度
3. **每步数据传输量小**: 仅 armature 位置变化，几何参数和线圈参数均为常量
4. **无缓存冲突**: GPU 共享内存/寄存器可替代 CPU 端的 `LRUCache<std::tuple<double,double,double>>`
5. **精度需求适中**: FP64 精度，NVIDIA A100/H100 提供良好 FP64 吞吐 (9.7 TFLOPS)
6. **Boost.Math 零额外依赖**: 项目已通过 CMake FetchContent 拉取 Boost.Math

## 主要障碍

1. **问题规模偏小**: S=25, N_fil=45 (1125 对/步) 对 GPU 而言偏小，kernel launch 开销可能占比大。若升级到高分辨率 (N_fil > 200, S > 50, n_nodes=16)，GPU 优势更显著
2. **消费级 GPU FP64 性能弱**: RTX 3090/4090 FP64:FP32 吞吐比约 1:32~1:64。需数据中心 GPU (A100/H100) 才能发挥 FP64 性能
3. **PCle 传输 + kernel launch 延迟**: 每步需传输 armature 位置到 GPU，取回 M1/dM1 矩阵。需衡量此开销是否超过 CPU OpenMP 的计算耗时
4. **Eigen 线性求解器不能跑 GPU**: LDLT 分解需保持在 CPU 端，但仅占总耗时的 <1%
5. **双代码路径维护成本**: 需同时维护 CPU (OpenMP) 和 GPU (CUDA) 两套实现

## 建议实现路径

### 方案: 混合 CPU/GPU (推荐)

```
每步循环:
┌─────────────────────────────────────────────────┐
│ CPU (OpenMP)                                    │
│   • 检查触发器 / 灭活静默台级                     │
│   • 构建 active_idx 列表                        │
├─────────────────────────────────────────────────┤
│ GPU (CUDA kernel)                               │
│   • 传输 armature_position → device             │
│   • kernel: flat_num 个 block, n²×n² 个 thread  │
│     每个 thread 计算 M_filament + dM_filament    │
│     块内归约完成 4D 积分                          │
│   • 传输 M1_mat, dM1_mat → host                  │
├─────────────────────────────────────────────────┤
│ CPU                                            │
│   • 构建 L_total 系统矩阵                        │
│   • Eigen LDLT 线性求解                          │
│   • 力计算 (OpenMP reduction)                    │
│   • 温度更新 (OpenMP)                            │
│   • 步进器推进                                   │
└─────────────────────────────────────────────────┘
```

### 需要修改的文件

| 文件 | 变更 |
|---|---|
| `include/coilgun/physics/elliptic.hpp` | 标注 `__host__ __device__` |
| `include/coilgun/physics/mutual_inductance.hpp` | 标注 `__host__ __device__` |
| `src/physics/mutual_inductance.cu` | 新增 CUDA kernel 实现 |
| `src/simulation/multi_stage_sim.cpp` | 添加 GPU 路径 |
| `CMakeLists.txt` | `enable_language(CUDA)`, `BOOST_MATH_ENABLE_CUDA` |

## 总体评估

| 维度 | 评分 | 说明 |
|---|---|---|
| 理论加速潜力 | ★★★☆☆ | 4D 积分可并行化，但受限于问题规模和 FP64 性能 |
| 工程复杂度 | ★★★☆☆ | 借助 Boost.Math 内置 CUDA 支持，无需重写特殊函数; 需写 kernel + 数据搬运 |
| ROI (投入产出比) | ★★☆☆☆ | 当前规模下 OpenMP 8-16 核已足够。ROI 随问题规模增长而升高 |
| 硬件兼容性 | ★★☆☆☆ | 需 A100/H100 才有好的 FP64 性能; 消费级 GPU FP64 极弱 |

## 最终结论

**当前不建议投入 GPU 加速。** 理由:

1. Boost.Math 的 CUDA 支持消除了重写椭圆积分的障碍，使 GPU 实现技术上可行
2. 但当前问题规模 (25 台级 × 45 丝元) 下，OpenMP 并行化 (已实现) 在 8-16 核 CPU 上提供足够好的性能
3. 引入 GPU 双路径将显著增加维护成本

**升级建议**: 若未来需高分辨率仿真 (N_fil > 500, S > 100, n_nodes ≥ 16)，GPU 加速的性价比将大幅提升。届时可直接利用 Boost.Math 的内置 CUDA 支持进行实施。
