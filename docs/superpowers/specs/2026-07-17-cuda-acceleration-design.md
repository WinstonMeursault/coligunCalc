# CUDA/GPU 加速 — 完整实现方案

> 日期: 2026-07-17 | 分支: `feature/cudaAcceleration` | 状态: 设计已确认，待用户审查

## 1. 设计目标

在现有 C++17 线圈炮仿真库中引入 CUDA 后端，GPU 加速 per-step 的 4D Gauss-Legendre 积分计算（占仿真耗时 >95%）。

目标场景：
- **参数扫描**: 同一模型扫描不同初始电压/电容/触发位置，10–100 个仿真并发
- **高分辨率仿真**: S > 50 级，N_fil > 200 丝元，让 GPU 优势实质化
- **工程实践**: 打通 C++17 + Boost + Eigen + CUDA 的集成管线

## 2. 架构决策

### 2.1 条件编译分支

- 新增 `libcoilgun_cuda.a` 独立 CUDA 静态库，依赖 `libcoilgun.a`
- CMake option `COILGUN_ENABLE_CUDA=ON/OFF` 控制是否构建
- CPU 路径完全不动，GPU 路径通过独立类提供

### 2.2 API 分离 — 两种模式

```
Batch 模式 (SimBatch):
  ├─ 所有仿真共享相同的线圈几何 + 丝元离散化
  ├─ 构造时校验 GeometryFingerprint 保证同构
  └─ 一次性批量 kernel: grid = (num_sims, n_active × F)

Stream 模式 (GpuSingleStageSim / GpuMultiStageSim):
  ├─ 单仿真独立 GPU 类
  ├─ 几何模型自由，不要求同构
  └─ 多个仿真通过独立 cudaStream 并发
```

### 2.3 内存管理 — 全设备端驻留

- 初始化时一次性上传所有不变数据（线圈几何、丝元信息、GL 节点）
- 预期显存占用: 100 仿真 × ~2.5MB ≈ 250MB，在 16GB 内绰绰有余
- 每步仅传 armature 位置向量 + 取回 M1/dM1 矩阵

### 2.4 GPU 优化等级

| 等级 | 距离截断 | 自适应 GL | 用途 |
|---|---|---|---|
| `Standard` | 否 | 否 (固定 n_nodes=9) | 调试 / 验证 |
| `Full` | 是 (>10×线圈长度) | 是 (近=9, 远=4) | 生产 |

原 CPU 路径的 T(q,p) 查表在构造时完成，GPU 路径自然继承，不需要单独优化等级。

## 3. Kernel 设计

### 3.1 单对积分 Kernel

```
kernel_mutual_coil_pair <<< grid, 512, shared_mem, stream >>>

  每个 block 处理 1 个 (stage, filament) 对:
    - 512 threads × ~13 次循环 = 覆盖全部 6561 个 4D 积分点
    - 每个点同时计算 M 和 dM (避免重复椭圆积分)
    - shared memory 树形归约 → tid=0 写入 global memory
```

### 3.2 Batch Kernel

```
grid = (num_sims, n_active × F_per_sim)
  blockIdx.y → 仿真 ID
  blockIdx.x → (stage, filament) 对索引
  从全局常量缓冲区按 sim_id 索引几何参数
```

### 3.3 设备内存布局

```
常量区 (一次上传):
  coil_geom[num_sims][n_stages]  = {ri, re, len, turns, pos}
  fil_geom[num_sims][N_fil]      = {ri, re, len, mean_radius}
  GL_nodes[9], GL_weights[9]

每步传输:
  H→D: separation[num_sims][n_stages × N_fil]   (位置向量)
  D→H: M_results[num_sims][n_stages × N_fil]
       dM_results[num_sims][n_stages × N_fil]
```

## 4. 公开 API

### 4.1 类型

```cpp
namespace coilgun::simulation::cuda {

enum class GpuOptLevel {
    Standard = 0,   // n_nodes=9 固定，不截断
    Full     = 1,   // 距离截断 + 自适应 n_nodes (近9远4)
};

struct GpuBackend {
    int  device_id        = 0;
    int  threads_per_block = 512;
    size_t max_batch_sims = 256;
    bool enable_profiling = false;
};

}
```

### 4.2 SimBatch（Batch 模式）

```cpp
template<typename SP>
class SimBatch {
public:
    SimBatch(std::vector<DrivingCoil> coils, Armature armature,
             int num_sims, double dt, const GpuBackend& = {});

    void set_excitations(int sim_id,
        std::vector<std::unique_ptr<Excitation>> excitations,
        std::vector<TriggerConfig> triggers);

    void run();
    void run(const TerminationPolicy&);
    const std::vector<MultiStageResult>& results() const;
    const MultiStageResult& result(int sim_id) const;

    int  num_sims() const;
};
```

### 4.3 GpuMultiStageSim / GpuSingleStageSim（Stream 模式）

```cpp
template<typename SP>
class GpuMultiStageSim {
public:
    GpuMultiStageSim(
        std::vector<DrivingCoil> coils, Armature armature,
        std::vector<std::unique_ptr<Excitation>> excitations,
        std::vector<TriggerConfig> triggers,
        double dt, bool enable_thermal = false,
        GpuOptLevel opt_level = GpuOptLevel::Full,
        const GpuBackend& = {});

    const MultiStageStep&  step();
    const MultiStageResult& run();
    const MultiStageResult& run(const TerminationPolicy&);
    void reset();

    const MultiStageResult& result() const;
    const MultiStageState&  state() const;
    double dt() const;
    int step_count() const;
    int num_stages() const;
};
```

API 与现有 CPU 版 `MultiStageSim` / `SingleStageSim` 保持一致，迁移仅需改类名和 include。

## 5. 头文件策略 (Plan B)

```
include/coilgun/physics/
  elliptic.hpp / mutual_inductance.hpp  ← 现存，CPU-only，不动
  elliptic.cuh / mutual_inductance.cuh  ← 新增，__host__ __device__ 内联声明

src/physics/
  elliptic.cpp / mutual_inductance.cpp  ← 现存

src/cuda/
  gpu_elliptic.cu           ← #include "elliptic.cuh"，device 端实现
  gpu_mutual_inductance.cu  ← kernel 实现
  gpu_adaptor.cu            ← 几何数据上传 / 驻留管理
  gpu_single_stage_sim.cu
  gpu_multi_stage_sim.cu
```

CUDA 编译时 `BOOST_MATH_ENABLE_CUDA` 宏启用 Boost.Math 的 device 端支持。

## 6. 构建系统

```cmake
# 顶层
project(coilgunCalc VERSION 1.2.0 LANGUAGES C CXX C)
option(COILGUN_ENABLE_CUDA "Build GPU-accelerated backend" OFF)

# src/cuda/CMakeLists.txt
add_library(coilgun_cuda STATIC ...)
target_link_libraries(coilgun_cuda
    PUBLIC coilgun Eigen3::Eigen Boost::headers
    PRIVATE CUDA::cudart)
target_compile_definitions(coilgun_cuda PRIVATE BOOST_MATH_ENABLE_CUDA)

# 目标架构: 默认 native 自动检测, 预设覆盖
set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING "CUDA arch")
# 预设中可写: "CMAKE_CUDA_ARCHITECTURES": "120"
```

CUDA 13.3 完整支持 Blackwell sm_120 (RTX 5080 Laptop)。

## 7. 测试策略

### 7.1 层级

| 层 | 测试 | 框架 | 精度要求 |
|---|---|---|---|
| L1 单元 | 现有 12 套 CPU 测试 | doctest | 不变 |
| L2 GPU 正确性 | `test_gpu_elliptic`, `test_gpu_filament`, `test_gpu_coil_pair`, `test_gpu_single_step` | doctest | 见下表 |
| L3 端到端 | `test_gpu_vs_cpu_single`, `test_gpu_vs_cpu_multi` | doctest | 见下表 |
| L4 Batch | `test_gpu_batch` (10 仿真) | doctest | 见下表 |

### 7.2 精度标准

| 测试 | 容差 | 理由 |
|---|---|---|
| GPU elliptic vs CPU | eps ≤ 1e-14 | 同算法 |
| GPU filament M vs CPU | rel ≤ 1e-13 | 浮点顺序差异 |
| GPU 4D 积分 vs CPU | rel ≤ 5e-7 | T(q,p) 查表量级 |
| 单步仿真 vs CPU | rel ≤ 1e-10 | 累加误差 |
| 端到端 v_muzzle vs CPU | rel ≤ 1e-5 | 多步误差累积 |

## 8. 实现阶段

| 阶段 | 内容 | 验证方式 |
|---|---|---|
| **M1** | 构建系统 + `elliptic.cuh` / `gpu_elliptic.cu` | `test_gpu_elliptic` |
| **M2** | `gpu_mutual_inductance.cu` — 单对 kernel + Stream | `test_gpu_filament` + `test_gpu_coil_pair` |
| **M3** | `GpuSingleStageSim` — 单级完整 CPU/GPU 混合循环 | `test_gpu_vs_cpu_single` |
| **M4** | `GpuMultiStageSim` — 多级仿真 | `test_gpu_vs_cpu_multi` |
| **M5** | `SimBatch` — Batch 模式容器 | `test_gpu_batch` |
| **M6** | 文档更新 (API.md / API_cn.md / README) | review |

## 9. 文件清单

### 新增文件

```
include/coilgun/physics/
  elliptic.cuh
  mutual_inductance.cuh
include/coilgun/simulation/cuda/
  gpu_backend.hpp
  gpu_adaptor.hpp
  gpu_single_stage_sim.hpp
  gpu_multi_stage_sim.hpp

src/cuda/
  CMakeLists.txt
  gpu_elliptic.cu
  gpu_mutual_inductance.cu
  gpu_adaptor.cu
  gpu_single_stage_sim.cu
  gpu_multi_stage_sim.cu

tests/
  test_gpu_elliptic.cu
  test_gpu_filament.cu
  test_gpu_coil_pair.cu
  test_gpu_single_step.cu
  test_gpu_vs_cpu_single.cu
  test_gpu_vs_cpu_multi.cu
  test_gpu_batch.cu
```

### 修改文件

```
CMakeLists.txt                        → 项目名/版本/COILGUN_ENABLE_CUDA
tests/CMakeLists.txt                  → 条件添加 GPU 测试
docs/API.md                           → GPU API 章节
docs/API_cn.md                        → GPU API 章节
docs/CUDA-feasibility.md              → 更新为实际结论
README.md / README_cn.md              → 增加 CUDA 构建说明
```

### 不变文件

所有 `src/physics/`, `src/components/`, `src/simulation/` 下的现有代码保持不动。

## 10. 边界条件与降级

| 场景 | 行为 |
|---|---|
| `n_active = 0` (所有级已结束) | 跳过 kernel launch，直接返回零矩阵 |
| `n_nodes = 4` (自适应远距) | 256 个积分点，512 threads → 1 轮循环即完成，半数 thread 空闲（可接受） |
| SimBatch 含不同几何 | 构造时 `GeometryFingerprint` 校验 → 抛出 `std::invalid_argument` |
| `COILGUN_ENABLE_CUDA=OFF` | 无 `.cuh`/GPU 头文件可用，编译报错（需用户感知） |
| CUDA toolkit < 12.8 | `sm_120` 不可用，CMake 报 warning；用户需降级架构或升级 CUDA |
| 显存不足 | `cudaMalloc` 失败 → 抛出 `std::runtime_error`，附带所请求的字节数 |
| 热模 (enable_thermal) | 温度更新在 CPU 端执行（<1% 耗时），GPU 只做 M/dM 积分 |

