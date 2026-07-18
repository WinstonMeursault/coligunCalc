# GPU 加速优化 — 第二阶段设计 (v2)

> 日期: 2026-07-18 | 分支: `feature/cudaAcceleration` | 状态: 设计已确认（第二轮细化）

## 变更记录

| 版本 | 日期 | 变更 |
|------|------|------|
| v1 | 2026-07-18 | 初版：持久化 kernel + P1-P7 后续方向 |
| v2 | 2026-07-18 | 第二轮细化：覆盖全部三个 GPU 类、GpuOptLevel 扩展三级、P7 提升优先级、FP32 融合 |

## 背景

第一阶段实现已完成：条件编译、GpuSingleStageSim、GpuMultiStageSim、SimBatch、API 文档、CMake 预设、19 套测试通过。当前性能：GPU RTX 5080 Laptop 仅比 16 核 CPU 快 1.2×（当前问题规模 S=2, N_fil=10）。

**核心瓶颈**：per-step 的 kernel launch overhead（20 次 launch/步 × 20000 步 = 400,000 次 launch）。

P0 batch kernel（grid 含 sim_id 维度）因 nvcc/nvlink 跨文件 struct-layout 问题阻塞。映射内存通路已验证可用——**持久化 kernel 方案可行**。

---

## 一、持久化 Kernel（P0 — 主方向）

### 1.1 覆盖范围（v2 扩展）

**v1 设计仅覆盖 SimBatch。v2 扩展为覆盖全部三个 GPU 仿真类**：

| 类 | 持久化模式 block 数 | 典型场景 |
|----|--------------------|----------|
| `GpuSingleStageSim` | N_fil | 1 stage × N 丝元 |
| `GpuMultiStageSim` | n_stages × N_fil | S 级 × N 丝元 |
| `SimBatch` | n_stages × N_fil | 同上，叠加 N_sim 个仿真 |

实现顺序：D2 Single → D3 Multi → D4 Batch（由简到繁）。

### 1.2 代码架构（v2 变更）

**v1**：kernel 定义在 `sim_batch.cu` 中。  
**v2**：kernel 定义在共享 `persistent_kernel.cuh` 头文件中，由三个 .cu 文件 `#include`。避免 CUDA 跨文件 extern `__global__` 的链接问题。

```
include/coilgun/simulation/cuda/
├── gpu_backend.hpp              ← + use_persistent (bool)
├── persistent_kernel.cuh        ← [新增] PersistentBuffers 结构 + persistent_batch_kernel 定义
├── gpu_single_stage_sim.hpp     ← + PersistentBuffers 成员 + init/free/launch 声明
├── gpu_multi_stage_sim.hpp      ← + PersistentBuffers 成员 + init/free/launch 声明
└── sim_batch.hpp                ← + PersistentBuffers 成员 + init/free/launch 声明

src/cuda/
├── gpu_single_stage_sim.cu      ← 修改 compute_M1_dM1 支持持久化路径
├── gpu_multi_stage_sim.cu       ← 修改 compute_M1_dM1 支持持久化路径
└── sim_batch.cu                 ← 修改 compute_all_M1_dM1 支持持久化路径
```

### 1.3 架构

```
Host (SimBatch::run):                Device (persistent kernel):
  ┌─ cudaHostAllocMapped ─┐         ┌───────────────────────────┐
  │  batch_id              │←──映射──│ ← block 读 → 新批次触发     │
  │  active_pairs          │←──映射──│ ← block 读 → 活跃对数目     │
  │  shutdown              │←──映射──│ ← block 读 → 退出标志       │
  │  doorbell[N_b]         │←──映射──│ ← host 写1 → block 读0      │
  │                        │         │ ← block 写0 → host 读0      │
  │  seps[N_b]             │←──映射──│ ← block 读 → 分离参数       │
  │  results_M[N_b]        │──映射──→│ ← block 写 → M 结果        │
  │  results_dM[N_b]       │──映射──→│ ← block 写 → dM 结果       │
  └────────────────────────┘         └───────────────────────────┘

N_b = n_stages × N_fil
```

### 1.4 Per-Step 数据流（适用于所有三个类）

```
Host side (step N)                    Device side (persistent kernel)
─────────────────                     ──────────────────────────────
1. check_triggers/extinguish          │
2. 填充 seps[0..K-1], K=活跃对数       │  while batch_id == last: spin
3. active_pairs = K                   │  last = batch_id
4. doorbell[0..K-1] = 1               │  while doorbell[idx]==0: spin
5. batch_id++                         │  if idx < active_pairs:
                                      │    4D GL → M, dM
                                      │  out_M[idx]=M, out_dM[idx]=dM
                                      │  doorbell[idx]=0
6. while any doorbell[i]!=0: spin     │
7. cudaStreamSynchronize              │
8. 读取 out_M/out_dM → M1_mat/dM1_mat │
9. LDLT solve → force → kinematics    │
```

**GpuSingleStageSim**：K = N_fil（1 stage，所有丝元始终活跃）。  
**GpuMultiStageSim**：K = n_active × N_fil（动态变化）。  
**SimBatch**：外层 for 循环逐 sim 触发，K 同上。

### 1.5 空闲 Block 的 SM 占用

所有 N_b 个 block 中，仅 K 个执行实际计算。剩余 (N_b - K) 个 block 只参与 doorbell 同步（读 doorbell[idx]，写 0）。

CUDA 硬件调度器限制同时驻留 SM 的 block 数：RTX 5080（84 SM × 4 block/SM）= **最多 336 block**。极端情况下（S=50, F=500, N_b=25000），仅 ~336 block 实际占用 SM，其中大部分空闲。spin 循环（load+compare+branch）极其轻量，<0.01% 总 GPU 吞吐。

**决策：接受 spin-wait 作为降低 launch overhead 的工程代价。**

### 1.6 Kernel 内部

```cuda
__global__ void persistent_batch_kernel(
        const CoilGeo* coils, const FilGeo* fils,
        const double* gl_nodes, const double* gl_weights,
        volatile int* batch_id,
        volatile int* active_pairs,
        volatile int* shutdown,
        volatile double* seps,
        volatile int* doorbell,
        volatile double* out_M,
        volatile double* out_dM,
        int n_stages, int N_fil, int n_nodes) {

    int idx = blockIdx.x;
    int si = idx / N_fil; int fi = idx % N_fil;
    int last_batch = -1;

    // Precompute constant geometry for this block
    const auto& coil = coils[si];
    const auto& fil  = fils[fi];
    const double ra_mid  = 0.5 * (coil.re + coil.ri);
    const double ra_half = 0.5 * (coil.re - coil.ri);
    const double rb_mid  = 0.5 * (fil.re + fil.ri);
    const double rb_half = 0.5 * (fil.re - fil.ri);
    const double la_half = 0.5 * coil.len;
    const double lb_half = 0.5 * fil.len;
    const double prefactor = coil.turns / 16.0;

    int n2 = n_nodes * n_nodes, n4 = n2 * n2;
    int tid = threadIdx.x, total = blockDim.x;

    while (*shutdown == 0) {
        // Wait for next batch
        int cur = *batch_id;
        if (cur == last_batch) continue;
        last_batch = cur;

        // Wait for host to fill doorbell
        while (doorbell[idx] == 0) {}

        double M = 0.0, dM = 0.0;

        if (idx < *active_pairs) {
            double separation = seps[idx];

            __shared__ double sM[512], sdM[512];
            sM[tid] = 0.0; sdM[tid] = 0.0;

            for (int pt = tid; pt < n4; pt += total) {
                int i1 = pt / (n_nodes*n_nodes*n_nodes);
                int rem1 = pt % (n_nodes*n_nodes*n_nodes);
                int j1 = rem1 / (n_nodes*n_nodes);
                int rem2 = rem1 % (n_nodes*n_nodes);
                int i2 = rem2 / n_nodes;
                int j2 = rem2 % n_nodes;

                double w = gl_weights[i1] * gl_weights[j1]
                         * gl_weights[i2] * gl_weights[j2];

                double ra = ra_mid + ra_half * gl_nodes[i1];
                double rb = rb_mid + rb_half * gl_nodes[i2];
                double za = la_half * gl_nodes[j1];
                double zb = separation + lb_half * gl_nodes[j2];

                double abs_sep = fabs(zb - za);
                sM[tid]  += w * mutual_inductance_filament_device(ra, rb, abs_sep);
                sdM[tid] += w * mutual_inductance_gradient_filament_device(ra, rb, zb - za);
            }

            __syncthreads();
            for (int stride = total/2; stride > 0; stride /= 2) {
                if (tid < stride) { sM[tid]+=sM[tid+stride]; sdM[tid]+=sdM[tid+stride]; }
                __syncthreads();
            }

            if (tid == 0) {
                M  = prefactor * sM[0];
                dM = prefactor * sdM[0];
            }
        }

        __threadfence();
        if (tid == 0) { out_M[idx] = M; out_dM[idx] = dM; }
        __threadfence();
        __syncthreads();
        if (tid == 0) doorbell[idx] = 0;
    }
}
```

### 1.7 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 同步方式 | 逐 block 门控 + batch_id | 无 atomic 竞争，无 false sharing |
| block 数 | 固定 N_b（含闲置 block）| 启动一次不重新 launch |
| 内存通道 | `cudaHostAllocMapped` | 绕过 P0 的跨文件 cudaMemcpy 问题 |
| per-step 屏障 | 逐 doorbell 轮询 + `cudaStreamSynchronize` | 保证 host 可见 device 完成 |
| 降级开关 | `GpuBackend::use_persistent` | debug 回退 |
| kernel 位置 | `.cuh` 头文件，#include 到各 `.cu` | 避免跨文件 extern `__global__` |
| 空闲 block | 参与 doorbell 同步，不计算 | 统一同步协议，简化实现 |

### 1.8 内存量

最坏情况 S=50, N_fil=500 → max_pairs=25000，单精度映射内存：
- doorbell: 25000 × 4 = 100 KB
- seps: 25000 × 8 = 200 KB
- out_M: 25000 × 8 = 200 KB
- out_dM: 25000 × 8 = 200 KB
- batch_id/active_pairs/shutdown: 忽略不计
- **总计: ~700 KB**，远低于 16 GB 限制

### 1.9 新增/修改文件

| 文件 | 变更 |
|------|------|
| `include/coilgun/simulation/cuda/persistent_kernel.cuh` | **[新增]** PersistentBuffers 结构 + kernel 定义 |
| `include/coilgun/simulation/cuda/gpu_backend.hpp` | GpuBackend 新增 `bool use_persistent = true`；GpuOptLevel 扩展为三级 |
| `include/coilgun/simulation/cuda/gpu_single_stage_sim.hpp` | 新增 PersistentBuffers 成员 + init/free/launch 声明 |
| `include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp` | 同上 |
| `include/coilgun/simulation/cuda/sim_batch.hpp` | 同上 |
| `src/cuda/gpu_single_stage_sim.cu` | 修改 compute_M1_dM1 支持持久化路径 |
| `src/cuda/gpu_multi_stage_sim.cu` | 同上 |
| `src/cuda/sim_batch.cu` | 修改 compute_all_M1_dM1 支持持久化路径 |
| `tests/test_gpu_vs_cpu_single.cpp` | 现有测试自动覆盖持久化（默认 on） |
| `tests/test_gpu_vs_cpu_multi.cpp` | 同上 |
| `tests/test_gpu_sim_batch.cpp` | 新增 persistent/fallback 一致性测试 |

---

## 二、GpuOptLevel 扩展（v2 变更）

### 2.1 三级优化等级

**v1** 为两级（Standard / Full）。**v2** 扩展为三级，将 FP32 混合精度融合为第三级：

```cpp
enum class GpuOptLevel {
    Standard   = 0,  // FP64, 无截断, n_nodes=9     — 验证/调试
    Full       = 1,  // FP64, 距离截断, n_nodes=9   — 当前生产默认
    Aggressive = 2,  // FP32 integrand + FP64 reduction, 距离截断, n_nodes=9 — 大规模扫描
};
```

| 等级 | 精度 | 距离截断 | 椭圆积分 | 场景 |
|:---:|:---:|:---:|:---:|------|
| Standard (0) | FP64 | 无 | Boost.Math FP64 | 验证/调试/reference |
| Full (1) | FP64 | 距离截断 (>10×) | Boost.Math FP64 | 生产默认 |
| Aggressive (2) | FP32 integrand + FP64 reduction | 距离截断 (>10×) | 手写 AGM FP32 | 大规模参数扫描 |

### 2.2 设计理由

- **累积语义**：与 CPU 端 `OptimizationLevel`（Reference → LookupTable → Full）一致
- **防止组合爆炸**：如果 `use_fp32` 是独立 bool，与 `use_cutoff`、`n_nodes` 组合产生 8 种配置，大量无意义组合
- **用户心智模型简单**：等级越高 = 计算越快 = 精度越低

### 2.3 FP32 椭圆积分（Aggressive 级别）

Boost.Math 椭圆积分仅支持 FP64 device-side。需要手写 AGM（Arithmetic-Geometric Mean）FP32 版本。

**AGM 方法**（NumericalModel.md §4.4.2）收敛速度快（~5-8 次迭代到机器精度），FP32 版本的 `float` mantissa 仅 23 bits（~7 位十进制），5 次 AGM 迭代足够。

内核改造：Kernel 在 Aggressive 模式下使用 `float` 类型做 4D 积分内层：
- GL 坐标变换、权重乘积 → `float`
- 椭圆积分 → 手写 AGM `float` 版本
- 互感 integrand → `float`
- 共享内存累加 → **`double`**（防止 FP32 归约积累误差）
- 最终结果（M, dM）→ `double`

预期加速：3-5×（RTX 5080: FP32 ~50 TFLOPS vs FP64 ~1.5 TFLOPS）。

### 2.4 注意

Aggressive 级别中，共享内存归约使用 `double` 累加器（非 float）。单 block 内 512 threads，sM[512] 和 sdM[512] 各 4096 字节，共 8192 字节，远低于 48 KB shared memory 限制。

---

## 三、GpuBackend API 变更

```cpp
// gpu_backend.hpp
enum class GpuOptLevel {
    Standard   = 0,  // FP64, no cutoff, n_nodes=9
    Full       = 1,  // FP64, distance cutoff, n_nodes=9
    Aggressive = 2,  // FP32 integrand + FP64 reduction, distance cutoff, n_nodes=9
};

struct GpuBackend {
    int     device_id         = 0;
    int     threads_per_block = 512;
    size_t  max_batch_sims    = 256;
    bool    enable_profiling  = false;
    bool    use_persistent    = true;   ///< [新增] 使用持久化 kernel。false 回退到 per-pair launch。
};
```

**降级策略**：
- `use_persistent = false` → 原 per-pair kernel launch（行为完全不变）
- `use_persistent = true` + `cudaHostAllocMapped` 成功 → 持久化 kernel
- `use_persistent = true` + `cudaHostAllocMapped` 失败 → 自动降级到 per-pair launch + stdout warning

**用户无感升级**：默认 `use_persistent = true`，现有用户代码无需修改。

---

## 四、其余优化方向（按优先级排列）

### P7：双缓冲 CPU/GPU 重叠（优先级: 高 — v2 提升）

**目标**：GPU 计算当前步 M1/dM1 的同时，CPU 并行执行上一步的 LDLT 求解。

**实现**：
- 两套 mapped 结果缓冲区 A/B 轮换
- 持久化 kernel 写入缓冲区 A → host 读取 B 做 LDLT
- 复用持久化 kernel 的 doorbell 机制
- 增加 `host_done` 标志

**预期加速**：1.3-1.5×。高分辨率（S>50, F>500）LDLT 占比更高，收益更大。

**风险**：中等。双缓冲区使状态管理复杂化。

---

### P2（原）：FP32 混合精度（优先级: 中 — 已融合到 GpuOptLevel::Aggressive）

**状态**：融入了三级 `GpuOptLevel`，作为 Persistent Kernel 完成后的第二优先实现项。

**目标**：利用 RTX 5080 的高 FP32 吞吐（~50 TFLOPS vs FP64 ~1.5 TFLOPS）。

**实现**：
- 4D 积分的 integrand 用 `float` 计算（包括 AGM 椭圆积分）
- 共享内存累加器和最终结果用 `double`（防止精度损失）
- 通过 `GpuOptLevel::Aggressive` 控制

**预期加速**：3-5×。

**风险**：中等。需手写 AGM 椭圆积分 FP32 版本。

---

### P1：CUDA Stream 并发（优先级: 中）

**目标**：多仿真并发时 GPU kernel 自动重叠。

**实现**：
- GpuMultiStageSim 构造时可选 `cudaStream_t`
- `cudaMemcpyAsync` 替代 `cudaMemcpy`
- 外部管理多个实例在不同 stream 上

**预期加速**：零额外 effort，RTX 5080 84 SM 可同时执行多个 kernel。

**风险**：低。需注意多个持久化 kernel 同时运行时的 SM 资源分配。

---

### P3：GL 节点 2D Block 并行（优先级: 低）

**目标**：提高 warp 占用率，减少 `__syncthreads` 开销。

**实现**：
- Block 维度从 `(512, 1)` 改为 `(9, 9)` = 81 threads
- `threadIdx.y` 映射外两层，`threadIdx.x` 循环内两层

**预期加速**：~1.2×。

**风险**：低。但 81 threads/block 远低于 warp scheduler 最佳 occupancy。

---

### P4：GPU LDLT 求解（优先级: 低）

**目标**：消除 M1/dM1 结果的 PCIe 传输，直接在 GPU 上求解。

**实现**：
- cuSOLVER `cusolverDnDpotrf` + `cusolverDnDpotrs` 替代 Eigen LDLT
- 仅传输电流向量回 host（~12 doubles）

**预期加速**：低分辨率忽略不计，S>50 时 2-3×。

**风险**：中等。需 cuSOLVER 依赖。要求正定矩阵。

---

### P5：CUDA Graph（优先级: 低）

**目标**：一次性录制 + 无限重放，零 launch overhead。

**实现**：
- 录制完整 step() 调用图
- 需固定 n_active（不支持动态触发/终止）

**预期加速**：~1.1-1.2×（相比持久化 kernel 的 0 launch overhead 仍有差距）。

**风险**：中等。不支持条件分支。

---

### P6（原）：增强型 Batch Kernel（优先级: 已取代）✅

**状态**：被持久化 kernel **完整取代**。持久化 kernel 在功能上完全覆盖原 batch kernel 设计目标，且消除了实现阻塞问题。不再作为独立方向。

---

### P8：Cooperative Groups 归约（优先级: 低）

**目标**：减少 shared memory 同步开销。

**实现**：
- `cooperative_groups::tiled_partition<32>` 做 warp 级归约
- 9 次 `__syncthreads()` → 1 次 cross-warp shuffle + 1 次 block sync

**预期加速**：~1.1×。

**风险**：低。

---

## 五、实现阶段（v2 重排）

| 阶段 | 内容 | 文件 | 验证 |
|------|------|------|------|
| **D1** | 共享 `persistent_kernel.cuh`：PersistentBuffers + kernel 定义 + init/free/launch 模板 | 新增 .cuh 头文件 + gpu_backend.hpp 修改 | 编译通过 |
| **D2** | `GpuSingleStageSim` 持久化集成 | gpu_single_stage_sim.{hpp,cu} | `test_gpu_vs_cpu_single` 通过 |
| **D3** | `GpuMultiStageSim` 持久化集成 | gpu_multi_stage_sim.{hpp,cu} | `test_gpu_vs_cpu_multi` 通过 |
| **D4** | `SimBatch` 持久化集成 | sim_batch.{hpp,cu} | `test_gpu_sim_batch` persistent/fallback 一致性 |
| **D5** | 双缓冲 CPU/GPU 重叠 | 各 GPU 类 | 性能对比测试 |
| **D6** | `GpuOptLevel::Aggressive`（FP32 混精 + AGM 椭圆积分）| persistent_kernel.cuh + mutual_inductance.cuh | 精度验证 + 性能对比 |
| **D7** | 文档更新 | API.md / API_cn.md | review |

---

## 六、测试策略

### 现有测试（自动覆盖持久化路径）

| 测试 | 覆盖 |
|------|------|
| `test_gpu_vs_cpu_single` | GpuSingleStageSim 持久化路径（默认 use_persistent=true） |
| `test_gpu_vs_cpu_multi` | GpuMultiStageSim 持久化路径 |
| `test_gpu_sim_batch` | SimBatch "vs individual runs" 一致性 |

### 新增/扩展测试

| 测试 | 内容 |
|------|------|
| `test_gpu_vs_cpu_single`（扩展） | D2 后增加：persistent vs fallback 结果一致性 CHECK（同一仿真配置，两次 run，结果相同） |
| `test_gpu_vs_cpu_multi`（扩展） | D3 后同上 |
| `test_gpu_sim_batch`（扩展） | persistent vs fallback 结果一致性对比 |
| `test_gpu_persistent_fallback` [新增] | `use_persistent=false` 显式降级所有三个类，验证 per-pair launch 仍正确 |
| `compute-sanitizer` | 所有 GPU 测试：mapped memory race/misaligned access 检查 |

### 精度验证（Aggressive 级别）

| 检查 | 方法 |
|------|------|
| 自一致性 | same params, Standard vs Aggressive, M/dM 相对误差 < 1e-4 |
| 参考一致性 | Aggressive GPU vs CPU Reference, 出口速度相对误差 < 0.1% |

---

## 七、边界条件与降级

| 场景 | 行为 |
|------|------|
| `use_persistent = false` | 使用原 per-pair kernel launch（行为完全不变） |
| `cudaHostAllocMapped` 失败 | 自动降级到 per-pair launch + stdout warning |
| GPU 不支持 mapped memory | 同上（cudaHostAllocMapped 返回错误） |
| 显存不足 | `cudaHostAlloc` 失败 → `std::runtime_error` |
| shutdown 信号 | kernel 检测 `*shutdown == 1` 后退出 while 循环 |
| `batch_id` 溢出 | `int` 计数器，2^31 步远大于仿真最大步数 (~10^5) |
| `doorbell` 死锁 | 超时未完成 → `cudaDeviceSynchronize` 捕获 kernel error |
| Aggressive 精度不足 | 默认 `GpuOptLevel::Full`，用户显式选择 Aggressive 即知晓 tradeoff |
