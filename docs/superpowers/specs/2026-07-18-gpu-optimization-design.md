# GPU 加速优化方向 — 第二阶段设计

> 日期: 2026-07-18 | 分支: `feature/cudaAcceleration` | 状态: 设计已确认

## 背景

第一阶段实现已完成：条件编译、GpuSingleStageSim、GpuMultiStageSim、SimBatch、API 文档、CMake 预设、19 套测试通过。当前性能：GPU RTX 5080 Laptop 仅比 16 核 CPU 快 1.2×（当前问题规模 S=2, N_fil=10）。

**核心瓶颈**：per-step 的 kernel launch overhead（20 次 launch/步 × 20000 步 = 400,000 次 launch）。

P0 batch kernel（grid 含 sim_id 维度）因 nvcc/nvlink 跨文件 struct-layout 问题阻塞。映射内存通路已验证可用——**持久化 kernel 方案可行**。

---

## 一、持久化 Kernel（主方向）

### 1.1 架构

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

### 1.2 每步流程

持久化 kernel 在一次 batch 内处理**单个仿真**的所有 (stage, filament) 对。SimBatch 外部循环遍历各仿真，逐次触发 kernel 批次。

```
Host side                              Device side (per block)
───────────                           ───────────────────────
for 每个 sim s:                        │
  1. 填充 seps[0..K-1] (sim s 的分    │
     离值，仅活跃 stage+fila 对)       │
  2. active_pairs = K                  │
  3. doorbell[0..K-1] = 1              │  while batch_id == last_id: spin
  4. batch_id++                        │  last_id = batch_id
                                       │  while doorbell[idx] == 0: spin
                                       │  if idx < active_pairs:
                                       │    读取 seps[idx]
                                       │    4D GL 积分 → M, dM
                                       │    __threadfence()
                                       │  results_M[idx]  = M
                                       │  results_dM[idx] = dM
                                       │  __threadfence()
                                       │  doorbell[idx] = 0
  5. for i in K: while doorbell[i]!=0: │
     spin                              │
  6. cudaStreamSynchronize             │
  7. 从 results 取值 → 填入            │
     batch_M1_(s, ...), batch_dM1_     │
──────────────────────────────────── ───────────────────────
结束所有 sim 后:
    每 sim LDLT → 力 → 运动 → 触发
    若终止: shutdown=1 → 下一大循环
```

**说明**：`K = n_active × N_fil`（当前步活跃的 stage 数 × 丝元数）。`N_b = n_stages × N_fil`（最大可能对数目）。闲置 block（idx ≥ K）同样参与 `doorbell[idx] = 0` 的同步周期，但不执行计算（`M=0, dM=0`）。

### 1.3 Kernel 内部

```cuda
__global__ void persistent_batch_kernel(
    const CoilGeo* coils, const FilGeo* fils,
    const double* gl_nodes, const double* gl_weights,
    volatile int* batch_id, volatile int* active_pairs,
    volatile int* shutdown,
    volatile double* seps, volatile int* doorbell,
    volatile double* out_M, volatile double* out_dM,
    int n_stages, int N_fil, int n_nodes) {

    int idx = blockIdx.x;
    int si = idx / N_fil;  int fi = idx % N_fil;
    int last_batch = -1;

    while (*shutdown == 0) {
        int cur = *batch_id;
        if (cur == last_batch) continue;
        last_batch = cur;

        while (doorbell[idx] == 0) {}

        double M = 0.0, dM = 0.0;
        if (idx < *active_pairs) {
            const auto& coil = coils[si];
            const auto& fil  = fils[fi];
            double ra_mid  = 0.5 * (coil.re + coil.ri);
            double ra_half = 0.5 * (coil.re - coil.ri);
            double rb_mid  = 0.5 * (fil.re + fil.ri);
            double rb_half = 0.5 * (fil.re - fil.ri);
            double la_half = 0.5 * coil.len;
            double lb_half = 0.5 * fil.len;
            double sep = seps[idx];

            // === 4D GL integration (same as single-pair kernel) ===
            int n2 = n_nodes * n_nodes, n4 = n2 * n2;
            double sM = 0.0, sdM = 0.0;
            // ... 512-thread shared-mem reduction ...
            M  = (coil.turns / 16.0) * sM;
            dM = (coil.turns / 16.0) * sdM;
        }

        out_M[idx]  = M;
        out_dM[idx] = dM;
        __threadfence();
        doorbell[idx] = 0;
    }
}
```

### 1.4 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 同步方式 | 逐 block 门控 + batch_id | 无 atomic 竞争，无 false sharing |
| block 数 | 固定 max_pairs（含闲置 block）| 启动一次不重新 launch |
| mapped memory | `cudaHostAllocMapped` | 绕过了 batch kernel 中 cudaMemcpy 的 nvcc 问题 |
| per-step 屏障 | `cudaStreamSynchronize` | 保证 host 可见 device 完成的写 |
| 降级开关 | `GpuBackend::use_persistent` | debug 时可回退到 per-sim 循环 |
| 核在同一 .cu 文件 | `sim_batch.cu` | 避免 P0 的跨文件 extern __global__ 问题 |

### 1.5 内存量

最坏情况 S=50, N_fil=500 → max_pairs=25000，单精度映射内存：
- doorbell: 25000 × 4 = 100KB
- seps: 25000 × 8 = 200KB
- results_M: 25000 × 8 = 200KB
- results_dM: 25000 × 8 = 200KB
- **总计: ~700KB**，远低于 16GB 限制

### 1.6 新增/修改文件

| 文件 | 变更 |
|---|---|
| `include/coilgun/simulation/cuda/sim_batch.hpp` | 新增 `PersistentBuffers` 结构，`init/free/launch` 方法 |
| `src/cuda/sim_batch.cu` | 添加 `persistent_batch_kernel`，替换 `compute_all_M1_dM1` 为 mapped 版本 |
| `include/coilgun/simulation/cuda/gpu_backend.hpp` | GpuBackend 新增 `bool use_persistent = true` |
| `tests/test_gpu_sim_batch.cpp` | 新增持久化模式测试（默认）和降级模式测试 |

---

## 二、其余优化方向（按优先级排列）

### P1：CUDA Stream 并发（优先级: 中）

**目标**：多仿真并发时 GPU kernel 和数据传输自动重叠。

**实现**：
- GpuMultiStageSim 构造时可选 `cudaStream_t`
- `compute_M1_dM1` 中 `<<<grid, block, 0, stream>>>`
- `cudaMemcpyAsync` 替代 `cudaMemcpy`
- 外部管理多个 GpuMultiStageSim 实例在不同 stream 上

**预期加速**：零额外 effort，多个 GPU 任务同时运行在 RTX 5080 上（84 SM 可同时执行多个 kernel）。

**风险**：低。stream 是 CUDA 标准特性。

---

### P2：FP32 混合精度（优先级: 中）

**目标**：利用 RTX 5080 的高 FP32 吞吐（~50 TFLOPS vs FP64 ~1.5 TFLOPS）。

**实现**：
- 4D 积分的 integrand 用 `float` 计算（椭圆积分计算的中间步骤）
- 最终归约累加器用 `double`（Kahan summation 或直接 double 累加）
- kernel 分两阶段：FP32 integrand → FP64 reduction

**预期加速**：3-5×（受限于 elliptic integral 内部仍有大量 double 运算）。

**风险**：中等。椭圆积分 Boost.Math 实现目前是 FP64，需要手写 FP32 版本或截断多项式逼近。

---

### P3：GL 节点 2D Block 并行（优先级: 低）

**目标**：提高 warp 占用率，减少 `__syncthreads` 开销。

**实现**：
- Block 维度从 `(512, 1)` 改为 `(9, 9)` = 81 threads
- `threadIdx.y` 映射外两层 (i1 → r_a, j1 → z_a)
- `threadIdx.x` 循环内两层 (i2, j2)，每次循环 81 次

**预期加速**：~1.2×（减少一次全局 __syncthreads）。

**风险**：低。但 81 threads/block 远低于 warp scheduler 最佳 occupancy。

---

### P4：GPU LDLT 求解（优先级: 低）

**目标**：消除 M1/dM1 结果的 PCIe 传输，直接在 GPU 上求解。

**实现**：
- cuSOLVER `cusolverDnDpotrf` + `cusolverDnDpotrs` 替代 Eigen LDLT
- 仅传输电流向量回 host（~12 doubles），而非 M1/dM1 矩阵（~S×F doubles）

**预期加速**：低分辨率下忽略不计（LDLT <1% 耗时），高分辨率 S>50, N_fil>500 时 2-3×。

**风险**：中等。需要 cuSOLVER 依赖。LDLT 用 Cholesky 替代要求正定矩阵（系统矩阵 [L-M] 是正定的）。

---

### P5：CUDA Graph（优先级: 低）

**目标**：一次性录制 + 无限重放，零 launch overhead。

**实现**：
- 录制一次完整的 step() 调用图（包含 kernel launch + cudaMemcpy）
- 使用 `cudaGraphInstantiate` + `cudaGraphLaunch` 重放
- 需固定 n_active（不支持提前终止分支）

**预期加速**：在当前规模下 ~1.1-1.2×（kernel launch 仅 5µs/步，相比持久化 kernel 的 0 仍有差距）。

**风险**：中等。不支持条件分支（需加 CUDA Graph conditional nodes，复杂度增加）。

---

### P6：增强型 Batch Kernel（优先级: 已取代） ✅

**状态**：被持久化 kernel **完整取代**。对比：

| 指标 | 原始 batch kernel (设计 §3.2) | 持久化 kernel (本文档主方向) |
|------|------|------|
| launch/步 | 1 次 `<<<>>>` | **0 次** |
| cudaMemcpy/步 | 1 次 H→D, 1 次 D→H | **0 次** (mapped memory 直访) |
| 数据通路 | kernel 参数 + cudaMemcpy | `cudaHostAllocMapped` load/store |
| nvcc 兼容性 | ❌ 跨文件 extern `__global__` + Eigen 导致 cudaMemcpy 失败 | ✅ 单 .cu 文件 kernel 定义，已验证 |
| 多 sim 支持 | grid=(num_sims, pairs) 单次内核 | 同一内核，`sim_id` 维度可轻松加入 |
| 实现状态 | 阻塞 | **待实现 (D1-D2)，设计已确认** |

**结论**：持久化 kernel 在功能上完全覆盖原始 batch kernel 设计目标，且消除了其实现阻塞问题。实现持久化 kernel 即视为 batch kernel 设计目标达成。

不再作为独立方向。

---

### P7：双缓冲 — CPU/GPU 重叠（优先级: 中）

**目标**：GPU 计算当前步 M1/dM1 的同时，CPU 并行执行上一步的 LDLT 求解。

**实现**：
- 两套 mapped 结果缓冲区 A/B 轮换
- 持久化 kernel 写入缓冲区 A → host 读取 B 做 LDLT
- 栅栏同步：host 等 GPU 完成 A，GPU 等 host 写完下一批 separations
- 可复用持久化 kernel 的 doorbell 机制，增加一个 `host_done` 标志

**预期加速**：1.3-1.5×。LDLT + 力学更新占 5-10% 总耗时（当前规模下）。高分辨率仿真下 LDLT 占比更高，收益更大。

**风险**：中等。双缓冲区使状态管理复杂化。

---

### P8：Cooperative Groups 归约（优先级: 低）

**目标**：减少 shared memory 同步开销。

**实现**：
- 替换 `for(int stride=total/2; stride>0; stride/=2){__syncthreads();}` 树形归约
- 使用 `cooperative_groups::tiled_partition<32>(this_thread_block())` 做 warp 级归约
- 9 次 `__syncthreads()` → 1 次跨 warp `shuffle_down` + 1 次 block-level `sync`

**预期加速**：~1.1×（归约开销本身 <1% 总计算时间，收益有限）。

**风险**：低。仅替换归约代码，不改变 kernel 结构。

---

## 三、实现阶段

| 阶段 | 内容 | 验证 |
|---|---|---|
| **D1** | 持久化 `persistent_batch_kernel` + mapped memory 集成 | 映射内存原型测试 |
| **D2** | SimBatch 集成：替换 compute_all_M1_dM1 | `test_gpu_sim_batch` 通过 |
| **D3** | GpuMultiStageSim 持久化模式（可选）| `test_gpu_vs_cpu_multi` 通过 |
| **D4** | 双缓冲 CPU/GPU 重叠 | 性能对比测试 |
| **D5** | CUDA Stream 并发 | 多仿真并发测试 |
| **D6** | 文档更新 | review |

---

## 四、边界条件与降级

| 场景 | 行为 |
|---|---|
| GPU 不支持 mapped memory | `GpuBackend::use_persistent = true` 构造时检测并 fallback 到 per-sim kernel 循环 |
| 显存不足 | `cudaHostAlloc` 失败 → `std::runtime_error` |
| shutdown 信号 | kernel 检测 `*shutdown == 1` 后退出 while 循环 → `cudaDeviceSynchronize` 等待退出 |
| `batch_id` 溢出 | `int` 计数器，2^31 步远大于仿真最大步数 (~10^5)，不处理 |
