# coligunCalc 性能审查记录

**审查范围：** 整个项目，包括 CPU/CUDA 物理实现、仿真编排、公共 API、构建配置、工具、测试和已有基准文档。

**审查目标：** 记录有代码或测量证据支持的性能优化机会；不以代码风格问题替代性能问题，不在本次任务中直接修改生产代码。

**终止条件：** 连续 5 轮独立复核没有新的可优化点。

## 状态

- 审查开始：2026-07-23
- 当前轮次：最终验证
- 连续无新发现：5 / 5
- 发现总数：52（CPU 16，CUDA 20，初始化/输出 3，构建/测试/基准/文档/CI/工具 11）

## 基线

- CPU Release 配置：成功。
- CPU Release 冷构建：44 个目标，约 8.13 秒。
- CPU Release 测试：20/20 通过，总耗时 52.40 秒。
- CUDA Release 构建：67 个编译/链接步骤成功；CUDA CTest 36/36 通过，总耗时 77.55 秒。
- CUDA Release `ctest -j 4`：36/36 通过，总耗时 31.88 秒，较串行约快 2.43 倍；该结果受机器线程和 OpenMP 资源影响。
- 主要测试耗时：`test_single_stage_sim` 约 22.28 秒，`test_multi_stage_sim` 约 31.33 秒。
- 本轮 GPU benchmark：RTX 5080 Laptop、CUDA Release；`medium-multi` Direct 5.710 ms（transfer 1.920 ms），`large-single` Direct 72.203 ms、Graph 13.624 ms，Persistent 682.086 ms 且 `gpu_executed=no`，batch 128 Direct 260.018 ms。

## 轮次记录

| 轮次 | 检查面 | 新发现 | 连续无新发现 |
|------|--------|--------|--------------|
| 0 | 基线建立 | 0 | 0 |
| 1 | CUDA kernel/传输/solver、wrapper、构建、测试反馈和 benchmark | 16 | 0 |
| A | 已知热点反向调用链、mask/active 生命周期、solver/事件/输出分配最终空白检查 | 0 | 1 |
| B | GPU/CPU 状态生产者到消费者、wrapper solver、thermal workspace 和 legacy adaptor | 5 | 0 |
| C | 几何元数据、公共 solver API、初始化和 summary 生命周期重复工作扫描 | 4 | 0 |
| D | geometry cache 生命周期、solver/thermal API 边界、summary/history、CUDA runtime 查询复核 | 9 | 0 |
| E | CPU 特殊函数/终止路径、CUDA kernel 并行度及反馈回路空白检查 | 4 | 0 |
| F | 全项目 CPU/CUDA/构建测试 benchmark 交叉去重空白检查 | 0 | 1 |
| G | CPU 数值/缓存、CUDA occupancy/内存、工程反馈链最终式空白检查 | 0 | 2 |
| H | CPU/CUDA 公共 API、错误路径、资源生命周期和工程验证最终空白检查 | 0 | 3 |
| I | CPU/CUDA/工程最终式逐项去重空白检查 | 1 | 0 |
| J | CPU/CUDA/工程最终空白复核 | 0 | 1 |
| K | CPU/CUDA/工程逐项终止前空白复核 | 0 | 2 |
| L | CPU/CUDA/工程严格去重空白复核 | 0 | 3 |
| M | CPU/CUDA/工程逐项终止前空白复核 | 0 | 4 |
| N | CPU/CUDA/工程最终终止前空白复核 | 0 | 5 |

## CPU 发现

### CPU-01 [P1] 合并互感与梯度的重复 4D 积分

**Symptom：** 单级和多级导数评估分别调用互感与梯度 coil-level API；两者都完整遍历 `n_nodes^4` 积分点并重复计算椭圆积分。位置：`src/simulation/single_stage_sim.cpp:153-162`、`src/simulation/multi_stage_sim.cpp:248-259`、`src/physics/mutual_inductance.cpp:184-222`。

**Source：** Ousterhout — Information Hiding；Code Complete — avoid redundant computation。

**Consequence：** `gprofng` 单级测试采样中 `ellint_1` 独占 54.65% CPU、`ellint_2` 独占 17.80%；两个 coil-level 函数包含时间分别为 40.72% 和 43.87%。

**Remedy：** 增加 fused kernel，一次积分同时得到 `M` 和 `dM/dz`，共享几何映射、`K/E` 和中间量；以数值回归和同机基准验证收益。

### CPU-02 [P1] 缓存 Gauss-Legendre 规则，避免热循环分配

**Symptom：** `integrate_4d` 每次调用 `gauss_legendre`，后者按值返回带两个 `std::vector<double>` 的结构。位置：`src/physics/mutual_inductance.cpp:144-145`、`src/physics/quadrature.cpp:17-177`。

**Source：** Fowler — Repeated Work；McConnell — Data Structure and Allocation Discipline。

**Consequence：** 每个 stage/filament 的每个导数评估都会重复构造节点与权重存储，产生不必要的小对象分配和初始化；尚需分配计数验证其占比。

**Remedy：** 使用进程内不可变静态数组和 `std::span`/固定容量视图；验证数值、线程安全和分配次数。

### CPU-03 [P2] 合并或自适应 OpenMP 并行区

**Symptom：** OpenMP `parallel for` 位于每个 stage 的导数计算内，每步/RK4 子步重复 fork/join。位置：`src/simulation/single_stage_sim.cpp:145`、`src/simulation/multi_stage_sim.cpp:230-234`。

**Source：** Hunt & Thomas — Orthogonality；Meszaros — Slow Tests/feedback loop。

**Consequence：** `gprofng` 中 `GOMP_parallel` 包含时间约占 4.33% CPU；低任务量时线程管理可能抵消并行收益。

**Remedy：** 展平 active stage/filament 任务或复用外层并行区，并按任务量自适应选择串行/并行；覆盖多 stage、filament 和线程数基准。

### CPU-04 [P2] 减少每步完整积分状态复制

**Symptom：** Euler 的 `step` 和 `make_trial` 连续复制 `IntegrationState`，其中 Eigen 向量和多态 excitation snapshot 重新分配/克隆。位置：`src/simulation/*_stage_sim.cpp` 的 `step`/`make_trial` 和 `src/simulation/integration_state.cpp:35-45`。

**Source：** Fowler — Repeated Work；Ousterhout — Strategic vs. Tactical Programming。

**Consequence：** 额外内存写入、堆分配和 cache 压力会随 filament 数和步数增长；当前 profile 说明它不是主热点，但仍是可测的次级优化点。

**Remedy：** 引入预分配 trial state；无事件 Euler 路径直接写入下一状态，RK4 复用存储，事件定位时才复制，并以 allocation/RSS/数值回归验证。

### CPU-05 [P2] 单级终止检查重复完整导数评估

`SingleStageSim::run` 在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:460) 每个 step 前调用 `check_termination`；达到 `velocity_decay_steps` 后，`:398-402` 又完整调用 `evaluate_derivatives`，尽管上一 `step()` 已计算当前物理导数。默认阈值后该额外互感/梯度遍历会放大 CPU-01 热点。建议缓存 acceleration/termination metric 或降低检查频率，并以导数调用数、wall 和最终状态验证。

### CPU-06 [P1] 每次导数评估重新组装并分解 dense 矩阵

`SingleStageSim` 在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:165) 每次导数评估重填矩阵并创建 `Eigen::LDLT`；多级对应路径在 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:263)。RK4 每步会执行四次。建议缓存不变 block 或评估 block/Schur 求解，按 `S/F=1/16、2/32、8/128` 分离 assembly/LDLT/solve 时间并回归残差。

### CPU-07 [P1] RK4 事件定位重复完整积分和线性求解

单级事件二分在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:343) 最多 64 次，每次运行完整 `advance_rk4_segment`，最终还会再算一次；多级在 `multi_stage_sim.cpp:523` 同样处理。单事件最坏约 66 个 segment、约 264 次导数评估。建议复用端点导数或使用 dense output/root interpolation，并以事件时间、最终状态、调用数和 wall 对比。

### CPU-08 [P2] RK4 cutoff 使用过期成员位置

`MultiStageSim::is_stage_within_range` 在 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:173) 读取成员 `integration_state_`，但 `evaluate_derivatives` 正在处理 s2/s3/s4 trial state。Full 优化级别可能错误地继续或跳过互感 work，同时影响 force。建议显式传入 trial position，在 cutoff 边界 RK4 case 上比较 active pair、wall、force 和最终速度。

## CUDA 发现

### CUDA-01 [P1] wrapper 重复 CPU 梯度积分

`GpuSingleStageSim::step` 在 [gpu_single_stage_sim.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_single_stage_sim.cu:228) 之后通过 `compute_force_at` 在每个 filament 调用 CPU `mutual_inductance_gradient_coil`（`:172-188`）；`GpuMultiStageSim::step` 在 [gpu_multi_stage_sim.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_multi_stage_sim.cu:504) 每步调用 `compute_pre_step_gradients`（`:284-313`）。引擎已经在 [gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:922) 下载本步 `dm1`，而 `SimBatch` 在 [sim_batch.cu](/mnt/data/Project/coligunCalc/src/cuda/sim_batch.cu:497) 直接复用边界梯度。该路径可能重新支付椭圆积分成本，当前 benchmark 未覆盖 wrapper；建议缓存初始边界梯度并复用引擎 `dm1`，以 wrapper/SimBatch 固定步数和数值回归验证。

### CUDA-02 [P1] 每步大量同步 H2D/D2H 往返

[gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:694) 逐项上传 masks、voltage 和 control 状态，`:907-955` 同步并逐项下载 status、互感、梯度、current、motion 和 derivative。本轮 `medium-multi` Direct 为 5.710 ms，其中 transfer 1.920 ms（约 33.6%）；`large-single` Direct 为 72.203 ms，其中 transfer 0.982 ms。建议把状态摘要合并为少量异步/固定 staging 回传，并让 force/status 尽量在设备侧归约。

### CUDA-03 [P1] batched solver 重复 D2D staging 和 status 同步

[gpu_solver.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_solver.cu:252) 每步把引擎 matrix/RHS 再复制进 solver workspace，`:273` 再复制 solution；`:292-299` 已下载 info/residual 并同步，[gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:907) 随后又调用 `validate_device_result`，在 [gpu_solver.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_solver.cu:340) 再次下载并同步。建议绑定 engine-owned workspace，并将验证延迟到单一 step barrier；覆盖 `D=9/32/129`、`B=1/8/128` 的 solver/transfer/wall benchmark。

### CUDA-04 [P1] Persistent backend 确定性 CPU fallback

[gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:554) 明确因缺少 dedicated control stream 抛异常，随后锁定 CPU fallback。本轮 `medium-multi` Persistent 为 682.086 ms，`gpu_executed=no`，而 Direct 为 5.710 ms。实现 control stream 协议前应在 capability/planner 中标为不可用，避免调用方误选；实现后再重新测量。

### CUDA-05 [P2] Direct/Graph 未按 workload 选择

[gpu_multi_stage_sim.hpp](/mnt/data/Project/coligunCalc/include/coilgun/simulation/cuda/gpu_multi_stage_sim.hpp:31) 默认 Direct；[gpu_execution_config.hpp](/mnt/data/Project/coligunCalc/include/coilgun/simulation/cuda/gpu_execution_config.hpp:117) 的 Auto 不读 timing，运行时 [gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:289) 仍选 Direct。本轮 `large-single` Direct 72.203 ms、Graph 13.624 ms，而 `small-single` Direct 6.035 ms、Graph 15.798 ms。应按 GPU、B/S/F/D、步数和 capture 成本建立阈值或 calibration，不能全局替换。

### CUDA-06 [P2] host state snapshot 和 boundary vector 每步复制/分配

`GpuEngine::step` 在 [gpu_engine.hpp](/mnt/data/Project/coligunCalc/include/coilgun/simulation/cuda/gpu_engine.hpp:375) 保存完整 `GpuEngineState`，pipeline 在 [gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:647) 再按值复制；每步还构造 active/compact-status vector。`SimBatch::configure_engine_boundary` 在 [sim_batch.cu](/mnt/data/Project/coligunCalc/src/cuda/sim_batch.cu:340) 创建多个 boundary vector，`:497` 复制完整 current，`:499` 复制 gradients。device-control 产生的 pair mask 下一步还会被 host active 重建覆盖并重新上传。大 batch/长仿真会累积 host bandwidth 和 allocator 成本。建议预分配固定容量并只保存回滚所需字段，需 allocation/RSS/B=1/8/128 验证。

### CUDA-07 [P3] kernel launch 重复查询已验证的 device buffers

[gpu_mutual_pipeline.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_mutual_pipeline.cu:141) 和 [gpu_state_kernels.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_state_kernels.cu:292) 在 launch 前重复调用 `cudaPointerGetAttributes`，部分路径还查询 device properties。建议保留公共 checked API，为 engine 内部提供初始化期校验后的 unchecked launch，并比较小 workload 的 runtime API 时间。

### CUDA-08 [P2] Graph key 过度包含运行时 mask

`graph_key()` 在 [gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:1049) 把 `stage_mask`/`mutual_stage_mask` 每个 bit 编入 key，导致 mask 变化时重新 capture；本机 8-step benchmark 出现 `small-single=1`、`medium-multi=2` 次 rebuild。若节点结构不变，应把 mask 作为 device input，只将真正改变拓扑的属性纳入 key，并用 mask 变化 benchmark 验证。

### CUDA-09 [P2] SimBatch 重算 GPU 已生成的 force

`SimBatch::run` 在 [sim_batch.cu](/mnt/data/Project/coligunCalc/src/cuda/sim_batch.cu:497) 复制 currents、`:499` 复制 gradients，并在 `:537-539` 做 `B*S*F` host force loop；[gpu_state_kernels.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_state_kernels.cu:11) 已执行 GPU force reduction，但 `d_force` 未作为 compact 结果回传。建议下载 per-batch/per-stage force 或提供 device-resident history，并测 B=1/8/32/128 的 transfer、host force 和 wall。

### CUDA-10 [P2] GpuMultiStageSim 执行无消费者的重复 force reduction

`GpuMultiStageSim::step` 在 [gpu_multi_stage_sim.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_multi_stage_sim.cu:517) 计算并写入 `applied_force_`，但该值无消费者；随后 `:535-539` 又为历史记录做 reduction。建议删除无效路径或让同一 device reduction 服务 applied/history，确认终止和历史数值不变后测 `S/F=2/32、8/128`。

### CUDA-11 [P2] stage completion 逐项触发 CUDA barrier

`GpuEngine::complete_stage` 在 [gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:183) 每次 `cudaMemsetAsync` 后立即同步并重新选择 graph variant；SimBatch 和 GpuMultiStageSim 在 `sim_batch.cu:526`、`gpu_multi_stage_sim.cu:528` 逐项调用。应收集 completion mask，一次性清零并只同步一次，用 B=1/8/128、S=1/2/8 API 计数和 wall 验证。

### CUDA-12 [P1] inactive batch row 仍支付完整 assembly/solve

固定容量 batch 不压缩 active rows；[gpu_state_kernels.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_state_kernels.cu:292) 仍按完整 B 组装，`GpuSolver::solve_device` 也按原始 B 调 cuBLAS。早结束比例高时，矩阵写入和 factorization 仍按原始 batch 增长。建议 active index list 或 active/inactive 分批处理，并测 inactive ratio、cuBLAS batch size 和 steps/s。

## 构建、测试和反馈回路发现

### BUILD-01 [P1] 76 MB 文本查表头拖慢增量编译

[lookup_table_data.hpp](/mnt/data/Project/coligunCalc/include/coilgun/physics/lookup_table_data.hpp:27) 为 76,415,075 bytes、583,348 行，约 2.9M 个 double，且被 `lookup_tables.cpp` 和 `self_inductance.cpp` 包含；object 当前约 23 MB，Ninja 日志记录单元编译约 10.0 秒（另一 CUDA build 记录约 33.3 秒）。建议改为紧凑二进制 blob/目标文件或校验后的运行时表，并比较增量构建与加载开销。

### TEST-01 [P1] 默认 CTest 串行且将长仿真标成 fast

[tests/CMakeLists.txt](/mnt/data/Project/coligunCalc/tests/CMakeLists.txt:18) 把通用测试统一标为 `fast`，但两个仿真测试约需 22.28/31.33 秒；[CMakePresets.json](/mnt/data/Project/coligunCalc/CMakePresets.json:91) 的 test preset 未配置 jobs。本轮串行 CUDA CTest 为 77.55 秒，`ctest -j 4` 为 31.88 秒，36/36 均通过。建议将长测试标为 slow/integration，并提供带资源策略的 quick/parallel preset。

### BENCH-01 [P2] GPU benchmark 只有单次冷态样本

[bench_gpu_engine.cu](/mnt/data/Project/coligunCalc/tests/bench_gpu_engine.cu:113) 和 `:135` 每个 workload 只运行一次，`:292-325` 没有 warm-up、重复次数、median/p95 或离散度；Graph capture 成本也混入首个样本。建议分离 cold/steady-state，提供 repetitions 和 percentile 输出，并区分 capture-inclusive/replay-only。

### DOC-01 [P2] benchmark 文档引用不存在的 preset

`AGENTS.md:21` 要求 `ninja-debug`/`ctest --preset debug`，中英文 [2026-07-19-unified-gpu-engine.md](/mnt/data/Project/coligunCalc/docs/benchmarks/2026-07-19-unified-gpu-engine.md:24) 又要求 `ninja-cuda-debug`，而当前只有 `cpu-debug`、`cpu-release`、`cuda-debug`、`cuda-release`。建议同步所有文档并添加命令 smoke check。

### CI-01 [P2] 性能 benchmark 未接入自动反馈

[src/cuda/CMakeLists.txt](/mnt/data/Project/coligunCalc/src/cuda/CMakeLists.txt:51) 将 benchmark 设为 `EXCLUDE_FROM_ALL`，仓库也没有 CI 配置或 benchmark 样本归档。固定 GPU runner 可用时，应增加非阻塞 performance job，记录设备/驱动/commit 和多次样本，以趋势或宽松阈值告警。

### BENCH-02 [P2] fallback 行仍输出数值 speedup

[bench_gpu_engine.cu](/mnt/data/Project/coligunCalc/tests/bench_gpu_engine.cu:262) 无条件计算并打印 `CPU/GPU speedup`，即使 `gpu_executed=false`。Persistent/Fallback 的 CPU 延迟可能被脚本误读为 GPU 性能。建议 fallback 输出 `N/A`，单独命名 fallback latency，并拒绝其作为 GPU baseline。

### TOOL-01 [P2] 查表校验仅抽样且未接入常规反馈

[verify_t_table.py](/mnt/data/Project/coligunCalc/scripts/verify_t_table.py:108) 默认随机抽样约 0.1% 表点并执行 Python 精确积分；CMake 只注册 generator worker-count 测试，没有固定边界点、元数据和 checksum 的快速校验。建议接入 CTest 快速子集，保留发布/夜间确定性抽样。

### BENCH-03 [P2] CPU/GPU 计时边界不对称

`benchmark_cpu` 在 [bench_gpu_engine.cu](/mnt/data/Project/coligunCalc/tests/bench_gpu_engine.cu:113) 把 `MultiStageSim` 构造纳入计时，而 `benchmark_gpu` 在 `:142-150` 先构造 `SimBatch`，只计 `batch.run()`。这会污染 speedup 对 setup/cold/steady-state 的解释。建议统一边界并分别报告构造、首步和 warm-up 后统计量。

### OUTPUT-01 [P2] 每步完整 history 输出并分配可变字段

`SingleStageSim::record_step` 在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:364) 每步 resize filament currents/temperatures；`MultiStageSim::record_step` 在 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:547) 还 resize cap/coil/stage fields。CUDA wrappers 在对应 `record_step` 中也逐步分配。`sample(every_n)` 只是运行结束后的复制，无法节省热循环的内存写入。建议提供记录步长/在线采样或连续存储，并用 allocation、RSS、wall 和结果一致性验证。

### CPU-09 [P2] 每次导数评估分配新的热态电阻和导数结果

`SingleStageSim::evaluate_derivatives` 在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:141) 返回新的 resistance 向量，并在 `:191-216` 组装动态 `DerivativeResult`；多级对应路径在 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:214) 和 `:312`。RK4 每步至少重复四次。thermal 或大 filament/stage 场景可通过预分配 workspace 复用这些输出，需以 allocation count、wall 和数值回归验证。

### CPU-10 [P2] 动生电势与力重复遍历同一梯度点积

单级在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:179) 计算动生电势点积，随后 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:201) 再调用 force 点积；多级对应 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:289) 和 `:322`。可在一次 active stage/filament 遍历中同时生成 RHS 所需 motional sum 与 force，并做轨迹回归。

### CUDA-13 [P1] GPU wrapper 固定 Eigen solver，Direct 路径每步退回 CPU 求解

`GpuSingleStageSim::make_config` 在 [gpu_single_stage_sim.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_single_stage_sim.cu:101) 与多级 [gpu_multi_stage_sim.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_multi_stage_sim.cu:111) 固定 `SolverMode::Eigen`。因此 [gpu_engine.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_engine.cu:852) 的 Direct 分支每步下载矩阵/RHS、执行 host Eigen LDLT、再上传 derivative。现有 SimBatch benchmark 的 Auto 路径解析为 cusolver，不能覆盖这个 wrapper 成本；应让 wrapper 按 workload 使用 Batched 或显式报告混合路径，并分别测 solver/transfer/wall。

### CUDA-14 [P1] 一次性 thermal helper 每次调用都重建 CUDA workspace

`update_thermal_batch` 在 [gpu_thermal.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_thermal.cu:363) 每次创建 `ThermalWorkspace`；其 `update` 在 `:192-210` 重新分配/上传输入、下载四个输出并同步，还在 `:182-190` 做设备属性查询和 BF 级校验。按时间步调用时，应该复用长期 workspace、分离不变量上传，并尽量只回传摘要；需以重复步数和 B*F 测 transfer、allocation、wall 和 CPU 数值一致性。

### CUDA-15 [P3] legacy GpuAdaptor 的 per-step API 固定同步并重分配输出

`GpuAdaptor` 在 [gpu_adaptor.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_adaptor.cu:207) 和 `:235` 的 upload/download API 使用同步 `cudaMemcpy`，下载还 resize 两个结果 vector。legacy/persistent 使用者无法复用异步 device output，批量 mutual 结果不能继续留在设备上。可增加异步 view/事件接口并复用 caller-owned capacity；当前仓库主要保留该路径给 persistent，因此标为条件性 P3。

### CPU-11 [P3] 导数热循环重复解析不变的 filament 几何元数据

`SingleStageSim::evaluate_derivatives` 在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:143) 每个 filament 重新计算 axial/radial 索引、相对轴向位置、filament length 和半径 accessor；多级 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:228) 对每个 active stage 重复同一元数据。构造期预计算连续 geometry/index 数组可减少每个 Euler/RK4 子步的除法、函数调用和 cache miss，需做小 workload/cutoff/RK4 数值与 wall 回归。

### CUDA-16 [P2] 公共 batched solver API 每次调用都做 host 转置和完整同步

`GpuSolver::solve_batch` 在 [gpu_solver.cu](/mnt/data/Project/coligunCalc/src/cuda/gpu_solver.cu:198) 将每个 `B*D*D` host row-major 矩阵转为 column-major，再在 `:202-224` 执行 H2D、factor/solve 状态回传、同步和 solution D2H，`:226` 继续 host residual 检查。该路径不同于 `solve_device` 的 D2D staging；应提供 device-view/column-major workspace 和延迟验证接口，并以 B/D 分解 transpose/transfer/sync/wall。

### INIT-01 [P2] 相同几何的多次 simulation 构造重复计算不变互感矩阵

`SingleStageSim` 构造在 [single_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/single_stage_sim.cpp:68) 每次重建 filament mutual，`MultiStageSim` 在 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:100) 还重建 stage mutual；GPU wrapper 传入空数组时，[gpu_engine.hpp](/mnt/data/Project/coligunCalc/include/coilgun/simulation/cuda/gpu_engine.hpp:228) 也为每个 engine 补算。对相同几何的参数扫描可共享 keyed immutable blocks，需测重复构造的 constructor wall、命中率、RSS 和数值。

### OUTPUT-02 [P3] 多级 summary 按 stage 重复扫描完整 history

`MultiStageSim::prepare_summary` 在 [multi_stage_sim.cpp](/mnt/data/Project/coligunCalc/src/simulation/multi_stage_sim.cpp:633) 对每个 stage 扫描全部 history，随后 `:655` 再做全局扫描；GPU multi 和 `SimBatch` 也有对应 per-stage pass。可在一次 history pass 或 record_step counters 中维护 peak/current/active 指标，按 S、B、H 比较 run 收尾 wall、cache 和 summary 数值。

### CPU-12 [P2] 规则 filament 网格重复求值相同的 filament-filament mutual

规则网格中 mutual 只取决于 radial pair 和 axial index 差，但 CPU 构造路径仍按全部 filament pair 调用积分 API，且关闭 cache（`single_stage_sim.cpp:93-99`、`multi_stage_sim.cpp:141-147`）。可按唯一几何组合预计算再映射回矩阵；需要对不同 radial/axial 配置测构造 wall 并回归矩阵和 force。

### CPU-13 [P3] inactive stage 导数评估仍全量清零 mutual workspace

多级导数评估每次先对完整 `S x F` mutual/gradient workspace `setZero()`，之后只计算 active/in-range stage（`multi_stage_sim.cpp:221-259`）。应按 active rows 或 generation 标记减少清零流量，覆盖 active ratio、Euler/RK4 和数值回归。

### CUDA-17 [P2] 非 GPU thermal 路径仍执行无效 pre-step D2D snapshot

`GpuEngine` 每步复制完整 current state 到 `d_pre_step_currents`，但该 buffer 只被 GPU thermal 使用（`gpu_engine.cu:762-804,877-892`）。thermal disabled/CPU thermal 仍支付 D2D bandwidth 和依赖；应仅在 GPU thermal 时捕获，并用 B/D、thermal mode 的 CUDA trace 和 steps/s 验证。

### CUDA-18 [P3] `GpuAdaptor::setup_batch` 额外保留单实例结果缓冲

`setup_batch()` 先调用 `setup()`，batch-only 使用也分配单实例 separation/result buffers，增加 3 次 `cudaMalloc` 和显存占用（`gpu_adaptor.cu:165-200`）。建议拆分或懒分配，并验证 setup 后单实例 API 兼容性。

### CUDA-19 [P3] resident thermal 初始化先零填充再覆盖不变量

`ThermalWorkspace::initialize()` 上传 `i/m/r0/material` 的临时零 vector，GPU engine 随后立即覆盖 `m/r0/material`，resident kernel 也不使用 `i`（`gpu_thermal.cu:142-163`、`gpu_engine.cu:544-551`）。建议按需分配或直接初始化 resident state，记录 H2D 和初始化 wall。

### CUDA-20 [P2] device control/status 尾部 kernel 在 batch 内串行扫描

device-control 每个 batch 线程串行遍历 `S`、`S*F`，compact-status 每个 batch 线程串行扫描 `D`，尽管 launch 使用 128-thread blocks（`gpu_state_kernels.cu:174-240,348-373`）。较大 S/F/D、B=1/小 batch 时可能成为固定尾部延迟；建议 block 内并行/reduction 并用 CUDA event/Nsight 测占比。

### DOC-02 [P2] Graph 性能范围的 API/实现文档口径不一致

实现中的 Graph capture 已覆盖 assembly、solver、state update、thermal、device control 和 compact status，但 `GpuMultiStageSim` 头文件仍称只捕获 mutual-inductance segment，API 主文档又称完整 step。应统一中英文 API、benchmark 记录和 execution report 的 Graph scope，避免 speedup 归因和优化优先级错误。

### BENCH-04 [P2] ExecutionReport 计时字段重叠且缺少控制/吞吐指标

`gpu_time_ms` 是完整 CUDA pipeline 的 host wall，solver/transfer/thermal 是其内部累计值，现有 benchmark 没有 control/sync、capture/replay、`steps/s` 或 `simulations/s`。应明确嵌套关系或改成互斥阶段计时，并补充吞吐指标。

### BUILD-02 [P2] 所有 preset 默认构建完整测试图

四个 CPU/CUDA preset 都打开 `COILGUN_BUILD_TESTS`，测试目标还各自编译 `test_main.cpp`。library-only、发布和性能迭代仍支付测试编译/链接成本。建议提供显式 tests-on 和 library-only preset，并用 ON/OFF clean/incremental build 对比验证。

### CPU-14 [P3] Struve 交叉区间同时计算幂级数和渐近分支

`struve_h0`/`struve_h1` 在 `8 <= x < 20` 先计算幂级数，再计算渐近分支并比较结果，只保留一个值（`src/physics/struve.cpp:133-161`）。exact self-inductance 输入若频繁落在交叉区间会支付两套特殊函数计算。应通过离线误差扫描缩窄区间或选择单一主分支，并用精度回归和构造 wall 验证。

### CPU-15 [P3] 多级终止判定重复扫描 stage 和 velocity history

`MultiStageSim::check_termination` 先调用 `check_all_finished()` 扫描 stages，随后又扫描 triggered/circuit 状态；达到 decay window 后还读取完整 velocity history 窗口（`multi_stage_sim.cpp:586-618`）。可维护 unfinished/active counters 和 rolling decay state，比较终止检查 wall 与停止条件回归。

### CPU-16 [P2] Struve 渐近路径重复构造 Gauss-Laguerre 规则

`k0_integral`/`k1_integral` 每次调用 `gauss_laguerre(30)`，后者按值构造节点和权重 vector（`src/physics/struve.cpp:78-96`、`src/physics/quadrature.cpp:220-287`）。exact self-inductance 的渐近 Struve 路径会反复支付该分配和填充成本；应缓存 15/30 阶 Laguerre 表，并以精度、分配计数和构造 wall 验证。

### CUDA-21 [P2] assembly kernel 的 RHS 动生项按行串行归约

`assembly_kernel` 只有 `column == 0` 的线程计算 RHS，stage row 串行遍历 `F`，filament row 串行遍历 `S` 并重复读取 gradient（`gpu_state_kernels.cu:119-157`）。`B` 小、`S/F` 大时可能形成 assembly 尾部瓶颈；建议 block 内归约并以 CUDA event 分离 RHS/assembly 时间，回归 residual 和数值。

### CUDA-22 [P3] state update kernel 每个 batch 串行更新全部 current 维度

state update 每个 batch 线程循环更新全部 `D` 个 current，launch 虽使用 256-thread block（`gpu_state_kernels.cu:43-72,406-431`）。`B=1/8`、D 较大时线程利用率不足；建议 batch x dimension 映射并比较 kernel 时间、吞吐和数值一致性。

## 否决项与限制

- `large-single` 的 Direct/Graph、Persistent fallback、medium transfer 和 CTest `-j 4` 已有本机实测，但单次 benchmark 受频率/环境影响，不能当作跨机器保证。
- wrapper 梯度重算、solver staging、host snapshot、pointer validation 和 history storage 的具体收益尚未拆分出独立微基准，当前记录为有代码路径支持的待验证候选。
- 未运行 `compute-sanitizer`；本次任务未修改生产代码，因此没有代码回归需要修复。
- 复核 B 的并行子代理首次和重启均按用户指定请求了 `gpt-5.6-luna` High，但服务端持续返回 403，故本轮结论由本地源码和定向测试支持，未改用 Terra/Sol。

## 验证记录

- CPU Release 构建：`cmake --build --preset cpu-release --parallel 2` 退出 0。
- CUDA Release 构建：`cmake --build --preset cuda-release --parallel 2` 退出 0。
- CPU Release 全量 CTest：20/20 通过，总耗时 56.46 秒。
- CUDA Release 全量 CTest：36/36 通过，总耗时 83.08 秒。
- 最终 `bench_gpu_engine`：GPU 可用；`medium-multi` Direct 6.099 ms（transfer 2.252 ms），`large-single` Direct 72.657 ms、Graph 14.666 ms，batch-128 Direct 261.396 ms；Persistent 为 CPU fallback 且 `gpu_executed=no`。
- 发现标题计数 52，连续 5 轮空白复核完成；审查文档行尾空白检查和 `git diff --check` 通过。
- 审查期间的轮次明细属于临时过程记录，完成汇总后已清理；本文件保留最终发现、证据和验证结果。

## B9 post-B8 多轮复核

B8 完成后重新执行了七轮六领域全项目 review（CPU physics、CUDA runtime、
simulation、architecture/debt、test/benchmark feedback、cross-check）。Round
2 新增一项 `ENGINE-01`：`GpuEngine::execute_cpu_physical_pipeline()` 的 CPU
fallback assembly 在 `gpu_engine.hpp:936-945` 仍分开调用 mutual/gradient
wrapper，未接入 fused pair。Round 3-7 连续五轮 `new_findings_count=0`，满足
终止条件；该项保留为未来独立 Task，B9 本身只读。

轮次报告和最终摘要属于已清理的临时过程记录；本节和本文件即为保留的最终
审查摘要。subagent 按用户指定的
`gpt-5.6-luna` High、Fast 关闭、禁止 Terra/Sol 尝试启动，但 controller 的
six-agent wave 遇到 thread limit，单 agent 重试收到 HTTP 403；没有使用其他
模型，最终结论由主 agent 的源码、测试和 benchmark 证据支持。
