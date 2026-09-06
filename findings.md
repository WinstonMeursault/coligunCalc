# CUDA 优化实施发现

## 2026-07-28 预检

- 当前 checkout 是普通仓库，不是 linked worktree。
- 当前分支已包含设计与执行计划三个文档提交，尚无实现 Task ledger。
- `docs/2026-07-24-GPU-and-Repository-Modernization-Assessment.md` 是既有未跟踪最终产物，必须保留。
- 计划中的 Task 并行是可选项；按 SDD 工作流将 implementer 串行派发，GPU 验证始终串行。
- 当前环境支持 `gpt-5.6-sol` 并可设置 High reasoning；不使用 Fast。

## 基线

- 隔离 worktree：`.worktrees/b0-t1`，分支 `task/cuda-opt-b0-t1`，base `89e9f36`。
- CPU Release build 成功，CTest `21/21` 通过，29.88 秒。
- CUDA Release build 成功，CTest `37/37` 通过，52.17 秒。
- CUDA 首次构建存在既有 `BOOST_MATH_ENABLE_CUDA` 重定义和 Eigen constexpr host/device warning；不是本轮回归。
- 首次并行 configure 的嵌套 session 没有被控制器正确等待，导致一次 CPU build 过早启动；串行恢复后基线正常。

## B0-T1 Review Context

- 当前 `bench_gpu_engine` 已实现 3 次独立 repeat，以及 setup、first-step/capture-inclusive、replay-only、warm-up、steady-state 阶段。
- 当前固定 workload 只有 `baseline` 和 `thermal`，尚未覆盖计划要求的 `small-single`、`medium-multi`、`large-single`、batch 128 矩阵。
- `ExecutionReport::gpu_time_ms` 是包含 host orchestration 与 transfer 的 CUDA-backed pipeline host wall；solver/thermal/transfer 是其内含字段，不能相加为总时间。
- 当前 schema 没有独立的 mask-update、capture/replay、control/status、sync 字段；Graph 仅通过 rebuild delta/total 间接观察。
- 当前 fallback row 已输出 `n/a` speedup，并有 `gpu_executed`/finite/fallback reason，必须保持。

## B0-T1 Implemented Evidence

- Task commit `35b3b16` 仅修改允许写集：benchmark、schema 文档和 standalone schema checker。
- 新 schema `gpu-benchmark/v2` 增加 measurement window、capture/replay、execution kind、active ratio、mask updates、simulations/s、stage timer/residual 占位和明确 timing relation。
- 当前生产 `ExecutionReport` 未暴露的 stage timer/residual 使用 `n/a`，没有伪造零值。
- 三次 final artifact 每次 165 GPU rows：165/165 finite、108 GPU execution、24 fallback、33 setup；fallback numeric speedup 为 0 行。
- Task 分支单一 commit、worktree clean、diff check 通过；等待独立 task reviewer verdict。

## B0-T1 Review Findings

- Important: runtime-mask GPU rows执行的工作与 CPU baseline 不同，却仍输出 numeric speedup；必须输出 `n/a` 并由 checker 强制。
- Important: checker 只用全局集合/布尔值验证矩阵，未要求 165 个唯一 composite rows、每个 request 的 3 repeats×5 windows，也未完整类型检查 numeric/enumerated fields。
- Minor: capture/replay 分类基于 requested backend，应该基于 resolved backend。
- Reviewer 无法从 diff 验证 raw artifacts、测试输出、RED/GREEN chronology、worktree clean 和 commit body；controller 将从 report、artifact 和 git 直接核验。

## B0-T1 Final Resolution

- runtime-mask speedup、165 个唯一 composite rows、完整类型/枚举校验和 resolved-backend capture/replay 分类均已修复。
- 最终 Task commit `f0f10a7` review clean，并以 `9e50677` 集成到执行分支。
- 三份最终 raw artifact 位于 `/tmp/coligunCalc-b0-t1-f0f10a7-run-{1,2,3}.md`，在 B0-T2 完成前保留。
- 后续所有新 subagent 必须显式使用 `gpt-5.6-sol`、`reasoning_effort=high`；不设置 Fast 或其他 service-tier 覆盖。

## B0-T2 Baseline Preflight

- worktree `.worktrees/b0-t2` 位于 base `9e50677bf04b...`，分支 `task/cuda-opt-b0-t2`。
- CPU Release 21/21、CUDA quick 30/30、CUDA integration 8/8，均为 0 failure。
- CUDA 构建 warning 与 B0-T1 记录一致，不构成新增失败。

## B0-T2 Final Baseline

- Final Task commit `ddc8ea0` 仅新增中文 baseline 文档，review clean。
- 当前-base 三份 artifact 各 165 GPU rows，合计 495/495 finite；324 GPU execution、99 setup、72 fallback。
- 每个精确请求有 9 个 steady-state 样本；文档保留 median、nearest-rank p95、population stddev 和原始 latency 样本。
- single-stage CPU/GPU velocity 均为 `0.00592475`；multi-stage 均为 `3.9989223864947752`；device solver residual 为 0（阈值 `<1e-12`）。
- 基线存在明显双峰机器状态方差；后续 before/after 必须交错运行并保留原始样本。

## B1-T1 Performance Evidence

- B1-T1 的精确受影响区段 microbenchmark 覆盖 `D={9,32,129}`、`B={1,8,32,128}`，
  每个 shape 进行 10 次 warm-up，并按 D 分别测量 300/120/20 次；base/after
  各 3 轮串行交错。
- 12 个 shape 的 `solve_device()` enqueue median 从 `45.1--1953.2 us` 降至
  `9.5--39.9 us`。包含 validation 的 total median：D=9 改善 20.5%--21.8%，
  D=32 改善 7.0%--7.6%，D=129 改善 0.5%--0.8%。最大 residual
  `1.97372982156e-16`，最大 solution error 为 0，所有 status 成功。
- B0 fixed workload 追加 base/after 各 5 份 artifact（每份 204 行 benchmark 输出，
  临时落盘因终止空行 `wc -l` 为 205；165 GPU rows，schema 全部通过）。设备状态存在与 B0 相同的双模：未受改动的
  small-single/direct/Eigen 控制行在 `0.554--1.021 ms` 间跳变；该幅度远高于
  单次 barrier 优化。15 个受影响 steady-state 样本的 raw median 变化为
  `large/direct +1.7%`、`large/graph -3.8%`、`medium/graph +0.8%`、
  `medium/graph/runtime-mask +0.1%`、`thermal 0.0%`，p95 仍受双模污染。
- B1-T1 结论：单一 validation barrier 的局部收益和正确性成立，但 fixed workload
  端到端证据没有达到 5% 门槛；本 Task 不改变默认 solver policy，是否接入 wrapper
  仍由 B1-T2 的独立 workload 证据决定。

## B1-T1 Final Integration

- Task commit `f0cebf5` 已集成为 `920778c`；集成分支定向增量构建 exit 0。
- 集成分支 focused CUDA tests：solver、resident pipeline、assembly、resource contracts 共 4/4 通过。
- B1-T2 必须单独证明 wrapper policy 接入后的数值正确性和性能门禁；不得把 B1-T1 的 microbenchmark 收益直接视为 wrapper 端到端收益。

## B1-T2 Preflight

- `GpuSingleStageSim` 与 `GpuMultiStageSim` 已在既有提交 `3ddd86e` 将 wrapper solver request
  从固定 Eigen 改为 Auto；large dimension 测试已验证 resolved Batched。B1-T2 不应重复
  修改 planner 阈值或硬编码 Batched。
- Single-stage 的 `check_termination()` 仍调用 `compute_force()`，后者逐 filament 在 CPU
  重算 `mutual_inductance_gradient_coil`；这是明确的 wrapper 热路径重复计算。
- Multi-stage applied/recorded force 已读取 `engine_->state().dm1`，没有同类 CPU mutual
  gradient 重算。它的 pre-step/post-step force 时序和 completion boundary 受现有测试保护。
- B1-T2 worktree focused baseline 3/3 通过；既有 CUDA warning 与前序基线一致。

## B1-T2 Implementation Evidence

- Single-stage velocity termination 旧路径每次按 filament 调用 CPU 9-point mutual gradient；
  60-step wrapper workload 中 Direct/Graph 分别从约 `141/146 ms` 降到 `30/40 ms`，
  对应 median 改善 `78.639%/72.799%`。关闭 velocity check 的 Direct control 仅变化
  `0.008%`，归因明确。
- small Direct wrapper 保持 Auto -> Eigen，small Graph 和 large Direct 为 Auto -> Batched；
  全部 benchmark 行 `gpu_executed=1`、无 fallback。
- 不受影响的 single/multi/thermal/large control 无超过 2% 负回退；multi thermal 的
  正向 6.889% 属机器状态观测，不归因于本 patch。
- 一对 Direct velocity 最终值差 `7.482404140723986e-9`，其余配对精确相等；focused
  与 integration 在既有 tolerance 下通过。
- B1-T2 review 补齐 empty-history committed-force guard：fresh run 在 acceleration
  termination 生效前至少提交一个 engine step。该 fix 不进入已提交 history 后的 benchmark
  热路径；final commit `d9106a8` review clean，0 Critical/Important/Minor。

## B1 Gate

- 集成 HEAD `3d254b6`：B1-T1 单一 validation barrier 与 B1-T2 committed wrapper force
  均已集成，两个 Task 都是唯一 final commit 且 review clean。
- 集成验证：focused wrapper/path 3/3、CUDA integration 8/8、CUDA quick 30/30。
- Solver microbenchmark 保留 D=9 约 21%、D=32 约 7%、D=129 约 0.5% 的 total
  改善；未把 D=129 或噪声端到端结果扩大为默认策略收益。
- Wrapper velocity-check Direct/Graph median 改善 `78.639%/72.799%`；disabled control
  `+0.008%`，归因明确。其余关键 workload 无超过 2% 负回退。
- B1 Gate 结论：通过。现有 Auto/Eigen/Batched/fallback policy 保留；不放宽 tolerance、
  不降低 precision、不修改公共 API。

## B2-T1 Preflight

- worktree `.worktrees/b2-t1` 位于 base `3d254b6`，分支 `task/cuda-opt-b2-t1`；focused baseline 4/4 通过。
- `GpuEngine::execute_physical_pipeline()` 当前每步执行 host `B*S*F` 三重循环，并把完整 `active_pairs` 复制到 device；该固定成本正是 B2-T1 的优化对象。
- `device_control_kernel()` 已在步尾基于 device-resident `active_mask`、`trigger_mask`、`mutual_stage_mask` 更新 `pair_active`。随后 host 在下一步开始时完整覆盖 `PairActive`，破坏了 device residency。
- host 当前 pair 表达式为 `active_mask && mutual_stage_mask`，device-control 表达式还包含 `trigger_mask`。实现必须先确认 mutual pipeline 的既有语义，不能让尚未 trigger 的 stage 是否参与 mutual 计算发生意外变化。
- 推荐的最小方向是在 mutual 执行前调用独立的 device mask-to-pair kernel；初始 host authoritative mask 仍需正确到达 device，device control 在步尾生成的变化不得在下一步被陈旧 host snapshot 覆盖。
- B2-T1 RED 应直接证明旧路径缺少“从当前 device masks 派生 pair activity”的契约，或证明 device-control 生成的 pair activity 会被下一步 host staging 覆盖；测试必须先在旧生产代码上按预期失败。
- `launch_device_step` 的 Graph/Direct 共用 capture body 当前顺序为 separation -> mutual -> assembly -> solve -> state -> device control；若增加 mask-to-pair kernel，必须在每次 mutual 前且位于 Graph capture body 内，不能只修非 Graph fallback 分支。
- device-control 路径在步尾把 active/trigger/stage/mutual/completion masks D2H 到 host，下一步又全部 H2D；B2-T1 review 应区分“host 确实修改了输入”与“host 只是镜像 device 状态”，防止仅删除 `active_pairs` 一份上传却继续无条件覆盖 device-authoritative masks。
- 非 Batched/Eigen fallback 路径有独立的 mutual launch，必须同样派生正确 pair activity；否则 Direct/Graph/Batched 可能出现语义分叉。
- 进一步核对 `GpuEngine::assemble_physical_system()` 的 canonical CPU fallback：pair 仅在 `mutual_stage_mask && trigger_mask` 时计算并写入矩阵；batch inactive 在外层跳过。因此 device-control 的表达式是正确契约，旧 CUDA host staging 遗漏 `trigger_mask` 本身是跨后端语义缺口，B2-T1 RED 应锁定该行为。
- 首轮自定义 active-matrix 在 100% active 的 B32/B128 显示约 47% 改善，远超本 Task 可解释范围；根 CUDA 静态库可能未包含已集成的 B1 solver barrier 改动。该 base artifact 暂时判为无效，必须重建 `3d254b6` 后重新链接/重测，不能用于性能结论。

## B2-T1 Implementation Evidence

- canonical CPU fallback 已确认 pair 契约为 `active batch && trigger_mask && mutual_stage_mask`；旧 CUDA host staging 遗漏 trigger，不只是性能冗余，也造成跨后端语义不一致。
- 最终选择把 pair activity 派生融合进既有 separation kernel，避免独立 mask kernel 的额外 launch；Graph、Direct/Batched 和 Eigen fallback 共用该 staging。
- 删除每步 host `B*S*F` active-pair 重建、完整 PairActive H2D，以及 device-control 步尾重复 `B*S*F` pair 写；没有新增同步或 precision/tolerance 变化。
- 有效 active-ratio matrix 使用 20 warm-up、20 measured、3 repeats，base/after 各 180 行；全部 mapping/finite/GPU execution 通过、0 fallback。最差 group median 为 B128/75%/low 的 `-0.7651%`，低于 2% 回退门禁。
- 各 batch 的 group-median 变化：B1 `+0.6405%`、B8 `+0.1325%`、B32 `+0.0683%`、B128 `-0.1482%`。没有证据支持当前引入 shape-specific fallback。
- 固定 workload 粗筛的 small Direct/Graph 典型 median 分别改善约 6.8%/5.8%，large Direct 约 2.5%；B128 Direct 基本持平，thermal Graph 无明确回退。
- 变化范围上传需要 dirty-generation lifecycle 和 Graph 边界 metadata，并与 B2-T3 重叠；本 Task 保留其他 host mask boundary，后续统一评估，不扩大 B2-T1 职责。
- B2-T1 controller review 0 Critical/Important、1 Minor：`StepWorkspace::active_pairs` 及其构造期 resize 已失去消费者，后续 lifecycle Task 应删除，避免保留无用 host 分配和失效状态。
- B2-T1 已集成为 `21c5332`；集成 build、focused 4/4 和 integration 8/8 通过。后续 B2-T2 可依赖 device-derived pair contract，不得重新引入 host PairActive staging。

## B5 CUDA Review Addendum

- R25 实现复核准备：原 `complete_stage()` 对每个完成 stage 立即执行 `cudaMemsetAsync` 并同步，导致多 stage completion 产生重复 host barrier。当前候选将位置收集到 host pending queue，在下一次 execution stream pipeline 开始处入队清零，并由后续统一 pipeline synchronization 提交结果。
- R25 性能取舍：没有为稀疏 completion 新增 device index buffer、H2D 上传和独立 clear kernel；这些固定开销可能大于少量同 stream `cudaMemsetAsync`。当前需重点审查 queue 去重、reset/fallback/step 异常边界、Graph capture 外部清零顺序以及 API 文档语义。

- R20 clean：审查 execution policy 的请求/解析/fallback 语义、`BackendSelectionReason` 与 `ExecutionReport` 合并规则、checked/trusted launch 边界和 resource-contract 覆盖。未发现新的正确性、生命周期或可安全合入的 CUDA 性能问题；当前工作树 CUDA Release 增量构建成功，定向测试 3/3 通过。
- R21 clean：审查 Graph topology key、runtime mask replay、capture/replay 异常锁定、resident allocation reuse/reset，以及 policy/report 交互。未发现新的正确性、生命周期或可安全合入的 CUDA 性能问题；定向测试 4/4 通过。
- R20/R21 均为修复后连续 clean，clean streak 正式更新为 `4/8`；尚不能宣称八轮 Review 完成。
- R22 clean：审查 physical pipeline、thermal workspace、reset/init、D2D snapshot 和所有 host/device boundary 的 execution-stream 归属，确认 D2H 前置 barrier 与异步源缓冲区生命周期成立。未发现新的正确性、生命周期或可安全合入的性能问题；single/thermal/resident 定向测试通过，resident memcheck 为 0 errors，clean streak 为 `5/8`。
- R23 clean：审查 mutual/state/assembly/control/status kernel 的边界与索引、尺寸溢出、shared-memory/线程限制、precision 分支及 checked/trusted launch 复用。未发现新的正确性、资源或可安全合入的性能问题；定向 CUDA 测试、precision、memcheck 和 `cuobjdump` 资源检查均通过，clean streak 为 `6/8`。
- R24 clean：审查 device solver 的 active compaction/scatter、pointer table、cuBLAS stream、Graph fixed-shape solver、mask-aware residual 和 status barrier。未发现新的正确性、生命周期或可安全合入的性能问题；solver/Graph/batch 回归、memcheck 和 racecheck 均通过，clean streak 为 `7/8`。
## B5 Final Findings

- 同环境 benchmark 采用 baseline `9e50677` 与 after `0adf191`，均为 CUDA Release、RTX 5080 Laptop、CUDA 13.3、driver 610.57.04；每个提交 3 份 artifact，每份 165 GPU rows，schema 全部 PASS。
- after 共 495 GPU rows：实际 GPU 324、fallback 72、setup 99；实际 GPU rows 全部 finite。fallback 不进入 speedup。
- steady-state median 变化：small Direct -2.37%、small Graph -1.25%、medium Direct -7.88%、medium Graph -4.71%、runtime-mask Graph -1.02%、large Direct -2.70%、large Graph -0.49%、batch128 -0.06%、GPU thermal -0.46%。统一 5% 主线门禁未通过，但无超过 2% 的关键 workload 未解释回退。
- 固定 whole-branch package `3ddd86e..0adf191` 的 controller review 结论为 0 Critical、0 Important、0 Minor。Luna reviewer 平台线程未交付，未把 running 状态视作 clean。
- B5 shell 错误 1：多行续行将变量声明拼入 benchmark binary 路径；未执行 benchmark，改用显式 `env` 单次调用。
- B5 shell 错误 2：zsh 将带空格的 artifact 列表作为单个 awk 文件参数；未读取 benchmark，改为显式传入每个文件。
- `/tmp` 下失败尝试和有效 raw artifact 均属于本轮临时证据；最终报告只保留摘要和 hash，清理时按精确路径删除。

## B5 Slow Validation Recheck

- 2026-08-19 当前集成分支重新执行 `ctest --preset cuda-slow --output-on-failure`，结果为 7/7、0 failure；但此前默认异步组合运行曾出现 `test_gpu_vs_cpu_multi` 与 `test_gpu_sim_batch` 的间歇性数值失败，因此单次通过不足以证明稳定性。
- `tests/CMakeLists.txt` 将 slow/integration 设为 `execution.jobs=1`，并为 GPU slow 测试设置 `RESOURCE_LOCK gpu`；当前没有并行 CTest 造成的直接资源竞争证据。
- `GpuEngine::execute_physical_pipeline()` 的 Batched 路径在 `validate_device_result()` 后读取 compact status；Eigen 路径在状态下载前调用 `context_->synchronize()`。设备上下文使用独立 non-blocking stream，cuBLAS/cuSOLVER 绑定同一 stream。
- `SimBatch::run()` 在 `engine_->step()` 后同步读取 engine state，再按 stable simulation ID 记录 host history；这条跨 host/device boundary 是后续复现和诊断的重点，但尚未证明是根因。
- 当前阶段不修改 tolerance、不宣称 slow Gate 稳定通过；需要重复默认顺序并与 `CUDA_LAUNCH_BLOCKING=1` 对照后再修正文档结论。

## 2026-08-19 Async Stream Regression Investigation

- 生产提交 `0adf191` 的 `GpuEngine::execute_physical_pipeline()` 将 mutual、assembly、state、control 和 status kernel 排入 `context_->stream()`；该 context 使用 `cudaStreamNonBlocking`。
- 同一 physical pipeline 的 boundary H2D 与 Eigen fallback derivative upload 仍通过封装后的同步 `cudaMemcpy`，即默认 stream；这没有形成默认 stream 与 execution stream 之间的显式依赖。
- 回归表现与该数据依赖一致：host solver derivative 已计算，但 state kernel 偶发读取旧的 device derivative，导致 GPU 首步 current/velocity 未更新；`CUDA_LAUNCH_BLOCKING=1` 改变时序后可消除现象。
- 当前诊断 worktree `/tmp/coligunCalc-regression` 的全设备同步和 `COILGUN_DEBUG_CUDA` 输出只能用于定位，不能作为修复或最终证据。
- 复现/修复要求：在未修改的 `0adf191` 上先取得 RED；将相关 H2D 改为 `cudaMemcpyAsync` 到 `context_->stream()`，补充最小稳定性回归测试，删除临时同步和诊断输出，再重复运行 focused/integration/sanitizer。
- 独立复现 worktree `/tmp/coligunCalc-repro` 的 CMake 配置因 FetchContent 下载 doctest 长时间无进展而中止；改用已有 CUDA Release 构建进行基线复现，不重复等待同一网络下载。
- 使用已有 CUDA Release binary 对未修复生产代码取得 RED：`ctest --preset cuda-slow -R '^test_gpu_vs_cpu_multi$' --repeat until-fail:5` 首次即失败，GPU 首步 `gpu_I0=0`、`gpu_v=4`，CPU 首步 `cpu_I0=139.499`；共 36 cases 中 35 passed、1 failed，59 assertions failed。
- 同一 binary 使用 `CUDA_LAUNCH_BLOCKING=1` 单独运行通过（1/1）。该对照支持默认 stream H2D 与 non-blocking execution stream 缺少顺序依赖的根因，不支持放宽 tolerance。
- regression worktree 的最小修复已将 physical pipeline、runtime reset、device assembly test、runtime initialization 和 solver pointer-table 的 H2D 排入 execution stream；初始化局部 host vectors 在 stream synchronize 后才离开作用域，step control time 改为持久 workspace scalar。
- 修复后的 `test_gpu_vs_cpu_multi` 连续 5 次通过（5/5，约 23.8 秒），未保留 `cudaDeviceSynchronize()`、环境调试输出或临时 MESSAGE；正式首步电流断言保留。

## CUDA Review Counter（修复后重新计数）

| 轮次 | 审查范围 | 新发现 | 结论 |
|---:|---|---:|---|
| R01 | execution stream 所有权、H2D ordering、初始化/reset/测试 host buffer 生命周期 | 0 | clean |
| R02 | solver enqueue/validation 边界、cuBLAS/cuSOLVER stream 绑定、device pointer table、active-row compaction、Graph capture 兼容性、资源释放顺序 | 0 | clean |
| R03-pre | state/mutual/assembly/control/status kernel 数据依赖、launch 参数和 mask 语义 | 1：`force_reduction_kernel` 的 shared reduction buffer 固定为 512，但 `valid_threads()` 允许 1024；已用 TDD 用例复现并修复为最多 512 | counter reset |
| R03-post | R03 修复后的 state-kernel launch contract、shared-memory 边界和既有 state-kernel 回归 | 0 | clean；修复后 clean streak 1/8 |
| R04-pre | Graph cache capture/replay、异常清理、runtime mask 和 backend fallback 语义 | 1：capture body 异常会遗留 stream capture 状态；已用 CUDA stream 状态 RED 用例复现并修复为结构化失败、EndCapture 和 graph 清理 | counter reset |
| R04-post | R04 Graph exception cleanup、EndCapture error cleanup、fallback lock 和现有 Graph pipeline | 0 | clean；修复后 clean streak 1/8 |
| R05 | device/host workspace ownership、reset/reuse、异常回滚、初始化 buffer 生命周期、析构顺序（当前单 device contract） | 0 | clean；修复后 clean streak 2/8 |
| R06-pre | steady-state launch validation、pointer attributes、device properties、kernel resource/occupancy、Graph/Direct 共用 launch path | 1：每步重复 `cudaPointerGetAttributes()`/`cudaGetDeviceProperties()`，固定 engine-owned device buffers 被重复检查；已增加 source-local trusted launch 入口，公共 checked API 保留 | counter reset |
| R06-post | trusted launch 调用边界、固定 buffer 生命周期、Graph/Direct/fallback 复用、资源报告和整数 launch contract | 0 | clean；修复后 clean streak 1/8 |
| R07 | device control、stage trigger/completion、active-mask termination、compact status、finite/residual barrier | 0 | clean；clean streak 2/8 |
| R08 | assembly row-major/RHS mapping、active-row compaction、cuBLAS pointer table、in-place factorization、solution scatter、residual barrier | 0 | clean；clean streak 3/8；D2D workspace fusion deferred pending independent design/benchmark |
| R09 | trusted launch 的公共/内部边界、整数与尺寸契约、Direct/Graph 共用路径、CUDA kernel 资源与线程配置 | 0 | clean；CUDA Release 定向 build/test 8/8；clean streak 4/8 |
| R10 | device solver active-row compaction、cuBLAS/cuSOLVER stream、pointer table、residual/status barrier、Graph capture solver 分支 | 0 | clean；solver/resident/Graph pipeline 3/3；clean streak 5/8 |
| R11-pre | Graph capture/replay、runtime mask 上传、Graph topology key、异常恢复和 fallback 生命周期 | 1：每步构造并复制不会参与 key/capture 的 runtime mask；已改为直接使用 topology overload，保留同一 stream mask upload | counter reset |
| R11-post | topology-only Graph selection、runtime mask replay、capture/replay/fallback cleanup | 0 | clean；Graph/multi/sim-batch 回归 4/4；修复后 clean streak 1/8 |
| R12 | Full/Standard/Aggressive mutual 数值语义、9 点索引、cutoff/clamp/AGM、shared reduction、寄存器/local memory/spill | 0 | clean；precision/mutual/elliptic/filament/engine physics 5/5；clean streak 2/8 |
| R13-pre | solver 测试 borrowed device view 的 H2D stream、跨 stream happens-before、residual/solution 生命周期 | 1：测试在 default stream 上传后交给 non-blocking execution stream；已改为 context stream，并补充双语 API 约束 | counter reset |
| R13-post | R13 测试 stream 修复、solver enqueue/validation 以及 integration 稳定性 | 0 | clean；CUDA quick 30/30、integration 8/8、solver 连续 8 次；clean streak 1/8 |
| R14 | physical pipeline/reset/assembly/init/solver 的 H2D/D2D/D2H stream 归属、异步拷贝源生命周期、context/BLAS/cuSOLVER 析构 | 0 | clean；针对性 CUDA 回归 5/5；clean streak 2/8 |
| R15 | GPU thermal workspace device state、thermal kernel/observation stream、inactive restore、CPU/GPU thermal 分支、reset/reuse、释放顺序 | 0 | clean；GPU single/multi thermal 对齐 2/2；clean streak 3/8 |
| R16 | Graph cache key/variant 生命周期、capture Begin/Body/End/Instantiate 清理、runtime mask replay、replay failure fallback、reset/shutdown | 0 | clean；Graph/Graph pipeline 2/2；clean streak 4/8 |
| R17 | mutual/separation、assembly、state/control/status kernel 的索引 guard、shape/乘法溢出、线程块/shared memory、checked/trusted launch | 0 | clean；kernel 回归 6/6；paths/assembly memcheck 0 errors；clean streak 5/8 |
| R18-pre | active-row compaction 后 residual kernel/validation 的 inactive 行处理与低 active ratio 工作量 | 1：全部 B 行仍执行 residual `D*D` 且 inactive 非零 RHS 会失败；新增 RED 后传入 mask，inactive residual 置零并跳过计算；counter reset |
| R18-post | mask-aware residual、active compaction/scatter、Graph fixed-shape solver、resident pipeline 和 solver status barrier | 0 | clean；focused 3/3；solver memcheck 0 errors；clean streak 1/8 |
| R19 | device-control 与 SimBatch/GPU wrapper 的 trigger/completion/active 终止、stable ID、host barrier、excitation 推进、complete-stage 更新 | 0 | clean；multi/SimBatch 2/2；已知 host boundary allocation backlog 不计新 finding；clean streak 2/8 |
