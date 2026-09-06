# CUDA 执行性能优化实施台账

## 目标

根据 `docs/2026-07-28-CUDA-Execution-Optimization-Design.md` 与
`docs/2026-07-28-CUDA-Execution-Optimization-Execution-Plan.md`，完成 B0-B5
全部 Task、review、Block Gate、最终效果评估和清理要求。

## 全局约束

- 项目固定使用 C++20。
- 后续 subagent 固定 `gpt-5.6-luna`、High reasoning、关闭 Fast；已完成 Task 的历史模型记录保持不变。
- Implementer 串行派发；GPU 测试与 benchmark 串行。
- 生产代码采用 TDD：先 RED，再 GREEN，再重构。
- 每个 Task 独立 worktree/branch、唯一 final commit、task review clean。
- 端到端代表性 workload median 默认至少改善 5%；关键 workload 不允许超过 2% 未解释回退。
- 不放宽数值 tolerance；fallback 不计入 GPU speedup。

## 阶段

- [x] 准备：隔离目录、ledger、干净基线
- [x] B0：benchmark schema 与基线
- [x] B1：device-resident solver 与 wrapper policy
- [x] B2：active dataflow、inactive rows、workspace 生命周期（开发 Task 完成；最终 benchmark 延后）
- [x] B3：mutual/state/reduction/residual kernels（B3-T1/T2/T3 已集成且 review clean；正式 benchmark/2% 性能门禁延后 B5）
- [x] B4：backend planner、Graph key、ExecutionReport/docs（B4-T1/T2 已集成且 review clean）
- [ ] B5：回归修复、8 轮 CUDA Review、全量验证、效果报告、whole-branch review、清理

## 当前状态

- 当前阶段：B5 最终验收与分步提交进行中
- 下一 Task：完成当前 HEAD 的文档提交、最终 review、精确清理和本地收尾提交
- B3-T2：Task commit `8163838`，集成 commit `7f0f2de`；controller fixed-package review clean、Brooks 100/100；focused 4/4、扩展回归 9/9、cuda-quick 30/30、memcheck 4/4 均通过，benchmark 按用户要求延后。
- B3-T3：Task commit `a8f5b12`，集成 commit `f1c128e`；controller fixed-package review clean、Brooks 100/100；full CUDA build 98/98、solver/resident/multi/precision 回归通过，memcheck 0 errors、racecheck 0 hazards，benchmark 按用户要求延后。
- 集成分支：`perf/cudaExecutionOptimization`
- B4-T1：Task commit `7dc9000`，集成 commit `420064f`；Luna/High/No Fast review approved，0 Critical/Important、1 个已记录 Minor；定向测试 5/5、集成 CUDA quick 30/30、CUDA Release build 101/101；integration 的 `test_gpu_sim_batch` 既有运行时波动延后 B5。
- B4-T2：Task commit `ea92436`，集成 commit `cff073c`；controller fixed-package review clean（0 findings）；CUDA Release build 110/110、test_gpu_graph 10/10、test_gpu_graph_pipeline 11/11、cuda-quick 30/30；integration 中 `test_gpu_vs_cpu_multi`/`test_gpu_sim_batch` 既有跨测试数值波动（隔离单跑 36/36、17/17 通过），延后 B5。
- 当前集成 HEAD：`ccd306a`（已分步集成 engine completion、trusted launch、kernel/solver 修复；CUDA 全量回归和最终 benchmark schema 已通过）
- B4 Gate：开发正确性与文档同步条件通过（Graph key 只含不可变拓扑/策略、mask 为 device input、backend_selection_reason 语义清晰、API.md/API_cn.md/benchmark schema 同步）；正式性能门禁延后 B5。
- B2-T2：已完成，Task commit `e123e9f`，集成 commit `9501faa`，review clean；focused CTest 3/3、
  19 cases/24,775 assertions、compute-sanitizer memcheck 0 errors。
- B2-T3：已完成，Task commit `f50d283`，集成 commit `e4c0f8f`，review clean；focused CTest 4/4、
  integration 8/8、CUDA Release 98/98、memcheck 0 errors。
- B2 Gate：开发正确性和生命周期条件通过；active-ratio/最终性能 benchmark 延后到 B5，未宣称性能门禁通过。
- B3-T1/T2/T3、B4-T1/T2 均已集成；B5 benchmark 已在全部开发 Task 完成后统一执行。

## B5 已有但待复核证据

- Final benchmark：before `9e50677` 与 after `0adf191` 在同一 RTX 5080 Laptop、CUDA 13.3、driver 610.57.04、CUDA Release 环境各运行 3 份 artifact；每份 165 GPU rows，六份均通过 schema checker。
- after benchmark 统计：495 GPU rows，其中实际 GPU 324、fallback 72、setup 99；实际 GPU rows `Finite=yes` 为 324/324，fallback 未计入 speedup。
- 主要结果：当前 HEAD 与基线各 3 份 artifact 均已生成并通过 schema；稳态样本受 GPU 时钟/运行时状态影响明显，报告不宣称统一固定 workload matrix 5% 主线门禁通过。
- 当前 HEAD 全仓测试证据：CPU Release 21/21、CUDA quick 30/30、CUDA parallel 30/30、CUDA slow 7/7、CUDA integration 8/8；solver/resident/Graph focused memcheck 均为 0 errors。
- Whole-branch review：固定范围 `3ddd86e..0adf191`；Luna 代理平台无 verdict，controller fallback review 结论为 0 Critical、0 Important、0 Minor。
- 最终报告：`docs/benchmarks/2026-07-28-CUDA-Execution-Optimization-Effectiveness.md`。
- 清理状态：未完成。当前仍保留 benchmark baseline worktree、多个历史 stale worktree 注册和 `/tmp` artifact，必须在最终验收后按精确路径清理。
- CUDA Review 计数：R03-pre 修复 shared-memory 边界，R04-pre 修复 Graph capture 异常生命周期，R06-pre 移除 steady-state 重复 launch validation，R11-pre 移除 Graph runtime mask 临时复制，R13-pre 修复测试 borrowed-view 跨 stream 准备问题，R18-pre 修复 active-row residual 无效工作/误报；R19-R24 clean。当前最新修复后的 clean streak 为 7/8，仍需连续完成 1 轮 clean。

## 当前回归阻塞

- `GpuEngine::execute_physical_pipeline()` 的 kernel 使用 `context_->stream()`，但边界 H2D 和 Eigen fallback derivative H2D 使用默认 stream 的同步 `cudaMemcpy`；默认 stream 与 non-blocking execution stream 没有明确 happens-before。
- 该问题可导致 GPU 多阶段首步读取旧 derivative，出现 GPU current/velocity 未更新而 CPU 正常更新的错误结果。
- `/tmp/coligunCalc-regression` 当前仅为诊断现场，包含临时 `cudaDeviceSynchronize()`、调试输出和测试 `MESSAGE`，均不能作为最终修复。
- 计划修复为所有 physical-pipeline H2D 使用 `cudaMemcpyAsync(..., context_->stream())`，并将 reset/test assembly 等同类边界一并审计；不保留全设备同步。

## 错误记录

| 错误 | 尝试 | 根因 | 处理 |
|---|---:|---|---|
| CPU build 缺少 CMakeCache/build.ninja | 2 | 并行配置的嵌套 exec 返回后台 session，控制器丢失句柄并过早启动 build | 保留 session，等待配置明确退出；串行重跑成功 |
| `review-package` permission denied | 1 | 技能脚本缺少 executable bit | 使用 `bash <script>` 调用，不修改技能文件 |
| B0-T2 brief 初次写入错误长 SHA | 1 | 手工补全短 SHA | 派发前用 `git rev-parse` 校验并更正为 `9e50677bf04b...` |
| B0-T2 首个 implementer 线程无进展 | 1 | agent 保持 running，但多次状态/中断请求均无响应且无文件或进程活动 | 关闭停滞线程；同模型、同 High/No Fast 配置重新派发，不切模型 |
| B0-T2 第二个 Sol High 线程同样无进展 | 2 | 新线程在 wait 和 interrupt 后仍无响应、无落盘或进程活动 | 依计划“平台阻塞时改为单 agent”；controller 本地完成，不尝试其他模型/Fast |
| B0-T2 `review-package` 内部 helper permission denied | 1 | `review-package` 在未给 OUTFILE 时直接执行无 executable bit 的 `sdd-workspace` | 显式提供唯一 OUTFILE，继续用 `bash review-package ... OUTFILE` |
| B0-T2 Sol High reviewer 无响应 | 1 | reviewer 在 wait/interrupt 后仍 running、无返回 | 关闭线程；依平台阻塞条款由 controller 按 task-reviewer rubric 审查并修复发现 |
| B1-T1 首次组合 patch 验证失败 | 1 | `findings.md` 含空 Update hunk，`apply_patch` 原子拒绝整组变更 | 删除空 hunk 后重新应用；未产生部分代码修改 |
| B1-T1 integration 首次 0/8 Not Run | 1 | 只构建了四个 focused target，integration 可执行文件尚不存在 | 先完成 `cuda-release` 全 target 构建，再重跑 integration preset |
| B1-T1 临时 solver benchmark 首次编译失败 | 1 | 独立 `c++` 命令遗漏 CUDA header include 路径 | 添加 `/opt/cuda/targets/x86_64-linux/include` 后重编译 base/after |
| B1-T1 Sol High reviewer 无响应 | 1 | reviewer 连续两个 wait 窗口及一次 interrupt 后仍为 `running`，无 verdict/文件/变更 | 关闭线程；依平台阻塞条款由 controller 对固定 package 执行 task-reviewer + Brooks rubric |
| B1-T1 集成构建会话恢复时句柄失效 | 1 | 上一会话的终端 session 已关闭，无法读取原始退出码 | 在集成分支重跑同一增量构建，`ninja: no work to do` 且 exit 0；再运行 focused CUDA tests 4/4 通过 |
| B1-T1 `/tmp` 清理命令受限 | 2 | 环境拒绝 `rm -f`，且 `/tmp` 挂载不支持 `gio trash` | 对预先核验的 15 个精确文件路径分别执行 `unlink`，15/15 删除成功 |
| B1-T2 Sol High implementer 无可用 handoff | 1 | 两个 wait 窗口无 worktree/report/process 活动；interrupt 后返回 BLOCKED，但稍后发现一个未提交测试草稿 | 关闭 agent；controller 核验并校准草稿后按 TDD 本地完成，不切模型/Fast |
| B1-T2 诊断 `CAPTURE` 编译失败 | 1 | doctest 2.5.3 的 `CAPTURE` 宏只接受单参数 | 拆分后完成测试工况诊断；最终删除延迟求值 capture，valid RED 在生产修改前取得 |
| B1-T2 临时 benchmark 首次编译缺 Eigen | 1 | 独立命令错误使用不存在的 Release `_deps/eigen-src` | 读取 CMakeCache，改用实际 `/build/ninja-debug/_deps/eigen-src` 后 base/after 均编译成功 |
| B1-T2 reviewer 两个窗口无初始回执 | 1 | reviewer 无文件/进程活动，直到 interrupt 状态请求后才完成初审 | 保留同模型/High/No Fast；按正式 finding 完成 TDD fix 与两轮 re-review，最终 0 findings |
| B2-T1 Sol High implementer 无 handoff | 1 | agent 先落盘 RED 测试和生产候选，但多个 wait 窗口及两次 interrupt 后仍无回执/report/进程活动并保持 running | 关闭线程；保留其测试先行现场，controller 按计划 fallback 核验 RED/GREEN、完成 benchmark/review，不切模型/Fast |
| B2-T1 coarse artifact schema 首次未运行 | 2 | CUDA Release 无 checker 可执行文件，且 CMake 未定义同名 target（先后返回 127 / unknown target） | 用项目既有 C++20 checker 源码直接编译临时 `/tmp` 工具后重跑；不修改 CMake，不视为 artifact 行为失败 |
| B2-T1 fused-kernel 首次 patch 未应用 | 1 | `launch_device_control` 指针校验区实际换行与组合 patch 上下文不一致 | `apply_patch` 原子拒绝、无部分修改；读取精确片段后拆分 patch |
| B2-T1 final helper 首次重链接失败 | 1 | `c++ -x c++` 作用域延伸到后续 archive，且首次遗漏 CUDA include，导致编译器把静态库当源码解析并产生巨量诊断 | 立即终止进程；改为先带 CUDA include 独立 `-c` helper、再无 `-x` 单独 link，编译和链接均成功 |
| B2-T1 review package 首次 revision range 无效 | 1 | 手工把短 SHA `47cee51` 补成了错误长 SHA | 使用 `git rev-parse 47cee51` 取得真实 `47cee513...`，修正 report 后重新生成；后续禁止手工补全 SHA |
| B2-T1 Sol High reviewer 无交付 | 1 | 180 秒 wait 及 interrupt 后 60 秒仍无 verdict，关闭时保持 `running` | 不切模型/Fast；按既定平台 fallback 对固定 package 执行 controller spec/quality + Brooks review |
| B2-T1 临时分支普通删除被 Git 拒绝 | 1 | cherry-pick 后 Task SHA 与集成 SHA 不在 ancestry，`git branch -d` 按保护规则拒绝 | 先用 `git diff --quiet 47cee51 21c5332` 验证树一致，再按用户清理要求精确执行 `git branch -D task/cuda-opt-b2-t1` |
| B2-T2 worktree 初次 build 缺少 `build/cuda-release` 配置 | 1 | 隔离 worktree 尚未生成 CUDA Release CMakeCache/build.ninja | 先串行执行 `cmake --preset cuda-release`，再执行完整 `cmake --build --preset cuda-release -j 2`；构建 150/150 成功 |
| B2-T2 集成首次被根目录用户改动拦截 | 1 | `docs/API_cn.md` 有未提交但与 Task 内容相同的用户改动 | 精确 stash 单文件、cherry-pick `e123e9f`、恢复并核对；未丢弃用户改动 |
