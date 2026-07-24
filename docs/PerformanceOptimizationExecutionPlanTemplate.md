# coligunCalc 性能优化执行计划框架

**状态：** 可复用框架，待后续具体优化轮次填充
**设计依据：** 对应轮次的性能优化设计框架或设计文档
**执行方式：** subagent-driven development，可按用户指定改为单 agent 执行
**适用分支：** 建议在独立 worktree / feature branch 中执行

本文档只保留执行框架，不记录某一轮具体 Block、Task、benchmark 结果或 review 轮次。后续启动新一轮性能优化时，应先复制本框架，再把具体设计目标、任务拆分、写集、验证命令和验收结果填入占位区域。

## 1. 执行目标

把性能优化设计拆成可恢复、可审查、可回滚的实现单元。执行计划必须回答：

- 哪些 Block 必须串行执行；
- 每个 Block 内哪些 Task 可以并行；
- 每个 Task 的写集、测试、benchmark 和提交边界；
- 如何调度 implementer、reviewer 和 fixer；
- 何时允许进入下一 Block；
- 如何评估本轮优化成效并启动后续 review。

### 1.1 Block / Task 关系

~~~text
Block A
  -> Block B
  -> Block C
  -> ...
  -> Validation / Effectiveness Review
  -> Follow-up Exploration
~~~

Block 之间只能串行：前一个 Block 的所有 Task、Task Review、必要 Fix、Block Gate 全部完成后，才能启动下一个 Block。

Block 内的 Task 只有在写集互不重叠、接口依赖已稳定、验证资源不冲突时，才允许并行执行。否则应拆成后续集成 Task。

## 2. 全局硬约束

### 2.1 性能优先级

1. 正确性和数值稳定性。
2. 仿真执行时间、steady-state latency、batch throughput 或本轮定义的主运行时指标。
3. 运行时内存、传输、同步、分配和启动成本。
4. 编译时间、测试反馈和 CI 便利性。

如果编译时间增加能换来经过 benchmark 证实的执行速度收益，保留该实现。任何编译/测试优化都不能以牺牲运行时速度、精度或覆盖率为代价。

### 2.2 Agent 与 API 限制模板

在启动 subagent 前填入本轮约束：

| 项目 | 约束 |
|---|---|
| 模型 | <model and reasoning effort> |
| 禁用模型 / 模式 | <forbidden models, fast modes, service tiers> |
| 最大并发 subagent | <max active agents> |
| API 并发 / RPM | <concurrency and rate limits> |
| 调度间隔 | <spawn/send/wait/close spacing> |
| 失败重试 | <retry policy without fallback drift> |

默认原则：

- 同一时刻保留 API 余量，不把平台上限当作可用并发。
- implementer wave、reviewer wave 和 fixer wave 分开调度。
- 启动失败时只按本轮允许的模型和模式重试；不得静默切换到未授权模型。
- 连续多次同类平台阻塞时暂停并报告，不扩大任务范围。

### 2.3 代码与分支约束

- 每个 Task 使用独立 worktree 和独立分支，或使用等价的隔离机制。
- 每个 Task 最终恰好对应一个实现 commit；fixer 修改同一 Task 时使用 amend，不新增第二个 Task commit。
- Task report 记录 base commit、initial task commit、final task commit、测试命令、benchmark 数据和 review 结论。
- Controller 按确定顺序逐个集成已 review clean 的 Task commit；不 squash、不合并多个 Task、不在 controller 工作树直接改生产代码。
- 若集成冲突暴露出共享写集或接口依赖，应停止当前 Block，新增 integration Task，而不是临时拼接补丁。
- Task 写集必须明确；未经重新划分，不得并行修改同一文件。
- Critical/Important finding 必须修复并 re-review clean；Minor finding 必须记录到 ledger 或 follow-up。
- Block Gate 失败时，不得启动下一个 Block。

## 3. 单 Task 生命周期

每个 Task 执行以下流程。implementer、reviewer 和 fixer 是否由 subagent 执行，按本轮资源约束决定。

1. Controller 记录当前 base commit，读取计划中的 Task 条目。
2. 生成 Task brief，明确目标、写集、禁止写入、测试、benchmark 和完成格式。
3. Implementer 只修改该 Task 写集，运行覆盖测试和必要 benchmark。
4. Implementer 创建唯一 Task commit，并输出 Task report。
5. Controller 生成 base..final_task_commit 的 review package。
6. Reviewer 检查 spec compliance、code quality、数值风险、性能证据和写集边界。
7. 若有 Critical/Important finding，fixer 在同一 Task 分支修复并 amend 原 Task commit。
8. Reviewer clean 后，Controller 记录 final commit、review verdict、测试和 benchmark 结果。
9. 关闭相关 agent / worktree，进入该 Block 的下一个可集成 Task。

### 3.1 Task report 模板

| 字段 | 内容 |
|---|---|
| Task | <Block-Task ID> |
| Base commit | <sha> |
| Initial task commit | <sha> |
| Final task commit | <sha> |
| Modified files | <files> |
| Summary | <implementation summary> |
| Tests | <commands and results> |
| Benchmarks | <commands, workload, before/after> |
| Numeric delta | <max error / tolerance> |
| Review verdict | <spec / quality> |
| Open concerns | <minor or deferred items> |

### 3.2 Commit ledger 模板

| Task | Base | Initial commit | Final task commit | Review | Integrated |
|---|---|---|---|---|---|
| <B?-T?> | <sha> | <sha> | <sha> | <spec/quality> | <feature-sha> |

一个 Task 行只能填一个 final task commit。Block Gate 只接受 ledger 完整且所有行均为 review clean 的状态。

## 4. Block 调度规则

每个 Block 使用如下波次：

~~~text
Block start
  -> 并行启动本 Block 的 implementer wave（不超过并发上限）
  -> 收集结果，关闭 terminal implementer
  -> 启动 reviewer wave
  -> 对失败 Task 启动 fixer wave
  -> 对 fixer 结果重新 review
  -> Controller 执行 Block Gate
  -> Gate 通过后进入下一个 Block
~~~

资源冲突处理：

- GPU benchmark、长测试、文件生成器或独占设备需要资源锁。
- 同一 Block 内若两个 Task 需要同一资源，应在该资源上排队，不必拆成不同 Block。
- 若一个 Task 依赖另一个 Task 的未完成接口，应移动到后续 Block 或 integration Task。

## 5. Block 总览模板

| Block | 目标 | Task 数 | 并行写集 | Gate |
|---|---|---:|---|---|
| <Block A> | <baseline / protocol / setup> | <n> | <artifact or source areas> | <exit criteria> |
| <Block B> | <first runtime area> | <n> | <non-overlapping files> | <tests + benchmark> |
| <Block C> | <integration area> | <n> | <non-overlapping files> | <numeric + perf gate> |
| <Block ?> | <validation / report / review> | <n> | <artifacts> | <final gate> |

## 6. Block 模板

复制本节为每个具体 Block 填写。

### Block <ID>：<Name>

**串行前置：** <previous block / prerequisite>
**并行原则：** <why tasks are independent>
**共享资源：** <GPU / test runner / generated data / none>

#### <ID>-T1 <Task name>

- **目标：** <what this task changes or measures>
- **写集：** <allowed files / directories>
- **禁止写入：** <forbidden files / boundaries>
- **输入：** <design sections, baseline data, APIs>
- **验证：** <tests, benchmark, numeric checks>
- **提交：** <one task commit>
- **退出：** <task-specific done condition>

#### <ID>-T2 <Task name>

- **目标：** <what this task changes or measures>
- **写集：** <allowed files / directories>
- **禁止写入：** <forbidden files / boundaries>
- **输入：** <design sections, baseline data, APIs>
- **验证：** <tests, benchmark, numeric checks>
- **提交：** <one task commit>
- **退出：** <task-specific done condition>

### Block <ID> Gate

- [ ] 所有 Task 有唯一 final task commit。
- [ ] 所有 Task review clean。
- [ ] Block-specific tests 通过或失败原因已记录。
- [ ] benchmark 原始样本和环境信息已保存。
- [ ] 数值差异在容差内，或回滚/延期结论明确。
- [ ] 写集没有越界。
- [ ] runtime performance gate 通过；compile/test feedback 只作次要记录。
- [ ] 所有相关 agent / worktree 已关闭或清理。

## 7. Benchmark 与效果评估计划

### 7.1 Baseline Block 要求

Baseline Block 必须在实现前建立：

- 固定 workload；
- setup / first-step / warm-up / steady-state 边界；
- repeat 数量和样本统计；
- 环境信息；
- 数值输出字段和容差；
- fallback / backend / solver / precision 语义；
- 原始 stdout 或机器可解析 artifact。

### 7.2 每个实现 Task 的性能证据

每个性能 Task 至少记录：

1. 优化前数据；
2. 优化后数据；
3. 局部指标；
4. 端到端指标；
5. 样本离散度；
6. 最大数值误差；
7. 是否命中目标；
8. 未命中时的保留、回滚或延期决定。

### 7.3 Final Effectiveness Block

最终效果评估应作为独立 Block 或独立 Task 完成：

- 汇总所有关键 workload 的 before/after；
- 区分净收益、混合收益和净回退；
- 列出未达标目标和原因；
- 记录环境敏感性；
- 给出下一轮 review 的基线 commit；
- 不以编译时间改善掩盖运行时回退。

## 8. Review 与后续探索计划

### 8.1 Task Review

每个 Task 的 reviewer 必须检查：

- 是否只修改允许写集；
- 是否满足 Task brief；
- 是否保持数值语义；
- 是否有足够 benchmark 证据；
- 是否引入公共 API / 文档同步义务；
- 是否需要 fixer 或 follow-up。

### 8.2 Whole-Branch Review

全部实现 Block 完成后，执行 whole-branch review：

- 设计/计划合规性；
- 生产代码质量和跨模块影响；
- 数值/性能回归风险；
- benchmark 解释是否可靠；
- 文档和报告是否同步；
- Critical/Important/Minor findings；
- 是否允许进入最终效果评估或 finishing workflow。

### 8.3 Follow-up Exploration

效果评估完成后，若用户要求继续探索，则启动只读 review loop：

1. 以优化后代码、效果报告和历史 findings 为基线。
2. 拆分多个独立 review scope。
3. 每轮综合去重，记录新发现数量。
4. 达到用户指定的停止条件后生成最终 review summary。
5. 新发现不在探索 Block 中隐式修复；需要新设计和新执行计划。

## 9. 恢复、失败和回滚

- 会话中断后，先读取本执行计划、设计模板、当前 final reports 和 git log，不重派已 clean Task。
- Implementer NEEDS_CONTEXT：补充 brief/context 后以同一模型和约束重新派发。
- Implementer BLOCKED：先诊断，必要时拆 Task；不得盲目重复相同 prompt。
- Reviewer 发现 Important 以上：当前 Block 停止，fixer 处理完整 findings list 后 re-review。
- Block Gate 性能未达标：保留数据，按设计文档决定回滚、延期或接受未达标。
- 数值或内存安全回归优先回滚该 Task。
- 不使用 destructive git 命令清理用户变更。
- Block 已集成后的回滚或新修复必须创建明确 follow-up Task 和 commit。

## 10. Block Gate 通用验收表

Controller 在每个 Gate 逐项记录：

- [ ] 本 Block 所有 Task implementer 已完成。
- [ ] 每个 Task 都只有一个 final task commit。
- [ ] Ledger 已记录 base / initial / final / integrated SHA。
- [ ] 每个 Task 都有独立 reviewer 的 spec / quality verdict。
- [ ] Critical / Important findings 已修复并 re-review clean。
- [ ] Minor findings 已写入 ledger 或 follow-up。
- [ ] Tests 命令和输出已保存。
- [ ] Benchmark 原始样本、环境和数值误差已保存。
- [ ] 写集没有越界，没有两个并行 Task 修改同一文件。
- [ ] Agent 已关闭，活动数回到调度允许状态。
- [ ] Runtime performance gate 通过，或未达标处理结论明确。

## 11. 完成定义

只有同时满足以下条件，才可以声称本轮执行完成：

1. 所有计划 Block 按顺序完成，所有 Block Gate 有记录。
2. 全部 Task 具有唯一 final task commit、测试报告和 review clean 记录。
3. 全量相关测试通过，或失败范围和延期原因已明确。
4. 关键 workload 的 benchmark 样本可复现。
5. 设计文档中的执行速度目标逐项有达标、延期或回滚结论。
6. Whole-branch review 无未处理 Critical/Important findings。
7. 生产代码、测试、benchmark、文档和规划文件的变更边界已审计。
8. 最终效果评估已完成，并为后续 review 提供基线。

## 12. Open Questions 模板

| 问题 | 影响 | 默认决策 | 需要谁确认 |
|---|---|---|---|
| <question> | <scope/risk/output> | <assumption> | <owner> |

答案明确的问题可由执行者自行决策。答案会改变交付物、风险边界、用户可见输出或外部状态的问题，应在执行前询问用户。
