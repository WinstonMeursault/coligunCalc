# 优化模块成果总结与下一阶段开发路线图

**日期：** 2026-09-16  
**基线：** `dev` / `131f549`（`feat(optimization): implement coilgun optimization framework`）  
**规划分支：** `feature/optimization-next-phase-roadmap`  
**状态：** 仅完成复盘与规划，尚未执行本路线图中的实现任务

## 1. 文档目的

本文档完成两件事：

1. 总结当前优化模块已经交付的能力、验证证据与已知边界；
2. 定义下一阶段开发方案，采用已确认的 **A/B 双轨策略**：
   - **A 轨：质量闭环与公共契约稳定化**，作为后续扩展的强制门禁；
   - **B 轨：真实 CUDA 批量优化评估器**，与 A 轨中不冲突的工作并行推进。

本阶段不引入新的优化算法。差分进化、CMA-ES、MOEA/D 等算法待核心契约和真实 GPU 评估链路稳定后再评估。

## 2. 当前成果总结

### 2.1 已完成的领域无关优化核心

当前模块已经形成独立于线圈炮物理模型的优化内核，主要能力如下。

| 能力 | 当前实现 |
|---|---|
| 变量建模 | 连续、整数、枚举三类变量；统一编码、解码与边界修复 |
| 目标建模 | 支持最大化/最小化方向、目标 ID 和多目标值 |
| 约束建模 | 支持硬/软约束、等式/范围等关系及归一化违反度 |
| 单目标搜索 | 遗传算法，含精英保留、锦标赛选择、SBX 与多项式变异 |
| 多目标搜索 | NSGA-II，含非支配排序、拥挤距离和约束支配 |
| 自动分流 | 一个目标自动选择单目标 GA；两个及以上目标自动选择 NSGA-II |
| 显式策略 | 调用方可显式选择策略；策略与目标数不匹配时返回配置错误 |
| 结果模型 | 默认返回 Pareto front；代表解必须由调用方显式选择 |
| 代表解选择 | 支持单目标极值、理想点距离、加权评分、字典序等选择器 |
| 可复现性 | 固定随机种子、终止原因、统计数据和诊断信息进入结果 |

### 2.2 已完成的评估基础设施

评估链路已经具备面向昂贵物理仿真的基础能力：

- 串行与批量评估接口；
- 批量失败后的逐候选隔离重试；
- 无效值、异常和非有限值统一转换为不可行候选；
- 内存缓存以及评估次数、命中次数、失败次数、fallback 次数和耗时统计；
- 批次结果保持输入顺序；
- GPU 回调注入边界和 CPU fallback 行为。

其中，GPU 边界目前只通过测试回调验证，尚无连接现有 CUDA 仿真引擎的生产实现。这是 B 轨的核心缺口。

### 2.3 已完成的线圈炮适配层

`CoilgunOptimizationProblem` 已将优化内核与现有仿真模型隔离，负责：

- 将候选变量映射到线圈、触发、电压等物理参数；
- 保留未参与优化的固定线圈语义；
- 支持完整几何规格映射；
- 计算末速度、峰值电流、峰值电压、最高温度、效率和损耗等指标；
- 将约束报告、失败诊断和仿真 metadata 返回优化器；
- 在 Full 路径完成搜索，并允许使用 Reference 路径复核结果。

适配层明确拒绝当前无法安全映射的 `ArmatureMass`，避免表面可配置、实际无效的参数绑定。

### 2.4 已完成的公共 API 与工程集成

- 优化 API 已由 `<coilgun/coilgun.hpp>` 总头文件导出；
- CMake 安装与包导出已覆盖优化头文件及实现；
- 外部 consumer 可通过 `find_package(coilgun)` 构建并运行；
- `docs/API.md` 与 `docs/API_cn.md` 已加入优化入口说明；
- 优化相关测试已纳入 CTest；
- CPU/CUDA Debug、Release preset 均已有历史验证记录。

### 2.5 已完成的端到端验证

现有基准 `docs/benchmarks/optimization-2026-09-08.md` 使用固定种子 `20260908`，覆盖受约束末速度优化、缓存、失败隔离、fallback 和 Reference 复核。

| 项目 | 记录值 |
|---|---:|
| 优化结果末速度 | `0.009764939958 m/s` |
| Reference 复核值 | `0.009764969509 m/s` |
| Full/Reference 误差 | `2.955109129e-08` |
| 可行性 | 通过 |

这证明 CPU 端完整工作流已经贯通，但该基准没有运行真实 GPU 优化评估器，因此不能用于宣称 GPU 优化吞吐收益。

## 3. 当前遗留问题与优先级

2026-09-11 的独立全分支审查结论为“有条件通过”：无 Critical，但仍有会产生静默错误或削弱扩展性的 Important 问题。对当前 `dev` 的抽查确认，以下关键问题仍存在。

### 3.1 P0：必须先解决的正确性与契约问题

| 编号 | 问题 | 风险 | 处理方向 |
|---|---|---|---|
| A1 | 缓存键不包含问题/评估器身份 | 跨问题共享缓存时可能静默返回其他问题的结果 | 给评估器引入稳定 namespace/version，并纳入缓存键 |
| A2 | `Population::initialize` 按值接收 RNG | 初始种群与后续算子复用随机流片段 | 改为引用传递，并补随机流推进回归测试 |
| A3 | 温度约束在未启用热模型时可静默失效 | 调用方以为约束生效，实际指标恒定 | 构造或配置阶段直接拒绝无效组合 |
| A4 | 多目标早停配置被静默忽略 | 预算和终止语义不可信 | 本阶段先显式拒绝；后续再引入超体积收敛指标 |
| A5 | 三条关键测试存在“名称覆盖、行为未覆盖” | 回归可能在测试全绿时进入主线 | 重写方向、异常隔离、物理参数绑定测试 |
| B0 | GPU 优化只有回调边界，没有生产评估器 | 设计中的批量优先承诺没有真正落地 | 实现连接 CUDA batch API 的正式评估器 |

### 3.2 P1：接口一致性与可维护性问题

| 编号 | 问题 | 处理方向 |
|---|---|---|
| A6 | `Penalty` 实现语义与设计文档不一致 | 明确采用加性罚分或更名现有策略；代码、设计、双语 API 文档保持一致 |
| A7 | 终止原因存在两套枚举，并依赖消息字符串反查 | 收敛为结构化终止原因；保留兼容转换但禁止文本解析 |
| A8 | selector 通过指针算术恢复下标 | 将选择器内部协议改为显式 index |
| A9 | `OptimizationProblem` 未拥有变量/目标/约束规格 | 引入稳定的 `ProblemSpec`/描述接口，分阶段迁移，避免一次破坏全部调用方 |
| A10 | 统计对象可变状态的并发语义不清晰 | 明确每次 run 的所有权与聚合边界，为 GPU/候选并行做准备 |
| A11 | 双语 API 只覆盖入口，未覆盖大部分公共符号 | 按类型、字段、默认值和错误条件补齐，并保持中英文等价 |

### 3.3 P2：性能与数值基线问题

现有性能分析显示，真实 GPU 优化集成不能只追求“能跑”，还必须建立可信的数值基线：

- `dM/dx` 在关键工作区间存在求积收敛风险；
- GPU 的优势在大批量吞吐，不在 batch=1 延迟；
- FP32 aggressive 路径有约 9.5–12 倍潜力，但必须通过长仿真误差门禁；
- 沿轨迹的 `M(x)` 制表与插值导数可能带来数量级 CPU 收益，但属于后续数值内核项目；
- 候选级并行、几何构造缓存与自感表范围扩展仍有较大收益空间。

因此，B 轨第一版以“固定几何、大批量、可证明正确”为目标，不同时引入求积算法重构、FP32 默认化和几何制表。

## 4. 下一阶段目标与非目标

### 4.1 目标

1. 消除可能产生静默错误结果的 P0 问题；
2. 冻结 evaluator、cache key、failure isolation、statistics 四组跨轨接口；
3. 提供真实 CUDA 批量优化评估器，并保留逐候选诊断和 CPU fallback；
4. 形成可重复的 CPU/GPU 吞吐、数值误差和失败统计基准；
5. 补齐双语公共 API 文档；
6. 通过完整 CPU/CUDA Debug/Release 验证与独立全分支审查。

### 4.2 非目标

- 不在本阶段新增 DE、CMA-ES、MOEA/D 或其他搜索算法；
- 不修改现有物理方程和结果语义；
- 不把 FP32 aggressive 设为默认精度；
- 不在第一版 GPU 评估器中支持任意候选几何重建；
- 不以单次仿真延迟作为 GPU 成功标准；
- 不承诺解决性能报告中的全部 CPU/GPU 内核问题。

## 5. A/B 双轨设计

### 5.1 A 轨：质量闭环与接口稳定化

#### A-T1 缓存身份与隔离

采用 evaluator-owned namespace：每个评估器提供稳定的 cache namespace 和 schema/version，缓存键至少包含：

```text
evaluator_namespace + evaluator_version + seed + fallback_mode + encoded_variables
```

验收重点：两个不同 evaluator 使用同一个缓存实例时不得互相命中；同 evaluator、同版本和同输入必须稳定命中。

#### A-T2 随机流与关键回归测试

- `Population::initialize` 消耗调用方持有的 RNG；
- 固定 seed 仍可完整复现；
- 初始化后的下一抽样不得与新 RNG 的第一抽样重复；
- 重写目标方向、未知异常隔离、匝数/触发绑定行为测试，使断言能够真实失败。

#### A-T3 约束与终止语义

- 配置 `MaximumTemperature` 时必须启用热模型，否则返回明确配置错误；
- 对 NSGA-II 的 `max_no_improvement_generations` 先采用 fail-fast；
- 决定并实现 `Penalty` 的唯一正式语义；建议采用与设计一致的加性罚分，同时要求尺度/归一化显式；
- 合并终止原因，消除对英文 message 的逻辑依赖。

#### A-T4 公共契约硬化

- selector 使用显式下标；
- 引入只读 `ProblemSpec`，聚合变量、目标、约束和 repair policy；
- 先提供新接口并保持旧构造方式兼容，再在后续主版本移除重复入口；
- 明确 statistics 是 run-local 快照，不允许多个并发 run 共享可变计数器。

#### A-T5 文档与工程闭环

- 同步补齐 `docs/API.md` 与 `docs/API_cn.md`；
- 修正旧计划中的无效 preset 名称和已知范围缩减；
- 将端到端数值摘要纳入可执行验证，而不是只保留手工 benchmark；
- 更新优化 SDD 台账和 deferred 列表。

### 5.2 B 轨：真实 CUDA 批量优化评估器

#### B-T1 数值基线与可行性探针

在修改生产接口前建立以下基线：

- 固定几何、固定时间步、多个电压/触发候选；
- CPU Reference、CPU Full、CUDA Full 的末速度、峰值电流、峰值电压和温度对比；
- batch size 至少覆盖 `1, 8, 32, 128`；
- 明确 `gpu_executed`、实际 backend 和 fallback 原因；
- GPU 运行必须证明不是 `BackendMode::Fallback`。

该任务可以与 A-T1/A-T2 并行，但只允许写测试探针、基准和设计记录，不得提前固化跨轨接口。

#### B-T2 CUDA 批量评估器生产实现

第一版只支持可共享构造数据的固定几何候选，优先允许变化：

- 激励电压；
- 触发阈值/触发时间；
- 其他不会导致线圈、电枢几何重建的运行时参数。

评估器职责：

1. 将一批 `CandidateVariables` 映射为 `SimBatch` 输入；
2. 一次调用现有 CUDA batch 引擎；
3. 将每个仿真行恢复为同序 `EvaluationResult`；
4. 保留逐候选 Invalid/Failed/diagnostics；
5. 对不支持的候选类型显式拒绝或受控 fallback，不得静默走 CPU；
6. 将实际 backend、GPU 执行标记和耗时写入统计。

#### B-T3 失败隔离、fallback 与统计集成

失败模型按以下层级处理：

| 失败范围 | 行为 |
|---|---|
| 单候选物理无效 | 仅该候选 Invalid/Failed，保留同批成功结果 |
| 单候选设备计算失败 | 仅对该候选进行允许的 CPU fallback，并递增对应计数 |
| 整批 CUDA 调用失败 | 按配置整体 CPU fallback 或终止；原因必须结构化记录 |
| 结果数量/顺序不一致 | 视为评估器协议错误，不接受部分静默修复 |

统计至少区分：请求候选数、GPU 实际执行数、GPU 成功/失败数、CPU fallback 数、缓存命中数、批次数、传输时间、内核时间和端到端时间。

#### B-T4 CPU/GPU 优化工作流基准

建立新的版本化 benchmark 文档，要求：

- 同一 seed、同一初始种群、同一物理配置；
- CPU 与 GPU 返回相同的候选顺序语义；
- Pareto front 或单目标最优解通过 CPU Reference 复核；
- 分别报告 batch=1 延迟和大批量吞吐；
- 不用理论估算替代实测；
- 记录 source revision、工作区状态、编译器、GPU、驱动、precision 和 backend；
- 所有性能结论同时给出数值误差门禁。

### 5.3 跨轨接口冻结门禁

B-T2 合入前，A 轨必须冻结下列契约：

| 契约 | 冻结条件 |
|---|---|
| Evaluator identity | 缓存 namespace/version 已定义并有跨 evaluator 测试 |
| Cache key | 编码稳定、覆盖上下文与 evaluator 身份、有版本升级策略 |
| Failure isolation | Invalid、Failed、batch protocol error、fallback 的结构化语义固定 |
| Statistics | run-local 所有权、CPU/GPU/fallback 计数口径固定 |

未通过门禁时，B 轨只能保留为实验分支，不得接入公共优化 API。

## 6. 依赖关系与并发开发波次

后续执行继续采用 TDD 与 subagent-driven development。每个实现任务遵循 Red → Green → Refactor，并由未参与实现的 reviewer 独立审查。

| 波次 | A 轨 | B 轨 | 合流条件 |
|---|---|---|---|
| W0 基线 | 建立问题复现测试、确认当前 CPU/CUDA 基线 | 建立真实 CUDA 可行性探针 | 所有失败均可重复，基线结果归档 |
| W1 可并行 | A-T1 缓存；A-T2 RNG/关键测试 | B-T1 数值基线与探针 | 不修改同一公共接口文件；测试数据可共享 |
| W2 A 轨门禁 | A-T3 约束/终止；A-T4 契约硬化 | 根据冻结草案调整实验适配器 | evaluator/cache/failure/statistics 审查通过 |
| W3 可并行 | A-T5 双语文档与工程闭环 | B-T2 生产 CUDA 批量评估器 | 公共接口不再漂移；各自测试通过 |
| W4 集成 | 审查 API 与统计口径 | B-T3 失败隔离；B-T4 基准 | CPU/GPU 端到端结果及 Reference 复核通过 |
| W5 收口 | 全需求审计 | 全 preset 验证与性能证据归档 | 无 Critical/Important 未处理项 |

并发约束：

- 每个 subagent 使用独立 worktree 和独立任务分支；
- 同一时间最多一个任务修改公共优化类型或 CMake 导出面；
- reviewer 不复用 implementer 的上下文结论；
- 发现共享接口变化时，先回到 A 轨契约任务，不允许在 B 轨临时扩展；
- 文档任务可以与实现并行，但双语 API 文件由单一 owner 串行修改。

## 7. 任务与 commit 拆分

每个 commit 只承担一个可独立审查、可独立回滚的逻辑变化。commit subject 使用简洁英文、首字母大写、祈使语气。

| 顺序 | 任务 | 计划 commit message | 主要产物 |
|---:|---|---|---|
| 1 | A-T1 | `Fix optimization cache identity` | evaluator namespace/version、缓存隔离测试 |
| 2 | A-T2a | `Advance optimization random stream` | RNG 引用语义和回归测试 |
| 3 | A-T2b | `Strengthen optimization regression tests` | 方向、异常隔离、物理绑定测试 |
| 4 | A-T3a | `Validate thermal optimization constraints` | 热模型配置校验 |
| 5 | A-T3b | `Align optimization constraint semantics` | Penalty 语义及测试 |
| 6 | A-T3c | `Unify optimization termination reasons` | 结构化终止原因、NSGA-II 配置校验 |
| 7 | A-T4a | `Harden optimization result selectors` | index-based selector |
| 8 | A-T4b | `Stabilize optimization problem contract` | `ProblemSpec` 与兼容迁移层 |
| 9 | A-T4c | `Define optimization statistics ownership` | run-local 统计与并发约束 |
| 10 | B-T1 | `Establish optimization numerical baseline` | CUDA/CPU 探针和数值证据 |
| 11 | B-T2 | `Add CUDA batch optimization evaluator` | 生产批量评估器 |
| 12 | B-T3 | `Integrate GPU optimization metrics` | 失败隔离、fallback、统计 |
| 13 | A-T5 | `Document optimization public API` | 同步更新 API.md/API_cn.md 与 SDD 台账 |
| 14 | B-T4 | `Benchmark GPU optimization workflow` | 版本化实测报告 |
| 15 | 收口 | `Validate next optimization phase` | 端到端验证、报告与门禁证据 |

如果某任务在 Red 阶段暴露更大的接口问题，应先修改计划并重新评审，不通过追加“临时修复 commit”掩盖任务边界。任务完成前将修复 amend/squash 到对应逻辑 commit；跨任务发现的问题另开明确任务。

## 8. TDD 与审查协议

每个任务必须留下三类证据：

1. **Red：** 新测试在旧实现上按预期失败，记录失败原因；
2. **Green：** 最小实现使目标测试通过；
3. **Refactor：** 清理实现后，相关测试和回归测试仍通过。

审查采用两级门禁：

- 任务级审查：核对需求、测试有效性、公共 API 影响和回滚能力；
- 波次级审查：核对 A/B 两轨契约是否一致，避免实验代码泄漏到公共层。

Critical 和 Important finding 必须在进入下一波次前清零，或由文档明确降级范围并获得新的设计确认。Minor finding 可以登记 deferred，但不得影响正确性、可复现性和诊断能力。

## 9. 验收标准

### 9.1 功能与正确性

- 单目标仍自动路由至 GA，多目标仍自动路由至 NSGA-II；
- 默认结果仍只返回 Pareto front，代表解仍需显式请求；
- 现有 CPU 行为和固定 seed 可复现性不退化；
- 不同 evaluator 共享缓存不会交叉污染；
- 无效热约束和无效多目标早停配置在运行前失败；
- GPU 批次内单候选失败不污染成功候选；
- GPU 结果顺序与输入候选顺序完全一致；
- 至少一个 GPU 优化结果通过 CPU Reference 复核；
- 请求 GPU 的验证必须断言 `gpu_executed == true` 且 backend 不是 Fallback。

### 9.2 工程验证

至少执行：

```bash
cmake --preset cpu-debug
cmake --build --preset cpu-debug
ctest --preset cpu-debug --output-on-failure

cmake --preset cpu-release
cmake --build --preset cpu-release
ctest --preset cpu-release --output-on-failure

cmake --preset cuda-debug
cmake --build --preset cuda-debug
ctest --preset cuda-debug --output-on-failure

cmake --preset cuda-release
cmake --build --preset cuda-release
ctest --preset cuda-release --output-on-failure
```

此外必须：

- 独立运行真实 GPU 优化集成测试；
- 记录 GPU 设备、驱动、CUDA Toolkit、precision 与 backend；
- 构建一个只依赖安装包和总头文件的外部 consumer；
- 对双语 API 文档执行结构与内容等价检查；
- 对 benchmark 输出执行 schema 校验；
- 最终运行全分支独立审查。

### 9.3 性能成功标准

第一阶段不预设夸张的加速倍数，采用可复现门槛：

- batch=32 和 batch=128 的 GPU 每候选吞吐必须优于同轮 CPU；
- 所有加速比使用同一轮、同一物理配置和同一精度口径；
- batch=1 单独报告，不与吞吐结论混合；
- 性能提升不能通过放宽既有数值容差获得；
- 若 FP32 被纳入实验，必须单列误差与适用范围，默认仍保留可靠精度路径。

## 10. 分支与 worktree 策略

本规划已在 `feature/optimization-next-phase-roadmap` 上建立，起点为 `dev` 的 `131f549`。

执行阶段建议使用：

```text
feature/optimization-quality-closure   # A 轨集成分支
feature/cuda-batch-optimization        # B 轨集成分支
task/optimization-<task-id>            # 独立任务分支/worktree
```

分支规则：

- 任务分支只从对应轨道最新已审查基线创建；
- A 轨接口冻结后，以明确 commit 为 B 轨 rebase/merge 基点；
- 禁止多个任务共享一个 worktree；
- 合并任务前要求工作区干净、报告完整、相关测试通过；
- 任务分支合入后及时删除，保留 A/B 集成分支直到最终收口；
- 不删除含未提交改动的 worktree，不用强制清理掩盖未知用户改动。

本次规划前已清理旧 OPT-T* 干净中间 worktree/分支；`fix/opt-t6-convergence` 因仍含未提交测试改动而保留，待其内容由所有者确认后再处理。

## 11. 风险、回滚与停止条件

| 风险 | 缓解与回滚 |
|---|---|
| `ProblemSpec` 迁移破坏现有调用方 | 新旧入口并存一个阶段；外部 consumer 作为门禁 |
| 缓存 key 变更导致旧缓存不可读 | 内存缓存直接版本隔离；未来持久缓存必须带 schema version |
| CUDA 批量接口无法表达逐候选状态 | 不压缩错误语义；必要时先扩展内部结果结构，再接公共 API |
| GPU 结果与 CPU 漂移 | 自动 Reference 复核；超阈值立即停止性能优化 |
| 固定几何限制过窄 | 第一版显式声明支持矩阵；不以静默 CPU fallback 伪装支持 |
| 并发开发造成接口漂移 | A 轨拥有公共契约；B 轨在冻结点前只做实验适配 |
| benchmark 受硬件噪声影响 | 多次运行、报告分布、同轮 CPU/GPU 对照，不在单元测试中写墙钟硬阈值 |

出现以下任一情况时停止 B 轨集成并回到设计评审：

- 无法可靠识别 GPU 是否实际执行；
- 无法保留逐候选失败隔离；
- CPU Reference 误差超过既有门槛；
- 需要修改现有物理方程或公共仿真语义；
- 需要让缓存、统计或 fallback 再次依赖消息文本或隐式全局状态。

## 12. 后续阶段候选项

完成本路线图后，再按证据选择下一阶段：

1. **算法扩展：** 为统一 optimizer 接口加入 DE/CMA-ES，或在多目标场景评估 MOEA/D；
2. **数值内核：** 修复 `dM/dx` 收敛基线，实施非对称求积阶数；
3. **轨迹制表：** 正式实现 `OptimizationLevel::LookupTable`；
4. **几何缓存：** 复用线圈/电枢构造和自感计算；
5. **精度分层：** 将 FP32 作为经过验证的可选吞吐模式；
6. **用户体验：** Pareto front 可视化、约束敏感性与运行恢复/检查点。

优先级原则是：先消除静默错误，再稳定契约；先证明真实批量收益，再扩展算法数量。
