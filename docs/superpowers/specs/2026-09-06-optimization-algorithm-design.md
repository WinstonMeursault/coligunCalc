# 优化模块设计：开放式约束遗传算法框架

**日期：** 2026-09-06  
**状态：** 设计已确认，待实现

## 1. 目标与范围

本模块用于在多约束条件下优化线圈炮配置，当前首要场景是最大化末速度。模块采用混合变量遗传算法，并将变量定义、物理仿真、目标、约束和遗传搜索解耦，以便后续增加约束类型、目标数量和其他优化算法。

本设计不改变现有线圈炮物理模型、仿真器数值语义或组件 API。优化模块通过问题适配器调用现有仿真接口。

## 2. 设计原则

- **领域无关的搜索核心：** 遗传算法不依赖线圈、温升或触发位置等物理概念。
- **多目标是一等能力：** 单目标是目标向量长度为 1 的特例。
- **约束独立演进：** 约束判断不写入交叉、变异和选择算法。
- **批量优先：** 优先批量评估一代候选，以对接 GPU batch；保留 CPU 串行兼容路径。
- **可复现与可诊断：** 随机种子、评估状态、失败原因和运行统计必须进入结果。
- **物理结果复核：** 优化候选须能使用更高精度 Reference 配置重新仿真。

## 3. 分层架构

```text
OptimizationProblem
  ├─ VariableSchema
  ├─ ObjectiveSpec[]
  ├─ ConstraintSpec[]
  ├─ RepairPolicy
  └─ Evaluator / BatchEvaluator

GeneticOptimizer
  ├─ Population
  ├─ Encoding and repair
  ├─ Selection strategy
  ├─ Crossover
  ├─ Mutation
  └─ Termination policy

OptimizationResult
  ├─ Pareto front
  ├─ Best candidate by objective
  └─ Statistics and diagnostics
```

优化核心只依赖 `OptimizationProblem` 和选择策略接口。`CoilgunOptimizationProblem` 负责将抽象变量解码为 `DrivingCoil`、`Armature`、`TriggerConfig` 和激励参数，调用 `MultiStageSim`，并从 `MultiStageResult` 提取目标与约束指标。

## 4. 变量模型

`VariableSchema` 为每个变量声明稳定的 ID、类型、边界和编码信息。第一版支持：

- 连续变量：实数编码，边界修复；
- 整数变量：整数编码，交叉/变异后进行边界修复；
- 枚举变量：类别编码，使用类别交叉和类别突变。

变量解码和修复由问题定义提供，遗传算法不假设变量对应何种物理量。变量数量和类型在一次优化运行内保持不变。

## 5. 目标与约束模型

每个目标通过 `ObjectiveSpec` 描述：

```text
id / name
direction: Minimize | Maximize
normalization policy
evaluator
```

每个约束通过 `ConstraintSpec` 描述，并生成 `ConstraintReport`：

```text
id / name
kind: Hard | Soft
relation: Equal | LessEqual | GreaterEqual | InRange
value
lower_bound / upper_bound
violation
normalized_violation
```

约束可以是静态变量约束，也可以依赖完整仿真结果的动态约束。新增温升、电流、机械应力、效率或能量损失约束不应修改 GA 核心。

约束比较采用可替换策略，第一版预留：

- `FeasibilityFirst`：硬约束可行性优先，其次比较归一化违反程度；
- `Penalty`：软约束按权重加入目标或评分；
- `Lexicographic`：按约束优先级逐级比较。

硬约束的默认语义是：可行解优于不可行解；两个不可行解比较总归一化违反程度。仿真失败、发散、非有限结果均视为不可行，并保留诊断信息。

## 6. 遗传算法与自动分流

遗传算法采用混合变量编码。连续变量使用 SBX 交叉和多项式变异；整数与枚举变量使用对应的离散交叉/突变及边界修复。种群规模、交叉率、变异率、精英保留、随机种子和终止条件均由配置提供。

选择策略提供：

```text
Auto
SingleObjective
NSGA2
```

`Auto` 在一次运行开始时读取固定的目标数量：

```text
1 个目标  -> SingleObjective
2 个或更多目标 -> NSGA-II
```

单目标路径使用标量比较器；其排序语义与一维 Pareto 支配等价，但避免非支配排序和拥挤距离的额外开销。显式选择 `SingleObjective` 处理多目标，或显式选择 `NSGA2` 处理单目标，应返回配置错误，而不是静默降级。

NSGA-II 使用父代与子代合并、约束优先的非支配排序和拥挤距离保留多样性。目标方向在比较前统一转换，归一化策略由 `ObjectiveSpec` 指定。

## 7. 评估执行模型

优化器以批量接口为主：

```text
BatchEvaluator::evaluate_batch(candidates)
```

输入和输出顺序严格一致。单个候选的失败只影响该候选，不影响同批次其他候选。串行 `Evaluator` 可由适配器包装为批量执行；CPU 实现可选择串行或 OpenMP 并行，GPU 实现直接对接现有 batch API。

评估器应支持稳定的缓存键，避免重复候选重复仿真。缓存命中、总评估数、失败数、GPU fallback 次数和各阶段耗时均进入统计信息。批量评估不得改变随机数种子语义。

## 8. 候选解与结果

```text
Candidate
  variables
  objectives: vector<ObjectiveValue>
  constraints: vector<ConstraintReport>
  evaluation_status: Success | Invalid | Failed
  diagnostics / metadata
```

结果模型为：

```text
OptimizationResult
  pareto_front
  best_by_objective
  statistics
  termination
```

多目标运行默认只返回最终 Pareto front，不自动选择代表解。代表解是显式的后处理操作：

```text
select_representative(RepresentativeSelector)
```

代表解策略作为独立接口，预留 `MaxObjective`、`MinConstraintViolationMargin`、`IdealPointDistance`、`WeightedScore`、`LexicographicObjectives` 和用户自定义回调。代表解选择不修改 Pareto front；空 Pareto front 必须返回明确错误或空结果。

单目标结果也统一填充 `pareto_front`，使调用方无需切换结果类型。

## 9. 线圈炮问题适配器

`CoilgunOptimizationProblem` 负责：

1. 校验并解码候选变量；
2. 构造或复用线圈、`Armature`、激励和触发配置；
3. 调用 CPU 或 GPU 多阶段仿真；
4. 提取末速度作为当前默认目标；
5. 提取温升、峰值电流、峰值电压、效率、能量损失等可选约束指标；
6. 将异常、非有限结果和终止状态转换为统一评估结果。

优化阶段可以使用 `OptimizationLevel::Full` 和 GPU batch；最终候选复核必须允许切换到更高精度、小时间步长和 `Reference` 路径。

## 10. 测试与验证

- 变量编码、边界修复和枚举/整数变异的单元测试；
- 等式、不等式、区间、硬/软约束和归一化违反程度测试；
- 单目标 Auto 分流及其与显式 `SingleObjective` 的等价性测试；
- 多目标 Auto 分流、非支配排序、拥挤距离和 Pareto front 测试；
- 代表解 selector 的确定性和空 front 行为测试；
- 批量/串行评估顺序、单候选失败隔离、缓存命中和随机种子复现测试；
- 使用固定线圈炮场景验证末速度目标、约束提取和 Reference 复核结果；
- 记录种群规模、评估次数、缓存命中、失败、fallback、终止原因和数值摘要。

## 11. 非目标与后续扩展

本阶段不改变现有仿真方程，不强制实现除遗传算法和 NSGA-II 以外的优化器，不自动把多目标压缩为加权单目标，也不在优化器内部决定工程代表解。

后续可在相同 `OptimizationProblem`、`Evaluator`、`ConstraintReport` 和 `OptimizationResult` 接口上增加 DE、CMA-ES、贝叶斯优化、MOEA/D、分布式评估和可视化工具。
