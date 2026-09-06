# coligunCalc 性能优化设计框架

**状态：** 可复用框架，待后续具体优化轮次填充
**用途：** 为后续性能优化设计提供统一结构、决策门槛、benchmark 证据要求和回滚规则。

本文档只保留设计框架，不记录某一轮具体优化方案。后续启动新一轮优化时，应先复制或扩展各章节中的占位内容，再把具体发现、目标数值、实现方案和验证结果写入对应位置。

## 1. 设计目标

性能优化设计必须先明确目标，再讨论实现。每轮设计至少填写以下内容：

- **主要目标：** 需要改善的执行路径、典型 workload、期望收益和优先级。
- **正确性目标：** 必须保持的数值语义、公共 API 语义、默认配置和兼容路径。
- **性能主指标：** 以仿真执行时间、吞吐、延迟、内存访问、传输或同步成本为主；编译时间和测试反馈为次级工程指标。
- **可观测性目标：** 每个优化必须能通过固定 workload、重复样本、阶段计时和数值回归验证。
- **回退目标：** 每个高风险优化都必须保留清晰的禁用、fallback 或回滚路径。

### 1.1 必须达成的原则

1. **正确性优先。** 不以降低精度、跳过必要计算或改变默认物理语义换取 benchmark 数字。
2. **执行速度优先。** 若编译时间增加能够换来经过测量确认的执行速度收益，可以接受更长编译时间。
3. **证据优先。** 不能只凭代码直觉合入优化；必须有优化前后同机同 workload 数据。
4. **边界清晰。** 每项优化必须说明影响的模块、写集、公共接口变化和兼容策略。
5. **可恢复。** 失败的优化应能独立回滚，不拖累已验证的优化成果。

### 1.2 目标指标模板

| 方向 | 当前证据 | 目标 | 验证方式 |
|---|---|---|---|
| <runtime path> | <profile/benchmark evidence> | <target threshold> | <benchmark/test command> |
| <memory/transfer path> | <observed cost> | <target threshold> | <measurement protocol> |
| <batch/throughput path> | <baseline throughput> | <target threshold> | <fixed workload> |
| <engineering feedback> | <build/test evidence> | <secondary target> | <smoke/CI check> |

目标指标是门槛，不是承诺。没有达到门槛时，设计必须给出保留、延期、缩小范围或回滚决定。

## 2. 范围与非目标

### 2.1 本轮范围模板

填写本轮允许修改或评估的区域：

- CPU runtime path：<components/functions/tests>
- CUDA runtime path：<components/functions/tests>
- 数据布局 / 缓存 / workspace：<owned state and lifecycle>
- benchmark / observability：<tools, reports, output schema>
- 构建 / 测试反馈：<presets, labels, CI artifacts>
- 文档同步：<docs that must be updated together>

### 2.2 非目标模板

每轮设计都应明确不做什么：

- 不改变的数值模型或物理假设：<list>
- 不引入的公共 API 破坏：<list>
- 不采用的近似算法或风险实现：<list>
- 不作为主门槛的工程指标：<list>
- 暂不处理但保留记录的候选项：<list>

## 3. 候选项记录模板

每个性能候选项进入设计前，至少补齐以下字段：

| 字段 | 内容 |
|---|---|
| ID | <CPU-xx / CUDA-xx / BUILD-xx / BENCH-xx> |
| 位置 | <file:line or component> |
| 现象 | <what is slow or repeated> |
| 机制 | <why it costs time/memory/sync> |
| 当前证据 | <profile, benchmark, trace, code path> |
| 预期收益 | <local and end-to-end hypothesis> |
| 风险 | <numeric, API, concurrency, memory, compatibility> |
| 设计方案 | <selected approach> |
| 替代方案 | <rejected options and why> |
| benchmark | <before/after command and workload> |
| 回滚方式 | <flag, fallback, revert scope> |

候选项只有在证据、风险和验证方式明确后，才能进入实现计划。

## 4. 方案选择框架

### 4.1 推荐方案结构

推荐方案应按依赖顺序拆分为阶段：

1. **测量基线。** 固定 workload、统一计时边界、记录环境和数值结果。
2. **局部热路径优化。** 先处理 profile 中证据最强、影响最集中的路径。
3. **数据生命周期优化。** 处理重复分配、重复传输、重复同步或重复扫描。
4. **批处理 / 并行 / 后端选择。** 处理规模变化、设备能力和 workload 自适应。
5. **反馈与文档。** 更新 benchmark schema、复现协议和必要的工程反馈入口。

每个阶段都应能独立验证和回滚。

### 4.2 替代方案评估模板

| 方案 | 优点 | 风险 | 是否采用 | 原因 |
|---|---|---|---|---|
| <option A> | <benefits> | <risks> | <yes/no> | <decision> |
| <option B> | <benefits> | <risks> | <yes/no> | <decision> |

不采用的方案也要记录，避免后续重复讨论同一设计分叉。

## 5. 总体架构框架

### 5.1 CPU 数据流占位

~~~text
<immutable inputs / cached metadata>
              |
              v
<hot-path evaluator or workspace>
              |
              v
<assembly / solve / state transition>
              |
              v
<history / summary / outputs>
~~~

填写要求：

- 哪些输入是不变量，哪些状态每步变化。
- 哪些计算可共享，哪些不能跨状态复用。
- workspace 的所有权、生命周期和线程安全边界。
- Reference/diagnostic 路径如何保留。

### 5.2 CUDA 数据流占位

~~~text
<host boundary>
      |
      v
<device-resident state/workspace>
      |
      +--> <kernel/solver/control/thermal path>
      |
      v
<compact outputs / synchronization boundary>
~~~

填写要求：

- 哪些数据必须留在 device，哪些数据必须回 host。
- backend/solver/precision/fallback 的选择规则。
- 同步点、transfer 点和 error/status 上报路径。
- 不同 GPU 可用性下的行为和报告字段。

## 6. 组件设计占位

### 6.1 CPU 组件

| 组件 | 设计要点 | 验证 |
|---|---|---|
| <component> | <cache/workspace/evaluator/lifecycle> | <numeric + perf check> |

### 6.2 CUDA 组件

| 组件 | 设计要点 | 验证 |
|---|---|---|
| <component> | <workspace/transfer/kernel/backend> | <numeric + perf check> |

### 6.3 Benchmark 与工具

| 工具 | 输出字段 | 用途 |
|---|---|---|
| <benchmark tool> | <phase/timing/status fields> | <decision supported> |

### 6.4 文档与反馈

| 文档/入口 | 更新内容 | 同步要求 |
|---|---|---|
| <doc> | <what changes> | <paired docs or smoke checks> |

## 7. Benchmark 与可观测性

### 7.1 统一测量协议

每个 workload 应记录：

- 构造 / 初始化时间；
- cold first step；
- warm-up 后 steady-state wall time；
- median、p95、标准差或样本离散度；
- 局部阶段计时；
- throughput 指标；
- backend、solver、precision、thermal、fallback 状态；
- 数值结果摘要和最大误差；
- commit、compiler、preset、OS、CPU/GPU、driver/runtime。

默认协议可按项目成本调整，但必须记录 warm-up、repeat 次数和 measured step 数。

### 7.2 Workload 模板

| 类别 | Workload | 目的 | 必要输出 |
|---|---|---|---|
| CPU | <single/multi/integration mode> | <hot path or end-to-end> | <timing + numeric> |
| CUDA | <backend/solver/batch shape> | <device path> | <timing + fallback> |
| Batch | <batch size / active ratio> | <throughput> | <steps/s or simulations/s> |
| 初始化 | <construct/reuse scenario> | <setup cost> | <constructor wall + memory> |
| 工程反馈 | <build/test/CI scenario> | <secondary signal> | <wall + status> |

### 7.3 回归规则

- 只比较同一 workload、同一输出语义、同一机器和同一构建配置。
- setup、first-step、warm-up 和 steady-state 不得混合解释。
- fallback latency 不得作为 GPU speedup。
- benchmark 行必须区分 finite、fallback reason 和实际执行后端。
- 数值容差沿用现有测试契约；不得为了性能优化放宽容差。

## 8. 实施阶段与退出条件

| 阶段 | 内容 | 退出条件 |
|---|---|---|
| 0 | 测量基线与协议 | 基线可重复，环境和计时边界明确 |
| 1 | 首批热路径优化 | 局部 benchmark 和端到端 benchmark 均有前后对比 |
| 2 | 数据生命周期优化 | 分配、传输、同步或扫描成本有可测变化 |
| 3 | 批处理/后端/并行策略 | 不同规模 workload 的选择规则明确且可回退 |
| 4 | 文档、反馈和最终 review | 报告、复现协议、未达标项和后续候选均记录完整 |

阶段之间应串行推进。任一阶段没有通过退出条件，不应启动依赖它的实现阶段。

## 9. 测试与验证矩阵

### 9.1 数值测试

- 核心数学/物理函数与旧实现或 reference oracle 对比。
- CPU 默认路径、优化路径、fallback 路径轨迹回归。
- CUDA 各实际可用 backend、solver 和 fallback 路径回归。
- event、termination、history、summary、reset/reuse 行为回归。
- 边界输入、空输入、inactive/active 状态和长运行稳定性回归。

### 9.2 工程验证

- CPU Release CTest 与 CUDA Release CTest 按项目可用 preset 运行。
- git diff --check 通过。
- benchmark 输出 schema 和字段解析通过。
- GPU 可用与不可用环境的 fallback 语义都被验证。
- 高风险 device lifecycle 或并发改动按需运行 sanitizer 或等价检查。

### 9.3 性能门禁

任何代码优化至少需要：

1. 优化前后同机同 workload 的原始数据；
2. 局部指标与端到端指标；
3. 数值回归结果和最大误差；
4. 未命中目标时的保留/回滚决定；
5. 对 unaffected workload 的回归检查。

## 10. 风险与回滚

| 风险 | 控制措施 | 回滚方式 |
|---|---|---|
| 数值语义漂移 | reference oracle、轨迹回归、最大误差记录 | 恢复 reference path 或禁用优化 |
| workspace 生命周期错误 | 明确 ownership、容量断言、sanitizer | 退回每次调用分配或旧生命周期 |
| 后端选择错误 | 显式配置、report 字段、fallback reason | 关闭 auto path 或固定旧 backend |
| batch/index 映射错误 | stable ID 测试、active/inactive 组合 | 禁用压缩或恢复固定 batch |
| 同步/传输边界错误 | trace、finite/status 检查 | 恢复旧同步边界 |
| benchmark 误读 | schema 校验、phase 分离、fallback 排除 | 作废该 benchmark 结论并重测 |

## 11. 完成定义

一轮性能优化设计只有满足以下条件，才能进入执行计划：

- 目标、范围、非目标和候选项列表完整。
- 每个候选项都有证据、风险、benchmark 和回滚方案。
- 阶段拆分满足依赖顺序，阶段内任务写集可并行或已明确串行。
- 测试矩阵覆盖数值、工程和性能门禁。
- 文档同步和最终报告位置明确。
- 未解决问题被列为 open questions 或 deferred items。

## 12. Open Questions 模板

| 问题 | 影响 | 默认决策 | 需要谁确认 |
|---|---|---|---|
| <question> | <scope/risk/output> | <assumption> | <owner> |

后续开始具体优化轮次时，应先填充本节；答案明确的问题可由执行者自行决策，答案会改变交付物或风险边界的问题再询问用户。
