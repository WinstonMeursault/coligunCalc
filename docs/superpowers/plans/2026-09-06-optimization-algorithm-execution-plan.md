# 优化模块执行计划

**设计依据：** [开放式约束遗传算法框架](/mnt/data/Project/coligunCalc/docs/superpowers/specs/2026-09-06-optimization-algorithm-design.md)  
**目标：** 实现支持混合变量、多约束、自动单/多目标分流、批量评估和显式代表解选择的遗传优化模块。  
**适用分支：** `feature/optimizationAlgorithm`  
**提交规则：** 每个 Task 一个独立 commit；commit message 使用简洁英文；Task 之间按依赖顺序集成，不 squash。

## 1. 全局约束

- 不修改现有物理方程、仿真结果语义和现有组件 API。
- 优化核心不得依赖 `DrivingCoil`、`Armature` 或具体约束名称。
- 目标数量在一次运行开始时固定；`Auto` 根据目标数量选择策略。
- 约束比较先于目标比较；失败、发散和非有限仿真结果视为不可行。
- 批量评估必须保持输入/输出顺序；单个失败不得影响同批候选。
- 所有随机过程接受显式 seed，并保证相同配置下可复现。
- 新增公共 API 时同步更新 `docs/API.md` 和 `docs/API_cn.md`。
- 不覆盖或恢复工作区中已有的用户删除和未提交变更。

## 2. Commit 总览

| 顺序 | Task | Commit message | 主要写集 | 依赖 |
|---:|---|---|---|---|
| 1 | OPT-T1 | `Add optimization domain types` | `include/coilgun/optimization/`, 基础源文件、CMake | 无 |
| 2 | OPT-T2 | `Implement variable encoding and repair` | 变量 schema、编码、修复实现和测试 | OPT-T1 |
| 3 | OPT-T3 | `Implement objectives and constraints` | 目标、约束、可行性比较和测试 | OPT-T1 |
| 4 | OPT-T4 | `Add batch evaluation interfaces` | evaluator、缓存、统计和测试 | OPT-T2, OPT-T3 |
| 5 | OPT-T5 | `Implement genetic operators` | population、selection、crossover、mutation 和测试 | OPT-T2, OPT-T3 |
| 6 | OPT-T6 | `Add single objective optimization` | GA 主循环、终止策略和测试 | OPT-T4, OPT-T5 |
| 7 | OPT-T7 | `Add NSGA-II selection` | 非支配排序、拥挤距离和测试 | OPT-T4, OPT-T5 |
| 8 | OPT-T8 | `Add automatic strategy routing` | Auto 分流、配置校验和集成测试 | OPT-T6, OPT-T7 |
| 9 | OPT-T9 | `Add coilgun optimization adapter` | 线圈炮问题适配器、仿真指标提取和测试 | OPT-T4, OPT-T8 |
| 10 | OPT-T10 | `Add Pareto result selectors` | 结果模型、代表解 selector 和测试 | OPT-T7, OPT-T9 |
| 11 | OPT-T11 | `Integrate optimization public API` | 总头文件、安装导出、双语 API 文档 | OPT-T10 |
| 12 | OPT-T12 | `Validate optimization workflow` | 端到端测试、benchmark、验证报告 | OPT-T11 |

## 3. Block A：基础领域模型

### OPT-T1：Add optimization domain types

- **目标：** 建立 `CandidateVariables`、`Candidate`、`EvaluationResult`、`OptimizationProblem`、`OptimizationConfig` 和状态/诊断类型的最小公共模型。
- **写集：** `include/coilgun/optimization/types.hpp`、`problem.hpp`、`config.hpp`；对应 `src/optimization/` 源文件；`src/CMakeLists.txt`。
- **禁止写入：** `src/simulation/`、`src/components/`、CUDA 实现。
- **实现要点：** 使用稳定的 ID、值向量和状态枚举；结果对象可携带目标、约束、诊断和统计；定义异常与无效结果的统一表达。
- **验证：** 类型构造、默认值、移动/复制语义和基本状态转换单元测试；CPU debug 构建。
- **退出条件：** 公共模型可被后续 Task 引用，测试通过，commit 唯一。
- **Commit：** `Add optimization domain types`

### OPT-T2：Implement variable encoding and repair

- **目标：** 实现连续、整数、枚举变量 schema，候选编码/解码和边界修复策略。
- **写集：** `include/coilgun/optimization/variables.hpp`、`src/optimization/variables.cpp`、`tests/test_optimization_variables.cpp`、测试 CMake 注册。
- **依赖：** OPT-T1。
- **实现要点：** 连续值按边界裁剪；整数采用确定性舍入和边界修复；枚举拒绝非法索引；schema 在运行期间不可变；非法 schema 在构造时报告。
- **验证：** 空 schema、单变量、多类型混合、上下界、NaN/Inf、非法枚举、重复 ID 和可复现修复测试。
- **退出条件：** 所有变量类型的编码和修复满足确定性及边界不变量。
- **Commit：** `Implement variable encoding and repair`

### OPT-T3：Implement objectives and constraints

- **目标：** 实现目标方向/归一化、硬软约束、等式/不等式/区间关系及可行性比较策略。
- **写集：** `include/coilgun/optimization/objective.hpp`、`constraint.hpp`、`comparator.hpp`；对应源文件和 `tests/test_optimization_constraints.cpp`。
- **依赖：** OPT-T1。
- **实现要点：** 明确 violation 和 normalized violation 公式；有限 scale 必须校验；目标方向转换为统一内部语义；`FeasibilityFirst`、`Penalty`、`Lexicographic` 策略通过接口隔离。
- **验证：** 可行/不可行两两比较、多个硬约束聚合、软约束不改变硬约束优先级、目标方向和不同量纲归一化测试。
- **退出条件：** 比较器不依赖物理类型，边界行为和错误信息明确。
- **Commit：** `Implement objectives and constraints`

### Block A Gate

- [ ] OPT-T1 至 OPT-T3 均已集成且各自只有一个 commit。
- [ ] 基础类型、变量边界和约束比较测试通过。
- [ ] 没有修改现有物理/仿真实现。

## 4. Block B：评估与遗传算子

### OPT-T4：Add batch evaluation interfaces

- **目标：** 提供串行 `Evaluator`、批量 `BatchEvaluator`、串行转批量适配器、缓存接口和评估统计。
- **写集：** `include/coilgun/optimization/evaluator.hpp`、`cache.hpp`、`statistics.hpp`；对应源文件和 `tests/test_optimization_evaluator.cpp`。
- **依赖：** OPT-T2、OPT-T3。
- **实现要点：** 输入/输出顺序严格一致；单候选失败隔离；缓存键由 evaluator 或问题定义提供；统计评估数、失败数、缓存命中、fallback 和耗时；seed 通过评估上下文传递。
- **验证：** 串行包装批量、顺序保持、部分失败、缓存命中、异常转换、空批次和重复 seed 测试。
- **退出条件：** 后续 GA 可以只依赖批量 evaluator，不感知 CPU/GPU 实现。
- **Commit：** `Add batch evaluation interfaces`

### OPT-T5：Implement genetic operators

- **目标：** 实现种群容器、初始化、选择基础设施、SBX/离散交叉、多项式/离散突变、精英保留和随机数上下文。
- **写集：** `include/coilgun/optimization/genetic_operators.hpp`、`population.hpp`；对应源文件和 `tests/test_optimization_operators.cpp`。
- **依赖：** OPT-T2、OPT-T3。
- **实现要点：** 操作符处理混合变量；操作后统一调用 repair；支持配置化概率和分布指数；不把目标/约束逻辑写入交叉和突变。
- **验证：** seed 可复现、概率边界、种群规模保持、变量类型保持、精英不被破坏、修复后候选始终合法。
- **退出条件：** 操作符可独立测试并可被单目标和 NSGA-II 共用。
- **Commit：** `Implement genetic operators`

### Block B Gate

- [ ] 批量 evaluator 和遗传算子测试通过。
- [ ] 相同 seed 产生相同初始种群及操作结果。
- [ ] 失败候选和缓存行为未污染同批次其他结果。

## 5. Block C：选择器与自动分流

### OPT-T6：Add single objective optimization

- **目标：** 实现单目标 GA 主循环、可行性优先选择、精英保留和终止策略。
- **写集：** `include/coilgun/optimization/genetic_optimizer.hpp`、`termination.hpp`；对应源文件和 `tests/test_optimization_single.cpp`。
- **依赖：** OPT-T4、OPT-T5。
- **实现要点：** 目标向量长度必须为 1；支持最大化/最小化；批量评估每一代；终止条件包括代数、评估预算、目标收敛和无改进代数。
- **验证：** 人工构造函数上的收敛、可行性优先、评估预算、早停、失败全体候选和 seed 复现测试。
- **退出条件：** 单目标结果稳定填充 `pareto_front` 和 `best_by_objective`。
- **Commit：** `Add single objective optimization`

### OPT-T7：Add NSGA-II selection

- **目标：** 实现约束优先的非支配排序、拥挤距离和父子代合并选择。
- **写集：** `include/coilgun/optimization/nsga2.hpp`、对应 `src/optimization/nsga2.cpp`、`tests/test_optimization_nsga2.cpp`。
- **依赖：** OPT-T4、OPT-T5。
- **实现要点：** 支持任意固定目标数量 >= 2；目标方向统一；边界点拥挤距离为无穷；重复目标值稳定处理；不可行候选按约束策略排序。
- **验证：** 已知二维 Pareto front、拥挤距离边界、重复点、约束优先、目标方向和种群截断测试。
- **退出条件：** NSGA-II 选择器不依赖线圈炮问题，且可独立替换。
- **Commit：** `Add NSGA-II selection`

### OPT-T8：Add automatic strategy routing

- **目标：** 实现 `SelectionStrategy::Auto` 和显式策略校验。
- **写集：** GA 配置/路由实现、`tests/test_optimization_routing.cpp`。
- **依赖：** OPT-T6、OPT-T7。
- **实现要点：** 运行开始时冻结目标数量；1 个目标路由到单目标选择器，>=2 个目标路由到 NSGA-II；显式策略与目标数量不匹配时报错；不得中途改变策略。
- **验证：** Auto/显式单目标等价、Auto 多目标进入 NSGA-II、非法组合报错、目标数量冻结测试。
- **退出条件：** 自动路由行为完全由配置和目标数量决定。
- **Commit：** `Add automatic strategy routing`

### Block C Gate

- [ ] 单目标和 NSGA-II 单元测试通过。
- [ ] Auto 路由和非法配置测试通过。
- [ ] 同一问题可在不改问题定义的情况下切换显式选择策略。

## 6. Block D：线圈炮适配与结果决策

### OPT-T9：Add coilgun optimization adapter

- **目标：** 将抽象候选解映射到现有线圈炮组件和多阶段仿真，提取末速度及可选约束指标。
- **写集：** `include/coilgun/optimization/coilgun_problem.hpp`、对应源文件、`tests/test_coilgun_optimization.cpp`。
- **依赖：** OPT-T4、OPT-T8。
- **实现要点：** 支持触发位置、几何和电路参数的可配置映射；复用或重建仿真对象；处理 `MultiStageResult` 终止状态；将温升、峰值电流、峰值电压、效率和能量损失暴露为指标；提供 GPU batch 和 CPU fallback 适配点。
- **验证：** 固定场景解码、末速度提取、静态/动态约束、仿真异常、非有限结果、GPU 不可用时 CPU fallback 测试。
- **退出条件：** 适配器之外不出现物理字段耦合，现有仿真测试不回归。
- **Commit：** `Add coilgun optimization adapter`

### OPT-T10：Add Pareto result selectors

- **目标：** 实现结果对象和显式代表解选择器。
- **写集：** `include/coilgun/optimization/result.hpp`、`selectors.hpp`；对应源文件和 `tests/test_optimization_selectors.cpp`。
- **依赖：** OPT-T7、OPT-T9。
- **实现要点：** 默认只返回 Pareto front；提供 `MaxObjective`、`MinConstraintViolationMargin`、`IdealPointDistance`、`WeightedScore`、`LexicographicObjectives` 和自定义 selector；选择不修改原结果；空 front 返回明确错误或空结果。
- **验证：** selector 确定性、方向处理、归一化、空 front、多目标折中和单目标统一结果模型测试。
- **退出条件：** API 不自动生成代表解，调用方可明确选择策略。
- **Commit：** `Add Pareto result selectors`

### Block D Gate

- [ ] 线圈炮适配器能够完成最小端到端仿真评估。
- [ ] Pareto front、约束报告和代表解 selector 行为通过测试。
- [ ] GPU fallback 和高精度 Reference 复核入口存在。

## 7. Block E：公共 API 与验证

### OPT-T11：Integrate optimization public API

- **目标：** 将优化模块接入总头文件、构建安装和双语 API 文档。
- **写集：** `include/coilgun/coilgun.hpp`、必要的 install/export CMake、`docs/API.md`、`docs/API_cn.md`、公共 API smoke test。
- **依赖：** OPT-T10。
- **实现要点：** 保持现有 include 兼容；所有新增公共符号在中英文文档中等价描述；不公开 detail 实现。
- **验证：** clean configure/build、公共头文件编译 smoke test、API 文档同步检查。
- **退出条件：** 使用者只包含总头文件即可构造问题、运行优化并读取 Pareto front。
- **Commit：** `Integrate optimization public API`

### OPT-T12：Validate optimization workflow

- **目标：** 完成端到端验证、性能基线和效果报告。
- **写集：** `tests/test_optimization_integration.cpp`、`tests/bench_optimization.cpp`、`docs/benchmarks/optimization-*.md`。
- **依赖：** OPT-T11。
- **验证命令：** `cmake --preset ninja-debug`、`cmake --build --preset ninja-debug`、`ctest --preset debug`；另行运行优化 benchmark。
- **验证内容：** 固定线圈炮 workload；记录 setup、first-step、warm-up、steady-state；比较 CPU/GPU batch；记录 seed、评估次数、缓存命中、失败、fallback、末速度和 Reference 复核误差。
- **退出条件：** 全量相关测试通过；结果可复现；数值差异在既有容差内；性能回退有明确接受、延期或回滚决定。
- **Commit：** `Validate optimization workflow`

### Block E Gate

- [ ] 构建、lint/typecheck（若项目提供）和相关测试通过。
- [ ] 双语 API 文档同步。
- [ ] benchmark 原始数据和环境信息已保存。
- [ ] 至少一个可行优化结果通过 Reference 复核。
- [ ] 没有未经记录的 Critical/Important 风险。

## 8. 集成顺序与回滚

Controller 按 OPT-T1 到 OPT-T12 顺序集成，每个 Task 只接受一个最终 commit。任一 Task 出现 Critical/Important review finding，先在该 Task 分支修复并 amend 原 commit，再重新 review；不得追加临时修复 commit。若跨 Task 冲突暴露出接口设计问题，暂停当前 Block，新增 integration Task 和独立英文 commit message。

任何数值语义、内存安全或现有仿真回归都优先回滚对应 Task；性能未达目标但正确性通过时，保留 benchmark 数据并记录延期或缩小范围决定，不以放宽数值容差解决。

## 9. 后续扩展入口

完成 OPT-T12 后，可在不修改 `OptimizationProblem`、`Evaluator`、`ConstraintReport` 和 `OptimizationResult` 的前提下增加 DE、CMA-ES、贝叶斯优化、MOEA/D、分布式评估和 Pareto 可视化。任何新算法都应新增独立选择/搜索策略和独立测试，不修改已验证的 GA/NSGA-II 默认路径。

## 10. Subagent-driven Development 调度

本计划允许使用 subagent-driven development，但本节描述的是未来执行协议，编写计划时不启动任何 agent。每个 Task 使用一个全新的 implementer；实现完成后使用独立 reviewer，Critical/Important finding 使用 fixer 修复并重新 review。Controller 负责记录 base commit、Task commit、测试结果和 review verdict。

### 10.1 资源限制

- 最大并发 agent 数：4（包含 implementer、reviewer、fixer，不把 Controller 计入）。
- 同一 Task 同时最多一个 implementer 和一个 reviewer；review 必须等待实现 commit。
- reviewer/fixer wave 与 implementer wave 分开，避免同一写集并发修改。
- GPU batch、CUDA 测试和长时间 benchmark 使用 `gpu` 资源锁，同一时刻最多一个 GPU workload。
- 同一文件或同一公共接口写集不得并行；若写集发生重叠，改为串行 Task 或新增 integration Task。
- implementer 使用能完成任务的最低模型档位；跨模块集成和最终 whole-branch review 使用标准或高能力模型。

### 10.2 并发波次

| 波次 | 前置条件 | 并行任务 | 说明 |
|---|---|---|---|
| W0 | 计划执行前 | Preflight | 检查 branch、`.git/index` 可写、构建 preset、测试依赖和 CUDA 状态；不改生产代码 |
| W1 | W0 通过 | OPT-T1 | 基础类型先串行完成，后续任务依赖其接口 |
| W2 | OPT-T1 review clean | OPT-T2 + OPT-T3 | 变量模型与目标/约束模型写集独立，可并行；各自先 review 再进入 W3 |
| W3 | OPT-T2、OPT-T3 clean | OPT-T4 + OPT-T5 | evaluator 与遗传算子写集独立，可并行；GPU 相关测试按资源锁排队 |
| W4 | OPT-T4、OPT-T5 clean | OPT-T6 + OPT-T7 | 单目标主循环与 NSGA-II 选择器独立，可并行 |
| W5 | OPT-T6、OPT-T7 clean | OPT-T8 | 自动路由依赖两个选择器，串行集成 |
| W6 | OPT-T8 clean | OPT-T9 | 线圈炮适配器依赖稳定的优化入口；仿真测试使用资源锁 |
| W7 | OPT-T9、OPT-T7 clean | OPT-T10 | 结果模型和 selector 依赖 Pareto 语义及物理评估结果 |
| W8 | OPT-T10 clean | OPT-T11 | 公共头文件、CMake 和双语文档集中集成，避免 API 写集冲突 |
| W9 | OPT-T11 clean | OPT-T12 | 端到端测试和 benchmark 串行执行，保证测量边界一致 |
| W10 | OPT-T12 clean | Whole-branch review | 一个高能力 reviewer 检查跨 Task 影响、数值和性能证据 |

W2、W3 和 W4 是主要并发窗口。每个窗口最多同时运行两个 implementer；其余槽位用于按完成顺序启动 reviewer。一个 reviewer 结束后才能为对应 Task 启动 fixer，不能因为有空槽位而提前启动依赖未完成的 Task。

### 10.3 每个 Task 的 subagent 生命周期

1. Controller 读取本计划对应 Task，记录当前 base commit，并生成只包含该 Task 要求的 brief 文件。
2. Dispatch 一个 fresh implementer，明确 brief 路径、允许写集、禁止写入、TDD 要求和 report 路径。
3. Implementer 完成 TDD、self-review，并创建唯一 Task commit；不得把无关格式化或其他 Task 修改混入 commit。
4. Controller 为该 Task 生成 review package，dispatch fresh reviewer 检查 spec compliance、代码质量、数值风险、写集边界和测试证据。
5. Reviewer 若发现 Critical/Important，停止该 Task 后续集成，dispatch 一个 fixer 在原 Task 分支修复并 amend 原 commit；fixer 必须重新运行覆盖测试并写入 report。
6. Reviewer re-review clean 后，Controller 更新 commit ledger 和 progress ledger，释放该 Task 的 agent/worktree。
7. 一个 Block 的所有 Task clean 且 Gate 通过后，才开启下一波次。

### 10.4 TDD 强制流程

每个实现 Task 的每项行为都必须遵循 Red → Green → Refactor：

1. **Red：** 先写最小行为测试；运行该测试并确认因功能缺失而失败，而非编译错误或测试错误。
2. **Green：** 写满足当前失败测试的最小生产代码；运行目标测试及已有相关测试，确认通过。
3. **Refactor：** 只在测试全绿后整理命名、重复逻辑和边界封装；不得借重构引入未测试行为。
4. 进入下一个行为前，保留前一个行为的绿色结果；Task commit 前运行该 Task 的完整覆盖测试。

TDD 测试写集属于对应 Task，不单独拆成“事后补测试” commit。测试必须验证公共行为和不变量，避免只验证 mock 调用次数。测试若一开始通过，implementer 必须确认是否测试了已有行为；若是，则重写为真正缺失的行为测试。

### 10.5 Review 与并发失败处理

- implementer `NEEDS_CONTEXT`：Controller 补充 brief 后由同一 agent 重新执行，不切换到无关 Task。
- implementer `BLOCKED`：先诊断是上下文、任务过大还是计划矛盾；必要时拆分 Task，不重复发送相同 prompt。
- reviewer 报告“无法从 diff 验证”的跨 Task 条目，由 Controller 在集成前自行核实；确认缺口则回到 implementer/fixer。
- 并发 Task 发生接口冲突时，停止受影响波次，新增 integration Task；不得在 Controller 工作树临时拼接补丁。
- 一个 Task 的 reviewer/fixer 不得修改另一个 Task 的写集；跨写集修复必须形成新 Task。

## 11. 执行前检查清单

- [ ] 用户已确认设计文档和本执行计划。
- [ ] 当前 branch 与设计适用 branch 一致。
- [ ] `.git/index` 可写且可创建 lock 文件；不能提交时不得启动实现 wave。
- [ ] `cmake --preset ninja-debug`、构建和基础 `ctest` 能运行，依赖状态已记录。
- [ ] CUDA 是否可用已记录；GPU workload 的资源锁策略已启用。
- [ ] Controller 已为每个 Task 准备 brief、report 和 review package 路径。
- [ ] 每个 Task 的英文 commit message 已固定，不在实现过程中临时改名。
- [ ] 首个实现 Task 从失败测试开始，未先写生产代码。
