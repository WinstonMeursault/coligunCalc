# 审计发现

## 初步确认

- `NumericalModel.md` 的核心 CPU 方程大部分有对应实现，但实现并非“完全覆盖”：`OptimizationLevel::LookupTable` 与 `OptimizationLevel::Reference` 在多级仿真运行时没有按文档区分自感计算模式；组件已在构造时决定自感路径，`opt_level` 实际只影响互感距离/GL策略。
- GPU 持久化 kernel 已进入三个类的代码，但 API 文档仍将 `SimBatch` 描述为逐 pair kernel 且“计划实现持久化”，与源码不符。
- GPU 单级实现保留 `dbl_buf_active_`/`dbl_M1_`，但没有真正完成独立的双缓冲流程；该成员会造成结果时序难以解释，且文档声称双缓冲仍计划中。需作为“未完成优化/实现风险”报告。
- `SimBatch` 构造函数没有验证 `num_sims > 0`、`num_sims <= backend.max_batch_sims`，也没有验证 `set_excitations` 已对所有 simulation 配置；API 文档却把这些写成约束。
- `SimBatch` 实际不实现 `TerminationPolicy.enable_velocity_check`，而单级/多级实现实现了该条件，属于行为不一致。
- `API.md` 与 `API_cn.md` 在 SimBatch 持久化描述、测试数量（API 末尾仍写 17）等处已确认不一致/过时，需要同步修正。
- `GpuAdaptor` 的公开 API 在 API 文档中列出，但 `coilgun_cuda.hpp` 和普通 `coilgun.hpp` 是否包含它需核对；其 accessor 的参数/所有权/调用前置条件需要补全。

## 复核后确认

- CPU 构建成功；配置实际包含 18 个 CTest 测试，而非 README/API 旧文档中的 17 个。
- CPU `test_multi_stage_sim` 通过。GPU 单级 fallback 测试通过，但耗时约 73 秒；完整 CTest 在 120 秒总工具超时前运行到 GPU 多级阶段，不能据此宣称全套通过。
- `test_gpu_batch` 在 90 秒诊断运行中被 SIGTERM，说明该测试仍然过慢或存在运行时问题，尚未证明 batch 目标达成。
- 设计文档 §1.4 要求 `batch_id` 协议，但实现仍使用 doorbell-only 协议；这不是设计的完全落地。当前源码的 `persistent_kernel.cu` 中，kernel 在启动后直接等待各 doorbell，主机每轮重置/置 1 的时序没有 batch generation token，存在跨轮次竞态风险。
- 设计文档 §1.2/§1.9 要求 `PersistentBuffers` 包含 `batch_id`，源码 `persistent_kernel.cuh` 仅有 `active_pairs/shutdown/doorbell/seps/results`，明确属于设计缩水。
- 设计文档 §1.4 期望 SimBatch 外层按 sim 触发，但 persistent 缓冲区只有一份；当前实现确实串行复用它，功能上可行但不等于设计图示的 per-batch protocol。
- `GpuOptLevel::Aggressive` 的源码 kernel 已存在 FP32 变体，但没有发现对应精度/性能测试；设计文档要求的 Standard vs Aggressive 自一致性和 CPU Reference <0.1% 验证未落地。
- API 文档已补充 persistent 默认行为、fallback、双缓冲未实现、SimBatch 正数约束和 18 个测试，但设计文档与源码的 batch_id 差距仍需在最终报告中明确列为未完成项。
