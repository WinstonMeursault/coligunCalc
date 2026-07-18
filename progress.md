# 进度日志

## 2026-07-18

- 已加载审计、规划和项目规则。
- 已发现仓库已有审计报告及近期文档修复提交，需以当前源码重新验证，不能直接信任历史结论。
- 已逐个检查公开头文件，API.md/API_cn.md 的章节结构一致，覆盖 physics、components、simulation、CUDA API。
- 修复了 persistent kernel 缺少 batch_id 的同步协议、固定索引空间错配、单级 stale double-buffer 路径、SimBatch 构造/配置验证和 velocity-decay termination 缺口。
- 同步修正 README 测试数量、CPU OptimizationLevel 实际语义、SimBatch persistent 描述、GpuBackend 约束、GPU 已知限制和 GpuAdaptor 前置条件。
- 验证：CPU 核心 10 项测试全部通过；GPU persistent `test_gpu_vs_cpu_single`、`test_gpu_vs_cpu_multi`、`test_gpu_sim_batch` 全部通过；完整配置共 18 项。`test_gpu_batch` 因运行时间过长仍未在 90 秒诊断窗口完成，不能宣称全套通过。
