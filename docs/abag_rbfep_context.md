# abag-rbfep 项目合并上下文

本文档合并以下两个历史会话的关键上下文，作为后续在本项目中继续工作的统一背景：

- `019e91de-a12b-7333-9ee4-e559e1d8c50a`：RBFEP 项目规划与架构设计
- `019eb5fb-0863-79c2-bd70-0ba71a9d11be`：项目深度扫描、验证、失败 case 和优化策略

## 项目目标与边界

`abag-rbfep` 是独立的抗体-抗原 RBFE 软件，运行目录为 `/mnt/data/liuchao/abag-rbfep`，不应在运行时依赖旧 `platform` 或 MM/GBSA 项目。

- V1：single-point mutation ddG
- V2：same-side double-point mutation ddG
- V2.1：charge-changing 与 cross-side mutation
- GROMACS 是动力学模拟底座
- `vendor/pmx` 负责突变、hybrid residue/topology 和 alchemical setup
- `gmx bar` 是正式 BAR 聚合主链路

## 工作流与代码结构

标准阶段为：

```text
ingest → prepare → mutate → build_legs
       → equilibrate → sample → bar → qc → report
```

- `src/abag_rbfe/`：CLI、计划、执行、GROMACS、QC、报告
- `src/abag_pmx/`：抗体-抗原突变逻辑
- `vendor/pmx/`：vendored pmx
- `benchmarks/ab_bind/`：AB-Bind 计划、watcher、rescue、calibration 和报告脚本
- `docs/`：架构、验证和状态文档
- `runs/`、`tmp/`：生成的运行产物，不作为源码维护

## 已完成的主要工作

- 建立独立项目和 `abag-rbfe` CLI
- 实现 system prepare、mutation validate、batch plan/run/report
- 支持 AB-Bind benchmark、split、calibration、rescue 和 merged report
- 修复多类 pmx insertion-code、forcefield、alchemy 及 `GMXLIB` 边界问题
- 增加 sidechain gap、clash、topology 和真实结构修复回归测试
- 增加 GPU 分配、stale job 恢复、target-specific watcher 和队列状态报告
- 报告链在校准汇总为空时可回退到 `validation_target_summary.json`
- 历史全量回归曾达到 `363 passed`

## 历史验证结论

某个历史快照中，accepted target-filtered Pearson `r` 为 `0.6073390390160122`，约 42 个 pair，排除 `1MLC`、`1CZ8`、`1BJ1`。该数值不是实时保证；后续引用指标前必须重新读取当前 JSON/Markdown 产物。

问题靶点的历史分类如下：

- `1MLC`：少数 outlier、采样方差和 overlap 不稳定
- `1BJ1`、`1CZ8`、`3NPS`：校准反转或 ranking instability
- `3HFM`：弱信号和采样不足，不只是简单离群点问题

## 后续优先级

1. 统一 prepare/mutate repair provenance JSON，记录原始问题、修复动作和残留问题。
2. 固化 `1YY9`、`2NZ9` 等真实失败样本回归。
3. 完善 insertion-code 映射或残基重编号，扩大真实样本覆盖面。
4. 对 inter-residue clash 分层，区分 backbone 和 sidechain 情况。
5. 针对 `1MLC`、`1BJ1`、`1CZ8`、`3HFM` 优先测试 targeted sampling、window spacing、约束释放和重复采样。

## 操作约定

使用项目 `.venv` 进行验证；以生成的验证报告和实时队列状态为准，不仅依赖 README。长时间运行时保持 GPU 并发约束，先分类失败原因，再决定是修输入、修报告、补采样还是调整校准。
