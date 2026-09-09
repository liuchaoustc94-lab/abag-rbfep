# 突变类型路由表（用户决策 2026-08-26，已固化进规划层）

> 来源：P1/P2/P3 靶向调研（`fep_strategy_research_p123_20260826.md`）的用户裁决版。
> 实现：`src/abag_rbfe/routing.py::classify_mutation_route`，接入 `planning.build_batch_plan`
> （jobs.csv 增加 `route`/`route_tags` 列 + `routing_summary.json`）；测试 `tests/test_routing.py`。

| 突变类型 | 推荐路径 | 路由标签 | 实现状态 |
|---|---|---|---|
| 普通中性 Ala 扫描 | 当前 12λ 基线 | `baseline_12lambda` | ✅ 现有协议 |
| Pro/Gly、侧链插入、柔性环区 | 端点系综 + 局部 REST2/RID | `endpoint_ensemble_rest` | 🟡 端点系综已有（40ns+帧注入）；局部 REST2 待建（P1-a） |
| Y/W/M 大空腔删除 | RID + 水占据分析；必要时 GCMC | `rid_hydration` | 🟡 水占据分析待建（P1-b OPES）；GCMC 仅作概念参照 |
| 电荷变化 | DSSB；失败时评估 co-alchemical ion | `dssb` | ✅ DSSB 已强制路由；co-alchemical ion 仅作失败兜底评估 |
| His/Asp/Glu/Lys 邻域 | PROPKA 多质子化状态 + pKa 加权 | `protonation_multi` | 🟡 PROPKA 隔离环境已有；多态加权待建（P2-b） |
| 高重复离散 | 增加 basin，而非单纯增加 λ 窗口 | `basin_diversity`（tag） | ✅ 已收敛研究支持（加深采样阴性） |

## 优先级与裁决规则

primary 路由按优先级单选：**dssb > rid_hydration > endpoint_ensemble_rest > protonation_multi > baseline**；
tags 独立累积（如 K→D 既是 `dssb` 又带 `titratable_neighbor` 标签）。

## 与用户原表的两处偏差（实现时扩写，已在代码注释标注）

1. "Y/W/M 大空腔删除"实现为 **Y/W/F/M → A/G**：补入 F（我们的 P1 实测误差类是 Y/W/F→A）。
2. "柔性环区"无法仅从突变判定——环区/CDR 归属需要系统级注释（留了 `loop_residues` 钩子，CDR 注释接入后自动生效）。

## 运维注记

- 路由是**建议性标签**；目前只有 dssb 在协议层强制执行（既有行为）。其余路径在对应引擎（REST2/OPES 帧源、PROPKA 加权）落地前仅作标记与人工分流依据。
- leg_topology 被 protocol 显式 pin 住时，路由保留为 advisory 并打 `leg_topology_pinned` 标签。
