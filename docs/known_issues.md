# abag-rbfep 已知问题追踪

维护说明：新发现的问题按 ISSUE-NNN 编号追加；状态 ∈ {open, in_progress, fixed, verified, wontfix}。
关联分析文档：`docs/target_rescue_analysis_1mlc_1bj1_1cz8_cn.md`。

---

## ISSUE-001 [verified 2026-08-10] pdb2gmx `-ignh` 导致杂化残基缺氢（全管线化学错误）

- **严重度**：critical —— 影响所有历史 ddG 结果
- **发现日期**：2026-08-06（Q89A 拓扑解剖中发现）
- **现象**：每个 job 的杂化残基（Q2A/N2A/S2A/Z2A…）在 processed.gro 与最终 itp 中只剩重原子。例：Q2A 应为 18 原子，实际 9 个；状态 A 电荷和 −1.31（应为 0），λ 两端净电荷差 +0.78e。
- **根因**：mutate 阶段 `pdb2gmx -ignh -missing` 的 `-ignh` 丢弃输入氢；杂化残基在 `aminoacids.hdb` 中无加氢规则（hdb 只覆盖标准残基），pmx 精心放置的氢与 alchemical dummy 被丢弃后无法重建。
- **历史影响**：极性/带电 WT 突变平均 |误差| 5.97 kcal/mol（31% >5），疏水 2.38（12%）；误差 top-8 全部为极性残基（N/Q/H/T）。Q89A −17.6/−26.8（exp ~1）、H90A −11/−14（exp 0）等灾难点均符合此指纹。所有历史指标（含 accepted R=0.607）均基于缺氢杂化。
- **修复**：`structure.py::strip_hydrogen_atoms()`（prepare 阶段统一剥氢 + provenance）；mutate 与 auto-repair 两处 pdb2gmx 移除 `-ignh`。探针验证 Q2A=18✓ / Z2A=21✓ / S2A=12✓ / GLH=16✓。回归 371/371。
- **遗留**：
  - [x] ~~pilot/1DVF 修复后数据完成后，头对头量化偏差~~ → **已验证（§4.10）：灾难点 MAE 14.9→3.3；1DVF 样本外 R=0.621**
  - [ ] 新增 QC：杂化残基完整性校验（已实现 informational 级，见 OBS-5；可升级 blocking）
  - [ ] 官方视图重建（fit pairs 9/11，3 个失败 job 重试中）
  - [ ] 历史结论（1MLC 排除、校准模型参数）全部需要重新评估

## ISSUE-002 [open] 2026-06-25 prune 删除原始 FEP 数据，官方报告无法再生

- **严重度**：high（工程/数据管理）
- **事件**：`prune_runs_keep_summaries` 删除全部非文本产物（1.39 TB，含 dhdl.xvg/gro/tpr/xtc），fit 侧 11 pair + validation 侧 80 pair 原始数据丢失；canonical `calibrated_validation_summary.json`（6-23 版）于 2026-08-06 被刷新脚本覆盖为 `insufficient_fit_pairs`，旧内容丢失。
- **缓解**：关键指标保存于 `docs/validation_target_summary/validation_target_summary.json` 与 `calibrations/.../summary.json`（6-23 冻结）；prune 工具已加报告引用保护（默认开启，~1 GB 保 dhdl/bar 输出 vs 释放 167 GB）。
- **遗留**：
  - [ ] 用修复后 pilot 的 11 个 fit pair 重建校准
  - [ ] 运维规则固化：先再生报告再 prune

## ISSUE-003 [open] 1BJ1/1CZ8 实验值检测限截断污染基准真值

- **严重度**：medium（基准统计）
- **现象**：1BJ1 实验 ddG 0.0×6、3.69×4（G92A/I83A/M81A/R82A）、0.82×2；1CZ8 0.0×5、4.1×2（G92A/M81A）、1.06×2。
- **文献核实（2026-08-06）**：数值源自 Muller et al. 1998 (Structure 6:1153, 1BJ1) 与 Chen et al. 1999 (Y0317/1CZ8) 的 VEGF 丙氨酸扫描（SPR, 37°C）。换算亲和力倍数（RT ln, 298K）：**3.69≈500×、4.10≈1000×** —— 整百/整千倍是典型检测饱和上报值；同一方法学体系（Muller 1997 PNAS, PMC23789）明确标注强破坏突变为 "NB (non-binder), >100-fold, cannot be precisely quantitated"。交叉印证：1BJ1 R82A=3.69 vs 1CZ8 R82A=0.82（同一突变两个抗体差 4.5×，符合"至少一个值是饱和上报"）。0.0×6/×5 为低端"无可测效应"。
- **结论**：**截断/饱和上报实锤**，这些点应以 censored 标注后合法剔除（理由 = experimental detection limit，区别于"难靶点"）。
- **遗留**：
  - [x] ~~指标层实现 censored 标注~~ → 已实现 `reporting.py::flag_censored_experimental_values()`（同复合物内重复≥2次且等于最大值=饱和）；效果：**80-pair 全量 calibrated R 0.200→0.292**（仅剔 5 点）；产物 `runs/analysis/target_rescue_20260805/censored_flagged_pairs.json`
  - [ ] 官方视图重建时将 censored 剔除接入正式报告链
  - [ ] 1BJ1/1CZ8 同一抗原同一突变面板（30/80 pairs 非独立）→ assay cluster 加权

## ISSUE-004 [open] strip 操作无动作级 provenance

- **严重度**：low（可审计性）
- **现象**：`strip_terminal_oxygen_atoms`、`strip_sidechain_atoms_for_residues` 不记录剥除了哪些原子/残基（`strip_hydrogen_atoms` 已带 provenance，可作范式）。
- **遗留**：给两个 strip 函数加返回值并写入 prepare_qc / mutate_qc。

## ISSUE-007 [in_progress] dt=0.002 + 严格分步调度在柔性/插入体系不稳定

- **严重度**：medium（协议鲁棒性）
- **现象（2026-08-12）**：新协议下 9 个 job sample 阶段失败（SIGSEGV / Too many LINCS / exit 1），全部为**插入型突变**（G/A→V/Q/A/T，在紧密堆积的环区：1AK4 CypA 环 487-496、1VFB CDR-H1、2JEL）。
- **根因定位**：失败窗口 = 严格分步调度的**交界窗口**（coul=1, vdw=0）——插入的侧链成为**全电荷零 vdW 的幽灵原子**，在紧密环区引发 LINCS 爆炸。dt=0.001 补丁不能根治（验证了错误假设）。
- **修复（2026-08-12）**：分步调度改为**重叠式**——coul 前 1/2 窗口 ramp、vdW 后 2/3 窗口 ramp，coul 满时 vdW 已达 0.25，消除幽灵态。7 个失败 job 已清空 sample 重跑验证中。
- **遗留**：
  - [x] ~~重试结果验证~~ → 演变为**方向感知调度**（见下）
  - [ ] 方向感知调度下 7 个插入 job 重跑验证（排队中）
  - [ ] 评估 sample 阶段自动 dt 减半回退机制的残留必要性

**调度演进（2026-08-12 三轮迭代）**：严格分步（修电荷变化 overlap，B2 验证）→ 重叠分步（幽灵态缓解不足，失败窗 λ003→λ005）→ **方向感知分步**（插入 vdW 先行、删除 coul 先行——插入型突变的侧链必须先在口袋中“长大”再带电，否则交界窗口的全电荷零 vdW 幽灵原子在紧密环区必然 LINCS 爆炸）。

## ISSUE-006 [open] 全链路正确性审计记录（2026-08-06，ISSUE-001 后续）

对 FEP 采样链路的系统审计结论——**未发现其他代码级 bug**，以下是设计层观察项（非 bug，可优化）：

**已验证健康的环节**：
- repeat 独立性：每 rep NVT `gen-vel=yes` + 按 job/leg/rep 哈希的稳定种子 ✓
- production dt/nsteps 一致（`_steps` 正确传入 dt=0.001）✓
- BAR 输入 λ 顺序正确（lambda_000 零填充字典序）✓
- gmx bar 与 alchemlyb BAR 交叉验证一致（见 §4.6 MBAR 报告）✓
- 修复后杂化残基电荷闭合实测：Q2A 18 原子、chargeA=chargeB=0.0000 ✓
- genion `-neutral -conc 0.15` 标准 ✓

**设计层观察项（可按需优化，非 bug）**：
- OBS-1：单一 `fep-lambdas` 调度（coulomb+vdW 同步 morph，sc-coul=on 缓解）。分步调度（先电荷后 vdW）收敛更好，但属协议改进而非错误。
- OBS-2：production `dt=0.001` 配合 `constraints=h-bonds` 过保守——dt=0.002 是标准做法，**可把全部采样成本减半**；变更需重基线。
- OBS-3：window_relax 较短（priority 下 0.2 ps MD），极端 λ 窗口可能欠松弛；deep rescue 已部分缓解（2x scale）。
- OBS-4（已验证关闭 2026-08-11）：`gmx bar` 未设 `-b` 截断——minibatch 数据实测截去前 20/40 ps 后 dG 变化仅 −0.12/+0.03 kcal/mol，window relax 预松弛充分，无需截断。

**补充（2026-08-11）**：
- 3ngb-g54s "盒子爆炸"非 bug：22 nm 超长复合物触发大盒阶梯回退（26.9 nm 立方盒）按设计工作；优化注记：细长体系可在装箱前做主轴对齐/旋转，水盒体积可减 ~50%。
- censored 感知指标已实现并应用（ISSUE-003）。
- OBS-5（已验证关闭 2026-08-06）：vendored pmx 容忍性补丁（alchemy 缺原子跳过）在 ISSUE-001 修复后**不再触发**——新代码产物实测 Q2A 18/18 原子、14 个 B-state 显式映射、完整键连网络（19 bonds/60 pairs/38 angles/116 dihedrals）。补丁保留作为不完整晶体残基场景的安全网，行为由 test_pmx_alchemy_patch.py 锁定。**防回归守卫已上线**：`gmx.py::validate_hybrid_topology_integrity()` 检查杂化残基原子完整性 + A/B 状态电荷整数闭合，已接入 mutate_qc（informational 级，pilot 验证后可升级为 blocking）。

## ISSUE-005 [open] 文档指标口径不一致

- **严重度**：low（文档）
- **现象**：README accepted R=0.672/32 pairs vs `validation_target_summary.json` R=0.607/42 pairs；且两者都基于 ISSUE-001 的缺氢数据，修复后均需重建。
- **遗留**：修复后数据就绪时统一以 `validation_target_summary.json` 为唯一权威来源重新生成全部文档数字。
