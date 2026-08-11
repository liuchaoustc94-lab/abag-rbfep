# 1MLC / 1BJ1 / 1CZ8 失败根因分析与 target-specific 拯救方案

生成日期：2026-08-05
状态：分析完成，待决策执行 Phase 1
关联文档：`docs/validation_status.md`、`docs/validation_target_summary/validation_target_summary.json`

## 0. TL;DR

> **2026-08-05 Phase 0 更新（颠覆性发现）**：1MLC 的"失败"主要是**报告合并伪影**——官方 validation 视图中 1MLC 全部 8 个值取自 quick preset（1rep/2λ/1ps），尽管 7/8 个 job 已有深采样结果。换用 best-available 值后 **R 从 −0.35 → 0.84，MAE 从 7.7 → 1.29**。详见 §4.1。拯救优先级因此重排：**P0 = 修正报告合并策略（偏向最深 protocol 值）**；质子化假设降级为仅解释 1BJ1/1CZ8 的 H90A/Q89A 极端离群点。

三个被排除靶点的失败**不是采样不足，而是系统性（模型/数据层面）误差**：

- 采样加深 8 倍（quick → 16λ/160ps/4rep），ddG 只移动 0.6–0.8 kcal/mol，而系统性偏差为 3.6–4.5 kcal/mol（1MLC）甚至 11–28 kcal/mol（1CZ8 极端离群点）。**继续堆采样是无效算力消耗，应立即停止对这三个靶点的 deep/ultra drain。**
- 1MLC 的 8 个抗体侧突变呈现**近乎恒定的 −3~−4 kcal/mol 偏移**；1BJ1/1CZ8（同一 VEGF 抗原、同一 15 位点 Ala 扫描面板）呈现**动态范围夸大**（热点预测过高、中性突变预测过负）。
- 最可疑根因：① `pdb2gmx -ignh -missing` 默认质子化态（无 pKa 计算、无 His 优化）；② AB-Bind 实验数据的检测限截断值（3.69/4.10 重复出现）污染基准真值。

## 1. 三个靶点的失败画像

| 靶点 | 体系 | pairs | raw R | MAE | sign acc | 失败模式 |
|---|---|---|---|---|---|---|
| 1MLC | D44.1 Fab – 溶菌酶 HEL，2.5 Å，结构完整无缺失 | 8 | **−0.35** | 7.7 | 0.75 | 恒定负偏移（−3~−4），排序错误 |
| 1BJ1 | Fab-12 – VEGF 二聚体，2.4 Å，H 链 138–143 环缺失 | 15 | 0.36 | 3.8 | 0.40 | 动态范围夸大 + 实验值疑似截断 |
| 1CZ8 | Y0317（1BJ1 亲和力成熟版）– 同一 VEGF | 15 | 0.46 | 4.7 | 0.33 | 同 1BJ1 + 极端离群点（H90A/Q89A 偏差 11–28） |

**重要事实：1BJ1 与 1CZ8 不是独立样本**——同一 VEGF 抗原、完全相同的 15 个突变位点、两个亲缘抗体。holdout 80 pairs 中 30 个来自同一套 assay，它们主导了全量 R=0.20。

## 2. 已排除的假设（有证据）

| 假设 | 结论 | 证据 |
|---|---|---|
| 采样不足 | ❌ 排除 | 8× 采样 ddG 仅移 0.6–0.8；BAR stderr 已降至 0.32–0.50，远小于 3.6+ 的偏差 |
| apo 参考态错误（VEGF 单体 vs 二聚体） | ❌ 排除 | `_leg_mutated_chains`（stages.py:147-152）返回全部 antigen_chains；apo manifest 实测 `chains_retained: [V, W]`，拓扑含双链 |
| 输入结构质量差（1MLC） | ❌ 排除 | 2.5 Å、无缺失残基、prepare_qc 全干净（incomplete=[]、无 clash、未触发修复） |
| BAR 统计误差 | ❌ 排除 | 从未触发 stderr>10 阈值；deep protocol 下 stderr ≤0.5 |
| 结晶水/SO4 污染 | ❌ 排除 | prepare 阶段只保留 ATOM 记录，HETATM 全部剔除；实测 apo/complex 体系 HETATM=0 |

## 3. 幸存假设（按可疑度排序）

### H1. 质子化态错误（最可疑）
- `pdb2gmx -ignh -missing` 不带 `-his/-glu/-asp/-lys` 标志 → 全部按默认质子化（His 中性、Asp/Glu 酸性、Lys/Arg 碱性），**无 pKa 计算**。
- 1MLC 突变全在抗体 CDR（T31/S57/N32/N92），界面埋藏可显著移动邻近 His/Asp/Glu 的 pKa；一个错误的质子化态足以产生 3–5 kcal/mol 的恒定偏移——与观测到的均匀偏移模式吻合。
- 1CZ8/1BJ1 面板含 **H86A、H90A 两个 His→Ala 突变**；WT 中 His 质子化态错误会直接毁掉这两个点（H90A 正是 11–28 kcal/mol 极端离群点之一）。
- 1MLC 的恒定偏移也可能意味着 apo Fab 腿的环境与 complex 腿差异被系统性误估（同一质子化错误在两条腿中不对称出现）。

### H2. 实验数据截断/质量问题（1BJ1/1CZ8）
- 实验 ddG 中 3.69 出现 2 次、4.10 出现 2 次、0.0 出现 5 次（1BJ1）——典型的**检测限截断**特征。被截断的点无论计算多准都无法相关，会同时拉低 R 和抬高 MAE。
- 这是"合法排除"的依据：不是"靶点难"，而是"真值被截断"，应从相关性指标中剔除截断点并**明确标注理由**。

### H3. 固定电荷力场对界面 Ala 扫描的动态范围夸大（1BJ1/1CZ8）
- 热点突变预测 +9~+11 vs 实验 +2.7~+4.1，中性突变预测 −1~−4 vs 实验 0~+1：预测方差系统性大于实验方差，是缺失极化/松弛响应的典型表现。校准（side_linear 斜率）可部分吸收，但逐靶点残差仍大。

### H4. 次要结构因素（1BJ1/1CZ8 独有）
- H 链 138–143 缺失环 + C 端截短：远离界面热点（W79–94），影响有限，但属于已知固有缺陷。
- VEGF 二聚体中未突变的 V 链随 apo 腿一起采样：W 突变若影响二聚体稳定性会被计入 ddG（模型假设，非 bug）。

## 4. QC 失败模式备忘（供 resume 策略参考）

- 65 个 job 跨 10 个 root：pass 7 / warning 58 / fail 0。
- warning 主因：**repeat spread >1.0 kcal/mol**（几乎 100% 触发），其次 overlap<0.2（仅浅采样 root）。complex 腿方差略大于 apo 腿。
- **2026-06-25 prune 伪影警告**：旧 root 的 `qc_report.json` 已全被重置为 fail（"observed 0 dhdl files"），分析必须基于 `stages/qc.json` message、`reports/merged/benchmark_pairs.csv` 和 `rescue_summary.json`。
- 保留完整 repeat 级原始数据的 root 仅：validation_plan、verify_preeq、minibatch_1mlc。
- 1MLC minibatch pilot（16λ/160ps/4rep）证实：QC 从 "overlap 失败" 改善到 "仅剩 repeat spread 边缘超阈"，但 ddG 几乎不动。**QC 指标与科学准确性已解耦——repeat spread 收紧不代表 ddG 变准。**

## 4.1 Phase 0 执行结果（2026-08-05）

产物目录：`runs/analysis/target_rescue_20260805/`（`official_pairs_3targets.csv`、`phase0_a1_a3_metrics.json`、`phase0_a3_1mlc_best_available.json`、`job_level_best_pairs.csv`）

### A3 结果：1MLC 官方指标被 quick-plan 值污染（最重要发现）
- 官方 `predict_pairs_calibrated.csv`（side_linear, 80 pairs）中，**1MLC 全部 8 个 raw 值来自 `abbind_1mlc_core_v1` quick 批次**（1rep/2λ/1ps），例如 n92a=+32.13（quick 垃圾值；minibatch 深采样为 −4.86）、s57v=−14.45（priority 为 −0.69）、t31a=+6.87（priority 为 +0.41）。
- 换用 best-available（每 job 取最深 protocol 值）后：**Pearson R = 0.84，MAE = 1.29 kcal/mol，残余偏移仅 −1.10**（官方视图：R=−0.35，MAE=7.69）。
- 1MLC 实质上**不是失败靶点**；"accepted 排除 1MLC"的必要性应重新评估。
- 1BJ1/1CZ8 的官方值大多已是 priority（8λ/20ps）级，污染主要影响 1MLC；但 g92a/q79a/i80a 等仍有更深的 rescue 值未被采用（如 1cz8-q79a 官方 −4.21 vs deep −2.33，exp 0.65）。

### A1 结果：实验数据截断确认（1BJ1/1CZ8）
- 1BJ1：0.0×6、3.69×3、0.82×2（15 值中 11 个重复/零值）；1CZ8：0.0×5、4.1×2、1.06×2。**典型的检测限截断/取整**，截断点无法参与相关性评估，应标注 censored 后剔除。
- 1MLC 实验值无重复，真值质量良好。

### 幸存问题（真实物理失败点）
- **H90A**：1BJ1 −11.03 / 1CZ8 −14.22（exp 均 0.0），两个抗体同一位置一致失败 → His 质子化假设仍适用于此点。
- **Q89A**：1BJ1 −17.60 / 1CZ8 −26.76（exp 1.75/1.06），同样两抗体一致失败 → 需单独检查该位点的杂化拓扑/局部几何（Q→A 不应产生如此能量，疑似 per-position 构建问题而非采样问题）。
- 这两点（4 个 job）是 1BJ1/1CZ8 剩余误差的主要贡献者；剔除后两靶点 MAE 大幅下降。

### A2 结果：PROPKA 质子化审计（2026-08-05，propka 3.5.1 @ tmp/propka_venv）

产物：`runs/analysis/target_rescue_20260805/propka/`（complex + apo 腿各 3 套 .pka）

**发现三个质子化耦合结合位点（complex↔apo pKa 位移 ≥2.5），与失败模式一一对应：**

| 位点 | complex pKa | apo pKa | ΔpKa | 后果 |
|---|---|---|---|---|
| 1BJ1/1CZ8 **W:H90** | 2.66 / 2.37 | 6.16 / 6.17 | **−3.5 / −3.8** | apo 腿中 H90 ~13% 质子化带电，pipeline 双腿均默认中性 → 解释 H90A 灾难性失败（−11/−14 vs exp 0） |
| 1MLC **H:E50**（CDR H2，邻接 S57） | 4.85 | **7.41** | −2.56 | apo Fab 中 E50 ~70% 质子化（中性），pipeline 默认带电 → 给全部 8 个抗体突变引入近似恒定偏移 |
| 1MLC **H:E35**（CDR H1，邻接 T31） | 2.42 | 5.15 | −2.73 | 同上，量级较小 |

- 1MLC 两个异常位点恰好都在突变所在的 CDR 环上 → 解释 best-available 后仍残余的 −1.10 kcal/mol 均匀偏移。
- 1BJ1/1CZ8 的 H86 中性边界（pKa 6.3，双腿一致）→ 影响小。
- **Q89A 质子化无法解释**（Gln 不可滴定）：−17.6/−26.8 的灾难值疑似 pmx 杂化拓扑/构建问题，需单独逐点检查。

**机理结论：三个靶点的真实物理瓶颈是"质子化耦合结合"（protonation-coupled binding）——实验 ddG 包含质子释放/摄取贡献，固定质子化 FEP 原理上无法完全复现；但可通过修正参考态质子化（apo 腿 E50→GLH、H90→HIP 敏感测试）大幅逼近。**

### 修正后的优先级
1. **P0：修正验证报告的合并策略**——每 job 取最深 protocol（或最高 QC）的值，禁止 quick 值覆盖深采样值；重新生成 validation target summary 并重评 accepted 排除清单。预计 1MLC 可从排除清单移除，全量 holdout R 显著改善。
2. **P1：质子化修正敏感测试**（GPU 小规模）：
   - 1BJ1/1CZ8 h90a：apo 腿 W:H90 改用 HIP（带电）重跑，验证灾难值是否收敛；
   - 1MLC 全 8 job：apo Fab 腿 H:E50 改用 GLH（中性）重跑，验证残余 −1.1 偏移是否消除；
   - 实现方式：prepare 阶段输入 PDB 重命名残基（GLU→GLH、HIS→HIP），pdb2gmx 会识别 amber 质子化变体名，无需改 pipeline 代码。
3. **P1：Q89A 构建检查**——检查 pmx 杂化拓扑/突变后结构是否有异常（两抗体同位点一致失败，疑似 per-position 构建缺陷）。
4. **P2：截断感知指标**（censored 标注）。
5. **长期**：考虑在 prepare 阶段集成可选 PROPKA 步骤（环境依赖：tmp/propka_venv 已验证可用）。

## 4.2 P0 实施记录（2026-08-06）

**已完成的代码修复**（全部测试通过）：

1. `benchmark.py::_merged_job_priority()`：采样深度（`repeats×λ×production_ps`）从第 10 级提升到第 8 级，位于 repeat spread / BAR stderr 之前。从此同一 job 的深采样行在同 QC/完成度等级下必然击败浅采样行。`merge_strategy.winner_priority` 元数据同步更新。此修改使代码与 `benchmarks/ab_bind/README.md`（L635-638 早已声称 effort 优先）一致。
2. `report_calibrated_validation.py::DEFAULT_EXTRA_PLAN_ROOTS`：新增 `abbind_core_v1_validation_plan`（此前 1MLC 深采样值所在的主验证 lane 不在合并范围内）。
3. `calibrate_ab_bind_plan()`：`predict_pairs_calibrated.csv` / raw pairs 新增 `source_plan_root` 列，最终产物可追溯每个值的来源 root。
4. 测试：`test_report_ab_bind_plan_prefers_lower_bar_stderr_rows_before_sampling_effort`（锁定旧 stderr 优先行为）改写为 `test_..._prefers_sampling_effort_before_bar_stderr`，断言深采样行胜出。

**best-available 修正视图**（`runs/analysis/target_rescue_20260805/phase0_best_available_3targets.json`）：

| 靶点 | 官方 R | best-available R | 官方 MAE | best-available MAE | 升级 job 数 |
|---|---|---|---|---|---|
| 1MLC | −0.349 | **+0.840** | 7.69 | **1.29** | 7/8 |
| 1BJ1 | 0.360 | 0.354 | 3.78 | 3.75 | 3/15 |
| 1CZ8 | 0.458 | 0.455 | 4.74 | 4.59 | 3/15 |

1BJ1/1CZ8 的官方值本就大多是 priority 级，修正合并对它们影响很小——它们的剩余误差来自 H90A/Q89A 物理失败点和实验值截断，见 §4.1 A1/A2。

**⚠️ 数据状态事件（需知晓）**：
- 2026-06-25 的 prune 删除了 fit 侧（calibration/development split）job 的 BAR 原始数据，导致官方报告**自 6-25 起无法完整再生**（`fit_pair_count=0`，`insufficient_fit_pairs`）。6-23 的校准产物（`calibrations/.../model.json`、`predict_pairs_calibrated.csv`、`summary.json`，含 accepted R=0.607 的全量指标）仍完好。
- 2026-08-06 本次执行刷新脚本时，canonical `calibrated_validation_summary.json`（6-23 版）被覆盖为 `insufficient_fit_pairs` 状态——这反映了 prune 后的真实数据状态，但意味着**该文件的 6-23 版内容已丢失**；其关键指标保存在 `docs/validation_target_summary/validation_target_summary.json` 与 `calibrations/.../summary.json` 中。
- 要让官方视图真正吃到 P0 修复的红利，需要重跑被 prune 的 fit 侧 job（或接受 6-23 校准产物为冻结基线 + 用 best-available 视图做靶点级修正）。

## 4.3 Fit 侧重跑（2026-08-06 启动）

- **背景**：2026-06-25 prune 删除了全部非文本产物（dhdl.xvg/gro/tpr/xtc 等，共 1.39 TB），fit 侧 11 个校准 pair 与 validation 侧 80 个 pair 的原始数据全部丢失，官方报告无法再生。唯一幸存深采样数据：minibatch_1mlc root（3 个 1MLC job，16λ/160ps，384 个 dhdl 文件完好）。
- **方案**：新建 plan root `runs/benchmarks/abbind_core_v1_refit_priority_20260806`，全部 18 个 batch / 316 个 core_v1 job，统一 `validation_priority` preset（8λ/3rep/20ps，effort=480），GPU 0+1、2 workers、自动 resume 循环。运行脚本：`benchmarks/ab_bind/run_refit_priority_20260806.sh`（nohup  detached，日志在 plan root 的 `logs/refit_runner.log`）。
- **与合并修复的协同**：refit 行 effort=480 > quick(2)/validation(16)，会自然成为合并 winner；minibatch 行 effort=10240 仍为 3 个 1MLC job 的 winner（修复按预期工作）。该 root 已加入 `DEFAULT_EXTRA_PLAN_ROOTS` 首位。
- **踩坑记录**：`plan-abbind --runs-root` 必须传**绝对路径**；相对路径会导致 stage 脚本以 workdir 为 CWD 解析失败（exit 127，mutate 全灭）。
- **预计耗时**：数天（316 jobs × ~1.2 ns，2×RTX 4090）。完成后执行 §5 Phase 2 的报告再生。

## 4.4 关键点重跑 + 质子化敏感测试 pilot（2026-08-06 启动）

全量 316-job refit 已按用户决策**停止**，改为 35 个关键 job 的定向重跑（~1 天）：

**Root 1：`abbind_keypoints_baseline_20260806`（23 jobs，priority preset）**
- 1MLC ×8（官方视图中的 quick 垃圾值的正名对照）
- 1BJ1/1CZ8 的 w-q89a + w-h90a ×2（4 个灾难离群点的可重复性验证）
- 11 个 fit pair job（恢复校准拟合能力，job 清单取自 6-23 fit_pairs.csv）

**Root 2：`abbind_protonation_pilot_20260806`（12 jobs，同 preset，质子化变体）**
- 1MLC ×8：apo 腿 H:E50 改名 GLU→GLH（中性化；PROPKA apo pKa 7.41）
- 1BJ1/1CZ8 q89a+h90a：apo 腿 W:H90 改名 HIS→HIP（带电；apo pKa 6.16 vs complex 2.66）

**质子化注入方法（已探针验证）**：plan 跑到 prepare → 改 apo 腿 `input.pdb` 残基名 → resume。pmx 正确识别 GLH（中性 Glu 拓扑含 HE2）；HIS→HIP 在突变位点也兼容——pmx 的 amber99sb-star-ildn-mut 库内置 `Z2A = Hip->Ala` 杂化残基。无需改 pipeline 代码。

**注意（机理修正）**：严格的质子化耦合结合热力学估计显示，pKa 6.16 的 apo H90 在 pH 7 的质子化系综贡献仅 ~0.1 kcal/mol——若 HIP 变体使 ddG 移动数 kcal/mol，说明问题在**固定质子化/构象耦合**（错误质子化态导致的结构畸变被采样锁定）而非平衡热力学贡献；若几乎不动，则 H90A/Q89A 失败应归因于杂化拓扑或力场端点效应。pilot 结果将直接裁决。

运行脚本：`benchmarks/ab_bind/run_keypoints_protonation_pilot_20260806.sh`；日志：`runs/benchmarks/abbind_keypoints_baseline_20260806/logs/pilot_runner.log`。

## 4.5 新靶点前瞻验证：1DVF（2026-08-06 排队）

**选择**：1DVF = D1.3 Fv（链 A,B）– 抗独特型抗体 E5.2 Fv（链 C,D），1.9 Å。AB-Bind 源数据但不在 core_v1，26 个单点突变；剔除 9 个电荷变化突变（V2.1 范围）后 **17 个电荷守恒 job**（实验 ddG 0.34–4.79，动态范围好）。A:Y49A 有两个测量值（2.05/1.75）取均值 1.90。

**案例定义**：`examples/real_cases/1dvf/`（system.yml / mutations.csv / experimental_ddg.csv），刻意与 core_v1 基准隔离，不污染官方指标。

**PROPKA 前瞻发现**：突变位点 **D:H33 是反向质子化耦合案例**——complex pKa 7.87（pH 7 下 ~74% 带电）vs apo pKa 5.09（中性），ΔpKa=+2.78。与 1BJ1 H90 相反：这里 default（中性 His）错在 **complex 腿**。形式化质子化耦合贡献估计 ~1.3 kcal/mol（1.42×[log10(1+10^0.87)−log10(1+10^-1.91)]），是可测量量级。另一个位点 A:H30 双腿均中性（ΔpKa −1.15 但 pKa 均 <6），作为阴性对照。

**运行**（pilot 完成后自动接力，GPU 0+1）：
- `runs/real_cases/1dvf_priority_20260806`：17 个 baseline job
- `runs/real_cases/1dvf_protonation_20260806`：仅 h33a，**complex 腿** D:H33 HIS→HIP
- 脚本：`benchmarks/ab_bind/run_1dvf_newtarget_20260806.sh`

**验证设计**：① 17 job baseline vs 实验 ddG → 管线真正的样本外 R/MAE（无任何校准）；② h33a baseline vs HIP 变体 → 前瞻检验质子化方案；③ 若 1MLC/1BJ1/1CZ8 回顾性结论在 1DVF 上重现，则质子化感知 protocol 可推广。

## 4.6 工程改进落地（2026-08-06）

**① prune 数据保留策略**（`tools/prune_runs_keep_summaries.py`）：新增报告引用保护（默认开启，`--no-protect-report-references` 可关）——凡被 `reports/calibrations/*/*.csv`、`reports/merged/plan_jobs.csv`（含 selections）、`real_cases/*/reports/*.csv` 引用的 job，其 `dhdl.xvg`（BAR/MBAR 原始输入）和 `bar/*.xvg`（QC 直方图）永久保留。真实 dry-run 验证：保护 646 个被引用 job 的 1,104 个文件（~1.0 GB），仍可删除 167 GB 轨迹类数据。测试：`tests/test_prune_runs_keep_summaries.py`（3 个用例）。**运维规则：先再生报告再 prune**，否则新跑完但未入报告的 job 不在保护范围。

**② MBAR vs BAR 对比**（`runs/analysis/target_rescue_20260805/mbar_comparison/`）：在 minibatch 深采样数据上，alchemlyb BAR 逐位复现 pipeline 的 gmx bar（实现正确性交叉验证通过），MBAR 与 BAR 差异 ≤0.22 kcal/mol（远小于 2.4–3.6 的系统性误差）。**结论：估计器不是误差来源，MBAR 不换入默认链路**，脚本保留为可选审计工具。

## 4.7 历史 backlog 三项改进落地（2026-08-06）

来自 `docs/abag_rbfep_context.md` 的优先级清单：

**① clash 分层**（structure.py）：`find_intra/inter_residue_heavy_atom_clashes` 的每条 clash 新增 `atom_a_class`/`atom_b_class`/`clash_class`（backbone_backbone / backbone_sidechain / sidechain_sidechain），issue 级新增 `clash_classes` 汇总；新增 `classify_clash_atom`/`classify_clash_pair` 单一事实来源，两个 partition 函数改为读它（行为不变）。纯增量字段，无下游破坏。

**② 真实失败样本回归**（`tests/test_real_structure_regressions.py`，纯 CPU）：
- 1YY9：23 个不完整残基**全部 sidechain-only**（0 个 backbone blocking）——锁定"当前代码可修复放行"的行为（历史 block 记录来自修复逻辑加入前的旧 run）；
- 2NZ9：WT 源结构无任何 pre-existing clash——H1064A 失败归因于突变后几何而非输入；
- 3BN9：CYS H140/H196 SG-SG 0.83 Å（二硫键邻近误报的典型）被正确分层为 sidechain_sidechain 且可修复。

**③ repair provenance schema 统一（最小版）**：新增 `empty_repair_summary()`——占位 schema 与触发后 schema 字段集完全一致（trigger/blocking/remaining 字段齐全，空列表填充），prepare 阶段两处占位已切换。遗留：strip 操作（terminal oxygen / sidechain strip）的动作级记录仍是缺口，列入后续。

回归：370/370 通过（新增 4 测试 + 2 个精确匹配测试同步更新字段）。

## 4.8 重大发现：`-ignh` 导致全管线杂化残基缺氢（2026-08-06）

**这可能是项目史上最重要的单个 bug。**

### 发现过程
Q89A 拓扑解剖（探索 #1）中发现：`1bj1-antigen-w-q89a` 的杂化残基 Q2A 在 mutant.pdb 里有 18 个原子，但 processed.gro / 最终 itp 里只剩 **9 个重原子**——所有氢（H/HA/HB1/HB2/HG1/HG2/HE21/HE22）和 alchemical dummy（HV）全部丢失。

### 机理
- pmx mutate 正确生成含氢杂化残基；
- 但 mutate 阶段的 `pdb2gmx -ignh -missing` 中 **`-ignh` 丢弃输入氢**，而杂化残基（Q2A/N2A/S2A/Z2A...）在 `aminoacids.hdb` 中**没有加氢规则**（hdb 只覆盖标准残基），氢永远无法重建；
- 后果：每个 job 的杂化残基都是缺氢的化学错误模型——Q2A 状态 A 电荷和为 **−1.31**（应为 0），状态 B −0.54，λ 两端净电荷差 +0.78e，缺 8 个 LJ 作用位点。

### 影响面
**所有历史 ddG 都受此影响**（prepare 注释显示 `-ignh` 本为"统一重建晶体结构氢"而加，但 AB-Bind 源结构实测 0 个氢原子，该理由不成立）。这统一解释了：动态范围夸大、极端离群点（Q89A −17.6/−26.8、H90A −11/−14）、以及"加深采样永远不收敛到实验值"。系统性残差如此之大约合理的部分原因是：缺氢缺陷对所有 job 一致，部分在校准和 ddG 相减中抵消。

### 修复（已落地）
- `structure.py` 新增 `strip_hydrogen_atoms()`（带 provenance 返回）；prepare 阶段对 prepared_input 统一剥氢；
- mutate 与 auto-repair 两处 `pdb2gmx` **移除 `-ignh`**；
- 探针验证（无 -ignh）：Q2A 18 原子✓、Z2A（带电 His→Ala）21 原子✓、S2A 12 原子✓、GLH 质子化变体 16 原子✓；
- 测试：test_planning 断言更新（`-ignh` 不得出现）+ `test_strip_hydrogen_atoms_*` 新增。

### 后续
- 旧 pilot roots 已清空，pilot + 1DVF 以修复后代码重启（2026-08-06）；
- 历史全部结果（含 accepted R=0.607 视图）都是在缺氢杂化下计算的，**修复后指标预期显著变化**，官方视图需在 pilot/1DVF 完成后整体重建；
- 此前"质子化假设"的优先级可能需要重排：先量化缺氢修复的影响，残余误差再归因子化。

## 4.9 收敛性研究 + 截断值文献核实（2026-08-06）

### 收敛性研究（产物：`runs/analysis/target_rescue_20260805/convergence/ddg_vs_effort.json`）

用 5 档协议在同批 job 上的历史 ddG（缺氢数据，收敛模式仍有指导意义）：

| 层次 | 结论 |
|---|---|
| quick/validation → priority | ddG 剧摆（中位 Δ5.18，最大 37.6）——浅层数据完全不可用，priority 是有证据的**最低**可用 preset |
| priority → 4-20× 更深 | ddG 漂移中位仅 **0.52**（最大 1.89）——priority 层已基本收敛；**deep/ultra rescue 是资源浪费**（4-20× 成本换 <0.9 kcal/mol） |
| 最深采样的 \|误差\| 中位 | 2.44 kcal/mol——收敛≠正确（缺氢 bug 主导残差） |

**对生产的指导**：priority（8λ/3rep/20ps）作为默认 preset 有了数据支撑；停止 deep/ultra 式深采样救援；QC 的 repeat-spread 阈值应在修复后数据上重新标定（旧阈值是与缺氢数据共同演化的）。

### 截断值文献核实（详见 known_issues.md ISSUE-003）
1BJ1/1CZ8 的 3.69/4.10 实测为 SPR 检测饱和上报值（500×/1000× 倍数），同一方法学体系明确标注 "NB, >100-fold, cannot be precisely quantitated"。**截断实锤**，合法剔除有文献依据。

## 4.10 修复后第一批完整数据（2026-08-09/10，决定性验证）

pilot 35 jobs + 1DVF 18 jobs 全部完成（3 个失败 job 已重试）。**这是项目第一批化学完整的 ddG。**

### ① ISSUE-001 修复效果：灾难点平均误差 14.9 → 3.3 kcal/mol ✅

| job | 缺氢(broken) | 修复(fixed) | 实验 |
|---|---|---|---|
| 1bj1-w-q89a | −17.60 | −3.94 | +1.75 |
| 1bj1-w-h90a | −11.03 | +3.23 | 0.00 |
| 1cz8-w-q89a | −26.77 | +6.35(var) | +1.06 |
| 1cz8-w-h90a | −14.22 | −0.89 | 0.00 |

**策略裁决：`-ignh` 修复对极性残基灾难点完全有效。**

### ② 1DVF 样本外验证：R=0.621, MAE=3.58, sign_acc=0.76（n=17，无校准）
- 首个化学完整管线的前瞻性验证，达到抗体界面 FEP 的文献级样本外水平。
- 剩余 6 个大误差（y102a/w52a/y49a/q104a/y32a/c-y49a）全部为 **Y/W 大芳香侧链删除**，杂化化学完好但 overlap 0.08-0.26、spread 1.8-3.1——**真正的采样受限类**，与极性灾难点（化学 bug）机理完全不同。

### ③ 1MLC：MAE 1.29→1.19 略改善，R 0.84→0.47（n=8 小样本噪声 + n32y/t31a 离群）
- 缺氢数据下 t31a 的"完美"（0.41 vs exp 0.45）实为 bug 误打误撞；修复后 t31a=2.25。
- 残差分布从 −4.9~+0.9 收窄到 −1.6~+3.1，排序尚不稳定，需更大样本评估。

### ④ 质子化变体：混合结果 + 一条重要边界教训
- 1MLC E50→GLH：n92a 改善（0.37→−0.87, exp −1.25）、s57v/n32y 改善，但 t31a/t31w/n32g 变差——**未证实**。
- **H33 HIP 教训（1DVF）**：PROPKA pKa 7.87（74% 带电）的 H33，纯 HIP 态 FEP 给出 ddG=−25.7（baseline 中性 −0.02，exp +1.88）。**纯质子化态 FEP 只在 pKa  decisively 偏离 pH（population >95%）时有效**；pKa 6-8 的中间态不能用纯态近似，默认中性反而更接近。这给 PROPKA 集成定了规则：仅当 |pKa−pH|>2 时才应用变体。

### ⑤ 误差结构的重构（方法论结论）
修复后误差分两类，需要不同对策：
1. ~~极性残基灾难~~ → 已被化学修复消灭（ISSUE-001）
2. **大芳香删除（Y/W→A）采样受限** → 需要更多 λ 窗口（12-16）而非更多 repeats——可做成按突变类型的自适应协议

### 遗留
- fit pairs 9/11（3ngb equilibrate 失败=盒子爆炸；1cz8-q89a sample=cuFFT 瞬时错误；3be1-y33a 未启动）——已重试
- 1DVF 的 Y/W→A job 可加 λ 窗口重跑验证自适应协议
- 用 fit pairs 重建校准（待 11/11）

## 4.11 腿级误差归因（2026-08-11，n=41 修复后 job）

产物：`runs/analysis/target_rescue_20260805/leg_error_attribution.json`

| 预测因子 | 与 \|ddG 误差\| 的相关性 |
|---|---|
| apo 腿 BAR stderr | +0.36 |
| apo 腿 repeat spread | +0.29 |
| complex 腿 repeat spread | +0.27 |
| overlap（负相关） | −0.17 ~ −0.27 |

- 所有 QC 指标与误差的相关性都**弱到中等**——现有 QC 体系对 ddG 准确性的预测力有限（阈值重标定的直接动机）。
- 高误差 job 的腿主导性混合（top-8 中 4 复杂腿 / 4 apo 腿），**不存在"界面腿系统性更差"的模式**——残差不是界面描述特有的问题。
- 值得注意：t31w 变体 complex 腿 spread=15.3（GLH 质子化变体引入的扰动远超采样噪声），质子化变体需要更谨慎的评估。

## 4.12 即时验证三则（2026-08-11，纯分析）

**① λ 布置分析（minibatch 16λ overlap 矩阵）**：相邻窗口 overlap 极均匀（0.216–0.283，中位 0.233），最弱点在 λ 0.53–0.60（vdW morph 区）但仅略低。**结论：均匀布置即可，收敛瓶颈在窗口数量而非位置**——支持 A2 的按数量自适应策略，否决非均匀 λ 布置的复杂性。

**② 1MLC 恒定偏移假设终结**：修复后 per-job 偏移 = [−0.75, −1.15, +1.8, −0.1, +0.53, +0.47, +3.12, +1.62]，均值 +0.69 但**不再均匀**（此前 best-available 为恒定 −1.10）。残差是逐位点的物理/采样问题，无单一参考态错误残留——质子化统一解释的必要性进一步降低。

**③ NPT 平衡充分性验证**：n92a apo 腿密度 1015.5 kg/m³（TIP3P@310K 正常值），前/后四分之一段漂移仅 +0.5——**20 ps NPT 充分**，无需延长。

## 5. 拯救方案

### Phase 0：纯分析（不跑 MD，数小时内完成）

| # | 任务 | 方法 | 判读 |
|---|---|---|---|
| A1 | 实验数据审计 | 核对 AB-Bind 原始表 3 靶点全部行，标记截断值（重复出现的相同 ddG）与重复测量；给出截断点剔除后的 per-target R | 若剔除截断点后 1BJ1/1CZ8 R 显著上升 → H2 成立，走"合法排除"路线 |
| A2 | 质子化审计 | 对 3 个 WT 复合物 + apo 跑 PROPKA（本地 .venv 无 propka 则在隔离 venv 安装），列出每个突变位点 8 Å 内 His/Asp/Glu/Lys 的预测 pKa 与默认指派差异 | 找出"默认指派 ≠ PROPKA 预测"的残基 → Phase 1 的修正目标 |
| A3 | 恒定偏移检验 | 对每个靶点做去均值后的相关性（per-target demeaned R） | 若 1MLC 去偏移后 R 很高 → 排序正确、仅参考态偏移，H1 可能性大；若仍低 → 逐点物理错误 |

### Phase 1：定向修正重跑（GPU，数天）

| # | 任务 | 范围 | 成功判据 |
|---|---|---|---|
| B1 | 质子化修正重跑 | 用 PROPKA 结果对 pdb2gmx 加 `-his/-glu/-asp` 显式指派，重跑 1MLC 全 8 job + 1CZ8/1BJ1 的 H86A/H90A（complex+apo 两腿，minibatch 级 protocol） | 1MLC 恒定偏移消除（MAE<2）；H90A/H86A 离群消失 |
| B2 | 截断感知指标 | 在 `report_calibrated_validation.py` / target summary 中实现 censored-value 标注与剔除（理由明确记录为 "experimental detection limit"，区别于 "bad target"） | 指标报告含 censored 标记，accepted 视图可审计 |
| B3 | 停止无效采样 | 关停/归档三靶点的 deep/ultra/rescue drain 队列（证据：采样收益已饱和） | 释放 GPU 给 B1 |

### Phase 2：决策门

- **B1 成功** → 把质子化步骤（可选 PROPKA stage）固化进默认 pipeline（prepare 阶段），全基准重跑受益面评估。
- **B1 失败（偏移依旧）** → 判定为力场/终点效应限制（H3），文档化结论，三靶点正式标注 `model_limited` 而非简单排除，转向 B2 的合法剔除 + 校准模型吸收。
- 1BJ1/1CZ8 的非独立性应在 benchmark 统计中永久标注（cluster-level 权重或合并为一个 assay cluster），避免 30/80 pairs 的单一 assay 主导全量指标。

## 6. 执行优先级

1. **A1 + A3**（今天可做，纯 pandas）
2. **A2**（PROPKA，需隔离环境安装）
3. **B3**（立即停损）
4. **B1**（A2 出结果后启动）
5. **B2**（与 B1 并行）
