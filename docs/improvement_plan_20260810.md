# 修复后时代改进计划（2026-08-10）

背景：ISSUE-001（`-ignh` 杂化缺氢）已修复并验证（灾难点 MAE 14.9→3.3；1DVF 样本外 R=0.621）。
修复后误差结构清晰化为两类：已消灭的极性灾难（化学 bug）与真正采样受限的大芳香删除（Y/W→A）。
关联：`docs/known_issues.md`、`docs/target_rescue_analysis_1mlc_1bj1_1cz8_cn.md` §4.10。

## RID 增强采样评估（run_abag_rid.sh）

`run_abag_rid.sh` 是 gromacs-abag-mmgbsa 仓库的 **RiD-kit 0.6.4（Reinforced Dynamics）** 封装：
4 walkers × N iterations、CV 集合（≤64，距离 CV ≤32）、basin 选择（≤3）、
无偏验证段（pilot 2 ns / production 10 ns）、专用 `gromacs-rid` (gmx_mpi) 工具链。
设计目标是探索结合/解离构象景观（绝对结合自由能方向）。

**对本项目（相对 FEP ddG）的适用性评估**：
- 潜在价值：为 WT/MUT 端点态生成多样构象 ensemble → 让 FEP repeats 真正独立、覆盖多个 basin，
  针对的是我们观测到的 repeat spread 超阈（不同 repeat 落入不同局部极小）。
- 不匹配点：RiD 用自己的 gmx_mpi 与 CV boosting 流程，与 alchemical λ-dhdl 采集不在一个框架内，
  整合成本高（需要把 RID 采样嵌入每个 λ 窗口或做端点 seeding 协议）。
- **结论**：先走便宜的自适应 λ（Phase A2/B1）；若 Y/W 删除类在 16λ 下 overlap 仍 <0.2，
  再评估 RID/REST2 做端点 ensemble seeding。**暂缓，作为 Phase B 失败时的 fallback。**

## Phase A：协议与 QC 代码改进（纯 CPU，当天）

- **A1 QC 升级与重标定准备**
  - 杂化完整性 QC 从 informational 升级为 **blocking**（mutate 阶段不过则 block job）
  - 收集 32 个修复后 job 的 overlap/spread 分布，作为阈值重标定数据集（先出分布报告，不急着改阈值）
- **A2 自适应 λ 协议**
  - plan 阶段按突变类型自动选择 λ 窗口数：Y/W/F→X 大删除 → 16λ；其余 → 8λ（priority 基线）
  - 实现为 job 级 protocol override（`lambda_windows` 字段），最小侵入
- **A3 分步 λ 调度（coul-lambdas + vdw-lambdas）**
  - 修 H33 HIP 类电荷变化转化的 overlap 崩塌（0.064），也是 V2.1 前置
  - 设计：前半 λ 步只变电荷，后半只变 vdW（沿用 sc-soft-core）；默认开启
- **A4 production dt 0.001 → 0.002**（OBS-2，采样成本减半；与 A3 同一次重基线）
- **A5 全量测试回归 + 文档同步**

## Phase B：GPU 验证（2 GPU，预计 1-2 天）

- **B1** 1DVF 6 个 Y/W→A 失败 job 用新协议（16λ + 分步调度 + dt 0.002）重跑 →
  判据：overlap_min ≥0.2 且误差下降（当前 4.4-9.9 kcal/mol）
- **B2** H33 HIP 变体用分步调度重跑 → 判据：complex 腿 overlap ≥0.2（当前 0.064）
- **B3** fit pairs 补齐 11/11（3 个失败 job 重试中）→ 重建 side_linear 校准 + 官方验证视图

## Phase D：RID 端点系综试点（2026-08-13 启动）

**依据**：basin 实验 spread 2.7-3.3 > 1.0 预设门槛（§4.14①）。
**设计**：`run_rid_pilot_1dvf_20260813.sh`
1. RID pilot 探索 1DVF WT 复合物（A,B|C,D；3 walkers × 3 iters，GPU 0-2）→ ≤3 个 basin（无偏弛豫后构象）
2. basin 蛋白坐标导出 PDB → 作为新输入结构走标准 FEP 全流程（prepare→mutate→…），y102a + y49a 各 1 job/basin
3. 对比：RID-basin ddG 均值 vs 自由 MD basin 均值（13.6/−6.9）vs 实验（4.79/1.90）

**判读规则**：RID basin 的 ddG 均值显著更接近实验（|误差| 减半以上）→ 端点系综策略成立，推广到 Y/W 删除类；仍不变 → Y/W 误差定位为力场/水合模型，停止采样方向投入。

**踩坑记录**：run_abag_rid.sh 的 GPU 选择器硬编码 --count 4；绕过方式 = --dry-run 生成 config → 改 RID_WALKERS=3 + RID_GPU_IDS + **RID_EXECUTE=1**（dry-run 会写 0）→ 直接调 run_rid_pipeline.sh。

## Phase C：分析与扩展（与 B 并行）

- **C1** 腿级误差归因（complex vs apo）——用 32 个新 job 的 bar_summary
- **C2** ~~扩展样本外靶点准备~~ → 已完成（2026-08-11）：1AK4（15 job，CycA/HIV 衣壳）、1KTZ（15 job，TGF-β 受体/配体）、3K2M（5 job，Monobody/SH2）。结构全部干净（0 缺失 0 clash）、编号核对通过、电荷守恒过滤完成、PROPKA 审计显示三个靶点均无决定性质子化位点（按 \|pKa−pH\|>2 规则无需变体）——干净的前瞻验证集。案例定义在 `examples/real_cases/{1ak4,1ktz,3k2m}/`，runner `benchmarks/ab_bind/run_outofsample_expansion_20260811.sh` 已排队（等 Phase B 完成自动接力，35 jobs）。注意 3K2M 有 4.0×3 重复值（将被 censored 规则捕获）。
- **C3** PROPKA 集成设计（规则：|pKa−pH|>2 才应用变体；纯态边界教训见 §4.10④）

## 验收标准

1. A 阶段全部测试通过（含新增自适应 λ / 分步调度 / blocking QC 的测试）
2. B1 的 6 个 job overlap_min 全部 ≥0.2 且 MAE 下降
3. B2 overlap ≥0.2（证明分步调度修好了电荷变化转化）
4. B3 官方视图用化学完整数据重建，README/文档数字更新
