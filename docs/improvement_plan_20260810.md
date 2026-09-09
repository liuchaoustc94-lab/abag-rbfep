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

**RID 执行进展（2026-08-13）**：RID 探索完成——3 个 basin 通过三重门控（接触保持率 0.75-0.87、uncertainty 达标）；无偏弛豫构象已导出为纯蛋白 PDB（6874 原子）；6 个 basin-FEP job（3 basin × y102a/y49a）已入队。

**踩坑记录**：① run_abag_rid.sh 的 GPU 选择器硬编码 --count 4；绕过 = --dry-run 生成 config → 改 RID_WALKERS=3 + RID_GPU_IDS + **RID_EXECUTE=1**（dry-run 会写 0）→ 直接调 run_rid_pipeline.sh。② trjconv 的 Protein 选择必须给 tpr（给 gro 会静默输出全溶剂体系）。

## Phase G：Sampson 式 outlier 分类 + 校正（2026-08-19 实施+裁决）

**实现**：`src/abag_rbfe/outliers.py`（结构感知五参数分类器 + 单参数收缩/偏移/封顶校正）+ `benchmarks/ab_bind/run_outlier_correction_analysis.py`（LOCO-CV 驱动）+ 测试。

**分类结果**（94 个化学完整验证点）：
| 类别 | n | MAE |
|---|---|---|
| buried_aromatic_large_deletion（Y/W/F 大删除） | 22 | **3.46**（主系统性偏差类） |
| large_deletion（其他大删除） | 39 | 1.84 |
| none | 33 | 2.90 |

**校正裁决（LOCO-CV，诚实阴性-阳性混合结果）**：
| 指标 | raw | 收缩校正 |
|---|---|---|
| MAE | 2.589 | **2.106（−19%）** |
| 强效应 MAE | 5.353 | **4.526** |
| 三分类 BA | 0.402 | 0.416 |
| Pearson R | 0.523 | **0.425（−19%）** |

**关键认识**：Sampson 校正在我们数据集上**用相关性换 MAE**——因为 Y/W 大删除类跨越整个实验动态范围，既承载信号（驱动 R）又带系统偏差（驱动 MAE）；任何均匀收缩在压掉过冲的同时也压掉了合法分布宽度。Sampson 有效是因为他们的 outlier 是真异常，我们的是系统性力场偏差。

**处置**：
- **官方指标一律报 raw**；分类标签（buried_aromatic_large_deletion 等）作为标准标注进报告
- 收缩校正仅作为 **opt-in 生产评分层**（MAE 导向的实验规划场景），不进基准报告
- 系统性偏差本身走物理路线：双力场（E6）+ 端点系综（已验证减半）

### E6 双力场 consensus（2026-08-20 裁决：阴性为主）

9 个最难 case 的 charmm36m-mut 臂完成（amber 臂取自现有数据）：

| 指标 | amber | charmm | consensus |
|---|---|---|---|
| Pearson R | 0.704 | **0.802** | 0.755 |
| Spearman | 0.517 | **0.567** | 0.533 |
| MAE | 6.92 | **6.06** | 6.49 |

**结论**：① charmm 单臂在三指标上略优于 amber；② **consensus 居中但不超 charmm 单臂**；③ 两力场对 Y/W 大删除的系统性高估**高度共享**（y102a：14.7 vs 16.2；w52a：11.8 vs 12.9）——**动态范围膨胀是力场无关的物理问题（采样/水合/端点），双力场无法解决**。⚠️ 注意：本子集是 cherry-picked 的最难点，指标不能与全集直接比。

**意义**：排除了一个假设，确认了 Y/W 问题的解决必须走系综/水合物理路线（端点系综、GCMC），而非力场切换。

## Phase F：调研驱动的零成本落地批次（2026-08-19 实施）

来源：`docs/fep_strategy_research_20260813.md` §四 的 1/2/3 项。

### F1 Petrov 2024 四指南审计（de Groot 组电荷伪影权威答案）✅
| 指南 | 我们的现状 | 处置 |
|---|---|---|
| 电中性盒 | genion -neutral 固定行为 ✓ | 无需改 |
| 溶剂缓冲 >1nm | 默认 1.0nm（边界） | **dssb_charge preset 提到 1.25nm**（电荷变化最敏感）；默认 preset 保持 1.0（成本） |
| 盐 0.1-0.15M | 0.15 ✓ | 无需改 |
| 路径保总电荷 | 方向感知调度 + DSSB 净电荷不变 ✓ | 无需改 |

### F2 cinnabar/MLE cycle-closure 分析层 —— 暂缓
当前数据全是星形图（所有突变从 WT 出发），无闭合环 → MLE 无可压缩的冗余。**等 V2 双突变链（WT→A、WT→B、A→AB、B→AB 成环）数据到位后再接入 cinnabar**（OpenFE 生态 pip 包）。已记录在案。

### F3 误差纪律升级 ✅
- 已实现：三分类 balanced accuracy（favorable/neutral/unfavorable，±1 kcal/mol 带）加入 `_benchmark_metrics_from_pairs`（benchmark.py）
- 现状基线：validation 面板 BA=**0.379**（中性召回 0.52、强去稳定 0.62、强稳定化 0/1）——Sampson 参照 0.69，诚实差距在案
- 独立 repeat SD 作为不确定性金标准（已有），MBAR/bootstrap 单轨迹误差不作报告口径

## Phase E：下一轮策略（调研驱动，2026-08-13 立项）

来源：`docs/fep_strategy_research_20260813.md` 的优先级清单。在当前队列（Wave/RID/REST2）完成后按序推进。

### E1 naked-charge λ 校验器（零成本，代码）✅ 已实现 2026-08-14
- `stages.py::check_lambda_schedule_naked_charge()`：逐窗口检查每类 alchemical 原子（deletion/insertion）是否存在"电荷>0.5 且 LJ<0.05"的幽灵态；原子类存在性由突变侧链重原子差推断
- 已接入 build_legs 的 `lambda_plan.json`（`lambda_schedule_naked_charge` 字段 + `lambda_distribution`）
- 验收达成：测试证明能对 ISSUE-007 的旧严格调度（coul 先行的插入型）预先报警；方向感知调度通过；另发现删除侧在插入调度下的残余风险也可被标记

### E2 Gapsys soft-core A/B（零成本，2 job）🟢 已入队 2026-08-14
- 前置重构：soft-core 参数从硬编码改为 protocol 字段（`sc_alpha/sc_sigma/sc_power/sc_coul`，默认保持 Beutler 兼容），3 个 renderer 已切换，测试通过
- 新增 `benchmarks/ab_bind/protocol.sc_gapsys.yml`（sc_alpha=0, sc_sigma=0.3）
- A/B job（队首）：`abbind_sc_gapsys_ab_20260814` 的 1bj1-q89a（极性删除）+ 1bj1-g88a（插入）；对照基线 = validation panels 的同 job Beutler 值
- 判读：overlap/spread/ddG 对比；顺带验证 grompp 对 sc_alpha=0 的接受性

### E3 Sigmoidal λ 分布可选参数（零成本，代码）✅ 已实现 2026-08-14
- `ProtocolConfig.lambda_distribution: linear|sigmoidal`（默认 linear 不变）；tanh 端点加密 warp（steepness=2.5），作用于方向感知调度的 ramp 进度
- 单测验证：端点步长更小、单调有界、端点精确 0/1

### E4 Patel 协议对齐实验（最高杠杆，GPU ~1-2 天）🟢 已启动 2026-08-17
- 范围修正：Patel 3HFM 8 突变中仅 **3 个电荷守恒**（y20f −0.48 / w98f +3.25 / y50l +4.39，覆盖完整动态范围）；电荷变化的 5 个（r21a/d101k/d32n/n31e/k97d）属 V2.1 范围，本轮不测
- 实现：`run_3hfm_patel_align_20260817.sh`——标准 job 跑到 equilibrate → 每腿 40ns 自由 MD（弃前 10ns）→ 12.5/25/37.5ns 三帧注入 rep01/02/03 → 标准 FEP（方向感知调度 + 自适应 λ）
- 判读预设：3HFM 切片 R/MAE 是否显著逼近 Patel R²=0.81；逼近 → 端点长平衡+多帧 seed 升级默认协议；不动 → 采样假设排除，定位力场/模型层

### E5（候选，V2.1 方向）
- 质子化多边 + pKa 加权合并（FEP+ Groups 式；我们的 GLH/HIP 注入机制已验证，只欠热力学权重合并）
- co-alchemical 水↔离子（电荷变化突变，Clark 2019/feflow 路径）
- charmm36-mut 双力场 consensus（Gapsys 验证过的系统性增益）
- **V2.1a 详细技术方案已完成：`docs/v2_1a_dssb_implementation_cn.md`**（DSSB 单腿架构、反向杂化、Patel 5 突变验证目标 MAE≤1.5）

### E6 双力场 consensus（已立项待启动，2026-08-19 记录）

**策略**：同 job 双臂跑 amber99sb-star-ildn-mut + charmm36m-mut，ddG_consensus = mean（两力场系统误差不相关，平均抵消）。

**可行性已验证（探针：1DVF y49a，charmm36m-mut）**：
- vendored pmx 自带 charmm36m-mut.ff（651 morph 条目 + 水模型 + 离子 + cmap 齐全）
- mutate 全链路通过；Y2A 杂化 24 原子、电荷 0/0 整数闭合
- equilibrate 正常推进

**实施拆分（待启动）**：
1. 双臂 runner（protocol 级 force_field 切换 + ddG 合并层）
2. 小规模验证（1DVF/1MLC 5-8 job 双臂 vs 实验）
3. 验收集 6 靶点双臂 → 双力场官方指标

**注意点**：charmm 质子化命名不同（HSD/HSE/HSP、GLUP/ASPP）——charmm 臂初期用默认质子化；成本 2×。
**预期收益**：MAE 2.3 → ~1.8（Gapsys 2016/2020 实证路径）；可能部分缓解强效应动态范围膨胀。

## Phase D2：REST2 端点系综试点（2026-08-14 启动）

**与 RID 试点构成同源三引擎对照**（同一评估问题：端点系综能否把 y102a/y49a 的 ddG 均值拉向实验 4.79/1.90）：
- 自由 MD 5ns：已完成（均值 13.6/−6.9，仍偏）
- RID（ML-CV 全域探索）：basin-FEP 排队中
- **REST2（内部实现，`run_abag_rest2.sh`，seed-residues 直接锁定突变位点热区——比 RID 更直接的机制对口）**：脚本 `run_rest2_pilot_1dvf_20260813.sh`，seed D:102/A:49，3 replicas（1.0/0.94/0.88），rep00 无偏轨迹取 3 远隔帧 → FEP

**踩坑**：① normalize_structure_to_pdb.py 需要 BioPython（装了 tmp/rest2_venv）；② `-s rest2md` 依赖 prepare 产物，需先 `-s prepare`。

## Phase C：分析与扩展（与 B 并行）

- **C1** 腿级误差归因（complex vs apo）——用 32 个新 job 的 bar_summary
- **C2** ~~扩展样本外靶点准备~~ → 已完成（2026-08-11）：1AK4（15 job，CycA/HIV 衣壳）、1KTZ（15 job，TGF-β 受体/配体）、3K2M（5 job，Monobody/SH2）。结构全部干净（0 缺失 0 clash）、编号核对通过、电荷守恒过滤完成、PROPKA 审计显示三个靶点均无决定性质子化位点（按 \|pKa−pH\|>2 规则无需变体）——干净的前瞻验证集。案例定义在 `examples/real_cases/{1ak4,1ktz,3k2m}/`，runner `benchmarks/ab_bind/run_outofsample_expansion_20260811.sh` 已排队（等 Phase B 完成自动接力，35 jobs）。注意 3K2M 有 4.0×3 重复值（将被 censored 规则捕获）。
- **C3** PROPKA 集成设计（规则：|pKa−pH|>2 才应用变体；纯态边界教训见 §4.10④）

## 验收标准

1. A 阶段全部测试通过（含新增自适应 λ / 分步调度 / blocking QC 的测试）
2. B1 的 6 个 job overlap_min 全部 ≥0.2 且 MAE 下降
3. B2 overlap ≥0.2（证明分步调度修好了电荷变化转化）
4. B3 官方视图用化学完整数据重建，README/文档数字更新

## Phase H：Y/W/F 大删除的水合与端点系综专项（2026-08-24，下一轮优先方向）

### H1 原理与问题定义

Y/W/F 大侧链删除会同时改变界面空腔、局部水占据和端点构象；固定起始构象的
单轨迹 FEP 往往把未充分采样的脱水/重排自由能误记成突变本身的结合能，表现为
强效应幅度系统性膨胀。已有 RID/REST2 结果说明 basin 多样性和水合重排是真实
误差源，因此继续单纯增加 lambda 窗口或更换常规 soft-core 不是主路线。

### H2 首批靶点与突变

- 1BJ1/1CZ8：G88A、G92A、Q89A、I80A
- 3K2M：W80A、M88A、Y36A
- 1DVF：全部 Y/W 大删除案例

3K2M 的重复检测限值仍保留作诊断，但不把 censored 点当作精确连续回归真值。
1BJ1/1CZ8 作为同一 VEGF assay cluster 报告，避免把相关面板当作独立靶点。

### H3 实验协议

1. WT 与 MUT 端点分别进行至少 3 个独立 basin/构象种子；
2. 先用 RID 或 REST2 生成端点构象，再执行标准 FEP；
3. 记录界面水占据、空腔体积和 basin 间 ddG spread；
4. 对代表性 Y/W 删除增加 GCMC、ghost-water 或 OpenMM grand 预平衡对照；
5. 高柔性环区保留 dt=0.002 成本路径，但 sample 失败自动回退 dt=0.001。

### H4 预注册验收门槛

- 单点绝对误差相对当前基线下降至少 50%；
- 端点 ensemble 的 basin 均值向实验移动，且 basin spread 不超过 2 kcal/mol；
- 未使用 OOS 实验值调参；
- 正式 OOS 指标同时满足 Pearson R、Spearman、MAE 和 censored-aware 指标要求。

该专项是下一轮的 P1 方向；力场切换仅作为对照，不作为主优化路线。
