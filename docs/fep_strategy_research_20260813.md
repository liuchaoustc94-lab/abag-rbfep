# 开源 FEP 实现策略调研与借鉴清单（2026-08-13）

调研范围：pmx/de Groot lab、OpenFE/feflow、QresFEP-2（Åqvist）、BFEE/BAT.py/CHARMM-GUI、FEP+ 公开协议（Clark 2019 / Sampson 2024）、Patel 2021 JCTC（3HFM）、LiveCoMS 最佳实践、Hui & de Groot 2026 GCMC/REST2 GROMACS 引擎。
目的：为 abag-rbfep 已定位的四层误差（化学/采样调度/系综/力场水合）寻找经过验证的解法。

## 一、与我们误差结构精确对应的可借鉴策略

### 层 1：化学/拓扑（已修 ISSUE-001，防线加固）
| 策略 | 来源 | 可行性 |
|---|---|---|
| **naked-charge 自动检查**：任何 λ 点不允许"有电荷无 LJ"的原子 | feflow lambda_protocol.py | ✅ 直接实现为我们的 λ 调度 QC（方向感知调度的正式校验器） |
| **feflow "quarters" 分段**：关旧 ele→关旧 sterics→开新 sterics→开新 ele | feflow | ✅ 可作我们方向感知调度的对照基线 |

### 层 2：λ 调度/采样协议（部分已做，可补强）
| 策略 | 来源 | 可行性 |
|---|---|---|
| **按突变类型分档窗口数**：中性 12 / proline 16 / 电荷变化 24 | FEP+（Sampson 2024, Sergeeva 2023） | ✅ 扩展我们已有的自适应 λ（加 proline/电荷分档） |
| **Sigmoidal λ 分布**（端点加密）：QresFEP-2 协议矩阵系统优化胜出 | Koenekoop 2025 | ✅ GROMACS 任意 λ 向量即可 |
| **Gapsys linearized soft-core**（社区已从 Beutler 迁移）：`sc-alpha=0, sc-power=1, sc-sigma=0.3, sc-coul=yes`；我们当前 Beutler α=0.3/σ=0.25 | de Groot / feflow 默认 / Schmidt 2025 | ✅ 纯 mdp 改动，值得 A/B |
| **transition 时长按体系标定**：先用最极端 ±突变做收敛扫描再推广（1BRS 5ns vs 3HFM 8ns） | Patel 2021 | ✅ 加一个 pilot 收敛测试阶段即可 |

### 层 3：端点系综（我们正在做的方向，文献给出具体参数）
| 策略 | 来源 | 可行性 |
|---|---|---|
| **端点系综标准参数**：WT/MUT 各 ≥20 ns 平衡，~300 ps 间隔抽 ~100 帧 | de Groot 标准协议 / Patel 2021 | ✅ 我们的 basin seeding 即此方向；参数可直接采用 |
| **REST2/lambda-hopping**（升温区只含突变残基） | FEP+ 生产协议 | ⚠️ 需 PLUMED 补丁版；或由 RID 覆盖 |
| **GCMC 水采样处理界面水重排**——Y/W 删除空腔水合问题的主流解法 | FEP+ GCMC / Zhang 2023 / Hui & de Groot 2026（GROMACS 开源 GCMC/water-swap/REST2 引擎，与 pmx 同实验室） | ⚠️ **Hui & de Groot 2026 引擎是最可行路径**（pmx 同源） |

### 层 4：系统性误差（力场/质子化/电荷）
| 策略 | 来源 | 可行性 |
|---|---|---|
| **质子化多边 + pKa 加权合并**：D/E/H/K 突变跑交替质子化状态边，按实验 pH 热力学加权 | FEP+ Groups / Sampson 2024（证明实质降低系统误差） | ⚠️ 中等——我们已有 propka_venv + GLH/HIP 注入机制（已探针验证），加权重合并即可 |
| **双力场 consensus**（amber99sb*-ildn-mut + charmm36-mut 平均）| Gapsys 2016/2020（系统性降误差） | ⚠️ 中等——pmx 官方支持 charmm36-mut，需 vendor 力场 |
| **电荷变化：co-alchemical 水↔离子** 维持净电荷中性 | Clark 2019 / feflow explicit_charge_correction（\|Δq\|=1, PME） | ⚠️ V2.1 主路径 |
| **电荷变化备选：double-system/single-box**（bound+unbound 同盒 30Å 分隔，天然电荷中性） | Patel 2021（pmx 原生） | ⚠️ 实现成本中等，盒大 |

### 分析层（零成本）
- 多 replica（≥3-10）报 SEM；forward/reverse hysteresis 图作为 QC（LiveCoMS 底线）
- MBAR 不确定度用 bootstrap(1000)（OpenFE 实践：比解析误差更可靠）
- 双向 work 分布 + CGI 交点（pmx analyze_crooks.py 现成）
- 报告双栏 raw + calibrated，τ/ρ 为主排名指标

## 二、精度锚点（设定验收标准用）

| 来源 | 体系 | 指标 |
|---|---|---|
| Patel 2021 (pmx, NEQ) | **3HFM 抗体-抗原** | **R²=0.81**，非电荷突变 ±0.5-1 kcal/mol |
| Sampson 2024 (FEP+) | PPI 大规模 186 例 | RMSE 1.03, R² 0.4；三分类 balanced acc 0.69 |
| Clark 2019 (FEP+) | 162 电荷变化 PPI 突变 | RMSE 1.41（非埋藏 1.23 / 埋藏 1.79） |
| QresFEP-2 2025 | 583 折叠 ddG | MAE 1.25（优于 FEP+ 1.38） |
| 实验自身噪声 | — | 0.4-0.9 kcal/mol（RMSE~1 已接近天花板） |

**关键参照**：Patel 2021 的 3HFM R²=0.81 用的就是 pmx + Amber99SB*ILDN——与我们同工具链；其数据在我们 `benchmarks/patel_2021_3hfm/` 已有对照集。我们 3HFM 目前全失败，Patel 的协议（40ns 平衡 + 100 帧 + NEQ 双向）是直接的改进模板。

## 三、优先级建议（性价比排序）

**零成本立即做**：
1. naked-charge λ 校验器（加固方向感知调度）
2. Gapsys soft-core A/B（mdp 改动）
3. Sigmoidal λ 分布（可选参数）
4. 精度锚点写入验收标准

**中等投入高回报**：
5. 质子化多边 + pKa 加权（机制已通，加权重合并）
6. Patel 2021 协议对齐测试：用我们的数据复算其 3HFM 协议的关键差异（平衡时长/帧数/NEQ）
7. Hui & de Groot 2026 GCMC/REST2 引擎评估（Y/W 删除水合问题）

**大投入（V2.1 时代）**：
8. co-alchemical 水↔离子（电荷变化突变）
9. charmm36-mut 双力场 consensus

---

## 四、第二轮靶向调研增补（2026-08-19，针对已定位问题）

三个主题（空腔水合/动态范围、精度前沿、质子化与柔性环区）的调研结论。

### A. 对我们具体问题的直接答案

**1. Y/W/F 大删除的动态范围膨胀**
- Sampson 2024 中有精确对应类别（"buried aromatic, large size change"），其**单参数经验校正模式可直接移植**（RMSE 1.35→1.03 的关键来源之一）
- GCMC 路线：GrandFEP/GCNCMC 全是 OpenMM 系——要么 prep-only 嫁接（OpenMM `grand` 只做预平衡），要么双引擎；GROMACS 内 alchemical ghost-water 耦合需自行开发（中期项）
- AMOEBA/Drude 极化 FEP：GROMACS 不支持，❌

**2. 弱效应排序差（R≈0.17）**
- 实验噪声地板 0.4-0.9 kcal/mol 下，弱效应子集的 R 理论上限仅 ~0.3-0.6——**我们离上限并不遥远**；提升路径 = 把计算 SE 压到 0.2 以下（4-5 repeats）+ 按动态范围校正报告指标，而非换新方法

**3. 质子化中间态灾难**
- **Sampson 2024 Groups 式多态边 + pKa 热力学加权是直接解药**（纯后处理+图扩展，我们已有 GLH/HIP 注入机制）；Duan 2020 证明单点 RT·ln10·(pKa−pH) 惩罚就能把 >5 kcal/mol 灾难点拉回 <1
- 不建议做真正的 constant-pH FEP（GROMACS 主线无支持，工程风险>>收益）

**4. 柔性环区爆炸**
- 我们的方向感知重叠调度 = 文献铁律（de Groot 组对 indel 的标准答案）✓
- **环区端点预建模 / MCFEP 式窗口播种**（Sampson 对 1DVF H3 的做法）：WT/MUT 环分别约束构象搜索后进 FEP；半 A 半 B 窗口播种在 7 对难例上 RMSE 0.62 vs 1.61——**与我们现有 λ 窗口架构兼容**（只是不同窗口用不同初始构象）

### B. 新增的高价值策略（第一轮未覆盖）

| 优先级 | 策略 | 收益 | 投入 |
|---|---|---|---|
| 1 | **Petrov 2024 四指南审计**：电中性盒 + >1nm 缓冲 + 0.15M 盐 + 路径保总电荷（我们基本已满足，盒缓冲 1.0nm 需微调） | 消除电荷伪影 | ~0 |
| 2 | **cinnabar/MLE cycle-closure 分析层**：多突变图逆方差加权 MLE，报 pairwise 指标（与 FEP+ 同口径） | 免费压缩噪声 | 1-2 天 |
| 3 | **误差纪律**：独立 repeat SD 替代 bootstrap/MBAR 单轨迹误差 + 三分类 balanced accuracy 指标 | 可信不确定性 | 1 天 |
| 4 | **Titratable-state 预筛/双态枚举加权**（Sampson Groups 等价物） | 文献单项最大收益（RMSE 2.33→1.48） | 中 |
| 5 | **ALS 式自适应 λ**（50ps pilot → overlap matrix → 按需插窗） | 省 ~40-50% GPU | 中（alchemlyb 可复用） |
| 6 | Rocklin 解析校正后处理（电荷突变边交叉验证） | 发表级严谨性 | 低-中 |
| 7 | **λ-HREX**（GROMACS 原生 `-replex`） | 难边收敛加速 | 低 |
| 8 | Sampson 式 outlier 五参数分类 + 带电异常值经验校正 | RMSE 1.35→1.03 | 低 |

### C. 精度天花板校准（重要认知）
- 实验 assay 间 pairwise RMSE ≈ **0.91 kcal/mol**（Ross 2023）；FEP+ ≈ 1.25
- **RMSE 1.0 已逼近实验可复现性极限**——我们 MAE 2.35 → 目标 1.3-1.5 是合理的物理上限，不是 0.5
- 超过该点后必须用分类指标（三分类 BA，Sampson 参照 0.69）而非继续压 RMSE

### D. 一句话结论
**从 R≈0.5 到 FEP+ 水平（RMSE~1.0）的收益排序：体系准备（电荷/盐/缓冲+质子化态）> 分析层（cycle closure MLE + 异常值分类校正 + 诚实误差棒）> 采样效率（ALS + λ-HREX）>> ML 力场。** 前两级零/低模拟成本。
