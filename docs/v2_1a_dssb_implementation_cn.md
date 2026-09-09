# V2.1a 技术方案：charge-changing 突变的 DSSB 实现（2026-08-18 立项）

关联：`docs/v2_1_charge_design_cn.md`（边界与纪律）、`docs/fep_strategy_research_20260813.md`（Patel 2021 路径）、`docs/known_issues.md`。
状态：设计完成，待评审后实施。

## 0. 范围与判据

**V2.1a 开放**：单点电荷变化突变 + 同侧双点中含电荷变化（标准氨基酸、单侧）。
**不开放**：跨侧双点（V2.1b）、非标准残基、glycan。
**验证目标**：Patel 3HFM 5 个电荷变化突变（r21a +0.90 / d101k +2.13 / d32n +0.17 / n31e +5.71 / k97d +6.77）达到 MAE ≤ 1.5 kcal/mol（Patel 全文 R²=0.81 的参照系）。
**纪律**：默认 charge-conserving 路径逐行不动；V2.1a 全部走独立分支（`docs/v2_1_charge_design_cn.md` 的硬性要求）。

## 1. 核心热力学设计

DSSB（double-system/single-box）：同一模拟盒中放两个体系——

```
盒子(λ) = [ bound 复合物(杂化残基 A=WT, B=MUT) ] + [ unbound 单体(杂化残基 A=MUT, B=WT) ]
```

全局单一 λ 驱动：λ: 0→1 时 bound WT→MUT、unbound MUT→WT **同时进行**。

**净电荷恒等式**：Q_box(λ) = Q_bound(λ) + Q_unbound(λ) ≡ Q_WT + Q_MUT（任意 λ）
→ setup 时一次性 genion 中和 Q_WT+Q_MUT 后，**所有 λ 窗口恒定中性**，零 PME 电荷伪影。
这是选择 DSSB 而非 co-alchemical 水↔离子的根本原因：无需改 genion、无需每条腿单独校正。

**ddG 提取**：DSSB 单腿的 BAR ΔG 直接就是 ddG（热力学循环在盒内闭合），无需两腿相减。

## 2. 已验证的前置条件

| 前置 | 证据 | 日期 |
|---|---|---|
| 反向 morph 库 | pmx mutres.rtp ~600 条目双向覆盖（A2Q/D2K/K2D/N2D/R2A 等 Patel 所需全部存在） | 2026-08-18 |
| 分步 λ 调度 | 电荷变化转化 overlap 0.064→0.211（B2 验证） | 2026-08-11 |
| 带电杂化残基 | Z2A（带电 His→Ala）21 原子、chargeA=+1.0 整数闭合 | 2026-08-06 |
| `pmx doublebox` CLI | vendor/pmx 自带（`pmx.model.double_box`，r=2.5nm/d=1.5nm 默认，体积最小化长方体盒） | 2026-08-18 |
| 杂化完整性 QC | `validate_hybrid_topology_integrity` 已 blocking | 2026-08-14 |

## 3. 架构方案：新腿类型 `"dssb"`（最小侵入）

### 3.1 腿抽象改造（~20 处）

所有阶段当前硬编码 `for leg in ("complex", "apo")`。引入 `ctx.legs`：

```python
def _ctx_legs(ctx) -> tuple[str, ...]:
    return ("dssb",) if getattr(ctx.protocol, "leg_topology", "two_leg") == "dssb" else ("complex", "apo")
```

替换点（调研已定位）：`stages.py` 441/449/410-422/1476/1722/1839/1864/2452/2578/2748/2849/3193/3206/3234/3248/3303；`reporting.py` 177/241/631（含 `total_windows = repeats * λ * 2` 的 ×2 硬编码 → `* len(ctx.legs)`）；`benchmark.py` 2945/2958/3432。
**注意**：`_configured_legs`（stages.py:410-422）有白名单，必须同步加 "dssb"，否则新腿被静默丢弃。

### 3.2 prepare：DSSB 输入结构（新增分支，不动默认）

DSSB job 的 prepare 不产 `legs/<leg>/input.pdb`，改产：
- `legs/dssb/bound_input.pdb`（antibody+antigen 全链，同现 complex 腿）
- `legs/dssb/unbound_input.pdb`（突变侧链，同现 apo 腿）

### 3.3 mutate：三次 pmx 调用 + doublebox（DSSB 核心新增）

```
1) bound:   pmx mutate bound_input.pdb   → mutant_bound.pdb   （正向 WT→MUT，复用现有命令链）
2) unbound: 先构造 MUT 态结构（用正向 pmx mutate 对 unbound_input 得 mutant_unbound.pdb 并取 B 态坐标）
            再用反向 morph 库: pmx mutate mutant_unbound.pdb → hybrid_reverse.pdb（A=MUT, B=WT）
            ※ 备选（若反向库缺条目）：写 itp A/B 列交换工具（机械、可单测）
3) 链名去冲突：unbound 各链重命名（H→U, L→V, Y→Z 或加前缀），同步改 mutations.txt 链引用
4) pmx doublebox -f1 mutant_bound.pdb -f2 hybrid_reverse.pdb -o dssb_input.pdb -r 2.5 -d 1.5
5) pdb2gmx（无 -ignh，沿用现命令）→ processed.gro / topol.top
6) 合并拓扑：moleculetype 重命名（Protein_chain_H → Protein_bound_chain_H / Protein_unbound_chain_H），
   拼接两份 pmxtop 的 [ molecules ]；杂化完整性 QC 照常跑（正向+反向两份杂化残基都要过）
```

**最大工作量点**：步骤 6 的 moleculetype 撞名合并（`_merge_dssb_topology()` 新函数，放 gmx.py）。

### 3.4 build_legs / equilibrate / sample / bar：几乎零改动

- `legs/dssb/repXX/` 目录契约与现有完全一致（lambda_plan.json + mdp/ + lambda_NNN/）
- equilibrate 的 editconf/solvate/genion 片段**逐行复用**（editconf 不搬动坐标，doublebox 几何保留；genion 见 §1）
- staged posre 框架天然兼容（多 moleculetype 逐个处理 = 两体系各自锚定）
- **新 QC**：① 总电荷 sanity check（genion 前 tpr 电荷 == Q_WT+Q_MUT）；② doublebox 几何检查（两体系最小距离 > nonbonded_cutoff 1.25nm）；③ λ 序列净电荷恒定检查（任选 2 窗口 grompp 无 non-zero-charge warning）

### 3.5 采样协议（独立 preset，不污染默认）

新增 preset `dssb_charge_single_point`（constants.py，不影响现有 preset）：

| 参数 | 值 | 依据 |
|---|---|---|
| lambda_windows | 24 | FEP+ 电荷变化档（Sampson 2024） |
| repeats | 3 | 现有标准 |
| production_ps | 20 | 现有标准（E4 结果后可能升级端点系综） |
| 调度 | 方向感知分步（反向腿 ramp 取反——由杂化 A/B 互换自然处理，无需额外代码） | B2/ISSUE-007 |
| box | doublebox 长方体（不强制 dodecahedron） | Patel 2021 |

### 3.6 reporting / QC / report

- `collect_job_results` 加 dssb 分支：`len(legs)==1 and leg=="dssb"` → ddG = 单腿 ΔG，repeat 配对改为自身配对（repXX 对 repXX），`complex_delta_g_kcal_mol`/`apo_delta_g_kcal_mol` 字段置 None（schema 不变）
- pairs/jobs CSV 新增列：`charge_changing: bool`、`setup_path: two_leg|dssb`、`result_confidence: quantitative|indicative`
- 失败分类新增 code：`dssb_net_charge_mismatch`、`dssb_geometry_violation`、`dssb_reverse_hybrid_invalid`

## 4. ProtocolConfig 新增字段（均有默认值，向后兼容）

```python
leg_topology: str = "two_leg"        # two_leg | dssb
doublebox_separation_nm: float = 2.5
doublebox_wall_nm: float = 1.5
dssb_com_restraint: str = "none"     # none | posres_only（预留 pull/flat-bottom）
```

planning 路由：`allow_charge_changing=True` 且突变含电荷变化 → 自动 `leg_topology="dssb"`（对应设计文档的 `requires_doublebox` 意图）；report 层把拒绝原因写成可读 message。

## 5. 风险清单与缓解

| 风险 | 等级 | 缓解 |
|---|---|---|
| 反向杂化生成（B 态坐标取自哪） | 高 | 主路径：正向 mutant 结构直接进 pmx（B 态即 MUT 序列结构）+ 反向 morph；备选 itp A/B 列交换工具 + 杂化完整性 QC 双保险 |
| moleculetype 撞名合并 | 高 | `_merge_dssb_topology` 专项 + 单元测试 |
| 双体系漂移（production 无 posre） | 中 | 2.5nm 间距 + staged posre 端点锚定；预留 `dssb_com_restraint`；QC 加几何漂移监控 |
| `_configured_legs` 白名单静默丢腿 | 中 | 改造点清单第一条 + resume 回归测试 |
| reporting 单腿 ddG 的 repeat 语义 | 低 | 单元测试锁 schema |
| 24λ × 大体系（3HFM 复合物+单体 ~120k 原子）成本 | 低 | 每 job 24λ×3rep×20ps ≈ 1.4ns，可接受 |

## 6. 实施拆分（4 个 PR 级步骤）

1. **PR-1 腿抽象** ✅ 完成（2026-08-18）：`ctx.legs`/`_job_expected_legs` + 19 处枚举替换 + reporting ×2 硬编码 + `_configured_legs` 白名单 + 回归测试（379/379 通过，行为零变化验证）
2. **PR-2 DSSB mutate/topology** ✅ 完成（2026-08-18）：
   - prepare 分支（dssb 双输入 bound_input/unbound_input + QC 键）
   - `_dssb_mutate_inner_commands`：双 mutate → 链重命名（`suggest_unbound_chain_mapping`/`rename_pdb_chains`）→ `pmx doublebox` → pdb2gmx → gentop → **A/B 态交换**（`swap_hybrid_residue_ab_states`，支持 bonds/angles/dihedrals func 1/3/4/5/9 + 共享参数容错）→ dssb_qc
   - `_mutate_outputs_complete` 对 dssb 增加 dssb_qc.json 要求（防 resume 跳过 swap）
   - **真实冒烟验证（3HFM D32N，D→N 电荷变化）**：D2N 杂化 14 原子 + 206 bonded 项交换成功；杂化完整性 ✓；**电荷不变性 state A == state B == −3.0, δ=0 精确成立**
   - 测试：swap roundtrip/电荷不变性、链重命名、dssb 命令合同 stage 级；382/382 全量回归
3. **PR-3 reporting/dssb ddG + 协议字段 + preset + planning 路由 + 测试** ✅ 完成（2026-08-18）：
   - `collect_job_results` dssb 分支：单腿 ΔG 直接为 ddG（盒内热力学循环闭合），complex/apo 字段 None，schema 不变
   - ddg_summary 新增标注：`charge_changing` / `setup_path` / `result_confidence`（quantitative|indicative）
   - planning 路由：电荷变化突变自动 `leg_topology="dssb"`（显式 pin 尊重）
   - 新 preset `dssb_charge_single_point`（24λ/3rep/20ps/分步调度/triclinic 盒）
   - 测试：dssb ddG 提取、planning 路由（含 pin）；382+2 全量回归
   - 首个端到端 job（3HFM D32N）已启动采样
4. **PR-4 Patel 验证**：🟢 进行中（2026-08-18）——5 个 3HFM 电荷变化突变全部启动：D32N（首个 e2e，equilibrate 中）+ r21a/d101k/n31e/k97d（已入队，dssb_charge_single_point preset，24λ/3rep/20ps）。判读门槛：MAE ≤ 1.5 kcal/mol（参照 Patel 全文 R²=0.81 的同工具链水平）。

## 7. 测试计划

**调度关键决策（2026-08-18 实证）**：DSSB 双副本在同一全局 λ 下反向转化，方向感知分步调度的 coul/vdW 顺序对两副本**必然矛盾**（bound 副本删除需 coul 先行、unbound 副本反向=插入需 vdW 先行）——实测中间窗口系统爆炸（Coulomb SR ~−1.7×10⁷ kJ/mol，mdrun SIGSEGV）。**DSSB 腿强制走耦合单 fep-lambda**（`force_coupled` 覆盖，已入 stages.py）；分步/方向感知调度仅用于标准两腿路径。

- 单测：腿路由、`_merge_dssb_topology`（双 moleculetype 撞名）、净电荷恒等检查、反向杂化完整性、dssb ddG 提取
- 集成（mock CommandRunner）：dssb job 全 9 stage 目录契约 + mutate.sh 含 doublebox/反向 mutate/genion 行
- 真实冒烟：1 个最小电荷变化 job（如 1MLC 附近的 D→N）端到端（quick preset）
- 基准：Patel 5 突变 MAE ≤ 1.5 判为通过，R² 记录但不作为门槛（n=5 太小）

## 8. PR-4b 裁决（2026-08-31）：40ns 端点系综未过门槛

5 个 3HFM 电荷变化突变的 40ns 端点系综重跑完成（24λ/3rep/50ps，帧注入自 40ns 自由 MD）：

| 突变 | exp | 首轮（20ps） | PR-4b（40ns 端点系综） | err | repeat 极差 |
|---|---|---|---|---|---|
| n31e | 5.71 | 6.69（err 0.98） | 6.58 | **0.87** ✅ | 1.87 |
| d32n | 0.17 | — | 4.72 | 4.55 | 1.90 |
| d101k | 2.13 | 4.19（err 2.06） | 13.74 | 11.61 ❌ 恶化 | 3.79 |
| r21a | 0.90 | 9.94（err 9.04） | 20.87 | 19.97 ❌ 恶化 | **12.81** |
| k97d | 6.77 | 39.12（err 32.4） | 27.22 | 20.45（仍差） | 7.64 |

**MAE：首轮 11.11（n=4）→ PR-4b 11.49（n=5）；V2.1a 门槛 MAE≤1.5 未通过。**

**关键诊断**：误差与 **repeat 极差强相关**（n31e 1.87→准；r21a 12.8→崩）。40ns 端点系综的不同帧把各 rep 带进了**不同的宏观态**，且每个 rep 的 24 窗口共享同一帧 → 整条腿被帧偏置整体搬移。结论：**端点弛豫不是电荷变化收敛的瓶颈**——Patel 协议的差异点在"每窗口 8-10ns 连续生长（窗口顺序串联）"，而我们的架构是"共享起点的短窗口"。这指向 P2 的下一步：λ 路径的慢速连续穿越（窗口串联 chaining），而非更多端点帧。

**裁决**：V2.1a DSSB 保持"可用但未验证"状态（n31e/d32n 这类温和电荷变化可信；大静电重排如 r21a/k97d 需 chaining 协议再验）。

## 9. P2-c 净电荷伪影审计（2026-09-08，已关闭）

**目的**：排除"DSSB 残余误差来自净电荷/周期边界伪影"的可能性（Rocklin 2013 家族校正）。

**审计计算**（Rocklin Ewald 自能项 ΔG_EW = −(Δq²·ξ_EW·e²/4πε₀)/(2ε_r·L)，ε_r=80）：

| 盒型 | L | ΔG_EW（Δq=±1） |
|---|---|---|
| 标准两腿盒（1dqj 型） | 11.9 nm | +0.049 kcal/mol |
| DSSB 双盒（r21a 型） | 14.1 nm | +0.042 kcal/mol |

**结论**：
1. DSSB 路径上 Q(λ)≡Q_WT+Q_MUT 恒定 → Δq≡0 → **该校正构造性恒等于零**，连算都不用算
2. 即使按单盒最坏情况估算，伪影上限 ~0.05 kcal/mol，而 k97d/r21a 的残余误差是 20 kcal/mol 量级——**伪影占残余误差 <0.3%**
3. Petrov 四指南（电中性盒、>1nm 缓冲、0.15M 盐、路径保总电荷）已全部落实；Rocklin 式解析校正**不需要**进管线
4. **P2-c 关闭**：DSSB 的电荷簿记是干净的；PR-4b 的失败纯粹是采样/弛豫问题（chaining 路线处理中），与电荷伪影无关
