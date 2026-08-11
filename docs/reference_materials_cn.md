# 参考材料速读（中文）

这份清单是给当前 `abag-rbfep` 路线用的，不是泛泛的自由能文献综述。它只服务于下面这条边界：

- `V1`: antibody 或 antigen 任一侧 `single-point` ddG
- `V2`: 同侧 `double-point` ddG
- `V2.1`: `charge-changing` 与跨侧 `double-point`
- 技术底座: `GROMACS + pmx`
- 当前正式结果源: `gmx bar`

如果你现在时间有限，先看“最短必读清单”；如果你要开始设计 `V2/V2.1`，直接跳到后面的“按版本阅读”。

## 最短必读清单

### 1. GROMACS 自由能主文档

- `Free energy calculations`
  - 作用: 先把 `lambda` 路径、自由能端态、DHDL 输出、积分/BAR 所在的位置看清楚。
  - 重点:
    - GROMACS 如何定义端态和路径
    - `dhdl`/Hamiltonian difference 从哪里来
    - 为什么大分子突变仍然是“端态 + 路径”问题，而不是单独的打分问题
- `gmx bar`
  - 作用: 当前仓库的正式分析源。
  - 重点:
    - BAR 的输入是什么
    - 相邻窗口如何汇总成总 `Delta G`
    - 结果里哪些字段值得进入 `qc/report`
- `Free energy interactions`
  - 作用: 以后调 `lambda` 组件、soft-core、bonded/non-bonded 插值时必须回头看这里。
  - 重点:
    - `lambda` 是向量，不同能量项可以独立推进
    - bonded 和 non-bonded 的插值行为不同

建议顺序:

1. `free-energy-calculations`
2. `gmx bar`
3. `free-energy-interactions`

## 2. pmx 科学核心

### Gapsys et al. 2015

- 题目: `pmx: Automated protein structure and topology generation for alchemical perturbations`
- DOI: `10.1002/jcc.23804`
- 作用:
  - 这是 `pmx` 的核心论文，也是你项目里 “mutation -> hybrid residue/topology -> FE setup” 这条链路最该引用的来源。
- 重点:
  - hybrid residue / topology 是怎么定义的
  - 为什么它适合做氨基酸突变的 alchemical setup
  - 它天生就不是只做单点突变的工具

### 官方文档和教程

- `pmx` 主页
- `tutorials/protein_mut`
- `examples/analysis`
- `scripts/gentop`
- `scripts/index`

先看什么:

1. `protein_mut`
2. `analysis`
3. `gentop`
4. `scripts/index`

重点:

- `protein_mut`: 理解突变定义、输入格式、hybrid 结构生成
- `analysis`: 看 `pmx` 如何把分析交给 `gmx bar`
- `gentop`: 理解 topology 完整化为什么必须保留在主链路
- `scripts/index`: 先记住 `doublebox` 的位置，给 `V2.1` 留钩子

## 3. 抗体抗原 benchmark 与问题定义

### AB-Bind

- 题目: `AB-Bind: Antibody binding mutational database for computational affinity predictions`
- DOI: `10.1002/pro.2829`
- 作用:
  - 这是当前项目最合适的外部主 benchmark 源。
- 重点:
  - 数据集规模和原始字段
  - 为什么真正能跑 RBFE 的样本一定比原始源集更少
  - 抗体抗原 `ddG` 不是普通 protein stability 问题

### AB-Bind 仓库

- 用途:
  - 对照原始表头、结构映射字段、mutant 标识和清洗过程。

## 4. 与当前工程最接近的蛋白-蛋白 alchemical 论文

### Patel et al. 2021

- 题目: `Implementing and Assessing an Alchemical Method for Calculating Protein-Protein Binding Free Energy`
- DOI: `10.1021/acs.jctc.0c01045`
- 作用:
  - 这是最接近你真实执行问题的一篇，不是泛理论论文。
- 重点:
  - protein-protein 突变自由能的执行难点
  - `pmx` hybrid setup 在大复合物里的现实成本
  - `charge-changing` 为什么要单独设计，不应和 `V2` 混着上

## 5. QC / 分析底座

### Klimovich, Shirts, Mobley 2015

- 题目: `Guidelines for the analysis of free energy calculations`
- DOI: `10.1007/s10822-015-9840-9`
- 作用:
  - 这是 `qc/report` 设计最实用的参考。
- 重点:
  - overlap 检查
  - 不确定度表达
  - 重复间一致性
  - 何时把一个结果判成 warning/fail

### Bennett 1976

- 题目: `Efficient estimation of free energy differences from Monte Carlo data`
- 作用:
  - BAR 的源头论文。
- 重点:
  - 为什么相邻状态 overlap 不够时，误差会直接变坏

### Shirts and Chodera 2008

- 题目: `Statistically optimal analysis of samples from multiple equilibrium states`
- DOI: `10.1063/1.2978177`
- 作用:
  - MBAR 的理论底座，先备着，不急着上生产主链路。

## 6. 你给的论文：DC-MBAR

### Jia, Ge, Mei 2021

- 题目: `Free energy change estimation: The Divide and Conquer MBAR method`
- DOI: `10.1002/jcc.26533`
- 作用:
  - 这篇最适合放在“未来分析 sidecar”位置，而不是现在替代 `gmx bar`。
- 先看哪几点:
  - 它如何定义相邻状态
  - overlap threshold 怎么设
  - local MBAR 如何拆分和重构总 free-energy profile
- 为什么现在先不替换 BAR:
  - 当前仓库正式口径已经固定在 `gmx bar`
  - overlap threshold 有明显系统依赖性
  - 论文主场景不是 antibody-antigen 的 pmx protein mutation RBFE

适合未来放到哪里:

- dense `lambda` 协议的事后重分析
- 低 overlap 窗口诊断
- `double_point` 保守协议下的分析 sidecar

## 7. V2 的外部多点参照

### mmCSM-AB

- 题目: `mmCSM-AB: guiding rational antibody engineering through multiple point mutations`
- DOI: `10.1093/nar/gkaa406`
- 作用:
  - 不是拿来复用算法。
  - 是拿来定义“抗体抗原 multiple-point mutation 是独立 benchmark 问题”的外部参照。
- 重点:
  - 多点突变的任务定义
  - 汇报格式
  - blind test / benchmark 的组织方式

## 8. 工作流层借鉴

### biobb_pmx

- 仓库: `bioexcel/biobb_pmx`
- 作用:
  - 不适合变成当前项目的运行时依赖。
  - 很适合借鉴它怎么把 `pmx` 封装成可复现、可组合、可追踪的 building blocks。

推荐只借鉴这些方面:

- idempotent stage 边界
- 参数显式化
- 产物路径和元数据组织

不建议直接照搬这些方面:

- 运行时依赖体系
- 项目整体对象模型

## 9. 这批补充调研对当前仓库的直接价值

这批资料是有帮助的，但不是每一条都该立刻变成默认协议。对当前
`abag-rbfep`，更合适的用法是分三层吸收。

### 已经被当前仓库吸收的结论

- `pmx + GROMACS + gmx bar` 作为主链路，这一点已经和当前仓库边界一致。
- `NVT -> NPT` 预平衡不是可选项，当前实现已经保留了 `nvt.mdp` 与
  `npt.mdp` 两段。
- 当前主链默认仍使用 `C-rescale`，这一点已经落实在生成的 `mdp` 里；
  但 Patel 2021 的 `3HFM` 正文和补充材料并不支持把它直接当成该论文的
  结论。该文正文描述的是 `Parrinello-Rahman`，补充 `mdp` 还出现了
  `Berendsen/Parrinello-Rahman` 混用，所以它更适合作为外部 regression
  参照，而不是直接替换主链默认值。
- `AB-Bind` 作为主 benchmark 源、`mmCSM-AB` 作为 `V2` 多点问题的外部
  framing、`DC-MBAR` 作为未来 sidecar 分析输入，这些判断已经进入仓库
  文档。

### 值得进入下一阶段 backlog 的结论

- `Patel et al. 2021` 里的 `3HFM` 应该补成一个外部校准案例。它不是替代
  `AB-Bind`，而是给协议级回归一个更直接的抗体-抗原 alchemical 参照。
- `Sampson et al. 2024` 的异常值分类，适合转成以后的 validation failure
  taxonomy，用来区分采样不足、结构问题、带电变化、局部重排等失败类型。
- `Clark 2017/2019` 适合放进 `V2.1` 设计文档，用来约束
  `charge-changing` 和更激进采样方案的边界，而不是提前混入 `V1/V2`
  默认流程。
- `biobb_pmx` 仍然适合只借鉴工作流分层和元数据组织，不要引成运行时依赖。

### 现在不应直接上默认协议的结论

- 文献里的 `40-100 ns`、`3-5 replicas`、`100` 个快照、长时间非平衡切换，
  可以作为 `robust/deep` 级协议上限参考，但当前仓库还不能在 validation
  主链路里直接把它们当成默认值。
- `REST2`、`GaMD`、`lambda-REMD` 更适合作为困难体系救援路线，不该在
  `V1` 主协议里默认启用。
- `DC-MBAR` 现在仍然更适合做 `gmx bar` 旁路重分析，而不是替换正式结果源。

### 当前最该转成工程动作的三件事

1. 增加一个 `3HFM` 外部校准运行与结果归档。
2. 在 validation 报表里加入更明确的 failure taxonomy。
3. `V2.1` 的 `charge-changing` 设计约束已经单独整理为
   [v2_1_charge_design_cn.md](/mnt/data/liuchao/abag-rbfep/docs/v2_1_charge_design_cn.md)，
   后续实现应以该文档为边界，避免和 `V2` same-side `double-point` 混线。

## 按版本阅读

### V1: single-point

只看下面这些就够了:

1. GROMACS `free-energy-calculations`
2. `gmx bar`
3. pmx 论文
4. `pmx/tutorials/protein_mut`
5. AB-Bind
6. Klimovich et al.

### V2: same-side double-point

在 V1 基础上加：

1. Patel et al. 2021
2. mmCSM-AB
3. `pmx/scripts/index` 里和多突变、`doublebox` 相关的脚本入口
4. DC-MBAR 作为后处理设计输入

### V2.1: charge-changing / 跨侧 double-point

在 V2 基础上再加：

1. `pmx` 的 `doublebox`
2. GROMACS `free-energy-interactions`
3. Klimovich et al. 重新看 overlap 和 uncertainty 那一部分

## 仓库内可以直接读的本地材料

这些都已经在当前仓库里，读起来最快：

- [reference_materials.md](/mnt/data/liuchao/abag-rbfep/docs/reference_materials.md)
- [architecture.md](/mnt/data/liuchao/abag-rbfep/docs/architecture.md)
- [upstream_pmx.md](/mnt/data/liuchao/abag-rbfep/docs/upstream_pmx.md)
- [protein_mut.rst](/mnt/data/liuchao/abag-rbfep/vendor/pmx/docs/tutorials/protein_mut.rst)
- [analysis.rst](/mnt/data/liuchao/abag-rbfep/vendor/pmx/docs/examples/analysis.rst)
- [gentop.rst](/mnt/data/liuchao/abag-rbfep/vendor/pmx/docs/scripts/gentop.rst)
- [index.rst](/mnt/data/liuchao/abag-rbfep/vendor/pmx/docs/scripts/index.rst)
- [large_scale_scans.rst](/mnt/data/liuchao/abag-rbfep/vendor/pmx/docs/tutorials/large_scale_scans.rst)

## 如果你只有半天

按这个顺序读：

1. GROMACS `free-energy-calculations`
2. `gmx bar`
3. pmx 论文
4. `pmx/tutorials/protein_mut`
5. AB-Bind
6. Patel et al. 2021

这六项读完，已经足够支撑 `V1` 和 `V2` 的大部分工程判断。

## 如果你下一步就要改代码

直接对应到当前项目：

- 改 `mutation/topology`：先看 pmx 论文 + `protein_mut` + `gentop`
- 改 `bar/qc/report`：先看 `gmx bar` + Klimovich + Bennett
- 设计 `double_point`：先看 Patel + mmCSM-AB
- 想评估 MBAR sidecar：最后再看 DC-MBAR + Shirts-Chodera
