# V2.1 设计约束：charge-changing 与跨侧 double-point

这份说明不是在提前实现 `V2.1`，而是把它和当前已经落地的 `V1/V2`
严格分开，避免把高风险电荷变化路径混进现在的默认 workflow。

## 为什么必须单独拆出来

当前仓库已经稳定下来的边界是：

- `V1`: 单点、标准残基、优先 `charge-conserving`
- `V2`: 同侧 `double-point`、标准残基、优先 `charge-conserving`
- 正式分析源: `gmx bar`
- 科学核心: `pmx`

`charge-changing` 和 `cross-side double-point` 不能直接沿用这条默认路径，
原因不是“再多几个窗口”这么简单，而是热力学循环、体系中性、有限尺寸误差、
约束设计和收敛诊断都会一起变复杂。

对当前项目，`V2.1` 的首要目标不是追求更广覆盖，而是避免污染已经能跑通的
`V1/V2` 主链路。

## 调研资料对项目的直接帮助

这批资料最有价值的地方不在于再证明一次 `pmx + GROMACS` 可行，而在于把
`V2.1` 的边界说清楚了：

- `Patel et al. 2021` 说明蛋白-蛋白突变自由能在大复合物里是可做的，但
  代价和稳定性要求明显高于小体系。
- `Clark 2019` 说明电荷变化突变确实值得做，但不应该和普通
  `charge-conserving` 突变共用同一默认 protocol。
- `pmx` 文档中的 `doublebox` 路径给了工程入口，说明 `V2.1` 不需要推翻
  当前 `pmx` 底座，而是要增加一条独立 setup 分支。
- `Klimovich et al.` 的分析规范给 `V2.1` 的 QC 设计提供了下限：更严格的
  overlap、重复一致性、不确定度表达不能省。

结论很直接：这些资料对项目有帮助，但帮助主要体现在“约束设计”和“避免误开
范围”，而不是立刻把默认采样时间调到更激进的数值。

## `V2.1` 的拆分策略

`V2.1` 继续拆成两个子阶段，不一次性全开：

### `V2.1a`: charge-changing single / same-side double

先只开放：

- 单点电荷变化突变
- 同侧 `double-point` 中包含电荷变化的情况
- 标准氨基酸
- 仍然限定 antibody 或 antigen 单侧变化

先不开放：

- 跨侧 `double-point`
- 非标准残基
- glycan 参与的 alchemical 路径

### `V2.1b`: cross-side double-point

只有在 `V2.1a` 已经完成真实案例、QC 和 benchmark 归档之后，再开放：

- 一个位点在 antibody，另一个位点在 antigen 的双点突变

这是因为跨侧双突变不只是多一个位点，而是会把结果汇总、缓存复用、热力学
解释和 failure taxonomy 一起变复杂。

## 对当前 workflow 的影响点

`V2.1` 不应重写全项目，而是对现有 stage 做受控增量扩展。

### 1. `mutation validate`

需要新增但暂不默认放开的信息：

- `net_charge_delta`
- `contains_charge_change`
- `requires_doublebox`
- `cross_side_double_point`

现阶段代码里已经有：

- `allow_charge_changing`
- `allow_cross_side_double_point`

但还缺：

- 明确写入 manifest / report 的电荷变化元数据
- 为什么被拒绝的可读原因
- `V2.1a` 与 `V2.1b` 的不同 gate

### 2. `mutate`

当前 `mutate` 是：

- `pmx mutate`
- `pdb2gmx`
- `pmx gentop`

`V2.1` 需要在这里新增独立分支，而不是改写现有默认分支：

- 默认 `charge-conserving` 继续走当前路径
- 电荷变化路径进入 `doublebox` / DSSB 风格 setup

这意味着 `mutate` 产物不能再只假设“一个混合蛋白 + 一个普通溶剂盒”。

### 3. `build_legs`

这里会是 `V2.1` 的主要变化点：

- 需要支持 double-system/single-box 风格的腿构建
- 需要记录 bound/unbound 两个体系在同一盒中的构型关系
- 需要把电荷中和策略作为 protocol 的显式字段，而不是隐含在脚本里

如果这一层没有独立建模，后面的 `equilibrate/sample/bar` 会继续假设当前
`V1/V2` 的目录结构，从而把 `V2.1` 弄成“看起来能跑，实际上语义不对”。

### 4. `equilibrate`

这里至少要新增三类约束：

- 体系总净电荷检查
- doublebox 几何检查
- 与电荷变化路径一致的 restraint 元数据

调研资料里关于 identical restraints、盒内分离距离、离子强度和稳定性检查
都应作为 `equilibrate` 的 QC 输入，而不是只留在文档里。

### 5. `sample / bar / qc / report`

正式结果源仍然保持 `gmx bar`，这一点不变。

但 `V2.1` 的报表必须额外回答下面这些问题：

- 这个作业是否为 `charge-changing`
- 使用了哪种 setup 路径
- 每条腿的 overlap 是否足够
- 重复间 spread 是否在可接受范围内
- 结果应被视为定量还是偏定性参考

换句话说，`V2.1` 不是换分析器，而是要让现有 `BAR + QC + report`
知道自己正在处理更困难的物理路径。

## 不应直接照搬到默认协议的调研结论

下面这些结论有参考价值，但现在不应该直接塞进 `validation` 主链路：

- `40-100 ns` 抗体-抗原级长采样
- `3-5` 个 replica 作为默认值
- `REST2` / `GaMD` / `lambda-REMD` 作为默认策略
- 用 `DC-MBAR` 替换 `gmx bar`

它们更适合：

- `robust` / `deep rescue` 的上限参考
- 难例专项救援
- 后处理 sidecar 分析

而不是现阶段 held-out validation 的默认入口。

## 推荐的工程落地顺序

### 第一批：只补边界和元数据

1. 在 mutation validation 输出里加入 `net_charge_delta` 与
   `contains_charge_change`
2. 在 `job_spec` / `report` / `plan_jobs.csv` 里显式暴露这些字段
3. 将 `V2.1` 的拒绝原因标准化，而不是只返回通用校验失败

### 第二批：新增独立 setup 分支

1. 为 `build_legs` 引入 `doublebox` 路径开关
2. 把 `charge-conserving` 与 `charge-changing` 的目录/产物语义分开
3. 为 `equilibrate` 增加 doublebox 几何和体系中性检查

### 第三批：真实案例和 benchmark

1. 先跑一个真实 `charge-changing single-point`
2. 再跑一个真实 same-side `charge-changing double-point`
3. 两者都能稳定 through `report` 后，再考虑跨侧 `double-point`

## 与当前项目目标的关系

当前 held-out validation 的主要瓶颈仍然是：

- 已完成可配对作业数太少
- 真实大体系 throughput 受限
- 少数作业被旧的 mutate/equilibrate 问题卡住

所以这些调研资料**有帮助**，但它们现在对项目最直接的作用不是立刻提高
`R > 0.6`，而是：

- 防止我们把 `V2.1` 误做成 `V2` 的小补丁
- 给后续 charge-changing 路径一个清晰的工程入口
- 让当前 `V1/V2` 验证与未来 `V2.1` 扩展保持解耦

## 当前结论

可以明确说，这些调研资料是有帮助的，而且帮助是高价值的。但它们更像：

- `V1/V2` 的边界校验器
- `V2.1` 的设计输入
- `3HFM` 外部校准和 failure taxonomy 的 backlog 依据

而不是今天就该把默认 protocol 改成更重采样、更复杂增强采样的理由。
