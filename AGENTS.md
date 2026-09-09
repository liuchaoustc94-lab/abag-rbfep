# AGENTS.md — abag-rbfep 项目工作约定

## 任务提交
- **不用 SLURM**（2026-09-08 用户决策）：任务直接 nohup/setsid 提交，GPU 用 ABAG_RBFE_VISIBLE_GPUS / CUDA_VISIBLE_DEVICES 显式指定。SLURM 虽然在体系里装了，但排队语义对我们的小批量交互式探索是负担。
- GPU3 长期被他人占用，项目用 GPU 0/1/2。

## 运维纪律（血泪史，全部有 ISSUE 记录）
- wipe 重跑三口径：① 删 lambda 目录必须连删 `stages/build_legs.json`；② 只删 sample.json 可幂等续跑（窗口产物保留时）；③ 强制重采样要连窗口产物（dhdl.xvg/md.gro/mutant 相关）一起删，否则恢复逻辑会拿旧数据重报
- `--runs-root` 必须绝对路径
- prune 前必须先再生报告
- grompp 在 GPU 污染环境用 `CUDA_VISIBLE_DEVICES=""` 强制 CPU（CUDA #700 连带事故）

## 评测口径
- 主指标：Pearson R / Spearman（排序能力）> MAE；分靶点+靶点中位数为主口径，pooled 为辅
- 验收集（2026-08-31 版）：剔除 2JEL/3NGB/3K2M/1JRH/1NMB/1MLC（数据质量原因，见 docs/validation_targets.md）
- 校正/过滤口径只作诊断，不进官方指标
