# Project Completion Status

Snapshot date: June 26, 2026.

This file can be regenerated with:
`python benchmarks/ab_bind/report_project_completion.py --root /mnt/data/liuchao/abag-rbfep`

This audit separates two questions that had been getting conflated in manual updates:
whether the accepted independent validation gate is already satisfied, and whether
the remaining AB-Bind and external 3HFM execution waves have fully drained.

## Independent Validation Gate

- summary file: `docs/validation_target_summary/validation_target_summary.json`
- generated at: `2026-06-25T06:40:58.525159+00:00`
- selected calibration model: `side_linear`
- accepted holdout view: `target_filtered`
- accepted excluded complexes: `1MLC, 1CZ8, 1BJ1`
- accepted calibrated `Pearson R = 0.6073390390160122`
- full unfiltered calibrated `Pearson R = 0.20005445478956882`
- held-out pair count: `80`
- accepted pair count: `42`
- accepted gate passed: `True`

## Real-Case Execution Gate

- required real-case checkpoints completed: `True`
- same-side double-point checkpoint completed: `True`

- 1VFB single-point quick validation:
  - file: `runs/real_cases/1vfb_y32f_quick/jobs/1vfb-antibody-b-y32f/results/ddg_summary.json`
  - generated at: `2026-06-05T04:21:40Z`
  - mutation: `B:Y32F@antibody`
  - mutation count / preset / side: `1` / `single_point` / `antibody`
  - `ddG = -3.5914742566415683 kcal/mol`
  - ddG BAR stderr: `5.6580069839823794 kcal/mol`
  - QC: `warning`
  - execution checkpoint passed: `True`
  - warnings:
    - complex:rep01 overlap score 0.091 below threshold 0.200
    - apo:rep01 overlap score 0.000 below threshold 0.200
- 1VFB same-side double-point quick validation:
  - file: `runs/real_cases/1vfb_y32f_v34i_quick/jobs/1vfb-antibody-b-y32f--b-v34i/results/ddg_summary.json`
  - generated at: `2026-06-10T10:52:59Z`
  - mutation: `B:Y32F@antibody__B:V34I@antibody`
  - mutation count / preset / side: `2` / `double_point` / `antibody`
  - `ddG = -2.11915462141458 kcal/mol`
  - ddG BAR stderr: `7.788108515900114 kcal/mol`
  - QC: `warning`
  - execution checkpoint passed: `True`
  - warnings:
    - complex:rep01 overlap score 0.148 below threshold 0.250
    - apo:rep01 overlap score 0.018 below threshold 0.250
- 4DN4 larger real-case quick validation:
  - file: `runs/real_cases/4dn4_v47i_quick/jobs/4dn4-antigen-m-v47i/results/ddg_summary.json`
  - generated at: `2026-06-10T11:50:54Z`
  - mutation: `M:V47I@antigen`
  - mutation count / preset / side: `1` / `single_point` / `antigen`
  - `ddG = 18.172982945270387 kcal/mol`
  - ddG BAR stderr: `55.64457153440632 kcal/mol`
  - QC: `warning`
  - execution checkpoint passed: `True`
  - warnings:
    - complex:rep01 overlap score 0.000 below threshold 0.200
    - complex mean BAR stderr 49.350 kcal/mol exceeds threshold 10.000
    - apo:rep01 overlap score 0.000 below threshold 0.200
    - apo mean BAR stderr 25.707 kcal/mol exceeds threshold 10.000
    - ddG BAR stderr 55.645 kcal/mol exceeds threshold 10.000

## External Regression Gate

- summary file: `runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference/reports/patel_2021_3hfm_summary.json`
- generated at: `2026-06-25T05:43:58Z`
- status: `insufficient_pairs`
- paired jobs: `0`
- incomplete jobs: `8`
- charge-conserving paired / incomplete: `0` / `3`
- charge-changing paired / incomplete: `0` / `5`
- external reference checkpoint passed: `False`
- note: `No completed Patel 2021 3HFM jobs are ready for external regression comparison yet. The charge-conserving subset is still in progress, while the charge-changing subset remains incomplete.`

## Live Execution State

- active core AB-Bind `gmx mdrun` processes: `2`
- active core AB-Bind `abag-rbfe resume` processes: `0`
- unique core active job ids (`resume` / `mdrun`): `0` / `2`
- orphaned core `resume` job ids: `none`
- stale core `gmx mdrun` processes (threshold `900` s): `0`
- active core benchmark roots: `abbind_core_v1_validation_target_specific_sampling_pilot_20260625`
- active untracked core benchmark roots: `none`
- active reference `gmx mdrun` processes: `0`
- active reference `abag-rbfe resume` processes: `0`
- orphaned reference `resume` job ids: `none`
- stale reference `gmx mdrun` processes: `0`
- active watcher processes: `0`

## Tracked Plan Drain State

- deep:
  - generated at: `2026-06-23T08:22:40Z`
  - selected / ddg-ready / paired: `19` / `3` / `3`
  - running sample / equilibrate: `9` / `0`
  - pending selected jobs: `16`
  - drained: `False`
- priority:
  - generated at: `2026-06-25T07:31:37Z`
  - selected / ddg-ready / paired: `80` / `0` / `0`
  - running sample / equilibrate: `0` / `0`
  - pending selected jobs: `80`
  - drained: `False`
- rescue:
  - generated at: `2026-06-23T08:22:24Z`
  - selected / ddg-ready / paired: `15` / `12` / `12`
  - running sample / equilibrate: `0` / `0`
  - pending selected jobs: `3`
  - drained: `False`
- robust:
  - generated at: `2026-06-23T08:22:46Z`
  - selected / ddg-ready / paired: `80` / `0` / `0`
  - running sample / equilibrate: `13` / `0`
  - pending selected jobs: `80`
  - drained: `False`
- sampling_qc:
  - generated at: `2026-06-25T07:30:00Z`
  - selected / ddg-ready / paired: `9` / `0` / `0`
  - running sample / equilibrate: `0` / `0`
  - pending selected jobs: `9`
  - drained: `False`
- target_specific_sampling_pilot_20260625:
  - generated at: `2026-06-26T05:32:36Z`
  - selected / ddg-ready / paired: `14` / `0` / `0`
  - running sample / equilibrate: `2` / `0`
  - pending selected jobs: `14`
  - drained: `False`
- targeted_lambda:
  - generated at: `2026-06-23T08:23:08Z`
  - selected / ddg-ready / paired: `4` / `4` / `4`
  - running sample / equilibrate: `0` / `0`
  - pending selected jobs: `0`
  - drained: `True`
- targeted_repeat:
  - generated at: `2026-06-23T08:23:23Z`
  - selected / ddg-ready / paired: `11` / `11` / `11`
  - running sample / equilibrate: `0` / `0`
  - pending selected jobs: `0`
  - drained: `True`
- ultra:
  - generated at: `2026-06-23T08:23:25Z`
  - selected / ddg-ready / paired: `12` / `0` / `0`
  - running sample / equilibrate: `9` / `0`
  - pending selected jobs: `12`
  - drained: `False`

## Completion Verdict

- accepted independent validation passed: `True`
- required real-case checkpoints completed: `True`
- same-side double-point checkpoint completed: `True`
- external 3HFM reference completed: `False`
- no live core processes remain: `False`
- no live reference processes remain: `True`
- tracked plan roots drained: `False`
- project complete: `False`

Current blockers:

- core AB-Bind execution is still active (mdrun=2, resume=0)
- Patel-like external 3HFM reference regression is incomplete (status=insufficient_pairs, paired=0, incomplete=8)
- tracked plan roots still report pending or running work: deep, priority, rescue, robust, sampling_qc, target_specific_sampling_pilot_20260625, ultra
