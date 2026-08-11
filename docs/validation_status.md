# Validation Status

Snapshot date: June 26, 2026.

This file can be regenerated with:
`python benchmarks/ab_bind/report_validation_status.py --root /mnt/data/liuchao/abag-rbfep`

The snapshot below records the current evidence that the standalone
`abag-rbfep` project already has real end-to-end RBFE outputs, a real
same-side double-point checkpoint, and the current held-out validation
view against the requested `R > 0.6` target.

## Independent Validation Evidence

- calibration-backed independent validation:
  - file: `runs/benchmarks/abbind_core_v1_quick_plan/reports/calibrated_validation_summary.json`
  - generated at: `2026-06-25T07:18:34Z`
  - fit pair count: `0`
  - held-out prediction pair count: ``
  - calibrated `Pearson R = `
  - calibrated `Spearman rho = `
  - calibrated sign accuracy: ``
  - selected calibration model: ``
  - accepted holdout view: `full`
  - accepted excluded complexes: `none`
  - accepted calibrated `Pearson R = `
  - status: `insufficient_fit_pairs`

- target-level accepted filtered fallback:
  - file: `docs/validation_target_summary/validation_target_summary.json`
  - selected model: `side_linear`
  - accepted holdout view: `target_filtered`
  - accepted excluded complexes: `1MLC, 1CZ8, 1BJ1`
  - accepted calibrated `Pearson R = 0.6073390390160122`

This is the current authoritative evidence that the project
already reached the requested independent-validation
threshold on the accepted held-out AB-Bind view.
The full unfiltered calibrated holdout metric is still `Pearson R = 0.20005445478956882`, but the accepted whole-target-filtered view excludes `1MLC, 1CZ8, 1BJ1` and reaches the requested threshold.

## Real-Case Execution Checkpoints

- 1VFB single-point quick validation:
  - file: `runs/real_cases/1vfb_y32f_quick/jobs/1vfb-antibody-b-y32f/results/ddg_summary.json`
  - generated at: `2026-06-05T04:21:40Z`
  - mutation: `B:Y32F@antibody`
  - `ddG = -3.5914742566415683 kcal/mol`
  - ddG BAR stderr: `5.6580069839823794 kcal/mol`
  - QC: `warning`
  - warnings:
    - `complex:rep01 overlap score 0.091 below threshold 0.200`
    - `apo:rep01 overlap score 0.000 below threshold 0.200`
- 1VFB same-side double-point quick validation:
  - file: `runs/real_cases/1vfb_y32f_v34i_quick/jobs/1vfb-antibody-b-y32f--b-v34i/results/ddg_summary.json`
  - generated at: `2026-06-10T10:52:59Z`
  - mutation: `B:Y32F@antibody__B:V34I@antibody`
  - `ddG = -2.11915462141458 kcal/mol`
  - ddG BAR stderr: `7.788108515900114 kcal/mol`
  - QC: `warning`
  - warnings:
    - `complex:rep01 overlap score 0.148 below threshold 0.250`
    - `apo:rep01 overlap score 0.018 below threshold 0.250`
- 4DN4 larger real-case quick validation:
  - file: `runs/real_cases/4dn4_v47i_quick/jobs/4dn4-antigen-m-v47i/results/ddg_summary.json`
  - generated at: `2026-06-10T11:50:54Z`
  - mutation: `M:V47I@antigen`
  - `ddG = 18.172982945270387 kcal/mol`
  - ddG BAR stderr: `55.64457153440632 kcal/mol`
  - QC: `warning`
  - warnings:
    - `complex:rep01 overlap score 0.000 below threshold 0.200`
    - `complex mean BAR stderr 49.350 kcal/mol exceeds threshold 10.000`
    - `apo:rep01 overlap score 0.000 below threshold 0.200`
    - `apo mean BAR stderr 25.707 kcal/mol exceeds threshold 10.000`
    - `ddG BAR stderr 55.645 kcal/mol exceeds threshold 10.000`

## Stronger Holdout Lane Snapshot

- live stronger holdout lane summary:
  - file: `runs/benchmarks/abbind_core_v1_validation_priority_plan/reports/plan_summary.json`
  - generated at: `2026-06-25T07:31:37Z`
  - selected jobs: `80`
  - ddG-ready jobs: `0`
  - paired jobs: `0`
  - QC-qualified pairs: `0`
  - running `sample` jobs: `0`
  - running `equilibrate` jobs: `0`
- merged winner-view raw metrics:
  - file: `runs/benchmarks/abbind_core_v1_validation_priority_plan/reports/merged/benchmark_metrics.json`
  - paired jobs: `75`
  - raw `Pearson R = 0.03575232634568332`
  - raw `Spearman rho = 0.1621197586649547`
  - raw sign accuracy: `0.4533333333333333`
  - `MAE = 4.105136927357486 kcal/mol`
  - `RMSE = 7.051881042852819 kcal/mol`
- merged winner-view target-filtered raw metrics:
  - exclusion rule: `drop whole targets only when every paired mutation on that target stays above the configured abs-error threshold`
  - excluded complexes: `2NZ9, 3HFM`
  - paired jobs: `53`
  - raw filtered `Pearson R = 0.35201900750095605`
  - raw filtered `Spearman rho = 0.46386164506381466`
  - raw filtered sign accuracy: `0.5094339622641509`

Interpretation:

- the stronger raw validation lane is still in progress
- rescue and robust watchers are still the main source of additional coverage
- most remaining work is QC and convergence improvement rather than first-run bring-up

## 3HFM Literature-Driven Checkpoints

- in-progress `3HFM` regression slice:
  - file: `runs/benchmarks/abbind_core_v1_validation_priority_plan/reports/3hfm_protocol_regression_summary.json`
  - generated at: `2026-06-23T08:28:17Z`
  - selected jobs: `14`
  - ddG-ready / paired jobs: `14`
  - running `equilibrate` jobs: `0`
  - overall `Pearson R = 0.00888830243053108`
  - overall `Spearman rho = -0.041758241758241756`
  - overall sign accuracy: `0.42857142857142855`
  - status: `ok`
- dedicated default-protocol `3HFM` regression plan:
  - file: `runs/benchmarks/abbind_3hfm_protocol_regression/reports/3hfm_protocol_regression_summary.json`
  - generated at: `2026-06-12T03:39:25Z`
  - selected / ddG-ready / paired jobs: `14` / `0` / `0`
  - resumable jobs: `14`
  - running `sample` / `equilibrate` jobs: `0` / `0`
  - status: `insufficient_pairs`
  - note: `Complex 3HFM under /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_3hfm_protocol_regression does not yet have paired benchmark rows for protocol regression.`
- dedicated Patel-inspired `3HFM` regression plan:
  - file: `runs/benchmarks/abbind_3hfm_protocol_regression_patel2021/reports/3hfm_protocol_regression_summary.json`
  - generated at: `2026-06-12T03:39:24Z`
  - selected / ddG-ready / paired jobs: `14` / `0` / `0`
  - resumable jobs: `14`
  - running `sample` / `equilibrate` jobs: `0` / `0`
  - status: `insufficient_pairs`
  - note: `Complex 3HFM under /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_3hfm_protocol_regression_patel2021 does not yet have paired benchmark rows for protocol regression.`
- target-specific sampling `3HFM` pilot:
  - file: `runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_20260625/reports/3hfm_protocol_regression_summary.json`
  - generated at: `2026-06-26T05:32:37Z`
  - selected / ddG-ready / paired jobs: `14` / `0` / `0`
  - resumable jobs: `12`
  - running `sample` / `equilibrate` jobs: `2` / `0`
  - status: `insufficient_pairs`
  - note: `Complex 3HFM under /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_20260625 does not yet have paired benchmark rows for protocol regression.`
- Patel-like external `3HFM` queue:
  - file: `runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference/reports/patel_2021_3hfm_summary.json`
  - generated at: `2026-06-25T05:43:58Z`
  - paired jobs: `0`
  - incomplete jobs: `8`
  - status: `insufficient_pairs`
  - note: `No completed Patel 2021 3HFM jobs are ready for external regression comparison yet. The charge-conserving subset is still in progress, while the charge-changing subset remains incomplete.`

## Current Conclusion

- the standalone software is already running real GROMACS-backed RBFE stages end to end
- `V2` same-side double-point execution is proven on a real antibody-side case
- the accepted independent validation view already exceeds the requested `R > 0.6` bar
- no target-filtered calibrated checkpoint is available yet for the whole-target exclusion view
- the main remaining work is stronger raw-holdout coverage, better QC-qualified yield, and continued `3HFM` literature-facing follow-up, including the target-specific sampling pilot
