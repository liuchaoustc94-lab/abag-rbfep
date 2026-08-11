# 1VFB Real Case

This example uses the D1.3 Fv - hen egg white lysozyme complex:

- `PDB`: `1VFB`
- antibody chains: `A`, `B`
- antigen chain: `C`
- quick test mutations:
  - `V1`: `B:Y32F@antibody`
  - `V2`: `B:Y32F@antibody + B:V34I@antibody`

`protocol.quick.yml` is a short workflow-validation preset. It is suitable for
end-to-end execution checks, not for production-quality ddG estimates.

`protocol.quick.double.yml` is the matching short preset for the same-side
double-point validation case.

Typical execution:

```bash
cd /mnt/data/liuchao/abag-rbfep
./examples/real_cases/1vfb/run_quick.sh
./examples/real_cases/1vfb/run_quick_double.sh
```

Outputs are written to:

```text
runs/real_cases/1vfb_y32f_quick/
runs/real_cases/1vfb_y32f_v34i_quick/
```

Current validation status:

- `V1` single-point:
  - stage path: `runs/real_cases/1vfb_y32f_quick/jobs/1vfb-antibody-b-y32f/stages/`
  - latest quick-run result: `ddG = -3.591 kcal/mol`
  - latest quick-run QC: `warning`
  - note: `ddG` is ready, but the short quick preset leaves BAR overlap below the QC threshold
- `V2` same-side double-point:
  - stage path: `runs/real_cases/1vfb_y32f_v34i_quick/jobs/1vfb-antibody-b-y32f--b-v34i/stages/`
  - latest quick-run result: `ddG = -2.119 kcal/mol`
  - latest quick-run QC: `warning`
  - note: `ddG` is ready, but the short quick preset leaves BAR overlap below the QC threshold
- formal outputs for both runs:
  - `results/ddg_summary.json`
  - `results/qc_report.json`
  - `report/summary.json`
