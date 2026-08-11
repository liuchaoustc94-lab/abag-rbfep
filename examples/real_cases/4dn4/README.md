# 4DN4 Real Case

This example uses the crystal structure of the CNTO888 Fab in complex with CCL2:

- `PDB`: `4DN4`
- antibody chains: `H`, `L`
- antigen chain: `M`
- quick test mutation: `M:V47I@antigen`

The bundled `protocol.quick.yml` is a short end-to-end validation preset. It is
only for workflow verification, not for production ddG reporting.

Typical execution:

```bash
cd /mnt/data/liuchao/abag-rbfep
./examples/real_cases/4dn4/run_quick.sh
```

Outputs are written to:

```text
runs/real_cases/4dn4_v47i_quick/
```

Current validation status:

- the current workflow no longer blocks during `prepare`
- sidechain-only atom gaps in the extracted structure are repaired during
  `prepare`, and the rerun batch completes `mutate -> equilibrate -> sample ->
  bar -> qc -> report`
- latest quick-run result bundle:
  - `ddG = 18.173 kcal/mol`
  - `QC = warning`
  - warning drivers are the expected quick-preset sampling limits for this
    larger complex: `overlap score = 0.000` in both legs and large BAR stderr
- current stage/QC/result files:
  - `jobs/4dn4-antigen-m-v47i/stages/prepare.json`
  - `jobs/4dn4-antigen-m-v47i/stages/mutate.json`
  - `jobs/4dn4-antigen-m-v47i/stages/equilibrate.json`
  - `jobs/4dn4-antigen-m-v47i/stages/sample.json`
  - `jobs/4dn4-antigen-m-v47i/stages/bar.json`
  - `jobs/4dn4-antigen-m-v47i/stages/qc.json`
  - `jobs/4dn4-antigen-m-v47i/stages/report.json`
  - `jobs/4dn4-antigen-m-v47i/artifacts/prepare_qc.json`
  - `jobs/4dn4-antigen-m-v47i/artifacts/mutate_qc.json`
  - `jobs/4dn4-antigen-m-v47i/results/ddg_summary.json`
  - `jobs/4dn4-antigen-m-v47i/results/qc_report.json`
  - `jobs/4dn4-antigen-m-v47i/report/summary.json`
