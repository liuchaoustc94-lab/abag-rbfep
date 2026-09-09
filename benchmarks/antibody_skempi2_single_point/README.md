# SKEMPI2-only antibody-side single-point panel

This is the experimental-only subset selected directly from the local
SKEMPI2 v2 table. It contains no Graphinity or other synthetic labels. The
authoritative counts are in [`target_manifest.csv`](target_manifest.csv), and
the PDB method/resolution audit is in [`structure_audit.csv`](structure_audit.csv).

| Target | Unique antibody-side singles | Status |
|---|---:|---|
| 3SE8, VRC03–gp120 | 28 | strict experimental |
| 3SE9, VRC-PG04–gp120 | 25 | strict experimental |
| 3N85, Fab37–HER2 | 9 | near threshold; one duplicate source row |
| 5C6T, 1G2–HCMV gB | 8 | near threshold; one combination row retained separately |

`mutations_single_unique.csv` is the recommended input for the V1 planner.
`mutations_single_all.csv` preserves repeated source measurements, and
`experimental_ddg.csv` contains the recomputed `RT ln(Kd_mut/Kd_wt)` labels,
source affinities, and censoring flags. All selected sites are on antibody
chains and passed the PDB residue/WT audit.

## Rebuild

```bash
python tools/build_skempi2_single_point_panel.py
```

For example, plan the 3SE8 experiment without executing GROMACS:

```bash
./.venv/bin/abag-rbfe batch plan \
  --system benchmarks/antibody_skempi2_single_point/targets/3SE8/system.yml \
  --mutations benchmarks/antibody_skempi2_single_point/targets/3SE8/mutations_single_unique.csv \
  --protocol benchmarks/ab_bind/protocol.validation.yml \
  --batch-id skempi2_single_3se8
```

Charge-changing substitutions require the charge-enabled DSSB protocol
preset. Do not pool the two near-threshold targets with the strict targets
without reporting the duplicate/low-count caveat.

An exhaustive untried SKEMPI2 candidate audit is available at
[`../antibody_new10_single_point/skempi2_untried_antibody_candidates.csv`](../antibody_new10_single_point/skempi2_untried_antibody_candidates.csv).
It shows that SKEMPI2 alone does not contain ten new independent targets with
at least ten unique antibody-side singles after excluding the project’s
already-used biology.
