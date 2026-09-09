# SKEMPI2 protein-protein single-point panel

This panel extends the validation pool beyond antibody-antigen complexes using
experimental SKEMPI2 v2 protein-protein systems. It contains ten new PDB
complexes, all with at least ten unique experimental single substitutions on
one selected partner. The counts and source labels are in
[`target_manifest.csv`](target_manifest.csv); structure quality is in
[`structure_audit.csv`](structure_audit.csv).

| PDB | Complex | Unique singles | Mutated partner |
|---|---|---:|---|
| 1C1Y | Rap1A–Raf RBD | 13 | Raf RBD |
| 1A4Y | ribonuclease inhibitor–angiogenin | 16 | ribonuclease inhibitor |
| 1C4Z | E6AP–UBCH7 | 22 | E6AP |
| 1DAN | factor VIIa–tissue factor | 85 | tissue factor |
| 1E50 | AML1–CBFβ | 11 | AML1 |
| 1EMV | colicin E9–Im9 | 34 | Im9 |
| 1F47 | FtsZ fragment–ZipA | 12 | FtsZ fragment |
| 1JTD | TEM-1–BLIP-II | 25 | BLIP-II |
| 3QHY | beta-lactamase–BLIP-II | 25 | BLIP-II |
| 1FFW | CheY–CheA | 11 | CheY |

## Schema compatibility

The current planner accepts the two sides under the historical
`antibody`/`antigen` names. For this panel, `antibody_chains` means the
selected mutated partner-A chains and `antigen_chains` means the binding
partner-B chains; `system.yml` records `complex_class: protein-protein` and
`mutated_partner: partner_a`. This is a schema alias only and does not claim
that these proteins are antibodies.

Use `mutations_single_unique.csv` for a single-point V1 plan. The source
`experimental_ddg.csv` retains duplicate measurements, multi-point rows,
affinities, and censor flags; do not mix the multi-point rows into the V1
single-point score.

## Rebuild and plan

```bash
python tools/build_skempi2_protein_protein_panel.py
```

Example plan (no GROMACS execution):

```bash
./.venv/bin/abag-rbfe batch plan \
  --system benchmarks/protein_protein_skempi2_single_point/targets/1C1Y/system.yml \
  --mutations benchmarks/protein_protein_skempi2_single_point/targets/1C1Y/mutations_single_unique.csv \
  --protocol benchmarks/ab_bind/protocol.validation.yml \
  --batch-id skempi2_pp_1c1y
```

The ten systems are experimental X-ray complexes. See the audit for
resolution and construct notes. Report this protein-protein stratum separately
from the antibody-antigen stratum when assessing generalization.

## Affinity-change statistics

Generate the per-target and per-mutation affinity tables with:

```bash
python tools/summarize_skempi2_affinity_changes.py
```

The resulting `affinity_change_report.md` uses unique single-point signatures;
`affinity_change_summary.csv` contains target-level ranges and gain/loss
counts, while `affinity_change_single_unique.csv` contains each mutation's
WT Kd, mutant Kd, fold change, ΔpKd, and ΔΔG.
