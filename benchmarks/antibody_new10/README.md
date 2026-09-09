# New antibody-side validation panel

This directory contains a conservative, provenance-tracked panel of ten
structures not present in the current formal/acceptance/materialized target
set. Mutations are restricted to antibody chains. The authoritative counts are
in [`target_manifest.csv`](target_manifest.csv); experimental method, resolution
and WT/construct status are recorded in [`structure_audit.csv`](structure_audit.csv).

## Important interpretation

The panel has two label domains and they must not be pooled into one Pearson
`R`:

| Tier | Targets | Use |
|---|---|---|
| Experimental, single-point ready | 3SE8 (28), 3SE9 (25) | Primary out-of-sample ddG validation |
| Experimental, near/exploratory | 3N85 (10 rows/9 unique), 3G6D (10), 3L5X (36), 2B2X (115) | 3N85 is a duplicate-measurement stress case; the other three are mainly multi-point and require V2/group-aware handling |
| Synthetic control | 3J1S (15), 4KUC (16), 5JHL (16), 4NNP (16) | Structure, mutation-routing, and sampling checks only; labels are Graphinity Rosetta Flex ddG, not measurements |

`2B2X` is a distinct structure of the AQC2–integrin-alpha-1 system and is
flagged as related to the already queued `1MHP`; do not describe it as an
independent biological target. The deposited 2B2X Fab is affinity-matured, so
its sequence must be reverted to a WT AQC2 baseline before an RBFE run.
`5C6T` is retained in
[`reserve_candidates.csv`](reserve_candidates.csv): it has nine antibody-side
records (eight unique singles), so it does not meet the ten-point rule yet.

## File layout

- `targets/<PDB>/system.yml`: chain-aware system manifest.
- `targets/<PDB>/experimental_ddg.csv`: one row per source mutation group,
  including the ddG label and censoring flag. For synthetic targets the same
  file stores the Rosetta control label and says so in `label_type`.
- `mutations_single_unique.csv`: deduplicated V1-compatible single-point sites.
- `mutations_all.csv`: all selected antibody-side sites, including combinations
  that are not accepted by the current V1 planner.
- `structures/`: cleaned SKEMPI2 structures for experimental targets and RCSB
  structures for the controls.
- `structure_audit.csv`: X-ray/EM method, resolution, and whether the deposited
  complex can be treated as a WT baseline.
- `source_metadata/sources.csv`: URLs, retrieval date, and SHA-256 checksums.

Experimental ddG is recomputed from SKEMPI2 affinities as
`RT ln(Kd_mut/Kd_wt)` using the row temperature. Missing/non-binder affinities
remain censored and are never silently converted to a numeric value.

For the synthetic controls, `source_mutation_tokens` retains Graphinity's
cleaned numbering while `mutation_tokens` and the site CSVs are remapped to the
downloaded RCSB PDB residue IDs by sequence alignment. A residue/WT audit was
run for every selected site (zero missing or WT-mismatch sites).

## Rebuild and run

Rebuild the selected artifacts after refreshing the two upstream files:

```bash
python tools/build_new_antibody_validation_panel.py
```

Plan the two primary experimental targets without executing GROMACS:

```bash
./.venv/bin/abag-rbfe batch plan \
  --system benchmarks/antibody_new10/targets/3SE8/system.yml \
  --mutations benchmarks/antibody_new10/targets/3SE8/mutations_single_unique.csv \
  --protocol benchmarks/ab_bind/protocol.validation.yml \
  --batch-id new10_3se8
```

Repeat with `3SE9`. Use `mutations_all.csv` only for a separately labelled
multi-point exploration; the current planner accepts single and double groups,
so the >2-site groups in 3G6D, 3L5X, and 2B2X must first be decomposed or
handled by a future multi-point extension. Do not mix them into the V1 score.

## Sources

- SKEMPI2 v2 CSV and cleaned PDB archive (Jankauskaite et al., 2019), URLs in
  `source_metadata/sources.csv`.
- Graphinity Flex ddG control CSV (Hummer et al., 2025), also recorded in the
  source metadata. Synthetic labels are explicitly excluded from experimental
  model-quality claims.
