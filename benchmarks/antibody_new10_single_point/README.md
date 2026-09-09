# Single-point antibody-side validation panel

This is the single-point-focused companion to `antibody_new10`. It contains
ten new antibody–antigen structures and only antibody-side one-residue
substitutions in the primary mutation files. Counts and provenance are in
[`target_manifest.csv`](target_manifest.csv); structure method, resolution,
and WT/construct caveats are in [`structure_audit.csv`](structure_audit.csv).

## What is suitable for model validation?

The public experimental pool is smaller than the requested ten-target panel:

| Tier | Targets | Mutation records | Use |
|---|---|---:|---|
| Strict experimental | 3SE8, 3SE9 | 28, 25 unique singles | Primary held-out ddG validation |
| Experimental near-threshold | 3N85 | 10 rows, 9 unique singles | Report separately; one duplicate source measurement |
| Synthetic controls | 3J1S, 4KUC, 5JHL, 4NNP, 4PLJ, 6A79, 5J56 | 15–16 each | Pipeline, residue mapping, and sampling checks only |

Synthetic labels are Graphinity Rosetta Flex ddG values, not experimental
affinity measurements. Never pool them with SKEMPI2 rows when calculating
Pearson `R`, MAE, or a publication claim. `6A79` contains a deposited P103A
scFv heavy-chain variant; restore the intended WT sequence before treating it
as an RBFE baseline. `3J1S` is an 8.5 Å EM model rather than a crystal
structure. The remaining entries are X-ray structures (see the audit file).

`mutations_single_unique.csv` is the recommended V1 input. The corresponding
`mutations_single_all.csv` preserves duplicate measurements, while
`experimental_ddg.csv` retains the source label, affinity values, and censor
flags. All selected sites passed a residue-name/WT audit against the copied
PDB files.

## Rebuild and plan

Refresh the local SKEMPI2/Graphinity inputs and rebuild with:

```bash
python tools/build_single_point_panel.py
```

Plan a primary experimental target without running GROMACS:

```bash
./.venv/bin/abag-rbfe batch plan \
  --system benchmarks/antibody_new10_single_point/targets/3SE8/system.yml \
  --mutations benchmarks/antibody_new10_single_point/targets/3SE8/mutations_single_unique.csv \
  --protocol benchmarks/ab_bind/protocol.validation.yml \
  --batch-id single_point_3se8
```

Charge-changing substitutions require the charge-enabled DSSB protocol
preset; do not silently use the ordinary validation protocol for those rows.
Keep synthetic controls in a separately labelled report namespace.

The reserve experimental candidate `5C6T` has only eight unique antibody-side
singles and is recorded in `reserve_candidates.csv` rather than counted as a
ten-point target.

## SKEMPI2-only search result

[`skempi2_untried_antibody_candidates.csv`](skempi2_untried_antibody_candidates.csv)
records the exhaustive untried candidate scan from the local SKEMPI2 v2 table,
using the PDB mutation-chain annotation and RCSB entity descriptions to keep
only antibody-side mutations. Under the current exclusions, only 3SE8 and
3SE9 satisfy the requirement of at least ten unique experimental singles;
3N85 and 5C6T are near-threshold. The other candidates either contain fewer
than ten singles, use a mutant/non-WT deposited baseline, or repeat biology
already used by the project. Therefore SKEMPI2 alone cannot supply ten new
strict targets under all the original constraints. The seven synthetic rows in
this directory are controls, not a substitute for missing SKEMPI2 experiments.
