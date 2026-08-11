# Vendored pmx

The `vendor/pmx` tree is a vendored snapshot of:

- Upstream: `https://github.com/deGrootLab/pmx`
- Branch: `develop`
- Resolved head at vendoring time: `e7b6a00c25c12cdd8fabe32929ec2dc2fee76ff3`

The project keeps pmx in a dedicated directory so future patches can be tracked
separately from the local orchestration code.

## Intended patch policy

- Antibody-antigen specific helper patches
- Batch execution and workflow-facing wrapper extensions
- GROMACS compatibility fixes needed by this project

Everything else should remain upstream-compatible.

## Current Project-Local Patches

The current real-case execution path required a small set of tolerance patches
for hybrid residues generated from crystal structures without a fully matching
hydrogen layout:

- `pmx/geometry.py`
  - `bb_super()` now copies only backbone atoms present in both residues, by
    atom name, instead of assuming identical main-chain atom lists
- `pmx/alchemy.py`
  - `_set_conformation()` now skips dihedral and anchor operations when the
    required atoms are absent
  - `_get_hybrid_residues()` now maps B-state topology attributes by atom name
    and skips atoms absent after `pdb2gmx -missing`
  - `_check_dih_ILDN_OPLS()` and `_find_predefined_dihedrals()` now skip
    predefined dihedral templates whose referenced atoms are absent

These patches were introduced to make the standalone workflow robust enough to
run the real `1VFB` quick cases end to end and to let `4DN4` advance from
sidechain-only gap repair through `prepare`/`mutate` and into a completed
end-to-end quick validation run instead of remaining blocked at input QC.
