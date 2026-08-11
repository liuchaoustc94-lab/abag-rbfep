"""Regression tests on real benchmark structures.

Locks the classification behaviour for documented real-world failure samples
(1YY9 sidechain-only incomplete residues, 3BN9 disulfide-adjacent clash) using
the actual PDB inputs under benchmarks/ab_bind/source/structures. These tests
are CPU-only and never invoke gmx/pmx.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from abag_rbfe.structure import (
    classify_incomplete_standard_residues,
    find_inter_residue_heavy_atom_clashes,
    find_intra_residue_heavy_atom_clashes,
    partition_inter_residue_sidechain_repairable_clashes,
)

STRUCTURES_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "ab_bind" / "source" / "structures"

PDB_1YY9 = STRUCTURES_DIR / "1YY9.pdb"
PDB_2NZ9 = STRUCTURES_DIR / "2NZ9.pdb"
PDB_3BN9 = STRUCTURES_DIR / "3BN9.pdb"


@pytest.mark.skipif(not PDB_1YY9.is_file(), reason="1YY9.pdb not available")
def test_1yy9_incomplete_residues_are_all_sidechain_only() -> None:
    """1YY9 was historically blocked at prepare, but all 23 incomplete standard
    residues are sidechain-only (repairable via PDBFixer), none block on
    missing backbone atoms."""
    incomplete = classify_incomplete_standard_residues(PDB_1YY9)

    assert len(incomplete) == 23
    blocking = [r for r in incomplete if r.get("blocking_prepare")]
    assert blocking == []
    for residue in incomplete:
        assert residue["missing_backbone_atoms"] == []
        assert residue["missing_sidechain_atoms"]


@pytest.mark.skipif(not PDB_2NZ9.is_file(), reason="2NZ9.pdb not available")
def test_2nz9_source_structure_has_no_preexisting_heavy_atom_clashes() -> None:
    """The 2NZ9 H1064A failure arises post-mutation; the WT source structure
    itself must be clash-free so the failure stays attributable to the mutant
    geometry, not the input."""
    assert find_intra_residue_heavy_atom_clashes(PDB_2NZ9) == []
    assert find_inter_residue_heavy_atom_clashes(PDB_2NZ9) == []


@pytest.mark.skipif(not PDB_3BN9.is_file(), reason="3BN9.pdb not available")
def test_3bn9_disulfide_adjacent_cys_clash_is_stratified_sidechain_sidechain() -> None:
    """3BN9 CYS H140 SG vs CYS H196 SG (0.83 A) is a disulfide-adjacent contact
    that must be stratified as sidechain_sidechain and remain repairable, not
    a backbone hard failure."""
    inter = find_inter_residue_heavy_atom_clashes(PDB_3BN9)

    assert len(inter) == 1
    issue = inter[0]
    assert (issue["chain_id"], issue["resseq"], issue["resname"]) == ("H", 140, "CYS")
    assert (issue["partner_chain_id"], issue["partner_resseq"], issue["partner_resname"]) == ("H", 196, "CYS")
    assert issue["clash_classes"] == ["sidechain_sidechain"]
    for clash in issue["clashes"]:
        assert clash["clash_class"] == "sidechain_sidechain"
        assert clash["atom_a_class"] == "sidechain"
        assert clash["atom_b_class"] == "sidechain"

    repairable, blocking = partition_inter_residue_sidechain_repairable_clashes(inter)
    assert len(repairable) == 1
    assert blocking == []
