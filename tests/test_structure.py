from pathlib import Path

from abag_rbfe.structure import (
    classify_incomplete_standard_residues,
    extract_pdb_chains,
    find_incomplete_standard_residues,
    find_inter_residue_heavy_atom_clashes,
    find_intra_residue_heavy_atom_clashes,
    partition_sidechain_repairable_clashes,
    pdb_chain_ids,
    repair_missing_atoms_with_pdbfixer,
    repair_incomplete_standard_residues_with_pdbfixer,
    repair_sidechain_only_incomplete_residues_with_pdbfixer,
    restore_incomplete_standard_residues_from_template,
    strip_sidechain_atoms_for_residues,
    strip_terminal_oxygen_atoms,
    write_inter_residue_heavy_atom_clash_report,
)


def test_structure_helpers_keep_only_protein_atom_records(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "HETATM    2  C1  GOL H 401      12.000  14.000   8.000  1.00 20.00           C",
                "HETATM    3  O   HOH H 501      13.000  14.500   8.500  1.00 20.00           O",
                "TER",
                "ATOM      4  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    assert pdb_chain_ids(input_path) == {"H", "A"}

    extract_pdb_chains(input_path, output_path, keep_chains=["H"])
    output = output_path.read_text(encoding="utf-8")
    assert "TYR H  32" in output
    assert "GOL" not in output
    assert "HOH" not in output
    assert "LYS A  58" not in output


def test_extract_pdb_chains_keeps_single_preferred_altloc_and_clears_altloc_label(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   LEU H  18      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA ALEU H  18      12.200  10.100   8.700  0.40 20.00           C",
                "ATOM      3  CA BLEU H  18      12.320  10.180   8.780  0.60 20.00           C",
                "ATOM      4  C   LEU H  18      13.410   9.280   8.960  1.00 20.00           C",
                "ATOM      5  O   LEU H  18      13.640   8.160   8.520  1.00 20.00           O",
                "ATOM      6  CB ALEU H  18      11.860  11.520   9.090  0.40 20.00           C",
                "ATOM      7  CB BLEU H  18      12.020  11.610   9.180  0.60 20.00           C",
                "ATOM      8  CG ALEU H  18      12.620  12.660   8.430  0.40 20.00           C",
                "ATOM      9  CG BLEU H  18      12.790  12.760   8.540  0.60 20.00           C",
                "ATOM     10 CD1ALEU H  18      12.330  14.040   9.000  0.40 20.00           C",
                "ATOM     11 CD1BLEU H  18      12.500  14.140   9.100  0.60 20.00           C",
                "ATOM     12 CD2ALEU H  18      14.070  12.450   8.640  0.40 20.00           C",
                "ATOM     13 CD2BLEU H  18      14.250  12.550   8.740  0.60 20.00           C",
                "TER",
                "ATOM     14  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    assert find_intra_residue_heavy_atom_clashes(input_path) == []

    extract_pdb_chains(input_path, output_path, keep_chains=["H"])
    output = output_path.read_text(encoding="utf-8")
    assert output.count(" CA  LEU H  18 ") == 1
    assert output.count(" CB  LEU H  18 ") == 1
    assert " CA ALEU H  18 " not in output
    assert " CA BLEU H  18 " not in output
    assert "  12.320  10.180   8.780" in output
    assert "  12.020  11.610   9.180" in output


def test_find_incomplete_standard_residues_reports_missing_heavy_atoms(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   GLU H   1      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  GLU H   1      12.936  14.622   9.024  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    issues = find_incomplete_standard_residues(input_path)
    assert issues == [
        {
            "chain_id": "H",
            "resseq": 1,
            "icode": "",
            "resname": "GLU",
            "normalized_resname": "GLU",
            "missing_atoms": ["CD", "CG", "OE1", "OE2"],
            "missing_backbone_atoms": [],
            "missing_sidechain_atoms": ["CD", "CG", "OE1", "OE2"],
            "blocking_prepare": False,
        }
    ]


def test_classify_incomplete_standard_residues_marks_backbone_missing_as_blocking(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  CB  GLU H   1      12.936  14.622   9.024  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    issues = classify_incomplete_standard_residues(input_path)
    assert issues == [
        {
            "chain_id": "H",
            "resseq": 1,
            "icode": "",
            "resname": "GLU",
            "normalized_resname": "GLU",
            "missing_atoms": ["CD", "CG", "O", "OE1", "OE2"],
            "missing_backbone_atoms": ["O"],
            "missing_sidechain_atoms": ["CD", "CG", "OE1", "OE2"],
            "blocking_prepare": True,
        }
    ]


def test_restore_incomplete_standard_residues_from_template_replaces_missing_sidechains(tmp_path: Path) -> None:
    template_path = tmp_path / "template.pdb"
    target_path = tmp_path / "target.pdb"
    restored_path = tmp_path / "restored.pdb"
    template_path.write_text(
        "\n".join(
            [
                "HEADER    TEMPLATE",
                "ATOM      1  N   GLU A   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU A   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU A   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   GLU A   1      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  GLU A   1      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  GLU A   1      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD  GLU A   1      12.402  16.521  10.579  1.00 20.00           C",
                "ATOM      8  OE1 GLU A   1      13.535  16.885  10.982  1.00 20.00           O",
                "ATOM      9  OE2 GLU A   1      11.503  17.327  10.495  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    target_path.write_text(
        "\n".join(
            [
                "HEADER    TARGET",
                "ATOM     99  N   GLU A   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM    100  CA  GLU A   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM    101  C   GLU A   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM    102  O   GLU A   1      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM    103  CB  GLU A   1      12.936  14.622   9.024  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = restore_incomplete_standard_residues_from_template(template_path, target_path, restored_path)
    restored_issues = classify_incomplete_standard_residues(restored_path)
    restored_text = restored_path.read_text(encoding="utf-8")

    assert summary["attempted"] is True
    assert summary["restored_residue_count"] == 1
    assert summary["missing_template_residue_count"] == 0
    assert summary["unresolved_residue_count"] == 0
    assert restored_issues == []
    assert "ATOM      1  N   GLU A   1" in restored_text
    assert "ATOM      9  OE2 GLU A   1" in restored_text


def test_strip_sidechain_atoms_for_residues_keeps_only_backbone(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   GLU H   1      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  GLU H   1      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  GLU H   1      12.086  15.074  10.191  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    strip_sidechain_atoms_for_residues(
        input_path,
        input_path,
        [{"chain_id": "H", "resseq": 1, "icode": ""}],
    )

    output = input_path.read_text(encoding="utf-8")
    assert " N   GLU H   1 " in output
    assert " CA  GLU H   1 " in output
    assert " C   GLU H   1 " in output
    assert " O   GLU H   1 " in output
    assert " CB  GLU H   1 " not in output
    assert " CG  GLU H   1 " not in output


def test_repair_missing_atoms_with_pdbfixer_restores_sidechain_only_gap(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   GLU H   1      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = repair_missing_atoms_with_pdbfixer(input_path, input_path)
    issues = classify_incomplete_standard_residues(input_path)

    assert summary["attempted"] is True
    assert summary["available"] is True
    assert summary["succeeded"] is True
    assert issues == []


def test_repair_sidechain_only_incomplete_residues_with_pdbfixer_strips_and_restores_gap(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   GLU H   1      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  GLU H   1      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  GLU H   1      12.086  15.074  10.191  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = repair_sidechain_only_incomplete_residues_with_pdbfixer(input_path, input_path)
    issues = classify_incomplete_standard_residues(input_path)

    assert summary["attempted"] is True
    assert summary["available"] is True
    assert summary["succeeded"] is True
    assert summary["trigger_residue_count"] == 1
    assert summary["remaining_incomplete_standard_residue_count"] == 0
    assert issues == []


def test_repair_incomplete_standard_residues_with_pdbfixer_restores_backbone_oxygen_gap(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  CB  GLU H   1      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      5  CG  GLU H   1      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      6  CD  GLU H   1      12.402  16.521  10.579  1.00 20.00           C",
                "ATOM      7  OE1 GLU H   1      13.535  16.885  10.982  1.00 20.00           O",
                "ATOM      8  OE2 GLU H   1      11.503  17.327  10.495  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = repair_incomplete_standard_residues_with_pdbfixer(input_path, input_path)
    issues = classify_incomplete_standard_residues(input_path)

    assert summary["attempted"] is True
    assert summary["available"] is True
    assert summary["succeeded"] is True
    assert summary["trigger_residue_count"] == 1
    assert summary["remaining_incomplete_standard_residue_count"] == 0
    assert issues == []


def test_strip_terminal_oxygen_atoms_removes_oxt_like_records(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   ARG H   1      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  ARG H   1      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   ARG H   1      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   ARG H   1      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5 OXT  ARG H   1      13.900  11.900   9.800  1.00 20.00           O",
                "ATOM      6 OT1  ARG H   1      14.100  11.700  10.200  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    strip_terminal_oxygen_atoms(input_path, input_path)
    output = input_path.read_text(encoding="utf-8")
    assert " OXT " not in output
    assert " OT1 " not in output
    assert " N   ARG H   1 " in output


def test_strip_terminal_oxygen_atoms_preserves_insertion_code_residues(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   PRO H 52A      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  PRO H 52A      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   PRO H 52A      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   PRO H 52A      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5 OXT  PRO H 52A      13.900  11.900   9.800  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    strip_terminal_oxygen_atoms(input_path, input_path)
    output = input_path.read_text(encoding="utf-8")
    assert " OXT " not in output
    assert " PRO H  52A" in output


def test_find_incomplete_standard_residues_preserves_insertion_code_labels(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H 52A      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  GLU H 52A      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   GLU H 52A      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   GLU H 52A      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  GLU H 52A      12.936  14.622   9.024  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    issues = find_incomplete_standard_residues(input_path)
    assert issues == [
        {
            "chain_id": "H",
            "resseq": 52,
            "icode": "A",
            "resname": "GLU",
            "normalized_resname": "GLU",
            "missing_atoms": ["CD", "CG", "OE1", "OE2"],
            "missing_backbone_atoms": [],
            "missing_sidechain_atoms": ["CD", "CG", "OE1", "OE2"],
            "blocking_prepare": False,
        }
    ]


def test_find_intra_residue_heavy_atom_clashes_reports_local_geometry_breaks(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  GLU H   1      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   GLU H   1      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   GLU H   1      13.350   9.250   9.700  1.00 20.00           O",
                "ATOM      5  CB  GLU H   1      12.400  11.500   9.200  1.00 20.00           C",
                "ATOM      6  CG  GLU H   1      13.400   9.350   9.100  1.00 20.00           C",
                "ATOM      7  CD  GLU H   1      14.300  12.800  10.000  1.00 20.00           C",
                "ATOM      8  OE1 GLU H   1      14.900  12.900  10.900  1.00 20.00           O",
                "ATOM      9  OE2 GLU H   1      13.900  13.600   9.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    clashes = find_intra_residue_heavy_atom_clashes(input_path)
    assert clashes == [
        {
            "chain_id": "H",
            "resseq": 1,
            "icode": "",
            "resname": "GLU",
            "normalized_resname": "GLU",
            "min_distance_angstrom": 0.2693,
            "clashes": [
                {
                    "atom_a": "C",
                    "atom_b": "O",
                    "distance_angstrom": 0.8031,
                    "atom_a_class": "backbone",
                    "atom_b_class": "backbone",
                    "clash_class": "backbone_backbone",
                },
                {
                    "atom_a": "C",
                    "atom_b": "CG",
                    "distance_angstrom": 0.2693,
                    "atom_a_class": "backbone",
                    "atom_b_class": "sidechain",
                    "clash_class": "backbone_sidechain",
                },
                {
                    "atom_a": "O",
                    "atom_b": "CG",
                    "distance_angstrom": 0.6103,
                    "atom_a_class": "backbone",
                    "atom_b_class": "sidechain",
                    "clash_class": "backbone_sidechain",
                },
            ],
            "clash_classes": ["backbone_backbone", "backbone_sidechain"],
            "blocking_prepare": True,
        }
    ]


def test_find_inter_residue_heavy_atom_clashes_reports_impossible_cross_residue_contacts(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   ARG H   1      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H   1      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H   1      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H   1      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H   1      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H   1      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H   1      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H   1      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H   1      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H   1      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H   1      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      22.500   8.900   8.300  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      14.950  12.250  11.150  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    clashes = find_inter_residue_heavy_atom_clashes(input_path)
    assert clashes == [
        {
            "chain_id": "A",
            "resseq": 4,
            "icode": "",
            "resname": "ILE",
            "normalized_resname": "ILE",
            "partner_chain_id": "H",
            "partner_resseq": 1,
            "partner_icode": "",
            "partner_resname": "ARG",
            "partner_normalized_resname": "ARG",
            "min_distance_angstrom": 0.0866,
            "clashes": [
                {
                    "atom_a": "CD1",
                    "atom_b": "CZ",
                    "distance_angstrom": 0.0866,
                    "atom_a_class": "sidechain",
                    "atom_b_class": "sidechain",
                    "clash_class": "sidechain_sidechain",
                }
            ],
            "clash_classes": ["sidechain_sidechain"],
            "blocking_prepare": True,
        }
    ]


def test_write_inter_residue_heavy_atom_clash_report_persists_payload(tmp_path: Path) -> None:
    input_path = tmp_path / "input.pdb"
    report_path = tmp_path / "report.json"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   ARG H   1      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H   1      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H   1      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H   1      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H   1      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H   1      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H   1      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H   1      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H   1      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H   1      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H   1      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      22.500   8.900   8.300  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      14.950  12.250  11.150  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    payload = write_inter_residue_heavy_atom_clash_report(input_path, report_path)

    assert report_path.exists()
    assert payload["blocking_inter_residue_heavy_atom_clashes"] is True
    assert len(payload["inter_residue_heavy_atom_clashes"]) == 1
    persisted = report_path.read_text(encoding="utf-8")
    assert '"blocking_inter_residue_heavy_atom_clashes": true' in persisted


def test_write_inter_residue_heavy_atom_clash_report_filters_preexisting_reference_clashes(tmp_path: Path) -> None:
    reference_path = tmp_path / "reference.pdb"
    candidate_path = tmp_path / "candidate.pdb"
    report_path = tmp_path / "report.json"
    pdb_text = (
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM      3  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM      4  O   ILE A   4      22.500   8.900   8.300  1.00 20.00           O",
                "ATOM      5  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM      6  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM      7  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM      8  CD1 ILE A   4      14.950  12.250  11.150  1.00 20.00           C",
                "TER",
                "ATOM      9  N   MET A   7      22.888   9.865  31.392  1.00 20.00           N",
                "ATOM     10  CA  MET A   7      23.269   9.735  32.793  1.00 20.00           C",
                "ATOM     11  C   MET A   7      22.692  10.873  33.629  1.00 20.00           C",
                "ATOM     12  O   MET A   7      23.409  11.521  34.392  1.00 20.00           O",
                "ATOM     13  CB  MET A   7      23.522   8.407  32.066  1.00 20.00           C",
                "ATOM     14  CG  MET A   7      23.226   9.713  10.100  1.00 20.00           C",
                "ATOM     15  SD  MET A   7      23.730   5.292  33.532  1.00 20.00           S",
                "ATOM     16  CE  MET A   7      22.617   6.441  34.657  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    reference_path.write_text(pdb_text, encoding="utf-8")
    candidate_path.write_text(pdb_text, encoding="utf-8")

    payload = write_inter_residue_heavy_atom_clash_report(
        candidate_path,
        report_path,
        reference_path=reference_path,
    )

    assert payload["blocking_inter_residue_heavy_atom_clashes"] is False
    assert payload["inter_residue_heavy_atom_clashes"] == []


def test_partition_sidechain_repairable_clashes_distinguishes_backbone_breaks() -> None:
    repairable, blocking = partition_sidechain_repairable_clashes(
        [
            {
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "resname": "GLU",
                "normalized_resname": "GLU",
                "clashes": [{"atom_a": "O", "atom_b": "CG", "distance_angstrom": 0.8}],
            },
            {
                "chain_id": "H",
                "resseq": 33,
                "icode": "",
                "resname": "GLY",
                "normalized_resname": "GLY",
                "clashes": [{"atom_a": "N", "atom_b": "CA", "distance_angstrom": 0.7}],
            },
            {
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "resname": "ASN",
                "normalized_resname": "ASN",
                "clashes": [{"atom_a": "N", "atom_b": "OD1", "distance_angstrom": 0.7}],
            },
        ]
    )

    assert [issue["resseq"] for issue in repairable] == [32, 34]
    assert [issue["resseq"] for issue in blocking] == [33]


def test_empty_repair_summary_matches_triggered_schema(tmp_path: Path) -> None:
    """The not-triggered placeholder must expose the same keys as a triggered
    repair summary so downstream consumers never branch on missing fields."""
    from abag_rbfe.structure import empty_repair_summary, repair_incomplete_standard_residues_with_pdbfixer

    complete_pdb = tmp_path / "complete.pdb"
    complete_pdb.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLY H   1      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      13.350   9.250   9.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    triggered = repair_incomplete_standard_residues_with_pdbfixer(complete_pdb, complete_pdb)
    assert set(empty_repair_summary()) == set(triggered)
    assert empty_repair_summary()["trigger_residues"] == []


def test_strip_hydrogen_atoms_removes_hydrogens_with_provenance(tmp_path: Path) -> None:
    from abag_rbfe.structure import strip_hydrogen_atoms

    input_path = tmp_path / "input.pdb"
    output_path = tmp_path / "output.pdb"
    input_path.write_text(
        "\n".join(
            [
                "HEADER    DEMO",
                "ATOM      1  N   GLU H   1      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  H   GLU H   1      10.500  10.200   8.800  1.00 20.00           H",
                "ATOM      3  CA  GLU H   1      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      4  HA  GLU H   1      12.500   9.600   9.600  1.00 20.00           H",
                "ATOM      5  C   GLU H   1      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      6  O   GLU H   1      13.350   9.250   9.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    summary = strip_hydrogen_atoms(input_path, output_path)

    assert summary["stripped_hydrogen_count"] == 2
    assert [a["atom_name"] for a in summary["stripped_hydrogen_atoms"]] == ["H", "HA"]
    remaining = [l for l in output_path.read_text().splitlines() if l.startswith("ATOM")]
    assert len(remaining) == 4
    assert all(not l[12:16].strip().startswith("H") for l in remaining)
