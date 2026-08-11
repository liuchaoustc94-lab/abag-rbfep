from pathlib import Path

from pmx.forcefield import TopolBase


def _write_minimal_itp(path: Path) -> None:
    path.write_text(
        "\n".join(
            [
                "[ moleculetype ]",
                "protein 3",
                "",
                "[ atoms ]",
                "1 amber99_0 52A PRO N 1 -0.300000 14.0100",
                "2 amber99_1 52A PRO CA 1 0.100000 12.0100",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_topolbase_reads_insertion_coded_topology_residue_labels(tmp_path: Path) -> None:
    itp_path = tmp_path / "insertion_code.itp"
    _write_minimal_itp(itp_path)

    topology = TopolBase(str(itp_path))

    assert [atom.resnr for atom in topology.atoms] == ["52A", "52A"]
    assert len(topology.residues) == 1
    assert topology.residues[0].id == "52A"
    assert all(atom.molecule is topology.residues[0] for atom in topology.atoms)


def test_topolbase_write_preserves_insertion_coded_topology_residue_labels(tmp_path: Path) -> None:
    itp_path = tmp_path / "insertion_code.itp"
    out_path = tmp_path / "written.itp"
    _write_minimal_itp(itp_path)

    topology = TopolBase(str(itp_path))
    topology.write(str(out_path), stateBonded="A", stateTypes="A", stateQ="A")

    atom_lines = [line for line in out_path.read_text(encoding="utf-8").splitlines() if line.strip().startswith(("1", "2"))]

    assert atom_lines
    assert atom_lines[0].split()[2] == "52A"
    assert atom_lines[1].split()[2] == "52A"
