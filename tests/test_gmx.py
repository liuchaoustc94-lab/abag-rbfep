from pathlib import Path

from abag_rbfe.gmx import (
    deduplicate_standard_topology_includes,
    ensure_local_gmxlib,
    generate_hybrid_topology,
    gro_file_is_valid,
    inspect_gro_file,
    materialize_staged_equilibration_restraints,
)


def _minimal_topology_text(
    *,
    molecule: str = "Protein_chain_A",
    residue: str = "ALA",
    atoms: tuple[str, ...] = ("N", "CA"),
    posre_include: str | None = None,
) -> str:
    lines = [
        "[ moleculetype ]",
        f"{molecule} 3",
        "",
        "[ atoms ]",
        "; nr type resnr residue atom cgnr charge mass",
    ]
    for index, atom_name in enumerate(atoms, start=1):
        lines.append(f"{index:5d}   CT      1    {residue:<3}   {atom_name:<4}   {index:5d}   0.000   12.0100")
    if len(atoms) >= 2:
        lines.extend(["", "[ bonds ]", "1 2"])
    if posre_include:
        lines.extend(["", "#ifdef POSRES", f'#include "{posre_include}"', "#endif"])
    return "\n".join(lines) + "\n"


def test_local_gmxlib_overlay_adds_hybrid_residues_to_residuetypes(tmp_path: Path) -> None:
    gmx_top = tmp_path / "gmx_top"
    gmx_top.mkdir()
    (gmx_top / "residuetypes.dat").write_text("ALA\tProtein\nCYS\tProtein\n", encoding="utf-8")
    (gmx_top / "specbond.dat").write_text("", encoding="utf-8")
    (gmx_top / "spc216.gro").write_text("water\n", encoding="utf-8")

    pmx_mutff_root = tmp_path / "mutff"
    ff_dir = pmx_mutff_root / "amber99sb-star-ildn-mut.ff"
    ff_dir.mkdir(parents=True)
    (pmx_mutff_root / "specbond.dat").write_text("pmx-specbond\n", encoding="utf-8")
    (ff_dir / "forcefield.itp").write_text("; ff\n", encoding="utf-8")
    (ff_dir / "mutres.rtp").write_text(
        "\n".join(
            [
                "[ bondedtypes ]",
                "[ C2A ] ; Cys -> Ala",
                "[ atoms ]",
                "  CA  C",
                "[ A2V ] ; Ala -> Val",
                "[ atoms ]",
                "  CA  C",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    gmxlib_dir = ensure_local_gmxlib(
        job_dir=tmp_path / "job",
        gmx_top_dir=gmx_top,
        pmx_mutff_root=pmx_mutff_root,
        force_field="amber99sb-star-ildn-mut",
    )

    residuetypes = (gmxlib_dir / "residuetypes.dat").read_text(encoding="utf-8")
    assert "C2A\tProtein" in residuetypes
    assert "A2V\tProtein" in residuetypes
    assert (gmxlib_dir / "specbond.dat").resolve() == (pmx_mutff_root / "specbond.dat").resolve()


def test_deduplicate_standard_topology_includes_keeps_first_water_include(tmp_path: Path) -> None:
    topology = tmp_path / "system.top"
    topology.write_text(
        "\n".join(
            [
                '#include "amber99sb-star-ildn-mut.ff/forcefield.itp"',
                "[ moleculetype ]",
                "Protein_chain_A 3",
                '#include "/tmp/gmxlib/amber99sb-star-ildn-mut.ff/tip3p.itp"',
                "#ifdef POSRES_WATER",
                "[ position_restraints ]",
                "1 1 1000 1000 1000",
                "#endif",
                '#include "/tmp/gmxlib/amber99sb-star-ildn-mut.ff/tip3p.itp"',
                '#include "/tmp/gmxlib/amber99sb-star-ildn-mut.ff/ions.itp"',
                "",
                "[ molecules ]",
                "Protein_chain_A 1",
                "SOL 100",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    removed = deduplicate_standard_topology_includes(topology)

    output = topology.read_text(encoding="utf-8")
    assert removed == {"tip3p.itp": 1}
    assert output.count("tip3p.itp") == 1
    assert output.count("ions.itp") == 1
    assert "#ifdef POSRES_WATER" in output
    assert "SOL 100" in output


def test_generate_hybrid_topology_converts_only_mutated_chain_itps(tmp_path: Path, monkeypatch) -> None:
    topology = tmp_path / "topol.top"
    topology.write_text(
        "\n".join(
            [
                '#include "amber99sb-star-ildn-mut.ff/forcefield.itp"',
                '#include "topol_Protein_chain_A.itp"',
                '#include "topol_Protein_chain_H.itp"',
                "",
                "[ system ]",
                "PMX MODEL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_minimal_topology_text(), encoding="utf-8")
    (tmp_path / "topol_Protein_chain_H.itp").write_text(
        _minimal_topology_text(molecule="Protein_chain_H"),
        encoding="utf-8",
    )

    calls: list[list[str]] = []

    def fake_run(command: list[str], check: bool, cwd: str) -> None:
        assert check is True
        calls.append(command)
        output_path = Path(command[command.index("-o") + 1])
        if not output_path.is_absolute():
            output_path = Path(cwd) / output_path
        output_path.write_text(_minimal_topology_text(), encoding="utf-8")

    monkeypatch.setattr("abag_rbfe.gmx.subprocess.run", fake_run)

    summary = generate_hybrid_topology(
        topology,
        tmp_path / "pmxtop.top",
        "amber99sb-star-ildn-mut",
        ["A"],
        ["pmx"],
    )

    output = (tmp_path / "pmxtop.top").read_text(encoding="utf-8")
    assert summary["mode"] == "per_chain"
    assert len(calls) == 1
    assert calls[0] == [
        "pmx",
        "gentop",
        "-p",
        str((tmp_path / "topol_Protein_chain_A.itp").resolve()),
        "-o",
        str((tmp_path / "pmx_topol_Protein_chain_A.itp").resolve()),
        "-ff",
        "amber99sb-star-ildn-mut",
        "--norecursive",
    ]
    assert '#include "pmx_topol_Protein_chain_A.itp"' in output
    assert '#include "topol_Protein_chain_H.itp"' in output
    assert (tmp_path / "pmx_topol_Protein_chain_A.itp").exists()


def test_generate_hybrid_topology_falls_back_to_recursive_topology_conversion(tmp_path: Path, monkeypatch) -> None:
    topology = tmp_path / "topol.top"
    topology.write_text(_minimal_topology_text(), encoding="utf-8")

    calls: list[list[str]] = []

    def fake_run(command: list[str], check: bool, cwd: str) -> None:
        assert check is True
        calls.append(command)
        output_path = Path(command[command.index("-o") + 1])
        if not output_path.is_absolute():
            output_path = Path(cwd) / output_path
        output_path.write_text(topology.read_text(encoding="utf-8"), encoding="utf-8")

    monkeypatch.setattr("abag_rbfe.gmx.subprocess.run", fake_run)

    summary = generate_hybrid_topology(
        topology,
        tmp_path / "pmxtop.top",
        "amber99sb-star-ildn-mut",
        ["A"],
        ["pmx"],
    )

    assert summary["mode"] == "recursive_fallback"
    assert summary["fallback_reason"] == "no_mutated_chain_itps_detected"
    assert calls == [
        [
            "pmx",
            "gentop",
            "-p",
            str(topology),
            "-o",
            str(tmp_path / "pmxtop.top"),
            "-ff",
            "amber99sb-star-ildn-mut",
        ]
    ]


def test_generate_hybrid_topology_reuses_existing_mutated_itps_when_repairs_are_noops(tmp_path: Path, monkeypatch) -> None:
    topology = tmp_path / "topol.top"
    topology.write_text(
        "\n".join(
            [
                '#include "amber99sb-star-ildn-mut.ff/forcefield.itp"',
                '#include "topol_Protein_chain_A.itp"',
                '#include "topol_Protein_chain_H.itp"',
                "",
                "[ system ]",
                "PMX MODEL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_minimal_topology_text(), encoding="utf-8")
    (tmp_path / "topol_Protein_chain_H.itp").write_text(
        _minimal_topology_text(molecule="Protein_chain_H"),
        encoding="utf-8",
    )
    (tmp_path / "pmx_topol_Protein_chain_A.itp").write_text(_minimal_topology_text(), encoding="utf-8")
    restore_summary = tmp_path / "mutant_standard_residue_repair.json"
    restore_summary.write_text(
        '{"attempted_residue_count": 0, "restored_residue_count": 0, "unresolved_residue_count": 0}\n',
        encoding="utf-8",
    )
    pdbfixer_summary = tmp_path / "mutant_pdbfixer_repair.json"
    pdbfixer_summary.write_text(
        '{"trigger_residue_count": 0, "blocking_residue_count": 0, "remaining_incomplete_standard_residue_count": 0}\n',
        encoding="utf-8",
    )

    def fail_run(*args, **kwargs) -> None:  # pragma: no cover - should never be called
        raise AssertionError("subprocess.run should not be called when reusing existing pmx outputs")

    monkeypatch.setattr("abag_rbfe.gmx.subprocess.run", fail_run)

    summary = generate_hybrid_topology(
        topology,
        tmp_path / "pmxtop.top",
        "amber99sb-star-ildn-mut",
        ["A"],
        ["pmx"],
        restore_summary_path=restore_summary,
        pdbfixer_summary_path=pdbfixer_summary,
        allow_reuse_existing=True,
    )

    output = (tmp_path / "pmxtop.top").read_text(encoding="utf-8")
    assert summary["mode"] == "reused_existing"
    assert '#include "pmx_topol_Protein_chain_A.itp"' in output
    assert '#include "topol_Protein_chain_H.itp"' in output


def test_generate_hybrid_topology_reuses_existing_monolithic_pmxtop_when_repairs_are_noops(
    tmp_path: Path, monkeypatch
) -> None:
    topology = tmp_path / "topol.top"
    topology.write_text(_minimal_topology_text(), encoding="utf-8")
    pmxtop = tmp_path / "pmxtop.top"
    pmxtop.write_text(_minimal_topology_text(), encoding="utf-8")
    restore_summary = tmp_path / "mutant_standard_residue_repair.json"
    restore_summary.write_text(
        '{"attempted_residue_count": 0, "restored_residue_count": 0, "unresolved_residue_count": 0}\n',
        encoding="utf-8",
    )
    pdbfixer_summary = tmp_path / "mutant_pdbfixer_repair.json"
    pdbfixer_summary.write_text(
        '{"trigger_residue_count": 0, "blocking_residue_count": 0, "remaining_incomplete_standard_residue_count": 0}\n',
        encoding="utf-8",
    )

    def fail_run(*args, **kwargs) -> None:  # pragma: no cover - should never be called
        raise AssertionError("subprocess.run should not be called when reusing monolithic pmx outputs")

    monkeypatch.setattr("abag_rbfe.gmx.subprocess.run", fail_run)

    summary = generate_hybrid_topology(
        topology,
        pmxtop,
        "amber99sb-star-ildn-mut",
        ["A"],
        ["pmx"],
        restore_summary_path=restore_summary,
        pdbfixer_summary_path=pdbfixer_summary,
        allow_reuse_existing=True,
    )

    assert summary["mode"] == "reused_existing"
    assert pmxtop.read_text(encoding="utf-8") == _minimal_topology_text()


def test_generate_hybrid_topology_regenerates_stale_existing_mutated_itps_even_when_repairs_are_noops(
    tmp_path: Path, monkeypatch
) -> None:
    topology = tmp_path / "topol.top"
    topology.write_text(
        "\n".join(
            [
                '#include "amber99sb-star-ildn-mut.ff/forcefield.itp"',
                '#include "topol_Protein_chain_A.itp"',
                "",
                "[ system ]",
                "PMX MODEL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    source_itp = _minimal_topology_text(atoms=("N", "CA", "CB"))
    stale_itp = _minimal_topology_text(atoms=("N", "CA"))
    (tmp_path / "topol_Protein_chain_A.itp").write_text(source_itp, encoding="utf-8")
    (tmp_path / "pmx_topol_Protein_chain_A.itp").write_text(stale_itp, encoding="utf-8")
    restore_summary = tmp_path / "mutant_standard_residue_repair.json"
    restore_summary.write_text(
        '{"attempted_residue_count": 0, "restored_residue_count": 0, "unresolved_residue_count": 0}\n',
        encoding="utf-8",
    )
    pdbfixer_summary = tmp_path / "mutant_pdbfixer_repair.json"
    pdbfixer_summary.write_text(
        '{"trigger_residue_count": 0, "blocking_residue_count": 0, "remaining_incomplete_standard_residue_count": 0}\n',
        encoding="utf-8",
    )

    calls: list[list[str]] = []

    def fake_run(command: list[str], check: bool, cwd: str) -> None:
        assert check is True
        calls.append(command)
        output_path = Path(command[command.index("-o") + 1])
        if not output_path.is_absolute():
            output_path = Path(cwd) / output_path
        output_path.write_text(source_itp, encoding="utf-8")

    monkeypatch.setattr("abag_rbfe.gmx.subprocess.run", fake_run)

    summary = generate_hybrid_topology(
        topology,
        tmp_path / "pmxtop.top",
        "amber99sb-star-ildn-mut",
        ["A"],
        ["pmx"],
        restore_summary_path=restore_summary,
        pdbfixer_summary_path=pdbfixer_summary,
        allow_reuse_existing=True,
    )

    assert summary["mode"] == "per_chain"
    assert len(calls) == 1
    assert (tmp_path / "pmx_topol_Protein_chain_A.itp").read_text(encoding="utf-8") == source_itp


def test_generate_hybrid_topology_regenerates_stale_existing_monolithic_pmxtop_even_when_repairs_are_noops(
    tmp_path: Path, monkeypatch
) -> None:
    topology = tmp_path / "topol.top"
    source_topology = _minimal_topology_text(atoms=("N", "CA", "CB"))
    topology.write_text(source_topology, encoding="utf-8")
    pmxtop = tmp_path / "pmxtop.top"
    pmxtop.write_text(_minimal_topology_text(atoms=("N", "CA")), encoding="utf-8")
    restore_summary = tmp_path / "mutant_standard_residue_repair.json"
    restore_summary.write_text(
        '{"attempted_residue_count": 0, "restored_residue_count": 0, "unresolved_residue_count": 0}\n',
        encoding="utf-8",
    )
    pdbfixer_summary = tmp_path / "mutant_pdbfixer_repair.json"
    pdbfixer_summary.write_text(
        '{"trigger_residue_count": 0, "blocking_residue_count": 0, "remaining_incomplete_standard_residue_count": 0}\n',
        encoding="utf-8",
    )

    calls: list[list[str]] = []

    def fake_run(command: list[str], check: bool, cwd: str) -> None:
        assert check is True
        calls.append(command)
        output_path = Path(command[command.index("-o") + 1])
        if not output_path.is_absolute():
            output_path = Path(cwd) / output_path
        output_path.write_text(source_topology, encoding="utf-8")

    monkeypatch.setattr("abag_rbfe.gmx.subprocess.run", fake_run)

    summary = generate_hybrid_topology(
        topology,
        pmxtop,
        "amber99sb-star-ildn-mut",
        ["A"],
        ["pmx"],
        restore_summary_path=restore_summary,
        pdbfixer_summary_path=pdbfixer_summary,
        allow_reuse_existing=True,
    )

    assert summary["mode"] == "recursive_fallback"
    assert calls == [
        [
            "pmx",
            "gentop",
            "-p",
            str(topology),
            "-o",
            str(pmxtop),
            "-ff",
            "amber99sb-star-ildn-mut",
        ]
    ]
    assert pmxtop.read_text(encoding="utf-8") == source_topology


def test_materialize_staged_equilibration_restraints_patches_chain_itps_and_writes_stage_files(tmp_path: Path) -> None:
    topology = tmp_path / "system.top"
    topology.write_text(
        "\n".join(
            [
                '#include "amber99sb-star-ildn-mut.ff/forcefield.itp"',
                '#include "topol_Protein_chain_A.itp"',
                "",
                "[ system ]",
                "PMX MODEL",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    chain_itp = tmp_path / "topol_Protein_chain_A.itp"
    chain_itp.write_text(
        _minimal_topology_text(
            atoms=("N", "CA", "CB", "H1"),
            posre_include="posre_Protein_chain_A.itp",
        ),
        encoding="utf-8",
    )

    summary = materialize_staged_equilibration_restraints(
        topology,
        heavy_force_constant=900.0,
        backbone_force_constant=150.0,
    )

    chain_text = chain_itp.read_text(encoding="utf-8")
    heavy_text = (tmp_path / "posre_Protein_chain_A_stage_heavy.itp").read_text(encoding="utf-8")
    backbone_text = (tmp_path / "posre_Protein_chain_A_stage_backbone.itp").read_text(encoding="utf-8")

    assert summary["status"] == "ok"
    assert summary["modified_topology_count"] == 1
    assert '#ifdef POSRES_STAGE_HEAVY' in chain_text
    assert '#include "posre_Protein_chain_A_stage_heavy.itp"' in chain_text
    assert '#include "posre_Protein_chain_A_stage_backbone.itp"' in chain_text
    assert "     3      1     900     900     900" in heavy_text
    assert "     2      1     150     150     150" in backbone_text
    assert "     3      1     150     150     150" not in backbone_text


def test_materialize_staged_equilibration_restraints_can_patch_monolithic_topology_file(tmp_path: Path) -> None:
    topology = tmp_path / "system.top"
    topology.write_text(
        _minimal_topology_text(
            molecule="Protein",
            atoms=("N", "CA", "CB"),
            posre_include="posre.itp",
        ),
        encoding="utf-8",
    )

    summary = materialize_staged_equilibration_restraints(topology)

    topology_text = topology.read_text(encoding="utf-8")
    assert summary["status"] == "ok"
    assert '#include "posre_stage_heavy.itp"' in topology_text
    assert '#include "posre_stage_backbone.itp"' in topology_text
    assert (tmp_path / "posre_stage_heavy.itp").is_file()
    assert (tmp_path / "posre_stage_backbone.itp").is_file()


def test_inspect_gro_file_accepts_minimal_valid_gro(tmp_path: Path) -> None:
    gro_path = tmp_path / "valid.gro"
    gro_path.write_text(
        "Mock GRO\n"
        "1\n"
        "    1ALA      N    1   0.000   0.000   0.000\n"
        "   1.00000   1.00000   1.00000\n",
        encoding="utf-8",
    )

    summary = inspect_gro_file(gro_path)

    assert summary["valid"] is True
    assert summary["reason"] == "ok"
    assert gro_file_is_valid(gro_path) is True


def test_inspect_gro_file_rejects_malformed_coordinate_field(tmp_path: Path) -> None:
    gro_path = tmp_path / "invalid_coords.gro"
    gro_path.write_text(
        "Mock GRO\n"
        "1\n"
        "    1ALA      N    1   0.000-409203.000-409203.000\n"
        "   1.00000   1.00000   1.00000\n",
        encoding="utf-8",
    )

    summary = inspect_gro_file(gro_path)

    assert summary["valid"] is False
    assert summary["reason"] == "invalid_coordinate"
    assert summary["line_number"] == 3
    assert gro_file_is_valid(gro_path) is False


def test_inspect_gro_file_rejects_absurd_coordinates_and_box(tmp_path: Path) -> None:
    gro_path = tmp_path / "invalid_range.gro"
    gro_path.write_text(
        "Mock GRO\n"
        "1\n"
        "    1ALA      N    1 200.000   0.000   0.000\n"
        "2000.00000   1.00000   1.00000\n",
        encoding="utf-8",
    )

    summary = inspect_gro_file(gro_path)

    assert summary["valid"] is False
    assert summary["reason"] == "coordinate_out_of_range"
    assert summary["line_number"] == 3


def test_inspect_gro_file_rejects_isolated_residue_hydrogen(tmp_path: Path) -> None:
    gro_path = tmp_path / "isolated_hydrogen.gro"
    gro_path.write_text(
        "Mock GRO\n"
        "5\n"
        "    1ALA      N    1   0.000   0.000   0.000\n"
        "    1ALA     CA    2   0.145   0.000   0.000\n"
        "    1ALA      C    3   0.290   0.000   0.000\n"
        "    1ALA      O    4   0.410   0.000   0.000\n"
        "    1ALA     HA    5   5.000   5.000   5.000\n"
        "   6.00000   6.00000   6.00000\n",
        encoding="utf-8",
    )

    summary = inspect_gro_file(gro_path)

    assert summary["valid"] is False
    assert summary["reason"] == "isolated_residue_hydrogen"
    assert summary["line_number"] == 7
    assert summary["residue_number"] == 1
    assert summary["residue_name"] == "ALA"
    assert summary["atom_name"] == "HA"
    assert summary["nearest_heavy_atom"] == "O"
    assert summary["nearest_heavy_distance_nm"] is not None
    assert float(summary["nearest_heavy_distance_nm"]) > 8.0
    assert gro_file_is_valid(gro_path) is False


def test_validate_hybrid_topology_integrity_detects_complete_and_broken(tmp_path: Path) -> None:
    from abag_rbfe.gmx import validate_hybrid_topology_integrity

    def write_case(root: Path, *, with_hydrogens: bool) -> None:
        atoms_pdb = [
            ("N", "N"), ("H", "H"), ("CA", "C"), ("HA", "H"), ("CB", "C"),
            ("HB1", "H"), ("HB2", "H"), ("CG", "C"), ("HG1", "H"), ("HG2", "H"),
            ("CD", "C"), ("OE1", "O"), ("NE2", "N"), ("HE21", "H"), ("HE22", "H"),
            ("C", "C"), ("O", "O"), ("HV", "H"),
        ]
        charges = {
            "N": ("N", -0.4157), "H": ("H", 0.2719), "CA": ("CT", -0.0031),
            "HA": ("H1", 0.0850), "CB": ("CT", -0.0036), "HB1": ("HC", 0.0171),
            "HB2": ("HC", 0.0171), "CG": ("CT", -0.0645), "HG1": ("HC", 0.0352),
            "HG2": ("HC", 0.0352), "CD": ("C", 0.6951), "OE1": ("O", -0.6086),
            "NE2": ("N", -0.9407), "HE21": ("H", 0.4251), "HE22": ("H", 0.4251),
            "C": ("C", 0.5973), "O": ("O", -0.5679), "HV": ("DUM_HC", 0.0),
        }
        pdb_lines = []
        for idx, (name, _elem) in enumerate(atoms_pdb, start=1):
            pdb_lines.append(
                f"ATOM  {idx:5d} {name:<4s} Q2A W  89      {1.0+idx*0.1:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00 10.00          {_elem:>2s}"
            )
        (root / "mutant.pdb").write_text("\n".join(pdb_lines) + "\nEND\n", encoding="utf-8")

        itp_atoms = atoms_pdb if with_hydrogens else [a for a in atoms_pdb if a[1] != "H"]
        lines = ["[ atoms ]"]
        for nr, (name, _e) in enumerate(itp_atoms, start=1):
            atype, charge = charges[name]
            lines.append(f"{nr:6d} {atype:>6s}     89    Q2A {name:>6s} {nr:6d} {charge:10.4f}    12.0100")
        (root / "pmx_topol_Protein_chain_W.itp").write_text("\n".join(lines) + "\n", encoding="utf-8")

    good = tmp_path / "good"
    good.mkdir()
    write_case(good, with_hydrogens=True)
    good_result = validate_hybrid_topology_integrity(good)
    assert good_result["checked"] is True
    assert good_result["ok"] is True
    assert good_result["issues"] == []
    assert good_result["residues"][0]["atom_count"] == 18
    assert good_result["residues"][0]["charge_a"] == 0.0

    broken = tmp_path / "broken"
    broken.mkdir()
    write_case(broken, with_hydrogens=False)
    broken_result = validate_hybrid_topology_integrity(broken)
    assert broken_result["checked"] is True
    assert broken_result["ok"] is False
    assert any("atoms" in issue for issue in broken_result["issues"])
    assert any("state A charge" in issue for issue in broken_result["issues"])


def test_swap_hybrid_residue_ab_states_roundtrip_and_charge_invariance(tmp_path: Path) -> None:
    from abag_rbfe.gmx import swap_hybrid_residue_ab_states, validate_dssb_charge_invariance

    itp = tmp_path / "pmx_topol_Protein_chain_U.itp"
    itp.write_text(
        "\n".join(
            [
                "[ atoms ]",
                "  1  N   89  Q2A  N   1  -0.415700  14.0100",
                "  2  CT  89  Q2A  CA  2  -0.003100  12.0100  CT  0.033700  12.0100",
                "  3  CT  89  Q2A  CB  3  -0.003600  12.0100  DUM_CT  0.000000  3.9633",
                "  4  C   89  Q2A  CD  4   0.695100  12.0100  DUM_C   0.000000  3.9633",
                "[ bonds ]",
                "  2  3  1  0.109000  284512.0  0.110000  250000.0",
                "[ angles ]",
                "  1  2  3  1  109.500000  418.400000  110.000000  400.000000 ; N CA CB",
                "[ dihedrals ]",
                "  1  2  3  4  9  0  8.368  1  0  1.12968  2 ; comment",
                "",
            ]
        ),
        encoding="utf-8",
    )
    summary = swap_hybrid_residue_ab_states(itp, set())
    assert summary["swapped_atoms"] == 3  # CA/CB/CD 三个带 B 态列的原子被交换（N 无 B 列保持不变）

    lines = itp.read_text().splitlines()
    atoms = [l.split() for l in lines if l.strip().startswith(("2", "3", "4"))]
    # CA: 交换后 type=CT, chargeA=0.0337, massA=12.01; chargeB=-0.0031
    ca = [l for l in atoms if l[0] == "2"][0]
    assert ca[6] == "0.033700" and ca[9] == "-0.003100"
    # CB: A 态变为 DUM（原 B），B 态变为实原子
    cb = [l for l in atoms if l[0] == "3"][0]
    assert cb[1] == "DUM_CT" and cb[8] == "CT" and cb[6] == "0.000000" and cb[9] == "-0.003600"
    # bonds: b0/kb 对换
    bond = [l.split() for l in lines if l.strip().startswith("2  3")][0]
    assert bond[3] == "0.110000" and bond[5] == "0.109000"
    # angles: theta/k 对换
    ang = [l.split() for l in lines if l.strip().startswith("1  2  3")][0]
    assert ang[4] == "110.000000" and ang[6] == "109.500000"
    # dihedrals func 9: 前后 6 参数对换
    dih = [l.split() for l in lines if l.strip().startswith("1  2  3  4  9")][0]
    assert dih[5:8] == ["0", "1.12968", "2"] and dih[8:11] == ["0", "8.368", "1"]

    # DSSB 电荷不变性：swap 前后 chargeA_sum 应等于 swap 前 chargeB_sum
    # 原拓扑 chargeA = -0.4157-0.0031-0.0036+0.6951 = 0.2727；chargeB = -0.4157+0.0337+0+0 = -0.382
    # swap 后 chargeA == 原 chargeB
    result = validate_dssb_charge_invariance([itp])
    assert result["total_charge_state_a"] == -0.382
    assert result["ok"] is False  # 单独 unbound 腿不满足不变性（需要 bound 腿配对）


def test_rename_pdb_chains_and_mapping_suggestion(tmp_path: Path) -> None:
    from abag_rbfe.structure import rename_pdb_chains, suggest_unbound_chain_mapping

    pdb = tmp_path / "in.pdb"
    pdb.write_text(
        "ATOM      1  N   GLY H   1      11.000  10.000   8.000  1.00 20.00           N\n"
        "ATOM      2  CA  GLY H   1      12.200  10.100   8.700  1.00 20.00           C\n"
        "ATOM      3  N   GLY L   2      15.000  10.000   8.000  1.00 20.00           N\n"
        "END\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.pdb"
    summary = rename_pdb_chains(pdb, out, {"H": "U", "L": "V"})
    assert summary["renamed_atom_count"] == 3
    text = out.read_text()
    assert "GLY U   1" in text and "GLY V   2" in text and "GLY H" not in text

    mapping = suggest_unbound_chain_mapping(["H", "L", "Y"], ["H", "L"])
    assert mapping == {"H": "A", "L": "B"}  # A/B 是未被占用的前两个字母


def test_outlier_classification_and_correction(tmp_path: Path) -> None:
    from abag_rbfe.outliers import (
        apply_outlier_correction,
        classify_outlier_features,
        fit_shrinkage_alpha,
        mutation_size_change,
        outlier_class_label,
    )

    assert mutation_size_change("Y", "A") == -7  # Tyr 8 - Ala 1
    assert mutation_size_change("G", "V") == 3
    feats = classify_outlier_features(wt="Y", mut="A")
    assert feats["is_aromatic_deletion"] and feats["is_large_deletion"]
    assert feats["is_buried_aromatic_large_deletion"]
    assert outlier_class_label(feats) == "buried_aromatic_large_deletion"

    feats2 = classify_outlier_features(wt="D", mut="K")
    assert outlier_class_label(feats2) == "charge_changing"

    feats3 = classify_outlier_features(wt="S", mut="T")
    assert outlier_class_label(feats3) == ""

    alpha = fit_shrinkage_alpha([10.0, 12.0], [5.0, 6.0], [True, True])
    assert 0.3 < alpha < 0.7
    corrected = apply_outlier_correction([10.0, 8.0], [True, False], 0.5)
    assert corrected == [5.0, 8.0]
