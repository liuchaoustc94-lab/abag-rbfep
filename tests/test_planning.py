import json
from pathlib import Path
import subprocess
import time

from abag_rbfe.execution import CommandOutcome, CommandRunner
from abag_rbfe.execution import discover_visible_gpu_devices
from abag_rbfe.planning import build_batch_plan, hydrate_protocol_config
from abag_rbfe.reporting import summarize_job
from abag_rbfe.stages import resume_job, run_job


def _write_valid_mock_gro(path: Path) -> None:
    path.write_text(
        "Mock GRO\n"
        "1\n"
        "    1ALA      N    1   0.000   0.000   0.000\n"
        "   1.00000   1.00000   1.00000\n",
        encoding="utf-8",
    )


def _read_mdp_parameter(path: Path, key: str) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" not in line:
            continue
        lhs, rhs = line.split("=", 1)
        if lhs.strip() == key:
            return line.split("=", 1)[1].strip()
    raise AssertionError(f"{key} not found in {path}")


def _build_seed_scope_demo_job(tmp_path: Path, *, batch_id: str) -> Path:
    case_dir = tmp_path / batch_id
    case_dir.mkdir(parents=True, exist_ok=True)

    system_path = case_dir / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: seed_scope_demo",
                f"input_structure: {case_dir / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    (case_dir / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED SCOPE DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     14  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM     15  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM     16  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "ATOM     17  CB  ILE A  10      21.500  11.900   9.100  1.00 20.00           C",
                "ATOM     18  CG1 ILE A  10      22.800  11.700   9.900  1.00 20.00           C",
                "ATOM     19  CG2 ILE A  10      20.400  12.600   9.900  1.00 20.00           C",
                "ATOM     20  CD1 ILE A  10      23.300  13.000  10.500  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = case_dir / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = case_dir / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id=batch_id,
        runs_root=tmp_path / "runs",
    )
    return Path(batch_plan.jobs[0].workdir)


def test_batch_plan_and_stage_generation(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "input.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
                "double_h_pair,H,52,,S,T,antibody",
                "double_h_pair,H,54,,N,Q,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 6",
                "repeats: 1",
                "nvt_ps: 100",
                "npt_ps: 500",
                "equilibration_restraint_schedule: staged_backbone_release",
                "equilibration_release_npt_ps: 250",
                "production_ps: 1000",
                "equilibrate_em_steps: 1000",
                "window_relax_em_steps: 250",
                "window_relax_md_ps: 0.2",
                "window_relax_md_dt_ps: 0.00025",
                "temperature_k: 310.0",
                "pressure_bar: 1.0",
                "allow_external_execute: false",
                "grompp_maxwarn_genion: 3",
                "grompp_maxwarn_equilibration: 2",
                "grompp_maxwarn_sampling: 4",
                "equilibration_pressure_coupling: Parrinello-Rahman",
                "equilibration_pressure_tau_ps: 2.0",
                "equilibration_refcoord_scaling: com",
                "sampling_pressure_coupling: Parrinello-Rahman",
                "sampling_pressure_tau_ps: 2.0",
                "sampling_refcoord_scaling: com",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(system_path, mutations_path, protocol_path, batch_id="unit_batch", runs_root=tmp_path / "runs")
    assert len(batch_plan.jobs) == 2
    assert batch_plan.jobs[0].protocol.preset == "single_point"
    assert batch_plan.jobs[1].protocol.preset == "double_point"
    assert batch_plan.jobs[1].protocol.lambda_windows == 6
    assert batch_plan.jobs[1].protocol.repeats == 1
    assert batch_plan.jobs[1].protocol.nvt_ps == 100
    assert batch_plan.jobs[1].protocol.npt_ps == 500
    assert batch_plan.jobs[1].protocol.production_ps == 1000
    assert batch_plan.jobs[1].protocol.production_dt_ps == 0.002
    assert batch_plan.jobs[1].protocol.equilibrate_em_steps == 1000
    assert batch_plan.jobs[1].protocol.window_relax_em_steps == 250
    assert batch_plan.jobs[1].protocol.window_relax_md_ps == 0.2
    assert batch_plan.jobs[1].protocol.window_relax_md_dt_ps == 0.00025
    assert batch_plan.jobs[1].protocol.equilibration_restraint_schedule == "staged_backbone_release"
    assert batch_plan.jobs[1].protocol.equilibration_release_npt_ps == 250
    assert batch_plan.jobs[1].protocol.nonbonded_cutoff_nm == 1.25
    assert batch_plan.jobs[1].protocol.vdw_switch_nm == 1.0
    assert batch_plan.jobs[1].protocol.grompp_maxwarn_genion == 3
    assert batch_plan.jobs[1].protocol.grompp_maxwarn_equilibration == 2
    assert batch_plan.jobs[1].protocol.grompp_maxwarn_sampling == 4

    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False)
    assert [item.stage for item in statuses] == [
        "ingest",
        "prepare",
        "mutate",
        "build_legs",
        "equilibrate",
        "sample",
        "bar",
        "qc",
        "report",
    ]
    summary = summarize_job(job_dir)
    assert summary["stage_count"] == 9
    assert summary["stages"][-1]["stage"] == "report"
    complex_input = (job_dir / "legs" / "complex" / "input.pdb").read_text(encoding="utf-8")
    apo_input = (job_dir / "legs" / "apo" / "input.pdb").read_text(encoding="utf-8")
    pre_relax_mdp = (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "pre_relax.mdp").read_text(encoding="utf-8")
    pre_md_mdp = (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "pre_md.mdp").read_text(encoding="utf-8")
    production_mdp = (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "production.mdp").read_text(encoding="utf-8")
    em_mdp = (job_dir / "legs" / "complex" / "rep01" / "mdp" / "em.mdp").read_text(encoding="utf-8")
    nvt_mdp = (job_dir / "legs" / "complex" / "rep01" / "mdp" / "nvt.mdp").read_text(encoding="utf-8")
    npt_mdp = (job_dir / "legs" / "complex" / "rep01" / "mdp" / "npt.mdp").read_text(encoding="utf-8")
    npt_release_mdp = (job_dir / "legs" / "complex" / "rep01" / "mdp" / "npt_release.mdp").read_text(encoding="utf-8")
    mutate_script = (job_dir / "artifacts" / "commands" / "mutate.sh").read_text(encoding="utf-8")
    equilibrate_script = (job_dir / "artifacts" / "commands" / "equilibrate.sh").read_text(encoding="utf-8")
    sample_script = (job_dir / "artifacts" / "commands" / "sample.sh").read_text(encoding="utf-8")
    assert "TYR H  32" in complex_input
    assert "LYS A  58" in complex_input
    assert "LYS A  58" not in apo_input
    assert "define                  = -DFLEXIBLE" in pre_relax_mdp
    assert "init-lambda-state       = 0" in pre_relax_mdp
    assert "nsteps                  = 250" in pre_relax_mdp
    assert "verlet-buffer-tolerance = -1" in pre_relax_mdp
    assert "rlist                   = 1.25" in pre_relax_mdp
    assert "rcoulomb                = 1.25" in pre_relax_mdp
    assert "rvdw                    = 1.25" in pre_relax_mdp
    assert "define                  = -DFLEXIBLE" in pre_md_mdp
    assert "dt                      = 0.00025" in pre_md_mdp
    assert "nsteps                  = 800" in pre_md_mdp
    assert "continuation            = no" in pre_md_mdp
    assert "verlet-buffer-tolerance = -1" in pre_md_mdp
    assert "rlist                   = 1.25" in pre_md_mdp
    assert "rcoulomb                = 1.25" in pre_md_mdp
    assert "rvdw                    = 1.25" in pre_md_mdp
    assert "init-lambda-state       = 0" in production_mdp
    assert "coul-lambdas            = 0.00000 0.50000 1.00000 1.00000 1.00000 1.00000" in production_mdp
    assert "vdw-lambdas             = 0.00000 0.00000 0.25000 0.50000 0.75000 1.00000" in production_mdp
    assert "bonded-lambdas" in production_mdp
    assert "mass-lambdas" in production_mdp
    assert "continuation            = yes" in production_mdp
    assert "gen-vel                 = no" in production_mdp
    assert "pcoupl                  = Parrinello-Rahman" in production_mdp
    assert "tau-p                   = 2.0" in production_mdp
    assert "refcoord-scaling        = com" in production_mdp
    assert "verlet-buffer-tolerance = -1" in production_mdp
    assert "rlist                   = 1.25" in production_mdp
    assert "rcoulomb                = 1.25" in production_mdp
    assert "rvdw                    = 1.25" in production_mdp
    assert "define                  = -DFLEXIBLE" in em_mdp
    assert "nsteps                  = 1000" in em_mdp
    assert "constraints             = none" in em_mdp
    assert "verlet-buffer-tolerance = -1" in em_mdp
    assert "rlist                   = 1.25" in em_mdp
    assert "rcoulomb                = 1.25" in em_mdp
    assert "rvdw                    = 1.25" in em_mdp
    assert "define                  = -DPOSRES_STAGE_HEAVY" in nvt_mdp
    assert "define                  = -DPOSRES_STAGE_HEAVY" in npt_mdp
    assert "pcoupl                  = Parrinello-Rahman" in npt_mdp
    assert "tau-p                   = 2.0" in npt_mdp
    assert "refcoord-scaling        = com" in npt_mdp
    assert "define                  = -DPOSRES_STAGE_BACKBONE" in npt_release_mdp
    assert "nsteps                  = 125000" in npt_release_mdp
    assert "pdb2gmx -missing" in mutate_script
    assert "pdb2gmx -ignh" not in mutate_script  # -ignh would drop pmx hybrid hydrogens
    assert "strip_terminal_oxygen_atoms" in mutate_script
    assert "restore_incomplete_standard_residues_from_template" in mutate_script
    assert "repair_sidechain_only_incomplete_residues_with_pdbfixer" in mutate_script
    assert "mutant_geometry_qc.json" in mutate_script
    assert "processed_gro_qc.json" in mutate_script
    assert "inspect_gro_file" in mutate_script
    assert "generate_hybrid_topology" in mutate_script
    assert "write_inter_residue_heavy_atom_clash_report" in mutate_script
    assert mutate_script.index("restore_incomplete_standard_residues_from_template") < mutate_script.index("repair_sidechain_only_incomplete_residues_with_pdbfixer")
    assert mutate_script.index("repair_sidechain_only_incomplete_residues_with_pdbfixer") < mutate_script.index("mutant_geometry_qc.json")
    assert mutate_script.index("restore_incomplete_standard_residues_from_template") < mutate_script.index("mutant_geometry_qc.json")
    assert mutate_script.index("mutant_geometry_qc.json") < mutate_script.index("pdb2gmx -missing")
    assert mutate_script.index("pdb2gmx -missing") < mutate_script.index("processed_gro_qc.json")
    assert "solvate" in equilibrate_script
    assert "-maxwarn 3" in equilibrate_script
    assert "-maxwarn 2" in equilibrate_script
    assert "mdrun" in equilibrate_script
    assert "inconsistent shifts over periodic boundaries" in equilibrate_script
    assert "force on at least one atom is not finite" in equilibrate_script
    assert "Maximum force     =            inf" in equilibrate_script
    assert "largest distance between excluded atoms" in equilibrate_script
    assert "'cubic|2.00|fallback'" in equilibrate_script
    assert "'cubic|5.00|expanded fallback'" in equilibrate_script
    assert "retrying equilibration with ${abag_next_label} box ${abag_next_box_type} and padding ${abag_next_box_padding} nm" in equilibrate_script
    assert "em.runtime.history.log" in equilibrate_script
    assert "equilibration_restraints.json" in equilibrate_script
    assert "npt_release.mdp" in equilibrate_script
    assert "skipping completed equilibrate repeat" in equilibrate_script
    assert "-maxwarn 4" in sample_script
    assert "pre_relax.mdp" in sample_script
    assert "pre_md.mdp" in sample_script
    assert "pre_md.cpt" in sample_script
    assert "skipping completed sample window" in sample_script
    assert "starting sample window" in sample_script
    assert "completed sample window" in sample_script
    assert "if [ -s" in sample_script
    assert "rm -f" in sample_script
    assert "md.log" in sample_script
    assert "-dhdl" in sample_script


def test_build_legs_uses_batch_scoped_deterministic_sampling_seeds(tmp_path: Path) -> None:
    first_job_dir = _build_seed_scope_demo_job(tmp_path, batch_id="seed_scope_batch_a")
    second_job_dir = _build_seed_scope_demo_job(tmp_path, batch_id="seed_scope_batch_b")

    run_job(first_job_dir, execute=False, to_stage="build_legs")
    run_job(second_job_dir, execute=False, to_stage="build_legs")

    assert json.loads((first_job_dir / "job_spec.json").read_text(encoding="utf-8"))["job_id"] == json.loads(
        (second_job_dir / "job_spec.json").read_text(encoding="utf-8")
    )["job_id"]

    first_nvt = first_job_dir / "legs" / "complex" / "rep01" / "mdp" / "nvt.mdp"
    first_pre_md = first_job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "pre_md.mdp"
    first_prod = first_job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "production.mdp"
    second_nvt = second_job_dir / "legs" / "complex" / "rep01" / "mdp" / "nvt.mdp"
    second_pre_md = second_job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "pre_md.mdp"
    second_prod = second_job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "production.mdp"

    first_nvt_seed = _read_mdp_parameter(first_nvt, "gen-seed")
    first_pre_md_seed = _read_mdp_parameter(first_pre_md, "ld-seed")
    first_prod_seed = _read_mdp_parameter(first_prod, "ld-seed")
    second_nvt_seed = _read_mdp_parameter(second_nvt, "gen-seed")
    second_pre_md_seed = _read_mdp_parameter(second_pre_md, "ld-seed")
    second_prod_seed = _read_mdp_parameter(second_prod, "ld-seed")

    assert first_nvt_seed != "-1"
    assert first_pre_md_seed != "-1"
    assert first_prod_seed != "-1"
    assert second_nvt_seed != "-1"
    assert second_pre_md_seed != "-1"
    assert second_prod_seed != "-1"
    assert first_nvt_seed != second_nvt_seed
    assert first_pre_md_seed != second_pre_md_seed
    assert first_prod_seed != second_prod_seed

    first_snapshot = {
        "nvt": first_nvt.read_text(encoding="utf-8"),
        "pre_md": first_pre_md.read_text(encoding="utf-8"),
        "production": first_prod.read_text(encoding="utf-8"),
    }
    run_job(first_job_dir, execute=False, from_stage="build_legs", to_stage="build_legs")

    assert first_snapshot["nvt"] == first_nvt.read_text(encoding="utf-8")
    assert first_snapshot["pre_md"] == first_pre_md.read_text(encoding="utf-8")
    assert first_snapshot["production"] == first_prod.read_text(encoding="utf-8")


def test_equilibrate_stage_seeds_completed_repeats_from_prior_plan_root(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: seed_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     14  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     15  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     16  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="seed_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text(f"; {leg} pmxtop\n", encoding="utf-8")
        (pmx_dir / "seed_support.itp").write_text(f"; {leg} support\n", encoding="utf-8")

    seed_root = tmp_path / "seed_runs"
    seed_job_dir = seed_root / "abbind_seed_demo" / "jobs" / job_dir.name
    seed_job_dir.mkdir(parents=True, exist_ok=True)
    (seed_job_dir / "job_spec.json").write_text((job_dir / "job_spec.json").read_text(encoding="utf-8"), encoding="utf-8")
    for leg in ("complex", "apo"):
        seed_repeat_dir = seed_job_dir / "legs" / leg / "rep01"
        (seed_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
        (seed_repeat_dir / "system.top").write_text(f"; seeded {leg} top\n", encoding="utf-8")
        (seed_repeat_dir / "seed_support.itp").write_text(f"; seeded {leg} itp\n", encoding="utf-8")
        (seed_repeat_dir / "equilibration" / "npt.gro").write_text(f"{leg} seeded npt\n", encoding="utf-8")

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "equilibrate commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(
        job_dir,
        execute=True,
        from_stage="equilibrate",
        to_stage="equilibrate",
        environment={"ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS": str(seed_root)},
    )

    assert [status.stage for status in statuses] == ["equilibrate"]
    assert statuses[0].state == "completed"

    equilibrate_script = (job_dir / "artifacts" / "commands" / "equilibrate.sh").read_text(encoding="utf-8")
    assert "skipping completed equilibrate repeat complex/rep01" in equilibrate_script
    assert "skipping completed equilibrate repeat apo/rep01" in equilibrate_script

    for leg in ("complex", "apo"):
        repeat_dir = job_dir / "legs" / leg / "rep01"
        assert (repeat_dir / "system.top").read_text(encoding="utf-8") == f"; seeded {leg} top\n"
        assert (repeat_dir / "seed_support.itp").read_text(encoding="utf-8") == f"; seeded {leg} itp\n"
        assert (repeat_dir / "equilibration" / "npt.gro").read_text(encoding="utf-8") == f"{leg} seeded npt\n"
        seed_source = json.loads((repeat_dir / "equilibration" / "seed_source.json").read_text(encoding="utf-8"))
        assert seed_source["seed_job_dir"] == str(seed_job_dir)
        assert seed_source["seed_repeat_dir"] == str(seed_job_dir / "legs" / leg / "rep01")
        assert seed_source["repeat_id"] == "rep01"


def test_equilibrate_stage_inherits_only_non_target_leg_from_rescue_source(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: targeted_seed_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    TARGETED SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="targeted_seed_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text(f"; {leg} pmxtop\n", encoding="utf-8")
        (pmx_dir / "seed_support.itp").write_text(f"; {leg} support\n", encoding="utf-8")

    source_job_dir = tmp_path / "targeted_source" / "batch" / "jobs" / job_dir.name
    source_job_dir.mkdir(parents=True, exist_ok=True)
    (source_job_dir / "job_spec.json").write_text((job_dir / "job_spec.json").read_text(encoding="utf-8"), encoding="utf-8")
    source_apo_repeat_dir = source_job_dir / "legs" / "apo" / "rep01"
    (source_apo_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
    (source_apo_repeat_dir / "system.top").write_text("; inherited apo top\n", encoding="utf-8")
    (source_apo_repeat_dir / "seed_support.itp").write_text("; inherited apo itp\n", encoding="utf-8")
    (source_apo_repeat_dir / "equilibration" / "npt.gro").write_text("apo inherited npt\n", encoding="utf-8")

    (job_dir / "config" / "rescue.json").write_text(
        json.dumps(
            {
                "mode": "targeted_primary_repeat_spread_leg",
                "source_job_dir": str(source_job_dir),
                "target_legs": ["complex"],
                "inherit_source_legs": ["apo"],
            }
        ),
        encoding="utf-8",
    )

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "equilibrate commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(job_dir, execute=True, from_stage="equilibrate", to_stage="equilibrate")

    assert [status.stage for status in statuses] == ["equilibrate"]
    assert statuses[0].state == "completed"

    equilibrate_script = (job_dir / "artifacts" / "commands" / "equilibrate.sh").read_text(encoding="utf-8")
    assert "skipping completed equilibrate repeat apo/rep01" in equilibrate_script

    apo_repeat_dir = job_dir / "legs" / "apo" / "rep01"
    complex_repeat_dir = job_dir / "legs" / "complex" / "rep01"
    assert (apo_repeat_dir / "system.top").read_text(encoding="utf-8") == "; inherited apo top\n"
    assert (apo_repeat_dir / "seed_support.itp").read_text(encoding="utf-8") == "; inherited apo itp\n"
    assert (apo_repeat_dir / "equilibration" / "npt.gro").read_text(encoding="utf-8") == "apo inherited npt\n"
    seed_source = json.loads((apo_repeat_dir / "equilibration" / "seed_source.json").read_text(encoding="utf-8"))
    assert seed_source["seed_job_dir"] == str(source_job_dir)
    assert not (complex_repeat_dir / "equilibration" / "seed_source.json").exists()
    assert not (complex_repeat_dir / "system.top").exists()
    assert not (complex_repeat_dir / "seed_support.itp").exists()
    assert not (complex_repeat_dir / "equilibration" / "npt.gro").exists()


def test_sample_stage_inherits_only_non_target_leg_from_rescue_source(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: targeted_sample_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    TARGETED SAMPLE DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="targeted_sample_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    for leg in ("complex", "apo"):
        repeat_dir = job_dir / "legs" / leg / "rep01"
        (repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
        (repeat_dir / "system.top").write_text(f"; current {leg} top\n", encoding="utf-8")
        (repeat_dir / "equilibration" / "npt.gro").write_text(f"{leg} current npt\n", encoding="utf-8")

    source_job_dir = tmp_path / "targeted_sample_source" / "batch" / "jobs" / job_dir.name
    source_job_dir.mkdir(parents=True, exist_ok=True)
    (source_job_dir / "job_spec.json").write_text((job_dir / "job_spec.json").read_text(encoding="utf-8"), encoding="utf-8")
    for window_name in ("lambda_000", "lambda_001"):
        source_window_dir = source_job_dir / "legs" / "apo" / "rep01" / window_name
        source_window_dir.mkdir(parents=True, exist_ok=True)
        (source_window_dir / "topol.tpr").write_text(f"{window_name} topol\n", encoding="utf-8")
        (source_window_dir / "dhdl.xvg").write_text(f"{window_name} dhdl\n", encoding="utf-8")
        (source_window_dir / "md.gro").write_text(f"{window_name} gro\n", encoding="utf-8")
        (source_window_dir / "md.log").write_text(f"{window_name} log\n", encoding="utf-8")

    (job_dir / "config" / "rescue.json").write_text(
        json.dumps(
            {
                "mode": "targeted_primary_repeat_spread_leg",
                "source_job_dir": str(source_job_dir),
                "target_legs": ["complex"],
                "inherit_source_legs": ["apo"],
            }
        ),
        encoding="utf-8",
    )

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "sample commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(job_dir, execute=True, from_stage="sample", to_stage="sample")

    assert [status.stage for status in statuses] == ["sample"]
    assert statuses[0].state == "completed"

    sample_script = (job_dir / "artifacts" / "commands" / "sample.sh").read_text(encoding="utf-8")
    assert "skipping completed sample window apo/rep01/lambda_000" in sample_script
    assert "skipping completed sample window apo/rep01/lambda_001" in sample_script

    for window_name in ("lambda_000", "lambda_001"):
        apo_window_dir = job_dir / "legs" / "apo" / "rep01" / window_name
        assert (apo_window_dir / "topol.tpr").read_text(encoding="utf-8") == f"{window_name} topol\n"
        assert (apo_window_dir / "dhdl.xvg").read_text(encoding="utf-8") == f"{window_name} dhdl\n"
        assert (apo_window_dir / "md.gro").read_text(encoding="utf-8") == f"{window_name} gro\n"
        assert (apo_window_dir / "md.log").read_text(encoding="utf-8") == f"{window_name} log\n"
        source_marker = json.loads((apo_window_dir / "sample_source.json").read_text(encoding="utf-8"))
        assert source_marker["seed_job_dir"] == str(source_job_dir)
        assert source_marker["window_id"] == window_name
        complex_window_dir = job_dir / "legs" / "complex" / "rep01" / window_name
        assert not (complex_window_dir / "sample_source.json").exists()
        assert not (complex_window_dir / "topol.tpr").exists()
        assert not (complex_window_dir / "dhdl.xvg").exists()
        assert not (complex_window_dir / "md.gro").exists()
        assert not (complex_window_dir / "md.log").exists()


def test_equilibrate_stage_seeding_overwrites_existing_topology_bundle_when_target_npt_is_missing(
    tmp_path: Path,
    monkeypatch,
) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'input.pdb'}",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "input.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "input.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "input.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "grp_seed,A,10,,A,V,antigen",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="seed_overwrite_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text(f"; {leg} pmxtop\n", encoding="utf-8")
        (pmx_dir / "seed_support.itp").write_text(f"; {leg} support\n", encoding="utf-8")
        repeat_dir = job_dir / "legs" / leg / "rep01"
        (repeat_dir / "system.top").write_text(f"; stale {leg} top\n", encoding="utf-8")
        (repeat_dir / "seed_support.itp").write_text(f"; stale {leg} itp\n", encoding="utf-8")

    seed_root = tmp_path / "seed_runs"
    seed_job_dir = seed_root / "abbind_seed_demo" / "jobs" / job_dir.name
    seed_job_dir.mkdir(parents=True, exist_ok=True)
    (seed_job_dir / "job_spec.json").write_text((job_dir / "job_spec.json").read_text(encoding="utf-8"), encoding="utf-8")
    for leg in ("complex", "apo"):
        seed_repeat_dir = seed_job_dir / "legs" / leg / "rep01"
        (seed_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
        (seed_repeat_dir / "system.top").write_text(f"; seeded {leg} top\n", encoding="utf-8")
        (seed_repeat_dir / "seed_support.itp").write_text(f"; seeded {leg} itp\n", encoding="utf-8")
        (seed_repeat_dir / "equilibration" / "npt.gro").write_text(f"{leg} seeded npt\n", encoding="utf-8")

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "equilibrate commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(
        job_dir,
        execute=True,
        from_stage="equilibrate",
        to_stage="equilibrate",
        environment={"ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS": str(seed_root)},
    )

    assert [status.stage for status in statuses] == ["equilibrate"]
    assert statuses[0].state == "completed"

    for leg in ("complex", "apo"):
        repeat_dir = job_dir / "legs" / leg / "rep01"
        assert (repeat_dir / "system.top").read_text(encoding="utf-8") == f"; seeded {leg} top\n"
        assert (repeat_dir / "seed_support.itp").read_text(encoding="utf-8") == f"; seeded {leg} itp\n"
        assert (repeat_dir / "equilibration" / "npt.gro").read_text(encoding="utf-8") == f"{leg} seeded npt\n"


def test_equilibrate_stage_reseeds_existing_bundle_when_seed_source_is_present(
    tmp_path: Path,
    monkeypatch,
) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'input.pdb'}",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "input.pdb").write_text(
        "\n".join(
            [
                "HEADER    SEED DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "grp_seed,A,10,,A,V,antigen",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="reseed_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    seed_root = tmp_path / "seed_runs"
    seed_job_dir = seed_root / "abbind_seed_demo" / "jobs" / job_dir.name
    seed_job_dir.mkdir(parents=True, exist_ok=True)
    (seed_job_dir / "job_spec.json").write_text((job_dir / "job_spec.json").read_text(encoding="utf-8"), encoding="utf-8")
    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text(f"; {leg} pmxtop\n", encoding="utf-8")
        (pmx_dir / "seed_support.itp").write_text(f"; {leg} support\n", encoding="utf-8")

        repeat_dir = job_dir / "legs" / leg / "rep01"
        (repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
        (repeat_dir / "system.top").write_text(f"; stale {leg} top\n", encoding="utf-8")
        (repeat_dir / "seed_support.itp").write_text(f"; stale {leg} itp\n", encoding="utf-8")
        (repeat_dir / "equilibration" / "npt.gro").write_text(f"{leg} stale npt\n", encoding="utf-8")

        seed_repeat_dir = seed_job_dir / "legs" / leg / "rep01"
        (seed_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
        (seed_repeat_dir / "system.top").write_text(f"; seeded {leg} top\n", encoding="utf-8")
        (seed_repeat_dir / "seed_support.itp").write_text(f"; seeded {leg} itp\n", encoding="utf-8")
        (seed_repeat_dir / "equilibration" / "npt.gro").write_text(f"{leg} seeded npt\n", encoding="utf-8")

        (repeat_dir / "equilibration" / "seed_source.json").write_text(
            json.dumps(
                {
                    "job_id": job_dir.name,
                    "leg": leg,
                    "repeat_id": "rep01",
                    "seed_job_dir": str(seed_job_dir),
                    "seed_repeat_dir": str(seed_repeat_dir),
                    "seeded_files": [],
                    "seeded_at": "2026-06-09T00:00:00Z",
                }
            ),
            encoding="utf-8",
        )

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "equilibrate commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(
        job_dir,
        execute=True,
        from_stage="equilibrate",
        to_stage="equilibrate",
        environment={"ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS": str(seed_root)},
    )

    assert [status.stage for status in statuses] == ["equilibrate"]
    assert statuses[0].state == "completed"

    for leg in ("complex", "apo"):
        repeat_dir = job_dir / "legs" / leg / "rep01"
        assert (repeat_dir / "system.top").read_text(encoding="utf-8") == f"; seeded {leg} top\n"
        assert (repeat_dir / "seed_support.itp").read_text(encoding="utf-8") == f"; seeded {leg} itp\n"
        assert (repeat_dir / "equilibration" / "npt.gro").read_text(encoding="utf-8") == f"{leg} seeded npt\n"
        seed_source = json.loads((repeat_dir / "equilibration" / "seed_source.json").read_text(encoding="utf-8"))
        assert seed_source["seed_repeat_dir"] == str(seed_job_dir / "legs" / leg / "rep01")


def test_equilibrate_stage_seeding_copies_complete_source_topology_bundle(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: source_bundle_demo",
                f"input_structure: {tmp_path / 'input.pdb'}",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "input.pdb").write_text(
        "\n".join(
            [
                "HEADER    SOURCE BUNDLE DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   ILE A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  ILE A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   ILE A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   ILE A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_a_i10v,A,10,,I,V,antigen",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="source_bundle_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text(f"; {leg} pmxtop\n", encoding="utf-8")
        (pmx_dir / "pmx_topol_Protein_chain_A.itp").write_text(f"; current {leg} chain A\n", encoding="utf-8")

    seed_root = tmp_path / "seed_runs"
    seed_job_dir = seed_root / "abbind_seed_demo" / "jobs" / job_dir.name
    seed_job_dir.mkdir(parents=True, exist_ok=True)
    (seed_job_dir / "job_spec.json").write_text((job_dir / "job_spec.json").read_text(encoding="utf-8"), encoding="utf-8")
    for leg in ("complex", "apo"):
        seed_repeat_dir = seed_job_dir / "legs" / leg / "rep01"
        (seed_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
        if leg == "complex":
            (seed_repeat_dir / "system.top").write_text(
                '\n'.join(
                    [
                        '#include "pmx_topol_Protein_chain_A.itp"',
                        '#include "pmx_topol_Protein_chain_B.itp"',
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (seed_repeat_dir / "pmx_topol_Protein_chain_A.itp").write_text("; seeded complex chain A\n", encoding="utf-8")
            (seed_repeat_dir / "pmx_topol_Protein_chain_B.itp").write_text("; seeded complex chain B\n", encoding="utf-8")
        else:
            (seed_repeat_dir / "system.top").write_text("; seeded apo top\n", encoding="utf-8")
            (seed_repeat_dir / "pmx_topol_Protein_chain_A.itp").write_text("; seeded apo chain A\n", encoding="utf-8")
        (seed_repeat_dir / "equilibration" / "npt.gro").write_text(f"{leg} seeded npt\n", encoding="utf-8")

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "equilibrate commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(
        job_dir,
        execute=True,
        from_stage="equilibrate",
        to_stage="equilibrate",
        environment={"ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS": str(seed_root)},
    )

    assert [status.stage for status in statuses] == ["equilibrate"]
    assert statuses[0].state == "completed"

    complex_repeat_dir = job_dir / "legs" / "complex" / "rep01"
    assert (complex_repeat_dir / "pmx_topol_Protein_chain_A.itp").read_text(encoding="utf-8") == "; seeded complex chain A\n"
    assert (complex_repeat_dir / "pmx_topol_Protein_chain_B.itp").read_text(encoding="utf-8") == "; seeded complex chain B\n"
    seed_source = json.loads((complex_repeat_dir / "equilibration" / "seed_source.json").read_text(encoding="utf-8"))
    assert str(complex_repeat_dir / "pmx_topol_Protein_chain_B.itp") in seed_source["seeded_files"]


def test_sample_stage_backfills_missing_seeded_repeat_support_itps(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: sample_backfill_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    SAMPLE BACKFILL DEMO",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 1",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="sample_backfill_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="build_legs")

    complex_repeat_dir = job_dir / "legs" / "complex" / "rep01"
    apo_repeat_dir = job_dir / "legs" / "apo" / "rep01"
    (complex_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
    (apo_repeat_dir / "equilibration").mkdir(parents=True, exist_ok=True)
    (complex_repeat_dir / "system.top").write_text(
        '\n'.join(
            [
                '#include "pmx_topol_Protein_chain_A.itp"',
                '#include "pmx_topol_Protein_chain_B.itp"',
                "",
            ]
        ),
        encoding="utf-8",
    )
    (complex_repeat_dir / "pmx_topol_Protein_chain_A.itp").write_text("; current chain A\n", encoding="utf-8")
    (complex_repeat_dir / "equilibration" / "npt.gro").write_text("complex current npt\n", encoding="utf-8")
    (apo_repeat_dir / "system.top").write_text("; apo top\n", encoding="utf-8")
    (apo_repeat_dir / "equilibration" / "npt.gro").write_text("apo current npt\n", encoding="utf-8")

    seed_repeat_dir = tmp_path / "seed_repeat" / "rep01"
    seed_repeat_dir.mkdir(parents=True, exist_ok=True)
    (seed_repeat_dir / "pmx_topol_Protein_chain_B.itp").write_text("; seeded chain B\n", encoding="utf-8")
    (complex_repeat_dir / "equilibration" / "seed_source.json").write_text(
        json.dumps(
            {
                "job_id": job_dir.name,
                "leg": "complex",
                "repeat_id": "rep01",
                "seed_repeat_dir": str(seed_repeat_dir),
                "seeded_files": [],
                "seeded_at": "2026-06-10T00:00:00Z",
            }
        ),
        encoding="utf-8",
    )

    from abag_rbfe import stages as stages_module

    gmxlib_dir = tmp_path / "gmxlib"
    gmxlib_dir.mkdir(parents=True, exist_ok=True)
    (gmxlib_dir / "spc216.gro").write_text("seed solvent\n", encoding="utf-8")
    monkeypatch.setattr(stages_module, "_stage_env", lambda _ctx: {"GMXLIB": str(gmxlib_dir)})

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("completed", "sample commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)

    statuses = run_job(job_dir, execute=True, from_stage="sample", to_stage="sample")

    assert [status.stage for status in statuses] == ["sample"]
    assert statuses[0].state == "completed"
    assert (complex_repeat_dir / "pmx_topol_Protein_chain_B.itp").read_text(encoding="utf-8") == "; seeded chain B\n"


def test_run_job_from_middle_stage_clears_stale_downstream_results(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 10",
                "npt_ps: 10",
                "production_ps: 10",
                "temperature_k: 310.0",
                "pressure_bar: 1.0",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(system_path, mutations_path, protocol_path, batch_id="unit_batch_rerun", runs_root=tmp_path / "runs")
    job_dir = Path(batch_plan.jobs[0].workdir)

    run_job(job_dir, execute=False)
    assert (job_dir / "results" / "ddg_summary.json").exists()
    assert (job_dir / "report" / "summary.json").exists()
    assert (job_dir / "stages" / "report.json").exists()

    rerun_statuses = run_job(job_dir, execute=False, from_stage="equilibrate", to_stage="sample")

    assert [item.stage for item in rerun_statuses] == ["equilibrate", "sample"]
    assert (job_dir / "stages" / "ingest.json").exists()
    assert (job_dir / "stages" / "equilibrate.json").exists()
    assert (job_dir / "stages" / "sample.json").exists()
    assert not (job_dir / "stages" / "bar.json").exists()
    assert not (job_dir / "stages" / "qc.json").exists()
    assert not (job_dir / "stages" / "report.json").exists()
    assert not (job_dir / "results" / "ddg_summary.json").exists()
    assert not (job_dir / "results" / "qc_report.json").exists()
    assert not (job_dir / "report" / "summary.json").exists()


def test_bar_stage_script_skips_completed_repeat_outputs(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(system_path, mutations_path, protocol_path, batch_id="unit_batch_bar", runs_root=tmp_path / "runs")
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False)

    for leg in ("complex", "apo"):
        repeat_dir = job_dir / "legs" / leg / "rep01"
        for window_index in (0, 1):
            window_dir = repeat_dir / f"lambda_{window_index:03d}"
            window_dir.mkdir(parents=True, exist_ok=True)
            (window_dir / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")

    run_job(job_dir, execute=False, from_stage="bar", to_stage="bar")
    bar_script = (job_dir / "artifacts" / "commands" / "bar.sh").read_text(encoding="utf-8")

    assert "skipping completed BAR repeat" in bar_script
    assert "if [ -s" in bar_script
    assert " bar -f " in bar_script


def test_run_job_marks_stage_running_before_dispatch(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(system_path, mutations_path, protocol_path, batch_id="unit_batch_running", runs_root=tmp_path / "runs")
    job_dir = Path(batch_plan.jobs[0].workdir)

    from abag_rbfe import stages as stages_module
    from abag_rbfe.io_utils import read_json

    original_ingest = stages_module.STAGE_DISPATCH["ingest"]

    def wrapped_ingest(ctx):
        stage_payload = read_json(ctx.job_dir / "stages" / "ingest.json")
        assert stage_payload["state"] == "running"
        assert stage_payload["message"] == "Stage execution started."
        return original_ingest(ctx)

    monkeypatch.setitem(stages_module.STAGE_DISPATCH, "ingest", wrapped_ingest)
    statuses = run_job(job_dir, execute=False, to_stage="ingest")

    assert [item.stage for item in statuses] == ["ingest"]
    assert statuses[0].state == "completed"


def test_run_job_exposes_external_stage_context_while_running(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="unit_batch_external_running",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    inspected: dict[str, dict] = {}

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        stage = script_path.stem
        payload = json.loads((job_dir / "stages" / f"{stage}.json").read_text(encoding="utf-8"))
        assert payload["stage"] == stage
        assert payload["state"] == "running"
        assert payload["completed_at"] is None
        assert payload["commands"] == commands
        assert str(script_path) in payload["artifacts"]
        assert str(script_path.with_suffix(".log")) in payload["artifacts"]
        inspected[stage] = payload
        self.write_script(script_path, commands, workdir, env=env)
        return CommandOutcome("planned", f"{stage} commands written")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)
    statuses = run_job(job_dir, execute=False, to_stage="sample")

    assert [item.stage for item in statuses] == [
        "ingest",
        "prepare",
        "mutate",
        "build_legs",
        "equilibrate",
        "sample",
    ]
    assert set(inspected) == {"mutate", "equilibrate", "sample"}
    for stage in inspected:
        payload = json.loads((job_dir / "stages" / f"{stage}.json").read_text(encoding="utf-8"))
        assert payload["state"] == "planned"
        assert str(job_dir / "artifacts" / "commands" / f"{stage}.sh") in payload["artifacts"]
        assert str(job_dir / "artifacts" / "commands" / f"{stage}.log") in payload["artifacts"]


def test_run_job_appends_env_requested_mdrun_pin_flags(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "mdrun_args: -ntmpi 1 -ntomp 4",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="unit_batch_pinning",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)

    run_job(
        job_dir,
        execute=False,
        to_stage="sample",
        environment={
            "ABAG_RBFE_MDRUN_PINOFFSET": "12",
            "ABAG_RBFE_MDRUN_PINSTRIDE": "1",
        },
    )

    equilibrate_script = (job_dir / "artifacts" / "commands" / "equilibrate.sh").read_text(encoding="utf-8")
    sample_script = (job_dir / "artifacts" / "commands" / "sample.sh").read_text(encoding="utf-8")

    assert "-ntmpi 1 -ntomp 4 -pin on -pinoffset 12 -pinstride 1" in equilibrate_script
    assert "-ntmpi 1 -ntomp 4 -pin on -pinoffset 12 -pinstride 1" in sample_script


def test_run_job_uses_mdrun_args_override_from_environment(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: unit_system_override",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "mdrun_args: -ntmpi 1 -ntomp 4",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="unit_batch_override",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)

    run_job(
        job_dir,
        execute=False,
        to_stage="sample",
        environment={
            "ABAG_RBFE_MDRUN_ARGS": "-ntmpi 1 -ntomp 2",
            "ABAG_RBFE_MDRUN_PINOFFSET": "6",
            "ABAG_RBFE_MDRUN_PINSTRIDE": "1",
        },
    )

    equilibrate_script = (job_dir / "artifacts" / "commands" / "equilibrate.sh").read_text(encoding="utf-8")
    sample_script = (job_dir / "artifacts" / "commands" / "sample.sh").read_text(encoding="utf-8")

    assert "-ntmpi 1 -ntomp 4 -pin on -pinoffset 6 -pinstride 1" not in equilibrate_script
    assert "-ntmpi 1 -ntomp 4 -pin on -pinoffset 6 -pinstride 1" not in sample_script
    assert "-ntmpi 1 -ntomp 2 -pin on -pinoffset 6 -pinstride 1" in equilibrate_script
    assert "-ntmpi 1 -ntomp 2 -pin on -pinoffset 6 -pinstride 1" in sample_script


def test_resume_job_prefers_sample_when_sample_artifacts_exist(tmp_path: Path, monkeypatch) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts" / "commands").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "apo" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "rep01" / "equilibration").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000").mkdir(parents=True, exist_ok=True)

    (job_dir / "stages" / "ingest.json").write_text(json.dumps({"stage": "ingest", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "prepare.json").write_text(json.dumps({"stage": "prepare", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "mutate.json").write_text(json.dumps({"stage": "mutate", "state": "planned"}), encoding="utf-8")
    (job_dir / "stages" / "build_legs.json").write_text(json.dumps({"stage": "build_legs", "state": "completed"}), encoding="utf-8")
    _write_valid_mock_gro(job_dir / "legs" / "complex" / "pmx" / "processed.gro")
    _write_valid_mock_gro(job_dir / "legs" / "apo" / "pmx" / "processed.gro")
    (job_dir / "legs" / "complex" / "pmx" / "pmxtop.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "apo" / "pmx" / "pmxtop.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "equilibration" / "npt.gro").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "system.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "topol.tpr").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "0"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["execute"] is True
    assert captured["from_stage"] == "sample"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "0"}


def test_resume_job_falls_back_to_equilibrate_when_sample_inputs_are_incomplete(tmp_path: Path, monkeypatch) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts" / "commands").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "apo" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "rep01" / "equilibration").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000").mkdir(parents=True, exist_ok=True)
    (job_dir / "legs" / "apo" / "rep01" / "equilibration").mkdir(parents=True)

    (job_dir / "stages" / "ingest.json").write_text(json.dumps({"stage": "ingest", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "prepare.json").write_text(json.dumps({"stage": "prepare", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "mutate.json").write_text(json.dumps({"stage": "mutate", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "build_legs.json").write_text(json.dumps({"stage": "build_legs", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "equilibrate.json").write_text(json.dumps({"stage": "equilibrate", "state": "running"}), encoding="utf-8")
    (job_dir / "stages" / "sample.json").write_text(json.dumps({"stage": "sample", "state": "running"}), encoding="utf-8")
    (job_dir / "artifacts" / "commands" / "sample.sh").write_text("#!/usr/bin/env bash\n", encoding="utf-8")
    (job_dir / "artifacts" / "commands" / "equilibrate.sh").write_text("#!/usr/bin/env bash\n", encoding="utf-8")

    for pmx_dir in (job_dir / "legs" / "complex" / "pmx", job_dir / "legs" / "apo" / "pmx"):
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "equilibration" / "npt.gro").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "system.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "topol.tpr").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")
    (job_dir / "legs" / "apo" / "rep01" / "system.top").write_text("mock\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "0"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["from_stage"] == "equilibrate"


def test_resume_job_prefers_build_legs_when_stage_files_are_missing_but_lambda_plan_exists(tmp_path: Path, monkeypatch) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "apo" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "rep01" / "mdp").mkdir(parents=True, exist_ok=True)

    (job_dir / "stages" / "ingest.json").write_text(json.dumps({"stage": "ingest", "state": "completed"}), encoding="utf-8")
    (job_dir / "stages" / "prepare.json").write_text(json.dumps({"stage": "prepare", "state": "completed"}), encoding="utf-8")
    _write_valid_mock_gro(job_dir / "legs" / "complex" / "pmx" / "processed.gro")
    _write_valid_mock_gro(job_dir / "legs" / "apo" / "pmx" / "processed.gro")
    (job_dir / "legs" / "complex" / "pmx" / "pmxtop.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "apo" / "pmx" / "pmxtop.top").write_text("mock\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "lambda_plan.json").write_text('{"leg":"complex"}\n', encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "mdp" / "nvt.mdp").write_text("mock\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "1"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["execute"] is True
    assert captured["from_stage"] == "build_legs"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "1"}


def test_resume_job_recovers_completed_stage_files_from_existing_outputs(tmp_path: Path, monkeypatch) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts" / "commands").mkdir(parents=True)

    job_spec = {
        "job_id": "demo-job",
        "batch_id": "demo-batch",
        "protocol": {
            "preset": "single_point",
            "repeats": 1,
            "lambda_windows": 1,
        },
    }
    (job_dir / "job_spec.json").write_text(json.dumps(job_spec), encoding="utf-8")

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "running"),
        ("build_legs", "running"),
        ("equilibrate", "failed"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps(
                {
                    "stage": stage,
                    "state": state,
                    "message": f"{stage} {state}",
                    "commands": [],
                    "artifacts": [],
                    "started_at": "2026-06-08T00:00:00Z",
                    "completed_at": None,
                }
            ),
            encoding="utf-8",
        )

    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text("mock\n", encoding="utf-8")

        repeat_dir = job_dir / "legs" / leg / "rep01"
        (repeat_dir / "mdp").mkdir(parents=True)
        (repeat_dir / "equilibration").mkdir(parents=True)
        lambda_dir = repeat_dir / "lambda_000"
        lambda_dir.mkdir(parents=True)

        (repeat_dir / "lambda_plan.json").write_text('{"leg": "mock"}\n', encoding="utf-8")
        for mdp_name in ("genion.mdp", "em.mdp", "nvt.mdp", "npt.mdp"):
            (repeat_dir / "mdp" / mdp_name).write_text("mock\n", encoding="utf-8")
        for mdp_name in ("pre_relax.mdp", "pre_md.mdp", "production.mdp"):
            (lambda_dir / mdp_name).write_text("mock\n", encoding="utf-8")
        (repeat_dir / "system.top").write_text("mock\n", encoding="utf-8")
        (repeat_dir / "equilibration" / "npt.gro").write_text("mock\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "2"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["from_stage"] == "sample"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "2"}

    mutate_payload = json.loads((job_dir / "stages" / "mutate.json").read_text(encoding="utf-8"))
    build_payload = json.loads((job_dir / "stages" / "build_legs.json").read_text(encoding="utf-8"))
    equilibrate_payload = json.loads((job_dir / "stages" / "equilibrate.json").read_text(encoding="utf-8"))

    assert mutate_payload["state"] == "completed"
    assert build_payload["state"] == "completed"
    assert equilibrate_payload["state"] == "completed"
    assert mutate_payload["message"] == "Recovered completed stage from existing pmx outputs."
    assert build_payload["message"] == "Recovered completed stage from existing lambda-plan and MDP outputs."
    assert equilibrate_payload["message"] == "Recovered completed stage from existing equilibrated repeat outputs."


def test_resume_job_invalidates_bad_processed_gro_and_restarts_from_mutate(tmp_path: Path, monkeypatch) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts" / "commands").mkdir(parents=True)

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "completed"),
        ("build_legs", "completed"),
        ("equilibrate", "running"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps({"stage": stage, "state": state}),
            encoding="utf-8",
        )

    complex_pmx = job_dir / "legs" / "complex" / "pmx"
    apo_pmx = job_dir / "legs" / "apo" / "pmx"
    complex_pmx.mkdir(parents=True)
    apo_pmx.mkdir(parents=True)
    (complex_pmx / "processed.gro").write_text("corrupt\n", encoding="utf-8")
    _write_valid_mock_gro(apo_pmx / "processed.gro")
    (complex_pmx / "pmxtop.top").write_text("mock\n", encoding="utf-8")
    (apo_pmx / "pmxtop.top").write_text("mock\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "5"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["from_stage"] == "mutate"
    assert not (job_dir / "stages" / "mutate.json").exists()
    assert not (job_dir / "stages" / "build_legs.json").exists()
    assert not (job_dir / "stages" / "equilibrate.json").exists()


def test_resume_job_restores_blocked_mutate_qc_and_discards_later_stages(
    tmp_path: Path, monkeypatch
) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts").mkdir(parents=True)

    job_spec = {
        "job_id": "demo-job",
        "batch_id": "demo-batch",
        "protocol": {
            "preset": "single_point",
            "repeats": 1,
            "lambda_windows": 1,
        },
    }
    (job_dir / "job_spec.json").write_text(json.dumps(job_spec), encoding="utf-8")

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "completed"),
        ("build_legs", "completed"),
        ("equilibrate", "running"),
        ("sample", "running"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps(
                {
                    "stage": stage,
                    "state": state,
                    "message": f"{stage} {state}",
                    "commands": [],
                    "artifacts": [],
                    "started_at": "2026-06-08T00:00:00Z",
                    "completed_at": None,
                }
            ),
            encoding="utf-8",
        )

    for leg in ("complex", "apo"):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True)
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text("mock\n", encoding="utf-8")

        repeat_dir = job_dir / "legs" / leg / "rep01"
        (repeat_dir / "mdp").mkdir(parents=True)
        (repeat_dir / "equilibration").mkdir(parents=True)
        lambda_dir = repeat_dir / "lambda_000"
        lambda_dir.mkdir(parents=True)

        (repeat_dir / "lambda_plan.json").write_text('{"leg": "mock"}\n', encoding="utf-8")
        for mdp_name in ("genion.mdp", "em.mdp", "nvt.mdp", "npt.mdp"):
            (repeat_dir / "mdp" / mdp_name).write_text("mock\n", encoding="utf-8")
        for mdp_name in ("pre_relax.mdp", "pre_md.mdp", "production.mdp"):
            (lambda_dir / mdp_name).write_text("mock\n", encoding="utf-8")
        (repeat_dir / "system.top").write_text("mock\n", encoding="utf-8")
        (repeat_dir / "equilibration" / "npt.gro").write_text("mock\n", encoding="utf-8")
        (lambda_dir / "dhdl.xvg").write_text("mock\n", encoding="utf-8")
        (lambda_dir / "md.gro").write_text("mock\n", encoding="utf-8")

    (job_dir / "results").mkdir(parents=True)
    (job_dir / "results" / "ddg_summary.json").write_text('{"ddg": 1.0}\n', encoding="utf-8")
    (job_dir / "report").mkdir(parents=True)
    (job_dir / "report" / "summary.json").write_text('{"status": "stale"}\n', encoding="utf-8")
    (job_dir / "artifacts" / "mutate_qc.json").write_text(
        json.dumps(
            {
                "job_id": "demo-job",
                "legs": {
                    "complex": {
                        "inter_residue_heavy_atom_clashes": [
                            {
                                "blocking_prepare": True,
                                "chain_id": "A",
                                "resseq": 565,
                                "icode": "",
                                "resname": "ARG",
                                "partner_chain_id": "A",
                                "partner_resseq": 749,
                                "partner_icode": "",
                                "partner_resname": "TYR",
                                "clashes": [
                                    {
                                        "atom_a": "CD",
                                        "atom_b": "CD2",
                                        "distance_angstrom": 1.0831,
                                    }
                                ],
                            }
                        ]
                    },
                    "apo": {"inter_residue_heavy_atom_clashes": []},
                },
            }
        ),
        encoding="utf-8",
    )

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "0"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["from_stage"] == "mutate"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "0"}

    mutate_payload = json.loads((job_dir / "stages" / "mutate.json").read_text(encoding="utf-8"))
    assert mutate_payload["state"] == "blocked_input"
    assert mutate_payload["message"].startswith(
        "Mutated structure contains impossible inter-residue heavy-atom clashes:"
    )
    assert not (job_dir / "stages" / "build_legs.json").exists()
    assert not (job_dir / "stages" / "equilibrate.json").exists()
    assert not (job_dir / "stages" / "sample.json").exists()
    assert not (job_dir / "results" / "ddg_summary.json").exists()
    assert not (job_dir / "report" / "summary.json").exists()


def test_resume_job_refreshes_stale_mutate_qc_against_reference_before_rerun(tmp_path: Path, monkeypatch) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts").mkdir(parents=True)
    (job_dir / "legs" / "complex" / "pmx").mkdir(parents=True)
    (job_dir / "legs" / "apo" / "pmx").mkdir(parents=True)

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "blocked_input"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps({"stage": stage, "state": state}),
            encoding="utf-8",
        )

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
    input_path = job_dir / "legs" / "complex" / "input.pdb"
    mutant_path = job_dir / "legs" / "complex" / "pmx" / "mutant.pdb"
    input_path.write_text(pdb_text, encoding="utf-8")
    mutant_path.write_text(pdb_text, encoding="utf-8")

    stale_clash = {
        "chain_id": "A",
        "resseq": 4,
        "icode": "",
        "resname": "ILE",
        "normalized_resname": "ILE",
        "partner_chain_id": "A",
        "partner_resseq": 7,
        "partner_icode": "",
        "partner_resname": "MET",
        "partner_normalized_resname": "MET",
        "min_distance_angstrom": 1.1261,
        "clashes": [{"atom_a": "C", "atom_b": "CG", "distance_angstrom": 1.1261}],
        "blocking_prepare": True,
    }
    (job_dir / "legs" / "complex" / "pmx" / "mutant_geometry_qc.json").write_text(
        json.dumps(
            {
                "input_structure": "mutant.pdb",
                "inter_residue_heavy_atom_clashes": [stale_clash],
                "blocking_inter_residue_heavy_atom_clashes": True,
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (job_dir / "artifacts" / "mutate_qc.json").write_text(
        json.dumps(
            {
                "job_id": "stale-job",
                "legs": {
                    "complex": {
                        "mutant_pdb": str(mutant_path),
                        "inter_residue_heavy_atom_clashes": [stale_clash],
                    },
                    "apo": {
                        "mutant_pdb": str(job_dir / "legs" / "apo" / "pmx" / "mutant.pdb"),
                        "inter_residue_heavy_atom_clashes": [],
                    },
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["from_stage"] = from_stage
        captured["execute"] = execute
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "7"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["from_stage"] == "mutate"
    refreshed_qc = json.loads((job_dir / "artifacts" / "mutate_qc.json").read_text(encoding="utf-8"))
    assert refreshed_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"] == []
    refreshed_geometry = json.loads(
        (job_dir / "legs" / "complex" / "pmx" / "mutant_geometry_qc.json").read_text(encoding="utf-8")
    )
    assert refreshed_geometry["reference_structure"] == str(input_path)
    assert refreshed_geometry["blocking_inter_residue_heavy_atom_clashes"] is False


def test_resume_job_synthesizes_missing_mutate_qc_and_restores_blocked_input(
    tmp_path: Path, monkeypatch
) -> None:
    job_dir = tmp_path / "job"
    (job_dir / "stages").mkdir(parents=True)
    (job_dir / "artifacts").mkdir(parents=True)

    job_spec = {
        "job_id": "demo-job",
        "batch_id": "demo-batch",
        "protocol": {
            "preset": "single_point",
            "repeats": 1,
            "lambda_windows": 1,
        },
    }
    (job_dir / "job_spec.json").write_text(json.dumps(job_spec), encoding="utf-8")

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "completed"),
        ("build_legs", "completed"),
        ("equilibrate", "completed"),
        ("sample", "running"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps(
                {
                    "stage": stage,
                    "state": state,
                    "message": f"{stage} {state}",
                    "commands": [],
                    "artifacts": [],
                    "started_at": "2026-06-08T00:00:00Z",
                    "completed_at": None,
                }
            ),
            encoding="utf-8",
        )

    complex_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
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
        + "\n"
    )
    apo_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    for leg, mutant_text in (("complex", complex_mutant), ("apo", apo_mutant)):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True)
        (pmx_dir / "mutant.pdb").write_text(mutant_text, encoding="utf-8")
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr("abag_rbfe.stages.run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "1"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["from_stage"] == "mutate"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "1"}

    mutate_payload = json.loads((job_dir / "stages" / "mutate.json").read_text(encoding="utf-8"))
    assert mutate_payload["state"] == "blocked_input"
    assert (job_dir / "artifacts" / "mutate_qc.json").exists()
    mutate_qc = json.loads((job_dir / "artifacts" / "mutate_qc.json").read_text(encoding="utf-8"))
    assert len(mutate_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"]) == 1
    assert mutate_qc["legs"]["apo"]["inter_residue_heavy_atom_clashes"] == []
    assert mutate_qc["legs"]["complex"]["mutant_pdbfixer_repair"] == {}
    assert mutate_qc["legs"]["apo"]["mutant_pdbfixer_repair"] == {}
    assert not (job_dir / "stages" / "build_legs.json").exists()
    assert not (job_dir / "stages" / "equilibrate.json").exists()
    assert not (job_dir / "stages" / "sample.json").exists()


def test_resume_job_repairs_non_mutated_mutate_sidechain_clashes_and_advances(
    tmp_path: Path, monkeypatch
) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: mutate_qc_repair_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      22.500   8.900   8.300  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      22.200  12.500  11.700  1.00 20.00           C",
                "TER",
                "ATOM     20  N   HIS A  50      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     21  CA  HIS A  50      31.100  10.500   8.800  1.00 20.00           C",
                "ATOM     22  C   HIS A  50      32.200   9.500   9.000  1.00 20.00           C",
                "ATOM     23  O   HIS A  50      32.400   8.700   8.200  1.00 20.00           O",
                "ATOM     24  CB  HIS A  50      30.700  11.400  10.000  1.00 20.00           C",
                "ATOM     25  CG  HIS A  50      31.700  12.300  10.600  1.00 20.00           C",
                "ATOM     26  ND1 HIS A  50      32.400  13.200   9.800  1.00 20.00           N",
                "ATOM     27  CD2 HIS A  50      32.100  12.500  11.900  1.00 20.00           C",
                "ATOM     28  CE1 HIS A  50      33.200  13.800  10.700  1.00 20.00           C",
                "ATOM     29  NE2 HIS A  50      33.100  13.300  11.900  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_a_h50a,A,50,,H,A,antigen",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="mutate_qc_repair_resume",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    (job_dir / "stages").mkdir(parents=True, exist_ok=True)
    (job_dir / "artifacts").mkdir(parents=True, exist_ok=True)

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "blocked_input"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps(
                {
                    "stage": stage,
                    "state": state,
                    "message": f"{stage} {state}",
                    "commands": [],
                    "artifacts": [],
                    "started_at": "2026-06-08T00:00:00Z",
                    "completed_at": "2026-06-08T00:05:00Z",
                }
            ),
            encoding="utf-8",
        )

    complex_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
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
                "ATOM     20  N   ALA A  50      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     21  CA  ALA A  50      31.100  10.500   8.800  1.00 20.00           C",
                "ATOM     22  C   ALA A  50      32.200   9.500   9.000  1.00 20.00           C",
                "ATOM     23  O   ALA A  50      32.400   8.700   8.200  1.00 20.00           O",
                "ATOM     24  CB  ALA A  50      30.700  11.400  10.000  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    apo_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM      3  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM      4  O   ILE A   4      22.500   8.900   8.300  1.00 20.00           O",
                "ATOM      5  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM      6  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM      7  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM      8  CD1 ILE A   4      22.200  12.500  11.700  1.00 20.00           C",
                "TER",
                "ATOM      9  N   ALA A  50      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     10  CA  ALA A  50      31.100  10.500   8.800  1.00 20.00           C",
                "ATOM     11  C   ALA A  50      32.200   9.500   9.000  1.00 20.00           C",
                "ATOM     12  O   ALA A  50      32.400   8.700   8.200  1.00 20.00           O",
                "ATOM     13  CB  ALA A  50      30.700  11.400  10.000  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    for leg, mutant_text in (("complex", complex_mutant), ("apo", apo_mutant)):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        (pmx_dir / "mutant.pdb").write_text(mutant_text, encoding="utf-8")
        if leg == "apo":
            _write_valid_mock_gro(pmx_dir / "processed.gro")
            (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")
            (pmx_dir / "topol.top").write_text("; generated\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        if not script_path.stem.startswith("mutate_repair_"):
            return CommandOutcome("planned", "No repair executed.")
        leg = script_path.stem.removeprefix("mutate_repair_")
        pmx_dir = job_dir / "legs" / leg / "pmx"
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "topol.top").write_text("; generated\n", encoding="utf-8")
        (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")
        return CommandOutcome("completed", "Mutant clash repair regenerated pmx outputs.")

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)
    from abag_rbfe import stages as stages_module

    monkeypatch.setattr(stages_module, "run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "2"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["execute"] is True
    assert captured["from_stage"] == "build_legs"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "2"}

    mutate_payload = json.loads((job_dir / "stages" / "mutate.json").read_text(encoding="utf-8"))
    assert mutate_payload["state"] == "completed"
    mutate_qc = json.loads((job_dir / "artifacts" / "mutate_qc.json").read_text(encoding="utf-8"))
    assert mutate_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"] == []
    repair_summary = mutate_qc["legs"]["complex"]["auto_repair_summary"]
    assert repair_summary["attempted"] is True
    assert repair_summary["succeeded"] is True


def test_resume_job_repairs_non_mutated_backbone_sidechain_mutate_clashes_and_advances(
    tmp_path: Path, monkeypatch
) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: mutate_qc_backbone_sidechain_repair_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      22.500   8.900   8.300  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      22.200  12.500  11.700  1.00 20.00           C",
                "TER",
                "ATOM     20  N   HIS A  50      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     21  CA  HIS A  50      31.100  10.500   8.800  1.00 20.00           C",
                "ATOM     22  C   HIS A  50      32.200   9.500   9.000  1.00 20.00           C",
                "ATOM     23  O   HIS A  50      32.400   8.700   8.200  1.00 20.00           O",
                "ATOM     24  CB  HIS A  50      30.700  11.400  10.000  1.00 20.00           C",
                "ATOM     25  CG  HIS A  50      31.700  12.300  10.600  1.00 20.00           C",
                "ATOM     26  ND1 HIS A  50      32.400  13.200   9.800  1.00 20.00           N",
                "ATOM     27  CD2 HIS A  50      32.100  12.500  11.900  1.00 20.00           C",
                "ATOM     28  CE1 HIS A  50      33.200  13.800  10.700  1.00 20.00           C",
                "ATOM     29  NE2 HIS A  50      33.100  13.300  11.900  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_a_h50a,A,50,,H,A,antigen",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="mutate_qc_backbone_sidechain_repair_resume",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    (job_dir / "stages").mkdir(parents=True, exist_ok=True)
    (job_dir / "artifacts").mkdir(parents=True, exist_ok=True)

    for stage, state in (
        ("ingest", "completed"),
        ("prepare", "completed"),
        ("mutate", "blocked_input"),
    ):
        (job_dir / "stages" / f"{stage}.json").write_text(
            json.dumps(
                {
                    "stage": stage,
                    "state": state,
                    "message": f"{stage} {state}",
                    "commands": [],
                    "artifacts": [],
                    "started_at": "2026-06-08T00:00:00Z",
                    "completed_at": "2026-06-08T00:05:00Z",
                }
            ),
            encoding="utf-8",
        )

    complex_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      14.950  12.250  11.150  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      22.200  12.500  11.700  1.00 20.00           C",
                "TER",
                "ATOM     20  N   ALA A  50      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     21  CA  ALA A  50      31.100  10.500   8.800  1.00 20.00           C",
                "ATOM     22  C   ALA A  50      32.200   9.500   9.000  1.00 20.00           C",
                "ATOM     23  O   ALA A  50      32.400   8.700   8.200  1.00 20.00           O",
                "ATOM     24  CB  ALA A  50      30.700  11.400  10.000  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    apo_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      20.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      22.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      14.950  12.250  11.150  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      20.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      21.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      18.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      22.200  12.500  11.700  1.00 20.00           C",
                "TER",
                "ATOM     20  N   ALA A  50      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     21  CA  ALA A  50      31.100  10.500   8.800  1.00 20.00           C",
                "ATOM     22  C   ALA A  50      32.200   9.500   9.000  1.00 20.00           C",
                "ATOM     23  O   ALA A  50      32.400   8.700   8.200  1.00 20.00           O",
                "ATOM     24  CB  ALA A  50      30.700  11.400  10.000  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    for leg, mutant_text in (("complex", complex_mutant), ("apo", apo_mutant)):
        pmx_dir = job_dir / "legs" / leg / "pmx"
        pmx_dir.mkdir(parents=True, exist_ok=True)
        (pmx_dir / "mutant.pdb").write_text(mutant_text, encoding="utf-8")
        if leg == "apo":
            _write_valid_mock_gro(pmx_dir / "processed.gro")
            (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")
            (pmx_dir / "topol.top").write_text("; generated\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        if not script_path.stem.startswith("mutate_repair_"):
            return CommandOutcome("planned", "No repair executed.")
        leg = script_path.stem.removeprefix("mutate_repair_")
        pmx_dir = job_dir / "legs" / leg / "pmx"
        _write_valid_mock_gro(pmx_dir / "processed.gro")
        (pmx_dir / "topol.top").write_text("; generated\n", encoding="utf-8")
        (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")
        return CommandOutcome("completed", "Mutant clash repair regenerated pmx outputs.")

    def fake_run_job(
        target_job_dir: Path,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[object]:
        captured["job_dir"] = target_job_dir
        captured["execute"] = execute
        captured["from_stage"] = from_stage
        captured["to_stage"] = to_stage
        captured["environment"] = environment
        return []

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)
    from abag_rbfe import stages as stages_module

    monkeypatch.setattr(stages_module, "run_job", fake_run_job)

    statuses = resume_job(job_dir, execute=True, environment={"CUDA_VISIBLE_DEVICES": "3"})

    assert statuses == []
    assert captured["job_dir"] == job_dir
    assert captured["execute"] is True
    assert captured["from_stage"] == "build_legs"
    assert captured["environment"] == {"CUDA_VISIBLE_DEVICES": "3"}

    mutate_payload = json.loads((job_dir / "stages" / "mutate.json").read_text(encoding="utf-8"))
    assert mutate_payload["state"] == "completed"
    mutate_qc = json.loads((job_dir / "artifacts" / "mutate_qc.json").read_text(encoding="utf-8"))
    assert mutate_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"] == []
    assert mutate_qc["legs"]["apo"]["inter_residue_heavy_atom_clashes"] == []
    repair_summary = mutate_qc["legs"]["complex"]["auto_repair_summary"]
    assert repair_summary["attempted"] is True
    assert repair_summary["succeeded"] is True


def test_run_job_marks_interrupted_stage_as_failed(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(system_path, mutations_path, protocol_path, batch_id="unit_batch_interrupt", runs_root=tmp_path / "runs")
    job_dir = Path(batch_plan.jobs[0].workdir)

    from abag_rbfe import stages as stages_module
    from abag_rbfe.io_utils import read_json

    def interrupted_ingest(_ctx):
        stage_payload = read_json(job_dir / "stages" / "ingest.json")
        assert stage_payload["state"] == "running"
        raise KeyboardInterrupt

    monkeypatch.setitem(stages_module.STAGE_DISPATCH, "ingest", interrupted_ingest)
    statuses = run_job(job_dir, execute=False, to_stage="ingest")

    assert [item.stage for item in statuses] == ["ingest"]
    assert statuses[0].state == "failed"
    assert statuses[0].message == "Stage execution interrupted."
    stage_payload = read_json(job_dir / "stages" / "ingest.json")
    assert stage_payload["state"] == "failed"
    assert stage_payload["message"] == "Stage execution interrupted."
    assert stage_payload["completed_at"] is not None


def test_double_point_uses_conservative_defaults_when_protocol_has_no_explicit_quick_overrides(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: demo_abag",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   SER H  52      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  SER H  52      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   SER H  52      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   SER H  52      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  SER H  52      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  OG  SER H  52      12.086  15.074  10.191  1.00 20.00           O",
                "TER",
                "ATOM      7  N   ASN H  54      11.104  16.207   8.154  1.00 20.00           N",
                "ATOM      8  CA  ASN H  54      12.560  16.311   8.321  1.00 20.00           C",
                "ATOM      9  C   ASN H  54      13.150  15.090   8.991  1.00 20.00           C",
                "ATOM     10  O   ASN H  54      12.811  13.945   8.710  1.00 20.00           O",
                "ATOM     11  CB  ASN H  54      12.936  17.622   9.024  1.00 20.00           C",
                "ATOM     12  CG  ASN H  54      12.086  18.074  10.191  1.00 20.00           C",
                "ATOM     13  OD1 ASN H  54      10.900  17.800  10.250  1.00 20.00           O",
                "ATOM     14  ND2 ASN H  54      12.700  18.800  11.130  1.00 20.00           N",
                "TER",
                "ATOM     15  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     16  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     17  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     18  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     19  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     20  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     21  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     22  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     23  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     24  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     25  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     26  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     27  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     28  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     29  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "double_h_pair,H,52,,S,T,antibody",
                "double_h_pair,H,54,,N,Q,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(system_path, mutations_path, protocol_path, batch_id="unit_batch_default", runs_root=tmp_path / "runs")
    protocol = batch_plan.jobs[0].protocol

    assert protocol.preset == "double_point"
    assert protocol.lambda_windows == 32
    assert protocol.repeats == 3
    assert protocol.nvt_ps == 200
    assert protocol.npt_ps == 1500
    assert protocol.production_ps == 3000
    assert protocol.production_dt_ps == 0.002


def test_hydrate_protocol_config_backfills_new_fields_for_legacy_job_configs() -> None:
    protocol = hydrate_protocol_config(
        {
            "preset": "single_point",
            "force_field": "amber99sb-star-ildn-mut",
            "water_model": "tip3p",
            "lambda_windows": 2,
            "repeats": 1,
            "nvt_ps": 1,
            "npt_ps": 1,
            "production_ps": 1,
            "temperature_k": 310.0,
            "pressure_bar": 1.0,
            "overlap_threshold": 0.2,
            "max_repeat_delta_kcal_mol": 1.0,
            "same_side_double_point": False,
            "allow_charge_changing": False,
            "allow_cross_side_double_point": False,
            "gmx_bin": "gmx",
            "pmx_bin": "pmx",
            "allow_external_execute": True,
            "box_type": "dodecahedron",
            "box_padding_nm": 1.0,
            "salt_concentration_m": 0.15,
            "mdrun_args": "-ntmpi 1 -ntomp 4",
        }
    )

    assert protocol.max_bar_stderr_kcal_mol == 10.0
    assert protocol.equilibrate_retry_box_padding_nm == 2.0
    assert protocol.equilibrate_emergency_box_padding_nm == 5.0
    assert protocol.equilibrate_em_steps == 5000
    assert protocol.production_dt_ps == 0.002
    assert protocol.window_relax_em_steps == 500
    assert protocol.window_relax_md_ps == 0.1
    assert protocol.window_relax_md_dt_ps == 0.0005
    assert protocol.grompp_maxwarn_genion == 2
    assert protocol.grompp_maxwarn_equilibration == 2
    assert protocol.grompp_maxwarn_sampling == 2
    assert protocol.equilibrate_retry_box_padding_nm == 2.0
    assert protocol.equilibration_restraint_schedule == "legacy_posres"
    assert protocol.equilibration_release_npt_ps == 0
    assert protocol.equilibration_heavy_posres_fc_kj_mol_nm2 == 1000.0
    assert protocol.equilibration_backbone_posres_fc_kj_mol_nm2 == 250.0
    assert protocol.equilibration_pressure_coupling == "C-rescale"
    assert protocol.equilibration_pressure_tau_ps == 5.0
    assert protocol.equilibration_refcoord_scaling == "com"
    assert protocol.sampling_pressure_coupling == "C-rescale"
    assert protocol.sampling_pressure_tau_ps == 5.0
    assert protocol.sampling_refcoord_scaling == "all"


def test_hydrate_protocol_config_supports_validation_priority_single_point_preset() -> None:
    protocol = hydrate_protocol_config({"preset": "validation_priority_single_point"})

    assert protocol.preset == "validation_priority_single_point"
    assert protocol.lambda_windows == 8
    assert protocol.repeats == 3
    assert protocol.nvt_ps == 10
    assert protocol.npt_ps == 20
    assert protocol.production_ps == 20
    assert protocol.window_relax_em_steps == 1000
    assert protocol.window_relax_md_ps == 0.2
    assert protocol.nonbonded_cutoff_nm == 1.25
    assert protocol.vdw_switch_nm == 1.0


def test_hydrate_protocol_config_supports_validation_robust_single_point_preset() -> None:
    protocol = hydrate_protocol_config({"preset": "validation_robust_single_point"})

    assert protocol.preset == "validation_robust_single_point"
    assert protocol.lambda_windows == 24
    assert protocol.repeats == 3
    assert protocol.nvt_ps == 50
    assert protocol.npt_ps == 200
    assert protocol.production_ps == 500
    assert protocol.window_relax_em_steps == 2000
    assert protocol.window_relax_md_ps == 1.0
    assert protocol.nonbonded_cutoff_nm == 1.25
    assert protocol.vdw_switch_nm == 1.0


def test_build_batch_plan_honors_allow_charge_changing_protocol_flag(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: charge_flag_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    CHARGE FLAG DEMO",
                "ATOM      1  N   ASP H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  ASP H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   ASP H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   ASP H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A  10      20.000  10.000   8.000  1.00 20.00           N",
                "ATOM      6  CA  GLY A  10      21.000  10.500   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY A  10      22.000   9.500   9.000  1.00 20.00           C",
                "ATOM      8  O   GLY A  10      21.800   8.300   8.900  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "charge_shift,H,32,,D,N,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: validation_robust_single_point",
                "allow_charge_changing: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path=system_path,
        mutations_path=mutations_path,
        protocol_path=protocol_path,
        batch_id="charge_demo",
        runs_root=tmp_path / "runs",
    )

    assert len(batch_plan.jobs) == 1
    assert batch_plan.jobs[0].protocol.allow_charge_changing is True
    assert batch_plan.jobs[0].mutation_group.charge_conserving is False


def test_prepare_allows_sidechain_only_incomplete_standard_residues(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: sidechain_incomplete_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   GLU A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     14  CA  GLU A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     15  C   GLU A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     16  O   GLU A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     17  CB  GLU A  58      20.300  11.700  10.450  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="sidechain_prepare_ok",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "completed"
    assert "Applied PDBFixer atom repair" in statuses[-1].message

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    complex_leg = prepare_qc["legs"]["complex"]
    assert len(complex_leg["blocking_incomplete_standard_residues"]) == 0
    assert complex_leg["repair_summary"]["attempted"] is True
    assert complex_leg["repair_summary"]["available"] is True
    assert complex_leg["repair_summary"]["succeeded"] is True
    assert len(complex_leg["sidechain_only_incomplete_standard_residues"]) == 0
    complex_input = (job_dir / "legs" / "complex" / "input.pdb").read_text(encoding="utf-8")
    assert " CB  GLU A  58 " in complex_input


def test_prepare_repairs_backbone_incomplete_standard_residues(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: backbone_incomplete_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      5  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      6  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      7  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      8  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM      9  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     10  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     11  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     12  N   GLU A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     13  CA  GLU A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     14  C   GLU A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     15  O   GLU A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     16  CB  GLU A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     17  CG  GLU A  58      22.300   9.750   8.350  1.00 20.00           C",
                "ATOM     18  CD  GLU A  58      20.800  13.000  11.000  1.00 20.00           C",
                "ATOM     19  OE1 GLU A  58      21.100  13.900  11.800  1.00 20.00           O",
                "ATOM     20  OE2 GLU A  58      20.000  13.700  10.200  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="backbone_prepare_ok",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "completed"
    assert "Applied PDBFixer atom repair" in statuses[-1].message

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    complex_leg = prepare_qc["legs"]["complex"]
    apo_leg = prepare_qc["legs"]["apo"]
    assert len(complex_leg["blocking_incomplete_standard_residues"]) == 0
    assert len(apo_leg["blocking_incomplete_standard_residues"]) == 0
    assert complex_leg["repair_summary"]["trigger_residue_count"] == 1
    assert apo_leg["repair_summary"]["trigger_residue_count"] == 1
    assert complex_leg["repair_summary"]["succeeded"] is True
    assert apo_leg["repair_summary"]["succeeded"] is True
    complex_input = (job_dir / "legs" / "complex" / "input.pdb").read_text(encoding="utf-8")
    apo_input = (job_dir / "legs" / "apo" / "input.pdb").read_text(encoding="utf-8")
    assert " O   TYR H  32 " in complex_input
    assert " O   TYR H  32 " in apo_input


def test_prepare_repairs_sidechain_involving_same_residue_heavy_atom_clashes(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: geometry_clash_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   GLU H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  GLU H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   GLU H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   GLU H  32      14.400   8.500   8.600  1.00 20.00           O",
                "ATOM      5  CB  GLU H  32      12.400  11.500   9.200  1.00 20.00           C",
                "ATOM      6  CG  GLU H  32      14.300   8.550   8.550  1.00 20.00           C",
                "ATOM      7  CD  GLU H  32      14.300  12.800  10.000  1.00 20.00           C",
                "ATOM      8  OE1 GLU H  32      14.900  12.900  10.900  1.00 20.00           O",
                "ATOM      9  OE2 GLU H  32      13.900  13.600   9.300  1.00 20.00           O",
                "TER",
                "ATOM     10  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     11  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     12  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     13  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_e32d,H,32,,E,D,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="geometry_clash_repaired",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "completed"
    assert "Repaired sidechain-involving same-residue clashes" in statuses[-1].message

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    assert prepare_qc["legs"]["complex"]["repairable_sidechain_clashes"] != []
    assert prepare_qc["legs"]["apo"]["repairable_sidechain_clashes"] != []
    assert prepare_qc["legs"]["complex"]["clash_repair_summary"]["succeeded"] is True
    assert prepare_qc["legs"]["apo"]["clash_repair_summary"]["succeeded"] is True
    assert prepare_qc["legs"]["complex"]["blocking_intra_residue_heavy_atom_clashes"] == []
    assert prepare_qc["legs"]["apo"]["blocking_intra_residue_heavy_atom_clashes"] == []


def test_prepare_defers_persistent_repairable_sidechain_clashes_on_non_mutated_residues(
    tmp_path: Path, monkeypatch
) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: persistent_geometry_clash_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    original_complex = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   GLU A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     14  CA  GLU A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     15  C   GLU A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     16  O   GLU A   1      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     17  CB  GLU A   1      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     18  CG  GLU A   1      22.300   9.750   8.350  1.00 20.00           C",
                "ATOM     19  CD  GLU A   1      20.800  13.000  11.000  1.00 20.00           C",
                "ATOM     20  OE1 GLU A   1      21.100  13.900  11.800  1.00 20.00           O",
                "ATOM     21  OE2 GLU A   1      20.000  13.700  10.200  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n"
    )
    (tmp_path / "complex.pdb").write_text(original_complex, encoding="utf-8")

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    from abag_rbfe import stages as stages_module

    def _reintroduce_original_clash(input_path: Path, output_path: Path) -> dict[str, object]:
        if "/legs/complex/" in output_path.as_posix():
            output_path.write_text(original_complex, encoding="utf-8")
        else:
            output_path.write_text(input_path.read_text(encoding="utf-8"), encoding="utf-8")
        return {
            "attempted": True,
            "available": True,
            "succeeded": True,
        }

    monkeypatch.setattr(stages_module, "repair_missing_atoms_with_pdbfixer", _reintroduce_original_clash)

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="persistent_geometry_clash_deferred",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "completed"
    assert "Deferred 1 repairable clash-sidechains as stripped residues" in statuses[-1].message

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    complex_leg = prepare_qc["legs"]["complex"]
    assert complex_leg["repairable_sidechain_clashes"] != []
    assert complex_leg["clash_repair_summary"]["succeeded"] is True
    assert complex_leg["blocking_intra_residue_heavy_atom_clashes"] == []
    assert len(complex_leg["deferred_sidechain_clash_residues"]) == 1
    assert any(item["resname"] == "GLU" for item in complex_leg["sidechain_only_incomplete_standard_residues"])


def test_prepare_blocks_backbone_same_residue_heavy_atom_clashes(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: backbone_clash_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   GLY H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  GLY H  32      11.300  10.100   8.100  1.00 20.00           C",
                "ATOM      3  C   GLY H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   GLY H  32      13.350   9.250   9.700  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_g32a,H,32,,G,A,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="backbone_clash_blocked",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "blocked_input"
    assert "same-residue heavy-atom clashes" in statuses[-1].message

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    assert prepare_qc["legs"]["complex"]["repairable_sidechain_clashes"] == []
    assert len(prepare_qc["legs"]["complex"]["blocking_intra_residue_heavy_atom_clashes"]) == 1
    assert prepare_qc["legs"]["apo"]["repairable_sidechain_clashes"] == []
    assert len(prepare_qc["legs"]["apo"]["blocking_intra_residue_heavy_atom_clashes"]) == 1


def test_prepare_blocks_impossible_inter_residue_heavy_atom_clashes(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: inter_residue_clash_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
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

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_r32k,H,32,,R,K,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="inter_residue_clash_blocked",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "blocked_input"
    assert "impossible inter-residue heavy-atom clashes" in statuses[-1].message

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    assert len(prepare_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"]) == 1
    assert prepare_qc["legs"]["apo"]["inter_residue_heavy_atom_clashes"] == []
    assert prepare_qc["legs"]["complex"]["blocking_intra_residue_heavy_atom_clashes"] == []
    assert prepare_qc["legs"]["apo"]["blocking_intra_residue_heavy_atom_clashes"] == []


def test_mutate_blocks_impossible_inter_residue_heavy_atom_clashes_in_mutant_pdb(
    tmp_path: Path, monkeypatch
) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: mutate_inter_residue_clash_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "ATOM     12  N   ILE A   4      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM     13  CA  ILE A   4      30.900  10.600   8.900  1.00 20.00           C",
                "ATOM     14  C   ILE A   4      32.100   9.700   9.100  1.00 20.00           C",
                "ATOM     15  O   ILE A   4      32.500   8.900   8.300  1.00 20.00           O",
                "ATOM     16  CB  ILE A   4      30.300  11.000  10.300  1.00 20.00           C",
                "ATOM     17  CG1 ILE A   4      31.000  12.200  11.000  1.00 20.00           C",
                "ATOM     18  CG2 ILE A   4      28.800  11.300  10.200  1.00 20.00           C",
                "ATOM     19  CD1 ILE A   4      31.600  12.500  12.350  1.00 20.00           C",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_r32k,H,32,,R,K,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="mutate_inter_residue_clash_blocked",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)

    complex_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
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
        + "\n"
    )
    apo_mutant = (
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   ARG H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  ARG H  32      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   ARG H  32      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   ARG H  32      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  ARG H  32      12.100  11.500   9.000  1.00 20.00           C",
                "ATOM      6  CG  ARG H  32      12.900  12.600   8.300  1.00 20.00           C",
                "ATOM      7  CD  ARG H  32      13.400  13.200   9.600  1.00 20.00           C",
                "ATOM      8  NE  ARG H  32      14.600  12.400  10.000  1.00 20.00           N",
                "ATOM      9  CZ  ARG H  32      15.000  12.300  11.200  1.00 20.00           C",
                "ATOM     10  NH1 ARG H  32      14.400  12.900  12.200  1.00 20.00           N",
                "ATOM     11  NH2 ARG H  32      15.900  11.500  11.300  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n"
    )

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        stage = script_path.stem
        if stage != "mutate":
            return CommandOutcome("planned", f"{stage} commands written")
        for leg, mutant_pdb_text in (("complex", complex_mutant), ("apo", apo_mutant)):
            pmx_dir = job_dir / "legs" / leg / "pmx"
            pmx_dir.mkdir(parents=True, exist_ok=True)
            (pmx_dir / "mutant.pdb").write_text(mutant_pdb_text, encoding="utf-8")
            _write_valid_mock_gro(pmx_dir / "processed.gro")
            (pmx_dir / "topol.top").write_text("; generated\n", encoding="utf-8")
            (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")
        return CommandOutcome("completed", "External commands completed. Log: fake.log")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)
    statuses = run_job(job_dir, execute=False, to_stage="mutate")

    assert [status.stage for status in statuses] == ["ingest", "prepare", "mutate"]
    assert statuses[-1].state == "blocked_input"
    assert "Mutated structure contains impossible inter-residue heavy-atom clashes" in statuses[-1].message

    mutate_qc = json.loads((job_dir / "artifacts" / "mutate_qc.json").read_text(encoding="utf-8"))
    assert len(mutate_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"]) == 1
    clash = mutate_qc["legs"]["complex"]["inter_residue_heavy_atom_clashes"][0]
    assert {clash["resname"], clash["partner_resname"]} == {"ARG", "ILE"}
    assert mutate_qc["legs"]["apo"]["inter_residue_heavy_atom_clashes"] == []


def test_prepare_ignores_altloc_duplicates_when_generating_geometry_qc(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: altloc_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   LEU H  32      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA ALEU H  32      12.200  10.100   8.700  0.40 20.00           C",
                "ATOM      3  CA BLEU H  32      12.320  10.180   8.780  0.60 20.00           C",
                "ATOM      4  C   LEU H  32      13.410   9.280   8.960  1.00 20.00           C",
                "ATOM      5  O   LEU H  32      13.640   8.160   8.520  1.00 20.00           O",
                "ATOM      6  CB ALEU H  32      11.860  11.520   9.090  0.40 20.00           C",
                "ATOM      7  CB BLEU H  32      12.020  11.610   9.180  0.60 20.00           C",
                "ATOM      8  CG ALEU H  32      12.620  12.660   8.430  0.40 20.00           C",
                "ATOM      9  CG BLEU H  32      12.790  12.760   8.540  0.60 20.00           C",
                "ATOM     10 CD1ALEU H  32      12.330  14.040   9.000  0.40 20.00           C",
                "ATOM     11 CD1BLEU H  32      12.500  14.140   9.100  0.60 20.00           C",
                "ATOM     12 CD2ALEU H  32      14.070  12.450   8.640  0.40 20.00           C",
                "ATOM     13 CD2BLEU H  32      14.250  12.550   8.740  0.60 20.00           C",
                "TER",
                "ATOM     14  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     15  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     16  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     17  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_l32i,H,32,,L,I,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="altloc_prepare_ok",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert [status.stage for status in statuses] == ["ingest", "prepare"]
    assert statuses[-1].state == "completed"

    prepare_qc = json.loads((job_dir / "artifacts" / "prepare_qc.json").read_text(encoding="utf-8"))
    assert prepare_qc["legs"]["complex"]["intra_residue_heavy_atom_clashes"] == []
    assert prepare_qc["legs"]["apo"]["intra_residue_heavy_atom_clashes"] == []


def test_run_job_skips_duplicate_execution_when_live_lock_exists(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: duplicate_lock_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="duplicate_lock_demo",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)

    sleeper = subprocess.Popen(["sleep", "60"])
    try:
        (job_dir / ".abag_job_execution_lock.json").write_text(
            json.dumps(
                {
                    "pid": sleeper.pid,
                    "started_at": "2026-06-08T00:00:00Z",
                    "from_stage": "ingest",
                    "to_stage": "prepare",
                    "execute": False,
                }
            )
            + "\n",
            encoding="utf-8",
        )

        statuses = run_job(job_dir, execute=False, to_stage="prepare")

        assert len(statuses) == 1
        assert statuses[0].stage == "ingest"
        assert statuses[0].state == "running"
        assert "already in progress" in statuses[0].message
        assert not any((job_dir / "stages").glob("*.json"))
    finally:
        sleeper.terminate()
        sleeper.wait(timeout=5)


def test_job_active_process_detects_live_job_script(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: active_process_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="active_process_demo",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    script = job_dir / "artifacts" / "commands" / "hold.sh"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text("#!/usr/bin/env bash\nwhile true; do sleep 60; done\n", encoding="utf-8")
    script.chmod(0o755)

    from abag_rbfe import stages as stages_module

    sleeper = subprocess.Popen(["bash", str(script)])
    try:
        active_process = None
        for _ in range(30):
            active_process = stages_module._job_active_process(job_dir)
            if active_process is not None:
                break
            time.sleep(0.1)

        assert active_process is not None
        assert active_process[0] == sleeper.pid
        assert str(script) in active_process[1]
    finally:
        sleeper.terminate()
        sleeper.wait(timeout=5)


def test_run_job_skips_duplicate_execution_when_live_job_process_exists(tmp_path: Path, monkeypatch) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: duplicate_process_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "TER",
                "ATOM      5  N   GLY A   1      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM      6  CA  GLY A   1      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM      7  C   GLY A   1      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM      8  O   GLY A   1      22.400   9.700   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="duplicate_process_demo",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)

    from abag_rbfe import stages as stages_module

    live_command = job_dir / "artifacts" / "commands" / "sample.sh"
    monkeypatch.setattr(
        stages_module,
        "_job_active_process",
        lambda *_args, **_kwargs: (43210, f"bash {live_command}"),
    )

    statuses = run_job(job_dir, execute=False, to_stage="prepare")

    assert len(statuses) == 1
    assert statuses[0].stage == "ingest"
    assert statuses[0].state == "running"
    assert "already in progress" in statuses[0].message
    assert "43210" in statuses[0].message
    assert str(live_command) in statuses[0].message
    assert not any((job_dir / "stages").glob("*.json"))


def test_stage_scripts_export_deterministic_gpu_assignment(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.setenv("ABAG_RBFE_VISIBLE_GPUS", "2,5")
    discover_visible_gpu_devices.cache_clear()

    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: gpu_env_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H, L]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   TYR H  32      11.104  13.207   8.154  1.00 20.00           N",
                "ATOM      2  CA  TYR H  32      12.560  13.311   8.321  1.00 20.00           C",
                "ATOM      3  C   TYR H  32      13.150  12.090   8.991  1.00 20.00           C",
                "ATOM      4  O   TYR H  32      12.811  10.945   8.710  1.00 20.00           O",
                "ATOM      5  CB  TYR H  32      12.936  14.622   9.024  1.00 20.00           C",
                "ATOM      6  CG  TYR H  32      12.086  15.074  10.191  1.00 20.00           C",
                "ATOM      7  CD1 TYR H  32      10.705  14.949  10.172  1.00 20.00           C",
                "ATOM      8  CD2 TYR H  32      12.671  15.624  11.341  1.00 20.00           C",
                "ATOM      9  CE1 TYR H  32       9.930  15.355  11.245  1.00 20.00           C",
                "ATOM     10  CE2 TYR H  32      11.904  16.038  12.414  1.00 20.00           C",
                "ATOM     11  CZ  TYR H  32      10.534  15.897  12.370  1.00 20.00           C",
                "ATOM     12  OH  TYR H  32       9.773  16.302  13.440  1.00 20.00           O",
                "TER",
                "ATOM     13  N   SER L   1      14.100  10.200   8.900  1.00 20.00           N",
                "ATOM     14  CA  SER L   1      15.200  10.800   9.200  1.00 20.00           C",
                "ATOM     15  C   SER L   1      16.200   9.800   9.800  1.00 20.00           C",
                "ATOM     16  O   SER L   1      16.000   8.600   9.700  1.00 20.00           O",
                "ATOM     17  CB  SER L   1      14.700  11.900  10.100  1.00 20.00           C",
                "ATOM     18  OG  SER L   1      15.600  12.900  10.400  1.00 20.00           O",
                "TER",
                "ATOM     19  N   LYS A  58      20.100  11.000   8.300  1.00 20.00           N",
                "ATOM     20  CA  LYS A  58      21.000  11.500   9.100  1.00 20.00           C",
                "ATOM     21  C   LYS A  58      22.200  10.600   9.100  1.00 20.00           C",
                "ATOM     22  O   LYS A  58      22.400   9.700   8.300  1.00 20.00           O",
                "ATOM     23  CB  LYS A  58      20.300  11.700  10.450  1.00 20.00           C",
                "ATOM     24  CG  LYS A  58      19.100  12.650  10.400  1.00 20.00           C",
                "ATOM     25  CD  LYS A  58      18.400  12.780  11.760  1.00 20.00           C",
                "ATOM     26  CE  LYS A  58      17.200  13.720  11.700  1.00 20.00           C",
                "ATOM     27  NZ  LYS A  58      16.500  13.840  13.020  1.00 20.00           N",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_y32f,H,32,,Y,F,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text(
        "\n".join(
            [
                "preset: single_point",
                "force_field: amber99sb-star-ildn-mut",
                "water_model: tip3p",
                "lambda_windows: 2",
                "repeats: 1",
                "nvt_ps: 1",
                "npt_ps: 1",
                "production_ps: 1",
                "allow_external_execute: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="gpu_env_batch",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)
    run_job(job_dir, execute=False, to_stage="sample")

    mutate_script = (job_dir / "artifacts" / "commands" / "mutate.sh").read_text(encoding="utf-8")
    equilibrate_script = (job_dir / "artifacts" / "commands" / "equilibrate.sh").read_text(encoding="utf-8")
    sample_script = (job_dir / "artifacts" / "commands" / "sample.sh").read_text(encoding="utf-8")

    assert "export CUDA_VISIBLE_DEVICES=" in mutate_script
    assert "export CUDA_VISIBLE_DEVICES=" in equilibrate_script
    assert "export CUDA_VISIBLE_DEVICES=" in sample_script
    assert "export GMX_MAXBACKUP=-1" in mutate_script
    assert "export GMX_MAXBACKUP=-1" in equilibrate_script
    assert "export GMX_MAXBACKUP=-1" in sample_script
    assert ("export CUDA_VISIBLE_DEVICES=2" in sample_script) or ("export CUDA_VISIBLE_DEVICES=5" in sample_script)

    discover_visible_gpu_devices.cache_clear()


def test_mutate_blocks_broken_hybrid_residue_integrity(tmp_path: Path, monkeypatch) -> None:
    """Mutate must block when the hybrid residue loses atoms vs mutant.pdb or
    its A/B-state charges are not integer (regression guard for the pdb2gmx
    -ignh heavy-atom-only hybrid defect)."""
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: hybrid_integrity_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text(
        "\n".join(
            [
                "HEADER    DEMO ABAG",
                "ATOM      1  N   GLN H  89      11.000  10.000   8.000  1.00 20.00           N",
                "ATOM      2  CA  GLN H  89      12.200  10.100   8.700  1.00 20.00           C",
                "ATOM      3  C   GLN H  89      13.300   9.200   8.900  1.00 20.00           C",
                "ATOM      4  O   GLN H  89      14.400   8.700   8.200  1.00 20.00           O",
                "ATOM      5  CB  GLN H  89      12.400  11.500   9.200  1.00 20.00           C",
                "TER",
                "ATOM      6  N   ILE A   4      30.000  10.000   8.000  1.00 20.00           N",
                "ATOM      7  CA  ILE A   4      30.900  10.600   8.900  1.00 20.00           C",
                "ATOM      8  C   ILE A   4      32.100   9.700   9.100  1.00 20.00           C",
                "ATOM      9  O   ILE A   4      32.500   8.900   8.300  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side\nsingle_h_q89a,H,89,,Q,A,antibody\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: single_point\n", encoding="utf-8")
    batch_plan = build_batch_plan(
        system_path,
        mutations_path,
        protocol_path,
        batch_id="mutate_hybrid_integrity_blocked",
        runs_root=tmp_path / "runs",
    )
    job_dir = Path(batch_plan.jobs[0].workdir)

    hybrid_pdb_atoms = [
        ("N", "N"), ("H", "H"), ("CA", "C"), ("HA", "H"), ("CB", "C"),
        ("HB1", "H"), ("HB2", "H"), ("CG", "C"), ("HG1", "H"), ("HG2", "H"),
        ("CD", "C"), ("OE1", "O"), ("NE2", "N"), ("HE21", "H"), ("HE22", "H"),
        ("C", "C"), ("O", "O"), ("HV", "H"),
    ]
    mutant_lines = ["HEADER    DEMO ABAG"]
    for idx, (name, element) in enumerate(hybrid_pdb_atoms, start=1):
        mutant_lines.append(
            f"ATOM  {idx:5d}  {name:<4s}Q2A H  89    {11.0 + idx * 0.1:8.3f}{10.0:8.3f}{8.0:8.3f}  1.00 20.00           {element}"
        )
    mutant_lines.extend(["TER", "END"])
    mutant_text = "\n".join(mutant_lines) + "\n"

    broken_itp_lines = ["[ atoms ]"]
    for nr, (name, _element) in enumerate(
        [a for a in hybrid_pdb_atoms if a[1] != "H"], start=1
    ):
        broken_itp_lines.append(
            f"{nr:6d}     CT     89    Q2A {name:>6s} {nr:6d}   -0.1500    12.0100"
        )

    def fake_run_script(
        self: CommandRunner,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        stage = script_path.stem
        if stage != "mutate":
            return CommandOutcome("planned", f"{stage} commands written")
        for leg in ("complex", "apo"):
            pmx_dir = job_dir / "legs" / leg / "pmx"
            pmx_dir.mkdir(parents=True, exist_ok=True)
            (pmx_dir / "mutant.pdb").write_text(mutant_text, encoding="utf-8")
            (pmx_dir / "pmx_topol_Protein_chain_H.itp").write_text(
                "\n".join(broken_itp_lines) + "\n", encoding="utf-8"
            )
            _write_valid_mock_gro(pmx_dir / "processed.gro")
            (pmx_dir / "topol.top").write_text("; generated\n", encoding="utf-8")
            (pmx_dir / "pmxtop.top").write_text("; generated\n", encoding="utf-8")
        return CommandOutcome("completed", "External commands completed. Log: fake.log")

    monkeypatch.setattr(CommandRunner, "run_script", fake_run_script)
    statuses = run_job(job_dir, execute=False, to_stage="mutate")

    assert [status.stage for status in statuses] == ["ingest", "prepare", "mutate"]
    assert statuses[-1].state == "blocked_input"
    assert "Hybrid residue integrity check failed" in statuses[-1].message
    assert "atoms" in statuses[-1].message


def test_adaptive_lambda_windows_raises_floor_by_mutation_size() -> None:
    from dataclasses import dataclass

    from abag_rbfe.planning import adaptive_lambda_windows

    @dataclass
    class Site:
        wt: str
        mut: str

    assert adaptive_lambda_windows(8, (Site("S", "A"),)) == 12
    assert adaptive_lambda_windows(8, (Site("Y", "A"),)) == 16
    assert adaptive_lambda_windows(8, (Site("A", "W"),)) == 16
    assert adaptive_lambda_windows(8, (Site("Y", "F"),)) == 12  # aromatic->aromatic is conservative
    assert adaptive_lambda_windows(24, (Site("Y", "A"),)) == 24  # never lowers explicit depth


def test_build_batch_plan_applies_adaptive_lambda_unless_pinned(tmp_path: Path) -> None:
    system_path = tmp_path / "system.yml"
    system_path.write_text(
        "\n".join(
            [
                "system_name: adaptive_lambda_demo",
                f"input_structure: {tmp_path / 'complex.pdb'}",
                "structure_source: experimental",
                "antibody_chains: [H]",
                "antigen_chains: [A]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (tmp_path / "complex.pdb").write_text("HEADER DEMO\nEND\n", encoding="utf-8")
    mutations_path = tmp_path / "mutations.csv"
    mutations_path.write_text(
        "\n".join(
            [
                "mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side",
                "single_h_s57a,H,57,,S,A,antibody",
                "single_h_y32a,H,32,,Y,A,antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    protocol_path = tmp_path / "protocol.yml"
    protocol_path.write_text("preset: validation_priority_single_point\n", encoding="utf-8")

    plan = build_batch_plan(
        system_path, mutations_path, protocol_path,
        batch_id="adaptive_lambda_demo", runs_root=tmp_path / "runs",
    )
    by_group = {job.mutation_group.mutation_group_id: job for job in plan.jobs}
    assert by_group["single_h_s57a"].protocol.lambda_windows == 12
    assert by_group["single_h_y32a"].protocol.lambda_windows == 16

    pinned_path = tmp_path / "protocol_pinned.yml"
    pinned_path.write_text(
        "preset: validation_priority_single_point\nlambda_windows: 8\n", encoding="utf-8"
    )
    plan_pinned = build_batch_plan(
        system_path, mutations_path, pinned_path,
        batch_id="adaptive_lambda_pinned", runs_root=tmp_path / "runs2",
    )
    assert all(job.protocol.lambda_windows == 8 for job in plan_pinned.jobs)
