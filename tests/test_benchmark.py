from math import exp, isclose, log1p
from pathlib import Path

import abag_rbfe.benchmark as benchmark_module
from abag_rbfe.benchmark import (
    calibrate_ab_bind_plan,
    curate_ab_bind,
    materialize_ab_bind_inputs,
    plan_ab_bind_batches,
    plan_ab_bind_rescues,
    report_ab_bind_plan,
    run_ab_bind_plan,
)
from abag_rbfe.io_utils import ensure_dir, read_csv_rows, read_json, read_yaml, write_csv_rows, write_json, write_yaml
from abag_rbfe.models import StageStatus


def _histogram_blocks(mode: str = "good") -> list[list[tuple[float, float]]]:
    if mode == "poor":
        forward = [(10.0, 1.0), (11.0, 1.0), (12.0, 1.0), (13.0, 1.0)]
        reverse = [(100.0, 1.0), (101.0, 1.0), (102.0, 1.0), (103.0, 1.0)]
    else:
        forward = [(10.0, 1.0), (11.0, 2.0), (12.0, 2.0), (13.0, 1.0)]
        reverse = [(11.0, 1.0), (12.0, 2.0), (13.0, 2.0), (14.0, 1.0)]
    return [
        [(1.0, 1.0), (2.0, 1.0), (3.0, 0.0), (4.0, 0.0)],
        [(0.0, 0.0), (0.0, 0.0), (0.0, 4.0), (0.0, 4.0)],
        forward,
        [(20.0, 1.0), (21.0, 1.0), (22.0, 0.0), (23.0, 0.0)],
        reverse,
        [(0.0, 0.0), (0.0, 0.0), (0.0, 4.0), (0.0, 4.0)],
    ]


def _write_bar_outputs(repeat_dir: Path, delta_kt: float, stderr_kt: float, *, overlap_mode: str = "good") -> None:
    ensure_dir(repeat_dir / "bar")
    ensure_dir(repeat_dir / "lambda_000")
    ensure_dir(repeat_dir / "lambda_001")
    (repeat_dir / "lambda_000" / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")
    (repeat_dir / "lambda_001" / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")
    (repeat_dir / "bar" / "bar.xvg").write_text(f"0.500 {delta_kt:.6f} {stderr_kt:.6f}\n", encoding="utf-8")
    (repeat_dir / "bar" / "barint.xvg").write_text(f"0 0.000000\n1 {delta_kt:.6f}\n", encoding="utf-8")
    histogram_lines = ["@ s0 legend \"mock\"", "@ s1 legend \"mock\"", "@ s2 legend \"mock\"", "@ s3 legend \"mock\"", "@ s4 legend \"mock\"", "@ s5 legend \"mock\""]
    for block in _histogram_blocks(overlap_mode):
        for x, y in block:
            histogram_lines.append(f"{x} {y}")
    (repeat_dir / "bar" / "histogram.xvg").write_text("\n".join(histogram_lines) + "\n", encoding="utf-8")


def _write_stage(job_dir: Path, stage: str) -> None:
    ensure_dir(job_dir / "stages")
    write_json(
        job_dir / "stages" / f"{stage}.json",
        {
            "stage": stage,
            "state": "completed",
            "message": f"{stage} completed",
            "commands": [],
            "artifacts": [],
        },
    )


def _write_running_stage(
    job_dir: Path,
    stage: str,
    *,
    commands: list[str] | None = None,
    artifacts: list[str] | None = None,
) -> None:
    ensure_dir(job_dir / "stages")
    write_json(
        job_dir / "stages" / f"{stage}.json",
        {
            "stage": stage,
            "state": "running",
            "message": "Stage execution started.",
            "commands": commands or [],
            "artifacts": artifacts or [],
        },
    )


def _write_sample_window_completed(lambda_dir: Path) -> None:
    ensure_dir(lambda_dir)
    (lambda_dir / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")
    (lambda_dir / "md.gro").write_text("mock\n", encoding="utf-8")


def _write_sample_window_started(lambda_dir: Path) -> None:
    ensure_dir(lambda_dir)
    (lambda_dir / "topol.tpr").write_text("mock\n", encoding="utf-8")


def _write_equilibrate_repeat_completed(repeat_dir: Path) -> None:
    ensure_dir(repeat_dir / "equilibration")
    (repeat_dir / "system.top").write_text("mock\n", encoding="utf-8")
    (repeat_dir / "equilibration" / "npt.gro").write_text("mock\n", encoding="utf-8")


def _write_equilibrate_repeat_started(repeat_dir: Path) -> None:
    ensure_dir(repeat_dir / "setup")
    (repeat_dir / "system.top").write_text("mock\n", encoding="utf-8")
    (repeat_dir / "setup" / "genion.tpr").write_text("mock\n", encoding="utf-8")


def test_validation_failure_taxonomy_classifies_core_categories() -> None:
    cases = [
        ({"benchmark_qc_qualified": True, "diagnostic_code": "qc_pass"}, ("benchmark_qc_qualified", "qc_pass")),
        (
            {"diagnostic_family": "qc", "diagnostic_code": "qc_low_overlap", "benchmark_qc_qualified": False},
            ("qc_sampling_issue", "qc_low_overlap"),
        ),
        (
            {"diagnostic_family": "input_structure", "diagnostic_code": "prepare_blocked_input"},
            ("input_structure_issue", "prepare_blocked_input"),
        ),
        (
            {
                "diagnostic_family": "mutate_setup",
                "diagnostic_code": "running_mutate",
                "current_invalid_mutate_output": True,
                "current_invalid_mutate_output_code": "mutate_processed_gro_invalid",
            },
            ("mutate_setup_issue", "mutate_processed_gro_invalid"),
        ),
        (
            {
                "diagnostic_family": "running",
                "diagnostic_code": "running_sample",
                "latest_stage_state": "running",
                "current_invalid_mutate_output": True,
                "current_invalid_mutate_output_code": "mutate_processed_gro_isolated_residue_hydrogen",
            },
            ("running_execution", "running_sample"),
        ),
        ({"diagnostic_code": "sample_missing_npt_gro"}, ("downstream_input_issue", "sample_missing_npt_gro")),
        ({"diagnostic_code": "sample_gmxlib_unavailable"}, ("external_dependency_issue", "sample_gmxlib_unavailable")),
        ({"diagnostic_code": "sample_failed"}, ("execution_failure", "sample_failed")),
        ({"diagnostic_family": "running", "diagnostic_code": "running_sample"}, ("running_execution", "running_sample")),
        ({"diagnostic_family": "pending", "diagnostic_code": "pending_equilibrate"}, ("pending_execution", "pending_equilibrate")),
        ({"diagnostic_family": "completed", "diagnostic_code": "reported"}, ("completed_unqualified", "reported")),
    ]

    for row, expected in cases:
        assert benchmark_module._validation_failure_taxonomy_category(row) == expected


def test_benchmark_target_metrics_bundle_flags_systematically_poor_targets() -> None:
    pair_rows = [
        {
            "complex_id": "BAD1",
            "predicted_ddg_kcal_mol": 0.0,
            "experimental_ddg_kcal_mol": 5.0,
            "ddg_error_kcal_mol": -5.0,
            "abs_error_kcal_mol": 5.0,
        },
        {
            "complex_id": "BAD1",
            "predicted_ddg_kcal_mol": 1.0,
            "experimental_ddg_kcal_mol": 6.0,
            "ddg_error_kcal_mol": -5.0,
            "abs_error_kcal_mol": 5.0,
        },
        {
            "complex_id": "BAD1",
            "predicted_ddg_kcal_mol": 2.0,
            "experimental_ddg_kcal_mol": 7.0,
            "ddg_error_kcal_mol": -5.0,
            "abs_error_kcal_mol": 5.0,
        },
        {
            "complex_id": "BAD1",
            "predicted_ddg_kcal_mol": 3.0,
            "experimental_ddg_kcal_mol": 8.0,
            "ddg_error_kcal_mol": -5.0,
            "abs_error_kcal_mol": 5.0,
        },
        {
            "complex_id": "GOOD1",
            "predicted_ddg_kcal_mol": 1.0,
            "experimental_ddg_kcal_mol": 1.0,
            "ddg_error_kcal_mol": 0.0,
            "abs_error_kcal_mol": 0.0,
        },
        {
            "complex_id": "GOOD1",
            "predicted_ddg_kcal_mol": 2.0,
            "experimental_ddg_kcal_mol": 2.0,
            "ddg_error_kcal_mol": 0.0,
            "abs_error_kcal_mol": 0.0,
        },
        {
            "complex_id": "GOOD1",
            "predicted_ddg_kcal_mol": 3.0,
            "experimental_ddg_kcal_mol": 3.0,
            "ddg_error_kcal_mol": 0.0,
            "abs_error_kcal_mol": 0.0,
        },
        {
            "complex_id": "GOOD1",
            "predicted_ddg_kcal_mol": 4.0,
            "experimental_ddg_kcal_mol": 4.0,
            "ddg_error_kcal_mol": 0.0,
            "abs_error_kcal_mol": 0.0,
        },
    ]

    bundle = benchmark_module._benchmark_target_metrics_bundle(pair_rows)

    assert bundle["excluded_target_ids"] == ["BAD1"]
    assert bundle["filtered_metrics"]["paired_job_count"] == 4
    assert isclose(bundle["filtered_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    metrics_by_complex = {row["complex_id"]: row for row in bundle["target_metrics"]}
    assert metrics_by_complex["BAD1"]["systematically_poor_target"] is True
    assert metrics_by_complex["BAD1"]["all_pairs_above_abs_error_threshold"] is True
    assert isclose(metrics_by_complex["BAD1"]["leave_one_out_pearson_r"], 1.0, rel_tol=1e-9)
    assert metrics_by_complex["GOOD1"]["systematically_poor_target"] is False
    assert metrics_by_complex["GOOD1"]["excluded_from_target_filtered_metrics"] is False


def test_benchmark_target_metrics_bundle_iteratively_excludes_targets_that_hurt_overall_correlation() -> None:
    pair_rows = []
    for complex_id, predicted_values in (
        ("GOOD1", [1.0, 2.0, 3.0, 4.0]),
        ("GOOD2", [1.5, 2.5, 3.5, 4.5]),
        ("BAD1", [1.5, 1.0, 0.5, 0.0]),
        ("BAD2", [2.0, 1.5, 1.0, 0.5]),
    ):
        for experimental_ddg, predicted_ddg in enumerate(predicted_values, start=1):
            error = predicted_ddg - float(experimental_ddg)
            pair_rows.append(
                {
                    "complex_id": complex_id,
                    "predicted_ddg_kcal_mol": predicted_ddg,
                    "experimental_ddg_kcal_mol": float(experimental_ddg),
                    "ddg_error_kcal_mol": error,
                    "abs_error_kcal_mol": abs(error),
                }
            )

    bundle = benchmark_module._benchmark_target_metrics_bundle(pair_rows)

    assert bundle["excluded_target_ids"] == ["BAD1", "BAD2"]
    assert bundle["filtered_metrics"]["paired_job_count"] == 8
    assert isclose(bundle["filtered_metrics"]["pearson_r"], 0.9759000729485332, rel_tol=1e-9)
    metrics_by_complex = {row["complex_id"]: row for row in bundle["target_metrics"]}
    assert metrics_by_complex["BAD1"]["systematically_poor_target"] is True
    assert metrics_by_complex["BAD1"]["systematically_poor_abs_error_target"] is False
    assert metrics_by_complex["BAD1"]["systematically_poor_correlation_target"] is True
    assert metrics_by_complex["BAD1"]["target_exclusion_iteration"] == 1
    assert metrics_by_complex["BAD1"]["systematically_poor_target_reason"] == "iterative_leave_one_out_gain"
    assert isclose(metrics_by_complex["BAD1"]["overall_pearson_r_at_exclusion"], 0.22032632461961588, rel_tol=1e-9)
    assert isclose(metrics_by_complex["BAD1"]["leave_one_out_pearson_gain_at_exclusion"], 0.23931223515730576, rel_tol=1e-9)
    assert metrics_by_complex["BAD2"]["systematically_poor_target"] is True
    assert metrics_by_complex["BAD2"]["systematically_poor_abs_error_target"] is False
    assert metrics_by_complex["BAD2"]["systematically_poor_correlation_target"] is True
    assert metrics_by_complex["BAD2"]["target_exclusion_iteration"] == 2
    assert isclose(metrics_by_complex["BAD2"]["leave_one_out_pearson_r_at_exclusion"], 0.9759000729485332, rel_tol=1e-9)
    assert isclose(metrics_by_complex["GOOD1"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert metrics_by_complex["GOOD1"]["systematically_poor_target"] is False


def test_ab_bind_curation_layers(tmp_path: Path) -> None:
    source = tmp_path / "ab_bind_source.csv"
    source.write_text(
        "\n".join(
            [
                "row_id,complex_id,complex_class,structure_source,ddg_kcal_mol,structure_mappable,stable_run_candidate,mutation_tokens",
                "1,1ABC,antibody-antigen,experimental,0.8,true,true,H:Y32F@antibody",
                "2,2BCD,antibody-antigen,experimental,1.1,true,true,H:S52T@antibody;H:N54Q@antibody",
                "3,3CDE,antibody-antigen,experimental,1.4,true,true,H:Y32F@antibody;A:T52S@antigen",
                "4,4DEF,antibody-antigen,homology_model,0.4,true,true,H:S31A@antibody",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    curate_ab_bind(source, tmp_path)
    core_v1 = read_csv_rows(tmp_path / "curated" / "ab_bind_rbfe_core_v1.csv")
    core_v2 = read_csv_rows(tmp_path / "curated" / "ab_bind_rbfe_core_v2.csv")
    summary = read_json(tmp_path / "summary.json")

    assert len(core_v1) == 1
    assert len(core_v2) == 2
    assert summary["core_v1"]["row_count"] == 1
    assert summary["core_v2"]["row_count"] == 2


def test_ab_bind_raw_schema_curation_with_annotations(tmp_path: Path) -> None:
    source = tmp_path / "AB-Bind_experimental_data.csv"
    source.write_text(
        "\n".join(
            [
                "#PDB,Partners(A_B),Protein-1,Protein-2,Mutation,ddG(kcal/mol)",
                "1VFB,C_HL,IGG1-Kappa D1.3 FV,Lysozyme,H:Y32F,0.8",
                "1VFB,C_HL,IGG1-Kappa D1.3 FV,Lysozyme,\"H:Y32F,H:V34I\",1.1",
                "1AK4,A_D,huCyc-A,HIV-1 CAPSID (N-Term),D:A488G,2.4",
                "HM_1YY9,A_HL,EGFR,hu225,H:S31A,0.4",
            ]
        )
        + "\n",
        encoding="latin-1",
    )

    annotations = tmp_path / "ab_bind_complex_annotations.csv"
    annotations.write_text(
        "\n".join(
            [
                "complex_id,complex_class,structure_source,antibody_chains,antigen_chains,structure_mappable,stable_run_candidate",
                "1VFB,antibody-antigen,experimental,HL,C,true,true",
                "1AK4,other,experimental,,,true,true",
                "HM_1YY9,antibody-antigen,homology_model,HL,A,true,true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    curate_ab_bind(source, tmp_path, annotations_path=annotations)
    registered = read_csv_rows(tmp_path / "curated" / "ab_bind_source_registered.csv")
    core_v1 = read_csv_rows(tmp_path / "curated" / "ab_bind_rbfe_core_v1.csv")
    core_v2 = read_csv_rows(tmp_path / "curated" / "ab_bind_rbfe_core_v2.csv")
    summary = read_json(tmp_path / "summary.json")

    assert len(registered) == 4
    assert len(core_v1) == 1
    assert len(core_v2) == 2
    assert core_v1[0]["mutation_tokens"] == "H:Y32F@antibody"
    assert core_v2[1]["mutation_tokens"] == "H:Y32F@antibody;H:V34I@antibody"
    assert "non_antibody_antigen" in registered[2]["exclusion_codes"]
    assert "non_experimental_structure" in registered[3]["exclusion_codes"]
    assert summary["source"]["row_count"] == 4
    assert summary["source"]["complex_count"] == 3
    assert summary["exclusion_counts"]["non_antibody_antigen"] == 1


def test_materialize_ab_bind_inputs_from_curated_rows(tmp_path: Path) -> None:
    source = tmp_path / "AB-Bind_experimental_data.csv"
    source.write_text(
        "\n".join(
            [
                "#PDB,Partners(A_B),Protein-1,Protein-2,Mutation,ddG(kcal/mol)",
                "1VFB,C_HL,IGG1-Kappa D1.3 FV,Lysozyme,H:Y32F,0.8",
                "1VFB,C_HL,IGG1-Kappa D1.3 FV,Lysozyme,\"H:Y32F,H:V34I\",1.1",
            ]
        )
        + "\n",
        encoding="latin-1",
    )

    annotations = tmp_path / "ab_bind_complex_annotations.csv"
    annotations.write_text(
        "\n".join(
            [
                "complex_id,complex_class,structure_source,antibody_chains,antigen_chains,structure_mappable,stable_run_candidate",
                "1VFB,antibody-antigen,experimental,HL,C,true,true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    source_dir = tmp_path / "source"
    source_dir.mkdir()
    (source_dir / "1VFB.pdb").write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    (source_dir / "ab_bind_complex_annotations.csv").write_text(annotations.read_text(encoding="utf-8"), encoding="utf-8")

    curate_ab_bind(source, tmp_path, annotations_path=annotations)
    summary = materialize_ab_bind_inputs(tmp_path, annotations_path=annotations)

    assert summary["generated"]["core_v1"]["complex_count"] == 1
    assert summary["generated"]["core_v2"]["complex_count"] == 1

    system_yml = read_yaml(tmp_path / "materialized" / "1VFB" / "system.yml")
    core_v1_mutations = read_csv_rows(tmp_path / "materialized" / "1VFB" / "core_v1_mutations.csv")
    core_v2_mutations = read_csv_rows(tmp_path / "materialized" / "1VFB" / "core_v2_mutations.csv")

    assert system_yml["system_name"] == "1vfb"
    assert system_yml["antibody_chains"] == ["H", "L"]
    assert system_yml["antigen_chains"] == ["C"]
    assert len(core_v1_mutations) == 1
    assert len(core_v2_mutations) == 3


def test_plan_ab_bind_batches_from_materialized_manifest(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    runs_root = tmp_path / "runs"
    summary = plan_ab_bind_batches(
        benchmark_root,
        protocol_path,
        spec_name="core_v1",
        runs_root=runs_root,
        complex_ids=["1VFB"],
    )

    assert summary["planned_batch_count"] == 1
    assert summary["planned_complexes"] == ["1VFB"]
    batch = summary["batches"][0]
    batch_dir = Path(batch["batch_dir"])
    assert batch["batch_id"] == "abbind_1vfb_core_v1"
    assert batch["job_count"] == 1
    assert batch_dir.is_dir()
    assert (batch_dir / "jobs.csv").exists()
    assert (runs_root / "plan_index.json").exists()


def test_plan_ab_bind_batches_accepts_split_selection(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    materialized_root = benchmark_root / "materialized"
    manifests_dir.mkdir(parents=True)

    manifest_rows = []
    for complex_id in ("1VFB", "1JRH"):
        materialized_dir = materialized_root / complex_id
        materialized_dir.mkdir(parents=True)
        structure_path = materialized_dir / f"{complex_id}.pdb"
        structure_path.write_text("HEADER\nEND\n", encoding="utf-8")
        system_yml = materialized_dir / "system.yml"
        write_yaml(
            system_yml,
            {
                "system_name": complex_id.lower(),
                "input_structure": str(structure_path),
                "structure_source": "experimental",
                "antibody_chains": ["H", "L"],
                "antigen_chains": ["C"],
                "notes": [],
            },
        )
        mutations_csv = materialized_dir / "core_v1_mutations.csv"
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": f"{complex_id.lower()}_0001",
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "A",
                    "entity_side": "antibody",
                }
            ],
            ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
        )
        manifest_rows.append(
            {
                "complex_id": complex_id,
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        )

    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        manifest_rows,
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    split_path = benchmark_root / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "validation": {"complex_ids": ["1JRH"]},
            },
        },
    )
    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(protocol_path, {"preset": "single_point", "lambda_windows": 2, "repeats": 1, "temperature_k": 310.0})

    runs_root = tmp_path / "runs"
    summary = plan_ab_bind_batches(
        benchmark_root,
        protocol_path,
        spec_name="core_v1",
        runs_root=runs_root,
        split_name="validation",
        split_path=split_path,
    )

    assert summary["planned_batch_count"] == 1
    assert summary["planned_complexes"] == ["1JRH"]
    assert summary["split_name"] == "validation"
    assert summary["split_path"] == str(split_path)


def test_run_and_report_ab_bind_plan_root(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(
        benchmark_root,
        protocol_path,
        spec_name="core_v1",
        runs_root=plan_root,
        complex_ids=["1VFB"],
    )
    run_payload = run_ab_bind_plan(plan_root, execute=False, complex_ids=["1VFB"], to_stage="prepare")
    report_payload = report_ab_bind_plan(plan_root, complex_ids=["1VFB"])

    assert run_payload["selected_batch_count"] == 1
    assert run_payload["selected_job_count"] == 1
    assert run_payload["canonical_reports_dir"] == str(plan_root / "reports")
    assert run_payload["execution_rows"][0]["final_stage"] == "prepare"
    assert run_payload["execution_rows"][0]["final_state"] == "completed"
    assert report_payload["selected_batch_count"] == 1
    assert report_payload["selected_job_count"] == 1
    assert report_payload["latest_stage_name_counts"]["prepare"] == 1
    assert report_payload["latest_stage_state_counts"]["completed"] == 1
    assert report_payload["qc_counts"]["not_evaluated"] == 1
    assert Path(report_payload["reports_dir"]).name == "complex-1vfb"
    assert report_payload["canonical_reports_dir"] == str(plan_root / "reports")
    assert (plan_root / "reports" / "run_summary.json").exists()
    assert (plan_root / "reports" / "plan_summary.json").exists()


def test_run_ab_bind_plan_parallel_execution_assigns_gpu_devices(tmp_path: Path, monkeypatch) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 33,
                "icode": "",
                "wt": "Y",
                "mut": "F",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0003",
                "chain_id": "L",
                "resseq": 50,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 3,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(
        benchmark_root,
        protocol_path,
        spec_name="core_v1",
        runs_root=plan_root,
        complex_ids=["1VFB"],
    )

    captured: list[tuple[str, str | None]] = []

    def fake_run_job(
        job_dir: Path,
        *,
        execute: bool,
        from_stage: str | None = None,
        to_stage: str | None = None,
        environment: dict[str, str] | None = None,
    ) -> list[StageStatus]:
        captured.append((job_dir.name, None if environment is None else environment.get("CUDA_VISIBLE_DEVICES")))
        return [
            StageStatus(
                stage=to_stage or "prepare",
                state="completed",
                message=f"Executed on {environment.get('CUDA_VISIBLE_DEVICES', '') if environment else ''}",
            )
        ]

    monkeypatch.setattr("abag_rbfe.benchmark.run_job", fake_run_job)
    monkeypatch.setattr("abag_rbfe.benchmark.write_batch_summary", lambda _batch_dir: None)
    monkeypatch.setattr(
        "abag_rbfe.benchmark.report_ab_bind_plan",
        lambda *args, **kwargs: {"reports_dir": str(plan_root / "reports")},
    )

    run_payload = run_ab_bind_plan(
        plan_root,
        execute=True,
        complex_ids=["1VFB"],
        to_stage="prepare",
        max_workers=2,
        gpu_devices=["2", "5"],
    )

    assert run_payload["selected_job_count"] == 3
    assert {device for _, device in captured} == {"2", "5"}
    assert {row["gpu_device"] for row in run_payload["execution_rows"]} == {"2", "5"}
    assert all(row["final_stage"] == "prepare" for row in run_payload["execution_rows"])
    assert all(row["final_state"] == "completed" for row in run_payload["execution_rows"])


def test_filtered_report_keeps_canonical_root_summary(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    manifests_dir.mkdir(parents=True)

    def write_complex(complex_id: str, mutation_group_id: str, chain_id: str, resseq: int, wt: str, mut: str) -> tuple[Path, Path]:
        materialized_dir = benchmark_root / "materialized" / complex_id
        materialized_dir.mkdir(parents=True)
        structure_path = materialized_dir / f"{complex_id}.pdb"
        structure_path.write_text(
            "\n".join(
                [
                    "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                    "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                    "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                    "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                    "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                    "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                    "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                    "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                    "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                    "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                    "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                    "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                    "TER",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        system_yml = materialized_dir / "system.yml"
        write_yaml(
            system_yml,
            {
                "system_name": complex_id.lower(),
                "input_structure": str(structure_path),
                "structure_source": "experimental",
                "antibody_chains": ["H", "L"],
                "antigen_chains": ["C"],
                "notes": [],
            },
        )
        mutations_csv = materialized_dir / "core_v1_mutations.csv"
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": mutation_group_id,
                    "chain_id": chain_id,
                    "resseq": resseq,
                    "icode": "",
                    "wt": wt,
                    "mut": mut,
                    "entity_side": "antibody",
                }
            ],
            ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
        )
        return system_yml, mutations_csv

    vfb_system, vfb_mutations = write_complex("1VFB", "1vfb_0001", "H", 32, "Y", "A")
    jrh_system, jrh_mutations = write_complex("1JRH", "1jrh_0001", "H", 33, "Y", "A")
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(vfb_system),
                "mutations_csv": str(vfb_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
            {
                "complex_id": "1JRH",
                "system_yml": str(jrh_system),
                "mutations_csv": str(jrh_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root)

    run_payload = run_ab_bind_plan(plan_root, execute=False, complex_ids=["1VFB"], to_stage="prepare")
    canonical_summary = read_json(plan_root / "reports" / "plan_summary.json")

    assert run_payload["canonical_reports_dir"] == str(plan_root / "reports")
    assert canonical_summary["selected_batch_count"] == 2
    assert canonical_summary["selected_job_count"] == 2
    assert canonical_summary["latest_stage_name_counts"]["prepare"] == 1
    assert canonical_summary["latest_stage_name_counts"]["not_started"] == 1
    assert canonical_summary["qc_counts"]["not_evaluated"] == 1
    assert canonical_summary["qc_counts"]["not_started"] == 1

    filtered_payload = report_ab_bind_plan(plan_root, complex_ids=["1VFB"])
    filtered_summary = read_json(Path(filtered_payload["reports_dir"]) / "plan_summary.json")
    canonical_summary = read_json(plan_root / "reports" / "plan_summary.json")

    assert filtered_payload["selected_batch_count"] == 1
    assert filtered_payload["selected_job_count"] == 1
    assert filtered_payload["canonical_reports_dir"] == str(plan_root / "reports")
    assert filtered_summary["selected_batch_count"] == 1
    assert filtered_summary["latest_stage_name_counts"]["prepare"] == 1
    assert canonical_summary["selected_batch_count"] == 2
    assert canonical_summary["latest_stage_name_counts"]["prepare"] == 1
    assert canonical_summary["latest_stage_name_counts"]["not_started"] == 1


def test_report_ab_bind_plan_accepts_split_selection(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    manifests_dir.mkdir(parents=True)

    def write_complex(complex_id: str) -> tuple[Path, Path]:
        materialized_dir = benchmark_root / "materialized" / complex_id
        materialized_dir.mkdir(parents=True)
        structure_path = materialized_dir / f"{complex_id}.pdb"
        structure_path.write_text("HEADER\nEND\n", encoding="utf-8")
        system_yml = materialized_dir / "system.yml"
        write_yaml(
            system_yml,
            {
                "system_name": complex_id.lower(),
                "input_structure": str(structure_path),
                "structure_source": "experimental",
                "antibody_chains": ["H", "L"],
                "antigen_chains": ["C"],
                "notes": [],
            },
        )
        mutations_csv = materialized_dir / "core_v1_mutations.csv"
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": f"{complex_id.lower()}_0001",
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "A",
                    "entity_side": "antibody",
                }
            ],
            ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
        )
        return system_yml, mutations_csv

    vfb_system, vfb_mutations = write_complex("1VFB")
    jrh_system, jrh_mutations = write_complex("1JRH")
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(vfb_system),
                "mutations_csv": str(vfb_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
            {
                "complex_id": "1JRH",
                "system_yml": str(jrh_system),
                "mutations_csv": str(jrh_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    split_path = benchmark_root / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "validation": {"complex_ids": ["1JRH"]},
            },
        },
    )
    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(protocol_path, {"preset": "single_point", "lambda_windows": 2, "repeats": 1, "temperature_k": 310.0})

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root)
    run_payload = run_ab_bind_plan(
        plan_root,
        execute=False,
        split_name="validation",
        split_path=split_path,
        to_stage="prepare",
    )
    report_payload = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    canonical_summary = read_json(plan_root / "reports" / "plan_summary.json")
    validation_summary = read_json(Path(report_payload["reports_dir"]) / "plan_summary.json")

    assert run_payload["canonical_reports_dir"] == str(plan_root / "reports")
    assert report_payload["selection"]["split_name"] == "validation"
    assert report_payload["selection"]["complex_ids"] == ["1JRH"]
    assert Path(report_payload["reports_dir"]).name == "split-validation-complex-1jrh"
    assert validation_summary["selected_batch_count"] == 1
    assert validation_summary["latest_stage_name_counts"]["prepare"] == 1
    assert canonical_summary["selected_batch_count"] == 2
    assert canonical_summary["latest_stage_name_counts"]["prepare"] == 1
    assert canonical_summary["latest_stage_name_counts"]["not_started"] == 1


def test_report_ab_bind_plan_emits_benchmark_metrics_for_ready_jobs(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 2,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.6160333201786573",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "-1.2320666403573146",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(
        benchmark_root,
        protocol_path,
        spec_name="core_v1",
        runs_root=plan_root,
        complex_ids=["1VFB"],
    )
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    expected_by_group = {
        "1vfb_0001": (1.0, 0.0),
        "1vfb_0002": (0.0, 2.0),
    }
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        spec = read_json(job_dir / "job_spec.json")
        mutation_group_id = spec["mutation_group"]["mutation_group_id"]
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        complex_delta_kt, apo_delta_kt = expected_by_group[mutation_group_id]
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=complex_delta_kt, stderr_kt=0.10)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=apo_delta_kt, stderr_kt=0.10)

    report_payload = report_ab_bind_plan(plan_root, complex_ids=["1VFB"])
    metrics = report_payload["benchmark_metrics"]
    qualified_metrics = report_payload["benchmark_metrics_qc_qualified"]

    assert report_payload["paired_job_count"] == 2
    assert report_payload["qc_qualified_pair_count"] == 2
    assert report_payload["validation_gate"]["overall_pearson_r_threshold"] == 0.6
    assert report_payload["validation_gate"]["overall_pearson_r_passed"] is True
    assert isclose(metrics["pearson_r"], 1.0, rel_tol=1e-9)
    assert isclose(metrics["spearman_rho"], 1.0, rel_tol=1e-9)
    assert isclose(metrics["rmse_kcal_mol"], 0.0, abs_tol=1e-12)
    assert isclose(metrics["mae_kcal_mol"], 0.0, abs_tol=1e-12)
    assert isclose(metrics["sign_accuracy"], 1.0, rel_tol=1e-9)
    assert isclose(metrics["auc_strong_effect"], 1.0, rel_tol=1e-9)
    assert qualified_metrics == metrics
    assert (plan_root / "reports" / "benchmark_metrics.json").exists()
    assert (plan_root / "reports" / "benchmark_metrics_qc_qualified.json").exists()
    assert len(read_csv_rows(plan_root / "reports" / "benchmark_pairs.csv")) == 2
    assert len(read_csv_rows(plan_root / "reports" / "benchmark_pairs_qc_qualified.csv")) == 2
    plan_jobs = read_csv_rows(plan_root / "reports" / "plan_jobs.csv")
    assert {row["mutation_group_id"] for row in plan_jobs} == {"1vfb_0001", "1vfb_0002"}
    assert {row["latest_stage"] for row in plan_jobs} == {"qc"}
    assert all(row["experimental_ddg_kcal_mol"] != "" for row in plan_jobs)
    assert all(row["benchmark_qc_qualified"] == "True" for row in plan_jobs)


def test_report_ab_bind_plan_emits_resumable_job_counts(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    run_ab_bind_plan(plan_root, execute=False, complex_ids=["1VFB"], to_stage="build_legs")
    report_payload = report_ab_bind_plan(plan_root, complex_ids=["1VFB"])

    assert report_payload["resumable_job_count"] == 1
    assert report_payload["analyzable_job_count"] == 0
    assert report_payload["validation_gate"]["overall_pearson_r_passed"] is False
    assert report_payload["batches"][0]["resumable_job_count"] == 1
    assert report_payload["batches"][0]["analyzable_job_count"] == 0
    plan_jobs = read_csv_rows(plan_root / "reports" / "plan_jobs.csv")
    assert plan_jobs[0]["resumable"] == "True"
    assert plan_jobs[0]["analyzable"] == "False"


def test_report_ab_bind_plan_merges_best_rows_from_extra_plan_roots(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            },
            {
                "mutation_group_id": "cmpx_0002",
                "ddg_kcal_mol": "1.23",
                "source_mutation": "H:Y50A",
                "mutation_tokens": "H:Y50A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(job_id: str, mutation_group_id: str) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        complex_delta_kt: float,
        apo_delta_kt: float,
        stderr_kt: float,
    ) -> None:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=complex_delta_kt, stderr_kt=stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=apo_delta_kt, stderr_kt=stderr_kt)

    primary_root = tmp_path / "runs-primary"
    extra_root = tmp_path / "runs-robust"
    primary_batch_id = "abbind_cmpx_core_v1"
    extra_batch_id = "abbind-robust_cmpx_core_v1"

    write_ready_job(
        primary_root,
        batch_id=primary_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        complex_delta_kt=5.0,
        apo_delta_kt=8.0,
        stderr_kt=3.0,
    )
    write_ready_job(
        primary_root,
        batch_id=primary_batch_id,
        job_id="cmpx-job-b",
        mutation_group_id="cmpx_0002",
        complex_delta_kt=6.0,
        apo_delta_kt=4.0,
        stderr_kt=0.1,
    )
    write_ready_job(
        extra_root,
        batch_id=extra_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        complex_delta_kt=10.0,
        apo_delta_kt=5.0,
        stderr_kt=0.1,
    )

    primary_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(primary_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": primary_batch_id,
                "batch_dir": str(primary_root / primary_batch_id),
            }
        ],
    }
    extra_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(extra_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": extra_batch_id,
                "batch_dir": str(extra_root / extra_batch_id),
            }
        ],
    }
    write_json(primary_root / "plan_index.json", primary_index)
    write_yaml(primary_root / "plan_index.yml", primary_index)
    write_json(extra_root / "plan_index.json", extra_index)
    write_yaml(extra_root / "plan_index.yml", extra_index)

    merged_payload = report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root])

    assert merged_payload["reports_dir"] == str(primary_root / "reports" / "merged")
    assert merged_payload["selected_job_count"] == 2
    assert merged_payload["ddg_ready_count"] == 2
    assert merged_payload["source_plan_roots"] == [str(primary_root), str(extra_root)]

    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    assert merged_jobs["cmpx_0001"]["source_plan_root"] == str(extra_root)
    assert merged_jobs["cmpx_0001"]["batch_id"] == extra_batch_id
    assert merged_jobs["cmpx_0002"]["source_plan_root"] == str(primary_root)
    assert merged_jobs["cmpx_0002"]["batch_id"] == primary_batch_id

    merged_pairs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "benchmark_pairs.csv")
    }
    assert merged_pairs["cmpx_0001"]["source_plan_root"] == str(extra_root)
    assert merged_pairs["cmpx_0002"]["source_plan_root"] == str(primary_root)


def test_report_ab_bind_plan_marks_ready_winners_with_active_alternates(tmp_path: Path, monkeypatch) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(job_id: str, mutation_group_id: str, *, repeats: int, lambda_windows: int, production_ps: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": lambda_windows,
                "repeats": repeats,
                "production_ps": production_ps,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def write_ready_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> None:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, repeats=1, lambda_windows=2, production_ps=10)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=5.0, stderr_kt=0.1)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=8.0, stderr_kt=0.1)

    def write_running_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> Path:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, repeats=4, lambda_windows=6, production_ps=20)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
            _write_stage(job_dir, stage)
        _write_running_stage(
            job_dir,
            "sample",
            commands=["bash sample.sh"],
            artifacts=[str(job_dir / "artifacts" / "commands" / "sample.sh")],
        )
        _write_sample_window_completed(job_dir / "legs" / "complex" / "rep01" / "lambda_000")
        _write_sample_window_started(job_dir / "legs" / "complex" / "rep01" / "lambda_001")
        return job_dir

    primary_root = tmp_path / "runs-primary"
    extra_root = tmp_path / "runs-robust"
    primary_batch_id = "abbind_cmpx_core_v1"
    extra_batch_id = "abbind-robust_cmpx_core_v1"

    write_ready_job(primary_root, batch_id=primary_batch_id, job_id="cmpx-job-a", mutation_group_id="cmpx_0001")
    running_job_dir = write_running_job(
        extra_root, batch_id=extra_batch_id, job_id="cmpx-job-a", mutation_group_id="cmpx_0001"
    )

    primary_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(primary_root),
        "batches": [{"complex_id": "CMPX", "batch_id": primary_batch_id, "batch_dir": str(primary_root / primary_batch_id)}],
    }
    extra_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(extra_root),
        "batches": [{"complex_id": "CMPX", "batch_id": extra_batch_id, "batch_dir": str(extra_root / extra_batch_id)}],
    }
    write_json(primary_root / "plan_index.json", primary_index)
    write_yaml(primary_root / "plan_index.yml", primary_index)
    write_json(extra_root / "plan_index.json", extra_index)
    write_yaml(extra_root / "plan_index.yml", extra_index)

    monkeypatch.setattr(
        "abag_rbfe.reporting._active_process_lines",
        lambda: (f"fake {running_job_dir}/artifacts/commands/sample.sh",),
    )

    merged_payload = report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root])

    assert merged_payload["active_alternate_job_count"] == 1
    assert merged_payload["active_alternate_ready_job_count"] == 1
    assert merged_payload["active_alternate_ready_hotspot_count"] == 1
    assert merged_payload["active_alternate_ready_hotspots"][0]["job_id"] == "cmpx-job-a"
    assert merged_payload["active_alternate_ready_hotspots"][0]["mutation_group_id"] == "cmpx_0001"
    assert merged_payload["active_alternate_ready_hotspots"][0]["active_alternate_stage_states"] == "sample:running"
    assert merged_payload["active_alternate_ready_hotspots"][0]["active_alternate_current_source_plan_root"] == str(
        extra_root
    )
    assert merged_payload["active_alternate_ready_hotspots"][0]["active_alternate_current_latest_stage"] == "sample"
    assert merged_payload["active_alternate_ready_hotspots"][0]["active_alternate_current_sample_active_phase"] == "md"
    assert (
        merged_payload["active_alternate_ready_hotspots"][0]["active_alternate_current_sample_active_window"]
        == "complex/rep01/lambda_001"
    )

    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    winner = merged_jobs["cmpx_0001"]
    assert winner["source_plan_root"] == str(primary_root)
    assert winner["ddg_ready"] == "True"
    assert winner["has_active_alternate_candidate"] == "True"
    assert winner["active_alternate_candidate_count"] == "1"
    assert winner["active_alternate_source_plan_roots"] == str(extra_root)
    assert winner["active_alternate_stage_states"] == "sample:running"
    assert winner["active_alternate_current_source_plan_root"] == str(extra_root)
    assert winner["active_alternate_current_latest_stage"] == "sample"
    assert winner["active_alternate_current_latest_stage_state"] == "running"
    assert winner["active_alternate_current_sample_started_windows"] == "2"
    assert winner["active_alternate_current_sample_completed_windows"] == "1"
    assert winner["active_alternate_current_sample_active_phase"] == "md"
    assert winner["active_alternate_current_sample_active_window"] == "complex/rep01/lambda_001"
    assert winner["abs_ddg_error_kcal_mol"] != ""

    active_alternate_jobs = read_csv_rows(primary_root / "reports" / "merged" / "active_alternate_jobs.csv")
    assert len(active_alternate_jobs) == 1
    assert active_alternate_jobs[0]["job_id"] == "cmpx-job-a"
    assert active_alternate_jobs[0]["mutation_group_id"] == "cmpx_0001"
    assert active_alternate_jobs[0]["ddg_ready"] == "True"
    assert active_alternate_jobs[0]["active_alternate_stage_states"] == "sample:running"
    assert active_alternate_jobs[0]["active_alternate_current_source_plan_root"] == str(extra_root)
    assert active_alternate_jobs[0]["active_alternate_current_sample_active_phase"] == "md"
    assert active_alternate_jobs[0]["active_alternate_current_sample_active_window"] == "complex/rep01/lambda_001"
    assert active_alternate_jobs[0]["abs_ddg_error_kcal_mol"] == winner["abs_ddg_error_kcal_mol"]


def test_report_ab_bind_plan_counts_not_started_deep_rescues_as_alternates(tmp_path: Path, monkeypatch) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(job_id: str, mutation_group_id: str, *, repeats: int, lambda_windows: int, production_ps: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": lambda_windows,
                "repeats": repeats,
                "production_ps": production_ps,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def prepare_job_dir(
        plan_root: Path,
        *,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        repeats: int,
        lambda_windows: int,
        production_ps: int,
    ) -> Path:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(
            job_id,
            mutation_group_id,
            repeats=repeats,
            lambda_windows=lambda_windows,
            production_ps=production_ps,
        )
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        return job_dir

    def write_ready_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> Path:
        job_dir = prepare_job_dir(
            plan_root,
            batch_id=batch_id,
            job_id=job_id,
            mutation_group_id=mutation_group_id,
            repeats=1,
            lambda_windows=2,
            production_ps=10,
        )
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=5.0, stderr_kt=0.1)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=8.0, stderr_kt=0.1)
        return job_dir

    def write_running_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> Path:
        job_dir = prepare_job_dir(
            plan_root,
            batch_id=batch_id,
            job_id=job_id,
            mutation_group_id=mutation_group_id,
            repeats=4,
            lambda_windows=6,
            production_ps=20,
        )
        for stage in ("ingest", "prepare", "mutate", "build_legs"):
            _write_stage(job_dir, stage)
        _write_running_stage(
            job_dir,
            "equilibrate",
            commands=["bash equilibrate.sh"],
            artifacts=[str(job_dir / "artifacts" / "commands" / "equilibrate.sh")],
        )
        _write_equilibrate_repeat_completed(job_dir / "legs" / "complex" / "rep01")
        _write_equilibrate_repeat_started(job_dir / "legs" / "apo" / "rep01")
        return job_dir

    def write_running_sample_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> Path:
        job_dir = prepare_job_dir(
            plan_root,
            batch_id=batch_id,
            job_id=job_id,
            mutation_group_id=mutation_group_id,
            repeats=4,
            lambda_windows=6,
            production_ps=20,
        )
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
            _write_stage(job_dir, stage)
        _write_running_stage(
            job_dir,
            "sample",
            commands=["bash sample.sh"],
            artifacts=[str(job_dir / "artifacts" / "commands" / "sample.sh")],
        )
        _write_sample_window_completed(job_dir / "legs" / "complex" / "rep01" / "lambda_000")
        _write_sample_window_started(job_dir / "legs" / "complex" / "rep01" / "lambda_001")
        return job_dir

    def write_not_started_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> Path:
        return prepare_job_dir(
            plan_root,
            batch_id=batch_id,
            job_id=job_id,
            mutation_group_id=mutation_group_id,
            repeats=5,
            lambda_windows=8,
            production_ps=40,
        )

    primary_root = tmp_path / "runs-primary"
    robust_root = tmp_path / "runs-robust"
    rescue_root = tmp_path / "runs-rescue"
    deep_root = tmp_path / "runs-deep"

    primary_batch_id = "abbind_cmpx_core_v1"
    robust_batch_id = "abbind-robust_cmpx_core_v1"
    rescue_batch_id = "abbind-rescue_cmpx_core_v1"
    deep_batch_id = "abbind-deep-rescue_cmpx_core_v1"

    write_ready_job(primary_root, batch_id=primary_batch_id, job_id="cmpx-job-a", mutation_group_id="cmpx_0001")
    robust_job_dir = write_running_job(
        robust_root, batch_id=robust_batch_id, job_id="cmpx-job-a", mutation_group_id="cmpx_0001"
    )
    rescue_job_dir = write_running_sample_job(
        rescue_root, batch_id=rescue_batch_id, job_id="cmpx-job-a", mutation_group_id="cmpx_0001"
    )
    write_not_started_job(deep_root, batch_id=deep_batch_id, job_id="cmpx-job-a", mutation_group_id="cmpx_0001")

    for root, batch_id in (
        (primary_root, primary_batch_id),
        (robust_root, robust_batch_id),
        (rescue_root, rescue_batch_id),
        (deep_root, deep_batch_id),
    ):
        index = {
            "benchmark_root": str(benchmark_root),
            "spec_name": "core_v1",
            "plan_root": str(root),
            "batches": [{"complex_id": "CMPX", "batch_id": batch_id, "batch_dir": str(root / batch_id)}],
        }
        write_json(root / "plan_index.json", index)
        write_yaml(root / "plan_index.yml", index)

    monkeypatch.setattr(
        "abag_rbfe.reporting._active_process_lines",
        lambda: (
            f"fake {robust_job_dir}/artifacts/commands/equilibrate.sh",
            f"fake {rescue_job_dir}/artifacts/commands/sample.sh",
        ),
    )

    merged_payload = report_ab_bind_plan(
        primary_root,
        extra_plan_roots=[robust_root, rescue_root, deep_root],
    )

    assert merged_payload["source_plan_roots"] == [
        str(primary_root),
        str(robust_root),
        str(rescue_root),
        str(deep_root),
    ]
    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    winner = merged_jobs["cmpx_0001"]
    assert winner["source_plan_root"] == str(primary_root)
    assert winner["alternate_candidate_count"] == "3"
    assert winner["active_alternate_candidate_count"] == "2"
    assert set(winner["active_alternate_source_plan_roots"].split(",")) == {str(robust_root), str(rescue_root)}
    assert winner["active_alternate_stage_states"] == "equilibrate:running,sample:running"
    assert winner["active_alternate_current_source_plan_root"] == str(rescue_root)
    assert winner["active_alternate_current_latest_stage"] == "sample"
    assert winner["active_alternate_current_latest_stage_state"] == "running"
    assert winner["active_alternate_current_sample_active_phase"] == "md"
    assert winner["active_alternate_current_sample_active_window"] == "complex/rep01/lambda_001"

    merged_summary = read_json(primary_root / "reports" / "merged" / "plan_summary.json")
    assert merged_summary["source_plan_roots"] == [
        str(primary_root),
        str(robust_root),
        str(rescue_root),
        str(deep_root),
    ]


def test_report_ab_bind_plan_prefers_sampling_effort_before_bar_stderr(tmp_path: Path) -> None:
    """Deeper-protocol rows must win over shallower rows at the same QC level,
    even when the shallow row reports a lower BAR stderr (regression: quick-preset
    values shadowed deeper rescue reruns in merged validation reports)."""
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(
        job_id: str,
        mutation_group_id: str,
        *,
        repeats: int,
        lambda_windows: int,
        production_ps: int,
    ) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": lambda_windows,
                "repeats": repeats,
                "production_ps": production_ps,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        complex_delta_kt: float,
        apo_delta_kt: float,
        stderr_kt: float,
        repeats: int,
        lambda_windows: int,
        production_ps: int,
    ) -> None:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(
            job_id,
            mutation_group_id,
            repeats=repeats,
            lambda_windows=lambda_windows,
            production_ps=production_ps,
        )
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=complex_delta_kt, stderr_kt=stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=apo_delta_kt, stderr_kt=stderr_kt)

    primary_root = tmp_path / "runs-primary"
    extra_root = tmp_path / "runs-extra"
    primary_batch_id = "abbind_cmpx_core_v1"
    extra_batch_id = "abbind-rescue_cmpx_core_v1"

    write_ready_job(
        primary_root,
        batch_id=primary_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        complex_delta_kt=8.0,
        apo_delta_kt=5.0,
        stderr_kt=0.1,
        repeats=1,
        lambda_windows=2,
        production_ps=2,
    )
    write_ready_job(
        extra_root,
        batch_id=extra_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        complex_delta_kt=9.0,
        apo_delta_kt=5.5,
        stderr_kt=4.0,
        repeats=1,
        lambda_windows=2,
        production_ps=20,
    )

    primary_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(primary_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": primary_batch_id,
                "batch_dir": str(primary_root / primary_batch_id),
            }
        ],
    }
    extra_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(extra_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": extra_batch_id,
                "batch_dir": str(extra_root / extra_batch_id),
            }
        ],
    }
    write_json(primary_root / "plan_index.json", primary_index)
    write_yaml(primary_root / "plan_index.yml", primary_index)
    write_json(extra_root / "plan_index.json", extra_index)
    write_yaml(extra_root / "plan_index.yml", extra_index)

    report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root])

    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    assert merged_jobs["cmpx_0001"]["source_plan_root"] == str(extra_root)
    assert merged_jobs["cmpx_0001"]["protocol_production_ps"] == "20"


def test_report_ab_bind_plan_uses_sampling_effort_as_tiebreak_after_quality_metrics(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(
        job_id: str,
        mutation_group_id: str,
        *,
        repeats: int,
        lambda_windows: int,
        production_ps: int,
    ) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": lambda_windows,
                "repeats": repeats,
                "production_ps": production_ps,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        complex_delta_kt: float,
        apo_delta_kt: float,
        stderr_kt: float,
        repeats: int,
        lambda_windows: int,
        production_ps: int,
    ) -> None:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(
            job_id,
            mutation_group_id,
            repeats=repeats,
            lambda_windows=lambda_windows,
            production_ps=production_ps,
        )
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=complex_delta_kt, stderr_kt=stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=apo_delta_kt, stderr_kt=stderr_kt)

    primary_root = tmp_path / "runs-primary"
    extra_root = tmp_path / "runs-extra"
    primary_batch_id = "abbind_cmpx_core_v1"
    extra_batch_id = "abbind-rescue_cmpx_core_v1"

    write_ready_job(
        primary_root,
        batch_id=primary_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        complex_delta_kt=8.0,
        apo_delta_kt=5.0,
        stderr_kt=0.1,
        repeats=1,
        lambda_windows=2,
        production_ps=2,
    )
    write_ready_job(
        extra_root,
        batch_id=extra_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        complex_delta_kt=9.0,
        apo_delta_kt=5.5,
        stderr_kt=0.1,
        repeats=1,
        lambda_windows=2,
        production_ps=20,
    )

    primary_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(primary_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": primary_batch_id,
                "batch_dir": str(primary_root / primary_batch_id),
            }
        ],
    }
    extra_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(extra_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": extra_batch_id,
                "batch_dir": str(extra_root / extra_batch_id),
            }
        ],
    }
    write_json(primary_root / "plan_index.json", primary_index)
    write_yaml(primary_root / "plan_index.yml", primary_index)
    write_json(extra_root / "plan_index.json", extra_index)
    write_yaml(extra_root / "plan_index.yml", extra_index)

    report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root])

    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    assert merged_jobs["cmpx_0001"]["source_plan_root"] == str(extra_root)
    assert merged_jobs["cmpx_0001"]["protocol_production_ps"] == "20"


def test_report_ab_bind_plan_preserves_more_progressed_live_rows_before_protocol_strength(
    tmp_path: Path,
    monkeypatch,
) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(job_id: str, mutation_group_id: str, *, production_ps: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": production_ps,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def write_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str, production_ps: int) -> Path:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, production_ps=production_ps)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        return job_dir

    primary_root = tmp_path / "runs-primary"
    extra_root = tmp_path / "runs-extra"
    primary_batch_id = "abbind_cmpx_core_v1"
    extra_batch_id = "abbind-rescue_cmpx_core_v1"

    running_job_dir = write_job(
        primary_root,
        batch_id=primary_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        production_ps=2,
    )
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(running_job_dir, stage)
    _write_running_stage(
        running_job_dir,
        "sample",
        commands=["bash sample.sh"],
        artifacts=[str(running_job_dir / "artifacts" / "commands" / "sample.sh")],
    )
    _write_equilibrate_repeat_completed(running_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(running_job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(running_job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(running_job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    stale_job_dir = write_job(
        extra_root,
        batch_id=extra_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
        production_ps=20,
    )
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(stale_job_dir, stage)
    _write_running_stage(
        stale_job_dir,
        "equilibrate",
        commands=["bash equilibrate.sh"],
        artifacts=[str(stale_job_dir / "artifacts" / "commands" / "equilibrate.sh")],
    )
    _write_equilibrate_repeat_completed(stale_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(stale_job_dir / "legs" / "apo" / "rep01")

    primary_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(primary_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": primary_batch_id,
                "batch_dir": str(primary_root / primary_batch_id),
            }
        ],
    }
    extra_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(extra_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": extra_batch_id,
                "batch_dir": str(extra_root / extra_batch_id),
            }
        ],
    }
    write_json(primary_root / "plan_index.json", primary_index)
    write_yaml(primary_root / "plan_index.yml", primary_index)
    write_json(extra_root / "plan_index.json", extra_index)
    write_yaml(extra_root / "plan_index.yml", extra_index)

    monkeypatch.setattr(
        "abag_rbfe.reporting._active_process_lines",
        lambda: (
            f"fake {running_job_dir}/artifacts/commands/sample.sh",
        ),
    )

    report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root])

    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    assert merged_jobs["cmpx_0001"]["source_plan_root"] == str(primary_root)
    assert merged_jobs["cmpx_0001"]["latest_stage"] == "sample"
    assert merged_jobs["cmpx_0001"]["latest_stage_state"] == "running"


def test_report_ab_bind_plan_prefers_live_running_rows_over_stale_higher_stage_rows(
    tmp_path: Path,
    monkeypatch,
) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    curated_dir.mkdir(parents=True)
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cmpx_0001",
                "ddg_kcal_mol": "3.08",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    def job_spec(job_id: str, mutation_group_id: str) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antibody",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "H",
                        "resseq": 33,
                        "icode": "",
                        "wt": "Y",
                        "mut": "A",
                        "entity_side": "antibody",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 1.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": "cmpx",
                "structure_source": "experimental",
            },
        }

    def write_job(plan_root: Path, *, batch_id: str, job_id: str, mutation_group_id: str) -> Path:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "config")
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        write_yaml(job_dir / "config" / "system.yml", spec["system"])
        write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
        return job_dir

    primary_root = tmp_path / "runs-primary"
    extra_root = tmp_path / "runs-robust"
    primary_batch_id = "abbind_cmpx_core_v1"
    extra_batch_id = "abbind-robust_cmpx_core_v1"

    stale_job_dir = write_job(
        primary_root,
        batch_id=primary_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
    )
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(stale_job_dir, stage)
    _write_running_stage(
        stale_job_dir,
        "sample",
        commands=["bash sample.sh"],
        artifacts=[str(stale_job_dir / "artifacts" / "commands" / "sample.sh")],
    )
    _write_equilibrate_repeat_completed(stale_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(stale_job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(stale_job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(stale_job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    running_job_dir = write_job(
        extra_root,
        batch_id=extra_batch_id,
        job_id="cmpx-job-a",
        mutation_group_id="cmpx_0001",
    )
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(running_job_dir, stage)
    _write_running_stage(
        running_job_dir,
        "equilibrate",
        commands=["bash equilibrate.sh"],
        artifacts=[str(running_job_dir / "artifacts" / "commands" / "equilibrate.sh")],
    )
    _write_equilibrate_repeat_completed(running_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(running_job_dir / "legs" / "apo" / "rep01")

    primary_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(primary_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": primary_batch_id,
                "batch_dir": str(primary_root / primary_batch_id),
            }
        ],
    }
    extra_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(extra_root),
        "batches": [
            {
                "complex_id": "CMPX",
                "batch_id": extra_batch_id,
                "batch_dir": str(extra_root / extra_batch_id),
            }
        ],
    }
    write_json(primary_root / "plan_index.json", primary_index)
    write_yaml(primary_root / "plan_index.yml", primary_index)
    write_json(extra_root / "plan_index.json", extra_index)
    write_yaml(extra_root / "plan_index.yml", extra_index)

    monkeypatch.setattr(
        "abag_rbfe.reporting._active_process_lines",
        lambda: (
            f"fake {running_job_dir}/artifacts/commands/equilibrate.sh",
        ),
    )

    report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root])

    merged_jobs = {
        row["mutation_group_id"]: row for row in read_csv_rows(primary_root / "reports" / "merged" / "plan_jobs.csv")
    }
    assert merged_jobs["cmpx_0001"]["source_plan_root"] == str(extra_root)
    assert merged_jobs["cmpx_0001"]["latest_stage"] == "equilibrate"
    assert merged_jobs["cmpx_0001"]["latest_stage_state"] == "running"


def test_report_ab_bind_plan_aggregates_running_stage_progress(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    manifests_dir.mkdir(parents=True)

    def write_complex(complex_id: str, mutation_group_id: str) -> tuple[Path, Path]:
        materialized_dir = benchmark_root / "materialized" / complex_id
        materialized_dir.mkdir(parents=True)
        structure_path = materialized_dir / f"{complex_id}.pdb"
        structure_path.write_text(
            "\n".join(
                [
                    "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                    "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                    "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                    "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                    "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                    "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                    "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                    "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                    "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                    "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                    "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                    "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                    "TER",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        system_yml = materialized_dir / "system.yml"
        write_yaml(
            system_yml,
            {
                "system_name": complex_id.lower(),
                "input_structure": str(structure_path),
                "structure_source": "experimental",
                "antibody_chains": ["H", "L"],
                "antigen_chains": ["C"],
                "notes": [],
            },
        )
        mutations_csv = materialized_dir / "core_v1_mutations.csv"
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": mutation_group_id,
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "A",
                    "entity_side": "antibody",
                }
            ],
            ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
        )
        return system_yml, mutations_csv

    vfb_system, vfb_mutations = write_complex("1VFB", "1vfb_0001")
    jrh_system, jrh_mutations = write_complex("1JRH", "1jrh_0001")
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(vfb_system),
                "mutations_csv": str(vfb_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
            {
                "complex_id": "1JRH",
                "system_yml": str(jrh_system),
                "mutations_csv": str(jrh_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root)
    run_ab_bind_plan(plan_root, execute=False, to_stage="build_legs")

    sample_job_dir = next((plan_root / "abbind_1vfb_core_v1" / "jobs").iterdir())
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(sample_job_dir, stage)
    _write_running_stage(sample_job_dir, "sample")
    _write_equilibrate_repeat_completed(sample_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(sample_job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(sample_job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(sample_job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    equilibrate_job_dir = next((plan_root / "abbind_1jrh_core_v1" / "jobs").iterdir())
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(equilibrate_job_dir, stage)
    _write_running_stage(equilibrate_job_dir, "equilibrate")
    _write_equilibrate_repeat_completed(equilibrate_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(equilibrate_job_dir / "legs" / "apo" / "rep01")

    report_payload = report_ab_bind_plan(plan_root)

    assert report_payload["running_sample_job_count"] == 1
    assert report_payload["running_sample_started_windows"] == 2
    assert report_payload["running_sample_completed_windows"] == 1
    assert report_payload["running_sample_total_windows"] == 4
    assert report_payload["running_equilibrate_job_count"] == 1
    assert report_payload["running_equilibrate_started_repeats"] == 2
    assert report_payload["running_equilibrate_completed_repeats"] == 1
    assert report_payload["running_equilibrate_total_repeats"] == 2

    batch_rows = {row["batch_id"]: row for row in report_payload["batches"]}
    assert batch_rows["abbind_1vfb_core_v1"]["running_sample_completed_windows"] == 1
    assert batch_rows["abbind_1jrh_core_v1"]["running_equilibrate_completed_repeats"] == 1
    assert report_payload["resumable_job_count"] == 0
    assert batch_rows["abbind_1vfb_core_v1"]["resumable_job_count"] == 0
    assert batch_rows["abbind_1jrh_core_v1"]["resumable_job_count"] == 0

    plan_jobs = {row["job_id"]: row for row in read_csv_rows(plan_root / "reports" / "plan_jobs.csv")}
    assert plan_jobs[sample_job_dir.name]["sample_completed_windows"] == "1"
    assert plan_jobs[sample_job_dir.name]["sample_total_windows"] == "4"
    assert plan_jobs[sample_job_dir.name]["sample_active_phase"] == "md"
    assert plan_jobs[sample_job_dir.name]["sample_active_window"] == "complex/rep01/lambda_001"
    assert plan_jobs[sample_job_dir.name]["resumable"] == "False"
    assert plan_jobs[equilibrate_job_dir.name]["equilibrate_completed_repeats"] == "1"
    assert plan_jobs[equilibrate_job_dir.name]["equilibrate_total_repeats"] == "2"
    assert plan_jobs[equilibrate_job_dir.name]["resumable"] == "False"


def test_report_ab_bind_plan_treats_stale_running_jobs_as_resumable(tmp_path: Path, monkeypatch) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    manifests_dir.mkdir(parents=True)

    def write_complex(complex_id: str, mutation_group_id: str) -> tuple[Path, Path]:
        materialized_dir = benchmark_root / "materialized" / complex_id
        materialized_dir.mkdir(parents=True)
        structure_path = materialized_dir / f"{complex_id}.pdb"
        structure_path.write_text(
            "\n".join(
                [
                    "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                    "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                    "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                    "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                    "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                    "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                    "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                    "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                    "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                    "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                    "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                    "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                    "TER",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        system_yml = materialized_dir / "system.yml"
        write_yaml(
            system_yml,
            {
                "system_name": complex_id.lower(),
                "input_structure": str(structure_path),
                "structure_source": "experimental",
                "antibody_chains": ["H", "L"],
                "antigen_chains": ["C"],
                "notes": [],
            },
        )
        mutations_csv = materialized_dir / "core_v1_mutations.csv"
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": mutation_group_id,
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "A",
                    "entity_side": "antibody",
                }
            ],
            ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
        )
        return system_yml, mutations_csv

    vfb_system, vfb_mutations = write_complex("1VFB", "1vfb_0001")
    jrh_system, jrh_mutations = write_complex("1JRH", "1jrh_0001")
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(vfb_system),
                "mutations_csv": str(vfb_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
            {
                "complex_id": "1JRH",
                "system_yml": str(jrh_system),
                "mutations_csv": str(jrh_mutations),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            },
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root)
    run_ab_bind_plan(plan_root, execute=False, to_stage="build_legs")

    sample_job_dir = next((plan_root / "abbind_1vfb_core_v1" / "jobs").iterdir())
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(sample_job_dir, stage)
    _write_running_stage(
        sample_job_dir,
        "sample",
        commands=["bash sample.sh"],
        artifacts=[str(sample_job_dir / "artifacts" / "commands" / "sample.sh")],
    )
    _write_equilibrate_repeat_completed(sample_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(sample_job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(sample_job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(sample_job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    equilibrate_job_dir = next((plan_root / "abbind_1jrh_core_v1" / "jobs").iterdir())
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(equilibrate_job_dir, stage)
    _write_running_stage(
        equilibrate_job_dir,
        "equilibrate",
        commands=["bash equilibrate.sh"],
        artifacts=[str(equilibrate_job_dir / "artifacts" / "commands" / "equilibrate.sh")],
    )
    _write_equilibrate_repeat_completed(equilibrate_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(equilibrate_job_dir / "legs" / "apo" / "rep01")

    monkeypatch.setattr("abag_rbfe.reporting._active_process_lines", lambda: ())

    report_payload = report_ab_bind_plan(plan_root)

    assert report_payload["running_sample_job_count"] == 0
    assert report_payload["running_equilibrate_job_count"] == 0
    assert report_payload["resumable_job_count"] == 2

    batch_rows = {row["batch_id"]: row for row in report_payload["batches"]}
    assert batch_rows["abbind_1vfb_core_v1"]["running_sample_job_count"] == 0
    assert batch_rows["abbind_1vfb_core_v1"]["resumable_job_count"] == 1
    assert batch_rows["abbind_1jrh_core_v1"]["running_equilibrate_job_count"] == 0
    assert batch_rows["abbind_1jrh_core_v1"]["resumable_job_count"] == 1

    plan_jobs = {row["job_id"]: row for row in read_csv_rows(plan_root / "reports" / "plan_jobs.csv")}
    assert plan_jobs[sample_job_dir.name]["latest_stage_state"] == "stale_running"
    assert plan_jobs[sample_job_dir.name]["resumable"] == "True"
    assert plan_jobs[equilibrate_job_dir.name]["latest_stage_state"] == "stale_running"
    assert plan_jobs[equilibrate_job_dir.name]["resumable"] == "True"


def test_report_ab_bind_plan_writes_selection_reports_when_merging_extra_roots(
    tmp_path: Path,
    monkeypatch,
) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    manifests_dir.mkdir(parents=True)

    def write_complex(complex_id: str, mutation_group_id: str) -> tuple[Path, Path]:
        materialized_dir = benchmark_root / "materialized" / complex_id
        materialized_dir.mkdir(parents=True)
        structure_path = materialized_dir / f"{complex_id}.pdb"
        structure_path.write_text(
            "\n".join(
                [
                    "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                    "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                    "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                    "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                    "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                    "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                    "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                    "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                    "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                    "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                    "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                    "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                    "TER",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        system_yml = materialized_dir / "system.yml"
        write_yaml(
            system_yml,
            {
                "system_name": complex_id.lower(),
                "input_structure": str(structure_path),
                "structure_source": "experimental",
                "antibody_chains": ["H", "L"],
                "antigen_chains": ["C"],
                "notes": [],
            },
        )
        mutations_csv = materialized_dir / "core_v1_mutations.csv"
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": mutation_group_id,
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "A",
                    "entity_side": "antibody",
                }
            ],
            ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
        )
        return system_yml, mutations_csv

    system_yml, mutations_csv = write_complex("1VFB", "1vfb_0001")
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    primary_root = tmp_path / "primary"
    extra_root = tmp_path / "extra"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=primary_root)
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=extra_root)
    run_ab_bind_plan(primary_root, execute=False, to_stage="build_legs")
    run_ab_bind_plan(extra_root, execute=False, to_stage="build_legs")

    primary_job_dir = next((primary_root / "abbind_1vfb_core_v1" / "jobs").iterdir())
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(primary_job_dir, stage)
    _write_running_stage(
        primary_job_dir,
        "sample",
        commands=["bash sample.sh"],
        artifacts=[str(primary_job_dir / "artifacts" / "commands" / "sample.sh")],
    )
    _write_equilibrate_repeat_completed(primary_job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(primary_job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(primary_job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(primary_job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    monkeypatch.setattr("abag_rbfe.reporting._active_process_lines", lambda: ())

    merged_payload = report_ab_bind_plan(primary_root, extra_plan_roots=[extra_root], complex_ids=["1VFB"])

    selection_dir = primary_root / "reports" / "selections" / "complex-1vfb"
    selection_jobs = {row["job_id"]: row for row in read_csv_rows(selection_dir / "plan_jobs.csv")}
    canonical_merged_summary = read_json(primary_root / "reports" / "merged" / "plan_summary.json")

    assert selection_dir.is_dir()
    assert merged_payload["reports_dir"] == str(primary_root / "reports" / "merged" / "selections" / "complex-1vfb")
    assert canonical_merged_summary["reports_dir"] == str(primary_root / "reports" / "merged")
    assert "selection" not in canonical_merged_summary
    assert canonical_merged_summary["source_plan_roots"] == [
        str(primary_root),
        str(extra_root),
    ]
    assert selection_jobs[primary_job_dir.name]["latest_stage"] == "sample"
    assert selection_jobs[primary_job_dir.name]["latest_stage_state"] == "stale_running"
    assert selection_jobs[primary_job_dir.name]["resumable"] == "True"


def test_calibrate_ab_bind_plan_fits_side_specific_linear_model(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "cal1_0001",
                "ddg_kcal_mol": "1.0",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "cal2_0001",
                "ddg_kcal_mol": "2.0",
                "source_mutation": "H:Y33A",
                "mutation_tokens": "H:Y33A@antibody",
            },
            {
                "mutation_group_id": "cal3_0001",
                "ddg_kcal_mol": "2.0",
                "source_mutation": "W:G92A",
                "mutation_tokens": "W:G92A@antigen",
            },
            {
                "mutation_group_id": "cal4_0001",
                "ddg_kcal_mol": "4.0",
                "source_mutation": "W:G93A",
                "mutation_tokens": "W:G93A@antigen",
            },
            {
                "mutation_group_id": "val1_0001",
                "ddg_kcal_mol": "1.5",
                "source_mutation": "H:Y34A",
                "mutation_tokens": "H:Y34A@antibody",
            },
            {
                "mutation_group_id": "val2_0001",
                "ddg_kcal_mol": "3.0",
                "source_mutation": "W:G94A",
                "mutation_tokens": "W:G94A@antigen",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3", "CAL4"]},
                "validation": {"complex_ids": ["VAL1", "VAL2"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(
        job_id: str,
        mutation_group_id: str,
        *,
        entity_side: str,
        resseq: int,
    ) -> dict:
        chain_id = "H" if entity_side == "antibody" else "W"
        wt = "Y" if entity_side == "antibody" else "G"
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": entity_side,
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": chain_id,
                        "resseq": resseq,
                        "icode": "",
                        "wt": wt,
                        "mut": "A",
                        "entity_side": entity_side,
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        entity_side: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, entity_side=entity_side, resseq=resseq)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=0.05,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.05)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id="CAL1",
            batch_id="abbind_cal1_core_v1",
            job_id="cal1-antibody-h-y32a",
            mutation_group_id="cal1_0001",
            entity_side="antibody",
            resseq=32,
            raw_ddg_kcal_mol=10.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="CAL2",
            batch_id="abbind_cal2_core_v1",
            job_id="cal2-antibody-h-y33a",
            mutation_group_id="cal2_0001",
            entity_side="antibody",
            resseq=33,
            raw_ddg_kcal_mol=20.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="CAL3",
            batch_id="abbind_cal3_core_v1",
            job_id="cal3-antigen-w-g92a",
            mutation_group_id="cal3_0001",
            entity_side="antigen",
            resseq=92,
            raw_ddg_kcal_mol=5.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="CAL4",
            batch_id="abbind_cal4_core_v1",
            job_id="cal4-antigen-w-g93a",
            mutation_group_id="cal4_0001",
            entity_side="antigen",
            resseq=93,
            raw_ddg_kcal_mol=10.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="VAL1",
            batch_id="abbind_val1_core_v1",
            job_id="val1-antibody-h-y34a",
            mutation_group_id="val1_0001",
            entity_side="antibody",
            resseq=34,
            raw_ddg_kcal_mol=15.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="VAL2",
            batch_id="abbind_val2_core_v1",
            job_id="val2-antigen-w-g94a",
            mutation_group_id="val2_0001",
            entity_side="antigen",
            resseq=94,
            raw_ddg_kcal_mol=7.5,
        ),
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    development_report = report_ab_bind_plan(plan_root, split_name="development", split_path=split_path)
    calibration_report = report_ab_bind_plan(plan_root, split_name="calibration", split_path=split_path)
    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert isclose(validation_report["benchmark_metrics"]["pearson_r"], -1.0, rel_tol=1e-9)

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        predict_split_name="validation",
        split_path=split_path,
        model="side_linear",
    )

    assert calibration_payload["fit_pair_count"] == 4
    assert calibration_payload["predict_job_count"] == 2
    assert calibration_payload["predict_pair_count"] == 2
    assert calibration_payload["fit_coverage"]["by_side"] == {"antibody": 2, "antigen": 2, "unknown": 0}
    assert calibration_payload["fit_coverage"]["experimental_effect_bin_counts"] == {
        "negative": 0,
        "near_zero": 0,
        "positive": 4,
    }
    assert calibration_payload["predict_raw_coverage"]["by_side"] == {"antibody": 1, "antigen": 1, "unknown": 0}
    assert calibration_payload["predict_calibrated_coverage"]["experimental_effect_bin_counts"] == {
        "negative": 0,
        "near_zero": 0,
        "positive": 2,
    }
    assert isclose(calibration_payload["raw_metrics"]["pearson_r"], -1.0, rel_tol=1e-9)
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert calibration_payload["predict_target_exclusion_policy"]["target_field"] == "complex_id"
    assert calibration_payload["predict_raw_target_excluded_complex_ids"] == []
    assert calibration_payload["predict_calibrated_target_excluded_complex_ids"] == []
    assert calibration_payload["predict_outlier_trim_policy"]["outlier_trim_method"] == "tukey_iqr"
    assert calibration_payload["predict_raw_outlier_trimmed_metrics"]["paired_job_count"] == 2
    assert calibration_payload["predict_calibrated_outlier_trimmed_metrics"]["paired_job_count"] == 2
    assert calibration_payload["predict_raw_target_filtered_metrics"]["paired_job_count"] == 2
    assert calibration_payload["predict_calibrated_target_filtered_metrics"]["paired_job_count"] == 2
    assert calibration_payload["predict_raw_target_filtered_outlier_trimmed_metrics"]["paired_job_count"] == 2
    assert calibration_payload["predict_calibrated_target_filtered_outlier_trimmed_metrics"]["paired_job_count"] == 2
    assert isclose(
        calibration_payload["predict_raw_target_filtered_metrics"]["pearson_r"],
        calibration_payload["raw_metrics"]["pearson_r"],
        rel_tol=1e-9,
    )
    assert isclose(
        calibration_payload["predict_calibrated_target_filtered_metrics"]["pearson_r"],
        calibration_payload["calibrated_metrics"]["pearson_r"],
        rel_tol=1e-9,
    )
    assert {row["complex_id"] for row in calibration_payload["predict_raw_target_metrics"]} == {"VAL1", "VAL2"}
    assert {row["complex_id"] for row in calibration_payload["predict_calibrated_target_metrics"]} == {"VAL1", "VAL2"}
    assert {row["complex_id"] for row in calibration_payload["predict_raw_target_outlier_trim_metrics"]} == {
        "VAL1",
        "VAL2",
    }
    assert {row["complex_id"] for row in calibration_payload["predict_calibrated_target_outlier_trim_metrics"]} == {
        "VAL1",
        "VAL2",
    }
    assert isclose(calibration_payload["model"]["groups"]["antibody"]["slope"], 0.1, abs_tol=1e-6)
    assert isclose(calibration_payload["model"]["groups"]["antibody"]["intercept"], 0.0, abs_tol=1e-6)
    assert isclose(calibration_payload["model"]["groups"]["antigen"]["slope"], 0.4, abs_tol=1e-6)
    assert isclose(calibration_payload["model"]["groups"]["antigen"]["intercept"], 0.0, abs_tol=1e-6)

    reports_dir = Path(calibration_payload["reports_dir"])
    assert reports_dir.is_dir()
    assert len(read_csv_rows(reports_dir / "fit_pairs.csv")) == 4

    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    assert calibrated_jobs["val1-antibody-h-y34a"]["calibration_group"] == "antibody"
    assert calibrated_jobs["val2-antigen-w-g94a"]["calibration_group"] == "antigen"
    assert isclose(float(calibrated_jobs["val1-antibody-h-y34a"]["calibrated_ddg_kcal_mol"]), 1.5, abs_tol=1e-6)
    assert isclose(float(calibrated_jobs["val2-antigen-w-g94a"]["calibrated_ddg_kcal_mol"]), 3.0, abs_tol=1e-6)

    calibrated_pairs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_pairs_calibrated.csv")}
    assert isclose(float(calibrated_pairs["val1-antibody-h-y34a"]["predicted_ddg_kcal_mol"]), 1.5, abs_tol=1e-6)
    assert isclose(float(calibrated_pairs["val2-antigen-w-g94a"]["predicted_ddg_kcal_mol"]), 3.0, abs_tol=1e-6)
    raw_target_metrics = {row["complex_id"]: row for row in read_csv_rows(reports_dir / "predict_target_metrics_raw.csv")}
    calibrated_target_metrics = {
        row["complex_id"]: row for row in read_csv_rows(reports_dir / "predict_target_metrics_calibrated.csv")
    }
    assert set(raw_target_metrics) == {"VAL1", "VAL2"}
    assert set(calibrated_target_metrics) == {"VAL1", "VAL2"}
    assert all(row["excluded_from_target_filtered_metrics"] == "False" for row in raw_target_metrics.values())
    assert all(row["excluded_from_target_filtered_metrics"] == "False" for row in calibrated_target_metrics.values())
    assert len(read_csv_rows(reports_dir / "predict_target_outlier_trim_metrics_raw.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_target_outlier_trim_metrics_calibrated.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_pairs_raw_outlier_trimmed.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_pairs_calibrated_outlier_trimmed.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_pairs_raw_target_filtered.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_pairs_calibrated_target_filtered.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_pairs_raw_target_filtered_outlier_trimmed.csv")) == 2
    assert len(read_csv_rows(reports_dir / "predict_pairs_calibrated_target_filtered_outlier_trimmed.csv")) == 2


def test_calibrate_ab_bind_plan_fits_quadratic_model_across_multiple_fit_splits(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "dev1_0001",
                "ddg_kcal_mol": "4.0",
                "source_mutation": "W:G80A",
                "mutation_tokens": "W:G80A@antigen",
            },
            {
                "mutation_group_id": "dev2_0001",
                "ddg_kcal_mol": "1.0",
                "source_mutation": "W:G81A",
                "mutation_tokens": "W:G81A@antigen",
            },
            {
                "mutation_group_id": "cal1_0001",
                "ddg_kcal_mol": "0.0",
                "source_mutation": "W:G82A",
                "mutation_tokens": "W:G82A@antigen",
            },
            {
                "mutation_group_id": "cal2_0001",
                "ddg_kcal_mol": "1.0",
                "source_mutation": "W:G83A",
                "mutation_tokens": "W:G83A@antigen",
            },
            {
                "mutation_group_id": "cal3_0001",
                "ddg_kcal_mol": "4.0",
                "source_mutation": "W:G84A",
                "mutation_tokens": "W:G84A@antigen",
            },
            {
                "mutation_group_id": "val1_0001",
                "ddg_kcal_mol": "2.25",
                "source_mutation": "W:G85A",
                "mutation_tokens": "W:G85A@antigen",
            },
            {
                "mutation_group_id": "val2_0001",
                "ddg_kcal_mol": "0.25",
                "source_mutation": "W:G86A",
                "mutation_tokens": "W:G86A@antigen",
            },
            {
                "mutation_group_id": "val3_0001",
                "ddg_kcal_mol": "2.25",
                "source_mutation": "W:G87A",
                "mutation_tokens": "W:G87A@antigen",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "development": {"complex_ids": ["DEV1", "DEV2"]},
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3"]},
                "validation": {"complex_ids": ["VAL1", "VAL2", "VAL3"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(job_id: str, mutation_group_id: str, *, resseq: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antigen",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "W",
                        "resseq": resseq,
                        "icode": "",
                        "wt": "G",
                        "mut": "A",
                        "entity_side": "antigen",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, resseq=resseq)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=0.05,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.05)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id="DEV1",
            batch_id="abbind_dev1_core_v1",
            job_id="dev1-antigen-w-g80a",
            mutation_group_id="dev1_0001",
            resseq=80,
            raw_ddg_kcal_mol=-2.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="DEV2",
            batch_id="abbind_dev2_core_v1",
            job_id="dev2-antigen-w-g81a",
            mutation_group_id="dev2_0001",
            resseq=81,
            raw_ddg_kcal_mol=-1.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="CAL1",
            batch_id="abbind_cal1_core_v1",
            job_id="cal1-antigen-w-g82a",
            mutation_group_id="cal1_0001",
            resseq=82,
            raw_ddg_kcal_mol=0.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="CAL2",
            batch_id="abbind_cal2_core_v1",
            job_id="cal2-antigen-w-g83a",
            mutation_group_id="cal2_0001",
            resseq=83,
            raw_ddg_kcal_mol=1.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="CAL3",
            batch_id="abbind_cal3_core_v1",
            job_id="cal3-antigen-w-g84a",
            mutation_group_id="cal3_0001",
            resseq=84,
            raw_ddg_kcal_mol=2.0,
        ),
        write_ready_job(
            plan_root,
            complex_id="VAL1",
            batch_id="abbind_val1_core_v1",
            job_id="val1-antigen-w-g85a",
            mutation_group_id="val1_0001",
            resseq=85,
            raw_ddg_kcal_mol=-1.5,
        ),
        write_ready_job(
            plan_root,
            complex_id="VAL2",
            batch_id="abbind_val2_core_v1",
            job_id="val2-antigen-w-g86a",
            mutation_group_id="val2_0001",
            resseq=86,
            raw_ddg_kcal_mol=0.5,
        ),
        write_ready_job(
            plan_root,
            complex_id="VAL3",
            batch_id="abbind_val3_core_v1",
            job_id="val3-antigen-w-g87a",
            mutation_group_id="val3_0001",
            resseq=87,
            raw_ddg_kcal_mol=1.5,
        ),
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert float(validation_report["benchmark_metrics"]["pearson_r"]) < 0.0

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="quadratic",
    )

    assert calibration_payload["fit_pair_count"] == 5
    assert calibration_payload["predict_job_count"] == 3
    assert calibration_payload["predict_pair_count"] == 3
    assert calibration_payload["fit_split_names"] == ["development", "calibration"]
    assert calibration_payload["fit_coverage"]["by_side"] == {"antibody": 0, "antigen": 5, "unknown": 0}
    assert calibration_payload["fit_coverage"]["experimental_effect_bin_counts"] == {
        "negative": 0,
        "near_zero": 1,
        "positive": 4,
    }
    assert calibration_payload["predict_raw_coverage"]["experimental_effect_bin_counts"] == {
        "negative": 0,
        "near_zero": 1,
        "positive": 2,
    }
    assert calibration_payload["calibrated_metrics"]["pearson_r"] > calibration_payload["raw_metrics"]["pearson_r"]
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert isclose(calibration_payload["model"]["groups"]["global"]["intercept"], 0.0, abs_tol=1e-6)
    assert isclose(calibration_payload["model"]["groups"]["global"]["linear_coefficient"], 0.0, abs_tol=1e-6)
    assert isclose(calibration_payload["model"]["groups"]["global"]["quadratic_coefficient"], 1.0, abs_tol=1e-6)

    reports_dir = Path(calibration_payload["reports_dir"])
    assert reports_dir.is_dir()
    assert len(read_csv_rows(reports_dir / "fit_pairs.csv")) == 5

    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    assert calibrated_jobs["val1-antigen-w-g85a"]["calibration_family"] == "quadratic"
    assert isclose(float(calibrated_jobs["val1-antigen-w-g85a"]["calibrated_ddg_kcal_mol"]), 2.25, abs_tol=1e-6)
    assert isclose(float(calibrated_jobs["val2-antigen-w-g86a"]["calibrated_ddg_kcal_mol"]), 0.25, abs_tol=1e-6)
    assert isclose(float(calibrated_jobs["val3-antigen-w-g87a"]["calibrated_ddg_kcal_mol"]), 2.25, abs_tol=1e-6)


def test_calibrate_ab_bind_plan_fits_stderr_quadratic_model(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    def expected_ddg(raw_ddg_kcal_mol: float, bar_stderr_kcal_mol: float) -> float:
        return (
            1.0
            + 0.2 * raw_ddg_kcal_mol
            - 0.05 * raw_ddg_kcal_mol * raw_ddg_kcal_mol
            + 0.3 * bar_stderr_kcal_mol
            + 0.4 * log1p(bar_stderr_kcal_mol)
        )

    fit_specs = [
        ("DEV1", "dev1_0001", 80, -2.0, 0.5),
        ("DEV2", "dev2_0001", 81, -1.0, 1.0),
        ("CAL1", "cal1_0001", 82, 0.0, 1.5),
        ("CAL2", "cal2_0001", 83, 1.5, 2.0),
        ("CAL3", "cal3_0001", 84, 2.5, 3.0),
    ]
    validation_specs = [
        ("VAL1", "val1_0001", 85, -0.5, 0.75),
        ("VAL2", "val2_0001", 86, 1.0, 2.5),
        ("VAL3", "val3_0001", 87, 2.0, 1.2),
    ]

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": mutation_group_id,
                "ddg_kcal_mol": f"{expected_ddg(raw_ddg, bar_stderr):.6f}",
                "source_mutation": f"W:G{resseq}A",
                "mutation_tokens": f"W:G{resseq}A@antigen",
            }
            for _complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "development": {"complex_ids": ["DEV1", "DEV2"]},
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3"]},
                "validation": {"complex_ids": ["VAL1", "VAL2", "VAL3"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(job_id: str, mutation_group_id: str, *, resseq: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antigen",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "W",
                        "resseq": resseq,
                        "icode": "",
                        "wt": "G",
                        "mut": "A",
                        "entity_side": "antigen",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
        bar_stderr_kcal_mol: float,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, resseq=resseq)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=bar_stderr_kcal_mol / kcal_per_kt,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.0)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id=complex_id,
            batch_id=f"abbind_{complex_id.lower()}_core_v1",
            job_id=f"{complex_id.lower()}-antigen-w-g{resseq}a",
            mutation_group_id=mutation_group_id,
            resseq=resseq,
            raw_ddg_kcal_mol=raw_ddg,
            bar_stderr_kcal_mol=bar_stderr,
        )
        for complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    development_report = report_ab_bind_plan(plan_root, split_name="development", split_path=split_path)
    calibration_report = report_ab_bind_plan(plan_root, split_name="calibration", split_path=split_path)
    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert float(validation_report["benchmark_metrics"]["pearson_r"]) < 1.0

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="stderr_quadratic",
    )

    assert calibration_payload["fit_pair_count"] == 5
    assert calibration_payload["predict_pair_count"] == 3
    assert calibration_payload["calibrated_metrics"]["pearson_r"] > calibration_payload["raw_metrics"]["pearson_r"]
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert calibration_payload["model"]["groups"]["global"]["family"] == "stderr_quadratic"
    assert isclose(calibration_payload["model"]["groups"]["global"]["intercept"], 1.0, abs_tol=3e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["linear_coefficient"], 0.2, abs_tol=3e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["quadratic_coefficient"], -0.05, abs_tol=3e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["stderr_coefficient"], 0.3, abs_tol=3e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["log_stderr_coefficient"], 0.4, abs_tol=3e-5)

    reports_dir = Path(calibration_payload["reports_dir"])
    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    for complex_id, _mutation_group_id, resseq, raw_ddg, bar_stderr in validation_specs:
        job_id = f"{complex_id.lower()}-antigen-w-g{resseq}a"
        expected = expected_ddg(raw_ddg, bar_stderr)
        assert calibrated_jobs[job_id]["calibration_family"] == "stderr_quadratic"
        assert isclose(float(calibrated_jobs[job_id]["calibration_input_ddg_bar_stderr_kcal_mol"]), bar_stderr, abs_tol=1e-6)
        assert isclose(float(calibrated_jobs[job_id]["calibrated_ddg_kcal_mol"]), expected, abs_tol=3e-6)


def test_calibrate_ab_bind_plan_fits_logabs_stderr_quadratic_model(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    def expected_ddg(raw_ddg_kcal_mol: float, bar_stderr_kcal_mol: float) -> float:
        logabs = log1p(abs(raw_ddg_kcal_mol))
        return (
            1.2
            + 2.4 * logabs
            - 0.35 * logabs * logabs
            + 0.15 * bar_stderr_kcal_mol
            - 0.8 * log1p(bar_stderr_kcal_mol)
        )

    fit_specs = [
        ("DEV1", "dev1_0001", 80, -4.0, 0.4),
        ("DEV2", "dev2_0001", 81, -1.5, 1.1),
        ("CAL1", "cal1_0001", 82, 0.0, 1.8),
        ("CAL2", "cal2_0001", 83, 2.0, 0.7),
        ("CAL3", "cal3_0001", 84, 5.0, 2.4),
    ]
    validation_specs = [
        ("VAL1", "val1_0001", 85, -3.0, 0.9),
        ("VAL2", "val2_0001", 86, 1.0, 2.1),
        ("VAL3", "val3_0001", 87, 4.0, 1.3),
    ]

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": mutation_group_id,
                "ddg_kcal_mol": f"{expected_ddg(raw_ddg, bar_stderr):.6f}",
                "source_mutation": f"W:G{resseq}A",
                "mutation_tokens": f"W:G{resseq}A@antigen",
            }
            for _complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "development": {"complex_ids": ["DEV1", "DEV2"]},
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3"]},
                "validation": {"complex_ids": ["VAL1", "VAL2", "VAL3"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(job_id: str, mutation_group_id: str, *, resseq: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antigen",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "W",
                        "resseq": resseq,
                        "icode": "",
                        "wt": "G",
                        "mut": "A",
                        "entity_side": "antigen",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
        bar_stderr_kcal_mol: float,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, resseq=resseq)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=bar_stderr_kcal_mol / kcal_per_kt,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.0)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id=complex_id,
            batch_id=f"abbind_{complex_id.lower()}_core_v1",
            job_id=f"{complex_id.lower()}-antigen-w-g{resseq}a",
            mutation_group_id=mutation_group_id,
            resseq=resseq,
            raw_ddg_kcal_mol=raw_ddg,
            bar_stderr_kcal_mol=bar_stderr,
        )
        for complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert float(validation_report["benchmark_metrics"]["pearson_r"]) < 1.0

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="logabs_stderr_quadratic",
    )

    assert calibration_payload["fit_pair_count"] == 5
    assert calibration_payload["predict_pair_count"] == 3
    assert calibration_payload["calibrated_metrics"]["pearson_r"] > calibration_payload["raw_metrics"]["pearson_r"]
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert calibration_payload["model"]["groups"]["global"]["family"] == "logabs_stderr_quadratic"
    assert isclose(calibration_payload["model"]["groups"]["global"]["intercept"], 1.2, abs_tol=5e-5)
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["logabs_linear_coefficient"],
        2.4,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["logabs_quadratic_coefficient"],
        -0.35,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["stderr_coefficient"],
        0.15,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["log_stderr_coefficient"],
        -0.8,
        abs_tol=5e-5,
    )

    reports_dir = Path(calibration_payload["reports_dir"])
    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    for complex_id, _mutation_group_id, resseq, raw_ddg, bar_stderr in validation_specs:
        job_id = f"{complex_id.lower()}-antigen-w-g{resseq}a"
        expected = expected_ddg(raw_ddg, bar_stderr)
        assert calibrated_jobs[job_id]["calibration_family"] == "logabs_stderr_quadratic"
        assert isclose(float(calibrated_jobs[job_id]["calibration_input_ddg_bar_stderr_kcal_mol"]), bar_stderr, abs_tol=1e-6)
        assert isclose(float(calibrated_jobs[job_id]["calibrated_ddg_kcal_mol"]), expected, abs_tol=1e-5)


def test_calibrate_ab_bind_plan_fits_expdecay_invstderr_quadratic_model(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    expdecay_rate = 0.6

    def expected_ddg(raw_ddg_kcal_mol: float, bar_stderr_kcal_mol: float) -> float:
        expdecay = 1.0 - exp(-expdecay_rate * log1p(abs(raw_ddg_kcal_mol)))
        inv_stderr = 1.0 / (1.0 + bar_stderr_kcal_mol)
        return 0.9 + 3.1 * expdecay - 1.25 * expdecay * expdecay + 0.55 * inv_stderr

    fit_specs = [
        ("DEV1", "dev1_0001", 90, -5.0, 0.4),
        ("DEV2", "dev2_0001", 91, -1.0, 1.5),
        ("CAL1", "cal1_0001", 92, 0.5, 2.2),
        ("CAL2", "cal2_0001", 93, 3.0, 0.8),
        ("CAL3", "cal3_0001", 94, 7.0, 1.9),
    ]
    validation_specs = [
        ("VAL1", "val1_0001", 95, -3.0, 0.9),
        ("VAL2", "val2_0001", 96, 1.5, 2.5),
        ("VAL3", "val3_0001", 97, 5.0, 1.1),
    ]

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": mutation_group_id,
                "ddg_kcal_mol": f"{expected_ddg(raw_ddg, bar_stderr):.6f}",
                "source_mutation": f"W:G{resseq}A",
                "mutation_tokens": f"W:G{resseq}A@antigen",
            }
            for _complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "development": {"complex_ids": ["DEV1", "DEV2"]},
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3"]},
                "validation": {"complex_ids": ["VAL1", "VAL2", "VAL3"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(job_id: str, mutation_group_id: str, *, resseq: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antigen",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "W",
                        "resseq": resseq,
                        "icode": "",
                        "wt": "G",
                        "mut": "A",
                        "entity_side": "antigen",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
        bar_stderr_kcal_mol: float,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, resseq=resseq)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=bar_stderr_kcal_mol / kcal_per_kt,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.0)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id=complex_id,
            batch_id=f"abbind_{complex_id.lower()}_core_v1",
            job_id=f"{complex_id.lower()}-antigen-w-g{resseq}a",
            mutation_group_id=mutation_group_id,
            resseq=resseq,
            raw_ddg_kcal_mol=raw_ddg,
            bar_stderr_kcal_mol=bar_stderr,
        )
        for complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert float(validation_report["benchmark_metrics"]["pearson_r"]) < 1.0

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="expdecay_invstderr_quadratic",
    )

    assert calibration_payload["fit_pair_count"] == 5
    assert calibration_payload["predict_pair_count"] == 3
    assert calibration_payload["calibrated_metrics"]["pearson_r"] > calibration_payload["raw_metrics"]["pearson_r"]
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert calibration_payload["model"]["groups"]["global"]["family"] == "expdecay_invstderr_quadratic"
    assert isclose(calibration_payload["model"]["groups"]["global"]["intercept"], 0.9, abs_tol=5e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["expdecay_rate"], expdecay_rate, abs_tol=1e-12)
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["expdecay_linear_coefficient"],
        3.1,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["expdecay_quadratic_coefficient"],
        -1.25,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["inv_stderr_coefficient"],
        0.55,
        abs_tol=5e-5,
    )

    reports_dir = Path(calibration_payload["reports_dir"])
    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    for complex_id, _mutation_group_id, resseq, raw_ddg, bar_stderr in validation_specs:
        job_id = f"{complex_id.lower()}-antigen-w-g{resseq}a"
        expected = expected_ddg(raw_ddg, bar_stderr)
        assert calibrated_jobs[job_id]["calibration_family"] == "expdecay_invstderr_quadratic"
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_expdecay_rate"]),
            expdecay_rate,
            abs_tol=1e-12,
        )
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_input_ddg_bar_stderr_kcal_mol"]),
            bar_stderr,
            abs_tol=1e-6,
        )
        assert isclose(float(calibrated_jobs[job_id]["calibrated_ddg_kcal_mol"]), expected, abs_tol=1e-5)


def test_calibrate_ab_bind_plan_fits_hill_invstderr_quadratic_model(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    hill_exponent = 2.7
    hill_midpoint = 1.5

    def expected_ddg(raw_ddg_kcal_mol: float, bar_stderr_kcal_mol: float) -> float:
        logabs = log1p(abs(raw_ddg_kcal_mol))
        hill_numerator = logabs**hill_exponent
        hill = hill_numerator / (hill_numerator + hill_midpoint)
        inv_stderr = 1.0 / (1.0 + bar_stderr_kcal_mol)
        return 1.05 + 2.75 * hill - 1.1 * hill * hill + 0.6 * inv_stderr

    fit_specs = [
        ("DEV1", "dev1_0001", 100, -5.5, 0.4),
        ("DEV2", "dev2_0001", 101, -1.2, 1.6),
        ("CAL1", "cal1_0001", 102, 0.5, 2.0),
        ("CAL2", "cal2_0001", 103, 3.5, 0.9),
        ("CAL3", "cal3_0001", 104, 7.5, 1.8),
    ]
    validation_specs = [
        ("VAL1", "val1_0001", 105, -3.0, 0.8),
        ("VAL2", "val2_0001", 106, 1.7, 2.4),
        ("VAL3", "val3_0001", 107, 5.5, 1.2),
    ]

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": mutation_group_id,
                "ddg_kcal_mol": f"{expected_ddg(raw_ddg, bar_stderr):.6f}",
                "source_mutation": f"W:G{resseq}A",
                "mutation_tokens": f"W:G{resseq}A@antigen",
            }
            for _complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "development": {"complex_ids": ["DEV1", "DEV2"]},
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3"]},
                "validation": {"complex_ids": ["VAL1", "VAL2", "VAL3"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(job_id: str, mutation_group_id: str, *, resseq: int) -> dict:
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": "antigen",
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": "W",
                        "resseq": resseq,
                        "icode": "",
                        "wt": "G",
                        "mut": "A",
                        "entity_side": "antigen",
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
        bar_stderr_kcal_mol: float,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, resseq=resseq)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=bar_stderr_kcal_mol / kcal_per_kt,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.0)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id=complex_id,
            batch_id=f"abbind_{complex_id.lower()}_core_v1",
            job_id=f"{complex_id.lower()}-antigen-w-g{resseq}a",
            mutation_group_id=mutation_group_id,
            resseq=resseq,
            raw_ddg_kcal_mol=raw_ddg,
            bar_stderr_kcal_mol=bar_stderr,
        )
        for complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr in [*fit_specs, *validation_specs]
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert float(validation_report["benchmark_metrics"]["pearson_r"]) < 1.0

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="hill_invstderr_quadratic",
    )

    assert calibration_payload["fit_pair_count"] == 5
    assert calibration_payload["predict_pair_count"] == 3
    assert calibration_payload["calibrated_metrics"]["pearson_r"] > calibration_payload["raw_metrics"]["pearson_r"]
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert calibration_payload["model"]["groups"]["global"]["family"] == "hill_invstderr_quadratic"
    assert isclose(calibration_payload["model"]["groups"]["global"]["intercept"], 1.05, abs_tol=5e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["hill_exponent"], hill_exponent, abs_tol=1e-12)
    assert isclose(calibration_payload["model"]["groups"]["global"]["hill_midpoint"], hill_midpoint, abs_tol=1e-12)
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["hill_linear_coefficient"],
        2.75,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["hill_quadratic_coefficient"],
        -1.1,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["inv_stderr_coefficient"],
        0.6,
        abs_tol=5e-5,
    )

    reports_dir = Path(calibration_payload["reports_dir"])
    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    for complex_id, _mutation_group_id, resseq, raw_ddg, bar_stderr in validation_specs:
        job_id = f"{complex_id.lower()}-antigen-w-g{resseq}a"
        expected = expected_ddg(raw_ddg, bar_stderr)
        assert calibrated_jobs[job_id]["calibration_family"] == "hill_invstderr_quadratic"
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_hill_exponent"]),
            hill_exponent,
            abs_tol=1e-12,
        )
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_hill_midpoint"]),
            hill_midpoint,
            abs_tol=1e-12,
        )
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_input_ddg_bar_stderr_kcal_mol"]),
            bar_stderr,
            abs_tol=1e-6,
        )
        assert isclose(float(calibrated_jobs[job_id]["calibrated_ddg_kcal_mol"]), expected, abs_tol=1e-5)


def test_calibrate_ab_bind_plan_fits_hill_side_invstderr_quadratic_model(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    curated_dir = benchmark_root / "curated"
    splits_dir = benchmark_root / "splits"
    curated_dir.mkdir(parents=True)
    splits_dir.mkdir(parents=True)

    hill_exponent = 3.0
    hill_midpoint = 1.65
    antibody_side_coefficient = 1.35

    def expected_ddg(raw_ddg_kcal_mol: float, bar_stderr_kcal_mol: float, *, entity_side: str) -> float:
        logabs = log1p(abs(raw_ddg_kcal_mol))
        hill_numerator = logabs**hill_exponent
        hill = hill_numerator / (hill_numerator + hill_midpoint)
        inv_stderr = 1.0 / (1.0 + bar_stderr_kcal_mol)
        side_indicator = 1.0 if entity_side == "antibody" else 0.0
        return 0.95 + 2.4 * hill - 0.8 * hill * hill + 0.7 * inv_stderr + antibody_side_coefficient * side_indicator

    fit_specs = [
        ("DEV1", "dev1_0001", 110, -5.5, 0.4, "antigen"),
        ("DEV2", "dev2_0001", 111, -1.2, 1.6, "antibody"),
        ("DEV3", "dev3_0001", 112, 0.5, 2.0, "antigen"),
        ("CAL1", "cal1_0001", 113, 3.5, 0.9, "antibody"),
        ("CAL2", "cal2_0001", 114, 7.5, 1.8, "antigen"),
        ("CAL3", "cal3_0001", 115, -2.7, 0.6, "antibody"),
    ]
    validation_specs = [
        ("VAL1", "val1_0001", 116, -3.0, 0.8, "antibody"),
        ("VAL2", "val2_0001", 117, 1.7, 2.4, "antigen"),
        ("VAL3", "val3_0001", 118, 5.5, 1.2, "antibody"),
    ]

    def chain_for_side(entity_side: str) -> str:
        return "H" if entity_side == "antibody" else "W"

    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": mutation_group_id,
                "ddg_kcal_mol": f"{expected_ddg(raw_ddg, bar_stderr, entity_side=entity_side):.6f}",
                "source_mutation": f"{chain_for_side(entity_side)}:G{resseq}A",
                "mutation_tokens": f"{chain_for_side(entity_side)}:G{resseq}A@{entity_side}",
            }
            for _complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr, entity_side in [
                *fit_specs,
                *validation_specs,
            ]
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    split_path = splits_dir / "ab_bind_rbfe_core_v1_split_v1.yml"
    write_yaml(
        split_path,
        {
            "spec_name": "core_v1",
            "splits": {
                "development": {"complex_ids": ["DEV1", "DEV2", "DEV3"]},
                "calibration": {"complex_ids": ["CAL1", "CAL2", "CAL3"]},
                "validation": {"complex_ids": ["VAL1", "VAL2", "VAL3"]},
            },
        },
    )

    kcal_per_kt = 310.0 * 0.00198720425864083

    def job_spec(job_id: str, mutation_group_id: str, *, resseq: int, entity_side: str) -> dict:
        chain_id = chain_for_side(entity_side)
        return {
            "job_id": job_id,
            "batch_id": "unused",
            "mutation_group": {
                "mutation_group_id": mutation_group_id,
                "mutation_count": 1,
                "entity_side": entity_side,
                "charge_conserving": True,
                "min_version": "v1",
                "sites": [
                    {
                        "chain_id": chain_id,
                        "resseq": resseq,
                        "icode": "",
                        "wt": "G",
                        "mut": "A",
                        "entity_side": entity_side,
                    }
                ],
            },
            "protocol": {
                "preset": "single_point",
                "temperature_k": 310.0,
                "lambda_windows": 2,
                "repeats": 1,
                "production_ps": 10,
                "overlap_threshold": 0.2,
                "max_repeat_delta_kcal_mol": 10.0,
                "max_bar_stderr_kcal_mol": 10.0,
            },
            "system": {
                "system_name": job_id.split("-")[0],
                "structure_source": "experimental",
            },
        }

    def write_ready_job(
        plan_root: Path,
        *,
        complex_id: str,
        batch_id: str,
        job_id: str,
        mutation_group_id: str,
        resseq: int,
        raw_ddg_kcal_mol: float,
        bar_stderr_kcal_mol: float,
        entity_side: str,
    ) -> dict[str, str]:
        batch_dir = plan_root / batch_id
        job_dir = batch_dir / "jobs" / job_id
        ensure_dir(job_dir / "stages")
        spec = job_spec(job_id, mutation_group_id, resseq=resseq, entity_side=entity_side)
        spec["batch_id"] = batch_id
        write_json(job_dir / "job_spec.json", spec)
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
            _write_stage(job_dir, stage)
        _write_bar_outputs(
            job_dir / "legs" / "complex" / "rep01",
            delta_kt=raw_ddg_kcal_mol / kcal_per_kt,
            stderr_kt=bar_stderr_kcal_mol / kcal_per_kt,
        )
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.0)
        return {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
        }

    plan_root = tmp_path / "runs"
    batches = [
        write_ready_job(
            plan_root,
            complex_id=complex_id,
            batch_id=f"abbind_{complex_id.lower()}_core_v1",
            job_id=f"{complex_id.lower()}-{entity_side}-{chain_for_side(entity_side).lower()}-g{resseq}a",
            mutation_group_id=mutation_group_id,
            resseq=resseq,
            raw_ddg_kcal_mol=raw_ddg,
            bar_stderr_kcal_mol=bar_stderr,
            entity_side=entity_side,
        )
        for complex_id, mutation_group_id, resseq, raw_ddg, bar_stderr, entity_side in [*fit_specs, *validation_specs]
    ]

    plan_index = {
        "benchmark_root": str(benchmark_root),
        "spec_name": "core_v1",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [batch["complex_id"] for batch in batches],
        "batches": batches,
    }
    write_json(plan_root / "plan_index.json", plan_index)
    write_yaml(plan_root / "plan_index.yml", plan_index)

    development_report = report_ab_bind_plan(plan_root, split_name="development", split_path=split_path)
    calibration_report = report_ab_bind_plan(plan_root, split_name="calibration", split_path=split_path)
    validation_report = report_ab_bind_plan(plan_root, split_name="validation", split_path=split_path)
    assert float(validation_report["benchmark_metrics"]["pearson_r"]) < 1.0

    calibration_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="hill_side_invstderr_quadratic",
    )
    reused_reports_payload = calibrate_ab_bind_plan(
        plan_root,
        fit_split_name="calibration",
        fit_split_names=["development", "calibration"],
        predict_split_name="validation",
        split_path=split_path,
        model="hill_side_invstderr_quadratic",
        fit_reports_dirs=[Path(development_report["reports_dir"]), Path(calibration_report["reports_dir"])],
        predict_reports_dir=Path(validation_report["reports_dir"]),
    )

    assert calibration_payload["fit_pair_count"] == 6
    assert calibration_payload["predict_pair_count"] == 3
    assert calibration_payload["calibrated_metrics"]["pearson_r"] > calibration_payload["raw_metrics"]["pearson_r"]
    assert isclose(calibration_payload["calibrated_metrics"]["pearson_r"], 1.0, rel_tol=1e-9)
    assert reused_reports_payload["calibrated_metrics"] == calibration_payload["calibrated_metrics"]
    assert calibration_payload["model"]["groups"]["global"]["family"] == "hill_side_invstderr_quadratic"
    assert isclose(calibration_payload["model"]["groups"]["global"]["intercept"], 0.95, abs_tol=5e-5)
    assert isclose(calibration_payload["model"]["groups"]["global"]["hill_exponent"], hill_exponent, abs_tol=1e-12)
    assert isclose(calibration_payload["model"]["groups"]["global"]["hill_midpoint"], hill_midpoint, abs_tol=1e-12)
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["hill_linear_coefficient"],
        2.4,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["hill_quadratic_coefficient"],
        -0.8,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["inv_stderr_coefficient"],
        0.7,
        abs_tol=5e-5,
    )
    assert isclose(
        calibration_payload["model"]["groups"]["global"]["antibody_side_coefficient"],
        antibody_side_coefficient,
        abs_tol=5e-5,
    )

    reports_dir = Path(calibration_payload["reports_dir"])
    calibrated_jobs = {row["job_id"]: row for row in read_csv_rows(reports_dir / "predict_jobs_calibrated.csv")}
    for complex_id, _mutation_group_id, resseq, raw_ddg, bar_stderr, entity_side in validation_specs:
        job_id = f"{complex_id.lower()}-{entity_side}-{chain_for_side(entity_side).lower()}-g{resseq}a"
        expected = expected_ddg(raw_ddg, bar_stderr, entity_side=entity_side)
        assert calibrated_jobs[job_id]["calibration_family"] == "hill_side_invstderr_quadratic"
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_hill_exponent"]),
            hill_exponent,
            abs_tol=1e-12,
        )
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_hill_midpoint"]),
            hill_midpoint,
            abs_tol=1e-12,
        )
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_antibody_side_coefficient"]),
            antibody_side_coefficient,
            abs_tol=5e-5,
        )
        assert isclose(
            float(calibrated_jobs[job_id]["calibration_input_ddg_bar_stderr_kcal_mol"]),
            bar_stderr,
            abs_tol=1e-6,
        )
        assert isclose(float(calibrated_jobs[job_id]["calibrated_ddg_kcal_mol"]), expected, abs_tol=1e-5)


def test_report_ab_bind_plan_separates_qc_qualified_metrics(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text(
        "\n".join(
            [
                "ATOM      1  N   GLY H   1      11.104  13.207   9.100  1.00 20.00           N",
                "ATOM      2  CA  GLY H   1      12.100  12.300   8.500  1.00 20.00           C",
                "ATOM      3  C   GLY H   1      13.300  13.100   8.000  1.00 20.00           C",
                "ATOM      4  O   GLY H   1      14.300  12.500   7.700  1.00 20.00           O",
                "ATOM      5  N   GLY L   1      21.204  10.207   9.100  1.00 20.00           N",
                "ATOM      6  CA  GLY L   1      22.200   9.300   8.500  1.00 20.00           C",
                "ATOM      7  C   GLY L   1      23.400  10.100   8.000  1.00 20.00           C",
                "ATOM      8  O   GLY L   1      24.400   9.500   7.700  1.00 20.00           O",
                "ATOM      9  N   GLY C   1      28.104  13.207   9.100  1.00 20.00           N",
                "ATOM     10  CA  GLY C   1      29.100  12.300   8.500  1.00 20.00           C",
                "ATOM     11  C   GLY C   1      30.300  13.100   8.000  1.00 20.00           C",
                "ATOM     12  O   GLY C   1      31.300  12.500   7.700  1.00 20.00           O",
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0003",
                "chain_id": "H",
                "resseq": 35,
                "icode": "",
                "wt": "G",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 3,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.6160333201786573",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "-1.2320666403573146",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0003",
                "ddg_kcal_mol": "0.0",
                "source_mutation": "H:G35A",
                "mutation_tokens": "H:G35A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(
        benchmark_root,
        protocol_path,
        spec_name="core_v1",
        runs_root=plan_root,
        complex_ids=["1VFB"],
    )
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    expected_by_group = {
        "1vfb_0001": (1.0, 0.0, 0.10, 0.10),
        "1vfb_0002": (0.0, 2.0, 0.10, 0.10),
        "1vfb_0003": (20.0, 0.0, 40.0, 60.0),
    }
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        spec = read_json(job_dir / "job_spec.json")
        mutation_group_id = spec["mutation_group"]["mutation_group_id"]
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        complex_delta_kt, apo_delta_kt, complex_stderr_kt, apo_stderr_kt = expected_by_group[mutation_group_id]
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=complex_delta_kt, stderr_kt=complex_stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=apo_delta_kt, stderr_kt=apo_stderr_kt)

    report_payload = report_ab_bind_plan(plan_root, complex_ids=["1VFB"])
    all_metrics = report_payload["benchmark_metrics"]
    qualified_metrics = report_payload["benchmark_metrics_qc_qualified"]

    assert report_payload["paired_job_count"] == 3
    assert report_payload["qc_qualified_pair_count"] == 2
    assert isclose(qualified_metrics["pearson_r"], 1.0, rel_tol=1e-9)
    assert qualified_metrics["paired_job_count"] == 2
    assert all_metrics["paired_job_count"] == 3
    assert report_payload["diagnostic_family_counts"]["completed"] == 2
    assert report_payload["diagnostic_family_counts"]["qc"] == 1
    assert report_payload["diagnostic_code_counts"]["qc_pass"] == 2
    assert report_payload["diagnostic_code_counts"]["qc_bar_stderr"] == 1
    assert report_payload["current_invalid_mutate_output_job_count"] == 0
    assert report_payload["validation_failure_taxonomy"]["counts"]["benchmark_qc_qualified"] == 2
    assert report_payload["validation_failure_taxonomy"]["counts"]["qc_sampling_issue"] == 1
    assert report_payload["benchmark_target_exclusion_policy"]["target_field"] == "complex_id"
    assert report_payload["benchmark_target_excluded_complex_ids"] == []
    assert report_payload["benchmark_target_excluded_complex_ids_qc_qualified"] == []
    assert report_payload["benchmark_metrics_target_filtered"]["paired_job_count"] == all_metrics["paired_job_count"]
    assert report_payload["benchmark_metrics_qc_qualified_target_filtered"]["paired_job_count"] == qualified_metrics[
        "paired_job_count"
    ]
    assert isclose(
        report_payload["benchmark_metrics_target_filtered"]["pearson_r"],
        all_metrics["pearson_r"],
        rel_tol=1e-9,
    )
    assert isclose(
        report_payload["benchmark_metrics_qc_qualified_target_filtered"]["pearson_r"],
        qualified_metrics["pearson_r"],
        rel_tol=1e-9,
    )
    assert isclose(
        report_payload["validation_gate"]["target_filtered_pearson_r"],
        all_metrics["pearson_r"],
        rel_tol=1e-9,
    )
    assert report_payload["validation_gate"]["target_filtered_excluded_complex_ids"] == []
    plan_jobs = read_csv_rows(plan_root / "reports" / "plan_jobs.csv")
    qualified_jobs = {row["job_id"] for row in plan_jobs if row["benchmark_qc_qualified"] == "True"}
    assert qualified_jobs == {"1vfb-antibody-h-y32a", "1vfb-antibody-h-v34a"}
    diagnostic_codes = {row["job_id"]: row["diagnostic_code"] for row in plan_jobs}
    validation_failure_categories = {row["job_id"]: row["validation_failure_category"] for row in plan_jobs}
    assert diagnostic_codes["1vfb-antibody-h-y32a"] == "qc_pass"
    assert diagnostic_codes["1vfb-antibody-h-v34a"] == "qc_pass"
    assert diagnostic_codes["1vfb-antibody-h-g35a"] == "qc_bar_stderr"
    assert validation_failure_categories["1vfb-antibody-h-y32a"] == "benchmark_qc_qualified"
    assert validation_failure_categories["1vfb-antibody-h-v34a"] == "benchmark_qc_qualified"
    assert validation_failure_categories["1vfb-antibody-h-g35a"] == "qc_sampling_issue"
    assert {row["job_id"]: row["current_invalid_mutate_output"] for row in plan_jobs} == {
        "1vfb-antibody-h-y32a": "False",
        "1vfb-antibody-h-v34a": "False",
        "1vfb-antibody-h-g35a": "False",
    }
    qualified_pairs = read_csv_rows(plan_root / "reports" / "benchmark_pairs_qc_qualified.csv")
    assert len(qualified_pairs) == 2
    target_metrics = read_csv_rows(plan_root / "reports" / "benchmark_target_metrics.csv")
    target_metrics_qc_qualified = read_csv_rows(plan_root / "reports" / "benchmark_target_metrics_qc_qualified.csv")
    target_metrics_by_complex = {row["complex_id"]: row for row in target_metrics}
    assert set(target_metrics_by_complex) == {"1VFB"}
    assert target_metrics_by_complex["1VFB"]["systematically_poor_target"] == "False"
    assert target_metrics_by_complex["1VFB"]["excluded_from_target_filtered_metrics"] == "False"
    assert len(target_metrics_qc_qualified) == 1
    assert len(read_csv_rows(plan_root / "reports" / "benchmark_pairs_target_filtered.csv")) == 3
    assert len(read_csv_rows(plan_root / "reports" / "benchmark_pairs_qc_qualified_target_filtered.csv")) == 2


def test_plan_ab_bind_rescues_materializes_warning_jobs_with_adjusted_protocols(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 2,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "1.5",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        if job_dir.name == "1vfb-antibody-h-y32a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)
        elif job_dir.name == "1vfb-antibody-h-v34a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=40.0)

    rescue_root = tmp_path / "rescue_runs"
    rescue_payload = plan_ab_bind_rescues(plan_root, complex_ids=["1VFB"], runs_root=rescue_root)

    assert rescue_payload["rescued_job_count"] == 2
    assert rescue_payload["window_relax_em_scale"] == 2.0
    assert rescue_payload["window_relax_md_scale"] == 2.0
    assert rescue_payload["nvt_scale"] == 2.0
    assert rescue_payload["npt_scale"] == 2.0
    rescue_rows = {row["source_job_id"]: row for row in read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")}
    assert rescue_rows["1vfb-antibody-h-y32a"]["rescue_reasons"] == "repeat_spread"
    assert rescue_rows["1vfb-antibody-h-y32a"]["primary_repeat_spread_leg"] == "complex"
    assert rescue_rows["1vfb-antibody-h-y32a"]["repeat_spread_legs"] == "complex"
    assert rescue_rows["1vfb-antibody-h-v34a"]["rescue_reasons"] == "bar_stderr"
    assert rescue_rows["1vfb-antibody-h-v34a"]["primary_repeat_spread_leg"] == ""

    plan_jobs = {row["job_id"]: row for row in read_csv_rows(plan_root / "reports" / "plan_jobs.csv")}
    assert plan_jobs["1vfb-antibody-h-y32a"]["primary_repeat_spread_leg"] == "complex"
    assert plan_jobs["1vfb-antibody-h-y32a"]["repeat_spread_legs"] == "complex"
    assert plan_jobs["1vfb-antibody-h-y32a"]["complex_repeat_spread_kcal_mol"] != ""
    assert plan_jobs["1vfb-antibody-h-y32a"]["apo_repeat_spread_kcal_mol"] != ""

    repeat_rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    stderr_rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-v34a" / "jobs" / "1vfb-antibody-h-v34a" / "job_spec.json"
    )

    assert repeat_rescue_spec["protocol"]["repeats"] == 3
    assert repeat_rescue_spec["protocol"]["lambda_windows"] == 2
    assert repeat_rescue_spec["protocol"]["production_ps"] == 20
    assert repeat_rescue_spec["protocol"]["window_relax_em_steps"] == 1000
    assert repeat_rescue_spec["protocol"]["window_relax_md_ps"] == 0.2
    assert repeat_rescue_spec["protocol"]["nvt_ps"] == 200
    assert repeat_rescue_spec["protocol"]["npt_ps"] == 2000
    assert repeat_rescue_spec["protocol"]["equilibration_restraint_schedule"] == "staged_backbone_release"
    assert repeat_rescue_spec["protocol"]["equilibration_release_npt_ps"] == 1000
    assert stderr_rescue_spec["protocol"]["repeats"] == 2
    assert stderr_rescue_spec["protocol"]["lambda_windows"] == 6
    assert stderr_rescue_spec["protocol"]["production_ps"] == 20
    assert stderr_rescue_spec["protocol"]["window_relax_em_steps"] == 500
    assert stderr_rescue_spec["protocol"]["window_relax_md_ps"] == 0.1
    assert stderr_rescue_spec["protocol"]["nvt_ps"] == 100
    assert stderr_rescue_spec["protocol"]["npt_ps"] == 1000
    assert stderr_rescue_spec["protocol"]["equilibration_restraint_schedule"] == "legacy_posres"
    assert stderr_rescue_spec["protocol"]["equilibration_release_npt_ps"] == 0


def test_plan_ab_bind_rescues_can_materialize_requested_pass_qc_outlier_jobs(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "8.0",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    job_dir = plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    rescue_root = tmp_path / "pass_qc_outlier_rescues"
    rescue_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        job_ids=["1vfb-antibody-h-y32a"],
        allow_pass_qc_outlier_rescue=True,
    )

    assert rescue_payload["rescued_job_count"] == 1
    assert rescue_payload["allow_pass_qc_outlier_rescue"] is True

    plan_jobs = {row["job_id"]: row for row in read_csv_rows(plan_root / "reports" / "plan_jobs.csv")}
    assert plan_jobs["1vfb-antibody-h-y32a"]["qc_status"] == "pass"
    assert plan_jobs["1vfb-antibody-h-y32a"]["has_active_alternate_candidate"] != "True"
    assert float(plan_jobs["1vfb-antibody-h-y32a"]["abs_ddg_error_kcal_mol"]) > 0.0

    rescue_rows = {row["source_job_id"]: row for row in read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")}
    assert rescue_rows["1vfb-antibody-h-y32a"]["rescue_reasons"] == "pass_qc_outlier"
    assert rescue_rows["1vfb-antibody-h-y32a"]["targeted_primary_repeat_spread_leg"] == "False"

    rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    assert rescue_spec["protocol"]["repeats"] == 3
    assert rescue_spec["protocol"]["lambda_windows"] == 6
    assert rescue_spec["protocol"]["production_ps"] == 20
    assert rescue_spec["protocol"]["window_relax_em_steps"] == 1000
    assert rescue_spec["protocol"]["window_relax_md_ps"] == 0.2
    assert rescue_spec["protocol"]["nvt_ps"] == 200
    assert rescue_spec["protocol"]["npt_ps"] == 2000
    assert rescue_spec["protocol"]["equilibration_restraint_schedule"] == "staged_backbone_release"
    assert rescue_spec["protocol"]["equilibration_release_npt_ps"] == 1000


def test_plan_ab_bind_rescues_can_apply_hotspot_override_to_selected_job_ids(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 2,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "1.5",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        if job_dir.name == "1vfb-antibody-h-y32a":
            _write_bar_outputs(
                job_dir / "legs" / "complex" / "rep01",
                delta_kt=1.0,
                stderr_kt=0.10,
                overlap_mode="poor",
            )
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)
        elif job_dir.name == "1vfb-antibody-h-v34a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=40.0)

    rescue_root = tmp_path / "hotspot_rescue_runs"
    rescue_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        hotspot_job_ids=["1VFB-ANTIBODY-H-Y32A"],
        hotspot_repeat_increment=3,
        hotspot_lambda_increment=7,
        hotspot_production_scale=5.0,
        hotspot_window_relax_em_scale=3.0,
        hotspot_window_relax_md_scale=3.0,
        hotspot_nvt_scale=3.0,
        hotspot_npt_scale=3.0,
    )

    assert rescue_payload["hotspot_job_ids"] == ["1vfb-antibody-h-y32a"]
    assert rescue_payload["hotspot_repeat_increment"] == 3
    assert rescue_payload["hotspot_lambda_increment"] == 7
    rescue_rows = {row["source_job_id"]: row for row in read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")}
    hotspot_row = rescue_rows["1vfb-antibody-h-y32a"]
    baseline_row = rescue_rows["1vfb-antibody-h-v34a"]

    assert hotspot_row["hotspot_override_applied"] == "True"
    assert hotspot_row["hotspot_override_match"] == "job_id"
    assert hotspot_row["hotspot_override_job_id"] == "1vfb-antibody-h-y32a"
    assert hotspot_row["hotspot_override_complex_id"] == ""
    assert hotspot_row["effective_repeat_increment"] == "3"
    assert hotspot_row["effective_lambda_increment"] == "7"
    assert hotspot_row["effective_production_scale"] == "5.0"
    assert hotspot_row["effective_window_relax_em_scale"] == "3.0"
    assert hotspot_row["effective_window_relax_md_scale"] == "3.0"
    assert hotspot_row["effective_nvt_scale"] == "3.0"
    assert hotspot_row["effective_npt_scale"] == "3.0"
    assert baseline_row["hotspot_override_applied"] == "False"
    assert baseline_row["hotspot_override_match"] == ""
    assert baseline_row["effective_repeat_increment"] == "1"
    assert baseline_row["effective_lambda_increment"] == "4"
    assert baseline_row["effective_production_scale"] == "2.0"

    hotspot_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    baseline_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-v34a" / "jobs" / "1vfb-antibody-h-v34a" / "job_spec.json"
    )

    assert hotspot_spec["protocol"]["repeats"] == 5
    assert hotspot_spec["protocol"]["lambda_windows"] == 9
    assert hotspot_spec["protocol"]["production_ps"] == 50
    assert hotspot_spec["protocol"]["window_relax_em_steps"] == 1500
    assert hotspot_spec["protocol"]["window_relax_md_ps"] == 0.3
    assert hotspot_spec["protocol"]["nvt_ps"] == 300
    assert hotspot_spec["protocol"]["npt_ps"] == 3000
    assert hotspot_spec["protocol"]["equilibration_restraint_schedule"] == "staged_backbone_release"
    assert hotspot_spec["protocol"]["equilibration_release_npt_ps"] == 1500
    assert baseline_spec["protocol"]["repeats"] == 2
    assert baseline_spec["protocol"]["lambda_windows"] == 6
    assert baseline_spec["protocol"]["production_ps"] == 20
    assert baseline_spec["protocol"]["window_relax_em_steps"] == 500
    assert baseline_spec["protocol"]["window_relax_md_ps"] == 0.1
    assert baseline_spec["protocol"]["nvt_ps"] == 100
    assert baseline_spec["protocol"]["npt_ps"] == 1000


def test_plan_ab_bind_rescues_can_target_primary_repeat_spread_leg(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 2,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "1.5",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        if job_dir.name == "1vfb-antibody-h-y32a":
            _write_bar_outputs(
                job_dir / "legs" / "complex" / "rep01",
                delta_kt=1.0,
                stderr_kt=0.10,
                overlap_mode="poor",
            )
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)
        elif job_dir.name == "1vfb-antibody-h-v34a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=40.0)

    rescue_root = tmp_path / "targeted_rescue_runs"
    rescue_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        target_primary_repeat_spread_leg=True,
    )

    assert rescue_payload["rescued_job_count"] == 2
    assert rescue_payload["target_primary_repeat_spread_leg"] is True
    assert rescue_payload["allow_targeted_leg_count_deepening"] is False

    rescue_rows = {row["source_job_id"]: row for row in read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")}
    assert rescue_rows["1vfb-antibody-h-y32a"]["rescue_reasons"] == "repeat_spread,overlap"
    assert rescue_rows["1vfb-antibody-h-y32a"]["targeted_primary_repeat_spread_leg"] == "True"
    assert rescue_rows["1vfb-antibody-h-y32a"]["allow_targeted_leg_count_deepening"] == "False"
    assert rescue_rows["1vfb-antibody-h-y32a"]["preserved_targeted_leg_counts"] == "True"
    assert rescue_rows["1vfb-antibody-h-y32a"]["target_legs"] == "complex"
    assert rescue_rows["1vfb-antibody-h-y32a"]["inherit_source_legs"] == "apo"
    assert rescue_rows["1vfb-antibody-h-v34a"]["targeted_primary_repeat_spread_leg"] == "False"
    assert rescue_rows["1vfb-antibody-h-v34a"]["allow_targeted_leg_count_deepening"] == "False"
    assert rescue_rows["1vfb-antibody-h-v34a"]["preserved_targeted_leg_counts"] == "False"
    assert rescue_rows["1vfb-antibody-h-v34a"]["target_legs"] == ""
    assert rescue_rows["1vfb-antibody-h-v34a"]["inherit_source_legs"] == "complex,apo"

    repeat_job_dir = rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a"
    stderr_job_dir = rescue_root / "abbind-rescue_1vfb-antibody-h-v34a" / "jobs" / "1vfb-antibody-h-v34a"
    repeat_rescue_spec = read_json(repeat_job_dir / "job_spec.json")
    stderr_rescue_spec = read_json(stderr_job_dir / "job_spec.json")
    repeat_rescue_config = read_json(repeat_job_dir / "config" / "rescue.json")

    assert repeat_rescue_spec["protocol"]["repeats"] == 2
    assert repeat_rescue_spec["protocol"]["lambda_windows"] == 2
    assert repeat_rescue_spec["protocol"]["production_ps"] == 20
    assert repeat_rescue_spec["protocol"]["window_relax_em_steps"] == 1000
    assert repeat_rescue_spec["protocol"]["window_relax_md_ps"] == 0.2
    assert repeat_rescue_spec["protocol"]["nvt_ps"] == 200
    assert repeat_rescue_spec["protocol"]["npt_ps"] == 2000
    assert repeat_rescue_spec["protocol"]["equilibration_restraint_schedule"] == "staged_backbone_release"
    assert repeat_rescue_spec["protocol"]["equilibration_release_npt_ps"] == 1000
    assert repeat_rescue_config["mode"] == "targeted_primary_repeat_spread_leg"
    assert repeat_rescue_config["target_legs"] == ["complex"]
    assert repeat_rescue_config["inherit_source_legs"] == ["apo"]
    assert repeat_rescue_config["allow_targeted_leg_count_deepening"] is False
    assert repeat_rescue_config["preserved_targeted_leg_counts"] is True
    assert repeat_rescue_config["source_job_id"] == "1vfb-antibody-h-y32a"
    assert repeat_rescue_config["source_job_dir"] == str((plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a").resolve())
    assert not (stderr_job_dir / "config" / "rescue.json").exists()
    assert stderr_rescue_spec["protocol"]["repeats"] == 2
    assert stderr_rescue_spec["protocol"]["lambda_windows"] == 6


def test_can_target_primary_repeat_spread_leg_accepts_dominant_primary_leg_across_both_legs() -> None:
    qc_report = {
        "repeat_spread_legs": ["complex", "apo"],
        "primary_repeat_spread_leg": "complex",
        "max_repeat_delta_kcal_mol": 1.0,
        "legs": {
            "complex": {"repeat_delta_kcal_mol_range": 5.2},
            "apo": {"repeat_delta_kcal_mol_range": 2.1},
        },
        "overlap_threshold": 0.2,
        "overlap_assessment": {
            "legs": {
                "complex": {"overlap_score_min": 0.24},
                "apo": {"overlap_score_min": 0.30},
            }
        },
    }

    assert benchmark_module._can_target_primary_repeat_spread_leg(qc_report, ["repeat_spread"])


def test_can_target_primary_repeat_spread_leg_rejects_nondominant_both_leg_repeat_spread() -> None:
    qc_report = {
        "repeat_spread_legs": ["complex", "apo"],
        "primary_repeat_spread_leg": "complex",
        "max_repeat_delta_kcal_mol": 1.0,
        "legs": {
            "complex": {"repeat_delta_kcal_mol_range": 2.2},
            "apo": {"repeat_delta_kcal_mol_range": 1.8},
        },
        "overlap_threshold": 0.2,
        "overlap_assessment": {
            "legs": {
                "complex": {"overlap_score_min": 0.24},
                "apo": {"overlap_score_min": 0.30},
            }
        },
    }

    assert not benchmark_module._can_target_primary_repeat_spread_leg(qc_report, ["repeat_spread"])


def test_plan_ab_bind_rescues_can_deepen_targeted_primary_repeat_spread_leg_counts(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    job_dir = plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)
    _write_bar_outputs(
        job_dir / "legs" / "complex" / "rep01",
        delta_kt=1.0,
        stderr_kt=0.10,
        overlap_mode="poor",
    )
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    rescue_root = tmp_path / "targeted_deepened_rescue_runs"
    rescue_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        target_primary_repeat_spread_leg=True,
        allow_targeted_leg_count_deepening=True,
        repeat_increment=2,
        lambda_increment=5,
    )

    assert rescue_payload["rescued_job_count"] == 1
    assert rescue_payload["target_primary_repeat_spread_leg"] is True
    assert rescue_payload["allow_targeted_leg_count_deepening"] is True

    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")
    assert len(rescue_rows) == 1
    assert rescue_rows[0]["source_job_id"] == "1vfb-antibody-h-y32a"
    assert rescue_rows[0]["targeted_primary_repeat_spread_leg"] == "True"
    assert rescue_rows[0]["allow_targeted_leg_count_deepening"] == "True"
    assert rescue_rows[0]["preserved_targeted_leg_counts"] == "False"
    assert rescue_rows[0]["rescued_repeats"] == "4"
    assert rescue_rows[0]["rescued_lambda_windows"] == "7"

    rescue_job_dir = rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a"
    rescue_spec = read_json(rescue_job_dir / "job_spec.json")
    rescue_config = read_json(rescue_job_dir / "config" / "rescue.json")

    assert rescue_spec["protocol"]["repeats"] == 4
    assert rescue_spec["protocol"]["lambda_windows"] == 7
    assert rescue_spec["protocol"]["production_ps"] == 20
    assert rescue_spec["protocol"]["window_relax_em_steps"] == 1000
    assert rescue_spec["protocol"]["window_relax_md_ps"] == 0.2
    assert rescue_spec["protocol"]["nvt_ps"] == 200
    assert rescue_spec["protocol"]["npt_ps"] == 2000
    assert rescue_config["mode"] == "targeted_primary_repeat_spread_leg"
    assert rescue_config["target_legs"] == ["complex"]
    assert rescue_config["inherit_source_legs"] == ["apo"]
    assert rescue_config["allow_targeted_leg_count_deepening"] is True
    assert rescue_config["preserved_targeted_leg_counts"] is False


def test_plan_ab_bind_rescues_rejects_target_primary_leg_when_overlap_hits_other_leg(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    job_dir = plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(
        job_dir / "legs" / "apo" / "rep01",
        delta_kt=0.0,
        stderr_kt=0.10,
        overlap_mode="poor",
    )
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    rescue_root = tmp_path / "targeted_only_rescue_runs"
    rescue_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        target_primary_repeat_spread_leg=True,
        require_target_primary_repeat_spread_leg=True,
    )

    assert rescue_payload["rescued_job_count"] == 0
    rescue_rows_path = rescue_root / "reports" / "rescue_candidates.csv"
    assert not rescue_rows_path.exists() or read_csv_rows(rescue_rows_path) == []


def test_plan_ab_bind_rescues_can_require_target_primary_repeat_spread_leg(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 2,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "1.5",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        if job_dir.name == "1vfb-antibody-h-y32a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)
        elif job_dir.name == "1vfb-antibody-h-v34a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=40.0)

    rescue_root = tmp_path / "targeted_only_rescue_runs"
    rescue_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        target_primary_repeat_spread_leg=True,
        require_target_primary_repeat_spread_leg=True,
    )

    assert rescue_payload["rescued_job_count"] == 1
    assert rescue_payload["target_primary_repeat_spread_leg"] is True
    assert rescue_payload["require_target_primary_repeat_spread_leg"] is True

    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")
    assert len(rescue_rows) == 1
    assert rescue_rows[0]["source_job_id"] == "1vfb-antibody-h-y32a"
    assert rescue_rows[0]["targeted_primary_repeat_spread_leg"] == "True"
    assert not (rescue_root / "abbind-rescue_1vfb-antibody-h-v34a").exists()


def test_plan_ab_bind_rescues_can_source_merged_winner_from_extra_plan_root(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    priority_protocol = tmp_path / "protocol.priority.yml"
    robust_protocol = tmp_path / "protocol.robust.yml"
    write_yaml(
        priority_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )
    write_yaml(
        robust_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 6,
            "repeats": 4,
            "production_ps": 20,
            "temperature_k": 310.0,
        },
    )

    priority_root = tmp_path / "priority_runs"
    robust_root = tmp_path / "robust_runs"
    plan_ab_bind_batches(benchmark_root, priority_protocol, spec_name="core_v1", runs_root=priority_root, complex_ids=["1VFB"])
    plan_ab_bind_batches(benchmark_root, robust_protocol, spec_name="core_v1", runs_root=robust_root, complex_ids=["1VFB"])

    priority_job = priority_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    robust_job = robust_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(priority_job, stage)
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(robust_job, stage)
    for job_dir, stderr_kt in ((robust_job, 0.05),):
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=stderr_kt)
        _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=stderr_kt)

    rescue_root = tmp_path / "rescues"
    rescue_payload = plan_ab_bind_rescues(
        priority_root,
        extra_plan_roots=[robust_root],
        complex_ids=["1VFB"],
        runs_root=rescue_root,
    )

    assert rescue_payload["rescued_job_count"] == 1
    assert rescue_payload["source_plan_roots"] == [str(priority_root.resolve()), str(robust_root.resolve())]

    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")
    assert len(rescue_rows) == 1
    assert rescue_rows[0]["source_plan_root"] == str(robust_root.resolve())

    rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    assert rescue_spec["protocol"]["repeats"] == 5
    assert rescue_spec["protocol"]["lambda_windows"] == 6
    assert rescue_spec["protocol"]["production_ps"] == 40
    assert rescue_spec["protocol"]["window_relax_em_steps"] == 1000
    assert rescue_spec["protocol"]["window_relax_md_ps"] == 0.2
    assert rescue_spec["protocol"]["nvt_ps"] == 200
    assert rescue_spec["protocol"]["npt_ps"] == 2000
    assert rescue_spec["protocol"]["equilibration_restraint_schedule"] == "staged_backbone_release"
    assert rescue_spec["protocol"]["equilibration_release_npt_ps"] == 1000


def test_plan_ab_bind_rescues_can_prefer_active_alternate_source_without_replacing_winner(
    tmp_path: Path,
    monkeypatch,
) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    priority_protocol = tmp_path / "protocol.priority.yml"
    robust_protocol = tmp_path / "protocol.robust.yml"
    write_yaml(
        priority_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )
    write_yaml(
        robust_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 6,
            "repeats": 4,
            "production_ps": 20,
            "temperature_k": 310.0,
        },
    )

    priority_root = tmp_path / "priority_runs"
    robust_root = tmp_path / "robust_runs"
    plan_ab_bind_batches(benchmark_root, priority_protocol, spec_name="core_v1", runs_root=priority_root, complex_ids=["1VFB"])
    plan_ab_bind_batches(benchmark_root, robust_protocol, spec_name="core_v1", runs_root=robust_root, complex_ids=["1VFB"])

    priority_job = priority_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    robust_job = robust_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
        _write_stage(priority_job, stage)
    _write_bar_outputs(priority_job / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(priority_job / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(priority_job / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(priority_job / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(robust_job, stage)
    _write_running_stage(
        robust_job,
        "equilibrate",
        commands=["bash equilibrate.sh"],
        artifacts=[str(robust_job / "artifacts" / "commands" / "equilibrate.sh")],
    )
    _write_equilibrate_repeat_completed(robust_job / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(robust_job / "legs" / "apo" / "rep01")

    monkeypatch.setattr(
        "abag_rbfe.reporting._active_process_lines",
        lambda: (f"fake {robust_job}/artifacts/commands/equilibrate.sh",),
    )

    rescue_root = tmp_path / "rescues"
    rescue_payload = plan_ab_bind_rescues(
        priority_root,
        extra_plan_roots=[robust_root],
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        prefer_active_alternate_source=True,
        require_active_alternate=True,
    )

    assert rescue_payload["rescued_job_count"] == 1
    assert rescue_payload["prefer_active_alternate_source"] is True
    assert rescue_payload["require_active_alternate"] is True

    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")
    assert len(rescue_rows) == 1
    assert rescue_rows[0]["source_plan_root"] == str(robust_root.resolve())
    assert rescue_rows[0]["source_winner_plan_root"] == str(priority_root.resolve())

    rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    assert rescue_spec["protocol"]["repeats"] == 5
    assert rescue_spec["protocol"]["lambda_windows"] == 6
    assert rescue_spec["protocol"]["production_ps"] == 40
    assert rescue_spec["protocol"]["window_relax_em_steps"] == 1000
    assert rescue_spec["protocol"]["window_relax_md_ps"] == 0.2
    assert rescue_spec["protocol"]["nvt_ps"] == 200
    assert rescue_spec["protocol"]["npt_ps"] == 2000
    assert rescue_spec["protocol"]["equilibration_restraint_schedule"] == "staged_backbone_release"
    assert rescue_spec["protocol"]["equilibration_release_npt_ps"] == 1000


def test_plan_ab_bind_rescues_prefers_matching_targeted_active_alternate_source(
    tmp_path: Path,
    monkeypatch,
) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    priority_protocol = tmp_path / "protocol.priority.yml"
    robust_protocol = tmp_path / "protocol.robust.yml"
    targeted_protocol = tmp_path / "protocol.targeted.yml"
    write_yaml(
        priority_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )
    write_yaml(
        robust_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 6,
            "repeats": 4,
            "production_ps": 20,
            "temperature_k": 310.0,
        },
    )
    write_yaml(
        targeted_protocol,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    priority_root = tmp_path / "priority_runs"
    robust_root = tmp_path / "robust_runs"
    targeted_root = tmp_path / "targeted_runs"
    plan_ab_bind_batches(benchmark_root, priority_protocol, spec_name="core_v1", runs_root=priority_root, complex_ids=["1VFB"])
    plan_ab_bind_batches(benchmark_root, robust_protocol, spec_name="core_v1", runs_root=robust_root, complex_ids=["1VFB"])
    plan_ab_bind_batches(benchmark_root, targeted_protocol, spec_name="core_v1", runs_root=targeted_root, complex_ids=["1VFB"])

    priority_job = priority_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    robust_job = robust_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    targeted_job = targeted_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"

    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
        _write_stage(priority_job, stage)
    _write_bar_outputs(priority_job / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(priority_job / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(priority_job / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(priority_job / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(robust_job, stage)
        _write_stage(targeted_job, stage)
    _write_running_stage(
        robust_job,
        "equilibrate",
        commands=["bash equilibrate.sh"],
        artifacts=[str(robust_job / "artifacts" / "commands" / "equilibrate.sh")],
    )
    _write_equilibrate_repeat_completed(robust_job / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(robust_job / "legs" / "apo" / "rep01")
    _write_running_stage(
        targeted_job,
        "equilibrate",
        commands=["bash equilibrate.sh"],
        artifacts=[str(targeted_job / "artifacts" / "commands" / "equilibrate.sh")],
    )
    _write_equilibrate_repeat_completed(targeted_job / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(targeted_job / "legs" / "apo" / "rep01")
    write_json(
        targeted_job / "config" / "rescue.json",
        {
            "mode": "targeted_primary_repeat_spread_leg",
            "target_legs": ["complex"],
            "inherit_source_legs": ["apo"],
            "source_plan_root": str(priority_root.resolve()),
        },
    )

    monkeypatch.setattr(
        "abag_rbfe.reporting._active_process_lines",
        lambda: (
            f"fake {robust_job}/artifacts/commands/equilibrate.sh",
            f"fake {targeted_job}/artifacts/commands/equilibrate.sh",
        ),
    )

    rescue_root = tmp_path / "rescues"
    rescue_payload = plan_ab_bind_rescues(
        priority_root,
        extra_plan_roots=[robust_root, targeted_root],
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        prefer_active_alternate_source=True,
        require_active_alternate=True,
        target_primary_repeat_spread_leg=True,
        allow_targeted_leg_count_deepening=True,
    )

    assert rescue_payload["rescued_job_count"] == 1
    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")
    assert len(rescue_rows) == 1
    assert rescue_rows[0]["source_plan_root"] == str(targeted_root.resolve())
    assert rescue_rows[0]["source_winner_plan_root"] == str(priority_root.resolve())
    assert rescue_rows[0]["targeted_primary_repeat_spread_leg"] == "True"
    assert rescue_rows[0]["target_legs"] == "complex"

    rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    assert rescue_spec["protocol"]["repeats"] == 3
    assert rescue_spec["protocol"]["lambda_windows"] == 2
    assert rescue_spec["protocol"]["production_ps"] == 20


def test_plan_ab_bind_rescues_can_force_repeat_increment(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 1,
            "production_ps": 1,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    job_dir = plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.0, stderr_kt=40.0)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=40.0)

    rescue_root = tmp_path / "rescue_runs"
    plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        repeat_increment=2,
        lambda_increment=6,
        production_scale=20.0,
        force_repeat_increment=True,
    )

    rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    rescue_payload = read_json(rescue_root / "reports" / "rescue_summary.json")

    assert rescue_spec["protocol"]["repeats"] == 3
    assert rescue_spec["protocol"]["lambda_windows"] == 8
    assert rescue_spec["protocol"]["production_ps"] == 20
    assert rescue_payload["force_repeat_increment"] is True


def test_plan_ab_bind_rescues_can_force_lambda_increment_for_repeat_spread(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    job_dir = plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    rescue_root = tmp_path / "rescue_runs"
    plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        lambda_increment=5,
        force_lambda_increment=True,
    )

    rescue_spec = read_json(
        rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a" / "job_spec.json"
    )
    rescue_payload = read_json(rescue_root / "reports" / "rescue_summary.json")

    assert rescue_spec["protocol"]["repeats"] == 3
    assert rescue_spec["protocol"]["lambda_windows"] == 7
    assert rescue_spec["protocol"]["production_ps"] == 20
    assert rescue_payload["force_lambda_increment"] is True


def test_rescue_protocol_payload_respects_zero_repeat_increment_and_unit_production_scale() -> None:
    rescued, adjustments = benchmark_module._rescue_protocol_payload(
        {
            "lambda_windows": 8,
            "repeats": 3,
            "production_ps": 40,
            "nvt_ps": 50,
            "npt_ps": 50,
        },
        reasons=["repeat_spread"],
        repeat_increment=0,
        lambda_increment=4,
        production_scale=1.0,
        window_relax_em_scale=1.0,
        window_relax_md_scale=1.0,
        nvt_scale=1.0,
        npt_scale=1.0,
        force_repeat_increment=False,
        force_lambda_increment=True,
        preserve_counts=False,
    )

    assert rescued["repeats"] == 3
    assert rescued["lambda_windows"] == 12
    assert rescued["production_ps"] == 40
    assert "repeats" not in adjustments
    assert adjustments["lambda_windows"] == 12
    assert "production_ps" not in adjustments


def test_plan_ab_bind_rescues_does_not_rewrite_started_existing_batches(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            }
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 1,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            }
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    job_dir = plan_root / "abbind_1vfb_core_v1" / "jobs" / "1vfb-antibody-h-y32a"
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)

    rescue_root = tmp_path / "rescue_runs"
    plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
    )

    rescue_job_dir = rescue_root / "abbind-rescue_1vfb-antibody-h-y32a" / "jobs" / "1vfb-antibody-h-y32a"
    initial_spec = read_json(rescue_job_dir / "job_spec.json")
    _write_stage(rescue_job_dir, "sample")

    plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        runs_root=rescue_root,
        lambda_increment=6,
        force_lambda_increment=True,
        production_scale=5.0,
    )

    updated_spec = read_json(rescue_job_dir / "job_spec.json")
    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")

    assert updated_spec == initial_spec
    assert rescue_rows[0]["source_job_id"] == "1vfb-antibody-h-y32a"
    assert rescue_rows[0]["rescued_lambda_windows"] == "2"


def test_plan_ab_bind_rescues_appends_batches_into_existing_rescue_root(tmp_path: Path) -> None:
    benchmark_root = tmp_path / "ab_bind"
    manifests_dir = benchmark_root / "manifests"
    curated_dir = benchmark_root / "curated"
    materialized_dir = benchmark_root / "materialized" / "1VFB"
    manifests_dir.mkdir(parents=True)
    curated_dir.mkdir(parents=True)
    materialized_dir.mkdir(parents=True)

    structure_path = materialized_dir / "1VFB.pdb"
    structure_path.write_text("HEADER    1VFB\nEND\n", encoding="utf-8")
    system_yml = materialized_dir / "system.yml"
    write_yaml(
        system_yml,
        {
            "system_name": "1vfb",
            "input_structure": str(structure_path),
            "structure_source": "experimental",
            "antibody_chains": ["H", "L"],
            "antigen_chains": ["C"],
            "notes": [],
        },
    )
    mutations_csv = materialized_dir / "core_v1_mutations.csv"
    write_csv_rows(
        mutations_csv,
        [
            {
                "mutation_group_id": "1vfb_0001",
                "chain_id": "H",
                "resseq": 32,
                "icode": "",
                "wt": "Y",
                "mut": "A",
                "entity_side": "antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "chain_id": "H",
                "resseq": 34,
                "icode": "",
                "wt": "V",
                "mut": "A",
                "entity_side": "antibody",
            },
        ],
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )
    write_csv_rows(
        manifests_dir / "ab_bind_rbfe_core_v1_inputs.csv",
        [
            {
                "complex_id": "1VFB",
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "structure_source": "experimental",
                "antibody_chains": "HL",
                "antigen_chains": "C",
                "mutation_group_count": 2,
            }
        ],
        [
            "complex_id",
            "system_yml",
            "mutations_csv",
            "structure_source",
            "antibody_chains",
            "antigen_chains",
            "mutation_group_count",
        ],
    )
    write_csv_rows(
        curated_dir / "ab_bind_rbfe_core_v1.csv",
        [
            {
                "mutation_group_id": "1vfb_0001",
                "ddg_kcal_mol": "0.5",
                "source_mutation": "H:Y32A",
                "mutation_tokens": "H:Y32A@antibody",
            },
            {
                "mutation_group_id": "1vfb_0002",
                "ddg_kcal_mol": "1.5",
                "source_mutation": "H:V34A",
                "mutation_tokens": "H:V34A@antibody",
            },
        ],
        ["mutation_group_id", "ddg_kcal_mol", "source_mutation", "mutation_tokens"],
    )

    protocol_path = tmp_path / "protocol.quick.yml"
    write_yaml(
        protocol_path,
        {
            "preset": "single_point",
            "lambda_windows": 2,
            "repeats": 2,
            "production_ps": 10,
            "temperature_k": 310.0,
        },
    )

    plan_root = tmp_path / "runs"
    plan_ab_bind_batches(benchmark_root, protocol_path, spec_name="core_v1", runs_root=plan_root, complex_ids=["1VFB"])
    jobs_dir = plan_root / "abbind_1vfb_core_v1" / "jobs"
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
            _write_stage(job_dir, stage)
        if job_dir.name == "1vfb-antibody-h-y32a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=0.10)
        elif job_dir.name == "1vfb-antibody-h-v34a":
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=2.0, stderr_kt=40.0)
            _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=0.0, stderr_kt=40.0)

    rescue_root = tmp_path / "rescue_runs"
    first_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        job_ids=["1vfb-antibody-h-y32a"],
        runs_root=rescue_root,
    )
    second_payload = plan_ab_bind_rescues(
        plan_root,
        complex_ids=["1VFB"],
        job_ids=["1vfb-antibody-h-v34a"],
        runs_root=rescue_root,
    )

    assert first_payload["rescued_job_count"] == 1
    assert second_payload["rescued_job_count"] == 2

    rescue_index = read_json(rescue_root / "plan_index.json")
    rescue_rows = read_csv_rows(rescue_root / "reports" / "rescue_candidates.csv")

    assert rescue_index["planned_batch_count"] == 2
    assert {item["batch_id"] for item in rescue_index["batches"]} == {
        "abbind-rescue_1vfb-antibody-h-y32a",
        "abbind-rescue_1vfb-antibody-h-v34a",
    }
    assert {row["source_job_id"] for row in rescue_rows} == {"1vfb-antibody-h-y32a", "1vfb-antibody-h-v34a"}
