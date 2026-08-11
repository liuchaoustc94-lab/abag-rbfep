from math import isclose
from pathlib import Path

from abag_rbfe.io_utils import ensure_dir, read_json, write_json, write_yaml
from abag_rbfe.reporting import build_qc_report, summarize_job, write_batch_summary, write_job_results, write_job_summary


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
    (repeat_dir / "bar" / "bar.xvg").write_text(f"0.500 {delta_kt:.2f} {stderr_kt:.2f}\n", encoding="utf-8")
    (repeat_dir / "bar" / "barint.xvg").write_text(f"0 0.00\n1 {delta_kt:.2f}\n", encoding="utf-8")
    histogram_lines = ["@ s0 legend \"mock\"", "@ s1 legend \"mock\"", "@ s2 legend \"mock\"", "@ s3 legend \"mock\"", "@ s4 legend \"mock\"", "@ s5 legend \"mock\""]
    for block in _histogram_blocks(overlap_mode):
        for x, y in block:
            histogram_lines.append(f"{x} {y}")
    (repeat_dir / "bar" / "histogram.xvg").write_text("\n".join(histogram_lines) + "\n", encoding="utf-8")


def _write_multilambda_bar_outputs(
    repeat_dir: Path,
    delta_kt: float,
    stderr_kt: float,
    *,
    lambda_values: tuple[float, ...] = (0.0, 0.3333, 0.6667, 1.0),
    overlap_mode: str = "good",
    signed_reverse: bool = False,
) -> None:
    ensure_dir(repeat_dir / "bar")
    for window_index in range(len(lambda_values)):
        ensure_dir(repeat_dir / f"lambda_{window_index:03d}")
        (repeat_dir / f"lambda_{window_index:03d}" / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")

    interval_value = delta_kt / max(len(lambda_values) - 1, 1)
    (repeat_dir / "bar" / "bar.xvg").write_text(
        "\n".join(
            f"{index + 0.5:.3f} {interval_value:.2f} {stderr_kt:.2f}"
            for index in range(len(lambda_values) - 1)
        )
        + "\n",
        encoding="utf-8",
    )
    (repeat_dir / "bar" / "barint.xvg").write_text(
        "\n".join(
            f"{index} {interval_value * index:.2f}"
            for index in range(len(lambda_values))
        )
        + "\n",
        encoding="utf-8",
    )

    derivative_series = [(1.0, 1.0), (2.0, 1.0), (3.0, 0.0), (4.0, 0.0)]
    if overlap_mode == "poor":
        forward_series = [(10.0, 1.0), (11.0, 1.0), (12.0, 1.0), (13.0, 1.0)]
        reverse_series = [(100.0, 1.0), (101.0, 1.0), (102.0, 1.0), (103.0, 1.0)]
    else:
        forward_series = [(10.0, 1.0), (11.0, 2.0), (12.0, 2.0), (13.0, 1.0)]
        reverse_series = [(11.0, 1.0), (12.0, 2.0), (13.0, 2.0), (14.0, 1.0)]
    empty_series = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)]

    histogram_lines = [
        '@    title "mock"',
        '@    xaxis  label "dH"',
        '@    yaxis  label "Samples"',
        "@TYPE xy",
    ]
    series_index = 0
    blocks: list[list[tuple[float, float]]] = []
    for current_index, current_lambda in enumerate(lambda_values):
        histogram_lines.append(f'@ s{series_index} legend "N(dH/dl | l={current_lambda})"')
        blocks.append(derivative_series)
        series_index += 1
        for target_index, target_lambda in enumerate(lambda_values):
            histogram_lines.append(
                f'@ s{series_index} legend "N(dH(l={target_lambda}) | l={current_lambda})"'
            )
            if target_index == current_index + 1:
                blocks.append(forward_series)
            elif target_index == current_index - 1:
                blocks.append([(-x, y) for x, y in reverse_series] if signed_reverse else reverse_series)
            else:
                blocks.append(empty_series)
            series_index += 1

    for block in blocks:
        histogram_lines.append("@")
        for x_value, y_value in block:
            histogram_lines.append(f"{x_value} {y_value}")
    (repeat_dir / "bar" / "histogram.xvg").write_text("\n".join(histogram_lines) + "\n", encoding="utf-8")


def test_write_job_results_extracts_ddg_and_qc(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "report")

    job_spec = {
        "job_id": "demo-job",
        "batch_id": "demo-batch",
        "mutation_group": {
            "mutation_group_id": "grp1",
            "mutation_count": 1,
            "entity_side": "antibody",
            "charge_conserving": True,
            "min_version": "v1",
            "sites": [
                {
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "F",
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
            "system_name": "demo-system",
            "structure_source": "experimental",
        },
    }
    write_json(job_dir / "job_spec.json", job_spec)
    write_yaml(job_dir / "config" / "system.yml", job_spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", job_spec["protocol"])

    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
        write_json(
            job_dir / "stages" / f"{stage}.json",
            {
                "stage": stage,
                "state": "completed",
                "message": f"{stage} ok",
                "commands": [],
                "artifacts": [],
            },
        )

    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=0.67)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=1.20)

    result_payload = write_job_results(job_dir)
    ddg_summary = result_payload["ddg_summary"]
    qc_report = result_payload["qc_report"]
    expected_ddg = (2.85 - 3.44) * 0.00198720425864083 * 310.0

    assert ddg_summary["ready"] is True
    assert isclose(ddg_summary["complex_delta_g_kcal_mol"], 2.85 * 0.00198720425864083 * 310.0, rel_tol=1e-6)
    assert isclose(ddg_summary["apo_delta_g_kcal_mol"], 3.44 * 0.00198720425864083 * 310.0, rel_tol=1e-6)
    assert isclose(ddg_summary["ddg_kcal_mol"], expected_ddg, rel_tol=1e-6)
    assert qc_report["status"] == "pass"
    assert Path(result_payload["paths"]["ddg_summary"]).exists()
    assert Path(result_payload["paths"]["qc_report"]).exists()

    summary = write_job_summary(job_dir)
    assert summary["results"]["ddg"]["ready"] is True
    assert summary["qc"]["status"] == "pass"
    assert Path(summary["result_files"]["ddg_summary_tsv"]).exists()


def _job_spec(job_id: str) -> dict:
    return {
        "job_id": job_id,
        "batch_id": "demo-batch",
        "mutation_group": {
            "mutation_group_id": f"{job_id}-grp",
            "mutation_count": 1,
            "entity_side": "antibody",
            "charge_conserving": True,
            "min_version": "v1",
            "sites": [
                {
                    "chain_id": "H",
                    "resseq": 32,
                    "icode": "",
                    "wt": "Y",
                    "mut": "F",
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
            "system_name": "demo-system",
            "structure_source": "experimental",
        },
    }


def _write_stage(
    job_dir: Path,
    stage: str,
    *,
    state: str = "completed",
    message: str | None = None,
    commands: list[str] | None = None,
    artifacts: list[str] | None = None,
) -> None:
    write_json(
        job_dir / "stages" / f"{stage}.json",
        {
            "stage": stage,
            "state": state,
            "message": message or f"{stage} {state}",
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


def _write_sample_window_pre_md(lambda_dir: Path) -> None:
    ensure_dir(lambda_dir)
    (lambda_dir / "pre_md.tpr").write_text("mock\n", encoding="utf-8")


def _write_equilibrate_repeat_completed(repeat_dir: Path) -> None:
    ensure_dir(repeat_dir / "equilibration")
    (repeat_dir / "system.top").write_text("mock\n", encoding="utf-8")
    (repeat_dir / "equilibration" / "npt.gro").write_text("mock\n", encoding="utf-8")


def _write_equilibrate_repeat_started(repeat_dir: Path) -> None:
    ensure_dir(repeat_dir / "setup")
    (repeat_dir / "system.top").write_text("mock\n", encoding="utf-8")
    (repeat_dir / "setup" / "genion.tpr").write_text("mock\n", encoding="utf-8")


def test_write_batch_summary_marks_unstarted_and_unevaluated_jobs(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    jobs_dir = batch_dir / "jobs"
    ensure_dir(jobs_dir)

    not_started_dir = jobs_dir / "job-not-started"
    ensure_dir(not_started_dir)
    write_json(not_started_dir / "job_spec.json", _job_spec("job-not-started"))

    partial_dir = jobs_dir / "job-partial"
    ensure_dir(partial_dir / "stages")
    write_json(partial_dir / "job_spec.json", _job_spec("job-partial"))
    _write_stage(partial_dir, "ingest")
    _write_stage(partial_dir, "prepare")

    sample_ready_dir = jobs_dir / "job-sample-ready"
    ensure_dir(sample_ready_dir / "stages")
    write_json(sample_ready_dir / "job_spec.json", _job_spec("job-sample-ready"))
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample"):
        _write_stage(sample_ready_dir, stage)

    completed_dir = jobs_dir / "job-completed"
    ensure_dir(completed_dir / "stages")
    write_json(completed_dir / "job_spec.json", _job_spec("job-completed"))
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(completed_dir, stage)
    _write_bar_outputs(completed_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=0.67)
    _write_bar_outputs(completed_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=1.20)

    payload = write_batch_summary(batch_dir)
    jobs = {job["job_id"]: job for job in payload["jobs"]}

    assert jobs["job-not-started"]["latest_stage_state"] == "not_started"
    assert jobs["job-not-started"]["latest_stage"] == ""
    assert jobs["job-not-started"]["qc_status"] == "not_started"
    assert jobs["job-not-started"]["benchmark_qc_qualified"] is False
    assert jobs["job-not-started"]["analyzable"] is False
    assert jobs["job-not-started"]["resumable"] is True
    assert jobs["job-not-started"]["diagnostic_family"] == "not_started"
    assert jobs["job-not-started"]["diagnostic_code"] == "not_started"

    assert jobs["job-partial"]["latest_stage"] == "prepare"
    assert jobs["job-partial"]["latest_stage_state"] == "completed"
    assert jobs["job-partial"]["qc_status"] == "not_evaluated"
    assert jobs["job-partial"]["benchmark_qc_qualified"] is False
    assert jobs["job-partial"]["analyzable"] is False
    assert jobs["job-partial"]["resumable"] is True
    assert jobs["job-partial"]["diagnostic_family"] == "pending"
    assert jobs["job-partial"]["diagnostic_code"] == "pending_mutate"

    assert jobs["job-sample-ready"]["latest_stage"] == "sample"
    assert jobs["job-sample-ready"]["latest_stage_state"] == "completed"
    assert jobs["job-sample-ready"]["analyzable"] is True
    assert jobs["job-sample-ready"]["resumable"] is False
    assert jobs["job-sample-ready"]["diagnostic_family"] == "pending"
    assert jobs["job-sample-ready"]["diagnostic_code"] == "pending_bar"

    assert jobs["job-completed"]["latest_stage"] == "qc"
    assert jobs["job-completed"]["latest_stage_state"] == "completed"
    assert jobs["job-completed"]["qc_status"] == "pass"
    assert jobs["job-completed"]["benchmark_qc_qualified"] is True
    assert jobs["job-completed"]["ddg_bar_stderr_kcal_mol"] is not None
    assert jobs["job-completed"]["analyzable"] is False
    assert jobs["job-completed"]["resumable"] is False
    assert jobs["job-completed"]["diagnostic_family"] == "completed"
    assert jobs["job-completed"]["diagnostic_code"] == "qc_pass"


def test_write_batch_summary_marks_failed_equilibrate_jobs_resumable(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-equilibrate-failed"
    ensure_dir(job_dir / "stages")
    write_json(job_dir / "job_spec.json", _job_spec("job-equilibrate-failed"))
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(job_dir, stage)
    _write_stage(job_dir, "equilibrate", state="failed", message="External command terminated by signal 15.")
    _write_equilibrate_repeat_started(job_dir / "legs" / "complex" / "rep01")

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "equilibrate"
    assert job_row["latest_stage_state"] == "failed"
    assert job_row["resumable"] is True
    assert job_row["equilibrate_started_repeats"] == 1
    assert job_row["equilibrate_completed_repeats"] == 0
    assert job_row["diagnostic_family"] == "equilibrate"
    assert job_row["diagnostic_code"] == "equilibrate_failed"


def test_write_batch_summary_marks_equilibrate_processed_gro_recovery_jobs_resumable(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-equilibrate-invalid-processed-gro"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    write_json(job_dir / "job_spec.json", _job_spec("job-equilibrate-invalid-processed-gro"))
    write_yaml(job_dir / "config" / "system.yml", _job_spec("job-equilibrate-invalid-processed-gro")["system"])
    write_yaml(job_dir / "config" / "protocol.yml", _job_spec("job-equilibrate-invalid-processed-gro")["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(job_dir, stage)
    _write_stage(
        job_dir,
        "equilibrate",
        state="blocked_input",
        message="Mutated coordinate file is invalid: /tmp/processed.gro (invalid_coordinate at line 42). Rerun mutate with --execute first.",
    )

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "equilibrate"
    assert job_row["latest_stage_state"] == "blocked_input"
    assert job_row["resumable"] is True
    assert job_row["diagnostic_family"] == "equilibrate"
    assert job_row["diagnostic_code"] == "equilibrate_invalid_processed_gro"


def test_write_batch_summary_classifies_prepare_and_mutate_input_failures(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    jobs_dir = batch_dir / "jobs"

    prepare_dir = jobs_dir / "job-prepare-blocked"
    ensure_dir(prepare_dir / "config")
    ensure_dir(prepare_dir / "stages")
    ensure_dir(prepare_dir / "artifacts")
    write_json(prepare_dir / "job_spec.json", _job_spec("job-prepare-blocked"))
    write_yaml(prepare_dir / "config" / "system.yml", _job_spec("job-prepare-blocked")["system"])
    write_yaml(prepare_dir / "config" / "protocol.yml", _job_spec("job-prepare-blocked")["protocol"])
    _write_stage(prepare_dir, "ingest")
    _write_stage(
        prepare_dir,
        "prepare",
        state="blocked_input",
        message="Prepared PDB contains backbone-incomplete standard residues: complex:H50 ASP missing backbone atoms C,O",
    )
    write_json(
        prepare_dir / "artifacts" / "prepare_qc.json",
        {
            "job_id": "job-prepare-blocked",
            "legs": {
                "complex": {
                    "blocking_incomplete_standard_residues": [
                        {"chain_id": "H", "resseq": 50, "icode": "", "resname": "ASP", "missing_backbone_atoms": ["C", "O"]}
                    ],
                    "blocking_intra_residue_heavy_atom_clashes": [],
                    "inter_residue_heavy_atom_clashes": [],
                },
                "apo": {
                    "blocking_incomplete_standard_residues": [],
                    "blocking_intra_residue_heavy_atom_clashes": [],
                    "inter_residue_heavy_atom_clashes": [],
                },
            },
        },
    )

    mutate_dir = jobs_dir / "job-mutate-invalid"
    ensure_dir(mutate_dir / "config")
    ensure_dir(mutate_dir / "stages")
    ensure_dir(mutate_dir / "artifacts")
    write_json(mutate_dir / "job_spec.json", _job_spec("job-mutate-invalid"))
    write_yaml(mutate_dir / "config" / "system.yml", _job_spec("job-mutate-invalid")["system"])
    write_yaml(mutate_dir / "config" / "protocol.yml", _job_spec("job-mutate-invalid")["protocol"])
    _write_stage(mutate_dir, "ingest")
    _write_stage(mutate_dir, "prepare")
    _write_stage(
        mutate_dir,
        "mutate",
        state="blocked_input",
        message="Mutated coordinate file is invalid: /tmp/processed.gro (invalid_coordinate at line 42).",
    )
    write_json(
        mutate_dir / "artifacts" / "mutate_qc.json",
        {
            "job_id": "job-mutate-invalid",
            "legs": {
                "complex": {
                    "inter_residue_heavy_atom_clashes": [],
                    "incomplete_standard_residues": [],
                    "processed_gro_qc": {"valid": False, "reason": "invalid_coordinate", "line_number": 42},
                },
                "apo": {
                    "inter_residue_heavy_atom_clashes": [],
                    "incomplete_standard_residues": [],
                    "processed_gro_qc": {"valid": True},
                },
            },
        },
    )

    payload = write_batch_summary(batch_dir)
    jobs = {job["job_id"]: job for job in payload["jobs"]}

    assert jobs["job-prepare-blocked"]["diagnostic_family"] == "input_structure"
    assert jobs["job-prepare-blocked"]["diagnostic_code"] == "input_backbone_incomplete"
    assert "backbone-incomplete" in jobs["job-prepare-blocked"]["diagnostic_detail"]
    assert jobs["job-mutate-invalid"]["diagnostic_family"] == "mutate_setup"
    assert jobs["job-mutate-invalid"]["diagnostic_code"] == "mutate_invalid_coordinate"
    assert "invalid_coordinate" in jobs["job-mutate-invalid"]["diagnostic_detail"]


def test_build_qc_report_marks_jobs_without_bar_data_as_not_evaluated(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    write_json(job_dir / "job_spec.json", _job_spec("job-empty"))
    write_yaml(job_dir / "config" / "system.yml", _job_spec("job-empty")["system"])
    write_yaml(job_dir / "config" / "protocol.yml", _job_spec("job-empty")["protocol"])

    qc_report = build_qc_report(job_dir)

    assert qc_report["status"] == "not_evaluated"
    assert qc_report["ddg_ready"] is False


def test_write_batch_summary_reclassifies_downstream_job_with_current_invalid_processed_gro(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    jobs_dir = batch_dir / "jobs"
    mutate_dir = jobs_dir / "job-mutated-output-drift"
    ensure_dir(mutate_dir / "artifacts")
    ensure_dir(mutate_dir / "config")
    ensure_dir(mutate_dir / "stages")
    ensure_dir(mutate_dir / "legs" / "complex" / "pmx")
    write_json(mutate_dir / "job_spec.json", _job_spec("job-mutated-output-drift"))
    write_yaml(mutate_dir / "config" / "system.yml", _job_spec("job-mutated-output-drift")["system"])
    write_yaml(mutate_dir / "config" / "protocol.yml", _job_spec("job-mutated-output-drift")["protocol"])
    _write_stage(mutate_dir, "ingest")
    _write_stage(mutate_dir, "prepare")
    _write_stage(mutate_dir, "mutate")
    _write_stage(mutate_dir, "build_legs")
    _write_stage(
        mutate_dir,
        "equilibrate",
        state="failed",
        message="Equilibrate stage failed after downstream setup.",
    )

    processed_gro = mutate_dir / "legs" / "complex" / "pmx" / "processed.gro"
    processed_gro.write_text(
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
    write_json(
        mutate_dir / "artifacts" / "mutate_qc.json",
        {
            "job_id": "job-mutated-output-drift",
            "legs": {
                "complex": {
                    "inter_residue_heavy_atom_clashes": [],
                    "incomplete_standard_residues": [],
                    "processed_gro": str(processed_gro),
                    "processed_gro_qc": {"valid": True, "reason": "ok"},
                }
            },
        },
    )

    payload = write_batch_summary(batch_dir)
    jobs = {job["job_id"]: job for job in payload["jobs"]}

    assert jobs["job-mutated-output-drift"]["diagnostic_family"] == "mutate_setup"
    assert jobs["job-mutated-output-drift"]["diagnostic_code"] == "mutate_processed_gro_isolated_residue_hydrogen"
    assert "ALA1 HA" in jobs["job-mutated-output-drift"]["diagnostic_detail"]
    assert jobs["job-mutated-output-drift"]["current_invalid_mutate_output"] is True
    assert jobs["job-mutated-output-drift"]["current_invalid_mutate_output_code"] == "mutate_processed_gro_isolated_residue_hydrogen"
    assert "ALA1 HA" in jobs["job-mutated-output-drift"]["current_invalid_mutate_output_detail"]


def test_write_batch_summary_marks_running_job_with_current_invalid_mutate_output(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-running-invalid-mutate-output"
    ensure_dir(job_dir / "artifacts")
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "legs" / "apo" / "pmx")
    write_json(job_dir / "job_spec.json", _job_spec("job-running-invalid-mutate-output"))
    write_yaml(job_dir / "config" / "system.yml", _job_spec("job-running-invalid-mutate-output")["system"])
    write_yaml(job_dir / "config" / "protocol.yml", _job_spec("job-running-invalid-mutate-output")["protocol"])
    _write_stage(job_dir, "ingest")
    _write_stage(job_dir, "prepare")
    _write_stage(job_dir, "mutate")
    _write_stage(job_dir, "build_legs")
    _write_stage(job_dir, "equilibrate")
    _write_stage(job_dir, "sample", state="running", message="Stage execution started. Script: sample.sh")

    processed_gro = job_dir / "legs" / "apo" / "pmx" / "processed.gro"
    processed_gro.write_text(
        "Mock GRO\n"
        "5\n"
        "    1ASP      N    1   0.000   0.000   0.000\n"
        "    1ASP     CA    2   0.145   0.000   0.000\n"
        "    1ASP      C    3   0.290   0.000   0.000\n"
        "    1ASP      O    4   0.410   0.000   0.000\n"
        "    1ASP     HA    5   5.000   5.000   5.000\n"
        "   6.00000   6.00000   6.00000\n",
        encoding="utf-8",
    )
    write_json(
        job_dir / "artifacts" / "mutate_qc.json",
        {
            "job_id": "job-running-invalid-mutate-output",
            "legs": {
                "apo": {
                    "inter_residue_heavy_atom_clashes": [],
                    "incomplete_standard_residues": [],
                    "processed_gro": str(processed_gro),
                    "processed_gro_qc": {"valid": True, "reason": "ok"},
                }
            },
        },
    )

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["diagnostic_family"] == "running"
    assert job_row["diagnostic_code"] == "running_sample"
    assert job_row["current_invalid_mutate_output"] is True
    assert job_row["current_invalid_mutate_output_code"] == "mutate_processed_gro_isolated_residue_hydrogen"
    assert "ASP1 HA" in job_row["current_invalid_mutate_output_detail"]


def test_build_qc_report_warns_on_large_bar_stderr(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "report")

    spec = _job_spec("job-noisy")
    spec["protocol"]["max_bar_stderr_kcal_mol"] = 10.0
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)

    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=40.0)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=60.0)

    qc_report = build_qc_report(job_dir)

    assert qc_report["status"] == "warning"
    assert any("BAR stderr" in warning for warning in qc_report["warnings"])


def test_build_qc_report_tracks_repeat_spread_legs_by_severity(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")

    spec = _job_spec("job-repeat-spread-legs")
    spec["protocol"]["repeats"] = 2
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)

    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=1.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=0.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "complex" / "rep02", delta_kt=4.0, stderr_kt=0.10)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep02", delta_kt=2.0, stderr_kt=0.10)

    qc_report = build_qc_report(job_dir)

    assert qc_report["status"] == "warning"
    assert qc_report["repeat_spread_legs"] == ["complex", "apo"]
    assert qc_report["primary_repeat_spread_leg"] == "complex"


def test_build_qc_report_warns_on_poor_overlap(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")

    spec = _job_spec("job-overlap")
    spec["protocol"]["overlap_threshold"] = 0.2
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)

    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=0.67, overlap_mode="poor")
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=1.20, overlap_mode="poor")

    qc_report = build_qc_report(job_dir)

    assert qc_report["status"] == "warning"
    assert any("overlap score" in warning for warning in qc_report["warnings"])


def test_build_qc_report_parses_multilambda_histogram_overlap(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")

    spec = _job_spec("job-multilambda")
    spec["protocol"]["lambda_windows"] = 4
    spec["protocol"]["overlap_threshold"] = 0.2
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)

    _write_multilambda_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=0.67)
    _write_multilambda_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=1.20)

    qc_report = build_qc_report(job_dir)

    assert qc_report["status"] == "pass"
    assert qc_report["overlap_assessment"]["legs"]["complex"]["overlap_score_min"] is not None
    assert qc_report["overlap_assessment"]["legs"]["apo"]["overlap_score_min"] is not None
    assert qc_report["overlap_assessment"]["legs"]["complex"]["overlap_score_min"] > 0.2
    assert qc_report["overlap_assessment"]["legs"]["apo"]["overlap_score_min"] > 0.2


def test_build_qc_report_reflects_signed_reverse_histogram_overlap(tmp_path: Path) -> None:
    job_dir = tmp_path / "job"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")

    spec = _job_spec("job-signed-multilambda")
    spec["protocol"]["lambda_windows"] = 4
    spec["protocol"]["overlap_threshold"] = 0.2
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc"):
        _write_stage(job_dir, stage)

    _write_multilambda_bar_outputs(
        job_dir / "legs" / "complex" / "rep01",
        delta_kt=2.85,
        stderr_kt=0.67,
        signed_reverse=True,
    )
    _write_multilambda_bar_outputs(
        job_dir / "legs" / "apo" / "rep01",
        delta_kt=3.44,
        stderr_kt=1.20,
        signed_reverse=True,
    )

    result_payload = write_job_results(job_dir)
    qc_report = result_payload["qc_report"]
    complex_repeat = result_payload["bar_summary"]["legs"]["complex"]["repeats"][0]
    apo_repeat = result_payload["bar_summary"]["legs"]["apo"]["repeats"][0]

    assert qc_report["status"] == "pass"
    assert complex_repeat["overlap_score"] is not None and complex_repeat["overlap_score"] > 0.2
    assert apo_repeat["overlap_score"] is not None and apo_repeat["overlap_score"] > 0.2
    assert complex_repeat["overlap_reverse_transform"] == "negate"
    assert apo_repeat["overlap_reverse_transform"] == "negate"


def test_write_batch_summary_refreshes_started_job_qc_files(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-noisy"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "report")

    spec = _job_spec("job-noisy")
    spec["protocol"]["max_bar_stderr_kcal_mol"] = 10.0
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
        _write_stage(job_dir, stage)

    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=40.0)
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=60.0)

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]
    qc_report_path = job_dir / "results" / "qc_report.json"

    assert job_row["qc_status"] == "warning"
    assert job_row["benchmark_qc_qualified"] is False
    assert qc_report_path.exists()
    assert read_json(qc_report_path)["status"] == "warning"
    assert job_row["diagnostic_family"] == "qc"
    assert job_row["diagnostic_code"] == "qc_bar_stderr"


def test_write_batch_summary_excludes_poor_overlap_from_benchmark_qc_qualified(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-overlap"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")

    spec = _job_spec("job-overlap")
    spec["protocol"]["overlap_threshold"] = 0.2
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"):
        _write_stage(job_dir, stage)

    _write_bar_outputs(job_dir / "legs" / "complex" / "rep01", delta_kt=2.85, stderr_kt=0.67, overlap_mode="poor")
    _write_bar_outputs(job_dir / "legs" / "apo" / "rep01", delta_kt=3.44, stderr_kt=1.20, overlap_mode="poor")

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["qc_status"] == "warning"
    assert job_row["benchmark_qc_qualified"] is False
    assert job_row["diagnostic_family"] == "qc"
    assert job_row["diagnostic_code"] == "qc_low_overlap"


def test_write_batch_summary_does_not_persist_transient_results_for_running_jobs(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-running"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "report")
    ensure_dir(job_dir / "results")

    spec = _job_spec("job-running")
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(job_dir, stage)
    _write_stage(job_dir, "sample", state="running", message="Stage execution started.")

    # Simulate stale transient artifacts left behind by an earlier partial refresh.
    write_json(job_dir / "results" / "ddg_summary.json", {"ready": False, "job_id": "job-running"})
    write_json(job_dir / "results" / "qc_report.json", {"status": "fail"})
    write_json(job_dir / "report" / "summary.json", {"job": {"job_id": "job-running"}})

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "sample"
    assert job_row["latest_stage_state"] == "running"
    assert job_row["qc_status"] == "not_evaluated"
    assert job_row["benchmark_qc_qualified"] is False
    assert job_row["analyzable"] is False
    assert job_row["resumable"] is False
    assert job_row["diagnostic_family"] == "running"
    assert job_row["diagnostic_code"] == "running_sample"
    assert not (job_dir / "results" / "ddg_summary.json").exists()
    assert not (job_dir / "results" / "qc_report.json").exists()
    assert not (job_dir / "report" / "summary.json").exists()


def test_summarize_job_tracks_running_sample_progress(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-running-sample"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    spec = _job_spec("job-running-sample")
    spec["protocol"]["lambda_windows"] = 2
    spec["protocol"]["repeats"] = 1
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(job_dir, stage)
    _write_stage(job_dir, "sample", state="running", message="Stage execution started.")

    _write_equilibrate_repeat_completed(job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    summary = summarize_job(job_dir)

    assert summary["progress"]["equilibrate"]["completed_repeats"] == 2
    assert summary["progress"]["equilibrate"]["total_repeats"] == 2
    assert summary["progress"]["sample"]["completed_windows"] == 1
    assert summary["progress"]["sample"]["started_windows"] == 2
    assert summary["progress"]["sample"]["total_windows"] == 4
    assert summary["progress"]["sample"]["active_leg"] == "complex"
    assert summary["progress"]["sample"]["active_repeat_id"] == "rep01"
    assert summary["progress"]["sample"]["active_lambda_id"] == "lambda_001"
    assert summary["progress"]["sample"]["active_lambda_index"] == 1
    assert summary["progress"]["sample"]["active_phase"] == "md"
    assert summary["progress"]["sample"]["active_window"] == "complex/rep01/lambda_001"

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]
    assert job_row["sample_completed_windows"] == 1
    assert job_row["sample_started_windows"] == 2
    assert job_row["sample_total_windows"] == 4
    assert job_row["sample_active_leg"] == "complex"
    assert job_row["sample_active_repeat_id"] == "rep01"
    assert job_row["sample_active_lambda_id"] == "lambda_001"
    assert job_row["sample_active_lambda_index"] == 1
    assert job_row["sample_active_phase"] == "md"
    assert job_row["sample_active_window"] == "complex/rep01/lambda_001"
    assert job_row["equilibrate_completed_repeats"] == 2
    assert job_row["equilibrate_total_repeats"] == 2


def test_summarize_job_tracks_pre_md_sample_phase(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-running-sample-pre-md"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    spec = _job_spec("job-running-sample-pre-md")
    spec["protocol"]["lambda_windows"] = 1
    spec["protocol"]["repeats"] = 1
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(job_dir, stage)
    _write_stage(job_dir, "sample", state="running", message="Stage execution started.")

    _write_equilibrate_repeat_completed(job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_pre_md(job_dir / "legs" / "apo" / "rep01" / "lambda_000")

    summary = summarize_job(job_dir)

    assert summary["progress"]["sample"]["started_windows"] == 1
    assert summary["progress"]["sample"]["completed_windows"] == 0
    assert summary["progress"]["sample"]["active_leg"] == "apo"
    assert summary["progress"]["sample"]["active_repeat_id"] == "rep01"
    assert summary["progress"]["sample"]["active_lambda_id"] == "lambda_000"
    assert summary["progress"]["sample"]["active_lambda_index"] == 0
    assert summary["progress"]["sample"]["active_phase"] == "pre_md"
    assert summary["progress"]["sample"]["active_window"] == "apo/rep01/lambda_000"


def test_write_batch_summary_tracks_running_equilibrate_progress(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-running-equilibrate"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    spec = _job_spec("job-running-equilibrate")
    spec["protocol"]["repeats"] = 2
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs"):
        _write_stage(job_dir, stage)
    _write_stage(job_dir, "equilibrate", state="running", message="Stage execution started.")

    _write_equilibrate_repeat_completed(job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_started(job_dir / "legs" / "complex" / "rep02")

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "equilibrate"
    assert job_row["latest_stage_state"] == "running"
    assert job_row["resumable"] is False
    assert job_row["equilibrate_completed_repeats"] == 1
    assert job_row["equilibrate_started_repeats"] == 2
    assert job_row["equilibrate_total_repeats"] == 4
    assert job_row["sample_total_windows"] == 8


def test_write_batch_summary_marks_stale_running_jobs_resumable(tmp_path: Path, monkeypatch) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-stale-sample"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    spec = _job_spec("job-stale-sample")
    spec["protocol"]["lambda_windows"] = 2
    spec["protocol"]["repeats"] = 1
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    for stage in ("ingest", "prepare", "mutate", "build_legs", "equilibrate"):
        _write_stage(job_dir, stage)
    _write_stage(
        job_dir,
        "sample",
        state="running",
        message="Stage execution started.",
        commands=["bash sample.sh"],
        artifacts=[str(job_dir / "artifacts" / "commands" / "sample.sh")],
    )

    _write_equilibrate_repeat_completed(job_dir / "legs" / "complex" / "rep01")
    _write_equilibrate_repeat_completed(job_dir / "legs" / "apo" / "rep01")
    _write_sample_window_completed(job_dir / "legs" / "complex" / "rep01" / "lambda_000")
    _write_sample_window_started(job_dir / "legs" / "complex" / "rep01" / "lambda_001")

    monkeypatch.setattr("abag_rbfe.reporting._active_process_lines", lambda: ())

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "sample"
    assert job_row["latest_stage_state"] == "stale_running"
    assert job_row["resumable"] is True
    assert job_row["sample_completed_windows"] == 1
    assert job_row["sample_started_windows"] == 2
    assert job_row["diagnostic_family"] == "running"
    assert job_row["diagnostic_code"] == "stale_running_sample"


def test_write_batch_summary_marks_repairable_blocked_mutate_qc_jobs_resumable(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-blocked-mutate-repairable"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "artifacts")
    spec = _job_spec("job-blocked-mutate-repairable")
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    _write_stage(job_dir, "ingest")
    _write_stage(job_dir, "prepare")
    _write_stage(job_dir, "mutate", state="blocked_input", message="repairable sidechain clash")
    write_json(
        job_dir / "artifacts" / "mutate_qc.json",
        {
            "job_id": job_dir.name,
            "legs": {
                "complex": {
                    "inter_residue_heavy_atom_clashes": [
                        {
                            "chain_id": "H",
                            "resseq": 50,
                            "icode": "",
                            "partner_chain_id": "H",
                            "partner_resseq": 51,
                            "partner_icode": "",
                            "clashes": [
                                {"atom_a": "CG", "atom_b": "CD", "distance_angstrom": 0.92},
                                {"atom_a": "NE", "atom_b": "CZ", "distance_angstrom": 1.01},
                            ],
                        }
                    ]
                },
                "apo": {
                    "inter_residue_heavy_atom_clashes": []
                },
            },
        },
    )

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "mutate"
    assert job_row["latest_stage_state"] == "blocked_input"
    assert job_row["resumable"] is True
    assert job_row["diagnostic_family"] == "mutate_setup"
    assert job_row["diagnostic_code"] == "mutate_inter_residue_clash"


def test_write_batch_summary_marks_backbone_sidechain_blocked_mutate_qc_jobs_resumable(tmp_path: Path) -> None:
    batch_dir = tmp_path / "batch"
    job_dir = batch_dir / "jobs" / "job-blocked-mutate-backbone-sidechain"
    ensure_dir(job_dir / "config")
    ensure_dir(job_dir / "stages")
    ensure_dir(job_dir / "artifacts")
    spec = _job_spec("job-blocked-mutate-backbone-sidechain")
    write_json(job_dir / "job_spec.json", spec)
    write_yaml(job_dir / "config" / "system.yml", spec["system"])
    write_yaml(job_dir / "config" / "protocol.yml", spec["protocol"])
    _write_stage(job_dir, "ingest")
    _write_stage(job_dir, "prepare")
    _write_stage(job_dir, "mutate", state="blocked_input", message="repairable backbone-sidechain clash")
    write_json(
        job_dir / "artifacts" / "mutate_qc.json",
        {
            "job_id": job_dir.name,
            "legs": {
                "complex": {
                    "inter_residue_heavy_atom_clashes": [
                        {
                            "chain_id": "H",
                            "resseq": 50,
                            "icode": "",
                            "partner_chain_id": "H",
                            "partner_resseq": 51,
                            "partner_icode": "",
                            "clashes": [
                                {"atom_a": "O", "atom_b": "CD", "distance_angstrom": 0.92},
                            ],
                        }
                    ]
                },
                "apo": {
                    "inter_residue_heavy_atom_clashes": []
                },
            },
        },
    )

    payload = write_batch_summary(batch_dir)
    job_row = payload["jobs"][0]

    assert job_row["latest_stage"] == "mutate"
    assert job_row["latest_stage_state"] == "blocked_input"
    assert job_row["resumable"] is True
    assert job_row["diagnostic_family"] == "mutate_setup"
    assert job_row["diagnostic_code"] == "mutate_inter_residue_clash"


def test_flag_censored_experimental_values_marks_saturated_max_duplicates() -> None:
    from abag_rbfe.reporting import flag_censored_experimental_values

    pairs = [
        {"complex_id": "1BJ1", "job_id": "j1", "experimental_ddg_kcal_mol": 3.69},
        {"complex_id": "1BJ1", "job_id": "j2", "experimental_ddg_kcal_mol": 3.69},
        {"complex_id": "1BJ1", "job_id": "j3", "experimental_ddg_kcal_mol": 3.69},
        {"complex_id": "1BJ1", "job_id": "j4", "experimental_ddg_kcal_mol": 0.82},
        {"complex_id": "1BJ1", "job_id": "j5", "experimental_ddg_kcal_mol": 0.0},
        {"complex_id": "1MLC", "job_id": "j6", "experimental_ddg_kcal_mol": -1.25},
        {"complex_id": "1MLC", "job_id": "j7", "experimental_ddg_kcal_mol": 0.53},
    ]
    flagged = flag_censored_experimental_values(pairs)
    by_job = {r["job_id"]: r for r in flagged}
    assert by_job["j1"]["censored_experimental"] is True
    assert by_job["j2"]["censored_experimental"] is True
    assert by_job["j3"]["censored_experimental"] is True
    assert "detection_limit" in by_job["j1"]["censored_reason"]
    assert by_job["j4"]["censored_experimental"] is False  # duplicated but not max
    assert by_job["j5"]["censored_experimental"] is False  # zero is not a saturation marker
    assert by_job["j6"]["censored_experimental"] is False  # 1MLC max is 0.53, unique
    assert by_job["j7"]["censored_experimental"] is False  # unique max is not censored
