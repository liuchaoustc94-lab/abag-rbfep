from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path


SCRIPT_PATH = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/report_project_completion.py")


def _load_script_module():
    spec = importlib.util.spec_from_file_location("report_project_completion", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def _write_real_case(
    root: Path,
    *,
    run_dir: str,
    job_dir: str,
    mutation_signature: str,
    mutation_count: int,
    protocol_preset: str,
    entity_side: str,
    ddg_kcal_mol: float,
    ddg_bar_stderr_kcal_mol: float,
    qc_status: str,
    warnings: list[str] | None = None,
) -> None:
    base = root / "runs" / "real_cases" / run_dir / "jobs" / job_dir
    _write_json(
        base / "results" / "ddg_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "ready": True,
            "mutation_signature": mutation_signature,
            "mutation_count": mutation_count,
            "protocol_preset": protocol_preset,
            "entity_side": entity_side,
            "ddg_kcal_mol": ddg_kcal_mol,
            "ddg_bar_stderr_kcal_mol": ddg_bar_stderr_kcal_mol,
        },
    )
    _write_json(
        base / "results" / "qc_report.json",
        {
            "status": qc_status,
            "warnings": warnings or [],
        },
    )
    _write_json(
        base / "job_spec.json",
        {
            "job_id": job_dir,
            "mutation_group": {
                "entity_side": entity_side,
                "sites": [{} for _ in range(mutation_count)],
            },
            "protocol": {
                "preset": protocol_preset,
            },
        },
    )


def _write_patellike_3hfm_summary(
    root: Path,
    *,
    status: str,
    paired_job_count: int,
    incomplete_job_count: int,
    generated_at: str = "2026-06-12T03:00:00Z",
    message: str = "",
    charge_conserving_paired: int = 0,
    charge_conserving_incomplete: int = 0,
    charge_changing_paired: int = 0,
    charge_changing_incomplete: int = 0,
) -> None:
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "patel_2021_3hfm"
        / "patel_2021_3hfm_reference"
        / "reports"
        / "patel_2021_3hfm_summary.json",
        {
            "generated_at": generated_at,
            "status": status,
            "paired_job_count": paired_job_count,
            "incomplete_job_count": incomplete_job_count,
            "message": message,
            "charge_class_summary": {
                "charge_conserving": {
                    "paired_job_count": charge_conserving_paired,
                    "incomplete_job_count": charge_conserving_incomplete,
                },
                "charge_changing": {
                    "paired_job_count": charge_changing_paired,
                    "incomplete_job_count": charge_changing_incomplete,
                },
            },
        },
    )


def test_default_output_paths_land_under_docs_and_runs() -> None:
    module = _load_script_module()
    root = Path("/tmp/abag-rbfep")

    assert module.default_summary_output_path(root) == root / "docs" / "project_completion_status.md"
    assert module.default_json_output_path(root) == root / "runs" / "benchmarks" / "project_completion_summary.json"
    assert (
        module.default_patellike_3hfm_summary_path(root)
        == root
        / "runs"
        / "benchmarks"
        / "patel_2021_3hfm"
        / "patel_2021_3hfm_reference"
        / "reports"
        / "patel_2021_3hfm_summary.json"
    )


def test_main_marks_project_incomplete_when_validation_passes_but_real_runs_are_still_active(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-12T02:15:22Z",
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_excluded_complex_ids": ["1MLC", "1CZ8", "1BJ1"],
            "accepted_calibrated_pearson_r": 0.6718905648684699,
            "accepted_calibrated_passed": True,
            "calibrated_pearson_r": 0.14533357582231624,
            "predict_pair_count": 61,
            "accepted_calibrated_pair_count": 32,
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-12T02:15:22Z",
            "selected_job_count": 80,
            "ddg_ready_count": 80,
            "paired_job_count": 80,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-12T02:11:02Z",
            "selected_job_count": 80,
            "ddg_ready_count": 59,
            "paired_job_count": 59,
            "running_sample_job_count": 19,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_robust_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-11T06:37:44Z",
            "selected_job_count": 80,
            "ddg_ready_count": 0,
            "paired_job_count": 0,
            "running_sample_job_count": 13,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_rescues" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-11T06:40:44Z",
            "selected_job_count": 15,
            "ddg_ready_count": 4,
            "paired_job_count": 4,
            "running_sample_job_count": 8,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_patellike_3hfm_summary(
        root,
        status="ok",
        paired_job_count=8,
        incomplete_job_count=0,
        charge_conserving_paired=3,
        charge_changing_paired=5,
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-3.59,
        ddg_bar_stderr_kcal_mol=5.66,
        qc_status="warning",
        warnings=["complex overlap low"],
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="double_point",
        entity_side="antibody",
        ddg_kcal_mol=-2.12,
        ddg_bar_stderr_kcal_mol=7.79,
        qc_status="warning",
        warnings=["apo overlap low"],
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=18.17,
        ddg_bar_stderr_kcal_mol=55.64,
        qc_status="warning",
        warnings=["ddG BAR stderr high"],
    )

    def fake_check_output(args, text=True):
        if args == ["ps", "-eo", "args"]:
            return "\n".join(
                [
                    "ARGS",
                    f"/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx mdrun -s {root}/runs/benchmarks/abbind_core_v1_validation_priority_plan/abbind_a/jobs/job-a/legs/complex/rep01/lambda_000/topol.tpr",
                    f"{root}/.venv/bin/python3.11 {root}/.venv/bin/abag-rbfe resume job-a --batch-dir {root}/runs/benchmarks/abbind_core_v1_validation_priority_plan/abbind_a --execute",
                    f"{root}/.venv/bin/python3.11 {root}/.venv/bin/abag-rbfe resume patel-a --batch-dir {root}/runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference --execute",
                ]
            )
        if args == ["ps", "-eo", "pid,etimes,args"]:
            return "\n".join(
                [
                    "PID ELAPSED COMMAND",
                    f"1234 300 /mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx mdrun -s {root}/runs/benchmarks/abbind_core_v1_validation_priority_plan/abbind_a/jobs/job-a/legs/complex/rep01/lambda_000/topol.tpr",
                ]
            )
        raise AssertionError(args)

    monkeypatch.setattr(module.subprocess, "check_output", fake_check_output)
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
            "--snapshot-date",
            "2026-06-12",
        ],
    )

    assert module.main() == 0

    payload = json.loads(json_output.read_text(encoding="utf-8"))
    assert payload["completion_gates"]["accepted_independent_validation_passed"] is True
    assert payload["completion_gates"]["required_real_case_checkpoints_completed"] is True
    assert payload["completion_gates"]["same_side_double_point_checkpoint_completed"] is True
    assert payload["completion_gates"]["external_3hfm_reference_completed"] is True
    assert payload["completion_gates"]["no_live_core_processes"] is False
    assert payload["completion_gates"]["no_live_reference_processes"] is False
    assert payload["completion_gates"]["tracked_plan_roots_drained"] is False
    assert payload["completion_gates"]["project_complete"] is False
    assert payload["live_processes"]["core_mdrun_process_count"] == 1
    assert payload["live_processes"]["core_resume_process_count"] == 1
    assert payload["live_processes"]["core_active_resume_job_count"] == 1
    assert payload["live_processes"]["core_active_mdrun_job_count"] == 1
    assert payload["live_processes"]["orphaned_core_resume_job_count"] == 0
    assert payload["live_processes"]["reference_resume_process_count"] == 1
    assert payload["tracked_plan_roots"]["priority"]["drained"] is False
    assert payload["tracked_plan_roots"]["robust"]["pending_selected_job_count"] == 80
    assert payload["real_case_checkpoints"][1]["execution_checkpoint_passed"] is True
    assert payload["real_case_checkpoints"][1]["protocol_preset"] == "double_point"
    text = summary_output.read_text(encoding="utf-8")
    assert "Snapshot date: June 12, 2026." in text
    assert "accepted excluded complexes: `1MLC, 1CZ8, 1BJ1`" in text
    assert "required real-case checkpoints completed: `True`" in text
    assert "same-side double-point checkpoint completed: `True`" in text
    assert "external reference checkpoint passed: `True`" in text
    assert "no live reference processes remain: `False`" in text
    assert "active core AB-Bind `gmx mdrun` processes: `1`" in text
    assert "unique core active job ids (`resume` / `mdrun`): `1` / `1`" in text
    assert "orphaned core `resume` job ids: `none`" in text
    assert "project complete: `False`" in text


def test_main_uses_validation_target_summary_when_calibrated_summary_lacks_accepted_metric(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-25T07:18:34Z",
            "status": "insufficient_fit_pairs",
            "fit_pair_count": 0,
        },
    )
    _write_json(
        root / "docs" / "validation_target_summary" / "validation_target_summary.json",
        {
            "generated_at": "2026-06-25T06:40:58.525159+00:00",
            "selected_model": "side_linear",
            "calibrated_pearson_r": 0.20005445478956882,
            "accepted_calibrated_pearson_r": 0.6073390390160122,
            "accepted_calibrated_excluded_complex_ids": ["1MLC", "1CZ8", "1BJ1"],
            "targets": [
                {
                    "complex_id": "GOOD1",
                    "pair_count": 20,
                    "calibrated_metrics": {"excluded_from_target_filtered_metrics": False},
                },
                {
                    "complex_id": "GOOD2",
                    "pair_count": 22,
                    "calibrated_metrics": {"excluded_from_target_filtered_metrics": False},
                },
                {
                    "complex_id": "BAD1",
                    "pair_count": 18,
                    "calibrated_metrics": {"excluded_from_target_filtered_metrics": True},
                },
            ],
        },
    )
    _write_patellike_3hfm_summary(
        root,
        status="insufficient_pairs",
        paired_job_count=0,
        incomplete_job_count=8,
    )
    for label in (
        "abbind_core_v1_validation_priority_plan",
        "abbind_core_v1_validation_robust_plan",
        "abbind_core_v1_validation_priority_rescues",
        "abbind_core_v1_validation_targeted_repeat_spread_rescues",
        "abbind_core_v1_validation_targeted_lambda_rescues",
        "abbind_core_v1_validation_sampling_qc_rescues",
        "abbind_core_v1_validation_deep_rescues",
        "abbind_core_v1_validation_ultra_rescues",
    ):
        _write_json(
            root / "runs" / "benchmarks" / label / "reports" / "plan_summary.json",
            {
                "generated_at": "2026-06-25T07:00:00Z",
                "selected_job_count": 0,
                "ddg_ready_count": 0,
                "paired_job_count": 0,
                "running_sample_job_count": 0,
                "running_equilibrate_job_count": 0,
            },
        )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-1.0,
        ddg_bar_stderr_kcal_mol=1.0,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="same_side_double_point",
        entity_side="antibody",
        ddg_kcal_mol=-1.0,
        ddg_bar_stderr_kcal_mol=1.0,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=1.0,
        ddg_bar_stderr_kcal_mol=1.0,
        qc_status="warning",
    )

    monkeypatch.setattr(module, "_ps_lines", lambda: [])
    monkeypatch.setattr(module, "_ps_status_rows", lambda: [])
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
            "--snapshot-date",
            "2026-06-25",
        ],
    )

    assert module.main() == 0
    summary = json.loads(json_output.read_text(encoding="utf-8"))
    assert summary["independent_validation"]["accepted_pearson_r"] == 0.6073390390160122
    assert summary["independent_validation"]["accepted_passed"] is True
    assert summary["independent_validation"]["accepted_pair_count"] == 42
    assert summary["independent_validation"]["summary_path"].endswith(
        "docs/validation_target_summary/validation_target_summary.json"
    )
    text = summary_output.read_text(encoding="utf-8")
    assert "accepted gate passed: `True`" in text
    assert "accepted independent validation passed: `True`" in text
    assert "Patel-like external 3HFM reference regression is incomplete" in text


def test_main_reports_active_untracked_core_benchmark_roots(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_pearson_r": 0.71,
            "accepted_calibrated_passed": True,
            "calibrated_pearson_r": 0.2,
            "predict_pair_count": 61,
            "accepted_calibrated_pair_count": 53,
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_job_count": 10,
            "ddg_ready_count": 10,
            "paired_job_count": 10,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_patellike_3hfm_summary(
        root,
        status="ok",
        paired_job_count=8,
        incomplete_job_count=0,
        charge_conserving_paired=3,
        charge_changing_paired=5,
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-3.59,
        ddg_bar_stderr_kcal_mol=5.66,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="double_point",
        entity_side="antibody",
        ddg_kcal_mol=-2.12,
        ddg_bar_stderr_kcal_mol=7.79,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=18.17,
        ddg_bar_stderr_kcal_mol=55.64,
        qc_status="warning",
    )

    untracked_root = "abbind_core_v1_validation_sampling_qc_rescues_verify_preeq_20260611_1"

    def fake_check_output(args, text=True):
        if args == ["ps", "-eo", "args"]:
            return "\n".join(
                [
                    "ARGS",
                    f"/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx mdrun -s {root}/runs/benchmarks/{untracked_root}/abbind_a/jobs/job-a/legs/complex/rep01/lambda_000/topol.tpr",
                    f"{root}/.venv/bin/python3.11 {root}/.venv/bin/abag-rbfe resume job-a --batch-dir {root}/runs/benchmarks/{untracked_root}/abbind_a --execute",
                ]
            )
        if args == ["ps", "-eo", "pid,etimes,args"]:
            return "\n".join(
                [
                    "PID ELAPSED COMMAND",
                    f"1234 300 /mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx mdrun -s {root}/runs/benchmarks/{untracked_root}/abbind_a/jobs/job-a/legs/complex/rep01/lambda_000/topol.tpr",
                ]
            )
        raise AssertionError(args)

    monkeypatch.setattr(module.subprocess, "check_output", fake_check_output)
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
        ],
    )

    assert module.main() == 0

    payload = json.loads(json_output.read_text(encoding="utf-8"))
    assert payload["completion_gates"]["accepted_independent_validation_passed"] is True
    assert payload["completion_gates"]["required_real_case_checkpoints_completed"] is True
    assert payload["completion_gates"]["external_3hfm_reference_completed"] is True
    assert payload["completion_gates"]["tracked_plan_roots_drained"] is True
    assert payload["completion_gates"]["project_complete"] is False
    assert payload["live_processes"]["active_untracked_core_benchmark_roots"] == [untracked_root]
    assert any("outside the tracked completion set" in item for item in payload["completion_blockers"])
    text = summary_output.read_text(encoding="utf-8")
    assert f"active untracked core benchmark roots: `{untracked_root}`" in text
    assert f"active core benchmark roots outside the tracked completion set: {untracked_root}" in text


def test_main_auto_tracks_active_core_benchmark_roots_with_plan_summary(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_pearson_r": 0.71,
            "accepted_calibrated_passed": True,
            "calibrated_pearson_r": 0.2,
            "predict_pair_count": 61,
            "accepted_calibrated_pair_count": 53,
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_job_count": 10,
            "ddg_ready_count": 10,
            "paired_job_count": 10,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_patellike_3hfm_summary(
        root,
        status="ok",
        paired_job_count=8,
        incomplete_job_count=0,
        charge_conserving_paired=3,
        charge_changing_paired=5,
    )
    extra_root = "abbind_core_v1_validation_sampling_qc_rescues_verify_preeq_20260611_1"
    _write_json(
        root / "runs" / "benchmarks" / extra_root / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-12T03:05:00Z",
            "selected_job_count": 3,
            "ddg_ready_count": 0,
            "paired_job_count": 0,
            "running_sample_job_count": 3,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-3.59,
        ddg_bar_stderr_kcal_mol=5.66,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="double_point",
        entity_side="antibody",
        ddg_kcal_mol=-2.12,
        ddg_bar_stderr_kcal_mol=7.79,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=18.17,
        ddg_bar_stderr_kcal_mol=55.64,
        qc_status="warning",
    )

    def fake_check_output(args, text=True):
        if args == ["ps", "-eo", "args"]:
            return "\n".join(
                [
                    "ARGS",
                    f"/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx mdrun -s {root}/runs/benchmarks/{extra_root}/abbind_a/jobs/job-a/legs/complex/rep01/lambda_000/topol.tpr",
                    f"{root}/.venv/bin/python3.11 {root}/.venv/bin/abag-rbfe resume job-a --batch-dir {root}/runs/benchmarks/{extra_root}/abbind_a --execute",
                ]
            )
        if args == ["ps", "-eo", "pid,etimes,args"]:
            return "\n".join(
                [
                    "PID ELAPSED COMMAND",
                    f"1234 300 /mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx mdrun -s {root}/runs/benchmarks/{extra_root}/abbind_a/jobs/job-a/legs/complex/rep01/lambda_000/topol.tpr",
                ]
            )
        raise AssertionError(args)

    monkeypatch.setattr(module.subprocess, "check_output", fake_check_output)
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
        ],
    )

    assert module.main() == 0

    payload = json.loads(json_output.read_text(encoding="utf-8"))
    assert payload["live_processes"]["active_untracked_core_benchmark_roots"] == []
    assert "sampling_qc_rescues_verify_preeq_20260611_1" in payload["tracked_plan_roots"]
    assert (
        payload["tracked_plan_roots"]["sampling_qc_rescues_verify_preeq_20260611_1"]["running_sample_job_count"]
        == 3
    )
    assert not any("outside the tracked completion set" in item for item in payload["completion_blockers"])
    text = summary_output.read_text(encoding="utf-8")
    assert "active untracked core benchmark roots: `none`" in text
    assert "- sampling_qc_rescues_verify_preeq_20260611_1:" in text


def test_main_marks_project_complete_when_validation_passes_and_no_live_work_remains(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_pearson_r": 0.71,
            "accepted_calibrated_passed": True,
            "calibrated_pearson_r": 0.2,
            "predict_pair_count": 61,
            "accepted_calibrated_pair_count": 53,
        },
    )
    for rel in (
        "abbind_core_v1_quick_plan",
        "abbind_core_v1_validation_priority_plan",
        "abbind_core_v1_validation_robust_plan",
    ):
        _write_json(
            root / "runs" / "benchmarks" / rel / "reports" / "plan_summary.json",
            {
                "generated_at": "2026-06-12T03:00:00Z",
                "selected_job_count": 10,
                "ddg_ready_count": 10,
                "paired_job_count": 10,
                "running_sample_job_count": 0,
                "running_equilibrate_job_count": 0,
            },
        )
    _write_patellike_3hfm_summary(
        root,
        status="ok",
        paired_job_count=8,
        incomplete_job_count=0,
        charge_conserving_paired=3,
        charge_changing_paired=5,
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-3.59,
        ddg_bar_stderr_kcal_mol=5.66,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="double_point",
        entity_side="antibody",
        ddg_kcal_mol=-2.12,
        ddg_bar_stderr_kcal_mol=7.79,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=18.17,
        ddg_bar_stderr_kcal_mol=55.64,
        qc_status="warning",
    )

    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: "ARGS\n")
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
        ],
    )

    assert module.main() == 0

    payload = json.loads(json_output.read_text(encoding="utf-8"))
    assert payload["completion_gates"]["accepted_independent_validation_passed"] is True
    assert payload["completion_gates"]["required_real_case_checkpoints_completed"] is True
    assert payload["completion_gates"]["same_side_double_point_checkpoint_completed"] is True
    assert payload["completion_gates"]["external_3hfm_reference_completed"] is True
    assert payload["completion_gates"]["no_live_core_processes"] is True
    assert payload["completion_gates"]["no_live_reference_processes"] is True
    assert payload["completion_gates"]["tracked_plan_roots_drained"] is True
    assert payload["completion_gates"]["project_complete"] is True
    text = summary_output.read_text(encoding="utf-8")
    assert "project complete: `True`" in text
    assert "Current blockers:" in text
    assert "- none" in text


def test_main_reports_stale_core_mdrun_processes_as_blocker(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_pearson_r": 0.71,
            "accepted_calibrated_passed": True,
            "calibrated_pearson_r": 0.2,
            "predict_pair_count": 61,
            "accepted_calibrated_pair_count": 53,
        },
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-3.59,
        ddg_bar_stderr_kcal_mol=5.66,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="double_point",
        entity_side="antibody",
        ddg_kcal_mol=-2.12,
        ddg_bar_stderr_kcal_mol=7.79,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=18.17,
        ddg_bar_stderr_kcal_mol=55.64,
        qc_status="warning",
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_job_count": 10,
            "ddg_ready_count": 10,
            "paired_job_count": 10,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
        },
    )
    _write_patellike_3hfm_summary(
        root,
        status="ok",
        paired_job_count=8,
        incomplete_job_count=0,
        charge_conserving_paired=3,
        charge_changing_paired=5,
    )

    deffnm = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "abbind_a"
        / "jobs"
        / "job-a"
        / "legs"
        / "complex"
        / "rep01"
        / "lambda_000"
        / "md"
    )
    deffnm.parent.mkdir(parents=True, exist_ok=True)
    log_path = deffnm.with_suffix(".log")
    log_path.write_text("old log\n", encoding="utf-8")
    stale_epoch = 1000.0
    os.utime(log_path, (stale_epoch, stale_epoch))

    def fake_check_output(args, text=True):
        if args == ["ps", "-eo", "args"]:
            return "\n".join(
                [
                    "ARGS",
                    f"/path/to/gmx mdrun -deffnm {deffnm}",
                ]
            )
        if args == ["ps", "-eo", "pid,etimes,args"]:
            return "\n".join(
                [
                    "PID ELAPSED COMMAND",
                    f"1234 1800 /path/to/gmx mdrun -deffnm {deffnm}",
                ]
            )
        raise AssertionError(args)

    monkeypatch.setattr(module.subprocess, "check_output", fake_check_output)
    monkeypatch.setattr(module.time, "time", lambda: stale_epoch + 5000.0)
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
        ],
    )

    assert module.main() == 0

    payload = json.loads(json_output.read_text(encoding="utf-8"))
    assert payload["live_processes"]["stale_core_mdrun_process_count"] == 1
    assert payload["completion_gates"]["project_complete"] is False
    assert any("appear stale" in item for item in payload["completion_blockers"])
    text = summary_output.read_text(encoding="utf-8")
    assert "stale core `gmx mdrun` processes" in text
    assert "Stale core processes:" in text
    assert "`job-a` pid=1234" in text


def test_main_reports_incomplete_patellike_3hfm_reference_as_blocker(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "project_completion_status.md"
    json_output = root / "runs" / "benchmarks" / "project_completion_summary.json"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-12T03:00:00Z",
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_pearson_r": 0.71,
            "accepted_calibrated_passed": True,
            "calibrated_pearson_r": 0.2,
            "predict_pair_count": 61,
            "accepted_calibrated_pair_count": 53,
        },
    )
    for rel in (
        "abbind_core_v1_quick_plan",
        "abbind_core_v1_validation_priority_plan",
        "abbind_core_v1_validation_robust_plan",
    ):
        _write_json(
            root / "runs" / "benchmarks" / rel / "reports" / "plan_summary.json",
            {
                "generated_at": "2026-06-12T03:00:00Z",
                "selected_job_count": 10,
                "ddg_ready_count": 10,
                "paired_job_count": 10,
                "running_sample_job_count": 0,
                "running_equilibrate_job_count": 0,
            },
        )
    _write_patellike_3hfm_summary(
        root,
        status="insufficient_pairs",
        paired_job_count=0,
        incomplete_job_count=8,
        message="No completed Patel 2021 3HFM jobs are ready for external regression comparison yet.",
        charge_conserving_incomplete=3,
        charge_changing_incomplete=5,
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_quick",
        job_dir="1vfb-antibody-b-y32f",
        mutation_signature="B:Y32F@antibody",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antibody",
        ddg_kcal_mol=-3.59,
        ddg_bar_stderr_kcal_mol=5.66,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="1vfb_y32f_v34i_quick",
        job_dir="1vfb-antibody-b-y32f--b-v34i",
        mutation_signature="B:Y32F@antibody__B:V34I@antibody",
        mutation_count=2,
        protocol_preset="double_point",
        entity_side="antibody",
        ddg_kcal_mol=-2.12,
        ddg_bar_stderr_kcal_mol=7.79,
        qc_status="warning",
    )
    _write_real_case(
        root,
        run_dir="4dn4_v47i_quick",
        job_dir="4dn4-antigen-m-v47i",
        mutation_signature="M:V47I@antigen",
        mutation_count=1,
        protocol_preset="single_point",
        entity_side="antigen",
        ddg_kcal_mol=18.17,
        ddg_bar_stderr_kcal_mol=55.64,
        qc_status="warning",
    )

    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: "ARGS\n")
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_project_completion.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--json-output",
            str(json_output),
        ],
    )

    assert module.main() == 0

    payload = json.loads(json_output.read_text(encoding="utf-8"))
    assert payload["completion_gates"]["accepted_independent_validation_passed"] is True
    assert payload["completion_gates"]["required_real_case_checkpoints_completed"] is True
    assert payload["completion_gates"]["external_3hfm_reference_completed"] is False
    assert payload["completion_gates"]["no_live_reference_processes"] is True
    assert payload["completion_gates"]["tracked_plan_roots_drained"] is True
    assert payload["completion_gates"]["project_complete"] is False
    assert any("Patel-like external 3HFM reference regression is incomplete" in item for item in payload["completion_blockers"])
    text = summary_output.read_text(encoding="utf-8")
    assert "status: `insufficient_pairs`" in text
    assert "external reference checkpoint passed: `False`" in text
    assert "no live reference processes remain: `True`" in text
