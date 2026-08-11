#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


VALIDATION_TARGET_R = 0.6
STALE_PROGRESS_THRESHOLD_SECONDS = 900
MDRUN_RE = re.compile(r"(?:^|\s)\S*gmx(?:_mpi)?\s+mdrun(?:\s|$)")
BENCHMARK_ROOT_RE = re.compile(r"/runs/benchmarks/([^/\s]+)/")
WATCHER_RE = re.compile(
    r"(?:run_validation_.*watcher\.sh|run_validation_watchers\.sh|watch_validation_priority\.py|"
    r"run_calibration_.*watcher\.sh|run_calibration_watchers\.sh|watch_calibration)"
)
CLI_JOB_RE = re.compile(r"\babag-rbfe\b\s+(?:run|resume|analyze)\s+([^\s]+)")
REAL_CASE_CHECKPOINTS = (
    {
        "slug": "1vfb_single",
        "label": "1VFB single-point quick validation",
        "ddg_relpath": (
            "runs/real_cases/1vfb_y32f_quick/jobs/1vfb-antibody-b-y32f/results/ddg_summary.json"
        ),
        "qc_relpath": (
            "runs/real_cases/1vfb_y32f_quick/jobs/1vfb-antibody-b-y32f/results/qc_report.json"
        ),
        "required": True,
        "expected_mutation_count": 1,
        "checkpoint_kind": "single_point",
    },
    {
        "slug": "1vfb_double",
        "label": "1VFB same-side double-point quick validation",
        "ddg_relpath": (
            "runs/real_cases/1vfb_y32f_v34i_quick/jobs/1vfb-antibody-b-y32f--b-v34i/results/ddg_summary.json"
        ),
        "qc_relpath": (
            "runs/real_cases/1vfb_y32f_v34i_quick/jobs/1vfb-antibody-b-y32f--b-v34i/results/qc_report.json"
        ),
        "required": True,
        "expected_mutation_count": 2,
        "checkpoint_kind": "same_side_double_point",
    },
    {
        "slug": "4dn4_single",
        "label": "4DN4 larger real-case quick validation",
        "ddg_relpath": (
            "runs/real_cases/4dn4_v47i_quick/jobs/4dn4-antigen-m-v47i/results/ddg_summary.json"
        ),
        "qc_relpath": (
            "runs/real_cases/4dn4_v47i_quick/jobs/4dn4-antigen-m-v47i/results/qc_report.json"
        ),
        "required": True,
        "expected_mutation_count": 1,
        "checkpoint_kind": "single_point_stress",
    },
)


def default_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_summary_output_path(root: Path) -> Path:
    return root / "docs" / "project_completion_status.md"


def default_json_output_path(root: Path) -> Path:
    return root / "runs" / "benchmarks" / "project_completion_summary.json"


def default_patellike_3hfm_summary_path(root: Path) -> Path:
    return (
        root
        / "runs"
        / "benchmarks"
        / "patel_2021_3hfm"
        / "patel_2021_3hfm_reference"
        / "reports"
        / "patel_2021_3hfm_summary.json"
    )


def default_validation_target_summary_path(root: Path) -> Path:
    return root / "docs" / "validation_target_summary" / "validation_target_summary.json"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="report_project_completion.py")
    parser.add_argument("--root", default=str(default_root()))
    parser.add_argument("--summary-output", default="")
    parser.add_argument("--json-output", default="")
    parser.add_argument("--snapshot-date", default="")
    return parser.parse_args(argv)


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _relpath(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _format_snapshot_date(raw_value: str) -> str:
    value = raw_value.strip()
    if not value:
        date_value = datetime.now(timezone.utc).date()
    else:
        date_value = datetime.fromisoformat(value).date()
    return date_value.strftime("%B %d, %Y").replace(" 0", " ")


def _format_target_list(value: Any) -> str:
    if isinstance(value, list):
        normalized = [str(item).strip() for item in value if str(item).strip()]
        if normalized:
            return ", ".join(normalized)
    return "none"


def _coerce_float(value: Any) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return None


def _coerce_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    if isinstance(value, str):
        try:
            return int(value.strip())
        except ValueError:
            return None
    return None


def _ps_lines() -> list[str]:
    output = subprocess.check_output(["ps", "-eo", "args"], text=True)
    return [line.strip() for line in output.splitlines()[1:] if line.strip()]


def _ps_status_rows() -> list[tuple[int, int, str]]:
    output = subprocess.check_output(["ps", "-eo", "pid,etimes,args"], text=True)
    rows: list[tuple[int, int, str]] = []
    for line in output.splitlines()[1:]:
        match = re.match(r"\s*(\d+)\s+(\d+)\s+(.*)$", line)
        if not match:
            continue
        rows.append((int(match.group(1)), int(match.group(2)), match.group(3)))
    return rows


def _file_age_seconds(path: str | None) -> float | None:
    if not path:
        return None
    try:
        return max(time.time() - os.path.getmtime(path), 0.0)
    except (FileNotFoundError, OSError):
        return None


def _parse_active_mdrun_status(pid: int, elapsed_seconds: int, command: str) -> dict[str, Any]:
    deffnm = ""
    dhdl_path = ""
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    for index, token in enumerate(tokens):
        if token == "-deffnm" and index + 1 < len(tokens):
            deffnm = tokens[index + 1]
        elif token == "-dhdl" and index + 1 < len(tokens):
            dhdl_path = tokens[index + 1]
    log_path = f"{deffnm}.log" if deffnm else ""
    progress_ages = [
        age
        for age in (
            _file_age_seconds(log_path),
            _file_age_seconds(dhdl_path),
        )
        if age is not None
    ]
    job_id_match = re.search(r"/jobs/([^/]+)/", deffnm or command)
    root_name_match = BENCHMARK_ROOT_RE.search(command)
    progress_age = min(progress_ages) if progress_ages else None
    return {
        "pid": pid,
        "elapsed_seconds": elapsed_seconds,
        "elapsed_minutes": round(elapsed_seconds / 60.0, 1),
        "job_id": job_id_match.group(1) if job_id_match else "",
        "root_name": root_name_match.group(1) if root_name_match else "",
        "deffnm": deffnm,
        "log_path": log_path,
        "dhdl_path": dhdl_path,
        "log_age_seconds": _file_age_seconds(log_path),
        "dhdl_age_seconds": _file_age_seconds(dhdl_path),
        "progress_age_seconds": progress_age,
    }


def _live_mdrun_statuses(root: Path) -> dict[str, list[dict[str, Any]]]:
    root_str = str(root.resolve())
    core_prefix = f"{root_str}/runs/benchmarks/abbind_core_v1_"
    patel_prefix = f"{root_str}/runs/benchmarks/patel_2021_3hfm/"
    core_statuses: list[dict[str, Any]] = []
    reference_statuses: list[dict[str, Any]] = []
    for pid, elapsed_seconds, command in _ps_status_rows():
        if not MDRUN_RE.search(command):
            continue
        status = _parse_active_mdrun_status(pid, elapsed_seconds, command)
        if core_prefix in command:
            core_statuses.append(status)
        elif patel_prefix in command:
            reference_statuses.append(status)
    return {
        "core": core_statuses,
        "reference": reference_statuses,
    }


def _stale_statuses(
    statuses: list[dict[str, Any]],
    *,
    warn_stale_mdrun_seconds: int = STALE_PROGRESS_THRESHOLD_SECONDS,
) -> list[dict[str, Any]]:
    if warn_stale_mdrun_seconds <= 0:
        return []
    stale: list[dict[str, Any]] = []
    for status in statuses:
        progress_age = status.get("progress_age_seconds")
        if progress_age is None:
            if int(status.get("elapsed_seconds") or 0) >= warn_stale_mdrun_seconds:
                stale.append(status)
            continue
        if float(progress_age) >= warn_stale_mdrun_seconds:
            stale.append(status)
    stale.sort(key=lambda item: (-float(item.get("progress_age_seconds") or 0.0), item.get("job_id") or ""))
    return stale


def _live_process_summary(root: Path) -> dict[str, Any]:
    root_str = str(root.resolve())
    core_prefix = f"{root_str}/runs/benchmarks/abbind_core_v1_"
    patel_prefix = f"{root_str}/runs/benchmarks/patel_2021_3hfm/"
    core_benchmark_roots: set[str] = set()
    reference_benchmark_roots: set[str] = set()
    counts = {
        "core_mdrun_process_count": 0,
        "core_resume_process_count": 0,
        "reference_mdrun_process_count": 0,
        "reference_resume_process_count": 0,
        "watcher_process_count": 0,
    }
    core_resume_job_ids: set[str] = set()
    reference_resume_job_ids: set[str] = set()
    for line in _ps_lines():
        benchmark_root_match = BENCHMARK_ROOT_RE.search(line)
        benchmark_root_name = benchmark_root_match.group(1) if benchmark_root_match else ""
        if WATCHER_RE.search(line) and root_str in line:
            counts["watcher_process_count"] += 1
        if MDRUN_RE.search(line):
            if core_prefix in line:
                counts["core_mdrun_process_count"] += 1
                if benchmark_root_name:
                    core_benchmark_roots.add(benchmark_root_name)
            elif patel_prefix in line:
                counts["reference_mdrun_process_count"] += 1
                if benchmark_root_name:
                    reference_benchmark_roots.add(benchmark_root_name)
        if "abag-rbfe resume" in line:
            cli_match = CLI_JOB_RE.search(line)
            resume_job_id = cli_match.group(1) if cli_match else ""
            if core_prefix in line:
                counts["core_resume_process_count"] += 1
                if resume_job_id:
                    core_resume_job_ids.add(resume_job_id)
                if benchmark_root_name:
                    core_benchmark_roots.add(benchmark_root_name)
            elif patel_prefix in line:
                counts["reference_resume_process_count"] += 1
                if resume_job_id:
                    reference_resume_job_ids.add(resume_job_id)
                if benchmark_root_name:
                    reference_benchmark_roots.add(benchmark_root_name)
    mdrun_statuses = _live_mdrun_statuses(root)
    stale_core_statuses = _stale_statuses(mdrun_statuses["core"])
    stale_reference_statuses = _stale_statuses(mdrun_statuses["reference"])
    core_mdrun_job_ids = {status["job_id"] for status in mdrun_statuses["core"] if status.get("job_id")}
    reference_mdrun_job_ids = {
        status["job_id"] for status in mdrun_statuses["reference"] if status.get("job_id")
    }
    orphaned_core_resume_job_ids = sorted(core_resume_job_ids - core_mdrun_job_ids)
    orphaned_reference_resume_job_ids = sorted(reference_resume_job_ids - reference_mdrun_job_ids)
    return {
        **counts,
        "stale_threshold_seconds": STALE_PROGRESS_THRESHOLD_SECONDS,
        "stale_core_mdrun_process_count": len(stale_core_statuses),
        "stale_reference_mdrun_process_count": len(stale_reference_statuses),
        "stale_core_mdrun_statuses": stale_core_statuses[:10],
        "stale_reference_mdrun_statuses": stale_reference_statuses[:10],
        "core_active_resume_job_count": len(core_resume_job_ids),
        "core_active_mdrun_job_count": len(core_mdrun_job_ids),
        "orphaned_core_resume_job_count": len(orphaned_core_resume_job_ids),
        "orphaned_core_resume_job_ids": orphaned_core_resume_job_ids[:20],
        "reference_active_resume_job_count": len(reference_resume_job_ids),
        "reference_active_mdrun_job_count": len(reference_mdrun_job_ids),
        "orphaned_reference_resume_job_count": len(orphaned_reference_resume_job_ids),
        "orphaned_reference_resume_job_ids": orphaned_reference_resume_job_ids[:20],
        "active_core_benchmark_roots": sorted(core_benchmark_roots),
        "active_reference_benchmark_roots": sorted(reference_benchmark_roots),
        "core_process_count": counts["core_mdrun_process_count"] + counts["core_resume_process_count"],
        "reference_process_count": counts["reference_mdrun_process_count"] + counts["reference_resume_process_count"],
    }


def _dynamic_tracked_root_label(root_name: str) -> str:
    prefix = "abbind_core_v1_validation_"
    if root_name.startswith(prefix):
        return root_name[len(prefix) :]
    return root_name


def _tracked_plan_summary_paths(
    root: Path,
    *,
    active_core_root_names: list[str] | None = None,
) -> dict[str, Path]:
    benchmarks_root = root / "runs" / "benchmarks"
    tracked_paths = {
        "priority": benchmarks_root / "abbind_core_v1_validation_priority_plan" / "reports" / "plan_summary.json",
        "robust": benchmarks_root / "abbind_core_v1_validation_robust_plan" / "reports" / "plan_summary.json",
        "rescue": benchmarks_root / "abbind_core_v1_validation_priority_rescues" / "reports" / "plan_summary.json",
        "targeted_repeat": benchmarks_root
        / "abbind_core_v1_validation_targeted_repeat_spread_rescues"
        / "reports"
        / "plan_summary.json",
        "targeted_lambda": benchmarks_root
        / "abbind_core_v1_validation_targeted_lambda_rescues"
        / "reports"
        / "plan_summary.json",
        "sampling_qc": benchmarks_root
        / "abbind_core_v1_validation_sampling_qc_rescues"
        / "reports"
        / "plan_summary.json",
        "deep": benchmarks_root / "abbind_core_v1_validation_deep_rescues" / "reports" / "plan_summary.json",
        "ultra": benchmarks_root / "abbind_core_v1_validation_ultra_rescues" / "reports" / "plan_summary.json",
    }
    seen_paths = {path.resolve() for path in tracked_paths.values()}
    for root_name in active_core_root_names or []:
        summary_path = benchmarks_root / root_name / "reports" / "plan_summary.json"
        resolved_summary_path = summary_path.resolve()
        if resolved_summary_path in seen_paths or not summary_path.is_file():
            continue
        label = _dynamic_tracked_root_label(root_name)
        if label in tracked_paths:
            label = root_name
        tracked_paths[label] = summary_path
        seen_paths.add(resolved_summary_path)
    return tracked_paths


def _tracked_plan_state(summary: dict[str, Any] | None) -> dict[str, Any] | None:
    if summary is None:
        return None
    selected_job_count = _coerce_int(summary.get("selected_job_count")) or 0
    ddg_ready_count = _coerce_int(summary.get("ddg_ready_count")) or 0
    paired_job_count = _coerce_int(summary.get("paired_job_count")) or 0
    running_sample_job_count = _coerce_int(summary.get("running_sample_job_count")) or 0
    running_equilibrate_job_count = _coerce_int(summary.get("running_equilibrate_job_count")) or 0
    pending_selected_job_count = max(selected_job_count - ddg_ready_count, 0)
    drained = (
        pending_selected_job_count == 0
        and running_sample_job_count == 0
        and running_equilibrate_job_count == 0
    )
    return {
        "generated_at": str(summary.get("generated_at", "")),
        "selected_job_count": selected_job_count,
        "ddg_ready_count": ddg_ready_count,
        "paired_job_count": paired_job_count,
        "running_sample_job_count": running_sample_job_count,
        "running_equilibrate_job_count": running_equilibrate_job_count,
        "pending_selected_job_count": pending_selected_job_count,
        "drained": drained,
    }


def _real_case_checkpoint_state(root: Path, spec: dict[str, Any]) -> dict[str, Any]:
    ddg_path = root / str(spec["ddg_relpath"])
    qc_path = root / str(spec["qc_relpath"])
    ddg_payload = _load_json(ddg_path) or {}
    qc_payload = _load_json(qc_path) or {}
    job_spec_path = ddg_path.parents[1] / "job_spec.json"
    job_spec_payload = _load_json(job_spec_path) or {}
    mutation_group = job_spec_payload.get("mutation_group", {}) if isinstance(job_spec_payload, dict) else {}
    mutation_count = (
        _coerce_int(ddg_payload.get("mutation_count"))
        or _coerce_int(mutation_group.get("mutation_count"))
        or len(mutation_group.get("sites", []))
    )
    protocol_preset = str(
        ddg_payload.get("protocol_preset")
        or ((job_spec_payload.get("protocol") or {}).get("preset") if isinstance(job_spec_payload.get("protocol"), dict) else "")
        or ""
    )
    ready = bool(ddg_payload.get("ready")) if ddg_payload else False
    qc_status = str(qc_payload.get("status") or "") if qc_payload else ""
    execution_checkpoint_passed = bool(ddg_payload and qc_payload and ready and qc_status in {"pass", "warning"})
    return {
        "slug": str(spec["slug"]),
        "label": str(spec["label"]),
        "required": bool(spec.get("required")),
        "checkpoint_kind": str(spec.get("checkpoint_kind", "")),
        "ddg_path": str(ddg_path),
        "qc_path": str(qc_path),
        "job_spec_path": str(job_spec_path),
        "available": bool(ddg_payload),
        "ready": ready,
        "generated_at": str(ddg_payload.get("generated_at", "")),
        "mutation_signature": str(ddg_payload.get("mutation_signature", "")),
        "mutation_count": mutation_count,
        "protocol_preset": protocol_preset,
        "entity_side": str(ddg_payload.get("entity_side", "") or mutation_group.get("entity_side", "")),
        "ddg_kcal_mol": _coerce_float(ddg_payload.get("ddg_kcal_mol")),
        "ddg_bar_stderr_kcal_mol": _coerce_float(ddg_payload.get("ddg_bar_stderr_kcal_mol")),
        "qc_status": qc_status,
        "qc_warning_count": len(qc_payload.get("warnings", [])) if isinstance(qc_payload.get("warnings"), list) else 0,
        "qc_warnings": list(qc_payload.get("warnings", []))[:5] if isinstance(qc_payload.get("warnings"), list) else [],
        "execution_checkpoint_passed": execution_checkpoint_passed,
        "matches_expected_mutation_count": mutation_count == _coerce_int(spec.get("expected_mutation_count")),
    }


def _external_reference_state(root: Path) -> dict[str, Any]:
    summary_path = default_patellike_3hfm_summary_path(root)
    payload = _load_json(summary_path) or {}
    charge_class_summary = payload.get("charge_class_summary", {})
    if not isinstance(charge_class_summary, dict):
        charge_class_summary = {}
    status = str(payload.get("status", "") or "")
    paired_job_count = _coerce_int(payload.get("paired_job_count")) or 0
    incomplete_job_count = _coerce_int(payload.get("incomplete_job_count")) or 0
    completed = bool(payload) and status == "ok" and paired_job_count > 0 and incomplete_job_count == 0
    return {
        "summary_path": str(summary_path),
        "available": bool(payload),
        "generated_at": str(payload.get("generated_at", "")),
        "status": status,
        "paired_job_count": paired_job_count,
        "incomplete_job_count": incomplete_job_count,
        "message": str(payload.get("message", "") or ""),
        "charge_class_summary": charge_class_summary,
        "completed": completed,
    }


def _validation_target_summary_state(root: Path) -> dict[str, Any]:
    summary_path = default_validation_target_summary_path(root)
    payload = _load_json(summary_path) or {}
    targets = payload.get("targets", [])
    if not isinstance(targets, list):
        targets = []
    predict_pair_count = 0
    accepted_pair_count = 0
    for item in targets:
        if not isinstance(item, dict):
            continue
        pair_count = _coerce_int(item.get("pair_count")) or _coerce_int(
            (item.get("calibrated_metrics") or {}).get("paired_job_count")
        ) or 0
        predict_pair_count += pair_count
        calibrated_metrics = item.get("calibrated_metrics") if isinstance(item.get("calibrated_metrics"), dict) else {}
        if not bool(calibrated_metrics.get("excluded_from_target_filtered_metrics")):
            accepted_pair_count += pair_count
    accepted_pearson_r = _coerce_float(payload.get("accepted_calibrated_pearson_r"))
    return {
        "summary_path": str(summary_path),
        "available": bool(payload),
        "generated_at": str(payload.get("generated_at", "")),
        "selected_model": str(payload.get("selected_model", "")),
        "accepted_view": "target_filtered" if accepted_pearson_r is not None else "",
        "accepted_pearson_r": accepted_pearson_r,
        "accepted_excluded_complex_ids": payload.get("accepted_calibrated_excluded_complex_ids", []),
        "accepted_passed": accepted_pearson_r is not None and accepted_pearson_r >= VALIDATION_TARGET_R,
        "full_calibrated_pearson_r": _coerce_float(payload.get("calibrated_pearson_r")),
        "predict_pair_count": predict_pair_count or None,
        "accepted_pair_count": accepted_pair_count or None,
    }


def build_completion_summary(*, root: Path) -> dict[str, Any]:
    calibrated_summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_quick_plan"
        / "reports"
        / "calibrated_validation_summary.json"
    )
    calibrated = _load_json(calibrated_summary_path) or {}
    validation_target_summary = _validation_target_summary_state(root)
    live_processes = _live_process_summary(root)
    real_case_checkpoints = [_real_case_checkpoint_state(root, spec) for spec in REAL_CASE_CHECKPOINTS]
    external_reference = _external_reference_state(root)
    tracked_plan_summary_paths = _tracked_plan_summary_paths(
        root,
        active_core_root_names=live_processes.get("active_core_benchmark_roots", []),
    )
    tracked_plan_summaries = {
        name: _tracked_plan_state(_load_json(path))
        for name, path in tracked_plan_summary_paths.items()
    }
    tracked_core_root_names = {
        path.parent.parent.name for path in tracked_plan_summary_paths.values()
    }
    active_untracked_core_benchmark_roots = sorted(
        name
        for name in live_processes.get("active_core_benchmark_roots", [])
        if name not in tracked_core_root_names
    )
    live_processes = dict(live_processes)
    live_processes["active_untracked_core_benchmark_roots"] = active_untracked_core_benchmark_roots
    undrained_plan_roots = sorted(
        name for name, state in tracked_plan_summaries.items() if state is not None and not state["drained"]
    )
    required_real_case_checkpoints = [item for item in real_case_checkpoints if item["required"]]
    incomplete_required_real_case_slugs = [
        item["slug"] for item in required_real_case_checkpoints if not item["execution_checkpoint_passed"]
    ]
    same_side_double_point_checkpoint_completed = any(
        item["execution_checkpoint_passed"] and item["checkpoint_kind"] == "same_side_double_point"
        for item in real_case_checkpoints
    )
    accepted_passed = bool(calibrated.get("accepted_calibrated_passed"))
    accepted_pearson_r = _coerce_float(calibrated.get("accepted_calibrated_pearson_r"))
    accepted_view = str(calibrated.get("accepted_calibrated_view", ""))
    accepted_excluded_complex_ids = calibrated.get("accepted_calibrated_excluded_complex_ids", [])
    selected_model = str(calibrated.get("selected_model", ""))
    full_calibrated_pearson_r = _coerce_float(calibrated.get("calibrated_pearson_r"))
    predict_pair_count = _coerce_int(calibrated.get("predict_pair_count"))
    accepted_pair_count = _coerce_int(calibrated.get("accepted_calibrated_pair_count")) or _coerce_int(
        calibrated.get("calibrated_target_filtered_pair_count")
    )
    independent_summary_path = calibrated_summary_path
    independent_generated_at = calibrated.get("generated_at", "")

    if accepted_pearson_r is None and validation_target_summary.get("available"):
        accepted_pearson_r = validation_target_summary.get("accepted_pearson_r")
        accepted_passed = bool(validation_target_summary.get("accepted_passed"))
        accepted_view = str(validation_target_summary.get("accepted_view", "")) or "target_filtered"
        accepted_excluded_complex_ids = validation_target_summary.get("accepted_excluded_complex_ids", [])
        selected_model = selected_model or str(validation_target_summary.get("selected_model", ""))
        full_calibrated_pearson_r = (
            full_calibrated_pearson_r
            if full_calibrated_pearson_r is not None
            else validation_target_summary.get("full_calibrated_pearson_r")
        )
        predict_pair_count = (
            predict_pair_count
            if predict_pair_count is not None
            else validation_target_summary.get("predict_pair_count")
        )
        accepted_pair_count = (
            accepted_pair_count
            if accepted_pair_count is not None
            else validation_target_summary.get("accepted_pair_count")
        )
        independent_summary_path = Path(str(validation_target_summary.get("summary_path")))
        independent_generated_at = str(validation_target_summary.get("generated_at", ""))

    blockers: list[str] = []
    if not accepted_passed:
        blockers.append("accepted independent validation is still below the requested R > 0.6 gate")
    if incomplete_required_real_case_slugs:
        blockers.append(
            "required real-case execution checkpoints are incomplete: "
            + ", ".join(incomplete_required_real_case_slugs)
        )
    if live_processes["stale_core_mdrun_process_count"] > 0:
        blockers.append(
            "some core AB-Bind mdrun processes appear stale "
            f"(count={live_processes['stale_core_mdrun_process_count']}, threshold={live_processes['stale_threshold_seconds']}s)"
        )
    if live_processes["stale_reference_mdrun_process_count"] > 0:
        blockers.append(
            "some Patel-like external 3HFM reference mdrun processes appear stale "
            f"(count={live_processes['stale_reference_mdrun_process_count']}, threshold={live_processes['stale_threshold_seconds']}s)"
        )
    if live_processes["core_process_count"] > 0:
        blockers.append(
            "core AB-Bind execution is still active "
            f"(mdrun={live_processes['core_mdrun_process_count']}, resume={live_processes['core_resume_process_count']})"
        )
    if live_processes["reference_process_count"] > 0:
        blockers.append(
            "Patel-like external 3HFM reference execution is still active "
            f"(mdrun={live_processes['reference_mdrun_process_count']}, resume={live_processes['reference_resume_process_count']})"
        )
    if not external_reference["completed"]:
        if external_reference["available"]:
            blockers.append(
                "Patel-like external 3HFM reference regression is incomplete "
                f"(status={external_reference['status'] or 'unknown'}, paired={external_reference['paired_job_count']}, incomplete={external_reference['incomplete_job_count']})"
            )
        else:
            blockers.append("Patel-like external 3HFM reference summary is not available yet")
    if active_untracked_core_benchmark_roots:
        blockers.append(
            "active core benchmark roots outside the tracked completion set: "
            + ", ".join(active_untracked_core_benchmark_roots)
        )
    if undrained_plan_roots:
        blockers.append(
            "tracked plan roots still report pending or running work: " + ", ".join(undrained_plan_roots)
        )
    project_complete = accepted_passed and not blockers
    return {
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "root": str(root.resolve()),
        "validation_target_r": VALIDATION_TARGET_R,
        "independent_validation": {
            "summary_path": str(independent_summary_path),
            "generated_at": independent_generated_at,
            "selected_model": selected_model,
            "accepted_view": accepted_view,
            "accepted_pearson_r": accepted_pearson_r,
            "accepted_excluded_complex_ids": accepted_excluded_complex_ids,
            "accepted_passed": accepted_passed,
            "full_calibrated_pearson_r": full_calibrated_pearson_r,
            "predict_pair_count": predict_pair_count,
            "accepted_pair_count": accepted_pair_count,
        },
        "real_case_checkpoints": real_case_checkpoints,
        "external_reference": external_reference,
        "live_processes": live_processes,
        "tracked_plan_roots": tracked_plan_summaries,
        "completion_gates": {
            "accepted_independent_validation_passed": accepted_passed,
            "required_real_case_checkpoints_completed": not incomplete_required_real_case_slugs,
            "same_side_double_point_checkpoint_completed": same_side_double_point_checkpoint_completed,
            "external_3hfm_reference_completed": external_reference["completed"],
            "no_live_core_processes": live_processes["core_process_count"] == 0,
            "no_live_reference_processes": live_processes["reference_process_count"] == 0,
            "tracked_plan_roots_drained": not undrained_plan_roots,
            "project_complete": project_complete,
        },
        "completion_blockers": blockers,
    }


def render_project_completion_status(
    *,
    root: Path,
    snapshot_date: str,
    summary: dict[str, Any],
) -> str:
    independent_validation = summary["independent_validation"]
    real_case_checkpoints = summary["real_case_checkpoints"]
    external_reference = summary["external_reference"]
    live_processes = summary["live_processes"]
    tracked_plan_roots = summary["tracked_plan_roots"]
    blockers = summary["completion_blockers"]
    lines = [
        "# Project Completion Status",
        "",
        f"Snapshot date: {_format_snapshot_date(snapshot_date)}.",
        "",
        "This file can be regenerated with:",
        f"`python benchmarks/ab_bind/report_project_completion.py --root {root}`",
        "",
        "This audit separates two questions that had been getting conflated in manual updates:",
        "whether the accepted independent validation gate is already satisfied, and whether",
        "the remaining AB-Bind and external 3HFM execution waves have fully drained.",
        "",
        "## Independent Validation Gate",
        "",
        f"- summary file: `{_relpath(root, Path(independent_validation['summary_path']))}`",
        f"- generated at: `{independent_validation.get('generated_at', '')}`",
        f"- selected calibration model: `{independent_validation.get('selected_model', '')}`",
        f"- accepted holdout view: `{independent_validation.get('accepted_view', '')}`",
        f"- accepted excluded complexes: `{_format_target_list(independent_validation.get('accepted_excluded_complex_ids'))}`",
        f"- accepted calibrated `Pearson R = {independent_validation.get('accepted_pearson_r', '')}`",
        f"- full unfiltered calibrated `Pearson R = {independent_validation.get('full_calibrated_pearson_r', '')}`",
        f"- held-out pair count: `{independent_validation.get('predict_pair_count', '')}`",
        f"- accepted pair count: `{independent_validation.get('accepted_pair_count', '')}`",
        f"- accepted gate passed: `{independent_validation.get('accepted_passed', False)}`",
        "",
        "## Real-Case Execution Gate",
        "",
        f"- required real-case checkpoints completed: `{summary['completion_gates']['required_real_case_checkpoints_completed']}`",
        f"- same-side double-point checkpoint completed: `{summary['completion_gates']['same_side_double_point_checkpoint_completed']}`",
        "",
    ]
    for item in real_case_checkpoints:
        lines.extend(
            [
                f"- {item['label']}:",
                f"  - file: `{_relpath(root, Path(item['ddg_path']))}`",
                f"  - generated at: `{item.get('generated_at', '')}`",
                f"  - mutation: `{item.get('mutation_signature', '')}`",
                f"  - mutation count / preset / side: `{item.get('mutation_count', '')}` / `{item.get('protocol_preset', '')}` / `{item.get('entity_side', '')}`",
                f"  - `ddG = {item.get('ddg_kcal_mol', '')} kcal/mol`",
                f"  - ddG BAR stderr: `{item.get('ddg_bar_stderr_kcal_mol', '')} kcal/mol`",
                f"  - QC: `{item.get('qc_status', '') or 'not available yet'}`",
                f"  - execution checkpoint passed: `{item.get('execution_checkpoint_passed', False)}`",
            ]
        )
        if item.get("qc_warnings"):
            lines.append("  - warnings:")
            for warning in item["qc_warnings"]:
                lines.append(f"    - {warning}")
    charge_class_summary = external_reference.get("charge_class_summary", {})
    charge_conserving = (
        charge_class_summary.get("charge_conserving", {}) if isinstance(charge_class_summary, dict) else {}
    )
    charge_changing = (
        charge_class_summary.get("charge_changing", {}) if isinstance(charge_class_summary, dict) else {}
    )
    lines.extend(
        [
            "",
            "## External Regression Gate",
            "",
            f"- summary file: `{_relpath(root, Path(external_reference['summary_path']))}`",
            f"- generated at: `{external_reference.get('generated_at', '')}`",
            f"- status: `{external_reference.get('status', '') or 'not available'}`",
            f"- paired jobs: `{external_reference.get('paired_job_count', 0)}`",
            f"- incomplete jobs: `{external_reference.get('incomplete_job_count', 0)}`",
            f"- charge-conserving paired / incomplete: `{_coerce_int(charge_conserving.get('paired_job_count')) or 0}` / `{_coerce_int(charge_conserving.get('incomplete_job_count')) or 0}`",
            f"- charge-changing paired / incomplete: `{_coerce_int(charge_changing.get('paired_job_count')) or 0}` / `{_coerce_int(charge_changing.get('incomplete_job_count')) or 0}`",
            f"- external reference checkpoint passed: `{external_reference.get('completed', False)}`",
        ]
    )
    if external_reference.get("message"):
        lines.append(f"- note: `{external_reference.get('message')}`")
    lines.extend(
        [
            "",
            "## Live Execution State",
            "",
        f"- active core AB-Bind `gmx mdrun` processes: `{live_processes.get('core_mdrun_process_count', 0)}`",
        f"- active core AB-Bind `abag-rbfe resume` processes: `{live_processes.get('core_resume_process_count', 0)}`",
        f"- unique core active job ids (`resume` / `mdrun`): `{live_processes.get('core_active_resume_job_count', 0)}` / `{live_processes.get('core_active_mdrun_job_count', 0)}`",
        f"- orphaned core `resume` job ids: `{_format_target_list(live_processes.get('orphaned_core_resume_job_ids'))}`",
        f"- stale core `gmx mdrun` processes (threshold `{live_processes.get('stale_threshold_seconds', STALE_PROGRESS_THRESHOLD_SECONDS)}` s): `{live_processes.get('stale_core_mdrun_process_count', 0)}`",
        f"- active core benchmark roots: `{_format_target_list(live_processes.get('active_core_benchmark_roots'))}`",
        f"- active untracked core benchmark roots: `{_format_target_list(live_processes.get('active_untracked_core_benchmark_roots'))}`",
        f"- active reference `gmx mdrun` processes: `{live_processes.get('reference_mdrun_process_count', 0)}`",
        f"- active reference `abag-rbfe resume` processes: `{live_processes.get('reference_resume_process_count', 0)}`",
        f"- orphaned reference `resume` job ids: `{_format_target_list(live_processes.get('orphaned_reference_resume_job_ids'))}`",
        f"- stale reference `gmx mdrun` processes: `{live_processes.get('stale_reference_mdrun_process_count', 0)}`",
        f"- active watcher processes: `{live_processes.get('watcher_process_count', 0)}`",
        "",
        "## Tracked Plan Drain State",
        "",
        ]
    )
    stale_core_statuses = live_processes.get("stale_core_mdrun_statuses", [])
    if stale_core_statuses:
        lines.append("Stale core processes:")
        for status in stale_core_statuses:
            lines.append(
                "- "
                f"`{status.get('job_id', '')}` pid={status.get('pid', '')} "
                f"elapsed={status.get('elapsed_seconds', '')}s "
                f"progress_age={status.get('progress_age_seconds', '')}s"
            )
        lines.append("")
    for name in sorted(tracked_plan_roots):
        state = tracked_plan_roots[name]
        if state is None:
            lines.append(f"- {name}: `not available`")
            continue
        lines.extend(
            [
                f"- {name}:",
                f"  - generated at: `{state.get('generated_at', '')}`",
                f"  - selected / ddg-ready / paired: `{state.get('selected_job_count', 0)}` / `{state.get('ddg_ready_count', 0)}` / `{state.get('paired_job_count', 0)}`",
                f"  - running sample / equilibrate: `{state.get('running_sample_job_count', 0)}` / `{state.get('running_equilibrate_job_count', 0)}`",
                f"  - pending selected jobs: `{state.get('pending_selected_job_count', 0)}`",
                f"  - drained: `{state.get('drained', False)}`",
            ]
        )
    lines.extend(
        [
            "",
            "## Completion Verdict",
            "",
            f"- accepted independent validation passed: `{summary['completion_gates']['accepted_independent_validation_passed']}`",
            f"- required real-case checkpoints completed: `{summary['completion_gates']['required_real_case_checkpoints_completed']}`",
            f"- same-side double-point checkpoint completed: `{summary['completion_gates']['same_side_double_point_checkpoint_completed']}`",
            f"- external 3HFM reference completed: `{summary['completion_gates']['external_3hfm_reference_completed']}`",
            f"- no live core processes remain: `{summary['completion_gates']['no_live_core_processes']}`",
            f"- no live reference processes remain: `{summary['completion_gates']['no_live_reference_processes']}`",
            f"- tracked plan roots drained: `{summary['completion_gates']['tracked_plan_roots_drained']}`",
            f"- project complete: `{summary['completion_gates']['project_complete']}`",
            "",
        ]
    )
    if blockers:
        lines.extend(
            [
                "Current blockers:",
                "",
                *[f"- {item}" for item in blockers],
            ]
        )
    else:
        lines.extend(
            [
                "Current blockers:",
                "",
                "- none",
            ]
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = Path(args.root).expanduser().resolve()
    summary_output = Path(args.summary_output).expanduser().resolve() if args.summary_output else default_summary_output_path(root)
    json_output = Path(args.json_output).expanduser().resolve() if args.json_output else default_json_output_path(root)
    summary = build_completion_summary(root=root)
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(
        render_project_completion_status(root=root, snapshot_date=args.snapshot_date, summary=summary),
        encoding="utf-8",
    )
    _write_json(json_output, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
