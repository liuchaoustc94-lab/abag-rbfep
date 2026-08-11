#!/mnt/data/liuchao/abag-rbfep/.venv/bin/python
"""Compare the same hotspot job IDs across multiple validation plan roots."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_JOB_IDS = [
    "3hfm-antibody-h-y33a",
    "3hfm-antibody-h-y50a",
    "3hfm-antibody-h-c95a",
    "3hfm-antigen-y-y20a",
]
DEFAULT_PLAN_ROOTS: list[tuple[str, Path]] = [
    ("priority", ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"),
    ("robust", ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_robust_plan"),
    (
        "targeted_repeat",
        ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_targeted_repeat_spread_rescues",
    ),
    ("targeted_lambda", ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_targeted_lambda_rescues"),
    ("sampling_qc_old", ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_sampling_qc_rescues"),
    (
        "pilot",
        ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_target_specific_sampling_pilot_20260625",
    ),
]
STAGE_ORDER = ["prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"]
PROTOCOL_KEYS = [
    "lambda_windows",
    "repeats",
    "production_ps",
    "window_relax_em_steps",
    "window_relax_md_ps",
    "nvt_ps",
    "npt_ps",
    "equilibration_restraint_schedule",
    "equilibration_release_npt_ps",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-id", action="append", default=[])
    parser.add_argument(
        "--plan-root",
        action="append",
        default=[],
        help="Optional labeled plan root in the form label=/abs/path. May be supplied multiple times.",
    )
    parser.add_argument("--json-output", default="")
    parser.add_argument("--md-output", default="")
    return parser.parse_args()


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve_job_dir(plan_root: Path, job_id: str) -> Path | None:
    matches = list(plan_root.glob(f"**/jobs/{job_id}"))
    if not matches:
        return None
    return matches[0]


def _stage_snapshot(job_dir: Path) -> dict[str, Any]:
    snapshots: list[dict[str, Any]] = []
    for stage in STAGE_ORDER:
        payload = _load_json(job_dir / "stages" / f"{stage}.json")
        if payload is None:
            continue
        snapshots.append(
            {
                "stage": stage,
                "state": payload.get("state", ""),
                "status": payload.get("status", ""),
                "message": payload.get("message", ""),
            }
        )
    latest = snapshots[-1] if snapshots else None
    return {
        "latest_stage": latest.get("stage", "") if latest else "",
        "latest_stage_state": latest.get("state", "") if latest else "",
        "latest_stage_status": latest.get("status", "") if latest else "",
        "latest_stage_message": latest.get("message", "") if latest else "",
        "stage_snapshots": snapshots,
    }


def _protocol_snapshot(job_dir: Path) -> dict[str, Any]:
    job_spec = _load_json(job_dir / "job_spec.json") or {}
    protocol = job_spec.get("protocol", {}) if isinstance(job_spec.get("protocol"), dict) else {}
    return {key: protocol.get(key) for key in PROTOCOL_KEYS}


def _result_snapshot(job_dir: Path) -> dict[str, Any]:
    ddg_summary = _load_json(job_dir / "results" / "ddg_summary.json") or {}
    qc_report = _load_json(job_dir / "results" / "qc_report.json") or {}
    return {
        "ddg_kcal_mol": ddg_summary.get("ddg_kcal_mol"),
        "ddg_bar_stderr_kcal_mol": ddg_summary.get("ddg_bar_stderr_kcal_mol"),
        "ddg_ready": ddg_summary.get("ddg_ready"),
        "qc_status": qc_report.get("status", ""),
        "qc_diagnostic_family": qc_report.get("diagnostic_family", ""),
        "qc_diagnostic_code": qc_report.get("diagnostic_code", ""),
    }


def collect_comparison(job_ids: list[str], plan_roots: list[tuple[str, Path]]) -> dict[str, Any]:
    payload: dict[str, Any] = {"job_ids": job_ids, "plan_roots": [], "jobs": {}}
    for label, path in plan_roots:
        payload["plan_roots"].append({"label": label, "path": str(path)})
    for job_id in job_ids:
        rows: list[dict[str, Any]] = []
        for label, plan_root in plan_roots:
            row: dict[str, Any] = {
                "label": label,
                "plan_root": str(plan_root),
                "present": False,
            }
            job_dir = _resolve_job_dir(plan_root, job_id)
            if job_dir is None:
                rows.append(row)
                continue
            row["present"] = True
            row["job_dir"] = str(job_dir)
            row.update(_stage_snapshot(job_dir))
            row["protocol"] = _protocol_snapshot(job_dir)
            row.update(_result_snapshot(job_dir))
            rows.append(row)
        payload["jobs"][job_id] = rows
    return payload


def render_markdown(payload: dict[str, Any]) -> str:
    lines = ["# Hotspot Root Comparison", ""]
    for job_id in payload["job_ids"]:
        lines.append(f"## {job_id}")
        lines.append("")
        lines.append("| Root | Present | Latest stage | State | ddG | stderr | QC | Protocol |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in payload["jobs"].get(job_id, []):
            protocol = row.get("protocol", {})
            protocol_summary = (
                f"rep={protocol.get('repeats')}, lam={protocol.get('lambda_windows')}, prod={protocol.get('production_ps')}, "
                f"preEM={protocol.get('window_relax_em_steps')}, preMD={protocol.get('window_relax_md_ps')}, "
                f"nvt={protocol.get('nvt_ps')}, npt={protocol.get('npt_ps')}, "
                f"restraint={protocol.get('equilibration_restraint_schedule')}, release={protocol.get('equilibration_release_npt_ps')}"
                if row.get("present")
                else ""
            )
            lines.append(
                "| {label} | {present} | {stage} | {state} | {ddg} | {stderr} | {qc} | {protocol} |".format(
                    label=row.get("label", ""),
                    present="yes" if row.get("present") else "no",
                    stage=row.get("latest_stage", ""),
                    state=row.get("latest_stage_state", ""),
                    ddg=row.get("ddg_kcal_mol", ""),
                    stderr=row.get("ddg_bar_stderr_kcal_mol", ""),
                    qc=row.get("qc_status", ""),
                    protocol=protocol_summary,
                )
            )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = parse_args()
    job_ids = [item.strip() for item in args.job_id if item.strip()] or list(DEFAULT_JOB_IDS)
    plan_roots = list(DEFAULT_PLAN_ROOTS)
    if args.plan_root:
        plan_roots = []
        for raw in args.plan_root:
            label, sep, path = raw.partition("=")
            if not sep or not label.strip() or not path.strip():
                raise SystemExit(f"Invalid --plan-root value: {raw!r}. Expected label=/abs/path")
            plan_roots.append((label.strip(), Path(path).expanduser().resolve()))
    payload = collect_comparison(job_ids, plan_roots)
    if args.json_output:
        Path(args.json_output).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    markdown = render_markdown(payload)
    if args.md_output:
        Path(args.md_output).write_text(markdown, encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
