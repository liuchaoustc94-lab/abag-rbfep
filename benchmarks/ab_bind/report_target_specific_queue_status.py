#!/mnt/data/liuchao/abag-rbfep/.venv/bin/python
"""Summarize the live target-specific queue state for 3HFM -> 1MLC handoff."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Any

from abag_rbfe.io_utils import utc_now, write_json, write_yaml

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_THREE_HFM_ROOT = (
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_target_specific_sampling_pilot_20260625"
)
DEFAULT_ONE_MLC_ROOT = (
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_target_specific_sampling_minibatch_1mlc_20260626"
)
DEFAULT_QUEUE_LOG = ROOT / "runs" / "benchmarks" / "queue_logs" / "target_specific_pilot_queue_20260626.log"
DEFAULT_JSON_OUTPUT = ROOT / "runs" / "benchmarks" / "queue_logs" / "target_specific_pilot_queue_status.json"
DEFAULT_MD_OUTPUT = ROOT / "runs" / "benchmarks" / "queue_logs" / "target_specific_pilot_queue_status.md"

THREE_HFM_JOB_IDS = [
    "3hfm-antibody-h-y33a",
    "3hfm-antibody-h-y50a",
    "3hfm-antibody-h-c95a",
    "3hfm-antigen-y-y20a",
]
ONE_MLC_JOB_IDS = [
    "1mlc-antibody-l-n92a",
    "1mlc-antibody-l-n32g",
    "1mlc-antibody-l-n32y",
]
STAGE_ORDER = ["prepare", "mutate", "build_legs", "equilibrate", "sample", "bar", "qc", "report"]
STEP_TIME_RE = re.compile(r"^\s*Step\s+Time\s*$")
STEP_VALUE_RE = re.compile(r"^\s*(\d+)\s+([0-9.]+)\s*$")
JOB_ID_RE = re.compile(r"/jobs/([^/]+)/")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--three-hfm-root", default=str(DEFAULT_THREE_HFM_ROOT))
    parser.add_argument("--one-mlc-root", default=str(DEFAULT_ONE_MLC_ROOT))
    parser.add_argument("--queue-log", default=str(DEFAULT_QUEUE_LOG))
    parser.add_argument("--json-output", default=str(DEFAULT_JSON_OUTPUT))
    parser.add_argument("--md-output", default=str(DEFAULT_MD_OUTPUT))
    return parser.parse_args()


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve_job_dir(plan_root: Path, job_id: str) -> Path | None:
    matches = list(plan_root.glob(f"**/jobs/{job_id}"))
    return matches[0] if matches else None


def _latest_stage(job_dir: Path) -> tuple[str, dict[str, Any] | None]:
    latest_name = ""
    latest_payload: dict[str, Any] | None = None
    for stage in STAGE_ORDER:
        payload = _load_json(job_dir / "stages" / f"{stage}.json")
        if payload is None:
            continue
        latest_name = stage
        latest_payload = payload
    return latest_name, latest_payload


def _tail_lines(path: Path, limit: int = 10) -> list[str]:
    if not path.is_file():
        return []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    return lines[-max(int(limit), 0) :]


def _parse_md_log_progress(md_log: Path) -> dict[str, Any]:
    if not md_log.is_file():
        return {}
    lines = md_log.read_text(encoding="utf-8", errors="replace").splitlines()
    for index in range(len(lines) - 1, -1, -1):
        if not STEP_TIME_RE.match(lines[index]):
            continue
        if index + 1 >= len(lines):
            continue
        match = STEP_VALUE_RE.match(lines[index + 1])
        if not match:
            continue
        step = int(match.group(1))
        time_ps = float(match.group(2))
        return {
            "md_log": str(md_log),
            "step": step,
            "time_ps": time_ps,
            "mtime": md_log.stat().st_mtime,
        }
    return {
        "md_log": str(md_log),
        "step": None,
        "time_ps": None,
        "mtime": md_log.stat().st_mtime,
    }


def _active_mdrun_progress(plan_root: Path) -> list[dict[str, Any]]:
    proc = subprocess.run(["ps", "-ef"], capture_output=True, text=True, check=True)
    rows: list[dict[str, Any]] = []
    root_text = str(plan_root)
    for line in proc.stdout.splitlines():
        if root_text not in line or "gmx mdrun" not in line:
            continue
        job_match = JOB_ID_RE.search(line)
        if not job_match:
            continue
        job_id = job_match.group(1)
        deffnm = ""
        tokens = line.split()
        for idx, token in enumerate(tokens):
            if token == "-deffnm" and idx + 1 < len(tokens):
                deffnm = tokens[idx + 1]
                break
        md_log = Path(f"{deffnm}.log") if deffnm else None
        progress = _parse_md_log_progress(md_log) if md_log is not None else {}
        rows.append(
            {
                "job_id": job_id,
                "command": line.strip(),
                "deffnm": deffnm,
                "progress": progress,
            }
        )
    rows.sort(key=lambda item: (item["job_id"], item.get("deffnm", "")))
    return rows


def _queue_processes() -> list[str]:
    proc = subprocess.run(["ps", "-ef"], capture_output=True, text=True, check=True)
    rows = []
    for line in proc.stdout.splitlines():
        if "run_target_specific_pilot_queue.sh" not in line:
            continue
        if "rg " in line or "grep " in line:
            continue
        rows.append(line.strip())
    return rows


def _gpu_snapshot() -> list[dict[str, Any]]:
    try:
        proc = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,utilization.gpu,memory.used,memory.total",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except Exception:
        return []
    rows: list[dict[str, Any]] = []
    for line in proc.stdout.splitlines():
        parts = [item.strip() for item in line.split(",")]
        if len(parts) != 5:
            continue
        rows.append(
            {
                "index": parts[0],
                "name": parts[1],
                "utilization_gpu": parts[2],
                "memory_used": parts[3],
                "memory_total": parts[4],
            }
        )
    return rows


def _job_status_rows(plan_root: Path, job_ids: list[str], active_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    active_by_job: dict[str, list[dict[str, Any]]] = {}
    for row in active_rows:
        active_by_job.setdefault(row["job_id"], []).append(row)

    rows: list[dict[str, Any]] = []
    for job_id in job_ids:
        row: dict[str, Any] = {"job_id": job_id, "present": False}
        job_dir = _resolve_job_dir(plan_root, job_id)
        if job_dir is None:
            rows.append(row)
            continue
        latest_stage_name, latest_stage_payload = _latest_stage(job_dir)
        row.update(
            {
                "present": True,
                "job_dir": str(job_dir),
                "latest_stage": latest_stage_name,
                "latest_stage_state": (latest_stage_payload or {}).get("state", ""),
                "latest_stage_message": (latest_stage_payload or {}).get("message", ""),
                "active_mdrun": active_by_job.get(job_id, []),
            }
        )
        ddg_summary = _load_json(job_dir / "results" / "ddg_summary.json") or {}
        qc_report = _load_json(job_dir / "results" / "qc_report.json") or {}
        if ddg_summary:
            row["ddg_ready"] = ddg_summary.get("ddg_ready")
            row["ddg_kcal_mol"] = ddg_summary.get("ddg_kcal_mol")
            row["ddg_bar_stderr_kcal_mol"] = ddg_summary.get("ddg_bar_stderr_kcal_mol")
        if qc_report:
            row["qc_status"] = qc_report.get("status", "")
            row["qc_diagnostic_code"] = qc_report.get("diagnostic_code", "")
        rows.append(row)
    return rows


def collect_status(args: argparse.Namespace) -> dict[str, Any]:
    three_hfm_root = Path(args.three_hfm_root).expanduser().resolve()
    one_mlc_root = Path(args.one_mlc_root).expanduser().resolve()
    queue_log = Path(args.queue_log).expanduser().resolve()

    three_hfm_active = _active_mdrun_progress(three_hfm_root)
    one_mlc_active = _active_mdrun_progress(one_mlc_root)
    payload: dict[str, Any] = {
        "generated_at": utc_now(),
        "three_hfm_root": str(three_hfm_root),
        "one_mlc_root": str(one_mlc_root),
        "queue_log": str(queue_log),
        "queue_processes": _queue_processes(),
        "queue_log_tail": _tail_lines(queue_log, limit=10),
        "gpus": _gpu_snapshot(),
        "three_hfm": {
            "active_mdrun_count": len(three_hfm_active),
            "active_mdruns": three_hfm_active,
            "jobs": _job_status_rows(three_hfm_root, THREE_HFM_JOB_IDS, three_hfm_active),
        },
        "one_mlc": {
            "active_mdrun_count": len(one_mlc_active),
            "active_mdruns": one_mlc_active,
            "jobs": _job_status_rows(one_mlc_root, ONE_MLC_JOB_IDS, one_mlc_active),
        },
    }
    return payload


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Target-Specific Queue Status",
        "",
        f"- Generated at: `{payload.get('generated_at', '')}`",
        f"- Queue log: `{payload.get('queue_log', '')}`",
        f"- Queue processes: `{len(payload.get('queue_processes', []))}`",
        f"- 3HFM active mdruns: `{payload.get('three_hfm', {}).get('active_mdrun_count', 0)}`",
        f"- 1MLC active mdruns: `{payload.get('one_mlc', {}).get('active_mdrun_count', 0)}`",
        "",
        "## GPU Snapshot",
        "",
        "| GPU | Util | Memory |",
        "| --- | --- | --- |",
    ]
    for gpu in payload.get("gpus", []):
        lines.append(
            f"| {gpu.get('index')} ({gpu.get('name')}) | {gpu.get('utilization_gpu')} | {gpu.get('memory_used')} / {gpu.get('memory_total')} |"
        )
    if not payload.get("gpus"):
        lines.append("| unavailable |  |  |")

    lines.extend(["", "## Queue Log Tail", ""])
    for line in payload.get("queue_log_tail", []):
        lines.append(f"- `{line}`")
    if not payload.get("queue_log_tail"):
        lines.append("- `queue log unavailable`")

    for section_key, title in (("three_hfm", "3HFM"), ("one_mlc", "1MLC")):
        section = payload.get(section_key, {})
        lines.extend(
            [
                "",
                f"## {title} Jobs",
                "",
                "| Job | Present | Latest stage | State | Active mdrun | Current ps |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in section.get("jobs", []):
            active_rows = row.get("active_mdrun", [])
            current_ps = ""
            if active_rows:
                progress = active_rows[0].get("progress", {})
                current_ps = progress.get("time_ps", "")
            lines.append(
                "| {job_id} | {present} | {stage} | {state} | {active} | {ps} |".format(
                    job_id=row.get("job_id", ""),
                    present="yes" if row.get("present") else "no",
                    stage=row.get("latest_stage", ""),
                    state=row.get("latest_stage_state", ""),
                    active="yes" if active_rows else "no",
                    ps=current_ps,
                )
            )
    return "\n".join(lines).rstrip() + "\n"


def main() -> int:
    args = parse_args()
    payload = collect_status(args)
    json_output = Path(args.json_output).expanduser().resolve()
    md_output = Path(args.md_output).expanduser().resolve()
    write_json(json_output, payload)
    write_yaml(json_output.with_suffix(".yml"), payload)
    md_output.write_text(render_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
