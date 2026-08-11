#!/usr/bin/env python3
"""Terminate active plan-root processes that fall outside an explicit allow-list."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
from typing import Any

JOB_RE = re.compile(r"/jobs/([^/]+)/")
CLI_JOB_RE = re.compile(r"\babag-rbfe\b\s+(?:run|resume|analyze)\s+([^\s]+)")
MDRUN_RE = re.compile(r"(?:^|\s)\S*gmx(?:_mpi)?\s+mdrun(?:\s|$)")
ACTIVE_PATH_MARKERS = ("/artifacts/commands/", "/legs/")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan-root", required=True, help="Plan root whose active processes should be audited.")
    parser.add_argument(
        "--signal",
        default="TERM",
        choices=("TERM", "KILL"),
        help="Signal used when --execute is set. Defaults to TERM.",
    )
    parser.add_argument("--execute", action="store_true", help="Actually send the selected signal to outside-scope PIDs.")
    parser.add_argument("job_id", nargs="+", help="Allow-listed job IDs that may continue running under the plan root.")
    return parser.parse_args()


def _canonical_job_id_from_text(text: str) -> str | None:
    path_match = JOB_RE.search(text)
    if path_match:
        return path_match.group(1)
    cli_match = CLI_JOB_RE.search(text)
    if cli_match:
        return cli_match.group(1)
    return None


def _process_kind(command: str) -> str:
    if MDRUN_RE.search(command):
        return "mdrun"
    if CLI_JOB_RE.search(command):
        return "controller"
    if any(marker in command for marker in ACTIVE_PATH_MARKERS):
        return "stage_script"
    return "job_process"


def read_ps_output() -> str:
    result = subprocess.run(["ps", "-ef"], check=False, capture_output=True, text=True)
    return result.stdout if result.returncode == 0 else ""


def list_plan_processes(plan_root: Path, *, ps_output: str | None = None) -> list[dict[str, Any]]:
    resolved_root = str(plan_root.expanduser().resolve())
    output = ps_output if ps_output is not None else read_ps_output()
    processes: list[dict[str, Any]] = []
    for line in output.splitlines():
        if resolved_root not in line:
            continue
        parts = line.split(None, 7)
        if len(parts) < 8 or not parts[1].isdigit():
            continue
        command = parts[7]
        job_id = _canonical_job_id_from_text(command)
        if not job_id:
            continue
        processes.append(
            {
                "pid": int(parts[1]),
                "job_id": job_id,
                "kind": _process_kind(command),
                "command": command,
            }
        )
    processes.sort(key=lambda item: (str(item["job_id"]), int(item["pid"])))
    return processes


def outside_scope_processes(
    plan_root: Path,
    *,
    allowed_job_ids: set[str],
    ps_output: str | None = None,
) -> list[dict[str, Any]]:
    allowed = {item.strip() for item in allowed_job_ids if item.strip()}
    return [item for item in list_plan_processes(plan_root, ps_output=ps_output) if item["job_id"] not in allowed]


def terminate_processes(
    processes: list[dict[str, Any]],
    *,
    sig: signal.Signals,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    terminated: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    for item in processes:
        try:
            os.kill(int(item["pid"]), sig)
        except ProcessLookupError:
            continue
        except OSError as exc:
            errors.append(
                {
                    "pid": item["pid"],
                    "job_id": item["job_id"],
                    "kind": item["kind"],
                    "error": str(exc),
                }
            )
            continue
        terminated.append(item)
    return terminated, errors


def main() -> int:
    args = parse_args()
    plan_root = Path(args.plan_root).expanduser().resolve()
    allowed_job_ids = {item.strip() for item in args.job_id if item.strip()}
    processes = outside_scope_processes(plan_root, allowed_job_ids=allowed_job_ids)
    payload = {
        "plan_root": str(plan_root),
        "allowed_job_ids": sorted(allowed_job_ids),
        "outside_scope_process_count": len(processes),
        "outside_scope_processes": processes,
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    if not args.execute or not processes:
        return 0

    sig = signal.SIGTERM if args.signal == "TERM" else signal.SIGKILL
    terminated, errors = terminate_processes(processes, sig=sig)
    result = {
        "signal": args.signal,
        "terminated_process_count": len(terminated),
        "terminated_processes": terminated,
        "error_count": len(errors),
        "errors": errors,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
