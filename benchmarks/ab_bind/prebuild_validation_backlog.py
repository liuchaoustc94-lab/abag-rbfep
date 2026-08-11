#!/usr/bin/env python3
"""Prebuild not-started validation jobs up to build_legs."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
ABAG_RBFE = ROOT / ".venv" / "bin" / "abag-rbfe"
RUNS_ROOT = ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan-root", default=str(RUNS_ROOT))
    parser.add_argument("--to-stage", default="build_legs")
    parser.add_argument("--chunk-size", type=int, default=int(os.environ.get("CHUNK_SIZE", "12")))
    parser.add_argument("--max-workers", type=int, default=int(os.environ.get("MAX_WORKERS", "4")))
    parser.add_argument("--gpu-devices", default=os.environ.get("GPU_DEVICES", ""))
    parser.add_argument("--complex-id", action="append", default=[])
    parser.add_argument("--limit-jobs", type=int)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def iter_job_dirs(plan_root: Path, complex_filters: set[str]) -> list[Path]:
    job_dirs: list[Path] = []
    for batch_dir in sorted(path for path in plan_root.glob("abbind_*_core_v1") if path.is_dir()):
        batch_name = batch_dir.name
        complex_id = batch_name.removeprefix("abbind_").removesuffix("_core_v1").upper()
        if complex_filters and complex_id not in complex_filters:
            continue
        jobs_root = batch_dir / "jobs"
        if not jobs_root.is_dir():
            continue
        for job_dir in sorted(path for path in jobs_root.iterdir() if path.is_dir()):
            job_dirs.append(job_dir)
    return job_dirs


def is_not_started(job_dir: Path) -> bool:
    stages_dir = job_dir / "stages"
    return not stages_dir.exists() or not any(stages_dir.glob("*.json"))


def chunked(values: list[str], chunk_size: int) -> list[list[str]]:
    return [values[index : index + chunk_size] for index in range(0, len(values), chunk_size)]


def main() -> int:
    args = parse_args()
    plan_root = Path(args.plan_root).expanduser().resolve()
    complex_filters = {value.strip().upper() for value in args.complex_id if value.strip()}
    job_ids = [job_dir.name for job_dir in iter_job_dirs(plan_root, complex_filters) if is_not_started(job_dir)]
    if args.limit_jobs is not None:
        job_ids = job_ids[: args.limit_jobs]

    if not job_ids:
        print("[prebuild] no not-started jobs matched the selection")
        return 0

    print(f"[prebuild] selected {len(job_ids)} not-started jobs")
    command_prefix = [
        str(ABAG_RBFE),
        "batch",
        "run-abbind",
        "--plan-root",
        str(plan_root),
        "--to-stage",
        args.to_stage,
        "--execute",
        "--max-workers",
        str(max(args.max_workers, 1)),
    ]
    if args.gpu_devices.strip():
        command_prefix.extend(["--gpu-devices", args.gpu_devices.strip()])

    for chunk_index, chunk in enumerate(chunked(job_ids, max(args.chunk_size, 1)), start=1):
        command = command_prefix + [token for job_id in chunk for token in ("--job-id", job_id)]
        print(f"[prebuild] chunk {chunk_index}: {len(chunk)} jobs")
        if args.dry_run:
            print(" ".join(command))
            continue
        result = subprocess.run(command, cwd=ROOT, check=False, text=True)
        if result.returncode != 0:
            print(f"[prebuild] chunk {chunk_index} failed with exit code {result.returncode}", file=sys.stderr)
            return result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
