#!/usr/bin/env python3
"""Selective cleanup of abag-rbfep run data (2026-08).

Deletes bulky MD intermediates that are NOT needed for resume, BAR recompute,
or QC audit:
  - *.xtc / *.trr / *.cpt everywhere under selected job dirs
  - lambda_*/pre_relax.* and lambda_*/pre_md.* (window relax intermediates)
  - equilibration em.* and nvt.* (npt.gro is kept as the seed-injection point)
  - setup/*.gro (regenerable solvation intermediates)

Keeps: dhdl.xvg (BAR input), md.gro/md.log/topol.tpr (window skip checks),
bar/*.xvg, all json/yaml/csv/txt/mdp/sh, npt.gro, edr.

Scope control:
  --completed-roots: entire roots treated as complete
  --active-roots:    only jobs with stages/report.json == completed are pruned
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path("/mnt/data/liuchao/abag-rbfep")

COMPLETED_ROOTS = [
    "runs/benchmarks/abbind_validation_panels_20260812",
    "runs/benchmarks/abbind_keypoints_baseline_20260806",
    "runs/benchmarks/abbind_protonation_pilot_20260806",
    "runs/benchmarks/abbind_sc_gapsys_ab_20260814",
    "runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_minibatch_1mlc_20260626",
]
ACTIVE_ROOTS = [
    "runs/benchmarks/abbind_fit_newprotocol_20260811",
    "runs/real_cases/1dvf_priority_20260806",
    "runs/real_cases/oos_1ak4_newprotocol_20260811",
    "runs/real_cases/oos_1ktz_newprotocol_20260811",
    "runs/real_cases/oos_3k2m_newprotocol_20260811",
]

DELETE_SUFFIXES = {".xtc", ".trr", ".cpt"}
DELETE_GLOBS = ("pre_relax.", "pre_md.")


def job_is_complete(job_dir: Path) -> bool:
    report = job_dir / "stages" / "report.json"
    if not report.is_file():
        return False
    try:
        return json.loads(report.read_text()).get("state") == "completed"
    except Exception:
        return False


def iter_job_dirs(root: Path, *, only_complete: bool):
    for job_dir in sorted(root.rglob("jobs/*")):
        if not job_dir.is_dir():
            continue
        if only_complete and not job_is_complete(job_dir):
            continue
        yield job_dir


def select_files(job_dir: Path):
    for path in job_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() in DELETE_SUFFIXES:
            yield path
            continue
        name = path.name
        if name.startswith(DELETE_GLOBS) and "lambda_" in str(path.parent):
            yield path
            continue
        if path.parent.name == "equilibration" and (name.startswith("em.") or name.startswith("nvt.")):
            yield path
            continue
        if path.parent.name == "setup" and path.suffix == ".gro":
            yield path
            continue


def main() -> int:
    dry_run = "--dry-run" in sys.argv
    total_bytes = 0
    total_files = 0
    per_root: dict[str, tuple[int, int]] = {}

    targets: list[tuple[Path, bool]] = []
    for rel in COMPLETED_ROOTS:
        targets.append((ROOT / rel, False))
    for rel in ACTIVE_ROOTS:
        targets.append((ROOT / rel, True))

    for root, only_complete in targets:
        root_files = 0
        root_bytes = 0
        for job_dir in iter_job_dirs(root, only_complete=only_complete):
            for path in select_files(job_dir):
                size = path.stat().st_size
                root_files += 1
                root_bytes += size
                if not dry_run:
                    path.unlink(missing_ok=True)
        per_root[str(root.relative_to(ROOT))] = (root_files, root_bytes)
        total_files += root_files
        total_bytes += root_bytes

    print(f"{'DRY-RUN ' if dry_run else ''}deleted {total_files} files, {total_bytes/1e9:.1f} GB")
    for root, (n, b) in sorted(per_root.items(), key=lambda kv: -kv[1][1]):
        print(f"  {b/1e9:7.1f} GB  {n:6d} files  {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
