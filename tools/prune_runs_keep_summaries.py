#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_KEEP_SUFFIXES = {
    ".json",
    ".yml",
    ".yaml",
    ".csv",
    ".tsv",
    ".md",
    ".txt",
    ".sh",
    ".mdp",
}

DEFAULT_KEEP_FILENAMES = {
    "LICENSE",
    "NOTICE.md",
    "README.md",
}

# Files that must survive pruning for jobs referenced by official reports
# (calibration fits, merged validation views). dhdl.xvg is the raw BAR/MBAR
# input; bar/*.xvg preserves QC histograms. Everything else stays prunable.
PROTECTED_RAW_FILENAME = "dhdl.xvg"
PROTECTED_BAR_DIR = "bar"

# Report CSVs whose (source_plan_root, batch_id, job_id) rows mark a job dir as
# protected. Merged/calibration reports are the canonical consumers of job data.
REPORT_CSV_PATTERNS = (
    "benchmarks/*/reports/calibrations/*/*.csv",
    "benchmarks/*/reports/merged/plan_jobs.csv",
    "benchmarks/*/reports/merged/selections/*/plan_jobs.csv",
    "benchmarks/*/reports/plan_jobs.csv",
    "real_cases/*/reports/*.csv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prune bulky intermediate files under runs/ while keeping summary, "
            "strategy, and result metadata."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("/mnt/data/liuchao/abag-rbfep/runs"),
        help="Runs directory to prune.",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=None,
        help="Directory where the prune report JSON will be written.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Inspect only. Do not delete files.",
    )
    parser.add_argument(
        "--sample-limit",
        type=int,
        default=5,
        help="How many example deleted files to keep per suffix in the report.",
    )
    parser.add_argument(
        "--no-protect-report-references",
        action="store_true",
        help=(
            "Disable protection of dhdl.xvg and bar/*.xvg for jobs referenced by "
            "official calibration/merged report CSVs."
        ),
    )
    return parser.parse_args()


def collect_protected_job_dirs(root: Path) -> set[Path]:
    """Resolve job directories referenced by official report CSVs.

    A job is protected when a report row carries its (source_plan_root,
    batch_id, job_id). Rows without a usable source_plan_root are ignored
    because the job dir cannot be resolved safely.
    """
    protected: set[Path] = set()
    for pattern in REPORT_CSV_PATTERNS:
        for csv_path in root.glob(pattern):
            try:
                with csv_path.open(newline="", encoding="utf-8") as handle:
                    for row in csv.DictReader(handle):
                        source_root = (row.get("source_plan_root") or "").strip()
                        batch_id = (row.get("batch_id") or "").strip()
                        job_id = (row.get("job_id") or "").strip()
                        if not source_root or not job_id:
                            continue
                        source_path = Path(source_root)
                        if not source_path.is_absolute():
                            source_path = (root / source_path)
                        job_dir = source_path / batch_id / "jobs" / job_id if batch_id else None
                        if job_dir is not None and job_dir.is_dir():
                            protected.add(job_dir.resolve())
            except (OSError, csv.Error, UnicodeDecodeError):
                continue
    return protected


def is_within(path: Path, parents: set[Path]) -> bool:
    current = path.parent
    while True:
        if current in parents:
            return True
        if current.parent == current:
            return False
        current = current.parent


def is_protected_raw_data(path: Path, protected_job_dirs: set[Path]) -> bool:
    if not protected_job_dirs:
        return False
    if path.name == PROTECTED_RAW_FILENAME:
        return is_within(path, protected_job_dirs)
    if path.suffix.lower() == ".xvg" and path.parent.name == PROTECTED_BAR_DIR:
        return is_within(path, protected_job_dirs)
    return False


def should_keep(path: Path, keep_suffixes: set[str], keep_filenames: set[str]) -> bool:
    if path.name in keep_filenames:
        return True
    if path.suffix.lower() in keep_suffixes:
        return True
    return False


def suffix_label(path: Path) -> str:
    if path.suffix:
        return path.suffix.lower()
    return "<noext>"


def iter_files(root: Path) -> Iterable[Path]:
    for path in root.rglob("*"):
        if path.is_file():
            yield path


def empty_dirs_under(root: Path) -> list[Path]:
    dirs = [path for path in root.rglob("*") if path.is_dir()]
    # Remove deepest directories first.
    dirs.sort(key=lambda path: len(path.parts), reverse=True)
    return [path for path in dirs if not any(path.iterdir())]


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    if not root.exists():
        raise SystemExit(f"runs root does not exist: {root}")

    report_dir = (
        args.report_dir.resolve()
        if args.report_dir is not None
        else (root / "cleanup_reports").resolve()
    )
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report_path = report_dir / f"prune_runs_keep_summaries_{timestamp}.json"

    keep_suffixes = set(DEFAULT_KEEP_SUFFIXES)
    keep_filenames = set(DEFAULT_KEEP_FILENAMES)
    keep_filenames.add(report_path.name)

    protected_job_dirs: set[Path] = set()
    if not args.no_protect_report_references:
        protected_job_dirs = collect_protected_job_dirs(root)
    protected_count = 0
    protected_bytes = 0

    to_delete: list[Path] = []
    kept_count = 0
    kept_bytes = 0
    deleted_counter: Counter[str] = Counter()
    deleted_bytes: defaultdict[str, int] = defaultdict(int)
    deleted_examples: defaultdict[str, list[str]] = defaultdict(list)

    for path in iter_files(root):
        if should_keep(path, keep_suffixes, keep_filenames):
            kept_count += 1
            kept_bytes += path.stat().st_size
            continue
        if is_protected_raw_data(path, protected_job_dirs):
            protected_count += 1
            protected_bytes += path.stat().st_size
            continue
        size = path.stat().st_size
        label = suffix_label(path)
        to_delete.append(path)
        deleted_counter[label] += 1
        deleted_bytes[label] += size
        if len(deleted_examples[label]) < args.sample_limit:
            deleted_examples[label].append(str(path.relative_to(root)))

    total_deleted_bytes = sum(deleted_bytes.values())
    empty_dir_candidates = empty_dirs_under(root) if not args.dry_run else []
    removed_dirs = 0

    if not args.dry_run:
        for path in to_delete:
            path.unlink(missing_ok=True)
        for path in empty_dir_candidates:
            try:
                path.rmdir()
                removed_dirs += 1
            except OSError:
                continue
        report_dir.mkdir(parents=True, exist_ok=True)

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "root": str(root),
        "dry_run": args.dry_run,
        "report_path": str(report_path),
        "keep_suffixes": sorted(keep_suffixes),
        "keep_filenames": sorted(keep_filenames),
        "report_reference_protection": not args.no_protect_report_references,
        "protected_job_dir_count": len(protected_job_dirs),
        "protected_file_count": protected_count,
        "protected_bytes": protected_bytes,
        "kept_file_count": kept_count,
        "kept_bytes": kept_bytes,
        "deleted_file_count": len(to_delete),
        "deleted_bytes": total_deleted_bytes,
        "removed_empty_dir_count": removed_dirs,
        "deleted_by_suffix": [
            {
                "suffix": suffix,
                "file_count": deleted_counter[suffix],
                "bytes": deleted_bytes[suffix],
                "sample_paths": deleted_examples[suffix],
            }
            for suffix in sorted(deleted_counter, key=lambda item: deleted_bytes[item], reverse=True)
        ],
    }

    if not args.dry_run:
        report_path.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n")

    print(json.dumps(report, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
