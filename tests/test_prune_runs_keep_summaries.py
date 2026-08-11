"""Tests for tools/prune_runs_keep_summaries.py report-reference protection."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
TOOL = REPO_ROOT / "tools" / "prune_runs_keep_summaries.py"


def _write_job(root: Path, plan_root_name: str, batch_id: str, job_id: str) -> Path:
    job_dir = root / plan_root_name / batch_id / "jobs" / job_id
    dhdl_dir = job_dir / "legs" / "complex" / "rep01" / "lambda_000"
    dhdl_dir.mkdir(parents=True)
    (dhdl_dir / "dhdl.xvg").write_text("# dhdl\n0 0\n", encoding="utf-8")
    bar_dir = job_dir / "legs" / "complex" / "rep01" / "bar"
    bar_dir.mkdir(parents=True)
    (bar_dir / "histogram.xvg").write_text("0 1\n", encoding="utf-8")
    (job_dir / "legs" / "complex" / "rep01" / "lambda_000" / "md.xtc").write_bytes(b"traj")
    (job_dir / "job_spec.json").write_text("{}", encoding="utf-8")
    return job_dir


def _write_pairs_csv(csv_path: Path, rows: list[dict[str, str]]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["complex_id", "batch_id", "source_plan_root", "job_id"]
        )
        writer.writeheader()
        writer.writerows(rows)


def _run_tool(root: Path, *extra: str) -> dict:
    result = subprocess.run(
        [sys.executable, str(TOOL), "--root", str(root), *extra],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


def test_prune_protects_report_referenced_job_raw_data(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    plan_root = root / "benchmarks" / "plan_a"
    protected_job = _write_job(root, "benchmarks/plan_a", "batch1", "job-protected")
    plain_job = _write_job(root, "benchmarks/plan_a", "batch1", "job-plain")

    _write_pairs_csv(
        plan_root / "reports" / "calibrations" / "fit-x" / "predict_pairs_calibrated.csv",
        [
            {
                "complex_id": "CMPX",
                "batch_id": "batch1",
                "source_plan_root": str(plan_root),
                "job_id": "job-protected",
            }
        ],
    )

    report = _run_tool(root)

    assert (protected_job / "legs" / "complex" / "rep01" / "lambda_000" / "dhdl.xvg").is_file()
    assert (protected_job / "legs" / "complex" / "rep01" / "bar" / "histogram.xvg").is_file()
    # Non-essential bulk inside a protected job is still pruned.
    assert not (protected_job / "legs" / "complex" / "rep01" / "lambda_000" / "md.xtc").exists()
    # Text summaries always survive.
    assert (protected_job / "job_spec.json").is_file()
    # Unreferenced jobs lose everything except text summaries.
    assert not (plain_job / "legs" / "complex" / "rep01" / "lambda_000" / "dhdl.xvg").exists()
    assert not (plain_job / "legs" / "complex" / "rep01" / "bar" / "histogram.xvg").exists()
    assert (plain_job / "job_spec.json").is_file()

    assert report["report_reference_protection"] is True
    assert report["protected_job_dir_count"] == 1
    assert report["protected_file_count"] == 2


def test_prune_no_protect_flag_disables_protection(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    plan_root = root / "benchmarks" / "plan_a"
    protected_job = _write_job(root, "benchmarks/plan_a", "batch1", "job-protected")
    _write_pairs_csv(
        plan_root / "reports" / "calibrations" / "fit-x" / "predict_pairs_calibrated.csv",
        [
            {
                "complex_id": "CMPX",
                "batch_id": "batch1",
                "source_plan_root": str(plan_root),
                "job_id": "job-protected",
            }
        ],
    )

    report = _run_tool(root, "--no-protect-report-references")

    assert not (protected_job / "legs" / "complex" / "rep01" / "lambda_000" / "dhdl.xvg").exists()
    assert report["report_reference_protection"] is False
    assert report["protected_file_count"] == 0


def test_prune_dry_run_counts_protection_without_deleting(tmp_path: Path) -> None:
    root = tmp_path / "runs"
    plan_root = root / "benchmarks" / "plan_a"
    protected_job = _write_job(root, "benchmarks/plan_a", "batch1", "job-protected")
    _write_pairs_csv(
        plan_root / "reports" / "merged" / "plan_jobs.csv",
        [
            {
                "complex_id": "CMPX",
                "batch_id": "batch1",
                "source_plan_root": str(plan_root),
                "job_id": "job-protected",
            }
        ],
    )

    report = _run_tool(root, "--dry-run")

    assert (protected_job / "legs" / "complex" / "rep01" / "lambda_000" / "dhdl.xvg").is_file()
    assert report["dry_run"] is True
    assert report["protected_file_count"] == 2
    assert report["deleted_file_count"] == 1  # only md.xtc
