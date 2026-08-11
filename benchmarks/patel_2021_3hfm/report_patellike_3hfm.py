#!/mnt/data/liuchao/abag-rbfep/.venv/bin/python
"""Summarize the Patel 2021 3HFM external regression batch."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from abag_pmx.mutations import is_charge_conserving
from abag_rbfe.benchmark import _benchmark_metrics_from_pairs
from abag_rbfe.io_utils import utc_now, write_json, write_yaml
from abag_rbfe.reporting import write_batch_summary

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BATCH_DIR = ROOT / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
DEFAULT_EXPERIMENTAL_CSV = ROOT / "benchmarks" / "patel_2021_3hfm" / "experimental_ddg.csv"


def default_summary_output_path(batch_dir: Path) -> Path:
    return batch_dir / "reports" / "patel_2021_3hfm_summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-dir", default=str(DEFAULT_BATCH_DIR))
    parser.add_argument("--experimental-csv", default=str(DEFAULT_EXPERIMENTAL_CSV))
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--summary-output", default="")
    parser.add_argument("--no-write-summary", action="store_true")
    return parser.parse_args()


def _safe_float(value: Any) -> float | None:
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_json(path: Path, default: Any) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


def write_summary_outputs(summary_path: Path, payload: dict[str, Any]) -> None:
    write_json(summary_path, payload)
    write_yaml(summary_path.with_suffix(".yml"), payload)


def _charge_class(charge_conserving: bool | None) -> str:
    if charge_conserving is True:
        return "charge_conserving"
    if charge_conserving is False:
        return "charge_changing"
    return "unknown"


def _job_charge_conserving(
    batch_dir: Path,
    job: dict[str, Any],
    experimental_row: dict[str, str],
) -> bool | None:
    if "charge_conserving" in job:
        return _safe_bool(job.get("charge_conserving"))

    job_id = str(job.get("job_id", "")).strip()
    if job_id:
        job_spec = _read_json(batch_dir / "jobs" / job_id / "job_spec.json", {})
        mutation_group = job_spec.get("mutation_group") or {}
        if "charge_conserving" in mutation_group:
            return bool(mutation_group.get("charge_conserving"))

    wt = str(experimental_row.get("wt", "")).strip().upper()
    mut = str(experimental_row.get("mut", "")).strip().upper()
    if wt and mut:
        try:
            return is_charge_conserving(wt, mut)
        except ValueError:
            return None
    return None


def _build_pair_rows(
    batch_dir: Path,
    jobs: list[dict[str, Any]],
    experimental_rows: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    experimental_by_job_id = {row["job_id"]: row for row in experimental_rows}
    job_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    qc_qualified_pair_rows: list[dict[str, Any]] = []
    incomplete_rows: list[dict[str, Any]] = []
    for job in jobs:
        experimental = experimental_by_job_id.get(str(job.get("job_id", "")))
        if experimental is None:
            continue
        charge_conserving = _job_charge_conserving(batch_dir, job, experimental)
        charge_class = _charge_class(charge_conserving)
        job_rows.append(
            {
                "job_id": job.get("job_id", ""),
                "mutation_group_id": job.get("mutation_group_id", ""),
                "charge_conserving": charge_conserving,
                "charge_class": charge_class,
                "ddg_ready": _safe_bool(job.get("ddg_ready")),
                "benchmark_qc_qualified": _safe_bool(job.get("benchmark_qc_qualified")),
                "qc_status": job.get("qc_status", ""),
                "diagnostic_code": job.get("diagnostic_code", ""),
                "latest_stage": job.get("latest_stage", ""),
                "latest_stage_state": job.get("latest_stage_state", ""),
            }
        )
        predicted_ddg = _safe_float(job.get("ddg_kcal_mol"))
        experimental_ddg = _safe_float(experimental.get("experimental_ddg_kcal_mol"))
        if predicted_ddg is None or experimental_ddg is None or not _safe_bool(job.get("ddg_ready")):
            incomplete_rows.append(
                {
                    "job_id": job.get("job_id", ""),
                    "mutation_group_id": job.get("mutation_group_id", ""),
                    "qc_status": job.get("qc_status", ""),
                    "diagnostic_code": job.get("diagnostic_code", ""),
                    "latest_stage": job.get("latest_stage", ""),
                    "latest_stage_state": job.get("latest_stage_state", ""),
                    "charge_conserving": charge_conserving,
                    "charge_class": charge_class,
                    "entity_side": experimental.get("entity_side", ""),
                    "wt": experimental.get("wt", ""),
                    "mut": experimental.get("mut", ""),
                }
            )
            continue
        error = predicted_ddg - experimental_ddg
        row = {
            "job_id": job.get("job_id", ""),
            "mutation_group_id": job.get("mutation_group_id", ""),
            "predicted_ddg_kcal_mol": predicted_ddg,
            "experimental_ddg_kcal_mol": experimental_ddg,
            "ddg_error_kcal_mol": error,
            "abs_error_kcal_mol": abs(error),
            "ddg_bar_stderr_kcal_mol": _safe_float(job.get("ddg_bar_stderr_kcal_mol")),
            "qc_status": job.get("qc_status", ""),
            "benchmark_qc_qualified": _safe_bool(job.get("benchmark_qc_qualified")),
            "chain_mapping_basis": experimental.get("chain_mapping_basis", ""),
            "charge_conserving": charge_conserving,
            "charge_class": charge_class,
            "entity_side": experimental.get("entity_side", ""),
        }
        pair_rows.append(row)
        if row["benchmark_qc_qualified"]:
            qc_qualified_pair_rows.append(row)
    return job_rows, pair_rows, qc_qualified_pair_rows, incomplete_rows


def _top_abs_error_pairs(rows: list[dict[str, Any]], *, limit: int) -> list[dict[str, Any]]:
    return sorted(rows, key=lambda item: (-float(item["abs_error_kcal_mol"]), str(item["job_id"])))[: max(int(limit), 0)]


def _count_by_field(rows: list[dict[str, Any]], field: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        key = str(row.get(field, "")).strip() or "unknown"
        counts[key] = counts.get(key, 0) + 1
    return {key: counts[key] for key in sorted(counts)}


def _count_by_charge_class_and_field(rows: list[dict[str, Any]], field: str) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for row in rows:
        charge_class = str(row.get("charge_class", "")).strip() or "unknown"
        field_key = str(row.get(field, "")).strip() or "unknown"
        class_counts = counts.setdefault(charge_class, {})
        class_counts[field_key] = class_counts.get(field_key, 0) + 1
    return {
        charge_class: {key: class_counts[key] for key in sorted(class_counts)}
        for charge_class, class_counts in sorted(counts.items())
    }


def _charge_class_summary(
    job_rows: list[dict[str, Any]],
    pair_rows: list[dict[str, Any]],
    qc_qualified_pair_rows: list[dict[str, Any]],
    incomplete_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for charge_class in ("charge_conserving", "charge_changing", "unknown"):
        class_jobs = [row for row in job_rows if row.get("charge_class") == charge_class]
        class_pairs = [row for row in pair_rows if row.get("charge_class") == charge_class]
        class_qc_pairs = [row for row in qc_qualified_pair_rows if row.get("charge_class") == charge_class]
        class_incomplete = [row for row in incomplete_rows if row.get("charge_class") == charge_class]
        summary[charge_class] = {
            "job_count": len(class_jobs),
            "paired_job_count": len(class_pairs),
            "qc_qualified_pair_count": len(class_qc_pairs),
            "incomplete_job_count": len(class_incomplete),
            "raw_metrics": _benchmark_metrics_from_pairs(class_pairs) if class_pairs else {},
            "qc_qualified_metrics": _benchmark_metrics_from_pairs(class_qc_pairs) if class_qc_pairs else {},
        }
    return summary


def _insufficient_pairs_message(charge_class_summary: dict[str, Any]) -> str:
    charge_conserving_incomplete = int(
        charge_class_summary.get("charge_conserving", {}).get("incomplete_job_count") or 0
    )
    charge_changing_incomplete = int(
        charge_class_summary.get("charge_changing", {}).get("incomplete_job_count") or 0
    )
    if charge_conserving_incomplete and charge_changing_incomplete:
        return (
            "No completed Patel 2021 3HFM jobs are ready for external regression comparison yet. "
            "The charge-conserving subset is still in progress, while the charge-changing subset remains incomplete."
        )
    if charge_conserving_incomplete:
        return (
            "No completed Patel 2021 3HFM jobs are ready for external regression comparison yet. "
            "The charge-conserving subset is still in progress."
        )
    if charge_changing_incomplete:
        return (
            "No completed Patel 2021 3HFM jobs are ready for external regression comparison yet. "
            "Only charge-changing Patel rows remain incomplete."
        )
    return "No completed Patel 2021 3HFM jobs are ready for external regression comparison yet."


def main() -> int:
    args = parse_args()
    batch_dir = Path(args.batch_dir).expanduser().resolve()
    experimental_csv = Path(args.experimental_csv).expanduser().resolve()
    summary_output = (
        Path(args.summary_output).expanduser().resolve()
        if str(args.summary_output).strip()
        else default_summary_output_path(batch_dir)
    )

    batch_summary = write_batch_summary(batch_dir)
    experimental_rows = _load_csv_rows(experimental_csv)
    job_rows, pair_rows, qc_qualified_pair_rows, incomplete_rows = _build_pair_rows(
        batch_dir,
        batch_summary.get("jobs", []),
        experimental_rows,
    )
    raw_metrics = _benchmark_metrics_from_pairs(pair_rows) if pair_rows else {}
    qc_qualified_metrics = _benchmark_metrics_from_pairs(qc_qualified_pair_rows) if qc_qualified_pair_rows else {}
    charge_class_summary = _charge_class_summary(job_rows, pair_rows, qc_qualified_pair_rows, incomplete_rows)

    summary: dict[str, Any] = {
        "status": "ok" if pair_rows else "insufficient_pairs",
        "generated_at": utc_now(),
        "batch_dir": str(batch_dir),
        "experimental_csv": str(experimental_csv),
        "experimental_row_count": len(experimental_rows),
        "job_count": len(batch_summary.get("jobs", [])),
        "paired_job_count": len(pair_rows),
        "qc_qualified_pair_count": len(qc_qualified_pair_rows),
        "incomplete_job_count": len(incomplete_rows),
        "raw_metrics": raw_metrics,
        "qc_qualified_metrics": qc_qualified_metrics,
        "charge_class_summary": charge_class_summary,
        "incomplete_diagnostic_code_counts": _count_by_field(incomplete_rows, "diagnostic_code"),
        "incomplete_diagnostic_code_counts_by_charge_class": _count_by_charge_class_and_field(
            incomplete_rows,
            "diagnostic_code",
        ),
        "incomplete_stage_state_counts": _count_by_field(incomplete_rows, "latest_stage_state"),
        "top_abs_error_pairs": _top_abs_error_pairs(pair_rows, limit=args.top_n),
        "incomplete_jobs": incomplete_rows,
    }
    if not pair_rows:
        summary["message"] = _insufficient_pairs_message(charge_class_summary)

    if not args.no_write_summary:
        write_summary_outputs(summary_output, summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if pair_rows else 2


if __name__ == "__main__":
    raise SystemExit(main())
