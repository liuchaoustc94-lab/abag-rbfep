#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
import sys
from pathlib import Path


DEFAULT_ROOT = Path(os.environ.get("ABAG_RBFE_ROOT", "/mnt/data/liuchao/abag-rbfep"))
DEFAULT_MERGED_SUMMARY_GLOB = (
    "runs/benchmarks/abbind_core_v1_*_plan/"
    "reports/merged/selections/*/plan_summary.json"
)
DEFAULT_PRIORITY_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_priority_plan"
)
DEFAULT_ROBUST_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_robust_plan"
)
DEFAULT_RESCUE_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_priority_rescues"
)
DEFAULT_TARGETED_REPEAT_SPREAD_RESCUE_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"
)
DEFAULT_TARGETED_LAMBDA_RESCUE_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"
)
DEFAULT_SAMPLING_QC_RESCUE_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"
)
DEFAULT_DEEP_RESCUE_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_deep_rescues"
)
DEFAULT_ULTRA_RESCUE_PLAN_ROOT = (
    "runs/benchmarks/abbind_core_v1_validation_ultra_rescues"
)
DEFAULT_QUEUE_EXCLUDED_JOB_IDS = (
    "2nz9-antigen-a-h1064a",
)
DEFAULT_PROACTIVE_ROBUST_JOB_IDS = (
    "1mlc-antibody-h-s57a",
    "1mlc-antibody-h-s57v",
    "1mlc-antibody-h-t31a",
    "1mlc-antibody-h-t31v",
    "1mlc-antibody-l-n92a",
)
DEFAULT_PREFERRED_SPLIT_NAME = "validation"
_HOTSPOT_QC_FIELDS = (
    "qc_status",
    "diagnostic_family",
    "diagnostic_code",
    "primary_repeat_spread_leg",
    "repeat_spread_legs",
    "complex_repeat_spread_kcal_mol",
    "apo_repeat_spread_kcal_mol",
    "ddg_repeat_range_kcal_mol",
    "ddg_bar_stderr_kcal_mol",
    "max_bar_stderr_kcal_mol",
)
_QC_REPORT_CACHE: dict[str, dict[str, object] | None] = {}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Refresh validation hotspot watchlists from the merged validation "
            "summary's active_alternate_ready_hotspots."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Repository root used to resolve default paths.",
    )
    parser.add_argument(
        "--merged-summary-glob",
        default=DEFAULT_MERGED_SUMMARY_GLOB,
        help="Glob, relative to --root unless absolute, for merged validation plan_summary.json files.",
    )
    parser.add_argument(
        "--robust-plan-root",
        default=DEFAULT_ROBUST_PLAN_ROOT,
        help="Robust validation plan root used to detect robust alternates.",
    )
    parser.add_argument(
        "--priority-plan-root",
        default=DEFAULT_PRIORITY_PLAN_ROOT,
        help="Priority validation plan root used to discover stale resumable jobs.",
    )
    parser.add_argument(
        "--candidate-plan-root",
        action="append",
        default=[],
        help=(
            "Additional validation plan root used alongside --priority-plan-root "
            "when discovering stale/gap resumable jobs. May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--rescue-plan-root",
        default=DEFAULT_RESCUE_PLAN_ROOT,
        help="Rescue validation plan root used to detect rescue alternates.",
    )
    parser.add_argument(
        "--targeted-repeat-spread-plan-root",
        default=DEFAULT_TARGETED_REPEAT_SPREAD_RESCUE_PLAN_ROOT,
        help="Targeted repeat-spread rescue plan root scanned for currently active alternates.",
    )
    parser.add_argument(
        "--targeted-lambda-plan-root",
        default=DEFAULT_TARGETED_LAMBDA_RESCUE_PLAN_ROOT,
        help="Targeted lambda rescue plan root scanned for currently active alternates.",
    )
    parser.add_argument(
        "--deep-rescue-plan-root",
        default=DEFAULT_DEEP_RESCUE_PLAN_ROOT,
        help="Deep rescue plan root scanned for currently active alternates.",
    )
    parser.add_argument(
        "--sampling-qc-plan-root",
        default=DEFAULT_SAMPLING_QC_RESCUE_PLAN_ROOT,
        help="Sampling/QC rescue plan root scanned for currently active alternates.",
    )
    parser.add_argument(
        "--ultra-rescue-plan-root",
        default=DEFAULT_ULTRA_RESCUE_PLAN_ROOT,
        help="Ultra rescue plan root scanned for currently active alternates.",
    )
    parser.add_argument(
        "--mode",
        choices=("robust", "rescue", "targeted", "sampling_qc", "stale", "gap", "backlog", "ultra", "hotspots", "all"),
        default="all",
        help="Which watchlist to emit on stdout.",
    )
    parser.add_argument(
        "--proactive-robust-job-id",
        action="append",
        default=list(DEFAULT_PROACTIVE_ROBUST_JOB_IDS),
        help=(
            "Additional proactive robust job IDs to append after hotspot-derived "
            "jobs. May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional JSON report path capturing the chosen summary and emitted watchlists.",
    )
    parser.add_argument(
        "--preferred-split-name",
        default=DEFAULT_PREFERRED_SPLIT_NAME,
        help=(
            "Prefer merged selection summaries whose embedded selection.split_name "
            "matches this value before falling back to the latest available summary. "
            "Set to an empty string to disable the split preference."
        ),
    )
    parser.add_argument(
        "--robust-pass-outlier-threshold",
        type=float,
        default=0.0,
        help=(
            "Optional |ddG error| gate for completed QC-pass rows that already "
            "have a robust alternate candidate materialized. Matching rows are "
            "inserted into the robust watchlist after active-alternate hotspots "
            "and before the proactive robust seed set. 0 disables this fallback."
        ),
    )
    parser.add_argument(
        "--ultra-pearson-gain-threshold",
        type=float,
        default=0.2,
        help=(
            "Minimum leave-one-complex-out Pearson gain required for a hotspot to be "
            "emitted in ultra mode."
        ),
    )
    parser.add_argument(
        "--ultra-abs-error-threshold",
        type=float,
        default=5.0,
        help=(
            "Optional job-level |ddG error| gate for ultra mode. Active-alternate "
            "hotspots meeting or exceeding this threshold are also emitted even when "
            "their complex-level Pearson gain stays below --ultra-pearson-gain-threshold. "
            "0 disables the outlier gate."
        ),
    )
    parser.add_argument(
        "--ultra-pass-outlier-threshold",
        type=float,
        default=0.0,
        help=(
            "Optional |ddG error| gate for completed QC-pass rows that are not part "
            "of the current hotspot taxonomy but already have at least one alternate "
            "candidate. Matching rows are appended after hotspot-driven ultra jobs. "
            "0 disables this fallback."
        ),
    )
    parser.add_argument(
        "--sampling-qc-complex-id",
        action="append",
        default=[],
        help=(
            "Optional complex IDs used to narrow sampling/QC hotspots. "
            "May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--sampling-qc-no-active-alt-abs-error-threshold",
        type=float,
        default=0.0,
        help=(
            "Optional |ddG error| gate for completed two-leg repeat-spread rows "
            "that belong on the broader sampling/QC lane but do not currently "
            "have any active alternate. Matching rows are appended after "
            "hotspot-derived sampling/QC jobs. 0 disables this fallback."
        ),
    )
    parser.add_argument(
        "--targeted-no-active-alt-abs-error-threshold",
        type=float,
        default=0.0,
        help=(
            "Optional |ddG error| gate for completed repeat-spread rows that can "
            "use the targeted primary-leg path but do not currently have any "
            "active alternate. Matching rows are appended after hotspot-derived "
            "targeted jobs. 0 disables this fallback."
        ),
    )
    parser.add_argument(
        "--queue-excluded-job-id",
        action="append",
        default=[],
        help=(
            "Job IDs removed from emitted live watchlists by default. "
            "May be supplied multiple times."
        ),
    )
    parser.add_argument(
        "--no-default-queue-excluded-job-ids",
        action="store_true",
        help=(
            "Disable the built-in queue exclusion list so explicitly blocked jobs "
            "can still be emitted for manual follow-up."
        ),
    )
    parser.add_argument(
        "--no-derived-invalid-mutate-output-exclusions",
        action="store_true",
        help=(
            "Do not automatically exclude job IDs whose queue plan rows already "
            "report current_invalid_mutate_output."
        ),
    )
    return parser.parse_args()


def resolve_path(root: Path, value: Path | str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return root / path


def resolve_glob(root: Path, pattern: str) -> list[Path]:
    pattern_path = Path(pattern)
    if pattern_path.is_absolute():
        return sorted(Path("/").glob(str(pattern_path)[1:]))
    return sorted(root.glob(pattern))


def _summary_matches_split_name(path: Path, preferred_split_name: str) -> int:
    split_name = preferred_split_name.strip()
    if not split_name:
        return 0
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 0
    selection = payload.get("selection")
    if not isinstance(selection, dict):
        return 0
    return 1 if str(selection.get("split_name", "") or "").strip() == split_name else 0


def _summary_source_plan_root_count(path: Path) -> int:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 0
    roots = payload.get("source_plan_roots")
    if not isinstance(roots, list):
        return 0
    return sum(1 for item in roots if str(item).strip())


def _summary_paired_job_count(path: Path) -> int:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 0
    benchmark_metrics = payload.get("benchmark_metrics")
    if not isinstance(benchmark_metrics, dict):
        return 0
    value = benchmark_metrics.get("paired_job_count")
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def choose_summary_path(
    paths: list[Path],
    *,
    preferred_path: Path | None = None,
    preferred_split_name: str = "",
) -> Path | None:
    candidates = list(paths)
    if preferred_path is not None and preferred_path.exists() and preferred_path not in candidates:
        candidates.append(preferred_path)
    if not candidates:
        return None
    preferred_path_str = str(preferred_path) if preferred_path is not None else ""
    return max(
        candidates,
        key=lambda path: (
            _summary_matches_split_name(path, preferred_split_name),
            _summary_source_plan_root_count(path),
            _summary_paired_job_count(path),
            path.stat().st_mtime_ns,
            1 if preferred_path_str and str(path) == preferred_path_str else 0,
            str(path),
        ),
    )


def normalize_job_ids(job_ids: list[str]) -> list[str]:
    normalized: list[str] = []
    seen: set[str] = set()
    for job_id in job_ids:
        value = job_id.strip()
        if not value or value in seen:
            continue
        seen.add(value)
        normalized.append(value)
    return normalized


def load_summary(summary_path: Path | None) -> dict[str, object]:
    if summary_path is None or not summary_path.exists():
        return {}
    return json.loads(summary_path.read_text(encoding="utf-8"))


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _path_mtime_utc(path: Path | None) -> str:
    if path is None or not path.exists():
        return ""
    return (
        datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


def _selected_split_name(summary: dict[str, object]) -> str:
    selection = summary.get("selection")
    if not isinstance(selection, dict):
        return ""
    return str(selection.get("split_name", "") or "").strip()


def _safe_float(value: object) -> float | None:
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def _csv_tokens(value: object) -> list[str]:
    if value in ("", None):
        return []
    tokens: list[str] = []
    seen: set[str] = set()
    for raw_token in str(value).split(","):
        token = raw_token.strip()
        if not token or token in seen:
            continue
        seen.add(token)
        tokens.append(token)
    return tokens


def _normalized_string_set(values: list[str] | tuple[str, ...] | None, *, upper: bool = False) -> set[str]:
    normalized: set[str] = set()
    for value in values or []:
        item = str(value or "").strip()
        if not item:
            continue
        normalized.add(item.upper() if upper else item)
    return normalized


def queue_excluded_job_ids(args: argparse.Namespace) -> set[str]:
    excluded = set() if args.no_default_queue_excluded_job_ids else set(DEFAULT_QUEUE_EXCLUDED_JOB_IDS)
    excluded.update(_normalized_string_set(list(args.queue_excluded_job_id)))
    return excluded


def derived_invalid_mutate_output_job_ids(
    rows: list[dict[str, object]],
) -> set[str]:
    derived: set[str] = set()
    for row in rows:
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id:
            continue
        if _safe_bool(row.get("current_invalid_mutate_output")):
            derived.add(job_id)
            continue
        if str(row.get("current_invalid_mutate_output_code", "") or "").strip():
            derived.add(job_id)
    return derived


def filter_queue_excluded_rows(
    rows: list[dict[str, object]],
    *,
    excluded_job_ids: set[str],
) -> list[dict[str, object]]:
    if not excluded_job_ids:
        return list(rows)
    return [
        row
        for row in rows
        if str(row.get("job_id", "") or "").strip() not in excluded_job_ids
    ]


def benchmark_pairs_path(summary_path: Path | None) -> Path | None:
    if summary_path is None:
        return None
    candidate = summary_path.with_name("benchmark_pairs.csv")
    if candidate.exists():
        return candidate
    return None


def load_pair_rows(summary_path: Path | None) -> list[dict[str, object]]:
    pairs_path = benchmark_pairs_path(summary_path)
    if pairs_path is None:
        return []
    with pairs_path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_plan_rows(plan_root: Path) -> list[dict[str, object]]:
    plan_jobs_path = plan_root / "reports" / "plan_jobs.csv"
    if not plan_jobs_path.exists():
        return []
    with plan_jobs_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row.setdefault("source_plan_root", str(plan_root))
    return rows


def collect_plan_rows(plan_roots: list[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for plan_root in plan_roots:
        rows.extend(load_plan_rows(plan_root))
    return rows


def current_active_alternate_roots_by_job_id(
    plan_rows: list[dict[str, object]],
    *,
    active_states: set[str] | None = None,
) -> dict[str, set[str]]:
    states = active_states or {"running", "stale_running"}
    grouped: dict[str, set[str]] = {}
    for row in plan_rows:
        latest_state = str(row.get("latest_stage_state", "") or "").strip()
        if latest_state not in states:
            continue
        job_id = str(row.get("job_id", "") or "").strip()
        source_plan_root = str(row.get("source_plan_root", "") or "").strip()
        if not job_id or not source_plan_root:
            continue
        grouped.setdefault(job_id, set()).add(source_plan_root)
    return grouped


def available_alternate_roots_by_job_id(plan_rows: list[dict[str, object]]) -> dict[str, set[str]]:
    grouped: dict[str, set[str]] = {}
    for row in plan_rows:
        job_id = str(row.get("job_id", "") or "").strip()
        source_plan_root = str(row.get("source_plan_root", "") or "").strip()
        if not job_id or not source_plan_root:
            continue
        grouped.setdefault(job_id, set()).add(source_plan_root)
    return grouped


def merged_plan_jobs_path(summary_path: Path | None) -> Path | None:
    if summary_path is None:
        return None
    merged_dir = summary_path.parents[2]
    candidate = merged_dir / "plan_jobs.csv"
    if candidate.exists():
        return candidate
    return None


def load_merged_plan_rows(summary_path: Path | None) -> list[dict[str, object]]:
    plan_jobs_path = merged_plan_jobs_path(summary_path)
    if plan_jobs_path is None:
        return []
    with plan_jobs_path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _pearson(x_values: list[float], y_values: list[float]) -> float | None:
    if len(x_values) < 2 or len(x_values) != len(y_values):
        return None
    x_mean = sum(x_values) / len(x_values)
    y_mean = sum(y_values) / len(y_values)
    x_delta = [value - x_mean for value in x_values]
    y_delta = [value - y_mean for value in y_values]
    denominator = (sum(value * value for value in x_delta) * sum(value * value for value in y_delta)) ** 0.5
    if denominator == 0:
        return None
    numerator = sum(x_value * y_value for x_value, y_value in zip(x_delta, y_delta))
    return numerator / denominator


def complex_impact_metrics(pair_rows: list[dict[str, object]]) -> dict[str, dict[str, float | int | None]]:
    usable_pairs: list[tuple[str, float, float]] = []
    for row in pair_rows:
        complex_id = str(row.get("complex_id", "") or "").strip()
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if not complex_id or predicted is None or experimental is None:
            continue
        usable_pairs.append((complex_id, predicted, experimental))

    overall_pearson = _pearson(
        [predicted for _complex_id, predicted, _experimental in usable_pairs],
        [experimental for _complex_id, _predicted, experimental in usable_pairs],
    )
    by_complex: dict[str, list[tuple[float, float]]] = {}
    for complex_id, predicted, experimental in usable_pairs:
        by_complex.setdefault(complex_id, []).append((predicted, experimental))

    impacts: dict[str, dict[str, float | int | None]] = {}
    for complex_id, values in by_complex.items():
        excluded = [item for item in usable_pairs if item[0] != complex_id]
        leave_one_out_pearson = _pearson(
            [predicted for _other_complex_id, predicted, _experimental in excluded],
            [experimental for _other_complex_id, _predicted, experimental in excluded],
        )
        pearson_gain = None
        if overall_pearson is not None and leave_one_out_pearson is not None:
            pearson_gain = leave_one_out_pearson - overall_pearson
        impacts[complex_id] = {
            "pair_count": len(values),
            "overall_pearson_r": overall_pearson,
            "leave_one_out_pearson_r": leave_one_out_pearson,
            "pearson_gain": pearson_gain,
        }
    return impacts


def pair_counts_by_complex(pair_rows: list[dict[str, object]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in pair_rows:
        complex_id = str(row.get("complex_id", "") or "").strip()
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if not complex_id or predicted is None or experimental is None:
            continue
        counts[complex_id] = counts.get(complex_id, 0) + 1
    return counts


def active_alternate_roots(row: dict[str, object]) -> set[str]:
    raw = str(row.get("active_alternate_source_plan_roots", "") or "")
    return {item.strip() for item in raw.split(",") if item.strip()}


def hotspot_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    raw = summary.get("active_alternate_ready_hotspots")
    if not isinstance(raw, list):
        return []
    rows: list[dict[str, object]] = []
    for item in raw:
        if isinstance(item, dict):
            rows.append(item)
    return rows


def taxonomy_hotspot_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    taxonomy = summary.get("validation_failure_taxonomy")
    if not isinstance(taxonomy, dict):
        return []
    raw = taxonomy.get("hotspots")
    if not isinstance(raw, list):
        return []
    rows: list[dict[str, object]] = []
    for item in raw:
        if isinstance(item, dict):
            rows.append(item)
    return rows


def merged_hotspot_rows(summary: dict[str, object]) -> list[dict[str, object]]:
    active_rows = hotspot_rows(summary)
    taxonomy_rows = taxonomy_hotspot_rows(summary)
    if not taxonomy_rows:
        return active_rows

    active_by_job_id = {
        str(row.get("job_id", "") or "").strip(): row
        for row in active_rows
        if str(row.get("job_id", "") or "").strip()
    }
    merged: list[dict[str, object]] = []
    seen: set[str] = set()

    for row in taxonomy_rows:
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id:
            continue
        seen.add(job_id)
        merged_row = dict(active_by_job_id.get(job_id, {}))
        for key, value in row.items():
            if value not in ("", None, [], {}):
                merged_row[key] = value
            elif key not in merged_row:
                merged_row[key] = value
        merged.append(merged_row)

    for row in active_rows:
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id or job_id in seen:
            continue
        merged.append(dict(row))
    return merged


def annotate_hotspots_with_plan_rows(
    hotspots: list[dict[str, object]],
    *,
    merged_plan_rows: list[dict[str, object]],
    queue_plan_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    merged_job_ids = {
        str(row.get("job_id", "") or "").strip()
        for row in merged_plan_rows
        if str(row.get("job_id", "") or "").strip()
    }
    plan_rows = list(merged_plan_rows) + [
        row
        for row in queue_plan_rows
        if str(row.get("job_id", "") or "").strip() not in merged_job_ids
    ]
    plan_rows_by_job_id = {
        str(row.get("job_id", "") or "").strip(): row
        for row in plan_rows
        if str(row.get("job_id", "") or "").strip()
    }
    annotated: list[dict[str, object]] = []
    for row in hotspots:
        item = dict(row)
        job_id = str(item.get("job_id", "") or "").strip()
        plan_row = plan_rows_by_job_id.get(job_id)
        if plan_row is not None:
            for field in _HOTSPOT_QC_FIELDS:
                if item.get(field) not in ("", None, [], {}):
                    continue
                value = plan_row.get(field)
                if value in ("", None, [], {}):
                    continue
                item[field] = value
        item["repeat_spread_legs"] = ",".join(_csv_tokens(item.get("repeat_spread_legs")))
        annotated.append(item)
    return annotated


def hotspot_sampling_qc_payload(hotspots: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    payload: dict[str, dict[str, object]] = {}
    for row in hotspots:
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id:
            continue
        payload[job_id] = {
            "qc_status": str(row.get("qc_status", "") or ""),
            "diagnostic_family": str(row.get("diagnostic_family", "") or ""),
            "diagnostic_code": str(row.get("diagnostic_code", "") or ""),
            "primary_repeat_spread_leg": str(row.get("primary_repeat_spread_leg", "") or ""),
            "repeat_spread_legs": _csv_tokens(row.get("repeat_spread_legs")),
            "complex_repeat_spread_kcal_mol": _safe_float(row.get("complex_repeat_spread_kcal_mol")),
            "apo_repeat_spread_kcal_mol": _safe_float(row.get("apo_repeat_spread_kcal_mol")),
            "ddg_repeat_range_kcal_mol": _row_ddg_repeat_range_kcal_mol(row),
            "ddg_bar_stderr_kcal_mol": _safe_float(row.get("ddg_bar_stderr_kcal_mol")),
            "max_bar_stderr_kcal_mol": _safe_float(row.get("max_bar_stderr_kcal_mol")),
            "prefer_targeted_primary_repeat_spread_leg": _prefer_targeted_primary_repeat_spread_leg(row),
        }
    return payload


def _qc_report_path_for_row(row: dict[str, object]) -> Path | None:
    source_plan_root = str(row.get("source_plan_root", "") or "").strip()
    batch_id = str(row.get("batch_id", "") or "").strip()
    job_id = str(row.get("job_id", "") or "").strip()
    if not source_plan_root or not batch_id or not job_id:
        return None
    return Path(source_plan_root) / batch_id / "jobs" / job_id / "results" / "qc_report.json"


def _qc_report_for_row(row: dict[str, object]) -> dict[str, object] | None:
    path = _qc_report_path_for_row(row)
    if path is None:
        return None
    cache_key = str(path.resolve())
    if cache_key in _QC_REPORT_CACHE:
        return _QC_REPORT_CACHE[cache_key]
    if not path.is_file():
        _QC_REPORT_CACHE[cache_key] = None
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = None
    _QC_REPORT_CACHE[cache_key] = payload
    return payload


def _qc_low_overlap_legs(qc_report: dict[str, object]) -> list[str]:
    overlap_threshold = _safe_float(qc_report.get("overlap_threshold"))
    if overlap_threshold is None:
        return []
    overlap_legs = qc_report.get("overlap_assessment", {}).get("legs", {})
    if not isinstance(overlap_legs, dict):
        return []

    failing_legs: list[str] = []
    for leg_name in ("complex", "apo"):
        payload = overlap_legs.get(leg_name, {})
        if not isinstance(payload, dict):
            continue
        score = _safe_float(payload.get("overlap_score_min"))
        if score is not None and score < overlap_threshold:
            failing_legs.append(leg_name)
    return failing_legs


def _qc_has_dominant_primary_repeat_spread_leg(qc_report: dict[str, object]) -> bool:
    repeat_spread_legs = list(qc_report.get("repeat_spread_legs", []))
    primary_repeat_spread_leg = str(qc_report.get("primary_repeat_spread_leg") or "")
    if primary_repeat_spread_leg not in {"complex", "apo"} or len(repeat_spread_legs) < 2:
        return False
    if primary_repeat_spread_leg not in repeat_spread_legs:
        return False

    secondary_legs = [leg for leg in repeat_spread_legs if leg != primary_repeat_spread_leg]
    if len(secondary_legs) != 1:
        return False

    legs = qc_report.get("legs", {})
    if not isinstance(legs, dict):
        return False
    primary_repeat_spread = _safe_float(
        legs.get(primary_repeat_spread_leg, {}).get("repeat_delta_kcal_mol_range")
    )
    secondary_repeat_spread = _safe_float(
        legs.get(secondary_legs[0], {}).get("repeat_delta_kcal_mol_range")
    )
    if primary_repeat_spread is None or secondary_repeat_spread is None:
        return False

    margin_threshold = max(_safe_float(qc_report.get("max_repeat_delta_kcal_mol")) or 0.0, 1.0)
    return (
        primary_repeat_spread >= secondary_repeat_spread + margin_threshold
        and primary_repeat_spread >= secondary_repeat_spread * 1.5
    )


def _qc_rescue_reason_codes(qc_report: dict[str, object]) -> list[str]:
    reasons: list[str] = []
    ddg_repeat_range = _safe_float(qc_report.get("ddg_repeat_range_kcal_mol"))
    max_repeat_delta = _safe_float(qc_report.get("max_repeat_delta_kcal_mol"))
    if ddg_repeat_range is not None and max_repeat_delta is not None and ddg_repeat_range > max_repeat_delta:
        reasons.append("repeat_spread")

    ddg_bar_stderr = _safe_float(qc_report.get("ddg_bar_stderr_kcal_mol"))
    max_bar_stderr = _safe_float(qc_report.get("max_bar_stderr_kcal_mol"))
    if ddg_bar_stderr is not None and max_bar_stderr is not None and ddg_bar_stderr > max_bar_stderr:
        reasons.append("bar_stderr")

    overlap_threshold = _safe_float(qc_report.get("overlap_threshold"))
    overlap_legs = qc_report.get("overlap_assessment", {}).get("legs", {})
    if overlap_threshold is not None and isinstance(overlap_legs, dict):
        for payload in overlap_legs.values():
            if not isinstance(payload, dict):
                continue
            score = _safe_float(payload.get("overlap_score_min"))
            if score is not None and score < overlap_threshold:
                reasons.append("overlap")
                break

    return reasons


def _qc_can_target_primary_repeat_spread_leg(qc_report: dict[str, object]) -> bool:
    rescue_reasons = _qc_rescue_reason_codes(qc_report)
    repeat_spread_legs = list(qc_report.get("repeat_spread_legs", []))
    primary_repeat_spread_leg = str(qc_report.get("primary_repeat_spread_leg") or "")
    if "repeat_spread" not in rescue_reasons:
        return False
    if primary_repeat_spread_leg not in {"complex", "apo"}:
        return False

    unexpected_reasons = set(rescue_reasons) - {"repeat_spread", "overlap"}
    if unexpected_reasons:
        return False

    if len(repeat_spread_legs) == 1:
        repeat_spread_targetable = True
    else:
        repeat_spread_targetable = _qc_has_dominant_primary_repeat_spread_leg(qc_report)
    if not repeat_spread_targetable:
        return False

    if "overlap" in rescue_reasons:
        low_overlap_legs = _qc_low_overlap_legs(qc_report)
        if set(low_overlap_legs) != {primary_repeat_spread_leg}:
            return False
    return True


def _row_has_dominant_primary_repeat_spread_leg(row: dict[str, object]) -> bool:
    repeat_spread_legs = _csv_tokens(row.get("repeat_spread_legs"))
    primary_repeat_spread_leg = str(row.get("primary_repeat_spread_leg", "") or "").strip()
    if primary_repeat_spread_leg not in {"complex", "apo"} or len(repeat_spread_legs) < 2:
        return False
    if primary_repeat_spread_leg not in repeat_spread_legs:
        return False

    secondary_legs = [leg for leg in repeat_spread_legs if leg != primary_repeat_spread_leg]
    if len(secondary_legs) != 1:
        return False

    primary_repeat_spread = _safe_float(row.get(f"{primary_repeat_spread_leg}_repeat_spread_kcal_mol"))
    secondary_repeat_spread = _safe_float(row.get(f"{secondary_legs[0]}_repeat_spread_kcal_mol"))
    if primary_repeat_spread is None or secondary_repeat_spread is None:
        return False

    margin_threshold = max(_safe_float(row.get("max_repeat_delta_kcal_mol")) or 0.0, 1.0)
    return (
        primary_repeat_spread >= secondary_repeat_spread + margin_threshold
        and primary_repeat_spread >= secondary_repeat_spread * 1.5
    )


def _row_repeat_spread_kcal_mol(row: dict[str, object], leg: str) -> float | None:
    normalized_leg = str(leg or "").strip().lower()
    if normalized_leg not in {"complex", "apo"}:
        return None
    value = _safe_float(row.get(f"{normalized_leg}_repeat_spread_kcal_mol"))
    if value is not None:
        return value
    qc_report = _qc_report_for_row(row)
    if not isinstance(qc_report, dict):
        return None
    legs = qc_report.get("legs", {})
    if not isinstance(legs, dict):
        return None
    payload = legs.get(normalized_leg, {})
    if not isinstance(payload, dict):
        return None
    return _safe_float(payload.get("repeat_delta_kcal_mol_range"))


def _row_ddg_repeat_range_kcal_mol(row: dict[str, object]) -> float | None:
    value = _safe_float(row.get("ddg_repeat_range_kcal_mol"))
    if value is not None:
        return value
    qc_report = _qc_report_for_row(row)
    if not isinstance(qc_report, dict):
        return None
    return _safe_float(qc_report.get("ddg_repeat_range_kcal_mol"))


def _descending_sort_value(value: float | None, *, missing: float = -1.0) -> float:
    return -(value if value is not None else missing)


def _targeted_primary_repeat_spread_sort_key(
    row: dict[str, object],
    *,
    index: int,
) -> tuple[float, float, int, float, float, int]:
    primary_repeat_spread_leg = str(row.get("primary_repeat_spread_leg", "") or "").strip().lower()
    repeat_spread_legs = _csv_tokens(row.get("repeat_spread_legs"))
    primary_repeat_spread = _row_repeat_spread_kcal_mol(row, primary_repeat_spread_leg)
    secondary_repeat_spreads = [
        _row_repeat_spread_kcal_mol(row, leg) or 0.0
        for leg in repeat_spread_legs
        if leg != primary_repeat_spread_leg
    ]
    secondary_repeat_spread = max(secondary_repeat_spreads, default=0.0)
    ddg_repeat_range = _row_ddg_repeat_range_kcal_mol(row)
    abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
    return (
        _descending_sort_value(primary_repeat_spread),
        _descending_sort_value(ddg_repeat_range),
        0 if len(repeat_spread_legs) == 1 else 1,
        _descending_sort_value(
            (primary_repeat_spread - secondary_repeat_spread) if primary_repeat_spread is not None else None
        ),
        _descending_sort_value(abs_error),
        index,
    )


def _sampling_qc_sort_key(
    row: dict[str, object],
    *,
    index: int,
) -> tuple[float, float, float, float, int]:
    complex_repeat_spread = _row_repeat_spread_kcal_mol(row, "complex")
    apo_repeat_spread = _row_repeat_spread_kcal_mol(row, "apo")
    observed_repeat_spreads = [
        value for value in (complex_repeat_spread, apo_repeat_spread) if value is not None
    ]
    floor_repeat_spread = min(observed_repeat_spreads) if observed_repeat_spreads else None
    total_repeat_spread = sum(observed_repeat_spreads) if observed_repeat_spreads else None
    ddg_repeat_range = _row_ddg_repeat_range_kcal_mol(row)
    abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
    return (
        _descending_sort_value(floor_repeat_spread),
        _descending_sort_value(total_repeat_spread),
        _descending_sort_value(ddg_repeat_range),
        _descending_sort_value(abs_error),
        index,
    )


def _prefer_targeted_primary_repeat_spread_leg(row: dict[str, object]) -> bool:
    if str(row.get("diagnostic_code", "") or "").strip() != "qc_repeat_spread":
        return False
    if len(_csv_tokens(row.get("repeat_spread_legs"))) < 2:
        return False
    qc_report = _qc_report_for_row(row)
    if isinstance(qc_report, dict):
        return _qc_can_target_primary_repeat_spread_leg(qc_report)
    return _row_has_dominant_primary_repeat_spread_leg(row)


def _is_targeted_primary_repeat_spread_candidate(row: dict[str, object]) -> bool:
    if str(row.get("qc_status", "") or "").strip() == "pass":
        return False
    qc_report = _qc_report_for_row(row)
    if isinstance(qc_report, dict):
        return _qc_can_target_primary_repeat_spread_leg(qc_report)

    if str(row.get("diagnostic_code", "") or "").strip() != "qc_repeat_spread":
        return False
    primary_repeat_spread_leg = str(row.get("primary_repeat_spread_leg", "") or "").strip()
    if primary_repeat_spread_leg not in {"complex", "apo"}:
        return False
    repeat_spread_legs = _csv_tokens(row.get("repeat_spread_legs"))
    if len(repeat_spread_legs) == 1:
        return primary_repeat_spread_leg in repeat_spread_legs
    if len(repeat_spread_legs) < 2:
        return False
    return _row_has_dominant_primary_repeat_spread_leg(row)


def targeted_primary_repeat_spread_candidates(
    hotspots: list[dict[str, object]],
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for row in hotspots:
        if not _is_targeted_primary_repeat_spread_candidate(row):
            continue
        selected.append(dict(row))
    return [
        row
        for _index, row in sorted(
            enumerate(selected),
            key=lambda item: _targeted_primary_repeat_spread_sort_key(item[1], index=item[0]),
        )
    ]


def targeted_no_active_alternate_outlier_candidates(
    plan_rows: list[dict[str, object]],
    *,
    hotspot_job_ids: set[str],
    current_active_alternate_roots: dict[str, set[str]],
    threshold: float,
    targeted_repeat_spread_plan_root: Path,
) -> list[dict[str, object]]:
    if threshold <= 0:
        return []

    targeted_plan_root_str = str(targeted_repeat_spread_plan_root)
    selected: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id or job_id in hotspot_job_ids:
            continue
        if not _safe_bool(row.get("ddg_ready")):
            continue
        source_plan_root = str(row.get("source_plan_root", "") or "").strip()
        if source_plan_root == targeted_plan_root_str:
            continue
        abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
        if abs_error is None or abs_error < threshold:
            continue
        if active_alternate_roots(row):
            continue
        if _has_current_active_alternate_elsewhere(
            row,
            current_active_alternate_roots=current_active_alternate_roots,
        ):
            continue
        if not _is_targeted_primary_repeat_spread_candidate(row):
            continue
        item = dict(row)
        item["_original_index"] = index
        item["targeted_outlier_reason"] = "no_active_alternate_abs_error"
        selected.append(item)

    ordered = sorted(
        selected,
        key=lambda item: _targeted_primary_repeat_spread_sort_key(item, index=int(item["_original_index"])),
    )
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def sampling_qc_candidates(
    hotspots: list[dict[str, object]],
    *,
    complex_ids: set[str],
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for row in hotspots:
        complex_id = str(row.get("complex_id", "") or "").strip().upper()
        if complex_ids and complex_id not in complex_ids:
            continue
        if str(row.get("qc_status", "") or "").strip() == "pass":
            continue
        if str(row.get("diagnostic_code", "") or "").strip() != "qc_repeat_spread":
            continue
        repeat_spread_legs = _csv_tokens(row.get("repeat_spread_legs"))
        if len(repeat_spread_legs) < 2:
            continue
        primary_repeat_spread_leg = str(row.get("primary_repeat_spread_leg", "") or "").strip()
        if primary_repeat_spread_leg not in {"complex", "apo"}:
            continue
        if _prefer_targeted_primary_repeat_spread_leg(row):
            continue
        selected.append(dict(row))
    return [
        row
        for _index, row in sorted(
            enumerate(selected),
            key=lambda item: _sampling_qc_sort_key(item[1], index=item[0]),
        )
    ]


def sampling_qc_no_active_alternate_outlier_candidates(
    plan_rows: list[dict[str, object]],
    *,
    hotspot_job_ids: set[str],
    current_active_alternate_roots: dict[str, set[str]],
    complex_ids: set[str],
    threshold: float,
    sampling_qc_plan_root: Path,
) -> list[dict[str, object]]:
    if threshold <= 0:
        return []

    sampling_qc_plan_root_str = str(sampling_qc_plan_root)
    selected: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id or job_id in hotspot_job_ids:
            continue
        if not _safe_bool(row.get("ddg_ready")):
            continue
        complex_id = str(row.get("complex_id", "") or "").strip().upper()
        if complex_ids and complex_id not in complex_ids:
            continue
        source_plan_root = str(row.get("source_plan_root", "") or "").strip()
        if source_plan_root == sampling_qc_plan_root_str:
            continue
        abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
        if abs_error is None or abs_error < threshold:
            continue
        if active_alternate_roots(row):
            continue
        if _has_current_active_alternate_elsewhere(
            row,
            current_active_alternate_roots=current_active_alternate_roots,
        ):
            continue
        if str(row.get("qc_status", "") or "").strip() == "pass":
            continue
        if str(row.get("diagnostic_code", "") or "").strip() != "qc_repeat_spread":
            continue
        repeat_spread_legs = _csv_tokens(row.get("repeat_spread_legs"))
        if len(repeat_spread_legs) < 2:
            continue
        primary_repeat_spread_leg = str(row.get("primary_repeat_spread_leg", "") or "").strip()
        if primary_repeat_spread_leg not in {"complex", "apo"}:
            continue
        if _prefer_targeted_primary_repeat_spread_leg(row):
            continue
        item = dict(row)
        item["_original_index"] = index
        item["sampling_qc_outlier_reason"] = "no_active_alternate_abs_error"
        selected.append(item)

    ordered = sorted(
        selected,
        key=lambda item: _sampling_qc_sort_key(item, index=int(item["_original_index"])),
    )
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def sort_hotspots(
    hotspots: list[dict[str, object]],
    *,
    pair_impacts: dict[str, dict[str, float | int | None]],
) -> list[dict[str, object]]:
    annotated: list[dict[str, object]] = []
    for index, row in enumerate(hotspots):
        item = dict(row)
        complex_id = str(item.get("complex_id", "") or "").strip()
        impact = pair_impacts.get(complex_id, {})
        item["complex_impact_pair_count"] = impact.get("pair_count")
        item["complex_impact_overall_pearson_r"] = impact.get("overall_pearson_r")
        item["complex_impact_leave_one_out_pearson_r"] = impact.get("leave_one_out_pearson_r")
        item["complex_impact_pearson_gain"] = impact.get("pearson_gain")
        item["_original_index"] = index
        annotated.append(item)

    def _sort_key(row: dict[str, object]) -> tuple[float, float, float, int]:
        pearson_gain = _safe_float(row.get("complex_impact_pearson_gain"))
        abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
        pair_count = _safe_float(row.get("complex_impact_pair_count"))
        return (
            -(pearson_gain if pearson_gain is not None else float("-inf")),
            -(pair_count if pair_count is not None else 0.0),
            -(abs_error if abs_error is not None else 0.0),
            int(row["_original_index"]),
        )

    ordered = sorted(annotated, key=_sort_key)
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def _has_current_active_alternate_elsewhere(
    row: dict[str, object],
    *,
    current_active_alternate_roots: dict[str, set[str]],
) -> bool:
    job_id = str(row.get("job_id", "") or "").strip()
    if not job_id:
        return False
    source_plan_root = str(row.get("source_plan_root", "") or "").strip()
    return bool(current_active_alternate_roots.get(job_id, set()) - {source_plan_root})


def ultra_hotspot_candidates(
    hotspots: list[dict[str, object]],
    *,
    current_active_alternate_roots: dict[str, set[str]],
    ultra_pearson_gain_threshold: float,
    ultra_abs_error_threshold: float,
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for row in hotspots:
        has_active_alternate = bool(active_alternate_roots(row))
        has_current_active_alternate = _has_current_active_alternate_elsewhere(
            row,
            current_active_alternate_roots=current_active_alternate_roots,
        )
        if not has_active_alternate and not has_current_active_alternate:
            continue
        pearson_gain = _safe_float(row.get("complex_impact_pearson_gain")) or 0.0
        abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol")) or 0.0
        if pearson_gain < ultra_pearson_gain_threshold and (
            ultra_abs_error_threshold <= 0 or abs_error < ultra_abs_error_threshold
        ):
            continue
        selected.append(dict(row))
    return selected


def ultra_pass_outlier_candidates(
    plan_rows: list[dict[str, object]],
    *,
    hotspot_job_ids: set[str],
    current_active_alternate_roots: dict[str, set[str]],
    threshold: float,
    ultra_rescue_plan_root: Path,
) -> list[dict[str, object]]:
    if threshold <= 0:
        return []

    selected: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id or job_id in hotspot_job_ids:
            continue
        if not _safe_bool(row.get("ddg_ready")):
            continue
        if str(row.get("qc_status", "") or "").strip() != "pass":
            continue
        source_plan_root = str(row.get("source_plan_root", "") or "").strip()
        if source_plan_root == str(ultra_rescue_plan_root):
            continue
        abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
        if abs_error is None or abs_error < threshold:
            continue
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental is None:
            continue
        alternate_candidate_count = int(_safe_float(row.get("alternate_candidate_count")) or 0)
        has_current_active_alternate = _has_current_active_alternate_elsewhere(
            row,
            current_active_alternate_roots=current_active_alternate_roots,
        )
        if alternate_candidate_count <= 0 and not has_current_active_alternate:
            continue
        item = dict(row)
        item["_original_index"] = index
        item["ultra_outlier_reason"] = "pass_qc_abs_error"
        item["ultra_has_current_active_alternate"] = has_current_active_alternate
        selected.append(item)

    ordered = sorted(
        selected,
        key=lambda row: (
            -(_safe_float(row.get("abs_ddg_error_kcal_mol")) or 0.0),
            0 if _safe_bool(row.get("ultra_has_current_active_alternate")) else 1,
            -(int(_safe_float(row.get("alternate_candidate_count")) or 0)),
            int(row["_original_index"]),
        ),
    )
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def robust_pass_outlier_candidates(
    plan_rows: list[dict[str, object]],
    *,
    hotspot_job_ids: set[str],
    alternate_roots_by_job_id: dict[str, set[str]],
    current_active_alternate_roots: dict[str, set[str]],
    threshold: float,
    robust_plan_root: Path,
) -> list[dict[str, object]]:
    if threshold <= 0:
        return []

    robust_plan_root_str = str(robust_plan_root)
    selected: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id or job_id in hotspot_job_ids:
            continue
        if not _safe_bool(row.get("ddg_ready")):
            continue
        if str(row.get("qc_status", "") or "").strip() != "pass":
            continue
        source_plan_root = str(row.get("source_plan_root", "") or "").strip()
        if source_plan_root == robust_plan_root_str:
            continue
        abs_error = _safe_float(row.get("abs_ddg_error_kcal_mol"))
        if abs_error is None or abs_error < threshold:
            continue
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental is None:
            continue
        alternate_roots = set(alternate_roots_by_job_id.get(job_id, set())) - {source_plan_root}
        if robust_plan_root_str not in alternate_roots:
            continue
        item = dict(row)
        item["_original_index"] = index
        item["robust_outlier_reason"] = "pass_qc_abs_error"
        item["robust_has_current_active_alternate"] = _has_current_active_alternate_elsewhere(
            row,
            current_active_alternate_roots=current_active_alternate_roots,
        )
        selected.append(item)

    ordered = sorted(
        selected,
        key=lambda row: (
            -(_safe_float(row.get("abs_ddg_error_kcal_mol")) or 0.0),
            0 if _safe_bool(row.get("robust_has_current_active_alternate")) else 1,
            int(row["_original_index"]),
        ),
    )
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def unique_rows_by_job_id(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    seen: set[str] = set()
    for row in rows:
        job_id = str(row.get("job_id", "") or "").strip()
        if not job_id or job_id in seen:
            continue
        seen.add(job_id)
        selected.append(row)
    return selected


def stale_candidates(
    plan_rows: list[dict[str, object]],
    *,
    pair_counts: dict[str, int],
    hotspot_job_ids: set[str],
) -> list[dict[str, object]]:
    stage_priority = {
        "sample": 0,
        "equilibrate": 1,
        "build_legs": 2,
        "mutate": 3,
        "prepare": 4,
        "ingest": 5,
    }
    candidates: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        latest_state = str(row.get("latest_stage_state", "") or "").strip()
        if latest_state != "stale_running":
            continue
        if not _safe_bool(row.get("resumable")):
            continue
        if _safe_bool(row.get("ddg_ready")):
            continue
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental is None:
            continue
        item = dict(row)
        complex_id = str(item.get("complex_id", "") or "").strip()
        job_id = str(item.get("job_id", "") or "").strip()
        item["current_pair_count"] = pair_counts.get(complex_id, 0)
        item["abs_experimental_ddg_kcal_mol"] = abs(experimental)
        item["is_hotspot_job"] = job_id in hotspot_job_ids
        item["_original_index"] = index
        candidates.append(item)

    def _sort_key(row: dict[str, object]) -> tuple[float, int, int, float, int]:
        pair_count = _safe_float(row.get("current_pair_count"))
        latest_stage = str(row.get("latest_stage", "") or "").strip()
        abs_experimental = _safe_float(row.get("abs_experimental_ddg_kcal_mol"))
        return (
            pair_count if pair_count is not None else 0.0,
            0 if _safe_bool(row.get("is_hotspot_job")) else 1,
            stage_priority.get(latest_stage, 99),
            -(abs_experimental if abs_experimental is not None else 0.0),
            int(row["_original_index"]),
        )

    ordered = sorted(candidates, key=_sort_key)
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def gap_candidates(
    plan_rows: list[dict[str, object]],
    *,
    pair_counts: dict[str, int],
) -> list[dict[str, object]]:
    stage_priority = {
        "sample": 0,
        "equilibrate": 1,
        "build_legs": 2,
        "mutate": 3,
        "prepare": 4,
        "ingest": 5,
        "": 6,
    }
    state_priority = {
        "stale_running": 0,
        "completed": 1,
        "not_started": 2,
        "": 3,
    }
    candidates: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        if _safe_bool(row.get("ddg_ready")):
            continue
        if not _safe_bool(row.get("resumable")):
            continue
        latest_state = str(row.get("latest_stage_state", "") or "").strip()
        if latest_state == "running":
            continue
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental is None:
            continue
        complex_id = str(row.get("complex_id", "") or "").strip()
        current_pair_count = pair_counts.get(complex_id, 0)
        if current_pair_count != 0:
            continue
        item = dict(row)
        item["current_pair_count"] = current_pair_count
        item["abs_experimental_ddg_kcal_mol"] = abs(experimental)
        item["_original_index"] = index
        candidates.append(item)

    def _sort_key(row: dict[str, object]) -> tuple[float, int, int, float, int]:
        pair_count = _safe_float(row.get("current_pair_count"))
        latest_stage = str(row.get("latest_stage", "") or "").strip()
        latest_state = str(row.get("latest_stage_state", "") or "").strip()
        abs_experimental = _safe_float(row.get("abs_experimental_ddg_kcal_mol"))
        return (
            pair_count if pair_count is not None else 0.0,
            stage_priority.get(latest_stage, 99),
            state_priority.get(latest_state, 99),
            -(abs_experimental if abs_experimental is not None else 0.0),
            int(row["_original_index"]),
        )

    ordered = sorted(candidates, key=_sort_key)
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def backlog_candidates(
    plan_rows: list[dict[str, object]],
    *,
    pair_counts: dict[str, int],
) -> list[dict[str, object]]:
    stage_priority = {
        "sample": 0,
        "equilibrate": 1,
        "build_legs": 2,
        "mutate": 3,
        "prepare": 4,
        "ingest": 5,
        "": 6,
    }
    state_priority = {
        "stale_running": 0,
        "running": 1,
        "completed": 2,
        "not_started": 3,
        "": 4,
    }
    candidates: list[dict[str, object]] = []
    for index, row in enumerate(plan_rows):
        if _safe_bool(row.get("ddg_ready")):
            continue
        if not _safe_bool(row.get("resumable", True)):
            continue
        latest_state = str(row.get("latest_stage_state", "") or "").strip()
        if latest_state == "blocked_input":
            continue
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental is None:
            continue
        active_alternates = int(_safe_float(row.get("active_alternate_candidate_count")) or 0)
        if active_alternates != 0:
            continue
        complex_id = str(row.get("complex_id", "") or "").strip()
        item = dict(row)
        item["current_pair_count"] = pair_counts.get(complex_id, 0)
        item["active_alternate_candidate_count"] = active_alternates
        item["abs_experimental_ddg_kcal_mol"] = abs(experimental)
        item["_original_index"] = index
        candidates.append(item)

    def _sort_key(row: dict[str, object]) -> tuple[float, float, int, int, float, int]:
        pair_count = _safe_float(row.get("current_pair_count"))
        active_alternates = _safe_float(row.get("active_alternate_candidate_count"))
        latest_stage = str(row.get("latest_stage", "") or "").strip()
        latest_state = str(row.get("latest_stage_state", "") or "").strip()
        abs_experimental = _safe_float(row.get("abs_experimental_ddg_kcal_mol"))
        return (
            pair_count if pair_count is not None else 0.0,
            active_alternates if active_alternates is not None else 0.0,
            stage_priority.get(latest_stage, 99),
            state_priority.get(latest_state, 99),
            -(abs_experimental if abs_experimental is not None else 0.0),
            int(row["_original_index"]),
        )

    ordered = sorted(candidates, key=_sort_key)
    for row in ordered:
        row.pop("_original_index", None)
    return ordered


def build_watchlists(
    *,
    hotspots: list[dict[str, object]],
    robust_pass_outlier_rows: list[dict[str, object]],
    targeted_rows: list[dict[str, object]],
    sampling_qc_rows: list[dict[str, object]],
    ultra_rows: list[dict[str, object]],
    stale_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    backlog_rows: list[dict[str, object]],
    current_active_alternate_roots: dict[str, set[str]],
    robust_plan_root: Path,
    rescue_plan_root: Path,
    proactive_robust_job_ids: list[str],
    queue_excluded_job_ids: set[str],
) -> dict[str, list[str]]:
    robust_plan_root_str = str(robust_plan_root)
    rescue_plan_root_str = str(rescue_plan_root)

    def _filter(job_ids: list[str]) -> list[str]:
        normalized = normalize_job_ids(job_ids)
        if not queue_excluded_job_ids:
            return normalized
        return [job_id for job_id in normalized if job_id not in queue_excluded_job_ids]

    robust_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in hotspots
            if robust_plan_root_str in active_alternate_roots(row)
        ]
        + [
            str(row.get("job_id", "") or "")
            for row in robust_pass_outlier_rows
        ]
        + proactive_robust_job_ids
    )
    rescue_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in hotspots
            if rescue_plan_root_str in active_alternate_roots(row)
        ]
    )
    targeted_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in targeted_rows
        ]
    )
    sampling_qc_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in sampling_qc_rows
        ]
    )
    ultra_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in ultra_rows
        ]
    )
    stale_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in stale_rows
        ]
    )
    gap_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in gap_rows
        ]
    )
    excluded_backlog_jobs = set(stale_jobs) | set(gap_jobs)
    backlog_jobs = _filter(
        [
            str(row.get("job_id", "") or "")
            for row in backlog_rows
            if str(row.get("job_id", "") or "") not in excluded_backlog_jobs
        ]
    )
    return {
        "robust": robust_jobs,
        "rescue": rescue_jobs,
        "targeted": targeted_jobs,
        "sampling_qc": sampling_qc_jobs,
        "stale": stale_jobs,
        "gap": gap_jobs,
        "backlog": backlog_jobs,
        "ultra": ultra_jobs,
    }


def write_summary_report(
    *,
    output_json: Path,
    summary: dict[str, object],
    summary_path: Path | None,
    active_alternate_hotspots: list[dict[str, object]],
    hotspots: list[dict[str, object]],
    robust_pass_outlier_rows: list[dict[str, object]],
    targeted_no_active_alternate_rows: list[dict[str, object]],
    sampling_qc_no_active_alternate_rows: list[dict[str, object]],
    stale_rows: list[dict[str, object]],
    gap_rows: list[dict[str, object]],
    backlog_rows: list[dict[str, object]],
    ultra_rows: list[dict[str, object]],
    watchlists: dict[str, list[str]],
    pair_impacts: dict[str, dict[str, float | int | None]],
    current_active_alternate_roots: dict[str, set[str]],
    robust_pass_outlier_threshold: float,
    sampling_qc_no_active_alt_abs_error_threshold: float,
    targeted_no_active_alt_abs_error_threshold: float,
    ultra_pearson_gain_threshold: float,
    ultra_abs_error_threshold: float,
    ultra_pass_outlier_threshold: float,
    queue_excluded_job_ids: set[str],
) -> None:
    output_json.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "watchlist_generated_at": _utc_now_iso(),
        "selected_summary_path": str(summary_path) if summary_path is not None else "",
        "selected_summary_generated_at": str(summary.get("generated_at", "") or "").strip(),
        "selected_summary_mtime_utc": _path_mtime_utc(summary_path),
        "selected_summary_selection_split_name": _selected_split_name(summary),
        "selected_benchmark_pairs_path": str(benchmark_pairs_path(summary_path) or ""),
        "active_alternate_ready_hotspot_count": len(active_alternate_hotspots),
        "active_alternate_ready_hotspot_job_ids": [
            str(row.get("job_id", "") or "")
            for row in active_alternate_hotspots
            if str(row.get("job_id", "") or "").strip()
        ],
        "robust_pass_outlier_threshold": robust_pass_outlier_threshold,
        "robust_pass_outlier_job_ids": [
            str(row.get("job_id", "") or "")
            for row in robust_pass_outlier_rows
            if str(row.get("robust_outlier_reason", "") or "").strip() == "pass_qc_abs_error"
            and str(row.get("job_id", "") or "").strip()
        ],
        "sampling_qc_no_active_alt_abs_error_threshold": sampling_qc_no_active_alt_abs_error_threshold,
        "sampling_qc_no_active_alt_outlier_job_ids": [
            str(row.get("job_id", "") or "")
            for row in sampling_qc_no_active_alternate_rows
            if str(row.get("sampling_qc_outlier_reason", "") or "").strip() == "no_active_alternate_abs_error"
            and str(row.get("job_id", "") or "").strip()
        ],
        "targeted_no_active_alt_abs_error_threshold": targeted_no_active_alt_abs_error_threshold,
        "targeted_no_active_alt_outlier_job_ids": [
            str(row.get("job_id", "") or "")
            for row in targeted_no_active_alternate_rows
            if str(row.get("targeted_outlier_reason", "") or "").strip() == "no_active_alternate_abs_error"
            and str(row.get("job_id", "") or "").strip()
        ],
        "validation_failure_hotspot_count": len(hotspots),
        "validation_failure_hotspot_job_ids": [
            str(row.get("job_id", "") or "") for row in hotspots if str(row.get("job_id", "") or "").strip()
        ],
        "hotspot_sampling_qc_by_job_id": hotspot_sampling_qc_payload(hotspots),
        "hotspot_complex_impact_pearson_gain": {
            str(row.get("job_id", "") or ""): _safe_float(row.get("complex_impact_pearson_gain"))
            for row in hotspots
            if str(row.get("job_id", "") or "").strip()
        },
        "current_active_alternate_roots_by_job_id": {
            job_id: sorted(roots)
            for job_id, roots in sorted(current_active_alternate_roots.items())
            if roots
        },
        "complex_impact_by_complex_id": pair_impacts,
        "robust_job_ids": watchlists["robust"],
        "rescue_job_ids": watchlists["rescue"],
        "targeted_job_ids": watchlists["targeted"],
        "sampling_qc_job_ids": watchlists["sampling_qc"],
        "sampling_qc_excluded_targeted_primary_repeat_spread_leg_job_ids": [
            str(row.get("job_id", "") or "")
            for row in hotspots
            if str(row.get("job_id", "") or "").strip()
            and _prefer_targeted_primary_repeat_spread_leg(row)
        ],
        "stale_job_ids": watchlists["stale"],
        "stale_job_current_pair_count": {
            str(row.get("job_id", "") or ""): int(_safe_float(row.get("current_pair_count")) or 0)
            for row in stale_rows
            if str(row.get("job_id", "") or "").strip()
        },
        "gap_job_ids": watchlists["gap"],
        "gap_job_current_pair_count": {
            str(row.get("job_id", "") or ""): int(_safe_float(row.get("current_pair_count")) or 0)
            for row in gap_rows
            if str(row.get("job_id", "") or "").strip()
        },
        "backlog_job_ids": watchlists["backlog"],
        "backlog_job_current_pair_count": {
            str(row.get("job_id", "") or ""): int(_safe_float(row.get("current_pair_count")) or 0)
            for row in backlog_rows
            if str(row.get("job_id", "") or "").strip()
        },
        "backlog_job_active_alternate_count": {
            str(row.get("job_id", "") or ""): int(_safe_float(row.get("active_alternate_candidate_count")) or 0)
            for row in backlog_rows
            if str(row.get("job_id", "") or "").strip()
        },
        "gap_complex_ids": sorted(
            {
                str(row.get("complex_id", "") or "").strip()
                for row in gap_rows
                if str(row.get("complex_id", "") or "").strip()
            }
        ),
        "ultra_pearson_gain_threshold": ultra_pearson_gain_threshold,
        "ultra_abs_error_threshold": ultra_abs_error_threshold,
        "ultra_pass_outlier_threshold": ultra_pass_outlier_threshold,
        "ultra_job_ids": watchlists["ultra"],
        "ultra_pass_outlier_job_ids": [
            str(row.get("job_id", "") or "")
            for row in ultra_rows
            if str(row.get("ultra_outlier_reason", "") or "").strip() == "pass_qc_abs_error"
            and str(row.get("job_id", "") or "").strip()
        ],
        "queue_excluded_job_ids": sorted(queue_excluded_job_ids),
    }
    output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def emit_watchlist(mode: str, watchlists: dict[str, list[str]]) -> int:
    if mode == "hotspots":
        for job_id in watchlists["hotspots"]:
            sys.stdout.write(f"{job_id}\n")
        return 0
    if mode == "all":
        json.dump(watchlists, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return 0
    for job_id in watchlists[mode]:
        sys.stdout.write(f"{job_id}\n")
    return 0


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    summary_paths = resolve_glob(root, args.merged_summary_glob)
    canonical_summary_path = resolve_path(root, args.priority_plan_root) / "reports" / "merged" / "plan_summary.json"
    preferred_summary_path = canonical_summary_path if args.merged_summary_glob == DEFAULT_MERGED_SUMMARY_GLOB else None
    summary_path = choose_summary_path(
        summary_paths,
        preferred_path=preferred_summary_path,
        preferred_split_name=str(args.preferred_split_name or ""),
    )
    summary = load_summary(summary_path)
    pair_rows = load_pair_rows(summary_path)
    merged_plan_rows = load_merged_plan_rows(summary_path)
    pair_impacts = complex_impact_metrics(pair_rows)
    pair_counts = pair_counts_by_complex(pair_rows)
    queue_plan_roots = [resolve_path(root, args.priority_plan_root)]
    robust_plan_root = resolve_path(root, args.robust_plan_root)
    rescue_plan_root = resolve_path(root, args.rescue_plan_root)
    targeted_repeat_spread_plan_root = resolve_path(root, args.targeted_repeat_spread_plan_root)
    targeted_lambda_plan_root = resolve_path(root, args.targeted_lambda_plan_root)
    sampling_qc_plan_root = resolve_path(root, args.sampling_qc_plan_root)
    deep_rescue_plan_root = resolve_path(root, args.deep_rescue_plan_root)
    ultra_rescue_plan_root = resolve_path(root, args.ultra_rescue_plan_root)
    for item in args.candidate_plan_root:
        resolved = resolve_path(root, item)
        if resolved not in queue_plan_roots:
            queue_plan_roots.append(resolved)
    queue_plan_rows = collect_plan_rows(queue_plan_roots)
    active_hotspots = annotate_hotspots_with_plan_rows(
        sort_hotspots(hotspot_rows(summary), pair_impacts=pair_impacts),
        merged_plan_rows=merged_plan_rows,
        queue_plan_rows=queue_plan_rows,
    )
    hotspots = annotate_hotspots_with_plan_rows(
        sort_hotspots(merged_hotspot_rows(summary), pair_impacts=pair_impacts),
        merged_plan_rows=merged_plan_rows,
        queue_plan_rows=queue_plan_rows,
    )
    hotspot_job_ids = {
        str(row.get("job_id", "") or "").strip()
        for row in hotspots
        if str(row.get("job_id", "") or "").strip()
    }
    alternate_scan_plan_roots = list(queue_plan_roots)
    for active_only_plan_root in (
        robust_plan_root,
        rescue_plan_root,
        targeted_repeat_spread_plan_root,
        targeted_lambda_plan_root,
        sampling_qc_plan_root,
        deep_rescue_plan_root,
        ultra_rescue_plan_root,
    ):
        if active_only_plan_root not in alternate_scan_plan_roots:
            alternate_scan_plan_roots.append(active_only_plan_root)
    alternate_scan_rows = collect_plan_rows(alternate_scan_plan_roots)
    current_active_alternates = current_active_alternate_roots_by_job_id(alternate_scan_rows)
    available_alternates = available_alternate_roots_by_job_id(alternate_scan_rows)
    queue_excluded = queue_excluded_job_ids(args)
    if not args.no_derived_invalid_mutate_output_exclusions:
        queue_excluded.update(derived_invalid_mutate_output_job_ids(queue_plan_rows))
    stale_rows = stale_candidates(
        queue_plan_rows,
        pair_counts=pair_counts,
        hotspot_job_ids=hotspot_job_ids,
    )
    stale_rows = filter_queue_excluded_rows(stale_rows, excluded_job_ids=queue_excluded)
    gap_rows = gap_candidates(
        queue_plan_rows,
        pair_counts=pair_counts,
    )
    gap_rows = filter_queue_excluded_rows(gap_rows, excluded_job_ids=queue_excluded)
    backlog_rows = backlog_candidates(
        merged_plan_rows or queue_plan_rows,
        pair_counts=pair_counts,
    )
    backlog_rows = filter_queue_excluded_rows(backlog_rows, excluded_job_ids=queue_excluded)
    targeted_rows = targeted_primary_repeat_spread_candidates(hotspots)
    targeted_no_active_alternate_rows = targeted_no_active_alternate_outlier_candidates(
        merged_plan_rows or queue_plan_rows,
        hotspot_job_ids=hotspot_job_ids,
        current_active_alternate_roots=current_active_alternates,
        threshold=float(args.targeted_no_active_alt_abs_error_threshold),
        targeted_repeat_spread_plan_root=targeted_repeat_spread_plan_root,
    )
    targeted_no_active_alternate_rows = filter_queue_excluded_rows(
        targeted_no_active_alternate_rows,
        excluded_job_ids=queue_excluded,
    )
    targeted_rows = unique_rows_by_job_id(targeted_rows + targeted_no_active_alternate_rows)
    targeted_rows = filter_queue_excluded_rows(targeted_rows, excluded_job_ids=queue_excluded)
    sampling_qc_rows = sampling_qc_candidates(
        hotspots,
        complex_ids=_normalized_string_set(list(args.sampling_qc_complex_id), upper=True),
    )
    sampling_qc_no_active_alternate_rows = sampling_qc_no_active_alternate_outlier_candidates(
        merged_plan_rows or queue_plan_rows,
        hotspot_job_ids=hotspot_job_ids,
        current_active_alternate_roots=current_active_alternates,
        complex_ids=_normalized_string_set(list(args.sampling_qc_complex_id), upper=True),
        threshold=float(args.sampling_qc_no_active_alt_abs_error_threshold),
        sampling_qc_plan_root=sampling_qc_plan_root,
    )
    sampling_qc_no_active_alternate_rows = filter_queue_excluded_rows(
        sampling_qc_no_active_alternate_rows,
        excluded_job_ids=queue_excluded,
    )
    sampling_qc_rows = unique_rows_by_job_id(sampling_qc_rows + sampling_qc_no_active_alternate_rows)
    sampling_qc_rows = filter_queue_excluded_rows(sampling_qc_rows, excluded_job_ids=queue_excluded)
    ultra_rows = ultra_hotspot_candidates(
        hotspots,
        current_active_alternate_roots=current_active_alternates,
        ultra_pearson_gain_threshold=float(args.ultra_pearson_gain_threshold),
        ultra_abs_error_threshold=float(args.ultra_abs_error_threshold),
    )
    ultra_pass_outlier_rows = ultra_pass_outlier_candidates(
        merged_plan_rows or queue_plan_rows,
        hotspot_job_ids=hotspot_job_ids,
        current_active_alternate_roots=current_active_alternates,
        threshold=float(args.ultra_pass_outlier_threshold),
        ultra_rescue_plan_root=ultra_rescue_plan_root,
    )
    ultra_rows = unique_rows_by_job_id(ultra_rows + ultra_pass_outlier_rows)
    ultra_rows = filter_queue_excluded_rows(ultra_rows, excluded_job_ids=queue_excluded)
    robust_pass_outlier_rows = robust_pass_outlier_candidates(
        merged_plan_rows or queue_plan_rows,
        hotspot_job_ids=hotspot_job_ids,
        alternate_roots_by_job_id=available_alternates,
        current_active_alternate_roots=current_active_alternates,
        threshold=float(args.robust_pass_outlier_threshold),
        robust_plan_root=robust_plan_root,
    )
    robust_pass_outlier_rows = filter_queue_excluded_rows(
        robust_pass_outlier_rows,
        excluded_job_ids=queue_excluded,
    )
    watchlists = build_watchlists(
        hotspots=hotspots,
        robust_pass_outlier_rows=robust_pass_outlier_rows,
        targeted_rows=targeted_rows,
        sampling_qc_rows=sampling_qc_rows,
        ultra_rows=ultra_rows,
        stale_rows=stale_rows,
        gap_rows=gap_rows,
        backlog_rows=backlog_rows,
        current_active_alternate_roots=current_active_alternates,
        robust_plan_root=robust_plan_root,
        rescue_plan_root=rescue_plan_root,
        proactive_robust_job_ids=normalize_job_ids(list(args.proactive_robust_job_id)),
        queue_excluded_job_ids=queue_excluded,
    )
    watchlists["hotspots"] = [
        job_id
        for job_id in normalize_job_ids([str(row.get("job_id", "") or "") for row in hotspots])
        if job_id not in queue_excluded
    ]

    if args.output_json is not None:
        output_json = resolve_path(root, args.output_json)
        write_summary_report(
            output_json=output_json,
            summary=summary,
            summary_path=summary_path,
            active_alternate_hotspots=active_hotspots,
            hotspots=hotspots,
            robust_pass_outlier_rows=robust_pass_outlier_rows,
            targeted_no_active_alternate_rows=targeted_no_active_alternate_rows,
            sampling_qc_no_active_alternate_rows=sampling_qc_no_active_alternate_rows,
            stale_rows=stale_rows,
            gap_rows=gap_rows,
            backlog_rows=backlog_rows,
            ultra_rows=ultra_rows,
            watchlists=watchlists,
            pair_impacts=pair_impacts,
            current_active_alternate_roots=current_active_alternates,
            robust_pass_outlier_threshold=float(args.robust_pass_outlier_threshold),
            sampling_qc_no_active_alt_abs_error_threshold=float(args.sampling_qc_no_active_alt_abs_error_threshold),
            targeted_no_active_alt_abs_error_threshold=float(args.targeted_no_active_alt_abs_error_threshold),
            ultra_pearson_gain_threshold=float(args.ultra_pearson_gain_threshold),
            ultra_abs_error_threshold=float(args.ultra_abs_error_threshold),
            ultra_pass_outlier_threshold=float(args.ultra_pass_outlier_threshold),
            queue_excluded_job_ids=queue_excluded,
        )

    return emit_watchlist(args.mode, watchlists)


if __name__ == "__main__":
    raise SystemExit(main())
