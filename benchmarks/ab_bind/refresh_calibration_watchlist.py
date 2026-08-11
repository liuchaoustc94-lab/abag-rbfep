#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

from abag_rbfe.benchmark import report_ab_bind_plan
from abag_rbfe.io_utils import read_csv_rows, read_json, write_json


DEFAULT_ROOT = Path(os.environ.get("ABAG_RBFE_ROOT", "/mnt/data/liuchao/abag-rbfep"))
DEFAULT_PLAN_ROOT = DEFAULT_ROOT / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
DEFAULT_SPLIT_FILE = DEFAULT_ROOT / "benchmarks" / "ab_bind" / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
DEFAULT_SUMMARY_PATH = DEFAULT_PLAN_ROOT / "reports" / "calibrated_validation_summary.json"
DEFAULT_OUTPUT_JSON = DEFAULT_PLAN_ROOT / "reports" / "watch" / "calibration_watchlist_refresh.json"
DEFAULT_TARGET_JOB_COUNT = 8
DEFAULT_FIT_SPLIT_NAMES = ("development", "calibration")
DEFAULT_FALLBACK_JOB_IDS = (
    "1dqj-antibody-h-y50a",
    "1dqj-antibody-h-y33a",
    "1n8z-antibody-h-w95a",
    "3be1-antibody-h-y33a",
    "1dqj-antigen-c-y20a",
    "2jel-antigen-p-s64t",
    "2nyy-antigen-a-f953a",
    "3bn9-antigen-a-f97a",
)
_ACTIVE_STATE_RANK = {
    "running": 0,
    "stale_running": 1,
    "failed": 4,
    "blocked_input": 5,
    "not_started": 3,
    "completed": 2,
}
_FIT_SPLIT_RANK = {"calibration": 0, "development": 1}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Refresh the calibration quick-watch watchlist so fit coverage keeps "
            "expanding across development and calibration splits."
        )
    )
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--plan-root", type=Path, default=DEFAULT_PLAN_ROOT)
    parser.add_argument("--split-file", type=Path, default=DEFAULT_SPLIT_FILE)
    parser.add_argument("--summary-path", type=Path, default=DEFAULT_SUMMARY_PATH)
    parser.add_argument("--fit-split-name", action="append", default=list(DEFAULT_FIT_SPLIT_NAMES))
    parser.add_argument("--target-job-count", type=int, default=DEFAULT_TARGET_JOB_COUNT)
    parser.add_argument("--near-zero-threshold", type=float, default=1.0)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT_JSON)
    return parser.parse_args()


def _safe_float(value: object) -> float | None:
    if value in ("", None):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _job_entity_side(job_id: str) -> str:
    if "-antibody-" in job_id:
        return "antibody"
    if "-antigen-" in job_id:
        return "antigen"
    return "unknown"


def _coverage_effect_bin(value: float, *, near_zero_threshold: float) -> str:
    if value < 0.0:
        return "negative"
    if value < near_zero_threshold:
        return "near_zero"
    return "positive"


def _stage_rank(row: dict[str, Any]) -> int:
    latest_stage_state = str(row.get("latest_stage_state", "not_started")).strip()
    return _ACTIVE_STATE_RANK.get(latest_stage_state, 9)


def _distance_rank(value: float, *, effect_bin: str) -> tuple[float, float]:
    if effect_bin == "negative":
        return (value, abs(value))
    if effect_bin == "near_zero":
        return (abs(value), value)
    return (value, abs(value))


def _load_fit_pair_counts(
    summary_path: Path,
    *,
    near_zero_threshold: float,
) -> tuple[dict[tuple[str, str], int], dict[str, int]]:
    bin_counts = {
        (side, effect_bin): 0
        for side in ("antibody", "antigen", "unknown")
        for effect_bin in ("negative", "near_zero", "positive")
    }
    side_counts = {side: 0 for side in ("antibody", "antigen", "unknown")}
    if not summary_path.exists():
        return bin_counts, side_counts

    try:
        summary = read_json(summary_path)
    except (OSError, json.JSONDecodeError):
        return bin_counts, side_counts
    reports_dir = Path(str(summary.get("reports_dir", "")).strip())
    fit_pairs_path = reports_dir / "fit_pairs.csv"
    if not fit_pairs_path.exists():
        return bin_counts, side_counts

    for row in read_csv_rows(fit_pairs_path):
        job_id = str(row.get("job_id", "")).strip()
        entity_side = _job_entity_side(job_id)
        experimental_ddg = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental_ddg is None:
            continue
        effect_bin = _coverage_effect_bin(experimental_ddg, near_zero_threshold=near_zero_threshold)
        bin_counts[(entity_side, effect_bin)] = bin_counts.get((entity_side, effect_bin), 0) + 1
        side_counts[entity_side] = side_counts.get(entity_side, 0) + 1
    return bin_counts, side_counts


def _candidate_job_rows(
    plan_root: Path,
    split_path: Path,
    *,
    fit_split_names: list[str],
    near_zero_threshold: float,
) -> list[dict[str, Any]]:
    candidates: dict[str, dict[str, Any]] = {}
    for split_name in fit_split_names:
        bundle = report_ab_bind_plan(plan_root, split_name=split_name, split_path=split_path)
        rows = read_csv_rows(Path(bundle["reports_dir"]) / "plan_jobs.csv")
        for row in rows:
            job_id = str(row.get("job_id", "")).strip()
            if not job_id or row.get("ddg_ready") == "True":
                continue
            latest_stage_state = str(row.get("latest_stage_state", "")).strip()
            if latest_stage_state in {"blocked_input", "failed"}:
                continue
            experimental_ddg = _safe_float(row.get("experimental_ddg_kcal_mol"))
            if experimental_ddg is None:
                continue
            entity_side = _job_entity_side(job_id)
            effect_bin = _coverage_effect_bin(experimental_ddg, near_zero_threshold=near_zero_threshold)
            candidate = dict(row)
            candidate["fit_split_name"] = split_name
            candidate["entity_side"] = entity_side
            candidate["effect_bin"] = effect_bin
            candidate["experimental_ddg_value"] = experimental_ddg
            candidates[job_id] = candidate
    return list(candidates.values())


def _fit_pair_count(current_side_counts: dict[str, int]) -> int:
    return sum(int(current_side_counts.get(side, 0)) for side in ("antibody", "antigen", "unknown"))


def _fit_pair_shortfall(current_bin_counts: dict[tuple[str, str], int]) -> int:
    tracked_bins = [
        (side, effect_bin)
        for side in ("antibody", "antigen")
        for effect_bin in ("negative", "near_zero", "positive")
    ]
    return sum(1 for key in tracked_bins if int(current_bin_counts.get(key, 0)) <= 0)


def _selection_priority_components(
    row: dict[str, Any],
    *,
    rolling_bin_counts: dict[tuple[str, str], int],
    rolling_side_counts: dict[str, int],
) -> dict[str, Any]:
    entity_side = str(row["entity_side"])
    effect_bin = str(row["effect_bin"])
    distance_primary, distance_secondary = _distance_rank(
        float(row["experimental_ddg_value"]),
        effect_bin=effect_bin,
    )
    fit_split_name = str(row.get("fit_split_name", "")).strip()
    fit_split_rank = _FIT_SPLIT_RANK.get(fit_split_name, 9)
    return {
        "stage_rank": _stage_rank(row),
        "fit_bin_count": rolling_bin_counts.get((entity_side, effect_bin), 0),
        "fit_side_count": rolling_side_counts.get(entity_side, 0),
        "fit_split_rank": fit_split_rank,
        "distance_primary": distance_primary,
        "distance_secondary": distance_secondary,
    }


def _select_watchlist_job_rows(
    candidates: list[dict[str, Any]],
    *,
    current_bin_counts: dict[tuple[str, str], int],
    current_side_counts: dict[str, int],
    target_job_count: int,
) -> list[dict[str, Any]]:
    remaining = [dict(row) for row in candidates]
    selected: list[dict[str, Any]] = []
    rolling_bin_counts = dict(current_bin_counts)
    rolling_side_counts = dict(current_side_counts)

    while remaining and len(selected) < target_job_count:
        priority_cache = {
            str(row.get("job_id", "")): _selection_priority_components(
                row,
                rolling_bin_counts=rolling_bin_counts,
                rolling_side_counts=rolling_side_counts,
            )
            for row in remaining
        }
        remaining.sort(
            key=lambda row: (
                priority_cache[str(row.get("job_id", ""))]["stage_rank"],
                priority_cache[str(row.get("job_id", ""))]["fit_bin_count"],
                priority_cache[str(row.get("job_id", ""))]["fit_side_count"],
                priority_cache[str(row.get("job_id", ""))]["fit_split_rank"],
                priority_cache[str(row.get("job_id", ""))]["distance_primary"],
                priority_cache[str(row.get("job_id", ""))]["distance_secondary"],
                str(row.get("job_id", "")),
            )
        )
        chosen = remaining.pop(0)
        chosen["selection_round"] = len(selected) + 1
        chosen["selection_priority"] = priority_cache.get(str(chosen.get("job_id", "")), {})
        chosen["priority_score"] = [
            chosen["selection_priority"].get("stage_rank"),
            chosen["selection_priority"].get("fit_bin_count"),
            chosen["selection_priority"].get("fit_side_count"),
            chosen["selection_priority"].get("fit_split_rank"),
            chosen["selection_priority"].get("distance_primary"),
            chosen["selection_priority"].get("distance_secondary"),
        ]
        selected.append(chosen)
        entity_side = str(chosen["entity_side"])
        effect_bin = str(chosen["effect_bin"])
        rolling_bin_counts[(entity_side, effect_bin)] = rolling_bin_counts.get((entity_side, effect_bin), 0) + 1
        rolling_side_counts[entity_side] = rolling_side_counts.get(entity_side, 0) + 1
    return selected


def _build_watchlist_payload(
    *,
    plan_root: Path,
    summary_path: Path,
    fit_split_names: list[str],
    near_zero_threshold: float,
    current_bin_counts: dict[tuple[str, str], int],
    current_side_counts: dict[str, int],
    target_job_count: int,
    candidates: list[dict[str, Any]],
    selected_rows: list[dict[str, Any]],
    selected_job_ids: list[str],
) -> dict[str, Any]:
    return {
        "plan_root": str(plan_root),
        "summary_path": str(summary_path),
        "fit_split_names": fit_split_names,
        "near_zero_threshold": near_zero_threshold,
        "target_job_count": target_job_count,
        "fit_pair_count": _fit_pair_count(current_side_counts),
        "fit_pair_shortfall": _fit_pair_shortfall(current_bin_counts),
        "current_fit_bin_counts": {
            f"{side}:{effect_bin}": count for (side, effect_bin), count in sorted(current_bin_counts.items())
        },
        "current_fit_side_counts": current_side_counts,
        "candidate_count": len(candidates),
        "selected_job_count": len(selected_rows),
        "selected_ready_job_count": sum(1 for row in selected_rows if str(row.get("ddg_ready", "")) == "True"),
        "selected_job_ids": selected_job_ids,
        "selected_rows": [
            {
                "job_id": row.get("job_id"),
                "complex_id": row.get("complex_id"),
                "fit_split_name": row.get("fit_split_name"),
                "fit_entity_side": row.get("entity_side"),
                "fit_effect_bin": row.get("effect_bin"),
                "experimental_ddg_kcal_mol": row.get("experimental_ddg_value"),
                "latest_stage": row.get("latest_stage"),
                "latest_stage_state": row.get("latest_stage_state"),
                "selection_round": row.get("selection_round"),
                "priority_score": row.get("priority_score"),
                "selection_priority": row.get("selection_priority"),
            }
            for row in selected_rows
        ],
    }


def main() -> int:
    args = parse_args()
    plan_root = args.plan_root.expanduser().resolve()
    split_path = args.split_file.expanduser().resolve()
    summary_path = args.summary_path.expanduser().resolve()
    output_json = args.output_json.expanduser().resolve() if args.output_json else None

    fit_split_names = [str(item).strip() for item in args.fit_split_name if str(item).strip()]
    if not fit_split_names:
        fit_split_names = list(DEFAULT_FIT_SPLIT_NAMES)

    current_bin_counts, current_side_counts = _load_fit_pair_counts(
        summary_path,
        near_zero_threshold=args.near_zero_threshold,
    )
    candidates = _candidate_job_rows(
        plan_root,
        split_path,
        fit_split_names=fit_split_names,
        near_zero_threshold=args.near_zero_threshold,
    )
    selected_rows = _select_watchlist_job_rows(
        candidates,
        current_bin_counts=current_bin_counts,
        current_side_counts=current_side_counts,
        target_job_count=max(int(args.target_job_count), 1),
    )
    job_ids = [str(row.get("job_id", "")).strip() for row in selected_rows if str(row.get("job_id", "")).strip()]
    if not job_ids:
        job_ids = list(DEFAULT_FALLBACK_JOB_IDS)

    if output_json is not None:
        write_json(
            output_json,
            _build_watchlist_payload(
                plan_root=plan_root,
                summary_path=summary_path,
                fit_split_names=fit_split_names,
                near_zero_threshold=float(args.near_zero_threshold),
                current_bin_counts=current_bin_counts,
                current_side_counts=current_side_counts,
                target_job_count=max(int(args.target_job_count), 1),
                candidates=candidates,
                selected_rows=selected_rows,
                selected_job_ids=job_ids,
            ),
        )

    for job_id in job_ids:
        print(job_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
