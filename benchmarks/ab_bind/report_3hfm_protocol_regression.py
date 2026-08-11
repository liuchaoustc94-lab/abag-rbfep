#!/mnt/data/liuchao/abag-rbfep/.venv/bin/python
"""Refresh and summarize the dedicated 3HFM protocol-regression benchmark view."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

from abag_rbfe.benchmark import report_ab_bind_plan
from abag_rbfe.io_utils import utc_now, write_json, write_yaml

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PLAN_ROOT = ROOT / "runs" / "benchmarks" / "abbind_3hfm_protocol_regression"
DEFAULT_COMPLEX_ID = "3HFM"


def default_summary_output_path(plan_root: Path) -> Path:
    return plan_root / "reports" / "3hfm_protocol_regression_summary.json"


def default_merged_summary_alias_path(plan_root: Path) -> Path:
    return plan_root / "reports" / "3hfm_protocol_regression_merged_summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan-root", default=str(DEFAULT_PLAN_ROOT))
    parser.add_argument("--extra-plan-root", action="append", default=[])
    parser.add_argument("--complex-id", default=DEFAULT_COMPLEX_ID)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--summary-output", default="")
    parser.add_argument("--no-write-summary", action="store_true")
    return parser.parse_args()


def resolve_extra_roots(args: argparse.Namespace, plan_root: Path) -> list[Path]:
    roots: list[Path] = []
    for item in args.extra_plan_root:
        resolved = Path(item).expanduser().resolve()
        if resolved == plan_root or resolved in roots or not resolved.is_dir():
            continue
        roots.append(resolved)
    return roots


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
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _top_abs_error_pairs(rows: list[dict[str, str]], *, limit: int) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    for row in rows:
        abs_error = _safe_float(row.get("abs_error_kcal_mol"))
        if abs_error is None:
            continue
        items.append(
            {
                "job_id": row.get("job_id", ""),
                "mutation_group_id": row.get("mutation_group_id", ""),
                "predicted_ddg_kcal_mol": _safe_float(row.get("predicted_ddg_kcal_mol")),
                "experimental_ddg_kcal_mol": _safe_float(row.get("experimental_ddg_kcal_mol")),
                "ddg_error_kcal_mol": _safe_float(row.get("ddg_error_kcal_mol")),
                "abs_error_kcal_mol": abs_error,
                "ddg_bar_stderr_kcal_mol": _safe_float(row.get("ddg_bar_stderr_kcal_mol")),
                "qc_status": row.get("qc_status", ""),
                "benchmark_qc_qualified": _safe_bool(row.get("benchmark_qc_qualified")),
            }
        )
    return sorted(items, key=lambda item: (-item["abs_error_kcal_mol"], item["job_id"]))[: max(int(limit), 0)]


def write_summary_outputs(summary_path: Path, payload: dict[str, Any]) -> None:
    write_json(summary_path, payload)
    write_yaml(summary_path.with_suffix(".yml"), payload)


def main() -> int:
    args = parse_args()
    plan_root = Path(args.plan_root).expanduser().resolve()
    extra_roots = resolve_extra_roots(args, plan_root)
    complex_id = str(args.complex_id).strip() or DEFAULT_COMPLEX_ID
    summary_output = (
        Path(args.summary_output).expanduser().resolve()
        if str(args.summary_output).strip()
        else default_summary_output_path(plan_root)
    )

    bundle = report_ab_bind_plan(
        plan_root,
        extra_plan_roots=extra_roots,
        complex_ids=[complex_id],
    )
    reports_dir = Path(bundle["reports_dir"]).expanduser().resolve()
    pair_rows = _load_csv_rows(reports_dir / "benchmark_pairs.csv")
    qc_pair_rows = _load_csv_rows(reports_dir / "benchmark_pairs_qc_qualified.csv")

    summary: dict[str, Any] = {
        "status": "ok" if pair_rows else "insufficient_pairs",
        "generated_at": utc_now(),
        "plan_root": str(plan_root),
        "reports_dir": str(reports_dir),
        "source_plan_roots": bundle.get("source_plan_roots", [str(plan_root)]),
        "complex_id": complex_id,
        "selected_job_count": bundle.get("selected_job_count"),
        "ddg_ready_count": bundle.get("ddg_ready_count"),
        "paired_job_count": len(pair_rows),
        "qc_qualified_pair_count": len(qc_pair_rows),
        "resumable_job_count": bundle.get("resumable_job_count"),
        "running_sample_job_count": bundle.get("running_sample_job_count"),
        "running_equilibrate_job_count": bundle.get("running_equilibrate_job_count"),
        "qc_counts": bundle.get("qc_counts", {}),
        "validation_failure_taxonomy": bundle.get("validation_failure_taxonomy", {}),
        "overall_pearson_r": bundle.get("benchmark_metrics", {}).get("pearson_r"),
        "overall_spearman_rho": bundle.get("benchmark_metrics", {}).get("spearman_rho"),
        "overall_sign_accuracy": bundle.get("benchmark_metrics", {}).get("sign_accuracy"),
        "overall_rmse_kcal_mol": bundle.get("benchmark_metrics", {}).get("rmse_kcal_mol"),
        "overall_mae_kcal_mol": bundle.get("benchmark_metrics", {}).get("mae_kcal_mol"),
        "qc_qualified_pearson_r": bundle.get("benchmark_metrics_qc_qualified", {}).get("pearson_r"),
        "qc_qualified_spearman_rho": bundle.get("benchmark_metrics_qc_qualified", {}).get("spearman_rho"),
        "qc_qualified_sign_accuracy": bundle.get("benchmark_metrics_qc_qualified", {}).get("sign_accuracy"),
        "qc_qualified_rmse_kcal_mol": bundle.get("benchmark_metrics_qc_qualified", {}).get("rmse_kcal_mol"),
        "qc_qualified_mae_kcal_mol": bundle.get("benchmark_metrics_qc_qualified", {}).get("mae_kcal_mol"),
        "top_abs_error_pairs": _top_abs_error_pairs(pair_rows, limit=args.top_n),
    }
    if not pair_rows:
        summary["message"] = (
            f"Complex {complex_id} under {plan_root} does not yet have paired benchmark rows for protocol regression."
        )

    if not args.no_write_summary:
        write_summary_outputs(summary_output, summary)
        if extra_roots and summary_output == default_summary_output_path(plan_root):
            write_summary_outputs(default_merged_summary_alias_path(plan_root), summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if pair_rows else 2


if __name__ == "__main__":
    raise SystemExit(main())
