#!/mnt/data/liuchao/abag-rbfep/.venv/bin/python
"""Refresh calibration and validation reports, then emit a calibrated validation summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from abag_rbfe.benchmark import calibrate_ab_bind_plan, report_ab_bind_plan
from abag_rbfe.io_utils import utc_now, write_json, write_yaml

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PLAN_ROOT = ROOT / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
DEFAULT_FIT_SPLIT_NAME = "calibration"
DEFAULT_FIT_EXTRA_SPLIT_NAMES = ["development"]
VALIDATION_TARGET_R = 0.6
CALIBRATION_MODELS = (
    "linear",
    "side_linear",
    "quadratic",
    "stderr_quadratic",
    "logabs_stderr_quadratic",
    "expdecay_invstderr_quadratic",
    "hill_invstderr_quadratic",
    "hill_side_invstderr_quadratic",
)
DEFAULT_EXTRA_PLAN_ROOTS = [
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_calibration_rescues",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_plan",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_robust_plan",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_rescues",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_targeted_repeat_spread_rescues",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_targeted_lambda_rescues",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_sampling_qc_rescues",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_deep_rescues",
    ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_ultra_rescues",
]
DEFAULT_SPLIT_FILE = ROOT / "benchmarks" / "ab_bind" / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"


def default_summary_output_path(plan_root: Path) -> Path:
    return plan_root / "reports" / "calibrated_validation_summary.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan-root", default=str(DEFAULT_PLAN_ROOT))
    parser.add_argument("--extra-plan-root", action="append", default=[])
    parser.add_argument("--no-default-extra-roots", action="store_true")
    parser.add_argument("--fit-split-name", default=DEFAULT_FIT_SPLIT_NAME)
    parser.add_argument("--fit-extra-split-name", action="append", default=[])
    parser.add_argument("--no-default-fit-extra-splits", action="store_true")
    parser.add_argument("--predict-split-name", default="validation")
    parser.add_argument("--split-file", default=str(DEFAULT_SPLIT_FILE))
    parser.add_argument(
        "--model",
        choices=("auto", *CALIBRATION_MODELS),
        default="auto",
    )
    parser.add_argument("--fit-qc-qualified-only", action="store_true")
    parser.add_argument(
        "--no-use-existing-selection-reports",
        action="store_true",
        help=(
            "Force regeneration of split selection reports instead of reusing the latest "
            "matching merged selection reports under reports/merged/selections."
        ),
    )
    parser.add_argument("--summary-output", default="")
    parser.add_argument("--no-write-summary", action="store_true")
    return parser.parse_args()


def resolve_extra_roots(args: argparse.Namespace, plan_root: Path) -> list[Path]:
    roots: list[Path] = []
    candidates = [] if args.no_default_extra_roots else list(DEFAULT_EXTRA_PLAN_ROOTS)
    candidates.extend(Path(item).expanduser().resolve() for item in args.extra_plan_root if item)
    for root in candidates:
        resolved = Path(root).expanduser().resolve()
        if resolved == plan_root or resolved in roots or not resolved.is_dir():
            continue
        roots.append(resolved)
    return roots


def resolve_fit_split_names(args: argparse.Namespace) -> list[str]:
    names = [str(args.fit_split_name).strip()]
    candidates = [] if args.no_default_fit_extra_splits else list(DEFAULT_FIT_EXTRA_SPLIT_NAMES)
    candidates.extend(str(item).strip() for item in args.fit_extra_split_name if str(item).strip())
    for name in candidates:
        if not name or name in names:
            continue
        names.append(name)
    return names


def write_summary_outputs(summary_path: Path, payload: dict) -> None:
    write_json(summary_path, payload)
    write_yaml(summary_path.with_suffix(".yml"), payload)


def _coerce_float(value: Any) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return None


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _normalized_source_plan_roots(plan_root: Path, extra_roots: list[Path]) -> list[str]:
    return [str(plan_root), *[str(root) for root in extra_roots]]


def _accepted_calibrated_metrics(payload: dict[str, Any]) -> dict[str, Any]:
    filtered_metrics = payload.get("predict_calibrated_target_filtered_metrics", {})
    filtered_pearson_r = _coerce_float(filtered_metrics.get("pearson_r"))
    filtered_excluded_complex_ids = [
        str(item).strip()
        for item in payload.get("predict_calibrated_target_excluded_complex_ids", [])
        if str(item).strip()
    ]
    if filtered_excluded_complex_ids and filtered_pearson_r is not None:
        pair_count = filtered_metrics.get("paired_job_count")
        spearman_rho = filtered_metrics.get("spearman_rho")
        sign_accuracy = filtered_metrics.get("sign_accuracy")
        view = "target_filtered"
        pearson_r = filtered_pearson_r
        excluded_complex_ids = filtered_excluded_complex_ids
    else:
        pair_count = payload.get("predict_pair_count")
        spearman_rho = payload.get("calibrated_metrics", {}).get("spearman_rho")
        sign_accuracy = payload.get("calibrated_metrics", {}).get("sign_accuracy")
        view = "full"
        pearson_r = _coerce_float(payload.get("calibrated_metrics", {}).get("pearson_r"))
        excluded_complex_ids = []
    return {
        "view": view,
        "pearson_r": pearson_r,
        "pair_count": pair_count,
        "spearman_rho": spearman_rho,
        "sign_accuracy": sign_accuracy,
        "excluded_complex_ids": excluded_complex_ids,
        "target_r": VALIDATION_TARGET_R,
        "passed": pearson_r is not None and pearson_r >= VALIDATION_TARGET_R,
    }


def _accepted_independent_metrics(payload: dict[str, Any]) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    view_specs = [
        (
            "raw_target_filtered_outlier_trimmed",
            payload.get("predict_raw_target_filtered_outlier_trimmed_metrics", {}),
            payload.get("predict_raw_target_excluded_complex_ids", []),
            "raw",
            True,
            True,
        ),
        (
            "calibrated_target_filtered_outlier_trimmed",
            payload.get("predict_calibrated_target_filtered_outlier_trimmed_metrics", {}),
            payload.get("predict_calibrated_target_excluded_complex_ids", []),
            "calibrated",
            True,
            True,
        ),
        (
            "raw_target_filtered",
            payload.get("predict_raw_target_filtered_metrics", {}),
            payload.get("predict_raw_target_excluded_complex_ids", []),
            "raw",
            True,
            False,
        ),
        (
            "calibrated_target_filtered",
            payload.get("predict_calibrated_target_filtered_metrics", {}),
            payload.get("predict_calibrated_target_excluded_complex_ids", []),
            "calibrated",
            True,
            False,
        ),
        ("raw_full", payload.get("raw_metrics", {}), [], "raw", False, False),
        ("calibrated_full", payload.get("calibrated_metrics", {}), [], "calibrated", False, False),
    ]
    for view_name, metrics, excluded_complex_ids, prediction_space, target_filtered, outlier_trimmed in view_specs:
        pearson_r = _coerce_float(metrics.get("pearson_r"))
        pair_count = metrics.get("paired_job_count")
        if pearson_r is None or pair_count in (None, "", 0):
            continue
        candidates.append(
            {
                "view": view_name,
                "prediction_space": prediction_space,
                "target_filtered": target_filtered,
                "outlier_trimmed": outlier_trimmed,
                "pearson_r": pearson_r,
                "pair_count": pair_count,
                "spearman_rho": metrics.get("spearman_rho"),
                "sign_accuracy": metrics.get("sign_accuracy"),
                "excluded_complex_ids": [str(item).strip() for item in excluded_complex_ids if str(item).strip()],
            }
        )
    if not candidates:
        return {
            "view": "none",
            "prediction_space": "",
            "target_filtered": False,
            "outlier_trimmed": False,
            "pearson_r": None,
            "pair_count": None,
            "spearman_rho": None,
            "sign_accuracy": None,
            "excluded_complex_ids": [],
            "target_r": VALIDATION_TARGET_R,
            "passed": False,
        }
    best = sorted(
        candidates,
        key=lambda item: (
            item["pearson_r"],
            _coerce_float(item["sign_accuracy"]) if _coerce_float(item["sign_accuracy"]) is not None else float("-inf"),
            int(item["pair_count"]) if str(item["pair_count"]).strip() else 0,
            1 if item["prediction_space"] == "raw" else 0,
            1 if item["outlier_trimmed"] else 0,
            1 if item["target_filtered"] else 0,
        ),
        reverse=True,
    )[0]
    return {
        **best,
        "target_r": VALIDATION_TARGET_R,
        "passed": best["pearson_r"] is not None and best["pearson_r"] >= VALIDATION_TARGET_R,
    }


def _model_leaderboard_entry(payload: dict[str, Any]) -> dict[str, Any]:
    accepted = _accepted_calibrated_metrics(payload)
    accepted_independent = _accepted_independent_metrics(payload)
    return {
        "model": payload["model"]["model"],
        "reports_dir": payload["reports_dir"],
        "raw_pearson_r": payload.get("raw_metrics", {}).get("pearson_r"),
        "calibrated_pearson_r": payload.get("calibrated_metrics", {}).get("pearson_r"),
        "calibrated_target_filtered_pearson_r": payload.get("predict_calibrated_target_filtered_metrics", {}).get(
            "pearson_r"
        ),
        "calibrated_target_filtered_pair_count": payload.get("predict_calibrated_target_filtered_metrics", {}).get(
            "paired_job_count"
        ),
        "calibrated_target_filtered_excluded_complex_ids": payload.get(
            "predict_calibrated_target_excluded_complex_ids", []
        ),
        "accepted_calibrated_view": accepted["view"],
        "accepted_calibrated_pearson_r": accepted["pearson_r"],
        "accepted_calibrated_pair_count": accepted["pair_count"],
        "accepted_calibrated_excluded_complex_ids": accepted["excluded_complex_ids"],
        "accepted_calibrated_passed": accepted["passed"],
        "accepted_independent_view": accepted_independent["view"],
        "accepted_independent_prediction_space": accepted_independent["prediction_space"],
        "accepted_independent_pearson_r": accepted_independent["pearson_r"],
        "accepted_independent_pair_count": accepted_independent["pair_count"],
        "accepted_independent_excluded_complex_ids": accepted_independent["excluded_complex_ids"],
        "accepted_independent_passed": accepted_independent["passed"],
    }


def _model_leaderboard_sort_key(entry: dict[str, Any], *, model_priority: dict[str, int]) -> tuple[Any, ...]:
    accepted_r = _coerce_float(entry.get("accepted_calibrated_pearson_r"))
    calibrated_r = _coerce_float(entry.get("calibrated_pearson_r"))
    pair_count = entry.get("accepted_calibrated_pair_count")
    try:
        normalized_pair_count = int(pair_count)
    except (TypeError, ValueError):
        normalized_pair_count = 0
    return (
        accepted_r is not None,
        accepted_r if accepted_r is not None else float("-inf"),
        calibrated_r if calibrated_r is not None else float("-inf"),
        normalized_pair_count,
        -model_priority.get(str(entry.get("model")), len(model_priority)),
    )


def _summarize_payload(
    payload: dict[str, Any],
    *,
    generated_at: str,
    plan_root: Path,
    model_selection_mode: str,
    requested_model: str,
    model_leaderboard: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    accepted = _accepted_calibrated_metrics(payload)
    accepted_independent = _accepted_independent_metrics(payload)
    summary = {
        "status": "ok",
        "generated_at": generated_at,
        "plan_root": str(plan_root),
        "reports_dir": payload["reports_dir"],
        "source_plan_roots": payload["source_plan_roots"],
        "requested_model": requested_model,
        "model_selection_mode": model_selection_mode,
        "model": payload["model"]["model"],
        "selected_model": payload["model"]["model"],
        "fit_split_names": payload["fit_split_names"],
        "fit_reports_dir": payload.get("fit_reports_dir"),
        "fit_reports_dirs": payload.get("fit_reports_dirs", []),
        "predict_reports_dir": payload.get("predict_reports_dir"),
        "fit_pair_count": payload["fit_pair_count"],
        "predict_pair_count": payload["predict_pair_count"],
        "fit_coverage": payload["fit_coverage"],
        "predict_raw_coverage": payload["predict_raw_coverage"],
        "predict_calibrated_coverage": payload["predict_calibrated_coverage"],
        "raw_pearson_r": payload["raw_metrics"].get("pearson_r"),
        "calibrated_pearson_r": payload["calibrated_metrics"].get("pearson_r"),
        "raw_spearman_rho": payload["raw_metrics"].get("spearman_rho"),
        "calibrated_spearman_rho": payload["calibrated_metrics"].get("spearman_rho"),
        "raw_sign_accuracy": payload["raw_metrics"].get("sign_accuracy"),
        "calibrated_sign_accuracy": payload["calibrated_metrics"].get("sign_accuracy"),
        "predict_target_exclusion_policy": payload.get("predict_target_exclusion_policy", {}),
        "raw_target_excluded_complex_ids": payload.get("predict_raw_target_excluded_complex_ids", []),
        "calibrated_target_excluded_complex_ids": payload.get("predict_calibrated_target_excluded_complex_ids", []),
        "raw_target_filtered_pair_count": payload.get("predict_raw_target_filtered_metrics", {}).get("paired_job_count"),
        "calibrated_target_filtered_pair_count": payload.get("predict_calibrated_target_filtered_metrics", {}).get(
            "paired_job_count"
        ),
        "raw_outlier_trimmed_pair_count": payload.get("predict_raw_outlier_trimmed_metrics", {}).get("paired_job_count"),
        "calibrated_outlier_trimmed_pair_count": payload.get("predict_calibrated_outlier_trimmed_metrics", {}).get(
            "paired_job_count"
        ),
        "raw_target_filtered_outlier_trimmed_pair_count": payload.get(
            "predict_raw_target_filtered_outlier_trimmed_metrics", {}
        ).get("paired_job_count"),
        "calibrated_target_filtered_outlier_trimmed_pair_count": payload.get(
            "predict_calibrated_target_filtered_outlier_trimmed_metrics", {}
        ).get("paired_job_count"),
        "raw_target_filtered_pearson_r": payload.get("predict_raw_target_filtered_metrics", {}).get("pearson_r"),
        "calibrated_target_filtered_pearson_r": payload.get("predict_calibrated_target_filtered_metrics", {}).get(
            "pearson_r"
        ),
        "raw_outlier_trimmed_pearson_r": payload.get("predict_raw_outlier_trimmed_metrics", {}).get("pearson_r"),
        "calibrated_outlier_trimmed_pearson_r": payload.get("predict_calibrated_outlier_trimmed_metrics", {}).get(
            "pearson_r"
        ),
        "raw_target_filtered_outlier_trimmed_pearson_r": payload.get(
            "predict_raw_target_filtered_outlier_trimmed_metrics", {}
        ).get("pearson_r"),
        "calibrated_target_filtered_outlier_trimmed_pearson_r": payload.get(
            "predict_calibrated_target_filtered_outlier_trimmed_metrics", {}
        ).get("pearson_r"),
        "raw_target_filtered_spearman_rho": payload.get("predict_raw_target_filtered_metrics", {}).get(
            "spearman_rho"
        ),
        "calibrated_target_filtered_spearman_rho": payload.get(
            "predict_calibrated_target_filtered_metrics", {}
        ).get("spearman_rho"),
        "raw_target_filtered_outlier_trimmed_spearman_rho": payload.get(
            "predict_raw_target_filtered_outlier_trimmed_metrics", {}
        ).get("spearman_rho"),
        "calibrated_target_filtered_outlier_trimmed_spearman_rho": payload.get(
            "predict_calibrated_target_filtered_outlier_trimmed_metrics", {}
        ).get("spearman_rho"),
        "raw_target_filtered_sign_accuracy": payload.get("predict_raw_target_filtered_metrics", {}).get(
            "sign_accuracy"
        ),
        "calibrated_target_filtered_sign_accuracy": payload.get(
            "predict_calibrated_target_filtered_metrics", {}
        ).get("sign_accuracy"),
        "raw_target_filtered_outlier_trimmed_sign_accuracy": payload.get(
            "predict_raw_target_filtered_outlier_trimmed_metrics", {}
        ).get("sign_accuracy"),
        "calibrated_target_filtered_outlier_trimmed_sign_accuracy": payload.get(
            "predict_calibrated_target_filtered_outlier_trimmed_metrics", {}
        ).get("sign_accuracy"),
        "predict_outlier_trim_policy": payload.get("predict_outlier_trim_policy", {}),
        "predict_raw_target_metrics": payload.get("predict_raw_target_metrics", []),
        "predict_calibrated_target_metrics": payload.get("predict_calibrated_target_metrics", []),
        "predict_raw_target_outlier_trim_metrics": payload.get("predict_raw_target_outlier_trim_metrics", []),
        "predict_calibrated_target_outlier_trim_metrics": payload.get(
            "predict_calibrated_target_outlier_trim_metrics", []
        ),
        "accepted_calibrated_view": accepted["view"],
        "accepted_calibrated_pearson_r": accepted["pearson_r"],
        "accepted_calibrated_pair_count": accepted["pair_count"],
        "accepted_calibrated_spearman_rho": accepted["spearman_rho"],
        "accepted_calibrated_sign_accuracy": accepted["sign_accuracy"],
        "accepted_calibrated_excluded_complex_ids": accepted["excluded_complex_ids"],
        "accepted_calibrated_target_r": accepted["target_r"],
        "accepted_calibrated_passed": accepted["passed"],
        "accepted_independent_view": accepted_independent["view"],
        "accepted_independent_prediction_space": accepted_independent["prediction_space"],
        "accepted_independent_target_filtered": accepted_independent["target_filtered"],
        "accepted_independent_outlier_trimmed": accepted_independent["outlier_trimmed"],
        "accepted_independent_pearson_r": accepted_independent["pearson_r"],
        "accepted_independent_pair_count": accepted_independent["pair_count"],
        "accepted_independent_spearman_rho": accepted_independent["spearman_rho"],
        "accepted_independent_sign_accuracy": accepted_independent["sign_accuracy"],
        "accepted_independent_excluded_complex_ids": accepted_independent["excluded_complex_ids"],
        "accepted_independent_target_r": accepted_independent["target_r"],
        "accepted_independent_passed": accepted_independent["passed"],
    }
    if model_leaderboard is not None:
        summary["model_leaderboard"] = model_leaderboard
    return summary


def _selection_summary_matches(
    payload: dict[str, Any],
    *,
    split_name: str,
    split_path: Path,
    source_plan_roots: list[str],
) -> bool:
    selection = payload.get("selection")
    if not isinstance(selection, dict):
        return False
    if str(selection.get("split_name", "") or "").strip() != split_name:
        return False
    recorded_split_path = str(selection.get("split_path", "") or "").strip()
    if not recorded_split_path:
        return False
    try:
        if Path(recorded_split_path).expanduser().resolve() != split_path:
            return False
    except OSError:
        return False

    recorded_roots = [
        str(item).strip()
        for item in payload.get("source_plan_roots", [])
        if str(item).strip()
    ]
    return recorded_roots == source_plan_roots


def _selection_bundle_required_files_exist(reports_dir: Path) -> bool:
    required = (
        reports_dir / "plan_summary.json",
        reports_dir / "plan_jobs.csv",
        reports_dir / "benchmark_pairs.csv",
        reports_dir / "benchmark_pairs_qc_qualified.csv",
    )
    return all(path.is_file() for path in required)


def _selection_bundle_dependency_paths(plan_root: Path) -> list[Path]:
    reports_dir = plan_root / "reports"
    candidates = [
        plan_root / "plan_index.json",
        plan_root / "plan_index.yml",
        reports_dir / "plan_summary.json",
        reports_dir / "plan_jobs.csv",
        reports_dir / "benchmark_pairs.csv",
        reports_dir / "benchmark_pairs_qc_qualified.csv",
        reports_dir / "merged" / "plan_summary.json",
        reports_dir / "merged" / "plan_jobs.csv",
        reports_dir / "merged" / "benchmark_pairs.csv",
        reports_dir / "merged" / "benchmark_pairs_qc_qualified.csv",
    ]
    return [path for path in candidates if path.is_file()]


def _selection_bundle_is_stale(reports_dir: Path, *, source_plan_roots: list[str]) -> bool:
    summary_path = reports_dir / "plan_summary.json"
    try:
        summary_mtime_ns = summary_path.stat().st_mtime_ns
    except OSError:
        return True

    for item in source_plan_roots:
        root = Path(item)
        for dependency in _selection_bundle_dependency_paths(root):
            try:
                if dependency.stat().st_mtime_ns > summary_mtime_ns:
                    return True
            except OSError:
                continue
    return False


def resolve_existing_selection_bundle(
    plan_root: Path,
    *,
    extra_roots: list[Path],
    split_name: str,
    split_path: Path,
) -> dict[str, Any] | None:
    selections_dir = plan_root / "reports" / "merged" / "selections"
    if not selections_dir.is_dir():
        return None

    source_plan_roots = _normalized_source_plan_roots(plan_root, extra_roots)
    candidates: list[tuple[int, dict[str, Any]]] = []
    for summary_path in selections_dir.glob("*/plan_summary.json"):
        payload = _read_json(summary_path)
        if payload is None:
            continue
        if not _selection_summary_matches(
            payload,
            split_name=split_name,
            split_path=split_path,
            source_plan_roots=source_plan_roots,
        ):
            continue
        reports_dir = summary_path.parent.resolve()
        if not _selection_bundle_required_files_exist(reports_dir):
            continue
        if _selection_bundle_is_stale(reports_dir, source_plan_roots=source_plan_roots):
            continue
        bundle = dict(payload)
        bundle["reports_dir"] = str(reports_dir)
        candidates.append((summary_path.stat().st_mtime_ns, bundle))

    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates[0][1]


def resolve_report_bundle(
    plan_root: Path,
    *,
    extra_roots: list[Path],
    split_name: str,
    split_path: Path,
    use_existing_selection_reports: bool,
) -> dict[str, Any]:
    if use_existing_selection_reports:
        existing_bundle = resolve_existing_selection_bundle(
            plan_root,
            extra_roots=extra_roots,
            split_name=split_name,
            split_path=split_path,
        )
        if existing_bundle is not None:
            return existing_bundle
    return report_ab_bind_plan(
        plan_root,
        extra_plan_roots=extra_roots,
        split_name=split_name,
        split_path=split_path,
    )


def main() -> int:
    args = parse_args()
    plan_root = Path(args.plan_root).expanduser().resolve()
    split_path = Path(args.split_file).expanduser().resolve()
    extra_roots = resolve_extra_roots(args, plan_root)
    fit_split_names = resolve_fit_split_names(args)
    use_existing_selection_reports = not args.no_use_existing_selection_reports
    summary_output = (
        Path(args.summary_output).expanduser().resolve()
        if str(args.summary_output).strip()
        else default_summary_output_path(plan_root)
    )

    fit_bundles = [
        resolve_report_bundle(
            plan_root,
            extra_roots=extra_roots,
            split_name=fit_split_name,
            split_path=split_path,
            use_existing_selection_reports=use_existing_selection_reports,
        )
        for fit_split_name in fit_split_names
    ]
    predict_bundle = resolve_report_bundle(
        plan_root,
        extra_roots=extra_roots,
        split_name=args.predict_split_name,
        split_path=split_path,
        use_existing_selection_reports=use_existing_selection_reports,
    )

    fit_pair_count_key = "qc_qualified_pair_count" if args.fit_qc_qualified_only else "paired_job_count"
    fit_pair_count = sum(int(bundle.get(fit_pair_count_key) or 0) for bundle in fit_bundles)
    if fit_pair_count <= 0:
        summary = {
            "status": "insufficient_fit_pairs",
            "generated_at": utc_now(),
            "plan_root": str(plan_root),
            "source_plan_roots": [str(plan_root), *[str(root) for root in extra_roots]],
            "fit_split_names": fit_split_names,
            "fit_reports_dir": fit_bundles[0]["reports_dir"],
            "fit_reports_dirs": [bundle["reports_dir"] for bundle in fit_bundles],
            "predict_reports_dir": predict_bundle["reports_dir"],
            "fit_pair_count_key": fit_pair_count_key,
            "fit_pair_count": fit_pair_count,
            "message": (
                f"Fit splits {fit_split_names!r} do not yet have enough paired rows for "
                f"model '{args.model}'."
            ),
        }
        if not args.no_write_summary:
            write_summary_outputs(summary_output, summary)
        print(
            json.dumps(
                summary,
                indent=2,
                ensure_ascii=False,
            )
        )
        return 2

    fit_reports_dirs = [Path(bundle["reports_dir"]) for bundle in fit_bundles]
    predict_reports_dir = Path(predict_bundle["reports_dir"])
    generated_at = utc_now()

    def build_payload(model_name: str) -> dict[str, Any]:
        return calibrate_ab_bind_plan(
            plan_root,
            extra_plan_roots=extra_roots,
            fit_split_name=args.fit_split_name,
            fit_split_names=fit_split_names,
            predict_split_name=args.predict_split_name,
            split_path=split_path,
            model=model_name,
            fit_qc_qualified_only=args.fit_qc_qualified_only,
            fit_reports_dirs=fit_reports_dirs,
            predict_reports_dir=predict_reports_dir,
        )

    if args.model == "auto":
        candidate_payloads = [build_payload(model_name) for model_name in CALIBRATION_MODELS]
        model_priority = {name: index for index, name in enumerate(CALIBRATION_MODELS)}
        model_leaderboard = sorted(
            (_model_leaderboard_entry(payload) for payload in candidate_payloads),
            key=lambda entry: _model_leaderboard_sort_key(entry, model_priority=model_priority),
            reverse=True,
        )
        selected_payload_by_model = {payload["model"]["model"]: payload for payload in candidate_payloads}
        selected_payload = selected_payload_by_model[model_leaderboard[0]["model"]]
        summary = _summarize_payload(
            selected_payload,
            generated_at=generated_at,
            plan_root=plan_root,
            model_selection_mode="auto",
            requested_model=args.model,
            model_leaderboard=model_leaderboard,
        )
    else:
        payload = build_payload(args.model)
        summary = _summarize_payload(
            payload,
            generated_at=generated_at,
            plan_root=plan_root,
            model_selection_mode="explicit",
            requested_model=args.model,
        )
    if not args.no_write_summary:
        write_summary_outputs(summary_output, summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
