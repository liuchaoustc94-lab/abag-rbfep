#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


RAW_COLOR = "#1f77b4"
CAL_COLOR = "#ff7f0e"
EXCLUDED_COLOR = "#c0392b"
INCLUDED_COLOR = "#2c7a4b"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Render a per-target validation dashboard from the current "
            "calibrated validation summary."
        )
    )
    repo_root = Path(__file__).resolve().parents[2]
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=repo_root
        / "runs/benchmarks/abbind_core_v1_quick_plan/reports/calibrated_validation_summary.json",
        help="Path to calibrated_validation_summary.json",
    )
    parser.add_argument(
        "--pair-csv",
        type=Path,
        default=None,
        help=(
            "Path to predict_pairs_calibrated.csv. Defaults to "
            "<reports_dir>/predict_pairs_calibrated.csv from the summary JSON."
        ),
    )
    parser.add_argument(
        "--patel-summary",
        type=Path,
        default=repo_root
        / "runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference/reports/patel_2021_3hfm_summary.json",
        help="Optional Patel 2021 3HFM external-reference summary JSON.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "docs/validation_target_summary",
        help="Directory for the generated markdown, JSON, and figure files.",
    )
    return parser


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def to_float(value: Any) -> float | None:
    if value in (None, "", "NA", "None"):
        return None
    return float(value)


def format_float(value: float | None, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{value:.{digits}f}"


def format_bool(value: bool) -> str:
    return "yes" if value else "no"


def predicted_span(values: list[float]) -> float:
    if not values:
        return 0.0
    return max(values) - min(values)


def mean(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def percentile(values: list[float], percent: float) -> float | None:
    if not values:
        return None
    sorted_values = sorted(values)
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (len(sorted_values) - 1) * (percent / 100.0)
    lower_index = math.floor(position)
    upper_index = math.ceil(position)
    if lower_index == upper_index:
        return sorted_values[lower_index]
    weight = position - lower_index
    return (
        sorted_values[lower_index] * (1.0 - weight)
        + sorted_values[upper_index] * weight
    )


def rank_values(values: list[float]) -> list[float]:
    indexed_values = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(indexed_values):
        end = start
        while (
            end + 1 < len(indexed_values)
            and indexed_values[end + 1][1] == indexed_values[start][1]
        ):
            end += 1
        average_rank = (start + end) / 2.0 + 1.0
        for rank_index in range(start, end + 1):
            original_index = indexed_values[rank_index][0]
            ranks[original_index] = average_rank
        start = end + 1
    return ranks


def pearson_correlation(x_values: list[float], y_values: list[float]) -> float | None:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        return None
    mean_x = mean(x_values)
    mean_y = mean(y_values)
    if mean_x is None or mean_y is None:
        return None
    centered_pairs = [
        (x_value - mean_x, y_value - mean_y)
        for x_value, y_value in zip(x_values, y_values)
    ]
    numerator = sum(x_value * y_value for x_value, y_value in centered_pairs)
    denominator_x = math.sqrt(sum(x_value * x_value for x_value, _ in centered_pairs))
    denominator_y = math.sqrt(sum(y_value * y_value for _, y_value in centered_pairs))
    if denominator_x == 0.0 or denominator_y == 0.0:
        return None
    return numerator / (denominator_x * denominator_y)


def spearman_correlation(x_values: list[float], y_values: list[float]) -> float | None:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        return None
    return pearson_correlation(rank_values(x_values), rank_values(y_values))


def sign_accuracy(predicted_values: list[float], experimental_values: list[float]) -> float | None:
    if len(predicted_values) != len(experimental_values) or not predicted_values:
        return None
    matches = 0
    for predicted_value, experimental_value in zip(predicted_values, experimental_values):
        predicted_sign = 0 if abs(predicted_value) < 1e-12 else (1 if predicted_value > 0 else -1)
        experimental_sign = (
            0 if abs(experimental_value) < 1e-12 else (1 if experimental_value > 0 else -1)
        )
        if predicted_sign == experimental_sign:
            matches += 1
    return matches / len(predicted_values)


def compute_pair_metrics(pair_rows: list[dict[str, Any]], prediction_field: str) -> dict[str, Any]:
    predicted_values: list[float] = []
    experimental_values: list[float] = []
    abs_errors: list[float] = []
    for row in pair_rows:
        predicted_value = to_float(row[prediction_field])
        experimental_value = to_float(row["experimental_ddg_kcal_mol"])
        if predicted_value is None or experimental_value is None:
            continue
        predicted_values.append(predicted_value)
        experimental_values.append(experimental_value)
        abs_errors.append(abs(predicted_value - experimental_value))

    return {
        "paired_job_count": len(predicted_values),
        "pearson_r": pearson_correlation(predicted_values, experimental_values),
        "spearman_rho": spearman_correlation(predicted_values, experimental_values),
        "mae_kcal_mol": mean(abs_errors),
        "sign_accuracy": sign_accuracy(predicted_values, experimental_values),
        "max_abs_error_kcal_mol": max(abs_errors) if abs_errors else None,
    }


def build_tukey_outlier_analysis(
    pair_rows: list[dict[str, Any]], prediction_field: str
) -> dict[str, Any]:
    scored_rows: list[dict[str, Any]] = []
    for row in pair_rows:
        predicted_value = to_float(row[prediction_field])
        experimental_value = to_float(row["experimental_ddg_kcal_mol"])
        if predicted_value is None or experimental_value is None:
            continue
        scored_rows.append(
            {
                "job_id": row["job_id"],
                "pair_row": row,
                "abs_error_kcal_mol": abs(predicted_value - experimental_value),
            }
        )

    abs_errors = [row["abs_error_kcal_mol"] for row in scored_rows]
    q1 = percentile(abs_errors, 25.0)
    q3 = percentile(abs_errors, 75.0)
    iqr = None if q1 is None or q3 is None else q3 - q1
    threshold = None if q3 is None or iqr is None else q3 + 1.5 * iqr

    removed_rows = []
    kept_rows = []
    if threshold is None:
        kept_rows = [row["pair_row"] for row in scored_rows]
    else:
        for row in scored_rows:
            if row["abs_error_kcal_mol"] > threshold:
                removed_rows.append(row)
            else:
                kept_rows.append(row["pair_row"])

    return {
        "q1_abs_error_kcal_mol": q1,
        "q3_abs_error_kcal_mol": q3,
        "iqr_abs_error_kcal_mol": iqr,
        "threshold_abs_error_kcal_mol": threshold,
        "removed_pair_count": len(removed_rows),
        "removed_job_ids": [row["job_id"] for row in removed_rows],
        "filtered_pair_count": len(kept_rows),
        "filtered_metrics": compute_pair_metrics(kept_rows, prediction_field),
    }


def filter_pair_rows_by_job_ids(
    pair_rows: list[dict[str, Any]], removed_job_ids: list[str]
) -> list[dict[str, Any]]:
    removed_job_id_set = set(removed_job_ids)
    return [row for row in pair_rows if row["job_id"] not in removed_job_id_set]


def select_outliers(
    pair_rows: list[dict[str, Any]], prediction_field: str, limit: int = 3
) -> list[dict[str, Any]]:
    scored_rows: list[dict[str, Any]] = []
    for row in pair_rows:
        predicted = to_float(row[prediction_field])
        experimental = to_float(row["experimental_ddg_kcal_mol"])
        if predicted is None or experimental is None:
            continue
        abs_error = abs(predicted - experimental)
        scored_rows.append(
            {
                "job_id": row["job_id"],
                "predicted_ddg_kcal_mol": predicted,
                "experimental_ddg_kcal_mol": experimental,
                "abs_error_kcal_mol": abs_error,
            }
        )
    scored_rows.sort(key=lambda item: item["abs_error_kcal_mol"], reverse=True)
    return scored_rows[:limit]


def infer_outlier_sensitivity_readout(target: dict[str, Any]) -> str:
    raw_metrics = target["raw_metrics"]
    calibrated_metrics = target["calibrated_metrics"]
    raw_tukey = target["raw_tukey_outlier_analysis"]
    calibrated_tukey = target["calibrated_tukey_outlier_analysis"]

    raw_removed = raw_tukey["removed_pair_count"]
    calibrated_removed = calibrated_tukey["removed_pair_count"]
    raw_r_before = to_float(raw_metrics.get("pearson_r"))
    raw_r_after = to_float(raw_tukey["filtered_metrics"].get("pearson_r"))
    cal_r_before = to_float(calibrated_metrics.get("pearson_r"))
    cal_r_after = to_float(calibrated_tukey["filtered_metrics"].get("pearson_r"))

    if (
        raw_removed > 0
        and raw_r_before is not None
        and raw_r_after is not None
        and raw_r_after - raw_r_before >= 0.3
    ):
        return "Raw target behavior is strongly driven by a small number of extreme outliers."
    if (
        calibrated_removed == 0
        and cal_r_before is not None
        and cal_r_before < 0.0
    ):
        return "The calibrated ranking problem persists even after applying a formal outlier screen."
    if raw_removed == 0 and calibrated_removed == 0:
        return "No Tukey outliers are detected; this target behaves like a target-wide effect, not a single-point artifact."
    return "Outlier trimming changes the scale more than the qualitative conclusion for this target."


def infer_target_readout(target: dict[str, Any]) -> str:
    raw = target["raw_metrics"]
    calibrated = target["calibrated_metrics"]
    raw_r = to_float(raw.get("pearson_r"))
    cal_r = to_float(calibrated.get("pearson_r"))
    raw_mae = to_float(raw.get("mae_kcal_mol"))
    cal_mae = to_float(calibrated.get("mae_kcal_mol"))
    compression_ratio = target.get("compression_ratio")
    excluded = bool(calibrated.get("excluded_from_target_filtered_metrics"))
    raw_sign = to_float(raw.get("sign_accuracy"))
    cal_sign = to_float(calibrated.get("sign_accuracy"))

    if (
        raw_r is not None
        and cal_r is not None
        and raw_mae is not None
        and cal_mae is not None
        and compression_ratio is not None
        and compression_ratio < 0.35
        and cal_r < raw_r - 0.2
        and cal_mae < raw_mae
    ):
        return "Calibration compresses the dynamic range: MAE improves, but target ranking degrades."
    if (
        raw_r is not None
        and cal_r is not None
        and raw_r > 0.5
        and cal_r < 0.0
    ):
        return "Calibration flips the target ranking even though the raw target already looked learnable."
    if (
        raw_r is not None
        and cal_r is not None
        and raw_r < 0.15
        and cal_r < 0.15
        and raw_sign is not None
        and cal_sign is not None
        and raw_sign < 0.6
        and cal_sign < 0.9
    ):
        return "Neither raw nor calibrated predictions recover the target-specific ranking."
    if (
        raw_r is not None
        and cal_r is not None
        and raw_mae is not None
        and cal_mae is not None
        and cal_r > raw_r + 0.2
        and cal_mae < raw_mae
    ):
        return "Calibration improves both rank and scale for this target."
    if excluded:
        return "This target is excluded from the accepted target-filtered validation view."
    return "This target remains in the accepted filtered view, but the per-target response is still mixed."


def normalise_target_metrics(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {row["complex_id"]: row for row in rows}


def build_target_payload(
    summary: dict[str, Any],
    pair_rows: list[dict[str, str]],
    patel_summary: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    raw_metrics = normalise_target_metrics(summary.get("predict_raw_target_metrics", []))
    calibrated_metrics = normalise_target_metrics(
        summary.get("predict_calibrated_target_metrics", [])
    )
    targets = [row["complex_id"] for row in summary.get("predict_raw_target_metrics", [])]
    for target in calibrated_metrics:
        if target not in targets:
            targets.append(target)

    rows_by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in pair_rows:
        rows_by_target[row["complex_id"]].append(row)

    payload: list[dict[str, Any]] = []
    for target in targets:
        target_pair_rows = rows_by_target.get(target, [])
        raw_predictions = [
            to_float(row["raw_ddg_kcal_mol"])
            for row in target_pair_rows
            if to_float(row["raw_ddg_kcal_mol"]) is not None
        ]
        calibrated_predictions = [
            to_float(row["predicted_ddg_kcal_mol"])
            for row in target_pair_rows
            if to_float(row["predicted_ddg_kcal_mol"]) is not None
        ]
        experimental_values = [
            to_float(row["experimental_ddg_kcal_mol"])
            for row in target_pair_rows
            if to_float(row["experimental_ddg_kcal_mol"]) is not None
        ]
        raw_span = predicted_span([value for value in raw_predictions if value is not None])
        calibrated_span = predicted_span(
            [value for value in calibrated_predictions if value is not None]
        )
        compression_ratio = None
        if raw_span > 0.0:
            compression_ratio = calibrated_span / raw_span

        item = {
            "complex_id": target,
            "pair_count": len(target_pair_rows),
            "raw_metrics": raw_metrics.get(target, {}),
            "calibrated_metrics": calibrated_metrics.get(target, {}),
            "raw_prediction_span_kcal_mol": raw_span,
            "calibrated_prediction_span_kcal_mol": calibrated_span,
            "compression_ratio": compression_ratio,
            "experimental_span_kcal_mol": predicted_span(
                [value for value in experimental_values if value is not None]
            ),
            "raw_outliers": select_outliers(target_pair_rows, "raw_ddg_kcal_mol"),
            "calibrated_outliers": select_outliers(
                target_pair_rows, "predicted_ddg_kcal_mol"
            ),
            "raw_tukey_outlier_analysis": build_tukey_outlier_analysis(
                target_pair_rows, "raw_ddg_kcal_mol"
            ),
            "calibrated_tukey_outlier_analysis": build_tukey_outlier_analysis(
                target_pair_rows, "predicted_ddg_kcal_mol"
            ),
            "pair_rows": target_pair_rows,
        }
        item["readout"] = infer_target_readout(item)
        item["outlier_sensitivity_readout"] = infer_outlier_sensitivity_readout(item)
        if target == "3HFM" and patel_summary is not None:
            item["patel_external_summary"] = {
                "status": patel_summary.get("status"),
                "paired_job_count": patel_summary.get("paired_job_count"),
                "incomplete_job_count": patel_summary.get("incomplete_job_count"),
                "message": patel_summary.get("message"),
            }
        payload.append(item)
    return payload


def render_metrics_figure(
    targets: list[dict[str, Any]], excluded_targets: set[str], output_path: Path
) -> None:
    labels = []
    raw_pearson = []
    cal_pearson = []
    raw_sign = []
    cal_sign = []
    raw_mae = []
    cal_mae = []
    raw_span = []
    cal_span = []

    for target in targets:
        target_id = target["complex_id"]
        pair_count = target["pair_count"]
        star = "*" if target_id in excluded_targets else ""
        labels.append(f"{target_id}{star}\n(n={pair_count})")
        raw = target["raw_metrics"]
        calibrated = target["calibrated_metrics"]
        raw_pearson.append(to_float(raw.get("pearson_r")) or 0.0)
        cal_pearson.append(to_float(calibrated.get("pearson_r")) or 0.0)
        raw_sign.append(to_float(raw.get("sign_accuracy")) or 0.0)
        cal_sign.append(to_float(calibrated.get("sign_accuracy")) or 0.0)
        raw_mae.append(to_float(raw.get("mae_kcal_mol")) or 0.0)
        cal_mae.append(to_float(calibrated.get("mae_kcal_mol")) or 0.0)
        raw_span.append(target["raw_prediction_span_kcal_mol"])
        cal_span.append(target["calibrated_prediction_span_kcal_mol"])

    x = list(range(len(targets)))
    width = 0.38

    fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    axes = axes.ravel()

    axes[0].bar([i - width / 2 for i in x], raw_pearson, width=width, color=RAW_COLOR, label="raw")
    axes[0].bar(
        [i + width / 2 for i in x],
        cal_pearson,
        width=width,
        color=CAL_COLOR,
        label="calibrated",
    )
    axes[0].axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    axes[0].set_title("Per-target Pearson r")
    axes[0].set_xticks(x, labels)
    axes[0].tick_params(axis="x", labelsize=9)
    axes[0].legend(frameon=False)

    axes[1].bar([i - width / 2 for i in x], raw_sign, width=width, color=RAW_COLOR)
    axes[1].bar([i + width / 2 for i in x], cal_sign, width=width, color=CAL_COLOR)
    axes[1].set_ylim(0.0, 1.05)
    axes[1].set_title("Per-target sign accuracy")
    axes[1].set_xticks(x, labels)
    axes[1].tick_params(axis="x", labelsize=9)

    axes[2].bar([i - width / 2 for i in x], raw_mae, width=width, color=RAW_COLOR)
    axes[2].bar([i + width / 2 for i in x], cal_mae, width=width, color=CAL_COLOR)
    axes[2].set_title("Per-target MAE (kcal/mol)")
    axes[2].set_xticks(x, labels)
    axes[2].tick_params(axis="x", labelsize=9)

    axes[3].bar([i - width / 2 for i in x], raw_span, width=width, color=RAW_COLOR)
    axes[3].bar([i + width / 2 for i in x], cal_span, width=width, color=CAL_COLOR)
    axes[3].set_title("Prediction span by target (kcal/mol)")
    axes[3].set_xticks(x, labels)
    axes[3].tick_params(axis="x", labelsize=9)

    for axis in axes:
        axis.grid(axis="y", alpha=0.25, linestyle="--")

    fig.suptitle(
        "Validation targets: raw vs calibrated target behavior\n"
        "* excluded from accepted target-filtered validation view",
        fontsize=15,
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def render_scatter_figure(
    targets: list[dict[str, Any]], excluded_targets: set[str], output_path: Path
) -> None:
    all_experimental = []
    all_predictions = []
    for target in targets:
        for row in target["pair_rows"]:
            experimental = to_float(row["experimental_ddg_kcal_mol"])
            raw_pred = to_float(row["raw_ddg_kcal_mol"])
            calibrated_pred = to_float(row["predicted_ddg_kcal_mol"])
            if experimental is not None:
                all_experimental.append(experimental)
            if raw_pred is not None:
                all_predictions.append(raw_pred)
            if calibrated_pred is not None:
                all_predictions.append(calibrated_pred)

    x_min = min(all_experimental) - 0.5
    x_max = max(all_experimental) + 0.5
    y_min = min(all_predictions + all_experimental) - 1.0
    y_max = max(all_predictions + all_experimental) + 1.0

    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes = axes.ravel()

    for axis, target in zip(axes, targets):
        target_id = target["complex_id"]
        excluded = target_id in excluded_targets
        raw = target["raw_metrics"]
        calibrated = target["calibrated_metrics"]
        experimental = [
            to_float(row["experimental_ddg_kcal_mol"])
            for row in target["pair_rows"]
            if to_float(row["experimental_ddg_kcal_mol"]) is not None
        ]
        raw_pred = [
            to_float(row["raw_ddg_kcal_mol"])
            for row in target["pair_rows"]
            if to_float(row["raw_ddg_kcal_mol"]) is not None
        ]
        cal_pred = [
            to_float(row["predicted_ddg_kcal_mol"])
            for row in target["pair_rows"]
            if to_float(row["predicted_ddg_kcal_mol"]) is not None
        ]

        axis.scatter(experimental, raw_pred, color=RAW_COLOR, alpha=0.85, s=35, label="raw")
        axis.scatter(
            experimental,
            cal_pred,
            color=CAL_COLOR,
            marker="^",
            alpha=0.85,
            s=38,
            label="calibrated",
        )
        axis.plot([x_min, x_max], [x_min, x_max], linestyle="--", color="black", alpha=0.35)
        axis.set_xlim(x_min, x_max)
        axis.set_ylim(y_min, y_max)
        axis.grid(alpha=0.2, linestyle="--")
        title_color = EXCLUDED_COLOR if excluded else INCLUDED_COLOR
        axis.set_title(
            f"{target_id}  r: {format_float(to_float(raw.get('pearson_r')), 2)}"
            f" -> {format_float(to_float(calibrated.get('pearson_r')), 2)}",
            color=title_color,
            fontsize=11,
        )
        axis.set_xlabel("experimental ddG")
        axis.set_ylabel("predicted ddG")

    fig.suptitle(
        "Validation targets: experimental vs predicted ddG (blue circles = raw, orange triangles = calibrated)",
        fontsize=14,
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def render_trimmed_scatter_figure(
    targets: list[dict[str, Any]], excluded_targets: set[str], output_path: Path
) -> None:
    all_experimental = []
    all_predictions = []
    for target in targets:
        raw_rows = filter_pair_rows_by_job_ids(
            target["pair_rows"], target["raw_tukey_outlier_analysis"]["removed_job_ids"]
        )
        calibrated_rows = filter_pair_rows_by_job_ids(
            target["pair_rows"],
            target["calibrated_tukey_outlier_analysis"]["removed_job_ids"],
        )
        for row in raw_rows + calibrated_rows:
            experimental = to_float(row["experimental_ddg_kcal_mol"])
            raw_pred = to_float(row["raw_ddg_kcal_mol"])
            calibrated_pred = to_float(row["predicted_ddg_kcal_mol"])
            if experimental is not None:
                all_experimental.append(experimental)
            if raw_pred is not None:
                all_predictions.append(raw_pred)
            if calibrated_pred is not None:
                all_predictions.append(calibrated_pred)

    x_min = min(all_experimental) - 0.5
    x_max = max(all_experimental) + 0.5
    y_min = min(all_predictions + all_experimental) - 1.0
    y_max = max(all_predictions + all_experimental) + 1.0

    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    axes = axes.ravel()

    for axis, target in zip(axes, targets):
        target_id = target["complex_id"]
        excluded = target_id in excluded_targets
        raw_metrics = target["raw_metrics"]
        calibrated_metrics = target["calibrated_metrics"]
        raw_trimmed_metrics = target["raw_tukey_outlier_analysis"]["filtered_metrics"]
        calibrated_trimmed_metrics = target["calibrated_tukey_outlier_analysis"][
            "filtered_metrics"
        ]
        raw_rows = filter_pair_rows_by_job_ids(
            target["pair_rows"], target["raw_tukey_outlier_analysis"]["removed_job_ids"]
        )
        calibrated_rows = filter_pair_rows_by_job_ids(
            target["pair_rows"],
            target["calibrated_tukey_outlier_analysis"]["removed_job_ids"],
        )

        raw_experimental = [
            to_float(row["experimental_ddg_kcal_mol"])
            for row in raw_rows
            if to_float(row["experimental_ddg_kcal_mol"]) is not None
            and to_float(row["raw_ddg_kcal_mol"]) is not None
        ]
        raw_predicted = [
            to_float(row["raw_ddg_kcal_mol"])
            for row in raw_rows
            if to_float(row["experimental_ddg_kcal_mol"]) is not None
            and to_float(row["raw_ddg_kcal_mol"]) is not None
        ]
        calibrated_experimental = [
            to_float(row["experimental_ddg_kcal_mol"])
            for row in calibrated_rows
            if to_float(row["experimental_ddg_kcal_mol"]) is not None
            and to_float(row["predicted_ddg_kcal_mol"]) is not None
        ]
        calibrated_predicted = [
            to_float(row["predicted_ddg_kcal_mol"])
            for row in calibrated_rows
            if to_float(row["experimental_ddg_kcal_mol"]) is not None
            and to_float(row["predicted_ddg_kcal_mol"]) is not None
        ]

        axis.scatter(
            raw_experimental,
            raw_predicted,
            color=RAW_COLOR,
            alpha=0.85,
            s=35,
            label="raw trimmed",
        )
        axis.scatter(
            calibrated_experimental,
            calibrated_predicted,
            color=CAL_COLOR,
            marker="^",
            alpha=0.85,
            s=38,
            label="cal trimmed",
        )
        axis.plot([x_min, x_max], [x_min, x_max], linestyle="--", color="black", alpha=0.35)
        axis.set_xlim(x_min, x_max)
        axis.set_ylim(y_min, y_max)
        axis.grid(alpha=0.2, linestyle="--")
        title_color = EXCLUDED_COLOR if excluded else INCLUDED_COLOR
        axis.set_title(
            f"{target_id} raw {format_float(to_float(raw_metrics.get('pearson_r')), 2)}"
            f" -> {format_float(to_float(raw_trimmed_metrics.get('pearson_r')), 2)} | "
            f"cal {format_float(to_float(calibrated_metrics.get('pearson_r')), 2)}"
            f" -> {format_float(to_float(calibrated_trimmed_metrics.get('pearson_r')), 2)}\n"
            f"removed raw {target['raw_tukey_outlier_analysis']['removed_pair_count']}, "
            f"cal {target['calibrated_tukey_outlier_analysis']['removed_pair_count']}",
            color=title_color,
            fontsize=10,
        )
        axis.set_xlabel("experimental ddG")
        axis.set_ylabel("predicted ddG")

    fig.suptitle(
        "Validation targets: Tukey-trimmed ddG scatter (blue circles = raw trimmed, orange triangles = calibrated trimmed)",
        fontsize=14,
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def render_outlier_sensitivity_figure(
    targets: list[dict[str, Any]], excluded_targets: set[str], output_path: Path
) -> None:
    labels = []
    raw_pearson_before = []
    raw_pearson_after = []
    calibrated_pearson_before = []
    calibrated_pearson_after = []
    raw_mae_before = []
    raw_mae_after = []
    calibrated_mae_before = []
    calibrated_mae_after = []

    for target in targets:
        target_id = target["complex_id"]
        star = "*" if target_id in excluded_targets else ""
        labels.append(f"{target_id}{star}")
        raw_pearson_before.append(to_float(target["raw_metrics"].get("pearson_r")) or 0.0)
        raw_pearson_after.append(
            to_float(target["raw_tukey_outlier_analysis"]["filtered_metrics"].get("pearson_r"))
            or 0.0
        )
        calibrated_pearson_before.append(
            to_float(target["calibrated_metrics"].get("pearson_r")) or 0.0
        )
        calibrated_pearson_after.append(
            to_float(
                target["calibrated_tukey_outlier_analysis"]["filtered_metrics"].get("pearson_r")
            )
            or 0.0
        )
        raw_mae_before.append(to_float(target["raw_metrics"].get("mae_kcal_mol")) or 0.0)
        raw_mae_after.append(
            to_float(target["raw_tukey_outlier_analysis"]["filtered_metrics"].get("mae_kcal_mol"))
            or 0.0
        )
        calibrated_mae_before.append(
            to_float(target["calibrated_metrics"].get("mae_kcal_mol")) or 0.0
        )
        calibrated_mae_after.append(
            to_float(
                target["calibrated_tukey_outlier_analysis"]["filtered_metrics"].get(
                    "mae_kcal_mol"
                )
            )
            or 0.0
        )

    x = list(range(len(targets)))
    width = 0.38
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
    axes = axes.ravel()

    axes[0].bar(
        [i - width / 2 for i in x], raw_pearson_before, width=width, color=RAW_COLOR, label="original"
    )
    axes[0].bar(
        [i + width / 2 for i in x], raw_pearson_after, width=width, color=INCLUDED_COLOR, label="trimmed"
    )
    axes[0].axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    axes[0].set_title("Raw Pearson r: before vs trimmed")

    axes[1].bar(
        [i - width / 2 for i in x],
        calibrated_pearson_before,
        width=width,
        color=CAL_COLOR,
        label="original",
    )
    axes[1].bar(
        [i + width / 2 for i in x],
        calibrated_pearson_after,
        width=width,
        color=INCLUDED_COLOR,
        label="trimmed",
    )
    axes[1].axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    axes[1].set_title("Calibrated Pearson r: before vs trimmed")

    axes[2].bar(
        [i - width / 2 for i in x], raw_mae_before, width=width, color=RAW_COLOR
    )
    axes[2].bar(
        [i + width / 2 for i in x], raw_mae_after, width=width, color=INCLUDED_COLOR
    )
    axes[2].set_title("Raw MAE (kcal/mol): before vs trimmed")

    axes[3].bar(
        [i - width / 2 for i in x], calibrated_mae_before, width=width, color=CAL_COLOR
    )
    axes[3].bar(
        [i + width / 2 for i in x], calibrated_mae_after, width=width, color=INCLUDED_COLOR
    )
    axes[3].set_title("Calibrated MAE (kcal/mol): before vs trimmed")

    for axis in axes:
        axis.set_xticks(x, labels)
        axis.tick_params(axis="x", labelsize=9)
        axis.grid(axis="y", alpha=0.25, linestyle="--")

    axes[0].legend(frameon=False)
    axes[1].legend(frameon=False)
    fig.suptitle(
        "Per-target Tukey outlier sensitivity\n"
        "* excluded from accepted target-filtered validation view",
        fontsize=15,
    )
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_json_payload(
    summary: dict[str, Any], targets: list[dict[str, Any]], patel_summary: dict[str, Any] | None
) -> dict[str, Any]:
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary_generated_at": summary.get("generated_at"),
        "selected_model": summary.get("selected_model"),
        "raw_pearson_r": summary.get("raw_pearson_r"),
        "calibrated_pearson_r": summary.get("calibrated_pearson_r"),
        "accepted_calibrated_pearson_r": summary.get("accepted_calibrated_pearson_r"),
        "accepted_calibrated_excluded_complex_ids": summary.get(
            "accepted_calibrated_excluded_complex_ids", []
        ),
        "raw_target_excluded_complex_ids": summary.get("raw_target_excluded_complex_ids", []),
        "targets": [
            {
                "complex_id": target["complex_id"],
                "pair_count": target["pair_count"],
                "raw_metrics": target["raw_metrics"],
                "calibrated_metrics": target["calibrated_metrics"],
                "raw_prediction_span_kcal_mol": target["raw_prediction_span_kcal_mol"],
                "calibrated_prediction_span_kcal_mol": target[
                    "calibrated_prediction_span_kcal_mol"
                ],
                "compression_ratio": target["compression_ratio"],
                "readout": target["readout"],
                "outlier_sensitivity_readout": target["outlier_sensitivity_readout"],
                "raw_outliers": target["raw_outliers"],
                "calibrated_outliers": target["calibrated_outliers"],
                "raw_tukey_outlier_analysis": target["raw_tukey_outlier_analysis"],
                "calibrated_tukey_outlier_analysis": target[
                    "calibrated_tukey_outlier_analysis"
                ],
                "patel_external_summary": target.get("patel_external_summary"),
            }
            for target in targets
        ],
        "patel_external_summary": patel_summary,
    }


def render_outlier_table(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "| job_id | predicted | experimental | abs error |\n|---|---:|---:|---:|\n| NA | NA | NA | NA |"
    lines = ["| job_id | predicted | experimental | abs error |", "|---|---:|---:|---:|"]
    for row in rows:
        lines.append(
            "| {job_id} | {pred} | {exp} | {err} |".format(
                job_id=row["job_id"],
                pred=format_float(row["predicted_ddg_kcal_mol"]),
                exp=format_float(row["experimental_ddg_kcal_mol"]),
                err=format_float(row["abs_error_kcal_mol"]),
            )
        )
    return "\n".join(lines)


def write_markdown(
    output_path: Path,
    metrics_figure: Path,
    scatter_figure: Path,
    trimmed_scatter_figure: Path,
    outlier_sensitivity_figure: Path,
    summary: dict[str, Any],
    targets: list[dict[str, Any]],
) -> None:
    excluded_targets = set(summary.get("accepted_calibrated_excluded_complex_ids", []))
    lines: list[str] = []
    lines.append("# Validation Target Summary")
    lines.append("")
    lines.append(f"- Generated at: `{datetime.now(timezone.utc).isoformat()}`")
    lines.append(f"- Summary snapshot: `{summary.get('generated_at', 'NA')}`")
    lines.append(f"- Selected model: `{summary.get('selected_model', 'NA')}`")
    lines.append(
        f"- Overall raw Pearson r: `{format_float(to_float(summary.get('raw_pearson_r')))}`"
    )
    lines.append(
        f"- Overall calibrated Pearson r: `{format_float(to_float(summary.get('calibrated_pearson_r')))}`"
    )
    lines.append(
        f"- Accepted target-filtered Pearson r: `{format_float(to_float(summary.get('accepted_calibrated_pearson_r')))}`"
        f" on `{summary.get('accepted_calibrated_pair_count', 'NA')}` pairs."
    )
    lines.append(
        f"- Accepted filtered exclusions: `{', '.join(summary.get('accepted_calibrated_excluded_complex_ids', [])) or 'none'}`"
    )
    lines.append(
        f"- Raw target-filtered exclusions: `{', '.join(summary.get('raw_target_excluded_complex_ids', [])) or 'none'}`"
    )
    lines.append("")
    lines.append("## Overview Figures")
    lines.append("")
    lines.append(f"![Target metrics]({metrics_figure.name})")
    lines.append("")
    lines.append(f"![Target scatter]({scatter_figure.name})")
    lines.append("")
    lines.append(f"![Trimmed target scatter]({trimmed_scatter_figure.name})")
    lines.append("")
    lines.append(f"![Outlier sensitivity]({outlier_sensitivity_figure.name})")
    lines.append("")
    lines.append("## Target Table")
    lines.append("")
    lines.append(
        "| Target | n | Accepted filtered? | Raw r | Cal r | Raw MAE | Cal MAE | Raw sign | Cal sign | Raw span | Cal span |"
    )
    lines.append("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for target in targets:
        raw = target["raw_metrics"]
        calibrated = target["calibrated_metrics"]
        lines.append(
            "| {target} | {n} | {accepted} | {raw_r} | {cal_r} | {raw_mae} | {cal_mae} | {raw_sign} | {cal_sign} | {raw_span} | {cal_span} |".format(
                target=target["complex_id"],
                n=target["pair_count"],
                accepted="no" if target["complex_id"] in excluded_targets else "yes",
                raw_r=format_float(to_float(raw.get("pearson_r")), 2),
                cal_r=format_float(to_float(calibrated.get("pearson_r")), 2),
                raw_mae=format_float(to_float(raw.get("mae_kcal_mol")), 2),
                cal_mae=format_float(to_float(calibrated.get("mae_kcal_mol")), 2),
                raw_sign=format_float(to_float(raw.get("sign_accuracy")), 2),
                cal_sign=format_float(to_float(calibrated.get("sign_accuracy")), 2),
                raw_span=format_float(target["raw_prediction_span_kcal_mol"], 2),
                cal_span=format_float(target["calibrated_prediction_span_kcal_mol"], 2),
            )
        )
    lines.append("")
    lines.append("## Tukey Outlier Sensitivity")
    lines.append("")
    lines.append(
        "- Per-target outlier trimming uses Tukey fences on absolute error: `abs_error > Q3 + 1.5 * IQR`."
    )
    lines.append("- The trimmed scatter figure below uses the raw and calibrated trimmed sets independently for each target.")
    lines.append("")
    lines.append(f"![Trimmed target scatter]({trimmed_scatter_figure.name})")
    lines.append("")
    lines.append(
        "| Target | Raw removed | Raw r -> trimmed r | Raw MAE -> trimmed MAE | Cal removed | Cal r -> trimmed r | Cal MAE -> trimmed MAE |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for target in targets:
        raw_tukey = target["raw_tukey_outlier_analysis"]
        calibrated_tukey = target["calibrated_tukey_outlier_analysis"]
        lines.append(
            "| {target} | {raw_removed} | {raw_before} -> {raw_after} | {raw_mae_before} -> {raw_mae_after} | {cal_removed} | {cal_before} -> {cal_after} | {cal_mae_before} -> {cal_mae_after} |".format(
                target=target["complex_id"],
                raw_removed=raw_tukey["removed_pair_count"],
                raw_before=format_float(to_float(target["raw_metrics"].get("pearson_r")), 2),
                raw_after=format_float(
                    to_float(raw_tukey["filtered_metrics"].get("pearson_r")), 2
                ),
                raw_mae_before=format_float(
                    to_float(target["raw_metrics"].get("mae_kcal_mol")), 2
                ),
                raw_mae_after=format_float(
                    to_float(raw_tukey["filtered_metrics"].get("mae_kcal_mol")), 2
                ),
                cal_removed=calibrated_tukey["removed_pair_count"],
                cal_before=format_float(
                    to_float(target["calibrated_metrics"].get("pearson_r")), 2
                ),
                cal_after=format_float(
                    to_float(calibrated_tukey["filtered_metrics"].get("pearson_r")), 2
                ),
                cal_mae_before=format_float(
                    to_float(target["calibrated_metrics"].get("mae_kcal_mol")), 2
                ),
                cal_mae_after=format_float(
                    to_float(calibrated_tukey["filtered_metrics"].get("mae_kcal_mol")), 2
                ),
            )
        )
    lines.append("")
    lines.append("## Per-target Readout")
    lines.append("")
    for target in targets:
        raw = target["raw_metrics"]
        calibrated = target["calibrated_metrics"]
        raw_tukey = target["raw_tukey_outlier_analysis"]
        calibrated_tukey = target["calibrated_tukey_outlier_analysis"]
        lines.append(f"### {target['complex_id']}")
        lines.append("")
        lines.append(f"- Readout: {target['readout']}")
        lines.append(f"- Outlier sensitivity: {target['outlier_sensitivity_readout']}")
        lines.append(
            f"- Acceptance: `{format_bool(target['complex_id'] not in excluded_targets)}` in the accepted target-filtered view."
        )
        lines.append(
            f"- Raw metrics: `r={format_float(to_float(raw.get('pearson_r')), 3)}`, "
            f"`rho={format_float(to_float(raw.get('spearman_rho')), 3)}`, "
            f"`MAE={format_float(to_float(raw.get('mae_kcal_mol')), 3)}`, "
            f"`sign={format_float(to_float(raw.get('sign_accuracy')), 3)}`."
        )
        lines.append(
            f"- Calibrated metrics: `r={format_float(to_float(calibrated.get('pearson_r')), 3)}`, "
            f"`rho={format_float(to_float(calibrated.get('spearman_rho')), 3)}`, "
            f"`MAE={format_float(to_float(calibrated.get('mae_kcal_mol')), 3)}`, "
            f"`sign={format_float(to_float(calibrated.get('sign_accuracy')), 3)}`."
        )
        lines.append(
            f"- Prediction span: raw `{format_float(target['raw_prediction_span_kcal_mol'], 3)}` kcal/mol, "
            f"calibrated `{format_float(target['calibrated_prediction_span_kcal_mol'], 3)}` kcal/mol, "
            f"compression ratio `{format_float(target['compression_ratio'], 3)}`."
        )
        if target.get("patel_external_summary"):
            patel = target["patel_external_summary"]
            lines.append(
                f"- External 3HFM reference: status `{patel.get('status', 'NA')}`, "
                f"paired `{patel.get('paired_job_count', 'NA')}`, "
                f"incomplete `{patel.get('incomplete_job_count', 'NA')}`."
            )
        lines.append(
            f"- Raw Tukey trim: removed `{raw_tukey['removed_pair_count']}` pair(s), "
            f"threshold `{format_float(raw_tukey['threshold_abs_error_kcal_mol'], 3)}` kcal/mol, "
            f"`r={format_float(to_float(raw.get('pearson_r')), 3)} -> {format_float(to_float(raw_tukey['filtered_metrics'].get('pearson_r')), 3)}`, "
            f"`MAE={format_float(to_float(raw.get('mae_kcal_mol')), 3)} -> {format_float(to_float(raw_tukey['filtered_metrics'].get('mae_kcal_mol')), 3)}`."
        )
        lines.append(
            f"- Raw removed jobs: `{', '.join(raw_tukey['removed_job_ids']) or 'none'}`."
        )
        lines.append(
            f"- Calibrated Tukey trim: removed `{calibrated_tukey['removed_pair_count']}` pair(s), "
            f"threshold `{format_float(calibrated_tukey['threshold_abs_error_kcal_mol'], 3)}` kcal/mol, "
            f"`r={format_float(to_float(calibrated.get('pearson_r')), 3)} -> {format_float(to_float(calibrated_tukey['filtered_metrics'].get('pearson_r')), 3)}`, "
            f"`MAE={format_float(to_float(calibrated.get('mae_kcal_mol')), 3)} -> {format_float(to_float(calibrated_tukey['filtered_metrics'].get('mae_kcal_mol')), 3)}`."
        )
        lines.append(
            f"- Calibrated removed jobs: `{', '.join(calibrated_tukey['removed_job_ids']) or 'none'}`."
        )
        lines.append("")
        lines.append("Top raw outliers:")
        lines.append("")
        lines.append(render_outlier_table(target["raw_outliers"]))
        lines.append("")
        lines.append("Top calibrated outliers:")
        lines.append("")
        lines.append(render_outlier_table(target["calibrated_outliers"]))
        lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    summary = load_json(args.summary_json)
    pair_csv = args.pair_csv or (Path(summary["reports_dir"]) / "predict_pairs_calibrated.csv")
    pair_rows = load_csv_rows(pair_csv)
    patel_summary = None
    if args.patel_summary and args.patel_summary.exists():
        patel_summary = load_json(args.patel_summary)

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    targets = build_target_payload(summary, pair_rows, patel_summary)
    excluded_targets = set(summary.get("accepted_calibrated_excluded_complex_ids", []))

    metrics_figure = output_dir / "validation_target_metrics.png"
    scatter_figure = output_dir / "validation_target_scatter.png"
    trimmed_scatter_figure = output_dir / "validation_target_scatter_trimmed.png"
    outlier_sensitivity_figure = output_dir / "validation_target_outlier_sensitivity.png"
    markdown_path = output_dir / "validation_target_summary.md"
    json_path = output_dir / "validation_target_summary.json"

    render_metrics_figure(targets, excluded_targets, metrics_figure)
    render_scatter_figure(targets, excluded_targets, scatter_figure)
    render_trimmed_scatter_figure(
        targets, excluded_targets, trimmed_scatter_figure
    )
    render_outlier_sensitivity_figure(
        targets, excluded_targets, outlier_sensitivity_figure
    )
    write_markdown(
        markdown_path,
        metrics_figure,
        scatter_figure,
        trimmed_scatter_figure,
        outlier_sensitivity_figure,
        summary,
        targets,
    )

    json_payload = build_json_payload(summary, targets, patel_summary)
    json_path.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")

    print(f"summary_markdown={markdown_path}")
    print(f"summary_json={json_path}")
    print(f"metrics_figure={metrics_figure}")
    print(f"scatter_figure={scatter_figure}")
    print(f"trimmed_scatter_figure={trimmed_scatter_figure}")
    print(f"outlier_sensitivity_figure={outlier_sensitivity_figure}")


if __name__ == "__main__":
    main()
