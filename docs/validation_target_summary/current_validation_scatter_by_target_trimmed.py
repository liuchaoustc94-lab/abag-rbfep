#!/usr/bin/env python3
"""Render target-level validation scatters after formal per-target outlier trimming.

Figure contract
---------------
Core conclusion: removing formally defined extreme absolute-error points helps
some raw target-level relationships, but does not resolve calibration-induced
ranking inversions in every target.
Archetype: quantitative grid.
Panels: one row per target; raw and side-linear calibrated scatters after the
existing per-target Tukey IQR screen has been applied independently.
Source data: selected side-linear prediction pairs and the recorded validation
target summary, which supplies the exact removed-job lists.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# Keep all SVG labels editable in vector-editing applications.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 6.4
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.65


REPO_ROOT = Path(__file__).resolve().parents[2]
PAIR_CSV = (
    REPO_ROOT
    / "runs/benchmarks/abbind_core_v1_quick_plan/reports/calibrations"
    / "fit-calibration-development-predict-validation-model-side-linear"
    / "predict_pairs_calibrated.csv"
)
SUMMARY_JSON = Path(__file__).resolve().parent / "validation_target_summary.json"
OUTPUT_DIR = Path(__file__).resolve().parent
OUTPUT_STEM = OUTPUT_DIR / "current_validation_scatter_by_target_trimmed"
TARGET_ORDER = ["1BJ1", "1CZ8", "1MLC", "2NZ9", "3HFM", "3NPS"]
EXCLUDED_TARGETS = {"1BJ1", "1CZ8", "1MLC"}
TARGET_COLORS = {
    "1BJ1": "#3775BA",
    "1CZ8": "#42949E",
    "1MLC": "#9A4D8E",
    "2NZ9": "#B64342",
    "3HFM": "#E28E2C",
    "3NPS": "#4F8A5B",
}


def read_pairs(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for field in ("raw_ddg_kcal_mol", "predicted_ddg_kcal_mol", "experimental_ddg_kcal_mol"):
            row[field] = float(str(row[field]))
    return rows


def read_trimmed_job_ids(path: Path) -> dict[str, dict[str, set[str]]]:
    summary = json.loads(path.read_text(encoding="utf-8"))
    trimmed: dict[str, dict[str, set[str]]] = {}
    for target in summary["targets"]:
        trimmed[str(target["complex_id"])] = {
            "raw": set(target["raw_tukey_outlier_analysis"]["removed_job_ids"]),
            "calibrated": set(target["calibrated_tukey_outlier_analysis"]["removed_job_ids"]),
        }
    return trimmed


def pearson_r(rows: list[dict[str, object]], prediction_field: str) -> float:
    predicted = np.array([float(row[prediction_field]) for row in rows])
    observed = np.array([float(row["experimental_ddg_kcal_mol"]) for row in rows])
    return float(np.corrcoef(predicted, observed)[0, 1])


def axis_limits(values: list[float], padding: float = 0.10) -> tuple[float, float]:
    lo, hi = min(values), max(values)
    span = max(hi - lo, 1.0)
    return lo - span * padding, hi + span * padding


def draw_target_panel(
    ax: plt.Axes,
    rows: list[dict[str, object]],
    target: str,
    prediction_field: str,
    title: str,
    show_y_label: bool,
    show_x_label: bool,
) -> None:
    predicted = [float(row[prediction_field]) for row in rows]
    observed = [float(row["experimental_ddg_kcal_mol"]) for row in rows]
    x_limits = axis_limits(predicted)
    y_limits = axis_limits(observed)
    color = TARGET_COLORS[target]
    ax.scatter(
        predicted,
        observed,
        s=18,
        facecolors=color,
        edgecolors=color,
        linewidths=0.6,
        alpha=0.82,
        zorder=3,
    )
    line_low = max(x_limits[0], y_limits[0])
    line_high = min(x_limits[1], y_limits[1])
    ax.plot([line_low, line_high], [line_low, line_high], color="#8A8A8A", lw=0.65, ls="--", zorder=1)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.grid(axis="y", color="#E7E7E7", linewidth=0.45, zorder=0)
    ax.set_title(title, loc="left", fontsize=6.8, fontweight="bold", pad=2)
    if show_y_label:
        ax.set_ylabel("Experimental ddG\n(kcal mol$^{-1}$)")
    if show_x_label:
        ax.set_xlabel("Predicted ddG (kcal mol$^{-1}$)")


def main() -> None:
    rows = read_pairs(PAIR_CSV)
    trimmed_ids = read_trimmed_job_ids(SUMMARY_JSON)
    fig, axes = plt.subplots(len(TARGET_ORDER), 2, figsize=(6.9, 12.4))
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.07, top=0.94, hspace=0.72, wspace=0.28)
    fig.text(0.25, 0.978, "Raw RBFE after Tukey trimming", ha="center", va="top", fontsize=9, fontweight="bold")
    fig.text(0.74, 0.978, "Calibrated prediction after Tukey trimming", ha="center", va="top", fontsize=9, fontweight="bold")

    for row_index, target in enumerate(TARGET_ORDER):
        target_rows = [row for row in rows if str(row["complex_id"]) == target]
        raw_rows = [row for row in target_rows if str(row["job_id"]) not in trimmed_ids[target]["raw"]]
        calibrated_rows = [
            row for row in target_rows if str(row["job_id"]) not in trimmed_ids[target]["calibrated"]
        ]
        status = "excluded" if target in EXCLUDED_TARGETS else "retained"
        draw_target_panel(
            axes[row_index, 0],
            raw_rows,
            target,
            "raw_ddg_kcal_mol",
            f"{target} | n={len(raw_rows)} (-{len(trimmed_ids[target]['raw'])}) | raw r={pearson_r(raw_rows, 'raw_ddg_kcal_mol'):.2f} | {status}",
            show_y_label=True,
            show_x_label=row_index == len(TARGET_ORDER) - 1,
        )
        draw_target_panel(
            axes[row_index, 1],
            calibrated_rows,
            target,
            "predicted_ddg_kcal_mol",
            f"{target} | n={len(calibrated_rows)} (-{len(trimmed_ids[target]['calibrated'])}) | calibrated r={pearson_r(calibrated_rows, 'predicted_ddg_kcal_mol'):.2f}",
            show_y_label=False,
            show_x_label=row_index == len(TARGET_ORDER) - 1,
        )
        axes[row_index, 0].text(
            -0.19,
            1.02,
            chr(ord("a") + row_index),
            transform=axes[row_index, 0].transAxes,
            fontsize=8.6,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    fig.text(
        0.5,
        0.015,
        "Per-target trimming rule: remove points with absolute error > Q3 + 1.5 x IQR. "
        "Raw and calibrated panels are trimmed independently; dashed line: identity.",
        ha="center",
        va="bottom",
        fontsize=6.2,
        color="#4D4D4D",
    )
    fig.savefig(OUTPUT_STEM.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
