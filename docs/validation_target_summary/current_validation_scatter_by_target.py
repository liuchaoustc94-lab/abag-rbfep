#!/usr/bin/env python3
"""Render one raw/calibrated scatter pair for every AB-Bind validation target.

Figure contract
---------------
Core conclusion: target-level heterogeneity, rather than a single global error
mode, explains why the accepted target-filtered validation result differs from
the full validation result.
Archetype: quantitative grid.
Panels: one row per target; raw prediction on the left and side-linear
calibrated prediction on the right.
Source data: predict_pairs_calibrated.csv from the selected side-linear run.
"""

from __future__ import annotations

import csv
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
OUTPUT_DIR = Path(__file__).resolve().parent
OUTPUT_STEM = OUTPUT_DIR / "current_validation_scatter_by_target"
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


def pearson_r(rows: list[dict[str, object]], prediction_field: str) -> float:
    if len(rows) < 2:
        return float("nan")
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
    else:
        ax.set_ylabel("")
    if show_x_label:
        ax.set_xlabel("Predicted ddG (kcal mol$^{-1}$)")
    else:
        ax.set_xlabel("")


def main() -> None:
    rows = read_pairs(PAIR_CSV)
    fig, axes = plt.subplots(len(TARGET_ORDER), 2, figsize=(6.9, 12.4))
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.07, top=0.94, hspace=0.72, wspace=0.28)
    fig.text(0.25, 0.978, "Raw RBFE prediction", ha="center", va="top", fontsize=9, fontweight="bold")
    fig.text(0.74, 0.978, "Side-linear calibrated prediction", ha="center", va="top", fontsize=9, fontweight="bold")

    for row_index, target in enumerate(TARGET_ORDER):
        target_rows = [row for row in rows if str(row["complex_id"]) == target]
        excluded = target in EXCLUDED_TARGETS
        status = "excluded" if excluded else "retained"
        raw_r = pearson_r(target_rows, "raw_ddg_kcal_mol")
        calibrated_r = pearson_r(target_rows, "predicted_ddg_kcal_mol")
        draw_target_panel(
            axes[row_index, 0],
            target_rows,
            target,
            "raw_ddg_kcal_mol",
            f"{target} | n={len(target_rows)} | raw r={raw_r:.2f} | {status}",
            show_y_label=True,
            show_x_label=row_index == len(TARGET_ORDER) - 1,
        )
        draw_target_panel(
            axes[row_index, 1],
            target_rows,
            target,
            "predicted_ddg_kcal_mol",
            f"{target} | calibrated r={calibrated_r:.2f}",
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
        "Each row uses target-specific axes. Dashed line: identity. "
        "Excluded/retained reflects the current accepted target-filtered validation view.",
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
