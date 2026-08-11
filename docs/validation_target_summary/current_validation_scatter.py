#!/usr/bin/env python3
"""Render current AB-Bind validation scatter plots from the selected side-linear run.

Figure contract
---------------
Core conclusion: the accepted target-filtered calibration view reaches r > 0.6,
while three target-level failures still dominate the full validation view.
Archetype: quantitative grid.
Panels: (a) raw RBFE prediction versus experiment; (b) side-linear calibrated
prediction versus experiment, with accepted targets filled and excluded targets open.
Source data: predict_pairs_calibrated.csv from the side-linear validation run.
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
plt.rcParams["font.size"] = 7
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.75


REPO_ROOT = Path(__file__).resolve().parents[2]
PAIR_CSV = (
    REPO_ROOT
    / "runs/benchmarks/abbind_core_v1_quick_plan/reports/calibrations"
    / "fit-calibration-development-predict-validation-model-side-linear"
    / "predict_pairs_calibrated.csv"
)
OUTPUT_DIR = Path(__file__).resolve().parent
OUTPUT_STEM = OUTPUT_DIR / "current_validation_scatter"
EXCLUDED_TARGETS = {"1MLC", "1CZ8", "1BJ1"}
TARGET_ORDER = ["1BJ1", "1CZ8", "1MLC", "2NZ9", "3HFM", "3NPS"]
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
    predicted = np.array([float(row[prediction_field]) for row in rows])
    observed = np.array([float(row["experimental_ddg_kcal_mol"]) for row in rows])
    return float(np.corrcoef(predicted, observed)[0, 1])


def axis_limits(values: list[float], padding: float = 0.08) -> tuple[float, float]:
    lo, hi = min(values), max(values)
    span = max(hi - lo, 1.0)
    return lo - span * padding, hi + span * padding


def draw_panel(
    ax: plt.Axes,
    rows: list[dict[str, object]],
    prediction_field: str,
    title: str,
    label: str,
    x_limits: tuple[float, float],
    observed_limits: tuple[float, float],
    show_legend: bool = False,
) -> None:
    for target in TARGET_ORDER:
        target_rows = [row for row in rows if row["complex_id"] == target]
        if not target_rows:
            continue
        x = [float(row[prediction_field]) for row in target_rows]
        y = [float(row["experimental_ddg_kcal_mol"]) for row in target_rows]
        excluded = target in EXCLUDED_TARGETS
        ax.scatter(
            x,
            y,
            s=20,
            marker="o",
            linewidths=0.85,
            facecolors="white" if excluded else TARGET_COLORS[target],
            edgecolors=TARGET_COLORS[target],
            alpha=0.9 if excluded else 0.82,
            label=target,
            zorder=3,
        )

    line_low = max(x_limits[0], observed_limits[0])
    line_high = min(x_limits[1], observed_limits[1])
    ax.plot([line_low, line_high], [line_low, line_high], color="#767676", lw=0.8, ls="--", zorder=1)
    ax.set_xlim(*x_limits)
    ax.set_ylim(*observed_limits)
    title_text = ax.set_title(title, loc="left", fontsize=8, fontweight="bold", pad=5)
    title_text.set_x(0.08)
    ax.set_xlabel("Predicted ddG (kcal mol$^{-1}$)")
    ax.grid(axis="y", color="#E7E7E7", linewidth=0.55, zorder=0)
    ax.text(
        -0.12,
        1.03,
        label,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    if show_legend:
        ax.legend(
            title="Target",
            loc="upper left",
            bbox_to_anchor=(1.01, 1.02),
            borderaxespad=0,
            labelspacing=0.4,
            handletextpad=0.4,
            markerscale=0.85,
            fontsize=6.2,
            title_fontsize=6.5,
        )


def annotate_raw_outliers(ax: plt.Axes, rows: list[dict[str, object]]) -> None:
    annotations = {
        "1bj1-antigen-w-q89a": ("1BJ1 Q89A", 3, -8),
        "1cz8-antigen-w-q89a": ("1CZ8 Q89A", 3, -8),
        "1mlc-antibody-l-n92a": ("1MLC N92A", 3, 5),
    }
    for row in rows:
        job_id = str(row["job_id"])
        if job_id not in annotations:
            continue
        label, dx, dy = annotations[job_id]
        ax.annotate(
            label,
            (float(row["raw_ddg_kcal_mol"]), float(row["experimental_ddg_kcal_mol"])),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=5.4,
            color="#4D4D4D",
            ha="left",
            va="bottom",
        )


def main() -> None:
    rows = read_pairs(PAIR_CSV)
    accepted_rows = [row for row in rows if str(row["complex_id"]) not in EXCLUDED_TARGETS]
    observed_values = [float(row["experimental_ddg_kcal_mol"]) for row in rows]
    raw_values = [float(row["raw_ddg_kcal_mol"]) for row in rows]
    calibrated_values = [float(row["predicted_ddg_kcal_mol"]) for row in rows]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35))
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.24, top=0.88, wspace=0.30)
    observed_limits = axis_limits(observed_values, padding=0.1)
    draw_panel(
        axes[0],
        rows,
        "raw_ddg_kcal_mol",
        f"Raw RBFE (all n={len(rows)}; r={pearson_r(rows, 'raw_ddg_kcal_mol'):.2f})",
        "a",
        axis_limits(raw_values, padding=0.07),
        observed_limits,
    )
    annotate_raw_outliers(axes[0], rows)
    axes[0].set_ylabel("Experimental ddG (kcal mol$^{-1}$)")

    draw_panel(
        axes[1],
        rows,
        "predicted_ddg_kcal_mol",
        (
            "Side-linear calibration "
            f"(full r={pearson_r(rows, 'predicted_ddg_kcal_mol'):.2f}; "
            f"accepted r={pearson_r(accepted_rows, 'predicted_ddg_kcal_mol'):.2f}, n={len(accepted_rows)})"
        ),
        "b",
        axis_limits(calibrated_values, padding=0.12),
        observed_limits,
        show_legend=True,
    )

    fig.text(
        0.5,
        0.055,
        "Filled: targets retained in accepted view (2NZ9, 3HFM, 3NPS). "
        "Open: excluded target-level failures (1BJ1, 1CZ8, 1MLC). Dashed line: identity.",
        ha="center",
        va="bottom",
        fontsize=6.4,
        color="#4D4D4D",
    )

    fig.savefig(OUTPUT_STEM.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".tiff"), dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
