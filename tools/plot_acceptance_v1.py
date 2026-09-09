#!/usr/bin/env python3
"""Plot the completed Acceptance V1 prediction/experiment pairs.

The script deliberately keeps incomplete jobs out of the statistics.  Experimental
rows without a completed ddg_summary.json remain in the exported manifest with a
missing prediction, so the completion gap is visible and auditable.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import OrderedDict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


TARGETS = ["1AHW", "1IAR", "1NMB", "1OGA", "2BDN", "4I77"]
COLORS = OrderedDict(
    [
        ("1AHW", "#3B82A0"),
        ("1IAR", "#D17A22"),
        ("1NMB", "#6B6E9F"),
        ("1OGA", "#5B8E7D"),
        ("2BDN", "#A06A86"),
        ("4I77", "#B54A45"),
    ]
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _latest_results(root: Path, target: str) -> dict[str, dict[str, Any]]:
    """Return one result per mutation, preferring the newest generated result."""

    target_dir = root / "runs" / "acceptance"
    paths = sorted(
        target_dir.glob(f"acc_{target.lower()}_v1_*/jobs/*/results/ddg_summary.json")
    )
    results: dict[str, dict[str, Any]] = {}
    timestamps: dict[str, str] = {}
    for path in paths:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        mutation_id = str(payload.get("mutation_group_id", "")).strip()
        if not mutation_id or "ddg_kcal_mol" not in payload:
            continue
        generated = str(payload.get("generated_at", ""))
        if mutation_id not in results or generated >= timestamps[mutation_id]:
            results[mutation_id] = payload
            timestamps[mutation_id] = generated
    return results


def collect_points(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for target in TARGETS:
        target_lc = target.lower()
        exp_path = root / "benchmarks" / "acceptance" / target_lc / "experimental_ddg_v1.csv"
        experiments = _read_csv(exp_path)
        results = _latest_results(root, target)
        for exp in experiments:
            mutation_id = exp["mutation_group_id"]
            payload = results.get(mutation_id)
            pred = None if payload is None else float(payload["ddg_kcal_mol"])
            repeat_sd = None
            confidence = "missing"
            generated_at = ""
            job_id = ""
            if payload is not None:
                repeat_sd = payload.get("ddg_repeat_stdev_kcal_mol")
                repeat_sd = None if repeat_sd is None else float(repeat_sd)
                confidence = str(payload.get("result_confidence", "unknown"))
                generated_at = str(payload.get("generated_at", ""))
                job_id = str(payload.get("job_id", ""))
            rows.append(
                {
                    "target": target,
                    "mutation_group_id": mutation_id,
                    "mutation_short_label": mutation_id.rsplit("_", 1)[-1],
                    "experimental_ddg_kcal_mol": float(exp["ddg_kcal_mol"]),
                    "predicted_ddg_kcal_mol": pred,
                    "prediction_repeat_sd_kcal_mol": repeat_sd,
                    "result_confidence": confidence,
                    "completed": payload is not None,
                    "experiment_note": exp.get("note", ""),
                    "generated_at": generated_at,
                    "job_id": job_id,
                }
            )

    points = pd.DataFrame(rows)
    completed = points.loc[points["completed"]].copy()
    metrics = []
    for target in TARGETS:
        sub = completed.loc[completed["target"] == target]
        x = sub["experimental_ddg_kcal_mol"].to_numpy(dtype=float)
        y = sub["predicted_ddg_kcal_mol"].to_numpy(dtype=float)
        # A two-point correlation is numerically defined but not scientifically
        # interpretable; keep it blank in the report and describe the panel as
        # descriptive-only instead.
        if len(sub) >= 3 and np.std(x) > 0 and np.std(y) > 0:
            r = float(pearsonr(x, y).statistic)
            rho = float(spearmanr(x, y).statistic)
        else:
            r = math.nan
            rho = math.nan
        metrics.append(
            {
                "target": target,
                "n_experimental": int((points["target"] == target).sum()),
                "n_completed": int(len(sub)),
                "n_missing": int((points["target"] == target).sum() - len(sub)),
                "pearson_r": r,
                "spearman_rho": rho,
                "mae_kcal_mol": float(np.mean(np.abs(y - x))) if len(sub) else math.nan,
                "rmse_kcal_mol": float(np.sqrt(np.mean((y - x) ** 2))) if len(sub) else math.nan,
                "sign_accuracy": float(np.mean(np.sign(y) == np.sign(x))) if len(sub) else math.nan,
                "within_1_kcal_mol": float(np.mean(np.abs(y - x) <= 1.0)) if len(sub) else math.nan,
                "within_2_kcal_mol": float(np.mean(np.abs(y - x) <= 2.0)) if len(sub) else math.nan,
            }
        )

    pooled = completed
    x = pooled["experimental_ddg_kcal_mol"].to_numpy(dtype=float)
    y = pooled["predicted_ddg_kcal_mol"].to_numpy(dtype=float)
    metrics.append(
        {
            "target": "POOLED",
            "n_experimental": int(len(points)),
            "n_completed": int(len(pooled)),
            "n_missing": int(len(points) - len(pooled)),
            "pearson_r": float(pearsonr(x, y).statistic),
            "spearman_rho": float(spearmanr(x, y).statistic),
            "mae_kcal_mol": float(np.mean(np.abs(y - x))),
            "rmse_kcal_mol": float(np.sqrt(np.mean((y - x) ** 2))),
            "sign_accuracy": float(np.mean(np.sign(y) == np.sign(x))),
            "within_1_kcal_mol": float(np.mean(np.abs(y - x) <= 1.0)),
            "within_2_kcal_mol": float(np.mean(np.abs(y - x) <= 2.0)),
        }
    )
    return points, pd.DataFrame(metrics)


def _format_metric(value: float) -> str:
    return "—" if not np.isfinite(value) else f"{value:.3f}"


def _style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.0,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 8.0,
            "axes.titlesize": 8.0,
            "xtick.labelsize": 7.0,
            "ytick.labelsize": 7.0,
            "legend.fontsize": 6.5,
            "legend.title_fontsize": 7.0,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def _limits(points: pd.DataFrame) -> tuple[float, float]:
    values = pd.concat(
        [points["experimental_ddg_kcal_mol"], points["predicted_ddg_kcal_mol"]], ignore_index=True
    ).dropna()
    lo = float(values.min())
    hi = float(values.max())
    span = max(hi - lo, 1.0)
    pad = max(0.35, 0.06 * span)
    return math.floor((lo - pad) * 2) / 2, math.ceil((hi + pad) * 2) / 2


def _add_identity(ax: Any, lo: float, hi: float) -> None:
    ax.plot([lo, hi], [lo, hi], color="#6E6E6E", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")


def _draw_points(ax: Any, sub: pd.DataFrame, color: str, *, pooled: bool = False) -> None:
    if sub.empty:
        return
    conf_alpha = np.where(sub["result_confidence"].eq("quantitative"), 0.92, 0.65)
    for row, alpha in zip(sub.itertuples(index=False), conf_alpha):
        x = float(row.experimental_ddg_kcal_mol)
        y = float(row.predicted_ddg_kcal_mol)
        sd = row.prediction_repeat_sd_kcal_mol
        if sd is not None and np.isfinite(sd) and sd > 0:
            ax.errorbar(
                x,
                y,
                yerr=float(sd),
                fmt="none",
                ecolor=color,
                elinewidth=0.55,
                capsize=1.3,
                alpha=0.32 if pooled else 0.42,
                zorder=2,
            )
        ax.scatter(
            [x],
            [y],
            s=24 if pooled else 25,
            color=color,
            edgecolor="white",
            linewidth=0.45,
            alpha=float(alpha),
            zorder=3,
        )


def plot_pooled(points: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    completed = points.loc[points["completed"]].copy()
    pooled = metrics.loc[metrics["target"] == "POOLED"].iloc[0]
    lo, hi = _limits(completed)
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for target, color in COLORS.items():
        _draw_points(ax, completed.loc[completed["target"] == target], color, pooled=True)
    _add_identity(ax, lo, hi)
    ax.set_xlabel(r"Experimental $\Delta\Delta G$ (kcal mol$^{-1}$)")
    ax.set_ylabel(r"Predicted $\Delta\Delta G$ (kcal mol$^{-1}$)")
    ax.set_title("Acceptance V1: pooled prediction versus experiment", loc="left", pad=8, fontweight="bold")
    ax.text(
        0.03,
        0.97,
        "n = {n} completed / {nt} experimental\n"
        "Pearson r = {r}\nSpearman ρ = {rho}\nMAE = {mae} kcal mol⁻¹\nRMSE = {rmse} kcal mol⁻¹".format(
            n=int(pooled["n_completed"]),
            nt=int(pooled["n_experimental"]),
            r=_format_metric(pooled["pearson_r"]),
            rho=_format_metric(pooled["spearman_rho"]),
            mae=_format_metric(pooled["mae_kcal_mol"]),
            rmse=_format_metric(pooled["rmse_kcal_mol"]),
        ),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=7.0,
        bbox={"facecolor": "white", "edgecolor": "#D0D0D0", "linewidth": 0.5, "pad": 4},
    )
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="white", markersize=5.5, label=target)
        for target, color in COLORS.items()
    ]
    ax.legend(handles=handles, title="Complex", ncol=3, loc="lower right", frameon=False, handletextpad=0.2, columnspacing=0.8)
    ax.grid(False)
    fig.tight_layout()
    _save_figure(fig, out)
    plt.close(fig)


def plot_target_grid(points: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    lo, hi = _limits(points.loc[points["completed"]])
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.7), sharex=True, sharey=True)
    for ax, target in zip(axes.flat, TARGETS):
        sub = points.loc[(points["target"] == target) & points["completed"]]
        row = metrics.loc[metrics["target"] == target].iloc[0]
        _add_identity(ax, lo, hi)
        _draw_points(ax, sub, COLORS[target])
        ax.set_title(target, loc="left", fontweight="bold", pad=4)
        if int(row["n_completed"]) == 0:
            ax.text(0.5, 0.5, "No completed\npredictions", ha="center", va="center", transform=ax.transAxes, color="#777777")
        else:
            small_n_note = "\nsmall n; descriptive only" if int(row["n_completed"]) < 3 else ""
            ax.text(
                0.04,
                0.95,
                "n = {n}/{nt}\nr = {r}\nρ = {rho}\nMAE = {mae}".format(
                    n=int(row["n_completed"]),
                    nt=int(row["n_experimental"]),
                    r=_format_metric(row["pearson_r"]),
                    rho=_format_metric(row["spearman_rho"]),
                    mae=_format_metric(row["mae_kcal_mol"]),
                ),
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=6.5,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2},
            )
            if small_n_note:
                ax.text(0.04, 0.60, small_n_note, transform=ax.transAxes, fontsize=6.0, color="#777777", va="top")
        ax.tick_params(length=3, width=0.6)
    axes[0, 0].set_ylabel(r"Predicted $\Delta\Delta G$ (kcal mol$^{-1}$)")
    axes[1, 0].set_ylabel(r"Predicted $\Delta\Delta G$ (kcal mol$^{-1}$)")
    for ax in axes[1, :]:
        ax.set_xlabel(r"Experimental $\Delta\Delta G$ (kcal mol$^{-1}$)")
    fig.suptitle("Acceptance V1: per-target prediction versus experiment", x=0.02, ha="left", y=0.995, fontweight="bold", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=1.0, w_pad=1.0)
    _save_figure(fig, out)
    plt.close(fig)


def _save_figure(fig: Any, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")


def write_notes(out: Path, points: pd.DataFrame, metrics: pd.DataFrame) -> None:
    pooled = metrics.loc[metrics["target"] == "POOLED"].iloc[0]
    missing = points.loc[~points["completed"]]
    (out / "figure_contract.md").write_text(
        "# Figure contract: Acceptance V1 scatter plots\n\n"
        "Core conclusion: Acceptance V1 shows moderate pooled agreement, with substantial heterogeneity across complexes and incomplete jobs kept outside the statistics.\n"
        "Figure archetype: quantitative grid (pooled hero panel plus per-target panels).\n"
        "Target output: reusable research/report figure; double-column width, editable SVG/PDF plus 600-dpi TIFF.\n"
        "Backend: Python/matplotlib only.\n"
        "Evidence: completed ddg_summary.json paired to experimental_ddg_v1.csv; missing predictions are retained in the source-data CSV.\n"
        "Statistics: Pearson r, Spearman rho, MAE, RMSE, sign accuracy, and within-1/2-kcal/mol fractions.\n"
        "Review risk: acceptance V1 is an execution/acceptance set, not an untouched scientific holdout; 2BDN has no completed predictions and must not be coded as zero.\n",
        encoding="utf-8",
    )
    (out / "qa_notes.md").write_text(
        f"# QA notes\n\n"
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}\n"
        f"Completed pairs: {int(pooled['n_completed'])}; experimental rows: {int(pooled['n_experimental'])}; missing predictions: {int(pooled['n_missing'])}.\n"
        f"Pooled Pearson r={pooled['pearson_r']:.6f}; Spearman rho={pooled['spearman_rho']:.6f}; MAE={pooled['mae_kcal_mol']:.6f}; RMSE={pooled['rmse_kcal_mol']:.6f}.\n"
        "Identity line is y=x. Prediction repeat SD is shown as a vertical error bar where available.\n"
        "All plotting, export, and visual QA renders were produced with Python/matplotlib.\n"
        "Missing rows:\n"
        + ("\n".join(f"- {r.target}: {r.mutation_group_id}" for r in missing.itertuples()) if not missing.empty else "- None")
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("reports/selections/acceptance_v1_scatter_20260826"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    _style()
    points, metrics = collect_points(root)
    args.output.mkdir(parents=True, exist_ok=True)
    points.to_csv(args.output / "acceptance_v1_points.csv", index=False)
    metrics.to_csv(args.output / "acceptance_v1_metrics.csv", index=False)
    (args.output / "acceptance_v1_metrics.json").write_text(metrics.to_json(orient="records", indent=2), encoding="utf-8")
    plot_pooled(points, metrics, args.output / "acceptance_v1_pooled_scatter")
    plot_target_grid(points, metrics, args.output / "acceptance_v1_target_grid")
    write_notes(args.output, points, metrics)
    print(metrics.to_string(index=False))
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
