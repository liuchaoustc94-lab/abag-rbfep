#!/usr/bin/env python3
"""Create current per-target RBFE validation scatter plots with Python.

The collector deliberately keeps protocol families separate.  The resulting CSV is
the source data for the figures; no stale aggregate report is used for the direct
result roots.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any

import matplotlib.pyplot as plt

# Required editable-text settings for the Python figure track.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7.2,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "legend.frameon": False,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "savefig.facecolor": "white",
    }
)

from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
SOURCE_REGISTERED = ROOT / "benchmarks/ab_bind/curated/ab_bind_source_registered.csv"
PATEL_EXP = ROOT / "benchmarks/patel_2021_3hfm/experimental_ddg.csv"

PALETTE = {
    "development_fit": "#3775BA",
    "independent_OOS": "#0F4D92",
    "validation_panel": "#767676",
    "stress_1DVF": "#D9853F",
    "Patel_3HFM": "#B64342",
    "CHARMM36_mut": "#9A4D8E",
}
MARKERS = {
    "development_fit": "o",
    "independent_OOS": "o",
    "validation_panel": "s",
    "stress_1DVF": "D",
    "Patel_3HFM": "^",
    "CHARMM36_mut": "P",
}
LABELS = {
    "development_fit": "development fit",
    "independent_OOS": "independent OOS",
    "validation_panel": "validation panel",
    "stress_1DVF": "1DVF stress",
    "Patel_3HFM": "3HFM Patel/DSSB",
    "CHARMM36_mut": "CHARMM36-mut exploratory",
}


def _float(value: Any) -> float | None:
    try:
        if value in (None, ""):
            return None
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _load_ab_bind_reference() -> dict[tuple[str, str], list[dict[str, str]]]:
    reference: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in _read_csv(SOURCE_REGISTERED):
        reference[(row["complex_id"].upper(), row["signature"].upper())].append(row)
    return reference


def _load_patel_reference() -> dict[str, float]:
    return {
        row["job_id"]: float(row["experimental_ddg_kcal_mol"])
        for row in _read_csv(PATEL_EXP)
    }


def _reference_for(
    complex_id: str,
    summary: dict[str, Any],
    reference: dict[tuple[str, str], list[dict[str, str]]],
) -> tuple[float | None, str, bool, float | None]:
    """Resolve AB-Bind references, retaining duplicate-source ambiguity metadata."""
    full = str(summary.get("mutation_signature", "")).upper()
    base = full.split("@", 1)[0]
    side = str(summary.get("entity_side", "")).upper()
    candidates = [full, base + (f"@{side}" if side else ""), base]
    matches: list[dict[str, str]] = []
    for key in candidates:
        matches = reference.get((complex_id.upper(), key), [])
        if matches:
            break
    if not matches:
        return None, "", False, None

    # 1DVF A:Y49A occurs twice in AB-Bind.  The single-point panel used 1.75;
    # preserve the alternative in the source-data table instead of hiding it.
    values = [float(row["ddg_kcal_mol"]) for row in matches]
    chosen = values[-1] if complex_id.upper() == "1DVF" and base == "A:Y49A" else values[0]
    alt = values[1] if len(values) > 1 else None
    source_id = ";".join(row.get("row_id", "") for row in matches)
    return chosen, source_id, len(values) > 1, alt


def _qc_status(job_dir: Path) -> tuple[str, int, int]:
    qc_path = job_dir / "results/qc_report.json"
    if not qc_path.exists():
        return "not_evaluated", 0, 0
    data = json.loads(qc_path.read_text(encoding="utf-8"))
    status = str(data.get("status", "not_evaluated"))
    if status not in {"pass", "warning", "fail"}:
        status = "not_evaluated"
    return status, len(data.get("warnings", [])), len(data.get("failures", []))


def _collect_root(
    root: Path,
    dataset: str,
    complex_id: str,
    reference: dict[tuple[str, str], list[dict[str, str]]],
    rows: list[dict[str, Any]],
) -> None:
    for summary_path in sorted(root.glob("**/results/ddg_summary.json")):
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if not summary.get("ready", True):
            continue
        predicted = _float(summary.get("ddg_kcal_mol"))
        if predicted is None:
            continue
        job_dir = summary_path.parents[1]
        qc, warning_count, failure_count = _qc_status(job_dir)
        experimental, source_id, ambiguous, alternative = _reference_for(
            complex_id, summary, reference
        )
        if experimental is None:
            continue
        rows.append(
            {
                "target": complex_id.upper(),
                "dataset": dataset,
                "dataset_label": LABELS[dataset],
                "job_id": str(summary.get("job_id", job_dir.name)),
                "mutation_signature": str(summary.get("mutation_signature", "")),
                "mutation_short_label": str(summary.get("mutation_short_label", "")),
                "experimental_ddg_kcal_mol": experimental,
                "predicted_ddg_kcal_mol": predicted,
                "error_kcal_mol": predicted - experimental,
                "repeat_stdev_kcal_mol": _float(summary.get("ddg_repeat_stdev_kcal_mol")),
                "bar_stderr_kcal_mol": _float(summary.get("ddg_bar_stderr_kcal_mol")),
                "qc_status": qc,
                "qc_warning_count": warning_count,
                "qc_failure_count": failure_count,
                "experimental_reference_id": source_id,
                "experimental_reference_ambiguous": ambiguous,
                "experimental_reference_alternative": alternative,
                "protocol_preset": str(summary.get("protocol_preset", "")),
                "protocol_lambda_windows": summary.get("protocol_lambda_windows", ""),
                "source_result": str(summary_path),
                "source_root": str(root),
            }
        )


def _collect_patel(rows: list[dict[str, Any]]) -> None:
    reference = _load_patel_reference()
    roots = [
        (ROOT / "runs/real_cases/3hfm_patel_align_20260817", "Patel_3HFM"),
        (ROOT / "runs/real_cases/3hfm_patel_dssb_pr4", "Patel_3HFM"),
    ]
    for root, dataset in roots:
        for summary_path in sorted(root.glob("**/results/ddg_summary.json")):
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if not summary.get("ready", True) or summary.get("job_id") not in reference:
                continue
            predicted = _float(summary.get("ddg_kcal_mol"))
            if predicted is None:
                continue
            job_dir = summary_path.parents[1]
            qc, warning_count, failure_count = _qc_status(job_dir)
            experimental = reference[summary["job_id"]]
            rows.append(
                {
                    "target": "3HFM",
                    "dataset": dataset,
                    "dataset_label": LABELS[dataset],
                    "job_id": summary["job_id"],
                    "mutation_signature": summary.get("mutation_signature", ""),
                    "mutation_short_label": summary.get("mutation_short_label", ""),
                    "experimental_ddg_kcal_mol": experimental,
                    "predicted_ddg_kcal_mol": predicted,
                    "error_kcal_mol": predicted - experimental,
                    "repeat_stdev_kcal_mol": _float(summary.get("ddg_repeat_stdev_kcal_mol")),
                    "bar_stderr_kcal_mol": _float(summary.get("ddg_bar_stderr_kcal_mol")),
                    "qc_status": qc,
                    "qc_warning_count": warning_count,
                    "qc_failure_count": failure_count,
                    "experimental_reference_id": "Patel 2021 Table 2",
                    "experimental_reference_ambiguous": False,
                    "experimental_reference_alternative": None,
                    "protocol_preset": str(summary.get("protocol_preset", "")),
                    "protocol_lambda_windows": summary.get("protocol_lambda_windows", ""),
                    "source_result": str(summary_path),
                    "source_root": str(root),
                }
            )


def collect_rows() -> list[dict[str, Any]]:
    reference = _load_ab_bind_reference()
    rows: list[dict[str, Any]] = []
    # The fit root contains target IDs in nested batch directories; use explicit
    # target roots rather than relying on a stale aggregate report.
    fit_root = ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811"
    fit_targets = {"1DQJ", "1JRH", "1N8Z", "1VFB", "2JEL", "3BE1", "3NGB"}
    for target in sorted(fit_targets):
        _collect_root(fit_root / f"abbind_{target.lower()}_core_v1", "development_fit", target, reference, rows)
    for target in ("1AK4", "1KTZ", "3K2M"):
        _collect_root(
            ROOT / f"runs/real_cases/oos_{target.lower()}_newprotocol_20260811",
            "independent_OOS",
            target,
            reference,
            rows,
        )
    for target in ("1BJ1", "1CZ8"):
        _collect_root(
            ROOT / "runs/benchmarks/abbind_validation_panels_20260812" / f"abbind_{target.lower()}_core_v1",
            "validation_panel",
            target,
            reference,
            rows,
        )
    _collect_root(
        ROOT / "runs/real_cases/1dvf_priority_20260806",
        "stress_1DVF",
        "1DVF",
        reference,
        rows,
    )
    for target in ("1AK4", "1BJ1", "1DVF", "1KTZ"):
        _collect_root(
            ROOT / f"runs/real_cases/e6_charmm_{target.lower()}",
            "CHARMM36_mut",
            target,
            reference,
            rows,
        )
    _collect_patel(rows)
    return rows


def _metric(values: list[dict[str, Any]]) -> dict[str, Any]:
    xs = [float(row["experimental_ddg_kcal_mol"]) for row in values]
    ys = [float(row["predicted_ddg_kcal_mol"]) for row in values]
    n = len(values)
    result: dict[str, Any] = {
        "n": n,
        "pearson_r": None,
        "spearman_rho": None,
        "mae_kcal_mol": None,
        "rmse_kcal_mol": None,
        "median_abs_error_kcal_mol": None,
        "sign_accuracy": None,
        "within_1_kcal_mol": None,
        "within_2_kcal_mol": None,
        "qc_pass": sum(row["qc_status"] == "pass" for row in values),
        "qc_warning": sum(row["qc_status"] == "warning" for row in values),
        "qc_fail": sum(row["qc_status"] == "fail" for row in values),
        "reference_ambiguous_n": sum(bool(row["experimental_reference_ambiguous"]) for row in values),
        "experimental_min": min(xs) if xs else None,
        "experimental_max": max(xs) if xs else None,
        "predicted_min": min(ys) if ys else None,
        "predicted_max": max(ys) if ys else None,
    }
    if n:
        errors = [y - x for x, y in zip(xs, ys, strict=True)]
        result["mae_kcal_mol"] = sum(abs(error) for error in errors) / n
        result["rmse_kcal_mol"] = math.sqrt(sum(error * error for error in errors) / n)
        result["median_abs_error_kcal_mol"] = median(abs(error) for error in errors)
        result["sign_accuracy"] = sum((x >= 0) == (y >= 0) for x, y in zip(xs, ys, strict=True)) / n
        result["within_1_kcal_mol"] = sum(abs(error) <= 1 for error in errors) / n
        result["within_2_kcal_mol"] = sum(abs(error) <= 2 for error in errors) / n
    if n >= 2 and len(set(xs)) > 1 and len(set(ys)) > 1:
        result["pearson_r"] = float(pearsonr(xs, ys).statistic)
        result["spearman_rho"] = float(spearmanr(xs, ys).statistic)
    return result


def _fmt(value: Any, digits: int = 2) -> str:
    return "—" if value is None else f"{float(value):.{digits}f}"


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _save_figure(fig: plt.Figure, stem: Path, dpi: int = 600) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=dpi, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def _limits(values: list[dict[str, Any]]) -> tuple[float, float]:
    numbers = [float(row["experimental_ddg_kcal_mol"]) for row in values]
    numbers += [float(row["predicted_ddg_kcal_mol"]) for row in values]
    lo, hi = min(numbers), max(numbers)
    span = max(hi - lo, 1.0)
    margin = max(0.12 * span, 0.5)
    return lo - margin, hi + margin


def _scatter_on_ax(
    ax: plt.Axes,
    target: str,
    values: list[dict[str, Any]],
    compact: bool = False,
) -> None:
    lo, hi = _limits(values)
    ax.plot([lo, hi], [lo, hi], color="#BDBDBD", lw=0.8, ls="--", zorder=0)
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in values:
        grouped[row["dataset"]].append(row)
    for dataset in sorted(grouped, key=lambda item: list(LABELS).index(item)):
        group = grouped[dataset]
        color = PALETTE[dataset]
        marker = MARKERS[dataset]
        x = [float(row["experimental_ddg_kcal_mol"]) for row in group]
        y = [float(row["predicted_ddg_kcal_mol"]) for row in group]
        exact = [not row["experimental_reference_ambiguous"] for row in group]
        ax.scatter(
            [xx for xx, is_exact in zip(x, exact, strict=True) if is_exact],
            [yy for yy, is_exact in zip(y, exact, strict=True) if is_exact],
            s=24 if not compact else 13,
            marker=marker,
            c=color,
            edgecolors="white",
            linewidths=0.35,
            alpha=0.9,
            label=LABELS[dataset],
            zorder=3,
        )
        ambiguous = [(xx, yy) for xx, yy, is_exact in zip(x, y, exact, strict=True) if not is_exact]
        if ambiguous:
            ax.scatter(
                [pair[0] for pair in ambiguous],
                [pair[1] for pair in ambiguous],
                s=28 if not compact else 16,
                marker=marker,
                facecolors="none",
                edgecolors=color,
                linewidths=0.9,
                zorder=4,
            )
        if len(group) >= 3 and len(set(x)) > 1:
            slope, intercept = _linear_fit(x, y)
            xx = [lo, hi]
            ax.plot(xx, [slope * value + intercept for value in xx], color=color, lw=0.9, alpha=0.75)
        if not compact and len(values) <= 12:
            for index, row in enumerate(group):
                label = _short_label(row)
                dx = 3 if index % 2 == 0 else -3
                dy = 3 if index % 3 else -4
                ax.annotate(
                    label,
                    (float(row["experimental_ddg_kcal_mol"]), float(row["predicted_ddg_kcal_mol"])),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    fontsize=5.4,
                    color=color,
                    alpha=0.82,
                )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Experimental ΔΔG (kcal mol⁻¹)")
    ax.set_ylabel("Predicted ΔΔG (kcal mol⁻¹)")
    ax.set_title(target, fontweight="bold", pad=4)
    ax.grid(False)


def _linear_fit(x: list[float], y: list[float]) -> tuple[float, float]:
    xbar = sum(x) / len(x)
    ybar = sum(y) / len(y)
    denominator = sum((value - xbar) ** 2 for value in x)
    if denominator == 0:
        return 0.0, ybar
    slope = sum((xx - xbar) * (yy - ybar) for xx, yy in zip(x, y, strict=True)) / denominator
    return slope, ybar - slope * xbar


def _short_label(row: dict[str, Any]) -> str:
    short = str(row.get("mutation_short_label", ""))
    if short:
        return short.replace("-", "")
    match = re.search(r":([A-Z])(\d+)([A-Z])", str(row.get("mutation_signature", "")))
    return f"{match.group(1)}{match.group(2)}{match.group(3)}" if match else "point"


def _metric_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(row["target"], row["dataset"])].append(row)
    output = []
    for (target, dataset), values in sorted(grouped.items()):
        output.append({"target": target, "dataset": dataset, **_metric(values)})
    return output


def _write_notes(output: Path, rows: list[dict[str, Any]], metrics: list[dict[str, Any]]) -> None:
    targets = sorted({row["target"] for row in rows})
    contract = (
        "# Figure contract\n\n"
        "Core conclusion: current completed RBFE results are heterogeneous across targets; "
        "the figures separate formal benchmark evidence from stress, DSSB, and force-field experiments.\n"
        "Figure archetype: quantitative grid plus per-target scatter panels.\n"
        "Backend: Python/matplotlib only.\n"
        "Output: editable SVG/PDF plus 600-dpi TIFF and PNG previews.\n"
        "Panel map: one panel per target; x = experimental ΔΔG, y = predicted ΔΔG; "
        "dashed line = identity; colored markers = protocol family.\n"
        "Statistics: n, Pearson R, Spearman rho, MAE, RMSE, median absolute error, "
        "sign accuracy, and QC counts are exported in metrics_by_target.csv.\n"
        "Reviewer risk: small n, repeated/censored assay values, protocol mixing, and "
        "incomplete 3HFM/DSSB jobs are explicitly separated.\n"
        f"Targets plotted: {', '.join(targets)}.\n"
    )
    output.joinpath("figure_contract.md").write_text(contract, encoding="utf-8")
    qa_notes = (
        "# Figure QA notes\n\n"
        "- All plot rendering, preview export, and vector/raster export used Python/matplotlib.\n"
        "- SVG text is kept editable (`svg.fonttype=none`); PDF uses TrueType fonts.\n"
        "- Every point in `all_target_points.csv` has a completed `ddg_summary.json` and a resolved experimental reference.\n"
        "- Open markers identify duplicate experimental-reference rows (currently the 1DVF A:Y49A source ambiguity).\n"
        "- Active jobs were not plotted: 3HFM D32N, W98F, and Y50L.\n"
        "- Acceptance smoke roots without completed paired results were not treated as validation targets.\n"
        f"- Source rows: {len(rows)}; metric groups: {len(metrics)}.\n"
    )
    output.joinpath("qa_notes.md").write_text(qa_notes, encoding="utf-8")


def make_figures(output: Path, rows: list[dict[str, Any]]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "all_target_points.csv", rows)
    metrics = _metric_rows(rows)
    _write_csv(output / "metrics_by_target.csv", metrics)
    _write_notes(output, rows, metrics)

    by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_target[row["target"]].append(row)

    # Individual target figures are the primary user-facing deliverable.
    for target, values in sorted(by_target.items()):
        fig, ax = plt.subplots(figsize=(3.45, 3.25))
        _scatter_on_ax(ax, target, values)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="upper left", fontsize=5.6, handletextpad=0.4, borderpad=0.25)
        target_metrics = [row for row in metrics if row["target"] == target]
        lines = []
        for metric in target_metrics:
            lines.append(
                f"{metric['dataset']}: n={metric['n']}, R={_fmt(metric['pearson_r'])}, "
                f"MAE={_fmt(metric['mae_kcal_mol'])}"
            )
        ax.text(
            0.98,
            0.02,
            "\n".join(lines),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=5.4,
            color="#333333",
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "#D0D0D0", "alpha": 0.88},
        )
        fig.tight_layout(pad=0.6)
        _save_figure(fig, output / f"scatter_{target.lower()}")

    formal_targets = sorted(
        target
        for target, values in by_target.items()
        if any(row["dataset"] in {"development_fit", "independent_OOS", "validation_panel"} for row in values)
    )
    special_targets = sorted(set(by_target) - set(formal_targets))
    for name, targets in (("formal_targets_grid", formal_targets), ("special_targets_grid", special_targets)):
        if not targets:
            continue
        ncols = 2 if len(targets) <= 2 else 4
        nrows = math.ceil(len(targets) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(7.2, 2.35 * nrows), squeeze=False)
        axes_flat = list(axes.ravel())
        for index, target in enumerate(targets):
            ax = axes_flat[index]
            _scatter_on_ax(ax, target, by_target[target], compact=True)
            ax.tick_params(labelsize=5.7)
            ax.xaxis.label.set_size(6.2)
            ax.yaxis.label.set_size(6.2)
            ax.title.set_size(7.4)
            if index % ncols:
                ax.set_ylabel("")
            if index < (nrows - 1) * ncols:
                ax.set_xlabel("")
            metric_lines = []
            for metric in [row for row in metrics if row["target"] == target]:
                metric_lines.append(f"{metric['dataset']}: {metric['n']}, R={_fmt(metric['pearson_r'])}")
            ax.text(0.98, 0.02, "\n".join(metric_lines), transform=ax.transAxes, ha="right", va="bottom", fontsize=4.6)
        for ax in axes_flat[len(targets) :]:
            ax.axis("off")
        handles, labels = axes_flat[0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 0.002),
                ncol=min(4, len(labels)),
                fontsize=6.2,
                frameon=False,
            )
        fig.suptitle(
            "Current completed RBFE results by target" if name == "formal_targets_grid" else "Special and exploratory target results",
            fontsize=9,
            fontweight="bold",
            y=0.998,
        )
        fig.tight_layout(rect=(0, 0.045, 1, 0.965), h_pad=1.0, w_pad=0.8)
        _save_figure(fig, output / name)

    output.joinpath("metrics_by_target.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "reports/selections/validation_by_target_20260824",
    )
    args = parser.parse_args()
    rows = collect_rows()
    if not rows:
        raise SystemExit("No completed paired results found")
    make_figures(args.output, rows)
    print(json.dumps({"output": str(args.output), "rows": len(rows), "targets": sorted({row['target'] for row in rows})}, indent=2))


if __name__ == "__main__":
    main()
