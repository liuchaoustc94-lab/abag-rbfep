#!/usr/bin/env python3
"""Collect the latest completed RBFE pairs and render per-target scatter plots.

The collector intentionally reads completed ``results/ddg_summary.json`` files
from the current run roots instead of stale aggregate reports.  The official
view follows the 2026-08-31 target policy; excluded and external targets are
kept in a separate diagnostic view.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "reports/selections/validation_latest_20260908"

# Current official policy from AGENTS.md / docs/validation_targets.md.
OFFICIAL_TARGETS = (
    "1AHW",
    "1AK4",
    "1BJ1",
    "1CZ8",
    "1DQJ",
    "1DVF",
    "1IAR",
    "1KTZ",
    "1N8Z",
    "1OGA",
    "1VFB",
    "3BE1",
    "4I77",
)
EXCLUDED_TARGETS = ("1JRH", "1MLC", "1NMB", "2JEL", "3K2M", "3NGB")
SPECIAL_TARGETS = ("2BDN", "3HFM")

AB_BIND_SOURCE = ROOT / "benchmarks/ab_bind/curated/ab_bind_source_registered.csv"
PATEL_SOURCE = ROOT / "benchmarks/patel_2021_3hfm/experimental_ddg.csv"

ROOT_SPECS = (
    ("1DQJ", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_1dqj_core_v1"),
    ("1JRH", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_1jrh_core_v1"),
    ("1N8Z", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_1n8z_core_v1"),
    ("1VFB", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_1vfb_core_v1"),
    ("2JEL", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_2jel_core_v1"),
    ("3BE1", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1"),
    ("3NGB", "development_fit", ROOT / "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3ngb_core_v1"),
    ("1AK4", "independent_OOS", ROOT / "runs/real_cases/oos_1ak4_newprotocol_20260811"),
    ("1KTZ", "independent_OOS", ROOT / "runs/real_cases/oos_1ktz_newprotocol_20260811"),
    ("3K2M", "independent_OOS", ROOT / "runs/real_cases/oos_3k2m_newprotocol_20260811"),
    ("1BJ1", "validation_panel", ROOT / "runs/benchmarks/abbind_validation_panels_20260812/abbind_1bj1_core_v1"),
    ("1CZ8", "validation_panel", ROOT / "runs/benchmarks/abbind_validation_panels_20260812/abbind_1cz8_core_v1"),
    ("1DVF", "stress_1DVF", ROOT / "runs/real_cases/1dvf_priority_20260806"),
)

ACCEPTANCE_TARGETS = ("1AHW", "1IAR", "1NMB", "1OGA", "2BDN", "4I77")

DATASET_LABELS = {
    "development_fit": "development fit",
    "independent_OOS": "independent OOS",
    "validation_panel": "validation panel",
    "stress_1DVF": "1DVF stress",
    "acceptance_V1": "acceptance V1",
    "Patel_3HFM": "Patel 3HFM / DSSB",
}

DATASET_COLORS = {
    "development_fit": "#3B6EA8",
    "independent_OOS": "#1C4E80",
    "validation_panel": "#6E7781",
    "stress_1DVF": "#C27836",
    "acceptance_V1": "#3D8B73",
    "Patel_3HFM": "#A64646",
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "font.size": 7.0,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def norm_signature(value: str) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


def load_ab_bind_reference() -> dict[tuple[str, str], list[dict[str, str]]]:
    reference: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_csv(AB_BIND_SOURCE):
        reference[(row["complex_id"].upper(), norm_signature(row["signature"]))].append(row)
    return reference


def resolve_ab_bind(
    target: str,
    summary: dict[str, Any],
    reference: dict[tuple[str, str], list[dict[str, str]]],
) -> tuple[float | None, bool, str]:
    full = norm_signature(summary.get("mutation_signature", ""))
    base = full.split("@", 1)[0]
    side = str(summary.get("entity_side", "")).upper()
    candidates = [full, base + (f"@{side}" if side else ""), base]
    matches: list[dict[str, str]] = []
    for candidate in candidates:
        matches = reference.get((target.upper(), candidate), [])
        if matches:
            break
    if not matches:
        return None, False, ""
    values = [float(row["ddg_kcal_mol"]) for row in matches]
    # Preserve the existing 1DVF A:Y49A panel choice while exposing ambiguity.
    chosen = values[-1] if target.upper() == "1DVF" and base == "A:Y49A" else values[0]
    source_id = ";".join(row.get("row_id", "") for row in matches)
    return chosen, len(values) > 1, source_id


def load_acceptance_reference(target: str) -> dict[str, tuple[float, str]]:
    path = ROOT / "benchmarks/acceptance" / target.lower() / "experimental_ddg_v1.csv"
    if not path.exists():
        return {}
    return {
        row["mutation_group_id"]: (float(row["ddg_kcal_mol"]), row.get("note", ""))
        for row in read_csv(path)
    }


def qc_status(job_dir: Path) -> str:
    path = job_dir / "results/qc_report.json"
    if not path.exists():
        return "not_evaluated"
    try:
        value = str(json.loads(path.read_text(encoding="utf-8")).get("status", "not_evaluated"))
    except (OSError, json.JSONDecodeError):
        return "not_evaluated"
    return value if value in {"pass", "warning", "fail"} else "not_evaluated"


def latest_summaries(
    root: Path, *, include_result_variants: bool = False
) -> dict[str, tuple[Path, dict[str, Any]]]:
    latest: dict[str, tuple[Path, dict[str, Any]]] = {}
    pattern = "**/results*/ddg_summary.json" if include_result_variants else "**/results/ddg_summary.json"
    for path in sorted(root.glob(pattern)):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not payload.get("ready", True) or finite(payload.get("ddg_kcal_mol")) is None:
            continue
        job_id = str(payload.get("job_id", path.parents[1].name))
        stamp = str(payload.get("generated_at", ""))
        previous = latest.get(job_id)
        if previous is None or stamp >= str(previous[1].get("generated_at", "")):
            latest[job_id] = (path, payload)
    return latest


def make_row(
    *,
    target: str,
    dataset: str,
    path: Path,
    summary: dict[str, Any],
    experimental: float,
    note: str = "",
    ambiguous: bool = False,
    reference_id: str = "",
    provisional: bool = False,
) -> dict[str, Any]:
    job_dir = path.parents[1]
    return {
        "target": target.upper(),
        "dataset": dataset,
        "dataset_label": DATASET_LABELS.get(dataset, dataset),
        "job_id": str(summary.get("job_id", job_dir.name)),
        "mutation_group_id": str(summary.get("mutation_group_id", "")),
        "mutation_signature": str(summary.get("mutation_signature", "")),
        "mutation_short_label": str(summary.get("mutation_short_label", "")),
        "experimental_ddg_kcal_mol": experimental,
        "predicted_ddg_kcal_mol": float(summary["ddg_kcal_mol"]),
        "repeat_sd_kcal_mol": finite(summary.get("ddg_repeat_stdev_kcal_mol")),
        "bar_stderr_kcal_mol": finite(summary.get("ddg_bar_stderr_kcal_mol")),
        "qc_status": qc_status(job_dir),
        "experimental_note": note,
        "experimental_reference_id": reference_id,
        "experimental_reference_ambiguous": ambiguous,
        "provisional": provisional,
        "generated_at": str(summary.get("generated_at", "")),
        "source_tag": str(path.relative_to(ROOT)),
    }


def collect_ab_bind_rows() -> list[dict[str, Any]]:
    reference = load_ab_bind_reference()
    rows: list[dict[str, Any]] = []
    for target, dataset, root in ROOT_SPECS:
        for job_id, (path, summary) in latest_summaries(root).items():
            experimental, ambiguous, reference_id = resolve_ab_bind(target, summary, reference)
            if experimental is None:
                continue
            rows.append(
                make_row(
                    target=target,
                    dataset=dataset,
                    path=path,
                    summary=summary,
                    experimental=experimental,
                    ambiguous=ambiguous,
                    reference_id=reference_id,
                )
            )
    return rows


def collect_acceptance_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target in ACCEPTANCE_TARGETS:
        reference = load_acceptance_reference(target)
        root = ROOT / "runs/acceptance" / f"acc_{target.lower()}_v1_20260818"
        for job_id, (path, summary) in latest_summaries(root).items():
            mutation_group_id = str(summary.get("mutation_group_id", ""))
            if mutation_group_id not in reference:
                continue
            experimental, note = reference[mutation_group_id]
            rows.append(
                make_row(
                    target=target,
                    dataset="acceptance_V1",
                    path=path,
                    summary=summary,
                    experimental=experimental,
                    note=note,
                    reference_id=mutation_group_id,
                )
            )
    return rows


def collect_patel_rows() -> list[dict[str, Any]]:
    reference = {row["job_id"]: row for row in read_csv(PATEL_SOURCE)}
    roots = (
        ROOT / "runs/real_cases/3hfm_patel_align_20260817",
        ROOT / "runs/real_cases/3hfm_patel_dssb_pr4",
    )
    latest: dict[str, tuple[Path, dict[str, Any]]] = {}
    for root in roots:
        for job_id, item in latest_summaries(root, include_result_variants=True).items():
            if job_id not in reference:
                continue
            previous = latest.get(job_id)
            if previous is None or str(item[1].get("generated_at", "")) >= str(previous[1].get("generated_at", "")):
                latest[job_id] = item
    rows = []
    for job_id, (path, summary) in latest.items():
        ref = reference[job_id]
        rows.append(
            make_row(
                target="3HFM",
                dataset="Patel_3HFM",
                path=path,
                summary=summary,
                experimental=float(ref["experimental_ddg_kcal_mol"]),
                note="Patel 2021 external reference; DSSB/chaining remains provisional",
                reference_id="Patel 2021 Table 2",
                provisional=True,
            )
        )
    return rows


def mark_censored(rows: list[dict[str, Any]]) -> None:
    by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_target[row["target"]].append(row)
    for target_rows in by_target.values():
        values = [row["experimental_ddg_kcal_mol"] for row in target_rows]
        counts: dict[float, int] = defaultdict(int)
        for value in values:
            counts[round(value, 8)] += 1
        lo, hi = min(values), max(values)
        for row in target_rows:
            value = row["experimental_ddg_kcal_mol"]
            row["experimental_censored_like"] = bool(
                counts[round(value, 8)] >= 2 and (value == lo or value == hi)
            )


def metric(rows: list[dict[str, Any]]) -> dict[str, Any]:
    x = np.asarray([row["experimental_ddg_kcal_mol"] for row in rows], dtype=float)
    y = np.asarray([row["predicted_ddg_kcal_mol"] for row in rows], dtype=float)
    result: dict[str, Any] = {
        "n": int(len(rows)),
        "pearson_r": None,
        "spearman_rho": None,
        "mae_kcal_mol": None,
        "rmse_kcal_mol": None,
        "median_abs_error_kcal_mol": None,
        "qc_pass": sum(row["qc_status"] == "pass" for row in rows),
        "qc_warning": sum(row["qc_status"] == "warning" for row in rows),
        "qc_fail": sum(row["qc_status"] == "fail" for row in rows),
        "censored_like_n": sum(bool(row.get("experimental_censored_like")) for row in rows),
    }
    if len(rows):
        error = y - x
        result.update(
            {
                "mae_kcal_mol": float(np.mean(np.abs(error))),
                "rmse_kcal_mol": float(np.sqrt(np.mean(error**2))),
                "median_abs_error_kcal_mol": float(np.median(np.abs(error))),
            }
        )
    if len(rows) >= 3 and np.std(x) > 0 and np.std(y) > 0:
        result["pearson_r"] = float(pearsonr(x, y).statistic)
        result["spearman_rho"] = float(spearmanr(x, y).statistic)
    return result


def fmt(value: Any) -> str:
    return "—" if value is None else f"{float(value):.2f}"


def short_label(row: dict[str, Any]) -> str:
    label = str(row.get("mutation_short_label", ""))
    if label:
        return label.replace("-", "")
    match = re.search(r":([A-Z])(\d+)([A-Z])", str(row.get("mutation_signature", "")))
    return "".join(match.groups()) if match else "point"


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def limits(rows: list[dict[str, Any]]) -> tuple[float, float]:
    values = [row["experimental_ddg_kcal_mol"] for row in rows] + [row["predicted_ddg_kcal_mol"] for row in rows]
    lo, hi = min(values), max(values)
    span = max(hi - lo, 1.0)
    pad = max(0.35, 0.07 * span)
    return math.floor((lo - pad) * 2) / 2, math.ceil((hi + pad) * 2) / 2


def draw_panel(ax: plt.Axes, target: str, rows: list[dict[str, Any]], *, labels: bool = False) -> None:
    lo, hi = limits(rows)
    ax.plot([lo, hi], [lo, hi], color="#8A8A8A", lw=0.75, ls=(0, (3, 2)), zorder=0)
    for row in rows:
        x = row["experimental_ddg_kcal_mol"]
        y = row["predicted_ddg_kcal_mol"]
        dataset = row["dataset"]
        color = DATASET_COLORS.get(dataset, "#4C78A8")
        status = row["qc_status"]
        marker = "x" if status == "fail" else "o"
        face = "none" if row.get("experimental_censored_like") or status == "warning" else color
        edge = "#B6483E" if status == "fail" else color
        if status == "fail":
            ax.scatter(x, y, s=30, marker="x", color=edge, linewidths=1.0, alpha=0.92, zorder=3)
        else:
            ax.scatter(x, y, s=25, marker="o", facecolors=face, edgecolors=edge, linewidths=0.8, alpha=0.88, zorder=3)
        if labels:
            ax.annotate(short_label(row), (x, y), xytext=(3, 3), textcoords="offset points", fontsize=5.4, color="#303030")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Experimental ΔΔG (kcal mol⁻¹)")
    ax.set_ylabel("Predicted ΔΔG (kcal mol⁻¹)")
    m = metric(rows)
    ax.set_title(f"{target}  n={m['n']}  R={fmt(m['pearson_r'])}  ρ={fmt(m['spearman_rho'])}", fontsize=8, fontweight="bold", pad=4)
    ax.grid(False)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_outputs(rows: list[dict[str, Any]]) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mark_censored(rows)
    by_target: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_target[row["target"]].append(row)

    write_csv(OUTPUT / "latest_completed_points.csv", rows)
    metrics: list[dict[str, Any]] = []
    for target in sorted(by_target):
        metrics.append({"target": target, "view": "official" if target in OFFICIAL_TARGETS else "diagnostic", **metric(by_target[target])})
    official_rows = [row for row in rows if row["target"] in OFFICIAL_TARGETS]
    metrics.append({"target": "POOLED_OFFICIAL", "view": "official", **metric(official_rows)})
    write_csv(OUTPUT / "latest_metrics_by_target.csv", metrics)
    (OUTPUT / "latest_metrics_by_target.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    official_targets = [target for target in OFFICIAL_TARGETS if target in by_target]
    ncols = 4
    nrows = math.ceil(len(official_targets) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(10.2, 2.75 * nrows), squeeze=False)
    for index, target in enumerate(official_targets):
        draw_panel(axes.ravel()[index], target, by_target[target])
    for ax in axes.ravel()[len(official_targets) :]:
        ax.axis("off")
    fig.suptitle("Latest completed official RBFE results by target", fontsize=11, fontweight="bold", y=0.998)
    fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=1.1, w_pad=0.8)
    save_figure(fig, OUTPUT / "official_target_grid")

    diagnostic_targets = [target for target in sorted(by_target) if target not in OFFICIAL_TARGETS]
    if diagnostic_targets:
        ncols = 3
        nrows = math.ceil(len(diagnostic_targets) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(8.0, 2.9 * nrows), squeeze=False)
        for index, target in enumerate(diagnostic_targets):
            draw_panel(axes.ravel()[index], target, by_target[target])
        for ax in axes.ravel()[len(diagnostic_targets) :]:
            ax.axis("off")
        fig.suptitle("Excluded and exploratory RBFE results (diagnostic only)", fontsize=10.5, fontweight="bold", y=0.998)
        fig.tight_layout(rect=(0, 0, 1, 0.985), h_pad=1.1, w_pad=0.8)
        save_figure(fig, OUTPUT / "diagnostic_target_grid")

    for target in sorted(by_target):
        fig, ax = plt.subplots(figsize=(3.55, 3.35))
        draw_panel(ax, target, by_target[target], labels=True)
        fig.tight_layout(pad=0.7)
        save_figure(fig, OUTPUT / f"scatter_{target.lower()}")

    pooled = official_rows
    fig, ax = plt.subplots(figsize=(4.0, 3.7))
    lo, hi = limits(pooled)
    ax.plot([lo, hi], [lo, hi], color="#8A8A8A", lw=0.8, ls=(0, (3, 2)))
    for target in official_targets:
        group = by_target[target]
        color = DATASET_COLORS.get(group[0]["dataset"], "#4C78A8")
        ax.scatter([r["experimental_ddg_kcal_mol"] for r in group], [r["predicted_ddg_kcal_mol"] for r in group], s=18, alpha=0.55, color=color, label=target)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Experimental ΔΔG (kcal mol⁻¹)")
    ax.set_ylabel("Predicted ΔΔG (kcal mol⁻¹)")
    m = metric(pooled)
    ax.set_title(f"Official pooled view  n={m['n']}  R={fmt(m['pearson_r'])}  ρ={fmt(m['spearman_rho'])}", fontsize=8.5, fontweight="bold")
    ax.legend(loc="upper left", fontsize=5.5, ncol=2, frameon=False, handletextpad=0.3, columnspacing=0.8)
    fig.tight_layout(pad=0.7)
    save_figure(fig, OUTPUT / "official_pooled_scatter")

    contract = """# Figure contract

Core conclusion: completed RBFE performance is heterogeneous across targets; the official view and diagnostic/external views must remain separate.
Archetype: quantitative grid plus one per-target scatter panel.
Backend: Python/matplotlib only.
Axes: x = experimental ΔΔG; y = predicted ΔΔG; dashed line = identity.
Official policy: 2026-08-31 acceptance set, excluding 2JEL, 3NGB, 3K2M, 1JRH, 1NMB and 1MLC.
Latest rule: only completed canonical result summaries available at render time are plotted; active P1 rollout jobs are excluded.
Statistics: Pearson R and Spearman rho are primary; MAE/RMSE, QC counts and censored-like flags are exported as supporting diagnostics.
"""
    (OUTPUT / "figure_contract.md").write_text(contract, encoding="utf-8")
    (OUTPUT / "qa_notes.md").write_text(
        """# Figure QA notes

- All rendering and preview/export used Python/matplotlib.
- SVG text remains editable and PDF uses TrueType fonts.
- Open markers indicate warning QC or repeated endpoint-like experimental values; red x marks QC failure.
- Official metrics are computed from the exact points in latest_completed_points.csv; no calibration or outlier filtering is applied.
- Active jobs from the 2026-09-08 P1 rollout and incomplete acceptance jobs are not treated as completed results.
""",
        encoding="utf-8",
    )

    pooled_metric = metric(official_rows)
    print(json.dumps({"output": str(OUTPUT), "rows": len(rows), "official_rows": len(official_rows), "targets": sorted(by_target), "pooled_official": pooled_metric}, indent=2))


def main() -> None:
    rows = collect_ab_bind_rows() + collect_acceptance_rows() + collect_patel_rows()
    if not rows:
        raise SystemExit("No completed paired results found")
    make_outputs(rows)


if __name__ == "__main__":
    main()
