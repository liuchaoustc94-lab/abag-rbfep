#!/usr/bin/env python3
"""Summarize SKEMPI2 affinity changes for a materialized target panel."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


def _median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def _fmt(value: float | None, digits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}g}"


def _read_unique_single_rows(path: Path) -> list[dict[str, str]]:
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        if row.get("single_point", "").strip().lower() != "true":
            continue
        if row.get("ddg_censored", "").strip().lower() == "true":
            continue
        signature = row.get("canonical_signature", "").strip()
        if not signature or signature in seen:
            continue
        try:
            wt = float(row["affinity_wt_M"])
            mut = float(row["affinity_mut_M"])
        except (KeyError, TypeError, ValueError):
            continue
        if wt <= 0 or mut <= 0:
            continue
        seen.add(signature)
        row["_wt"] = str(wt)
        row["_mut"] = str(mut)
        selected.append(row)
    return selected


def build(panel: Path) -> tuple[Path, Path, Path]:
    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    manifest = {
        row["pdb_id"]: row
        for row in csv.DictReader(
            (panel / "target_manifest.csv").open(encoding="utf-8", newline="")
        )
    }

    for target_dir in sorted((panel / "targets").iterdir()):
        if not target_dir.is_dir():
            continue
        rows = _read_unique_single_rows(target_dir / "experimental_ddg.csv")
        if not rows:
            continue
        target_id = target_dir.name
        wt_values = [float(row["_wt"]) for row in rows]
        mut_values = [float(row["_mut"]) for row in rows]
        kd_folds = [mut / wt for mut, wt in zip(mut_values, wt_values)]
        affinity_folds = [1.0 / fold for fold in kd_folds]
        delta_pkd = [-math.log10(mut / wt) for mut, wt in zip(mut_values, wt_values)]
        ddg = [float(row["ddg_kcal_mol"]) for row in rows]
        gain = sum(value > 0 for value in delta_pkd)
        loss = sum(value < 0 for value in delta_pkd)
        neutral = len(delta_pkd) - gain - loss
        for row, kd_fold, affinity_fold, dpka in zip(
            rows, kd_folds, affinity_folds, delta_pkd
        ):
            detail_rows.append(
                {
                    "pdb_id": target_id,
                    "target_name": manifest.get(target_id, {}).get("name", ""),
                    "source_row_id": row["source_row_id"],
                    "mutation_tokens": row["mutation_tokens"],
                    "canonical_signature": row["canonical_signature"],
                    "temperature_K": row["temperature_K"],
                    "affinity_wt_M": _fmt(float(row["_wt"]), 8),
                    "affinity_mut_M": _fmt(float(row["_mut"]), 8),
                    "kd_fold_mut_over_wt": _fmt(kd_fold, 8),
                    "affinity_fold_mut_over_wt": _fmt(affinity_fold, 8),
                    "delta_pKd": _fmt(dpka, 8),
                    "ddg_kcal_mol": row["ddg_kcal_mol"],
                    "reference": row["reference"],
                    "method": row["method"],
                }
            )
        summary_rows.append(
            {
                "pdb_id": target_id,
                "target_name": manifest.get(target_id, {}).get("name", ""),
                "n_unique_single": len(rows),
                "n_affinity_gain": gain,
                "n_affinity_loss": loss,
                "n_affinity_neutral": neutral,
                "wt_kd_median_M": _fmt(_median(wt_values), 8),
                "wt_pKd_median": _fmt(_median([-math.log10(v) for v in wt_values]), 6),
                "mut_kd_median_M": _fmt(_median(mut_values), 8),
                "mut_kd_min_M": _fmt(min(mut_values), 8),
                "mut_kd_max_M": _fmt(max(mut_values), 8),
                "kd_fold_median_mut_over_wt": _fmt(_median(kd_folds), 8),
                "kd_fold_min_mut_over_wt": _fmt(min(kd_folds), 8),
                "kd_fold_max_mut_over_wt": _fmt(max(kd_folds), 8),
                "affinity_fold_median_mut_over_wt": _fmt(_median(affinity_folds), 8),
                "delta_pKd_median": _fmt(_median(delta_pkd), 6),
                "delta_pKd_mean": _fmt(statistics.mean(delta_pkd), 6),
                "delta_pKd_min": _fmt(min(delta_pkd), 6),
                "delta_pKd_max": _fmt(max(delta_pkd), 6),
                "ddg_median_kcal_mol": _fmt(_median(ddg), 6),
                "ddg_mean_kcal_mol": _fmt(statistics.mean(ddg), 6),
                "ddg_min_kcal_mol": _fmt(min(ddg), 6),
                "ddg_max_kcal_mol": _fmt(max(ddg), 6),
            }
        )

    detail_path = panel / "affinity_change_single_unique.csv"
    summary_path = panel / "affinity_change_summary.csv"
    report_path = panel / "affinity_change_report.md"
    detail_fields = list(detail_rows[0].keys()) if detail_rows else []
    summary_fields = list(summary_rows[0].keys()) if summary_rows else []
    for path, rows, fields in (
        (detail_path, detail_rows, detail_fields),
        (summary_path, summary_rows, summary_fields),
    ):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# SKEMPI2 affinity-change summary\n\n")
        handle.write(
            "Statistics use the first occurrence of each unique single-point "
            "signature. `Kd_mut/Kd_wt > 1` means weaker binding; positive "
            "`delta_pKd` means stronger binding; positive `ddG` means weaker "
            "binding. Values are descriptive and retain source-specific WT "
            "baselines.\n\n"
        )
        handle.write(
            "| PDB | n | WT Kd (M) | mutant Kd median (M) | Kd fold median | "
            "affinity fold median | median ΔpKd | median ΔΔG (kcal/mol) | gain/loss |\n"
        )
        handle.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in summary_rows:
            handle.write(
                f"| {row['pdb_id']} | {row['n_unique_single']} | "
                f"{row['wt_kd_median_M']} | {row['mut_kd_median_M']} | "
                f"{row['kd_fold_median_mut_over_wt']} | "
                f"{row['affinity_fold_median_mut_over_wt']} | "
                f"{row['delta_pKd_median']} | {row['ddg_median_kcal_mol']} | "
                f"{row['n_affinity_gain']}/{row['n_affinity_loss']} |\n"
            )
        handle.write(
            "\nThe per-mutation values, ranges, source references, and assay methods "
            "are in `affinity_change_single_unique.csv`; the full summary table "
            "is in `affinity_change_summary.csv`.\n"
        )
    return summary_path, detail_path, report_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        type=Path,
        default=Path("benchmarks/protein_protein_skempi2_single_point"),
    )
    args = parser.parse_args()
    paths = build(args.panel)
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
