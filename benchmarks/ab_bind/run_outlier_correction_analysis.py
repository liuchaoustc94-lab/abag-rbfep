#!/usr/bin/env python3
"""Sampson-style outlier classification + LOCO-CV empirical correction driver.

Assembles all fixed-chemistry validation points (validation panels + out-of-sample
+ acceptance), classifies each mutation with structure-informed features, fits the
single shrinkage parameter with leave-one-complex-out CV, and reports
before/after metrics.
"""
from __future__ import annotations

import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "/mnt/data/liuchao/abag-rbfep/src")
from abag_rbfe.outliers import apply_outlier_correction, classify_outlier_features, fit_shrinkage_alpha, outlier_class_label

ROOT = Path("/mnt/data/liuchao/abag-rbfep")
OUT = ROOT / "runs/analysis/target_rescue_20260805/outlier_correction"
OUT.mkdir(parents=True, exist_ok=True)

STRUCT = ROOT / "benchmarks/ab_bind/source/structures"
ACC_STRUCT = ROOT / "benchmarks/acceptance/structures"


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    return cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else float("nan")


def ba3(pred, exp, t=1.0):
    def c(v):
        return 0 if v <= -t else (2 if v >= t else 1)
    recalls = []
    for k in (0, 1, 2):
        m = [(p, e) for p, e in zip(pred, exp) if c(e) == k]
        if m:
            recalls.append(sum(1 for p, _ in m if c(p) == k) / len(m))
    return sum(recalls) / len(recalls) if recalls else float("nan")


def metrics(pred, exp):
    return {
        "n": len(pred),
        "pearson_r": round(pearson(pred, exp), 3),
        "mae": round(sum(abs(a - b) for a, b in zip(pred, exp)) / len(pred), 3),
        "ba3": round(ba3(pred, exp), 3),
        "strong_mae": round(
            sum(abs(a - b) for a, b in zip(pred, exp) if abs(b) >= 2.5) / max(1, sum(1 for b in exp if abs(b) >= 2.5)), 3
        ),
    }


MUT_TOK = re.compile(r"([a-z])-([a-z])(\d+)([a-z])$")


def parse_job_mutation(job_id: str) -> tuple[str, int, str, str] | None:
    # e.g. 1bj1-antigen-w-q89a / 1dvf-d13-e52-antigen-d-y102a / acc-4i77-antibody-h-y100a
    m = MUT_TOK.search(job_id)
    if not m:
        return None
    chain, wt, resseq, mut = m.groups()
    return chain.upper(), int(resseq), wt.upper(), mut.upper()


def structure_for(complex_id: str) -> Path | None:
    for base in (STRUCT, ACC_STRUCT):
        p = base / f"{complex_id.upper()}.pdb"
        if p.is_file():
            return p
    return None


def main() -> int:
    # 1) 汇总全部点（job, complex, pred, exp）
    points = []
    for r in json.load(open(ROOT / "runs/analysis/target_rescue_20260805/validation_panels_newprotocol_pairs.json")):
        points.append({"complex": r.get("complex", r["job"][:4].upper()), "job": r["job"], "pred": r["pred"], "exp": r["exp"], "src": "validation_panels"})
    for case, pattern in [
        ("1DVF", "runs/real_cases/1dvf_priority_20260806/jobs/*/results/ddg_summary.json"),
        ("1AK4", "runs/real_cases/oos_1ak4_newprotocol_20260811/jobs/*/results/ddg_summary.json"),
        ("1KTZ", "runs/real_cases/oos_1ktz_newprotocol_20260811/jobs/*/results/ddg_summary.json"),
        ("3K2M", "runs/real_cases/oos_3k2m_newprotocol_20260811/jobs/*/results/ddg_summary.json"),
    ]:
        exp_map = {
            r["mutation_group_id"].replace("_", "-"): float(r["ddg_kcal_mol"])
            for r in csv.DictReader(open(ROOT / f"examples/real_cases/{case.lower()}/experimental_ddg.csv"))
        }
        for f in Path(ROOT).glob(pattern):
            j = f.parts[-3]
            suffix = j.split("-out-of-sample-")[-1] if "-out-of-sample-" in j else j.replace(f"{case.lower()}-d13-e52-", "")
            d = json.loads(f.read_text())
            if d.get("ready") and suffix in exp_map:
                points.append({"complex": case, "job": j, "pred": d["ddg_kcal_mol"], "exp": exp_map[suffix], "src": "oos"})
    # 验收 4I77
    exp77 = {r["mutation_group_id"].replace("_", "-"): float(r["ddg_kcal_mol"]) for r in csv.DictReader(open(ROOT / "benchmarks/acceptance/4i77/experimental_ddg_v1.csv"))}
    for f in Path(ROOT).glob("runs/acceptance/acc_4i77_v1_20260818/jobs/*/results/ddg_summary.json"):
        j = f.parts[-3]
        d = json.loads(f.read_text())
        suffix = j.replace("acc-4i77-", "")
        if d.get("ready") and suffix in exp77:
            points.append({"complex": "4I77", "job": j, "pred": d["ddg_kcal_mol"], "exp": exp77[suffix], "src": "acc"})

    print(f"assembled {len(points)} points")

    # 2) 逐点特征
    for pt in points:
        mut = parse_job_mutation(pt["job"])
        if mut is None:
            pt["features"] = {}
            continue
        chain, resseq, wt, mut_aa = mut
        pt["mut"] = {"chain": chain, "resseq": resseq, "wt": wt, "mut": mut_aa}
        pt["features"] = classify_outlier_features(
            wt=wt, mut=mut_aa, structure_path=structure_for(pt["complex"]), chain_id=chain, resseq=resseq
        )
        pt["outlier_class"] = outlier_class_label(pt["features"])

    # 类别统计
    class_counts = defaultdict(list)
    for pt in points:
        class_counts[pt.get("outlier_class", "")].append(pt)
    print("\n=== outlier 分类 ===")
    for cls, pts in sorted(class_counts.items(), key=lambda kv: -len(kv[1])):
        mae = sum(abs(p["pred"] - p["exp"]) for p in pts) / len(pts)
        print(f"  [{cls or 'none'}] n={len(pts)} MAE={mae:.2f}")

    # 3) LOCO-CV：主类 buried_aromatic_large_deletion 的收缩因子
    flags_all = [pt.get("outlier_class") == "buried_aromatic_large_deletion" for pt in points]
    print(f"\nburied_aromatic_large_deletion 点数: {sum(flags_all)}")
    complexes = sorted({pt["complex"] for pt in points})
    corrected = []
    for pt, flag in zip(points, flags_all):
        if not flag:
            corrected.append(pt["pred"])
            continue
        fit_idx = [i for i, q in enumerate(points) if q["complex"] != pt["complex"]]
        alpha = fit_shrinkage_alpha(
            [points[i]["pred"] for i in fit_idx],
            [points[i]["exp"] for i in fit_idx],
            [flags_all[i] for i in fit_idx],
        )
        corrected.append(alpha * pt["pred"])

    pred_raw = [pt["pred"] for pt in points]
    exp = [pt["exp"] for pt in points]
    result = {
        "raw": metrics(pred_raw, exp),
        "corrected_loco": metrics(corrected, exp),
        "class_counts": {cls: {"n": len(pts), "mae": round(sum(abs(p["pred"] - p["exp"]) for p in pts) / len(pts), 3)} for cls, pts in class_counts.items()},
        "points": points,
    }
    print("\n=== 指标对比（LOCO-CV）===")
    print("raw      :", result["raw"])
    print("corrected:", result["corrected_loco"])

    (OUT / "outlier_correction_loco.json").write_text(json.dumps(result, indent=1, default=str))
    print(f"\nsaved {OUT / 'outlier_correction_loco.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
