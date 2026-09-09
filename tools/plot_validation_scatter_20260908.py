#!/usr/bin/env python
"""Validation scatter report (2026-09-08): per-target + pooled pred-vs-exp.

- Main set: standard two-leg protocol, judged evaluation set (13 targets;
  2JEL/3NGB/3K2M/1JRH/1NMB/1MLC excluded per user decisions)
- Seeded P1 pilot points (w50a/y50a/w52a) shown as green arrows old->new.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/mnt/data/liuchao/abag-rbfep")
OUT = ROOT / "runs/analysis/validation_scatter_20260908"
OUT.mkdir(parents=True, exist_ok=True)

EXCLUDED = {"2jel", "3ngb", "3k2m", "1jrh", "1nmb", "1mlc", "3hfm"}
SEEDED = {  # (complex, mut): new seeded value
    ("3be1", "w50a"): 5.93,
    ("1dqj", "y50a"): 5.03,
    ("1vfb", "w52a"): -1.88,
}

rows = []
for r in csv.DictReader(open(ROOT / "runs/analysis/protonation_reweight_20260831/pairs_protonation.csv")):
    if r["complex"] in EXCLUDED:
        continue
    rows.append({"cx": r["complex"], "mut": r["mut"],
                 "pred": float(r["pred"]), "exp": float(r["exp"])})


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    vx = sum((a - mx) ** 2 for a in x)
    vy = sum((b - my) ** 2 for b in y)
    vp = vx * vy
    return cov / math.sqrt(vp) if vp > 0 else float("nan")


def spearman(x, y):
    def rk(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    return pearson(rk(x), rk(y))


by = defaultdict(list)
for r in rows:
    by[r["cx"]].append(r)

targets = sorted(by, key=lambda c: -pearson([p["pred"] for p in by[c]], [p["exp"] for p in by[c]]))

# ---------- figure ----------
fig = plt.figure(figsize=(16, 13))
gs = fig.add_gridspec(4, 5, hspace=0.6, wspace=0.35)

# pooled panel (left, spans 2x2)
ax = fig.add_subplot(gs[0:2, 0:2])
allp = [r["pred"] for r in rows]
alle = [r["exp"] for r in rows]
R = pearson(allp, alle)
Sp = spearman(allp, alle)
mae = sum(abs(a - b) for a, b in zip(allp, alle)) / len(allp)
cmap = plt.get_cmap("tab20")
colors = {cx: cmap(i % 20) for i, cx in enumerate(sorted(by))}
for cx in sorted(by):
    v = by[cx]
    ax.scatter([p["exp"] for p in v], [p["pred"] for p in v], s=22, alpha=0.8,
               color=colors[cx], label=cx.upper(), edgecolors="none")
lo = min(min(alle), min(allp)) - 1
hi = max(max(alle), max(allp)) + 1
ax.plot([lo, hi], [lo, hi], "k-", lw=0.8)
ax.plot([lo, hi], [lo + 1, hi + 1], "k:", lw=0.6)
ax.plot([lo, hi], [lo - 1, hi - 1], "k:", lw=0.6)
ax.set_xlabel("Experimental ΔΔG (kcal/mol)")
ax.set_ylabel("Predicted ΔΔG (kcal/mol)")
ax.set_title(f"Pooled (judged set, n={len(rows)})\nR={R:.2f}  Sp={Sp:.2f}  MAE={mae:.2f}", fontsize=10)
ax.legend(fontsize=6, ncol=2, frameon=False, loc="upper left")

# per-target mini panels
mini_positions = [(0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4),
                  (2, 0), (2, 1), (2, 2), (2, 3), (2, 4),
                  (3, 0), (3, 1), (3, 2), (3, 3), (3, 4)]
for (rr, cc), cx in zip(mini_positions, targets[: len(mini_positions)]):
    axm = fig.add_subplot(gs[rr, cc])
    v = by[cx]
    p = [x["pred"] for x in v]
    e = [x["exp"] for x in v]
    axm.scatter(e, p, s=16, color=colors[cx], edgecolors="none", alpha=0.85)
    lo = min(min(e), min(p)) - 0.5
    hi = max(max(e), max(p)) + 0.5
    axm.plot([lo, hi], [lo, hi], "k-", lw=0.6)
    r = pearson(p, e)
    s = spearman(p, e)
    m = sum(abs(a - b) for a, b in zip(p, e)) / len(p)
    axm.set_title(f"{cx.upper()}  n={len(v)}\nR={r:.2f} Sp={s:.2f} MAE={m:.1f}", fontsize=8)
    axm.tick_params(labelsize=6)
    # seeded pilot arrows
    for x in v:
        key = (cx, x["mut"])
        if key in SEEDED:
            axm.annotate("", xy=(x["exp"], SEEDED[key]), xytext=(x["exp"], x["pred"]),
                         arrowprops=dict(arrowstyle="->", color="green", lw=1.4))
            axm.scatter([x["exp"]], [SEEDED[key]], color="green", s=28, marker="*", zorder=5)

# remaining targets note
if len(targets) > len(mini_positions):
    fig.text(0.02, 0.02, f"(+{len(targets)-len(mini_positions)} more targets not shown)", fontsize=7)

fig.text(0.02, 0.965, "abag-rbfep validation scatter — judged set (13 targets, standard two-leg protocol)",
         fontsize=12, weight="bold")
fig.text(0.02, 0.945, "green star/arrow: P1 hydration-seeded pilot points (w50a 9.73→5.93, y50a 9.95→5.03, w52a 4.58→−1.88)",
         fontsize=8, color="green")
fig.savefig(OUT / "validation_scatter_20260908.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "validation_scatter_20260908.pdf", bbox_inches="tight")
print("saved", OUT)

# ---------- summary table ----------
print(f"\n{'靶点':<6}{'n':>4}{'R':>8}{'Sp':>8}{'MAE':>7}  评级")
import statistics
rs = []
for cx in targets:
    v = by[cx]
    p = [x["pred"] for x in v]; e = [x["exp"] for x in v]
    r = pearson(p, e); s = spearman(p, e)
    m = sum(abs(a - b) for a, b in zip(p, e)) / len(p)
    rs.append(r)
    grade = "✅强" if r >= 0.7 else ("✅良好" if r >= 0.5 else ("⚠️中" if r >= 0.4 else "❌弱"))
    print(f"{cx.upper():<6}{len(v):>4}{r:>8.3f}{s:>8.3f}{m:>7.2f}  {grade}")
print(f"\npooled n={len(rows)}: R={R:.3f} Sp={Sp:.3f} MAE={mae:.2f}")
print(f"靶点中位 R={statistics.median(rs):.3f}")
