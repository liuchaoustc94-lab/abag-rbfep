"""Mutation-type → protocol-path routing (user policy table 2026-08-26).

Classifies each mutation group into a primary protocol route plus advisory
tags. Routes are advisory labels carried into batch plans (jobs.csv /
routing_summary.json); only paths with existing infrastructure are additionally
enforced at the protocol level (charge-changing → DSSB, already wired in
planning.build_batch_plan).

Policy table (user decision 2026-08-26):

| 突变类型 | 推荐路径 |
|---|---|
| 普通中性 Ala 扫描 | 12λ 基线 |
| Pro/Gly、侧链插入、柔性环区 | 端点系综 + 局部 REST2/RID |
| Y/W/M 大空腔删除 | RID + 水占据分析；必要时 GCMC |
| 电荷变化 | DSSB；失败时评估 co-alchemical ion |
| His/Asp/Glu/Lys 邻域 | PROPKA 多质子化状态 + pKa 加权 |
| 高重复离散 | 增加 basin，而非单纯增加 λ 窗口 |

Notes:
- "Y/W/M 大空腔删除" is implemented as wt ∈ {Y,W,F,M} → mut ∈ {A,G}: the
  user's Y/W/M plus F (our measured P1 error class is Y/W/F→A; M is the
  large hydrophobe the user added). Deviations are flagged in the summary.
- "柔性环区" cannot be derived from the mutation alone; loop/CDR membership
  should come from system-level annotation (hook left open, see
  LOOP_ANNOTATION_HOOK).
- "高重复离散 → 增加 basin" is a QC-time advisory, not a mutation route; it is
  documented here for completeness and surfaces as the `basin_diversity` tag
  when a mutation is large-deletion class (where repeat spread is likeliest).
"""

from __future__ import annotations

from dataclasses import dataclass

# Sidechain heavy-atom counts (mirrors stages._SIDECHAIN_HEAVY_ATOMS; kept local
# to avoid a planning -> stages import cycle).
_SIDECHAIN_HEAVY_ATOMS = {
    "G": 0, "A": 1, "S": 2, "C": 2, "T": 3, "V": 3, "P": 3,
    "N": 4, "D": 4, "L": 4, "I": 4, "M": 4, "K": 5, "Q": 5,
    "E": 5, "H": 6, "F": 7, "R": 7, "Y": 8, "W": 10,
}

_LARGE_CAVITY_RESIDUES = frozenset({"Y", "W", "F", "M"})
_SMALL_DELETION_TARGETS = frozenset({"A", "G"})
_TITRATABLE_RESIDUES = frozenset({"H", "D", "E", "K"})

ROUTE_BASELINE = "baseline_12lambda"
ROUTE_ENDPOINT_ENSEMBLE = "endpoint_ensemble_rest"  # 端点系综 + 局部 REST2/RID
ROUTE_RID_HYDRATION = "rid_hydration"  # RID + 水占据分析（必要时 GCMC）
ROUTE_DSSB = "dssb"  # 电荷变化双盒单盒
ROUTE_PROTONATION = "protonation_multi"  # PROPKA 多态 + pKa 加权

#: Loop/CDR membership must come from system annotation; placeholder hook.
LOOP_ANNOTATION_HOOK = "loop_membership_requires_system_annotation"


@dataclass(frozen=True)
class MutationRoute:
    primary: str
    tags: tuple[str, ...]
    reasons: tuple[str, ...]


def _site_upper(site: object) -> tuple[str, str]:
    return str(getattr(site, "wt", "")).upper(), str(getattr(site, "mut", "")).upper()


def classify_mutation_route(
    sites: tuple,
    *,
    charge_conserving: bool,
    loop_residues: frozenset[int] | None = None,
) -> MutationRoute:
    """Classify a mutation group into the policy-table route.

    Priority (first match wins for `primary`):
      1. charge-changing           -> dssb
      2. large cavity deletion     -> rid_hydration   (Y/W/F/M -> A/G)
      3. Pro/Gly or insertion      -> endpoint_ensemble_rest
      4. titratable neighborhood   -> protonation_multi (H/D/E/K involved)
      5. otherwise                 -> baseline_12lambda
    Tags accumulate independently of the primary route.
    """
    tags: list[str] = []
    reasons: list[str] = []

    wt_mut = [(_site_upper(site)) for site in sites]

    involves_pro_gly = any("P" in pair or "G" in pair for pair in wt_mut)
    net_growth = sum(
        _SIDECHAIN_HEAVY_ATOMS.get(mut, 0) - _SIDECHAIN_HEAVY_ATOMS.get(wt, 0)
        for wt, mut in wt_mut
    )
    is_insertion = net_growth > 0
    large_cavity_deletion = any(
        wt in _LARGE_CAVITY_RESIDUES and mut in _SMALL_DELETION_TARGETS
        for wt, mut in wt_mut
    )
    titratable = any(wt in _TITRATABLE_RESIDUES or mut in _TITRATABLE_RESIDUES for wt, mut in wt_mut)
    in_loop = False
    if loop_residues:
        in_loop = any(int(getattr(site, "resseq", -1)) in loop_residues for site in sites)

    if titratable:
        tags.append("titratable_neighbor")
        reasons.append("involves H/D/E/K: PROPKA multi-state + pKa weighting advised")
    if in_loop:
        tags.append("flexible_loop")
        reasons.append("site in annotated loop/CDR region")
    if large_cavity_deletion:
        tags.append("basin_diversity")
        reasons.append("large-deletion class: prefer more basins over more lambda windows")

    if not charge_conserving:
        reasons.insert(0, "net charge change: DSSB double-system single-box")
        return MutationRoute(ROUTE_DSSB, tuple(tags), tuple(reasons))
    if large_cavity_deletion:
        reasons.insert(0, "Y/W/F/M -> A/G large cavity deletion: RID + water occupancy analysis")
        return MutationRoute(ROUTE_RID_HYDRATION, tuple(tags), tuple(reasons))
    if involves_pro_gly or is_insertion or in_loop:
        reasons.insert(0, "Pro/Gly / insertion / loop: endpoint ensemble + local REST2/RID")
        return MutationRoute(ROUTE_ENDPOINT_ENSEMBLE, tuple(tags), tuple(reasons))
    if titratable:
        reasons.insert(0, "titratable neighborhood without size/charge class: protonation path")
        return MutationRoute(ROUTE_PROTONATION, tuple(tags), tuple(reasons))
    return MutationRoute(ROUTE_BASELINE, tuple(tags), ("ordinary neutral mutation: 12-lambda baseline",))
