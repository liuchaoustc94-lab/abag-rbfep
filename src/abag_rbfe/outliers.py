"""Sampson-2024-style outlier classification and empirical correction.

Classifies mutation ddG outliers using structure-informed features and applies
a single-parameter shrinkage correction to the dominant systematic class
(buried large aromatic deletions, e.g. Y/W/F->Ala), whose alchemical
transformations systematically overestimate |ddG| in fixed-charge force fields.

Reference: Sampson et al. 2024, JCTC (RMSE 1.35 -> 1.03 via outlier flagging
+ single-parameter charged-outlier correction).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

SIDECHAIN_HEAVY_ATOMS: dict[str, int] = {
    "G": 0, "A": 1, "S": 2, "C": 2, "T": 3, "V": 3, "P": 3,
    "N": 4, "D": 4, "L": 4, "I": 4, "M": 4, "K": 5, "Q": 5,
    "E": 5, "H": 6, "F": 7, "R": 7, "Y": 8, "W": 10,
}
AROMATIC = frozenset({"F", "Y", "W"})
CHARGED = frozenset({"D", "E", "K", "R", "H"})
NOMINAL_CHARGE: dict[str, int] = {"R": 1, "K": 1, "H": 0, "D": -1, "E": -1}
LARGE_DELETION_THRESHOLD = 3  # |sidechain heavy atoms removed| >= 3


def nominal_charge(residue: str) -> int:
    return NOMINAL_CHARGE.get(residue.upper(), 0)


def mutation_size_change(wt: str, mut: str) -> int:
    """Sidechain heavy-atom change (mut - wt); negative for deletions."""
    return SIDECHAIN_HEAVY_ATOMS.get(mut.upper(), 0) - SIDECHAIN_HEAVY_ATOMS.get(wt.upper(), 0)


def _ca_neighbors(structure_path: Path, chain_id: str, resseq: int, cutoff_angstrom: float) -> list[tuple[str, int, str, float]]:
    """C-alpha neighbors of a residue within cutoff (structure-informed burial proxy)."""
    target: tuple[float, float, float] | None = None
    others: list[tuple[str, int, str, float, float, float, float]] = []
    with open(structure_path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if line[21] == chain_id and int(line[22:26]) == resseq:
                target = (x, y, z)
            else:
                others.append((line[21], int(line[22:26]), line[17:20].strip(), x, y, z))
    if target is None:
        return []
    neighbors = []
    for chain, seq, resname, x, y, z in others:
        dist = math.dist(target, (x, y, z))
        if dist <= cutoff_angstrom:
            neighbors.append((chain, seq, resname, dist))
    return neighbors


def classify_outlier_features(
    *,
    wt: str,
    mut: str,
    structure_path: Path | None = None,
    chain_id: str = "",
    resseq: int = 0,
) -> dict[str, Any]:
    """Sampson-style outlier features for one mutation point."""
    wt = wt.upper()
    mut = mut.upper()
    size_change = mutation_size_change(wt, mut)
    features: dict[str, Any] = {
        "size_change": size_change,
        "is_aromatic_deletion": wt in AROMATIC and size_change < 0,
        "is_large_deletion": abs(size_change) >= LARGE_DELETION_THRESHOLD and size_change < 0,
        "wt_charged": wt in CHARGED,
        "mut_charged": mut in CHARGED,
        "is_charge_changing": nominal_charge(wt) != nominal_charge(mut),
        "buried_neighbor_count": None,
        "charged_neighbor_count": None,
    }
    if structure_path is not None and Path(structure_path).is_file() and chain_id and resseq:
        neighbors = _ca_neighbors(Path(structure_path), chain_id, int(resseq), 10.0)
        features["buried_neighbor_count"] = len(neighbors)
        features["charged_neighbor_count"] = sum(1 for *_rest, resname, _d in neighbors if resname in {"ASP", "GLU", "LYS", "ARG", "HIS"})
    features["is_buried_aromatic_large_deletion"] = bool(
        features["is_aromatic_deletion"] and features["is_large_deletion"]
    )
    return features


def outlier_class_label(features: dict[str, Any]) -> str:
    """Sampson-style class label for reporting; empty string = no outlier."""
    if features.get("is_buried_aromatic_large_deletion"):
        return "buried_aromatic_large_deletion"
    if features.get("is_charge_changing"):
        return "charge_changing"
    if features.get("is_large_deletion"):
        return "large_deletion"
    return ""


def fit_shrinkage_alpha(
    predicted: list[float],
    experimental: list[float],
    flags: list[bool],
) -> float:
    """Fit the single shrinkage parameter alpha for flagged points: choose alpha
    minimizing MAE(alpha * pred) over the flagged subset."""
    flagged = [(p, e) for p, e, f in zip(predicted, experimental, flags) if f]
    if not flagged:
        return 1.0
    best_alpha, best_mae = 1.0, float("inf")
    for alpha_i in range(10, 101):
        alpha = alpha_i / 100.0
        mae = sum(abs(alpha * p - e) for p, e in flagged) / len(flagged)
        if mae < best_mae:
            best_mae, best_alpha = mae, alpha
    return best_alpha


def apply_outlier_correction(
    predicted: list[float],
    flags: list[bool],
    alpha: float,
) -> list[float]:
    """Apply the single-parameter shrinkage to flagged points only."""
    return [alpha * p if f else p for p, f in zip(predicted, flags)]
