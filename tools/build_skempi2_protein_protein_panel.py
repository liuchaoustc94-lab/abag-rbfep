#!/usr/bin/env python3
"""Materialize an experimental SKEMPI2 protein-protein single-point panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import build_new_antibody_validation_panel as base


ROOT = Path(__file__).resolve().parents[1]
base.OUT = ROOT / "benchmarks" / "protein_protein_skempi2_single_point"
base.TARGETS = base.PROTEIN_PROTEIN_TARGETS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skempi", type=Path, default=base.SKEMPI_DEFAULT)
    args = parser.parse_args()
    base.build(args.skempi, base.GRAPHINITY_DEFAULT, base.GRAPHINITY_SEQUENCES_DEFAULT)
    print(f"Wrote SKEMPI2 protein-protein panel to {base.OUT}")


if __name__ == "__main__":
    main()
