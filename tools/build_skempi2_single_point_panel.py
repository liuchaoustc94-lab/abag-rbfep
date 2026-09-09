#!/usr/bin/env python3
"""Materialize the SKEMPI2-only antibody-side single-point reserve panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import build_new_antibody_validation_panel as base


ROOT = Path(__file__).resolve().parents[1]
base.OUT = ROOT / "benchmarks" / "antibody_skempi2_single_point"
base.TARGETS = base.SKEMPI_SINGLE_POINT_TARGETS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skempi", type=Path, default=base.SKEMPI_DEFAULT)
    args = parser.parse_args()
    # Graphinity paths are accepted by the shared builder for API stability;
    # the all-experimental target set makes the builder skip those inputs.
    base.build(args.skempi, base.GRAPHINITY_DEFAULT, base.GRAPHINITY_SEQUENCES_DEFAULT)
    print(f"Wrote SKEMPI2-only panel to {base.OUT}")


if __name__ == "__main__":
    main()
