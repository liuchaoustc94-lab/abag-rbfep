#!/usr/bin/env python3
"""Build the single-point-focused companion panel.

The broad panel remains under benchmarks/antibody_new10.  This companion
panel intentionally excludes all combination-mutant rows from its primary
targets and uses clearly labelled Graphinity Flex-ddG controls to reach ten
structure-level systems.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import build_new_antibody_validation_panel as base


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmarks" / "antibody_new10_single_point"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skempi", type=Path, default=base.SKEMPI_DEFAULT)
    parser.add_argument("--graphinity", type=Path, default=base.GRAPHINITY_DEFAULT)
    parser.add_argument(
        "--graphinity-sequences",
        type=Path,
        default=base.GRAPHINITY_SEQUENCES_DEFAULT,
    )
    args = parser.parse_args()

    base.OUT = OUT
    base.TARGETS = base.SINGLE_POINT_TARGETS
    base.build(args.skempi, args.graphinity, args.graphinity_sequences)
    print(f"Wrote single-point-focused panel to {OUT}")


if __name__ == "__main__":
    main()
