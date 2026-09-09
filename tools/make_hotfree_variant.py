#!/usr/bin/env python
"""Build the hot-region-unrestrained variant of a hydration MetaD run (A/B arm B).

Takes a job where tools/setup_metad_hydration.py has already run (hydra.mdp,
hydra_groups.ndx, system_hydra.top, hydra.tpr exist) and produces:
  posre_Protein_chain_*.hotfree.itp  (posre minus atoms within --radius of the
                                      hybrid residue)
  system_hotfree.top                 (includes the hotfree posre variants)
  hydra_hotfree.tpr                  (grompp with the same hydra.mdp)

Physics: everything frozen except the cavity neighborhood, so the MetaD water
CV is coupled to local sidechain/backbone gating motions (arm A = fully frozen
isolates hydration alone; arm B tests the coupled water+conformation ensemble).

Usage:
  tools/make_hotfree_variant.py <job_dir> --resname W2A --ref-atom CZ3 [--radius 0.8]
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys

GMX = "/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"


def run(cmd: list[str], **kw) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        raise RuntimeError("failed: " + " ".join(cmd) + "\n" + proc.stdout[-800:] + "\n" + proc.stderr[-800:])


def select_indices(tpr: Path, selection: str) -> set[int]:
    out = Path("/tmp/hotfree_sel.ndx")
    run([GMX, "select", "-s", str(tpr), "-select", selection, "-on", str(out)])
    idx: set[int] = set()
    for line in open(out):
        if not line.startswith("["):
            idx.update(int(t) for t in line.split())
    return idx


def itp_atom_count(itp: Path) -> int:
    in_atoms = False
    n = 0
    for line in open(itp):
        s = line.strip()
        if s.startswith("["):
            if in_atoms:
                break  # next section
            in_atoms = s.lower().startswith("[ atoms")
            continue
        if not in_atoms or not s or s.startswith(";"):
            continue
        parts = s.split()
        if parts and parts[0].isdigit():
            n = max(n, int(parts[0]))
    return n


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("job_dir", type=Path)
    ap.add_argument("--resname", required=True)
    ap.add_argument("--ref-atom", required=True)
    ap.add_argument("--radius", type=float, default=0.8)
    args = ap.parse_args()
    job_dir = args.job_dir.resolve()
    rep = job_dir / "legs" / "complex" / "rep01"

    # hot atoms: global indices within radius of the hybrid residue
    hot = select_indices(
        rep / "hydra.tpr",
        f"not resname SOL NA CL and within {args.radius} of resname {args.resname}",
    )
    print(f"hot atoms (global): {len(hot)}")

    # chain boundaries from system.top molecule order + per-chain atom counts
    top = (rep / "system.top").read_text()
    mol_section = top.split("[ molecules ]")[-1]
    chains: list[tuple[str, int]] = []  # (itp stem, natoms)
    for line in mol_section.splitlines():
        s = line.strip()
        if not s or s.startswith(";"):
            continue
        mol = s.split()[0]
        if not mol.startswith("Protein_chain_"):
            continue
        for cand in (rep / f"pmx_topol_{mol}.itp", rep / f"topol_{mol}.itp"):
            if cand.is_file():
                chains.append((mol, itp_atom_count(cand)))
                break
    # global offsets
    offset = 0
    chain_of: dict[int, tuple[str, int]] = {}  # global idx -> (mol, local idx)
    for mol, natoms in chains:
        for local in range(1, natoms + 1):
            chain_of[offset + local] = (mol, local)
        offset += natoms

    hot_by_chain: dict[str, set[int]] = {}
    for g in hot:
        if g in chain_of:
            mol, local = chain_of[g]
            hot_by_chain.setdefault(mol, set()).add(local)

    # write hotfree posre variants
    for mol, _ in chains:
        src = rep / f"posre_{mol}.itp"
        if not src.is_file():
            continue
        drop = hot_by_chain.get(mol, set())
        kept = 0
        out_lines = []
        for line in open(src):
            s = line.strip()
            if s and not s.startswith((";", "[")) :
                parts = s.split()
                if parts and parts[0].isdigit() and int(parts[0]) in drop:
                    continue
                if parts and parts[0].isdigit():
                    kept += 1
            out_lines.append(line)
        (rep / f"posre_{mol}.hotfree.itp").write_text("".join(out_lines))
        print(f"posre_{mol}: kept {kept}, dropped {len(drop)} (hot)")

    # system_hotfree.top: swap posre includes to hotfree variants
    src = (rep / "system_hydra.top").read_text()
    src = re.sub(r'#include "(posre_Protein_chain_\w+)\.itp"', r'#include "\1.hotfree.itp"', src)
    (rep / "system_hotfree.top").write_text(src)

    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = ""
    gmxlib = job_dir / "artifacts" / "gmxlib"
    if gmxlib.is_dir():
        env["GMXLIB"] = str(gmxlib)
    npt = rep / "equilibration" / "npt.gro"
    run([GMX, "grompp", "-f", str(rep / "hydra.mdp"), "-c", str(npt), "-r", str(npt),
         "-p", str(rep / "system_hotfree.top"), "-o", str(rep / "hydra_hotfree.tpr"),
         "-maxwarn", "3"], cwd=rep, env=env)
    print(f"OK -> {rep / 'hydra_hotfree.tpr'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
