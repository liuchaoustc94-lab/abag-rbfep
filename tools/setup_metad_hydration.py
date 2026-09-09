#!/usr/bin/env python
"""Set up a cavity-hydration MetaD run at the mutant endpoint (lambda=last)
for a completed FEP job (P1-b rollout, generalized from the w50a pilot).

Creates in <job>/legs/complex/rep01/:
  hydra_groups.ndx  (waterO + cavity groups; simple names, PLUMED-safe)
  system_hydra.top  (system.top + posre includes inside moleculetype scope)
  hydra.mdp         (last-window production settings, 25 ns, POSRES, fresh vel)
  hydra_plumed.dat  (WT-MetaD on COORDINATION(cavity COM, water O) with NLIST)
  hydra.tpr         (grompp -c/-r equilibration/npt.gro)

Usage:
  tools/setup_metad_hydration.py <job_dir> --resname W2A --ref-atom CZ3
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys

GMX = "/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"

PLUMED_TEMPLATE = """# cavity hydration MetaD (mutant endpoint). NLIST = 3x speedup vs naive CV.
UNITS LENGTH=nm ENERGY=kj/mol TIME=ps
wat: GROUP NDX_FILE=hydra_groups.ndx NDX_GROUP=waterO
cavsel: GROUP NDX_FILE=hydra_groups.ndx NDX_GROUP=cavity
cav: COM ATOMS=cavsel
cv: COORDINATION GROUPA=cav GROUPB=wat R_0=0.45 NN=8 MM=16 NLIST NL_CUTOFF=0.9 NL_STRIDE=10
metad: METAD ARG=cv PACE=500 HEIGHT=2.0 SIGMA=1.5 TEMP=310 BIASFACTOR=15 FILE=HILLS
PRINT ARG=cv,metad.bias STRIDE=500 FILE=colvar.dat
"""


def run(cmd: list[str], **kw) -> None:
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        raise RuntimeError(
            "failed: " + " ".join(cmd)
            + "\nSTDOUT:\n" + proc.stdout[-1200:]
            + "\nSTDERR:\n" + proc.stderr[-1200:]
        )


def select_indices(tpr: Path, selection: str) -> list[int]:
    out = Path("/tmp/hydra_sel.ndx")
    run([GMX, "select", "-s", str(tpr), "-select", selection, "-on", str(out)])
    idx: list[int] = []
    for line in open(out):
        if line.startswith("["):
            continue
        idx.extend(int(t) for t in line.split())
    return idx


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("job_dir", type=Path)
    ap.add_argument("--resname", required=True)
    ap.add_argument("--ref-atom", required=True)
    ap.add_argument("--ns", type=float, default=25.0)
    args = ap.parse_args()
    args.job_dir = args.job_dir.resolve()

    rep = args.job_dir / "legs" / "complex" / "rep01"
    lambda_dirs = sorted(rep.glob("lambda_[0-9][0-9][0-9]"))
    if not lambda_dirs:
        raise SystemExit(f"no lambda dirs under {rep}")
    last = lambda_dirs[-1]
    tpr0 = lambda_dirs[0] / "topol.tpr"
    if not tpr0.is_file():
        tpr0 = last / "topol.tpr"

    water_idx = select_indices(tpr0, "name OW")
    cavity_idx = select_indices(
        tpr0, f"not resname SOL NA CL and within 0.6 of (resname {args.resname} and name {args.ref_atom})"
    )
    print(f"water O: {len(water_idx)}, cavity atoms: {len(cavity_idx)}")
    with open(rep / "hydra_groups.ndx", "w") as f:
        f.write("[ waterO ]\n")
        for i in range(0, len(water_idx), 15):
            f.write(" ".join(str(x) for x in water_idx[i : i + 15]) + "\n")
        f.write("[ cavity ]\n")
        for i in range(0, len(cavity_idx), 15):
            f.write(" ".join(str(x) for x in cavity_idx[i : i + 15]) + "\n")

    # system_hydra.top: posre includes inside moleculetype scope
    src = (rep / "system.top").read_text()
    posres = {p.name for p in rep.glob("posre_Protein_chain_*.itp")}
    out = []
    for line in src.splitlines():
        out.append(line)
        m = re.search(r'#include "(?:pmx_)?topol_Protein_chain_(\w+)\.itp"', line)
        if m:
            pf = f"posre_Protein_chain_{m.group(1)}.itp"
            if pf in posres:
                out += ["", "#ifdef POSRES", f'#include "{pf}"', "#endif", ""]
    (rep / "system_hydra.top").write_text("\n".join(out) + "\n")

    # hydra.mdp from the last window's production.mdp
    mdp = (last / "production.mdp").read_text()
    lines = []
    for line in mdp.splitlines():
        key = line.split("=")[0].strip() if "=" in line else ""
        if key == "nsteps":
            lines.append(f"nsteps                  = {int(args.ns * 1000 / 0.002)}")
        elif key == "ld-seed":
            lines.append("ld-seed                 = 26090101")
        elif key == "nstxout-compressed":
            lines.append("nstxout-compressed      = 12500")
        elif key == "continuation":
            lines.append("continuation            = no")
        elif key == "gen-vel":
            lines.append("gen-vel                 = yes\ngen-temp                = 310.0")
        else:
            lines.append(line)
    (rep / "hydra.mdp").write_text("define                  = -DPOSRES\n" + "\n".join(lines) + "\n")

    (rep / "hydra_plumed.dat").write_text(PLUMED_TEMPLATE)

    npt = rep / "equilibration" / "npt.gro"
    env = dict(os.environ)
    gmxlib = args.job_dir / "artifacts" / "gmxlib"
    env["CUDA_VISIBLE_DEVICES"] = ""
    if gmxlib.is_dir():
        env["GMXLIB"] = str(gmxlib)
    env["CUDA_VISIBLE_DEVICES"] = ""  # grompp: force CPU (CUDA #700 collateral fix)
    run([GMX, "grompp", "-f", str(rep / "hydra.mdp"), "-c", str(npt), "-r", str(npt),
         "-p", str(rep / "system_hydra.top"), "-o", str(rep / "hydra.tpr"), "-maxwarn", "3"],
        cwd=rep, env=env)
    print(f"OK -> {rep / 'hydra.tpr'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
