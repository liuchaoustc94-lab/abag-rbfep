#!/usr/bin/env python
"""Alchembed-style lambda sub-step pre-growth for a crashing FEP window (P3-a pilot).

Rationale: insertion/deletion windows crash because the alchemical sidechain is
forced into (or out of) a tightly packed pocket in one lambda jump. Alchembed
(Jakubec/Vymetal) instead grows the sidechain gradually: a series of short
soft-core EM runs at sub-lambda values interpolating the previous window's
state into the target window's state, each starting from the previous EM
output. The relaxed structure is then a physically adapted starting point for
the window's standard pre_md/production chain.

Usage:
  python tools/pregrow_window.py \
      --window-dir <job>/legs/complex/rep01/lambda_005 \
      --start-gro <job>/legs/complex/rep01/equilibration/npt.gro \
      --repeat-top <job>/legs/complex/rep01/system.top \
      --gmx /path/to/gmx --gmxlib <job>/artifacts/gmxlib \
      --substeps 4 --em-steps 1000
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys


def _parse_lambda_vectors(mdp_text: str) -> dict[str, list[float]]:
    vectors: dict[str, list[float]] = {}
    for key in ("coul-lambdas", "vdw-lambdas", "bonded-lambdas", "mass-lambdas", "fep-lambdas"):
        m = re.search(rf"^{re.escape(key)}\s*=\s*(.+)$", mdp_text, flags=re.MULTILINE)
        if m:
            vectors[key] = [float(x) for x in m.group(1).split()]
    return vectors


def _window_state_index(mdp_text: str) -> int:
    m = re.search(r"^init-lambda-state\s*=\s*(\d+)", mdp_text, flags=re.MULTILINE)
    if not m:
        raise ValueError("init-lambda-state not found in mdp")
    return int(m.group(1))


def _render_substep_mdp(base_mdp: str, coul: float, vdw: float, *, em_steps: int) -> str:
    """Rewrite the window pre_relax mdp into a single-state sub-step EM mdp."""
    text = base_mdp
    text = re.sub(r"^init-lambda-state\s*=\s*\d+", "init-lambda-state       = 0", text, flags=re.MULTILINE)
    if "coul-lambdas" in text:
        text = re.sub(r"^coul-lambdas\s*=.*$", f"coul-lambdas            = {coul:.5f}", text, flags=re.MULTILINE)
        text = re.sub(r"^vdw-lambdas\s*=.*$", f"vdw-lambdas             = {vdw:.5f}", text, flags=re.MULTILINE)
        text = re.sub(r"^bonded-lambdas\s*=.*$", f"bonded-lambdas          = {vdw:.5f}", text, flags=re.MULTILINE)
        text = re.sub(r"^mass-lambdas\s*=.*$", f"mass-lambdas            = {vdw:.5f}", text, flags=re.MULTILINE)
    else:
        # coupled schedule: interpolate the single fep vector end value directly
        text = re.sub(r"^fep-lambdas\s*=.*$", f"fep-lambdas             = {vdw:.5f}", text, flags=re.MULTILINE)
    text = re.sub(r"^nsteps\s*=\s*\d+", f"nsteps                  = {em_steps}", text, flags=re.MULTILINE)
    return text


def _run(cmd: list[str], *, env: dict[str, str], cwd: Path, log) -> None:
    proc = subprocess.run(cmd, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window-dir", required=True, type=Path)
    ap.add_argument("--start-gro", required=True, type=Path, help="structure to start the sub-step chain from")
    ap.add_argument("--repeat-top", required=True, type=Path)
    ap.add_argument("--gmx", required=True)
    ap.add_argument("--gmxlib", default=None)
    ap.add_argument("--substeps", type=int, default=4)
    ap.add_argument("--em-steps", type=int, default=1000)
    ap.add_argument("--mdrun-args", default="-ntmpi 1 -ntomp 4")
    args = ap.parse_args()

    window_dir = args.window_dir
    # Prefer the couple-intramol=yes variant as the EM base: the exclusionchecker
    # fatal ("perturbed excluded pairs beyond pair-list cutoff") is a hard crash
    # for spread hybrids under couple-intramol=no, and pre-growth EM is a
    # structure-relaxation phase only (production mdps are untouched).
    retry2_mdp = window_dir / "pre_relax.retry2.mdp"
    base_mdp_path = retry2_mdp if retry2_mdp.is_file() else window_dir / "pre_relax.mdp"
    base_mdp = base_mdp_path.read_text(encoding="utf-8")
    if "couple-intramol" in base_mdp:
        base_mdp = re.sub(r"^couple-intramol\s*=\s*no", "couple-intramol         = yes", base_mdp, flags=re.MULTILINE)
    else:
        base_mdp += "\ncouple-intramol         = yes\n"
    vectors = _parse_lambda_vectors(base_mdp)
    k = _window_state_index(base_mdp)
    if k == 0:
        raise SystemExit("window 0 needs no pre-growth (it is the physical start state)")

    if "coul-lambdas" in vectors:
        prev_coul, tgt_coul = vectors["coul-lambdas"][k - 1], vectors["coul-lambdas"][k]
        prev_vdw, tgt_vdw = vectors["vdw-lambdas"][k - 1], vectors["vdw-lambdas"][k]
    else:
        fep = vectors["fep-lambdas"]
        prev_coul = prev_vdw = fep[k - 1]
        tgt_coul = tgt_vdw = fep[k]

    print(f"[pregrow] window {window_dir.name} (state {k}): "
          f"coul {prev_coul:.3f}->{tgt_coul:.3f}, vdw {prev_vdw:.3f}->{tgt_vdw:.3f}, "
          f"{args.substeps} substeps x {args.em_steps} EM steps")

    env = dict(os.environ)
    if args.gmxlib:
        env["GMXLIB"] = str(args.gmxlib)
    env.setdefault("GMX_MAXBACKUP", "-1")

    out_dir = window_dir / "pregrow"
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "pregrow.log"
    current_gro = args.start_gro
    with log_path.open("w", encoding="utf-8") as log:
        for i in range(1, args.substeps + 1):
            frac = i / args.substeps
            coul = prev_coul + frac * (tgt_coul - prev_coul)
            vdw = prev_vdw + frac * (tgt_vdw - prev_vdw)
            tag = f"sub{i:02d}"
            mdp = out_dir / f"{tag}.mdp"
            tpr = out_dir / f"{tag}.tpr"
            gro = out_dir / f"{tag}.gro"
            mdp.write_text(_render_substep_mdp(base_mdp, coul, vdw, em_steps=args.em_steps), encoding="utf-8")
            print(f"[pregrow] {tag}: coul={coul:.3f} vdw={vdw:.3f}")
            _run([args.gmx, "grompp", "-f", str(mdp), "-c", str(current_gro), "-p", str(args.repeat_top),
                  "-o", str(tpr), "-maxwarn", "2"], env=env, cwd=window_dir, log=log)
            _run([args.gmx, "mdrun", "-s", str(tpr), "-deffnm", str(out_dir / tag),
                  *args.mdrun_args.split()], env=env, cwd=window_dir, log=log)
            if not gro.is_file() or gro.stat().st_size == 0:
                raise RuntimeError(f"sub-step {tag} produced no gro")
            current_gro = gro

    final_gro = window_dir / "pregrown.gro"
    final_gro.write_bytes(current_gro.read_bytes())
    print(f"[pregrow] OK -> {final_gro}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
