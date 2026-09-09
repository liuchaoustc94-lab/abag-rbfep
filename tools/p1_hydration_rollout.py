#!/usr/bin/env python
"""P1 hydration rollout driver (2026-09-08): per-window hydrated seeding for
all remaining aromatic-deletion (Y/W/F -> A) benchmark points.

Per target job:
  1. setup_metad_hydration.py  (build 25ns MetaD system at mutant endpoint)
  2. run MetaD (GPU)
  3. analyze colvar.dat -> pick 2 hydrated frames (sustained cv>=8 epochs)
  4. extract frames as rep02/rep03 hydrated seeds
  5. write config/window_seeds.json (last 3 windows of complex rep02/rep03)
  6. wipe those windows' outputs + sample/bar/qc/report stages + results, resume
  7. record before/after ddG in the rollout state file

State file makes the campaign resumable. Two workers on GPU1/GPU2 by default.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import time

ROOT = Path("/mnt/data/liuchao/abag-rbfep")
PY = ROOT / ".venv/bin/python"
GMX = "/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
PLUMED_KERNEL = "/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/plumed-mpi/lib/libplumedKernel.so"
STATE_DIR = ROOT / "runs/analysis/p1_rollout_20260908"
LOG = STATE_DIR / "rollout.log"
STATE = STATE_DIR / "state.json"
PAIRS = ROOT / "runs/analysis/protonation_reweight_20260831/pairs_protonation.csv"
DONE = {("3be1", "w50a"), ("1dqj", "y50a"), ("1vfb", "w52a")}
REF_ATOM = {"W": "CZ3", "Y": "CZ", "F": "CZ"}

_lock = threading.Lock()


def log(msg: str) -> None:
    with _lock:
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        with open(LOG, "a") as f:
            f.write(line + "\n")


def load_state() -> dict:
    if STATE.is_file():
        return json.load(open(STATE))
    return {}


def save_state(state: dict) -> None:
    with _lock:
        STATE.write_text(json.dumps(state, indent=1))


def run(cmd, *, log_prefix, env=None, check=True, cwd=None):
    e = dict(os.environ)
    e.update(env or {})
    with open(STATE_DIR / f"{log_prefix}.log", "a") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=e, cwd=cwd)
    if check and proc.returncode != 0:
        raise RuntimeError(f"{log_prefix} failed rc={proc.returncode}")
    return proc.returncode


def pick_frames(colvar: Path, n: int = 2):
    ts, cv = [], []
    for line in open(colvar):
        if line.startswith("#"):
            continue
        p = line.split()
        ts.append(float(p[0])); cv.append(float(p[1]))
    runs, cur = [], []
    for t, v in zip(ts, cv):
        if t > 5000 and v >= 8.0:
            cur.append(t)
        else:
            if len(cur) >= 300:
                runs.append(cur)
            cur = []
    if len(cur) >= 300:
        runs.append(cur)
    picks = [int(r[len(r) // 2]) for r in runs[:n]]
    return picks


def process(job_dir: Path, complex_id: str, mut: str, wt: str, gpu: str, state: dict) -> None:
    tag = f"{complex_id}-{mut}"
    state.setdefault(tag, {})
    st = state[tag]
    if st.get("stage") == "done":
        log(f"{tag}: already done, skip")
        return
    try:
        rep1 = job_dir / "legs/complex/rep01"
        colvar = rep1 / "colvar.dat"
        if not (rep1 / "hydra.tpr").is_file():
            st["stage"] = "setup"
            log(f"{tag}: setup MetaD")
            run([str(PY), str(ROOT / "tools/setup_metad_hydration.py"), str(job_dir),
                 "--resname", f"{wt}2A", "--ref-atom", REF_ATOM[wt]], log_prefix=tag)
            save_state(state)
        if not st.get("metad_done"):
            st["stage"] = "metad"
            log(f"{tag}: MetaD 25ns on gpu{gpu}")
            env = {"GMXLIB": str(job_dir / "artifacts/gmxlib"), "GMX_MAXBACKUP": "-1",
                   "PLUMED_KERNEL": PLUMED_KERNEL, "CUDA_VISIBLE_DEVICES": gpu}
            # 断点续跑（2026-09-09）：有 cpt 就接着跑，PLUMED 用 RESTART 读回 HILLS
            cpt = rep1 / "hydra.cpt"
            hills = rep1 / "HILLS"
            if cpt.is_file() and hills.is_file():
                restart_dat = rep1 / "hydra_plumed_restart.dat"
                txt = (rep1 / "hydra_plumed.dat").read_text()
                txt = txt.replace("metad: METAD ", "metad: METAD RESTART=YES ")
                txt = "RESTART\n" + txt
                restart_dat.write_text(txt)
                log(f"{tag}: 从 cpt 续跑（已有 {sum(1 for l in open(rep1/'colvar.dat') if not l.startswith('#'))} ps）")
                run([GMX, "mdrun", "-s", "hydra.tpr", "-cpi", "hydra.cpt", "-deffnm", "hydra",
                     "-plumed", "hydra_plumed_restart.dat", "-ntmpi", "1", "-ntomp", "8"],
                    log_prefix=tag, env=env, cwd=rep1)
            else:
                run([GMX, "mdrun", "-s", "hydra.tpr", "-deffnm", "hydra",
                     "-plumed", "hydra_plumed.dat", "-ntmpi", "1", "-ntomp", "8"],
                    log_prefix=tag, env=env, cwd=rep1)
            st["metad_done"] = True
            save_state(state)
        # frames
        picks = pick_frames(colvar)
        if len(picks) < 2:
            st["stage"] = "no_hydrated_frames"
            st["picks"] = picks
            save_state(state)
            log(f"{tag}: 持续水合段不足（{len(picks)} 帧）——记录为欠水合失败")
            return
        st["stage"] = "seeding"
        seeds = {}
        for rep, pick in (("rep02", picks[0]), ("rep03", picks[1])):
            tgt = job_dir / "legs/complex" / rep / "equilibration"
            if not tgt.is_dir():
                continue
            seed = tgt / "npt_seed_hydrated.gro"
            proc = subprocess.run(
                f"printf 'System\\nSystem\\n' | {GMX} trjconv -f {rep1/'hydra.xtc'} -s {rep1/'hydra.tpr'} -o {seed} -dump {pick} -pbc mol -ur compact",
                shell=True, cwd=job_dir, env={**os.environ, "CUDA_VISIBLE_DEVICES": ""},
                capture_output=True)
            if not seed.is_file() or seed.stat().st_size == 0:
                raise RuntimeError(f"frame extraction failed {rep}@{pick}")
            seeds[rep] = seed
        # window seeds: last 3 windows of complex rep02/rep03; wipe their outputs
        # so the sample stage reruns exactly those from the hydrated seed.
        seed_map = {}
        for rep, seed in seeds.items():
            rep_dir = job_dir / "legs/complex" / rep
            for w in sorted(rep_dir.glob("lambda_[0-9][0-9][0-9]"))[-3:]:
                seed_map[f"complex/{rep}/{w.name}"] = str(seed)
                for fn in ("dhdl.xvg", "md.gro", "md.log", "topol.tpr"):
                    f = w / fn
                    if f.exists():
                        f.unlink()
        (job_dir / "config/window_seeds.json").write_text(json.dumps(seed_map, indent=1))
        old = None
        ddg_path = job_dir / "results/ddg_summary.json"
        if ddg_path.is_file():
            old = json.load(open(ddg_path)).get("ddg_kcal_mol")
        for s in ("sample", "bar", "qc", "report"):
            f = job_dir / "stages" / f"{s}.json"
            if f.is_file():
                f.unlink()
        for f in (job_dir / "results").glob("*.json"):
            f.unlink()
        for bd in (job_dir / "legs").glob("*/rep*/bar"):
            subprocess.run(["rm", "-rf", str(bd)])
        log(f"{tag}: reseeded resume (old ddG={old})")
        env = {"PATH": f"{ROOT}/.venv/bin:" + os.environ["PATH"],
               "ABAG_RBFE_VISIBLE_GPUS": gpu}
        run([str(ROOT / ".venv/bin/abag-rbfe"), "resume", job_dir.name,
             "--batch-dir", str(job_dir.parent.parent), "--execute"],
            log_prefix=tag, env=env, check=False)
        new = None
        if ddg_path.is_file():
            new = json.load(open(ddg_path)).get("ddg_kcal_mol")
        st["stage"] = "done"
        st["old_ddg"] = old
        st["new_ddg"] = new
        st["picks"] = picks
        save_state(state)
        log(f"{tag}: DONE old={old} new={new}")
    except Exception as exc:  # noqa: BLE001
        st["stage"] = "error"
        st["error"] = str(exc)
        save_state(state)
        log(f"{tag}: ERROR {exc}")


def worker(name: str, gpu: str, queue: list, state: dict) -> None:
    for job_dir, cx, mut, wt in queue:
        process(job_dir, cx, mut, wt, gpu, state)


def main() -> int:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    targets = []
    for r in csv.DictReader(open(PAIRS)):
        mut = r["mut"]
        m = re.match(r"^([ywf])\d+a$", mut)
        if not m:
            continue
        key = (r["complex"], mut)
        if key in DONE or r["complex"] == "3hfm":  # 3HFM = DSSB track, not this pipeline
            continue
        err = abs(float(r["pred"]) - float(r["exp"]))
        targets.append((err, Path(r["job_dir"]), r["complex"], mut, m.group(1).upper()))
    targets.sort(key=lambda t: -t[0])
    targets = [(jd, cx, mut, wt) for _, jd, cx, mut, wt in targets]
    print(f"rollout targets: {len(targets)}")
    for t in targets:
        print("  ", t[1], t[2])
    state = load_state()
    half = (len(targets) + 1) // 2
    shards = [targets[:half], targets[half:]]
    threads = [threading.Thread(target=worker, args=(f"w{i}", str(i + 1), shards[i], state))
               for i in range(2) if shards[i]]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    log("rollout finished")
    return 0


if __name__ == "__main__":
    sys.exit(main())
