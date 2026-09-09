#!/usr/bin/env bash
# 1BJ1 endpoint-ensemble fix (2026-08-19): rerun the 4 worst 1BJ1 jobs with
# 40 ns endpoint ensemble + frame injection (E4 pattern).
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
ROOTDIR="${ROOT}/runs/benchmarks/abbind_validation_panels_20260812/abbind_1bj1_core_v1"
export PATH="${ROOT}/.venv/bin:${PATH}"
JOBS=(1bj1-antigen-w-g92a 1bj1-antigen-w-g88a 1bj1-antigen-w-q89a 1bj1-antigen-w-i91a)

for job in "${JOBS[@]}"; do
  for leg in complex apo; do
    rep="${ROOTDIR}/jobs/${job}/legs/${leg}/rep01"
    eq="${rep}/equilibration"
    [ -f "${eq}/npt.gro" ] || { echo "[1bj1] skip ${job}/${leg}: no npt.gro"; continue; }
    [ -f "${eq}/freemd40.xtc" ] && continue
    cat > "${eq}/freemd.mdp" <<'MDP'
integrator = sd
dt = 0.002
nsteps = 20000000
ld-seed = 52011
cutoff-scheme = Verlet
verlet-buffer-tolerance = -1
rlist = 1.25
rcoulomb = 1.25
rvdw = 1.25
vdw-type = Cut-off
vdw-modifier = Potential-switch
rvdw-switch = 1.0
coulombtype = PME
pme-order = 4
fourierspacing = 0.12
DispCorr = EnerPres
tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = 310.0
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 2.0
compressibility = 4.5e-5
ref-p = 1.0
constraints = h-bonds
constraint-algorithm = lincs
pbc = xyz
nstxout-compressed = 12500
compressed-x-grps = System
MDP
    "${GMX}" grompp -f "${eq}/freemd.mdp" -c "${eq}/npt.gro" -p "${rep}/system.top" \
      -o "${eq}/freemd40.tpr" -maxwarn 2 > "${eq}/freemd40.log" 2>&1 && \
    ABAG_RBFE_VISIBLE_GPUS="${BJ1_GPU:-0}" "${GMX}" mdrun -s "${eq}/freemd40.tpr" \
      -deffnm "${eq}/freemd40" -ntmpi 1 -ntomp 4 >> "${eq}/freemd40.log" 2>&1 &
    echo "[1bj1] launched 40ns ${job}/${leg}"
  done
done
wait
echo "[1bj1] free-MD done; injecting frames"
for job in "${JOBS[@]}"; do
  for leg in complex apo; do
    eq="${ROOTDIR}/jobs/${job}/legs/${leg}/rep01/equilibration"
    [ -f "${eq}/freemd40.xtc" ] || { echo "[1bj1] no traj ${job}/${leg}"; continue; }
    i=1
    for frame_ps in 12500 25000 37500; do
      target="${ROOTDIR}/jobs/${job}/legs/${leg}/rep0${i}/equilibration"
      printf 'System\nSystem\n' | "${GMX}" trjconv -f "${eq}/freemd40.xtc" -s "${eq}/freemd40.tpr" \
        -o "${target}/npt_seed.gro" -dump "${frame_ps}" -pbc mol -ur compact > /dev/null 2>&1
      if [ -s "${target}/npt_seed.gro" ]; then
        cp "${target}/npt.gro" "${target}/npt.gro.orig" 2>/dev/null || true
        cp "${target}/npt_seed.gro" "${target}/npt.gro"
      fi
      i=$((i+1))
    done
  done
done
"${ROOT}/.venv/bin/python" - <<'PYEOF'
import glob, os, shutil
rootdir = "/mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_validation_panels_20260812/abbind_1bj1_core_v1"
for job in ["1bj1-antigen-w-g92a", "1bj1-antigen-w-g88a", "1bj1-antigen-w-q89a", "1bj1-antigen-w-i91a"]:
    jd = f"{rootdir}/jobs/{job}"
    for pat in ["legs/*/rep*/lambda_*", "legs/*/rep*/bar"]:
        for p in glob.glob(f"{jd}/{pat}"):
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
    for s in ["sample", "bar", "qc", "report"]:
        f = f"{jd}/stages/{s}.json"
        if os.path.isfile(f):
            os.remove(f)
    for f in glob.glob(f"{jd}/results/*.json"):
        os.remove(f)
    print("wiped", job)
PYEOF
for job in "${JOBS[@]}"; do
  echo "[1bj1] FEP ${job}"
  "${ABAG_RBFE}" resume "${job}" --batch-dir "${ROOTDIR}" --execute 2>&1 | tail -1
done
echo "[1bj1] $(date --iso-8601=seconds) all done"
