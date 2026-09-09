#!/usr/bin/env bash
# PR-4b: charge-changing 3HFM mutations with E4-style deep endpoint ensemble.
# Rationale (Patel Table 2): charge flips need ns-scale electrostatic relaxation;
# our 20ps windows are in the unrelaxed regime (k97d +39 vs exp +6.77).
# Protocol: 40 ns free-MD per DSSB box from npt endpoint -> 3 frames
# (12.5/25/37.5 ns) -> rep01/02/03 seeds -> 24 lambda x 3 rep x 50 ps production.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
export PATH="${ROOT}/.venv/bin:${PATH}"

declare -A JOBS=(
  ["dssb_smoke_d32n"]="3hfm-patel-2021-antibody-h-d32n"
)
BATCH2="${ROOT}/runs/real_cases/3hfm_patel_dssb_pr4"
JOBS2=(3hfm-patel-2021-antigen-y-r21a 3hfm-patel-2021-antigen-y-d101k 3hfm-patel-2021-antibody-l-n31e 3hfm-patel-2021-antigen-y-k97d)

# 0) 提升 production_ps 到 50（协议级弛豫增强）
"${ROOT}/.venv/bin/python" - <<'PYEOF'
import json, re
jobs = [
    "runs/real_cases/dssb_smoke_d32n/jobs/3hfm-patel-2021-antibody-h-d32n",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-r21a",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-d101k",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antibody-l-n31e",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-k97d",
]
for jd in jobs:
    spec = json.load(open(f"{jd}/job_spec.json"))
    spec["protocol"]["production_ps"] = 50
    json.dump(spec, open(f"{jd}/job_spec.json", "w"), indent=2)
    yp = f"{jd}/config/protocol.yml"
    t = open(yp).read()
    t = re.sub(r"production_ps: [\d.]+", "production_ps: 50", t)
    open(yp, "w").write(t)
    print("patched production_ps=50:", jd.split("/")[-1])
PYEOF

# 1) 每个 job 的 dssb 腿 rep01 起 40ns 自由 MD
launch_freemd() {
  local jd="$1"
  local rep="${jd}/legs/dssb/rep01"
  local eq="${rep}/equilibration"
  [ -f "${eq}/npt.gro" ] || { echo "[pr4b] skip ${jd}: no npt.gro"; return 0; }
  [ -f "${eq}/freemd40.xtc" ] && { echo "[pr4b] freemd exists ${jd##*/}"; return 0; }
  cat > "${eq}/freemd.mdp" <<'MDP'
integrator = sd
dt = 0.002
nsteps = 20000000
ld-seed = 52001
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
  ABAG_RBFE_VISIBLE_GPUS="${PR4B_GPU:-0}" "${GMX}" mdrun -s "${eq}/freemd40.tpr" \
    -deffnm "${eq}/freemd40" -ntmpi 1 -ntomp 4 >> "${eq}/freemd40.log" 2>&1 &
  echo "[pr4b] launched 40ns free-MD ${jd##*/}"
}

GPU_I=0
for jd in "${ROOT}/runs/real_cases/dssb_smoke_d32n/jobs/3hfm-patel-2021-antibody-h-d32n" \
          "${BATCH2}/jobs/3hfm-patel-2021-antigen-y-r21a" \
          "${BATCH2}/jobs/3hfm-patel-2021-antigen-y-d101k" \
          "${BATCH2}/jobs/3hfm-patel-2021-antibody-l-n31e" \
          "${BATCH2}/jobs/3hfm-patel-2021-antigen-y-k97d"; do
  PR4B_GPU=$((GPU_I % 3)) launch_freemd "${jd}"
  GPU_I=$((GPU_I + 1))
done
wait
echo "[pr4b] free-MD done; injecting frames"

# 2) 帧注入 rep01/02/03
for jd in "${ROOT}/runs/real_cases/dssb_smoke_d32n/jobs/3hfm-patel-2021-antibody-h-d32n" \
          "${BATCH2}/jobs/3hfm-patel-2021-antigen-y-r21a" \
          "${BATCH2}/jobs/3hfm-patel-2021-antigen-y-d101k" \
          "${BATCH2}/jobs/3hfm-patel-2021-antibody-l-n31e" \
          "${BATCH2}/jobs/3hfm-patel-2021-antigen-y-k97d"; do
  eq="${jd}/legs/dssb/rep01/equilibration"
  [ -f "${eq}/freemd40.xtc" ] || { echo "[pr4b] no trajectory ${jd##*/}"; continue; }
  i=1
  for frame_ps in 12500 25000 37500; do
    target="${jd}/legs/dssb/rep0${i}/equilibration"
    printf 'System\nSystem\n' | "${GMX}" trjconv -f "${eq}/freemd40.xtc" -s "${eq}/freemd40.tpr" \
      -o "${target}/npt_seed.gro" -dump "${frame_ps}" -pbc mol -ur compact > /dev/null 2>&1
    if [ -s "${target}/npt_seed.gro" ]; then
      cp "${target}/npt.gro" "${target}/npt.gro.orig" 2>/dev/null || true
      cp "${target}/npt_seed.gro" "${target}/npt.gro"
      echo "[pr4b] rep0${i} ${jd##*/} <- ${frame_ps} ps"
    fi
    i=$((i+1))
  done
done

# 3) 清 sample 下游并 resume
"${ROOT}/.venv/bin/python" - <<'PYEOF'
import glob, os, shutil
jobs = [
    "runs/real_cases/dssb_smoke_d32n/jobs/3hfm-patel-2021-antibody-h-d32n",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-r21a",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-d101k",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antibody-l-n31e",
    "runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-k97d",
]
for jd in jobs:
    for pat in ["legs/dssb/rep*/lambda_*", "legs/dssb/rep*/bar"]:
        for p in glob.glob(f"{jd}/{pat}"):
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
    for s in ["build_legs", "sample", "bar", "qc", "report"]:
        # ISSUE-009: wiping lambda_* REQUIRES deleting build_legs.json too,
        # otherwise the stage is skipped and lambda mdp files never regenerate.
        f = f"{jd}/stages/{s}.json"
        if os.path.isfile(f):
            os.remove(f)
    for f in glob.glob(f"{jd}/results/*.json"):
        os.remove(f)
    print("wiped sample+:", jd.split("/")[-1])
PYEOF

for jd_job in "dssb_smoke_d32n:3hfm-patel-2021-antibody-h-d32n" \
              "3hfm_patel_dssb_pr4:3hfm-patel-2021-antigen-y-r21a" \
              "3hfm_patel_dssb_pr4:3hfm-patel-2021-antigen-y-d101k" \
              "3hfm_patel_dssb_pr4:3hfm-patel-2021-antibody-l-n31e" \
              "3hfm_patel_dssb_pr4:3hfm-patel-2021-antigen-y-k97d"; do
  batch="${jd_job%%:*}"; job="${jd_job##*:}"
  echo "[pr4b] FEP ${job}"
  "${ABAG_RBFE}" resume "${job}" --batch-dir "${ROOT}/runs/real_cases/${batch}" --execute 2>&1 | tail -1
done
echo "[pr4b] $(date --iso-8601=seconds) all done"
