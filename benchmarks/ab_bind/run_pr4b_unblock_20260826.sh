#!/usr/bin/env bash
# PR-4b unblock (2026-08-26): the original run_pr4b_charge_ensemble_20260819.sh
# is stuck at `wait` behind d32n's 40ns free-MD (81%, ETA ~1.5d, still healthy).
# The 4 real PR-4 mutants finished their 40ns runs on 08-24/25, so inject frames
# and launch their FEP now. Also fixes the ISSUE-009 trap: the original script's
# wipe forgot stages/build_legs.json (lambda mdp files would never regenerate).
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
BATCH2="${ROOT}/runs/real_cases/3hfm_patel_dssb_pr4"
export PATH="${ROOT}/.venv/bin:${PATH}"
LOG="${ROOT}/runs/manual/pr4b_unblock_20260826.log"

log() { echo "[pr4b-unblock] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

JOBS=(3hfm-patel-2021-antigen-y-r21a 3hfm-patel-2021-antigen-y-d101k 3hfm-patel-2021-antibody-l-n31e 3hfm-patel-2021-antigen-y-k97d)

# 1) frame injection rep01/02/03 from the finished 40ns trajectories
for job in "${JOBS[@]}"; do
  jd="${BATCH2}/jobs/${job}"
  eq="${jd}/legs/dssb/rep01/equilibration"
  [ -f "${eq}/freemd40.xtc" ] || { log "skip ${job}: no freemd40.xtc"; continue; }
  i=1
  for frame_ps in 12500 25000 37500; do
    target="${jd}/legs/dssb/rep0${i}/equilibration"
    printf 'System\nSystem\n' | "${GMX}" trjconv -f "${eq}/freemd40.xtc" -s "${eq}/freemd40.tpr" \
      -o "${target}/npt_seed.gro" -dump "${frame_ps}" -pbc mol -ur compact > /dev/null 2>&1
    if [ -s "${target}/npt_seed.gro" ]; then
      cp "${target}/npt.gro" "${target}/npt.gro.orig" 2>/dev/null || true
      cp "${target}/npt_seed.gro" "${target}/npt.gro"
      log "injected ${job} rep0${i} <- ${frame_ps} ps"
    else
      log "WARN ${job} rep0${i} injection failed"
    fi
    i=$((i+1))
  done
done

# 2) wipe sample downstream INCLUDING build_legs.json (ISSUE-009 lesson)
"${ROOT}/.venv/bin/python" - <<'PYEOF'
import glob, os, shutil
base = "runs/real_cases/3hfm_patel_dssb_pr4/jobs"
jobs = ["3hfm-patel-2021-antigen-y-r21a", "3hfm-patel-2021-antigen-y-d101k",
        "3hfm-patel-2021-antibody-l-n31e", "3hfm-patel-2021-antigen-y-k97d"]
for job in jobs:
    jd = f"{base}/{job}"
    for pat in ["legs/dssb/rep*/lambda_*", "legs/dssb/rep*/bar"]:
        for p in glob.glob(f"{jd}/{pat}"):
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
    for s in ["build_legs", "sample", "bar", "qc", "report"]:
        f = f"{jd}/stages/{s}.json"
        if os.path.isfile(f):
            os.remove(f)
    for f in glob.glob(f"{jd}/results/*.json"):
        os.remove(f)
    print("wiped(build_legs+sample+)", job)
PYEOF

# 3) resume FEP on GPU1/GPU2 (GPU0 carries d32n freemd40)
i=0
for job in "${JOBS[@]}"; do
  gpu=$(( i % 2 + 1 ))
  log "FEP ${job} on gpu${gpu}"
  ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" --batch-dir "${BATCH2}" --execute \
    > "${ROOT}/runs/manual/pr4b_unblock_${job}.log" 2>&1 &
  i=$((i+1))
done
wait
log "ALL DONE"
