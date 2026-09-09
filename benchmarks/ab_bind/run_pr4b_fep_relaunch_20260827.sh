#!/usr/bin/env bash
# PR-4b FEP relaunch (2026-08-27): the original orchestrator woke up when d32n
# free-MD finished, wiped lambda dirs mid-flight (racing the unblock run), then
# exited leaving all 5 jobs with stale build_legs.json and no lambda dirs.
# This script does the correct full wipe (incl. build_legs.json, ISSUE-009) and
# relaunches all 5 jobs on idle GPUs.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
export PATH="${ROOT}/.venv/bin:${PATH}"
LOG="${ROOT}/runs/manual/pr4b_relaunch_20260827.log"
log() { echo "[pr4b-relaunch] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

declare -A BATCH_OF=(
  [3hfm-patel-2021-antibody-h-d32n]="dssb_smoke_d32n"
  [3hfm-patel-2021-antigen-y-r21a]="3hfm_patel_dssb_pr4"
  [3hfm-patel-2021-antigen-y-d101k]="3hfm_patel_dssb_pr4"
  [3hfm-patel-2021-antibody-l-n31e]="3hfm_patel_dssb_pr4"
  [3hfm-patel-2021-antigen-y-k97d]="3hfm_patel_dssb_pr4"
)

for job in "${!BATCH_OF[@]}"; do
  jd="${ROOT}/runs/real_cases/${BATCH_OF[$job]}/jobs/${job}"
  # full wipe: lambda dirs + bar + all downstream stage state + results
  find "${jd}/legs" -maxdepth 3 -name "lambda_*" -type d -exec rm -rf {} + 2>/dev/null
  find "${jd}/legs" -maxdepth 3 -name "bar" -type d -exec rm -rf {} + 2>/dev/null
  rm -f "${jd}"/stages/{build_legs,sample,bar,qc,report}.json
  rm -f "${jd}"/results/*.json
  # sanity: frames must be injected already (npt_seed.gro per rep)
  nseed=$(find "${jd}/legs" -name "npt_seed.gro" | wc -l)
  log "wiped ${job} (npt_seed=${nseed}/3)"
done

GPUS=(0 1 2 1 2)
i=0
for job in 3hfm-patel-2021-antigen-y-r21a 3hfm-patel-2021-antigen-y-d101k 3hfm-patel-2021-antibody-l-n31e 3hfm-patel-2021-antigen-y-k97d 3hfm-patel-2021-antibody-h-d32n; do
  gpu="${GPUS[$i]}"
  log "FEP ${job} on gpu${gpu}"
  ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" \
    --batch-dir "${ROOT}/runs/real_cases/${BATCH_OF[$job]}" --execute \
    > "${ROOT}/runs/manual/pr4b_relaunch_${job}.log" 2>&1 &
  i=$((i+1))
done
wait
log "ALL DONE"
