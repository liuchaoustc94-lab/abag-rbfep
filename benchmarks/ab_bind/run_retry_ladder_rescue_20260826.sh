#!/usr/bin/env bash
# ISSUE-007 window-crash rescue queue (2026-08-26): rerun the 12 MD-crash jobs
# with the new per-window retry ladder (standard -> reseed -> half-dt robust).
# QC-blocked jobs (1OGA crystal defects: n55a/q52a/q155a/v67a/v76a) are
# intentionally NOT requeued — the blocking is correct behavior.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
export PATH="${ROOT}/.venv/bin:${PATH}"
LOG="${ROOT}/runs/manual/pool_retry_ladder_20260826.log"
QDIR="${ROOT}/runs/manual/retry_ladder_20260826"
mkdir -p "${QDIR}"

log() { echo "[pool] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

QUEUE="${QDIR}/queue.txt"
: > "${QUEUE}"

E4="${ROOT}/runs/real_cases/3hfm_patel_align_20260817"
echo "${E4}|3hfm-patel-2021-antibody-h-w98f" >> "${QUEUE}"
echo "${E4}|3hfm-patel-2021-antibody-l-y50l" >> "${QUEUE}"
echo "${ROOT}/runs/real_cases/oos_1ak4_newprotocol_20260811|1ak4-out-of-sample-antigen-d-h487q" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_4i77_v1_20260818|acc-4i77-antibody-h-w52f" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_2bdn_v1_20260818|acc-2bdn-antibody-h-d31e" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_2bdn_v1_20260818|acc-2bdn-antibody-h-n28q" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-h-d56e" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-l-l94v" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-l-t93f" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-l-t93w" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_1oga_v1_20260818|acc-1oga-antibody-e-i53l" >> "${QUEUE}"
echo "${ROOT}/runs/acceptance/acc_1oga_v1_20260818|acc-1oga-antibody-e-i53n" >> "${QUEUE}"

# clear the failed sample stage so the ladder-equipped sample stage re-runs
while IFS='|' read -r batch_dir job; do
  jd="${batch_dir}/jobs/${job}"
  rm -f "${jd}/stages/sample.json" "${jd}/stages/bar.json" "${jd}/stages/qc.json" "${jd}/stages/report.json"
  log "cleared sample stage: ${job}"
done < "${QUEUE}"

pop_job() {
  flock "${QUEUE}.lock" bash -c '
    q="$1"
    [ -s "$q" ] || exit 0
    head -1 "$q"
    tail -n +2 "$q" > "$q.tmp" && mv "$q.tmp" "$q"
  ' _ "${QUEUE}"
}

worker() {
  local wid="$1" gpu="$2"
  while true; do
    local line
    line="$(pop_job)"
    [ -n "${line}" ] || break
    local batch_dir="${line%%|*}" job="${line##*|}"
    log "[w${wid}][gpu${gpu}] start ${job}"
    ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" \
      --batch-dir "${batch_dir}" --execute > "${QDIR}/w${wid}_${job}.log" 2>&1
    rc=$?
    if [ -f "${batch_dir}/jobs/${job}/results/ddg_summary.json" ]; then
      log "[w${wid}][gpu${gpu}] done ${job} rc=${rc} ddg=READY"
    else
      log "[w${wid}][gpu${gpu}] done ${job} rc=${rc} ddg=MISSING"
    fi
  done
  log "[w${wid}][gpu${gpu}] queue drained"
}

worker 0 1 &
worker 1 2 &
worker 2 2 &
wait
log "ALL DONE"
