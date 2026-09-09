#!/usr/bin/env bash
# Retry-ladder round 2 (2026-08-26): requeue the round-1 failures with the
# 4-rung ladder (standard -> reseed -> couple-intramol -> ci+half-dt robust).
# Round 1 (3-rung) showed the dominant residual failure is the GROMACS
# exclusionchecker fatal ("perturbed excluded pairs beyond pair-list cutoff"),
# whose official remedy is couple-intramol=yes (now rung 3).
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
export PATH="${ROOT}/.venv/bin:${PATH}"
LOG="${ROOT}/runs/manual/pool_retry_ladder2_20260826.log"
QDIR="${ROOT}/runs/manual/retry_ladder2_20260826"
mkdir -p "${QDIR}"
log() { echo "[pool] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

# wait for round-1 pool to fully drain (locks would block us anyway)
while pgrep -f "run_retry_ladder_rescue_20260826" > /dev/null; do sleep 60; done

E4="${ROOT}/runs/real_cases/3hfm_patel_align_20260817"
QUEUE="${QDIR}/queue.txt"
cat > "${QUEUE}" <<EOF
${E4}|3hfm-patel-2021-antibody-h-w98f
${E4}|3hfm-patel-2021-antibody-l-y50l
${ROOT}/runs/real_cases/oos_1ak4_newprotocol_20260811|1ak4-out-of-sample-antigen-d-h487q
${ROOT}/runs/acceptance/acc_4i77_v1_20260818|acc-4i77-antibody-h-w52f
${ROOT}/runs/acceptance/acc_2bdn_v1_20260818|acc-2bdn-antibody-h-d31e
${ROOT}/runs/acceptance/acc_2bdn_v1_20260818|acc-2bdn-antibody-h-n28q
${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-h-d56e
${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-l-l94v
${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-l-t93f
${ROOT}/runs/acceptance/acc_1nmb_v1_20260818|acc-1nmb-antibody-l-t93w
${ROOT}/runs/acceptance/acc_1oga_v1_20260818|acc-1oga-antibody-e-i53l
${ROOT}/runs/acceptance/acc_1oga_v1_20260818|acc-1oga-antibody-e-i53n
EOF

# requeue only jobs that still lack a result; clear sample stage again
TMPQ="${QDIR}/queue.active.txt"; : > "${TMPQ}"
while IFS='|' read -r batch_dir job; do
  jd="${batch_dir}/jobs/${job}"
  if [ -f "${jd}/results/ddg_summary.json" ]; then
    log "skip ${job}: already has ddg"
    continue
  fi
  rm -f "${jd}/stages/sample.json" "${jd}/stages/bar.json" "${jd}/stages/qc.json" "${jd}/stages/report.json"
  echo "${batch_dir}|${job}" >> "${TMPQ}"
  log "requeue ${job}"
done < "${QUEUE}"
mv "${TMPQ}" "${QUEUE}"
n=$(wc -l < "${QUEUE}")
log "round-2 queue: ${n} jobs"
[ "${n}" -eq 0 ] && { log "nothing to do"; exit 0; }

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
