#!/usr/bin/env bash
# Stalled-work queue (2026-08-24). GPU capacity freed: 3 worker slots on GPU1/GPU2
# (GPU0 kept for the two running PR-4b 40ns free-MD jobs; GPU3 belongs to others).
#
# Contents:
#   Wave A (priority): 1BJ1 endpoint-ensemble FEP rerun (4 jobs)  — needs build_legs.json
#                       deletion so lambda mdp files regenerate (the classic wipe trap:
#                       lambda dirs were wiped on 08-19 but build_legs.json survived, so
#                       sample failed with missing pre_relax.mdp on 08-24).
#                       E4 tail: w98f/y50l sample resume (per-window idempotent skip).
#                       h487q: 24-lambda retry with rlist 1.25->1.50 patched into
#                       pre_relax/pre_md mdp files (mitigation for the exclusionchecker
#                       fatal at complex/rep01/lambda_011; production.mdp untouched).
#   Wave B: acceptance V1 remaining jobs (2bdn 2, 1nmb 6, 1ahw 8, 1iar 11, 1oga 31,
#           4i77 leftover 5) that were planned 08-18 but never executed.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
export PATH="${ROOT}/.venv/bin:${PATH}"
LOG="${ROOT}/runs/manual/pool_stalled_20260824.log"
QDIR="${ROOT}/runs/manual/stalled_queue_20260824"
mkdir -p "${QDIR}"

log() { echo "[pool] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

# ---------------------------------------------------------------- Phase 0: fixes
BJ1_DIR="${ROOT}/runs/benchmarks/abbind_validation_panels_20260812/abbind_1bj1_core_v1"
for job in 1bj1-antigen-w-g92a 1bj1-antigen-w-g88a 1bj1-antigen-w-q89a 1bj1-antigen-w-i91a; do
  jd="${BJ1_DIR}/jobs/${job}"
  rm -f "${jd}/stages/build_legs.json" "${jd}/stages/sample.json" \
        "${jd}/stages/bar.json" "${jd}/stages/qc.json" "${jd}/stages/report.json"
  log "fix 1bj1 ${job}: removed stale build_legs/sample stage state"
done

E4_DIR="${ROOT}/runs/real_cases/3hfm_patel_align_20260817"
for job in 3hfm-patel-2021-antibody-h-w98f 3hfm-patel-2021-antibody-l-y50l; do
  jd="${E4_DIR}/jobs/${job}"
  rm -f "${jd}/stages/sample.json" "${jd}/stages/bar.json" "${jd}/stages/qc.json" "${jd}/stages/report.json"
  log "fix e4 ${job}: cleared failed sample stage (completed windows are skipped idempotently)"
done

H487Q_DIR="${ROOT}/runs/real_cases/oos_1ak4_newprotocol_20260811"
H487Q_JOB="1ak4-out-of-sample-antigen-d-h487q"
jd="${H487Q_DIR}/jobs/${H487Q_JOB}"
patched=0
for mdp in "${jd}"/legs/*/rep*/lambda_*/pre_relax.mdp "${jd}"/legs/*/rep*/lambda_*/pre_md.mdp; do
  [ -f "${mdp}" ] || continue
  sed -i 's/^rlist\s*=\s*1\.25/rlist                   = 1.50/' "${mdp}"
  patched=$((patched+1))
done
rm -f "${jd}/stages/sample.json" "${jd}/stages/bar.json" "${jd}/stages/qc.json" "${jd}/stages/report.json"
log "fix h487q: patched rlist->1.50 in ${patched} pre_relax/pre_md mdp files, cleared sample stage"

# ---------------------------------------------------------------- Phase 1: queue
QUEUE="${QDIR}/queue.txt"
: > "${QUEUE}"
for job in 1bj1-antigen-w-g92a 1bj1-antigen-w-g88a 1bj1-antigen-w-q89a 1bj1-antigen-w-i91a; do
  echo "${BJ1_DIR}|${job}" >> "${QUEUE}"
done
for job in 3hfm-patel-2021-antibody-h-w98f 3hfm-patel-2021-antibody-l-y50l; do
  echo "${E4_DIR}|${job}" >> "${QUEUE}"
done
echo "${H487Q_DIR}|${H487Q_JOB}" >> "${QUEUE}"

for acc in acc_1bj1_skip acc_2bdn_v1_20260818 acc_1nmb_v1_20260818 acc_1ahw_v1_20260818 acc_1iar_v1_20260818 acc_1oga_v1_20260818 acc_4i77_v1_20260818; do
  [ "${acc}" = "acc_1bj1_skip" ] && continue
  bd="${ROOT}/runs/acceptance/${acc}"
  [ -d "${bd}/jobs" ] || continue
  for jd in "${bd}/jobs"/*/; do
    job="$(basename "${jd}")"
    [ -f "${jd}/results/ddg_summary.json" ] && continue
    echo "${bd}|${job}" >> "${QUEUE}"
  done
done
total=$(wc -l < "${QUEUE}")
log "queue built: ${total} jobs"

# ---------------------------------------------------------------- Phase 2: workers
# Shared queue with flock-atomic pop; 8 workers: gpu0 x1, gpu1 x3, gpu2 x4.
# (gpu0 also carries two PR-4b 40ns free-MD jobs; gpu1 carries one d32n lambda job;
#  gpu3 belongs to other users and stays untouched.)
GPU_MAP=(0 1 1 1 2 2 2 2)
NW=${#GPU_MAP[@]}

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
    log "[w${wid}][gpu${gpu}] start ${job} ($(basename "${batch_dir}"))"
    ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" \
      --batch-dir "${batch_dir}" --execute > "${QDIR}/w${wid}_last.log" 2>&1
    rc=$?
    if [ -f "${batch_dir}/jobs/${job}/results/ddg_summary.json" ]; then
      log "[w${wid}][gpu${gpu}] done ${job} rc=${rc} ddg=READY"
    else
      log "[w${wid}][gpu${gpu}] done ${job} rc=${rc} ddg=MISSING (see ${QDIR}/w${wid}_last.log)"
    fi
  done
  log "[w${wid}][gpu${gpu}] queue drained"
}

for wid in "${!GPU_MAP[@]}"; do
  worker "${wid}" "${GPU_MAP[${wid}]}" &
done
wait
log "ALL DONE"
