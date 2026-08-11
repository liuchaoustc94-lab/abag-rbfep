#!/usr/bin/env bash
# Unified parallel pool runner (2026-08-11): all pending jobs across batches,
# 6 workers = 3 GPUs x 2 concurrent jobs each (CPU: 64 cores, mdrun -ntomp 4).
# Queue lines: "<batch_dir>|<job_id>"
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
QUEUE="${ROOT}/runs/manual/parallel_queue_20260811.txt"
LOCK="${QUEUE}.lock"
DONE_LOG="${ROOT}/runs/manual/parallel_pool_20260811.done.log"
export PATH="${ROOT}/.venv/bin:${PATH}"

worker() {
  local gpu="$1"
  while true; do
    local line
    line=$(flock "${LOCK}" -c "sed -n '1p' '${QUEUE}' && sed -i '1d' '${QUEUE}'" 2>/dev/null | head -1)
    [ -z "${line}" ] && break
    local batch_dir="${line%%|*}"
    local job_id="${line##*|}"
    echo "[pool][gpu${gpu}] $(date --iso-8601=seconds) start ${job_id} (${batch_dir##*/})"
    ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job_id}" --batch-dir "${batch_dir}" --execute >> "${DONE_LOG}" 2>&1
    echo "[pool][gpu${gpu}] $(date --iso-8601=seconds) done ${job_id} rc=$?" | tee -a "${DONE_LOG}"
  done
}

mkdir -p "${ROOT}/runs/manual"
: > "${DONE_LOG}"
worker 0 & worker 0 & worker 1 & worker 1 & worker 2 & worker 2 &
wait
echo "[pool] $(date --iso-8601=seconds) queue drained"
