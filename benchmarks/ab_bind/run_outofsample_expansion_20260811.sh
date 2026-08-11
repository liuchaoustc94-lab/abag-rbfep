#!/usr/bin/env bash
# Out-of-sample expansion runner (2026-08-11): 1AK4(15) + 1KTZ(15) + 3K2M(5) = 35 jobs.
# Waits for Phase B (run_phase_b_newprotocol_20260811) to finish, then runs all
# three cases under the new protocol (adaptive lambda + decoupled schedule + dt 2fs).
# Usage: nohup bash benchmarks/ab_bind/run_outofsample_expansion_20260811.sh > log 2>&1 &
set -uo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
PROTOCOL="${ROOT}/benchmarks/ab_bind/protocol.validation_priority.yml"
PHASE_B_PATTERN="run_phase_b_newprotocol_20260811"
export PATH="${ROOT}/.venv/bin:${PATH}"

echo "[oos] $(date --iso-8601=seconds) waiting for Phase B..."
while pgrep -f "${PHASE_B_PATTERN}" > /dev/null; do sleep 120; done
echo "[oos] $(date --iso-8601=seconds) Phase B done; starting out-of-sample batches"

for case in 1ak4 1ktz 3k2m; do
  batch="oos_${case}_newprotocol_20260811"
  batch_dir="${ROOT}/runs/real_cases/${batch}"
  "${ABAG_RBFE}" batch plan \
    --system "${ROOT}/examples/real_cases/${case}/system.yml" \
    --mutations "${ROOT}/examples/real_cases/${case}/mutations.csv" \
    --protocol "${PROTOCOL}" \
    --batch-id "${batch}" \
    --runs-root "${ROOT}/runs/real_cases"
  mapfile -t jobs < <(ls "${batch_dir}/jobs")
  half=$(( (${#jobs[@]} + 1) / 2 ))
  for i in "${!jobs[@]}"; do
    gpu=$(( i < half ? 0 : 1 ))
    echo "[oos][gpu${gpu}] $(date --iso-8601=seconds) run ${jobs[$i]}"
    ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${jobs[$i]}" --batch-dir "${batch_dir}" --execute 2>&1 | tail -1 &
    # 每 GPU 串行、双 GPU 并行：当后台任务数达到 2 时等待任意一个完成
    while [ "$(jobs -r | wc -l)" -ge 2 ]; do sleep 30; done
  done
  wait
done

echo "[oos] $(date --iso-8601=seconds) all done"
