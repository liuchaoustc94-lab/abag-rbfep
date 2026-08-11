#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"
PROTOCOL_PATH="${BENCHMARK_ROOT}/protocol.validation_robust.yml"
WAIT_FOR_PID="${WAIT_FOR_PID:-}"
MAX_WORKERS="${MAX_WORKERS:-1}"
GPU_DEVICES="${GPU_DEVICES:-}"

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=(
    "1cz8-antigen-w-g92a"
    "3hfm-antibody-h-c95a"
    "3hfm-antibody-h-y33a"
    "3hfm-antibody-h-y50a"
  )
fi

"${ROOT}/benchmarks/ab_bind/refresh_curated.sh" >/dev/null
"${ROOT}/benchmarks/ab_bind/materialize_inputs.sh" >/dev/null

"${ABAG_RBFE}" batch plan-abbind \
  --benchmark-root "${BENCHMARK_ROOT}" \
  --protocol "${PROTOCOL_PATH}" \
  --spec core_v1 \
  --runs-root "${RUNS_ROOT}" \
  --split-name validation \
  --split-file "${SPLIT_FILE}"

RUN_ARGS=(
  batch run-abbind
  --plan-root "${RUNS_ROOT}"
  --split-name validation
  --split-file "${SPLIT_FILE}"
  --resume
  --execute
)

if [ -n "${GPU_DEVICES}" ]; then
  RUN_ARGS+=(--gpu-devices "${GPU_DEVICES}")
fi

if [ "${MAX_WORKERS}" -gt 1 ]; then
  RUN_ARGS+=(--max-workers "${MAX_WORKERS}")
fi

for job_id in "${JOB_IDS[@]}"; do
  RUN_ARGS+=(--job-id "${job_id}")
done

"${ABAG_RBFE}" "${RUN_ARGS[@]}"

"${ABAG_RBFE}" batch report-abbind \
  --plan-root "${RUNS_ROOT}" \
  --split-name validation \
  --split-file "${SPLIT_FILE}"

echo "Finished. Inspect ${RUNS_ROOT}/reports/selections/split-validation-complex-1bj1-1cz8-1mlc-2nz9-3hfm-3nps/benchmark_metrics.json"
