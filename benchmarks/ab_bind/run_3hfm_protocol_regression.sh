#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_3hfm_protocol_regression}"
PROTOCOL_PATH="${PROTOCOL_PATH:-${BENCHMARK_ROOT}/protocol.validation_robust.yml}"
COMPLEX_ID="${COMPLEX_ID:-3HFM}"
MAX_WORKERS="${MAX_WORKERS:-1}"
GPU_DEVICES="${GPU_DEVICES:-}"

"${ROOT}/benchmarks/ab_bind/refresh_curated.sh" >/dev/null
"${ROOT}/benchmarks/ab_bind/materialize_inputs.sh" >/dev/null

"${ABAG_RBFE}" batch plan-abbind \
  --benchmark-root "${BENCHMARK_ROOT}" \
  --protocol "${PROTOCOL_PATH}" \
  --spec core_v1 \
  --runs-root "${RUNS_ROOT}" \
  --complex-id "${COMPLEX_ID}"

RUN_ARGS=(
  batch run-abbind
  --plan-root "${RUNS_ROOT}"
  --complex-id "${COMPLEX_ID}"
  --resume
  --execute
)

if [ -n "${GPU_DEVICES}" ]; then
  RUN_ARGS+=(--gpu-devices "${GPU_DEVICES}")
fi

if [ "${MAX_WORKERS}" -gt 1 ]; then
  RUN_ARGS+=(--max-workers "${MAX_WORKERS}")
fi

if [ "$#" -gt 0 ]; then
  for job_id in "$@"; do
    RUN_ARGS+=(--job-id "${job_id}")
  done
fi

"${ABAG_RBFE}" "${RUN_ARGS[@]}"

"${ROOT}/benchmarks/ab_bind/report_3hfm_protocol_regression.py" \
  --plan-root "${RUNS_ROOT}" \
  --complex-id "${COMPLEX_ID}"

echo "Finished. Inspect ${RUNS_ROOT}/reports/3hfm_protocol_regression_summary.json"
