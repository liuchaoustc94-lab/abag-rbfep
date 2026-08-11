#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_quick_plan"

COMPLEX_ID="${COMPLEX_ID:-1MLC}"

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=(
    "1mlc-antibody-h-s57a"
    "1mlc-antibody-h-t31a"
    "1mlc-antibody-h-t31v"
    "1mlc-antibody-l-n92a"
  )
fi

"${ROOT}/benchmarks/ab_bind/refresh_curated.sh" >/dev/null
"${ROOT}/benchmarks/ab_bind/materialize_inputs.sh" >/dev/null

"${ABAG_RBFE}" batch plan-abbind \
  --benchmark-root "${BENCHMARK_ROOT}" \
  --protocol "${BENCHMARK_ROOT}/protocol.quick.yml" \
  --spec core_v1 \
  --runs-root "${RUNS_ROOT}" \
  --split-name validation \
  --split-file "${SPLIT_FILE}"

for job_id in "${JOB_IDS[@]}"; do
  "${ABAG_RBFE}" batch run-abbind \
    --plan-root "${RUNS_ROOT}" \
    --split-name validation \
    --split-file "${SPLIT_FILE}" \
    --complex-id "${COMPLEX_ID}" \
    --job-id "${job_id}" \
    --execute
done

"${ABAG_RBFE}" batch report-abbind \
  --plan-root "${RUNS_ROOT}" \
  --split-name validation \
  --split-file "${SPLIT_FILE}"

echo "Finished. Inspect ${RUNS_ROOT}/reports/selections/split-validation-complex-1bj1-1cz8-1mlc-2nz9-3hfm-3nps/benchmark_metrics.json"
