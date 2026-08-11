#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_1jrh_runner_truequick"
JOB_ID="1jrh-antigen-i-m25l--i-i28v"

"${ROOT}/benchmarks/ab_bind/refresh_curated.sh" >/dev/null
"${ROOT}/benchmarks/ab_bind/materialize_inputs.sh" >/dev/null

"${ABAG_RBFE}" batch plan-abbind \
  --benchmark-root "${BENCHMARK_ROOT}" \
  --protocol "${BENCHMARK_ROOT}/protocol.quick.yml" \
  --spec core_v2 \
  --runs-root "${RUNS_ROOT}" \
  --complex-id 1JRH

"${ABAG_RBFE}" batch run-abbind \
  --plan-root "${RUNS_ROOT}" \
  --complex-id 1JRH \
  --job-id "${JOB_ID}" \
  --execute

"${ABAG_RBFE}" batch report-abbind \
  --plan-root "${RUNS_ROOT}" \
  --complex-id 1JRH

echo "Finished. Inspect ${RUNS_ROOT}/abbind_1jrh_core_v2/jobs/${JOB_ID}/results/ddg_summary.json"
