#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
RUNS_ROOT="${ROOT}/runs/benchmarks"
BATCH_ID="abbind_1vfb_core_v1_quick"
JOB_ID="1vfb-antibody-h-y32a"

"${ROOT}/benchmarks/ab_bind/refresh_curated.sh" >/dev/null
"${ROOT}/benchmarks/ab_bind/materialize_inputs.sh" >/dev/null

"${ABAG_RBFE}" batch plan \
  --system "${BENCHMARK_ROOT}/materialized/1VFB/system.yml" \
  --mutations "${BENCHMARK_ROOT}/materialized/1VFB/core_v1_mutations.csv" \
  --protocol "${BENCHMARK_ROOT}/protocol.quick.yml" \
  --batch-id "${BATCH_ID}" \
  --runs-root "${RUNS_ROOT}"

"${ABAG_RBFE}" run "${JOB_ID}" \
  --batch-dir "${RUNS_ROOT}/${BATCH_ID}" \
  --execute

"${ABAG_RBFE}" report "${RUNS_ROOT}/${BATCH_ID}"

echo "Finished. Inspect ${RUNS_ROOT}/${BATCH_ID}/jobs/${JOB_ID}/results/ddg_summary.json"
