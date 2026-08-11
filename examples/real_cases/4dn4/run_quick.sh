#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
BATCH_ID="4dn4_v47i_quick"
RUNS_ROOT="${ROOT}/runs/real_cases"
SYSTEM="${ROOT}/examples/real_cases/4dn4/system.yml"
MUTATIONS="${ROOT}/examples/real_cases/4dn4/mutations.csv"
PROTOCOL="${ROOT}/examples/real_cases/4dn4/protocol.quick.yml"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"

mkdir -p "${RUNS_ROOT}"

"${ABAG_RBFE}" batch plan \
  --system "${SYSTEM}" \
  --mutations "${MUTATIONS}" \
  --protocol "${PROTOCOL}" \
  --batch-id "${BATCH_ID}" \
  --runs-root "${RUNS_ROOT}"

JOB_ID="$(tail -n +2 "${RUNS_ROOT}/${BATCH_ID}/jobs.csv" | cut -d, -f1)"

"${ABAG_RBFE}" run "${JOB_ID}" \
  --batch-dir "${RUNS_ROOT}/${BATCH_ID}" \
  --execute

"${ABAG_RBFE}" report "${RUNS_ROOT}/${BATCH_ID}"

echo "Finished. Inspect ${RUNS_ROOT}/${BATCH_ID}/reports/batch_summary.json"
