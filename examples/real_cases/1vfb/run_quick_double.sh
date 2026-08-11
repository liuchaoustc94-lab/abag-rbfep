#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
cd "${ROOT}"

RUNS_ROOT="${ROOT}/runs/real_cases"
SYSTEM="${ROOT}/examples/real_cases/1vfb/system.yml"
MUTATIONS="${ROOT}/examples/real_cases/1vfb/mutations.double.csv"
PROTOCOL="${ROOT}/examples/real_cases/1vfb/protocol.quick.double.yml"
BATCH_ID="1vfb_y32f_v34i_quick"

abag-rbfe mutation validate \
  --mutations "${MUTATIONS}" \
  --output "${ROOT}/tmp/1vfb_double_validate.json"

abag-rbfe batch plan \
  --system "${SYSTEM}" \
  --mutations "${MUTATIONS}" \
  --protocol "${PROTOCOL}" \
  --batch-id "${BATCH_ID}" \
  --runs-root "${RUNS_ROOT}"

JOB_ID="$(tail -n +2 "${RUNS_ROOT}/${BATCH_ID}/jobs.csv" | cut -d, -f1)"

abag-rbfe run \
  "${JOB_ID}" \
  --batch-dir "${RUNS_ROOT}/${BATCH_ID}" \
  --execute
