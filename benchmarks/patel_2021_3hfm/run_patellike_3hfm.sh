#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/patel_2021_3hfm"
SYSTEM_PATH="${BENCHMARK_ROOT}/system.yml"
MUTATIONS_PATH="${BENCHMARK_ROOT}/mutations.csv"
PROTOCOL_PATH="${PROTOCOL_PATH:-${ROOT}/benchmarks/ab_bind/protocol.3hfm_patel2021_reference.yml}"
RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/patel_2021_3hfm}"
BATCH_ID="${BATCH_ID:-patel_2021_3hfm_reference}"
MAX_WORKERS="${MAX_WORKERS:-1}"
BATCH_DIR="${RUNS_ROOT}/${BATCH_ID}"

mkdir -p "${RUNS_ROOT}"

"${ABAG_RBFE}" batch plan \
  --system "${SYSTEM_PATH}" \
  --mutations "${MUTATIONS_PATH}" \
  --protocol "${PROTOCOL_PATH}" \
  --batch-id "${BATCH_ID}" \
  --runs-root "${RUNS_ROOT}" >/dev/null

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  mapfile -t JOB_IDS < <(BATCH_DIR_ENV="${BATCH_DIR}" python - <<'PY'
import json
import os
from pathlib import Path
batch_plan = Path(os.environ["BATCH_DIR_ENV"]) / "batch_plan.json"
payload = json.loads(batch_plan.read_text())
for job in payload["jobs"]:
    print(job["job_id"])
PY
)
fi

for job_id in "${JOB_IDS[@]}"; do
  "${ABAG_RBFE}" resume "${job_id}" --batch-dir "${BATCH_DIR}" --execute
done

"${ROOT}/benchmarks/patel_2021_3hfm/report_patellike_3hfm.py" \
  --batch-dir "${BATCH_DIR}"

echo "Finished. Inspect ${BATCH_DIR}/reports/patel_2021_3hfm_summary.json"
