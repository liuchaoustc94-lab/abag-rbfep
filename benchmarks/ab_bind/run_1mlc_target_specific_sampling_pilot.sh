#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
PROTOCOL_PATH="${PROTOCOL_PATH:-${BENCHMARK_ROOT}/protocol.1mlc_target_specific_sampling_pilot.yml}"
RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_1mlc_20260626}"
BATCH_PREFIX="${BATCH_PREFIX:-abbind-target-specific-sampling-pilot-1mlc}"
COMPLEX_ID="${COMPLEX_ID:-1MLC}"
MAX_WORKERS="${MAX_WORKERS:-2}"
GPU_DEVICES="${GPU_DEVICES:-0,1}"

DEFAULT_JOB_IDS=(
  "1mlc-antibody-l-n92a"
  "1mlc-antibody-h-s57v"
  "1mlc-antibody-h-t31a"
)

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=("${DEFAULT_JOB_IDS[@]}")
fi

PLAN_ARGS=(
  batch plan-abbind
  --benchmark-root "${BENCHMARK_ROOT}"
  --protocol "${PROTOCOL_PATH}"
  --spec core_v1
  --runs-root "${RUNS_ROOT}"
  --batch-prefix "${BATCH_PREFIX}"
  --complex-id "${COMPLEX_ID}"
)

"${ABAG_RBFE}" "${PLAN_ARGS[@]}"

RUN_ARGS=(
  batch run-abbind
  --plan-root "${RUNS_ROOT}"
  --complex-id "${COMPLEX_ID}"
  --resume
  --execute
  --max-workers "${MAX_WORKERS}"
)

if [ -n "${GPU_DEVICES}" ]; then
  RUN_ARGS+=(--gpu-devices "${GPU_DEVICES}")
fi

for job_id in "${JOB_IDS[@]}"; do
  RUN_ARGS+=(--job-id "${job_id}")
done

"${ABAG_RBFE}" "${RUN_ARGS[@]}"

"${ABAG_RBFE}" batch report-abbind \
  --plan-root "${RUNS_ROOT}" \
  --complex-id "${COMPLEX_ID}" \
  --split-name validation \
  --split-file "${SPLIT_FILE}"

"${ROOT}/.venv/bin/python" "${BENCHMARK_ROOT}/report_hotspot_root_comparison.py" \
  --job-id 1mlc-antibody-l-n92a \
  --job-id 1mlc-antibody-l-n32g \
  --job-id 1mlc-antibody-l-n32y \
  --plan-root "priority=${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan" \
  --plan-root "robust=${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan" \
  --plan-root "sampling_qc=${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues" \
  --plan-root "deep=${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues" \
  --plan-root "ultra=${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues" \
  --plan-root "pilot_1mlc=${RUNS_ROOT}" \
  --json-output "${RUNS_ROOT}/reports/hotspot_root_comparison.json" \
  --md-output "${RUNS_ROOT}/reports/hotspot_root_comparison.md" \
  || true

echo "Finished. Inspect ${RUNS_ROOT}/reports/plan_summary.json"
