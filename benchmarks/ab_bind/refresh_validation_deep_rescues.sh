#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
SOURCE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
ROBUST_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"
RESCUE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"
TARGETED_REPEAT_SPREAD_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"
SAMPLING_QC_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"
BATCH_PREFIX="${BATCH_PREFIX:-abbind-deep-rescue}"
REPEAT_INCREMENT="${REPEAT_INCREMENT:-1}"
LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-4}"
PRODUCTION_SCALE="${PRODUCTION_SCALE:-4.0}"
WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-2.0}"
WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-2.0}"
NVT_SCALE="${NVT_SCALE:-2.0}"
NPT_SCALE="${NPT_SCALE:-2.0}"
REQUIRE_ACTIVE_ALTERNATE="${REQUIRE_ACTIVE_ALTERNATE:-0}"
HOTSPOT_COMPLEX_IDS="${HOTSPOT_COMPLEX_IDS:-3HFM}"
HOTSPOT_REPEAT_INCREMENT="${HOTSPOT_REPEAT_INCREMENT:-2}"
HOTSPOT_LAMBDA_INCREMENT="${HOTSPOT_LAMBDA_INCREMENT:-6}"
HOTSPOT_PRODUCTION_SCALE="${HOTSPOT_PRODUCTION_SCALE:-6.0}"
HOTSPOT_WINDOW_RELAX_EM_SCALE="${HOTSPOT_WINDOW_RELAX_EM_SCALE:-3.0}"
HOTSPOT_WINDOW_RELAX_MD_SCALE="${HOTSPOT_WINDOW_RELAX_MD_SCALE:-3.0}"
HOTSPOT_NVT_SCALE="${HOTSPOT_NVT_SCALE:-3.0}"
HOTSPOT_NPT_SCALE="${HOTSPOT_NPT_SCALE:-3.0}"

CMD=(
  "${ABAG_RBFE}"
  batch
  rescue-abbind
  --plan-root "${SOURCE_PLAN_ROOT}"
  --extra-plan-root "${TARGETED_REPEAT_SPREAD_PLAN_ROOT}"
  --extra-plan-root "${SAMPLING_QC_PLAN_ROOT}"
  --extra-plan-root "${ROBUST_PLAN_ROOT}"
  --extra-plan-root "${RESCUE_PLAN_ROOT}"
  --split-name validation
  --split-file "${SPLIT_FILE}"
  --runs-root "${RUNS_ROOT}"
  --batch-prefix "${BATCH_PREFIX}"
  --repeat-increment "${REPEAT_INCREMENT}"
  --lambda-increment "${LAMBDA_INCREMENT}"
  --production-scale "${PRODUCTION_SCALE}"
  --window-relax-em-scale "${WINDOW_RELAX_EM_SCALE}"
  --window-relax-md-scale "${WINDOW_RELAX_MD_SCALE}"
  --nvt-scale "${NVT_SCALE}"
  --npt-scale "${NPT_SCALE}"
  --force-repeat-increment
  --prefer-active-alternate-source
  --target-primary-repeat-spread-leg
  --allow-targeted-leg-count-deepening
)

if [ "${REQUIRE_ACTIVE_ALTERNATE}" = "1" ]; then
  CMD+=(--require-active-alternate)
fi

if [ -n "${HOTSPOT_COMPLEX_IDS}" ]; then
  IFS=',' read -r -a HOTSPOT_COMPLEX_ID_VALUES <<< "${HOTSPOT_COMPLEX_IDS}"
  for complex_id in "${HOTSPOT_COMPLEX_ID_VALUES[@]}"; do
    complex_id="${complex_id//[[:space:]]/}"
    if [ -n "${complex_id}" ]; then
      CMD+=(--hotspot-complex-id "${complex_id}")
    fi
  done
  CMD+=(
    --hotspot-repeat-increment "${HOTSPOT_REPEAT_INCREMENT}"
    --hotspot-lambda-increment "${HOTSPOT_LAMBDA_INCREMENT}"
    --hotspot-production-scale "${HOTSPOT_PRODUCTION_SCALE}"
    --hotspot-window-relax-em-scale "${HOTSPOT_WINDOW_RELAX_EM_SCALE}"
    --hotspot-window-relax-md-scale "${HOTSPOT_WINDOW_RELAX_MD_SCALE}"
    --hotspot-nvt-scale "${HOTSPOT_NVT_SCALE}"
    --hotspot-npt-scale "${HOTSPOT_NPT_SCALE}"
  )
fi

if [ "$#" -gt 0 ]; then
  for job_id in "$@"; do
    CMD+=(--job-id "${job_id}")
  done
fi

exec "${CMD[@]}"
