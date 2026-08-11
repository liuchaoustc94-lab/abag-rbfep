#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
SOURCE_PLAN_ROOT="${SOURCE_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues}"
RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues}"
BATCH_PREFIX="${BATCH_PREFIX:-abbind-targeted-lambda-rescue}"
REPEAT_INCREMENT="${REPEAT_INCREMENT:-0}"
LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-4}"
PRODUCTION_SCALE="${PRODUCTION_SCALE:-1.0}"
WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-1.0}"
WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-1.0}"
NVT_SCALE="${NVT_SCALE:-1.0}"
NPT_SCALE="${NPT_SCALE:-1.0}"
FORCE_LAMBDA_INCREMENT="${FORCE_LAMBDA_INCREMENT:-1}"
ALLOW_TARGETED_LEG_COUNT_DEEPENING="${ALLOW_TARGETED_LEG_COUNT_DEEPENING:-1}"

CMD=(
  "${ABAG_RBFE}"
  batch
  rescue-abbind
  --plan-root "${SOURCE_PLAN_ROOT}"
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
  --target-primary-repeat-spread-leg
  --require-target-primary-repeat-spread-leg
)

if [ "${FORCE_LAMBDA_INCREMENT}" = "1" ]; then
  CMD+=(--force-lambda-increment)
fi

if [ "${ALLOW_TARGETED_LEG_COUNT_DEEPENING}" = "1" ]; then
  CMD+=(--allow-targeted-leg-count-deepening)
fi

if [ "$#" -gt 0 ]; then
  for job_id in "$@"; do
    CMD+=(--job-id "${job_id}")
  done
fi

exec "${CMD[@]}"
