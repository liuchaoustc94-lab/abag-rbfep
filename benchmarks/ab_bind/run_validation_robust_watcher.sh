#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"
PRIORITY_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"
TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"
TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"
SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"
DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"
ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"
CALIBRATED_VALIDATION_REFRESH="${CALIBRATED_VALIDATION_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh}"
VALIDATION_POST_REPORT_REFRESH="${VALIDATION_POST_REPORT_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh}"

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-6}"
MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"
MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-56}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-2}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-180}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${PRIORITY_RUNS_ROOT}}"
MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PRIORITY_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"
WATCH_ONCE="${WATCH_ONCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"
WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_watchlist_refresh.json}"
ROBUST_PASS_OUTLIER_THRESHOLD="${ROBUST_PASS_OUTLIER_THRESHOLD:-5.0}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"
export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

DEFAULT_JOB_IDS=(
  "3hfm-antibody-h-y50a"
  "3hfm-antibody-h-y33a"
  "3hfm-antibody-h-c95a"
  "3hfm-antigen-y-y20a"
  "1cz8-antigen-w-g92a"
  "1bj1-antigen-w-g92a"
  "1mlc-antibody-h-s57a"
  "1mlc-antibody-h-s57v"
  "1mlc-antibody-h-t31a"
  "1mlc-antibody-h-t31v"
  "1mlc-antibody-l-n92a"
)

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
  AUTO_REFRESH_WATCHLIST=0
else
  JOB_IDS=("${DEFAULT_JOB_IDS[@]}")
  AUTO_REFRESH_WATCHLIST="${REFRESH_WATCHLIST_EACH_PASS}"
fi

refresh_job_ids() {
  local refreshed=()
  mapfile -t refreshed < <(
    python3 -u "${BENCHMARK_ROOT}/refresh_validation_watchlists.py" \
      --root "${ROOT}" \
      --mode robust \
      --robust-plan-root "${RUNS_ROOT}" \
      --rescue-plan-root "${RESCUE_RUNS_ROOT}" \
      --robust-pass-outlier-threshold "${ROBUST_PASS_OUTLIER_THRESHOLD}" \
      --output-json "${WATCHLIST_REFRESH_JSON}"
  )
  if [ "${#refreshed[@]}" -gt 0 ]; then
    JOB_IDS=("${refreshed[@]}")
  fi
}

run_watch_pass() {
  if [ "${AUTO_REFRESH_WATCHLIST}" = "1" ]; then
    refresh_job_ids
  fi

  WATCH_ARGS=(
    python3
    -u
    "${BENCHMARK_ROOT}/watch_validation_priority.py"
    --plan-root "${RUNS_ROOT}"
    --split-name validation
    --split-file "${SPLIT_FILE}"
    --only-listed
    --poll-seconds "${POLL_SECONDS}"
    --max-compute-apps-per-gpu "${MAX_COMPUTE_APPS_PER_GPU}"
    --min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"
    --max-gpu-utilization "${MAX_GPU_UTILIZATION}"
    --max-load-per-core "${MAX_LOAD_PER_CORE}"
    --max-active-mdrun-threads "${MAX_ACTIVE_MDRUN_THREADS}"
    --max-launches-per-pass "${MAX_LAUNCHES_PER_PASS}"
    --launch-cooldown-seconds "${LAUNCH_COOLDOWN_SECONDS}"
    --warn-stale-mdrun-seconds "${WARN_STALE_MDRUN_SECONDS}"
    --mdrun-args-override "${MDRUN_ARGS_OVERRIDE}"
    --once
  )

  if [ -n "${MERGED_PLAN_ROOT}" ]; then
    WATCH_ARGS+=(--merged-plan-root "${MERGED_PLAN_ROOT}")
  fi

  if [ -n "${MERGED_EXTRA_PLAN_ROOT_1}" ]; then
    WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_1}")
  fi

  if [ -n "${MERGED_EXTRA_PLAN_ROOT_2}" ]; then
    WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_2}")
  fi

if [ -n "${MERGED_EXTRA_PLAN_ROOT_3}" ]; then
  WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_3}")
fi

if [ -n "${MERGED_EXTRA_PLAN_ROOT_4}" ]; then
  WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_4}")
fi

if [ -n "${MERGED_EXTRA_PLAN_ROOT_5}" ]; then
  WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")
fi

if [ -n "${MERGED_EXTRA_PLAN_ROOT_6}" ]; then
  WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")
fi

if [ -n "${MERGED_EXTRA_PLAN_ROOT_7}" ]; then
  WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")
fi

  if [ -n "${GPU_DEVICES}" ]; then
    WATCH_ARGS+=(--gpu-devices "${GPU_DEVICES}")
  fi

  if [ -n "${POST_REFRESH_COMMAND}" ]; then
    WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")
  fi

  if [ "${DRY_RUN}" = "1" ]; then
    WATCH_ARGS+=(--dry-run)
  fi

  WATCH_ARGS+=("${JOB_IDS[@]}")
  "${WATCH_ARGS[@]}"
}

run_watch_pass
if [ "${WATCH_ONCE}" = "1" ]; then
  exit 0
fi

while true; do
  SLEEP_SECONDS="${POLL_SECONDS}"
  if [ "${SLEEP_SECONDS}" -lt 5 ]; then
    SLEEP_SECONDS=5
  fi
  sleep "${SLEEP_SECONDS}"
  run_watch_pass
done
