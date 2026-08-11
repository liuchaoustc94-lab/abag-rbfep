#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
ROBUST_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"
RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"
TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"
TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"
SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"
DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"
ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"
CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"
DISABLE_EXIT_CODE="${WATCH_SUPERVISOR_NO_RESTART_CODE:-75}"

if [ "${ABAG_ENABLE_STALE_RECOVERY_WATCHER:-0}" != "1" ]; then
  echo "Stale recovery watcher is disabled by default. Set ABAG_ENABLE_STALE_RECOVERY_WATCHER=1 to enable it." >&2
  exit "${DISABLE_EXIT_CODE}"
fi

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-14}"
MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10000}"
MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-100}"
MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-2}"
THREAD_BUDGET_PLAN_ROOT_1="${THREAD_BUDGET_PLAN_ROOT_1:-${RUNS_ROOT}}"
THREAD_BUDGET_PLAN_ROOT_2="${THREAD_BUDGET_PLAN_ROOT_2:-${ROBUST_RUNS_ROOT}}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-2}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-60}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${ROBUST_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${DEEP_RESCUE_RUNS_ROOT}:${ULTRA_RESCUE_RUNS_ROOT}}"
MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"
WATCH_ONCE="${WATCH_ONCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"
WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_stale_watchlist_refresh.json}"
INCLUDE_GAP_JOBS="${INCLUDE_GAP_JOBS:-0}"
GAP_WATCHLIST_REFRESH_JSON="${GAP_WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_gap_watchlist_refresh.json}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"
export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

DEFAULT_JOB_IDS=(
  "1bj1-antigen-v-f17a"
  "1bj1-antigen-v-y21a"
  "1bj1-antigen-w-h86a"
  "1bj1-antigen-w-h90a"
  "1bj1-antigen-w-i80a"
  "1mlc-antibody-h-s57a"
  "1mlc-antibody-h-s57v"
  "1mlc-antibody-h-t31a"
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
  local gap_refreshed=()
  local combined=()
  local seen=()
  local job_id=""
  mapfile -t refreshed < <(
    python3 -u "${BENCHMARK_ROOT}/refresh_validation_watchlists.py" \
      --root "${ROOT}" \
      --mode stale \
      --priority-plan-root "${RUNS_ROOT}" \
      --candidate-plan-root "${ROBUST_RUNS_ROOT}" \
      --candidate-plan-root "${RESCUE_RUNS_ROOT}" \
      --candidate-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}" \
      --candidate-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}" \
      --candidate-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}" \
      --candidate-plan-root "${DEEP_RESCUE_RUNS_ROOT}" \
      --candidate-plan-root "${ULTRA_RESCUE_RUNS_ROOT}" \
      --output-json "${WATCHLIST_REFRESH_JSON}"
  )
  if [ "${INCLUDE_GAP_JOBS}" = "1" ]; then
    mapfile -t gap_refreshed < <(
      python3 -u "${BENCHMARK_ROOT}/refresh_validation_watchlists.py" \
        --root "${ROOT}" \
        --mode gap \
        --priority-plan-root "${RUNS_ROOT}" \
        --candidate-plan-root "${ROBUST_RUNS_ROOT}" \
        --candidate-plan-root "${RESCUE_RUNS_ROOT}" \
        --candidate-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}" \
        --candidate-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}" \
        --candidate-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}" \
        --candidate-plan-root "${DEEP_RESCUE_RUNS_ROOT}" \
        --candidate-plan-root "${ULTRA_RESCUE_RUNS_ROOT}" \
        --output-json "${GAP_WATCHLIST_REFRESH_JSON}"
    )
  fi
  combined=()
  seen=()
  for job_id in "${refreshed[@]}" "${gap_refreshed[@]}"; do
    if [ -z "${job_id}" ]; then
      continue
    fi
    if [[ " ${seen[*]} " == *" ${job_id} "* ]]; then
      continue
    fi
    combined+=("${job_id}")
    seen+=("${job_id}")
  done
  if [ "${#combined[@]}" -gt 0 ]; then
    JOB_IDS=("${combined[@]}")
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
    --watch-tag stale
    --allow-active-elsewhere-job-ids
    --poll-seconds "${POLL_SECONDS}"
    --max-compute-apps-per-gpu "${MAX_COMPUTE_APPS_PER_GPU}"
    --min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"
    --max-gpu-utilization "${MAX_GPU_UTILIZATION}"
    --max-load-per-core "${MAX_LOAD_PER_CORE}"
    --max-active-mdrun-threads "${MAX_ACTIVE_MDRUN_THREADS}"
    --max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"
    --max-launches-per-pass "${MAX_LAUNCHES_PER_PASS}"
    --launch-cooldown-seconds "${LAUNCH_COOLDOWN_SECONDS}"
    --warn-stale-mdrun-seconds "${WARN_STALE_MDRUN_SECONDS}"
    --mdrun-args-override "${MDRUN_ARGS_OVERRIDE}"
    --once
  )

  if [ -n "${THREAD_BUDGET_PLAN_ROOT_1}" ]; then
    WATCH_ARGS+=(--thread-budget-plan-root "${THREAD_BUDGET_PLAN_ROOT_1}")
  fi

  if [ -n "${THREAD_BUDGET_PLAN_ROOT_2}" ]; then
    WATCH_ARGS+=(--thread-budget-plan-root "${THREAD_BUDGET_PLAN_ROOT_2}")
  fi

  WATCH_ARGS+=(--extra-plan-root "${ROBUST_RUNS_ROOT}")
  WATCH_ARGS+=(--extra-plan-root "${RESCUE_RUNS_ROOT}")
  WATCH_ARGS+=(--extra-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}")
  WATCH_ARGS+=(--extra-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}")
  WATCH_ARGS+=(--extra-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}")
  WATCH_ARGS+=(--extra-plan-root "${DEEP_RESCUE_RUNS_ROOT}")
  WATCH_ARGS+=(--extra-plan-root "${ULTRA_RESCUE_RUNS_ROOT}")

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
