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

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-12}"
MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10500}"
MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-96}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-2}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-60}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
WATCH_ONCE="${WATCH_ONCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"
WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_backlog_watchlist_refresh.json}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=(
    "2nz9-antigen-a-f953a"
    "2nz9-antigen-a-t1063a"
    "2nz9-antigen-a-l919a"
    "2nz9-antigen-a-n918a"
    "3hfm-antibody-l-y50a"
    "3hfm-antibody-l-y50l"
    "3hfm-antibody-h-c95f"
    "1cz8-antigen-w-g88a"
  )
  AUTO_REFRESH_WATCHLIST="${REFRESH_WATCHLIST_EACH_PASS}"
fi

if [ "$#" -gt 0 ]; then
  AUTO_REFRESH_WATCHLIST=0
fi

refresh_job_ids() {
  local refreshed=()
  mapfile -t refreshed < <(
    python3 -u "${BENCHMARK_ROOT}/refresh_validation_watchlists.py" \
      --root "${ROOT}" \
      --mode backlog \
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
    --watch-tag backlog
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
    --merged-plan-root "${RUNS_ROOT}"
    --merged-extra-plan-root "${ROBUST_RUNS_ROOT}"
    --merged-extra-plan-root "${RESCUE_RUNS_ROOT}"
    --merged-extra-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}"
    --merged-extra-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}"
    --merged-extra-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}"
    --merged-extra-plan-root "${DEEP_RESCUE_RUNS_ROOT}"
    --merged-extra-plan-root "${ULTRA_RESCUE_RUNS_ROOT}"
  )

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
