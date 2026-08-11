#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_calibration_rescues"
QUICK_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_quick_plan"
CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
DISABLE_EXIT_CODE="${WATCH_SUPERVISOR_NO_RESTART_CODE:-75}"

if [ "${ABAG_ENABLE_CALIBRATION_RESCUE_WATCHER:-0}" != "1" ]; then
  echo "Calibration rescue watcher is disabled by default. Use run_calibration_watchers.sh start rescue or set ABAG_ENABLE_CALIBRATION_RESCUE_WATCHER=1." >&2
  exit "${DISABLE_EXIT_CODE}"
fi

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-6}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-48}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-300}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${QUICK_RUNS_ROOT}}"
APPEND_REST="${APPEND_REST:-0}"
WATCH_ONCE="${WATCH_ONCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
REFRESH_RESCUES_EACH_PASS="${REFRESH_RESCUES_EACH_PASS:-1}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${CALIBRATED_VALIDATION_REFRESH}}"
export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=(
    "1dqj-antibody-h-y50a"
    "1dqj-antibody-h-y33a"
    "1n8z-antibody-h-w95a"
    "3be1-antibody-h-y33a"
    "1dqj-antigen-c-y20a"
    "2jel-antigen-p-s64t"
    "2nyy-antigen-a-f953a"
    "3bn9-antigen-a-f97a"
  )
fi

run_watch_pass() {
  if [ "${REFRESH_RESCUES_EACH_PASS}" = "1" ]; then
    "${ROOT}/benchmarks/ab_bind/refresh_calibration_rescues.sh" "${JOB_IDS[@]}" >/dev/null
  fi

  WATCH_ARGS=(
    python3
    -u
    "${BENCHMARK_ROOT}/watch_validation_priority.py"
    --plan-root "${RUNS_ROOT}"
    --split-name calibration
    --split-file "${SPLIT_FILE}"
    --allow-active-elsewhere-job-ids
    --poll-seconds "${POLL_SECONDS}"
    --max-compute-apps-per-gpu "${MAX_COMPUTE_APPS_PER_GPU}"
    --max-load-per-core "${MAX_LOAD_PER_CORE}"
    --max-active-mdrun-threads "${MAX_ACTIVE_MDRUN_THREADS}"
    --max-launches-per-pass "${MAX_LAUNCHES_PER_PASS}"
    --launch-cooldown-seconds "${LAUNCH_COOLDOWN_SECONDS}"
    --warn-stale-mdrun-seconds "${WARN_STALE_MDRUN_SECONDS}"
    --mdrun-args-override "${MDRUN_ARGS_OVERRIDE}"
    --merged-plan-root "${RUNS_ROOT}"
    --merged-extra-plan-root "${QUICK_RUNS_ROOT}"
    --once
  )

  if [ "${APPEND_REST}" = "1" ]; then
    WATCH_ARGS+=(--append-rest)
  else
    WATCH_ARGS+=(--only-listed)
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
