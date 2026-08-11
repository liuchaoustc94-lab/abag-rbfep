#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
PYTHON_BIN="${PYTHON_BIN:-${ROOT}/.venv/bin/python}"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_quick_plan"
PROTOCOL_PATH="${BENCHMARK_ROOT}/protocol.quick.yml"
CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
REFRESH_CALIBRATION_WATCHLIST="${ROOT}/benchmarks/ab_bind/refresh_calibration_watchlist.py"
DISABLE_EXIT_CODE="${WATCH_SUPERVISOR_NO_RESTART_CODE:-75}"

if [ "${ABAG_ENABLE_QUICK_WATCHER:-0}" != "1" ]; then
  echo "Calibration quick watcher is disabled by default. Use run_calibration_watchers.sh start or set ABAG_ENABLE_QUICK_WATCHER=1." >&2
  exit "${DISABLE_EXIT_CODE}"
fi

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-6}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-64}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-300}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
APPEND_REST="${APPEND_REST:-0}"
WATCH_ONCE="${WATCH_ONCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"
WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/calibration_watchlist_refresh.json}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${CALIBRATED_VALIDATION_REFRESH}}"

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

DEFAULT_JOB_IDS=(
  "1dqj-antibody-h-y50a"
  "1dqj-antibody-h-y33a"
  "1n8z-antibody-h-w95a"
  "3be1-antibody-h-y33a"
  "1dqj-antigen-c-y20a"
  "2jel-antigen-p-s64t"
  "2nyy-antigen-a-f953a"
  "3bn9-antigen-a-f97a"
)

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
  AUTO_REFRESH_WATCHLIST=0
else
  JOB_IDS=("${DEFAULT_JOB_IDS[@]}")
  AUTO_REFRESH_WATCHLIST="${REFRESH_WATCHLIST_EACH_PASS}"
fi

"${ROOT}/benchmarks/ab_bind/refresh_curated.sh" >/dev/null
"${ROOT}/benchmarks/ab_bind/materialize_inputs.sh" >/dev/null

"${ABAG_RBFE}" batch plan-abbind \
  --benchmark-root "${BENCHMARK_ROOT}" \
  --protocol "${PROTOCOL_PATH}" \
  --spec core_v1 \
  --runs-root "${RUNS_ROOT}"

refresh_job_ids() {
  local refreshed=()
  mapfile -t refreshed < <(
    "${PYTHON_BIN}" -u "${REFRESH_CALIBRATION_WATCHLIST}" \
      --root "${ROOT}" \
      --plan-root "${RUNS_ROOT}" \
      --split-file "${SPLIT_FILE}" \
      --output-json "${WATCHLIST_REFRESH_JSON}"
  )
  if [ "${#refreshed[@]}" -gt 0 ]; then
    JOB_IDS=("${refreshed[@]}")
  fi
}

if [ "${AUTO_REFRESH_WATCHLIST}" = "1" ]; then
  refresh_job_ids
fi

run_watch_pass() {
  if [ "${AUTO_REFRESH_WATCHLIST}" = "1" ]; then
    refresh_job_ids
  fi

  WATCH_ARGS=(
    "${PYTHON_BIN}"
    -u
    "${BENCHMARK_ROOT}/watch_validation_priority.py"
    --plan-root "${RUNS_ROOT}"
    --split-name calibration
    --split-file "${SPLIT_FILE}"
    --poll-seconds "${POLL_SECONDS}"
    --max-compute-apps-per-gpu "${MAX_COMPUTE_APPS_PER_GPU}"
    --max-load-per-core "${MAX_LOAD_PER_CORE}"
    --max-active-mdrun-threads "${MAX_ACTIVE_MDRUN_THREADS}"
    --max-launches-per-pass "${MAX_LAUNCHES_PER_PASS}"
    --launch-cooldown-seconds "${LAUNCH_COOLDOWN_SECONDS}"
    --warn-stale-mdrun-seconds "${WARN_STALE_MDRUN_SECONDS}"
    --mdrun-args-override "${MDRUN_ARGS_OVERRIDE}"
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
