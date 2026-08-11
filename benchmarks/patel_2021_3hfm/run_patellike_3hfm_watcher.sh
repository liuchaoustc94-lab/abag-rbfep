#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
BENCHMARK_ROOT="${ROOT}/benchmarks/patel_2021_3hfm"
BATCH_DIR="${BATCH_DIR:-${ROOT}/runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference}"

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-13}"
MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"
MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-0}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-2}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-180}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-}"
APPEND_REST="${APPEND_REST:-0}"
SKIP_CHARGE_CHANGING="${SKIP_CHARGE_CHANGING:-1}"

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

WATCH_ARGS=(
  python3
  -u
  "${BENCHMARK_ROOT}/watch_patellike_3hfm.py"
  --batch-dir "${BATCH_DIR}"
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
)

if [ -n "${GPU_DEVICES}" ]; then
  WATCH_ARGS+=(--gpu-devices "${GPU_DEVICES}")
fi

if [ "${SKIP_CHARGE_CHANGING}" = "1" ]; then
  WATCH_ARGS+=(--skip-charge-changing)
fi

if [ -n "${POST_REFRESH_COMMAND}" ]; then
  WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")
fi

if [ "$#" -gt 0 ]; then
  if [ "${APPEND_REST}" = "1" ]; then
    WATCH_ARGS+=(--append-rest)
  else
    WATCH_ARGS+=(--only-listed)
  fi
  WATCH_ARGS+=("$@")
fi

exec "${WATCH_ARGS[@]}"
