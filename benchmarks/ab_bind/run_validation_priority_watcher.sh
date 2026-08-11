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

if [ "${ABAG_ENABLE_PRIORITY_WATCHER:-0}" != "1" ]; then
  echo "Priority watcher is disabled by default. Use run_validation_watchers.sh start priority|all or set ABAG_ENABLE_PRIORITY_WATCHER=1." >&2
  exit "${DISABLE_EXIT_CODE}"
fi

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-5}"
MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"
MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-80}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-4}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-60}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOTS="${MERGED_EXTRA_PLAN_ROOTS:-}"
APPEND_REST="${APPEND_REST:-1}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

WATCH_ARGS=(
  python3
  -u
  "${BENCHMARK_ROOT}/watch_validation_priority.py"
  --plan-root "${RUNS_ROOT}"
  --split-name validation
  --split-file "${SPLIT_FILE}"
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

if [ -n "${MERGED_EXTRA_PLAN_ROOTS}" ]; then
  IFS=':' read -r -a MERGED_EXTRA_PLAN_ROOT_VALUES <<< "${MERGED_EXTRA_PLAN_ROOTS}"
  for merged_extra_plan_root in "${MERGED_EXTRA_PLAN_ROOT_VALUES[@]}"; do
    merged_extra_plan_root="${merged_extra_plan_root//[[:space:]]/}"
    if [ -n "${merged_extra_plan_root}" ]; then
      WATCH_ARGS+=(--merged-extra-plan-root "${merged_extra_plan_root}")
    fi
  done
fi

if [ -n "${GPU_DEVICES}" ]; then
  WATCH_ARGS+=(--gpu-devices "${GPU_DEVICES}")
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
