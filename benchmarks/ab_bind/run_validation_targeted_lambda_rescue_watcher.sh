#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"
PRIORITY_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"
ROBUST_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"
RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"
TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"
SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"
DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"
ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"
CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"
DISABLE_EXIT_CODE="${WATCH_SUPERVISOR_NO_RESTART_CODE:-75}"

if [ "${ABAG_ENABLE_TARGETED_LAMBDA_RESCUE_WATCHER:-0}" != "1" ]; then
  echo "Targeted lambda rescue watcher is disabled by default. Set ABAG_ENABLE_TARGETED_LAMBDA_RESCUE_WATCHER=1 to enable it." >&2
  exit "${DISABLE_EXIT_CODE}"
fi

WAIT_FOR_PID="${WAIT_FOR_PID:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
GPU_DEVICES="${GPU_DEVICES:-}"
MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-4}"
MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"
MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"
MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"
MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-24}"
MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"
MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"
LAUNCH_COOLDOWN_SECONDS="${LAUNCH_COOLDOWN_SECONDS:-180}"
WARN_STALE_MDRUN_SECONDS="${WARN_STALE_MDRUN_SECONDS:-900}"
MDRUN_ARGS_OVERRIDE="${MDRUN_ARGS_OVERRIDE:--ntmpi 1 -ntomp 2}"
ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${PRIORITY_RUNS_ROOT}}"
MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PRIORITY_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${DEEP_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${ULTRA_RESCUE_RUNS_ROOT}}"
MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"
WATCH_ONCE="${WATCH_ONCE:-0}"
DRY_RUN="${DRY_RUN:-0}"
REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"
REFRESH_TARGETED_LAMBDA_RESCUES_EACH_PASS="${REFRESH_TARGETED_LAMBDA_RESCUES_EACH_PASS:-1}"
WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_targeted_lambda_watchlist_refresh.json}"
POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"
export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS

if [ -n "${WAIT_FOR_PID}" ]; then
  while kill -0 "${WAIT_FOR_PID}" 2>/dev/null; do
    sleep 30
  done
fi

DEFAULT_JOB_IDS=(
  "1bj1-antigen-w-g88a"
  "1cz8-antigen-w-g92a"
  "1cz8-antigen-w-m81a"
  "3hfm-antibody-h-c95a"
  "3hfm-antibody-h-y33a"
  "3nps-antigen-a-h138a"
)

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
  AUTO_REFRESH_WATCHLIST=0
else
  JOB_IDS=("${DEFAULT_JOB_IDS[@]}")
  AUTO_REFRESH_WATCHLIST="${REFRESH_WATCHLIST_EACH_PASS}"
fi

append_materialized_job_ids() {
  local plan_jobs_csv="${RUNS_ROOT}/reports/plan_jobs.csv"
  if [ ! -f "${plan_jobs_csv}" ]; then
    return
  fi
  local materialized=()
  mapfile -t materialized < <(
    python3 - "${plan_jobs_csv}" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.exists():
    raise SystemExit(0)

active_states = {"running", "stale_running"}
seen = set()
with path.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        job_id = (row.get("job_id") or "").strip()
        latest_stage_state = (row.get("latest_stage_state") or "").strip()
        if not job_id or latest_stage_state not in active_states or job_id in seen:
            continue
        seen.add(job_id)
        print(job_id)
PY
  )
  if [ "${#materialized[@]}" -eq 0 ]; then
    return
  fi
  local merged=()
  local seen=()
  local job_id
  for job_id in "${JOB_IDS[@]}" "${materialized[@]}"; do
    if [ -z "${job_id}" ]; then
        continue
    fi
    if [[ " ${seen[*]} " == *" ${job_id} "* ]]; then
      continue
    fi
    seen+=("${job_id}")
    merged+=("${job_id}")
  done
  JOB_IDS=("${merged[@]}")
}

refresh_job_ids() {
  local refreshed=()
  mapfile -t refreshed < <(
    python3 -u "${BENCHMARK_ROOT}/refresh_validation_watchlists.py" \
      --root "${ROOT}" \
      --mode targeted \
      --robust-plan-root "${ROBUST_RUNS_ROOT}" \
      --rescue-plan-root "${RESCUE_RUNS_ROOT}" \
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

  if [ "${REFRESH_TARGETED_LAMBDA_RESCUES_EACH_PASS}" = "1" ]; then
    "${ROOT}/benchmarks/ab_bind/refresh_validation_targeted_lambda_rescues.sh" "${JOB_IDS[@]}" >/dev/null
  fi

  append_materialized_job_ids

  WATCH_ARGS=(
    python3
    -u
    "${BENCHMARK_ROOT}/watch_validation_priority.py"
    --plan-root "${RUNS_ROOT}"
    --split-name validation
    --split-file "${SPLIT_FILE}"
    --only-listed
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
