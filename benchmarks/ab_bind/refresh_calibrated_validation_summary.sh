#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
PLAN_ROOT="${PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_quick_plan}"
REPORT_PLAN_ROOT="${REPORT_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_quick_plan}"
WATCH_DIR="${WATCH_DIR:-${PLAN_ROOT}/reports/watch}"
# Keep watcher locks/stamps under the caller's PLAN_ROOT, but refresh the
# canonical summary file under the report-generating root unless overridden.
SUMMARY_OUTPUT="${SUMMARY_OUTPUT:-${REPORT_PLAN_ROOT}/reports/calibrated_validation_summary.json}"
PYTHON_BIN="${PYTHON_BIN:-${ROOT}/.venv/bin/python}"
REPORT_SCRIPT="${REPORT_SCRIPT:-${ROOT}/benchmarks/ab_bind/report_calibrated_validation.py}"
LOCK_PATH="${LOCK_PATH:-${WATCH_DIR}/calibrated_validation_summary.lock}"
STAMP_PATH="${STAMP_PATH:-${WATCH_DIR}/calibrated_validation_summary.last_run}"
LOG_PATH="${LOG_PATH:-${WATCH_DIR}/calibrated_validation_summary_refresh.log}"
MIN_INTERVAL_SECONDS="${MIN_INTERVAL_SECONDS:-300}"
FAIL_ON_REFRESH_ERROR="${FAIL_ON_REFRESH_ERROR:-0}"
EXTRA_PLAN_ROOTS="${EXTRA_PLAN_ROOTS:-${MERGED_EXTRA_PLAN_ROOTS:-}}"

mkdir -p "${WATCH_DIR}"

if command -v flock >/dev/null 2>&1; then
  exec 9>"${LOCK_PATH}"
  if ! flock -n 9; then
    exit 0
  fi
fi

NOW_EPOCH="$(date +%s)"
LAST_RUN_EPOCH=""
if [ -f "${STAMP_PATH}" ]; then
  LAST_RUN_EPOCH="$(tr -d '[:space:]' < "${STAMP_PATH}" || true)"
fi

if [[ "${LAST_RUN_EPOCH}" =~ ^[0-9]+$ ]]; then
  ELAPSED_SECONDS="$((NOW_EPOCH - LAST_RUN_EPOCH))"
  if [ "${ELAPSED_SECONDS}" -lt "${MIN_INTERVAL_SECONDS}" ]; then
    exit 0
  fi
fi

status=0
extra_args=()
if [ -n "${EXTRA_PLAN_ROOTS}" ]; then
  OLD_IFS="${IFS}"
  IFS=':'
  read -r -a extra_roots <<< "${EXTRA_PLAN_ROOTS}"
  IFS="${OLD_IFS}"
  for extra_root in "${extra_roots[@]}"; do
    if [ -z "${extra_root}" ] || [ "${extra_root}" = "${REPORT_PLAN_ROOT}" ]; then
      continue
    fi
    extra_args+=(--extra-plan-root "${extra_root}")
  done
fi
{
  printf '[refresh] %s start min_interval=%s command=%q script=%q\n' \
    "$(date --iso-8601=seconds)" \
    "${MIN_INTERVAL_SECONDS}" \
    "${PYTHON_BIN}" \
    "${REPORT_SCRIPT}"
  if "${PYTHON_BIN}" -u "${REPORT_SCRIPT}" --plan-root "${REPORT_PLAN_ROOT}" --summary-output "${SUMMARY_OUTPUT}" "${extra_args[@]}" "$@"; then
    date +%s > "${STAMP_PATH}"
    printf '[refresh] %s success\n' "$(date --iso-8601=seconds)"
  else
    status=$?
    printf '[refresh] %s failed rc=%s\n' "$(date --iso-8601=seconds)" "${status}"
  fi
} >> "${LOG_PATH}" 2>&1

if [ "${status}" -ne 0 ] && [ "${FAIL_ON_REFRESH_ERROR}" = "1" ]; then
  exit "${status}"
fi

exit 0
