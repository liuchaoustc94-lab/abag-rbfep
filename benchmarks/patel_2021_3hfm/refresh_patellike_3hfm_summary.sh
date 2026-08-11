#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
BATCH_DIR="${BATCH_DIR:-${ROOT}/runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference}"
WATCH_DIR="${WATCH_DIR:-${BATCH_DIR}/reports/watch}"
SUMMARY_OUTPUT="${SUMMARY_OUTPUT:-${BATCH_DIR}/reports/patel_2021_3hfm_summary.json}"
PYTHON_BIN="${PYTHON_BIN:-${ROOT}/.venv/bin/python}"
REPORT_SCRIPT="${REPORT_SCRIPT:-${ROOT}/benchmarks/patel_2021_3hfm/report_patellike_3hfm.py}"
LOCK_PATH="${LOCK_PATH:-${WATCH_DIR}/patellike_3hfm_summary.lock}"
STAMP_PATH="${STAMP_PATH:-${WATCH_DIR}/patellike_3hfm_summary.last_run}"
LOG_PATH="${LOG_PATH:-${WATCH_DIR}/patellike_3hfm_summary_refresh.log}"
MIN_INTERVAL_SECONDS="${MIN_INTERVAL_SECONDS:-300}"
FAIL_ON_REFRESH_ERROR="${FAIL_ON_REFRESH_ERROR:-0}"

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
{
  printf '[refresh] %s start min_interval=%s command=%q script=%q\n' \
    "$(date --iso-8601=seconds)" \
    "${MIN_INTERVAL_SECONDS}" \
    "${PYTHON_BIN}" \
    "${REPORT_SCRIPT}"
  if "${PYTHON_BIN}" -u "${REPORT_SCRIPT}" --batch-dir "${BATCH_DIR}" --summary-output "${SUMMARY_OUTPUT}" "$@"; then
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
