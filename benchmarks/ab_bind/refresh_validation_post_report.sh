#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
PLAN_ROOT="${PLAN_ROOT:-${MERGED_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan}}"
MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PLAN_ROOT}}"
WATCH_DIR="${WATCH_DIR:-${PLAN_ROOT}/reports/watch}"
LOG_PATH="${LOG_PATH:-${WATCH_DIR}/validation_post_report_refresh.log}"
LOCK_PATH="${LOCK_PATH:-${WATCH_DIR}/validation_post_report_refresh.lock}"
PYTHON_BIN="${PYTHON_BIN:-${ROOT}/.venv/bin/python}"
CALIBRATED_VALIDATION_REFRESH="${CALIBRATED_VALIDATION_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh}"
THREE_HFM_REPORT_SCRIPT="${THREE_HFM_REPORT_SCRIPT:-${ROOT}/benchmarks/ab_bind/report_3hfm_protocol_regression.py}"
PATELLIKE_3HFM_REFRESH="${PATELLIKE_3HFM_REFRESH:-${ROOT}/benchmarks/patel_2021_3hfm/refresh_patellike_3hfm_summary.sh}"
PATELLIKE_3HFM_BATCH_DIR="${PATELLIKE_3HFM_BATCH_DIR:-${ROOT}/runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference}"
PATELLIKE_3HFM_SUMMARY_OUTPUT="${PATELLIKE_3HFM_SUMMARY_OUTPUT:-${PATELLIKE_3HFM_BATCH_DIR}/reports/patel_2021_3hfm_summary.json}"
PATELLIKE_3HFM_MIN_INTERVAL_SECONDS="${PATELLIKE_3HFM_MIN_INTERVAL_SECONDS:-0}"
PATELLIKE_3HFM_FAIL_ON_REFRESH_ERROR="${PATELLIKE_3HFM_FAIL_ON_REFRESH_ERROR:-1}"
VALIDATION_STATUS_REPORT_SCRIPT="${VALIDATION_STATUS_REPORT_SCRIPT:-${ROOT}/benchmarks/ab_bind/report_validation_status.py}"
VALIDATION_STATUS_OUTPUT="${VALIDATION_STATUS_OUTPUT:-${ROOT}/docs/validation_status.md}"
PROJECT_COMPLETION_REPORT_SCRIPT="${PROJECT_COMPLETION_REPORT_SCRIPT:-${ROOT}/benchmarks/ab_bind/report_project_completion.py}"
PROJECT_COMPLETION_OUTPUT="${PROJECT_COMPLETION_OUTPUT:-${ROOT}/docs/project_completion_status.md}"
PROJECT_COMPLETION_JSON_OUTPUT="${PROJECT_COMPLETION_JSON_OUTPUT:-${ROOT}/runs/benchmarks/project_completion_summary.json}"
THREE_HFM_COMPLEX_ID="${THREE_HFM_COMPLEX_ID:-3HFM}"
THREE_HFM_PLAN_ROOT="${THREE_HFM_PLAN_ROOT:-${MERGED_PLAN_ROOT}}"
THREE_HFM_EXTRA_PLAN_ROOTS="${THREE_HFM_EXTRA_PLAN_ROOTS:-${MERGED_EXTRA_PLAN_ROOTS:-}}"
FAIL_ON_REFRESH_ERROR="${FAIL_ON_REFRESH_ERROR:-0}"

mkdir -p "${WATCH_DIR}"

if command -v flock >/dev/null 2>&1; then
  exec 9>"${LOCK_PATH}"
  if ! flock -n 9; then
    exit 0
  fi
fi

status=0
{
  printf '[refresh] %s start plan_root=%q merged_plan_root=%q complex_id=%q\n' \
    "$(date --iso-8601=seconds)" \
    "${PLAN_ROOT}" \
    "${MERGED_PLAN_ROOT}" \
    "${THREE_HFM_COMPLEX_ID}"

  if "${CALIBRATED_VALIDATION_REFRESH}" "$@"; then
    printf '[refresh] %s calibrated_validation success\n' "$(date --iso-8601=seconds)"
  else
    status=$?
    printf '[refresh] %s calibrated_validation failed rc=%s\n' "$(date --iso-8601=seconds)" "${status}"
  fi

  extra_args=()
  if [ -n "${THREE_HFM_EXTRA_PLAN_ROOTS}" ]; then
    OLD_IFS="${IFS}"
    IFS=':'
    read -r -a extra_roots <<< "${THREE_HFM_EXTRA_PLAN_ROOTS}"
    IFS="${OLD_IFS}"
    for extra_root in "${extra_roots[@]}"; do
      if [ -z "${extra_root}" ] || [ "${extra_root}" = "${THREE_HFM_PLAN_ROOT}" ]; then
        continue
      fi
      extra_args+=(--extra-plan-root "${extra_root}")
    done
  fi

  regression_status=0
  if "${PYTHON_BIN}" -u "${THREE_HFM_REPORT_SCRIPT}" \
    --plan-root "${THREE_HFM_PLAN_ROOT}" \
    --complex-id "${THREE_HFM_COMPLEX_ID}" \
    "${extra_args[@]}"; then
    printf '[refresh] %s 3hfm_protocol_regression success\n' "$(date --iso-8601=seconds)"
  else
    regression_status=$?
    printf '[refresh] %s 3hfm_protocol_regression failed rc=%s\n' \
      "$(date --iso-8601=seconds)" \
      "${regression_status}"
    if [ "${regression_status}" -ne 2 ] && [ "${status}" -eq 0 ]; then
      status="${regression_status}"
    fi
  fi

  patellike_3hfm_refresh_status=0
  if BATCH_DIR="${PATELLIKE_3HFM_BATCH_DIR}" \
    SUMMARY_OUTPUT="${PATELLIKE_3HFM_SUMMARY_OUTPUT}" \
    MIN_INTERVAL_SECONDS="${PATELLIKE_3HFM_MIN_INTERVAL_SECONDS}" \
    FAIL_ON_REFRESH_ERROR="${PATELLIKE_3HFM_FAIL_ON_REFRESH_ERROR}" \
    "${PATELLIKE_3HFM_REFRESH}"; then
    printf '[refresh] %s patellike_3hfm success\n' "$(date --iso-8601=seconds)"
  else
    patellike_3hfm_refresh_status=$?
    printf '[refresh] %s patellike_3hfm failed rc=%s\n' \
      "$(date --iso-8601=seconds)" \
      "${patellike_3hfm_refresh_status}"
    if [ "${patellike_3hfm_refresh_status}" -ne 2 ] && [ "${status}" -eq 0 ]; then
      status="${patellike_3hfm_refresh_status}"
    fi
  fi

  validation_status_refresh_status=0
  if "${PYTHON_BIN}" -u "${VALIDATION_STATUS_REPORT_SCRIPT}" \
    --root "${ROOT}" \
    --summary-output "${VALIDATION_STATUS_OUTPUT}"; then
    printf '[refresh] %s validation_status success\n' "$(date --iso-8601=seconds)"
  else
    validation_status_refresh_status=$?
    printf '[refresh] %s validation_status failed rc=%s\n' \
      "$(date --iso-8601=seconds)" \
      "${validation_status_refresh_status}"
    if [ "${status}" -eq 0 ]; then
      status="${validation_status_refresh_status}"
    fi
  fi

  project_completion_refresh_status=0
  if "${PYTHON_BIN}" -u "${PROJECT_COMPLETION_REPORT_SCRIPT}" \
    --root "${ROOT}" \
    --summary-output "${PROJECT_COMPLETION_OUTPUT}" \
    --json-output "${PROJECT_COMPLETION_JSON_OUTPUT}"; then
    printf '[refresh] %s project_completion success\n' "$(date --iso-8601=seconds)"
  else
    project_completion_refresh_status=$?
    printf '[refresh] %s project_completion failed rc=%s\n' \
      "$(date --iso-8601=seconds)" \
      "${project_completion_refresh_status}"
    if [ "${status}" -eq 0 ]; then
      status="${project_completion_refresh_status}"
    fi
  fi

  printf '[refresh] %s done status=%s\n' "$(date --iso-8601=seconds)" "${status}"
} >> "${LOG_PATH}" 2>&1

if [ "${status}" -ne 0 ] && [ "${FAIL_ON_REFRESH_ERROR}" = "1" ]; then
  exit "${status}"
fi

exit 0
