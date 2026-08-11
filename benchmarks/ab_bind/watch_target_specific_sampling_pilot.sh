#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
RUN_SCRIPT="${ROOT}/benchmarks/ab_bind/run_target_specific_sampling_pilot.sh"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SUMMARY_PATH="${ROOT}/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_20260625/reports/3hfm_protocol_regression_summary.json"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_20260625"
POLL_SECONDS="${POLL_SECONDS:-300}"
WATCH_ONCE="${WATCH_ONCE:-0}"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
MAX_WORKERS="${MAX_WORKERS:-2}"

DEFAULT_JOB_IDS=(
  "3hfm-antibody-h-y33a"
  "3hfm-antibody-h-y50a"
  "3hfm-antibody-h-c95a"
  "3hfm-antigen-y-y20a"
)

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=("${DEFAULT_JOB_IDS[@]}")
fi

refresh_reports() {
  "${ROOT}/.venv/bin/python" "${BENCHMARK_ROOT}/report_3hfm_protocol_regression.py" \
    --plan-root "${RUNS_ROOT}" \
    --complex-id 3HFM \
    --summary-output "${SUMMARY_PATH}" \
    >/dev/null 2>&1 || true

  "${ROOT}/.venv/bin/python" "${BENCHMARK_ROOT}/report_hotspot_root_comparison.py" \
    --json-output "${RUNS_ROOT}/reports/hotspot_root_comparison.json" \
    --md-output "${RUNS_ROOT}/reports/hotspot_root_comparison.md" \
    >/dev/null 2>&1 || true
}

has_active_pilot_execution() {
  ps -ef | grep -F "${RUNS_ROOT}" | grep -E 'abag-rbfe batch run-abbind|gmx mdrun' | grep -v grep >/dev/null 2>&1
}

run_pass() {
  if has_active_pilot_execution; then
    echo "Active pilot execution detected under ${RUNS_ROOT}; refreshing reports only."
    refresh_reports
  else
    GPU_DEVICES="${GPU_DEVICES}" MAX_WORKERS="${MAX_WORKERS}" "${RUN_SCRIPT}" "${JOB_IDS[@]}" || true
  fi

  if [ -f "${SUMMARY_PATH}" ]; then
    python3 - <<'PY' "${SUMMARY_PATH}"
import json, sys
path = sys.argv[1]
data = json.load(open(path, "r", encoding="utf-8"))
print(
    json.dumps(
        {
            "generated_at": data.get("generated_at"),
            "status": data.get("status"),
            "selected_job_count": data.get("selected_job_count"),
            "ddg_ready_count": data.get("ddg_ready_count"),
            "paired_job_count": data.get("paired_job_count"),
            "resumable_job_count": data.get("resumable_job_count"),
            "running_equilibrate_job_count": data.get("running_equilibrate_job_count"),
            "running_sample_job_count": data.get("running_sample_job_count"),
        },
        ensure_ascii=False,
    )
)
PY
  fi
}

run_pass
if [ "${WATCH_ONCE}" = "1" ]; then
  exit 0
fi

while true; do
  sleep "${POLL_SECONDS}"
  run_pass
done
