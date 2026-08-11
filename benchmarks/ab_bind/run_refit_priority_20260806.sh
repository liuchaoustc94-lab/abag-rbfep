#!/usr/bin/env bash
# Refit runner: re-run all core_v1 jobs (fit + validation splits) at
# validation_priority preset (8 lambda / 3 repeats / 20 ps) to replace the
# 2026-06-25 pruned quick-preset data and let the fixed merge precedence
# (sampling effort > spread/stderr) upgrade the official validation view.
#
# Usage: nohup benchmarks/ab_bind/run_refit_priority_20260806.sh > logs 2>&1 &
set -uo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_refit_priority_20260806"
PROTOCOL_PATH="${BENCHMARK_ROOT}/protocol.validation_priority.yml"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
MAX_WORKERS="${MAX_WORKERS:-2}"
MAX_PASSES="${MAX_PASSES:-100}"
SLEEP_BETWEEN_PASSES="${SLEEP_BETWEEN_PASSES:-60}"

export PATH="${ROOT}/.venv/bin:${PATH}"

pass=0
while [ "${pass}" -lt "${MAX_PASSES}" ]; do
  pass=$((pass + 1))
  echo "[refit] $(date --iso-8601=seconds) pass=${pass} start"

  "${ABAG_RBFE}" batch run-abbind \
    --plan-root "${RUNS_ROOT}" \
    --resume \
    --execute \
    --max-workers "${MAX_WORKERS}" \
    --gpu-devices "${GPU_DEVICES}"
  rc=$?
  echo "[refit] $(date --iso-8601=seconds) pass=${pass} run-abbind rc=${rc}"

  "${ABAG_RBFE}" batch report-abbind --plan-root "${RUNS_ROOT}" >/dev/null 2>&1 || true

  remaining=$("${ROOT}/.venv/bin/python" - "${RUNS_ROOT}" <<'PYEOF'
import json, sys
from pathlib import Path
summary = Path(sys.argv[1]) / "reports" / "plan_summary.json"
if not summary.is_file():
    print("unknown")
    raise SystemExit(0)
payload = json.loads(summary.read_text())
jobs = payload.get("jobs") or payload.get("job_rows") or []
if not jobs:
    # fall back to run_summary.csv based counting handled by caller
    print("unknown")
    raise SystemExit(0)
incomplete = [j for j in jobs if str(j.get("latest_stage_state")) != "completed" or str(j.get("latest_stage")) != "report"]
print(len(incomplete))
PYEOF
)
  echo "[refit] $(date --iso-8601=seconds) pass=${pass} remaining=${remaining}"
  if [ "${remaining}" = "0" ]; then
    echo "[refit] all jobs completed report stage; stopping."
    break
  fi
  sleep "${SLEEP_BETWEEN_PASSES}"
done

echo "[refit] $(date --iso-8601=seconds) done after ${pass} passes"
