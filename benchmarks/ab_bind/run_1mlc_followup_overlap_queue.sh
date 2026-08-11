#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
RUN_SCRIPT="${ROOT}/benchmarks/ab_bind/run_1mlc_overlap_focus_pilot.sh"
POLL_SECONDS="${POLL_SECONDS:-300}"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
MAX_WORKERS="${MAX_WORKERS:-2}"
QUEUE_ONCE="${QUEUE_ONCE:-0}"

PRIMARY_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_1mlc_20260626"
FOLLOWUP_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_overlap_focus_pilot_1mlc_20260626"

has_active_execution() {
  ps -ef | grep -E "${PRIMARY_ROOT}|${FOLLOWUP_ROOT}|run_1mlc_target_specific_sampling_pilot\\.sh|run_1mlc_overlap_focus_pilot\\.sh|abag-rbfe batch run-abbind" | grep -v grep | grep -q "1MLC"
}

primary_finished() {
  python3 - <<'PY' "${PRIMARY_ROOT}"
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
job_ids = [
    "1mlc-antibody-l-n92a",
    "1mlc-antibody-h-s57v",
    "1mlc-antibody-h-t31a",
]

def latest_payload(job_dir: Path):
    for rel in (
        "results/qc_report.json",
        "results/ddg_summary.json",
        "stages/bar.json",
        "stages/sample.json",
        "stages/equilibrate.json",
        "stages/build_legs.json",
        "stages/mutate.json",
        "stages/prepare.json",
    ):
        path = job_dir / rel
        if path.exists():
            return rel, json.loads(path.read_text())
    return "", {}

all_done = True
for job_id in job_ids:
    matches = list(root.glob(f"**/jobs/{job_id}"))
    if not matches:
        all_done = False
        continue
    rel, payload = latest_payload(matches[0])
    state = str(payload.get("state") or payload.get("status") or "")
    if rel in {"results/qc_report.json", "results/ddg_summary.json"}:
        continue
    if rel == "stages/bar.json" and state == "completed":
        continue
    all_done = False

print("done" if all_done else "waiting")
PY
}

followup_finished() {
  python3 - <<'PY' "${FOLLOWUP_ROOT}"
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
job_ids = [
    "1mlc-antibody-l-n92a",
    "1mlc-antibody-h-s57v",
]

done = True
for job_id in job_ids:
    matches = list(root.glob(f"**/jobs/{job_id}"))
    if not matches:
        done = False
        continue
    ok = False
    for rel in ("results/qc_report.json", "results/ddg_summary.json", "stages/bar.json"):
        path = matches[0] / rel
        if path.exists():
            payload = json.loads(path.read_text())
            state = str(payload.get("state") or payload.get("status") or "")
            if rel in {"results/qc_report.json", "results/ddg_summary.json"}:
                ok = True
                break
            if rel == "stages/bar.json" and state == "completed":
                ok = True
                break
    if not ok:
        done = False

print("done" if done else "waiting")
PY
}

main_loop() {
  while true; do
    if has_active_execution; then
      echo "[queue] 1MLC lane still active; waiting."
    elif [ "$(primary_finished)" != "done" ]; then
      echo "[queue] primary 1MLC pilot not finished yet; waiting."
    elif [ "$(followup_finished)" = "done" ]; then
      echo "[queue] overlap-focus 1MLC follow-up already finished."
      break
    else
      echo "[queue] launching overlap-focus 1MLC follow-up."
      GPU_DEVICES="${GPU_DEVICES}" MAX_WORKERS="${MAX_WORKERS}" "${RUN_SCRIPT}" || true
      break
    fi

    if [ "${QUEUE_ONCE}" = "1" ]; then
      break
    fi
    sleep "${POLL_SECONDS}"
  done
}

main_loop
