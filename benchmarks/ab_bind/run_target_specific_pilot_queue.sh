#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
THREE_HFM_SCRIPT="${ROOT}/benchmarks/ab_bind/run_target_specific_sampling_pilot.sh"
ONE_MLC_SCRIPT="${ROOT}/benchmarks/ab_bind/run_1mlc_target_specific_sampling_minibatch.sh"
POLL_SECONDS="${POLL_SECONDS:-300}"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
MAX_WORKERS="${MAX_WORKERS:-2}"
QUEUE_ONCE="${QUEUE_ONCE:-0}"

THREE_HFM_PATTERN="abbind_core_v1_validation_target_specific_sampling_pilot_20260625"
THREE_HFM_JOB_IDS=(
  "3hfm-antibody-h-y33a"
  "3hfm-antibody-h-y50a"
  "3hfm-antibody-h-c95a"
  "3hfm-antigen-y-y20a"
)

_active_three_hfm_count() {
  python3 - <<'PY' "${THREE_HFM_PATTERN}"
import subprocess
import sys

pattern = sys.argv[1]
proc = subprocess.run(["ps", "-ef"], capture_output=True, text=True, check=True)
count = 0
for line in proc.stdout.splitlines():
    if pattern not in line:
        continue
    if "run_target_specific_pilot_queue.sh" in line:
        continue
    if "python - <<'PY'" in line or "python3 - <<'PY'" in line:
        continue
    if "grep -E" in line or "rg " in line:
        continue
    if (
        " batch run-abbind " not in line
        and "/artifacts/commands/sample.sh" not in line
        and "gmx mdrun" not in line
    ):
        continue
    count += 1
print(count)
PY
}

_run_three_hfm_pass() {
  GPU_DEVICES="${GPU_DEVICES}" MAX_WORKERS="${MAX_WORKERS}" "${THREE_HFM_SCRIPT}" || true
}

_run_one_mlc_pass() {
  GPU_DEVICES="${GPU_DEVICES}" MAX_WORKERS="${MAX_WORKERS}" "${ONE_MLC_SCRIPT}" || true
}

_queue_log() {
  printf '[queue][%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

_three_hfm_next_action() {
  python3 - <<'PY' "${ROOT}" "${THREE_HFM_JOB_IDS[@]}"
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
job_ids = sys.argv[2:]
plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_target_specific_sampling_pilot_20260625"

def latest_state(job_dir: Path) -> tuple[str, str]:
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
            payload = json.loads(path.read_text())
            return rel, str(payload.get("state") or payload.get("status") or "")
    return "", ""

needs_more = False
for job_id in job_ids:
    matches = list(plan_root.glob(f"**/jobs/{job_id}"))
    if not matches:
        needs_more = True
        continue
    rel, state = latest_state(matches[0])
    if rel == "results/qc_report.json":
        continue
    if rel == "results/ddg_summary.json":
        ddg_payload = json.loads((matches[0] / rel).read_text())
        if ddg_payload.get("ddg_ready"):
            continue
    needs_more = True

print("rerun_3hfm" if needs_more else "launch_1mlc")
PY
}

main_loop() {
  while true; do
    active_count="$(_active_three_hfm_count)"
    if [ "${active_count}" -gt 0 ]; then
      _queue_log "3HFM lane still active (${active_count} process matches); waiting."
    else
      next_action="$(_three_hfm_next_action)"
      if [ "${next_action}" = "rerun_3hfm" ]; then
        _queue_log "3HFM lane idle but unfinished; resuming 3HFM pilot."
        _run_three_hfm_pass
      else
        _queue_log "3HFM hotspot set finished; launching 1MLC pilot."
        _run_one_mlc_pass
        break
      fi
    fi

    if [ "${QUEUE_ONCE}" = "1" ]; then
      break
    fi
    sleep "${POLL_SECONDS}"
  done
}

main_loop
