#!/usr/bin/env bash
# 1DVF new-target validation runner (2026-08-06).
# Waits for the keypoints/protonation pilot to release GPUs, then runs:
#   - 17 baseline jobs (priority preset) in runs/real_cases/1dvf_priority_20260806
#   - 1 protonation variant: antigen-d-h33a with COMPLEX-leg D:H33 HIS->HIP
#     (PROPKA: complex pKa 7.87 vs apo 5.09; charged in complex, neutral in apo)
# Usage: nohup bash benchmarks/ab_bind/run_1dvf_newtarget_20260806.sh > log 2>&1 &
set -uo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BASELINE="${ROOT}/runs/real_cases/1dvf_priority_20260806"
VARIANT="${ROOT}/runs/real_cases/1dvf_protonation_20260806"
PILOT_PATTERN="run_keypoints_protonation_pilot_20260806"
export PATH="${ROOT}/.venv/bin:${PATH}"

echo "[1dvf] $(date --iso-8601=seconds) waiting for pilot to finish..."
while pgrep -f "${PILOT_PATTERN}" > /dev/null; do
  sleep 120
done
echo "[1dvf] $(date --iso-8601=seconds) pilot done; starting 1DVF baseline"

JOBS=(
  1dvf-d13-e52-antibody-a-h30a 1dvf-d13-e52-antibody-a-s93a 1dvf-d13-e52-antibody-a-w92a
  1dvf-d13-e52-antibody-a-y32a 1dvf-d13-e52-antibody-a-y49a 1dvf-d13-e52-antibody-a-y50a
  1dvf-d13-e52-antibody-b-n56a 1dvf-d13-e52-antibody-b-t30a 1dvf-d13-e52-antibody-b-w52a
  1dvf-d13-e52-antibody-b-y101f 1dvf-d13-e52-antibody-b-y32a 1dvf-d13-e52-antigen-c-y49a
  1dvf-d13-e52-antigen-d-h33a 1dvf-d13-e52-antigen-d-i101a 1dvf-d13-e52-antigen-d-n55a
  1dvf-d13-e52-antigen-d-q104a 1dvf-d13-e52-antigen-d-y102a
)

run_queue() {
  local gpu="$1"; shift
  for job in "$@"; do
    echo "[1dvf][gpu${gpu}] $(date --iso-8601=seconds) run ${job}"
    ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" \
      --batch-dir "${BASELINE}" --execute 2>&1 | tail -1
  done
}

queue0=(); queue1=()
for i in "${!JOBS[@]}"; do
  if [ $((i % 2)) -eq 0 ]; then queue0+=("${JOBS[$i]}"); else queue1+=("${JOBS[$i]}"); fi
done

run_queue 0 "${queue0[@]}" &
Q0=$!
run_queue 1 "${queue1[@]}" &
Q1=$!
wait "${Q0}" "${Q1}"

echo "[1dvf] $(date --iso-8601=seconds) baseline done; preparing H33HIP variant"
"${ABAG_RBFE}" run "1dvf-d13-e52-antigen-d-h33a" \
  --batch-dir "${VARIANT}" --execute --to-stage prepare 2>&1 | tail -1

"${ROOT}/.venv/bin/python" - <<'PYEOF'
from pathlib import Path
pdb = Path("/mnt/data/liuchao/abag-rbfep/runs/real_cases/1dvf_protonation_20260806/jobs/1dvf-d13-e52-antigen-d-h33a/legs/complex/input.pdb")
lines = pdb.read_text().splitlines(keepends=True)
out, n = [], 0
for line in lines:
    if line.startswith("ATOM") and line[21] == "D" and int(line[22:26]) == 33 and line[17:20].strip() == "HIS":
        line = line[:17] + "HIP" + line[20:]
        n += 1
    out.append(line)
pdb.write_text("".join(out))
print(f"[patch] complex-leg D:H33 HIS->HIP ({n} atoms)")
PYEOF

echo "[1dvf] $(date --iso-8601=seconds) running H33HIP variant on gpu0"
ABAG_RBFE_VISIBLE_GPUS=0 "${ABAG_RBFE}" resume "1dvf-d13-e52-antigen-d-h33a" \
  --batch-dir "${VARIANT}" --execute 2>&1 | tail -1

echo "[1dvf] $(date --iso-8601=seconds) all done"
