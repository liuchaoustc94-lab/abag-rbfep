#!/usr/bin/env bash
# Key-points rerun + protonation sensitivity pilot (2026-08-06).
#
# Root 1 (baseline, priority preset 8lambda/3rep/20ps):
#   - 1MLC x8 antibody jobs (official view was quick-plan garbage)
#   - 1BJ1/1CZ8 w-q89a + w-h90a (catastrophic outliers)
#   - 11 fit-split jobs from the 2026-06-23 calibration fit_pairs.csv
# Root 2 (protonation variants, same preset):
#   - 1MLC x8 with apo-leg H:E50 renamed GLU->GLH  (PROPKA apo pKa 7.41)
#   - 1BJ1/1CZ8 w-q89a + w-h90a with apo-leg W:H90 renamed HIS->HIP
#     (PROPKA apo pKa 6.16 vs complex 2.66/2.37)
#
# Usage: nohup bash benchmarks/ab_bind/run_keypoints_protonation_pilot_20260806.sh > log 2>&1 &
set -uo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BASELINE="${ROOT}/runs/benchmarks/abbind_keypoints_baseline_20260806"
VARIANT="${ROOT}/runs/benchmarks/abbind_protonation_pilot_20260806"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
MAX_WORKERS="${MAX_WORKERS:-2}"
export PATH="${ROOT}/.venv/bin:${PATH}"

BASELINE_JOBS=(
  1mlc-antibody-h-s57a 1mlc-antibody-h-s57v 1mlc-antibody-h-t31a 1mlc-antibody-h-t31v
  1mlc-antibody-h-t31w 1mlc-antibody-l-n32g 1mlc-antibody-l-n32y 1mlc-antibody-l-n92a
  1bj1-antigen-w-q89a 1bj1-antigen-w-h90a 1cz8-antigen-w-q89a 1cz8-antigen-w-h90a
  1dqj-antibody-h-y33a 1dqj-antibody-h-y50a 1dqj-antigen-c-y20a 1n8z-antibody-h-w95a
  2jel-antigen-p-s64t 3be1-antibody-h-y33a 1jrh-antigen-i-t14v 1vfb-antibody-h-g31a
  1vfb-antibody-h-y32a 1vfb-antigen-c-y23a 3ngb-antibody-h-g54s
)
VARIANT_JOBS=(
  1mlc-antibody-h-s57a 1mlc-antibody-h-s57v 1mlc-antibody-h-t31a 1mlc-antibody-h-t31v
  1mlc-antibody-h-t31w 1mlc-antibody-l-n32g 1mlc-antibody-l-n32y 1mlc-antibody-l-n92a
  1bj1-antigen-w-q89a 1bj1-antigen-w-h90a 1cz8-antigen-w-q89a 1cz8-antigen-w-h90a
)

run_root() {
  local plan_root="$1"; shift
  local extra_args=("$@")
  local pass=0
  while [ "${pass}" -lt 40 ]; do
    pass=$((pass + 1))
    echo "[pilot] $(date --iso-8601=seconds) root=${plan_root##*/} pass=${pass}"
    local args=(batch run-abbind --plan-root "${plan_root}" --resume --execute
                --max-workers "${MAX_WORKERS}" --gpu-devices "${GPU_DEVICES}")
    for j in "${extra_args[@]}"; do args+=(--job-id "${j}"); done
    "${ABAG_RBFE}" "${args[@]}" && break
    echo "[pilot] run-abbind rc=$?; retrying after 60s"
    sleep 60
  done
}

echo "[pilot] $(date --iso-8601=seconds) phase 1: baseline root"
run_root "${BASELINE}" "${BASELINE_JOBS[@]}"

echo "[pilot] $(date --iso-8601=seconds) phase 2: variant prepare (CPU only)"
PREP_ARGS=(batch run-abbind --plan-root "${VARIANT}" --resume --execute --to-stage prepare)
for j in "${VARIANT_JOBS[@]}"; do PREP_ARGS+=(--job-id "${j}"); done
"${ABAG_RBFE}" "${PREP_ARGS[@]}"

echo "[pilot] $(date --iso-8601=seconds) phase 3: patch apo-leg protonation variants"
"${ROOT}/.venv/bin/python" - <<'PYEOF'
from pathlib import Path

ROOT = Path("/mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_protonation_pilot_20260806")

def rename_residue(pdb: Path, chain: str, resseq: int, old: str, new: str) -> int:
    lines = pdb.read_text().splitlines(keepends=True)
    count = 0
    out = []
    for line in lines:
        if (
            line.startswith("ATOM")
            and line[21] == chain
            and int(line[22:26]) == resseq
            and line[17:20].strip() == old
        ):
            line = line[:17] + new.ljust(3) + line[20:]
            count += 1
        out.append(line)
    pdb.write_text("".join(out))
    return count

patched = []
for job_dir in sorted(ROOT.glob("*/jobs/*")):
    job_id = job_dir.name
    apo_input = job_dir / "legs" / "apo" / "input.pdb"
    if not apo_input.is_file():
        continue
    if job_id.startswith("1mlc-"):
        n = rename_residue(apo_input, "H", 50, "GLU", "GLH")
        patched.append((job_id, "H:E50->GLH", n))
    elif job_id in ("1bj1-antigen-w-q89a", "1bj1-antigen-w-h90a",
                    "1cz8-antigen-w-q89a", "1cz8-antigen-w-h90a"):
        n = rename_residue(apo_input, "W", 90, "HIS", "HIP")
        patched.append((job_id, "W:H90->HIP", n))
for job_id, what, n in patched:
    print(f"[patch] {job_id}: {what} ({n} atoms)")
PYEOF

echo "[pilot] $(date --iso-8601=seconds) phase 4: variant full run"
run_root "${VARIANT}" "${VARIANT_JOBS[@]}"

echo "[pilot] $(date --iso-8601=seconds) all done"
