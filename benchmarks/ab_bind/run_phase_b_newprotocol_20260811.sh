#!/usr/bin/env bash
# Phase B: new-protocol (adaptive lambda + decoupled schedule + dt=2fs) validation.
# B3: 11 fit-pair jobs (fresh, unified protocol)
# B1: 1DVF 6 worst Y/W->A jobs (adaptive 16-lambda test)
# B2: 1DVF H33 HIP variant (decoupled-schedule overlap test)
# Usage: nohup bash benchmarks/ab_bind/run_phase_b_newprotocol_20260811.sh > log 2>&1 &
set -uo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
PROTOCOL="${BENCHMARK_ROOT}/protocol.validation_priority.yml"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
MAX_WORKERS="${MAX_WORKERS:-2}"
export PATH="${ROOT}/.venv/bin:${PATH}"

FIT_ROOT="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811"
B1_BATCH="${ROOT}/runs/real_cases/1dvf_adaptive_20260811"
B2_BATCH="${ROOT}/runs/real_cases/1dvf_h33hip_newprotocol_20260811"

FIT_JOBS=(
  1dqj-antibody-h-y33a 1dqj-antibody-h-y50a 1dqj-antigen-c-y20a 1n8z-antibody-h-w95a
  2jel-antigen-p-s64t 3be1-antibody-h-y33a 1jrh-antigen-i-t14v 1vfb-antibody-h-g31a
  1vfb-antibody-h-y32a 1vfb-antigen-c-y23a 3ngb-antibody-h-g54s
)
B1_JOBS=(
  1dvf-d13-e52-antigen-d-y102a 1dvf-d13-e52-antibody-b-w52a 1dvf-d13-e52-antibody-a-y49a
  1dvf-d13-e52-antigen-c-y49a 1dvf-d13-e52-antigen-d-q104a 1dvf-d13-e52-antibody-a-y32a
)

# ---------- B3: fit pairs ----------
echo "[B3] $(date --iso-8601=seconds) plan fit root"
COMPLEX_ARGS=()
for c in 1DQJ 1N8Z 2JEL 3BE1 1JRH 1VFB 3NGB; do COMPLEX_ARGS+=(--complex-id "${c}"); done
"${ABAG_RBFE}" batch plan-abbind \
  --benchmark-root "${BENCHMARK_ROOT}" --protocol "${PROTOCOL}" --spec core_v1 \
  --runs-root "${FIT_ROOT}" "${COMPLEX_ARGS[@]}"

FIT_ARGS=(batch run-abbind --plan-root "${FIT_ROOT}" --resume --execute
          --max-workers "${MAX_WORKERS}" --gpu-devices "${GPU_DEVICES}")
for j in "${FIT_JOBS[@]}"; do FIT_ARGS+=(--job-id "${j}"); done
echo "[B3] $(date --iso-8601=seconds) run 11 fit jobs"
"${ABAG_RBFE}" "${FIT_ARGS[@]}"

# ---------- B1: 1DVF adaptive lambda ----------
echo "[B1] $(date --iso-8601=seconds) plan 1dvf adaptive batch"
cat > /tmp/1dvf_b1_mutations.csv <<'CSV'
mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side
antigen_d_y102a,D,102,,Y,A,antigen
antibody_b_w52a,B,52,,W,A,antibody
antibody_a_y49a,A,49,,Y,A,antibody
antigen_c_y49a,C,49,,Y,A,antigen
antigen_d_q104a,D,104,,Q,A,antigen
antibody_a_y32a,A,32,,Y,A,antibody
CSV
"${ABAG_RBFE}" batch plan \
  --system "${ROOT}/examples/real_cases/1dvf/system.yml" \
  --mutations /tmp/1dvf_b1_mutations.csv \
  --protocol "${PROTOCOL}" \
  --batch-id "$(basename ${B1_BATCH})" \
  --runs-root "${ROOT}/runs/real_cases"

run_queue() {
  local batch="$1" gpu="$2"; shift 2
  for job in "$@"; do
    echo "[B1][gpu${gpu}] $(date --iso-8601=seconds) run ${job}"
    ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" --batch-dir "${batch}" --execute 2>&1 | tail -1
  done
}
run_queue "${B1_BATCH}" 0 "${B1_JOBS[0]}" "${B1_JOBS[2]}" "${B1_JOBS[4]}" &
Q0=$!
run_queue "${B1_BATCH}" 1 "${B1_JOBS[1]}" "${B1_JOBS[3]}" "${B1_JOBS[5]}" &
Q1=$!
wait "${Q0}" "${Q1}"

# ---------- B2: H33 HIP with decoupled schedule ----------
echo "[B2] $(date --iso-8601=seconds) plan h33hip batch"
cat > /tmp/1dvf_b2_mutations.csv <<'CSV'
mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side
antigen_d_h33a,D,33,,H,A,antigen
CSV
"${ABAG_RBFE}" batch plan \
  --system "${ROOT}/examples/real_cases/1dvf/system.yml" \
  --mutations /tmp/1dvf_b2_mutations.csv \
  --protocol "${PROTOCOL}" \
  --batch-id "$(basename ${B2_BATCH})" \
  --runs-root "${ROOT}/runs/real_cases"

"${ABAG_RBFE}" run 1dvf-d13-e52-antigen-d-h33a --batch-dir "${B2_BATCH}" --execute --to-stage prepare 2>&1 | tail -1
"${ROOT}/.venv/bin/python" - <<'PYEOF'
from pathlib import Path
pdb = Path("/mnt/data/liuchao/abag-rbfep/runs/real_cases/1dvf_h33hip_newprotocol_20260811/jobs/1dvf-d13-e52-antigen-d-h33a/legs/complex/input.pdb")
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
ABAG_RBFE_VISIBLE_GPUS=0 "${ABAG_RBFE}" resume 1dvf-d13-e52-antigen-d-h33a --batch-dir "${B2_BATCH}" --execute 2>&1 | tail -1

echo "[phaseB] $(date --iso-8601=seconds) all done"
