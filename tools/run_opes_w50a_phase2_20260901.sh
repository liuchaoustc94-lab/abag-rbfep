#!/usr/bin/env bash
# P1-b phase 2 (2026-09-01): per-window hydrated-start pilot on 3be1 w50a.
# After phase 1 (rep02/rep03 rerun all-standard) completes, replace the high
# vdw-lambda windows (13-15, vdw >= 0.82, W ring ~dummy) of rep02/rep03 with
# runs started from the hydrated MetaD frames, then re-BAR. rep01 and the apo
# leg stay untouched as the standard reference.
# Verdict logic: ddG moves 9.73 -> toward exp 1.40 with small rep spread =
# hydration hypothesis + per-window delivery both validated.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
JOB="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1/jobs/3be1-antibody-l-w50a"
BATCH="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1"
LOG="${ROOT}/runs/manual/opes_w50a_phase2_20260901.log"
log() { echo "[p1b-phase2] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

# wait for phase 1 (all-standard rerun) to finish: ddg_summary regenerated
while [ ! -f "${JOB}/results/ddg_summary.json" ]; do sleep 300; done
# extra guard: phase-1 resume process must be gone
while pgrep -f "resume 3be1-antibody-l-w50a" > /dev/null; do sleep 60; done
log "phase 1 done; replacing high-lambda windows with hydrated starts"

for rep in rep02 rep03; do
  rep_dir="${JOB}/legs/complex/${rep}"
  start_gro="${rep_dir}/equilibration/npt_seed_hydrated.gro"
  [ -s "${start_gro}" ] || { log "no hydrated seed for ${rep}"; continue; }
  for w in 013 014 015; do
    wdir="${rep_dir}/lambda_${w}"
    [ -d "${wdir}" ] || continue
    log "rerun ${rep}/lambda_${w} from hydrated frame"
    (
      set -e
      export GMXLIB="${JOB}/artifacts/gmxlib" GMX_MAXBACKUP="-1" CUDA_VISIBLE_DEVICES=0
      cd "${JOB}"
      rm -f "${wdir}/pre_relax.tpr" "${wdir}/pre_relax.gro" "${wdir}/pre_md.tpr" "${wdir}/pre_md.gro" \
            "${wdir}/pre_md.cpt" "${wdir}/topol.tpr" "${wdir}/dhdl.xvg" "${wdir}/md.gro" "${wdir}/md.log"
      "${GMX}" grompp -f "${wdir}/pre_relax.mdp" -c "${start_gro}" -p "${rep_dir}/system.top" \
        -o "${wdir}/pre_relax.tpr" -maxwarn 2
      "${GMX}" mdrun -s "${wdir}/pre_relax.tpr" -deffnm "${wdir}/pre_relax" -ntmpi 1 -ntomp 4
      "${GMX}" grompp -f "${wdir}/pre_md.mdp" -c "${wdir}/pre_relax.gro" -p "${rep_dir}/system.top" \
        -o "${wdir}/pre_md.tpr" -maxwarn 2
      "${GMX}" mdrun -s "${wdir}/pre_md.tpr" -deffnm "${wdir}/pre_md" -ntmpi 1 -ntomp 4
      "${GMX}" grompp -f "${wdir}/production.mdp" -c "${wdir}/pre_md.gro" -t "${wdir}/pre_md.cpt" \
        -p "${rep_dir}/system.top" -o "${wdir}/topol.tpr" -maxwarn 2
      "${GMX}" mdrun -s "${wdir}/topol.tpr" -deffnm "${wdir}/md" -dhdl "${wdir}/dhdl.xvg" -ntmpi 1 -ntomp 4
    ) >> "${LOG}" 2>&1
    rc=$?
    if [ ${rc} -eq 0 ] && [ -s "${wdir}/dhdl.xvg" ]; then
      echo "hydrated start from ${start_gro}" > "${wdir}/HYDRATED_START.txt"
      log "  ${rep}/lambda_${w} OK"
    else
      log "  ${rep}/lambda_${w} FAILED rc=${rc} (window left incomplete)"
    fi
  done
done

# re-BAR with the mixed ensemble (rep01 standard; rep02/03 windows 0-12 standard + 13-15 hydrated)
rm -f "${JOB}"/stages/{bar,qc,report}.json "${JOB}"/results/*.json
find "${JOB}/legs" -maxdepth 3 -name bar -type d -exec rm -rf {} +
log "re-BAR"
ABAG_RBFE_VISIBLE_GPUS=0 "${ABAG_RBFE}" resume 3be1-antibody-l-w50a --batch-dir "${BATCH}" --execute \
  >> "${ROOT}/runs/manual/opes_w50a_phase2_rebar.log" 2>&1
if [ -f "${JOB}/results/ddg_summary.json" ]; then
  "${ROOT}/.venv/bin/python" -c "
import json
d = json.load(open('${JOB}/results/ddg_summary.json'))
print(f\"PHASE2 ddG={d['ddg_kcal_mol']:.2f} rep_range={d.get('ddg_repeat_range_kcal_mol',-1):.2f} (all-standard 9.73, shared-hydrated -3.87, exp 1.40)\")
" | tee -a "${LOG}"
else
  log "rebar produced no ddg_summary"
fi
log "phase2 done"
