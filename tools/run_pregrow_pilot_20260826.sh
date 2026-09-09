#!/usr/bin/env bash
# P3-a pilot (2026-08-26): alchembed-style lambda sub-step pre-growth on the
# three canary crash windows, then run the window's standard chain manually
# starting from the pregrown structure. If the completion artifacts appear,
# the window is rescued and a normal resume can finish the job.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
PY="${ROOT}/.venv/bin/python"
LOG="${ROOT}/runs/manual/pregrow_pilot_20260826.log"

log() { echo "[pregrow-pilot] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

# canary: job_dir|leg|rep|window|gpu
CANARIES=(
  "${ROOT}/runs/real_cases/3hfm_patel_align_20260817/jobs/3hfm-patel-2021-antibody-l-y50l|complex|rep02|lambda_007|1"
  "${ROOT}/runs/acceptance/acc_2bdn_v1_20260818/jobs/acc-2bdn-antibody-h-d31e|complex|rep01|lambda_002|2"
  "${ROOT}/runs/real_cases/p1c_stepdown_3be1_20260826/jobs/3be1-antibody-l-w50f|complex|rep01|lambda_005|2"
)

run_canary() {
  local spec="$1"
  local jd leg rep win gpu
  IFS='|' read -r jd leg rep win gpu <<< "${spec}"
  local tag="$(basename "${jd}")/${leg}/${rep}/${win}"
  local rep_dir="${jd}/legs/${leg}/${rep}"
  local wdir="${rep_dir}/${win}"
  local env_prefix="GMXLIB=${jd}/artifacts/gmxlib GMX_MAXBACKUP=-1 CUDA_VISIBLE_DEVICES=${gpu}"

  log "== ${tag}: pre-growth"
  if ! env GMXLIB="${jd}/artifacts/gmxlib" CUDA_VISIBLE_DEVICES="${gpu}" "${PY}" "${ROOT}/tools/pregrow_window.py" \
      --window-dir "${wdir}" --start-gro "${rep_dir}/equilibration/npt.gro" \
      --repeat-top "${rep_dir}/system.top" --gmx "${GMX}" --gmxlib "${jd}/artifacts/gmxlib" \
      --substeps 4 --em-steps 1000 >> "${LOG}" 2>&1; then
    log "== ${tag}: PREGROW FAILED"
    return 1
  fi

  # standard chain from the pregrown structure (skip pre_relax)
  log "== ${tag}: standard chain from pregrown.gro"
  (
    set -e
    export GMXLIB="${jd}/artifacts/gmxlib" GMX_MAXBACKUP="-1" CUDA_VISIBLE_DEVICES="${gpu}"
    cd "${jd}"
    rm -f "${wdir}/pre_md.tpr" "${wdir}/pre_md.gro" "${wdir}/pre_md.cpt" "${wdir}/topol.tpr" \
          "${wdir}/dhdl.xvg" "${wdir}/md.gro" "${wdir}/md.log"
    "${GMX}" grompp -f "${wdir}/pre_md.mdp" -c "${wdir}/pregrown.gro" -p "${rep_dir}/system.top" \
      -o "${wdir}/pre_md.tpr" -maxwarn 2
    "${GMX}" mdrun -s "${wdir}/pre_md.tpr" -deffnm "${wdir}/pre_md" -ntmpi 1 -ntomp 4
    "${GMX}" grompp -f "${wdir}/production.mdp" -c "${wdir}/pre_md.gro" -t "${wdir}/pre_md.cpt" \
      -p "${rep_dir}/system.top" -o "${wdir}/topol.tpr" -maxwarn 2
    "${GMX}" mdrun -s "${wdir}/topol.tpr" -deffnm "${wdir}/md" -dhdl "${wdir}/dhdl.xvg" -ntmpi 1 -ntomp 4
  ) >> "${LOG}" 2>&1
  rc=$?
  if [ ${rc} -eq 0 ] && [ -s "${wdir}/dhdl.xvg" ] && [ -s "${wdir}/md.gro" ] && [ -s "${wdir}/md.log" ] && [ -s "${wdir}/topol.tpr" ]; then
    log "== ${tag}: RESCUED (window complete)"
    return 0
  fi
  log "== ${tag}: CHAIN FAILED rc=${rc}"
  return 1
}

pids=()
for spec in "${CANARIES[@]}"; do
  run_canary "${spec}" &
  pids+=($!)
done
rc_all=0
for p in "${pids[@]}"; do wait "${p}" || rc_all=1; done
log "pilot done, overall_rc=${rc_all}"
