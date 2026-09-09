#!/usr/bin/env bash
# P3-a pilot round 2 (2026-08-26): rerun the two failed canaries.
# - d31e: pregrow now uses the couple-intramol=yes EM base (fixes sub02 fatal)
# - y50l: pregrown.gro already exists from round 1; just retry the standard
#   chain (round-1 failure was a transient CUDA #700 in grompp)
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
PY="${ROOT}/.venv/bin/python"
LOG="${ROOT}/runs/manual/pregrow_pilot2_20260826.log"
log() { echo "[pregrow-pilot2] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

run_chain() {
  local jd="$1" leg="$2" rep="$3" win="$4" gpu="$5"
  local tag="$(basename "${jd}")/${leg}/${rep}/${win}"
  local rep_dir="${jd}/legs/${leg}/${rep}"
  local wdir="${rep_dir}/${win}"
  [ -s "${wdir}/pregrown.gro" ] || { log "== ${tag}: no pregrown.gro"; return 1; }
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
    log "== ${tag}: RESCUED"
    return 0
  fi
  log "== ${tag}: CHAIN FAILED rc=${rc}"
  return 1
}

# d31e: fresh pregrow (ci=yes base) + chain
jd_d="${ROOT}/runs/acceptance/acc_2bdn_v1_20260818/jobs/acc-2bdn-antibody-h-d31e"
log "== d31e: pregrow with ci=yes base"
if env GMXLIB="${jd_d}/artifacts/gmxlib" CUDA_VISIBLE_DEVICES=2 "${PY}" "${ROOT}/tools/pregrow_window.py" \
    --window-dir "${jd_d}/legs/complex/rep01/lambda_002" \
    --start-gro "${jd_d}/legs/complex/rep01/equilibration/npt.gro" \
    --repeat-top "${jd_d}/legs/complex/rep01/system.top" \
    --gmx "${GMX}" --gmxlib "${jd_d}/artifacts/gmxlib" \
    --substeps 4 --em-steps 1000 >> "${LOG}" 2>&1; then
  run_chain "${jd_d}" complex rep01 lambda_002 2 &
else
  log "== d31e: PREGROW FAILED again"
fi

# y50l: chain retry only
jd_y="${ROOT}/runs/real_cases/3hfm_patel_align_20260817/jobs/3hfm-patel-2021-antibody-l-y50l"
run_chain "${jd_y}" complex rep02 lambda_007 1 &

wait
log "pilot2 done"
