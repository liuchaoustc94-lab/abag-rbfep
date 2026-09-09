#!/usr/bin/env bash
# P3-a pilot round 3: targeted fixes per canary.
# y50l: grompp on CPU (CUDA_VISIBLE_DEVICES="" for grompp only) to dodge the
#       transient CUDA #700; production unchanged on GPU.
# w50f: pregrown.gro exists; chain with retry3 (half-dt) mdps.
# d31e: 8 substeps x 2000 EM steps (slower ramp).
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
PY="${ROOT}/.venv/bin/python"
LOG="${ROOT}/runs/manual/pregrow_pilot3_20260826.log"
log() { echo "[pilot3] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

run_chain() {  # jd leg rep win gpu pre_md_mdp production_mdp
  local jd="$1" leg="$2" rep="$3" win="$4" gpu="$5" pmmdp="$6" prodmdp="$7"
  local tag="$(basename "${jd}")/${win}"
  local rep_dir="${jd}/legs/${leg}/${rep}"
  local wdir="${rep_dir}/${win}"
  (
    set -e
    export GMXLIB="${jd}/artifacts/gmxlib" GMX_MAXBACKUP="-1"
    cd "${jd}"
    rm -f "${wdir}/pre_md.tpr" "${wdir}/pre_md.gro" "${wdir}/pre_md.cpt" "${wdir}/topol.tpr" \
          "${wdir}/dhdl.xvg" "${wdir}/md.gro" "${wdir}/md.log"
    # grompp strictly on CPU (CUDA #700 collateral fix)
    CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${pmmdp}" -c "${wdir}/pregrown.gro" -p "${rep_dir}/system.top" \
      -o "${wdir}/pre_md.tpr" -maxwarn 2
    CUDA_VISIBLE_DEVICES="${gpu}" "${GMX}" mdrun -s "${wdir}/pre_md.tpr" -deffnm "${wdir}/pre_md" -ntmpi 1 -ntomp 4
    CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${prodmdp}" -c "${wdir}/pre_md.gro" -t "${wdir}/pre_md.cpt" \
      -p "${rep_dir}/system.top" -o "${wdir}/topol.tpr" -maxwarn 2
    CUDA_VISIBLE_DEVICES="${gpu}" "${GMX}" mdrun -s "${wdir}/topol.tpr" -deffnm "${wdir}/md" -dhdl "${wdir}/dhdl.xvg" -ntmpi 1 -ntomp 4
  ) >> "${LOG}" 2>&1
  rc=$?
  if [ ${rc} -eq 0 ] && [ -s "${wdir}/dhdl.xvg" ] && [ -s "${wdir}/md.gro" ] && [ -s "${wdir}/md.log" ] && [ -s "${wdir}/topol.tpr" ]; then
    log "== ${tag}: RESCUED"; return 0
  fi
  log "== ${tag}: CHAIN FAILED rc=${rc}"; return 1
}

jd_y="${ROOT}/runs/real_cases/3hfm_patel_align_20260817/jobs/3hfm-patel-2021-antibody-l-y50l"
run_chain "${jd_y}" complex rep02 lambda_007 1 \
  "${jd_y}/legs/complex/rep02/lambda_007/pre_md.mdp" \
  "${jd_y}/legs/complex/rep02/lambda_007/production.mdp" &

jd_w="${ROOT}/runs/real_cases/p1c_stepdown_3be1_20260826/jobs/3be1-antibody-l-w50f"
run_chain "${jd_w}" complex rep01 lambda_005 2 \
  "${jd_w}/legs/complex/rep01/lambda_005/pre_md.retry3.mdp" \
  "${jd_w}/legs/complex/rep01/lambda_005/production.retry3.mdp" &

jd_d="${ROOT}/runs/acceptance/acc_2bdn_v1_20260818/jobs/acc-2bdn-antibody-h-d31e"
(
  if env GMXLIB="${jd_d}/artifacts/gmxlib" CUDA_VISIBLE_DEVICES=2 "${PY}" "${ROOT}/tools/pregrow_window.py" \
      --window-dir "${jd_d}/legs/complex/rep01/lambda_002" \
      --start-gro "${jd_d}/legs/complex/rep01/equilibration/npt.gro" \
      --repeat-top "${jd_d}/legs/complex/rep01/system.top" \
      --gmx "${GMX}" --gmxlib "${jd_d}/artifacts/gmxlib" \
      --substeps 8 --em-steps 2000 >> "${LOG}" 2>&1; then
    run_chain "${jd_d}" complex rep01 lambda_002 2 \
      "${jd_d}/legs/complex/rep01/lambda_002/pre_md.mdp" \
      "${jd_d}/legs/complex/rep01/lambda_002/production.mdp"
  else
    log "== d31e lambda_002: PREGROW FAILED (8 substeps)"
  fi
) &

wait
log "pilot3 done"
