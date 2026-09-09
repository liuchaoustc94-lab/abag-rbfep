#!/usr/bin/env bash
# Per-window hydrated seeding rollout (2026-09-08): 1dqj y50a + 1vfb w52a.
# Protocol (validated on 3be1 w50a, 9.73 -> 5.93): only high vdw-lambda windows
# (13-15, vdw>=0.82) get hydrated-start reruns; all other windows untouched.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
export PATH="${ROOT}/.venv/bin:${PATH}"
LOG="${ROOT}/runs/manual/perwindow_seed_rollout_20260908.log"
log() { echo "[seed-rollout] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

run_job() {
  local batch="$1" job="$2" gpu="$3" pick2="$4" pick3="$5"
  local jd="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811/${batch}/jobs/${job}"
  local rep1="${jd}/legs/complex/rep01"
  log "== ${job}: 提帧 (rep02<-${pick2}ps, rep03<-${pick3}ps)"
  for spec in "rep02 ${pick2}" "rep03 ${pick3}"; do
    set -- ${spec}
    local tgt="${jd}/legs/complex/$1/equilibration"
    printf 'System\nSystem\n' | "${GMX}" trjconv -f "${rep1}/hydra_hotfree.xtc" -s "${rep1}/hydra_hotfree.tpr" \
      -o "${tgt}/npt_seed_hydrated.gro" -dump "$2" -pbc mol -ur compact > /dev/null 2>&1
    [ -s "${tgt}/npt_seed_hydrated.gro" ] || { log "== ${job}: 帧提取失败 $1"; return 1; }
  done
  # 高 vdw-λ 窗口换水合起点重跑
  for rep in rep02 rep03; do
    local rep_dir="${jd}/legs/complex/${rep}"
    local seed="${rep_dir}/equilibration/npt_seed_hydrated.gro"
    for w in 013 014 015; do
      local wdir="${rep_dir}/lambda_${w}"
      [ -d "${wdir}" ] || continue
      (
        set -e
        export GMXLIB="${jd}/artifacts/gmxlib" GMX_MAXBACKUP="-1" CUDA_VISIBLE_DEVICES="${gpu}"
        cd "${jd}"
        rm -f "${wdir}/pre_relax.tpr" "${wdir}/pre_relax.gro" "${wdir}/pre_md.tpr" "${wdir}/pre_md.gro" \
              "${wdir}/pre_md.cpt" "${wdir}/topol.tpr" "${wdir}/dhdl.xvg" "${wdir}/md.gro" "${wdir}/md.log"
        CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${wdir}/pre_relax.mdp" -c "${seed}" -p "${rep_dir}/system.top" -o "${wdir}/pre_relax.tpr" -maxwarn 2
        "${GMX}" mdrun -s "${wdir}/pre_relax.tpr" -deffnm "${wdir}/pre_relax" -ntmpi 1 -ntomp 4
        CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${wdir}/pre_md.mdp" -c "${wdir}/pre_relax.gro" -p "${rep_dir}/system.top" -o "${wdir}/pre_md.tpr" -maxwarn 2
        "${GMX}" mdrun -s "${wdir}/pre_md.tpr" -deffnm "${wdir}/pre_md" -ntmpi 1 -ntomp 4
        CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${wdir}/production.mdp" -c "${wdir}/pre_md.gro" -t "${wdir}/pre_md.cpt" -p "${rep_dir}/system.top" -o "${wdir}/topol.tpr" -maxwarn 2
        "${GMX}" mdrun -s "${wdir}/topol.tpr" -deffnm "${wdir}/md" -dhdl "${wdir}/dhdl.xvg" -ntmpi 1 -ntomp 4
      ) >> "${LOG}" 2>&1
      rc=$?
      if [ ${rc} -eq 0 ] && [ -s "${wdir}/dhdl.xvg" ]; then
        echo "hydrated start (B-arm hotfree frame)" > "${wdir}/HYDRATED_START.txt"
        log "== ${job} ${rep}/lambda_${w} OK"
      else
        log "== ${job} ${rep}/lambda_${w} FAILED rc=${rc}"
        return 1
      fi
    done
  done
  # 重算 bar
  rm -f "${jd}"/stages/{bar,qc,report}.json "${jd}"/results/*.json
  find "${jd}/legs" -maxdepth 3 -name bar -type d -exec rm -rf {} + 2>/dev/null
  ABAG_RBFE_VISIBLE_GPUS="${gpu}" "${ABAG_RBFE}" resume "${job}" --batch-dir "${jd%/*}" --execute >> "${LOG}" 2>&1
  if [ -f "${jd}/results/ddg_summary.json" ]; then
    "${ROOT}/.venv/bin/python" -c "
import json
d=json.load(open('${jd}/results/ddg_summary.json'))
print(f\"SEEDED ${job}: ddG={d['ddg_kcal_mol']:.2f} rep_range={d.get('ddg_repeat_range_kcal_mol',-1):.2f}\")
" | tee -a "${LOG}"
  fi
}

run_job abbind_1dqj_core_v1 1dqj-antibody-l-y50a 1 9350 14450 &
run_job abbind_1vfb_core_v1 1vfb-antibody-h-w52a 2 6200 10725 &
wait
log "rollout done (对照: 1dqj y50a 旧 9.95/exp 2.70; 1vfb w52a 旧 4.58/exp 0.90; w50a 参照 9.73->5.93/exp 1.40)"
