#!/usr/bin/env bash
# Restrained-growth pilot v2 (2026-08-31): full couple-intramol=yes consistency.
# v1 mixed a ci=yes pregrown EM structure with ci=no pre_md/production mdps
# (different Hamiltonian -> LINCS storm). v2: build the restrained pre_md from
# pre_md.retry2.mdp (ci=yes) and use production.retry2.mdp (ci=yes).
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
PY="${ROOT}/.venv/bin/python"
LOG="${ROOT}/runs/manual/restrained_growth_pilot2_20260831.log"
log() { echo "[rg2] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

CANARIES=(
  "${ROOT}/runs/real_cases/3hfm_patel_align_20260817/jobs/3hfm-patel-2021-antibody-l-y50l|complex|rep02|lambda_007|1"
  "${ROOT}/runs/acceptance/acc_2bdn_v1_20260818/jobs/acc-2bdn-antibody-h-d31e|complex|rep01|lambda_002|2"
  "${ROOT}/runs/real_cases/p1c_stepdown_3be1_20260826/jobs/3be1-antibody-l-w50f|complex|rep01|lambda_005|0"
)

run_canary() {
  local spec="$1"
  local jd leg rep win gpu
  IFS='|' read -r jd leg rep win gpu <<< "${spec}"
  local tag="$(basename "${jd}")/${win}"
  local rep_dir="${jd}/legs/${leg}/${rep}"
  local wdir="${rep_dir}/${win}"
  local pre_src="${wdir}/pre_md.retry2.mdp"
  local prod_src="${wdir}/production.retry2.mdp"
  [ -s "${wdir}/pregrown.gro" ] || { log "== ${tag}: no pregrown.gro"; return 1; }
  [ -f "${pre_src}" ] || { log "== ${tag}: no pre_md.retry2.mdp"; return 1; }
  [ -f "${prod_src}" ] || { log "== ${tag}: no production.retry2.mdp"; return 1; }
  [ -f "${rep_dir}/system_rg.top" ] || { log "== ${tag}: no system_rg.top"; return 1; }

  "${PY}" - "${wdir}" "${pre_src}" <<'PYEOF'
import sys
from pathlib import Path
w = Path(sys.argv[1]); src = Path(sys.argv[2]).read_text()
lines = []
for line in src.splitlines():
    key = line.split("=")[0].strip() if "=" in line else ""
    if key == "define":
        lines.append("define                  = -DFLEXIBLE -DPOSRES")
    elif key == "dt":
        lines.append("dt                      = 0.001")
    elif key == "nsteps":
        lines.append("nsteps                  = 2000")
    elif key == "constraints":
        lines.append("constraints             = h-bonds")
    elif key == "gen-vel":
        lines.append("gen-vel                 = yes")
    elif key == "continuation":
        lines.append("continuation            = no")
    else:
        lines.append(line)
text = "\n".join(lines)
if "gen-temp" not in text:
    text = text.replace("gen-vel                 = yes", "gen-vel                 = yes\ngen-temp                = 310.0")
(w / "pre_md_posres_ci.mdp").write_text(text + "\n")
PYEOF

  log "== ${tag}: ci-consistent restrained chain"
  (
    set -e
    export GMXLIB="${jd}/artifacts/gmxlib" GMX_MAXBACKUP="-1"
    cd "${jd}"
    rm -f "${wdir}/pre_md.tpr" "${wdir}/pre_md.gro" "${wdir}/pre_md.cpt" "${wdir}/topol.tpr" \
          "${wdir}/dhdl.xvg" "${wdir}/md.gro" "${wdir}/md.log"
    CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${wdir}/pre_md_posres_ci.mdp" -c "${wdir}/pregrown.gro" \
      -r "${wdir}/pregrown.gro" -p "${rep_dir}/system_rg.top" -o "${wdir}/pre_md.tpr" -maxwarn 2
    CUDA_VISIBLE_DEVICES="${gpu}" "${GMX}" mdrun -s "${wdir}/pre_md.tpr" -deffnm "${wdir}/pre_md" -ntmpi 1 -ntomp 4
    CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${prod_src}" -c "${wdir}/pre_md.gro" \
      -t "${wdir}/pre_md.cpt" -p "${rep_dir}/system.top" -o "${wdir}/topol.tpr" -maxwarn 2
    CUDA_VISIBLE_DEVICES="${gpu}" "${GMX}" mdrun -s "${wdir}/topol.tpr" -deffnm "${wdir}/md" \
      -dhdl "${wdir}/dhdl.xvg" -ntmpi 1 -ntomp 4
  ) >> "${LOG}" 2>&1
  rc=$?
  if [ ${rc} -eq 0 ] && [ -s "${wdir}/dhdl.xvg" ] && [ -s "${wdir}/md.gro" ] && [ -s "${wdir}/md.log" ] && [ -s "${wdir}/topol.tpr" ]; then
    log "== ${tag}: RESCUED (ci-consistent restrained growth)"
    return 0
  fi
  log "== ${tag}: FAILED rc=${rc}"
  return 1
}

pids=()
for spec in "${CANARIES[@]}"; do run_canary "${spec}" & pids+=($!); done
rc_all=0
for p in "${pids[@]}"; do wait "${p}" || rc_all=1; done
log "pilot2 done, overall_rc=${rc_all}"
