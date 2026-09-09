#!/usr/bin/env bash
# Restrained-growth pilot (2026-08-31): P3-a follow-up after the plain
# pre-growth pilot showed EM-level relaxation is insufficient (LINCS storms
# from step 0 of pre_md). Root suspicion: the standard pre_md is a brutal cold
# start (0.2 ps, constraints=none, zero velocities inherited from EM output).
# This pilot runs, per canary window:
#   pregrown.gro -> pre_md_posres (POSRES on protein heavy atoms fc=1000,
#                    ref=pregrown.gro, gen-vel 310K, constraints=h-bonds,
#                    dt=0.001, 2 ps) -> standard production (unrestrained)
# Success criterion: window completion artifacts (dhdl.xvg, md.gro, md.log,
# topol.tpr) appear with rc=0.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
PY="${ROOT}/.venv/bin/python"
LOG="${ROOT}/runs/manual/restrained_growth_pilot_20260831.log"
log() { echo "[rg-pilot] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

# canary: job_dir|leg|rep|window|gpu
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
  [ -s "${wdir}/pregrown.gro" ] || { log "== ${tag}: no pregrown.gro"; return 1; }

  # 1) window-local top with posre includes in moleculetype scope
  "${PY}" - "${rep_dir}" <<'PYEOF'
import sys, re
from pathlib import Path
rep = Path(sys.argv[1])
src = (rep / "system.top").read_text()
posres = {p.name for p in rep.glob("posre_Protein_chain_*.itp")}
out = []
for line in src.splitlines():
    out.append(line)
    m = re.search(r'#include "(?:pmx_)?topol_Protein_chain_(\w+)\.itp"', line)
    if m:
        pf = f"posre_Protein_chain_{m.group(1)}.itp"
        if pf in posres:
            out += ["", "#ifdef POSRES", f'#include "{pf}"', "#endif", ""]
(rep / "system_rg.top").write_text("\n".join(out) + "\n")
PYEOF

  # 2) restrained pre_md mdp derived from the window's pre_md.mdp
  "${PY}" - "${wdir}" <<'PYEOF'
import sys
from pathlib import Path
w = Path(sys.argv[1])
src = (w / "pre_md.mdp").read_text()
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
(w / "pre_md_posres.mdp").write_text(text + "\n")
PYEOF

  # 3) restrained pre_md -> standard production
  log "== ${tag}: restrained pre_md"
  (
    set -e
    export GMXLIB="${jd}/artifacts/gmxlib" GMX_MAXBACKUP="-1"
    cd "${jd}"
    rm -f "${wdir}/pre_md.tpr" "${wdir}/pre_md.gro" "${wdir}/pre_md.cpt" "${wdir}/topol.tpr" \
          "${wdir}/dhdl.xvg" "${wdir}/md.gro" "${wdir}/md.log"
    CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${wdir}/pre_md_posres.mdp" -c "${wdir}/pregrown.gro" \
      -r "${wdir}/pregrown.gro" -p "${rep_dir}/system_rg.top" -o "${wdir}/pre_md.tpr" -maxwarn 2
    CUDA_VISIBLE_DEVICES="${gpu}" "${GMX}" mdrun -s "${wdir}/pre_md.tpr" -deffnm "${wdir}/pre_md" -ntmpi 1 -ntomp 4
    CUDA_VISIBLE_DEVICES="" "${GMX}" grompp -f "${wdir}/production.mdp" -c "${wdir}/pre_md.gro" \
      -t "${wdir}/pre_md.cpt" -p "${rep_dir}/system.top" -o "${wdir}/topol.tpr" -maxwarn 2
    CUDA_VISIBLE_DEVICES="${gpu}" "${GMX}" mdrun -s "${wdir}/topol.tpr" -deffnm "${wdir}/md" \
      -dhdl "${wdir}/dhdl.xvg" -ntmpi 1 -ntomp 4
  ) >> "${LOG}" 2>&1
  rc=$?
  if [ ${rc} -eq 0 ] && [ -s "${wdir}/dhdl.xvg" ] && [ -s "${wdir}/md.gro" ] && [ -s "${wdir}/md.log" ] && [ -s "${wdir}/topol.tpr" ]; then
    log "== ${tag}: RESCUED (restrained growth)"
    return 0
  fi
  log "== ${tag}: FAILED rc=${rc}"
  return 1
}

pids=()
for spec in "${CANARIES[@]}"; do
  run_canary "${spec}" &
  pids+=($!)
done
rc_all=0
for p in "${pids[@]}"; do wait "${p}" || rc_all=1; done
log "restrained-growth pilot done, overall_rc=${rc_all}"
