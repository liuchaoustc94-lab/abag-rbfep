#!/usr/bin/env bash
# E4: Patel 2021 protocol-alignment experiment on 3HFM (2026-08-17).
# Question: does Patel-style deep endpoint equilibration (40 ns + multi-frame
# seeds) close the gap to Patel's reported R^2=0.81 on the same toolchain?
# Design: 3 charge-conserving 3HFM mutations (y20f/w98f/y50l, exp -0.48/3.25/4.39),
# 40 ns free-MD per leg from the standard equilibration endpoint, frames at
# 12.5/25/37.5 ns injected as rep01/02/03 equilibration endpoints, then the
# standard FEP chain (direction-aware decoupled schedule, adaptive lambda).
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
BATCH_DIR="${ROOT}/runs/real_cases/3hfm_patel_align_20260817"
export PATH="${ROOT}/.venv/bin:${PATH}"

echo "[E4] plan 3HFM jobs"
if [ ! -d "${BATCH_DIR}" ]; then
  printf 'mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side\npatel_3hfm_y20f,Y,20,,Y,F,antigen\npatel_3hfm_w98f,H,98,,W,F,antibody\npatel_3hfm_y50l,L,50,,Y,L,antibody\n' > /tmp/patel3_e4.csv
  "${ABAG_RBFE}" batch plan --system "${ROOT}/benchmarks/patel_2021_3hfm/system.yml" \
    --mutations /tmp/patel3_e4.csv --protocol "${ROOT}/benchmarks/ab_bind/protocol.validation_priority.yml" \
    --batch-id 3hfm_patel_align_20260817 --runs-root "${ROOT}/runs/real_cases"
fi

JOBS=(3hfm-patel-2021-antigen-y-y20f 3hfm-patel-2021-antibody-h-w98f 3hfm-patel-2021-antibody-l-y50l)

echo "[E4] run to equilibrate"
for job in "${JOBS[@]}"; do
  "${ABAG_RBFE}" run "${job}" --batch-dir "${BATCH_DIR}" \
    --execute --to-stage equilibrate 2>&1 | tail -1
done

echo "[E4] 40 ns free-MD per leg"
for job in "${JOBS[@]}"; do
  for leg in complex apo; do
    rep="${BATCH_DIR}/jobs/${job}/legs/${leg}/rep01"
    eq="${rep}/equilibration"
    [ -f "${eq}/npt.gro" ] || { echo "[E4] skip ${job}/${leg}"; continue; }
    if [ ! -f "${eq}/freemd40.xtc" ]; then
      cat > "${eq}/freemd.mdp" <<'MDP'
integrator = sd
dt = 0.002
nsteps = 20000000
ld-seed = 51703
cutoff-scheme = Verlet
verlet-buffer-tolerance = -1
rlist = 1.25
rcoulomb = 1.25
rvdw = 1.25
vdw-type = Cut-off
vdw-modifier = Potential-switch
rvdw-switch = 1.0
coulombtype = PME
pme-order = 4
fourierspacing = 0.12
DispCorr = EnerPres
tcoupl = v-rescale
tc-grps = System
tau-t = 1.0
ref-t = 310.0
pcoupl = C-rescale
pcoupltype = isotropic
tau-p = 2.0
compressibility = 4.5e-5
ref-p = 1.0
constraints = h-bonds
constraint-algorithm = lincs
pbc = xyz
nstxout-compressed = 12500
compressed-x-grps = System
MDP
      "${GMX}" grompp -f "${eq}/freemd.mdp" -c "${eq}/npt.gro" -p "${rep}/system.top" \
        -o "${eq}/freemd40.tpr" -maxwarn 2 > "${eq}/freemd40.log" 2>&1 && \
      ABAG_RBFE_VISIBLE_GPUS="${E4_GPU:-0}" "${GMX}" mdrun -s "${eq}/freemd40.tpr" \
        -deffnm "${eq}/freemd40" -ntmpi 1 -ntomp 4 >> "${eq}/freemd40.log" 2>&1 &
      echo "[E4] launched free-MD ${job}/${leg}"
    fi
  done
done
wait
echo "[E4] free-MD done; injecting frames"
for job in "${JOBS[@]}"; do
  for leg in complex apo; do
    rep1="${BATCH_DIR}/jobs/${job}/legs/${leg}/rep01"
    eq="${rep1}/equilibration"
    [ -f "${eq}/freemd40.xtc" ] || { echo "[E4] no trajectory ${job}/${leg}"; continue; }
    i=1
    for frame_ps in 12500 25000 37500; do
      target="${BATCH_DIR}/jobs/${job}/legs/${leg}/rep0${i}/equilibration"
      printf 'System\nSystem\n' | "${GMX}" trjconv -f "${eq}/freemd40.xtc" -s "${eq}/freemd40.tpr" \
        -o "${target}/npt_seed.gro" -dump "${frame_ps}" -pbc mol -ur compact > /dev/null 2>&1
      if [ -s "${target}/npt_seed.gro" ]; then
        cp "${target}/npt.gro" "${target}/npt.gro.orig" 2>/dev/null || true
        cp "${target}/npt_seed.gro" "${target}/npt.gro"
        echo "[E4] rep0${i} ${job}/${leg} <- ${frame_ps} ps frame"
      fi
      i=$((i+1))
    done
  done
done

echo "[E4] run FEP stages"
for job in "${JOBS[@]}"; do
  "${ABAG_RBFE}" resume "${job}" --batch-dir "${BATCH_DIR}" --execute 2>&1 | tail -1
done
echo "[E4] $(date --iso-8601=seconds) all done"
