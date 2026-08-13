#!/usr/bin/env bash
# Basin-diversity control experiment (2026-08-12).
# Question: does endpoint conformational basin diversity move ddG for the worst
# aromatic-deletion jobs (1DVF y102a, y49a)? Method: replace the 3 repeats'
# equilibration endpoints with 3 distant frames (1.67/3.33/5 ns) from a free
# 5 ns WT trajectory per leg, then run the standard FEP. If per-repeat ddG
# spread >> 1 kcal/mol, basin diversity matters -> RID/REST2 worth it.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
BATCH="${ROOT}/runs/real_cases/1dvf_basin_test_20260812"
export PATH="${ROOT}/.venv/bin:${PATH}"

echo "[basin] plan batch"
cat > /tmp/1dvf_basin_mutations.csv <<'CSV'
mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side
antigen_d_y102a,D,102,,Y,A,antigen
antibody_a_y49a,A,49,,Y,A,antibody
CSV
"${ABAG_RBFE}" batch plan \
  --system "${ROOT}/examples/real_cases/1dvf/system.yml" \
  --mutations /tmp/1dvf_basin_mutations.csv \
  --protocol "${ROOT}/benchmarks/ab_bind/protocol.validation_priority.yml" \
  --batch-id 1dvf_basin_test_20260812 \
  --runs-root "${ROOT}/runs/real_cases"

for job in 1dvf-d13-e52-antigen-d-y102a 1dvf-d13-e52-antibody-a-y49a; do
  echo "[basin] equilibrate ${job}"
  "${ABAG_RBFE}" run "${job}" --batch-dir "${BATCH}" --execute --to-stage equilibrate 2>&1 | tail -1
done

# 5 ns free-MD extension per leg from rep01 NPT endpoint; frames -> rep01/02/03
for job in 1dvf-d13-e52-antigen-d-y102a 1dvf-d13-e52-antibody-a-y49a; do
  for leg in complex apo; do
    rep="${BATCH}/jobs/${job}/legs/${leg}/rep01"
    eq="${rep}/equilibration"
    [ -f "${eq}/npt.gro" ] || { echo "[basin] skip ${job}/${leg}: no npt.gro"; continue; }
    echo "[basin] free-MD 5ns ${job}/${leg}"
    cat > "${eq}/freemd.mdp" <<EOF
integrator = sd
dt = 0.002
nsteps = 2500000
ld-seed = 42001
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
nstxout-compressed = 2500
compressed-x-grps = System
EOF
    "${GMX}" grompp -f "${eq}/freemd.mdp" -c "${eq}/npt.gro" -p "${rep}/system.top" \
      -o "${eq}/freemd.tpr" -maxwarn 2 > "${eq}/freemd_grompp.log" 2>&1
    ABAG_RBFE_VISIBLE_GPUS="${RID_GPU:-0}" "${GMX}" mdrun -s "${eq}/freemd.tpr" \
      -deffnm "${eq}/freemd" -ntmpi 1 -ntomp 4 >> "${eq}/freemd_grompp.log" 2>&1
    # 提取 3 帧并替换三个 rep 的 NPT 终点
    i=1
    for frame_ps in 1667 3333 5000; do
      target="${BATCH}/jobs/${job}/legs/${leg}/rep0${i}/equilibration"
      echo "System" | "${GMX}" trjconv -f "${eq}/freemd.xtc" -s "${eq}/freemd.tpr" \
        -o "${target}/npt_seed.gro" -dump "${frame_ps}" -pbc mol -ur compact -center > /dev/null 2>&1
      if [ -s "${target}/npt_seed.gro" ]; then
        cp "${target}/npt.gro" "${target}/npt.gro.orig" 2>/dev/null || true
        cp "${target}/npt_seed.gro" "${target}/npt.gro"
        echo "[basin] rep0${i} ${job}/${leg} seeded from ${frame_ps} ps"
      fi
      i=$((i+1))
    done
  done
done

echo "[basin] run FEP stages for both jobs"
for job in 1dvf-d13-e52-antigen-d-y102a 1dvf-d13-e52-antibody-a-y49a; do
  "${ABAG_RBFE}" resume "${job}" --batch-dir "${BATCH}" --execute 2>&1 | tail -1
done
echo "[basin] all done"
