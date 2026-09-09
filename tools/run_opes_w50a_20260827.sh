#!/usr/bin/env bash
# P1-b pilot (2026-08-27): OPES on cavity water count at the 3be1 W50A mutant
# endpoint (complex leg, lambda state 15). 25 ns, posres on protein heavy atoms,
# PLUMED OPES_METAD on COORDINATION(cavity COM, water O). Goal: hydrated frames
# for endpoint-ensemble injection.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
JOB="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1/jobs/3be1-antibody-l-w50a"
REP="${JOB}/legs/complex/rep01"
cd "${REP}"
export GMXLIB="${JOB}/artifacts/gmxlib" GMX_MAXBACKUP=-1
export PLUMED_KERNEL="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/plumed-mpi/lib/libplumedKernel.so"
CUDA_VISIBLE_DEVICES=0 "${GMX}" mdrun -s opes.tpr -deffnm opes -plumed opes_plumed.dat -ntmpi 1 -ntomp 8 >> opes_run.log 2>&1
rc=$?
# 只在成功时写完成标记（2026-08-31 bug 修复：此前 rc=1 也写，导致 stage2
# watcher 被早期失败启动的残留标记提前触发、读到只有 1 帧的 colvar 而误判）
if [ ${rc} -eq 0 ]; then
  echo "OPES done rc=0 $(date --iso-8601=seconds)" >> "${ROOT}/runs/manual/opes_w50a.status"
else
  echo "OPES failed rc=${rc} $(date --iso-8601=seconds)" >> "${ROOT}/runs/manual/opes_w50a.status"
fi
