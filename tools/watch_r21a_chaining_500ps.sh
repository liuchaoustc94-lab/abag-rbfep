#!/usr/bin/env bash
# Watcher: r21a chaining+500ps pilot -> when all 3 rep scripts finish, re-BAR.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
JOB="${ROOT}/runs/real_cases/3hfm_patel_dssb_pr4/jobs/3hfm-patel-2021-antigen-y-r21a"
BATCH="${ROOT}/runs/real_cases/3hfm_patel_dssb_pr4"
LOG="${ROOT}/runs/manual/r21a_chaining_500ps_watch.log"
export PATH="${ROOT}/.venv/bin:${PATH}"
while pgrep -f "sample_dssb_rep0" > /dev/null; do sleep 300; done
echo "$(date --iso-8601=seconds) rep scripts finished" >> "${LOG}"
n=$(find "${JOB}/legs" -path "*lambda*" -name dhdl.xvg | wc -l)
echo "dhdl windows: ${n}/72" >> "${LOG}"
rm -f "${JOB}"/stages/{bar,qc,report}.json "${JOB}"/results/*.json 2>/dev/null
find "${JOB}/legs" -maxdepth 3 -name bar -type d -exec rm -rf {} + 2>/dev/null
ABAG_RBFE_VISIBLE_GPUS=0 "${ROOT}/.venv/bin/abag-rbfe" resume 3hfm-patel-2021-antigen-y-r21a \
  --batch-dir "${BATCH}" --execute >> "${LOG}" 2>&1
if [ -f "${JOB}/results/ddg_summary.json" ]; then
  "${ROOT}/.venv/bin/python" -c "
import json
d=json.load(open('${JOB}/results/ddg_summary.json'))
print(f\"CHAINING+500ps ddG={d['ddg_kcal_mol']:.2f} rep_range={d.get('ddg_repeat_range_kcal_mol',-1):.2f} (chaining+50ps: 8.43/7.37; shared-start: 20.87/12.81; exp 0.90)\")
" >> "${LOG}" 2>&1
fi
echo "$(date --iso-8601=seconds) done" >> "${LOG}"
