#!/usr/bin/env bash
# P1-b pilot stage 2 (2026-08-27): when the 25ns MetaD water-CV run finishes,
# extract 2 hydrated frames, inject them as rep02/rep03 seeds of the 3be1 w50a
# job, rerun sample+bar, and log old vs new ddG (old: 9.73, exp: 1.40).
# ISSUE-009 compliance: lambda dirs are KEPT (no rebuild needed) -> only the
# sample/bar/qc/report stage jsons and results are removed.
set -uo pipefail
ROOT="/mnt/data/liuchao/abag-rbfep"
GMX="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa/tools/gromacs-gpu/bin/gmx"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
JOB="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1/jobs/3be1-antibody-l-w50a"
BATCH="${ROOT}/runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1"
REP="${JOB}/legs/complex/rep01"
LOG="${ROOT}/runs/manual/opes_w50a_stage2.log"
log() { echo "[opes-stage2] $(date --iso-8601=seconds) $*" | tee -a "${LOG}"; }

# wait for the metad run: status must report rc=0 AND colvar must reach the
# planned duration (guards against the premature-trigger bug of 2026-08-27,
# where rc=1 residue from failed launches fired the watcher early).
while true; do
  ok_status=0
  grep -q "OPES done rc=0" "${ROOT}/runs/manual/opes_w50a.status" 2>/dev/null && ok_status=1
  ok_cv=0
  if [ -f "${REP}/colvar.dat" ]; then
    last_t=$(awk '!/^#!/ {t=$1} END {print int(t)}' "${REP}/colvar.dat" 2>/dev/null || echo 0)
    [ "${last_t:-0}" -ge 24900 ] && ok_cv=1
  fi
  [ "${ok_status}" -eq 1 ] && [ "${ok_cv}" -eq 1 ] && break
  sleep 300
done
log "metad run finished (verified rc=0 + colvar complete), analyzing colvar"

# pick 2 hydrated frames from distinct epochs
"${ROOT}/.venv/bin/python" - <<'PYEOF'
import json, subprocess
colvar = "runs/benchmarks/abbind_fit_newprotocol_20260811/abbind_3be1_core_v1/jobs/3be1-antibody-l-w50a/legs/complex/rep01/colvar.dat"
ts, cv = [], []
for line in open(colvar):
    if line.startswith("#!") or not line.strip():
        continue
    p = line.split()
    ts.append(float(p[0])); cv.append(float(p[1]))
import statistics
mx = max(cv); med = statistics.median(cv)
thresh = max(med + 3.0, 0.7 * mx)
print(f"cv: median={med:.1f} max={mx:.1f} -> hydrated threshold {thresh:.1f}")
# hydrated epochs separated by >= 2000 ps
picks = []
for t, v in zip(ts, cv):
    if t < 5000:  # skip metad equilibration phase (bug fix 2026-09-01)
        continue
    if v >= thresh and all(abs(t - p) > 2000 for p in picks):
        picks.append(t)
    if len(picks) == 2:
        break
json.dump({"threshold": thresh, "median": med, "max": mx, "picks_ps": picks},
          open("runs/manual/opes_w50a_frame_picks.json", "w"))
print("picks:", picks)
PYEOF

picks=$("${ROOT}/.venv/bin/python" -c "import json; print(' '.join(str(p) for p in json.load(open('runs/manual/opes_w50a_frame_picks.json'))['picks_ps']))")
if [ -z "${picks}" ]; then
  log "no hydrated frames found (cv never crossed threshold) -- WATER DID NOT ENTER; stop here"
  exit 0
fi

# inject into rep02 / rep03 (complex leg)
i=2
for t in ${picks}; do
  target="${JOB}/legs/complex/rep0${i}/equilibration"
  [ -d "${target}" ] || { log "no rep0${i}, skip"; continue; }
  printf 'System\nSystem\n' | "${GMX}" trjconv -f "${REP}/opes.xtc" -s "${REP}/opes.tpr" \
    -o "${target}/npt_seed_hydrated.gro" -dump "${t}" -pbc mol -ur compact > /dev/null 2>&1
  if [ -s "${target}/npt_seed_hydrated.gro" ]; then
    cp "${target}/npt.gro" "${target}/npt.gro.pre_hydrated" 2>/dev/null || true
    cp "${target}/npt_seed_hydrated.gro" "${target}/npt.gro"
    log "injected hydrated frame ${t} ps -> complex/rep0${i}"
  fi
  i=$((i+1))
done

# wipe sample downstream only (lambda dirs kept -> no build_legs wipe needed)
# full window-output wipe: recovery marks sample complete if dhdl.xvg remain
find "${JOB}/legs" -path "*lambda*" \( -name dhdl.xvg -o -name md.gro -o -name md.log -o -name topol.tpr \) -delete
find "${JOB}/legs" -maxdepth 3 -name bar -type d -exec rm -rf {} +
rm -f "${JOB}"/stages/{sample,bar,qc,report}.json "${JOB}"/results/*.json
log "resuming job with hydrated ensemble"
ABAG_RBFE_VISIBLE_GPUS=1 "${ABAG_RBFE}" resume 3be1-antibody-l-w50a --batch-dir "${BATCH}" --execute \
  > "${ROOT}/runs/manual/opes_w50a_refep.log" 2>&1
if [ -f "${JOB}/results/ddg_summary.json" ]; then
  new=$(python3 -c "import json; print(json.load(open('${JOB}/results/ddg_summary.json'))['ddg_kcal_mol'])")
  log "NEW ddg=${new} (old 9.73, exp 1.40)"
else
  log "refep did not produce ddg_summary"
fi
log "stage2 done"
