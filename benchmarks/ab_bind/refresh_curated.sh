#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
SOURCE="${ROOT}/benchmarks/ab_bind/source/AB-Bind_experimental_data.csv"
ANNOTATIONS="${ROOT}/benchmarks/ab_bind/source/ab_bind_complex_annotations.csv"
OUTDIR="${ROOT}/benchmarks/ab_bind"

"${ABAG_RBFE}" batch curate-abbind \
  --source-csv "${SOURCE}" \
  --annotations "${ANNOTATIONS}" \
  --output-dir "${OUTDIR}"

cat "${OUTDIR}/summary.json"
