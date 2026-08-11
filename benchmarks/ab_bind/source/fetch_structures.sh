#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
TMPDIR="${ROOT}/tmp/ab_bind_structures_sync"
OUTDIR="${ROOT}/benchmarks/ab_bind/source/structures"

rm -rf "${TMPDIR}"
git clone --depth 1 https://github.com/sarahsirin/AB-Bind-Database "${TMPDIR}" >/tmp/ab_bind_structures_clone.log 2>&1

mkdir -p "${OUTDIR}"
find "${TMPDIR}" -maxdepth 1 -type f -name '*.pdb' -exec cp {} "${OUTDIR}/" \;

echo "Fetched structures into ${OUTDIR}"
