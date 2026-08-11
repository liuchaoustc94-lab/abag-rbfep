#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"

"${ABAG_RBFE}" batch materialize-abbind --benchmark-root "${BENCHMARK_ROOT}"

cat "${BENCHMARK_ROOT}/materialized/summary.json"
