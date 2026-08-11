#!/usr/bin/env bash
set -euo pipefail

ROOT="/mnt/data/liuchao/abag-rbfep"
OUT="${ROOT}/benchmarks/ab_bind/source/AB-Bind_experimental_data.csv"

curl -L \
  -o "${OUT}" \
  "https://raw.githubusercontent.com/sarahsirin/AB-Bind-Database/master/AB-Bind_experimental_data.csv"

echo "Fetched ${OUT}"
