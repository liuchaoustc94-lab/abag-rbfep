#!/usr/bin/env bash
set -euo pipefail

ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"
ABAG_RBFE="${ROOT}/.venv/bin/abag-rbfe"
BENCHMARK_ROOT="${ROOT}/benchmarks/ab_bind"
SPLIT_FILE="${BENCHMARK_ROOT}/splits/ab_bind_rbfe_core_v1_split_v1.yml"
PROTOCOL_PATH="${PROTOCOL_PATH:-${BENCHMARK_ROOT}/protocol.1mlc_target_specific_sampling_pilot.yml}"
RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_minibatch_1mlc_20260626}"
BATCH_ID="${BATCH_ID:-abbind-target-specific-sampling-minibatch-1mlc_1mlc_core_v1}"
COMPLEX_ID="${COMPLEX_ID:-1MLC}"
MAX_WORKERS="${MAX_WORKERS:-2}"
GPU_DEVICES="${GPU_DEVICES:-0,1}"
PLAN_ONLY="${PLAN_ONLY:-0}"
INPUT_ROOT="${RUNS_ROOT}/inputs/${COMPLEX_ID}"
SYSTEM_PATH="${BENCHMARK_ROOT}/materialized/${COMPLEX_ID}/system.yml"
SOURCE_MUTATIONS="${BENCHMARK_ROOT}/materialized/${COMPLEX_ID}/core_v1_mutations.csv"
FILTERED_MUTATIONS="${INPUT_ROOT}/core_v1_mutations.minibatch.csv"
JOB_MAP_CSV="${JOB_MAP_CSV:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_target_specific_sampling_pilot_1mlc_20260626/abbind-target-specific-sampling-pilot-1mlc_1mlc_core_v1/jobs.csv}"

DEFAULT_JOB_IDS=(
  "1mlc-antibody-l-n92a"
  "1mlc-antibody-l-n32g"
  "1mlc-antibody-l-n32y"
)

if [ "$#" -gt 0 ]; then
  JOB_IDS=("$@")
else
  JOB_IDS=("${DEFAULT_JOB_IDS[@]}")
fi

mkdir -p "${INPUT_ROOT}"

# Build a mutations CSV that only contains the requested mini-batch jobs.
"${ROOT}/.venv/bin/python" - <<'PY' "${SOURCE_MUTATIONS}" "${FILTERED_MUTATIONS}" "${JOB_MAP_CSV}" "${JOB_IDS[@]}"
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
dest = Path(sys.argv[2])
job_map_csv = Path(sys.argv[3])
requested_job_ids = sys.argv[4:]
fallback_map = {
    "1mlc-antibody-l-n32g": "1mlc_0319",
    "1mlc-antibody-l-n32y": "1mlc_0320",
    "1mlc-antibody-l-n92a": "1mlc_0321",
}

job_to_group = dict(fallback_map)
if job_map_csv.is_file():
    with job_map_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            job_id = row.get("job_id", "").strip()
            mutation_group_id = row.get("mutation_group_id", "").strip()
            if job_id and mutation_group_id:
                job_to_group[job_id] = mutation_group_id

requested_groups = []
missing_job_ids = []
for job_id in requested_job_ids:
    mutation_group_id = job_to_group.get(job_id)
    if mutation_group_id is None:
        missing_job_ids.append(job_id)
    else:
        requested_groups.append(mutation_group_id)

if missing_job_ids:
    raise SystemExit(f"Missing job_id to mutation_group_id mapping: {', '.join(sorted(missing_job_ids))}")

requested = set(requested_groups)

with source.open(newline="") as handle:
    reader = csv.DictReader(handle)
    rows = [row for row in reader if row.get("mutation_group_id", "").strip() in requested]
    fieldnames = list(reader.fieldnames or [])

missing = sorted(requested.difference(row["mutation_group_id"] for row in rows))
if missing:
    raise SystemExit(f"Missing requested mutation_group_id entries: {', '.join(missing)}")

dest.parent.mkdir(parents=True, exist_ok=True)
with dest.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY

"${ABAG_RBFE}" batch plan \
  --system "${SYSTEM_PATH}" \
  --mutations "${FILTERED_MUTATIONS}" \
  --protocol "${PROTOCOL_PATH}" \
  --batch-id "${BATCH_ID}" \
  --runs-root "${RUNS_ROOT}"

"${ROOT}/.venv/bin/python" - <<'PY' "${RUNS_ROOT}" "${BENCHMARK_ROOT}" "${PROTOCOL_PATH}" "${BATCH_ID}" "${SYSTEM_PATH}" "${FILTERED_MUTATIONS}" "${COMPLEX_ID}"
import json
import sys
from pathlib import Path

try:
    import yaml
except ImportError as exc:
    raise SystemExit(f"PyYAML is required to materialize plan_index.yml: {exc}")

runs_root = Path(sys.argv[1])
benchmark_root = Path(sys.argv[2])
protocol_path = Path(sys.argv[3])
batch_id = sys.argv[4]
system_path = Path(sys.argv[5])
mutations_path = Path(sys.argv[6])
complex_id = sys.argv[7]

batch_dir = runs_root / batch_id
jobs_csv = batch_dir / "jobs.csv"
job_count = 0
if jobs_csv.is_file():
    job_count = max(sum(1 for _ in jobs_csv.open()) - 1, 0)

payload = {
    "benchmark_root": str(benchmark_root),
    "spec_name": "core_v1",
    "protocol_path": str(protocol_path),
    "plan_root": str(runs_root),
    "split_name": "",
    "split_path": "",
    "planned_batch_count": 1,
    "planned_complexes": [complex_id],
    "batches": [
        {
            "complex_id": complex_id,
            "batch_id": batch_id,
            "batch_dir": str(batch_dir),
            "system_yml": str(system_path),
            "mutations_csv": str(mutations_path),
            "job_count": job_count,
            "mutation_group_count": job_count,
            "structure_source": "experimental",
            "antibody_chains": "HL",
            "antigen_chains": "E",
        }
    ],
}

(runs_root / "plan_index.json").write_text(json.dumps(payload, indent=2) + "\n")
(runs_root / "plan_index.yml").write_text(yaml.safe_dump(payload, sort_keys=False))
PY

if [ "${PLAN_ONLY}" = "1" ]; then
  echo "Planned mini-batch only. Inspect ${RUNS_ROOT}/${BATCH_ID}"
  exit 0
fi

RUN_ARGS=(
  batch run-abbind
  --plan-root "${RUNS_ROOT}"
  --batch-id "${BATCH_ID}"
  --complex-id "${COMPLEX_ID}"
  --resume
  --execute
  --max-workers "${MAX_WORKERS}"
)

if [ -n "${GPU_DEVICES}" ]; then
  RUN_ARGS+=(--gpu-devices "${GPU_DEVICES}")
fi

for job_id in "${JOB_IDS[@]}"; do
  RUN_ARGS+=(--job-id "${job_id}")
done

"${ABAG_RBFE}" "${RUN_ARGS[@]}"

"${ABAG_RBFE}" batch report-abbind \
  --plan-root "${RUNS_ROOT}" \
  --batch-id "${BATCH_ID}" \
  --complex-id "${COMPLEX_ID}" \
  --split-name validation \
  --split-file "${SPLIT_FILE}"

"${ROOT}/.venv/bin/python" "${BENCHMARK_ROOT}/report_hotspot_root_comparison.py" \
  --job-id 1mlc-antibody-l-n92a \
  --job-id 1mlc-antibody-l-n32g \
  --job-id 1mlc-antibody-l-n32y \
  --plan-root "priority=${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan" \
  --plan-root "robust=${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan" \
  --plan-root "sampling_qc=${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues" \
  --plan-root "deep=${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues" \
  --plan-root "ultra=${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues" \
  --plan-root "minibatch_1mlc=${RUNS_ROOT}" \
  --json-output "${RUNS_ROOT}/reports/hotspot_root_comparison.json" \
  --md-output "${RUNS_ROOT}/reports/hotspot_root_comparison.md" \
  || true

printf '%s\n' "${JOB_IDS[@]}" > "${RUNS_ROOT}/reports/minibatch_job_ids.txt"

echo "Finished. Inspect ${RUNS_ROOT}/reports/plan_summary.json"
