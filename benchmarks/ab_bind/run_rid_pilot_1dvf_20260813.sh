#!/usr/bin/env bash
# RID pilot: endpoint-ensemble exploration for 1DVF -> FEP re-run from basins.
# Decision context: basin-diversity test showed per-basin ddG spread 2.7-3.3
# kcal/mol for y102a/y49a (worst aromatic deletions). This pilot tests whether
# RID-sampled basins shift the MEAN toward experiment (4.79 / 1.90).
#
# Boundary discipline: rid_kit only produces conformations; abag-rbfep only
# consumes protein coordinates (fresh FEP jobs from basin PDBs, standard stages).
set -uo pipefail
RID_REPO="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa"
RBFEP="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${RBFEP}/.venv/bin/abag-rbfe"
GMX="${RID_REPO}/tools/gromacs-gpu/bin/gmx"
STAMP="20260813"
RID_OUT="${RID_REPO}/rid_kit/runs/1dvf_rid_pilot_${STAMP}"
export PATH="${RBFEP}/.venv/bin:${PATH}"

# ---------- Stage 1: RID pilot on WT 1DVF complex ----------
if [ ! -f "${RID_OUT}/results/basin_manifest.tsv" ]; then
  echo "[rid-pilot] $(date --iso-8601=seconds) starting RID pilot (1DVF, A,B|C,D)"
  cd "${RID_REPO}"
  bash run_abag_rid.sh -i "${RBFEP}/benchmarks/ab_bind/source/structures/1DVF.pdb" \
    --antibody-chains A,B --antigen-chains C,D \
    --profile pilot --gpu-ids 0,1,2 --run-dir "${RID_OUT}"
fi
echo "[rid-pilot] basins:"
cat "${RID_OUT}/results/basin_manifest.tsv" 2>/dev/null || { echo "no basins produced"; exit 1; }

# ---------- Stage 2: basin -> protein-only PDB ----------
BASIN_DIR="${RBFEP}/runs/real_cases/rid_basins_1dvf_${STAMP}"
mkdir -p "${BASIN_DIR}"
for gro in "${RID_OUT}"/unbiased/basin_*/rep_*/md/unbiased_npt.gro; do
  [ -f "${gro}" ] || continue
  tag=$(echo "${gro}" | sed 's/.*\(basin_[0-9]*\).*/\1/')
  out="${BASIN_DIR}/${tag}.pdb"
  if [ ! -f "${out}" ]; then
    echo "Protein" | "${GMX}" trjconv -f "${gro}" -s "${gro}" -o "${out}" -pbc mol 2>/dev/null || \
    echo "Protein" | "${GMX}" editconf -f "${gro}" -o "${out}" 2>/dev/null
  fi
  echo "[rid-pilot] basin pdb: ${out} ($(grep -c '^ATOM' "${out}" 2>/dev/null || echo 0) atoms)"
done

# ---------- Stage 3: per-basin FEP jobs (y102a + y49a), standard pipeline ----------
cat > /tmp/1dvf_rid_mutations.csv <<'CSV'
mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side
antigen_d_y102a,D,102,,Y,A,antigen
antibody_a_y49a,A,49,,Y,A,antibody
CSV
for pdb in "${BASIN_DIR}"/basin_*.pdb; do
  [ -f "${pdb}" ] || continue
  tag=$(basename "${pdb}" .pdb)
  sysyml="${BASIN_DIR}/system_${tag}.yml"
  cat > "${sysyml}" <<EOF
system_name: 1dvf_rid_${tag}
input_structure: ${pdb}
structure_source: experimental
antibody_chains: [A, B]
antigen_chains: [C, D]
notes:
  - 1DVF conformation from RID basin ${tag} (unbiased-relaxed), pilot ${STAMP}.
EOF
  batch="1dvf_rid_${tag}_${STAMP}"
  if [ ! -d "${RBFEP}/runs/real_cases/${batch}" ]; then
    "${ABAG_RBFE}" batch plan \
      --system "${sysyml}" \
      --mutations /tmp/1dvf_rid_mutations.csv \
      --protocol "${RBFEP}/benchmarks/ab_bind/protocol.validation_priority.yml" \
      --batch-id "${batch}" \
      --runs-root "${RBFEP}/runs/real_cases"
  fi
  for job in "1dvf-rid-${tag}-antigen-d-y102a" "1dvf-rid-${tag}-antibody-a-y49a"; do
    echo "[rid-pilot] FEP ${job}"
    "${ABAG_RBFE}" resume "${job}" --batch-dir "${RBFEP}/runs/real_cases/${batch}" --execute 2>&1 | tail -1
  done
done
echo "[rid-pilot] $(date --iso-8601=seconds) all done"
