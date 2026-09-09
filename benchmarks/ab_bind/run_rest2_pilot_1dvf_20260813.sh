#!/usr/bin/env bash
# REST2-pilot: endpoint ensemble via internal REST2 (hot region = mutation site)
# for the two worst aromatic-deletion jobs (1DVF y102a / y49a).
# Same evaluation question as the RID pilot: does engine-sampled endpoint
# diversity move the FEP ddG mean toward experiment (4.79 / 1.90)?
set -uo pipefail
REPO="/mnt/data/liuchao/platform/gromacs-abag-mmgbsa"
RBFEP="/mnt/data/liuchao/abag-rbfep"
ABAG_RBFE="${RBFEP}/.venv/bin/abag-rbfe"
GMX="${REPO}/tools/gromacs-gpu/bin/gmx"
export PATH="${RBFEP}/tmp/rest2_venv/bin:${RBFEP}/.venv/bin:${PATH}"

run_one() {
  local seed_chain="$1" seed_res="$2" job_suffix="$3"
  local tag="${seed_chain}${seed_res}"
  local rest2_dir="${REPO}/runs/1dvf_rest2_seed_${tag}_20260813"
  echo "[rest2-pilot] $(date --iso-8601=seconds) REST2 seed ${seed_chain}:${seed_res}"
  if [ ! -d "${rest2_dir}" ]; then
    ( cd "${REPO}" && bash run_abag_rest2.sh \
        -i "${RBFEP}/benchmarks/ab_bind/source/structures/1DVF.pdb" \
        --seed-residues "${seed_chain}:${seed_res}" \
        --antibody-chains A,B --antigen-chains C,D \
        --replicas 3 --scales 1.0,0.94,0.88 --gpu-ids 0,1,2 \
        --run-dir "${rest2_dir}" -s prepare )
  fi
  if [ ! -f "${rest2_dir}/rest2/rep00/rest2_prod.xtc" ]; then
    ( cd "${REPO}" && bash run_abag_rest2.sh \
        -i "${RBFEP}/benchmarks/ab_bind/source/structures/1DVF.pdb" \
        --seed-residues "${seed_chain}:${seed_res}" \
        --antibody-chains A,B --antigen-chains C,D \
        --replicas 3 --scales 1.0,0.94,0.88 --gpu-ids 0,1,2 \
        --run-dir "${rest2_dir}" -s rest2md )
  fi
  # 从 rep00 轨迹取 3 个远隔帧作为 FEP seed（无偏物理 replica）
  local prod="${rest2_dir}/rest2/rep00"
  [ -f "${prod}/rest2_prod.xtc" ] || { echo "[rest2-pilot] no rep00 xtc for ${tag}"; return 1; }
  local frames=(1000 3000 5000)
  local i=1
  local seedpdb_dir="${RBFEP}/runs/real_cases/rest2_seeds_1dvf_20260813"
  mkdir -p "${seedpdb_dir}"
  for fps in "${frames[@]}"; do
    local out="${seedpdb_dir}/${tag}_frame${fps}.pdb"
    printf 'Protein\nProtein\nProtein\n' | "${GMX}" trjconv -f "${prod}/rest2_prod.xtc" -s "${prod}/rest2_prod.tpr" \
      -o "${out}" -dump "${fps}" -pbc mol -ur compact -center >/dev/null 2>&1 || true
    [ -s "${out}" ] && echo "[rest2-pilot] seed pdb ${out} ($(grep -c '^ATOM' "${out}") atoms)"
    i=$((i+1))
  done
  # 每个 seed 构象一个 FEP job（对应突变）
  cat > /tmp/rest2_mut_${tag}.csv <<CSV
mutation_group_id,chain_id,resseq,icode,wt,mut,entity_side
${job_suffix},${seed_chain},${seed_res},,Y,A,$( [ "${seed_chain}" = "A" ] || [ "${seed_chain}" = "B" ] && echo antibody || echo antigen )
CSV
  for pdb in "${seedpdb_dir}/${tag}_frame"*.pdb; do
    [ -f "${pdb}" ] || continue
    local ftag=$(basename "${pdb}" .pdb)
    cat > /tmp/sys_${ftag}.yml <<EOF
system_name: 1dvf_rest2_${ftag}
input_structure: ${pdb}
structure_source: experimental
antibody_chains: [A, B]
antigen_chains: [C, D]
EOF
    local batch="1dvf_rest2_${ftag}_20260813"
    if [ ! -d "${RBFEP}/runs/real_cases/${batch}" ]; then
      "${ABAG_RBFE}" batch plan --system "/tmp/sys_${ftag}.yml" \
        --mutations "/tmp/rest2_mut_${tag}.csv" \
        --protocol "${RBFEP}/benchmarks/ab_bind/protocol.validation_priority.yml" \
        --batch-id "${batch}" --runs-root "${RBFEP}/runs/real_cases" >/dev/null
    fi
    for jd in "${RBFEP}/runs/real_cases/${batch}/jobs"/*; do
      "${ABAG_RBFE}" resume "$(basename ${jd})" --batch-dir "${RBFEP}/runs/real_cases/${batch}" --execute 2>&1 | tail -1
    done
  done
}

run_one D 102 antigen_d_y102a
run_one A 49 antibody_a_y49a
echo "[rest2-pilot] $(date --iso-8601=seconds) all done"
