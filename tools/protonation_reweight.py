#!/usr/bin/env python
"""P2-b: protonation multi-state (Wyman linkage) ddG reweighting.

Our FEP models titratable residues in a single protonation state (His neutral,
Asp/Glu charged, Lys charged). The physical WT ensemble contains the other
state with weight set by (pKa - pH). For a mutation that REMOVES a titratable
site (His->Ala etc.), the measured ddG differs from the single-state FEP ddG by

    ddG_true = ddG_FEP + RT * ln( f_complex / f_apo )

    f = 1 + 10^(pKa - pH)   for bases (His, Lys: protonation adds charge)
    f = 1 + 10^(pH - pKa)   for acids (Asp, Glu: deprotonation adds charge)

For a mutation CREATING a titratable site (Ala->His) the sign flips.
pKa values come from PROPKA on the complex and on the unbound entity that
carries the site (antibody-only / antigen-only chain-filtered PDB).

Diagnostic only (opt-in): corrections never enter official metrics by default
(project decision: Sampson-style corrections trade R for MAE; report both).

Usage:
  tmp/propka_venv/bin/python tools/protonation_reweight.py \
      --out runs/analysis/protonation_reweight_20260831 [--ph 7.4]
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

R_KCAL = 0.0019872  # kcal/(mol K)
TITRATABLE = {"H": "HIS", "D": "ASP", "E": "GLU", "K": "LYS", "R": "ARG"}
BASES = {"HIS", "LYS", "ARG"}
GOOD_BATCH_HINTS = (
    "abbind_validation_panels_20260812", "abbind_fit_newprotocol_20260811",
    "runs/acceptance/", "runs/real_cases/oos_", "3hfm_patel_align", "e6_charmm",
    "runs/real_cases/1dvf_", "abbind_core_v1_validation_priority_plan", "abbind_core_v1_validation_plan",
    "3hfm_patel_dssb", "dssb_smoke",
)


def collect_pairs(root: Path) -> list[dict]:
    """Best-available (deepest sampling) ddg per (complex, mutation)."""
    exp: dict[tuple[str, str], float] = {}
    # full registered set first (covers SKEMPI out-of-sample targets like
    # 1AK4/1KTZ/1MLC/1DVF/3K2M that the core subsets filter out)
    for r in csv.DictReader(open(root / "benchmarks/ab_bind/curated/ab_bind_source_registered.csv")):
        m = re.search(r"([A-Z]+):([A-Z])(\d+)([A-Z])", r["mutation_tokens"])
        if m and r["ddg_kcal_mol"] not in ("", "NA"):
            exp[(r["complex_id"].lower(), f"{m.group(2)}{m.group(3)}{m.group(4)}".lower())] = float(r["ddg_kcal_mol"])
    for f in ("benchmarks/ab_bind/curated/ab_bind_rbfe_core_v1.csv",
              "benchmarks/ab_bind/curated/ab_bind_rbfe_core_v2.csv"):
        for r in csv.DictReader(open(root / f)):
            m = re.match(r"([A-Z]+):([A-Z])(\d+)([A-Z])@", r["mutation_tokens"])
            if m:
                exp[(r["complex_id"].lower(), f"{m.group(2)}{m.group(3)}{m.group(4)}".lower())] = float(r["ddg_kcal_mol"])
    panels = json.load(open(root / "benchmarks/acceptance/skedmi_mutation_panels.json"))
    for tgt, pdata in panels.items():
        for m in pdata["muts"]:
            exp[(tgt.lower(), f"{m['wt']}{m['resseq']}{m['mut']}".lower())] = m["ddg"]
    # Patel 2021 3HFM reference set (DSSB charge-changing jobs)
    patel = root / "benchmarks/patel_2021_3hfm/experimental_ddg.csv"
    if patel.is_file():
        for r in csv.DictReader(open(patel)):
            vals = list(r.values())
            # columns: job_id, group, chain, resseq, wt, mut, side, ddg, ...
            exp[("3hfm", f"{vals[4]}{vals[3]}{vals[5]}".lower())] = float(vals[7])

    cands: dict[tuple[str, str], tuple[int, float, dict, str]] = {}
    for f in glob.glob(str(root / "runs/**/jobs/*/results/ddg_summary.json"), recursive=True):
        if not any(h in f for h in GOOD_BATCH_HINTS):
            continue
        job_dir = Path(f).parent.parent
        job_id = job_dir.name
        m = re.search(r"-([a-z]+)-([a-z])(\d+)([a-z])$", job_id)
        if not m:
            continue
        mut = f"{m.group(2)}{m.group(3)}{m.group(4)}".lower()
        cx = None
        tokens = job_id.split("-")
        for c in {k[0] for k in exp}:
            if tokens[0] == c or (tokens[0] == "acc" and len(tokens) > 1 and tokens[1] == c):
                cx = c
                break
        if not cx or (cx, mut) not in exp:
            continue
        d = json.load(open(f))
        pred = d.get("ddg_kcal_mol")
        if pred is None:
            continue
        spec = json.load(open(job_dir / "job_spec.json"))
        prot = spec["protocol"]
        score = int(prot.get("lambda_windows", 0)) * int(prot.get("repeats", 0))
        key = (cx, mut)
        if key not in cands or score > cands[key][0]:
            cands[key] = (score, float(pred), spec, str(job_dir))
    return [
        {"complex": cx, "mut": mut, "pred": v[1], "exp": exp[(cx, mut)],
         "job_dir": v[3], "sites": v[2]["mutation_group"]["sites"]}
        for (cx, mut), v in sorted(cands.items())
    ]


def build_entity_pdb(complex_pdb: Path, chains: set[str], out: Path) -> Path:
    if out.is_file():
        return out
    lines = []
    for l in open(complex_pdb):
        if not l.startswith(("ATOM", "HETATM")) or l[21] not in chains:
            continue
        # propka 3.5 chokes on explicit hydrogens (e.g. N-terminal H1/H2);
        # strip them — propka builds its own protonation model.
        if l[76:78].strip() == "H" or l[12:16].strip().startswith(("H", "1H", "2H", "3H")):
            continue
        lines.append(l)
    out.write_text("".join(lines) + "END\n", encoding="utf-8")
    return out


def run_propka(pdb: Path, out_dir: Path, propka_python: Path) -> Path:
    """Run propka 3.5 CLI (writes <input_stem>.pka next to the input pdb)."""
    pka = out_dir / (pdb.stem + ".pka")
    if pka.is_file():
        return pka
    pdb_abs = pdb.resolve()
    subprocess.run([str(propka_python), "-m", "propka", str(pdb_abs)],
                   check=True, capture_output=True, cwd=out_dir)
    produced = pdb_abs.with_suffix(".pka")
    if produced.is_file() and produced != pka:
        produced.rename(pka)
    if not pka.is_file():
        raise RuntimeError(f"propka failed for {pdb}")
    return pka


def parse_pka(pka_path: Path) -> dict[tuple[str, int, str], float]:
    """(resname, resnum, chain) -> pKa from the SUMMARY table."""
    table = {}
    in_summary = False
    for line in open(pka_path):
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if in_summary:
            m = re.match(r"\s+([A-Z]{3})\s+(\d+)\s+(\w)\s+([\d.]+)", line)
            if m:
                table[(m.group(1), int(m.group(2)), m.group(3))] = float(m.group(4))
            elif line.strip().startswith("---") is False and table and not line.strip():
                break
            if "Free energy" in line or "------------" in line and table:
                break
    return table


def factor(resname: str, pka: float, ph: float) -> float:
    if resname in BASES:
        return 1.0 + 10 ** (pka - ph)
    return 1.0 + 10 ** (ph - pka)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--temp", type=float, default=310.0)
    ap.add_argument("--propka-python", default="tmp/propka_venv/bin/python")
    args = ap.parse_args()
    root = Path(".").resolve()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    propka_python = root / args.propka_python
    rt = R_KCAL * args.temp

    pairs = collect_pairs(root)
    print(f"collected {len(pairs)} best-available pairs")

    # per-complex propka: complex + both entities
    pka_cache: dict[str, dict[str, dict]] = {}
    results = []
    n_corr = 0
    for pair in pairs:
        cx = pair["complex"]
        if cx not in pka_cache:
            sys_yml = None
            for cand in glob.glob(str(root / f"benchmarks/ab_bind/materialized/{cx.upper()}/system.yml")) + \
                        glob.glob(str(root / f"benchmarks/acceptance/{cx}/system.yml")):
                sys_yml = cand
                break
            if not sys_yml:
                # fall back to the job's own config
                jd = pair["job_dir"]
                sys_yml = str(Path(jd) / "config" / "system.yml")
            import yaml
            sysconf = yaml.safe_load(open(sys_yml))
            complex_pdb = Path(sysconf["input_structure"])
            ab_chains = set(sysconf["antibody_chains"])
            ag_chains = set(sysconf["antigen_chains"])
            cdir = out / cx
            cdir.mkdir(exist_ok=True)
            local_complex = build_entity_pdb(complex_pdb, ab_chains | ag_chains, cdir / "complex.pdb")
            ab_pdb = build_entity_pdb(complex_pdb, ab_chains, cdir / "apo_antibody.pdb")
            ag_pdb = build_entity_pdb(complex_pdb, ag_chains, cdir / "apo_antigen.pdb")
            pka_cache[cx] = {
                "complex": parse_pka(run_propka(local_complex, cdir, propka_python)),
                "antibody": parse_pka(run_propka(ab_pdb, cdir, propka_python)),
                "antigen": parse_pka(run_propka(ag_pdb, cdir, propka_python)),
            }

        site = pair["sites"][0]
        wt, mut = str(site["wt"]).upper(), str(site["mut"]).upper()
        chain, resseq = site["chain_id"], int(site["resseq"])
        correction = 0.0
        note = ""
        if wt in TITRATABLE or mut in TITRATABLE:
            if wt in TITRATABLE and mut in TITRATABLE:
                # titratable -> titratable (e.g. K->R): the site is not removed,
                # the linkage factors nearly cancel; v1 skips (charge-changing
                # swaps like K->D belong to the DSSB path anyway).
                note = "titratable-to-titratable swap: no correction in v1"
                results.append({**pair, "protonation_correction": 0.0, "note": note})
                continue
            # the titratable residue exists in WT (wt side) or in MUT (mut side)
            side = "antibody" if site.get("entity_side") == "antibody" else "antigen"
            resname_wt = TITRATABLE.get(wt)
            resname_mut = TITRATABLE.get(mut)
            if resname_wt:
                pka_cx = pka_cache[cx]["complex"].get((resname_wt, resseq, chain))
                pka_apo = pka_cache[cx][side].get((resname_wt, resseq, chain))
                sign = +1.0
                resname = resname_wt
            else:
                # mutating INTO a titratable residue: pKa of the mutant residue
                # is not measurable on the WT structure -> flag only
                resname = resname_mut
                sign = -1.0
                pka_cx = pka_apo = None
            if pka_cx is not None and pka_apo is not None:
                correction = sign * rt * math.log(factor(resname, pka_cx, args.ph) / factor(resname, pka_apo, args.ph))
                note = f"{resname}{chain}{resseq} pKa cx={pka_cx:.2f} apo={pka_apo:.2f}"
                n_corr += 1
            else:
                note = f"{resname}{chain}{resseq} pKa lookup miss (resnumbering?)"
        results.append({**pair, "protonation_correction": correction, "note": note})

    # metrics before/after
    def metrics(preds, exps):
        n = len(preds)
        mp, me = sum(preds) / n, sum(exps) / n
        cov = sum((a - mp) * (b - me) for a, b in zip(preds, exps))
        vp = sum((a - mp) ** 2 for a in preds)
        ve = sum((b - me) ** 2 for b in exps)
        r = cov / math.sqrt(vp * ve) if vp > 0 and ve > 0 else float("nan")
        mae = sum(abs(a - b) for a, b in zip(preds, exps)) / n
        return r, mae

    pred_raw = [p["pred"] for p in results]
    pred_corr = [p["pred"] + p["protonation_correction"] for p in results]
    exps = [p["exp"] for p in results]
    r0, mae0 = metrics(pred_raw, exps)
    r1, mae1 = metrics(pred_corr, exps)

    with open(out / "pairs_protonation.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["complex", "mut", "pred", "exp", "protonation_correction", "note", "job_dir"])
        w.writeheader()
        for p in results:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in p.items() if k in w.fieldnames})

    corrected = [p for p in results if abs(p["protonation_correction"]) > 0.01]
    summary = {
        "n_pairs": len(results),
        "n_titratable_corrected": n_corr,
        "n_correction_gt_0.01": len(corrected),
        "ph": args.ph,
        "raw": {"pearson": round(r0, 4), "mae": round(mae0, 4)},
        "corrected": {"pearson": round(r1, 4), "mae": round(mae1, 4)},
        "top_corrections": sorted(
            [{"complex": p["complex"], "mut": p["mut"], "corr": round(p["protonation_correction"], 3),
              "pred": round(p["pred"], 2), "exp": p["exp"], "note": p["note"]} for p in corrected],
            key=lambda d: -abs(d["corr"]))[:15],
    }
    json.dump(summary, open(out / "summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1, ensure_ascii=False)[:2000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
