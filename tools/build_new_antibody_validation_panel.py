#!/usr/bin/env python3
"""Build a conservative, provenance-tracked panel of new Ab--Ag systems.

The panel deliberately keeps experimental SKEMPI2 labels separate from the
Rosetta Flex ddG controls distributed with Graphinity.  This script only
materializes selected rows; the full upstream benchmark files are not copied
into the repository.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import date
import hashlib
import math
from pathlib import Path
import re
import shutil
import ssl
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmarks" / "antibody_new10"
SKEMPI_DEFAULT = Path("/tmp/skempi_new10/skempi_v2.csv")
SKEMPI_PDB_ROOT = Path("/tmp/skempi_new10/PDBs")
GRAPHINITY_DEFAULT = Path(
    "/tmp/Graphinity/data/ddg_synthetic/Flex_ddG/Synthetic_FlexddG_ddG_20829.csv"
)
GRAPHINITY_SEQUENCES_DEFAULT = Path(
    "/tmp/ab-ag-ddg/ddg_synthetic/Flex_ddG/Synthetic_FlexddG_ddG_20829_with_sequences.csv"
)

R_KCAL = 0.001987204258
TOKEN_RE = re.compile(
    r"^(?P<wt>[A-Z])(?P<chain>[A-Za-z0-9])(?P<resseq>-?\d+)(?P<icode>[A-Za-z]?)(?P<mut>[A-Z])$"
)

# This is intentionally conservative: a target is excluded if it has already
# appeared in a formal report, acceptance queue, historical run, or local
# materialized/queued benchmark assets.
TRIED_TARGETS = {
    "1AK4",
    "1AHW",
    "1BJ1",
    "1CZ8",
    "1DQJ",
    "1DVF",
    "1IAR",
    "1JRH",
    "1KTZ",
    "1MLC",
    "1MHP",
    "1N8Z",
    "1NMB",
    "1OGA",
    "1VFB",
    "1YY9",
    "2BDN",
    "2JEL",
    "2NYY",
    "2NRZ",  # harmless guard for a historical spelling typo
    "2NY7",
    "2NZ9",
    "3BDY",
    "3BE1",
    "3BN9",
    "3HFM",
    "3K2M",
    "3NGB",
    "3NPS",
    "4DN4",
    "4I77",
}


@dataclass(frozen=True)
class TargetSpec:
    pdb: str
    name: str
    antibody_chains: tuple[str, ...]
    antigen_chains: tuple[str, ...]
    label_type: str
    source_dataset: str
    notes: str = ""
    related_to_existing: str = ""
    allow_below_minimum: bool = False
    complex_class: str = "antibody-antigen"
    mutated_partner: str = "antibody"


TARGETS = (
    TargetSpec(
        "3SE8",
        "VRC03 Fab–HIV-1 gp120",
        ("H", "L"),
        ("G",),
        "experimental",
        "SKEMPI2",
        "28 antibody-side single substitutions; strongest strict new experimental target.",
    ),
    TargetSpec(
        "3SE9",
        "VRC-PG04–HIV-1 gp120",
        ("H", "L"),
        ("G",),
        "experimental",
        "SKEMPI2",
        "25 antibody-side single substitutions; independent antibody on the same antigen family.",
    ),
    TargetSpec(
        "3N85",
        "Fab37–HER2",
        ("H", "L"),
        ("A",),
        "experimental",
        "SKEMPI2",
        "10 source rows, but WH100F is duplicated; only 9 unique single signatures.",
    ),
    TargetSpec(
        "3G6D",
        "CNTO607 Fab–IL-13",
        ("H", "L"),
        ("A",),
        "experimental",
        "SKEMPI2",
        "10 antibody-side rows, dominated by same-side combination mutations; exploratory only.",
    ),
    TargetSpec(
        "3L5X",
        "H2L6 Fab–IL-13",
        ("H", "L"),
        ("A",),
        "experimental",
        "SKEMPI2",
        "36 antibody-side rows, mostly framework-adaptation combinations; V2/multi-point stress test.",
    ),
    TargetSpec(
        "2B2X",
        "AQC2 Fab–integrin alpha-1 (alternate structure)",
        ("H", "L"),
        ("A",),
        "experimental",
        "SKEMPI2",
        "115 antibody-side rows, but mostly combinations; structural replicate of the already queued 1MHP biology. The deposited 2B2X complex contains an affinity-matured Fab, so revert to a WT AQC2 baseline before RBFE.",
        "1MHP",
    ),
    TargetSpec(
        "3J1S",
        "A20 Fab–AAV2",
        ("H", "L"),
        ("A",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "15 Rosetta Flex ddG antibody substitutions (one of the 16 source rows is antigen-side); do not pool with experimental R metrics.",
    ),
    TargetSpec(
        "4KUC",
        "6C2 antibody–ricin A chain",
        ("H", "D"),
        ("A",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "16 Rosetta Flex ddG antibody substitutions; useful for pipeline/structure stress testing.",
    ),
    TargetSpec(
        "5JHL",
        "flavivirus broadly protective antibody–Zika E",
        ("H", "L"),
        ("A",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "16 Rosetta Flex ddG antibody substitutions; do not treat labels as experimental affinity.",
    ),
    TargetSpec(
        "4NNP",
        "antagonistic Fab–MntC",
        ("H", "L"),
        ("A",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "16 Rosetta Flex ddG antibody substitutions on a bacterial transporter target.",
    ),
)

# A single-point-focused companion panel.  The public experimental pool has
# only two clearly independent new targets with >=10 antibody-side singles;
# the remaining slots are explicitly labelled Rosetta controls rather than
# being filled with combination mutants.
SINGLE_POINT_TARGETS = (
    TARGETS[0],
    TARGETS[1],
    TARGETS[2],
    TARGETS[6],
    TARGETS[7],
    TARGETS[8],
    TARGETS[9],
    TargetSpec(
        "4PLJ",
        "8G12 antibody–hepatitis E virus E2",
        ("D", "C"),
        ("B",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "16 Rosetta Flex ddG antibody-side single substitutions.",
    ),
    TargetSpec(
        "6A79",
        "B5209B scFv–Robo1 Ig5",
        ("H", "L"),
        ("A",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "16 Rosetta Flex ddG antibody-side single substitutions; source paper reports only a smaller experimental mutant panel.",
    ),
    TargetSpec(
        "5J56",
        "V1C7 single-domain antibody–ricin A chain",
        ("B",),
        ("A",),
        "synthetic_control",
        "Graphinity-Flex-ddG",
        "16 Rosetta Flex ddG antibody-side single substitutions.",
    ),
)

# The SKEMPI2-only materialized subset.  5C6T is retained as an explicitly
# labelled near-threshold reserve (eight unique antibody-side singles); it is
# useful for source curation but is not formal-ready under the ten-point rule.
SKEMPI_SINGLE_POINT_TARGETS = (
    TARGETS[0],
    TARGETS[1],
    TARGETS[2],
    TargetSpec(
        "5C6T",
        "1G2 Fab–HCMV glycoprotein B",
        ("H", "L"),
        ("A",),
        "experimental",
        "SKEMPI2",
        "Eight unique antibody-side singles; near-threshold reserve.",
        allow_below_minimum=True,
    ),
)

# Experimental protein-protein extension.  The current CLI and mutation
# schema use ``antibody``/``antigen`` as the two-side compatibility names; for
# these entries ``antibody_chains`` means the mutated partner-A chains and the
# system metadata records the true protein-protein class.
PROTEIN_PROTEIN_TARGETS = (
    TargetSpec(
        "1C1Y",
        "Rap1A–Raf RBD",
        ("B",),
        ("A",),
        "experimental",
        "SKEMPI2",
        "13 unique single substitutions on the Raf RBD chain.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1A4Y",
        "ribonuclease inhibitor–angiogenin",
        ("A",),
        ("B",),
        "experimental",
        "SKEMPI2",
        "16 unique single substitutions on ribonuclease inhibitor chain A.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1C4Z",
        "E6AP–UBCH7",
        ("A", "B", "C"),
        ("D",),
        "experimental",
        "SKEMPI2",
        "22 unique single substitutions on E6AP copies A–C.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1DAN",
        "factor VIIa–tissue factor",
        ("T", "U"),
        ("H", "L"),
        "experimental",
        "SKEMPI2",
        "97 single records / 85 unique signatures on tissue-factor copies T/U.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1E50",
        "AML1 Runt domain–CBFβ",
        ("A", "C", "E", "G", "Q", "R"),
        ("B", "D", "F", "H"),
        "experimental",
        "SKEMPI2",
        "15 single records / 11 unique signatures on AML1 copies.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1EMV",
        "colicin E9 DNase–immunity protein Im9",
        ("A",),
        ("B",),
        "experimental",
        "SKEMPI2",
        "34 unique single substitutions on immunity protein Im9.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1F47",
        "FtsZ fragment–ZipA",
        ("A",),
        ("B",),
        "experimental",
        "SKEMPI2",
        "12 unique single substitutions on the FtsZ fragment.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1JTD",
        "TEM-1 beta-lactamase–BLIP-II",
        ("B",),
        ("A",),
        "experimental",
        "SKEMPI2",
        "25 unique single substitutions on BLIP-II.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "3QHY",
        "class-A beta-lactamase–BLIP-II",
        ("B",),
        ("A",),
        "experimental",
        "SKEMPI2",
        "25 unique single substitutions on BLIP-II; independent complex from 1JTD.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
    TargetSpec(
        "1FFW",
        "CheY–CheA",
        ("A",),
        ("B",),
        "experimental",
        "SKEMPI2",
        "11 unique single substitutions on the CheY chain.",
        complex_class="protein-protein",
        mutated_partner="partner_a",
    ),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "abag-rbfep-new10/1.0"})
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except ImportError as exc:  # pragma: no cover - only used on minimal hosts
        raise RuntimeError("certifi is required for structure downloads") from exc
    with urlopen(request, timeout=60, context=context) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def _read_csv(path: Path, delimiter: str = ",", encoding: str = "utf-8") -> list[dict[str, str]]:
    with path.open("r", encoding=encoding, newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def _write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _parse_token(token: str) -> dict[str, object] | None:
    match = TOKEN_RE.match(token.strip())
    if not match:
        return None
    values = match.groupdict()
    return {
        "wt": values["wt"].upper(),
        "chain_id": values["chain"].upper(),
        "resseq": int(values["resseq"]),
        "icode": values["icode"].upper(),
        "mut": values["mut"].upper(),
    }


def _split_tokens(value: str) -> list[str]:
    return [token.strip() for token in value.split(",") if token.strip()]


def _canonical_signature(sites: list[dict[str, object]]) -> str:
    ordered = sorted(
        sites,
        key=lambda site: (
            str(site["chain_id"]),
            int(site["resseq"]),
            str(site["icode"]),
            str(site["wt"]),
            str(site["mut"]),
        ),
    )
    return ";".join(
        f"{site['chain_id']}:{site['wt']}{site['resseq']}{site['icode']}{site['mut']}"
        for site in ordered
    )


def _temperature(row: dict[str, str]) -> float:
    match = re.search(r"[-+]?\d+(?:\.\d+)?", row.get("Temperature", ""))
    return float(match.group(0)) if match else 298.0


def _ddg_from_affinity(row: dict[str, str]) -> tuple[float | None, str]:
    mut = row.get("Affinity_mut_parsed", "").strip()
    wt = row.get("Affinity_wt_parsed", "").strip()
    try:
        mut_value = float(mut)
        wt_value = float(wt)
        if mut_value <= 0 or wt_value <= 0:
            raise ValueError
    except ValueError:
        return None, "missing_or_censored_affinity"
    return R_KCAL * _temperature(row) * math.log(mut_value / wt_value), ""


def _experimental_rows(spec: TargetSpec, source_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for index, row in enumerate(source_rows, start=1):
        pdb_key = row.get("#Pdb", "").strip().upper()
        if not pdb_key.startswith(spec.pdb):
            continue
        raw_tokens = _split_tokens(row.get("Mutation(s)_cleaned", ""))
        sites = [_parse_token(token) for token in raw_tokens]
        if not raw_tokens or any(site is None for site in sites):
            continue
        parsed_sites = [site for site in sites if site is not None]
        if not all(str(site["chain_id"]) in spec.antibody_chains for site in parsed_sites):
            continue
        ddg, censor_reason = _ddg_from_affinity(row)
        signature = _canonical_signature(parsed_sites)
        output.append(
            {
                "source_row_id": f"skempi2_{spec.pdb.lower()}_{index:04d}",
                "pdb_id": spec.pdb,
                "mutation_tokens": ";".join(
                    f"{site['chain_id']}:{site['wt']}{site['resseq']}{site['icode']}{site['mut']}"
                    for site in parsed_sites
                ),
                "source_mutation_tokens": row.get("Mutation(s)_cleaned", ""),
                "source_mut_index": "",
                "structure_mapped": True,
                "canonical_signature": signature,
                "mutation_count": len(parsed_sites),
                "single_point": len(parsed_sites) == 1,
                "antibody_chains": ",".join(spec.antibody_chains),
                "antigen_chains": ",".join(spec.antigen_chains),
                "ddg_kcal_mol": "" if ddg is None else f"{ddg:.8f}",
                "ddg_censored": ddg is None,
                "censor_reason": censor_reason,
                "temperature_K": f"{_temperature(row):.2f}",
                "affinity_mut_M": row.get("Affinity_mut_parsed", ""),
                "affinity_wt_M": row.get("Affinity_wt_parsed", ""),
                "reference": row.get("Reference", ""),
                "method": row.get("Method", ""),
                "source_dataset": "SKEMPI2",
                "label_type": "experimental_ddg",
            }
        )
    return output


THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def _pdb_chain_sequences(path: Path) -> dict[str, tuple[str, list[tuple[int, str]]]]:
    chains: dict[str, tuple[list[str], list[tuple[int, str]]]] = {}
    seen: set[tuple[str, int, str]] = set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = THREE_TO_ONE.get(line[17:20].strip().upper())
        if residue is None:
            continue
        chain = line[21].strip()
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            continue
        icode = line[26].strip().upper()
        key = (chain, resseq, icode)
        if key in seen:
            continue
        seen.add(key)
        letters, identifiers = chains.setdefault(chain, ([], []))
        letters.append(residue)
        identifiers.append((resseq, icode))
    return {chain: ("".join(letters), ids) for chain, (letters, ids) in chains.items()}


def _global_sequence_map(query: str, target: str) -> dict[int, int]:
    """Return query-index -> target-index for a Needleman-Wunsch alignment."""
    match_score, mismatch_score, gap_score = 2, -1, -5
    n, m = len(query), len(target)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = i * gap_score
        trace[i][0] = "U"
    for j in range(1, m + 1):
        scores[0][j] = j * gap_score
        trace[0][j] = "L"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i - 1][j - 1] + (
                match_score if query[i - 1] == target[j - 1] else mismatch_score
            )
            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score
            best = max(diagonal, up, left)
            scores[i][j] = best
            trace[i][j] = "D" if best == diagonal else ("U" if best == up else "L")
    mapping: dict[int, int] = {}
    i, j = n, m
    while i or j:
        direction = trace[i][j]
        if direction == "D":
            if query[i - 1] == target[j - 1]:
                mapping[i - 1] = j - 1
            i -= 1
            j -= 1
        elif direction == "U":
            i -= 1
        elif direction == "L":
            j -= 1
        else:  # pragma: no cover - defensive guard for malformed traces
            break
    return mapping


def _synthetic_structure_mapping(
    spec: TargetSpec,
    sequence_rows: list[dict[str, str]],
    structure_path: Path,
) -> dict[tuple[str, int], tuple[int, str, str]]:
    matching = [row for row in sequence_rows if row.get("pdb", "").strip().upper() == spec.pdb]
    if not matching:
        raise ValueError(f"No Graphinity sequence row found for synthetic target {spec.pdb}")
    row = matching[0]
    chain_sequences = _pdb_chain_sequences(structure_path)
    graph_sequences: dict[str, str] = {}
    antibody_chains = list(spec.antibody_chains)
    graph_sequences[antibody_chains[0]] = row.get("ab_chain1_seq", "").strip()
    if len(antibody_chains) > 1:
        graph_sequences[antibody_chains[1]] = row.get("ab_chain2_seq", "").strip()
    graph_sequences[spec.antigen_chains[0]] = row.get("ag_chain_seq", "").strip()

    mapping: dict[tuple[str, int], tuple[int, str, str]] = {}
    for chain, query in graph_sequences.items():
        if chain not in chain_sequences or not query:
            raise ValueError(f"Cannot map {spec.pdb} chain {chain} to the selected RCSB structure")
        target, identifiers = chain_sequences[chain]
        index_map = _global_sequence_map(query, target)
        for query_index, target_index in index_map.items():
            resseq, icode = identifiers[target_index]
            mapping[(chain, query_index)] = (resseq, icode, target[target_index])
    return mapping


def _synthetic_rows(
    spec: TargetSpec,
    source_rows: list[dict[str, str]],
    sequence_rows: list[dict[str, str]],
    structure_mapping: dict[tuple[str, int], tuple[int, str, str]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    sequence_by_complex = {
        row.get("complex", "").strip(): row for row in sequence_rows
    }
    for index, row in enumerate(source_rows, start=1):
        if row.get("pdb", "").strip().upper() != spec.pdb:
            continue
        complex_value = row.get("complex", "").strip()
        token = complex_value.split("_")[-1]
        site = _parse_token(token)
        if site is None or str(site["chain_id"]) not in spec.antibody_chains:
            continue
        sequence_row = sequence_by_complex.get(complex_value, {})
        source_index = int(sequence_row.get("mut_index", row.get("mut_index", "-1")))
        mapped = structure_mapping.get((str(site["chain_id"]), source_index))
        if mapped is None:
            raise ValueError(f"Cannot map {spec.pdb} {token} from Graphinity numbering to PDB numbering")
        mapped_resseq, mapped_icode, mapped_wt = mapped
        if mapped_wt != str(site["wt"]):
            raise ValueError(
                f"WT mismatch while mapping {spec.pdb} {token}: sequence has {mapped_wt}"
            )
        source_token = f"{site['chain_id']}:{site['wt']}{site['resseq']}{site['icode']}{site['mut']}"
        site["resseq"] = mapped_resseq
        site["icode"] = mapped_icode
        signature = _canonical_signature([site])
        output.append(
            {
                "source_row_id": f"graphinity_flexddg_{spec.pdb.lower()}_{index:04d}",
                "pdb_id": spec.pdb,
                "mutation_tokens": f"{site['chain_id']}:{site['wt']}{site['resseq']}{site['icode']}{site['mut']}",
                "source_mutation_tokens": source_token,
                "source_mut_index": source_index,
                "structure_mapped": True,
                "canonical_signature": signature,
                "mutation_count": 1,
                "single_point": True,
                "antibody_chains": ",".join(spec.antibody_chains),
                "antigen_chains": ",".join(spec.antigen_chains),
                "ddg_kcal_mol": row.get("labels", ""),
                "ddg_censored": False,
                "censor_reason": "",
                "temperature_K": "",
                "affinity_mut_M": "",
                "affinity_wt_M": "",
                "reference": "",
                "method": "Rosetta Flex ddG (Graphinity control)",
                "source_dataset": "Graphinity-Flex-ddG",
                "label_type": "synthetic_rosetta_flexddg",
            }
        )
    return output


def _site_rows(group_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in group_rows:
        group_id = str(row["source_row_id"])
        for token in str(row["mutation_tokens"]).split(";"):
            if ":" in token:
                chain, compact = token.split(":", 1)
                token = f"{compact[0]}{chain}{compact[1:]}"
            site = _parse_token(token)
            if site is None:
                continue
            rows.append(
                {
                    "mutation_group_id": group_id,
                    "chain_id": site["chain_id"],
                    "resseq": site["resseq"],
                    "icode": site["icode"],
                    "wt": site["wt"],
                    "mut": site["mut"],
                    "entity_side": "antibody",
                }
            )
    return rows


def _write_system(spec: TargetSpec, structure_path: Path) -> None:
    target_dir = OUT / "targets" / spec.pdb
    target_dir.mkdir(parents=True, exist_ok=True)
    text = "\n".join(
        [
            f"system_name: {spec.pdb.lower()}",
            f"input_structure: {structure_path.relative_to(ROOT)}",
            "structure_source: experimental",
            f"complex_class: {spec.complex_class}",
            f"mutated_partner: {spec.mutated_partner}",
            "antibody_chains:",
            *[f"- {chain}" for chain in spec.antibody_chains],
            "antigen_chains:",
            *[f"- {chain}" for chain in spec.antigen_chains],
            "notes:",
            f"- {spec.name}",
            f"- {spec.notes}",
        ]
    )
    (target_dir / "system.yml").write_text(text + "\n", encoding="utf-8")


def _copy_structure(spec: TargetSpec) -> Path:
    destination = OUT / "structures" / f"{spec.pdb}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        return destination
    source = SKEMPI_PDB_ROOT / f"{spec.pdb}.pdb"
    if not source.exists():
        source = Path("/tmp/skempi_new10") / f"{spec.pdb}.pdb"
    if source.exists():
        shutil.copy2(source, destination)
        return destination
    url = f"https://files.rcsb.org/download/{spec.pdb}.pdb"
    _download(url, destination)
    return destination


def _write_source_metadata(
    skempi: Path,
    graphinity: Path | None = None,
    graphinity_sequences: Path | None = None,
) -> None:
    rows = [
        {
            "source": "SKEMPI2 CSV",
            "url": "https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv",
            "local_input": str(skempi),
            "sha256": _sha256(skempi),
            "retrieved_or_verified": str(date.today()),
        },
        {
            "source": "SKEMPI2 cleaned PDB archive",
            "url": "https://life.bsc.es/pid/skempi2/database/download/SKEMPI2_PDBs.tgz",
            "local_input": "/tmp/skempi_new10/SKEMPI2_PDBs.tgz",
            "sha256": _sha256(Path("/tmp/skempi_new10/SKEMPI2_PDBs.tgz"))
            if Path("/tmp/skempi_new10/SKEMPI2_PDBs.tgz").exists()
            else "",
            "retrieved_or_verified": str(date.today()),
        },
    ]
    if graphinity is not None and graphinity_sequences is not None:
        rows.extend(
            [
                {
                    "source": "Graphinity Flex ddG CSV",
                    "url": "https://github.com/amhummer/Graphinity/blob/main/data/ddg_synthetic/Flex_ddG/Synthetic_FlexddG_ddG_20829.csv",
                    "local_input": str(graphinity),
                    "sha256": _sha256(graphinity),
                    "retrieved_or_verified": str(date.today()),
                },
                {
                    "source": "Graphinity Flex ddG sequence mapping CSV",
                    "url": "https://github.com/amhummer/Graphinity/blob/main/data/ddg_synthetic/Flex_ddG/Synthetic_FlexddG_ddG_20829_with_sequences.csv",
                    "local_input": str(graphinity_sequences),
                    "sha256": _sha256(graphinity_sequences),
                    "retrieved_or_verified": str(date.today()),
                },
            ]
        )
    _write_csv(
        OUT / "source_metadata" / "sources.csv",
        rows,
        ["source", "url", "local_input", "sha256", "retrieved_or_verified"],
    )


def build(skempi_path: Path, graphinity_path: Path, graphinity_sequences_path: Path) -> None:
    if not skempi_path.exists():
        raise FileNotFoundError(
            f"SKEMPI2 CSV not found at {skempi_path}; download it from the URL in README.md"
        )
    uses_graphinity = any(spec.label_type != "experimental" for spec in TARGETS)
    if uses_graphinity and not graphinity_path.exists():
        raise FileNotFoundError(
            f"Graphinity Flex ddG CSV not found at {graphinity_path}; it is only needed for control targets"
        )
    if uses_graphinity and not graphinity_sequences_path.exists():
        raise FileNotFoundError(
            f"Graphinity sequence mapping CSV not found at {graphinity_sequences_path}; it is required to map controls to RCSB residue IDs"
        )

    skempi_rows = _read_csv(skempi_path, delimiter=";", encoding="latin-1")
    graphinity_rows = _read_csv(graphinity_path) if uses_graphinity else []
    graphinity_sequence_rows = _read_csv(graphinity_sequences_path) if uses_graphinity else []
    _write_source_metadata(
        skempi_path,
        graphinity_path if uses_graphinity else None,
        graphinity_sequences_path if uses_graphinity else None,
    )

    manifest_rows: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    all_site_rows: list[dict[str, object]] = []

    for spec in TARGETS:
        if spec.pdb in TRIED_TARGETS:
            raise ValueError(f"Target {spec.pdb} is in the conservative tried-target exclusion set")
        target_dir = OUT / "targets" / spec.pdb
        target_dir.mkdir(parents=True, exist_ok=True)
        structure = _copy_structure(spec)
        _write_system(spec, structure)

        if spec.label_type == "experimental":
            rows = _experimental_rows(spec, skempi_rows)
        else:
            structure_mapping = _synthetic_structure_mapping(
                spec, graphinity_sequence_rows, structure
            )
            rows = _synthetic_rows(
                spec, graphinity_rows, graphinity_sequence_rows, structure_mapping
            )
        if len(rows) < 10 and not spec.allow_below_minimum:
            raise ValueError(f"{spec.pdb} yielded only {len(rows)} selected mutation rows")

        single_rows = [row for row in rows if bool(row["single_point"])]
        unique_signatures = {str(row["canonical_signature"]) for row in rows}
        unique_single_signatures = {
            str(row["canonical_signature"]) for row in single_rows
        }
        multi_rows = [row for row in rows if not bool(row["single_point"])]
        censored_rows = [row for row in rows if bool(row["ddg_censored"])]
        formal_ready = (
            spec.label_type == "experimental"
            and len(unique_single_signatures) >= 10
            and not censored_rows
        )
        if formal_ready:
            quality_tier = "experimental_single_point_ready"
        elif spec.label_type == "experimental" and (
            len(single_rows) >= 10 or spec.allow_below_minimum
        ):
            quality_tier = "experimental_single_point_near_threshold"
        elif spec.label_type == "experimental":
            quality_tier = "experimental_multi_point_exploratory"
        else:
            quality_tier = "synthetic_control_only"

        group_fields = [
            "source_row_id",
            "pdb_id",
            "mutation_tokens",
            "canonical_signature",
            "source_mutation_tokens",
            "source_mut_index",
            "structure_mapped",
            "mutation_count",
            "single_point",
            "antibody_chains",
            "antigen_chains",
            "ddg_kcal_mol",
            "ddg_censored",
            "censor_reason",
            "temperature_K",
            "affinity_mut_M",
            "affinity_wt_M",
            "reference",
            "method",
            "source_dataset",
            "label_type",
        ]
        _write_csv(target_dir / "experimental_ddg.csv", rows, group_fields)
        _write_csv(target_dir / "mutations_all.csv", _site_rows(rows), [
            "mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"
        ])
        _write_csv(target_dir / "mutations_single_all.csv", _site_rows(single_rows), [
            "mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"
        ])
        unique_single_rows: list[dict[str, object]] = []
        seen: set[str] = set()
        for row in single_rows:
            signature = str(row["canonical_signature"])
            if signature in seen:
                continue
            seen.add(signature)
            unique_single_rows.append(row)
        _write_csv(target_dir / "mutations_single_unique.csv", _site_rows(unique_single_rows), [
            "mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"
        ])

        for row in rows:
            all_rows.append(row)
        all_site_rows.extend(_site_rows(rows))
        manifest_rows.append(
            {
                "target_id": spec.pdb,
                "pdb_id": spec.pdb,
                "name": spec.name,
                "structure_path": str(structure.relative_to(ROOT)),
                "antibody_chains": ",".join(spec.antibody_chains),
                "antigen_chains": ",".join(spec.antigen_chains),
                "label_type": spec.label_type,
                "source_dataset": spec.source_dataset,
                "complex_class": spec.complex_class,
                "mutated_partner": spec.mutated_partner,
                "n_total_records": len(rows),
                "n_unique_signatures": len(unique_signatures),
                "n_single_records": len(single_rows),
                "n_unique_single_records": len(unique_single_signatures),
                "n_multi_records": len(multi_rows),
                "n_censored_records": len(censored_rows),
                "formal_ready": formal_ready,
                "quality_tier": quality_tier,
                "related_to_existing": spec.related_to_existing,
                "notes": spec.notes,
            }
        )

    manifest_fields = list(manifest_rows[0].keys())
    _write_csv(OUT / "target_manifest.csv", manifest_rows, manifest_fields)
    all_fields = list(all_rows[0].keys())
    _write_csv(OUT / "selected_mutations.csv", all_rows, all_fields)
    _write_csv(
        OUT / "selected_mutation_sites.csv",
        all_site_rows,
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )

    # A reserve list is useful for future source expansion, but it is not
    # counted among the ten because it has fewer than ten antibody-side rows.
    reserve = [
        {
            "pdb_id": "5C6T",
            "name": "1G2 Fab–HCMV glycoprotein B",
            "antibody_rows": 9,
            "unique_single_rows": 8,
            "reason": "best near-threshold experimental candidate; one additional independent antibody measurement is needed",
        },
        {
            "pdb_id": "4NM8",
            "name": "CR8043 Fab–influenza HA",
            "antibody_rows": 0,
            "unique_single_rows": 0,
            "reason": "SKEMPI mutations are on HA chain B (antigen), so excluded by the antibody-side requirement",
        },
    ]
    _write_csv(
        OUT / "reserve_candidates.csv",
        reserve,
        ["pdb_id", "name", "antibody_rows", "unique_single_rows", "reason"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skempi", type=Path, default=SKEMPI_DEFAULT)
    parser.add_argument("--graphinity", type=Path, default=GRAPHINITY_DEFAULT)
    parser.add_argument(
        "--graphinity-sequences", type=Path, default=GRAPHINITY_SEQUENCES_DEFAULT
    )
    args = parser.parse_args()
    build(args.skempi, args.graphinity, args.graphinity_sequences)
    print(f"Wrote new antibody validation panel to {OUT}")


if __name__ == "__main__":
    main()
