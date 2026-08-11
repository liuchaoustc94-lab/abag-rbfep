"""Mutation parsing and validation for antibody-antigen RBFE."""

from __future__ import annotations

from hashlib import sha1
from pathlib import Path
import re

from abag_rbfe.io_utils import read_csv_rows
from abag_rbfe.models import MutationGroup, MutationSite

STANDARD_AA_CODES = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

NOMINAL_CHARGE = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0,
}

TOKEN_RE = re.compile(
    r"^(?:(?P<chain>[A-Za-z0-9]):)?(?P<wt>[A-Z])(?P<resseq>-?\d+)(?P<icode>[A-Za-z]?)(?P<mut>[A-Z])(?:@(?P<entity_side>antibody|antigen))?$"
)


def _require_standard_residue(code: str) -> str:
    code = code.strip().upper()
    if code not in STANDARD_AA_CODES:
        raise ValueError(f"Unsupported residue code: {code}")
    return code


def nominal_charge(code: str) -> int:
    return NOMINAL_CHARGE.get(_require_standard_residue(code), 0)


def is_charge_conserving(wt: str, mut: str) -> bool:
    return nominal_charge(wt) == nominal_charge(mut)


def parse_mutation_token(token: str, default_entity_side: str | None = None) -> MutationSite:
    match = TOKEN_RE.match(token.strip())
    if not match:
        raise ValueError(f"Unsupported mutation token: {token}")
    entity_side = match.group("entity_side") or default_entity_side
    if entity_side not in {"antibody", "antigen"}:
        raise ValueError(f"Mutation token must resolve to antibody or antigen side: {token}")
    return MutationSite(
        chain_id=(match.group("chain") or "").upper(),
        resseq=int(match.group("resseq")),
        icode=(match.group("icode") or "").upper(),
        wt=_require_standard_residue(match.group("wt")),
        mut=_require_standard_residue(match.group("mut")),
        entity_side=entity_side,
    )


def parse_mutation_group_tokens(tokens: str, default_entity_side: str | None = None) -> list[MutationSite]:
    return [
        parse_mutation_token(token, default_entity_side=default_entity_side)
        for token in tokens.split(";")
        if token.strip()
    ]


def parse_mutation_site_dict(row: dict[str, str]) -> MutationSite:
    required = ("chain_id", "resseq", "wt", "mut", "entity_side")
    missing = [field for field in required if row.get(field, "").strip() == ""]
    if missing:
        raise ValueError(f"Missing mutation fields: {', '.join(missing)}")
    return MutationSite(
        chain_id=row["chain_id"].strip().upper(),
        resseq=int(row["resseq"]),
        icode=row.get("icode", "").strip().upper(),
        wt=_require_standard_residue(row["wt"]),
        mut=_require_standard_residue(row["mut"]),
        entity_side=row["entity_side"].strip().lower(),
    )


def _normalize_sites(sites: list[MutationSite]) -> tuple[MutationSite, ...]:
    ordered = sorted(sites, key=lambda item: (item.entity_side, item.chain_id, item.resseq, item.icode, item.wt, item.mut))
    dedupe_key = [(site.entity_side, site.chain_id, site.resseq, site.icode, site.wt, site.mut) for site in ordered]
    if len(set(dedupe_key)) != len(dedupe_key):
        raise ValueError("Duplicate mutation sites are not allowed inside a mutation group")
    return tuple(ordered)


def build_mutation_group(
    mutation_group_id: str,
    sites: list[MutationSite],
    allow_double_same_side: bool,
    allow_charge_change: bool,
) -> MutationGroup:
    normalized = _normalize_sites(sites)
    if not normalized:
        raise ValueError("Mutation groups must contain at least one mutation site")
    if len(normalized) > 2:
        raise ValueError("Only single-point and double-point mutation groups are supported")
    sides = {site.entity_side for site in normalized}
    if not sides.issubset({"antibody", "antigen"}):
        raise ValueError("Mutation sites must be assigned to antibody or antigen")
    if len(normalized) == 2 and not allow_double_same_side:
        raise ValueError("Double-point mutation groups are not enabled")
    if len(normalized) == 2 and len(sides) != 1:
        raise ValueError("V2 only supports same-side double-point mutation groups")
    charge_conserving = all(is_charge_conserving(site.wt, site.mut) for site in normalized)
    if not charge_conserving and not allow_charge_change:
        raise ValueError("Charge-changing mutations are not enabled in this protocol")
    min_version = "v1" if len(normalized) == 1 else "v2"
    return MutationGroup(
        mutation_group_id=mutation_group_id,
        sites=normalized,
        mutation_count=len(normalized),
        entity_side=normalized[0].entity_side if len(sides) == 1 else "mixed",
        charge_conserving=charge_conserving,
        min_version=min_version,
    )


def load_mutation_groups_from_csv(path: Path) -> list[dict[str, object]]:
    rows = read_csv_rows(path)
    if not rows:
        return []

    if "mutation_tokens" in {field.lower() for field in rows[0].keys()}:
        output = []
        for index, row in enumerate(rows, start=1):
            group_id = row.get("mutation_group_id") or row.get("job_id") or f"group_{index}"
            default_side = row.get("entity_side", "").strip().lower() or None
            output.append(
                {
                    "mutation_group_id": group_id,
                    "sites": parse_mutation_group_tokens(row["mutation_tokens"], default_entity_side=default_side),
                }
            )
        return output

    grouped: dict[str, list[MutationSite]] = {}
    for index, row in enumerate(rows, start=1):
        group_id = row.get("mutation_group_id") or row.get("job_id") or f"group_{index}"
        grouped.setdefault(group_id, []).append(parse_mutation_site_dict(row))
    return [{"mutation_group_id": group_id, "sites": sites} for group_id, sites in grouped.items()]


def build_job_id(system_name: str, mutation_group: MutationGroup) -> str:
    prefix = system_name.lower().replace("_", "-")
    readable = f"{prefix}-{mutation_group.entity_side}-{mutation_group.short_label()}"
    readable = re.sub(r"[^a-z0-9-]+", "-", readable).strip("-")
    if len(readable) > 72:
        digest = sha1(mutation_group.signature().encode("utf-8")).hexdigest()[:8]
        readable = f"{readable[:63].rstrip('-')}-{digest}"
    return readable


def mutation_script_lines(sites: tuple[MutationSite, ...]) -> list[str]:
    lines = []
    for site in sites:
        if site.icode:
            raise ValueError("pmx mutation script generation does not support insertion codes without prior residue mapping")
        lines.append(f"{site.chain_id} {site.resseq} {site.mut}")
    return lines
