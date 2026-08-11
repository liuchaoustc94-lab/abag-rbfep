"""Minimal structure helpers for independent project-local preprocessing."""

from __future__ import annotations

from collections import defaultdict
import json
from math import dist
from pathlib import Path
import re

STRUCTURE_HEADER_PREFIXES = (
    "MODEL ",
    "ENDMDL",
    "END",
    "CRYST1",
    "HEADER",
    "TITLE ",
    "REMARK",
    "COMPND",
    "SOURCE",
    "KEYWDS",
    "EXPDTA",
    "AUTHOR",
)

RESIDUE_ALIASES = {
    "ASH": "ASP",
    "CYM": "CYS",
    "CYX": "CYS",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "LYN": "LYS",
}

EXPECTED_HEAVY_ATOMS = {
    "ALA": {"N", "CA", "C", "O", "CB"},
    "ARG": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "ASP": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "CYS": {"N", "CA", "C", "O", "CB", "SG"},
    "GLN": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "GLY": {"N", "CA", "C", "O"},
    "HIS": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "LYS": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "MET": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "PHE": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "SER": {"N", "CA", "C", "O", "CB", "OG"},
    "THR": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "TRP": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}

BACKBONE_HEAVY_ATOMS = frozenset({"N", "CA", "C", "O"})
BACKBONE_RETAINED_ATOMS = BACKBONE_HEAVY_ATOMS.union({"OXT", "OT1", "OT2"})
CLASH_CLASS_BACKBONE_BACKBONE = "backbone_backbone"
CLASH_CLASS_BACKBONE_SIDECHAIN = "backbone_sidechain"
CLASH_CLASS_SIDECHAIN_SIDECHAIN = "sidechain_sidechain"
TERMINAL_OXYGEN_ATOMS = frozenset({"OXT", "OT1", "OT2"})
INTRA_RESIDUE_HEAVY_ATOM_MIN_DISTANCE_ANGSTROM = 1.05
INTER_RESIDUE_HEAVY_ATOM_MIN_DISTANCE_ANGSTROM = 1.15
PDB_RESIDUE_LOCATION_RE = re.compile(r"^\s*(?P<resseq>-?\d+)\s*(?P<icode>[A-Za-z]?)\s*$")


def classify_clash_atom(atom_name: object) -> str:
    """Classify a PDB atom name as backbone or sidechain for clash stratification."""
    return "backbone" if str(atom_name).strip().upper() in BACKBONE_RETAINED_ATOMS else "sidechain"


def classify_clash_pair(atom_a: object, atom_b: object) -> str:
    """Stratify a clash pair: backbone_backbone (hard fail), backbone_sidechain,
    or sidechain_sidechain (both repairable by sidechain rebuild)."""
    a_backbone = classify_clash_atom(atom_a) == "backbone"
    b_backbone = classify_clash_atom(atom_b) == "backbone"
    if a_backbone and b_backbone:
        return CLASH_CLASS_BACKBONE_BACKBONE
    if a_backbone or b_backbone:
        return CLASH_CLASS_BACKBONE_SIDECHAIN
    return CLASH_CLASS_SIDECHAIN_SIDECHAIN


def _is_protein_atom_record(line: str) -> bool:
    return line.startswith("ATOM  ") and len(line) >= 22


def _slice(line: str, start: int, end: int) -> str:
    if start >= len(line):
        return ""
    return line[start:end]


def _normalize_pdb_atom_line(
    line: str,
    *,
    clear_altloc: bool = False,
    serial_override: int | None = None,
) -> str:
    record = _slice(line, 0, 6) or "ATOM  "
    serial = serial_override if serial_override is not None else int((_slice(line, 6, 11) or "0").strip())
    atom_field = (_slice(line, 12, 16) or "    ")[:4].ljust(4)
    altloc = " " if clear_altloc else (_slice(line, 16, 17) or " ")[:1]
    resname = (_slice(line, 17, 20) or "   ")[:3].rjust(3)
    chain = (_slice(line, 21, 22) or " ")[:1]
    resseq, parsed_icode = _parse_residue_location(line)
    icode = parsed_icode if parsed_icode else " "
    x = float((_slice(line, 30, 38) or "0").strip())
    y = float((_slice(line, 38, 46) or "0").strip())
    z = float((_slice(line, 46, 54) or "0").strip())
    occupancy_raw = _slice(line, 54, 60).strip()
    bfactor_raw = _slice(line, 60, 66).strip()
    occupancy = float(occupancy_raw) if occupancy_raw else 1.0
    bfactor = float(bfactor_raw) if bfactor_raw else 0.0
    element = _slice(line, 76, 78).strip()
    charge = _slice(line, 78, 80).strip()
    return (
        f"{record:<6}{serial:>5d} {atom_field}{altloc}{resname} {chain}{resseq:>4d}{icode}"
        f"   {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}"
        f"          {element:>2}{charge:>2}\n"
    )


def _parse_residue_location(line: str) -> tuple[int, str]:
    raw = _slice(line, 22, 27)
    match = PDB_RESIDUE_LOCATION_RE.match(raw)
    if match is not None:
        return int(match.group("resseq")), (match.group("icode") or "").upper()

    fallback_resseq = (_slice(line, 22, 26) or "").strip()
    fallback_icode = (_slice(line, 26, 27) or "").strip().upper()
    if fallback_resseq:
        return int(fallback_resseq), fallback_icode
    return 0, fallback_icode


def _parse_resseq(line: str) -> int | None:
    try:
        resseq, _icode = _parse_residue_location(line)
        return resseq
    except ValueError:
        return None


def _parse_icode(line: str) -> str:
    try:
        _resseq, icode = _parse_residue_location(line)
    except ValueError:
        return (_slice(line, 26, 27) or "").strip().upper()
    return icode


def _residue_selector_key(item: dict[str, object]) -> tuple[str, int, str]:
    return (
        str(item.get("chain_id") or "").strip(),
        int(item.get("resseq") or 0),
        str(item.get("icode") or "").strip().upper(),
    )


def _preferred_altloc_choice_key(line: str) -> tuple[int, float, str, int]:
    altloc = (_slice(line, 16, 17) or " ").strip().upper()
    occupancy_raw = _slice(line, 54, 60).strip()
    occupancy = float(occupancy_raw) if occupancy_raw else 0.0
    serial_raw = _slice(line, 6, 11).strip()
    serial = int(serial_raw) if serial_raw else 0
    # Prefer blank altlocs, then higher occupancy, then deterministic serial order.
    return (0 if not altloc else 1, -occupancy, altloc or " ", serial)


def _preferred_altloc_atom_indices(lines: list[str], *, keep_chains: set[str] | None = None) -> set[int]:
    atom_groups: dict[tuple[str, int, str, str, str], list[tuple[int, str]]] = defaultdict(list)
    for index, line in enumerate(lines):
        if not _is_protein_atom_record(line):
            continue
        chain_id = _slice(line, 21, 22).strip()
        if keep_chains is not None and chain_id not in keep_chains:
            continue
        try:
            resseq, icode = _parse_residue_location(line)
        except ValueError:
            continue
        resname = _slice(line, 17, 20).strip().upper()
        atom_name = _slice(line, 12, 16).strip().upper()
        atom_groups[(chain_id, resseq, icode, resname, atom_name)].append((index, line))

    selected_indices: set[int] = set()
    for candidates in atom_groups.values():
        chosen_index, _ = min(candidates, key=lambda item: _preferred_altloc_choice_key(item[1]))
        selected_indices.add(chosen_index)
    return selected_indices


def pdb_chain_ids(path: Path) -> set[str]:
    chains: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if _is_protein_atom_record(line):
                chain = line[21].strip()
                if chain:
                    chains.add(chain)
    return chains


def extract_pdb_chains(input_path: Path, output_path: Path, keep_chains: list[str]) -> None:
    keep = {chain.strip() for chain in keep_chains if chain.strip()}
    if not keep:
        raise ValueError("At least one chain must be retained")

    source_lines = input_path.read_text(encoding="utf-8").splitlines(keepends=True)
    selected_indices = _preferred_altloc_atom_indices(source_lines, keep_chains=keep)
    output_lines: list[str] = []
    segment_has_kept_atoms = False
    for index, line in enumerate(source_lines):
        if _is_protein_atom_record(line):
            chain = line[21].strip()
            if chain in keep and index in selected_indices:
                output_lines.append(_normalize_pdb_atom_line(line, clear_altloc=True))
                segment_has_kept_atoms = True
            continue
        if line.startswith("TER"):
            if segment_has_kept_atoms:
                output_lines.append(line)
            segment_has_kept_atoms = False
            continue
        if line.startswith(STRUCTURE_HEADER_PREFIXES):
            output_lines.append(line)

    if not output_lines or not any(line.startswith("ATOM  ") for line in output_lines):
        raise ValueError(f"No retained atoms remain after filtering chains {sorted(keep)} from {input_path}")

    if not any(line.startswith("END") for line in output_lines):
        output_lines.append("END\n")
    output_path.write_text("".join(output_lines), encoding="utf-8")


def _renumber_protein_atom_serials(lines: list[str]) -> list[str]:
    output_lines: list[str] = []
    next_serial = 1
    for line in lines:
        if _is_protein_atom_record(line):
            output_lines.append(_normalize_pdb_atom_line(line, serial_override=next_serial))
            next_serial += 1
            continue
        output_lines.append(line)
    return output_lines


def _collect_preferred_residue_blocks(
    lines: list[str],
) -> dict[tuple[str, int, str], dict[str, object]]:
    selected_indices = _preferred_altloc_atom_indices(lines)
    residue_blocks: dict[tuple[str, int, str], dict[str, object]] = {}
    for index, line in enumerate(lines):
        if not _is_protein_atom_record(line) or index not in selected_indices:
            continue
        try:
            resseq, icode = _parse_residue_location(line)
        except ValueError:
            continue
        chain_id = _slice(line, 21, 22).strip()
        resname = _slice(line, 17, 20).strip().upper()
        normalized_resname = RESIDUE_ALIASES.get(resname, resname)
        key = (chain_id, resseq, icode)
        block = residue_blocks.setdefault(
            key,
            {
                "chain_id": chain_id,
                "resseq": resseq,
                "icode": icode,
                "resname": resname,
                "normalized_resname": normalized_resname,
                "lines": [],
            },
        )
        block["lines"].append(_normalize_pdb_atom_line(line, clear_altloc=True))
    return residue_blocks


def restore_incomplete_standard_residues_from_template(
    template_path: Path,
    target_path: Path,
    output_path: Path | None = None,
    *,
    exclude_residues: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    output_path = output_path or target_path
    target_lines = target_path.read_text(encoding="utf-8").splitlines(keepends=True)
    target_issues = classify_incomplete_standard_residues(target_path)
    excluded_keys = {_residue_selector_key(item) for item in (exclude_residues or [])}

    attempted_issues = [issue for issue in target_issues if _residue_selector_key(issue) not in excluded_keys]
    summary: dict[str, object] = {
        "attempted": True,
        "attempted_residue_count": len(attempted_issues),
        "skipped_excluded_residue_count": len(target_issues) - len(attempted_issues),
        "skipped_excluded_residues": [
            issue for issue in target_issues if _residue_selector_key(issue) in excluded_keys
        ],
        "restored_residue_count": 0,
        "restored_residues": [],
        "missing_template_residue_count": 0,
        "missing_template_residues": [],
        "unresolved_residue_count": 0,
        "unresolved_residues": [],
    }

    if not attempted_issues:
        if output_path != target_path:
            output_path.write_text("".join(target_lines), encoding="utf-8")
        return summary

    template_blocks = _collect_preferred_residue_blocks(template_path.read_text(encoding="utf-8").splitlines(keepends=True))
    replacements: dict[tuple[str, int, str], list[str]] = {}

    for issue in attempted_issues:
        key = _residue_selector_key(issue)
        template_block = template_blocks.get(key)
        if template_block is None:
            summary["missing_template_residues"].append(issue)
            summary["unresolved_residues"].append(issue)
            continue
        target_normalized = str(issue.get("normalized_resname") or issue.get("resname") or "").strip().upper()
        template_normalized = str(template_block.get("normalized_resname") or "").strip().upper()
        if target_normalized != template_normalized:
            unresolved_issue = dict(issue)
            unresolved_issue["template_resname"] = template_block.get("resname")
            unresolved_issue["template_normalized_resname"] = template_block.get("normalized_resname")
            summary["unresolved_residues"].append(unresolved_issue)
            continue
        replacements[key] = list(template_block["lines"])
        summary["restored_residues"].append(
            {
                "chain_id": key[0],
                "resseq": key[1],
                "icode": key[2],
                "resname": issue.get("resname"),
                "normalized_resname": issue.get("normalized_resname"),
                "missing_atoms": list(issue.get("missing_atoms", [])),
            }
        )

    if replacements:
        rebuilt_lines: list[str] = []
        index = 0
        while index < len(target_lines):
            line = target_lines[index]
            if not _is_protein_atom_record(line):
                if not line.startswith("END"):
                    rebuilt_lines.append(line)
                index += 1
                continue
            try:
                residue_key = (line[21].strip(), *_parse_residue_location(line))
            except ValueError:
                rebuilt_lines.append(line)
                index += 1
                continue
            residue_lines: list[str] = []
            while index < len(target_lines):
                current_line = target_lines[index]
                if not _is_protein_atom_record(current_line):
                    break
                try:
                    current_key = (current_line[21].strip(), *_parse_residue_location(current_line))
                except ValueError:
                    residue_lines.append(current_line)
                    index += 1
                    continue
                if current_key != residue_key:
                    break
                residue_lines.append(current_line)
                index += 1
            rebuilt_lines.extend(replacements.get(residue_key, residue_lines))
        final_lines = _renumber_protein_atom_serials(rebuilt_lines)
    else:
        final_lines = target_lines

    if not any(line.startswith("END") for line in final_lines):
        final_lines.append("END\n")
    output_path.write_text("".join(final_lines), encoding="utf-8")

    summary["restored_residue_count"] = len(summary["restored_residues"])
    summary["missing_template_residue_count"] = len(summary["missing_template_residues"])
    summary["unresolved_residue_count"] = len(summary["unresolved_residues"])
    return summary


def strip_sidechain_atoms_for_residues(
    input_path: Path,
    output_path: Path,
    residues: list[dict[str, object]],
) -> None:
    targets = {
        (
            str(item["chain_id"]).strip(),
            int(item["resseq"]),
            str(item.get("icode", "")).strip().upper(),
        )
        for item in residues
    }
    if not targets:
        if input_path != output_path:
            output_path.write_text(input_path.read_text(encoding="utf-8"), encoding="utf-8")
        return

    output_lines: list[str] = []
    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if _is_protein_atom_record(line):
                chain_id = line[21].strip()
                try:
                    resseq, icode = _parse_residue_location(line)
                except ValueError:
                    output_lines.append(line)
                    continue
                atom_name = line[12:16].strip().upper()
                if (chain_id, resseq, icode) in targets and atom_name not in BACKBONE_RETAINED_ATOMS:
                    continue
                output_lines.append(_normalize_pdb_atom_line(line))
                continue
            output_lines.append(line)

    if not any(line.startswith("END") for line in output_lines):
        output_lines.append("END\n")
    output_path.write_text("".join(output_lines), encoding="utf-8")


def empty_repair_summary() -> dict[str, object]:
    """Field-complete placeholder for a repair operation that was not triggered.

    Keeps the schema identical to a triggered repair summary so downstream
    consumers never have to branch on missing keys.
    """
    return {
        "attempted": False,
        "available": False,
        "succeeded": False,
        "trigger_residue_count": 0,
        "trigger_residues": [],
        "blocking_residue_count": 0,
        "blocking_residues": [],
        "remaining_incomplete_standard_residue_count": 0,
        "remaining_incomplete_standard_residues": [],
    }


def repair_missing_atoms_with_pdbfixer(input_path: Path, output_path: Path) -> dict[str, object]:
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:  # pragma: no cover - optional dependency guard
        return {
            "attempted": False,
            "available": False,
            "succeeded": False,
            "error": f"{type(exc).__name__}: {exc}",
        }

    try:
        fixer = PDBFixer(filename=str(input_path))
        fixer.findMissingResidues()
        missing_residue_count = sum(len(items) for items in fixer.missingResidues.values())
        fixer.missingResidues = {}
        fixer.findNonstandardResidues()
        nonstandard_residue_count = len(fixer.nonstandardResidues)
        if fixer.nonstandardResidues:
            fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        missing_atom_count = sum(len(items) for items in fixer.missingAtoms.values())
        missing_terminal_atom_count = sum(len(items) for items in fixer.missingTerminals.values())
        if missing_atom_count or missing_terminal_atom_count or nonstandard_residue_count:
            fixer.addMissingAtoms()
        with output_path.open("w", encoding="utf-8") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        return {
            "attempted": True,
            "available": True,
            "succeeded": True,
            "missing_residue_count_detected": missing_residue_count,
            "missing_atom_count_detected": missing_atom_count,
            "missing_terminal_atom_count_detected": missing_terminal_atom_count,
            "nonstandard_residue_count_detected": nonstandard_residue_count,
        }
    except Exception as exc:  # pragma: no cover - runtime repair fallback
        return {
            "attempted": True,
            "available": True,
            "succeeded": False,
            "error": f"{type(exc).__name__}: {exc}",
        }


def _repair_incomplete_standard_residues_with_pdbfixer(
    input_path: Path,
    output_path: Path,
    repairable_issues: list[dict[str, object]],
    *,
    issues_before: list[dict[str, object]] | None = None,
    blocking_issues: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    issues_before = list(issues_before) if issues_before is not None else classify_incomplete_standard_residues(input_path)
    blocking_issues = (
        list(blocking_issues)
        if blocking_issues is not None
        else [issue for issue in issues_before if issue not in repairable_issues]
    )
    summary: dict[str, object] = {
        "attempted": False,
        "available": False,
        "succeeded": False,
        "trigger_residue_count": len(repairable_issues),
        "trigger_residues": repairable_issues,
        "blocking_residue_count": len(blocking_issues),
        "blocking_residues": blocking_issues,
        "remaining_incomplete_standard_residue_count": len(issues_before),
        "remaining_incomplete_standard_residues": issues_before,
    }

    if not repairable_issues:
        if input_path != output_path:
            output_path.write_text(input_path.read_text(encoding="utf-8"), encoding="utf-8")
        return summary

    strip_sidechain_atoms_for_residues(input_path, output_path, repairable_issues)
    summary.update(repair_missing_atoms_with_pdbfixer(output_path, output_path))
    if summary.get("succeeded"):
        strip_terminal_oxygen_atoms(output_path, output_path)
    issues_after = classify_incomplete_standard_residues(output_path)
    summary["remaining_incomplete_standard_residue_count"] = len(issues_after)
    summary["remaining_incomplete_standard_residues"] = issues_after
    return summary


def repair_incomplete_standard_residues_with_pdbfixer(
    input_path: Path,
    output_path: Path,
) -> dict[str, object]:
    issues_before = classify_incomplete_standard_residues(input_path)
    blocking_issues = [issue for issue in issues_before if issue["blocking_prepare"]]
    return _repair_incomplete_standard_residues_with_pdbfixer(
        input_path,
        output_path,
        issues_before,
        issues_before=issues_before,
        blocking_issues=blocking_issues,
    )


def repair_sidechain_only_incomplete_residues_with_pdbfixer(
    input_path: Path,
    output_path: Path,
) -> dict[str, object]:
    issues_before = classify_incomplete_standard_residues(input_path)
    blocking_issues = [issue for issue in issues_before if issue["blocking_prepare"]]
    sidechain_only_issues = [issue for issue in issues_before if not issue["blocking_prepare"]]
    return _repair_incomplete_standard_residues_with_pdbfixer(
        input_path,
        output_path,
        sidechain_only_issues,
        issues_before=issues_before,
        blocking_issues=blocking_issues,
    )


def strip_hydrogen_atoms(input_path: Path, output_path: Path) -> dict[str, object]:
    """Remove all hydrogen atoms from protein ATOM records.

    Hydrogens are stripped at prepare time so pdb2gmx can run WITHOUT -ignh:
    normal residues get hydrogens rebuilt from the force-field .hdb rules,
    while pmx hybrid residues keep the explicit hydrogens (and alchemical
    dummy atoms) written by pmx mutate. Using -ignh instead would silently
    discard the hybrid-residue hydrogens because hybrid residues have no
    .hdb entries, producing chemically broken heavy-atom-only hybrids.
    Returns a small provenance summary.
    """
    output_lines: list[str] = []
    removed: list[dict[str, object]] = []
    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if _is_protein_atom_record(line):
                atom_name = line[12:16].strip()
                element = line[76:78].strip().upper() if len(line) >= 78 else ""
                is_hydrogen = element == "H" or (
                    not element and atom_name.upper().startswith("H")
                )
                if is_hydrogen:
                    removed.append(
                        {
                            "chain_id": line[21].strip(),
                            "resseq": int(line[22:26]),
                            "resname": line[17:20].strip(),
                            "atom_name": atom_name,
                        }
                    )
                    continue
                output_lines.append(_normalize_pdb_atom_line(line))
                continue
            output_lines.append(line)

    if not any(line.startswith("END") for line in output_lines):
        output_lines.append("END\n")
    output_path.write_text("".join(output_lines), encoding="utf-8")
    return {
        "stripped_hydrogen_count": len(removed),
        "stripped_hydrogen_atoms": removed,
    }


def strip_terminal_oxygen_atoms(input_path: Path, output_path: Path) -> None:
    output_lines: list[str] = []
    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if _is_protein_atom_record(line):
                atom_name = line[12:16].strip().upper()
                if atom_name in TERMINAL_OXYGEN_ATOMS:
                    continue
                output_lines.append(_normalize_pdb_atom_line(line))
                continue
            output_lines.append(line)

    if not any(line.startswith("END") for line in output_lines):
        output_lines.append("END\n")
    output_path.write_text("".join(output_lines), encoding="utf-8")


def classify_incomplete_standard_residues(path: Path) -> list[dict[str, object]]:
    residue_atoms: dict[tuple[str, int, str, str], set[str]] = defaultdict(set)
    residue_names: dict[tuple[str, int, str, str], tuple[str, str]] = {}

    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    selected_indices = _preferred_altloc_atom_indices(lines)
    for index, line in enumerate(lines):
        if not _is_protein_atom_record(line) or index not in selected_indices:
            continue
        atom_name = line[12:16].strip()
        if atom_name.startswith("H"):
            continue
        chain_id = line[21].strip()
        resname = line[17:20].strip().upper()
        normalized = RESIDUE_ALIASES.get(resname, resname)
        expected = EXPECTED_HEAVY_ATOMS.get(normalized)
        if expected is None:
            continue
        try:
            resseq, icode = _parse_residue_location(line)
        except ValueError:
            continue
        key = (chain_id, resseq, icode, normalized)
        residue_names[key] = (resname, normalized)
        residue_atoms[key].add(atom_name)

    incomplete: list[dict[str, object]] = []
    for key in sorted(residue_atoms, key=lambda item: (item[0], item[1], item[2], item[3])):
        expected = EXPECTED_HEAVY_ATOMS[key[3]]
        missing_atoms = sorted(expected - residue_atoms[key])
        if not missing_atoms:
            continue
        missing_backbone_atoms = sorted(BACKBONE_HEAVY_ATOMS.intersection(missing_atoms))
        missing_sidechain_atoms = sorted(set(missing_atoms) - BACKBONE_HEAVY_ATOMS)
        resname, normalized = residue_names[key]
        incomplete.append(
            {
                "chain_id": key[0],
                "resseq": key[1],
                "icode": key[2],
                "resname": resname,
                "normalized_resname": normalized,
                "missing_atoms": missing_atoms,
                "missing_backbone_atoms": missing_backbone_atoms,
                "missing_sidechain_atoms": missing_sidechain_atoms,
                "blocking_prepare": bool(missing_backbone_atoms),
            }
        )
    return incomplete


def find_incomplete_standard_residues(path: Path) -> list[dict[str, object]]:
    return classify_incomplete_standard_residues(path)


def _collect_standard_residue_heavy_atoms(
    path: Path,
) -> tuple[
    dict[tuple[str, int, str, str], list[tuple[str, tuple[float, float, float]]]],
    dict[tuple[str, int, str, str], tuple[str, str]],
]:
    residue_atoms: dict[tuple[str, int, str, str], list[tuple[str, tuple[float, float, float]]]] = defaultdict(list)
    residue_names: dict[tuple[str, int, str, str], tuple[str, str]] = {}

    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    selected_indices = _preferred_altloc_atom_indices(lines)
    for index, line in enumerate(lines):
        if not _is_protein_atom_record(line) or index not in selected_indices:
            continue
        atom_name = line[12:16].strip().upper()
        if atom_name.startswith("H"):
            continue
        resname = line[17:20].strip().upper()
        normalized = RESIDUE_ALIASES.get(resname, resname)
        if normalized not in EXPECTED_HEAVY_ATOMS:
            continue
        chain_id = line[21].strip()
        try:
            resseq, icode = _parse_residue_location(line)
        except ValueError:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        key = (chain_id, resseq, icode, normalized)
        residue_names[key] = (resname, normalized)
        residue_atoms[key].append((atom_name, (x, y, z)))
    return residue_atoms, residue_names


def find_intra_residue_heavy_atom_clashes(path: Path) -> list[dict[str, object]]:
    residue_atoms, residue_names = _collect_standard_residue_heavy_atoms(path)

    clashes: list[dict[str, object]] = []
    threshold = INTRA_RESIDUE_HEAVY_ATOM_MIN_DISTANCE_ANGSTROM
    for key in sorted(residue_atoms, key=lambda item: (item[0], item[1], item[2], item[3])):
        atoms = residue_atoms[key]
        residue_clashes: list[dict[str, object]] = []
        min_distance: float | None = None
        for index, (atom_a, coords_a) in enumerate(atoms):
            for atom_b, coords_b in atoms[index + 1 :]:
                distance = dist(coords_a, coords_b)
                if distance >= threshold:
                    continue
                residue_clashes.append(
                    {
                        "atom_a": atom_a,
                        "atom_b": atom_b,
                        "distance_angstrom": round(distance, 4),
                        "atom_a_class": classify_clash_atom(atom_a),
                        "atom_b_class": classify_clash_atom(atom_b),
                        "clash_class": classify_clash_pair(atom_a, atom_b),
                    }
                )
                min_distance = distance if min_distance is None else min(min_distance, distance)
        if not residue_clashes:
            continue
        resname, normalized = residue_names[key]
        clashes.append(
            {
                "chain_id": key[0],
                "resseq": key[1],
                "icode": key[2],
                "resname": resname,
                "normalized_resname": normalized,
                "min_distance_angstrom": round(min_distance or 0.0, 4),
                "clashes": residue_clashes,
                "clash_classes": sorted({str(clash["clash_class"]) for clash in residue_clashes}),
                "blocking_prepare": True,
            }
        )
    return clashes


def find_inter_residue_heavy_atom_clashes(path: Path) -> list[dict[str, object]]:
    residue_atoms, residue_names = _collect_standard_residue_heavy_atoms(path)

    clashes: list[dict[str, object]] = []
    threshold = INTER_RESIDUE_HEAVY_ATOM_MIN_DISTANCE_ANGSTROM
    residue_keys = sorted(residue_atoms, key=lambda item: (item[0], item[1], item[2], item[3]))
    for index, key_a in enumerate(residue_keys):
        atoms_a = residue_atoms[key_a]
        for key_b in residue_keys[index + 1 :]:
            atoms_b = residue_atoms[key_b]
            pair_clashes: list[dict[str, object]] = []
            min_distance: float | None = None
            for atom_a, coords_a in atoms_a:
                for atom_b, coords_b in atoms_b:
                    distance = dist(coords_a, coords_b)
                    if distance >= threshold:
                        continue
                    pair_clashes.append(
                        {
                            "atom_a": atom_a,
                            "atom_b": atom_b,
                            "distance_angstrom": round(distance, 4),
                            "atom_a_class": classify_clash_atom(atom_a),
                            "atom_b_class": classify_clash_atom(atom_b),
                            "clash_class": classify_clash_pair(atom_a, atom_b),
                        }
                    )
                    min_distance = distance if min_distance is None else min(min_distance, distance)
            if not pair_clashes:
                continue
            resname_a, normalized_a = residue_names[key_a]
            resname_b, normalized_b = residue_names[key_b]
            clashes.append(
                {
                    "chain_id": key_a[0],
                    "resseq": key_a[1],
                    "icode": key_a[2],
                    "resname": resname_a,
                    "normalized_resname": normalized_a,
                    "partner_chain_id": key_b[0],
                    "partner_resseq": key_b[1],
                    "partner_icode": key_b[2],
                    "partner_resname": resname_b,
                    "partner_normalized_resname": normalized_b,
                    "min_distance_angstrom": round(min_distance or 0.0, 4),
                    "clashes": pair_clashes,
                    "clash_classes": sorted({str(clash["clash_class"]) for clash in pair_clashes}),
                    "blocking_prepare": True,
                }
            )
    return clashes


def _inter_residue_clash_entry_key(issue: dict[str, object], clash: dict[str, object]) -> tuple[object, ...]:
    return (
        str(issue.get("chain_id") or "").strip(),
        int(issue.get("resseq") or 0),
        str(issue.get("icode") or "").strip().upper(),
        str(issue.get("partner_chain_id") or "").strip(),
        int(issue.get("partner_resseq") or 0),
        str(issue.get("partner_icode") or "").strip().upper(),
        str(clash.get("atom_a") or "").strip().upper(),
        str(clash.get("atom_b") or "").strip().upper(),
    )


def filter_preexisting_inter_residue_heavy_atom_clashes(
    candidate_clashes: list[dict[str, object]],
    reference_clashes: list[dict[str, object]],
) -> list[dict[str, object]]:
    reference_keys = {
        _inter_residue_clash_entry_key(issue, clash)
        for issue in reference_clashes
        for clash in issue.get("clashes", [])
        if isinstance(clash, dict)
    }

    filtered: list[dict[str, object]] = []
    for issue in candidate_clashes:
        retained_clashes = [
            clash
            for clash in issue.get("clashes", [])
            if isinstance(clash, dict) and _inter_residue_clash_entry_key(issue, clash) not in reference_keys
        ]
        if not retained_clashes:
            continue
        filtered_issue = dict(issue)
        filtered_issue["clashes"] = retained_clashes
        filtered_issue["min_distance_angstrom"] = min(
            float(clash.get("distance_angstrom", 0.0))
            for clash in retained_clashes
        )
        filtered.append(filtered_issue)
    return filtered


def write_inter_residue_heavy_atom_clash_report(
    input_path: Path,
    output_path: Path,
    *,
    reference_path: Path | None = None,
) -> dict[str, object]:
    clashes = find_inter_residue_heavy_atom_clashes(input_path)
    if reference_path is not None:
        reference_clashes = find_inter_residue_heavy_atom_clashes(reference_path)
        clashes = filter_preexisting_inter_residue_heavy_atom_clashes(clashes, reference_clashes)
    payload = {
        "input_structure": str(input_path),
        "reference_structure": str(reference_path) if reference_path is not None else None,
        "inter_residue_heavy_atom_clashes": clashes,
        "blocking_inter_residue_heavy_atom_clashes": bool(clashes),
    }
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return payload


def partition_sidechain_repairable_clashes(
    clashes: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    repairable: list[dict[str, object]] = []
    blocking: list[dict[str, object]] = []
    for issue in clashes:
        normalized_resname = str(issue.get("normalized_resname") or issue.get("resname") or "").strip().upper()
        clash_entries = issue.get("clashes", [])
        if not clash_entries or normalized_resname == "GLY":
            blocking.append(issue)
            continue

        repairable_issue = True
        for clash in clash_entries:
            # Only sidechain-involving clashes are repairable by stripping and
            # rebuilding the sidechain; backbone-backbone collisions remain hard failures.
            if classify_clash_pair(clash.get("atom_a", ""), clash.get("atom_b", "")) == CLASH_CLASS_BACKBONE_BACKBONE:
                repairable_issue = False
                break
        if repairable_issue:
            repairable.append(issue)
        else:
            blocking.append(issue)
    return repairable, blocking


def partition_inter_residue_sidechain_repairable_clashes(
    clashes: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    repairable: list[dict[str, object]] = []
    blocking: list[dict[str, object]] = []
    for issue in clashes:
        clash_entries = issue.get("clashes", [])
        if not clash_entries:
            blocking.append(issue)
            continue

        repairable_issue = True
        for clash in clash_entries:
            # Inter-residue clashes remain repairable as long as at least one sidechain
            # atom can be stripped and rebuilt; backbone-backbone collisions still hard-fail.
            if classify_clash_pair(clash.get("atom_a", ""), clash.get("atom_b", "")) == CLASH_CLASS_BACKBONE_BACKBONE:
                repairable_issue = False
                break
        if repairable_issue:
            repairable.append(issue)
        else:
            blocking.append(issue)
    return repairable, blocking
