"""Helpers for working with local GROMACS installations."""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from abag_rbfe.io_utils import ensure_dir
from abag_rbfe.structure import EXPECTED_HEAVY_ATOMS, RESIDUE_ALIASES

WATER_COORDINATE_FILES = {
    "spc": "spc216.gro",
    "spce": "spc216.gro",
    "tip3p": "spc216.gro",
}

RTP_SECTION_EXCLUSIONS = {
    "angles",
    "atoms",
    "bonds",
    "bondedtypes",
    "cmap",
    "dihedrals",
    "exclusions",
    "impropers",
    "pairs",
}

DEDUP_INCLUDE_BASENAMES = {
    "ions.itp",
    "spc.itp",
    "spce.itp",
    "tip3p.itp",
    "tip4p.itp",
    "tip4pew.itp",
    "tip5p.itp",
}

CHAIN_ITP_BASENAME_RE = re.compile(r"^topol_Protein_chain_(?P<chain>.+)\.itp$", re.IGNORECASE)
PROTEIN_BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}


def force_field_dir_name(force_field: str) -> str:
    if force_field.endswith(".ff"):
        return force_field
    return f"{force_field}.ff"


def resolve_gmx_binary(gmx_bin: str) -> Path | None:
    resolved = shutil.which(gmx_bin)
    if resolved:
        return Path(resolved).resolve()
    candidate = Path(gmx_bin).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    return None


def discover_gmx_top_dir(gmx_bin: str) -> Path | None:
    binary = resolve_gmx_binary(gmx_bin)
    if binary is not None:
        top_dir = binary.parents[1] / "share" / "gromacs" / "top"
        if top_dir.is_dir():
            return top_dir

    gmxdata = os.environ.get("GMXDATA")
    if gmxdata:
        candidate = Path(gmxdata) / "top"
        if candidate.is_dir():
            return candidate.resolve()

    gmxlib = os.environ.get("GMXLIB")
    if gmxlib:
        candidate = Path(gmxlib)
        if candidate.is_dir():
            return candidate.resolve()

    return None


def ensure_local_gmxlib(job_dir: Path, gmx_top_dir: Path, pmx_mutff_root: Path, force_field: str) -> Path:
    """Create a local GMXLIB overlay for this job.

    The overlay mirrors the upstream GROMACS top directory via symlinks and
    injects the pmx mutation force field so ``pdb2gmx`` and subsequent GROMACS
    tools can resolve both standard data files and hybrid residue definitions.
    """

    gmxlib_dir = ensure_dir(job_dir / "artifacts" / "gmxlib")
    ff_dir_name = force_field_dir_name(force_field)
    pmx_ff_dir = pmx_mutff_root / ff_dir_name
    if not pmx_ff_dir.is_dir():
        raise FileNotFoundError(f"pmx force-field directory not found: {pmx_ff_dir}")

    for source in sorted(gmx_top_dir.iterdir()):
        target = gmxlib_dir / source.name
        if target.exists() or target.is_symlink():
            continue
        if source.name in {"residuetypes.dat", "specbond.dat"}:
            continue
        target.symlink_to(source)

    residuetypes_path = gmxlib_dir / "residuetypes.dat"
    _write_residuetypes_overlay(
        source_path=gmx_top_dir / "residuetypes.dat",
        output_path=residuetypes_path,
        additions=_hybrid_residue_types(pmx_ff_dir),
    )

    target_ff = gmxlib_dir / ff_dir_name
    if not target_ff.exists():
        target_ff.symlink_to(pmx_ff_dir)

    pmx_specbond = pmx_mutff_root / "specbond.dat"
    if pmx_specbond.is_file():
        target_specbond = gmxlib_dir / "specbond.dat"
        if target_specbond.exists() or target_specbond.is_symlink():
            target_specbond.unlink()
        target_specbond.symlink_to(pmx_specbond)

    return gmxlib_dir


def water_coordinate_path(gmxlib_dir: Path, water_model: str) -> Path:
    model = water_model.strip().lower()
    filename = WATER_COORDINATE_FILES.get(model)
    if filename is None:
        raise ValueError(f"Unsupported water model for automatic solvent coordinates: {water_model}")
    path = gmxlib_dir / filename
    if not path.is_file():
        raise FileNotFoundError(f"Water coordinate file not found in GMXLIB overlay: {path}")
    return path


def inspect_gro_file(
    path: Path,
    coordinate_abs_threshold_nm: float = 100.0,
    box_abs_threshold_nm: float = 1000.0,
    residue_hydrogen_heavy_distance_threshold_nm: float = 0.25,
) -> dict[str, object]:
    summary: dict[str, object] = {
        "path": str(path),
        "valid": False,
        "reason": "",
        "line_number": None,
        "declared_atom_count": None,
        "atom_count": 0,
        "max_abs_coordinate_nm": None,
        "box_vectors_nm": [],
        "residue_number": None,
        "residue_name": "",
        "atom_name": "",
        "nearest_heavy_atom": "",
        "nearest_heavy_distance_nm": None,
    }

    if not path.is_file():
        summary["reason"] = "missing_file"
        return summary
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        summary["reason"] = f"unreadable_file:{type(exc).__name__}"
        return summary
    if not text.strip():
        summary["reason"] = "empty_file"
        return summary

    lines = text.splitlines()
    if len(lines) < 3:
        summary["reason"] = "too_few_lines"
        return summary

    try:
        declared_atom_count = int(lines[1].strip())
    except ValueError:
        summary["reason"] = "invalid_declared_atom_count"
        summary["line_number"] = 2
        return summary
    if declared_atom_count < 0:
        summary["reason"] = "negative_declared_atom_count"
        summary["line_number"] = 2
        return summary
    summary["declared_atom_count"] = declared_atom_count

    expected_line_count = declared_atom_count + 3
    if len(lines) != expected_line_count:
        summary["reason"] = "atom_count_mismatch"
        summary["atom_count"] = max(len(lines) - 3, 0)
        return summary

    atom_lines = lines[2 : 2 + declared_atom_count]
    summary["atom_count"] = len(atom_lines)
    max_abs_coordinate = 0.0
    protein_residue_atoms: list[list[dict[str, object]]] = []
    current_residue_signature: tuple[int, str, str] | None = None
    for atom_index, atom_line in enumerate(atom_lines, start=1):
        if len(atom_line) < 44:
            summary["reason"] = "atom_line_too_short"
            summary["line_number"] = atom_index + 2
            return summary
        try:
            x = float(atom_line[20:28])
            y = float(atom_line[28:36])
            z = float(atom_line[36:44])
        except ValueError:
            summary["reason"] = "invalid_coordinate"
            summary["line_number"] = atom_index + 2
            return summary
        if not all(math.isfinite(value) for value in (x, y, z)):
            summary["reason"] = "nonfinite_coordinate"
            summary["line_number"] = atom_index + 2
            return summary
        max_abs_coordinate = max(max_abs_coordinate, abs(x), abs(y), abs(z))
        if max_abs_coordinate > coordinate_abs_threshold_nm:
            summary["reason"] = "coordinate_out_of_range"
            summary["line_number"] = atom_index + 2
            summary["max_abs_coordinate_nm"] = max_abs_coordinate
            return summary

        try:
            residue_number = int(atom_line[0:5].strip())
        except ValueError:
            current_residue_signature = None
            continue
        residue_name = atom_line[5:10].strip().upper()
        atom_name = atom_line[10:15].strip().upper()
        normalized_residue_name = RESIDUE_ALIASES.get(residue_name, residue_name)
        if normalized_residue_name not in EXPECTED_HEAVY_ATOMS:
            current_residue_signature = None
            continue
        residue_signature = (residue_number, residue_name, normalized_residue_name)
        if residue_signature != current_residue_signature:
            protein_residue_atoms.append([])
            current_residue_signature = residue_signature
        protein_residue_atoms[-1].append(
            {
                "line_number": atom_index + 2,
                "residue_number": residue_number,
                "residue_name": residue_name,
                "normalized_residue_name": normalized_residue_name,
                "atom_name": atom_name,
                "coordinates_nm": (x, y, z),
            }
        )
    summary["max_abs_coordinate_nm"] = max_abs_coordinate

    box_line = lines[2 + declared_atom_count]
    try:
        box_values = [float(token) for token in box_line.split()]
    except ValueError:
        summary["reason"] = "invalid_box"
        summary["line_number"] = declared_atom_count + 3
        return summary
    if len(box_values) < 3:
        summary["reason"] = "invalid_box"
        summary["line_number"] = declared_atom_count + 3
        return summary
    if not all(math.isfinite(value) for value in box_values):
        summary["reason"] = "nonfinite_box"
        summary["line_number"] = declared_atom_count + 3
        return summary
    if any(value <= 0.0 for value in box_values[:3]) or any(abs(value) > box_abs_threshold_nm for value in box_values):
        summary["reason"] = "box_out_of_range"
        summary["line_number"] = declared_atom_count + 3
        summary["box_vectors_nm"] = box_values
        return summary

    for residue_atoms in protein_residue_atoms:
        heavy_atoms = [item for item in residue_atoms if not str(item["atom_name"]).startswith("H")]
        if not heavy_atoms:
            continue
        for atom in residue_atoms:
            atom_name = str(atom["atom_name"])
            if not atom_name.startswith("H"):
                continue
            hydrogen_coordinates = tuple(float(value) for value in atom["coordinates_nm"])
            nearest_heavy_atom = min(
                heavy_atoms,
                key=lambda item: math.dist(hydrogen_coordinates, tuple(float(value) for value in item["coordinates_nm"])),
            )
            nearest_heavy_distance_nm = math.dist(
                hydrogen_coordinates,
                tuple(float(value) for value in nearest_heavy_atom["coordinates_nm"]),
            )
            if nearest_heavy_distance_nm > residue_hydrogen_heavy_distance_threshold_nm:
                summary["reason"] = "isolated_residue_hydrogen"
                summary["line_number"] = atom["line_number"]
                summary["residue_number"] = atom["residue_number"]
                summary["residue_name"] = atom["residue_name"]
                summary["atom_name"] = atom_name
                summary["nearest_heavy_atom"] = nearest_heavy_atom["atom_name"]
                summary["nearest_heavy_distance_nm"] = nearest_heavy_distance_nm
                return summary

    summary["box_vectors_nm"] = box_values
    summary["valid"] = True
    summary["reason"] = "ok"
    return summary


def gro_file_is_valid(
    path: Path,
    coordinate_abs_threshold_nm: float = 100.0,
    box_abs_threshold_nm: float = 1000.0,
    residue_hydrogen_heavy_distance_threshold_nm: float = 0.25,
) -> bool:
    return bool(
        inspect_gro_file(
            path,
            coordinate_abs_threshold_nm=coordinate_abs_threshold_nm,
            box_abs_threshold_nm=box_abs_threshold_nm,
            residue_hydrogen_heavy_distance_threshold_nm=residue_hydrogen_heavy_distance_threshold_nm,
        ).get("valid")
    )


def deduplicate_standard_topology_includes(topology_path: Path) -> dict[str, int]:
    """Remove duplicate standard water/ion includes from a GROMACS topology.

    `pmx gentop` can emit a duplicate water-model include for monolithic
    single-chain topologies. Keeping the first occurrence matches the standard
    GROMACS topology layout and avoids `moleculetype SOL is redefined` fatals
    during `grompp`.
    """

    original_lines = topology_path.read_text(encoding="utf-8").splitlines()
    deduped_lines: list[str] = []
    seen_basenames: set[str] = set()
    removed_counts: dict[str, int] = {}

    for raw_line in original_lines:
        include_target = _parse_include_target(raw_line)
        if include_target is None:
            deduped_lines.append(raw_line)
            continue
        basename = Path(include_target).name.lower()
        if basename not in DEDUP_INCLUDE_BASENAMES:
            deduped_lines.append(raw_line)
            continue
        if basename in seen_basenames:
            removed_counts[basename] = removed_counts.get(basename, 0) + 1
            continue
        seen_basenames.add(basename)
        deduped_lines.append(raw_line)

    if deduped_lines != original_lines:
        topology_path.write_text("\n".join(deduped_lines) + "\n", encoding="utf-8")
    return removed_counts


def generate_hybrid_topology(
    topology_path: Path,
    output_path: Path,
    force_field: str,
    mutated_chain_ids: list[str] | tuple[str, ...],
    pmx_command: list[str] | None = None,
    restore_summary_path: Path | None = None,
    pdbfixer_summary_path: Path | None = None,
    allow_reuse_existing: bool = False,
) -> dict[str, object]:
    """Generate a pmx hybrid topology, preferring per-chain conversion when possible.

    For protein-complex topologies produced by ``pdb2gmx``, only the chain itp
    files that actually contain hybrid residues need B-state filling. Running
    ``pmx gentop`` recursively over every included chain is correct but
    unnecessarily expensive for large antibody-antigen systems. This helper
    converts only the mutated chain itps and rewrites the top-level include file
    to point at the generated ``pmx_*.itp`` files.

    If the topology does not expose chain itp includes, or if the mutated chain
    itps cannot be resolved, it falls back to the standard recursive
    ``pmx gentop`` path.
    """

    topology_path = Path(topology_path)
    output_path = Path(output_path)
    normalized_mutated_chains = {
        str(chain_id).strip().upper()
        for chain_id in mutated_chain_ids
        if str(chain_id).strip()
    }
    pmx_cli = list(pmx_command) if pmx_command else [sys.executable, "-m", "pmx.scripts.cli"]
    summary: dict[str, object] = {
        "mode": "recursive_fallback",
        "topology_path": str(topology_path),
        "output_path": str(output_path),
        "mutated_chain_ids": sorted(normalized_mutated_chains),
        "converted_itps": [],
        "reused_itps": [],
        "fallback_reason": "",
    }

    lines = topology_path.read_text(encoding="utf-8").splitlines()
    mutated_targets: list[dict[str, object]] = []
    for raw_line in lines:
        include_target = _parse_include_target(raw_line)
        if include_target is None:
            continue
        basename = Path(include_target).name
        match = CHAIN_ITP_BASENAME_RE.match(basename)
        if match is None:
            continue
        chain_id = str(match.group("chain") or "").strip().upper()
        if normalized_mutated_chains and chain_id not in normalized_mutated_chains:
            summary["reused_itps"].append(include_target)
            continue
        input_itp = _resolve_include_path(topology_path.parent, include_target)
        output_include_target = str(Path(include_target).with_name(f"pmx_{basename}"))
        output_itp = _resolve_include_path(topology_path.parent, output_include_target)
        mutated_targets.append(
            {
                "chain_id": chain_id,
                "include_target": include_target,
                "output_include_target": output_include_target,
                "input_itp": input_itp,
                "output_itp": output_itp,
            }
        )

    if not mutated_targets:
        if allow_reuse_existing and _topology_repairs_were_noops(restore_summary_path, pdbfixer_summary_path):
            if _hybrid_topology_matches_source(topology_path, output_path):
                summary["mode"] = "reused_existing"
                summary["fallback_reason"] = ""
                return summary
        summary["fallback_reason"] = "no_mutated_chain_itps_detected"
        _run_pmx_gentop(topology_path, output_path, force_field, pmx_cli, recursive=True)
        _ensure_hybrid_topology_matches_source(topology_path, output_path)
        return summary

    for target in mutated_targets:
        input_itp = target["input_itp"]
        if not isinstance(input_itp, Path) or not input_itp.is_file():
            summary["fallback_reason"] = f"missing_mutated_chain_itp:{target['include_target']}"
            _run_pmx_gentop(topology_path, output_path, force_field, pmx_cli, recursive=True)
            _ensure_hybrid_topology_matches_source(topology_path, output_path)
            return summary

    replacements = {
        str(target["include_target"]): str(target["output_include_target"])
        for target in mutated_targets
    }
    if allow_reuse_existing and _topology_repairs_were_noops(restore_summary_path, pdbfixer_summary_path):
        if all(
            isinstance(target["input_itp"], Path)
            and isinstance(target["output_itp"], Path)
            and _hybrid_topology_matches_source(target["input_itp"], target["output_itp"])
            for target in mutated_targets
        ):
            _write_rewritten_topology(lines, output_path, replacements)
            summary["mode"] = "reused_existing"
            summary["fallback_reason"] = ""
            summary["converted_itps"] = [
                {
                    "chain_id": target["chain_id"],
                    "input": str(target["input_itp"]),
                    "output": str(target["output_itp"]),
                }
                for target in mutated_targets
            ]
            return summary

    for target in mutated_targets:
        output_itp = target["output_itp"]
        if isinstance(output_itp, Path):
            output_itp.parent.mkdir(parents=True, exist_ok=True)
        _run_pmx_gentop(
            target["input_itp"],
            target["output_itp"],
            force_field,
            pmx_cli,
            recursive=False,
        )
        _ensure_hybrid_topology_matches_source(target["input_itp"], target["output_itp"])
        summary["converted_itps"].append(
            {
                "chain_id": target["chain_id"],
                "input": str(target["input_itp"]),
                "output": str(target["output_itp"]),
            }
        )

    _write_rewritten_topology(lines, output_path, replacements)
    summary["mode"] = "per_chain"
    summary["fallback_reason"] = ""
    return summary


def materialize_staged_equilibration_restraints(
    topology_path: Path,
    *,
    heavy_force_constant: float = 1000.0,
    backbone_force_constant: float = 250.0,
) -> dict[str, object]:
    topology_path = Path(topology_path).resolve()
    summary: dict[str, object] = {
        "topology_path": str(topology_path),
        "modified_topology_count": 0,
        "generated_file_count": 0,
        "files": [],
    }
    if not topology_path.is_file():
        summary["status"] = "missing_topology"
        return summary

    candidate_paths: list[Path] = [topology_path]
    seen: set[Path] = {topology_path}
    for raw_line in topology_path.read_text(encoding="utf-8").splitlines():
        include_target = _parse_include_target(raw_line)
        if include_target is None or not include_target.lower().endswith(".itp"):
            continue
        include_path = _resolve_include_path(topology_path.parent, include_target)
        if include_path.parent != topology_path.parent or include_path in seen or not include_path.is_file():
            continue
        candidate_paths.append(include_path)
        seen.add(include_path)

    for candidate in candidate_paths:
        payload = _materialize_staged_posres_for_topology_file(
            candidate,
            heavy_force_constant=heavy_force_constant,
            backbone_force_constant=backbone_force_constant,
        )
        if payload is None:
            continue
        summary["files"].append(payload)
        if payload.get("topology_modified"):
            summary["modified_topology_count"] = int(summary["modified_topology_count"]) + 1
        generated_files = payload.get("generated_posre_files", [])
        summary["generated_file_count"] = int(summary["generated_file_count"]) + len(generated_files)

    summary["status"] = "ok" if summary["files"] else "no_position_restraint_targets"
    return summary


def _hybrid_residue_types(force_field_dir: Path) -> dict[str, str]:
    additions: dict[str, str] = {}
    for filename, residue_type in (
        ("mutres.rtp", "Protein"),
        ("mutres_dna.rtp", "DNA"),
        ("mutres_rna.rtp", "RNA"),
    ):
        rtp_path = force_field_dir / filename
        if not rtp_path.is_file():
            continue
        for residue_name in _parse_rtp_residue_names(rtp_path):
            additions[residue_name] = residue_type
    return additions


def _parse_rtp_residue_names(path: Path) -> list[str]:
    residue_names: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split(";", 1)[0].strip()
        if not line.startswith("[") or "]" not in line:
            continue
        section = line[1 : line.index("]")].strip()
        if not section or section.lower() in RTP_SECTION_EXCLUSIONS:
            continue
        residue_names.append(section)
    return residue_names


def _write_residuetypes_overlay(source_path: Path, output_path: Path, additions: dict[str, str]) -> None:
    if output_path.is_symlink() or output_path.is_file():
        output_path.unlink()
    existing_lines = source_path.read_text(encoding="utf-8").splitlines()
    existing_names = {
        line.split()[0]
        for line in existing_lines
        if line.strip() and not line.lstrip().startswith(";")
    }
    merged_lines = list(existing_lines)
    if merged_lines and merged_lines[-1].strip():
        merged_lines.append("")
    for residue_name in sorted(additions):
        if residue_name not in existing_names:
            merged_lines.append(f"{residue_name}\t{additions[residue_name]}")
    output_path.write_text("\n".join(merged_lines) + "\n", encoding="utf-8")


def _parse_include_target(raw_line: str) -> str | None:
    line = raw_line.strip()
    if not line.startswith("#include"):
        return None
    remainder = line[len("#include") :].strip()
    if not remainder:
        return None
    opener = remainder[0]
    closer = '"' if opener == '"' else ">" if opener == "<" else ""
    if not closer:
        return None
    end_index = remainder.find(closer, 1)
    if end_index == -1:
        return None
    return remainder[1:end_index]


def _materialize_staged_posres_for_topology_file(
    topology_path: Path,
    *,
    heavy_force_constant: float,
    backbone_force_constant: float,
) -> dict[str, object] | None:
    lines = topology_path.read_text(encoding="utf-8").splitlines()
    posres_include_target = _posres_include_target(lines)
    if posres_include_target is None:
        return None

    atom_entries = _topology_atoms(lines)
    if not atom_entries:
        return None

    heavy_indices = [atom_index for atom_index, atom_name in atom_entries if not atom_name.startswith("H")]
    backbone_indices = [atom_index for atom_index, atom_name in atom_entries if atom_name in PROTEIN_BACKBONE_ATOM_NAMES]
    if not heavy_indices or not backbone_indices:
        return None

    include_path = Path(posres_include_target)
    suffix = include_path.suffix or ".itp"
    stem = include_path.stem or "posre"
    heavy_include_target = str(include_path.with_name(f"{stem}_stage_heavy{suffix}"))
    backbone_include_target = str(include_path.with_name(f"{stem}_stage_backbone{suffix}"))
    heavy_path = _resolve_include_path(topology_path.parent, heavy_include_target)
    backbone_path = _resolve_include_path(topology_path.parent, backbone_include_target)
    heavy_path.parent.mkdir(parents=True, exist_ok=True)
    backbone_path.parent.mkdir(parents=True, exist_ok=True)
    heavy_path.write_text(_render_position_restraints(heavy_indices, heavy_force_constant), encoding="utf-8")
    backbone_path.write_text(_render_position_restraints(backbone_indices, backbone_force_constant), encoding="utf-8")

    topology_modified = False
    joined_lines = "\n".join(lines)
    if "#ifdef POSRES_STAGE_HEAVY" not in joined_lines and "#ifdef POSRES_STAGE_BACKBONE" not in joined_lines:
        rewritten_lines = _insert_staged_posres_blocks(
            lines,
            heavy_include_target=heavy_include_target,
            backbone_include_target=backbone_include_target,
        )
        topology_path.write_text("\n".join(rewritten_lines) + "\n", encoding="utf-8")
        topology_modified = rewritten_lines != lines

    return {
        "topology_path": str(topology_path),
        "position_restraint_include": posres_include_target,
        "heavy_atom_restraint_count": len(heavy_indices),
        "backbone_atom_restraint_count": len(backbone_indices),
        "generated_posre_files": [str(heavy_path), str(backbone_path)],
        "heavy_include_target": heavy_include_target,
        "backbone_include_target": backbone_include_target,
        "topology_modified": topology_modified,
    }


def _posres_include_target(lines: list[str]) -> str | None:
    in_posres_block = False
    for raw_line in lines:
        line = raw_line.strip()
        if line == "#ifdef POSRES":
            in_posres_block = True
            continue
        if not in_posres_block:
            continue
        if line == "#endif":
            return None
        include_target = _parse_include_target(raw_line)
        if include_target is not None:
            return include_target
    return None


def _topology_atoms(lines: list[str]) -> list[tuple[int, str]]:
    atom_entries: list[tuple[int, str]] = []
    in_atoms = False
    for raw_line in lines:
        line = raw_line.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            if "]" not in line:
                return []
            section = line[1 : line.index("]")].strip().lower()
            if section == "atoms":
                in_atoms = True
                continue
            if in_atoms:
                break
            continue
        if not in_atoms:
            continue
        fields = line.split()
        if len(fields) < 5:
            return []
        try:
            atom_index = int(fields[0])
        except ValueError:
            return []
        atom_entries.append((atom_index, fields[4].strip().upper()))
    return atom_entries


def _render_position_restraints(atom_indices: list[int], force_constant: float) -> str:
    formatted_force = f"{float(force_constant):g}"
    lines = [
        "[ position_restraints ]",
        "; ai  funct  fcx  fcy  fcz",
    ]
    for atom_index in atom_indices:
        lines.append(
            f"{atom_index:6d}{1:7d}{formatted_force:>8}{formatted_force:>8}{formatted_force:>8}"
        )
    return "\n".join(lines) + "\n"


def _insert_staged_posres_blocks(
    lines: list[str],
    *,
    heavy_include_target: str,
    backbone_include_target: str,
) -> list[str]:
    insert_at = len(lines)
    for index, raw_line in enumerate(lines):
        if raw_line.strip() != "#ifdef POSRES":
            continue
        for end_index in range(index + 1, len(lines)):
            if lines[end_index].strip() == "#endif":
                insert_at = end_index + 1
                return [
                    *lines[:insert_at],
                    "",
                    "#ifdef POSRES_STAGE_HEAVY",
                    f'#include "{heavy_include_target}"',
                    "#endif",
                    "",
                    "#ifdef POSRES_STAGE_BACKBONE",
                    f'#include "{backbone_include_target}"',
                    "#endif",
                    *lines[insert_at:],
                ]
        break
    return lines


def _nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _topology_atom_signatures(path: Path) -> tuple[tuple[str, str, str], ...] | None:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None

    signatures: list[tuple[str, str, str]] = []
    in_atoms = False
    for raw_line in lines:
        line = raw_line.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            if "]" not in line:
                return None
            section = line[1 : line.index("]")].strip().lower()
            if section == "atoms":
                in_atoms = True
                continue
            if in_atoms:
                break
            continue
        if not in_atoms:
            continue
        fields = line.split()
        if len(fields) < 5:
            return None
        signatures.append((fields[2], fields[3], fields[4]))
    if not signatures:
        return None
    return tuple(signatures)


def _hybrid_topology_matches_source(source_path: Path, output_path: Path) -> bool:
    if not _nonempty_file(output_path):
        return False
    source_signatures = _topology_atom_signatures(source_path)
    output_signatures = _topology_atom_signatures(output_path)
    return source_signatures is not None and source_signatures == output_signatures


def _ensure_hybrid_topology_matches_source(source_path: Path, output_path: Path) -> None:
    if _hybrid_topology_matches_source(source_path, output_path):
        return
    source_signatures = _topology_atom_signatures(source_path) or ()
    output_signatures = _topology_atom_signatures(output_path) or ()
    missing = len(set(source_signatures) - set(output_signatures))
    extra = len(set(output_signatures) - set(source_signatures))
    raise RuntimeError(
        "Generated hybrid topology does not match the source topology atoms: "
        f"source={source_path} output={output_path} "
        f"source_atom_count={len(source_signatures)} output_atom_count={len(output_signatures)} "
        f"missing_atoms={missing} extra_atoms={extra}"
    )


def _read_optional_summary(path: Path | None) -> dict[str, object]:
    if path is None or not path.is_file():
        return {}
    try:
        import json

        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _topology_repairs_were_noops(
    restore_summary_path: Path | None,
    pdbfixer_summary_path: Path | None,
) -> bool:
    restore_summary = _read_optional_summary(restore_summary_path)
    pdbfixer_summary = _read_optional_summary(pdbfixer_summary_path)
    if not restore_summary and not pdbfixer_summary:
        return False

    restore_changed = any(
        int(restore_summary.get(key) or 0) > 0
        for key in ("attempted_residue_count", "restored_residue_count", "unresolved_residue_count")
    )
    pdbfixer_changed = any(
        int(pdbfixer_summary.get(key) or 0) > 0
        for key in ("trigger_residue_count", "blocking_residue_count", "remaining_incomplete_standard_residue_count")
    )
    return not restore_changed and not pdbfixer_changed


def _write_rewritten_topology(lines: list[str], output_path: Path, replacements: dict[str, str]) -> None:
    rewritten_lines: list[str] = []
    for raw_line in lines:
        include_target = _parse_include_target(raw_line)
        replacement = replacements.get(include_target or "")
        if include_target is not None and replacement is not None:
            rewritten_lines.append(raw_line.replace(include_target, replacement, 1))
            continue
        rewritten_lines.append(raw_line)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(rewritten_lines) + "\n", encoding="utf-8")


def _resolve_include_path(base_dir: Path, include_target: str) -> Path:
    include_path = Path(include_target)
    if include_path.is_absolute():
        return include_path
    return (base_dir / include_path).resolve()


def _run_pmx_gentop(
    topology_path: Path,
    output_path: Path,
    force_field: str,
    pmx_command: list[str],
    *,
    recursive: bool,
) -> None:
    command = [*pmx_command, "gentop", "-p", str(topology_path), "-o", str(output_path), "-ff", force_field]
    if not recursive:
        command.append("--norecursive")
    subprocess.run(
        command,
        check=True,
        cwd=str(topology_path.parent.resolve()),
    )


_HYBRID_RESNAME_RE = re.compile(r"^[A-Z]{1,3}2[A-Z]{1,2}$")
_CHARGE_INTEGER_TOLERANCE = 0.05


def _parse_itp_hybrid_residues(itp_path: Path) -> dict[int, dict[str, object]]:
    """Collect per-residue atom counts and A/B-state charge sums for hybrid
    residues (pmx names like Q2A, Z2A) from a pmx topology itp file."""
    residues: dict[int, dict[str, object]] = {}
    in_atoms = False
    for line in itp_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("["):
            in_atoms = stripped.lower().startswith("[ atoms ]")
            continue
        if not in_atoms or not stripped or stripped.startswith(";"):
            continue
        parts = stripped.split()
        if len(parts) < 8 or not parts[0].lstrip("-").isdigit():
            continue
        resname = parts[3]
        if not _HYBRID_RESNAME_RE.match(resname):
            continue
        resnr = int(parts[2])
        entry = residues.setdefault(
            resnr,
            {"resnr": resnr, "resname": resname, "atom_count": 0, "charge_a": 0.0, "charge_b": 0.0},
        )
        entry["atom_count"] = int(entry["atom_count"]) + 1
        charge_a = float(parts[6])
        entry["charge_a"] = float(entry["charge_a"]) + charge_a
        # B-state columns (typeB chargeB massB) exist only for morphing atoms.
        charge_b = float(parts[9]) if len(parts) >= 11 else charge_a
        entry["charge_b"] = float(entry["charge_b"]) + charge_b
    for entry in residues.values():
        entry["charge_a"] = round(float(entry["charge_a"]), 4)
        entry["charge_b"] = round(float(entry["charge_b"]), 4)
    return residues


def _mutant_pdb_hybrid_atom_counts(mutant_pdb_path: Path) -> dict[tuple[str, int], int]:
    """Count atoms per hybrid residue (chain, resseq) in a pmx mutant PDB."""
    counts: dict[tuple[str, int], int] = {}
    if not mutant_pdb_path.is_file():
        return counts
    for line in mutant_pdb_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM"):
            continue
        resname = line[17:20].strip()
        if not _HYBRID_RESNAME_RE.match(resname):
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        key = (line[21].strip(), resseq)
        counts[key] = counts.get(key, 0) + 1
    return counts


def _charge_is_integer(charge: float) -> bool:
    return abs(charge - round(charge)) <= _CHARGE_INTEGER_TOLERANCE


def validate_hybrid_topology_integrity(pmx_dir: Path) -> dict[str, object]:
    """Validate hybrid-residue chemistry in a leg's pmx directory.

    Guards against heavy-atom-only hybrids (e.g. the historical pdb2gmx -ignh
    defect): every hybrid residue in pmx_topol_*.itp must keep all atoms that
    pmx wrote into mutant.pdb, and its A/B-state charges must each sum to an
    integer (0 for charge-conserving mutations, +/-1 for protonation variants).
    """
    pmx_dir = Path(pmx_dir)
    mutant_counts = _mutant_pdb_hybrid_atom_counts(pmx_dir / "mutant.pdb")
    chain_from_itp = re.compile(r"pmx_topol_Protein_chain_(?P<chain>.+)\.itp$")

    residues: list[dict[str, object]] = []
    issues: list[str] = []
    for itp_path in sorted(pmx_dir.glob("pmx_topol_Protein_chain_*.itp")):
        match = chain_from_itp.match(itp_path.name)
        chain = match.group("chain") if match else ""
        for resnr, entry in sorted(_parse_itp_hybrid_residues(itp_path).items()):
            expected = mutant_counts.get((chain, resnr))
            record = dict(entry)
            record["chain_id"] = chain
            record["expected_atom_count"] = expected
            record["charge_a_integer"] = _charge_is_integer(float(entry["charge_a"]))
            record["charge_b_integer"] = _charge_is_integer(float(entry["charge_b"]))
            record["atoms_complete"] = expected is None or int(entry["atom_count"]) == expected
            residues.append(record)
            label = f"{chain}:{entry['resname']}{resnr}"
            if not record["atoms_complete"]:
                issues.append(
                    f"{label}: topology has {entry['atom_count']} atoms but mutant.pdb has {expected}"
                )
            if not record["charge_a_integer"]:
                issues.append(f"{label}: state A charge sum {entry['charge_a']} is not an integer")
            if not record["charge_b_integer"]:
                issues.append(f"{label}: state B charge sum {entry['charge_b']} is not an integer")

    return {
        "checked": bool(residues),
        "ok": not issues,
        "residues": residues,
        "issues": issues,
    }
