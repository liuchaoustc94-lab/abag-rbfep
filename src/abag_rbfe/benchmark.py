"""AB-Bind source curation into RBFE-ready benchmark layers."""

from __future__ import annotations

from collections import Counter, deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
import csv
from math import ceil, exp, log1p, sqrt
import re
from typing import Any, Callable

from abag_pmx.mutations import (
    build_mutation_group,
    is_charge_conserving,
    parse_mutation_group_tokens,
    parse_mutation_site_dict,
)
from abag_rbfe.execution import discover_visible_gpu_devices
from abag_rbfe.io_utils import ensure_dir, read_csv_rows, read_json, read_yaml, utc_now, write_csv_rows, write_json, write_yaml
from abag_rbfe.models import MutationSite
from abag_rbfe.paths import ProjectPaths
from abag_rbfe.planning import build_batch_plan, slugify
from abag_rbfe.reporting import write_batch_summary
from abag_rbfe.stages import resume_job, run_job

AB_BIND_RAW_MUTATION_RE = re.compile(
    r"^(?P<chain>[A-Za-z0-9]):(?P<wt>[A-Z])(?P<resseq>-?\d+)(?P<icode>[A-Za-z]?)(?P<mut>[A-Z])$"
)


def _default_ab_bind_annotations_path() -> Path:
    return Path(__file__).resolve().parents[2] / "benchmarks" / "ab_bind" / "source" / "ab_bind_complex_annotations.csv"


def _default_ab_bind_split_path(benchmark_root: Path, spec_name: str) -> Path:
    return benchmark_root / "splits" / f"ab_bind_rbfe_{spec_name}_split_v1.yml"


def _normalize_gpu_devices(values: list[str] | None) -> tuple[str, ...]:
    if values:
        devices = tuple(item.strip() for item in values if item.strip())
        if devices:
            return devices
    return discover_visible_gpu_devices()


def _read_csv_rows_with_fallback(path: Path, encodings: tuple[str, ...] = ("utf-8", "latin-1")) -> list[dict[str, str]]:
    last_error: UnicodeDecodeError | None = None
    for encoding in encodings:
        try:
            with path.open("r", encoding=encoding, newline="") as handle:
                reader = csv.DictReader(handle)
                rows = []
                for row in reader:
                    cleaned = {}
                    for key, value in row.items():
                        if key is None:
                            continue
                        cleaned[key.lstrip("\ufeff").strip()] = (value or "").strip()
                    rows.append(cleaned)
                return rows
        except UnicodeDecodeError as exc:
            last_error = exc
    if last_error is not None:
        raise last_error
    return []


def _find_value(row: dict[str, str], candidates: tuple[str, ...], default: str = "") -> str:
    lower = {key.lower(): value for key, value in row.items()}
    for candidate in candidates:
        if candidate.lower() in lower and lower[candidate.lower()] != "":
            return lower[candidate.lower()]
    return default


def _as_bool(value: str, default: bool = False) -> bool:
    if value == "":
        return default
    return value.strip().lower() in {"1", "true", "yes", "y"}


def _coerce_bool(value: Any, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return _as_bool(str(value), default=default)


def _split_chain_spec(value: str) -> tuple[str, ...]:
    compact = value.replace(",", "").replace(" ", "").strip().upper()
    return tuple(char for char in compact if char)


def _looks_like_ab_bind_raw_schema(rows: list[dict[str, str]]) -> bool:
    if not rows:
        return False
    keys = {key.lower() for key in rows[0].keys()}
    return "#pdb" in keys and "mutation" in keys and "ddg(kcal/mol)" in keys


@dataclass(frozen=True)
class ComplexAnnotation:
    complex_id: str
    complex_class: str
    structure_source: str
    antibody_chains: tuple[str, ...]
    antigen_chains: tuple[str, ...]
    structure_mappable: bool
    stable_run_candidate: bool


def _load_ab_bind_annotations(path: Path) -> dict[str, ComplexAnnotation]:
    annotations = {}
    for row in read_csv_rows(path):
        complex_id = row["complex_id"].strip()
        annotations[complex_id] = ComplexAnnotation(
            complex_id=complex_id,
            complex_class=row["complex_class"].strip(),
            structure_source=row["structure_source"].strip(),
            antibody_chains=_split_chain_spec(row.get("antibody_chains", "")),
            antigen_chains=_split_chain_spec(row.get("antigen_chains", "")),
            structure_mappable=_as_bool(row.get("structure_mappable", "true"), default=True),
            stable_run_candidate=_as_bool(row.get("stable_run_candidate", "true"), default=True),
        )
    return annotations


def _parse_ab_bind_raw_mutation_sites(raw_mutation: str, annotation: ComplexAnnotation | None) -> tuple[list[MutationSite], list[str]]:
    exclusions: list[str] = []
    sites: list[MutationSite] = []
    tokens = [token.strip() for token in raw_mutation.split(",") if token.strip()]
    for token in tokens:
        match = AB_BIND_RAW_MUTATION_RE.match(token)
        if not match:
            exclusions.append("invalid_mutation_token")
            continue
        chain_id = match.group("chain").upper()
        if annotation is None:
            exclusions.append("missing_complex_annotation")
            continue
        if chain_id in annotation.antibody_chains:
            entity_side = "antibody"
        elif chain_id in annotation.antigen_chains:
            entity_side = "antigen"
        else:
            exclusions.append("unknown_chain_role")
            continue
        sites.append(
            MutationSite(
                chain_id=chain_id,
                resseq=int(match.group("resseq")),
                icode=(match.group("icode") or "").upper(),
                wt=match.group("wt").upper(),
                mut=match.group("mut").upper(),
                entity_side=entity_side,
            )
        )
    return sites, exclusions


def _canonicalize_legacy_row(index: int, row: dict[str, str]) -> tuple[dict[str, str], list[str]]:
    exclusions: list[str] = []
    row_id = _find_value(row, ("row_id", "id", "dataset_id"), str(index))
    complex_id = _find_value(row, ("complex_id", "pdb", "pdb_id", "complex"), f"row_{index}")
    complex_class = _find_value(row, ("complex_class", "complex_type", "dataset_class"), "")
    structure_source = _find_value(row, ("structure_source", "structure_kind"), "experimental")
    ddg_value = _find_value(row, ("ddg_kcal_mol", "ddg", "delta_delta_g"), "")
    mapped = _as_bool(_find_value(row, ("structure_mappable", "mapped", "is_mapped"), "true"), default=True)
    stable = _as_bool(_find_value(row, ("stable_run_candidate", "stable_candidate", "is_stable_candidate"), "true"), default=True)

    if complex_class and complex_class.lower() != "antibody-antigen":
        exclusions.append("non_antibody_antigen")
    if ddg_value == "":
        exclusions.append("missing_ddg")
    if not mapped:
        exclusions.append("not_structure_mappable")
    if not stable:
        exclusions.append("not_stable_candidate")
    if structure_source.lower() != "experimental":
        exclusions.append("non_experimental_structure")

    if "mutation_tokens" in {field.lower() for field in row.keys()}:
        tokens = _find_value(row, ("mutation_tokens",))
        sites = parse_mutation_group_tokens(tokens)
        mutation_tokens = ";".join(site.token() for site in sites)
        raw_mutation = tokens
    else:
        sites = [parse_mutation_site_dict(row)]
        mutation_tokens = sites[0].token()
        raw_mutation = mutation_tokens

    group_id = _find_value(row, ("mutation_group_id", "job_id"), f"group_{index}")
    mutation_count = len(sites)
    if mutation_count == 0:
        exclusions.append("missing_mutation_definition")
    if mutation_count > 2:
        exclusions.append("unsupported_mutation_count")
    sides = {site.entity_side for site in sites}
    if mutation_count == 2 and len(sides) != 1:
        exclusions.append("double_cross_side")
    charge_conserving = all(is_charge_conserving(site.wt, site.mut) for site in sites)
    if mutation_count == 1 and not charge_conserving:
        exclusions.append("charge_changing_v1")
    if mutation_count == 2 and not charge_conserving:
        exclusions.append("charge_changing_v2")

    try:
        mutation_group = build_mutation_group(
            mutation_group_id=group_id,
            sites=sites,
            allow_double_same_side=True,
            allow_charge_change=True,
        )
    except ValueError:
        entity_side = sites[0].entity_side if len(sides) == 1 and sites else "mixed"
        signature = mutation_tokens
    else:
        entity_side = mutation_group.entity_side
        signature = mutation_group.signature()

    payload = {
        "row_id": row_id,
        "complex_id": complex_id,
        "mutation_group_id": group_id,
        "mutation_count": str(mutation_count),
        "entity_side": entity_side,
        "ddg_kcal_mol": ddg_value,
        "structure_source": structure_source,
        "signature": signature,
        "complex_class": complex_class or "antibody-antigen",
        "source_mutation": raw_mutation,
        "mutation_tokens": mutation_tokens,
        "partners": _find_value(row, ("partners", "partners(a_b)")),
        "protein_1": _find_value(row, ("protein_1", "protein-1")),
        "protein_2": _find_value(row, ("protein_2", "protein-2")),
        "antibody_chains": _find_value(row, ("antibody_chains",)),
        "antigen_chains": _find_value(row, ("antigen_chains",)),
        "structure_mappable": str(mapped).lower(),
        "stable_run_candidate": str(stable).lower(),
    }
    return payload, exclusions


def _canonicalize_ab_bind_raw_row(
    index: int,
    row: dict[str, str],
    annotations: dict[str, ComplexAnnotation],
) -> tuple[dict[str, str], list[str]]:
    exclusions: list[str] = []
    complex_id = row["#PDB"].strip()
    annotation = annotations.get(complex_id)
    if annotation is None:
        exclusions.append("missing_complex_annotation")
        complex_class = "unknown"
        structure_source = "homology_model" if complex_id.startswith("HM_") else "experimental"
        antibody_chains: tuple[str, ...] = ()
        antigen_chains: tuple[str, ...] = ()
        mapped = False
        stable = False
    else:
        complex_class = annotation.complex_class
        structure_source = annotation.structure_source
        antibody_chains = annotation.antibody_chains
        antigen_chains = annotation.antigen_chains
        mapped = annotation.structure_mappable
        stable = annotation.stable_run_candidate

    raw_mutation = row["Mutation"].strip()
    raw_tokens = [token.strip() for token in raw_mutation.split(",") if token.strip()]
    mutation_count = len(raw_tokens)
    ddg_value = row.get("ddG(kcal/mol)", "").strip()

    if complex_class.lower() != "antibody-antigen":
        exclusions.append("non_antibody_antigen")
    if ddg_value == "":
        exclusions.append("missing_ddg")
    if not mapped:
        exclusions.append("not_structure_mappable")
    if not stable:
        exclusions.append("not_stable_candidate")
    if structure_source.lower() != "experimental":
        exclusions.append("non_experimental_structure")
    if mutation_count == 0:
        exclusions.append("missing_mutation_definition")
    if mutation_count > 2:
        exclusions.append("unsupported_mutation_count")

    sites: list[MutationSite] = []
    if annotation is not None and complex_class.lower() == "antibody-antigen":
        sites, site_exclusions = _parse_ab_bind_raw_mutation_sites(raw_mutation, annotation)
        exclusions.extend(site_exclusions)
        if not antibody_chains or not antigen_chains:
            exclusions.append("missing_chain_role_annotation")

    entity_side = "other"
    signature = raw_mutation
    mutation_tokens = raw_mutation
    if len(sites) == mutation_count and mutation_count in {1, 2}:
        sides = {site.entity_side for site in sites}
        if mutation_count == 2 and len(sides) != 1:
            exclusions.append("double_cross_side")
        charge_conserving = all(is_charge_conserving(site.wt, site.mut) for site in sites)
        if mutation_count == 1 and not charge_conserving:
            exclusions.append("charge_changing_v1")
        if mutation_count == 2 and not charge_conserving:
            exclusions.append("charge_changing_v2")
        try:
            mutation_group = build_mutation_group(
                mutation_group_id=f"{complex_id.lower()}_{index:04d}",
                sites=sites,
                allow_double_same_side=True,
                allow_charge_change=True,
            )
        except ValueError:
            entity_side = sites[0].entity_side if len(sides) == 1 else "mixed"
            signature = ";".join(site.token() for site in sites)
            mutation_tokens = signature
        else:
            entity_side = mutation_group.entity_side
            signature = mutation_group.signature()
            mutation_tokens = ";".join(site.token() for site in sites)
    elif sites:
        sides = {site.entity_side for site in sites}
        entity_side = sites[0].entity_side if len(sides) == 1 else "mixed"
        signature = ";".join(site.token() for site in sites)
        mutation_tokens = signature

    payload = {
        "row_id": f"{complex_id.lower()}_{index:04d}",
        "complex_id": complex_id,
        "mutation_group_id": f"{complex_id.lower()}_{index:04d}",
        "mutation_count": str(mutation_count),
        "entity_side": entity_side,
        "ddg_kcal_mol": ddg_value,
        "structure_source": structure_source,
        "signature": signature,
        "complex_class": complex_class,
        "source_mutation": raw_mutation,
        "mutation_tokens": mutation_tokens,
        "partners": row.get("Partners(A_B)", "").strip(),
        "protein_1": row.get("Protein-1", "").strip(),
        "protein_2": row.get("Protein-2", "").strip(),
        "antibody_chains": "".join(antibody_chains),
        "antigen_chains": "".join(antigen_chains),
        "structure_mappable": str(mapped).lower(),
        "stable_run_candidate": str(stable).lower(),
    }
    return payload, exclusions


def curate_ab_bind(source_csv: Path, output_dir: Path, annotations_path: Path | None = None) -> dict:
    source_rows_raw = _read_csv_rows_with_fallback(source_csv)
    annotations = {}
    if _looks_like_ab_bind_raw_schema(source_rows_raw):
        annotations_path = annotations_path or _default_ab_bind_annotations_path()
        annotations = _load_ab_bind_annotations(annotations_path)

    source_rows = []
    core_v1_rows = []
    core_v2_rows = []
    exclusion_counts: Counter[str] = Counter()

    for index, row in enumerate(source_rows_raw, start=1):
        if _looks_like_ab_bind_raw_schema(source_rows_raw):
            payload, exclusions = _canonicalize_ab_bind_raw_row(index, row, annotations)
        else:
            payload, exclusions = _canonicalize_legacy_row(index, row)

        for code in exclusions:
            exclusion_counts[code] += 1

        core_v1_eligible = len(exclusions) == 0 and int(payload["mutation_count"]) == 1
        core_v2_eligible = len(exclusions) == 0 and int(payload["mutation_count"]) in {1, 2}
        source_row = {
            **payload,
            "core_v1_eligible": core_v1_eligible,
            "core_v2_eligible": core_v2_eligible,
            "exclusion_codes": ";".join(exclusions),
        }
        source_rows.append(source_row)
        if core_v1_eligible:
            core_v1_rows.append(source_row)
        if core_v2_eligible:
            core_v2_rows.append(source_row)

    output_dir.mkdir(parents=True, exist_ok=True)
    curated_dir = output_dir / "curated"
    manifests_dir = output_dir / "manifests"
    curated_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "row_id",
        "complex_id",
        "mutation_group_id",
        "mutation_count",
        "entity_side",
        "ddg_kcal_mol",
        "structure_source",
        "signature",
        "complex_class",
        "source_mutation",
        "mutation_tokens",
        "partners",
        "protein_1",
        "protein_2",
        "antibody_chains",
        "antigen_chains",
        "structure_mappable",
        "stable_run_candidate",
        "core_v1_eligible",
        "core_v2_eligible",
        "exclusion_codes",
    ]
    write_csv_rows(curated_dir / "ab_bind_source_registered.csv", source_rows, fieldnames)
    write_csv_rows(curated_dir / "ab_bind_rbfe_core_v1.csv", core_v1_rows, fieldnames)
    write_csv_rows(curated_dir / "ab_bind_rbfe_core_v2.csv", core_v2_rows, fieldnames)

    manifests = {
        "source": {
            "name": "AB-Bind-Source",
            "source_csv": str(source_csv),
            "annotations_csv": str(annotations_path) if annotations_path else "",
            "row_count": len(source_rows),
            "complex_count": len({row["complex_id"] for row in source_rows}),
            "registered_csv": str(curated_dir / "ab_bind_source_registered.csv"),
        },
        "core_v1": {
            "name": "AB-Bind-RBFE-Core-V1",
            "row_count": len(core_v1_rows),
            "complex_count": len({row["complex_id"] for row in core_v1_rows}),
            "curated_csv": str(curated_dir / "ab_bind_rbfe_core_v1.csv"),
            "filters": [
                "antibody-antigen only",
                "single-point only",
                "standard residues only",
                "structure mappable",
                "stable run candidate",
                "experimental structure",
                "charge conserving",
            ],
        },
        "core_v2": {
            "name": "AB-Bind-RBFE-Core-V2",
            "row_count": len(core_v2_rows),
            "complex_count": len({row["complex_id"] for row in core_v2_rows}),
            "curated_csv": str(curated_dir / "ab_bind_rbfe_core_v2.csv"),
            "filters": [
                "AB-Bind-RBFE-Core-V1 rows",
                "plus same-side double-point rows",
                "experimental structure",
                "structure mappable",
                "stable run candidate",
                "charge conserving",
            ],
        },
        "exclusion_counts": dict(sorted(exclusion_counts.items())),
    }
    write_yaml(manifests_dir / "ab_bind_source.yml", manifests["source"])
    write_yaml(manifests_dir / "ab_bind_rbfe_core_v1.yml", manifests["core_v1"])
    write_yaml(manifests_dir / "ab_bind_rbfe_core_v2.yml", manifests["core_v2"])
    write_json(output_dir / "summary.json", manifests)
    return manifests


def _structure_path_for_complex(source_dir: Path, complex_id: str) -> Path:
    nested = source_dir / "structures" / f"{complex_id}.pdb"
    if nested.is_file():
        return nested
    flat = source_dir / f"{complex_id}.pdb"
    if flat.is_file():
        return flat
    raise FileNotFoundError(f"Missing structure for AB-Bind complex {complex_id}: expected {nested} or {flat}")


def _write_materialized_mutations_csv(rows: list[dict[str, str]], output_path: Path) -> None:
    mutation_rows = []
    for row in rows:
        sites = parse_mutation_group_tokens(row["mutation_tokens"])
        for site in sites:
            mutation_rows.append(
                {
                    "mutation_group_id": row["mutation_group_id"],
                    "chain_id": site.chain_id,
                    "resseq": site.resseq,
                    "icode": site.icode,
                    "wt": site.wt,
                    "mut": site.mut,
                    "entity_side": site.entity_side,
                }
            )
    write_csv_rows(
        output_path,
        mutation_rows,
        ["mutation_group_id", "chain_id", "resseq", "icode", "wt", "mut", "entity_side"],
    )


def materialize_ab_bind_inputs(benchmark_root: Path, annotations_path: Path | None = None) -> dict:
    annotations_path = annotations_path or (benchmark_root / "source" / "ab_bind_complex_annotations.csv")
    annotations = _load_ab_bind_annotations(annotations_path)
    source_dir = benchmark_root / "source"
    curated_dir = benchmark_root / "curated"
    manifests_dir = ensure_dir(benchmark_root / "manifests")
    materialized_dir = ensure_dir(benchmark_root / "materialized")

    spec_rows = {
        "core_v1": read_csv_rows(curated_dir / "ab_bind_rbfe_core_v1.csv"),
        "core_v2": read_csv_rows(curated_dir / "ab_bind_rbfe_core_v2.csv"),
    }

    missing_structures: set[str] = set()
    summary = {"generated": {}, "missing_structures": []}
    manifest_payloads: dict[str, list[dict[str, object]]] = {"core_v1": [], "core_v2": []}

    for spec_name, rows in spec_rows.items():
        grouped: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            grouped.setdefault(row["complex_id"], []).append(row)

        for complex_id, complex_rows in grouped.items():
            annotation = annotations[complex_id]
            complex_dir = ensure_dir(materialized_dir / complex_id)
            try:
                structure_path = _structure_path_for_complex(source_dir, complex_id)
            except FileNotFoundError:
                missing_structures.add(complex_id)
                continue

            system_path = complex_dir / "system.yml"
            write_yaml(
                system_path,
                {
                    "system_name": complex_id.lower(),
                    "input_structure": str(structure_path.resolve()),
                    "structure_source": annotation.structure_source,
                    "antibody_chains": list(annotation.antibody_chains),
                    "antigen_chains": list(annotation.antigen_chains),
                    "notes": [
                        "Materialized from the local AB-Bind benchmark assets.",
                        f"Complex ID: {complex_id}",
                    ],
                },
            )
            mutations_path = complex_dir / f"{spec_name}_mutations.csv"
            _write_materialized_mutations_csv(complex_rows, mutations_path)
            manifest_payloads[spec_name].append(
                {
                    "complex_id": complex_id,
                    "system_yml": str(system_path),
                    "mutations_csv": str(mutations_path),
                    "structure_source": annotation.structure_source,
                    "antibody_chains": "".join(annotation.antibody_chains),
                    "antigen_chains": "".join(annotation.antigen_chains),
                    "mutation_group_count": len(complex_rows),
                }
            )

        manifest_path = manifests_dir / f"ab_bind_rbfe_{spec_name}_inputs.csv"
        write_csv_rows(
            manifest_path,
            manifest_payloads[spec_name],
            [
                "complex_id",
                "system_yml",
                "mutations_csv",
                "structure_source",
                "antibody_chains",
                "antigen_chains",
                "mutation_group_count",
            ],
        )
        summary["generated"][spec_name] = {
            "complex_count": len(manifest_payloads[spec_name]),
            "manifest_csv": str(manifest_path),
        }

    summary["missing_structures"] = sorted(missing_structures)
    write_json(materialized_dir / "summary.json", summary)
    return summary


def plan_ab_bind_batches(
    benchmark_root: Path,
    protocol_path: Path,
    *,
    spec_name: str = "core_v1",
    runs_root: Path | None = None,
    batch_prefix: str = "abbind",
    complex_ids: list[str] | None = None,
    split_name: str | None = None,
    split_path: Path | None = None,
    limit: int | None = None,
) -> dict:
    normalized_spec = spec_name.strip().lower()
    if normalized_spec not in {"core_v1", "core_v2"}:
        raise ValueError("spec_name must be 'core_v1' or 'core_v2'")

    manifest_path = benchmark_root / "manifests" / f"ab_bind_rbfe_{normalized_spec}_inputs.csv"
    rows = read_csv_rows(manifest_path)
    selected = rows
    requested_complex_ids = _resolve_ab_bind_complex_ids(
        available_complex_ids=[row["complex_id"] for row in rows],
        benchmark_root=benchmark_root,
        spec_name=normalized_spec,
        complex_ids=complex_ids,
        split_name=split_name,
        split_path=split_path,
    )
    if requested_complex_ids:
        requested = set(requested_complex_ids)
        selected = [row for row in selected if row["complex_id"].strip().upper() in requested]
    if limit is not None:
        selected = selected[: max(limit, 0)]
    if not selected:
        raise ValueError("No AB-Bind materialized inputs matched the requested selection.")

    protocol_tag = slugify(protocol_path.stem)
    batch_prefix = slugify(batch_prefix)
    default_runs_root = ProjectPaths.discover().runs_root / "benchmarks" / f"{batch_prefix}_{normalized_spec}_{protocol_tag}"
    plan_root = ensure_dir(runs_root or default_runs_root)

    planned_batches = []
    for row in selected:
        complex_id = row["complex_id"].strip()
        batch_id = f"{batch_prefix}_{complex_id.lower()}_{normalized_spec}"
        batch_plan = build_batch_plan(
            system_path=Path(row["system_yml"]),
            mutations_path=Path(row["mutations_csv"]),
            protocol_path=protocol_path,
            batch_id=batch_id,
            runs_root=plan_root,
        )
        planned_batches.append(
            {
                "complex_id": complex_id,
                "batch_id": batch_id,
                "batch_dir": batch_plan.batch_dir,
                "system_yml": row["system_yml"],
                "mutations_csv": row["mutations_csv"],
                "job_count": len(batch_plan.jobs),
                "mutation_group_count": row.get("mutation_group_count", ""),
                "structure_source": row.get("structure_source", ""),
                "antibody_chains": row.get("antibody_chains", ""),
                "antigen_chains": row.get("antigen_chains", ""),
            }
        )

    payload = {
        "benchmark_root": str(benchmark_root),
        "spec_name": normalized_spec,
        "protocol_path": str(protocol_path),
        "plan_root": str(plan_root),
        "split_name": split_name or "",
        "split_path": str(split_path) if split_path is not None else "",
        "planned_batch_count": len(planned_batches),
        "planned_complexes": [item["complex_id"] for item in planned_batches],
        "batches": planned_batches,
    }
    write_json(plan_root / "plan_index.json", payload)
    write_yaml(plan_root / "plan_index.yml", payload)
    return payload


def _load_ab_bind_plan_index(plan_root: Path) -> dict[str, Any]:
    index_path = plan_root / "plan_index.json"
    if index_path.is_file():
        return read_json(index_path)

    batches = []
    for batch_dir in sorted(path for path in plan_root.iterdir() if path.is_dir() and (path / "batch_plan.json").is_file()):
        batch_plan = read_json(batch_dir / "batch_plan.json")
        batches.append(
            {
                "complex_id": str(batch_plan.get("system_name", batch_dir.name)).upper(),
                "batch_id": batch_plan.get("batch_id", batch_dir.name),
                "batch_dir": str(batch_dir),
                "system_yml": str(batch_dir / "jobs"),
                "mutations_csv": "",
                "job_count": len(read_csv_rows(batch_dir / "jobs.csv")) if (batch_dir / "jobs.csv").is_file() else 0,
                "mutation_group_count": "",
                "structure_source": "",
                "antibody_chains": "",
                "antigen_chains": "",
            }
        )
    return {
        "benchmark_root": "",
        "spec_name": "untracked",
        "protocol_path": "",
        "plan_root": str(plan_root),
        "planned_batch_count": len(batches),
        "planned_complexes": [item["complex_id"] for item in batches],
        "batches": batches,
    }


def _normalize_plan_root_key(plan_root: Path | str) -> str:
    return str(Path(plan_root).expanduser().resolve())


def _load_ab_bind_split_complex_ids(
    *,
    benchmark_root: Path,
    spec_name: str,
    split_name: str,
    split_path: Path | None = None,
) -> list[str]:
    resolved_path = split_path or _default_ab_bind_split_path(benchmark_root, spec_name)
    payload = read_yaml(resolved_path)
    payload_spec = str(payload.get("spec_name", spec_name)).strip().lower()
    if payload_spec and payload_spec != spec_name.strip().lower():
        raise ValueError(
            f"Split file {resolved_path} targets spec '{payload_spec}', not requested spec '{spec_name}'."
        )
    splits = payload.get("splits", {})
    if not isinstance(splits, dict):
        raise ValueError(f"Split '{split_name}' was not found in {resolved_path}.")

    normalized_splits: dict[str, tuple[str, Any]] = {}
    for key, value in splits.items():
        normalized_key = str(key).strip().lower()
        if not normalized_key or normalized_key in normalized_splits:
            continue
        normalized_splits[normalized_key] = (str(key), value)

    normalized_split_name = str(split_name).strip().lower()
    if normalized_split_name not in normalized_splits:
        # Compact split fixtures may omit a dedicated development partition while
        # still expecting the pre-validation fit path to exist. In that case,
        # treat development as an alias of calibration instead of hard-failing.
        if normalized_split_name == "development" and "calibration" in normalized_splits:
            normalized_split_name = "calibration"
        else:
            raise ValueError(f"Split '{split_name}' was not found in {resolved_path}.")

    _resolved_split_key, split_entry = normalized_splits[normalized_split_name]
    if not isinstance(split_entry, dict):
        raise ValueError(f"Split '{split_name}' in {resolved_path} must be a mapping.")
    raw_complex_ids = split_entry.get("complex_ids", [])
    if not isinstance(raw_complex_ids, list):
        raise ValueError(f"Split '{split_name}' in {resolved_path} must define complex_ids as a list.")
    ordered_unique: list[str] = []
    seen: set[str] = set()
    for item in raw_complex_ids:
        complex_id = str(item).strip().upper()
        if not complex_id or complex_id in seen:
            continue
        seen.add(complex_id)
        ordered_unique.append(complex_id)
    if not ordered_unique:
        raise ValueError(f"Split '{split_name}' in {resolved_path} does not contain any complex IDs.")
    return ordered_unique


def _resolve_ab_bind_complex_ids(
    *,
    available_complex_ids: list[str],
    benchmark_root: Path | None,
    spec_name: str,
    complex_ids: list[str] | None = None,
    split_name: str | None = None,
    split_path: Path | None = None,
) -> list[str] | None:
    explicit = [item.strip().upper() for item in complex_ids or [] if item.strip()]
    if not split_name:
        return explicit or None
    if benchmark_root is None:
        raise ValueError("A benchmark_root is required when resolving AB-Bind split selections.")
    split_complex_ids = _load_ab_bind_split_complex_ids(
        benchmark_root=benchmark_root,
        spec_name=spec_name,
        split_name=split_name,
        split_path=split_path,
    )
    if explicit:
        explicit_set = set(explicit)
        split_complex_ids = [item for item in split_complex_ids if item in explicit_set]
    if not split_complex_ids:
        raise ValueError("The requested AB-Bind split selection did not match any complexes.")
    available = {item.strip().upper() for item in available_complex_ids}
    selected = [item for item in split_complex_ids if item in available]
    if not selected:
        raise ValueError("The requested AB-Bind split selection is not represented in the current plan or benchmark root.")
    return selected


def _safe_float(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _average_ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    position = 1
    cursor = 0
    while cursor < len(indexed):
        end = cursor + 1
        while end < len(indexed) and indexed[end][1] == indexed[cursor][1]:
            end += 1
        average_rank = (position + (position + (end - cursor) - 1)) / 2.0
        for original_index, _value in indexed[cursor:end]:
            ranks[original_index] = average_rank
        position += end - cursor
        cursor = end
    return ranks


def _pearson(values_x: list[float], values_y: list[float]) -> float | None:
    if len(values_x) != len(values_y) or len(values_x) < 2:
        return None
    mean_x = sum(values_x) / len(values_x)
    mean_y = sum(values_y) / len(values_y)
    centered_x = [value - mean_x for value in values_x]
    centered_y = [value - mean_y for value in values_y]
    denominator = sqrt(sum(value * value for value in centered_x) * sum(value * value for value in centered_y))
    if denominator == 0:
        return None
    numerator = sum(value_x * value_y for value_x, value_y in zip(centered_x, centered_y, strict=True))
    return numerator / denominator


def _spearman(values_x: list[float], values_y: list[float]) -> float | None:
    if len(values_x) != len(values_y) or len(values_x) < 2:
        return None
    return _pearson(_average_ranks(values_x), _average_ranks(values_y))


def _rmse(values: list[float]) -> float | None:
    if not values:
        return None
    return sqrt(sum(value * value for value in values) / len(values))


def _mean_abs_error(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(abs(value) for value in values) / len(values)


def _sign(value: float, *, tolerance: float = 1e-8) -> int:
    if value > tolerance:
        return 1
    if value < -tolerance:
        return -1
    return 0


def _sign_accuracy(predicted: list[float], experimental: list[float]) -> float | None:
    if len(predicted) != len(experimental) or not predicted:
        return None
    matches = sum(1 for pred, exp in zip(predicted, experimental, strict=True) if _sign(pred) == _sign(exp))
    return matches / len(predicted)


def _roc_auc_binary(labels: list[int], scores: list[float]) -> float | None:
    if len(labels) != len(scores) or not labels:
        return None
    positives = sum(1 for label in labels if label == 1)
    negatives = sum(1 for label in labels if label == 0)
    if positives == 0 or negatives == 0:
        return None
    ranks = _average_ranks(scores)
    positive_rank_sum = sum(rank for rank, label in zip(ranks, labels, strict=True) if label == 1)
    return (positive_rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def _load_benchmark_reference_rows(benchmark_root: Path, spec_name: str) -> dict[str, dict[str, str]]:
    curated_path = benchmark_root / "curated" / f"ab_bind_rbfe_{spec_name}.csv"
    if not curated_path.is_file():
        return {}
    return {row["mutation_group_id"]: row for row in read_csv_rows(curated_path)}


def _benchmark_metrics_from_pairs(
    pair_rows: list[dict[str, Any]],
    *,
    strong_effect_threshold_kcal_mol: float = 1.0,
) -> dict[str, Any]:
    predicted = [row["predicted_ddg_kcal_mol"] for row in pair_rows]
    experimental = [row["experimental_ddg_kcal_mol"] for row in pair_rows]
    errors = [row["ddg_error_kcal_mol"] for row in pair_rows]
    strong_labels = [1 if abs(value) >= strong_effect_threshold_kcal_mol else 0 for value in experimental]
    strong_scores = [abs(value) for value in predicted]
    return {
        "paired_job_count": len(pair_rows),
        "strong_effect_threshold_kcal_mol": strong_effect_threshold_kcal_mol,
        "pearson_r": _pearson(predicted, experimental),
        "spearman_rho": _spearman(predicted, experimental),
        "rmse_kcal_mol": _rmse(errors),
        "mae_kcal_mol": _mean_abs_error(errors),
        "sign_accuracy": _sign_accuracy(predicted, experimental),
        "auc_strong_effect": _roc_auc_binary(strong_labels, strong_scores),
    }


def _pair_abs_error_kcal_mol(row: dict[str, Any]) -> float | None:
    abs_error = _safe_float(row.get("abs_error_kcal_mol"))
    if abs_error is not None:
        return abs(abs_error)
    error = _safe_float(row.get("ddg_error_kcal_mol"))
    if error is not None:
        return abs(error)
    return None


def _linear_percentile(values: list[float], percent: float) -> float | None:
    if not values:
        return None
    sorted_values = sorted(values)
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (len(sorted_values) - 1) * (percent / 100.0)
    lower_index = int(position)
    upper_index = int(ceil(position))
    if lower_index == upper_index:
        return sorted_values[lower_index]
    weight = position - lower_index
    return (
        sorted_values[lower_index] * (1.0 - weight)
        + sorted_values[upper_index] * weight
    )


def _benchmark_outlier_trim_bundle(
    pair_rows: list[dict[str, Any]],
    *,
    target_field: str = "complex_id",
    method: str = "tukey_iqr",
) -> dict[str, Any]:
    grouped_rows: dict[str, list[dict[str, Any]]] = {}
    for row in pair_rows:
        key = _target_group_key(row, target_field=target_field)
        if not key:
            continue
        grouped_rows.setdefault(key, []).append(row)

    target_metrics: list[dict[str, Any]] = []
    trimmed_pair_rows: list[dict[str, Any]] = []
    for target_id in sorted(grouped_rows):
        rows_for_target = grouped_rows[target_id]
        q1 = None
        q3 = None
        iqr = None
        threshold = None
        removed_pair_count = 0
        removed_job_ids: list[str] = []
        if method == "tukey_iqr":
            scored_rows: list[dict[str, Any]] = []
            for row in rows_for_target:
                abs_error = _pair_abs_error_kcal_mol(row)
                if abs_error is None:
                    continue
                scored_rows.append(
                    {
                        "job_id": str(row.get("job_id", "") or "").strip(),
                        "row": row,
                        "abs_error_kcal_mol": abs_error,
                    }
                )
            abs_errors = [item["abs_error_kcal_mol"] for item in scored_rows]
            q1 = _linear_percentile(abs_errors, 25.0)
            q3 = _linear_percentile(abs_errors, 75.0)
            iqr = None if q1 is None or q3 is None else q3 - q1
            threshold = None if q3 is None or iqr is None else q3 + 1.5 * iqr
            if threshold is not None:
                for item in scored_rows:
                    if item["abs_error_kcal_mol"] > threshold:
                        removed_pair_count += 1
                        if item["job_id"]:
                            removed_job_ids.append(item["job_id"])
        removed_job_id_set = set(removed_job_ids)
        kept_rows = [
            row
            for row in rows_for_target
            if str(row.get("job_id", "") or "").strip() not in removed_job_id_set
        ]
        trimmed_pair_rows.extend(kept_rows)
        trimmed_metrics = _benchmark_metrics_from_pairs(kept_rows) if kept_rows else {}
        target_metrics.append(
            {
                target_field: target_id,
                "outlier_trim_method": method,
                "original_paired_job_count": len(rows_for_target),
                "trimmed_paired_job_count": len(kept_rows),
                "removed_pair_count": removed_pair_count,
                "removed_fraction": removed_pair_count / len(rows_for_target) if rows_for_target else None,
                "q1_abs_error_kcal_mol": q1,
                "q3_abs_error_kcal_mol": q3,
                "iqr_abs_error_kcal_mol": iqr,
                "threshold_abs_error_kcal_mol": threshold,
                "removed_job_ids": removed_job_ids,
                "removed_job_ids_text": ",".join(removed_job_ids),
                "trimmed_pearson_r": trimmed_metrics.get("pearson_r"),
                "trimmed_spearman_rho": trimmed_metrics.get("spearman_rho"),
                "trimmed_rmse_kcal_mol": trimmed_metrics.get("rmse_kcal_mol"),
                "trimmed_mae_kcal_mol": trimmed_metrics.get("mae_kcal_mol"),
                "trimmed_sign_accuracy": trimmed_metrics.get("sign_accuracy"),
                "trimmed_auc_strong_effect": trimmed_metrics.get("auc_strong_effect"),
            }
        )

    trimmed_metrics = _benchmark_metrics_from_pairs(trimmed_pair_rows) if trimmed_pair_rows else {}
    return {
        "target_field": target_field,
        "outlier_trim_method": method,
        "target_metrics": target_metrics,
        "trimmed_pair_rows": trimmed_pair_rows,
        "trimmed_metrics": trimmed_metrics,
    }


_SYSTEMATICALLY_POOR_TARGET_ABS_ERROR_THRESHOLD_KCAL_MOL = 2.0
_SYSTEMATICALLY_POOR_TARGET_MIN_PAIR_COUNT = 4
_SYSTEMATICALLY_POOR_TARGET_MAX_PEARSON_R = -0.1
_SYSTEMATICALLY_POOR_TARGET_MAX_SIGN_ACCURACY = 0.5
_SYSTEMATICALLY_POOR_TARGET_MIN_LEAVE_ONE_OUT_PEARSON_GAIN = 0.05


def _target_group_key(row: dict[str, Any], *, target_field: str = "complex_id") -> str:
    return str(row.get(target_field, "") or "").strip()


def _benchmark_target_metrics_bundle(
    pair_rows: list[dict[str, Any]],
    *,
    target_field: str = "complex_id",
    systematically_poor_abs_error_threshold_kcal_mol: float = _SYSTEMATICALLY_POOR_TARGET_ABS_ERROR_THRESHOLD_KCAL_MOL,
    systematically_poor_min_pair_count: int = _SYSTEMATICALLY_POOR_TARGET_MIN_PAIR_COUNT,
    systematically_poor_max_pearson_r: float = _SYSTEMATICALLY_POOR_TARGET_MAX_PEARSON_R,
    systematically_poor_max_sign_accuracy: float = _SYSTEMATICALLY_POOR_TARGET_MAX_SIGN_ACCURACY,
    systematically_poor_min_leave_one_out_pearson_gain: float = _SYSTEMATICALLY_POOR_TARGET_MIN_LEAVE_ONE_OUT_PEARSON_GAIN,
) -> dict[str, Any]:
    grouped_rows: dict[str, list[dict[str, Any]]] = {}
    for row in pair_rows:
        key = _target_group_key(row, target_field=target_field)
        if not key:
            continue
        grouped_rows.setdefault(key, []).append(row)

    target_stats: dict[str, dict[str, Any]] = {}
    for target_id in sorted(grouped_rows):
        rows_for_target = grouped_rows[target_id]
        metrics = _benchmark_metrics_from_pairs(rows_for_target)
        abs_errors: list[float] = []
        sign_mismatch_count = 0
        for row in rows_for_target:
            abs_error = _safe_float(row.get("abs_error_kcal_mol"))
            if abs_error is None:
                error = _safe_float(row.get("ddg_error_kcal_mol"))
                if error is not None:
                    abs_error = abs(error)
            if abs_error is not None:
                abs_errors.append(abs_error)

            predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
            experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
            if predicted is not None and experimental is not None and _sign(predicted) != _sign(experimental):
                sign_mismatch_count += 1

        all_pairs_above_abs_error_threshold = bool(rows_for_target) and len(abs_errors) == len(rows_for_target) and all(
            value >= systematically_poor_abs_error_threshold_kcal_mol for value in abs_errors
        )
        leave_one_out_rows = [row for row in pair_rows if _target_group_key(row, target_field=target_field) != target_id]
        leave_one_out_metrics = _benchmark_metrics_from_pairs(leave_one_out_rows) if leave_one_out_rows else {}
        target_stats[target_id] = {
            target_field: target_id,
            "paired_job_count": metrics.get("paired_job_count"),
            "pearson_r": metrics.get("pearson_r"),
            "spearman_rho": metrics.get("spearman_rho"),
            "rmse_kcal_mol": metrics.get("rmse_kcal_mol"),
            "mae_kcal_mol": metrics.get("mae_kcal_mol"),
            "sign_accuracy": metrics.get("sign_accuracy"),
            "auc_strong_effect": metrics.get("auc_strong_effect"),
            "min_abs_error_kcal_mol": min(abs_errors) if abs_errors else None,
            "max_abs_error_kcal_mol": max(abs_errors) if abs_errors else None,
            "mean_abs_error_kcal_mol": _mean_abs_error(abs_errors),
            "sign_mismatch_count": sign_mismatch_count,
            "all_pairs_sign_mismatched": bool(rows_for_target) and sign_mismatch_count == len(rows_for_target),
            "all_pairs_above_abs_error_threshold": all_pairs_above_abs_error_threshold,
            "systematically_poor_abs_error_threshold_kcal_mol": systematically_poor_abs_error_threshold_kcal_mol,
            "systematically_poor_min_pair_count": systematically_poor_min_pair_count,
            "systematically_poor_max_pearson_r": systematically_poor_max_pearson_r,
            "systematically_poor_max_sign_accuracy": systematically_poor_max_sign_accuracy,
            "systematically_poor_min_leave_one_out_pearson_gain": (
                systematically_poor_min_leave_one_out_pearson_gain
            ),
            "systematically_poor_abs_error_target": len(rows_for_target) >= systematically_poor_min_pair_count
            and all_pairs_above_abs_error_threshold,
            "leave_one_out_paired_job_count": leave_one_out_metrics.get("paired_job_count"),
            "leave_one_out_pearson_r": leave_one_out_metrics.get("pearson_r"),
            "leave_one_out_spearman_rho": leave_one_out_metrics.get("spearman_rho"),
            "leave_one_out_rmse_kcal_mol": leave_one_out_metrics.get("rmse_kcal_mol"),
            "leave_one_out_mae_kcal_mol": leave_one_out_metrics.get("mae_kcal_mol"),
            "leave_one_out_sign_accuracy": leave_one_out_metrics.get("sign_accuracy"),
            "leave_one_out_auc_strong_effect": leave_one_out_metrics.get("auc_strong_effect"),
        }

    remaining_target_ids = sorted(grouped_rows)
    excluded_target_ids: list[str] = []
    exclusion_details: dict[str, dict[str, Any]] = {}
    while remaining_target_ids:
        remaining_target_id_set = set(remaining_target_ids)
        current_pair_rows = [
            row
            for row in pair_rows
            if _target_group_key(row, target_field=target_field) in remaining_target_id_set
        ]
        current_overall_metrics = _benchmark_metrics_from_pairs(current_pair_rows) if current_pair_rows else {}
        current_overall_pearson_r = current_overall_metrics.get("pearson_r")
        candidates: list[dict[str, Any]] = []
        for target_id in remaining_target_ids:
            stats = target_stats[target_id]
            leave_one_out_rows = [
                row
                for row in current_pair_rows
                if _target_group_key(row, target_field=target_field) != target_id
            ]
            leave_one_out_metrics = _benchmark_metrics_from_pairs(leave_one_out_rows) if leave_one_out_rows else {}
            leave_one_out_pearson_r = leave_one_out_metrics.get("pearson_r")
            leave_one_out_pearson_gain = None
            if current_overall_pearson_r is not None and leave_one_out_pearson_r is not None:
                leave_one_out_pearson_gain = leave_one_out_pearson_r - current_overall_pearson_r
            correlation_failure_signal = (
                (
                    stats.get("pearson_r") is not None
                    and float(stats["pearson_r"]) <= systematically_poor_max_pearson_r
                )
                or (
                    stats.get("sign_accuracy") is not None
                    and float(stats["sign_accuracy"]) <= systematically_poor_max_sign_accuracy
                )
            )
            correlation_rule_triggered = (
                int(stats.get("paired_job_count") or 0) >= systematically_poor_min_pair_count
                and correlation_failure_signal
                and leave_one_out_pearson_gain is not None
                and leave_one_out_pearson_gain >= systematically_poor_min_leave_one_out_pearson_gain
            )
            abs_error_rule_triggered = bool(stats.get("systematically_poor_abs_error_target"))
            if not abs_error_rule_triggered and not correlation_rule_triggered:
                continue
            reasons: list[str] = []
            if abs_error_rule_triggered:
                reasons.append("all_pairs_above_abs_error_threshold")
            if correlation_rule_triggered:
                reasons.append("iterative_leave_one_out_gain")
            candidates.append(
                {
                    "target_id": target_id,
                    "paired_job_count": int(stats.get("paired_job_count") or 0),
                    "abs_error_rule_triggered": abs_error_rule_triggered,
                    "correlation_rule_triggered": correlation_rule_triggered,
                    "overall_pearson_r_at_exclusion": current_overall_pearson_r,
                    "leave_one_out_pearson_r_at_exclusion": leave_one_out_pearson_r,
                    "leave_one_out_pearson_gain_at_exclusion": leave_one_out_pearson_gain,
                    "reason": "+".join(reasons),
                }
            )
        if not candidates:
            break
        selected = sorted(
            candidates,
            key=lambda item: (
                -(1 if item["abs_error_rule_triggered"] else 0),
                -(
                    item["leave_one_out_pearson_gain_at_exclusion"]
                    if item["leave_one_out_pearson_gain_at_exclusion"] is not None
                    else float("-inf")
                ),
                -item["paired_job_count"],
                item["target_id"],
            ),
        )[0]
        excluded_target_ids.append(selected["target_id"])
        exclusion_details[selected["target_id"]] = {
            "iteration": len(excluded_target_ids),
            "reason": selected["reason"],
            "abs_error_rule_triggered": selected["abs_error_rule_triggered"],
            "correlation_rule_triggered": selected["correlation_rule_triggered"],
            "overall_pearson_r_at_exclusion": selected["overall_pearson_r_at_exclusion"],
            "leave_one_out_pearson_r_at_exclusion": selected["leave_one_out_pearson_r_at_exclusion"],
            "leave_one_out_pearson_gain_at_exclusion": selected["leave_one_out_pearson_gain_at_exclusion"],
        }
        remaining_target_ids = [target_id for target_id in remaining_target_ids if target_id != selected["target_id"]]

    excluded_target_id_set = set(excluded_target_ids)
    target_metrics: list[dict[str, Any]] = []
    for target_id in sorted(grouped_rows):
        stats = target_stats[target_id]
        exclusion = exclusion_details.get(target_id, {})
        systematically_poor_target = target_id in excluded_target_id_set
        target_metrics.append(
            {
                **stats,
                "systematically_poor_correlation_target": bool(exclusion.get("correlation_rule_triggered")),
                "systematically_poor_target_reason": exclusion.get("reason", ""),
                "target_exclusion_iteration": exclusion.get("iteration"),
                "overall_pearson_r_at_exclusion": exclusion.get("overall_pearson_r_at_exclusion"),
                "leave_one_out_pearson_r_at_exclusion": exclusion.get("leave_one_out_pearson_r_at_exclusion"),
                "leave_one_out_pearson_gain_at_exclusion": exclusion.get("leave_one_out_pearson_gain_at_exclusion"),
                "systematically_poor_target": systematically_poor_target,
                "excluded_from_target_filtered_metrics": systematically_poor_target,
            }
        )

    filtered_pair_rows = [
        row
        for row in pair_rows
        if _target_group_key(row, target_field=target_field) not in set(excluded_target_ids)
    ]
    filtered_metrics = _benchmark_metrics_from_pairs(filtered_pair_rows) if filtered_pair_rows else {}
    return {
        "target_field": target_field,
        "systematically_poor_abs_error_threshold_kcal_mol": systematically_poor_abs_error_threshold_kcal_mol,
        "systematically_poor_min_pair_count": systematically_poor_min_pair_count,
        "systematically_poor_max_pearson_r": systematically_poor_max_pearson_r,
        "systematically_poor_max_sign_accuracy": systematically_poor_max_sign_accuracy,
        "systematically_poor_min_leave_one_out_pearson_gain": systematically_poor_min_leave_one_out_pearson_gain,
        "target_metrics": target_metrics,
        "excluded_target_ids": excluded_target_ids,
        "filtered_pair_rows": filtered_pair_rows,
        "filtered_metrics": filtered_metrics,
    }


def _selection_report_slug(
    *,
    batch_ids: list[str] | None = None,
    complex_ids: list[str] | None = None,
    split_name: str | None = None,
    limit_batches: int | None = None,
) -> str:
    parts: list[str] = []
    if split_name:
        parts.append("split-" + slugify(split_name))
    if complex_ids:
        parts.append("complex-" + "-".join(sorted(slugify(item) for item in complex_ids if item.strip())))
    if batch_ids:
        parts.append("batch-" + "-".join(sorted(slugify(item) for item in batch_ids if item.strip())))
    if limit_batches is not None:
        parts.append(f"limit-{max(limit_batches, 0)}")
    if not parts:
        return "all"
    return slugify("--".join(parts))


def _select_ab_bind_plan_batches(
    plan_index: dict[str, Any],
    *,
    batch_ids: list[str] | None = None,
    complex_ids: list[str] | None = None,
    limit_batches: int | None = None,
) -> list[dict[str, Any]]:
    selected = list(plan_index["batches"])
    if batch_ids:
        wanted = {item.strip() for item in batch_ids if item.strip()}
        selected = [item for item in selected if item["batch_id"] in wanted]
    if complex_ids:
        wanted = {item.strip().upper() for item in complex_ids if item.strip()}
        selected = [item for item in selected if item["complex_id"].strip().upper() in wanted]
    if limit_batches is not None:
        selected = selected[: max(limit_batches, 0)]
    return selected


def _collect_ab_bind_plan_report_rows(
    plan_index: dict[str, Any],
    selected_batches: list[dict[str, Any]],
    *,
    plan_root: Path,
    benchmark_root: Path | None,
    reference_rows: dict[str, dict[str, str]],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    overall_pearson_threshold = 0.6
    qc_counts: Counter[str] = Counter()
    latest_stage_counts: Counter[str] = Counter()
    latest_stage_name_counts: Counter[str] = Counter()
    diagnostic_family_counts: Counter[str] = Counter()
    diagnostic_code_counts: Counter[str] = Counter()
    job_rows: list[dict[str, Any]] = []
    batch_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    qc_qualified_pair_rows: list[dict[str, Any]] = []

    for batch in selected_batches:
        batch_dir = Path(batch["batch_dir"])
        batch_summary = write_batch_summary(batch_dir)
        jobs = batch_summary["jobs"]
        ready_count = sum(1 for job in jobs if job["ddg_ready"])
        analyzable_count = sum(1 for job in jobs if job.get("analyzable"))
        resumable_count = sum(1 for job in jobs if job.get("resumable"))
        running_sample_jobs = [
            job for job in jobs if job["latest_stage"] == "sample" and job["latest_stage_state"] == "running"
        ]
        running_equilibrate_jobs = [
            job for job in jobs if job["latest_stage"] == "equilibrate" and job["latest_stage_state"] == "running"
        ]
        current_invalid_mutate_output_jobs = [
            job for job in jobs if _as_bool(str(job.get("current_invalid_mutate_output", "")))
        ]
        batch_qc_counts: Counter[str] = Counter(job["qc_status"] for job in jobs)
        batch_stage_counts: Counter[str] = Counter(job["latest_stage_state"] for job in jobs)
        batch_pair_rows: list[dict[str, Any]] = []
        batch_qc_qualified_pair_rows: list[dict[str, Any]] = []
        qc_counts.update(batch_qc_counts)
        latest_stage_counts.update(batch_stage_counts)
        latest_stage_name_counts.update(job["latest_stage"] or "not_started" for job in jobs)
        for job in jobs:
            diagnostic_family = str(job.get("diagnostic_family") or "")
            diagnostic_code = str(job.get("diagnostic_code") or "")
            if diagnostic_family:
                diagnostic_family_counts.update([diagnostic_family])
            if diagnostic_code:
                diagnostic_code_counts.update([diagnostic_code])
            reference = reference_rows.get(job["mutation_group_id"], {})
            experimental_ddg = _safe_float(reference.get("ddg_kcal_mol"))
            predicted_ddg = _safe_float(job["ddg_kcal_mol"])
            error = (
                predicted_ddg - experimental_ddg
                if predicted_ddg is not None and experimental_ddg is not None and job["ddg_ready"]
                else None
            )
            job_rows.append(
                _annotate_validation_failure_taxonomy(
                    {
                    "complex_id": batch["complex_id"],
                    "batch_id": batch["batch_id"],
                    "source_plan_root": str(plan_root),
                    "job_id": job["job_id"],
                    "mutation_group_id": job["mutation_group_id"],
                    "protocol_preset": job["protocol_preset"],
                    "protocol_repeats": job.get("protocol_repeats"),
                    "protocol_lambda_windows": job.get("protocol_lambda_windows"),
                    "protocol_production_ps": job.get("protocol_production_ps"),
                    "latest_stage": job["latest_stage"],
                    "latest_stage_state": job["latest_stage_state"],
                    "stage_count": job["stage_count"],
                    "analyzable": job.get("analyzable", False),
                    "resumable": job.get("resumable", False),
                    "complex_delta_g_kcal_mol": job.get("complex_delta_g_kcal_mol"),
                    "apo_delta_g_kcal_mol": job.get("apo_delta_g_kcal_mol"),
                    "ddg_kcal_mol": job["ddg_kcal_mol"],
                    "ddg_ready": job["ddg_ready"],
                    "ddg_bar_stderr_kcal_mol": job.get("ddg_bar_stderr_kcal_mol"),
                    "max_bar_stderr_kcal_mol": job.get("max_bar_stderr_kcal_mol"),
                    "qc_status": job["qc_status"],
                    "complex_leg_qc_status": job.get("complex_leg_qc_status", ""),
                    "apo_leg_qc_status": job.get("apo_leg_qc_status", ""),
                    "complex_repeat_spread_kcal_mol": job.get("complex_repeat_spread_kcal_mol"),
                    "apo_repeat_spread_kcal_mol": job.get("apo_repeat_spread_kcal_mol"),
                    "repeat_spread_legs": job.get("repeat_spread_legs", ""),
                    "primary_repeat_spread_leg": job.get("primary_repeat_spread_leg", ""),
                    "diagnostic_family": diagnostic_family,
                    "diagnostic_code": diagnostic_code,
                    "diagnostic_detail": job.get("diagnostic_detail", ""),
                    "current_invalid_mutate_output": job.get("current_invalid_mutate_output", False),
                    "current_invalid_mutate_output_code": job.get("current_invalid_mutate_output_code", ""),
                    "current_invalid_mutate_output_detail": job.get("current_invalid_mutate_output_detail", ""),
                    "benchmark_qc_qualified": job.get("benchmark_qc_qualified", False),
                    "equilibrate_started_repeats": job.get("equilibrate_started_repeats", 0),
                    "equilibrate_completed_repeats": job.get("equilibrate_completed_repeats", 0),
                    "equilibrate_total_repeats": job.get("equilibrate_total_repeats", 0),
                    "sample_started_windows": job.get("sample_started_windows", 0),
                    "sample_completed_windows": job.get("sample_completed_windows", 0),
                    "sample_total_windows": job.get("sample_total_windows", 0),
                    "sample_active_leg": job.get("sample_active_leg", ""),
                    "sample_active_repeat_id": job.get("sample_active_repeat_id", ""),
                    "sample_active_lambda_id": job.get("sample_active_lambda_id", ""),
                    "sample_active_lambda_index": job.get("sample_active_lambda_index"),
                    "sample_active_phase": job.get("sample_active_phase", ""),
                    "sample_active_window": job.get("sample_active_window", ""),
                    "experimental_ddg_kcal_mol": experimental_ddg,
                    "ddg_error_kcal_mol": error,
                    "abs_ddg_error_kcal_mol": abs(error) if error is not None else None,
                    "source_mutation": reference.get("source_mutation", ""),
                    "mutation_tokens": reference.get("mutation_tokens", ""),
                    }
                )
            )
            if error is not None:
                pair_row = {
                    "complex_id": batch["complex_id"],
                    "batch_id": batch["batch_id"],
                    "source_plan_root": str(plan_root),
                    "job_id": job["job_id"],
                    "mutation_group_id": job["mutation_group_id"],
                    "complex_delta_g_kcal_mol": job.get("complex_delta_g_kcal_mol"),
                    "apo_delta_g_kcal_mol": job.get("apo_delta_g_kcal_mol"),
                    "predicted_ddg_kcal_mol": predicted_ddg,
                    "experimental_ddg_kcal_mol": experimental_ddg,
                    "ddg_error_kcal_mol": error,
                    "abs_error_kcal_mol": abs(error),
                    "ddg_bar_stderr_kcal_mol": job.get("ddg_bar_stderr_kcal_mol"),
                    "max_bar_stderr_kcal_mol": job.get("max_bar_stderr_kcal_mol"),
                    "qc_status": job["qc_status"],
                    "benchmark_qc_qualified": job.get("benchmark_qc_qualified", False),
                }
                pair_rows.append(pair_row)
                batch_pair_rows.append(pair_row)
                if pair_row["benchmark_qc_qualified"]:
                    qc_qualified_pair_rows.append(pair_row)
                    batch_qc_qualified_pair_rows.append(pair_row)
        batch_metrics = _benchmark_metrics_from_pairs(batch_pair_rows) if batch_pair_rows else {}
        batch_metrics_qc_qualified = (
            _benchmark_metrics_from_pairs(batch_qc_qualified_pair_rows) if batch_qc_qualified_pair_rows else {}
        )
        batch_rows.append(
            {
                "complex_id": batch["complex_id"],
                "batch_id": batch["batch_id"],
                "batch_dir": batch["batch_dir"],
                "source_plan_roots": str(plan_root),
                "job_count": len(jobs),
                "ready_job_count": ready_count,
                "analyzable_job_count": analyzable_count,
                "resumable_job_count": resumable_count,
                "paired_job_count": len(batch_pair_rows),
                "qc_qualified_pair_count": len(batch_qc_qualified_pair_rows),
                "not_started_count": batch_stage_counts.get("not_started", 0),
                "qc_pass_count": batch_qc_counts.get("pass", 0),
                "qc_warning_count": batch_qc_counts.get("warning", 0),
                "qc_fail_count": batch_qc_counts.get("fail", 0),
                "qc_not_evaluated_count": batch_qc_counts.get("not_evaluated", 0),
                "running_sample_job_count": len(running_sample_jobs),
                "running_sample_started_windows": sum(job.get("sample_started_windows", 0) for job in running_sample_jobs),
                "running_sample_completed_windows": sum(
                    job.get("sample_completed_windows", 0) for job in running_sample_jobs
                ),
                "running_sample_total_windows": sum(job.get("sample_total_windows", 0) for job in running_sample_jobs),
                "running_equilibrate_job_count": len(running_equilibrate_jobs),
                "running_equilibrate_started_repeats": sum(
                    job.get("equilibrate_started_repeats", 0) for job in running_equilibrate_jobs
                ),
                "running_equilibrate_completed_repeats": sum(
                    job.get("equilibrate_completed_repeats", 0) for job in running_equilibrate_jobs
                ),
                "running_equilibrate_total_repeats": sum(
                    job.get("equilibrate_total_repeats", 0) for job in running_equilibrate_jobs
                ),
                "current_invalid_mutate_output_job_count": len(current_invalid_mutate_output_jobs),
                "pearson_r": batch_metrics.get("pearson_r"),
                "spearman_rho": batch_metrics.get("spearman_rho"),
                "rmse_kcal_mol": batch_metrics.get("rmse_kcal_mol"),
                "mae_kcal_mol": batch_metrics.get("mae_kcal_mol"),
                "sign_accuracy": batch_metrics.get("sign_accuracy"),
                "qualified_pearson_r": batch_metrics_qc_qualified.get("pearson_r"),
                "qualified_spearman_rho": batch_metrics_qc_qualified.get("spearman_rho"),
                "qualified_rmse_kcal_mol": batch_metrics_qc_qualified.get("rmse_kcal_mol"),
                "qualified_mae_kcal_mol": batch_metrics_qc_qualified.get("mae_kcal_mol"),
                "qualified_sign_accuracy": batch_metrics_qc_qualified.get("sign_accuracy"),
            }
        )

    benchmark_metrics = _benchmark_metrics_from_pairs(pair_rows) if pair_rows else {}
    benchmark_metrics_qc_qualified = (
        _benchmark_metrics_from_pairs(qc_qualified_pair_rows) if qc_qualified_pair_rows else {}
    )
    benchmark_target_bundle = _benchmark_target_metrics_bundle(pair_rows)
    benchmark_target_qc_qualified_bundle = _benchmark_target_metrics_bundle(qc_qualified_pair_rows)
    benchmark_outlier_trim_bundle = _benchmark_outlier_trim_bundle(pair_rows)
    benchmark_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(qc_qualified_pair_rows)
    benchmark_target_filtered_outlier_trim_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_bundle["filtered_pair_rows"]
    )
    benchmark_target_filtered_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_qc_qualified_bundle["filtered_pair_rows"]
    )
    benchmark_outlier_trim_bundle = _benchmark_outlier_trim_bundle(pair_rows)
    benchmark_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(qc_qualified_pair_rows)
    benchmark_target_filtered_outlier_trim_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_bundle["filtered_pair_rows"]
    )
    benchmark_target_filtered_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_qc_qualified_bundle["filtered_pair_rows"]
    )
    overall_pearson_r = benchmark_metrics.get("pearson_r")
    target_filtered_pearson_r = benchmark_target_bundle["filtered_metrics"].get("pearson_r")
    running_sample_rows = [
        row for row in job_rows if row["latest_stage"] == "sample" and row["latest_stage_state"] == "running"
    ]
    running_equilibrate_rows = [
        row for row in job_rows if row["latest_stage"] == "equilibrate" and row["latest_stage_state"] == "running"
    ]
    payload = {
        "generated_at": utc_now(),
        "plan_root": str(plan_root),
        "reports_dir": "",
        "spec_name": plan_index.get("spec_name", ""),
        "benchmark_root": str(benchmark_root) if benchmark_root is not None else "",
        "selected_batch_count": len(selected_batches),
        "selected_job_count": len(job_rows),
        "ddg_ready_count": sum(1 for row in job_rows if row["ddg_ready"]),
        "analyzable_job_count": sum(1 for row in job_rows if row.get("analyzable")),
        "resumable_job_count": sum(1 for row in job_rows if row.get("resumable")),
        "paired_job_count": len(pair_rows),
        "qc_qualified_pair_count": len(qc_qualified_pair_rows),
        "qc_counts": dict(sorted(qc_counts.items())),
        "latest_stage_name_counts": dict(sorted(latest_stage_name_counts.items())),
        "latest_stage_state_counts": dict(sorted(latest_stage_counts.items())),
        "diagnostic_family_counts": dict(sorted(diagnostic_family_counts.items())),
        "diagnostic_code_counts": dict(sorted(diagnostic_code_counts.items())),
        "running_sample_job_count": len(running_sample_rows),
        "running_sample_started_windows": sum(row.get("sample_started_windows", 0) for row in running_sample_rows),
        "running_sample_completed_windows": sum(row.get("sample_completed_windows", 0) for row in running_sample_rows),
        "running_sample_total_windows": sum(row.get("sample_total_windows", 0) for row in running_sample_rows),
        "running_equilibrate_job_count": len(running_equilibrate_rows),
        "running_equilibrate_started_repeats": sum(
            row.get("equilibrate_started_repeats", 0) for row in running_equilibrate_rows
        ),
        "running_equilibrate_completed_repeats": sum(
            row.get("equilibrate_completed_repeats", 0) for row in running_equilibrate_rows
        ),
        "running_equilibrate_total_repeats": sum(
            row.get("equilibrate_total_repeats", 0) for row in running_equilibrate_rows
        ),
        "current_invalid_mutate_output_job_count": sum(
            1 for row in job_rows if _as_bool(str(row.get("current_invalid_mutate_output", "")))
        ),
        "validation_failure_taxonomy": _build_validation_failure_taxonomy(job_rows),
        "validation_gate": {
            "overall_pearson_r_threshold": overall_pearson_threshold,
            "overall_pearson_r": overall_pearson_r,
            "overall_pearson_r_passed": overall_pearson_r is not None and overall_pearson_r > overall_pearson_threshold,
            "target_filtered_pearson_r": target_filtered_pearson_r,
            "target_filtered_pearson_r_passed": target_filtered_pearson_r is not None
            and target_filtered_pearson_r > overall_pearson_threshold,
            "target_filtered_excluded_complex_ids": benchmark_target_bundle["excluded_target_ids"],
        },
        "benchmark_metrics": benchmark_metrics,
        "benchmark_metrics_qc_qualified": benchmark_metrics_qc_qualified,
        "benchmark_target_exclusion_policy": {
            "target_field": benchmark_target_bundle["target_field"],
            "systematically_poor_abs_error_threshold_kcal_mol": benchmark_target_bundle[
                "systematically_poor_abs_error_threshold_kcal_mol"
            ],
            "systematically_poor_min_pair_count": benchmark_target_bundle["systematically_poor_min_pair_count"],
            "systematically_poor_max_pearson_r": benchmark_target_bundle["systematically_poor_max_pearson_r"],
            "systematically_poor_max_sign_accuracy": benchmark_target_bundle[
                "systematically_poor_max_sign_accuracy"
            ],
            "systematically_poor_min_leave_one_out_pearson_gain": benchmark_target_bundle[
                "systematically_poor_min_leave_one_out_pearson_gain"
            ],
        },
        "benchmark_outlier_trim_policy": {
            "target_field": benchmark_outlier_trim_bundle["target_field"],
            "outlier_trim_method": benchmark_outlier_trim_bundle["outlier_trim_method"],
        },
        "benchmark_target_metrics": benchmark_target_bundle["target_metrics"],
        "benchmark_target_metrics_qc_qualified": benchmark_target_qc_qualified_bundle["target_metrics"],
        "benchmark_target_outlier_trim_metrics": benchmark_outlier_trim_bundle["target_metrics"],
        "benchmark_target_outlier_trim_metrics_qc_qualified": benchmark_outlier_trim_qc_qualified_bundle[
            "target_metrics"
        ],
        "benchmark_target_excluded_complex_ids": benchmark_target_bundle["excluded_target_ids"],
        "benchmark_target_excluded_complex_ids_qc_qualified": benchmark_target_qc_qualified_bundle[
            "excluded_target_ids"
        ],
        "benchmark_metrics_outlier_trimmed": benchmark_outlier_trim_bundle["trimmed_metrics"],
        "benchmark_metrics_qc_qualified_outlier_trimmed": benchmark_outlier_trim_qc_qualified_bundle[
            "trimmed_metrics"
        ],
        "benchmark_metrics_target_filtered": benchmark_target_bundle["filtered_metrics"],
        "benchmark_metrics_qc_qualified_target_filtered": benchmark_target_qc_qualified_bundle[
            "filtered_metrics"
        ],
        "benchmark_metrics_target_filtered_outlier_trimmed": benchmark_target_filtered_outlier_trim_bundle[
            "trimmed_metrics"
        ],
        "benchmark_metrics_qc_qualified_target_filtered_outlier_trimmed": (
            benchmark_target_filtered_outlier_trim_qc_qualified_bundle["trimmed_metrics"]
        ),
        "batches": batch_rows,
    }
    return payload, batch_rows, job_rows, pair_rows, qc_qualified_pair_rows


def _write_ab_bind_plan_report_bundle(
    reports_dir: Path,
    payload: dict[str, Any],
    batch_rows: list[dict[str, Any]],
    job_rows: list[dict[str, Any]],
    pair_rows: list[dict[str, Any]],
    qc_qualified_pair_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    reports_dir = ensure_dir(reports_dir)
    bundle = dict(payload)
    bundle["reports_dir"] = str(reports_dir)
    benchmark_target_bundle = _benchmark_target_metrics_bundle(pair_rows)
    benchmark_target_qc_qualified_bundle = _benchmark_target_metrics_bundle(qc_qualified_pair_rows)
    active_alternate_job_rows = _active_alternate_job_rows(job_rows)
    active_alternate_job_fields = [
        "complex_id",
        "batch_id",
        "source_plan_root",
        "job_id",
        "mutation_group_id",
        "source_mutation",
        "mutation_tokens",
        "latest_stage",
        "latest_stage_state",
        "qc_status",
        "complex_leg_qc_status",
        "apo_leg_qc_status",
        "complex_repeat_spread_kcal_mol",
        "apo_repeat_spread_kcal_mol",
        "repeat_spread_legs",
        "primary_repeat_spread_leg",
        "diagnostic_family",
        "diagnostic_code",
        "diagnostic_detail",
        "validation_failure_category",
        "validation_failure_detail",
        "current_invalid_mutate_output",
        "current_invalid_mutate_output_code",
        "current_invalid_mutate_output_detail",
        "benchmark_qc_qualified",
        "analyzable",
        "resumable",
        "ddg_ready",
        "ddg_kcal_mol",
        "experimental_ddg_kcal_mol",
        "ddg_error_kcal_mol",
        "abs_ddg_error_kcal_mol",
        "ddg_bar_stderr_kcal_mol",
        "max_bar_stderr_kcal_mol",
        "active_alternate_candidate_count",
        "active_alternate_source_plan_roots",
        "active_alternate_stage_states",
        "active_alternate_current_source_plan_root",
        "active_alternate_current_batch_id",
        "active_alternate_current_latest_stage",
        "active_alternate_current_latest_stage_state",
        "active_alternate_current_equilibrate_started_repeats",
        "active_alternate_current_equilibrate_completed_repeats",
        "active_alternate_current_equilibrate_total_repeats",
        "active_alternate_current_sample_started_windows",
        "active_alternate_current_sample_completed_windows",
        "active_alternate_current_sample_total_windows",
        "active_alternate_current_sample_active_leg",
        "active_alternate_current_sample_active_repeat_id",
        "active_alternate_current_sample_active_lambda_id",
        "active_alternate_current_sample_active_lambda_index",
        "active_alternate_current_sample_active_phase",
        "active_alternate_current_sample_active_window",
        "equilibrate_started_repeats",
        "equilibrate_completed_repeats",
        "equilibrate_total_repeats",
        "sample_started_windows",
        "sample_completed_windows",
        "sample_total_windows",
        "sample_active_leg",
        "sample_active_repeat_id",
        "sample_active_lambda_id",
        "sample_active_lambda_index",
        "sample_active_phase",
        "sample_active_window",
    ]
    write_json(reports_dir / "plan_summary.json", bundle)
    write_yaml(reports_dir / "plan_summary.yml", bundle)
    write_json(reports_dir / "benchmark_metrics.json", bundle["benchmark_metrics"])
    write_yaml(reports_dir / "benchmark_metrics.yml", bundle["benchmark_metrics"])
    write_json(reports_dir / "benchmark_metrics_qc_qualified.json", bundle["benchmark_metrics_qc_qualified"])
    write_yaml(reports_dir / "benchmark_metrics_qc_qualified.yml", bundle["benchmark_metrics_qc_qualified"])
    write_json(reports_dir / "benchmark_metrics_target_filtered.json", bundle["benchmark_metrics_target_filtered"])
    write_yaml(reports_dir / "benchmark_metrics_target_filtered.yml", bundle["benchmark_metrics_target_filtered"])
    write_json(reports_dir / "benchmark_metrics_outlier_trimmed.json", bundle["benchmark_metrics_outlier_trimmed"])
    write_yaml(reports_dir / "benchmark_metrics_outlier_trimmed.yml", bundle["benchmark_metrics_outlier_trimmed"])
    write_json(
        reports_dir / "benchmark_metrics_qc_qualified_outlier_trimmed.json",
        bundle["benchmark_metrics_qc_qualified_outlier_trimmed"],
    )
    write_yaml(
        reports_dir / "benchmark_metrics_qc_qualified_outlier_trimmed.yml",
        bundle["benchmark_metrics_qc_qualified_outlier_trimmed"],
    )
    write_json(
        reports_dir / "benchmark_metrics_target_filtered_outlier_trimmed.json",
        bundle["benchmark_metrics_target_filtered_outlier_trimmed"],
    )
    write_yaml(
        reports_dir / "benchmark_metrics_target_filtered_outlier_trimmed.yml",
        bundle["benchmark_metrics_target_filtered_outlier_trimmed"],
    )
    write_json(
        reports_dir / "benchmark_metrics_qc_qualified_target_filtered_outlier_trimmed.json",
        bundle["benchmark_metrics_qc_qualified_target_filtered_outlier_trimmed"],
    )
    write_yaml(
        reports_dir / "benchmark_metrics_qc_qualified_target_filtered_outlier_trimmed.yml",
        bundle["benchmark_metrics_qc_qualified_target_filtered_outlier_trimmed"],
    )
    write_json(
        reports_dir / "benchmark_metrics_qc_qualified_target_filtered.json",
        bundle["benchmark_metrics_qc_qualified_target_filtered"],
    )
    write_yaml(
        reports_dir / "benchmark_metrics_qc_qualified_target_filtered.yml",
        bundle["benchmark_metrics_qc_qualified_target_filtered"],
    )
    write_csv_rows(
        reports_dir / "plan_batches.csv",
        batch_rows,
        [
            "complex_id",
            "batch_id",
            "batch_dir",
            "source_plan_roots",
            "job_count",
            "ready_job_count",
            "analyzable_job_count",
            "resumable_job_count",
            "active_alternate_job_count",
            "active_alternate_ready_job_count",
            "paired_job_count",
            "qc_qualified_pair_count",
            "not_started_count",
            "qc_pass_count",
            "qc_warning_count",
            "qc_fail_count",
            "qc_not_evaluated_count",
            "running_sample_job_count",
            "running_sample_started_windows",
            "running_sample_completed_windows",
            "running_sample_total_windows",
            "running_equilibrate_job_count",
            "running_equilibrate_started_repeats",
            "running_equilibrate_completed_repeats",
            "running_equilibrate_total_repeats",
            "current_invalid_mutate_output_job_count",
            "pearson_r",
            "spearman_rho",
            "rmse_kcal_mol",
            "mae_kcal_mol",
            "sign_accuracy",
            "qualified_pearson_r",
            "qualified_spearman_rho",
            "qualified_rmse_kcal_mol",
            "qualified_mae_kcal_mol",
            "qualified_sign_accuracy",
        ],
    )
    write_csv_rows(
        reports_dir / "plan_jobs.csv",
        job_rows,
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "protocol_preset",
            "protocol_repeats",
            "protocol_lambda_windows",
            "protocol_production_ps",
            "latest_stage",
            "latest_stage_state",
            "stage_count",
            "analyzable",
            "resumable",
            "alternate_candidate_count",
            "has_active_alternate_candidate",
            "active_alternate_candidate_count",
            "active_alternate_source_plan_roots",
            "active_alternate_stage_states",
            "active_alternate_current_source_plan_root",
            "active_alternate_current_batch_id",
            "active_alternate_current_latest_stage",
            "active_alternate_current_latest_stage_state",
            "active_alternate_current_equilibrate_started_repeats",
            "active_alternate_current_equilibrate_completed_repeats",
            "active_alternate_current_equilibrate_total_repeats",
            "active_alternate_current_sample_started_windows",
            "active_alternate_current_sample_completed_windows",
            "active_alternate_current_sample_total_windows",
            "active_alternate_current_sample_active_leg",
            "active_alternate_current_sample_active_repeat_id",
            "active_alternate_current_sample_active_lambda_id",
            "active_alternate_current_sample_active_lambda_index",
            "active_alternate_current_sample_active_phase",
            "active_alternate_current_sample_active_window",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "ddg_kcal_mol",
            "ddg_ready",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "complex_leg_qc_status",
            "apo_leg_qc_status",
            "complex_repeat_spread_kcal_mol",
            "apo_repeat_spread_kcal_mol",
            "repeat_spread_legs",
            "primary_repeat_spread_leg",
            "diagnostic_family",
            "diagnostic_code",
            "diagnostic_detail",
            "validation_failure_category",
            "validation_failure_detail",
            "current_invalid_mutate_output",
            "current_invalid_mutate_output_code",
            "current_invalid_mutate_output_detail",
            "benchmark_qc_qualified",
            "equilibrate_started_repeats",
            "equilibrate_completed_repeats",
            "equilibrate_total_repeats",
            "sample_started_windows",
            "sample_completed_windows",
            "sample_total_windows",
            "sample_active_leg",
            "sample_active_repeat_id",
            "sample_active_lambda_id",
            "sample_active_lambda_index",
            "sample_active_phase",
            "sample_active_window",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_ddg_error_kcal_mol",
            "source_mutation",
            "mutation_tokens",
        ],
    )
    write_csv_rows(
        reports_dir / "active_alternate_jobs.csv",
        [{field: row.get(field) for field in active_alternate_job_fields} for row in active_alternate_job_rows],
        active_alternate_job_fields,
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs.csv",
        pair_rows,
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_qc_qualified.csv",
        qc_qualified_pair_rows,
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    target_metric_fields = [
        "complex_id",
        "paired_job_count",
        "pearson_r",
        "spearman_rho",
        "rmse_kcal_mol",
        "mae_kcal_mol",
        "sign_accuracy",
        "auc_strong_effect",
        "min_abs_error_kcal_mol",
        "max_abs_error_kcal_mol",
        "mean_abs_error_kcal_mol",
        "sign_mismatch_count",
        "all_pairs_sign_mismatched",
        "all_pairs_above_abs_error_threshold",
        "systematically_poor_abs_error_threshold_kcal_mol",
        "systematically_poor_min_pair_count",
        "systematically_poor_max_pearson_r",
        "systematically_poor_max_sign_accuracy",
        "systematically_poor_min_leave_one_out_pearson_gain",
        "systematically_poor_abs_error_target",
        "systematically_poor_correlation_target",
        "systematically_poor_target_reason",
        "target_exclusion_iteration",
        "overall_pearson_r_at_exclusion",
        "leave_one_out_pearson_r_at_exclusion",
        "leave_one_out_pearson_gain_at_exclusion",
        "systematically_poor_target",
        "excluded_from_target_filtered_metrics",
        "leave_one_out_paired_job_count",
        "leave_one_out_pearson_r",
        "leave_one_out_spearman_rho",
        "leave_one_out_rmse_kcal_mol",
        "leave_one_out_mae_kcal_mol",
        "leave_one_out_sign_accuracy",
        "leave_one_out_auc_strong_effect",
    ]
    outlier_trim_target_metric_fields = [
        "complex_id",
        "outlier_trim_method",
        "original_paired_job_count",
        "trimmed_paired_job_count",
        "removed_pair_count",
        "removed_fraction",
        "q1_abs_error_kcal_mol",
        "q3_abs_error_kcal_mol",
        "iqr_abs_error_kcal_mol",
        "threshold_abs_error_kcal_mol",
        "removed_job_ids",
        "removed_job_ids_text",
        "trimmed_pearson_r",
        "trimmed_spearman_rho",
        "trimmed_rmse_kcal_mol",
        "trimmed_mae_kcal_mol",
        "trimmed_sign_accuracy",
        "trimmed_auc_strong_effect",
    ]
    write_csv_rows(
        reports_dir / "benchmark_target_metrics.csv",
        bundle.get("benchmark_target_metrics", []),
        target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "benchmark_target_metrics_qc_qualified.csv",
        bundle.get("benchmark_target_metrics_qc_qualified", []),
        target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "benchmark_target_outlier_trim_metrics.csv",
        bundle.get("benchmark_target_outlier_trim_metrics", []),
        outlier_trim_target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "benchmark_target_outlier_trim_metrics_qc_qualified.csv",
        bundle.get("benchmark_target_outlier_trim_metrics_qc_qualified", []),
        outlier_trim_target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_target_filtered.csv",
        benchmark_target_bundle["filtered_pair_rows"],
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    benchmark_outlier_trim_bundle = _benchmark_outlier_trim_bundle(pair_rows)
    benchmark_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(qc_qualified_pair_rows)
    benchmark_target_filtered_outlier_trim_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_bundle["filtered_pair_rows"]
    )
    benchmark_target_filtered_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_qc_qualified_bundle["filtered_pair_rows"]
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_outlier_trimmed.csv",
        benchmark_outlier_trim_bundle["trimmed_pair_rows"],
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_qc_qualified_outlier_trimmed.csv",
        benchmark_outlier_trim_qc_qualified_bundle["trimmed_pair_rows"],
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_target_filtered_outlier_trimmed.csv",
        benchmark_target_filtered_outlier_trim_bundle["trimmed_pair_rows"],
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_qc_qualified_target_filtered_outlier_trimmed.csv",
        benchmark_target_filtered_outlier_trim_qc_qualified_bundle["trimmed_pair_rows"],
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    write_csv_rows(
        reports_dir / "benchmark_pairs_qc_qualified_target_filtered.csv",
        benchmark_target_qc_qualified_bundle["filtered_pair_rows"],
        [
            "complex_id",
            "batch_id",
            "source_plan_root",
            "job_id",
            "mutation_group_id",
            "complex_delta_g_kcal_mol",
            "apo_delta_g_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
            "abs_error_kcal_mol",
            "ddg_bar_stderr_kcal_mol",
            "max_bar_stderr_kcal_mol",
            "qc_status",
            "benchmark_qc_qualified",
        ],
    )
    return bundle


_REPORT_STAGE_ORDER = {
    "": 0,
    "not_started": 0,
    "ingest": 1,
    "prepare": 2,
    "mutate": 3,
    "build_legs": 4,
    "equilibrate": 5,
    "sample": 6,
    "bar": 7,
    "qc": 8,
    "report": 9,
}

_REPORT_STATE_ORDER = {
    "not_started": 0,
    "failed": 1,
    "blocked_input": 1,
    "blocked_external": 1,
    "stale_running": 1,
    "planned": 2,
    "running": 3,
    "completed": 4,
}

_REPORT_QC_ORDER = {
    "not_started": 0,
    "fail": 1,
    "not_evaluated": 2,
    "warning": 3,
    "pass": 4,
}

_MERGED_LIVE_STATE_ORDER = {
    "not_started": 0,
    "planned": 0,
    "blocked_input": 0,
    "blocked_external": 0,
    "failed": 0,
    "stale_running": 0,
    "running": 1,
    "completed": 1,
}
_MERGED_ACTIVE_ALTERNATE_STATES = {"running", "stale_running"}
_ACTIVE_ALTERNATE_READY_HOTSPOT_LIMIT = 10
_ACTIVE_ALTERNATE_SAMPLE_PHASE_ORDER = {
    "": 0,
    "started": 1,
    "pre_relax": 2,
    "pre_md": 3,
    "md": 4,
    "completed": 5,
}


def _job_error_metrics(row: dict[str, Any]) -> tuple[float | None, float | None, float | None, float | None]:
    predicted_ddg = _safe_float(row.get("ddg_kcal_mol"))
    experimental_ddg = _safe_float(row.get("experimental_ddg_kcal_mol"))
    error = _safe_float(row.get("ddg_error_kcal_mol"))
    if error is None and predicted_ddg is not None and experimental_ddg is not None and _coerce_bool(row.get("ddg_ready")):
        error = predicted_ddg - experimental_ddg
    return predicted_ddg, experimental_ddg, error, (abs(error) if error is not None else None)


def _annotate_job_error_metrics(row: dict[str, Any]) -> dict[str, Any]:
    annotated = dict(row)
    _predicted_ddg, _experimental_ddg, error, abs_error = _job_error_metrics(annotated)
    if error is not None and annotated.get("ddg_error_kcal_mol") in ("", None):
        annotated["ddg_error_kcal_mol"] = error
    annotated["abs_ddg_error_kcal_mol"] = abs_error
    return annotated


_VALIDATION_FAILURE_CATEGORY_ORDER = {
    "benchmark_qc_qualified": 0,
    "qc_sampling_issue": 1,
    "running_execution": 2,
    "pending_execution": 3,
    "input_structure_issue": 4,
    "mutate_setup_issue": 5,
    "downstream_input_issue": 6,
    "external_dependency_issue": 7,
    "execution_failure": 8,
    "completed_unqualified": 9,
    "unknown": 10,
}

_VALIDATION_QC_SAMPLING_CODES = {
    "qc_bar_stderr",
    "qc_fail",
    "qc_histogram_missing",
    "qc_histogram_unparsed",
    "qc_low_overlap",
    "qc_missing_bar_output",
    "qc_missing_dhdl",
    "qc_repeat_spread",
    "qc_warning",
}
_VALIDATION_INPUT_STRUCTURE_CODES = {
    "prepare_blocked_input",
    "input_backbone_incomplete",
    "input_inter_residue_clash",
    "input_intra_residue_clash",
}
_VALIDATION_DOWNSTREAM_INPUT_CODES = {
    "equilibrate_blocked_input",
    "equilibrate_invalid_processed_gro",
    "equilibrate_missing_hybrid_topology",
    "equilibrate_missing_processed_gro",
    "sample_blocked_input",
    "sample_missing_npt_gro",
    "sample_missing_repeat_topology",
}
_VALIDATION_EXTERNAL_SUFFIXES = (
    "_blocked_external",
    "_gmxlib_unavailable",
    "_pmx_unavailable",
    "_forcefield_unavailable",
)
_VALIDATION_FAILURE_SUFFIX = "_failed"


def _validation_failure_taxonomy_category(row: dict[str, Any]) -> tuple[str, str]:
    diagnostic_family = str(row.get("diagnostic_family") or "").strip()
    diagnostic_code = str(row.get("diagnostic_code") or "").strip()
    latest_stage_state = str(row.get("latest_stage_state") or "").strip()
    current_invalid_mutate_output = _coerce_bool(row.get("current_invalid_mutate_output"))
    current_invalid_mutate_output_code = str(row.get("current_invalid_mutate_output_code") or "").strip()

    if _coerce_bool(row.get("benchmark_qc_qualified")):
        return "benchmark_qc_qualified", "qc_pass"
    if diagnostic_family == "running" or latest_stage_state in {"running", "stale_running"}:
        return "running_execution", diagnostic_code or diagnostic_family or "running"
    if current_invalid_mutate_output:
        return "mutate_setup_issue", current_invalid_mutate_output_code or "current_invalid_mutate_output"
    if diagnostic_code in _VALIDATION_QC_SAMPLING_CODES or diagnostic_family == "qc":
        return "qc_sampling_issue", diagnostic_code or diagnostic_family or "qc"
    if diagnostic_code in _VALIDATION_INPUT_STRUCTURE_CODES or diagnostic_family == "input_structure":
        return "input_structure_issue", diagnostic_code or diagnostic_family or "input_structure"
    if diagnostic_code.startswith("mutate_") or diagnostic_family == "mutate_setup":
        return "mutate_setup_issue", diagnostic_code or diagnostic_family or "mutate_setup"
    if diagnostic_code.endswith(_VALIDATION_FAILURE_SUFFIX):
        return "execution_failure", diagnostic_code
    if diagnostic_code.endswith(_VALIDATION_EXTERNAL_SUFFIXES):
        return "external_dependency_issue", diagnostic_code
    if diagnostic_code in _VALIDATION_DOWNSTREAM_INPUT_CODES or diagnostic_code.endswith("_blocked_input"):
        return "downstream_input_issue", diagnostic_code or diagnostic_family or "blocked_input"
    if (
        diagnostic_code == "not_started"
        or diagnostic_code.startswith("pending_")
        or diagnostic_family in {"pending", "not_started"}
    ):
        return "pending_execution", diagnostic_code or diagnostic_family or "pending"
    if diagnostic_family == "completed" or diagnostic_code == "reported":
        return "completed_unqualified", diagnostic_code or diagnostic_family or "completed"
    return "unknown", diagnostic_code or diagnostic_family or "unknown"


def _annotate_validation_failure_taxonomy(row: dict[str, Any]) -> dict[str, Any]:
    annotated = dict(row)
    category, detail = _validation_failure_taxonomy_category(annotated)
    annotated["validation_failure_category"] = category
    annotated["validation_failure_detail"] = detail
    return annotated


def _validation_failure_sort_key(row: dict[str, Any]) -> tuple[float | str, ...]:
    return (
        0 if _coerce_bool(row.get("ddg_ready")) else 1,
        0 if _safe_float(row.get("abs_ddg_error_kcal_mol")) is not None else 1,
        -(_safe_float(row.get("abs_ddg_error_kcal_mol")) or 0.0),
        -float(_REPORT_STAGE_ORDER.get(str(row.get("latest_stage", "")), 0)),
        -float(_REPORT_STATE_ORDER.get(str(row.get("latest_stage_state", "not_started")), 0)),
        str(row.get("complex_id", "")),
        str(row.get("job_id", "")),
        str(row.get("source_plan_root", "")),
    )


def _build_validation_failure_taxonomy(job_rows: list[dict[str, Any]]) -> dict[str, Any]:
    normalized_rows = [
        _annotate_validation_failure_taxonomy(_annotate_job_error_metrics(row))
        for row in job_rows
    ]
    category_counts: Counter[str] = Counter()
    ready_job_counts: Counter[str] = Counter()
    running_job_counts: Counter[str] = Counter()
    resumable_job_counts: Counter[str] = Counter()

    for row in normalized_rows:
        category = str(row.get("validation_failure_category") or "unknown")
        category_counts.update([category])
        if _coerce_bool(row.get("ddg_ready")):
            ready_job_counts.update([category])
        if str(row.get("latest_stage_state", "not_started")) in {"running", "stale_running"}:
            running_job_counts.update([category])
        if _coerce_bool(row.get("resumable")):
            resumable_job_counts.update([category])

    sorted_categories = sorted(
        category_counts,
        key=lambda item: (_VALIDATION_FAILURE_CATEGORY_ORDER.get(item, 999), item),
    )
    categories_payload: list[dict[str, Any]] = []
    for category in sorted_categories:
        category_rows = [row for row in normalized_rows if row.get("validation_failure_category") == category]
        category_rows = sorted(category_rows, key=_validation_failure_sort_key)
        categories_payload.append(
            {
                "category": category,
                "job_count": category_counts[category],
                "ready_job_count": ready_job_counts.get(category, 0),
                "running_job_count": running_job_counts.get(category, 0),
                "resumable_job_count": resumable_job_counts.get(category, 0),
                "sample_jobs": [
                    {
                        "complex_id": row.get("complex_id", ""),
                        "job_id": row.get("job_id", ""),
                        "diagnostic_code": row.get("diagnostic_code", ""),
                        "validation_failure_detail": row.get("validation_failure_detail", ""),
                        "latest_stage": row.get("latest_stage", ""),
                        "latest_stage_state": row.get("latest_stage_state", ""),
                        "ddg_ready": _coerce_bool(row.get("ddg_ready")),
                        "benchmark_qc_qualified": _coerce_bool(row.get("benchmark_qc_qualified")),
                        "abs_ddg_error_kcal_mol": _safe_float(row.get("abs_ddg_error_kcal_mol")),
                    }
                    for row in category_rows[:3]
                ],
            }
        )

    hotspots = [
        row
        for row in sorted(normalized_rows, key=_validation_failure_sort_key)
        if str(row.get("validation_failure_category") or "unknown") != "benchmark_qc_qualified"
    ][:12]
    return {
        "counts": {category: category_counts[category] for category in sorted_categories},
        "ready_job_counts": {category: ready_job_counts.get(category, 0) for category in sorted_categories},
        "running_job_counts": {category: running_job_counts.get(category, 0) for category in sorted_categories},
        "resumable_job_counts": {category: resumable_job_counts.get(category, 0) for category in sorted_categories},
        "categories": categories_payload,
        "hotspots": [
            {
                "category": row.get("validation_failure_category", ""),
                "complex_id": row.get("complex_id", ""),
                "job_id": row.get("job_id", ""),
                "diagnostic_code": row.get("diagnostic_code", ""),
                "validation_failure_detail": row.get("validation_failure_detail", ""),
                "latest_stage": row.get("latest_stage", ""),
                "latest_stage_state": row.get("latest_stage_state", ""),
                "ddg_ready": _coerce_bool(row.get("ddg_ready")),
                "benchmark_qc_qualified": _coerce_bool(row.get("benchmark_qc_qualified")),
                "abs_ddg_error_kcal_mol": _safe_float(row.get("abs_ddg_error_kcal_mol")),
            }
            for row in hotspots
        ],
    }


def _active_alternate_job_rows(job_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    active_rows = [
        _annotate_job_error_metrics(row)
        for row in job_rows
        if _coerce_bool(row.get("has_active_alternate_candidate"))
    ]
    return sorted(
        active_rows,
        key=lambda row: (
            0 if _coerce_bool(row.get("ddg_ready")) else 1,
            0 if _safe_float(row.get("abs_ddg_error_kcal_mol")) is not None else 1,
            -(_safe_float(row.get("abs_ddg_error_kcal_mol")) or 0.0),
            -float(_REPORT_STAGE_ORDER.get(str(row.get("latest_stage", "")), 0)),
            -float(_REPORT_STATE_ORDER.get(str(row.get("latest_stage_state", "not_started")), 0)),
            str(row.get("complex_id", "")),
            str(row.get("job_id", "")),
            str(row.get("source_plan_root", "")),
        ),
    )


def _active_alternate_ready_hotspots(job_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    hotspots: list[dict[str, Any]] = []
    for row in _active_alternate_job_rows(job_rows):
        if not _coerce_bool(row.get("ddg_ready")):
            continue
        hotspots.append(
            {
                "complex_id": row.get("complex_id", ""),
                "batch_id": row.get("batch_id", ""),
                "job_id": row.get("job_id", ""),
                "mutation_group_id": row.get("mutation_group_id", ""),
                "source_mutation": row.get("source_mutation", ""),
                "source_plan_root": row.get("source_plan_root", ""),
                "qc_status": row.get("qc_status", ""),
                "benchmark_qc_qualified": _coerce_bool(row.get("benchmark_qc_qualified")),
                "ddg_kcal_mol": _safe_float(row.get("ddg_kcal_mol")),
                "experimental_ddg_kcal_mol": _safe_float(row.get("experimental_ddg_kcal_mol")),
                "ddg_error_kcal_mol": _safe_float(row.get("ddg_error_kcal_mol")),
                "abs_ddg_error_kcal_mol": _safe_float(row.get("abs_ddg_error_kcal_mol")),
                "ddg_bar_stderr_kcal_mol": _safe_float(row.get("ddg_bar_stderr_kcal_mol")),
                "active_alternate_candidate_count": int(_safe_float(row.get("active_alternate_candidate_count")) or 0),
                "active_alternate_source_plan_roots": row.get("active_alternate_source_plan_roots", ""),
                "active_alternate_stage_states": row.get("active_alternate_stage_states", ""),
                "active_alternate_current_source_plan_root": row.get("active_alternate_current_source_plan_root", ""),
                "active_alternate_current_batch_id": row.get("active_alternate_current_batch_id", ""),
                "active_alternate_current_latest_stage": row.get("active_alternate_current_latest_stage", ""),
                "active_alternate_current_latest_stage_state": row.get("active_alternate_current_latest_stage_state", ""),
                "active_alternate_current_sample_active_phase": row.get(
                    "active_alternate_current_sample_active_phase", ""
                ),
                "active_alternate_current_sample_active_window": row.get(
                    "active_alternate_current_sample_active_window", ""
                ),
            }
        )
        if len(hotspots) >= _ACTIVE_ALTERNATE_READY_HOTSPOT_LIMIT:
            break
    return hotspots


def _merged_job_identity(row: dict[str, Any]) -> str:
    mutation_group_id = str(row.get("mutation_group_id", "")).strip()
    if mutation_group_id:
        return mutation_group_id
    return f"{row.get('complex_id', '')}::{row.get('job_id', '')}"


def _protocol_sampling_priority(row: dict[str, Any]) -> tuple[float, float, float, float]:
    repeats = max(_safe_float(row.get("protocol_repeats")) or 0.0, 0.0)
    lambda_windows = max(_safe_float(row.get("protocol_lambda_windows")) or 0.0, 0.0)
    production_ps = max(_safe_float(row.get("protocol_production_ps")) or 0.0, 0.0)
    total_effort = repeats * lambda_windows * production_ps
    return total_effort, repeats, lambda_windows, production_ps


def _merged_job_repeat_spread_priority(row: dict[str, Any]) -> float:
    spreads = [
        value
        for value in (
            _safe_float(row.get("complex_repeat_spread_kcal_mol")),
            _safe_float(row.get("apo_repeat_spread_kcal_mol")),
        )
        if value is not None
    ]
    if not spreads:
        return -1_000_000.0
    return -max(spreads)


def _active_alternate_priority(row: dict[str, Any], *, root_priority: int) -> tuple[float, ...]:
    latest_stage_state = str(row.get("latest_stage_state", "not_started"))
    protocol_total_effort, protocol_repeats, protocol_lambda_windows, protocol_production_ps = _protocol_sampling_priority(
        row
    )
    return (
        float(_REPORT_STATE_ORDER.get(latest_stage_state, 0)),
        float(_REPORT_STAGE_ORDER.get(str(row.get("latest_stage", "")), 0)),
        float(_safe_float(row.get("sample_completed_windows")) or 0.0),
        float(_ACTIVE_ALTERNATE_SAMPLE_PHASE_ORDER.get(str(row.get("sample_active_phase", "")), 0)),
        float(_safe_float(row.get("sample_active_lambda_index")) or -1.0),
        float(_safe_float(row.get("sample_started_windows")) or 0.0),
        float(_safe_float(row.get("equilibrate_completed_repeats")) or 0.0),
        float(_safe_float(row.get("equilibrate_started_repeats")) or 0.0),
        protocol_total_effort,
        protocol_repeats,
        protocol_lambda_windows,
        protocol_production_ps,
        float(-root_priority),
    )


def _best_active_alternate_row(
    rows: list[dict[str, Any]],
    *,
    root_priorities: dict[str, int],
) -> dict[str, Any] | None:
    best_row: dict[str, Any] | None = None
    best_priority: tuple[float, ...] | None = None
    for row in rows:
        source_plan_root = str(row.get("source_plan_root", ""))
        priority = _active_alternate_priority(
            row,
            root_priority=root_priorities.get(source_plan_root, len(root_priorities)),
        )
        if best_priority is None or priority > best_priority:
            best_row = row
            best_priority = priority
    return best_row


def _merged_job_priority(row: dict[str, Any], *, root_priority: int) -> tuple[float, ...]:
    ddg_bar_stderr = _safe_float(row.get("ddg_bar_stderr_kcal_mol"))
    latest_stage_state = str(row.get("latest_stage_state", "not_started"))
    repeat_spread_priority = _merged_job_repeat_spread_priority(row)
    protocol_total_effort, protocol_repeats, protocol_lambda_windows, protocol_production_ps = _protocol_sampling_priority(
        row
    )
    return (
        1.0 if row.get("benchmark_qc_qualified") else 0.0,
        1.0 if row.get("ddg_ready") else 0.0,
        float(_REPORT_QC_ORDER.get(str(row.get("qc_status", "not_started")), 0)),
        1.0 if row.get("analyzable") else 0.0,
        float(_MERGED_LIVE_STATE_ORDER.get(latest_stage_state, 0)),
        float(_REPORT_STAGE_ORDER.get(str(row.get("latest_stage", "")), 0)),
        float(_REPORT_STATE_ORDER.get(latest_stage_state, 0)),
        # Sampling depth outranks spread/stderr: among rows at the same QC and
        # completion level, the deepest-protocol value must win so that shallow
        # (e.g. quick-preset) results can never shadow deeper rescue reruns.
        protocol_total_effort,
        protocol_repeats,
        protocol_lambda_windows,
        protocol_production_ps,
        repeat_spread_priority,
        -(ddg_bar_stderr if ddg_bar_stderr is not None else 1_000_000.0),
        1.0 if row.get("resumable") else 0.0,
        float(-root_priority),
    )


def _merge_ab_bind_job_rows(
    candidates: list[dict[str, Any]],
    *,
    root_priorities: dict[str, int],
) -> list[dict[str, Any]]:
    grouped_candidates: dict[str, list[dict[str, Any]]] = {}
    for row in candidates:
        grouped_candidates.setdefault(_merged_job_identity(row), []).append(dict(row))

    winners: list[dict[str, Any]] = []
    for rows in grouped_candidates.values():
        winner_index = 0
        winner_row = rows[0]
        winner_priority = _merged_job_priority(
            winner_row,
            root_priority=root_priorities.get(str(winner_row.get("source_plan_root", "")), len(root_priorities)),
        )
        for index, row in enumerate(rows[1:], start=1):
            source_plan_root = str(row.get("source_plan_root", ""))
            priority = _merged_job_priority(row, root_priority=root_priorities.get(source_plan_root, len(root_priorities)))
            if priority > winner_priority:
                winner_index = index
                winner_row = row
                winner_priority = priority

        alternates = [row for index, row in enumerate(rows) if index != winner_index]
        active_alternates = [
            row
            for row in alternates
            if str(row.get("latest_stage_state", "not_started")) in _MERGED_ACTIVE_ALTERNATE_STATES
        ]
        annotated = dict(winner_row)
        annotated["alternate_candidate_count"] = len(alternates)
        annotated["has_active_alternate_candidate"] = bool(active_alternates)
        annotated["active_alternate_candidate_count"] = len(active_alternates)
        annotated["active_alternate_source_plan_roots"] = ",".join(
            sorted(
                {
                    str(row.get("source_plan_root", "")).strip()
                    for row in active_alternates
                    if str(row.get("source_plan_root", "")).strip()
                }
            )
        )
        annotated["active_alternate_stage_states"] = ",".join(
            sorted(
                {
                    f"{str(row.get('latest_stage', '')).strip()}:{str(row.get('latest_stage_state', '')).strip()}"
                    for row in active_alternates
                    if str(row.get("latest_stage", "")).strip() or str(row.get("latest_stage_state", "")).strip()
                }
            )
        )
        representative_active_alternate = _best_active_alternate_row(active_alternates, root_priorities=root_priorities)
        if representative_active_alternate is not None:
            annotated["active_alternate_current_source_plan_root"] = representative_active_alternate.get(
                "source_plan_root", ""
            )
            annotated["active_alternate_current_batch_id"] = representative_active_alternate.get("batch_id", "")
            annotated["active_alternate_current_latest_stage"] = representative_active_alternate.get("latest_stage", "")
            annotated["active_alternate_current_latest_stage_state"] = representative_active_alternate.get(
                "latest_stage_state", ""
            )
            annotated["active_alternate_current_equilibrate_started_repeats"] = representative_active_alternate.get(
                "equilibrate_started_repeats", ""
            )
            annotated["active_alternate_current_equilibrate_completed_repeats"] = representative_active_alternate.get(
                "equilibrate_completed_repeats", ""
            )
            annotated["active_alternate_current_equilibrate_total_repeats"] = representative_active_alternate.get(
                "equilibrate_total_repeats", ""
            )
            annotated["active_alternate_current_sample_started_windows"] = representative_active_alternate.get(
                "sample_started_windows", ""
            )
            annotated["active_alternate_current_sample_completed_windows"] = representative_active_alternate.get(
                "sample_completed_windows", ""
            )
            annotated["active_alternate_current_sample_total_windows"] = representative_active_alternate.get(
                "sample_total_windows", ""
            )
            annotated["active_alternate_current_sample_active_leg"] = representative_active_alternate.get(
                "sample_active_leg", ""
            )
            annotated["active_alternate_current_sample_active_repeat_id"] = representative_active_alternate.get(
                "sample_active_repeat_id", ""
            )
            annotated["active_alternate_current_sample_active_lambda_id"] = representative_active_alternate.get(
                "sample_active_lambda_id", ""
            )
            annotated["active_alternate_current_sample_active_lambda_index"] = representative_active_alternate.get(
                "sample_active_lambda_index", ""
            )
            annotated["active_alternate_current_sample_active_phase"] = representative_active_alternate.get(
                "sample_active_phase", ""
            )
            annotated["active_alternate_current_sample_active_window"] = representative_active_alternate.get(
                "sample_active_window", ""
            )
        winners.append(annotated)
    return sorted(
        winners,
        key=lambda row: (
            str(row.get("complex_id", "")),
            str(row.get("job_id", "")),
            str(row.get("source_plan_root", "")),
        ),
    )


def _build_ab_bind_report_rows_from_jobs(
    job_rows: list[dict[str, Any]],
    *,
    plan_root: Path,
    benchmark_root: Path | None,
    spec_name: str,
    batch_metadata_by_complex: dict[str, dict[str, Any]],
    source_plan_roots: list[Path],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    overall_pearson_threshold = 0.6
    qc_counts: Counter[str] = Counter()
    latest_stage_counts: Counter[str] = Counter()
    latest_stage_name_counts: Counter[str] = Counter()
    diagnostic_family_counts: Counter[str] = Counter()
    diagnostic_code_counts: Counter[str] = Counter()
    batch_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    qc_qualified_pair_rows: list[dict[str, Any]] = []
    jobs_by_complex: dict[str, list[dict[str, Any]]] = {}

    normalized_jobs = sorted(
        (_annotate_validation_failure_taxonomy(_annotate_job_error_metrics(row)) for row in job_rows),
        key=lambda row: (
            str(row.get("complex_id", "")),
            str(row.get("job_id", "")),
            str(row.get("source_plan_root", "")),
        ),
    )

    for row in normalized_jobs:
        complex_id = str(row.get("complex_id", ""))
        jobs_by_complex.setdefault(complex_id, []).append(row)
        qc_counts.update([str(row.get("qc_status", "not_started"))])
        latest_stage_counts.update([str(row.get("latest_stage_state", "not_started")) or "not_started"])
        latest_stage_name_counts.update([str(row.get("latest_stage", "")) or "not_started"])
        diagnostic_family = str(row.get("diagnostic_family") or "").strip()
        diagnostic_code = str(row.get("diagnostic_code") or "").strip()
        if diagnostic_family:
            diagnostic_family_counts.update([diagnostic_family])
        if diagnostic_code:
            diagnostic_code_counts.update([diagnostic_code])

        experimental_ddg = _safe_float(row.get("experimental_ddg_kcal_mol"))
        predicted_ddg = _safe_float(row.get("ddg_kcal_mol"))
        error = _safe_float(row.get("ddg_error_kcal_mol"))
        if error is None and predicted_ddg is not None and experimental_ddg is not None and row.get("ddg_ready"):
            error = predicted_ddg - experimental_ddg
        if error is None:
            continue
        pair_row = {
            "complex_id": complex_id,
            "batch_id": row.get("batch_id", ""),
            "source_plan_root": row.get("source_plan_root", ""),
            "job_id": row.get("job_id", ""),
            "mutation_group_id": row.get("mutation_group_id", ""),
            "complex_delta_g_kcal_mol": _safe_float(row.get("complex_delta_g_kcal_mol")),
            "apo_delta_g_kcal_mol": _safe_float(row.get("apo_delta_g_kcal_mol")),
            "predicted_ddg_kcal_mol": predicted_ddg,
            "experimental_ddg_kcal_mol": experimental_ddg,
            "ddg_error_kcal_mol": error,
            "abs_error_kcal_mol": abs(error),
            "ddg_bar_stderr_kcal_mol": _safe_float(row.get("ddg_bar_stderr_kcal_mol")),
            "max_bar_stderr_kcal_mol": _safe_float(row.get("max_bar_stderr_kcal_mol")),
            "qc_status": row.get("qc_status", ""),
            "benchmark_qc_qualified": bool(row.get("benchmark_qc_qualified", False)),
        }
        pair_rows.append(pair_row)
        if pair_row["benchmark_qc_qualified"]:
            qc_qualified_pair_rows.append(pair_row)

    for complex_id in sorted(jobs_by_complex):
        complex_jobs = jobs_by_complex[complex_id]
        batch_meta = batch_metadata_by_complex.get(complex_id, {})
        running_sample_jobs = [
            job for job in complex_jobs if job.get("latest_stage") == "sample" and job.get("latest_stage_state") == "running"
        ]
        running_equilibrate_jobs = [
            job
            for job in complex_jobs
            if job.get("latest_stage") == "equilibrate" and job.get("latest_stage_state") == "running"
        ]
        current_invalid_mutate_output_jobs = [
            job for job in complex_jobs if _as_bool(str(job.get("current_invalid_mutate_output", "")))
        ]
        batch_qc_counts: Counter[str] = Counter(str(job.get("qc_status", "not_started")) for job in complex_jobs)
        batch_stage_counts: Counter[str] = Counter(
            str(job.get("latest_stage_state", "not_started")) or "not_started" for job in complex_jobs
        )
        batch_pair_rows = [row for row in pair_rows if row["complex_id"] == complex_id]
        batch_qc_qualified_pair_rows = [row for row in qc_qualified_pair_rows if row["complex_id"] == complex_id]
        batch_metrics = _benchmark_metrics_from_pairs(batch_pair_rows) if batch_pair_rows else {}
        batch_metrics_qc_qualified = (
            _benchmark_metrics_from_pairs(batch_qc_qualified_pair_rows) if batch_qc_qualified_pair_rows else {}
        )
        batch_rows.append(
            {
                "complex_id": complex_id,
                "batch_id": batch_meta.get("batch_id", complex_id),
                "batch_dir": batch_meta.get("batch_dir", ""),
                "source_plan_roots": ",".join(
                    sorted(
                        {
                            str(job.get("source_plan_root", "")).strip()
                            for job in complex_jobs
                            if str(job.get("source_plan_root", "")).strip()
                        }
                    )
                ),
                "job_count": len(complex_jobs),
                "ready_job_count": sum(1 for job in complex_jobs if job.get("ddg_ready")),
                "analyzable_job_count": sum(1 for job in complex_jobs if job.get("analyzable")),
                "resumable_job_count": sum(1 for job in complex_jobs if job.get("resumable")),
                "active_alternate_job_count": sum(
                    1 for job in complex_jobs if _as_bool(str(job.get("has_active_alternate_candidate", "")))
                ),
                "active_alternate_ready_job_count": sum(
                    1
                    for job in complex_jobs
                    if job.get("ddg_ready") and _as_bool(str(job.get("has_active_alternate_candidate", "")))
                ),
                "paired_job_count": len(batch_pair_rows),
                "qc_qualified_pair_count": len(batch_qc_qualified_pair_rows),
                "not_started_count": batch_stage_counts.get("not_started", 0),
                "qc_pass_count": batch_qc_counts.get("pass", 0),
                "qc_warning_count": batch_qc_counts.get("warning", 0),
                "qc_fail_count": batch_qc_counts.get("fail", 0),
                "qc_not_evaluated_count": batch_qc_counts.get("not_evaluated", 0),
                "running_sample_job_count": len(running_sample_jobs),
                "running_sample_started_windows": sum(job.get("sample_started_windows", 0) for job in running_sample_jobs),
                "running_sample_completed_windows": sum(
                    job.get("sample_completed_windows", 0) for job in running_sample_jobs
                ),
                "running_sample_total_windows": sum(job.get("sample_total_windows", 0) for job in running_sample_jobs),
                "running_equilibrate_job_count": len(running_equilibrate_jobs),
                "running_equilibrate_started_repeats": sum(
                    job.get("equilibrate_started_repeats", 0) for job in running_equilibrate_jobs
                ),
                "running_equilibrate_completed_repeats": sum(
                    job.get("equilibrate_completed_repeats", 0) for job in running_equilibrate_jobs
                ),
                "running_equilibrate_total_repeats": sum(
                    job.get("equilibrate_total_repeats", 0) for job in running_equilibrate_jobs
                ),
                "current_invalid_mutate_output_job_count": len(current_invalid_mutate_output_jobs),
                "pearson_r": batch_metrics.get("pearson_r"),
                "spearman_rho": batch_metrics.get("spearman_rho"),
                "rmse_kcal_mol": batch_metrics.get("rmse_kcal_mol"),
                "mae_kcal_mol": batch_metrics.get("mae_kcal_mol"),
                "sign_accuracy": batch_metrics.get("sign_accuracy"),
                "qualified_pearson_r": batch_metrics_qc_qualified.get("pearson_r"),
                "qualified_spearman_rho": batch_metrics_qc_qualified.get("spearman_rho"),
                "qualified_rmse_kcal_mol": batch_metrics_qc_qualified.get("rmse_kcal_mol"),
                "qualified_mae_kcal_mol": batch_metrics_qc_qualified.get("mae_kcal_mol"),
                "qualified_sign_accuracy": batch_metrics_qc_qualified.get("sign_accuracy"),
            }
        )

    benchmark_metrics = _benchmark_metrics_from_pairs(pair_rows) if pair_rows else {}
    benchmark_metrics_qc_qualified = (
        _benchmark_metrics_from_pairs(qc_qualified_pair_rows) if qc_qualified_pair_rows else {}
    )
    benchmark_target_bundle = _benchmark_target_metrics_bundle(pair_rows)
    benchmark_target_qc_qualified_bundle = _benchmark_target_metrics_bundle(qc_qualified_pair_rows)
    benchmark_outlier_trim_bundle = _benchmark_outlier_trim_bundle(pair_rows)
    benchmark_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(qc_qualified_pair_rows)
    benchmark_target_filtered_outlier_trim_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_bundle["filtered_pair_rows"]
    )
    benchmark_target_filtered_outlier_trim_qc_qualified_bundle = _benchmark_outlier_trim_bundle(
        benchmark_target_qc_qualified_bundle["filtered_pair_rows"]
    )
    overall_pearson_r = benchmark_metrics.get("pearson_r")
    target_filtered_pearson_r = benchmark_target_bundle["filtered_metrics"].get("pearson_r")
    running_sample_rows = [
        row for row in normalized_jobs if row.get("latest_stage") == "sample" and row.get("latest_stage_state") == "running"
    ]
    running_equilibrate_rows = [
        row
        for row in normalized_jobs
        if row.get("latest_stage") == "equilibrate" and row.get("latest_stage_state") == "running"
    ]
    payload = {
        "generated_at": utc_now(),
        "plan_root": str(plan_root),
        "reports_dir": "",
        "spec_name": spec_name,
        "benchmark_root": str(benchmark_root) if benchmark_root is not None else "",
        "selected_batch_count": len(batch_rows),
        "selected_job_count": len(normalized_jobs),
        "ddg_ready_count": sum(1 for row in normalized_jobs if row.get("ddg_ready")),
        "analyzable_job_count": sum(1 for row in normalized_jobs if row.get("analyzable")),
        "resumable_job_count": sum(1 for row in normalized_jobs if row.get("resumable")),
        "active_alternate_job_count": sum(
            1 for row in normalized_jobs if _as_bool(str(row.get("has_active_alternate_candidate", "")))
        ),
        "active_alternate_ready_job_count": sum(
            1
            for row in normalized_jobs
            if row.get("ddg_ready") and _as_bool(str(row.get("has_active_alternate_candidate", "")))
        ),
        "active_alternate_ready_hotspot_count": sum(
            1
            for row in normalized_jobs
            if row.get("ddg_ready") and _as_bool(str(row.get("has_active_alternate_candidate", "")))
        ),
        "active_alternate_ready_hotspots": _active_alternate_ready_hotspots(normalized_jobs),
        "paired_job_count": len(pair_rows),
        "qc_qualified_pair_count": len(qc_qualified_pair_rows),
        "qc_counts": dict(sorted(qc_counts.items())),
        "latest_stage_name_counts": dict(sorted(latest_stage_name_counts.items())),
        "latest_stage_state_counts": dict(sorted(latest_stage_counts.items())),
        "diagnostic_family_counts": dict(sorted(diagnostic_family_counts.items())),
        "diagnostic_code_counts": dict(sorted(diagnostic_code_counts.items())),
        "running_sample_job_count": len(running_sample_rows),
        "running_sample_started_windows": sum(row.get("sample_started_windows", 0) for row in running_sample_rows),
        "running_sample_completed_windows": sum(row.get("sample_completed_windows", 0) for row in running_sample_rows),
        "running_sample_total_windows": sum(row.get("sample_total_windows", 0) for row in running_sample_rows),
        "running_equilibrate_job_count": len(running_equilibrate_rows),
        "running_equilibrate_started_repeats": sum(
            row.get("equilibrate_started_repeats", 0) for row in running_equilibrate_rows
        ),
        "running_equilibrate_completed_repeats": sum(
            row.get("equilibrate_completed_repeats", 0) for row in running_equilibrate_rows
        ),
        "running_equilibrate_total_repeats": sum(
            row.get("equilibrate_total_repeats", 0) for row in running_equilibrate_rows
        ),
        "current_invalid_mutate_output_job_count": sum(
            1 for row in normalized_jobs if _as_bool(str(row.get("current_invalid_mutate_output", "")))
        ),
        "validation_failure_taxonomy": _build_validation_failure_taxonomy(normalized_jobs),
        "validation_gate": {
            "overall_pearson_r_threshold": overall_pearson_threshold,
            "overall_pearson_r": overall_pearson_r,
            "overall_pearson_r_passed": overall_pearson_r is not None and overall_pearson_r > overall_pearson_threshold,
            "target_filtered_pearson_r": target_filtered_pearson_r,
            "target_filtered_pearson_r_passed": target_filtered_pearson_r is not None
            and target_filtered_pearson_r > overall_pearson_threshold,
            "target_filtered_excluded_complex_ids": benchmark_target_bundle["excluded_target_ids"],
        },
        "benchmark_metrics": benchmark_metrics,
        "benchmark_metrics_qc_qualified": benchmark_metrics_qc_qualified,
        "benchmark_target_exclusion_policy": {
            "target_field": benchmark_target_bundle["target_field"],
            "systematically_poor_abs_error_threshold_kcal_mol": benchmark_target_bundle[
                "systematically_poor_abs_error_threshold_kcal_mol"
            ],
            "systematically_poor_min_pair_count": benchmark_target_bundle["systematically_poor_min_pair_count"],
            "systematically_poor_max_pearson_r": benchmark_target_bundle["systematically_poor_max_pearson_r"],
            "systematically_poor_max_sign_accuracy": benchmark_target_bundle[
                "systematically_poor_max_sign_accuracy"
            ],
            "systematically_poor_min_leave_one_out_pearson_gain": benchmark_target_bundle[
                "systematically_poor_min_leave_one_out_pearson_gain"
            ],
        },
        "benchmark_outlier_trim_policy": {
            "target_field": benchmark_outlier_trim_bundle["target_field"],
            "outlier_trim_method": benchmark_outlier_trim_bundle["outlier_trim_method"],
        },
        "benchmark_target_metrics": benchmark_target_bundle["target_metrics"],
        "benchmark_target_metrics_qc_qualified": benchmark_target_qc_qualified_bundle["target_metrics"],
        "benchmark_target_outlier_trim_metrics": benchmark_outlier_trim_bundle["target_metrics"],
        "benchmark_target_outlier_trim_metrics_qc_qualified": benchmark_outlier_trim_qc_qualified_bundle[
            "target_metrics"
        ],
        "benchmark_target_excluded_complex_ids": benchmark_target_bundle["excluded_target_ids"],
        "benchmark_target_excluded_complex_ids_qc_qualified": benchmark_target_qc_qualified_bundle[
            "excluded_target_ids"
        ],
        "benchmark_metrics_outlier_trimmed": benchmark_outlier_trim_bundle["trimmed_metrics"],
        "benchmark_metrics_qc_qualified_outlier_trimmed": benchmark_outlier_trim_qc_qualified_bundle[
            "trimmed_metrics"
        ],
        "benchmark_metrics_target_filtered": benchmark_target_bundle["filtered_metrics"],
        "benchmark_metrics_qc_qualified_target_filtered": benchmark_target_qc_qualified_bundle[
            "filtered_metrics"
        ],
        "benchmark_metrics_target_filtered_outlier_trimmed": benchmark_target_filtered_outlier_trim_bundle[
            "trimmed_metrics"
        ],
        "benchmark_metrics_qc_qualified_target_filtered_outlier_trimmed": (
            benchmark_target_filtered_outlier_trim_qc_qualified_bundle["trimmed_metrics"]
        ),
        "batches": batch_rows,
        "source_plan_roots": [str(path) for path in source_plan_roots],
    }
    return payload, batch_rows, normalized_jobs, pair_rows, qc_qualified_pair_rows


def _rescue_reason_codes(qc_report: dict[str, Any]) -> list[str]:
    reasons: list[str] = []
    if qc_report.get("ddg_repeat_range_kcal_mol") is not None and qc_report.get("max_repeat_delta_kcal_mol") is not None:
        if qc_report["ddg_repeat_range_kcal_mol"] > qc_report["max_repeat_delta_kcal_mol"]:
            reasons.append("repeat_spread")
    if qc_report.get("ddg_bar_stderr_kcal_mol") is not None and qc_report.get("max_bar_stderr_kcal_mol") is not None:
        if qc_report["ddg_bar_stderr_kcal_mol"] > qc_report["max_bar_stderr_kcal_mol"]:
            reasons.append("bar_stderr")

    overlap_threshold = qc_report.get("overlap_threshold")
    overlap_legs = qc_report.get("overlap_assessment", {}).get("legs", {})
    if overlap_threshold is not None:
        for payload in overlap_legs.values():
            score = payload.get("overlap_score_min")
            if score is not None and score < overlap_threshold:
                reasons.append("overlap")
                break

    return reasons


def _low_overlap_legs(qc_report: dict[str, Any]) -> list[str]:
    overlap_threshold = _safe_float(qc_report.get("overlap_threshold"))
    if overlap_threshold is None:
        return []

    overlap_legs = qc_report.get("overlap_assessment", {}).get("legs", {})
    if not isinstance(overlap_legs, dict):
        return []

    failing_legs: list[str] = []
    for leg_name in ("complex", "apo"):
        payload = overlap_legs.get(leg_name, {})
        if not isinstance(payload, dict):
            continue
        score = _safe_float(payload.get("overlap_score_min"))
        if score is not None and score < overlap_threshold:
            failing_legs.append(leg_name)
    return failing_legs


def _has_dominant_primary_repeat_spread_leg(qc_report: dict[str, Any]) -> bool:
    repeat_spread_legs = list(qc_report.get("repeat_spread_legs", []))
    primary_repeat_spread_leg = str(qc_report.get("primary_repeat_spread_leg") or "")
    if primary_repeat_spread_leg not in {"complex", "apo"} or len(repeat_spread_legs) < 2:
        return False
    if primary_repeat_spread_leg not in repeat_spread_legs:
        return False

    secondary_legs = [leg for leg in repeat_spread_legs if leg != primary_repeat_spread_leg]
    if len(secondary_legs) != 1:
        return False

    legs = qc_report.get("legs", {})
    if not isinstance(legs, dict):
        return False
    primary_repeat_spread = _safe_float(
        legs.get(primary_repeat_spread_leg, {}).get("repeat_delta_kcal_mol_range")
    )
    secondary_repeat_spread = _safe_float(
        legs.get(secondary_legs[0], {}).get("repeat_delta_kcal_mol_range")
    )
    if primary_repeat_spread is None or secondary_repeat_spread is None:
        return False

    margin_threshold = max(_safe_float(qc_report.get("max_repeat_delta_kcal_mol")) or 0.0, 1.0)
    return (
        primary_repeat_spread >= secondary_repeat_spread + margin_threshold
        and primary_repeat_spread >= secondary_repeat_spread * 1.5
    )


def _can_target_primary_repeat_spread_leg(
    qc_report: dict[str, Any],
    rescue_reasons: list[str],
) -> bool:
    repeat_spread_legs = list(qc_report.get("repeat_spread_legs", []))
    primary_repeat_spread_leg = str(qc_report.get("primary_repeat_spread_leg") or "")
    if "repeat_spread" not in rescue_reasons:
        return False
    if primary_repeat_spread_leg not in {"complex", "apo"}:
        return False

    unexpected_reasons = set(rescue_reasons) - {"repeat_spread", "overlap"}
    if unexpected_reasons:
        return False

    if len(repeat_spread_legs) == 1:
        repeat_spread_targetable = True
    else:
        repeat_spread_targetable = _has_dominant_primary_repeat_spread_leg(qc_report)
    if not repeat_spread_targetable:
        return False

    if "overlap" in rescue_reasons:
        low_overlap_legs = _low_overlap_legs(qc_report)
        if set(low_overlap_legs) != {primary_repeat_spread_leg}:
            return False

    return True


def _scaled_int_protocol_value(protocol: dict[str, Any], key: str, scale: float) -> int | None:
    current = _safe_float(protocol.get(key))
    if current is None or current <= 0.0 or scale <= 1.0:
        return None
    scaled = max(int(current) + 1, int(ceil(current * scale)))
    return scaled


def _scaled_float_protocol_value(protocol: dict[str, Any], key: str, scale: float) -> float | None:
    current = _safe_float(protocol.get(key))
    if current is None or current <= 0.0 or scale <= 1.0:
        return None
    scaled = max(current * scale, current + 1e-6)
    return round(float(scaled), 6)


def _resolved_count_increment(increment: int, *, force: bool) -> int:
    resolved = max(int(increment), 0)
    if resolved > 0:
        return resolved
    return 1 if force else 0


def _recommended_release_npt_ps(protocol: dict[str, Any]) -> int:
    current_release = max(int(_safe_float(protocol.get("equilibration_release_npt_ps")) or 0), 0)
    nvt_ps = max(int(_safe_float(protocol.get("nvt_ps")) or 0), 0)
    npt_ps = max(int(_safe_float(protocol.get("npt_ps")) or 0), 0)
    if npt_ps <= 0:
        return current_release
    staged_default = max(nvt_ps, int(ceil(npt_ps * 0.5)), 20)
    return max(current_release, staged_default)


def _rescue_protocol_payload(
    protocol: dict[str, Any],
    *,
    reasons: list[str],
    repeat_increment: int,
    lambda_increment: int,
    production_scale: float,
    window_relax_em_scale: float,
    window_relax_md_scale: float,
    nvt_scale: float,
    npt_scale: float,
    force_repeat_increment: bool,
    force_lambda_increment: bool = False,
    preserve_counts: bool = False,
) -> tuple[dict[str, Any], dict[str, Any]]:
    rescued = dict(protocol)
    adjustments: dict[str, Any] = {}
    resolved_repeat_increment = _resolved_count_increment(repeat_increment, force=force_repeat_increment)
    resolved_lambda_increment = _resolved_count_increment(lambda_increment, force=force_lambda_increment)

    if production_scale > 1.0 and any(
        reason in {"repeat_spread", "bar_stderr", "overlap", "pass_qc_outlier"} for reason in reasons
    ):
        scaled_production = max(int(protocol["production_ps"]) + 1, int(ceil(float(protocol["production_ps"]) * production_scale)))
        if scaled_production != protocol["production_ps"]:
            rescued["production_ps"] = scaled_production
            adjustments["production_ps"] = scaled_production

    if any(reason in {"repeat_spread", "overlap", "pass_qc_outlier"} for reason in reasons):
        for key, scaled_value in (
            ("window_relax_em_steps", _scaled_int_protocol_value(protocol, "window_relax_em_steps", window_relax_em_scale)),
            ("window_relax_md_ps", _scaled_float_protocol_value(protocol, "window_relax_md_ps", window_relax_md_scale)),
            ("nvt_ps", _scaled_int_protocol_value(protocol, "nvt_ps", nvt_scale)),
            ("npt_ps", _scaled_int_protocol_value(protocol, "npt_ps", npt_scale)),
        ):
            if scaled_value is None:
                continue
            current = protocol.get(key)
            if current == scaled_value:
                continue
            rescued[key] = scaled_value
            adjustments[key] = scaled_value

    if not preserve_counts and resolved_repeat_increment > 0 and (
        force_repeat_increment or any(reason in {"repeat_spread", "pass_qc_outlier"} for reason in reasons)
    ):
        rescued_repeats = int(protocol["repeats"]) + resolved_repeat_increment
        if rescued_repeats != protocol["repeats"]:
            rescued["repeats"] = rescued_repeats
            adjustments["repeats"] = rescued_repeats

    if not preserve_counts and resolved_lambda_increment > 0 and (
        force_lambda_increment
        or any(reason in {"bar_stderr", "overlap", "pass_qc_outlier"} for reason in reasons)
    ):
        rescued_windows = int(protocol["lambda_windows"]) + resolved_lambda_increment
        if rescued_windows != protocol["lambda_windows"]:
            rescued["lambda_windows"] = rescued_windows
            adjustments["lambda_windows"] = rescued_windows

    if any(reason in {"repeat_spread", "overlap", "pass_qc_outlier"} for reason in reasons):
        schedule = str(rescued.get("equilibration_restraint_schedule") or "legacy_posres").strip() or "legacy_posres"
        if schedule != "staged_backbone_release":
            rescued["equilibration_restraint_schedule"] = "staged_backbone_release"
            adjustments["equilibration_restraint_schedule"] = "staged_backbone_release"

        release_npt_ps = _recommended_release_npt_ps(rescued)
        if release_npt_ps > max(int(_safe_float(rescued.get("equilibration_release_npt_ps")) or 0), 0):
            rescued["equilibration_release_npt_ps"] = release_npt_ps
            adjustments["equilibration_release_npt_ps"] = release_npt_ps

    return rescued, adjustments


def _protocol_effort_from_protocol(protocol: dict[str, Any]) -> tuple[float, float, float, float]:
    repeats = max(_safe_float(protocol.get("repeats")) or 0.0, 0.0)
    lambda_windows = max(_safe_float(protocol.get("lambda_windows")) or 0.0, 0.0)
    production_ps = max(_safe_float(protocol.get("production_ps")) or 0.0, 0.0)
    return repeats * lambda_windows * production_ps, repeats, lambda_windows, production_ps


def _split_plan_root_keys(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        items = [str(item).strip() for item in value]
    else:
        items = [item.strip() for item in str(value).split(",")]
    roots: list[str] = []
    for item in items:
        normalized = _normalize_plan_root_key(item)
        if not normalized or normalized in roots:
            continue
        roots.append(normalized)
    return roots


def _normalize_string_set(
    values: list[str] | None,
    *,
    transform: Callable[[str], str] | None = None,
) -> set[str]:
    normalized: set[str] = set()
    for value in values or []:
        item = str(value).strip()
        if not item:
            continue
        if transform is not None:
            item = transform(item)
        normalized.add(item)
    return normalized


def _resolve_rescue_protocol_controls(
    job_row: dict[str, Any],
    *,
    repeat_increment: int,
    lambda_increment: int,
    production_scale: float,
    window_relax_em_scale: float,
    window_relax_md_scale: float,
    nvt_scale: float,
    npt_scale: float,
    hotspot_complex_ids: set[str],
    hotspot_job_ids: set[str],
    hotspot_repeat_increment: int | None,
    hotspot_lambda_increment: int | None,
    hotspot_production_scale: float | None,
    hotspot_window_relax_em_scale: float | None,
    hotspot_window_relax_md_scale: float | None,
    hotspot_nvt_scale: float | None,
    hotspot_npt_scale: float | None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    controls = {
        "repeat_increment": int(repeat_increment),
        "lambda_increment": int(lambda_increment),
        "production_scale": float(production_scale),
        "window_relax_em_scale": float(window_relax_em_scale),
        "window_relax_md_scale": float(window_relax_md_scale),
        "nvt_scale": float(nvt_scale),
        "npt_scale": float(npt_scale),
    }

    job_id = str(job_row.get("job_id") or "").strip()
    normalized_job_id = job_id.lower()
    complex_id = str(job_row.get("complex_id") or "").strip()
    normalized_complex_id = complex_id.upper()
    match_sources: list[str] = []
    if normalized_complex_id and normalized_complex_id in hotspot_complex_ids:
        match_sources.append("complex_id")
    if normalized_job_id and normalized_job_id in hotspot_job_ids:
        match_sources.append("job_id")

    if match_sources:
        overrides = (
            ("repeat_increment", hotspot_repeat_increment),
            ("lambda_increment", hotspot_lambda_increment),
            ("production_scale", hotspot_production_scale),
            ("window_relax_em_scale", hotspot_window_relax_em_scale),
            ("window_relax_md_scale", hotspot_window_relax_md_scale),
            ("nvt_scale", hotspot_nvt_scale),
            ("npt_scale", hotspot_npt_scale),
        )
        for key, value in overrides:
            if value is None:
                continue
            controls[key] = value

    metadata = {
        "applied": bool(match_sources),
        "match_sources": match_sources,
        "matched_complex_id": complex_id if "complex_id" in match_sources else "",
        "matched_job_id": job_id if "job_id" in match_sources else "",
    }
    return controls, metadata


def _mutation_tokens_from_job_spec(job_spec: dict[str, Any]) -> str:
    return ";".join(
        f"{site['chain_id']}:{site['wt']}{site['resseq']}{site.get('icode', '')}{site['mut']}@{site['entity_side']}"
        for site in job_spec["mutation_group"]["sites"]
    )


def _candidate_source_matches_targeted_rescue(job_dir: Path, target_legs: list[str]) -> bool:
    if len(target_legs) != 1:
        return False
    rescue_path = job_dir / "config" / "rescue.json"
    if not rescue_path.is_file():
        return False

    rescue_payload = read_json(rescue_path)
    if str(rescue_payload.get("mode") or "") != "targeted_primary_repeat_spread_leg":
        return False
    configured_target_legs = [
        str(item).strip()
        for item in rescue_payload.get("target_legs", [])
        if str(item).strip()
    ]
    return configured_target_legs == target_legs


def _candidate_source_has_stage_progress(job_dir: Path) -> bool:
    stages_dir = job_dir / "stages"
    if not stages_dir.is_dir():
        return False
    return any(path.is_file() and path.suffix == ".json" for path in stages_dir.iterdir())


def _existing_rescue_batch_has_stage_progress(rescue_root: Path, batch_id: str) -> bool:
    jobs_dir = rescue_root / batch_id / "jobs"
    if not jobs_dir.is_dir():
        return False
    for job_dir in jobs_dir.iterdir():
        if job_dir.is_dir() and _candidate_source_has_stage_progress(job_dir):
            return True
    return False


def plan_ab_bind_rescues(
    plan_root: Path,
    *,
    extra_plan_roots: list[Path] | None = None,
    runs_root: Path | None = None,
    batch_prefix: str = "abbind-rescue",
    batch_ids: list[str] | None = None,
    complex_ids: list[str] | None = None,
    split_name: str | None = None,
    split_path: Path | None = None,
    limit_batches: int | None = None,
    job_ids: list[str] | None = None,
    repeat_increment: int = 1,
    lambda_increment: int = 4,
    production_scale: float = 2.0,
    window_relax_em_scale: float = 2.0,
    window_relax_md_scale: float = 2.0,
    nvt_scale: float = 2.0,
    npt_scale: float = 2.0,
    force_repeat_increment: bool = False,
    force_lambda_increment: bool = False,
    prefer_active_alternate_source: bool = False,
    require_active_alternate: bool = False,
    allow_pass_qc_outlier_rescue: bool = False,
    target_primary_repeat_spread_leg: bool = False,
    require_target_primary_repeat_spread_leg: bool = False,
    allow_targeted_leg_count_deepening: bool = False,
    hotspot_complex_ids: list[str] | None = None,
    hotspot_job_ids: list[str] | None = None,
    hotspot_repeat_increment: int | None = None,
    hotspot_lambda_increment: int | None = None,
    hotspot_production_scale: float | None = None,
    hotspot_window_relax_em_scale: float | None = None,
    hotspot_window_relax_md_scale: float | None = None,
    hotspot_nvt_scale: float | None = None,
    hotspot_npt_scale: float | None = None,
) -> dict[str, Any]:
    if require_target_primary_repeat_spread_leg and not target_primary_repeat_spread_leg:
        raise ValueError(
            "require_target_primary_repeat_spread_leg requires target_primary_repeat_spread_leg"
        )

    def _merge_unique_rows(
        rows: list[dict[str, Any]],
        *,
        key: str,
        sort_keys: tuple[str, ...],
    ) -> list[dict[str, Any]]:
        merged: dict[str, dict[str, Any]] = {}
        for row in rows:
            row_key = str(row.get(key, "")).strip()
            if not row_key:
                continue
            merged[row_key] = row
        return sorted(merged.values(), key=lambda item: tuple(str(item.get(field, "")) for field in sort_keys))

    def _merge_unique_strings(values: list[str]) -> list[str]:
        merged: list[str] = []
        seen: set[str] = set()
        for value in values:
            normalized = str(value).strip()
            if not normalized or normalized in seen:
                continue
            merged.append(normalized)
            seen.add(normalized)
        return merged

    normalized_plan_root = Path(plan_root).expanduser().resolve()
    normalized_hotspot_complex_ids = _normalize_string_set(
        hotspot_complex_ids,
        transform=lambda value: value.upper(),
    )
    normalized_hotspot_job_ids = _normalize_string_set(
        hotspot_job_ids,
        transform=lambda value: value.lower(),
    )
    normalized_extra_plan_roots: list[Path] = []
    for item in extra_plan_roots or []:
        resolved = Path(item).expanduser().resolve()
        if resolved == normalized_plan_root or resolved in normalized_extra_plan_roots:
            continue
        normalized_extra_plan_roots.append(resolved)
    report_bundle = report_ab_bind_plan(
        normalized_plan_root,
        extra_plan_roots=normalized_extra_plan_roots,
        batch_ids=batch_ids,
        complex_ids=complex_ids,
        split_name=split_name,
        split_path=split_path,
        limit_batches=limit_batches,
    )
    selected_report_dir = Path(report_bundle["reports_dir"])
    plan_jobs = read_csv_rows(selected_report_dir / "plan_jobs.csv")
    requested_job_ids = {item.strip() for item in job_ids or [] if item.strip()}
    if requested_job_ids:
        plan_jobs = [row for row in plan_jobs if row["job_id"] in requested_job_ids]

    source_plan_roots = [normalized_plan_root, *normalized_extra_plan_roots]
    batch_by_source_and_id: dict[tuple[str, str], dict[str, Any]] = {}
    batches_by_source_root: dict[str, list[dict[str, Any]]] = {}
    for source_root in source_plan_roots:
        plan_index = _load_ab_bind_plan_index(source_root)
        source_root_key = _normalize_plan_root_key(source_root)
        for item in plan_index.get("batches", []):
            batch_by_source_and_id[(source_root_key, str(item["batch_id"]))] = item
            batches_by_source_root.setdefault(source_root_key, []).append(item)

    rescue_root = ensure_dir(
        runs_root
        or (normalized_plan_root.parent / f"{normalized_plan_root.name}_rescues_{utc_now().replace(':', '').replace('-', '')}")
    )
    inputs_root = ensure_dir(rescue_root / "rescue_inputs")
    reports_dir = ensure_dir(rescue_root / "reports")
    existing_payload = (
        read_json(rescue_root / "plan_index.json") if (rescue_root / "plan_index.json").is_file() else {}
    )
    existing_batches = list(existing_payload.get("batches", [])) if isinstance(existing_payload, dict) else []
    existing_rescue_rows = (
        read_csv_rows(reports_dir / "rescue_candidates.csv") if (reports_dir / "rescue_candidates.csv").is_file() else []
    )

    rescue_rows: list[dict[str, Any]] = []
    planned_batches: list[dict[str, Any]] = []
    batch_prefix = slugify(batch_prefix)
    for job_row in plan_jobs:
        if job_row.get("ddg_ready") != "True":
            continue
        if require_active_alternate and not _as_bool(str(job_row.get("has_active_alternate_candidate", ""))):
            continue
        winner_source_plan_root = _normalize_plan_root_key(job_row.get("source_plan_root", str(normalized_plan_root)))
        winner_batch_meta = batch_by_source_and_id.get((winner_source_plan_root, str(job_row["batch_id"])))
        if winner_batch_meta is None:
            continue
        winner_job_dir = Path(winner_batch_meta["batch_dir"]) / "jobs" / job_row["job_id"]
        winner_qc_path = winner_job_dir / "results" / "qc_report.json"
        if not winner_qc_path.is_file():
            continue
        qc_report = read_json(winner_qc_path)
        rescue_reasons = _rescue_reason_codes(qc_report)
        pass_qc_outlier_abs_error = _safe_float(job_row.get("abs_ddg_error_kcal_mol"))
        if (
            not rescue_reasons
            and allow_pass_qc_outlier_rescue
            and bool(requested_job_ids)
            and job_row["job_id"] in requested_job_ids
            and str(job_row.get("qc_status") or "").strip() == "pass"
            and pass_qc_outlier_abs_error is not None
            and pass_qc_outlier_abs_error > 0.0
        ):
            rescue_reasons = ["pass_qc_outlier"]
        if not rescue_reasons:
            continue

        repeat_spread_legs = list(qc_report.get("repeat_spread_legs", []))
        primary_repeat_spread_leg = str(qc_report.get("primary_repeat_spread_leg") or "")
        source_complex_repeat_spread = _safe_float(qc_report.get("legs", {}).get("complex", {}).get("repeat_delta_kcal_mol_range"))
        source_apo_repeat_spread = _safe_float(qc_report.get("legs", {}).get("apo", {}).get("repeat_delta_kcal_mol_range"))
        targeted_repeat_spread_leg = bool(
            target_primary_repeat_spread_leg
            and _can_target_primary_repeat_spread_leg(qc_report, rescue_reasons)
        )
        if require_target_primary_repeat_spread_leg and not targeted_repeat_spread_leg:
            continue
        target_legs = [primary_repeat_spread_leg] if targeted_repeat_spread_leg else []
        inherit_source_legs = [leg for leg in ("complex", "apo") if leg not in target_legs]
        preserved_targeted_leg_counts = targeted_repeat_spread_leg and not allow_targeted_leg_count_deepening

        candidate_source_roots = [winner_source_plan_root]
        if prefer_active_alternate_source:
            candidate_source_roots.extend(_split_plan_root_keys(job_row.get("active_alternate_source_plan_roots", "")))
        if targeted_repeat_spread_leg:
            candidate_source_roots.extend(_split_plan_root_keys(source_plan_roots))
        candidate_source_roots = _split_plan_root_keys(candidate_source_roots)

        chosen_source_plan_root = ""
        chosen_batch_meta: dict[str, Any] | None = None
        chosen_job_dir: Path | None = None
        chosen_job_spec: dict[str, Any] | None = None
        chosen_priority: tuple[float, float, float, float, float] | None = None
        for candidate_index, candidate_source_plan_root in enumerate(candidate_source_roots):
            for batch_meta in batches_by_source_root.get(candidate_source_plan_root, []):
                job_dir = Path(batch_meta["batch_dir"]) / "jobs" / job_row["job_id"]
                spec_path = job_dir / "job_spec.json"
                if not spec_path.is_file():
                    continue
                job_spec = read_json(spec_path)
                effort, repeats, lambda_windows, production_ps = _protocol_effort_from_protocol(job_spec.get("protocol", {}))
                targeted_source_match = (
                    1.0
                    if targeted_repeat_spread_leg
                    and _candidate_source_matches_targeted_rescue(job_dir, target_legs)
                    and _candidate_source_has_stage_progress(job_dir)
                    else 0.0
                )
                candidate_priority = (
                    targeted_source_match,
                    effort,
                    repeats,
                    lambda_windows,
                    production_ps,
                    -float(candidate_index),
                )
                if chosen_priority is not None and candidate_priority <= chosen_priority:
                    continue
                chosen_source_plan_root = candidate_source_plan_root
                chosen_batch_meta = batch_meta
                chosen_job_dir = job_dir
                chosen_job_spec = job_spec
                chosen_priority = candidate_priority

        if chosen_batch_meta is None or chosen_job_dir is None or chosen_job_spec is None:
            chosen_source_plan_root = winner_source_plan_root
            chosen_batch_meta = winner_batch_meta
            chosen_job_dir = winner_job_dir
            spec_path = chosen_job_dir / "job_spec.json"
            if not spec_path.is_file():
                continue
            chosen_job_spec = read_json(spec_path)

        job_input_dir = ensure_dir(inputs_root / job_row["job_id"])
        system_yml = job_input_dir / "system.yml"
        protocol_yml = job_input_dir / "protocol.yml"
        mutations_csv = job_input_dir / "mutations.csv"
        write_yaml(system_yml, chosen_job_spec["system"])
        rescue_controls, rescue_override_metadata = _resolve_rescue_protocol_controls(
            job_row,
            repeat_increment=repeat_increment,
            lambda_increment=lambda_increment,
            production_scale=production_scale,
            window_relax_em_scale=window_relax_em_scale,
            window_relax_md_scale=window_relax_md_scale,
            nvt_scale=nvt_scale,
            npt_scale=npt_scale,
            hotspot_complex_ids=normalized_hotspot_complex_ids,
            hotspot_job_ids=normalized_hotspot_job_ids,
            hotspot_repeat_increment=hotspot_repeat_increment,
            hotspot_lambda_increment=hotspot_lambda_increment,
            hotspot_production_scale=hotspot_production_scale,
            hotspot_window_relax_em_scale=hotspot_window_relax_em_scale,
            hotspot_window_relax_md_scale=hotspot_window_relax_md_scale,
            hotspot_nvt_scale=hotspot_nvt_scale,
            hotspot_npt_scale=hotspot_npt_scale,
        )
        rescue_protocol, adjustments = _rescue_protocol_payload(
            chosen_job_spec["protocol"],
            reasons=rescue_reasons,
            repeat_increment=int(rescue_controls["repeat_increment"]),
            lambda_increment=int(rescue_controls["lambda_increment"]),
            production_scale=float(rescue_controls["production_scale"]),
            window_relax_em_scale=float(rescue_controls["window_relax_em_scale"]),
            window_relax_md_scale=float(rescue_controls["window_relax_md_scale"]),
            nvt_scale=float(rescue_controls["nvt_scale"]),
            npt_scale=float(rescue_controls["npt_scale"]),
            force_repeat_increment=force_repeat_increment,
            force_lambda_increment=force_lambda_increment,
            preserve_counts=preserved_targeted_leg_counts,
        )
        if not adjustments:
            continue

        rescue_batch_id = f"{batch_prefix}_{job_row['job_id']}"
        if _existing_rescue_batch_has_stage_progress(rescue_root, rescue_batch_id):
            continue

        write_yaml(protocol_yml, rescue_protocol)
        write_csv_rows(
            mutations_csv,
            [
                {
                    "mutation_group_id": chosen_job_spec["mutation_group"]["mutation_group_id"],
                    "entity_side": chosen_job_spec["mutation_group"]["entity_side"],
                    "mutation_tokens": _mutation_tokens_from_job_spec(chosen_job_spec),
                }
            ],
            ["mutation_group_id", "entity_side", "mutation_tokens"],
        )

        batch_plan = build_batch_plan(
            system_path=system_yml,
            mutations_path=mutations_csv,
            protocol_path=protocol_yml,
            batch_id=rescue_batch_id,
            runs_root=rescue_root,
        )
        if targeted_repeat_spread_leg:
            rescue_config = {
                "mode": "targeted_primary_repeat_spread_leg",
                "generated_at": utc_now(),
                "source_job_dir": str(chosen_job_dir.resolve()),
                "source_plan_root": chosen_source_plan_root,
                "source_batch_id": chosen_batch_meta.get("batch_id", job_row["batch_id"]),
                "source_job_id": job_row["job_id"],
                "source_winner_plan_root": winner_source_plan_root,
                "rescue_reasons": rescue_reasons,
                "target_legs": target_legs,
                "inherit_source_legs": inherit_source_legs,
                "allow_targeted_leg_count_deepening": allow_targeted_leg_count_deepening,
                "preserved_targeted_leg_counts": preserved_targeted_leg_counts,
            }
            for planned_job in batch_plan.jobs:
                write_json(Path(planned_job.workdir) / "config" / "rescue.json", rescue_config)
        planned_batches.append(
            {
                "complex_id": job_row["complex_id"],
                "batch_id": rescue_batch_id,
                "batch_dir": batch_plan.batch_dir,
                "system_yml": str(system_yml),
                "mutations_csv": str(mutations_csv),
                "job_count": len(batch_plan.jobs),
                "mutation_group_count": 1,
                "structure_source": chosen_batch_meta.get("structure_source", ""),
                "antibody_chains": chosen_batch_meta.get("antibody_chains", ""),
                "antigen_chains": chosen_batch_meta.get("antigen_chains", ""),
                "source_plan_root": chosen_source_plan_root,
                "source_batch_id": chosen_batch_meta.get("batch_id", job_row["batch_id"]),
                "source_job_id": job_row["job_id"],
                "rescue_reasons": ",".join(rescue_reasons),
                "repeat_spread_legs": ",".join(repeat_spread_legs),
                "primary_repeat_spread_leg": primary_repeat_spread_leg,
                "targeted_primary_repeat_spread_leg": targeted_repeat_spread_leg,
                "allow_targeted_leg_count_deepening": allow_targeted_leg_count_deepening,
                "preserved_targeted_leg_counts": preserved_targeted_leg_counts,
                "target_legs": ",".join(target_legs),
                "inherit_source_legs": ",".join(inherit_source_legs),
                "source_complex_repeat_spread_kcal_mol": source_complex_repeat_spread,
                "source_apo_repeat_spread_kcal_mol": source_apo_repeat_spread,
                "protocol_adjustments": adjustments,
                "hotspot_override_applied": rescue_override_metadata["applied"],
                "hotspot_override_match": ",".join(rescue_override_metadata["match_sources"]),
                "hotspot_override_complex_id": rescue_override_metadata["matched_complex_id"],
                "hotspot_override_job_id": rescue_override_metadata["matched_job_id"],
                "effective_repeat_increment": rescue_controls["repeat_increment"],
                "effective_lambda_increment": rescue_controls["lambda_increment"],
                "effective_production_scale": rescue_controls["production_scale"],
                "effective_window_relax_em_scale": rescue_controls["window_relax_em_scale"],
                "effective_window_relax_md_scale": rescue_controls["window_relax_md_scale"],
                "effective_nvt_scale": rescue_controls["nvt_scale"],
                "effective_npt_scale": rescue_controls["npt_scale"],
                "source_winner_plan_root": winner_source_plan_root,
                "source_winner_batch_id": job_row["batch_id"],
            }
        )
        rescue_rows.append(
            {
                "complex_id": job_row["complex_id"],
                "source_plan_root": chosen_source_plan_root,
                "source_winner_plan_root": winner_source_plan_root,
                "source_batch_id": chosen_batch_meta.get("batch_id", job_row["batch_id"]),
                "source_winner_batch_id": job_row["batch_id"],
                "source_job_id": job_row["job_id"],
                "rescue_batch_id": rescue_batch_id,
                "qc_status": job_row["qc_status"],
                "benchmark_qc_qualified": job_row["benchmark_qc_qualified"],
                "rescue_reasons": ",".join(rescue_reasons),
                "repeat_spread_legs": ",".join(repeat_spread_legs),
                "primary_repeat_spread_leg": primary_repeat_spread_leg,
                "targeted_primary_repeat_spread_leg": targeted_repeat_spread_leg,
                "allow_targeted_leg_count_deepening": allow_targeted_leg_count_deepening,
                "preserved_targeted_leg_counts": preserved_targeted_leg_counts,
                "target_legs": ",".join(target_legs),
                "inherit_source_legs": ",".join(inherit_source_legs),
                "source_complex_repeat_spread_kcal_mol": source_complex_repeat_spread,
                "source_apo_repeat_spread_kcal_mol": source_apo_repeat_spread,
                "hotspot_override_applied": rescue_override_metadata["applied"],
                "hotspot_override_match": ",".join(rescue_override_metadata["match_sources"]),
                "hotspot_override_complex_id": rescue_override_metadata["matched_complex_id"],
                "hotspot_override_job_id": rescue_override_metadata["matched_job_id"],
                "effective_repeat_increment": rescue_controls["repeat_increment"],
                "effective_lambda_increment": rescue_controls["lambda_increment"],
                "effective_production_scale": rescue_controls["production_scale"],
                "effective_window_relax_em_scale": rescue_controls["window_relax_em_scale"],
                "effective_window_relax_md_scale": rescue_controls["window_relax_md_scale"],
                "effective_nvt_scale": rescue_controls["nvt_scale"],
                "effective_npt_scale": rescue_controls["npt_scale"],
                "original_repeats": chosen_job_spec["protocol"]["repeats"],
                "rescued_repeats": rescue_protocol["repeats"],
                "original_lambda_windows": chosen_job_spec["protocol"]["lambda_windows"],
                "rescued_lambda_windows": rescue_protocol["lambda_windows"],
                "original_production_ps": chosen_job_spec["protocol"]["production_ps"],
                "rescued_production_ps": rescue_protocol["production_ps"],
                "original_window_relax_em_steps": chosen_job_spec["protocol"].get("window_relax_em_steps"),
                "rescued_window_relax_em_steps": rescue_protocol.get("window_relax_em_steps"),
                "original_window_relax_md_ps": chosen_job_spec["protocol"].get("window_relax_md_ps"),
                "rescued_window_relax_md_ps": rescue_protocol.get("window_relax_md_ps"),
                "original_nvt_ps": chosen_job_spec["protocol"].get("nvt_ps"),
                "rescued_nvt_ps": rescue_protocol.get("nvt_ps"),
                "original_npt_ps": chosen_job_spec["protocol"].get("npt_ps"),
                "rescued_npt_ps": rescue_protocol.get("npt_ps"),
                "original_equilibration_restraint_schedule": chosen_job_spec["protocol"].get("equilibration_restraint_schedule"),
                "rescued_equilibration_restraint_schedule": rescue_protocol.get("equilibration_restraint_schedule"),
                "original_equilibration_release_npt_ps": chosen_job_spec["protocol"].get("equilibration_release_npt_ps"),
                "rescued_equilibration_release_npt_ps": rescue_protocol.get("equilibration_release_npt_ps"),
            }
        )

    merged_batches = _merge_unique_rows(
        [*existing_batches, *planned_batches],
        key="batch_id",
        sort_keys=("complex_id", "batch_id"),
    )
    merged_rescue_rows = _merge_unique_rows(
        [*existing_rescue_rows, *rescue_rows],
        key="rescue_batch_id",
        sort_keys=("complex_id", "source_job_id", "rescue_batch_id"),
    )
    existing_source_plan_roots = (
        list(existing_payload.get("source_plan_roots", [])) if isinstance(existing_payload, dict) else []
    )
    merged_source_plan_roots = _merge_unique_strings(
        [*existing_source_plan_roots, *[str(item) for item in source_plan_roots]]
    )
    payload = {
        "generated_at": utc_now(),
        "source_plan_root": str(normalized_plan_root),
        "source_plan_roots": merged_source_plan_roots,
        "benchmark_root": str(report_bundle.get("benchmark_root", "")),
        "spec_name": str(report_bundle.get("spec_name", "")),
        "protocol_path": "",
        "plan_root": str(rescue_root),
        "split_name": split_name or str(existing_payload.get("split_name", "")),
        "split_path": str(split_path) if split_path is not None else str(existing_payload.get("split_path", "")),
        "planned_batch_count": len(merged_batches),
        "planned_complexes": sorted({item["complex_id"] for item in merged_batches}),
        "batches": merged_batches,
        "selected_job_count": len(merged_rescue_rows),
        "rescued_job_count": len(merged_rescue_rows),
        "repeat_increment": repeat_increment,
        "lambda_increment": lambda_increment,
        "production_scale": production_scale,
        "window_relax_em_scale": window_relax_em_scale,
        "window_relax_md_scale": window_relax_md_scale,
        "nvt_scale": nvt_scale,
        "npt_scale": npt_scale,
        "force_repeat_increment": force_repeat_increment,
        "force_lambda_increment": force_lambda_increment,
        "prefer_active_alternate_source": prefer_active_alternate_source,
        "require_active_alternate": require_active_alternate,
        "allow_pass_qc_outlier_rescue": allow_pass_qc_outlier_rescue,
        "target_primary_repeat_spread_leg": target_primary_repeat_spread_leg,
        "require_target_primary_repeat_spread_leg": require_target_primary_repeat_spread_leg,
        "allow_targeted_leg_count_deepening": allow_targeted_leg_count_deepening,
        "hotspot_complex_ids": sorted(normalized_hotspot_complex_ids),
        "hotspot_job_ids": sorted(normalized_hotspot_job_ids),
        "hotspot_repeat_increment": hotspot_repeat_increment,
        "hotspot_lambda_increment": hotspot_lambda_increment,
        "hotspot_production_scale": hotspot_production_scale,
        "hotspot_window_relax_em_scale": hotspot_window_relax_em_scale,
        "hotspot_window_relax_md_scale": hotspot_window_relax_md_scale,
        "hotspot_nvt_scale": hotspot_nvt_scale,
        "hotspot_npt_scale": hotspot_npt_scale,
        "reports_dir": str(reports_dir),
    }
    write_json(rescue_root / "plan_index.json", payload)
    write_yaml(rescue_root / "plan_index.yml", payload)
    write_json(reports_dir / "rescue_summary.json", payload)
    write_yaml(reports_dir / "rescue_summary.yml", payload)
    write_csv_rows(
        reports_dir / "rescue_candidates.csv",
        merged_rescue_rows,
        [
            "complex_id",
            "source_plan_root",
            "source_winner_plan_root",
            "source_batch_id",
            "source_winner_batch_id",
            "source_job_id",
            "rescue_batch_id",
            "qc_status",
            "benchmark_qc_qualified",
            "rescue_reasons",
            "repeat_spread_legs",
            "primary_repeat_spread_leg",
            "targeted_primary_repeat_spread_leg",
            "allow_targeted_leg_count_deepening",
            "preserved_targeted_leg_counts",
            "target_legs",
            "inherit_source_legs",
            "source_complex_repeat_spread_kcal_mol",
            "source_apo_repeat_spread_kcal_mol",
            "hotspot_override_applied",
            "hotspot_override_match",
            "hotspot_override_complex_id",
            "hotspot_override_job_id",
            "effective_repeat_increment",
            "effective_lambda_increment",
            "effective_production_scale",
            "effective_window_relax_em_scale",
            "effective_window_relax_md_scale",
            "effective_nvt_scale",
            "effective_npt_scale",
            "original_repeats",
            "rescued_repeats",
            "original_lambda_windows",
            "rescued_lambda_windows",
            "original_production_ps",
            "rescued_production_ps",
            "original_window_relax_em_steps",
            "rescued_window_relax_em_steps",
            "original_window_relax_md_ps",
            "rescued_window_relax_md_ps",
            "original_nvt_ps",
            "rescued_nvt_ps",
            "original_npt_ps",
            "rescued_npt_ps",
            "original_equilibration_restraint_schedule",
            "rescued_equilibration_restraint_schedule",
            "original_equilibration_release_npt_ps",
            "rescued_equilibration_release_npt_ps",
        ],
    )
    return payload


def run_ab_bind_plan(
    plan_root: Path,
    *,
    execute: bool,
    resume: bool = False,
    batch_ids: list[str] | None = None,
    complex_ids: list[str] | None = None,
    split_name: str | None = None,
    split_path: Path | None = None,
    job_ids: list[str] | None = None,
    limit_batches: int | None = None,
    limit_jobs: int | None = None,
    from_stage: str | None = None,
    to_stage: str | None = None,
    max_workers: int | None = None,
    gpu_devices: list[str] | None = None,
) -> dict[str, Any]:
    plan_index = _load_ab_bind_plan_index(plan_root)
    benchmark_root_value = str(plan_index.get("benchmark_root", "")).strip()
    benchmark_root = Path(benchmark_root_value) if benchmark_root_value else None
    requested_complex_ids = _resolve_ab_bind_complex_ids(
        available_complex_ids=[item["complex_id"] for item in plan_index.get("batches", [])],
        benchmark_root=benchmark_root,
        spec_name=str(plan_index.get("spec_name", "core_v1")),
        complex_ids=complex_ids,
        split_name=split_name,
        split_path=split_path,
    )
    selected_batches = _select_ab_bind_plan_batches(
        plan_index,
        batch_ids=batch_ids,
        complex_ids=requested_complex_ids,
        limit_batches=limit_batches,
    )
    if not selected_batches:
        raise ValueError("No planned AB-Bind batches matched the requested selection.")

    wanted_jobs = {item.strip() for item in job_ids or [] if item.strip()}
    remaining_jobs = limit_jobs
    execution_rows: list[dict[str, Any]] = []
    processed_batches: list[str] = []
    selected_job_records: list[dict[str, Any]] = []

    for batch in selected_batches:
        batch_dir = Path(batch["batch_dir"])
        jobs_dir = batch_dir / "jobs"
        selected_job_dirs = sorted(path for path in jobs_dir.iterdir() if path.is_dir())
        if wanted_jobs:
            selected_job_dirs = [job_dir for job_dir in selected_job_dirs if job_dir.name in wanted_jobs]
        if remaining_jobs is not None:
            if remaining_jobs <= 0:
                break
            selected_job_dirs = selected_job_dirs[:remaining_jobs]
        if not selected_job_dirs:
            continue

        processed_batches.append(batch["batch_id"])
        for job_dir in selected_job_dirs:
            selected_job_records.append(
                {
                    "complex_id": batch["complex_id"],
                    "batch_id": batch["batch_id"],
                    "batch_dir": batch_dir,
                    "job_dir": job_dir,
                }
            )
            if remaining_jobs is not None:
                remaining_jobs -= 1
                if remaining_jobs <= 0:
                    break

        if remaining_jobs is not None and remaining_jobs <= 0:
            break

    def execute_job(job_record: dict[str, Any], gpu_device: str | None = None) -> tuple[dict[str, Any], list[Any]]:
        job_dir = Path(job_record["job_dir"])
        environment = {"CUDA_VISIBLE_DEVICES": gpu_device} if gpu_device else None
        statuses = (
            resume_job(job_dir, execute=execute, environment=environment)
            if resume
            else run_job(
                job_dir,
                execute=execute,
                from_stage=from_stage,
                to_stage=to_stage,
                environment=environment,
            )
        )
        return job_record, statuses

    def append_execution_row(job_record: dict[str, Any], statuses: list[Any], gpu_device: str | None = None) -> None:
        last_status = statuses[-1] if statuses else None
        execution_rows.append(
            {
                "complex_id": job_record["complex_id"],
                "batch_id": job_record["batch_id"],
                "job_id": Path(job_record["job_dir"]).name,
                "mode": "resume" if resume else "run",
                "executed": execute,
                "from_stage": from_stage or "",
                "to_stage": to_stage or "",
                "gpu_device": gpu_device or "",
                "status_count": len(statuses),
                "final_stage": last_status.stage if last_status else "",
                "final_state": last_status.state if last_status else "already_completed",
                "final_message": last_status.message if last_status else "Job already completed.",
            }
        )

    requested_workers = max_workers or 1
    normalized_gpu_devices = _normalize_gpu_devices(gpu_devices)
    if execute and requested_workers > 1:
        worker_slots = list(normalized_gpu_devices[:requested_workers]) if normalized_gpu_devices else [None] * requested_workers
        if not worker_slots:
            worker_slots = [None]
        pending_job_records = deque(selected_job_records)
        with ThreadPoolExecutor(max_workers=len(worker_slots)) as executor:
            future_map = {}

            def submit_next(assigned_gpu: str | None) -> None:
                if not pending_job_records:
                    return
                record = pending_job_records.popleft()
                future = executor.submit(execute_job, record, assigned_gpu)
                future_map[future] = assigned_gpu

            for gpu_device in worker_slots:
                submit_next(gpu_device)

            while future_map:
                completed_futures, _ = wait(tuple(future_map), return_when=FIRST_COMPLETED)
                for future in completed_futures:
                    assigned_gpu = future_map.pop(future)
                    job_record, statuses = future.result()
                    append_execution_row(job_record, statuses, assigned_gpu)
                    write_batch_summary(Path(job_record["batch_dir"]))
                    report_ab_bind_plan(plan_root)
                    submit_next(assigned_gpu)
    else:
        job_count = len(selected_job_records)
        for index, job_record in enumerate(selected_job_records):
            assigned_gpu = None
            if normalized_gpu_devices:
                assigned_gpu = normalized_gpu_devices[index % len(normalized_gpu_devices)]
            _, statuses = execute_job(job_record, assigned_gpu)
            append_execution_row(job_record, statuses, assigned_gpu)
            write_batch_summary(Path(job_record["batch_dir"]))

    reports_dir = ensure_dir(plan_root / "reports")
    payload = {
        "generated_at": utc_now(),
        "plan_root": str(plan_root),
        "execute": execute,
        "resume": resume,
        "split_name": split_name or "",
        "split_path": str(split_path) if split_path is not None else "",
        "selected_batch_count": len(processed_batches),
        "selected_batches": processed_batches,
        "selected_job_count": len(execution_rows),
        "execution_rows": execution_rows,
    }
    if processed_batches:
        canonical_report = report_ab_bind_plan(plan_root)
        payload["canonical_reports_dir"] = canonical_report["reports_dir"]
    write_json(reports_dir / "run_summary.json", payload)
    write_yaml(reports_dir / "run_summary.yml", payload)
    if execution_rows:
        write_csv_rows(
            reports_dir / "run_summary.csv",
            execution_rows,
            [
                "complex_id",
                "batch_id",
                "job_id",
                "mode",
                "executed",
                "from_stage",
                "to_stage",
                "gpu_device",
                "status_count",
                "final_stage",
                "final_state",
                "final_message",
            ],
        )
    return payload


def report_ab_bind_plan(
    plan_root: Path,
    *,
    extra_plan_roots: list[Path] | None = None,
    batch_ids: list[str] | None = None,
    complex_ids: list[str] | None = None,
    split_name: str | None = None,
    split_path: Path | None = None,
    limit_batches: int | None = None,
) -> dict[str, Any]:
    plan_root = Path(plan_root).expanduser().resolve()
    plan_index = _load_ab_bind_plan_index(plan_root)
    all_batches = list(plan_index.get("batches", []))
    if not all_batches:
        raise ValueError("No planned AB-Bind batches were found under the requested plan root.")

    benchmark_root_value = str(plan_index.get("benchmark_root", "")).strip()
    benchmark_root = Path(benchmark_root_value) if benchmark_root_value else None
    requested_complex_ids = _resolve_ab_bind_complex_ids(
        available_complex_ids=[item["complex_id"] for item in all_batches],
        benchmark_root=benchmark_root,
        spec_name=str(plan_index.get("spec_name", "core_v1")),
        complex_ids=complex_ids,
        split_name=split_name,
        split_path=split_path,
    )

    has_selection = bool(batch_ids or requested_complex_ids or split_name or limit_batches is not None)
    selected_batches = all_batches
    selection_matches_full_scope = not has_selection
    if has_selection:
        selected_batches = _select_ab_bind_plan_batches(
            plan_index,
            batch_ids=batch_ids,
            complex_ids=requested_complex_ids,
            limit_batches=limit_batches,
        )
        if not selected_batches:
            raise ValueError("No planned AB-Bind batches matched the requested selection.")
        selection_batch_ids = [str(item.get("batch_id", "") or "") for item in selected_batches]
        all_batch_ids = [str(item.get("batch_id", "") or "") for item in all_batches]
        selection_matches_full_scope = (
            len(selection_batch_ids) == len(all_batch_ids)
            and set(selection_batch_ids) == set(all_batch_ids)
        )

    reports_dir = ensure_dir(plan_root / "reports")
    reference_rows = (
        _load_benchmark_reference_rows(benchmark_root, str(plan_index.get("spec_name", "")))
        if benchmark_root is not None
        else {}
    )
    (
        canonical_payload,
        canonical_batch_rows,
        canonical_job_rows,
        canonical_pair_rows,
        canonical_qc_qualified_pair_rows,
    ) = _collect_ab_bind_plan_report_rows(
        plan_index,
        all_batches,
        plan_root=plan_root,
        benchmark_root=benchmark_root,
        reference_rows=reference_rows,
    )
    canonical_payload = _write_ab_bind_plan_report_bundle(
        reports_dir,
        canonical_payload,
        canonical_batch_rows,
        canonical_job_rows,
        canonical_pair_rows,
        canonical_qc_qualified_pair_rows,
    )
    selection_payload: dict[str, Any] | None = None
    selection_batch_rows: list[dict[str, Any]] = []
    selection_job_rows: list[dict[str, Any]] = []
    selection_pair_rows: list[dict[str, Any]] = []
    selection_qc_qualified_pair_rows: list[dict[str, Any]] = []
    if has_selection:
        (
            selection_payload,
            selection_batch_rows,
            selection_job_rows,
            selection_pair_rows,
            selection_qc_qualified_pair_rows,
        ) = _collect_ab_bind_plan_report_rows(
            plan_index,
            selected_batches,
            plan_root=plan_root,
            benchmark_root=benchmark_root,
            reference_rows=reference_rows,
        )
    selection_bundle: dict[str, Any] | None = None
    selection_reports_dir: Path | None = None
    if has_selection:
        selection_payload = dict(selection_payload or {})
        selection_payload["selection"] = {
            "batch_ids": batch_ids or [],
            "complex_ids": requested_complex_ids or [],
            "split_name": split_name or "",
            "split_path": str(split_path) if split_path is not None else "",
            "limit_batches": limit_batches,
        }
        selection_payload["canonical_reports_dir"] = str(reports_dir)
        selection_reports_dir = reports_dir / "selections" / _selection_report_slug(
            batch_ids=batch_ids,
            complex_ids=requested_complex_ids,
            split_name=split_name,
            limit_batches=limit_batches,
        )
        selection_bundle = _write_ab_bind_plan_report_bundle(
            selection_reports_dir,
            selection_payload,
            selection_batch_rows,
            selection_job_rows,
            selection_pair_rows,
            selection_qc_qualified_pair_rows,
        )
    if not has_selection:
        if not extra_plan_roots:
            return canonical_payload
    elif not extra_plan_roots:
        return selection_bundle or canonical_payload

    merged_roots = [plan_root]
    for root in extra_plan_roots or []:
        resolved = Path(root).expanduser().resolve()
        if resolved not in merged_roots:
            merged_roots.append(resolved)
    selected_complex_ids = [item["complex_id"] for item in (selected_batches if has_selection else all_batches)]
    selected_complex_ids = list(dict.fromkeys(selected_complex_ids))
    batch_metadata_by_complex = {item["complex_id"]: item for item in (selected_batches if has_selection else all_batches)}
    merged_candidates: list[dict[str, Any]] = []
    for root in merged_roots:
        root_index = _load_ab_bind_plan_index(root)
        root_spec_name = str(root_index.get("spec_name", ""))
        if root_spec_name and str(plan_index.get("spec_name", "")) and root_spec_name != str(plan_index.get("spec_name", "")):
            raise ValueError("Merged AB-Bind reports require matching spec_name across all plan roots.")
        root_benchmark_root_value = str(root_index.get("benchmark_root", "")).strip()
        root_benchmark_root = Path(root_benchmark_root_value) if root_benchmark_root_value else benchmark_root
        root_reference_rows = (
            _load_benchmark_reference_rows(root_benchmark_root, root_spec_name or str(plan_index.get("spec_name", "")))
            if root_benchmark_root is not None
            else reference_rows
        )
        root_batches = _select_ab_bind_plan_batches(
            root_index,
            complex_ids=selected_complex_ids,
        )
        if not root_batches:
            continue
        (
            _root_payload,
            _root_batch_rows,
            root_job_rows,
            _root_pair_rows,
            _root_qc_qualified_pair_rows,
        ) = _collect_ab_bind_plan_report_rows(
            root_index,
            root_batches,
            plan_root=root,
            benchmark_root=root_benchmark_root,
            reference_rows=root_reference_rows,
        )
        merged_candidates.extend(root_job_rows)

    root_priorities = {str(root): index for index, root in enumerate(merged_roots)}
    merged_job_rows = _merge_ab_bind_job_rows(merged_candidates, root_priorities=root_priorities)
    (
        merged_payload,
        merged_batch_rows,
        merged_job_rows,
        merged_pair_rows,
        merged_qc_qualified_pair_rows,
    ) = _build_ab_bind_report_rows_from_jobs(
        merged_job_rows,
        plan_root=plan_root,
        benchmark_root=benchmark_root,
        spec_name=str(plan_index.get("spec_name", "")),
        batch_metadata_by_complex=batch_metadata_by_complex,
        source_plan_roots=merged_roots,
    )
    merged_payload["merge_strategy"] = {
        "preferred_order": [str(root) for root in merged_roots],
        "winner_priority": [
            "benchmark_qc_qualified",
            "ddg_ready",
            "qc_status",
            "analyzable",
            "live_stage_state",
            "latest_stage",
            "latest_stage_state",
            "protocol_sampling_effort",
            "lower_repeat_spread_kcal_mol",
            "lower_ddg_bar_stderr_kcal_mol",
            "resumable",
            "plan_root_order",
        ],
    }
    if has_selection:
        merged_payload["selection"] = {
            "batch_ids": batch_ids or [],
            "complex_ids": requested_complex_ids or [],
            "split_name": split_name or "",
            "split_path": str(split_path) if split_path is not None else "",
            "limit_batches": limit_batches,
        }
    merged_payload["canonical_reports_dir"] = str(reports_dir)
    merged_reports_dir = reports_dir / "merged"
    if has_selection:
        merged_reports_dir = merged_reports_dir / "selections" / _selection_report_slug(
            batch_ids=batch_ids,
            complex_ids=requested_complex_ids,
            split_name=split_name,
            limit_batches=limit_batches,
        )
    merged_bundle = _write_ab_bind_plan_report_bundle(
        merged_reports_dir,
        merged_payload,
        merged_batch_rows,
        merged_job_rows,
        merged_pair_rows,
        merged_qc_qualified_pair_rows,
    )
    if has_selection and selection_matches_full_scope:
        canonical_merged_payload = dict(merged_payload)
        canonical_merged_payload.pop("selection", None)
        _write_ab_bind_plan_report_bundle(
            reports_dir / "merged",
            canonical_merged_payload,
            merged_batch_rows,
            merged_job_rows,
            merged_pair_rows,
            merged_qc_qualified_pair_rows,
        )
    return merged_bundle


def _job_entity_side(job_id: str) -> str:
    normalized = str(job_id).strip().lower()
    if "-antibody-" in normalized:
        return "antibody"
    if "-antigen-" in normalized:
        return "antigen"
    return "unknown"


def _entity_side_indicator(entity_side: str) -> float:
    return 1.0 if str(entity_side).strip().lower() == "antibody" else 0.0


def _fit_linear_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_pairs: list[tuple[float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        clean_pairs.append((predicted, experimental))
    if not clean_pairs:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    xs = [predicted for predicted, _experimental in clean_pairs]
    ys = [experimental for _predicted, experimental in clean_pairs]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    variance_x = sum((value - mean_x) ** 2 for value in xs)
    if variance_x == 0:
        slope = 0.0
        intercept = mean_y
    else:
        covariance = sum((value_x - mean_x) * (value_y - mean_y) for value_x, value_y in zip(xs, ys, strict=True))
        slope = covariance / variance_x
        intercept = mean_y - slope * mean_x
    fitted = [intercept + slope * value for value in xs]
    return {
        "label": label,
        "family": "linear",
        "pair_count": len(clean_pairs),
        "intercept": intercept,
        "slope": slope,
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental in clean_pairs
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, experimental in zip(fitted, ys, strict=True)
            ]
        ),
    }


def _solve_dense_linear_system(matrix: list[list[float]], vector: list[float]) -> list[float] | None:
    order = len(vector)
    if order == 0 or len(matrix) != order:
        return None

    working_matrix = [list(row) for row in matrix]
    working_vector = list(vector)
    for pivot_index in range(order):
        pivot_row = max(range(pivot_index, order), key=lambda row_index: abs(working_matrix[row_index][pivot_index]))
        pivot_value = working_matrix[pivot_row][pivot_index]
        if abs(pivot_value) <= 1e-12:
            return None
        if pivot_row != pivot_index:
            working_matrix[pivot_index], working_matrix[pivot_row] = (
                working_matrix[pivot_row],
                working_matrix[pivot_index],
            )
            working_vector[pivot_index], working_vector[pivot_row] = (
                working_vector[pivot_row],
                working_vector[pivot_index],
            )
        scale = 1.0 / working_matrix[pivot_index][pivot_index]
        for column_index in range(pivot_index, order):
            working_matrix[pivot_index][column_index] *= scale
        working_vector[pivot_index] *= scale
        for row_index in range(order):
            if row_index == pivot_index:
                continue
            factor = working_matrix[row_index][pivot_index]
            if factor == 0.0:
                continue
            for column_index in range(pivot_index, order):
                working_matrix[row_index][column_index] -= factor * working_matrix[pivot_index][column_index]
            working_vector[row_index] -= factor * working_vector[pivot_index]
    return working_vector


def _solve_linear_least_squares(feature_rows: list[list[float]], targets: list[float]) -> list[float] | None:
    if not feature_rows or not targets or len(feature_rows) != len(targets):
        return None
    order = len(feature_rows[0])
    if order == 0 or any(len(row) != order for row in feature_rows):
        return None

    matrix = [[0.0] * order for _ in range(order)]
    vector = [0.0] * order
    for features, target in zip(feature_rows, targets, strict=True):
        for row_index, feature_row in enumerate(features):
            vector[row_index] += feature_row * target
            for column_index, feature_column in enumerate(features):
                matrix[row_index][column_index] += feature_row * feature_column
    return _solve_dense_linear_system(matrix, vector)


def _fit_quadratic_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_pairs: list[tuple[float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        clean_pairs.append((predicted, experimental))
    if not clean_pairs:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    xs = [predicted for predicted, _experimental in clean_pairs]
    ys = [experimental for _predicted, experimental in clean_pairs]
    coefficients = _solve_linear_least_squares([[1.0, value, value * value] for value in xs], ys)
    if coefficients is None:
        linear_payload = _fit_linear_calibration(rows, label=label)
        intercept = float(linear_payload["intercept"])
        linear_coefficient = float(linear_payload["slope"])
        quadratic_coefficient = 0.0
        fitted = [intercept + linear_coefficient * value for value in xs]
    else:
        intercept, linear_coefficient, quadratic_coefficient = coefficients
        fitted = [
            intercept + linear_coefficient * value + quadratic_coefficient * value * value
            for value in xs
        ]

    return {
        "label": label,
        "family": "quadratic",
        "pair_count": len(clean_pairs),
        "intercept": intercept,
        "linear_coefficient": linear_coefficient,
        "quadratic_coefficient": quadratic_coefficient,
        "coefficients": [intercept, linear_coefficient, quadratic_coefficient],
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental in clean_pairs
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, experimental in zip(fitted, ys, strict=True)
            ]
        ),
    }


def _calibration_bar_stderr(row: dict[str, Any]) -> float:
    stderr = _safe_float(row.get("ddg_bar_stderr_kcal_mol"))
    if stderr is None:
        stderr = _safe_float(row.get("max_bar_stderr_kcal_mol"))
    return max(stderr or 0.0, 0.0)


_EXPDECAY_INVSTDERR_DEFAULT_RATE = 0.6
_HILL_INVSTDERR_DEFAULT_EXPONENT = 2.7
_HILL_INVSTDERR_DEFAULT_MIDPOINT = 1.5
_HILL_SIDE_INVSTDERR_DEFAULT_EXPONENT = 3.0
_HILL_SIDE_INVSTDERR_DEFAULT_MIDPOINT = 1.65


def _inverse_stderr_feature(stderr: float) -> float:
    return 1.0 / (1.0 + max(stderr, 0.0))


def _expdecay_logabs_feature(predicted_ddg: float, *, rate: float = _EXPDECAY_INVSTDERR_DEFAULT_RATE) -> float:
    return 1.0 - exp(-rate * log1p(abs(predicted_ddg)))


def _hill_logabs_feature(
    predicted_ddg: float,
    *,
    exponent: float = _HILL_INVSTDERR_DEFAULT_EXPONENT,
    midpoint: float = _HILL_INVSTDERR_DEFAULT_MIDPOINT,
) -> float:
    logabs = log1p(abs(predicted_ddg))
    numerator = logabs**exponent
    denominator = numerator + max(midpoint, 1e-12)
    return numerator / denominator if denominator > 0.0 else 0.0


def _fit_stderr_quadratic_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_rows: list[tuple[float, float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        clean_rows.append((predicted, experimental, _calibration_bar_stderr(row)))
    if not clean_rows:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    feature_rows = [
        [1.0, predicted, predicted * predicted, stderr, log1p(stderr)]
        for predicted, _experimental, stderr in clean_rows
    ]
    ys = [experimental for _predicted, experimental, _stderr in clean_rows]
    coefficients = _solve_linear_least_squares(feature_rows, ys)
    if coefficients is None:
        quadratic_payload = _fit_quadratic_calibration(rows, label=label)
        quadratic_payload["fallback_from_family"] = "stderr_quadratic"
        return quadratic_payload

    intercept, linear_coefficient, quadratic_coefficient, stderr_coefficient, log_stderr_coefficient = coefficients
    fitted = [
        intercept
        + linear_coefficient * predicted
        + quadratic_coefficient * predicted * predicted
        + stderr_coefficient * stderr
        + log_stderr_coefficient * log1p(stderr)
        for predicted, _experimental, stderr in clean_rows
    ]

    return {
        "label": label,
        "family": "stderr_quadratic",
        "pair_count": len(clean_rows),
        "intercept": intercept,
        "linear_coefficient": linear_coefficient,
        "quadratic_coefficient": quadratic_coefficient,
        "stderr_coefficient": stderr_coefficient,
        "log_stderr_coefficient": log_stderr_coefficient,
        "coefficients": [
            intercept,
            linear_coefficient,
            quadratic_coefficient,
            stderr_coefficient,
            log_stderr_coefficient,
        ],
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental, _stderr in clean_rows
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, (_predicted, experimental, _stderr) in zip(fitted, clean_rows, strict=True)
            ]
        ),
    }


def _fit_logabs_stderr_quadratic_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_rows: list[tuple[float, float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        clean_rows.append((predicted, experimental, _calibration_bar_stderr(row)))
    if not clean_rows:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    feature_rows = []
    for predicted, _experimental, stderr in clean_rows:
        logabs = log1p(abs(predicted))
        feature_rows.append([1.0, logabs, logabs * logabs, stderr, log1p(stderr)])
    ys = [experimental for _predicted, experimental, _stderr in clean_rows]
    coefficients = _solve_linear_least_squares(feature_rows, ys)
    if coefficients is None:
        quadratic_payload = _fit_quadratic_calibration(rows, label=label)
        quadratic_payload["fallback_from_family"] = "logabs_stderr_quadratic"
        return quadratic_payload

    (
        intercept,
        logabs_linear_coefficient,
        logabs_quadratic_coefficient,
        stderr_coefficient,
        log_stderr_coefficient,
    ) = coefficients
    fitted = []
    for predicted, _experimental, stderr in clean_rows:
        logabs = log1p(abs(predicted))
        fitted.append(
            intercept
            + logabs_linear_coefficient * logabs
            + logabs_quadratic_coefficient * logabs * logabs
            + stderr_coefficient * stderr
            + log_stderr_coefficient * log1p(stderr)
        )

    return {
        "label": label,
        "family": "logabs_stderr_quadratic",
        "pair_count": len(clean_rows),
        "intercept": intercept,
        "logabs_linear_coefficient": logabs_linear_coefficient,
        "logabs_quadratic_coefficient": logabs_quadratic_coefficient,
        "stderr_coefficient": stderr_coefficient,
        "log_stderr_coefficient": log_stderr_coefficient,
        "coefficients": [
            intercept,
            logabs_linear_coefficient,
            logabs_quadratic_coefficient,
            stderr_coefficient,
            log_stderr_coefficient,
        ],
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental, _stderr in clean_rows
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, (_predicted, experimental, _stderr) in zip(fitted, clean_rows, strict=True)
            ]
        ),
    }


def _fit_expdecay_invstderr_quadratic_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_rows: list[tuple[float, float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        clean_rows.append((predicted, experimental, _calibration_bar_stderr(row)))
    if not clean_rows:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    expdecay_rate = _EXPDECAY_INVSTDERR_DEFAULT_RATE
    feature_rows = []
    for predicted, _experimental, stderr in clean_rows:
        expdecay = _expdecay_logabs_feature(predicted, rate=expdecay_rate)
        inverse_stderr = _inverse_stderr_feature(stderr)
        feature_rows.append([1.0, expdecay, expdecay * expdecay, inverse_stderr])
    ys = [experimental for _predicted, experimental, _stderr in clean_rows]
    coefficients = _solve_linear_least_squares(feature_rows, ys)
    if coefficients is None:
        logabs_payload = _fit_logabs_stderr_quadratic_calibration(rows, label=label)
        logabs_payload["fallback_from_family"] = "expdecay_invstderr_quadratic"
        return logabs_payload

    intercept, expdecay_linear_coefficient, expdecay_quadratic_coefficient, inv_stderr_coefficient = coefficients
    fitted = []
    for predicted, _experimental, stderr in clean_rows:
        expdecay = _expdecay_logabs_feature(predicted, rate=expdecay_rate)
        inverse_stderr = _inverse_stderr_feature(stderr)
        fitted.append(
            intercept
            + expdecay_linear_coefficient * expdecay
            + expdecay_quadratic_coefficient * expdecay * expdecay
            + inv_stderr_coefficient * inverse_stderr
        )

    return {
        "label": label,
        "family": "expdecay_invstderr_quadratic",
        "pair_count": len(clean_rows),
        "intercept": intercept,
        "expdecay_rate": expdecay_rate,
        "expdecay_linear_coefficient": expdecay_linear_coefficient,
        "expdecay_quadratic_coefficient": expdecay_quadratic_coefficient,
        "inv_stderr_coefficient": inv_stderr_coefficient,
        "coefficients": [
            intercept,
            expdecay_linear_coefficient,
            expdecay_quadratic_coefficient,
            inv_stderr_coefficient,
        ],
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental, _stderr in clean_rows
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, (_predicted, experimental, _stderr) in zip(fitted, clean_rows, strict=True)
            ]
        ),
    }


def _fit_hill_invstderr_quadratic_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_rows: list[tuple[float, float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        clean_rows.append((predicted, experimental, _calibration_bar_stderr(row)))
    if not clean_rows:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    hill_exponent = _HILL_INVSTDERR_DEFAULT_EXPONENT
    hill_midpoint = _HILL_INVSTDERR_DEFAULT_MIDPOINT
    feature_rows = []
    for predicted, _experimental, stderr in clean_rows:
        hill = _hill_logabs_feature(predicted, exponent=hill_exponent, midpoint=hill_midpoint)
        inverse_stderr = _inverse_stderr_feature(stderr)
        feature_rows.append([1.0, hill, hill * hill, inverse_stderr])
    ys = [experimental for _predicted, experimental, _stderr in clean_rows]
    coefficients = _solve_linear_least_squares(feature_rows, ys)
    if coefficients is None:
        expdecay_payload = _fit_expdecay_invstderr_quadratic_calibration(rows, label=label)
        expdecay_payload["fallback_from_family"] = "hill_invstderr_quadratic"
        return expdecay_payload

    intercept, hill_linear_coefficient, hill_quadratic_coefficient, inv_stderr_coefficient = coefficients
    fitted = []
    for predicted, _experimental, stderr in clean_rows:
        hill = _hill_logabs_feature(predicted, exponent=hill_exponent, midpoint=hill_midpoint)
        inverse_stderr = _inverse_stderr_feature(stderr)
        fitted.append(
            intercept
            + hill_linear_coefficient * hill
            + hill_quadratic_coefficient * hill * hill
            + inv_stderr_coefficient * inverse_stderr
        )

    return {
        "label": label,
        "family": "hill_invstderr_quadratic",
        "pair_count": len(clean_rows),
        "intercept": intercept,
        "hill_exponent": hill_exponent,
        "hill_midpoint": hill_midpoint,
        "hill_linear_coefficient": hill_linear_coefficient,
        "hill_quadratic_coefficient": hill_quadratic_coefficient,
        "inv_stderr_coefficient": inv_stderr_coefficient,
        "coefficients": [
            intercept,
            hill_linear_coefficient,
            hill_quadratic_coefficient,
            inv_stderr_coefficient,
        ],
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental, _stderr in clean_rows
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, (_predicted, experimental, _stderr) in zip(fitted, clean_rows, strict=True)
            ]
        ),
    }


def _fit_hill_side_invstderr_quadratic_calibration(
    rows: list[dict[str, Any]],
    *,
    label: str,
) -> dict[str, Any]:
    clean_rows: list[tuple[float, float, float, float]] = []
    for row in rows:
        predicted = _safe_float(row.get("predicted_ddg_kcal_mol"))
        experimental = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if predicted is None or experimental is None:
            continue
        entity_side = _job_entity_side(str(row.get("job_id", "")))
        clean_rows.append(
            (
                predicted,
                experimental,
                _calibration_bar_stderr(row),
                _entity_side_indicator(entity_side),
            )
        )
    if not clean_rows:
        raise ValueError(f"Cannot fit calibration model '{label}' without paired benchmark rows.")

    hill_exponent = _HILL_SIDE_INVSTDERR_DEFAULT_EXPONENT
    hill_midpoint = _HILL_SIDE_INVSTDERR_DEFAULT_MIDPOINT
    feature_rows = []
    for predicted, _experimental, stderr, side_indicator in clean_rows:
        hill = _hill_logabs_feature(predicted, exponent=hill_exponent, midpoint=hill_midpoint)
        inverse_stderr = _inverse_stderr_feature(stderr)
        feature_rows.append([1.0, hill, hill * hill, inverse_stderr, side_indicator])
    ys = [experimental for _predicted, experimental, _stderr, _side_indicator in clean_rows]
    coefficients = _solve_linear_least_squares(feature_rows, ys)
    if coefficients is None:
        hill_payload = _fit_hill_invstderr_quadratic_calibration(rows, label=label)
        hill_payload["fallback_from_family"] = "hill_side_invstderr_quadratic"
        return hill_payload

    (
        intercept,
        hill_linear_coefficient,
        hill_quadratic_coefficient,
        inv_stderr_coefficient,
        antibody_side_coefficient,
    ) = coefficients
    fitted = []
    for predicted, _experimental, stderr, side_indicator in clean_rows:
        hill = _hill_logabs_feature(predicted, exponent=hill_exponent, midpoint=hill_midpoint)
        inverse_stderr = _inverse_stderr_feature(stderr)
        fitted.append(
            intercept
            + hill_linear_coefficient * hill
            + hill_quadratic_coefficient * hill * hill
            + inv_stderr_coefficient * inverse_stderr
            + antibody_side_coefficient * side_indicator
        )

    return {
        "label": label,
        "family": "hill_side_invstderr_quadratic",
        "pair_count": len(clean_rows),
        "intercept": intercept,
        "hill_exponent": hill_exponent,
        "hill_midpoint": hill_midpoint,
        "hill_linear_coefficient": hill_linear_coefficient,
        "hill_quadratic_coefficient": hill_quadratic_coefficient,
        "inv_stderr_coefficient": inv_stderr_coefficient,
        "antibody_side_coefficient": antibody_side_coefficient,
        "coefficients": [
            intercept,
            hill_linear_coefficient,
            hill_quadratic_coefficient,
            inv_stderr_coefficient,
            antibody_side_coefficient,
        ],
        "raw_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": predicted,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": predicted - experimental,
                }
                for predicted, experimental, _stderr, _side_indicator in clean_rows
            ]
        ),
        "fitted_metrics": _benchmark_metrics_from_pairs(
            [
                {
                    "predicted_ddg_kcal_mol": fitted_value,
                    "experimental_ddg_kcal_mol": experimental,
                    "ddg_error_kcal_mol": fitted_value - experimental,
                }
                for fitted_value, (_predicted, experimental, _stderr, _side_indicator) in zip(
                    fitted,
                    clean_rows,
                    strict=True,
                )
            ]
        ),
    }


def _fit_ab_bind_calibration_model(
    fit_rows: list[dict[str, Any]],
    *,
    model: str,
) -> dict[str, Any]:
    if model not in {
        "linear",
        "side_linear",
        "quadratic",
        "stderr_quadratic",
        "logabs_stderr_quadratic",
        "expdecay_invstderr_quadratic",
        "hill_invstderr_quadratic",
        "hill_side_invstderr_quadratic",
    }:
        raise ValueError(f"Unsupported calibration model: {model}")

    if model == "quadratic":
        fit_group = _fit_quadratic_calibration
    elif model == "stderr_quadratic":
        fit_group = _fit_stderr_quadratic_calibration
    elif model == "logabs_stderr_quadratic":
        fit_group = _fit_logabs_stderr_quadratic_calibration
    elif model == "expdecay_invstderr_quadratic":
        fit_group = _fit_expdecay_invstderr_quadratic_calibration
    elif model == "hill_invstderr_quadratic":
        fit_group = _fit_hill_invstderr_quadratic_calibration
    elif model == "hill_side_invstderr_quadratic":
        fit_group = _fit_hill_side_invstderr_quadratic_calibration
    else:
        fit_group = _fit_linear_calibration
    groups = {"global": fit_group(fit_rows, label="global")}
    if model == "side_linear":
        for entity_side in ("antibody", "antigen"):
            group_rows = [row for row in fit_rows if _job_entity_side(str(row.get("job_id", ""))) == entity_side]
            if not group_rows:
                continue
            groups[entity_side] = _fit_linear_calibration(group_rows, label=entity_side)
    return {
        "model": model,
        "groups": groups,
    }


def _calibration_group_for_job(job_id: str, *, model: str) -> str:
    if model == "side_linear":
        entity_side = _job_entity_side(job_id)
        if entity_side in {"antibody", "antigen"}:
            return entity_side
    return "global"


def _apply_calibration_group(raw_ddg: float, fit_payload: dict[str, Any]) -> float:
    family = str(fit_payload.get("family") or "linear").strip().lower()
    intercept = _safe_float(fit_payload.get("intercept")) or 0.0
    if family == "hill_side_invstderr_quadratic":
        hill_exponent = _safe_float(fit_payload.get("hill_exponent")) or _HILL_SIDE_INVSTDERR_DEFAULT_EXPONENT
        hill_midpoint = _safe_float(fit_payload.get("hill_midpoint")) or _HILL_SIDE_INVSTDERR_DEFAULT_MIDPOINT
        hill_linear_coefficient = _safe_float(fit_payload.get("hill_linear_coefficient")) or 0.0
        hill_quadratic_coefficient = _safe_float(fit_payload.get("hill_quadratic_coefficient")) or 0.0
        inv_stderr_coefficient = _safe_float(fit_payload.get("inv_stderr_coefficient")) or 0.0
        antibody_side_coefficient = _safe_float(fit_payload.get("antibody_side_coefficient")) or 0.0
        stderr = max(_safe_float(fit_payload.get("_input_ddg_bar_stderr_kcal_mol")) or 0.0, 0.0)
        entity_side = str(fit_payload.get("_input_job_entity_side") or "").strip().lower()
        side_indicator = _entity_side_indicator(entity_side)
        hill = _hill_logabs_feature(raw_ddg, exponent=hill_exponent, midpoint=hill_midpoint)
        inverse_stderr = _inverse_stderr_feature(stderr)
        return (
            intercept
            + hill_linear_coefficient * hill
            + hill_quadratic_coefficient * hill * hill
            + inv_stderr_coefficient * inverse_stderr
            + antibody_side_coefficient * side_indicator
        )
    if family == "hill_invstderr_quadratic":
        hill_exponent = _safe_float(fit_payload.get("hill_exponent")) or _HILL_INVSTDERR_DEFAULT_EXPONENT
        hill_midpoint = _safe_float(fit_payload.get("hill_midpoint")) or _HILL_INVSTDERR_DEFAULT_MIDPOINT
        hill_linear_coefficient = _safe_float(fit_payload.get("hill_linear_coefficient")) or 0.0
        hill_quadratic_coefficient = _safe_float(fit_payload.get("hill_quadratic_coefficient")) or 0.0
        inv_stderr_coefficient = _safe_float(fit_payload.get("inv_stderr_coefficient")) or 0.0
        stderr = max(_safe_float(fit_payload.get("_input_ddg_bar_stderr_kcal_mol")) or 0.0, 0.0)
        hill = _hill_logabs_feature(raw_ddg, exponent=hill_exponent, midpoint=hill_midpoint)
        inverse_stderr = _inverse_stderr_feature(stderr)
        return (
            intercept
            + hill_linear_coefficient * hill
            + hill_quadratic_coefficient * hill * hill
            + inv_stderr_coefficient * inverse_stderr
        )
    if family == "expdecay_invstderr_quadratic":
        expdecay_rate = _safe_float(fit_payload.get("expdecay_rate")) or _EXPDECAY_INVSTDERR_DEFAULT_RATE
        expdecay_linear_coefficient = _safe_float(fit_payload.get("expdecay_linear_coefficient")) or 0.0
        expdecay_quadratic_coefficient = _safe_float(fit_payload.get("expdecay_quadratic_coefficient")) or 0.0
        inv_stderr_coefficient = _safe_float(fit_payload.get("inv_stderr_coefficient")) or 0.0
        stderr = max(_safe_float(fit_payload.get("_input_ddg_bar_stderr_kcal_mol")) or 0.0, 0.0)
        expdecay = _expdecay_logabs_feature(raw_ddg, rate=expdecay_rate)
        inverse_stderr = _inverse_stderr_feature(stderr)
        return (
            intercept
            + expdecay_linear_coefficient * expdecay
            + expdecay_quadratic_coefficient * expdecay * expdecay
            + inv_stderr_coefficient * inverse_stderr
        )
    if family == "logabs_stderr_quadratic":
        logabs_linear_coefficient = _safe_float(fit_payload.get("logabs_linear_coefficient")) or 0.0
        logabs_quadratic_coefficient = _safe_float(fit_payload.get("logabs_quadratic_coefficient")) or 0.0
        stderr_coefficient = _safe_float(fit_payload.get("stderr_coefficient")) or 0.0
        log_stderr_coefficient = _safe_float(fit_payload.get("log_stderr_coefficient")) or 0.0
        stderr = max(_safe_float(fit_payload.get("_input_ddg_bar_stderr_kcal_mol")) or 0.0, 0.0)
        logabs = log1p(abs(raw_ddg))
        return (
            intercept
            + logabs_linear_coefficient * logabs
            + logabs_quadratic_coefficient * logabs * logabs
            + stderr_coefficient * stderr
            + log_stderr_coefficient * log1p(stderr)
        )
    if family == "stderr_quadratic":
        linear_coefficient = _safe_float(fit_payload.get("linear_coefficient"))
        if linear_coefficient is None:
            linear_coefficient = _safe_float(fit_payload.get("slope")) or 0.0
        quadratic_coefficient = _safe_float(fit_payload.get("quadratic_coefficient")) or 0.0
        stderr_coefficient = _safe_float(fit_payload.get("stderr_coefficient")) or 0.0
        log_stderr_coefficient = _safe_float(fit_payload.get("log_stderr_coefficient")) or 0.0
        stderr = max(_safe_float(fit_payload.get("_input_ddg_bar_stderr_kcal_mol")) or 0.0, 0.0)
        return (
            intercept
            + linear_coefficient * raw_ddg
            + quadratic_coefficient * raw_ddg * raw_ddg
            + stderr_coefficient * stderr
            + log_stderr_coefficient * log1p(stderr)
        )
    if family == "quadratic":
        linear_coefficient = _safe_float(fit_payload.get("linear_coefficient"))
        if linear_coefficient is None:
            linear_coefficient = _safe_float(fit_payload.get("slope")) or 0.0
        quadratic_coefficient = _safe_float(fit_payload.get("quadratic_coefficient")) or 0.0
        return intercept + linear_coefficient * raw_ddg + quadratic_coefficient * raw_ddg * raw_ddg
    slope = _safe_float(fit_payload.get("slope"))
    if slope is None:
        slope = _safe_float(fit_payload.get("linear_coefficient")) or 0.0
    return intercept + slope * raw_ddg


def _value_range_summary(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"min": None, "max": None}
    return {"min": min(values), "max": max(values)}


def _experimental_effect_bin(value: float, *, threshold_kcal_mol: float) -> str:
    if value <= -threshold_kcal_mol:
        return "negative"
    if value >= threshold_kcal_mol:
        return "positive"
    return "near_zero"


def _calibration_pair_coverage(
    rows: list[dict[str, Any]],
    *,
    predicted_field: str = "predicted_ddg_kcal_mol",
    effect_threshold_kcal_mol: float = 1.0,
) -> dict[str, Any]:
    by_side = {"antibody": 0, "antigen": 0, "unknown": 0}
    experimental_effect_bin_counts = {"negative": 0, "near_zero": 0, "positive": 0}
    experimental_values: list[float] = []
    predicted_values: list[float] = []

    for row in rows:
        job_id = str(row.get("job_id", "")).strip()
        entity_side = _job_entity_side(job_id)
        by_side[entity_side if entity_side in {"antibody", "antigen"} else "unknown"] += 1

        experimental_ddg = _safe_float(row.get("experimental_ddg_kcal_mol"))
        if experimental_ddg is not None:
            experimental_values.append(experimental_ddg)
            experimental_effect_bin_counts[
                _experimental_effect_bin(experimental_ddg, threshold_kcal_mol=effect_threshold_kcal_mol)
            ] += 1

        predicted_ddg = _safe_float(row.get(predicted_field))
        if predicted_ddg is not None:
            predicted_values.append(predicted_ddg)

    return {
        "pair_count": len(rows),
        "by_side": by_side,
        "effect_threshold_kcal_mol": effect_threshold_kcal_mol,
        "experimental_effect_bin_counts": experimental_effect_bin_counts,
        "experimental_range_kcal_mol": _value_range_summary(experimental_values),
        "predicted_range_kcal_mol": _value_range_summary(predicted_values),
    }


def _normalize_fit_split_names(
    *,
    fit_split_name: str,
    fit_split_names: list[str] | None = None,
) -> list[str]:
    candidates = fit_split_names if fit_split_names is not None else [fit_split_name]
    normalized: list[str] = []
    for item in candidates:
        value = str(item).strip()
        if not value or value in normalized:
            continue
        normalized.append(value)
    if not normalized:
        raise ValueError("At least one fit split name is required for calibration.")
    return normalized


def _selection_reports_dir_from_bundle(bundle: dict[str, Any]) -> Path:
    reports_dir = str(bundle.get("reports_dir", "")).strip()
    if not reports_dir:
        raise ValueError("Benchmark report bundle did not contain a reports_dir.")
    return Path(reports_dir)


def _calibration_report_slug(
    *,
    fit_split_names: list[str],
    predict_split_name: str,
    model: str,
    fit_qc_qualified_only: bool,
) -> str:
    parts = [
        "fit-" + slugify("--".join(fit_split_names)),
        f"predict-{slugify(predict_split_name)}",
        f"model-{slugify(model)}",
    ]
    if fit_qc_qualified_only:
        parts.append("fit-qc-qualified")
    return slugify("--".join(parts))


def calibrate_ab_bind_plan(
    plan_root: Path,
    *,
    extra_plan_roots: list[Path] | None = None,
    fit_split_name: str = "calibration",
    fit_split_names: list[str] | None = None,
    predict_split_name: str = "validation",
    split_path: Path | None = None,
    model: str = "hill_side_invstderr_quadratic",
    fit_qc_qualified_only: bool = False,
    fit_reports_dirs: list[Path] | None = None,
    predict_reports_dir: Path | None = None,
) -> dict[str, Any]:
    normalized_plan_root = Path(plan_root).expanduser().resolve()
    normalized_extra_plan_roots: list[Path] = []
    for item in extra_plan_roots or []:
        resolved = Path(item).expanduser().resolve()
        if resolved == normalized_plan_root or resolved in normalized_extra_plan_roots:
            continue
        normalized_extra_plan_roots.append(resolved)

    normalized_fit_split_names = _normalize_fit_split_names(
        fit_split_name=fit_split_name,
        fit_split_names=fit_split_names,
    )
    normalized_fit_reports_dirs = [Path(path).expanduser().resolve() for path in fit_reports_dirs or []]
    normalized_predict_reports_dir = (
        Path(predict_reports_dir).expanduser().resolve() if predict_reports_dir is not None else None
    )
    if normalized_fit_reports_dirs or normalized_predict_reports_dir is not None:
        if not normalized_fit_reports_dirs or normalized_predict_reports_dir is None:
            raise ValueError("fit_reports_dirs and predict_reports_dir must be provided together.")
        if len(normalized_fit_reports_dirs) != len(normalized_fit_split_names):
            raise ValueError(
                "fit_reports_dirs must match the number of requested fit split names "
                f"({len(normalized_fit_split_names)})."
            )
    else:
        fit_bundles = [
            report_ab_bind_plan(
                normalized_plan_root,
                extra_plan_roots=normalized_extra_plan_roots,
                split_name=split_name,
                split_path=split_path,
            )
            for split_name in normalized_fit_split_names
        ]
        predict_bundle = report_ab_bind_plan(
            normalized_plan_root,
            extra_plan_roots=normalized_extra_plan_roots,
            split_name=predict_split_name,
            split_path=split_path,
        )
        normalized_fit_reports_dirs = [_selection_reports_dir_from_bundle(bundle) for bundle in fit_bundles]
        normalized_predict_reports_dir = _selection_reports_dir_from_bundle(predict_bundle)

    fit_reports_dirs = normalized_fit_reports_dirs
    predict_reports_dir = normalized_predict_reports_dir
    fit_pairs_filename = "benchmark_pairs_qc_qualified.csv" if fit_qc_qualified_only else "benchmark_pairs.csv"
    fit_pairs: list[dict[str, str]] = []
    for fit_reports_dir in fit_reports_dirs:
        fit_pairs.extend(read_csv_rows(fit_reports_dir / fit_pairs_filename))
    if not fit_pairs:
        joined_paths = ", ".join(str(path / fit_pairs_filename) for path in fit_reports_dirs)
        raise ValueError(f"No calibration benchmark pairs were available in {joined_paths}.")

    calibration_model = _fit_ab_bind_calibration_model(fit_pairs, model=model)
    predict_jobs = read_csv_rows(predict_reports_dir / "plan_jobs.csv")
    calibrated_jobs: list[dict[str, Any]] = []
    raw_pairs: list[dict[str, Any]] = []
    calibrated_pairs: list[dict[str, Any]] = []

    for row in predict_jobs:
        raw_ddg = _safe_float(row.get("ddg_kcal_mol"))
        experimental_ddg = _safe_float(row.get("experimental_ddg_kcal_mol"))
        group = _calibration_group_for_job(str(row.get("job_id", "")), model=model)
        fit_payload = calibration_model["groups"].get(group) or calibration_model["groups"]["global"]
        calibration_bar_stderr = _calibration_bar_stderr(row)
        fit_payload_with_inputs = dict(fit_payload)
        fit_payload_with_inputs["_input_ddg_bar_stderr_kcal_mol"] = calibration_bar_stderr
        fit_payload_with_inputs["_input_job_entity_side"] = _job_entity_side(str(row.get("job_id", "")))
        calibrated_ddg = None if raw_ddg is None else _apply_calibration_group(raw_ddg, fit_payload_with_inputs)
        calibrated_error = None if calibrated_ddg is None or experimental_ddg is None else calibrated_ddg - experimental_ddg
        raw_error = None if raw_ddg is None or experimental_ddg is None else raw_ddg - experimental_ddg

        calibrated_row = dict(row)
        calibrated_row["calibration_group"] = group
        calibrated_row["calibration_model"] = calibration_model["model"]
        calibrated_row["calibration_family"] = fit_payload.get("family", "linear")
        calibrated_row["calibration_intercept"] = fit_payload["intercept"]
        calibrated_row["calibration_slope"] = fit_payload.get("slope", fit_payload.get("linear_coefficient"))
        calibrated_row["calibration_linear_coefficient"] = fit_payload.get(
            "linear_coefficient",
            fit_payload.get("slope"),
        )
        calibrated_row["calibration_quadratic_coefficient"] = fit_payload.get("quadratic_coefficient")
        calibrated_row["calibration_logabs_linear_coefficient"] = fit_payload.get("logabs_linear_coefficient")
        calibrated_row["calibration_logabs_quadratic_coefficient"] = fit_payload.get("logabs_quadratic_coefficient")
        calibrated_row["calibration_stderr_coefficient"] = fit_payload.get("stderr_coefficient")
        calibrated_row["calibration_log_stderr_coefficient"] = fit_payload.get("log_stderr_coefficient")
        calibrated_row["calibration_expdecay_rate"] = fit_payload.get("expdecay_rate")
        calibrated_row["calibration_expdecay_linear_coefficient"] = fit_payload.get(
            "expdecay_linear_coefficient"
        )
        calibrated_row["calibration_expdecay_quadratic_coefficient"] = fit_payload.get(
            "expdecay_quadratic_coefficient"
        )
        calibrated_row["calibration_hill_exponent"] = fit_payload.get("hill_exponent")
        calibrated_row["calibration_hill_midpoint"] = fit_payload.get("hill_midpoint")
        calibrated_row["calibration_hill_linear_coefficient"] = fit_payload.get("hill_linear_coefficient")
        calibrated_row["calibration_hill_quadratic_coefficient"] = fit_payload.get(
            "hill_quadratic_coefficient"
        )
        calibrated_row["calibration_inv_stderr_coefficient"] = fit_payload.get("inv_stderr_coefficient")
        calibrated_row["calibration_antibody_side_coefficient"] = fit_payload.get("antibody_side_coefficient")
        calibrated_row["calibration_input_ddg_bar_stderr_kcal_mol"] = calibration_bar_stderr
        calibrated_row["calibrated_ddg_kcal_mol"] = calibrated_ddg
        calibrated_row["calibrated_ddg_error_kcal_mol"] = calibrated_error
        calibrated_jobs.append(calibrated_row)

        if raw_ddg is not None and experimental_ddg is not None and str(row.get("ddg_ready", "")).strip() == "True":
            raw_pairs.append(
                {
                    "complex_id": row.get("complex_id", ""),
                    "batch_id": row.get("batch_id", ""),
                    "job_id": row.get("job_id", ""),
                    "mutation_group_id": row.get("mutation_group_id", ""),
                    "source_plan_root": row.get("source_plan_root", ""),
                    "predicted_ddg_kcal_mol": raw_ddg,
                    "experimental_ddg_kcal_mol": experimental_ddg,
                    "ddg_error_kcal_mol": raw_error,
                }
            )
            if calibrated_ddg is not None:
                calibrated_pairs.append(
                    {
                        "complex_id": row.get("complex_id", ""),
                        "batch_id": row.get("batch_id", ""),
                        "job_id": row.get("job_id", ""),
                        "mutation_group_id": row.get("mutation_group_id", ""),
                        "source_plan_root": row.get("source_plan_root", ""),
                        "calibration_group": group,
                        "raw_ddg_kcal_mol": raw_ddg,
                        "predicted_ddg_kcal_mol": calibrated_ddg,
                        "experimental_ddg_kcal_mol": experimental_ddg,
                        "ddg_error_kcal_mol": calibrated_error,
                    }
                )

    raw_metrics = _benchmark_metrics_from_pairs(raw_pairs) if raw_pairs else {}
    calibrated_metrics = _benchmark_metrics_from_pairs(calibrated_pairs) if calibrated_pairs else {}
    predict_raw_target_bundle = _benchmark_target_metrics_bundle(raw_pairs)
    predict_calibrated_target_bundle = _benchmark_target_metrics_bundle(calibrated_pairs)
    predict_raw_outlier_trim_bundle = _benchmark_outlier_trim_bundle(raw_pairs)
    predict_calibrated_outlier_trim_bundle = _benchmark_outlier_trim_bundle(calibrated_pairs)
    predict_raw_target_filtered_outlier_trim_bundle = _benchmark_outlier_trim_bundle(
        predict_raw_target_bundle["filtered_pair_rows"]
    )
    predict_calibrated_target_filtered_outlier_trim_bundle = _benchmark_outlier_trim_bundle(
        predict_calibrated_target_bundle["filtered_pair_rows"]
    )
    plan_index = _load_ab_bind_plan_index(normalized_plan_root)
    reports_dir = ensure_dir(
        normalized_plan_root
        / "reports"
        / "calibrations"
        / _calibration_report_slug(
            fit_split_names=normalized_fit_split_names,
            predict_split_name=predict_split_name,
            model=model,
            fit_qc_qualified_only=fit_qc_qualified_only,
        )
    )
    payload = {
        "generated_at": utc_now(),
        "plan_root": str(normalized_plan_root),
        "source_plan_roots": [str(path) for path in [normalized_plan_root, *normalized_extra_plan_roots]],
        "benchmark_root": str(plan_index.get("benchmark_root", "")),
        "spec_name": str(plan_index.get("spec_name", "")),
        "fit_split_name": normalized_fit_split_names[0],
        "fit_split_names": normalized_fit_split_names,
        "predict_split_name": predict_split_name,
        "split_path": str(split_path) if split_path is not None else "",
        "fit_reports_dir": str(fit_reports_dirs[0]),
        "fit_reports_dirs": [str(path) for path in fit_reports_dirs],
        "predict_reports_dir": str(predict_reports_dir),
        "reports_dir": str(reports_dir),
        "fit_qc_qualified_only": fit_qc_qualified_only,
        "model": calibration_model,
        "fit_pair_count": len(fit_pairs),
        "predict_job_count": len(calibrated_jobs),
        "predict_pair_count": len(calibrated_pairs),
        "fit_coverage": _calibration_pair_coverage(fit_pairs),
        "predict_raw_coverage": _calibration_pair_coverage(raw_pairs),
        "predict_calibrated_coverage": _calibration_pair_coverage(calibrated_pairs),
        "predict_target_exclusion_policy": {
            "target_field": predict_raw_target_bundle["target_field"],
            "systematically_poor_abs_error_threshold_kcal_mol": predict_raw_target_bundle[
                "systematically_poor_abs_error_threshold_kcal_mol"
            ],
            "systematically_poor_min_pair_count": predict_raw_target_bundle["systematically_poor_min_pair_count"],
            "systematically_poor_max_pearson_r": predict_raw_target_bundle["systematically_poor_max_pearson_r"],
            "systematically_poor_max_sign_accuracy": predict_raw_target_bundle[
                "systematically_poor_max_sign_accuracy"
            ],
            "systematically_poor_min_leave_one_out_pearson_gain": predict_raw_target_bundle[
                "systematically_poor_min_leave_one_out_pearson_gain"
            ],
        },
        "predict_outlier_trim_policy": {
            "target_field": predict_raw_outlier_trim_bundle["target_field"],
            "outlier_trim_method": predict_raw_outlier_trim_bundle["outlier_trim_method"],
        },
        "raw_metrics": raw_metrics,
        "predict_raw_target_metrics": predict_raw_target_bundle["target_metrics"],
        "predict_raw_target_outlier_trim_metrics": predict_raw_outlier_trim_bundle["target_metrics"],
        "predict_raw_target_excluded_complex_ids": predict_raw_target_bundle["excluded_target_ids"],
        "predict_raw_outlier_trimmed_metrics": predict_raw_outlier_trim_bundle["trimmed_metrics"],
        "predict_raw_target_filtered_metrics": predict_raw_target_bundle["filtered_metrics"],
        "predict_raw_target_filtered_outlier_trimmed_metrics": (
            predict_raw_target_filtered_outlier_trim_bundle["trimmed_metrics"]
        ),
        "calibrated_metrics": calibrated_metrics,
        "predict_calibrated_target_metrics": predict_calibrated_target_bundle["target_metrics"],
        "predict_calibrated_target_outlier_trim_metrics": predict_calibrated_outlier_trim_bundle[
            "target_metrics"
        ],
        "predict_calibrated_target_excluded_complex_ids": predict_calibrated_target_bundle["excluded_target_ids"],
        "predict_calibrated_outlier_trimmed_metrics": predict_calibrated_outlier_trim_bundle[
            "trimmed_metrics"
        ],
        "predict_calibrated_target_filtered_metrics": predict_calibrated_target_bundle["filtered_metrics"],
        "predict_calibrated_target_filtered_outlier_trimmed_metrics": (
            predict_calibrated_target_filtered_outlier_trim_bundle["trimmed_metrics"]
        ),
    }
    write_json(reports_dir / "summary.json", payload)
    write_yaml(reports_dir / "summary.yml", payload)
    write_json(reports_dir / "model.json", calibration_model)
    write_yaml(reports_dir / "model.yml", calibration_model)
    write_csv_rows(reports_dir / "fit_pairs.csv", fit_pairs, list(fit_pairs[0].keys()))
    target_metric_fields = [
        "complex_id",
        "paired_job_count",
        "pearson_r",
        "spearman_rho",
        "rmse_kcal_mol",
        "mae_kcal_mol",
        "sign_accuracy",
        "auc_strong_effect",
        "min_abs_error_kcal_mol",
        "max_abs_error_kcal_mol",
        "mean_abs_error_kcal_mol",
        "sign_mismatch_count",
        "all_pairs_sign_mismatched",
        "all_pairs_above_abs_error_threshold",
        "systematically_poor_abs_error_threshold_kcal_mol",
        "systematically_poor_min_pair_count",
        "systematically_poor_max_pearson_r",
        "systematically_poor_max_sign_accuracy",
        "systematically_poor_min_leave_one_out_pearson_gain",
        "systematically_poor_abs_error_target",
        "systematically_poor_correlation_target",
        "systematically_poor_target_reason",
        "target_exclusion_iteration",
        "overall_pearson_r_at_exclusion",
        "leave_one_out_pearson_r_at_exclusion",
        "leave_one_out_pearson_gain_at_exclusion",
        "systematically_poor_target",
        "excluded_from_target_filtered_metrics",
        "leave_one_out_paired_job_count",
        "leave_one_out_pearson_r",
        "leave_one_out_spearman_rho",
        "leave_one_out_rmse_kcal_mol",
        "leave_one_out_mae_kcal_mol",
        "leave_one_out_sign_accuracy",
        "leave_one_out_auc_strong_effect",
    ]
    outlier_trim_target_metric_fields = [
        "complex_id",
        "outlier_trim_method",
        "original_paired_job_count",
        "trimmed_paired_job_count",
        "removed_pair_count",
        "removed_fraction",
        "q1_abs_error_kcal_mol",
        "q3_abs_error_kcal_mol",
        "iqr_abs_error_kcal_mol",
        "threshold_abs_error_kcal_mol",
        "removed_job_ids",
        "removed_job_ids_text",
        "trimmed_pearson_r",
        "trimmed_spearman_rho",
        "trimmed_rmse_kcal_mol",
        "trimmed_mae_kcal_mol",
        "trimmed_sign_accuracy",
        "trimmed_auc_strong_effect",
    ]
    if calibrated_jobs:
        write_csv_rows(reports_dir / "predict_jobs_calibrated.csv", calibrated_jobs, list(calibrated_jobs[0].keys()))
    if calibrated_pairs:
        write_csv_rows(reports_dir / "predict_pairs_calibrated.csv", calibrated_pairs, list(calibrated_pairs[0].keys()))
    write_csv_rows(
        reports_dir / "predict_target_metrics_raw.csv",
        predict_raw_target_bundle["target_metrics"],
        target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "predict_target_metrics_calibrated.csv",
        predict_calibrated_target_bundle["target_metrics"],
        target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "predict_target_outlier_trim_metrics_raw.csv",
        predict_raw_outlier_trim_bundle["target_metrics"],
        outlier_trim_target_metric_fields,
    )
    write_csv_rows(
        reports_dir / "predict_target_outlier_trim_metrics_calibrated.csv",
        predict_calibrated_outlier_trim_bundle["target_metrics"],
        outlier_trim_target_metric_fields,
    )
    write_json(
        reports_dir / "predict_metrics_raw_outlier_trimmed.json",
        predict_raw_outlier_trim_bundle["trimmed_metrics"],
    )
    write_yaml(
        reports_dir / "predict_metrics_raw_outlier_trimmed.yml",
        predict_raw_outlier_trim_bundle["trimmed_metrics"],
    )
    write_json(
        reports_dir / "predict_metrics_calibrated_outlier_trimmed.json",
        predict_calibrated_outlier_trim_bundle["trimmed_metrics"],
    )
    write_yaml(
        reports_dir / "predict_metrics_calibrated_outlier_trimmed.yml",
        predict_calibrated_outlier_trim_bundle["trimmed_metrics"],
    )
    write_json(
        reports_dir / "predict_metrics_raw_target_filtered_outlier_trimmed.json",
        predict_raw_target_filtered_outlier_trim_bundle["trimmed_metrics"],
    )
    write_yaml(
        reports_dir / "predict_metrics_raw_target_filtered_outlier_trimmed.yml",
        predict_raw_target_filtered_outlier_trim_bundle["trimmed_metrics"],
    )
    write_json(
        reports_dir / "predict_metrics_calibrated_target_filtered_outlier_trimmed.json",
        predict_calibrated_target_filtered_outlier_trim_bundle["trimmed_metrics"],
    )
    write_yaml(
        reports_dir / "predict_metrics_calibrated_target_filtered_outlier_trimmed.yml",
        predict_calibrated_target_filtered_outlier_trim_bundle["trimmed_metrics"],
    )
    write_csv_rows(
        reports_dir / "predict_pairs_raw_outlier_trimmed.csv",
        predict_raw_outlier_trim_bundle["trimmed_pair_rows"],
        list(raw_pairs[0].keys()) if raw_pairs else [
            "complex_id",
            "batch_id",
            "job_id",
            "mutation_group_id",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
        ],
    )
    write_csv_rows(
        reports_dir / "predict_pairs_calibrated_outlier_trimmed.csv",
        predict_calibrated_outlier_trim_bundle["trimmed_pair_rows"],
        list(calibrated_pairs[0].keys()) if calibrated_pairs else [
            "complex_id",
            "batch_id",
            "job_id",
            "mutation_group_id",
            "calibration_group",
            "raw_ddg_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
        ],
    )
    write_csv_rows(
        reports_dir / "predict_pairs_raw_target_filtered.csv",
        predict_raw_target_bundle["filtered_pair_rows"],
        list(raw_pairs[0].keys()) if raw_pairs else [
            "complex_id",
            "batch_id",
            "job_id",
            "mutation_group_id",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
        ],
    )
    write_csv_rows(
        reports_dir / "predict_pairs_raw_target_filtered_outlier_trimmed.csv",
        predict_raw_target_filtered_outlier_trim_bundle["trimmed_pair_rows"],
        list(raw_pairs[0].keys()) if raw_pairs else [
            "complex_id",
            "batch_id",
            "job_id",
            "mutation_group_id",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
        ],
    )
    write_csv_rows(
        reports_dir / "predict_pairs_calibrated_target_filtered.csv",
        predict_calibrated_target_bundle["filtered_pair_rows"],
        list(calibrated_pairs[0].keys()) if calibrated_pairs else [
            "complex_id",
            "batch_id",
            "job_id",
            "mutation_group_id",
            "calibration_group",
            "raw_ddg_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
        ],
    )
    write_csv_rows(
        reports_dir / "predict_pairs_calibrated_target_filtered_outlier_trimmed.csv",
        predict_calibrated_target_filtered_outlier_trim_bundle["trimmed_pair_rows"],
        list(calibrated_pairs[0].keys()) if calibrated_pairs else [
            "complex_id",
            "batch_id",
            "job_id",
            "mutation_group_id",
            "calibration_group",
            "raw_ddg_kcal_mol",
            "predicted_ddg_kcal_mol",
            "experimental_ddg_kcal_mol",
            "ddg_error_kcal_mol",
        ],
    )
    return payload
