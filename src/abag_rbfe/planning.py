"""System/protocol/mutation loading and batch planning."""

from __future__ import annotations

from pathlib import Path
from dataclasses import replace
import re

from abag_pmx.mutations import (
    build_job_id,
    build_mutation_group,
    load_mutation_groups_from_csv,
)
from abag_rbfe.constants import preset_copy
from abag_rbfe.io_utils import ensure_dir, read_yaml, utc_now, write_csv_rows, write_json, write_yaml
from abag_rbfe.models import BatchPlan, JobSpec, ProtocolConfig, SystemConfig, dataclass_to_dict
from abag_rbfe.paths import ProjectPaths


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "batch"


def load_system_config(path: Path) -> SystemConfig:
    data = read_yaml(path)
    for key in ("system_name", "input_structure", "antibody_chains", "antigen_chains"):
        if key not in data:
            raise ValueError(f"Missing required system field: {key}")
    notes = tuple(str(item) for item in data.get("notes", []))
    return SystemConfig(
        system_name=str(data["system_name"]),
        input_structure=str(Path(str(data["input_structure"])).expanduser().resolve()),
        structure_source=str(data.get("structure_source", "experimental")),
        antibody_chains=tuple(str(item) for item in data["antibody_chains"]),
        antigen_chains=tuple(str(item) for item in data["antigen_chains"]),
        notes=notes,
    )


def load_protocol_config(path: Path, preset_name: str | None = None) -> ProtocolConfig:
    raw = read_yaml(path)
    return hydrate_protocol_config(raw, preset_name=preset_name)


def hydrate_protocol_config(raw: dict, preset_name: str | None = None) -> ProtocolConfig:
    preset_key = preset_name or str(raw.get("preset", "single_point"))
    data = preset_copy(preset_key)
    data.update(raw)
    return ProtocolConfig(**data)


_LARGE_AROMATIC_RESIDUES = frozenset({"Y", "W", "F"})
# Evidence-based floors (2026-08 fixed-chemistry overlap study): at 8 lambda
# windows the chemically complete hybrids show median overlap ~0.16 (<0.2
# threshold for 56% of legs); 16-window runs reached 0.28-0.54.
ADAPTIVE_LAMBDA_MIN_DEFAULT = 12
ADAPTIVE_LAMBDA_MIN_LARGE_AROMATIC = 16


def adaptive_lambda_windows(base_windows: int, sites: tuple) -> int:
    """Raise lambda windows based on mutation size: aromatic<->non-aromatic
    transformations (large sidechain deletion/insertion) need denser lambda
    spacing. Applies only when the protocol did not pin lambda_windows."""
    large_aromatic = any(
        (site.wt in _LARGE_AROMATIC_RESIDUES) != (site.mut in _LARGE_AROMATIC_RESIDUES)
        for site in sites
    )
    floor = ADAPTIVE_LAMBDA_MIN_LARGE_AROMATIC if large_aromatic else ADAPTIVE_LAMBDA_MIN_DEFAULT
    return max(int(base_windows), floor)


def choose_protocol_for_group(
    base_protocol: ProtocolConfig,
    mutation_count: int,
    *,
    explicit_overrides: set[str] | None = None,
) -> ProtocolConfig:
    if mutation_count == 1:
        return base_protocol
    if mutation_count == 2 and base_protocol.preset != "double_point":
        data = preset_copy("double_point")
        base_data = dataclass_to_dict(base_protocol)
        always_inherit = {
            "force_field",
            "water_model",
            "gmx_bin",
            "pmx_bin",
            "allow_external_execute",
            "equilibration_restraint_schedule",
            "equilibration_release_npt_ps",
            "equilibration_heavy_posres_fc_kj_mol_nm2",
            "equilibration_backbone_posres_fc_kj_mol_nm2",
            "box_type",
            "box_padding_nm",
            "nonbonded_cutoff_nm",
            "vdw_switch_nm",
            "equilibrate_retry_box_padding_nm",
            "equilibrate_emergency_box_padding_nm",
            "salt_concentration_m",
            "max_bar_stderr_kcal_mol",
            "production_dt_ps",
            "window_relax_em_steps",
            "window_relax_md_ps",
            "window_relax_md_dt_ps",
            "grompp_maxwarn_genion",
            "grompp_maxwarn_equilibration",
            "grompp_maxwarn_sampling",
            "mdrun_args",
            "equilibration_pressure_coupling",
            "equilibration_pressure_tau_ps",
            "equilibration_refcoord_scaling",
            "sampling_pressure_coupling",
            "sampling_pressure_tau_ps",
            "sampling_refcoord_scaling",
            "equilibrate_em_steps",
        }
        for key in always_inherit:
            data[key] = base_data[key]

        for key in explicit_overrides or set():
            if key == "preset":
                continue
            if key in base_data:
                data[key] = base_data[key]
        return ProtocolConfig(**data)
    return base_protocol


def build_batch_plan(
    system_path: Path,
    mutations_path: Path,
    protocol_path: Path,
    batch_id: str | None = None,
    runs_root: Path | None = None,
) -> BatchPlan:
    system = load_system_config(system_path)
    raw_protocol = read_yaml(protocol_path)
    base_protocol = load_protocol_config(protocol_path)
    explicit_protocol_fields = set(raw_protocol)
    project_paths = ProjectPaths.discover()
    runs_root = runs_root or project_paths.runs_root
    ensure_dir(runs_root)

    mutation_groups = []
    for group in load_mutation_groups_from_csv(mutations_path):
        validated = build_mutation_group(
            mutation_group_id=group["mutation_group_id"],
            sites=group["sites"],
            allow_double_same_side=True,
            allow_charge_change=bool(base_protocol.allow_charge_changing),
        )
        mutation_groups.append(validated)

    batch_id = batch_id or f"{slugify(system.system_name)}_{base_protocol.preset}_{utc_now().replace(':', '').replace('-', '')}"
    batch_dir = ensure_dir(Path(runs_root) / batch_id)
    jobs_dir = ensure_dir(batch_dir / "jobs")
    ensure_dir(batch_dir / "artifacts")
    ensure_dir(batch_dir / "reports")

    jobs: list[JobSpec] = []
    job_rows: list[dict[str, object]] = []
    for mutation_group in mutation_groups:
        protocol = choose_protocol_for_group(
            base_protocol,
            mutation_group.mutation_count,
            explicit_overrides=explicit_protocol_fields,
        )
        if "lambda_windows" not in explicit_protocol_fields:
            adapted = adaptive_lambda_windows(protocol.lambda_windows, mutation_group.sites)
            if adapted != protocol.lambda_windows:
                protocol = replace(protocol, lambda_windows=adapted)
        job_id = build_job_id(system.system_name, mutation_group)
        workdir = ensure_dir(jobs_dir / job_id)
        job = JobSpec(
            job_id=job_id,
            mutation_group=mutation_group,
            protocol=protocol,
            system=system,
            batch_id=batch_id,
            workdir=str(workdir),
        )
        jobs.append(job)
        ensure_dir(workdir / "config")
        ensure_dir(workdir / "stages")
        ensure_dir(workdir / "artifacts" / "commands")
        ensure_dir(workdir / "results")
        ensure_dir(workdir / "report")
        write_yaml(workdir / "config" / "system.yml", dataclass_to_dict(system))
        write_yaml(workdir / "config" / "protocol.yml", dataclass_to_dict(protocol))
        write_json(workdir / "config" / "mutation_group.json", dataclass_to_dict(mutation_group))
        write_json(workdir / "job_spec.json", dataclass_to_dict(job))
        job_rows.append(
            {
                "job_id": job_id,
                "mutation_group_id": mutation_group.mutation_group_id,
                "mutation_count": mutation_group.mutation_count,
                "entity_side": mutation_group.entity_side,
                "min_version": mutation_group.min_version,
                "protocol_preset": protocol.preset,
                "workdir": workdir,
            }
        )

    batch_plan = BatchPlan(
        batch_id=batch_id,
        system_name=system.system_name,
        batch_dir=str(batch_dir),
        jobs=tuple(jobs),
    )
    write_json(batch_dir / "batch_plan.json", dataclass_to_dict(batch_plan))
    write_yaml(batch_dir / "batch_plan.yml", dataclass_to_dict(batch_plan))
    write_csv_rows(
        batch_dir / "jobs.csv",
        job_rows,
        ["job_id", "mutation_group_id", "mutation_count", "entity_side", "min_version", "protocol_preset", "workdir"],
    )
    return batch_plan
