"""Datamodels for manifests, plans, and stage state."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class MutationSite:
    chain_id: str
    resseq: int
    wt: str
    mut: str
    entity_side: str
    icode: str = ""

    def token(self) -> str:
        insertion = self.icode or ""
        return f"{self.chain_id}:{self.wt}{self.resseq}{insertion}{self.mut}@{self.entity_side}"

    def short_label(self) -> str:
        insertion = self.icode or ""
        return f"{self.chain_id.lower()}-{self.wt.lower()}{self.resseq}{insertion.lower()}{self.mut.lower()}"


@dataclass(frozen=True)
class MutationGroup:
    mutation_group_id: str
    sites: tuple[MutationSite, ...]
    mutation_count: int
    entity_side: str
    charge_conserving: bool
    min_version: str

    def signature(self) -> str:
        return "__".join(site.token() for site in self.sites)

    def short_label(self) -> str:
        return "--".join(site.short_label() for site in self.sites)


@dataclass(frozen=True)
class SystemConfig:
    system_name: str
    input_structure: str
    structure_source: str
    antibody_chains: tuple[str, ...]
    antigen_chains: tuple[str, ...]
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class ProtocolConfig:
    preset: str
    force_field: str
    water_model: str
    lambda_windows: int
    repeats: int
    nvt_ps: int
    npt_ps: int
    production_ps: int
    production_dt_ps: float
    temperature_k: float
    pressure_bar: float
    overlap_threshold: float
    max_repeat_delta_kcal_mol: float
    max_bar_stderr_kcal_mol: float
    window_relax_em_steps: int
    window_relax_md_ps: float
    window_relax_md_dt_ps: float
    same_side_double_point: bool
    allow_charge_changing: bool
    allow_cross_side_double_point: bool
    equilibration_restraint_schedule: str = "legacy_posres"
    equilibration_release_npt_ps: int = 0
    equilibration_heavy_posres_fc_kj_mol_nm2: float = 1000.0
    equilibration_backbone_posres_fc_kj_mol_nm2: float = 250.0
    gmx_bin: str = "gmx"
    pmx_bin: str = "pmx"
    allow_external_execute: bool = False
    box_type: str = "dodecahedron"
    box_padding_nm: float = 1.0
    nonbonded_cutoff_nm: float = 1.25
    vdw_switch_nm: float = 1.0
    equilibrate_retry_box_padding_nm: float = 2.0
    equilibrate_emergency_box_padding_nm: float = 5.0
    salt_concentration_m: float = 0.15
    grompp_maxwarn_genion: int = 1
    grompp_maxwarn_equilibration: int = 1
    grompp_maxwarn_sampling: int = 1
    mdrun_args: str = "-ntmpi 1 -ntomp 1"
    equilibration_pressure_coupling: str = "C-rescale"
    equilibration_pressure_tau_ps: float = 5.0
    equilibration_refcoord_scaling: str = "com"
    sampling_pressure_coupling: str = "C-rescale"
    sampling_pressure_tau_ps: float = 5.0
    sampling_refcoord_scaling: str = "all"
    equilibrate_em_steps: int = 5000


@dataclass(frozen=True)
class JobSpec:
    job_id: str
    mutation_group: MutationGroup
    protocol: ProtocolConfig
    system: SystemConfig
    batch_id: str
    workdir: str


@dataclass(frozen=True)
class BatchPlan:
    batch_id: str
    system_name: str
    batch_dir: str
    jobs: tuple[JobSpec, ...]


@dataclass
class StageStatus:
    stage: str
    state: str
    message: str
    commands: list[str] = field(default_factory=list)
    artifacts: list[str] = field(default_factory=list)
    started_at: str | None = None
    completed_at: str | None = None


def dataclass_to_dict(value: Any) -> Any:
    """Convert dataclasses recursively into plain serializable structures."""

    if hasattr(value, "__dataclass_fields__"):
        return asdict(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [dataclass_to_dict(item) for item in value]
    if isinstance(value, list):
        return [dataclass_to_dict(item) for item in value]
    if isinstance(value, dict):
        return {key: dataclass_to_dict(item) for key, item in value.items()}
    return value
