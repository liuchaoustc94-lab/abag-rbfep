"""Stage orchestration for abag-rbfep jobs."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha1
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

from abag_pmx.mutations import mutation_script_lines
from abag_rbfe.constants import STAGES
from abag_rbfe.execution import CommandOutcome, CommandRunner, discover_visible_gpu_devices
from abag_rbfe.gmx import (
    discover_gmx_top_dir,
    ensure_local_gmxlib,
    gro_file_is_valid,
    inspect_gro_file,
    materialize_staged_equilibration_restraints,
    resolve_gmx_binary,
    validate_hybrid_topology_integrity,
    water_coordinate_path,
)
from abag_rbfe.io_utils import read_json, read_yaml, utc_now, write_json
from abag_rbfe.models import JobSpec, MutationGroup, MutationSite, ProtocolConfig, StageStatus, SystemConfig, dataclass_to_dict
from abag_rbfe.paths import ProjectPaths
from abag_rbfe.planning import hydrate_protocol_config
from abag_rbfe.reporting import write_job_results, write_job_summary
from abag_rbfe.structure import (
    classify_incomplete_standard_residues,
    empty_repair_summary,
    extract_pdb_chains,
    find_inter_residue_heavy_atom_clashes,
    find_intra_residue_heavy_atom_clashes,
    pdb_chain_ids,
    partition_inter_residue_sidechain_repairable_clashes,
    partition_sidechain_repairable_clashes,
    repair_missing_atoms_with_pdbfixer,
    repair_incomplete_standard_residues_with_pdbfixer,
    restore_incomplete_standard_residues_from_template,
    strip_hydrogen_atoms,
    strip_terminal_oxygen_atoms,
    strip_sidechain_atoms_for_residues,
    write_inter_residue_heavy_atom_clash_report,
)

JOB_EXECUTION_LOCK = ".abag_job_execution_lock.json"
TOPOLOGY_INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"')


@dataclass(frozen=True)
class StageContext:
    job_dir: Path
    job: JobSpec
    system: SystemConfig
    protocol: ProtocolConfig
    mutation_group: MutationGroup
    runner: CommandRunner
    project_paths: ProjectPaths
    rescue_config: dict[str, object]
    env_overrides: dict[str, str] | None = None


def _load_optional_json_dict(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        payload = read_json(path)
    except OSError:
        return {}
    return payload if isinstance(payload, dict) else {}


def _load_context(job_dir: Path, execute: bool, environment: dict[str, str] | None = None) -> StageContext:
    spec = read_json(job_dir / "job_spec.json")
    system = SystemConfig(**read_yaml(job_dir / "config" / "system.yml"))
    protocol = hydrate_protocol_config(read_yaml(job_dir / "config" / "protocol.yml"))
    rescue_config = _load_optional_json_dict(job_dir / "config" / "rescue.json")
    mutation_group = MutationGroup(
        mutation_group_id=spec["mutation_group"]["mutation_group_id"],
        mutation_count=spec["mutation_group"]["mutation_count"],
        entity_side=spec["mutation_group"]["entity_side"],
        charge_conserving=spec["mutation_group"]["charge_conserving"],
        min_version=spec["mutation_group"]["min_version"],
        sites=tuple(MutationSite(**site) for site in spec["mutation_group"]["sites"]),
    )
    job = JobSpec(
        job_id=spec["job_id"],
        mutation_group=mutation_group,
        protocol=protocol,
        system=system,
        batch_id=spec["batch_id"],
        workdir=spec["workdir"],
    )
    return StageContext(
        job_dir=job_dir,
        job=job,
        system=system,
        protocol=protocol,
        mutation_group=mutation_group,
        runner=CommandRunner(execute=execute),
        project_paths=ProjectPaths.discover(),
        rescue_config=rescue_config,
        env_overrides=dict(environment) if environment else None,
    )


def _write_stage_status(job_dir: Path, status: StageStatus) -> StageStatus:
    payload = dataclass_to_dict(status)
    write_json(job_dir / "stages" / f"{status.stage}.json", payload)
    return status


def _external_stage_artifacts(script_file: Path, artifacts: list[str]) -> list[str]:
    ordered: list[str] = []
    for candidate in [str(script_file), str(script_file.with_suffix(".log")), *artifacts]:
        if candidate not in ordered:
            ordered.append(candidate)
    return ordered


def _write_external_stage_running(
    job_dir: Path,
    *,
    stage: str,
    script_file: Path,
    commands: list[str],
    artifacts: list[str],
    started_at: str,
    message: str | None = None,
) -> StageStatus:
    return _write_stage_status(
        job_dir,
        StageStatus(
            stage=stage,
            state="running",
            message=message or f"Stage execution started. Script: {script_file.name}",
            commands=commands,
            artifacts=_external_stage_artifacts(script_file, artifacts),
            started_at=started_at,
            completed_at=None,
        ),
    )


def _leg_mutated_chains(ctx: StageContext, leg: str) -> list[str]:
    if leg == "complex":
        return list(ctx.system.antibody_chains + ctx.system.antigen_chains)
    if ctx.mutation_group.entity_side == "antibody":
        return list(ctx.system.antibody_chains)
    return list(ctx.system.antigen_chains)


def _pmx_mutff_root(ctx: StageContext) -> Path:
    return ctx.project_paths.vendor_pmx_root / "src" / "pmx" / "data" / "mutff"


def _resolve_pmx_command(ctx: StageContext) -> str | None:
    pmx_argv = _resolve_pmx_argv(ctx)
    if pmx_argv is None:
        return None
    return " ".join(sh_quote(token) for token in pmx_argv)


def _resolve_pmx_argv(ctx: StageContext) -> list[str] | None:
    project_python = ctx.project_paths.repo_root / ".venv" / "bin" / "python"
    if project_python.is_file() and _python_has_pmx(project_python):
        return [str(project_python), "-m", "pmx.scripts.cli"]

    pmx_binary = shutil.which(ctx.protocol.pmx_bin)
    if pmx_binary:
        return [pmx_binary]

    candidate = Path(ctx.protocol.pmx_bin).expanduser()
    if candidate.is_file():
        return [str(candidate.resolve())]
    return None


def _resolve_gmx_command(ctx: StageContext) -> str:
    binary = resolve_gmx_binary(ctx.protocol.gmx_bin)
    if binary is not None:
        return sh_quote(str(binary))
    return sh_quote(ctx.protocol.gmx_bin)


def _resolve_project_python_command(ctx: StageContext) -> str:
    project_python = ctx.project_paths.repo_root / ".venv" / "bin" / "python"
    if project_python.is_file():
        return sh_quote(str(project_python))
    fallback = shutil.which("python3") or shutil.which("python") or "python3"
    return sh_quote(fallback)


def _prepare_local_gmxlib(ctx: StageContext) -> Path | None:
    gmx_top_dir = discover_gmx_top_dir(ctx.protocol.gmx_bin)
    if gmx_top_dir is None:
        return None
    return ensure_local_gmxlib(
        job_dir=ctx.job_dir,
        gmx_top_dir=gmx_top_dir,
        pmx_mutff_root=_pmx_mutff_root(ctx),
        force_field=ctx.protocol.force_field,
    )


def _stage_env(ctx: StageContext) -> dict[str, str] | None:
    env: dict[str, str] = dict(ctx.env_overrides or {})
    gmxlib_dir = _prepare_local_gmxlib(ctx)
    if gmxlib_dir is not None:
        env["GMXLIB"] = str(gmxlib_dir)

    if "GMX_MAXBACKUP" not in env and "GMX_MAXBACKUP" not in os.environ:
        env["GMX_MAXBACKUP"] = "-1"

    if "CUDA_VISIBLE_DEVICES" not in env and "CUDA_VISIBLE_DEVICES" not in os.environ:
        visible_devices = discover_visible_gpu_devices()
        if visible_devices:
            selected_device = visible_devices[(_job_seed(ctx, "gpu-device") - 1) % len(visible_devices)]
            env["CUDA_VISIBLE_DEVICES"] = selected_device

    return env or None


def _python_has_pmx(python_executable: Path) -> bool:
    result = subprocess.run(
        [
            str(python_executable),
            "-c",
            "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('pmx.scripts.cli') else 1)",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def _lambda_values(window_count: int) -> list[float]:
    if window_count <= 1:
        return [0.0]
    return [index / (window_count - 1) for index in range(window_count)]


_DECOUPLED_LAMBDA_MIN_WINDOWS = 6


def _lambda_component_schedules(window_count: int) -> tuple[list[float], list[float]] | None:
    """Coulomb-first / vdW-second decoupled lambda schedule.

    Decoupling electrostatics from vdW avoids the overlap collapse observed
    for charge-changing transformations (e.g. charged-His hybrids showed
    overlap 0.064 with a coupled single-lambda schedule). Returns
    (coul_lambdas, vdw_lambdas); bonded/mass follow the vdW leg. Small window
    counts (smoke presets) keep the coupled single-vector behaviour."""
    if window_count < _DECOUPLED_LAMBDA_MIN_WINDOWS:
        return None
    coul_windows = max(2, window_count // 3)
    coul: list[float] = []
    vdw: list[float] = []
    for index in range(window_count):
        coul.append(min(1.0, index / (coul_windows - 1)))
        if index < coul_windows - 1:
            vdw.append(0.0)
        else:
            vdw.append((index - coul_windows + 1) / (window_count - coul_windows))
    return coul, vdw


def _format_lambda_schedule_lines(lambda_values: list[float]) -> list[str]:
    schedules = _lambda_component_schedules(len(lambda_values))
    if schedules is None:
        return [f"fep-lambdas             = {_format_lambda_values(lambda_values)}"]
    coul, vdw = schedules
    vdw_text = _format_lambda_values(vdw)
    return [
        f"coul-lambdas            = {_format_lambda_values(coul)}",
        f"vdw-lambdas             = {vdw_text}",
        f"bonded-lambdas          = {vdw_text}",
        f"mass-lambdas            = {vdw_text}",
    ]


def _format_lambda_values(values: list[float]) -> str:
    return " ".join(f"{value:.5f}" for value in values)


def _steps(duration_ps: float, dt_ps: float = 0.002) -> int:
    return max(int(round(duration_ps / dt_ps)), 1)


def _stable_seed(*tokens: object) -> int:
    digest = sha1("|".join(str(token) for token in tokens).encode("utf-8")).hexdigest()
    return int(digest[:8], 16) % 2147483646 + 1


def _job_seed(ctx: StageContext, *tokens: object) -> int:
    return _stable_seed(ctx.job.batch_id, ctx.job.job_id, *tokens)


def _context_env_value(ctx: StageContext, key: str) -> str | None:
    if ctx.env_overrides and key in ctx.env_overrides:
        return ctx.env_overrides[key]
    return os.environ.get(key)


def _configured_legs(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, (list, tuple, set)):
        items = [str(item).strip().lower() for item in value]
    else:
        items = [item.strip().lower() for item in str(value).split(",")]
    legs: list[str] = []
    for item in items:
        if item not in {"complex", "apo"} or item in legs:
            continue
        legs.append(item)
    return tuple(legs)


def _rescue_source_job_dir(ctx: StageContext) -> Path | None:
    raw_source_job_dir = str(ctx.rescue_config.get("source_job_dir") or "").strip()
    if not raw_source_job_dir:
        return None
    source_job_dir = Path(raw_source_job_dir).expanduser()
    try:
        source_job_dir = source_job_dir.resolve()
    except OSError:
        pass
    if not source_job_dir.is_dir():
        return None
    return source_job_dir


def _job_target_legs(ctx: StageContext) -> tuple[str, ...]:
    target_legs = _configured_legs(ctx.rescue_config.get("target_legs"))
    return target_legs or ("complex", "apo")


def _job_inherit_source_legs(ctx: StageContext) -> tuple[str, ...]:
    inherit_source_legs = _configured_legs(ctx.rescue_config.get("inherit_source_legs"))
    if inherit_source_legs:
        return inherit_source_legs
    targeted = set(_job_target_legs(ctx))
    return tuple(leg for leg in ("complex", "apo") if leg not in targeted)


def _should_inherit_leg_from_source(ctx: StageContext, leg: str) -> bool:
    return _rescue_source_job_dir(ctx) is not None and leg in _job_inherit_source_legs(ctx)


def _should_seed_equilibrated_repeat_from_source(ctx: StageContext, leg: str) -> bool:
    if _rescue_source_job_dir(ctx) is not None:
        return leg in _job_inherit_source_legs(ctx)
    return bool(_equilibration_seed_job_dirs(ctx))


def _split_path_env(value: str | None) -> list[Path]:
    if value is None:
        return []
    ordered: list[Path] = []
    seen: set[Path] = set()
    normalized = value.replace("\n", os.pathsep)
    for raw_token in normalized.split(os.pathsep):
        token = raw_token.strip()
        if not token:
            continue
        candidate = Path(token).expanduser()
        try:
            candidate = candidate.resolve()
        except OSError:
            pass
        if candidate in seen:
            continue
        seen.add(candidate)
        ordered.append(candidate)
    return ordered


def _seed_job_dir_matches_target(ctx: StageContext, seed_job_dir: Path) -> bool:
    spec_path = seed_job_dir / "job_spec.json"
    if not spec_path.is_file():
        return False
    try:
        payload = read_json(spec_path)
    except OSError:
        return False
    if str(payload.get("job_id") or "") != ctx.job.job_id:
        return False
    mutation_group = payload.get("mutation_group", {})
    if isinstance(mutation_group, dict):
        seed_group_id = str(mutation_group.get("mutation_group_id") or "").strip()
        if seed_group_id and seed_group_id != ctx.mutation_group.mutation_group_id:
            return False
    return True


def _equilibration_seed_job_dirs(ctx: StageContext) -> list[Path]:
    ordered: list[Path] = []
    seen: set[Path] = set()

    def register(candidate: Path) -> None:
        if candidate in seen or not candidate.is_dir():
            return
        if not _seed_job_dir_matches_target(ctx, candidate):
            return
        seen.add(candidate)
        ordered.append(candidate)

    rescue_source_job_dir = _rescue_source_job_dir(ctx)
    if rescue_source_job_dir is not None:
        register(rescue_source_job_dir)

    for candidate in _split_path_env(_context_env_value(ctx, "ABAG_RBFE_EQUILIBRATION_SEED_JOB_DIRS")):
        register(candidate)

    for root in _split_path_env(_context_env_value(ctx, "ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS")):
        register(root / "jobs" / ctx.job.job_id)
        if root.is_dir():
            for candidate in sorted(root.glob(f"*/jobs/{ctx.job.job_id}")):
                register(candidate)

    return ordered


def _split_shell_args(value: str) -> list[str]:
    if not value.strip():
        return []
    try:
        return shlex.split(value)
    except ValueError:
        return value.split()


def _truthy_env(value: str | None) -> bool:
    if value is None:
        return False
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _mdrun_suffix(ctx: StageContext) -> str:
    args_override = (_context_env_value(ctx, "ABAG_RBFE_MDRUN_ARGS") or "").strip()
    args = args_override or ctx.protocol.mdrun_args.strip()
    tokens = _split_shell_args(args)
    extra_tokens: list[str] = []

    pin_offset = (_context_env_value(ctx, "ABAG_RBFE_MDRUN_PINOFFSET") or "").strip()
    pin_stride = (_context_env_value(ctx, "ABAG_RBFE_MDRUN_PINSTRIDE") or "").strip()
    pin_requested = _truthy_env(_context_env_value(ctx, "ABAG_RBFE_MDRUN_PIN")) or bool(pin_offset) or bool(pin_stride)

    if pin_requested and "-pin" not in tokens:
        extra_tokens.extend(["-pin", "on"])
    if pin_offset and "-pinoffset" not in tokens:
        extra_tokens.extend(["-pinoffset", pin_offset])
    if pin_stride and "-pinstride" not in tokens:
        extra_tokens.extend(["-pinstride", pin_stride])

    if extra_tokens:
        args = " ".join(part for part in [args, *[shlex.quote(token) for token in extra_tokens]] if part)
    return f" {args}" if args else ""


def _shell_all_nonempty(paths: list[Path]) -> str:
    return " && ".join(f"[ -s {sh_quote(str(path))} ]" for path in paths)


def _parse_topology_include_target(raw_line: str) -> str | None:
    line = raw_line.split(";", 1)[0].strip()
    if not line:
        return None
    match = TOPOLOGY_INCLUDE_RE.match(line)
    if match is None:
        return None
    return str(match.group(1) or "").strip() or None


def _local_topology_support_paths(topology_path: Path) -> list[Path]:
    if not topology_path.is_file():
        return []
    topology_dir = topology_path.parent.resolve()
    support_paths: list[Path] = []
    seen: set[Path] = set()
    for raw_line in topology_path.read_text(encoding="utf-8").splitlines():
        include_target = _parse_topology_include_target(raw_line)
        if include_target is None or not include_target.lower().endswith(".itp"):
            continue
        include_path = Path(include_target)
        resolved = include_path if include_path.is_absolute() else (topology_dir / include_path).resolve()
        if resolved.parent != topology_dir or resolved in seen:
            continue
        seen.add(resolved)
        support_paths.append(resolved)
    return support_paths


def _render_nonbonded_block(ctx: StageContext) -> list[str]:
    cutoff_nm = max(float(ctx.protocol.nonbonded_cutoff_nm), 0.1)
    switch_nm = max(min(float(ctx.protocol.vdw_switch_nm), cutoff_nm - 0.01), 0.0)
    return [
        "cutoff-scheme           = Verlet",
        "verlet-buffer-tolerance = -1",
        "nstlist                 = 20",
        f"rlist                   = {cutoff_nm:.2f}",
        "coulombtype             = PME",
        f"rcoulomb                = {cutoff_nm:.2f}",
        "vdw-modifier            = Force-switch",
        f"rvdw-switch             = {switch_nm:.2f}",
        f"rvdw                    = {cutoff_nm:.2f}",
    ]


def _render_genion_mdp(ctx: StageContext) -> str:
    return "\n".join(
        [
            "integrator              = steep",
            "nsteps                  = 200",
            "emtol                   = 1000",
            *_render_nonbonded_block(ctx),
            "pbc                     = xyz",
            "",
        ]
    )


def _render_em_mdp(ctx: StageContext) -> str:
    return "\n".join(
        [
            "define                  = -DFLEXIBLE",
            "integrator              = steep",
            f"nsteps                  = {ctx.protocol.equilibrate_em_steps}",
            "emtol                   = 1000",
            *_render_nonbonded_block(ctx),
            "constraints             = none",
            "pbc                     = xyz",
            "",
        ]
    )


def _equilibration_restraint_schedule(ctx: StageContext) -> str:
    normalized = str(ctx.protocol.equilibration_restraint_schedule or "legacy_posres").strip().lower()
    if normalized in {"staged_backbone_release", "staged_backbone"}:
        return "staged_backbone_release"
    return "legacy_posres"


def _uses_staged_equilibration_restraints(ctx: StageContext) -> bool:
    return _equilibration_restraint_schedule(ctx) == "staged_backbone_release"


def _equilibration_release_npt_ps(ctx: StageContext) -> int:
    return max(int(ctx.protocol.equilibration_release_npt_ps), 0)


def _has_equilibration_release_stage(ctx: StageContext) -> bool:
    return _equilibration_release_npt_ps(ctx) > 0


def _equilibration_nvt_define(ctx: StageContext) -> str:
    if _uses_staged_equilibration_restraints(ctx):
        return "-DPOSRES_STAGE_HEAVY"
    return "-DPOSRES"


def _equilibration_npt_define(ctx: StageContext) -> str:
    if _uses_staged_equilibration_restraints(ctx):
        return "-DPOSRES_STAGE_HEAVY"
    return "-DPOSRES"


def _equilibration_release_define(ctx: StageContext) -> str:
    if _uses_staged_equilibration_restraints(ctx):
        return "-DPOSRES_STAGE_BACKBONE"
    return "-DPOSRES"


def _render_nvt_mdp(ctx: StageContext, seed: int) -> str:
    return "\n".join(
        [
            f"define                  = {_equilibration_nvt_define(ctx)}",
            "integrator              = md",
            "dt                      = 0.002",
            f"nsteps                  = {_steps(ctx.protocol.nvt_ps)}",
            *_render_nonbonded_block(ctx),
            "DispCorr                = EnerPres",
            "tcoupl                  = v-rescale",
            "tc-grps                 = System",
            "tau-t                   = 1.0",
            f"ref-t                   = {ctx.protocol.temperature_k}",
            "pcoupl                  = no",
            "constraints             = h-bonds",
            "constraint-algorithm    = lincs",
            "continuation            = no",
            "gen-vel                 = yes",
            f"gen-temp                = {ctx.protocol.temperature_k}",
            f"gen-seed                = {seed}",
            "pbc                     = xyz",
            "nstenergy               = 500",
            "nstlog                  = 500",
            "nstxout-compressed      = 1000",
            "",
        ]
    )


def _render_npt_stage_mdp(ctx: StageContext, *, define: str, duration_ps: int) -> str:
    pressure_coupling = str(ctx.protocol.equilibration_pressure_coupling).strip() or "C-rescale"
    return "\n".join(
        [
            f"define                  = {define}",
            "integrator              = md",
            "dt                      = 0.002",
            f"nsteps                  = {_steps(duration_ps)}",
            *_render_nonbonded_block(ctx),
            "DispCorr                = EnerPres",
            "tcoupl                  = v-rescale",
            "tc-grps                 = System",
            "tau-t                   = 1.0",
            f"ref-t                   = {ctx.protocol.temperature_k}",
            f"pcoupl                  = {pressure_coupling}",
            "pcoupltype              = isotropic",
            f"tau-p                   = {ctx.protocol.equilibration_pressure_tau_ps}",
            "compressibility         = 4.5e-5",
            f"ref-p                   = {ctx.protocol.pressure_bar}",
            "constraints             = h-bonds",
            "constraint-algorithm    = lincs",
            "continuation            = yes",
            "gen-vel                 = no",
            "pbc                     = xyz",
            "nstenergy               = 500",
            "nstlog                  = 500",
            "nstxout-compressed      = 1000",
            f"refcoord-scaling        = {ctx.protocol.equilibration_refcoord_scaling}",
            "",
        ]
    )


def _render_npt_mdp(ctx: StageContext) -> str:
    return _render_npt_stage_mdp(
        ctx,
        define=_equilibration_npt_define(ctx),
        duration_ps=ctx.protocol.npt_ps,
    )


def _render_npt_release_mdp(ctx: StageContext) -> str:
    return _render_npt_stage_mdp(
        ctx,
        define=_equilibration_release_define(ctx),
        duration_ps=_equilibration_release_npt_ps(ctx),
    )


def _write_equilibration_mdps(ctx: StageContext, *, leg: str, repeat_index: int, repeat_dir: Path) -> list[str]:
    mdp_dir = repeat_dir / "mdp"
    mdp_dir.mkdir(parents=True, exist_ok=True)
    artifacts = [
        mdp_dir / "genion.mdp",
        mdp_dir / "em.mdp",
        mdp_dir / "nvt.mdp",
        mdp_dir / "npt.mdp",
    ]
    (mdp_dir / "genion.mdp").write_text(_render_genion_mdp(ctx), encoding="utf-8")
    (mdp_dir / "em.mdp").write_text(_render_em_mdp(ctx), encoding="utf-8")
    (mdp_dir / "nvt.mdp").write_text(
        _render_nvt_mdp(ctx, seed=_job_seed(ctx, leg, repeat_index, "nvt")),
        encoding="utf-8",
    )
    (mdp_dir / "npt.mdp").write_text(_render_npt_mdp(ctx), encoding="utf-8")

    release_path = mdp_dir / "npt_release.mdp"
    if _has_equilibration_release_stage(ctx):
        release_path.write_text(_render_npt_release_mdp(ctx), encoding="utf-8")
        artifacts.append(release_path)
    elif release_path.exists():
        release_path.unlink()

    return [str(path) for path in artifacts]


def _equilibrate_retry_attempts(ctx: StageContext) -> list[tuple[str, float, str]]:
    attempts: list[tuple[str, float, str]] = [
        (ctx.protocol.box_type, float(ctx.protocol.box_padding_nm), "requested"),
    ]
    retry_candidates = [
        ("cubic", max(float(ctx.protocol.box_padding_nm), float(ctx.protocol.equilibrate_retry_box_padding_nm)), "fallback"),
        (
            "cubic",
            max(float(ctx.protocol.box_padding_nm), float(ctx.protocol.equilibrate_emergency_box_padding_nm)),
            "expanded fallback",
        ),
    ]
    seen = {(ctx.protocol.box_type, round(float(ctx.protocol.box_padding_nm), 3))}
    for box_type, padding_nm, label in retry_candidates:
        signature = (box_type, round(padding_nm, 3))
        if signature in seen:
            continue
        attempts.append((box_type, padding_nm, label))
        seen.add(signature)
    return attempts


def _equilibrate_repeat_snippet(
    *,
    ctx: StageContext,
    gmx_command: str,
    mdrun_suffix: str,
    processed_gro: Path,
    pmx_top: Path,
    support_itps: list[Path],
    repeat_dir: Path,
    setup_dir: Path,
    equil_dir: Path,
    mdp_dir: Path,
    solvent_coordinates: Path,
) -> str:
    project_python_command = _resolve_project_python_command(ctx)
    use_release_stage = _has_equilibration_release_stage(ctx)
    repeat_top = repeat_dir / "system.top"
    repeat_support_itps = [repeat_dir / support_itp.name for support_itp in support_itps]
    boxed_gro = setup_dir / "boxed.gro"
    solvated_gro = setup_dir / "solvated.gro"
    genion_tpr = setup_dir / "genion.tpr"
    ions_gro = setup_dir / "ions.gro"
    em_tpr = equil_dir / "em.tpr"
    nvt_tpr = equil_dir / "nvt.tpr"
    npt_tpr = equil_dir / "npt.tpr"
    npt_stage1_tpr = equil_dir / "npt_stage1.tpr"
    npt_release_tpr = equil_dir / "npt_release.tpr"
    em_deffnm = equil_dir / "em"
    nvt_deffnm = equil_dir / "nvt"
    npt_deffnm = equil_dir / "npt"
    npt_stage1_deffnm = equil_dir / "npt_stage1"
    npt_release_deffnm = equil_dir / "npt_release"
    em_runtime_log = equil_dir / "em.runtime.log"
    em_runtime_history_log = equil_dir / "em.runtime.history.log"
    restraint_summary_path = repeat_dir / "equilibration_restraints.json"
    completion_check = _shell_all_nonempty([repeat_top, *repeat_support_itps, npt_deffnm.with_suffix(".gro")])
    repeat_tag = f"{repeat_dir.parent.name}/{repeat_dir.name}"
    retry_attempts = _equilibrate_retry_attempts(ctx)
    snippet_name = f"{repeat_dir.parent.name}_{repeat_dir.name}".replace("-", "_")
    shift_pattern = "inconsistent shifts over periodic boundaries"
    unstable_force_pattern = "force on at least one atom is not finite"
    infinite_force_pattern = "Maximum force     =            inf"
    excluded_distance_pattern = "largest distance between excluded atoms"
    staged_restraints_command = ""
    if _uses_staged_equilibration_restraints(ctx):
        staged_restraints_command = (
            "  "
            + _materialize_staged_equilibration_restraints_command(
                project_python_command,
                repeat_topology_path=repeat_top,
                summary_path=restraint_summary_path,
                heavy_force_constant=ctx.protocol.equilibration_heavy_posres_fc_kj_mol_nm2,
                backbone_force_constant=ctx.protocol.equilibration_backbone_posres_fc_kj_mol_nm2,
            )
            + " || return $?"
        )
    first_npt_tpr = npt_stage1_tpr if use_release_stage else npt_tpr
    first_npt_deffnm = npt_stage1_deffnm if use_release_stage else npt_deffnm

    copy_lines = [f"cp {sh_quote(str(pmx_top))} {sh_quote(str(repeat_top))} || return $?"]
    copy_lines.extend(
        f"cp {sh_quote(str(support_itp))} {sh_quote(str(repeat_dir / support_itp.name))} || return $?"
        for support_itp in support_itps
    )
    copy_block = "\n".join(f"  {line}" for line in copy_lines)
    sanitize_topology_command = (
        f"  {project_python_command} -c "
        + sh_quote(
            "from pathlib import Path; "
            "from abag_rbfe.gmx import deduplicate_standard_topology_includes; "
            f"deduplicate_standard_topology_includes(Path({str(repeat_top)!r}))"
        )
        + " || return $?"
    )
    cleanup_targets = [
        boxed_gro,
        solvated_gro,
        genion_tpr,
        ions_gro,
        em_tpr,
        nvt_tpr,
        npt_tpr,
        npt_stage1_tpr,
        npt_release_tpr,
        em_deffnm.with_suffix(".gro"),
        em_deffnm.with_suffix(".edr"),
        em_deffnm.with_suffix(".log"),
        em_deffnm.with_suffix(".trr"),
        nvt_deffnm.with_suffix(".gro"),
        nvt_deffnm.with_suffix(".edr"),
        nvt_deffnm.with_suffix(".log"),
        nvt_deffnm.with_suffix(".trr"),
        npt_deffnm.with_suffix(".gro"),
        npt_deffnm.with_suffix(".edr"),
        npt_deffnm.with_suffix(".log"),
        npt_deffnm.with_suffix(".trr"),
        npt_stage1_deffnm.with_suffix(".gro"),
        npt_stage1_deffnm.with_suffix(".edr"),
        npt_stage1_deffnm.with_suffix(".log"),
        npt_stage1_deffnm.with_suffix(".trr"),
        npt_release_deffnm.with_suffix(".gro"),
        npt_release_deffnm.with_suffix(".edr"),
        npt_release_deffnm.with_suffix(".log"),
        npt_release_deffnm.with_suffix(".trr"),
        em_runtime_log,
        restraint_summary_path,
    ]
    cleanup_line = "rm -f " + " ".join(sh_quote(str(path)) for path in cleanup_targets)
    attempt_specs = [
        f"  {sh_quote(f'{box_type}|{padding_nm:.2f}|{label}')}"
        for box_type, padding_nm, label in retry_attempts
    ]
    attempt_specs_block = "\n".join(attempt_specs)
    snippet_lines = [
        f"run_equilibrate_em_{snippet_name}() {{",
        '  local abag_box_type="$1"',
        '  local abag_box_padding="$2"',
        f"  {cleanup_line}",
        f"  : > {sh_quote(str(em_runtime_log))}",
        copy_block,
        sanitize_topology_command,
        *([staged_restraints_command] if staged_restraints_command else []),
        f"  {gmx_command} editconf -f {sh_quote(str(processed_gro))} -o {sh_quote(str(boxed_gro))} -bt \"$abag_box_type\" -d \"$abag_box_padding\" || return $?",
        f"  {gmx_command} solvate -cp {sh_quote(str(boxed_gro))} -cs {sh_quote(str(solvent_coordinates))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(solvated_gro))} || return $?",
        f"  {gmx_command} grompp -f {sh_quote(str(mdp_dir / 'genion.mdp'))} -c {sh_quote(str(solvated_gro))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(genion_tpr))} -maxwarn {ctx.protocol.grompp_maxwarn_genion} || return $?",
        f"  printf 'SOL\\n' | {gmx_command} genion -s {sh_quote(str(genion_tpr))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(ions_gro))} -neutral -conc {ctx.protocol.salt_concentration_m} || return $?",
        f"  {gmx_command} grompp -f {sh_quote(str(mdp_dir / 'em.mdp'))} -c {sh_quote(str(ions_gro))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(em_tpr))} -maxwarn {ctx.protocol.grompp_maxwarn_equilibration} || return $?",
        f"  printf '[abag-rbfep] EM attempt box=%s padding=%s nm\\n' \"$abag_box_type\" \"$abag_box_padding\" > {sh_quote(str(em_runtime_log))}",
        "  set +e",
        f"  {gmx_command} mdrun -s {sh_quote(str(em_tpr))} -deffnm {sh_quote(str(em_deffnm))}{mdrun_suffix} >> {sh_quote(str(em_runtime_log))} 2>&1",
        "  local abag_em_rc=$?",
        "  set -e",
        f"  if grep -q {sh_quote(unstable_force_pattern)} {sh_quote(str(em_runtime_log))} || grep -q {sh_quote(infinite_force_pattern)} {sh_quote(str(em_runtime_log))}; then",
        "    abag_em_rc=1",
        "  fi",
        '  if [ "$abag_em_rc" -eq 0 ]; then',
        "    set +e",
        f"    {gmx_command} grompp -f {sh_quote(str(mdp_dir / 'nvt.mdp'))} -c {sh_quote(str(em_deffnm.with_suffix('.gro')))} -r {sh_quote(str(em_deffnm.with_suffix('.gro')))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(nvt_tpr))} -maxwarn {ctx.protocol.grompp_maxwarn_equilibration} >> {sh_quote(str(em_runtime_log))} 2>&1",
        "    abag_em_rc=$?",
        "    set -e",
        "  fi",
        '  if [ "$abag_em_rc" -eq 0 ]; then',
        "    set +e",
        f"    {gmx_command} mdrun -s {sh_quote(str(nvt_tpr))} -deffnm {sh_quote(str(nvt_deffnm))}{mdrun_suffix} >> {sh_quote(str(em_runtime_log))} 2>&1",
        "    abag_em_rc=$?",
        "    set -e",
        "  fi",
        '  if [ "$abag_em_rc" -eq 0 ]; then',
        "    set +e",
        f"    {gmx_command} grompp -f {sh_quote(str(mdp_dir / 'npt.mdp'))} -c {sh_quote(str(nvt_deffnm.with_suffix('.gro')))} -r {sh_quote(str(nvt_deffnm.with_suffix('.gro')))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(first_npt_tpr))} -maxwarn {ctx.protocol.grompp_maxwarn_equilibration} >> {sh_quote(str(em_runtime_log))} 2>&1",
        "    abag_em_rc=$?",
        "    set -e",
        "  fi",
        '  if [ "$abag_em_rc" -eq 0 ]; then',
        "    set +e",
        f"    {gmx_command} mdrun -s {sh_quote(str(first_npt_tpr))} -deffnm {sh_quote(str(first_npt_deffnm))}{mdrun_suffix} >> {sh_quote(str(em_runtime_log))} 2>&1",
        "    abag_em_rc=$?",
        "    set -e",
        "  fi",
        *(
            [
                '  if [ "$abag_em_rc" -eq 0 ]; then',
                "    set +e",
                f"    {gmx_command} grompp -f {sh_quote(str(mdp_dir / 'npt_release.mdp'))} -c {sh_quote(str(npt_stage1_deffnm.with_suffix('.gro')))} -r {sh_quote(str(npt_stage1_deffnm.with_suffix('.gro')))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(npt_tpr))} -maxwarn {ctx.protocol.grompp_maxwarn_equilibration} >> {sh_quote(str(em_runtime_log))} 2>&1",
                "    abag_em_rc=$?",
                "    set -e",
                "  fi",
                '  if [ "$abag_em_rc" -eq 0 ]; then',
                "    set +e",
                f"    {gmx_command} mdrun -s {sh_quote(str(npt_tpr))} -deffnm {sh_quote(str(npt_deffnm))}{mdrun_suffix} >> {sh_quote(str(em_runtime_log))} 2>&1",
                "    abag_em_rc=$?",
                "    set -e",
                "  fi",
            ]
            if use_release_stage
            else []
        ),
        f"  printf '[abag-rbfep] EM attempt box=%s padding=%s nm exit=%s\\n' \"$abag_box_type\" \"$abag_box_padding\" \"$abag_em_rc\" >> {sh_quote(str(em_runtime_history_log))}",
        f"  cat {sh_quote(str(em_runtime_log))} >> {sh_quote(str(em_runtime_history_log))}",
        f"  printf '\\n' >> {sh_quote(str(em_runtime_history_log))}",
        "  return $abag_em_rc",
        "}",
        "",
        f"rm -f {sh_quote(str(em_runtime_history_log))}",
        f"run_equilibrate_attempts_{snippet_name}() {{",
        "  local attempt_specs=(",
        attempt_specs_block,
        "  )",
        "  local attempt_idx=0",
        "  local attempt_total=${#attempt_specs[@]}",
        "  while true; do",
        '    IFS="|" read -r abag_box_type abag_box_padding abag_label <<< "${attempt_specs[$attempt_idx]}"',
        f"    if run_equilibrate_em_{snippet_name} \"$abag_box_type\" \"$abag_box_padding\"; then",
        "      return 0",
        "    fi",
        f"    if ! grep -q {sh_quote(shift_pattern)} {sh_quote(str(em_runtime_log))} && ! grep -q {sh_quote(unstable_force_pattern)} {sh_quote(str(em_runtime_log))} && ! grep -q {sh_quote(infinite_force_pattern)} {sh_quote(str(em_runtime_log))} && ! grep -q {sh_quote(excluded_distance_pattern)} {sh_quote(str(em_runtime_log))}; then",
        "      return 1",
        "    fi",
        "    attempt_idx=$((attempt_idx + 1))",
        '    if [ "$attempt_idx" -ge "$attempt_total" ]; then',
        "      return 1",
        "    fi",
        '    IFS="|" read -r abag_next_box_type abag_next_box_padding abag_next_label <<< "${attempt_specs[$attempt_idx]}"',
        f"    echo \"[abag-rbfep] retrying equilibration with ${'{'}abag_next_label{'}'} box ${'{'}abag_next_box_type{'}'} and padding ${'{'}abag_next_box_padding{'}'} nm for {repeat_dir.parent.name}/{repeat_dir.name}\" >&2",
        "  done",
        "}",
        "",
        f"run_equilibrate_attempts_{snippet_name}",
    ]
    return "\n".join(
        [
            f"if {completion_check}; then",
            f"  echo \"[abag-rbfep] skipping completed equilibrate repeat {repeat_tag}\"",
            "else",
            *[(f"  {line}" if line else "") for line in snippet_lines],
            "fi",
        ]
    )


def _seed_equilibrated_repeat_from_sources(
    ctx: StageContext,
    *,
    leg: str,
    repeat_dir: Path,
    repeat_top: Path,
    support_itps: list[Path],
) -> list[str]:
    seed_source_path = repeat_dir / "equilibration" / "seed_source.json"
    required_outputs = [
        repeat_top,
        repeat_dir / "equilibration" / "npt.gro",
        *[repeat_dir / support_itp.name for support_itp in support_itps],
    ]
    if required_outputs and all(_stage_artifact_exists(path) for path in required_outputs) and not _stage_artifact_exists(seed_source_path):
        return []

    preferred_seed_repeat_dir: Path | None = None
    if _stage_artifact_exists(seed_source_path):
        try:
            seed_source_payload = read_json(seed_source_path)
        except OSError:
            seed_source_payload = {}
        raw_seed_repeat_dir = str(seed_source_payload.get("seed_repeat_dir") or "").strip()
        if raw_seed_repeat_dir:
            preferred_seed_repeat_dir = Path(raw_seed_repeat_dir).expanduser()
            try:
                preferred_seed_repeat_dir = preferred_seed_repeat_dir.resolve()
            except OSError:
                pass

    target_npt = repeat_dir / "equilibration" / "npt.gro"
    sync_seed_bundle = not _stage_artifact_exists(target_npt) or _stage_artifact_exists(seed_source_path)
    repeat_id = repeat_dir.name
    for seed_job_dir in _equilibration_seed_job_dirs(ctx):
        seed_repeat_dir = seed_job_dir / "legs" / leg / repeat_id
        if preferred_seed_repeat_dir is not None:
            try:
                resolved_seed_repeat_dir = seed_repeat_dir.resolve()
            except OSError:
                resolved_seed_repeat_dir = seed_repeat_dir
            if resolved_seed_repeat_dir != preferred_seed_repeat_dir:
                continue
        seed_top = seed_repeat_dir / "system.top"
        seed_npt = seed_repeat_dir / "equilibration" / "npt.gro"
        seed_itps = sorted(seed_repeat_dir.glob("*.itp"))
        if not _stage_artifact_exists(seed_top) or not _stage_artifact_exists(seed_npt):
            continue
        if any(not _stage_artifact_exists(path) for path in seed_itps):
            continue

        copied: list[str] = []
        repeat_top.parent.mkdir(parents=True, exist_ok=True)
        if sync_seed_bundle or not _stage_artifact_exists(repeat_top):
            shutil.copy2(seed_top, repeat_top)
            copied.append(str(repeat_top))

        for seed_itp in seed_itps:
            target_itp = repeat_dir / seed_itp.name
            if sync_seed_bundle or not _stage_artifact_exists(target_itp):
                target_itp.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(seed_itp, target_itp)
                copied.append(str(target_itp))

        if sync_seed_bundle or not _stage_artifact_exists(target_npt):
            target_npt.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(seed_npt, target_npt)
            copied.append(str(target_npt))

        write_json(
            seed_source_path,
            {
                "job_id": ctx.job.job_id,
                "leg": leg,
                "repeat_id": repeat_id,
                "seed_job_dir": str(seed_job_dir),
                "seed_repeat_dir": str(seed_repeat_dir),
                "seeded_files": copied,
                "seeded_at": utc_now(),
            },
        )
        copied.append(str(seed_source_path))
        return copied
    return []


def _backfill_repeat_support_itps_from_seed_source(repeat_dir: Path, repeat_top: Path) -> list[str]:
    if not repeat_top.is_file():
        return []
    seed_source_path = repeat_dir / "equilibration" / "seed_source.json"
    if not seed_source_path.is_file():
        return []
    try:
        seed_source_payload = read_json(seed_source_path)
    except OSError:
        return []
    raw_seed_repeat_dir = str(seed_source_payload.get("seed_repeat_dir") or "").strip()
    if not raw_seed_repeat_dir:
        return []
    seed_repeat_dir = Path(raw_seed_repeat_dir).expanduser()
    try:
        seed_repeat_dir = seed_repeat_dir.resolve()
    except OSError:
        pass
    if not seed_repeat_dir.is_dir():
        return []

    copied: list[str] = []
    for target_support_path in _local_topology_support_paths(repeat_top):
        if _stage_artifact_exists(target_support_path):
            continue
        source_support_path = seed_repeat_dir / target_support_path.name
        if not _stage_artifact_exists(source_support_path):
            continue
        target_support_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_support_path, target_support_path)
        copied.append(str(target_support_path))
    return copied


def _seed_sample_window_from_source(
    ctx: StageContext,
    *,
    leg: str,
    repeat_dir: Path,
    window_dir: Path,
) -> list[str]:
    if not _should_inherit_leg_from_source(ctx, leg):
        return []

    source_job_dir = _rescue_source_job_dir(ctx)
    if source_job_dir is None:
        return []

    source_window_dir = source_job_dir / "legs" / leg / repeat_dir.name / window_dir.name
    source_marker_path = window_dir / "sample_source.json"
    target_files = [
        window_dir / "topol.tpr",
        window_dir / "dhdl.xvg",
        window_dir / "md.gro",
        window_dir / "md.log",
    ]
    if all(_stage_artifact_exists(path) for path in target_files) and not _stage_artifact_exists(source_marker_path):
        return []

    source_files = {
        "topol.tpr": source_window_dir / "topol.tpr",
        "dhdl.xvg": source_window_dir / "dhdl.xvg",
        "md.gro": source_window_dir / "md.gro",
        "md.log": source_window_dir / "md.log",
    }
    if any(not _stage_artifact_exists(path) for path in source_files.values()):
        return []

    copied: list[str] = []
    for filename, source_path in source_files.items():
        target_path = window_dir / filename
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target_path)
        copied.append(str(target_path))

    write_json(
        source_marker_path,
        {
            "job_id": ctx.job.job_id,
            "leg": leg,
            "repeat_id": repeat_dir.name,
            "window_id": window_dir.name,
            "seed_job_dir": str(source_job_dir),
            "seed_window_dir": str(source_window_dir),
            "seeded_files": copied,
            "seeded_at": utc_now(),
        },
    )
    copied.append(str(source_marker_path))
    return copied


def _render_lambda_mdp(ctx: StageContext, lambda_values: list[float], window_index: int, seed: int) -> str:
    pressure_coupling = str(ctx.protocol.sampling_pressure_coupling).strip() or "C-rescale"
    return "\n".join(
        [
            "integrator              = sd",
            f"dt                      = {ctx.protocol.production_dt_ps:.3f}",
            f"nsteps                  = {_steps(ctx.protocol.production_ps, ctx.protocol.production_dt_ps)}",
            f"ld-seed                 = {seed}",
            *_render_nonbonded_block(ctx),
            "DispCorr                = EnerPres",
            "tcoupl                  = v-rescale",
            "tc-grps                 = System",
            "tau-t                   = 1.0",
            f"ref-t                   = {ctx.protocol.temperature_k}",
            f"pcoupl                  = {pressure_coupling}",
            "pcoupltype              = isotropic",
            f"tau-p                   = {ctx.protocol.sampling_pressure_tau_ps}",
            "compressibility         = 4.5e-5",
            f"ref-p                   = {ctx.protocol.pressure_bar}",
            "constraints             = h-bonds",
            "constraint-algorithm    = lincs",
            "pbc                     = xyz",
            "free-energy             = yes",
            f"init-lambda-state       = {window_index}",
            *_format_lambda_schedule_lines(lambda_values),
            "calc-lambda-neighbors   = -1",
            "sc-alpha                = 0.3",
            "sc-sigma                = 0.25",
            "sc-power                = 1",
            "sc-coul                 = yes",
            "nstdhdl                 = 100",
            "dhdl-print-energy       = total",
            f"refcoord-scaling        = {ctx.protocol.sampling_refcoord_scaling}",
            "nstenergy               = 100",
            "nstlog                  = 100",
            "nstxout-compressed      = 1000",
            "continuation            = yes",
            "gen-vel                 = no",
            "",
        ]
    )


def _render_window_relax_em_mdp(ctx: StageContext, lambda_values: list[float], window_index: int) -> str:
    return "\n".join(
        [
            "define                  = -DFLEXIBLE",
            "integrator              = steep",
            f"nsteps                  = {ctx.protocol.window_relax_em_steps}",
            "emtol                   = 1000",
            *_render_nonbonded_block(ctx),
            "constraints             = none",
            "pbc                     = xyz",
            "free-energy             = yes",
            f"init-lambda-state       = {window_index}",
            *_format_lambda_schedule_lines(lambda_values),
            "calc-lambda-neighbors   = -1",
            "sc-alpha                = 0.3",
            "sc-sigma                = 0.25",
            "sc-power                = 1",
            "sc-coul                 = yes",
            "nstdhdl                 = 100",
            "dhdl-print-energy       = total",
            "refcoord-scaling        = all",
            "nstenergy               = 100",
            "nstlog                  = 100",
            "",
        ]
    )


def _render_window_relax_md_mdp(
    ctx: StageContext,
    lambda_values: list[float],
    window_index: int,
    seed: int,
) -> str:
    return "\n".join(
        [
            "define                  = -DFLEXIBLE",
            "integrator              = sd",
            f"dt                      = {ctx.protocol.window_relax_md_dt_ps:.5f}",
            f"nsteps                  = {_steps(ctx.protocol.window_relax_md_ps, ctx.protocol.window_relax_md_dt_ps)}",
            f"ld-seed                 = {seed}",
            *_render_nonbonded_block(ctx),
            "DispCorr                = EnerPres",
            "tcoupl                  = v-rescale",
            "tc-grps                 = System",
            "tau-t                   = 1.0",
            f"ref-t                   = {ctx.protocol.temperature_k}",
            "pcoupl                  = no",
            "constraints             = none",
            "pbc                     = xyz",
            "free-energy             = yes",
            f"init-lambda-state       = {window_index}",
            *_format_lambda_schedule_lines(lambda_values),
            "calc-lambda-neighbors   = -1",
            "sc-alpha                = 0.3",
            "sc-sigma                = 0.25",
            "sc-power                = 1",
            "sc-coul                 = yes",
            "nstdhdl                 = 100",
            "dhdl-print-energy       = total",
            "refcoord-scaling        = all",
            "nstenergy               = 100",
            "nstlog                  = 100",
            "nstxout-compressed      = 1000",
            "continuation            = no",
            "gen-vel                 = no",
            "",
        ]
    )


def _sample_window_snippet(
    *,
    gmx_command: str,
    mdrun_suffix: str,
    repeat_top: Path,
    start_gro: Path,
    window_dir: Path,
    grompp_maxwarn_sampling: int,
) -> str:
    pre_relax_mdp = window_dir / "pre_relax.mdp"
    pre_relax_tpr = window_dir / "pre_relax.tpr"
    pre_relax_deffnm = window_dir / "pre_relax"
    pre_md_mdp = window_dir / "pre_md.mdp"
    pre_md_tpr = window_dir / "pre_md.tpr"
    pre_md_deffnm = window_dir / "pre_md"
    production_mdp = window_dir / "production.mdp"
    tpr_path = window_dir / "topol.tpr"
    deffnm = window_dir / "md"
    dhdl_path = window_dir / "dhdl.xvg"
    window_tag = f"{window_dir.parent.parent.name}/{window_dir.parent.name}/{window_dir.name}"
    completion_check = _shell_all_nonempty(
        [
            dhdl_path,
            deffnm.with_suffix(".gro"),
            deffnm.with_suffix(".log"),
            tpr_path,
        ]
    )
    cleanup_targets = [
        pre_relax_tpr,
        pre_relax_deffnm.with_suffix(".edr"),
        pre_relax_deffnm.with_suffix(".gro"),
        pre_relax_deffnm.with_suffix(".log"),
        pre_relax_deffnm.with_suffix(".trr"),
        pre_md_tpr,
        pre_md_deffnm.with_suffix(".cpt"),
        pre_md_deffnm.with_suffix(".edr"),
        pre_md_deffnm.with_suffix(".gro"),
        pre_md_deffnm.with_suffix(".log"),
        pre_md_deffnm.with_suffix(".trr"),
        pre_md_deffnm.with_suffix(".xtc"),
        pre_md_deffnm.with_suffix(".xvg"),
        tpr_path,
        dhdl_path,
        deffnm.with_suffix(".cpt"),
        deffnm.with_suffix(".edr"),
        deffnm.with_suffix(".gro"),
        deffnm.with_suffix(".log"),
        deffnm.with_suffix(".xtc"),
        deffnm.with_suffix(".trr"),
    ]
    cleanup_line = "rm -f " + " ".join(sh_quote(str(path)) for path in cleanup_targets)
    return "\n".join(
        [
            f"if {completion_check}; then",
            f"  echo \"[abag-rbfep] skipping completed sample window {window_tag}\"",
            "else",
            f"  echo \"[abag-rbfep] starting sample window {window_tag}\"",
            f"  {cleanup_line}",
            f"  {gmx_command} grompp -f {sh_quote(str(pre_relax_mdp))} -c {sh_quote(str(start_gro))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(pre_relax_tpr))} -maxwarn {grompp_maxwarn_sampling}",
            f"  {gmx_command} mdrun -s {sh_quote(str(pre_relax_tpr))} -deffnm {sh_quote(str(pre_relax_deffnm))}{mdrun_suffix}",
            f"  {gmx_command} grompp -f {sh_quote(str(pre_md_mdp))} -c {sh_quote(str(pre_relax_deffnm.with_suffix('.gro')))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(pre_md_tpr))} -maxwarn {grompp_maxwarn_sampling}",
            f"  {gmx_command} mdrun -s {sh_quote(str(pre_md_tpr))} -deffnm {sh_quote(str(pre_md_deffnm))}{mdrun_suffix}",
            f"  {gmx_command} grompp -f {sh_quote(str(production_mdp))} -c {sh_quote(str(pre_md_deffnm.with_suffix('.gro')))} -t {sh_quote(str(pre_md_deffnm.with_suffix('.cpt')))} -p {sh_quote(str(repeat_top))} -o {sh_quote(str(tpr_path))} -maxwarn {grompp_maxwarn_sampling}",
            f"  {gmx_command} mdrun -s {sh_quote(str(tpr_path))} -deffnm {sh_quote(str(deffnm))} -dhdl {sh_quote(str(dhdl_path))}{mdrun_suffix}",
            f"  echo \"[abag-rbfep] completed sample window {window_tag}\"",
            "fi",
        ]
    )


def _bar_repeat_snippet(
    *,
    gmx_command: str,
    output_dir: Path,
    dhdl_files: list[Path],
) -> str:
    bar_path = output_dir / "bar.xvg"
    barint_path = output_dir / "barint.xvg"
    histogram_path = output_dir / "histogram.xvg"
    completion_check = _shell_all_nonempty([bar_path, barint_path, histogram_path])
    cleanup_line = "rm -f " + " ".join(
        sh_quote(str(path))
        for path in [
            bar_path,
            barint_path,
            histogram_path,
        ]
    )
    input_files = " ".join(sh_quote(str(path)) for path in dhdl_files)
    repeat_tag = f"{output_dir.parent.parent.name}/{output_dir.parent.name}"
    return "\n".join(
        [
            f"if {completion_check}; then",
            f"  echo \"[abag-rbfep] skipping completed BAR repeat {repeat_tag}\"",
            "else",
            f"  {cleanup_line}",
            f"  {gmx_command} bar -f {input_files} -o {sh_quote(str(bar_path))} -oi {sh_quote(str(barint_path))} -oh {sh_quote(str(histogram_path))}",
            "fi",
        ]
    )


def _stage_ingest(ctx: StageContext) -> StageStatus:
    started = utc_now()
    status = StageStatus(
        stage="ingest",
        state="completed",
        message="Config manifests are present in the job workspace.",
        artifacts=[
            str(ctx.job_dir / "config" / "system.yml"),
            str(ctx.job_dir / "config" / "protocol.yml"),
            str(ctx.job_dir / "config" / "mutation_group.json"),
        ],
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_prepare(ctx: StageContext) -> StageStatus:
    started = utc_now()
    input_path = Path(ctx.system.input_structure)
    if input_path.suffix.lower() != ".pdb":
        status = StageStatus(
            stage="prepare",
            state="blocked_input",
            message="Current implementation only prepares PDB inputs. Convert CIF/mmCIF to PDB first.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    available = pdb_chain_ids(input_path)
    required = set(ctx.system.antibody_chains + ctx.system.antigen_chains)
    missing = sorted(required - available)
    if missing:
        status = StageStatus(
            stage="prepare",
            state="blocked_input",
            message=f"Input PDB is missing required chains: {', '.join(missing)}",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    artifacts = []
    qc_payload: dict[str, object] = {
        "job_id": ctx.job.job_id,
        "source_input_structure": ctx.system.input_structure,
        "legs": {},
    }
    mutated_site_keys = {(site.chain_id, site.resseq, site.icode or "") for site in ctx.mutation_group.sites}
    for leg in ("complex", "apo"):
        leg_dir = ctx.job_dir / "legs" / leg
        leg_dir.mkdir(parents=True, exist_ok=True)
        prepared_input = leg_dir / "input.pdb"
        extract_pdb_chains(input_path, prepared_input, keep_chains=_leg_mutated_chains(ctx, leg))
        # Strip any input hydrogens so pdb2gmx can run without -ignh: normal
        # residues are re-protonated from .hdb rules, while pmx hybrid residues
        # keep their explicit pmx hydrogens/dummies (-ignh would silently drop
        # them and produce heavy-atom-only, fractionally charged hybrids).
        hydrogen_strip_summary = strip_hydrogen_atoms(prepared_input, prepared_input)
        incomplete_residues = classify_incomplete_standard_residues(prepared_input)
        repair_summary: dict[str, object] = empty_repair_summary()
        if incomplete_residues:
            repair_summary = repair_incomplete_standard_residues_with_pdbfixer(prepared_input, prepared_input)
            if repair_summary.get("succeeded"):
                strip_terminal_oxygen_atoms(prepared_input, prepared_input)
            incomplete_residues = classify_incomplete_standard_residues(prepared_input)
        blocking_incomplete_residues = [item for item in incomplete_residues if item["blocking_prepare"]]
        sidechain_only_incomplete_residues = [item for item in incomplete_residues if not item["blocking_prepare"]]
        heavy_atom_clashes = find_intra_residue_heavy_atom_clashes(prepared_input)
        inter_residue_heavy_atom_clashes = find_inter_residue_heavy_atom_clashes(prepared_input)
        clash_repair_summary: dict[str, object] = empty_repair_summary()
        repairable_sidechain_clashes, blocking_heavy_atom_clashes = partition_sidechain_repairable_clashes(heavy_atom_clashes)
        deferred_sidechain_clash_residues: list[dict[str, object]] = []
        if repairable_sidechain_clashes:
            strip_sidechain_atoms_for_residues(prepared_input, prepared_input, repairable_sidechain_clashes)
            clash_repair_summary = repair_missing_atoms_with_pdbfixer(prepared_input, prepared_input)
            if clash_repair_summary.get("succeeded"):
                strip_terminal_oxygen_atoms(prepared_input, prepared_input)
            incomplete_residues = classify_incomplete_standard_residues(prepared_input)
            blocking_incomplete_residues = [item for item in incomplete_residues if item["blocking_prepare"]]
            sidechain_only_incomplete_residues = [item for item in incomplete_residues if not item["blocking_prepare"]]
            heavy_atom_clashes = find_intra_residue_heavy_atom_clashes(prepared_input)
            remaining_repairable_sidechain_clashes, blocking_heavy_atom_clashes = partition_sidechain_repairable_clashes(
                heavy_atom_clashes
            )
            deferred_candidates = [
                issue
                for issue in remaining_repairable_sidechain_clashes
                if (issue["chain_id"], issue["resseq"], issue["icode"]) not in mutated_site_keys
            ]
            if deferred_candidates and len(deferred_candidates) == len(remaining_repairable_sidechain_clashes):
                # If rebuilding reintroduces only sidechain-local clashes on non-mutated residues,
                # leave those sidechains stripped and let downstream pdb2gmx -missing rebuild them.
                strip_sidechain_atoms_for_residues(prepared_input, prepared_input, deferred_candidates)
                strip_terminal_oxygen_atoms(prepared_input, prepared_input)
                deferred_sidechain_clash_residues = deferred_candidates
                incomplete_residues = classify_incomplete_standard_residues(prepared_input)
                blocking_incomplete_residues = [item for item in incomplete_residues if item["blocking_prepare"]]
                sidechain_only_incomplete_residues = [item for item in incomplete_residues if not item["blocking_prepare"]]
                heavy_atom_clashes = find_intra_residue_heavy_atom_clashes(prepared_input)
                _, blocking_heavy_atom_clashes = partition_sidechain_repairable_clashes(heavy_atom_clashes)
        inter_residue_heavy_atom_clashes = find_inter_residue_heavy_atom_clashes(prepared_input)
        manifest = {
            "leg": leg,
            "mutated_entity_side": ctx.mutation_group.entity_side,
            "chains_retained": _leg_mutated_chains(ctx, leg),
            "input_structure": str(prepared_input),
            "source_input_structure": ctx.system.input_structure,
            "structure_source": ctx.system.structure_source,
            "incomplete_standard_residue_count": len(incomplete_residues),
            "blocking_incomplete_standard_residue_count": len(blocking_incomplete_residues),
            "sidechain_only_incomplete_standard_residue_count": len(sidechain_only_incomplete_residues),
            "intra_residue_heavy_atom_clash_count": len(heavy_atom_clashes),
            "inter_residue_heavy_atom_clash_count": len(inter_residue_heavy_atom_clashes),
        }
        write_json(leg_dir / "manifest.json", manifest)
        qc_payload["legs"][leg] = {
            "chains_retained": _leg_mutated_chains(ctx, leg),
            "incomplete_standard_residues": incomplete_residues,
            "blocking_incomplete_standard_residues": blocking_incomplete_residues,
            "sidechain_only_incomplete_standard_residues": sidechain_only_incomplete_residues,
            "intra_residue_heavy_atom_clashes": heavy_atom_clashes,
            "blocking_intra_residue_heavy_atom_clashes": blocking_heavy_atom_clashes,
            "inter_residue_heavy_atom_clashes": inter_residue_heavy_atom_clashes,
            "repair_summary": repair_summary,
            "hydrogen_strip_summary": hydrogen_strip_summary,
            "repairable_sidechain_clashes": repairable_sidechain_clashes,
            "clash_repair_summary": clash_repair_summary,
            "deferred_sidechain_clash_residues": deferred_sidechain_clash_residues,
        }
        artifacts.append(str(prepared_input))
        artifacts.append(str(leg_dir / "manifest.json"))
    qc_path = ctx.job_dir / "artifacts" / "prepare_qc.json"
    write_json(qc_path, qc_payload)
    artifacts.append(str(qc_path))
    blocking_legs = {
        leg: payload["blocking_incomplete_standard_residues"]
        for leg, payload in qc_payload["legs"].items()
        if payload["blocking_incomplete_standard_residues"]
    }
    if blocking_legs:
        samples = []
        for leg, issues in blocking_legs.items():
            for issue in issues[:2]:
                insertion = issue["icode"] or ""
                missing_atoms = ",".join(issue["missing_backbone_atoms"])
                samples.append(
                    f"{leg}:{issue['chain_id']}{issue['resseq']}{insertion} {issue['resname']} missing backbone atoms {missing_atoms}"
                )
        extra_count = sum(len(issues) for issues in blocking_legs.values()) - len(samples)
        suffix = f" (+{extra_count} more)" if extra_count > 0 else ""
        status = StageStatus(
            stage="prepare",
            state="blocked_input",
            message="Prepared PDB contains backbone-incomplete standard residues: " + "; ".join(samples) + suffix,
            artifacts=artifacts,
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    blocking_clash_legs = {
        leg: payload["blocking_intra_residue_heavy_atom_clashes"]
        for leg, payload in qc_payload["legs"].items()
        if payload["blocking_intra_residue_heavy_atom_clashes"]
    }
    if blocking_clash_legs:
        samples = []
        for leg, issues in blocking_clash_legs.items():
            for issue in issues[:2]:
                insertion = issue["icode"] or ""
                first_clash = issue["clashes"][0]
                samples.append(
                    f"{leg}:{issue['chain_id']}{issue['resseq']}{insertion} {issue['resname']} "
                    f"{first_clash['atom_a']}-{first_clash['atom_b']} {first_clash['distance_angstrom']:.3f} A"
                )
        extra_count = sum(len(issues) for issues in blocking_clash_legs.values()) - len(samples)
        suffix = f" (+{extra_count} more)" if extra_count > 0 else ""
        status = StageStatus(
            stage="prepare",
            state="blocked_input",
            message="Prepared PDB contains same-residue heavy-atom clashes: " + "; ".join(samples) + suffix,
            artifacts=artifacts,
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    blocking_inter_residue_clash_legs = {
        leg: payload["inter_residue_heavy_atom_clashes"]
        for leg, payload in qc_payload["legs"].items()
        if payload["inter_residue_heavy_atom_clashes"]
    }
    if blocking_inter_residue_clash_legs:
        samples = []
        for leg, issues in blocking_inter_residue_clash_legs.items():
            for issue in issues[:2]:
                insertion = issue["icode"] or ""
                partner_insertion = issue["partner_icode"] or ""
                first_clash = issue["clashes"][0]
                samples.append(
                    f"{leg}:{issue['chain_id']}{issue['resseq']}{insertion} {issue['resname']} vs "
                    f"{issue['partner_chain_id']}{issue['partner_resseq']}{partner_insertion} {issue['partner_resname']} "
                    f"{first_clash['atom_a']}-{first_clash['atom_b']} {first_clash['distance_angstrom']:.3f} A"
                )
        extra_count = sum(len(issues) for issues in blocking_inter_residue_clash_legs.values()) - len(samples)
        suffix = f" (+{extra_count} more)" if extra_count > 0 else ""
        status = StageStatus(
            stage="prepare",
            state="blocked_input",
            message="Prepared PDB contains impossible inter-residue heavy-atom clashes: " + "; ".join(samples) + suffix,
            artifacts=artifacts,
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    sidechain_warning_count = sum(
        len(payload["sidechain_only_incomplete_standard_residues"]) for payload in qc_payload["legs"].values()
    )
    message = "Leg-specific PDB inputs, manifests, and structure QC generated."
    if sidechain_warning_count:
        message += (
            f" Detected {sidechain_warning_count} sidechain-incomplete standard residues;"
            " stripped their remaining sidechain atoms so downstream pdb2gmx -missing can attempt reconstruction."
        )
    if any(payload["repair_summary"].get("succeeded") for payload in qc_payload["legs"].values()):
        message += " Applied PDBFixer atom repair for incomplete residues."
    if any(payload["clash_repair_summary"].get("succeeded") for payload in qc_payload["legs"].values()):
        message += " Repaired sidechain-involving same-residue clashes by stripping and rebuilding affected sidechains."
    deferred_sidechain_clash_count = sum(
        len(payload.get("deferred_sidechain_clash_residues", [])) for payload in qc_payload["legs"].values()
    )
    if deferred_sidechain_clash_count:
        message += (
            f" Deferred {deferred_sidechain_clash_count} repairable clash-sidechains as stripped residues"
            " so downstream pdb2gmx -missing can rebuild them."
        )
    status = StageStatus(
        stage="prepare",
        state="completed",
        message=message,
        artifacts=artifacts,
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_mutate(ctx: StageContext) -> StageStatus:
    started = utc_now()
    if any(site.icode for site in ctx.mutation_group.sites):
        status = StageStatus(
            stage="mutate",
            state="blocked_input",
            message="Insertion codes are tracked in manifests but not yet mapped to pmx integer residue IDs.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    pmx_command_argv = _resolve_pmx_argv(ctx)
    pmx_command = " ".join(sh_quote(token) for token in pmx_command_argv) if pmx_command_argv is not None else None
    if pmx_command is None and ctx.runner.execute:
        status = StageStatus(
            stage="mutate",
            state="blocked_external",
            message="pmx is not available. Create the project .venv or install pmx in PATH.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    if pmx_command is None:
        pmx_command = sh_quote(ctx.protocol.pmx_bin)

    gmx_command = _resolve_gmx_command(ctx)
    project_python_command = _resolve_project_python_command(ctx)
    env = _stage_env(ctx)
    if ctx.runner.execute and env is None:
        status = StageStatus(
            stage="mutate",
            state="blocked_external",
            message="Could not determine a usable GROMACS topology library for the pmx mutation force field.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    stage_commands: list[str] = []
    artifacts: list[str] = []
    leg_qc_paths: dict[str, Path] = {}
    leg_restore_paths: dict[str, Path] = {}
    leg_pdbfixer_paths: dict[str, Path] = {}
    leg_processed_gro_qc_paths: dict[str, Path] = {}
    mutated_chain_ids = sorted({str(site.chain_id).strip() for site in ctx.mutation_group.sites if str(site.chain_id).strip()})
    gmxlib = env.get("GMXLIB") if env else None
    if gmxlib:
        artifacts.append(gmxlib)
    for leg in ("complex", "apo"):
        prepared_input = ctx.job_dir / "legs" / leg / "input.pdb"
        if not prepared_input.is_file():
            status = StageStatus(
                stage="mutate",
                state="blocked_input",
                message=f"Prepared leg input is missing: {prepared_input}",
                started_at=started,
                completed_at=utc_now(),
            )
            return _write_stage_status(ctx.job_dir, status)
        leg_dir = ctx.job_dir / "legs" / leg / "pmx"
        leg_dir.mkdir(parents=True, exist_ok=True)
        script_path = leg_dir / "mutations.txt"
        script_path.write_text("\n".join(mutation_script_lines(ctx.mutation_group.sites)) + "\n", encoding="utf-8")
        mutant_pdb = leg_dir / "mutant.pdb"
        processed_gro = leg_dir / "processed.gro"
        topology_top = leg_dir / "topol.top"
        pmx_top = leg_dir / "pmxtop.top"
        mutant_geometry_qc = leg_dir / "mutant_geometry_qc.json"
        mutant_standard_residue_repair = leg_dir / "mutant_standard_residue_repair.json"
        mutant_pdbfixer_repair = leg_dir / "mutant_pdbfixer_repair.json"
        processed_gro_qc = leg_dir / "processed_gro_qc.json"
        leg_qc_paths[leg] = mutant_geometry_qc
        leg_restore_paths[leg] = mutant_standard_residue_repair
        leg_pdbfixer_paths[leg] = mutant_pdbfixer_repair
        leg_processed_gro_qc_paths[leg] = processed_gro_qc
        artifacts.extend(
            [
                str(script_path),
                str(mutant_pdb),
                str(processed_gro),
                str(topology_top),
                str(pmx_top),
                str(mutant_standard_residue_repair),
                str(mutant_pdbfixer_repair),
                str(processed_gro_qc),
            ]
        )
        strip_mutant_terminal_oxygen_command = (
            f"{project_python_command} -c "
            + sh_quote(
                "from pathlib import Path; "
                "from abag_rbfe.structure import strip_terminal_oxygen_atoms; "
                "strip_terminal_oxygen_atoms(Path('mutant.pdb'), Path('mutant.pdb'))"
            )
        )
        mutant_standard_residue_restore_command = _mutant_standard_residue_restore_command(
            project_python_command,
            template_path=prepared_input,
            target_path=Path(mutant_pdb.name),
            summary_path=Path(mutant_standard_residue_repair.name),
        )
        mutant_pdbfixer_repair_command = _mutant_pdbfixer_sidechain_repair_command(
            project_python_command,
            target_path=Path(mutant_pdb.name),
            summary_path=Path(mutant_pdbfixer_repair.name),
        )
        mutant_inter_residue_qc_command = (
            f"{project_python_command} -c "
            + sh_quote(
                "from pathlib import Path; "
                "from abag_rbfe.structure import write_inter_residue_heavy_atom_clash_report; "
                "import sys; "
                f"payload = write_inter_residue_heavy_atom_clash_report(Path('mutant.pdb'), Path({mutant_geometry_qc.name!r}), reference_path=Path({str(prepared_input)!r})); "
                "sys.exit(2 if payload.get('blocking_inter_residue_heavy_atom_clashes') else 0)"
            )
        )
        processed_gro_validation_command = _processed_gro_validation_command(
            project_python_command,
            gro_path=Path(processed_gro.name),
            summary_path=Path(processed_gro_qc.name),
        )
        hybrid_topology_generation_command = _hybrid_topology_generation_command(
            project_python_command,
            topology_path=Path(topology_top.name),
            output_path=Path(pmx_top.name),
            force_field=ctx.protocol.force_field,
            mutated_chain_ids=mutated_chain_ids,
            pmx_command_argv=pmx_command_argv or [ctx.protocol.pmx_bin],
            restore_summary_path=Path(mutant_standard_residue_repair.name),
            pdbfixer_summary_path=Path(mutant_pdbfixer_repair.name),
        )
        inner_commands = " && ".join(
            [
                f"cd {sh_quote(str(leg_dir))}",
                f"{pmx_command} mutate -f {sh_quote(str(prepared_input))} -o {sh_quote(str(mutant_pdb))} -ff {ctx.protocol.force_field} --script {sh_quote(str(script_path))} --keep_resid",
                strip_mutant_terminal_oxygen_command,
                mutant_standard_residue_restore_command,
                mutant_pdbfixer_repair_command,
                mutant_inter_residue_qc_command,
                # No -ignh here: hydrogens were already stripped at prepare time,
                # and pmx hybrid residues must keep their explicit pmx hydrogens
                # (hybrid residues have no .hdb entries to rebuild them from).
                f"{gmx_command} pdb2gmx -missing -f {sh_quote(str(mutant_pdb))} -o {sh_quote(str(processed_gro))} -p {sh_quote(str(topology_top))} -ff {ctx.protocol.force_field} -water {ctx.protocol.water_model}",
                processed_gro_validation_command,
                hybrid_topology_generation_command,
            ]
        )
        stage_commands.append(f"bash -lc {sh_quote(inner_commands)}")
    script_file = ctx.job_dir / "artifacts" / "commands" / "mutate.sh"
    _write_external_stage_running(
        ctx.job_dir,
        stage="mutate",
        script_file=script_file,
        commands=stage_commands,
        artifacts=artifacts,
        started_at=started,
    )
    outcome = ctx.runner.run_script(script_file, stage_commands, ctx.job_dir, env=env)
    qc_artifacts: list[str] = []
    qc_payload: dict[str, object] = {
        "job_id": ctx.job.job_id,
        "legs": {},
    }
    qc_path = ctx.job_dir / "artifacts" / "mutate_qc.json"
    missing_mutant_pdbs: list[str] = []
    for leg in ("complex", "apo"):
        leg_payload, leg_artifacts, missing_mutant = _collect_mutate_leg_qc(
            ctx.job_dir,
            leg,
            geometry_qc_path=leg_qc_paths[leg],
            restore_summary_path=leg_restore_paths[leg],
            pdbfixer_summary_path=leg_pdbfixer_paths[leg],
            processed_gro_qc_path=leg_processed_gro_qc_paths[leg],
        )
        qc_payload["legs"][leg] = leg_payload
        qc_artifacts.extend(leg_artifacts)
        if missing_mutant and outcome.state == "completed":
            missing_mutant_pdbs.append(str(ctx.job_dir / "legs" / leg / "pmx" / "mutant.pdb"))
    blocking_inter_residue_clash_legs = _blocking_mutate_qc_legs_from_payload(qc_payload)
    if blocking_inter_residue_clash_legs and ctx.runner.execute:
        repaired_any = False
        for leg, issues in list(blocking_inter_residue_clash_legs.items()):
            repaired_issues, repair_summary = _attempt_repair_mutant_sidechain_clashes(ctx, leg=leg, issues=issues)
            qc_payload["legs"][leg]["inter_residue_heavy_atom_clashes"] = repaired_issues
            qc_payload["legs"][leg]["auto_repair_summary"] = repair_summary
            if repaired_issues != issues:
                repaired_any = True
        if repaired_any:
            qc_artifacts = []
            missing_mutant_pdbs = []
            for leg in ("complex", "apo"):
                leg_payload, leg_artifacts, missing_mutant = _collect_mutate_leg_qc(
                    ctx.job_dir,
                    leg,
                    geometry_qc_path=leg_qc_paths[leg],
                    restore_summary_path=leg_restore_paths[leg],
                    pdbfixer_summary_path=leg_pdbfixer_paths[leg],
                    processed_gro_qc_path=leg_processed_gro_qc_paths[leg],
                )
                auto_repair_summary = qc_payload["legs"].get(leg, {}).get("auto_repair_summary")
                if auto_repair_summary:
                    leg_payload["auto_repair_summary"] = auto_repair_summary
                qc_payload["legs"][leg] = leg_payload
                qc_artifacts.extend(leg_artifacts)
                if missing_mutant and outcome.state == "completed":
                    missing_mutant_pdbs.append(str(ctx.job_dir / "legs" / leg / "pmx" / "mutant.pdb"))
            blocking_inter_residue_clash_legs = _blocking_mutate_qc_legs_from_payload(qc_payload)
            if not blocking_inter_residue_clash_legs:
                outcome = CommandOutcome(
                    "completed",
                    "Mutant sidechain clashes were auto-repaired and mutate outputs were regenerated.",
                )
    write_json(qc_path, qc_payload)
    qc_artifacts.append(str(qc_path))
    blocking_incomplete_standard_residue_legs = {
        str(leg): [
            issue
            for issue in leg_payload.get("incomplete_standard_residues", [])
            if isinstance(issue, dict)
        ]
        for leg, leg_payload in qc_payload["legs"].items()
        if isinstance(leg_payload, dict) and leg_payload.get("incomplete_standard_residues")
    }
    invalid_processed_gro_legs = {
        str(leg): dict(leg_payload.get("processed_gro_qc", {}))
        for leg, leg_payload in qc_payload["legs"].items()
        if isinstance(leg_payload, dict)
        and isinstance(leg_payload.get("processed_gro_qc"), dict)
        and leg_payload["processed_gro_qc"]
        and not bool(leg_payload["processed_gro_qc"].get("valid"))
    }
    invalid_hybrid_integrity_legs = {
        str(leg): list(leg_payload["hybrid_integrity"].get("issues", []))
        for leg, leg_payload in qc_payload["legs"].items()
        if isinstance(leg_payload, dict)
        and isinstance(leg_payload.get("hybrid_integrity"), dict)
        and leg_payload["hybrid_integrity"].get("checked")
        and not leg_payload["hybrid_integrity"].get("ok")
    }
    if invalid_hybrid_integrity_legs:
        samples = [
            f"{leg}:{issue}"
            for leg, issues in invalid_hybrid_integrity_legs.items()
            for issue in issues[:2]
        ]
        status = StageStatus(
            stage="mutate",
            state="blocked_input",
            message="Hybrid residue integrity check failed: " + "; ".join(samples),
            commands=stage_commands,
            artifacts=_external_stage_artifacts(script_file, artifacts + qc_artifacts),
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    if blocking_inter_residue_clash_legs:
        samples = []
        for leg, issues in blocking_inter_residue_clash_legs.items():
            for issue in issues[:2]:
                insertion = issue["icode"] or ""
                partner_insertion = issue["partner_icode"] or ""
                first_clash = issue["clashes"][0]
                samples.append(
                    f"{leg}:{issue['chain_id']}{issue['resseq']}{insertion} {issue['resname']} vs "
                    f"{issue['partner_chain_id']}{issue['partner_resseq']}{partner_insertion} {issue['partner_resname']} "
                    f"{first_clash['atom_a']}-{first_clash['atom_b']} {first_clash['distance_angstrom']:.3f} A"
                )
        extra_count = sum(len(issues) for issues in blocking_inter_residue_clash_legs.values()) - len(samples)
        suffix = f" (+{extra_count} more)" if extra_count > 0 else ""
        status = StageStatus(
            stage="mutate",
            state="blocked_input",
            message="Mutated structure contains impossible inter-residue heavy-atom clashes: "
            + "; ".join(samples)
            + suffix,
            commands=stage_commands,
            artifacts=_external_stage_artifacts(script_file, artifacts + qc_artifacts),
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    if outcome.state != "completed" and blocking_incomplete_standard_residue_legs:
        status = StageStatus(
            stage="mutate",
            state="blocked_input",
            message=_incomplete_standard_residue_message(blocking_incomplete_standard_residue_legs),
            commands=stage_commands,
            artifacts=_external_stage_artifacts(script_file, artifacts + qc_artifacts),
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    if invalid_processed_gro_legs:
        status = StageStatus(
            stage="mutate",
            state="blocked_input",
            message=_invalid_processed_gro_message(invalid_processed_gro_legs),
            commands=stage_commands,
            artifacts=_external_stage_artifacts(script_file, artifacts + qc_artifacts),
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    if outcome.state == "completed" and missing_mutant_pdbs:
        status = StageStatus(
            stage="mutate",
            state="blocked_input",
            message=(
                "Mutate completed but expected mutant coordinate files are missing: "
                + ", ".join(missing_mutant_pdbs)
            ),
            commands=stage_commands,
            artifacts=_external_stage_artifacts(script_file, artifacts + qc_artifacts),
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)
    status = StageStatus(
        stage="mutate",
        state=outcome.state,
        message=outcome.message + (" Mutated geometry QC passed." if outcome.state == "completed" else ""),
        commands=stage_commands,
        artifacts=_external_stage_artifacts(script_file, artifacts + qc_artifacts),
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _mutation_site_keys(ctx: StageContext) -> set[tuple[str, int, str]]:
    return {
        (
            str(site.chain_id).strip(),
            int(site.resseq),
            str(site.icode or "").strip().upper(),
        )
        for site in ctx.mutation_group.sites
    }


def _issue_residue_key(issue: dict[str, object], *, partner: bool = False) -> tuple[str, int, str]:
    prefix = "partner_" if partner else ""
    return (
        str(issue.get(f"{prefix}chain_id") or "").strip(),
        int(issue.get(f"{prefix}resseq") or 0),
        str(issue.get(f"{prefix}icode") or "").strip().upper(),
    )


def _issue_residue_payload(issue: dict[str, object], *, partner: bool = False) -> dict[str, object]:
    prefix = "partner_" if partner else ""
    return {
        "chain_id": str(issue.get(f"{prefix}chain_id") or "").strip(),
        "resseq": int(issue.get(f"{prefix}resseq") or 0),
        "icode": str(issue.get(f"{prefix}icode") or "").strip().upper(),
    }


def _read_json_dict_if_exists(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        payload = read_json(path)
    except OSError:
        return {}
    return payload if isinstance(payload, dict) else {}


def _path_mtime_ns(path: Path) -> int | None:
    try:
        return path.stat().st_mtime_ns
    except OSError:
        return None


def _should_refresh_geometry_qc(
    geometry_qc_path: Path,
    geometry_qc: dict[str, object],
    *,
    mutant_pdb: Path,
    reference_path: Path | None,
) -> bool:
    if not geometry_qc_path.is_file() or not geometry_qc:
        return True
    if "inter_residue_heavy_atom_clashes" not in geometry_qc:
        return True

    expected_reference = str(reference_path) if reference_path is not None else None
    recorded_reference = geometry_qc.get("reference_structure")
    if recorded_reference != expected_reference:
        return True

    geometry_mtime = _path_mtime_ns(geometry_qc_path)
    mutant_mtime = _path_mtime_ns(mutant_pdb)
    if geometry_mtime is None or mutant_mtime is None:
        return True
    if mutant_mtime > geometry_mtime:
        return True
    if reference_path is not None:
        reference_mtime = _path_mtime_ns(reference_path)
        if reference_mtime is not None and reference_mtime > geometry_mtime:
            return True
    return False


def _mutant_standard_residue_restore_command(
    project_python_command: str,
    *,
    template_path: Path,
    target_path: Path,
    summary_path: Path,
    exclude_residues: list[dict[str, object]] | None = None,
) -> str:
    exclude_payload = json.dumps(exclude_residues or [], ensure_ascii=False, separators=(",", ":"))
    return (
        f"{project_python_command} -c "
        + sh_quote(
            "import json; "
            "from pathlib import Path; "
            "from abag_rbfe.structure import restore_incomplete_standard_residues_from_template; "
            f"summary = restore_incomplete_standard_residues_from_template(Path({str(template_path)!r}), Path({str(target_path)!r}), Path({str(target_path)!r}), exclude_residues=json.loads({exclude_payload!r})); "
            f"Path({str(summary_path)!r}).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\\n', encoding='utf-8')"
        )
    )


def _mutant_pdbfixer_sidechain_repair_command(
    project_python_command: str,
    *,
    target_path: Path,
    summary_path: Path,
) -> str:
    return (
        f"{project_python_command} -c "
        + sh_quote(
            "import json; "
            "from pathlib import Path; "
            "from abag_rbfe.structure import repair_sidechain_only_incomplete_residues_with_pdbfixer; "
            f"summary = repair_sidechain_only_incomplete_residues_with_pdbfixer(Path({str(target_path)!r}), Path({str(target_path)!r})); "
            f"Path({str(summary_path)!r}).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\\n', encoding='utf-8')"
        )
    )


def _processed_gro_validation_command(
    project_python_command: str,
    *,
    gro_path: Path,
    summary_path: Path,
) -> str:
    return (
        f"{project_python_command} -c "
        + sh_quote(
            "import json, sys; "
            "from pathlib import Path; "
            "from abag_rbfe.gmx import inspect_gro_file; "
            f"summary = inspect_gro_file(Path({str(gro_path)!r})); "
            f"Path({str(summary_path)!r}).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\\n', encoding='utf-8'); "
            "sys.exit(0 if summary.get('valid') else 2)"
        )
    )


def _hybrid_topology_generation_command(
    project_python_command: str,
    *,
    topology_path: Path,
    output_path: Path,
    force_field: str,
    mutated_chain_ids: list[str],
    pmx_command_argv: list[str],
    restore_summary_path: Path,
    pdbfixer_summary_path: Path,
) -> str:
    return (
        f"{project_python_command} -c "
        + sh_quote(
            "from pathlib import Path; "
            "from abag_rbfe.gmx import generate_hybrid_topology; "
            f"generate_hybrid_topology(Path({str(topology_path)!r}), Path({str(output_path)!r}), {force_field!r}, {mutated_chain_ids!r}, {pmx_command_argv!r}, "
            f"restore_summary_path=Path({str(restore_summary_path)!r}), pdbfixer_summary_path=Path({str(pdbfixer_summary_path)!r}), allow_reuse_existing=True)"
        )
    )


def _materialize_staged_equilibration_restraints_command(
    project_python_command: str,
    *,
    repeat_topology_path: Path,
    summary_path: Path,
    heavy_force_constant: float,
    backbone_force_constant: float,
) -> str:
    return (
        f"{project_python_command} -c "
        + sh_quote(
            "import json; "
            "from pathlib import Path; "
            "from abag_rbfe.gmx import materialize_staged_equilibration_restraints; "
            f"summary = materialize_staged_equilibration_restraints(Path({str(repeat_topology_path)!r}), "
            f"heavy_force_constant={heavy_force_constant!r}, backbone_force_constant={backbone_force_constant!r}); "
            f"Path({str(summary_path)!r}).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\\n', encoding='utf-8')"
        )
    )


def _collect_mutate_leg_qc(
    job_dir: Path,
    leg: str,
    *,
    geometry_qc_path: Path,
    restore_summary_path: Path,
    pdbfixer_summary_path: Path,
    processed_gro_qc_path: Path,
) -> tuple[dict[str, object], list[str], bool]:
    pmx_dir = job_dir / "legs" / leg / "pmx"
    mutant_pdb = pmx_dir / "mutant.pdb"
    processed_gro = pmx_dir / "processed.gro"
    reference_path = job_dir / "legs" / leg / "input.pdb"

    artifacts: list[str] = []
    for candidate in (geometry_qc_path, restore_summary_path, pdbfixer_summary_path, processed_gro_qc_path):
        if candidate.is_file():
            artifacts.append(str(candidate))

    geometry_qc = _read_json_dict_if_exists(geometry_qc_path)
    if mutant_pdb.is_file() and _should_refresh_geometry_qc(
        geometry_qc_path,
        geometry_qc,
        mutant_pdb=mutant_pdb,
        reference_path=reference_path if reference_path.is_file() else None,
    ):
        qc_payload = write_inter_residue_heavy_atom_clash_report(
            mutant_pdb,
            geometry_qc_path,
            reference_path=reference_path if reference_path.is_file() else None,
        )
        inter_residue_heavy_atom_clashes = list(qc_payload.get("inter_residue_heavy_atom_clashes", []))
        geometry_qc = qc_payload
        artifacts.append(str(geometry_qc_path))
    elif geometry_qc:
        inter_residue_heavy_atom_clashes = list(geometry_qc.get("inter_residue_heavy_atom_clashes", []))
    else:
        inter_residue_heavy_atom_clashes = []

    mutant_standard_residue_repair = _read_json_dict_if_exists(restore_summary_path)
    mutant_pdbfixer_repair = _read_json_dict_if_exists(pdbfixer_summary_path)
    processed_gro_qc = _read_json_dict_if_exists(processed_gro_qc_path)
    if processed_gro.is_file():
        processed_gro_qc = inspect_gro_file(processed_gro)
        write_json(processed_gro_qc_path, processed_gro_qc)
        if str(processed_gro_qc_path) not in artifacts:
            artifacts.append(str(processed_gro_qc_path))

    incomplete_standard_residues = classify_incomplete_standard_residues(mutant_pdb) if mutant_pdb.is_file() else []
    hybrid_integrity = validate_hybrid_topology_integrity(pmx_dir)
    return (
        {
            "mutant_pdb": str(mutant_pdb),
            "processed_gro": str(processed_gro),
            "inter_residue_heavy_atom_clashes": inter_residue_heavy_atom_clashes,
            "incomplete_standard_residues": incomplete_standard_residues,
            "mutant_standard_residue_repair": mutant_standard_residue_repair,
            "mutant_pdbfixer_repair": mutant_pdbfixer_repair,
            "processed_gro_qc": processed_gro_qc,
            "hybrid_integrity": hybrid_integrity,
        },
        artifacts,
        not mutant_pdb.is_file(),
    )


def _format_incomplete_standard_residue_issue(issue: dict[str, object]) -> str:
    insertion = str(issue.get("icode") or "").strip()
    missing_atoms = ",".join(str(atom) for atom in issue.get("missing_atoms", []))
    return f"{issue.get('chain_id', '?')}{issue.get('resseq', '?')}{insertion} {issue.get('resname', '?')} missing {missing_atoms}"


def _incomplete_standard_residue_message(
    legs: dict[str, list[dict[str, object]]],
) -> str:
    samples: list[str] = []
    for leg, issues in legs.items():
        for issue in issues[:2]:
            samples.append(f"{leg}:{_format_incomplete_standard_residue_issue(issue)}")
    extra_count = sum(len(issues) for issues in legs.values()) - len(samples)
    suffix = f" (+{extra_count} more)" if extra_count > 0 else ""
    return "Mutant PDB retains incomplete standard residues after mutation: " + "; ".join(samples) + suffix


def _invalid_processed_gro_message(legs: dict[str, dict[str, object]]) -> str:
    samples: list[str] = []
    for leg, payload in legs.items():
        reason = str(payload.get("reason") or "invalid_processed_gro")
        detail = reason
        line_number = payload.get("line_number")
        if line_number is not None:
            detail += f" at line {line_number}"
        max_abs_coordinate = payload.get("max_abs_coordinate_nm")
        if isinstance(max_abs_coordinate, (int, float)) and max_abs_coordinate > 0:
            detail += f" (max |coord| {max_abs_coordinate:.3f} nm)"
        residue_number = payload.get("residue_number")
        residue_name = str(payload.get("residue_name") or "").strip()
        atom_name = str(payload.get("atom_name") or "").strip()
        nearest_heavy_atom = str(payload.get("nearest_heavy_atom") or "").strip()
        nearest_heavy_distance_nm = payload.get("nearest_heavy_distance_nm")
        if residue_number is not None and residue_name and atom_name:
            residue_detail = f"{residue_name}{residue_number} {atom_name}"
            if nearest_heavy_atom and isinstance(nearest_heavy_distance_nm, (int, float)):
                residue_detail += f" vs {nearest_heavy_atom} {nearest_heavy_distance_nm:.3f} nm"
            detail += f" ({residue_detail})"
        samples.append(f"{leg}:{detail}")
    return "Mutate generated an invalid processed.gro: " + "; ".join(samples)


def _mutate_repair_commands(
    ctx: StageContext,
    *,
    prepared_input: Path,
    leg_dir: Path,
    mutant_pdb: Path,
    processed_gro: Path,
    topology_top: Path,
    pmx_top: Path,
    exclude_residues: list[dict[str, object]],
) -> list[str] | None:
    pmx_command_argv = _resolve_pmx_argv(ctx)
    pmx_command = " ".join(sh_quote(token) for token in pmx_command_argv) if pmx_command_argv is not None else None
    if pmx_command is None and ctx.runner.execute:
        return None
    if pmx_command is None:
        pmx_command = sh_quote(ctx.protocol.pmx_bin)
    gmx_command = _resolve_gmx_command(ctx)
    project_python_command = _resolve_project_python_command(ctx)
    restore_summary_path = leg_dir / "mutant_standard_residue_repair.json"
    pdbfixer_summary_path = leg_dir / "mutant_pdbfixer_repair.json"
    processed_gro_qc_path = leg_dir / "processed_gro_qc.json"
    restore_command = _mutant_standard_residue_restore_command(
        project_python_command,
        template_path=prepared_input,
        target_path=Path(mutant_pdb.name),
        summary_path=Path(restore_summary_path.name),
        exclude_residues=exclude_residues,
    )
    pdbfixer_repair_command = _mutant_pdbfixer_sidechain_repair_command(
        project_python_command,
        target_path=Path(mutant_pdb.name),
        summary_path=Path(pdbfixer_summary_path.name),
    )
    processed_gro_validation_command = _processed_gro_validation_command(
        project_python_command,
        gro_path=Path(processed_gro.name),
        summary_path=Path(processed_gro_qc_path.name),
    )
    hybrid_topology_generation_command = _hybrid_topology_generation_command(
        project_python_command,
        topology_path=Path(topology_top.name),
        output_path=Path(pmx_top.name),
        force_field=ctx.protocol.force_field,
        mutated_chain_ids=sorted({str(site.chain_id).strip() for site in ctx.mutation_group.sites if str(site.chain_id).strip()}),
        pmx_command_argv=pmx_command_argv or [ctx.protocol.pmx_bin],
        restore_summary_path=Path(restore_summary_path.name),
        pdbfixer_summary_path=Path(pdbfixer_summary_path.name),
    )
    inner_commands = " && ".join(
        [
            f"cd {sh_quote(str(leg_dir))}",
            restore_command,
            pdbfixer_repair_command,
            # See the mutate-stage note: no -ignh so hybrid hydrogens survive.
            f"{gmx_command} pdb2gmx -missing -f {sh_quote(str(mutant_pdb))} -o {sh_quote(str(processed_gro))} -p {sh_quote(str(topology_top))} -ff {ctx.protocol.force_field} -water {ctx.protocol.water_model}",
            processed_gro_validation_command,
            hybrid_topology_generation_command,
        ]
    )
    return [f"bash -lc {sh_quote(inner_commands)}"]


def _attempt_repair_mutant_sidechain_clashes(
    ctx: StageContext,
    *,
    leg: str,
    issues: list[dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    summary: dict[str, object] = {
        "attempted": False,
        "succeeded": False,
        "repair_target_residues": [],
    }
    repairable_issues, blocking_issues = partition_inter_residue_sidechain_repairable_clashes(issues)
    mutated_site_keys = _mutation_site_keys(ctx)
    targets: list[dict[str, object]] = []
    seen_targets: set[tuple[str, int, str]] = set()
    for issue in repairable_issues:
        residue_keys = (
            _issue_residue_key(issue, partner=False),
            _issue_residue_key(issue, partner=True),
        )
        if any(key in mutated_site_keys for key in residue_keys):
            blocking_issues.append(issue)
            continue
        for partner in (False, True):
            key = _issue_residue_key(issue, partner=partner)
            if key in seen_targets:
                continue
            seen_targets.add(key)
            targets.append(_issue_residue_payload(issue, partner=partner))

    summary["repair_target_residues"] = targets
    if blocking_issues or not targets or not ctx.runner.execute:
        return issues, summary

    mutant_pdb = ctx.job_dir / "legs" / leg / "pmx" / "mutant.pdb"
    processed_gro = ctx.job_dir / "legs" / leg / "pmx" / "processed.gro"
    topology_top = ctx.job_dir / "legs" / leg / "pmx" / "topol.top"
    pmx_top = ctx.job_dir / "legs" / leg / "pmx" / "pmxtop.top"
    leg_dir = mutant_pdb.parent
    prepared_input = ctx.job_dir / "legs" / leg / "input.pdb"
    script_file = ctx.job_dir / "artifacts" / "commands" / f"mutate_repair_{leg}.sh"
    summary["attempted"] = True
    summary["repair_script"] = str(script_file)
    summary["repair_log"] = str(script_file.with_suffix(".log"))

    strip_sidechain_atoms_for_residues(mutant_pdb, mutant_pdb, targets)
    strip_terminal_oxygen_atoms(mutant_pdb, mutant_pdb)
    stripped_clashes = find_inter_residue_heavy_atom_clashes(mutant_pdb)
    summary["post_strip_clash_count"] = len(stripped_clashes)
    if stripped_clashes:
        return stripped_clashes, summary

    env = _stage_env(ctx)
    if env is None:
        summary["regeneration_state"] = "blocked_external"
        summary["regeneration_message"] = "Could not determine a usable GROMACS topology library for clash repair."
        return issues, summary
    commands = _mutate_repair_commands(
        ctx,
        prepared_input=prepared_input,
        leg_dir=leg_dir,
        mutant_pdb=mutant_pdb,
        processed_gro=processed_gro,
        topology_top=topology_top,
        pmx_top=pmx_top,
        exclude_residues=targets,
    )
    if commands is None:
        summary["regeneration_state"] = "blocked_external"
        summary["regeneration_message"] = "pmx is not available for mutant clash repair."
        return issues, summary

    outcome = ctx.runner.run_script(script_file, commands, ctx.job_dir, env=env)
    summary["regeneration_state"] = outcome.state
    summary["regeneration_message"] = outcome.message
    if outcome.state != "completed":
        return issues, summary

    qc_payload = write_inter_residue_heavy_atom_clash_report(
        mutant_pdb,
        leg_dir / "mutant_geometry_qc.json",
        reference_path=prepared_input if prepared_input.is_file() else None,
    )
    remaining = list(qc_payload.get("inter_residue_heavy_atom_clashes", []))
    summary["post_regeneration_clash_count"] = len(remaining)
    summary["succeeded"] = not remaining
    return remaining, summary


def _stage_build_legs(ctx: StageContext) -> StageStatus:
    started = utc_now()
    lambda_values = _lambda_values(ctx.protocol.lambda_windows)
    env = _stage_env(ctx)
    artifacts: list[str] = []
    gmxlib = env.get("GMXLIB") if env else None
    if gmxlib:
        artifacts.append(gmxlib)
    command_lines: list[str] = []
    gmxlib_dir = Path(gmxlib) if gmxlib else None

    for leg in ("complex", "apo"):
        pmx_dir = ctx.job_dir / "legs" / leg / "pmx"
        processed_gro = pmx_dir / "processed.gro"
        pmx_top = pmx_dir / "pmxtop.top"
        for repeat_index in range(1, ctx.protocol.repeats + 1):
            repeat_dir = ctx.job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            setup_dir = repeat_dir / "setup"
            equil_dir = repeat_dir / "equilibration"
            setup_dir.mkdir(parents=True, exist_ok=True)
            equil_dir.mkdir(parents=True, exist_ok=True)
            mdp_artifacts = _write_equilibration_mdps(
                ctx,
                leg=leg,
                repeat_index=repeat_index,
                repeat_dir=repeat_dir,
            )
            mdp_dir = repeat_dir / "mdp"

            manifest = {
                "leg": leg,
                "repeat": repeat_index,
                "lambda_values": [round(value, 5) for value in lambda_values],
                "processed_gro": str(processed_gro),
                "topology": str(pmx_top),
                "repeat_topology": str(repeat_dir / "system.top"),
                "equilibration_restraint_schedule": _equilibration_restraint_schedule(ctx),
                "equilibration_release_npt_ps": _equilibration_release_npt_ps(ctx),
                "gmxlib": str(gmxlib_dir) if gmxlib_dir else None,
                "water_coordinates": str(water_coordinate_path(gmxlib_dir, ctx.protocol.water_model)) if gmxlib_dir else None,
            }
            write_json(repeat_dir / "lambda_plan.json", manifest)
            artifacts.extend(
                [
                    str(repeat_dir / "lambda_plan.json"),
                    *mdp_artifacts,
                ]
            )

            for window_index, _lambda_value in enumerate(lambda_values):
                window_dir = repeat_dir / f"lambda_{window_index:03d}"
                window_dir.mkdir(parents=True, exist_ok=True)
                pre_relax_mdp = window_dir / "pre_relax.mdp"
                pre_md_mdp = window_dir / "pre_md.mdp"
                mdp_path = window_dir / "production.mdp"
                pre_relax_mdp.write_text(
                    _render_window_relax_em_mdp(ctx, lambda_values, window_index),
                    encoding="utf-8",
                )
                pre_md_mdp.write_text(
                    _render_window_relax_md_mdp(
                        ctx,
                        lambda_values,
                        window_index,
                        _job_seed(ctx, leg, repeat_index, window_index, "pre-md"),
                    ),
                    encoding="utf-8",
                )
                mdp_path.write_text(
                    _render_lambda_mdp(
                        ctx,
                        lambda_values,
                        window_index,
                        _job_seed(ctx, leg, repeat_index, window_index, "prod"),
                    ),
                    encoding="utf-8",
                )
                artifacts.extend([str(pre_relax_mdp), str(pre_md_mdp), str(mdp_path)])
                command_lines.append(
                    f"# {leg} rep{repeat_index:02d} lambda_{window_index:03d}: window relax will start from {equil_dir / 'npt.gro'} before production"
                )

    script_file = ctx.job_dir / "artifacts" / "commands" / "build_legs.sh"
    ctx.runner.write_script(script_file, command_lines, ctx.job_dir, env=env)
    status = StageStatus(
        stage="build_legs",
        state="completed",
        message="Per-leg directories, GMXLIB overlay, and MDP files generated.",
        artifacts=artifacts + [str(script_file)],
        commands=command_lines,
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_equilibrate(ctx: StageContext) -> StageStatus:
    started = utc_now()
    env = _stage_env(ctx)
    if ctx.runner.execute and env is None:
        status = StageStatus(
            stage="equilibrate",
            state="blocked_external",
            message="Could not determine a usable GMXLIB overlay for solvation and equilibration.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    gmxlib = env.get("GMXLIB") if env else None
    gmxlib_dir = Path(gmxlib) if gmxlib else None
    if ctx.runner.execute and gmxlib_dir is not None:
        solvent_coordinates = water_coordinate_path(gmxlib_dir, ctx.protocol.water_model)
    else:
        solvent_coordinates = (gmxlib_dir / "spc216.gro") if gmxlib_dir is not None else Path("spc216.gro")

    gmx_command = _resolve_gmx_command(ctx)
    mdrun_suffix = _mdrun_suffix(ctx)
    commands: list[str] = []
    artifacts: list[str] = []
    if gmxlib:
        artifacts.append(gmxlib)

    for leg in ("complex", "apo"):
        pmx_dir = ctx.job_dir / "legs" / leg / "pmx"
        processed_gro = pmx_dir / "processed.gro"
        pmx_top = pmx_dir / "pmxtop.top"
        for repeat_index in range(1, ctx.protocol.repeats + 1):
            repeat_dir = ctx.job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            setup_dir = repeat_dir / "setup"
            equil_dir = repeat_dir / "equilibration"
            mdp_dir = repeat_dir / "mdp"
            setup_dir.mkdir(parents=True, exist_ok=True)
            equil_dir.mkdir(parents=True, exist_ok=True)
            mdp_artifacts = _write_equilibration_mdps(
                ctx,
                leg=leg,
                repeat_index=repeat_index,
                repeat_dir=repeat_dir,
            )
            repeat_top = repeat_dir / "system.top"
            support_itps = sorted(pmx_dir.glob("*.itp"))
            if ctx.runner.execute and not processed_gro.is_file():
                status = StageStatus(
                    stage="equilibrate",
                    state="blocked_input",
                    message=f"Mutated coordinate file is missing: {processed_gro}. Run mutate with --execute first.",
                    started_at=started,
                    completed_at=utc_now(),
                )
                return _write_stage_status(ctx.job_dir, status)
            if ctx.runner.execute and processed_gro.is_file():
                processed_gro_qc = inspect_gro_file(processed_gro)
                if not processed_gro_qc.get("valid"):
                    reason = str(processed_gro_qc.get("reason") or "invalid_processed_gro")
                    line_number = processed_gro_qc.get("line_number")
                    line_suffix = f" at line {line_number}" if line_number is not None else ""
                    status = StageStatus(
                        stage="equilibrate",
                        state="blocked_input",
                        message=(
                            f"Mutated coordinate file is invalid: {processed_gro} ({reason}{line_suffix}). "
                            "Rerun mutate with --execute first."
                        ),
                        started_at=started,
                        completed_at=utc_now(),
                    )
                    return _write_stage_status(ctx.job_dir, status)
            if ctx.runner.execute and not pmx_top.is_file():
                status = StageStatus(
                    stage="equilibrate",
                    state="blocked_input",
                    message=f"Hybrid topology is missing: {pmx_top}. Run mutate with --execute first.",
                    started_at=started,
                    completed_at=utc_now(),
                )
                return _write_stage_status(ctx.job_dir, status)

            boxed_gro = setup_dir / "boxed.gro"
            solvated_gro = setup_dir / "solvated.gro"
            genion_tpr = setup_dir / "genion.tpr"
            ions_gro = setup_dir / "ions.gro"
            em_tpr = equil_dir / "em.tpr"
            nvt_tpr = equil_dir / "nvt.tpr"
            npt_tpr = equil_dir / "npt.tpr"
            npt_stage1_tpr = equil_dir / "npt_stage1.tpr"
            npt_release_tpr = equil_dir / "npt_release.tpr"
            em_deffnm = equil_dir / "em"
            nvt_deffnm = equil_dir / "nvt"
            npt_deffnm = equil_dir / "npt"
            npt_stage1_deffnm = equil_dir / "npt_stage1"
            npt_release_deffnm = equil_dir / "npt_release"
            em_runtime_log = equil_dir / "em.runtime.log"
            restraint_summary_path = repeat_dir / "equilibration_restraints.json"
            release_artifacts = []
            if _has_equilibration_release_stage(ctx):
                release_artifacts.extend(
                    [
                        str(npt_stage1_tpr),
                        str(npt_release_tpr),
                        str(npt_stage1_deffnm.with_suffix(".gro")),
                        str(npt_release_deffnm.with_suffix(".gro")),
                    ]
                )
            restraint_artifacts = [str(restraint_summary_path)] if _uses_staged_equilibration_restraints(ctx) else []

            artifacts.extend(
                [
                    *mdp_artifacts,
                    str(repeat_top),
                    str(boxed_gro),
                    str(solvated_gro),
                    str(genion_tpr),
                    str(ions_gro),
                    str(em_tpr),
                    str(nvt_tpr),
                    str(npt_tpr),
                    *release_artifacts,
                    str(em_runtime_log),
                    str(npt_deffnm.with_suffix(".gro")),
                    *restraint_artifacts,
                ]
            )
            artifacts.extend(str(repeat_dir / support_itp.name) for support_itp in support_itps)
            if ctx.runner.execute and _should_seed_equilibrated_repeat_from_source(ctx, leg):
                artifacts.extend(
                    _seed_equilibrated_repeat_from_sources(
                        ctx,
                        leg=leg,
                        repeat_dir=repeat_dir,
                        repeat_top=repeat_top,
                        support_itps=support_itps,
                    )
                )
            commands.append(
                _equilibrate_repeat_snippet(
                    ctx=ctx,
                    gmx_command=gmx_command,
                    mdrun_suffix=mdrun_suffix,
                    processed_gro=processed_gro,
                    pmx_top=pmx_top,
                    support_itps=support_itps,
                    repeat_dir=repeat_dir,
                    setup_dir=setup_dir,
                    equil_dir=equil_dir,
                    mdp_dir=mdp_dir,
                    solvent_coordinates=solvent_coordinates,
                )
            )

    script_file = ctx.job_dir / "artifacts" / "commands" / "equilibrate.sh"
    _write_external_stage_running(
        ctx.job_dir,
        stage="equilibrate",
        script_file=script_file,
        commands=commands,
        artifacts=artifacts,
        started_at=started,
    )
    outcome = ctx.runner.run_script(script_file, commands, ctx.job_dir, env=env)
    status = StageStatus(
        stage="equilibrate",
        state=outcome.state,
        message=outcome.message,
        commands=commands,
        artifacts=_external_stage_artifacts(script_file, artifacts),
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_sample(ctx: StageContext) -> StageStatus:
    started = utc_now()
    env = _stage_env(ctx)
    if ctx.runner.execute and env is None:
        status = StageStatus(
            stage="sample",
            state="blocked_external",
            message="Could not determine a usable GMXLIB overlay for lambda sampling.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    gmx_command = _resolve_gmx_command(ctx)
    mdrun_suffix = _mdrun_suffix(ctx)
    commands: list[str] = []
    artifacts: list[str] = []
    gmxlib = env.get("GMXLIB") if env else None
    if gmxlib:
        artifacts.append(gmxlib)

    for leg in ("complex", "apo"):
        for repeat_index in range(1, ctx.protocol.repeats + 1):
            repeat_dir = ctx.job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            repeat_top = repeat_dir / "system.top"
            start_gro = repeat_dir / "equilibration" / "npt.gro"
            if ctx.runner.execute and not start_gro.is_file():
                status = StageStatus(
                    stage="sample",
                    state="blocked_input",
                    message=f"Equilibrated starting structure is missing: {start_gro}. Run equilibrate with --execute first.",
                    started_at=started,
                    completed_at=utc_now(),
                )
                return _write_stage_status(ctx.job_dir, status)
            if ctx.runner.execute and not repeat_top.is_file():
                status = StageStatus(
                    stage="sample",
                    state="blocked_input",
                    message=f"Repeat-specific topology is missing: {repeat_top}. Run equilibrate with --execute first.",
                    started_at=started,
                    completed_at=utc_now(),
                )
                return _write_stage_status(ctx.job_dir, status)
            if ctx.runner.execute:
                artifacts.extend(_backfill_repeat_support_itps_from_seed_source(repeat_dir, repeat_top))

            for window_index, _lambda_value in enumerate(_lambda_values(ctx.protocol.lambda_windows)):
                window_dir = repeat_dir / f"lambda_{window_index:03d}"
                if ctx.runner.execute and _should_inherit_leg_from_source(ctx, leg):
                    artifacts.extend(
                        _seed_sample_window_from_source(
                            ctx,
                            leg=leg,
                            repeat_dir=repeat_dir,
                            window_dir=window_dir,
                        )
                    )
                artifacts.extend(
                    [
                        str(window_dir / "pre_relax.tpr"),
                        str((window_dir / "pre_relax").with_suffix(".gro")),
                        str(window_dir / "pre_md.tpr"),
                        str((window_dir / "pre_md").with_suffix(".gro")),
                        str((window_dir / "pre_md").with_suffix(".cpt")),
                        str(window_dir / "topol.tpr"),
                        str(window_dir / "dhdl.xvg"),
                    ]
                )
                commands.append(
                    _sample_window_snippet(
                        gmx_command=gmx_command,
                        mdrun_suffix=mdrun_suffix,
                        repeat_top=repeat_top,
                        start_gro=start_gro,
                        window_dir=window_dir,
                        grompp_maxwarn_sampling=ctx.protocol.grompp_maxwarn_sampling,
                    )
                )

    script_file = ctx.job_dir / "artifacts" / "commands" / "sample.sh"
    _write_external_stage_running(
        ctx.job_dir,
        stage="sample",
        script_file=script_file,
        commands=commands,
        artifacts=artifacts,
        started_at=started,
    )
    outcome = ctx.runner.run_script(script_file, commands, ctx.job_dir, env=env)
    status = StageStatus(
        stage="sample",
        state=outcome.state,
        message=outcome.message,
        commands=commands,
        artifacts=_external_stage_artifacts(script_file, artifacts),
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_bar(ctx: StageContext) -> StageStatus:
    started = utc_now()
    env = _stage_env(ctx)
    if ctx.runner.execute and env is None:
        status = StageStatus(
            stage="bar",
            state="blocked_external",
            message="Could not determine a usable GMXLIB overlay for BAR analysis.",
            started_at=started,
            completed_at=utc_now(),
        )
        return _write_stage_status(ctx.job_dir, status)

    gmx_command = _resolve_gmx_command(ctx)
    commands: list[str] = []
    artifacts: list[str] = []
    gmxlib = env.get("GMXLIB") if env else None
    if gmxlib:
        artifacts.append(gmxlib)

    for leg in ("complex", "apo"):
        for repeat_index in range(1, ctx.protocol.repeats + 1):
            repeat_dir = ctx.job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            dhdl_files = sorted(str(path) for path in repeat_dir.glob("lambda_*/dhdl.xvg"))
            output_dir = repeat_dir / "bar"
            output_dir.mkdir(parents=True, exist_ok=True)
            if dhdl_files:
                commands.append(
                    _bar_repeat_snippet(
                        gmx_command=gmx_command,
                        output_dir=output_dir,
                        dhdl_files=[Path(path) for path in dhdl_files],
                    )
                )
            artifacts.append(str(output_dir))
    script_file = ctx.job_dir / "artifacts" / "commands" / "bar.sh"
    if commands:
        _write_external_stage_running(
            ctx.job_dir,
            stage="bar",
            script_file=script_file,
            commands=commands,
            artifacts=artifacts,
            started_at=started,
        )
    outcome = ctx.runner.run_script(script_file, commands, ctx.job_dir, env=env)
    status = StageStatus(
        stage="bar",
        state=outcome.state if commands else "planned",
        message=outcome.message if commands else "No dhdl.xvg files found yet; BAR stage remains planned.",
        commands=commands,
        artifacts=_external_stage_artifacts(script_file, artifacts),
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_qc(ctx: StageContext) -> StageStatus:
    started = utc_now()
    result_payload = write_job_results(ctx.job_dir)
    qc_path = Path(result_payload["paths"]["qc_report"])
    status = StageStatus(
        stage="qc",
        state="completed",
        message=f"QC report written with status {result_payload['qc_report']['status']}.",
        artifacts=[
            str(qc_path),
            result_payload["paths"]["bar_summary"],
            result_payload["paths"]["ddg_summary"],
            result_payload["paths"]["ddg_summary_tsv"],
        ],
        started_at=started,
        completed_at=utc_now(),
    )
    return _write_stage_status(ctx.job_dir, status)


def _stage_report(ctx: StageContext) -> StageStatus:
    started = utc_now()
    status = StageStatus(
        stage="report",
        state="completed",
        message="Job summary written.",
        artifacts=[str(ctx.job_dir / "report" / "summary.json"), str(ctx.job_dir / "report" / "summary.yml")],
        started_at=started,
        completed_at=utc_now(),
    )
    _write_stage_status(ctx.job_dir, status)
    write_job_summary(ctx.job_dir)
    return status


STAGE_DISPATCH = {
    "ingest": _stage_ingest,
    "prepare": _stage_prepare,
    "mutate": _stage_mutate,
    "build_legs": _stage_build_legs,
    "equilibrate": _stage_equilibrate,
    "sample": _stage_sample,
    "bar": _stage_bar,
    "qc": _stage_qc,
    "report": _stage_report,
}


def sh_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _invalidate_job_outputs(job_dir: Path, start_index: int) -> None:
    stages_dir = job_dir / "stages"
    for stage in STAGES[start_index:]:
        stage_path = stages_dir / f"{stage}.json"
        if stage_path.exists():
            stage_path.unlink()

    results_dir = job_dir / "results"
    if results_dir.exists():
        for result_name in ("bar_summary.json", "ddg_summary.json", "ddg_summary.tsv", "qc_report.json"):
            result_path = results_dir / result_name
            if result_path.exists():
                result_path.unlink()

    report_dir = job_dir / "report"
    if report_dir.exists():
        for summary_name in ("summary.json", "summary.yml"):
            summary_path = report_dir / summary_name
            if summary_path.exists():
                summary_path.unlink()


def _write_stage_running(job_dir: Path, stage: str) -> None:
    _write_stage_status(
        job_dir,
        StageStatus(
            stage=stage,
            state="running",
            message="Stage execution started.",
            started_at=utc_now(),
            completed_at=None,
        ),
    )


def _interrupted_stage_status(job_dir: Path, stage: str, message: str) -> StageStatus:
    stage_path = job_dir / "stages" / f"{stage}.json"
    started_at = None
    if stage_path.exists():
        started_at = read_json(stage_path).get("started_at")
    return StageStatus(
        stage=stage,
        state="failed",
        message=message,
        started_at=started_at or utc_now(),
        completed_at=utc_now(),
    )


def _job_execution_lock_path(job_dir: Path) -> Path:
    return job_dir / JOB_EXECUTION_LOCK


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _iter_process_commands() -> list[tuple[int, str]]:
    proc_root = Path("/proc")
    if proc_root.is_dir():
        commands: list[tuple[int, str]] = []
        for proc_dir in proc_root.iterdir():
            if not proc_dir.name.isdigit():
                continue
            pid = int(proc_dir.name)
            try:
                raw_cmdline = (proc_dir / "cmdline").read_bytes()
            except OSError:
                continue
            command = raw_cmdline.replace(b"\x00", b" ").decode("utf-8", errors="ignore").strip()
            if command:
                commands.append((pid, command))
        if commands:
            return commands

    try:
        result = subprocess.run(
            ["ps", "-eo", "pid,args"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError:
        return []
    if result.returncode != 0:
        return []

    commands = []
    for raw_line in result.stdout.splitlines()[1:]:
        line = raw_line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2 or not parts[0].isdigit():
            continue
        commands.append((int(parts[0]), parts[1]))
    return commands


def _job_active_process(job_dir: Path) -> tuple[int, str] | None:
    job_dir_text = str(job_dir)
    command_markers = (
        str(job_dir / "artifacts" / "commands"),
        str(job_dir / "legs"),
    )

    current_pid = os.getpid()
    for pid, command in _iter_process_commands():
        if pid == current_pid:
            continue
        if job_dir_text not in command:
            continue
        if not any(marker in command for marker in command_markers):
            continue
        return pid, command
    return None


def _job_execution_conflict(job_dir: Path) -> str | None:
    lock_path = _job_execution_lock_path(job_dir)
    if lock_path.exists():
        try:
            payload = json.loads(lock_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            payload = {}
        pid = int(payload.get("pid") or 0)
        started_at = str(payload.get("started_at") or "").strip()
        if pid and _pid_is_alive(pid):
            detail = f" by PID {pid}"
            if started_at:
                detail += f" since {started_at}"
            return "Job execution already in progress" + detail + "."
        try:
            lock_path.unlink()
        except OSError:
            pass

    active_process = _job_active_process(job_dir)
    if active_process is not None:
        pid, command = active_process
        return f"Job execution already in progress by PID {pid}: {command}"
    return None


def _acquire_job_execution_lock(
    job_dir: Path,
    *,
    from_stage: str | None,
    to_stage: str | None,
    execute: bool,
) -> tuple[Path | None, str | None]:
    conflict = _job_execution_conflict(job_dir)
    if conflict is not None:
        return None, conflict

    lock_path = _job_execution_lock_path(job_dir)
    payload = {
        "pid": os.getpid(),
        "started_at": utc_now(),
        "from_stage": from_stage or "",
        "to_stage": to_stage or "",
        "execute": execute,
    }
    while True:
        try:
            fd = os.open(lock_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        except FileExistsError:
            conflict = _job_execution_conflict(job_dir)
            if conflict is not None:
                return None, conflict
            continue
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        return lock_path, None


def _release_job_execution_lock(lock_path: Path | None) -> None:
    if lock_path is None:
        return
    try:
        if lock_path.exists():
            lock_path.unlink()
    except OSError:
        pass


def _duplicate_execution_status(job_dir: Path, stage: str, message: str) -> StageStatus:
    return StageStatus(
        stage=stage,
        state="running",
        message=message,
        started_at=utc_now(),
        completed_at=utc_now(),
        artifacts=[str(_job_execution_lock_path(job_dir))] if _job_execution_lock_path(job_dir).exists() else [],
    )


def _job_has_any(job_dir: Path, *patterns: str) -> bool:
    for pattern in patterns:
        if any(job_dir.glob(pattern)):
            return True
    return False


def _job_protocol_counts(job_dir: Path) -> tuple[int, int]:
    spec_path = job_dir / "job_spec.json"
    if not spec_path.exists():
        return 0, 0
    try:
        payload = read_json(spec_path)
    except OSError:
        return 0, 0
    protocol = payload.get("protocol", {})
    try:
        repeats = max(int(protocol.get("repeats", 0)), 0)
    except (TypeError, ValueError):
        repeats = 0
    try:
        lambda_windows = max(int(protocol.get("lambda_windows", 0)), 0)
    except (TypeError, ValueError):
        lambda_windows = 0
    return repeats, lambda_windows


def _existing_stage_payload(job_dir: Path, stage: str) -> dict[str, object]:
    stage_path = job_dir / "stages" / f"{stage}.json"
    if not stage_path.exists():
        return {}
    try:
        payload = read_json(stage_path)
    except OSError:
        return {}
    if isinstance(payload, dict):
        return payload
    return {}


def _stage_artifact_exists(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _mutate_outputs_complete(job_dir: Path) -> bool:
    pmx_dirs = [job_dir / "legs" / "complex" / "pmx", job_dir / "legs" / "apo" / "pmx"]
    return all(
        _stage_artifact_exists(pmx_dir / "processed.gro")
        and gro_file_is_valid(pmx_dir / "processed.gro")
        and _stage_artifact_exists(pmx_dir / "pmxtop.top")
        for pmx_dir in pmx_dirs
    )


def _build_legs_outputs_complete(job_dir: Path) -> bool:
    repeats, lambda_windows = _job_protocol_counts(job_dir)
    if repeats <= 0 or lambda_windows <= 0:
        return False
    for leg in ("complex", "apo"):
        for repeat_index in range(1, repeats + 1):
            repeat_dir = job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            required_files = [
                repeat_dir / "lambda_plan.json",
                repeat_dir / "mdp" / "genion.mdp",
                repeat_dir / "mdp" / "em.mdp",
                repeat_dir / "mdp" / "nvt.mdp",
                repeat_dir / "mdp" / "npt.mdp",
            ]
            if not all(_stage_artifact_exists(path) for path in required_files):
                return False
            for window_index in range(lambda_windows):
                lambda_dir = repeat_dir / f"lambda_{window_index:03d}"
                lambda_files = [
                    lambda_dir / "pre_relax.mdp",
                    lambda_dir / "pre_md.mdp",
                    lambda_dir / "production.mdp",
                ]
                if not all(_stage_artifact_exists(path) for path in lambda_files):
                    return False
    return True


def _equilibrate_outputs_complete(job_dir: Path) -> bool:
    repeats, _lambda_windows = _job_protocol_counts(job_dir)
    if repeats <= 0:
        return False
    for leg in ("complex", "apo"):
        for repeat_index in range(1, repeats + 1):
            repeat_dir = job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            if not _stage_artifact_exists(repeat_dir / "system.top"):
                return False
            if not _stage_artifact_exists(repeat_dir / "equilibration" / "npt.gro"):
                return False
    return True


def _sample_outputs_complete(job_dir: Path) -> bool:
    repeats, lambda_windows = _job_protocol_counts(job_dir)
    if repeats <= 0 or lambda_windows <= 0:
        return False
    for leg in ("complex", "apo"):
        for repeat_index in range(1, repeats + 1):
            repeat_dir = job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            for window_index in range(lambda_windows):
                lambda_dir = repeat_dir / f"lambda_{window_index:03d}"
                if not _stage_artifact_exists(lambda_dir / "dhdl.xvg"):
                    return False
                if not _stage_artifact_exists(lambda_dir / "md.gro"):
                    return False
    return True


def _recovered_stage_status(
    job_dir: Path,
    stage: str,
    *,
    message: str,
) -> StageStatus:
    payload = _existing_stage_payload(job_dir, stage)
    commands = payload.get("commands")
    artifacts = payload.get("artifacts")
    started_at = str(payload.get("started_at") or "").strip() or utc_now()
    completed_at = str(payload.get("completed_at") or "").strip() or utc_now()
    return StageStatus(
        stage=stage,
        state="completed",
        message=message,
        commands=list(commands) if isinstance(commands, list) else [],
        artifacts=list(artifacts) if isinstance(artifacts, list) else [],
        started_at=started_at,
        completed_at=completed_at,
    )


def _job_id_for_recovery(job_dir: Path) -> str:
    spec_path = job_dir / "job_spec.json"
    if spec_path.exists():
        try:
            payload = read_json(spec_path)
        except OSError:
            payload = {}
        job_id = str(payload.get("job_id") or "").strip()
        if job_id:
            return job_id
    return job_dir.name


def _ensure_mutate_qc_payload(job_dir: Path) -> dict[str, object]:
    qc_path = job_dir / "artifacts" / "mutate_qc.json"
    existing_payload = _read_json_dict_if_exists(qc_path)
    payload: dict[str, object] = {
        "job_id": _job_id_for_recovery(job_dir),
        "legs": {},
    }
    found_any_mutant = False
    for leg in ("complex", "apo"):
        leg_payload, _artifacts, missing_mutant = _collect_mutate_leg_qc(
            job_dir,
            leg,
            geometry_qc_path=job_dir / "legs" / leg / "pmx" / "mutant_geometry_qc.json",
            restore_summary_path=job_dir / "legs" / leg / "pmx" / "mutant_standard_residue_repair.json",
            pdbfixer_summary_path=job_dir / "legs" / leg / "pmx" / "mutant_pdbfixer_repair.json",
            processed_gro_qc_path=job_dir / "legs" / leg / "pmx" / "processed_gro_qc.json",
        )
        payload["legs"][leg] = leg_payload
        if not missing_mutant:
            found_any_mutant = True
    if found_any_mutant:
        write_json(qc_path, payload)
        return payload
    if existing_payload:
        return existing_payload
    if qc_path.exists():
        try:
            qc_path.unlink()
        except OSError:
            pass
    return {}


def _blocking_mutate_qc_legs_from_payload(payload: dict[str, object]) -> dict[str, list[dict[str, object]]]:
    legs = payload.get("legs", {})
    if not isinstance(legs, dict):
        return {}
    blocking: dict[str, list[dict[str, object]]] = {}
    for leg, leg_payload in legs.items():
        if not isinstance(leg_payload, dict):
            continue
        issues = leg_payload.get("inter_residue_heavy_atom_clashes", [])
        if not isinstance(issues, list):
            continue
        filtered = [
            issue
            for issue in issues
            if isinstance(issue, dict) and bool(issue.get("blocking_prepare", True))
        ]
        if filtered:
            blocking[str(leg)] = filtered
    return blocking


def _blocking_mutate_qc_legs(job_dir: Path) -> dict[str, list[dict[str, object]]]:
    return _blocking_mutate_qc_legs_from_payload(_ensure_mutate_qc_payload(job_dir))


def _repair_blocking_mutate_qc(
    job_dir: Path,
    *,
    execute: bool,
    environment: dict[str, str] | None = None,
) -> dict[str, list[dict[str, object]]]:
    payload = _ensure_mutate_qc_payload(job_dir)
    if not payload:
        return {}
    blocking = _blocking_mutate_qc_legs_from_payload(payload)
    if not blocking or not execute:
        return blocking

    try:
        ctx = _load_context(job_dir, execute=True, environment=environment)
    except (FileNotFoundError, KeyError, TypeError, ValueError):
        return blocking
    legs = payload.get("legs", {})
    if not isinstance(legs, dict):
        return blocking

    updated = False
    for leg, issues in list(blocking.items()):
        leg_payload = legs.get(leg)
        if not isinstance(leg_payload, dict):
            continue
        repaired_issues, repair_summary = _attempt_repair_mutant_sidechain_clashes(ctx, leg=leg, issues=issues)
        leg_payload["inter_residue_heavy_atom_clashes"] = repaired_issues
        leg_payload["auto_repair_summary"] = repair_summary
        if repaired_issues != issues or repair_summary.get("attempted"):
            updated = True
    if updated:
        write_json(job_dir / "artifacts" / "mutate_qc.json", payload)
    return _blocking_mutate_qc_legs_from_payload(payload)


def _blocking_mutate_qc_message(blocking_legs: dict[str, list[dict[str, object]]]) -> str:
    samples: list[str] = []
    for leg, issues in blocking_legs.items():
        for issue in issues[:2]:
            insertion = issue.get("icode", "") or ""
            partner_insertion = issue.get("partner_icode", "") or ""
            clashes = issue.get("clashes", [])
            if not isinstance(clashes, list) or not clashes:
                continue
            first_clash = clashes[0]
            if not isinstance(first_clash, dict):
                continue
            samples.append(
                f"{leg}:{issue.get('chain_id', '?')}{issue.get('resseq', '?')}{insertion} {issue.get('resname', '?')} vs "
                f"{issue.get('partner_chain_id', '?')}{issue.get('partner_resseq', '?')}{partner_insertion} {issue.get('partner_resname', '?')} "
                f"{first_clash.get('atom_a', '?')}-{first_clash.get('atom_b', '?')} {float(first_clash.get('distance_angstrom', 0.0)):.3f} A"
            )
    extra_count = sum(len(issues) for issues in blocking_legs.values()) - len(samples)
    suffix = f" (+{extra_count} more)" if extra_count > 0 else ""
    return (
        "Mutated structure contains impossible inter-residue heavy-atom clashes: "
        + "; ".join(samples)
        + suffix
    )


def _recovered_blocked_mutate_status(
    job_dir: Path,
    blocking_legs: dict[str, list[dict[str, object]]],
) -> StageStatus:
    payload = _existing_stage_payload(job_dir, "mutate")
    commands = payload.get("commands")
    artifacts = payload.get("artifacts")
    started_at = str(payload.get("started_at") or "").strip() or utc_now()
    completed_at = utc_now()
    return StageStatus(
        stage="mutate",
        state="blocked_input",
        message=_blocking_mutate_qc_message(blocking_legs),
        commands=list(commands) if isinstance(commands, list) else [],
        artifacts=list(artifacts) if isinstance(artifacts, list) else [],
        started_at=started_at,
        completed_at=completed_at,
    )


def _recover_completed_stage_files(
    job_dir: Path,
    stage_states: dict[str, str],
    *,
    execute: bool,
    environment: dict[str, str] | None = None,
) -> dict[str, str]:
    mutate_state = stage_states.get("mutate")
    has_downstream_stage_state = any(stage in stage_states for stage in STAGES[STAGES.index("build_legs") :])
    if not _mutate_outputs_complete(job_dir) and (mutate_state != "blocked_input" or has_downstream_stage_state):
        _invalidate_job_outputs(job_dir, STAGES.index("mutate"))
        for stage in STAGES[STAGES.index("mutate") :]:
            stage_states.pop(stage, None)
        return stage_states

    blocking_mutate_legs = _repair_blocking_mutate_qc(job_dir, execute=execute, environment=environment)
    if blocking_mutate_legs:
        _invalidate_job_outputs(job_dir, STAGES.index("build_legs"))
        _write_stage_status(
            job_dir,
            _recovered_blocked_mutate_status(job_dir, blocking_mutate_legs),
        )
        stage_states["mutate"] = "blocked_input"
        for stage in STAGES[STAGES.index("build_legs") :]:
            stage_states.pop(stage, None)
        return stage_states

    recovered_checks = [
        (
            "mutate",
            _mutate_outputs_complete(job_dir),
            "Recovered completed stage from existing pmx outputs.",
        ),
        (
            "build_legs",
            _build_legs_outputs_complete(job_dir),
            "Recovered completed stage from existing lambda-plan and MDP outputs.",
        ),
        (
            "equilibrate",
            _equilibrate_outputs_complete(job_dir),
            "Recovered completed stage from existing equilibrated repeat outputs.",
        ),
        (
            "sample",
            _sample_outputs_complete(job_dir),
            "Recovered completed stage from existing lambda sampling outputs.",
        ),
    ]
    for stage, complete, message in recovered_checks:
        if not complete or stage_states.get(stage) == "completed":
            continue
        _write_stage_status(
            job_dir,
            _recovered_stage_status(
                job_dir,
                stage,
                message=message,
            ),
        )
        stage_states[stage] = "completed"
    return stage_states


def _job_repeat_dirs(job_dir: Path) -> list[Path]:
    return sorted(path for path in job_dir.glob("legs/*/rep*") if path.is_dir())


def _job_pmx_dirs(job_dir: Path) -> list[Path]:
    return sorted(path for path in job_dir.glob("legs/*/pmx") if path.is_dir())


def _recoverable_resume_stage(job_dir: Path, stage_states: dict[str, str]) -> str | None:
    sample_outputs_exist = _job_has_any(
        job_dir,
        "legs/*/rep*/lambda_*/dhdl.xvg",
        "legs/*/rep*/lambda_*/topol.tpr",
        "legs/*/rep*/lambda_*/pre_md.tpr",
        "legs/*/rep*/lambda_*/pre_relax.tpr",
    )
    equilibrate_outputs_exist = _job_has_any(
        job_dir,
        "legs/*/rep*/equilibration/npt.gro",
        "legs/*/rep*/equilibration/nvt.gro",
        "legs/*/rep*/equilibration/em.gro",
        "legs/*/rep*/equilibration/*.tpr",
        "legs/*/rep*/system.top",
        "legs/*/rep*/setup/*.gro",
        "legs/*/rep*/setup/*.tpr",
    )
    build_outputs_exist = _job_has_any(
        job_dir,
        "legs/*/rep*/lambda_plan.json",
        "legs/*/rep*/mdp/*.mdp",
    )
    mutate_outputs_exist = _job_has_any(
        job_dir,
        "legs/*/pmx/processed.gro",
        "legs/*/pmx/pmxtop.top",
    )
    mutate_outputs_complete = _mutate_outputs_complete(job_dir)
    repeat_dirs = _job_repeat_dirs(job_dir)
    pmx_dirs = _job_pmx_dirs(job_dir)
    sample_prereqs_complete = bool(repeat_dirs) and all((repeat_dir / "equilibration" / "npt.gro").is_file() for repeat_dir in repeat_dirs) and all(
        (repeat_dir / "system.top").is_file() for repeat_dir in repeat_dirs
    )
    equilibrate_prereqs_complete = bool(pmx_dirs) and mutate_outputs_complete

    if stage_states.get("sample") != "completed" and sample_outputs_exist:
        if sample_prereqs_complete:
            return "sample"

    if stage_states.get("equilibrate") != "completed" and equilibrate_outputs_exist:
        if equilibrate_prereqs_complete:
            return "equilibrate"

    if stage_states.get("build_legs") != "completed" and build_outputs_exist:
        if mutate_outputs_complete:
            return "build_legs"

    if stage_states.get("mutate") != "completed" and mutate_outputs_complete:
        return "mutate"

    return None


def run_job(
    job_dir: Path,
    execute: bool,
    from_stage: str | None = None,
    to_stage: str | None = None,
    environment: dict[str, str] | None = None,
) -> list[StageStatus]:
    ctx = _load_context(job_dir, execute=execute, environment=environment)
    start_index = STAGES.index(from_stage) if from_stage else 0
    stop_index = STAGES.index(to_stage) if to_stage else len(STAGES) - 1
    if stop_index < start_index:
        raise ValueError("to-stage must be after from-stage")
    lock_path, conflict = _acquire_job_execution_lock(
        job_dir,
        from_stage=from_stage,
        to_stage=to_stage,
        execute=execute,
    )
    if conflict is not None:
        return [_duplicate_execution_status(job_dir, STAGES[start_index], conflict)]
    try:
        _invalidate_job_outputs(job_dir, start_index)
        statuses = []
        for stage in STAGES[start_index : stop_index + 1]:
            _write_stage_running(ctx.job_dir, stage)
            try:
                status = STAGE_DISPATCH[stage](ctx)
            except KeyboardInterrupt:
                status = _write_stage_status(
                    ctx.job_dir,
                    _interrupted_stage_status(ctx.job_dir, stage, "Stage execution interrupted."),
                )
            statuses.append(status)
            if status.state in {"failed", "blocked_external", "blocked_input"}:
                break
        return statuses
    finally:
        _release_job_execution_lock(lock_path)


def resume_job(job_dir: Path, execute: bool, environment: dict[str, str] | None = None) -> list[StageStatus]:
    stage_states: dict[str, str] = {}
    for stage_file in (job_dir / "stages").glob("*.json"):
        payload = read_json(stage_file)
        state = str(payload.get("state", ""))
        stage_states[stage_file.stem] = state
    if _job_execution_conflict(job_dir) is None:
        stage_states = _recover_completed_stage_files(
            job_dir,
            stage_states,
            execute=execute,
            environment=environment,
        )

    blocked_stage = next((stage for stage in STAGES if stage_states.get(stage) == "blocked_input"), None)
    completed = {stage for stage, state in stage_states.items() if state == "completed"}
    next_stage = None
    for stage in STAGES:
        if stage not in completed:
            next_stage = stage
            break
    if blocked_stage is not None:
        next_stage = blocked_stage
    recoverable_stage = _recoverable_resume_stage(job_dir, stage_states)
    if blocked_stage is None and recoverable_stage is not None and next_stage is not None:
        if STAGES.index(recoverable_stage) > STAGES.index(next_stage):
            next_stage = recoverable_stage
    if next_stage is None:
        return []
    return run_job(job_dir, execute=execute, from_stage=next_stage, environment=environment)
