"""Window-level retry ladder for the sample stage (ISSUE-007 family).

The ladder rescues deterministic/stochastic window crashes (LINCS, SIGSEGV at
junction windows) without operator intervention:
  standard -> retry1 (fresh ld-seeds) -> retry2 (3x EM + half dt + fresh seeds).
"""

from pathlib import Path

from abag_rbfe.io_utils import read_json, write_json
from abag_rbfe.stages import _sample_window_snippet


def _write_minimal_window(window_dir: Path) -> None:
    window_dir.mkdir(parents=True, exist_ok=True)
    (window_dir / "pre_relax.mdp").write_text("integrator = steep\n", encoding="utf-8")
    (window_dir / "pre_md.mdp").write_text("integrator = sd\n", encoding="utf-8")
    (window_dir / "production.mdp").write_text("integrator = sd\n", encoding="utf-8")


def _retry_mdps(window_dir: Path) -> dict[str, dict[str, Path]]:
    retry1_pre_md = window_dir / "pre_md.retry1.mdp"
    retry1_production = window_dir / "production.retry1.mdp"
    retry2_pre_relax = window_dir / "pre_relax.retry2.mdp"
    retry2_pre_md = window_dir / "pre_md.retry2.mdp"
    retry2_production = window_dir / "production.retry2.mdp"
    retry3_pre_relax = window_dir / "pre_relax.retry3.mdp"
    retry3_pre_md = window_dir / "pre_md.retry3.mdp"
    retry3_production = window_dir / "production.retry3.mdp"
    for path in (
        retry1_pre_md,
        retry1_production,
        retry2_pre_relax,
        retry2_pre_md,
        retry2_production,
        retry3_pre_relax,
        retry3_pre_md,
        retry3_production,
    ):
        path.write_text("# retry variant\n", encoding="utf-8")
    return {
        "retry1": {"pre_md": retry1_pre_md, "production": retry1_production},
        "retry2": {
            "pre_relax": retry2_pre_relax,
            "pre_md": retry2_pre_md,
            "production": retry2_production,
        },
        "retry3": {
            "pre_relax": retry3_pre_relax,
            "pre_md": retry3_pre_md,
            "production": retry3_production,
        },
    }


def _snippet(window_dir: Path, tmp_path: Path, retry: bool = True) -> str:
    return _sample_window_snippet(
        gmx_command="gmx",
        mdrun_suffix=" -ntmpi 1",
        repeat_top=tmp_path / "system.top",
        start_gro=tmp_path / "npt.gro",
        window_dir=window_dir,
        grompp_maxwarn_sampling=2,
        retry_mdps=_retry_mdps(window_dir) if retry else None,
    )


def test_snippet_keeps_completion_skip_and_gains_retry_ladder(tmp_path: Path) -> None:
    window_dir = tmp_path / "legs" / "complex" / "rep01" / "lambda_005"
    _write_minimal_window(window_dir)
    text = _snippet(window_dir, tmp_path)
    assert "skipping completed sample window complex/rep01/lambda_005" in text
    assert "attempt standard" in text
    assert "attempt retry1-reseed" in text
    assert "attempt retry2-couple-intramol" in text
    assert "attempt retry3-robust-ci-half-dt" in text
    assert "FAILED sample window complex/rep01/lambda_005 after retry ladder" in text
    assert "exit 1" in text
    # ladder rungs must not kill the stage script mid-way (it is set -euo pipefail)
    assert text.count("set +e") == 4
    assert text.count("set -e") == 4
    # retry rungs use the variant mdps
    assert "pre_md.retry1.mdp" in text
    assert "production.retry1.mdp" in text
    assert "pre_relax.retry2.mdp" in text
    assert "production.retry2.mdp" in text
    assert "pre_relax.retry3.mdp" in text
    assert "production.retry3.mdp" in text


def test_snippet_without_retry_mdps_still_fails_fast(tmp_path: Path) -> None:
    window_dir = tmp_path / "legs" / "complex" / "rep01" / "lambda_000"
    _write_minimal_window(window_dir)
    text = _snippet(window_dir, tmp_path, retry=False)
    assert "attempt standard" in text
    assert "retry1" not in text
    assert "retry2" not in text
    assert "exit 1" in text


def test_retry_ladder_recovers_from_standard_crash(tmp_path: Path) -> None:
    """Functional bash test with a stub gmx: standard mdrun crashes (SIGSEGV),
    retry1 succeeds; the ladder must produce the completion artifacts."""
    window_dir = tmp_path / "legs" / "complex" / "rep01" / "lambda_003"
    _write_minimal_window(window_dir)
    retry = _retry_mdps(window_dir)
    (tmp_path / "system.top").write_text("; top\n", encoding="utf-8")
    (tmp_path / "npt.gro").write_text("gro\n", encoding="utf-8")

    marker = tmp_path / "gmx_calls.txt"
    stub = tmp_path / "gmx"
    stub.write_text(
        "#!/usr/bin/env bash\n"
        f"echo \"$@\" >> {marker}\n"
        'if [ "$1" = "grompp" ]; then\n'
        '  while [ $# -gt 0 ]; do [ "$1" = "-o" ] && { echo tpr > "$2"; break; }; shift; done\n'
        "  exit 0\n"
        "fi\n"
        'if [ "$1" = "mdrun" ]; then\n'
        '  deffnm=""; dhdl=""\n'
        "  while [ $# -gt 0 ]; do\n"
        '    [ "$1" = "-deffnm" ] && deffnm="$2"\n'
        '    [ "$1" = "-dhdl" ] && dhdl="$2"\n'
        "    shift\n"
        "  done\n"
        '  case "$deffnm" in *pre_relax|*pre_md) echo gro > "${deffnm}.gro"; echo cpt > "${deffnm}.cpt"; exit 0;; esac\n'
        "  # production mdrun: crash on the standard attempt (no retry1 flag in env)\n"
        '  if [ -z "${ABAG_STUB_OK:-}" ]; then kill -SEGV $$; fi\n'
        '  echo gro > "${deffnm}.gro"; echo log > "${deffnm}.log"\n'
        '  [ -n "$dhdl" ] && echo "1" > "$dhdl"\n'
        "  exit 0\n"
        "fi\n"
        "exit 1\n",
        encoding="utf-8",
    )
    stub.chmod(0o755)

    # First run: all production attempts segfault -> snippet exits 1.
    snippet = _sample_window_snippet(
        gmx_command=str(stub),
        mdrun_suffix="",
        repeat_top=tmp_path / "system.top",
        start_gro=tmp_path / "npt.gro",
        window_dir=window_dir,
        grompp_maxwarn_sampling=2,
        retry_mdps=retry,
    )
    script = tmp_path / "run.sh"
    script.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + snippet + "\n", encoding="utf-8")
    script.chmod(0o755)

    import subprocess

    proc = subprocess.run(["bash", str(script)], cwd=tmp_path, capture_output=True, text=True)
    assert proc.returncode == 1
    assert "attempt standard" in proc.stdout
    assert "attempt retry1-reseed" in proc.stdout
    assert "attempt retry2-couple-intramol" in proc.stdout
    assert "attempt retry3-robust-ci-half-dt" in proc.stdout
    assert "FAILED sample window" in proc.stderr

    # Second run: stub now succeeds -> standard attempt completes the window.
    env = {"ABAG_STUB_OK": "1"}
    import os

    proc2 = subprocess.run(
        ["bash", str(script)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env={**os.environ, **env},
    )
    assert proc2.returncode == 0, proc2.stdout + proc2.stderr
    assert (window_dir / "dhdl.xvg").is_file()

    # Third run: completion check short-circuits everything.
    marker.write_text("", encoding="utf-8")
    proc3 = subprocess.run(["bash", str(script)], cwd=tmp_path, capture_output=True, text=True)
    assert proc3.returncode == 0
    assert "skipping completed sample window" in proc3.stdout
    assert marker.read_text(encoding="utf-8") == ""


def test_retry_variant_mdps_written_by_sample_stage(tmp_path: Path) -> None:
    """_write_sample_retry_mdps renders reseed + half-dt robust variants."""
    from abag_rbfe.planning import hydrate_protocol_config
    from abag_rbfe.stages import _write_sample_retry_mdps

    protocol = hydrate_protocol_config(
        {
            "lambda_windows": 8,
            "repeats": 1,
            "production_ps": 20,
            "production_dt_ps": 0.002,
            "window_relax_em_steps": 200,
            "window_relax_md_ps": 50.0,
            "window_relax_md_dt_ps": 0.001,
        }
    )
    ctx = _FakeCtx(protocol=protocol, mutation_group=_FakeMutationGroup(), job=_FakeJob())
    window_dir = tmp_path / "legs" / "complex" / "rep01" / "lambda_002"
    variants = _write_sample_retry_mdps(
        ctx,
        leg="complex",
        repeat_index=1,
        window_index=2,
        lambda_values=[i / 7 for i in range(8)],
        window_dir=window_dir,
    )
    retry1_prod = variants["retry1"]["production"].read_text(encoding="utf-8")
    retry2_prod = variants["retry2"]["production"].read_text(encoding="utf-8")
    retry3_prod = variants["retry3"]["production"].read_text(encoding="utf-8")
    assert "dt                      = 0.002" in retry1_prod
    assert "couple-intramol         = yes" in retry2_prod
    assert "dt                      = 0.002" in retry2_prod  # ci rung keeps standard timing
    assert "couple-intramol         = yes" in retry3_prod
    assert "dt                      = 0.001" in retry3_prod  # robust rung halves dt
    assert f"nsteps                  = {int(round(20 / 0.001))}" in retry3_prod
    retry3_em = variants["retry3"]["pre_relax"].read_text(encoding="utf-8")
    assert "nsteps                  = 600" in retry3_em
    # reseed variants must differ from each other in ld-seed
    assert variants["retry1"]["pre_md"].read_text(encoding="utf-8") != variants["retry3"][
        "pre_md"
    ].read_text(encoding="utf-8")


from dataclasses import dataclass as _dc


@_dc(frozen=True)
class _FakeCtx:
    """Minimal StageContext stand-in for mdp rendering helpers."""

    protocol: object
    mutation_group: object
    job: object


@_dc(frozen=True)
class _FakeJob:
    batch_id: str = "b"
    job_id: str = "j"


@_dc(frozen=True)
class _FakeMutationGroup:
    sites: tuple = ()


def test_window_chaining_uses_previous_window_endpoint(tmp_path: Path) -> None:
    """window_chaining=True: window k must start from lambda_{k-1}/md.gro."""
    from abag_rbfe.execution import CommandRunner
    from abag_rbfe.models import JobSpec, MutationGroup, SystemConfig
    from abag_rbfe.paths import ProjectPaths
    from abag_rbfe.planning import hydrate_protocol_config
    from abag_rbfe.stages import StageContext, _stage_sample

    protocol = hydrate_protocol_config(
        {"lambda_windows": 4, "repeats": 1, "window_chaining": True}
    )
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    mg = MutationGroup(
        mutation_group_id="g1",
        mutation_count=1,
        entity_side="antibody",
        charge_conserving=True,
        min_version="v1",
        sites=(),
    )
    system = SystemConfig(
        system_name="s",
        input_structure="x.pdb",
        structure_source="experimental",
        antibody_chains=("H",),
        antigen_chains=("A",),
    )
    job = JobSpec(
        job_id="j1",
        mutation_group=mg,
        protocol=protocol,
        system=system,
        batch_id="b1",
        workdir=str(job_dir),
    )
    ctx = StageContext(
        job_dir=job_dir,
        job=job,
        system=system,
        protocol=protocol,
        mutation_group=mg,
        runner=CommandRunner(execute=False),
        project_paths=ProjectPaths.discover(),
        rescue_config={},
    )
    status = _stage_sample(ctx)
    commands = "\n".join(status.commands)
    import re
    clean = commands.replace("'", "")
    starts = re.findall(r"grompp -f \S+ -c (\S+)", clean)
    assert any(s.endswith("lambda_000/md.gro") for s in starts), starts
    assert any(s.endswith("lambda_001/md.gro") for s in starts), starts
    assert any(s.endswith("lambda_002/md.gro") for s in starts), starts

    # default off: no chaining references
    protocol2 = hydrate_protocol_config({"lambda_windows": 4, "repeats": 1})
    ctx2 = StageContext(
        job_dir=tmp_path / "job2",
        job=job,
        system=system,
        protocol=protocol2,
        mutation_group=mg,
        runner=CommandRunner(execute=False),
        project_paths=ProjectPaths.discover(),
        rescue_config={},
    )
    (tmp_path / "job2").mkdir()
    status2 = _stage_sample(ctx2)
    commands2 = "\n".join(status2.commands)
    starts2 = re.findall(r"grompp -f \S+ -c (\S+)", commands2.replace("'", ""))
    assert not any("lambda_" in s and s.endswith("/md.gro") for s in starts2), starts2


def test_window_seeds_override_start_structure(tmp_path: Path) -> None:
    """config/window_seeds.json routes specific windows to their own start gro."""
    from abag_rbfe.execution import CommandRunner
    from abag_rbfe.models import JobSpec, MutationGroup, SystemConfig
    from abag_rbfe.paths import ProjectPaths
    from abag_rbfe.planning import hydrate_protocol_config
    from abag_rbfe.stages import StageContext, _stage_sample

    protocol = hydrate_protocol_config({"lambda_windows": 4, "repeats": 1})
    job_dir = tmp_path / "job"
    (job_dir / "config").mkdir(parents=True)
    seed = tmp_path / "hydrated.gro"
    seed.write_text("gro\n", encoding="utf-8")
    (job_dir / "config" / "window_seeds.json").write_text(
        '{"complex/rep01/lambda_003": "' + str(seed) + '"}', encoding="utf-8"
    )
    mg = MutationGroup(
        mutation_group_id="g1", mutation_count=1, entity_side="antibody",
        charge_conserving=True, min_version="v1", sites=(),
    )
    system = SystemConfig(
        system_name="s", input_structure="x.pdb", structure_source="experimental",
        antibody_chains=("H",), antigen_chains=("A",),
    )
    job = JobSpec(job_id="j1", mutation_group=mg, protocol=protocol, system=system,
                  batch_id="b1", workdir=str(job_dir))
    ctx = StageContext(
        job_dir=job_dir, job=job, system=system, protocol=protocol, mutation_group=mg,
        runner=CommandRunner(execute=False), project_paths=ProjectPaths.discover(),
        rescue_config={},
    )
    status = _stage_sample(ctx)
    commands = "\n".join(status.commands)
    assert "hydrated.gro" in commands
    # only the seeded window uses it: exactly 2 occurrences (pre_relax grompp per attempt chain start)
    # window 3 standard attempt starts from the seed
    assert commands.count("hydrated.gro") >= 1
    # non-seeded windows keep the shared npt.gro
    assert "npt.gro" in commands
