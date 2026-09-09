#!/usr/bin/env python
"""Chaining + long-window pilot driver for r21a (2026-09-07).

Generates three rep-specific sample scripts (chained starts, retry ladder
included) so the 3 reps run in parallel on 3 GPUs (the stock sample.sh is
sequential). Requires the job's protocol already carrying
window_chaining=true + the target production_ps (build_legs must have been
regenerated with it).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from abag_rbfe.stages import (
    _job_legs,
    _lambda_values,
    _load_context,
    _mdrun_suffix,
    _resolve_gmx_command,
    _sample_window_snippet,
    _stage_env,
    _write_sample_retry_mdps,
)
from abag_rbfe.execution import CommandRunner


def main() -> int:
    job_dir = Path(sys.argv[1]).resolve()
    ctx = _load_context(job_dir, execute=False)
    if not getattr(ctx.protocol, "window_chaining", False):
        raise SystemExit("protocol.window_chaining is not true for this job")
    gmx_command = _resolve_gmx_command(ctx)
    mdrun_suffix = _mdrun_suffix(ctx)
    env = _stage_env(ctx)
    lambda_values = _lambda_values(ctx.protocol.lambda_windows)

    out_dir = job_dir / "artifacts" / "commands"
    out_dir.mkdir(parents=True, exist_ok=True)
    runner = CommandRunner(execute=False)

    for leg in _job_legs(ctx):
        for repeat_index in range(1, ctx.protocol.repeats + 1):
            repeat_dir = job_dir / "legs" / leg / f"rep{repeat_index:02d}"
            repeat_top = repeat_dir / "system.top"
            npt_gro = repeat_dir / "equilibration" / "npt.gro"
            commands: list[str] = []
            start_gro = npt_gro
            for window_index, _ in enumerate(lambda_values):
                window_dir = repeat_dir / f"lambda_{window_index:03d}"
                if window_index > 0:
                    start_gro = repeat_dir / f"lambda_{window_index - 1:03d}" / "md.gro"
                retry_mdps = _write_sample_retry_mdps(
                    ctx,
                    leg=leg,
                    repeat_index=repeat_index,
                    window_index=window_index,
                    lambda_values=lambda_values,
                    window_dir=window_dir,
                )
                commands.append(
                    _sample_window_snippet(
                        gmx_command=gmx_command,
                        mdrun_suffix=mdrun_suffix,
                        repeat_top=repeat_top,
                        start_gro=start_gro,
                        window_dir=window_dir,
                        grompp_maxwarn_sampling=ctx.protocol.grompp_maxwarn_sampling,
                        retry_mdps=retry_mdps,
                    )
                )
            script = out_dir / f"sample_{leg}_rep{repeat_index:02d}.sh"
            runner.write_script(script, commands, job_dir, env=env)
            print(f"wrote {script}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
