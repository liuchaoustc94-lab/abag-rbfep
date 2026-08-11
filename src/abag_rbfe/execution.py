"""External command planning and optional execution."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from functools import lru_cache
from dataclasses import dataclass
from pathlib import Path

from abag_rbfe.io_utils import ensure_dir


@dataclass(frozen=True)
class CommandOutcome:
    state: str
    message: str


class CommandRunner:
    """Write shell command scripts and optionally execute them."""

    def __init__(self, execute: bool):
        self.execute = execute

    def write_script(
        self,
        path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> None:
        ensure_dir(path.parent)
        lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            f"cd {shlex.quote(str(workdir))}",
            "",
        ]
        if env:
            for key in sorted(env):
                lines.append(f"export {key}={shlex.quote(env[key])}")
            lines.append("")
        lines.extend(commands or [":"])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        path.chmod(0o755)

    def ensure_binaries(self, commands: list[str]) -> list[str]:
        missing: list[str] = []
        for command in commands:
            # Multi-line shell snippets and control-flow blocks are executed by
            # bash directly; preflight token inspection is too lossy there.
            if "\n" in command or any(token in command for token in {"{", "}", "if ", "fi", "then", "else"}):
                continue
            parts = shlex.split(command)
            if not parts:
                continue
            binary = parts[0]
            if binary in {"bash", "python", "python3"}:
                continue
            if shutil.which(binary) is None:
                missing.append(binary)
        return sorted(set(missing))

    def run_script(
        self,
        script_path: Path,
        commands: list[str],
        workdir: Path,
        env: dict[str, str] | None = None,
    ) -> CommandOutcome:
        self.write_script(script_path, commands, workdir, env=env)
        if not commands:
            return CommandOutcome("completed", "No external commands were required.")
        if not self.execute:
            return CommandOutcome("planned", "Stage commands written; external execution not requested.")
        missing = self.ensure_binaries(commands)
        if missing:
            return CommandOutcome("blocked_external", f"Missing external binaries: {', '.join(missing)}")
        process_env = os.environ.copy()
        if env:
            process_env.update(env)
        log_path = script_path.with_suffix(".log")
        with log_path.open("w", encoding="utf-8") as log_handle:
            process = subprocess.Popen(
                ["/usr/bin/env", "bash", str(script_path)],
                cwd=workdir,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                env=process_env,
            )
            return_code = process.wait()
        if return_code != 0:
            if return_code < 0:
                signal_number = -return_code
                return CommandOutcome(
                    "failed",
                    f"External command terminated by signal {signal_number}. See {log_path}",
                )
            if return_code >= 128:
                signal_number = return_code - 128
                return CommandOutcome(
                    "failed",
                    f"External command terminated by signal {signal_number}. See {log_path}",
                )
            return CommandOutcome(
                "failed",
                f"External command failed with exit code {return_code}. See {log_path}",
            )
        return CommandOutcome("completed", f"External commands completed. Log: {log_path}")


@lru_cache(maxsize=1)
def discover_visible_gpu_devices() -> tuple[str, ...]:
    override = os.environ.get("ABAG_RBFE_VISIBLE_GPUS", "").strip()
    if override:
        devices = tuple(token.strip() for token in override.split(",") if token.strip())
        if devices:
            return devices

    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi is None:
        return ()

    result = subprocess.run(
        [nvidia_smi, "-L"],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return ()

    devices: list[str] = []
    for raw_line in result.stdout.splitlines():
        line = raw_line.strip()
        if not line.startswith("GPU "):
            continue
        prefix = line.split(":", 1)[0]
        parts = prefix.split()
        if len(parts) < 2:
            continue
        index = parts[1]
        if index.isdigit():
            devices.append(index)
    return tuple(devices)
