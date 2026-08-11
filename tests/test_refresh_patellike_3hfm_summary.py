from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path


SCRIPT = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/patel_2021_3hfm/refresh_patellike_3hfm_summary.sh"
)


def _write_executable(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _fake_root(tmp_path: Path, *, exit_code: int = 0) -> tuple[Path, Path]:
    root = tmp_path / "abag-rbfep"
    python_log = tmp_path / "python.log"
    fake_python = root / ".venv" / "bin" / "python"
    _write_executable(
        fake_python,
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "printf '%s\\n' \"$*\" >> \"${FAKE_PYTHON_LOG}\"\n"
        f"exit {exit_code}\n",
    )
    report_script = root / "benchmarks" / "patel_2021_3hfm" / "report_patellike_3hfm.py"
    report_script.parent.mkdir(parents=True, exist_ok=True)
    report_script.write_text("print('stub report')\n", encoding="utf-8")
    return root, python_log


def test_refresh_wrapper_throttles_successful_runs(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    batch_dir = root / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    summary_output = batch_dir / "reports" / "patel_2021_3hfm_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "MIN_INTERVAL_SECONDS": "3600",
    }

    first = subprocess.run(
        ["bash", str(SCRIPT)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    second = subprocess.run(
        ["bash", str(SCRIPT)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert first.returncode == 0
    assert second.returncode == 0
    assert python_log.read_text(encoding="utf-8").splitlines() == [
        (
            f"-u {root}/benchmarks/patel_2021_3hfm/report_patellike_3hfm.py "
            f"--batch-dir {batch_dir} --summary-output {summary_output}"
        )
    ]
    stamp_path = batch_dir / "reports" / "watch" / "patellike_3hfm_summary.last_run"
    assert stamp_path.exists()


def test_refresh_wrapper_swallows_report_failures_by_default(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=2)
    batch_dir = root / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    summary_output = batch_dir / "reports" / "patel_2021_3hfm_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "MIN_INTERVAL_SECONDS": "0",
    }

    result = subprocess.run(
        ["bash", str(SCRIPT)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0
    assert python_log.read_text(encoding="utf-8").splitlines() == [
        (
            f"-u {root}/benchmarks/patel_2021_3hfm/report_patellike_3hfm.py "
            f"--batch-dir {batch_dir} --summary-output {summary_output}"
        )
    ]
    log_path = batch_dir / "reports" / "watch" / "patellike_3hfm_summary_refresh.log"
    log_text = log_path.read_text(encoding="utf-8")
    assert "failed rc=2" in log_text


def test_refresh_wrapper_honors_batch_dir_override(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    batch_dir = root / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference_custom"
    summary_output = batch_dir / "reports" / "patel_2021_3hfm_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "BATCH_DIR": str(batch_dir),
        "MIN_INTERVAL_SECONDS": "0",
    }

    result = subprocess.run(
        ["bash", str(SCRIPT)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0
    assert python_log.read_text(encoding="utf-8").splitlines() == [
        (
            f"-u {root}/benchmarks/patel_2021_3hfm/report_patellike_3hfm.py "
            f"--batch-dir {batch_dir} --summary-output {summary_output}"
        )
    ]
    assert (batch_dir / "reports" / "watch" / "patellike_3hfm_summary.last_run").exists()

