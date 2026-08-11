from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path


SCRIPT = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
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
    report_script = root / "benchmarks" / "ab_bind" / "report_calibrated_validation.py"
    report_script.parent.mkdir(parents=True, exist_ok=True)
    report_script.write_text("print('stub report')\n", encoding="utf-8")
    return root, python_log


def test_refresh_wrapper_throttles_successful_runs(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    summary_output = plan_root / "reports" / "calibrated_validation_summary.json"
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
            f"-u {root}/benchmarks/ab_bind/report_calibrated_validation.py "
            f"--plan-root {plan_root} --summary-output {summary_output}"
        )
    ]
    stamp_path = plan_root / "reports" / "watch" / "calibrated_validation_summary.last_run"
    assert stamp_path.exists()


def test_refresh_wrapper_swallows_report_failures_by_default(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=2)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    summary_output = plan_root / "reports" / "calibrated_validation_summary.json"
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
            f"-u {root}/benchmarks/ab_bind/report_calibrated_validation.py "
            f"--plan-root {plan_root} --summary-output {summary_output}"
        )
    ]
    log_path = plan_root / "reports" / "watch" / "calibrated_validation_summary_refresh.log"
    log_text = log_path.read_text(encoding="utf-8")
    assert "failed rc=2" in log_text


def test_refresh_wrapper_honors_plan_root_env_override(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    report_plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    summary_output = report_plan_root / "reports" / "calibrated_validation_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "PLAN_ROOT": str(plan_root),
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
            f"-u {root}/benchmarks/ab_bind/report_calibrated_validation.py "
            f"--plan-root {report_plan_root} --summary-output {summary_output}"
        )
    ]
    assert (plan_root / "reports" / "watch" / "calibrated_validation_summary.last_run").exists()


def test_refresh_wrapper_honors_report_plan_root_override(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    report_plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_rescues"
    summary_output = report_plan_root / "reports" / "calibrated_validation_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "PLAN_ROOT": str(plan_root),
        "REPORT_PLAN_ROOT": str(report_plan_root),
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
            f"-u {root}/benchmarks/ab_bind/report_calibrated_validation.py "
            f"--plan-root {report_plan_root} --summary-output {summary_output}"
        )
    ]


def test_refresh_wrapper_propagates_merged_extra_plan_roots(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    report_plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_rescues"
    extra_root_1 = root / "runs" / "benchmarks" / "abbind_core_v1_validation_sampling_qc_rescues_verify_preeq_20260611_1"
    extra_root_2 = root / "runs" / "benchmarks" / "abbind_core_v1_validation_deep_rescues"
    summary_output = report_plan_root / "reports" / "calibrated_validation_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "PLAN_ROOT": str(plan_root),
        "REPORT_PLAN_ROOT": str(report_plan_root),
        "MERGED_EXTRA_PLAN_ROOTS": os.pathsep.join([str(extra_root_1), str(report_plan_root), "", str(extra_root_2)]),
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
            f"-u {root}/benchmarks/ab_bind/report_calibrated_validation.py "
            f"--plan-root {report_plan_root} --summary-output {summary_output} "
            f"--extra-plan-root {extra_root_1} --extra-plan-root {extra_root_2}"
        )
    ]


def test_refresh_wrapper_honors_explicit_summary_output_override(tmp_path: Path) -> None:
    root, python_log = _fake_root(tmp_path, exit_code=0)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    report_plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_rescues"
    summary_output = plan_root / "reports" / "custom" / "calibrated_validation_summary.json"
    env = {
        **os.environ,
        "ABAG_RBFE_ROOT": str(root),
        "FAKE_PYTHON_LOG": str(python_log),
        "PLAN_ROOT": str(plan_root),
        "REPORT_PLAN_ROOT": str(report_plan_root),
        "SUMMARY_OUTPUT": str(summary_output),
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
            f"-u {root}/benchmarks/ab_bind/report_calibrated_validation.py "
            f"--plan-root {report_plan_root} --summary-output {summary_output}"
        )
    ]
