from __future__ import annotations

import subprocess
from pathlib import Path


CALIBRATION_MANAGER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_calibration_watchers.sh"
)
CALIBRATION_QUICK_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_calibration_quick_watcher.sh"
)
CALIBRATION_RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_calibration_rescue_watcher.sh"
)
CALIBRATION_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_calibration_rescues.sh"
)
CALIBRATED_VALIDATION_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
)


def test_calibration_quick_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(CALIBRATION_QUICK_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Calibration quick watcher is disabled by default" in result.stderr


def test_calibration_manager_starts_quick_with_explicit_enable_env() -> None:
    text = CALIBRATION_MANAGER.read_text()

    assert "env ABAG_ENABLE_QUICK_WATCHER=1" in text
    assert '"${ROOT}/benchmarks/ab_bind/run_calibration_quick_watcher.sh"' in text


def test_calibration_quick_wrapper_uses_relaxed_gpu_slot_default() -> None:
    text = CALIBRATION_QUICK_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-6}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${CALIBRATED_VALIDATION_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text


def test_calibration_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(CALIBRATION_RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Calibration rescue watcher is disabled by default" in result.stderr


def test_calibration_manager_starts_rescue_with_explicit_enable_env() -> None:
    text = CALIBRATION_MANAGER.read_text()

    assert "env ABAG_ENABLE_CALIBRATION_RESCUE_WATCHER=1" in text
    assert '"${ROOT}/benchmarks/ab_bind/run_calibration_rescue_watcher.sh"' in text


def test_calibration_rescue_refresh_uses_priority_strength_defaults() -> None:
    text = CALIBRATION_RESCUE_REFRESH.read_text()

    assert 'REPEAT_INCREMENT="${REPEAT_INCREMENT:-2}"' in text
    assert 'LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-6}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-20.0}"' in text
    assert 'WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-2.0}"' in text
    assert 'WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-2.0}"' in text
    assert 'NVT_SCALE="${NVT_SCALE:-2.0}"' in text
    assert 'NPT_SCALE="${NPT_SCALE:-2.0}"' in text
    assert "--window-relax-em-scale" in text
    assert "--window-relax-md-scale" in text
    assert "--nvt-scale" in text
    assert "--npt-scale" in text
    assert "rescue-abbind" in text
    assert "--force-repeat-increment" in text


def test_calibration_rescue_wrapper_refreshes_rescues_every_pass() -> None:
    text = CALIBRATION_RESCUE_WRAPPER.read_text()

    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${CALIBRATED_VALIDATION_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${QUICK_RUNS_ROOT}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'REFRESH_RESCUES_EACH_PASS="${REFRESH_RESCUES_EACH_PASS:-1}"' in text
    assert '"${ROOT}/benchmarks/ab_bind/refresh_calibration_rescues.sh" "${JOB_IDS[@]}" >/dev/null' in text
    assert "--once" in text
