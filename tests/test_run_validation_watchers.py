from __future__ import annotations

import json
import os
import stat
import subprocess
import time
from pathlib import Path


SCRIPT = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_watchers.sh")
ROBUST_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_robust_watcher.sh"
)
PRIORITY_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_priority_watcher.sh"
)
RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_rescue_watcher.sh"
)
TARGETED_REPEAT_SPREAD_RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_targeted_repeat_spread_rescue_watcher.sh"
)
TARGETED_LAMBDA_RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_targeted_lambda_rescue_watcher.sh"
)
SAMPLING_QC_RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_sampling_qc_rescue_watcher.sh"
)
DEEP_RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_deep_rescue_watcher.sh"
)
ULTRA_RESCUE_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_ultra_rescue_watcher.sh"
)
STALE_RECOVERY_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_stale_recovery_watcher.sh"
)
GAP_RECOVERY_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_gap_recovery_watcher.sh"
)
PRIORITY_BACKLOG_WRAPPER = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_validation_priority_backlog_watcher.sh"
)
VALIDATION_WATCHLIST_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_watchlists.py"
)
VALIDATION_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_rescues.sh"
)
VALIDATION_TARGETED_REPEAT_SPREAD_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_targeted_repeat_spread_rescues.sh"
)
VALIDATION_TARGETED_LAMBDA_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_targeted_lambda_rescues.sh"
)
VALIDATION_SAMPLING_QC_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_sampling_qc_rescues.sh"
)
VALIDATION_DEEP_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_deep_rescues.sh"
)
VALIDATION_ULTRA_RESCUE_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_ultra_rescues.sh"
)
CALIBRATED_VALIDATION_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"
)
VALIDATION_POST_REPORT_REFRESH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_post_report.sh"
)
PERSISTENT_WATCH = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/run_persistent_watch.sh"
)


def _write_executable(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _fake_root(tmp_path: Path) -> Path:
    root = tmp_path / "abag-rbfep"
    bench = root / "benchmarks" / "ab_bind"
    sleep_script = "#!/usr/bin/env bash\nexec >/dev/null 2>&1 </dev/null\nsleep 30\n"
    _write_executable(bench / "run_persistent_watch.sh", sleep_script)
    _write_executable(bench / "run_validation_robust_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_priority_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_rescue_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_targeted_repeat_spread_rescue_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_targeted_lambda_rescue_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_sampling_qc_rescue_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_deep_rescue_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_ultra_rescue_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_stale_recovery_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_gap_recovery_watcher.sh", sleep_script)
    _write_executable(bench / "run_validation_priority_backlog_watcher.sh", sleep_script)
    return root


def _run_watchers(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["ABAG_RBFE_ROOT"] = str(root)
    return subprocess.run(
        ["bash", str(SCRIPT), *args],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )


def _pidfile(root: Path, name: str) -> Path:
    reports = {
        "robust": "abbind_core_v1_validation_robust_plan",
        "priority": "abbind_core_v1_validation_priority_plan",
        "rescue": "abbind_core_v1_validation_priority_rescues",
        "targeted": "abbind_core_v1_validation_targeted_repeat_spread_rescues",
        "lambda": "abbind_core_v1_validation_targeted_lambda_rescues",
        "sampling-qc": "abbind_core_v1_validation_sampling_qc_rescues",
        "deep": "abbind_core_v1_validation_deep_rescues",
        "ultra": "abbind_core_v1_validation_ultra_rescues",
        "stale": "abbind_core_v1_validation_priority_plan",
        "gap": "abbind_core_v1_validation_priority_plan",
        "backlog": "abbind_core_v1_validation_priority_plan",
    }
    filenames = {
        "robust": "watch_validation_robust.pid",
        "priority": "watch_validation_priority.pid",
        "rescue": "watch_validation_rescue.pid",
        "targeted": "watch_validation_targeted_repeat_spread_rescue.pid",
        "lambda": "watch_validation_targeted_lambda_rescue.pid",
        "sampling-qc": "watch_validation_sampling_qc_rescue.pid",
        "deep": "watch_validation_deep_rescue.pid",
        "ultra": "watch_validation_ultra_rescue.pid",
        "stale": "watch_validation_stale_recovery.pid",
        "gap": "watch_validation_gap_recovery.pid",
        "backlog": "watch_validation_priority_backlog.pid",
    }
    return root / "runs" / "benchmarks" / reports[name] / "reports" / "watch" / filenames[name]


def _read_pid(root: Path, name: str) -> int:
    return int(_pidfile(root, name).read_text().strip())


def _spawn_dummy_priority_process(root: Path) -> subprocess.Popen[bytes]:
    process = subprocess.Popen(
        [
            "python3",
            "-c",
            "import time; time.sleep(30)",
            str(root / "benchmarks" / "ab_bind" / "run_validation_priority_watcher.sh"),
        ]
    )
    time.sleep(0.2)
    return process


def test_status_reports_all_watchers_by_default(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run_watchers(root, "status")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "robust\tstatus=stopped",
        "priority\tstatus=stopped",
        "rescue\tstatus=stopped",
        "targeted\tstatus=stopped",
        "lambda\tstatus=stopped",
        "sampling-qc\tstatus=stopped",
        "deep\tstatus=stopped",
        "ultra\tstatus=stopped",
        "stale\tstatus=stopped",
        "gap\tstatus=stopped",
        "backlog\tstatus=stopped",
    ]


def test_status_reports_priority_running_without_pidfile_when_matching_process_exists(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    process = _spawn_dummy_priority_process(root)

    try:
        result = _run_watchers(root, "status", "priority")

        assert result.returncode == 0
        assert result.stdout.strip() == f"priority\tpid={process.pid}\tstatus=running(no-pidfile)"
    finally:
        process.terminate()
        process.wait(timeout=5)


def test_start_priority_reuses_matching_process_and_writes_pidfile(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    process = _spawn_dummy_priority_process(root)

    try:
        start = _run_watchers(root, "start", "priority")

        assert start.returncode == 0
        assert start.stdout.strip() == f"priority already running pid={process.pid}"
        assert _read_pid(root, "priority") == process.pid
    finally:
        process.terminate()
        process.wait(timeout=5)


def test_start_defaults_to_robust_only(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    start = _run_watchers(root, "start")
    try:
        assert start.returncode == 0
        assert "started robust pid=" in start.stdout
        assert _pidfile(root, "robust").exists()
        assert not _pidfile(root, "priority").exists()
        assert not _pidfile(root, "rescue").exists()
        assert not _pidfile(root, "targeted").exists()
        assert not _pidfile(root, "lambda").exists()
        assert not _pidfile(root, "sampling-qc").exists()
        assert not _pidfile(root, "deep").exists()
        assert not _pidfile(root, "ultra").exists()
        assert not _pidfile(root, "stale").exists()
        assert not _pidfile(root, "gap").exists()
        assert not _pidfile(root, "backlog").exists()
    finally:
        stop = _run_watchers(root, "stop", "all")
        assert stop.returncode == 0


def test_start_all_requires_explicit_target(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    start = _run_watchers(root, "start", "all")
    try:
        assert start.returncode == 0
        assert _pidfile(root, "robust").exists()
        assert _pidfile(root, "priority").exists()
        assert _pidfile(root, "rescue").exists()
        assert _pidfile(root, "targeted").exists()
        assert _pidfile(root, "lambda").exists()
        assert _pidfile(root, "sampling-qc").exists()
        assert _pidfile(root, "deep").exists()
        assert _pidfile(root, "ultra").exists()
        assert _pidfile(root, "stale").exists()
        assert _pidfile(root, "gap").exists()
        assert _pidfile(root, "backlog").exists()
    finally:
        stop = _run_watchers(root, "stop", "all")
        assert stop.returncode == 0


def test_restart_defaults_to_robust_only(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    start = _run_watchers(root, "start", "all")
    assert start.returncode == 0
    initial_pids = {
        name: _read_pid(root, name)
        for name in ("robust", "priority", "rescue", "targeted", "lambda", "sampling-qc", "deep", "ultra", "stale", "gap", "backlog")
    }

    try:
        restart = _run_watchers(root, "restart")

        assert restart.returncode == 0
        assert _read_pid(root, "robust") != initial_pids["robust"]
        assert _read_pid(root, "priority") == initial_pids["priority"]
        assert _read_pid(root, "rescue") == initial_pids["rescue"]
        assert _read_pid(root, "targeted") == initial_pids["targeted"]
        assert _read_pid(root, "lambda") == initial_pids["lambda"]
        assert _read_pid(root, "sampling-qc") == initial_pids["sampling-qc"]
        assert _read_pid(root, "deep") == initial_pids["deep"]
        assert _read_pid(root, "ultra") == initial_pids["ultra"]
        assert _read_pid(root, "stale") == initial_pids["stale"]
        assert _read_pid(root, "gap") == initial_pids["gap"]
        assert _read_pid(root, "backlog") == initial_pids["backlog"]
    finally:
        stop = _run_watchers(root, "stop", "all")
        assert stop.returncode == 0


def test_restart_all_replaces_all_watcher_pidfiles(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    start = _run_watchers(root, "start", "all")
    assert start.returncode == 0
    initial_pids = {
        name: _read_pid(root, name)
        for name in ("robust", "priority", "rescue", "targeted", "lambda", "sampling-qc", "deep", "ultra", "stale", "gap", "backlog")
    }

    try:
        restart = _run_watchers(root, "restart", "all")

        assert restart.returncode == 0
        for name in ("robust", "priority", "rescue", "targeted", "lambda", "sampling-qc", "deep", "ultra", "stale", "gap", "backlog"):
            assert _read_pid(root, name) != initial_pids[name]
    finally:
        stop = _run_watchers(root, "stop", "all")
        assert stop.returncode == 0


def test_start_priority_propagates_env_appended_merged_extra_plan_roots(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    bench = root / "benchmarks" / "ab_bind"
    observed_env = root / "priority_watcher_env.txt"

    _write_executable(
        bench / "run_persistent_watch.sh",
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                "logfile=\"$1\"",
                "shift",
                "\"$@\" &",
                "child=$!",
                "wait \"$child\"",
            ]
        )
        + "\n",
    )
    _write_executable(
        bench / "run_validation_priority_watcher.sh",
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f"env | sort > {str(observed_env)!r}",
                "sleep 30",
            ]
        )
        + "\n",
    )

    extra_roots = ":".join(
        [
            str(root / "runs" / "benchmarks" / "extra_one"),
            str(root / "runs" / "benchmarks" / "extra_two"),
        ]
    )
    env = os.environ.copy()
    env["ABAG_RBFE_ROOT"] = str(root)
    env["MERGED_EXTRA_PLAN_ROOTS"] = extra_roots

    start = subprocess.run(
        ["bash", str(SCRIPT), "start", "priority"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    try:
        assert start.returncode == 0
        assert "started priority pid=" in start.stdout

        deadline = time.time() + 5
        while not observed_env.exists() and time.time() < deadline:
            time.sleep(0.1)

        assert observed_env.exists()
        recorded_env = observed_env.read_text(encoding="utf-8")
        assert f"MERGED_EXTRA_PLAN_ROOTS={extra_roots}" in recorded_env
    finally:
        stop = _run_watchers(root, "stop", "all")
        assert stop.returncode == 0


def test_robust_wrapper_default_allow_list_excludes_2nz9() -> None:
    text = ROBUST_WRAPPER.read_text()

    assert '"2nz9-antigen-a-h1064a"' not in text


def test_robust_wrapper_merges_rescue_root_into_priority_view() -> None:
    text = ROBUST_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-6}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${CALIBRATED_VALIDATION_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh}"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${VALIDATION_POST_REPORT_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh}"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"' in text
    assert 'ROBUST_PASS_OUTLIER_THRESHOLD="${ROBUST_PASS_OUTLIER_THRESHOLD:-5.0}"' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${PRIORITY_RUNS_ROOT}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'refresh_validation_watchlists.py' in text
    assert "--mode robust" in text
    assert '--robust-pass-outlier-threshold "${ROBUST_PASS_OUTLIER_THRESHOLD}"' in text
    assert "--once" in text
    assert '"1bj1-antigen-w-g92a"' in text
    assert '"1mlc-antibody-h-s57a"' in text
    assert '"1mlc-antibody-h-s57v"' in text
    assert '"1mlc-antibody-h-t31a"' in text
    assert '"1mlc-antibody-h-t31v"' in text
    assert '"1mlc-antibody-l-n92a"' in text
    assert '"3hfm-antigen-y-y20a"' in text
    assert '"3hfm-antibody-l-n31a"' not in text
    assert '"3hfm-antibody-l-n32a"' not in text
    assert 'RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_2}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_3}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_4}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text


def test_priority_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(PRIORITY_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Priority watcher is disabled by default" in result.stderr


def test_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Rescue watcher is disabled by default" in result.stderr


def test_targeted_repeat_spread_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(TARGETED_REPEAT_SPREAD_RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Targeted repeat-spread rescue watcher is disabled by default" in result.stderr


def test_targeted_lambda_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(TARGETED_LAMBDA_RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Targeted lambda rescue watcher is disabled by default" in result.stderr


def test_sampling_qc_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(SAMPLING_QC_RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Sampling/QC rescue watcher is disabled by default" in result.stderr


def test_deep_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(DEEP_RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Deep rescue watcher is disabled by default" in result.stderr


def test_ultra_rescue_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(ULTRA_RESCUE_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Ultra rescue watcher is disabled by default" in result.stderr


def test_stale_recovery_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(STALE_RECOVERY_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Stale recovery watcher is disabled by default" in result.stderr


def test_gap_recovery_wrapper_requires_explicit_enable() -> None:
    result = subprocess.run(
        ["bash", str(GAP_RECOVERY_WRAPPER)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 75
    assert "Gap recovery watcher is disabled by default" in result.stderr


def test_rescue_wrapper_allows_cross_root_same_job_ids() -> None:
    text = RESCUE_WRAPPER.read_text()

    assert "--allow-active-elsewhere-job-ids" in text


def test_rescue_wrapper_targets_current_warning_followups() -> None:
    text = RESCUE_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-6}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${CALIBRATED_VALIDATION_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh}"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${VALIDATION_POST_REPORT_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh}"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"' in text
    assert 'REFRESH_RESCUES_EACH_PASS="${REFRESH_RESCUES_EACH_PASS:-1}"' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${PRIORITY_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'refresh_validation_rescues.sh' in text
    assert 'refresh_validation_watchlists.py' in text
    assert "--mode hotspots" in text
    assert "--once" in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_3}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_4}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text
    assert '"1bj1-antigen-w-g92a"' in text
    assert '"1cz8-antigen-w-g92a"' in text
    assert '"3hfm-antibody-h-c95a"' in text
    assert '"3hfm-antibody-h-y33a"' in text
    assert '"3hfm-antibody-h-y50a"' in text
    assert '"3hfm-antigen-y-y20a"' in text


def test_validation_rescue_refresh_uses_validation_priority_source_defaults() -> None:
    text = VALIDATION_RESCUE_REFRESH.read_text()

    assert 'SOURCE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"' in text
    assert 'RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"' in text
    assert 'BATCH_PREFIX="${BATCH_PREFIX:-abbind-rescue}"' in text
    assert 'REPEAT_INCREMENT="${REPEAT_INCREMENT:-1}"' in text
    assert 'LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-4}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-2.0}"' in text
    assert 'WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-2.0}"' in text
    assert 'WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-2.0}"' in text
    assert 'NVT_SCALE="${NVT_SCALE:-2.0}"' in text
    assert 'NPT_SCALE="${NPT_SCALE:-2.0}"' in text
    assert 'REQUIRE_ACTIVE_ALTERNATE="${REQUIRE_ACTIVE_ALTERNATE:-1}"' in text
    assert "--window-relax-em-scale" in text
    assert "--window-relax-md-scale" in text
    assert "--nvt-scale" in text
    assert "--npt-scale" in text
    assert "--split-name validation" in text
    assert "rescue-abbind" in text


def test_validation_targeted_repeat_spread_refresh_targets_primary_leg_only() -> None:
    text = VALIDATION_TARGETED_REPEAT_SPREAD_RESCUE_REFRESH.read_text()

    assert 'SOURCE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"' in text
    assert 'RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'BATCH_PREFIX="${BATCH_PREFIX:-abbind-targeted-repeat-spread-rescue}"' in text
    assert 'REPEAT_INCREMENT="${REPEAT_INCREMENT:-0}"' in text
    assert 'LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-0}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-2.0}"' in text
    assert "--target-primary-repeat-spread-leg" in text
    assert "--require-target-primary-repeat-spread-leg" in text
    assert "--split-name validation" in text
    assert "rescue-abbind" in text


def test_validation_targeted_lambda_rescue_refresh_densifies_single_leg_path() -> None:
    text = VALIDATION_TARGETED_LAMBDA_RESCUE_REFRESH.read_text()

    assert (
        'SOURCE_PLAN_ROOT="${SOURCE_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues}"'
        in text
    )
    assert (
        'RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues}"'
        in text
    )
    assert 'BATCH_PREFIX="${BATCH_PREFIX:-abbind-targeted-lambda-rescue}"' in text
    assert 'REPEAT_INCREMENT="${REPEAT_INCREMENT:-0}"' in text
    assert 'LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-4}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-1.0}"' in text
    assert 'WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-1.0}"' in text
    assert 'WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-1.0}"' in text
    assert 'NVT_SCALE="${NVT_SCALE:-1.0}"' in text
    assert 'NPT_SCALE="${NPT_SCALE:-1.0}"' in text
    assert 'FORCE_LAMBDA_INCREMENT="${FORCE_LAMBDA_INCREMENT:-1}"' in text
    assert 'ALLOW_TARGETED_LEG_COUNT_DEEPENING="${ALLOW_TARGETED_LEG_COUNT_DEEPENING:-1}"' in text
    assert "--force-lambda-increment" in text
    assert "--allow-targeted-leg-count-deepening" in text
    assert "--target-primary-repeat-spread-leg" in text
    assert "--require-target-primary-repeat-spread-leg" in text
    assert "--split-name validation" in text
    assert "rescue-abbind" in text


def test_validation_sampling_qc_rescue_refresh_sources_sampling_hotspots_from_rescue_root() -> None:
    text = VALIDATION_SAMPLING_QC_RESCUE_REFRESH.read_text()

    assert (
        'SOURCE_PLAN_ROOT="${SOURCE_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues}"'
        in text
    )
    assert (
        'PRIORITY_PLAN_ROOT="${PRIORITY_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan}"'
        in text
    )
    assert (
        'TARGETED_REPEAT_SPREAD_PLAN_ROOT="${TARGETED_REPEAT_SPREAD_PLAN_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues}"'
        in text
    )
    assert (
        'RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues}"'
        in text
    )
    assert 'BATCH_PREFIX="${BATCH_PREFIX:-abbind-sampling-qc-rescue}"' in text
    assert 'REPEAT_INCREMENT="${REPEAT_INCREMENT:-1}"' in text
    assert 'LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-4}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-1.0}"' in text
    assert 'WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-3.0}"' in text
    assert 'WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-3.0}"' in text
    assert 'NVT_SCALE="${NVT_SCALE:-3.0}"' in text
    assert 'NPT_SCALE="${NPT_SCALE:-3.0}"' in text
    assert 'FORCE_LAMBDA_INCREMENT="${FORCE_LAMBDA_INCREMENT:-1}"' in text
    assert 'PREFER_ACTIVE_ALTERNATE_SOURCE="${PREFER_ACTIVE_ALTERNATE_SOURCE:-1}"' in text
    assert 'HOTSPOT_COMPLEX_IDS="${HOTSPOT_COMPLEX_IDS:-3HFM}"' in text
    assert 'HOTSPOT_REPEAT_INCREMENT="${HOTSPOT_REPEAT_INCREMENT:-1}"' in text
    assert 'HOTSPOT_LAMBDA_INCREMENT="${HOTSPOT_LAMBDA_INCREMENT:-6}"' in text
    assert 'HOTSPOT_PRODUCTION_SCALE="${HOTSPOT_PRODUCTION_SCALE:-1.0}"' in text
    assert 'HOTSPOT_WINDOW_RELAX_EM_SCALE="${HOTSPOT_WINDOW_RELAX_EM_SCALE:-4.0}"' in text
    assert 'HOTSPOT_WINDOW_RELAX_MD_SCALE="${HOTSPOT_WINDOW_RELAX_MD_SCALE:-4.0}"' in text
    assert 'HOTSPOT_NVT_SCALE="${HOTSPOT_NVT_SCALE:-4.0}"' in text
    assert 'HOTSPOT_NPT_SCALE="${HOTSPOT_NPT_SCALE:-4.0}"' in text
    assert "--prefer-active-alternate-source" in text
    assert "--force-lambda-increment" in text
    assert "--hotspot-complex-id" in text
    assert "--hotspot-repeat-increment" in text
    assert "--hotspot-lambda-increment" in text
    assert "--hotspot-production-scale" in text
    assert "--hotspot-window-relax-em-scale" in text
    assert "--hotspot-window-relax-md-scale" in text
    assert "--hotspot-nvt-scale" in text
    assert "--hotspot-npt-scale" in text
    assert "--target-primary-repeat-spread-leg" not in text
    assert "--split-name validation" in text
    assert "rescue-abbind" in text


def test_validation_deep_rescue_refresh_prefers_sampling_rich_alternates() -> None:
    text = VALIDATION_DEEP_RESCUE_REFRESH.read_text()

    assert 'SOURCE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"' in text
    assert 'TARGETED_REPEAT_SPREAD_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'ROBUST_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"' in text
    assert 'RESCUE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"' in text
    assert 'RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'BATCH_PREFIX="${BATCH_PREFIX:-abbind-deep-rescue}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-4.0}"' in text
    assert 'WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-2.0}"' in text
    assert 'WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-2.0}"' in text
    assert 'NVT_SCALE="${NVT_SCALE:-2.0}"' in text
    assert 'NPT_SCALE="${NPT_SCALE:-2.0}"' in text
    assert 'REQUIRE_ACTIVE_ALTERNATE="${REQUIRE_ACTIVE_ALTERNATE:-0}"' in text
    assert 'HOTSPOT_COMPLEX_IDS="${HOTSPOT_COMPLEX_IDS:-3HFM}"' in text
    assert 'HOTSPOT_REPEAT_INCREMENT="${HOTSPOT_REPEAT_INCREMENT:-2}"' in text
    assert 'HOTSPOT_LAMBDA_INCREMENT="${HOTSPOT_LAMBDA_INCREMENT:-6}"' in text
    assert 'HOTSPOT_PRODUCTION_SCALE="${HOTSPOT_PRODUCTION_SCALE:-6.0}"' in text
    assert 'HOTSPOT_WINDOW_RELAX_EM_SCALE="${HOTSPOT_WINDOW_RELAX_EM_SCALE:-3.0}"' in text
    assert 'HOTSPOT_WINDOW_RELAX_MD_SCALE="${HOTSPOT_WINDOW_RELAX_MD_SCALE:-3.0}"' in text
    assert 'HOTSPOT_NVT_SCALE="${HOTSPOT_NVT_SCALE:-3.0}"' in text
    assert 'HOTSPOT_NPT_SCALE="${HOTSPOT_NPT_SCALE:-3.0}"' in text
    assert "--window-relax-em-scale" in text
    assert "--window-relax-md-scale" in text
    assert "--nvt-scale" in text
    assert "--npt-scale" in text
    assert '--extra-plan-root "${SAMPLING_QC_PLAN_ROOT}"' in text
    assert "--hotspot-complex-id" in text
    assert "--hotspot-repeat-increment" in text
    assert "--hotspot-lambda-increment" in text
    assert "--hotspot-production-scale" in text
    assert "--hotspot-window-relax-em-scale" in text
    assert "--hotspot-window-relax-md-scale" in text
    assert "--hotspot-nvt-scale" in text
    assert "--hotspot-npt-scale" in text
    assert "--force-repeat-increment" in text
    assert "--prefer-active-alternate-source" in text
    assert "--target-primary-repeat-spread-leg" in text
    assert "--allow-targeted-leg-count-deepening" in text
    assert 'if [ "${REQUIRE_ACTIVE_ALTERNATE}" = "1" ]; then' in text
    assert 'CMD+=(--require-active-alternate)' in text


def test_validation_ultra_rescue_refresh_prefers_deepest_active_alternates() -> None:
    text = VALIDATION_ULTRA_RESCUE_REFRESH.read_text()

    assert 'SOURCE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_plan"' in text
    assert 'TARGETED_REPEAT_SPREAD_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'ROBUST_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_robust_plan"' in text
    assert 'RESCUE_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"' in text
    assert 'DEEP_PLAN_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"' in text
    assert 'BATCH_PREFIX="${BATCH_PREFIX:-abbind-ultra-rescue}"' in text
    assert 'HOTSPOT_COMPLEX_IDS="${HOTSPOT_COMPLEX_IDS:-3HFM}"' in text
    assert 'HOTSPOT_REPEAT_INCREMENT="${HOTSPOT_REPEAT_INCREMENT:-3}"' in text
    assert 'HOTSPOT_LAMBDA_INCREMENT="${HOTSPOT_LAMBDA_INCREMENT:-12}"' in text
    assert 'HOTSPOT_PRODUCTION_SCALE="${HOTSPOT_PRODUCTION_SCALE:-3.0}"' in text
    assert 'HOTSPOT_WINDOW_RELAX_EM_SCALE="${HOTSPOT_WINDOW_RELAX_EM_SCALE:-3.0}"' in text
    assert 'HOTSPOT_WINDOW_RELAX_MD_SCALE="${HOTSPOT_WINDOW_RELAX_MD_SCALE:-3.0}"' in text
    assert 'HOTSPOT_NVT_SCALE="${HOTSPOT_NVT_SCALE:-3.0}"' in text
    assert 'HOTSPOT_NPT_SCALE="${HOTSPOT_NPT_SCALE:-3.0}"' in text
    assert 'ALLOW_PASS_QC_OUTLIER_RESCUE="${ALLOW_PASS_QC_OUTLIER_RESCUE:-0}"' in text
    assert 'REPEAT_INCREMENT="${REPEAT_INCREMENT:-2}"' in text
    assert 'LAMBDA_INCREMENT="${LAMBDA_INCREMENT:-8}"' in text
    assert 'PRODUCTION_SCALE="${PRODUCTION_SCALE:-2.0}"' in text
    assert 'WINDOW_RELAX_EM_SCALE="${WINDOW_RELAX_EM_SCALE:-2.0}"' in text
    assert 'WINDOW_RELAX_MD_SCALE="${WINDOW_RELAX_MD_SCALE:-2.0}"' in text
    assert 'NVT_SCALE="${NVT_SCALE:-2.0}"' in text
    assert 'NPT_SCALE="${NPT_SCALE:-2.0}"' in text
    assert "--window-relax-em-scale" in text
    assert "--window-relax-md-scale" in text
    assert "--nvt-scale" in text
    assert "--npt-scale" in text
    assert '--extra-plan-root "${SAMPLING_QC_PLAN_ROOT}"' in text
    assert "--hotspot-complex-id" in text
    assert "--hotspot-repeat-increment" in text
    assert "--hotspot-lambda-increment" in text
    assert "--hotspot-production-scale" in text
    assert "--hotspot-window-relax-em-scale" in text
    assert "--hotspot-window-relax-md-scale" in text
    assert "--hotspot-nvt-scale" in text
    assert "--hotspot-npt-scale" in text
    assert "--force-repeat-increment" in text
    assert "--prefer-active-alternate-source" in text
    assert 'if [ "${REQUIRE_ACTIVE_ALTERNATE}" = "1" ]; then' in text
    assert 'CMD+=(--require-active-alternate)' in text
    assert 'if [ "${ALLOW_PASS_QC_OUTLIER_RESCUE}" = "1" ]; then' in text
    assert 'CMD+=(--allow-pass-qc-outlier-rescue)' in text
    assert "--target-primary-repeat-spread-leg" in text
    assert "--allow-targeted-leg-count-deepening" in text


def test_validation_watchlist_refresh_has_default_proactive_1mlc_jobs() -> None:
    text = VALIDATION_WATCHLIST_REFRESH.read_text()

    assert "DEFAULT_PROACTIVE_ROBUST_JOB_IDS" in text
    assert '"1mlc-antibody-h-s57a"' in text
    assert '"1mlc-antibody-h-s57v"' in text
    assert '"1mlc-antibody-h-t31a"' in text
    assert '"1mlc-antibody-h-t31v"' in text
    assert '"1mlc-antibody-l-n92a"' in text
    assert '"hotspots"' in text
    assert '"stale"' in text
    assert '"gap"' in text
    assert '"backlog"' in text
    assert '"ultra"' in text


def test_priority_backlog_wrapper_targets_non_overlapping_partial_backlog() -> None:
    text = PRIORITY_BACKLOG_WRAPPER.read_text()

    assert "--only-listed" in text
    assert '--watch-tag backlog' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10500}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-96}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"' in text
    assert 'WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_backlog_watchlist_refresh.json}"' in text
    assert '--mode backlog' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"' in text
    assert '--candidate-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}"' in text
    assert '--merged-extra-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}"' in text
    assert '--merged-extra-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}"' in text
    assert '--merged-extra-plan-root "${RESCUE_RUNS_ROOT}"' in text
    assert '--merged-extra-plan-root "${DEEP_RESCUE_RUNS_ROOT}"' in text
    assert '--merged-extra-plan-root "${ULTRA_RESCUE_RUNS_ROOT}"' in text
    assert '"2nz9-antigen-a-f953a"' in text
    assert '"2nz9-antigen-a-t1063a"' in text
    assert '"2nz9-antigen-a-l919a"' in text
    assert '"2nz9-antigen-a-n918a"' in text
    assert '"3hfm-antibody-l-y50a"' in text
    assert '"3hfm-antibody-l-y50l"' in text
    assert '"3hfm-antibody-h-c95f"' in text
    assert '"1cz8-antigen-w-g88a"' in text


def test_priority_wrapper_merges_rescue_root_into_priority_view() -> None:
    text = PRIORITY_WRAPPER.read_text()

    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_priority_rescues"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_2}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_3}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_4}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text


def test_priority_wrapper_accepts_env_appended_merged_extra_plan_roots() -> None:
    text = PRIORITY_WRAPPER.read_text()

    assert 'MERGED_EXTRA_PLAN_ROOTS="${MERGED_EXTRA_PLAN_ROOTS:-}"' in text
    assert "IFS=':' read -r -a MERGED_EXTRA_PLAN_ROOT_VALUES <<< \"${MERGED_EXTRA_PLAN_ROOTS}\"" in text
    assert 'merged_extra_plan_root="${merged_extra_plan_root//[[:space:]]/}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${merged_extra_plan_root}")' in text


def test_validation_watchers_merge_targeted_lambda_root_into_merged_reports() -> None:
    expected = {
        ROBUST_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        PRIORITY_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        RESCUE_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        TARGETED_REPEAT_SPREAD_RESCUE_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        SAMPLING_QC_RESCUE_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${TARGETED_LAMBDA_RESCUE_RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues}"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        DEEP_RESCUE_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        ULTRA_RESCUE_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        STALE_RECOVERY_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        GAP_RECOVERY_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}}"',
            'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")',
        ],
        PRIORITY_BACKLOG_WRAPPER: [
            'TARGETED_LAMBDA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues"',
            '--merged-extra-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}"',
        ],
    }

    for wrapper, snippets in expected.items():
        text = wrapper.read_text()
        for snippet in snippets:
            assert snippet in text, f"{wrapper.name} missing snippet: {snippet}"


def test_duplicate_tolerant_validation_rescue_wrappers_cap_active_copies_per_job_id() -> None:
    expected = {
        RESCUE_WRAPPER: [
            'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"',
            '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"',
        ],
        TARGETED_REPEAT_SPREAD_RESCUE_WRAPPER: [
            'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"',
            '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"',
        ],
        TARGETED_LAMBDA_RESCUE_WRAPPER: [
            'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"',
            '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"',
        ],
        SAMPLING_QC_RESCUE_WRAPPER: [
            'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"',
            '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"',
        ],
        DEEP_RESCUE_WRAPPER: [
            'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"',
            '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"',
        ],
        ULTRA_RESCUE_WRAPPER: [
            'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-3}"',
            '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"',
        ],
    }

    for wrapper, snippets in expected.items():
        text = wrapper.read_text()
        for snippet in snippets:
            assert snippet in text, f"{wrapper.name} missing snippet: {snippet}"


def test_priority_wrapper_appends_env_merged_extra_plan_roots_to_watch_args(tmp_path: Path) -> None:
    root = tmp_path / "abag-rbfep"
    bench = root / "benchmarks" / "ab_bind"
    observed_args = root / "observed_watch_args.json"

    _write_executable(
        bench / "watch_validation_priority.py",
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json, sys",
                "from pathlib import Path",
                f"Path({str(observed_args)!r}).write_text(json.dumps(sys.argv[1:]), encoding='utf-8')",
            ]
        )
        + "\n",
    )
    _write_executable(bench / "refresh_validation_post_report.sh", "#!/usr/bin/env bash\nexit 0\n")

    wrapper_path = bench / "run_validation_priority_watcher.sh"
    wrapper_text = PRIORITY_WRAPPER.read_text(encoding="utf-8").replace(
        'ROOT="/mnt/data/liuchao/abag-rbfep"',
        f'ROOT="{root}"',
    )
    _write_executable(wrapper_path, wrapper_text)

    extra_root_one = root / "runs" / "benchmarks" / "extra_one"
    extra_root_two = root / "runs" / "benchmarks" / "extra_two"

    env = os.environ.copy()
    env["ABAG_ENABLE_PRIORITY_WATCHER"] = "1"
    env["MERGED_EXTRA_PLAN_ROOTS"] = f"{extra_root_one}: {extra_root_two}::"

    result = subprocess.run(
        ["bash", str(wrapper_path)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0
    observed = json.loads(observed_args.read_text(encoding="utf-8"))
    merged_extra_plan_roots = [
        observed[index + 1]
        for index, arg in enumerate(observed[:-1])
        if arg == "--merged-extra-plan-root"
    ]
    assert str(extra_root_one) in merged_extra_plan_roots
    assert str(extra_root_two) in merged_extra_plan_roots


def test_queue_view_watchers_scan_targeted_lambda_alternates() -> None:
    expected = {
        STALE_RECOVERY_WRAPPER: [
            '--candidate-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}"',
            'WATCH_ARGS+=(--extra-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}")',
        ],
        GAP_RECOVERY_WRAPPER: [
            '--candidate-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}"',
            'WATCH_ARGS+=(--extra-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}")',
        ],
        PRIORITY_BACKLOG_WRAPPER: [
            '--candidate-plan-root "${TARGETED_LAMBDA_RESCUE_RUNS_ROOT}"',
        ],
    }

    for wrapper, snippets in expected.items():
        text = wrapper.read_text()
        for snippet in snippets:
            assert snippet in text, f"{wrapper.name} missing snippet: {snippet}"


def test_watchers_script_uses_distinct_priority_and_stale_patterns() -> None:
    text = SCRIPT.read_text()

    assert 'PRIORITY_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_priority_watcher.sh"' in text
    assert 'STALE_RECOVERY_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_stale_recovery_watcher.sh"' in text
    assert 'GAP_RECOVERY_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_gap_recovery_watcher.sh"' in text
    assert 'BACKLOG_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_priority_backlog_watcher.sh"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_targeted_repeat_spread_rescue_watcher.sh"' in text
    assert 'TARGETED_LAMBDA_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_targeted_lambda_rescue_watcher.sh"' in text
    assert 'SAMPLING_QC_RESCUE_PATTERN="${ROOT}/benchmarks/ab_bind/run_validation_sampling_qc_rescue_watcher.sh"' in text
    assert 'stop_one "priority" "${PRIORITY_PIDFILE}" "${PRIORITY_PATTERN}"' in text
    assert 'status_line "priority" "${PRIORITY_PIDFILE}" "${PRIORITY_PATTERN}"' in text
    assert 'stop_one "targeted" "${TARGETED_REPEAT_SPREAD_RESCUE_PIDFILE}" "${TARGETED_REPEAT_SPREAD_RESCUE_PATTERN}"' in text
    assert 'status_line "targeted" "${TARGETED_REPEAT_SPREAD_RESCUE_PIDFILE}" "${TARGETED_REPEAT_SPREAD_RESCUE_PATTERN}"' in text
    assert 'stop_one "lambda" "${TARGETED_LAMBDA_RESCUE_PIDFILE}" "${TARGETED_LAMBDA_RESCUE_PATTERN}"' in text
    assert 'status_line "lambda" "${TARGETED_LAMBDA_RESCUE_PIDFILE}" "${TARGETED_LAMBDA_RESCUE_PATTERN}"' in text
    assert 'stop_one "sampling-qc" "${SAMPLING_QC_RESCUE_PIDFILE}" "${SAMPLING_QC_RESCUE_PATTERN}"' in text
    assert 'status_line "sampling-qc" "${SAMPLING_QC_RESCUE_PIDFILE}" "${SAMPLING_QC_RESCUE_PATTERN}"' in text
    assert 'stop_one "stale" "${STALE_RECOVERY_PIDFILE}" "${STALE_RECOVERY_PATTERN}"' in text
    assert 'status_line "stale" "${STALE_RECOVERY_PIDFILE}" "${STALE_RECOVERY_PATTERN}"' in text
    assert 'stop_one "gap" "${GAP_RECOVERY_PIDFILE}" "${GAP_RECOVERY_PATTERN}"' in text
    assert 'status_line "gap" "${GAP_RECOVERY_PIDFILE}" "${GAP_RECOVERY_PATTERN}"' in text
    assert 'stop_one "backlog" "${BACKLOG_PIDFILE}" "${BACKLOG_PATTERN}"' in text
    assert 'status_line "backlog" "${BACKLOG_PIDFILE}" "${BACKLOG_PATTERN}"' in text


def test_targeted_repeat_spread_rescue_wrapper_targets_single_leg_repeat_spread_path() -> None:
    text = TARGETED_REPEAT_SPREAD_RESCUE_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-4}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-24}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"' in text
    assert 'REFRESH_TARGETED_RESCUES_EACH_PASS="${REFRESH_TARGETED_RESCUES_EACH_PASS:-1}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'TARGETED_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD="${TARGETED_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD:-3.0}"' in text
    assert 'MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PRIORITY_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${RESCUE_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${PRIORITY_RUNS_ROOT}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'refresh_validation_targeted_repeat_spread_rescues.sh' in text
    assert 'reports/plan_jobs.csv' in text
    assert 'append_materialized_job_ids' in text
    assert 'latest_stage_state' in text
    assert '"running", "stale_running"' in text
    assert "--mode targeted" in text
    assert '--targeted-no-active-alt-abs-error-threshold "${TARGETED_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD}"' in text
    assert "--allow-active-elsewhere-job-ids" in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text
    assert '"1bj1-antigen-w-g88a"' in text
    assert '"3hfm-antibody-h-c95a"' in text
    assert '"3nps-antigen-a-h138a"' in text


def test_targeted_lambda_rescue_wrapper_targets_single_leg_lambda_densification_path() -> None:
    text = TARGETED_LAMBDA_RESCUE_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-4}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-24}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"' in text
    assert 'REFRESH_TARGETED_LAMBDA_RESCUES_EACH_PASS="${REFRESH_TARGETED_LAMBDA_RESCUES_EACH_PASS:-1}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert (
        'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${PRIORITY_RUNS_ROOT}}"'
        in text
    )
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PRIORITY_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_7="${MERGED_EXTRA_PLAN_ROOT_7:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'refresh_validation_targeted_lambda_rescues.sh' in text
    assert 'reports/plan_jobs.csv' in text
    assert 'append_materialized_job_ids' in text
    assert 'latest_stage_state' in text
    assert '"running", "stale_running"' in text
    assert "--mode targeted" in text
    assert "--allow-active-elsewhere-job-ids" in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_7}")' in text
    assert '"1bj1-antigen-w-g88a"' in text
    assert '"3hfm-antibody-h-c95a"' in text
    assert '"3nps-antigen-a-h138a"' in text


def test_sampling_qc_rescue_wrapper_targets_two_leg_sampling_qc_path() -> None:
    text = SAMPLING_QC_RESCUE_WRAPPER.read_text()

    assert 'ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"' in text
    assert 'ABAG_RBFE="${ABAG_RBFE:-${ROOT}/.venv/bin/abag-rbfe}"' in text
    assert 'RUNS_ROOT="${RUNS_ROOT:-${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues}"' in text
    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-4}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-24}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"' in text
    assert 'SAMPLING_QC_COMPLEX_IDS="${SAMPLING_QC_COMPLEX_IDS:-}"' in text
    assert 'SAMPLING_QC_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD="${SAMPLING_QC_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD:-2.5}"' in text
    assert 'REFRESH_SAMPLING_QC_RESCUES_EACH_PASS="${REFRESH_SAMPLING_QC_RESCUES_EACH_PASS:-1}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${CALIBRATED_VALIDATION_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh}"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${VALIDATION_POST_REPORT_REFRESH:-${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh}"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PRIORITY_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${RESCUE_RUNS_ROOT}:${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${PRIORITY_RUNS_ROOT}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'refresh_validation_sampling_qc_rescues.sh' in text
    assert "report-abbind" in text
    assert 'reports/plan_jobs.csv' in text
    assert 'refresh_validation_watchlists.py' in text
    assert 'latest_stage_state' in text
    assert '"running", "stale_running"' in text
    assert "--mode sampling_qc" in text
    assert "--sampling-qc-complex-id" in text
    assert '--sampling-qc-no-active-alt-abs-error-threshold "${SAMPLING_QC_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD}"' in text
    assert "--allow-active-elsewhere-job-ids" in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text
    assert '"1mlc-antibody-l-n92a"' in text
    assert '"3hfm-antibody-h-c95a"' in text
    assert '"3hfm-antibody-h-y33a"' in text
    assert '"3hfm-antibody-h-y50a"' in text
    assert '"3hfm-antibody-l-n32a"' in text
    assert '"3hfm-antigen-y-y20a"' in text


def test_sampling_qc_rescue_wrapper_only_appends_active_materialized_jobs(tmp_path: Path) -> None:
    root = tmp_path / "abag-rbfep"
    bench = root / "benchmarks" / "ab_bind"
    runs_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_sampling_qc_rescues"
    reports_dir = runs_root / "reports"
    watch_dir = reports_dir / "watch"
    watch_dir.mkdir(parents=True, exist_ok=True)

    (reports_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                "job_id,latest_stage_state",
                "active-job,running",
                "stale-job,stale_running",
                "finished-job,completed",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (reports_dir / "rescue_candidates.csv").write_text(
        "\n".join(
            [
                "source_job_id",
                "active-job",
                "finished-job",
                "legacy-only-job",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    observed_args = root / "observed_watch_args.json"
    _write_executable(
        bench / "refresh_validation_watchlists.py",
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json, sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "if '--output-json' in args:",
                "    output = Path(args[args.index('--output-json') + 1])",
                "    output.parent.mkdir(parents=True, exist_ok=True)",
                "    output.write_text(json.dumps({'ok': True}), encoding='utf-8')",
                "print('fresh-job')",
            ]
        )
        + "\n",
    )
    _write_executable(
        bench / "watch_validation_priority.py",
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json, sys",
                "from pathlib import Path",
                f"Path({str(observed_args)!r}).write_text(json.dumps(sys.argv[1:]), encoding='utf-8')",
            ]
        )
        + "\n",
    )
    _write_executable(bench / "refresh_validation_sampling_qc_rescues.sh", "#!/usr/bin/env bash\nexit 0\n")
    _write_executable(bench / "refresh_validation_post_report.sh", "#!/usr/bin/env bash\nexit 0\n")

    wrapper_path = bench / "run_validation_sampling_qc_rescue_watcher.sh"
    wrapper_text = SAMPLING_QC_RESCUE_WRAPPER.read_text(encoding="utf-8").replace(
        'ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"',
        f'ROOT="{root}"',
    )
    _write_executable(wrapper_path, wrapper_text)

    env = os.environ.copy()
    env["ABAG_ENABLE_SAMPLING_QC_RESCUE_WATCHER"] = "1"
    env["WATCH_ONCE"] = "1"
    env["REFRESH_SAMPLING_QC_RESCUES_EACH_PASS"] = "0"
    env["POST_REFRESH_COMMAND"] = "/bin/true"

    result = subprocess.run(
        ["bash", str(wrapper_path)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0
    observed = json.loads(observed_args.read_text(encoding="utf-8"))
    assert observed[-3:] == ["fresh-job", "active-job", "stale-job"]
    assert "finished-job" not in observed
    assert "legacy-only-job" not in observed


def test_sampling_qc_rescue_wrapper_bootstraps_missing_plan_reports(tmp_path: Path) -> None:
    root = tmp_path / "abag-rbfep"
    bench = root / "benchmarks" / "ab_bind"
    runs_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_sampling_qc_rescues"
    watch_dir = runs_root / "reports" / "watch"
    watch_dir.mkdir(parents=True, exist_ok=True)

    observed_args = root / "observed_watch_args.json"
    report_bootstrap_sentinel = root / "report_bootstrap.txt"

    _write_executable(
        root / ".venv" / "bin" / "abag-rbfe",
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                "if [ \"$1\" = \"batch\" ] && [ \"$2\" = \"report-abbind\" ]; then",
                "  shift 2",
                "  plan_root=\"\"",
                "  while [ \"$#\" -gt 0 ]; do",
                "    case \"$1\" in",
                "      --plan-root)",
                "        shift",
                "        plan_root=\"$1\"",
                "        ;;",
                "    esac",
                "    shift",
                "  done",
                "  mkdir -p \"${plan_root}/reports\"",
                "  printf 'job_id,latest_stage_state\\nfresh-job,not_started\\n' > \"${plan_root}/reports/plan_jobs.csv\"",
                "  printf 'batch_id\\nfresh-batch\\n' > \"${plan_root}/reports/plan_batches.csv\"",
                f"  printf 'bootstrapped\\n' > {str(report_bootstrap_sentinel)!r}",
                "  exit 0",
                "fi",
                "exit 1",
            ]
        )
        + "\n",
    )
    _write_executable(
        bench / "watch_validation_priority.py",
        "\n".join(
            [
                "#!/usr/bin/env python3",
                "import json, sys",
                "from pathlib import Path",
                f"Path({str(observed_args)!r}).write_text(json.dumps(sys.argv[1:]), encoding='utf-8')",
            ]
        )
        + "\n",
    )
    _write_executable(bench / "refresh_validation_sampling_qc_rescues.sh", "#!/usr/bin/env bash\nexit 0\n")
    _write_executable(bench / "refresh_validation_post_report.sh", "#!/usr/bin/env bash\nexit 0\n")

    wrapper_path = bench / "run_validation_sampling_qc_rescue_watcher.sh"
    wrapper_text = SAMPLING_QC_RESCUE_WRAPPER.read_text(encoding="utf-8").replace(
        'ROOT="${ABAG_RBFE_ROOT:-/mnt/data/liuchao/abag-rbfep}"',
        f'ROOT="{root}"',
    )
    _write_executable(wrapper_path, wrapper_text)

    env = os.environ.copy()
    env["ABAG_ENABLE_SAMPLING_QC_RESCUE_WATCHER"] = "1"
    env["WATCH_ONCE"] = "1"
    env["REFRESH_SAMPLING_QC_RESCUES_EACH_PASS"] = "0"
    env["POST_REFRESH_COMMAND"] = "/bin/true"

    result = subprocess.run(
        ["bash", str(wrapper_path), "fresh-job"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0
    assert report_bootstrap_sentinel.read_text(encoding="utf-8").strip() == "bootstrapped"
    observed = json.loads(observed_args.read_text(encoding="utf-8"))
    assert observed[-1] == "fresh-job"


def test_deep_rescue_wrapper_targets_hotspot_deepening_path() -> None:
    text = DEEP_RESCUE_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-4}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-24}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"' in text
    assert 'REFRESH_DEEP_RESCUES_EACH_PASS="${REFRESH_DEEP_RESCUES_EACH_PASS:-1}"' in text
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert 'WATCH_ARGS+=(--post-refresh-command "${POST_REFRESH_COMMAND}")' in text
    assert 'MERGED_PLAN_ROOT="${MERGED_PLAN_ROOT:-${PRIORITY_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${RUNS_ROOT}}"' in text
    assert 'ULTRA_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_ultra_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${PRIORITY_RUNS_ROOT}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'refresh_validation_deep_rescues.sh' in text
    assert "--mode hotspots" in text
    assert "--allow-active-elsewhere-job-ids" in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_4}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text
    assert '"3hfm-antibody-l-n31a"' in text
    assert '"3hfm-antigen-y-y20a"' in text


def test_ultra_rescue_wrapper_targets_high_impact_complexes() -> None:
    text = ULTRA_RESCUE_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-2}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-12000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_LOAD_PER_CORE="${MAX_LOAD_PER_CORE:-0}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-20}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"' in text
    assert 'REFRESH_ULTRA_RESCUES_EACH_PASS="${REFRESH_ULTRA_RESCUES_EACH_PASS:-1}"' in text
    assert 'ULTRA_PEARSON_GAIN_THRESHOLD="${ULTRA_PEARSON_GAIN_THRESHOLD:-0.2}"' in text
    assert 'ULTRA_ABS_ERROR_THRESHOLD="${ULTRA_ABS_ERROR_THRESHOLD:-5.0}"' in text
    assert 'ULTRA_PASS_OUTLIER_THRESHOLD="${ULTRA_PASS_OUTLIER_THRESHOLD:-${ULTRA_ABS_ERROR_THRESHOLD}}"' in text
    assert 'ULTRA_PASS_OUTLIER_ALLOW_INACTIVE_ALTERNATE="${ULTRA_PASS_OUTLIER_ALLOW_INACTIVE_ALTERNATE:-1}"' in text
    assert 'DEEP_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_deep_rescues"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${DEEP_RESCUE_RUNS_ROOT}:${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${ROBUST_RUNS_ROOT}:${PRIORITY_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'refresh_validation_ultra_rescues.sh' in text
    assert '--mode ultra' in text
    assert '--ultra-abs-error-threshold "${ULTRA_ABS_ERROR_THRESHOLD}"' in text
    assert '--ultra-pass-outlier-threshold "${ULTRA_PASS_OUTLIER_THRESHOLD}"' in text
    assert 'ULTRA_PASS_OUTLIER_JOB_IDS=()' in text
    assert 'payload.get("ultra_pass_outlier_job_ids", [])' in text
    assert 'ALLOW_PASS_QC_OUTLIER_RESCUE=1 REQUIRE_ACTIVE_ALTERNATE=0 "${ROOT}/benchmarks/ab_bind/refresh_validation_ultra_rescues.sh" "${ULTRA_PASS_OUTLIER_JOB_IDS[@]}" >/dev/null' in text
    assert '--allow-active-elsewhere-job-ids' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_5}")' in text
    assert 'WATCH_ARGS+=(--merged-extra-plan-root "${MERGED_EXTRA_PLAN_ROOT_6}")' in text
    assert '"3hfm-antibody-h-y33a"' in text
    assert '"3hfm-antibody-l-n32a"' in text


def test_stale_recovery_wrapper_targets_stale_validation_jobs() -> None:
    text = STALE_RECOVERY_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-14}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-100}"' in text
    assert 'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-2}"' in text
    assert 'THREAD_BUDGET_PLAN_ROOT_1="${THREAD_BUDGET_PLAN_ROOT_1:-${RUNS_ROOT}}"' in text
    assert 'THREAD_BUDGET_PLAN_ROOT_2="${THREAD_BUDGET_PLAN_ROOT_2:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-2}"' in text
    assert 'REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"' in text
    assert 'WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_stale_watchlist_refresh.json}"' in text
    assert 'INCLUDE_GAP_JOBS="${INCLUDE_GAP_JOBS:-0}"' in text
    assert 'GAP_WATCHLIST_REFRESH_JSON="${GAP_WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_gap_watchlist_refresh.json}"' in text
    assert (
        'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${ROBUST_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${DEEP_RESCUE_RUNS_ROOT}:${ULTRA_RESCUE_RUNS_ROOT}}"'
        in text
    )
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'refresh_validation_watchlists.py' in text
    assert '--mode stale' in text
    assert '--mode gap' in text
    assert '--candidate-plan-root "${ROBUST_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${DEEP_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${ULTRA_RESCUE_RUNS_ROOT}"' in text
    assert '--watch-tag stale' in text
    assert '--allow-active-elsewhere-job-ids' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"' in text
    assert 'WATCH_ARGS+=(--thread-budget-plan-root "${THREAD_BUDGET_PLAN_ROOT_1}")' in text
    assert 'WATCH_ARGS+=(--thread-budget-plan-root "${THREAD_BUDGET_PLAN_ROOT_2}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${ROBUST_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${DEEP_RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${ULTRA_RESCUE_RUNS_ROOT}")' in text
    assert '"1bj1-antigen-v-f17a"' in text
    assert '"1mlc-antibody-h-t31a"' in text


def test_gap_recovery_wrapper_targets_missing_complex_jobs() -> None:
    text = GAP_RECOVERY_WRAPPER.read_text()

    assert 'MAX_COMPUTE_APPS_PER_GPU="${MAX_COMPUTE_APPS_PER_GPU:-12}"' in text
    assert 'MIN_FREE_GPU_MEMORY_MB="${MIN_FREE_GPU_MEMORY_MB:-10000}"' in text
    assert 'MAX_GPU_UTILIZATION="${MAX_GPU_UTILIZATION:-60}"' in text
    assert 'MAX_ACTIVE_MDRUN_THREADS="${MAX_ACTIVE_MDRUN_THREADS:-100}"' in text
    assert 'MAX_ACTIVE_COPIES_PER_JOB_ID="${MAX_ACTIVE_COPIES_PER_JOB_ID:-2}"' in text
    assert 'THREAD_BUDGET_PLAN_ROOT_1="${THREAD_BUDGET_PLAN_ROOT_1:-${RUNS_ROOT}}"' in text
    assert 'THREAD_BUDGET_PLAN_ROOT_2="${THREAD_BUDGET_PLAN_ROOT_2:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MAX_LAUNCHES_PER_PASS="${MAX_LAUNCHES_PER_PASS:-1}"' in text
    assert 'REFRESH_WATCHLIST_EACH_PASS="${REFRESH_WATCHLIST_EACH_PASS:-1}"' in text
    assert 'WATCHLIST_REFRESH_JSON="${WATCHLIST_REFRESH_JSON:-${RUNS_ROOT}/reports/watch/validation_gap_watchlist_refresh.json}"' in text
    assert (
        'ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS="${ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS:-${ROBUST_RUNS_ROOT}:${RESCUE_RUNS_ROOT}:${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}:${DEEP_RESCUE_RUNS_ROOT}:${ULTRA_RESCUE_RUNS_ROOT}}"'
        in text
    )
    assert 'CALIBRATED_VALIDATION_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_calibrated_validation_summary.sh"' in text
    assert 'VALIDATION_POST_REPORT_REFRESH="${ROOT}/benchmarks/ab_bind/refresh_validation_post_report.sh"' in text
    assert 'POST_REFRESH_COMMAND="${POST_REFRESH_COMMAND:-${VALIDATION_POST_REPORT_REFRESH}}"' in text
    assert "export ABAG_RBFE_EQUILIBRATION_SEED_PLAN_ROOTS" in text
    assert 'MERGED_EXTRA_PLAN_ROOT_1="${MERGED_EXTRA_PLAN_ROOT_1:-${ROBUST_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_2="${MERGED_EXTRA_PLAN_ROOT_2:-${RESCUE_RUNS_ROOT}}"' in text
    assert 'TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_3="${MERGED_EXTRA_PLAN_ROOT_3:-${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_4="${MERGED_EXTRA_PLAN_ROOT_4:-${DEEP_RESCUE_RUNS_ROOT}}"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_5="${MERGED_EXTRA_PLAN_ROOT_5:-${ULTRA_RESCUE_RUNS_ROOT}}"' in text
    assert 'SAMPLING_QC_RESCUE_RUNS_ROOT="${ROOT}/runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues"' in text
    assert 'MERGED_EXTRA_PLAN_ROOT_6="${MERGED_EXTRA_PLAN_ROOT_6:-${SAMPLING_QC_RESCUE_RUNS_ROOT}}"' in text
    assert 'refresh_validation_watchlists.py' in text
    assert '--mode gap' in text
    assert '--candidate-plan-root "${ROBUST_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${DEEP_RESCUE_RUNS_ROOT}"' in text
    assert '--candidate-plan-root "${ULTRA_RESCUE_RUNS_ROOT}"' in text
    assert '--watch-tag gap' in text
    assert '--allow-active-elsewhere-job-ids' in text
    assert '--min-free-gpu-memory-mb "${MIN_FREE_GPU_MEMORY_MB}"' in text
    assert '--max-gpu-utilization "${MAX_GPU_UTILIZATION}"' in text
    assert '--max-active-copies-per-job-id "${MAX_ACTIVE_COPIES_PER_JOB_ID}"' in text
    assert 'WATCH_ARGS+=(--thread-budget-plan-root "${THREAD_BUDGET_PLAN_ROOT_1}")' in text
    assert 'WATCH_ARGS+=(--thread-budget-plan-root "${THREAD_BUDGET_PLAN_ROOT_2}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${ROBUST_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${TARGETED_REPEAT_SPREAD_RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${SAMPLING_QC_RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${DEEP_RESCUE_RUNS_ROOT}")' in text
    assert 'WATCH_ARGS+=(--extra-plan-root "${ULTRA_RESCUE_RUNS_ROOT}")' in text
    assert '"1mlc-antibody-h-s57a"' in text
    assert '"1mlc-antibody-l-n92a"' in text
    assert '"2nz9-antigen-a-f953a"' in text
    assert '"2nz9-antigen-a-t1063a"' in text


def test_persistent_watch_exits_without_restart_on_disable_code(tmp_path: Path) -> None:
    logfile = tmp_path / "watch.log"
    child = tmp_path / "disabled_child.sh"
    _write_executable(
        child,
        "#!/usr/bin/env bash\n"
        "echo 'disabled child'\n"
        "exit 75\n",
    )

    result = subprocess.run(
        ["bash", str(PERSISTENT_WATCH), str(logfile), str(child)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    log_text = logfile.read_text()
    assert "child exited code=75" in log_text
    assert "not restarting after disable code=75" in log_text
