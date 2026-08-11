from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace


def _load_watch_module():
    module_path = (
        Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/watch_validation_priority.py")
    )
    spec = importlib.util.spec_from_file_location("watch_validation_priority", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_available_gpu_devices_prefers_lower_loaded_gpus(monkeypatch) -> None:
    module = _load_watch_module()
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 2, "1": 3, "2": 0, "3": 1})

    available = module.available_gpu_devices(
        ["0", "1", "2", "3"],
        gpu_counts=module.gpu_compute_counts(),
        max_compute_apps_per_gpu=3,
    )

    assert available == ["2", "3", "0"]


def test_available_gpu_devices_respects_zero_threshold(monkeypatch) -> None:
    module = _load_watch_module()
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0, "1": 0})

    available = module.available_gpu_devices(
        ["0", "1"],
        gpu_counts=module.gpu_compute_counts(),
        max_compute_apps_per_gpu=0,
    )

    assert available == []


def test_gpu_device_stats_parses_nvidia_smi_output(monkeypatch) -> None:
    module = _load_watch_module()
    monkeypatch.setattr(module, "shutil_which", lambda binary: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(
            returncode=0,
            stdout="\n".join(
                [
                    "0, 7274 MiB, 24564 MiB, 44 %, 3 %",
                    "1, 8837 MiB, 24564 MiB, 47 %, 3 %",
                ]
            ),
        ),
    )

    stats = module.gpu_device_stats()

    assert stats == {
        "0": {
            "memory_used_mb": 7274,
            "memory_total_mb": 24564,
            "gpu_utilization_percent": 44,
            "memory_utilization_percent": 3,
        },
        "1": {
            "memory_used_mb": 8837,
            "memory_total_mb": 24564,
            "gpu_utilization_percent": 47,
            "memory_utilization_percent": 3,
        },
    }


def test_available_gpu_devices_can_use_gpu_headroom_override_when_count_gate_is_full() -> None:
    module = _load_watch_module()

    available = module.available_gpu_devices(
        ["0", "1"],
        gpu_counts={"0": 14, "1": 14},
        max_compute_apps_per_gpu=4,
        gpu_stats={
            "0": {
                "memory_used_mb": 7274,
                "memory_total_mb": 24564,
                "gpu_utilization_percent": 44,
            },
            "1": {
                "memory_used_mb": 15064,
                "memory_total_mb": 24564,
                "gpu_utilization_percent": 44,
            },
        },
        min_free_gpu_memory_mb=12000,
        max_gpu_utilization=60,
    )

    assert available == ["0"]


def test_launch_allowed_by_cpu_blocks_when_load_is_high(monkeypatch) -> None:
    module = _load_watch_module()
    monkeypatch.setattr(module.os, "cpu_count", lambda: 8)
    monkeypatch.setattr(module.os, "getloadavg", lambda: (8.0, 7.0, 6.0))

    allowed, load_per_core = module.launch_allowed_by_cpu(max_load_per_core=0.9)

    assert allowed is False
    assert load_per_core == 1.0


def test_active_canonical_job_copy_counts_deduplicates_multiple_processes_per_copy(monkeypatch) -> None:
    module = _load_watch_module()
    ps_stdout = "\n".join(
        [
            "liuchao 101 1 0 00:00 ? 00:00:00 /mnt/data/liuchao/abag-rbfep/.venv/bin/abag-rbfe resume 3hfm-antibody-h-y33a --batch-dir /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_core_v1_validation_ultra_rescues/abbind-ultra-rescue_3hfm-antibody-h-y33a --execute",
            "liuchao 102 1 0 00:00 ? 00:00:00 bash /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_core_v1_validation_ultra_rescues/abbind-ultra-rescue_3hfm-antibody-h-y33a/jobs/3hfm-antibody-h-y33a/artifacts/commands/sample.sh",
            "liuchao 103 1 0 00:00 ? 00:00:00 gmx mdrun -deffnm /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues/abbind-targeted-lambda-rescue_3hfm-antibody-h-y33a/jobs/3hfm-antibody-h-y33a/legs/complex/repeat_0/window_00/md",
            "liuchao 104 1 0 00:00 ? 00:00:00 bash /mnt/data/liuchao/abag-rbfep/runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues/abbind-targeted-lambda-rescue_3hfm-antibody-h-y33a/jobs/3hfm-antibody-h-y33a/artifacts/commands/sample.sh",
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_stdout),
    )

    counts = module.active_canonical_job_copy_counts()

    assert counts == {"3hfm-antibody-h-y33a": 2}


def test_launch_allowed_by_thread_budget_can_scope_to_selected_plan_roots(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    priority_root = tmp_path / "priority_plan"
    robust_root = tmp_path / "robust_plan"
    deep_root = tmp_path / "deep_rescues"
    ps_stdout = "\n".join(
        [
            f"liuchao 101 1 0 00:00 ? 00:00:00 gmx mdrun -ntmpi 1 -ntomp 2 -deffnm {priority_root}/abbind_1/jobs/job_a/legs/complex/repeat_0/window_00/md",
            f"liuchao 102 1 0 00:00 ? 00:00:00 gmx mdrun -ntmpi 1 -ntomp 2 -deffnm {robust_root}/abbind_2/jobs/job_b/legs/complex/repeat_0/window_00/md",
            f"liuchao 103 1 0 00:00 ? 00:00:00 gmx mdrun -ntmpi 1 -ntomp 4 -deffnm {deep_root}/abbind_3/jobs/job_c/legs/complex/repeat_0/window_00/md",
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_stdout),
    )

    allowed, active_threads, process_count = module.launch_allowed_by_thread_budget(
        max_active_mdrun_threads=5,
        plan_roots=[priority_root, robust_root],
    )

    assert allowed is True
    assert active_threads == 4
    assert process_count == 2


def test_merged_priority_job_rows_prefers_validation_split_selection_over_stale_or_narrower_candidates(
    monkeypatch,
    tmp_path,
) -> None:
    module = _load_watch_module()
    merged_root = tmp_path / "priority_plan"
    merged_dir = merged_root / "reports" / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)
    canonical_plan_jobs = merged_dir / "plan_jobs.csv"
    canonical_plan_jobs.write_text(
        "\n".join(
            [
                "job_id,batch_id,abs_ddg_error_kcal_mol",
                "canonical-job,canonical-batch,1.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    validation_dir = merged_dir / "selections" / "split-validation-fake"
    validation_dir.mkdir(parents=True, exist_ok=True)
    (validation_dir / "plan_summary.json").write_text(
        json.dumps({"selection": {"split_name": "validation"}}),
        encoding="utf-8",
    )
    (validation_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                "job_id,batch_id,abs_ddg_error_kcal_mol",
                "validation-job,validation-batch,5.5",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    complex_dir = merged_dir / "selections" / "complex-1bj1"
    complex_dir.mkdir(parents=True, exist_ok=True)
    (complex_dir / "plan_summary.json").write_text(
        json.dumps({"selection": {"split_name": ""}}),
        encoding="utf-8",
    )
    (complex_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                "job_id,batch_id,abs_ddg_error_kcal_mol",
                "complex-job,complex-batch,9.9",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    validation_time = canonical_plan_jobs.stat().st_mtime + 5
    complex_time = validation_time + 5
    os.utime(validation_dir / "plan_summary.json", (validation_time, validation_time))
    os.utime(validation_dir / "plan_jobs.csv", (validation_time, validation_time))
    os.utime(complex_dir / "plan_summary.json", (complex_time, complex_time))
    os.utime(complex_dir / "plan_jobs.csv", (complex_time, complex_time))

    monkeypatch.setattr(module, "MERGED_PLAN_ROOT", merged_root)
    monkeypatch.setattr(module, "SPLIT_NAME", "validation")

    rows = module.merged_priority_job_rows()

    assert list(rows) == ["validation-job"]
    assert rows["validation-job"]["batch_id"] == "validation-batch"


def test_active_mdrun_threads_sums_ntmpi_times_ntomp(monkeypatch) -> None:
    module = _load_watch_module()
    plan_root = str(module.RUNS_ROOT)
    ps_output = "\n".join(
        [
            (
                "liuchao  1111  1 99 00:00 ? 00:00:00 "
                f"/path/to/gmx mdrun -s {plan_root}/abbind_1cz8_core_v1/jobs/1cz8-antigen-w-g92a/legs/complex/rep01/lambda_001/topol.tpr "
                "-ntmpi 2 -ntomp 3"
            ),
            (
                "liuchao  2222  1 99 00:00 ? 00:00:00 "
                f"/path/to/gmx mdrun -s {plan_root}/abbind_3nps_core_v1/jobs/3nps-antigen-a-h138a/legs/complex/rep01/lambda_001/topol.tpr "
                "-ntomp 4"
            ),
            (
                "liuchao  3333  1 99 00:00 ? 00:00:00 "
                f"/path/to/gmx mdrun -s {plan_root}/abbind_3nps_core_v1/jobs/3nps-antigen-a-y141a/legs/complex/rep01/lambda_000/topol.tpr"
            ),
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_output),
    )

    active_threads, process_count = module.active_mdrun_threads()

    assert active_threads == 11
    assert process_count == 3


def test_stale_mdrun_statuses_uses_recent_progress_signal(monkeypatch) -> None:
    module = _load_watch_module()
    monkeypatch.setattr(
        module,
        "active_mdrun_statuses",
        lambda: [
            {
                "job_id": "1bj1-antigen-v-f17a",
                "pid": 1111,
                "elapsed_seconds": 1800,
                "elapsed_minutes": 30.0,
                "deffnm_tail": "complex/rep01/equilibration/npt",
                "progress_age_seconds": 950.0,
            },
            {
                "job_id": "1cz8-antigen-w-g92a",
                "pid": 2222,
                "elapsed_seconds": 1200,
                "elapsed_minutes": 20.0,
                "deffnm_tail": "apo/rep03/lambda_004/md",
                "progress_age_seconds": 12.0,
            },
        ],
    )

    stale = module.stale_mdrun_statuses(warn_stale_mdrun_seconds=900)

    assert stale == [
        {
            "job_id": "1bj1-antigen-v-f17a",
            "pid": 1111,
            "elapsed_seconds": 1800,
            "elapsed_minutes": 30.0,
            "deffnm_tail": "complex/rep01/equilibration/npt",
            "progress_age_seconds": 950.0,
        }
    ]


def test_refresh_reports_uses_configured_plan_root_and_split(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path / "abbind_core_v1_validation_robust_plan"
    module.PLAN_ROOTS = None
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    captured: list[list[str]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append(command) or SimpleNamespace(returncode=0, stdout=""),
    )

    module.refresh_reports(dry_run=False)

    assert captured == [
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(module.RUNS_ROOT),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ]
    ]


def test_refresh_reports_uses_all_configured_plan_roots(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    primary = tmp_path / "rescues"
    extra = tmp_path / "main"
    module.PLAN_ROOTS = [primary, extra]
    module.RUNS_ROOT = primary
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.REPORTS_DIR = primary / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    captured: list[list[str]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append(command) or SimpleNamespace(returncode=0, stdout=""),
    )

    module.refresh_reports(dry_run=False)

    assert captured == [
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(primary),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ],
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(extra),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ],
    ]


def test_refresh_reports_runs_post_refresh_command_after_report_updates(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path / "validation"
    module.PLAN_ROOTS = None
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.POST_REPORT_REFRESH_COMMAND = ["/tmp/refresh-calibrated-validation.sh"]
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    captured: list[tuple[list[str], dict[str, str] | None]] = []
    lock_states: list[bool] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append((command, kwargs.get("env")))
        or lock_states.append(module._report_refresh_lock_path().exists())
        or SimpleNamespace(returncode=0, stdout="", stderr=""),
    )

    module.refresh_reports(dry_run=False)

    assert captured == [
        (
            [
                str(module.ABAG_RBFE),
                "batch",
                "report-abbind",
                "--plan-root",
                str(module.RUNS_ROOT),
                "--split-name",
                "validation",
                "--split-file",
                str(module.SPLIT_FILE),
            ],
            None,
        ),
        (
            ["/tmp/refresh-calibrated-validation.sh"],
            {"PLAN_ROOT": str(module.RUNS_ROOT)},
        ),
    ]
    assert lock_states == [True, False]
    assert not module._report_refresh_lock_path().exists()


def test_refresh_reports_logs_progress_for_report_and_post_refresh(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path / "validation"
    module.PLAN_ROOTS = None
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.POST_REPORT_REFRESH_COMMAND = ["/tmp/refresh-calibrated-validation.sh"]
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout="", stderr=""),
    )

    module.refresh_reports(dry_run=False)

    out = capsys.readouterr().out
    assert f"[watch] report refresh plan_root={module.RUNS_ROOT} start:" in out
    assert f"[watch] report refresh plan_root={module.RUNS_ROOT} done rc=0" in out
    assert "[watch] post-refresh command env=PLAN_ROOT=" in out
    assert "[watch] post-refresh command env=PLAN_ROOT=" in out
    assert "done rc=0" in out


def test_refresh_reports_dry_run_logs_progress_without_calling_run_command(
    monkeypatch, tmp_path, capsys
) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path / "validation"
    module.PLAN_ROOTS = None
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.POST_REPORT_REFRESH_COMMAND = ["/tmp/refresh-calibrated-validation.sh"]
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    monkeypatch.setattr(
        module,
        "run_command",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("run_command should not be called in dry-run")),
    )

    module.refresh_reports(dry_run=True)

    out = capsys.readouterr().out
    assert f"[watch] dry-run report refresh plan_root={module.RUNS_ROOT}:" in out
    assert "[watch] dry-run post-refresh command env=PLAN_ROOT=" in out


def test_refresh_reports_passes_merged_plan_root_to_post_refresh_command(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    primary = tmp_path / "priority"
    robust = tmp_path / "robust"
    rescue = tmp_path / "rescue"
    module.RUNS_ROOT = robust
    module.PLAN_ROOTS = [robust]
    module.MERGED_PLAN_ROOT = primary
    module.MERGED_EXTRA_PLAN_ROOTS = [robust, rescue]
    module.POST_REPORT_REFRESH_COMMAND = ["/tmp/refresh-calibrated-validation.sh"]
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    captured: list[tuple[list[str], dict[str, str] | None]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append((command, kwargs.get("env")))
        or SimpleNamespace(returncode=0, stdout="", stderr=""),
    )

    module.refresh_reports(dry_run=False)

    assert captured[-1] == (
        ["/tmp/refresh-calibrated-validation.sh"],
        {
            "PLAN_ROOT": str(primary),
            "MERGED_PLAN_ROOT": str(primary),
            "MERGED_EXTRA_PLAN_ROOTS": f"{robust}{module.os.pathsep}{rescue}",
        },
    )


def test_refresh_reports_can_refresh_merged_report(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    primary = tmp_path / "priority"
    robust = tmp_path / "robust"
    rescue = tmp_path / "rescue"
    module.PLAN_ROOTS = [robust]
    module.RUNS_ROOT = robust
    module.MERGED_PLAN_ROOT = primary
    module.MERGED_EXTRA_PLAN_ROOTS = [robust, rescue]
    module.REPORTS_DIR = robust / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    captured: list[list[str]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append(command) or SimpleNamespace(returncode=0, stdout=""),
    )

    module.refresh_reports(dry_run=False)

    assert captured == [
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(robust),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ],
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(rescue),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ],
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(primary),
            "--extra-plan-root",
            str(robust),
            "--extra-plan-root",
            str(rescue),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ],
    ]


def test_refresh_reports_skips_when_global_refresh_lock_is_held(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path / "validation"
    module.PLAN_ROOTS = None
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.POST_REPORT_REFRESH_COMMAND = []
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    lock_path = module._report_refresh_lock_path()
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path.write_text(json.dumps({"pid": 4242, "started_at": "2026-06-09T00:00:00Z"}) + "\n", encoding="utf-8")
    monkeypatch.setattr(module, "_pid_is_alive", lambda pid: pid == 4242)

    captured: list[list[str]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append(command) or SimpleNamespace(returncode=0, stdout=""),
    )

    module.refresh_reports(dry_run=False)

    out = capsys.readouterr().out
    assert "skipping report refresh" in out
    assert captured == []
    assert lock_path.exists()


def test_refresh_reports_reclaims_stale_global_refresh_lock(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path / "validation"
    module.PLAN_ROOTS = None
    module.MERGED_PLAN_ROOT = None
    module.MERGED_EXTRA_PLAN_ROOTS = []
    module.POST_REPORT_REFRESH_COMMAND = []
    module.REPORTS_DIR = module.RUNS_ROOT / "reports"
    module.LOG_DIR = module.REPORTS_DIR / "watch"
    module.SPLIT_NAME = "validation"
    module.SPLIT_FILE = tmp_path / "custom_split.yml"

    lock_path = module._report_refresh_lock_path()
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path.write_text(json.dumps({"pid": 4242, "started_at": "2026-06-09T00:00:00Z"}) + "\n", encoding="utf-8")
    monkeypatch.setattr(module, "_pid_is_alive", lambda pid: False)

    captured: list[list[str]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append(command) or SimpleNamespace(returncode=0, stdout=""),
    )

    module.refresh_reports(dry_run=False)

    assert captured == [
        [
            str(module.ABAG_RBFE),
            "batch",
            "report-abbind",
            "--plan-root",
            str(module.RUNS_ROOT),
            "--split-name",
            "validation",
            "--split-file",
            str(module.SPLIT_FILE),
        ]
    ]
    assert not lock_path.exists()


def test_job_queue_only_listed_skips_full_plan_fallback(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.RUNS_ROOT = tmp_path
    listed_job = tmp_path / "abbind_1cz8_core_v1" / "jobs" / "1cz8-antigen-w-g92a"
    other_job = tmp_path / "abbind_3hfm_core_v1" / "jobs" / "3hfm-antibody-h-c95a"
    listed_job.mkdir(parents=True)
    other_job.mkdir(parents=True)

    queue = module.job_queue([listed_job.name], only_listed=True)

    assert queue == [listed_job]


def test_main_defaults_explicit_job_ids_to_listed_only(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    args = SimpleNamespace(
        plan_root=str(tmp_path / "plan"),
        extra_plan_root=[],
        merged_plan_root="",
        merged_extra_plan_root=[],
        split_name="validation",
        split_file=str(tmp_path / "split.yml"),
        poll_seconds=60,
        gpu_devices="",
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        launch_cooldown_seconds=180,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=900,
        mdrun_args_override="-ntmpi 1 -ntomp 2",
        max_launches_per_pass=2,
        once=True,
        dry_run=True,
        only_listed=False,
        append_rest=False,
        wait_for_pid=0,
        job_id=["3hfm-antibody-h-y50a", "1cz8-antigen-w-g92a"],
    )
    observed: dict[str, object] = {}
    monkeypatch.setattr(module, "parse_args", lambda: args)
    monkeypatch.setattr(module, "normalize_gpu_devices", lambda raw: [])
    monkeypatch.setattr(
        module,
        "pass_once",
        lambda job_ids, gpu_devices, **kwargs: observed.update(
            {
                "job_ids": job_ids,
                "gpu_devices": gpu_devices,
                "only_listed": kwargs["only_listed"],
            }
        ),
    )

    exit_code = module.main()

    assert exit_code == 0
    assert observed == {
        "job_ids": ["3hfm-antibody-h-y50a", "1cz8-antigen-w-g92a"],
        "gpu_devices": [],
        "only_listed": True,
    }
    out = capsys.readouterr().out
    assert "mode=listed-only" in out


def test_main_allows_append_rest_opt_in_for_explicit_job_ids(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    args = SimpleNamespace(
        plan_root=str(tmp_path / "plan"),
        extra_plan_root=[],
        merged_plan_root="",
        merged_extra_plan_root=[],
        split_name="validation",
        split_file=str(tmp_path / "split.yml"),
        poll_seconds=60,
        gpu_devices="",
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        launch_cooldown_seconds=180,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=900,
        mdrun_args_override="-ntmpi 1 -ntomp 2",
        max_launches_per_pass=2,
        once=True,
        dry_run=True,
        only_listed=False,
        append_rest=True,
        wait_for_pid=0,
        job_id=["3hfm-antibody-h-y50a", "1cz8-antigen-w-g92a"],
    )
    observed: dict[str, object] = {}
    monkeypatch.setattr(module, "parse_args", lambda: args)
    monkeypatch.setattr(module, "normalize_gpu_devices", lambda raw: [])
    monkeypatch.setattr(
        module,
        "pass_once",
        lambda job_ids, gpu_devices, **kwargs: observed.update(
            {
                "job_ids": job_ids,
                "gpu_devices": gpu_devices,
                "only_listed": kwargs["only_listed"],
            }
        ),
    )

    exit_code = module.main()

    assert exit_code == 0
    assert observed == {
        "job_ids": ["3hfm-antibody-h-y50a", "1cz8-antigen-w-g92a"],
        "gpu_devices": [],
        "only_listed": False,
    }
    out = capsys.readouterr().out
    assert "mode=listed-plus-rest" in out


def test_main_refuses_duplicate_watcher_lock(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    plan_root = tmp_path / "plan"
    watch_dir = plan_root / "reports" / "watch"
    watch_dir.mkdir(parents=True)
    (watch_dir / module.WATCH_LOCK_NAME).write_text(
        json.dumps(
            {
                "pid": 4242,
                "started_at": "2026-06-08T10:00:00Z",
                "queue_mode": "listed-only",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    args = SimpleNamespace(
        plan_root=str(plan_root),
        extra_plan_root=[],
        merged_plan_root="",
        merged_extra_plan_root=[],
        split_name="validation",
        split_file=str(tmp_path / "split.yml"),
        poll_seconds=60,
        gpu_devices="",
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        launch_cooldown_seconds=180,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=900,
        mdrun_args_override="-ntmpi 1 -ntomp 2",
        max_launches_per_pass=2,
        once=True,
        dry_run=True,
        only_listed=True,
        append_rest=False,
        wait_for_pid=0,
        job_id=["3hfm-antibody-h-y50a"],
    )
    monkeypatch.setattr(module, "parse_args", lambda: args)
    monkeypatch.setattr(module, "normalize_gpu_devices", lambda raw: [])
    monkeypatch.setattr(module, "_pid_is_alive", lambda pid: pid == 4242)
    monkeypatch.setattr(module, "pass_once", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("pass_once should not run")))

    exit_code = module.main()

    assert exit_code == 75
    err = capsys.readouterr().err
    assert "Watcher already running" in err
    assert "PID 4242" in err


def test_main_reclaims_stale_watcher_lock(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    plan_root = tmp_path / "plan"
    watch_dir = plan_root / "reports" / "watch"
    watch_dir.mkdir(parents=True)
    lock_path = watch_dir / module.WATCH_LOCK_NAME
    lock_path.write_text(json.dumps({"pid": 4242}) + "\n", encoding="utf-8")
    args = SimpleNamespace(
        plan_root=str(plan_root),
        extra_plan_root=[],
        merged_plan_root="",
        merged_extra_plan_root=[],
        split_name="validation",
        split_file=str(tmp_path / "split.yml"),
        poll_seconds=60,
        gpu_devices="",
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        launch_cooldown_seconds=180,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=900,
        mdrun_args_override="-ntmpi 1 -ntomp 2",
        max_launches_per_pass=2,
        once=True,
        dry_run=True,
        only_listed=False,
        append_rest=False,
        wait_for_pid=0,
        job_id=["3hfm-antibody-h-y50a"],
    )
    observed: dict[str, object] = {}
    monkeypatch.setattr(module, "parse_args", lambda: args)
    monkeypatch.setattr(module, "normalize_gpu_devices", lambda raw: [])
    monkeypatch.setattr(module, "_pid_is_alive", lambda pid: False)
    monkeypatch.setattr(
        module,
        "pass_once",
        lambda job_ids, gpu_devices, **kwargs: observed.update(
            {
                "job_ids": job_ids,
                "gpu_devices": gpu_devices,
                "only_listed": kwargs["only_listed"],
            }
        ),
    )

    exit_code = module.main()

    assert exit_code == 0
    assert observed == {
        "job_ids": ["3hfm-antibody-h-y50a"],
        "gpu_devices": [],
        "only_listed": True,
    }
    assert not lock_path.exists()


def test_main_allows_distinct_watch_tag_lock(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    plan_root = tmp_path / "plan"
    watch_dir = plan_root / "reports" / "watch"
    watch_dir.mkdir(parents=True)
    (watch_dir / module.WATCH_LOCK_NAME).write_text(
        json.dumps(
            {
                "pid": 4242,
                "started_at": "2026-06-09T06:10:45Z",
                "queue_mode": "listed-only",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    args = SimpleNamespace(
        plan_root=str(plan_root),
        extra_plan_root=[],
        merged_plan_root="",
        merged_extra_plan_root=[],
        split_name="validation",
        split_file=str(tmp_path / "split.yml"),
        poll_seconds=60,
        gpu_devices="",
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        launch_cooldown_seconds=180,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=900,
        mdrun_args_override="-ntmpi 1 -ntomp 2",
        max_launches_per_pass=2,
        once=True,
        dry_run=True,
        only_listed=True,
        append_rest=False,
        wait_for_pid=0,
        watch_tag="gap",
        job_id=["2nz9-antigen-a-f953a"],
    )
    observed: dict[str, object] = {}
    monkeypatch.setattr(module, "parse_args", lambda: args)
    monkeypatch.setattr(module, "normalize_gpu_devices", lambda raw: [])
    monkeypatch.setattr(module, "_pid_is_alive", lambda pid: pid == 4242)
    monkeypatch.setattr(
        module,
        "pass_once",
        lambda job_ids, gpu_devices, **kwargs: observed.update(
            {
                "job_ids": job_ids,
                "gpu_devices": gpu_devices,
                "only_listed": kwargs["only_listed"],
            }
        ),
    )

    exit_code = module.main()

    assert exit_code == 0
    assert observed == {
        "job_ids": ["2nz9-antigen-a-f953a"],
        "gpu_devices": [],
        "only_listed": True,
    }
    assert (watch_dir / "watch.gap.lock.json").exists() is False


def test_active_job_ids_detects_equilibrate_bar_and_cli_processes(monkeypatch) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    plan_root = str(module.RUNS_ROOT)
    ps_output = "\n".join(
        [
            (
                "liuchao  1234     1  0 00:00 ? 00:00:00 "
                f"{plan_root}/.venv/bin/python3.11 {module.ABAG_RBFE} resume 1bj1-antigen-w-i83a "
                f"--batch-dir {plan_root}/abbind_1bj1_core_v1 --execute"
            ),
            (
                "liuchao  2345  1234  0 00:00 ? 00:00:00 bash "
                f"{plan_root}/abbind_1bj1_core_v1/jobs/1bj1-antigen-v-f17a/artifacts/commands/equilibrate.sh"
            ),
            (
                "liuchao  3456  1234  0 00:00 ? 00:00:00 bash "
                f"{plan_root}/abbind_3nps_core_v1/jobs/3nps-antigen-a-h138a/artifacts/commands/bar.sh"
            ),
            (
                "liuchao  4567  1234 99 00:00 ? 00:00:00 gmx mdrun -s "
                f"{plan_root}/abbind_1cz8_core_v1/jobs/1cz8-antigen-w-g92a/legs/complex/rep01/lambda_001/topol.tpr"
            ),
            (
                "liuchao  5678  1234  0 00:00 ? 00:00:00 python3 -c "
                f"\"print('{plan_root}/abbind_1bj1_core_v1/jobs/1bj1-antigen-w-g92a/stages/sample.json')\""
            ),
            (
                "liuchao  6789  1234  0 00:00 ? 00:00:00 "
                f"{module.ABAG_RBFE} batch report-abbind --plan-root {plan_root}"
            ),
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_output),
    )

    active = module.active_job_ids()

    assert active == {
        "1bj1-antigen-w-i83a",
        "1bj1-antigen-v-f17a",
        "3nps-antigen-a-h138a",
        "1cz8-antigen-w-g92a",
    }


def test_active_job_ids_ignores_shell_template_job_literals(monkeypatch) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    plan_root = str(module.RUNS_ROOT)
    ps_output = "\n".join(
        [
            (
                "liuchao  1234  1  0 00:00 ? 00:00:00 /bin/bash -c "
                f"for job in 1mlc-antibody-h-s57a; do {module.ABAG_RBFE} run "
                f"{plan_root}/abbind_1mlc_core_v1/jobs/$job --from-stage equilibrate --execute & done; wait"
            ),
            (
                "liuchao  2345  1  0 00:00 ? 00:00:00 bash "
                f"{plan_root}/abbind_1mlc_core_v1/jobs/1mlc-antibody-h-s57a/artifacts/commands/equilibrate.sh"
            ),
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_output),
    )

    active = module.active_job_ids()

    assert active == {"1mlc-antibody-h-s57a"}


def test_active_canonical_job_ids_detects_jobs_outside_current_plan_root(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.ROOT = tmp_path
    module.PLAN_ROOTS = [tmp_path / "runs" / "benchmarks" / "robust"]
    other_root = tmp_path / "runs" / "benchmarks" / "priority"
    current_root = tmp_path / "runs" / "benchmarks" / "robust"
    ps_output = "\n".join(
        [
            (
                "liuchao  1234     1 99 00:00 ? 00:00:00 gmx mdrun -s "
                f"{other_root}/abbind_3hfm_core_v1/jobs/3hfm-antibody-l-n31a/legs/complex/rep01/lambda_000/topol.tpr"
            ),
            (
                "liuchao  2345     1  0 00:00 ? 00:00:00 "
                f"{current_root}/.venv/bin/python3.11 {module.ABAG_RBFE} resume 2nz9-antigen-a-h1064a "
                f"--batch-dir {current_root}/abbind_2nz9_core_v1 --execute"
            ),
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_output),
    )

    active = module.active_canonical_job_ids()

    assert active == {
        "3hfm-antibody-l-n31a",
        "2nz9-antigen-a-h1064a",
    }


def test_default_queue_excludes_known_blocked_live_job() -> None:
    module = _load_watch_module()

    assert "2nz9-antigen-a-h1064a" not in module.DEFAULT_QUEUE


def test_invalid_mutate_output_canonical_job_ids_detects_canonical_job_ids() -> None:
    module = _load_watch_module()

    blocked = module.invalid_mutate_output_canonical_job_ids(
        {
            "priority/abbind_2nz9_core_v1/2nz9-antigen-a-h1064a": {
                "job_id": "2nz9-antigen-a-h1064a",
                "current_invalid_mutate_output": "True",
                "current_invalid_mutate_output_code": "mutate_processed_gro_isolated_residue_hydrogen",
            },
            "robust/abbind_2nz9_core_v1/2nz9-antigen-a-f953a": {
                "job_id": "2nz9-antigen-a-f953a",
            },
            "rescues/abbind_misc/3hfm-antibody-h-y33a": {
                "current_invalid_mutate_output_code": "mutate_processed_gro_invalid",
            },
        }
    )

    assert blocked == {
        "2nz9-antigen-a-h1064a",
        "3hfm-antibody-h-y33a",
    }


def test_active_canonical_job_ids_ignores_shell_template_job_literals(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.ROOT = tmp_path
    module.PLAN_ROOTS = [tmp_path / "runs" / "benchmarks" / "robust"]
    current_root = tmp_path / "runs" / "benchmarks" / "robust"
    ps_output = "\n".join(
        [
            (
                "liuchao  1234  1  0 00:00 ? 00:00:00 /bin/bash -c "
                f"for job in 1mlc-antibody-h-s57a; do {module.ABAG_RBFE} run "
                f"{current_root}/abbind_1mlc_core_v1/jobs/$job --from-stage equilibrate --execute & done; wait"
            ),
            (
                "liuchao  2345  1 99 00:00 ? 00:00:00 gmx mdrun -s "
                f"{current_root}/abbind_1mlc_core_v1/jobs/1mlc-antibody-h-s57a/legs/complex/rep01/lambda_000/topol.tpr"
            ),
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_output),
    )

    active = module.active_canonical_job_ids()

    assert active == {"1mlc-antibody-h-s57a"}


def test_resume_launch_environment_assigns_pin_offset_from_gpu_slot(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    job_dir = tmp_path / "abbind_1bj1_core_v1" / "jobs" / "1bj1-antigen-w-g88a"
    job_dir.mkdir(parents=True)
    (job_dir / "job_spec.json").write_text(
        json.dumps({"protocol": {"mdrun_args": "-ntmpi 1 -ntomp 4"}}),
        encoding="utf-8",
    )
    monkeypatch.setattr(module.os, "cpu_count", lambda: 64)

    env = module.resume_launch_environment(
        job_dir,
        "2",
        gpu_devices=["0", "1", "2", "3"],
        gpu_counts={"0": 1, "1": 0, "2": 2, "3": 0},
        max_compute_apps_per_gpu=3,
    )

    assert env == {
        "CUDA_VISIBLE_DEVICES": "2",
        "ABAG_RBFE_MDRUN_PINOFFSET": "32",
        "ABAG_RBFE_MDRUN_PINSTRIDE": "1",
    }


def test_resume_launch_environment_uses_mdrun_args_override_for_thread_count(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    job_dir = tmp_path / "abbind_1bj1_core_v1" / "jobs" / "1bj1-antigen-w-g88a"
    job_dir.mkdir(parents=True)
    (job_dir / "job_spec.json").write_text(
        json.dumps({"protocol": {"mdrun_args": "-ntmpi 1 -ntomp 4"}}),
        encoding="utf-8",
    )
    monkeypatch.setattr(module.os, "cpu_count", lambda: 64)

    env = module.resume_launch_environment(
        job_dir,
        "2",
        gpu_devices=["0", "1", "2", "3"],
        gpu_counts={"0": 1, "1": 0, "2": 2, "3": 0},
        max_compute_apps_per_gpu=3,
        mdrun_args_override="-ntmpi 1 -ntomp 2",
    )

    assert env == {
        "CUDA_VISIBLE_DEVICES": "2",
        "ABAG_RBFE_MDRUN_ARGS": "-ntmpi 1 -ntomp 2",
        "ABAG_RBFE_MDRUN_PINOFFSET": "16",
        "ABAG_RBFE_MDRUN_PINSTRIDE": "1",
    }


def test_desired_active_mdrun_affinities_spreads_slots_by_gpu_position() -> None:
    module = _load_watch_module()

    affinity = module.desired_active_mdrun_affinities(
        [
            {"pid": 101, "job_id": "job-a", "gpu_device": "0", "thread_count": 4},
            {"pid": 102, "job_id": "job-b", "gpu_device": "0", "thread_count": 4},
            {"pid": 201, "job_id": "job-c", "gpu_device": "1", "thread_count": 4},
        ],
        gpu_devices=["0", "1"],
        max_compute_apps_per_gpu=3,
        cpu_count=32,
    )

    assert affinity == {
        101: (0, 1, 2, 3),
        102: (4, 5, 6, 7),
        201: (12, 13, 14, 15),
    }


def test_rebalance_active_mdrun_affinity_updates_mismatched_processes(monkeypatch) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    monkeypatch.setattr(module.os, "cpu_count", lambda: 32)
    monkeypatch.setattr(
        module,
        "active_mdrun_processes",
        lambda: [
            {"pid": 101, "job_id": "job-a", "command": "gmx mdrun ...", "thread_count": 4},
            {"pid": 201, "job_id": "job-b", "command": "gmx mdrun ...", "thread_count": 4},
        ],
    )
    monkeypatch.setattr(module, "gpu_pid_to_index", lambda: {101: "0", 201: "1"})
    monkeypatch.setattr(module, "process_affinity", lambda pid: (0, 1, 2, 3) if pid == 101 else (0, 1, 2, 3))
    applied: list[tuple[int, tuple[int, ...]]] = []
    monkeypatch.setattr(module, "set_process_affinity", lambda pid, cpus: applied.append((pid, cpus)) or True)

    changes = module.rebalance_active_mdrun_affinity(
        gpu_devices=["0", "1"],
        max_compute_apps_per_gpu=3,
        dry_run=False,
    )

    assert applied == [(201, (12, 13, 14, 15))]
    assert changes == [
        {
            "pid": 201,
            "job_id": "job-b",
            "gpu_device": "1",
            "thread_count": 4,
            "current_cpus": [0, 1, 2, 3],
            "target_cpus": [12, 13, 14, 15],
        }
    ]


def test_job_queue_keeps_duplicate_job_ids_distinct_across_plan_roots(tmp_path) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    main_root = tmp_path / "main"
    rescue_job = rescue_root / "abbind-rescue_dup" / "jobs" / "dup-job"
    main_job = main_root / "abbind_main_dup" / "jobs" / "dup-job"
    other_job = main_root / "abbind_main_other" / "jobs" / "other-job"
    rescue_job.mkdir(parents=True)
    main_job.mkdir(parents=True)
    other_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root, main_root]

    queue = module.job_queue(["dup-job"])

    assert queue == [rescue_job, main_job, other_job]
    assert module.job_cache_key(rescue_job) != module.job_cache_key(main_job)


def test_pass_once_prioritizes_more_advanced_and_unpaired_jobs(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch_a = tmp_path / "abbind_3hfm_core_v1" / "jobs"
    batch_b = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    batch_c = tmp_path / "abbind_1cz8_core_v1" / "jobs"
    job_build = batch_a / "3hfm-antibody-l-y50a"
    job_equilibrated = batch_b / "1bj1-antigen-w-g88a"
    job_build_unpaired = batch_c / "1cz8-antigen-w-h86a"
    jobs = [job_build, job_equilibrated, job_build_unpaired]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0, "1": 0})
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed" if job_dir == job_equilibrated else "",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "3hfm-antibody-l-y50a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "0",
                    "equilibrate_total_repeats": "6",
                },
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
                "1cz8-antigen-w-h86a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "0",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_3hfm_core_v1": {"paired_job_count": "3"},
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
                "abbind_1cz8_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0", "1"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [
        ("1bj1-antigen-w-g88a", "0"),
        ("1cz8-antigen-w-h86a", "1"),
    ]
    assert refresh_called == [True]


def test_pass_once_blocks_jobs_with_invalid_mutate_output(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_2nz9_core_v1" / "jobs"
    job = batch / "2nz9-antigen-a-h1064a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "priority/abbind_2nz9_core_v1/2nz9-antigen-a-h1064a": {
                    "job_id": "2nz9-antigen-a-h1064a",
                    "current_invalid_mutate_output": "True",
                    "current_invalid_mutate_output_code": "mutate_processed_gro_isolated_residue_hydrogen",
                },
            },
            {},
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []
    assert refresh_called == []
    out = capsys.readouterr().out
    assert '"blocked": ["2nz9-antigen-a-h1064a"]' in out


def test_pass_once_defers_resume_launches_when_cpu_gate_blocks(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    job = batch / "1bj1-antigen-w-g88a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (False, 1.05))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    recent_launches: dict[str, float] = {}
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches=recent_launches,
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []
    assert refresh_called == []
    assert recent_launches == {}


def test_pass_once_skips_launch_when_launch_coordination_lock_conflicts(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    job = batch / "1bj1-antigen-w-g88a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    monkeypatch.setattr(module, "merged_priority_job_rows", lambda: {})
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    released: list[Path | None] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))
    monkeypatch.setattr(module, "acquire_launch_coordination_lock", lambda: (None, "busy"))
    monkeypatch.setattr(module, "release_launch_coordination_lock", lambda lock_path: released.append(lock_path))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []
    assert refresh_called == []
    assert released == []


def test_pass_once_rechecks_thread_budget_under_launch_coordination_lock(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    job = batch / "1bj1-antigen-w-g88a"
    jobs = [job]
    thread_budget_calls = {"count": 0}
    active_calls = {"count": 0}
    release_calls: list[Path | None] = []
    lock_path = tmp_path / "launch.lock.json"

    def launch_allowed_by_thread_budget(**kwargs):
        thread_budget_calls["count"] += 1
        if thread_budget_calls["count"] == 1:
            return True, 0, 0
        return False, 60, 15

    def active_job_ids():
        active_calls["count"] += 1
        return set()

    monkeypatch.setattr(module, "active_job_ids", active_job_ids)
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", launch_allowed_by_thread_budget)
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    monkeypatch.setattr(module, "merged_priority_job_rows", lambda: {})
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))
    monkeypatch.setattr(module, "acquire_launch_coordination_lock", lambda: (lock_path, None))
    monkeypatch.setattr(module, "release_launch_coordination_lock", lambda path: release_calls.append(path))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=56,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []
    assert refresh_called == []
    assert thread_budget_calls["count"] == 2
    assert active_calls["count"] == 2
    assert release_calls == [lock_path]


def test_pass_once_can_launch_using_gpu_headroom_override(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_3hfm_core_v1" / "jobs"
    job = batch / "3hfm-antibody-h-c95a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 14})
    monkeypatch.setattr(
        module,
        "gpu_device_stats",
        lambda: {
            "0": {
                "memory_used_mb": 7274,
                "memory_total_mb": 24564,
                "gpu_utilization_percent": 44,
            }
        },
    )
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(
        module,
        "launch_allowed_by_thread_budget",
        lambda **kwargs: (True, 0, 0),
    )
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "3hfm-antibody-h-c95a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_3hfm_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=4,
        min_free_gpu_memory_mb=12000,
        max_gpu_utilization=60,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [("3hfm-antibody-h-c95a", "0")]
    assert refresh_called == [True]


def test_pass_once_prioritizes_higher_current_validation_error_from_merged_rows(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    low_batch = tmp_path / "abbind_3nps_core_v1" / "jobs"
    high_batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    low_job = low_batch / "3nps-antigen-a-h138a"
    high_job = high_batch / "1bj1-antigen-w-g88a"
    jobs = [low_job, high_job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(
        module,
        "launch_allowed_by_thread_budget",
        lambda **kwargs: (True, 0, 0),
    )
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "3nps-antigen-a-h138a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "4",
                    "equilibrate_total_repeats": "4",
                },
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "4",
                    "equilibrate_total_repeats": "4",
                },
            },
            {
                "abbind_3nps_core_v1": {"paired_job_count": "2"},
                "abbind_1bj1_core_v1": {"paired_job_count": "4"},
            },
        ),
    )
    monkeypatch.setattr(
        module,
        "merged_priority_job_rows",
        lambda: {
            "3nps-antigen-a-h138a": {"abs_ddg_error_kcal_mol": "1.194273489694477"},
            "1bj1-antigen-w-g88a": {"abs_ddg_error_kcal_mol": "5.503060267996389"},
        },
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [("1bj1-antigen-w-g88a", "0")]
    assert refresh_called == [True]


def test_pass_once_refreshes_reports_while_jobs_remain_active(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    job = batch / "1bj1-antigen-w-g88a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: {"1bj1-antigen-w-g88a"})
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(
        module,
        "launch_allowed_by_thread_budget",
        lambda **kwargs: (True, 2, 1),
    )
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "running",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "3",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []
    assert refresh_called == [True]


def test_pass_once_skips_recently_launched_job_until_cooldown_expires(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    job_dir = tmp_path / "abbind_1bj1_core_v1" / "jobs" / "1bj1-antigen-w-g92a"
    jobs = [job_dir]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(
        module,
        "launch_allowed_by_thread_budget",
        lambda **kwargs: (True, 0, 0),
    )
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g92a": {
                    "sample_completed_windows": "30",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                }
            },
            {"abbind_1bj1_core_v1": {"paired_job_count": "0"}},
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    recent_launches = {"1bj1-antigen-w-g92a": 1000.0}
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches=recent_launches,
        dry_run=False,
        now_ts=1100.0,
    )
    assert launches == []
    assert refresh_called == []
    assert recent_launches == {"1bj1-antigen-w-g92a": 1000.0}

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches=recent_launches,
        dry_run=False,
        now_ts=1185.0,
    )
    assert launches == [("1bj1-antigen-w-g92a", "0")]
    assert refresh_called == [True]
    assert recent_launches == {"1bj1-antigen-w-g92a": 1185.0}


def test_pass_once_defers_resume_launches_when_thread_gate_blocks(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    job = batch / "1bj1-antigen-w-g88a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.6))
    monkeypatch.setattr(
        module,
        "launch_allowed_by_thread_budget",
        lambda **kwargs: (False, 56, 14),
    )
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    refresh_called: list[bool] = []
    recent_launches: dict[str, float] = {}
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: refresh_called.append(True))

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=56,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches=recent_launches,
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []
    assert refresh_called == []
    assert recent_launches == {}


def test_pass_once_reports_stale_active_jobs(monkeypatch, tmp_path, capsys) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_1bj1_core_v1" / "jobs"
    job = batch / "1bj1-antigen-w-g88a"
    jobs = [job]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(
        module,
        "launch_allowed_by_thread_budget",
        lambda **kwargs: (False, 56, 14),
    )
    monkeypatch.setattr(
        module,
        "stale_mdrun_statuses",
        lambda **kwargs: [
            {
                "job_id": "1bj1-antigen-v-f17a",
                "pid": 1111,
                "elapsed_minutes": 30.0,
                "progress_age_seconds": 950.0,
                "deffnm_tail": "complex/rep01/equilibration/npt",
            }
        ],
    )
    monkeypatch.setattr(module, "job_queue", lambda job_ids, **kwargs: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                "1bj1-antigen-w-g88a": {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                "abbind_1bj1_core_v1": {"paired_job_count": "0"},
            },
        ),
    )
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "launch_resume", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=56,
        warn_stale_mdrun_seconds=900,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    out = capsys.readouterr().out
    assert '"stale_active": ["1bj1-antigen-v-f17a"]' in out
    assert "stale active mdrun processes" in out


def test_pass_once_uses_root_sensitive_recent_launch_keys(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    main_root = tmp_path / "main"
    rescue_job = rescue_root / "abbind-rescue_dup" / "jobs" / "dup-job"
    main_job = main_root / "abbind_main_dup" / "jobs" / "dup-job"
    rescue_job.mkdir(parents=True)
    main_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root, main_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, only_listed=False: [rescue_job, main_job])
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                module.job_cache_key(rescue_job): {
                    "sample_completed_windows": "24",
                    "sample_total_windows": "64",
                    "equilibrate_completed_repeats": "8",
                    "equilibrate_total_repeats": "8",
                },
                module.job_cache_key(main_job): {
                    "sample_completed_windows": "24",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                module.batch_cache_key(rescue_job.parent.parent): {"paired_job_count": "1"},
                module.batch_cache_key(main_job.parent.parent): {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((module.job_cache_key(job_dir), gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={module.job_cache_key(rescue_job): 1000.0},
        dry_run=False,
        now_ts=1100.0,
    )

    assert launches == [(module.job_cache_key(main_job), "0")]


def test_pass_once_can_launch_not_started_rescue_job(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    rescue_job = rescue_root / "abbind-rescue_dup" / "jobs" / "dup-job"
    rescue_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, only_listed=False: [rescue_job])
    monkeypatch.setattr(module, "read_stage_states", lambda job_dir: {})
    monkeypatch.setattr(module, "report_priority_data", lambda: ({}, {}))
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [("dup-job", "0")]


def test_pass_once_prioritizes_primary_plan_root_jobs(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    main_root = tmp_path / "main"
    rescue_job = rescue_root / "abbind-rescue_dup" / "jobs" / "dup-job"
    main_job = main_root / "abbind_main_other" / "jobs" / "other-job"
    rescue_job.mkdir(parents=True)
    main_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root, main_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, only_listed=False: [rescue_job, main_job])
    monkeypatch.setattr(module, "read_stage_states", lambda job_dir: {})
    monkeypatch.setattr(module, "report_priority_data", lambda: ({}, {}))
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((module.job_cache_key(job_dir), gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [(module.job_cache_key(rescue_job), "0")]


def test_pass_once_limits_resume_launches_per_pass(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    module.PLAN_ROOTS = None
    batch = tmp_path / "abbind_3hfm_core_v1" / "jobs"
    jobs = [
        batch / "3hfm-antibody-h-y50a",
        batch / "3hfm-antibody-h-y33a",
        batch / "3hfm-antibody-h-c95a",
    ]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0, "1": 0, "2": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.2))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, only_listed=False: jobs)
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                job.name: {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "144",
                    "equilibrate_completed_repeats": "0",
                    "equilibrate_total_repeats": "6",
                }
                for job in jobs
            },
            {"abbind_3hfm_core_v1": {"paired_job_count": "0"}},
        ),
    )
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)
    monkeypatch.setattr(
        module,
        "acquire_launch_coordination_lock",
        lambda: (tmp_path / "launch.lock.json", None),
    )
    monkeypatch.setattr(module, "release_launch_coordination_lock", lambda lock_path: None)

    module.pass_once(
        [],
        ["0", "1", "2"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=1.2,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=2,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [
        ("3hfm-antibody-h-y50a", "0"),
        ("3hfm-antibody-h-y33a", "1"),
    ]


def test_pass_once_skips_job_active_in_other_plan_root(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    robust_root = tmp_path / "robust"
    priority_root = tmp_path / "priority"
    robust_job = robust_root / "abbind_3hfm_core_v1" / "jobs" / "3hfm-antibody-l-n31a"
    priority_job = priority_root / "abbind_3hfm_core_v1" / "jobs" / "3hfm-antibody-l-n31a"
    robust_job.mkdir(parents=True)
    priority_job.mkdir(parents=True)
    module.PLAN_ROOTS = [robust_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: {"3hfm-antibody-l-n31a"})
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, only_listed=False: [robust_job])
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                module.job_cache_key(robust_job): {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                module.batch_cache_key(robust_job.parent.parent): {"paired_job_count": "0"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=0,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == []


def test_pass_once_can_launch_job_active_in_other_plan_root_when_allowed(monkeypatch, tmp_path) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    rescue_job = rescue_root / "abbind_3hfm_core_v1" / "jobs" / "3hfm-antibody-l-n31a"
    rescue_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: {"3hfm-antibody-l-n31a"})
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "job_queue", lambda job_ids, only_listed=False: [rescue_job])
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                module.job_cache_key(rescue_job): {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "64",
                    "equilibrate_completed_repeats": "8",
                    "equilibrate_total_repeats": "8",
                },
            },
            {
                module.batch_cache_key(rescue_job.parent.parent): {"paired_job_count": "1"},
            },
        ),
    )
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        allow_active_elsewhere_job_ids=True,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [("3hfm-antibody-l-n31a", "0")]


def test_pass_once_prefers_unique_job_before_active_elsewhere_duplicate_when_allowed(
    monkeypatch, tmp_path
) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    duplicate_job = rescue_root / "abbind_1mlc_core_v1" / "jobs" / "1mlc-antibody-l-n92a"
    unique_job = rescue_root / "abbind_1cz8_core_v1" / "jobs" / "1cz8-antigen-w-q79a"
    duplicate_job.mkdir(parents=True)
    unique_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: {"1mlc-antibody-l-n92a"})
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(
        module,
        "job_queue",
        lambda job_ids, only_listed=False: [duplicate_job, unique_job],
    )
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                module.job_cache_key(duplicate_job): {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
                module.job_cache_key(unique_job): {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                module.batch_cache_key(duplicate_job.parent.parent): {"paired_job_count": "0"},
                module.batch_cache_key(unique_job.parent.parent): {"paired_job_count": "0"},
            },
        ),
    )
    monkeypatch.setattr(module, "merged_priority_job_rows", lambda: {})
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)
    monkeypatch.setattr(
        module,
        "acquire_launch_coordination_lock",
        lambda: (tmp_path / "launch.lock.json", None),
    )
    monkeypatch.setattr(module, "release_launch_coordination_lock", lambda lock_path: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        allow_active_elsewhere_job_ids=True,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [("1cz8-antigen-w-q79a", "0")]


def test_pass_once_respects_active_copy_cap_when_duplicates_are_allowed(
    monkeypatch, tmp_path
) -> None:
    module = _load_watch_module()
    rescue_root = tmp_path / "rescues"
    duplicate_job = rescue_root / "abbind_1mlc_core_v1" / "jobs" / "1mlc-antibody-l-n92a"
    unique_job = rescue_root / "abbind_1cz8_core_v1" / "jobs" / "1cz8-antigen-w-q79a"
    duplicate_job.mkdir(parents=True)
    unique_job.mkdir(parents=True)
    module.PLAN_ROOTS = [rescue_root]

    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "active_canonical_job_ids", lambda: {"1mlc-antibody-l-n92a"})
    monkeypatch.setattr(module, "active_canonical_job_copy_counts", lambda: {"1mlc-antibody-l-n92a": 3})
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(
        module,
        "job_queue",
        lambda job_ids, only_listed=False: [duplicate_job, unique_job],
    )
    monkeypatch.setattr(
        module,
        "read_stage_states",
        lambda job_dir: {
            "report": "",
            "prepare": "completed",
            "build_legs": "completed",
            "equilibrate": "completed",
            "sample": "",
        },
    )
    monkeypatch.setattr(
        module,
        "report_priority_data",
        lambda: (
            {
                module.job_cache_key(duplicate_job): {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
                module.job_cache_key(unique_job): {
                    "sample_completed_windows": "0",
                    "sample_total_windows": "48",
                    "equilibrate_completed_repeats": "6",
                    "equilibrate_total_repeats": "6",
                },
            },
            {
                module.batch_cache_key(duplicate_job.parent.parent): {"paired_job_count": "0"},
                module.batch_cache_key(unique_job.parent.parent): {"paired_job_count": "0"},
            },
        ),
    )
    monkeypatch.setattr(module, "merged_priority_job_rows", lambda: {})
    launches: list[tuple[str, str]] = []
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "rebalance_active_mdrun_affinity", lambda **kwargs: [])
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append((job_dir.name, gpu_device)),
    )
    monkeypatch.setattr(module, "refresh_reports", lambda dry_run: None)
    monkeypatch.setattr(
        module,
        "acquire_launch_coordination_lock",
        lambda: (tmp_path / "launch.lock.json", None),
    )
    monkeypatch.setattr(module, "release_launch_coordination_lock", lambda lock_path: None)

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        allow_active_elsewhere_job_ids=True,
        max_active_copies_per_job_id=3,
        max_compute_apps_per_gpu=1,
        max_load_per_core=0.95,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=1,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=1000.0,
    )

    assert launches == [("1cz8-antigen-w-q79a", "0")]
