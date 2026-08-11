from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace


def _load_watch_module():
    module_path = (
        Path("/mnt/data/liuchao/abag-rbfep/benchmarks/patel_2021_3hfm/watch_patellike_3hfm.py")
    )
    spec = importlib.util.spec_from_file_location("watch_patellike_3hfm", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_job_spec(job_dir: Path, *, charge_conserving: bool, mdrun_args: str = "-ntmpi 1 -ntomp 4") -> None:
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "job_spec.json").write_text(
        json.dumps(
            {
                "job_id": job_dir.name,
                "protocol": {"mdrun_args": mdrun_args},
                "mutation_group": {"charge_conserving": charge_conserving},
            }
        ),
        encoding="utf-8",
    )


def test_available_gpu_devices_prefers_lower_loaded_gpus(monkeypatch) -> None:
    module = _load_watch_module()
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 2, "1": 3, "2": 0, "3": 1})

    available = module.available_gpu_devices(
        ["0", "1", "2", "3"],
        gpu_counts=module.gpu_compute_counts(),
        max_compute_apps_per_gpu=3,
    )

    assert available == ["2", "3", "0"]


def test_available_gpu_devices_uses_headroom_override_when_compute_slots_full() -> None:
    module = _load_watch_module()

    available = module.available_gpu_devices(
        ["0", "1", "2"],
        gpu_counts={"0": 3, "1": 3, "2": 1},
        max_compute_apps_per_gpu=3,
        gpu_stats={
            "0": {
                "memory_used_mb": 1000,
                "memory_total_mb": 24000,
                "gpu_utilization_percent": 20,
                "memory_utilization_percent": 5,
            },
            "1": {
                "memory_used_mb": 18000,
                "memory_total_mb": 24000,
                "gpu_utilization_percent": 90,
                "memory_utilization_percent": 75,
            },
            "2": {
                "memory_used_mb": 12000,
                "memory_total_mb": 24000,
                "gpu_utilization_percent": 30,
                "memory_utilization_percent": 50,
            },
        },
        min_free_gpu_memory_mb=12000,
        max_gpu_utilization=60,
    )

    assert available == ["2", "0"]


def test_job_queue_uses_batch_plan_order_and_respects_only_listed(tmp_path: Path) -> None:
    module = _load_watch_module()
    batch_dir = tmp_path / "patel" / "patel_2021_3hfm_reference"
    jobs_dir = batch_dir / "jobs"
    for job_id in ("job-a", "job-b", "job-c"):
        (jobs_dir / job_id).mkdir(parents=True, exist_ok=True)
    (batch_dir / "batch_plan.json").write_text(
        json.dumps({"jobs": [{"job_id": "job-a"}, {"job_id": "job-b"}, {"job_id": "job-c"}]}),
        encoding="utf-8",
    )
    module.BATCH_DIR = batch_dir

    only_listed = module.job_queue(["job-c", "job-a"], only_listed=True)
    append_rest = module.job_queue(["job-c"], only_listed=False)

    assert [path.name for path in only_listed] == ["job-c", "job-a"]
    assert [path.name for path in append_rest] == ["job-c", "job-a", "job-b"]


def test_resumable_priority_prefers_charge_conserving_jobs(tmp_path: Path) -> None:
    module = _load_watch_module()
    charge_job = tmp_path / "jobs" / "charge"
    neutral_job = tmp_path / "jobs" / "neutral"
    _write_job_spec(charge_job, charge_conserving=False)
    _write_job_spec(neutral_job, charge_conserving=True)
    row = {
        "latest_stage": "build_legs",
        "latest_stage_state": "completed",
        "sample_completed_windows": 0,
        "sample_total_windows": 144,
        "equilibrate_completed_repeats": 0,
        "equilibrate_total_repeats": 6,
    }

    neutral_priority = module.resumable_priority(neutral_job, row, queue_position=0)
    charge_priority = module.resumable_priority(charge_job, row, queue_position=1)

    assert neutral_priority < charge_priority


def test_pass_once_defers_charge_changing_while_charge_conserving_pending(
    tmp_path: Path,
    monkeypatch,
) -> None:
    module = _load_watch_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    jobs_dir = batch_dir / "jobs"
    conserving_job = jobs_dir / "3hfm-patel-2021-antigen-y-y20f"
    charge_job = jobs_dir / "3hfm-patel-2021-antibody-h-d32n"
    _write_job_spec(conserving_job, charge_conserving=True)
    _write_job_spec(charge_job, charge_conserving=False)
    (batch_dir / "batch_plan.json").write_text(
        json.dumps({"jobs": [{"job_id": conserving_job.name}, {"job_id": charge_job.name}]}),
        encoding="utf-8",
    )
    module.BATCH_DIR = batch_dir
    module.LOG_DIR = batch_dir / "reports" / "watch"

    launches: list[str] = []
    refresh_calls: list[bool] = []
    monkeypatch.setattr(
        module,
        "refresh_batch_summary",
        lambda *, dry_run: refresh_calls.append(dry_run)
        or {
            "jobs": [
                {
                    "job_id": conserving_job.name,
                    "ddg_ready": False,
                    "latest_stage": "build_legs",
                    "latest_stage_state": "completed",
                    "analyzable": False,
                    "resumable": True,
                    "sample_completed_windows": 0,
                    "sample_total_windows": 144,
                    "equilibrate_completed_repeats": 0,
                    "equilibrate_total_repeats": 6,
                },
                {
                    "job_id": charge_job.name,
                    "ddg_ready": False,
                    "latest_stage": "build_legs",
                    "latest_stage_state": "completed",
                    "analyzable": False,
                    "resumable": True,
                    "sample_completed_windows": 0,
                    "sample_total_windows": 144,
                    "equilibrate_completed_repeats": 0,
                    "equilibrate_total_repeats": 6,
                },
            ]
        },
    )
    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "gpu_device_stats", lambda: {})
    monkeypatch.setattr(module, "available_gpu_devices", lambda *args, **kwargs: ["0"])
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "run_post_refresh", lambda **kwargs: None)
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append(job_dir.name),
    )

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        skip_charge_changing=True,
        max_compute_apps_per_gpu=1,
        min_free_gpu_memory_mb=0,
        max_gpu_utilization=0,
        max_load_per_core=0.0,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=2,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=0.0,
    )

    assert launches == [conserving_job.name]
    assert len(refresh_calls) == 2


def test_pass_once_releases_charge_changing_after_conserving_subset_finishes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    module = _load_watch_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    jobs_dir = batch_dir / "jobs"
    conserving_job = jobs_dir / "3hfm-patel-2021-antigen-y-y20f"
    charge_job = jobs_dir / "3hfm-patel-2021-antibody-h-d32n"
    _write_job_spec(conserving_job, charge_conserving=True)
    _write_job_spec(charge_job, charge_conserving=False)
    (batch_dir / "batch_plan.json").write_text(
        json.dumps({"jobs": [{"job_id": conserving_job.name}, {"job_id": charge_job.name}]}),
        encoding="utf-8",
    )
    module.BATCH_DIR = batch_dir
    module.LOG_DIR = batch_dir / "reports" / "watch"

    launches: list[str] = []
    refresh_calls: list[bool] = []
    monkeypatch.setattr(
        module,
        "refresh_batch_summary",
        lambda *, dry_run: refresh_calls.append(dry_run)
        or {
            "jobs": [
                {
                    "job_id": conserving_job.name,
                    "ddg_ready": True,
                    "latest_stage": "report",
                    "latest_stage_state": "completed",
                    "analyzable": False,
                    "resumable": False,
                    "sample_completed_windows": 144,
                    "sample_total_windows": 144,
                    "equilibrate_completed_repeats": 6,
                    "equilibrate_total_repeats": 6,
                },
                {
                    "job_id": charge_job.name,
                    "ddg_ready": False,
                    "latest_stage": "build_legs",
                    "latest_stage_state": "completed",
                    "analyzable": False,
                    "resumable": True,
                    "sample_completed_windows": 0,
                    "sample_total_windows": 144,
                    "equilibrate_completed_repeats": 0,
                    "equilibrate_total_repeats": 6,
                },
            ]
        },
    )
    monkeypatch.setattr(module, "active_job_ids", lambda: set())
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "gpu_device_stats", lambda: {})
    monkeypatch.setattr(module, "available_gpu_devices", lambda *args, **kwargs: ["0"])
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "run_post_refresh", lambda **kwargs: None)
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append(job_dir.name),
    )

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        skip_charge_changing=True,
        max_compute_apps_per_gpu=1,
        min_free_gpu_memory_mb=0,
        max_gpu_utilization=0,
        max_load_per_core=0.0,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=2,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=0.0,
    )

    assert launches == [charge_job.name]
    assert len(refresh_calls) == 2


def test_pass_once_allows_charge_changing_once_conserving_subset_is_in_flight(
    tmp_path: Path,
    monkeypatch,
) -> None:
    module = _load_watch_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    jobs_dir = batch_dir / "jobs"
    conserving_job = jobs_dir / "3hfm-patel-2021-antigen-y-y20f"
    charge_job = jobs_dir / "3hfm-patel-2021-antibody-h-d32n"
    _write_job_spec(conserving_job, charge_conserving=True)
    _write_job_spec(charge_job, charge_conserving=False)
    (batch_dir / "batch_plan.json").write_text(
        json.dumps({"jobs": [{"job_id": conserving_job.name}, {"job_id": charge_job.name}]}),
        encoding="utf-8",
    )
    module.BATCH_DIR = batch_dir
    module.LOG_DIR = batch_dir / "reports" / "watch"

    launches: list[str] = []
    refresh_calls: list[bool] = []
    monkeypatch.setattr(
        module,
        "refresh_batch_summary",
        lambda *, dry_run: refresh_calls.append(dry_run)
        or {
            "jobs": [
                {
                    "job_id": conserving_job.name,
                    "ddg_ready": False,
                    "latest_stage": "sample",
                    "latest_stage_state": "running",
                    "analyzable": False,
                    "resumable": True,
                    "sample_completed_windows": 36,
                    "sample_total_windows": 144,
                    "equilibrate_completed_repeats": 6,
                    "equilibrate_total_repeats": 6,
                },
                {
                    "job_id": charge_job.name,
                    "ddg_ready": False,
                    "latest_stage": "build_legs",
                    "latest_stage_state": "completed",
                    "analyzable": False,
                    "resumable": True,
                    "sample_completed_windows": 0,
                    "sample_total_windows": 144,
                    "equilibrate_completed_repeats": 0,
                    "equilibrate_total_repeats": 6,
                },
            ]
        },
    )
    monkeypatch.setattr(module, "active_job_ids", lambda: {conserving_job.name})
    monkeypatch.setattr(module, "gpu_compute_counts", lambda: {"0": 0})
    monkeypatch.setattr(module, "gpu_device_stats", lambda: {})
    monkeypatch.setattr(module, "available_gpu_devices", lambda *args, **kwargs: ["0"])
    monkeypatch.setattr(module, "launch_allowed_by_cpu", lambda **kwargs: (True, 0.0))
    monkeypatch.setattr(module, "launch_allowed_by_thread_budget", lambda **kwargs: (True, 0, 0))
    monkeypatch.setattr(module, "stale_mdrun_statuses", lambda **kwargs: [])
    monkeypatch.setattr(module, "analyze_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(module, "run_post_refresh", lambda **kwargs: None)
    monkeypatch.setattr(
        module,
        "launch_resume",
        lambda job_dir, gpu_device, dry_run, **kwargs: launches.append(job_dir.name),
    )

    module.pass_once(
        [],
        ["0"],
        only_listed=False,
        skip_charge_changing=True,
        max_compute_apps_per_gpu=1,
        min_free_gpu_memory_mb=0,
        max_gpu_utilization=0,
        max_load_per_core=0.0,
        max_active_mdrun_threads=0,
        warn_stale_mdrun_seconds=0,
        mdrun_args_override="",
        max_launches_per_pass=2,
        launch_cooldown_seconds=180,
        recent_launches={},
        dry_run=False,
        now_ts=0.0,
    )

    assert launches == [charge_job.name]
    assert len(refresh_calls) == 2


def test_active_mdrun_threads_sums_repo_jobs(monkeypatch) -> None:
    module = _load_watch_module()
    root_runs = str(module.ROOT / "runs")
    ps_output = "\n".join(
        [
            (
                "liuchao  1111  1 99 00:00 ? 00:00:00 "
                f"/path/to/gmx mdrun -s {root_runs}/benchmarks/patel/jobs/job-a/legs/complex/rep01/lambda_001/topol.tpr "
                "-ntmpi 2 -ntomp 3"
            ),
            (
                "liuchao  2222  1 99 00:00 ? 00:00:00 "
                f"/path/to/gmx mdrun -s {root_runs}/benchmarks/abbind/jobs/job-b/legs/complex/rep01/lambda_001/topol.tpr "
                "-ntomp 4"
            ),
        ]
    )
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: SimpleNamespace(returncode=0, stdout=ps_output),
    )

    active_threads, process_count = module.active_mdrun_threads()

    assert active_threads == 10
    assert process_count == 2


def test_refresh_batch_summary_runs_report_and_post_refresh(monkeypatch, tmp_path: Path) -> None:
    module = _load_watch_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    reports_dir = batch_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / "batch_summary.json").write_text(json.dumps({"jobs": [{"job_id": "job-a"}]}), encoding="utf-8")
    module.BATCH_DIR = batch_dir
    module.LOG_DIR = reports_dir / "watch"
    module.POST_REPORT_REFRESH_COMMAND = ["/tmp/refresh-patel.sh"]

    captured: list[tuple[list[str], dict[str, str] | None]] = []
    monkeypatch.setattr(
        module,
        "run_command",
        lambda command, **kwargs: captured.append((command, kwargs.get("env")))
        or SimpleNamespace(returncode=0, stdout="", stderr=""),
    )

    payload = module.refresh_batch_summary(dry_run=False)
    module.run_post_refresh(dry_run=False)

    assert payload == {"jobs": [{"job_id": "job-a"}]}
    assert captured == [
        (
            [str(module.PYTHON_BIN), "-u", str(module.REPORT_SCRIPT), "--batch-dir", str(batch_dir)],
            None,
        ),
        (
            ["/tmp/refresh-patel.sh"],
            {"BATCH_DIR": str(batch_dir), "PLAN_ROOT": str(batch_dir.parent)},
        ),
    ]
