#!/usr/bin/env python3
"""Monitor the Patel 2021 3HFM batch and keep a small execution lane active."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
ABAG_RBFE = ROOT / ".venv" / "bin" / "abag-rbfe"
PYTHON_BIN = ROOT / ".venv" / "bin" / "python"
BENCHMARK_ROOT = ROOT / "benchmarks" / "patel_2021_3hfm"
REPORT_SCRIPT = BENCHMARK_ROOT / "report_patellike_3hfm.py"
BATCH_DIR = ROOT / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
LOG_DIR = BATCH_DIR / "reports" / "watch"
POST_REPORT_REFRESH_COMMAND: list[str] = []
WATCH_LOCK_NAME = "watch.lock.json"

JOB_RE = re.compile(r"/jobs/([^/]+)/")
CLI_JOB_RE = re.compile(r"\babag-rbfe\b\s+(?:run|resume|analyze)\s+([^\s]+)")
BATCH_DIR_RE = re.compile(r"--batch-dir\s+([^\s]+)")
ACTIVE_PATH_MARKERS = ("/artifacts/commands/", "/legs/")
NTOMP_RE = re.compile(r"(?:^|\s)-ntomp\s+(\d+)(?:\s|$)")
NTMPI_RE = re.compile(r"(?:^|\s)-ntmpi\s+(\d+)(?:\s|$)")
MDRUN_RE = re.compile(r"(?:^|\s)\S*gmx(?:_mpi)?\s+mdrun(?:\s|$)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_id", nargs="*")
    parser.add_argument(
        "--batch-dir",
        default=os.environ.get("BATCH_DIR", str(BATCH_DIR)),
        help="Patel 2021 3HFM batch directory to monitor.",
    )
    parser.add_argument(
        "--post-refresh-command",
        default=os.environ.get("POST_REFRESH_COMMAND", ""),
        help="Optional command executed after the Patel summary refresh.",
    )
    parser.add_argument("--poll-seconds", type=int, default=int(os.environ.get("POLL_SECONDS", "60")))
    parser.add_argument("--gpu-devices", default=os.environ.get("GPU_DEVICES", ""))
    parser.add_argument(
        "--max-compute-apps-per-gpu",
        type=int,
        default=int(os.environ.get("MAX_COMPUTE_APPS_PER_GPU", "13")),
    )
    parser.add_argument(
        "--min-free-gpu-memory-mb",
        type=int,
        default=int(os.environ.get("MIN_FREE_GPU_MEMORY_MB", "0")),
    )
    parser.add_argument(
        "--max-gpu-utilization",
        type=int,
        default=int(os.environ.get("MAX_GPU_UTILIZATION", "0")),
    )
    parser.add_argument(
        "--max-load-per-core",
        type=float,
        default=float(os.environ.get("MAX_LOAD_PER_CORE", "0")),
    )
    parser.add_argument(
        "--max-active-mdrun-threads",
        type=int,
        default=int(os.environ.get("MAX_ACTIVE_MDRUN_THREADS", "0")),
    )
    parser.add_argument(
        "--launch-cooldown-seconds",
        type=int,
        default=int(os.environ.get("LAUNCH_COOLDOWN_SECONDS", "180")),
    )
    parser.add_argument(
        "--warn-stale-mdrun-seconds",
        type=int,
        default=int(os.environ.get("WARN_STALE_MDRUN_SECONDS", "900")),
    )
    parser.add_argument(
        "--mdrun-args-override",
        default=os.environ.get("ABAG_RBFE_MDRUN_ARGS", os.environ.get("MDRUN_ARGS_OVERRIDE", "")),
    )
    parser.add_argument(
        "--max-launches-per-pass",
        type=int,
        default=int(os.environ.get("MAX_LAUNCHES_PER_PASS", "2")),
    )
    parser.add_argument(
        "--skip-charge-changing",
        action="store_true",
        default=os.environ.get("SKIP_CHARGE_CHANGING", "1").strip().lower() in {"1", "true", "yes", "on"},
        help="Skip charge-changing jobs and prioritize the charge-conserving Patel subset.",
    )
    parser.add_argument(
        "--watch-tag",
        default=os.environ.get("WATCH_TAG", "").strip(),
        help="Optional watcher tag used only for lock diagnostics.",
    )
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    queue_scope = parser.add_mutually_exclusive_group()
    queue_scope.add_argument(
        "--only-listed",
        action="store_true",
        help="Restrict the queue to explicitly listed job IDs.",
    )
    queue_scope.add_argument(
        "--append-rest",
        action="store_true",
        help="After explicit job IDs, append the remainder of the Patel batch queue.",
    )
    return parser.parse_args()


def _bool_env(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _read_json(path: Path, default: Any) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return default


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


def _watch_lock_path(watch_tag: str = "") -> Path:
    normalized_tag = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(watch_tag or "").strip()).strip("-.")
    if not normalized_tag:
        return LOG_DIR / WATCH_LOCK_NAME
    return LOG_DIR / f"watch.{normalized_tag}.lock.json"


def acquire_watch_lock(job_ids: list[str], *, only_listed: bool, watch_tag: str = "") -> tuple[Path | None, str | None]:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    lock_path = _watch_lock_path(watch_tag)
    payload = {
        "pid": os.getpid(),
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "batch_dir": str(BATCH_DIR),
        "job_ids": list(job_ids),
        "queue_mode": "listed-only" if only_listed else "listed-plus-rest",
        "watch_tag": watch_tag,
    }
    for _attempt in range(2):
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            existing = _read_json(lock_path, {})
            existing_pid = int(existing.get("pid") or 0)
            if existing_pid and existing_pid != os.getpid() and _pid_is_alive(existing_pid):
                queue_mode = str(existing.get("queue_mode") or "").strip()
                existing_watch_tag = str(existing.get("watch_tag") or "").strip()
                started_at = str(existing.get("started_at") or "").strip()
                details = f" by PID {existing_pid}"
                if queue_mode:
                    details += f" ({queue_mode})"
                if existing_watch_tag:
                    details += f" tag={existing_watch_tag}"
                if started_at:
                    details += f" since {started_at}"
                return None, f"Watcher already running for {BATCH_DIR}{details}."
            try:
                lock_path.unlink()
            except OSError:
                pass
            continue
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        return lock_path, None
    return None, f"Could not acquire watcher lock: {lock_path}"


def release_watch_lock(lock_path: Path | None) -> None:
    if lock_path is None or not lock_path.exists():
        return
    payload = _read_json(lock_path, {})
    lock_pid = int(payload.get("pid") or 0)
    if lock_pid not in (0, os.getpid()):
        return
    try:
        lock_path.unlink()
    except OSError:
        pass


def run_command(
    command: list[str],
    *,
    env: dict[str, str] | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    process_env = os.environ.copy()
    if env:
        process_env.update(env)
    return subprocess.run(
        command,
        cwd=ROOT,
        check=check,
        capture_output=True,
        text=True,
        env=process_env,
    )


def spawn_command(
    label: str,
    command: list[str],
    *,
    env: dict[str, str] | None = None,
    log_dir: Path | None = None,
) -> int:
    selected_log_dir = log_dir or LOG_DIR
    selected_log_dir.mkdir(parents=True, exist_ok=True)
    log_path = selected_log_dir / f"{label}.log"
    process_env = os.environ.copy()
    if env:
        process_env.update(env)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(f"$ {' '.join(command)}\n")
        handle.flush()
        process = subprocess.Popen(
            command,
            cwd=ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            env=process_env,
            start_new_session=True,
        )
    return process.pid


def shutil_which(binary: str) -> str | None:
    for directory in os.environ.get("PATH", "").split(os.pathsep):
        candidate = Path(directory) / binary
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def _parse_nvidia_metric(raw_value: str) -> int | None:
    token = str(raw_value).strip()
    if not token or token.lower() in {"n/a", "[not supported]"}:
        return None
    match = re.search(r"-?\d+", token)
    if match is None:
        return None
    return int(match.group(0))


def normalize_gpu_devices(raw_value: str) -> list[str]:
    tokens = [item.strip() for item in raw_value.split(",") if item.strip()]
    if tokens:
        return tokens
    if not shutil_which("nvidia-smi"):
        return []
    result = run_command(["nvidia-smi", "--query-gpu=index", "--format=csv,noheader"], check=False)
    if result.returncode != 0:
        return []
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def gpu_compute_counts() -> dict[str, int]:
    if not shutil_which("nvidia-smi"):
        return {}
    gpu_result = run_command(
        ["nvidia-smi", "--query-gpu=index,uuid", "--format=csv,noheader"],
        check=False,
    )
    if gpu_result.returncode != 0:
        return {}
    uuid_to_index: dict[str, str] = {}
    for raw_line in gpu_result.stdout.splitlines():
        if not raw_line.strip():
            continue
        index, uuid = [item.strip() for item in raw_line.split(",", 1)]
        uuid_to_index[uuid] = index
    app_result = run_command(
        ["nvidia-smi", "--query-compute-apps=gpu_uuid,pid", "--format=csv,noheader"],
        check=False,
    )
    if app_result.returncode != 0:
        return {index: 0 for index in uuid_to_index.values()}
    counts: dict[str, int] = {index: 0 for index in uuid_to_index.values()}
    for raw_line in app_result.stdout.splitlines():
        if not raw_line.strip():
            continue
        uuid = raw_line.split(",", 1)[0].strip()
        index = uuid_to_index.get(uuid)
        if index is not None:
            counts[index] = counts.get(index, 0) + 1
    return counts


def gpu_device_stats() -> dict[str, dict[str, int | None]]:
    if not shutil_which("nvidia-smi"):
        return {}
    result = run_command(
        [
            "nvidia-smi",
            "--query-gpu=index,memory.used,memory.total,utilization.gpu,utilization.memory",
            "--format=csv,noheader",
        ],
        check=False,
    )
    if result.returncode != 0:
        return {}
    stats: dict[str, dict[str, int | None]] = {}
    for raw_line in result.stdout.splitlines():
        if not raw_line.strip():
            continue
        parts = [item.strip() for item in raw_line.split(",")]
        if len(parts) < 5:
            continue
        index = parts[0]
        stats[index] = {
            "memory_used_mb": _parse_nvidia_metric(parts[1]),
            "memory_total_mb": _parse_nvidia_metric(parts[2]),
            "gpu_utilization_percent": _parse_nvidia_metric(parts[3]),
            "memory_utilization_percent": _parse_nvidia_metric(parts[4]),
        }
    return stats


def _gpu_headroom_allows_launch(
    stats: dict[str, int | None],
    *,
    min_free_gpu_memory_mb: int,
    max_gpu_utilization: int,
) -> bool:
    checks: list[bool] = []
    if min_free_gpu_memory_mb > 0:
        memory_used_mb = stats.get("memory_used_mb")
        memory_total_mb = stats.get("memory_total_mb")
        if memory_used_mb is None or memory_total_mb is None:
            return False
        checks.append(max(memory_total_mb - memory_used_mb, 0) >= min_free_gpu_memory_mb)
    if max_gpu_utilization > 0:
        gpu_utilization_percent = stats.get("gpu_utilization_percent")
        if gpu_utilization_percent is None:
            return False
        checks.append(gpu_utilization_percent <= max_gpu_utilization)
    return bool(checks) and all(checks)


def available_gpu_devices(
    gpu_devices: list[str],
    *,
    gpu_counts: dict[str, int],
    max_compute_apps_per_gpu: int,
    gpu_stats: dict[str, dict[str, int | None]] | None = None,
    min_free_gpu_memory_mb: int = 0,
    max_gpu_utilization: int = 0,
) -> list[str]:
    if max_compute_apps_per_gpu <= 0:
        return []
    indexed_devices = list(enumerate(gpu_devices))
    stats_by_device = gpu_stats or {}
    available: list[tuple[tuple[int, int, int, int, int], str]] = []
    for position, device in indexed_devices:
        count = gpu_counts.get(device, 0)
        stats = stats_by_device.get(device, {})
        count_allowed = count < max_compute_apps_per_gpu
        headroom_allowed = False
        if not count_allowed:
            headroom_allowed = _gpu_headroom_allows_launch(
                stats,
                min_free_gpu_memory_mb=max(min_free_gpu_memory_mb, 0),
                max_gpu_utilization=max(max_gpu_utilization, 0),
            )
        if not count_allowed and not headroom_allowed:
            continue
        memory_used_mb = stats.get("memory_used_mb")
        memory_total_mb = stats.get("memory_total_mb")
        if memory_used_mb is not None and memory_total_mb is not None:
            free_memory_mb = max(memory_total_mb - memory_used_mb, 0)
        else:
            free_memory_mb = -1
        gpu_utilization_percent = stats.get("gpu_utilization_percent")
        util_sort = gpu_utilization_percent if gpu_utilization_percent is not None else 10**9
        availability_class = 0 if count_allowed else 1
        available.append(
            (
                (
                    availability_class,
                    count,
                    -free_memory_mb,
                    util_sort,
                    position,
                ),
                device,
            )
        )
    available.sort(key=lambda item: item[0])
    return [device for _sort_key, device in available]


def cpu_load_per_core() -> float | None:
    cpu_count = os.cpu_count() or 0
    if cpu_count <= 0 or not hasattr(os, "getloadavg"):
        return None
    try:
        load1, _load5, _load15 = os.getloadavg()
    except OSError:
        return None
    return load1 / cpu_count


def launch_allowed_by_cpu(*, max_load_per_core: float) -> tuple[bool, float | None]:
    load_per_core = cpu_load_per_core()
    if max_load_per_core <= 0 or load_per_core is None:
        return True, load_per_core
    return load_per_core < max_load_per_core, load_per_core


def _command_matches_repo_job(command: str) -> bool:
    return str(ROOT / "runs") in command and "/jobs/" in command


def _command_matches_batch(command: str) -> bool:
    return str(BATCH_DIR) in command


def _job_identity_from_text(text: str) -> dict[str, str] | None:
    if str(BATCH_DIR) not in text and "abag-rbfe" not in text:
        return None
    path_match = JOB_RE.search(text)
    if path_match:
        return {"job_id": path_match.group(1)}
    cli_match = CLI_JOB_RE.search(text)
    batch_dir_match = BATCH_DIR_RE.search(text)
    if cli_match and batch_dir_match:
        try:
            parsed_batch_dir = Path(batch_dir_match.group(1)).expanduser().resolve()
        except OSError:
            return None
        if parsed_batch_dir == BATCH_DIR:
            return {"job_id": cli_match.group(1)}
    return None


def _parse_mdrun_thread_count(mdrun_args: str) -> int:
    nt_total = 0
    ntmpi = 1
    ntomp = 1
    try:
        tokens = shlex.split(mdrun_args)
    except ValueError:
        tokens = mdrun_args.split()
    for index, token in enumerate(tokens):
        if token not in {"-nt", "-ntmpi", "-ntomp"}:
            continue
        if index + 1 >= len(tokens):
            continue
        try:
            value = int(tokens[index + 1])
        except ValueError:
            continue
        if value <= 0:
            continue
        if token == "-nt":
            nt_total = value
        elif token == "-ntmpi":
            ntmpi = value
        elif token == "-ntomp":
            ntomp = value
    if nt_total > 0:
        return nt_total
    return max(ntmpi * ntomp, 1)


def job_mdrun_threads(job_dir: Path, *, mdrun_args_override: str = "") -> int:
    override = (mdrun_args_override or "").strip()
    if override:
        return _parse_mdrun_thread_count(override)
    payload = _read_json(job_dir / "job_spec.json", {})
    protocol = payload.get("protocol") or {}
    return _parse_mdrun_thread_count(str(protocol.get("mdrun_args") or ""))


def _mdrun_threads_from_command(command: str) -> int:
    ntomp_match = NTOMP_RE.search(command)
    ntmpi_match = NTMPI_RE.search(command)
    ntomp = int(ntomp_match.group(1)) if ntomp_match else 1
    ntmpi = int(ntmpi_match.group(1)) if ntmpi_match else 1
    return max(ntomp, 1) * max(ntmpi, 1)


def active_mdrun_threads() -> tuple[int, int]:
    result = run_command(["ps", "-ef"], check=False)
    if result.returncode != 0:
        return 0, 0
    total_threads = 0
    process_count = 0
    for line in result.stdout.splitlines():
        if not _command_matches_repo_job(line) or not MDRUN_RE.search(line):
            continue
        total_threads += _mdrun_threads_from_command(line)
        process_count += 1
    return total_threads, process_count


def launch_allowed_by_thread_budget(*, max_active_mdrun_threads: int) -> tuple[bool, int, int]:
    active_threads, process_count = active_mdrun_threads()
    if max_active_mdrun_threads <= 0:
        return True, active_threads, process_count
    return active_threads < max_active_mdrun_threads, active_threads, process_count


def _active_mdrun_rows() -> list[tuple[int, int, str]]:
    result = run_command(["ps", "-eo", "pid,etimes,args"], check=False)
    if result.returncode != 0:
        return []
    rows: list[tuple[int, int, str]] = []
    for line in result.stdout.splitlines()[1:]:
        match = re.match(r"\s*(\d+)\s+(\d+)\s+(.*)$", line)
        if not match:
            continue
        pid = int(match.group(1))
        elapsed_seconds = int(match.group(2))
        command = match.group(3)
        if not _command_matches_repo_job(command) or not MDRUN_RE.search(command):
            continue
        rows.append((pid, elapsed_seconds, command))
    return rows


def _file_age_seconds(path: str) -> float | None:
    try:
        return max(time.time() - os.path.getmtime(path), 0.0)
    except FileNotFoundError:
        return None
    except OSError:
        return None


def _parse_mdrun_status(pid: int, elapsed_seconds: int, command: str) -> dict[str, Any]:
    deffnm = ""
    dhdl_path: str | None = None
    try:
        tokens = shlex.split(command)
    except ValueError:
        tokens = command.split()
    for index, token in enumerate(tokens):
        if token == "-deffnm" and index + 1 < len(tokens):
            deffnm = tokens[index + 1]
        elif token == "-dhdl" and index + 1 < len(tokens):
            dhdl_path = tokens[index + 1]
    log_path = f"{deffnm}.log" if deffnm else None
    progress_ages = [
        age
        for age in (
            _file_age_seconds(log_path) if log_path else None,
            _file_age_seconds(dhdl_path) if dhdl_path else None,
        )
        if age is not None
    ]
    progress_age = min(progress_ages) if progress_ages else None
    identity = _job_identity_from_text(deffnm or command)
    return {
        "pid": pid,
        "elapsed_seconds": elapsed_seconds,
        "elapsed_minutes": round(elapsed_seconds / 60.0, 1),
        "job_id": identity["job_id"] if identity is not None else "",
        "deffnm": deffnm,
        "deffnm_tail": "/".join(deffnm.split("/")[-4:]) if deffnm else "",
        "progress_age_seconds": progress_age,
    }


def active_mdrun_statuses() -> list[dict[str, Any]]:
    return [_parse_mdrun_status(pid, elapsed_seconds, command) for pid, elapsed_seconds, command in _active_mdrun_rows()]


def stale_mdrun_statuses(*, warn_stale_mdrun_seconds: int) -> list[dict[str, Any]]:
    if warn_stale_mdrun_seconds <= 0:
        return []
    stale: list[dict[str, Any]] = []
    for status in active_mdrun_statuses():
        progress_age = status.get("progress_age_seconds")
        if progress_age is None:
            if status["elapsed_seconds"] >= warn_stale_mdrun_seconds:
                stale.append(status)
            continue
        if progress_age >= warn_stale_mdrun_seconds:
            stale.append(status)
    stale.sort(key=lambda item: (-float(item.get("progress_age_seconds") or 0.0), item.get("job_id") or ""))
    return stale


def active_job_ids() -> set[str]:
    result = run_command(["ps", "-ef"], check=False)
    if result.returncode != 0:
        return set()
    jobs: set[str] = set()
    for line in result.stdout.splitlines():
        if not _command_matches_batch(line):
            continue
        identity = _job_identity_from_text(line)
        if identity is not None and any(marker in line for marker in ACTIVE_PATH_MARKERS):
            jobs.add(identity["job_id"])
            continue
        if identity is not None and CLI_JOB_RE.search(line):
            jobs.add(identity["job_id"])
    return jobs


def batch_plan_job_ids(batch_dir: Path) -> list[str]:
    payload = _read_json(batch_dir / "batch_plan.json", {})
    jobs = payload.get("jobs")
    if isinstance(jobs, list):
        ordered = [str(job.get("job_id") or "").strip() for job in jobs]
        return [job_id for job_id in ordered if job_id]
    return sorted(path.name for path in (batch_dir / "jobs").iterdir() if path.is_dir())


def queue_only_listed(job_ids: list[str], *, only_listed: bool, append_rest: bool) -> bool:
    if not job_ids:
        return False
    if append_rest:
        return False
    return True if only_listed or job_ids else False


def job_queue(job_ids: list[str], *, only_listed: bool = False) -> list[Path]:
    ordered_job_ids = batch_plan_job_ids(BATCH_DIR)
    explicit = [job_id for job_id in job_ids if (BATCH_DIR / "jobs" / job_id).is_dir()]
    if explicit and only_listed:
        return [BATCH_DIR / "jobs" / job_id for job_id in explicit]
    queue: list[str] = []
    seen: set[str] = set()
    for job_id in explicit:
        if job_id in seen:
            continue
        seen.add(job_id)
        queue.append(job_id)
    for job_id in ordered_job_ids:
        if job_id in seen:
            continue
        seen.add(job_id)
        queue.append(job_id)
    return [BATCH_DIR / "jobs" / job_id for job_id in queue if (BATCH_DIR / "jobs" / job_id).is_dir()]


def _job_charge_conserving(job_dir: Path) -> bool:
    payload = _read_json(job_dir / "job_spec.json", {})
    mutation_group = payload.get("mutation_group") or {}
    return bool(mutation_group.get("charge_conserving", False))


def _default_summary_row(job_id: str) -> dict[str, Any]:
    return {
        "job_id": job_id,
        "ddg_ready": False,
        "latest_stage": "",
        "latest_stage_state": "not_started",
        "analyzable": False,
        "resumable": True,
        "sample_completed_windows": 0,
        "sample_total_windows": 0,
        "equilibrate_completed_repeats": 0,
        "equilibrate_total_repeats": 0,
    }


def _job_has_remaining_progress(
    job_id: str,
    job_row: dict[str, Any],
    *,
    active_job_ids: set[str],
    recent_launches: dict[str, float],
) -> bool:
    latest_stage = str(job_row.get("latest_stage") or "").strip()
    latest_state = str(job_row.get("latest_stage_state") or "").strip()
    if _bool_env(job_row.get("ddg_ready")) or (latest_stage == "report" and latest_state == "completed"):
        return False
    if latest_state == "blocked_input":
        return False
    if job_id in active_job_ids or job_id in recent_launches:
        return True
    if _bool_env(job_row.get("analyzable")):
        return True
    return _bool_env(job_row.get("resumable", True))


def pending_charge_conserving_job_ids(
    queue: list[Path],
    summary_rows: dict[str, dict[str, Any]],
    *,
    active_job_ids: set[str],
    recent_launches: dict[str, float],
) -> list[str]:
    pending: list[str] = []
    for job_dir in queue:
        if not _job_charge_conserving(job_dir):
            continue
        job_id = job_dir.name
        row = summary_rows.get(job_id, _default_summary_row(job_id))
        if _job_has_remaining_progress(
            job_id,
            row,
            active_job_ids=active_job_ids,
            recent_launches=recent_launches,
        ):
            pending.append(job_id)
    return pending


def charge_conserving_launch_backlog_job_ids(
    queue: list[Path],
    summary_rows: dict[str, dict[str, Any]],
    *,
    active_job_ids: set[str],
    recent_launches: dict[str, float],
) -> list[str]:
    backlog: list[str] = []
    for job_dir in queue:
        if not _job_charge_conserving(job_dir):
            continue
        job_id = job_dir.name
        row = summary_rows.get(job_id, _default_summary_row(job_id))
        latest_stage = str(row.get("latest_stage") or "").strip()
        latest_state = str(row.get("latest_stage_state") or "").strip()
        if _bool_env(row.get("ddg_ready")) or (latest_stage == "report" and latest_state == "completed"):
            continue
        if latest_state == "blocked_input":
            continue
        if job_id in active_job_ids or job_id in recent_launches:
            continue
        if _bool_env(row.get("analyzable")):
            continue
        if _bool_env(row.get("resumable", True)):
            backlog.append(job_id)
    return backlog


def resumable_priority(job_dir: Path, job_row: dict[str, Any], *, queue_position: int) -> tuple[int, int, float, int, str]:
    charge_penalty = 0 if _job_charge_conserving(job_dir) else 1
    latest_stage = str(job_row.get("latest_stage") or "").strip()
    latest_state = str(job_row.get("latest_stage_state") or "").strip()
    sample_completed = int(job_row.get("sample_completed_windows") or 0)
    sample_total = max(int(job_row.get("sample_total_windows") or 0), 1)
    equilibrate_completed = int(job_row.get("equilibrate_completed_repeats") or 0)
    equilibrate_total = max(int(job_row.get("equilibrate_total_repeats") or 0), 1)

    if latest_stage == "sample":
        return (
            charge_penalty,
            0,
            -(sample_completed / sample_total),
            queue_position,
            job_dir.name,
        )
    if latest_stage == "equilibrate":
        return (
            charge_penalty,
            1,
            -(equilibrate_completed / equilibrate_total),
            queue_position,
            job_dir.name,
        )
    if latest_stage == "build_legs":
        return (charge_penalty, 2, 0.0, queue_position, job_dir.name)
    if latest_stage == "mutate":
        return (charge_penalty, 3, 0.0, queue_position, job_dir.name)
    if latest_stage == "prepare":
        return (charge_penalty, 4, 0.0, queue_position, job_dir.name)
    if latest_stage == "ingest":
        return (charge_penalty, 5, 0.0, queue_position, job_dir.name)
    if latest_state == "not_started":
        return (charge_penalty, 6, 0.0, queue_position, job_dir.name)
    return (charge_penalty, 7, 0.0, queue_position, job_dir.name)


def refresh_batch_summary(*, dry_run: bool) -> dict[str, Any]:
    command = [str(PYTHON_BIN), "-u", str(REPORT_SCRIPT), "--batch-dir", str(BATCH_DIR)]
    if dry_run:
        print("[watch] dry-run Patel refresh:", " ".join(command))
    else:
        result = run_command(command, check=False)
        if result.returncode not in {0, 2}:
            sys.stderr.write(f"[watch] Patel summary refresh exited code={result.returncode}: {' '.join(command)}\n")
            sys.stderr.write(result.stdout)
            sys.stderr.write(result.stderr)
    return _read_json(BATCH_DIR / "reports" / "batch_summary.json", {"jobs": []})


def run_post_refresh(*, dry_run: bool) -> None:
    if not POST_REPORT_REFRESH_COMMAND:
        return
    env = {
        "BATCH_DIR": str(BATCH_DIR),
        "PLAN_ROOT": str(BATCH_DIR.parent),
    }
    if dry_run:
        env_summary = " ".join(f"{key}={value}" for key, value in sorted(env.items()))
        print("[watch] dry-run post-refresh command:", env_summary, " ".join(POST_REPORT_REFRESH_COMMAND))
        return
    result = run_command(POST_REPORT_REFRESH_COMMAND, env=env, check=False)
    if result.returncode != 0:
        sys.stderr.write(
            f"[watch] post-refresh command exited code={result.returncode}: {' '.join(POST_REPORT_REFRESH_COMMAND)}\n"
        )
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)


def analyze_job(job_dir: Path, dry_run: bool) -> None:
    command = [str(ABAG_RBFE), "analyze", job_dir.name, "--batch-dir", str(BATCH_DIR), "--execute"]
    if dry_run:
        print("[watch] dry-run analyze:", " ".join(command))
        return
    print(f"[watch] analyze {job_dir.name}")
    result = run_command(command, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)


def resume_launch_environment(
    job_dir: Path,
    gpu_device: str,
    *,
    gpu_devices: list[str],
    gpu_counts: dict[str, int],
    max_compute_apps_per_gpu: int,
    mdrun_args_override: str = "",
) -> dict[str, str]:
    thread_count = max(job_mdrun_threads(job_dir, mdrun_args_override=mdrun_args_override), 1)
    env = {"CUDA_VISIBLE_DEVICES": gpu_device}
    try:
        gpu_position = gpu_devices.index(gpu_device)
    except ValueError:
        gpu_position = 0
    gpu_slot = max(gpu_counts.get(gpu_device, 0), 0)
    slots_per_gpu = max(max_compute_apps_per_gpu, 1)
    base_slot = gpu_position * slots_per_gpu + gpu_slot
    cpu_count = os.cpu_count() or 0
    slot_capacity = max(cpu_count // thread_count, 1) if cpu_count > 0 else 0
    if slot_capacity > 0:
        base_slot %= slot_capacity
    env["ABAG_RBFE_MDRUN_PINOFFSET"] = str(base_slot * thread_count)
    env["ABAG_RBFE_MDRUN_PINSTRIDE"] = "1"
    override = (mdrun_args_override or "").strip()
    if override:
        env["ABAG_RBFE_MDRUN_ARGS"] = override
    return env


def launch_resume(
    job_dir: Path,
    gpu_device: str,
    dry_run: bool,
    *,
    gpu_devices: list[str],
    gpu_counts: dict[str, int],
    max_compute_apps_per_gpu: int,
    mdrun_args_override: str = "",
) -> None:
    command = [str(ABAG_RBFE), "resume", job_dir.name, "--batch-dir", str(BATCH_DIR), "--execute"]
    launch_env = resume_launch_environment(
        job_dir,
        gpu_device,
        gpu_devices=gpu_devices,
        gpu_counts=gpu_counts,
        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
        mdrun_args_override=mdrun_args_override,
    )
    effective_threads = job_mdrun_threads(job_dir, mdrun_args_override=mdrun_args_override)
    if dry_run:
        print(
            f"[watch] dry-run resume gpu={gpu_device} pin_offset={launch_env['ABAG_RBFE_MDRUN_PINOFFSET']} "
            f"threads={effective_threads} mdrun_args_override={launch_env.get('ABAG_RBFE_MDRUN_ARGS', '')!r}:",
            " ".join(command),
        )
        return
    pid = spawn_command(
        f"{job_dir.name}_{int(time.time())}",
        command,
        env=launch_env,
        log_dir=LOG_DIR,
    )
    print(
        f"[watch] launched {job_dir.name} on GPU {gpu_device} "
        f"pin_offset={launch_env['ABAG_RBFE_MDRUN_PINOFFSET']} "
        f"threads={effective_threads} "
        f"mdrun_args_override={launch_env.get('ABAG_RBFE_MDRUN_ARGS', '')!r} (pid={pid})"
    )


def prune_recent_launches(recent_launches: dict[str, float], *, now_ts: float, cooldown_seconds: int) -> None:
    if cooldown_seconds <= 0:
        recent_launches.clear()
        return
    expired = [
        job_id
        for job_id, launched_at in recent_launches.items()
        if now_ts - launched_at >= cooldown_seconds
    ]
    for job_id in expired:
        recent_launches.pop(job_id, None)


def pass_once(
    job_ids: list[str],
    gpu_devices: list[str],
    *,
    only_listed: bool,
    skip_charge_changing: bool,
    max_compute_apps_per_gpu: int,
    min_free_gpu_memory_mb: int,
    max_gpu_utilization: int,
    max_load_per_core: float,
    max_active_mdrun_threads: int,
    warn_stale_mdrun_seconds: int,
    mdrun_args_override: str,
    max_launches_per_pass: int,
    launch_cooldown_seconds: int,
    recent_launches: dict[str, float],
    dry_run: bool,
    now_ts: float | None = None,
) -> None:
    current_time = time.time() if now_ts is None else now_ts
    prune_recent_launches(
        recent_launches,
        now_ts=current_time,
        cooldown_seconds=max(launch_cooldown_seconds, 0),
    )
    summary_payload = refresh_batch_summary(dry_run=dry_run)
    summary_rows = {
        str(row.get("job_id") or "").strip(): row
        for row in summary_payload.get("jobs", [])
        if str(row.get("job_id") or "").strip()
    }
    active = active_job_ids()
    gpu_counts = gpu_compute_counts()
    gpu_stats = gpu_device_stats()
    gpu_free_gpus = available_gpu_devices(
        gpu_devices,
        gpu_counts=gpu_counts,
        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
        gpu_stats=gpu_stats,
        min_free_gpu_memory_mb=max(min_free_gpu_memory_mb, 0),
        max_gpu_utilization=max(max_gpu_utilization, 0),
    )
    cpu_launch_allowed, load_per_core = launch_allowed_by_cpu(max_load_per_core=max_load_per_core)
    thread_launch_allowed, active_threads, active_mdrun_processes = launch_allowed_by_thread_budget(
        max_active_mdrun_threads=max_active_mdrun_threads,
    )
    free_gpus = gpu_free_gpus if cpu_launch_allowed and thread_launch_allowed else []
    print(
        f"[watch] batch={BATCH_DIR.name} active_jobs={sorted(active)} gpu_counts={gpu_counts} "
        f"gpu_stats={gpu_stats} "
        f"max_compute_apps_per_gpu={max_compute_apps_per_gpu} "
        f"min_free_gpu_memory_mb={min_free_gpu_memory_mb} "
        f"max_gpu_utilization={max_gpu_utilization} "
        f"max_load_per_core={max_load_per_core} cpu_load_per_core={load_per_core} "
        f"max_active_mdrun_threads={max_active_mdrun_threads} "
        f"active_mdrun_threads={active_threads} active_mdrun_processes={active_mdrun_processes} "
        f"free_gpus={free_gpus} skip_charge_changing={skip_charge_changing}"
    )
    if not cpu_launch_allowed and gpu_free_gpus:
        print(
            f"[watch] CPU launch gate active: load/core={load_per_core:.3f} "
            f">= threshold {max_load_per_core:.3f}; deferring new resumes."
        )
    if not thread_launch_allowed and gpu_free_gpus:
        print(
            f"[watch] mdrun thread gate active: active_threads={active_threads} "
            f">= threshold {max_active_mdrun_threads}; deferring new resumes."
        )
    stale_statuses = stale_mdrun_statuses(
        warn_stale_mdrun_seconds=max(warn_stale_mdrun_seconds, 0),
    )
    if stale_statuses:
        compact = [
            {
                "job_id": status["job_id"],
                "pid": status["pid"],
                "deffnm_tail": status["deffnm_tail"],
                "elapsed_minutes": status["elapsed_minutes"],
                "progress_age_seconds": round(float(status["progress_age_seconds"] or 0.0), 1),
            }
            for status in stale_statuses
        ]
        print(
            f"[watch] stale active mdrun processes (threshold={warn_stale_mdrun_seconds}s): "
            + json.dumps(compact, ensure_ascii=False)
        )

    analyzable: list[Path] = []
    resumable: list[tuple[Path, dict[str, Any]]] = []
    blocked: list[str] = []
    completed: list[str] = []
    cooling: list[str] = []
    skipped_charge: list[str] = []
    active_batch: list[str] = []

    queue = job_queue(job_ids, only_listed=only_listed)
    queue_positions = {job_dir.name: index for index, job_dir in enumerate(queue)}
    conserving_pending = pending_charge_conserving_job_ids(
        queue,
        summary_rows,
        active_job_ids=active,
        recent_launches=recent_launches,
    )
    conserving_launch_backlog = charge_conserving_launch_backlog_job_ids(
        queue,
        summary_rows,
        active_job_ids=active,
        recent_launches=recent_launches,
    )
    defer_charge_changing = skip_charge_changing and bool(conserving_launch_backlog)

    for job_dir in queue:
        job_id = job_dir.name
        row = summary_rows.get(job_id, _default_summary_row(job_id))
        if defer_charge_changing and not _job_charge_conserving(job_dir):
            skipped_charge.append(job_id)
            continue
        latest_stage = str(row.get("latest_stage") or "").strip()
        latest_state = str(row.get("latest_stage_state") or "").strip()
        if _bool_env(row.get("ddg_ready")) or (latest_stage == "report" and latest_state == "completed"):
            completed.append(job_id)
            continue
        if latest_state == "blocked_input":
            blocked.append(job_id)
            continue
        if job_id in active:
            active_batch.append(job_id)
            continue
        if job_id in recent_launches:
            cooling.append(job_id)
            continue
        if _bool_env(row.get("analyzable")):
            analyzable.append(job_dir)
            continue
        if not _bool_env(row.get("resumable", True)):
            blocked.append(job_id)
            continue
        resumable.append((job_dir, row))

    resumable.sort(key=lambda item: resumable_priority(item[0], item[1], queue_position=queue_positions.get(item[0].name, 9999)))

    print(
        "[watch] queue",
        json.dumps(
            {
                "completed": completed,
                "blocked": blocked,
                "active": active_batch,
                "analyzable": [job_dir.name for job_dir in analyzable],
                "cooling": cooling,
                "defer_charge_changing": defer_charge_changing,
                "charge_conserving_pending": conserving_pending,
                "charge_conserving_launch_backlog": conserving_launch_backlog,
                "skipped_charge_changing": skipped_charge,
                "resumable": [job_dir.name for job_dir, _row in resumable],
            },
            ensure_ascii=False,
        ),
    )

    changed = False
    for job_dir in analyzable:
        analyze_job(job_dir, dry_run=dry_run)
        changed = True

    launch_pairs = list(zip(free_gpus, [job_dir for job_dir, _row in resumable]))
    if max_launches_per_pass > 0:
        launch_pairs = launch_pairs[:max_launches_per_pass]

    for gpu_device, job_dir in launch_pairs:
        launch_resume(
            job_dir,
            gpu_device,
            dry_run=dry_run,
            gpu_devices=gpu_devices,
            gpu_counts=gpu_counts,
            max_compute_apps_per_gpu=max_compute_apps_per_gpu,
            mdrun_args_override=mdrun_args_override,
        )
        recent_launches[job_dir.name] = current_time
        changed = True

    if changed or active_batch:
        refresh_batch_summary(dry_run=dry_run)
        run_post_refresh(dry_run=dry_run)


def main() -> int:
    global BATCH_DIR, LOG_DIR, POST_REPORT_REFRESH_COMMAND
    args = parse_args()
    BATCH_DIR = Path(args.batch_dir).expanduser().resolve()
    LOG_DIR = BATCH_DIR / "reports" / "watch"
    post_refresh_command = str(getattr(args, "post_refresh_command", "") or "").strip()
    POST_REPORT_REFRESH_COMMAND = [post_refresh_command] if post_refresh_command else []

    gpu_devices = normalize_gpu_devices(args.gpu_devices)
    if not gpu_devices:
        print("[watch] warning: no GPU devices detected; resumable jobs will not be launched.")
    only_listed = queue_only_listed(
        args.job_id,
        only_listed=args.only_listed,
        append_rest=args.append_rest,
    )
    if args.job_id:
        queue_mode = "listed-only" if only_listed else "listed-plus-rest"
        print(
            f"[watch] explicit job scope mode={queue_mode} "
            f"job_ids={json.dumps(args.job_id, ensure_ascii=False)}"
        )
    recent_launches: dict[str, float] = {}
    lock_path, conflict = acquire_watch_lock(
        args.job_id,
        only_listed=only_listed,
        watch_tag=str(getattr(args, "watch_tag", "") or "").strip(),
    )
    if conflict is not None:
        print(f"[watch] {conflict}", file=sys.stderr)
        return 1

    try:
        while True:
            pass_once(
                args.job_id,
                gpu_devices,
                only_listed=only_listed,
                skip_charge_changing=bool(args.skip_charge_changing),
                max_compute_apps_per_gpu=max(args.max_compute_apps_per_gpu, 0),
                min_free_gpu_memory_mb=max(args.min_free_gpu_memory_mb, 0),
                max_gpu_utilization=max(args.max_gpu_utilization, 0),
                max_load_per_core=max(args.max_load_per_core, 0.0),
                max_active_mdrun_threads=max(args.max_active_mdrun_threads, 0),
                warn_stale_mdrun_seconds=max(args.warn_stale_mdrun_seconds, 0),
                mdrun_args_override=args.mdrun_args_override,
                max_launches_per_pass=max(args.max_launches_per_pass, 0),
                launch_cooldown_seconds=max(args.launch_cooldown_seconds, 0),
                recent_launches=recent_launches,
                dry_run=args.dry_run,
            )
            if args.once:
                return 0
            time.sleep(max(args.poll_seconds, 5))
    finally:
        release_watch_lock(lock_path)


if __name__ == "__main__":
    raise SystemExit(main())
