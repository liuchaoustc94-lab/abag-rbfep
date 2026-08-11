#!/usr/bin/env python3
"""Monitor the validation priority queue and keep GPUs filled."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
from types import SimpleNamespace
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
ABAG_RBFE = ROOT / ".venv" / "bin" / "abag-rbfe"
BENCHMARK_ROOT = ROOT / "benchmarks" / "ab_bind"
SPLIT_FILE = BENCHMARK_ROOT / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
SPLIT_NAME = "validation"
RUNS_ROOT = ROOT / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
LOG_DIR = RUNS_ROOT / "reports" / "watch"
REPORTS_DIR = RUNS_ROOT / "reports"
PLAN_ROOTS: list[Path] | None = None
MERGED_PLAN_ROOT: Path | None = None
MERGED_EXTRA_PLAN_ROOTS: list[Path] = []
POST_REPORT_REFRESH_COMMAND: list[str] = []
WATCH_LOCK_NAME = "watch.lock.json"
REPORT_REFRESH_LOCK_NAME = "report_refresh.lock.json"
LAUNCH_COORDINATION_LOCK_NAME = "launch_coordination.lock.json"
WATCH_SUPERVISOR_NO_RESTART_CODE = int(os.environ.get("WATCH_SUPERVISOR_NO_RESTART_CODE", "75"))

DEFAULT_QUEUE = [
    "3hfm-antibody-h-y50a",
    "3hfm-antibody-h-y33a",
    "3hfm-antibody-h-c95a",
    "1cz8-antigen-w-g92a",
    "1bj1-antigen-w-g92a",
    "3nps-antigen-a-h138a",
    "3hfm-antibody-l-n31a",
    "3hfm-antibody-l-n32a",
    "3hfm-antigen-y-y20a",
    "1cz8-antigen-w-m81a",
    "1bj1-antigen-w-i83a",
    "3nps-antigen-a-y141a",
]

JOB_RE = re.compile(r"/jobs/([^/]+)/")
CLI_JOB_RE = re.compile(r"\babag-rbfe\b\s+(?:run|resume|analyze)\s+([^\s]+)")
BATCH_DIR_RE = re.compile(r"--batch-dir\s+([^\s]+)")
ACTIVE_PATH_MARKERS = ("/artifacts/commands/", "/legs/")
NTOMP_RE = re.compile(r"(?:^|\s)-ntomp\s+(\d+)(?:\s|$)")
NTMPI_RE = re.compile(r"(?:^|\s)-ntmpi\s+(\d+)(?:\s|$)")
MDRUN_RE = re.compile(r"(?:^|\s)\S*gmx(?:_mpi)?\s+mdrun(?:\s|$)")
VALID_QUEUE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
JOB_DIR_PATH_RE = re.compile(r"(/[^\s]+/jobs/[^/\s]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_id", nargs="*")
    parser.add_argument(
        "--plan-root",
        default=os.environ.get("PLAN_ROOT", str(RUNS_ROOT)),
        help="AB-Bind plan root to monitor. Defaults to the validation priority plan root.",
    )
    parser.add_argument(
        "--extra-plan-root",
        action="append",
        default=[],
        help="Additional plan root to monitor with lower priority than the primary --plan-root.",
    )
    parser.add_argument(
        "--merged-plan-root",
        default=os.environ.get("MERGED_PLAN_ROOT", ""),
        help="Optional primary plan root for a merged report refresh after per-root refreshes.",
    )
    parser.add_argument(
        "--merged-extra-plan-root",
        action="append",
        default=[],
        help="Additional plan roots passed to the merged report refresh.",
    )
    parser.add_argument(
        "--post-refresh-command",
        default=os.environ.get("POST_REFRESH_COMMAND", ""),
        help="Optional command executed after report refreshes.",
    )
    parser.add_argument(
        "--split-name",
        default=os.environ.get("SPLIT_NAME", SPLIT_NAME),
        help="Selection split name passed to report-abbind refreshes.",
    )
    parser.add_argument(
        "--split-file",
        default=os.environ.get("SPLIT_FILE", str(SPLIT_FILE)),
        help="Selection split file passed to report-abbind refreshes.",
    )
    parser.add_argument("--poll-seconds", type=int, default=int(os.environ.get("POLL_SECONDS", "60")))
    parser.add_argument("--gpu-devices", default=os.environ.get("GPU_DEVICES", ""))
    parser.add_argument(
        "--max-compute-apps-per-gpu",
        type=int,
        default=int(os.environ.get("MAX_COMPUTE_APPS_PER_GPU", "1")),
    )
    parser.add_argument(
        "--min-free-gpu-memory-mb",
        type=int,
        default=int(os.environ.get("MIN_FREE_GPU_MEMORY_MB", "0")),
        help=(
            "Optional GPU headroom override: when a device already reached "
            "--max-compute-apps-per-gpu, still allow launches if the device keeps "
            "at least this much free GPU memory. 0 disables the override."
        ),
    )
    parser.add_argument(
        "--max-gpu-utilization",
        type=int,
        default=int(os.environ.get("MAX_GPU_UTILIZATION", "0")),
        help=(
            "Optional GPU headroom override companion gate: when set, the "
            "headroom override only applies if current GPU utilization is at or "
            "below this percentage. 0 disables the utilization check."
        ),
    )
    parser.add_argument(
        "--max-load-per-core",
        type=float,
        default=float(os.environ.get("MAX_LOAD_PER_CORE", "0.95")),
    )
    parser.add_argument(
        "--launch-cooldown-seconds",
        type=int,
        default=int(os.environ.get("LAUNCH_COOLDOWN_SECONDS", "180")),
    )
    parser.add_argument(
        "--max-active-mdrun-threads",
        type=int,
        default=int(os.environ.get("MAX_ACTIVE_MDRUN_THREADS", "0")),
        help="Optional hard cap on the summed -ntmpi * -ntomp of active GROMACS mdrun processes. 0 disables the gate.",
    )
    parser.add_argument(
        "--thread-budget-plan-root",
        action="append",
        default=[],
        help=(
            "Optional plan root scope for the mdrun thread budget. When omitted, the budget "
            "counts active mdrun threads across every monitored plan root."
        ),
    )
    parser.add_argument(
        "--warn-stale-mdrun-seconds",
        type=int,
        default=int(os.environ.get("WARN_STALE_MDRUN_SECONDS", "900")),
        help="Warn when an active mdrun has not updated its log/dhdl for this many seconds. 0 disables stale warnings.",
    )
    parser.add_argument(
        "--mdrun-args-override",
        default=os.environ.get("ABAG_RBFE_MDRUN_ARGS", os.environ.get("MDRUN_ARGS_OVERRIDE", "")),
        help="Optional mdrun argument override injected into future resumed jobs, for example '-ntmpi 1 -ntomp 2'.",
    )
    parser.add_argument(
        "--max-launches-per-pass",
        type=int,
        default=int(os.environ.get("MAX_LAUNCHES_PER_PASS", "0")),
        help="Maximum number of resumable jobs to launch per watcher pass. 0 disables the cap.",
    )
    parser.add_argument(
        "--allow-active-elsewhere-job-ids",
        action="store_true",
        default=os.environ.get("ALLOW_ACTIVE_ELSEWHERE_JOB_IDS", "0").strip().lower() in {"1", "true", "yes", "on"},
        help=(
            "Allow launching a job even when the same canonical job ID is already active in a different plan root. "
            "Same-root duplicate protection still uses root-sensitive job keys."
        ),
    )
    parser.add_argument(
        "--max-active-copies-per-job-id",
        type=int,
        default=int(os.environ.get("MAX_ACTIVE_COPIES_PER_JOB_ID", "0")),
        help=(
            "Optional cap on the number of concurrently active copies of the same canonical "
            "job_id across the repository. 0 disables the cap."
        ),
    )
    parser.add_argument(
        "--watch-tag",
        default=os.environ.get("WATCH_TAG", "").strip(),
        help="Optional watcher tag used only for process identification and lock diagnostics.",
    )
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    queue_scope = parser.add_mutually_exclusive_group()
    queue_scope.add_argument(
        "--only-listed",
        action="store_true",
        help=(
            "Restrict the queue to the explicitly listed job IDs. This is the default "
            "whenever positional job IDs are supplied."
        ),
    )
    queue_scope.add_argument(
        "--append-rest",
        action="store_true",
        help="After the explicitly listed job IDs, append the remainder of the plan root queue.",
    )
    parser.add_argument("--wait-for-pid", type=int, default=int(os.environ.get("WAIT_FOR_PID", "0") or "0"))
    return parser.parse_args()


def all_plan_roots() -> list[Path]:
    roots = PLAN_ROOTS or [RUNS_ROOT]
    unique: list[Path] = []
    for root in roots:
        path = Path(root)
        if path not in unique:
            unique.append(path)
    return unique


def _resolve_path_list(raw_items: list[str]) -> list[Path]:
    resolved_paths: list[Path] = []
    for item in raw_items:
        token = str(item or "").strip()
        if not token:
            continue
        resolved = Path(token).expanduser().resolve()
        if resolved not in resolved_paths:
            resolved_paths.append(resolved)
    return resolved_paths


def queue_only_listed(job_ids: list[str], *, only_listed: bool, append_rest: bool) -> bool:
    if not job_ids:
        return False
    if append_rest:
        return False
    return True if only_listed or job_ids else False


def multi_plan_mode() -> bool:
    return len(all_plan_roots()) > 1


def _watch_lock_path(watch_tag: str = "") -> Path:
    normalized_tag = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(watch_tag or "").strip()).strip("-.")
    if not normalized_tag:
        return LOG_DIR / WATCH_LOCK_NAME
    return LOG_DIR / f"watch.{normalized_tag}.lock.json"


def _report_refresh_lock_path() -> Path:
    anchor = Path(MERGED_PLAN_ROOT or RUNS_ROOT).parent
    return anchor / REPORT_REFRESH_LOCK_NAME


def _launch_coordination_lock_path() -> Path:
    roots = all_plan_roots()
    anchor_root = Path(roots[0] if roots else RUNS_ROOT)
    return anchor_root.parent / LAUNCH_COORDINATION_LOCK_NAME


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


def acquire_report_refresh_lock() -> tuple[Path | None, str | None]:
    lock_path = _report_refresh_lock_path()
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "pid": os.getpid(),
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_root": str(RUNS_ROOT),
        "merged_plan_root": str(MERGED_PLAN_ROOT) if MERGED_PLAN_ROOT is not None else "",
        "split_name": SPLIT_NAME,
    }
    for _attempt in range(2):
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            try:
                existing = json.loads(lock_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                existing = {}
            existing_pid = int(existing.get("pid") or 0)
            if existing_pid and existing_pid != os.getpid() and _pid_is_alive(existing_pid):
                started_at = str(existing.get("started_at") or "").strip()
                details = f"report refresh already running by PID {existing_pid}"
                if started_at:
                    details += f" since {started_at}"
                return None, details
            try:
                lock_path.unlink()
            except OSError:
                pass
            continue
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        return lock_path, None
    return None, f"Could not acquire report-refresh lock: {lock_path}"


def release_report_refresh_lock(lock_path: Path | None) -> None:
    if lock_path is None or not lock_path.exists():
        return
    try:
        payload = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = {}
    lock_pid = int(payload.get("pid") or 0)
    if lock_pid not in (0, os.getpid()):
        return
    try:
        lock_path.unlink()
    except OSError:
        pass


def acquire_launch_coordination_lock() -> tuple[Path | None, str | None]:
    lock_path = _launch_coordination_lock_path()
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "pid": os.getpid(),
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_roots": [str(root) for root in all_plan_roots()],
        "split_name": SPLIT_NAME,
    }
    for _attempt in range(2):
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            try:
                existing = json.loads(lock_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                existing = {}
            existing_pid = int(existing.get("pid") or 0)
            if existing_pid and existing_pid != os.getpid() and _pid_is_alive(existing_pid):
                started_at = str(existing.get("started_at") or "").strip()
                details = f"launch coordination already held by PID {existing_pid}"
                if started_at:
                    details += f" since {started_at}"
                return None, details
            try:
                lock_path.unlink()
            except OSError:
                pass
            continue
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        return lock_path, None
    return None, f"Could not acquire launch-coordination lock: {lock_path}"


def release_launch_coordination_lock(lock_path: Path | None) -> None:
    if lock_path is None or not lock_path.exists():
        return
    try:
        payload = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = {}
    lock_pid = int(payload.get("pid") or 0)
    if lock_pid not in (0, os.getpid()):
        return
    try:
        lock_path.unlink()
    except OSError:
        pass


def acquire_watch_lock(
    job_ids: list[str],
    *,
    only_listed: bool,
    watch_tag: str = "",
) -> tuple[Path | None, str | None]:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    lock_path = _watch_lock_path(watch_tag)
    payload = {
        "pid": os.getpid(),
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "plan_root": str(RUNS_ROOT),
        "job_ids": list(job_ids),
        "queue_mode": "listed-only" if only_listed else "listed-plus-rest",
        "watch_tag": watch_tag,
    }
    for _attempt in range(2):
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            try:
                existing = json.loads(lock_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                existing = {}
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
                return None, f"Watcher already running for {RUNS_ROOT}{details}."
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
    try:
        payload = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = {}
    lock_pid = int(payload.get("pid") or 0)
    if lock_pid not in (0, os.getpid()):
        return
    try:
        lock_path.unlink()
    except OSError:
        pass


def plan_root_for_job_dir(job_dir: Path) -> Path:
    for root in all_plan_roots():
        try:
            job_dir.relative_to(root)
        except ValueError:
            continue
        return root
    return job_dir.parent.parent.parent


def batch_key_from_parts(plan_root: Path, batch_id: str) -> str:
    if not multi_plan_mode():
        return batch_id
    return f"{plan_root.name}/{batch_id}"


def job_key_from_parts(plan_root: Path, batch_id: str, job_id: str) -> str:
    if not multi_plan_mode():
        return job_id
    return f"{plan_root.name}/{batch_id}/{job_id}"


def batch_cache_key(batch_dir: Path) -> str:
    return batch_key_from_parts(batch_dir.parent, batch_dir.name)


def job_cache_key(job_dir: Path) -> str:
    batch_dir = job_dir.parent.parent
    return job_key_from_parts(plan_root_for_job_dir(job_dir), batch_dir.name, job_dir.name)


def job_display_name(job_dir: Path) -> str:
    return job_cache_key(job_dir)


def _command_matches_plan_root(command: str) -> bool:
    return _command_matches_plan_roots(command)


def _command_matches_plan_roots(command: str, plan_roots: list[Path] | None = None) -> bool:
    roots = all_plan_roots() if plan_roots is None else list(plan_roots)
    return any(str(root) in command for root in roots)


def _command_matches_repo_job(command: str) -> bool:
    return str(ROOT / "runs") in command and "/jobs/" in command


def _canonical_job_id_from_text(text: str) -> str | None:
    path_match = JOB_RE.search(text)
    if path_match:
        job_id = path_match.group(1)
        if VALID_QUEUE_ID_RE.fullmatch(job_id):
            return job_id
    cli_match = CLI_JOB_RE.search(text)
    if cli_match:
        job_id = cli_match.group(1)
        if VALID_QUEUE_ID_RE.fullmatch(job_id):
            return job_id
    return None


def _active_copy_key_from_text(text: str) -> str | None:
    job_dir_match = JOB_DIR_PATH_RE.search(text)
    if job_dir_match:
        return job_dir_match.group(1)
    cli_match = CLI_JOB_RE.search(text)
    batch_dir_match = BATCH_DIR_RE.search(text)
    if cli_match and batch_dir_match:
        batch_dir = batch_dir_match.group(1).strip()
        job_id = cli_match.group(1).strip()
        if batch_dir and job_id:
            return f"{batch_dir}/jobs/{job_id}"
    return None


def _job_identity_from_text(text: str) -> dict[str, Any] | None:
    for plan_root in all_plan_roots():
        root_str = str(plan_root)
        if root_str not in text:
            continue
        suffix = text.split(root_str, 1)[1]
        path_match = re.search(r"/([^/\s]+)/jobs/([^/\s]+)", suffix)
        if path_match:
            batch_id, job_id = path_match.groups()
            if not (
                VALID_QUEUE_ID_RE.fullmatch(batch_id)
                and VALID_QUEUE_ID_RE.fullmatch(job_id)
            ):
                continue
            return {
                "plan_root": plan_root,
                "batch_id": batch_id,
                "job_id": job_id,
                "job_key": job_key_from_parts(plan_root, batch_id, job_id),
            }
        cli_match = CLI_JOB_RE.search(text)
        batch_dir_match = BATCH_DIR_RE.search(text)
        if cli_match and batch_dir_match:
            batch_dir = Path(batch_dir_match.group(1))
            try:
                batch_dir.relative_to(plan_root)
            except ValueError:
                continue
            job_id = cli_match.group(1)
            if not (
                VALID_QUEUE_ID_RE.fullmatch(batch_dir.name)
                and VALID_QUEUE_ID_RE.fullmatch(job_id)
            ):
                continue
            return {
                "plan_root": plan_root,
                "batch_id": batch_dir.name,
                "job_id": job_id,
                "job_key": job_key_from_parts(plan_root, batch_dir.name, job_id),
            }
    cli_match = CLI_JOB_RE.search(text)
    if cli_match and not multi_plan_mode():
        job_id = cli_match.group(1)
        if not VALID_QUEUE_ID_RE.fullmatch(job_id):
            return None
        return {
            "plan_root": RUNS_ROOT,
            "batch_id": "",
            "job_id": job_id,
            "job_key": job_id,
        }
    return None


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


def gpu_uuid_to_index() -> dict[str, str]:
    if not shutil_which("nvidia-smi"):
        return {}
    result = run_command(
        ["nvidia-smi", "--query-gpu=index,uuid", "--format=csv,noheader"],
        check=False,
    )
    if result.returncode != 0:
        return {}
    mapping: dict[str, str] = {}
    for raw_line in result.stdout.splitlines():
        if not raw_line.strip():
            continue
        index, uuid = [item.strip() for item in raw_line.split(",", 1)]
        mapping[uuid] = index
    return mapping


def gpu_pid_to_index() -> dict[int, str]:
    uuid_to_index = gpu_uuid_to_index()
    if not uuid_to_index:
        return {}
    result = run_command(
        ["nvidia-smi", "--query-compute-apps=gpu_uuid,pid", "--format=csv,noheader"],
        check=False,
    )
    if result.returncode != 0:
        return {}
    mapping: dict[int, str] = {}
    for raw_line in result.stdout.splitlines():
        if not raw_line.strip():
            continue
        uuid, raw_pid = [item.strip() for item in raw_line.split(",", 1)]
        if not raw_pid.isdigit():
            continue
        index = uuid_to_index.get(uuid)
        if index is not None:
            mapping[int(raw_pid)] = index
    return mapping


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


def _mdrun_threads_from_command(command: str) -> int:
    ntomp_match = NTOMP_RE.search(command)
    ntmpi_match = NTMPI_RE.search(command)
    ntomp = int(ntomp_match.group(1)) if ntomp_match else 1
    ntmpi = int(ntmpi_match.group(1)) if ntmpi_match else 1
    return max(ntomp, 1) * max(ntmpi, 1)


def active_mdrun_threads(*, plan_roots: list[Path] | None = None) -> tuple[int, int]:
    result = run_command(["ps", "-ef"], check=False)
    if result.returncode != 0:
        return 0, 0
    total_threads = 0
    process_count = 0
    for line in result.stdout.splitlines():
        if not _command_matches_plan_roots(line, plan_roots) or not MDRUN_RE.search(line):
            continue
        total_threads += _mdrun_threads_from_command(line)
        process_count += 1
    return total_threads, process_count


def active_mdrun_processes() -> list[dict[str, Any]]:
    result = run_command(["ps", "-ef"], check=False)
    if result.returncode != 0:
        return []
    processes: list[dict[str, Any]] = []
    for line in result.stdout.splitlines():
        if not _command_matches_plan_root(line) or not MDRUN_RE.search(line):
            continue
        parts = line.split(None, 7)
        if len(parts) < 8 or not parts[1].isdigit():
            continue
        pid = int(parts[1])
        command = parts[7]
        identity = _job_identity_from_text(command)
        processes.append(
            {
                "pid": pid,
                "command": command,
                "job_id": identity["job_id"] if identity is not None else "",
                "job_key": identity["job_key"] if identity is not None else "",
                "thread_count": _mdrun_threads_from_command(command),
            }
        )
    processes.sort(key=lambda item: (str(item.get("job_key") or item.get("job_id") or ""), int(item.get("pid") or 0)))
    return processes


def launch_allowed_by_thread_budget(
    *,
    max_active_mdrun_threads: int,
    plan_roots: list[Path] | None = None,
) -> tuple[bool, int, int]:
    if max_active_mdrun_threads <= 0:
        return True, 0, 0
    active_threads, process_count = active_mdrun_threads(plan_roots=plan_roots)
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
        if not _command_matches_plan_root(command) or not MDRUN_RE.search(command):
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
        "job_key": identity["job_key"] if identity is not None else "",
        "deffnm": deffnm,
        "deffnm_tail": "/".join(deffnm.split("/")[-4:]) if deffnm else "",
        "log_path": log_path,
        "log_age_seconds": _file_age_seconds(log_path) if log_path else None,
        "dhdl_path": dhdl_path,
        "dhdl_age_seconds": _file_age_seconds(dhdl_path) if dhdl_path else None,
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
        if not _command_matches_plan_root(line):
            continue
        identity = _job_identity_from_text(line)
        if identity is not None and any(marker in line for marker in ACTIVE_PATH_MARKERS):
            jobs.add(identity["job_key"])
            continue
        if identity is not None and CLI_JOB_RE.search(line):
            jobs.add(identity["job_key"])
    return jobs


def active_canonical_job_ids() -> set[str]:
    result = run_command(["ps", "-ef"], check=False)
    if result.returncode != 0:
        return set()
    jobs: set[str] = set()
    for line in result.stdout.splitlines():
        if not (_command_matches_repo_job(line) or "abag-rbfe" in line):
            continue
        if not (any(marker in line for marker in ACTIVE_PATH_MARKERS) or CLI_JOB_RE.search(line)):
            continue
        canonical_job_id = _canonical_job_id_from_text(line)
        if canonical_job_id:
            jobs.add(canonical_job_id)
    return jobs


def active_canonical_job_copy_counts() -> dict[str, int]:
    result = run_command(["ps", "-ef"], check=False)
    if result.returncode != 0:
        return {}
    copy_keys_by_job_id: dict[str, set[str]] = {}
    for line in result.stdout.splitlines():
        if not (_command_matches_repo_job(line) or "abag-rbfe" in line):
            continue
        if not (any(marker in line for marker in ACTIVE_PATH_MARKERS) or CLI_JOB_RE.search(line)):
            continue
        canonical_job_id = _canonical_job_id_from_text(line)
        copy_key = _active_copy_key_from_text(line)
        if not canonical_job_id or not copy_key:
            continue
        copy_keys_by_job_id.setdefault(canonical_job_id, set()).add(copy_key)
    return {job_id: len(copy_keys) for job_id, copy_keys in copy_keys_by_job_id.items()}


def resolve_job_dir(job_id: str) -> Path | None:
    matches = resolve_job_dirs(job_id)
    return matches[0] if matches else None


def resolve_job_dirs(job_id: str) -> list[Path]:
    matches: list[Path] = []
    for plan_root in all_plan_roots():
        for batch_dir in sorted(path for path in plan_root.iterdir() if path.is_dir() and (path / "jobs").is_dir()):
            candidate = batch_dir / "jobs" / job_id
            if candidate.is_dir():
                matches.append(candidate)
    return matches


def all_job_dirs() -> list[Path]:
    jobs: list[Path] = []
    for plan_root in all_plan_roots():
        for batch_dir in sorted(path for path in plan_root.iterdir() if path.is_dir() and (path / "jobs").is_dir()):
            jobs.extend(sorted(path for path in (batch_dir / "jobs").iterdir() if path.is_dir()))
    return jobs


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def report_priority_data() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    job_rows: dict[str, dict[str, str]] = {}
    batch_rows: dict[str, dict[str, str]] = {}
    for plan_root in all_plan_roots():
        reports_dir = plan_root / "reports"
        for row in _read_csv_rows(reports_dir / "plan_jobs.csv"):
            job_id = row.get("job_id", "")
            batch_id = row.get("batch_id", "")
            if not job_id or not batch_id:
                continue
            job_rows[job_key_from_parts(plan_root, batch_id, job_id)] = row
        for row in _read_csv_rows(reports_dir / "plan_batches.csv"):
            batch_id = row.get("batch_id", "")
            if not batch_id:
                continue
            batch_rows[batch_key_from_parts(plan_root, batch_id)] = row
    return job_rows, batch_rows


def merged_priority_job_rows() -> dict[str, dict[str, str]]:
    if MERGED_PLAN_ROOT is None:
        return {}
    merged_dir = MERGED_PLAN_ROOT / "reports" / "merged"
    candidates: list[Path] = []
    canonical_plan_jobs = merged_dir / "plan_jobs.csv"
    if canonical_plan_jobs.exists():
        candidates.append(canonical_plan_jobs)
    selections_dir = merged_dir / "selections"
    if selections_dir.exists():
        candidates.extend(sorted(selections_dir.glob("*/plan_jobs.csv")))
    if not candidates:
        return {}

    def _split_match_score(plan_jobs_path: Path) -> int:
        split_name = str(SPLIT_NAME or "").strip()
        if not split_name:
            return 0
        summary_path = plan_jobs_path.with_name("plan_summary.json")
        if not summary_path.exists():
            return 0
        try:
            payload = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return 0
        selection = payload.get("selection")
        if not isinstance(selection, dict):
            return 0
        return 1 if str(selection.get("split_name", "") or "").strip() == split_name else 0

    merged_plan_jobs = max(
        candidates,
        key=lambda path: (
            _split_match_score(path),
            path.stat().st_mtime_ns,
            1 if path == canonical_plan_jobs else 0,
            str(path),
        ),
    )
    merged_rows: dict[str, dict[str, str]] = {}
    for row in _read_csv_rows(merged_plan_jobs):
        job_id = row.get("job_id", "")
        if not job_id:
            continue
        merged_rows[job_id] = row
    return merged_rows


def read_stage_states(job_dir: Path) -> dict[str, str]:
    states: dict[str, str] = {}
    for stage_path in (job_dir / "stages").glob("*.json"):
        try:
            payload = json.loads(stage_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            continue
        states[stage_path.stem] = str(payload.get("state", ""))
    return states


def refresh_reports(dry_run: bool) -> None:
    def _format_command(command: list[str]) -> str:
        return " ".join(shlex.quote(str(token)) for token in command)

    def _run_refresh_command(
        command: list[str],
        *,
        description: str,
        env: dict[str, str] | None = None,
    ):
        rendered = _format_command(command)
        if dry_run:
            print(f"[watch] dry-run {description}: {rendered}")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        print(f"[watch] {description} start: {rendered}")
        result = run_command(command, env=env, check=False)
        print(f"[watch] {description} done rc={result.returncode}")
        return result

    report_refresh_lock: Path | None = None
    if not dry_run:
        report_refresh_lock, lock_message = acquire_report_refresh_lock()
        if report_refresh_lock is None:
            print(f"[watch] skipping report refresh: {lock_message}")
            return
    refreshed_roots: list[Path] = []
    try:
        for plan_root in all_plan_roots():
            path = Path(plan_root)
            if path in refreshed_roots:
                continue
            refreshed_roots.append(path)
            command = [
                str(ABAG_RBFE),
                "batch",
                "report-abbind",
                "--plan-root",
                str(path),
            ]
            if SPLIT_NAME:
                command.extend(["--split-name", SPLIT_NAME])
            if str(SPLIT_FILE):
                command.extend(["--split-file", str(SPLIT_FILE)])
            _run_refresh_command(command, description=f"report refresh plan_root={path}")
        if MERGED_PLAN_ROOT is not None:
            for extra_root in MERGED_EXTRA_PLAN_ROOTS:
                path = Path(extra_root)
                if path in refreshed_roots:
                    continue
                refreshed_roots.append(path)
                command = [
                    str(ABAG_RBFE),
                    "batch",
                    "report-abbind",
                    "--plan-root",
                    str(path),
                ]
                if SPLIT_NAME:
                    command.extend(["--split-name", SPLIT_NAME])
                if str(SPLIT_FILE):
                    command.extend(["--split-file", str(SPLIT_FILE)])
                _run_refresh_command(command, description=f"extra report refresh plan_root={path}")
            command = [
                str(ABAG_RBFE),
                "batch",
                "report-abbind",
                "--plan-root",
                str(MERGED_PLAN_ROOT),
            ]
            for extra_root in MERGED_EXTRA_PLAN_ROOTS:
                command.extend(["--extra-plan-root", str(extra_root)])
            if SPLIT_NAME:
                command.extend(["--split-name", SPLIT_NAME])
            if str(SPLIT_FILE):
                command.extend(["--split-file", str(SPLIT_FILE)])
            extra_roots_text = os.pathsep.join(str(Path(root)) for root in MERGED_EXTRA_PLAN_ROOTS)
            _run_refresh_command(
                command,
                description=(
                    f"merged report refresh plan_root={Path(MERGED_PLAN_ROOT)}"
                    + (f" extra_plan_roots={extra_roots_text}" if extra_roots_text else "")
                ),
            )
    finally:
        if not dry_run:
            release_report_refresh_lock(report_refresh_lock)
            report_refresh_lock = None

    if not POST_REPORT_REFRESH_COMMAND:
        return
    post_refresh_root = Path(MERGED_PLAN_ROOT or RUNS_ROOT)
    post_refresh_env = {"PLAN_ROOT": str(post_refresh_root)}
    if MERGED_PLAN_ROOT is not None:
        post_refresh_env["MERGED_PLAN_ROOT"] = str(MERGED_PLAN_ROOT)
    if MERGED_EXTRA_PLAN_ROOTS:
        post_refresh_env["MERGED_EXTRA_PLAN_ROOTS"] = os.pathsep.join(str(Path(root)) for root in MERGED_EXTRA_PLAN_ROOTS)
    env_summary = " ".join(f"{key}={value}" for key, value in sorted(post_refresh_env.items()))
    try:
        result = _run_refresh_command(
            POST_REPORT_REFRESH_COMMAND,
            description=f"post-refresh command env={env_summary}",
            env=post_refresh_env,
        )
    except OSError as exc:
        sys.stderr.write(
            f"[watch] post-refresh command failed to start: {' '.join(POST_REPORT_REFRESH_COMMAND)} ({exc})\n"
        )
        return
    if result.returncode != 0:
        sys.stderr.write(
            f"[watch] post-refresh command exited code={result.returncode}: "
            f"{' '.join(POST_REPORT_REFRESH_COMMAND)}\n"
        )
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)


def analyze_job(job_dir: Path, dry_run: bool) -> None:
    command = [str(ABAG_RBFE), "analyze", job_dir.name, "--batch-dir", str(job_dir.parent.parent), "--execute"]
    if dry_run:
        print("[watch] dry-run analyze:", " ".join(command))
        return
    print("[watch] analyze", job_display_name(job_dir))
    result = run_command(command, check=False)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)


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
        value = _parse_int(tokens[index + 1])
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
    spec_path = job_dir / "job_spec.json"
    if not spec_path.exists():
        return 1
    try:
        payload = json.loads(spec_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 1
    protocol = payload.get("protocol") or {}
    return _parse_mdrun_thread_count(str(protocol.get("mdrun_args") or ""))


def desired_active_mdrun_affinities(
    processes: list[dict[str, Any]],
    *,
    gpu_devices: list[str],
    max_compute_apps_per_gpu: int,
    cpu_count: int,
) -> dict[int, tuple[int, ...]]:
    if cpu_count <= 0 or not gpu_devices:
        return {}
    eligible = [process for process in processes if str(process.get("gpu_device") or "") in gpu_devices]
    if not eligible:
        return {}

    slot_width = max(max(int(process.get("thread_count") or 1), 1) for process in eligible)
    slot_capacity = cpu_count // slot_width
    if slot_capacity <= 0:
        return {}

    grouped: dict[str, list[dict[str, Any]]] = {device: [] for device in gpu_devices}
    for process in eligible:
        grouped[str(process["gpu_device"])].append(process)
    for group in grouped.values():
        group.sort(key=lambda item: (int(item.get("pid") or 0), str(item.get("job_id") or "")))

    observed_max = max((len(group) for group in grouped.values()), default=0)
    slots_per_gpu = max(max_compute_apps_per_gpu, observed_max, 1)

    affinity: dict[int, tuple[int, ...]] = {}
    for gpu_position, gpu_device in enumerate(gpu_devices):
        for local_slot, process in enumerate(grouped.get(gpu_device, [])):
            pid = int(process.get("pid") or 0)
            if pid <= 0:
                continue
            thread_count = max(int(process.get("thread_count") or 1), 1)
            global_slot = (gpu_position * slots_per_gpu + local_slot) % slot_capacity
            start_cpu = global_slot * slot_width
            affinity[pid] = tuple(range(start_cpu, min(start_cpu + thread_count, cpu_count)))
    return affinity


def process_affinity(pid: int) -> tuple[int, ...] | None:
    if pid <= 0:
        return None
    if hasattr(os, "sched_getaffinity"):
        try:
            return tuple(sorted(os.sched_getaffinity(pid)))
        except OSError:
            return None
    taskset = shutil_which("taskset")
    if taskset is None:
        return None
    result = run_command([taskset, "-pc", str(pid)], check=False)
    if result.returncode != 0:
        return None
    match = re.search(r"affinity list:\s*(.+)$", result.stdout.strip())
    if not match:
        return None
    cpu_values: list[int] = []
    for token in match.group(1).split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            start, end = token.split("-", 1)
            if start.isdigit() and end.isdigit():
                cpu_values.extend(range(int(start), int(end) + 1))
        elif token.isdigit():
            cpu_values.append(int(token))
    return tuple(sorted(set(cpu_values))) if cpu_values else None


def set_process_affinity(pid: int, cpu_ids: tuple[int, ...]) -> bool:
    if pid <= 0 or not cpu_ids:
        return False
    if hasattr(os, "sched_setaffinity"):
        try:
            os.sched_setaffinity(pid, set(cpu_ids))
            return True
        except OSError:
            return False
    taskset = shutil_which("taskset")
    if taskset is None:
        return False
    cpu_list = ",".join(str(cpu_id) for cpu_id in cpu_ids)
    result = run_command([taskset, "-pc", cpu_list, str(pid)], check=False)
    return result.returncode == 0


def rebalance_active_mdrun_affinity(
    *,
    gpu_devices: list[str],
    max_compute_apps_per_gpu: int,
    dry_run: bool,
) -> list[dict[str, Any]]:
    cpu_count = os.cpu_count() or 0
    processes = active_mdrun_processes()
    if cpu_count <= 0 or not processes or not gpu_devices:
        return []

    pid_to_gpu = gpu_pid_to_index()
    enriched: list[dict[str, Any]] = []
    for process in processes:
        enriched.append(
            {
                **process,
                "gpu_device": pid_to_gpu.get(int(process.get("pid") or 0), ""),
            }
        )

    desired_affinity = desired_active_mdrun_affinities(
        enriched,
        gpu_devices=gpu_devices,
        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
        cpu_count=cpu_count,
    )
    if not desired_affinity:
        return []

    changes: list[dict[str, Any]] = []
    for process in enriched:
        pid = int(process.get("pid") or 0)
        target = desired_affinity.get(pid)
        if pid <= 0 or target is None:
            continue
        current = process_affinity(pid)
        if current == target:
            continue
        change = {
            "pid": pid,
            "job_id": str(process.get("job_id") or ""),
            "gpu_device": str(process.get("gpu_device") or ""),
            "thread_count": int(process.get("thread_count") or 1),
            "current_cpus": list(current) if current is not None else [],
            "target_cpus": list(target),
        }
        if dry_run or set_process_affinity(pid, target):
            changes.append(change)
    return changes


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
    command = [str(ABAG_RBFE), "resume", job_dir.name, "--batch-dir", str(job_dir.parent.parent), "--execute"]
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
        f"{job_display_name(job_dir).replace('/', '_')}_{int(time.time())}",
        command,
        env=launch_env,
        log_dir=plan_root_for_job_dir(job_dir) / "reports" / "watch",
    )
    print(
        f"[watch] launched {job_display_name(job_dir)} on GPU {gpu_device} "
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


def _parse_int(value: str | None) -> int:
    if value is None or value == "":
        return 0
    try:
        return int(value)
    except ValueError:
        try:
            return int(float(value))
        except ValueError:
            return 0


def _parse_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value or "").strip().lower() in {"1", "true", "yes", "on"}


def canonical_job_id_from_report_key(job_key: str) -> str:
    token = str(job_key or "").strip()
    if not token:
        return ""
    if "/" not in token:
        return token
    return token.rsplit("/", 1)[-1].strip()


def invalid_mutate_output_canonical_job_ids(
    job_rows: dict[str, dict[str, str]],
) -> set[str]:
    blocked: set[str] = set()
    for job_key, row in job_rows.items():
        if not (
            _parse_bool(row.get("current_invalid_mutate_output"))
            or str(row.get("current_invalid_mutate_output_code", "") or "").strip()
        ):
            continue
        job_id = str(row.get("job_id", "") or "").strip() or canonical_job_id_from_report_key(job_key)
        if job_id:
            blocked.add(job_id)
    return blocked


def resumable_priority(
    job_dir: Path,
    states: dict[str, str],
    *,
    queue_positions: dict[str, int],
    job_rows: dict[str, dict[str, str]],
    batch_rows: dict[str, dict[str, str]],
    merged_job_rows: dict[str, dict[str, str]],
    active_canonical_job_ids: set[str] | None = None,
) -> tuple[int, int, int, float, float, int, int, str]:
    try:
        root_priority = all_plan_roots().index(plan_root_for_job_dir(job_dir))
    except ValueError:
        root_priority = len(all_plan_roots())
    job_row = job_rows.get(job_cache_key(job_dir), {})
    merged_job_row = merged_job_rows.get(job_dir.name, {})
    batch_row = batch_rows.get(batch_cache_key(job_dir.parent.parent), {})
    batch_pairs = _parse_int(batch_row.get("paired_job_count"))
    queue_position = queue_positions.get(job_cache_key(job_dir), len(queue_positions) + 10_000)
    sample_completed = _parse_int(job_row.get("sample_completed_windows"))
    sample_total = max(_parse_int(job_row.get("sample_total_windows")), 1)
    equilibrate_completed = _parse_int(job_row.get("equilibrate_completed_repeats"))
    equilibrate_total = max(_parse_int(job_row.get("equilibrate_total_repeats")), 1)
    abs_ddg_error = _parse_float(job_row.get("abs_ddg_error_kcal_mol"))
    if abs_ddg_error is None:
        abs_ddg_error = _parse_float(merged_job_row.get("abs_ddg_error_kcal_mol"))
    error_priority = -(abs_ddg_error if abs_ddg_error is not None else -1.0)
    # When a lane explicitly allows the same canonical job to stay active in another
    # plan root, prefer unique coverage first and only backfill duplicates after that.
    active_elsewhere_priority = (
        1
        if active_canonical_job_ids is not None and job_dir.name in active_canonical_job_ids
        else 0
    )

    if states.get("sample") in {"failed", "blocked_external", "running"}:
        return (
            root_priority,
            active_elsewhere_priority,
            0,
            -(sample_completed / sample_total),
            error_priority,
            batch_pairs,
            queue_position,
            job_cache_key(job_dir),
        )
    if states.get("equilibrate") in {"failed", "blocked_external", "running"}:
        return (
            root_priority,
            active_elsewhere_priority,
            1,
            -(equilibrate_completed / equilibrate_total),
            error_priority,
            batch_pairs,
            queue_position,
            job_cache_key(job_dir),
        )
    if states.get("equilibrate") == "completed":
        return (
            root_priority,
            active_elsewhere_priority,
            2,
            -1.0,
            error_priority,
            batch_pairs,
            queue_position,
            job_cache_key(job_dir),
        )
    if states.get("build_legs") == "completed":
        return (
            root_priority,
            active_elsewhere_priority,
            3,
            0.0,
            error_priority,
            batch_pairs,
            queue_position,
            job_cache_key(job_dir),
        )
    return (
        root_priority,
        active_elsewhere_priority,
        4,
        0.0,
        error_priority,
        batch_pairs,
        queue_position,
        job_cache_key(job_dir),
    )


def job_queue(job_ids: list[str], *, only_listed: bool = False) -> list[Path]:
    queue = job_ids or DEFAULT_QUEUE
    resolved: list[Path] = []
    seen: set[str] = set()
    for job_id in queue:
        for job_dir in resolve_job_dirs(job_id):
            cache_key = job_cache_key(job_dir)
            if cache_key in seen:
                continue
            seen.add(cache_key)
            resolved.append(job_dir)
    if only_listed and job_ids:
        return resolved
    for job_dir in all_job_dirs():
        cache_key = job_cache_key(job_dir)
        if cache_key in seen:
            continue
        seen.add(cache_key)
        resolved.append(job_dir)
    return resolved


def pass_once(
    job_ids: list[str],
    gpu_devices: list[str],
    *,
    only_listed: bool,
    allow_active_elsewhere_job_ids: bool = False,
    max_active_copies_per_job_id: int = 0,
    thread_budget_plan_roots: list[Path] | None = None,
    max_compute_apps_per_gpu: int,
    min_free_gpu_memory_mb: int = 0,
    max_gpu_utilization: int = 0,
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
    active = active_job_ids()
    active_canonical = active_canonical_job_ids()
    active_canonical_counts = active_canonical_job_copy_counts()
    gpu_counts = gpu_compute_counts()
    gpu_stats = (
        gpu_device_stats()
        if max(min_free_gpu_memory_mb, 0) > 0 or max(max_gpu_utilization, 0) > 0
        else {}
    )
    gpu_free_gpus = available_gpu_devices(
        gpu_devices,
        gpu_counts=gpu_counts,
        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
        gpu_stats=gpu_stats,
        min_free_gpu_memory_mb=min_free_gpu_memory_mb,
        max_gpu_utilization=max_gpu_utilization,
    )
    affinity_changes = rebalance_active_mdrun_affinity(
        gpu_devices=gpu_devices,
        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
        dry_run=dry_run,
    )
    cpu_launch_allowed, load_per_core = launch_allowed_by_cpu(max_load_per_core=max_load_per_core)
    thread_launch_allowed, active_threads, active_mdrun_processes = launch_allowed_by_thread_budget(
        max_active_mdrun_threads=max_active_mdrun_threads,
        plan_roots=thread_budget_plan_roots,
    )
    free_gpus = gpu_free_gpus if cpu_launch_allowed and thread_launch_allowed else []
    thread_budget_scope = [str(root) for root in thread_budget_plan_roots or []]
    gpu_status = {}
    for device in gpu_devices:
        entry: dict[str, Any] = {"count": gpu_counts.get(device, 0)}
        if device in gpu_stats:
            entry.update(
                {
                    "memory_used_mb": gpu_stats[device].get("memory_used_mb"),
                    "memory_total_mb": gpu_stats[device].get("memory_total_mb"),
                    "gpu_utilization_percent": gpu_stats[device].get("gpu_utilization_percent"),
                }
            )
        gpu_status[device] = entry
    print(
        f"[watch] active_jobs={sorted(active)} gpu_status={gpu_status} "
        f"max_compute_apps_per_gpu={max_compute_apps_per_gpu} "
        f"min_free_gpu_memory_mb={min_free_gpu_memory_mb} "
        f"max_gpu_utilization={max_gpu_utilization} "
        f"max_load_per_core={max_load_per_core} cpu_load_per_core={load_per_core} "
        f"max_active_mdrun_threads={max_active_mdrun_threads} "
        f"active_mdrun_threads={active_threads} active_mdrun_processes={active_mdrun_processes} "
        f"max_active_copies_per_job_id={max_active_copies_per_job_id} "
        f"thread_budget_plan_roots={thread_budget_scope or ['<all>']} "
        f"free_gpus={free_gpus}"
    )
    if affinity_changes:
        compact = [
            {
                "job_id": change["job_id"],
                "pid": change["pid"],
                "gpu_device": change["gpu_device"],
                "target_cpus": change["target_cpus"],
            }
            for change in affinity_changes
        ]
        print("[watch] rebalanced active mdrun affinity:", json.dumps(compact, ensure_ascii=False))
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
    resumable: list[Path] = []
    blocked: list[str] = []
    completed: list[str] = []
    cooling: list[str] = []
    active_elsewhere: list[str] = []
    active_copy_limited: list[str] = []
    resumable_states: dict[str, dict[str, str]] = {}
    queue = job_queue(job_ids, only_listed=only_listed)
    queue_positions = {job_cache_key(job_dir): index for index, job_dir in enumerate(queue)}
    job_rows, batch_rows = report_priority_data()
    merged_job_rows = merged_priority_job_rows()
    invalid_mutate_output_job_ids = invalid_mutate_output_canonical_job_ids(job_rows)

    for job_dir in queue:
        states = read_stage_states(job_dir)
        job_key = job_cache_key(job_dir)
        job_label = job_display_name(job_dir)
        if states.get("report") == "completed":
            completed.append(job_label)
            continue
        if states.get("prepare") == "blocked_input":
            blocked.append(job_label)
            continue
        if job_dir.name in invalid_mutate_output_job_ids:
            blocked.append(job_label)
            continue
        if job_key in active:
            continue
        if not allow_active_elsewhere_job_ids and job_dir.name in active_canonical:
            active_elsewhere.append(job_label)
            continue
        if (
            allow_active_elsewhere_job_ids
            and max_active_copies_per_job_id > 0
            and active_canonical_counts.get(job_dir.name, 0) >= max_active_copies_per_job_id
        ):
            active_copy_limited.append(job_label)
            continue
        if states.get("sample") == "completed":
            analyzable.append(job_dir)
            continue
        if job_key in recent_launches:
            cooling.append(job_label)
            continue
        resumable.append(job_dir)
        resumable_states[job_key] = states

    resumable.sort(
        key=lambda job_dir: resumable_priority(
            job_dir,
            resumable_states.get(job_cache_key(job_dir), {}),
            queue_positions=queue_positions,
            job_rows=job_rows,
            batch_rows=batch_rows,
            merged_job_rows=merged_job_rows,
            active_canonical_job_ids=active_canonical if allow_active_elsewhere_job_ids else None,
        )
    )

    print(
        "[watch] queue",
        json.dumps(
            {
                "completed": completed,
                "blocked": blocked,
                "analyzable": [job_display_name(job) for job in analyzable],
                "cooling": cooling,
                "active_elsewhere": active_elsewhere,
                "active_copy_limited": active_copy_limited,
                "stale_active": [status.get("job_key") or status["job_id"] for status in stale_statuses],
                "resumable": [job_display_name(job) for job in resumable],
            },
            ensure_ascii=False,
        ),
    )

    changed = False
    for job_dir in analyzable:
        analyze_job(job_dir, dry_run=dry_run)
        changed = True

    launch_pairs: list[tuple[str, Path]] = []
    launch_lock_path: Path | None = None
    if free_gpus and resumable:
        launch_lock_path, launch_lock_message = acquire_launch_coordination_lock()
        if launch_lock_path is None:
            print(f"[watch] launch coordination gate active: {launch_lock_message}; deferring new resumes.")
        else:
            try:
                refreshed_active = active_job_ids()
                refreshed_active_canonical = active_canonical_job_ids()
                refreshed_active_canonical_counts = active_canonical_job_copy_counts()
                refreshed_gpu_counts = gpu_compute_counts()
                refreshed_gpu_stats = (
                    gpu_device_stats()
                    if max(min_free_gpu_memory_mb, 0) > 0 or max(max_gpu_utilization, 0) > 0
                    else {}
                )
                refreshed_cpu_launch_allowed, refreshed_load_per_core = launch_allowed_by_cpu(
                    max_load_per_core=max_load_per_core
                )
                (
                    refreshed_thread_launch_allowed,
                    refreshed_active_threads,
                    _refreshed_active_mdrun_processes,
                ) = launch_allowed_by_thread_budget(
                    max_active_mdrun_threads=max_active_mdrun_threads,
                    plan_roots=thread_budget_plan_roots,
                )
                refreshed_free_gpus = (
                    available_gpu_devices(
                        gpu_devices,
                        gpu_counts=refreshed_gpu_counts,
                        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
                        gpu_stats=refreshed_gpu_stats,
                        min_free_gpu_memory_mb=min_free_gpu_memory_mb,
                        max_gpu_utilization=max_gpu_utilization,
                    )
                    if refreshed_cpu_launch_allowed and refreshed_thread_launch_allowed
                    else []
                )
                coordinated_resumable: list[Path] = []
                for job_dir in resumable:
                    job_key = job_cache_key(job_dir)
                    if job_key in refreshed_active:
                        continue
                    if not allow_active_elsewhere_job_ids and job_dir.name in refreshed_active_canonical:
                        continue
                    if (
                        allow_active_elsewhere_job_ids
                        and max_active_copies_per_job_id > 0
                        and refreshed_active_canonical_counts.get(job_dir.name, 0) >= max_active_copies_per_job_id
                    ):
                        continue
                    if job_key in recent_launches:
                        continue
                    coordinated_resumable.append(job_dir)

                if not refreshed_cpu_launch_allowed and refreshed_free_gpus:
                    print(
                        f"[watch] coordinated CPU launch gate active: load/core={refreshed_load_per_core:.3f} "
                        f">= threshold {max_load_per_core:.3f}; deferring new resumes."
                    )
                if not refreshed_thread_launch_allowed and refreshed_free_gpus:
                    print(
                        f"[watch] coordinated mdrun thread gate active: active_threads={refreshed_active_threads} "
                        f">= threshold {max_active_mdrun_threads}; deferring new resumes."
                    )

                launch_pairs = list(zip(refreshed_free_gpus, coordinated_resumable))
                if max_launches_per_pass > 0:
                    launch_pairs = launch_pairs[:max_launches_per_pass]

                for gpu_device, job_dir in launch_pairs:
                    launch_resume(
                        job_dir,
                        gpu_device,
                        dry_run=dry_run,
                        gpu_devices=gpu_devices,
                        gpu_counts=refreshed_gpu_counts,
                        max_compute_apps_per_gpu=max_compute_apps_per_gpu,
                        mdrun_args_override=mdrun_args_override,
                    )
                    recent_launches[job_cache_key(job_dir)] = current_time
                    changed = True
            finally:
                release_launch_coordination_lock(launch_lock_path)

    # Keep plan reports aligned with live process state while jobs are running.
    if changed or active:
        refresh_reports(dry_run=dry_run)


def main() -> int:
    global RUNS_ROOT, REPORTS_DIR, LOG_DIR, SPLIT_FILE, SPLIT_NAME, PLAN_ROOTS, MERGED_PLAN_ROOT, MERGED_EXTRA_PLAN_ROOTS, POST_REPORT_REFRESH_COMMAND
    args = parse_args()
    roots = [Path(args.plan_root).expanduser().resolve()]
    roots.extend(Path(item).expanduser().resolve() for item in args.extra_plan_root if item)
    PLAN_ROOTS = []
    for root in roots:
        if root not in PLAN_ROOTS:
            PLAN_ROOTS.append(root)
    RUNS_ROOT = PLAN_ROOTS[0]
    REPORTS_DIR = RUNS_ROOT / "reports"
    LOG_DIR = REPORTS_DIR / "watch"
    SPLIT_FILE = Path(args.split_file).expanduser().resolve() if args.split_file else Path()
    SPLIT_NAME = args.split_name.strip()
    MERGED_PLAN_ROOT = Path(args.merged_plan_root).expanduser().resolve() if args.merged_plan_root else None
    merged_extra_roots: list[Path] = []
    for item in args.merged_extra_plan_root:
        if not item:
            continue
        resolved = Path(item).expanduser().resolve()
        if MERGED_PLAN_ROOT is not None and resolved == MERGED_PLAN_ROOT:
            continue
        if resolved not in merged_extra_roots:
            merged_extra_roots.append(resolved)
    MERGED_EXTRA_PLAN_ROOTS = merged_extra_roots
    post_refresh_command = str(getattr(args, "post_refresh_command", "") or "").strip()
    POST_REPORT_REFRESH_COMMAND = [post_refresh_command] if post_refresh_command else []
    if args.wait_for_pid > 0:
        while Path(f"/proc/{args.wait_for_pid}").exists():
            time.sleep(30)

    gpu_devices = normalize_gpu_devices(args.gpu_devices)
    if not gpu_devices:
        print("[watch] warning: no GPU devices detected; resumable jobs will not be launched.")
    only_listed = queue_only_listed(
        args.job_id,
        only_listed=args.only_listed,
        append_rest=args.append_rest,
    )
    thread_budget_plan_roots = _resolve_path_list(getattr(args, "thread_budget_plan_root", []))
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
        return WATCH_SUPERVISOR_NO_RESTART_CODE

    try:
        allow_active_elsewhere_job_ids = bool(getattr(args, "allow_active_elsewhere_job_ids", False))
        while True:
            pass_once(
                args.job_id,
                gpu_devices,
                only_listed=only_listed,
                allow_active_elsewhere_job_ids=allow_active_elsewhere_job_ids,
                max_compute_apps_per_gpu=max(args.max_compute_apps_per_gpu, 0),
                min_free_gpu_memory_mb=max(getattr(args, "min_free_gpu_memory_mb", 0), 0),
                max_gpu_utilization=max(getattr(args, "max_gpu_utilization", 0), 0),
                max_load_per_core=max(args.max_load_per_core, 0.0),
                max_active_mdrun_threads=max(args.max_active_mdrun_threads, 0),
                max_active_copies_per_job_id=max(getattr(args, "max_active_copies_per_job_id", 0), 0),
                thread_budget_plan_roots=thread_budget_plan_roots or None,
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
