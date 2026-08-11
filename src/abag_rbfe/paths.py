"""Project path discovery helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ProjectPaths:
    repo_root: Path
    runs_root: Path
    benchmarks_root: Path
    vendor_pmx_root: Path

    @classmethod
    def discover(cls) -> "ProjectPaths":
        repo_root = Path(__file__).resolve().parents[2]
        return cls(
            repo_root=repo_root,
            runs_root=repo_root / "runs",
            benchmarks_root=repo_root / "benchmarks" / "ab_bind",
            vendor_pmx_root=repo_root / "vendor" / "pmx",
        )


def resolve_job_dir(identifier: str, batch_dir: str | None) -> Path:
    candidate = Path(identifier)
    if candidate.is_dir():
        return candidate.resolve()
    if batch_dir is None:
        raise ValueError("Provide --batch-dir when IDENTIFIER is not a direct job path")
    job_dir = Path(batch_dir).resolve() / "jobs" / identifier
    if not job_dir.is_dir():
        raise FileNotFoundError(f"Job directory not found: {job_dir}")
    return job_dir

