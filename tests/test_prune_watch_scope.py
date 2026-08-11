from __future__ import annotations

import importlib.util
from pathlib import Path
import signal


def _load_module():
    module_path = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/prune_watch_scope.py")
    spec = importlib.util.spec_from_file_location("prune_watch_scope", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_outside_scope_processes_filters_to_plan_root_and_allow_list(tmp_path: Path) -> None:
    module = _load_module()
    plan_root = tmp_path / "robust_plan"
    ps_output = "\n".join(
        [
            (
                "liuchao  1111     1 99 00:00 ? 00:00:00 "
                f"/path/to/gmx mdrun -s {plan_root}/abbind_1bj1_core_v1/jobs/1bj1-antigen-w-g92a/legs/complex/rep01/lambda_000/topol.tpr"
            ),
            (
                "liuchao  2222     1  0 00:00 ? 00:00:00 "
                f"{plan_root}/.venv/bin/python {tmp_path}/.venv/bin/abag-rbfe resume 3hfm-antibody-h-y50a "
                f"--batch-dir {plan_root}/abbind_3hfm_core_v1 --execute"
            ),
            (
                "liuchao  3333     1  0 00:00 ? 00:00:00 bash "
                f"{plan_root}/abbind_1bj1_core_v1/jobs/1bj1-antigen-w-g92a/artifacts/commands/equilibrate.sh"
            ),
            (
                "liuchao  4444     1  0 00:00 ? 00:00:00 "
                f"/other/root/abbind_1bj1_core_v1/jobs/1bj1-antigen-w-g92a/artifacts/commands/equilibrate.sh"
            ),
        ]
    )

    outside = module.outside_scope_processes(
        plan_root,
        allowed_job_ids={"3hfm-antibody-h-y50a"},
        ps_output=ps_output,
    )

    assert outside == [
        {
            "pid": 1111,
            "job_id": "1bj1-antigen-w-g92a",
            "kind": "mdrun",
            "command": (
                f"/path/to/gmx mdrun -s {plan_root}/abbind_1bj1_core_v1/jobs/1bj1-antigen-w-g92a/"
                "legs/complex/rep01/lambda_000/topol.tpr"
            ),
        },
        {
            "pid": 3333,
            "job_id": "1bj1-antigen-w-g92a",
            "kind": "stage_script",
            "command": f"bash {plan_root}/abbind_1bj1_core_v1/jobs/1bj1-antigen-w-g92a/artifacts/commands/equilibrate.sh",
        },
    ]


def test_terminate_processes_sends_requested_signal(monkeypatch) -> None:
    module = _load_module()
    seen: list[tuple[int, signal.Signals]] = []
    monkeypatch.setattr(module.os, "kill", lambda pid, sig: seen.append((pid, sig)))

    terminated, errors = module.terminate_processes(
        [
            {"pid": 1111, "job_id": "job-a", "kind": "mdrun", "command": "gmx mdrun ..."},
            {"pid": 2222, "job_id": "job-b", "kind": "controller", "command": "abag-rbfe resume ..."},
        ],
        sig=signal.SIGTERM,
    )

    assert seen == [(1111, signal.SIGTERM), (2222, signal.SIGTERM)]
    assert terminated == [
        {"pid": 1111, "job_id": "job-a", "kind": "mdrun", "command": "gmx mdrun ..."},
        {"pid": 2222, "job_id": "job-b", "kind": "controller", "command": "abag-rbfe resume ..."},
    ]
    assert errors == []
