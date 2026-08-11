from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path


SCRIPT = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_post_report.sh"
)


def _write_executable(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _fake_root(tmp_path: Path) -> Path:
    root = tmp_path / "abag-rbfep"
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    (plan_root / "reports" / "watch").mkdir(parents=True, exist_ok=True)
    (root / "calls").mkdir(parents=True, exist_ok=True)

    _write_executable(
        root / "benchmarks" / "ab_bind" / "refresh_calibrated_validation_summary.sh",
        """#!/usr/bin/env bash
set -euo pipefail
printf 'calibrated\\n' >> "${ABAG_RBFE_ROOT}/calls/calibrated.log"
""",
    )
    _write_executable(
        root / "benchmarks" / "ab_bind" / "report_3hfm_protocol_regression.py",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

root = pathlib.Path(os.environ["ABAG_RBFE_ROOT"])
payload = {
    "argv": sys.argv[1:],
    "plan_root_env": os.environ.get("PLAN_ROOT", ""),
    "merged_plan_root_env": os.environ.get("MERGED_PLAN_ROOT", ""),
    "merged_extra_plan_roots_env": os.environ.get("MERGED_EXTRA_PLAN_ROOTS", ""),
}
(root / "calls" / "three_hfm.json").write_text(json.dumps(payload), encoding="utf-8")
raise SystemExit(int(os.environ.get("THREE_HFM_EXIT_CODE", "0")))
""",
    )
    _write_executable(
        root / "benchmarks" / "patel_2021_3hfm" / "refresh_patellike_3hfm_summary.sh",
        """#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PY'
import json
import os
from pathlib import Path

root = Path(os.environ["ABAG_RBFE_ROOT"])
payload = {
    "argv": [],
    "batch_dir_env": os.environ.get("BATCH_DIR", ""),
    "summary_output_env": os.environ.get("SUMMARY_OUTPUT", ""),
    "min_interval_seconds_env": os.environ.get("MIN_INTERVAL_SECONDS", ""),
    "fail_on_refresh_error_env": os.environ.get("FAIL_ON_REFRESH_ERROR", ""),
}
(root / "calls" / "patellike_3hfm.json").write_text(json.dumps(payload), encoding="utf-8")
PY
exit "${PATELLIKE_3HFM_EXIT_CODE:-0}"
""",
    )
    _write_executable(
        root / "benchmarks" / "ab_bind" / "report_validation_status.py",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

root = pathlib.Path(os.environ["ABAG_RBFE_ROOT"])
payload = {"argv": sys.argv[1:]}
(root / "calls" / "validation_status.json").write_text(json.dumps(payload), encoding="utf-8")
output_path = pathlib.Path(sys.argv[sys.argv.index("--summary-output") + 1])
output_path.parent.mkdir(parents=True, exist_ok=True)
output_path.write_text("# mock validation status\\n", encoding="utf-8")
""",
    )
    _write_executable(
        root / "benchmarks" / "ab_bind" / "report_project_completion.py",
        """#!/usr/bin/env python3
import json
import os
import pathlib
import sys

root = pathlib.Path(os.environ["ABAG_RBFE_ROOT"])
payload = {"argv": sys.argv[1:]}
(root / "calls" / "project_completion.json").write_text(json.dumps(payload), encoding="utf-8")
summary_output = pathlib.Path(sys.argv[sys.argv.index("--summary-output") + 1])
summary_output.parent.mkdir(parents=True, exist_ok=True)
summary_output.write_text("# mock project completion\\n", encoding="utf-8")
json_output = pathlib.Path(sys.argv[sys.argv.index("--json-output") + 1])
json_output.parent.mkdir(parents=True, exist_ok=True)
json_output.write_text('{"project_complete": false}\\n', encoding="utf-8")
""",
    )
    return root


def _run(root: Path, **env_overrides: str) -> subprocess.CompletedProcess[str]:
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    env = os.environ.copy()
    env.update(
        {
            "ABAG_RBFE_ROOT": str(root),
            "PLAN_ROOT": str(plan_root),
            "MERGED_PLAN_ROOT": str(plan_root),
            "PYTHON_BIN": "python3",
        }
    )
    env.update(env_overrides)
    return subprocess.run(
        ["bash", str(SCRIPT)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )


def test_refresh_validation_post_report_propagates_merged_extra_roots(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    extra_root_1 = root / "runs" / "benchmarks" / "abbind_core_v1_validation_robust_plan"
    extra_root_2 = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_rescues"

    result = _run(
        root,
        MERGED_EXTRA_PLAN_ROOTS=os.pathsep.join(
            [str(extra_root_1), str(plan_root), "", str(extra_root_2)]
        ),
    )

    assert result.returncode == 0
    assert (root / "calls" / "calibrated.log").read_text(encoding="utf-8").splitlines() == [
        "calibrated"
    ]
    payload = json.loads((root / "calls" / "three_hfm.json").read_text(encoding="utf-8"))
    assert payload["argv"] == [
        "--plan-root",
        str(plan_root),
        "--complex-id",
        "3HFM",
        "--extra-plan-root",
        str(extra_root_1),
        "--extra-plan-root",
        str(extra_root_2),
    ]
    patellike_payload = json.loads((root / "calls" / "patellike_3hfm.json").read_text(encoding="utf-8"))
    assert patellike_payload == {
        "argv": [],
        "batch_dir_env": str(
            root / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
        ),
        "summary_output_env": str(
            root
            / "runs"
            / "benchmarks"
            / "patel_2021_3hfm"
            / "patel_2021_3hfm_reference"
            / "reports"
            / "patel_2021_3hfm_summary.json"
        ),
        "min_interval_seconds_env": "0",
        "fail_on_refresh_error_env": "1",
    }
    validation_status_payload = json.loads(
        (root / "calls" / "validation_status.json").read_text(encoding="utf-8")
    )
    assert validation_status_payload["argv"] == [
        "--root",
        str(root),
        "--summary-output",
        str(root / "docs" / "validation_status.md"),
    ]
    assert (root / "docs" / "validation_status.md").read_text(encoding="utf-8") == "# mock validation status\n"
    project_completion_payload = json.loads(
        (root / "calls" / "project_completion.json").read_text(encoding="utf-8")
    )
    assert project_completion_payload["argv"] == [
        "--root",
        str(root),
        "--summary-output",
        str(root / "docs" / "project_completion_status.md"),
        "--json-output",
        str(root / "runs" / "benchmarks" / "project_completion_summary.json"),
    ]
    assert (root / "docs" / "project_completion_status.md").read_text(encoding="utf-8") == "# mock project completion\n"
    log_text = (
        plan_root / "reports" / "watch" / "validation_post_report_refresh.log"
    ).read_text(encoding="utf-8")
    assert "3hfm_protocol_regression success" in log_text
    assert "patellike_3hfm success" in log_text
    assert "validation_status success" in log_text
    assert "project_completion success" in log_text


def test_refresh_validation_post_report_tolerates_insufficient_pairs_exit_code(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"

    result = _run(root, THREE_HFM_EXIT_CODE="2")

    assert result.returncode == 0
    log_text = (
        plan_root / "reports" / "watch" / "validation_post_report_refresh.log"
    ).read_text(encoding="utf-8")
    assert "3hfm_protocol_regression failed rc=2" in log_text


def test_refresh_validation_post_report_tolerates_patellike_insufficient_pairs_exit_code(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"

    result = _run(root, PATELLIKE_3HFM_EXIT_CODE="2")

    assert result.returncode == 0
    log_text = (
        plan_root / "reports" / "watch" / "validation_post_report_refresh.log"
    ).read_text(encoding="utf-8")
    assert "patellike_3hfm failed rc=2" in log_text


def test_refresh_validation_post_report_propagates_extra_roots_to_calibrated_refresh(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    extra_root_1 = root / "runs" / "benchmarks" / "abbind_core_v1_validation_sampling_qc_rescues_verify_preeq_20260611_1"
    extra_root_2 = root / "runs" / "benchmarks" / "abbind_core_v1_validation_deep_rescues"
    _write_executable(
        root / "benchmarks" / "ab_bind" / "refresh_calibrated_validation_summary.sh",
        """#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PY'
import json
import os
from pathlib import Path

root = Path(os.environ["ABAG_RBFE_ROOT"])
payload = {
    "plan_root_env": os.environ.get("PLAN_ROOT", ""),
    "merged_plan_root_env": os.environ.get("MERGED_PLAN_ROOT", ""),
    "merged_extra_plan_roots_env": os.environ.get("MERGED_EXTRA_PLAN_ROOTS", ""),
}
(root / "calls" / "calibrated_env.json").write_text(json.dumps(payload), encoding="utf-8")
PY
""",
    )

    result = _run(
        root,
        MERGED_EXTRA_PLAN_ROOTS=os.pathsep.join(
            [str(extra_root_1), str(plan_root), "", str(extra_root_2)]
        ),
    )

    assert result.returncode == 0
    payload = json.loads((root / "calls" / "calibrated_env.json").read_text(encoding="utf-8"))
    assert payload == {
        "plan_root_env": str(plan_root),
        "merged_plan_root_env": str(plan_root),
        "merged_extra_plan_roots_env": os.pathsep.join(
            [str(extra_root_1), str(plan_root), "", str(extra_root_2)]
        ),
    }


def test_refresh_validation_post_report_can_fail_hard_when_requested(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, THREE_HFM_EXIT_CODE="7", FAIL_ON_REFRESH_ERROR="1")

    assert result.returncode == 7
