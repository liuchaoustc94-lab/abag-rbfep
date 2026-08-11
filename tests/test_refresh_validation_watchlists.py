from __future__ import annotations

from datetime import datetime
import importlib.util
import json
import os
import subprocess
from pathlib import Path


SCRIPT = Path(
    "/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_validation_watchlists.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location("refresh_validation_watchlists", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_root(tmp_path: Path) -> Path:
    root = tmp_path / "abag-rbfep"
    summary_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
    )
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "generated_at": "2026-06-10T09:23:29Z",
        "active_alternate_ready_hotspots": [
            {
                "job_id": "1cz8-antigen-w-g92a",
                "complex_id": "1CZ8",
                "active_alternate_source_plan_roots": ",".join(
                    [
                        str(
                            root
                            / "runs"
                            / "benchmarks"
                            / "abbind_core_v1_validation_priority_rescues"
                        ),
                    ]
                ),
                "abs_ddg_error_kcal_mol": 5.59,
            },
            {
                "job_id": "3hfm-antibody-h-y33a",
                "complex_id": "3HFM",
                "active_alternate_source_plan_roots": ",".join(
                    [
                        str(
                            root
                            / "runs"
                            / "benchmarks"
                            / "abbind_core_v1_validation_priority_rescues"
                        ),
                        str(
                            root
                            / "runs"
                            / "benchmarks"
                            / "abbind_core_v1_validation_robust_plan"
                        ),
                    ]
                ),
                "abs_ddg_error_kcal_mol": 16.61,
            },
            {
                "job_id": "3hfm-antigen-y-y20a",
                "complex_id": "3HFM",
                "active_alternate_source_plan_roots": str(
                    root
                    / "runs"
                    / "benchmarks"
                    / "abbind_core_v1_validation_robust_plan"
                ),
                "abs_ddg_error_kcal_mol": 6.45,
            },
        ]
    }
    (summary_dir / "plan_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    (summary_dir / "benchmark_pairs.csv").write_text(
        "\n".join(
            [
                "complex_id,job_id,predicted_ddg_kcal_mol,experimental_ddg_kcal_mol",
                "1BJ1,1bj1-antigen-w-g92a,3.7,3.7",
                "1CZ8,1cz8-antigen-w-g92a,4.1,4.1",
                "3HFM,3hfm-antibody-h-y33a,-10.6,6.0",
                "3HFM,3hfm-antigen-y-y20a,-1.6,4.9",
                "3NPS,3nps-antigen-a-h138a,1.9,1.9",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    priority_reports_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
    )
    priority_reports_dir.mkdir(parents=True, exist_ok=True)
    (priority_reports_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,job_id,latest_stage,latest_stage_state,resumable,"
                    "ddg_ready,experimental_ddg_kcal_mol,qc_status,diagnostic_family,"
                    "diagnostic_code,complex_repeat_spread_kcal_mol,apo_repeat_spread_kcal_mol,"
                    "repeat_spread_legs,primary_repeat_spread_leg,ddg_bar_stderr_kcal_mol,"
                    "max_bar_stderr_kcal_mol,current_invalid_mutate_output,current_invalid_mutate_output_code,"
                    "benchmark_qc_qualified,abs_ddg_error_kcal_mol"
                ),
                "1MLC,1mlc-antibody-h-s57a,equilibrate,stale_running,True,False,-0.38,,,,,,,,,",
                "1MLC,1mlc-antibody-l-n92a,build_legs,completed,True,False,-1.25,,,,,,,,,",
                "1BJ1,1bj1-antigen-w-i80a,sample,stale_running,True,False,0.82,,,,,,,,,",
                "1BJ1,1bj1-antigen-v-f17a,equilibrate,stale_running,True,False,0.0,,,,,,,,,",
                "2NZ9,2nz9-antigen-a-f953a,build_legs,completed,True,False,3.34,,,,,,,,,",
                "2NZ9,2nz9-antigen-a-t1063a,build_legs,completed,True,False,2.37,,,,,,,,,",
                "2NZ9,2nz9-antigen-a-h1064a,mutate,blocked_input,False,False,7.42,,,,,,,,,",
                "3HFM,3hfm-antibody-l-y50a,build_legs,completed,True,False,4.57,,,,,,,,,",
                "3HFM,3hfm-antibody-h-y33a,report,completed,False,True,6.0,warning,qc,qc_repeat_spread,3.2,2.1,\"complex,apo\",complex,2.5,10.0,,,False,16.61",
                "3HFM,3hfm-antigen-y-y20a,report,completed,False,True,4.9,warning,qc,qc_overlap,0.8,2.6,apo,apo,1.8,10.0,,,False,6.45",
                "1CZ8,1cz8-antigen-w-g92a,report,completed,False,True,4.1,warning,qc,qc_bar_stderr,1.4,1.3,complex,complex,5.6,10.0,,,False,5.59",
                "1CZ8,1cz8-antigen-w-h90a,report,completed,False,True,1.0,pass,qc,qc_pass,0.4,0.55,\"complex,apo\",complex,0.6,10.0,,,True,14.22",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    robust_reports_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
    )
    robust_reports_dir.mkdir(parents=True, exist_ok=True)
    (robust_reports_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,job_id,latest_stage,latest_stage_state,resumable,"
                    "ddg_ready,experimental_ddg_kcal_mol,qc_status,diagnostic_family,"
                    "diagnostic_code,complex_repeat_spread_kcal_mol,apo_repeat_spread_kcal_mol,"
                    "repeat_spread_legs,primary_repeat_spread_leg,ddg_bar_stderr_kcal_mol,"
                    "max_bar_stderr_kcal_mol,current_invalid_mutate_output,current_invalid_mutate_output_code"
                ),
                "1BJ1,1bj1-antigen-w-g88a,equilibrate,stale_running,True,False,2.76,,,,,,,,,",
                "1BJ1,1bj1-antigen-w-h86a,equilibrate,stale_running,True,False,0.0,,,,,,,,,",
                "1CZ8,1cz8-antigen-w-h90a,,not_started,True,False,1.0,,,,,,,,,",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return root


def _run(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["python3", str(SCRIPT), "--root", str(root), *args],
        check=False,
        capture_output=True,
        text=True,
    )


def test_choose_summary_path_prefers_richer_validation_bundle_over_newer_priority_bundle(tmp_path: Path) -> None:
    module = _load_module()
    priority = tmp_path / "priority" / "plan_summary.json"
    quick = tmp_path / "quick" / "plan_summary.json"
    priority.parent.mkdir(parents=True, exist_ok=True)
    quick.parent.mkdir(parents=True, exist_ok=True)

    priority.write_text(
        json.dumps(
            {
                "selection": {"split_name": "validation"},
                "source_plan_roots": ["a", "b", "c", "d", "e", "f"],
                "benchmark_metrics": {"paired_job_count": 41},
            }
        ),
        encoding="utf-8",
    )
    quick.write_text(
        json.dumps(
            {
                "selection": {"split_name": "validation"},
                "source_plan_roots": ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"],
                "benchmark_metrics": {"paired_job_count": 45},
            }
        ),
        encoding="utf-8",
    )
    os.utime(quick, (1_700_000_000, 1_700_000_000))
    os.utime(priority, (1_800_000_000, 1_800_000_000))

    chosen = module.choose_summary_path(
        [priority, quick],
        preferred_path=priority,
        preferred_split_name="validation",
    )

    assert chosen == quick


def test_refresh_validation_watchlists_selects_robust_and_rescue_jobs(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    robust = _run(
        root,
        "--mode",
        "robust",
        "--proactive-robust-job-id",
        "1mlc-antibody-h-s57a",
        "--proactive-robust-job-id",
        "1mlc-antibody-h-s57v",
    )
    rescue = _run(root, "--mode", "rescue")

    assert robust.returncode == 0
    assert robust.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1mlc-antibody-h-s57a",
        "1mlc-antibody-h-s57v",
        "1mlc-antibody-h-t31a",
        "1mlc-antibody-h-t31v",
        "1mlc-antibody-l-n92a",
    ]
    assert rescue.returncode == 0
    assert rescue.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "1cz8-antigen-w-g92a",
    ]


def test_refresh_validation_watchlists_robust_mode_can_append_pass_qc_outlier_fallback(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    robust = _run(
        root,
        "--mode",
        "robust",
        "--robust-pass-outlier-threshold",
        "5.0",
        "--output-json",
        "runs/benchmarks/abbind_core_v1_validation_robust_plan/reports/watch/validation_watchlist_refresh.json",
    )

    assert robust.returncode == 0
    assert robust.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-h90a",
        "1mlc-antibody-h-s57a",
        "1mlc-antibody-h-s57v",
        "1mlc-antibody-h-t31a",
        "1mlc-antibody-h-t31v",
        "1mlc-antibody-l-n92a",
    ]
    output_json = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
        / "watch"
        / "validation_watchlist_refresh.json"
    )
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    assert payload["robust_pass_outlier_threshold"] == 5.0
    assert payload["robust_pass_outlier_job_ids"] == ["1cz8-antigen-w-h90a"]


def test_refresh_validation_watchlists_writes_summary_report(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    output_json = root / "runs" / "benchmarks" / "abbind_core_v1_validation_robust_plan" / "reports" / "watch" / "validation_watchlist_refresh.json"
    summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["selection"] = {"split_name": "validation"}
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = _run(root, "--mode", "all", "--output-json", str(output_json))

    assert result.returncode == 0
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    assert payload["watchlist_generated_at"].endswith("Z")
    assert datetime.fromisoformat(payload["watchlist_generated_at"].replace("Z", "+00:00"))
    assert payload["selected_summary_generated_at"] == "2026-06-10T09:23:29Z"
    assert payload["selected_summary_selection_split_name"] == "validation"
    assert payload["selected_summary_mtime_utc"].endswith("Z")
    assert payload["active_alternate_ready_hotspot_count"] == 3
    assert payload["selected_benchmark_pairs_path"].endswith("benchmark_pairs.csv")
    assert payload["robust_pass_outlier_threshold"] == 0.0
    assert payload["robust_pass_outlier_job_ids"] == []
    assert payload["robust_job_ids"][:2] == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
    ]
    assert payload["rescue_job_ids"] == [
        "3hfm-antibody-h-y33a",
        "1cz8-antigen-w-g92a",
    ]
    assert payload["targeted_job_ids"] == [
        "3hfm-antibody-h-y33a",
    ]
    assert payload["sampling_qc_job_ids"] == []
    assert payload["sampling_qc_excluded_targeted_primary_repeat_spread_leg_job_ids"] == [
        "3hfm-antibody-h-y33a",
    ]
    assert payload["stale_job_ids"] == [
        "1mlc-antibody-h-s57a",
        "1bj1-antigen-w-i80a",
        "1bj1-antigen-v-f17a",
    ]
    assert payload["gap_job_ids"] == [
        "1mlc-antibody-h-s57a",
        "2nz9-antigen-a-f953a",
        "2nz9-antigen-a-t1063a",
        "1mlc-antibody-l-n92a",
    ]
    assert payload["backlog_job_ids"] == [
        "3hfm-antibody-l-y50a",
    ]
    assert payload["gap_complex_ids"] == [
        "1MLC",
        "2NZ9",
    ]
    assert payload["ultra_pearson_gain_threshold"] == 0.2
    assert payload["ultra_abs_error_threshold"] == 5.0
    assert payload["ultra_pass_outlier_threshold"] == 0.0
    assert payload["ultra_job_ids"] == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
    ]
    assert payload["ultra_pass_outlier_job_ids"] == []
    assert payload["queue_excluded_job_ids"] == [
        "2nz9-antigen-a-h1064a",
    ]
    assert payload["hotspot_complex_impact_pearson_gain"]["3hfm-antibody-h-y33a"] is not None


def test_refresh_validation_watchlists_writes_hotspot_sampling_qc_payload(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    output_json = root / "runs" / "benchmarks" / "abbind_core_v1_validation_robust_plan" / "reports" / "watch" / "validation_watchlist_refresh.json"

    result = _run(root, "--mode", "all", "--output-json", str(output_json))

    assert result.returncode == 0
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    y33a = payload["hotspot_sampling_qc_by_job_id"]["3hfm-antibody-h-y33a"]
    assert y33a == {
        "qc_status": "warning",
        "diagnostic_family": "qc",
        "diagnostic_code": "qc_repeat_spread",
        "primary_repeat_spread_leg": "complex",
        "repeat_spread_legs": ["complex", "apo"],
        "complex_repeat_spread_kcal_mol": 3.2,
        "apo_repeat_spread_kcal_mol": 2.1,
        "ddg_repeat_range_kcal_mol": None,
        "ddg_bar_stderr_kcal_mol": 2.5,
        "max_bar_stderr_kcal_mol": 10.0,
        "prefer_targeted_primary_repeat_spread_leg": True,
    }
    y20a = payload["hotspot_sampling_qc_by_job_id"]["3hfm-antigen-y-y20a"]
    assert y20a["primary_repeat_spread_leg"] == "apo"
    assert y20a["repeat_spread_legs"] == ["apo"]
    assert y20a["diagnostic_code"] == "qc_overlap"
    assert y20a["prefer_targeted_primary_repeat_spread_leg"] is False


def test_refresh_validation_watchlists_emits_hotspots_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "hotspots")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
    ]


def test_refresh_validation_watchlists_emits_sampling_qc_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "sampling_qc")

    assert result.returncode == 0
    assert result.stdout.splitlines() == []


def test_refresh_validation_watchlists_sampling_qc_mode_can_append_no_active_alternate_outlier_fallback(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    priority_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    priority_plan_jobs.write_text(
        priority_plan_jobs.read_text(encoding="utf-8")
        + (
            "1MLC,1mlc-antibody-l-n32g,report,completed,False,True,-1.25,warning,qc,qc_repeat_spread,"
            "1.324,1.164,\"complex,apo\",complex,2.2,10.0,,,False,3.24\n"
        ),
        encoding="utf-8",
    )
    output_json = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_sampling_qc_rescues"
        / "reports"
        / "watch"
        / "validation_watchlist_refresh.json"
    )

    result = _run(
        root,
        "--mode",
        "sampling_qc",
        "--sampling-qc-no-active-alt-abs-error-threshold",
        "3.0",
        "--output-json",
        str(output_json),
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "1mlc-antibody-l-n32g",
    ]
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    assert payload["sampling_qc_job_ids"] == [
        "1mlc-antibody-l-n32g",
    ]
    assert payload["sampling_qc_no_active_alt_abs_error_threshold"] == 3.0
    assert payload["sampling_qc_no_active_alt_outlier_job_ids"] == [
        "1mlc-antibody-l-n32g",
    ]


def test_refresh_validation_watchlists_emits_targeted_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "targeted")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
    ]


def test_refresh_validation_watchlists_targeted_mode_can_append_no_active_alternate_outlier_fallback(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    priority_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    priority_plan_jobs.write_text(
        priority_plan_jobs.read_text(encoding="utf-8")
        + (
            "3NPS,3nps-antigen-a-f50a,report,completed,False,True,1.25,warning,qc,qc_repeat_spread,"
            "2.76,1.44,\"complex,apo\",complex,2.0,10.0,,,False,0.70\n"
        ),
        encoding="utf-8",
    )
    output_json = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_targeted_repeat_spread_rescues"
        / "reports"
        / "watch"
        / "validation_watchlist_refresh.json"
    )

    result = _run(
        root,
        "--mode",
        "targeted",
        "--targeted-no-active-alt-abs-error-threshold",
        "0.5",
        "--output-json",
        str(output_json),
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3nps-antigen-a-f50a",
    ]
    payload = json.loads(output_json.read_text(encoding="utf-8"))
    assert payload["targeted_job_ids"] == [
        "3hfm-antibody-h-y33a",
        "3nps-antigen-a-f50a",
    ]
    assert payload["targeted_no_active_alt_abs_error_threshold"] == 0.5
    assert payload["targeted_no_active_alt_outlier_job_ids"] == [
        "3nps-antigen-a-f50a",
    ]


def test_sampling_qc_candidates_excludes_dominant_primary_leg_repeat_spread_rows() -> None:
    module = _load_module()

    dominant = {
        "job_id": "3hfm-antibody-h-y33a",
        "qc_status": "warning",
        "diagnostic_code": "qc_repeat_spread",
        "repeat_spread_legs": "complex,apo",
        "primary_repeat_spread_leg": "complex",
        "complex_repeat_spread_kcal_mol": 3.2,
        "apo_repeat_spread_kcal_mol": 2.1,
    }
    balanced = {
        "job_id": "1bj1-antigen-v-y21a",
        "qc_status": "warning",
        "diagnostic_code": "qc_repeat_spread",
        "repeat_spread_legs": "complex,apo",
        "primary_repeat_spread_leg": "complex",
        "complex_repeat_spread_kcal_mol": 2.9,
        "apo_repeat_spread_kcal_mol": 2.4,
    }

    selected = module.sampling_qc_candidates([dominant, balanced], complex_ids=set())

    assert module._prefer_targeted_primary_repeat_spread_leg(dominant) is True
    assert module._prefer_targeted_primary_repeat_spread_leg(balanced) is False
    assert [row["job_id"] for row in selected] == ["1bj1-antigen-v-y21a"]


def test_targeted_primary_repeat_spread_candidates_include_single_leg_and_dominant_rows() -> None:
    module = _load_module()

    dominant = {
        "job_id": "3hfm-antibody-h-y33a",
        "qc_status": "warning",
        "diagnostic_code": "qc_repeat_spread",
        "repeat_spread_legs": "complex,apo",
        "primary_repeat_spread_leg": "complex",
        "complex_repeat_spread_kcal_mol": 3.2,
        "apo_repeat_spread_kcal_mol": 2.1,
    }
    single_leg = {
        "job_id": "1vfb-antibody-h-y32a",
        "qc_status": "warning",
        "diagnostic_code": "qc_repeat_spread",
        "repeat_spread_legs": "complex",
        "primary_repeat_spread_leg": "complex",
    }
    balanced = {
        "job_id": "1bj1-antigen-v-y21a",
        "qc_status": "warning",
        "diagnostic_code": "qc_repeat_spread",
        "repeat_spread_legs": "complex,apo",
        "primary_repeat_spread_leg": "complex",
        "complex_repeat_spread_kcal_mol": 2.9,
        "apo_repeat_spread_kcal_mol": 2.4,
    }

    selected = module.targeted_primary_repeat_spread_candidates([dominant, single_leg, balanced])

    assert [row["job_id"] for row in selected] == [
        "3hfm-antibody-h-y33a",
        "1vfb-antibody-h-y32a",
    ]


def test_targeted_primary_repeat_spread_candidates_prioritize_primary_leg_repeat_spread() -> None:
    module = _load_module()

    rows = [
        {
            "job_id": "dominant-low",
            "qc_status": "warning",
            "diagnostic_code": "qc_repeat_spread",
            "repeat_spread_legs": "complex,apo",
            "primary_repeat_spread_leg": "complex",
            "complex_repeat_spread_kcal_mol": 3.2,
            "apo_repeat_spread_kcal_mol": 2.0,
            "ddg_repeat_range_kcal_mol": 4.2,
        },
        {
            "job_id": "single-mid",
            "qc_status": "warning",
            "diagnostic_code": "qc_repeat_spread",
            "repeat_spread_legs": "complex",
            "primary_repeat_spread_leg": "complex",
            "complex_repeat_spread_kcal_mol": 4.1,
            "ddg_repeat_range_kcal_mol": 4.4,
        },
        {
            "job_id": "dominant-high",
            "qc_status": "warning",
            "diagnostic_code": "qc_repeat_spread",
            "repeat_spread_legs": "complex,apo",
            "primary_repeat_spread_leg": "apo",
            "complex_repeat_spread_kcal_mol": 2.0,
            "apo_repeat_spread_kcal_mol": 4.6,
            "ddg_repeat_range_kcal_mol": 5.1,
        },
    ]

    selected = module.targeted_primary_repeat_spread_candidates(rows)

    assert [row["job_id"] for row in selected] == [
        "dominant-high",
        "single-mid",
        "dominant-low",
    ]


def test_sampling_qc_candidates_prioritize_broad_two_leg_repeat_spread() -> None:
    module = _load_module()

    rows = [
        {
            "job_id": "broad-mid",
            "complex_id": "3HFM",
            "qc_status": "warning",
            "diagnostic_code": "qc_repeat_spread",
            "repeat_spread_legs": "complex,apo",
            "primary_repeat_spread_leg": "complex",
            "complex_repeat_spread_kcal_mol": 4.3,
            "apo_repeat_spread_kcal_mol": 3.1,
            "ddg_repeat_range_kcal_mol": 5.0,
        },
        {
            "job_id": "broad-strong",
            "complex_id": "3HFM",
            "qc_status": "warning",
            "diagnostic_code": "qc_repeat_spread",
            "repeat_spread_legs": "complex,apo",
            "primary_repeat_spread_leg": "complex",
            "complex_repeat_spread_kcal_mol": 4.2,
            "apo_repeat_spread_kcal_mol": 3.9,
            "ddg_repeat_range_kcal_mol": 5.8,
        },
        {
            "job_id": "broad-weak",
            "complex_id": "1MLC",
            "qc_status": "warning",
            "diagnostic_code": "qc_repeat_spread",
            "repeat_spread_legs": "complex,apo",
            "primary_repeat_spread_leg": "complex",
            "complex_repeat_spread_kcal_mol": 3.6,
            "apo_repeat_spread_kcal_mol": 2.9,
            "ddg_repeat_range_kcal_mol": 4.6,
        },
    ]

    selected = module.sampling_qc_candidates(rows, complex_ids=set())

    assert [row["job_id"] for row in selected] == [
        "broad-strong",
        "broad-mid",
        "broad-weak",
    ]


def test_sampling_qc_candidates_keeps_rows_when_qc_report_has_overlap_on_nonprimary_leg(tmp_path: Path) -> None:
    module = _load_module()
    root = tmp_path / "abag-rbfep"
    plan_root = root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan"
    qc_report_path = (
        plan_root
        / "abbind_3hfm_core_v1"
        / "jobs"
        / "3hfm-antibody-l-n32a"
        / "results"
        / "qc_report.json"
    )
    qc_report_path.parent.mkdir(parents=True, exist_ok=True)
    qc_report_path.write_text(
        json.dumps(
            {
                "job_id": "3hfm-antibody-l-n32a",
                "ddg_repeat_range_kcal_mol": 4.996030226648912,
                "max_repeat_delta_kcal_mol": 1.0,
                "ddg_bar_stderr_kcal_mol": 0.8575962466989301,
                "max_bar_stderr_kcal_mol": 10.0,
                "repeat_spread_legs": ["complex", "apo"],
                "primary_repeat_spread_leg": "complex",
                "legs": {
                    "complex": {"repeat_delta_kcal_mol_range": 7.55256850539034},
                    "apo": {"repeat_delta_kcal_mol_range": 2.556538278741428},
                },
                "overlap_threshold": 0.2,
                "overlap_assessment": {
                    "legs": {
                        "complex": {"overlap_score_min": 0.20637600241960083},
                        "apo": {"overlap_score_min": 0.17190950889142617},
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    row = {
        "source_plan_root": str(plan_root),
        "batch_id": "abbind_3hfm_core_v1",
        "job_id": "3hfm-antibody-l-n32a",
        "qc_status": "warning",
        "diagnostic_code": "qc_repeat_spread",
        "repeat_spread_legs": "complex,apo",
        "primary_repeat_spread_leg": "complex",
        "complex_repeat_spread_kcal_mol": 7.55256850539034,
        "apo_repeat_spread_kcal_mol": 2.556538278741428,
    }

    assert module._prefer_targeted_primary_repeat_spread_leg(row) is False
    assert [item["job_id"] for item in module.sampling_qc_candidates([row], complex_ids=set())] == [
        "3hfm-antibody-l-n32a"
    ]


def test_refresh_validation_watchlists_prefers_taxonomy_hotspots_when_present(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["validation_failure_taxonomy"] = {
        "hotspots": [
            {
                "category": "qc_sampling_issue",
                "complex_id": "3HFM",
                "job_id": "3hfm-antibody-h-y33a",
                "diagnostic_code": "qc_repeat_spread",
                "abs_ddg_error_kcal_mol": 16.61,
                "ddg_ready": True,
            },
            {
                "category": "qc_sampling_issue",
                "complex_id": "3HFM",
                "job_id": "3hfm-antibody-h-c95a",
                "diagnostic_code": "qc_repeat_spread",
                "abs_ddg_error_kcal_mol": 7.91,
                "ddg_ready": True,
            },
        ]
    }
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = _run(root, "--mode", "hotspots")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antibody-h-c95a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
    ]


def test_refresh_validation_watchlists_emits_ultra_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "ultra")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
    ]


def test_refresh_validation_watchlists_ultra_mode_requires_active_alternate_for_abs_error_outliers(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["validation_failure_taxonomy"] = {
        "hotspots": [
            {
                "category": "qc_sampling_issue",
                "complex_id": "1BJ1",
                "job_id": "1bj1-antigen-w-g92a",
                "diagnostic_code": "qc_repeat_spread",
                "abs_ddg_error_kcal_mol": 7.98,
                "ddg_ready": True,
            }
        ]
    }
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = _run(root, "--mode", "ultra")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
    ]


def test_refresh_validation_watchlists_ultra_mode_can_append_pass_qc_outlier_fallback(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    priority_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    lines = priority_plan_jobs.read_text(encoding="utf-8").splitlines()
    lines[0] += ",alternate_candidate_count"
    updated = [lines[0]]
    for line in lines[1:]:
        if line.startswith("1CZ8,1cz8-antigen-w-h90a,"):
            updated.append(f"{line},1")
        else:
            updated.append(f"{line},")
    priority_plan_jobs.write_text("\n".join(updated) + "\n", encoding="utf-8")

    result = _run(
        root,
        "--mode",
        "ultra",
        "--ultra-pass-outlier-threshold",
        "10.0",
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
        "1cz8-antigen-w-h90a",
    ]


def test_refresh_validation_watchlists_ultra_mode_uses_current_active_alternates_from_robust_root(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["validation_failure_taxonomy"] = {
        "hotspots": [
            {
                "category": "qc_sampling_issue",
                "complex_id": "1BJ1",
                "job_id": "1bj1-antigen-w-g88a",
                "diagnostic_code": "qc_repeat_spread",
                "abs_ddg_error_kcal_mol": 5.50,
                "ddg_ready": True,
            }
        ]
    }
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = _run(root, "--mode", "ultra")

    assert result.returncode == 0
    lines = result.stdout.splitlines()
    assert set(lines) == {
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1cz8-antigen-w-g92a",
        "1bj1-antigen-w-g88a",
    }


def test_refresh_validation_watchlists_ultra_mode_uses_current_active_alternates_from_deep_root(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    robust_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    robust_lines = robust_plan_jobs.read_text(encoding="utf-8").splitlines()
    robust_plan_jobs.write_text(
        "\n".join(
            [
                robust_lines[0],
                robust_lines[2],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    deep_reports_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_deep_rescues"
        / "reports"
    )
    deep_reports_dir.mkdir(parents=True, exist_ok=True)
    (deep_reports_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,job_id,latest_stage,latest_stage_state,resumable,"
                    "ddg_ready,experimental_ddg_kcal_mol"
                ),
                "1BJ1,1bj1-antigen-w-g88a,equilibrate,running,True,False,2.76",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["validation_failure_taxonomy"] = {
        "hotspots": [
            {
                "category": "qc_sampling_issue",
                "complex_id": "1BJ1",
                "job_id": "1bj1-antigen-w-g88a",
                "diagnostic_code": "qc_repeat_spread",
                "abs_ddg_error_kcal_mol": 5.50,
                "ddg_ready": True,
            }
        ]
    }
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = _run(root, "--mode", "ultra")

    assert result.returncode == 0
    assert "1bj1-antigen-w-g88a" in result.stdout.splitlines()


def test_refresh_validation_watchlists_emits_stale_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "stale")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "1mlc-antibody-h-s57a",
        "1bj1-antigen-w-i80a",
        "1bj1-antigen-v-f17a",
    ]


def test_refresh_validation_watchlists_emits_stale_mode_across_candidate_plan_roots(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    robust_root = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
    )
    summary_path = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    hotspots = list(summary["active_alternate_ready_hotspots"])
    hotspots.extend(
        [
            {
                "job_id": "1bj1-antigen-w-g88a",
                "complex_id": "1BJ1",
                "active_alternate_source_plan_roots": str(robust_root),
                "abs_ddg_error_kcal_mol": 5.50,
            },
            {
                "job_id": "1bj1-antigen-w-h86a",
                "complex_id": "1BJ1",
                "active_alternate_source_plan_roots": str(robust_root),
                "abs_ddg_error_kcal_mol": 0.77,
            },
        ]
    )
    summary["active_alternate_ready_hotspots"] = hotspots
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = _run(
        root,
        "--mode",
        "stale",
        "--candidate-plan-root",
        str(robust_root),
    )

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "1mlc-antibody-h-s57a",
        "1bj1-antigen-w-g88a",
        "1bj1-antigen-w-h86a",
        "1bj1-antigen-w-i80a",
        "1bj1-antigen-v-f17a",
    ]


def test_refresh_validation_watchlists_emits_gap_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "gap")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "1mlc-antibody-h-s57a",
        "2nz9-antigen-a-f953a",
        "2nz9-antigen-a-t1063a",
        "1mlc-antibody-l-n92a",
    ]


def test_refresh_validation_watchlists_emits_backlog_mode(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)

    result = _run(root, "--mode", "backlog")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-l-y50a",
    ]


def test_refresh_validation_watchlists_gap_mode_excludes_default_queue_excluded_jobs(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    plan_jobs.write_text(
        plan_jobs.read_text(encoding="utf-8").replace(
            "2NZ9,2nz9-antigen-a-h1064a,mutate,blocked_input,False,False,7.42",
            "2NZ9,2nz9-antigen-a-h1064a,mutate,blocked_input,True,False,7.42",
        ),
        encoding="utf-8",
    )

    result = _run(root, "--mode", "gap")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "1mlc-antibody-h-s57a",
        "2nz9-antigen-a-f953a",
        "2nz9-antigen-a-t1063a",
        "1mlc-antibody-l-n92a",
    ]


def test_refresh_validation_watchlists_gap_mode_can_include_resumable_blocked_mutate_jobs_when_defaults_disabled(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    plan_jobs.write_text(
        plan_jobs.read_text(encoding="utf-8").replace(
            "2NZ9,2nz9-antigen-a-h1064a,mutate,blocked_input,False,False,7.42",
            "2NZ9,2nz9-antigen-a-h1064a,mutate,blocked_input,True,False,7.42",
        ),
        encoding="utf-8",
    )

    result = _run(root, "--mode", "gap", "--no-default-queue-excluded-job-ids")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "1mlc-antibody-h-s57a",
        "2nz9-antigen-a-f953a",
        "2nz9-antigen-a-t1063a",
        "1mlc-antibody-l-n92a",
        "2nz9-antigen-a-h1064a",
    ]


def test_refresh_validation_watchlists_gap_mode_excludes_derived_invalid_mutate_output_jobs(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    priority_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    robust_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    priority_plan_jobs.write_text(
        priority_plan_jobs.read_text(encoding="utf-8")
        + (
            "2NZ9,2nz9-antigen-a-q915a,mutate,blocked_input,False,False,1.15,"
            ",,,,,,,,,True,mutate_processed_gro_isolated_residue_hydrogen\n"
        ),
        encoding="utf-8",
    )
    robust_plan_jobs.write_text(
        robust_plan_jobs.read_text(encoding="utf-8")
        + "2NZ9,2nz9-antigen-a-q915a,build_legs,completed,True,False,1.15,,,,,,,,,,\n",
        encoding="utf-8",
    )

    result = _run(
        root,
        "--mode",
        "gap",
        "--candidate-plan-root",
        "runs/benchmarks/abbind_core_v1_validation_robust_plan",
    )

    assert result.returncode == 0
    assert "2nz9-antigen-a-q915a" not in result.stdout.splitlines()


def test_refresh_validation_watchlists_gap_mode_can_include_derived_invalid_mutate_output_jobs_when_opted_out(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    priority_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    robust_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    priority_plan_jobs.write_text(
        priority_plan_jobs.read_text(encoding="utf-8")
        + (
            "2NZ9,2nz9-antigen-a-q915a,mutate,blocked_input,False,False,1.15,"
            ",,,,,,,,,True,mutate_processed_gro_isolated_residue_hydrogen\n"
        ),
        encoding="utf-8",
    )
    robust_plan_jobs.write_text(
        robust_plan_jobs.read_text(encoding="utf-8")
        + "2NZ9,2nz9-antigen-a-q915a,build_legs,completed,True,False,1.15,,,,,,,,,,\n",
        encoding="utf-8",
    )

    result = _run(
        root,
        "--mode",
        "gap",
        "--candidate-plan-root",
        "runs/benchmarks/abbind_core_v1_validation_robust_plan",
        "--no-derived-invalid-mutate-output-exclusions",
    )

    assert result.returncode == 0
    assert "2nz9-antigen-a-q915a" in result.stdout.splitlines()


def test_refresh_validation_watchlists_stale_mode_excludes_default_queue_excluded_jobs(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    robust_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    robust_plan_jobs.write_text(
        robust_plan_jobs.read_text(encoding="utf-8")
        + "2NZ9,2nz9-antigen-a-h1064a,sample,stale_running,True,False,7.42,,,,,,,,,\n",
        encoding="utf-8",
    )

    result = _run(
        root,
        "--mode",
        "stale",
        "--candidate-plan-root",
        "runs/benchmarks/abbind_core_v1_validation_robust_plan",
    )

    assert result.returncode == 0
    assert "2nz9-antigen-a-h1064a" not in result.stdout.splitlines()


def test_refresh_validation_watchlists_stale_mode_can_include_queue_excluded_jobs_when_defaults_disabled(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    robust_plan_jobs = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_robust_plan"
        / "reports"
        / "plan_jobs.csv"
    )
    robust_plan_jobs.write_text(
        robust_plan_jobs.read_text(encoding="utf-8")
        + "2NZ9,2nz9-antigen-a-h1064a,sample,stale_running,True,False,7.42,,,,,,,,,\n",
        encoding="utf-8",
    )

    result = _run(
        root,
        "--mode",
        "stale",
        "--candidate-plan-root",
        "runs/benchmarks/abbind_core_v1_validation_robust_plan",
        "--no-default-queue-excluded-job-ids",
    )

    assert result.returncode == 0
    assert "2nz9-antigen-a-h1064a" in result.stdout.splitlines()


def test_refresh_validation_watchlists_prefers_canonical_merged_summary_when_present(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    canonical_merged_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
    )
    canonical_merged_dir.mkdir(parents=True, exist_ok=True)
    (canonical_merged_dir / "plan_summary.json").write_text(
        json.dumps(
            {
                "active_alternate_ready_hotspots": [
                    {
                        "job_id": "2nz9-antigen-a-f953a",
                        "complex_id": "2NZ9",
                        "active_alternate_source_plan_roots": str(
                            root
                            / "runs"
                            / "benchmarks"
                            / "abbind_core_v1_validation_robust_plan"
                        ),
                        "abs_ddg_error_kcal_mol": 3.34,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    (canonical_merged_dir / "benchmark_pairs.csv").write_text(
        "\n".join(
            [
                "complex_id,job_id,predicted_ddg_kcal_mol,experimental_ddg_kcal_mol",
                "2NZ9,2nz9-antigen-a-f953a,3.3,3.3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (canonical_merged_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,job_id,latest_stage,latest_stage_state,resumable,"
                    "ddg_ready,experimental_ddg_kcal_mol,active_alternate_candidate_count"
                ),
                "2NZ9,2nz9-antigen-a-f953a,build_legs,completed,True,False,3.34,0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run(root, "--mode", "robust")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "2nz9-antigen-a-f953a",
        "1mlc-antibody-h-s57a",
        "1mlc-antibody-h-s57v",
        "1mlc-antibody-h-t31a",
        "1mlc-antibody-h-t31v",
        "1mlc-antibody-l-n92a",
    ]


def test_refresh_validation_watchlists_prefers_newer_selection_over_stale_canonical(tmp_path: Path) -> None:
    root = _fake_root(tmp_path)
    canonical_merged_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
    )
    canonical_merged_dir.mkdir(parents=True, exist_ok=True)
    canonical_summary = canonical_merged_dir / "plan_summary.json"
    canonical_pairs = canonical_merged_dir / "benchmark_pairs.csv"
    canonical_plan_jobs = canonical_merged_dir / "plan_jobs.csv"
    canonical_summary.write_text(
        json.dumps(
            {
                "active_alternate_ready_hotspots": [
                    {
                        "job_id": "2nz9-antigen-a-f953a",
                        "complex_id": "2NZ9",
                        "active_alternate_source_plan_roots": str(
                            root
                            / "runs"
                            / "benchmarks"
                            / "abbind_core_v1_validation_robust_plan"
                        ),
                        "abs_ddg_error_kcal_mol": 3.34,
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    canonical_pairs.write_text(
        "\n".join(
            [
                "complex_id,job_id,predicted_ddg_kcal_mol,experimental_ddg_kcal_mol",
                "2NZ9,2nz9-antigen-a-f953a,3.3,3.3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    canonical_plan_jobs.write_text(
        "\n".join(
            [
                (
                    "complex_id,job_id,latest_stage,latest_stage_state,resumable,"
                    "ddg_ready,experimental_ddg_kcal_mol,active_alternate_candidate_count"
                ),
                "2NZ9,2nz9-antigen-a-f953a,build_legs,completed,True,False,3.34,0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    selection_summary = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    selection_mtime = selection_summary.stat().st_mtime_ns
    older_mtime_seconds = (selection_mtime - 1_000_000_000) / 1_000_000_000
    for path in (canonical_summary, canonical_pairs, canonical_plan_jobs):
        os.utime(path, (older_mtime_seconds, older_mtime_seconds))

    result = _run(root, "--mode", "robust")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1mlc-antibody-h-s57a",
        "1mlc-antibody-h-s57v",
        "1mlc-antibody-h-t31a",
        "1mlc-antibody-h-t31v",
        "1mlc-antibody-l-n92a",
    ]


def test_refresh_validation_watchlists_prefers_validation_split_selection_over_newer_complex_selection(
    tmp_path: Path,
) -> None:
    root = _fake_root(tmp_path)
    validation_summary = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "split-validation-fake"
        / "plan_summary.json"
    )
    validation_payload = json.loads(validation_summary.read_text(encoding="utf-8"))
    validation_payload["selection"] = {"split_name": "validation"}
    validation_summary.write_text(json.dumps(validation_payload), encoding="utf-8")

    complex_selection_dir = (
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "selections"
        / "complex-1bj1"
    )
    complex_selection_dir.mkdir(parents=True, exist_ok=True)
    (complex_selection_dir / "plan_summary.json").write_text(
        json.dumps(
            {
                "selection": {"split_name": ""},
                "active_alternate_ready_hotspots": [
                    {
                        "job_id": "2nz9-antigen-a-f953a",
                        "complex_id": "2NZ9",
                        "active_alternate_source_plan_roots": str(
                            root
                            / "runs"
                            / "benchmarks"
                            / "abbind_core_v1_validation_robust_plan"
                        ),
                        "abs_ddg_error_kcal_mol": 3.34,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    (complex_selection_dir / "benchmark_pairs.csv").write_text(
        "\n".join(
            [
                "complex_id,job_id,predicted_ddg_kcal_mol,experimental_ddg_kcal_mol",
                "2NZ9,2nz9-antigen-a-f953a,3.3,3.3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (complex_selection_dir / "plan_jobs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,job_id,latest_stage,latest_stage_state,resumable,"
                    "ddg_ready,experimental_ddg_kcal_mol,active_alternate_candidate_count"
                ),
                "2NZ9,2nz9-antigen-a-f953a,build_legs,completed,True,False,3.34,0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    newer_seconds = validation_summary.stat().st_mtime + 10
    for path in (
        complex_selection_dir / "plan_summary.json",
        complex_selection_dir / "benchmark_pairs.csv",
        complex_selection_dir / "plan_jobs.csv",
    ):
        os.utime(path, (newer_seconds, newer_seconds))

    result = _run(root, "--mode", "robust")

    assert result.returncode == 0
    assert result.stdout.splitlines() == [
        "3hfm-antibody-h-y33a",
        "3hfm-antigen-y-y20a",
        "1mlc-antibody-h-s57a",
        "1mlc-antibody-h-s57v",
        "1mlc-antibody-h-t31a",
        "1mlc-antibody-h-t31v",
        "1mlc-antibody-l-n92a",
    ]
