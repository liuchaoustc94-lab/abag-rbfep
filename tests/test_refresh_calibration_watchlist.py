from __future__ import annotations

import importlib.util
from pathlib import Path

from abag_rbfe.io_utils import write_csv_rows, write_json


def _load_module():
    module_path = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/refresh_calibration_watchlist.py")
    spec = importlib.util.spec_from_file_location("refresh_calibration_watchlist", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_load_fit_pair_counts_reads_existing_summary(tmp_path: Path) -> None:
    module = _load_module()
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir()
    fit_pairs = reports_dir / "fit_pairs.csv"
    write_csv_rows(
        fit_pairs,
        [
            {
                "job_id": "a-antibody-h-y32a",
                "predicted_ddg_kcal_mol": "2.0",
                "experimental_ddg_kcal_mol": "-0.2",
            },
            {
                "job_id": "b-antibody-h-y33a",
                "predicted_ddg_kcal_mol": "3.0",
                "experimental_ddg_kcal_mol": "0.6",
            },
            {
                "job_id": "c-antigen-w-g92a",
                "predicted_ddg_kcal_mol": "4.0",
                "experimental_ddg_kcal_mol": "2.4",
            },
        ],
        ["job_id", "predicted_ddg_kcal_mol", "experimental_ddg_kcal_mol"],
    )
    summary_path = tmp_path / "calibrated_validation_summary.json"
    write_json(summary_path, {"reports_dir": str(reports_dir)})

    bin_counts, side_counts = module._load_fit_pair_counts(summary_path, near_zero_threshold=1.0)

    assert side_counts == {"antibody": 2, "antigen": 1, "unknown": 0}
    assert bin_counts[("antibody", "negative")] == 1
    assert bin_counts[("antibody", "near_zero")] == 1
    assert bin_counts[("antigen", "positive")] == 1


def test_select_watchlist_job_rows_balances_sparse_bins_and_keeps_active_jobs() -> None:
    module = _load_module()
    candidates = [
        {
            "job_id": "job-running-antibody-positive",
            "entity_side": "antibody",
            "effect_bin": "positive",
            "experimental_ddg_value": 2.2,
            "fit_split_name": "calibration",
            "latest_stage_state": "running",
        },
        {
            "job_id": "job-antigen-negative",
            "entity_side": "antigen",
            "effect_bin": "negative",
            "experimental_ddg_value": -0.4,
            "fit_split_name": "calibration",
            "latest_stage_state": "not_started",
        },
        {
            "job_id": "job-antigen-near-zero",
            "entity_side": "antigen",
            "effect_bin": "near_zero",
            "experimental_ddg_value": 0.2,
            "fit_split_name": "development",
            "latest_stage_state": "not_started",
        },
        {
            "job_id": "job-antibody-positive",
            "entity_side": "antibody",
            "effect_bin": "positive",
            "experimental_ddg_value": 3.1,
            "fit_split_name": "calibration",
            "latest_stage_state": "not_started",
        },
    ]
    current_bin_counts = {
        ("antibody", "negative"): 1,
        ("antibody", "near_zero"): 3,
        ("antibody", "positive"): 7,
        ("antigen", "negative"): 0,
        ("antigen", "near_zero"): 1,
        ("antigen", "positive"): 4,
        ("unknown", "negative"): 0,
        ("unknown", "near_zero"): 0,
        ("unknown", "positive"): 0,
    }
    current_side_counts = {"antibody": 11, "antigen": 5, "unknown": 0}

    selected = module._select_watchlist_job_rows(
        candidates,
        current_bin_counts=current_bin_counts,
        current_side_counts=current_side_counts,
        target_job_count=3,
    )

    assert [row["job_id"] for row in selected] == [
        "job-running-antibody-positive",
        "job-antigen-negative",
        "job-antigen-near-zero",
    ]
    assert selected[0]["selection_round"] == 1
    assert selected[0]["priority_score"][:4] == [0, 7, 11, 0]
    assert selected[1]["selection_round"] == 2
    assert selected[1]["selection_priority"]["fit_bin_count"] == 0
    assert selected[2]["selection_priority"]["fit_side_count"] == 6


def test_build_watchlist_payload_reports_fit_shortfall_and_selected_counts(tmp_path: Path) -> None:
    module = _load_module()
    payload = module._build_watchlist_payload(
        plan_root=tmp_path / "runs",
        summary_path=tmp_path / "summary.json",
        fit_split_names=["calibration", "development"],
        near_zero_threshold=1.0,
        current_bin_counts={
            ("antibody", "negative"): 0,
            ("antibody", "near_zero"): 1,
            ("antibody", "positive"): 2,
            ("antigen", "negative"): 0,
            ("antigen", "near_zero"): 0,
            ("antigen", "positive"): 1,
            ("unknown", "negative"): 0,
            ("unknown", "near_zero"): 0,
            ("unknown", "positive"): 0,
        },
        current_side_counts={"antibody": 3, "antigen": 1, "unknown": 0},
        target_job_count=4,
        candidates=[{"job_id": "job-a"}, {"job_id": "job-b"}],
        selected_rows=[
            {
                "job_id": "job-a",
                "complex_id": "1ABC",
                "fit_split_name": "calibration",
                "entity_side": "antibody",
                "effect_bin": "negative",
                "experimental_ddg_value": -0.4,
                "latest_stage": "",
                "latest_stage_state": "not_started",
                "selection_round": 1,
                "priority_score": [3, 0, 3, 0, -0.4, 0.4],
                "selection_priority": {"stage_rank": 3, "fit_bin_count": 0, "fit_side_count": 3},
                "ddg_ready": "False",
            }
        ],
        selected_job_ids=["job-a"],
    )

    assert payload["fit_pair_count"] == 4
    assert payload["fit_pair_shortfall"] == 3
    assert payload["selected_job_count"] == 1
    assert payload["selected_ready_job_count"] == 0
    assert payload["selected_rows"][0]["fit_entity_side"] == "antibody"
    assert payload["selected_rows"][0]["fit_effect_bin"] == "negative"
