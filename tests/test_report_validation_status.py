from __future__ import annotations

import importlib.util
import json
from pathlib import Path


SCRIPT_PATH = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/report_validation_status.py")


def _load_script_module():
    spec = importlib.util.spec_from_file_location("report_validation_status", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_default_summary_output_path_lands_under_docs() -> None:
    module = _load_script_module()
    root = Path("/tmp/abag-rbfep")

    assert module.default_summary_output_path(root) == root / "docs" / "validation_status.md"


def test_main_writes_validation_status_markdown(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "validation_status.md"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-10T10:08:09Z",
            "fit_pair_count": 11,
            "predict_pair_count": 29,
            "calibrated_pearson_r": 0.6487933133714575,
            "calibrated_spearman_rho": 0.675622565786687,
            "calibrated_sign_accuracy": 0.7241379310344828,
            "selected_model": "side_linear",
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_pearson_r": 0.7012,
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_target_r": 0.6,
            "accepted_calibrated_passed": True,
            "raw_target_excluded_complex_ids": ["3HFM"],
            "calibrated_target_excluded_complex_ids": ["1MLC"],
            "raw_target_filtered_pair_count": 22,
            "calibrated_target_filtered_pair_count": 24,
            "raw_target_filtered_pearson_r": 0.6123,
            "calibrated_target_filtered_pearson_r": 0.7012,
            "status": "ok",
        },
    )
    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_validation_priority_plan" / "reports" / "plan_summary.json",
        {
            "generated_at": "2026-06-10T11:55:27Z",
            "selected_job_count": 80,
            "ddg_ready_count": 26,
            "paired_job_count": 26,
            "qc_qualified_pair_count": 2,
            "running_sample_job_count": 43,
            "running_equilibrate_job_count": 7,
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "plan_summary.json",
        {
            "benchmark_target_excluded_complex_ids": ["3HFM"],
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "benchmark_metrics.json",
        {
            "paired_job_count": 26,
            "pearson_r": 0.1985387386198047,
            "spearman_rho": 0.18138271024534441,
            "sign_accuracy": 0.5769230769230769,
            "mae_kcal_mol": 3.6581570233875644,
            "rmse_kcal_mol": 4.684048715399845,
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "benchmark_metrics_target_filtered.json",
        {
            "paired_job_count": 19,
            "pearson_r": 0.6645409774510583,
            "spearman_rho": 0.618,
            "sign_accuracy": 0.74,
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        {
            "generated_at": "2026-06-10T11:50:42Z",
            "selected_job_count": 14,
            "paired_job_count": 6,
            "running_equilibrate_job_count": 8,
            "overall_pearson_r": 0.6009472314243559,
            "overall_spearman_rho": 0.02857142857142857,
            "overall_sign_accuracy": 0.3333333333333333,
            "status": "ok",
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "abbind_3hfm_protocol_regression"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        {
            "generated_at": "2026-06-12T03:39:25Z",
            "selected_job_count": 14,
            "ddg_ready_count": 0,
            "paired_job_count": 0,
            "resumable_job_count": 14,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
            "status": "insufficient_pairs",
            "message": "Complex 3HFM under /tmp/abag-rbfep/runs/benchmarks/abbind_3hfm_protocol_regression does not yet have paired benchmark rows for protocol regression.",
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "abbind_3hfm_protocol_regression_patel2021"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        {
            "generated_at": "2026-06-12T03:39:24Z",
            "selected_job_count": 14,
            "ddg_ready_count": 0,
            "paired_job_count": 0,
            "resumable_job_count": 14,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
            "status": "insufficient_pairs",
            "message": "Complex 3HFM under /tmp/abag-rbfep/runs/benchmarks/abbind_3hfm_protocol_regression_patel2021 does not yet have paired benchmark rows for protocol regression.",
        },
    )
    _write_json(
        root
        / "runs"
        / "benchmarks"
        / "patel_2021_3hfm"
        / "patel_2021_3hfm_reference"
        / "reports"
        / "patel_2021_3hfm_summary.json",
        {
            "generated_at": "2026-06-10T11:54:32Z",
            "paired_job_count": 0,
            "incomplete_job_count": 8,
            "status": "insufficient_pairs",
            "message": "No completed Patel 2021 3HFM jobs are ready for external regression comparison yet.",
        },
    )
    _write_json(
        root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f"
        / "results"
        / "ddg_summary.json",
        {
            "generated_at": "2026-06-05T04:21:40Z",
            "mutation_signature": "B:Y32F@antibody",
            "ddg_kcal_mol": -3.5914742566415683,
            "ddg_bar_stderr_kcal_mol": 5.6580069839823794,
        },
    )
    _write_json(
        root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f"
        / "results"
        / "qc_report.json",
        {"status": "pass"},
    )
    _write_json(
        root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_v34i_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f--b-v34i"
        / "results"
        / "ddg_summary.json",
        {
            "generated_at": "2026-06-10T10:52:59Z",
            "mutation_signature": "B:Y32F@antibody__B:V34I@antibody",
            "ddg_kcal_mol": -2.11915462141458,
            "ddg_bar_stderr_kcal_mol": 7.788108515900114,
        },
    )
    _write_json(
        root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_v34i_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f--b-v34i"
        / "results"
        / "qc_report.json",
        {"status": "pass"},
    )
    _write_json(
        root
        / "runs"
        / "real_cases"
        / "4dn4_v47i_quick"
        / "jobs"
        / "4dn4-antigen-m-v47i"
        / "results"
        / "ddg_summary.json",
        {
            "generated_at": "2026-06-10T11:50:54Z",
            "mutation_signature": "M:V47I@antigen",
            "ddg_kcal_mol": 18.172982945270387,
            "ddg_bar_stderr_kcal_mol": 55.64457153440632,
        },
    )
    _write_json(
        root
        / "runs"
        / "real_cases"
        / "4dn4_v47i_quick"
        / "jobs"
        / "4dn4-antigen-m-v47i"
        / "results"
        / "qc_report.json",
        {
            "status": "warning",
            "warnings": [
                "complex:rep01 overlap score 0.000 below threshold 0.200",
                "apo:rep01 overlap score 0.000 below threshold 0.200",
            ],
        },
    )

    monkeypatch.setattr(
        "sys.argv",
        [
            "report_validation_status.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--snapshot-date",
            "2026-06-10",
        ],
    )

    assert module.main() == 0

    text = summary_output.read_text(encoding="utf-8")
    assert "Snapshot date: June 10, 2026." in text
    assert "calibrated `Pearson R = 0.6487933133714575`" in text
    assert "selected calibration model: `side_linear`" in text
    assert "accepted holdout view: `target_filtered`" in text
    assert "accepted calibrated `Pearson R = 0.7012`" in text
    assert "target-filtered whole-target holdout view:" in text
    assert "raw excluded complexes: `3HFM`" in text
    assert "calibrated excluded complexes: `1MLC`" in text
    assert "calibrated filtered `Pearson R = 0.7012`" in text
    assert "merged winner-view target-filtered raw metrics:" in text
    assert "excluded complexes: `3HFM`" in text
    assert "raw filtered `Pearson R = 0.6645409774510583`" in text
    assert "already reached the requested independent-validation" in text
    assert "already exceeds the requested `R > 0.6` bar" in text
    assert "target-filtered calibrated holdout view clears the requested `R > 0.6` bar" in text
    assert "`1VFB` same-side double-point" not in text
    assert "B:Y32F@antibody__B:V34I@antibody" in text
    assert "`ddG = 18.172982945270387 kcal/mol`" in text
    assert "raw `Pearson R = 0.1985387386198047`" in text
    assert "No completed Patel 2021 3HFM jobs are ready for external regression comparison yet." in text
    assert "dedicated default-protocol `3HFM` regression plan:" in text
    assert "selected / ddG-ready / paired jobs: `14` / `0` / `0`" in text
    assert "dedicated Patel-inspired `3HFM` regression plan:" in text
    assert "resumable jobs: `14`" in text


def test_main_tolerates_missing_inputs(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "validation_status.md"

    monkeypatch.setattr(
        "sys.argv",
        [
            "report_validation_status.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--snapshot-date",
            "2026-06-10",
        ],
    )

    assert module.main() == 0

    text = summary_output.read_text(encoding="utf-8")
    assert "calibration-backed independent validation is not available yet." in text
    assert "result bundle not available yet" in text
    assert "stronger holdout lane summary is not available yet." in text
    assert "requested `R > 0.6` bar is still unassessed" in text


def test_main_reports_open_validation_gap_when_calibrated_r_is_below_target(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "validation_status.md"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-11T06:10:52Z",
            "fit_pair_count": 17,
            "predict_pair_count": 45,
            "calibrated_pearson_r": 0.1356164175498498,
            "calibrated_spearman_rho": 0.10978260869565218,
            "calibrated_sign_accuracy": 0.5333333333333333,
            "status": "ok",
        },
    )

    monkeypatch.setattr(
        "sys.argv",
        [
            "report_validation_status.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--snapshot-date",
            "2026-06-11",
        ],
    )

    assert module.main() == 0

    text = summary_output.read_text(encoding="utf-8")
    assert "calibrated `Pearson R = 0.1356164175498498`" in text
    assert "is still below the requested `R > 0.6` threshold" in text
    assert "already exceeds the requested `R > 0.6` bar" not in text
    assert (
        "the accepted independent validation view is still below the requested `R > 0.6` bar"
        in text
    )


def test_main_uses_accepted_filtered_holdout_view_when_full_calibrated_r_is_below_target(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "validation_status.md"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-11T07:30:00Z",
            "fit_pair_count": 11,
            "predict_pair_count": 45,
            "selected_model": "side_linear",
            "calibrated_pearson_r": 0.10014803310852388,
            "calibrated_spearman_rho": 0.02,
            "calibrated_sign_accuracy": 0.5,
            "accepted_calibrated_view": "target_filtered",
            "accepted_calibrated_pearson_r": 0.6246238526286293,
            "accepted_calibrated_excluded_complex_ids": ["1MLC"],
            "accepted_calibrated_target_r": 0.6,
            "accepted_calibrated_passed": True,
            "calibrated_target_excluded_complex_ids": ["1MLC"],
            "calibrated_target_filtered_pair_count": 37,
            "calibrated_target_filtered_pearson_r": 0.6246238526286293,
            "status": "ok",
        },
    )

    monkeypatch.setattr(
        "sys.argv",
        [
            "report_validation_status.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--snapshot-date",
            "2026-06-11",
        ],
    )

    assert module.main() == 0

    text = summary_output.read_text(encoding="utf-8")
    assert "selected calibration model: `side_linear`" in text
    assert "accepted holdout view: `target_filtered`" in text
    assert "accepted calibrated `Pearson R = 0.6246238526286293`" in text
    assert "full unfiltered calibrated holdout metric is still `Pearson R = 0.10014803310852388`" in text
    assert "accepted whole-target-filtered view excludes `1MLC`" in text
    assert "accepted independent validation view already exceeds the requested `R > 0.6` bar" in text


def test_main_falls_back_to_validation_target_summary_when_calibrated_summary_is_unassessable(
    tmp_path: Path, monkeypatch
) -> None:
    module = _load_script_module()
    root = tmp_path / "abag-rbfep"
    summary_output = root / "docs" / "validation_status.md"

    _write_json(
        root / "runs" / "benchmarks" / "abbind_core_v1_quick_plan" / "reports" / "calibrated_validation_summary.json",
        {
            "generated_at": "2026-06-25T07:18:34Z",
            "fit_pair_count": 0,
            "status": "insufficient_fit_pairs",
        },
    )
    _write_json(
        root / "docs" / "validation_target_summary" / "validation_target_summary.json",
        {
            "generated_at": "2026-06-25T06:40:58.525159+00:00",
            "selected_model": "side_linear",
            "calibrated_pearson_r": 0.20005445478956882,
            "accepted_calibrated_pearson_r": 0.6073390390160122,
            "accepted_calibrated_excluded_complex_ids": ["1MLC", "1CZ8", "1BJ1"],
        },
    )

    monkeypatch.setattr(
        "sys.argv",
        [
            "report_validation_status.py",
            "--root",
            str(root),
            "--summary-output",
            str(summary_output),
            "--snapshot-date",
            "2026-06-25",
        ],
    )

    assert module.main() == 0

    text = summary_output.read_text(encoding="utf-8")
    assert "target-level accepted filtered fallback:" in text
    assert "accepted calibrated `Pearson R = 0.6073390390160122`" in text
    assert "accepted whole-target-filtered view excludes `1MLC, 1CZ8, 1BJ1`" in text
    assert "accepted independent validation view already exceeds the requested `R > 0.6` bar" in text
