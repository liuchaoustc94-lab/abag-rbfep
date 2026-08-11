from __future__ import annotations

import importlib.util
import json
import os
from datetime import datetime
from pathlib import Path


SCRIPT_PATH = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/report_calibrated_validation.py")


def _load_script_module():
    spec = importlib.util.spec_from_file_location("report_calibrated_validation", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_default_extra_roots_include_calibration_rescues() -> None:
    module = _load_script_module()

    roots = [Path(item) for item in module.DEFAULT_EXTRA_PLAN_ROOTS]

    assert SCRIPT_PATH.resolve().parents[2] / "runs" / "benchmarks" / "abbind_core_v1_calibration_rescues" in roots
    assert (
        SCRIPT_PATH.resolve().parents[2]
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_targeted_repeat_spread_rescues"
        in roots
    )
    assert (
        SCRIPT_PATH.resolve().parents[2]
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_targeted_lambda_rescues"
        in roots
    )
    assert (
        SCRIPT_PATH.resolve().parents[2]
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_sampling_qc_rescues"
        in roots
    )
    assert SCRIPT_PATH.resolve().parents[2] / "runs" / "benchmarks" / "abbind_core_v1_validation_deep_rescues" in roots
    assert SCRIPT_PATH.resolve().parents[2] / "runs" / "benchmarks" / "abbind_core_v1_validation_ultra_rescues" in roots


def test_default_fit_split_defaults_include_development() -> None:
    module = _load_script_module()

    assert module.DEFAULT_FIT_SPLIT_NAME == "calibration"
    assert module.DEFAULT_FIT_EXTRA_SPLIT_NAMES == ["development"]


def test_default_summary_output_path_lands_under_plan_root_reports() -> None:
    module = _load_script_module()
    plan_root = Path("/tmp/abag-rbfep/runs/benchmarks/abbind_core_v1_quick_plan")

    assert module.default_summary_output_path(plan_root) == plan_root / "reports" / "calibrated_validation_summary.json"


def test_parse_args_defaults_to_auto(monkeypatch) -> None:
    module = _load_script_module()

    monkeypatch.setattr("sys.argv", ["report_calibrated_validation.py"])

    args = module.parse_args()

    assert args.model == "auto"
    assert args.no_use_existing_selection_reports is False


def test_main_writes_model_into_summary(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    reports_dir = plan_root / "reports"
    summary_output = reports_dir / "calibrated_validation_summary.json"
    plan_root.mkdir(parents=True)

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir),
            "paired_job_count": 6,
            "qc_qualified_pair_count": 0,
        },
    )
    monkeypatch.setattr(
        module,
        "calibrate_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir / "calibrations" / "mock"),
            "source_plan_roots": [str(plan_root)],
            "fit_split_names": ["calibration", "development"],
            "fit_pair_count": 6,
            "predict_pair_count": 3,
            "fit_coverage": {"pair_count": 6},
            "predict_raw_coverage": {"pair_count": 3},
            "predict_calibrated_coverage": {"pair_count": 3},
            "raw_metrics": {"pearson_r": 0.1, "spearman_rho": 0.2, "sign_accuracy": 0.3},
            "calibrated_metrics": {"pearson_r": 0.7, "spearman_rho": 0.6, "sign_accuracy": 0.8},
            "predict_target_exclusion_policy": {
                "target_field": "complex_id",
                "systematically_poor_abs_error_threshold_kcal_mol": 2.0,
                "systematically_poor_min_pair_count": 4,
            },
            "predict_raw_target_excluded_complex_ids": ["BAD1"],
            "predict_calibrated_target_excluded_complex_ids": [],
            "predict_raw_target_filtered_metrics": {
                "paired_job_count": 2,
                "pearson_r": 0.55,
                "spearman_rho": 0.45,
                "sign_accuracy": 0.5,
            },
            "predict_raw_target_filtered_outlier_trimmed_metrics": {
                "paired_job_count": 2,
                "pearson_r": 0.61,
                "spearman_rho": 0.51,
                "sign_accuracy": 0.6,
            },
            "predict_calibrated_target_filtered_metrics": {
                "paired_job_count": 3,
                "pearson_r": 0.72,
                "spearman_rho": 0.62,
                "sign_accuracy": 0.82,
            },
            "predict_calibrated_target_filtered_outlier_trimmed_metrics": {
                "paired_job_count": 3,
                "pearson_r": 0.74,
                "spearman_rho": 0.64,
                "sign_accuracy": 0.84,
            },
            "predict_raw_outlier_trimmed_metrics": {
                "paired_job_count": 3,
                "pearson_r": 0.21,
            },
            "predict_calibrated_outlier_trimmed_metrics": {
                "paired_job_count": 3,
                "pearson_r": 0.71,
            },
            "predict_outlier_trim_policy": {"outlier_trim_method": "tukey_iqr"},
            "predict_raw_target_metrics": [
                {"complex_id": "BAD1", "systematically_poor_target": True},
                {"complex_id": "GOOD1", "systematically_poor_target": False},
            ],
            "predict_raw_target_outlier_trim_metrics": [
                {"complex_id": "BAD1", "removed_pair_count": 1},
                {"complex_id": "GOOD1", "removed_pair_count": 0},
            ],
            "predict_calibrated_target_metrics": [
                {"complex_id": "GOOD1", "systematically_poor_target": False},
            ],
            "predict_calibrated_target_outlier_trim_metrics": [
                {"complex_id": "GOOD1", "removed_pair_count": 0},
            ],
            "model": {"model": "stderr_quadratic", "groups": {"global": {"family": "stderr_quadratic"}}},
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_calibrated_validation.py",
            "--plan-root",
            str(plan_root),
            "--model",
            "stderr_quadratic",
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 0

    summary = json.loads(summary_output.read_text(encoding="utf-8"))

    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["requested_model"] == "stderr_quadratic"
    assert summary["model_selection_mode"] == "explicit"
    assert summary["model"] == "stderr_quadratic"
    assert summary["selected_model"] == "stderr_quadratic"
    assert summary["calibrated_pearson_r"] == 0.7
    assert summary["predict_target_exclusion_policy"]["target_field"] == "complex_id"
    assert summary["raw_target_excluded_complex_ids"] == ["BAD1"]
    assert summary["calibrated_target_excluded_complex_ids"] == []
    assert summary["raw_target_filtered_pair_count"] == 2
    assert summary["calibrated_target_filtered_pair_count"] == 3
    assert summary["raw_target_filtered_outlier_trimmed_pair_count"] == 2
    assert summary["calibrated_target_filtered_outlier_trimmed_pair_count"] == 3
    assert summary["raw_target_filtered_pearson_r"] == 0.55
    assert summary["calibrated_target_filtered_pearson_r"] == 0.72
    assert summary["raw_target_filtered_outlier_trimmed_pearson_r"] == 0.61
    assert summary["calibrated_target_filtered_outlier_trimmed_pearson_r"] == 0.74
    assert summary["predict_outlier_trim_policy"]["outlier_trim_method"] == "tukey_iqr"
    assert summary["predict_raw_target_metrics"][0]["complex_id"] == "BAD1"
    assert summary["predict_raw_target_outlier_trim_metrics"][0]["complex_id"] == "BAD1"
    assert summary["predict_calibrated_target_metrics"][0]["complex_id"] == "GOOD1"
    assert summary["predict_calibrated_target_outlier_trim_metrics"][0]["complex_id"] == "GOOD1"
    assert summary["accepted_calibrated_view"] == "full"
    assert summary["accepted_calibrated_pearson_r"] == 0.7
    assert summary["accepted_calibrated_excluded_complex_ids"] == []
    assert summary["accepted_calibrated_passed"] is True
    assert summary["accepted_independent_view"] == "calibrated_target_filtered_outlier_trimmed"
    assert summary["accepted_independent_pearson_r"] == 0.74
    assert summary["accepted_independent_passed"] is True


def test_main_reuses_existing_selection_reports(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    reports_dir = plan_root / "reports"
    summary_output = reports_dir / "calibrated_validation_summary.json"
    split_path = tmp_path / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
    split_path.parent.mkdir(parents=True)
    split_path.write_text("splits: {}\n", encoding="utf-8")
    plan_root.mkdir(parents=True)

    def write_selection_bundle(name: str, split_name: str, paired_job_count: int) -> Path:
        bundle_dir = reports_dir / "merged" / "selections" / name
        bundle_dir.mkdir(parents=True)
        payload = {
            "reports_dir": str(bundle_dir),
            "paired_job_count": paired_job_count,
            "qc_qualified_pair_count": 0,
            "selection": {
                "split_name": split_name,
                "split_path": str(split_path.resolve()),
            },
            "source_plan_roots": [str(plan_root.resolve())],
        }
        (bundle_dir / "plan_summary.json").write_text(json.dumps(payload), encoding="utf-8")
        (bundle_dir / "plan_jobs.csv").write_text("job_id\nmock\n", encoding="utf-8")
        (bundle_dir / "benchmark_pairs.csv").write_text("job_id\nmock\n", encoding="utf-8")
        (bundle_dir / "benchmark_pairs_qc_qualified.csv").write_text("job_id\n", encoding="utf-8")
        return bundle_dir

    development_dir = write_selection_bundle("split-development", "development", 5)
    calibration_dir = write_selection_bundle("split-calibration", "calibration", 6)
    validation_dir = write_selection_bundle("split-validation", "validation", 3)

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("report_ab_bind_plan should not be called")),
    )
    monkeypatch.setattr(
        module,
        "calibrate_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir / "calibrations" / "mock"),
            "source_plan_roots": [str(plan_root.resolve())],
            "fit_split_names": ["calibration", "development"],
            "fit_pair_count": 11,
            "predict_pair_count": 3,
            "fit_coverage": {"pair_count": 11},
            "predict_raw_coverage": {"pair_count": 3},
            "predict_calibrated_coverage": {"pair_count": 3},
            "raw_metrics": {"pearson_r": 0.1, "spearman_rho": 0.2, "sign_accuracy": 0.3},
            "calibrated_metrics": {"pearson_r": 0.7, "spearman_rho": 0.6, "sign_accuracy": 0.8},
            "model": {"model": "stderr_quadratic", "groups": {"global": {"family": "stderr_quadratic"}}},
            "fit_reports_dirs": [str(calibration_dir), str(development_dir)],
            "predict_reports_dir": str(validation_dir),
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_calibrated_validation.py",
            "--plan-root",
            str(plan_root),
            "--split-file",
            str(split_path),
            "--model",
            "stderr_quadratic",
            "--no-default-extra-roots",
            "--no-default-fit-extra-splits",
            "--fit-extra-split-name",
            "development",
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 0

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["fit_pair_count"] == 11
    assert summary["predict_pair_count"] == 3
    assert summary["fit_reports_dirs"] == [str(calibration_dir), str(development_dir)]
    assert summary["predict_reports_dir"] == str(validation_dir)


def test_resolve_existing_selection_bundle_rejects_stale_source_reports(tmp_path: Path) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    reports_dir = plan_root / "reports"
    split_path = tmp_path / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
    split_path.parent.mkdir(parents=True)
    split_path.write_text("splits: {}\n", encoding="utf-8")
    plan_root.mkdir(parents=True)

    bundle_dir = reports_dir / "merged" / "selections" / "split-validation"
    bundle_dir.mkdir(parents=True)
    payload = {
        "reports_dir": str(bundle_dir),
        "paired_job_count": 3,
        "qc_qualified_pair_count": 0,
        "selection": {
            "split_name": "validation",
            "split_path": str(split_path.resolve()),
        },
        "source_plan_roots": [str(plan_root.resolve())],
    }
    required_files = [
        bundle_dir / "plan_summary.json",
        bundle_dir / "plan_jobs.csv",
        bundle_dir / "benchmark_pairs.csv",
        bundle_dir / "benchmark_pairs_qc_qualified.csv",
    ]
    required_files[0].write_text(json.dumps(payload), encoding="utf-8")
    required_files[1].write_text("job_id\nmock\n", encoding="utf-8")
    required_files[2].write_text("job_id\nmock\n", encoding="utf-8")
    required_files[3].write_text("job_id\n", encoding="utf-8")
    for path in required_files:
        os.utime(path, (1, 1))

    source_report = reports_dir / "plan_jobs.csv"
    source_report.parent.mkdir(parents=True, exist_ok=True)
    source_report.write_text("job_id\nnewer\n", encoding="utf-8")

    bundle = module.resolve_existing_selection_bundle(
        plan_root,
        extra_roots=[],
        split_name="validation",
        split_path=split_path.resolve(),
    )

    assert bundle is None


def test_main_can_force_report_refresh(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    reports_dir = plan_root / "reports"
    summary_output = reports_dir / "calibrated_validation_summary.json"
    split_path = tmp_path / "splits" / "ab_bind_rbfe_core_v1_split_v1.yml"
    split_path.parent.mkdir(parents=True)
    split_path.write_text("splits: {}\n", encoding="utf-8")
    plan_root.mkdir(parents=True)

    (reports_dir / "merged" / "selections" / "split-validation").mkdir(parents=True)

    report_calls: list[str] = []

    def fake_report_ab_bind_plan(_plan_root, **kwargs):
        split_name = kwargs["split_name"]
        report_calls.append(split_name)
        return {
            "reports_dir": str(reports_dir / "generated" / split_name),
            "paired_job_count": 4 if split_name != "validation" else 2,
            "qc_qualified_pair_count": 0,
        }

    monkeypatch.setattr(module, "report_ab_bind_plan", fake_report_ab_bind_plan)
    monkeypatch.setattr(
        module,
        "calibrate_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir / "calibrations" / "mock"),
            "source_plan_roots": [str(plan_root.resolve())],
            "fit_split_names": ["calibration", "development"],
            "fit_pair_count": 8,
            "predict_pair_count": 2,
            "fit_coverage": {"pair_count": 8},
            "predict_raw_coverage": {"pair_count": 2},
            "predict_calibrated_coverage": {"pair_count": 2},
            "raw_metrics": {"pearson_r": 0.1, "spearman_rho": 0.2, "sign_accuracy": 0.3},
            "calibrated_metrics": {"pearson_r": 0.7, "spearman_rho": 0.6, "sign_accuracy": 0.8},
            "model": {"model": "stderr_quadratic", "groups": {"global": {"family": "stderr_quadratic"}}},
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_calibrated_validation.py",
            "--plan-root",
            str(plan_root),
            "--split-file",
            str(split_path),
            "--model",
            "stderr_quadratic",
            "--no-default-extra-roots",
            "--no-default-fit-extra-splits",
            "--fit-extra-split-name",
            "development",
            "--no-use-existing-selection-reports",
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 0
    assert report_calls == ["calibration", "development", "validation"]


def test_main_auto_selects_best_model_by_accepted_calibrated_pearson(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    reports_dir = plan_root / "reports"
    summary_output = reports_dir / "calibrated_validation_summary.json"
    plan_root.mkdir(parents=True)

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir / kwargs["split_name"]),
            "paired_job_count": 6,
            "qc_qualified_pair_count": 0,
        },
    )

    def fake_payload(model_name: str) -> dict:
        accepted = {
            "linear": (0.22, [], 45),
            "side_linear": (0.6246238526286293, ["1MLC"], 37),
            "quadratic": (0.11, ["3HFM"], 38),
            "stderr_quadratic": (-0.2, ["3HFM"], 38),
            "logabs_stderr_quadratic": (-0.1, ["3HFM"], 38),
            "expdecay_invstderr_quadratic": (0.05, ["3HFM"], 38),
            "hill_invstderr_quadratic": (0.31, ["1BJ1"], 39),
            "hill_side_invstderr_quadratic": (0.40, ["1MLC"], 37),
        }[model_name]
        accepted_r, excluded, pair_count = accepted
        return {
            "reports_dir": str(reports_dir / "calibrations" / model_name),
            "source_plan_roots": [str(plan_root)],
            "fit_split_names": ["calibration", "development"],
            "fit_reports_dir": str(reports_dir / "calibration"),
            "fit_reports_dirs": [str(reports_dir / "calibration"), str(reports_dir / "development")],
            "predict_reports_dir": str(reports_dir / "validation"),
            "fit_pair_count": 11,
            "predict_pair_count": 45,
            "fit_coverage": {"pair_count": 11},
            "predict_raw_coverage": {"pair_count": 45},
            "predict_calibrated_coverage": {"pair_count": pair_count},
            "raw_metrics": {"pearson_r": 0.04, "spearman_rho": 0.2, "sign_accuracy": 0.55},
            "calibrated_metrics": {"pearson_r": 0.1 if model_name == "side_linear" else 0.08, "spearman_rho": 0.2, "sign_accuracy": 0.5},
            "predict_target_exclusion_policy": {
                "target_field": "complex_id",
                "systematically_poor_abs_error_threshold_kcal_mol": 2.0,
                "systematically_poor_min_pair_count": 4,
            },
            "predict_raw_target_excluded_complex_ids": ["3HFM"],
            "predict_calibrated_target_excluded_complex_ids": excluded,
            "predict_raw_target_filtered_metrics": {
                "paired_job_count": 38,
                "pearson_r": 0.2493,
                "spearman_rho": 0.47,
                "sign_accuracy": 0.60,
            },
            "predict_raw_target_filtered_outlier_trimmed_metrics": {
                "paired_job_count": 37,
                "pearson_r": 0.31,
                "spearman_rho": 0.49,
                "sign_accuracy": 0.62,
            },
            "predict_calibrated_target_filtered_metrics": {
                "paired_job_count": pair_count,
                "pearson_r": accepted_r,
                "spearman_rho": 0.5,
                "sign_accuracy": 0.6,
            },
            "predict_calibrated_target_filtered_outlier_trimmed_metrics": {
                "paired_job_count": pair_count - 1,
                "pearson_r": accepted_r + 0.02,
                "spearman_rho": 0.52,
                "sign_accuracy": 0.62,
            },
            "predict_raw_outlier_trimmed_metrics": {"paired_job_count": 44, "pearson_r": 0.12},
            "predict_calibrated_outlier_trimmed_metrics": {
                "paired_job_count": max(pair_count - 1, 0),
                "pearson_r": accepted_r + 0.01,
            },
            "predict_outlier_trim_policy": {"outlier_trim_method": "tukey_iqr"},
            "predict_raw_target_metrics": [],
            "predict_raw_target_outlier_trim_metrics": [],
            "predict_calibrated_target_metrics": [],
            "predict_calibrated_target_outlier_trim_metrics": [],
            "model": {"model": model_name, "groups": {"global": {"family": model_name}}},
        }

    monkeypatch.setattr(
        module,
        "calibrate_ab_bind_plan",
        lambda *args, **kwargs: fake_payload(kwargs["model"]),
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_calibrated_validation.py",
            "--plan-root",
            str(plan_root),
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 0

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["requested_model"] == "auto"
    assert summary["model_selection_mode"] == "auto"
    assert summary["selected_model"] == "side_linear"
    assert summary["model"] == "side_linear"
    assert summary["accepted_calibrated_view"] == "target_filtered"
    assert summary["accepted_calibrated_pearson_r"] == 0.6246238526286293
    assert summary["accepted_calibrated_excluded_complex_ids"] == ["1MLC"]
    assert summary["accepted_calibrated_passed"] is True
    assert summary["calibrated_target_filtered_pearson_r"] == 0.6246238526286293
    assert summary["calibrated_target_filtered_outlier_trimmed_pearson_r"] == 0.6446238526286293
    assert summary["accepted_independent_view"] == "calibrated_target_filtered_outlier_trimmed"
    assert summary["accepted_independent_pearson_r"] == 0.6446238526286293
    assert summary["accepted_independent_passed"] is True
    assert summary["model_leaderboard"][0]["model"] == "side_linear"
    assert summary["model_leaderboard"][0]["accepted_calibrated_passed"] is True
    assert (
        summary["model_leaderboard"][0]["accepted_independent_view"]
        == "calibrated_target_filtered_outlier_trimmed"
    )


def test_main_writes_generated_at_for_insufficient_fit_pairs(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_core_v1_quick_plan"
    reports_dir = plan_root / "reports"
    summary_output = reports_dir / "calibrated_validation_summary.json"
    plan_root.mkdir(parents=True)

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir),
            "paired_job_count": 0,
            "qc_qualified_pair_count": 0,
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_calibrated_validation.py",
            "--plan-root",
            str(plan_root),
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 2

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["status"] == "insufficient_fit_pairs"
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
