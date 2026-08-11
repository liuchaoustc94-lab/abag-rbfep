from __future__ import annotations

import importlib.util
import json
from datetime import datetime
from pathlib import Path


SCRIPT_PATH = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/ab_bind/report_3hfm_protocol_regression.py")


def _load_script_module():
    spec = importlib.util.spec_from_file_location("report_3hfm_protocol_regression", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_default_summary_output_path_lands_under_plan_root_reports() -> None:
    module = _load_script_module()
    plan_root = Path("/tmp/abag-rbfep/runs/benchmarks/abbind_3hfm_protocol_regression")

    assert module.default_summary_output_path(plan_root) == plan_root / "reports" / "3hfm_protocol_regression_summary.json"
    assert (
        module.default_merged_summary_alias_path(plan_root)
        == plan_root / "reports" / "3hfm_protocol_regression_merged_summary.json"
    )


def test_main_writes_insufficient_pairs_summary(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_3hfm_protocol_regression"
    reports_dir = plan_root / "reports" / "selections" / "complex-3hfm"
    summary_output = plan_root / "reports" / "3hfm_protocol_regression_summary.json"
    reports_dir.mkdir(parents=True)
    plan_root.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir),
            "source_plan_roots": [str(plan_root)],
            "selected_job_count": 14,
            "ddg_ready_count": 0,
            "resumable_job_count": 14,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
            "benchmark_metrics": {},
            "benchmark_metrics_qc_qualified": {},
            "validation_failure_taxonomy": {"counts": {"pending_execution": 14}},
            "qc_counts": {"not_started": 14},
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_3hfm_protocol_regression.py",
            "--plan-root",
            str(plan_root),
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 2

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["status"] == "insufficient_pairs"
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["paired_job_count"] == 0
    assert summary["validation_failure_taxonomy"]["counts"]["pending_execution"] == 14


def test_main_summarizes_top_abs_error_pairs(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_3hfm_protocol_regression"
    reports_dir = plan_root / "reports" / "selections" / "complex-3hfm"
    summary_output = plan_root / "reports" / "3hfm_protocol_regression_summary.json"
    reports_dir.mkdir(parents=True)
    plan_root.mkdir(parents=True, exist_ok=True)
    (reports_dir / "benchmark_pairs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,batch_id,source_plan_root,job_id,mutation_group_id,"
                    "complex_delta_g_kcal_mol,apo_delta_g_kcal_mol,predicted_ddg_kcal_mol,"
                    "experimental_ddg_kcal_mol,ddg_error_kcal_mol,abs_error_kcal_mol,"
                    "ddg_bar_stderr_kcal_mol,max_bar_stderr_kcal_mol,qc_status,benchmark_qc_qualified"
                ),
                "3HFM,abbind_3hfm_core_v1,/tmp/root,3hfm-antibody-h-y50a,3hfm_0925,0,0,7.0,8.0,-1.0,1.0,0.5,10.0,warning,False",
                "3HFM,abbind_3hfm_core_v1,/tmp/root,3hfm-antibody-h-c95a,3hfm_0919,0,0,2.0,5.52,-3.52,3.52,0.4,10.0,warning,False",
                "3HFM,abbind_3hfm_core_v1,/tmp/root,3hfm-antigen-y-y20a,3hfm_0940,0,0,4.6,4.87,-0.27,0.27,0.3,10.0,pass,True",
            ]
        ),
        encoding="utf-8",
    )
    (reports_dir / "benchmark_pairs_qc_qualified.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,batch_id,source_plan_root,job_id,mutation_group_id,"
                    "complex_delta_g_kcal_mol,apo_delta_g_kcal_mol,predicted_ddg_kcal_mol,"
                    "experimental_ddg_kcal_mol,ddg_error_kcal_mol,abs_error_kcal_mol,"
                    "ddg_bar_stderr_kcal_mol,max_bar_stderr_kcal_mol,qc_status,benchmark_qc_qualified"
                ),
                "3HFM,abbind_3hfm_core_v1,/tmp/root,3hfm-antigen-y-y20a,3hfm_0940,0,0,4.6,4.87,-0.27,0.27,0.3,10.0,pass,True",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir),
            "source_plan_roots": [str(plan_root)],
            "selected_job_count": 14,
            "ddg_ready_count": 3,
            "resumable_job_count": 11,
            "running_sample_job_count": 1,
            "running_equilibrate_job_count": 0,
            "benchmark_metrics": {
                "pearson_r": 0.65,
                "spearman_rho": 0.5,
                "sign_accuracy": 2 / 3,
                "rmse_kcal_mol": 2.1,
                "mae_kcal_mol": 1.6,
            },
            "benchmark_metrics_qc_qualified": {
                "pearson_r": 1.0,
                "spearman_rho": 1.0,
                "sign_accuracy": 1.0,
                "rmse_kcal_mol": 0.27,
                "mae_kcal_mol": 0.27,
            },
            "validation_failure_taxonomy": {"counts": {"qc_sampling_issue": 2, "benchmark_qc_qualified": 1}},
            "qc_counts": {"warning": 2, "pass": 1},
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_3hfm_protocol_regression.py",
            "--plan-root",
            str(plan_root),
            "--summary-output",
            str(summary_output),
            "--top-n",
            "2",
        ],
    )

    assert module.main() == 0

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["status"] == "ok"
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["paired_job_count"] == 3
    assert summary["qc_qualified_pair_count"] == 1
    assert summary["overall_pearson_r"] == 0.65
    assert [item["job_id"] for item in summary["top_abs_error_pairs"]] == [
        "3hfm-antibody-h-c95a",
        "3hfm-antibody-h-y50a",
    ]
    assert summary["top_abs_error_pairs"][0]["abs_error_kcal_mol"] == 3.52


def test_main_writes_merged_alias_when_extra_roots_are_used(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    plan_root = tmp_path / "runs" / "benchmarks" / "abbind_3hfm_protocol_regression"
    extra_root = tmp_path / "runs" / "benchmarks" / "abbind_3hfm_protocol_regression_extra"
    reports_dir = plan_root / "reports" / "selections" / "complex-3hfm"
    summary_output = plan_root / "reports" / "3hfm_protocol_regression_summary.json"
    merged_alias = plan_root / "reports" / "3hfm_protocol_regression_merged_summary.json"
    reports_dir.mkdir(parents=True)
    extra_root.mkdir(parents=True)
    plan_root.mkdir(parents=True, exist_ok=True)
    (reports_dir / "benchmark_pairs.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,batch_id,source_plan_root,job_id,mutation_group_id,"
                    "complex_delta_g_kcal_mol,apo_delta_g_kcal_mol,predicted_ddg_kcal_mol,"
                    "experimental_ddg_kcal_mol,ddg_error_kcal_mol,abs_error_kcal_mol,"
                    "ddg_bar_stderr_kcal_mol,max_bar_stderr_kcal_mol,qc_status,benchmark_qc_qualified"
                ),
                "3HFM,abbind_3hfm_core_v1,/tmp/root,3hfm-antigen-y-y20a,3hfm_0940,0,0,4.6,4.87,-0.27,0.27,0.3,10.0,pass,True",
            ]
        ),
        encoding="utf-8",
    )
    (reports_dir / "benchmark_pairs_qc_qualified.csv").write_text(
        "\n".join(
            [
                (
                    "complex_id,batch_id,source_plan_root,job_id,mutation_group_id,"
                    "complex_delta_g_kcal_mol,apo_delta_g_kcal_mol,predicted_ddg_kcal_mol,"
                    "experimental_ddg_kcal_mol,ddg_error_kcal_mol,abs_error_kcal_mol,"
                    "ddg_bar_stderr_kcal_mol,max_bar_stderr_kcal_mol,qc_status,benchmark_qc_qualified"
                ),
                "3HFM,abbind_3hfm_core_v1,/tmp/root,3hfm-antigen-y-y20a,3hfm_0940,0,0,4.6,4.87,-0.27,0.27,0.3,10.0,pass,True",
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        module,
        "report_ab_bind_plan",
        lambda *args, **kwargs: {
            "reports_dir": str(reports_dir),
            "source_plan_roots": [str(plan_root), str(extra_root)],
            "selected_job_count": 1,
            "ddg_ready_count": 1,
            "resumable_job_count": 0,
            "running_sample_job_count": 0,
            "running_equilibrate_job_count": 0,
            "benchmark_metrics": {
                "pearson_r": 1.0,
                "spearman_rho": 1.0,
                "sign_accuracy": 1.0,
                "rmse_kcal_mol": 0.27,
                "mae_kcal_mol": 0.27,
            },
            "benchmark_metrics_qc_qualified": {
                "pearson_r": 1.0,
                "spearman_rho": 1.0,
                "sign_accuracy": 1.0,
                "rmse_kcal_mol": 0.27,
                "mae_kcal_mol": 0.27,
            },
            "validation_failure_taxonomy": {"counts": {"benchmark_qc_qualified": 1}},
            "qc_counts": {"pass": 1},
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_3hfm_protocol_regression.py",
            "--plan-root",
            str(plan_root),
            "--extra-plan-root",
            str(extra_root),
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 0

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    merged_summary = json.loads(merged_alias.read_text(encoding="utf-8"))
    assert merged_summary == summary
