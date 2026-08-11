from __future__ import annotations

import importlib.util
import json
from datetime import datetime
from pathlib import Path


SCRIPT_PATH = Path("/mnt/data/liuchao/abag-rbfep/benchmarks/patel_2021_3hfm/report_patellike_3hfm.py")


def _load_script_module():
    spec = importlib.util.spec_from_file_location("report_patellike_3hfm", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_default_summary_output_path_lands_under_batch_reports() -> None:
    module = _load_script_module()
    batch_dir = Path("/tmp/abag-rbfep/runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference")

    assert module.default_summary_output_path(batch_dir) == batch_dir / "reports" / "patel_2021_3hfm_summary.json"


def test_main_writes_insufficient_pairs_summary(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    reports_dir = batch_dir / "reports"
    reports_dir.mkdir(parents=True)
    experimental_csv = tmp_path / "experimental_ddg.csv"
    summary_output = reports_dir / "patel_2021_3hfm_summary.json"
    experimental_csv.write_text(
        "\n".join(
            [
                "job_id,mutation_group_id,chain_id,resseq,wt,mut,entity_side,experimental_ddg_kcal_mol,chain_mapping_basis,source_reference",
                "3hfm-patel-2021-antigen-y-y20f,patel_3hfm_y20f,Y,20,Y,F,antigen,-0.48,explicit,Patel 2021 Table 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        module,
        "write_batch_summary",
        lambda _batch_dir: {
            "batch_dir": str(batch_dir),
            "jobs": [
                {
                    "job_id": "3hfm-patel-2021-antigen-y-y20f",
                    "mutation_group_id": "patel_3hfm_y20f",
                    "ddg_kcal_mol": None,
                    "ddg_ready": False,
                    "benchmark_qc_qualified": False,
                    "ddg_bar_stderr_kcal_mol": None,
                    "qc_status": "not_evaluated",
                    "diagnostic_code": "pending_equilibrate",
                    "latest_stage": "build_legs",
                    "latest_stage_state": "completed",
                }
            ],
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_patellike_3hfm.py",
            "--batch-dir",
            str(batch_dir),
            "--experimental-csv",
            str(experimental_csv),
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 2

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["status"] == "insufficient_pairs"
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["incomplete_job_count"] == 1
    assert summary["paired_job_count"] == 0
    assert summary["charge_class_summary"]["charge_conserving"]["job_count"] == 1
    assert summary["charge_class_summary"]["charge_conserving"]["incomplete_job_count"] == 1
    assert summary["charge_class_summary"]["charge_changing"]["job_count"] == 0
    assert summary["incomplete_diagnostic_code_counts"] == {"pending_equilibrate": 1}
    assert summary["incomplete_diagnostic_code_counts_by_charge_class"] == {
        "charge_conserving": {"pending_equilibrate": 1}
    }
    assert "charge-conserving subset is still in progress" in summary["message"]


def test_main_summarizes_ready_pairs(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    reports_dir = batch_dir / "reports"
    reports_dir.mkdir(parents=True)
    experimental_csv = tmp_path / "experimental_ddg.csv"
    summary_output = reports_dir / "patel_2021_3hfm_summary.json"
    experimental_csv.write_text(
        "\n".join(
            [
                "job_id,mutation_group_id,chain_id,resseq,wt,mut,entity_side,experimental_ddg_kcal_mol,chain_mapping_basis,source_reference",
                "3hfm-patel-2021-antigen-y-y20f,patel_3hfm_y20f,Y,20,Y,F,antigen,-0.48,explicit,Patel 2021 Table 2",
                "3hfm-patel-2021-antibody-h-w98f,patel_3hfm_w98f,H,98,W,F,antibody,3.25,inferred,Patel 2021 Table 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        module,
        "write_batch_summary",
        lambda _batch_dir: {
            "batch_dir": str(batch_dir),
            "jobs": [
                {
                    "job_id": "3hfm-patel-2021-antigen-y-y20f",
                    "mutation_group_id": "patel_3hfm_y20f",
                    "ddg_kcal_mol": -0.34,
                    "ddg_ready": True,
                    "benchmark_qc_qualified": True,
                    "ddg_bar_stderr_kcal_mol": 0.3,
                    "qc_status": "pass",
                    "diagnostic_code": "qc_pass",
                    "latest_stage": "report",
                    "latest_stage_state": "completed",
                },
                {
                    "job_id": "3hfm-patel-2021-antibody-h-w98f",
                    "mutation_group_id": "patel_3hfm_w98f",
                    "ddg_kcal_mol": 2.23,
                    "ddg_ready": True,
                    "benchmark_qc_qualified": False,
                    "ddg_bar_stderr_kcal_mol": 0.5,
                    "qc_status": "warning",
                    "diagnostic_code": "qc_repeat_spread",
                    "latest_stage": "report",
                    "latest_stage_state": "completed",
                },
            ],
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_patellike_3hfm.py",
            "--batch-dir",
            str(batch_dir),
            "--experimental-csv",
            str(experimental_csv),
            "--summary-output",
            str(summary_output),
            "--top-n",
            "1",
        ],
    )

    assert module.main() == 0

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["status"] == "ok"
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["paired_job_count"] == 2
    assert summary["qc_qualified_pair_count"] == 1
    assert summary["raw_metrics"]["pearson_r"] is not None
    assert summary["top_abs_error_pairs"][0]["job_id"] == "3hfm-patel-2021-antibody-h-w98f"
    assert summary["charge_class_summary"]["charge_conserving"]["job_count"] == 2
    assert summary["charge_class_summary"]["charge_conserving"]["paired_job_count"] == 2
    assert summary["charge_class_summary"]["charge_conserving"]["qc_qualified_pair_count"] == 1
    assert summary["charge_class_summary"]["charge_changing"]["job_count"] == 0


def test_main_classifies_charge_changing_incomplete_jobs(tmp_path: Path, monkeypatch) -> None:
    module = _load_script_module()
    batch_dir = tmp_path / "runs" / "benchmarks" / "patel_2021_3hfm" / "patel_2021_3hfm_reference"
    reports_dir = batch_dir / "reports"
    reports_dir.mkdir(parents=True)
    experimental_csv = tmp_path / "experimental_ddg.csv"
    summary_output = reports_dir / "patel_2021_3hfm_summary.json"
    experimental_csv.write_text(
        "\n".join(
            [
                "job_id,mutation_group_id,chain_id,resseq,wt,mut,entity_side,experimental_ddg_kcal_mol,chain_mapping_basis,source_reference",
                "3hfm-patel-2021-antibody-h-d32n,patel_3hfm_d32n,H,32,D,N,antibody,0.17,explicit,Patel 2021 Table 2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        module,
        "write_batch_summary",
        lambda _batch_dir: {
            "batch_dir": str(batch_dir),
            "jobs": [
                {
                    "job_id": "3hfm-patel-2021-antibody-h-d32n",
                    "mutation_group_id": "patel_3hfm_d32n",
                    "ddg_kcal_mol": None,
                    "ddg_ready": False,
                    "benchmark_qc_qualified": False,
                    "ddg_bar_stderr_kcal_mol": None,
                    "qc_status": "not_evaluated",
                    "diagnostic_code": "pending_equilibrate",
                    "latest_stage": "build_legs",
                    "latest_stage_state": "completed",
                }
            ],
        },
    )
    monkeypatch.setattr(
        "sys.argv",
        [
            "report_patellike_3hfm.py",
            "--batch-dir",
            str(batch_dir),
            "--experimental-csv",
            str(experimental_csv),
            "--summary-output",
            str(summary_output),
        ],
    )

    assert module.main() == 2

    summary = json.loads(summary_output.read_text(encoding="utf-8"))
    assert summary["generated_at"].endswith("Z")
    assert datetime.fromisoformat(summary["generated_at"].replace("Z", "+00:00"))
    assert summary["charge_class_summary"]["charge_changing"]["job_count"] == 1
    assert summary["charge_class_summary"]["charge_changing"]["incomplete_job_count"] == 1
    assert summary["incomplete_diagnostic_code_counts_by_charge_class"] == {
        "charge_changing": {"pending_equilibrate": 1}
    }
    assert summary["charge_class_summary"]["charge_conserving"]["job_count"] == 0
    assert "charge-changing Patel rows remain incomplete" in summary["message"]
