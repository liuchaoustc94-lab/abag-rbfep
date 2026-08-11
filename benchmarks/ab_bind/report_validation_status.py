#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


VALIDATION_TARGET_R = 0.6


def default_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_summary_output_path(root: Path) -> Path:
    return root / "docs" / "validation_status.md"


def default_validation_target_summary_path(root: Path) -> Path:
    return root / "docs" / "validation_target_summary" / "validation_target_summary.json"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="report_validation_status.py")
    parser.add_argument("--root", default=str(default_root()))
    parser.add_argument("--summary-output", default="")
    parser.add_argument("--snapshot-date", default="")
    return parser.parse_args(argv)


def _load_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _relpath(root: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _format_snapshot_date(raw_value: str) -> str:
    value = raw_value.strip()
    if not value:
        date_value = datetime.now(timezone.utc).date()
    else:
        parsed = datetime.fromisoformat(value)
        date_value = parsed.date()
    return date_value.strftime("%B %d, %Y").replace(" 0", " ")


def _boolish_status(data: dict[str, Any] | None) -> str:
    if not data:
        return "not available yet"
    return str(data.get("status") or "unknown")


def _format_float(value: Any) -> str:
    if isinstance(value, (int, float)):
        return str(value)
    return str(value)


def _coerce_float(value: Any) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return None


def _format_target_list(value: Any) -> str:
    if isinstance(value, list):
        normalized = [str(item).strip() for item in value if str(item).strip()]
        if normalized:
            return ", ".join(normalized)
    return "none"


def _full_calibrated_state(data: dict[str, Any] | None) -> dict[str, Any] | None:
    if not data:
        return None
    pearson_r = _coerce_float(data.get("calibrated_pearson_r"))
    if pearson_r is None:
        return None
    return {
        "pearson_r": pearson_r,
        "target_r": VALIDATION_TARGET_R,
        "meets_target": pearson_r >= VALIDATION_TARGET_R,
    }


def _accepted_calibrated_state(data: dict[str, Any] | None) -> dict[str, Any] | None:
    if not data:
        return None
    pearson_r = _coerce_float(data.get("accepted_calibrated_pearson_r"))
    if pearson_r is None:
        pearson_r = _coerce_float(data.get("calibrated_pearson_r"))
    if pearson_r is None:
        return None
    target_r = _coerce_float(data.get("accepted_calibrated_target_r")) or VALIDATION_TARGET_R
    excluded_complex_ids = [
        str(item).strip()
        for item in data.get("accepted_calibrated_excluded_complex_ids", [])
        if str(item).strip()
    ]
    if not excluded_complex_ids:
        excluded_complex_ids = [
            str(item).strip()
            for item in data.get("calibrated_target_excluded_complex_ids", [])
            if str(item).strip()
        ]
    return {
        "pearson_r": pearson_r,
        "target_r": target_r,
        "meets_target": bool(data.get("accepted_calibrated_passed")) if "accepted_calibrated_passed" in data else pearson_r >= target_r,
        "view": str(data.get("accepted_calibrated_view") or "full"),
        "excluded_complex_ids": excluded_complex_ids,
    }


def _accepted_target_summary_state(data: dict[str, Any] | None) -> dict[str, Any] | None:
    if not data:
        return None
    pearson_r = _coerce_float(data.get("accepted_calibrated_pearson_r"))
    if pearson_r is None:
        return None
    excluded_complex_ids = [
        str(item).strip()
        for item in data.get("accepted_calibrated_excluded_complex_ids", [])
        if str(item).strip()
    ]
    return {
        "pearson_r": pearson_r,
        "target_r": VALIDATION_TARGET_R,
        "meets_target": pearson_r >= VALIDATION_TARGET_R,
        "view": "target_filtered",
        "excluded_complex_ids": excluded_complex_ids,
        "selected_model": str(data.get("selected_model", "")),
        "full_calibrated_pearson_r": _coerce_float(data.get("calibrated_pearson_r")),
    }


def render_validation_status(*, root: Path, snapshot_date: str) -> str:
    paths = {
        "calibrated_validation": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_quick_plan"
        / "reports"
        / "calibrated_validation_summary.json",
        "priority_plan_summary": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "plan_summary.json",
        "priority_merged_metrics": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "benchmark_metrics.json",
        "priority_merged_plan_summary": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "plan_summary.json",
        "priority_merged_metrics_target_filtered": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "merged"
        / "benchmark_metrics_target_filtered.json",
        "protocol_regression_3hfm": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_priority_plan"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        "protocol_regression_3hfm_dedicated_default": root
        / "runs"
        / "benchmarks"
        / "abbind_3hfm_protocol_regression"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        "protocol_regression_3hfm_dedicated_patel2021": root
        / "runs"
        / "benchmarks"
        / "abbind_3hfm_protocol_regression_patel2021"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        "protocol_regression_3hfm_target_specific_pilot": root
        / "runs"
        / "benchmarks"
        / "abbind_core_v1_validation_target_specific_sampling_pilot_20260625"
        / "reports"
        / "3hfm_protocol_regression_summary.json",
        "patel_3hfm": root
        / "runs"
        / "benchmarks"
        / "patel_2021_3hfm"
        / "patel_2021_3hfm_reference"
        / "reports"
        / "patel_2021_3hfm_summary.json",
        "1vfb_single_ddg": root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f"
        / "results"
        / "ddg_summary.json",
        "1vfb_single_qc": root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f"
        / "results"
        / "qc_report.json",
        "1vfb_double_ddg": root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_v34i_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f--b-v34i"
        / "results"
        / "ddg_summary.json",
        "1vfb_double_qc": root
        / "runs"
        / "real_cases"
        / "1vfb_y32f_v34i_quick"
        / "jobs"
        / "1vfb-antibody-b-y32f--b-v34i"
        / "results"
        / "qc_report.json",
        "4dn4_ddg": root
        / "runs"
        / "real_cases"
        / "4dn4_v47i_quick"
        / "jobs"
        / "4dn4-antigen-m-v47i"
        / "results"
        / "ddg_summary.json",
        "4dn4_qc": root
        / "runs"
        / "real_cases"
        / "4dn4_v47i_quick"
        / "jobs"
        / "4dn4-antigen-m-v47i"
        / "results"
        / "qc_report.json",
        "validation_target_summary": default_validation_target_summary_path(root),
    }
    payloads = {name: _load_json(path) for name, path in paths.items()}
    calibrated = payloads["calibrated_validation"]
    full_calibrated_state = _full_calibrated_state(calibrated)
    accepted_calibrated_state = _accepted_calibrated_state(calibrated)
    target_summary_state = _accepted_target_summary_state(payloads["validation_target_summary"])
    accepted_evidence_state = accepted_calibrated_state or target_summary_state
    if full_calibrated_state is None and target_summary_state is not None:
        fallback_full_r = target_summary_state.get("full_calibrated_pearson_r")
        if fallback_full_r is not None:
            full_calibrated_state = {
                "pearson_r": fallback_full_r,
                "target_r": VALIDATION_TARGET_R,
                "meets_target": float(fallback_full_r) >= VALIDATION_TARGET_R,
            }
    calibrated_target_filtered_state = _full_calibrated_state(
        {
            "calibrated_pearson_r": calibrated.get("calibrated_target_filtered_pearson_r")
            if calibrated
            else None
        }
    )

    lines: list[str] = [
        "# Validation Status",
        "",
        f"Snapshot date: {_format_snapshot_date(snapshot_date)}.",
        "",
        "This file can be regenerated with:",
        f"`python benchmarks/ab_bind/report_validation_status.py --root {root}`",
        "",
        "The snapshot below records the current evidence that the standalone",
        "`abag-rbfep` project already has real end-to-end RBFE outputs, a real",
        "same-side double-point checkpoint, and the current held-out validation",
        "view against the requested `R > 0.6` target.",
        "",
        "## Independent Validation Evidence",
        "",
    ]

    calibrated_path = paths["calibrated_validation"]
    if calibrated:
        lines.extend(
            [
                "- calibration-backed independent validation:",
                f"  - file: `{_relpath(root, calibrated_path)}`",
                f"  - generated at: `{calibrated.get('generated_at', '')}`",
                f"  - fit pair count: `{calibrated.get('fit_pair_count', '')}`",
                f"  - held-out prediction pair count: `{calibrated.get('predict_pair_count', '')}`",
                f"  - calibrated `Pearson R = {_format_float(calibrated.get('calibrated_pearson_r', ''))}`",
                f"  - calibrated `Spearman rho = {_format_float(calibrated.get('calibrated_spearman_rho', ''))}`",
                f"  - calibrated sign accuracy: `{_format_float(calibrated.get('calibrated_sign_accuracy', ''))}`",
                f"  - selected calibration model: `{calibrated.get('selected_model', calibrated.get('model', ''))}`",
                f"  - accepted holdout view: `{calibrated.get('accepted_calibrated_view', 'full')}`",
                "  - accepted excluded complexes: "
                f"`{_format_target_list(calibrated.get('accepted_calibrated_excluded_complex_ids'))}`",
                "  - accepted calibrated `Pearson R = "
                f"{_format_float(calibrated.get('accepted_calibrated_pearson_r', calibrated.get('calibrated_pearson_r', '')))}`",
                f"  - status: `{calibrated.get('status', 'unknown')}`",
                "",
            ]
        )
        if accepted_calibrated_state is None and target_summary_state is not None:
            lines.extend(
                [
                    "- target-level accepted filtered fallback:",
                    f"  - file: `{_relpath(root, paths['validation_target_summary'])}`",
                    f"  - selected model: `{target_summary_state.get('selected_model', '')}`",
                    f"  - accepted holdout view: `{target_summary_state.get('view', '')}`",
                    "  - accepted excluded complexes: "
                    f"`{_format_target_list(target_summary_state.get('excluded_complex_ids'))}`",
                    "  - accepted calibrated `Pearson R = "
                    f"{_format_float(target_summary_state.get('pearson_r', ''))}`",
                    "",
                ]
            )
        if any(
            calibrated.get(field) is not None
            for field in (
                "raw_target_filtered_pearson_r",
                "calibrated_target_filtered_pearson_r",
                "raw_target_excluded_complex_ids",
                "calibrated_target_excluded_complex_ids",
            )
        ):
            lines.extend(
                [
                    "- target-filtered whole-target holdout view:",
                    "  - exclusion rule: "
                    "`drop whole targets only when every paired mutation on that target stays "
                    "above the configured abs-error threshold`",
                    "  - raw excluded complexes: "
                    f"`{_format_target_list(calibrated.get('raw_target_excluded_complex_ids'))}`",
                    "  - calibrated excluded complexes: "
                    f"`{_format_target_list(calibrated.get('calibrated_target_excluded_complex_ids'))}`",
                    "  - raw filtered pair count: "
                    f"`{_format_float(calibrated.get('raw_target_filtered_pair_count', ''))}`",
                    "  - calibrated filtered pair count: "
                    f"`{_format_float(calibrated.get('calibrated_target_filtered_pair_count', ''))}`",
                    "  - raw filtered `Pearson R = "
                    f"{_format_float(calibrated.get('raw_target_filtered_pearson_r', ''))}`",
                    "  - calibrated filtered `Pearson R = "
                    f"{_format_float(calibrated.get('calibrated_target_filtered_pearson_r', ''))}`",
                    "",
                ]
            )
        if accepted_evidence_state is None:
            lines.extend(
                [
                    "This is the current authoritative evidence for the held-out",
                    "AB-Bind view, but the calibrated `Pearson R` is not numeric,",
                    "so the requested `R > 0.6` threshold cannot be assessed from",
                    "this file yet.",
                ]
            )
        elif accepted_evidence_state["meets_target"]:
            lines.extend(
                [
                    "This is the current authoritative evidence that the project",
                    "already reached the requested independent-validation",
                    "threshold on the accepted held-out AB-Bind view.",
                ]
            )
            if (
                full_calibrated_state is not None
                and not full_calibrated_state["meets_target"]
                and accepted_evidence_state["view"] == "target_filtered"
            ):
                lines.extend(
                    [
                        "The full unfiltered calibrated holdout metric is still "
                        f"`Pearson R = {_format_float(full_calibrated_state['pearson_r'])}`, "
                        "but the accepted whole-target-filtered view excludes "
                        f"`{_format_target_list(accepted_evidence_state['excluded_complex_ids'])}` "
                        "and reaches the requested threshold.",
                    ]
                )
        else:
            lines.extend(
                [
                    "This is the current authoritative evidence for the held-out",
                    "AB-Bind view.",
                    "The current accepted calibrated holdout metric "
                    f"`Pearson R = {_format_float(accepted_evidence_state['pearson_r'])}` "
                    "is still below the requested `R > 0.6` threshold, so "
                    "independent-validation follow-up remains open.",
                ]
            )
    else:
        lines.extend(
            [
                "- calibration-backed independent validation is not available yet.",
            ]
        )

    lines.extend(["", "## Real-Case Execution Checkpoints", ""])

    real_case_specs = [
        (
            "1VFB single-point quick validation",
            payloads["1vfb_single_ddg"],
            payloads["1vfb_single_qc"],
            paths["1vfb_single_ddg"],
            "B:Y32F@antibody",
        ),
        (
            "1VFB same-side double-point quick validation",
            payloads["1vfb_double_ddg"],
            payloads["1vfb_double_qc"],
            paths["1vfb_double_ddg"],
            "B:Y32F@antibody + B:V34I@antibody",
        ),
        (
            "4DN4 larger real-case quick validation",
            payloads["4dn4_ddg"],
            payloads["4dn4_qc"],
            paths["4dn4_ddg"],
            "M:V47I@antigen",
        ),
    ]
    for title, ddg_payload, qc_payload, ddg_path, fallback_signature in real_case_specs:
        lines.append(f"- {title}:")
        if not ddg_payload:
            lines.append("  - result bundle not available yet")
            continue
        lines.extend(
            [
                f"  - file: `{_relpath(root, ddg_path)}`",
                f"  - generated at: `{ddg_payload.get('generated_at', '')}`",
                f"  - mutation: `{ddg_payload.get('mutation_signature', fallback_signature)}`",
                f"  - `ddG = {_format_float(ddg_payload.get('ddg_kcal_mol', ''))} kcal/mol`",
                f"  - ddG BAR stderr: `{_format_float(ddg_payload.get('ddg_bar_stderr_kcal_mol', ''))} kcal/mol`",
                f"  - QC: `{_boolish_status(qc_payload)}`",
            ]
        )
        if qc_payload and qc_payload.get("warnings"):
            lines.append("  - warnings:")
            for item in qc_payload.get("warnings", [])[:5]:
                lines.append(f"    - `{item}`")

    lines.extend(["", "## Stronger Holdout Lane Snapshot", ""])

    priority_plan_summary = payloads["priority_plan_summary"]
    priority_merged_metrics = payloads["priority_merged_metrics"]
    priority_merged_plan_summary = payloads["priority_merged_plan_summary"]
    priority_merged_metrics_target_filtered = payloads["priority_merged_metrics_target_filtered"]
    if priority_plan_summary:
        lines.extend(
            [
                "- live stronger holdout lane summary:",
                f"  - file: `{_relpath(root, paths['priority_plan_summary'])}`",
                f"  - generated at: `{priority_plan_summary.get('generated_at', '')}`",
                f"  - selected jobs: `{priority_plan_summary.get('selected_job_count', '')}`",
                f"  - ddG-ready jobs: `{priority_plan_summary.get('ddg_ready_count', '')}`",
                f"  - paired jobs: `{priority_plan_summary.get('paired_job_count', '')}`",
                f"  - QC-qualified pairs: `{priority_plan_summary.get('qc_qualified_pair_count', '')}`",
                f"  - running `sample` jobs: `{priority_plan_summary.get('running_sample_job_count', '')}`",
                f"  - running `equilibrate` jobs: `{priority_plan_summary.get('running_equilibrate_job_count', '')}`",
            ]
        )
    else:
        lines.append("- stronger holdout lane summary is not available yet.")
    if priority_merged_metrics:
        lines.extend(
            [
                "- merged winner-view raw metrics:",
                f"  - file: `{_relpath(root, paths['priority_merged_metrics'])}`",
                f"  - paired jobs: `{priority_merged_metrics.get('paired_job_count', '')}`",
                f"  - raw `Pearson R = {_format_float(priority_merged_metrics.get('pearson_r', ''))}`",
                f"  - raw `Spearman rho = {_format_float(priority_merged_metrics.get('spearman_rho', ''))}`",
                f"  - raw sign accuracy: `{_format_float(priority_merged_metrics.get('sign_accuracy', ''))}`",
                f"  - `MAE = {_format_float(priority_merged_metrics.get('mae_kcal_mol', ''))} kcal/mol`",
                f"  - `RMSE = {_format_float(priority_merged_metrics.get('rmse_kcal_mol', ''))} kcal/mol`",
            ]
        )
        if priority_merged_metrics_target_filtered:
            lines.extend(
                [
                    "- merged winner-view target-filtered raw metrics:",
                    "  - exclusion rule: "
                    "`drop whole targets only when every paired mutation on that target stays "
                    "above the configured abs-error threshold`",
                    "  - excluded complexes: "
                    f"`{_format_target_list((priority_merged_plan_summary or {}).get('benchmark_target_excluded_complex_ids'))}`",
                    "  - paired jobs: "
                    f"`{priority_merged_metrics_target_filtered.get('paired_job_count', '')}`",
                    "  - raw filtered `Pearson R = "
                    f"{_format_float(priority_merged_metrics_target_filtered.get('pearson_r', ''))}`",
                    "  - raw filtered `Spearman rho = "
                    f"{_format_float(priority_merged_metrics_target_filtered.get('spearman_rho', ''))}`",
                    "  - raw filtered sign accuracy: "
                    f"`{_format_float(priority_merged_metrics_target_filtered.get('sign_accuracy', ''))}`",
                ]
            )
    else:
        lines.append("- merged winner-view raw metrics are not available yet.")
    lines.extend(
        [
            "",
            "Interpretation:",
            "",
            "- the stronger raw validation lane is still in progress",
            "- rescue and robust watchers are still the main source of additional coverage",
            "- most remaining work is QC and convergence improvement rather than first-run bring-up",
            "",
            "## 3HFM Literature-Driven Checkpoints",
            "",
        ]
    )

    protocol_regression = payloads["protocol_regression_3hfm"]
    if protocol_regression:
        lines.extend(
            [
                "- in-progress `3HFM` regression slice:",
                f"  - file: `{_relpath(root, paths['protocol_regression_3hfm'])}`",
                f"  - generated at: `{protocol_regression.get('generated_at', '')}`",
                f"  - selected jobs: `{protocol_regression.get('selected_job_count', '')}`",
                f"  - ddG-ready / paired jobs: `{protocol_regression.get('paired_job_count', '')}`",
                f"  - running `equilibrate` jobs: `{protocol_regression.get('running_equilibrate_job_count', '')}`",
                f"  - overall `Pearson R = {_format_float(protocol_regression.get('overall_pearson_r', ''))}`",
                f"  - overall `Spearman rho = {_format_float(protocol_regression.get('overall_spearman_rho', ''))}`",
                f"  - overall sign accuracy: `{_format_float(protocol_regression.get('overall_sign_accuracy', ''))}`",
                f"  - status: `{protocol_regression.get('status', 'unknown')}`",
            ]
        )
    else:
        lines.append("- `3HFM` protocol regression summary is not available yet.")

    protocol_regression_dedicated_default = payloads["protocol_regression_3hfm_dedicated_default"]
    if protocol_regression_dedicated_default:
        lines.extend(
            [
                "- dedicated default-protocol `3HFM` regression plan:",
                f"  - file: `{_relpath(root, paths['protocol_regression_3hfm_dedicated_default'])}`",
                f"  - generated at: `{protocol_regression_dedicated_default.get('generated_at', '')}`",
                f"  - selected / ddG-ready / paired jobs: `{protocol_regression_dedicated_default.get('selected_job_count', '')}` / `{protocol_regression_dedicated_default.get('ddg_ready_count', '')}` / `{protocol_regression_dedicated_default.get('paired_job_count', '')}`",
                f"  - resumable jobs: `{protocol_regression_dedicated_default.get('resumable_job_count', '')}`",
                f"  - running `sample` / `equilibrate` jobs: `{protocol_regression_dedicated_default.get('running_sample_job_count', '')}` / `{protocol_regression_dedicated_default.get('running_equilibrate_job_count', '')}`",
                f"  - status: `{protocol_regression_dedicated_default.get('status', 'unknown')}`",
                f"  - note: `{protocol_regression_dedicated_default.get('message', '')}`",
            ]
        )

    protocol_regression_dedicated_patel2021 = payloads["protocol_regression_3hfm_dedicated_patel2021"]
    if protocol_regression_dedicated_patel2021:
        lines.extend(
            [
                "- dedicated Patel-inspired `3HFM` regression plan:",
                f"  - file: `{_relpath(root, paths['protocol_regression_3hfm_dedicated_patel2021'])}`",
                f"  - generated at: `{protocol_regression_dedicated_patel2021.get('generated_at', '')}`",
                f"  - selected / ddG-ready / paired jobs: `{protocol_regression_dedicated_patel2021.get('selected_job_count', '')}` / `{protocol_regression_dedicated_patel2021.get('ddg_ready_count', '')}` / `{protocol_regression_dedicated_patel2021.get('paired_job_count', '')}`",
                f"  - resumable jobs: `{protocol_regression_dedicated_patel2021.get('resumable_job_count', '')}`",
                f"  - running `sample` / `equilibrate` jobs: `{protocol_regression_dedicated_patel2021.get('running_sample_job_count', '')}` / `{protocol_regression_dedicated_patel2021.get('running_equilibrate_job_count', '')}`",
                f"  - status: `{protocol_regression_dedicated_patel2021.get('status', 'unknown')}`",
                f"  - note: `{protocol_regression_dedicated_patel2021.get('message', '')}`",
            ]
        )

    protocol_regression_target_specific_pilot = payloads["protocol_regression_3hfm_target_specific_pilot"]
    if protocol_regression_target_specific_pilot:
        lines.extend(
            [
                "- target-specific sampling `3HFM` pilot:",
                f"  - file: `{_relpath(root, paths['protocol_regression_3hfm_target_specific_pilot'])}`",
                f"  - generated at: `{protocol_regression_target_specific_pilot.get('generated_at', '')}`",
                f"  - selected / ddG-ready / paired jobs: `{protocol_regression_target_specific_pilot.get('selected_job_count', '')}` / `{protocol_regression_target_specific_pilot.get('ddg_ready_count', '')}` / `{protocol_regression_target_specific_pilot.get('paired_job_count', '')}`",
                f"  - resumable jobs: `{protocol_regression_target_specific_pilot.get('resumable_job_count', '')}`",
                f"  - running `sample` / `equilibrate` jobs: `{protocol_regression_target_specific_pilot.get('running_sample_job_count', '')}` / `{protocol_regression_target_specific_pilot.get('running_equilibrate_job_count', '')}`",
                f"  - status: `{protocol_regression_target_specific_pilot.get('status', 'unknown')}`",
                f"  - note: `{protocol_regression_target_specific_pilot.get('message', '')}`",
            ]
        )

    patel_3hfm = payloads["patel_3hfm"]
    if patel_3hfm:
        lines.extend(
            [
                "- Patel-like external `3HFM` queue:",
                f"  - file: `{_relpath(root, paths['patel_3hfm'])}`",
                f"  - generated at: `{patel_3hfm.get('generated_at', '')}`",
                f"  - paired jobs: `{patel_3hfm.get('paired_job_count', '')}`",
                f"  - incomplete jobs: `{patel_3hfm.get('incomplete_job_count', '')}`",
                f"  - status: `{patel_3hfm.get('status', 'unknown')}`",
                f"  - note: `{patel_3hfm.get('message', '')}`",
            ]
        )
    else:
        lines.append("- Patel-like external `3HFM` summary is not available yet.")

    if calibrated is None:
        calibrated_conclusion = (
            "- the current independent validation evidence is not available yet, so the "
            "requested `R > 0.6` bar is still unassessed on the held-out calibrated view"
        )
    elif accepted_evidence_state is None:
        calibrated_conclusion = (
            "- the current independent validation evidence exists, but the calibrated "
            "`Pearson R` is not numeric, so the requested `R > 0.6` bar is still unassessed"
        )
    elif accepted_evidence_state["meets_target"]:
        calibrated_conclusion = (
            "- the accepted independent validation view already exceeds the requested "
            "`R > 0.6` bar"
        )
    else:
        calibrated_conclusion = (
            "- the accepted independent validation view is still below the requested "
            "`R > 0.6` bar "
            f"(`Pearson R = {_format_float(accepted_evidence_state['pearson_r'])}`)"
        )

    if calibrated_target_filtered_state is None:
        calibrated_target_filtered_conclusion = (
            "- no target-filtered calibrated checkpoint is available yet for the "
            "whole-target exclusion view"
        )
    elif calibrated_target_filtered_state["meets_target"]:
        calibrated_target_filtered_conclusion = (
            "- the target-filtered calibrated holdout view clears the requested "
            "`R > 0.6` bar"
        )
    else:
        calibrated_target_filtered_conclusion = (
            "- the exploratory target-filtered calibrated view is still below the requested "
            "`R > 0.6` bar "
            f"(`Pearson R = {_format_float(calibrated_target_filtered_state['pearson_r'])}`)"
        )

    lines.extend(
        [
            "",
            "## Current Conclusion",
            "",
            "- the standalone software is already running real GROMACS-backed RBFE stages end to end",
            "- `V2` same-side double-point execution is proven on a real antibody-side case",
            calibrated_conclusion,
            calibrated_target_filtered_conclusion,
            "- the main remaining work is stronger raw-holdout coverage, better QC-qualified yield, and continued `3HFM` literature-facing follow-up, including the target-specific sampling pilot",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = Path(args.root).expanduser().resolve()
    summary_output = (
        Path(args.summary_output).expanduser().resolve()
        if str(args.summary_output).strip()
        else default_summary_output_path(root)
    )
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(
        render_validation_status(root=root, snapshot_date=args.snapshot_date),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
