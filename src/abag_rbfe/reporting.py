"""Result extraction and report generation."""

from __future__ import annotations

from functools import lru_cache
from math import sqrt
from pathlib import Path
import re
from statistics import mean, stdev
import subprocess
from typing import Any

from abag_rbfe.constants import STAGES
from abag_rbfe.gmx import inspect_gro_file
from abag_rbfe.io_utils import ensure_dir, read_json, utc_now, write_csv_rows, write_json, write_yaml
from abag_rbfe.structure import partition_inter_residue_sidechain_repairable_clashes

KCAL_PER_MOL_PER_K = 0.00198720425864083
HISTOGRAM_LEGEND_RE = re.compile(r'^@ s(?P<index>\d+) legend "(?P<legend>.*)"$')
HISTOGRAM_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d+|\d+)(?:[eE][-+]?\d+)?")
ACTIVE_STAGE_PATH_MARKERS = ("/artifacts/commands/", "/legs/")
SAMPLE_PHASE_ARTIFACTS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (
        "md",
        (
            "topol.tpr",
            "md.cpt",
            "md.edr",
            "md.gro",
            "md.log",
            "md.trr",
            "md.xtc",
            "dhdl.xvg",
        ),
    ),
    (
        "pre_md",
        (
            "pre_md.tpr",
            "pre_md.cpt",
            "pre_md.edr",
            "pre_md.gro",
            "pre_md.log",
            "pre_md.trr",
            "pre_md.xtc",
            "pre_md.xvg",
        ),
    ),
    (
        "pre_relax",
        (
            "pre_relax.tpr",
            "pre_relax.edr",
            "pre_relax.gro",
            "pre_relax.log",
            "pre_relax.trr",
        ),
    ),
)


@lru_cache(maxsize=1)
def _active_process_lines() -> tuple[str, ...]:
    try:
        result = subprocess.run(["ps", "-ef"], check=False, capture_output=True, text=True)
    except OSError:
        return ()
    if result.returncode != 0:
        return ()
    return tuple(result.stdout.splitlines())


def _job_has_live_process(job_dir: Path) -> bool:
    job_dir_str = str(job_dir)
    batch_dir_str = str(job_dir.parent.parent)
    job_id_pattern = re.compile(rf"\b{re.escape(job_dir.name)}\b")
    for line in _active_process_lines():
        if job_dir_str in line and any(marker in line for marker in ACTIVE_STAGE_PATH_MARKERS):
            return True
        if "abag-rbfe" in line and batch_dir_str in line and job_id_pattern.search(line):
            return True
    return False


def _normalize_running_stage_records(job_dir: Path, stage_records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not stage_records:
        return stage_records
    latest_stage = stage_records[-1]
    if latest_stage.get("state") != "running":
        return stage_records
    if not latest_stage.get("commands") and not latest_stage.get("artifacts"):
        return stage_records
    if _job_has_live_process(job_dir):
        return stage_records

    normalized = list(stage_records[:-1])
    stale_stage = dict(latest_stage)
    stale_stage["reported_state"] = latest_stage.get("state")
    stale_stage["state"] = "stale_running"
    stale_stage["stale_reason"] = "no_active_process_detected"
    normalized.append(stale_stage)
    return normalized


def _read_stage_records(job_dir: Path) -> list[dict[str, Any]]:
    stage_records = []
    for stage in STAGES:
        stage_file = job_dir / "stages" / f"{stage}.json"
        if stage_file.exists():
            stage_records.append(read_json(stage_file))
    return _normalize_running_stage_records(job_dir, stage_records)


def _load_job_spec(job_dir: Path) -> dict[str, Any]:
    spec_path = job_dir / "job_spec.json"
    if not spec_path.exists():
        return {}
    return read_json(spec_path)


def _nonempty_file(path: Path) -> bool:
    try:
        return path.is_file() and path.stat().st_size > 0
    except OSError:
        return False


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _fraction(numerator: int, denominator: int) -> float | None:
    if denominator <= 0:
        return None
    return numerator / denominator


def _sample_window_started(lambda_dir: Path) -> bool:
    if not lambda_dir.is_dir():
        return False
    for child in lambda_dir.iterdir():
        if child.name in {"pre_relax.mdp", "pre_md.mdp", "production.mdp"}:
            continue
        return True
    return False


def _sample_window_completed(lambda_dir: Path) -> bool:
    return _nonempty_file(lambda_dir / "dhdl.xvg") and _nonempty_file(lambda_dir / "md.gro")


def _sample_window_phase(lambda_dir: Path) -> str | None:
    if not lambda_dir.is_dir():
        return None
    if _sample_window_completed(lambda_dir):
        return "completed"
    for phase, artifact_names in SAMPLE_PHASE_ARTIFACTS:
        if any(_nonempty_file(lambda_dir / artifact_name) for artifact_name in artifact_names):
            return phase
    if _sample_window_started(lambda_dir):
        return "started"
    return None


def _sample_progress(job_dir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    protocol = spec.get("protocol", {})
    repeats = max(_safe_int(protocol.get("repeats"), 0), 0)
    lambda_windows = max(_safe_int(protocol.get("lambda_windows"), 0), 0)
    total_windows = repeats * lambda_windows * 2
    started_windows = 0
    completed_windows = 0
    active_window: dict[str, Any] | None = None

    for leg_name in ("complex", "apo"):
        for repeat_index in range(1, repeats + 1):
            repeat_id = f"rep{repeat_index:02d}"
            repeat_dir = job_dir / "legs" / leg_name / repeat_id
            for window_index in range(lambda_windows):
                lambda_id = f"lambda_{window_index:03d}"
                lambda_dir = repeat_dir / lambda_id
                phase = _sample_window_phase(lambda_dir)
                if phase == "completed":
                    started_windows += 1
                    completed_windows += 1
                    continue
                if phase is not None:
                    started_windows += 1
                    if active_window is None:
                        active_window = {
                            "leg": leg_name,
                            "repeat_id": repeat_id,
                            "repeat_index": repeat_index,
                            "lambda_id": lambda_id,
                            "lambda_index": window_index,
                            "phase": phase,
                            "window": f"{leg_name}/{repeat_id}/{lambda_id}",
                        }

    return {
        "completed_windows": completed_windows,
        "started_windows": started_windows,
        "remaining_windows": max(total_windows - completed_windows, 0),
        "total_windows": total_windows,
        "fraction_complete": _fraction(completed_windows, total_windows),
        "active_leg": active_window["leg"] if active_window is not None else None,
        "active_repeat_id": active_window["repeat_id"] if active_window is not None else None,
        "active_repeat_index": active_window["repeat_index"] if active_window is not None else None,
        "active_lambda_id": active_window["lambda_id"] if active_window is not None else None,
        "active_lambda_index": active_window["lambda_index"] if active_window is not None else None,
        "active_phase": active_window["phase"] if active_window is not None else None,
        "active_window": active_window["window"] if active_window is not None else None,
    }


def _equilibrate_repeat_started(repeat_dir: Path) -> bool:
    if _nonempty_file(repeat_dir / "system.top"):
        return True
    for directory in (repeat_dir / "setup", repeat_dir / "equilibration"):
        if not directory.is_dir():
            continue
        for child in directory.iterdir():
            if child.is_file():
                return True
    return False


def _equilibrate_repeat_completed(repeat_dir: Path) -> bool:
    return _nonempty_file(repeat_dir / "equilibration" / "npt.gro")


def _equilibrate_progress(job_dir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    protocol = spec.get("protocol", {})
    repeats = max(_safe_int(protocol.get("repeats"), 0), 0)
    total_repeats = repeats * 2
    started_repeats = 0
    completed_repeats = 0

    for leg_name in ("complex", "apo"):
        for repeat_index in range(1, repeats + 1):
            repeat_dir = job_dir / "legs" / leg_name / f"rep{repeat_index:02d}"
            if _equilibrate_repeat_completed(repeat_dir):
                started_repeats += 1
                completed_repeats += 1
            elif _equilibrate_repeat_started(repeat_dir):
                started_repeats += 1

    return {
        "completed_repeats": completed_repeats,
        "started_repeats": started_repeats,
        "remaining_repeats": max(total_repeats - completed_repeats, 0),
        "total_repeats": total_repeats,
        "fraction_complete": _fraction(completed_repeats, total_repeats),
    }


def _numeric_lines(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped[0] in {"#", "@"}:
            continue
        try:
            rows.append([float(item) for item in stripped.split()])
        except ValueError:
            continue
    return rows


def _safe_mean(values: list[float]) -> float | None:
    if not values:
        return None
    return mean(values)


def _safe_stdev(values: list[float]) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    return stdev(values)


def _safe_range(values: list[float]) -> float | None:
    if not values:
        return None
    return max(values) - min(values)


def _kt_to_kcal(delta_g_kt: float | None, temperature_k: float) -> float | None:
    if delta_g_kt is None:
        return None
    return delta_g_kt * KCAL_PER_MOL_PER_K * temperature_k


def _histogram_legend_metadata(legend: str) -> dict[str, Any]:
    numeric_tokens = [float(token) for token in HISTOGRAM_FLOAT_RE.findall(legend)]
    current_lambda = numeric_tokens[-1] if numeric_tokens else None
    target_lambda = numeric_tokens[0] if len(numeric_tokens) >= 2 else None
    return {
        "legend": legend,
        "kind": "dhdl" if "dH/d" in legend else "delta_h",
        "current_lambda": current_lambda,
        "target_lambda": target_lambda,
    }


def _series_payload(
    *,
    index: int,
    block: list[tuple[float, float]],
    legend: str,
) -> dict[str, Any]:
    return {
        "series_index": index,
        "series": block,
        **_histogram_legend_metadata(legend),
    }


def _parse_histogram_series(path: Path) -> list[dict[str, Any]] | None:
    if not path.exists():
        return None

    legends: dict[int, str] = {}
    blocks: list[list[tuple[float, float]]] = []
    current_block: list[tuple[float, float]] = []
    saw_separator = False

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        legend_match = HISTOGRAM_LEGEND_RE.match(stripped)
        if legend_match:
            legends[int(legend_match.group("index"))] = legend_match.group("legend")
            continue

        if stripped == "@":
            saw_separator = True
            if current_block:
                blocks.append(current_block)
                current_block = []
            continue

        if stripped.startswith("@"):
            continue

        try:
            x_value, y_value = [float(item) for item in stripped.split()[:2]]
        except ValueError:
            continue
        current_block.append((x_value, y_value))

    if current_block:
        blocks.append(current_block)

    if saw_separator and blocks:
        return [
            _series_payload(index=index, block=block, legend=legends.get(index, ""))
            for index, block in enumerate(blocks)
        ]

    rows = _numeric_lines(path)
    if not rows:
        return None

    block_count = len(legends) or 6
    if block_count <= 0 or len(rows) % block_count != 0:
        return None

    block_size = len(rows) // block_count
    if block_size == 0:
        return None

    return [
        _series_payload(
            index=index,
            block=[(row[0], row[1]) for row in rows[index * block_size : (index + 1) * block_size]],
            legend=legends.get(index, ""),
        )
        for index in range(block_count)
    ]


def _lambda_key(value: float | None) -> float | None:
    if value is None:
        return None
    return round(value, 6)


def _histogram_overlap_payload(path: Path) -> dict[str, Any] | None:
    histogram_series = _parse_histogram_series(path)
    if not histogram_series:
        return None

    lambda_order = [
        lambda_value
        for lambda_value in (_lambda_key(series["current_lambda"]) for series in histogram_series if series["kind"] == "dhdl")
        if lambda_value is not None
    ]
    delta_h_series = {
        (_lambda_key(series["current_lambda"]), _lambda_key(series["target_lambda"])): series["series"]
        for series in histogram_series
        if series["kind"] == "delta_h"
        and _lambda_key(series["current_lambda"]) is not None
        and _lambda_key(series["target_lambda"]) is not None
    }

    interval_payloads = []
    for lower_lambda, upper_lambda in zip(lambda_order, lambda_order[1:], strict=False):
        forward = delta_h_series.get((lower_lambda, upper_lambda))
        reverse = delta_h_series.get((upper_lambda, lower_lambda))
        if forward is None or reverse is None:
            continue
        overlap = _support_overlap_fraction(forward, reverse)
        if overlap is None:
            continue
        interval_payloads.append(
            {
                "from_lambda": lower_lambda,
                "to_lambda": upper_lambda,
                **overlap,
            }
        )

    if not interval_payloads and len(histogram_series) >= 5:
        overlap = _support_overlap_fraction(histogram_series[2]["series"], histogram_series[4]["series"])
        if overlap is not None:
            interval_payloads.append(
                {
                    "from_lambda": _lambda_key(histogram_series[2]["current_lambda"]),
                    "to_lambda": _lambda_key(histogram_series[4]["current_lambda"]),
                    **overlap,
                }
            )

    if not interval_payloads:
        return {
            "series_count": len(histogram_series),
            "interval_count": 0,
            "intervals": [],
            "overlap_score_mean": None,
            "overlap_score_min": None,
        }

    scores = [item["overlap_fraction"] for item in interval_payloads if item["overlap_fraction"] is not None]
    worst_interval = min(interval_payloads, key=lambda item: item["overlap_fraction"] if item["overlap_fraction"] is not None else 1.0)
    return {
        "series_count": len(histogram_series),
        "interval_count": len(interval_payloads),
        "intervals": interval_payloads,
        "overlap_score_mean": _safe_mean(scores),
        "overlap_score_min": min(scores) if scores else None,
        "worst_interval": worst_interval,
    }


def _support_bounds(
    series: list[tuple[float, float]],
    *,
    negate_x: bool = False,
) -> tuple[float, float, float] | None:
    support = [(-x if negate_x else x) for x, count in series if count > 0]
    if not support:
        return None
    return min(support), max(support), sum(count for _, count in series if count > 0)


def _support_overlap_score(
    forward_bounds: tuple[float, float] | None,
    reverse_bounds: tuple[float, float] | None,
) -> float | None:
    if forward_bounds is None or reverse_bounds is None:
        return None
    forward_min, forward_max = forward_bounds
    reverse_min, reverse_max = reverse_bounds
    union_min = min(forward_min, reverse_min)
    union_max = max(forward_max, reverse_max)
    intersection_min = max(forward_min, reverse_min)
    intersection_max = min(forward_max, reverse_max)
    union_width = union_max - union_min
    if union_width == 0:
        return 1.0 if intersection_max >= intersection_min else 0.0
    return max(intersection_max - intersection_min, 0.0) / union_width


def _support_overlap_fraction(
    forward: list[tuple[float, float]],
    reverse: list[tuple[float, float]],
) -> dict[str, float | None] | None:
    forward_support = _support_bounds(forward)
    reverse_support = _support_bounds(reverse)
    reverse_reflected_support = _support_bounds(reverse, negate_x=True)
    if forward_support is None or reverse_support is None or reverse_reflected_support is None:
        return None

    forward_min, forward_max, forward_sample_count = forward_support
    reverse_raw_min, reverse_raw_max, reverse_sample_count = reverse_support
    reverse_reflected_min, reverse_reflected_max, _reverse_reflected_sample_count = reverse_reflected_support
    raw_score = _support_overlap_score((forward_min, forward_max), (reverse_raw_min, reverse_raw_max))
    reflected_score = _support_overlap_score(
        (forward_min, forward_max),
        (reverse_reflected_min, reverse_reflected_max),
    )
    reverse_transform = "identity"
    reverse_min = reverse_raw_min
    reverse_max = reverse_raw_max
    score = raw_score
    if reflected_score is not None and (score is None or reflected_score > score):
        reverse_transform = "negate"
        reverse_min = reverse_reflected_min
        reverse_max = reverse_reflected_max
        score = reflected_score
    return {
        "forward_support_min_kj_mol": forward_min,
        "forward_support_max_kj_mol": forward_max,
        "reverse_raw_support_min_kj_mol": reverse_raw_min,
        "reverse_raw_support_max_kj_mol": reverse_raw_max,
        "reverse_reflected_support_min_kj_mol": reverse_reflected_min,
        "reverse_reflected_support_max_kj_mol": reverse_reflected_max,
        "reverse_support_min_kj_mol": reverse_min,
        "reverse_support_max_kj_mol": reverse_max,
        "reverse_transform": reverse_transform,
        "overlap_fraction": score,
        "overlap_fraction_raw": raw_score,
        "overlap_fraction_reflected": reflected_score,
        "forward_sample_count": forward_sample_count,
        "reverse_sample_count": reverse_sample_count,
    }


def _build_job_metadata(spec: dict[str, Any], job_dir: Path) -> dict[str, Any]:
    mutation_group = spec.get("mutation_group", {})
    system = spec.get("system", {})
    protocol = spec.get("protocol", {})
    sites = mutation_group.get("sites", [])
    return {
        "job_id": spec.get("job_id", job_dir.name),
        "batch_id": spec.get("batch_id"),
        "job_dir": str(job_dir),
        "system_name": system.get("system_name"),
        "structure_source": system.get("structure_source"),
        "protocol_preset": protocol.get("preset"),
        "protocol_repeats": protocol.get("repeats"),
        "protocol_lambda_windows": protocol.get("lambda_windows"),
        "protocol_production_ps": protocol.get("production_ps"),
        "mutation_group_id": mutation_group.get("mutation_group_id"),
        "mutation_count": mutation_group.get("mutation_count"),
        "entity_side": mutation_group.get("entity_side"),
        "mutation_signature": "__".join(
            f"{site.get('chain_id')}:{site.get('wt')}{site.get('resseq')}{site.get('icode', '')}{site.get('mut')}@{site.get('entity_side')}"
            for site in sites
        ),
        "mutation_short_label": "--".join(
            f"{site.get('chain_id', '').lower()}-{site.get('wt', '').lower()}{site.get('resseq')}{site.get('icode', '').lower()}{site.get('mut', '').lower()}"
            for site in sites
        ),
        "sites": sites,
    }


def _parse_repeat_bar(repeat_dir: Path, *, leg: str, temperature_k: float) -> dict[str, Any]:
    dhdl_files = sorted(repeat_dir.glob("lambda_*/dhdl.xvg"))
    bar_dir = repeat_dir / "bar"
    interval_rows = _numeric_lines(bar_dir / "bar.xvg")
    integral_rows = _numeric_lines(bar_dir / "barint.xvg")
    overlap_payload = _histogram_overlap_payload(bar_dir / "histogram.xvg")
    worst_interval = overlap_payload.get("worst_interval") if overlap_payload is not None else None
    total_delta_kt = integral_rows[-1][1] if integral_rows else None
    stderr_kt = sqrt(sum(row[2] ** 2 for row in interval_rows if len(row) >= 3)) if interval_rows else None
    return {
        "leg": leg,
        "repeat_id": repeat_dir.name,
        "repeat_index": int(repeat_dir.name.replace("rep", "")) if repeat_dir.name.startswith("rep") else None,
        "repeat_dir": str(repeat_dir),
        "bar_dir": str(bar_dir),
        "observed_dhdl_count": len(dhdl_files),
        "observed_bar_intervals": len(interval_rows),
        "bar_complete": (bar_dir / "bar.xvg").exists() and (bar_dir / "barint.xvg").exists(),
        "histogram_present": (bar_dir / "histogram.xvg").exists(),
        "histogram_parsed": overlap_payload is not None,
        "histogram_series_count": overlap_payload["series_count"] if overlap_payload is not None else None,
        "overlap_interval_count": overlap_payload["interval_count"] if overlap_payload is not None else None,
        "overlap_score": overlap_payload["overlap_score_min"] if overlap_payload is not None else None,
        "overlap_score_mean": overlap_payload["overlap_score_mean"] if overlap_payload is not None else None,
        "overlap_forward_support_min_kj_mol": worst_interval["forward_support_min_kj_mol"] if worst_interval is not None else None,
        "overlap_forward_support_max_kj_mol": worst_interval["forward_support_max_kj_mol"] if worst_interval is not None else None,
        "overlap_reverse_support_min_kj_mol": worst_interval["reverse_support_min_kj_mol"] if worst_interval is not None else None,
        "overlap_reverse_support_max_kj_mol": worst_interval["reverse_support_max_kj_mol"] if worst_interval is not None else None,
        "overlap_reverse_raw_support_min_kj_mol": (
            worst_interval["reverse_raw_support_min_kj_mol"] if worst_interval is not None else None
        ),
        "overlap_reverse_raw_support_max_kj_mol": (
            worst_interval["reverse_raw_support_max_kj_mol"] if worst_interval is not None else None
        ),
        "overlap_reverse_reflected_support_min_kj_mol": (
            worst_interval["reverse_reflected_support_min_kj_mol"] if worst_interval is not None else None
        ),
        "overlap_reverse_reflected_support_max_kj_mol": (
            worst_interval["reverse_reflected_support_max_kj_mol"] if worst_interval is not None else None
        ),
        "overlap_reverse_transform": worst_interval["reverse_transform"] if worst_interval is not None else None,
        "overlap_score_raw": worst_interval["overlap_fraction_raw"] if worst_interval is not None else None,
        "overlap_score_reflected": worst_interval["overlap_fraction_reflected"] if worst_interval is not None else None,
        "overlap_forward_sample_count": worst_interval["forward_sample_count"] if worst_interval is not None else None,
        "overlap_reverse_sample_count": worst_interval["reverse_sample_count"] if worst_interval is not None else None,
        "overlap_intervals": overlap_payload["intervals"] if overlap_payload is not None else [],
        "delta_g_kt": total_delta_kt,
        "stderr_kt": stderr_kt,
        "delta_g_kcal_mol": _kt_to_kcal(total_delta_kt, temperature_k),
        "stderr_kcal_mol": _kt_to_kcal(stderr_kt, temperature_k),
    }


def collect_job_results(job_dir: Path) -> dict[str, Any]:
    spec = _load_job_spec(job_dir)
    protocol = spec.get("protocol", {})
    temperature_k = float(protocol.get("temperature_k", 310.0))
    expected_repeats = int(protocol.get("repeats", 0))
    expected_lambda_windows = int(protocol.get("lambda_windows", 0))
    max_repeat_delta = float(protocol.get("max_repeat_delta_kcal_mol", 0.0))

    legs: dict[str, Any] = {}
    repeat_results_by_leg: dict[str, list[dict[str, Any]]] = {}
    for leg in ("complex", "apo"):
        leg_root = job_dir / "legs" / leg
        repeats = [
            _parse_repeat_bar(repeat_dir, leg=leg, temperature_k=temperature_k)
            for repeat_dir in sorted(path for path in leg_root.glob("rep*") if path.is_dir())
        ]
        repeat_results_by_leg[leg] = repeats
        delta_values = [item["delta_g_kcal_mol"] for item in repeats if item["delta_g_kcal_mol"] is not None]
        overlap_scores = [item["overlap_score"] for item in repeats if item["overlap_score"] is not None]
        complete_repeats = [
            item
            for item in repeats
            if item["bar_complete"] and item["observed_dhdl_count"] == expected_lambda_windows and item["delta_g_kcal_mol"] is not None
        ]
        repeat_range = _safe_range(delta_values)
        legs[leg] = {
            "expected_repeats": expected_repeats,
            "expected_lambda_windows": expected_lambda_windows,
            "observed_repeats": len(repeats),
            "complete_repeats": len(complete_repeats),
            "delta_g_kcal_mol_mean": _safe_mean(delta_values),
            "delta_g_kcal_mol_stdev": _safe_stdev(delta_values),
            "delta_g_kcal_mol_range": repeat_range,
            "delta_g_kt_mean": _safe_mean([item["delta_g_kt"] for item in repeats if item["delta_g_kt"] is not None]),
            "bar_stderr_kcal_mol_mean": _safe_mean(
                [item["stderr_kcal_mol"] for item in repeats if item["stderr_kcal_mol"] is not None]
            ),
            "overlap_score_mean": _safe_mean(overlap_scores),
            "overlap_score_min": min(overlap_scores) if overlap_scores else None,
            "repeat_within_threshold": repeat_range <= max_repeat_delta if repeat_range is not None else None,
            "repeats": repeats,
        }

    complex_by_repeat = {item["repeat_id"]: item for item in repeat_results_by_leg["complex"] if item["delta_g_kcal_mol"] is not None}
    apo_by_repeat = {item["repeat_id"]: item for item in repeat_results_by_leg["apo"] if item["delta_g_kcal_mol"] is not None}
    common_repeat_ids = sorted(set(complex_by_repeat) & set(apo_by_repeat))
    ddg_repeats = []
    for repeat_id in common_repeat_ids:
        complex_item = complex_by_repeat[repeat_id]
        apo_item = apo_by_repeat[repeat_id]
        stderr_kcal = None
        if complex_item["stderr_kcal_mol"] is not None and apo_item["stderr_kcal_mol"] is not None:
            stderr_kcal = sqrt(complex_item["stderr_kcal_mol"] ** 2 + apo_item["stderr_kcal_mol"] ** 2)
        ddg_kcal = complex_item["delta_g_kcal_mol"] - apo_item["delta_g_kcal_mol"]
        ddg_repeats.append(
            {
                "repeat_id": repeat_id,
                "complex_delta_g_kcal_mol": complex_item["delta_g_kcal_mol"],
                "apo_delta_g_kcal_mol": apo_item["delta_g_kcal_mol"],
                "ddg_kcal_mol": ddg_kcal,
                "propagated_stderr_kcal_mol": stderr_kcal,
            }
        )

    ddg_values = [item["ddg_kcal_mol"] for item in ddg_repeats]
    ddg_range = _safe_range(ddg_values)
    ddg_summary = {
        **_build_job_metadata(spec, job_dir),
        "generated_at": utc_now(),
        "temperature_k": temperature_k,
        "kT_to_kcal_mol": KCAL_PER_MOL_PER_K * temperature_k,
        "ready": bool(ddg_repeats),
        "paired_repeat_count": len(ddg_repeats),
        "complex_delta_g_kcal_mol": legs["complex"]["delta_g_kcal_mol_mean"],
        "apo_delta_g_kcal_mol": legs["apo"]["delta_g_kcal_mol_mean"],
        "ddg_kcal_mol": _safe_mean(ddg_values),
        "ddg_repeat_stdev_kcal_mol": _safe_stdev(ddg_values),
        "ddg_repeat_range_kcal_mol": ddg_range,
        "ddg_bar_stderr_kcal_mol": _safe_mean(
            [item["propagated_stderr_kcal_mol"] for item in ddg_repeats if item["propagated_stderr_kcal_mol"] is not None]
        ),
        "repeat_within_threshold": ddg_range <= max_repeat_delta if ddg_range is not None else None,
        "repeats": ddg_repeats,
    }
    return {
        "metadata": _build_job_metadata(spec, job_dir),
        "bar_summary": {
            **_build_job_metadata(spec, job_dir),
            "generated_at": utc_now(),
            "temperature_k": temperature_k,
            "kT_to_kcal_mol": KCAL_PER_MOL_PER_K * temperature_k,
            "legs": legs,
            "ddg": ddg_summary,
        },
        "ddg_summary": ddg_summary,
    }


def build_qc_report(job_dir: Path, results: dict[str, Any] | None = None) -> dict[str, Any]:
    spec = _load_job_spec(job_dir)
    protocol = spec.get("protocol", {})
    expected_lambda_windows = int(protocol.get("lambda_windows", 0))
    max_repeat_delta = float(protocol.get("max_repeat_delta_kcal_mol", 0.0))
    max_bar_stderr = float(protocol.get("max_bar_stderr_kcal_mol", 10.0))
    overlap_threshold = float(protocol.get("overlap_threshold", 0.0))
    results = results or collect_job_results(job_dir)
    bar_summary = results["bar_summary"]

    warnings: list[str] = []
    failures: list[str] = []
    leg_reports: dict[str, Any] = {}
    overlap_assessment: dict[str, Any] = {"threshold": overlap_threshold, "legs": {}}
    repeat_spread_legs: list[tuple[str, float]] = []
    for leg, payload in bar_summary["legs"].items():
        repeat_reports = []
        for repeat in payload["repeats"]:
            repeat_warnings: list[str] = []
            repeat_failures: list[str] = []
            if repeat["observed_dhdl_count"] != expected_lambda_windows:
                repeat_failures.append(
                    f"{leg}:{repeat['repeat_id']} observed {repeat['observed_dhdl_count']} dhdl files, expected {expected_lambda_windows}"
                )
            if not repeat["bar_complete"]:
                repeat_failures.append(f"{leg}:{repeat['repeat_id']} missing BAR output files")
            if not repeat["histogram_present"]:
                repeat_warnings.append(f"{leg}:{repeat['repeat_id']} missing histogram.xvg")
            elif not repeat["histogram_parsed"]:
                repeat_warnings.append(f"{leg}:{repeat['repeat_id']} histogram.xvg could not be parsed for overlap assessment")
            elif repeat["overlap_score"] is not None and repeat["overlap_score"] < overlap_threshold:
                repeat_warnings.append(
                    f"{leg}:{repeat['repeat_id']} overlap score {repeat['overlap_score']:.3f} below threshold {overlap_threshold:.3f}"
                )
            repeat_reports.append(
                {
                    "repeat_id": repeat["repeat_id"],
                    "status": "fail" if repeat_failures else ("warning" if repeat_warnings else "pass"),
                    "warnings": repeat_warnings,
                    "failures": repeat_failures,
                    "overlap_score": repeat["overlap_score"],
                }
            )
            warnings.extend(repeat_warnings)
            failures.extend(repeat_failures)

        if payload["repeat_within_threshold"] is False:
            repeat_spread_legs.append((leg, float(payload["delta_g_kcal_mol_range"])))
            warnings.append(
                f"{leg} repeat spread {payload['delta_g_kcal_mol_range']:.3f} kcal/mol exceeds threshold {max_repeat_delta:.3f}"
            )
        if payload["bar_stderr_kcal_mol_mean"] is not None and payload["bar_stderr_kcal_mol_mean"] > max_bar_stderr:
            warnings.append(
                f"{leg} mean BAR stderr {payload['bar_stderr_kcal_mol_mean']:.3f} kcal/mol exceeds threshold {max_bar_stderr:.3f}"
            )
        leg_reports[leg] = {
            "status": "fail" if any(item["status"] == "fail" for item in repeat_reports) else (
                "warning"
                if any(item["status"] == "warning" for item in repeat_reports)
                or payload["repeat_within_threshold"] is False
                or (payload["bar_stderr_kcal_mol_mean"] is not None and payload["bar_stderr_kcal_mol_mean"] > max_bar_stderr)
                else "pass"
            ),
            "repeat_reports": repeat_reports,
            "expected_repeats": payload["expected_repeats"],
            "observed_repeats": payload["observed_repeats"],
            "complete_repeats": payload["complete_repeats"],
            "repeat_delta_kcal_mol_range": payload["delta_g_kcal_mol_range"],
            "repeat_within_threshold": payload["repeat_within_threshold"],
            "bar_stderr_kcal_mol_mean": payload["bar_stderr_kcal_mol_mean"],
            "overlap_score_mean": payload["overlap_score_mean"],
            "overlap_score_min": payload["overlap_score_min"],
        }
        overlap_assessment["legs"][leg] = {
            "status": leg_reports[leg]["status"],
            "overlap_score_mean": payload["overlap_score_mean"],
            "overlap_score_min": payload["overlap_score_min"],
            "repeat_scores": [
                {
                    "repeat_id": repeat["repeat_id"],
                    "overlap_score": repeat["overlap_score"],
                    "histogram_present": repeat["histogram_present"],
                    "histogram_parsed": repeat["histogram_parsed"],
                }
                for repeat in payload["repeats"]
            ],
        }

    if bar_summary["ddg"]["repeat_within_threshold"] is False:
        warnings.append(
            f"ddG repeat spread {bar_summary['ddg']['ddg_repeat_range_kcal_mol']:.3f} kcal/mol exceeds threshold {max_repeat_delta:.3f}"
        )
    if (
        bar_summary["ddg"]["ddg_bar_stderr_kcal_mol"] is not None
        and bar_summary["ddg"]["ddg_bar_stderr_kcal_mol"] > max_bar_stderr
    ):
        warnings.append(
            f"ddG BAR stderr {bar_summary['ddg']['ddg_bar_stderr_kcal_mol']:.3f} kcal/mol exceeds threshold {max_bar_stderr:.3f}"
        )

    observed_repeats = sum(payload["observed_repeats"] for payload in bar_summary["legs"].values())
    status = "pass"
    if observed_repeats == 0 and not bar_summary["ddg"]["ready"]:
        status = "not_evaluated"
    elif failures:
        status = "fail"
    elif warnings:
        status = "warning"

    sorted_repeat_spread_legs = [
        leg for leg, _range in sorted(repeat_spread_legs, key=lambda item: (-item[1], item[0]))
    ]
    primary_repeat_spread_leg = sorted_repeat_spread_legs[0] if sorted_repeat_spread_legs else None

    return {
        **_build_job_metadata(spec, job_dir),
        "generated_at": utc_now(),
        "status": status,
        "overlap_threshold": overlap_threshold,
        "overlap_assessment": overlap_assessment,
        "max_repeat_delta_kcal_mol": max_repeat_delta,
        "max_bar_stderr_kcal_mol": max_bar_stderr,
        "legs": leg_reports,
        "repeat_spread_legs": sorted_repeat_spread_legs,
        "primary_repeat_spread_leg": primary_repeat_spread_leg,
        "ddg_ready": bar_summary["ddg"]["ready"],
        "ddg_repeat_range_kcal_mol": bar_summary["ddg"]["ddg_repeat_range_kcal_mol"],
        "ddg_bar_stderr_kcal_mol": bar_summary["ddg"]["ddg_bar_stderr_kcal_mol"],
        "warnings": warnings,
        "failures": failures,
    }


def write_job_results(job_dir: Path) -> dict[str, Any]:
    results_dir = ensure_dir(job_dir / "results")
    results = collect_job_results(job_dir)
    qc_report = build_qc_report(job_dir, results=results)

    ddg_summary = results["ddg_summary"]
    bar_summary = results["bar_summary"]
    ddg_row = {
        "job_id": ddg_summary["job_id"],
        "mutation_group_id": ddg_summary["mutation_group_id"],
        "system_name": ddg_summary["system_name"],
        "protocol_preset": ddg_summary["protocol_preset"],
        "entity_side": ddg_summary["entity_side"],
        "mutation_count": ddg_summary["mutation_count"],
        "complex_delta_g_kcal_mol": ddg_summary["complex_delta_g_kcal_mol"],
        "apo_delta_g_kcal_mol": ddg_summary["apo_delta_g_kcal_mol"],
        "ddg_kcal_mol": ddg_summary["ddg_kcal_mol"],
        "ddg_repeat_stdev_kcal_mol": ddg_summary["ddg_repeat_stdev_kcal_mol"],
        "ddg_bar_stderr_kcal_mol": ddg_summary["ddg_bar_stderr_kcal_mol"],
        "paired_repeat_count": ddg_summary["paired_repeat_count"],
        "ready": ddg_summary["ready"],
    }

    bar_summary_path = results_dir / "bar_summary.json"
    ddg_summary_path = results_dir / "ddg_summary.json"
    ddg_tsv_path = results_dir / "ddg_summary.tsv"
    qc_report_path = results_dir / "qc_report.json"
    write_json(bar_summary_path, bar_summary)
    write_json(ddg_summary_path, ddg_summary)
    write_csv_rows(ddg_tsv_path, [ddg_row], list(ddg_row.keys()))
    write_json(qc_report_path, qc_report)
    return {
        "bar_summary": bar_summary,
        "ddg_summary": ddg_summary,
        "qc_report": qc_report,
        "paths": {
            "bar_summary": str(bar_summary_path),
            "ddg_summary": str(ddg_summary_path),
            "ddg_summary_tsv": str(ddg_tsv_path),
            "qc_report": str(qc_report_path),
        },
    }


def summarize_job(job_dir: Path) -> dict:
    stage_records = _read_stage_records(job_dir)
    spec = _load_job_spec(job_dir)
    results = collect_job_results(job_dir)
    qc_report = build_qc_report(job_dir, results=results)
    summary = {
        "generated_at": utc_now(),
        "job": _build_job_metadata(spec, job_dir),
        "job_dir": str(job_dir),
        "stage_count": len(stage_records),
        "stages": stage_records,
        "results": {
            "ddg": results["ddg_summary"],
            "bar": results["bar_summary"],
        },
        "progress": {
            "equilibrate": _equilibrate_progress(job_dir, spec),
            "sample": _sample_progress(job_dir, spec),
        },
        "qc": qc_report,
    }
    return summary


def write_job_summary(job_dir: Path) -> dict:
    result_payload = write_job_results(job_dir)
    summary = summarize_job(job_dir)
    summary["result_files"] = result_payload["paths"]
    write_json(job_dir / "report" / "summary.json", summary)
    write_yaml(job_dir / "report" / "summary.yml", summary)
    return summary


def _can_persist_job_summary(summary: dict[str, Any]) -> bool:
    if not summary["stages"]:
        return False
    latest_stage = summary["stages"][-1]
    return latest_stage.get("state") == "completed" and latest_stage.get("stage") in {"bar", "qc", "report"}


def _remove_transient_job_outputs(job_dir: Path) -> None:
    for path in [
        job_dir / "results" / "bar_summary.json",
        job_dir / "results" / "ddg_summary.json",
        job_dir / "results" / "ddg_summary.tsv",
        job_dir / "results" / "qc_report.json",
        job_dir / "report" / "summary.json",
        job_dir / "report" / "summary.yml",
    ]:
        if path.exists():
            path.unlink()


def _batch_qc_status(summary: dict[str, Any]) -> str:
    if not summary["stages"]:
        return "not_started"
    qc_stage = next((record for record in summary["stages"] if record.get("stage") == "qc"), None)
    if qc_stage is None or qc_stage.get("state") != "completed":
        return "not_evaluated"
    return summary["qc"]["status"]


def _batch_analyzable(summary: dict[str, Any]) -> bool:
    if not summary["stages"]:
        return False
    latest_stage = summary["stages"][-1]
    return latest_stage.get("stage") == "sample" and latest_stage.get("state") == "completed"


def _mutation_site_keys_from_spec(spec: dict[str, Any]) -> set[tuple[str, int, str]]:
    mutation_group = spec.get("mutation_group")
    if not isinstance(mutation_group, dict):
        return set()
    raw_sites = mutation_group.get("sites")
    if not isinstance(raw_sites, list):
        return set()
    site_keys: set[tuple[str, int, str]] = set()
    for site in raw_sites:
        if not isinstance(site, dict):
            continue
        try:
            resseq = int(site.get("resseq") or 0)
        except (TypeError, ValueError):
            continue
        site_keys.add(
            (
                str(site.get("chain_id") or "").strip(),
                resseq,
                str(site.get("icode") or "").strip().upper(),
            )
        )
    return site_keys


def _issue_residue_key(issue: dict[str, Any], *, partner: bool = False) -> tuple[str, int, str]:
    prefix = "partner_" if partner else ""
    try:
        resseq = int(issue.get(f"{prefix}resseq") or 0)
    except (TypeError, ValueError):
        resseq = 0
    return (
        str(issue.get(f"{prefix}chain_id") or "").strip(),
        resseq,
        str(issue.get(f"{prefix}icode") or "").strip().upper(),
    )


def _blocked_mutate_qc_is_auto_repairable(job_dir: Path, summary: dict[str, Any]) -> bool:
    if not summary["stages"]:
        return False
    latest_stage = summary["stages"][-1]
    if latest_stage.get("stage") != "mutate" or latest_stage.get("state") != "blocked_input":
        return False

    mutate_qc_path = job_dir / "artifacts" / "mutate_qc.json"
    if not mutate_qc_path.is_file():
        return False

    try:
        mutate_qc_payload = read_json(mutate_qc_path)
    except (OSError, ValueError, TypeError):
        return False

    site_keys = _mutation_site_keys_from_spec(_load_job_spec(job_dir))
    if not site_keys:
        return False

    legs = mutate_qc_payload.get("legs")
    if not isinstance(legs, dict):
        return False

    repairable_issue_found = False
    for leg_payload in legs.values():
        if not isinstance(leg_payload, dict):
            continue
        issues = leg_payload.get("inter_residue_heavy_atom_clashes")
        if not isinstance(issues, list) or not issues:
            continue
        auto_repair_summary = leg_payload.get("auto_repair_summary")
        if isinstance(auto_repair_summary, dict) and auto_repair_summary.get("attempted") and not auto_repair_summary.get("succeeded"):
            return False
        repairable_issues, blocking_issues = partition_inter_residue_sidechain_repairable_clashes(issues)
        if blocking_issues:
            return False
        for issue in repairable_issues:
            residue_keys = (
                _issue_residue_key(issue, partner=False),
                _issue_residue_key(issue, partner=True),
            )
            if any(key in site_keys for key in residue_keys):
                return False
        if repairable_issues:
            repairable_issue_found = True
    return repairable_issue_found


def _downstream_mutate_recovery_is_resumable(job_dir: Path, summary: dict[str, Any]) -> bool:
    if not summary["stages"]:
        return False
    latest_stage = summary["stages"][-1]
    stage = str(latest_stage.get("stage") or "").strip()
    state = str(latest_stage.get("state") or "").strip()
    if state != "blocked_input" or stage != "equilibrate":
        return False

    stage_states = {
        str(record.get("stage") or "").strip(): str(record.get("state") or "").strip()
        for record in summary["stages"]
        if str(record.get("stage") or "").strip()
    }
    if not any(stage_name in stage_states for stage_name in ("mutate", "build_legs")):
        return False

    diagnostic_code, _diagnostic_detail = _blocked_input_diagnostic(
        job_dir,
        stage,
        str(latest_stage.get("message") or "").strip(),
    )
    return diagnostic_code in {
        "equilibrate_missing_processed_gro",
        "equilibrate_invalid_processed_gro",
        "equilibrate_missing_hybrid_topology",
    }


def _batch_resumable(job_dir: Path, summary: dict[str, Any]) -> bool:
    if not summary["stages"]:
        return True
    latest_stage = summary["stages"][-1]
    latest_name = latest_stage.get("stage")
    latest_state = latest_stage.get("state")
    if latest_state == "completed":
        return latest_name in {"ingest", "prepare", "mutate", "build_legs", "equilibrate"}
    if latest_state == "blocked_input":
        return (latest_name == "mutate" and _blocked_mutate_qc_is_auto_repairable(job_dir, summary)) or (
            _downstream_mutate_recovery_is_resumable(job_dir, summary)
        )
    if latest_state in {"failed", "blocked_external", "stale_running"}:
        return latest_name in {"ingest", "prepare", "mutate", "build_legs", "equilibrate", "sample"}
    return False


def _batch_benchmark_qc_qualified(summary: dict[str, Any]) -> bool:
    qc_status = _batch_qc_status(summary)
    if qc_status in {"not_started", "not_evaluated", "fail"}:
        return False
    ddg_summary = summary["results"]["ddg"]
    if not ddg_summary["ready"]:
        return False
    ddg_bar_stderr = summary["qc"].get("ddg_bar_stderr_kcal_mol")
    max_bar_stderr = summary["qc"].get("max_bar_stderr_kcal_mol")
    if ddg_bar_stderr is None or max_bar_stderr is None:
        return False
    if ddg_bar_stderr > max_bar_stderr:
        return False

    ddg_repeat_range = summary["qc"].get("ddg_repeat_range_kcal_mol")
    max_repeat_delta = summary["qc"].get("max_repeat_delta_kcal_mol")
    if ddg_repeat_range is None or max_repeat_delta is None or ddg_repeat_range > max_repeat_delta:
        return False

    overlap_threshold = summary["qc"].get("overlap_threshold")
    overlap_legs = summary["qc"].get("overlap_assessment", {}).get("legs", {})
    leg_reports = summary["qc"].get("legs", {})
    for leg_name, overlap_payload in overlap_legs.items():
        if leg_reports.get(leg_name, {}).get("repeat_within_threshold") is False:
            return False
        repeat_scores = overlap_payload.get("repeat_scores", [])
        if not repeat_scores:
            return False
        if any(not item.get("histogram_present") or not item.get("histogram_parsed") for item in repeat_scores):
            return False
        overlap_score_min = overlap_payload.get("overlap_score_min")
        if overlap_score_min is None or overlap_threshold is None or overlap_score_min < overlap_threshold:
            return False

    return True


_PENDING_STAGE_CODES = {
    "ingest": "pending_prepare",
    "prepare": "pending_mutate",
    "mutate": "pending_build_legs",
    "build_legs": "pending_equilibrate",
    "equilibrate": "pending_sample",
    "sample": "pending_bar",
    "bar": "pending_qc",
    "qc": "pending_report",
}


def _read_optional_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        payload = read_json(path)
    except (OSError, ValueError, TypeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def _prepare_blocked_diagnostic(job_dir: Path, message: str) -> tuple[str, str]:
    prepare_qc = _read_optional_json(job_dir / "artifacts" / "prepare_qc.json")
    legs = prepare_qc.get("legs")
    if isinstance(legs, dict):
        if any(isinstance(payload, dict) and payload.get("blocking_incomplete_standard_residues") for payload in legs.values()):
            return "input_backbone_incomplete", message or "Prepared structure has backbone-incomplete residues."
        if any(
            isinstance(payload, dict) and payload.get("blocking_intra_residue_heavy_atom_clashes")
            for payload in legs.values()
        ):
            return "input_intra_residue_clash", message or "Prepared structure has same-residue heavy-atom clashes."
        if any(isinstance(payload, dict) and payload.get("inter_residue_heavy_atom_clashes") for payload in legs.values()):
            return "input_inter_residue_clash", message or "Prepared structure has inter-residue heavy-atom clashes."
    return "prepare_blocked_input", message or "Prepare stage is blocked by input structure issues."


def _mutate_blocked_diagnostic(job_dir: Path, message: str) -> tuple[str, str]:
    lowered = message.lower()
    if "insertion codes" in lowered:
        return "mutate_insertion_code_unsupported", message
    if "expected mutant coordinate files are missing" in lowered:
        return "mutate_missing_mutant_pdb", message
    if "prepared leg input is missing" in lowered:
        return "mutate_missing_leg_input", message

    mutate_qc = _read_optional_json(job_dir / "artifacts" / "mutate_qc.json")
    legs = mutate_qc.get("legs")
    if isinstance(legs, dict):
        if any(isinstance(payload, dict) and payload.get("inter_residue_heavy_atom_clashes") for payload in legs.values()):
            return "mutate_inter_residue_clash", message or "Mutated structure has inter-residue heavy-atom clashes."
        if any(isinstance(payload, dict) and payload.get("incomplete_standard_residues") for payload in legs.values()):
            return (
                "mutate_incomplete_standard_residue",
                message or "Mutated structure still contains incomplete standard residues.",
            )
        for payload in legs.values():
            if not isinstance(payload, dict):
                continue
            processed_gro_qc = payload.get("processed_gro_qc")
            if not isinstance(processed_gro_qc, dict) or not processed_gro_qc or processed_gro_qc.get("valid", True):
                continue
            reason = str(processed_gro_qc.get("reason") or "").strip().lower()
            if reason == "invalid_coordinate":
                return "mutate_invalid_coordinate", message or "Mutated processed.gro contains invalid coordinates."
            if reason == "isolated_residue_hydrogen":
                return (
                    "mutate_processed_gro_isolated_residue_hydrogen",
                    message or "Mutated processed.gro contains residue-local hydrogens far from their heavy-atom frame.",
                )
            if reason:
                return f"mutate_processed_gro_{reason}", message or f"Mutated processed.gro is invalid ({reason})."
            return "mutate_processed_gro_invalid", message or "Mutated processed.gro is invalid."
    return "mutate_blocked_input", message or "Mutate stage is blocked by input or topology issues."


def _format_processed_gro_invalid_detail(leg: str, gro_path: Path, summary: dict[str, Any]) -> str:
    reason = str(summary.get("reason") or "invalid_processed_gro").strip()
    detail = f"{leg} processed.gro is invalid ({reason})"
    line_number = summary.get("line_number")
    if line_number is not None:
        detail += f" at line {line_number}"
    residue_number = summary.get("residue_number")
    residue_name = str(summary.get("residue_name") or "").strip()
    atom_name = str(summary.get("atom_name") or "").strip()
    nearest_heavy_atom = str(summary.get("nearest_heavy_atom") or "").strip()
    nearest_heavy_distance_nm = summary.get("nearest_heavy_distance_nm")
    if residue_number is not None and residue_name and atom_name:
        detail += f": {residue_name}{residue_number} {atom_name}"
        if nearest_heavy_atom and isinstance(nearest_heavy_distance_nm, (int, float)):
            detail += f" vs {nearest_heavy_atom} {nearest_heavy_distance_nm:.3f} nm"
    return f"{detail}. Source: {gro_path}"


def _current_mutate_output_diagnostic(job_dir: Path) -> tuple[str, str] | None:
    mutate_qc = _read_optional_json(job_dir / "artifacts" / "mutate_qc.json")
    legs = mutate_qc.get("legs")
    if not isinstance(legs, dict):
        return None

    for leg, payload in legs.items():
        if not isinstance(payload, dict):
            continue
        raw_processed_path = str(payload.get("processed_gro") or "").strip()
        gro_path = Path(raw_processed_path) if raw_processed_path else job_dir / "legs" / str(leg) / "pmx" / "processed.gro"
        if not gro_path.is_file():
            continue
        current_summary = inspect_gro_file(gro_path)
        if current_summary.get("valid", False):
            continue
        reason = str(current_summary.get("reason") or "").strip().lower()
        detail = _format_processed_gro_invalid_detail(str(leg), gro_path, current_summary)
        if reason == "invalid_coordinate":
            return "mutate_invalid_coordinate", detail
        if reason:
            return f"mutate_processed_gro_{reason}", detail
        return "mutate_processed_gro_invalid", detail
    return None


def _current_mutate_output_status(job_dir: Path) -> dict[str, Any]:
    diagnostic = _current_mutate_output_diagnostic(job_dir)
    if diagnostic is None:
        return {
            "current_invalid_mutate_output": False,
            "current_invalid_mutate_output_code": "",
            "current_invalid_mutate_output_detail": "",
        }
    diagnostic_code, diagnostic_detail = diagnostic
    return {
        "current_invalid_mutate_output": True,
        "current_invalid_mutate_output_code": diagnostic_code,
        "current_invalid_mutate_output_detail": diagnostic_detail,
    }


def _qc_diagnostic(summary: dict[str, Any]) -> tuple[str, str]:
    qc_report = summary.get("qc", {})
    if not isinstance(qc_report, dict):
        return "qc_warning", "QC reported an unspecified issue."

    items = [
        *[str(item).strip() for item in qc_report.get("failures", []) if str(item).strip()],
        *[str(item).strip() for item in qc_report.get("warnings", []) if str(item).strip()],
    ]
    for item in items:
        lowered = item.lower()
        if "missing bar output files" in lowered:
            return "qc_missing_bar_output", item
        if "observed" in lowered and "dhdl files" in lowered:
            return "qc_missing_dhdl", item
        if "missing histogram.xvg" in lowered:
            return "qc_histogram_missing", item
        if "could not be parsed for overlap assessment" in lowered:
            return "qc_histogram_unparsed", item
        if "overlap score" in lowered and "below threshold" in lowered:
            return "qc_low_overlap", item
        if "repeat spread" in lowered:
            return "qc_repeat_spread", item
        if "bar stderr" in lowered:
            return "qc_bar_stderr", item

    status = str(qc_report.get("status") or "").strip().lower()
    if status == "fail":
        return "qc_fail", items[0] if items else "QC failed without a classified reason."
    return "qc_warning", items[0] if items else "QC reported a warning without a classified reason."


def _blocked_external_diagnostic(stage: str, message: str) -> tuple[str, str]:
    lowered = message.lower()
    if stage == "mutate":
        if "pmx is not available" in lowered:
            return "mutate_pmx_unavailable", message
        if "gromacs topology library" in lowered:
            return "mutate_forcefield_unavailable", message
        return "mutate_blocked_external", message or "Mutate stage is blocked by missing external dependencies."
    if stage == "equilibrate":
        if "gmxlib overlay" in lowered:
            return "equilibrate_gmxlib_unavailable", message
        return "equilibrate_blocked_external", message or "Equilibrate stage is blocked by missing external dependencies."
    if stage == "sample":
        if "gmxlib overlay" in lowered:
            return "sample_gmxlib_unavailable", message
        return "sample_blocked_external", message or "Sample stage is blocked by missing external dependencies."
    if stage == "bar":
        return "bar_blocked_external", message or "BAR analysis is blocked by missing external dependencies."
    return f"{stage or 'workflow'}_blocked_external", message or "Stage is blocked by missing external dependencies."


def _blocked_input_diagnostic(job_dir: Path, stage: str, message: str) -> tuple[str, str]:
    if stage == "prepare":
        return _prepare_blocked_diagnostic(job_dir, message)
    if stage == "mutate":
        return _mutate_blocked_diagnostic(job_dir, message)
    lowered = message.lower()
    if stage == "equilibrate":
        if "coordinate file is missing" in lowered:
            return "equilibrate_missing_processed_gro", message
        if "coordinate file is invalid" in lowered:
            return "equilibrate_invalid_processed_gro", message
        if "hybrid topology is missing" in lowered:
            return "equilibrate_missing_hybrid_topology", message
        return "equilibrate_blocked_input", message or "Equilibrate stage is blocked by missing inputs."
    if stage == "sample":
        if "equilibrated starting structure is missing" in lowered:
            return "sample_missing_npt_gro", message
        if "repeat-specific topology is missing" in lowered:
            return "sample_missing_repeat_topology", message
        return "sample_blocked_input", message or "Sample stage is blocked by missing inputs."
    return f"{stage or 'workflow'}_blocked_input", message or "Stage is blocked by missing required inputs."


def _failed_diagnostic(stage: str, message: str) -> tuple[str, str]:
    if stage == "equilibrate":
        return "equilibrate_failed", message or "Equilibrate stage failed."
    if stage == "sample":
        return "sample_failed", message or "Sample stage failed."
    if stage == "bar":
        return "bar_failed", message or "BAR analysis failed."
    return f"{stage or 'workflow'}_failed", message or "Stage failed."


def _job_diagnostic(job_dir: Path, summary: dict[str, Any]) -> dict[str, str]:
    stage_records = summary.get("stages", [])
    if not stage_records:
        return {
            "diagnostic_family": "not_started",
            "diagnostic_code": "not_started",
            "diagnostic_detail": "No stage records exist for this job.",
        }

    qc_status = _batch_qc_status(summary)
    if qc_status in {"warning", "fail"}:
        diagnostic_code, diagnostic_detail = _qc_diagnostic(summary)
        return {
            "diagnostic_family": "qc",
            "diagnostic_code": diagnostic_code,
            "diagnostic_detail": diagnostic_detail,
        }
    if qc_status == "pass":
        return {
            "diagnostic_family": "completed",
            "diagnostic_code": "qc_pass",
            "diagnostic_detail": "QC passed.",
        }

    latest_stage = stage_records[-1]
    stage = str(latest_stage.get("stage") or "").strip()
    state = str(latest_stage.get("state") or "").strip()
    message = str(latest_stage.get("message") or "").strip()

    if state in {"running", "stale_running"}:
        return {
            "diagnostic_family": "running",
            "diagnostic_code": f"{state}_{stage or 'unknown'}",
            "diagnostic_detail": message or f"{state} {stage or 'workflow'}",
        }

    current_mutate_output_diagnostic = _current_mutate_output_diagnostic(job_dir)
    if current_mutate_output_diagnostic is not None:
        diagnostic_code, diagnostic_detail = current_mutate_output_diagnostic
        return {
            "diagnostic_family": "mutate_setup",
            "diagnostic_code": diagnostic_code,
            "diagnostic_detail": diagnostic_detail,
        }
    if state == "blocked_input":
        diagnostic_code, diagnostic_detail = _blocked_input_diagnostic(job_dir, stage, message)
        diagnostic_family = "input_structure" if stage == "prepare" else (stage or "workflow")
        if stage == "mutate":
            diagnostic_family = "mutate_setup"
        return {
            "diagnostic_family": diagnostic_family,
            "diagnostic_code": diagnostic_code,
            "diagnostic_detail": diagnostic_detail,
        }
    if state == "blocked_external":
        diagnostic_code, diagnostic_detail = _blocked_external_diagnostic(stage, message)
        diagnostic_family = "mutate_setup" if stage == "mutate" else (stage or "external")
        return {
            "diagnostic_family": diagnostic_family,
            "diagnostic_code": diagnostic_code,
            "diagnostic_detail": diagnostic_detail,
        }
    if state == "failed":
        diagnostic_code, diagnostic_detail = _failed_diagnostic(stage, message)
        return {
            "diagnostic_family": stage or "workflow",
            "diagnostic_code": diagnostic_code,
            "diagnostic_detail": diagnostic_detail,
        }
    if state == "completed":
        return {
            "diagnostic_family": "pending" if stage != "report" else "completed",
            "diagnostic_code": "reported" if stage == "report" else _PENDING_STAGE_CODES.get(stage, f"completed_{stage or 'workflow'}"),
            "diagnostic_detail": message or f"{stage or 'workflow'} completed.",
        }

    return {
        "diagnostic_family": stage or "workflow",
        "diagnostic_code": f"{stage or 'workflow'}_{state or 'unknown'}",
        "diagnostic_detail": message or "Stage reached an unclassified state.",
    }


def write_batch_summary(batch_dir: Path) -> dict:
    jobs_dir = batch_dir / "jobs"
    job_summaries = []
    for job_dir in sorted(path for path in jobs_dir.iterdir() if path.is_dir()):
        existing_summary = summarize_job(job_dir)
        if _can_persist_job_summary(existing_summary):
            summary = write_job_summary(job_dir)
        else:
            _remove_transient_job_outputs(job_dir)
            summary = existing_summary
        diagnostic = _job_diagnostic(job_dir, summary)
        current_mutate_output_status = _current_mutate_output_status(job_dir)
        ddg_summary = summary["results"]["ddg"]
        job_summaries.append(
            {
                "job_id": ddg_summary["job_id"],
                "mutation_group_id": ddg_summary["mutation_group_id"],
                "protocol_preset": ddg_summary["protocol_preset"],
                "protocol_repeats": ddg_summary.get("protocol_repeats"),
                "protocol_lambda_windows": ddg_summary.get("protocol_lambda_windows"),
                "protocol_production_ps": ddg_summary.get("protocol_production_ps"),
                "latest_stage": summary["stages"][-1]["stage"] if summary["stages"] else "",
                "latest_stage_state": summary["stages"][-1]["state"] if summary["stages"] else "not_started",
                "stage_count": summary["stage_count"],
                "analyzable": _batch_analyzable(summary),
                "resumable": _batch_resumable(job_dir, summary),
                "ddg_kcal_mol": ddg_summary["ddg_kcal_mol"],
                "complex_delta_g_kcal_mol": ddg_summary.get("complex_delta_g_kcal_mol"),
                "apo_delta_g_kcal_mol": ddg_summary.get("apo_delta_g_kcal_mol"),
                "ddg_ready": ddg_summary["ready"],
                "ddg_bar_stderr_kcal_mol": summary["qc"].get("ddg_bar_stderr_kcal_mol"),
                "max_bar_stderr_kcal_mol": summary["qc"].get("max_bar_stderr_kcal_mol"),
                "qc_status": _batch_qc_status(summary),
                "benchmark_qc_qualified": _batch_benchmark_qc_qualified(summary),
                "complex_leg_qc_status": summary["qc"].get("legs", {}).get("complex", {}).get("status", ""),
                "apo_leg_qc_status": summary["qc"].get("legs", {}).get("apo", {}).get("status", ""),
                "complex_repeat_spread_kcal_mol": summary["qc"].get("legs", {}).get("complex", {}).get("repeat_delta_kcal_mol_range"),
                "apo_repeat_spread_kcal_mol": summary["qc"].get("legs", {}).get("apo", {}).get("repeat_delta_kcal_mol_range"),
                "repeat_spread_legs": ",".join(summary["qc"].get("repeat_spread_legs", [])),
                "primary_repeat_spread_leg": summary["qc"].get("primary_repeat_spread_leg", ""),
                "equilibrate_started_repeats": summary["progress"]["equilibrate"]["started_repeats"],
                "equilibrate_completed_repeats": summary["progress"]["equilibrate"]["completed_repeats"],
                "equilibrate_total_repeats": summary["progress"]["equilibrate"]["total_repeats"],
                "sample_started_windows": summary["progress"]["sample"]["started_windows"],
                "sample_completed_windows": summary["progress"]["sample"]["completed_windows"],
                "sample_total_windows": summary["progress"]["sample"]["total_windows"],
                "sample_active_leg": summary["progress"]["sample"].get("active_leg"),
                "sample_active_repeat_id": summary["progress"]["sample"].get("active_repeat_id"),
                "sample_active_lambda_id": summary["progress"]["sample"].get("active_lambda_id"),
                "sample_active_lambda_index": summary["progress"]["sample"].get("active_lambda_index"),
                "sample_active_phase": summary["progress"]["sample"].get("active_phase"),
                "sample_active_window": summary["progress"]["sample"].get("active_window"),
                "diagnostic_family": diagnostic["diagnostic_family"],
                "diagnostic_code": diagnostic["diagnostic_code"],
                "diagnostic_detail": diagnostic["diagnostic_detail"],
                "current_invalid_mutate_output": current_mutate_output_status["current_invalid_mutate_output"],
                "current_invalid_mutate_output_code": current_mutate_output_status["current_invalid_mutate_output_code"],
                "current_invalid_mutate_output_detail": current_mutate_output_status["current_invalid_mutate_output_detail"],
            }
        )
    payload = {
        "batch_dir": str(batch_dir),
        "generated_at": utc_now(),
        "jobs": job_summaries,
    }
    write_json(batch_dir / "reports" / "batch_summary.json", payload)
    write_yaml(batch_dir / "reports" / "batch_summary.yml", payload)
    return payload


def flag_censored_experimental_values(
    pairs: list[dict[str, Any]],
    *,
    complex_id_key: str = "complex_id",
    experimental_key: str = "experimental_ddg_kcal_mol",
) -> list[dict[str, Any]]:
    """Flag assay-saturation (censored) experimental ddG values.

    Rule: within one complex, an experimental value that (a) appears two or
    more times AND (b) equals the complex maximum is almost certainly an
    upper-detection-limit report (e.g. 1BJ1's repeated 3.69 ~ 500-fold and
    1CZ8's 4.10 ~ 1000-fold from the VEGF alanine-scan SPR assays; documented
    in docs/known_issues.md ISSUE-003). Such points cannot constrain
    correlation and should be excluded with an explicit, auditable reason.

    Returns the input rows annotated with ``censored_experimental`` (bool) and
    ``censored_reason`` (str).
    """
    by_complex: dict[str, list[dict[str, Any]]] = {}
    for row in pairs:
        by_complex.setdefault(str(row.get(complex_id_key, "")), []).append(row)

    annotated: list[dict[str, Any]] = []
    for complex_id, rows in by_complex.items():
        values = [float(r[experimental_key]) for r in rows if r.get(experimental_key) not in (None, "")]
        if not values:
            annotated.extend({**r, "censored_experimental": False, "censored_reason": ""} for r in rows)
            continue
        max_value = max(values)
        counts = {v: values.count(v) for v in set(values)}
        censored_values = {v for v, n in counts.items() if n >= 2 and v == max_value and v > 0.0}
        for row in rows:
            value = row.get(experimental_key)
            is_censored = value not in (None, "") and float(value) in censored_values
            annotated.append(
                {
                    **row,
                    "censored_experimental": is_censored,
                    "censored_reason": (
                        f"experimental_detection_limit (value {value} repeated "
                        f"{counts.get(float(value), 0)}x at complex maximum)"
                    )
                    if is_censored
                    else "",
                }
            )
    return annotated
