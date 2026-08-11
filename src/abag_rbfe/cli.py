"""Command line entrypoint for abag-rbfep."""

from __future__ import annotations

import argparse
from pathlib import Path

from abag_pmx.mutations import build_mutation_group, load_mutation_groups_from_csv
from abag_rbfe.benchmark import (
    calibrate_ab_bind_plan,
    curate_ab_bind,
    materialize_ab_bind_inputs,
    plan_ab_bind_batches,
    plan_ab_bind_rescues,
    report_ab_bind_plan,
    run_ab_bind_plan,
)
from abag_rbfe.io_utils import write_json, write_yaml
from abag_rbfe.models import dataclass_to_dict
from abag_rbfe.paths import ProjectPaths, resolve_job_dir
from abag_rbfe.planning import build_batch_plan
from abag_rbfe.reporting import write_batch_summary, write_job_summary
from abag_rbfe.stages import resume_job, run_job


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="abag-rbfe")
    subparsers = parser.add_subparsers(dest="command", required=True)

    system_parser = subparsers.add_parser("system")
    system_subparsers = system_parser.add_subparsers(dest="system_command", required=True)
    prepare = system_subparsers.add_parser("prepare")
    prepare.add_argument("--name", required=True)
    prepare.add_argument("--input-structure", required=True)
    prepare.add_argument("--antibody-chains", required=True, help="Comma-separated chain IDs")
    prepare.add_argument("--antigen-chains", required=True, help="Comma-separated chain IDs")
    prepare.add_argument("--structure-source", default="experimental")
    prepare.add_argument("--note", action="append", default=[])
    prepare.add_argument("--output", required=True)

    mutation_parser = subparsers.add_parser("mutation")
    mutation_subparsers = mutation_parser.add_subparsers(dest="mutation_command", required=True)
    validate = mutation_subparsers.add_parser("validate")
    validate.add_argument("--mutations", required=True)
    validate.add_argument("--output", required=True)

    batch_parser = subparsers.add_parser("batch")
    batch_subparsers = batch_parser.add_subparsers(dest="batch_command", required=True)
    curate = batch_subparsers.add_parser("curate-abbind")
    curate.add_argument("--source-csv", required=True)
    curate.add_argument("--annotations")
    curate.add_argument("--output-dir", default=str(ProjectPaths.discover().benchmarks_root))

    materialize = batch_subparsers.add_parser("materialize-abbind")
    materialize.add_argument("--benchmark-root", default=str(ProjectPaths.discover().benchmarks_root))
    materialize.add_argument("--annotations")

    plan_abbind = batch_subparsers.add_parser("plan-abbind")
    plan_abbind.add_argument("--benchmark-root", default=str(ProjectPaths.discover().benchmarks_root))
    plan_abbind.add_argument("--protocol", required=True)
    plan_abbind.add_argument("--spec", choices=("core_v1", "core_v2"), default="core_v1")
    plan_abbind.add_argument("--runs-root")
    plan_abbind.add_argument("--batch-prefix", default="abbind")
    plan_abbind.add_argument("--complex-id", action="append", default=[])
    plan_abbind.add_argument("--split-name")
    plan_abbind.add_argument("--split-file")
    plan_abbind.add_argument("--limit", type=int)

    run_abbind = batch_subparsers.add_parser("run-abbind")
    run_abbind.add_argument("--plan-root", required=True)
    run_abbind.add_argument("--batch-id", action="append", default=[])
    run_abbind.add_argument("--complex-id", action="append", default=[])
    run_abbind.add_argument("--split-name")
    run_abbind.add_argument("--split-file")
    run_abbind.add_argument("--job-id", action="append", default=[])
    run_abbind.add_argument("--limit-batches", type=int)
    run_abbind.add_argument("--limit-jobs", type=int)
    run_abbind.add_argument("--from-stage")
    run_abbind.add_argument("--to-stage")
    run_abbind.add_argument("--max-workers", type=int)
    run_abbind.add_argument("--gpu-devices", help="Comma-separated GPU device IDs for execution scheduling.")
    run_abbind.add_argument("--execute", action="store_true")
    run_abbind.add_argument("--resume", action="store_true")

    report_abbind = batch_subparsers.add_parser("report-abbind")
    report_abbind.add_argument("--plan-root", required=True)
    report_abbind.add_argument("--extra-plan-root", action="append", default=[])
    report_abbind.add_argument("--batch-id", action="append", default=[])
    report_abbind.add_argument("--complex-id", action="append", default=[])
    report_abbind.add_argument("--split-name")
    report_abbind.add_argument("--split-file")
    report_abbind.add_argument("--limit-batches", type=int)

    calibrate_abbind = batch_subparsers.add_parser("calibrate-abbind")
    calibrate_abbind.add_argument("--plan-root", required=True)
    calibrate_abbind.add_argument("--extra-plan-root", action="append", default=[])
    calibrate_abbind.add_argument("--fit-split-name", default="calibration")
    calibrate_abbind.add_argument("--fit-extra-split-name", action="append", default=[])
    calibrate_abbind.add_argument("--predict-split-name", default="validation")
    calibrate_abbind.add_argument("--split-file")
    calibrate_abbind.add_argument(
        "--model",
        choices=(
            "linear",
            "side_linear",
            "quadratic",
            "stderr_quadratic",
            "logabs_stderr_quadratic",
            "expdecay_invstderr_quadratic",
            "hill_invstderr_quadratic",
            "hill_side_invstderr_quadratic",
        ),
        default="hill_side_invstderr_quadratic",
    )
    calibrate_abbind.add_argument("--fit-qc-qualified-only", action="store_true")

    rescue_abbind = batch_subparsers.add_parser("rescue-abbind")
    rescue_abbind.add_argument("--plan-root", required=True)
    rescue_abbind.add_argument("--extra-plan-root", action="append", default=[])
    rescue_abbind.add_argument("--runs-root")
    rescue_abbind.add_argument("--batch-prefix", default="abbind-rescue")
    rescue_abbind.add_argument("--batch-id", action="append", default=[])
    rescue_abbind.add_argument("--complex-id", action="append", default=[])
    rescue_abbind.add_argument("--split-name")
    rescue_abbind.add_argument("--split-file")
    rescue_abbind.add_argument("--limit-batches", type=int)
    rescue_abbind.add_argument("--job-id", action="append", default=[])
    rescue_abbind.add_argument("--repeat-increment", type=int, default=1)
    rescue_abbind.add_argument("--lambda-increment", type=int, default=4)
    rescue_abbind.add_argument("--production-scale", type=float, default=2.0)
    rescue_abbind.add_argument("--window-relax-em-scale", type=float, default=2.0)
    rescue_abbind.add_argument("--window-relax-md-scale", type=float, default=2.0)
    rescue_abbind.add_argument("--nvt-scale", type=float, default=2.0)
    rescue_abbind.add_argument("--npt-scale", type=float, default=2.0)
    rescue_abbind.add_argument("--force-repeat-increment", action="store_true")
    rescue_abbind.add_argument("--force-lambda-increment", action="store_true")
    rescue_abbind.add_argument("--prefer-active-alternate-source", action="store_true")
    rescue_abbind.add_argument("--require-active-alternate", action="store_true")
    rescue_abbind.add_argument("--allow-pass-qc-outlier-rescue", action="store_true")
    rescue_abbind.add_argument("--target-primary-repeat-spread-leg", action="store_true")
    rescue_abbind.add_argument("--require-target-primary-repeat-spread-leg", action="store_true")
    rescue_abbind.add_argument("--allow-targeted-leg-count-deepening", action="store_true")
    rescue_abbind.add_argument("--hotspot-complex-id", action="append", default=[])
    rescue_abbind.add_argument("--hotspot-job-id", action="append", default=[])
    rescue_abbind.add_argument("--hotspot-repeat-increment", type=int)
    rescue_abbind.add_argument("--hotspot-lambda-increment", type=int)
    rescue_abbind.add_argument("--hotspot-production-scale", type=float)
    rescue_abbind.add_argument("--hotspot-window-relax-em-scale", type=float)
    rescue_abbind.add_argument("--hotspot-window-relax-md-scale", type=float)
    rescue_abbind.add_argument("--hotspot-nvt-scale", type=float)
    rescue_abbind.add_argument("--hotspot-npt-scale", type=float)

    plan = batch_subparsers.add_parser("plan")
    plan.add_argument("--system", required=True)
    plan.add_argument("--mutations", required=True)
    plan.add_argument("--protocol", required=True)
    plan.add_argument("--batch-id")
    plan.add_argument("--runs-root")

    run = subparsers.add_parser("run")
    run.add_argument("identifier")
    run.add_argument("--batch-dir")
    run.add_argument("--from-stage")
    run.add_argument("--to-stage")
    run.add_argument("--execute", action="store_true")

    resume = subparsers.add_parser("resume")
    resume.add_argument("identifier")
    resume.add_argument("--batch-dir")
    resume.add_argument("--execute", action="store_true")

    analyze = subparsers.add_parser("analyze")
    analyze.add_argument("identifier")
    analyze.add_argument("--batch-dir")
    analyze.add_argument("--execute", action="store_true")

    report = subparsers.add_parser("report")
    report.add_argument("identifier")
    report.add_argument("--batch-dir")

    return parser


def _split_csv_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "system" and args.system_command == "prepare":
        payload = {
            "system_name": args.name,
            "input_structure": str(Path(args.input_structure).expanduser().resolve()),
            "structure_source": args.structure_source,
            "antibody_chains": _split_csv_list(args.antibody_chains),
            "antigen_chains": _split_csv_list(args.antigen_chains),
            "notes": args.note,
        }
        write_yaml(Path(args.output), payload)
        return 0

    if args.command == "mutation" and args.mutation_command == "validate":
        raw_groups = load_mutation_groups_from_csv(Path(args.mutations))
        payload = {
            "mutation_group_count": len(raw_groups),
            "mutation_groups": [
                dataclass_to_dict(
                    build_mutation_group(
                        mutation_group_id=item["mutation_group_id"],
                        sites=item["sites"],
                        allow_double_same_side=True,
                        allow_charge_change=True,
                    )
                )
                for item in raw_groups
            ],
        }
        write_json(Path(args.output), payload)
        return 0

    if args.command == "batch" and args.batch_command == "curate-abbind":
        curate_ab_bind(
            Path(args.source_csv),
            Path(args.output_dir),
            annotations_path=Path(args.annotations) if args.annotations else None,
        )
        return 0

    if args.command == "batch" and args.batch_command == "materialize-abbind":
        materialize_ab_bind_inputs(
            Path(args.benchmark_root),
            annotations_path=Path(args.annotations) if args.annotations else None,
        )
        return 0

    if args.command == "batch" and args.batch_command == "plan-abbind":
        payload = plan_ab_bind_batches(
            Path(args.benchmark_root),
            Path(args.protocol),
            spec_name=args.spec,
            runs_root=Path(args.runs_root) if args.runs_root else None,
            batch_prefix=args.batch_prefix,
            complex_ids=args.complex_id,
            split_name=args.split_name,
            split_path=Path(args.split_file) if args.split_file else None,
            limit=args.limit,
        )
        print(payload["plan_root"])
        return 0

    if args.command == "batch" and args.batch_command == "run-abbind":
        payload = run_ab_bind_plan(
            Path(args.plan_root),
            execute=args.execute,
            resume=args.resume,
            batch_ids=args.batch_id,
            complex_ids=args.complex_id,
            split_name=args.split_name,
            split_path=Path(args.split_file) if args.split_file else None,
            job_ids=args.job_id,
            limit_batches=args.limit_batches,
            limit_jobs=args.limit_jobs,
            from_stage=args.from_stage,
            to_stage=args.to_stage,
            max_workers=args.max_workers,
            gpu_devices=_split_csv_list(args.gpu_devices) if args.gpu_devices else None,
        )
        print(payload["plan_root"])
        return 0

    if args.command == "batch" and args.batch_command == "report-abbind":
        payload = report_ab_bind_plan(
            Path(args.plan_root),
            extra_plan_roots=[Path(item) for item in args.extra_plan_root],
            batch_ids=args.batch_id,
            complex_ids=args.complex_id,
            split_name=args.split_name,
            split_path=Path(args.split_file) if args.split_file else None,
            limit_batches=args.limit_batches,
        )
        print(payload["plan_root"])
        return 0

    if args.command == "batch" and args.batch_command == "calibrate-abbind":
        payload = calibrate_ab_bind_plan(
            Path(args.plan_root),
            extra_plan_roots=[Path(item) for item in args.extra_plan_root],
            fit_split_name=args.fit_split_name,
            fit_split_names=[args.fit_split_name, *args.fit_extra_split_name],
            predict_split_name=args.predict_split_name,
            split_path=Path(args.split_file) if args.split_file else None,
            model=args.model,
            fit_qc_qualified_only=args.fit_qc_qualified_only,
        )
        print(payload["reports_dir"])
        return 0

    if args.command == "batch" and args.batch_command == "rescue-abbind":
        payload = plan_ab_bind_rescues(
            Path(args.plan_root),
            extra_plan_roots=[Path(item) for item in args.extra_plan_root],
            runs_root=Path(args.runs_root) if args.runs_root else None,
            batch_prefix=args.batch_prefix,
            batch_ids=args.batch_id,
            complex_ids=args.complex_id,
            split_name=args.split_name,
            split_path=Path(args.split_file) if args.split_file else None,
            limit_batches=args.limit_batches,
            job_ids=args.job_id,
            repeat_increment=args.repeat_increment,
            lambda_increment=args.lambda_increment,
            production_scale=args.production_scale,
            window_relax_em_scale=args.window_relax_em_scale,
            window_relax_md_scale=args.window_relax_md_scale,
            nvt_scale=args.nvt_scale,
            npt_scale=args.npt_scale,
            force_repeat_increment=args.force_repeat_increment,
            force_lambda_increment=args.force_lambda_increment,
            prefer_active_alternate_source=args.prefer_active_alternate_source,
            require_active_alternate=args.require_active_alternate,
            allow_pass_qc_outlier_rescue=args.allow_pass_qc_outlier_rescue,
            target_primary_repeat_spread_leg=args.target_primary_repeat_spread_leg,
            require_target_primary_repeat_spread_leg=args.require_target_primary_repeat_spread_leg,
            allow_targeted_leg_count_deepening=args.allow_targeted_leg_count_deepening,
            hotspot_complex_ids=args.hotspot_complex_id,
            hotspot_job_ids=args.hotspot_job_id,
            hotspot_repeat_increment=args.hotspot_repeat_increment,
            hotspot_lambda_increment=args.hotspot_lambda_increment,
            hotspot_production_scale=args.hotspot_production_scale,
            hotspot_window_relax_em_scale=args.hotspot_window_relax_em_scale,
            hotspot_window_relax_md_scale=args.hotspot_window_relax_md_scale,
            hotspot_nvt_scale=args.hotspot_nvt_scale,
            hotspot_npt_scale=args.hotspot_npt_scale,
        )
        print(payload["plan_root"])
        return 0

    if args.command == "batch" and args.batch_command == "plan":
        batch_plan = build_batch_plan(
            system_path=Path(args.system),
            mutations_path=Path(args.mutations),
            protocol_path=Path(args.protocol),
            batch_id=args.batch_id,
            runs_root=Path(args.runs_root) if args.runs_root else None,
        )
        print(batch_plan.batch_dir)
        return 0

    if args.command == "run":
        job_dir = resolve_job_dir(args.identifier, args.batch_dir)
        run_job(job_dir, execute=args.execute, from_stage=args.from_stage, to_stage=args.to_stage)
        return 0

    if args.command == "resume":
        job_dir = resolve_job_dir(args.identifier, args.batch_dir)
        resume_job(job_dir, execute=args.execute)
        return 0

    if args.command == "analyze":
        job_dir = resolve_job_dir(args.identifier, args.batch_dir)
        run_job(job_dir, execute=args.execute, from_stage="bar", to_stage="report")
        return 0

    if args.command == "report":
        candidate = Path(args.identifier)
        if candidate.is_dir():
            if (candidate / "jobs").is_dir():
                write_batch_summary(candidate)
            else:
                write_job_summary(candidate)
            return 0
        if args.batch_dir:
            job_dir = resolve_job_dir(args.identifier, args.batch_dir)
            write_job_summary(job_dir)
            return 0
        batch_dir = ProjectPaths.discover().runs_root / args.identifier
        if batch_dir.is_dir():
            write_batch_summary(batch_dir)
            return 0
        parser.error("Identifier is neither a batch dir nor a job dir. Provide --batch-dir for job reports.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
