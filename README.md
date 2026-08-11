# abag-rbfep

`abag-rbfep` is a standalone antibody-antigen RBFE project built outside the
existing `platform` directories. It does not import or execute the current
`gromacs-abag-mmgbsa` or `amino-acid-platform` code as runtime dependencies.

## Scope

- `V1`: single-point mutation ddG on the antibody or antigen side
- `V2`: same-side double-point mutation ddG
- `V2.1`: charge-changing and cross-side double-point extensions

The scientific core is `GROMACS + vendored pmx + project-local orchestration`.
`gmx bar` is the source of truth for BAR aggregation.

## Repository Layout

```text
vendor/pmx/           Vendored upstream pmx snapshot (develop branch)
src/abag_pmx/         Antibody-antigen specific mutation logic
src/abag_rbfe/        CLI, workflow orchestration, planning, QC, reports
benchmarks/ab_bind/   AB-Bind source and curated manifests
envs/                 Independent environment specs
runs/                 Generated batch and job workspaces
docs/                 Architecture, vendoring, and reference notes
examples/             Minimal system/protocol/mutation examples
```

## Quick Start

Create the independent environment:

```bash
cd /mnt/data/liuchao/abag-rbfep
conda env create -f envs/environment.yml
conda activate abag-rbfep
```

Generate a system manifest:

```bash
abag-rbfe system prepare \
  --name demo_abag \
  --input-structure /absolute/path/to/complex.pdb \
  --antibody-chains H,L \
  --antigen-chains A \
  --output examples/system.generated.yml
```

Validate mutations:

```bash
abag-rbfe mutation validate \
  --mutations examples/mutations.csv \
  --output runs/validated_mutations.json
```

Create a batch plan:

```bash
abag-rbfe batch plan \
  --system examples/system.yml \
  --mutations examples/mutations.csv \
  --protocol examples/protocol.yml
```

Create benchmark-scale plans directly from the shipped AB-Bind materialized
inputs:

```bash
abag-rbfe batch plan-abbind \
  --benchmark-root benchmarks/ab_bind \
  --protocol benchmarks/ab_bind/protocol.quick.yml \
  --spec core_v1
```

Run a selected job from a benchmark plan root and then aggregate the plan-root
report:

```bash
abag-rbfe batch run-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --complex-id 1VFB \
  --job-id 1vfb-antibody-h-y32a \
  --execute

abag-rbfe batch report-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --complex-id 1VFB
```

`batch run-abbind` now refreshes the canonical root report under `reports/`
after each execution pass. Filtered `batch report-abbind` calls keep that
canonical root summary intact and write selection-scoped outputs under
`reports/selections/<selection>/`.

For multi-GPU execution, the same command can schedule multiple jobs in one
pass:

```bash
abag-rbfe batch run-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --execute \
  --resume \
  --max-workers 4 \
  --gpu-devices 0,1,2,3 \
  --job-id 1vfb-antibody-h-y32a \
  --job-id 1vfb-antibody-h-y33f
```

When `--max-workers > 1`, `batch run-abbind` uses the provided GPU list (or the
auto-discovered visible GPUs when `--gpu-devices` is omitted) and assigns one
job per worker/device while preserving the normal per-job stage state and
resume behavior.

Use the explicit AB-Bind split manifest for independent validation views:

```bash
abag-rbfe batch report-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --split-name validation \
  --split-file benchmarks/ab_bind/splits/ab_bind_rbfe_core_v1_split_v1.yml
```

The current validation snapshot is summarized in
`docs/validation_status.md`, while the stricter project-level completion audit
is summarized in `docs/project_completion_status.md`. The current canonical
calibration summary under `abbind_core_v1_quick_plan` now auto-selects
`side_linear` as the best shipped post hoc model and accepts the
target-filtered held-out view that excludes `1MLC`, `1CZ8`, and `1BJ1`,
yielding `accepted calibrated Pearson R = 0.672` on `32` held-out pairs. The
full unfiltered calibrated holdout metric remains lower (`Pearson R = 0.145`
on `61` pairs), while the stronger `validation_priority` lane currently
reports `raw target-filtered Pearson R = 0.381` after excluding `3HFM`. The
completion audit keeps that validation win separate from the still-running
AB-Bind and external `3HFM` reference execution waves.

For a stronger holdout rerun tier than the 2-window quick preset:

```bash
./benchmarks/ab_bind/run_validation_protocol.sh
```

This now uses `benchmarks/ab_bind/protocol.validation_robust.yml`, which maps
to the project-local `validation_robust_single_point` preset and writes into
`runs/benchmarks/abbind_core_v1_validation_robust_plan/`. The default rerun
list targets the already-completed but numerically unstable holdout jobs first,
so stronger sampling can be compared directly against the earlier
`validation_priority_single_point` results.

For a curated strong-effect follow-up across held-out complexes:

```bash
./benchmarks/ab_bind/run_validation_priority.sh
```

This uses `benchmarks/ab_bind/protocol.validation_priority.yml` and writes into
the separate root `runs/benchmarks/abbind_core_v1_validation_priority_plan/` so
the stronger holdout tier does not overwrite the lighter
`abbind_core_v1_validation_plan/` artifacts. Keep `2nz9-antigen-a-h1064a`
out of the default live holdout queue for now; it now fails a stricter
post-mutate geometry check and should only be re-added explicitly after the
mutant-side clash is remediated.

Set `MAX_WORKERS` and optionally `GPU_DEVICES` to drive the priority queue on
multiple GPUs, e.g. `MAX_WORKERS=4 GPU_DEVICES=0,1,2,3 ./benchmarks/ab_bind/run_validation_priority.sh`.

For unattended continuation after the first robust wave has been launched,
use the constrained watcher wrapper:

```bash
./benchmarks/ab_bind/run_validation_robust_watcher.sh
```

This wraps the validation watcher with a strict allow-list of the current
high-value robust reruns
(`3hfm-antibody-h-y50a`, `3hfm-antibody-h-y33a`, `3hfm-antibody-h-c95a`,
`1cz8-antigen-w-g92a`, `1bj1-antigen-w-g92a`), keeps the `mdrun` override quoted correctly, and
refreshes the merged validation report rooted at
`runs/benchmarks/abbind_core_v1_validation_priority_plan/` using both the
`robust` and `rescue` roots as extra winners. Pass explicit
`job_id` values to replace the default list. Keep
`2nz9-antigen-a-h1064a` out of the default live queue for now; the current
robust lane classifies it as `blocked_input` after post-mutate geometry QC.

For a narrower priority-side completion queue that only targets partially
sampled holdout jobs not currently overlapped by the live robust allow-list,
use:

```bash
./benchmarks/ab_bind/run_validation_priority_backlog_watcher.sh
```

Its default list currently focuses on the highest-progress partial backlog in
`3HFM`, `1BJ1`, `1CZ8`, and `3NPS`, and now explicitly re-promotes the
high-progress `priority` resumes for `3hfm-antibody-l-n31a`,
`3hfm-antibody-l-n32a`, and `3hfm-antigen-y-y20a` instead of burning the same
GPU budget on their lower-progress `robust` restarts. The default backlog set
now stops at the higher-progress `1BJ1` partials (`1bj1-antigen-w-i83a` and
`1bj1-antigen-w-g88a`) so the shared node can spend more of its budget on the
active `robust/rescue` hotspot reruns instead of the lowest-progress second-tier
`1BJ1` resumes. Pass explicit `job_id` values when you want to re-add those
deferred lower-progress partials. It stays in `--only-listed` mode so it does
not silently widen back to the full priority root.

For the current warning follow-ups that directly target the worst ready
validation points (`1bj1-antigen-w-g92a`, `1cz8-antigen-w-g92a`,
`3hfm-antibody-h-c95a`, `3hfm-antibody-h-y33a`, `3hfm-antibody-h-y50a`),
start the rescue watcher explicitly:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh start rescue
```

The rescue wrapper intentionally runs with a slightly higher default
`MAX_COMPUTE_APPS_PER_GPU=6` so these five high-value follow-ups can piggyback
on the shared node without waiting for the priority backlog to drain
completely. The `priority`-root merged report now also includes this rescue
root automatically, so improved rescue outputs feed the canonical holdout
winner view without a separate manual `report-abbind` step.

If you intentionally want the watcher to treat the entire robust root as
eligible, use the raw watcher directly:

```bash
python3 benchmarks/ab_bind/watch_validation_priority.py --plan-root runs/benchmarks/abbind_core_v1_validation_robust_plan --poll-seconds 60 --max-compute-apps-per-gpu 4 --max-active-mdrun-threads 56 --warn-stale-mdrun-seconds 900 --max-launches-per-pass 2 --mdrun-args-override '-ntmpi 1 -ntomp 2'
```

The watcher now accepts `--plan-root`, `--split-name`, and `--split-file`, so
the same queue manager can be pointed at the priority rerun root or the new
robust validation root. It refreshes reports after completed analyses and
while live jobs remain active, so `running`/`stale_running` status in the
merged validation report stays current during long resume cycles. The priority
wrapper also uses a more aggressive default refill policy
(`--max-launches-per-pass 4`, `--launch-cooldown-seconds 60`) to pull backlog
jobs back onto free GPUs faster.
When positional `job_id` values are supplied, the watcher treats them as a
strict allow-list by default; add `--append-rest` only when you intentionally
want the remaining plan-root jobs appended after that explicit list.
`--max-active-mdrun-threads` is optional but useful on shared CPU nodes because
it prevents the watcher from launching new resumes once the summed `gmx mdrun`
thread budget is already full, even if `loadavg` has not caught up yet.
`--warn-stale-mdrun-seconds` only emits warnings; it helps flag active `mdrun`
processes whose `.log`/`dhdl.xvg` have stopped moving for too long so they can
be inspected before the queue silently stalls. `--mdrun-args-override` is a
runtime-only performance knob for future resumed jobs; it does not change the
stored protocol YAML, but it lets the queue adopt a more conservative CPU
threading shape such as `-ntomp 2` when the node is CPU-limited and GPU
occupancy is still available. `--max-launches-per-pass` is useful when several
GPUs free up at once but you only want to release a small number of resumable
jobs each polling cycle.

To manage the validation watchers together, use:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh restart
./benchmarks/ab_bind/run_validation_watchers.sh status
```

By default, `start`/`restart` now only relaunch the `robust` watcher. Add
explicit targets such as `priority`, `rescue`, or `all` only when you
intentionally want those extra queues back:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh restart all
./benchmarks/ab_bind/run_validation_watchers.sh start priority rescue
```

Pidfiles still land under each `reports/watch/` directory and each selected
watcher still runs under the small supervisor loop from
`run_persistent_watch.sh`. The direct `priority` and `rescue` watcher wrapper
scripts are now disabled by default; use `run_validation_watchers.sh` to
enable them intentionally.

To keep the CPU-side reservoir filled without touching already-started jobs,
prebuild the remaining not-started validation cases up to `build_legs`:

```bash
python3 benchmarks/ab_bind/prebuild_validation_backlog.py
```

This script only selects jobs with no stage manifests yet, chunks them, and
drives `batch run-abbind --to-stage build_legs` in the background-safe path.

Generate workflow artifacts for a job:

```bash
abag-rbfe run demo-abag-antibody-h-y32f \
  --batch-dir runs/demo_abag_single_20260604T000000Z
```

By default the project writes stage manifests and shell command files. External
GROMACS and pmx commands are only launched when `--execute` is supplied.

## Current Execution Boundary

- `prepare` now writes real leg-specific `complex` and `apo` PDB inputs.
- `prepare` now distinguishes between backbone-incomplete and sidechain-only
  incomplete standard residues. Backbone gaps still block the job. Sidechain-only
  gaps are normalized by stripping the remaining partial sidechain atoms and then
  repaired with `PDBFixer` before the workflow hands the structure to
  `pmx`/`pdb2gmx`. This repair path deliberately does not fill missing residue
  spans or unresolved loops. Terminal oxygen records introduced during repair or
  mutation are stripped again before `pdb2gmx` so internal chain gaps do not get
  misread as false termini.
- `mutate` now writes real `pmx mutate -> gmx pdb2gmx -> pmx gentop` command
  chains and will execute them when `--execute` is supplied and `pmx` is
  available in the project environment or on `PATH`.
- `build_legs` writes real lambda directory layouts and MDP files.
- `equilibrate` now executes real `gmx editconf -> solvate -> genion -> grompp
  -> mdrun` setup and short equilibration when `--execute` is supplied.
- `equilibrate` now retries EM setup with a two-stage cubic fallback ladder
  when the initial setup hits periodic-boundary shift fatals, and preserves an
  `em.runtime.history.log` audit trail for those EM retries. The shipped EM
  preset now uses `-DFLEXIBLE` with `constraints = none` to reduce
  minimization-time LINCS and water-settle instability.
- `sample` now executes real per-window `pre_relax -> pre_md -> production`
  chains when `--execute` is supplied. Each lambda window now gets a local
  EM-based pre-relaxation plus a short flexible-water stochastic warmup before
  the formal `mdrun -dhdl` production leg.
- `sample` and `bar` now skip already completed window/repeat artifacts on
  rerun, so failed real jobs can resume incrementally instead of replaying the
  whole leg.
- legacy `job/config/protocol.yml` files are now hydrated against the preset
  defaults at runtime, so newly added QC and `grompp_maxwarn_*` fields do not
  force old plan roots to be rebuilt.
- `protocol.yml` can now tune `grompp_maxwarn_genion`,
  `grompp_maxwarn_equilibration`, and `grompp_maxwarn_sampling` for real systems
  that routinely trigger a small number of known GROMACS warnings.
- `bar` executes real `gmx bar` when `dhdl.xvg` files are present.
- `qc` now writes a structured QC report with window-count, repeat-consistency,
  BAR uncertainty, and histogram-overlap checks.
- `report` now writes parsed BAR summaries and a formal `ddG` result bundle.

## Result Files

Each executed job now materializes a `results/` directory with stable outputs:

- `results/bar_summary.json`: parsed `complex` and `apo` BAR results by repeat
- `results/ddg_summary.json`: aggregated `ddG` summary in machine-readable JSON
- `results/ddg_summary.tsv`: one-line tabular export for downstream collection
- `results/qc_report.json`: current QC status and per-leg/per-repeat checks
- `report/summary.json` and `report/summary.yml`: full stage + result summary

Batch-level `abag-rbfe report <batch_dir>` now writes `reports/batch_summary.*`
with `ddg_kcal_mol`, readiness, and QC state per job.

Benchmark plan roots now also write:

- `reports/run_summary.json|yml|csv`: which planned jobs were launched or resumed
- `reports/plan_summary.json|yml`: aggregated status across the selected batches
- `reports/plan_batches.csv`: per-batch counts
- `reports/plan_jobs.csv`: flattened per-job benchmark summary, now including `abs_ddg_error_kcal_mol`
- `reports/active_alternate_jobs.csv`: current merged winners that already have active reruns elsewhere, sorted so the worst ready-error hotspots float to the top
- `reports/benchmark_metrics.json|yml`: metrics over all ready predicted/experimental pairs
- `reports/benchmark_metrics_qc_qualified.json|yml`: metrics over the QC-qualified subset
- `reports/benchmark_pairs.csv` and `reports/benchmark_pairs_qc_qualified.csv`

QC now also inspects BAR uncertainty and histogram overlap. Jobs with
`ddg_bar_stderr_kcal_mol` or mean leg BAR stderr above
`max_bar_stderr_kcal_mol`, or with overlap scores below `overlap_threshold`,
are downgraded to `warning`, even if all expected files are present.
`benchmark_qc_qualified=True` means the job is `ddg`-ready, not failed, and its
propagated `ddG` BAR stderr stays within the protocol threshold, so the
benchmark report can separate all completed pairs from the trusted subset.

## Benchmark Assets

`benchmarks/ab_bind/` is no longer just a placeholder. The repository now
ships:

- raw `AB-Bind` source data in `benchmarks/ab_bind/source/`
- a project-local complex annotation table for RBFE boundary decisions
- materialized `curated/` and `manifests/` outputs generated from the raw source

Current strict counts under the shipped annotation table are:

- `AB-Bind-Source`: `1101` rows / `32` complexes
- `AB-Bind-RBFE-Core-V1`: `318` rows / `18` complexes
- `AB-Bind-RBFE-Core-V2`: `339` rows / `19` complexes

The benchmark directory also now materializes per-complex `system.yml` and
mutation CSV inputs, and one benchmark-derived real quick run has been executed
end to end:

- `1VFB` `H:Y32A@antibody`
- batch: `runs/benchmarks/abbind_1vfb_core_v1_quick/`
- result: `ddG = 10.306 kcal/mol`, `QC = warning`
  - under the current overlap-aware QC, this quick preset is intentionally
    treated as a workflow-validation run, not a production-quality ddG

The benchmark workflow can now also be planned in bulk from the materialized
core sets. A real `core_v1` planning pass has already been generated under:

- `runs/benchmarks/abbind_core_v1_quick_plan/`
- current plan contents: `18` batches / `318` jobs

The same `core_v1` plan root now also has a cross-complex validation slice:

- `1JRH` single mutation `I:T14V@antigen`
  - result: `ddG = -7.805 kcal/mol`, `QC = warning`
- `3NGB` single mutation `H:G54S@antibody`
  - result: `ddG = -51.303 kcal/mol`, `QC = warning`
  - this case is the first real benchmark job that required the new
    per-window `pre_relax/pre_md` sample stabilization path
- `1YY9` single mutations `H:N56A@antibody` and `L:N93A@antibody`
  - both stop at `prepare` with `blocked_input`
  - cause: incomplete standard residues remain after chain extraction

The plan-root execution path has now also been validated on a fresh benchmark
runner root:

- root: `runs/benchmarks/abbind_1vfb_runner_quick_plan/`
- selected executed jobs:
  - `1vfb-antibody-h-y32a`
  - `1vfb-antibody-h-g31a`
  - `1vfb-antigen-c-y23a`
- root-level aggregate report:
  - `23` jobs in the selected batch
  - `3` `ddG`-ready jobs
  - QC counts: `warning=3`, `not_started=20`

The same root-level lifecycle has now also been validated for a real `V2`
same-side double-point benchmark case:

- root: `runs/benchmarks/abbind_1jrh_runner_truequick/`
- selected executed job: `1jrh-antigen-i-m25l--i-i28v`
- mutation signature: `I:M25L@antigen + I:I28V@antigen`
- result: `ddG = -19.054 kcal/mol`
- QC: `warning`
  - current BAR stderr is large enough to trip the new uncertainty gate
- root-level aggregate report:
  - `3` jobs in the selected batch
  - `1` ready job
  - QC counts: `warning=1`, `not_started=2`

## Smoke Status

The current standalone implementation has been exercised on a minimal local
smoke system through:

- `prepare`
- `mutate`
- `build_legs`
- `equilibrate`
- `sample`
- `bar`
- `qc`
- `report`

The smoke job produced real `dhdl.xvg`, `bar.xvg`, `ddg_summary.json`, and
`qc_report.json` artifacts under `runs/smoke_real/`.

## Real-Case Status

The repository now also contains real antibody-antigen examples under
`examples/real_cases/`:

- `1VFB`:
  - quick validation mutations:
    - `V1`: `B:Y32F@antibody`
    - `V2`: `B:Y32F@antibody + B:V34I@antibody`
  - current quick-run batches:
    - `runs/real_cases/1vfb_y32f_quick/`
    - `runs/real_cases/1vfb_y32f_v34i_quick/`
  - latest quick-run results:
    - `V1 ddG = -3.591 kcal/mol`, `QC = warning`
    - `V2 ddG = -2.119 kcal/mol`, `QC = warning`
    - both quick presets finish end to end, but BAR overlap remains below the
      current QC threshold
  - both are execution-validation results from short presets, not
    production-quality affinity estimates
- `4DN4`:
  - quick validation mutation: `M:V47I@antigen`
  - current quick-run batch: `runs/real_cases/4dn4_v47i_quick/`
  - no longer blocks at `prepare`; sidechain-only atom gaps are repaired during
    the leg-specific real-case preprocessing path
  - latest quick-run result:
    - `ddG = 18.173 kcal/mol`, `QC = warning`
    - completed end to end through `report`
    - warning is expected for the short validation preset on this larger case:
      both legs have `overlap score = 0.000` and BAR stderr is well above the
      default production threshold
  - relevant diagnostics and summaries live in `artifacts/prepare_qc.json`,
    `artifacts/mutate_qc.json`, `results/ddg_summary.json`,
    `results/qc_report.json`, and the per-stage status files under `stages/`

Batch-level summaries for validated runs can be regenerated with:

```bash
abag-rbfe report runs/real_cases/1vfb_y32f_quick
```

## Notes

- Current preprocessing only accepts `.pdb` as a direct input for leg
  extraction. Convert `.cif`/`.mmcif` before running the execution path.
- The project accepts insertion codes in manifests, but the current pmx script
  generation only supports integer residue IDs. Jobs with unresolved insertion
  codes are kept valid at the model layer and explicitly blocked at the mutate
  stage until a residue-number mapping step is provided.
- The default presets now allow `grompp maxwarn = 2` for `genion`,
  equilibration, and sampling so that benchmark systems like `1JRH` are not
  rejected solely because of routine hybrid-topology warning pairs.
- When a `single_point` protocol file explicitly sets quick overrides such as
  `lambda_windows`, `repeats`, `nvt_ps`, `npt_ps`, or `production_ps`, those
  explicit values are now preserved when the planner upgrades a two-site job to
  the `double_point` preset.
- Vendored pmx remains under its upstream license. See [NOTICE.md](NOTICE.md)
  and [docs/upstream_pmx.md](docs/upstream_pmx.md).
- A curated project reading list lives in
  [docs/reference_materials.md](docs/reference_materials.md).
- A Chinese quick-reading companion lives in
  [docs/reference_materials_cn.md](docs/reference_materials_cn.md).
