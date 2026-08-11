# Architecture

## Isolation

`abag-rbfep` is designed as a fresh project. It does not import modules from the
existing `platform` directories and does not write results into any pre-existing
project path.

## Workflow

The workflow stages are fixed and persisted per job:

1. `ingest`
2. `prepare`
3. `mutate`
4. `build_legs`
5. `equilibrate`
6. `sample`
7. `bar`
8. `qc`
9. `report`

Each stage writes:

- `stages/<stage>.json`: machine-readable status and message
- `artifacts/commands/<stage>.sh`: generated shell commands
- stage-specific artifacts under `artifacts/`, `results/`, and `report/`

Benchmark plan roots add another orchestration layer on top of per-batch
directories:

- `plan_index.json|yml`: all planned per-complex batches from `batch plan-abbind`
- `reports/run_summary.*`: selected job launches/resumes from `batch run-abbind`
- `reports/plan_summary.*`: aggregated per-job and per-batch benchmark status
- `reports/plan_jobs.csv` and `reports/plan_batches.csv`: flat exports for
  downstream analysis
- `reports/benchmark_metrics_qc_qualified.*` and
  `reports/benchmark_pairs_qc_qualified.csv`: QC-qualified benchmark views

Current implementation status:

- `prepare`: implemented with project-local PDB chain filtering
- `prepare`: backbone-incomplete residues still block the run, while
  sidechain-only missing heavy atoms are normalized and repaired with
  `PDBFixer` before the workflow proceeds to `pmx`/`pdb2gmx`; missing residue
  spans remain out of scope for v1
- `mutate`: implemented as real `pmx/gmx/pmx` command generation and optional execution
- `build_legs`: implemented for lambda directory and MDP generation
- `equilibrate`: implemented as real `gmx editconf/solvate/genion/grompp/mdrun`
- `equilibrate`: EM now uses `-DFLEXIBLE` and `constraints = none`, and the
  setup auto-retries with a two-stage larger cubic box ladder when the initial
  box hits a periodic-boundary shift fatal
- `sample`: implemented as real per-window `pre_relax -> pre_md -> grompp/mdrun -dhdl`
- `sample`/`bar`: reruns now skip completed window/repeat artifacts so resume
  can be incremental instead of whole-leg replay
- `bar`: real command execution when DHDL files exist
- `batch run-abbind`: refreshes canonical root benchmark reports after each run
- `batch run-abbind`: supports multi-worker scheduling with explicit GPU device
  assignment while keeping per-job stage manifests and resume semantics
- filtered `batch report-abbind`: writes selection-scoped reports under
  `reports/selections/` without overwriting the canonical root summary
- `batch rescue-abbind`: refreshes QC, selects completed `warning` jobs with
  actionable repeat-spread / BAR-stderr / overlap issues, and materializes a
  separate rescue plan root with more conservative protocol settings; the QC
  payload now also records structured per-leg repeat-spread diagnostics so
  rescue/watch tooling can distinguish `complex`-dominated from
  `apo`-dominated instability
- `AB-Bind` now has an explicit complex-level split manifest for
  `development/calibration/validation` at
  `benchmarks/ab_bind/splits/ab_bind_rbfe_core_v1_split_v1.yml`
- `qc`: implemented as structured output from observed BAR artifacts
- `report`: implemented as `BAR -> leg dG -> ddG` extraction plus job summaries
- benchmark orchestration: implemented for plan-root selection, job execution,
  and aggregated reporting across planned AB-Bind batches
- benchmark reporting: now separates all completed pairs from the
  `benchmark_qc_qualified` subset using propagated `ddG` BAR uncertainty
- benchmark reporting: now annotates each job with
  `diagnostic_family`/`diagnostic_code`/`diagnostic_detail` and aggregates
  `diagnostic_*_counts` into `plan_summary.*`, so validation failures can be
  grouped into structure-input, mutate-setup, sampling/QC, and in-progress
  buckets instead of being treated as one undifferentiated backlog
- protocol controls: implemented for `grompp_maxwarn_genion`,
  `grompp_maxwarn_equilibration`, and `grompp_maxwarn_sampling`
- legacy job protocol configs are hydrated against preset defaults at runtime,
  so new QC fields can be introduced without forcing existing plan roots to be
  replanned

## Result Model

Per job, the project now writes:

- `results/bar_summary.json`
  - parsed BAR outputs for `complex` and `apo`
  - per-repeat `delta_g_kt`, `delta_g_kcal_mol`, BAR stderr, and observed DHDL count
- `results/ddg_summary.json`
  - `complex_delta_g_kcal_mol`
  - `apo_delta_g_kcal_mol`
  - `ddg_kcal_mol`
  - repeat spread and propagated BAR stderr
- `results/ddg_summary.tsv`
  - flattened one-row export for batch collection
- `results/qc_report.json`
  - window-count checks
  - BAR artifact presence
  - repeat spread checks against protocol thresholds
  - BAR uncertainty checks
  - histogram-overlap assessment per leg and repeat using the better of raw and
    sign-reflected reverse support alignment

`report/summary.json` embeds stage history, QC, parsed BAR summaries, and the
formal `ddG` payload. Batch reports aggregate `ddg_kcal_mol`, readiness, and QC
state for each job. Plan-root benchmark reports aggregate the same fields across
multiple per-complex batches.

QC status semantics are now explicit:

- `not_started`: no stage files for the job
- `not_evaluated`: no BAR-backed free-energy result exists yet
- `pass` / `warning` / `fail`: QC evaluated against observed BAR artifacts

Current warning triggers include:

- missing optional BAR histogram output
- histogram present but unparseable for overlap assessment
- overlap score below `overlap_threshold`
- repeat spread above `max_repeat_delta_kcal_mol`
- mean leg BAR stderr above `max_bar_stderr_kcal_mol`
- propagated `ddG` BAR stderr above `max_bar_stderr_kcal_mol`

## Validated Real Cases

The current repository state has two real-case checkpoints:

- `examples/real_cases/1vfb`
  - clean crystal structure
  - validated end-to-end through `report` for both:
    - `V1` single-point `B:Y32F`
    - `V2` same-side double-point `B:Y32F + B:V34I`
  - current quick preset result bundles live under:
    - `runs/real_cases/1vfb_y32f_quick/`
    - `runs/real_cases/1vfb_y32f_v34i_quick/`
- `examples/real_cases/4dn4`
  - real crystal structure with sidechain-only atom gaps after extraction
  - validated as a positive real-case input after project-local prepare repair
  - the current quick rerun now completes through `report` for
    `M:V47I@antigen`; the repaired-input audit trail still lives in
    `artifacts/prepare_qc.json` plus `artifacts/mutate_qc.json`
  - the bundled quick preset remains an execution-validation preset only:
    the latest run finishes with `QC = warning` because overlap and BAR
    uncertainty stay outside production thresholds on this larger system
- `benchmarks/abbind`
  - `1VFB` benchmark-derived single mutations
    `H:Y32A@antibody`, `H:G31A@antibody`, and `C:Y23A@antigen`
  - validated through the new plan-root lifecycle:
    - `batch plan-abbind`
    - `batch run-abbind`
    - `batch report-abbind`
  - current root-level validation bundle lives under:
    - `runs/benchmarks/abbind_1vfb_runner_quick_plan/`
  - `1JRH` benchmark-derived same-side double-point mutation
    `I:M25L + I:I28V @ antigen`
  - validated through the same root-level lifecycle under:
    - `runs/benchmarks/abbind_1jrh_runner_truequick/`
  - `abbind_core_v1_quick_plan` also now has a cross-complex validation slice:
    - `1JRH` single mutation `I:T14V @ antigen` completed through `report`
    - `3NGB` single mutation `H:G54S @ antibody` completed through `report`
      after exercising the window-local `pre_relax/pre_md` sample stabilization
      path
    - `1YY9` antibody-side single mutations block in `prepare` because the
      extracted structures still contain incomplete standard residues
  - held-out validation now has two separate roots by design:
    - `abbind_core_v1_validation_plan` for the lighter validation tier
    - `abbind_core_v1_validation_priority_plan` for stronger single-point
      reruns on curated high-effect holdout cases
  - sidechain-only missing heavy atoms are now repaired before mutation with a
    local `PDBFixer` pass; backbone-incomplete residues still block early, and
    terminal oxygen records are stripped before `pdb2gmx` to avoid false
    internal-terminus handling across real chain gaps

## Scientific Boundaries

- `V1`: single-point, same thermodynamic cycle for antibody-side and
  antigen-side mutations
- `V2`: same-side double-point, same `complex - apo` cycle
- `V2.1`: charge-changing and cross-side double-point, separate design

See `docs/v2_1_charge_design_cn.md` for the current `V2.1` design constraints.

`gmx bar` is the formal BAR aggregator. pmx is used for mutation definition,
hybrid structure generation, and hybrid topology completion.
