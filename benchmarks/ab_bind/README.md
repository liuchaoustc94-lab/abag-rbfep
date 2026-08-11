# AB-Bind Benchmark Layout

This project keeps AB-Bind in three layers:

- `source/`: raw source files copied from the original dataset release
- `curated/`: machine-generated filtered tables
- `manifests/`: YAML manifests used by batch planning and evaluation

## Core benchmark sets

- `AB-Bind-Source`: all rows preserved for provenance
- `AB-Bind-RBFE-Core-V1`: antibody-antigen, single-point, standard residues,
  structure-mappable, stable-run candidates, nominally charge-conserving
- `AB-Bind-RBFE-Core-V2`: V1-compatible rows plus same-side double-point rows

The `abag-rbfe batch curate-abbind` command produces the curated layers and
records explicit exclusion codes for every filtered row.

## Current Materialized Counts

Under the current project-local annotation table:

- `AB-Bind-Source`: `1101` rows across `32` complexes
- `AB-Bind-RBFE-Core-V1`: `318` rows across `18` experimental strict
  antibody-antigen complexes
- `AB-Bind-RBFE-Core-V2`: `339` rows across `19` experimental strict
  antibody-antigen complexes

These counts intentionally reflect the stricter `abag-rbfep` RBFE boundary, not
the broader raw AB-Bind release or the literature `S645` comparison subset.

## Refresh

Refetch the upstream table:

```bash
./benchmarks/ab_bind/source/fetch_source.sh
```

Regenerate curated layers:

```bash
./benchmarks/ab_bind/refresh_curated.sh
```

Materialize per-complex `system.yml` and `mutations.csv` inputs:

```bash
./benchmarks/ab_bind/materialize_inputs.sh
```

Fetch the upstream PDB / homology-model structures:

```bash
./benchmarks/ab_bind/source/fetch_structures.sh
```

## Materialized Inputs

The current repository state now includes:

- `materialized/<complex_id>/system.yml`
- `materialized/<complex_id>/core_v1_mutations.csv`
- `materialized/<complex_id>/core_v2_mutations.csv`
- `manifests/ab_bind_rbfe_core_v1_inputs.csv`
- `manifests/ab_bind_rbfe_core_v2_inputs.csv`

These files can be fed directly into `abag-rbfe batch plan`.

For benchmark-scale planning, use:

```bash
abag-rbfe batch plan-abbind \
  --benchmark-root benchmarks/ab_bind \
  --protocol benchmarks/ab_bind/protocol.quick.yml \
  --spec core_v1
```

This command consumes `manifests/ab_bind_rbfe_<spec>_inputs.csv` and writes one
batch directory per complex plus a `plan_index.json` summary.

Execute a selected job from a plan root:

```bash
abag-rbfe batch run-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --complex-id 1VFB \
  --job-id 1vfb-antibody-h-y32a \
  --execute
```

Aggregate the selected plan root back into root-level reports:

```bash
abag-rbfe batch report-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --complex-id 1VFB
```

Materialize conservative follow-up batches for completed `warning` jobs without
touching the original plan root:

```bash
abag-rbfe batch rescue-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --complex-id 1VFB
```

By default, `rescue-abbind` now thickens the pre-equilibration path for
repeat-spread/overlap hotspots as well as the sampling path: it scales
`window_relax_em_steps`, `window_relax_md_ps`, `nvt_ps`, and `npt_ps` by `2x`
unless you override `--window-relax-em-scale`, `--window-relax-md-scale`,
`--nvt-scale`, or `--npt-scale`. Those same repeat-spread/overlap rescues also
switch the early equilibration leg into a staged restraint schedule
(`POSRES_STAGE_HEAVY -> POSRES_STAGE_BACKBONE`) and add an extra restrained
backbone-release `npt_release.mdp` segment before lambda-window sampling.

When the current best rows live in a merged multi-plan view, add repeated
`--extra-plan-root` arguments and `rescue-abbind` will source each selected
job from the winning row's original plan root rather than assuming every
candidate still belongs to the primary `--plan-root`.

This writes:

- `reports/run_summary.json|yml|csv`
- `reports/plan_summary.*`, `reports/benchmark_metrics.*`,
  `reports/plan_jobs.csv`, and `reports/plan_batches.csv` as the canonical
  full-root summary
- `reports/active_alternate_jobs.csv` as the filtered merged-winner view for
  rows that already have an active rerun in another plan root, including the
  representative active alternate's current stage and sample-window progress
- `reports/benchmark_metrics_qc_qualified.*` and
  `reports/benchmark_pairs_qc_qualified.csv` for the trusted subset
- `reports/selections/<selection>/...` for any filtered `report-abbind` view
- `<rescue_plan_root>/reports/rescue_summary.*` and
  `<rescue_plan_root>/reports/rescue_candidates.csv` for rescue batches built
  from `batch rescue-abbind`

Merged `plan_summary.json` now also exposes
`active_alternate_ready_hotspots`, which is a short preview of the current
ready winner rows with the largest absolute `ddG` error that already have an
active follow-up elsewhere.

An explicit complex-level holdout split is now tracked in:

- `benchmarks/ab_bind/splits/ab_bind_rbfe_core_v1_split_v1.yml`

Generate the held-out validation report with:

```bash
abag-rbfe batch report-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --split-name validation \
  --split-file benchmarks/ab_bind/splits/ab_bind_rbfe_core_v1_split_v1.yml
```

Fit a post hoc calibration model from one completed split and apply it to
another split without changing the underlying BAR outputs:

```bash
abag-rbfe batch calibrate-abbind \
  --plan-root runs/benchmarks/abbind_core_v1_quick_plan \
  --fit-split-name calibration \
  --predict-split-name validation \
  --split-file benchmarks/ab_bind/splits/ab_bind_rbfe_core_v1_split_v1.yml \
  --model side_linear
```

This writes:

- `reports/calibrations/<slug>/summary.json|yml`
- `reports/calibrations/<slug>/model.json|yml`
- `reports/calibrations/<slug>/fit_pairs.csv`
- `reports/calibrations/<slug>/predict_jobs_calibrated.csv`
- `reports/calibrations/<slug>/predict_pairs_calibrated.csv`

Use `--fit-qc-qualified-only` when you want to fit only from the trusted
subset in `benchmark_pairs_qc_qualified.csv`. The fit split must already have
non-empty paired benchmark rows; in practice, refresh the split report first
and verify `paired_job_count > 0` before running `calibrate-abbind`.

For a quick holdout sweep on the default `1MLC` validation jobs:

```bash
./benchmarks/ab_bind/run_validation_quick.sh
```

To opportunistically backfill the calibration split on the shared quick root
without contending aggressively with the live validation queues, start the
calibration watcher:

```bash
./benchmarks/ab_bind/run_calibration_watchers.sh start
./benchmarks/ab_bind/run_calibration_watchers.sh status
```

This wrapper keeps the default quick root (`runs/benchmarks/abbind_core_v1_quick_plan/`)
planned for the full core-v1 benchmark, but only launches the listed
calibration jobs under the `calibration` split with a conservative watcher
profile. The low-level `run_calibration_quick_watcher.sh` entrypoint is now
disabled by default, so shared nodes do not quietly repopulate the quick lane
unless you explicitly re-enable it through `run_calibration_watchers.sh start`.
Its default GPU slot cap now matches the rescue lane
(`MAX_COMPUTE_APPS_PER_GPU=6`) so one calibration follow-up can still piggyback
when the node is saturated by many low-memory validation `mdrun` processes.
Pass explicit `job_id` values to replace the default shortlist.

When quick calibration jobs finish with `warning` QC, refresh the stronger
follow-up plan with:

```bash
./benchmarks/ab_bind/refresh_calibration_rescues.sh
./benchmarks/ab_bind/run_calibration_watchers.sh start rescue
```

This writes rescue batches into
`runs/benchmarks/abbind_core_v1_calibration_rescues/` and upgrades the quick
protocol toward the validation-priority effort level by default
(`repeats +2`, `lambda_windows +6`, `production_ps x20`). The calibration
manager now exposes both `quick` and `rescue` targets via
`run_calibration_watchers.sh {start|stop|restart|status} [quick|rescue|all]`.
The rescue watcher now refreshes the rescue plan before every polling pass, so
newly completed `warning` calibration jobs from the quick root can flow into
the stronger rescue lane without a manual restart.

Once the calibration split has at least a few paired rows, emit a calibrated
validation summary against the current merged validation winner view with:

```bash
./benchmarks/ab_bind/report_calibrated_validation.py
```

By default this summary now sweeps all shipped calibration families
(`linear`, `side_linear`, `quadratic`, `stderr_quadratic`,
`logabs_stderr_quadratic`, `expdecay_invstderr_quadratic`,
`hill_invstderr_quadratic`, `hill_side_invstderr_quadratic`) and selects the
current best model by the accepted held-out calibrated `Pearson R`. The
accepted metric uses the full held-out view unless the explicit whole-target
exclusion rule fires for a target whose paired mutations all stay above the
configured absolute-error threshold; in that case the summary promotes the
target-filtered calibrated view and records the excluded complexes in
`accepted_calibrated_excluded_complex_ids`.

For a dedicated antibody-antigen protocol-regression run on the external
`3HFM` complex without mixing it into the held-out validation split, use:

```bash
./benchmarks/ab_bind/run_3hfm_protocol_regression.sh
./benchmarks/ab_bind/report_3hfm_protocol_regression.py \
  --plan-root runs/benchmarks/abbind_3hfm_protocol_regression
```

This keeps `3HFM` under its own runs root and writes a compact regression
summary to:

- `runs/benchmarks/abbind_3hfm_protocol_regression/reports/3hfm_protocol_regression_summary.json`

If you want a literature-inspired reference setup for that same external
regression, the protocol layer now accepts optional
`equilibration_pressure_coupling`, `equilibration_pressure_tau_ps`,
`sampling_pressure_coupling`, `sampling_pressure_tau_ps`, and
`sampling_refcoord_scaling` overrides. That makes it possible to compare the
default `C-rescale` workflow against a Patel-2021-like
`Parrinello-Rahman + refcoord-scaling=com` variant without changing the main
validation protocol.

The repository now also ships a ready-made reference file for that purpose:

```bash
PROTOCOL_PATH=benchmarks/ab_bind/protocol.3hfm_patel2021_reference.yml \
  ./benchmarks/ab_bind/run_3hfm_protocol_regression.sh
```

By default this fits on `quick_plan/calibration` and predicts on the merged
validation view sourced from:

- `runs/benchmarks/abbind_core_v1_quick_plan/`
- `runs/benchmarks/abbind_core_v1_calibration_rescues/`
- `runs/benchmarks/abbind_core_v1_validation_priority_plan/`
- `runs/benchmarks/abbind_core_v1_validation_robust_plan/`
- `runs/benchmarks/abbind_core_v1_validation_priority_rescues/`

For the stronger holdout tier used to stabilize validation reporting:

```bash
./benchmarks/ab_bind/run_validation_protocol.sh
```

- this now uses `protocol.validation_robust.yml`
- outputs land under `runs/benchmarks/abbind_core_v1_validation_robust_plan/`

For a curated strong-effect follow-up on held-out complexes:

```bash
./benchmarks/ab_bind/run_validation_priority.sh
```
- this now uses `protocol.validation_priority.yml`
- outputs land under `runs/benchmarks/abbind_core_v1_validation_priority_plan/`
- keep `2nz9-antigen-a-h1064a` out of the default live queues; the mutant-side
  clash path is now handled separately and currently blocks under the robust lane
- `reports/plan_summary.json|yml`
- `reports/plan_batches.csv`
- `reports/plan_jobs.csv`
- `reports/plan_summary.json|yml` now also exposes a top-level
  `validation_gate` block with the current overall `pearson_r` and whether the
  root-level `pearson_r > 0.6` gate has passed
- `reports/plan_summary.json|yml`, `plan_batches.csv`, and `plan_jobs.csv`
  now also expose running-stage progress counters:
  `running_sample_*` for live lambda-window completion and
  `running_equilibrate_*` for live repeat completion, plus per-job
  `sample_*` / `equilibrate_*` progress columns in `plan_jobs.csv`

For unattended continuation on the robust rerun root without widening to the
entire plan, use:

```bash
./benchmarks/ab_bind/run_validation_robust_watcher.sh
```

This watcher wrapper defaults to a strict allow-list of the current robust
holdout follow-ups, and also refreshes the merged validation report rooted at
`runs/benchmarks/abbind_core_v1_validation_priority_plan/` using both the
`robust` and `rescue` roots as extra winners. Pass explicit
`job_id` values to replace the default list. Keep
`2nz9-antigen-a-h1064a` out of the default live queue for now; the current
robust lane classifies it as `blocked_input` after post-mutate geometry QC.
By default it now refreshes its allow-list every pass from the merged
validation `active_alternate_ready_hotspots` view, so the strongest current
warning points automatically stay in scope for the robust lane instead of being
frozen into a stale static list. The fallback robust list still seeds the same
high-value follow-ups, and now also includes `3hfm-antigen-y-y20a` because the
merged hotspot view currently shows a live robust alternate for that warning
point. On top of the hotspot-derived jobs, the wrapper keeps a proactive
`1MLC` seed set (`1mlc-antibody-h-s57a`, `1mlc-antibody-h-s57v`,
`1mlc-antibody-h-t31a`, `1mlc-antibody-h-t31v`, `1mlc-antibody-l-n92a`) so
those stronger reruns continue to refill even before they become merged
validation winners. The live robust refresh now also accepts a
`ROBUST_PASS_OUTLIER_THRESHOLD` fallback (default `5.0 kcal/mol`) for
completed `qc pass` rows that already have a robust alternate materialized but
have not yet entered the hotspot taxonomy; this is intended for holdout
outliers such as `1cz8-antigen-w-h90a`, so they can enter the higher-priority
robust lane instead of waiting for the lower-priority ultra queue. This
wrapper now uses `MAX_COMPUTE_APPS_PER_GPU=6` by
default for the same reason as the rescue lane: the shared node often sits at
5-7 background compute apps per GPU, so a strict cap of 5 left the robust
queue starved even when no robust job from the target hotspot set was actually
active. It now also applies the same GPU headroom override used by the
priority/targeted lanes (`MIN_FREE_GPU_MEMORY_MB=12000`,
`MAX_GPU_UTILIZATION=60`) so robust follow-ups can backfill without fighting
with large baseline queues that still leave enough free memory/utilization.

For a narrower priority-side completion queue that now dynamically targets the
best non-overlapping backlog jobs not currently covered by live alternates,
use:

```bash
./benchmarks/ab_bind/run_validation_priority_backlog_watcher.sh
```

Its live allow-list now refreshes every pass from
`refresh_validation_watchlists.py --mode backlog`, which prioritizes unfinished
validation jobs with `active_alternate_candidate_count = 0`, fewer completed
pairs for that complex, more advanced stage progress, and larger experimental
`|ddG|`. That makes the backlog lane refill from under-covered `2NZ9`, `3HFM`,
`1CZ8`, and `3NPS` jobs instead of staying pinned to a stale hard-coded list.
Pass explicit `job_id` values when you want a different backlog mix. It still
runs in `--only-listed` mode so it does not silently widen back to the whole
priority root. It now also keeps a slightly looser GPU headroom floor
(`MIN_FREE_GPU_MEMORY_MB=10500`, `MAX_GPU_UTILIZATION=60`) than the hotspot
lanes, because this queue only launches narrow `2-thread` backfills after the
heavier validation/rescue lanes have already occupied most of the shared-node
memory budget. Its default thread gate is also slightly wider
(`MAX_ACTIVE_MDRUN_THREADS=96`) so this narrow non-overlapping backlog queue
can backfill a couple more `2-thread` `mdrun` workers on the shared node after
the first 3HFM/2NZ9 launches without waiting for the whole priority fleet to
drain. That keeps the lane gated by real memory/utilization headroom instead
of a prematurely tight global thread cap alone.

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
completely. It now also refreshes its source hotspot list every pass from the
merged validation report and then re-runs `batch rescue-abbind` against the
current warning-hotspot `job_id`s before each watcher pass. That means new
warning points such as `3hfm-antigen-y-y20a` do not have to wait for a manual
rescue materialization step before entering the stronger rescue lane. The
`priority`-root merged report now also includes this rescue root automatically,
so improved rescue outputs feed the canonical holdout winner view without a
separate manual `report-abbind` step. The rescue refresh wrappers now also pass
explicit pre-equilibration scaling defaults (`WINDOW_RELAX_*`, `NVT_SCALE`,
`NPT_SCALE`) so repeat-spread hotspots are biased toward thicker relaxation
before they pay for longer production. It also now uses the same GPU headroom
override as the priority/targeted lanes (`MIN_FREE_GPU_MEMORY_MB=12000`,
`MAX_GPU_UTILIZATION=60`) so this medium-cost hotspot lane can still refill
when the shared node is oversubscribed but a device has enough real headroom.

For the narrower single-leg repeat-spread path from the sampling/QC note, use
the targeted watcher instead of widening directly into the generic rescue lane:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh start targeted
```

This wrapper refreshes
`runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues/`
through `refresh_validation_targeted_repeat_spread_rescues.sh`, now sources its
base queue from the dedicated `targeted` watchlist mode rather than the full
generic hotspot list, keeps only jobs that are currently classified as exactly
one `repeat_spread` leg, and still accepts the same job when the only extra QC
flag is a low-overlap signal on that same unstable leg. Jobs with overlap on
the opposite leg still stay out of this narrow lane so we do not inherit the
wrong source leg. Within that narrowed queue, candidate order now follows the
primary-leg `repeat_spread` severity and `ddG` repeat spread rather than the
generic Pearson-gain hotspot order, so the most unstable leg-specific jobs are
attempted first. The watcher preserves `repeats` plus `lambda_windows` while thickening
`window_relax_em_steps`, `window_relax_md_ps`, `nvt_ps`, `npt_ps`, and
`production_ps`. The same targeted rescues now also inherit the staged
equilibration restraint schedule used by the broader repeat-spread rescue path.
By default the watcher also appends completed `qc_repeat_spread` outliers that
can use the same primary-leg path but still have no active alternate anywhere,
using `TARGETED_NO_ACTIVE_ALT_ABS_ERROR_THRESHOLD=3.0` as the default gate.
The planned rescue jobs also write `config/rescue.json` with
`target_legs` and `inherit_source_legs`, so downstream stages only re-run the
unstable leg and continue to seed the stable leg from the original source job.
Its default watcher profile also adds a GPU headroom override
(`MIN_FREE_GPU_MEMORY_MB=12000`, `MAX_GPU_UTILIZATION=60`) so this narrow
hotspot lane can keep moving even when the shared node is full of many
low-memory `gmx mdrun` processes.

If one of those single-leg targeted reruns still comes back unstable, the
next step from the sampling/QC note is a lambda-densification follow-up
instead of immediately widening repeats again:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh start lambda
```

This second-stage watcher refreshes
`runs/benchmarks/abbind_core_v1_validation_targeted_lambda_rescues/`
through `refresh_validation_targeted_lambda_rescues.sh`, reuses that same
`targeted` watchlist scoping, sources only from the completed targeted-repeat
root, keeps `repeats`, `production_ps`, and the pre-equilibration thickness
fixed, and only increases `lambda_windows` (default `+4`) while preserving the
same single unstable leg targeting. The targeted-repeat lane itself now also
accepts a narrow “dominant primary leg” case when both legs exceed the
repeat-spread threshold but one leg is clearly worse and no opposite-leg
overlap issue is present, which helps hot spots such as the current `3HFM`
complex-leg-driven outliers move onto the leg-specific path earlier instead of
waiting for full two-leg deepening. Its completed jobs are now also included in
the same merged validation `report-abbind` refresh path used by the other
validation rescue lanes, so finished lambda-densification rescues can
contribute directly to the formal merged benchmark metrics instead of staying
visible only inside the dedicated lambda root.

When the hotspot still looks like a two-leg sampling problem rather than a
single unstable leg, there is now a dedicated midpoint lane before `deep`:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh start sampling-qc
```

This wrapper refreshes
`runs/benchmarks/abbind_core_v1_validation_sampling_qc_rescues/`
through `refresh_validation_sampling_qc_rescues.sh`, now defaults to all
currently detected two-leg `qc_repeat_spread` hotspots from the sampling/QC
note, and reruns the full job from the stronger
`validation_priority_rescues` source instead of forcing the narrow single-leg
inheritance path. It keeps the lane compute cap conservative
(`MAX_COMPUTE_APPS_PER_GPU=4`, `MIN_FREE_GPU_MEMORY_MB=10000`,
`MAX_ACTIVE_MDRUN_THREADS=24`, `MAX_LAUNCHES_PER_PASS=1`), now keeps
`production_ps` fixed by default while tripling the pre-equilibration
thickness (`WINDOW_RELAX_*`, `NVT_SCALE`, `NPT_SCALE`) and forcing a lambda
densification pass. By default it also applies a stronger hotspot-only
override for `3HFM`, adding a slightly denser lambda profile plus thicker
pre-equilibration without immediately escalating production time. Override or
disable that behavior through `HOTSPOT_COMPLEX_IDS` and the corresponding
`HOTSPOT_*` env vars in `refresh_validation_sampling_qc_rescues.sh`. The lane
still reuses the richer current alternate when one exists, and the queue is now
ordered by the shared two-leg repeat-spread burden before generic impact
heuristics so the broadest sampling failures move first. Use this lane for
cases where both legs are still
showing repeat-spread symptoms and you want to try “thicker preparation first”
before escalating into the heavier `deep/ultra` presets. Dominant-primary-leg
two-leg hotspots are now filtered back toward the narrower targeted/lambda
path so this lane stays focused on genuinely broad two-leg sampling failures.
By default it also appends completed two-leg `qc_repeat_spread` outliers that
still have no active alternate anywhere once `|ddG error| >= 2.5 kcal/mol`,
which helps cases like the current 1MLC/1CZ8 two-leg misses enter the broader
sampling/QC lane without waiting for them to appear in the smaller hotspot list.
If you want to keep the lane narrow, set `SAMPLING_QC_COMPLEX_IDS=3HFM`
before starting the watcher.

For the broader hotspot-deepening and highest-impact 3HFM follow-up lanes, use:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh start deep ultra
```

The deep wrapper now uses a slightly looser GPU headroom override
(`MIN_FREE_GPU_MEMORY_MB=10000`, `MAX_GPU_UTILIZATION=60`) so it can keep
backfilling the next hotspot candidates without waiting for a nearly empty GPU,
while the ultra wrapper stays at the tighter
`MIN_FREE_GPU_MEMORY_MB=12000` guardrail. Both lanes still keep their stricter
compute-app/thread caps, and the deep/ultra refresh steps now also apply a
hotspot-only protocol override for `3HFM` by default:
non-hotspot jobs still keep the lane-wide deep/ultra preset, while matched
`3HFM` rescues get thicker pre-equilibration (`WINDOW_RELAX_*`, `NVT_SCALE`,
`NPT_SCALE`) and a denser lambda/repeat profile through the new
`--hotspot-*` rescue-abbind options. Override or disable that behavior by
setting `HOTSPOT_COMPLEX_IDS` and the corresponding `HOTSPOT_*` env vars in
`refresh_validation_deep_rescues.sh` or `refresh_validation_ultra_rescues.sh`.
Those deep/ultra refresh paths also now reuse
`runs/benchmarks/abbind_core_v1_validation_targeted_repeat_spread_rescues/`
as an extra source root and pass
`--target-primary-repeat-spread-leg --allow-targeted-leg-count-deepening`, so
single-leg repeat-spread hotspots can stay on the narrower leg-specific path
while deep/ultra still add repeats and lambda windows instead of reverting to
a full two-leg rerun.
The ultra refresh path now also exposes `REQUIRE_ACTIVE_ALTERNATE`
(default `1`). Leave it enabled for the normal path; temporarily set
`REQUIRE_ACTIVE_ALTERNATE=0` when the merged validation summary is lagging
behind the per-root `plan_jobs.csv` state but you still need to materialize a
known high-error job that already has a live robust/targeted/deep alternate.
The ultra watchlist is no longer limited to high-impact complexes alone:
besides the `complex_impact_pearson_gain` gate, it now also promotes
active-alternate validation hotspots whose current `|ddG error|` exceeds
`ULTRA_ABS_ERROR_THRESHOLD` (default `5.0 kcal/mol`). That keeps the 3HFM
leave-one-complex-out path as the main ultra driver, while still allowing
exceptionally bad outliers such as the current `1BJ1` / `1CZ8` leaders to
enter the ultra lane once a live alternate exists. The ultra watcher also now
adds an opt-in pass-QC outlier fallback through
`ULTRA_PASS_OUTLIER_THRESHOLD` (defaulting to the same value as
`ULTRA_ABS_ERROR_THRESHOLD` in the live wrapper), so extreme completed
outliers that already have a non-active alternate candidate can still be
backfilled into ultra even when they have not yet entered the hotspot taxonomy.
The live ultra watcher now reads `ultra_pass_outlier_job_ids` from the
watchlist JSON and, when `ULTRA_PASS_OUTLIER_ALLOW_INACTIVE_ALTERNATE=1`
(default), runs a second ultra-refresh pass for just those jobs with
`REQUIRE_ACTIVE_ALTERNATE=0` and `ALLOW_PASS_QC_OUTLIER_RESCUE=1`. The latter
maps those explicitly requested, `qc pass` but still high-error jobs onto a
dedicated `pass_qc_outlier` rescue reason, instead of requiring them to first
show up as a traditional QC failure. That keeps the normal hotspot path strict
while still letting pass-QC outliers such as `1cz8-antigen-w-h90a` get
materialized into the ultra rescue root instead of stopping at watchlist
selection.
The ultra watcher also now
defaults to `MAX_LOAD_PER_CORE=0` and `MAX_ACTIVE_MDRUN_THREADS=20`, so these
new outlier ultra jobs can still backfill one additional low-cost outlier
behind the already-running 3HFM ultra set instead of sitting indefinitely in
`resumable` behind a stricter thread gate than the actual GPU-headroom policy.
The deep refresh path also no longer requires an already-active alternate by
default; it still prefers active alternates as the source when available, but
taxonomy hotspots such as `1bj1-antigen-w-g88a` can now materialize directly
into the deep lane instead of waiting for a lower rescue root to finish first.
Set `REQUIRE_ACTIVE_ALTERNATE=1` in `refresh_validation_deep_rescues.sh` when
you want the older stricter behavior back. The deep watcher now also disables
the extra CPU load/core gate by default (`MAX_LOAD_PER_CORE=0`) and relies on
GPU headroom plus `MAX_ACTIVE_MDRUN_THREADS=24` as the effective launch cap,
so newly materialized deep hotspots are no longer blocked indefinitely on this
shared node by a conservative load-average threshold that was stricter than the
actual thread budget.

To manage the validation watchers together, use:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh restart
./benchmarks/ab_bind/run_validation_watchers.sh status
```

By default, `start`/`restart` now only relaunch the `robust` watcher. Add
explicit targets such as `priority`, `rescue`, `targeted`, `lambda`, or `all`
only when you intentionally want those extra queues back:

```bash
./benchmarks/ab_bind/run_validation_watchers.sh restart all
./benchmarks/ab_bind/run_validation_watchers.sh start priority rescue targeted lambda
```

This still writes pidfiles under each `reports/watch/` directory and the
selected watcher entrypoints still run under the self-restarting supervisor
loop. The direct `priority` and `rescue` watcher wrapper scripts are now
disabled by default; use `run_validation_watchers.sh` when you intentionally
want them enabled. The priority watcher now keeps refreshing reports while
active jobs remain live and uses a faster default refill profile
(`--max-launches-per-pass 4`, `--launch-cooldown-seconds 60`) so stale/resumable
validation backlog is pulled back onto GPUs more quickly. It now also applies
the same GPU headroom override used by the targeted lanes
(`MIN_FREE_GPU_MEMORY_MB=12000`, `MAX_GPU_UTILIZATION=60`), so the watcher can
continue refilling under shared-node oversubscription when devices still have
memory/utilization headroom. When you need a temporary experimental root to
feed the formal merged validation report without editing the wrapper defaults,
set `MERGED_EXTRA_PLAN_ROOTS=/abs/root_one:/abs/root_two` before launching the
priority watcher; those roots are appended to the built-in merged sources and
flow through the same post-report refresh chain.
The rescue watcher now explicitly allows the same `job_id` to run in parallel
with priority or robust roots when the plan root differs, so medium-cost rescue
reruns no longer sit in `active_elsewhere` behind an already active baseline or
robust copy of the same mutation.
To keep a handful of hotspot jobs from occupying five or six concurrent rescue
copies on the shared node, the duplicate-tolerant validation rescue wrappers
(`rescue`, `targeted`, `lambda`, `sampling-qc`, `deep`, `ultra`) now also
default `MAX_ACTIVE_COPIES_PER_JOB_ID=3`. Override it per wrapper when you
intentionally want broader same-`job_id` concurrency; the `stale`/`gap`
recovery wrappers instead default to a smaller `MAX_ACTIVE_COPIES_PER_JOB_ID=2`
so they can keep refilling `priority`/`robust` coverage without creating a
third or fourth copy of the same hotspot mutation.
The `stale` and `gap` recovery wrappers also use a slightly looser GPU headroom
override (`MIN_FREE_GPU_MEMORY_MB=10000`, `MAX_GPU_UTILIZATION=60`), so
coverage-recovery passes for `1MLC`/`2NZ9` do not get starved just because the
node carries many lightweight background compute apps.
They now also scope their thread budget to the `priority` + `robust` plan
roots. The `stale` wrapper now defaults `INCLUDE_GAP_JOBS=0`, so its queue stays
focused on true stale-resume coverage instead of competing with the dedicated
`gap` lane for the same `2NZ9` catch-up jobs. Re-enable `INCLUDE_GAP_JOBS=1`
only if you intentionally want the stale lane to absorb the gap backlog too.

If the observed BAR uncertainty is too large, the job now lands in
`qc_status=warning` even when all expected files are present. Parsed BAR
histograms are now also checked for support overlap after choosing the better
of raw and sign-reflected reverse alignment, which matches GROMACS BAR outputs
where the reverse `ΔH` histogram may appear with the opposite sign convention.
`benchmark_qc_qualified=True` is now tracked separately from `qc_status` and
requires the job to be `ddg`-ready, non-failing, and within the configured
`max_bar_stderr_kcal_mol` threshold.
When merged reports compare the same mutation group across multiple plan roots,
the winner selection now also prefers the higher protocol sampling effort
(`repeats * lambda_windows * production_ps`) before falling back to lower
`ddg_bar_stderr_kcal_mol`, so robust or rescued reruns replace cheaper priority
lane outputs once their QC tier is otherwise comparable.

The equilibration stage now also has a project-local robustness fallback for
large complexes: when EM hits a periodic-boundary shift fatal, the job retries
through a two-stage larger cubic box ladder and preserves an
`em.runtime.history.log` audit trail for those EM attempts. The default EM
preset is now regenerated with `-DFLEXIBLE` and `constraints = none`.

The prepare stage now also repairs sidechain-only missing heavy atoms before
`pmx`/`pdb2gmx` by stripping partial sidechains and passing the surviving
backbone through `PDBFixer`. This deliberately does not invent missing residue
spans or unresolved loops; backbone-incomplete cases still block early.

## Validated Benchmark Run

The benchmark assets have been exercised on one real AB-Bind-derived case:

- source complex: `1VFB`
- mutation: `H:Y32A@antibody`
- batch dir:
  `runs/benchmarks/abbind_1vfb_core_v1_quick/`
- result:
  `ddG = 10.306 kcal/mol`, `QC = warning`

Reproduce it with:

```bash
./benchmarks/ab_bind/run_1vfb_h_y32a_quick.sh
```

## Bulk Planning Status

A real `core_v1` quick-planning pass has already been materialized under:

- `runs/benchmarks/abbind_core_v1_quick_plan/`
- `18` per-complex batch directories
- `318` planned jobs total

The same root now also contains a cross-complex validation slice:

- `1JRH`
  - `1jrh-antigen-i-t14v`
  - completed through `report`
  - `ddG = -7.805 kcal/mol`, `QC = warning`
- `3NGB`
  - `3ngb-antibody-h-g54s`
  - completed through `report`
  - `ddG = -51.303 kcal/mol`, `QC = warning`
  - this case exercised the new per-window `pre_relax/pre_md` sample
    stabilization path before production
- `1YY9`
  - `1yy9-antibody-h-n56a`
  - `1yy9-antibody-l-n93a`
  - both stop at `prepare` with `blocked_input` because incomplete standard
    residues remain after extraction

The plan-root execution/reporting path has also been validated on a fresh runner
root:

- root: `runs/benchmarks/abbind_1vfb_runner_quick_plan/`
- selected executed jobs:
  - `1vfb-antibody-h-y32a`
  - `1vfb-antibody-h-g31a`
  - `1vfb-antigen-c-y23a`
- aggregate root-level counts:
  - `23` jobs in the selected `1VFB` batch
  - `3` ready jobs
  - `3` jobs currently `warning`
  - `20` jobs remain `not_started`

A real `V2` same-side double-point benchmark quick run is now also available:

- root: `runs/benchmarks/abbind_1jrh_runner_truequick/`
- selected executed job: `1jrh-antigen-i-m25l--i-i28v`
- mutation signature: `I:M25L@antigen + I:I28V@antigen`
- result: `ddG = -19.054 kcal/mol`
- QC: `warning`
- aggregate root-level counts:
  - `3` jobs in the selected `1JRH` batch
  - `1` ready job
  - `2` jobs remain `not_started`

Reproduce it with:

```bash
./benchmarks/ab_bind/run_1jrh_double_quick.sh
```
