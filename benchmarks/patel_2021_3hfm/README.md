# Patel 2021 3HFM External Regression Set

This directory captures the eight single-point `3HFM` mutations reported in:

- Patel, Patel, Ytreberg. `Implementing and Assessing an Alchemical Method for Calculating Protein-Protein Binding Free Energy`. JCTC 2021.

It is intentionally separate from the main `AB-Bind` benchmark flow.

Use it for:

- protocol regression against a published antibody-antigen FEP case
- charge-changing single-point experiments outside the current held-out validation lane
- comparing the default equilibrium BAR workflow with a Patel-inspired reference protocol

## Files

- `system.yml`: `3HFM` system definition
- `mutations.csv`: eight literature mutations
- `experimental_ddg.csv`: experimental `ddG` values from Table 2

## Chain mapping notes

The paper table uses residue-style labels such as `W98F` and `D101K` without
chain letters. Most rows map directly from unique residue identity in the
`3HFM` structure. Two rows are inferred:

- `H:W98F`: `TRP 98` is unique to heavy chain `H`
- `Y:D101K`: both heavy chain `H` and antigen chain `Y` contain `ASP 101`, but
  antigen `Y:101` is the interface-exposed candidate consistent with the paper's
  reported mutation set

## Run

```bash
./benchmarks/patel_2021_3hfm/run_patellike_3hfm.sh
./benchmarks/patel_2021_3hfm/report_patellike_3hfm.py \
  --batch-dir runs/benchmarks/patel_2021_3hfm/patel_2021_3hfm_reference
```

For sustained execution on a shared GPU host, use the watcher instead of the
one-shot runner:

```bash
./benchmarks/patel_2021_3hfm/run_patellike_3hfm_watcher.sh
./benchmarks/patel_2021_3hfm/run_patellike_3hfm_watcher.sh 3hfm-patel-2021-antigen-y-y20f
```

Watcher defaults:

- refreshes `batch_summary.json` and `patel_2021_3hfm_summary.json` every pass
- launches up to two new Patel jobs per pass
- defaults to `SKIP_CHARGE_CHANGING=1`, so the charge-conserving literature subset
  is still launched first; once every conserving row is already active,
  recently launched, analyzable, or completed, the watcher can use spare launch
  capacity for the deferred charge-changing rows instead of idling free GPUs
- defaults to `MAX_COMPUTE_APPS_PER_GPU=13` and `MDRUN_ARGS_OVERRIDE='-ntmpi 1 -ntomp 2'`
  to avoid competing too aggressively with the main AB-Bind validation lanes
- now also applies the same GPU headroom override used by the validation
  hotspot lanes (`MIN_FREE_GPU_MEMORY_MB=12000`, `MAX_GPU_UTILIZATION=60`), so
  the watcher can still prefer devices with usable memory/utilization headroom
  even when the shared node is saturated by many low-memory `mdrun` processes

The Patel summary report now breaks out `charge_conserving` versus
`charge_changing` rows explicitly. When the default watcher is left at
`SKIP_CHARGE_CHANGING=1`, the summary should therefore show the
charge-conserving subset as launching/running first, while the charge-changing
rows remain deferred until there is no conserving launch backlog.
