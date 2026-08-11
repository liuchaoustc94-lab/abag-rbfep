# Repository Guidelines

## Project Structure & Module Organization

`src/abag_rbfe/` contains the CLI, workflow stages, GROMACS integration, planning, QC, and reporting. Mutation helpers live in `src/abag_pmx/`; upstream `pmx` is pinned under `vendor/pmx/`. Put tests in `tests/`. Small runnable inputs belong in `examples/`; benchmark protocols, manifests, and scripts belong in `benchmarks/`. Documentation lives in `docs/`. Treat `runs/` and `tmp/` as generated workspaces.

## Build, Test, and Development Commands

Create the supported Python 3.13 environment and editable installs with:

```bash
conda env create -f envs/environment.yml
conda activate abag-rbfep
```

For an existing environment, mirror CI explicitly:

```bash
python -m pip install -e ./vendor/pmx
python -m pip install -e '.[dev]'
python -m pytest
```

Run a focused test during development with `python -m pytest tests/test_planning.py -q`. Exercise the CLI without launching simulations using commands such as `abag-rbfe mutation validate --mutations examples/mutations.csv --output /tmp/mutations.json`. Add `--execute` only when a real GROMACS run is intended.

## Coding Style & Naming Conventions

Use four-space indentation, PEP 8 layout, type hints on public interfaces, and `pathlib.Path` for filesystem work. Use `snake_case` for modules, functions, variables, and tests; use `CapWords` for classes. Keep stage and report payloads backward-compatible because existing runs are resumable. No formatter or linter is configured, so match nearby code.

## Testing Guidelines

Pytest is the required framework; CI runs the complete suite. Name files `test_<area>.py` and tests `test_<behavior>`. Use `tmp_path`, monkeypatching, and mocked command runners for deterministic tests; do not require GPUs or long molecular-dynamics jobs in unit tests. Changes to report wording, JSON fields, watcher behavior, or resume logic must update the corresponding regression tests.

## Commit & Pull Request Guidelines

This checkout contains no usable Git history, so no repository-specific commit pattern can be verified. Use short imperative subjects, optionally scoped, for example `reporting: preserve QC diagnostics`. Keep commits focused. Pull requests should explain the scientific or workflow impact, list tests run, identify changed CLI or artifact schemas, and link the relevant issue. Include compact before/after report excerpts when outputs change; never commit generated trajectories, `runs/`, caches, or local benchmark source data.

## Configuration & Safety

Keep the project isolated from older platform repositories. Preserve the vendored `pmx` boundary, use project-relative configuration, and review generated shell commands before external execution. Do not overwrite canonical benchmark reports with filtered selections; selection-specific output belongs under `reports/selections/`.
