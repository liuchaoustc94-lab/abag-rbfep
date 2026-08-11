# Reference Materials

This note collects papers and official documentation that are worth reading
before extending `abag-rbfep`. It is intentionally biased toward the current
project boundary:

- `V1`: antibody- or antigen-side single-point ddG
- `V2`: same-side double-point ddG
- `GROMACS + pmx` as the execution stack
- `gmx bar` as the current formal result source

This is not a general free-energy bibliography. It is an implementation-facing
reading list for the current repository state.

For a shorter Chinese reading map tied to the same project boundary, see
`docs/reference_materials_cn.md`.

## Suggested Reading Order

If you only read a short list, use this order:

1. GROMACS free-energy documentation
2. GROMACS `gmx bar` manual
3. GROMACS free-energy interaction functions
4. Gapsys et al. `pmx` paper
5. pmx tutorials for protein mutation and analysis
6. AB-Bind paper
7. Patel et al. protein-protein alchemical assessment
8. Klimovich et al. analysis guidelines
9. Jia et al. DC-MBAR paper
10. mmCSM-AB paper

## Reading Tracks By Milestone

### V1: single-point, charge-conserving

Read these before touching protocol defaults or benchmark conclusions:

1. GROMACS free-energy docs
2. `gmx bar` manual
3. pmx paper
4. pmx mutation tutorial
5. AB-Bind paper
6. Klimovich et al.

### V2: same-side double-point

Read these before changing the mutation data model or the default lambda
schedule for double-point runs:

1. Re-read the pmx paper with the "one or multiple mutations" setup in mind
2. Patel et al. for protein-protein alchemical execution realities
3. mmCSM-AB for the antibody multiple-mutation benchmark framing
4. Jia et al. DC-MBAR as an analysis-scale reference, not as a default
   replacement for BAR

### V2.1: charge-changing and cross-side double-point

Read these before opening the charge-changing scope:

1. pmx scripts/docs around `doublebox`
2. GROMACS free-energy interaction functions
3. Klimovich et al. overlap / uncertainty guidance

This is the point where BAR-vs-MBAR sidecar analysis becomes more interesting,
but it still should not displace the formal BAR path without dedicated
cross-validation.

## User-Provided Paper: DC-MBAR

### Jia, Ge, Mei (2021)

- Title: `Free energy change estimation: The Divide and Conquer MBAR method`
- Journal: `Journal of Computational Chemistry`
- DOI: `10.1002/jcc.26533`

Why this paper matters for `abag-rbfep`:

- It proposes a scalable alternative to full MBAR when the number of
  alchemical states becomes large.
- It starts from pairwise overlap estimates, defines adjacent states with an
  overlap threshold, solves smaller local MBAR problems, and reconstructs the
  total free-energy profile from adjacent increments.
- The pairwise decomposition is naturally parallel, so it maps well to the
  existing multi-job / multi-device orchestration style of this project.

Why it should not replace the current default analysis path yet:

- The repository currently treats `gmx bar` as the formal source of truth for
  reported ddG values, so any MBAR-based alternative should start as an
  optional sidecar analysis path.
- The paper itself notes a real practical weakness: the overlap threshold is
  system-dependent.
- The published examples are umbrella-sampling systems, not antibody-antigen
  RBFE on pmx-generated hybrid protein topologies.

Best use in this project:

- future post-hoc reanalysis of DHDL data for dense lambda schedules
- future QC diagnostics for low-overlap lambda windows
- future large-window `double_point` protocols where full multistate analysis
  may become a bottleneck

Questions to ask before investing implementation time here:

- do we already have enough completed lambda windows for MBAR to add value
- is overlap failure actually the dominant error source, or is sampling still
  the main problem
- can the sidecar analysis be benchmarked against the existing `gmx bar` path
  on the same completed jobs before it is exposed to users

## Core Theory

### Bennett (1976)

- Citation: Bennett CH. `Efficient estimation of free energy differences from
  Monte Carlo data.` J Comput Phys. 1976.
- Why read it:
  - foundation of BAR
  - overlap intuition still matters when diagnosing bad lambda spacing
  - explains why weak overlap turns directly into noisy ddG estimates

### Shirts and Chodera (2008)

- Citation: Shirts MR, Chodera JD. `Statistically optimal analysis of samples
  from multiple equilibrium states.` J Chem Phys. 2008.
- DOI: `10.1063/1.2978177`
- Why read it:
  - formal MBAR reference
  - useful when deciding whether a future `abag-rbfep` analysis layer should
    go beyond BAR
  - gives the right conceptual frame for combining multiple lambda states and
    for uncertainty handling

### Klimovich, Shirts, Mobley (2015)

- Citation: Klimovich PV, Shirts MR, Mobley DL. `Guidelines for the analysis of
  free energy calculations.` J Comput Aided Mol Des. 2015.
- DOI: `10.1007/s10822-015-9840-9`
- Why read it:
  - practical best-practice paper for overlap checks, uncertainty assessment,
    convergence inspection, and analysis hygiene
  - directly relevant to `qc` and `report` stage design

## Stack-Specific References

### GROMACS official documentation

- Free-energy algorithms:
  `https://manual.gromacs.org/current/reference-manual/algorithms/free-energy-calculations.html`
- `gmx bar` command reference:
  `https://manual.gromacs.org/current/onlinehelp/gmx-bar.html`
- Free-energy interaction functions:
  `https://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html`

Why these matter:

- They define the engine-level assumptions behind lambda states, DHDL output,
  and BAR aggregation.
- `gmx bar` already supports combining adjacent free-energy contributions into a
  total estimate, which is exactly why it remains the formal analysis path in
  this repository.
- The interaction-functions section is especially useful when adjusting
  alchemical soft-core, bonded interpolation, or other lambda-dependent
  behavior that can quietly destabilize large protein mutation jobs.

### Gapsys et al. (2015) pmx paper

- Citation: Gapsys V, Michielssens S, Seeliger D, de Groot BL. `pmx:
  Automated protein structure and topology generation for alchemical
  perturbations.` J Comput Chem. 2015.
- DOI: `10.1002/jcc.23804`
- Why read it:
  - this is the scientific core behind hybrid residue and hybrid topology
    generation
  - it explains the exact class of mutation problems `pmx` was built to solve
  - it is the right paper to cite when describing the mutation/topology part of
    the `abag-rbfep` workflow

### pmx official docs and tutorials

- Project docs: `https://degrootlab.github.io/pmx/`
- Tutorials index: `https://degrootlab.github.io/pmx/tutorials/index.html`
- Protein mutation tutorial:
  `https://degrootlab.github.io/pmx/tutorials/protein_mut.html`
- Analysis example:
  `https://degrootlab.github.io/pmx/examples/analysis.html`
- Scripts index:
  `https://degrootlab.github.io/pmx/scripts/index.html`

Why these matter:

- They are the most direct reference for the current `vendor/pmx` execution
  model.
- They provide concrete command-level expectations for mutation setup,
  topology generation, and analysis handoff.
- The analysis example is especially useful because it explicitly points to
  `gmx bar` for equilibrium free-energy estimation.
- The scripts index is worth bookmarking because it exposes forward-looking
  tools like `doublebox`, which are directly relevant to the planned
  `V2.1` charge-changing path.

## Antibody-Antigen Problem Framing

### AB-Bind

- Citation: Sirin S, Apgar JR, Bennett EM, Keating AE. `AB-Bind: Antibody
  binding mutational database for computational affinity predictions.` Protein
  Sci. 2016.
- DOI: `10.1002/pro.2829`
- Why read it:
  - source benchmark for this project
  - defines the mutation-data landscape and why structure-mapped antibody-
    antigen ddG prediction is hard
  - useful when explaining why curated RBFE-ready subsets are smaller than the
    full literature set

### AB-Bind repository

- Repository: `https://github.com/sarahsirin/AB-Bind-Database`
- Why read it:
  - practical source for understanding how the benchmark was distributed
  - useful when reproducing curation decisions outside the paper text
  - helpful for tracking field names and original identifiers during data
    ingest and audit

### Patel et al. (2021)

- Citation: Patel D, Patel JS, Ytreberg FM. `Implementing and Assessing an
  Alchemical Method for Calculating Protein-Protein Binding Free Energy.` J
  Chem Theory Comput. 2021.
- DOI: `10.1021/acs.jctc.0c01045`
- Why read it:
  - directly about mutation-induced protein-protein binding free-energy changes
  - uses `pmx` for hybrid structures/topologies
  - highlights the engineering pain points that still show up in practice:
    hybrid setup, neutral net charge handling, and large-complex execution cost
  - the closest paper in this list to the operational problems seen in this
    repository during real antibody-antigen runs

### mmCSM-AB

- Citation: Myung Y, Pires DEV, Ascher DB. `mmCSM-AB: guiding rational antibody
  engineering through multiple point mutations.` Nucleic Acids Res. 2020.
- DOI: `10.1093/nar/gkaa406`
- Why read it:
  - not for runtime reuse
  - useful as an external multiple-mutation reference when framing `V2`
    benchmark design and report format
  - confirms that same-side multiple-point mutation assessment is a real
    antibody-engineering use case, not a synthetic extension

## Optional Companion References

These are not required to keep the current BAR-first path moving, but they are
good design inputs for future modules.

### biobb_pmx

- Project: `https://github.com/bioexcel/biobb_pmx`
- Why read it:
  - useful workflow-layer reference for wrapping pmx steps in reproducible
    Python building blocks
  - a good place to borrow ideas for idempotent stage interfaces and metadata
    capture without inheriting a runtime dependency

### Shirts and Chodera (2008)

- Citation: Shirts MR, Chodera JD. `Statistically optimal analysis of samples
  from multiple equilibrium states.` J Chem Phys. 2008.
- DOI: `10.1063/1.2978177`
- Why keep it nearby:
  - formal MBAR foundation for any future sidecar analysis path
  - useful when reasoning about uncertainty propagation beyond adjacent BAR
    pairs

## How To Map The Reading Back To This Repo

- `vendor/pmx/` and `src/abag_pmx/`
  - read Gapsys et al. plus the pmx docs
- `src/abag_rbfe/stages.py`, `src/abag_rbfe/gmx.py`
  - read the GROMACS free-energy docs, interaction functions, and `gmx bar`
    manual
- `src/abag_rbfe/reporting.py`
  - read Bennett, Shirts-Chodera, and Klimovich et al.
- `benchmarks/ab_bind/`
  - read AB-Bind first, then the AB-Bind repository, then mmCSM-AB for
    multi-point framing
- future optional MBAR sidecar
  - read Jia et al. DC-MBAR after the BAR path is already stable

## How The New Survey Should Influence This Repo

The additional survey is useful, but only part of it should drive immediate
repository changes.

### Already aligned with the current implementation

- It reinforces the existing project boundary: `GROMACS + pmx` for setup and
  `gmx bar` as the formal result source.
- It supports the current decision to keep explicit `NVT -> NPT`
  equilibration instead of relying on a thinner pmx-style setup path.
- The current workflow still uses `C-rescale` by default in generated
  equilibration and sampling `mdp` files, but Patel 2021 should not be cited as
  direct support for that choice. The main text describes
  `Parrinello-Rahman`, while the supporting `mdp` bundle mixes
  `Berendsen/Parrinello-Rahman`, so this paper is better treated as an
  external regression reference than as a literal default-protocol template.
- It is consistent with the current roles of `AB-Bind`, `mmCSM-AB`, and
  `DC-MBAR` in this repository: benchmark source, `V2` framing reference, and
  future analysis-sidecar input respectively.

### Best next-use items

- `Patel et al. 2021` makes `3HFM` a strong external calibration candidate for
  protocol regression outside the main AB-Bind plan roots.
- `Sampson et al. 2024` is useful for designing a clearer validation failure
  taxonomy around outliers and hard protein-protein cases.
- `Clark 2017/2019` are more useful as design constraints for `V2.1`
  charge-changing scope than as immediate `V1/V2` default-protocol inputs.
- `biobb_pmx` remains a workflow-pattern reference, not a runtime dependency.

### Useful, but not ready to become defaults

- Literature-scale production settings such as long equilibrations, multiple
  replicas, and large snapshot counts are better treated as upper-bound design
  inputs for robust/deep presets, not immediate validation defaults.
- Enhanced-sampling routes such as `REST2`, `GaMD`, and lambda-exchange are
  better treated as rescue paths for difficult systems.
- `DC-MBAR` should remain an optional sidecar analysis path until it is
  benchmarked head-to-head against the existing BAR pipeline on completed jobs.

### Immediate engineering actions suggested by the survey

1. Add a `3HFM` external calibration run and archive its outputs.
2. Add a clearer validation failure taxonomy to reports.
3. Write the `V2.1` charge-changing design constraints down explicitly so they
   do not leak into `V2` same-side double-point work.

## Practical Recommendation For The Current Project State

For the next implementation cycle, the best reading sequence is:

1. keep `gmx bar` and the GROMACS free-energy manual as the operational base
2. use the `pmx` paper and tutorials to guide hybrid-topology robustness work
3. use AB-Bind as the benchmark ground truth source
4. use Klimovich et al. to tighten QC and reporting criteria
5. keep DC-MBAR as a design input for a later optional analysis module, not as
   the first-line reported result path
6. use Patel et al. and mmCSM-AB when extending from `V1` to `V2`, because
   they are the most relevant references for large-complex protein-protein
   execution cost and antibody multiple-mutation framing
