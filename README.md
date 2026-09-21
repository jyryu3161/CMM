# CMM — Constraint-based Metabolic Modeling

[![CI](https://github.com/jyryu3161/CMM/actions/workflows/ci.yml/badge.svg)](https://github.com/jyryu3161/CMM/actions/workflows/ci.yml)

CMM is a solver-aware Python platform for reproducible constraint-based metabolic engineering.
It proposes and tests genetic interventions for target-metabolite production, integrates
expression data, analyses perturbation responses, and records the evidence needed to reproduce
a computational study. The optional Qt desktop application calls the same numerical services;
Python and the thin CLI are the publication workflow boundaries.

> **Research software status:** CMM 0.5.0 is beta software being prepared for journal
> publication. Passing the validation suite establishes implementation and artifact
> integrity, not biological validation. Treat every predicted intervention as an *in silico*
> hypothesis until it is tested experimentally.

![CMM platform](docs/images/overview.png)

## Reproducible production workflow

The canonical production workflow keeps proposal methods and forward checks distinct:

| Scientific question | Methods | Evidence reported |
|---|---|---|
| Which single-gene deletions change the phenotype? | MOMA-L2 and ROOM | Separate growth-versus-product screens, GPR-resolved reaction information, and flux response plus paired sampling for every unique D1–D5 candidate |
| Which deletion sets couple production to growth? | OptKnock and RobustKnock | Deterministically seeded MILP searches; maximum and guaranteed product; designs ranked by the guarantee |
| Which fluxes are amplification hypotheses? | FSEOF and FVSEOF | Independent top-10 method rankings and trajectories, visible loop diagnostics, and a flux-response record for every candidate in their union |
| Can another researcher audit the claims? | Provenance, manifest, R renderer, run validator | Exact model/config, raw CSVs, metadata, 300-DPI PNG, editable PDF/SVG, linked and standalone HTML |

One UTF-8 JSON config specifies the exact SBML file, product exchange, medium, substrate
uptake, oxygen bounds, solver, search limits, the strain-design seed, and sampling seed. A
complete run is:

```bash
uv run --frozen --all-extras python examples/production-targets/prepare.py
uv run --frozen --all-extras cmm production-targets --config examples/production-targets/config.json
uv run --frozen --all-extras cmm report validate results/example-production-succinate --json
```

This [checked-in succinate example](examples/production-targets/README.md) first saves the
existing COBRApy textbook model. The `production-targets` command then performs preflight,
numerical analysis, publication rendering, and final
validation. Use `--analysis-only` to separate the solver run from R rendering. A run is not
complete until validation succeeds; unavailable methods and infeasible solutions remain
visible rather than being silently omitted.

See [Building or customizing a CMM workflow](docs/building-custom-workflows.md) for a complete
generic config, the equivalent Python API, and the boundary between a downstream study and an
installed workflow. Contributors can follow the separate
[canonical-workflow tutorial](docs/tutorials/adding-a-canonical-workflow.md). The scientific
sequence of the production workflow is specified in
[SC-01 production target discovery](docs/scenarios/SC-01-production-target-discovery.md).

## Reproducible transformation workflow (MTA/rMTA)

The second canonical workflow ranks gene or reaction knockouts that move a **source**
metabolic state toward a **target** state. It uses published MTA or rMTA, reports the MOMA
baseline, and can repeat the ranking over an explicit epsilon sensitivity grid.

```bash
uv run --frozen --all-extras python examples/transformation-targets/prepare.py
uv run --frozen --all-extras cmm transformation-targets --config examples/transformation-targets/config.json
uv run --frozen --all-extras cmm report validate results/example-transformation-mta --json
```

This [checked-in MTA/rMTA example](examples/transformation-targets/README.md) exports CMM's
existing synthetic disease → healthy fixture. It also explains how to provide Recon inputs;
the historical Recon1 timing below has no bundled input pair or run config. Both runnable
workflows are indexed in [examples](examples/README.md).

For a new study, set the exact model, source/target expression files, condition and epsilon.
Expression inputs are finite, non-negative **linear**
measurements: reference inference uses linear replicate means, while direction tests compare
`log2(value + 1)` per replicate. Use the t-test for replicated measurements and explicitly
select fold-change thresholds for single measurements. Confirm source → target before running.

Both methods require an MIQP-capable Gurobi or CPLEX installation. The default reference is
E-Flux2, not the original paper's iMAT-plus-sampling reference, and epsilon must be chosen for
the data; these departures are recorded in the report. See the
[SC-02 contract](docs/scenarios/SC-02-transformation-target-discovery.md) and
[config/API reference](docs/agent-reference.md#complete-transformation-workflow--cmmworkflowstransformation).

### A decision model plays the model

Both workflows above are deterministic optimisations. The third entry point is not: it puts
TypeSafe's JEV decision model in a loop and lets it search the design space, one move at a
time, while CMM enforces the rules.

```bash
export OPENROUTER_API_KEY=...   # the only part of CMM that needs a credential
uv run --frozen --all-extras cmm jev-design --config examples/jev-design/config.json
```

JEV is a *System One* model: it generates no text and returns a typed answer chosen from
criteria CMM supplied, so it cannot name a reaction the model does not contain and there is no
output to parse. Each tick, CMM renders the metabolic state as compact JSON — fluxes, distance
to the product, ATP and redox balance, what each active intervention actually bought — and JEV
picks one move: a knockout, a knockdown, an amplification, switching an unused reaction on,
withdrawing an earlier move, running an FSEOF or essentiality scan, or ending the round. CMM
applies it, re-solves with pFBA and MOMA, and redraws the flux map. In the desktop
application's *JEV Agent* tab the map moves while the agent is still playing.

**CMM owns the rules and JEV owns the strategy.** A move that makes the model infeasible or
pushes growth below the configured floor is reverted by CMM whatever the agent predicted, and
the reversal is a recorded row. A move that is merely unhelpful is kept for the agent to
withdraw itself.

**Every run scores itself against the deterministic methods**, on the same model and the
same growth floor, every design applied and solved the same way. On the shipped anaerobic
succinate example the best single gene deletion of 71 reaches 0.211 mmol gDW⁻¹ h⁻¹, OptKnock
and RobustKnock reach 9.911, and the agent reaches 10.761 — but only because it is handed
OptKnock's own proven design and adds the one move OptKnock cannot express, forcing flux
through the glyoxylate shunt, at a real cost in growth. Left to search on its own it does not
match 9.911 at all: the winning design deletes reactions carrying no flux in the wild type,
and a board built from where the flux is today cannot see them. Over ten runs of one
configuration, nine reached 10.761 and one stopped at 9.911, for $0.02 and 66 seconds in
total.

Measured on `e_coli_core`: two calls per move, about 0.6 s and $0.00016 each, so a full run
costs a fraction of a cent. **The CMM solves in a run are deterministic; JEV's decisions are
not.** Every request and response is saved to `04_agent/transcript.jsonl` so a single run can
be audited, which is not the same as the method being reproducible. See the
[SC-03 contract](docs/scenarios/SC-03-jev-agent-design.md) and the
[succinate example](examples/jev-design/README.md).
The command writes raw tables, archives both expression inputs, renders R figures and both
HTML reports, and validates the bundle. `--analysis-only` omits R rendering;
`scripts/reproduce.py` replays the archived inputs after moving the run directory.

## Availability and implementation

Source code, documentation, test data, and reproducibility workflows are available
at <https://github.com/jyryu3161/CMM> under the [MIT License](LICENSE). CMM is implemented in
Python and supports Python 3.10–3.12 on Linux, macOS, and Windows. It provides both a Python
API and a Qt desktop interface; installation and a complete test run require no registration.
Tagged releases and their test data will remain available for at least two years after
publication, with issue reporting through the repository's
[GitHub Issues](https://github.com/jyryu3161/CMM/issues).

Before journal submission, the exact `v0.5.0` release must additionally be archived in
Zenodo or an equivalent long-term repository and its DOI added to this section,
`CITATION.cff`, and the manuscript's Availability and Implementation statement.

## Implemented methods

- Simulation: FBA, pFBA, FVA, editable conditions, and growth-media presets.
- Perturbations: reaction/gene/multiple knockouts, L1/L2 MOMA, ROOM with both published
  tolerance pairs, and batch screens. The reported distance is a distance — Segrè et al.
  Eq. (4)'s Euclidean value for L2 MOMA — kept separate from the raw solver objective.
- Omics: LAD, strict two-stage E-Flux2, multi-condition flux prediction, and log2 changes.
- Production: theoretical yield with carbon/CO₂ disclosure (media presets close CO₂ *uptake*
  by default; secretion stays free, and the fraction of product carbon supplied by any residual
  CO₂ uptake is reported and warned on), production envelopes, FSEOF, and FVA-based FVSEOF
  classifying on Park et al.'s nine types, with optional caller-supplied linear flux couplings.
- Response and sampling: flux-response scans reporting the exact LP shadow price and its phase
  boundaries, and seeded random flux sampling both uniform and constrained around a reference
  flux state.
- Strain design: distinct OptKnock and three-level RobustKnock modules through StrainDesign,
  with an explicit deterministic seed forwarded to the MILP backend, followed by independent
  maximum/guaranteed-product evaluation. The default seed is `0`; publication configs record it
  rather than accepting a hidden backend-generated seed. Knockout candidates are restricted to
  gene-associated internal reactions, so a returned design is buildable as a gene deletion:
  exchanges are excluded, and so are reactions whose only gene is COBRA's `s0001`
  spontaneous-reaction placeholder.
- MTA/rMTA: published MTA MIQP, published rMTA best/MOMA/worst scoring, and an explicitly
  labeled legacy continuous heuristic. Reached as `revert_targets` and `transformation_targets`
  in Python and as the *Revert Metabolism* and *Transform (A→B)* tabs in the application.
- Visualization: flux maps of any current flux state — FBA, pFBA, an E-Flux2/LAD prediction,
  or a MOMA/ROOM redistribution, each labelled with the method that produced it — on a curated
  Escher layout — Escher's *E. coli* core map is bundled
  and offered automatically to any model containing at least half its reactions, including
  genome-scale ones — or a dependency-free schematic of the highest-flux reactions for models
  no map fits. Any Escher JSON can be loaded from the GUI.
- JEV agent: a decision model drives CMM's own services in a round-based loop to search for
  a production design, with viability enforced by the solver rather than predicted by the
  agent. Reached as `run_jev_design` in Python, `cmm jev-design` on the command line, and the
  *JEV Agent* tab in the application. Needs `OPENROUTER_API_KEY`; nothing else in CMM does.
- Auditability: deterministic model fingerprints and solver/package/parameter provenance on
  numerical results.
- Production workflow: the concrete `production_target_workflow` composes MOMA-L2/ROOM single
  knockouts, OptKnock/RobustKnock, FSEOF/FVSEOF, loop diagnostics, flux response, paired
  sampling, recommendations, and the fixed schema-v2 artifact tree. FSEOF and FVSEOF contribute
  independent method-specific rankings; overlap is reported but is not required for a
  hypothesis to receive its own forward validation. Every unique MOMA/ROOM display-ranked
  D1–D5 knockout candidate receives matched paired sampling. A representable single-reaction
  signature also receives a pre-deletion wild-type reaction-to-product response scan;
  multi-reaction signatures remain explicitly unavailable rather than being reduced to an
  arbitrary axis. Every candidate in the independent FSEOF/FVSEOF top-10 union receives a
  flux-response scan. Loop-flagged or unresolved amplification targets are still scanned but
  remain ineligible for support/recommendation; non-runnable or failed analyses remain visible
  with status and reason rather than shortening those candidate sets. The workflow is available
  through `cmm production-targets --config CONFIG` and the Python API.
- Transformation workflow: `cmm transformation-targets --config CONFIG` composes expression
  integration, direction testing, candidate construction, MTA/rMTA ranking, MOMA comparison,
  optional epsilon sensitivity, archived-input replay, R reporting and validation. Its Python
  boundary is `cmm.workflows.transformation.run_transformation_target_discovery`.
- Publication reporting: the concrete `publication_reporting` service validates
  manifest-declared source tables and renders deterministic English HTML plus 300-DPI PNG and editable
  PDF/SVG figures through R. Generic scenario-template and scenario-file-format engines remain
  explicitly excluded from the shipped feature manifest.

## Reproducible installation

CMM requires Python 3.10–3.12. Exact Python dependencies are recorded in `uv.lock`; Python 3.12
is the preferred cross-platform interpreter. The source installer is the shortest editable setup:
it resolves compatible dependencies from the project metadata, does not rely on the operating
system's `python3`, and locates a newly bootstrapped `uv` directly even before the current shell's
`PATH` is refreshed. Use the separate frozen command below when the lock must be enforced.

```bash
git clone https://github.com/jyryu3161/CMM.git
cd CMM
./install.sh --dev
.venv/bin/cmm --version
.venv/bin/python -m cmm.app
```

On Windows PowerShell, run `.\install.ps1 -Dev` and launch with
`.\.venv\Scripts\python.exe -m cmm.app`. The installers accept an existing `uv >= 0.8.0` and
print the version and executable they use. Only when `uv` is absent do they bootstrap the tested
0.12.5 release. An older installed `uv` is rejected with pinned-upgrade instructions.

For an exactly frozen CI or manuscript environment, use `uv` 0.12.5 with `uv.lock`. A fresh
Astral bootstrap normally writes `uv` and a shell environment file under `$HOME/.local/bin`; load
that file (or prepend the directory) before invoking `uv` in the same shell:

```bash
curl -LsSf https://astral.sh/uv/0.12.5/install.sh | sh
[ -f "$HOME/.local/bin/env" ] && . "$HOME/.local/bin/env"
export PATH="$HOME/.local/bin:$PATH"
uv --version  # must report uv 0.12.5 for this frozen path
uv python install 3.12
uv sync --python 3.12 --frozen --all-extras
uv run cmm --version
uv run python -m cmm.app
```

For the equivalent manual PowerShell flow, bootstrap with
`irm https://astral.sh/uv/0.12.5/install.ps1 | iex`, invoke
`$env:USERPROFILE\.local\bin\uv.exe`, and confirm it reports 0.12.5 before running the
`python install` and `sync --frozen` subcommands. R is optional for the scientific services,
workflow analysis, CLI validation, and desktop application.

Publication rendering additionally requires R. In a source checkout, `renv.lock` pins R 4.3.2
and every renderer package. Restore it, then expose the project library to the current shell;
the export is needed because CMM deliberately launches `Rscript --vanilla`.

macOS/Linux/WSL:

```bash
Rscript --vanilla -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(prompt=FALSE)'
CMM_RENV_LIBRARY="$(Rscript --vanilla -e 'renv::load(); cat(.libPaths()[1])')"
export R_LIBS_USER="$CMM_RENV_LIBRARY"
uv run cmm report render RUN_DIR
```

Windows PowerShell:

```powershell
Rscript --vanilla -e "if (!requireNamespace('renv', quietly=TRUE)) install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(prompt=FALSE)"
$env:R_LIBS_USER = ((Rscript --vanilla -e "renv::load(); cat(.libPaths()[1])") -join "").Trim()
uv run cmm report render RUN_DIR
```

The Python wheel and sdist contain the checked-in R renderer. The sdist also contains
`renv.lock`; a wheel-only installation does not. Wheel-only users can install the renderer's
runtime packages directly:

```bash
Rscript --vanilla -e "install.packages(c('jsonlite','ggplot2','ggrepel','patchwork','svglite','ragg'), repos='https://cloud.r-project.org')"
```

That follows current CRAN rather than the publication lock. For an exactly reproducible
manuscript render, restore `renv.lock` from the matching Git checkout, sdist, or GitHub source
archive.
At runtime the renderer verifies that `Rscript` and its named packages can be loaded; those
packages enforce their declared compatible dependency minima. The runtime records the actual
versions in the figure manifest but does not reinterpret the repository lock. Exact versions
come from `renv::restore()` and are compared with the complete lock by CI on all three OSes.

The current lock records exact R/package versions and CRAN sources but has no per-package
`Hash` fields. Do not add hashes by hand: regenerate them only from a fully materialized R
4.3.2 library when a canonical `renv::snapshot()` diff preserves the complete package set and
every recorded version.

CI restores that lock and executes the R package smoke test plus publication-report tests on
Ubuntu, macOS, and Windows. `ggplot2` itself is pure R; native compilation risk comes from
graphics dependencies such as `ragg`, `systemfonts`, `textshaping`, and `svglite`. The CI job
uses available CRAN/Posit binaries, installs Linux and macOS build libraries, and provisions
Rtools43 on Windows. Those toolchains prepare and support source compilation when an archived
R 4.3 binary is unavailable. The normal matrix does not force a source-only restore; it tests
the exact restored versions and renderer result regardless of whether each package arrived as
a binary or was compiled from source.

For a conventional editable installation, invoke a supported interpreter explicitly (on Windows,
use `py -3.12` or the full path to Python 3.10–3.12):

```bash
python3.12 -m pip install -e ".[desktop,design,solver-gurobi]"
python3.12 -m cmm.app
```

The cross-platform source installers create new environments with a uv-managed Python 3.12, so
an older system Python such as macOS 3.9.6 is not selected. They use an installed `uv >= 0.8.0`
or bootstrap 0.12.5 when none is available. They install the desktop, strain-design, and Gurobi
extras by default:

```bash
./install.sh                 # macOS / Linux / WSL
# .\install.ps1             # Windows PowerShell
./install.sh --dev           # also install test, coverage, lint, and type-check tools
./install.sh --no-gurobi     # GLPK LP/MILP only
CMM_PYTHON=3.11 ./install.sh # ask uv to install/use another supported version
./install.sh --python /path/to/python3.12  # command-line override wins over CMM_PYTHON
```

PowerShell accepts the corresponding `-Dev`, `-NoGurobi`, `-CoreOnly`, `-Python`, and
`-VenvDir` options; its environment override is `$env:CMM_PYTHON = "3.11"`. A numeric override
asks uv to install that managed version. A command or path must already identify an executable.
Both installers reject Python outside 3.10–3.12 and never silently replace an existing virtual
environment that conflicts with an explicit override. Reuse also requires `pyvenv.cfg` and an
interpreter that reports `sys.prefix != sys.base_prefix`; an ordinary directory containing a
Python executable is not trusted as a venv.

If a checkout already has a legacy Python 3.9 `.venv`, preserve it while recovering. Either move
it aside and let the installer create a fresh `.venv`, or install into a new path:

```bash
mv .venv .venv-python39-backup
./install.sh
# Alternatively, leave .venv untouched:
./install.sh --venv .venv312
```

PowerShell equivalents are `Rename-Item .venv .venv-python39-backup; .\install.ps1` or
`.\install.ps1 -VenvDir .venv312`. The installers do not delete or overwrite the old environment.

Tagged wheels and source archives are published on the
[GitHub Releases page](https://github.com/jyryu3161/CMM/releases). The current source version
is 0.5.0; 0.4.0 was a **breaking** release — see `CHANGELOG.md`.

## Solver requirements

Solver capability is checked before a solve; a method never silently changes formulation.

| Class | Methods |
|---|---|
| LP | FBA, pFBA, FVA, LAD, yield, envelope, FSEOF, FVSEOF, flux response, flux sampling |
| MILP | ROOM, OptKnock, RobustKnock |
| QP | L2 MOMA, E-Flux2, `rmta_continuous` |
| MIQP | published MTA and rMTA |

GLPK supports LP/MILP. Gurobi and CPLEX support the full table, subject to their licenses and
model-size limits. The restricted Gurobi license bundled with `pip install gurobipy` is limited
to 2000 variables and 2000 constraints, dropping to **200 variables for any model containing
quadratic terms** — the cap counts variables, not quadratic terms. Because a COBRA LP uses
about two variables per reaction, that restricted license covers CMM's QP/MIQP *validation*
models (all ≤190 variables) but only permits QP/MIQP on networks of roughly **100 reactions or
fewer** — and fewer still for L2 MOMA, which adds one variable per reaction: **L2 MOMA on
`e_coli_core` is 286 variables and already fails** under the restricted license, though it is
neither genome-scale nor mixed-integer. Genome-scale work of any kind, including plain FBA on
iJO1366 (5166 variables), requires a full academic or commercial license.

The `rmta_continuous` QP path is exercised on a small deterministic network and explicitly
labels its CSV as a continuous heuristic. It is not published rMTA. See
[Known limits](docs/VALIDATION.md#known-limits).

### Expect MTA/rMTA on a genome-scale model to take hours

A slow MIQP screen is the normal cost of these methods, not a hang or a misconfiguration. One
measured reference point: published rMTA over Recon1 (3741 reactions, 1905 genes) at the gene
level, 413 candidates after the blocked-reaction reduction, with Gurobi 13 on an 8-core Apple
M4, took **about 95 s per candidate and 11.5 h end to end**. rMTA solves three MIQPs per
candidate — best case, MOMA and worst case — against MTA's one, and candidates are evaluated
sequentially, so wall time scales with the candidate count and extra cores do not shorten it
much. Each additional value in `validation.epsilon_sweep` re-ranks every candidate and costs
another full pass.

Run these screens in the background rather than interrupting and retrying with smaller
parameters, and size a run from its candidate count before launching it: cost is roughly the
candidate count times the per-candidate solve time, so a smaller model with far fewer
candidates is correspondingly cheaper.

## Python quick start

```python
from cobra.io import load_model
from cmm.core import apply_medium, fba, pfba
from cmm.features import fseof, theoretical_yield

model = load_model("textbook")
apply_medium(model, "glucose_anaerobic")     # the condition, set once

growth = fba(model)
minimal = pfba(model)
yield_result = theoretical_yield(model, "EX_succ_e")
scan = fseof(model, "EX_succ_e", n_steps=10)

print(growth.objective_value, minimal.status)
print(yield_result.molar_yield, yield_result.metadata["model_sha256"])
print(scan.amplification_targets())
```

Calling the response and sampling building blocks directly:

```python
from cmm.features import flux_response, random_flux_sampling

response = flux_response(model, "PGI", "EX_succ_e", biomass_fraction=0.3)
print(response.optimum(), response.limit.found, response.feasible_range())

ensemble = random_flux_sampling(model, n=1000, seed=0)
print(ensemble.statistics().loc["EX_succ_e"])
```

The sampling call above characterizes the model as passed; it is not paired knockout
validation by itself. The canonical production workflow applies every unique MOMA/ROOM D1–D5
deletion and exports matched wild-type/knockout ensembles for each candidate, while applying
flux response to those knockout candidates and the complete FSEOF/FVSEOF top-10 union.

Expression integration:

```python
import pandas as pd
from cmm.omics import flux_log_change, predict_condition_fluxes

expression = pd.read_csv("expression.csv").set_index("gene")
predicted = predict_condition_fluxes(model, expression, method="eflux2")
change = flux_log_change(
    predicted.fluxes("condition_A"),
    predicted.fluxes("condition_B"),
)
```

## Validation

The suite includes direct COBRApy cross-checks, a non-optional iJO1366 genome-scale test,
the official COBRA Toolbox MTA test topology, E-Flux2's independent two-stage QP, target
sensitivity checks, offscreen GUI workflows, static analysis, coverage, and distribution
build verification.

For focused offscreen and workflow checks with the same locked dependencies:

```bash
QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest -q -ra \
  tests/test_app_smoke.py tests/test_scenarios.py
uv run --frozen --all-extras pytest -q -ra \
  tests/test_agent_contract.py tests/test_production_workflow.py \
  tests/test_transformation_workflow.py tests/test_publication_reporting.py \
  tests/test_transformation_reporting.py
```

These tests include numerical checks, orchestration fixtures, and report validation. To keep
the actual GUI captures, run the three commands in the
[scenario-figure manifest](docs/scenario-figures.md). Solver-specific tests require the listed
capabilities and a sufficient license; report tests require the R environment described above.

```bash
uv sync --frozen --all-extras
QT_QPA_PLATFORM=offscreen uv run pytest -q -ra --strict-markers \
  --durations=10 --cov=cmm --cov-branch --cov-report=term-missing --cov-fail-under=80
uv run ruff check src tests
uv run ruff format --check src tests
uv run mypy src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
uvx --from cffconvert==2.0.0 cffconvert --validate
uv build && uvx twine check dist/*
```

CI keeps every Linux, Windows, and macOS test shard required, but evaluates the 80% branch-
coverage policy once over the combined Linux 3.12 evidence from the solver-neutral suite, the
restricted-license QP/MIQP suite, and the locked-R publication suite. Capability-specific
deselection or an unavailable R installation therefore cannot turn a successful platform test
shard into a misleading coverage failure, while the combined public and publication surface
must still meet the same threshold.

The exact evidence, method contracts, references, provenance schema, and limitations are in
[Scientific validation and reproducibility](docs/VALIDATION.md). Passing these checks
supports implementation correctness; it does not constitute wet-lab validation of a new
biological prediction.

## Documentation

- [Build or customize a reproducible workflow](docs/building-custom-workflows.md)
- [Contributor tutorial: add a canonical workflow](docs/tutorials/adding-a-canonical-workflow.md)
- [Desktop and Python tutorial](docs/TUTORIAL.md)
- [Scientific validation and reproducibility](docs/VALIDATION.md)
- [MTA/rMTA design and equations](docs/design-revert-metabolism.md)
- [Architecture and solver contracts](docs/architecture.md)
- [Contributor scenario-figure manifest](docs/scenario-figures.md) — GUI regression captures and regeneration commands, not publication evidence
- [AI-assisted use and disclosure](docs/AI-USAGE.md)
- [Release changes](CHANGELOG.md)

For driving CMM from a repository-aware AI coding tool:

- [Agent operating instructions](AGENTS.md) — the scenario router, solver gate, and run contract
- [Metabolic-engineering scenarios](docs/scenarios/README.md) — step-by-step pipelines
- [Function reference for agents](docs/agent-reference.md) — signatures and result objects

Two tested Codex/OpenAI-compatible repository skills route requests to the corresponding CLI:

| Goal | Repository skill | Public command |
|---|---|---|
| Production / knockout / amplification targets | [cmm-production-engineering](.agents/skills/cmm-production-engineering/SKILL.md) | `cmm production-targets` |
| Move a source state toward a target with MTA/rMTA | [cmm-transformation-engineering](.agents/skills/cmm-transformation-engineering/SKILL.md) | `cmm transformation-targets` |

Skill auto-discovery is host-dependent and is not
claimed for other agent hosts; they can still follow `AGENTS.md`, the scenarios, and the public
API directly. Both skills ship in the GitHub/source checkout and sdist, but are not Python
runtime dependencies or installed-wheel features. Installing the wheel alone provides the
workflow APIs and CLI, not the repository-level skills.

The documentation and agent contract are reproducibility/source materials rather than numerical
runtime dependencies. They ship in the sdist; the wheel contains the installed Python runtime,
R renderer, and licensed map resources. `.github/workflows/` is GitHub-only maintenance
automation: it is not imported by CMM, but a green run for the archived commit is part of the
software-validation evidence.

## Citation and license

CMM is open source under the [MIT License](LICENSE); external inspiration and redistributed
material are recorded in [Third-party notices](THIRD_PARTY_NOTICES.md). Citation metadata are
machine-readable in [CITATION.cff](CITATION.cff). The file intentionally still contains an
organization-only placeholder: replace it with the final authors and ORCIDs, archive the exact
release, and add the DOI before journal submission. Until then, cite the repository URL, version,
and commit.
A manuscript must also cite the original paper for every method it uses, as listed in
[Scientific validation and reproducibility](docs/VALIDATION.md).

Three points that are easy to get wrong, all detailed there:

- **OptKnock and RobustKnock need three citations, not two** — Burgard et al. (2003) and
  Tepper & Shlomi (2010) for the formulations, *and* Schneider et al. (2022) for the
  `straindesign` package that actually solves them, which carries no citation of its own.
- **`transformation_targets` is not a CMM invention**; both of its paths map to published
  Yizhak et al. (2013) methods and should be cited to that paper.
- **The bundled Escher map is Escher's, not CMM's.** `src/cmm/resources/` redistributes the
  *E. coli* core map under Escher's MIT license; work using that layout should cite King et al.
  (2015). Provenance, digest, and license are in `src/cmm/resources/ATTRIBUTION.md`.
- **`flux_log_change` has no published source.** It is a CMM utility and must not be cited to
  any paper. CMM's FSEOF selection rule and FVSEOF's `robust_targets()` flag are likewise
  CMM's own and must not be attributed to Choi et al. or Park et al.

## Release process

Every push to `main` and every pull request installs the frozen lockfile and runs the
cross-platform quality gates. A tag must exactly match `pyproject.toml`; the release workflow reruns all checks,
builds the wheel and sdist, validates them, installs the wheel in a clean environment, and
then attaches the artifacts to a GitHub Release.

```bash
git tag v0.5.0
git push origin v0.5.0
```
