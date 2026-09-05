# Building or customizing a CMM workflow

CMM ships **two canonical workflows**: SC-01 production-target discovery and SC-02
transformation-target discovery. Each binds a scientific question to an exact model, explicit
conditions, solver requirements, typed numerical results, provenance, a report renderer, and an
artifact directory. An analysis recipe composed over public services is a legitimate way to
work, but it is neither an installed workflow command nor a validated run schema.

This guide separates two legitimate extension paths:

- **Track A — downstream use:** configure SC-01 or SC-02, or compose a private reproducible study from
  public CMM services without changing CMM itself.
- **Track B — upstream contribution:** add a new canonical workflow, schema, reporter,
  validator, CLI boundary, documentation, and tests to CMM.

The desktop application is useful for exploration, but a workflow intended for a paper should
run through Python or the thin CLI so that its inputs and outputs can be replayed.

## 1. Choose the smallest correct boundary

There are three ways to use CMM. Choose by scientific scope, not by preferred interface.

| Need | Boundary | Reproducible output |
|---|---|---|
| Complete production-target study | `ProductionWorkflowConfig` and `run_production_target_discovery` | Canonical schema-v2 run, R report, validator |
| Which knockout moves one metabolic state toward another | `TransformationWorkflowConfig` and `run_transformation_target_discovery` | Canonical schema-v2 run, R report, validator |
| Same study with different thresholds, candidate counts, search seed, or sampling settings | Change the canonical config | Same canonical schema and validation contract |
| A different scientific question or a single analysis | Track A: compose documented functions in `cmm.core`, `cmm.features`, or `cmm.omics` | A downstream study with its own declared outputs |
| Another reusable, installed CMM workflow | Track B: contribute a workflow-specific API, schema, reporter, validator, CLI, and tests | A separately versioned canonical run contract |

Do not create a second SC-01 implementation just to change parameters. The existing production
workflow already exposes stage switches, method limits, the number of method-specific FSEOF and
FVSEOF targets, validation settings, a strain-design seed, and a sampling seed. Use direct
service calls when the question itself is different—for example, a standalone FVA study or a
condition comparison.

Start with the [scenario router](scenarios/README.md), then use the relevant signatures in the
[function reference](agent-reference.md). The [architecture](architecture.md) explains which
layer owns numerical work, workflow composition, rendering, and validation.

## 2. Write the scientific contract first

Before writing code or JSON, state the following in a short protocol:

1. **Question and decision rule.** What result changes the conclusion, and which thresholds are
   declared before looking at the ranking?
2. **Exact inputs.** Identify the SBML path, product exchange where applicable, expression
   files, and any id mapping. A model name is not a model version.
3. **Declared conditions.** For every state the question compares, record the medium, substrate
   exchange and uptake, oxygen exchange and bounds, objective, and every other changed bound.
   A single-condition workflow still declares exactly one complete condition.
4. **Method roles.** Separate inverse target discovery from forward phenotype prediction. Do
   not merge unlike scores into one ranking.
5. **Capabilities.** State the required LP, MILP, QP, or MIQP support and any optional package
   requirement. Define whether a missing capability stops or explicitly narrows the run.
6. **Artifacts.** Name each table, metadata sidecar, figure, and narrative claim that the run
   must produce.
7. **Completion gate.** Define the validator and tests that must pass. A generated HTML file is
   not itself proof that the workflow completed correctly.

The shared [preflight](scenarios/_preflight.md) and
[reporting contract](scenarios/_reporting.md) provide the detailed checks.

## 3. Track A1 — customize the canonical production workflow

The recommended entry point for a production-engineering study is a UTF-8 JSON config. The
following generic example is intentionally not tied to one organism. Replace every reaction id
with an id present in the exact SBML model.

```json
{
  "model_path": "../model/organism.xml",
  "product": "EX_product_e",
  "output_dir": "../results/product_aerobic",
  "solver": "gurobi",
  "substrate": "EX_substrate_e",
  "biomass": "BIOMASS",
  "medium": {
    "mode": "model_as_loaded"
  },
  "condition": {
    "name": "defined_aerobic",
    "objective": {
      "coefficients": {
        "BIOMASS": 1.0
      },
      "direction": "max"
    },
    "bounds": [
      {
        "reaction_id": "EX_substrate_e",
        "lower_bound": -10.0,
        "upper_bound": 1000.0
      },
      {
        "reaction_id": "EX_o2_e",
        "lower_bound": -20.0,
        "upper_bound": 1000.0
      }
    ],
    "notes": "Model-as-loaded medium explicitly accepted after inspection; substrate and oxygen bounds overridden here"
  },
  "top_single_knockouts_per_method": 5,
  "top_amplification_targets_per_method": 10,
  "max_knockouts": 3,
  "strain_design_seed": 0,
  "fseof_steps": 10,
  "fvseof_steps": 10,
  "amplification_loop_diagnostic_top_n": 20,
  "validation": {
    "enabled": true,
    "max_flux_response_targets": 30,
    "flux_response_steps": 20,
    "flux_response_biomass_fraction": 0.3,
    "sampling_growth_fraction": 0.1,
    "sampling": {
      "enabled": true,
      "n": 1000,
      "method": "achr",
      "thinning": 100,
      "processes": 1,
      "seed": 0,
      "store_raw_samples": true
    }
  }
}
```

The generic example uses `model_as_loaded` because a biologically complete medium cannot be
invented from only substrate and oxygen. Use this mode only after inspecting and explicitly
accepting every uptake bound in the archived SBML; the condition then overrides the named
substrate and oxygen bounds and explicitly replaces the loaded objective with maximum biomass.
The selected full-capability solver must be installed and available to COBRApy, and the full
strain-design stage also requires importable `straindesign`. For a publication run, an explicit
medium is preferable and its positive uptake mapping must list **every required nutrient**, not
only carbon and oxygen. Relative `model_path` and `output_dir` values resolve from the config
file's directory.

Keep `strain_design_seed` explicit even though its default is `0`. It is forwarded unchanged to
both OptKnock and RobustKnock and must be an integer in `0..2_000_000_000`; booleans and floats
are invalid. It is separate from `validation.sampling.seed`. Omitting explicit seed forwarding
would allow `straindesign` to invent a new random seed per call, changing MILP search paths,
runtime, or finite solution pools in otherwise identical runs.

Run analysis, rendering, and validation together:

```bash
uv run cmm production-targets --config workflows/product.json
uv run cmm report validate results/product_aerobic --json
```

Or separate numerical analysis from the R renderer:

```bash
uv run cmm production-targets --config workflows/product.json --analysis-only
uv run cmm report render results/product_aerobic
uv run cmm report validate results/product_aerobic --json
```

The equivalent Python boundary is useful when the workflow is called by another application:

```python
from cmm.reporting import render_production_report, validate_production_run
from cmm.workflows.production import (
    ProductionWorkflowConfig,
    ProductionWorkflowResult,
    run_production_target_discovery,
)

config = ProductionWorkflowConfig.from_json("workflows/product.json")
result: ProductionWorkflowResult = run_production_target_discovery(config)
if result.run_directory is None:
    raise RuntimeError("output_dir is required for a reportable run")

render_production_report(result.run_directory, renderer="nature-r")
validation = validate_production_run(result.run_directory)
validation.raise_for_errors()
```

FSEOF and FVSEOF remain independent methods. The workflow exports and plots the top configured
number from **each** ranking; overlap is descriptive, not an admission rule. Every
method-specific trajectory retains its diagnostic status, and every unique reaction in the two
report-visible lists proceeds to its own flux-response check. Loop-flagged or unresolved
candidates are still scanned but cannot be promoted to supported recommendations.
Consequently, the diagnostic cap must be at least twice the per-method amplification count, and
the response capacity must also accommodate both per-method D1–D5 single-knockout candidate
lists. With 10 amplification targets and five single knockouts per method, those minima are 20
and 30. These limits are capacity guards before deduplication, not top-*N* selectors. Config
validation rejects smaller enabled capacities, and execution never silently drops the tail of
a method.

The canonical knockout validation universe is the unique blocked-reaction signatures
represented by MOMA D1–D5 and ROOM D1–D5. Every representable single-reaction candidate receives
a pre-deletion target-reaction-to-product response scan, and every candidate receives matched
wild-type/knockout sampling. Multi-reaction signatures remain explicitly unavailable for the
one-axis response scan rather than being reduced to an arbitrary reaction. Equivalent
signatures are simulated once, with all represented gene ids retained as candidate-id
provenance. Every candidate also remains in the response and sampling indexes if its analysis
is infeasible, unavailable, skipped, or failed.

### Customize the transformation workflow (SC-02)

For source → target studies, use `TransformationWorkflowConfig` from
`cmm.workflows.transformation` and the `cmm transformation-targets` command. A completed run's
`00_config.json` records the full configuration; `scripts/transformation_config.json` is its
replay configuration. When copying either file into a new study directory, update its relative
paths because they resolve from that file's directory.

Start with the [runnable MTA/rMTA example](../examples/transformation-targets/README.md), which
also lists the inputs needed to adapt it to Recon. The
[succinate example](../examples/production-targets/README.md) provides the matching SC-01 entry
point; neither requires copying numerical workflow code.

1. Set `model_path`, `source_expression_path`, `target_expression_path` and a new `output_dir`.
   Confirm which state to move away from and toward; this cannot be inferred from the model.
2. Set the medium and condition, including substrate and oxygen bounds, for that exact model.
   Inspect exchange ids and verify growth; a Recon1 condition cannot be transferred unchanged
   to another Recon release with different ids or reactions.
3. Provide finite, non-negative linear expression with matching model gene ids, one row per
   gene and one column per replicate. Keep the original measurements: SC-02 performs the
   `log2(value + 1)` comparison itself. For replicated data use `direction.significance="ttest"`;
   single measurements require an explicit fold-change configuration and disclosure.
4. Choose `method` (`mta` or `rmta`), `perturbation`, the reference method, `epsilon`, the
   changed-reaction cutoff, and any `validation.epsilon_sweep`. Epsilon is a flux magnitude;
   do not inherit it from a different model without checking the reference flux scale.
5. Run the same public boundary and completion gate:

   ```bash
   uv run --frozen --all-extras cmm transformation-targets --config studies/my-transformation.json
   uv run --frozen --all-extras cmm report validate results/my-transformation --json
   ```

The workflow archives the model and both expression files, constructs candidates, ranks them,
compares the MOMA baseline, and writes R figures plus both HTML reports. Preserve explicit
failure/status rows and tie warnings. A new dataset or condition is a new **configuration** of
SC-02; it does not need copied solver code or a third installed workflow. A different scientific
question that needs its own analyses and artifact contract belongs in Track B below.

See the [SC-02 contract](scenarios/SC-02-transformation-target-discovery.md) and
[public signatures](agent-reference.md#complete-transformation-workflow--cmmworkflowstransformation)
for the complete fields and their meanings.

## 4. Understand what the canonical run owns

The production workflow writes the exact source and conditioned models, resolved config,
provenance, raw CSVs, per-analysis metadata, and an authoritative manifest before reporting.
Its high-level flow is:

```text
config + exact SBML
        │
        ▼
preflight ─► yield/reference ─► MOMA + ROOM single deletions
        │                           │
        ├────────► OptKnock + RobustKnock
        │
        └────────► FSEOF + FVSEOF (independent rankings)
                              │
                              ▼
             all-candidate response / loop / KO sampling checks
                              │
                              ▼
                  manifest ─► R figures + HTML ─► validator
```

`00_manifest.json` is the authority for artifacts; the renderer does not discover files by
filename patterns. `render_production_report(..., renderer="nature-r")` reads only the declared
source tables and writes 300-DPI PNG plus editable PDF/SVG figures and linked/standalone HTML.
`validate_production_run` recomputes declared sizes and SHA-256 values, checks required methods
and columns, and verifies figure/report coverage without changing results.

The generated `scripts/production_config.json`, `reproduce.py`, `render.py`, and `validate.py`
are the run-local replay boundary. Preserve them with the exact model, `uv.lock`, `renv.lock`,
CMM version or commit, and the raw tables used by the manuscript.

## 5. Track A2 — compose a downstream study from public services

When the question is outside both canonical workflows, compose the documented public services and keep their typed
results intact. This small example asks whether one FSEOF hypothesis has a supportive
target-to-product response. It is a narrow in-memory study, **not** an installed CMM workflow,
CLI command, artifact schema, or report contract.

```python
from dataclasses import asdict, dataclass
from typing import Mapping

from cobra.io import read_sbml_model

from cmm.core import (
    Condition,
    FluxSolution,
    Medium,
    ObjectiveSpec,
    apply_medium,
    fba,
    require,
    run_provenance,
)
from cmm.features import (
    FluxResponseResult,
    FseofResult,
    ProductionYield,
    flux_response,
    fseof,
    theoretical_yield,
)


@dataclass(frozen=True)
class TargetResponseStudy:
    baseline: FluxSolution
    ceiling: ProductionYield
    ranking: FseofResult
    target: str
    response: FluxResponseResult
    metadata: Mapping[str, object]


def run_target_response_study(
    model_path: str,
    *,
    product: str,
    target: str,
    biomass: str,
    medium: Medium,
    substrate_exchange: str,
    oxygen_exchange: str,
) -> TargetResponseStudy:
    model = read_sbml_model(model_path)
    applied = apply_medium(model, medium)
    missing = {substrate_exchange, oxygen_exchange} - set(applied)
    if missing:
        raise ValueError(
            "the loaded model could not apply the declared substrate/oxygen exchanges: "
            f"{sorted(missing)}"
        )
    condition = Condition(
        name=medium.name,
        objective=ObjectiveSpec(coefficients={biomass: 1.0}, direction="max"),
        notes="Complete caller-supplied medium and explicit biomass objective",
    )
    require("LP", model.solver.interface, feature="target-response study")

    baseline = fba(model, condition=condition)
    if baseline.status != "optimal" or baseline.objective_value <= 1e-6:
        raise RuntimeError("wild type does not grow in the declared condition")

    ceiling = theoretical_yield(
        model,
        product,
        substrate=substrate_exchange,
        condition=condition,
    )
    if ceiling.molar_yield <= 1e-6:
        raise RuntimeError("product is unreachable in the declared condition")

    ranking = fseof(model, product, biomass=biomass, condition=condition, n_steps=10)
    if target not in ranking.amplification_targets():
        raise RuntimeError("declared target is not an actionable FSEOF hypothesis")

    response = flux_response(
        model,
        target,
        product,
        biomass=biomass,
        condition=condition,
        biomass_fraction=0.3,
        n_steps=20,
    )
    return TargetResponseStudy(
        baseline=baseline,
        ceiling=ceiling,
        ranking=ranking,
        target=target,
        response=response,
        metadata=run_provenance(
            model,
            method="target_response_study",
            product=product,
            target=target,
            biomass=biomass,
            condition=asdict(condition),
            applied_medium=applied.to_provenance(),
            substrate=substrate_exchange,
            substrate_uptake=applied[substrate_exchange],
            oxygen_exchange=oxygen_exchange,
            oxygen_uptake=applied[oxygen_exchange],
        ),
    )
```

The caller can now export `baseline.to_frame()`, `ceiling.to_frame()`, `ranking.to_frame()`, and
`response.to_frame()` without reconstructing scientific results from selected scalars. The
`substrate_exchange` and `oxygen_exchange` arguments must be resolved reaction ids from the
loaded model; checking them against `MediumApplication` prevents a dropped medium component
from being reported as applied. The
study still lacks a config file, manifest, replay boundary, renderer, validator, and CLI, so do
not present it as a canonical workflow. Downstream code may add its own declared artifact
contract, but it must not claim compatibility with the production schema merely because the
directory contains similarly named files.

## 6. Track B — contribute a canonical workflow to CMM

Adding a reusable workflow to CMM is a larger product boundary than composing service calls.
Follow the complete contributor tutorial:

- [Adding a canonical workflow to CMM](tutorials/adding-a-canonical-workflow.md)

That tutorial uses SC-01 and SC-02 as its two shipped references and generic `MyWorkflow...`
names for non-installed skeletons. It covers the scientific contract, typed public services, config and
result types, orchestration, schema and manifest, dedicated validation and R reporting, CLI and
package exports, documentation, tests, and release criteria.

Do not import or copy private helpers from `cmm.workflows.production`, and treat
`cmm.workflows._bundle` the same way — it is shared between the two shipped workflows, but the
leading underscore means it carries no stability promise for anything outside the package. The
current artifact contracts, `validate_run` implementation, report narrative, figure order, and
R script are production-specific even where a name appears generic. A new workflow needs its own schema id
and version, semantic roles, scientific invariants, renderer inputs, and validation boundary.
Only genuinely reusable infrastructure should be extracted, named neutrally, and tested before
a second workflow depends on it.

## 7. Test a study or workflow as a research instrument

A new workflow is not complete until tests cover both successful results and scientific
failure modes.

- **Config tests:** JSON round-trip, path resolution, invalid ranges, and explicit stage
  disablement.
- **Preflight tests:** zero growth, missing exchange, zero yield, and missing solver capability
  must stop or produce the documented unavailable state.
- **Composition tests:** assert service call order, exact condition propagation, complete
  declared coverage without runtime slicing, capability gates, and deterministic parameters.
- **Scientific invariants:** define assertions for the workflow's own question. SC-01's FSEOF,
  FVSEOF, MOMA, ROOM, response, sampling, and RobustKnock invariants are examples, not universal
  requirements for another workflow.
- **Artifact tests:** required roles, metadata sidecars, hashes, path containment, replay scripts,
  and no hand-transcribed report values.
- **Rendering tests:** source-data traceability, non-empty PNG/PDF/SVG siblings, HTML figure
  coverage, and visual inspection for clipped labels.
- **Regression tests:** run the repository's full offscreen quality gate, not only the new test
  module.

```bash
QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest -q \
  --cov=cmm --cov-branch --cov-fail-under=80
uv run --frozen --all-extras ruff check src tests
uv run --frozen --all-extras ruff format --check src tests
uv run --frozen --all-extras mypy \
  src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
```

## 8. Prepare the publication handoff

Archive the following together:

- exact release tag or commit, `uv.lock`, `renv.lock`, and solver/version information;
- exact input and conditioned models, including their redistribution terms;
- resolved config, provenance, manifest, raw tables, validation result, and replay scripts;
- 300-DPI raster and editable vector figures with their source-data tables;
- CMM's citation plus the original paper for every method used;
- declared medium, aeration, substrate uptake, thresholds, search limits, sample count,
  thinning, strain-design seed, and sampling seed;
- limitations covering alternative optima, solver tolerances, GPR resolution, sampler
  convergence, model scope, and experimental validation.

Passing CMM's tests supports implementation and artifact integrity. It does not establish the
biological validity of a predicted intervention. Present targets as *in silico* hypotheses
until they have been tested experimentally.

Further reading:

- [Contributor tutorial: add a canonical workflow](tutorials/adding-a-canonical-workflow.md)
- [SC-01 production target discovery](scenarios/SC-01-production-target-discovery.md)
- [Reporting and artifact contract](scenarios/_reporting.md)
- [Scientific validation and method references](VALIDATION.md)
- [Public function and result-object reference](agent-reference.md)
- [Repository architecture](architecture.md)
