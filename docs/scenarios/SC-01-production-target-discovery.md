---
id: SC-01
title: Production target discovery
goal: Find and verify genetic interventions that increase production of a target metabolite
when_to_use:
  - "increase production of X"
  - "which genes should I over-express or knock out to make more X"
  - "design a strain where production is guaranteed, not merely possible"
  - "find production-enhancing targets"
  - "design a growth-coupled strain"
role: spine
requires:
  model: "path to the exact COBRA model"
  product: "exchange reaction id, for example EX_succ_e"
  condition: "explicit medium, substrate uptake, oxygen/aeration bounds, and other changed bounds"
optional_inputs_from:
  - "a resolved biological condition from an omics comparison"
solver:
  canonical_workflow: "QP for MOMA-L2; MILP for ROOM and strain design"
  strain_design: "importable straindesign; any additional backend requirement is surfaced"
renderer: "Rscript plus loadable renderer packages; renv.lock and CI supply exact versions"
steps: [preflight, yield, reference, single_knockout, strain_design, amplification, validation, report]
runtime: "minutes on a core model; hours on a genome-scale model"
---

# SC-01 — Production target discovery

## Objective

Produce a reproducible *in silico* strain-engineering report that keeps four claims distinct:

1. MOMA and ROOM predict the phenotype of individual deletions.
2. OptKnock and RobustKnock search multi-reaction designs that couple product to growth.
3. FSEOF and FVSEOF nominate flux-amplification targets.
4. Flux-response indexes cover every canonical knockout and amplification candidate; numeric
   scans cover every amplification and every single-reaction knockout signature, while
   knockout-conditioned sampling covers every canonical single-knockout candidate. Each result
   may support, contradict, or leave the nomination inconclusive.

Predictions prioritize experiments; they do not establish a wet-lab phenotype.

## Canonical execution boundary

Use `.agents/skills/cmm-production-engineering/` for a complete production request. The
checked-in config and thin CLI are the reproducible boundary:

```bash
cmm production-targets --config CONFIG
cmm production-targets --config CONFIG --analysis-only
cmm report render RUN_DIR
cmm report validate RUN_DIR --json
```

The equivalent Python API is:

```python
from cmm.workflows.production import (
    ProductionWorkflowConfig,
    ProductionWorkflowResult,
    run_production_target_discovery,
)
from cmm.reporting import render_production_report, validate_production_run

config = ProductionWorkflowConfig(...)  # use the public fields in agent-reference.md
result: ProductionWorkflowResult = run_production_target_discovery(config)
assert result.run_directory is not None  # set when config.output_dir is supplied
render_production_report(result.run_directory, renderer="nature-r")
validation = validate_production_run(result.run_directory)
```

Do not replace this with a bespoke script that re-ranks CSVs or invents private thresholds.
Direct `cmm.features` calls remain appropriate for a request explicitly limited to one method.

## Resolve one biological condition before running

The run definition contains:

- exact model path; preflight-computed model id and fingerprint;
- product exchange reaction id;
- medium and every applied exchange bound;
- substrate exchange and uptake bound;
- oxygen exchange and aeration bound.

If the prompt or config already makes these values explicit and consistent, continue without
asking for redundant confirmation. Before asking, inspect the model and local capabilities and
resolve facts that are unique. If a remaining user decision is missing or ambiguous, follow
[`_preflight.md` step P0](_preflight.md#step-p0--resolve-consequential-ambiguity): ask only the
current decision, include a clearly labeled recommended option with its evidence and effect on
the scientific claim, and resolve dependent decisions in order. Do not infer conditions from a
filename, a prior run, or unconfirmed model bounds. Aerobic and anaerobic analyses are different
runs, not panels silently mixed into one run.

If clarification was required, show the exact model, product exchange, condition, method scope or
substitutions, and relevant experimental constraints in a final run-definition summary. Start the
canonical workflow only after the user confirms it. An initially complete request needs no extra
confirmation. Save the exact input model and resolved config.

## Pipeline and artifacts

| Step | Question | Required methods | Primary artifacts |
|---|---|---|---|
| 0 | Is the request and environment valid? | model/condition/solver/R preflight | `01_preflight/` |
| 1 | What is the production ceiling and growth trade-off? | theoretical yield, envelope | `02_yield/` |
| 2 | What is the wild-type reference state? | pFBA reference | `03_reference/` |
| 3 | Which single deletions improve product while retaining growth? | GPR mapping, MOMA-L2 and ROOM | `04_single_knockout/` |
| 4 | Which knockout sets couple product to growth? | OptKnock and RobustKnock | `05_strain_design/` |
| 5 | Which fluxes rise as product is enforced? | FSEOF and FVSEOF | `06_amplification/` |
| 6 | Do forward analyses support the candidates? | loopless FVA diagnostic, flux response and paired sampling | `07_validation/` |
| 7 | Can a reader audit every claim? | R rendering and run validation | `figures/`, reports, manifest |

The numbered directories are fixed. A skipped or unavailable method is recorded explicitly;
directories and tables are not silently renamed to hide it.

## Step 0 — Preflight

Run [`_preflight.md`](_preflight.md) under the confirmed condition before any search.

The complete workflow stops when:

- the model cannot grow under that condition;
- the product is not an exchange reaction or has zero theoretical yield;
- MOMA-L2 lacks QP support or ROOM lacks MILP support;
- OptKnock/RobustKnock were requested but MILP or importable `straindesign` is unavailable;
- the selected strain-design backend reports an unmet additional requirement;
- `renderer="nature-r"` was requested but `Rscript` or a required package is unavailable.

There is no silent MOMA-L1 substitution in the canonical workflow. A narrower, explicitly
approved run may use a different method, but its report must name the substitution and the
scientific claim that was lost.

Save preflight results, capabilities, package versions, model fingerprints before and after
condition application, and the fully resolved condition.

## Step 1 — Yield and production envelope

`theoretical_yield` establishes reachability and the molar ceiling. `production_envelope`
shows the feasible growth/product trade-off. Report carbon uptake, CO2 uptake, substrate
uptake, and whether the carbon-balance diagnostic raises a warning; a large product flux is
not interpretable without these quantities.

Do not continue when yield is zero. Infeasible envelope points remain in the source table.
Use the same condition object and resolved model state that every later step uses.

## Step 2 — Wild-type reference

Use the configured FBA or deterministic pFBA reference and record wild-type growth, product
flux, and the complete state. MOMA, ROOM, and paired sampling must all refer to this same state
and model fingerprint. An omics comparison may supply the resolved biological condition, but the current
`ProductionWorkflowConfig` does not accept an LAD/E-Flux2 flux state as its reference.

An expression-derived reference would change the biological question and requires a separate
typed workflow boundary before it can be called canonical. Do not mix pFBA, LAD, and E-Flux2
rows in one ranking as though they were replicate measurements.

## Step 3 — Single-knockout predictions

Evaluate the same single-deletion universe independently with:

- MOMA-L2, which minimizes squared adjustment from the wild-type flux state; and
- exact ROOM with the documented flux-prediction tolerance pair, which minimizes the number
  of substantial flux changes.

Keep method-specific tables separate. Every row includes target id and kind, solver status,
growth, product flux, the method-specific objective, and coverage/GPR information. Preserve
lethal and infeasible deletions. A gene with no blocked reaction is *uninformative*, not a
validated neutral knockout.

The primary figure is a two-panel scatter plot with growth rate on the x-axis and target
production on the y-axis, one panel for MOMA and one for ROOM. Use identical limits and units;
mark wild type and the declared viability floor; label the five highest-ranked candidates in
each method panel when five eligible rows are available. A consensus table records agreement
and disagreement without averaging the two models into a fictitious score.

Those method-specific display ranks define the canonical single-knockout candidate universe:
take D1–D5 from MOMA and D1–D5 from ROOM, then deduplicate only equivalent blocked-reaction
signatures. The result contains at most ten candidates. It is a coverage definition, not a
benefit or recommendation filter. Equivalent model phenotypes are simulated once, while every
represented gene id remains in the validation index provenance. Every unique candidate in this
universe proceeds to matched wild-type/knockout sampling in step 6. A candidate with exactly one
blocked reaction also receives a pre-deletion wild-type reaction→product scan: reference↔zero
when reference flux is nonzero, otherwise the full feasible reaction domain as exploratory
response. A multi-reaction signature remains explicit unavailable/skipped because it has no
unambiguous single flux-response axis. Coverage does not depend on which method proposed a
candidate or whether it later passes the positive evidence policy.

Reaction and gene ids are different identifier spaces. Map through the model GPR before
presenting an experimental deletion. For `or` rules, all isozyme genes may need deletion; for
`and` rules, deleting one required subunit may block the reaction. Preserve the rule in the
report instead of reducing it to an unexplained gene list. The canonical export records this
audit trail in `04_single_knockout/gene_knockout_mapping.csv`, including inert genes rather
than quietly treating them as validated neutral deletions. The MOMA and ROOM result tables in
the report join this mapping so each displayed gene includes blocked reaction id(s), reaction
name(s), and the relevant GPR, not only a model-specific gene id.

## Step 4 — Growth-coupled strain design

Run OptKnock and RobustKnock as inverse design methods and export their results separately.
Each design reports knockout set, growth, maximum product, guaranteed product, coupling
verdict, solver, condition, search limits, and the strain-design seed. Set
`strain_design_seed` explicitly in the workflow config (default `0`); the workflow forwards the
same value to both methods. Direct calls likewise pass `seed=` to `optknock` and `robustknock`.
The supported range is the Gurobi-compatible integer range `0..2_000_000_000`.

The seed is part of the computational method, not an incidental runtime setting.
`straindesign` otherwise creates a seed internally, which can send identical MILP formulations
through different search paths and change runtime or the returned solution pool. Never compare
or publish unseeded searches as though they were exact reruns. Record the seed in each method's
metadata and the run-level config/provenance.

Each search also has a wall-clock budget, `strain_design_max_seconds` (default 3600 s;
`null` = none). `straindesign` has no limit of its own and can re-solve one MILP forever,
so an unbounded search is not merely slow. A search cut at the budget is recorded as
incomplete (`straindesign_time_limit_reached`, `straindesign_search_complete=false`), the
artifact is `partial`, and its designs are not ranked as a result. Raising the budget is a
choice to record in the protocol; it does not make a stalled search finish.

Rank by `guaranteed_product`, never by `max_product`. Maximum product is a cooperative upper
bound; guaranteed product is the worst production among growth-optimal states. A design with
zero guaranteed product is not growth-coupled even if its maximum product is large.

Do not treat the number of designs returned by a solver pool as a stable biological result.
Report the configured search limits and ranked designs. Map reaction designs to GPR-resolved
gene interventions before recommending experiments.

This step and the single-knockout step answer different questions. OptKnock/RobustKnock do not
replace MOMA/ROOM, and a beneficial MOMA/ROOM deletion does not prove growth coupling.

## Step 5 — Amplification targets

Run both FSEOF and FVSEOF under the confirmed condition. Export reaction classifications,
ranking quantities, and tidy enforced-product trajectories. Show the top 10 actionable FSEOF
targets and the top 10 actionable FVSEOF targets as independent method-specific lists; do not
filter either list to their intersection. An amplification candidate is a hypothesis about
flux direction, not a prescription for a gene-expression fold change.

For each method-specific candidate, report its own rank, direction, GPR actionability, loop
diagnostic, and step 6 response. FVSEOF robustness and agreement between methods are useful
descriptive columns, but neither substitutes for target-level forward validation and neither is
an admission rule. Keep reversible-reaction direction and boundary flags visible; do not hide
them inside a single opaque score.

For every candidate in the independent top-10 lists, compare standard and fastSNP loopless FVA
capacity under the enforced product and biomass floors. Export the complete diagnostic to
`07_validation/amplification_loop_diagnostic.csv`. The union of both lists is the canonical
amplification candidate universe; every member receives a flux-response index row in step 6.
Keep every method-specific top-10 candidate in its publication trajectory and visibly mark
complete, flagged, failed, or inconclusive diagnostics. A loop-flagged or unresolved candidate
still receives its flux-response scan and retains its diagnostic/eligibility columns, but it is
not promoted to supported evidence or a recommendation and must not be interpreted as
cycle-free. The raw FSEOF/FVSEOF rankings remain unchanged for auditability.

## Step 6 — Forward validation

### Flux response

Create flux-response coverage for the complete canonical candidate universes, with no
downstream selection or recommendation filter:

- every unique MOMA/ROOM D1–D5 single-knockout candidate; and
- every unique reaction in the independent FSEOF top 10 and FVSEOF top 10.

Execute those candidates as follows:

- For amplification targets, force the nominated reaction across its feasible range and use
  product exchange as the response while retaining the configured growth floor.
- For a knockout candidate with one blocked reaction, keep the pre-deletion wild-type
  background and use product exchange as the response under the same growth floor. Scan the
  closed interval between nonzero reference flux and zero. If reference flux is already zero,
  scan the full feasible target-reaction domain and label it exploratory; it is not causal
  support for deletion. Neither scan simulates the complete gene deletion.
- For a knockout signature with multiple blocked reactions, retain an explicit
  unavailable/skipped index row stating that no unambiguous single x-axis exists. Never select
  one blocked reaction silently.

Figure 5 uses the standard response coordinates for every completed scan: enforced candidate-
reaction flux on the x-axis (`target_flux`) and target-product flux on the y-axis
(`response_flux`). Growth is the configured minimum-growth constraint and secondary
`biomass_flux` output, not a Figure 5 axis. MOMA/ROOM and matched knockout sampling, rather than
the pre-deletion reaction scan, establish the complete-knockout phenotype. Show zero-reference
full-domain scans as exploratory, not causal deletion evidence.

Every candidate has a row in `flux_response_index.csv`, and loop flags do not suppress the
numerical scan. Preserve loop status as separate eligibility evidence for later synthesis. If
an unavailable capability, lethal/infeasible background, or analysis failure prevents a
meaningful scan, keep the candidate as an explicit non-complete status row with a reason. It
must not vanish from the index or be replaced by a smaller top subset.
`max_flux_response_targets` is checked before execution as a worst-case capacity guard; it must
never be used to slice this candidate list at runtime.

Report whether the response supports, contradicts, or is inconclusive for the proposed
direction. Preserve infeasible points and phase boundaries. A flat curve does not support an
intervention; an optimum at an artificial bound or a loop-dominated range is a limitation,
not a high-confidence target.

### Paired random sampling for single knockouts

Sampling answers whether a canonical deletion candidate changes the feasible flux
distribution. For every unique MOMA/ROOM D1–D5 candidate, compare a wild-type ensemble with an
ensemble from the knockout model using the same condition, objective conditioning, sampler,
seed policy, sample count, thinning, and reaction set. Do not limit sampling to beneficial,
consensus, recommended, or visually selected candidates.

Export sample-level or auditable summary data for both states and report medians, intervals,
and distribution shifts for product, biomass, the blocked reaction(s), and the mechanistically
relevant reactions. Do not interpret correlated Markov-chain draws as independent biological
replicates or headline an unqualified p-value. A sampler failure or infeasible knockout is a
reported validation result.

Every candidate has a row in `random_sampling_index.csv`. A lethal/infeasible knockout,
unavailable sampler, or failed sampling run remains as an explicit status row with its reason.
Random sampling is targeted to the canonical single-knockout candidate universe in this
workflow. It is not a generic decorative violin plot and does not replace the MOMA/ROOM
prediction.

### Evidence synthesis

The consensus table keeps claim types separate:

| Intervention | Proposal method | Forward check | Growth retained | Product effect | Verdict |
|---|---|---|---|---|---|
| single knockout | MOMA and/or ROOM | MOMA/ROOM + paired sampling; WT reaction scan when single-reaction | explicit | explicit | support / contradict / inconclusive |
| knockout set | OptKnock/RobustKnock | perturbed response | explicit | guaranteed and maximum | coupled / uncoupled / unavailable |
| amplification | FSEOF or FVSEOF | loopless diagnostic + target-to-product response | explicit | direction and range | support / contradict / inconclusive |

Do not collapse these rows into a single score whose meaning changes by method.

`07_validation/recommendations.csv` is a strict positive-evidence subset, not another ranking:
a single-gene deletion needs a selected beneficial MOMA-L2/ROOM result, positive paired-
sampling product shift, and retained growth. A nonzero-reference WT reaction titration may add
mechanistic support but does not replace those complete-knockout analyses; a zero-reference
full-domain scan is exploratory and never causal support for deletion. An
amplification target may originate from either method and needs a supporting response plus a
complete non-flagged loopless diagnostic; a multi-knockout recommendation needs a RobustKnock
design with positive guaranteed product and retained growth. OptKnock remains a required design
table but is not by itself worst-case evidence. Never recommend a combined
knockout-amplification strain because this workflow has not simulated that combined
intervention. Rejected or unvalidated hypotheses remain in their source tables.

## Step 7 — Report and validation

Follow [`_reporting.md`](_reporting.md). Render with:

```python
render_production_report(run_dir, renderer="nature-r")
validation = validate_production_run(run_dir)
```

The run is complete only when validation reports no errors. The report must include:

- the confirmed condition and solver/R capability table;
- method-specific MOMA and ROOM results and their 2D growth/product figure;
- separate OptKnock and RobustKnock tables ranked by guaranteed product;
- independent top-10 FSEOF and top-10 FVSEOF trajectories and classifications;
- loopless-capacity diagnostics, target-indexed flux-response, and paired-sampling evidence;
- links from every numerical claim and figure caption to source CSVs.

The report presents these method-specific results without a **Recommended targets and strain
proposal** section and without a standalone **Limitations** section. It must not promote a row in
the summary, tables, or figures. Put method-specific assumptions and interpretation boundaries
beside the analysis they qualify, and leave intervention selection to the user. Keep
`07_validation/recommendations.csv` as a machine-readable validation artifact for compatibility;
it is not a publication-report section or visual category.

The manuscript renderer writes English labels, 300-DPI raster files, and editable PDF/SVG
counterparts. It fails visibly on missing R packages, a nonzero renderer exit, missing source
data, or incomplete figures.

## Cross-checks

- The config, provenance, pinned model, and all result metadata have the same condition and
  model fingerprint.
- MOMA and ROOM cover the same deletion universe; class counts include lethal/infeasible and
  uninformative rows.
- Every plotted point is present in a CSV, and every D1–D5 display label resolves to one row.
- OptKnock/RobustKnock tables include both maximum and guaranteed product and are ranked by the
  latter.
- The resolved `strain_design_seed` is identical in config, run provenance, and both
  OptKnock/RobustKnock result metadata records; reruns do not rely on a backend-generated seed.
- Every canonical MOMA/ROOM D1–D5 knockout candidate has both flux-response and paired-sampling
  index entries, independent of recommendation status. Single-reaction signatures have a
  pre-deletion WT scan: reference↔zero when reference flux is nonzero, otherwise the full
  feasible domain marked exploratory. Multi-reaction signatures have an explicit
  unavailable/skipped response reason rather than a silently chosen reaction.
- Every candidate in the independent FSEOF/FVSEOF top-10 union has a flux-response index entry;
  loop-flagged/unresolved candidates are still scanned, and any non-runnable case carries an
  explicit status and reason.
- Every recommended amplification target is traceable to its own FSEOF or FVSEOF ranking, has
  a supporting response, and has a completed loopless diagnostic with no artifact flag.
- Paired sampling uses matched settings and the knockout is actually applied in the perturbed
  ensemble.
- Every authoritative analysis CSV links to its method/provenance metadata sidecar, and the
  run-local reproduction, rendering, and validation scripts are manifest-declared.
- Every figure has PNG and PDF/SVG output, a caption, units, and declared source data.
- `00_manifest.json` lists every authoritative artifact and `validate_production_run` has no
  errors.

## Do not

- Do not guess the product exchange, medium, substrate, or oxygen state.
- Do not use MOMA/ROOM as inverse target-finding algorithms; they predict consequences of a
  supplied deletion universe.
- Do not describe an OptKnock maximum as guaranteed production.
- Do not omit the strain-design seed or allow the backend to choose one implicitly.
- Do not merge MOMA, ROOM, OptKnock, RobustKnock, FSEOF, and FVSEOF into one undifferentiated
  ranking.
- Do not sample only wild type and claim that it validates a knockout.
- Do not turn a flux ratio into a wet-lab over-expression fold without an expression/control
  model and experimental design.
- Do not omit unavailable methods, infeasible points, lethal knockouts, or contradictory
  validation from the report.
