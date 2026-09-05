---
id: SC-02
title: Transformation target discovery
goal: Find knockouts that move a source metabolic state toward a target state
when_to_use:
  - "what should I knock out to make the diseased state look like the healthy one"
  - "which perturbation explains the difference between these two expression profiles"
  - "rank drug targets that revert a metabolic phenotype"
  - "find the gene whose deletion turns state A into state B"
role: complete study
requires:
  model: "path to the exact COBRA model"
  source_expression: "expression for the state to move AWAY from"
  target_expression: "expression for the state to move TOWARD"
  condition: "explicit medium, substrate uptake, oxygen/aeration bounds, and other changed bounds"
optional_inputs_from:
  - "a resolved biological condition from an omics comparison"
solver:
  canonical_workflow: "MIQP for MTA and rMTA; there is no LP or QP substitute"
renderer: "nature-r, the same R/ggplot2 path SC-01 uses; PNG/SVG/PDF figures plus linked and self-contained HTML"
steps: [preflight, reference, direction, candidates, transformation, validation]
completion_gate: "validate_transformation_run; a rendered page is not a finished run"
runtime: "minutes on a core model; hours on a genome-scale model, and about 3x that for rMTA"
---

# SC-02 — Transformation target discovery

## Objective

Given expression for two metabolic states, rank the single knockouts most likely to move the
first state toward the second. This is the inverse of SC-01: instead of pushing flux toward a
product, it asks which intervention reproduces an observed difference.

The methods are published:

- **MTA** — Yizhak K, Gabay O, Cohen H, Ruppin E (2013) *Nat Commun* **4**:2632. One MIQP per
  candidate, ranked by their Equation 10 transformation score.
- **rMTA** — Valcárcel LV *et al.* (2019) *Bioinformatics* **35**(21):4350–4355. Three solves
  per candidate — best case, MOMA, and the direction reversed for the worst case — combined by
  their Equation 9 so a candidate that scores well *whichever way it is pushed* is demoted.

Neither originates in CMM. Cite them for any reported ranking.

## Source and target are not interchangeable

**The single most consequential input is which file is which.** Nothing in the model can
detect a swap, and the same pair answers two different questions:

| source | target | the question being asked |
|---|---|---|
| wild type | knockout | which perturbation turns wild type into the knockout? |
| disease | healthy | which perturbation reverts the disease? |

A run with the two exchanged is not a worse answer to the intended question — it is a correct
answer to a different one. Confirm the direction with the user; do not infer it from file
names.

## Canonical execution boundary

```bash
cmm transformation-targets --config CONFIG
```

```python
from cmm.workflows.transformation import (
    TransformationWorkflowConfig,
    TransformationWorkflowResult,
    run_transformation_target_discovery,
)
```

The workflow composes existing CMM services and adds no new numerical method. Use the
individual functions directly only when the user narrows the request to one of them.

## Steps

### 1 — Preflight

Checks the expression files parse, how many replicates each state has, what fraction of the
model's genes the data covers, that the model grows under the requested condition, and that
the solver can run MIQP.

Both expression files must contain finite, non-negative linear measurements. Log-transformed
values are not the workflow input format.

**The MIQP gate is a stop.** `rmta_continuous` is a QP heuristic and explicitly not published
rMTA; it must never stand in for a method the solver cannot run. Report the method as
unavailable instead.

### 2 — Reference state (v_ref)

The source expression's linear replicate means become a flux distribution through **E-Flux2**
or **LAD**.

Everything downstream is measured against v_ref: the sign flip that puts expression labels
into flux-value space, the MIQP's success thresholds `v_ref ± ε`, and the denominator of the
transformation score. Choosing the estimator is therefore a scientific decision, not a
default.

> **Yizhak et al. use iMAT.** CMM implements no iMAT. The substitution is not cosmetic: iMAT
> places no objective on growth, whereas E-Flux2 at `objective_fraction=1.0` forces a
> growth-maximal state, so the ranking is conditioned on whichever estimator is used. State the
> substitution in any report. An externally computed reference state can be supplied through
> the Python API to restore the published pipeline.

### 3 — Direction map

Compares the two states gene by gene, resolves the result through each reaction's GPR (AND
requires all subunits, OR at least one, mixed evidence yields unchanged), and converts it to a
desired direction of flux change using the sign of v_ref.

Each replicate is transformed with `log2(value + 1)` before the comparison. Fold-change
thresholds therefore apply in log2 space; the source reference still uses linear values.
The transform and pseudocount are recorded in provenance.

**With replicates, use the t-test.** That is what the paper specifies, and CMM implements it
as `gene_directions_from_replicates`. Fall back to a fold-change cut only when the data has
one measurement per gene, and say so.

The changed set is then cut to the most differentially expressed reactions — the paper keeps
100–200. The cut is not cosmetic: every changed reaction adds one binary variable to the MIQP,
so on a genome-scale model it decides whether the run finishes.

### 4 — Candidate universe

Dead-end reactions out, essential reactions out (the paper: growth reduced by more than 80%),
then one member per coupled set.

**The candidate count is the denominator of any "ranked in the top *N*%" statement,** so how
it was built belongs in the report. A linear pathway's three reactions are one intervention;
counting them separately inflates the denominator.

Coupled sets are defined on reactions and therefore apply to a reaction-level run. A
gene-level run deduplicates genes that block the same reaction signature, as SC-01's
single-knockout screen already does. Signatures include every reaction blocked by the gene,
including reactions excluded as individual candidates. Essentiality is tested for the full
gene deletion, including joint loss of otherwise nonessential pathways. The construction
metadata preserves all equivalent gene ids under `represented_genes`.

> CMM computes **full** coupling from the null space of S, not the paper's **partial**
> coupling, which needs O(*n*²) linear programmes. Full coupling is stronger, so the grouping
> is conservative — it can split one of the paper's sets but never merge two, which means the
> candidate count it yields is an upper bound on theirs.

The null space always uses the full model. Filtering knockout candidates does not remove
those reactions from the network used to determine coupling.

### 5 — Transformation search

One optimisation per candidate, ranked by the published transformation score. rMTA also
reports its three components (`bTS`, `mTS`, `wTS`) per candidate, because Equation 9's branch
depends on their signs and a reader cannot reconstruct it from `rTS` alone.

**Check the tie structure before quoting a top-*k*.** Ties are broken on target id, so a slice
taken inside a tie block is alphabetical rather than meaningful. The ranking's metadata carries
`n_distinct_scores` and `largest_tie_block`.

### 6 — Validation

**The MOMA baseline.** Yizhak et al. compare MTA against MOMA and report it as *markedly
inferior*; reproducing that contrast is what shows a ranking's signal comes from the method
rather than from the inputs. A ranking whose MOMA baseline agrees with it has not demonstrated
that, whatever its top candidate.

Unsuccessful MOMA solves remain in the baseline CSV with their solver `status` and a `-inf`
sentinel score, after all successful solves. They are listed in the report and excluded from
the rank-comparison figure; an infeasible deletion is a recorded outcome, not a valid score.

**Epsilon sensitivity, when configured.** ε is a flux magnitude with no derivable value here,
so reporting how the ranking moves across ε is the honest substitute for the paper's
derivation. Both source papers report such an analysis.

## Report

`cmm transformation-targets` renders the report unless `--analysis-only` is given; the same
page is available as `cmm.reporting.render_transformation_report(run_dir, highlight=…)`.

Figures come from `render_transformation_figures.R`, the same checked-in R/ggplot2 path SC-01
uses; the two workflows draw different panels, so each has its own script. Three panels are
drawn — score against rank, transformation rank against the MOMA baseline, and rank against
epsilon when a sweep was configured — each as PNG, SVG and PDF, with two copies of the page
beside them:

| file | figures | use |
|---|---|---|
| `report.html` | linked as `figures/*.png` | reading the run in place |
| `report_standalone.html` | embedded as data URIs | sending to someone |

**Send the standalone copy.** The linked one renders every figure as blank space once it leaves
the run directory, and nothing on the page says so.

### The completion gate

```bash
cmm report validate RUN_DIR --json
```

`validate_transformation_run` is what says the run finished. A page that opens is not evidence:
it can point at figures that are gone, quote a CSV that was edited after the run, or be stale
with respect to the figures beside it. Each of those looks like success in a browser. The gate
checks every declared artifact against its recorded hash and size, that the ranking is ordered
1..N by descending score, that a skipped stage records why, that every rendered figure has a
non-empty PNG, SVG and PDF, and that the standalone page carries every image it references.

The bundle includes byte-for-byte copies of both expression inputs under `inputs/` and the
source model under `model/`. Both `00_config.json` and `scripts/transformation_config.json`
refer to these files with relative paths. The gate checks input hashes and these references;
`scripts/reproduce.py` can therefore replay a moved bundle without the original external
files. Older bundles missing these inputs must be regenerated before they pass this gate.

`cmm report render` and `cmm report validate` read the run's own `00_manifest.json` to decide
which workflow it belongs to, so neither needs to be told.

**R is required to render, not to analyse.** `Rscript` and the renderer packages are needed for
the report; `--analysis-only` produces the full run bundle without them.

## What the report must state

The renderer writes all of these into the page body rather than leaving them in a provenance
file. Beyond the standard contract in [`_reporting.md`](_reporting.md):

1. **Which file was source and which was target**, in words, not just as paths.
2. **How v_ref was computed**, and that it is not the published iMAT-plus-sampling state.
3. **ε, and that it was chosen rather than derived.** `TransformationWorkflowConfig.suggest_epsilon`
   reads percentiles off v_ref so the choice has a basis.
4. **The candidate construction and its count** — the denominator of any percentile claim.
5. **The tie structure** of the ranking.
6. That predictions are *in silico* hypotheses requiring experimental validation. rMTA is a
   prioritisation robustified by its worst-case term, not a proof.

## Interpretation limits

- A high-ranking knockout is a hypothesis about which intervention *could* produce the observed
  difference, not evidence that it *did*.
- The ranking is conditioned on one medium and one reference state. Treat a different condition
  as a separate run.
- Reactions ranked near the true answer are often the ones sharing its consequences. That is
  the method working, not a failure, but it means the top of the list is a neighbourhood rather
  than a single call.
