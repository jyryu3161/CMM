# Scientific validation and reproducibility

This document defines what CMM 0.5.0 has been validated against, how to reproduce the
evidence, and what the tests do not establish. A green test suite supports numerical and
implementation correctness; it is not a substitute for experimental validation of a new
biological target.

CMM source code, documentation, and test data are freely available under the MIT License at
<https://github.com/jyryu3161/CMM>. For a Bioinformatics Application Note submission, archive
the exact tagged release in Zenodo or an equivalent long-term repository and add its DOI to
the README, `CITATION.cff`, and manuscript Availability and Implementation statement.

The two optional repository skills are interfaces to the same numerical workflows, not sources of
solver evidence. When AI assistance materially affects development, analysis, interpretation,
or manuscript preparation, record and disclose it separately as described in
[`AI-USAGE.md`](AI-USAGE.md).

## Reproduce the publication environment

The repository tracks `uv.lock`, including hashes and resolutions for every supported Python
version. From a clean checkout:

```bash
uv sync --frozen --all-extras
QT_QPA_PLATFORM=offscreen uv run pytest -q -ra --strict-markers \
  --durations=10 --cov=cmm --cov-branch --cov-report=term-missing --cov-fail-under=80
uv run ruff check src tests
uv run ruff format --check src tests
uv run mypy src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
Rscript --vanilla -e 'if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org"); renv::restore(prompt=FALSE)'
R_LIBS_USER="$(Rscript --vanilla -e 'renv::load(); cat(.libPaths()[1])')" Rscript --vanilla -e 'p <- c("jsonlite","ggplot2","ggrepel","patchwork","svglite","ragg"); stopifnot(all(vapply(p, requireNamespace, logical(1), quietly=TRUE)))'
uvx --from cffconvert==2.0.0 cffconvert --validate
uv export --frozen --all-extras --no-emit-project -o /tmp/cmm-requirements-audit.txt
uvx --from pip-audit==2.10.1 pip-audit \
  -r /tmp/cmm-requirements-audit.txt --disable-pip --require-hashes
uv build
uvx twine check dist/*
```

The R commands are optional for Python-only use and require the matching `renv.lock`. Git
checkouts, GitHub source archives, and the sdist include that lock; the Python wheel does not.
Both Python distributions contain the R renderer itself. A wheel-only user can install compatible
renderer packages, but exact publication reproduction requires the lock from the same release.

The general Python matrix executes the locked solver-neutral/GLPK environment on Linux
(Python 3.10, 3.11, and 3.12), Windows 3.12, and macOS 3.12. It deselects exactly the tests
marked `requires_qp` or `requires_miqp`; a separate Linux 3.12 job installs Gurobi plus the
desktop dependencies and runs that complete marked universe with offscreen Qt. It also runs the
exact unmarked MILP node that solves both OptKnock and RobustKnock and checks deterministic seed
forwarding. The current audited solver collection therefore contains 73 marked tests plus that
one explicit MILP test, with no exception. A collection/JUnit identity guard fails if an audited
node is omitted, skipped, duplicated, or replaced without review.

The solver job also restores locked R because the SC-02 tests render and validate real MIQP
run bundles. Those workflow-contract tests use the six-reaction branched model, keeping both
MIQP and MOMA within the restricted license's 200-variable quadratic limit. Genome-scale
numerical checks remain separate from these orchestration and report tests.

A separate publication-report matrix restores `renv.lock`, prints and checks every locked
package version, and runs `test_publication_reporting.py` and `test_transformation_reporting.py`
on Ubuntu, macOS, and Windows. The
lock records R 4.3.2, so archived package versions may compile when a matching binary is
unavailable: Linux and macOS jobs install the graphics build libraries, and the Windows job
provisions matching Rtools43. `ggplot2` itself declares `NeedsCompilation: no`; the relevant
compiled packages are `ragg`, `systemfonts`, `textshaping`, and `svglite`. The matrix permits
both repository binaries and source builds and tests the restored result; it does not force a
source-only installation. The build toolchains prepare and support source fallback without
claiming that every matrix run actually compiled a package.

The 80% branch-coverage requirement applies to the union of the public runtime surfaces, not to
each capability-limited shard in isolation. The Linux 3.12 general job, restricted-license
solver job, and locked-R publication job each upload their coverage data only after their tests
pass. A final coverage job combines those three datasets and enforces `--fail-under=80`. It also
depends on the complete cross-platform matrices, so combining coverage does not make an OS test
optional. This partition preserves the solver-neutral matrix and exact R lock checks while
preventing expected QP/MIQP deselection or R-package skips from being misreported as a product-
wide coverage regression.

The current `renv.lock` records exact versions and CRAN sources but has no per-package `Hash`
fields. A clean renv 1.1.4 restore/snapshot audit did not provide a safe Hash-only round trip:
the canonical snapshot omitted records that were not fully materialized in its isolated
library. The lock therefore remains hashless rather than carrying hand-generated hashes.
Regenerate hashes only from a complete R 4.3.2 project library and accept the result only when
the package set and every R/package version are unchanged; until then, the three-OS complete
version comparison is the explicit CI guard. That guard also parses every `Depends`, `Imports`,
and `LinkingTo` entry, requires each hard dependency to be locked or supplied by base/recommended
R, and requires the restored project's non-base package set to match the lock exactly. `renv`
itself is locked, so CI permits no additional unpinned bootstrap package.

QP/MIQP checks use small models that fit the bundled Gurobi restricted license
(largest QP: 190 variables / 73 constraints, against that license's 200-variable limit for
quadratic models); all 60 capability-marked nodes are required tests, not environment-dependent
skips. The additional audited MILP node exercises OptKnock, RobustKnock, and explicit seed
forwarding under the same license. That 10-variable QP margin is deliberate and narrow —
enlarging the E-Flux2 test model beyond 100 reactions would make these tests require a full
license. Genome-scale LP validation uses GLPK. Pytest's ten slowest durations are retained in
every CI log as a lightweight runtime regression record; solver/model/license details remain
part of result provenance.

## Evidence matrix

| Area | Independent/reference evidence | Automated checks |
|---|---|---|
| FBA, pFBA, FVA | Direct COBRApy calls on `e_coli_core` and iJO1366 | `test_validation.py`, `test_genome_scale.py` |
| MOMA | Direct COBRApy L2 MOMA on an independently constructed branched network | `test_validation.py` |
| E-Flux2 | Independent two-stage QP: maximize biological objective, fix it, minimize total squared flux | `test_validation.py`, `test_expression.py` |
| Yield and media | Direct objective/bound calculations in COBRApy | `test_validation.py`, `test_production.py`, `test_media.py` |
| FSEOF/FVSEOF | Enforced product levels, biomass optimization/FVA, regression and boundary filtering | `test_production.py`; scan-resolution sensitivity in `test_scientific_sensitivity.py` |
| MTA/rMTA | Official COBRA Toolbox test topology, expected score signs, published Equation 9 | `test_revert.py`; alpha sensitivity in `test_scientific_sensitivity.py` |
| OptKnock/RobustKnock | Distinct StrainDesign module types (`OPTKNOCK` and three-level `ROBUSTKNOCK`), explicit seed forwarding/provenance, plus post-solve max/min product evaluation | `test_strain_design.py` |
| Flux response | Scanned optimum reproduces the known `e_coli_core` oxygen growth optimum; competing-branch model gives an exact −1 response gradient | `test_response.py` |
| Flux sampling | Drawn samples satisfy `S · v = 0` and every model bound; identical seeds reproduce identical ensembles | `test_sampling.py` |
| Production-target workflow | Nested JSON config/path resolution, deterministic unique-gene ranking, strict QP gate, complete candidate-universe response/sampling coverage, candidate-reaction response semantics with explicit multi-reaction unavailability, paired-sampling distance, and schema-v2 export | `test_production_workflow.py` |
| Publication reporting | Schema-v2 path/column validation, path-escape rejection, deterministic standalone HTML, R-rendered PNG/PDF/SVG, and explicit optional-panel status | `test_publication_reporting.py` |
| GUI/state | Real offscreen Qt workflows, invalid-file rejection, model reload state invalidation | `test_app_smoke.py`, `test_scenarios.py` |
| Provenance | Deterministic SHA-256 fingerprint changes with model bounds and accompanies numerical results | `test_core_primitives.py` and feature tests |

## Genome-scale check

The non-optional genome-scale fixture loads `iJO1366` through COBRApy's model repository and
asserts the imported artifact's exact dimensions (2,583 reactions and 1,367 genes) before any
comparison. CMM FBA, pFBA, and selected-reaction FVA must match direct COBRApy results; a
succinate FSEOF scan must complete without failed levels and must not return biomass or the
product exchange as an actionable intervention.

The reconstruction derives from Orth et al. (2011),
<https://doi.org/10.1038/msb.2011.65>. Exact imported dimensions include repository-level
boundary and annotation representation, so the dimension assertion also detects a changed
upstream artifact.

## Method contracts

### Gene-protein-reaction resolution

Two different GPR rules ship, because the source papers specify two different operations, and
each result records which one produced it.

- **Continuous expression values** (`gene_to_reaction_weights`, and therefore E-Flux2 and LAD):
  `AND` is `min`, `OR` is `or_rule`, defaulting to `"sum"`. Kim et al. (2016) take *"the sum of
  the expression values of the associated genes"* for isozymes; Lee et al. (2012) state *"the
  total capacity is given by the sum of its components"*. `or_rule="max"` is CMM's pre-0.4.0
  behaviour, matches no source paper, and is retained only to reproduce an old result. The rule
  is archived as `gpr_or_rule` in the result's provenance.
- **Ternary direction labels** (`reaction_directions`, and therefore MTA/rMTA): Yizhak et al.'s
  rule — all genes elevated or all reduced under `AND`, at least one under `OR`, and **mixed ⇒
  unchanged**. It is implemented as a two-pass binary decomposition, not by reusing the numeric
  `min`/`max` helper, and is recorded as `gpr_rule = "yizhak2013_two_pass_binary"`. A sum is
  not type-correct over direction labels, so the two rules are deliberately not shared.

### E-Flux2

CMM normalizes non-negative expression-derived reaction weights, constrains reaction bounds,
maximizes the model's biological objective, enforces the requested fraction of that optimum,
and minimizes the L2 norm of all fluxes. The default objective fraction is `1.0`. Without QP
support the method raises; the opt-in `allow_l1_fallback=True` result is labeled
`eflux2_l1_fallback` and must not be reported as E-Flux2.

Reference: Kim MK, Lane A, Kelley JJ, Lun DS (2016), *PLoS ONE* 11(6):e0157101,
<https://doi.org/10.1371/journal.pone.0157101>. **Cite that author list.** "Kim, Woo, Choi" is
a wrong author list for this paper and must not appear in any manuscript or docstring.
E-Flux2 builds on E-Flux, Colijn et al. (2009), *PLoS Comput Biol* 5(8):e1000489,
<https://doi.org/10.1371/journal.pcbi.1000489>.

### LAD

LAD fits fluxes to expression-derived reaction targets by minimizing absolute deviation, an LP.
It is the method named "Lee-12" in the Machado & Herrgård (2014) survey of expression-integration
methods.

The residual is taken against the **absolute** flux, using COBRApy's forward/reverse split
(`v = f − b`, `f, b ≥ 0`, with `b` pinned to 0 when `lb ≥ 0` and `f` pinned to 0 when `ub ≤ 0`),
so the constraint is `f + b − d⁺ + d⁻ = target`. Fitting the signed flux instead — which CMM did
before 0.4.0 — assessed a maximum penalty on a correct answer: a reversible reaction at `v = −5`
with target 5 scores deviation 0 in the source formulation and 10 in the signed one.

Two deviations from Lee et al., both recorded in `metadata["cmm_deviations"]` and neither
silent. (i) The per-reaction `1/σ` weights of Eq. 2/3 are **off by default**; supply them with
`reaction_sigma=`, and `metadata["sigma_weighted"]` records which was used. (ii) Fitting `f + b`
in a single LP does not force `f · b = 0`, so a reversible reaction may satisfy a target by
carrying flux in both directions; that is a property of the reference formulation too, which
avoids it with an FVA-driven iteration CMM does not implement. `weight_threshold` defaults to
`0.0` from 0.4.0 (previously 0.01), so a low-expression reaction is driven toward zero flux
rather than excluded from the objective.

Reference: Lee D, Smallbone K, Dunn WB, Murabito E, Winder CL, Kell DB, Mendes P, Swainston N
(2012), *BMC Systems Biology* 6:73, <https://doi.org/10.1186/1752-0509-6-73>. Survey placement:
Machado & Herrgård (2014), *PLoS Comput Biol* 10(4):e1003580,
<https://doi.org/10.1371/journal.pcbi.1003580>.

### FSEOF and FVSEOF

Scans begin at product flux observed under maximum growth and end at a selected fraction of
the theoretical maximum. FSEOF fixes product flux and maximizes biomass at each level. FVSEOF
fixes both product and an exact fraction of maximum biomass, computes reaction FVA, and reports
midpoint (Park's `V_avg`), forced-minimum magnitude, range width (Park's `l_sol`), and their
slopes. `n_steps` defaults to 10, Park et al.'s stated minimum. Diagnostic tables contain all
requested reactions; actionable lists exclude boundary/objective/no-GPR reactions.

FVSEOF classification follows Park et al.'s **nine types**, taken from the joint sign of ΔV_avg
and Δl_sol: types 1–3 (V_avg rising) are the amplification candidates, 4–6 the knockdown
candidates, 7–9 neither. `amplification_targets()` returns types 1–3 ordered by ascending mean
`l_sol`, which is Park's stated priority (*"reactions with smaller values of l_sol received
higher priorities"*). The type index is exported as `park_type`. The index *within* each band
is CMM's convention and is labelled as such: the paper's own numbering is in a figure that
could not be resolved from the text, and no numbering was invented for it.

The optional `linear_flux_couplings` keyword (renamed from `group_constraints` in 0.4.0) takes
caller-supplied linear equalities `Σc·v = 0`.

Two selection rules here are CMM's own and must not be attributed to the source papers.
FSEOF's criterion (endpoint difference, positive linear slope, no sign reversal, baselined at
the 10% scan level) is stricter than Choi et al.'s `|v_j|max > |v_j^initial|` ∧
`v_j^max · v_j^min ≥ 0`. FVSEOF's `robust_targets()` flag — forced FVA minimum rising
monotonically — has no counterpart in Park et al. and is not their variability criterion; Park's
variability signal is Δl_sol, which the nine types encode. FVSEOF's `linear_flux_couplings` are
caller-supplied linear equalities `Σc·v = 0`, a different object from Park et al.'s
STRING-derived on/off reaction pairs with `|v1−mean| + |v2−mean| ≤ δ` on uptake-normalised
fluxes, none of which CMM implements. FVSEOF's `slope` regresses `|V_avg|`, matching CMM's
magnitude convention for reverse-running reactions, not Park's signed `q_slope`.

References: Choi HS, Lee SY, Kim TY, Woo HM (2010), *Appl Environ Microbiol* 76(10):3097–3105,
<https://doi.org/10.1128/AEM.00115-10>; Park JM, Park HM, Kim WJ, Kim HU, Kim TY, Lee SY (2012),
*BMC Systems Biology* 6:106, <https://doi.org/10.1186/1752-0509-6-106>.

### OptKnock and RobustKnock

OptKnock calls StrainDesign's optimistic two-level module. RobustKnock calls its distinct
three-level module, which guards the minimum product flux among growth-optimal states. CMM
then independently fixes each design at its optimal growth and reports both maximum and
minimum product flux. A nonzero `guaranteed_product` is therefore an evaluated property, not
an alias of the optimistic value.

**Candidate sets are restricted by default.** `actionable_only=True` limits knockout candidates
to gene-associated internal reactions, because a design deleting a boundary exchange with no GPR
is not realisable as a gene deletion; Burgard et al. likewise restrict candidates to central
metabolism. On `e_coli_core` this excludes 26 of 95 reactions and removes 10 of 18 returned
designs without changing the top design or its numbers. `actionable_only=False` restores the
unrestricted set. `max_solutions` caps MILP solutions, not distinct designs, and designs are
de-duplicated by knockout set before being returned. Both functions take `condition=` as of
0.4.0; before that they accepted no context at all.

**The MILP search seed is explicit.** `optknock(..., seed=0)` and
`robustknock(..., seed=0)` forward the integer to `straindesign`; accepted values are
`0..2_000_000_000`, while booleans and floats are rejected. Without an explicit value,
`straindesign` generates a new random seed, which can alter the MILP path, runtime, and returned
solution pool between calls. CMM therefore defaults to `0` and records the value in
`metadata["seed"]`, `metadata["parameters"]["seed"]`, and
`metadata["parameters"]["strain_design_seed"]`. Reproduction still requires the recorded
solver/version/platform because a seed alone does not guarantee identical behavior across
different MILP implementations or versions.

**Three citations are required together, not two.** The bilevel and three-level searches are
solved by the `straindesign` package, which carries no literature citation of its own; citing
only the original method papers would attribute that package's candidate set, network
compression, decompression, and solution enumeration to the original authors.

References: Burgard AP, Pharkya P, Maranas CD (2003), *Biotechnol Bioeng* 84(6):647–657,
<https://doi.org/10.1002/bit.10803>; Tepper N, Shlomi T (2010), *Bioinformatics* 26(4):536–543,
<https://doi.org/10.1093/bioinformatics/btp704>; Schneider P, Bekiaris PS, von Kamp A, Klamt S
(2022), *Bioinformatics* 38(21):4981–4983, <https://doi.org/10.1093/bioinformatics/btac632>
(the StrainDesign package itself).

### MOMA and ROOM

Both are delegated to COBRApy; nothing is reimplemented, so their numerical behaviour is the
ecosystem's — including ROOM's constraint capping the perturbed objective at the reference
objective, which is COBRApy's addition and is not in Shlomi et al.'s Eqs 1–3.

**L2 MOMA is Segrè et al. (2002); the L1 variant is not.** That paper contains no linear
variant, so `moma_l1` must be described as the linear variant as implemented in COBRApy, with
no citation to Segrè et al., until a quotable source for it is obtained.

ROOM's δ = 0.03 / ε = 0.001 are Shlomi et al.'s **flux-prediction** pair; the paper gives
δ = 0.1 / ε = 0.01 for **lethality prediction**. Both ship as `cmm.features.ROOM_TOLERANCES` and
are selected by `use_case` / `room_use_case`: `room()` and `knockout_comparison()` default to
flux prediction, `batch_comparison` to lethality, and every result records the pair it used.
Explicit `delta`/`epsilon` override the preset. The choice is not cosmetic — 531 against 401
total switches over the same 35 gene knockouts, a 24% shift in the ranking score.

**Every result in this family carries the full `run_provenance` block**, alongside the
reference identity and (for ROOM) the tolerance pair. `moma`, `room` and `knockout_comparison`
fingerprint the model *as handed to the solver*, so a knockout's record is the fingerprint of
the knocked-out model; `batch_comparison` carries one block for the screen on its
`BatchComparisonResult` container, fingerprinting the model before any knockout — the one
model every row derives from — together with what the enumeration dropped
(`n_inert_dropped`). `seed` is recorded as `null`, because MOMA and ROOM are deterministic and
CMM does not invent a seed where a method has none.

**The reported quantity is separated from the solver objective.** `objective_value` is the raw
objective and differs in kind per method: `Σd²` for `moma_l2`, `Σ|d|` for `moma_l1`, and a
*count of significantly changed reactions* for ROOM. `distance` is a distance and only a
distance — the Euclidean `√(Σd²)` of Segrè et al.'s Eq. (4) for `moma_l2`, the L1 sum for
`moma_l1`, and `None` for ROOM, whose count is exported as `n_changed_reactions`.
`distance_kind` records which. Before 0.4.0 one field held all three; on the SC-01 design it
reported 1303.99 where the Euclidean distance is 36.11.

**The reference state is pFBA in both the Python API and the GUI.** Segrè et al. specify FBA,
and `reference_flux(model, "fba")` remains available; pFBA is chosen for determinism and is a
disclosed deviation. Measured: FBA vs pFBA changes MOMA growth by at most 0.00096 h⁻¹ (0.098%
of wild type) with 0 of 24 knockouts reclassified on *i*JO1366 and byte-identical vectors on
`e_coli_core`, while 8 identical `reference_flux(..., "fba")` calls on *i*JO1366 returned 3
distinct vectors. Every result records the reference state it was measured against.

References: Segrè D, Vitkup D, Church GM (2002), *PNAS* 99(23):15112–15117,
<https://doi.org/10.1073/pnas.232349399>; Shlomi T, Berkman O, Ruppin E (2005), *PNAS*
102(21):7695–7700, <https://doi.org/10.1073/pnas.0406346102>.

### Flux response

`flux_response` fixes the target reaction at each point of a linear scan and maximizes the
response reaction, recording the response flux, biomass flux, and solver status per point.
The default scan range is the target's full feasible interval (FVA at a zero fraction of the
optimum); an explicit range may exceed the reaction's declared bounds and is then flagged
`range_outside_bounds` in provenance. `response` defaults to the objective reaction, giving
the robustness reading; naming a product exchange gives the production reading. Because
maximizing a product without a growth floor returns non-growing solutions, `biomass_fraction`
holds biomass at a fraction of the wild-type optimum across the scan. Infeasible points are
retained with NaN fluxes and their solver status: the flux at which a scan stops being
solvable is a reported result, not an error. `feasible_range()` is computed by FVA on the
scanned reaction under the applied growth floor, not read off the scan grid.

Sensitivity is reported as the **shadow price** `d(response)/d(target)`, taken exactly from the
LP dual, together with the phase boundaries at which it changes. The response curve is the
optimal-value function of an LP parameterized in one bound and is therefore concave piecewise
linear, so a finite-difference "steepest decline" locates the edge of the scan grid rather than
a property of the network. The `bottleneck` field that reported it was removed in 0.4.0: its
location moved by up to 29.53 flux units and its `found` flag inverted as `n_steps` went from 6
to 160, and no published criterion in constraint-based modelling defines a bottleneck as the
argmin of a finite-difference slope. Regions of constant shadow price and the boundaries
between them are the published objects.

References: Lee KH, Park JH, Kim TY, Kim HU, Lee SY (2007), *Mol Syst Biol* 3:149,
<https://doi.org/10.1038/msb4100196> (PMC2174629, PMID 18059444 — **not** PMC1911197, which is
Feist et al. 2007), for in silico flux response analysis; Edwards JS, Palsson BØ (2000),
*Biotechnol Prog* 16(6):927–939, <https://doi.org/10.1021/bp0000712>, for the biomass case
(robustness analysis); Edwards JS, Ramakrishna R, Palsson BØ (2002), *Biotechnol Bioeng*
77(1):27–36, <https://doi.org/10.1002/bit.10047>, "Characterizing the metabolic phenotype: a
phenotype phase plane analysis", for phases as regions of constant shadow price. **Orth, Thiele
& Palsson (2010) contains no definition of shadow price or reduced cost and must not be cited
for either.**

### Random flux sampling

`random_flux_sampling` draws flux distributions with COBRApy's OptGP or ACHR samplers under
the model and condition as constrained. Every run takes an explicit `seed` and defaults to
`processes=1`, because parallel chains are seeded independently and a multi-process run is
therefore not bit-for-bit reproducible; provenance records `reproducible` accordingly.
`reference_constrained_sampling` first narrows each reaction to `[v * min_fraction,
v * max_fraction]` around a reference flux state (mirrored for negative `v`, and
`[-zero_window, zero_window]` for reactions at essentially zero flux), intersected with the
existing bounds. Since `min_fraction <= 1 <= max_fraction`, every window contains its
reference flux, so an empty window can only mean the reference violates the model's bounds —
which raises rather than silently sampling a different space.

Samples satisfy the steady-state and bound constraints, and `test_sampling.py` asserts
`S · v = 0` on every drawn sample. CMM deliberately provides no post-hoc noise addition:
perturbing sampled fluxes independently would break both constraints, so the perturbed
vectors would no longer be flux distributions.

References: Megchelenbrink et al. (2014), <https://doi.org/10.1371/journal.pone.0086587>
(OptGP); Kaufman and Smith (1998), <https://doi.org/10.1287/opre.46.1.84> (ACHR).

### MTA and rMTA

The optimization, transformation score, and robust Equation 9 are documented in
[design-revert-metabolism.md](design-revert-metabolism.md). The GUI uses deterministic
E-Flux2 source-state preprocessing, which differs from the original contextualization plus
sampling protocol and must be disclosed in manuscripts.

Two of Yizhak et al.'s preprocessing parameters are exposed rather than assumed, and a run
that intends to follow the paper must set both:

- **ε**, the flux change that counts as significant. `revert.DEFAULT_EPSILON` is a fixed
  scalar; the paper derives ε per data set from a sampled reference distribution. State the
  value used.
- **the changed-set size.** The paper keeps only the top 100–200 most differentially expressed
  reactions as "changed"; pass `differential_expression(top_n_changed=…)` or call
  `restrict_to_top_changed()`. Left unset, every reaction clearing the fold-change thresholds
  stays changed, which is not the published protocol and adds one MIQP binary variable each.

The published significance test is a Student's t-test over replicates, implemented as
`gene_directions_from_replicates`. It is what the SC-02 workflow and the GUI's Revert
Metabolism tab run whenever both expression files carry at least two replicate columns. Data
that cannot support it fall back to a fold-change cut — `gene_directions` on a pair of vectors,
or `gene_directions_by_fold_change` on the same replicate frames — which is not the published
test, and a run that used it must say so.

SC-02 accepts finite, non-negative linear expression measurements. It keeps linear replicate
means for reference inference and applies `log2(value + 1)` to each replicate before direction
tests, matching the desktop comparison. Regression tests distinguish a small relative change
(`100 → 102`) from a two-unit log2 change and check the replicate-level t-test transform.

References: Yizhak K, Gabay O, Cohen H, Ruppin E (2013), *Nat Commun* 4:2632,
<https://doi.org/10.1038/ncomms3632>; Valcárcel LV, Torrano V, Tobalina L, Carracedo A,
Planes FJ (2019), *Bioinformatics* 35(21):4350–4355,
<https://doi.org/10.1093/bioinformatics/btz231>.

### Transformation targets

`transformation_targets` is **not a CMM invention.** Both of its paths map to published Yizhak
et al. (2013) methods: the `mta` path is that paper's MIQP formulation applied to an arbitrary
source→target pair rather than a disease→healthy one, and the `moma` path is the
distance-reduction scoring the same paper uses as its comparison method. Cite Yizhak et al.
(2013) for both; do not describe either as a method originating in CMM.

Note that Yizhak et al. report the MOMA-style scoring as *markedly inferior* to MTA, and it is
CMM's default. State which path produced a ranking.

Regression tests cover failed MOMA solves with negative-scoring viable alternatives: failures
retain their status and a `-inf` sentinel, follow successful targets, and are excluded from the
report's rank comparison. Candidate tests retain the full network for null-space coupling and
the full GPR deletion signature for gene grouping, including jointly essential deletions.
Bundle tests remove the original expression files, relocate the run, replay it from archived
inputs, and check identical rankings. Input corruption, missing input artifacts, and external
paths in replay configs fail validation.

### Methods with no published source

`flux_log_change` is a **CMM utility**: a log2 ratio of two flux vectors with a pseudocount. It
implements no published method and **must not be cited to any paper.** Describe it as CMM's own
reporting convenience, and state the pseudocount, which determines the value returned for
switch-on and switch-off reactions.

`network_flux_map`'s schematic layout is likewise **CMM's own**: it selects the highest-`|flux|`
reactions, drops a fixed currency-metabolite list, and runs a seeded force-directed layout. It
is a flux overview, not a curated map, and implements no published layout algorithm. The
curated alternative is not CMM's either — `escher_flux_map` renders a map authored in Escher,
and work using that layout should cite King et al. (2015); the bundled map's provenance,
SHA-256 and license are recorded in `src/cmm/resources/ATTRIBUTION.md` and asserted by
`tests/test_resources.py`.

### JEV agent design

`run_jev_design` implements **no published method** and must not be cited to one. It is a
search loop whose proposal step is TypeSafe's JEV decision model and whose every measurement
is an ordinary CMM solve. What is validated is the loop, not the agent:

- Each applied move becomes bounds on one reaction, and the resulting state is a `pfba` solve
  and (unless disabled) a `moma` solve against the wild-type pFBA reference. Both carry the
  contracts documented above; the agent adds no numerics of its own.
- **Viability is enforced by CMM, not predicted.** A move leaving the model infeasible, or
  growth below `growth_floor`, is reverted and recorded. `tests/test_jev.py` asserts this
  against a scripted agent that repeatedly proposes a lethal knockout.
- **The move vocabulary is down-regulation only**: `knockout` (bounds forced to zero) and
  `knockdown_50` (magnitude capped at half the wild-type flux). Both are caps, which is what a
  gene deletion or a weakened promoter imposes. Amplification was removed because a forced
  lower bound makes the solver carry a flux by whatever route is cheapest, which is not what
  raising an enzyme's expression does; an in-silico gain obtained that way overstates what a
  strain would do. Asserted by `test_nothing_in_the_vocabulary_can_force_flux_up`.
- **The cost of that restriction is measured, not assumed.** Every comparison carries a
  `best amplification on top of this design` row: FSEOF ranked on the design being scored, its
  top ten targets forced one at a time, the best one clearing the growth floor reported. It is
  excluded from the "best deterministic method" the agent is scored against, because scoring
  an agent against a move it was forbidden to play is not a comparison. The forced level uses
  the reaction's **loopless** feasible maximum when it carries no flux: the plain LP maximum
  is a loop artifact on `e_coli_core` (`FRD7` returns its 1000 bound through the
  `FRD7`/`SUCDi` cycle at 10 mmol gDW⁻¹ h⁻¹ glucose uptake).
- **MOMA's reference is the wild-type pFBA state for every design, at every step**, never the
  previous step's state. MOMA's premise is that a freshly perturbed cell keeps its parent's
  regulatory setpoints, and the parent of the design a run produces is the wild type — the
  design is built and characterised as one strain, and its intermediate states are search
  positions rather than organisms. Chaining the reference would make each step's distance
  describe only the last edit, so a five-edit design would appear as easy to build as a
  one-edit design. SC-03 states this in the run's own documentation.
- Candidates are restricted to reactions carrying a gene association, so a proposal is an
  intervention a laboratory could make rather than a bound edit that only exists in silico.
- **Knockouts and knockdowns are budgeted separately** (`max_knockouts`, `max_knockdowns`),
  because they are different things for a laboratory to build and a single shared cap is not
  a quantity anyone plans against.

**Reproducibility is partial and must be stated as such.** Every CMM solve in a run is
deterministic. JEV's decisions are not guaranteed to repeat: repeated runs of the same
configuration have produced designs differing several-fold in product flux.
`04_agent/transcript.jsonl` records every request and response so a single run can be
audited; it does not make the method reproducible. A single run is evidence about that run.
Claiming agent performance requires repeated runs and a stated distribution, which this
release does not provide, and comparing it against FSEOF or OptKnock requires running those
on the same model and condition, which SC-03 does not do.

Predictions are computational hypotheses for experimental test, exactly as elsewhere in CMM.

### Production workflow and publication-report contract

`run_production_target_discovery` is validated as orchestration, not as a new metabolic
algorithm. Its numerical evidence comes from the method contracts above. The workflow-level
contract is that it:

- requires a file-backed SBML model, product exchange, and serializable resolved config;
- applies one condition before every analysis and records source and conditioned model
  fingerprints;
- gates QP for MOMA-L2 and MILP for exact ROOM before the matched single-gene screens, with no
  L1 or relaxed-ROOM fallback;
- screens the same GPR-resolved gene universe with both methods, retains non-optimal rows, and
  records separate deterministic D1–D5 display ranks and downstream positive-benefit ranks;
- gates MILP plus importable `straindesign` for OptKnock/RobustKnock and surfaces any
  additional requirement reported by the selected backend;
- forwards the explicit workflow `strain_design_seed` to both design methods and records it in
  resolved config, run provenance, and method metadata; it never permits a hidden
  backend-generated random seed;
- exports independent top-10 FSEOF and top-10 FVSEOF rankings and tidy trajectory tables
  separately, without using their intersection as a selection gate;
- compares standard and fastSNP loopless FVA capacity for publication amplification candidates,
  preserves all independent top-10 rows in the tidy publication trajectories with visible
  diagnostic status, runs flux response even for flagged or unresolved targets, and withholds
  them only from forward-validation support or recommendations;
- defines the canonical single-knockout validation universe as the unique blocked-reaction
  signatures represented by MOMA D1–D5 and ROOM D1–D5, then applies matched wild-type/knockout
  sampling to every candidate in that universe; a signature with one blocked reaction receives
  a pre-deletion wild-type reference↔zero scan when reference flux is nonzero, or a full feasible-
  domain exploratory scan when it is already zero; a multi-reaction signature remains explicit
  unavailable/skipped because no unambiguous single response axis exists, and no blocked
  reaction is selected silently;
- applies flux response to every unique candidate in the independent report-visible FSEOF
  top-10 and FVSEOF top-10 lists, with no recommendation or beneficial-selection gate;
- uses `max_flux_response_targets` only as a preflight capacity guard for those complete
  universes and never as a runtime slicing limit; unavailable, infeasible, skipped, and failed
  analyses remain explicit index rows with reasons;
- records candidate identity/scope, all equivalent knockout gene ids, loop eligibility, and
  expected/attempted/completed/failed/skipped coverage in the exported indexes, summary, and
  provenance so full-set validation is machine-checkable;
- renders every completed Figure 5 scan with enforced candidate-reaction flux on x
  (`target_flux`) and target-product flux on y (`response_flux`); amplification is a wild-type
  candidate→product scan; a single-reaction knockout uses pre-deletion wild-type
  reference↔zero when nonzero or the full feasible domain marked exploratory when already zero;
  growth remains the configured minimum-growth constraint plus secondary `biomass_flux` output
  rather than an axis;
- writes the fixed `01_preflight`–`07_validation` schema, resolved config, byte-for-byte source
  SBML, conditioned SBML, per-analysis metadata sidecars, replay/render/validate scripts,
  provenance, summary, and authoritative manifest;
- emits recommendations only under the declared evidence policy: validated beneficial
  single-gene deletion; a method-specific FSEOF or FVSEOF hypothesis with supporting response
  and a complete non-loop diagnostic; or RobustKnock positive guaranteed product. It does not
  require FSEOF/FVSEOF overlap and does not invent combined knockout-amplification
  recommendations.

The publication layer has a separate contract. `validate_production_run` validates source
artifacts, recomputes declared SHA-256/byte sizes, requires metadata sidecars for complete CSV
artifacts, and returns all structural issues without changing data. `nature-r` is the sole
publication backend: it invokes the checked-in R script with `--vanilla`, fails on missing R
or a required renderer package, reads only manifest-declared CSVs, and writes non-empty 300-DPI
PNG plus editable PDF/SVG figures. Package dependency metadata supplies compatible minimum
versions at runtime; the exact versions are supplied by `renv.lock`, asserted in CI, and
recorded in the figure manifest. HTML construction is deterministic; the standalone copy
embeds every figure.
Optional response/sampling panels must be declared as rendered, skipped, or failed, so absence
cannot be mistaken for negative evidence.

These tests establish deterministic composition and artifact integrity. They do not validate
a predicted knockout or amplification target experimentally, and the Nature
Genetics-inspired theme is not journal approval.

## Run provenance

Numerical result metadata includes:

- deterministic model SHA-256 over reaction bounds, objective coefficients, GPR rules, and
  stoichiometry;
- model id and active solver, **and the solver version**;
- the run timestamp in UTC (`timestamp_utc`, ISO-8601, never local time) and the run `seed`
  (`null` when the method has none — a statement that this run had no seed, not an omission;
  OptKnock/RobustKnock record the explicit strain-design seed, default `0`);
- the platform, machine and processor;
- Python, CMM, COBRApy, NumPy, pandas, and SciPy versions;
- method parameters, source-state provenance, and counts of non-optimal scan levels where
  applicable;
- for the production family, `applied_condition`: the medium as applied (with the components
  the model could not express, under `dropped`), the oxygen exchange bounds, and every carbon
  uptake with its rate — so a result file states its own conditions rather than leaving a reader
  to reconstruct them from a fingerprint.

Save result tables together with their metadata, the exact input model and omics files, the
Git commit, and the `uv.lock` used for the analysis. Model fingerprints detect scientific
model changes; they do not replace archival storage of the model file.

## Known limits

- Under the bundled restricted Gurobi license, L2 MOMA and `rmta_continuous` are exercised on
  small deterministic networks, including a direct QP solve for `rmta_continuous`. L2 MOMA on
  `e_coli_core` (286 variables, because COBRApy adds one distance variable per reaction)
  exceeds the restricted license and requires a full one.
- `fva` resolves `processes` to 1 below 500 analysed reactions, but leaves the pool to cobra's
  configuration at genome scale. A full *i*JO1366 FVA on macOS therefore still needs the caller
  to wrap it in `if __name__ == "__main__":`, or to pass `processes=1`. That residual is
  cobra's multiprocessing behaviour, not CMM's.
- No current automated test validates novel target predictions against new wet-lab data.
- MIQP/MILP genome-scale runtime and feasibility depend strongly on the solver license,
  tolerances, candidate set, and network compression.
- rMTA source-state preprocessing and FVSEOF's `linear_flux_couplings` must be described
  exactly; CMM does not silently claim the full preprocessing protocol of either paper.
- MTA/rMTA transformation scores floor the published denominator at the run's own `epsilon`, so
  a "perfect transformation" receives a large finite score rather than `+∞`. The scale is a
  documented function of a documented parameter. Check `n_distinct_scores`,
  `largest_tie_block` and `score_resolution` in the ranking's metadata before quoting a top-k:
  `TargetRanking.sorted` breaks ties on `target_id`, so a large tie block is an alphabetical
  slice, not a ranking.
- Floating-point optima can vary within solver tolerance. Manuscripts should report the
  solver, tolerance, model fingerprint, parameters, and acceptance tolerance.
- Sampler convergence is not asserted. OptGP needs large sample counts (roughly >1000) to mix;
  a small or heavily thinned run can under-represent the feasible space without failing. Report
  the sampler, sample count, thinning, and seed, and check stability across seeds before
  drawing conclusions from sampled means.
