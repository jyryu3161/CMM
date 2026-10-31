# Changelog

## 0.5.0 — unreleased

The flux map, which shipped in 0.4.0 but could not be reached from the application,
is now a working part of it. Numerical APIs and result semantics are additive in this release;
obsolete internal planning documents are removed from the public source tree.

### Added

- **A second canonical workflow: transformation-target discovery.** `cmm transformation-targets
  --config CONFIG` ranks the knockouts that move a source metabolic state toward a target one,
  from two gene-expression profiles, using the published MTA (Yizhak et al. 2013) or rMTA
  (Valcárcel et al. 2019). It composes existing CMM services and adds no numerical method: the
  MIQP, transformation score and Equation 9 are `revert_targets`, the MOMA baseline is
  `transformation_targets`, and the reference state is E-Flux2 or LAD. Six stages write the
  same schema-v2 run bundle SC-01 does, so one manifest format and one path-discovery surface
  serve both. Shipped with [`SC-02`](docs/scenarios/SC-02-transformation-target-discovery.md)
  and the `cmm-transformation-engineering` skill, whose interview confirms on every run which
  file is the source — nothing in the model can detect a swap, and the reversed run is a
  correct answer to a different question.
- **A report renderer and a completion gate for transformation runs.**
  `cmm.reporting.render_transformation_report` draws four panels — best-case against worst-case
  component score, score against rank, transformation rank against the MOMA baseline, and rank
  against epsilon — through
  `render_transformation_figures.R`, the same checked-in R/ggplot2 path SC-01 uses, and writes
  the same artifact pair: a linked `report.html` and a `report_standalone.html` carrying every
  figure as a data URI, with 300-DPI PNG plus editable SVG and PDF for each panel.
  `validate_transformation_run` is the completion gate: it checks each declared artifact
  against its recorded hash and size, that the ranking is ordered 1..N by descending score,
  that a skipped stage records why, and that the standalone page carries every image it
  references. A page that opens is not evidence a run finished, and each of those failures
  looks like success in a browser. `cmm report render` and `cmm report validate` read the
  run's own manifest to tell the two workflows apart, so neither has to be told which it is.
  The page states in its own body what a reader would otherwise have to dig out of the
  provenance: that the reference state is not the published iMAT one, that epsilon was chosen
  rather than derived, that the candidate
  count is the denominator of any percentile claim, and that source and target are an input
  rather than a finding.
- **`cmm.omics.gene_directions_from_replicates`** implements the Student's t-test Yizhak et al.
  specify for selecting changed genes. `gene_directions` cuts on fold change, which is all a
  single measurement per gene supports; the published route needs replicates and had no
  implementation. `restrict_to_top_changed` now also accepts `gene_p_values`, ranking the
  changed set by the strongest evidence rather than the largest fold change, and records which
  ordering was used.
- **`cmm.features.coupled_reaction_sets`** groups reactions whose deletion has the same
  consequence, from the null space of S. A ranking over knockout candidates is only meaningful
  if each candidate is a distinct intervention: three reactions of an unbranched pathway are
  one intervention, and counting them separately inflates the denominator of any "top *N*%"
  claim.
- **The Revert Metabolism tab now says what its parameters are for, and lets you choose the
  source-state estimator.** The tab drives a three-stage pipeline whose parameters mean nothing
  out of that context — ε and α belong to different stages and answer different questions — so
  the controls are grouped by stage, each group opens with a sentence naming that stage's job,
  and every field carries a tooltip giving its published value and what happens if it is wrong.
  The source flux state was hard-wired to E-Flux2; **LAD is now selectable**, and the tab states
  that neither is the iMAT-plus-sampling state Yizhak et al. use. The fold-change threshold that
  decides which genes count as changed was likewise fixed at ±1 with no way to see or set it.
- **MTA's two published preprocessing parameters are now reachable, including from the GUI.**
  The significant-flux-change threshold ε was fixed at its default on every graphical run —
  the Revert Metabolism tab never passed it — and Yizhak et al.'s cut that keeps only the most
  differentially expressed reactions in the changed set had no implementation at all. Both are
  now controls on the tab and arguments in the API: `differential_expression(top_n_changed=…)`
  and the new `restrict_to_top_changed()`. The cut is not cosmetic: each changed reaction adds
  one binary variable to the MIQP, so it decides whether a genome-scale run is tractable.
  Defaults are unchanged: no cut, and ε at `revert.DEFAULT_EPSILON`.
- **Source-checkout installation now uses a uv-managed Python 3.12 by default.** The shell and
  PowerShell installers no longer inherit an unsupported operating-system Python such as the
  Python 3.9.6 supplied by older macOS releases. They accept an explicit Python 3.10–3.12
  override, validate existing virtual environments before reuse, and fail with recovery
  guidance instead of silently replacing one.
- **Publication and agent-interface provenance is explicit.** The source distribution now
  includes the repository skill, installation scripts, Python and R lock files, third-party
  notices, and AI-use guidance. The independently written clarification interview records its
  public behavioral inspiration, and the complete Escher map license ships beside the map.
- **The Flux Map tab is always available, and needs no file to get started.** It previously
  appeared only when the window was constructed with `map_path=`, which `python -m cmm.app`
  never passes — so the feature existed but no ordinary launch could reach it. The tab now
  offers two layouts: a curated Escher map, and a dependency-free schematic of the highest
  `|flux|` reactions for models no map fits. The caption states which one is on screen.
- **Escher's *E. coli* core map is bundled** (`src/cmm/resources/`) and offered automatically
  to any model containing at least half its reactions — `e_coli_core` (94 of 95) and also
  genome-scale reconstructions such as iJO1366, since viewing a genome-scale model on a
  central-metabolism map is how Escher is used. Redistributed byte-for-byte under Escher's
  MIT license; provenance, SHA-256 and the citation to King et al. (2015) are in
  `src/cmm/resources/ATTRIBUTION.md`, and a test asserts the digest so the attribution cannot
  silently go stale.
- **The flux map draws expression-derived and knockout flux states, not just FBA/pFBA.** The
  Omics tab gained a condition selector, and both it and the Comparison tab gained a **Show on
  flux map** button. These are the two results a map is most useful for — a distribution
  over the whole network is not something to read as a 95-row table — and they were the two it
  could not show.
- **The Flux Map tab runs FBA or pFBA itself** (`draw: FBA | pFBA`), and has no separate
  render button: changing the layout or the reaction count redraws the loaded flux state, and
  those controls never solve. A render button invited pressing **FBA** to redraw, which
  silently replaced a drawn omics or knockout flux state with a fresh FBA solve. pFBA was reachable only by
  running it on the Simulation tab first, which nothing on the map suggested; the tab silently
  drew FBA and looked like that was all it could do.
- **The figure title names the method behind the numbers** (`e_coli_core — LAD · condB`). The
  map previously drew whichever of FBA or pFBA had run last while its title said only "flux
  map", so a pFBA map was indistinguishable from an FBA one.
- **A map drawing can be laid under the flux colouring** (**Background…** on the Flux Map tab).
  CMM reuses an Escher map's layout but not Escher's rendering, so a picture exported from the
  same map — SVG, PNG or JPEG — supplies the drawing while CMM supplies the flux. It is placed
  in the map's canvas coordinates, so an export of the loaded map needs no adjustment, and CMM
  drops its own labels and node markers while one is shown rather than doubling them. SVG is
  rasterised with Qt in `cmm.app`; `cmm.visualization` takes an array and stays free of Qt.
- **`File ▸ Open Escher Map…`** (and **Load map…** on the tab) loads any Escher map JSON. A map
  whose reactions are absent from the loaded model is refused with a message instead of being
  drawn as a blank grey network.
- `cmm.resources.bundled_map_for(model)` returns the bundled map's path when it suits a model,
  and `None` when nothing does.
- **The Revert Metabolism tab takes replicates and runs the published t-test.** Its expression
  files were limited to one value per gene, so the only direction test the GUI could offer was
  a fold-change cut — Yizhak et al.'s (2013) step 2 needed the Python API. Both files now read
  a gene id column followed by one column per replicate, and with at least two on each side
  the tab runs `gene_directions_from_replicates` (Student's t-test at a settable P cutoff)
  through the same `reaction_directions` / `restrict_to_top_changed` path the SC-02 workflow
  uses. With fewer, the t-test entry is disabled and a note states that the published test was
  not applied. Values stay linear expression, as everywhere else in the app: they are
  log2(x + 1) transformed per replicate for the comparison, and the E-Flux2/LAD source state
  gets the replicate mean as measured. The summary line now records the test, the replicate
  counts, how many genes were called changed and how many reactions were labelled. A
  one-column file scores exactly as before.
- `cmm.omics.gene_directions_by_fold_change` — the fold-change counterpart of
  `gene_directions_from_replicates`, returning the same evidence frame so a caller can swap
  tests without also changing what the numbers mean. The SC-02 workflow's fold-change branch
  is now this function.

### Fixed

- **OptKnock and RobustKnock can no longer run forever.** `straindesign` 1.18 finds the optimal
  design, then solves a second MILP to minimise the interventions. When that solve is not
  optimal - on a 1,440-reaction model Gurobi reported "Model is infeasible" straight after
  "max constraint violation 4.75e-09 exceeds tolerance" - the inner loop breaks without
  excluding the solution it just found, and the outer loop re-solves the identical MILP.
  `max_solutions` is never reached because nothing is added, and CMM passed no time limit, so
  two production runs sat in it for 24 hours with no output and no way to see where. Whether
  a run falls in depends on solver numerics (it moved with the thread count), not on the
  model alone. `optknock`/`robustknock` now take `max_seconds` and the workflow takes
  `strain_design_max_seconds`, both defaulting to 3600 s, forwarded as `straindesign`'s
  `time_limit`. The regression test drives the real `straindesign` loop on a five-reaction
  model with the second solve forced non-optimal.
- **A search the clock ended is no longer reported as complete.** `straindesign` tests the
  remaining time only before a solve and keeps the status of the last one, so it can return
  `optimal` with no designs, or fewer than a finished search would have found. CMM took that
  at its word (`straindesign_search_complete=True`). It now also compares the elapsed time
  with the limit and records `straindesign_time_limit_reached`,
  `straindesign_search_elapsed_s` and `max_seconds`; a cut search emits a `UserWarning` and
  the workflow reports its artifact as `partial`, as it already did for `time_limit`. The
  elapsed time includes preprocessing, so a search that finished just inside its limit can be
  counted as cut - the safe error.
- **The desktop forms no longer lay themselves out differently on macOS.** `QFormLayout` takes
  its alignment and field-growth policy from the platform style, and the styles disagree: most
  give `AlignLeft` with `AllNonFixedFieldsGrow`, macOS gives `AlignHCenter` with
  `FieldsStayAtSizeHint`. Under the macOS pair every group box centred its own label-and-field
  block at that block's natural width, so a column of boxes landed at as many different x
  positions as it had boxes and no two inputs shared a width — while the same code rendered
  correctly under every other style, CI included. Both policies are now set explicitly. Three
  smaller defects went with it: the stage boxes' label columns are pinned to one width
  measured from the resolved font on first show (a minimum width is only a floor, so a label
  that rendered wider still dragged its own box out of line); combo boxes and spin boxes are
  styled to one 26px height, so a form row's pitch no longer depends on which control it
  holds; and the inputs stop growing at a readable width instead of running the full width of
  the panel on a large display.
- **MTA/rMTA candidates are once again scored against a common yardstick.** The
  impossible-change mask was evaluated inside each candidate's knockout context, so it saw
  that candidate's modified bounds. A reaction could be masked for one knockout and not for
  another, which meant the steady set — the denominator of the transformation score — differed
  between candidates whose scores were then ranked against each other. The mask is now applied
  once against the unperturbed model, as in Yizhak et al.'s preprocessing, and the resulting
  count is reported as a single `n_impossible_masked` in the ranking's provenance rather than
  being a per-candidate quantity that never surfaced. Rankings on models where the mask fires
  will change; the toy-network test suite is unaffected. Preparing the direction maps once
  also removes two map constructions per candidate.
- The source-tree fallback for `cmm.__version__` now matches the `0.5.0` package metadata.
- **The flux map no longer drifts right, clip its title, or draw its colorbar as a hairline**
  when the GUI stretches the figure to a wide panel. `colorbar` re-anchors its parent axes to
  the right, so every bit of the shrinkage `set_aspect("equal")` applies was taken off the left
  edge; `tight_layout` solved the margins once at the authored size and they no longer fitted
  once stretched. The figure now anchors centre and re-solves its layout on every draw.
- **The schematic layout folds long chains into rows instead of collapsing them.** At 25
  reactions the carbon backbone is one ~30-node chain, and laying it along a single row put
  nodes at *zero* separation — markers and labels merged into a smear. Rows now use one fixed
  node spacing, so a two-node branch no longer stretches across the whole panel either.

### Changed

- **OptKnock and RobustKnock now have a one-hour budget by default** (`max_seconds=3600.0`,
  workflow `strain_design_max_seconds`). A search that legitimately needs longer returns its
  designs so far, marked incomplete, instead of running to the end; pass `None`/`null` or a
  larger value to restore the unbounded behaviour. Where it terminated, a 1,440-reaction model
  took 4-5 minutes per search.
- Internal planning ledgers, completion notes, and machine-specific analysis notes are no
  longer part of the public source tree or source distribution. Generated SC-01 run bundles
  remain user-owned analysis artifacts rather than release inputs.
- CMM's public expansion is now **Constraint-based Metabolic Modeling** across package metadata,
  the README, citation metadata, documentation, and desktop title/header.
- **The schematic carries a colorbar** instead of an `∝ |flux|` formula in the margin, drawn
  from the same truncated colormap span the arrows are coloured from.
- The schematic's reaction count is capped at 25. Beyond that the cross-row arrows dominate and
  the figure stops being readable — that is a curated map's job, not a fallback schematic's.
- `test_data/` is removed. Its last remaining file was an unattributed copy of the same Escher
  map; the scenario harness now uses the bundled, attributed one.

## 0.4.0 — unreleased · **BREAKING**

A fidelity release. Five independent audits compared every implemented method against its
original publication; this release acts on their findings. It changes reported numbers, removes
one public parameter and one result field, and corrects documentation that overstated what the
code does. **Results produced with 0.3.0 are not directly comparable and must be regenerated.**

The scenario *documents* have been corrected. Any result produced by 0.3.0 — a saved run
directory, an exported CSV, a figure — is stale and must be regenerated: the changes below move
yields, design rankings, expression-derived fluxes and one result field's meaning.

### Breaking

- **`aerobic=True|False` is removed** from `theoretical_yield`, `production_envelope`, `fseof`
  and `fvseof`, in favour of `condition=`. There is no deprecation period. It was a second,
  redundant way of saying what the medium already says, and it was actively dangerous:
  `optknock` and `robustknock` never accepted it, so a caller who set `aerobic=False` got an
  aerobic design with no warning. `condition=` is now the one convention across
  `fba`/`pfba`/`fva`, the production family, `flux_response`, the samplers and — new here —
  `optknock`/`robustknock`, which previously took no context at all and depended silently on
  the caller's model state.
- **`theoretical_yield` returns different numbers.** Media presets now close CO₂ *uptake*
  (secretion stays free), because an open CO₂ uptake inflated yields with no guard firing: the
  only check tests the carbon ceiling, and a CO₂-inflated yield below the ceiling passes it. On
  anaerobic `e_coli_core`/`EX_succ_e` the reported molar yield falls from **1.3906 to 1.2000**
  against a ceiling of 1.5 — the old figure was 15.9% high, obtained by taking up 6.95 mmol of
  CO₂ that a closed anaerobic fermentation does not supply. Closing it costs 0.000% of aerobic
  growth and ≤0.25% anaerobic, and removes 8 of 8 carbon-ceiling violations. The error also
  propagated to every FSEOF/FVSEOF scan level and to the production envelope (~13.7% high).
  `theoretical_yield` additionally raises on a non-boundary reaction, computes the carbon
  ceiling from every carbon-containing uptake flux rather than one nominated substrate, and
  multiplies element counts by their stoichiometric coefficients.
- **The flux-response `bottleneck` field is removed and replaced by a shadow-price measure.**
  `ResponseBottleneck` and `.bottleneck` are gone; `FluxResponseResult` gains `.phases`
  (intervals of constant shadow price), `.limit` (the phase boundary past which the response
  falls faster than a stated threshold), `.shadow_price_at()` and `.phases_frame()`. The old
  field reported the steepest decline of a finite-difference gradient. The response curve is an
  LP optimal-value function parameterized in one bound and is therefore concave piecewise
  linear, so that argmin locates the edge of the scan grid, not a property of the network: the
  reported location moved by up to 29.53 flux units and the `found` flag inverted on `PGI`,
  `TPI` and `EX_o2_e` as `n_steps` went from 6 to 160, where the replacement moves by
  ≤2.6 × 10⁻¹³. No published criterion defines a bottleneck as the argmin of a
  finite-difference slope. The word "bottleneck" is dropped from the shipped surface.
- `ComparisonResult` splits its overloaded `distance` field. The raw solver objective is kept
  as `objective_value`; `distance` is Euclidean for `moma_l2` (the square root of the QP
  objective — the exported value was `Σd²`, e.g. 1303.99 where the true distance is 36.11),
  L1 for `moma_l1`, and `None` for ROOM, whose value was never a distance but a count of
  switched reactions. `distance_kind` records which. Frame columns are renamed accordingly.
- `fva` returns a result object carrying `run_provenance` instead of a bare
  `dict[str, FluxRange]`, forwards `loopless` so loopless FVA is reachable, and exposes
  `processes`.
- Media presets are renamed so none implies M9 or another standard, and every component is
  documented. A missing growth-limiting component now raises instead of silently producing a
  different experiment; missing minerals warn and are recorded in provenance under `dropped`.
  `Medium.apply_to` (and therefore `apply_medium`) returns a `MediumApplication` rather than a
  plain dict; it still behaves as the mapping of applied uptakes, and adds `dropped` and
  `to_provenance()`.
- `fvseof`'s `group_constraints` keyword is renamed **`linear_flux_couplings`**, and its
  `n_group_constraints` metadata key to `n_linear_flux_couplings`. Park et al.'s
  grouping-reaction constraints are STRING-derived on/off pairs with a normalised-flux
  inequality; CMM takes caller-supplied linear equalities, which is a different object and now
  says so in its name.
- `batch_comparison(method="room")` defaults to the **lethality** tolerance pair (δ=0.1,
  ε=0.01) because a screen asks which deletions are lethal, while `room()` and
  `knockout_comparison()` keep the flux-prediction pair. A ROOM screen re-run therefore reports
  about 24% fewer switches than an archived pre-0.4.0 one; pass
  `room_use_case="flux_prediction"` to reproduce the old numbers.
- `lad`'s `weight_threshold` default changes from **0.01 to 0.0**, so a low-expression reaction
  is driven toward zero flux as Lee et al. intend instead of being dropped from the objective,
  and the test is now `weight >= weight_threshold`. `lad` also gains `reaction_sigma` for Lee
  et al.'s per-reaction `1/sigma` weights; the unweighted default is a documented deviation.
- `gene_perturbations`, `grouped_gene_perturbations` and `reaction_perturbations` return a
  `PerturbationList` (a `list` subclass) carrying `inert_dropped`, `n_inert_dropped` and
  `provenance()`. Iteration, `len()`, indexing and slicing are unchanged, so no caller breaks.
- The scenario/report layer renames its knockout-screen column. What was
  `moma_l2_distance_to_reference` held the QP objective `Σd²`; it is now written as
  `moma_l2_solver_objective`, and `moma_l2_distance_to_reference` holds the Euclidean distance
  — a factor of about 36 apart. `distance_kind` and `n_changed_reactions` are written
  alongside, and both SC-01 and SC-03 reports state which quantity each column holds.

- `BatchComparisonResult.to_frame()` omits the column the screen's method cannot fill: a MOMA
  screen has no switch count, and a ROOM screen has no distance since ROOM reports a count
  rather than a norm. Both were previously written as a column of `NaN`, which in an exported
  CSV is indistinguishable from a run that failed. `distance_kind` survives on a ROOM screen —
  `"none"` is a statement about ROOM, not a missing value — and
  `to_frame(drop_empty_method_columns=False)` returns the full schema for a caller
  concatenating screens run under different methods.

### Scientific correctness

- **Ranked targets carry the model's own name for the thing they name.** `TargetScore` gains
  `target_name`, populated from the SBML the run already stores — `model.genes[…].name` for a
  gene-level target, `model.reactions[…].name` for a reaction-level one, never from an imported
  expression table, which is untrusted data. Report tables and figure labels read `6510_AT1` or
  `b1602` without it, which no reader can interpret. The name sits *beside* `target_id` and never
  replaces it: names are neither unique nor always present, whereas the id is what provenance and
  reproduction are written against. The export column appears only when something was actually
  named, so a ranking over a model that names nothing exports exactly the table it did before.
  Shared by FSEOF, OptKnock and revert, so both workflows gain it at once.
- **The transformation report leads with the component scores rather than the shape of the list.**
  The old first panel plotted transformation score against rank; both axes are artifacts of the
  ranking, and under rMTA a single amplified candidate can be an order of magnitude above the
  rest, flattening the other 400 onto the axis. The new `fig01_component_scores` plots best-case
  (bTS) against worst-case (wTS), the space Equation 9 actually branches on, marking the
  candidates the robust score amplifies and labelling them. Score against rank is kept as
  `fig02_transformation_ranking` and is now the required panel, because it is the only one every
  method can draw: under MTA the three component scores are one number, so the component panel
  records why it cannot be drawn instead of failing the run. The MOMA and epsilon panels shift to
  `fig03` and `fig04`.
- **The coupling disclosure no longer claims the source papers reduce candidates by partial
  coupling.** Yizhak et al. (2013) apply it only in their reaction-level validation analyses
  — *"in the validation analyses, the set of simulated knockouts is composed of a member from
  each partially coupled set"* — and Valcárcel et al. (2019) rank gene knockouts with no
  coupling reduction at all. `cmm.features.coupling`, the transformation workflow docstring,
  the `cmm-transformation-engineering` disclosure list and the rendered report nonetheless
  presented full-versus-partial coupling as one of three structural departures from "the
  papers", which overstates what either source does and reads as a deviation from rMTA that
  does not exist. Each now states what a reaction-level CMM run computes, names the Yizhak
  scope, and says rMTA applies no such reduction. A gene-level run — the default, and the
  perturbation level of both published human analyses — never reaches the coupling path. No
  numerical behaviour changes.
- **GPR `OR` resolution is per method, matching each source paper.** CMM's global `max` matched
  none of them. `gene_to_reaction_weights` gains an `or_rule` parameter defaulting to `"sum"`,
  which is what both Kim et al. 2016 (E-Flux2) and Lee et al. 2012 (LAD) specify; `AND = min`
  is unchanged and was already correct. MTA/rMTA direction sets are rewritten separately using
  Yizhak et al.'s ternary rule — all subunits changed (AND), at least one changed (OR), **mixed
  ⇒ unchanged** — which the previous signed `min`/`max` over `{−1, 0, +1}` got wrong in 4 of 7
  label cases with a directional bias. The chosen rule is recorded in provenance. This changes
  every E-Flux2, LAD, MTA and rMTA result.
- LAD now fits **absolute** flux as Lee et al. specify, via a forward/backward split
  (`v = f − b`, residual against `f + b`), instead of the signed residual that systematically
  penalized reverse flux: a reversible reaction at `v = −5` with target 5 scored deviation 0 in
  the source and 10 in CMM. The problem remains an LP.
- OptKnock and RobustKnock apply the existing `_actionable_reaction` filter to their candidate
  set by default, with an opt-out. Designs deleting boundary reactions with no GPR
  (`EX_co2_e`, `EX_ac_e`, `EX_for_e`, `EX_etoh_e`, `EX_lac__D_e`) are not realisable as gene
  deletions and are no longer proposed.
- The same filter now also rejects reactions whose only gene is COBRA's `s0001`
  spontaneous-reaction placeholder, which is not a gene anyone can delete. On `e_coli_core`
  that removes `ACALDt`, `CO2t` and `O2t`, taking the candidate set from 69 to 66 — and with it
  every design the filter had been letting through under a gene it could not name. The
  top-ranked anaerobic succinate design changes accordingly, and both arms of the production
  scenario now return a design that is buildable as written.
- FVSEOF selects on Park et al.'s **joint** sign of ΔV_avg and Δl_sol rather than V_avg alone;
  the capacity slope was already computed and simply never used for selection. Its default
  `n_steps` rises from 8 to 10, Park's stated minimum.
- ROOM exposes both published δ/ε presets and selects by use case — 0.03/0.001 for flux
  prediction, 0.1/0.01 for lethality — instead of using the flux-prediction pair everywhere.
  Measured effect on ranking: 531 vs 401 switches over 35 genes, a 24% shift.
- MOMA's reference protocol is pFBA everywhere. The Python API and the GUI previously
  disagreed (pFBA vs FBA), which was the actual defect; the numerical difference is negligible
  (≤0.00096 h⁻¹ on iJO1366, 0 of 24 knockouts reclassified, byte-identical on `e_coli_core`)
  while FBA's alternate optima are not — 8 identical FBA solves returned 3 distinct vectors.
- `feasible_range()` is computed by FVA on the scanned reaction rather than from the scan grid.
  At the documented default `n_steps=20` it returned (−37.3684, 6.8421) against a true range of
  (−38.0997, 9.9463), understating `PGI`'s headroom by 31%. The bias appeared only when a
  growth floor was applied.
- `run_provenance` records the run timestamp, the seed, the solver version and the platform,
  none of which it emitted before — a provenance record without the seed and solver version
  does not let anyone reproduce a sampling or MILP result. It also records the applied
  condition in full: medium name, oxygen exchange bounds, substrate and uptake rate, and the
  medium components applied and dropped.
- **The whole MOMA/ROOM family now carries that block, and it did not before.** `moma`,
  `room`, `knockout_comparison` and `batch_comparison` returned results whose only metadata
  was `reference` / `reference_provenance` / `reference_method` (plus ROOM's tolerances) —
  no `timestamp_utc`, no `seed`, no `solver`, no `solver_version`, no `platform` and **no
  `model_sha256`** — so the perturbation-response engine, which is what a knockout screen is
  made of, failed `AGENTS.md` rule 4. All four now carry the full `run_provenance` block
  beside the reference keys they already had. `seed` is `null`: the methods are
  deterministic and no seed is invented. For `moma`, `room` and `knockout_comparison` the
  fingerprint is of the model **as handed to the solver**, so a `knockout_comparison` record
  fingerprints the knocked-out model, not the wild type, and `parameters["knockouts"]` names
  the reactions forced to zero.
- **`batch_comparison` returns a `BatchComparisonResult`**, a `list` subclass with the same
  iteration, `len()`, indexing and `sorted()` behaviour as the plain list of rows it replaces,
  adding `.to_frame()` and `.metadata`. The screen's provenance lives on the container rather
  than the row because every row of one screen shares one model, reference, method, tolerance
  pair, solver and machine; per-row duplication would copy a 16-key block across the 1,367
  gene knockouts of a genome-scale screen and compute 1,367 model fingerprints to do it.
  `metadata["model_sha256"]` is the screened model before any knockout — the one model common
  to every row. The container's provenance also records what the enumeration left out
  (`n_perturbations`, `n_inert_dropped`, `n_candidates_considered`): `gene_perturbations`
  drops genes whose deletion blocks no reaction, 66 of 137 on `e_coli_core`, and the screen no
  longer understates its coverage silently.
- `predict_condition_fluxes` returns a `ConditionFluxes` carrying a `metadata` block for the
  multi-condition job — the model fingerprint every condition was solved on, the integration
  method, the condition names and how many solves were non-optimal. Each condition's own
  `OmicsFluxResult` keeps its own block for when its numbers are lifted out alone.
- `tests/test_provenance_surface.py` parametrises over every public service that returns
  numbers and fails if any of them drops the block, so the gap cannot be reintroduced
  silently. 23 services, 23 carrying it.
- `transformation_targets` records provenance on its `moma` path, which had none, and records
  the target-state identity on both paths so an `mta` run is distinguishable from
  `revert_targets`. Its default stays `method="moma"` and is now labelled as such, with the
  paper's own "markedly inferior" verdict quoted in the docstring.
- **MTA/rMTA transformation scores no longer return `±inf`.** The published denominator (the
  steady-set L1 deviation) is floored at the run's own `epsilon`: a steady-set deviation
  smaller than the change the method itself calls significant cannot be resolved from zero.
  On the SC-02 condition pair this removes a 38-gene `+∞` tie block (distinct scores 11 → 14,
  largest tie block 38 → 18), so an `mta` top-k is a ranking rather than an alphabetical slice
  of a tie. Rankings also now carry `n_distinct_scores`, `largest_tie_block` and
  `score_resolution` so a reader can tell the difference. Every MTA/rMTA artefact must be
  regenerated — for `rmta` the effect of the GPR rule alone is ρ = 0.361 against the old
  ranking.
- `rmta_continuous` results carry a `continuous_heuristic_not_rmta` column in `to_frame()`, so
  an exported CSV cannot be mistaken for a published rMTA result. `rmta`/`mta` runs do not get
  the column.
- E-Flux2's L1 fallback records `eflux2_l1_fallback` in its archived provenance instead of
  mislabelling itself as the QP path.

### Documentation

- **`docs/scenarios/SC-01` presented an aerobic result as the anaerobic answer.** The design
  `{CO2t, FORti, PGI}` at `guaranteed_product` 10.4063 **does not grow anaerobically at all**;
  it was obtained under an aerobic medium while `aerobic=False` was passed to the functions
  that had the parameter. The correct anaerobic design is `{ACALD, D_LACt2, THD2}`, guaranteed
  succinate **9.910758** at growth **0.090648**, growth-coupled, with an enforced-flux scan on
  the deleted model bounded below at **4.794286**. The scenario's condition-setting procedure
  is rewritten so the medium and aeration are set once, before the first solve, and every later
  step visibly inherits them.
- Provenance labels that overstated what is tested are corrected: FSEOF's `criterion` string no
  longer claims a monotonicity test the code does not perform; FVSEOF's `robust_targets()` flag
  and CMM's FSEOF selection rule are labelled as CMM's own rather than Choi et al.'s or Park
  et al.'s; FVSEOF's grouping constraints are renamed to reflect that they are caller-supplied
  linear flux couplings, not Park et al.'s STRING-derived reaction pairs;
  `production_envelope`'s docstring no longer calls its output a phenotypic phase plane; and
  `max_solutions` is documented as a cap on MILP solutions, not on distinct designs.
- Citations added: **Lee D, Smallbone K, Dunn WB et al. (2012)** *BMC Syst Biol* 6:73 for LAD;
  **Lee KH, Park JH, Kim TY, Kim HU, Lee SY (2007)** *Mol Syst Biol* 3:149 for `flux_response`,
  with **Edwards & Palsson (2000)** for its biomass (robustness) case; **Schneider et al.
  (2022)** for the `straindesign` package, which carries no citation of its own and must be
  cited alongside Burgard et al. (2003) and Tepper & Shlomi (2010).
- Citations corrected: Kim et al. (2016)'s authors are **Kim MK, Lane A, Kelley JJ, Lun DS**;
  Edwards et al. (2002)'s title is "Characterizing the metabolic phenotype: a phenotype phase
  plane analysis"; **Orth, Thiele & Palsson (2010) contains no definition of shadow price or
  reduced cost** and must not be cited for either. `transformation_targets` is documented as
  *not* a CMM invention — both paths map to published Yizhak et al. (2013) methods — and
  `flux_log_change` is marked explicitly as a CMM utility with no published source that must
  not be cited to any paper. L2 MOMA is attributed to Segrè et al. (2002); the L1 variant,
  which that paper does not contain, is described as COBRApy's linear variant with no citation
  until a quotable source is obtained.
- The README's Gurobi statement is corrected: the restricted license caps quadratic problems at
  **200 variables** (not quadratic terms), a practical ceiling of roughly **100 reactions**, and
  **L2 MOMA on `e_coli_core` (286 variables) already fails** under it — neither genome-scale nor
  mixed-integer. `docs/VALIDATION.md` records the measured margin: the largest QP is E-Flux2 on
  `e_coli_core` at **190 variables against the 200-variable cap**. `rmta_continuous` has no test
  that solves its QP, so its README QP row is noted as unverified.
- `docs/architecture.md` no longer lists flux sampling and flux-response analysis as roadmap
  items; both ship, with method contracts, GUI tabs and figures. Dynamic FBA and
  enzyme-constrained modeling remain roadmap.
- Internal development notes no longer imply that unavailable private expression/model fixtures
  are part of the public repository. No substitute data has been invented.
- `CITATION.cff` carries an explicit `TODO` for the manuscript's final authors, ORCIDs and
  release DOI. The organization-only author entry is not acceptable in the published record and
  no names have been invented in its place.

### Removed

- `test_data/LAD.py`, `test_data/Simulator.py`, `test_data/saved_map.png` and
  `test_data/saved_map.svg`. They carried no author, licence or provenance header, their origin
  is unknown, and `docs/clean-room-policy.md` forbids both copied source files and "copied
  screenshots, icons, SVG maps, or bundled examples". Verified unreferenced before removal: the
  two Python files import only each other and nothing in `src/` or `tests/` imports either; the
  two map files have zero references anywhere in the repository; and `MANIFEST.in` ships only
  `test_data/*.json`, so none was in the distribution. `test_data/e_coli_core.Core
  metabolism.json` is kept — it is the Escher map used by `src/cmm/app/succinate_scenario.py`
  and documented in `docs/TUTORIAL.md`. The files remain in git history; the LAD attribution
  they supported is independently confirmed by Lee et al. (2012) and by Machado & Herrgård
  (2014).

### Added

- `flux_response`: scan an enforced flux through one reaction and maximize a response
  reaction at each point. Omitting the response gives the robustness reading (growth vs the
  target); naming a product exchange gives the production reading, with biomass recorded
  throughout and an optional `biomass_fraction` growth floor so the curve describes a viable
  strain rather than a non-growing ceiling. Reports the feasible window (by FVA on the scanned
  reaction), the optimum, the phases of constant shadow price and the response limit;
  infeasible scan points are returned with their solver status instead of raising.
- `random_flux_sampling` and `reference_constrained_sampling`: seeded, single-process-by-
  default flux sampling (OptGP/ACHR) over the feasible space, or narrowed to a window around
  a reference flux state. Results carry per-reaction statistics, a correlation matrix over
  varying reactions, and a bridge to `FluxState` for use as a MOMA/ROOM/MTA reference.
- Two GUI tabs, Flux Response and Sampling, each with a plot beside its result table, CSV
  export, and 300 DPI figure export. The Sampling tab additionally exports the raw ensemble
  (one row per sample, one column per reaction) separately from the per-reaction summary
  table, and clears it on model reload so a stale export cannot be attributed to a new model.
  The Flux Response tab's scan range is always shown and editable, filled with the target's
  detected feasible interval on every target change and resettable with "Detect range", so a
  scan never starts from a placeholder or another reaction's numbers.
- `flux_response_figure` and `sampling_figure` publication figures: the response curve with
  its infeasible range, wild-type marker, optimum, phase-boundary dividers, response-limit
  marker, and a secondary growth axis for product responses; and per-reaction sampled-flux
  violins against a reference solution.
- Simulation: FBA and pFBA fluxes are shown in separate columns (Reaction / FBA flux /
  pFBA flux / FVA range) so running pFBA no longer overwrites the FBA result; pFBA's minimal
  total flux is shown directly under the objective value.
- Production: a result table beside each plot — FSEOF/FVSEOF amplify/knockdown targets with
  their low/high enforced-level fluxes, and the production-envelope growth range per product
  flux — with a "Show all reactions" toggle that also lists unchanged reactions.
- Comparison (single run): a "Significant change ≥ X % of reference" threshold (default 3%)
  replacing the fixed 1e-6 cutoff, so alternate-optimum drift (notably ROOM) no longer reads
  as a knockout response; the solve is cached and re-filters without re-solving.
- Comparison (batch): the table reports wild-type and post-knockout biomass, an essentiality
  flag, and — when a target product is selected — that product's wild-type and post-knockout
  flux columns.
- Omics: multi-condition expression tables are supported in one tab, computing one predicted-
  flux column per selected condition, with a "Show all reactions" toggle.
- `AGENTS.md` and `CLAUDE.md`: agent operating instructions covering the scenario router, the
  goal→function table, the solver gate, the rules that keep results reportable, the run
  contract, and when to stop and ask.
- `docs/agent-reference.md`: signatures and result objects for every shipped service.
- `docs/scenarios/`: step-by-step metabolic-engineering pipelines (`SC-01`–`SC-03`) with a
  shared preflight and reporting contract, each step stating its preconditions, call,
  artifacts, decision rule, and failure handling.

### Changed

- **The desktop app is brought in line with the library.** The Comparison tab's reference
  template defaults to `pfba`, matching `reference_flux`'s own default — the two entry points
  disagreed, which was the actual defect behind the reference-state decision. The tab gains a
  ROOM tolerance-pair selector and reports the pair a count was produced under. The Flux
  Response tab presents phases, the response limit and the shadow price at the wild-type flux,
  with a phase table beside the scan table, in place of the removed bottleneck sentence. The
  Production tab passes a `Condition` instead of the removed `aerobic=` flag, reports the
  fraction of product carbon supplied by CO₂ uptake and any residual carbon imbalance, shows
  Park's nine-type index per FVSEOF row, labels the `robust` flag as CMM's own, and no longer
  forces `n_steps=8` under FVSEOF's Park-minimum default of 10. The Strain Design tab gains its
  own aeration selector and passes it to `optknock`/`robustknock` as a `Condition`. Applying a
  medium reports the preset's full display name and everything the loaded model could not
  express.
- **A GUI crash on the ROOM path is fixed.** The Comparison summary formatted
  `result.distance` with `:.4g` unconditionally; after the `ComparisonResult` split that value
  is `None` for ROOM, so selecting ROOM raised `TypeError: unsupported format string passed to
  NoneType.__format__`. The summary now names the quantity each method actually reports, and
  the offscreen smoke tests cover ROOM as well as MOMA.
- `src/cmm/app/genome_scale_scenario.py` no longer closes `EX_co2_e` by hand; the CO₂ policy
  lives in `cmm.core.media` and the scenario applies the preset instead. Verified equivalent on
  the textbook model: zero differing exchange bounds and identical growth, 0.8739215069684303.
- **Package surface.** `cmm` exports `FvaResult`, `MediumApplication`, `ResponsePhase` and
  `ResponseLimit` alongside `pfba` and `Medium`; `cmm.features` exports `ROOM_TOLERANCES`,
  `PerturbationList`, `perturbation_provenance`, `tie_structure`, `BatchComparisonResult`,
  `CarbonUptake`,
  `ProductionYield`, `ProductionEnvelope`, `FseofResult` and `FvseofResult`; `cmm.omics`
  exports `EFLUX2_DEVIATIONS` and `LAD_DEVIATIONS`. `INCLUDED_FEATURES` and `PLANNED_FEATURES`
  are unchanged: nothing new ships in 0.4.0, and `flux_response_analysis` is still the name of
  a shipped feature whose internals changed.
- Tabs are reordered so each analysis sits next to the one asking a similar question:
  Simulation and Sampling (what the model can do with no intervention), Comparison and Flux
  Response (the consequence of one — discrete and continuous respectively), then Production and
  Strain Design (proposing interventions), then the omics-driven tabs. Names and behavior are
  unchanged.
- Disabled combo boxes, spin boxes, and their labels are now visibly inert (muted fill and
  text). Previously a disabled control was almost indistinguishable from an active one, so the
  Sampling tab's reference options read as editable in `uniform` mode, where they do nothing.
- Figure swapping on the plotting tabs goes through one shared helper, so Production, Flux
  Response, and Sampling no longer each repeat the canvas/toolbar teardown.
- The Sampling tab can sample a deletion strain: an optional knockout picker (the Comparison
  tab's two-panel selector, now shared by both tabs) applies reaction or GPR-resolved gene
  deletions as a scoped condition, leaving the loaded model untouched. In "around a reference"
  mode the reference is built under the same knockouts, since a wild-type reference would put
  every deleted reaction outside its own sampling window. A knockout set that leaves no
  feasible space is reported as probably lethal, and the previous ensemble is dropped so a
  stale result cannot be exported under the new settings.
- Flux Response and Sampling figures re-run their layout on every draw, so axis labels and
  titles stay clear when the window resizes rather than only fitting the size they were
  authored at.
- `SC-02` is folded into `SC-01`. Finding knockout targets is an *inverse* problem, so the
  design step now uses OptKnock/RobustKnock rather than a MOMA/ROOM single-deletion screen, and
  MOMA/ROOM moves to the verification step where it belongs — predicting the built strain's
  immediate phenotype. A gap between the MOMA prediction and the design's guaranteed product is
  an estimate of the adaptive evolution required, not a refutation of the coupling proof.
  Without MILP the design step falls back to a single-deletion screen and the report must state
  that only single deletions were examined and coupling was not established. The reference-state
  step is correspondingly narrower: OptKnock does not consult a reference, so multiple reference
  states no longer strengthen the design, only the interpretation.
- The `cmm-guide` project-local skill is replaced by `AGENTS.md` + `docs/agent-reference.md`.
  Skills are only read by Claude Code; plain repository documents are read by any agent CLI,
  and the split keeps decision-critical material in the always-loaded entry point while the
  function reference loads on demand.
- Comparison: a two-panel knockout picker (searchable catalogue on the left, chosen knockout
  set on the right) replaces Ctrl/Shift-click selection, making the selection visible and
  clearable.
- Comparison: LAD/E-Flux2 reference templates are disabled (greyed out, not selectable) until
  their integration method has actually been computed on the Omics tab.
- Omics: loading an expression file only stores it and shows its filename; a separate Compute
  button runs the selected method, so the method can change without reloading.
- Revert / Transform: the loaded source/target expression filename is shown next to each input.
- Production: removed the redundant Run FBA button (duplicated the Simulation tab and produced
  no Production-tab output).
- GUI: the main window opens at a narrower default (1160×760). The wide single-row control
  bars on the Production, Comparison, Strain Design, and Omics tabs were split across rows so
  the content minimum width no longer forces the window far wider (~1574 → ~1146 px).
- GUI: combo popup entries are left-aligned.

### Fixed

- The `ruff check` and `ruff format --check` release gates pass again. `main_window.py` had a
  module-level assignment above its imports, which made every import in the file an `E402`
  error (17 in total), and two files were unformatted.
- GUI: combo-box and spin-box arrows now render (drawn from bundled SVG assets); the previous
  CSS border-triangle never drew in Qt and showed a grey box.
- GUI: tab labels no longer clip or overflow into a scroll button — the stylesheet font-weight
  that Qt's tab sizing ignored was removed and tab padding trimmed so all tabs fit.

- Comparison batch (MOMA/ROOM) no longer aborts the whole run when a lethal knockout makes the
  model infeasible; such a knockout is recorded as infeasible and the run continues.
- Revert / MTA (`_mta_miqp`, `_mta_qp`) no longer crash under Gurobi when a lethal knockout is
  infeasible (backend "Unable to retrieve attribute 'X'"); feasibility is probed before reading
  primals, so infeasible knockouts are skipped and the ranking completes.

### Known deviations, now disclosed rather than implied

- pFBA minimises total flux over **all** reactions, following COBRApy, rather than only
  gene-associated reactions as in Lewis et al. (11 reactions differ on iJO1366, none on
  `e_coli_core`).
- MOMA, ROOM and both samplers are pure delegation to COBRApy, so their numerical behaviour is
  the ecosystem's — including ROOM's constraint capping the perturbed objective at the reference
  objective, which is COBRApy's addition and not in Shlomi et al.'s Eqs 1–3.
- MOMA-L1 and ROOM return different flux vectors from COBRApy (12/95 and 54/95 reactions on
  `e_coli_core`); these are proven alternate optima, identical in objective to 4e-14. ROOM's
  predicted growth on iJO1366 differs by 1.235e-4 at an identical objective, so ROOM growth on
  a genome-scale model must not be quoted to four decimals. FSEOF is reproducible to solver
  tolerance, not bit-exact.
- MTA's source state is deterministic E-Flux2 at `objective_fraction=1.0` rather than
  contextualization plus sampling, and rMTA's ε is a fixed scalar rather than derived from a
  sampled reference distribution. Both are determinism choices; the external-`FluxState` route
  remains available.
- E-Flux2's per-condition normalisation means a uniform global expression shift produces zero
  predicted change in `predict_condition_fluxes`.
- Sampling runs `processes=1` by design, which disables optGpSampler's parallel-chain design;
  `thinning=100` is below what genome-scale sampling needs.
- `gene_perturbations` drops inert genes (71 of 137 on `e_coli_core`); the dropped count is now
  reported in provenance rather than being silent.
- The media presets are not M9 and cannot be made M9 — a literal M9 element set gives growth 0.0
  on iJO1366. `glucose_aerobic` matches iJO1366's shipped default medium exactly; that artefact,
  not a published composition table, is what is cited.

## 0.3.0 — 2026-07-10

### Scientific correctness

- Connected `robustknock` to StrainDesign's actual three-level `ROBUSTKNOCK` module and
  evaluated maximum and guaranteed product at exactly optimal growth.
- Reimplemented `mta` and `rmta` with the published MIQP/MOMA/MIQP workflow, L1
  transformation score, and rMTA Equation 9. Renamed the historical QP approximation to
  `rmta_continuous`.
- Made E-Flux2 a strict two-stage QP with a full-objective default and an explicitly labeled,
  opt-in L1 fallback.
- Aligned FSEOF/FVSEOF scan origins, enforced product/biomass constraints, trend handling,
  and actionable-target filtering with their documented method contracts.

### Reliability and validation

- Added finite/range/completeness checks across conditions, media, expression, flux states,
  perturbations, and scan parameters.
- Added deterministic model fingerprints and solver/runtime/parameter provenance to
  numerical result objects.
- Added direct COBRApy comparisons, the official MTA test topology, non-optional iJO1366
  checks, scientific sensitivity tests, and malformed-input/GUI-state regressions.
- Added an 80% branch-coverage gate, Ruff formatting/linting, mypy checks for scientific
  services, citation validation, and locked-dependency vulnerability auditing.

### Reproducibility and distribution

- Bumped the package to 0.3.0, constrained supported dependency/API ranges, and committed
  `uv.lock` for Python 3.10–3.12.
- Added `CITATION.cff`, an MIT `LICENSE`, scientific validation notes, and
  corrected solver/method documentation.
- Release builds now verify tag/version agreement, rerun all quality gates, validate wheel
  metadata, and install the built wheel in a clean environment.

### Behavior changes

- Code that relied on E-Flux2 silently falling back to pFBA must now request
  `allow_l1_fallback=True` and handle the `eflux2_l1_fallback` method label.
- `rmta` now requires MIQP and denotes the published robust workflow. Use
  `rmta_continuous` only when the historical QP heuristic is intentionally desired.
- Non-optimal FBA results return no objective value or flux vector instead of exposing
  solver-generated invalid numbers.
