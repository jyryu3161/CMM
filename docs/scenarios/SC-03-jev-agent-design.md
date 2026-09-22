# SC-03 — JEV agent design

What the numbers in a `cmm jev-design` run mean, what the run does and does not establish,
and where it differs from every other scenario in this directory.

## What this is

CMM's other target-discovery methods are deterministic optimisations. FSEOF scans enforced
product flux; OptKnock solves a bilevel MILP; MOMA finds the nearest feasible state. Each
answers a precise mathematical question and answers it the same way every time.

SC-03 is not one of those. It puts a **decision model in the loop**: TypeSafe's JEV, reached
through OpenRouter, is shown the metabolic state and picks one move at a time from a fixed
vocabulary, many times over. CMM executes each move, re-solves, and shows it the result.

JEV does not generate text. Every answer is drawn from criteria CMM supplied, so it cannot
name a reaction the model does not contain — there is no output to parse and nothing to
sanitise. Each answer also carries a probability for every option it was offered, so one call
ranks the entire board.

## The division of labour

This is the only reason a run is worth reading.

**CMM owns the rules.** Every number on the screen is a solve. A move that makes the model
infeasible, or that pushes growth below the configured floor, is reverted by CMM whatever the
agent predicted, and the reversal is a row in `02_game/ticks.csv` with its reason.

**JEV owns the strategy.** A move that is legal but unhelpful is kept and flagged in the
history the agent reads next tick, for it to withdraw with `undo_last`. Auto-reverting
everything that failed to improve the product would make the loop greedy hill climbing and
would make a two-step manoeuvre impossible to express.

## Rounds and steps

**A round is one independent attempt.** It starts from the wild type with an empty design,
plays until its steps run out or the agent ends it, and is scored on its own. What carries
across rounds is knowledge, not bounds: the agent is shown what each earlier round reached
and with which design, so a later round can go after something different, or head for what
already worked with `restore_best_design`.

Rounds used to continue one another, and the effect was not subtle — the first round filled
the design and the rest had nothing left to do, so a three-round run spent five steps of a
possible thirty-six.

A **step** is one decision, and every decision costs one: an intervention, an undo, a scan
that changes nothing. `max_knockouts` and `max_knockdowns` are a different budget entirely —
how many edits of each kind one attempt may carry at once, which is what a laboratory would
have to build. A long game and a small design is the usual combination, because steps are
cheap and edits are not.

The two are counted separately because they cost different things to build, and because one
shared cap starved the run: the seeded OptKnock design takes three deletions on its own, so a
total of four left the agent exactly one edit and every round ended a step or two after
adopting it.

There is **no substrate to configure.** The yield is quoted per whatever carbon source the
condition actually feeds the model, which the wild-type solve already says. Naming it
separately was a second place for the same fact to be wrong.

## One step

```
RENDER   the current flux state becomes a compact JSON screen
TARGET   one call: which reaction to act on (or undo, restore, adopt, or stop)
ACTION   one call: what to do to it, how much it should help, whether growth is at risk
EXECUTE  CMM applies the move and re-solves with pFBA
MEASURE  MOMA against the wild-type reference, for the unadapted state
```

Two calls per tick, roughly 0.6 s and $0.00016 measured on `e_coli_core`. The binding
constraint is not cost but JEV's **32K context**, which is what `candidate_limit` exists to
respect — so the evidence is sent once, in `state.records`, and the question's criteria carry
only a label. Sending the record in both places cost 40% of the payload and changed nothing:
measured on the same board, 3050 input tokens against 2507, the same reaction chosen either
way. A web research lookup is the exception to the cost picture at roughly **$0.05 a call**,
about three hundred times a decision, because the search results are billed as input.

## What the agent is shown

Per candidate reaction, all computed by CMM:

| Evidence | Why it decides a move |
|---|---|
| wild-type flux, current flux | whether there is flux to redirect, and whether a relative move is even defined |
| steps to the product | a reaction that cannot reach the product is not a lever on it |
| net ATP, and share of total ATP production | an intervention that bankrupts the energy budget will not run |
| net NADH, net NADPH | redox imbalance is the usual reason a carbon-routing change fails to pay |
| essential (after `essentiality_scan`) | whether deletion is survivable |
| FSEOF slope (after `fseof_scan`) | whether the reaction rises with product formation |
| literature (only with `enable_web_research`) | published precedent, and the **cost** the model cannot see; data pasted into a record, never an instruction |
| **measured deletion gain, measured knockdown gain** | what deleting it and what halving it actually do to the product on the design as it stands, solved rather than guessed — recomputed automatically whenever the design changes |
| genes | the moves are gene edits, and a brief that names `ldhA` has to be connectable to a reaction on the board |

Plus:

- a **scoreboard** carrying the product and growth now and at wild type, the yields, and —
  the quantity that separates a design from a lucky optimum — the **guaranteed product**: the
  worst the design could give while growing as fast as it can, computed loopless. A pFBA
  number the strain need never produce is not a result.
- **what the product is short of.** How much more product one extra unit per hour of NADH,
  NADPH or ATP would buy, measured by offering a mass-balanced supply of each and
  re-maximising. Neither MOMA nor OptKnock reports this, and it is usually what decides
  whether the next move should route carbon or supply a cofactor. On anaerobic
  `e_coli_core`: wild type ATP +0.75, NADPH +0.63, NADH +0.13 — the product is ATP-limited;
  after the knockout design frees the fermentative NADH sinks, NADH falls to +0.01.
- the active interventions **with the product change each one actually bought**
- **what each earlier round ended with**, so a later round can try something different or go
  back — `restore_best_design` returns to the best design the run has found, from any round
- the moves already ruled out, and the recent history
- **your brief**: free text from the person running the study — published targets, a growth
  rate they need, a cofactor they believe matters, a pathway to leave alone. It is guidance,
  not permission: it cannot widen the move vocabulary, name a reaction outside the model, or
  lift the growth floor CMM enforces, so the worst a mistaken brief can do is waste steps.

## Another organism, another product

Nothing in the agent's vocabulary is specific to *E. coli* or to succinate. The moves are
metabolic — carbon, reducing power, ATP, growth — and the product, the growth floor and the
board are all parameters. Two things had to be fixed before that was actually true:

**Cofactor pools are found by formula, not by id.** Yeast-GEM calls ATP `s_0434`; AGORA
models differ again. An implementation that matches `atp_c` does not *fail* on those — it
silently finds nothing, drops the cofactor reading, and treats every hub metabolite as a
carbon carrier so the graph distances become noise. Pools are identified instead by the part
of a formula that does not change with protonation or naming: the counts of carbon, nitrogen,
phosphorus and sulfur. Sulfur is in the key because coenzyme A and NADP share `(21, 7, 3)`
and differ only by it, and within a redox pair the reduced member is the one carrying one
more hydrogen — matching the skeleton alone makes NAD and NADH cancel, so a reaction that
produces one NADH reads as redox-neutral.

On `e_coli_core` the formula route reproduces the id route exactly: same net coefficient on
every reaction, same graph distances, same turnover figures. A pool that cannot be identified
is **named on the screen** under `not_measured`, because an absent row reads as a zero.

**The organism is required for a literature lookup and has no default.** Asking the published
record about the wrong species returns an answer that is confident and wrong, so
`enable_web_research` without `organism` is refused rather than guessed at.

## The board

Composed from four slates rather than one ranking, because one ranking demonstrably fails:

- **pinned** — reactions a deterministic strain designer named, after `strain_design_scan` or
  the seeding that runs it before the first move. These lead the board, because a proof that
  deleting a reaction forces the product outranks every other kind of evidence on it.
- **near** — closest to the product. These open the route.
- **competing** — reactions feeding the carbon byproducts the strain is currently secreting.
  These close the competition. On anaerobic *E. coli* they are the formate, acetate and
  ethanol branches, and no proximity or flux ranking surfaces them.
- **carriers** — the remaining largest flux carriers.

**Only reactions with a gene association are candidates.** A move must be one a laboratory
could make; an exchange reaction, a biomass pseudo-reaction or an ATP maintenance term has no
gene to delete or over-express.

## The moves

The vocabulary is **down-regulation only**: delete a gene, or weaken it to half.

| Move | Meaning |
|---|---|
| `knockout` | delete the reaction's genes: flux forced to zero |
| `knockdown_50` | weaken them: capped at half the **wild-type** flux magnitude |
| `undo_last` | withdraw the most recent intervention |
| `fseof_scan`, `essentiality_scan`, `envelope_probe` | run a CMM analysis; the model is unchanged |
| `state_distance_check` | MOMA and ROOM on the current design: how far the cell has to move, and how many reactions have to change |
| `strain_design_scan` | run OptKnock and RobustKnock and put the reactions they name on the board |
| `adopt_best_design` | apply a proven knockout set in one move, since its deletions only pay off together |
| `restore_best_design` | discard the current design and go back to the best the run has found |
| `end_round` | stop spending steps |

There is **one knockdown strength**, not two. A second, deeper cap mostly bought a second
rejection of the same idea, and roughly halving an activity is the level a promoter swap or an
RBS change can actually aim at.

A `knockdown_50` needs a non-zero wild-type flux to be half of, and is not offered without
one. A `knockout` is always offered, including on a reaction carrying nothing today — which is
not a pointless move: OptKnock's most valuable deletions close routes the cell would only
switch to once the obvious ones are shut.

### Why there is no amplification

The vocabulary used to include `amplify_2x`, `amplify_5x` and a `force_on` move that switched
a zero-flux reaction on at a fraction of its feasible maximum, and that is what produced this
project's best succinate design: **10.7613** against OptKnock's 9.9108. It was removed anyway.

A lower bound on a reaction is not what over-expression does. Forcing `v ≥ x` tells the solver
the cell **must** carry that flux, and the solver will satisfy it by whatever route is
cheapest — including one the enzyme has nothing to do with. Stronger expression of an enzyme
raises a *capacity*; the cell still decides whether to use it. The in-silico gain from a
forced lower bound is therefore an upper bound on an upper bound, and it is the kind of number
that survives review and fails in a flask. A deletion and a knockdown are caps: they say what
the cell **cannot** do, which is what deleting a gene or weakening its promoter achieves.

The restriction costs product, and every run prices it rather than arguing about it. The
baseline table carries a
`best amplification on top of this design (outside the vocabulary)` row: FSEOF is ranked on
the design being scored, its top ten targets are forced one at a time, and the best one that
still clears the growth floor is reported. The verdict excludes that row from "best
deterministic method" — scoring the agent against a move it was forbidden to play is not a
comparison — and names it in the sentence instead.

It is measured per run, not quoted as a constant, because it is not one. Ranked on the wild
type, FSEOF's top amplification target for anaerobic succinate buys **nothing**. Ranked on a
design that already deletes `ACALD`, `D_LACt2` and `THD2`, the same method puts the glyoxylate
shunt fourth and forcing it reaches **10.76** — but at a growth rate of 0.041, so at a floor
of 0.05 that move does not exist, and the best one that does reaches **10.04** against the
design's 9.95. One design, three defensible numbers; a headline percentage would have been
true of one of them.

### The screen: a measurement, not a move

One measurement is deliberately **not** a move. What deleting and what halving each candidate
would do to the product is recomputed whenever the design changes — two pFBA solves per
candidate, about two seconds for a board of 24 on `e_coli_core` — because it is a fact and not
a decision. Left as a move the agent could choose, it was skipped: given the cofactor reading
it would infer a plausible answer and act on the inference instead of the measurement. Partial
information displacing measurement is worse than no information.

Both moves are measured because the pair is the decision: a reaction whose deletion is lethal
and whose knockdown pays is exactly what the knockdown exists for, and screening deletions
alone would hide it. Essentiality comes free with the deletion solve, so `essentiality_scan`
is withheld once the screen has run — its answer is already on the board.

A move refused for dropping growth below the floor says **"too strong, not wrong"** when a
gentler version of it is still available on that reaction. There is exactly one such pair,
`knockout → knockdown_50`, and it is the one that matters: a gene the cell cannot live without
can very often live at half.

### Which state MOMA measures from

MOMA's reference here is the **wild type**, for every design at every step, and that is a
choice worth stating because the loop makes its edits one at a time and the obvious
alternative is to re-reference each step to the design before it.

MOMA's premise ([Segrè et al. 2002](https://www.pnas.org/doi/full/10.1073/pnas.232349399)) is
that a freshly perturbed cell keeps the regulatory setpoints of the cell it was made from, so
the reference must be the *parent strain*. The parent of the design this run produces is the
wild type: the design is built and characterised as one strain, and the intermediate designs
the agent passes through on the way are search positions, not organisms anyone will culture.
The published convention agrees — epistasis maps built with MOMA reference single and double
mutants alike to the wild type, not each single mutant to its own parent.

Re-referencing each step to the previous design models a different experiment — an edit
introduced into a strain that has already been grown up — and it is a real protocol. It is
still the wrong number to report here, for two reasons:

- **It destroys the quantity.** Chained, each step's MOMA distance describes only the last
  edit, so a five-edit design looks exactly as easy to build as a one-edit design. The whole
  point of the number is that it does not.
- **It would not be MOMA-from-MOMA in any case.** A strain you can make a second edit in is a
  strain you have cultured, and cultured knockout strains move *away* from the MOMA state
  toward the FBA optimum — Fong & Palsson evolved knockout strains to within 10% of the
  predicted optimum in 38 of 50 cases. So the honest parent state would be the previous
  design's **pFBA**, which is the quantity the loop already reports as the score.

If you want the sequential reading, it is the difference between consecutive `product_flux`
rows in `02_game/ticks.csv`; the `moma_*` columns are deliberately not that.

## Reading a result

`00_summary.json` and the CLI output give the headline. `best_product_flux` is the product
flux at the **pFBA** optimum of the best design that held the growth floor — the adapted
strain, and the right score, because a design that only pays off if the cell chooses to make
the product is not a design. `moma_product_flux` in `02_game/ticks.csv` is the unadapted
state immediately after the change, against the wild type.

Read these together with the tick table:

- an **`applied`** row whose reason says *the forced flux was consumed inside the network* is
  a move that changed nothing useful: the constraint was met by an internal cycle. The move
  still occupies one of the design's places.
- a **`reverted_growth_floor`** row is not a failure of the method. It is the rule working.
- an **`infeasible`** row is data. A lethal intervention is a finding.
- `02_game/candidate_rankings.csv` holds the agent's probability over every option at every
  tick. The runner-up is often the more interesting row.

## Measured against the deterministic methods

Every JEV run scores itself against the methods it sits beside, on the same model, the same
condition and the same growth floor, with every design evaluated the same way — apply it,
solve pFBA, read the product and the growth. The table is `05_baseline/comparison.csv` and
the verdict is in `00_summary.json`.

On the shipped anaerobic succinate example:

| Method | Best succinate (mmol gDW⁻¹ h⁻¹) | Growth (h⁻¹) | Time | Deterministic |
|---|---:|---:|---:|:---:|
| Wild type | 0.000 | 0.212 | — | — |
| Best single gene deletion, MOMA-L2, 71 genes | 0.211 | 0.165 | 3.4 s | yes |
| **OptKnock**, 3 knockouts | **9.911** | 0.091 | 0.8 s | yes |
| **RobustKnock**, 3 knockouts | **9.911** | 0.091 | 0.9 s | yes |
| **JEV agent**, 3 deletions + 1 knockdown | **9.946** | 0.055 | 29 steps, $0.009 | **no** |
| *Best amplification on top of that design* | *10.039* | *0.053* | — | *outside the vocabulary* |

Read this carefully, because the obvious reading is wrong in several directions.

**A single gene deletion cannot solve this problem at all.** The best of 71 reaches 0.211.
Anaerobic succinate needs several routes closed at once, and no single-deletion screen —
however it is scored — can find that.

**OptKnock is not beaten by an agent searching on its own.** Before the deterministic designer
was wired into the loop, JEV runs landed between 1.9 and 9.1 and never matched 9.911. The
reason is structural and worth stating: OptKnock's winning design deletes `LDH_D` and `THD2`,
neither of which carries any flux in the wild type. They are escape routes the cell would
switch to once the obvious ones are shut, and a board built from where the flux is *today*
cannot see them. No amount of play fixes that; it is a blindness in what the agent is shown.

**What the agent adds is a move OptKnock cannot express.** OptKnock's variables are
present-or-absent: a 50% cap is not a constraint its formulation can write down. Handed its
own proven design, the agent halved acetate kinase (`ACKr`) on top of it and reached 9.946.
That is **+0.4%** — a small margin honestly won, on a problem where the deterministic method
is already close to the ceiling, and bought with growth (0.055 against 0.091, both above the
floor).

**The margin used to be 8.6%, and the vocabulary change is why.** Forcing flux through the
glyoxylate shunt on top of OptKnock's design reaches 10.761, and that is what earlier runs
did. Amplification was removed anyway, for the reason in *Why there is no amplification*
above. The headroom row keeps the consequence visible in every run instead of leaving it as
an argument.

**The agent found the optimum of the space it was given.** An exhaustive screen of every
gene-associated deletion and every 50% knockdown on top of the OptKnock set reaches 9.9461 —
the same number, via `ACKr` or the equivalent `PTAr`. That is worth more than the margin: it
says the loop is searching its space properly, not that this space is the best one.

Four things made the difference, and each was a failure before it was a fix:

- **Seeding.** `seed_with_strain_design` runs the designer once before the first move and puts
  the reactions it names on the board with the guaranteed product they buy.
- **Adopting a design as one move.** A design's deletions only pay off together — `ACALD`
  alone buys almost nothing — so an agent judging each move by the product change it causes
  abandons the design after the first deletion. `adopt_best_design` applies the set.
- **Measuring instead of guessing.** Which branch competes with the product depends on the
  whole network at the current bounds, not on a reaction's own stoichiometry. CMM solves for
  it — two pFBA solves per candidate — and puts the measured change on the board.
- **Making that measurement automatic.** While it was still a move the agent could choose, it
  stopped choosing it once the cofactor reading gave it something plausible to infer from.

**This is one problem on one small model.** It is evidence that the loop can add something to
a deterministic optimum on a problem where the two intervention classes are complementary. It
is not evidence that it will on yours. Run the comparison; it is on by default.

## What a run does **not** establish

- **It is not reproducible in the agent's choices.** The CMM solves are deterministic; JEV's
  answers are not guaranteed to repeat, and two runs of the same config have produced designs
  differing several-fold in product flux. `04_agent/transcript.jsonl` records every request
  and response so that *one* run can be audited. Never present a single run as the method's
  performance.
- **It is not a benchmark against FSEOF or OptKnock.** Comparing them needs
  `cmm production-targets` on the same model and condition, and a comparison of the two run
  directories. That is a separate piece of work and this scenario does not do it.
- **It is not a wet-lab claim.** Every design here is a computational hypothesis to test
  experimentally, exactly as in SC-01 and SC-02.

## Stop and ask

In addition to the shared preflight rules in `_preflight.md`:

- **No growth floor is agreed.** It is the only rule CMM enforces against the agent, so it
  decides which designs can exist. It is not a default to inherit silently.
- **Theoretical yield is zero for the product.** No sequence of moves can change that; report
  it and ask about the medium, substrate or aeration.
- **The user wants a reproducible result.** Say plainly that the agent's choices are not
  guaranteed to repeat, and offer `cmm production-targets` instead, which is deterministic.

## Run contract

```
<run directory>/
  00_config.json         resolved inputs, paths rewritten to be bundle-relative
  00_provenance.json     model fingerprint, solver, versions, question-set version,
                         the JEV model id that actually served the run
  00_summary.json        headline result and the best design
  00_manifest.json       authoritative artifact inventory (schema_version 2)
  model/<model-id>.xml   the model that was solved
  01_wild_type/          reference pFBA fluxes and the starting scoreboard
  02_game/               ticks.csv, rounds.csv, candidate_rankings.csv
  03_design/             best_design.csv, final_design.csv, flux_trajectory.csv
  04_agent/              transcript.jsonl (every request and response), usage.json
```

Units are CMM's throughout: fluxes in mmol gDW⁻¹ h⁻¹, growth in h⁻¹, molar yield in mol/mol.

There is **no publication renderer and no completion gate** for SC-03. `cmm report validate`
refuses a JEV run rather than measuring it against the production contract it was never meant
to satisfy.
