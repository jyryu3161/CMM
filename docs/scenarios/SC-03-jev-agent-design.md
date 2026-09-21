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

A **round** is one game. A **step** is one decision, and every decision costs one — an
intervention, an undo, a scan that changes nothing. `steps_per_round` is how long the agent
may play; a round ends when its steps run out or the agent chooses `end_round`.

`max_interventions` is a different budget entirely: how many changes may be active at once,
which is the number a laboratory would have to build. A long game and a small design is the
usual combination, because steps are cheap and edits are not.

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
| measured amplification gain (after `amplification_screen`) | what forcing flux through it actually does to the product, solved rather than guessed |

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

| Move | Meaning |
|---|---|
| `knockout` | flux forced to zero |
| `knockdown_50`, `knockdown_25` | capped at that fraction of the **wild-type** flux |
| `amplify_2x`, `amplify_5x` | at least that multiple of the wild-type flux, direction preserved |
| `force_on_low`, `force_on_high` | for a reaction carrying **no** wild-type flux: 25% or 60% of its loop-free maximum |
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

One measurement is deliberately **not** a move. What forcing flux through each candidate would
do to the product is recomputed whenever the design changes, because it is a fact and not a
decision. Left as a move the agent could choose, it was skipped: given the cofactor reading it
would infer a plausible answer — NADPH is short, so over-express an NADPH-producing enzyme —
and act on the inference instead of the measurement. Partial information displacing
measurement is worse than no information.

A move refused for dropping growth below the floor says **"too strong, not wrong"** when a
gentler version of it is still available on that reaction. Without that line, a refused
`force_on_high` on the glyoxylate shunt sent the agent to a different reaction and left behind
the 8% that `force_on_low` on the same one collects.

The relative moves and the `force_on` moves are offered to disjoint sets of reactions, so no
move ever means two things.

### Why `force_on` uses a loop-free maximum

A plain LP maximisation of `FRD7` on anaerobic `e_coli_core` returns its 1000 bound, reached
through the thermodynamically infeasible `FRD7`/`SUCDi` cycle, on a model taking up
10 mmol gDW⁻¹ h⁻¹ of glucose. Sixty per cent of that would be a physically meaningless target
the solver would satisfy with a futile cycle. The loopless range gives 13.6, and costs less
time than the plain one because the tighter problem solves faster.

## Reading a result

`00_summary.json` and the CLI output give the headline. `best_product_flux` is the product
flux at the **pFBA** optimum of the best design that held the growth floor — the adapted
strain, and the right score, because a design that only pays off if the cell chooses to make
the product is not a design. `moma_product_flux` in `02_game/ticks.csv` is the unadapted
state immediately after the change.

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
| **JEV agent**, 4 interventions | **10.761** | 0.068 | 7 s, $0.002 | **no** |

Read this carefully, because the obvious reading is wrong in both directions.

**A single gene deletion cannot solve this problem at all.** The best of 71 reaches 0.211.
Anaerobic succinate needs several routes closed at once, and no single-deletion screen —
however it is scored — can find that.

**OptKnock is not beaten by an agent searching on its own.** Before the deterministic designer
was wired into the loop, JEV runs landed between 1.9 and 9.1 and never matched 9.911. The
reason is structural and worth stating: OptKnock's winning design deletes `LDH_D` and `THD2`,
neither of which carries any flux in the wild type. They are escape routes the cell would
switch to once the obvious ones are shut, and a board built from where the flux is *today*
cannot see them. No amount of play fixes that; it is a blindness in what the agent is shown.

**What the agent adds is a move OptKnock cannot express.** OptKnock searches deletions only.
Handed its own proven design and one place left in the budget, the agent forces flux through
the glyoxylate shunt (`ICL`) and reaches 10.761 — 8.6% above the deterministic optimum, at a
real cost in growth (0.068 against 0.091, both above the floor). Verified independently: the
result lies inside the **loop-free** feasible succinate range for that design, so it is not a
thermodynamic artifact.

**Over eight independent runs of one configuration**: all eight reached 10.761, in six steps,
five seconds and $0.0011 each. Standard deviation 0.0. The agent does not do worse than the
deterministic method because it starts from it.

Five things made the difference, and each was a failure before it was a fix:

- **Seeding.** `seed_with_strain_design` runs the designer once before the first move and puts
  the reactions it names on the board with the guaranteed product they buy.
- **Adopting a design as one move.** A design's deletions only pay off together — `ACALD`
  alone buys almost nothing — so an agent judging each move by the product change it causes
  abandons the design after the first deletion. `adopt_best_design` applies the set.
- **Measuring instead of guessing.** Asked which reaction to amplify, the agent reached for
  fumarate reductase, the direct product-forming step, which is already saturated and buys
  nothing. CMM now solves for the answer — one pFBA per candidate, about a second for a board
  of 24 — and puts the measured change on the board.
- **Making that measurement automatic.** While it was still a move the agent could choose, it
  stopped choosing it once the cofactor reading gave it something plausible to infer from.
- **Saying when a move was too strong rather than wrong.** `force_on_high` on `ICL` breaches
  the floor and `force_on_low` on the same reaction is the 8%. Naming the gentler move in the
  rejection took the outcome from nine runs in ten to ten in ten.

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
