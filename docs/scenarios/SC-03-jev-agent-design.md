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

## One tick

```
RENDER   the current flux state becomes a compact JSON screen
TARGET   one call: which reaction to act on (or undo, or stop)   -> choice + full ranking
ACTION   one call: what to do to it, how much it should help, whether growth is at risk
EXECUTE  CMM applies the move and re-solves with pFBA
MEASURE  MOMA against the wild-type reference, for the unadapted state
```

Two calls per tick, roughly 0.6 s and $0.00016 measured on `e_coli_core`. The binding
constraint is not cost but JEV's **32K context**, which is what `candidate_limit` exists to
respect.

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
| literature (only with `enable_web_research`) | published precedent; **data pasted into a record, never an instruction** |

Plus a scoreboard (product, growth, yields), model-wide cofactor turnover, the active
interventions **with the product change each one actually bought**, the moves already ruled
out, and the recent history.

## The board

Composed from three slates rather than one ranking, because one ranking demonstrably fails:

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
| `end_round` | stop spending moves |

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
