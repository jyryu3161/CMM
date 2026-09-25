# SC-03 on a genome-scale model: D-lactate in *E. coli* iJO1366

This is the run the `e_coli_core` succinate example cannot be: a problem where **no single
move works at all**, so a search that takes the best move each step has nothing to take.

```bash
export OPENROUTER_API_KEY=...   # or save it once: see "The key" below
uv run --frozen --all-extras python examples/jev-lactate-genome-scale/prepare.py
uv run --frozen --all-extras cmm jev-design --config examples/jev-lactate-genome-scale/config.json
```

`prepare.py` unpacks the iJO1366 that COBRApy ships (2583 reactions, 1367 genes) rather than
committing a second copy, the same arrangement the SC-01 and SC-02 examples use.

## Why this problem

Measured on this model and condition, before any agent was involved:

| | best single move | best pair | gain from the second |
|---|---:|---:|---:|
| succinate, `e_coli_core`, on OptKnock's design | 9.9457 | 9.9457 | +0.000000 |
| succinate, `iJO1366`, on OptKnock's design | 12.4015 (`PGI`) | 12.4022 | +0.000727 |
| **D-lactate, `iJO1366`, from the wild type** | **0.0000** | **17.5858** | **+17.585806** |

Every one of 2526 single knockouts and knockdowns leaves D-lactate at exactly zero. The pair
`ATPS4rpp` + `ALCD2x` reaches 17.59. Neither edit pays alone; together they are the design.

Confirmed against the whole board, not a shortlist: anchoring on each of those two moves and
trying **every** partner among 2123 gene-associated reactions in both actions — about 1790
partners give product per anchor — the best guaranteed pair is still
`knockout:ALCD2x` + `knockout:ATPS4rpp` at **17.5858**. What that sweep also turns up is a
portfolio, with the growth trade-off laid out:

| design | guaranteed | growth |
|---|---:|---:|
| `ALCD2x` ✗ + `ATPS4rpp` ✗ | 17.5858 | 0.1625 |
| `ALCD2x` ↓50% + `ATPS4rpp` ✗ | 13.5470 | 0.1625 |
| `ALCD2x` ↓50% + `PFL` ✗ | 13.1526 | 0.1892 |
| `ALCD2x` ↓50% + `ACKr` ✗ (or `PTAr` ✗) | 12.3947 | 0.1972 |
| `ATPS4rpp` ↓50% + `ACALD` ✗ | 11.8282 | 0.1948 |

Four of those five **mix a knockdown with a knockout**, which is exactly the shape OptKnock's
present-or-absent formulation cannot write down. Trading 5.2 of product for 0.035 h⁻¹ of growth
is the kind of choice a laboratory makes, and it only exists because the vocabulary has both
moves in it.

**The deterministic designer does not solve this problem.** OptKnock on this model and
product, `max_knockouts=3`, `max_solutions=5`, `min_growth=0.05`, returned after **2773 s**
with five designs, and every one of them guarantees **0.0000**:

| OptKnock's answer | growth | pFBA | guaranteed |
|---|---:|---:|---:|
| `ATPS4rpp`, `PPC`, `TPI` | 0.0504 | 9.8924 | **0.0000** |
| `GLCptspp`, `PFL`, `TPI` | 0.0539 | 9.8790 | **0.0000** |
| `PFL`, `PPC`, `TPI` | 0.0545 | 9.6982 | **0.0000** |
| *(two more, same)* | | | **0.0000** |
| **`ALCD2x`, `ATPS4rpp`** — which it did not return | **0.1625** | 17.5974 | **17.5858** |

That last row is **two** knockouts, inside OptKnock's own three-knockout budget, and it is
better on every axis: 17.59 guaranteed against 0, at three times the growth. Each row above was
re-scored here independently of the designer, and the result is the same applied as bare
reaction knockouts or as the gene edits that achieve them.

Why the designer missed it is *not* established here — a node or time limit inside
`straindesign`, a compression artifact, or something in the formulation are all live
possibilities, and none of them was checked. What is established is that on this problem the
method the comparison treats as the standard to beat returns nothing usable after 46 minutes,
while a two-edit design that guarantees 17.59 sits inside its stated search space. Read the
`OptKnock` and `RobustKnock` rows of `05_baseline/comparison.csv` on this run with that in mind.

**And the guarantee is doing the work.** Ranked by pFBA the top partners for `ATPS4rpp` are
`GLCptspp`, `PGI`, `CHTBSptspp` and friends at 19.02 — every one of which has a guaranteed
product of **0.0000**. The strain could make 19 and need never make any. A screen ranked on the
pFBA optimum would have reported four spectacular designs that are not designs.

That is the whole point of running this one. A greedy search — which is what the
`best deterministic design + N knockdowns` control row is, and what a hill-climbing agent is —
stops dead on that plateau and reports nothing. The control says so in its own note when it
happens. So on this problem the control is a floor and not a ceiling, and the question "does
judgement buy anything the search would not find" finally has room to be answered either way.

The `brief` in the config states the plateau to the agent in words, because that is what a
person who had run the screen would tell it. It is guidance, not permission: it cannot widen
the vocabulary or lift the growth floor.

## What happened when it was run

It was run twice over, because the first attempt found a defect in CMM rather than a design.

**First attempt — the agent returned nothing, ten times out of ten.** `best_design: []`, and
every replicate took the identical path: `PFL`✗ then `GLCptspp`✗, first move 10/10 the same,
nine distinct reactions touched across all ten runs out of 2123. Two causes, both now fixed:

1. *It was told the wrong number.* The verdict after each move quoted the change in the pFBA
   product — "product rose by +17.2" — while the run ranks designs on the guarantee, which never
   left zero. It ended all six rounds satisfied. The verdict now states the scored quantity.
2. *The moves it needed were never on the board.* `ALCD2x` and `ATPS4rpp` — the only pair that
   guarantees any D-lactate — were **never offered once in ten runs**. The board's slates are
   built on the carbon graph and on flux magnitude; `ATPS4rpp` is ATP synthase and has no path
   to the product at all, so no board size could reach it. Meanwhile the run's own
   `cofactor_limitation` was reporting ATP as limiting by +3.0 on every tick. There is a
   cofactor slate now, and this example sets `candidate_limit: 32`.

**Second attempt — it beats everything deterministic here, in nine runs out of nine that
finished.** Ten replicates through `evals/jev_replicates.py`:

| method | guaranteed D-lactate | edits | growth | cost |
|---|---:|---:|---:|---|
| OptKnock, 3 knockouts | 0.0000 | 3 | 0.050 | 2773 s |
| exhaustive **pair** sweep | 17.5858 | 2 | 0.1625 | ~350 s per anchor |
| **agent**, median of 9 | **17.9200** | 8 | 0.1340 | ~570 s, $0.077 |
| **agent**, best of 9 | **19.1361** | 6 | 0.0581 | 546 s, $0.070 |

| outcome | runs |
|---|---:|
| 17.9200 | 6 |
| 18.0150 | 1 |
| 19.1199 | 1 |
| 19.1361 | 1 |
| timed out at 1200 s | 1 |

**9 of 9 completions exceeded 17.5858**, the best design an exhaustive pair search can find. Six
distinct designs, $0.70 for the study. Read at one growth rate the gap is not a trade: held at
the best run's 0.0581 its design guarantees 19.14 where the pair guarantees 0.0000. Every design
quoted here was re-scored independently of the run that produced it.

None of them is the pair. The best is `ALCD2x`, `GLUDy`, `ATPS4rpp`, `GLCptspp`, `PPK`, `FADRx`
— four of those six are reactions only the cofactor slate can put on a board.

**What this is evidence of, and what it is not.** The exhaustive search that found 17.5858 was
over *pairs*, because pairs are where exhaustive search stops being affordable — 496 of them on
the reduced board here, 3.2 million over the whole model. A six- or eight-edit design was never
in its space. So this is not "the agent beat exhaustive search"; it is **the agent returning a
better design at a depth exhaustive search cannot reach**, which is the only place judgement
could have paid. It is ten runs of one configuration on one product and one model, and it says
nothing about the next one.

**One run in ten hung** and was recorded as a timeout rather than being allowed to stall the
study. See *What it costs*.

## What to read afterwards

- `00_summary.json` — `best_guaranteed_product` is the number that matters. Anything above
  zero means the agent crossed a plateau that one-move-at-a-time search cannot.
- `05_baseline/comparison.csv` — the control row's note will say whether its greedy search hit
  the plateau. If it did, the agent is being compared against a floor.
- `02_game/ticks.csv` — `runner_up` and `decided_by` beside every move. On a plateau every
  first move measures zero gain, so what the agent chose *on no measured evidence* is the
  interesting column.
- `report.html` — the whole run on one page.

## What it costs

Genome scale is not `e_coli_core`, and the deterministic parts dominate:

| step | `e_coli_core` | `iJO1366` |
|---|---:|---:|
| one pFBA | 0.01 s | 0.33 s |
| one loopless guarantee | 0.1 s | 0.8 s |
| OptKnock, 3 knockouts, 3 solutions | 1.2 s | **1417 s** |
| the intervention screen, 24 candidates | ~2 s | ~16 s |

`seed_with_strain_design` is **off** in this config. Seeding is what makes the agent reachable
on succinate, but here it would cost 1417 s before the first move and would hand the agent a
deterministic answer on the very problem meant to test whether it can find one itself.

`run_baseline_comparison` is on and is the slow part of the tail: OptKnock and RobustKnock at
~1400 s each on succinate and 2773 s on D-lactate, plus a MOMA-L2 screen over 1367 genes. Expect
the comparison to take longer than the game. Set it to `false` for a first look — the two runs
reported above did.

**And one LOOK move can dominate everything else.** Measured with baselines off: run 1 took
588 s over 56 steps, run 2 took **17272 s — 4.8 hours — over 86 steps**. The difference is not
the extra steps and it is not the design getting harder to solve. Profiled at every depth of
run 2's own design, the per-tick solver cost is flat at about 9.5 s:

| depth | pFBA | guarantee | cofactors | MOMA | screen (32) | total |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.13 | 0.40 | 0.07 | 1.95 | 6.69 | 9.48 |
| 3 | 0.11 | 0.39 | 0.08 | 2.37 | 6.13 | 9.38 |
| 6 | 0.11 | 0.45 | 0.06 | 2.43 | 6.21 | 9.50 |

86 steps of that is 820 s, and the decision calls added 66 s. The other **4.6 hours** were two
`strain_design_scan` LOOK moves: on this model that is OptKnock *and* RobustKnock over 2583
reactions, about 2.3 hours each. Run 1 made none, which is the whole of the difference. Neither
call put a single reaction into the design run 2 reported — the reactions it pinned were
`ACACT3r`, `ACOAD1f`, `PGCD` and friends, and `ATPS4rpp` had already been applied at R1T2.

It also should not have been able to make two. A scan is recorded as done, but `ScanCache`
snapshots are taken before a move and restored when it is withdrawn, which rolled that record
back. Permanent scans survive a restore now, and every LOOK move is timed, reported in its own
`reason`, and flagged on the board when it runs long — so a run can no longer spend 2.3 hours
without saying so.

The rest is honest cost. `screen_interventions` is 6 s of the 9.5, and turning it off or
shrinking `steps_per_round` changes what the agent sees, so say which you used.

**Two looks are unaffordable on this model and this config refuses them.** Measured here:
`fseof_scan` 1.0 s, `envelope_probe` 0.2 s, `strain_design_scan` about 2.3 hours, and
`state_distance_check` — which runs ROOM, a MILP over 2583 binaries — had not returned after
25 minutes. All five are sub-second on `e_coli_core`, which is the model the LOOK vocabulary was
designed against. `disabled_look_actions` names the two, which is also the honest setting when
`seed_with_strain_design` is off: calling the designer through a look is seeding by another
route. `max_scan_seconds` is the backstop for a model nobody has measured yet — a look that
overruns it is withheld for the rest of the run, so it is paid for once rather than every time
the agent asks.

**And `max_run_seconds` bounds the run.** Ten replicates of this config ran 10 minutes, 4.8
hours, and longer again on identical settings, because the paths differ; a study cannot be
planned against that. Stopping on the clock is the same answer the stop button gives — keep what
was played, score it, say so.

*It is checked between steps, not inside one*, so a single runaway step overruns it — and about
one run in ten does exactly that. It is not the design getting deeper (flat, above), not either
disabled look, not the loopless guarantee (0.43 s worst over 40 random designs) and not MOMA-L2
(2.6 s worst over the same). It remains unidentified, and rare enough that six traced runs in a
row completed cleanly without reproducing it.

So `evals/jev_replicates.py` runs each replicate in its own subprocess under `--timeout` and
records a hang as a data point. That is what makes a ten-run study finishable while the cause is
still open: run 4 of ten timed out at 1200 s and the other nine reported normally.
`02_game/ticks.csv` carries `elapsed_s` and `solver_s` per step, so a completed run says where
its time went without a sampler attached to a live process.

## The key

The agent needs `OPENROUTER_API_KEY`; nothing else in CMM does. Either export it for one shell,
or save it once so later runs pick it up:

```bash
mkdir -p ~/.config/cmm && printf '%s\n' 'sk-or-...' > ~/.config/cmm/openrouter.key
chmod 600 ~/.config/cmm/openrouter.key
```

The file is plain text at `0600`; `cmm.jev.credentials.clear_key()` removes it.

## What a run here does not establish

Everything SC-03 says still applies, and one thing more. The ten-run study quoted above is of
the **old** board: 0 of 10, 20 to 26 steps, $0.245 in total. On the fixed board only two runs
have been made, both successful but one of them 4.8 hours long, so there is no distribution here
yet — two successes are not a success rate. Repeat it before quoting it:

```bash
uv run --frozen --all-extras python evals/jev_replicates.py \
    examples/jev-lactate-genome-scale/config.json --runs 10 --target 17.5858
```
