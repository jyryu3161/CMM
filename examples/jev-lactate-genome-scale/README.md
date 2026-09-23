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

Live, `typesafe/jev-1.13-20260917`, six rounds, 26 steps, 43 calls, $0.027, about three
minutes:

| method | guaranteed D-lactate | cost |
|---|---:|---|
| OptKnock, 3 knockouts | **0.0000** | 2773 s |
| **JEV agent**, 6 rounds | **0.0000** | 43 calls, $0.027 |
| greedy over the same vocabulary | 0.0000 — stops on the plateau | seconds |
| **exhaustive pair sweep** | **17.5858** (`ALCD2x`✗ + `ATPS4rpp`✗) | ~350 s per anchor |

**The agent did not cross the plateau.** It reached pFBA 17.6 repeatedly — `PFL`✗ then
`GLCptspp`✗ — and every one of those designs guarantees zero, so none was ever promoted and
the run ends with `best_design: []`. It never reached `ALCD2x` + `ATPS4rpp`.

So on this problem every method that reasons fails and only brute force succeeds. That is a
result, not a failed run, and it is the first problem here sharp enough to produce one.

One thing the run changed about CMM. On the first attempt the agent was told *"product rose by
+17.2"* after each move, because the verdict quoted the pFBA product while the run ranks on the
guarantee — a design climbing impressively on a number nobody was scoring. It ended all six
rounds satisfied. The verdict now states the scored quantity, and on the same problem the agent
visibly changes behaviour: it starts *withdrawing* moves it had applied, which it never did when
the feedback flattered it. It still does not find the pair, but it is now failing against an
honest signal.

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
~1400 s each, plus a MOMA-L2 screen over 1367 genes. Expect the comparison to take longer than
the game. Set it to `false` for a first look.

## The key

The agent needs `OPENROUTER_API_KEY`; nothing else in CMM does. Either export it for one shell,
or save it once so later runs pick it up:

```bash
mkdir -p ~/.config/cmm && printf '%s\n' 'sk-or-...' > ~/.config/cmm/openrouter.key
chmod 600 ~/.config/cmm/openrouter.key
```

The file is plain text at `0600`; `cmm.jev.credentials.clear_key()` removes it.

## What a run here does not establish

Everything SC-03 says still applies, and one thing more. The agent's choices are not
guaranteed to repeat, so a single run that crosses the plateau is evidence about that run, not
about the method. Crossing it once is interesting; crossing it in most of *n* runs, with the
distribution stated, is a result. This example does not provide that.
