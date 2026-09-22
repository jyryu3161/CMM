# Agent: succinate in anaerobic *E. coli*

A decision model plays `e_coli_core` to raise succinate export, one move at a time. The
condition is the same as the [SC-01 succinate example](../production-targets/README.md) —
glucose uptake 10 mmol gDW⁻¹ h⁻¹, oxygen closed — so the two runs are directly comparable.

## Run it

```bash
export OPENROUTER_API_KEY=...   # use `read -rs OPENROUTER_API_KEY` to keep it out of history
uv run --frozen --all-extras python examples/production-targets/prepare.py
uv run --frozen --all-extras cmm jev-design --config examples/jev-design/config.json
```

`prepare.py` exports the model COBRApy supplies; this example reuses it rather than shipping
a second copy. Moves print as they are played:

```
R1T8 knockout on D_LACt2 → product 9.381, growth 0.08559
R1T9 JEV ended the round
R2T1 knockout on PFL → product 0.9015, growth 0.1776
R4T1 adopt_best_design on adopt_best_design → product 9.911, growth 0.09065
R4T2 knockdown_50 on ACKr → product 9.946, growth 0.05472
R4T3 JEV ended the round
```

Every round starts again from the wild type — round 2 opens on `PFL` at 0.678, not on round
1's 9.381 — and carries forward only what earlier rounds learned.

## A portfolio, not one design

Each round ends with a diagnosis, printed and written to `02_game/rounds.csv`, and each is
asked a *different* question: one reaction of every design already found is withheld, so a
later round has to reach somewhere else.

```
What each round asked, reached, and left undone:
  round 1 — best design available: 9.381 at growth 0.08559; the agent judged no remaining
    move worth making; the product was still NADPH-limited: one more unit per hour would have
    bought +1.06 more; it finished with 2 unused knockdown(s)
  round 4 — best design without PYK: 9.946 at growth 0.05472; …
  5 distinct design(s) over 6 round(s). The headline is the best of them; the rest are what a
  laboratory that cannot build it would use.
```

Telling the agent that a design had already been found was measured and was not enough: six
rounds gave **two** distinct designs and four exact repeats, because every round starts from
the same wild type and plays the same game. The cut — the device OptKnock uses to enumerate
alternative designs — took the same run to **five** distinct designs with the same headline.

## Every target, for and against

`06_targets/targets.csv` has one row per reaction the run acted on **or merely measured**:

| reaction | genes | measured | for | against |
|---|---|---|---|---|
| `ACKr` | purT, ackA, tdcD | halving it 0 | in the best design; the cell grows without it | needs 3 genes (isozymes), not one; CMM measured halving it changing the product by 0 |
| `SUCCt2_2` | dctA | deleting it −9.9 | — | a single-gene edit, but the same gene runs `FUMt2_2` and `MALt2_2`, which the edit stops too |

They are not scored against each other. "Raises the product by 0.035" and "needs three
isozymes deleted" are not commensurable, and the trade belongs to whoever builds the strain.

Note what `ACKr` reads here. The tick that applied it recorded `+0.03537`, because that was the
change against the design standing at that moment; the targets table re-measures every gain
against the design the run *ended* with, where halving it changes nothing. Both are true of
different designs, which is why each number says which design it belongs to rather than being
averaged into one.

## What it is measured against

Every run scores itself against the deterministic methods on the same problem, every design
applied and solved the same way, and prints the table:

```
  wild type                                guaranteed  0.0000  growth  0.2117
  best single gene deletion                guaranteed  0.0000  growth  0.2076
  OptKnock                                 guaranteed  9.9098  growth  0.0906
  RobustKnock                              guaranteed  9.9098  growth  0.0906
  best deterministic design + one knockdown (exhaustive)
                                           guaranteed  9.9475  growth  0.0527
  best amplification (outside the vocabulary)
                                           guaranteed 10.0310  growth  0.0532
  JEV agent                                guaranteed  9.9457  growth  0.0547  (not deterministic)
```

Every row is the **guaranteed** product: the least that design can make while growing as fast
as it can. A pFBA number is the best case, and a design whose worst case is zero is one the
strain may grow just as fast without ever using.

Read it honestly, in three steps.

*A single gene deletion cannot solve this problem* — the best of 71 reaches nothing once it is
re-optimised. Anaerobic succinate needs several routes closed at once.

*OptKnock proves 9.910 in under a second*, and an agent searching on its own does not match it:
the winning design deletes `LDH_D` and `THD2`, which carry no flux in the wild type, and a board
built from where the flux is today cannot see them. So the run hands the agent that design and
asks it to improve on it — which means the agent's design **contains** OptKnock's, and the
verdict says so.

*What the agent adds is coupling, not ceiling.* Halving acetate kinase on top of OptKnock's
design moves the guarantee from 9.910 to 9.946, and the growth rate from 0.091 to 0.055. Quoted
as "+0.4% product" that is meaningless — the two designs sit at different points of the same
trade-off, and any design can buy product by spending growth. Held at one growth rate the real
difference shows: OptKnock's design is free to make anything from 7.60 to 11.55, and the
agent's is pinned at 9.946. The knockdown does not raise what the strain can make; it removes
the strain's freedom to make less.

*And enumeration finds the same thing.* The `+ one knockdown (exhaustive)` row is the same
proven design plus the best of every `knockdown_50` the agent could have played, tried one at a
time: 9.9475 through `ATPS4r`, in under a second, with no judgement anywhere in it. The agent
landed on `ACKr` at 9.9459. Any claim that the decision model contributes something has to
clear that row, not OptKnock — which is why every run now computes it.

The `best amplification` row is deliberately a move the agent may **not** make. Forcing flux
through a reaction is not what over-expression does to a cell, so it is excluded from the
vocabulary — and the row prices that decision in every run rather than leaving it as an
argument. It is measured on the design being scored, not on the wild type, because on the
wild type FSEOF's top amplification target for succinate buys nothing at all.

This is one run, on one problem, on one small model. The agent's choices are not guaranteed to
repeat — two runs of this configuration have produced designs differing several-fold — so
nothing here is a claim about the method, let alone about your model. Claiming agent performance
needs repeated runs and a stated distribution, which this example does not provide.

## What it costs

Measured on this configuration: two calls per move, about 0.6 s and $0.00016 per move. A
typical run here is 10 to 13 moves, about 7 seconds and $0.002. The limit that matters is
JEV's 32K context, which `candidate_limit` keeps the state inside — the evidence is sent once,
in `state.records`, and the question's criteria carry only a label, because sending the record
in both places cost 40% of the payload and changed nothing.

`enable_web_research` is the exception: a lookup costs roughly **$0.05**, about three hundred
times a decision, because the search results are billed as input. `max_research_calls` bounds
it.

## What comes out

```
results/example-jev-succinate/
  00_summary.json        the headline and the best design
  02_game/ticks.csv      every move, why it was chosen, what it did
  02_game/candidate_rankings.csv   the agent's probability over every option, every tick
  03_design/best_design.csv        the design that scored best while holding the growth floor
  04_agent/transcript.jsonl        every request and response
  04_agent/literature.csv          what the web lookup returned, in full, with its sources
  05_baseline/comparison.csv       what the deterministic methods give on the same problem
```

## Read the result honestly

**The agent's choices are not guaranteed to repeat.** The CMM solves behind them are
deterministic; the decisions are not. Repeated runs of this exact config have produced designs
differing several-fold in succinate flux. The transcript is what makes a single run auditable;
it is not what makes the method reproducible. Do not quote one run as the method's
performance, and do not compare a single JEV run against a deterministic method as though the
two were the same kind of measurement.

A run that fails to beat the wild type is a real result. Report it as one.

`cmm report validate` does not accept a JEV run: there is no publication renderer or
completion gate for this workflow, and validating it against the production contract would
report a long list of artifacts it was never meant to contain.

See [SC-03](../../docs/scenarios/SC-03-jev-agent-design.md) for what each number means, how
the board of candidate reactions is built, and what the run does not establish.

## Telling the agent what you know

`brief` in the config — and the text box on the *JEV Agent* tab — is for what the model does
not contain:

```
- The published targets for succinate in anaerobic E. coli are ldhA, pflB, ptsG and pta/ackA.
- Growth has to stay above 0.1 per hour for this strain to be useful.
- NADPH supply is the cofactor I expect to be limiting.
- Leave the pentose phosphate pathway alone; we cannot engineer it here.
```

It is guidance, not permission. It cannot widen the move vocabulary, name a reaction the model
does not have, or lift the growth floor CMM enforces — the agent still answers only with the
criteria CMM supplies — so the worst a mistaken brief can do is spend steps. It is recorded in
`00_provenance.json`, because a result that was steered deserves to say so.

## Changing it

Copy this directory and edit `config.json`. The product must be an exchange reaction that
exists in the model, and the condition must state the medium, substrate uptake and aeration
explicitly — the agent optimises against whatever condition it is given, so an unstated one
is an unstated assumption in the result.

**Another organism works.** The cofactor pools are found by chemical formula rather than by
id, so a model using Yeast-GEM or AGORA naming is read correctly rather than silently
mis-read, and a pool that cannot be identified is named on the screen instead of omitted.
`organism` is required as soon as `enable_web_research` is on, and has no default: a
literature lookup about the wrong species returns an answer that is confident and wrong.

`enable_web_research: true` adds **one** literature lookup, through OpenRouter's web plugin,
before the first move and never inside a step. What it returns is shown in every state under
`published_evidence` and kept whole in `04_agent/literature.csv` with its sources. It used to
be a lookup per candidate reaction, run between the two calls of a step: correct, and unusable
— a web search takes tens of seconds, the loop waits on it with nothing to do, and eight of
them turned a run that plays in a minute into one that takes ten. It is data, never an
instruction: JEV can still only answer with the moves this package defines.
