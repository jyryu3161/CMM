# JEV agent: succinate in anaerobic *E. coli*

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
R2T1 knockout on PFL → product 0.6781, growth 0.18
R2T2 restore_best_design on restore_best_design → product 9.381, growth 0.08559
R2T3 knockdown_50 on PFL → product 9.423, growth 0.07983
R3T1 adopt_best_design on adopt_best_design → product 9.911, growth 0.09065
R3T2 knockdown_50 on ACKr → product 9.946, growth 0.05472
```

Every round starts again from the wild type — round 2 opens on `PFL` at 0.678, not on round
1's 9.381 — and carries forward only the record of what earlier rounds reached, which is why
`restore_best_design` exists.

## What it is measured against

Every run scores itself against the deterministic methods on the same problem, every design
applied and solved the same way, and prints the table:

```
  wild type                              product    0.0000  growth  0.2117
  best single gene deletion (MOMA-L2)    product    0.2114  growth  0.1646
  OptKnock                               product    9.9108  growth  0.0906
  RobustKnock                            product    9.9108  growth  0.0906
  best amplification on top of this design (outside the vocabulary)
                                         product   10.0387  growth  0.0532
  JEV agent                              product    9.9461  growth  0.0547  (not deterministic)
```

Read it honestly. A single gene deletion cannot solve this problem — the best of 71 reaches
0.211. OptKnock proves 9.911 in under a second, and an agent searching on its own does not
match it: the winning design deletes `LDH_D` and `THD2`, which carry no flux in the wild type,
and a board built from where the flux is today cannot see them. What the agent adds is the
move OptKnock's formulation cannot express — a **partial** knockdown, its own variables being
present-or-absent. Halving acetate kinase on top of OptKnock's design takes 9.911 to 9.946,
at a cost in growth (0.055 against 0.091). That +0.4% is the whole margin, and it is an
honest one: an exhaustive screen of every deletion and every 50% knockdown on top of that
design reaches the same 9.9461, so the agent found the best its vocabulary allows.

The `best amplification` row is deliberately a move the agent may **not** make. Forcing flux
through a reaction is not what over-expression does to a cell, so it is excluded from the
vocabulary — and the row prices that decision in every run rather than leaving it as an
argument. It is measured on the design being scored, not on the wild type, because on the
wild type FSEOF's top amplification target for succinate buys nothing at all.

Over eight runs of this configuration: all eight reached 10.761, in six steps and five
seconds each. That is one problem on one small model, and it is not a claim about yours.

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

`enable_web_research: true` adds a literature lookup for each candidate through OpenRouter's
web plugin. What it returns is pasted into that candidate's record as evidence for the agent
to weigh. It is data, never an instruction: JEV can still only answer with the moves this
package defines.
