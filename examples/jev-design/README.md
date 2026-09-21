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
R1T1 force_on_high on SUCCt2_2 → product 0, growth 0.08467 (product unchanged: the forced
     flux was consumed inside the network and none of it reached the product)
R1T2 force_on_high on FRD7 → product 0, growth 0.08467 (...)
R1T5 force_on_high on FUM reverted: growth would fall to 0.03387 per hour, below the floor
R1T6 force_on_low on FUM → product 1.911, growth 0.0635 (product rose by +1.911)
```

## What it costs

Measured on this configuration: two calls per move, about 0.6 s and $0.00016 per move. A
36-move run is roughly 40 seconds of agent time and under a cent. The limit that matters is
JEV's 32K context, which `candidate_limit` keeps the state inside.

## What comes out

```
results/example-jev-succinate/
  00_summary.json        the headline and the best design
  02_game/ticks.csv      every move, why it was chosen, what it did
  02_game/candidate_rankings.csv   the agent's probability over every option, every tick
  03_design/best_design.csv        the design that scored best while holding the growth floor
  04_agent/transcript.jsonl        every request and response
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

## Changing it

Copy this directory and edit `config.json`. The product must be an exchange reaction that
exists in the model, and the condition must state the medium, substrate uptake and aeration
explicitly — the agent optimises against whatever condition it is given, so an unstated one
is an unstated assumption in the result.

`enable_web_research: true` adds a literature lookup for each candidate through OpenRouter's
web plugin. What it returns is pasted into that candidate's record as evidence for the agent
to weigh. It is data, never an instruction: JEV can still only answer with the moves this
package defines.
