# `evals` — agent task contract checks

## Purpose

This directory owns the evidence that an agent driving CMM produced a run which **honours the
contracts in [`AGENTS.md`](../AGENTS.md)** — not that its biology is right.

Each task in [`tasks/`](tasks) names a representative request and the machine-checkable
properties its output must have: the run validates, provenance carries the fingerprint and
solver, the condition is explicit rather than inherited, infeasible rows survive in the tables,
and the shipped/planned/excluded inventories stay disjoint so that "CMM does not do this" and
"CMM does this" cannot both be true of one capability.

What it deliberately does **not** measure:

- **Biological validity.** A target passing every check here is still a hypothesis to test
  experimentally ([`AGENTS.md`](../AGENTS.md) §5).
- **The conversation.** Only the artifacts an agent produced are inspected, never its reasoning
  or its tool calls.
- **Whether the science is the best available.** `validate_production_run` is a completion gate,
  not a peer review.

## The JEV replicate harness

[`jev_replicates.py`](jev_replicates.py) answers a different question from the task checks
above: not "does this run honour the contracts" but "what does this configuration do when you
run it more than once". A JEV design quoted from one run says nothing about the method, because
the decision model's choices are not guaranteed to repeat.

```bash
uv run --frozen --all-extras python evals/jev_replicates.py \
    examples/jev-design/config.json --runs 10 --target 9.945652
```

`--target` is what an exhaustive search of the same move vocabulary reached, when that is known,
so the report can say how often the agent got there. It reports the **guaranteed product**,
which is what a run ranks designs on. Needs `OPENROUTER_API_KEY`; costs about $0.015 a run on
`e_coli_core` and $0.025 on `iJO1366`.

Measured, ten runs each: the succinate configuration reaches 9.945652 in 10 of 10 with one
distinct design, and the genome-scale D-lactate configuration returns no design in 10 of 10.
Both outcomes are degenerate and opposite; what varies is the path.

## Quick commands

```bash
# check the bundled example runs
uv run --frozen --all-extras python evals/run_evals.py

# check your own run and record the result outside the repository
uv run --frozen --all-extras python evals/run_evals.py \
  --task evals/tasks/sc01_production.json --run-dir results/my-run \
  --json ../cmm-audit/agent-results.json
```

Exit status is `1` when any check fails, so this runs in CI unchanged.

## Reading an agent's own scripts

`run_evals.py` inspects a finished run bundle, which by definition came from CMM. The other
risk is an agent that never produced one: asked for something CMM does not ship, it writes the
method inline and reports numbers that carry no provenance (rule 12). `check_service_use.py`
reads such scripts and reports where they solve, sample, or mutate a model outside CMM:

```bash
uv run --frozen --all-extras python evals/check_service_use.py analysis.py
uv run --frozen --all-extras python evals/check_service_use.py --strict scripts/*.py
```

Two severities. **bypass** is a solve, a cobra analysis call, or a direct solver import — each
message names the CMM function that owns it. **review** is hand-set `objective`, `medium`, or
bounds, which is sometimes legitimate but leaves the change out of the condition record.
Exit status is `1` on any bypass, and `--strict` fails on a review too.

Point it at agent-authored scripts only: CMM's own source calls the solver legitimately.
**A clean report means "no bypass in the files given", never "CMM produced every number"** — an
agent that works in a REPL or deletes its script leaves nothing to read. This is detection
after the fact; nothing here can prevent the bypass.

## How to add a task

A task file declares `id`, `kind` (`production` or `transformation`), the `request` an agent was
given, a default `run_dir`, and an `expect` block:

| Key in `expect` | Checks |
|---|---|
| `validates` | the run passes its own workflow validator |
| `provenance_fields` | each named field exists and is non-empty in `00_provenance.json` |
| `condition_bounds_include` | `00_config.json` pins these exchange bounds explicitly |
| `tables_retain_status` | each CSV has a `status` column; infeasible rows are counted, not dropped |

Add the request you actually care about rather than a synthetic one — a task nobody asks is not
evidence.

## Non-obvious patterns

- **Why: a check that raises is a failure, not a skip.** A contract that cannot be evaluated is
  not a contract that passed, so exceptions are caught and recorded as `FAIL`.
- **Note: results belong outside the repository.** They are run-specific artifacts, and
  `MANIFEST.in` sweeps `docs/**/*.md` into the sdist — a stray result file ships to users.
- **Gotcha: passing every check does not mean the run is good.** These are necessary conditions
  drawn from `AGENTS.md`, not sufficient ones.

## Cross-module dependencies

`run_evals.py` imports only the public surfaces it validates — `cmm.reporting` for
`validate_production_run`, `cmm.reporting.transformation` for `validate_transformation_run`, and
`cmm.features` for the feature-inventory check. It never imports a private helper and never
invokes a solver, so a failing check always points at the run bundle or at those public
contracts, never at this harness re-implementing the science.

## See also

- [`docs/AI-USAGE.md`](../docs/AI-USAGE.md) — why this exists: run provenance alone does not
  establish agent reliability, and a paper relying on an AI-assisted interface needs more
- [`AGENTS.md`](../AGENTS.md) — the contracts these tasks check
- [`docs/VALIDATION.md`](../docs/VALIDATION.md) — per-method scientific contracts
