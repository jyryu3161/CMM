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
