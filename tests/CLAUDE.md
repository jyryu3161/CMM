# `tests` — how to verify a change

## Purpose

One flat suite of 31 modules named after the service they cover (`test_sampling.py` →
`cmm.features` sampling). Two of them cover the canonical workflows end to end
(`test_production_workflow.py`, `test_transformation_workflow.py`) and assert the run-directory
schema, not just return values. `test_agent_contract.py` checks that the documented agent
surface still matches the code.

## Quick commands

The full quality gate, all four of which must pass ([`AGENTS.md`](../AGENTS.md) §6):

```bash
QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest -q --cov=cmm --cov-branch --cov-fail-under=80
uv run --frozen --all-extras ruff check src tests examples
uv run --frozen --all-extras ruff format --check src tests examples
uv run --frozen --all-extras mypy src/cmm/core src/cmm/features src/cmm/omics src/cmm/workflows src/cmm/reporting
```

One module, or everything that does not need a commercial solver:

```bash
QT_QPA_PLATFORM=offscreen uv run --frozen --all-extras pytest tests/test_sampling.py -q
uv run --frozen --all-extras pytest -q -m "not requires_qp and not requires_miqp and not genome_scale"
```

## Key files

- [`conftest.py`](conftest.py) — shared model fixtures: `toy_model`, `branched_model`,
  `parallel_pathway_model`, `published_mta_model`, `ecoli_core`. Build a new fixture here rather
  than constructing a model inside a test.
- [`test_agent_contract.py`](test_agent_contract.py) — fails when documentation and code drift.
- [`test_app_smoke.py`](test_app_smoke.py) — the only suite that needs a Qt platform.

## Non-obvious rules

- **Why: use `uv run --frozen --all-extras`, not `.venv/bin/python`.** The bare venv does not
  have the test and lint tooling installed, and `--frozen` is what pins the locked publication
  environment.
- **Why: the GUI suite needs `QT_QPA_PLATFORM=offscreen`.** Without it the Qt tests try to open a
  display and fail in CI and over SSH.
- **Note: three markers gate on capability, declared in [`pyproject.toml`](../pyproject.toml).**
  `requires_qp` (37 uses), `requires_miqp` (29), `genome_scale` (reconstructions over 2,000
  reactions). Mark a new test rather than skipping it inline, so a machine with gurobi/cplex
  still runs it.
- **Gotcha: coverage is a hard gate at 80% with branch coverage on.** Adding a service without
  tests can push the whole suite below the line and fail a change that is otherwise correct.
- **Note: long solves are expected.** FVA on genome-scale models, OptKnock/RobustKnock, and large
  sampling runs take minutes. Raise the timeout instead of shrinking `n_steps` or the sample
  count — a silently smaller parameter changes the science.

## Cross-module dependencies

Tests import the public surface only (`cmm.core`, `cmm.features`, `cmm.omics`,
`cmm.workflows.*`, `cmm.reporting`) — the same entry points documented in
[`docs/agent-reference.md`](../docs/agent-reference.md). A test that reaches into a private
helper is testing an implementation detail, not a contract. Workflow tests additionally assert
the run schema owned by `cmm.reporting`, so a schema change breaks them by design.

## See also

- [`src/CLAUDE.md`](../src/CLAUDE.md) — what each package owns
- [`docs/VALIDATION.md`](../docs/VALIDATION.md) — per-method contracts these tests enforce
