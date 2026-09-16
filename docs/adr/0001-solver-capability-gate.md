# ADR-0001 — A missing solver capability fails the run

**Status:** Accepted · recorded 2026-09-16 (decision predates this record)

## Context

The cobra default solver, GLPK, does LP and MILP only. Several shipped methods need more:
MOMA-L2, E-Flux2, and `rmta_continuous` need QP; published MTA and rMTA need MIQP;
OptKnock/RobustKnock need MILP plus an importable `straindesign`.

Every one of these has a cheaper neighbour that GLPK *can* run — MOMA-L1 for MOMA-L2, LAD for
E-Flux2, `rmta_continuous` for rMTA. Substituting one for another is a one-line change and
produces numbers that look entirely normal. Nothing in the output distinguishes a result from
the requested method from a result from its substitute.

## Decision

Check capability **before** the call, via `cmm.core.solver_status` / `cmm.core.supports`.

For a narrow request, an LP-capable substitute is permitted and **must be named in the report**
along with the reason. Silence is not an option in either direction: never produce nothing, and
never downgrade without saying so.

For the canonical SC-01 workflow the bar is higher. Its single-knockout comparison requires
MOMA-L2 *and* ROOM so the two requested methods stay comparable, so it **fails its capability
gate** rather than replacing MOMA-L2 with MOMA-L1.

## Consequences

- A user without gurobi or cplex sees an explicit failure naming the missing capability, not a
  quietly different analysis.
- CMM cannot claim "works with the default solver" for the full method set, and the install
  documentation has to carry that caveat.
- Test suites mark capability-gated cases (`requires_qp`, `requires_miqp`) rather than skipping
  inline, so a machine with a commercial solver still exercises them.
- `rmta_continuous` is a QP heuristic and is **not** published rMTA. It must never stand in for
  a method the solver cannot run; report the method as unavailable instead.

## Source

[`AGENTS.md`](../../AGENTS.md) §2 and rule 4; the method table in
[`docs/VALIDATION.md`](../VALIDATION.md).
