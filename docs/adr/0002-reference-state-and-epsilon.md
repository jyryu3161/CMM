# ADR-0002 — `v_ref` is an E-Flux2/LAD state, and ε has no default

**Status:** Accepted · recorded 2026-09-16 (decision predates this record)

## Context

MTA and rMTA measure everything against a reference flux state `v_ref`: the sign flip that puts
expression labels into flux-value space, the MIQP's success thresholds `v_ref ± ε`, and the
denominator of the transformation score.

The source papers derive `v_ref` from iMAT followed by sampling. CMM does not implement iMAT.
It derives `v_ref` from the source expression through E-Flux2 or LAD instead.

ε is a flux magnitude. Its correct value depends on the model, the medium, and the units of the
data at hand. There is no value that is safe across models, and the papers derive theirs for
their own setting.

## Decision

Use E-Flux2 (or LAD, when QP is unavailable) for `v_ref`, and **state in every run's provenance
that this is not the published iMAT-plus-sampling state**.

Do not ship a default ε. `suggest_epsilon(reference_fluxes)` returns percentiles of `|v_ref|` so
a value can be chosen against the model at hand, and the chosen value is recorded. Where a run
configures it, report how the ranking moves across an ε grid — that sensitivity analysis is the
honest substitute for the papers' derivation, and both source papers report one.

## Consequences

- A CMM rMTA ranking is **not** a reproduction of the published pipeline, and no report may
  imply that it is. The departure is disclosed, not buried.
- Two runs with different ε are not comparable, which is why ε belongs in provenance rather than
  in a constant.
- `result.summary()["candidate_construction"]` carries the count that is the denominator of any
  percentile claim, so a "top 5%" statement can be checked.
- Implementing iMAT-plus-sampling would reopen this record.

## Source

[`docs/agent-reference.md`](../agent-reference.md), transformation section;
[`docs/scenarios/SC-02-transformation-target-discovery.md`](../scenarios/SC-02-transformation-target-discovery.md)
steps 2 and 5; [`docs/VALIDATION.md`](../VALIDATION.md).
