# ADR-0004 — The publication report proposes no strain

**Status:** Accepted · recorded 2026-09-16 (decision predates this record)

## Context

An SC-01 run produces several independent rankings: FSEOF and FVSEOF amplification targets,
OptKnock and RobustKnock designs, MOMA and ROOM single-knockout screens, plus flux-response and
sampling verification for each candidate. `07_validation/recommendations.csv` collects the
candidates that passed verification.

The obvious next step — merging those into one ranked list of "recommended targets" or a single
proposed strain — is what a reader expects a report to do. It is also the step CMM cannot
justify: there is no validated method for ranking a FSEOF target against an OptKnock design, and
the source papers do not provide one. A merged list would invent an ordering and present it with
the authority of the numbers underneath it.

## Decision

Keep `recommendations.csv` as a **machine-readable validation artifact**. The canonical
publication report must **not** synthesize it into recommended targets, a strain proposal, a
summary promotion, or a figure category.

Present each method's results separately and leave intervention selection to the user.

FSEOF and FVSEOF keep independent top-10 rankings. Membership in both is useful provenance, not
a prerequisite for validation or promotion.

## Consequences

- The report answers "what did each method find, and did it survive verification?" — not "what
  should I build?". That is a deliberate limit on what the document claims.
- Readers wanting a single recommendation must apply their own criteria, which they can, because
  every method's ranking and its verification status are present.
- A downstream study may synthesize across methods; it then owns that claim rather than
  inheriting it from CMM.
- CMM's tests establish implementation correctness, not biological validity. Every predicted
  target is a hypothesis to test experimentally, which is a second reason not to present one as
  a recommendation.

## Source

[`AGENTS.md`](../../AGENTS.md) rule 11 and §5; [`docs/scenarios/_reporting.md`](../scenarios/_reporting.md).
