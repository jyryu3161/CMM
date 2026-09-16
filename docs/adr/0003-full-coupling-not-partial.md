# ADR-0003 — Coupled sets come from the null space, not partial coupling

**Status:** Accepted · recorded 2026-09-16 (decision predates this record)

## Context

A transformation run deduplicates reactions that cannot vary independently, so that one
biological intervention is counted once rather than several times. Counting coupled reactions
separately inflates the denominator of every percentile claim.

The source paper groups reactions by **partial** coupling, which requires O(*n*²) linear
programmes — prohibitive on a genome-scale reconstruction inside an interactive workflow.

## Decision

Compute **full** coupling from the null space of the stoichiometric matrix S.

The null space always uses the **full model**, even when only a subset of reactions is eligible
for deletion. Filtering knockout candidates does not remove those reactions from the network
that determines coupling.

`collapse_coupled_sets=None` follows the perturbation level: on for reactions, off for genes,
because coupled sets are defined on reactions. Asking for them on a gene-level run is rejected
rather than silently ignored. A gene-level run instead deduplicates genes that block the same
reaction signature.

## Consequences

- Full coupling is a **stronger** condition than partial coupling, so the grouping is
  conservative: it can split one of the paper's sets but never merge two. The candidate count it
  yields is an **upper bound** on the paper's.
- Percentile claims computed against that count are therefore conservative in the same
  direction, which is the safe direction for a claim.
- CMM's candidate counts will not match the paper's exactly, and a comparison must say why
  rather than treating the difference as a bug.
- This is a tractability trade-off. If O(*n*²) LPs become affordable, reopen the record.

## Source

[`docs/scenarios/SC-02-transformation-target-discovery.md`](../scenarios/SC-02-transformation-target-discovery.md)
step 4; [`docs/agent-reference.md`](../agent-reference.md), transformation section.
