"""Flux coupling — grouping reactions whose deletion has the same consequence.

A ranking over knockout candidates is only meaningful if each candidate is a distinct
intervention. Three reactions of an unbranched pathway always carry the same flux, so deleting
any one of them stops all three; counting them as three candidates counts one intervention three
times and inflates the denominator of any "top *N*%" claim.

Yizhak et al. (2013) reduce their candidate set this way before their reaction-level validation
analyses — dead-end reactions out, essential reactions out, then *"the set of simulated
knockouts is composed of a member from each partially coupled set (including singleton sets)"*.
Valcárcel et al. (2019) rank gene knockouts and apply no coupling reduction.

**This module computes full coupling.** Two reactions are fully coupled when ``v_i / v_j`` is
the same constant in every steady state; partially coupled when each is non-zero exactly when
the other is, with a ratio free to vary. Full coupling is strictly stronger, so the grouping is
conservative: it can split a partially coupled set but can never merge two, so the candidate
count it yields is an upper bound on the partial one. Full coupling also falls out of one
null-space computation, whereas partial coupling needs flux coupling analysis — O(n^2) linear
programmes, impractical at genome scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence

import numpy as np
from cobra import Model
from scipy.linalg import null_space

#: Two null-space rows count as parallel when their normalised difference falls below this.
#: Also the magnitude under which a row is treated as carrying no steady-state flux at all.
PARALLEL_TOLERANCE = 1e-6


@dataclass(frozen=True)
class CoupledSets:
    """Reactions grouped by full flux coupling, with one representative per group."""

    #: reaction id -> group index. Reactions with no steady-state flux each get their own.
    groups: Mapping[str, int]
    #: One reaction id per group, chosen as the alphabetically first so reruns agree.
    representatives: tuple[str, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.representatives)

    def members(self, reaction_id: str) -> tuple[str, ...]:
        """Every reaction coupled to ``reaction_id``, itself included."""

        group = self.groups[reaction_id]
        return tuple(sorted(rid for rid, gid in self.groups.items() if gid == group))

    def representative_for(self, reaction_id: str) -> str:
        """The representative standing in for ``reaction_id``'s group."""

        group = self.groups[reaction_id]
        for rid in self.representatives:
            if self.groups[rid] == group:
                return rid
        raise KeyError(reaction_id)

    def to_provenance(self) -> dict[str, object]:
        return {
            "coupling": "full",
            "n_reactions": len(self.groups),
            "n_sets": len(self.representatives),
            **dict(self.metadata),
        }


def coupled_reaction_sets(
    model: Model,
    reactions: Iterable[str] | None = None,
    *,
    tolerance: float = PARALLEL_TOLERANCE,
) -> CoupledSets:
    """Group ``reactions`` into fully coupled sets from the null space of S.

    Every steady state satisfies ``Sv = 0``, so ``v = N c`` for ``N = null_space(S)`` and some
    coefficient vector ``c``. Reaction *i*'s flux is then the inner product of ``N``'s *i*-th
    row with ``c``. If two rows are parallel, the ratio of their fluxes is that constant of
    proportionality in *every* steady state — which is full coupling.

    Rows are normalised and sign-canonicalised before hashing, so ``v_i = -3 v_j`` groups with
    ``v_i = 3 v_j``: the deletion consequence is the same either way.

    ``reactions`` selects the candidates to return, not the network to solve. Coupling is
    always computed from the full model: reactions excluded as knockout candidates still
    carry flux and must remain in S.
    """

    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    rxn_ids: Sequence[str] = (
        list(reactions) if reactions is not None else [r.id for r in model.reactions]
    )
    if not rxn_ids:
        return CoupledSets({}, (), {"tolerance": tolerance})
    missing = [rid for rid in rxn_ids if rid not in model.reactions]
    if missing:
        raise KeyError(f"reactions absent from the model: {sorted(missing)[:5]}")

    index = {reaction.id: i for i, reaction in enumerate(model.reactions)}
    met_index = {met.id: i for i, met in enumerate(model.metabolites)}
    stoich = np.zeros((len(met_index), len(model.reactions)))
    for reaction in model.reactions:
        for met, coeff in reaction.metabolites.items():
            stoich[met_index[met.id], index[reaction.id]] = coeff

    basis = null_space(stoich)
    norms = (
        np.linalg.norm(basis, axis=1) if basis.size else np.zeros(len(model.reactions))
    )

    groups: dict[str, int] = {}
    signatures: dict[tuple, int] = {}
    n_no_flux = 0
    next_singleton = -1
    for rid in rxn_ids:
        row = index[rid]
        if basis.size == 0 or norms[row] < tolerance:
            # No steady state carries flux through it, so it is coupled to nothing. Give it a
            # distinct negative id rather than pooling every such reaction into one group.
            groups[rid] = next_singleton
            next_singleton -= 1
            n_no_flux += 1
            continue
        unit = basis[row] / norms[row]
        lead = next((x for x in unit if abs(x) > tolerance), 1.0)
        if lead < 0:
            unit = -unit
        key = tuple(np.round(unit / tolerance).astype(np.int64))
        groups[rid] = signatures.setdefault(key, len(signatures))

    representatives: list[str] = []
    seen: set[int] = set()
    for rid in sorted(rxn_ids):
        group = groups[rid]
        if group < 0:
            representatives.append(rid)
            continue
        if group in seen:
            continue
        seen.add(group)
        representatives.append(rid)

    sizes = [
        sum(1 for gid in groups.values() if gid == g) for g in set(groups.values())
    ]
    return CoupledSets(
        groups=groups,
        representatives=tuple(representatives),
        metadata={
            "tolerance": tolerance,
            "n_reactions_without_steady_state_flux": n_no_flux,
            "largest_set": max(sizes) if sizes else 0,
            "n_singleton_sets": sum(1 for size in sizes if size == 1),
        },
    )
