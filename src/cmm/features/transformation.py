"""Find knockouts that transform one flux state (condition A) toward another (condition B).

Given a source flux distribution (e.g. condition A, possibly omics-predicted) and a target
distribution (condition B), rank candidate gene/reaction knockouts by how close they bring
the flux to the target. Two methods, **both of them published**:

- ``moma``: for each knockout, the cell minimally adjusts from A (MOMA); rank by how much that
  adjusted state moves toward B (reduction in distance-to-B). This is the MOMA-based
  transforming-perturbation baseline of Yizhak K, Gabay O, Cohen H & Ruppin E (2013),
  "Model-based identification of drug targets that revert disrupted metabolism and its
  application to ageing", *Nat Commun* 4:2632, Methods.
- ``mta``: derive a per-reaction desired direction from B − A and run the MTA MIQP of the same
  paper, via :func:`~cmm.features.revert.revert_targets`.

Neither path is CMM's invention; what is CMM's is deriving the direction set from two explicit
flux states rather than from expression through a GPR, so omics-predicted condition fluxes can
be used directly as source/target.

**The default is ``method="moma"``, and Yizhak et al. report that path as "markedly inferior"
to MTA.** It is kept as the default anyway, for two measured reasons, and labelled here rather
than silently: it needs only QP instead of MIQP and runs about 8× faster per gene at genome
scale (2.4 s against 19.3 s on *i*JO1366, 55 min against 7.3 h for a full gene scan); and on
CMM's own SC-02 condition pair the ``mta`` ranking is the more degenerate of the two — 11
distinct scores against MOMA's 21 before the score fix in :mod:`cmm.features.revert`. The two
rankings share no gene in their top 20 (Spearman ρ = −0.068 over 71 genes), so the choice is
not a mild bias: report both when the answer matters, and read the ``score_resolution`` and
``largest_tie_block`` entries of the result metadata before trusting a top-k.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Literal

from cobra import Model

from cmm.core import solvers
from cmm.core.flux_state import FluxState
from cmm.core.provenance import run_provenance
from cmm.core.results import TargetRanking, TargetScore
from cmm.features._perturbation import (
    apply_perturbation,
    gene_perturbations,
    perturbation_provenance,
    reaction_perturbations,
    target_display_name,
)
from cmm.features.comparison import moma
from cmm.features.revert import revert_targets, tie_structure
from cmm.omics.differential import DirectionMap

TransformationMethod = Literal["moma", "mta"]


def direction_from_states(
    source: FluxState,
    target: FluxState,
    *,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-6,
) -> DirectionMap:
    """Per-reaction desired flux-value direction to go from source to target (B − A)."""

    keys = (
        set(reactions)
        if reactions is not None
        else set(source.fluxes) | set(target.fluxes)
    )
    directions: dict[str, int] = {}
    for rid in keys:
        delta = target.get(rid) - source.get(rid)
        directions[rid] = 1 if delta > tol else (-1 if delta < -tol else 0)
    return DirectionMap(
        directions=directions,
        metadata={
            "from": source.name,
            "to": target.name,
            # Not the published expression/GPR route: the direction set comes from the
            # difference of two flux vectors, which is CMM's generalization.
            "direction_rule": "flux_state_difference",
            "tol": tol,
        },
    )


def _state_identity(source: FluxState, target: FluxState) -> dict[str, object]:
    """Which two states a transformation ranking was computed between.

    Without this a ``transformation_targets`` provenance record is indistinguishable from a
    :func:`~cmm.features.revert.revert_targets` one — same method family, same formulation
    tag, no trace of the target state that defined the direction.
    """

    return {
        "source": source.name,
        "source_provenance": source.provenance,
        "target": target.name,
        "target_provenance": target.provenance,
    }


def transformation_targets(
    model: Model,
    source_state: FluxState,
    target_state: FluxState,
    *,
    method: TransformationMethod = "moma",
    perturbation: Literal["gene", "reaction"] = "gene",
    targets: Iterable[str] | None = None,
    order: Literal[1, 2] = 2,
    alpha: float = 0.66,
) -> TargetRanking:
    """Rank knockouts that move the source flux state toward the target flux state.

    ``method="moma"`` (the default) is Yizhak et al.'s MOMA baseline; ``method="mta"`` is
    their MTA MIQP. See the module docstring for why the inferior-by-the-paper's-own-account
    method is the default and what to check before acting on either ranking.
    """

    if method not in {"moma", "mta"}:
        raise ValueError(f"unknown method {method!r}; use 'moma' or 'mta'")
    if perturbation not in {"gene", "reaction"}:
        raise ValueError("perturbation must be 'gene' or 'reaction'")
    source_state.validate()
    target_state.validate()

    if method == "mta":
        direction = direction_from_states(source_state, target_state)
        ranking = revert_targets(
            model,
            None,
            source_state,
            direction,
            targets=targets,
            method="mta",
            alpha=alpha,
            perturbation=perturbation,
        )
        return TargetRanking(
            method="transform_mta",
            targets=ranking.targets,
            metadata={
                **ranking.metadata,
                **_state_identity(source_state, target_state),
                # The optimisation is Yizhak's MTA MIQP; the direction set is not derived the
                # published way, so the inherited tag would overstate the correspondence.
                "formulation": "yizhak_2013_mta_with_flux_state_directions",
                "transformation_method": "mta",
            },
        )

    solvers.require("QP", model.solver.interface, feature="MOMA transformation targets")
    compare_rxns = [r.id for r in model.reactions]
    baseline = source_state.distance(target_state, reactions=compare_rxns, order=order)

    provenance = run_provenance(
        model,
        method="transformation_targets",
        transformation_method="moma",
        perturbation=perturbation,
        order=order,
        alpha=alpha,
    )

    if perturbation == "gene":
        perts = gene_perturbations(model, targets)
    else:
        perts = reaction_perturbations(model, targets)

    scored: list[TargetScore] = []
    nonoptimal = 0
    for pert in perts:
        with apply_perturbation(model, pert):
            result = moma(model, source_state, linear=False)
        if result.status != "optimal":
            nonoptimal += 1
            scored.append(
                TargetScore(
                    pert.target_id,
                    float("-inf"),
                    target_name=target_display_name(model, pert.target_id, pert.kind),
                    status=result.status,
                )
            )
            continue
        predicted = FluxState(result.fluxes, name="perturbed")
        distance_to_target = predicted.distance(
            target_state, reactions=compare_rxns, order=order
        )
        scored.append(
            TargetScore(
                pert.target_id,
                baseline - distance_to_target,
                target_name=target_display_name(model, pert.target_id, pert.kind),
                status=result.status,
            )
        )

    return TargetRanking(
        method="transform_moma",
        targets=tuple(scored),
        metadata={
            **provenance,
            **_state_identity(source_state, target_state),
            "formulation": "yizhak_2013_moma_baseline",
            "transformation_method": "moma",
            "perturbation": perturbation,
            "baseline_distance": baseline,
            "distance_order": order,
            "n_targets": len(scored),
            "n_nonoptimal": nonoptimal,
            **perturbation_provenance(perts),
            **tie_structure(scored),
        },
    ).sorted()
