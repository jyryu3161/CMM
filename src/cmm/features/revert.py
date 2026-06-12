"""Revert-metabolism: normalization-target prediction via robust MTA (rMTA).

Predicts gene/reaction knockouts that move a perturbed (e.g. disease) metabolic state back
toward a reference (e.g. healthy) state. This is the inverse of the production-design
services: instead of pushing flux toward a product, it ranks interventions by how well they
transform the source flux distribution in the direction implied by two-state expression.

Method (robust MTA, in the spirit of Valcarcel 2019, continuous form):
  For each knockout, solve a MOMA-based QP three times with different transformation weights
  - best  (alpha): reward moving target reactions in the desired direction
  - moma  (0):      pure minimal adjustment, neutral
  - worst (-alpha): reward moving target reactions the WRONG way (adversarial)
  Each solution gets a transformation score TS = (correct - wrong movement of target
  reactions) / (1 + steady-reaction disturbance). The robust score rTS combines the three so
  a knockout ranks high only when it transforms well in the best case and is not harmful in
  the worst case.

  NOTE — deliberate deviations from the published rMTA, documented and validated in
  ``docs/design-revert-metabolism.md``: the robust score uses the additive worst-case form
  ``mTS + wTS`` (gated on bTS/wTS) rather than the textbook multiplicative gate
  ``mTS * (bTS - wTS)``, and the TS denominator is the L2 steady disturbance with a ``+1``
  regularizer rather than the published L1 steady deviation. The textbook gate (which fires
  only when wTS < 0) zeroes out knockouts that *force* a correct reroute; the additive form
  is the one validated to rank such forced-reversion targets correctly, so it is the default.
  Treat the scores as a self-consistent ranking heuristic, not numerically identical to the
  cobratoolbox rMTA.m output.

The continuous rMTA runs on any QP solver. An optional original-MTA MIQP mode is gated on a
MIQP-capable solver.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from cobra import Model

from cmm.core import solvers
from cmm.core.condition import Condition
from cmm.core.flux_state import FluxState
from cmm.core.results import TargetRanking
from cmm.features._perturbation import (
    Perturbation,
    apply_perturbation,
    gene_perturbations,
    reaction_perturbations,
)
from cmm.omics.differential import DirectionMap

RevertMethod = Literal["rmta", "mta", "mta_miqp"]


@dataclass(frozen=True)
class _Scores:
    moma: float
    best: float
    worst: float
    robust: float
    status: str


def _sign(x: float, tol: float = 1e-9) -> int:
    if x > tol:
        return 1
    if x < -tol:
        return -1
    return 0


def _target_reactions(model: Model, reference: FluxState) -> list[str]:
    present = {r.id for r in model.reactions}
    return [rid for rid in reference.reactions() if rid in present]


def _mta_qp(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    alpha: float,
) -> tuple[dict[str, float], str]:
    """Solve the MOMA-based transformation QP and return (fluxes, status).

    Objective (minimize): sum_i (v_i - v_ref_i)^2  -  alpha * sum_{nonsteady} d_i (v_i - v_ref_i)
    The quadratic term is full MOMA (keeps the solve bounded/unique); the linear term rewards
    moving target reactions in their desired direction (alpha<0 rewards the wrong way).
    """

    prob = model.problem
    saved_objective = model.objective
    saved_direction = model.objective_direction
    target_set = set(target_rxns)
    try:
        expr = 0.0
        # Full MOMA quadratic over EVERY model reaction (reactions absent from the reference
        # contribute (v_i - 0)^2), matching comparison.py's reference solution and keeping the
        # QP bounded and unique. The linear transformation reward applies only to directional
        # reactions that are part of the reference target set.
        for rxn in model.reactions:
            rid = rxn.id
            vref = reference.get(rid)
            diff = rxn.flux_expression - vref
            expr = expr + diff * diff
            if alpha != 0 and rid in target_set:
                d = direction.get(rid, 0)
                if d != 0:
                    expr = expr - alpha * d * diff
        model.objective = prob.Objective(expr, direction="min")
        solution = model.optimize()
        fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
        return fluxes, solution.status
    finally:
        model.objective = saved_objective
        model.objective_direction = saved_direction


def _mta_miqp(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    weight: float,
    eps: float,
) -> tuple[dict[str, float], str]:
    """Original-MTA MIQP solve: maximize the count of correctly transformed reactions.

    Binary z_i is 1 only when nonsteady reaction i moves past its reference by ``eps`` in the
    desired direction (big-M indicator). Minimizes ``weight * steady_distance - sum z_i``,
    i.e. maximizes successful transformations while keeping steady reactions near reference.
    Requires a MIQP-capable solver (quadratic objective + binary variables).
    """

    from optlang.symbolics import Zero

    # Work on a copy: isolates the binary vars from the caller's model/history and avoids
    # cross-knockout name collisions. The copy captures the current (KO-applied) bounds.
    work = model.copy()
    prob = work.problem

    # Phase 1: add the binary variables and FLUSH, so the constraints that reference them in
    # phase 2 do not make optlang re-queue the same variable (the duplicate-add trap).
    switches: list = []
    specs: list[tuple[str, int, float, float]] = []
    for rid in target_rxns:
        d = direction.get(rid, 0)
        if d == 0:
            continue
        rxn = work.reactions.get_by_id(rid)
        vref = reference.get(rid)
        big_m = abs(rxn.lower_bound) + abs(rxn.upper_bound) + abs(vref) + 1.0
        z = prob.Variable(f"mta_z_{rid}", type="binary")
        switches.append(z)
        specs.append((rid, d, vref, big_m))
    if switches:
        work.add_cons_vars(switches)
        work.solver.update()

    # Phase 2: constraints referencing the now-existing binaries.
    constraints = []
    for (rid, d, vref, big_m), z in zip(specs, switches, strict=True):
        flux = work.reactions.get_by_id(rid).flux_expression
        if d > 0:
            c = prob.Constraint(flux - big_m * z, lb=vref + eps - big_m, name=f"mta_up_{rid}")
        else:
            c = prob.Constraint(flux + big_m * z, ub=vref - eps + big_m, name=f"mta_dn_{rid}")
        constraints.append(c)
    if constraints:
        work.add_cons_vars(constraints)
        work.solver.update()

    # Objective: minimize weighted steady distance minus the count of transformed reactions.
    quad = Zero
    for rid in target_rxns:
        if direction.get(rid, 0) == 0:
            diff = work.reactions.get_by_id(rid).flux_expression - reference.get(rid)
            quad = quad + weight * diff * diff
    work.objective = prob.Objective(quad - sum(switches) if switches else quad, direction="min")

    solution = work.optimize()
    fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
    return fluxes, solution.status


def _transformation_score(
    fluxes: dict[str, float],
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
) -> float:
    """Net correct movement of target reactions, normalized by steady disturbance.

    TS = (sum |delta| over correctly-moved targets - sum |delta| over wrongly-moved targets)
         / (1 + sum delta^2 over steady reactions).
    The L1 numerator matches the published transformation score; the denominator is the L2
    steady disturbance plus a +1 regularizer (a deliberate, documented choice — see the module
    docstring — that keeps the score finite and bounded rather than the published L1 form).
    """

    correct = 0.0
    wrong = 0.0
    steady_dev = 0.0
    for rid in target_rxns:
        vref = reference.get(rid)
        delta = fluxes.get(rid, 0.0) - vref
        d = direction.get(rid, 0)
        if d == 0:
            steady_dev += delta * delta
        else:
            moved = _sign(delta)
            if moved == d:
                correct += abs(delta)
            elif moved == -d:
                wrong += abs(delta)
    return (correct - wrong) / (1.0 + steady_dev)


def _score_knockout(
    model: Model,
    reference: FluxState,
    direction: DirectionMap,
    target_rxns: list[str],
    method: RevertMethod,
    alpha: float,
) -> _Scores:
    if method == "mta_miqp":
        flux, status = _mta_miqp(model, reference, direction, target_rxns, weight=1.0, eps=1e-3)
        if status != "optimal":
            return _Scores(moma=0.0, best=0.0, worst=0.0, robust=0.0, status=status)
        ts = _transformation_score(flux, reference, direction, target_rxns)
        return _Scores(moma=ts, best=ts, worst=ts, robust=ts, status="optimal")

    moma_flux, status = _mta_qp(model, reference, direction, target_rxns, 0.0)
    if status != "optimal":
        return _Scores(moma=0.0, best=0.0, worst=0.0, robust=0.0, status=status)
    moma_ts = _transformation_score(moma_flux, reference, direction, target_rxns)
    if method == "mta":
        best_flux, b_status = _mta_qp(model, reference, direction, target_rxns, alpha)
        best_ts = (
            _transformation_score(best_flux, reference, direction, target_rxns)
            if b_status == "optimal"
            else moma_ts
        )
        return _Scores(moma=moma_ts, best=best_ts, worst=moma_ts, robust=best_ts, status="optimal")

    best_flux, b_status = _mta_qp(model, reference, direction, target_rxns, alpha)
    worst_flux, w_status = _mta_qp(model, reference, direction, target_rxns, -alpha)
    best_ts = (
        _transformation_score(best_flux, reference, direction, target_rxns)
        if b_status == "optimal"
        else moma_ts
    )
    worst_ts = (
        _transformation_score(worst_flux, reference, direction, target_rxns)
        if w_status == "optimal"
        else moma_ts
    )
    robust = _robust_score(moma_ts, best_ts, worst_ts)
    return _Scores(moma=moma_ts, best=best_ts, worst=worst_ts, robust=robust, status="optimal")


def _robust_score(moma_ts: float, best_ts: float, worst_ts: float) -> float:
    """Combine the three transformation scores into a robust score.

    A knockout is a robust normalization target when it transforms toward the target in the
    best case (best_ts > 0) and remains beneficial — not adversarially reversible — in the
    worst case (worst_ts >= 0). The robust score rewards the guaranteed (worst-case)
    transformation, anchored by the neutral MOMA score.
    """

    if best_ts > 0 and worst_ts >= 0:
        return moma_ts + worst_ts
    if best_ts > 0:
        # Helps in the best case but is adversarially reversible: keep only a damped signal.
        return max(0.0, moma_ts + worst_ts)
    return 0.0


def revert_targets(
    model: Model,
    source_condition: Condition | None,
    reference_state: FluxState,
    direction: DirectionMap,
    *,
    targets: Iterable[str] | None = None,
    method: RevertMethod = "rmta",
    alpha: float = 0.9,
    perturbation: Literal["gene", "reaction"] = "gene",
) -> TargetRanking:
    """Rank knockouts by how well they revert the source state toward the target.

    Parameters mirror the design doc: ``reference_state`` is the source/disease reference
    flux distribution, ``direction`` is the per-reaction desired flux direction from
    two-state differential expression, and ``alpha`` is the transformation weight.
    """

    if method == "mta_miqp":
        solvers.require("MIQP", model.solver.interface, feature="original MTA (MIQP)")
    else:
        solvers.require("QP", model.solver.interface, feature="rMTA revert-metabolism")

    with model:
        if source_condition is not None:
            source_condition.apply_to(model)

        target_rxns = _target_reactions(model, reference_state)

        if perturbation == "gene":
            perts: list[Perturbation] = gene_perturbations(model, targets)
        else:
            perts = reaction_perturbations(model, targets)

        scored: list[tuple[str, float]] = []
        for pert in perts:
            with apply_perturbation(model, pert):
                scores = _score_knockout(
                    model, reference_state, direction, target_rxns, method, alpha
                )
            scored.append((pert.target_id, scores.robust))

    return TargetRanking.from_scores(
        method=f"revert_{method}",
        scores=scored,
        metadata={
            "alpha": alpha,
            "perturbation": perturbation,
            "reference": reference_state.name,
            "reference_provenance": reference_state.provenance,
            "n_targets": len(scored),
        },
    )
