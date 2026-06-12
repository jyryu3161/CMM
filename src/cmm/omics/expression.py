"""Single-state omics integration: map gene expression to fluxes (E-Flux2, LAD).

These methods predict an intracellular flux distribution from one transcriptomic sample by
turning gene expression into reaction-level weights through the GPR, then either scaling
reaction bounds by expression (E-Flux2, Kim 2016) or fitting fluxes to expression-derived
targets (LAD). They complement the two-state ``differential`` module used by revert-metabolism.
"""

from __future__ import annotations

import ast
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

import numpy as np
from cobra import Model

from cmm.core import solvers
from cmm.core.flux_state import FluxState


@dataclass(frozen=True)
class OmicsFluxResult:
    """Flux distribution predicted from expression data."""

    method: str
    status: str
    objective_value: float | None
    fluxes: dict[str, float] = field(default_factory=dict)
    detail: str = ""

    def to_flux_state(self, name: str = "omics") -> FluxState:
        return FluxState(self.fluxes, name=name, provenance="imported")


def _eval_gpr_value(node: ast.AST | None, expr: Mapping[str, float]) -> float | None:
    """Continuous GPR evaluation: AND (complex) -> min, OR (isozymes) -> max.

    A missing gene contributes ``None``; an AND with any missing subunit is ``None`` (the
    complex cannot form); an OR keeps the available alternatives.
    """

    if node is None:
        return None
    if isinstance(node, ast.Module):
        body = getattr(node, "body", None)
        if isinstance(body, list):
            body = body[0] if body else None
        return _eval_gpr_value(body, expr)
    if isinstance(node, ast.Expression):
        return _eval_gpr_value(node.body, expr)
    if isinstance(node, ast.Name):
        return expr.get(node.id)
    if isinstance(node, ast.BoolOp):
        values = [_eval_gpr_value(v, expr) for v in node.values]
        if isinstance(node.op, ast.And):
            return None if any(v is None for v in values) else min(values)
        available = [v for v in values if v is not None]
        return max(available) if available else None
    return None


def gene_to_reaction_weights(
    model: Model, gene_expression: Mapping[str, float]
) -> dict[str, float]:
    """Map gene expression to reaction weights through each reaction's GPR."""

    lookup: dict[str, float] = {}
    for gene in model.genes:
        for key in (gene.id, gene.id.upper(), gene.id.lower()):
            if key in gene_expression:
                lookup[gene.id] = float(gene_expression[key])
                break

    weights: dict[str, float] = {}
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        value = _eval_gpr_value(rxn.gpr, lookup)
        if value is not None:
            weights[rxn.id] = float(value)
    return weights


def _is_excluded(rxn) -> bool:
    return rxn.boundary or "biomass" in rxn.id.lower() or "growth" in rxn.id.lower()


def eflux2(
    model: Model,
    reaction_weights: Mapping[str, float],
    *,
    objective_fraction: float = 0.99,
    min_scale: float = 1e-3,
    weight_threshold: float = 0.0,
    normalization_percentile: float = 100.0,
    exclude_exchange: bool = True,
    excluded_reactions: Iterable[str] | None = None,
) -> OmicsFluxResult:
    """E-Flux2: scale reaction bounds by normalized expression, then minimize the L2 flux.

    Scales each non-excluded reaction's bounds by ``clip(weight / denom, min_scale, 1)`` where
    ``denom`` is the ``normalization_percentile`` of the above-threshold weights; runs FBA for
    the optimum, holds the objective at ``objective_fraction`` of it, then minimizes the
    L2-norm of the flux vector (QP) for a unique distribution. Falls back to pFBA (L1) when no
    QP solver is available.
    """

    valid = {rid: abs(w) for rid, w in reaction_weights.items() if rid in model.reactions}
    above = {rid: w for rid, w in valid.items() if w > weight_threshold}
    if not above:
        return OmicsFluxResult("eflux2", "no_targets", None, {})

    denom = float(np.percentile(list(above.values()), normalization_percentile)) or 1.0
    excluded = set(excluded_reactions or ())

    with model:
        for rid, weight in above.items():
            rxn = model.reactions.get_by_id(rid)
            if rid in excluded or (exclude_exchange and _is_excluded(rxn)):
                continue
            ratio = min(1.0, max(min_scale, weight / denom))
            lb, ub = rxn.lower_bound, rxn.upper_bound
            rxn.bounds = (lb * ratio if lb < 0 else lb, ub * ratio if ub > 0 else ub)
        # Reactions below threshold are squeezed to the noise floor (E-Flux semantics).
        for rid in valid.keys() - above.keys():
            rxn = model.reactions.get_by_id(rid)
            if rid in excluded or (exclude_exchange and _is_excluded(rxn)):
                continue
            lb, ub = rxn.lower_bound, rxn.upper_bound
            rxn.bounds = (lb * min_scale if lb < 0 else lb, ub * min_scale if ub > 0 else ub)

        optimum = model.slim_optimize()
        if optimum is None or not np.isfinite(optimum):
            return OmicsFluxResult("eflux2", "infeasible", None, {})

        prob = model.problem
        objective_expr = model.objective.expression
        # Hold the biological objective near its optimum. For a maximized objective that is a
        # lower bound; for a minimized one it is an upper bound. The min-direction band must
        # relax UPWARD from the optimum regardless of sign — dividing a *negative* minimized
        # optimum by objective_fraction tightens past the optimum and makes the QP infeasible,
        # so relax by a fraction of |optimum| instead (equals optimum/fraction for optimum>0).
        if model.objective_direction == "min":
            slack = abs(optimum) * (1.0 / objective_fraction - 1.0)
            floor = prob.Constraint(objective_expr, ub=optimum + slack, name="_eflux2_floor")
        else:
            floor = prob.Constraint(objective_expr, lb=objective_fraction * optimum, name="_eflux2_floor")
        model.add_cons_vars([floor])

        if solvers.supports("QP", model.solver.interface):
            l2 = 0.0
            for rxn in model.reactions:
                l2 = l2 + rxn.flux_expression * rxn.flux_expression
            model.objective = prob.Objective(l2, direction="min")
            solution = model.optimize()
            detail = "qp"
        else:
            from cobra.flux_analysis import pfba

            solution = pfba(model)
            detail = "pfba"

        fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
        return OmicsFluxResult("eflux2", solution.status, float(optimum), fluxes, detail)


def lad(
    model: Model,
    reaction_weights: Mapping[str, float],
    *,
    scaling_factor: float = 1.0,
    weight_threshold: float = 0.01,
) -> OmicsFluxResult:
    """LAD: fit fluxes to expression-derived targets by minimizing total absolute deviation.

    For each target reaction, ``flux = target + d+ - d-`` with ``target = |weight| * scaling``;
    minimizes ``sum(d+ + d-)`` (an LP, runs on any solver).
    """

    targets = {
        rid: abs(w)
        for rid, w in reaction_weights.items()
        if rid in model.reactions and abs(w) > weight_threshold
    }
    if not targets:
        return OmicsFluxResult("lad", "no_targets", None, {})

    with model:
        prob = model.problem
        deviations = []
        constraints = []
        terms = []
        for rid, weight in targets.items():
            rxn = model.reactions.get_by_id(rid)
            target_flux = weight * scaling_factor
            d_pos = prob.Variable(f"lad_pos_{rid}", lb=0)
            d_neg = prob.Variable(f"lad_neg_{rid}", lb=0)
            deviations.extend([d_pos, d_neg])
            terms.extend([d_pos, d_neg])
            constraints.append(
                prob.Constraint(rxn.flux_expression - d_pos + d_neg,
                                lb=target_flux, ub=target_flux, name=f"lad_{rid}")
            )
        # Add variables first, then the constraints that reference them (optlang-safe order).
        model.add_cons_vars(deviations)
        model.solver.update()
        model.add_cons_vars(constraints)
        model.objective = prob.Objective(sum(terms), direction="min")

        solution = model.optimize()
        fluxes = {rid: float(v) for rid, v in solution.fluxes.items()}
        objective = (
            float(solution.objective_value) if solution.objective_value is not None else None
        )
        return OmicsFluxResult("lad", solution.status, objective, fluxes)


def integrate_expression(
    model: Model,
    gene_expression: Mapping[str, float],
    *,
    method: str = "eflux2",
    **kwargs,
) -> OmicsFluxResult:
    """Convenience: gene expression -> reaction weights -> predicted fluxes (E-Flux2 or LAD)."""

    weights = gene_to_reaction_weights(model, gene_expression)
    if method == "eflux2":
        return eflux2(model, weights, **kwargs)
    if method == "lad":
        return lad(model, weights, **kwargs)
    raise ValueError(f"unknown omics method {method!r}; use 'eflux2' or 'lad'")
