"""Measure a JEV run against the deterministic methods it sits beside.

An agent that proposes a design is worth nothing unless someone asks what the established
methods give on the same problem. This module asks, and the answer on the shipped succinate
example is not flattering in the way one might hope: OptKnock reaches 9.91 mmol gDW⁻¹ h⁻¹ in
under a second, deterministically and provably, and an unaided JEV run does not match it.

What the agent adds is a different thing, and the comparison is what makes it visible.
OptKnock searches knockouts only — its formulation cannot express "force more flux through
this reaction". Handed OptKnock's own proven design and one place left in the budget, JEV
added an amplification on the glyoxylate shunt and reached 10.76, above the deterministic
optimum, at a real cost in growth. Neither number means anything without the other.

Every design here is evaluated the same way: apply it, solve pFBA, read the product and the
growth. A comparison in which each method reports its own favourite quantity is not a
comparison.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import time

import pandas as pd
from cobra import Model

from cmm.core.simulation import FluxSolution, pfba
from cmm.jev.actions import Intervention


@dataclass(frozen=True)
class BaselineRow:
    """One method's best design, scored the same way as every other method's."""

    method: str
    design: tuple[str, ...]
    product_flux: float
    growth: float
    seconds: float
    deterministic: bool
    status: str = "optimal"
    note: str = ""

    def to_row(self) -> dict[str, object]:
        return {
            "method": self.method,
            "n_interventions": len(self.design),
            "design": "; ".join(self.design),
            "product_flux": self.product_flux,
            "growth": self.growth,
            "seconds": round(self.seconds, 3),
            "deterministic": self.deterministic,
            "status": self.status,
            "note": self.note,
        }


def _evaluate(model: Model, bounds: Mapping[str, tuple[float, float]]) -> FluxSolution:
    """Score one design: apply its bounds, solve, restore. Same lens for every method."""

    with model:
        for reaction_id, (lower, upper) in bounds.items():
            model.reactions.get_by_id(reaction_id).bounds = (lower, upper)
        from cmm.jev.engine import _solve

        return _solve(model)


def compare_with_baselines(
    model: Model,
    *,
    product: str,
    biomass: str,
    growth_floor: float,
    jev_interventions: Sequence[Intervention] = (),
    max_knockouts: int = 3,
    max_solutions: int = 5,
    seed: int = 0,
    run_single_gene_screen: bool = True,
) -> tuple[BaselineRow, ...]:
    """Score the JEV design and the deterministic methods on one problem, identically.

    ``model`` must already carry the condition the JEV run used; every method is given the
    same one. The row order is the order they are reported in, wild type first so every other
    number has something to be a change from.
    """

    rows: list[BaselineRow] = []

    wild = _evaluate(model, {})
    rows.append(
        BaselineRow(
            method="wild type",
            design=(),
            product_flux=float(wild.fluxes.get(product, 0.0)),
            growth=float(wild.fluxes.get(biomass, 0.0)),
            seconds=0.0,
            deterministic=True,
            status=wild.status,
            note="no intervention",
        )
    )

    if run_single_gene_screen:
        rows.append(_single_gene_row(model, product, biomass, growth_floor, seed=seed))

    for label, solver in (("OptKnock", "optknock"), ("RobustKnock", "robustknock")):
        rows.append(
            _strain_design_row(
                model,
                label,
                solver,
                product,
                biomass,
                growth_floor,
                max_knockouts=max_knockouts,
                max_solutions=max_solutions,
                seed=seed,
            )
        )

    rows.append(_fseof_row(model, product, biomass, growth_floor))

    if jev_interventions:
        bounds = {
            i.reaction_id: (i.lower_bound, i.upper_bound) for i in jev_interventions
        }
        solution = _evaluate(model, bounds)
        rows.append(
            BaselineRow(
                method="JEV agent",
                design=tuple(i.describe() for i in jev_interventions),
                product_flux=float(solution.fluxes.get(product, 0.0)),
                growth=float(solution.fluxes.get(biomass, 0.0)),
                seconds=float("nan"),
                deterministic=False,
                status=solution.status,
                note=(
                    "the agent's choices are not guaranteed to repeat; this is one run, not "
                    "the method's performance"
                ),
            )
        )

    return tuple(rows)


def _single_gene_row(
    model: Model, product: str, biomass: str, growth_floor: float, *, seed: int
) -> BaselineRow:
    """Best single-gene deletion by product flux, among those that stay viable."""

    from cmm.features import batch_comparison
    from cmm.features._perturbation import gene_perturbations
    from cmm.features.comparison import reference_flux

    started = time.perf_counter()
    try:
        reference = reference_flux(model, "pfba")
        screen = batch_comparison(
            model,
            reference,
            gene_perturbations(model),
            method="moma_l2",
            objective_reaction=biomass,
            product_reaction=product,
        )
    except Exception as error:
        return BaselineRow(
            method="best single gene deletion (MOMA-L2)",
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=time.perf_counter() - started,
            deterministic=True,
            status="failed",
            note=f"the screen could not run: {error}",
        )

    viable = [
        row
        for row in screen
        if row.status == "optimal" and row.objective >= growth_floor
    ]
    viable.sort(key=lambda row: (-row.product_flux, row.target_id))
    elapsed = time.perf_counter() - started
    if not viable:
        return BaselineRow(
            method="best single gene deletion (MOMA-L2)",
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=elapsed,
            deterministic=True,
            status="none viable",
            note=f"no single deletion of {len(screen)} held the growth floor",
        )
    best = viable[0]
    return BaselineRow(
        method="best single gene deletion (MOMA-L2)",
        design=(best.target_id,),
        product_flux=float(best.product_flux),
        growth=float(best.objective),
        seconds=elapsed,
        deterministic=True,
        note=(
            f"best of {len(screen)} genes screened; scored at the minimal-adjustment state, "
            "not the re-optimised one"
        ),
    )


def _strain_design_row(
    model: Model,
    label: str,
    solver: str,
    product: str,
    biomass: str,
    growth_floor: float,
    *,
    max_knockouts: int,
    max_solutions: int,
    seed: int,
) -> BaselineRow:
    """The designer's best design by guaranteed product, re-scored under pFBA."""

    from cmm.features import strain_design as sd

    started = time.perf_counter()
    try:
        result = getattr(sd, solver)(
            model,
            product,
            biomass=biomass,
            max_knockouts=max_knockouts,
            max_solutions=max_solutions,
            min_growth=growth_floor,
            seed=seed,
        )
    except Exception as error:
        return BaselineRow(
            method=label,
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=time.perf_counter() - started,
            deterministic=True,
            status="failed",
            note=f"{label} could not run: {error}",
        )
    elapsed = time.perf_counter() - started
    if not result.designs:
        return BaselineRow(
            method=label,
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=elapsed,
            deterministic=True,
            status="no design",
            note=f"{label} returned no design under these bounds",
        )

    # CMM's own rule: rank designs by guaranteed, not maximum, product.
    best = max(result.designs, key=lambda design: design.guaranteed_product)
    solution = _evaluate(model, {rid: (0.0, 0.0) for rid in best.knockouts})
    return BaselineRow(
        method=label,
        design=tuple(best.knockouts),
        product_flux=float(solution.fluxes.get(product, 0.0)),
        growth=float(solution.fluxes.get(biomass, 0.0)),
        seconds=elapsed,
        deterministic=True,
        status=solution.status,
        note=(
            f"best of {len(result.designs)} designs by guaranteed product "
            f"({best.guaranteed_product:.4g}); knockouts only — this formulation cannot "
            "express an amplification"
        ),
    )


def _fseof_row(
    model: Model, product: str, biomass: str, growth_floor: float
) -> BaselineRow:
    """FSEOF's top amplification target, forced on by the same rule the agent uses.

    FSEOF ranks reactions; it does not state how hard to push one. Scoring its top target at
    the level ``force_on_low`` would use makes the two comparable, and the choice is stated
    rather than left implicit.
    """

    from cmm.features.production import fseof
    from cmm.jev.actions import (
        ACTION_CATALOGUE,
        ActionNotApplicable,
        build_intervention,
    )

    started = time.perf_counter()
    try:
        targets = fseof(model, product, biomass).amplification_targets()
    except Exception as error:
        return BaselineRow(
            method="FSEOF top amplification target",
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=time.perf_counter() - started,
            deterministic=True,
            status="failed",
            note=f"FSEOF could not run: {error}",
        )
    if not targets:
        return BaselineRow(
            method="FSEOF top amplification target",
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=time.perf_counter() - started,
            deterministic=True,
            status="no target",
            note="FSEOF selected no amplification target",
        )

    reference = pfba(model).fluxes
    target = targets[0]
    action = ACTION_CATALOGUE[
        "force_on_low"
        if abs(float(reference.get(target, 0.0))) <= 1e-9
        else "amplify_2x"
    ]
    try:
        intervention = build_intervention(
            model, target, action, reference_flux=float(reference.get(target, 0.0))
        )
    except ActionNotApplicable as error:
        return BaselineRow(
            method="FSEOF top amplification target",
            design=(target,),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=time.perf_counter() - started,
            deterministic=True,
            status="not applicable",
            note=str(error),
        )
    solution = _evaluate(
        model, {target: (intervention.lower_bound, intervention.upper_bound)}
    )
    return BaselineRow(
        method="FSEOF top amplification target",
        design=(intervention.describe(),),
        product_flux=float(solution.fluxes.get(product, 0.0)),
        growth=float(solution.fluxes.get(biomass, 0.0)),
        seconds=time.perf_counter() - started,
        deterministic=True,
        status=solution.status,
        note=(
            f"top of {len(targets)} ranked targets, applied with the agent's own "
            f"{action.name} rule so the two are comparable"
        ),
    )


def comparison_frame(rows: Sequence[BaselineRow]) -> pd.DataFrame:
    """The comparison as a table, in reporting order."""

    return pd.DataFrame([row.to_row() for row in rows])


def comparison_summary(
    rows: Sequence[BaselineRow], *, product: str
) -> dict[str, object]:
    """The one paragraph a reader needs, stated without flattering the agent."""

    scored = [
        row
        for row in rows
        if row.status == "optimal" and row.product_flux == row.product_flux
    ]
    if not scored:
        return {"product": product, "verdict": "no method produced a scorable design"}
    best = max(scored, key=lambda row: row.product_flux)
    deterministic = [row for row in scored if row.deterministic and row.design]
    best_deterministic = (
        max(deterministic, key=lambda row: row.product_flux) if deterministic else None
    )
    agent = next((row for row in scored if row.method == "JEV agent"), None)

    verdict: str
    if agent is None:
        verdict = "no agent design was scored"
    elif best_deterministic is None:
        verdict = "no deterministic method produced a design to compare against"
    elif agent.product_flux > best_deterministic.product_flux + 1e-6:
        margin = agent.product_flux / best_deterministic.product_flux - 1.0
        verdict = (
            f"the agent's design beats the best deterministic one by {margin:.1%}, at a "
            f"growth rate of {agent.growth:.4g} against {best_deterministic.growth:.4g}. "
            "One run is not the method's performance."
        )
    elif agent.product_flux < best_deterministic.product_flux - 1e-6:
        margin = 1.0 - agent.product_flux / best_deterministic.product_flux
        verdict = (
            f"the agent's design is {margin:.1%} below {best_deterministic.method}, which "
            "found its answer deterministically and in less time."
        )
    else:
        verdict = f"the agent matched {best_deterministic.method}."

    return {
        "product": product,
        "best_method": best.method,
        "best_product_flux": best.product_flux,
        "best_growth": best.growth,
        "best_deterministic_method": (
            best_deterministic.method if best_deterministic else None
        ),
        "best_deterministic_product_flux": (
            best_deterministic.product_flux if best_deterministic else None
        ),
        "agent_product_flux": agent.product_flux if agent else None,
        "verdict": verdict,
    }


__all__ = [
    "BaselineRow",
    "comparison_frame",
    "comparison_summary",
    "compare_with_baselines",
]
