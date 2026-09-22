"""Measure a JEV run against the deterministic methods it sits beside.

An agent that proposes a design is worth nothing unless someone asks what the established
methods give on the same problem. This module asks, and the answer on the shipped succinate
example is not flattering in the way one might hope: OptKnock reaches 9.91 mmol gDW⁻¹ h⁻¹ in
under a second, deterministically and provably, and an unaided JEV run does not match it.

What the agent adds is a different thing, and the comparison is what makes it visible.
OptKnock's variables are present-or-absent: its formulation cannot express "carry half as
much flux through this reaction". A knockdown on top of a proven design is therefore a move
the designer could not have considered, and is where an agent has something to add.

One row here is deliberately a move the agent may **not** make. The vocabulary is deletions
and knockdowns only, on the grounds that a forced lower bound is not what over-expression
does to a cell, and that restriction costs real product: on anaerobic succinate, forcing the
glyoxylate shunt on reaches 10.76 against 9.95 for the best deletion-and-knockdown design.
The FSEOF row prices that decision in every run rather than leaving it as an argument.

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

from cmm.core.simulation import FluxSolution
from cmm.jev.actions import Intervention

#: One label, used by every branch of the FSEOF row so the comparison table cannot end up
#: with two spellings of the same method depending on whether it succeeded.
_HEADROOM_LABEL = "best amplification on top of this design (outside the vocabulary)"


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

    # The agent's design, if there is one, is what the headroom row probes on top of; with no
    # design it probes the wild type, which is the honest thing to compare a wild type to.
    bounds = {i.reaction_id: (i.lower_bound, i.upper_bound) for i in jev_interventions}
    rows.append(
        _amplification_headroom_row(model, product, biomass, growth_floor, bounds)
    )

    if jev_interventions:
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
            f"({best.guaranteed_product:.4g}); complete deletions only — this formulation "
            "cannot express a partial knockdown"
        ),
    )


#: Fraction of a reaction's loop-free feasible maximum the amplification probe forces through
#: it when it carries no flux. Well short of the ceiling on purpose: forcing a reaction to its
#: own maximum leaves the rest of the network no freedom and almost always collapses growth,
#: which is not a fair reading of what over-expressing that enzyme would do.
_FORCE_ON_FRACTION = 0.25

#: How many of FSEOF's ranked targets the probe tries. The answer is almost always in the
#: first few, and each one costs a solve.
_HEADROOM_TARGETS = 10


def _amplification_headroom_row(
    model: Model,
    product: str,
    biomass: str,
    growth_floor: float,
    design: Mapping[str, tuple[float, float]],
) -> BaselineRow:
    """The best amplification available **on top of the design being scored**.

    This row is the price tag on a policy decision, measured rather than argued. The agent's
    vocabulary is deletions and knockdowns only, because a lower bound on a flux is not what
    over-expressing an enzyme does to a cell: it tells the solver the flux *must* be carried,
    by whatever route is cheapest, while stronger expression only raises a capacity the cell
    may decline to use. That restriction has a cost, and the honest place to show it is beside
    the agent's own result.

    **On top of the design, not on the wild type**, because that is where the question lives
    and the two answers are nothing alike. FSEOF's top amplification target on wild-type
    anaerobic ``e_coli_core`` buys exactly nothing. Ranked instead on a design that already
    deletes ``ACALD``, ``D_LACt2`` and ``THD2``, the same method puts the glyoxylate shunt
    fourth, and forcing it reaches 10.76 against the design's 9.95.

    It is also why this is computed per run instead of quoted as a constant. Whether that
    10.76 is *available* depends on how much growth the design has already spent: at a floor
    of 0.05 it is not, because it leaves growth at 0.041, and the best amplification that
    keeps the strain alive reaches 10.04. A headline percentage would have been true of one
    design and wrong about the next.

    The loop-free range is not optional for a reaction at zero. A plain LP maximisation of
    ``FRD7`` on anaerobic ``e_coli_core`` returns its 1000 bound through the thermodynamically
    infeasible ``FRD7``/``SUCDi`` cycle, on a model taking up 10 mmol gDW^-1 h^-1 of glucose.
    A quarter of that would be physically meaningless, and would flatter this row.
    """

    from cmm.features.production import fseof

    started = time.perf_counter()

    def failed(status: str, note: str, found: tuple[str, ...] = ()) -> BaselineRow:
        return BaselineRow(
            method=_HEADROOM_LABEL,
            design=found,
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=time.perf_counter() - started,
            deterministic=True,
            status=status,
            note=note,
        )

    with model:
        for reaction_id, applied in design.items():
            model.reactions.get_by_id(reaction_id).bounds = applied
        base = _evaluate(model, {})
        if base.status != "optimal":
            return failed("infeasible", "the design itself does not solve")
        before = float(base.fluxes.get(product, 0.0))
        try:
            ranked = fseof(model, product, biomass).amplification_targets()
        except Exception as error:
            return failed("failed", f"FSEOF could not run on this design: {error}")
        if not ranked:
            return failed("no target", "FSEOF ranked no amplification target")

        best: tuple[str, float, float, str] | None = None
        refused = 0
        for target in ranked[:_HEADROOM_TARGETS]:
            probe = _forced_bounds(model, target, float(base.fluxes.get(target, 0.0)))
            if probe is None:
                continue
            solution = _evaluate(model, {target: probe[0]})
            growth = float(solution.fluxes.get(biomass, 0.0))
            if solution.status != "optimal" or growth < growth_floor:
                refused += 1
                continue
            flux = float(solution.fluxes.get(product, 0.0))
            if best is None or flux > best[1]:
                best = (target, flux, growth, probe[1])

    tail = (
        f"; {refused} of the ranked targets were refused for dropping growth below the floor"
        if refused
        else ""
    )
    if best is None:
        return failed(
            "no viable target",
            f"no amplification among FSEOF's top {_HEADROOM_TARGETS} keeps the strain "
            f"above the growth floor on this design{tail}",
        )
    target, flux, growth, level = best
    return BaselineRow(
        method=_HEADROOM_LABEL,
        design=(f"{target}: forced to {level}",),
        product_flux=flux,
        growth=growth,
        seconds=time.perf_counter() - started,
        deterministic=True,
        status="optimal",
        note=(
            f"best of FSEOF's top {_HEADROOM_TARGETS} on this design, worth "
            f"{flux - before:+.4g} over it. This move is OUTSIDE the agent's vocabulary, "
            f"which is deletions and knockdowns only; the row prices that restriction rather "
            f"than competing with it{tail}"
        ),
    )


def _forced_bounds(
    model: Model, reaction_id: str, reference: float
) -> tuple[tuple[float, float], str] | None:
    """Bounds that force flux through a reaction, and a phrase describing the level.

    Twice its current flux when it carries one; a quarter of its loop-free maximum when it
    does not. ``None`` when the reaction cannot carry flux in either direction, which is not
    an error — there is simply nothing to switch on.
    """

    reaction = model.reactions.get_by_id(reaction_id)
    lower, upper = float(reaction.lower_bound), float(reaction.upper_bound)
    if abs(reference) > 1e-9:
        target = 2.0 * reference
        if reference > 0:
            return (
                min(target, upper),
                upper,
            ), f"at least {abs(target):.4g} (2x its flux)"
        return (lower, max(target, lower)), f"at least {abs(target):.4g} (2x its flux)"

    low, high = _loopless_extreme(model, reaction_id)
    forward, reverse = max(high, 0.0), min(low, 0.0)
    if max(forward, -reverse) <= 1e-9:
        return None
    if forward >= -reverse:
        target = min(_FORCE_ON_FRACTION * forward, upper)
        bounds = (target, upper)
    else:
        target = max(_FORCE_ON_FRACTION * reverse, lower)
        bounds = (lower, target)
    return bounds, (
        f"at least {abs(target):.4g} ({_FORCE_ON_FRACTION:.0%} of its loop-free maximum, "
        "it carrying no flux here)"
    )


def _loopless_extreme(model: Model, reaction_id: str) -> tuple[float, float]:
    """The largest negative and positive flux a reaction can carry, free of internal loops."""

    from cmm.core.simulation import fva

    try:
        ranges = fva(
            model,
            reactions=[reaction_id],
            fraction_of_optimum=0.0,
            loopless="fastSNP",
            processes=1,
        )
    except (
        Exception
    ):  # pragma: no cover - solver-specific; the plain range still bounds it
        ranges = fva(
            model, reactions=[reaction_id], fraction_of_optimum=0.0, processes=1
        )
    flux_range = ranges[reaction_id]
    return float(flux_range.minimum), float(flux_range.maximum)


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
    # The FSEOF row applies an amplification, which the agent is not permitted to make. It
    # belongs in the table — it is what the vocabulary restriction costs — but scoring the
    # agent against it would be scoring it on a move it was forbidden to play.
    deterministic = [
        row
        for row in scored
        if row.deterministic and row.design and row.method != _HEADROOM_LABEL
    ]
    best_deterministic = (
        max(deterministic, key=lambda row: row.product_flux) if deterministic else None
    )
    agent = next((row for row in scored if row.method == "JEV agent"), None)
    amplification = next((row for row in scored if row.method == _HEADROOM_LABEL), None)

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

    # State the price of the restriction in the same breath as the verdict, so a reader is
    # never left to infer that deletions and knockdowns are all there was.
    if agent is not None and amplification is not None:
        gap = amplification.product_flux - agent.product_flux
        if gap > 1e-6 and amplification.design:
            verdict += (
                f" Forcing flux through {amplification.design[0].split(':')[0]} on top of "
                f"this design reaches {amplification.product_flux:.4g}, {gap:+.4g} more — "
                "but that is an amplification, which this run does not allow itself, because "
                "a forced lower bound is not what over-expressing an enzyme does to a cell."
            )

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
