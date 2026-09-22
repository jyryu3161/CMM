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
from dataclasses import dataclass, replace
import time

import pandas as pd
from cobra import Model

from cmm.core.simulation import FluxSolution
from cmm.jev.actions import Intervention

#: One label, used by every branch of the FSEOF row so the comparison table cannot end up
#: with two spellings of the same method depending on whether it succeeded.
_HEADROOM_LABEL = "best amplification (outside the vocabulary)"

#: Likewise for the other rows whose label is referred to from more than one place.
_SINGLE_GENE_LABEL = "best single gene deletion"
_AGENT_LABEL = "JEV agent"

#: The control the agent has to beat to have contributed anything: the best deterministic
#: design plus one knockdown, chosen by exhausting every knockdown the agent could have made.
#: Same vocabulary, same growth floor, same gene resolution — no judgement anywhere in it.
_SWEEP_LABEL = "best deterministic design + one knockdown (exhaustive)"


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
    #: The worst product this design could make while growing as fast as it can, measured
    #: loopless. This is what the table is ranked on. ``product_flux`` beside it is the best
    #: case: a design whose guarantee is zero is one the strain may grow just as fast without
    #: ever using, however good its pFBA number looks.
    guaranteed_product: float | None = None
    #: The *worst* product this design could make with growth held at the common rate defined
    #: in :func:`compare_with_baselines` — the lowest maximum growth among the compared
    #: designs. Without it a row that trades growth for product looks like a better method
    #: rather than the same method read at a different point of the same trade-off. It is the
    #: worst and not the best for the same reason ``guaranteed_product`` is: read at its
    #: maximum, the unmodified wild type scores higher here than any design, because a network
    #: that is free to do anything is free to make the product too.
    guaranteed_at_matched_growth: float | None = None
    #: The method whose design this one contains, if any. A row seeded with a deterministic
    #: answer and then extended is not an independent result, and the table has to say so.
    contains_design: str | None = None
    #: The reactions this design constrains once its gene edits are resolved, so two rows can
    #: be checked for being the same kind of object rather than assumed to be.
    constrained: tuple[str, ...] = ()

    def to_row(self) -> dict[str, object]:
        return {
            "method": self.method,
            "n_interventions": len(self.design),
            "design": "; ".join(self.design),
            "guaranteed_product": self.guaranteed_product,
            "product_flux": self.product_flux,
            "growth": self.growth,
            "guaranteed_at_matched_growth": self.guaranteed_at_matched_growth,
            "contains_design": self.contains_design,
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


def _guarantee(
    model: Model,
    bounds: Mapping[str, tuple[float, float]],
    *,
    product: str,
    biomass: str,
) -> float | None:
    """The worst product this design could make while growing as fast as it can."""

    from cmm.jev.state import guaranteed_product

    with model:
        for reaction_id, (lower, upper) in bounds.items():
            model.reactions.get_by_id(reaction_id).bounds = (lower, upper)
        measured = guaranteed_product(model, product=product, biomass=biomass)
    return None if measured is None else float(measured[0])


def _guarantee_at_growth(
    model: Model,
    bounds: Mapping[str, tuple[float, float]],
    *,
    product: str,
    biomass: str,
    growth: float,
) -> float | None:
    """The least product this design must make with growth *held* at ``growth``.

    The column that makes the table a comparison. Two designs read at their own maximum growth
    rates are two points on two different trade-off curves, and the difference between them
    says as much about where each one sits on its curve as about which curve is better. Pinning
    growth reads every design at one operating point.

    The minimum and not the maximum, for the reason the guarantee is the minimum everywhere
    else here: held at a low growth rate the *unmodified* model can reach more product than any
    design, because nothing stops it — and a column the wild type wins is not measuring design.
    """

    with model:
        for reaction_id, (lower, upper) in bounds.items():
            model.reactions.get_by_id(reaction_id).bounds = (lower, upper)
        reaction = model.reactions.get_by_id(biomass)
        reaction.bounds = (growth, growth)
        model.objective = product
        model.objective_direction = "min"
        value = model.slim_optimize()
    return None if value != value else float(value)


def _gene_resolved_bounds(
    model: Model, knockouts: Sequence[str]
) -> dict[str, tuple[float, float]]:
    """A deletion set as the gene edit that achieves it, collateral included.

    The deterministic designers name reactions, and this used to score them by zeroing those
    reactions alone. The agent's own row has always carried the full consequence of the gene
    edit — isozymes deleted together, a shared gene taking its other reactions with it — so
    the two rows described different kinds of object in the same column. They are both the
    strain that would be built now.
    """

    from cmm.jev.actions import (
        ACTION_CATALOGUE,
        ActionNotApplicable,
        build_intervention,
    )

    reference = {r.id: 0.0 for r in model.reactions}
    bounds: dict[str, tuple[float, float]] = {}
    for reaction_id in knockouts:
        try:
            intervention = build_intervention(
                model, reaction_id, ACTION_CATALOGUE["knockout"], reference
            )
        except (ActionNotApplicable, KeyError):
            # A reaction with no gene association is an honest bound edit and nothing more.
            bounds[reaction_id] = (0.0, 0.0)
            continue
        for rid, lower, upper in intervention.bounds:
            bounds[rid] = (lower, upper)
    return bounds


def _gene_deletion_bounds(model: Model, gene_id: str) -> dict[str, tuple[float, float]]:
    """Every reaction one gene's loss stops, as bounds."""

    from cobra.manipulation import knock_out_model_genes

    if (
        gene_id not in model.genes
    ):  # pragma: no cover - the screen names the model's genes
        return {}
    with model:
        stopped = knock_out_model_genes(model, [gene_id])
        return {reaction.id: (0.0, 0.0) for reaction in stopped}


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
    #: Each row's design as bounds, so the comparison columns below can re-read every design
    #: through the same lens instead of each row measuring itself its own way.
    design_bounds: dict[str, dict[str, tuple[float, float]]] = {}

    def keep(
        row: BaselineRow, bounds: Mapping[str, tuple[float, float]] | None = None
    ) -> BaselineRow:
        rows.append(row)
        design_bounds[row.method] = dict(bounds or {})
        return row

    wild = _evaluate(model, {})
    keep(
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
        keep(*_single_gene_row(model, product, biomass, growth_floor))

    for label, solver in (("OptKnock", "optknock"), ("RobustKnock", "robustknock")):
        keep(
            *_strain_design_row(
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

    # The control. Built on whichever deterministic design scored best so far, because that is
    # the design the agent itself was seeded with and started from.
    proven = [
        row
        for row in rows
        if row.deterministic and row.design and row.status == "optimal"
    ]
    if proven:
        base = max(proven, key=lambda row: row.product_flux)
        keep(
            *_knockdown_sweep_row(
                model,
                product=product,
                biomass=biomass,
                growth_floor=growth_floor,
                base_label=base.method,
                base_design=base.design,
                base_bounds=design_bounds[base.method],
            )
        )

    # The agent's design, if there is one, is what the headroom row probes on top of; with no
    # design it probes the wild type, which is the honest thing to compare a wild type to.
    # Every reaction each gene edit constrains, not only the one the agent named: the row has
    # to score the strain that would actually be built.
    bounds = {
        rid: (low, high)
        for intervention in jev_interventions
        for rid, low, high in intervention.bounds
    }
    keep(*_amplification_headroom_row(model, product, biomass, growth_floor, bounds))

    if jev_interventions:
        solution = _evaluate(model, bounds)
        keep(
            BaselineRow(
                method=_AGENT_LABEL,
                design=tuple(i.describe() for i in jev_interventions),
                product_flux=float(solution.fluxes.get(product, 0.0)),
                growth=float(solution.fluxes.get(biomass, 0.0)),
                seconds=float("nan"),
                deterministic=False,
                status=solution.status,
                constrained=tuple(sorted(bounds)),
                note=(
                    "the agent's choices are not guaranteed to repeat; this is one run, not "
                    "the method's performance"
                ),
            ),
            bounds,
        )

    return _fill_comparison_columns(
        rows, design_bounds, model=model, product=product, biomass=biomass
    )


def _fill_comparison_columns(
    rows: Sequence[BaselineRow],
    design_bounds: Mapping[str, Mapping[str, tuple[float, float]]],
    *,
    model: Model,
    product: str,
    biomass: str,
) -> tuple[BaselineRow, ...]:
    """Measure every scorable design on the two columns that make the table a comparison.

    ``guaranteed_product`` is what each design must make; ``product_at_matched_growth`` is what
    each one can make at a single shared growth rate. Both are needed because a design bought by
    spending growth reads as a better method in the pFBA column alone, and a design bought by
    coupling reads as a worse one.

    The shared rate is the *lowest* maximum growth among the scorable designs, which is the only
    rate every design in the table can reach.
    """

    scorable = [
        row
        for row in rows
        if row.status == "optimal" and row.growth == row.growth and row.design
    ]
    matched_growth = min((row.growth for row in scorable), default=None)

    filled: list[BaselineRow] = []
    for row in rows:
        if row.status != "optimal" or row.product_flux != row.product_flux:
            filled.append(row)
            continue
        bounds = design_bounds.get(row.method, {})
        guaranteed = _guarantee(model, bounds, product=product, biomass=biomass)
        matched = (
            _guarantee_at_growth(
                model, bounds, product=product, biomass=biomass, growth=matched_growth
            )
            if matched_growth is not None and row.growth >= matched_growth - 1e-9
            else None
        )
        contains = row.contains_design or _containing_design(row, rows)
        filled.append(
            replace(
                row,
                guaranteed_product=guaranteed,
                guaranteed_at_matched_growth=matched,
                contains_design=contains,
            )
        )
    return tuple(filled)


def _containing_design(row: BaselineRow, rows: Sequence[BaselineRow]) -> str | None:
    """The deterministic design this row's own design is a superset of, if any.

    A run seeded with the deterministic answer and given a move that applies it whole can end
    on that answer plus one edit, and then be scored against it. That is not two methods
    disagreeing; it is one method and an increment. The table has to say which.
    """

    if row.method in (_HEADROOM_LABEL, "wild type") or not row.constrained:
        return None
    mine = set(row.constrained)
    contained = [
        other
        for other in rows
        if other is not row
        and other.deterministic
        and other.constrained
        and other.method not in (_HEADROOM_LABEL, _SWEEP_LABEL)
        and set(other.constrained) < mine
    ]
    # Only the largest ones. OptKnock's design contains the best single deletion here, so
    # naming both says the same thing twice and buries the one that matters.
    maximal = [
        other
        for other in contained
        if not any(
            set(other.constrained) < set(bigger.constrained) for bigger in contained
        )
    ]
    return ", ".join(sorted(row.method for row in maximal)) or None


def _single_gene_row(
    model: Model, product: str, biomass: str, growth_floor: float
) -> tuple[BaselineRow, dict[str, tuple[float, float]]]:
    """Best single-gene deletion, chosen by the MOMA screen and scored under pFBA.

    The screen is MOMA-L2 because the question it answers — what does the cell do immediately
    after one gene goes — is MOMA's question. The *score* in the table is not: the row used to
    report the minimal-adjustment product straight from the screen while every other row
    reported a re-optimised pFBA product, so one column held two quantities. The gene's
    deletion is re-applied and re-solved here like everything else, and the MOMA number it was
    chosen on is kept in the note.
    """

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
        return (
            BaselineRow(
                method=_SINGLE_GENE_LABEL,
                design=(),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=time.perf_counter() - started,
                deterministic=True,
                status="failed",
                note=f"the screen could not run: {error}",
            ),
            {},
        )

    viable = [
        row
        for row in screen
        if row.status == "optimal" and row.objective >= growth_floor
    ]
    viable.sort(key=lambda row: (-row.product_flux, row.target_id))
    elapsed = time.perf_counter() - started
    if not viable:
        return (
            BaselineRow(
                method=_SINGLE_GENE_LABEL,
                design=(),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=elapsed,
                deterministic=True,
                status="none viable",
                note=f"no single deletion of {len(screen)} held the growth floor",
            ),
            {},
        )
    best = viable[0]
    bounds = _gene_deletion_bounds(model, best.target_id)
    solution = _evaluate(model, bounds)
    return (
        BaselineRow(
            method=_SINGLE_GENE_LABEL,
            design=(best.target_id,),
            product_flux=float(solution.fluxes.get(product, 0.0)),
            growth=float(solution.fluxes.get(biomass, 0.0)),
            seconds=elapsed,
            deterministic=True,
            status=solution.status,
            constrained=tuple(sorted(bounds)),
            note=(
                f"best of {len(screen)} genes screened, chosen on the MOMA-L2 screen "
                f"({best.product_flux:.4g} product at the minimal-adjustment state) and "
                "re-scored here under pFBA like every other row"
            ),
        ),
        bounds,
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
) -> tuple[BaselineRow, dict[str, tuple[float, float]]]:
    """The designer's best design by guaranteed product, re-scored under pFBA.

    The designer names reactions; the row scores the **gene edit** that achieves them, so the
    collateral a laboratory would actually get is in the number. Scoring the named reactions
    alone put a different kind of object in the same column as the agent's row, which has
    carried its gene edits' full consequence from the start.
    """

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
        return (
            BaselineRow(
                method=label,
                design=(),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=time.perf_counter() - started,
                deterministic=True,
                status="failed",
                note=f"{label} could not run: {error}",
            ),
            {},
        )
    elapsed = time.perf_counter() - started
    if not result.designs:
        return (
            BaselineRow(
                method=label,
                design=(),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=elapsed,
                deterministic=True,
                status="no design",
                note=f"{label} returned no design under these bounds",
            ),
            {},
        )

    # CMM's own rule: rank designs by guaranteed, not maximum, product.
    best = max(result.designs, key=lambda design: design.guaranteed_product)
    bounds = _gene_resolved_bounds(model, best.knockouts)
    solution = _evaluate(model, bounds)
    collateral = sorted(set(bounds) - set(best.knockouts))
    return (
        BaselineRow(
            method=label,
            design=tuple(best.knockouts),
            product_flux=float(solution.fluxes.get(product, 0.0)),
            growth=float(solution.fluxes.get(biomass, 0.0)),
            seconds=elapsed,
            deterministic=True,
            status=solution.status,
            constrained=tuple(sorted(bounds)),
            note=(
                f"best of {len(result.designs)} designs by guaranteed product "
                f"({best.guaranteed_product:.4g}); complete deletions only — this formulation "
                "cannot express a partial knockdown"
                + (
                    "; scored as the gene edit, which also stops "
                    + ", ".join(collateral)
                    if collateral
                    else ""
                )
            ),
        ),
        bounds,
    )


def _knockdown_sweep_row(
    model: Model,
    *,
    product: str,
    biomass: str,
    growth_floor: float,
    base_label: str,
    base_design: Sequence[str],
    base_bounds: Mapping[str, tuple[float, float]],
) -> tuple[BaselineRow, dict[str, tuple[float, float]]]:
    """The best deterministic design plus the single best knockdown, found by exhausting them.

    This is the control that says what the decision model contributed. The agent's claim to add
    something the designers cannot is specific: OptKnock's variables are present-or-absent, so a
    knockdown on top of a proven design is a move it could not have considered. True — but
    "could not have considered" is not "requires judgement to find". There are only as many such
    moves as there are reactions carrying flux, each costs one solve, and trying all of them
    takes about a second on ``e_coli_core``.

    So the row is built the way the agent would have to build it: the same proven design, the
    same ``knockdown_50``, resolved through the same GPR, held to the same growth floor, ranked
    on the same guaranteed product. Whatever the agent's row has over this one is what the
    judgement bought; whatever it does not is what exhaustive enumeration was already going to
    find.
    """

    from cmm.jev.actions import (
        ACTION_CATALOGUE,
        ActionNotApplicable,
        build_intervention,
    )
    from cmm.jev.engine import _solve

    started = time.perf_counter()
    if not base_design:
        return (
            BaselineRow(
                method=_SWEEP_LABEL,
                design=(),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=time.perf_counter() - started,
                deterministic=True,
                status="no design",
                note="no deterministic design was found to build on",
            ),
            {},
        )

    with model:
        for reaction_id, (lower, upper) in base_bounds.items():
            model.reactions.get_by_id(reaction_id).bounds = (lower, upper)
        reference = _solve(model)
        if (
            reference.status != "optimal"
        ):  # pragma: no cover - the base row already solved
            return (
                BaselineRow(
                    method=_SWEEP_LABEL,
                    design=(),
                    product_flux=float("nan"),
                    growth=float("nan"),
                    seconds=time.perf_counter() - started,
                    deterministic=True,
                    status=reference.status,
                    note=f"the {base_label} design does not solve, so nothing can be added",
                ),
                {},
            )
        # "Half of wild type" means half of what the reaction carries in the state the
        # knockdown is applied to, which here is the proven design — exactly what the agent's
        # own screen measures against.
        reference_fluxes = dict(reference.fluxes)
        best: tuple[float, str, dict[str, tuple[float, float]], str] | None = None
        n_tried = 0
        for reaction in list(model.reactions):
            # The same universe the agent plays on: a reaction with no gene association is a
            # bound edit nobody can build, and offering the control a move the agent was never
            # offered would stop it being a control.
            if (
                reaction.id in base_bounds
                or reaction.id == product
                or not reaction.genes
                or reaction.objective_coefficient != 0
            ):
                continue
            try:
                trial = build_intervention(
                    model,
                    reaction.id,
                    ACTION_CATALOGUE["knockdown_50"],
                    reference_fluxes,
                )
            except (ActionNotApplicable, KeyError):
                continue
            n_tried += 1
            saved = [
                (rid, model.reactions.get_by_id(rid).bounds)
                for rid, _, _ in trial.bounds
            ]
            for rid, lower, upper in trial.bounds:
                model.reactions.get_by_id(rid).bounds = (lower, upper)
            solution = _solve(model)
            growth = float(solution.fluxes.get(biomass, 0.0))
            guaranteed = (
                _guarantee(model, {}, product=product, biomass=biomass)
                if solution.status == "optimal" and growth >= growth_floor
                else None
            )
            for rid, bounds in saved:
                model.reactions.get_by_id(rid).bounds = bounds
            if guaranteed is None:
                continue
            if best is None or guaranteed > best[0]:
                combined = dict(base_bounds)
                for rid, lower, upper in trial.bounds:
                    combined[rid] = (lower, upper)
                best = (guaranteed, reaction.id, combined, trial.describe())

    elapsed = time.perf_counter() - started
    if best is None:
        return (
            BaselineRow(
                method=_SWEEP_LABEL,
                design=tuple(base_design),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=elapsed,
                deterministic=True,
                status="none viable",
                note=(
                    f"none of {n_tried} knockdowns on the {base_label} design held the "
                    "growth floor"
                ),
            ),
            {},
        )

    _, reaction_id, combined, described = best
    solution = _evaluate(model, combined)
    return (
        BaselineRow(
            method=_SWEEP_LABEL,
            design=(*base_design, described),
            product_flux=float(solution.fluxes.get(product, 0.0)),
            growth=float(solution.fluxes.get(biomass, 0.0)),
            seconds=elapsed,
            deterministic=True,
            status=solution.status,
            contains_design=base_label,
            constrained=tuple(sorted(combined)),
            note=(
                f"the {base_label} design plus the best of {n_tried} knockdowns, ranked on "
                f"guaranteed product; {reaction_id} won. No judgement anywhere in this row — "
                "it is what exhausting the agent's own vocabulary on the agent's own starting "
                "point already gives"
            ),
        ),
        combined,
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
) -> tuple[BaselineRow, dict[str, tuple[float, float]]]:
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

    def failed(
        status: str, note: str, found: tuple[str, ...] = ()
    ) -> tuple[BaselineRow, dict[str, tuple[float, float]]]:
        return (
            BaselineRow(
                method=_HEADROOM_LABEL,
                design=found,
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=time.perf_counter() - started,
                deterministic=True,
                status=status,
                note=note,
            ),
            {},
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

        best: tuple[str, float, float, str, tuple[float, float]] | None = None
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
                best = (target, flux, growth, probe[1], probe[0])

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
    target, flux, growth, level, forced = best
    # The design *plus* the forced reaction. This row is scored on its own bounds like every
    # other: handing it the design's alone gave it the design's guarantee, so the row that
    # exists to price the restriction reported a number identical to the design it was pricing.
    combined = {**dict(design), target: forced}
    return (
        BaselineRow(
            method=_HEADROOM_LABEL,
            design=(f"{target}: forced to {level}",),
            product_flux=flux,
            growth=growth,
            seconds=time.perf_counter() - started,
            deterministic=True,
            status="optimal",
            constrained=tuple(sorted(combined)),
            note=(
                f"best of FSEOF's top {_HEADROOM_TARGETS} on this design, worth "
                f"{flux - before:+.4g} over it. This move is OUTSIDE the agent's vocabulary, "
                f"which is deletions and knockdowns only; the row prices that restriction "
                f"rather than competing with it{tail}"
            ),
        ),
        combined,
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


def _rank(row: BaselineRow) -> float:
    """What a design is worth: its guarantee, or its pFBA product when it has none."""

    if row.guaranteed_product is not None:
        return row.guaranteed_product
    return row.product_flux


def comparison_summary(
    rows: Sequence[BaselineRow], *, product: str
) -> dict[str, object]:
    """The one paragraph a reader needs, stated without flattering the agent.

    Three things changed here, all of them because the old verdict could report a win that was
    not one. It ranked on the pFBA product, which is the best a design could do rather than the
    worst it must; it compared designs read at their own maximum growth rates, so a design that
    simply spent growth read as a better method; and it never said when the agent's design
    contained the design it was being scored against.
    """

    scored = [
        row
        for row in rows
        if row.status == "optimal" and row.product_flux == row.product_flux
    ]
    if not scored:
        return {"product": product, "verdict": "no method produced a scorable design"}
    best = max(scored, key=_rank)
    # The FSEOF row applies an amplification, which the agent is not permitted to make. It
    # belongs in the table — it is what the vocabulary restriction costs — but scoring the
    # agent against it would be scoring it on a move it was forbidden to play.
    # The published methods. The exhaustive sweep is deterministic too, but it is the control
    # for what the agent's judgement added rather than a method anyone would cite, so it gets
    # its own sentence instead of competing for this one.
    deterministic = [
        row
        for row in scored
        if row.deterministic
        and row.design
        and row.method not in (_HEADROOM_LABEL, _SWEEP_LABEL)
    ]
    best_deterministic = max(deterministic, key=_rank) if deterministic else None
    agent = next((row for row in scored if row.method == _AGENT_LABEL), None)
    amplification = next((row for row in scored if row.method == _HEADROOM_LABEL), None)
    sweep = next((row for row in scored if row.method == _SWEEP_LABEL), None)

    quantity = (
        "guaranteed product"
        if agent is not None and agent.guaranteed_product is not None
        else "pFBA product"
    )

    verdict: str
    if agent is None:
        verdict = "no agent design was scored"
    elif best_deterministic is None:
        verdict = "no deterministic method produced a design to compare against"
    else:
        mine, theirs = _rank(agent), _rank(best_deterministic)
        # Both numbers with both growth rates, never one percentage. The designs sit at
        # different points of the same trade-off, and a single ratio hides which.
        verdict = (
            f"the agent's design reaches {mine:.4g} {quantity} at a growth rate of "
            f"{agent.growth:.4g} per hour; {best_deterministic.method} reaches {theirs:.4g} "
            f"at {best_deterministic.growth:.4g}."
        )
        matched = agent.guaranteed_at_matched_growth
        theirs_matched = best_deterministic.guaranteed_at_matched_growth
        if (
            matched is not None
            and theirs_matched is not None
            and abs(agent.growth - best_deterministic.growth) > 1e-6
        ):
            # The comparison that survives the growth difference: both designs pinned to the
            # one growth rate every design in the table can hold.
            verb = (
                "still ahead of"
                if matched > theirs_matched + 1e-6
                else (
                    "level with" if abs(matched - theirs_matched) <= 1e-6 else "behind"
                )
            )
            verdict += (
                f" Those are different operating points, so read them held at one: with growth "
                f"fixed at the rate every design here can hold, the agent's design makes "
                f"{matched:.4g} and {best_deterministic.method}'s makes {theirs_matched:.4g} — "
                f"{verb} it. A design that buys product by spending growth has moved along the "
                "trade-off, not beaten it."
            )
        verdict += " One run is not the method's performance."

    if agent is not None and agent.contains_design:
        verdict += (
            f" The agent's design contains {agent.contains_design}'s, which it was seeded with "
            "and could adopt in one move, so what is its own is the increment on top and not "
            "the whole design."
        )

    # The control, stated in the same breath as the verdict: exhausting the agent's own
    # vocabulary on the agent's own starting point is not a method it gets to beat quietly.
    if agent is not None and sweep is not None:
        gap = _rank(agent) - _rank(sweep)
        if gap > 1e-6:
            verdict += (
                f" Exhaustively trying every knockdown on the same proven design reaches "
                f"{_rank(sweep):.4g}, which the agent beat by {gap:+.4g}."
            )
        else:
            verdict += (
                f" Exhaustively trying every knockdown on the same proven design reaches "
                f"{_rank(sweep):.4g} in {sweep.seconds:.3g} s and deterministically, "
                f"{-gap:+.4g} against the agent — so on this problem the judgement bought "
                "nothing enumeration was not already going to find."
            )

    # State the price of the restriction in the same breath as the verdict, so a reader is
    # never left to infer that deletions and knockdowns are all there was.
    if agent is not None and amplification is not None:
        gap = _rank(amplification) - _rank(agent)
        if gap > 1e-6 and amplification.design:
            verdict += (
                f" Forcing flux through {amplification.design[0].split(':')[0]} on top of "
                f"this design reaches {_rank(amplification):.4g}, {gap:+.4g} more — "
                "but that is an amplification, which this run does not allow itself, because "
                "a forced lower bound is not what over-expressing an enzyme does to a cell."
            )

    return {
        "product": product,
        "ranked_on": quantity,
        "best_method": best.method,
        "best_product_flux": best.product_flux,
        "best_guaranteed_product": best.guaranteed_product,
        "best_growth": best.growth,
        "agent_guaranteed_at_matched_growth": (
            agent.guaranteed_at_matched_growth if agent is not None else None
        ),
        "best_deterministic_method": (
            best_deterministic.method if best_deterministic else None
        ),
        "best_deterministic_product_flux": (
            best_deterministic.product_flux if best_deterministic else None
        ),
        "best_deterministic_guaranteed_product": (
            best_deterministic.guaranteed_product if best_deterministic else None
        ),
        "agent_product_flux": agent.product_flux if agent else None,
        "agent_guaranteed_product": agent.guaranteed_product if agent else None,
        "agent_contains_design": agent.contains_design if agent else None,
        "exhaustive_sweep_guaranteed_product": (
            sweep.guaranteed_product if sweep else None
        ),
        "verdict": verdict,
    }


__all__ = [
    "BaselineRow",
    "comparison_frame",
    "comparison_summary",
    "compare_with_baselines",
]
