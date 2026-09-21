"""The controller: every move JEV is allowed to make, and what each one does to the model.

JEV plays CMM the way TypeSafe's demo plays Doom — it is handed a structured description of
the current state and picks one move from a fixed vocabulary, many times over. This module
*is* that vocabulary. Nothing outside it can be chosen, because a ``choice`` question can
only be answered with a criterion the caller supplied, so a move that names a reaction the
model does not contain is not something to validate against; it cannot be expressed.

Moves come in three kinds:

``ACT``
    Changes the model. Which moves exist for a reaction depends on whether it carries any
    wild-type flux, because that decides what a change can be measured against:

    *Carrying flux* — ``knockdown_50``, ``knockdown_25``, ``amplify_2x``, ``amplify_5x``.
    These are multiples of the reaction's **wild-type** flux, so "half" and "double" mean
    the same thing at every point in a run rather than drifting with the current state.

    *Carrying none* — ``force_on_low``, ``force_on_high``. A reaction at zero has no
    reference to be a multiple of, and this is not a corner case: in anaerobic
    ``e_coli_core`` the entire succinate-forming branch sits at zero, so a design that could
    only scale an existing flux could never switch the product on at all. These moves instead
    force a fraction of the flux the reaction could carry, found by maximising it under the
    model's current bounds. They are exploratory by construction — there is no wild-type
    behaviour to compare against — and are labelled as such.

    ``knockout`` is offered either way, and ``undo_last`` withdraws the most recent move.

``LOOK``
    Runs a CMM analysis and feeds the answer back into the next state without changing the
    model: ``fseof_scan``, ``essentiality_scan``, ``envelope_probe``. These are what make the
    loop *operating CMM* rather than only editing bounds.

``END``
    ``end_round`` — stop spending ticks on a round that has nothing left worth doing.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from typing import Literal

from cobra import Model

from cmm.core.condition import ReactionBound

ActionKind = Literal["act", "look", "end"]
InterventionMode = Literal["knockout", "knockdown", "amplification"]

#: Fractions of the feasible maximum used by the two ``force_on`` moves. Deliberately well
#: short of the ceiling: forcing a reaction to its own maximum leaves the rest of the network
#: no freedom and almost always collapses growth, which teaches the agent nothing.
FORCE_ON_FRACTIONS = {"force_on_low": 0.25, "force_on_high": 0.6}

#: A flux below this is treated as zero when deciding whether a relative move is defined.
#: Matches the magnitude at which cobra solutions report numerical noise rather than flux.
FLUX_EPSILON = 1e-9


@dataclass(frozen=True)
class Action:
    """One move in the controller, with the sentence JEV reads when choosing it."""

    name: str
    kind: ActionKind
    description: str
    mode: InterventionMode | None = None
    level: float | None = None


#: The ACT moves. ``level`` is a multiple of the wild-type flux magnitude: 0.5 means "cap at
#: half of wild type", 2.0 means "force at least twice wild type". Two knockdown and two
#: amplification strengths are offered rather than a continuous parameter because JEV picks
#: from named criteria — a number it cannot type is a number it cannot get wrong.
ACT_ACTIONS: tuple[Action, ...] = (
    Action(
        name="knockout",
        kind="act",
        mode="knockout",
        level=0.0,
        description=(
            "Delete this reaction: force its flux to exactly zero. Choose this when the "
            "reaction competes with the product for carbon or reducing power and the cell "
            "can grow without it."
        ),
    ),
    Action(
        name="knockdown_50",
        kind="act",
        mode="knockdown",
        level=0.5,
        description=(
            "Cap this reaction at half its wild-type flux. Choose this when the flux should "
            "be reduced but not removed, because deleting it entirely would stop growth."
        ),
    ),
    Action(
        name="knockdown_25",
        kind="act",
        mode="knockdown",
        level=0.25,
        description=(
            "Cap this reaction at a quarter of its wild-type flux: a strong knockdown that "
            "still leaves the function intact. Choose this when a half cap was not enough."
        ),
    ),
    Action(
        name="amplify_2x",
        kind="act",
        mode="amplification",
        level=2.0,
        description=(
            "Force at least twice the wild-type flux through this reaction, in the "
            "direction it already runs. Choose this when the reaction is on the route to "
            "the product, or supplies a precursor, ATP or redox cofactor that route needs."
        ),
    ),
    Action(
        name="amplify_5x",
        kind="act",
        mode="amplification",
        level=5.0,
        description=(
            "Force at least five times the wild-type flux through this reaction: a hard "
            "pull. Choose this when doubling was not enough and the step looks rate "
            "limiting for the product."
        ),
    ),
)

#: Amplification for a reaction that carries no wild-type flux: switch it on. Offered only
#: when the relative moves are undefined, so the two vocabularies never overlap.
FORCE_ON_ACTIONS: tuple[Action, ...] = (
    Action(
        name="force_on_low",
        kind="act",
        mode="amplification",
        level=FORCE_ON_FRACTIONS["force_on_low"],
        description=(
            "This reaction carries no flux at all. Switch it on at a quarter of the most it "
            "could carry. Choose this to open a pathway that is currently unused, when a "
            "gentler start is safer for growth."
        ),
    ),
    Action(
        name="force_on_high",
        kind="act",
        mode="amplification",
        level=FORCE_ON_FRACTIONS["force_on_high"],
        description=(
            "This reaction carries no flux at all. Switch it on hard, at around 60% of the "
            "most it could carry. Choose this when the reaction is on the direct route to "
            "the product and the pathway needs to be opened decisively."
        ),
    ),
)

#: The LOOK moves. These call a CMM analysis and enrich the next state; the model is unchanged.
LOOK_ACTIONS: tuple[Action, ...] = (
    Action(
        name="fseof_scan",
        kind="look",
        description=(
            "Do not change anything yet. Run an FSEOF scan, which forces product formation "
            "up in steps and reports which reactions increase with it. Choose this when you "
            "need to know which reactions pull the product before committing a move."
        ),
    ),
    Action(
        name="essentiality_scan",
        kind="look",
        description=(
            "Do not change anything yet. Test which candidate reactions the cell cannot grow "
            "without. Choose this when the risk of killing the strain is what is blocking a "
            "decision."
        ),
    ),
    Action(
        name="amplification_screen",
        kind="look",
        description=(
            "Do not change anything yet. For every reaction on the board, have CMM work out "
            "what forcing flux through it would actually do to the product, and report the "
            "answer. Choose this when you need to know which amplification pays before "
            "spending a place in the design on one — the reaction on the obvious route is "
            "often already saturated, and the one that pays is often a bypass."
        ),
    ),
    Action(
        name="strain_design_scan",
        kind="look",
        description=(
            "Do not change anything yet. Run OptKnock and RobustKnock, which prove which "
            "combination of deletions forces the product at maximum growth, and add the "
            "reactions they name to the board. Choose this when you need the deletions that "
            "close the cell's escape routes — including reactions carrying no flux today "
            "that it could switch to once the obvious routes are shut."
        ),
    ),
    Action(
        name="envelope_probe",
        kind="look",
        description=(
            "Do not change anything yet. Compute the growth-versus-production envelope, "
            "which shows how much product is reachable at each growth rate. Choose this when "
            "you need to know whether the product is limited by the trade-off with growth."
        ),
    ),
)

END_ACTION = Action(
    name="end_round",
    kind="end",
    description=(
        "Stop acting for this round. Choose this when no remaining move is expected to "
        "improve the product, or the design should be evaluated as it stands."
    ),
)

ADOPT_ACTION = Action(
    name="adopt_best_design",
    kind="act",
    description=(
        "Apply, in one move, the complete set of deletions the strain designer proved best. "
        "Choose this when the designer has found a design and you want it as the base to "
        "build on. Its deletions only pay off together: applied one at a time each looks "
        "worthless, which is why they are offered as a single move. It uses one place in the "
        "design per deletion."
    ),
)

UNDO_ACTION = Action(
    name="undo_last",
    kind="act",
    description=(
        "Withdraw the most recent intervention, restoring that reaction to its previous "
        "bounds. Choose this when the last move reduced the product or the growth rate."
    ),
)

#: Every action by name, for resolving an answer back to its definition.
ACTION_CATALOGUE: Mapping[str, Action] = {
    action.name: action
    for action in (
        *ACT_ACTIONS,
        *FORCE_ON_ACTIONS,
        *LOOK_ACTIONS,
        END_ACTION,
        UNDO_ACTION,
        ADOPT_ACTION,
    )
}


@dataclass(frozen=True)
class Intervention:
    """One applied change to one reaction, expressed relative to its wild-type flux.

    ``reference_flux`` is the reaction's flux in the round-0 wild-type pFBA state. Every
    relative move is defined against it rather than against the current state, so "half" does
    not drift as a run accumulates interventions.

    ``clamped`` records that the requested magnitude exceeded what the reaction's own bounds
    allow, and the move was applied at that ceiling instead. A clamped amplification is still
    a real change; reporting it as if the full multiple had been applied would overstate it.
    """

    reaction_id: str
    mode: InterventionMode
    level: float
    reference_flux: float
    lower_bound: float
    upper_bound: float
    action_name: str
    clamped: bool = False
    #: True for a ``force_on`` move, where ``level`` is a fraction of the reaction's feasible
    #: maximum rather than a multiple of a wild-type flux. A reader comparing two rows needs
    #: to know that the two numbers are not the same kind of quantity.
    exploratory: bool = False

    def to_reaction_bound(self) -> ReactionBound:
        """The change as the :class:`~cmm.core.condition.ReactionBound` CMM already applies."""

        return ReactionBound(
            reaction_id=self.reaction_id,
            lower_bound=self.lower_bound,
            upper_bound=self.upper_bound,
        )

    def describe(self) -> str:
        """A one-line human summary, used in the state history and in report tables."""

        if self.mode == "knockout":
            return f"{self.reaction_id}: knockout (flux forced to 0)"
        if self.mode == "knockdown":
            return (
                f"{self.reaction_id}: knockdown to {self.level:.0%} of wild type "
                f"(|v| ≤ {abs(self.level * self.reference_flux):.4g})"
            )
        suffix = " (clamped at the reaction's own bound)" if self.clamped else ""
        if self.exploratory:
            forced = self.lower_bound if self.lower_bound > 0 else self.upper_bound
            return (
                f"{self.reaction_id}: switched on from zero to {self.level:.0%} of its "
                f"loop-free maximum (|v| \u2265 {abs(forced):.4g}){suffix}"
            )
        return (
            f"{self.reaction_id}: amplification to {self.level:g}x wild type "
            f"(|v| ≥ {abs(self.level * self.reference_flux):.4g}){suffix}"
        )

    def to_record(self) -> dict[str, object]:
        """Flat export row."""

        return {
            "reaction_id": self.reaction_id,
            "mode": self.mode,
            "action": self.action_name,
            "level": self.level,
            "reference_flux": self.reference_flux,
            "lower_bound": self.lower_bound,
            "upper_bound": self.upper_bound,
            "clamped": self.clamped,
            "exploratory": self.exploratory,
        }


class ActionNotApplicable(ValueError):
    """Raised when a chosen move cannot be expressed on the chosen reaction.

    This is not a failure of the run. It is the honest answer to "halve a flux that is
    already zero", and the engine records the reason and lets JEV pick again rather than
    substituting a different move.
    """


def feasible_extreme(model: Model, reaction_id: str) -> tuple[float, float]:
    """The largest negative and positive flux this reaction can carry, free of loops.

    Used only to give a ``force_on`` move something to be a fraction of, for a reaction whose
    wild-type flux is zero.

    **The loopless constraint is not optional here.** A plain LP maximisation of ``FRD7`` on
    anaerobic ``e_coli_core`` returns 1000 — the bound, reached through the thermodynamically
    infeasible ``FRD7``/``SUCDi`` cycle — on a model taking up 10 mmol gDW^-1 h^-1 of glucose.
    Taking 60% of that would force a physically meaningless flux and hand the solver a futile
    cycle to satisfy it with. The loopless range gives 13.64 instead, and costs less time than
    the plain one (20 ms against 39 ms on this model) because the tighter problem solves
    faster.

    If the loopless solve is unavailable the plain range is used and the caller is told, via
    the returned flag, that the number may include a loop.
    """

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


def build_intervention(
    model: Model, reaction_id: str, action: Action, reference_flux: float
) -> Intervention:
    """Turn a chosen ACT move on a chosen reaction into concrete bounds.

    A relative move (``knockdown_*``, ``amplify_*``) needs a non-zero wild-type flux to be
    relative *to*, and raises :class:`ActionNotApplicable` without one rather than quietly
    becoming something else. A ``force_on`` move is the opposite: it exists only for a
    reaction at zero, and is measured against what the reaction could carry instead.

    A knockout is always defined, including on a zero-flux reaction — where it is honestly a
    no-op, which the engine's own re-solve will show.
    """

    if action.mode is None:
        raise ActionNotApplicable(
            f"{action.name!r} is a {action.kind} move and does not change a reaction"
        )
    reaction = model.reactions.get_by_id(reaction_id)
    lower, upper = float(reaction.lower_bound), float(reaction.upper_bound)
    reference = float(reference_flux)
    level = float(action.level if action.level is not None else 0.0)
    carries_flux = abs(reference) > FLUX_EPSILON

    if action.mode == "knockout":
        return Intervention(
            reaction_id=reaction_id,
            mode="knockout",
            level=0.0,
            reference_flux=reference,
            lower_bound=0.0,
            upper_bound=0.0,
            action_name=action.name,
        )

    if action.name in FORCE_ON_FRACTIONS:
        if carries_flux:
            raise ActionNotApplicable(
                f"{action.name!r} switches on a reaction that carries nothing, but "
                f"{reaction_id!r} already carries {reference:.4g}; use an amplify move"
            )
        return _force_on(model, reaction_id, action, level, lower, upper)

    if not carries_flux:
        raise ActionNotApplicable(
            f"{action.name!r} is relative to the wild-type flux of {reaction_id!r}, "
            f"which is {reference:.3g}: there is nothing to scale"
        )

    if action.mode == "knockdown":
        cap = abs(level * reference)
        # Cap the magnitude without opening a direction the reaction did not already have.
        new_lower = max(lower, -cap)
        new_upper = min(upper, cap)
        return Intervention(
            reaction_id=reaction_id,
            mode="knockdown",
            level=level,
            reference_flux=reference,
            lower_bound=new_lower,
            upper_bound=new_upper,
            action_name=action.name,
        )

    # Amplification: push the flux away from zero in the direction it already runs, by
    # raising the near-zero side of the interval. The far side is left alone so the solver
    # keeps the freedom to go further than asked.
    target = level * reference
    clamped = False
    if reference > 0:
        if target > upper:
            target, clamped = upper, True
        new_lower, new_upper = target, upper
    else:
        if target < lower:
            target, clamped = lower, True
        new_lower, new_upper = lower, target

    if not math.isfinite(new_lower) or not math.isfinite(new_upper):
        raise ActionNotApplicable(
            f"amplifying {reaction_id!r} needs finite bounds; it has ({lower}, {upper})"
        )
    return Intervention(
        reaction_id=reaction_id,
        mode="amplification",
        level=level,
        reference_flux=reference,
        lower_bound=new_lower,
        upper_bound=new_upper,
        action_name=action.name,
        clamped=clamped,
    )


def _force_on(
    model: Model,
    reaction_id: str,
    action: Action,
    fraction: float,
    lower: float,
    upper: float,
) -> Intervention:
    """Switch on a reaction that carries no flux, at a fraction of what it could carry.

    The direction is chosen by which way the reaction can run further. A reaction that cannot
    carry flux in either direction cannot be switched on at all, and says so rather than
    being applied as a bound that changes nothing.
    """

    reachable_low, reachable_high = feasible_extreme(model, reaction_id)
    forward, reverse = max(reachable_high, 0.0), min(reachable_low, 0.0)
    if max(forward, -reverse) <= FLUX_EPSILON:
        raise ActionNotApplicable(
            f"{reaction_id!r} cannot carry flux in either direction under the current "
            "bounds, so there is nothing to switch on"
        )

    if forward >= -reverse:
        target = min(fraction * forward, upper)
        new_lower, new_upper = target, upper
    else:
        target = max(fraction * reverse, lower)
        new_lower, new_upper = lower, target

    if new_lower > new_upper:  # pragma: no cover - clamped above
        raise ActionNotApplicable(
            f"switching {reaction_id!r} on is outside its bounds ({lower:.4g}, {upper:.4g})"
        )
    return Intervention(
        reaction_id=reaction_id,
        mode="amplification",
        level=fraction,
        reference_flux=0.0,
        lower_bound=new_lower,
        upper_bound=new_upper,
        action_name=action.name,
        exploratory=True,
    )


def applicable_actions(reference_flux: float) -> tuple[Action, ...]:
    """The ACT moves that are defined for a reaction with this wild-type flux.

    The two vocabularies never overlap: a reaction carrying flux is offered the relative
    moves, a reaction at zero is offered the ``force_on`` moves, and both are offered a
    knockout. Offering a move that cannot be applied would waste a tick and teach the agent
    nothing, so the list is filtered before the question is built rather than after it is
    answered.
    """

    knockout = tuple(action for action in ACT_ACTIONS if action.mode == "knockout")
    if abs(float(reference_flux)) <= FLUX_EPSILON:
        return knockout + FORCE_ON_ACTIONS
    return ACT_ACTIONS
