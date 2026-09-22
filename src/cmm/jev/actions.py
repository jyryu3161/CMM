"""The controller: every move JEV is allowed to make, and what each one does to the model.

JEV plays CMM the way TypeSafe's demo plays Doom — it is handed a structured description of
the current state and picks one move from a fixed vocabulary, many times over. This module
*is* that vocabulary. Nothing outside it can be chosen, because a ``choice`` question can
only be answered with a criterion the caller supplied, so a move that names a reaction the
model does not contain is not something to validate against; it cannot be expressed.

**The vocabulary is down-regulation only: deletion and partial knockdown.** It used to
include amplification, both as a multiple of an existing flux and as a ``force_on`` move that
switched a zero-flux reaction on at a fraction of its feasible maximum, and that is what
produced this project's best succinate design — 10.7613 against OptKnock's 9.9108. It was
removed deliberately, and the reason is not that it scored badly:

A lower bound on a reaction is not what over-expression does. Forcing ``v >= x`` tells the
solver the cell *must* carry that flux, and the solver will satisfy it through whatever route
is cheapest, including one the enzyme has nothing to do with. Stronger expression of an
enzyme raises a *capacity*; the cell still decides whether to use it. So the in-silico gain
from a forced lower bound is an upper bound on an upper bound, and it is the kind of number
that survives review and fails in a flask. A deletion and a 50% knockdown, by contrast, are
caps — they say what the cell *cannot* do, which is exactly what deleting a gene or
weakening its promoter achieves.

Measured cost of the restriction on anaerobic succinate in ``e_coli_core``: the best design
reachable with deletions and knockdowns is 9.9461 against 10.7613 with amplification, about
7.6%. The run records that, so the trade is visible rather than implied.

Moves come in three kinds:

``ACT``
    Changes the model.

    ``knockout`` forces the flux to exactly zero and is always defined, including on a
    reaction carrying nothing today — which is not a pointless move: OptKnock's most valuable
    deletions are escape routes the cell would switch to once the obvious ones are shut.

    ``knockdown_50`` caps the magnitude at half the reaction's **wild-type** flux, so "half"
    means the same thing at every point in a run rather than drifting with the current state.
    It needs a non-zero wild-type flux to be half *of*, and says so rather than quietly
    becoming something else. There is one knockdown strength, not two: a second, deeper cap
    mostly bought a second rejection of the same idea, and roughly halving an activity is the
    level a laboratory can actually aim at with a promoter swap or an RBS change.

    ``undo_last`` withdraws the most recent move.

``LOOK``
    Runs a CMM analysis and feeds the answer back into the next state without changing the
    model: ``fseof_scan``, ``essentiality_scan``, ``envelope_probe``, ``state_distance_check``,
    ``strain_design_scan``. These are what make the loop *operating CMM* rather than only
    editing bounds.

    One measurement is deliberately **not** a move. What deleting or halving each candidate
    would do to the product is recomputed by CMM whenever the design changes, because it is a
    fact and not a decision — and because leaving it as a move meant the agent skipped it.
    Given the cofactor reading it would infer a plausible answer and act on the inference
    instead of the measurement. Partial information displacing measurement is worse than no
    information.

``END``
    ``end_round`` — stop spending steps on a round that has nothing left worth doing.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

from cobra import Model

from cmm.core.condition import ReactionBound
from cmm.jev.genes import gene_names, resolve_gene_edit

ActionKind = Literal["act", "look", "end"]
InterventionMode = Literal["knockout", "knockdown"]

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


#: The ACT moves. ``level`` is a fraction of the wild-type flux magnitude: 0.0 means "hold at
#: zero", 0.5 means "cap at half of wild type". A named strength rather than a continuous
#: parameter because JEV picks from named criteria — a number it cannot type is a number it
#: cannot get wrong.
ACT_ACTIONS: tuple[Action, ...] = (
    Action(
        name="knockout",
        kind="act",
        mode="knockout",
        level=0.0,
        description=(
            "Delete the gene or genes behind this reaction, forcing its flux to zero. The "
            "record says which genes that is; if they are isozymes, all of them go, and if "
            "one of them also serves another reaction, that reaction stops too and CMM "
            "applies that as part of the move. Choose this when the reaction competes with "
            "the product for carbon or reducing power and the cell can grow without it, or "
            "when it is an alternative route the cell would escape down once the obvious "
            "ones are shut."
        ),
    ),
    Action(
        name="knockdown_50",
        kind="act",
        mode="knockdown",
        level=0.5,
        description=(
            "Weaken the gene or genes behind this reaction so it carries at most half its "
            "wild-type flux \u2014 a promoter swap or an RBS change rather than a deletion. "
            "Any other reaction those genes alone carry is weakened with it. Choose this "
            "when the flux should be reduced but not removed, because deleting it entirely "
            "would stop growth or cost more product than it frees."
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
            "up in steps and reports which reactions rise and which fall with it. A reaction "
            "whose flux FALLS as the product is forced up is one the product pathway does "
            "not need, and is therefore a candidate for deletion or knockdown. Choose this "
            "when you need to know which fluxes compete with the product."
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
        name="state_distance_check",
        kind="look",
        description=(
            "Do not change anything yet. Run MOMA and ROOM on the current design against the "
            "wild type, and report how far the cell has to move and how many reactions have "
            "to change. Choose this when you need to know whether the design is a small "
            "rewiring or a wholesale one — a design needing forty reactions to change is a "
            "harder strain to build than one needing five, at the same product flux."
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
        "worthless, which is why they are offered as a single move. It uses one deletion "
        "from the knockout budget per reaction."
    ),
)

RESTORE_ACTION = Action(
    name="restore_best_design",
    kind="act",
    description=(
        "Throw away the current design and go back to the best one this run has found, "
        "whichever round it came from. Choose this when the last rounds have made things "
        "worse and there is nothing left to learn from where the design is now."
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

#: The weaker version of a move, for when the stronger one was rejected for being too much.
#: A move turned down because it dropped growth below the floor has an obvious next thing to
#: try — the same move, gentler — and the agent will not find it on its own: watching a run,
#: a refused move sent it to a different reaction entirely, leaving on the table what the
#: gentler version of the same move would have collected.
GENTLER_ALTERNATIVE: Mapping[str, str] = {"knockout": "knockdown_50"}

#: Every action by name, for resolving an answer back to its definition.
ACTION_CATALOGUE: Mapping[str, Action] = {
    action.name: action
    for action in (
        *ACT_ACTIONS,
        *LOOK_ACTIONS,
        END_ACTION,
        UNDO_ACTION,
        ADOPT_ACTION,
        RESTORE_ACTION,
    )
}


@dataclass(frozen=True)
class Intervention:
    """One gene edit, and every reaction bound that follows from it.

    The agent chooses a *reaction*; what gets built is a *gene* edit; and the two are not the
    same thing, so this object holds both. ``reaction_id`` is what was chosen and ``genes`` is
    the minimal set whose loss achieves it \u2014 for ``ACKr`` that is all three of *ackA*,
    *tdcD* and *purT*, because deleting the textbook one alone leaves two isozymes and changes
    nothing at all. ``bounds`` is then every reaction the edit actually constrains, which for a
    shared gene is more than one: *dctA* runs ``SUCCt2_2``, ``FUMt2_2`` and ``MALt2_2``, and a
    design that edits it edits all three whether anyone intended that or not.

    ``reference_flux`` is the chosen reaction's flux in the round-0 wild-type pFBA state. The
    knockdown is defined against it rather than against the current state, so "half" does not
    drift as a run accumulates interventions.
    """

    reaction_id: str
    mode: InterventionMode
    level: float
    reference_flux: float
    action_name: str
    #: The minimal gene set whose deletion achieves this move. Empty only for a reaction with
    #: no gene association, which the board does not offer.
    genes: tuple[str, ...] = ()
    #: ``(reaction id, lower, upper)`` for every reaction the edit constrains, the chosen one
    #: included. This is what the engine applies; applying it to the chosen reaction alone
    #: would model a strain nobody can build.
    bounds: tuple[tuple[str, float, float], ...] = ()
    #: Reactions the gene edit also hits but whose share of the change cannot be expressed \u2014
    #: a knockdown of a reaction carrying no wild-type flux has nothing to be half of. Named
    #: rather than silently forced to zero, which would turn a knockdown into a knockout of
    #: something the agent never chose.
    unmodelled: tuple[str, ...] = ()
    #: Readable gene names where the model carries them, so a brief naming *ldhA* and a board
    #: naming ``b1380`` are connectable.
    gene_names: tuple[str, ...] = ()

    @property
    def lower_bound(self) -> float:
        """The chosen reaction's new lower bound."""

        return next(
            (low for rid, low, _ in self.bounds if rid == self.reaction_id), 0.0
        )

    @property
    def upper_bound(self) -> float:
        """The chosen reaction's new upper bound."""

        return next(
            (high for rid, _, high in self.bounds if rid == self.reaction_id), 0.0
        )

    @property
    def side_effects(self) -> tuple[str, ...]:
        """Reactions constrained that the agent did not choose."""

        return tuple(rid for rid, _, _ in self.bounds if rid != self.reaction_id)

    def to_reaction_bounds(self) -> tuple[ReactionBound, ...]:
        """The change as the :class:`~cmm.core.condition.ReactionBound` CMM already applies."""

        return tuple(
            ReactionBound(reaction_id=rid, lower_bound=low, upper_bound=high)
            for rid, low, high in self.bounds
        )

    def gene_phrase(self) -> str:
        """How the edit is named to a reader: gene names where the model has them."""

        if not self.genes:
            return "no gene association"
        return ", ".join(self.gene_names or self.genes)

    def describe(self) -> str:
        """A one-line human summary, used in the state history and in report tables."""

        head = (
            self.reaction_id
            if not self.genes
            else f"{self.reaction_id} ({self.gene_phrase()})"
        )
        if self.mode == "knockout":
            line = f"{head}: delete the gene(s), flux forced to 0"
        else:
            line = (
                f"{head}: weaken the gene(s) to 50% of wild type "
                f"(|v| \u2264 {abs(self.level * self.reference_flux):.4g})"
            )
        if self.side_effects:
            line += f" \u2014 also constrains {', '.join(self.side_effects)}"
        if self.unmodelled:
            line += (
                f"; the edit also touches {', '.join(self.unmodelled)}, which carry no "
                "wild-type flux, so the weakening cannot be expressed on them"
            )
        return line

    def to_record(self) -> dict[str, object]:
        """Flat export row."""

        return {
            "reaction_id": self.reaction_id,
            "genes": ";".join(self.genes),
            "gene_names": ";".join(self.gene_names),
            "mode": self.mode,
            "action": self.action_name,
            "level": self.level,
            "reference_flux": self.reference_flux,
            "lower_bound": self.lower_bound,
            "upper_bound": self.upper_bound,
            "reactions_constrained": ";".join(rid for rid, _, _ in self.bounds),
            "side_effects": ";".join(self.side_effects),
        }


class ActionNotApplicable(ValueError):
    """Raised when a chosen move cannot be expressed on the chosen reaction.

    This is not a failure of the run. It is the honest answer to "halve a flux that is
    already zero", and the engine records the reason and lets JEV pick again rather than
    substituting a different move.
    """


def build_intervention(
    model: Model,
    reaction_id: str,
    action: Action,
    reference_fluxes: Mapping[str, float],
) -> Intervention:
    """Turn a chosen ACT move into the gene edit that achieves it, and its bounds.

    The move names a reaction; what is applied is the consequence of deleting or weakening the
    minimal gene set behind it, which is frequently more than that one reaction. Resolving it
    here rather than at report time is what makes CMM's viability rule a rule about strains
    that can be built: a deletion whose collateral kills the cell is reverted by the engine
    like any other unviable move, instead of surviving to the end and failing in a flask.

    A knockdown needs a non-zero wild-type flux to be half *of*, and raises
    :class:`ActionNotApplicable` without one rather than quietly becoming something else. A
    knockout is always defined, including on a zero-flux reaction \u2014 where it changes
    nothing today, which the engine's own re-solve will show, but closes a route the cell
    could otherwise escape down once another deletion lands.
    """

    if action.mode is None:
        raise ActionNotApplicable(
            f"{action.name!r} is a {action.kind} move and does not change a reaction"
        )
    edit = resolve_gene_edit(model, reaction_id)
    names = gene_names(model, edit.genes)
    reference = float(reference_fluxes.get(reaction_id, 0.0))
    level = float(action.level if action.level is not None else 0.0)
    readable = tuple(names.get(gene, gene) for gene in edit.genes)

    if action.mode == "knockout":
        return Intervention(
            reaction_id=reaction_id,
            mode="knockout",
            level=0.0,
            reference_flux=reference,
            action_name=action.name,
            genes=edit.genes,
            gene_names=readable,
            bounds=tuple((rid, 0.0, 0.0) for rid in edit.blocks),
        )

    if abs(reference) <= FLUX_EPSILON:
        raise ActionNotApplicable(
            f"{action.name!r} is relative to the wild-type flux of {reaction_id!r}, "
            f"which is {reference:.3g}: there is nothing to halve"
        )

    # Weakening the gene set weakens every reaction that set alone carries. The cap is each
    # reaction's own wild-type flux, because that is the quantity the level is a fraction of.
    bounds: list[tuple[str, float, float]] = []
    unmodelled: list[str] = []
    for rid in edit.blocks:
        affected = model.reactions.get_by_id(rid)
        own = float(reference_fluxes.get(rid, 0.0))
        if abs(own) <= FLUX_EPSILON:
            # Half of nothing is nothing, and writing zero here would knock the reaction out
            # rather than weaken it. Say so instead of doing the wrong thing quietly.
            if rid != reaction_id:
                unmodelled.append(rid)
            continue
        cap = abs(level * own)
        # Cap the magnitude without opening a direction the reaction did not already have.
        lower = max(float(affected.lower_bound), -cap)
        upper = min(float(affected.upper_bound), cap)
        if lower > upper:
            # The reaction is already held above the cap — a maintenance demand pinned at a
            # floor, or a bound an earlier edit tightened. Weakening it to half is not a
            # constraint that exists, and writing the crossed bounds would raise from inside
            # cobra with nothing to say which move caused it.
            raise ActionNotApplicable(
                f"{action.name!r} would cap {rid!r} at {cap:.4g}, below the lower bound of "
                f"{affected.lower_bound:.4g} it is already held at; the reaction cannot carry "
                "half its flux"
            )
        bounds.append((rid, lower, upper))

    return Intervention(
        reaction_id=reaction_id,
        mode="knockdown",
        level=level,
        reference_flux=reference,
        action_name=action.name,
        genes=edit.genes,
        gene_names=readable,
        bounds=tuple(bounds),
        unmodelled=tuple(unmodelled),
    )


def applicable_actions(reference_flux: float) -> tuple[Action, ...]:
    """The ACT moves that are defined for a reaction with this wild-type flux.

    A knockout is always defined. A knockdown is not defined for a reaction carrying nothing,
    because there is no flux for it to be half of. Offering a move that cannot be applied
    would waste a step and teach the agent nothing, so the list is filtered before the
    question is built rather than after it is answered.
    """

    if abs(float(reference_flux)) <= FLUX_EPSILON:
        return tuple(action for action in ACT_ACTIONS if action.mode == "knockout")
    return ACT_ACTIONS
