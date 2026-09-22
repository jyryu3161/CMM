"""How JEV is asked. The wording here is the scientific content of the feature.

JEV cannot say anything that is not offered to it, so the quality of a run is decided almost
entirely by two things: what the state shows (:mod:`cmm.jev.state`) and how the question is
put. That makes the wording a first-class, versioned artifact rather than a string literal
buried in the loop — a run records which :class:`QuestionSet` produced it, so a change to the
phrasing is a change that can be compared against earlier runs instead of silently
reinterpreting them.

One tick is two calls:

1. **target** — a ``choice`` over the shortlisted reactions plus ``undo_last`` and
   ``end_round``. The answer's ``probabilities`` cover every offered key, so this single call
   also yields a complete ranking of the board, which the engine records.
2. **action** — a ``choice`` over the moves that are actually applicable to the chosen
   reaction, alongside a ``score`` for the expected gain and a ``noul`` for the growth risk.
   "Applicable" now also means "within budget": knockouts and knockdowns are limited
   separately, so a move whose own budget is spent is never put on the list.

Two calls rather than one because the product of the two vocabularies (every reaction times
every move) would be a criteria set too large to read and too large for a 32K context.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass

from cmm.jev.actions import (
    ADOPT_ACTION,
    END_ACTION,
    RESTORE_ACTION,
    LOOK_ACTIONS,
    UNDO_ACTION,
    Action,
    applicable_actions,
)
from cmm.jev.state import CandidateEvidence
from cmm.jev._transport import choice_question, noul_question, score_question


class NoAvailableAction(RuntimeError):
    """Raised when a candidate has no move left that has not already been tried.

    The engine treats this as a reason to drop the reaction from the board, not as an error:
    a reaction whose every move has failed is a reaction there is nothing more to learn from.
    """


#: Answer keys. Fixed names so the engine, the transcript and the tests agree.
TARGET_KEY = "target"
ACTION_KEY = "action"
BENEFIT_KEY = "benefit"
RISK_KEY = "growth_risk"


def available_actions(
    candidate: CandidateEvidence,
    *,
    allow_look: bool,
    exclude: Collection[str] = (),
    allowed_modes: Collection[str] = ("knockout", "knockdown"),
) -> tuple[Action, ...]:
    """Every move that could actually be executed on this candidate right now.

    Four things narrow the list, and each one exists because leaving it out cost real steps:

    * a move undefined for the reaction (a knockdown of a flux that is already zero);
    * a move already tried and rejected on this reaction, in ``exclude``;
    * a move whose budget is spent — ``allowed_modes`` is how the separate knockout and
      knockdown limits reach the question, so a design that may take no further deletion is
      never offered one;
    * a scan whose answer is already in the record, which would spend a step to learn nothing.

    The engine uses this directly as well as through the question: when exactly one move is
    left there is nothing to decide, so it is taken without paying for a second call.
    """

    modes = set(allowed_modes)
    blocked = set(exclude)
    actions: list[Action] = [
        action
        for action in applicable_actions(candidate.reference_flux)
        if action.name not in blocked and action.mode in modes
    ]
    if allow_look:
        actions.extend(
            action
            for action in LOOK_ACTIONS
            if action.name not in blocked
            and not (
                action.name == "essentiality_scan" and candidate.essential is not None
            )
            and not (action.name == "fseof_scan" and candidate.fseof_slope is not None)
        )
    return tuple(actions)


@dataclass(frozen=True)
class QuestionSet:
    """A named, versioned way of asking. ``version`` goes into every run's provenance."""

    version: str
    description: str

    def target_question(
        self,
        candidates: Sequence[CandidateEvidence],
        *,
        product: str,
        growth_floor: float,
        allow_undo: bool,
        allow_look: bool,
        design_full: bool = False,
        proven_design: str = "",
        best_design: str = "",
        rounds_found: int = 0,
    ) -> dict[str, Mapping[str, object]]:
        """Stage one: which reaction to act on next, or stop, undo, adopt or go back."""

        raise NotImplementedError

    def action_question(
        self,
        candidate: CandidateEvidence,
        *,
        product: str,
        growth_floor: float,
        allow_look: bool,
        exclude: Collection[str] = (),
        allowed_modes: Collection[str] = ("knockout", "knockdown"),
    ) -> dict[str, Mapping[str, object]]:
        """Stage two: what to do to the reaction stage one chose."""

        raise NotImplementedError


class ProductionV1(QuestionSet):
    """The first production-maximisation question set.

    Design choices worth stating, because each one was a real fork:

    * The target question asks for **one** reaction, not a ranked list. JEV returns a
      probability for every offered key regardless, so the ranking arrives anyway and the
      question stays a question a person could answer.
    * ``end_round`` and ``undo_last`` sit in the same vocabulary as the reactions rather than
      in a separate "should you continue?" call. Stopping is a move, and making it one costs
      nothing and removes a whole call from every tick.
    * The benefit ``score`` is asked next to the action ``choice`` in the same call, so the
      grade applies to the move that was actually chosen.
    * The growth risk is asked as a ``noul`` for the record, but the engine **does not act on
      it**: whether growth survives is something CMM solves for, and a prediction is no
      substitute for the solve. It is kept because the gap between what JEV expected and what
      the solver found is the interesting quantity in a post-hoc analysis.
    """

    def __init__(self) -> None:
        super().__init__(
            version="production_v1",
            description=(
                "Two-stage tick for raising a product exchange flux under a growth floor: "
                "choose a reaction, then choose what to do to it."
            ),
        )

    def target_question(
        self,
        candidates: Sequence[CandidateEvidence],
        *,
        product: str,
        growth_floor: float,
        allow_undo: bool,
        allow_look: bool,
        design_full: bool = False,
        proven_design: str = "",
        best_design: str = "",
        rounds_found: int = 0,
    ) -> dict[str, Mapping[str, object]]:
        if not candidates and not design_full:
            raise ValueError(
                "the target question needs at least one candidate reaction"
            )

        criteria: dict[str, str] = {}
        if not design_full:
            # A full design cannot take another intervention, so listing the reactions would
            # be offering a move that is certain to be refused. Observed without this guard:
            # the same knockout proposed on ten consecutive ticks, rejected every time for
            # want of room, until the run ended.
            criteria.update(
                {
                    candidate.reaction_id: candidate.to_label()
                    for candidate in candidates
                }
            )
        if proven_design and not design_full:
            criteria[ADOPT_ACTION.name] = (
                f"{ADOPT_ACTION.description} The design available is: {proven_design}."
            )
        if best_design:
            criteria[RESTORE_ACTION.name] = (
                f"{RESTORE_ACTION.description} The best design so far: {best_design}."
            )
        if allow_undo:
            criteria[UNDO_ACTION.name] = UNDO_ACTION.description
        criteria[END_ACTION.name] = END_ACTION.description

        if design_full:
            return {
                TARGET_KEY: choice_question(
                    (
                        "The design already carries the maximum number of interventions, so "
                        "nothing more can be added. Either withdraw the most recent "
                        "intervention to make room for a different one, or end the round and "
                        "let the design be measured as it stands. The product being "
                        f"maximised is {product}."
                    ),
                    criteria,
                )
            }

        instructions = (
            f"Choose the single reaction to act on next so that flux through {product} "
            f"increases, while the growth rate stays at or above {growth_floor} per hour. "
            "The only moves available are deleting the reaction's genes or weakening them to "
            "half their wild-type activity; nothing can be over-expressed, so look for flux "
            "that competes with the product rather than flux that feeds it. "
            "The evidence for each option is the record with the same id in the state "
            "above: its wild-type flux, its current flux, how many steps it sits from the "
            "product, what it does to the ATP and redox pools, and what editing it costs. "
            "That last part matters: a move is made on genes, not on reactions, so a record "
            "saying the edit also stops other reactions means those will be stopped too, and "
            "a record naming three isozymes means all three have to be deleted for the move "
            "to do anything at all. "
            "A record beginning MEASURED reports what CMM found when it actually made the "
            "move on the design as it stands, and outranks every other line, including your "
            "own reasoning about the stoichiometry. "
            "A record saying OptKnock or RobustKnock deletes the reaction is the next "
            "strongest: those are proofs that the deletion forces the product at maximum "
            "growth, and several of them name reactions carrying no flux today, which the "
            "cell would switch to once the obvious routes are shut. "
        )
        if allow_look:
            instructions += (
                "If the records do not yet say what you need, you may choose a reaction and "
                "then ask for a scan instead of intervening. "
            )
        if rounds_found:
            instructions += (
                "Earlier rounds of this run already found the designs listed under "
                "designs_already_found, and each round's entry under previous_rounds says "
                "what it left undone. Rebuilding a design that is already there adds nothing "
                "to the run: prefer a route those rounds did not take, or go after something "
                "one of them explicitly left on the table. Returning to the best design so "
                "far is still the right move when the alternatives have been tried and were "
                "worse \u2014 but say so by choosing restore_best_design rather than by "
                "reassembling it move by move. "
            )
        instructions += (
            f"Choose {END_ACTION.name} only when no remaining move is expected to help."
        )
        return {TARGET_KEY: choice_question(instructions, criteria)}

    def action_question(
        self,
        candidate: CandidateEvidence,
        *,
        product: str,
        growth_floor: float,
        allow_look: bool,
        exclude: Collection[str] = (),
        allowed_modes: Collection[str] = ("knockout", "knockdown"),
    ) -> dict[str, Mapping[str, object]]:
        actions = available_actions(
            candidate,
            allow_look=allow_look,
            exclude=exclude,
            allowed_modes=allowed_modes,
        )
        criteria = {action.name: action.description for action in actions}
        if len(criteria) < 2:
            # One option is not a question. The engine takes a lone move without asking; it
            # only reaches here with nothing left at all.
            raise NoAvailableAction(
                f"every move on {candidate.reaction_id!r} has already been tried or is "
                "undefined for it"
            )

        rid = candidate.reaction_id
        missing = []
        if candidate.essential is None:
            missing.append("whether it is essential")
        if candidate.fseof_slope is None:
            missing.append("whether its flux falls as the product is forced up")
        gap = (
            f" What is still unknown about it: {', and '.join(missing)}."
            if missing and allow_look
            else ""
        )

        return {
            ACTION_KEY: choice_question(
                (
                    f'Choose what to do to the genes behind reaction "{rid}" so that flux '
                    f"through {product} increases while growth stays at or above "
                    f"{growth_floor} per hour. Its evidence is the record with id {rid} in "
                    f"the state above; a MEASURED line there is what CMM found when it made "
                    f"that exact move on the current design.{gap}"
                ),
                criteria,
            ),
            BENEFIT_KEY: score_question(
                (
                    f"Rate how much you expect flux through {product} to increase as a "
                    f'result of the move you just chose for "{rid}". If its record reports '
                    "a MEASURED change for that move, grade against that number rather "
                    "than against what the stoichiometry suggests."
                ),
                [
                    "No gain at all, or the cell stops growing.",
                    "A marginal gain, or a gain you are not confident in.",
                    "A moderate gain, by a mechanism you can point to in the record.",
                    "A large gain; the record makes the mechanism explicit.",
                    "A decisive gain; this step is what is limiting the product.",
                ],
            ),
            RISK_KEY: noul_question(
                (
                    f'Will this move on "{rid}" push the growth rate below '
                    f"{growth_floor} per hour?"
                ),
                true_means=(
                    "Growth very likely falls below the floor, making the strain unusable."
                ),
                false_means="Growth very likely stays at or above the floor.",
            ),
        }


#: Every shipped question set, by version. A config names one; provenance records it.
QUESTION_SETS: Mapping[str, QuestionSet] = {"production_v1": ProductionV1()}

DEFAULT_QUESTION_SET = "production_v1"


def get_question_set(version: str = DEFAULT_QUESTION_SET) -> QuestionSet:
    """Look a question set up by version, naming the alternatives when it is not found."""

    try:
        return QUESTION_SETS[version]
    except KeyError:
        known = ", ".join(sorted(QUESTION_SETS))
        raise ValueError(
            f"unknown JEV question set {version!r}; available sets: {known}"
        ) from None


def research_query(candidate: CandidateEvidence, product: str, organism: str) -> str:
    """The web-research prompt for one candidate.

    Deliberately narrow, and deliberately two-sided. An early version asked only what raises
    the product, and a source that answers only that question is worse than none: the model
    already predicts the yield, and what it cannot predict is the fitness cost, the
    regulatory response or the byproduct that a paper would report. So the prompt asks for
    the penalty as explicitly as the benefit, and asks for silence to be reported as silence
    rather than filled in.
    """

    gene_hint = f" (genes {', '.join(candidate.genes)})" if candidate.genes else ""
    return (
        f"In metabolic engineering of {organism}, what has been published about modifying "
        f"the reaction {candidate.reaction_id}{gene_hint} \u2014 {candidate.name} \u2014 to "
        f"increase production of {product}?\n"
        "Answer these three things separately and briefly:\n"
        "1. Has knockout, knockdown or overexpression of this step been reported, and what "
        "happened to the product?\n"
        "2. What was the cost? State any reported growth defect, fitness burden, reduced "
        "biomass yield, or dependence on a supplement or a specific medium.\n"
        "3. Were there side effects a flux model would not predict \u2014 byproduct "
        "accumulation, regulatory compensation, protein burden, toxicity?\n"
        "If a point has no published work behind it, say so for that point rather than "
        "generalising from a related enzyme or a different organism."
    )
