"""The game loop: JEV plays CMM, one tick at a time.

Each tick is a frame. CMM renders the current metabolic state, JEV picks a reaction and a
move, CMM executes the move and re-solves, and the new flux distribution becomes the next
frame. A round is a run of ticks with a checkpoint at the end; a run is a sequence of rounds.

The division of labour is deliberate and is the reason the loop is trustworthy:

**CMM owns the rules.** Every number on the screen is a solve, not a claim. Viability is
enforced here, not predicted: a move that makes the model infeasible or drops growth below
the configured floor is undone by CMM regardless of how confident JEV was, and the reason is
recorded.

**JEV owns the strategy.** A move that is legal but unhelpful is kept, flagged in the history
JEV reads on the next tick, and left for JEV to withdraw with ``undo_last``. Auto-reverting
every move that failed to improve the product would turn the loop into greedy hill climbing
and make a two-step manoeuvre — sacrifice a branch now, gain from it after the next move —
impossible to express.

Two flux states are computed after every intervention, because they answer different
questions and a report that quoted only one would mislead:

``pFBA``
    The re-optimised state. What the strain does once it has adapted to the change. This is
    the score, because a design that only pays off if the cell chooses to make the product is
    not a design.
``MOMA``
    The minimal-adjustment state against the wild-type reference. What the strain does
    immediately after the change, before any adaptation.

**Which state MOMA measures from is a choice, and it is the wild type here.** MOMA's premise
is that a freshly perturbed cell keeps the regulatory setpoints of the cell it was made from,
so the reference must be the parent strain — and the parent of the design this run produces
is the wild type, because the design is built and characterised as one strain, not as a
sequence of strains. The intermediate designs the agent passes through on its way there are
search positions, not organisms anyone will culture.

The alternative — re-referencing each step to the previous design's flux state — models a
different experiment: an edit introduced into a strain that has already been grown up. It is
a real protocol, but it is the wrong one to report here for two reasons. Chaining the
reference makes each step's MOMA distance describe only the last edit, so a five-edit design
would look exactly as easy to build as a one-edit design, which is the opposite of what the
number is for. And it would not be MOMA from the previous *MOMA* state in any case: a strain
you can make a second edit in is a strain you have cultured, and cultured knockout strains
move toward the FBA optimum (Fong & Palsson 2004), so the honest parent state would be the
previous design's pFBA — which is what re-referencing would have to use, and what it would
then be comparing everything against instead of the organism.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import pandas as pd
from cobra import Model
from cobra.exceptions import OptimizationError

from cmm.core.condition import Condition
from cmm.core.flux_state import FluxState
from cmm.core.media import Medium
from cmm.core.provenance import run_provenance
from cmm.core.simulation import FluxSolution, pfba
from cmm.core.solvers import solver_status, supports
from cmm.jev._transport import DecisionResult, JevClient, JevTransportError
from cmm.jev.actions import (
    ACTION_CATALOGUE,
    ADOPT_ACTION,
    END_ACTION,
    GENTLER_ALTERNATIVE,
    RESTORE_ACTION,
    UNDO_ACTION,
    ActionNotApplicable,
    Intervention,
    applicable_actions,
    build_intervention,
)
from cmm.jev.questions import (
    ACTION_KEY,
    BENEFIT_KEY,
    DEFAULT_QUESTION_SET,
    RISK_KEY,
    TARGET_KEY,
    available_actions,
    get_question_set,
    literature_briefing,
)

if TYPE_CHECKING:  # the benchmark imports this module, so the type is compile-time only
    from cmm.jev.benchmark import BaselineRow

from cmm.jev.state import (
    CandidateEvidence,
    GameState,
    ScanCache,
    build_candidates,
    cofactor_balance,
    cofactor_limitation,
    guaranteed_product,
    resolve_pools,
)

#: How many times the identical move may land the identical way before the round is cut
#: short. Three is enough to be sure it is a loop and cheap enough not to matter if it is not.
_MAX_IDENTICAL_MOVES = 3

#: How many consecutive rounds may end without the agent changing anything before the run
#: stops. Every round now starts from the same wild type, so a round that applies nothing at
#: all means the agent has stopped finding the board worth acting on; two in a row is enough
#: to believe it. A round that *tries* and is refused is not idle.
_MAX_IDLE_ROUNDS = 2

TickOutcome = Literal[
    "applied",
    "undone_by_agent",
    "reverted_infeasible",
    "reverted_growth_floor",
    "scan",
    "end_round",
    "not_applicable",
    "budget_exhausted",
]

#: Reasons a run can stop before playing every round, for the summary to state plainly.
STOP_REQUESTED = "the run was stopped from the interface"


class JevWorkflowError(RuntimeError):
    """Raised when the run cannot proceed on scientific or configuration grounds."""


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JevConfig:
    """A complete, serializable invocation of a JEV design run.

    **A round is one independent attempt**, not a phase of a longer one. Every round starts
    from the wild type with an empty design, plays until its steps run out or the agent ends
    it, and is scored on its own. What carries across is knowledge, not bounds: the agent is
    shown what each earlier round reached and with which design, so a later round can go after
    something different or head for what already worked.

    Rounds used to continue one another, and the effect was not subtle — the first round
    filled the design and the rest had nothing left to do, so a three-round run spent five
    steps of a possible thirty-six.

    Two kinds of budget, and they mean different things. ``steps_per_round`` bounds how long
    one attempt may *play*: every decision costs one step, including an undo and including a
    scan that changes nothing. ``max_knockouts`` and ``max_knockdowns`` bound how many edits
    of each kind an attempt may carry at once, which is the quantity a wet-lab reader cares
    about — a design needing twelve edits is not the same proposal as one needing three, and
    six deletions is a different project from six promoter swaps.

    They are counted separately rather than as one total because they cost different things
    to build, and because one shared cap starved the run: the seeded OptKnock design takes
    three deletions on its own, so a total of four left exactly one edit for the agent and
    every round ended within a step or two of adopting it.
    """

    model_path: str | Path
    product: str
    output_dir: str | Path | None = None
    #: The exchange the yield is quoted per. Leave it unset: the substrate is whichever
    #: carbon source the model is actually taking up, which the wild-type solve already says,
    #: and the medium or condition is what decides that. Naming it separately is a second
    #: place for the same fact to be wrong.
    substrate: str | None = None
    biomass: str | None = None
    solver: str | None = None
    medium: Medium | str | None = None
    condition: Condition | None = None
    #: The organism the model represents. Only used for the literature lookup, and required
    #: when that is on: a default would ask the published record about the wrong species and
    #: return an answer that is confident and wrong. Left empty when no lookup is done,
    #: because guessing it from a model id is not a thing that works.
    organism: str = ""
    #: What the person running this knows and the model does not: published targets for this
    #: product, a growth rate the strain has to hold, a cofactor they believe is decisive, a
    #: reaction they want left alone. Shown to the agent under ``your_brief`` on every step.
    #:
    #: It is **guidance, not permission**. It cannot widen the move vocabulary, name a
    #: reaction outside the model, or lift the growth floor CMM enforces — the agent still
    #: answers only with the criteria this package supplies, so the worst a mistaken brief can
    #: do is waste steps.
    brief: str = ""

    # -- the game -----------------------------------------------------------
    rounds: int = 5
    #: Steps the agent may spend in one round. **Every** decision costs one — an intervention,
    #: an undo, a scan — so this is the length of the game, not a count of edits. The design
    #: is bounded separately by ``max_knockouts`` and ``max_knockdowns``, which are the
    #: numbers a laboratory would have to build. A long round is cheap: a step is at most two
    #: calls, about 0.6 s and $0.00016.
    steps_per_round: int = 40
    #: How many gene deletions the design may carry, and how many 50% knockdowns. The
    #: deletion budget has to clear the seeded strain design (three deletions by default)
    #: with room to spare, or the agent inherits a full design and has nothing to play.
    max_knockouts: int = 6
    max_knockdowns: int = 3
    growth_floor: float = 0.05
    candidate_limit: int = 24
    allow_look_actions: bool = True
    run_moma: bool = True
    #: Bounds on the ``strain_design_scan`` LOOK move. It runs the same OptKnock and
    #: RobustKnock services the SC-01 workflow uses, so its cost is theirs.
    design_max_knockouts: int = 3
    design_max_solutions: int = 5
    #: Run the deterministic strain designer once before the first move and put the reactions
    #: it names on the board from the start. On anaerobic succinate the designer's answer
    #: deletes reactions carrying no flux at all, which no flux-derived board can surface, so
    #: without this the agent cannot reach the known optimum however long it plays. Seeding it
    #: is not cheating: the agent is being handed CMM's best deterministic result and asked to
    #: improve on it, which is the only comparison that means anything.
    seed_with_strain_design: bool = True
    #: Require each round to reach a design no earlier round already found, by withholding one
    #: reaction from each design already in hand — the integer cut OptKnock itself uses to
    #: enumerate alternative designs.
    #:
    #: Telling the agent that a design had been found was measured and was not enough: six
    #: rounds produced two distinct designs and four exact repeats, because every round starts
    #: from the same wild type and sees the same board, so it plays the same game. A run asked
    #: for several attempts should return several answers.
    #:
    #: It never loses the best design: the global best is tracked across rounds and reported
    #: whichever round found it. What a later round cannot do is find it again.
    require_distinct_rounds: bool = True
    #: After the game, score the agent's design against the deterministic methods on the same
    #: problem, every design evaluated the same way. An agent result with nothing to compare
    #: it to is not a result.
    run_baseline_comparison: bool = True
    #: Measure, every tick, how much more product one extra unit of NADH, NADPH or ATP would
    #: buy. Three extra LPs. This is the reading neither MOMA nor OptKnock reports, and it is
    #: usually what decides whether the next move should route carbon or supply a cofactor.
    measure_cofactor_limits: bool = True
    #: Measure, every tick, the worst product the design could give while growing as fast as
    #: it can. A design whose worst case is zero is not a design, however good its pFBA
    #: number looks. Loopless, so it costs one FVA solve.
    measure_guaranteed_product: bool = True
    #: Recompute, whenever the design changes, what deleting and what halving each reaction
    #: on the board would do to the product. Two pFBA solves per candidate — about two
    #: seconds for a board of 24 on ``e_coli_core``, and the single most useful thing on the
    #: screen, because it is the measured consequence of the exact moves the agent may make.
    #: Turn it off for a genome-scale model where a board of solves is not cheap.
    screen_interventions: bool = True

    # -- the agent ----------------------------------------------------------
    jev_model: str = "typesafe/jev-1.13"
    question_set: str = DEFAULT_QUESTION_SET
    enable_web_research: bool = False
    research_model: str = "openai/gpt-5.6-luna"
    #: Genes, reactions or subsystems the run may not touch, whatever the agent thinks. The
    #: brief is *guidance* — the agent reads it and may still disagree with it — and this is
    #: *enforcement*: anything named here is taken off the board before the agent sees it, and
    #: no proven design containing it is offered. Matches a reaction id, a gene id, a gene
    #: name, or a subsystem name, case-insensitively.
    off_limits: tuple[str, ...] = ()
    max_decisions: int = 2000
    max_cost_usd: float = 5.0
    request_timeout_s: float = 60.0
    seed: int = 0

    overwrite: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "model_path", Path(self.model_path))
        if self.output_dir is not None:
            object.__setattr__(self, "output_dir", Path(self.output_dir))
        self.validate()

    def validate(self) -> None:
        if not str(self.model_path):
            raise ValueError("model_path must not be empty")
        if not self.product:
            raise ValueError("product must not be empty")
        if self.rounds < 1:
            raise ValueError("rounds must be at least 1")
        if self.steps_per_round < 1:
            raise ValueError("steps_per_round must be at least 1")
        if self.max_knockouts < 0:
            raise ValueError("max_knockouts must be non-negative")
        if self.max_knockdowns < 0:
            raise ValueError("max_knockdowns must be non-negative")
        if self.max_interventions < 1:
            raise ValueError(
                "max_knockouts and max_knockdowns cannot both be zero: the agent would "
                "have no move it is allowed to make"
            )
        if self.growth_floor < 0:
            raise ValueError("growth_floor must be non-negative")
        if self.candidate_limit < 2:
            raise ValueError("candidate_limit must be at least 2 for a choice question")
        if self.max_decisions < 1:
            raise ValueError("max_decisions must be at least 1")
        if self.max_cost_usd <= 0:
            raise ValueError("max_cost_usd must be positive")

        if self.enable_web_research and not self.organism.strip():
            raise ValueError(
                "enable_web_research needs organism: a literature lookup asks about a named "
                "species, and the wrong one returns evidence that is confident and wrong"
            )
        get_question_set(self.question_set)

    @property
    def max_interventions(self) -> int:
        """Edits the design may carry in total. Derived: the two budgets are the inputs."""

        return self.max_knockouts + self.max_knockdowns

    def room_for(self, used_knockouts: int, used_knockdowns: int) -> tuple[str, ...]:
        """Which intervention modes still have budget, in the order the vocabulary lists them."""

        modes = []
        if used_knockouts < self.max_knockouts:
            modes.append("knockout")
        if used_knockdowns < self.max_knockdowns:
            modes.append("knockdown")
        return tuple(modes)

    @classmethod
    def from_json(cls, path: str | Path) -> "JevConfig":
        """Load from a UTF-8 JSON object, resolving relative paths against the file."""

        config_path = Path(path).expanduser().resolve()
        payload = json.loads(config_path.read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise ValueError("JEV workflow JSON must contain an object")
        values = dict(payload)
        for name in ("model_path", "output_dir"):
            raw = values.get(name)
            if raw is None:
                continue
            candidate = Path(str(raw)).expanduser()
            if not candidate.is_absolute():
                candidate = config_path.parent / candidate
            values[name] = candidate.resolve()
        return cls.from_mapping(values)

    @classmethod
    def from_mapping(cls, payload: Mapping[str, object]) -> "JevConfig":
        from cmm.workflows.production import (
            _condition_from_payload,
            _medium_from_payload,
        )

        values = dict(payload)
        if "max_interventions" in values:
            raise ValueError(
                "max_interventions is no longer a setting: deletions and knockdowns are "
                "budgeted separately. Replace it with max_knockouts and max_knockdowns, "
                "which say what a laboratory would actually have to build."
            )
        if "max_research_calls" in values:
            raise ValueError(
                "max_research_calls is no longer a setting: the run reads the literature "
                "once before the first move rather than once per candidate, because a web "
                "search inside a step made the loop wait tens of seconds for it"
            )
        raw_limits = values.get("off_limits")
        if isinstance(raw_limits, (list, tuple)):
            values["off_limits"] = tuple(str(name) for name in raw_limits)
        elif isinstance(raw_limits, str):
            # A single name, written as a string rather than a one-element list. Accepting it
            # beats failing on the shape of a constraint that is otherwise perfectly clear.
            values["off_limits"] = (raw_limits,)
        if "screen_amplifications" in values:
            raise ValueError(
                "screen_amplifications is no longer a setting: the agent cannot amplify. "
                "Use screen_interventions, which measures deletions and knockdowns instead."
            )
        values["medium"] = _medium_from_payload(values.get("medium"))
        values["condition"] = _condition_from_payload(values.get("condition"))
        return cls(**values)  # type: ignore[arg-type]

    def to_provenance(self) -> dict[str, object]:
        """Every parameter that can change the outcome, flat, for the provenance record."""

        return {
            "product": self.product,
            "substrate": self.substrate,
            "biomass": self.biomass,
            "organism": self.organism,
            "brief": self.brief,
            "rounds": self.rounds,
            "steps_per_round": self.steps_per_round,
            "max_knockouts": self.max_knockouts,
            "max_knockdowns": self.max_knockdowns,
            "max_interventions": self.max_interventions,
            "intervention_vocabulary": "knockout, knockdown_50",
            "growth_floor": self.growth_floor,
            "candidate_limit": self.candidate_limit,
            "allow_look_actions": self.allow_look_actions,
            "design_max_knockouts": self.design_max_knockouts,
            "design_max_solutions": self.design_max_solutions,
            "seed_with_strain_design": self.seed_with_strain_design,
            "require_distinct_rounds": self.require_distinct_rounds,
            "run_baseline_comparison": self.run_baseline_comparison,
            "measure_cofactor_limits": self.measure_cofactor_limits,
            "measure_guaranteed_product": self.measure_guaranteed_product,
            "screen_interventions": self.screen_interventions,
            "run_moma": self.run_moma,
            "jev_model_requested": self.jev_model,
            "question_set": self.question_set,
            "enable_web_research": self.enable_web_research,
            "research_model": self.research_model if self.enable_web_research else None,
            "off_limits": list(self.off_limits),
            "max_decisions": self.max_decisions,
            "max_cost_usd": self.max_cost_usd,
            "seed": self.seed,
            "medium": self.medium if isinstance(self.medium, str) else None,
            "condition": self.condition.name if self.condition is not None else None,
            "agent_determinism": (
                "The CMM solves in this run are deterministic. The JEV decisions are not "
                "guaranteed to be: re-running may choose different moves. The full "
                "request/response transcript is saved so any single run can be audited."
            ),
        }


# ---------------------------------------------------------------------------
# records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TickRecord:
    """One frame of the game: what was shown, what was chosen, what happened."""

    round_index: int
    tick_index: int
    target: str
    target_confidence: float | None
    target_ranking: tuple[tuple[str, float], ...]
    action: str | None
    action_confidence: float | None
    benefit_score: float | None
    predicted_growth_risk: float | None
    outcome: TickOutcome
    reason: str
    intervention: Intervention | None
    product_flux: float
    growth: float
    status: str
    moma_product_flux: float | None = None
    moma_growth: float | None = None
    moma_distance: float | None = None
    n_active_interventions: int = 0
    decision_cost_usd: float = 0.0
    decision_latency_s: float = 0.0
    #: What the agent weighed at the second stage, and by how much. Kept for the same reason
    #: as the first stage's ranking: the runner-up is often the interesting row.
    action_ranking: tuple[tuple[str, float], ...] = ()

    def to_row(self) -> dict[str, object]:
        return {
            "round": self.round_index,
            "tick": self.tick_index,
            "target": self.target,
            "target_confidence": self.target_confidence,
            "action": self.action,
            "action_confidence": self.action_confidence,
            "benefit_score": self.benefit_score,
            "predicted_growth_risk": self.predicted_growth_risk,
            "outcome": self.outcome,
            "reason": self.reason,
            "reaction_id": (
                self.intervention.reaction_id if self.intervention else None
            ),
            "mode": self.intervention.mode if self.intervention else None,
            "status": self.status,
            "product_flux": self.product_flux,
            "growth": self.growth,
            "moma_product_flux": self.moma_product_flux,
            "moma_growth": self.moma_growth,
            "moma_distance": self.moma_distance,
            "n_active_interventions": self.n_active_interventions,
            "decision_cost_usd": self.decision_cost_usd,
            "decision_latency_s": self.decision_latency_s,
        }

    def headline(self) -> str:
        """The one line the GUI status bar and the history show for this tick."""

        if self.outcome == "end_round":
            return f"R{self.round_index}T{self.tick_index} JEV ended the round"
        if self.outcome == "scan":
            return f"R{self.round_index}T{self.tick_index} ran {self.action} on {self.target}"
        if self.outcome == "undone_by_agent":
            return f"R{self.round_index}T{self.tick_index} JEV withdrew {self.reason}"
        if self.outcome.startswith("reverted"):
            return (
                f"R{self.round_index}T{self.tick_index} {self.action} on {self.target} "
                f"reverted: {self.reason}"
            )
        if self.outcome == "not_applicable":
            return (
                f"R{self.round_index}T{self.tick_index} {self.action} on {self.target} "
                f"not applicable: {self.reason}"
            )
        return (
            f"R{self.round_index}T{self.tick_index} {self.action} on {self.target} "
            f"→ product {self.product_flux:.4g}, growth {self.growth:.4g}"
        )


@dataclass(frozen=True)
class RoundRecord:
    """The checkpoint at the end of a round: what it reached, and what it did not.

    A round that only records its score teaches the next round nothing. The shortfall is the
    useful half — the moves measured to still pay that this round never took, the cofactor
    still limiting the product when it stopped, the budget it left unspent — because that is
    what a later round can act on.
    """

    round_index: int
    n_ticks: int
    product_flux: float
    growth: float
    interventions: tuple[str, ...]
    improved: bool
    ended_early: bool
    #: Why the round stopped: the agent ended it, the steps ran out, nothing was left to try.
    stopped_because: str = ""
    #: What this round left on the table, one clause each, in the order they matter.
    shortfall: tuple[str, ...] = ()
    #: The design as an order-independent key, so a round that rediscovers an earlier one can
    #: be told apart from one that found something new.
    signature: tuple[str, ...] = ()
    #: The earlier round this one reproduced, if any.
    repeated: int | None = None
    #: Reactions this round was not allowed to use, and therefore the question it answers:
    #: "the best design that does not use these". A portfolio a laboratory can choose from
    #: needs that question answered, not six attempts at the same one.
    withheld: tuple[str, ...] = ()

    def to_row(self) -> dict[str, object]:
        return {
            "round": self.round_index,
            "n_ticks": self.n_ticks,
            "product_flux": self.product_flux,
            "growth": self.growth,
            "n_interventions": len(self.interventions),
            "interventions": "; ".join(self.interventions),
            "improved": self.improved,
            "ended_early": self.ended_early,
            "stopped_because": self.stopped_because,
            "shortfall": "; ".join(self.shortfall),
            "design_signature": "; ".join(self.signature),
            "repeated_round": self.repeated,
            "withheld": "; ".join(self.withheld),
            "question_answered": (
                "best design available"
                if not self.withheld
                else "best design that does not use " + ", ".join(self.withheld)
            ),
        }


@dataclass(frozen=True)
class JevResult:
    """Everything one JEV run produced."""

    config: JevConfig
    provenance: Mapping[str, object]
    wild_type_product_flux: float
    wild_type_growth: float
    theoretical_max_yield: float | None
    ticks: tuple[TickRecord, ...]
    rounds: tuple[RoundRecord, ...]
    best_product_flux: float
    best_growth: float
    best_interventions: tuple[Intervention, ...]
    final_interventions: tuple[Intervention, ...]
    flux_frames: tuple[Mapping[str, float], ...] = ()
    transcript: tuple[Mapping[str, object], ...] = ()
    baselines: tuple["BaselineRow", ...] = ()
    #: What the one web search returned, in full, and where it came from. The agent was shown
    #: this on every step, so a reader checking a design needs to be able to read it too.
    literature_brief: str = ""
    literature_sources: tuple[str, ...] = ()
    #: The last evidence the board held for every reaction the run measured, which is what
    #: :func:`~cmm.jev.targets.build_target_reports` reads to state each target's case.
    candidates_seen: tuple[CandidateEvidence, ...] = ()
    notes: tuple[str, ...] = ()
    usage: Mapping[str, object] = field(default_factory=dict)
    run_directory: Path | None = None

    def summary(self) -> dict[str, object]:
        improvement = self.best_product_flux - self.wild_type_product_flux
        return {
            "product": self.config.product,
            "wild_type_product_flux": self.wild_type_product_flux,
            "wild_type_growth": self.wild_type_growth,
            "best_product_flux": self.best_product_flux,
            "best_growth": self.best_growth,
            "absolute_improvement": improvement,
            "fold_improvement": (
                self.best_product_flux / self.wild_type_product_flux
                if abs(self.wild_type_product_flux) > 1e-9
                else None
            ),
            "beat_wild_type": improvement > 1e-9,
            "n_rounds": len(self.rounds),
            "n_ticks": len(self.ticks),
            "n_best_interventions": len(self.best_interventions),
            "best_design": [i.describe() for i in self.best_interventions],
            "question_set": self.config.question_set,
            "baseline_comparison": self.baseline_summary(),
            "targets": self.targets_summary(),
            "notes": list(self.notes),
            "usage": dict(self.usage),
        }

    def targets_summary(self) -> dict[str, object]:
        """How many targets the run weighed, and how many of them pay."""

        from cmm.jev.targets import targets_summary

        return dict(targets_summary(self.targets()))

    def baseline_summary(self) -> dict[str, object] | None:
        """How the agent's design scored against the deterministic methods, or None."""

        if not self.baselines:
            return None
        from cmm.jev.benchmark import comparison_summary

        return comparison_summary(self.baselines, product=self.config.product)

    def baselines_frame(self) -> pd.DataFrame:
        from cmm.jev.benchmark import comparison_frame

        return comparison_frame(self.baselines)

    def literature_frame(self) -> pd.DataFrame:
        """What the web search returned, in full, with its sources."""

        rows = (
            [
                {
                    "product": self.config.product,
                    "organism": self.config.organism,
                    "briefing": self.literature_brief,
                    "n_sources": len(self.literature_sources),
                    "sources": "; ".join(self.literature_sources),
                }
            ]
            if self.literature_brief
            else []
        )
        return pd.DataFrame(
            rows,
            columns=["product", "organism", "briefing", "n_sources", "sources"],
        )

    def targets(self):
        """Every target the run touched, with the case for and against each one."""

        from cmm.jev.targets import build_target_reports

        return build_target_reports(self)

    def targets_frame(self) -> pd.DataFrame:
        """The per-target pros-and-cons table."""

        from cmm.jev.targets import targets_frame

        return targets_frame(self.targets())

    def ticks_frame(self) -> pd.DataFrame:
        return pd.DataFrame([tick.to_row() for tick in self.ticks])

    def rounds_frame(self) -> pd.DataFrame:
        return pd.DataFrame([record.to_row() for record in self.rounds])

    def interventions_frame(self) -> pd.DataFrame:
        return pd.DataFrame([i.to_record() for i in self.best_interventions])


# ---------------------------------------------------------------------------
# the loop
# ---------------------------------------------------------------------------


@dataclass
class _Board:
    """Mutable run state. Not part of the public surface."""

    model: Model
    reference: FluxState
    product: str
    biomass: str
    interventions: list[Intervention] = field(default_factory=list)
    previous_bounds: list[tuple[tuple[str, float, float], ...]] = field(
        default_factory=list
    )
    history: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    scans: ScanCache = field(default_factory=ScanCache)
    #: reaction id -> the moves already tried on it that did not stick. A move is recorded
    #: here the moment it is reverted, and is never offered on that reaction again. Without
    #: this the agent re-proposed an identical rejected move every tick, because the state it
    #: sees afterwards is the same state it saw before.
    failed_moves: dict[str, set[str]] = field(default_factory=dict)
    #: Scan results as they stood before each applied intervention, so an undo gives them back.
    scan_stack: list[ScanCache] = field(default_factory=list)
    #: Lines shown to the agent under ``notes`` in the state. Not history and not evidence
    #: about one reaction: facts about the situation it is now in.
    state_notes: list[str] = field(default_factory=list)
    #: One line per completed round, shown to the agent so a later round can either try
    #: something different or go back to what worked. Without it every round starts blind to
    #: what the ones before it achieved.
    round_log: list[str] = field(default_factory=list)
    #: The distinct designs completed rounds ended on, in the order found. Shown to the agent
    #: so a later round can deliberately go elsewhere: six rounds returning the same design
    #: have produced one result, not six.
    designs_found: list[str] = field(default_factory=list)
    #: reaction id -> why later rounds may not use it. One member of each design already found,
    #: so the next round has to reach somewhere else. Stated to the agent rather than applied
    #: invisibly: a reaction that vanishes from the board with no explanation is a reaction the
    #: agent will waste steps looking for.
    round_bans: dict[str, str] = field(default_factory=dict)
    #: Reactions the run definition put off limits. Unlike ``round_bans`` these never lift:
    #: the person running this said not to touch them, and that is not the agent's call.
    forbidden: frozenset[str] = frozenset()
    #: The best design the run has seen, and what it scored, so ``restore_best_design`` has
    #: somewhere to go back to.
    best_snapshot: tuple[Intervention, ...] = ()
    best_score: float = float("-inf")
    #: True when the design has changed since the intervention gains were last measured.
    screen_stale: bool = True
    #: reaction id -> the most recent evidence the board carried for it, kept so the run can
    #: report on every target it measured and not only on the ones it played. A reaction CMM
    #: checked and rejected is a result; dropping it would make the report a record of the
    #: agent's attention rather than of the evidence.
    candidates_seen: dict[str, CandidateEvidence] = field(default_factory=dict)
    #: reaction id -> the change in product flux measured when that intervention was applied.
    #: Shown next to each active intervention so dead weight is visible: an intervention that
    #: bought nothing still occupies one of the design's places, and withdrawing it is a real
    #: move the agent can only choose if it can see the cost.
    contribution: dict[str, float] = field(default_factory=dict)

    def apply(self, intervention: Intervention) -> None:
        # Every reaction the gene edit constrains, not only the one the agent named. A
        # shared gene takes its other reactions with it, and applying the bound to the chosen
        # reaction alone would model a strain that cannot be built.
        self.previous_bounds.append(
            tuple(
                (
                    rid,
                    float(self.model.reactions.get_by_id(rid).lower_bound),
                    float(self.model.reactions.get_by_id(rid).upper_bound),
                )
                for rid, _, _ in intervention.bounds
            )
        )
        self.scan_stack.append(self.scans.snapshot())
        for rid, lower, upper in intervention.bounds:
            self.model.reactions.get_by_id(rid).bounds = (lower, upper)
        self.interventions.append(intervention)
        self.scans.invalidate()
        self.screen_stale = True

    def undo(self) -> Intervention | None:
        if not self.interventions:
            return None
        intervention = self.interventions.pop()
        self.contribution.pop(intervention.reaction_id, None)
        for reaction_id, lower, upper in self.previous_bounds.pop():
            self.model.reactions.get_by_id(reaction_id).bounds = (lower, upper)
        # The model is back where it was, so the scans taken before the change are valid
        # again. Clearing them here made the agent re-run the same scan every tick.
        if self.scan_stack:
            self.scans.restore(self.scan_stack.pop())
        self.screen_stale = True
        return intervention

    def record_failure(self, reaction_id: str, action_name: str) -> None:
        self.failed_moves.setdefault(reaction_id, set()).add(action_name)

    def clear_failures(self) -> None:
        """Forget every rejected move, because the design they were rejected against is gone.

        Whether a move is feasible depends on the bounds it is applied to, so a rejection is
        a fact about one design and not about the move. Keeping it forever was measurably
        wrong: ``force_on_high`` on ``FUM`` breaches the growth floor on the wild type and
        succeeds once ``FRD7`` is carrying flux, and a run that had tried it early could
        never reach the better design afterwards.

        A rejection *is* still valid while the design does not change, which is the case that
        matters — a reverted move leaves the model exactly as it was, so its rejection
        survives until something else sticks.
        """

        self.failed_moves.clear()

    def exhausted(self, reaction_id: str, reference_flux: float) -> bool:
        """True when every move defined for this reaction has already been tried."""

        tried = self.failed_moves.get(reaction_id)
        if not tried:
            return False
        return all(
            action.name in tried for action in applicable_actions(reference_flux)
        )

    def restore(self, design: Sequence[Intervention]) -> None:
        """Replace the current design with ``design``, wholesale."""

        while self.interventions:
            self.undo()
        for intervention in design:
            self.apply(intervention)
        self.contribution.clear()
        self.clear_failures()

    def active_labels(self) -> tuple[str, ...]:
        labels = []
        for intervention in self.interventions:
            delta = self.contribution.get(intervention.reaction_id)
            if delta is None:
                labels.append(intervention.describe())
            elif abs(delta) <= 1e-9:
                labels.append(
                    f"{intervention.describe()} \u2014 bought no product at all when applied"
                )
            else:
                labels.append(
                    f"{intervention.describe()} \u2014 changed the product by {delta:+.4g} "
                    "when applied"
                )
        return tuple(labels)


def run_jev_design(
    config: JevConfig,
    *,
    client: JevClient | None = None,
    on_tick: Callable[[TickRecord, Mapping[str, float]], None] | None = None,
    should_stop: Callable[[], bool] | None = None,
) -> JevResult:
    """Run the JEV game loop and, when ``output_dir`` is set, write a run bundle.

    ``client`` is injectable so a test can play the whole game against a scripted agent with
    no network. ``on_tick`` is called after every frame with the tick record and the flux
    distribution it produced, which is how the desktop app redraws the flux map mid-round.

    ``should_stop`` is polled once per step, from the thread the loop runs on. A run that is
    stopped this way still returns everything it played and still reports its best design —
    stopping is an answer about how long to look, not a reason to throw the answer away. What
    it does skip is the baseline comparison, which costs an OptKnock solve and is not what
    someone pressing stop is waiting for; the run says so in its notes.
    """

    from cobra.io import read_sbml_model

    config.validate()
    # Before anything is loaded or solved: an agent run without a credential fails, and it
    # should fail in the first second rather than after a genome-scale model and a yield
    # calculation have been paid for.
    agent = client or JevClient(
        model=config.jev_model,
        research_model=config.research_model,
        timeout_s=config.request_timeout_s,
        seed=config.seed,
    )

    model = read_sbml_model(str(config.model_path))
    if config.solver:
        model.solver = config.solver

    notes: list[str] = []

    # -- condition ----------------------------------------------------------
    if config.medium is not None:
        from cmm.core.media import apply_medium, preset_medium

        medium = (
            preset_medium(config.medium)
            if isinstance(config.medium, str)
            else config.medium
        )
        application = apply_medium(model, medium)
        missing = getattr(application, "missing", ()) or ()
        if missing:
            notes.append(
                f"medium {medium.name!r} names {len(missing)} exchange(s) this model does "
                f"not contain; they were skipped: {', '.join(sorted(missing)[:8])}"
            )
    if config.condition is not None:
        config.condition.apply_to(model)

    product = _resolve_exchange(model, config.product)
    biomass = config.biomass or _objective_reaction_id(model)
    if biomass is None:
        raise JevWorkflowError(
            "the model has no objective reaction and config.biomass was not set; "
            "growth cannot be scored"
        )

    # -- solver gate, before anything expensive ------------------------------
    use_linear_moma = False
    if config.run_moma and not supports("QP", model.solver.interface):
        use_linear_moma = True
        status = solver_status(model)
        notes.append(
            f"MOMA was run as the LP L1 variant, not the published L2 variant: the active "
            f"solver {status.name!r} provides {', '.join(status.capabilities)} and L2 MOMA "
            "needs QP. Install gurobi, cplex or osqp for the L2 distance."
        )

    # -- wild type ----------------------------------------------------------
    wild_type = _solve(model)
    if wild_type.status != "optimal":
        raise JevWorkflowError(
            f"the wild-type pFBA solve is {wild_type.status}: the model does not grow under "
            "this condition, so nothing downstream would mean anything"
        )
    wild_product = float(wild_type.fluxes.get(product, 0.0))
    wild_growth = float(wild_type.fluxes.get(biomass, 0.0))
    if wild_growth < config.growth_floor:
        notes.append(
            f"wild-type growth ({wild_growth:.4g} per hour) is already below the configured "
            f"floor ({config.growth_floor}); every move will be reverted on the floor rule"
        )
    reference = FluxState(
        wild_type.fluxes,
        name="wild_type",
        provenance="pfba",
        metadata={"source_method": "pfba"},
    )

    substrate = config.substrate or detect_substrate(model, wild_type.fluxes)
    if config.substrate is None and substrate is not None:
        notes.append(
            f"the yield is quoted per {substrate}, the carbon source this condition actually "
            "feeds the model"
        )
    theoretical = _theoretical_max_yield(model, product, substrate, notes)

    # Resolved once: the cofactor pools and the currency metabolites of *this* model, found
    # by formula so a model that does not use BiGG ids is handled rather than silently
    # mis-read. What could not be found is reported, never left implicit.
    pools = resolve_pools(model)
    if pools.missing:
        notes.append(
            "these cofactor pools could not be identified in this model, so the run makes "
            f"no claim about them: {', '.join(pools.missing)}"
            + (
                ""
                if pools.by_formula
                else " (the model carries no formulas, so identification fell back to BiGG "
                "id stems)"
            )
        )

    forbidden, unmatched = _resolve_off_limits(model, config.off_limits)
    if unmatched:
        raise JevWorkflowError(
            "these off-limits names match nothing in this model, so the constraint would "
            f"be silently ignored: {', '.join(sorted(unmatched))}. They are matched against "
            "reaction ids, gene ids, gene names and subsystem names."
        )
    if forbidden:
        notes.append(
            f"{len(forbidden)} reaction(s) were held off limits by the run definition and "
            f"never offered to the agent: {', '.join(sorted(forbidden))}"
        )

    board = _Board(model=model, reference=reference, product=product, biomass=biomass)
    board.forbidden = forbidden

    if config.seed_with_strain_design:
        # Before the first move: hand the agent what the deterministic designer already knows.
        # The reactions it names go on the board with the guaranteed product they buy, which
        # is the only way an escape route carrying no flux today can ever be considered.
        seeded = _run_scan(board, "strain_design_scan", (), config)
        board.scans.completed.add("strain_design_scan")
        notes.append(f"strain design seeded the board before the first move: {seeded}")

    # One web search, before the game, and none inside it. It used to be a lookup per
    # candidate reaction run between the two calls of a step: correct, and unusable — a web
    # search takes tens of seconds, the loop waits on it with nothing to do, and eight of them
    # turned a run that plays in a minute into one that takes ten. Read once, play with it in
    # hand.
    literature_brief = ""
    literature_sources: tuple[str, ...] = ()
    if config.enable_web_research:
        try:
            literature_brief, urls = agent.web_research(
                literature_briefing(config.product, config.organism)
            )
            literature_sources = tuple(urls)
        except JevTransportError as error:
            notes.append(
                f"the literature briefing could not be fetched, so the run played without "
                f"it: {error}"
            )
        else:
            notes.append(
                f"a literature briefing on {config.product} in {config.organism} was read "
                f"once before the first move, from {len(literature_sources)} source(s), and "
                "shown to the agent on every step. It is data the agent weighs, not an "
                "instruction: it cannot name a move outside this package's vocabulary."
            )

    question_set = get_question_set(config.question_set)
    ticks: list[TickRecord] = []
    round_records: list[RoundRecord] = []
    flux_frames: list[Mapping[str, float]] = [dict(wild_type.fluxes)]
    transcript: list[Mapping[str, object]] = []

    idle_rounds = 0
    stopped_by_user = False
    best_product = wild_product
    best_growth = wild_growth
    best_interventions: tuple[Intervention, ...] = ()

    current_product, current_growth = wild_product, wild_growth
    stop_run = False

    for round_index in range(1, config.rounds + 1):
        # A fresh attempt: the bounds go back to the wild type, and what the agent keeps from
        # the last round is the record of what it reached, not the design that reached it.
        while board.interventions:
            board.undo()
        board.contribution.clear()
        board.clear_failures()
        board.scans.invalidate()
        board.screen_stale = True
        board.state_notes.clear()
        current_product, current_growth = wild_product, wild_growth

        round_start_product = current_product
        ended_early = False
        ticks_this_round = 0
        last_signature: tuple[object, ...] | None = None
        repeats = 0

        for tick_index in range(1, config.steps_per_round + 1):
            if should_stop is not None and should_stop():
                notes.append(
                    f"{STOP_REQUESTED} during round {round_index}, step {tick_index}; "
                    "everything played up to that point is kept"
                )
                stopped_by_user = True
                stop_run = True
                break
            if _budget_exhausted(agent, config):
                notes.append(
                    f"stopped in round {round_index}: the agent budget was reached "
                    f"({agent.usage.calls} decisions, ${agent.usage.cost_usd:.4f})"
                )
                stop_run = True
                break

            ticks_this_round += 1
            solution = _solve(board.model)
            if solution.status != "optimal":
                # Every move that breaks feasibility is reverted as it happens, so reaching
                # here means the board itself is unusable. Stop and say so rather than
                # building a screen out of an empty flux vector.
                notes.append(
                    f"round {round_index} tick {tick_index}: the model is "
                    f"{solution.status} before any move was made this tick; the run stopped"
                )
                stop_run = True
                break
            current_fluxes = dict(solution.fluxes)
            current_product = float(current_fluxes.get(product, 0.0))
            current_growth = float(current_fluxes.get(biomass, 0.0))

            exhausted = [
                reaction_id
                for reaction_id in board.failed_moves
                if board.exhausted(
                    reaction_id, float(reference.fluxes.get(reaction_id, 0.0))
                )
            ]
            candidates = board.scans.apply(
                build_candidates(
                    board.model,
                    product_reaction_id=product,
                    reference_fluxes=reference.fluxes,
                    current_fluxes=current_fluxes,
                    excluded=[
                        *(i.reaction_id for i in board.interventions),
                        *exhausted,
                        *board.round_bans,
                        *board.forbidden,
                    ],
                    pinned=board.scans.design_notes,
                    limit=config.candidate_limit,
                    pools=pools,
                )
            )
            if not candidates:
                notes.append(
                    f"round {round_index} tick {tick_index}: no candidate reactions are "
                    "left to act on"
                )
                ended_early = True
                break

            if config.screen_interventions and board.screen_stale:
                _run_scan(board, "intervention_screen", candidates, config)
                board.screen_stale = False
                candidates = board.scans.apply(candidates)

            used_knockouts = sum(1 for i in board.interventions if i.mode == "knockout")
            used_knockdowns = len(board.interventions) - used_knockouts
            for evidence in candidates:
                board.candidates_seen[evidence.reaction_id] = evidence
            state = GameState(
                product_reaction_id=product,
                product_flux=current_product,
                growth=current_growth,
                wild_type_product_flux=wild_product,
                wild_type_growth=wild_growth,
                theoretical_max_yield=theoretical,
                molar_yield=None,
                growth_floor=config.growth_floor,
                balance=cofactor_balance(board.model, current_fluxes, pools),
                unresolved_pools=pools.missing,
                cofactor_limits=(
                    cofactor_limitation(
                        board.model,
                        product=product,
                        biomass=biomass,
                        growth_floor=config.growth_floor,
                        pools=pools,
                    )
                    if config.measure_cofactor_limits
                    else {}
                ),
                guaranteed=(
                    guaranteed_product(board.model, product=product, biomass=biomass)
                    if config.measure_guaranteed_product
                    else None
                ),
                previous_rounds=tuple(board.round_log),
                designs_found=tuple(board.designs_found),
                brief=config.brief,
                literature_brief=literature_brief,
                literature_sources=literature_sources,
                candidates=candidates,
                active_interventions=board.active_labels(),
                ruled_out=tuple(
                    f"{reaction_id}: {move}"
                    for reaction_id, moves in sorted(board.failed_moves.items())
                    for move in sorted(moves)
                ),
                history=tuple(board.history),
                notes=tuple(
                    [
                        *(
                            [
                                "The run definition puts these off limits and they are not "
                                "on the board at all: "
                                + ", ".join(sorted(board.forbidden))
                                + ". That is the person running this telling you what they "
                                "will not build, and it is not open to argument."
                            ]
                            if board.forbidden
                            else []
                        ),
                        *(
                            [
                                "This round may not use "
                                + ", ".join(sorted(board.round_bans))
                                + ", so that it reaches a design an earlier round did not. "
                                + "; ".join(
                                    f"{rid}: {why}"
                                    for rid, why in sorted(board.round_bans.items())
                                )
                            ]
                            if board.round_bans
                            else []
                        ),
                        *board.state_notes,
                        *(
                            [board.scans.envelope_note]
                            if board.scans.envelope_note
                            else []
                        ),
                    ]
                ),
                round_index=round_index,
                tick_index=tick_index,
                ticks_left=config.steps_per_round - tick_index,
                knockouts_used=used_knockouts,
                max_knockouts=config.max_knockouts,
                knockdowns_used=used_knockdowns,
                max_knockdowns=config.max_knockdowns,
            )

            try:
                tick = _play_tick(
                    board=board,
                    agent=agent,
                    config=config,
                    question_set=question_set,
                    state=state,
                    candidates=candidates,
                    round_index=round_index,
                    tick_index=tick_index,
                    use_linear_moma=use_linear_moma,
                    transcript=transcript,
                    allowed_modes=config.room_for(used_knockouts, used_knockdowns),
                    rounds_found=len(board.designs_found),
                )
            except JevTransportError as error:
                # The network failed after the retries. Everything played up to here is real
                # and was paid for in solver time, so the run ends the way a stop request
                # ends it rather than throwing the work away. Observed: a read timeout on the
                # second call of a six-round run discarded the whole thing.
                notes.append(
                    f"the run stopped in round {round_index}, step {tick_index}: {error}. "
                    "Everything played before that is kept and scored; the design is "
                    "whatever the run had reached, not whatever it would have reached."
                )
                stop_run = True
                stopped_by_user = True
                break
            ticks.append(tick)
            board.history.append(tick.headline())

            # Circuit breaker. Every dead end found so far has been closed at its source, but
            # an agent that cannot make progress should end its round rather than spend the
            # remaining ticks discovering the same refusal. A repeat is the same reaction and
            # the same move landing the same way.
            signature = (tick.target, tick.action, tick.outcome)
            repeats = repeats + 1 if signature == last_signature else 0
            last_signature = signature
            if repeats >= _MAX_IDENTICAL_MOVES:
                board.history.append(
                    f"R{round_index}T{tick_index} the same move repeated "
                    f"{repeats + 1} times; the round was ended"
                )
                ended_early = True
                break

            post = _solve(board.model)
            post_fluxes = dict(post.fluxes)
            flux_frames.append(post_fluxes)
            current_product = float(post_fluxes.get(product, 0.0))
            current_growth = float(post_fluxes.get(biomass, 0.0))

            if on_tick is not None:
                on_tick(tick, post_fluxes)

            if (
                current_growth >= config.growth_floor
                and current_product > best_product + 1e-9
            ):
                best_product = current_product
                best_growth = current_growth
                best_interventions = tuple(board.interventions)
                board.best_snapshot = best_interventions
                board.best_score = current_product

            if tick.outcome == "end_round":
                ended_early = True
                break

        design = (
            "; ".join(i.describe() for i in board.interventions) or "no interventions"
        )
        design_key = _design_signature(board.interventions)
        repeated = next(
            (
                record.round_index
                for record in round_records
                if record.signature == design_key and design_key
            ),
            None,
        )
        round_ticks = ticks[-ticks_this_round:] if ticks_this_round else []
        withheld = tuple(sorted(board.round_bans))
        shortfall = _round_shortfall(
            board,
            config,
            pools=pools,
            ticks=round_ticks,
            product=current_product,
            growth=current_growth,
        )
        stopped = _why_the_round_stopped(
            round_ticks, ended_early=ended_early, steps=config.steps_per_round
        )

        board.round_log.append(
            (
                ""
                if not withheld
                else f"(round {round_index} could not use {', '.join(withheld)}, so it "
                "answers: what is the best design without them?) "
            )
            + f"round {round_index} reached {current_product:.4g} product at "
            f"{current_growth:.4g} growth with {design}"
            + (
                f" (best so far: {board.best_score:.4g})"
                if board.best_score > current_product + 1e-9
                else " (the best round so far)"
            )
            + f". It stopped because {stopped}."
            + (
                " What it left undone: " + "; ".join(shortfall) + "."
                if shortfall
                else " Nothing measurable was left undone."
            )
            + (
                f" NOTE: this is the same design round {repeated} already found, so the two "
                "rounds produced one result between them, not two."
                if repeated is not None
                else ""
            )
        )

        acted = any(
            tick.outcome in ("applied", "undone_by_agent") for tick in round_ticks
        )
        idle_rounds = 0 if acted else idle_rounds + 1

        if config.require_distinct_rounds and design_key:
            _ban_one_member(board, config, round_index)
        if design_key and repeated is None:
            board.designs_found.append(
                f"round {round_index}: {design} \u2014 {current_product:.4g} product at "
                f"{current_growth:.4g} growth"
            )
        round_records.append(
            RoundRecord(
                round_index=round_index,
                n_ticks=ticks_this_round,
                product_flux=current_product,
                growth=current_growth,
                interventions=board.active_labels(),
                improved=current_product > round_start_product + 1e-9,
                ended_early=ended_early,
                stopped_because=stopped,
                shortfall=shortfall,
                signature=design_key,
                repeated=repeated,
                withheld=withheld,
            )
        )
        if idle_rounds >= _MAX_IDLE_ROUNDS:
            notes.append(
                f"the run stopped after round {round_index}: {idle_rounds} consecutive "
                "rounds ended without the agent changing anything, and the state it is "
                "shown does not change between them"
            )
            break
        if stop_run:
            break

    provenance = {
        **run_provenance(model, method="jev_target_design", seed=config.seed),
        **config.to_provenance(),
        "jev_model_served": ", ".join(agent.served_models) or None,
        "question_set_version": question_set.version,
        "moma_variant": (
            None
            if not config.run_moma
            else ("moma_l1" if use_linear_moma else "moma_l2")
        ),
        "wild_type_reference": "pfba",
        "cofactor_pools_resolved_by": "formula"
        if pools.by_formula
        else "bigg_id_stems",
        "cofactor_pools_not_found": list(pools.missing),
    }

    # The design the run ended on, captured before the comparison below strips the board
    # back to the unmodified model. Reading it afterwards returned an empty design.
    final_interventions = tuple(board.interventions)

    baselines: tuple["BaselineRow", ...] = ()
    if stopped_by_user and config.run_baseline_comparison:
        notes.append(
            "the baseline comparison was skipped because the run was stopped; the agent's "
            "own numbers are complete, but there is nothing here to compare them against"
        )
    if config.run_baseline_comparison and not stopped_by_user:
        from cmm.jev.benchmark import compare_with_baselines

        # Every method's design is applied by the comparison itself, within its own reverting
        # context, so the board's edits must not still be standing while it runs.
        for _ in list(board.interventions):
            board.undo()
        try:
            baselines = compare_with_baselines(
                model,
                product=product,
                biomass=biomass,
                growth_floor=config.growth_floor,
                jev_interventions=best_interventions,
                max_knockouts=config.design_max_knockouts,
                max_solutions=config.design_max_solutions,
                seed=config.seed,
            )
        except Exception as error:  # a comparison that fails must not lose the run
            notes.append(f"the baseline comparison could not run: {error}")

    result = JevResult(
        config=config,
        provenance=provenance,
        wild_type_product_flux=wild_product,
        wild_type_growth=wild_growth,
        theoretical_max_yield=theoretical,
        ticks=tuple(ticks),
        rounds=tuple(round_records),
        best_product_flux=best_product,
        best_growth=best_growth,
        best_interventions=best_interventions,
        final_interventions=final_interventions,
        flux_frames=tuple(flux_frames),
        transcript=tuple(transcript),
        baselines=baselines,
        literature_brief=literature_brief,
        literature_sources=literature_sources,
        candidates_seen=tuple(board.scans.apply(tuple(board.candidates_seen.values()))),
        notes=tuple([*notes, *board.notes]),
        usage=agent.usage.to_dict(),
    )

    if config.output_dir is None:
        return result

    from cmm.jev.artifacts import export_run

    return export_run(result, model=model, reference=reference)


# ---------------------------------------------------------------------------
# one tick
# ---------------------------------------------------------------------------


def _play_tick(
    *,
    board: _Board,
    agent: JevClient,
    config: JevConfig,
    question_set,
    state: GameState,
    candidates: Sequence[CandidateEvidence],
    round_index: int,
    tick_index: int,
    use_linear_moma: bool,
    transcript: list[Mapping[str, object]],
    allowed_modes: tuple[str, ...],
    rounds_found: int = 0,
) -> TickRecord:
    """Ask, execute, measure. Returns the frame; the caller re-solves and redraws."""

    payload = state.to_payload()
    allow_undo = bool(board.interventions)
    allow_look = config.allow_look_actions

    # A proven design is all deletions, so it is the deletion budget it has to fit inside.
    knockout_room = config.max_knockouts - sum(
        1 for i in board.interventions if i.mode == "knockout"
    )
    best_design = next(
        (
            design
            for design in board.scans.designs
            if len(design[1]) <= knockout_room
            and not set(design[1]) & {i.reaction_id for i in board.interventions}
            # A proven design that reaches through a withheld reaction would walk straight
            # back into the design the cut exists to move away from.
            and not set(design[1]) & set(board.round_bans)
            and not set(design[1]) & board.forbidden
        ),
        None,
    )
    target_questions = question_set.target_question(
        candidates,
        product=board.product,
        growth_floor=config.growth_floor,
        allow_undo=allow_undo,
        allow_look=allow_look,
        design_full=not allowed_modes,
        best_design=(
            f"{len(board.best_snapshot)} interventions reaching {board.best_score:.4g}"
            if board.best_snapshot
            and tuple(board.best_snapshot) != tuple(board.interventions)
            and not {i.reaction_id for i in board.best_snapshot} & set(board.round_bans)
            else ""
        ),
        proven_design=(
            f"{best_design[0]} deletes {', '.join(best_design[1])} for "
            f"{best_design[2]:.3g} guaranteed product"
            if best_design
            else ""
        ),
        rounds_found=rounds_found,
    )
    stage1 = agent.decide(payload, target_questions)
    _record(transcript, round_index, tick_index, "target", payload, stage1)
    target_answer = stage1[TARGET_KEY]
    target = target_answer.choice
    ranking = target_answer.ranked()
    cost = stage1.cost_usd
    latency = stage1.latency_s

    action_ranking: tuple[tuple[str, float], ...] = ()

    def frame(
        *,
        action: str | None,
        outcome: TickOutcome,
        reason: str,
        intervention: Intervention | None = None,
        action_confidence: float | None = None,
        benefit: float | None = None,
        risk: float | None = None,
        moma: tuple[float | None, float | None, float | None] = (None, None, None),
        status: str = "optimal",
        product_flux: float = float("nan"),
        growth: float = float("nan"),
    ) -> TickRecord:
        return TickRecord(
            round_index=round_index,
            tick_index=tick_index,
            target=target,
            target_confidence=target_answer.confidence,
            target_ranking=ranking,
            action_ranking=action_ranking,
            action=action,
            action_confidence=action_confidence,
            benefit_score=benefit,
            predicted_growth_risk=risk,
            outcome=outcome,
            reason=reason,
            intervention=intervention,
            product_flux=product_flux,
            growth=growth,
            status=status,
            moma_product_flux=moma[0],
            moma_growth=moma[1],
            moma_distance=moma[2],
            n_active_interventions=len(board.interventions),
            decision_cost_usd=cost,
            decision_latency_s=latency,
        )

    # -- the two moves that need no second call -----------------------------
    if target == END_ACTION.name:
        return frame(
            action=END_ACTION.name,
            outcome="end_round",
            reason="the agent judged no remaining move worthwhile",
            product_flux=state.product_flux,
            growth=state.growth,
        )

    if target == ADOPT_ACTION.name:
        return _adopt_design(
            board=board,
            config=config,
            design=best_design,
            frame=frame,
            use_linear_moma=use_linear_moma,
            state=state,
        )

    if target == RESTORE_ACTION.name:
        board.restore(board.best_snapshot)
        solution = _solve(board.model)
        return frame(
            action=RESTORE_ACTION.name,
            outcome="applied"
            if solution.status == "optimal"
            else "reverted_infeasible",
            reason=(
                f"went back to the best design this run has found ({board.best_score:.4g} "
                f"product, {len(board.best_snapshot)} interventions)"
            ),
            status=solution.status,
            product_flux=float(solution.fluxes.get(board.product, 0.0)),
            growth=float(solution.fluxes.get(board.biomass, 0.0)),
        )

    if target == UNDO_ACTION.name:
        removed = board.undo()
        board.clear_failures()
        return frame(
            action=UNDO_ACTION.name,
            outcome="undone_by_agent",
            reason=(removed.describe() if removed else "there was nothing to undo"),
            product_flux=state.product_flux,
            growth=state.growth,
        )

    candidate = next((c for c in candidates if c.reaction_id == target), None)
    if (
        candidate is None
    ):  # pragma: no cover - JEV can only answer inside the vocabulary
        return frame(
            action=None,
            outcome="not_applicable",
            reason=f"{target!r} is not on the board",
            product_flux=state.product_flux,
            growth=state.growth,
        )

    # -- stage two ----------------------------------------------------------
    # Moves already ruled out on this reaction, plus every model-wide scan already run
    # against the current bounds. Both would spend a tick to arrive back where we are.
    blocked = set(board.failed_moves.get(candidate.reaction_id, ()))
    blocked |= board.scans.completed
    offered = available_actions(
        candidate,
        allow_look=allow_look,
        exclude=blocked,
        allowed_modes=allowed_modes,
    )
    benefit: float | None = None
    risk: float | None = None
    action_confidence: float | None = None

    if not offered:
        # Nothing is left to try here. Mark the whole reaction spent so the next tick's board
        # is built without it, instead of offering a dead end again.
        for spent in applicable_actions(candidate.reference_flux):
            board.record_failure(candidate.reaction_id, spent.name)
        return frame(
            action=None,
            outcome="not_applicable",
            reason=(
                f"every move on {candidate.reaction_id!r} has already been tried, is "
                "undefined for it, or is out of budget"
            ),
            product_flux=state.product_flux,
            growth=state.growth,
        )

    if len(offered) == 1:
        # One option is not a decision. Asking would spend a call and a step to be told the
        # only thing that could be said, which matters now that the budgets are separate:
        # once the knockdown budget is gone, a zero-flux reaction has exactly one legal move.
        action = offered[0]
        action_ranking = ((action.name, 1.0),)
    else:
        action_questions = question_set.action_question(
            candidate,
            product=board.product,
            growth_floor=config.growth_floor,
            allow_look=allow_look,
            exclude=blocked,
            allowed_modes=allowed_modes,
        )
        stage2 = agent.decide(
            {**payload, "selected_reaction": candidate.reaction_id}, action_questions
        )
        _record(
            transcript, round_index, tick_index, "action", candidate.reaction_id, stage2
        )
        cost += stage2.cost_usd
        latency += stage2.latency_s

        action_answer = stage2[ACTION_KEY]
        action_ranking = action_answer.ranked()
        action_confidence = action_answer.confidence
        benefit = _safe_score(stage2, BENEFIT_KEY)
        risk = _safe_noul(stage2, RISK_KEY)
        resolved = ACTION_CATALOGUE.get(action_answer.choice)
        if resolved is None:  # pragma: no cover - vocabulary is closed
            return frame(
                action=action_answer.choice,
                outcome="not_applicable",
                reason=f"{action_answer.choice!r} is not a known move",
                benefit=benefit,
                risk=risk,
                action_confidence=action_confidence,
                product_flux=state.product_flux,
                growth=state.growth,
            )
        action = resolved

    # -- LOOK: run a CMM analysis, change nothing ---------------------------
    if action.kind == "look":
        reason = _run_scan(board, action.name, candidates, config)
        board.scans.completed.add(action.name)
        return frame(
            action=action.name,
            outcome="scan",
            reason=reason,
            benefit=benefit,
            risk=risk,
            action_confidence=action_confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    # -- ACT ----------------------------------------------------------------
    # No budget check here: a move whose budget is spent was never offered, which is the
    # point of gating the question rather than the answer. The alternative — offering it and
    # refusing it — was measured, and it cost ten consecutive steps on one run.
    try:
        intervention = build_intervention(
            board.model, candidate.reaction_id, action, board.reference.fluxes
        )
    except ActionNotApplicable as error:
        board.record_failure(candidate.reaction_id, action.name)
        return frame(
            action=action.name,
            outcome="not_applicable",
            reason=str(error),
            benefit=benefit,
            risk=risk,
            action_confidence=action_confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    board.apply(intervention)
    solution = _solve(board.model)
    # From here a revert leaves the design untouched, so its rejection stays on the record;
    # a move that sticks changes the design and every earlier rejection with it.

    # CMM enforces viability. A move that makes the model infeasible, or drops growth below
    # the floor, is undone here whatever the agent predicted — an infeasible design is not a
    # design. A legal but unhelpful move is kept for the agent to withdraw itself.
    if solution.status != "optimal":
        board.undo()
        board.record_failure(candidate.reaction_id, action.name)
        return frame(
            action=action.name,
            outcome="reverted_infeasible",
            reason=f"the model became {solution.status} with this intervention applied",
            intervention=intervention,
            benefit=benefit,
            risk=risk,
            action_confidence=action_confidence,
            status=solution.status,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    new_product = float(solution.fluxes.get(board.product, 0.0))
    new_growth = float(solution.fluxes.get(board.biomass, 0.0))
    if new_growth < config.growth_floor:
        board.undo()
        board.record_failure(candidate.reaction_id, action.name)
        # The move was too strong, not wrong. Say that the same move has a gentler version
        # still on offer for this reaction, because the agent does not infer it: watching a
        # run, a refused ``force_on_high`` sent it to a different reaction and left behind
        # what ``force_on_low`` on the same one would have collected.
        gentler = GENTLER_ALTERNATIVE.get(action.name)
        still_open = (
            gentler is not None
            and gentler not in board.failed_moves.get(candidate.reaction_id, set())
            and gentler
            in {
                available.name
                for available in applicable_actions(candidate.reference_flux)
            }
        )
        hint = (
            f"; the move was too strong, not wrong \u2014 {gentler} on the same reaction is "
            "still available"
            if still_open
            else ""
        )
        return frame(
            action=action.name,
            outcome="reverted_growth_floor",
            reason=(
                f"growth would fall to {new_growth:.4g} per hour, below the floor of "
                f"{config.growth_floor}{hint}"
            ),
            intervention=intervention,
            benefit=benefit,
            risk=risk,
            action_confidence=action_confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    moma = _moma_snapshot(board, config, use_linear_moma)
    delta = new_product - state.product_flux
    board.contribution[intervention.reaction_id] = delta
    if delta > 1e-9:
        verdict = f"product rose by {delta:+.4g}"
    elif delta < -1e-9:
        verdict = f"product FELL by {delta:+.4g}"
    else:
        verdict = (
            "product unchanged: this edit costs a place in the design and has bought "
            "nothing so far"
        )
    # The move stuck, so the design that every earlier rejection was measured against is
    # gone, and with it the grounds for the rejection.
    board.clear_failures()
    return frame(
        action=action.name,
        outcome="applied",
        reason=verdict,
        intervention=intervention,
        benefit=benefit,
        risk=risk,
        action_confidence=action_confidence,
        moma=moma,
        product_flux=new_product,
        growth=new_growth,
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _solve(model: Model) -> FluxSolution:
    """pFBA that *reports* infeasibility instead of raising it.

    An intervention that kills the strain is a result, not a crash — CMM's own rule is that
    an infeasible point is data and stays in the table. COBRApy's ``pfba`` raises
    ``Infeasible`` from inside ``fix_objective_as_constraint`` before any status can be read,
    so feasibility is probed first with a slim solve, the same way
    :mod:`cmm.features.comparison` does before a MOMA solve.
    """

    if not math.isfinite(model.slim_optimize(error_value=float("nan"))):
        return FluxSolution(
            status="infeasible", objective_value=None, fluxes={}, metadata={}
        )
    try:
        return pfba(model)
    except OptimizationError as error:  # pragma: no cover - the probe above catches it
        return FluxSolution(
            status="infeasible",
            objective_value=None,
            fluxes={},
            metadata={"error": str(error)},
        )


def _adopt_design(
    *,
    board: _Board,
    config: JevConfig,
    design: tuple[str, tuple[str, ...], float] | None,
    frame,
    use_linear_moma: bool,
    state: GameState,
) -> TickRecord:
    """Apply a whole proven knockout set as one move.

    A design's deletions are only worth anything together. ``ACALD``, ``LDH_D`` and ``THD2``
    reach 9.9 mmol gDW^-1 h^-1 of succinate as a set, while ``ACALD`` alone reaches almost
    nothing — so an agent that judges each move by the product change it causes will abandon
    the design after the first deletion. Offering the set as one move is what makes the
    deterministic result reachable, and it leaves the interesting question open: whether a
    partial knockdown on top of a proven design beats the design alone, which is a question
    the designer itself cannot ask, its variables being present-or-absent.

    The whole set is reverted together if it turns out infeasible or unviable, because a
    partially applied design is not the thing that was proven.
    """

    if design is None:  # pragma: no cover - the move is not offered without one
        return frame(
            action=ADOPT_ACTION.name,
            outcome="not_applicable",
            reason="no proven design is available; run a strain design scan first",
            product_flux=state.product_flux,
            growth=state.growth,
        )

    method, knockouts, guaranteed = design
    applied: list[Intervention] = []
    for reaction_id in knockouts:
        try:
            intervention = build_intervention(
                board.model,
                reaction_id,
                ACTION_CATALOGUE["knockout"],
                board.reference.fluxes,
            )
        except (ActionNotApplicable, KeyError) as error:
            for _ in applied:
                board.undo()
            return frame(
                action=ADOPT_ACTION.name,
                outcome="not_applicable",
                reason=f"the design could not be applied: {error}",
                product_flux=state.product_flux,
                growth=state.growth,
            )
        board.apply(intervention)
        applied.append(intervention)

    solution = _solve(board.model)
    if solution.status != "optimal":
        for _ in applied:
            board.undo()
        return frame(
            action=ADOPT_ACTION.name,
            outcome="reverted_infeasible",
            reason=f"the {method} design left the model {solution.status}",
            intervention=applied[-1] if applied else None,
            status=solution.status,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    new_product = float(solution.fluxes.get(board.product, 0.0))
    new_growth = float(solution.fluxes.get(board.biomass, 0.0))
    if new_growth < config.growth_floor:
        for _ in applied:
            board.undo()
        return frame(
            action=ADOPT_ACTION.name,
            outcome="reverted_growth_floor",
            reason=(
                f"the {method} design drops growth to {new_growth:.4g} per hour, below the "
                f"floor of {config.growth_floor}"
            ),
            intervention=applied[-1] if applied else None,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    board.clear_failures()
    # The one thing the designer cannot have considered, stated plainly because it is the
    # whole remaining opportunity. Without a line like this, eight runs out of eight adopted
    # the design and ended the round on the spot.
    note = (
        f"The active design came from {method}, whose formulation searches complete "
        "deletions only: every gene is either present or absent. A knockdown to half of "
        "wild type is a constraint it cannot express, so a knockdown on top of this design "
        "is a move it could not have considered, and is the kind of move most likely to "
        "improve on a proven one. The measured knockdown gains on the board are for the "
        "design as it now stands, with these deletions applied."
    )
    if note not in board.state_notes:
        board.state_notes.append(note)
    for intervention in applied:
        board.contribution[intervention.reaction_id] = 0.0
    # The set bought the change, so the whole change is credited to its last member rather
    # than split arbitrarily between deletions that mean nothing on their own.
    if applied:
        board.contribution[applied[-1].reaction_id] = new_product - state.product_flux
    return frame(
        action=ADOPT_ACTION.name,
        outcome="applied",
        reason=(
            f"adopted the {method} design ({', '.join(knockouts)}), "
            f"proven for {guaranteed:.3g} guaranteed product"
        ),
        intervention=applied[-1] if applied else None,
        moma=_moma_snapshot(board, config, use_linear_moma),
        product_flux=new_product,
        growth=new_growth,
    )


def _ban_one_member(board: _Board, config: JevConfig, round_index: int) -> None:
    """Withhold one reaction of the design just found, so the next round reaches elsewhere.

    This is the integer cut OptKnock uses to enumerate alternative designs, in the only form
    a board of single reactions can express: forbid one member, and every design containing it
    becomes unreachable. Which member decides how far the next round is pushed, so it is the
    one that did the most — the largest measured contribution to the product — because
    banning a member that bought nothing would leave the same design one substitution away.

    Bans accumulate across rounds, which is what makes a six-round run an enumeration rather
    than six attempts at the same answer. They stop accumulating when the board would be left
    too thin to play on; a run that cannot find anything new should say so rather than play
    rounds with nothing on the table.
    """

    members = [i.reaction_id for i in board.interventions]
    if not members:
        return
    remaining = sum(
        1
        for reaction in board.model.reactions
        if reaction.genes and reaction.id not in board.round_bans
    )
    if remaining <= config.candidate_limit:
        board.notes.append(
            f"round {round_index} found a design but nothing further was withheld: the "
            "board would have been left too thin to play on, so later rounds may repeat it"
        )
        return

    ranked = sorted(
        members,
        key=lambda rid: (-abs(board.contribution.get(rid, 0.0)), rid),
    )
    banned = ranked[0]
    board.round_bans[banned] = (
        f"it carried the design round {round_index} found. Withholding it is what makes the "
        "next round answer a different question \u2014 what is the best design without it? "
        "\u2014 which is the question a laboratory that cannot edit it needs answered"
    )


def _resolve_off_limits(
    model: Model, names: Sequence[str]
) -> tuple[frozenset[str], set[str]]:
    """Turn what a person wrote into the reactions CMM will refuse to touch.

    A name may be a reaction id, a gene id, a gene name or a subsystem, because those are the
    four ways people describe the thing they will not build: "don't touch ``PFL``", "leave
    *ldhA* alone", "the pentose phosphate pathway is off the table". Matching is
    case-insensitive; nothing else is inferred.

    A name matching nothing comes back in the second set, and the caller fails the run on it.
    Silently ignoring a constraint is the one behaviour that must not happen here: the person
    would get a design built on exactly what they said they could not do, with nothing on
    screen to say the instruction was dropped.
    """

    # Keyed by the folded name, valued by what the person actually typed: an error that
    # reports "notageneorreaction" back at someone who wrote "notAGeneOrReaction" is an error
    # they have to squint at to find in their own input.
    wanted = {
        name.strip().casefold(): name.strip() for name in names if name and name.strip()
    }
    if not wanted:
        return frozenset(), set()

    matched: set[str] = set()
    seen: set[str] = set()
    for reaction in model.reactions:
        keys = {str(reaction.id).casefold()}
        subsystem = str(getattr(reaction, "subsystem", "") or "").strip().casefold()
        if subsystem:
            keys.add(subsystem)
        for gene in reaction.genes:
            keys.add(str(gene.id).casefold())
            name = str(getattr(gene, "name", "") or "").strip().casefold()
            if name:
                keys.add(name)
        hit = keys & set(wanted)
        if hit:
            matched.add(str(reaction.id))
            seen |= hit
    return frozenset(matched), {
        original for folded, original in wanted.items() if folded not in seen
    }


def _design_signature(interventions: Sequence[Intervention]) -> tuple[str, ...]:
    """A design as an order-independent key, so two rounds can be compared.

    Order does not make a different strain. Without sorting, the same three deletions applied
    in two sequences read as two distinct results, which is exactly the illusion of diversity
    a multi-round run should not produce.
    """

    return tuple(sorted(f"{i.reaction_id}:{i.action_name}" for i in interventions))


def _why_the_round_stopped(
    ticks: Sequence[TickRecord], *, ended_early: bool, steps: int
) -> str:
    """One clause naming what ended the round, for the next round to read."""

    if not ticks:
        return "it never got a move in"
    last = ticks[-1]
    if last.outcome == "end_round":
        return "the agent judged no remaining move worth making"
    if len(ticks) >= steps:
        return f"it used all {steps} of its steps"
    if last.outcome == "not_applicable":
        return f"nothing was left it could do ({last.reason})"
    return "it ran out of moves that changed anything"


def _round_shortfall(
    board: _Board,
    config: JevConfig,
    *,
    pools,
    ticks: Sequence[TickRecord],
    product: float,
    growth: float,
) -> tuple[str, ...]:
    """What this round left on the table, measured on the design it ended with.

    Four questions, in the order a next round would ask them. Each is a measurement against
    the *final* design rather than a memory of what the board said earlier, because a gain
    measured three moves ago is a gain against a design that no longer exists.
    """

    lines: list[str] = []

    # 1. Moves this round could still have made and did not. The sharpest of the four: a
    #    reaction the screen says pays, sitting untouched when the round ended.
    if config.screen_interventions:
        board.scans.deletion_gains.clear()
        board.scans.knockdown_gains.clear()
        current = board.scans.apply(
            build_candidates(
                board.model,
                product_reaction_id=board.product,
                reference_fluxes=board.reference.fluxes,
                current_fluxes=dict(_solve(board.model).fluxes),
                excluded=[i.reaction_id for i in board.interventions],
                pinned=board.scans.design_notes,
                limit=config.candidate_limit,
                pools=pools,
            )
        )
        if current:
            _run_scan(board, "intervention_screen", current, config)
            untaken = sorted(
                (
                    (gain, reaction_id, move)
                    for store, move in (
                        (board.scans.deletion_gains, "deleting"),
                        (board.scans.knockdown_gains, "halving"),
                    )
                    for reaction_id, gain in store.items()
                    if gain > 1e-6
                ),
                reverse=True,
            )
            if untaken:
                named = ", ".join(
                    f"{move} {reaction_id} ({gain:+.4g})"
                    for gain, reaction_id, move in untaken[:3]
                )
                lines.append(f"moves still worth making that it did not take: {named}")

    # 2. What the product was short of when it stopped. Carbon routing and cofactor supply
    #    call for different moves, and nothing else in the run tells them apart.
    if config.measure_cofactor_limits:
        limits = cofactor_limitation(
            board.model,
            product=board.product,
            biomass=board.biomass,
            growth_floor=config.growth_floor,
            pools=pools,
        )
        binding = sorted(
            ((value, name) for name, value in limits.items() if value > 1e-6),
            reverse=True,
        )
        if binding:
            value, name = binding[0]
            lines.append(
                f"the product was still {name}-limited: one more unit per hour would have "
                f"bought {value:+.3g} more"
            )

    # 3. Budget it never spent. A round that ends with edits in hand stopped for a reason
    #    other than the budget, and the next round should not assume otherwise.
    used_knockouts = sum(1 for i in board.interventions if i.mode == "knockout")
    spare = [
        f"{config.max_knockouts - used_knockouts} unused deletion(s)"
        if used_knockouts < config.max_knockouts
        else "",
        f"{config.max_knockdowns - (len(board.interventions) - used_knockouts)} "
        "unused knockdown(s)"
        if len(board.interventions) - used_knockouts < config.max_knockdowns
        else "",
    ]
    spare = [item for item in spare if item]
    if spare:
        lines.append("it finished with " + " and ".join(spare))

    # 4. Moves the rules took away from it, which a later round may reach by another route:
    #    a deletion refused on the growth floor here can be affordable on a design that
    #    spends its growth differently.
    refused = [
        f"{tick.target} ({tick.action})"
        for tick in ticks
        if tick.outcome.startswith("reverted")
    ]
    if refused:
        # dict.fromkeys keeps first-seen order while dropping repeats, so a move refused on
        # five consecutive steps is named once.
        lines.append("the rules refused " + ", ".join(list(dict.fromkeys(refused))[:4]))

    if growth - config.growth_floor > 1e-6 and product > 1e-9:
        lines.append(
            f"it ended {growth - config.growth_floor:.4g} per hour above the growth floor, "
            "which is room a bolder design could have spent"
        )
    return tuple(lines)


def _run_scan(
    board: _Board,
    name: str,
    candidates: Sequence[CandidateEvidence],
    config: JevConfig,
) -> str:
    """Execute a LOOK move: a real CMM analysis whose answer enriches the next frame."""

    if name == "essentiality_scan":
        n_essential = 0
        for candidate in candidates:
            reaction = board.model.reactions.get_by_id(candidate.reaction_id)
            saved = reaction.bounds
            reaction.bounds = (0.0, 0.0)
            growth = board.model.slim_optimize(error_value=float("nan"))
            reaction.bounds = saved
            essential = math.isnan(growth) or growth < config.growth_floor
            board.scans.essential[candidate.reaction_id] = bool(essential)
            n_essential += int(essential)
        return f"tested {len(candidates)} reactions; {n_essential} are essential under this floor"

    if name == "fseof_scan":
        from cmm.features.production import fseof

        try:
            scan = fseof(board.model, board.product, board.biomass)
        except Exception as error:  # a scan that cannot run is data, not a crash
            return f"FSEOF could not run: {error}"
        trends = scan.trends
        if "slope" in trends:
            for reaction_id, slope in trends["slope"].items():
                value = float(slope)
                if math.isfinite(value):
                    board.scans.fseof_slopes[str(reaction_id)] = value
        return f"FSEOF scanned {len(trends)} reactions at {len(scan.enforced_levels)} levels"

    if name == "state_distance_check":
        from cmm.features.comparison import knockout_comparison, moma

        if not board.interventions:
            return "there is no design yet; MOMA and ROOM would compare the wild type to itself"
        lines = []
        try:
            adjusted = moma(
                board.model, board.reference, linear=config.run_moma is False
            )
            if adjusted.status == "optimal" and adjusted.distance is not None:
                lines.append(
                    f"MOMA: the cell has to move {adjusted.distance:.4g} from the wild-type "
                    f"flux state, and makes {adjusted.fluxes.get(board.product, 0.0):.4g} "
                    "product before adapting"
                )
        except Exception as error:
            lines.append(f"MOMA could not run: {error}")
        try:
            switched = knockout_comparison(
                board.model,
                board.reference,
                [i.reaction_id for i in board.interventions if i.mode == "knockout"]
                or [board.interventions[0].reaction_id],
                method="room",
            )
            if (
                switched.status == "optimal"
                and switched.n_changed_reactions is not None
            ):
                lines.append(
                    f"ROOM: {switched.n_changed_reactions:.0f} reactions have to change "
                    "their flux for this design to work"
                )
        except Exception as error:
            lines.append(f"ROOM could not run: {error}")
        return "; ".join(lines) or "neither MOMA nor ROOM could be run on this design"

    if name == "intervention_screen":
        # Two solves per candidate, against the design as it stands: what happens if this
        # reaction is deleted, and what happens if it is capped at half its wild-type flux.
        # This is CMM answering the question the agent would otherwise guess at, and the
        # guess is systematically wrong in a way worth naming — which branch competes with
        # the product depends on the whole network at the current bounds, not on the
        # reaction's own stoichiometry, and a deletion that pays on the wild type can be
        # worthless once three other deletions are standing.
        #
        # Both moves are measured because the pair is the decision. A reaction whose deletion
        # is lethal and whose knockdown pays is exactly what the knockdown move exists for,
        # and screening deletions alone would hide it.
        baseline = _solve(board.model)
        if baseline.status != "optimal":
            return "the screen needs a feasible starting point; the model is not"
        before = float(baseline.fluxes.get(board.product, 0.0))
        measured = 0
        best: tuple[str, str, float] | None = None
        for candidate in candidates:
            for action_name, store in (
                ("knockout", board.scans.deletion_gains),
                ("knockdown_50", board.scans.knockdown_gains),
            ):
                try:
                    trial = build_intervention(
                        board.model,
                        candidate.reaction_id,
                        ACTION_CATALOGUE[action_name],
                        board.reference.fluxes,
                    )
                except ActionNotApplicable:
                    # A knockdown of a flux that is already zero. Not a gap in the screen:
                    # the move does not exist, and the agent is never offered it.
                    continue
                # The gene edit's whole consequence, so the measured gain is the gain of the
                # move as it would be built. Screening the chosen reaction alone would report
                # a number the strain never produces.
                saved = [
                    (rid, board.model.reactions.get_by_id(rid).bounds)
                    for rid, _, _ in trial.bounds
                ]
                for rid, low, high in trial.bounds:
                    board.model.reactions.get_by_id(rid).bounds = (low, high)
                solution = _solve(board.model)
                for rid, bounds in saved:
                    board.model.reactions.get_by_id(rid).bounds = bounds
                viable = solution.status == "optimal" and (
                    float(solution.fluxes.get(board.biomass, 0.0))
                    >= config.growth_floor
                )
                if action_name == "knockout":
                    # Essentiality comes free with the deletion solve, so the agent never has
                    # to spend a step on an essentiality scan to learn it.
                    board.scans.essential[candidate.reaction_id] = not viable
                if not viable:
                    # Recorded as no gain rather than omitted: a move that kills the strain is
                    # not a move, and a blank would read as "not yet measured".
                    store[candidate.reaction_id] = 0.0
                    measured += 1
                    continue
                gain = float(solution.fluxes.get(board.product, 0.0)) - before
                store[candidate.reaction_id] = gain
                measured += 1
                if best is None or gain > best[2]:
                    best = (candidate.reaction_id, action_name, gain)
        board.scans.completed.add("essentiality_scan")
        if best is None or best[2] <= 1e-9:
            return f"measured {measured} moves; none of them raises the product"
        return (
            f"measured {measured} moves; the best is {best[1]} on {best[0]} at "
            f"{best[2]:+.4g} product flux"
        )

    if name == "strain_design_scan":
        from cmm.features.strain_design import optknock, robustknock

        named: dict[str, list[str]] = {}
        summaries: list[str] = []
        for method, solve in (("OptKnock", optknock), ("RobustKnock", robustknock)):
            try:
                proven = solve(
                    board.model,
                    board.product,
                    biomass=board.biomass,
                    max_knockouts=config.design_max_knockouts,
                    max_solutions=config.design_max_solutions,
                    min_growth=config.growth_floor,
                    seed=config.seed,
                )
            except Exception as error:
                # straindesign missing, or no MILP solver: a scan that cannot run is data.
                summaries.append(f"{method} could not run: {error}")
                continue
            designs = proven.designs[: config.design_max_solutions]
            summaries.append(f"{method} returned {len(proven.designs)} designs")
            for design in designs:
                board.scans.designs.append(
                    (method, tuple(design.knockouts), float(design.guaranteed_product))
                )
                for reaction_id in design.knockouts:
                    named.setdefault(reaction_id, []).append(
                        f"{method} deletes it in a {len(design.knockouts)}-knockout design "
                        f"reaching {design.guaranteed_product:.3g} guaranteed product"
                    )
        for reaction_id, notes in named.items():
            board.scans.design_notes[reaction_id] = notes[0]
        # Best guaranteed product first: that is the one ``adopt_best_design`` applies, and
        # ranking by guaranteed rather than maximum product is CMM's own rule for designs.
        board.scans.designs.sort(key=lambda item: (-item[2], item[0], item[1]))
        if named:
            summaries.append(
                "these reactions were added to the board: " + ", ".join(sorted(named))
            )
        return "; ".join(summaries)

    if name == "envelope_probe":
        from cmm.features.production import production_envelope

        try:
            envelope = production_envelope(board.model, board.product, points=12)
        except Exception as error:
            return f"the envelope could not be computed: {error}"
        note = (
            f"production envelope: the product reaches at most "
            f"{envelope.max_product:.4g} and growth at most {envelope.max_growth:.4g}"
        )
        board.scans.envelope_note = note
        return note

    return f"{name} is not a known scan"  # pragma: no cover - vocabulary is closed


def _moma_snapshot(
    board: _Board, config: JevConfig, use_linear: bool
) -> tuple[float | None, float | None, float | None]:
    """The minimal-adjustment state of the current design against the wild type."""

    if not config.run_moma:
        return (None, None, None)
    from cmm.features.comparison import moma

    try:
        result = moma(board.model, board.reference, linear=use_linear)
    except Exception as error:  # pragma: no cover - surfaced as a note, never fatal
        board.notes.append(f"MOMA failed on this design: {error}")
        return (None, None, None)
    if result.status != "optimal":
        return (None, None, None)
    return (
        float(result.fluxes.get(board.product, 0.0)),
        float(result.fluxes.get(board.biomass, 0.0)),
        result.distance,
    )


def _budget_exhausted(agent: JevClient, config: JevConfig) -> bool:
    return (
        agent.usage.calls >= config.max_decisions
        or agent.usage.cost_usd >= config.max_cost_usd
    )


def _safe_score(result: DecisionResult, key: str) -> float | None:
    answer = result.answers.get(key)
    if answer is None or answer.type != "score":
        return None
    try:
        return answer.score
    except (TypeError, ValueError):  # pragma: no cover - defensive
        return None


def _safe_noul(result: DecisionResult, key: str) -> float | None:
    answer = result.answers.get(key)
    if answer is None or answer.type != "noul":
        return None
    try:
        return float(answer.value)  # type: ignore[arg-type]
    except (TypeError, ValueError):  # pragma: no cover - defensive
        return None


def _record(
    transcript: list[Mapping[str, object]],
    round_index: int,
    tick_index: int,
    stage: str,
    request: object,
    result: DecisionResult,
) -> None:
    """Append one request/response pair to the audit transcript.

    The state payload is kept for the target stage only; repeating it for the action stage
    would double the transcript for no information, since the action stage sends the same
    state plus one field.
    """

    transcript.append(
        {
            "round": round_index,
            "tick": tick_index,
            "stage": stage,
            "request": request,
            "served_model": result.model,
            "request_id": result.request_id,
            "latency_s": round(result.latency_s, 4),
            "cost_usd": result.cost_usd,
            "answers": {
                key: {
                    "type": answer.type,
                    "value": answer.value,
                    "confidence": answer.confidence,
                    "probabilities": dict(answer.probabilities),
                }
                for key, answer in result.answers.items()
            },
        }
    )


def _resolve_exchange(model: Model, product: str) -> str:
    """Check the product names a real reaction, and say what is available when it does not."""

    if product in model.reactions:
        return product
    exchanges = sorted(r.id for r in model.exchanges)
    if not exchanges:
        raise JevWorkflowError(
            "this model has no exchange reactions, so production design is unavailable"
        )
    near = [rid for rid in exchanges if product.lower() in rid.lower()]
    hint = f" Did you mean one of: {', '.join(near[:8])}?" if near else ""
    raise JevWorkflowError(
        f"product {product!r} is not a reaction in this model.{hint}"
    )


def _objective_reaction_id(model: Model) -> str | None:
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            return str(reaction.id)
    return None


def detect_substrate(model: Model, fluxes: Mapping[str, float]) -> str | None:
    """The carbon source the model is actually consuming, from the wild-type solve.

    The yield needs a denominator, and asking for one separately is asking the user to repeat
    a fact the medium already fixed — in a second place, where it can disagree. The substrate
    is the organic exchange carrying the largest uptake; CO2 is excluded because a model
    fixing carbon dioxide is not being fed on it.
    """

    best: tuple[float, str] | None = None
    for reaction in model.exchanges:
        uptake = -float(fluxes.get(reaction.id, 0.0))
        if uptake <= 1e-9:
            continue
        carbon = 0
        for metabolite in reaction.metabolites:
            if str(getattr(metabolite, "formula", "") or "") == "CO2":
                carbon = 0
                break
            carbon = max(
                carbon, int((getattr(metabolite, "elements", None) or {}).get("C", 0))
            )
        if carbon <= 0:
            continue
        if best is None or uptake > best[0]:
            best = (uptake, reaction.id)
    return best[1] if best else None


def _theoretical_max_yield(
    model: Model, product: str, substrate: str | None, notes: list[str]
) -> float | None:
    """The yield ceiling for the scoreboard, or ``None`` with a note when it cannot be had."""

    if substrate is None:
        return None
    from cmm.features.production import theoretical_yield

    try:
        result = theoretical_yield(model, product, substrate)
    except Exception as error:
        notes.append(f"the theoretical yield could not be computed: {error}")
        return None
    if result.status != "optimal":
        notes.append(f"the theoretical yield solve is {result.status}")
        return None
    if abs(result.molar_yield) <= 1e-12:
        notes.append(
            f"the theoretical maximum yield of {product} is zero under this condition: "
            "the product is unreachable, and no intervention can change that"
        )
    return float(result.molar_yield)


__all__ = [
    "JevConfig",
    "JevResult",
    "JevWorkflowError",
    "RoundRecord",
    "TickRecord",
    "run_jev_design",
]
