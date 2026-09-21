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
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
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
    NoAvailableAction,
    get_question_set,
    research_query,
)

if TYPE_CHECKING:  # the benchmark imports this module, so the type is compile-time only
    from cmm.jev.benchmark import BaselineRow

from cmm.jev.state import (
    CandidateEvidence,
    GameState,
    ScanCache,
    build_candidates,
    cofactor_balance,
)

#: How many times the identical move may land the identical way before the round is cut
#: short. Three is enough to be sure it is a loop and cheap enough not to matter if it is not.
_MAX_IDENTICAL_MOVES = 3

#: How many consecutive rounds may end without the agent changing anything before the run
#: stops. The state it is shown is identical each time, so a third round would ask the same
#: question and get the same answer at the same cost.
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


class JevWorkflowError(RuntimeError):
    """Raised when the run cannot proceed on scientific or configuration grounds."""


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JevConfig:
    """A complete, serializable invocation of a JEV design run.

    ``rounds`` x ``ticks_per_round`` bounds the number of moves; ``max_interventions`` bounds
    how many changes may be active at once, which is the quantity a wet-lab reader cares
    about — a design needing twelve edits is not the same proposal as one needing three.
    """

    model_path: str | Path
    product: str
    output_dir: str | Path | None = None
    substrate: str | None = None
    biomass: str | None = None
    solver: str | None = None
    medium: Medium | str | None = None
    condition: Condition | None = None
    organism: str = "Escherichia coli"

    # -- the game -----------------------------------------------------------
    rounds: int = 5
    ticks_per_round: int = 6
    max_interventions: int = 4
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
    #: After the game, score the agent's design against the deterministic methods on the same
    #: problem, every design evaluated the same way. An agent result with nothing to compare
    #: it to is not a result.
    run_baseline_comparison: bool = True

    # -- the agent ----------------------------------------------------------
    jev_model: str = "typesafe/jev-1.13"
    question_set: str = DEFAULT_QUESTION_SET
    enable_web_research: bool = False
    research_model: str = "openai/gpt-5.6-luna"
    #: A web lookup costs roughly $0.05 — about three hundred times a JEV decision, because
    #: the search results themselves are billed as input. Eight is a few dollars at most and
    #: still enough to cover the candidates that matter.
    max_research_calls: int = 8
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
        if self.ticks_per_round < 1:
            raise ValueError("ticks_per_round must be at least 1")
        if self.max_interventions < 1:
            raise ValueError("max_interventions must be at least 1")
        if self.growth_floor < 0:
            raise ValueError("growth_floor must be non-negative")
        if self.candidate_limit < 2:
            raise ValueError("candidate_limit must be at least 2 for a choice question")
        if self.max_decisions < 1:
            raise ValueError("max_decisions must be at least 1")
        if self.max_cost_usd <= 0:
            raise ValueError("max_cost_usd must be positive")
        if self.max_research_calls < 0:
            raise ValueError("max_research_calls must be non-negative")
        get_question_set(self.question_set)

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
            "rounds": self.rounds,
            "ticks_per_round": self.ticks_per_round,
            "max_interventions": self.max_interventions,
            "growth_floor": self.growth_floor,
            "candidate_limit": self.candidate_limit,
            "allow_look_actions": self.allow_look_actions,
            "design_max_knockouts": self.design_max_knockouts,
            "design_max_solutions": self.design_max_solutions,
            "seed_with_strain_design": self.seed_with_strain_design,
            "run_baseline_comparison": self.run_baseline_comparison,
            "run_moma": self.run_moma,
            "jev_model_requested": self.jev_model,
            "question_set": self.question_set,
            "enable_web_research": self.enable_web_research,
            "research_model": self.research_model if self.enable_web_research else None,
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
    """The checkpoint at the end of a round."""

    round_index: int
    n_ticks: int
    product_flux: float
    growth: float
    interventions: tuple[str, ...]
    improved: bool
    ended_early: bool

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
    #: reaction id -> (full literature answer, citation urls). The board carries only an
    #: excerpt; this is the whole thing, so a reader can check what the agent was told.
    literature: Mapping[str, tuple[str, tuple[str, ...]]] = field(default_factory=dict)
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
            "notes": list(self.notes),
            "usage": dict(self.usage),
        }

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
        """What the web lookup returned, in full, with its sources."""

        return pd.DataFrame(
            [
                {
                    "reaction_id": reaction_id,
                    "evidence": text,
                    "n_sources": len(urls),
                    "sources": "; ".join(urls),
                }
                for reaction_id, (text, urls) in sorted(self.literature.items())
            ],
            columns=["reaction_id", "evidence", "n_sources", "sources"],
        )

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
    previous_bounds: list[tuple[str, float, float]] = field(default_factory=list)
    history: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    scans: ScanCache = field(default_factory=ScanCache)
    research_calls: int = 0
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
    #: reaction id -> the change in product flux measured when that intervention was applied.
    #: Shown next to each active intervention so dead weight is visible: an intervention that
    #: bought nothing still occupies one of the design's places, and withdrawing it is a real
    #: move the agent can only choose if it can see the cost.
    contribution: dict[str, float] = field(default_factory=dict)

    def apply(self, intervention: Intervention) -> None:
        reaction = self.model.reactions.get_by_id(intervention.reaction_id)
        self.previous_bounds.append(
            (
                intervention.reaction_id,
                float(reaction.lower_bound),
                float(reaction.upper_bound),
            )
        )
        self.scan_stack.append(self.scans.snapshot())
        reaction.bounds = (intervention.lower_bound, intervention.upper_bound)
        self.interventions.append(intervention)
        self.scans.invalidate()

    def undo(self) -> Intervention | None:
        if not self.interventions:
            return None
        intervention = self.interventions.pop()
        self.contribution.pop(intervention.reaction_id, None)
        reaction_id, lower, upper = self.previous_bounds.pop()
        self.model.reactions.get_by_id(reaction_id).bounds = (lower, upper)
        # The model is back where it was, so the scans taken before the change are valid
        # again. Clearing them here made the agent re-run the same scan every tick.
        if self.scan_stack:
            self.scans.restore(self.scan_stack.pop())
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
) -> JevResult:
    """Run the JEV game loop and, when ``output_dir`` is set, write a run bundle.

    ``client`` is injectable so a test can play the whole game against a scripted agent with
    no network. ``on_tick`` is called after every frame with the tick record and the flux
    distribution it produced, which is how the desktop app redraws the flux map mid-round.
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

    theoretical = _theoretical_max_yield(model, product, config.substrate, notes)

    board = _Board(model=model, reference=reference, product=product, biomass=biomass)

    if config.seed_with_strain_design:
        # Before the first move: hand the agent what the deterministic designer already knows.
        # The reactions it names go on the board with the guaranteed product they buy, which
        # is the only way an escape route carrying no flux today can ever be considered.
        seeded = _run_scan(board, "strain_design_scan", (), config)
        board.scans.completed.add("strain_design_scan")
        notes.append(f"strain design seeded the board before the first move: {seeded}")

    question_set = get_question_set(config.question_set)
    ticks: list[TickRecord] = []
    round_records: list[RoundRecord] = []
    flux_frames: list[Mapping[str, float]] = [dict(wild_type.fluxes)]
    transcript: list[Mapping[str, object]] = []

    idle_rounds = 0
    best_product = wild_product
    best_growth = wild_growth
    best_interventions: tuple[Intervention, ...] = ()

    current_product, current_growth = wild_product, wild_growth
    stop_run = False

    for round_index in range(1, config.rounds + 1):
        round_start_product = current_product
        ended_early = False
        ticks_this_round = 0
        last_signature: tuple[object, ...] | None = None
        repeats = 0

        for tick_index in range(1, config.ticks_per_round + 1):
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
                    ],
                    pinned=board.scans.design_notes,
                    limit=config.candidate_limit,
                )
            )
            if not candidates:
                notes.append(
                    f"round {round_index} tick {tick_index}: no candidate reactions are "
                    "left to act on"
                )
                ended_early = True
                break

            state = GameState(
                product_reaction_id=product,
                product_flux=current_product,
                growth=current_growth,
                wild_type_product_flux=wild_product,
                wild_type_growth=wild_growth,
                theoretical_max_yield=theoretical,
                molar_yield=None,
                growth_floor=config.growth_floor,
                balance=cofactor_balance(board.model, current_fluxes),
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
                ticks_left=config.ticks_per_round - tick_index,
                interventions_used=len(board.interventions),
                max_interventions=config.max_interventions,
            )

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
            )
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

            if tick.outcome == "end_round":
                ended_early = True
                break

        acted = any(
            tick.outcome in ("applied", "undone_by_agent")
            for tick in ticks[-ticks_this_round:]
        )
        idle_rounds = 0 if acted else idle_rounds + 1

        round_records.append(
            RoundRecord(
                round_index=round_index,
                n_ticks=ticks_this_round,
                product_flux=current_product,
                growth=current_growth,
                interventions=board.active_labels(),
                improved=current_product > round_start_product + 1e-9,
                ended_early=ended_early,
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
    }

    # The design the run ended on, captured before the comparison below strips the board
    # back to the unmodified model. Reading it afterwards returned an empty design.
    final_interventions = tuple(board.interventions)

    baselines: tuple["BaselineRow", ...] = ()
    if config.run_baseline_comparison:
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
        literature=dict(board.scans.literature),
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
) -> TickRecord:
    """Ask, execute, measure. Returns the frame; the caller re-solves and redraws."""

    payload = state.to_payload()
    allow_undo = bool(board.interventions)
    allow_look = config.allow_look_actions

    room = config.max_interventions - len(board.interventions)
    best_design = next(
        (
            design
            for design in board.scans.designs
            if len(design[1]) <= room
            and not set(design[1]) & {i.reaction_id for i in board.interventions}
        ),
        None,
    )
    target_questions = question_set.target_question(
        candidates,
        product=board.product,
        growth_floor=config.growth_floor,
        allow_undo=allow_undo,
        allow_look=allow_look,
        design_full=room <= 0,
        proven_design=(
            f"{best_design[0]} deletes {', '.join(best_design[1])} for "
            f"{best_design[2]:.3g} guaranteed product"
            if best_design
            else ""
        ),
    )
    stage1 = agent.decide(payload, target_questions)
    _record(transcript, round_index, tick_index, "target", payload, stage1)
    target_answer = stage1[TARGET_KEY]
    target = target_answer.choice
    ranking = target_answer.ranked()
    cost = stage1.cost_usd
    latency = stage1.latency_s

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

    # -- optional literature lookup, before the action is chosen -------------
    if (
        config.enable_web_research
        and board.research_calls < config.max_research_calls
        and candidate.reaction_id not in board.scans.literature
    ):
        board.research_calls += 1
        try:
            text, urls = agent.web_research(
                research_query(candidate, board.product, config.organism)
            )
        except JevTransportError as error:
            board.notes.append(
                f"web research for {candidate.reaction_id} failed: {error}"
            )
        else:
            board.scans.literature[candidate.reaction_id] = (text, tuple(urls))
            candidate = replace(candidate, literature=text, citations=tuple(urls))

    # -- stage two ----------------------------------------------------------
    # Moves already ruled out on this reaction, plus every model-wide scan already run
    # against the current bounds. Both would spend a tick to arrive back where we are.
    blocked = set(board.failed_moves.get(candidate.reaction_id, ()))
    blocked |= board.scans.completed
    try:
        action_questions = question_set.action_question(
            candidate,
            product=board.product,
            growth_floor=config.growth_floor,
            allow_look=allow_look,
            exclude=blocked,
        )
    except NoAvailableAction as error:
        # Nothing is left to try here. Mark the whole reaction spent so the next tick's board
        # is built without it, instead of offering a dead end again.
        for spent in applicable_actions(candidate.reference_flux):
            board.record_failure(candidate.reaction_id, spent.name)
        return frame(
            action=None,
            outcome="not_applicable",
            reason=str(error),
            product_flux=state.product_flux,
            growth=state.growth,
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
    action_name = action_answer.choice
    benefit = _safe_score(stage2, BENEFIT_KEY)
    risk = _safe_noul(stage2, RISK_KEY)
    action = ACTION_CATALOGUE.get(action_name)
    if action is None:  # pragma: no cover - vocabulary is closed
        return frame(
            action=action_name,
            outcome="not_applicable",
            reason=f"{action_name!r} is not a known move",
            benefit=benefit,
            risk=risk,
            action_confidence=action_answer.confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

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
            action_confidence=action_answer.confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    # -- ACT ----------------------------------------------------------------
    if len(board.interventions) >= config.max_interventions:
        # Not recorded as a failure: the move may be fine, there is simply no room for it.
        return frame(
            action=action.name,
            outcome="not_applicable",
            reason=(
                f"the design already carries the maximum of {config.max_interventions} "
                "interventions; withdraw one before adding another"
            ),
            benefit=benefit,
            risk=risk,
            action_confidence=action_answer.confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    try:
        intervention = build_intervention(
            board.model,
            candidate.reaction_id,
            action,
            reference_flux=candidate.reference_flux,
        )
    except ActionNotApplicable as error:
        board.record_failure(candidate.reaction_id, action.name)
        return frame(
            action=action.name,
            outcome="not_applicable",
            reason=str(error),
            benefit=benefit,
            risk=risk,
            action_confidence=action_answer.confidence,
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
            action_confidence=action_answer.confidence,
            status=solution.status,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    new_product = float(solution.fluxes.get(board.product, 0.0))
    new_growth = float(solution.fluxes.get(board.biomass, 0.0))
    if new_growth < config.growth_floor:
        board.undo()
        board.record_failure(candidate.reaction_id, action.name)
        return frame(
            action=action.name,
            outcome="reverted_growth_floor",
            reason=(
                f"growth would fall to {new_growth:.4g} per hour, below the floor of "
                f"{config.growth_floor}"
            ),
            intervention=intervention,
            benefit=benefit,
            risk=risk,
            action_confidence=action_answer.confidence,
            product_flux=state.product_flux,
            growth=state.growth,
        )

    moma = _moma_snapshot(board, config, use_linear_moma)
    delta = new_product - state.product_flux
    board.contribution[intervention.reaction_id] = delta
    if delta > 1e-9:
        verdict = f"product rose by {delta:+.4g}"
        board.clear_failures()
    elif delta < -1e-9:
        verdict = f"product FELL by {delta:+.4g}"
        board.clear_failures()
    elif intervention.mode == "amplification":
        # The forced flux was met without any of it reaching the product — the network
        # satisfied the constraint internally, typically through a cycle. Naming this is the
        # difference between the agent learning something and repeating the move elsewhere.
        verdict = (
            "product unchanged: the forced flux was consumed inside the network and none "
            "of it reached the product"
        )
        board.record_failure(candidate.reaction_id, action.name)
    else:
        verdict = "product unchanged"
    board.clear_failures()
    return frame(
        action=action.name,
        outcome="applied",
        reason=verdict,
        intervention=intervention,
        benefit=benefit,
        risk=risk,
        action_confidence=action_answer.confidence,
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
    deterministic result reachable, and it leaves the interesting question open: whether an
    amplification on top of a proven design beats the design alone, which is a question the
    designer itself cannot answer.

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
                reference_flux=float(board.reference.get(reaction_id, 0.0)),
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
    # whole remaining opportunity: OptKnock and RobustKnock search deletions only. Without
    # this line, eight runs out of eight adopted the design and ended the round immediately,
    # while forcing flux through the glyoxylate shunt on top of it reaches 10.76 against the
    # design's 9.91.
    note = (
        f"The active design came from {method}, which searches deletions only — its "
        "formulation cannot express forcing more flux through a reaction. An amplification "
        "or a knockdown on top of it is a move it could not have considered, and is the only "
        "kind of move left that might improve on a proven design."
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

    if name == "amplification_screen":
        # One solve per candidate, against the design as it stands. This is CMM answering
        # the question the agent would otherwise have to guess at — and the guess is
        # systematically wrong in a way worth naming: on anaerobic succinate the agent
        # reaches for fumarate reductase, the direct product-forming step, which is already
        # saturated and buys nothing, while the glyoxylate shunt buys 8.6%.
        baseline = _solve(board.model)
        if baseline.status != "optimal":
            return "the screen needs a feasible starting point; the model is not"
        before = float(baseline.fluxes.get(board.product, 0.0))
        measured = 0
        best: tuple[str, float] | None = None
        for candidate in candidates:
            action = (
                ACTION_CATALOGUE["force_on_low"]
                if abs(candidate.reference_flux) <= 1e-9
                else ACTION_CATALOGUE["amplify_2x"]
            )
            try:
                trial = build_intervention(
                    board.model,
                    candidate.reaction_id,
                    action,
                    reference_flux=candidate.reference_flux,
                )
            except ActionNotApplicable:
                continue
            reaction = board.model.reactions.get_by_id(candidate.reaction_id)
            saved = reaction.bounds
            reaction.bounds = (trial.lower_bound, trial.upper_bound)
            solution = _solve(board.model)
            reaction.bounds = saved
            if solution.status != "optimal":
                continue
            if float(solution.fluxes.get(board.biomass, 0.0)) < config.growth_floor:
                # Reported as no gain rather than omitted: a move that kills the strain is
                # not a move, and leaving the row blank would read as "not yet measured".
                board.scans.amplification_gains[candidate.reaction_id] = 0.0
                measured += 1
                continue
            gain = float(solution.fluxes.get(board.product, 0.0)) - before
            board.scans.amplification_gains[candidate.reaction_id] = gain
            measured += 1
            if best is None or gain > best[1]:
                best = (candidate.reaction_id, gain)
        if best is None:
            return f"measured {measured} reactions; none of them raises the product"
        return (
            f"measured {measured} reactions; the best is {best[0]} at {best[1]:+.4g} "
            "product flux"
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
