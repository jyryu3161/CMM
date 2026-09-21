"""The screen: what JEV sees of the metabolic model before every move.

TypeSafe's Doom demo does not feed JEV pixels — it feeds structured facts (enemy positions,
distances, angles) and asks "what now?" ten times a second. This module is the equivalent for
a metabolic model: it turns a cobra model and a flux distribution into a compact JSON state
and a candidate list, each candidate carrying the handful of numbers that actually decide a
metabolic engineering move.

Those numbers are, per candidate reaction:

* **where the flux is** — its wild-type flux, and what it is carrying right now
* **how close it sits to the product** — shortest path through the metabolite graph to the
  target product, ignoring currency metabolites so the answer is not "2" for everything
* **what it does to the ATP pool** — net ATP stoichiometry, and the share of the model's
  total ATP production it accounts for
* **what it does to redox** — net NADH and NADPH stoichiometry
* **whether removing it kills the cell** — filled in once an ``essentiality_scan`` has run
* **whether it pulls the product** — FSEOF slope, filled in once an ``fseof_scan`` has run
* **what the literature says** — filled in only when web research is enabled

The JEV decision model has a 32K context, so size is the binding constraint, not cost. Every
record is one line, the candidate list is capped, and the budget is checked before a call
rather than discovered as a truncation afterwards.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace

from cobra import Model, Reaction

#: Metabolites that participate in most reactions and would collapse every graph distance to
#: 2 if left in. Matched on the id stem, so ``atp_c``/``atp_p``/``atp[c]`` all resolve.
CURRENCY_STEMS = frozenset(
    {
        "h",
        "h2o",
        "atp",
        "adp",
        "amp",
        "pi",
        "ppi",
        "nad",
        "nadh",
        "nadp",
        "nadph",
        "co2",
        "o2",
        "nh4",
        "pe",
        "coa",
        "q8",
        "q8h2",
        "fad",
        "fadh2",
        "so4",
        "hco3",
    }
)

#: Cofactor pools reported per reaction. The value is the stem of the "charged" member of the
#: pair, whose stoichiometric coefficient gives the net production of that pool.
ATP_STEM = "atp"
NADH_STEM = "nadh"
NADPH_STEM = "nadph"

#: A flux smaller than this reads as zero on the screen.
DISPLAY_EPSILON = 1e-9

#: How the board is divided between the three slates (see :func:`build_candidates`). The
#: remainder goes to the plain flux carriers.
#: How much of a literature answer goes onto the board. A web lookup returns around a
#: thousand characters; four of those would be a quarter of JEV's whole 32K context spent on
#: prose. The full text and its citations are kept in the run bundle, where length costs
#: nothing and a reader can check the source.
LITERATURE_EXCERPT_CHARS = 420

NEAR_SHARE = 0.45
COMPETING_SHARE = 0.35

#: How far from a secreted byproduct a reaction still counts as part of that branch.
BYPRODUCT_RADIUS = 2


def _flux(value: float) -> str:
    """Render a flux for the screen, without a sign on a zero."""

    return "0" if abs(value) <= DISPLAY_EPSILON else f"{value:+.3g}"


def _stem(metabolite_id: str) -> str:
    """The compartment-free stem of a metabolite id (``atp_c`` -> ``atp``, ``atp[c]`` -> ``atp``)."""

    base = metabolite_id.split("[")[0]
    if "_" in base:
        head, _, tail = base.rpartition("_")
        # Only strip a trailing single-letter compartment tag, not part of the name.
        if head and len(tail) <= 2:
            return head
    return base


def _net_coefficient(reaction: Reaction, stem: str) -> float:
    """Net stoichiometric coefficient of a cofactor pool in one reaction.

    Summed across compartments, because a reaction that consumes cytosolic ATP and produces
    periplasmic ATP is not ATP-neutral for the purposes of this screen.
    """

    return float(
        sum(
            coefficient
            for metabolite, coefficient in reaction.metabolites.items()
            if _stem(metabolite.id) == stem
        )
    )


@dataclass(frozen=True)
class CandidateEvidence:
    """Everything the screen says about one candidate reaction.

    Fields that are ``None`` are genuinely unknown rather than zero — ``essential`` is unknown
    until an ``essentiality_scan`` move has run, and ``fseof_slope`` until an ``fseof_scan``
    has. Rendering them as ``None`` is what lets JEV choose a LOOK move to find out.
    """

    reaction_id: str
    name: str
    subsystem: str
    genes: tuple[str, ...]
    reference_flux: float
    current_flux: float
    lower_bound: float
    upper_bound: float
    distance_to_product: int | None
    net_atp: float
    net_nadh: float
    net_nadph: float
    atp_production_share: float
    essential: bool | None = None
    fseof_slope: float | None = None
    #: Change in product flux CMM measured when this reaction was forced on, from an
    #: ``amplification_screen``. ``None`` until that scan has run.
    amplification_gain: float | None = None
    design_note: str = ""
    literature: str = ""
    citations: tuple[str, ...] = ()

    def to_record(self) -> str:
        """The one line JEV reads for this candidate.

        Written as short clauses rather than a JSON blob: the same content costs fewer tokens
        and the spike showed JEV reads it correctly. Every number carries its meaning, because
        a bare figure in a list is a figure that can be read as the wrong quantity.
        """

        parts = [
            f"wild-type flux {_flux(self.reference_flux)}",
            f"now {_flux(self.current_flux)}",
        ]
        if self.distance_to_product is None:
            parts.append("no pathway connects it to the product")
        elif self.distance_to_product == 1:
            parts.append("directly on the product's own exchange step")
        else:
            parts.append(f"{self.distance_to_product} steps from the product")
        if abs(self.net_atp) > 1e-9:
            share = (
                f", {self.atp_production_share:.0%} of total ATP production"
                if self.atp_production_share > 0.01
                else ""
            )
            verb = "makes" if self.net_atp > 0 else "spends"
            parts.append(f"{verb} {abs(self.net_atp):g} ATP per turnover{share}")
        if abs(self.net_nadh) > 1e-9:
            verb = "makes" if self.net_nadh > 0 else "consumes"
            parts.append(f"{verb} {abs(self.net_nadh):g} NADH")
        if abs(self.net_nadph) > 1e-9:
            verb = "makes" if self.net_nadph > 0 else "consumes"
            parts.append(f"{verb} {abs(self.net_nadph):g} NADPH")
        if self.essential is True:
            parts.append("ESSENTIAL: deleting it stops growth")
        elif self.essential is False:
            parts.append("not essential for growth")
        if self.fseof_slope is not None:
            direction = "rises" if self.fseof_slope > 0 else "falls"
            parts.append(
                f"FSEOF: flux {direction} ({self.fseof_slope:+.3g}) as product is forced up"
            )
        if self.amplification_gain is not None:
            if self.amplification_gain > 1e-6:
                parts.append(
                    f"MEASURED: forcing flux through it raises the product by "
                    f"{self.amplification_gain:+.4g}"
                )
            elif self.amplification_gain < -1e-6:
                parts.append(
                    f"measured: forcing flux through it LOWERS the product by "
                    f"{self.amplification_gain:+.4g}"
                )
            else:
                parts.append(
                    "measured: forcing flux through it changes the product by 0"
                )
        if self.design_note:
            parts.append(self.design_note)
        if self.subsystem:
            parts.append(f"subsystem {self.subsystem}")
        if self.literature:
            excerpt = " ".join(self.literature.split())
            if len(excerpt) > LITERATURE_EXCERPT_CHARS:
                excerpt = excerpt[: LITERATURE_EXCERPT_CHARS - 1].rstrip() + "\u2026"
            cited = f" [{len(self.citations)} sources]" if self.citations else ""
            parts.append(f"published evidence{cited}: {excerpt}")
        head = f"{self.name or self.reaction_id}"
        return f"{head} — " + "; ".join(parts)

    def to_label(self) -> str:
        """The short form that names this option in the answer space.

        A ``choice`` question's criteria define what may be answered; the evidence for
        weighing the options belongs in ``state.records``, which is the pattern the Decisions
        API is built around. Sending the whole record in both places duplicated 40% of the
        payload for nothing: measured against the live service on the same board, 3050 input
        tokens with the record repeated against 2500 with only a label, the same reaction
        chosen and the same shape of distribution either way. The saving is not the point on
        its own — the 32K context is the binding constraint on how many reactions fit on the
        board, so it buys candidates.
        """

        label = self.name or self.reaction_id
        return (
            f"{label}."
            if not self.design_note
            else f"{label}; the strain designer names it."
        )

    def to_row(self) -> dict[str, object]:
        """Flat export row for the run bundle."""

        return {
            "reaction_id": self.reaction_id,
            "name": self.name,
            "subsystem": self.subsystem,
            "genes": ";".join(self.genes),
            "reference_flux": self.reference_flux,
            "current_flux": self.current_flux,
            "distance_to_product": self.distance_to_product,
            "net_atp": self.net_atp,
            "net_nadh": self.net_nadh,
            "net_nadph": self.net_nadph,
            "atp_production_share": self.atp_production_share,
            "essential": self.essential,
            "fseof_slope": self.fseof_slope,
            "n_citations": len(self.citations),
        }


@dataclass(frozen=True)
class CofactorBalance:
    """Model-wide ATP and redox accounting at one flux distribution.

    Production and consumption are summed separately because their difference is ~0 in any
    balanced solution; what tells a metabolic engineer something is the *turnover* and which
    reaction dominates each side.
    """

    atp_production: float
    atp_consumption: float
    nadh_production: float
    nadh_consumption: float
    nadph_production: float
    nadph_consumption: float
    largest_atp_source: tuple[str, float] | None = None
    largest_nadh_sink: tuple[str, float] | None = None

    def to_payload(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "atp_turnover": round(self.atp_production, 4),
            "nadh_turnover": round(self.nadh_production, 4),
            "nadph_turnover": round(self.nadph_production, 4),
        }
        if self.largest_atp_source:
            rid, value = self.largest_atp_source
            payload["largest_atp_source"] = f"{rid} ({value:+.3g})"
        if self.largest_nadh_sink:
            rid, value = self.largest_nadh_sink
            payload["largest_nadh_sink"] = f"{rid} ({value:+.3g})"
        return payload


def cofactor_balance(model: Model, fluxes: Mapping[str, float]) -> CofactorBalance:
    """ATP/NADH/NADPH turnover at ``fluxes``, and the reaction dominating each side."""

    totals = {
        "atp_p": 0.0,
        "atp_c": 0.0,
        "nadh_p": 0.0,
        "nadh_c": 0.0,
        "nadph_p": 0.0,
        "nadph_c": 0.0,
    }
    best_atp: tuple[str, float] | None = None
    worst_nadh: tuple[str, float] | None = None

    for reaction in model.reactions:
        flux = float(fluxes.get(reaction.id, 0.0))
        if abs(flux) <= DISPLAY_EPSILON:
            continue
        for stem, key in (
            (ATP_STEM, "atp"),
            (NADH_STEM, "nadh"),
            (NADPH_STEM, "nadph"),
        ):
            rate = _net_coefficient(reaction, stem) * flux
            if rate > 0:
                totals[f"{key}_p"] += rate
            else:
                totals[f"{key}_c"] += -rate
            if (
                stem == ATP_STEM
                and rate > 0
                and (best_atp is None or rate > best_atp[1])
            ):
                best_atp = (reaction.id, rate)
            if (
                stem == NADH_STEM
                and rate < 0
                and (worst_nadh is None or rate < worst_nadh[1])
            ):
                worst_nadh = (reaction.id, rate)

    return CofactorBalance(
        atp_production=totals["atp_p"],
        atp_consumption=totals["atp_c"],
        nadh_production=totals["nadh_p"],
        nadh_consumption=totals["nadh_c"],
        nadph_production=totals["nadph_p"],
        nadph_consumption=totals["nadph_c"],
        largest_atp_source=best_atp,
        largest_nadh_sink=worst_nadh,
    )


def product_distances(model: Model, product_reaction_id: str) -> dict[str, int]:
    """Shortest number of reaction steps from each reaction to the product exchange.

    Breadth-first over the reaction/metabolite bipartite graph, with currency metabolites
    removed as edges. Keeping them in would make almost every reaction two steps from almost
    every other through the ATP or proton pool, which is true and useless.

    A reaction absent from the result is not connected to the product at all once currency
    metabolites are excluded.
    """

    product = model.reactions.get_by_id(product_reaction_id)
    by_metabolite: dict[str, list[str]] = {}
    for reaction in model.reactions:
        for metabolite in reaction.metabolites:
            if _stem(metabolite.id) in CURRENCY_STEMS:
                continue
            by_metabolite.setdefault(metabolite.id, []).append(reaction.id)

    distances = {product.id: 0}
    queue: deque[tuple[str, int]] = deque([(product.id, 0)])
    while queue:
        reaction_id, depth = queue.popleft()
        reaction = model.reactions.get_by_id(reaction_id)
        for metabolite in reaction.metabolites:
            if _stem(metabolite.id) in CURRENCY_STEMS:
                continue
            for neighbour in by_metabolite.get(metabolite.id, ()):
                if neighbour not in distances:
                    distances[neighbour] = depth + 1
                    queue.append((neighbour, depth + 1))
    return distances


def carbon_byproducts(
    model: Model, fluxes: Mapping[str, float], product_reaction_id: str
) -> tuple[str, ...]:
    """Exchange reactions secreting a carbon compound other than the product.

    These are where the carbon is going instead of into the product, which makes the
    reactions feeding them the competition. CO2 is excluded: it is the unavoidable end of
    oxidative metabolism, not a branch that can be redirected.
    """

    secreted: list[tuple[float, str]] = []
    for reaction in model.exchanges:
        if reaction.id == product_reaction_id:
            continue
        flux = float(fluxes.get(reaction.id, 0.0))
        if flux <= DISPLAY_EPSILON:  # uptake or idle
            continue
        carbon = 0
        for metabolite in reaction.metabolites:
            elements = getattr(metabolite, "elements", None) or {}
            carbon = max(carbon, int(elements.get("C", 0)))
            if _stem(metabolite.id) == "co2":
                carbon = 0
                break
        if carbon > 0:
            secreted.append((-flux, reaction.id))
    return tuple(rid for _, rid in sorted(secreted))


def build_candidates(
    model: Model,
    *,
    product_reaction_id: str,
    reference_fluxes: Mapping[str, float],
    current_fluxes: Mapping[str, float],
    excluded: Iterable[str] = (),
    pinned: Iterable[str] = (),
    limit: int = 24,
) -> tuple[CandidateEvidence, ...]:
    """Build the board for one tick: the reactions JEV may act on, shortlisted.

        A shortlist is needed because the state has to fit a 32K context, and it is computed by
        CMM rather than asked of JEV so the same model and flux state always produce the same
        board.

        **The board is composed from three slates, not one ranking.** A single blended score does
        not work, and the failure is instructive: ranked on proximity alone, an anaerobic
        succinate board contains only the succinate branch, and every move is "switch on another
        step of a pathway that has no reason to carry flux". Ranked on flux magnitude alone, it
        contains the respiratory chain and glycolysis, none of which can reach the product. Real
        strain design needs both halves — open the route *and* close what competes with it — so
        the board reserves places for each:

        ``near``
            Closest to the product through the metabolite graph. These are the moves that open
            the route.
        ``competing``
            Reactions feeding the carbon byproducts the strain is currently secreting, found by
            the same graph walk run backwards from each byproduct exchange. In anaerobic
            *E. coli* these are the ethanol, acetate and formate branches, and closing them is the
            textbook way to push carbon into succinate. Neither proximity nor flux magnitude
            surfaces them: ``ALCD2x`` is five steps from ``EX_succ_e`` and carries less flux than
            glycolysis, which cannot be touched at all.
        ``carriers``
            The remaining largest flux carriers, so the board is never blind to where the carbon
            actually is.

    ``pinned``
            Reactions named by a deterministic strain designer after a ``strain_design_scan``.
            These take their places first, because the slates above cannot reach them: the
            winning anaerobic succinate design deletes ``LDH_D`` and ``THD2``, neither of which
            carries any flux in the wild type, so neither appears near the product, near a
            secreted byproduct, or among the flux carriers. They are escape routes, and a board
            built from where the flux is today cannot see them.

        The slates are filled in that order, deduplicated, and truncated to ``limit``.

        **Only reactions with a gene association are candidates**, whenever the model carries GPRs
        at all. A move has to be something a laboratory could actually make, and an exchange
        reaction, a biomass pseudo-reaction, an ATP maintenance term or a passive diffusion step
        has no gene to delete or over-express. On ``e_coli_core`` this removes 26 of 95 reactions
        — every exchange, ``ATPM``, the biomass reaction and four gene-less transporters — and
        keeps every enzyme. Without the filter the board fills with "knock out acetate exchange",
        which raises the product on paper and cannot be built.
    """

    if limit < 2:
        raise ValueError("limit must be at least 2 for a choice question")
    excluded_ids = set(excluded)
    excluded_ids.add(product_reaction_id)
    genes_available = bool(model.genes)
    for reaction in model.reactions:
        if reaction.objective_coefficient != 0:
            excluded_ids.add(reaction.id)
        elif genes_available and not reaction.genes:
            excluded_ids.add(reaction.id)

    distances = product_distances(model, product_reaction_id)
    balance = cofactor_balance(model, reference_fluxes)
    atp_total = balance.atp_production or 1.0

    evidence: dict[str, CandidateEvidence] = {}
    for reaction in model.reactions:
        if reaction.id in excluded_ids:
            continue
        reference = float(reference_fluxes.get(reaction.id, 0.0))
        net_atp = _net_coefficient(reaction, ATP_STEM)
        evidence[reaction.id] = CandidateEvidence(
            reaction_id=reaction.id,
            name=reaction.name or reaction.id,
            subsystem=str(getattr(reaction, "subsystem", "") or ""),
            genes=tuple(sorted(gene.id for gene in reaction.genes)),
            reference_flux=reference,
            current_flux=float(current_fluxes.get(reaction.id, 0.0)),
            lower_bound=float(reaction.lower_bound),
            upper_bound=float(reaction.upper_bound),
            distance_to_product=distances.get(reaction.id),
            net_atp=net_atp,
            net_nadh=_net_coefficient(reaction, NADH_STEM),
            net_nadph=_net_coefficient(reaction, NADPH_STEM),
            atp_production_share=max(net_atp * reference, 0.0) / atp_total,
        )
    if not evidence:
        return ()

    def by_flux(candidate: CandidateEvidence) -> tuple[float, str]:
        return (-abs(candidate.reference_flux), candidate.reaction_id)

    near = sorted(
        (c for c in evidence.values() if c.distance_to_product is not None),
        key=lambda c: (c.distance_to_product, -abs(c.reference_flux), c.reaction_id),
    )
    # The competing slate: everything within two steps of a secreted carbon byproduct that
    # is not already close to the product. Two steps reaches the branch itself without
    # dragging in the shared trunk that feeds every branch.
    competing_ids: dict[str, int] = {}
    for byproduct in carbon_byproducts(model, reference_fluxes, product_reaction_id):
        for reaction_id, depth in product_distances(model, byproduct).items():
            if depth > BYPRODUCT_RADIUS:
                continue
            near_product = distances.get(reaction_id)
            if near_product is not None and near_product <= BYPRODUCT_RADIUS:
                continue
            competing_ids[reaction_id] = min(
                competing_ids.get(reaction_id, depth), depth
            )
    competing = sorted(
        (
            evidence[reaction_id]
            for reaction_id in competing_ids
            if reaction_id in evidence
        ),
        key=by_flux,
    )
    carriers = sorted(
        (c for c in evidence.values() if abs(c.reference_flux) > DISPLAY_EPSILON),
        key=by_flux,
    )

    pinned_slate = [
        evidence[reaction_id]
        for reaction_id in dict.fromkeys(pinned)
        if reaction_id in evidence
    ]
    quota_near = max(1, round(limit * NEAR_SHARE))
    quota_competing = max(1, round(limit * COMPETING_SHARE))
    board: dict[str, CandidateEvidence] = {}
    for slate, quota in (
        (pinned_slate, limit),
        (near, quota_near),
        (competing, quota_competing),
        (carriers, limit),
    ):
        taken = 0
        for candidate in slate:
            if taken >= quota or len(board) >= limit:
                break
            if candidate.reaction_id in board:
                continue
            board[candidate.reaction_id] = candidate
            taken += 1

    # Presentation order is part of the question, so it is deterministic and it leads with
    # the strongest evidence. A reaction a deterministic designer has *proved* forces the
    # product when deleted outranks one whose only recommendation is sitting two steps away:
    # ordered by distance alone, the designer's reactions landed at positions 12 to 24 of 24
    # and the agent, reading from the top, never reached them.
    pinned_ids = {candidate.reaction_id for candidate in pinned_slate}
    return tuple(
        sorted(
            board.values(),
            key=lambda c: (
                0 if c.reaction_id in pinned_ids else 1,
                c.distance_to_product if c.distance_to_product is not None else 10**6,
                -abs(c.reference_flux),
                c.reaction_id,
            ),
        )
    )


@dataclass(frozen=True)
class GameState:
    """One frame: the scoreboard, the cofactor accounting, the history and the board.

    :meth:`to_payload` is what goes over the wire as the Decisions ``state``.
    """

    product_reaction_id: str
    product_flux: float
    growth: float
    wild_type_product_flux: float
    wild_type_growth: float
    theoretical_max_yield: float | None
    molar_yield: float | None
    growth_floor: float
    balance: CofactorBalance
    candidates: tuple[CandidateEvidence, ...]
    active_interventions: tuple[str, ...] = ()
    #: Moves already tried that did not stick, as "REACTION: move" lines. Shown so the agent
    #: does not spend a tick re-proposing something the rules have already rejected.
    ruled_out: tuple[str, ...] = ()
    history: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    round_index: int = 1
    tick_index: int = 1
    ticks_left: int = 0
    interventions_used: int = 0
    max_interventions: int = 4

    def to_payload(self) -> dict[str, object]:
        """The compact JSON object sent as the Decisions ``state``."""

        scoreboard: dict[str, object] = {
            "product_reaction": self.product_reaction_id,
            "product_flux_now": round(self.product_flux, 5),
            "product_flux_wild_type": round(self.wild_type_product_flux, 5),
            "growth_now_per_h": round(self.growth, 5),
            "growth_wild_type_per_h": round(self.wild_type_growth, 5),
            "growth_floor_per_h": self.growth_floor,
        }
        if self.molar_yield is not None:
            scoreboard["molar_yield_mol_per_mol"] = round(self.molar_yield, 5)
        if self.theoretical_max_yield is not None:
            scoreboard["theoretical_max_yield_mol_per_mol"] = round(
                self.theoretical_max_yield, 5
            )

        payload: dict[str, object] = {
            "description": (
                "Metabolic engineering game state. The goal is to raise the product "
                f"exchange flux {self.product_reaction_id} as high as possible while "
                f"keeping the growth rate at or above {self.growth_floor} per hour. "
                "Each record below is one reaction you may act on."
            ),
            "scoreboard": scoreboard,
            "cofactor_balance": self.balance.to_payload(),
            "budget": {
                "round": self.round_index,
                "tick": self.tick_index,
                "ticks_left_this_round": self.ticks_left,
                "interventions_active": self.interventions_used,
                "max_interventions": self.max_interventions,
            },
            "active_interventions": list(self.active_interventions),
            "already_ruled_out": list(self.ruled_out),
            "history": list(self.history[-10:]),
            "records": [
                {"id": candidate.reaction_id, "record": candidate.to_record()}
                for candidate in self.candidates
            ],
        }
        if self.notes:
            payload["notes"] = list(self.notes)
        return payload

    def with_candidates(self, candidates: Sequence[CandidateEvidence]) -> "GameState":
        return replace(self, candidates=tuple(candidates))


@dataclass
class ScanCache:
    """Results of LOOK moves, carried forward so a scan is not paid for twice.

    A scan describes the model, not the tick, so its answer stays valid until an intervention
    changes the model — at which point :meth:`invalidate` drops the parts that moved.
    """

    essential: dict[str, bool] = field(default_factory=dict)
    fseof_slopes: dict[str, float] = field(default_factory=dict)
    #: reaction id -> product-flux change CMM measured when it was forced on, against the
    #: current design. This is the evidence an agent needs to choose an amplification, and it
    #: is measured rather than reasoned about: the reaction on the direct route to the product
    #: is often already saturated, and the one that pays is often a bypass no one would guess.
    amplification_gains: dict[str, float] = field(default_factory=dict)
    envelope_note: str = ""
    literature: dict[str, tuple[str, tuple[str, ...]]] = field(default_factory=dict)
    #: reaction id -> what the deterministic strain designer found about it. These are the
    #: reactions OptKnock and RobustKnock name, and they are forced onto the board because the
    #: ordinary slates cannot reach them: the winning anaerobic succinate design deletes
    #: ``LDH_D`` and ``THD2``, which carry no flux at all in the wild type. They are escape
    #: routes the cell would switch to once the obvious ones are shut, and a board built from
    #: where the flux is today is structurally blind to them.
    design_notes: dict[str, str] = field(default_factory=dict)
    #: The complete designs, best guaranteed product first, as
    #: ``(method, knockout reaction ids, guaranteed product)``. Kept whole because a design's
    #: deletions only pay off together — applied one at a time each looks worthless, and an
    #: agent judging a move by the product change it causes will never assemble one.
    designs: list[tuple[str, tuple[str, ...], float]] = field(default_factory=list)
    #: LOOK moves already run against the model as it currently stands. Every scan here is
    #: model-wide, not per-reaction — an FSEOF scan reports on the whole network — so "have I
    #: already asked this?" cannot be answered from one candidate's record. Tracking it per
    #: candidate let the agent re-run ``envelope_probe`` twenty ticks in a row.
    completed: set[str] = field(default_factory=set)

    def apply(
        self, candidates: Sequence[CandidateEvidence]
    ) -> tuple[CandidateEvidence, ...]:
        """Fold every cached scan result into a fresh candidate list."""

        enriched = []
        for candidate in candidates:
            literature, citations = self.literature.get(candidate.reaction_id, ("", ()))
            enriched.append(
                replace(
                    candidate,
                    essential=self.essential.get(candidate.reaction_id),
                    fseof_slope=self.fseof_slopes.get(candidate.reaction_id),
                    amplification_gain=self.amplification_gains.get(
                        candidate.reaction_id
                    ),
                    design_note=self.design_notes.get(candidate.reaction_id, ""),
                    literature=literature,
                    citations=citations,
                )
            )
        return tuple(enriched)

    def invalidate(self) -> None:
        """Forget the model-dependent scans after the model has been changed.

        Essentiality and FSEOF slopes are properties of the current bounds, so an
        intervention makes them stale. Literature is a property of the reaction and survives.
        """

        self.essential.clear()
        self.fseof_slopes.clear()
        self.amplification_gains.clear()
        self.envelope_note = ""
        self.design_notes.clear()
        self.designs.clear()
        self.completed.clear()

    def snapshot(self) -> "ScanCache":
        """A copy, so a change that is later undone can give the scans back.

        Without this, applying an intervention and immediately reverting it — which is what
        happens every time a move breaches the growth floor — destroys scan results that are
        still perfectly valid, because the model ended up back where it started. Observed
        effect: the agent re-ran the same two scans every tick and never progressed.
        """

        return ScanCache(
            essential=dict(self.essential),
            fseof_slopes=dict(self.fseof_slopes),
            amplification_gains=dict(self.amplification_gains),
            envelope_note=self.envelope_note,
            literature=dict(self.literature),
            design_notes=dict(self.design_notes),
            designs=list(self.designs),
            completed=set(self.completed),
        )

    def restore(self, snapshot: "ScanCache") -> None:
        """Adopt a previous snapshot's contents in place."""

        self.essential = dict(snapshot.essential)
        self.fseof_slopes = dict(snapshot.fseof_slopes)
        self.amplification_gains = dict(snapshot.amplification_gains)
        self.envelope_note = snapshot.envelope_note
        self.literature = dict(snapshot.literature)
        self.design_notes = dict(snapshot.design_notes)
        self.designs = list(snapshot.designs)
        self.completed = set(snapshot.completed)

    def knows(self, reaction_id: str) -> tuple[bool, bool]:
        """Whether essentiality and the FSEOF slope are already known for this reaction."""

        return (reaction_id in self.essential, reaction_id in self.fseof_slopes)
