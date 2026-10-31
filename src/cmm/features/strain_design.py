"""Computational strain design: OptKnock and RobustKnock.

OptKnock (Burgard 2003) searches for a small reaction-knockout set that couples a target
product to growth — at maximum growth the cell is *forced* to make the product. RobustKnock
(Tepper & Shlomi 2010) additionally requires the coupling to hold under the worst case (the
guaranteed minimum product at maximum growth, not the optimistic maximum).

The bilevel MILP search is delegated to the ``straindesign`` package; each returned design is
then evaluated for its optimistic and guaranteed product at maximum growth so OptKnock and
RobustKnock can be ranked on the right quantity.
"""

from __future__ import annotations

import contextlib
import io
import math
import time
import warnings
from collections.abc import Iterable
from dataclasses import dataclass
from dataclasses import field
from importlib.metadata import version as package_version
from numbers import Integral, Real
from typing import Literal, NamedTuple, cast

from cobra import Model
import pandas as pd

from cmm.core import solvers
from cmm.core.condition import Condition
from cmm.core.provenance import run_provenance
from cmm.features.production import (
    _actionable_reaction,
    _apply_condition,
    _condition_provenance,
    _record_condition,
)


# Gurobi's Seed parameter accepts integers through 2,000,000,000.  Using the same
# conservative range for the solver-neutral API also remains valid for the other
# straindesign backends while preventing a backend error after expensive preprocessing.
STRAIN_DESIGN_SEED_MAX = 2_000_000_000

# Wall-clock budget for one OptKnock/RobustKnock search.  straindesign takes no time limit
# unless it is given one, and its default 'best' approach can re-solve one MILP forever (see
# ``_search_designs``), so "no limit" is not a search that is merely slow.  Where the search
# terminated, a 1,440-reaction model took 4-5 minutes.
STRAIN_DESIGN_MAX_SECONDS = 3600.0

# What straindesign reports when its own clock ended the search.
_TIME_LIMIT_STATUSES = frozenset({"time_limit", "time_limit_w_sols"})


def validate_max_seconds(value: object, name: str = "max_seconds") -> None:
    """Accept ``None`` (no limit) or a positive, finite number of seconds."""

    if value is None:
        return
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ValueError(
            f"{name} must be a positive, finite number of seconds, or None for no limit"
        )


@dataclass(frozen=True)
class StrainDesign:
    """One reaction-knockout design and the product it forces at maximum growth.

    The reaction ids are a computational design, not a gene recipe.  Experimental use
    requires resolving every reaction's GPR, including isoenzymes and multisubunit enzymes.
    """

    knockouts: tuple[str, ...]
    growth: float
    max_product: float  # optimistic: max product at max growth
    guaranteed_product: float  # robust: min product at max growth

    @property
    def growth_coupled(self) -> bool:
        return self.guaranteed_product > 1e-6


@dataclass(frozen=True)
class StrainDesignResult:
    method: str
    product: str
    designs: tuple[StrainDesign, ...]
    metadata: dict[str, object] = field(default_factory=dict)

    def best(self) -> StrainDesign | None:
        return self.designs[0] if self.designs else None

    def to_frame(self) -> pd.DataFrame:
        """One row per design, in the method's documented rank order.

        ``knockouts`` is stored as a semicolon-delimited stable string for CSV portability;
        ``n_knockouts`` makes filtering independent of parsing that display field.
        """

        return pd.DataFrame(
            [
                {
                    "rank": rank,
                    "method": self.method,
                    "product": self.product,
                    "knockouts": ";".join(design.knockouts),
                    "n_knockouts": len(design.knockouts),
                    "growth": design.growth,
                    "max_product": design.max_product,
                    "guaranteed_product": design.guaranteed_product,
                    "growth_coupled": design.growth_coupled,
                }
                for rank, design in enumerate(self.designs, start=1)
            ],
            columns=[
                "rank",
                "method",
                "product",
                "knockouts",
                "n_knockouts",
                "growth",
                "max_product",
                "guaranteed_product",
                "growth_coupled",
            ],
        )


def _objective_reaction(model: Model) -> str:
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    raise ValueError("model has no objective reaction")


def _evaluate_design(
    model: Model, knockouts: Iterable[str], product: str, biomass: str
) -> StrainDesign | None:
    """Growth and optimistic/guaranteed product for a knockout set, at maximum growth."""

    kos = tuple(knockouts)
    with model:
        for rid in kos:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        growth = model.slim_optimize()
        if (
            growth is None or growth != growth or growth < 1e-6
        ):  # NaN/infeasible/no growth
            return None
        biomass_rxn = model.reactions.get_by_id(biomass)
        # Evaluate the product range on the *growth-optimal face*: exactly Tepper & Shlomi's
        # Omega', the set of growth-maximal flux states, over which RobustKnock minimises.
        # Relaxing the pin to a 99.9% face does not overstate the guarantee - measured, it
        # *lowers* the reported minimum every time, because a slightly sub-optimal face
        # contains extra flux states with less product. Pinning by equality is therefore the
        # faithful reading and the less conservative-looking one; both are reasons to keep it.
        biomass_rxn.bounds = (growth, growth)
        product_rxn = model.reactions.get_by_id(product)
        model.objective = product_rxn
        model.objective_direction = "max"
        max_product = model.slim_optimize()
        model.objective_direction = "min"
        guaranteed = model.slim_optimize()
    return StrainDesign(
        knockouts=kos,
        growth=float(growth),
        max_product=float(max_product) if max_product is not None else 0.0,
        guaranteed_product=float(guaranteed) if guaranteed is not None else 0.0,
    )


class _SearchOutcome(NamedTuple):
    ko_sets: list[tuple[str, ...]]
    n_candidates: int
    n_decompressed_designs: int
    status: str
    elapsed_s: float
    time_limit_reached: bool


def _search_designs(
    model: Model,
    product: str,
    biomass: str,
    *,
    method: Literal["optknock", "robustknock"],
    max_knockouts: int,
    max_solutions: int,
    min_growth: float,
    actionable_only: bool,
    seed: int,
    max_seconds: float | None,
) -> _SearchOutcome:
    try:
        import straindesign as sd
    except (
        ImportError
    ) as exc:  # pragma: no cover - straindesign is a declared dependency
        raise RuntimeError("strain design requires the 'straindesign' package") from exc

    # A product that cannot carry flux in the current medium has no coupling design — and
    # straindesign prunes such reactions and then crashes on the missing outer objective.
    with model:
        model.objective = model.reactions.get_by_id(product)
        model.objective_direction = "max"
        if (model.slim_optimize() or 0.0) < 1e-6:
            return _SearchOutcome([], 0, 0, "not_run_product_unreachable", 0.0, False)

    module_type = sd.OPTKNOCK if method == "optknock" else sd.ROBUSTKNOCK
    module = sd.SDModule(
        model,
        module_type,
        inner_objective=biomass,
        outer_objective=product,
        constraints=[f"{biomass} >= {min_growth}"],
    )
    search_kwargs: dict[str, object] = {}
    n_candidates = len(model.reactions)
    if actionable_only:
        # Burgard et al. restrict the candidate set to central metabolism. Without any
        # restriction the search returns designs that delete boundary reactions with no GPR
        # (EX_co2_e, EX_ac_e, EX_for_e, EX_etoh_e, EX_lac__D_e) — not realisable as gene
        # deletions. straindesign's ko_cost restricts candidates to the listed reactions.
        candidates = {
            rxn.id: 1.0
            for rxn in model.reactions
            if _actionable_reaction(model, rxn.id, product, biomass)
        }
        n_candidates = len(candidates)
        search_kwargs["ko_cost"] = candidates
    if max_seconds is not None:
        # Without this a search whose MILP stalls never returns.  straindesign 1.18 finds the
        # optimum, then solves a second MILP to minimise the interventions; when that one is
        # not optimal (Gurobi: "Model is infeasible", right after "max constraint violation
        # 4.75e-09 exceeds tolerance") the inner loop breaks without excluding the solution
        # it just found, the outer loop solves the identical MILP again, and ``max_solutions``
        # is never reached because nothing is ever added.  The limit is the only exit.
        search_kwargs["time_limit"] = float(max_seconds)
    started = time.monotonic()
    # straindesign is chatty; silence its solver logs.
    with (
        contextlib.redirect_stdout(io.StringIO()),
        contextlib.redirect_stderr(io.StringIO()),
    ):
        solutions = sd.compute_strain_designs(
            model,
            sd_modules=[module],
            max_cost=max_knockouts,
            max_solutions=max_solutions,
            compress=True,
            seed=seed,
            **search_kwargs,
        )
    elapsed_s = time.monotonic() - started
    reaction_sd = getattr(solutions, "reaction_sd", None) or []
    search_status = str(getattr(solutions, "status", "unknown"))
    # straindesign checks the remaining time only before a solve and keeps the status of the
    # last one, so the clock can end a search that still says 'optimal' - with no designs, or
    # fewer than a finished search would have found.  The elapsed time is the only witness.
    # It includes preprocessing, which straindesign's own clock does not, so a search that
    # finished just inside its limit can be counted as cut; that error is the safe one.
    time_limit_reached = search_status in _TIME_LIMIT_STATUSES or (
        max_seconds is not None
        and search_status == "optimal"
        and elapsed_s >= max_seconds
    )
    n_decompressed_designs = len(reaction_sd)
    designs: list[tuple[str, ...]] = []
    for design in reaction_sd:
        kos = tuple(sorted(rid for rid, coeff in design.items() if coeff <= 0))
        if kos and kos not in designs:
            designs.append(kos)
    return _SearchOutcome(
        designs,
        n_candidates,
        n_decompressed_designs,
        search_status,
        elapsed_s,
        time_limit_reached,
    )


def _strain_design(
    model: Model,
    product: str,
    *,
    method: Literal["optknock", "robustknock"],
    biomass: str | None,
    max_knockouts: int,
    max_solutions: int,
    min_growth: float,
    condition: Condition | None,
    actionable_only: bool,
    seed: int,
    max_seconds: float | None,
) -> StrainDesignResult:
    if max_knockouts < 1:
        raise ValueError("max_knockouts must be at least 1")
    if max_solutions < 1:
        raise ValueError("max_solutions must be at least 1")
    if min_growth < 0:
        raise ValueError("min_growth must be non-negative")
    validate_max_seconds(max_seconds)
    if (
        isinstance(seed, bool)
        or not isinstance(seed, Integral)
        or not 0 <= seed <= STRAIN_DESIGN_SEED_MAX
    ):
        raise ValueError(
            "seed must be an integer in the Gurobi-compatible range "
            f"[0, {STRAIN_DESIGN_SEED_MAX}]"
        )
    seed = int(seed)
    solvers.require("MILP", model.solver.interface, feature=method)
    biomass = biomass or _objective_reaction(model)
    provenance = run_provenance(
        model,
        seed=seed,
        method=method,
        product=product,
        biomass=biomass,
        max_knockouts=max_knockouts,
        max_solutions=max_solutions,
        min_growth=min_growth,
        condition=condition.name if condition is not None else None,
        actionable_only=actionable_only,
        strain_design_seed=seed,
        max_seconds=max_seconds,
    )
    with model:
        _apply_condition(model, condition=condition)
        _record_condition(provenance, _condition_provenance(model, condition))
        search = _search_designs(
            model,
            product,
            biomass,
            method=method,
            max_knockouts=max_knockouts,
            max_solutions=max_solutions,
            min_growth=min_growth,
            actionable_only=actionable_only,
            seed=seed,
            max_seconds=max_seconds,
        )
        evaluated = [
            d
            for d in (
                _evaluate_design(model, kos, product, biomass) for kos in search.ko_sets
            )
            if d is not None
        ]
    if search.time_limit_reached:
        limit = f"the {max_seconds:g} s time limit" if max_seconds else "a time limit"
        warnings.warn(
            f"{method}: straindesign was stopped by {limit} (backend status "
            f"'{search.status}', {search.elapsed_s:.0f} s elapsed) after finding "
            f"{len(search.ko_sets)} design(s); this is not a complete search. A search "
            "that found none may be re-solving one MILP forever (straindesign 1.18) or "
            "may only need a larger max_seconds.",
            UserWarning,
            stacklevel=3,
        )
    n_evaluated_designs = len(evaluated)
    if method == "robustknock":
        evaluated = [d for d in evaluated if d.growth_coupled]
    # A high optimistic maximum is not a coupling guarantee.  Both methods use the same
    # publication rank: worst-case product first, then optimistic product and growth, with
    # the reaction set as a deterministic final tie-break.
    evaluated.sort(
        key=lambda d: (
            -d.guaranteed_product,
            -d.max_product,
            -d.growth,
            d.knockouts,
        )
    )
    deduplicated = _deduplicate_designs(evaluated)
    provenance["parameters"] = {
        **cast(dict[str, object], provenance["parameters"]),
        "n_knockout_candidates": search.n_candidates,
        "strain_design_backend": "straindesign",
        "straindesign_version": package_version("straindesign"),
        "straindesign_search_status": search.status,
        "straindesign_search_elapsed_s": round(search.elapsed_s, 3),
        "straindesign_time_limit_reached": search.time_limit_reached,
        "straindesign_search_complete": (
            search.status == "optimal" and not search.time_limit_reached
        ),
        "straindesign_compress": True,
        "max_compressed_milp_solutions_requested": max_solutions,
        "max_solutions_semantics": (
            "caps compressed MILP solutions before straindesign network decompression; "
            "the returned reaction-design count can be larger"
        ),
        "n_compressed_milp_solutions": None,
        "n_compressed_milp_solutions_reason": (
            "straindesign does not expose the generated compressed-solution count"
        ),
        "n_decompressed_designs": search.n_decompressed_designs,
        "n_unique_decompressed_designs": len(search.ko_sets),
        "n_evaluated_designs": n_evaluated_designs,
        "n_method_eligible_designs": len(evaluated),
        "n_returned_designs_after_deduplication": len(deduplicated),
        "intervention_level": "reaction",
        "requires_gpr_resolution": True,
        # Backward-compatible historical keys.  ``n_milp_designs`` was a misnomer: it
        # counted evaluated decompressed designs, not compressed MILP solutions.
        "n_milp_designs": len(evaluated),
        "n_milp_designs_semantics": (
            "legacy alias for n_method_eligible_designs; not a compressed MILP count"
        ),
        "n_designs_after_deduplication": len(deduplicated),
    }
    return StrainDesignResult(
        method=method,
        product=product,
        designs=tuple(deduplicated),
        metadata=provenance,
    )


def _deduplicate_designs(designs: Iterable[StrainDesign]) -> list[StrainDesign]:
    """Keep the first design for each distinct knockout *set*, preserving rank order.

    ``max_solutions`` caps MILP solutions, not designs: straindesign decompresses one MILP
    solution into every member of the compressed reaction's coupled set, so the same
    intervention comes back repeatedly under different member ids. Collapsing on the
    knockout set as a set (order-insensitive) leaves one row per genuinely distinct design.
    """

    seen: set[frozenset[str]] = set()
    unique: list[StrainDesign] = []
    for design in designs:
        key = frozenset(design.knockouts)
        if key in seen:
            continue
        seen.add(key)
        unique.append(design)
    return unique


def optknock(
    model: Model,
    product: str,
    *,
    biomass: str | None = None,
    max_knockouts: int = 3,
    max_solutions: int = 5,
    min_growth: float = 0.05,
    condition: Condition | None = None,
    actionable_only: bool = True,
    seed: int = 0,
    max_seconds: float | None = STRAIN_DESIGN_MAX_SECONDS,
) -> StrainDesignResult:
    """OptKnock: knockout sets that maximize the product at maximum growth (optimistic).

    ``condition`` states the medium and aeration the search runs under, and is recorded in
    full in the provenance. Passing it is strongly preferred to conditioning the model
    beforehand: these two functions previously took no condition at all and depended
    silently on the caller's model state, which is how an aerobic design came to be
    published as the anaerobic answer.

    ``max_solutions`` caps **MILP solutions, not designs.** straindesign decompresses each
    MILP solution into every member of the compressed reaction's coupled set, so the number
    of designs returned is normally larger — ``max_solutions=5`` has returned 23. Designs
    are de-duplicated by knockout set before being returned; the raw and de-duplicated
    counts are both in the provenance.

    ``actionable_only`` restricts knockout candidates to reactions that have a GPR and are
    not boundary reactions, following Burgard et al.'s restriction to central metabolism.
    Set it to False to allow exchange knockouts, which are not realisable as gene deletions.

    ``seed`` is forwarded explicitly to straindesign's MILP backend.  Keeping it fixed pins
    the solver's search perturbation for reproducible reruns on the same solver/platform;
    omitting it no longer lets straindesign invent a different seed on every call.

    ``max_seconds`` is the wall-clock budget for the search (default one hour; ``None`` is
    no limit).  straindesign takes none by itself, and a search whose MILP stalls then never
    returns - a 1,440-reaction model sat in one for 24 hours.  When the budget ends the
    search, the designs found so far are returned with a ``UserWarning``,
    ``straindesign_time_limit_reached=True`` and ``straindesign_search_complete=False`` in
    the provenance, and the production workflow reports the artifact as partial.  A search
    cut at the limit is not a ranking: raise the budget or accept that the list is open.
    """

    return _strain_design(
        model,
        product,
        method="optknock",
        biomass=biomass,
        max_knockouts=max_knockouts,
        max_solutions=max_solutions,
        min_growth=min_growth,
        condition=condition,
        actionable_only=actionable_only,
        seed=seed,
        max_seconds=max_seconds,
    )


def robustknock(
    model: Model,
    product: str,
    *,
    biomass: str | None = None,
    max_knockouts: int = 3,
    max_solutions: int = 8,
    min_growth: float = 0.05,
    condition: Condition | None = None,
    actionable_only: bool = True,
    seed: int = 0,
    max_seconds: float | None = STRAIN_DESIGN_MAX_SECONDS,
) -> StrainDesignResult:
    """RobustKnock: maximize guaranteed product over all growth-maximal flux states.

    ``condition``, ``max_solutions``, ``actionable_only``, ``seed`` and ``max_seconds``
    behave exactly as in :func:`optknock`; in particular ``max_solutions`` caps MILP
    solutions, not designs.
    """

    return _strain_design(
        model,
        product,
        method="robustknock",
        biomass=biomass,
        max_knockouts=max_knockouts,
        max_solutions=max_solutions,
        min_growth=min_growth,
        condition=condition,
        actionable_only=actionable_only,
        seed=seed,
        max_seconds=max_seconds,
    )
