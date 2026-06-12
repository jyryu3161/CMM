"""Production-design services: theoretical yield, production envelope, FSEOF.

These answer "how do we make more of a target product?" — the metabolic-engineering side of
the workbench. They run on plain cobra models and return solver-neutral, serializable
results so the visualization layer and the GUI can render them without re-solving.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from cobra import Model
from cobra.flux_analysis import flux_variability_analysis as _cobra_fva
from cobra.flux_analysis import production_envelope as _cobra_envelope

_CO2_IDS = ("EX_co2_e", "EX_co2(e)", "EX_co2_e_")


@dataclass(frozen=True)
class ProductionYield:
    """Maximum product flux and yield per substrate at a fixed substrate uptake."""

    product: str
    substrate: str
    product_flux: float
    substrate_uptake: float
    molar_yield: float
    status: str
    aerobic: bool
    carbon_ceiling: float | None = None  # max mol product / mol substrate from substrate C
    co2_exchange: float = 0.0  # net CO2 exchange (negative = fixation)

    @property
    def co2_fixed(self) -> bool:
        """Whether net CO2 is being fixed (consumed), which can lift yield above the ceiling."""

        return self.co2_exchange < -1e-6

    @property
    def exceeds_carbon_ceiling(self) -> bool:
        if self.carbon_ceiling is None:
            return False
        return self.molar_yield > self.carbon_ceiling + 1e-6


@dataclass(frozen=True)
class EnvelopePoint:
    product_flux: float
    growth_min: float
    growth_max: float


@dataclass(frozen=True)
class ProductionEnvelope:
    product: str
    objective: str
    points: tuple[EnvelopePoint, ...]

    @property
    def max_growth(self) -> float:
        return max((p.growth_max for p in self.points), default=0.0)

    @property
    def max_product(self) -> float:
        return max((p.product_flux for p in self.points), default=0.0)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [(p.product_flux, p.growth_min, p.growth_max) for p in self.points],
            columns=["product_flux", "growth_min", "growth_max"],
        )


@dataclass(frozen=True)
class FseofResult:
    """Flux Scanning based on Enforced Objective Flux (Choi 2010)."""

    product: str
    biomass: str
    enforced_levels: tuple[float, ...]
    trends: pd.DataFrame  # index reaction_id; columns = enforced levels; + 'classification'
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        # frozen=True only blocks reassignment; copy the DataFrame so a caller mutating it
        # cannot corrupt this (or another caller's) cached result.
        object.__setattr__(self, "trends", self.trends.copy())

    def amplification_targets(self) -> list[str]:
        return self.trends.index[self.trends["classification"] == "amplify"].tolist()

    def knockout_targets(self) -> list[str]:
        return self.trends.index[self.trends["classification"] == "knockdown"].tolist()


@dataclass(frozen=True)
class FvseofResult:
    """Flux Variability Scanning based on Enforced Objective Flux.

    Like FSEOF but runs FVA at each enforced product level, so each reaction has a flux
    *range* per step. A reaction is a *robust* amplification target when its forced-minimum
    flux magnitude rises with enforced product — it cannot avoid carrying more flux — which
    FSEOF (single FBA) cannot distinguish.
    """

    product: str
    biomass: str
    enforced_levels: tuple[float, ...]
    mean: pd.DataFrame   # reaction x level: midpoint flux
    forced: pd.DataFrame  # reaction x level: forced minimum |flux| (0 if the range spans 0)
    classification: pd.Series  # reaction -> amplify / knockdown / none (on |mean| trend)
    robust: pd.Series  # reaction -> bool (forced-min |flux| monotonically rises)
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "mean", self.mean.copy())
        object.__setattr__(self, "forced", self.forced.copy())

    def amplification_targets(self) -> list[str]:
        return self.classification.index[self.classification == "amplify"].tolist()

    def knockout_targets(self) -> list[str]:
        return self.classification.index[self.classification == "knockdown"].tolist()

    def robust_targets(self) -> list[str]:
        """Amplification targets that are robustly forced (FVA minimum also rises)."""

        amp = set(self.amplification_targets())
        return [rid for rid in self.robust.index[self.robust] if rid in amp]


def _carbon_count(rxn) -> int:
    """Total carbon atoms across the reaction's metabolites (0 if none/unknown)."""

    total = 0
    for met in rxn.metabolites:
        try:
            total += int(met.elements.get("C", 0))
        except (AttributeError, TypeError):
            continue
    return total


def _detect_substrate(model: Model, unlimited: float = -999.0) -> str:
    """Pick the main carbon source: the limiting *carbon-containing* uptake exchange.

    The substrate is the carbon exchange with an explicitly restricted uptake bound (e.g.
    glucose at -10), not an unconstrained one and not a non-carbon nutrient (phosphate,
    sulfate, ...). Exchanges left at the default ~-1000 (CO2, H2O, NH4, O2) are availability
    constraints, not the carbon source.
    """

    carbon_limited = [
        rxn
        for rxn in model.exchanges
        if unlimited < rxn.lower_bound < 0.0 and _carbon_count(rxn) > 0
    ]
    if carbon_limited:
        return min(carbon_limited, key=lambda r: r.lower_bound).id
    # Fall back to the carbon exchange actually taken up at the wild-type optimum,
    # inside a context guard so the probe does not disturb the caller's model state.
    with model:
        solution = model.optimize()
    uptaken = [
        (rxn.id, solution.fluxes[rxn.id])
        for rxn in model.exchanges
        if solution.fluxes[rxn.id] < 0 and _carbon_count(rxn) > 0
    ]
    if not uptaken:
        raise ValueError("no carbon uptake exchange found; specify substrate explicitly")
    return min(uptaken, key=lambda x: x[1])[0]


def _set_anaerobic(model: Model) -> None:
    for oxygen_id in ("EX_o2_e", "EX_o2(e)", "EX_o2_e_"):
        if oxygen_id in model.reactions:
            model.reactions.get_by_id(oxygen_id).lower_bound = 0.0
            return


def _net_co2(solution, model: Model) -> float:
    for co2_id in _CO2_IDS:
        if co2_id in model.reactions:
            return float(solution.fluxes[co2_id])
    return 0.0


def theoretical_yield(
    model: Model,
    product: str,
    substrate: str | None = None,
    *,
    aerobic: bool = True,
) -> ProductionYield:
    """Maximum molar yield of ``product`` per ``substrate`` taken up.

    Maximizes the product exchange with the substrate uptake fixed exactly (so the yield
    denominator is deterministic). The reported yield can exceed the substrate's own carbon
    ceiling when the model fixes external CO2; ``carbon_ceiling`` and ``co2_exchange`` expose
    that so the number is never read as carbon-from-substrate alone.
    """

    with model:
        if not aerobic:
            _set_anaerobic(model)
        substrate = substrate or _detect_substrate(model)
        substrate_rxn = model.reactions.get_by_id(substrate)
        substrate_uptake = abs(substrate_rxn.lower_bound)
        if substrate_uptake <= 1e-9:
            # A closed substrate (lower_bound 0) gives a 0 denominator and a NaN yield. This
            # is plausible (e.g. a genome-scale medium with glucose shut, see rule R14), so
            # fail loudly with an actionable message instead of returning a silent NaN.
            raise ValueError(
                f"substrate '{substrate}' has no uptake capacity (lower_bound "
                f"{substrate_rxn.lower_bound}); open its uptake or choose another substrate "
                f"before computing a theoretical yield."
            )
        substrate_rxn.bounds = (-substrate_uptake, -substrate_uptake)
        model.objective = model.reactions.get_by_id(product)
        model.objective_direction = "max"
        solution = model.optimize()
        product_flux = float(solution.fluxes[product])
        molar_yield = product_flux / substrate_uptake

        substrate_carbon = _carbon_count(substrate_rxn)
        product_carbon = _carbon_count(model.reactions.get_by_id(product))
        carbon_ceiling = (
            substrate_carbon / product_carbon if product_carbon > 0 and substrate_carbon > 0 else None
        )
        return ProductionYield(
            product=product,
            substrate=substrate,
            product_flux=product_flux,
            substrate_uptake=substrate_uptake,
            molar_yield=molar_yield,
            status=solution.status,
            aerobic=aerobic,
            carbon_ceiling=carbon_ceiling,
            co2_exchange=_net_co2(solution, model),
        )


def production_envelope(
    model: Model,
    product: str,
    *,
    objective: str | None = None,
    substrate: str | None = None,
    aerobic: bool = True,
    points: int = 20,
) -> ProductionEnvelope:
    """Growth-vs-product envelope (phenotypic phase plane) for the product."""

    objective = objective or _objective_reaction(model)
    with model:
        if not aerobic:
            _set_anaerobic(model)
        substrate = substrate or _detect_substrate(model)
        frame = _cobra_envelope(
            model, [product], objective=objective, carbon_sources=substrate, points=points
        )
    envelope_points = tuple(
        EnvelopePoint(
            product_flux=float(row[product]),
            growth_min=float(row["flux_minimum"]),
            growth_max=float(row["flux_maximum"]),
        )
        for _, row in frame.iterrows()
    )
    return ProductionEnvelope(product=product, objective=objective, points=envelope_points)


def fseof(
    model: Model,
    product: str,
    biomass: str | None = None,
    *,
    n_steps: int = 10,
    fraction_min: float = 0.1,
    fraction_max: float = 0.9,
    aerobic: bool = True,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-3,
) -> FseofResult:
    """Scan enforced product flux and classify each reaction's flux trend (Choi 2010).

    At each enforced product level (a fraction of the theoretical maximum) the product flux
    is fixed and biomass is maximized. Reactions whose flux *magnitude* rises monotonically
    with enforced product are amplification targets; those whose magnitude falls are
    knockdown/knockout targets. Classifying on magnitude (not signed flux) keeps reactions
    that operate in the reverse direction — e.g. the reductive succinate pathway — correct.
    """

    biomass = biomass or _objective_reaction(model)
    yield_result = theoretical_yield(model, product, aerobic=aerobic)
    max_product = yield_result.product_flux
    levels = np.linspace(fraction_min, fraction_max, n_steps) * max_product

    rxn_ids = list(reactions) if reactions is not None else [r.id for r in model.reactions]
    columns: dict[float, dict[str, float]] = {}
    with model:
        if not aerobic:
            _set_anaerobic(model)
        product_rxn = model.reactions.get_by_id(product)
        model.objective = model.reactions.get_by_id(biomass)
        model.objective_direction = "max"
        for level in levels:
            product_rxn.bounds = (float(level), float(level))  # fix flux per Choi 2010
            solution = model.optimize()
            if solution.status != "optimal":
                columns[float(level)] = {rid: float("nan") for rid in rxn_ids}
                continue
            columns[float(level)] = {rid: float(solution.fluxes[rid]) for rid in rxn_ids}

    trends = pd.DataFrame(columns, index=rxn_ids)
    trends["classification"] = [
        _classify_trend(trends.loc[rid, list(columns.keys())].to_numpy(), tol)
        for rid in rxn_ids
    ]
    return FseofResult(
        product=product,
        biomass=biomass,
        enforced_levels=tuple(float(x) for x in levels),
        trends=trends,
        metadata={"max_product": max_product, "aerobic": aerobic},
    )


def _forced_min_magnitude(minimum: float, maximum: float, tol: float = 1e-9) -> float:
    """The smallest |flux| the reaction must carry given its feasible [min, max] range."""

    if minimum > tol:
        return minimum
    if maximum < -tol:
        return -maximum
    return 0.0  # the range spans zero, so the reaction can be off


def fvseof(
    model: Model,
    product: str,
    biomass: str | None = None,
    *,
    n_steps: int = 8,
    fraction_min: float = 0.1,
    fraction_max: float = 0.9,
    biomass_fraction: float = 0.95,
    aerobic: bool = True,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-3,
) -> FvseofResult:
    """Flux Variability Scanning based on Enforced Objective Flux (FSEOF + FVA per step).

    At each enforced product level the product flux is fixed and biomass is held at
    ``biomass_fraction`` of its maximum; FVA then gives every reaction's flux range. Reactions
    whose mean flux magnitude rises are amplification targets; those whose forced minimum also
    rises are *robust* targets (the reaction is compelled to carry more flux).
    """

    biomass = biomass or _objective_reaction(model)
    yield_result = theoretical_yield(model, product, aerobic=aerobic)
    max_product = yield_result.product_flux
    levels = np.linspace(fraction_min, fraction_max, n_steps) * max_product

    rxn_ids = list(reactions) if reactions is not None else [r.id for r in model.reactions]
    mean_cols: dict[float, dict[str, float]] = {}
    forced_cols: dict[float, dict[str, float]] = {}
    with model:
        if not aerobic:
            _set_anaerobic(model)
        product_rxn = model.reactions.get_by_id(product)
        model.objective = model.reactions.get_by_id(biomass)
        model.objective_direction = "max"
        for level in levels:
            level = float(level)
            with model:
                product_rxn.bounds = (level, level)
                try:
                    table = _cobra_fva(
                        model, reaction_list=rxn_ids,
                        fraction_of_optimum=biomass_fraction, processes=1,
                    )
                except Exception:
                    mean_cols[level] = {rid: float("nan") for rid in rxn_ids}
                    forced_cols[level] = {rid: float("nan") for rid in rxn_ids}
                    continue
            mean_cols[level] = {
                rid: 0.5 * (float(row["minimum"]) + float(row["maximum"]))
                for rid, row in table.iterrows()
            }
            forced_cols[level] = {
                rid: _forced_min_magnitude(float(row["minimum"]), float(row["maximum"]))
                for rid, row in table.iterrows()
            }

    mean = pd.DataFrame(mean_cols, index=rxn_ids)
    forced = pd.DataFrame(forced_cols, index=rxn_ids)
    levels_key = list(mean_cols.keys())
    classification = pd.Series(
        [_classify_trend(mean.loc[rid, levels_key].to_numpy(), tol) for rid in rxn_ids],
        index=rxn_ids,
    )
    robust = pd.Series(
        [_rises_monotonically(forced.loc[rid, levels_key].to_numpy(), tol) for rid in rxn_ids],
        index=rxn_ids,
    )
    return FvseofResult(
        product=product,
        biomass=biomass,
        # Use the actual scan keys so enforced_levels always matches the result columns
        # (a zero-yield product collapses all levels to a single 0.0).
        enforced_levels=tuple(float(x) for x in levels_key),
        mean=mean,
        forced=forced,
        classification=classification,
        robust=robust,
        metadata={"max_product": max_product, "biomass_fraction": biomass_fraction},
    )


def _rises_monotonically(values: np.ndarray, tol: float) -> bool:
    finite = values[np.isfinite(values)]
    if len(finite) < 2:
        return False
    return bool(finite[-1] - finite[0] > tol and np.all(np.diff(finite) >= -tol))


def _classify_trend(values: np.ndarray, tol: float) -> str:
    """Classify a reaction's flux-magnitude trend across the enforced-product scan."""

    finite = np.abs(values[np.isfinite(values)])
    if len(finite) < 2:
        return "none"
    start, end = finite[0], finite[-1]
    diffs = np.diff(finite)
    if end - start > tol and np.all(diffs >= -tol):
        return "amplify"
    if start - end > tol and np.all(diffs <= tol):
        return "knockdown"
    return "none"


def _objective_reaction(model: Model) -> str:
    for rxn in model.reactions:
        if rxn.objective_coefficient != 0:
            return rxn.id
    raise ValueError("model has no objective reaction")
