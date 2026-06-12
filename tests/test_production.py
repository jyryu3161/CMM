from __future__ import annotations

import pytest
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)

GLC = "EX_glc__D_e"
SUCC = "EX_succ_e"
BIOMASS = "Biomass_Ecoli_core"


def test_theoretical_yield_succinate_anaerobic(ecoli_core):
    result = theoretical_yield(ecoli_core, SUCC, aerobic=False)
    assert result.status == "optimal"
    assert result.substrate == GLC
    assert result.substrate_uptake == pytest.approx(10.0, abs=1e-6)
    # Known anaerobic succinate yield on e_coli_core ~1.39 mol/mol.
    assert result.molar_yield == pytest.approx(1.391, abs=0.02)


def test_theoretical_yield_aerobic_higher(ecoli_core):
    aerobic = theoretical_yield(ecoli_core, SUCC, aerobic=True)
    anaerobic = theoretical_yield(ecoli_core, SUCC, aerobic=False)
    assert aerobic.molar_yield == pytest.approx(1.638, abs=0.02)
    assert aerobic.molar_yield > anaerobic.molar_yield


def test_theoretical_yield_does_not_mutate_model(ecoli_core):
    before = ecoli_core.reactions.EX_o2_e.lower_bound
    theoretical_yield(ecoli_core, SUCC, aerobic=False)
    assert ecoli_core.reactions.EX_o2_e.lower_bound == before
    # Objective restored to biomass.
    assert ecoli_core.slim_optimize() == pytest.approx(0.8739, abs=1e-3)


def test_theoretical_yield_raises_on_closed_substrate(ecoli_core):
    # A closed substrate (lower_bound 0) would divide by zero -> NaN; fail loudly instead.
    ecoli_core.reactions.get_by_id(GLC).bounds = (0.0, 1000.0)
    with pytest.raises(ValueError, match="no uptake capacity"):
        theoretical_yield(ecoli_core, SUCC, substrate=GLC, aerobic=True)


def test_production_envelope_is_growth_coupled(ecoli_core):
    envelope = production_envelope(ecoli_core, SUCC, points=15)
    assert envelope.max_growth == pytest.approx(0.8739, abs=1e-2)
    frame = envelope.to_frame()
    # Growth falls as enforced succinate rises (growth-coupled tradeoff).
    assert frame["growth_max"].iloc[0] > frame["growth_max"].iloc[-1]
    assert frame["product_flux"].iloc[-1] > frame["product_flux"].iloc[0]


def test_production_envelope_respects_anaerobic_flag(ecoli_core):
    aerobic = production_envelope(ecoli_core, SUCC, aerobic=True, points=8)
    before_o2 = ecoli_core.reactions.EX_o2_e.lower_bound
    anaerobic = production_envelope(ecoli_core, SUCC, aerobic=False, points=8)

    assert anaerobic.max_growth < aerobic.max_growth
    assert anaerobic.max_product < aerobic.max_product
    assert ecoli_core.reactions.EX_o2_e.lower_bound == before_o2


def test_fseof_finds_known_succinate_targets(ecoli_core):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, aerobic=False)
    amplification = set(result.amplification_targets())
    # FRD7 (fumarate reductase) and PPC (PEP carboxylase) are textbook succinate
    # amplification targets and must appear.
    assert "FRD7" in amplification
    assert "PPC" in amplification
    assert SUCC in amplification


def test_fseof_trends_shape(ecoli_core):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=6, aerobic=False)
    assert result.trends.shape[0] == len(ecoli_core.reactions)
    assert "classification" in result.trends.columns
    assert len(result.enforced_levels) == 6


_FVSEOF_RXNS = ["FRD7", "PPC", "MDH", "FUM", "SUCCt3", "EX_succ_e", "CS", "ACALD"]


def test_fvseof_finds_robust_succinate_targets(ecoli_core):
    result = fvseof(
        ecoli_core, SUCC, BIOMASS, n_steps=4, aerobic=False, reactions=_FVSEOF_RXNS
    )
    amplification = set(result.amplification_targets())
    # The reductive succinate pathway must be amplification targets...
    assert {"FRD7", "FUM", "MDH"} <= amplification
    # ...and robustly forced (FVA minimum rises, not just the mean).
    robust = set(result.robust_targets())
    assert "FRD7" in robust
    assert robust <= amplification  # robust targets are a subset of amplification targets


def test_fvseof_ranges_shape_and_columns(ecoli_core):
    result = fvseof(ecoli_core, SUCC, BIOMASS, n_steps=3, aerobic=False, reactions=_FVSEOF_RXNS)
    assert result.mean.shape == (len(_FVSEOF_RXNS), 3)
    assert result.forced.shape == (len(_FVSEOF_RXNS), 3)
    # The forced minimum magnitude is never negative (it is a |flux| lower bound).
    assert (result.forced.to_numpy() >= -1e-9).all()
    assert len(result.enforced_levels) == 3


def test_fvseof_zero_yield_product_is_consistent(ecoli_core):
    # A blocked product collapses every enforced level to 0; result stays shape-consistent.
    ecoli_core.reactions.EX_succ_e.bounds = (0.0, 0.0)
    result = fvseof(ecoli_core, SUCC, BIOMASS, n_steps=4, aerobic=False, reactions=["FRD7", "PPC"])
    assert result.metadata["max_product"] == pytest.approx(0.0, abs=1e-9)
    assert len(result.enforced_levels) == result.mean.shape[1]  # no level/column mismatch


def test_fvseof_does_not_mutate_model(ecoli_core):
    objective = str(ecoli_core.objective.expression)
    growth = ecoli_core.slim_optimize()
    fvseof(ecoli_core, SUCC, BIOMASS, n_steps=3, aerobic=False, reactions=_FVSEOF_RXNS)
    assert str(ecoli_core.objective.expression) == objective
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_fseof_classifies_reductive_pathway_as_amplify(ecoli_core):
    # MDH/FUM operate in reverse (negative flux) whose magnitude rises with enforced
    # succinate; classifying on magnitude (not signed flux) must label them amplify, and
    # the ethanol pathway (magnitude falls to zero) as knockdown.
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, aerobic=False)
    assert result.trends.loc["MDH", "classification"] == "amplify"
    assert result.trends.loc["FUM", "classification"] == "amplify"
    assert result.trends.loc["ALCD2x", "classification"] == "knockdown"


def test_fseof_knockout_targets_are_disjoint_and_nonempty(ecoli_core):
    result = fseof(ecoli_core, SUCC, BIOMASS, n_steps=8, aerobic=False)
    amplify = set(result.amplification_targets())
    knockdown = set(result.knockout_targets())
    assert knockdown  # there are reactions to turn down
    assert amplify.isdisjoint(knockdown)


def test_fseof_does_not_mutate_model(ecoli_core):
    objective = str(ecoli_core.objective.expression)
    direction = ecoli_core.objective_direction
    succ_bounds = ecoli_core.reactions.EX_succ_e.bounds
    growth = ecoli_core.slim_optimize()
    fseof(ecoli_core, SUCC, BIOMASS, n_steps=5, aerobic=False)
    assert str(ecoli_core.objective.expression) == objective
    assert ecoli_core.objective_direction == direction
    assert ecoli_core.reactions.EX_succ_e.bounds == succ_bounds
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_production_envelope_does_not_mutate_model(ecoli_core):
    growth = ecoli_core.slim_optimize()
    production_envelope(ecoli_core, SUCC, points=8)
    assert ecoli_core.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_fseof_results_hold_independent_dataframes(ecoli_core):
    # FseofResult.__post_init__ copies the trends frame, so mutating one result's frame
    # cannot corrupt another's (the false-immutability finding).
    r1 = fseof(ecoli_core, SUCC, BIOMASS, n_steps=5, aerobic=False)
    r2 = fseof(ecoli_core, SUCC, BIOMASS, n_steps=5, aerobic=False)
    r1.trends.loc["FRD7", "classification"] = "TAMPERED"
    assert r2.trends.loc["FRD7", "classification"] != "TAMPERED"


def test_theoretical_yield_discloses_co2_fixation(ecoli_core):
    aerobic = theoretical_yield(ecoli_core, SUCC, aerobic=True)
    # Aerobic succinate yield (1.64) exceeds the glucose carbon ceiling (1.5) only via CO2
    # fixation; the result must expose that rather than presenting 1.64 as carbon-from-glucose.
    assert aerobic.carbon_ceiling == pytest.approx(1.5, abs=1e-6)
    assert aerobic.exceeds_carbon_ceiling
    assert aerobic.co2_fixed
    assert aerobic.co2_exchange < 0


def test_detect_substrate_ignores_non_carbon_nutrient(ecoli_core):
    from cmm.features.production import _detect_substrate

    ecoli_core.reactions.EX_pi_e.lower_bound = -15  # tighter than glucose but no carbon
    assert _detect_substrate(ecoli_core) == GLC
