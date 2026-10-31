from __future__ import annotations

import time

import pytest
from cobra import Metabolite, Model, Reaction

# OptKnock/RobustKnock delegate the bilevel MILP to the optional ``straindesign`` package
# (the ``design`` extra). Skip this module gracefully where it is not installed so CI on a
# platform without it stays green instead of erroring at call time.
pytest.importorskip("straindesign")

from cmm.core.condition import Condition, ReactionBound  # noqa: E402
from cmm.core.solvers import SolverCapabilityError  # noqa: E402
from cmm.features.strain_design import (  # noqa: E402
    STRAIN_DESIGN_MAX_SECONDS,
    optknock,
    robustknock,
)

SUCC = "EX_succ_e"


@pytest.fixture
def anaerobic_ecoli(ecoli_core):
    ecoli_core.reactions.EX_o2_e.lower_bound = 0.0
    ecoli_core.reactions.EX_glc__D_e.lower_bound = -10.0
    return ecoli_core


@pytest.fixture
def restricted_license_strain_model():
    """Five-reaction coupling problem, deliberately far below bundled-solver limits."""

    model = Model("restricted_license_strain_design")
    substrate = Metabolite("substrate_c", compartment="c")
    precursor = Metabolite("precursor_c", compartment="c")
    product_metabolite = Metabolite("product_c", compartment="c")

    def reaction(reaction_id, stoichiometry, gene_rule="", upper_bound=1000.0):
        item = Reaction(reaction_id, lower_bound=0.0, upper_bound=upper_bound)
        item.add_metabolites(stoichiometry)
        item.gene_reaction_rule = gene_rule
        return item

    model.add_reactions(
        [
            reaction("SUPPLY", {substrate: 1.0}, upper_bound=10.0),
            reaction(
                "DIRECT",
                {substrate: -1.0, precursor: 2.0},
                "g_direct",
            ),
            reaction(
                "COUPLED",
                {substrate: -1.0, precursor: 1.0, product_metabolite: 1.0},
                "g_coupled",
            ),
            reaction("BIOMASS", {precursor: -1.0}, "g_biomass"),
            reaction("PRODUCT", {product_metabolite: -1.0}),
        ]
    )
    model.objective = "BIOMASS"
    return model


def test_restricted_license_strain_design_forwards_seed(
    restricted_license_strain_model,
):
    """OptKnock and RobustKnock solve a tiny MILP and preserve the requested seed."""

    for search in (optknock, robustknock):
        result = search(
            restricted_license_strain_model,
            "PRODUCT",
            max_knockouts=1,
            max_solutions=2,
            seed=17,
        )
        assert result.best() is not None
        assert result.best().knockouts == ("DIRECT",)
        assert result.best().guaranteed_product == pytest.approx(10.0)
        assert result.metadata["seed"] == 17
        assert result.metadata["parameters"]["strain_design_seed"] == 17


def test_optknock_finds_growth_coupled_succinate_design(anaerobic_ecoli):
    result = optknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=4)
    assert result.method == "optknock"
    assert len(result.designs) >= 1
    best = result.best()
    assert 1 <= len(best.knockouts) <= 3
    assert best.max_product > 5.0  # succinate is forced to a high flux
    assert best.growth > 0.0
    guaranteed = [design.guaranteed_product for design in result.designs]
    assert guaranteed == sorted(guaranteed, reverse=True)
    assert result.to_frame()["guaranteed_product"].tolist() == guaranteed


def test_robustknock_returns_only_guaranteed_designs(anaerobic_ecoli):
    result = robustknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=8)
    assert result.method == "robustknock"
    assert len(result.designs) >= 1
    # Every robust design must guarantee product at maximum growth (worst case > 0).
    assert all(d.guaranteed_product > 1e-6 for d in result.designs)
    assert all(d.growth_coupled for d in result.designs)
    # Ranked by guaranteed product (descending).
    guaranteed = [d.guaranteed_product for d in result.designs]
    assert guaranteed == sorted(guaranteed, reverse=True)

    exported = result.to_frame()
    assert exported["rank"].tolist() == list(range(1, len(result.designs) + 1))
    assert exported["guaranteed_product"].tolist() == guaranteed
    assert exported["growth_coupled"].all()
    assert exported.loc[0, "knockouts"] == ";".join(result.designs[0].knockouts)


def test_optknock_and_robustknock_use_distinct_nested_searches(
    anaerobic_ecoli, monkeypatch
):
    """Prevent RobustKnock from regressing to post-filtered OptKnock candidates."""

    import straindesign as sd

    seen: list[str] = []
    real_module = sd.SDModule

    def recording_module(model, module_type, *args, **kwargs):
        seen.append(module_type)
        return real_module(model, module_type, *args, **kwargs)

    monkeypatch.setattr(sd, "SDModule", recording_module)
    search_kwargs: list[dict[str, object]] = []

    def no_solutions(*args, **kwargs):
        del args
        search_kwargs.append(kwargs)
        return None

    monkeypatch.setattr(sd, "compute_strain_designs", no_solutions)

    opt_result = optknock(
        anaerobic_ecoli,
        SUCC,
        max_knockouts=1,
        max_solutions=1,
        seed=17,
    )
    robust_result = robustknock(
        anaerobic_ecoli,
        SUCC,
        max_knockouts=1,
        max_solutions=1,
        seed=23,
    )

    assert seen == [sd.OPTKNOCK, sd.ROBUSTKNOCK]
    assert [kwargs["compress"] for kwargs in search_kwargs] == [True, True]
    assert [kwargs["seed"] for kwargs in search_kwargs] == [17, 23]
    assert opt_result.metadata["seed"] == 17
    assert opt_result.metadata["parameters"]["strain_design_seed"] == 17
    assert robust_result.metadata["seed"] == 23
    assert robust_result.metadata["parameters"]["strain_design_seed"] == 23


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"max_knockouts": 0}, "max_knockouts"),
        ({"max_solutions": 0}, "max_solutions"),
        ({"min_growth": -0.1}, "min_growth"),
        ({"seed": -1}, "seed"),
        ({"seed": 2_000_000_001}, "seed"),
        ({"seed": 1.5}, "seed"),
    ],
)
def test_strain_design_validates_search_parameters(anaerobic_ecoli, kwargs, message):
    with pytest.raises(ValueError, match=message):
        optknock(anaerobic_ecoli, SUCC, **kwargs)


def test_strain_design_does_not_mutate_model(anaerobic_ecoli):
    growth = anaerobic_ecoli.slim_optimize()
    optknock(anaerobic_ecoli, SUCC, max_knockouts=2, max_solutions=2)
    assert anaerobic_ecoli.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_strain_design_handles_uncouplable_product(anaerobic_ecoli):
    # A product that cannot be produced in the medium yields no designs, not a crash.
    anaerobic_ecoli.reactions.EX_succ_e.bounds = (0.0, 0.0)
    result = optknock(anaerobic_ecoli, SUCC, max_knockouts=2, max_solutions=2)
    assert result.designs == ()
    assert result.best() is None


def test_strain_design_requires_milp(anaerobic_ecoli):
    anaerobic_ecoli.solver = "glpk_exact"  # LP only, no MILP
    with pytest.raises(SolverCapabilityError) as exc:
        optknock(anaerobic_ecoli, SUCC, max_knockouts=2)
    assert exc.value.capability == "MILP"


def test_optknock_candidate_filter_excludes_unrealisable_knockouts(anaerobic_ecoli):
    """Burgard restrict candidates to central metabolism; CMM did not.

    Without the filter the search returns designs deleting boundary reactions with no GPR
    (``EX_co2_e``, ``EX_ac_e``, ``EX_for_e``, ``EX_lac__D_e``), which cannot be realised as
    gene deletions. ``_actionable_reaction`` already existed in ``production.py`` and was
    simply not applied here.
    """

    result = optknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=5)
    assert result.designs
    for design in result.designs:
        for rid in design.knockouts:
            reaction = anaerobic_ecoli.reactions.get_by_id(rid)
            assert not reaction.boundary, f"{rid} is an exchange reaction"
            assert reaction.genes, f"{rid} has no GPR"
    assert result.metadata["parameters"]["actionable_only"] is True
    assert result.metadata["parameters"]["n_knockout_candidates"] < len(
        anaerobic_ecoli.reactions
    )


def test_optknock_candidate_filter_can_be_opted_out(anaerobic_ecoli):
    unfiltered = optknock(
        anaerobic_ecoli,
        SUCC,
        max_knockouts=3,
        max_solutions=5,
        actionable_only=False,
    )
    filtered = optknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=5)
    assert unfiltered.metadata["parameters"]["n_knockout_candidates"] == len(
        anaerobic_ecoli.reactions
    )
    # The opt-out is what restores the exchange knockouts, so it must return strictly more.
    assert len(unfiltered.designs) > len(filtered.designs)
    assert any(
        anaerobic_ecoli.reactions.get_by_id(rid).boundary
        for design in unfiltered.designs
        for rid in design.knockouts
    )


def test_optknock_designs_are_deduplicated_by_knockout_set(anaerobic_ecoli):
    """``max_solutions`` caps MILP solutions, not designs - so designs must be collapsed."""

    result = optknock(anaerobic_ecoli, SUCC, max_knockouts=3, max_solutions=5)
    knockout_sets = [frozenset(d.knockouts) for d in result.designs]
    assert len(knockout_sets) == len(set(knockout_sets))
    parameters = result.metadata["parameters"]
    assert parameters["n_designs_after_deduplication"] == len(result.designs)
    assert parameters["n_milp_designs"] >= parameters["n_designs_after_deduplication"]
    assert parameters["strain_design_backend"] == "straindesign"
    assert parameters["straindesign_version"]
    assert parameters["straindesign_search_status"] in {"optimal", "time_limit"}
    assert parameters["straindesign_search_complete"] is (
        parameters["straindesign_search_status"] == "optimal"
    )
    assert parameters["straindesign_compress"] is True
    assert parameters["max_compressed_milp_solutions_requested"] == 5
    assert "compressed MILP solutions" in parameters["max_solutions_semantics"]
    assert parameters["n_compressed_milp_solutions"] is None
    assert parameters["n_decompressed_designs"] >= len(result.designs)
    assert parameters["n_unique_decompressed_designs"] >= len(result.designs)
    assert parameters["n_returned_designs_after_deduplication"] == len(result.designs)
    assert parameters["n_milp_designs_semantics"].startswith("legacy alias")
    assert parameters["intervention_level"] == "reaction"
    assert parameters["requires_gpr_resolution"] is True


def test_strain_design_accepts_and_records_a_condition(ecoli_core):
    """optknock/robustknock took no condition and depended on the caller's model state."""

    condition = Condition(
        name="anaerobic",
        bounds=(ReactionBound("EX_o2_e", lower_bound=0.0),),
    )
    before = ecoli_core.reactions.EX_o2_e.lower_bound
    result = optknock(
        ecoli_core, SUCC, max_knockouts=3, max_solutions=3, condition=condition
    )
    applied = result.metadata["parameters"]["applied_condition"]
    assert applied["condition"] == "anaerobic"
    assert applied["aerobic"] is False
    assert applied["oxygen_exchange"]["lower_bound"] == pytest.approx(0.0)
    assert ecoli_core.reactions.EX_o2_e.lower_bound == before  # model left untouched
    assert result.designs


class _Solutions:
    """The part of straindesign's ``SDSolutions`` that CMM reads."""

    def __init__(self, status: str, reaction_sd: list | None = None) -> None:
        self.status = status
        self.reaction_sd = reaction_sd or []


def test_strain_design_forwards_a_time_limit_to_straindesign(
    anaerobic_ecoli, monkeypatch
):
    """straindesign has no time limit unless it is given one, and then it can never return."""

    import straindesign as sd

    assert STRAIN_DESIGN_MAX_SECONDS == 3600.0
    seen: list[object] = []

    def record(*args, **kwargs):
        del args
        seen.append(kwargs.get("time_limit", "absent"))
        return _Solutions("optimal")

    monkeypatch.setattr(sd, "compute_strain_designs", record)

    default = optknock(anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1)
    explicit = robustknock(
        anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1, max_seconds=120
    )
    unbounded = optknock(
        anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1, max_seconds=None
    )

    assert seen == [3600.0, 120.0, "absent"]
    assert default.metadata["parameters"]["max_seconds"] == 3600.0
    assert explicit.metadata["parameters"]["max_seconds"] == 120
    assert unbounded.metadata["parameters"]["max_seconds"] is None


@pytest.mark.parametrize(
    "value", [0, -1.0, float("nan"), float("inf"), True, "60"], ids=repr
)
def test_strain_design_rejects_a_meaningless_time_limit(anaerobic_ecoli, value):
    with pytest.raises(ValueError, match="max_seconds"):
        optknock(anaerobic_ecoli, SUCC, max_seconds=value)


def test_a_search_cut_by_the_clock_is_not_reported_complete(
    anaerobic_ecoli, monkeypatch
):
    """straindesign leaves status 'optimal' when the clock runs out between two solves.

    Its outer loop tests the remaining time only before a solve, and ``status`` is the last
    solve's. So the limit can end a search that says 'optimal' with fewer designs than a
    finished search would have - which CMM then reported as a complete search.
    """

    import straindesign as sd

    def runs_past_its_limit(*args, **kwargs):
        del args
        time.sleep(kwargs["time_limit"] + 0.05)
        return _Solutions("optimal")

    monkeypatch.setattr(sd, "compute_strain_designs", runs_past_its_limit)

    with pytest.warns(UserWarning, match="time limit"):
        result = optknock(
            anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1, max_seconds=0.05
        )

    parameters = result.metadata["parameters"]
    assert parameters["straindesign_search_status"] == "optimal"  # backend's own word
    assert parameters["straindesign_time_limit_reached"] is True
    assert parameters["straindesign_search_complete"] is False
    assert parameters["straindesign_search_elapsed_s"] >= 0.05


@pytest.mark.parametrize("status", ["time_limit", "time_limit_w_sols"])
def test_a_backend_reported_time_limit_is_never_complete(
    anaerobic_ecoli, monkeypatch, status
):
    import straindesign as sd

    monkeypatch.setattr(
        sd, "compute_strain_designs", lambda *a, **k: _Solutions(status)
    )
    with pytest.warns(UserWarning, match="time limit"):
        result = robustknock(
            anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1, max_seconds=600
        )
    parameters = result.metadata["parameters"]
    assert parameters["straindesign_search_status"] == status
    assert parameters["straindesign_time_limit_reached"] is True
    assert parameters["straindesign_search_complete"] is False


def test_a_search_that_finishes_inside_its_limit_stays_complete(
    anaerobic_ecoli, monkeypatch, recwarn
):
    import straindesign as sd

    monkeypatch.setattr(
        sd, "compute_strain_designs", lambda *a, **k: _Solutions("optimal")
    )
    unbounded = optknock(
        anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1, max_seconds=None
    )
    bounded = optknock(
        anaerobic_ecoli, SUCC, max_knockouts=1, max_solutions=1, max_seconds=600
    )
    for result in (unbounded, bounded):
        parameters = result.metadata["parameters"]
        assert parameters["straindesign_time_limit_reached"] is False
        assert parameters["straindesign_search_complete"] is True
    assert not [w for w in recwarn if "time limit" in str(w.message)]


def test_a_stalled_straindesign_search_returns_at_the_time_limit(
    restricted_license_strain_model, monkeypatch
):
    """Regression: straindesign 1.18 ``compute_optimal`` re-solved the same MILP forever.

    After the MILP finds its optimum it runs a second solve to minimise the number of
    interventions. When that solve is not optimal (Gurobi reported 'Model is infeasible'
    straight after 'max constraint violation 4.75e-09 exceeds tolerance' on a 1,440-reaction
    model) the inner loop breaks WITHOUT adding an exclusion constraint, the outer loop
    solves the identical MILP again, and ``max_solutions`` is never reached because nothing
    is ever added. A production run sat in it for 24 hours with no output.

    The second solve of every outer iteration is forced non-optimal here; everything else
    is straindesign's real code path.
    """

    from straindesign.names import INFEASIBLE
    from straindesign.strainDesignMILP import SDMILP

    real_solve = SDMILP.solveZ
    solves = {"n": 0}

    def minimisation_never_optimal(self):
        solves["n"] += 1
        # Without the time limit this loop does not end. Fail instead of hanging CI.
        assert solves["n"] < 200_000, "straindesign kept re-solving past the time limit"
        z, x, opt, status = real_solve(self)
        if solves["n"] % 2 == 0:  # the "minimise interventions" solve
            status = INFEASIBLE
        return z, x, opt, status

    monkeypatch.setattr(SDMILP, "solveZ", minimisation_never_optimal)

    started = time.monotonic()
    with pytest.warns(UserWarning, match="time limit"):
        result = optknock(
            restricted_license_strain_model,
            "PRODUCT",
            max_knockouts=1,
            max_solutions=2,
            seed=17,
            max_seconds=2.0,
        )

    assert time.monotonic() - started < 60
    parameters = result.metadata["parameters"]
    assert parameters["straindesign_time_limit_reached"] is True
    assert parameters["straindesign_search_complete"] is False
    assert parameters["max_seconds"] == 2.0
    assert not result.designs  # nothing was recorded: the search never made progress
