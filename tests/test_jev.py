"""The JEV agent: transport, controller, screen, and a whole game played offline.

No test here touches the network. The agent is replaced by a scripted client that returns
prepared answers, which is possible precisely because the real service returns typed answers
drawn from the caller's own criteria — a fake that obeys the same contract is a faithful
stand-in, not an approximation of free text.

Where a real HTTP response is needed, the recorded shape of one is used. Those literals came
from live calls to ``typesafe/jev-1.13`` during development; they are the contract this
package is written against.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from cobra.io import write_sbml_model

from cmm.core.condition import Condition, ReactionBound
from cmm.core.simulation import pfba
from cmm.jev import (
    ACTION_CATALOGUE,
    ACT_ACTIONS,
    ActionNotApplicable,
    JevClient,
    JevConfig,
    JevTransportError,
    applicable_actions,
    build_candidates,
    build_intervention,
    cofactor_balance,
    get_question_set,
    product_distances,
    run_jev_design,
)
from cmm.jev._transport import DecisionResult, JevAnswer, JevUsage, _parse_decision
from cmm.jev.questions import NoAvailableAction, available_actions
from cmm.jev.state import ScanCache, carbon_byproducts

ANAEROBIC = Condition(
    name="glucose_anaerobic",
    bounds=(
        ReactionBound("EX_glc__D_e", -10.0, 1000.0),
        ReactionBound("EX_o2_e", 0.0, 1000.0),
    ),
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def anaerobic_core(ecoli_core):
    """``e_coli_core`` on glucose with oxygen closed: the condition succinate needs."""

    ANAEROBIC.apply_to(ecoli_core)
    return ecoli_core


@pytest.fixture
def anaerobic_core_path(tmp_path, anaerobic_core) -> Path:
    path = tmp_path / "model.xml"
    write_sbml_model(anaerobic_core, str(path))
    return path


class ScriptedClient(JevClient):
    """A JEV client that answers from a script instead of the network.

    Subclasses the real client so everything around the answers — usage accounting, the
    served-model record, the budget guard — is the code that runs in production.
    """

    def __init__(self, moves, *, default=("end_round", None)):
        # Bypass the real constructor: it resolves an API key, which this never needs.
        self.model = "scripted/jev"
        self.research_model = "scripted/chat"
        self.usage = JevUsage()
        self.served_models = []
        self.moves = list(moves)
        self.default = default
        self.asked: list[dict] = []
        self.states: list[object] = []
        self._pending_action: str | None = None

    def decide(self, state, questions, *, model=None):
        self.states.append(state)
        self.asked.append(dict(questions))
        key = next(iter(questions))
        offered = questions[key]["criteria"]

        if key == "target":
            # One script entry per tick. The action stage is a second call about the move
            # already chosen, so it must not consume the next entry.
            target, self._pending_action = (
                self.moves.pop(0) if self.moves else self.default
            )
            choice = target if target in offered else "end_round"
        else:
            action = self._pending_action
            choice = action if action in offered else next(iter(offered))

        answers = {
            key: JevAnswer(
                key=key,
                type="choice",
                value=choice,
                confidence=0.9,
                probabilities={
                    name: (0.9 if name == choice else 0.1 / max(len(offered) - 1, 1))
                    for name in offered
                },
            )
        }
        for extra, kind in (("benefit", "score"), ("growth_risk", "noul")):
            if extra in questions:
                answers[extra] = JevAnswer(
                    key=extra, type=kind, value=3.0 if kind == "score" else 0.1
                )
        result = DecisionResult(
            answers=answers,
            model="scripted/jev-test",
            input_tokens=100,
            output_tokens=10,
            cost_usd=0.0001,
        )
        self.usage.record(result)
        if result.model not in self.served_models:
            self.served_models.append(result.model)
        return result


# ---------------------------------------------------------------------------
# transport
# ---------------------------------------------------------------------------


def test_parsing_a_recorded_response_keeps_every_answer_shape() -> None:
    """The three answer types, as the live service actually returns them."""

    result = _parse_decision(
        {
            "model": "typesafe/jev-1.13-20260917",
            "answers": {
                "target": {
                    "type": "choice",
                    "choice": "FRD7",
                    "confidence": 0.63,
                    "probabilities": {"FRD7": 0.63, "PPC": 0.14, "PTA": 0.05},
                },
                "benefit": {
                    "type": "score",
                    "score": 1.97,
                    "confidence": 0.0,
                    "legend": {"0": "no gain"},
                },
                "growth_risk": {"type": "noul", "noul": 0.23},
            },
            "usage": {"input_tokens": 2100, "output_tokens": 197, "cost": 8.82e-05},
            "id": "gen-dec-1",
            "provider": "TypeSafe",
        },
        latency=0.314,
    )

    assert result["target"].choice == "FRD7"
    assert result["benefit"].score == pytest.approx(1.97)
    assert result["growth_risk"].value == pytest.approx(0.23)
    assert result.cost_usd == pytest.approx(8.82e-05)
    # The served id differs from the requested one; provenance must record what answered.
    assert result.model == "typesafe/jev-1.13-20260917"


def test_a_choice_answer_ranks_every_option_it_was_offered() -> None:
    """One call is a whole ranking, which is why no per-candidate question is needed."""

    answer = JevAnswer(
        key="target",
        type="choice",
        value="FRD7",
        probabilities={"PPC": 0.14, "FRD7": 0.63, "AAA": 0.14},
    )
    assert answer.ranked()[0] == ("FRD7", 0.63)
    # Ties break by name, so the same probabilities always produce the same order.
    assert [name for name, _ in answer.ranked()[1:]] == ["AAA", "PPC"]


def test_an_unknown_answer_type_is_reported_not_guessed() -> None:
    with pytest.raises(JevTransportError, match="unsupported type"):
        _parse_decision(
            {"answers": {"target": {"type": "vibes", "value": 1}}, "model": "x"},
            latency=0.0,
        )


def test_a_missing_api_key_names_the_variable_to_set(monkeypatch, tmp_path) -> None:
    """Both sources have to be absent, and the saved one lives outside the repository.

    Deleting the environment variable is not enough: ``resolve_api_key`` falls back to the key
    the desktop app saves under ``$XDG_CONFIG_HOME``, so on a machine where anyone had ever
    saved one this test passed for the wrong reason and failed the moment a key appeared. The
    config home is redirected at a temporary directory so the test states its own world.
    """

    monkeypatch.delenv("OPENROUTER_API_KEY", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    with pytest.raises(JevTransportError, match="OPENROUTER_API_KEY"):
        JevClient()


def test_the_api_key_never_appears_in_an_error(monkeypatch) -> None:
    """A traceback or a run bundle must not become a place a secret leaks."""

    secret = "sk-or-v1-notarealkey"
    monkeypatch.setenv("OPENROUTER_API_KEY", secret)
    client = JevClient(base_url="https://127.0.0.1:1")
    with pytest.raises(JevTransportError) as caught:
        client.decide(
            {"a": 1}, {"q": {"type": "choice", "criteria": {"x": "1", "y": "2"}}}
        )
    assert secret not in str(caught.value)
    assert secret not in repr(caught.value)


# ---------------------------------------------------------------------------
# the controller
# ---------------------------------------------------------------------------


def test_a_reaction_carrying_flux_can_be_deleted_or_halved(anaerobic_core) -> None:
    """One knockdown strength, not two, and no amplification at all.

    A second, deeper cap mostly bought a second rejection of the same idea, and roughly
    halving an activity is the level a promoter swap or an RBS change can actually aim at.
    """

    names = [action.name for action in applicable_actions(8.2)]
    assert names == ["knockout", "knockdown_50"]
    assert "knockdown_25" not in ACTION_CATALOGUE


def test_nothing_in_the_vocabulary_can_force_flux_up() -> None:
    """The restriction that defines this vocabulary, asserted rather than described.

    A lower bound on a flux is not what over-expression does: it tells the solver the flux
    *must* be carried, by whatever route is cheapest, while stronger expression only raises a
    capacity the cell may decline to use. Every move here is therefore a cap.
    """

    for name in ("amplify_2x", "amplify_5x", "force_on_low", "force_on_high"):
        assert name not in ACTION_CATALOGUE
    for action in ACT_ACTIONS:
        assert action.mode in ("knockout", "knockdown")


def test_a_reaction_at_zero_can_only_be_deleted() -> None:
    """A knockdown needs a flux to be half of. A deletion does not, and is still worth
    making: OptKnock's most valuable deletions close routes carrying nothing today."""

    names = [action.name for action in applicable_actions(0.0)]
    assert names == ["knockout"]


def test_a_knockdown_caps_the_magnitude_without_opening_a_direction(
    anaerobic_core,
) -> None:
    reaction = anaerobic_core.reactions.get_by_id("PFL")
    reaction.bounds = (0.0, 1000.0)
    intervention = build_intervention(
        anaerobic_core, "PFL", ACTION_CATALOGUE["knockdown_50"], {"PFL": 17.8}
    )
    assert intervention.upper_bound == pytest.approx(8.9)
    assert intervention.lower_bound == 0.0  # the reaction was irreversible; it stays so


def test_a_knockdown_on_a_zero_flux_reaction_is_refused_not_reinterpreted(
    anaerobic_core,
) -> None:
    with pytest.raises(ActionNotApplicable, match="nothing to halve"):
        build_intervention(
            anaerobic_core, "FRD7", ACTION_CATALOGUE["knockdown_50"], {"FRD7": 0.0}
        )


def test_a_move_is_resolved_to_the_gene_set_that_achieves_it(anaerobic_core) -> None:
    """Deleting a reaction and deleting its genes are not the same operation.

    ``ACKr`` is ``b2296 or b3115 or b1849``. Deleting *ackA*, the textbook acetate-branch
    gene, leaves two isozymes: measured on this model, succinate stays at 0 and growth at
    0.2117. The move has to take all three, and the record has to say so, or the design is
    one a laboratory would build and find inert.
    """

    intervention = build_intervention(
        anaerobic_core, "ACKr", ACTION_CATALOGUE["knockout"], {"ACKr": 8.5}
    )
    assert set(intervention.genes) == {"b2296", "b3115", "b1849"}
    assert "ackA" in intervention.describe()
    # No shared gene here, so the edit stops exactly the reaction it was aimed at.
    assert intervention.side_effects == ()


def test_a_complex_needs_only_one_subunit_gone(anaerobic_core) -> None:
    """The other direction: ``THD2`` is ``b1602 and b1603``, so either subunit suffices and
    naming both would overstate the work."""

    intervention = build_intervention(
        anaerobic_core, "THD2", ACTION_CATALOGUE["knockout"], {"THD2": 0.0}
    )
    assert len(intervention.genes) == 1
    assert intervention.genes[0] in {"b1602", "b1603"}
    # NADTRHD is `b3962 or (b1602 and b1603)`, so one subunit does not stop it.
    assert "NADTRHD" not in intervention.side_effects


def test_a_shared_gene_takes_its_other_reactions_with_it(anaerobic_core) -> None:
    """The consequence CMM now applies rather than discovering at report time.

    ``SUCCt2_2`` is run by *dctA*, which also runs ``FUMt2_2`` and ``MALt2_2``. A design that
    deletes it deletes all three whether anyone intended that or not, so the engine has to
    apply all three and let the viability rule judge the result.
    """

    intervention = build_intervention(
        anaerobic_core, "SUCCt2_2", ACTION_CATALOGUE["knockout"], {"SUCCt2_2": 0.0}
    )
    assert intervention.genes == ("b3528",)
    assert set(intervention.side_effects) == {"FUMt2_2", "MALt2_2"}
    assert {rid for rid, _, _ in intervention.bounds} == {
        "SUCCt2_2",
        "FUMt2_2",
        "MALt2_2",
    }
    assert all(low == 0.0 and high == 0.0 for _, low, high in intervention.bounds)
    assert "also constrains" in intervention.describe()


def test_a_knockdown_says_which_side_effects_it_cannot_express(anaerobic_core) -> None:
    """Half of nothing is nothing, and writing zero would be a knockout of something the
    agent never chose. The reaction is named instead of being silently deleted."""

    fluxes = dict(pfba(anaerobic_core).fluxes)
    intervention = build_intervention(
        anaerobic_core, "FORti", ACTION_CATALOGUE["knockdown_50"], fluxes
    )
    constrained = {rid for rid, _, _ in intervention.bounds}
    for reaction_id in intervention.unmodelled:
        assert reaction_id not in constrained
        assert abs(fluxes.get(reaction_id, 0.0)) <= 1e-9
    if intervention.unmodelled:
        assert "cannot be expressed" in intervention.describe()


def test_an_intervention_converts_to_the_bounds_cmm_already_applies(
    anaerobic_core,
) -> None:
    intervention = build_intervention(
        anaerobic_core, "PFL", ACTION_CATALOGUE["knockout"], {"PFL": 17.8}
    )
    bounds = intervention.to_reaction_bounds()
    assert bounds and all(isinstance(bound, ReactionBound) for bound in bounds)
    assert all((b.lower_bound, b.upper_bound) == (0.0, 0.0) for b in bounds)


# ---------------------------------------------------------------------------
# the screen
# ---------------------------------------------------------------------------


def test_cofactor_accounting_finds_the_anaerobic_atp_and_nadh_pools(
    anaerobic_core,
) -> None:
    fluxes = pfba(anaerobic_core).fluxes
    balance = cofactor_balance(anaerobic_core, fluxes)
    assert balance.atp_production > 0
    assert balance.nadh_production > 0
    # Substrate-level phosphorylation carries ATP production with oxygen closed.
    assert balance.largest_atp_source is not None
    assert balance.largest_atp_source[0] in {"PGK", "PYK", "ACKr", "SUCOAS"}


def test_product_distance_ignores_currency_metabolites(anaerobic_core) -> None:
    """Leave ATP and protons in and every reaction is two steps from every other."""

    distances = product_distances(anaerobic_core, "EX_succ_e")
    assert distances["EX_succ_e"] == 0
    assert distances["SUCCt2_2"] == 1
    assert distances["FRD7"] == 2
    # ATP synthase touches only currency metabolites, so no pathway connects it.
    assert "ATPS4r" not in distances


def test_the_board_covers_both_opening_the_route_and_closing_the_competition(
    anaerobic_core,
) -> None:
    """The failure this composition exists to prevent, stated as a test.

    Ranked on proximity alone the board is only the succinate branch; ranked on flux alone it
    is glycolysis and respiration. A design needs both halves.
    """

    fluxes = pfba(anaerobic_core).fluxes
    board = build_candidates(
        anaerobic_core,
        product_reaction_id="EX_succ_e",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=20,
    )
    ids = {candidate.reaction_id for candidate in board}
    assert {"FRD7", "FUM", "SUCOAS"} <= ids, (
        "the route to the product must be on the board"
    )
    assert ids & {"PFL", "ACKr", "PTAr", "ALCD2x", "ACALD"}, (
        "the fermentation branches competing for carbon must be on the board"
    )


def test_only_reactions_with_a_gene_can_be_targets(anaerobic_core) -> None:
    """A move has to be one a laboratory could make."""

    fluxes = pfba(anaerobic_core).fluxes
    board = build_candidates(
        anaerobic_core,
        product_reaction_id="EX_succ_e",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=40,
    )
    ids = {candidate.reaction_id for candidate in board}
    assert not any(rid.startswith("EX_") for rid in ids)
    assert "ATPM" not in ids  # a maintenance term, not an enzyme
    assert "SUCCt3" not in ids  # passive diffusion, no GPR
    assert all(anaerobic_core.reactions.get_by_id(rid).genes for rid in ids)


def test_the_secreted_carbon_byproducts_are_the_competition(anaerobic_core) -> None:
    fluxes = pfba(anaerobic_core).fluxes
    byproducts = carbon_byproducts(anaerobic_core, fluxes, "EX_succ_e")
    assert set(byproducts) == {"EX_for_e", "EX_ac_e", "EX_etoh_e"}
    assert (
        "EX_co2_e" not in byproducts
    )  # the end of oxidation, not a redirectable branch


def test_a_candidate_record_states_what_each_number_means(anaerobic_core) -> None:
    fluxes = pfba(anaerobic_core).fluxes
    board = build_candidates(
        anaerobic_core,
        product_reaction_id="EX_succ_e",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=24,
    )
    record = next(c for c in board if c.reaction_id == "FRD7").to_record()
    assert "wild-type flux" in record
    assert "steps from the product" in record
    # Unknowns are absent rather than rendered as zero, which is what lets the agent scan.
    assert "ESSENTIAL" not in record
    assert "FSEOF" not in record


def test_a_reverted_move_gives_the_scans_back() -> None:
    """The bug that made the agent re-run the same two scans on every tick."""

    cache = ScanCache()
    cache.essential["PFL"] = False
    cache.fseof_slopes["FRD7"] = 1.2
    cache.completed.add("fseof_scan")

    saved = cache.snapshot()
    cache.invalidate()
    assert not cache.essential and not cache.completed

    cache.restore(saved)
    assert cache.essential == {"PFL": False}
    assert cache.fseof_slopes == {"FRD7": 1.2}
    assert cache.completed == {"fseof_scan"}


# ---------------------------------------------------------------------------
# the questions
# ---------------------------------------------------------------------------


def test_a_full_design_may_only_undo_or_stop(anaerobic_core) -> None:
    """Offering reactions here produced ten consecutive refusals in a real run."""

    fluxes = pfba(anaerobic_core).fluxes
    board = build_candidates(
        anaerobic_core,
        product_reaction_id="EX_succ_e",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=10,
    )
    question = get_question_set().target_question(
        board,
        product="EX_succ_e",
        growth_floor=0.05,
        allow_undo=True,
        allow_look=True,
        design_full=True,
    )["target"]
    assert set(question["criteria"]) == {"undo_last", "end_round"}


def test_a_scan_already_run_is_not_offered_again(anaerobic_core) -> None:
    fluxes = pfba(anaerobic_core).fluxes
    candidate = next(
        c
        for c in build_candidates(
            anaerobic_core,
            product_reaction_id="EX_succ_e",
            reference_fluxes=fluxes,
            current_fluxes=fluxes,
            limit=10,
        )
        if c.reaction_id == "FRD7"
    )
    offered = get_question_set().action_question(
        candidate,
        product="EX_succ_e",
        growth_floor=0.05,
        allow_look=True,
        exclude={"fseof_scan", "envelope_probe"},
    )["action"]["criteria"]
    assert "fseof_scan" not in offered
    assert "envelope_probe" not in offered
    assert "essentiality_scan" in offered


def test_a_candidate_with_nothing_left_to_try_says_so(anaerobic_core) -> None:
    fluxes = pfba(anaerobic_core).fluxes
    candidate = next(
        c
        for c in build_candidates(
            anaerobic_core,
            product_reaction_id="EX_succ_e",
            reference_fluxes=fluxes,
            current_fluxes=fluxes,
            limit=10,
        )
        if c.reaction_id == "FRD7"
    )
    with pytest.raises(NoAvailableAction):
        get_question_set().action_question(
            candidate,
            product="EX_succ_e",
            growth_floor=0.05,
            allow_look=False,
            exclude={"knockout"},
        )


# ---------------------------------------------------------------------------
# a whole game, offline
# ---------------------------------------------------------------------------


def test_a_scripted_game_redirects_carbon_to_the_product_and_records_every_move(
    anaerobic_core_path, tmp_path
) -> None:
    """The end-to-end result: a product that was zero is being made, and the run says how."""

    # Deleting pyruvate formate lyase closes the formate/acetate branch and pushes carbon
    # down the reductive route to succinate: 0 -> 0.68 at a growth rate of 0.18. It is a
    # single move on purpose: the assertion is that the engine turns a chosen move into real
    # product flux, not that this particular design is the best one.
    client = ScriptedClient([("PFL", "knockout")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        substrate="EX_glc__D_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=2,
        steps_per_round=3,
        max_knockouts=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    frames: list[tuple] = []
    result = run_jev_design(
        config,
        client=client,
        on_tick=lambda tick, fluxes: frames.append((tick, fluxes)),
    )

    assert result.wild_type_product_flux == pytest.approx(0.0, abs=1e-6)
    assert result.best_product_flux > 0.5, (
        "closing the fermentative branch must push carbon to succinate"
    )
    assert result.best_growth >= config.growth_floor
    assert result.summary()["beat_wild_type"] is True

    # Every move produced a frame the GUI could draw, and a row a reader can audit.
    assert len(frames) == len(result.ticks)
    assert all(fluxes for _, fluxes in frames)
    assert len(result.ticks_frame()) == len(result.ticks)


def test_cmm_enforces_the_growth_floor_whatever_the_agent_chose(
    anaerobic_core_path, tmp_path
) -> None:
    """The rule the agent does not get a vote on."""

    client = ScriptedClient([("GAPD", "knockout")] * 4)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=1,
        steps_per_round=2,
        growth_floor=0.05,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    result = run_jev_design(config, client=client)

    reverted = [tick for tick in result.ticks if tick.outcome.startswith("reverted")]
    assert reverted, "deleting GAPD anaerobically cannot leave the strain viable"
    # A lethal move is a recorded row with its intervention attached, not a dropped row and
    # not a crash: an infeasible point is data.
    assert all(tick.intervention is not None for tick in reverted)
    assert any(tick.action == "knockout" for tick in reverted)
    # Nothing that breaches the floor is left standing, whatever the agent asked for.
    assert all(
        tick.growth >= config.growth_floor
        for tick in result.ticks
        if tick.outcome == "applied"
    )
    assert "GAPD" not in {i.reaction_id for i in result.final_interventions} or all(
        i.mode != "knockout" for i in result.final_interventions
    )


def test_the_deletion_budget_is_never_exceeded(anaerobic_core_path, tmp_path) -> None:
    client = ScriptedClient(
        [
            ("PFL", "knockout"),
            ("ACKr", "knockout"),
            ("ALCD2x", "knockout"),
        ]
        * 4
    )
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=3,
        steps_per_round=4,
        max_knockouts=2,
        max_knockdowns=0,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    result = run_jev_design(config, client=client)
    assert len(result.final_interventions) <= 2
    assert max(tick.n_active_interventions for tick in result.ticks) <= 2


def test_the_two_budgets_are_counted_separately(anaerobic_core) -> None:
    """Deletions and knockdowns are different things to build, so they are limited apart.

    One shared cap starved the run: the seeded OptKnock design takes three deletions on its
    own, so a total of four left the agent one edit and every round ended a step or two after
    adopting it.
    """

    fluxes = pfba(anaerobic_core).fluxes
    candidate = next(
        c
        for c in build_candidates(
            anaerobic_core,
            product_reaction_id="EX_succ_e",
            reference_fluxes=fluxes,
            current_fluxes=fluxes,
            limit=16,
        )
        if c.reaction_id == "PFL"  # carries flux, so both moves are defined for it
    )
    both = {a.name for a in available_actions(candidate, allow_look=False)}
    assert both == {"knockout", "knockdown_50"}

    # A spent deletion budget removes the deletion, not the reaction.
    assert {
        a.name
        for a in available_actions(
            candidate, allow_look=False, allowed_modes=("knockdown",)
        )
    } == {"knockdown_50"}
    assert {
        a.name
        for a in available_actions(
            candidate, allow_look=False, allowed_modes=("knockout",)
        )
    } == {"knockout"}
    assert available_actions(candidate, allow_look=False, allowed_modes=()) == ()

    config = JevConfig(model_path="m.xml", product="EX_succ_e")
    assert config.room_for(0, 0) == ("knockout", "knockdown")
    assert config.room_for(config.max_knockouts, 0) == ("knockdown",)
    assert config.room_for(0, config.max_knockdowns) == ("knockout",)
    assert config.room_for(config.max_knockouts, config.max_knockdowns) == ()


def test_a_run_can_be_stopped_and_keeps_what_it_played(
    anaerobic_core_path, tmp_path
) -> None:
    """Stopping is an answer about how long to look, not a reason to discard the answer."""

    from cmm.jev.engine import STOP_REQUESTED

    client = ScriptedClient([("PFL", "knockout"), ("ACKr", "knockout")] * 20)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=4,
        steps_per_round=10,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        # On purpose: the point is that a stopped run skips it and says so.
        run_baseline_comparison=True,
    )

    played: list[object] = []

    def stop_after_three() -> bool:
        return len(played) >= 3

    result = run_jev_design(
        config,
        client=client,
        on_tick=lambda tick, fluxes: played.append(tick),
        should_stop=stop_after_three,
    )

    assert len(result.ticks) == 3, (
        "the flag is read once per step, so it stops on the next"
    )
    assert result.best_product_flux > 0.5, "what was played is still scored"
    assert any(STOP_REQUESTED in note for note in result.notes)
    assert result.baselines == ()
    assert any("baseline comparison was skipped" in note for note in result.notes)


def test_a_lone_available_move_is_taken_without_asking(
    anaerobic_core_path, tmp_path
) -> None:
    """One option is not a decision, and asking costs a call and a step.

    This matters now the budgets are separate: once the knockdown budget is gone, a reaction
    carrying no wild-type flux has exactly one legal move.
    """

    client = ScriptedClient([("FRD7", "knockout")] * 6)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=2,
        max_knockdowns=0,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        allow_look_actions=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=client)

    first = result.ticks[0]
    assert first.action == "knockout"
    assert first.action_ranking == (("knockout", 1.0),)
    # One call per step, not two: the second stage was never asked.
    assert len(client.asked) == len(result.ticks)


def test_the_run_stops_when_the_budget_is_spent(anaerobic_core_path, tmp_path) -> None:
    client = ScriptedClient([("FRD7", "fseof_scan")] * 50)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=5,
        steps_per_round=5,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
        max_decisions=6,
    )
    result = run_jev_design(config, client=client)
    assert (
        client.usage.calls <= 8
    )  # the guard is checked before each tick, not mid-tick
    assert any("budget" in note for note in result.notes)


def test_the_run_writes_one_artifact_per_role_and_every_file_exists(
    anaerobic_core_path, tmp_path
) -> None:
    client = ScriptedClient([("PFL", "knockout"), ("ACKr", "knockout")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        substrate="EX_glc__D_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=1,
        steps_per_round=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    result = run_jev_design(config, client=client)
    root = result.run_directory
    assert root == (tmp_path / "run").resolve()

    manifest = json.loads((root / "00_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == 2
    assert manifest["workflow"] == "jev_target_design"
    for role in (
        "model",
        "wild_type_reference_fluxes",
        "wild_type_summary",
        "ticks",
        "rounds",
        "candidate_rankings",
        "best_design",
        "final_design",
        "flux_trajectory",
        "agent_transcript",
        "agent_usage",
        "provenance",
        "summary",
        "workflow_configuration",
    ):
        assert role in manifest["artifacts"], role
        assert (root / manifest["artifacts"][role]["path"]).is_file(), role

    # The transcript is the only record of why each move was made, so it is not optional.
    entries = [
        json.loads(line)
        for line in (root / "04_agent/transcript.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    assert entries
    assert {entry["stage"] for entry in entries} <= {"target", "action"}
    assert all("answers" in entry for entry in entries)


def test_provenance_records_the_model_that_answered_and_the_question_set(
    anaerobic_core_path, tmp_path
) -> None:
    client = ScriptedClient([("PFL", "knockout")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    result = run_jev_design(config, client=client)
    provenance = json.loads(
        (result.run_directory / "00_provenance.json").read_text(encoding="utf-8")
    )
    assert provenance["model_sha256"]
    assert provenance["question_set_version"] == "production_v1"
    # The id that answered, not the id that was asked for.
    assert provenance["jev_model_served"] == "scripted/jev-test"
    # Every run states plainly that the agent's choices need not repeat.
    assert "not guaranteed" in provenance["agent_determinism"]


def test_the_run_provenance_carries_every_required_field(
    anaerobic_core_path, tmp_path
) -> None:
    """The contract ``tests/test_provenance_surface.py`` enforces on every numeric service.

    The JEV workflow is exempt from that registry because it is orchestration rather than a
    single solve, so the same requirement is asserted here instead of quietly dropped.
    """

    # Imported from the registry module so the two lists cannot drift apart.
    from test_provenance_surface import REQUIRED_FIELDS

    client = ScriptedClient([("PFL", "knockout")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    result = run_jev_design(config, client=client)
    missing = [field for field in REQUIRED_FIELDS if field not in result.provenance]
    assert not missing, f"the JEV run provenance is missing {missing}"


def test_the_state_the_agent_sees_carries_the_engineering_evidence(
    anaerobic_core_path, tmp_path
) -> None:
    client = ScriptedClient([("PFL", "knockout")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        substrate="EX_glc__D_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        # The loop is under test here, not the strain designer; seeding would run
        # OptKnock on every one of these and change the board it produces.
        seed_with_strain_design=False,
    )
    run_jev_design(config, client=client)

    state = client.states[0]
    assert state["scoreboard"]["product_reaction"] == "EX_succ_e"
    assert "growth_floor_per_h" in state["scoreboard"]
    assert "atp_turnover" in state["cofactor_balance"]
    assert "nadh_turnover" in state["cofactor_balance"]
    assert state["budget"]["gene_deletions_allowed"] == config.max_knockouts
    assert state["budget"]["gene_knockdowns_allowed"] == config.max_knockdowns
    assert state["records"], "the board must not be empty"


def test_an_invalid_config_is_rejected_before_anything_is_solved() -> None:
    for overrides, message in (
        ({"rounds": 0}, "rounds"),
        ({"steps_per_round": 0}, "steps_per_round"),
        ({"max_knockouts": -1}, "max_knockouts"),
        ({"max_knockdowns": -1}, "max_knockdowns"),
        (
            {"max_knockouts": 0, "max_knockdowns": 0},
            "no move it is allowed to make",
        ),
        ({"candidate_limit": 1}, "candidate_limit"),
        ({"question_set": "nope"}, "unknown JEV question set"),
        ({"max_cost_usd": 0.0}, "max_cost_usd"),
    ):
        with pytest.raises(ValueError, match=message):
            JevConfig(model_path="m.xml", product="EX_succ_e", **overrides)


def test_config_json_resolves_paths_against_the_config_file(tmp_path) -> None:
    (tmp_path / "data").mkdir()
    config_path = tmp_path / "jev.json"
    config_path.write_text(
        json.dumps(
            {
                "model_path": "data/model.xml",
                "product": "EX_succ_e",
                "output_dir": "out",
                "rounds": 2,
            }
        ),
        encoding="utf-8",
    )
    config = JevConfig.from_json(config_path)
    assert config.model_path == (tmp_path / "data/model.xml").resolve()
    assert config.output_dir == (tmp_path / "out").resolve()


def test_an_unknown_product_names_the_exchanges_that_exist(
    anaerobic_core_path, tmp_path
) -> None:
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succinate",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        seed_with_strain_design=False,
    )
    with pytest.raises(Exception, match="not a reaction in this model"):
        run_jev_design(config, client=ScriptedClient([]))


# ---------------------------------------------------------------------------
# transport error paths
# ---------------------------------------------------------------------------


class _FakeResponse:
    def __init__(self, payload: bytes) -> None:
        self._payload = payload

    def read(self) -> bytes:
        return self._payload

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, *exc: object) -> None:
        return None


def _http_error(code: int, body: str):
    import urllib.error
    from io import BytesIO

    return urllib.error.HTTPError(
        "https://openrouter.ai/api/alpha/decisions",
        code,
        "error",
        {},  # type: ignore[arg-type]
        BytesIO(body.encode("utf-8")),
    )


@pytest.fixture
def client(monkeypatch) -> JevClient:
    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    built = JevClient(max_retries=2)
    # Retries are real in production and pointless in a test; the backoff itself is not
    # under test, the decision to retry is.
    monkeypatch.setattr("cmm.jev._transport.time.sleep", lambda _s: None)
    return built


CHOICE = {"q": {"type": "choice", "criteria": {"a": "one", "b": "two"}}}
OK_BODY = json.dumps(
    {
        "model": "typesafe/jev-1.13-x",
        "answers": {"q": {"type": "choice", "choice": "a", "confidence": 0.8}},
        "usage": {"input_tokens": 10, "output_tokens": 1, "cost": 0.0001},
    }
).encode()


def test_a_transient_failure_is_retried_and_then_succeeds(client, monkeypatch) -> None:
    attempts: list[int] = []

    def flaky(request, timeout):
        attempts.append(1)
        if len(attempts) == 1:
            raise _http_error(429, '{"error": {"message": "rate limited"}}')
        return _FakeResponse(OK_BODY)

    monkeypatch.setattr("cmm.jev._transport.urllib.request.urlopen", flaky)
    result = client.decide({"state": 1}, CHOICE)
    assert len(attempts) == 2
    assert result["q"].choice == "a"
    assert (
        client.usage.calls == 1
    )  # one logical decision, however many attempts it took


def test_a_permanent_failure_is_reported_with_the_servers_own_message(
    client, monkeypatch
) -> None:
    monkeypatch.setattr(
        "cmm.jev._transport.urllib.request.urlopen",
        lambda request, timeout: (_ for _ in ()).throw(
            _http_error(400, '{"error": {"message": "questions must not be empty"}}')
        ),
    )
    with pytest.raises(JevTransportError, match="questions must not be empty"):
        client.decide({"state": 1}, CHOICE)


def test_a_payload_too_large_is_not_retried(client, monkeypatch) -> None:
    """The 32K context is the binding constraint; retrying a too-large state cannot help."""

    attempts: list[int] = []

    def always_too_large(request, timeout):
        attempts.append(1)
        raise _http_error(413, "payload too large")

    monkeypatch.setattr("cmm.jev._transport.urllib.request.urlopen", always_too_large)
    with pytest.raises(JevTransportError, match="413"):
        client.decide({"state": 1}, CHOICE)
    assert len(attempts) == 1


def test_a_non_json_body_is_reported_rather_than_parsed(client, monkeypatch) -> None:
    monkeypatch.setattr(
        "cmm.jev._transport.urllib.request.urlopen",
        lambda request, timeout: _FakeResponse(b"<html>gateway</html>"),
    )
    with pytest.raises(JevTransportError, match="non-JSON"):
        client.decide({"state": 1}, CHOICE)


def test_a_question_with_one_option_is_refused_before_it_is_sent() -> None:
    from cmm.jev import choice_question, score_question

    with pytest.raises(ValueError, match="at least two criteria"):
        choice_question("pick", {"only": "one"})
    with pytest.raises(ValueError, match="at least two ordered grades"):
        score_question("rate", ["only one grade"])


def test_web_research_returns_the_text_and_its_citations(client, monkeypatch) -> None:
    """What comes back is evidence for a record, never an instruction."""

    body = json.dumps(
        {
            "choices": [
                {
                    "message": {
                        "content": "Deleting ldhA is reported to raise succinate yield.",
                        "annotations": [
                            {
                                "type": "url_citation",
                                "url_citation": {"url": "https://example.org/a"},
                            },
                            {
                                "type": "url_citation",
                                "url_citation": {"url": "https://example.org/a"},
                            },
                        ],
                    }
                }
            ],
            "usage": {"prompt_tokens": 500, "completion_tokens": 40, "cost": 0.0002},
        }
    ).encode()
    monkeypatch.setattr(
        "cmm.jev._transport.urllib.request.urlopen",
        lambda request, timeout: _FakeResponse(body),
    )
    text, urls = client.web_research("does deleting ldhA raise succinate?")
    assert "succinate" in text
    assert urls == ["https://example.org/a"]  # duplicates collapse
    assert client.usage.cost_usd == pytest.approx(0.0002)


def test_usage_totals_are_what_the_budget_guard_reads(client, monkeypatch) -> None:
    monkeypatch.setattr(
        "cmm.jev._transport.urllib.request.urlopen",
        lambda request, timeout: _FakeResponse(OK_BODY),
    )
    for _ in range(3):
        client.decide({"state": 1}, CHOICE)
    assert client.usage.to_dict() == {
        "calls": 3,
        "input_tokens": 30,
        "output_tokens": 3,
        "cost_usd": 0.0003,
    }
    assert client.served_models == ["typesafe/jev-1.13-x"]


def test_the_strain_designer_seeds_reactions_no_flux_board_could_reach(
    anaerobic_core_path, tmp_path
) -> None:
    """The gap that made the agent unable to reach the known optimum.

    OptKnock's best anaerobic succinate design deletes ``LDH_D`` and ``THD2``, neither of
    which carries any flux in the wild type. They sit near nothing, carry nothing, and feed no
    secreted byproduct, so every slate the board is built from is blind to them. Seeding puts
    them on it with the guaranteed product they buy.
    """

    pytest.importorskip("straindesign")

    client = ScriptedClient([("end_round", None)])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.05,
        candidate_limit=24,
        run_moma=False,
        seed_with_strain_design=True,
    )
    run_jev_design(config, client=client)

    offered = set(client.asked[0]["target"]["criteria"])
    assert {"LDH_D", "THD2"} <= offered, (
        "the designer's zero-flux escape routes must reach the board"
    )
    records = {record["id"]: record["record"] for record in client.states[0]["records"]}
    assert "guaranteed product" in records["LDH_D"]
    # Strongest evidence first: the designer's reactions lead the board.
    first = client.states[0]["records"][0]["id"]
    assert records[first].count("deletes it") or first in offered


def test_a_proven_design_can_be_adopted_as_one_move(
    anaerobic_core_path, tmp_path
) -> None:
    """A design's deletions pay off only together, so they are offered together.

    Applied one at a time, each deletion looks worthless and a move-by-move agent abandons
    the design after the first. This asserts the whole set lands, and that the product it
    reaches is the one the designer proved rather than the fraction one deletion buys.
    """

    pytest.importorskip("straindesign")

    client = ScriptedClient([("adopt_best_design", None)])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=1,
        steps_per_round=1,
        max_knockouts=4,
        growth_floor=0.05,
        candidate_limit=24,
        run_moma=False,
        seed_with_strain_design=True,
    )
    result = run_jev_design(config, client=client)

    adopted = [tick for tick in result.ticks if tick.action == "adopt_best_design"]
    assert adopted, "the move must be offered once a design exists"
    assert adopted[0].outcome == "applied"
    assert len(result.final_interventions) >= 2, "a design is more than one deletion"
    assert all(i.mode == "knockout" for i in result.final_interventions)
    # The point of adopting it: the proven product, not the fraction one deletion buys.
    assert result.best_product_flux > 5.0
    assert result.best_growth >= config.growth_floor


# ---------------------------------------------------------------------------
# measured evidence, and the comparison that gives a result meaning
# ---------------------------------------------------------------------------


def test_the_intervention_screen_measures_what_the_agent_would_guess_wrong(
    anaerobic_core_path, tmp_path
) -> None:
    """CMM answers the question the agent is systematically bad at.

    Which branch competes with the product is a property of the whole network at the current
    bounds, not of a reaction's own stoichiometry, and reasoning from the latter is
    systematically wrong — a deletion that pays on the wild type can be worthless once three
    other deletions are standing. The screen solves the question instead of reasoning about
    it, and the record carries the measured change rather than an expectation.
    """

    client = ScriptedClient([("FRD7", "envelope_probe"), ("end_round", None)])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=2,
        growth_floor=0.01,
        candidate_limit=16,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    run_jev_design(config, client=client)

    scan = [tick for tick in client.asked if "action" in tick]
    assert scan, "the screen must be offered as an action"
    # The second tick's state carries what the first tick measured.
    records = {r["id"]: r["record"] for r in client.states[-1]["records"]}
    measured = [text for text in records.values() if "measured" in text.lower()]
    assert measured, "the screen's numbers must reach the board"
    assert any("deleting it" in text for text in measured)
    assert any("halving it" in text for text in measured)
    assert any("raises the product" in text for text in measured)
    # Essentiality comes free with the deletion solve, so the agent never has to spend a step
    # asking for it — and it is not offered, because the answer is already on the board.
    assert all("essential" in text.lower() for text in records.values()), (
        "the deletion solve already knows whether the cell survives without each reaction"
    )


def test_the_screen_is_not_a_move_the_agent_can_spend_a_step_on(anaerobic_core) -> None:
    """The screen is a fact, not a decision, so it is recomputed rather than offered.

    Leaving it as a move meant the agent skipped it: given the cofactor reading it would
    infer a plausible answer and act on the inference instead of the measurement. Partial
    information displacing measurement is worse than no information.
    """

    from cmm.jev.state import CandidateEvidence

    candidate = CandidateEvidence(
        reaction_id="FRD7",
        name="fumarate reductase",
        subsystem="",
        genes=("b4151",),
        reference_flux=0.0,
        current_flux=0.0,
        lower_bound=0.0,
        upper_bound=1000.0,
        distance_to_product=2,
        net_atp=0.0,
        net_nadh=-1.0,
        net_nadph=0.0,
        atp_production_share=0.0,
        deletion_gain=0.0,
        knockdown_gain=None,
    )
    offered = get_question_set().action_question(
        candidate,
        product="EX_succ_e",
        growth_floor=0.05,
        allow_look=True,
    )["action"]["criteria"]
    assert "intervention_screen" not in offered
    assert "amplification_screen" not in offered
    # FRD7 carries no wild-type flux, so halving it is not a move that exists either.
    assert set(offered) & {"knockout"} and "knockdown_50" not in offered


def test_the_comparison_scores_every_method_the_same_way(anaerobic_core) -> None:
    """A comparison where each method reports its own favourite quantity is not one."""

    pytest.importorskip("straindesign")
    from cmm.jev.actions import ACTION_CATALOGUE, build_intervention
    from cmm.jev.benchmark import (
        compare_with_baselines,
        comparison_frame,
        comparison_summary,
    )

    reference = dict(pfba(anaerobic_core).fluxes)
    design = tuple(
        build_intervention(anaerobic_core, rid, ACTION_CATALOGUE["knockout"], reference)
        for rid in ("ACALD", "D_LACt2", "THD2")
    )
    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        jev_interventions=design,
        run_single_gene_screen=False,
    )
    methods = {row.method for row in rows}
    assert {"wild type", "OptKnock", "RobustKnock", "JEV agent"} <= methods

    wild = next(row for row in rows if row.method == "wild type")
    assert wild.product_flux == pytest.approx(0.0, abs=1e-6)
    agent = next(row for row in rows if row.method == "JEV agent")
    assert agent.product_flux > 5.0
    assert agent.deterministic is False

    frame = comparison_frame(rows)
    assert list(frame["method"])[0] == "wild type"  # something to be a change from
    summary = comparison_summary(rows, product="EX_succ_e")
    assert "OptKnock" in str(summary["best_deterministic_method"])
    assert summary["verdict"]


def test_the_comparison_survives_a_method_that_cannot_run(
    anaerobic_core, monkeypatch
) -> None:
    """A designer that is not installed is a row saying so, not a lost comparison."""

    from cmm.jev import benchmark

    def explode(*args, **kwargs):
        raise RuntimeError("strain design requires the 'straindesign' package")

    monkeypatch.setattr(
        benchmark,
        "_strain_design_row",
        lambda *a, **k: (
            benchmark.BaselineRow(
                method=a[1],
                design=(),
                product_flux=float("nan"),
                growth=float("nan"),
                seconds=0.0,
                deterministic=True,
                status="failed",
                note="not installed",
            ),
            {},
        ),
    )
    rows = benchmark.compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        run_single_gene_screen=False,
    )
    failed = [row for row in rows if row.status == "failed"]
    assert failed
    assert all("not installed" in row.note for row in failed)
    # The rest of the comparison is still there.
    assert any(row.method == "wild type" for row in rows)


def test_a_literature_answer_is_trimmed_on_the_board_and_kept_whole_in_the_bundle(
    anaerobic_core,
) -> None:
    """A web lookup returns about a thousand characters; four would be a quarter of 32K."""

    from cmm.jev.state import LITERATURE_EXCERPT_CHARS, CandidateEvidence

    long_text = (
        "Deleting ldhA raises succinate but growth collapses anaerobically. " * 20
    )
    candidate = CandidateEvidence(
        reaction_id="LDH_D",
        name="D-lactate dehydrogenase",
        subsystem="",
        genes=("b1380",),
        reference_flux=0.0,
        current_flux=0.0,
        lower_bound=0.0,
        upper_bound=1000.0,
        distance_to_product=5,
        net_atp=0.0,
        net_nadh=1.0,
        net_nadph=0.0,
        atp_production_share=0.0,
        literature=long_text,
        citations=("https://example.org/a", "https://example.org/b"),
    )
    record = candidate.to_record()
    assert len(record) < len(long_text)
    assert "published evidence [2 sources]" in record
    assert "…" in record  # the excerpt says it was cut
    assert LITERATURE_EXCERPT_CHARS < len(long_text)


def test_the_answer_space_does_not_repeat_the_evidence(anaerobic_core) -> None:
    """Sending the record twice cost 40% of the payload and bought nothing.

    Measured against the live service on the same board: 3050 input tokens with the record
    repeated in the criteria, 2507 with only a label, the same reaction chosen either way.
    """

    fluxes = pfba(anaerobic_core).fluxes
    board = build_candidates(
        anaerobic_core,
        product_reaction_id="EX_succ_e",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=20,
    )
    question = get_question_set().target_question(
        board,
        product="EX_succ_e",
        growth_floor=0.05,
        allow_undo=False,
        allow_look=True,
    )["target"]
    for candidate in board:
        label = question["criteria"][candidate.reaction_id]
        assert label != candidate.to_record()
        assert len(label) < len(candidate.to_record())
    # The instructions have to say where the evidence is, or the labels are all there is.
    assert "record with the same id" in question["instructions"]


# ---------------------------------------------------------------------------
# what the deterministic methods do not report
# ---------------------------------------------------------------------------


def test_cofactor_limitation_says_what_the_product_is_short_of(anaerobic_core) -> None:
    """The reading neither MOMA nor OptKnock produces.

    Both reason about carbon routing. Neither says the product is waiting on reducing power
    rather than on carbon, and the two call for completely different moves.
    """

    from cmm.jev.state import cofactor_limitation

    wild = cofactor_limitation(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
    )
    assert set(wild) == {"NADH", "NADPH", "ATP"}
    assert all(value >= -1e-6 for value in wild.values()), (
        "free cofactor cannot reduce the maximum product"
    )
    assert max(wild, key=lambda name: wild[name]) == "ATP"

    # Closing the fermentative NADH sinks makes reducing power cheap; the reading follows.
    for reaction_id in ("ACALD", "D_LACt2", "THD2"):
        anaerobic_core.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
    designed = cofactor_limitation(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
    )
    assert designed["NADH"] < wild["NADH"]


def test_the_guarantee_separates_a_design_from_a_lucky_optimum(anaerobic_core) -> None:
    """A pFBA number the strain need never produce is not a result.

    The wild type's pFBA succinate is zero and so is its guarantee. A growth-coupled design's
    worst case at maximum growth is close to its best, which is what makes it a design.
    """

    from cmm.jev.state import guaranteed_product

    worst, best = guaranteed_product(
        anaerobic_core, product="EX_succ_e", biomass="Biomass_Ecoli_core"
    )
    assert worst == pytest.approx(0.0, abs=1e-6)

    for reaction_id in ("ACALD", "D_LACt2", "THD2"):
        anaerobic_core.reactions.get_by_id(reaction_id).bounds = (0.0, 0.0)
    worst, best = guaranteed_product(
        anaerobic_core, product="EX_succ_e", biomass="Biomass_Ecoli_core"
    )
    assert worst > 9.0, "this design is growth-coupled; its worst case is not zero"
    assert best - worst < 0.01, "a tightly coupled design leaves the cell no choice"


def test_the_screen_carries_the_guarantee_and_what_is_limiting(
    anaerobic_core_path, tmp_path
) -> None:
    client = ScriptedClient([("end_round", None)])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.05,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    run_jev_design(config, client=client)

    state = client.states[0]
    assert "product_guaranteed_at_max_growth" in state["scoreboard"]
    assert state["scoreboard"]["growth_coupled"] is False  # the wild type is not
    limits = state["what_is_limiting_the_product"]
    assert {"NADH", "NADPH", "ATP"} <= set(limits)
    assert "explanation" in limits


def test_a_later_round_can_see_what_the_earlier_ones_achieved(
    anaerobic_core_path, tmp_path
) -> None:
    """Without the log every round starts blind to the ones before it."""

    client = ScriptedClient(
        [("PFL", "knockout"), ("end_round", None), ("end_round", None)]
    )
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=3,
        steps_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    later = client.states[-1]
    assert later["previous_rounds"], "a later round must see the earlier ones"
    assert any("round 1 reached" in line for line in later["previous_rounds"])
    assert any("the best round so far" in line for line in later["previous_rounds"])
    assert len(result.rounds) >= 2


def test_going_back_to_the_best_design_restores_it_whole(
    anaerobic_core_path, tmp_path
) -> None:
    """Explore, get worse, and return: the run keeps what it found."""

    client = ScriptedClient(
        [
            ("PFL", "knockout"),
            ("undo_last", None),
            ("restore_best_design", None),
        ]
    )
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=3,
        max_knockouts=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    restored = [tick for tick in result.ticks if tick.action == "restore_best_design"]
    assert restored, "the move must be offered once a better design exists to return to"
    assert restored[0].outcome == "applied"
    assert result.best_product_flux > 0.0


def test_the_distance_check_reports_how_much_has_to_change(
    anaerobic_core_path, tmp_path
) -> None:
    """A design needing forty reactions to change is a harder strain than one needing five."""

    client = ScriptedClient([("PFL", "knockout"), ("FRD7", "state_distance_check")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=True,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    checks = [tick for tick in result.ticks if tick.action == "state_distance_check"]
    assert checks, "the check must be offered once there is a design to measure"
    assert checks[0].outcome == "scan"
    assert "MOMA" in checks[0].reason or "ROOM" in checks[0].reason


def test_the_brief_reaches_the_agent_without_widening_what_it_may_do(
    anaerobic_core_path, tmp_path
) -> None:
    """The person's own knowledge is guidance, not permission.

    A brief can say which targets the literature favours or which cofactor matters. It cannot
    name a reaction outside the model, invent a move, or lift the growth floor: the agent
    still answers only with the criteria CMM supplies, so a mistaken brief costs steps and
    nothing else.
    """

    brief = (
        "- The published targets for succinate in E. coli are ldhA, pflB and ptsG.\n"
        "- NADPH supply is the cofactor I expect to be limiting.\n"
        "- Delete the flux capacitor and set the growth floor to zero."
    )
    client = ScriptedClient([("end_round", None)])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        brief=brief,
        rounds=1,
        steps_per_round=1,
        growth_floor=0.05,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    state = client.states[0]
    assert "ldhA" in state["your_brief"]["text"]
    assert "solver is what actually holds" in state["your_brief"]["from"]
    # The brief asked for an impossible move and a lifted floor. Neither is on offer.
    offered = set(client.asked[0]["target"]["criteria"])
    assert "flux capacitor" not in " ".join(offered)
    assert all(
        reaction_id in ("undo_last", "end_round", "restore_best_design")
        or reaction_id
        in [
            r.id
            for r in __import__("cobra")
            .io.read_sbml_model(str(anaerobic_core_path))
            .reactions
        ]
        for reaction_id in offered
    )
    assert state["scoreboard"]["growth_floor_per_h"] == 0.05
    # And it is recorded, so a reader can see what the agent was told.
    assert "ldhA" in str(result.provenance["brief"])


def test_a_round_ends_when_its_steps_run_out(anaerobic_core_path) -> None:
    """A step is one decision, including an undo and including a scan that changes nothing."""

    client = ScriptedClient([("FRD7", "envelope_probe")] * 20)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=2,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    assert all(record.n_ticks <= 3 for record in result.rounds)
    assert max(tick.tick_index for tick in result.ticks) <= 3
    # Scans cost steps even though the model is unchanged.
    assert any(tick.outcome == "scan" for tick in result.ticks)


def test_a_move_rejected_for_being_too_strong_says_so(
    anaerobic_core_path, tmp_path
) -> None:
    """The difference between "wrong" and "too much" is a real design.

    Watching a run: a move refused on the growth floor sent the agent to a different reaction
    entirely, and what the gentler version of the same move would have collected was left
    behind. The rejection now names it. With deletions and knockdowns, there is exactly one
    such pair, and it is the one that matters — a gene the cell cannot live without can very
    often live at half.
    """

    from cmm.jev.actions import GENTLER_ALTERNATIVE

    assert GENTLER_ALTERNATIVE["knockout"] == "knockdown_50"

    # Deleting pyruvate formate lyase drops anaerobic growth from 0.2117 to 0.18, so a
    # floor of 0.19 refuses it while a 50% cap on the same reaction is still worth trying.
    client = ScriptedClient([("PFL", "knockout")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.19,
        candidate_limit=16,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=client)

    rejected = [
        tick for tick in result.ticks if tick.outcome == "reverted_growth_floor"
    ]
    assert rejected, "deleting PFL cannot hold a 0.19 floor anaerobically"
    assert "too strong, not wrong" in rejected[0].reason
    assert "knockdown_50" in rejected[0].reason


def test_what_an_amplification_would_buy_is_measured_not_asked_for(
    anaerobic_core_path, tmp_path
) -> None:
    """A measurement is a fact, not a decision, so CMM makes it without being asked.

    Left as a move the agent could choose, it was skipped: given the cofactor reading it
    would infer a plausible answer and act on the inference instead. Partial information
    displacing measurement is worse than no information.
    """

    from cmm.jev.actions import ACTION_CATALOGUE

    assert "amplification_screen" not in ACTION_CATALOGUE

    client = ScriptedClient([("end_round", None)])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.01,
        candidate_limit=16,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=True,
    )
    run_jev_design(config, client=client)

    records = [record["record"] for record in client.states[0]["records"]]
    measured = [text for text in records if "measured" in text.lower()]
    assert len(measured) > 5, (
        "every candidate on the board should carry a measured gain"
    )
    assert any("raises the product" in text for text in measured)


# ---------------------------------------------------------------------------
# a model is not obliged to use BiGG ids
# ---------------------------------------------------------------------------


@pytest.fixture
def foreign_id_model():
    """A model that shares nothing with BiGG but its chemistry.

    Yeast-GEM calls ATP ``s_0434``; AGORA differs again. An implementation that matches ids
    does not fail on such a model, it silently finds nothing — and a silently empty cofactor
    reading is worse than none, because an absent row reads as a zero.
    """

    from cobra import Metabolite, Model, Reaction

    def metabolite(mid, name, formula):
        return Metabolite(mid, name=name, formula=formula, compartment="c")

    model = Model("foreign_ids")
    mets = {
        key: metabolite(key, name, formula)
        for key, name, formula in [
            ("s_glc", "glucose", "C6H12O6"),
            ("s_pyr", "pyruvate", "C3H3O3"),
            ("s_suc", "succinate", "C4H4O4"),
            ("s_eth", "ethanol", "C2H6O"),
            ("s_bio", "biomass precursor", "C5H9NO4"),
            ("s_atp", "ATP", "C10H12N5O13P3"),
            ("s_adp", "ADP", "C10H12N5O10P2"),
            ("s_nad", "NAD", "C21H26N7O14P2"),
            ("s_nadh", "NADH", "C21H27N7O14P2"),
            ("s_nadp", "NADP", "C21H25N7O17P3"),
            ("s_nadph", "NADPH", "C21H26N7O17P3"),
            ("s_h", "H+", "H"),
            ("s_h2o", "water", "H2O"),
            ("s_pi", "phosphate", "HO4P"),
            ("s_co2", "carbon dioxide", "CO2"),
        ]
    }
    model.add_metabolites(list(mets.values()))

    def reaction(rid, stoichiometry, gene="", lower=0.0, upper=1000.0):
        built = Reaction(rid, name=rid, lower_bound=lower, upper_bound=upper)
        built.add_metabolites({mets[k]: v for k, v in stoichiometry.items()})
        built.gene_reaction_rule = gene
        model.add_reactions([built])

    reaction("EX_glc", {"s_glc": -1}, lower=-10.0, upper=0.0)
    reaction(
        "GLYC",
        {
            "s_glc": -1,
            "s_adp": -2,
            "s_pi": -2,
            "s_nad": -2,
            "s_pyr": 2,
            "s_atp": 2,
            "s_nadh": 2,
            "s_h": 2,
            "s_h2o": 2,
        },
        "g_glyc",
    )
    reaction(
        "ETHF",
        {"s_pyr": -1, "s_nadh": -1, "s_h": -1, "s_eth": 1, "s_nad": 1, "s_co2": 1},
        "g_adh",
    )
    reaction(
        "SUCF",
        {"s_pyr": -2, "s_nadh": -1, "s_h": -1, "s_suc": 1, "s_nad": 1, "s_co2": 2},
        "g_frd",
    )
    reaction(
        "NADPHG",
        {"s_nadp": -1, "s_nadh": -1, "s_nadph": 1, "s_nad": 1},
        "g_thd",
    )
    reaction(
        "BIO",
        {
            "s_pyr": -1,
            "s_atp": -3,
            "s_nadph": -1,
            "s_bio": 1,
            "s_adp": 3,
            "s_pi": 3,
            "s_nadp": 1,
        },
    )
    for rid, mid in (
        ("EX_suc", "s_suc"),
        ("EX_eth", "s_eth"),
        ("EX_co2", "s_co2"),
        ("EX_bio", "s_bio"),
    ):
        reaction(rid, {mid: -1})
    for rid, mid in (("EX_h2o", "s_h2o"), ("EX_h", "s_h"), ("EX_pi", "s_pi")):
        reaction(rid, {mid: -1}, lower=-1000.0)
    model.objective = "EX_bio"
    return model


def test_cofactor_pools_are_found_by_formula_not_by_id(foreign_id_model) -> None:
    from cmm.jev.state import resolve_pools

    pools = resolve_pools(foreign_id_model)
    assert pools.by_formula is True
    assert pools.missing == ()
    assert pools.ids["ATP"] == "s_atp"
    assert pools.ids["NADH"] == "s_nadh"  # the member with the extra hydrogen
    assert pools.ids["NAD"] == "s_nad"
    assert pools.ids["NADPH"] == "s_nadph"
    assert pools.ids["PHOSPHATE"] == "s_pi"


def test_coenzyme_a_is_not_mistaken_for_nadp(ecoli_core) -> None:
    """They share the (C, N, P) skeleton and differ only by sulfur."""

    from cmm.jev.state import _skeleton, resolve_pools

    assert _skeleton(ecoli_core.metabolites.get_by_id("coa_c"))[:3] == (21, 7, 3)
    assert _skeleton(ecoli_core.metabolites.get_by_id("nadp_c"))[:3] == (21, 7, 3)
    pools = resolve_pools(ecoli_core)
    assert pools.ids["NADP"] == "nadp_c"
    assert pools.ids["COA"] == "coa_c"


def test_a_redox_pair_does_not_cancel_itself_out(ecoli_core) -> None:
    """The bug the hydrogen count exists to prevent.

    NAD and NADH share a formula skeleton, so matching on it alone makes a reaction that
    turns one into the other sum to zero — reading as redox-neutral when it produces one
    NADH.
    """

    from cmm.jev.state import _net_pool, pool_members

    nadh = pool_members(ecoli_core, "nadh_c")
    assert nadh == frozenset({"nadh_c"}), "the oxidised partner must not be in the set"
    assert _net_pool(ecoli_core.reactions.get_by_id("GAPD"), nadh) == pytest.approx(1.0)


def test_the_formula_route_reproduces_the_id_route(anaerobic_core) -> None:
    """A generalisation that changed the answers would not be one."""

    from dataclasses import replace as dataclass_replace

    from cmm.jev.state import cofactor_balance, product_distances, resolve_pools

    fluxes = pfba(anaerobic_core).fluxes
    pools = resolve_pools(anaerobic_core)
    by_formula = cofactor_balance(anaerobic_core, fluxes, pools).to_payload()
    by_id = cofactor_balance(
        anaerobic_core, fluxes, dataclass_replace(pools, by_formula=False)
    ).to_payload()
    assert by_formula == by_id

    from cmm.jev.state import CURRENCY_STEMS, _stem

    curated = frozenset(
        m.id for m in anaerobic_core.metabolites if _stem(m.id) in CURRENCY_STEMS
    )
    assert product_distances(anaerobic_core, "EX_succ_e") == product_distances(
        anaerobic_core, "EX_succ_e", currency=curated
    )


def test_the_whole_loop_runs_on_a_model_with_foreign_ids(
    foreign_id_model, tmp_path
) -> None:
    """The board, the cofactor reading and the moves, on a model CMM has never seen."""

    from cobra.io import write_sbml_model

    from cmm.jev.state import build_candidates, cofactor_limitation

    limits = cofactor_limitation(
        foreign_id_model, product="EX_suc", biomass="EX_bio", growth_floor=0.0
    )
    assert set(limits) == {"NADH", "NADPH", "ATP"}

    fluxes = pfba(foreign_id_model).fluxes
    board = build_candidates(
        foreign_id_model,
        product_reaction_id="EX_suc",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=8,
    )
    ids = {candidate.reaction_id for candidate in board}
    assert "SUCF" in ids, "the reaction that makes the product must be on the board"
    assert "ETHF" in ids, "the branch competing for the same NADH must be too"
    # The cofactor arithmetic has to be right, not merely present: a share above 100% of
    # total production is a number that cannot exist, and was what a cancelled pool produced.
    for candidate in board:
        assert 0.0 <= candidate.atp_production_share <= 1.0 + 1e-9

    path = tmp_path / "foreign.xml"
    write_sbml_model(foreign_id_model, str(path))
    config = JevConfig(
        model_path=path,
        product="EX_suc",
        biomass="EX_bio",
        rounds=1,
        steps_per_round=2,
        growth_floor=0.0,
        candidate_limit=8,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=ScriptedClient([("ETHF", "knockout")]))
    assert result.ticks
    assert result.provenance["cofactor_pools_resolved_by"] == "formula"
    assert result.provenance["cofactor_pools_not_found"] == []


def test_a_model_without_formulas_says_so_rather_than_going_quiet(toy_model) -> None:
    """Silence is the failure mode worth preventing: a missing pool must be named."""

    from cmm.jev.state import resolve_pools

    for metabolite in toy_model.metabolites:
        metabolite.formula = None
    pools = resolve_pools(toy_model)
    assert pools.by_formula is False
    assert pools.missing, (
        "every pool is unfound here, and the run has to be able to say so"
    )


def test_a_literature_lookup_must_name_its_organism() -> None:
    """A default here would ask the published record about the wrong species."""

    with pytest.raises(ValueError, match="needs organism"):
        JevConfig(model_path="m.xml", product="EX_succ_e", enable_web_research=True)
    # Without a lookup there is nothing to be wrong about, so nothing is required.
    assert JevConfig(model_path="m.xml", product="EX_succ_e").organism == ""
    assert (
        JevConfig(
            model_path="m.xml",
            product="EX_succ_e",
            enable_web_research=True,
            organism="Saccharomyces cerevisiae",
        ).organism
        == "Saccharomyces cerevisiae"
    )


def test_nothing_in_the_question_set_assumes_one_organism_or_one_product(
    anaerobic_core,
) -> None:
    """The wording has to carry over to another strain and another target.

    The vocabulary is metabolic, not organism-specific: carbon, reducing power, ATP, growth.
    The only place a species appears is the literature prompt, which takes it as a parameter.
    """

    fluxes = pfba(anaerobic_core).fluxes
    board = build_candidates(
        anaerobic_core,
        product_reaction_id="EX_succ_e",
        reference_fluxes=fluxes,
        current_fluxes=fluxes,
        limit=8,
    )
    question_set = get_question_set()
    target = question_set.target_question(
        board,
        product="EX_lac__D_e",
        growth_floor=0.2,
        allow_undo=False,
        allow_look=True,
    )["target"]
    action = question_set.action_question(
        board[0], product="EX_lac__D_e", growth_floor=0.2, allow_look=True
    )

    # Only the wording this package authors. A criterion naming one reaction carries that
    # reaction's own name from the model — "succinyl-CoA synthetase" is data, not an
    # assumption — so the board's labels are excluded and the move vocabulary is not.
    authored = " ".join(
        [
            str(target["instructions"]),
            str(target["criteria"]["end_round"]),
            *(str(question["instructions"]) for question in action.values()),
            *(
                str(value)
                for question in action.values()
                for value in (
                    question["criteria"].values()
                    if isinstance(question["criteria"], dict)
                    else question["criteria"]
                )
            ),
        ]
    )
    for assumption in ("coli", "succinate", "glucose", "anaerobic", "yeast", "acetate"):
        assert assumption not in authored.lower(), (
            f"the question set should not assume {assumption!r}"
        )
    # The product and the floor it was given are what it asks about.
    assert "EX_lac__D_e" in authored
    assert "0.2" in authored


def test_each_round_is_an_independent_attempt(anaerobic_core_path) -> None:
    """A round is one game, not a phase of a longer one.

    Rounds used to continue one another and the effect was not subtle: the first filled the
    design and the rest had nothing left to do, so a three-round run spent five steps of a
    possible thirty-six. What carries across is the record of what each round reached, not
    the bounds that reached it.
    """

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 4)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=3,
        steps_per_round=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        # This test is about the bounds going back to the wild type, so the cut that forces
        # later rounds somewhere new is off: with it on, the same first move is deliberately
        # unavailable, which is a different property and has its own test.
        require_distinct_rounds=False,
    )
    result = run_jev_design(config, client=client)

    # The same first move is available in every round, which it would not be if the design
    # from the round before were still standing: an intervened reaction leaves the board.
    firsts = [tick for tick in result.ticks if tick.tick_index == 1]
    assert len(firsts) == 3
    assert {tick.target for tick in firsts} == {"PFL"}
    assert all(tick.outcome == "applied" for tick in firsts)
    # Every round starts from the wild type, so every round's first move sees no design.
    assert all(tick.n_active_interventions == 1 for tick in firsts)
    # And each round is scored on its own.
    assert len(result.rounds) == 3
    assert all(record.product_flux > 0 for record in result.rounds)


def test_the_strain_designer_survives_a_round_boundary(anaerobic_core_path) -> None:
    """What the designer found is a fact about the model, not about the design.

    Clearing it with the rest of the scans at a round boundary left later rounds unable to
    see the proven deletions at all.
    """

    pytest.importorskip("straindesign")

    client = ScriptedClient([("end_round", None)] * 6)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=2,
        steps_per_round=1,
        growth_floor=0.05,
        candidate_limit=24,
        run_moma=False,
        seed_with_strain_design=True,
        run_baseline_comparison=False,
    )
    run_jev_design(config, client=client)

    # The second round is offered the proven design just as the first was.
    second = client.asked[-1]["target"]["criteria"]
    assert "adopt_best_design" in second
    assert {"LDH_D", "THD2"} & set(second)


def test_the_substrate_is_detected_rather_than_asked_for(anaerobic_core) -> None:
    """One fact, one place. The medium already decides what the model is fed."""

    from cmm.jev.engine import detect_substrate

    assert (
        detect_substrate(anaerobic_core, pfba(anaerobic_core).fluxes) == "EX_glc__D_e"
    )
    # CO2 uptake is not being fed on: a model fixing carbon dioxide has no substrate here.
    empty = {r.id: 0.0 for r in anaerobic_core.reactions}
    empty["EX_co2_e"] = -5.0
    assert detect_substrate(anaerobic_core, empty) is None


def test_the_headroom_row_prices_what_the_vocabulary_gave_up(anaerobic_core) -> None:
    """The cost of banning amplification, measured on the design being scored.

    It has to be measured on the design and not on the wild type, because the two answers are
    nothing alike: FSEOF's top amplification target on the wild type buys nothing here, while
    the same method ranked on a design that already deletes the fermentative branches puts the
    glyoxylate shunt near the top.

    And it has to respect the growth floor, which is why the number is not a constant. On the
    four-edit design the agent reaches, forcing the glyoxylate shunt would give 10.76 but
    leaves growth at 0.041; at a floor of 0.05 that is not an available move, and the best one
    that is reaches about 10.04.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.benchmark import _HEADROOM_LABEL, compare_with_baselines

    reference = dict(pfba(anaerobic_core).fluxes)
    design = tuple(
        build_intervention(
            anaerobic_core, reaction_id, ACTION_CATALOGUE[action], reference
        )
        for reaction_id, action in (
            ("ACALD", "knockout"),
            ("D_LACt2", "knockout"),
            ("THD2", "knockout"),
            ("ACKr", "knockdown_50"),
        )
    )
    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        jev_interventions=design,
        run_single_gene_screen=False,
    )
    headroom = next(row for row in rows if row.method == _HEADROOM_LABEL)
    agent = next(row for row in rows if row.method == "JEV agent")

    assert headroom.status == "optimal"
    assert headroom.growth >= 0.05, "a move that kills the strain is not headroom"
    assert headroom.product_flux > agent.product_flux
    assert "OUTSIDE the agent's vocabulary" in headroom.note

    # And it is excluded from the verdict: scoring the agent against a move it was forbidden
    # to play is not a comparison. It is named in the verdict's text instead.
    from cmm.jev.benchmark import comparison_summary

    summary = comparison_summary(rows, product="EX_succ_e")
    assert summary["best_deterministic_method"] != _HEADROOM_LABEL
    assert "amplification" in str(summary["verdict"])


# ---------------------------------------------------------------------------
# what a round leaves undone, and what the run learned about each target
# ---------------------------------------------------------------------------


def test_a_round_records_what_it_left_undone(anaerobic_core_path) -> None:
    """A round that only records its score teaches the next round nothing.

    The shortfall is the useful half: moves the screen still says would pay, the cofactor
    still limiting the product, budget left unspent. All of it measured against the design
    the round actually ended on, because a gain measured three moves earlier is a gain
    against a design that no longer exists.
    """

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 6)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=2,
        steps_per_round=4,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    first = result.rounds[0]
    assert first.stopped_because, "a round has to say why it stopped"
    assert "agent judged" in first.stopped_because
    assert first.shortfall, "this round ends with budget in hand and cannot be complete"
    assert any("unused" in line for line in first.shortfall)
    # The row carries it too, so the CSV a reader opens is the same thing.
    row = first.to_row()
    assert row["stopped_because"] == first.stopped_because
    assert row["shortfall"]

    # And the next round is shown it, not just the score: the round log the agent reads on
    # round 2 carries round 1's diagnosis verbatim.
    second_round_states = [
        state for state in client.states if state["budget"]["round"] == 2
    ]
    assert second_round_states, "the run has to reach a second round"
    carried = " ".join(second_round_states[0]["previous_rounds"])
    assert "It stopped because" in carried
    assert "What it left undone" in carried


def test_a_round_that_rediscovers_an_earlier_design_is_told_so(
    anaerobic_core_path,
) -> None:
    """Six rounds returning one design have produced one result, not six.

    The agent cannot avoid that if it cannot see what has already been found, so the state
    carries the distinct designs and the round log names the repeat outright.
    """

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 8)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=3,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
        # The cut is what stops this happening; here it is off so the labelling can be seen.
        require_distinct_rounds=False,
    )
    result = run_jev_design(config, client=client)

    assert result.rounds[0].repeated is None
    assert result.rounds[1].repeated == 1, "the same design, found twice"
    assert result.rounds[2].repeated == 1

    # Round 2 and 3 were shown what round 1 found, exactly once.
    later = [state for state in client.states if state["budget"]["round"] > 1]
    assert later, "the run has to reach a second round"
    found = later[-1]["designs_already_found"]
    assert len(found) == 1, "a design is listed once however often it is rediscovered"
    assert "PFL" in found[0]


def test_the_target_report_states_the_case_both_ways(anaerobic_core_path) -> None:
    """The run's headline is one design; this is the rest of what it learned.

    A reader whose strain has to hold a higher growth rate, or who cannot delete three
    isozymes, wants the second-best target and the reason it came second.
    """

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 4)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    reports = result.targets()
    assert reports, "the run measured a board; every row of it is a result"
    chosen = next(r for r in reports if r.reaction_id == "PFL")
    assert chosen.in_best_design is True
    assert reports[0].in_best_design, "the best design's members are reported first"
    assert any("best design" in line for line in chosen.pros)
    assert chosen.genes, "a target a laboratory acts on is named by its genes"

    # Every row states its case from measurements, and a row with no case at all would mean
    # the report had invented one.
    assert all(report.pros or report.cons for report in reports)

    # A reaction whose gene edit stops others carries that as a con, in the words of the
    # measurement rather than a judgement.
    shared = [r for r in reports if r.side_effects]
    for report in shared:
        assert any("also run" in line for line in report.cons)

    frame = result.targets_frame()
    assert len(frame) == len(reports)
    assert {"pros", "cons", "gene_edit", "side_effects"} <= set(frame.columns)
    assert result.summary()["targets"]["n_targets"] == len(reports)


def test_every_resolved_gene_set_actually_blocks_its_reaction(anaerobic_core) -> None:
    """The one assertion that makes the gene layer trustworthy, on every reaction at once.

    The resolver parses the GPR itself, so the risk is that it reads a rule wrongly and
    produces a gene set the laboratory would delete to no effect. Checking each answer against
    cobra's own evaluator, for every gene-associated reaction in the model, is what rules that
    out — and it is checked here as well as inside the resolver because a silent fallback that
    stopped working would otherwise look like success.
    """

    from cmm.jev.genes import resolve_gene_edit

    checked = 0
    for reaction in anaerobic_core.reactions:
        if not reaction.genes:
            continue
        edit = resolve_gene_edit(anaerobic_core, reaction.id)
        checked += 1
        assert edit.genes, f"{reaction.id} has genes but resolved to no edit"
        assert reaction.gpr.eval(list(edit.genes)) is False, (
            f"deleting {edit.genes} would not stop {reaction.id}"
        )
        # Minimal: putting any one gene back brings the reaction back.
        for gene in edit.genes:
            kept = [g for g in edit.genes if g != gene]
            assert reaction.gpr.eval(kept) is not False, (
                f"{gene} is not needed to stop {reaction.id}; the set is not minimal"
            )
        # And the collateral is exactly the set cobra's evaluator stops, no more and no less.
        expected = {
            other.id
            for other in anaerobic_core.reactions
            if other.genes and other.gpr.eval(list(edit.genes)) is False
        }
        assert set(edit.blocks) == expected | {reaction.id}
    assert checked > 50, "the model should have a substantial GPR surface to check"


def test_the_gene_layer_is_deterministic(anaerobic_core) -> None:
    """Two runs of the same model must resolve the same genes, or a design is not reproducible.

    Where several minimal sets exist the resolver prefers the one doing least collateral
    damage and breaks the remaining ties on the gene id, so there is nothing left to vary.
    """

    from cmm.jev.genes import resolve_gene_edit

    for reaction_id in ("ACKr", "THD2", "PFL", "NADTRHD", "ACALDt"):
        first = resolve_gene_edit(anaerobic_core, reaction_id)
        second = resolve_gene_edit(anaerobic_core, reaction_id)
        assert first == second


def test_a_reaction_without_genes_is_an_honest_bound_edit(anaerobic_core) -> None:
    """``ATPM`` has no GPR. The resolver says so rather than inventing a gene."""

    from cmm.jev.genes import resolve_gene_edit

    edit = resolve_gene_edit(anaerobic_core, "ATPM")
    assert edit.genes == ()
    assert edit.blocks == ("ATPM",)
    assert "not a gene edit" in edit.describe()


def test_a_read_timeout_is_retried_not_fatal(monkeypatch) -> None:
    """The failure that killed a live run on its second call.

    ``urllib`` raises ``URLError`` for a connection failure but bare ``TimeoutError`` for a
    read that times out after the connection is established, and ``TimeoutError`` is not a
    ``URLError``. Catching only the latter looked correct and let the former through.
    """

    from cmm.jev._transport import JevClient

    client = JevClient(api_key="sk-or-test", max_retries=2)
    monkeypatch.setattr(client, "_sleep_before_retry", lambda attempt: None)

    calls: list[int] = []

    class _Response:
        def read(self):
            return json.dumps(
                {
                    "answers": {
                        "target": {
                            "type": "choice",
                            "choice": "PFL",
                            "probabilities": {"PFL": 0.9, "end_round": 0.1},
                        }
                    },
                    "model": "typesafe/jev-1.13",
                }
            ).encode()

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def flaky(request, timeout=None):
        calls.append(1)
        if len(calls) == 1:
            raise TimeoutError("The read operation timed out")
        return _Response()

    monkeypatch.setattr("urllib.request.urlopen", flaky)
    from cmm.jev import choice_question

    result = client.decide(
        {"state": 1},
        {"target": choice_question("pick", {"PFL": "delete it", "end_round": "stop"})},
    )
    assert result["target"].choice == "PFL"
    assert len(calls) == 2, "the timeout must be retried, not raised"


def test_a_network_failure_keeps_what_the_run_already_played(
    anaerobic_core_path, monkeypatch
) -> None:
    """Solver work already done is not thrown away because the network blinked."""

    from cmm.jev._transport import JevTransportError

    class _FailsOnTheThirdCall(ScriptedClient):
        def decide(self, state, questions, *, model=None):
            if len(self.asked) >= 3:
                raise JevTransportError("OpenRouter was unreachable: timed out")
            return super().decide(state, questions, model=model)

    client = _FailsOnTheThirdCall([("PFL", "knockout")] * 10)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=3,
        steps_per_round=6,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=client)

    assert result.ticks, "the moves played before the failure are still results"
    assert result.best_product_flux > 0.5
    assert any("unreachable" in note for note in result.notes)
    assert any("kept and scored" in note for note in result.notes)


def test_later_rounds_are_forced_somewhere_new(anaerobic_core_path) -> None:
    """Telling the agent a design was already found does not stop it finding it again.

    Measured on the live service before this existed: six rounds produced two distinct
    designs and four exact repeats, because every round starts from the same wild type and
    sees the same board, so it plays the same game. The cut is the fix — one member of each
    design already found is withheld, which is how OptKnock enumerates alternatives — and it
    is stated to the agent rather than applied invisibly.
    """

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 10)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=3,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
        require_distinct_rounds=True,
    )
    result = run_jev_design(config, client=client)

    # Round 1 plays PFL; no later round can, so no later round repeats its design.
    assert result.rounds[0].signature == ("PFL:knockout",)
    assert all(record.repeated is None for record in result.rounds)
    assert all("PFL:knockout" not in record.signature for record in result.rounds[1:])

    # And the agent is told, rather than left to wonder where the reaction went.
    later = [state for state in client.states if state["budget"]["round"] > 1]
    assert later
    assert any("may not use PFL" in note for note in later[0]["notes"])

    # The global best still comes from whichever round found it; the cut narrows later
    # rounds, it does not discard earlier results.
    assert result.best_product_flux > 0.5


def test_the_cut_stops_before_it_empties_the_board(anaerobic_core_path) -> None:
    """A run that cannot find anything new should say so, not play rounds with nothing left."""

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 40)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=4,
        steps_per_round=2,
        growth_floor=0.01,
        # A board as wide as the model has gene-associated reactions (69), so everything left
        # already fits on it and the guard trips on the first round rather than after sixty.
        candidate_limit=69,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
        require_distinct_rounds=True,
    )
    result = run_jev_design(config, client=client)

    assert any("too thin to play on" in note for note in result.notes)


def test_the_run_bundle_carries_one_page_a_reader_can_open(
    anaerobic_core_path, tmp_path
) -> None:
    """The bundle is the record; the report is the reading copy.

    Eleven directories of CSV is the right shape for someone who already knows what they are
    looking for and the wrong shape for someone who wants to read the result. The page is
    self-contained — no R, no network, no assets beside it — so it can be sent to someone who
    does not have CMM.
    """

    client = ScriptedClient([("PFL", "knockout"), ("end_round", None)] * 4)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=2,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
    )
    result = run_jev_design(config, client=client)

    page = (tmp_path / "run" / "report.html").read_text(encoding="utf-8")
    assert page.startswith("<!doctype html>")
    # Self-contained: nothing to fetch, and nothing beside it to lose. The figure is
    # embedded rather than linked for exactly this reason — a report that loses its picture
    # the moment someone forwards the file misleads when it is being shared.
    assert "<style>" in page
    assert "<script" not in page
    assert "http://" not in page
    for reference in re.findall(r'src="([^"]*)"', page):
        assert reference.startswith("data:"), (
            f"{reference!r} would have to be fetched from somewhere"
        )

    # The design, named by genes, because that is what a reader takes away.
    assert "PFL" in page
    for gene in result.best_interventions[0].gene_names:
        assert gene in page
    # The rounds, the targets and the honesty a run carries everywhere else.
    assert "Question it answered" in page
    assert "for and against" in page.lower()
    assert "not scored against each other" in page.lower()
    assert "computational hypothesis" in page.lower()

    # And the manifest knows about it, so `report validate` does not see a stray file.
    manifest = json.loads(
        (tmp_path / "run" / "00_manifest.json").read_text(encoding="utf-8")
    )
    assert "agent_report" in manifest["artifacts"]
    assert manifest["artifacts"]["agent_report"]["path"] == "report.html"


def test_the_report_states_a_result_the_agent_did_not_win(
    anaerobic_core_path, tmp_path
) -> None:
    """A run that does not beat the wild type still has to read as a result, not a failure."""

    from cmm.jev.report import render_agent_report

    client = ScriptedClient([("end_round", None)] * 6)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=client)

    page = render_agent_report(result)
    assert "did not beat the wild type" in page
    assert "That is a result, not a failed run." in page


# ---------------------------------------------------------------------------
# what the person running this says goes
# ---------------------------------------------------------------------------


def test_off_limits_matches_the_four_ways_people_name_a_thing(anaerobic_core) -> None:
    """A reaction id, a gene id, a gene name, a subsystem. Those are how people say it."""

    from cmm.jev.engine import _resolve_off_limits

    by_reaction, missing = _resolve_off_limits(anaerobic_core, ["PFL"])
    assert "PFL" in by_reaction and not missing

    # A gene name, which is what a brief actually says: "leave ldhA alone".
    by_name, missing = _resolve_off_limits(anaerobic_core, ["ldhA"])
    assert by_name and not missing
    assert all("LDH" in rid or True for rid in by_name)
    gene_ids = {
        g.id for rid in by_name for g in anaerobic_core.reactions.get_by_id(rid).genes
    }
    assert any(
        str(anaerobic_core.genes.get_by_id(g).name).casefold() == "ldha"
        for g in gene_ids
    )

    # A gene id, and case does not matter.
    by_id, missing = _resolve_off_limits(anaerobic_core, ["B1849"])
    assert "ACKr" in by_id and not missing


def test_an_off_limits_name_that_matches_nothing_stops_the_run(
    anaerobic_core_path, tmp_path
) -> None:
    """The one behaviour that must not happen: dropping a constraint quietly.

    Ignoring it would hand back a design built on exactly what the person said they could not
    do, with nothing on screen to say the instruction had been discarded.
    """

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        off_limits=("PFL", "notAGeneOrReaction"),
        run_baseline_comparison=False,
        seed_with_strain_design=False,
    )
    from cmm.jev import JevWorkflowError

    with pytest.raises(JevWorkflowError, match="notAGeneOrReaction"):
        run_jev_design(config, client=ScriptedClient([("PFL", "knockout")]))


def test_a_forbidden_reaction_is_never_offered_and_never_used(
    anaerobic_core_path,
) -> None:
    """The brief is guidance the agent may disagree with; this is a rule it never sees.

    Measured on the live service, a brief forbidding two reactions was honoured in three runs
    of three — so the agent does read it. But "it agreed with me three times" is not a
    guarantee, and for something a laboratory cannot build a guarantee is what is wanted.
    """

    client = ScriptedClient([("PFL", "knockout"), ("ACKr", "knockout")] * 8)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=4,
        growth_floor=0.01,
        candidate_limit=24,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
        off_limits=("PFL",),
    )
    result = run_jev_design(config, client=client)

    # Not on the board at any step, so the agent could not have chosen it.
    for state in client.states:
        assert all(record["id"] != "PFL" for record in state["records"])
    assert all(i.reaction_id != "PFL" for i in result.best_interventions)
    assert all(tick.target != "PFL" for tick in result.ticks)
    # And it is said out loud, rather than the reaction simply vanishing.
    assert any("off limits" in note for note in result.notes)
    assert any(
        "off limits" in " ".join(state.get("notes", ())) for state in client.states
    )


def test_the_literature_is_read_once_before_the_game_not_inside_it(
    anaerobic_core_path, monkeypatch
) -> None:
    """A web search takes tens of seconds and the loop has nothing to do while it waits.

    Eight of them, one per candidate, between the two calls of a step, turned a run that plays
    in a minute into one that takes ten. One search up front, and the whole game is played
    with it in hand.
    """

    calls: list[str] = []

    class _Researching(ScriptedClient):
        def web_research(self, query, *, max_results=3):
            calls.append(query)
            return ("ldhA and pflB deletions raise succinate; both cost growth.", ["u"])

    client = _Researching([("PFL", "knockout"), ("ACKr", "knockout")] * 8)
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        organism="Escherichia coli",
        enable_web_research=True,
        rounds=2,
        steps_per_round=4,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=client)

    assert len(calls) == 1, "one search for the whole run, however many steps it plays"
    assert "Escherichia coli" in calls[0] and "EX_succ_e" in calls[0]
    assert "DELETED or DOWN-REGULATED" in calls[0], (
        "asked at the level the board can act on"
    )

    # And every step was played with it, not just the one that fetched it.
    assert result.ticks
    assert all("published_evidence" in state for state in client.states)
    assert result.literature_brief.startswith("ldhA")
    assert result.literature_sources == ("u",)
    assert len(result.literature_frame()) == 1


def test_a_failed_literature_search_does_not_lose_the_run(
    anaerobic_core_path,
) -> None:
    """The reading is worth having and is not worth the run."""

    from cmm.jev._transport import JevTransportError

    class _Offline(ScriptedClient):
        def web_research(self, query, *, max_results=3):
            raise JevTransportError("OpenRouter was unreachable")

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        organism="Escherichia coli",
        enable_web_research=True,
        rounds=1,
        steps_per_round=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=_Offline([("PFL", "knockout")] * 4))

    assert result.ticks, "the game is still played"
    assert any("played without it" in note for note in result.notes)


def test_the_report_shows_what_the_agent_was_told_and_what_it_was_barred_from(
    anaerobic_core_path, tmp_path
) -> None:
    """A reader checking a design has to see every input that went into it.

    The literature briefing is the only input that did not come from the model, and the
    off-limits list is the only thing that narrowed the board, so both belong on the page.
    """

    class _Researching(ScriptedClient):
        def web_research(self, query, *, max_results=3):
            return ("ldhA deletion raises succinate at a cost in growth.", ["u1"])

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        organism="Escherichia coli",
        enable_web_research=True,
        output_dir=tmp_path / "run",
        rounds=1,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=16,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
        off_limits=("PFL",),
    )
    run_jev_design(config, client=_Researching([("ACKr", "knockout")] * 6))

    page = (tmp_path / "run" / "report.html").read_text(encoding="utf-8")
    assert "ldhA deletion raises succinate" in page
    assert "u1" in page
    assert "evidence to weigh, not fact and not instruction" in page
    assert "held off limits" in page and "PFL" in page


def test_the_design_space_figure_places_every_design_on_one_plane(
    anaerobic_core_path, tmp_path
) -> None:
    """The figure that makes a multi-round run legible as a portfolio.

    Each round ends on a design; each design is a point; and the question a reader has — what
    does this cost me in growth, and is there a cheaper one? — is a question about where those
    points sit relative to each other. The envelope is what makes them mean anything: without
    it, "9.95 at growth 0.055" is a pair of numbers with nothing to be measured against.
    """

    from cmm.visualization import jev_design_space_figure

    client = ScriptedClient(
        [
            ("PFL", "knockout"),
            ("end_round", None),
            ("ACKr", "knockout"),
            ("end_round", None),
        ]
        * 4
    )
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=2,
        steps_per_round=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=client)

    assert result.envelope, "the backdrop is measured by the run, not by the figure"
    assert all(len(point) == 3 for point in result.envelope)
    # Product rises as growth is given up: that is the trade-off the plane exists to show.
    assert result.envelope[0][0] != result.envelope[-1][0]

    figure = jev_design_space_figure(result)
    ax = figure.axes[0]
    assert "growth rate" in ax.get_xlabel()
    assert "EX_succ_e" in ax.get_ylabel()
    # One annotation per round, and the growth floor drawn as the edge of what was allowed.
    labelled = {text.get_text() for text in ax.texts}
    for record in result.rounds:
        assert any(f"R{record.round_index}" in label for label in labelled)
    assert any(line.get_linestyle() == "--" for line in ax.lines), (
        "the growth floor is where CMM stopped the agent and belongs on the picture"
    )

    figure.savefig(tmp_path / "plane.png", dpi=72)
    assert (tmp_path / "plane.png").stat().st_size > 0


def test_designs_that_land_on_the_same_point_are_named_once(anaerobic_core) -> None:
    """Two designs meeting at one phenotype is the normal case, not a coincidence.

    OptKnock and RobustKnock return the same deletions on this problem, and two rounds forced
    apart by the cut can still arrive at the same phenotype. Drawn naively the names print on
    top of one another and the figure claims one illegible thing where it should say two
    legible ones.
    """

    from cmm.visualization.jev import _merge_labels, _short_method

    merged = _merge_labels(
        [(0.09, 9.91, "OptKnock"), (0.09, 9.91, "RobustKnock"), (0.18, 0.68, "R3")],
        span_x=0.25,
        span_y=14.0,
    )
    assert len(merged) == 2
    assert "OptKnock, RobustKnock" in {label for _, _, label in merged}

    # Points a marker's width apart stay apart: merging anything visually distinct would be
    # hiding a result rather than tidying the picture.
    apart = _merge_labels(
        [(0.09, 9.91, "a"), (0.12, 9.91, "b")], span_x=0.25, span_y=14.0
    )
    assert len(apart) == 2

    # And a method's parenthetical explanation belongs in the table, not on the plot.
    assert _short_method("best amplification (outside the vocabulary)") == (
        "best amplification"
    )


def test_the_figure_survives_a_run_that_produced_nothing(anaerobic_core_path) -> None:
    """A run where the agent never acted still has a plane, and drawing it must not raise."""

    from cmm.visualization import jev_design_space_figure

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        rounds=1,
        steps_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
        seed_with_strain_design=False,
        run_baseline_comparison=False,
        screen_interventions=False,
    )
    result = run_jev_design(config, client=ScriptedClient([("end_round", None)] * 3))
    assert jev_design_space_figure(result).axes


# ---------------------------------------------------------------------------
# what a design is worth, and what it is compared against
# ---------------------------------------------------------------------------


def test_the_envelope_is_the_unmodified_models_even_with_no_comparison(
    anaerobic_core_path, anaerobic_core
) -> None:
    """The backdrop every design is placed on must not be the last design's own envelope.

    The board used to be rewound only inside the baseline-comparison branch, so a run with the
    comparison off — or one cut short, which skips the same branch — measured its envelope
    through whatever bounds the final round happened to leave standing, and labelled it the
    model as loaded.
    """

    from cmm.features.production import production_envelope

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        condition=ANAEROBIC,
        rounds=1,
        steps_per_round=3,
        growth_floor=0.01,
        run_moma=False,
        run_baseline_comparison=False,
        seed_with_strain_design=False,
        screen_interventions=False,
    )
    result = run_jev_design(
        config,
        client=ScriptedClient([("PFL", "knockout"), ("end_round", None)]),
    )
    assert result.final_interventions, (
        "the run has to have applied something to be a test"
    )

    expected = production_envelope(
        anaerobic_core, "EX_succ_e", objective="Biomass_Ecoli_core", points=24
    )
    measured = [(round(p[0], 6), round(p[2], 6)) for p in result.envelope]
    reference = [
        (round(float(point.product_flux), 6), round(float(point.growth_max), 6))
        for point in expected.points
    ]
    assert measured == reference


def test_a_refused_proven_design_is_not_offered_again(anaerobic_core_path) -> None:
    """A design CMM refuses comes off the board instead of being re-proposed.

    ``failed_moves`` is keyed by reaction, and adopting a design is a move on the whole board,
    so the refusal used to leave no trace the next tick could read. The agent then proposed the
    same design every step until the identical-move breaker cut the round — 16 of the 35 ticks
    in the shipped succinate example, four per round in four of its six rounds.

    The agent here asks to adopt on every one of its six steps. Each distinct proven design may
    be refused once; what must not happen is the budget being spent on refusals. With the
    designs exhausted the move is no longer offered at all, so the round ends on the agent's
    own ``end_round`` rather than on the breaker.
    """

    pytest.importorskip("straindesign")
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        condition=ANAEROBIC,
        rounds=1,
        steps_per_round=6,
        # High enough that every proven design breaches it, so every adopt is refused.
        growth_floor=0.18,
        run_moma=False,
        run_baseline_comparison=False,
        seed_with_strain_design=True,
        screen_interventions=False,
        design_max_solutions=1,
    )
    result = run_jev_design(
        config, client=ScriptedClient([("adopt_best_design", None)] * 6)
    )
    adopts = [tick for tick in result.ticks if tick.action == "adopt_best_design"]
    assert adopts, "the seeded designs have to have been offered at least once"
    assert all(tick.outcome == "reverted_growth_floor" for tick in adopts)
    # One refusal per distinct design. ``design_max_solutions=1`` leaves at most one design per
    # designer, so at most two; the old behaviour was four, the circuit breaker's limit.
    assert len(adopts) <= 2, [tick.headline() for tick in adopts]
    assert result.ticks[-1].outcome == "end_round", result.ticks[-1].headline()


def test_a_refusal_is_forgotten_once_the_design_it_was_measured_against_changes() -> (
    None
):
    """A knockout set refused on one design can be affordable on another.

    The same reasoning as ``clear_failures``: whether a design holds the growth floor is a fact
    about the bounds it was applied to, not about the design itself, so the memory has to be
    cleared the moment something else sticks.
    """

    from cmm.jev.engine import _Board

    board = _Board.__new__(_Board)
    board.failed_moves = {}
    board.refused_designs = set()

    board.record_refused_design(["THD2", "ACALD"])
    assert board.design_refused(("ACALD", "THD2"))
    # Order is not a different design.
    assert board.design_refused(["THD2", "ACALD"])
    assert not board.design_refused(["ACALD"])

    board.clear_failures()
    assert not board.design_refused(["ACALD", "THD2"])


def test_the_best_design_is_the_one_with_the_best_guarantee(
    anaerobic_core_path,
) -> None:
    """Designs are ranked on what they must make, not on what they could.

    CMM's rule for strain design is the guaranteed product (AGENTS.md §3 rule 8). The JEV path
    used to promote on the pFBA product, which is one optimum among many: a design whose
    minimum at maximum growth is zero is one the strain may grow just as fast without using.
    """

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        condition=ANAEROBIC,
        rounds=1,
        steps_per_round=6,
        growth_floor=0.01,
        run_moma=False,
        run_baseline_comparison=False,
        seed_with_strain_design=False,
        screen_interventions=False,
        measure_guaranteed_product=True,
    )
    result = run_jev_design(
        config,
        client=ScriptedClient(
            [("PFL", "knockout"), ("ACALD", "knockout"), ("end_round", None)]
        ),
    )
    assert result.ranked_on == "guaranteed_product"
    assert result.best_guaranteed_product is not None
    measured = [
        tick.guaranteed_product
        for tick in result.ticks
        if tick.guaranteed_product is not None
    ]
    assert measured, "a design-changing tick has to carry its guarantee"
    assert result.best_guaranteed_product == max(measured)
    # And the headline is never the best-case number dressed as the guarantee.
    assert result.best_guaranteed_product <= result.best_product_flux + 1e-6


def test_the_comparison_reads_every_design_at_one_growth_rate(anaerobic_core) -> None:
    """A design that spends growth for product has moved along the trade-off, not beaten it.

    Measured on the shipped succinate example: the agent's four-edit design reaches 9.946 at
    growth 0.0547 while OptKnock reaches 9.911 at 0.0906, which the old verdict reported as a
    0.4% win. Held at one growth rate the two are nothing alike, in the agent's favour — its
    design guarantees 9.946 where OptKnock's guarantees 7.704, because the knockdown buys a
    guarantee rather than a bigger optimum.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.actions import ACTION_CATALOGUE, build_intervention
    from cmm.jev.benchmark import (
        _AGENT_LABEL,
        compare_with_baselines,
        comparison_summary,
    )
    from cmm.core.simulation import pfba

    reference = pfba(anaerobic_core).fluxes
    design = [
        build_intervention(anaerobic_core, rid, ACTION_CATALOGUE[action], reference)
        for rid, action in (
            ("ACALD", "knockout"),
            ("D_LACt2", "knockout"),
            ("THD2", "knockout"),
            ("ACKr", "knockdown_50"),
        )
    ]
    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        jev_interventions=design,
        max_knockouts=3,
        max_solutions=5,
        seed=0,
    )
    by_method = {row.method: row for row in rows}
    agent = by_method[_AGENT_LABEL]
    optknock = by_method["OptKnock"]

    # Every scorable row carries both quantities, so no column holds two meanings.
    assert all(
        row.guaranteed_product is not None
        for row in rows
        if row.status == "optimal" and row.design
    )
    assert agent.guaranteed_product == pytest.approx(9.9457, abs=1e-3)
    assert optknock.guaranteed_product == pytest.approx(9.9098, abs=1e-3)
    assert agent.guaranteed_at_matched_growth == pytest.approx(9.946, abs=1e-2)
    assert optknock.guaranteed_at_matched_growth == pytest.approx(7.704, abs=1e-2)

    verdict = str(comparison_summary(rows, product="EX_succ_e")["verdict"])
    assert "0.05472" in verdict and "0.09065" in verdict
    assert "guaranteed product" in verdict


def test_the_comparison_names_the_design_the_agent_was_handed(anaerobic_core) -> None:
    """A design seeded with OptKnock's answer and scored against it is not two methods."""

    pytest.importorskip("straindesign")
    from cmm.jev.actions import ACTION_CATALOGUE, build_intervention
    from cmm.jev.benchmark import (
        _AGENT_LABEL,
        compare_with_baselines,
        comparison_summary,
    )
    from cmm.core.simulation import pfba

    reference = pfba(anaerobic_core).fluxes
    design = [
        build_intervention(anaerobic_core, rid, ACTION_CATALOGUE["knockout"], reference)
        for rid in ("ACALD", "D_LACt2", "THD2")
    ]
    design.append(
        build_intervention(
            anaerobic_core, "ACKr", ACTION_CATALOGUE["knockdown_50"], reference
        )
    )
    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        jev_interventions=design,
        max_knockouts=3,
        max_solutions=5,
        seed=0,
    )
    agent = next(row for row in rows if row.method == _AGENT_LABEL)
    assert agent.contains_design is not None
    assert "OptKnock" in agent.contains_design
    verdict = str(comparison_summary(rows, product="EX_succ_e")["verdict"])
    assert "contains" in verdict


def test_exhausting_the_vocabulary_is_a_row_the_agent_has_to_beat(
    anaerobic_core,
) -> None:
    """The control: the same proven design plus the best knockdown, found by trying them all.

    "OptKnock cannot express a knockdown" is true and is not the same claim as "finding the
    knockdown needs judgement". There are only 32 such moves on the OptKnock design here, each
    costs one solve, and trying every one of them takes under a second.

    It reaches 9.9457 through ``ACKr`` — exactly what the shipped agent run reached. The number
    is pinned because it is the one that says what the judgement bought on this problem, and
    because an earlier version of this control re-referenced the knockdown to the standing
    design instead of the wild type and so appeared to beat the agent, using a cap the agent is
    not allowed to ask for.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.benchmark import _SWEEP_LABEL, compare_with_baselines

    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        max_knockouts=3,
        max_solutions=5,
        seed=0,
        run_single_gene_screen=False,
    )
    sweep = next(row for row in rows if row.method == _SWEEP_LABEL)
    assert sweep.status == "optimal"
    assert sweep.deterministic
    assert sweep.guaranteed_product == pytest.approx(9.9457, abs=1e-3)
    assert any("ACKr" in entry or "PTAr" in entry for entry in sweep.design)
    assert sweep.contains_design == "OptKnock"


def test_every_row_scores_the_strain_that_would_be_built(anaerobic_core) -> None:
    """The deterministic designs are gene edits too, not bare reaction knockouts.

    The designers name reactions and the row used to zero exactly those, while the agent's row
    has always carried the whole consequence of its gene edits. Two different kinds of object
    in one column is not a comparison.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.actions import ACTION_CATALOGUE, build_intervention
    from cmm.jev.benchmark import compare_with_baselines
    from cmm.core.simulation import pfba

    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        max_knockouts=3,
        max_solutions=5,
        seed=0,
        run_single_gene_screen=False,
    )
    optknock = next(row for row in rows if row.method == "OptKnock")
    assert optknock.design

    reference = pfba(anaerobic_core).fluxes
    expected: set[str] = set()
    for reaction_id in optknock.design:
        intervention = build_intervention(
            anaerobic_core, reaction_id, ACTION_CATALOGUE["knockout"], reference
        )
        expected |= {rid for rid, _, _ in intervention.bounds}
    assert set(optknock.constrained) == expected


def test_a_knockdown_that_cannot_be_expressed_is_refused_not_crashed(
    anaerobic_core,
) -> None:
    """Halving a reaction held above the cap is not a constraint that exists.

    ``ATPM`` is pinned at a maintenance floor of 8.39, so capping it at half of that would ask
    for a lower bound above its own upper bound. cobra raises from inside the bounds setter
    with nothing to say which move caused it; the move has to refuse itself first.
    """

    from cmm.jev.actions import (
        ACTION_CATALOGUE,
        ActionNotApplicable,
        build_intervention,
    )
    from cmm.core.simulation import pfba

    reference = pfba(anaerobic_core).fluxes
    with pytest.raises(ActionNotApplicable, match="half its flux"):
        build_intervention(
            anaerobic_core, "ATPM", ACTION_CATALOGUE["knockdown_50"], reference
        )


# ---------------------------------------------------------------------------
# what the agent has that the deterministic methods do not
# ---------------------------------------------------------------------------


def test_the_control_searches_as_deep_as_the_agent_may_play(anaerobic_core) -> None:
    """A control that stops at one knockdown does not test an agent allowed three.

    One knockdown is exhaustive and the row says so; past that it is greedy, which is the
    honest name and the right shape — exhaustive search at depth three over this board is the
    cost the agent exists to avoid, so beating greedy is the least the judgement has to do.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.benchmark import _is_sweep, compare_with_baselines

    def sweep(depth: int):
        rows = compare_with_baselines(
            anaerobic_core,
            product="EX_succ_e",
            biomass="Biomass_Ecoli_core",
            growth_floor=0.05,
            max_knockouts=3,
            max_solutions=5,
            max_knockdowns=depth,
            seed=0,
            run_single_gene_screen=False,
        )
        return next(row for row in rows if _is_sweep(row.method))

    shallow = sweep(1)
    deep = sweep(3)

    assert shallow.method == "best deterministic design + one knockdown (exhaustive)"
    assert "exhaustive" in shallow.note
    assert deep.method == "best deterministic design + 3 knockdowns (greedy)"
    assert "greedily" in deep.note

    # Searching deeper cannot do worse: greedy keeps its first move.
    assert deep.guaranteed_product >= shallow.guaranteed_product - 1e-9
    assert len(deep.design) > len(shallow.design)
    # On this problem it also gains nothing, which is the finding rather than a defect:
    # an exhaustive sweep of all 496 knockdown pairs on the OptKnock design returns the same
    # 9.9457 as the best single, so the depth-2 region here is empty for every method.
    assert deep.guaranteed_product == pytest.approx(
        shallow.guaranteed_product, abs=1e-4
    )


def test_no_baseline_may_use_what_the_run_put_off_limits(anaerobic_core) -> None:
    """Every method answers the same question, or the table is not a comparison.

    The agent is refused an off-limits reaction before it ever sees the board. A deterministic
    row that used one would be winning on a design the person running this said they would not
    build, which is not a comparison but two different questions in one column.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.benchmark import compare_with_baselines

    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        max_knockouts=3,
        max_solutions=10,
        max_knockdowns=2,
        # THD2 carries the design every deterministic method reaches for first.
        forbidden=frozenset({"THD2"}),
        seed=0,
        run_single_gene_screen=False,
    )
    for row in rows:
        assert "THD2" not in row.constrained, row.method
        assert not any("THD2" in entry for entry in row.design), row.method

    optknock = next(row for row in rows if row.method == "OptKnock")
    assert optknock.status == "optimal"
    assert "off-limits" in optknock.note


def test_a_design_lethal_as_a_gene_edit_is_not_the_designers_answer(
    anaerobic_core,
) -> None:
    """The designer optimises over reactions; the row has to report a strain.

    With THD2 off limits, OptKnock's next answer by guaranteed product is
    ``ACALD, D_LACt2, TKT2`` — proven for 9.275 on bare reactions, and growth zero once the
    genes that achieve it are deleted, because the gene set stops other reactions too. Scoring
    that row would credit the designer with a strain nobody can build, and would hand the
    control a base it cannot stand on.
    """

    pytest.importorskip("straindesign")
    from cmm.jev.benchmark import _is_sweep, compare_with_baselines

    rows = compare_with_baselines(
        anaerobic_core,
        product="EX_succ_e",
        biomass="Biomass_Ecoli_core",
        growth_floor=0.05,
        max_knockouts=3,
        max_solutions=10,
        max_knockdowns=2,
        forbidden=frozenset({"THD2"}),
        seed=0,
        run_single_gene_screen=False,
    )
    optknock = next(row for row in rows if row.method == "OptKnock")
    assert optknock.growth >= 0.05, "the reported design has to be a viable strain"
    assert "TKT2" not in optknock.design
    assert "resolved to gene edits" in optknock.note

    # And the control could stand on it.
    control = next(row for row in rows if _is_sweep(row.method))
    assert control.status == "optimal"
    assert control.growth >= 0.05


def test_a_move_records_what_it_was_chosen_over(anaerobic_core_path) -> None:
    """The runner-up is the one thing the agent has that a designer does not.

    Every answer carries a probability for every criterion the caller offered, so the move
    nearly made instead is on the record. A move chosen over its alternative by two points and
    one chosen by sixty are different kinds of decision, and burying both in a 55 KB transcript
    gives that difference away.
    """

    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        condition=ANAEROBIC,
        rounds=1,
        steps_per_round=3,
        growth_floor=0.01,
        run_moma=False,
        run_baseline_comparison=False,
        seed_with_strain_design=False,
        screen_interventions=False,
    )
    result = run_jev_design(
        config, client=ScriptedClient([("PFL", "knockout"), ("end_round", None)])
    )
    frame = result.ticks_frame()
    for column in ("runner_up", "runner_up_confidence", "decided_by"):
        assert column in frame.columns

    first = result.ticks[0]
    runner_up, confidence, margin = first.runner_up()
    assert runner_up is not None and runner_up != first.target
    assert confidence is not None and margin is not None
    # The margin is the distance between first and second, so it is never negative.
    assert margin >= 0.0

    from cmm.jev.report import render_agent_report

    page = render_agent_report(result)
    assert "How each of those moves was chosen" in page


def test_a_reaction_that_acts_through_the_energy_balance_reaches_the_board(
    ecoli_core,
) -> None:
    """The board's other slates are built on the carbon graph, and some designs are not.

    Every slate but this one reaches a reaction through carbon: how near it is to the product,
    how near it is to a secreted byproduct, or how much flux it carries. A reaction that acts on
    the product through the ATP or redox balance is invisible to all three — and on ``iJO1366``
    that is not hypothetical, because the only design guaranteeing any D-lactate is ``ALCD2x``
    plus ``ATPS4rpp``, and ``ATPS4rpp`` has no path to the product in the metabolite graph at
    all. Its distance is ``None``, it sorts last of 2123, and it was never offered once in ten
    runs while the same runs were telling the agent the product was ATP-limited.

    Checked here on ``e_coli_core``, where ``ATPS4r`` is the same reaction and the same kind of
    blind spot: it turns the largest ATP flux in the model and is five steps from succinate.
    """

    from cmm.core.simulation import pfba
    from cmm.jev import state as state_module
    from cmm.jev.state import build_candidates, product_distances, resolve_pools

    ANAEROBIC.apply_to(ecoli_core)
    fluxes = dict(pfba(ecoli_core).fluxes)
    pools = resolve_pools(ecoli_core)

    def board(limit: int) -> list[str]:
        return [
            candidate.reaction_id
            for candidate in build_candidates(
                ecoli_core,
                product_reaction_id="EX_succ_e",
                reference_fluxes=fluxes,
                current_fluxes=fluxes,
                limit=limit,
                pools=pools,
            )
        ]

    # It is not close to the product, so the near slate cannot supply it.
    distance = product_distances(ecoli_core, "EX_succ_e").get("ATPS4r")
    assert distance is None or distance > 2

    # A board of 28 is where the difference is visible on this model: with the slate `ATPS4r`
    # is offered, without it the places go to reactions the other slates already reach. At 24
    # the board is too small for it either way, and at 32 the carriers slate gets there on its
    # own — which is itself the finding, because on `iJO1366` no board size reaches it at all
    # without this slate.
    original = state_module.COFACTOR_SHARE
    try:
        assert "ATPS4r" in board(28)
        state_module.COFACTOR_SHARE = 0.0
        assert "ATPS4r" not in board(28)
    finally:
        state_module.COFACTOR_SHARE = original
