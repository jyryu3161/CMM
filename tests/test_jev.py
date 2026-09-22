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
from pathlib import Path

import pytest
from cobra.io import write_sbml_model

from cmm.core.condition import Condition, ReactionBound
from cmm.core.simulation import pfba
from cmm.jev import (
    ACTION_CATALOGUE,
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
from cmm.jev.actions import FORCE_ON_ACTIONS, feasible_extreme
from cmm.jev.questions import NoAvailableAction
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


def test_a_missing_api_key_names_the_variable_to_set(monkeypatch) -> None:
    monkeypatch.delenv("OPENROUTER_API_KEY", raising=False)
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


def test_a_reaction_carrying_flux_gets_the_relative_moves(anaerobic_core) -> None:
    """One knockdown strength, not two.

    A second, deeper cap mostly bought a second rejection of the same idea, and roughly
    halving an activity is the level a promoter swap or an RBS change can actually aim at.
    """

    names = [action.name for action in applicable_actions(8.2)]
    assert names == ["knockout", "knockdown_50", "amplify_2x", "amplify_5x"]
    assert "knockdown_25" not in ACTION_CATALOGUE


def test_a_reaction_at_zero_gets_the_switch_on_moves_instead() -> None:
    """Without these the agent could never start a pathway that is off."""

    names = [action.name for action in applicable_actions(0.0)]
    assert names == ["knockout", "force_on_low", "force_on_high"]
    assert not {"amplify_2x", "knockdown_50"} & set(names)


def test_a_knockdown_caps_the_magnitude_without_opening_a_direction(
    anaerobic_core,
) -> None:
    reaction = anaerobic_core.reactions.get_by_id("PFL")
    reaction.bounds = (0.0, 1000.0)
    intervention = build_intervention(
        anaerobic_core, "PFL", ACTION_CATALOGUE["knockdown_50"], reference_flux=17.8
    )
    assert intervention.upper_bound == pytest.approx(8.9)
    assert intervention.lower_bound == 0.0  # the reaction was irreversible; it stays so


def test_a_relative_move_on_a_zero_flux_reaction_is_refused_not_reinterpreted(
    anaerobic_core,
) -> None:
    with pytest.raises(ActionNotApplicable, match="nothing to scale"):
        build_intervention(
            anaerobic_core, "FRD7", ACTION_CATALOGUE["amplify_2x"], reference_flux=0.0
        )


def test_switching_a_reaction_on_uses_its_loop_free_maximum(anaerobic_core) -> None:
    """The headline reason the extreme must be loopless.

    A plain LP maximisation of ``FRD7`` returns its 1000 bound through the
    ``FRD7``/``SUCDi`` cycle, on a model taking up 10 mmol gDW^-1 h^-1 of glucose. Sixty per
    cent of that would be a physically meaningless target that the solver would satisfy with
    a futile cycle.
    """

    low, high = feasible_extreme(anaerobic_core, "FRD7")
    assert low == pytest.approx(0.0, abs=1e-6)
    assert 5.0 < high < 30.0, (
        "a loop-free maximum must be on the scale of the carbon input"
    )

    intervention = build_intervention(
        anaerobic_core, "FRD7", ACTION_CATALOGUE["force_on_high"], reference_flux=0.0
    )
    assert intervention.exploratory is True
    assert intervention.lower_bound == pytest.approx(0.6 * high)
    assert "loop-free maximum" in intervention.describe()


def test_a_switch_on_move_is_refused_on_a_reaction_that_already_carries_flux(
    anaerobic_core,
) -> None:
    with pytest.raises(ActionNotApplicable, match="already carries"):
        build_intervention(
            anaerobic_core, "PFL", FORCE_ON_ACTIONS[0], reference_flux=17.8
        )


def test_an_intervention_converts_to_the_bound_cmm_already_applies(
    anaerobic_core,
) -> None:
    intervention = build_intervention(
        anaerobic_core, "PFL", ACTION_CATALOGUE["knockout"], reference_flux=17.8
    )
    bound = intervention.to_reaction_bound()
    assert isinstance(bound, ReactionBound)
    assert (bound.lower_bound, bound.upper_bound) == (0.0, 0.0)


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
            exclude={"knockout", "force_on_low", "force_on_high"},
        )


# ---------------------------------------------------------------------------
# a whole game, offline
# ---------------------------------------------------------------------------


def test_a_scripted_game_switches_the_pathway_on_and_records_every_move(
    anaerobic_core_path, tmp_path
) -> None:
    """The end-to-end result: a product that was zero is being made, and the run says how."""

    # Switching succinyl-CoA synthetase on opens the reductive route to succinate. It is a
    # single move on purpose: the assertion is that the engine turns a chosen move into real
    # product flux, not that this particular design is the best one.
    client = ScriptedClient([("SUCOAS", "force_on_high")])
    config = JevConfig(
        model_path=anaerobic_core_path,
        product="EX_succ_e",
        substrate="EX_glc__D_e",
        biomass="Biomass_Ecoli_core",
        output_dir=tmp_path / "run",
        rounds=2,
        steps_per_round=3,
        max_interventions=3,
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
    assert result.best_product_flux > 1.0, (
        "forcing the reductive branch must make succinate"
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


def test_the_intervention_cap_is_never_exceeded(anaerobic_core_path, tmp_path) -> None:
    client = ScriptedClient(
        [
            ("FRD7", "force_on_high"),
            ("FUM", "force_on_high"),
            ("SUCOAS", "force_on_high"),
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
        max_interventions=2,
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
    client = ScriptedClient([("FRD7", "force_on_high"), ("FUM", "force_on_low")])
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
    client = ScriptedClient([("FRD7", "force_on_high")])
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

    client = ScriptedClient([("FRD7", "force_on_high")])
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
    client = ScriptedClient([("FRD7", "force_on_high")])
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
    assert state["budget"]["max_interventions"] == config.max_interventions
    assert state["records"], "the board must not be empty"


def test_an_invalid_config_is_rejected_before_anything_is_solved() -> None:
    for overrides, message in (
        ({"rounds": 0}, "rounds"),
        ({"steps_per_round": 0}, "steps_per_round"),
        ({"max_interventions": 0}, "max_interventions"),
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
        max_interventions=4,
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


def test_the_amplification_screen_measures_what_the_agent_would_guess_wrong(
    anaerobic_core_path, tmp_path
) -> None:
    """CMM answers the question the agent is systematically bad at.

    Asked which reaction to amplify for succinate, the agent reaches for fumarate reductase —
    the direct product-forming step, which is already saturated and buys nothing. The screen
    solves the question instead of reasoning about it, and the record carries the measured
    change rather than an expectation.
    """

    client = ScriptedClient([("FRD7", "amplification_screen"), ("end_round", None)])
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
    assert any("raises the product" in text for text in measured)


def test_a_screened_reaction_is_not_screened_again(anaerobic_core) -> None:
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
        amplification_gain=0.0,
    )
    offered = get_question_set().action_question(
        candidate,
        product="EX_succ_e",
        growth_floor=0.05,
        allow_look=True,
    )["action"]["criteria"]
    assert "amplification_screen" not in offered


def test_the_comparison_scores_every_method_the_same_way(anaerobic_core) -> None:
    """A comparison where each method reports its own favourite quantity is not one."""

    pytest.importorskip("straindesign")
    from cmm.jev.actions import ACTION_CATALOGUE, build_intervention
    from cmm.jev.benchmark import (
        compare_with_baselines,
        comparison_frame,
        comparison_summary,
    )

    design = tuple(
        build_intervention(
            anaerobic_core, rid, ACTION_CATALOGUE["knockout"], reference_flux=0.0
        )
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
        lambda *a, **k: benchmark.BaselineRow(
            method=a[1],
            design=(),
            product_flux=float("nan"),
            growth=float("nan"),
            seconds=0.0,
            deterministic=True,
            status="failed",
            note="not installed",
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
        [("SUCOAS", "force_on_high"), ("end_round", None), ("end_round", None)]
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
            ("SUCOAS", "force_on_high"),
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
        max_interventions=3,
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

    client = ScriptedClient(
        [("SUCOAS", "force_on_high"), ("FRD7", "state_distance_check")]
    )
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

    client = ScriptedClient([("FRD7", "essentiality_scan")] * 20)
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
    """The difference between "wrong" and "too much" is worth 8% of the product.

    Watching a run: ``force_on_high`` on the glyoxylate shunt was refused on the growth floor,
    the agent moved to a different reaction, and what ``force_on_low`` on that same reaction
    would have collected was left behind. The rejection now names the gentler move.
    """

    from cmm.jev.actions import GENTLER_ALTERNATIVE

    assert GENTLER_ALTERNATIVE["force_on_high"] == "force_on_low"
    assert GENTLER_ALTERNATIVE["amplify_5x"] == "amplify_2x"

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
        screen_amplifications=False,
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
        screen_amplifications=True,
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
        board[0], product="EX_lac__D_e", growth_floor=0.2, allow_look=False
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

    client = ScriptedClient([("SUCOAS", "force_on_high"), ("end_round", None)] * 4)
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
    )
    result = run_jev_design(config, client=client)

    # The same first move is available in every round, which it would not be if the design
    # from the round before were still standing: an intervened reaction leaves the board.
    firsts = [tick for tick in result.ticks if tick.tick_index == 1]
    assert len(firsts) == 3
    assert {tick.target for tick in firsts} == {"SUCOAS"}
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
