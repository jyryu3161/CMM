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
    names = [action.name for action in applicable_actions(8.2)]
    assert names == [
        "knockout",
        "knockdown_50",
        "knockdown_25",
        "amplify_2x",
        "amplify_5x",
    ]


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
        ticks_per_round=3,
        max_interventions=3,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=2,
        growth_floor=0.05,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=4,
        max_interventions=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=5,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=2,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ticks_per_round=1,
        growth_floor=0.01,
        candidate_limit=12,
        run_moma=False,
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
        ({"ticks_per_round": 0}, "ticks_per_round"),
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
        ticks_per_round=1,
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
