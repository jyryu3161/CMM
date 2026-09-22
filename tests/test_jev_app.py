"""The JEV Agent tab, driven headlessly with a scripted agent.

The point of the tab is that the flux map moves *while* the agent is still playing, so the
test that matters is the one counting redraws during the run rather than after it.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest

pytest.importorskip("qtpy")

from cobra.io import load_model  # noqa: E402

from cmm.core.condition import Condition, ReactionBound  # noqa: E402
from test_jev import ScriptedClient  # noqa: E402

ANAEROBIC = Condition(
    name="glucose_anaerobic",
    bounds=(
        ReactionBound("EX_glc__D_e", -10.0, 1000.0),
        ReactionBound("EX_o2_e", 0.0, 1000.0),
    ),
)


@pytest.fixture(scope="module")
def app():
    from qtpy.QtWidgets import QApplication

    instance = QApplication.instance()
    if instance is None:
        try:
            instance = QApplication([])
        except Exception:  # pragma: no cover - no Qt platform available
            pytest.skip("no Qt platform available")
    return instance


@pytest.fixture
def window(app):
    from cmm.app.main_window import CmmMainWindow

    model = load_model("textbook")
    ANAEROBIC.apply_to(model)
    return CmmMainWindow(model)


def test_the_tab_exists_and_sits_last(window) -> None:
    assert window._tab_index("Agent") == window.tabs.count() - 1


def test_without_an_api_key_the_tab_explains_itself_and_nothing_else_breaks(
    window, monkeypatch
) -> None:
    """The whole feature degrades to one disabled button and a sentence."""

    monkeypatch.delenv("OPENROUTER_API_KEY", raising=False)
    window._refresh_jev_inputs()

    assert not window.jev_run_btn.isEnabled()
    assert "OPENROUTER_API_KEY" in window.jev_summary.text()
    # Every other tab is untouched by the missing credential.
    assert window._tab_index("Simulation") is not None
    assert window.tabs.count() == 11


def test_the_product_combo_is_filled_from_the_model(window, monkeypatch) -> None:
    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()

    products = [
        window.jev_product_combo.itemText(i)
        for i in range(window.jev_product_combo.count())
    ]
    assert "EX_succ_e" in products
    assert window.jev_run_btn.isEnabled()
    # There is no substrate to choose: the yield is quoted per whatever carbon source the
    # condition actually feeds the model, which the wild-type solve already says.
    assert not hasattr(window, "jev_substrate_combo")


def test_the_map_is_redrawn_once_per_move_while_the_run_is_still_going(
    window, monkeypatch
) -> None:
    """The reason the tab uses a queued signal instead of a modal progress dialog."""

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)

    client = ScriptedClient([("PFL", "knockout")])
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient", lambda **kwargs: client, raising=True
    )

    redraws: list[str] = []
    original = window._draw_jev_map
    monkeypatch.setattr(
        window,
        "_draw_jev_map",
        lambda tick, fluxes, previous=None: (
            redraws.append(tick.target),
            original(tick, fluxes, previous),
        ),
    )

    window.run_jev_agent()

    assert redraws, "the map must be redrawn during the run, not only at the end"
    assert len(redraws) == window.jev_table.rowCount()
    assert window._active_figure() is not None
    assert window._active_table() is window.jev_table
    assert "EX_succ_e" in window.jev_summary.text()
    # The tab always says a run is a hypothesis, whatever the numbers were.
    assert "hypothesis" in window.jev_summary.text()
    # Controls come back regardless of how the run ended.
    assert window.jev_run_btn.isEnabled()
    assert window.jev_controls.isEnabled()


def test_the_progress_and_decision_figures_render_after_a_run(
    window, monkeypatch
) -> None:
    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("PFL", "knockout")]),
        raising=True,
    )

    window.run_jev_agent()
    window.show_jev_progress()
    assert window._active_figure() is not None
    window.show_jev_last_decision()
    assert window._active_figure() is not None


def test_asking_for_a_figure_before_a_run_says_so_rather_than_failing(window) -> None:
    window.show_jev_progress()
    assert "Run the agent first" in window.jev_summary.text()
    window.show_jev_last_decision()
    assert "Run the agent first" in window.jev_summary.text()


def test_the_progress_bars_count_rounds_and_steps_separately(
    window, monkeypatch
) -> None:
    """Two clocks, because a run has two and they answer different questions.

    Neither maximum is a target: the agent may end a round on its first step, so a bar is
    allowed to stop short rather than being stretched to look complete.
    """

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(2)
    window.jev_ticks_spin.setValue(3)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("PFL", "knockout")]),
        raising=True,
    )

    window.run_jev_agent()

    assert window.jev_round_progress.maximum() == 2  # rounds
    assert window.jev_progress.maximum() == 3  # steps within one round
    assert window.jev_progress.value() <= window.jev_progress.maximum()
    assert window.jev_round_progress.value() <= window.jev_round_progress.maximum()
    assert "finished after" in window.jev_progress.format()
    assert "played" in window.jev_round_progress.format()


def test_the_brief_box_reaches_the_run(window, monkeypatch) -> None:
    """What the person types is what the agent is shown."""

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(1)
    window.jev_brief.setPlainText(
        "- Known targets for succinate are ldhA and pflB.\n- NADPH is the limiting cofactor."
    )

    client = ScriptedClient([("end_round", None)])
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient", lambda **kwargs: client, raising=True
    )
    window.run_jev_agent()

    assert client.states, "the run must have reached the agent"
    assert "ldhA" in client.states[0]["your_brief"]["text"]


def test_the_step_ceiling_allows_a_long_game(window) -> None:
    """Steps are decisions, not edits: an undo and a scan each cost one, so the game is long.

    The design stays small — that budget is the pair of spin boxes below it, because the
    number a laboratory has to build is a different quantity from the number of moves played,
    and a deletion is a different quantity from a promoter swap.
    """

    assert window.jev_ticks_spin.maximum() >= 1000
    assert window.jev_knockouts_spin.maximum() <= 30
    assert window.jev_knockdowns_spin.maximum() <= 30
    # The seeded strain design takes three deletions on its own, so the default has to leave
    # the agent something to play with after adopting it.
    assert window.jev_knockouts_spin.value() > 3


def test_the_flux_map_keeps_its_space(window, monkeypatch) -> None:
    """The map is the point of this tab, so nothing else may squeeze it out.

    Both directions went wrong once. Six stretched table columns claimed enough width between
    them to leave the map ninety pixels across, and a finished run's summary grew until it
    took the window's height and left the map a strip.
    """

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("PFL", "knockout")]),
        raising=True,
    )

    window.resize(1500, 950)
    window.run_jev_agent()

    canvas = window._jev_canvas
    assert canvas is not None
    assert canvas.width() >= 400, "the flux map must not be squeezed to a strip"
    assert canvas.height() >= 200


def test_a_move_that_is_not_about_a_reaction_does_not_name_one(
    window, monkeypatch
) -> None:
    """``adopt_best_design`` is a move, not a reaction id, in the table and in the title."""

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(1)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("end_round", None)]),
        raising=True,
    )

    window.run_jev_agent()

    assert window.jev_table.rowCount() == 1
    assert window.jev_table.item(0, 1).text() == "—"  # no reaction
    assert window.jev_table.item(0, 2).text() == "end_round"
    assert window.jev_table.item(0, 0).text().startswith("R1S")  # steps, not ticks


def test_the_dashboard_shows_what_the_agent_weighed(window, monkeypatch) -> None:
    """A choice answer carries a probability for every option, so show the whole ranking.

    "The agent picked ICL" is an assertion; "it picked ICL over FUM by 0.40 to 0.16" is
    evidence, and the runner-up is often the interesting row.
    """

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(1)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("PFL", "knockout")]),
        raising=True,
    )

    window.run_jev_agent()

    assert window.jev_weights.rowCount() > 0, "the deliberation must be shown"
    labels = [
        window.jev_weights.item(row, 0).text()
        for row in range(window.jev_weights.rowCount())
    ]
    assert any("which target" in label for label in labels)
    assert any("what to do" in label for label in labels)
    # Every weight row carries a bar, not just a number.
    assert window.jev_weights.cellWidget(0, 1) is not None
    text = window.jev_thinking.text()
    assert "chose" in text and "confidence" in text


def test_the_key_can_be_saved_and_cleared(window, monkeypatch, tmp_path) -> None:
    """The key is the user's to keep or remove, and the tab follows it.

    The menu handlers around this are three lines each and end in a modal box, which a test
    run cannot answer; what is worth asserting is the credential logic and the tab's
    response to it.
    """

    from cmm.jev import credentials

    monkeypatch.delenv("OPENROUTER_API_KEY", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

    window._refresh_jev_inputs()
    assert not window.jev_run_btn.isEnabled()
    assert credentials.key_source() == "none"
    assert "Agent menu" in window.jev_summary.text()

    credentials.save_key("sk-or-v1-testkey0000")
    window._refresh_jev_inputs()
    assert window.jev_run_btn.isEnabled()
    assert credentials.key_source() == "saved"
    # Readable only by this user, and never shown in full.
    assert credentials.key_path().stat().st_mode & 0o077 == 0
    shown = credentials.masked(credentials.stored_key())
    assert shown.endswith("0000")
    assert "testkey" not in shown

    # The environment always wins, so a key exported for one session is never overridden.
    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-v1-fromenv")
    assert credentials.key_source() == "environment"
    monkeypatch.delenv("OPENROUTER_API_KEY")

    assert credentials.clear_key() is True
    assert credentials.clear_key() is False  # nothing left to remove
    window._refresh_jev_inputs()
    assert credentials.key_source() == "none"
    assert not window.jev_run_btn.isEnabled()


def test_the_key_menu_items_exist(window) -> None:
    """The three things a user needs to do with a credential are one click each."""

    jev_menu = next(
        action.menu()
        for action in window.menuBar().actions()
        if action.text().replace("&", "") == "Agent"
    )
    labels = [action.text() for action in jev_menu.actions() if action.text()]
    assert "Set OpenRouter API Key…" in labels
    assert "Clear OpenRouter API Key" in labels
    assert "Where is my key?" in labels


def test_the_run_reports_what_it_has_spent(window, monkeypatch) -> None:
    """A decision is about $0.00016, which is exactly why the number is worth showing.

    Left off the screen, the only available estimate is whatever order of magnitude a person
    assumes an agent loop costs, and that assumption is wrong by three or four of them.
    """

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("PFL", "knockout")]),
        raising=True,
    )

    assert window.jev_cost.text() == "$0.0000"
    window.run_jev_agent()

    text = window.jev_cost.text()
    assert text.startswith("$")
    assert "decisions" in text


def test_a_step_that_changed_no_flux_says_so(window) -> None:
    """The answer to "why does the map not move?".

    It does not move because it did not change: a scan changes nothing by definition, a move
    that breached the growth floor was reverted before the frame was taken, and deleting a
    reaction already carrying nothing is a real move with no immediate consequence. Only a
    move that stuck and mattered redraws differently. An unchanged picture that does not
    admit it reads as a broken redraw.
    """

    same = {"A": 1.0, "B": -2.0}
    assert "No flux changed" in window._flux_change(same, dict(same))
    assert "Wild-type flux distribution" in window._flux_change(same, None)

    moved = window._flux_change({"A": 4.0, "B": -2.0}, same)
    assert "1 reaction changed flux" in moved
    assert "A 1→4" in moved


def test_stopping_is_asked_for_once_and_keeps_the_run(window, monkeypatch) -> None:
    """Not a kill: the worker is inside a solver call for most of its life."""

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")

    assert not window.jev_stop_btn.isEnabled(), "nothing to stop before a run"
    window._jev_stop_requested = False
    window.jev_stop_btn.setEnabled(True)

    window.stop_jev_agent()
    assert window._jev_stop_requested is True
    assert not window.jev_stop_btn.isEnabled()
    assert "Stopping" in window.jev_stop_btn.text()

    # Asking twice is not an error and does not undo anything.
    window.stop_jev_agent()
    assert window._jev_stop_requested is True


def test_a_played_move_can_be_put_back_on_the_map(window, monkeypatch) -> None:
    """Clicking a row in the move log redraws that move's flux distribution."""

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(3)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("PFL", "knockout")]),
        raising=True,
    )
    window.run_jev_agent()

    assert window.jev_table.rowCount() >= 2
    drawn: list[object] = []
    original = window._draw_jev_map
    monkeypatch.setattr(
        window,
        "_draw_jev_map",
        lambda tick, fluxes, previous=None: (
            drawn.append(tick.tick_index),
            original(tick, fluxes, previous),
        ),
    )
    window.jev_table.selectRow(0)
    assert drawn == [window._jev_frames[0][0].tick_index]


def test_the_key_is_reachable_without_finding_the_menu(
    window, monkeypatch, tmp_path
) -> None:
    """A credential that lives only in a menu is a credential nobody finds.

    It is also a prerequisite for the run button working at all, so it belongs beside the
    run button and it has to say whether one is in force.
    """

    from cmm.jev import credentials

    monkeypatch.delenv("OPENROUTER_API_KEY", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    window._refresh_jev_inputs()

    assert window.jev_key_btn.text() == "Set key…"
    assert "no key" in window.jev_key_label.text()
    assert not window.jev_run_btn.isEnabled()
    # And the reason is on screen, not in a tooltip.
    assert not window.jev_hint.isHidden()
    assert "Set key" in window.jev_hint.text()

    credentials.save_key("sk-or-v1-testkey0000")
    window._refresh_jev_inputs()
    assert window.jev_key_btn.text() == "Change key…"
    assert window.jev_key_label.text().endswith("0000)")
    assert window.jev_run_btn.isEnabled()
    assert window.jev_hint.isHidden()
    credentials.clear_key()


def test_the_web_lookup_says_what_it_needs_instead_of_going_quiet(
    window, monkeypatch
) -> None:
    """Ticking the box greyed the run button out with no visible reason.

    The requirement is real — a literature lookup aimed at the wrong species returns an
    answer that is confident and wrong — but a disabled button that does not say why reads
    as a broken button rather than as a question.
    """

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("Agent")
    assert window.jev_run_btn.isEnabled()
    assert window.jev_hint.isHidden()

    window.jev_web_check.setChecked(True)
    assert not window.jev_run_btn.isEnabled()
    assert not window.jev_hint.isHidden()
    assert "organism" in window.jev_hint.text()
    assert "Untick" in window.jev_hint.text(), "the way out has to be named too"

    # One field, and it runs.
    window.jev_organism.setText("Escherichia coli")
    assert window.jev_run_btn.isEnabled()
    assert window.jev_hint.isHidden()

    # Or untick it, and it runs without the lookup.
    window.jev_organism.setText("")
    assert not window.jev_run_btn.isEnabled()
    window.jev_web_check.setChecked(False)
    assert window.jev_run_btn.isEnabled()
