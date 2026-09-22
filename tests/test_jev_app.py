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
    assert window._tab_index("JEV Agent") == window.tabs.count() - 1


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
    assert window.jev_substrate_combo.currentText() == "EX_glc__D_e"
    assert window.jev_run_btn.isEnabled()


def test_the_map_is_redrawn_once_per_move_while_the_run_is_still_going(
    window, monkeypatch
) -> None:
    """The reason the tab uses a queued signal instead of a modal progress dialog."""

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("JEV Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)

    client = ScriptedClient([("SUCOAS", "force_on_high")])
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient", lambda **kwargs: client, raising=True
    )

    redraws: list[str] = []
    original = window._draw_jev_map
    monkeypatch.setattr(
        window,
        "_draw_jev_map",
        lambda tick, fluxes: (redraws.append(tick.target), original(tick, fluxes)),
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
    window._goto_tab("JEV Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("SUCOAS", "force_on_high")]),
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
    window._goto_tab("JEV Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(2)
    window.jev_ticks_spin.setValue(3)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("SUCOAS", "force_on_high")]),
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
    window._goto_tab("JEV Agent")
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

    The design stays small — that budget is the spin box next to it, because the number a
    laboratory has to build is a different quantity from the number of moves played.
    """

    assert window.jev_ticks_spin.maximum() >= 1000
    assert window.jev_targets_spin.maximum() <= 20


def test_the_flux_map_keeps_its_space(window, monkeypatch) -> None:
    """The map is the point of this tab, so nothing else may squeeze it out.

    Both directions went wrong once. Six stretched table columns claimed enough width between
    them to leave the map ninety pixels across, and a finished run's summary grew until it
    took the window's height and left the map a strip.
    """

    monkeypatch.setenv("OPENROUTER_API_KEY", "sk-or-test")
    window._refresh_jev_inputs()
    window._goto_tab("JEV Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(2)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("SUCOAS", "force_on_high")]),
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
    window._goto_tab("JEV Agent")
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
    window._goto_tab("JEV Agent")
    window.jev_product_combo.setCurrentText("EX_succ_e")
    window.jev_rounds_spin.setValue(1)
    window.jev_ticks_spin.setValue(1)
    monkeypatch.setattr(
        "cmm.jev.engine.JevClient",
        lambda **kwargs: ScriptedClient([("SUCOAS", "force_on_high")]),
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
    assert "JEV menu" in window.jev_summary.text()

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
        if action.text().replace("&", "") == "JEV"
    )
    labels = [action.text() for action in jev_menu.actions() if action.text()]
    assert "Set API Key…" in labels
    assert "Clear API Key" in labels
    assert "Where is my key?" in labels
