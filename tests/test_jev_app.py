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


def test_the_progress_bar_counts_moves_against_the_budget(window, monkeypatch) -> None:
    """A determinate bar: the move budget is known before the first call.

    Its maximum is an upper bound, not a target — the agent may end a round early — so the
    bar is allowed to stop short rather than being stretched to look complete.
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

    assert window.jev_progress.maximum() == 6  # 2 rounds x 3 moves
    assert window.jev_progress.value() == len(window._jev_frames)
    assert window.jev_progress.value() <= window.jev_progress.maximum()
    assert "finished after" in window.jev_progress.format()


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
