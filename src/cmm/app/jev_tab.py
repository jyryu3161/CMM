"""The JEV Agent tab: watch a decision model play the model, one move at a time.

Kept in its own module and mixed into :class:`~cmm.app.main_window.CmmMainWindow` so the
whole feature costs the existing 4,700-line window four lines. Nothing here runs unless the
tab is opened, and nothing here changes how any other tab behaves.

**Why the live redraw needs a signal.** The engine runs on a worker thread and calls back
once per move from it. Touching a widget or a matplotlib figure from a worker thread is not
allowed, so the callback does the one thing it may do — emit a Qt signal. The connection is
queued, so the slot runs on the UI thread, and the nested event loop the run sits in delivers
those queued slots *during* the run rather than after it. That is what makes the flux map move
while the agent is still playing.

The model is serialized to a temporary SBML file before the worker starts, for the same
reason :meth:`~cmm.app.main_window.CmmMainWindow._run_model_in_background` rebuilds one: a
``cobra.Model`` owns a solver object that belongs to the thread that made it.
"""

from __future__ import annotations

import html
import os
from pathlib import Path
import tempfile

from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtGui import QColor
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

JEV_TAB_NAME = "JEV Agent"

#: Row colours, matching the progress figure's outcome palette.
_OUTCOME_COLOUR = {
    "applied": QColor("#dce9d6"),
    "undone_by_agent": QColor("#fdf3d8"),
    "reverted_infeasible": QColor("#f8dcd6"),
    "reverted_growth_floor": QColor("#f8dcd6"),
    "not_applicable": QColor("#eeeeee"),
    "scan": QColor("#e3edf5"),
    "end_round": QColor("#eeeeee"),
}

_NO_KEY_MESSAGE = (
    "The JEV agent needs an OpenRouter API key. Set OPENROUTER_API_KEY in the environment "
    "and restart CMM. Every other tab works without one."
)


class _TickBridge(QObject):
    """Carries one played move from the worker thread to the UI thread."""

    played = Signal(object, object)


class JevTabMixin:
    """The JEV Agent tab. Mixed into the main window; holds no state of its own until built."""

    # -- construction -------------------------------------------------------

    def _build_jev_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        self._jev_result = None
        self._jev_frames: list[tuple[object, dict]] = []
        self._jev_canvas = None
        self._jev_toolbar = None

        controls = QGroupBox("JEV agent — a decision model plays this model")
        # Held so the whole control block can be disabled while a run is in flight.
        self.jev_controls = controls
        form = QFormLayout(controls)

        self.jev_product_combo = QComboBox()
        self.jev_product_combo.setToolTip(
            "The exchange reaction whose flux the agent is trying to raise."
        )
        form.addRow("Product exchange:", self.jev_product_combo)

        self.jev_substrate_combo = QComboBox()
        self.jev_substrate_combo.setToolTip(
            "Used only for the theoretical yield shown on the agent's scoreboard."
        )
        form.addRow("Substrate exchange:", self.jev_substrate_combo)

        budget_row = QHBoxLayout()
        self.jev_rounds_spin = QSpinBox()
        self.jev_rounds_spin.setRange(1, 50)
        self.jev_rounds_spin.setValue(5)
        self.jev_ticks_spin = QSpinBox()
        self.jev_ticks_spin.setRange(1, 40)
        self.jev_ticks_spin.setValue(6)
        self.jev_targets_spin = QSpinBox()
        self.jev_targets_spin.setRange(1, 20)
        self.jev_targets_spin.setValue(4)
        budget_row.addWidget(QLabel("rounds"))
        budget_row.addWidget(self.jev_rounds_spin)
        budget_row.addWidget(QLabel("moves per round"))
        budget_row.addWidget(self.jev_ticks_spin)
        budget_row.addWidget(QLabel("max interventions"))
        budget_row.addWidget(self.jev_targets_spin)
        budget_row.addStretch(1)
        form.addRow("Budget:", budget_row)

        limits_row = QHBoxLayout()
        self.jev_growth_spin = QDoubleSpinBox()
        self.jev_growth_spin.setDecimals(4)
        self.jev_growth_spin.setRange(0.0, 10.0)
        self.jev_growth_spin.setSingleStep(0.01)
        self.jev_growth_spin.setValue(0.05)
        self.jev_growth_spin.setToolTip(
            "A move that pushes growth below this is reverted by CMM, whatever the agent "
            "predicted."
        )
        self.jev_board_spin = QSpinBox()
        self.jev_board_spin.setRange(4, 60)
        self.jev_board_spin.setValue(24)
        self.jev_board_spin.setToolTip(
            "How many reactions the agent is shown each move. The decision model has a 32K "
            "context, so this is what keeps the state inside it."
        )
        limits_row.addWidget(QLabel("growth floor (1/h)"))
        limits_row.addWidget(self.jev_growth_spin)
        limits_row.addWidget(QLabel("reactions on the board"))
        limits_row.addWidget(self.jev_board_spin)
        limits_row.addStretch(1)
        form.addRow("Rules:", limits_row)

        self.jev_web_check = QCheckBox(
            "Look up published evidence on the web for each candidate (slower, costs more)"
        )
        form.addRow("", self.jev_web_check)

        self.jev_run_btn = QPushButton("Let JEV play")
        self.jev_run_btn.clicked.connect(self.run_jev_agent)
        form.addRow("", self.jev_run_btn)
        layout.addWidget(controls)

        self.jev_summary = QLabel(
            "JEV is a decision model: it is shown the metabolic state and picks one move at "
            "a time from a fixed set. CMM executes each move, re-solves, and redraws the "
            "flux map — so you watch the design happen."
        )
        self.jev_summary.setWordWrap(True)
        layout.addWidget(self.jev_summary)

        results = QSplitter(Qt.Horizontal)
        self.jev_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.jev_canvas_holder)
        results.addWidget(holder)

        self.jev_table = QTableWidget(0, 6)
        self.jev_table.setHorizontalHeaderLabels(
            ["Move", "Reaction", "Action", "Outcome", "Product", "Growth"]
        )
        self.jev_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.jev_table.verticalHeader().setVisible(False)
        self.jev_table.setAlternatingRowColors(True)
        results.addWidget(self.jev_table)
        results.setStretchFactor(0, 5)
        results.setStretchFactor(1, 4)
        layout.addWidget(results, 1)

        self._jev_bridge = _TickBridge()
        # Queued so the slot runs on the UI thread even though the engine emits from the
        # worker. This is what lets the map redraw mid-run instead of once at the end.
        self._jev_bridge.played.connect(self._on_jev_tick, Qt.QueuedConnection)

        self._refresh_jev_inputs()
        return tab

    # -- state --------------------------------------------------------------

    def _refresh_jev_inputs(self) -> None:
        """Fill the exchange combos from the loaded model and gate the Run button."""

        if not hasattr(self, "jev_product_combo"):
            return
        exchanges = sorted(r.id for r in self.model.exchanges)
        for combo in (self.jev_product_combo, self.jev_substrate_combo):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(exchanges)
            combo.blockSignals(False)
        if self._default_product and self._default_product in exchanges:
            self.jev_product_combo.setCurrentText(self._default_product)
        for guess in ("EX_glc__D_e", "EX_glc_e"):
            if guess in exchanges:
                self.jev_substrate_combo.setCurrentText(guess)
                break

        has_key = bool(os.environ.get("OPENROUTER_API_KEY", "").strip())
        ready = bool(exchanges) and has_key
        self.jev_run_btn.setEnabled(ready)
        if not exchanges:
            self.jev_summary.setText(
                "This model has no exchange reactions, so there is no product flux to raise."
            )
            self.jev_run_btn.setToolTip("The model has no exchange reactions.")
        elif not has_key:
            self.jev_summary.setText(_NO_KEY_MESSAGE)
            self.jev_run_btn.setToolTip(_NO_KEY_MESSAGE)
        else:
            self.jev_run_btn.setToolTip("")

        self.jev_table.setRowCount(0)
        self._jev_frames = []
        self._jev_result = None

    # -- the run ------------------------------------------------------------

    def run_jev_agent(self) -> None:
        """Play a full run, redrawing the flux map after every move."""

        from cmm.jev import JevConfig, run_jev_design

        product = self.jev_product_combo.currentText()
        if not product:
            self.jev_summary.setText("Choose a product exchange reaction first.")
            return

        self.jev_table.setRowCount(0)
        self._jev_frames = []
        self.jev_summary.setText(f"JEV is playing for {html.escape(product)}…")

        # Serialize on the UI thread: the worker must not touch this model's solver object.
        from cobra.io import write_sbml_model

        scratch = Path(tempfile.mkdtemp(prefix="cmm-jev-"))
        model_path = scratch / "model.xml"
        write_sbml_model(self.model, str(model_path))

        config = JevConfig(
            model_path=model_path,
            product=product,
            substrate=self.jev_substrate_combo.currentText() or None,
            rounds=self.jev_rounds_spin.value(),
            ticks_per_round=self.jev_ticks_spin.value(),
            max_interventions=self.jev_targets_spin.value(),
            growth_floor=self.jev_growth_spin.value(),
            candidate_limit=self.jev_board_spin.value(),
            enable_web_research=self.jev_web_check.isChecked(),
        )
        bridge = self._jev_bridge

        def _compute():
            return run_jev_design(
                config,
                on_tick=lambda tick, fluxes: bridge.played.emit(tick, dict(fluxes)),
            )

        try:
            result = self._run_jev_visibly(_compute)
        except Exception as exc:
            self.jev_summary.setText(f"The JEV run failed: {html.escape(str(exc))}")
            self.status_label.setText("JEV run failed.")
            return

        self._jev_result = result
        self._show_jev_summary(result)

    def _run_jev_visibly(self, compute):
        """Run ``compute`` off the UI thread with the window left visible.

        The other tabs use ``_run_in_background``, which raises a modal progress dialog over
        the window. That is right for a solve whose only output is a final table, and wrong
        here: the dialog would cover the flux map this tab exists to show moving.

        Dropping the dialog is safe for this run specifically, because the engine loads its
        own model from the temporary SBML written above and never touches ``self.model``.
        The only thing the worker shares with the window is the signal it emits. The controls
        are disabled for the duration so a second run cannot be started on top of the first.

        The nested event loop keeps the call synchronous for the caller — and, more to the
        point, is what delivers the queued per-move signals while the worker is still running.
        """

        from qtpy.QtCore import QEventLoop, QThread, QTimer

        from cmm.app.main_window import _AnalysisWorker

        worker = _AnalysisWorker(compute)
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)

        loop = QEventLoop()
        worker.finished.connect(loop.quit)

        self.jev_run_btn.setEnabled(False)
        self.jev_controls.setEnabled(False)
        QTimer.singleShot(0, thread.start)
        try:
            loop.exec_()
        finally:
            thread.quit()
            thread.wait()
            worker.deleteLater()
            self.jev_controls.setEnabled(True)
            self.jev_run_btn.setEnabled(True)

        if worker.error is not None:
            raise worker.error
        return worker.result

    def _on_jev_tick(self, tick, fluxes) -> None:
        """One move arrived from the worker. Append it and redraw the map."""

        self._jev_frames.append((tick, fluxes))
        self._append_jev_row(tick)
        self._draw_jev_map(tick, fluxes)
        self.status_label.setText(tick.headline())

    def _append_jev_row(self, tick) -> None:
        row = self.jev_table.rowCount()
        self.jev_table.insertRow(row)
        values = [
            f"R{tick.round_index}T{tick.tick_index}",
            tick.target,
            tick.action or "—",
            tick.outcome.replace("_", " "),
            f"{tick.product_flux:.4g}",
            f"{tick.growth:.4g}",
        ]
        colour = _OUTCOME_COLOUR.get(tick.outcome)
        for column, value in enumerate(values):
            item = QTableWidgetItem(value)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            if colour is not None:
                item.setBackground(colour)
            self.jev_table.setItem(row, column, item)
        self.jev_table.scrollToBottom()

    def _draw_jev_map(self, tick, fluxes) -> None:
        """Redraw the flux map for this move, on the curated map when the model has one."""

        from cmm.visualization import escher_flux_map, network_flux_map

        title = (
            f"R{tick.round_index}T{tick.tick_index}  {tick.action or 'no move'} "
            f"on {tick.target}  —  product {tick.product_flux:.4g}"
        )
        try:
            if self._map_path:
                figure = escher_flux_map(self._map_path, dict(fluxes), title=title)
            else:
                figure = network_flux_map(self.model, dict(fluxes), title=title)
        except Exception as exc:  # a map that cannot be drawn must not stop the game
            self.status_label.setText(f"{tick.headline()} (map not drawn: {exc})")
            return
        self._set_figure(self.jev_canvas_holder, "jev", figure)

    def _show_jev_summary(self, result) -> None:
        summary = result.summary()
        product = html.escape(str(summary["product"]))
        wild = summary["wild_type_product_flux"]
        best = summary["best_product_flux"]
        usage = summary["usage"]

        if summary["beat_wild_type"]:
            fold = summary["fold_improvement"]
            gain = (
                f" ({fold:.1f}× wild type)"
                if isinstance(fold, (int, float))
                else " (wild type made none)"
            )
            headline = (
                f"<b>{product} rose from {wild:.4g} to {best:.4g} "
                f"mmol gDW⁻¹ h⁻¹{gain}</b>, at a growth rate of "
                f"{summary['best_growth']:.4g} h⁻¹."
            )
        else:
            headline = (
                f"<b>The agent did not beat the wild type.</b> {product} stayed at "
                f"{wild:.4g} mmol gDW⁻¹ h⁻¹. That is a real result, not "
                "a failed run."
            )

        design = (
            "".join(f"<li>{html.escape(line)}</li>" for line in summary["best_design"])
            or "<li>no interventions survived the rules</li>"
        )
        notes = "".join(
            f"<li>{html.escape(str(note))}</li>" for note in summary["notes"]
        )
        note_block = f"<br>Notes:<ul>{notes}</ul>" if notes else ""

        self.jev_summary.setText(
            f"{headline}<br>Design:<ul>{design}</ul>"
            f"{summary['n_ticks']} moves over {summary['n_rounds']} rounds; "
            f"{usage['calls']} agent decisions costing ${usage['cost_usd']:.4f}."
            f"{note_block}"
            "<br><i>This is a computational hypothesis. The agent's choices are not "
            "guaranteed to repeat on a re-run; the CMM solves behind them are.</i>"
        )
        self.status_label.setText(
            f"JEV run complete: {summary['n_ticks']} moves, best {product} "
            f"{best:.4g} mmol gDW-1 h-1."
        )

    def show_jev_progress(self) -> None:
        """Swap the map for the run's progress curve."""

        if self._jev_result is None:
            self.jev_summary.setText(
                "Run the agent first; there is no progress to plot."
            )
            return
        from cmm.visualization import jev_progress_figure

        self._set_figure(
            self.jev_canvas_holder, "jev", jev_progress_figure(self._jev_result)
        )
        self.status_label.setText("Showing the JEV run's progress.")

    def show_jev_last_decision(self) -> None:
        """Swap the map for the probability the agent put on every option of the last move."""

        if not self._jev_frames:
            self.jev_summary.setText(
                "Run the agent first; there is no decision to show."
            )
            return
        from cmm.visualization import jev_decision_figure

        tick, _ = self._jev_frames[-1]
        self._set_figure(self.jev_canvas_holder, "jev", jev_decision_figure(tick))
        self.status_label.setText("Showing the agent's ranking for the last move.")
