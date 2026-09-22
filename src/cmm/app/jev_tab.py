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
from pathlib import Path
import tempfile

from qtpy.QtCore import QObject, Qt, Signal
from qtpy.QtGui import QColor
from qtpy.QtWidgets import (
    QCheckBox,
    QInputDialog,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QScrollArea,
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

#: How many options of each stage the deliberation panel shows. Enough to see the runner-up
#: and the shape of the tail, short enough to read at a glance while the run is moving.
_WEIGHT_ROWS = 6

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
    "The JEV agent needs an OpenRouter API key. Add one from the JEV menu, or set "
    "OPENROUTER_API_KEY in the environment. Every other tab works without one."
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

        budget_row = QHBoxLayout()
        self.jev_rounds_spin = QSpinBox()
        self.jev_rounds_spin.setRange(1, 50)
        self.jev_rounds_spin.setValue(5)
        self.jev_ticks_spin = QSpinBox()
        # A step is one decision, not one edit: an undo and a scan each cost one. Steps are
        # cheap (about 0.6 s and $0.00016), so the ceiling is generous and the design is
        # bounded separately by the intervention count next to it.
        self.jev_ticks_spin.setRange(1, 1000)
        self.jev_ticks_spin.setValue(20)
        self.jev_targets_spin = QSpinBox()
        self.jev_targets_spin.setRange(1, 20)
        self.jev_targets_spin.setValue(4)
        budget_row.addWidget(QLabel("rounds"))
        budget_row.addWidget(self.jev_rounds_spin)
        budget_row.addWidget(QLabel("steps per round"))
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

        # What the person running this knows and the model does not. Guidance, not
        # permission: it cannot widen the move vocabulary or lift the growth floor, because
        # the agent still answers only with the criteria CMM supplies.
        self.jev_brief = QPlainTextEdit()
        self.jev_brief.setPlaceholderText(
            "What should the agent know that the model does not? One point per line, for "
            "example:\n"
            "- The published targets for succinate in E. coli are ldhA, pflB and ptsG.\n"
            "- Growth has to stay above 0.1 per hour for this strain to be useful.\n"
            "- NADPH supply is the cofactor I expect to be limiting.\n"
            "- Leave the pentose phosphate pathway alone; we cannot engineer it here."
        )
        self.jev_brief.setMaximumHeight(78)
        form.addRow("Brief for the agent:", self.jev_brief)

        web_row = QHBoxLayout()
        self.jev_web_check = QCheckBox(
            "Look up published evidence on the web (slower, about $0.05 a lookup)"
        )
        self.jev_organism = QLineEdit()
        self.jev_organism.setPlaceholderText("organism, e.g. Escherichia coli")
        self.jev_organism.setToolTip(
            "Required for the literature lookup. There is no default: asking the published "
            "record about the wrong species returns an answer that is confident and wrong."
        )
        self.jev_organism.setEnabled(False)
        self.jev_web_check.toggled.connect(self.jev_organism.setEnabled)
        self.jev_web_check.toggled.connect(lambda _: self._update_jev_run_state())
        self.jev_organism.textChanged.connect(lambda _: self._update_jev_run_state())
        web_row.addWidget(self.jev_web_check)
        web_row.addWidget(self.jev_organism, 1)
        form.addRow("", web_row)

        self.jev_run_btn = QPushButton("Let JEV play")
        self.jev_run_btn.clicked.connect(self.run_jev_agent)
        form.addRow("", self.jev_run_btn)
        layout.addWidget(controls)

        # Two determinate bars, because a run has two clocks and they answer different
        # questions: how far through the game, and how far through this round. Both maxima
        # are known before the first call.
        bars = QHBoxLayout()
        bars.setSpacing(10)
        self.jev_round_progress = QProgressBar()
        self.jev_round_progress.setFormat("round \u2014")
        self.jev_round_progress.setRange(0, 1)
        self.jev_round_progress.setValue(0)
        self.jev_progress = QProgressBar()
        self.jev_progress.setFormat("idle")
        self.jev_progress.setRange(0, 1)
        self.jev_progress.setValue(0)
        # Tall enough to read across the room, and two colours so a glance tells which clock
        # is which: the round is the game, the step is the move inside it.
        for bar, chunk in (
            (self.jev_round_progress, "#2f5d8a"),
            (self.jev_progress, "#4c7a34"),
        ):
            bar.setMinimumHeight(26)
            bar.setTextVisible(True)
            bar.setStyleSheet(
                "QProgressBar { border: 1px solid #c2ccd6; border-radius: 4px; "
                "background: #f2f5f8; text-align: center; font-weight: bold; "
                "color: #23313f; } "
                f"QProgressBar::chunk {{ background: {chunk}; border-radius: 3px; }}"
            )
        bars.addWidget(self.jev_round_progress, 2)
        bars.addWidget(self.jev_progress, 3)
        layout.addLayout(bars)

        self.jev_summary = QLabel(
            "JEV is a decision model: it is shown the metabolic state and picks one move at "
            "a time from a fixed set. CMM executes each move, re-solves, and redraws the "
            "flux map — so you watch the design happen."
        )
        self.jev_summary.setWordWrap(True)
        self.jev_summary.setAlignment(Qt.AlignTop)
        # Capped and scrollable. A finished run's summary lists every intervention and every
        # note, and left to grow it took the window's whole height and squeezed the flux map
        # into a ninety-pixel strip — the one thing this tab exists to show.
        summary_area = QScrollArea()
        summary_area.setWidget(self.jev_summary)
        summary_area.setWidgetResizable(True)
        summary_area.setFrameShape(QScrollArea.NoFrame)
        summary_area.setMaximumHeight(120)
        summary_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        layout.addWidget(summary_area)

        results = QSplitter(Qt.Horizontal)
        map_box = QVBoxLayout()
        self.jev_canvas_holder = QVBoxLayout()
        map_box.addLayout(self.jev_canvas_holder, 1)
        # The picture needs to say what it is. Asked what they were looking at, the honest
        # answer was "a curated Escher map, or a schematic when no curated map fits the
        # model", and nothing on screen said either.
        self.jev_map_caption = QLabel("")
        self.jev_map_caption.setWordWrap(True)
        self.jev_map_caption.setStyleSheet("color: #5a6b7c; font-size: 11px;")
        map_box.addWidget(self.jev_map_caption)
        holder = QWidget()
        holder.setLayout(map_box)
        # The map is the point of this tab, so it gets a floor no stretch factor can take
        # away. Stretch factors alone lost: six stretched table columns claim a minimum width
        # between them that squeezed the map to ninety pixels of unreadable smear.
        holder.setMinimumWidth(620)
        results.addWidget(holder)

        self.jev_table = QTableWidget(0, 6)
        self.jev_table.setHorizontalHeaderLabels(
            ["Move", "Reaction", "Action", "Outcome", "Product", "Growth"]
        )
        header = self.jev_table.horizontalHeader()
        # Only the outcome column carries a sentence; the rest are short and fixed, so sizing
        # them to their contents leaves the width for the map instead of spreading it evenly.
        for column in range(6):
            header.setSectionResizeMode(
                column,
                QHeaderView.Stretch if column == 3 else QHeaderView.ResizeToContents,
            )
        self.jev_table.verticalHeader().setVisible(False)
        self.jev_table.setAlternatingRowColors(True)
        self.jev_table.setMinimumWidth(500)

        # What the agent weighed, above what it did. A ``choice`` answer carries a
        # probability for every option it was offered, so each move records a complete
        # ranking of the board — showing it turns "the agent picked ICL" into "it picked ICL
        # over FUM by 0.40 to 0.16", which is the difference between an assertion and
        # evidence, and is often where the interesting runner-up is.
        self.jev_thinking = QLabel("JEV has not been asked anything yet.")
        self.jev_thinking.setWordWrap(True)
        self.jev_thinking.setObjectName("jevthinking")
        self.jev_thinking.setAlignment(Qt.AlignTop)
        # Capped so a long reason cannot eat the weights below it, which are the part that
        # has to stay readable while the run moves.
        self.jev_thinking.setMaximumHeight(96)

        self.jev_weights = QTableWidget(0, 2)
        self.jev_weights.setHorizontalHeaderLabels(["Option", "How strongly"])
        weights_header = self.jev_weights.horizontalHeader()
        weights_header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        weights_header.setSectionResizeMode(1, QHeaderView.Stretch)
        self.jev_weights.verticalHeader().setVisible(False)
        self.jev_weights.setShowGrid(False)
        self.jev_weights.setMinimumHeight(150)

        panel = QWidget()
        panel_layout = QVBoxLayout(panel)
        panel_layout.setContentsMargins(0, 0, 0, 0)
        panel_layout.addWidget(self.jev_thinking)
        panel_layout.addWidget(self.jev_weights, 1)

        right = QSplitter(Qt.Vertical)
        right.addWidget(panel)
        right.addWidget(self.jev_table)
        right.setSizes([260, 200])
        right.setMinimumWidth(500)
        results.addWidget(right)
        results.setStretchFactor(0, 3)
        results.setStretchFactor(1, 2)
        results.setSizes([780, 520])
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
        self.jev_product_combo.blockSignals(True)
        self.jev_product_combo.clear()
        self.jev_product_combo.addItems(exchanges)
        self.jev_product_combo.blockSignals(False)
        if self._default_product and self._default_product in exchanges:
            self.jev_product_combo.setCurrentText(self._default_product)

        from cmm.jev import credentials

        has_key = credentials.key_source() != "none"
        self._jev_has_key = has_key
        self._jev_has_exchanges = bool(exchanges)
        self._update_jev_run_state()
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

    def _update_jev_run_state(self) -> None:
        """Enable the run only when everything it needs is present, and say what is missing."""

        if not hasattr(self, "jev_run_btn"):
            return
        needs_organism = (
            self.jev_web_check.isChecked() and not self.jev_organism.text().strip()
        )
        ready = (
            getattr(self, "_jev_has_exchanges", False)
            and getattr(self, "_jev_has_key", False)
            and not needs_organism
        )
        self.jev_run_btn.setEnabled(ready)
        if needs_organism and ready is False and getattr(self, "_jev_has_key", False):
            self.jev_run_btn.setToolTip(
                "Name the organism: a literature lookup about the wrong species returns an "
                "answer that is confident and wrong."
            )

    # -- the credential -----------------------------------------------------

    def set_jev_api_key(self) -> None:
        """Ask for an OpenRouter key and save it, stating plainly where it will live.

        A secret written to disk in plain text is the user's decision to make, so the dialog
        names the file before they type, rather than after.
        """

        from cmm.jev import credentials

        location = credentials.key_path()
        key, accepted = QInputDialog.getText(
            self,
            "OpenRouter API key",
            "The JEV agent is the only part of CMM that needs a key.\n\n"
            f"It will be saved in plain text at:\n{location}\n"
            "readable only by you, and removable from this menu at any time.\n"
            f"Setting {credentials.ENV_VAR} in the environment overrides it.\n\n"
            "Key:",
            QLineEdit.Password,
        )
        if not accepted:
            return
        try:
            saved = credentials.save_key(key)
        except ValueError:
            QMessageBox.warning(
                self, "OpenRouter API key", "No key was entered; nothing was saved."
            )
            return
        except OSError as error:
            QMessageBox.warning(
                self, "OpenRouter API key", f"The key could not be saved: {error}"
            )
            return
        self._refresh_jev_inputs()
        QMessageBox.information(
            self,
            "OpenRouter API key",
            f"Saved to {saved}.\nThe JEV Agent tab is ready.",
        )

    def clear_jev_api_key(self) -> None:
        """Delete the saved key, and say what is still in force if anything is."""

        from cmm.jev import credentials

        removed = credentials.clear_key()
        self._refresh_jev_inputs()
        if not removed:
            QMessageBox.information(
                self, "OpenRouter API key", "There was no saved key to remove."
            )
            return
        still_set = credentials.key_source() == "environment"
        message = f"Removed {credentials.key_path()}."
        if still_set:
            message += (
                f"\n\n{credentials.ENV_VAR} is still set in this environment, so the agent "
                "can still run. Unset it and restart CMM to stop that too."
            )
        QMessageBox.information(self, "OpenRouter API key", message)

    def show_jev_api_key_status(self) -> None:
        """Where the key in force came from, without showing the key."""

        from cmm.jev import credentials

        source = credentials.key_source()
        if source == "environment":
            import os

            shown = credentials.masked(os.environ[credentials.ENV_VAR])
            text = (
                f"In force: the {credentials.ENV_VAR} environment variable ({shown}).\n"
                "It takes precedence over any saved key."
            )
        elif source == "saved":
            text = (
                f"In force: the key saved at {credentials.key_path()} "
                f"({credentials.masked(credentials.stored_key())})."
            )
        else:
            text = (
                "No key is set. The JEV Agent tab is disabled until one is.\n\n"
                f"Either set {credentials.ENV_VAR} in the environment, or use "
                "JEV \u25b8 Set API Key."
            )
        QMessageBox.information(self, "OpenRouter API key", text)

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
        self.jev_round_progress.setRange(0, self.jev_rounds_spin.value())
        self.jev_round_progress.setValue(0)
        self.jev_round_progress.setFormat(f"round 0 of {self.jev_rounds_spin.value()}")
        self.jev_progress.setRange(0, self.jev_ticks_spin.value())
        self.jev_progress.setValue(0)
        self.jev_progress.setFormat(
            f"step 0 of up to {self.jev_ticks_spin.value()} this round"
        )
        self.jev_weights.setRowCount(0)
        self.jev_thinking.setText("Asking JEV for its first move\u2026")
        self.jev_summary.setText(f"JEV is playing for {html.escape(product)}…")

        # Serialize on the UI thread: the worker must not touch this model's solver object.
        from cobra.io import write_sbml_model

        scratch = Path(tempfile.mkdtemp(prefix="cmm-jev-"))
        model_path = scratch / "model.xml"
        write_sbml_model(self.model, str(model_path))

        config = JevConfig(
            model_path=model_path,
            product=product,
            rounds=self.jev_rounds_spin.value(),
            steps_per_round=self.jev_ticks_spin.value(),
            brief=self.jev_brief.toPlainText(),
            max_interventions=self.jev_targets_spin.value(),
            growth_floor=self.jev_growth_spin.value(),
            candidate_limit=self.jev_board_spin.value(),
            enable_web_research=self.jev_web_check.isChecked(),
            organism=self.jev_organism.text().strip(),
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
            played = len(self._jev_frames)
            self.jev_progress.setFormat(
                f"finished after {played} step{'' if played == 1 else 's'}"
            )
            self.jev_round_progress.setFormat(
                f"{self.jev_round_progress.value()} round"
                f"{'' if self.jev_round_progress.value() == 1 else 's'} played"
            )

        if worker.error is not None:
            raise worker.error
        return worker.result

    def _on_jev_tick(self, tick, fluxes) -> None:
        """One move arrived from the worker. Append it and redraw the map."""

        self._jev_frames.append((tick, fluxes))
        self._append_jev_row(tick)
        self._draw_jev_map(tick, fluxes)
        self._show_jev_thinking(tick)
        self._advance_jev_progress(tick)
        self.status_label.setText(tick.headline())

    def _advance_jev_progress(self, tick) -> None:
        """Move both bars. Neither maximum is a target.

        The agent may end a round on its first step, so a bar is allowed to stop short rather
        than being stretched to look complete.
        """

        best = max(
            (float(frame.product_flux) for frame, _ in self._jev_frames),
            default=0.0,
        )
        self.jev_round_progress.setValue(
            min(tick.round_index, self.jev_round_progress.maximum())
        )
        self.jev_round_progress.setFormat(
            f"round {tick.round_index} of {self.jev_round_progress.maximum()} "
            f"\u00b7 best {best:.4g}"
        )
        self.jev_progress.setValue(min(tick.tick_index, self.jev_progress.maximum()))
        self.jev_progress.setFormat(
            f"step {tick.tick_index} of up to {self.jev_progress.maximum()} this round"
        )

    def _append_jev_row(self, tick) -> None:
        row = self.jev_table.rowCount()
        self.jev_table.insertRow(row)
        # A move that is not about one reaction has no reaction to name. Repeating the move
        # in both columns was wrong as well as wide: adopt_best_design is not a reaction id.
        reaction = "—" if tick.target == tick.action else tick.target
        values = [
            f"R{tick.round_index}S{tick.tick_index}",
            reaction,
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

    def _show_jev_thinking(self, tick) -> None:
        """What the agent weighed on this step, and how strongly.

        Both stages are shown: which reaction, then what to do to it. A ``choice`` answer
        carries a probability for every option it was offered, so this is the whole
        deliberation rather than a report of the winner.
        """

        chosen = tick.action or "no move"
        target = "\u2014" if tick.target == tick.action else tick.target
        headline = (
            f"<b>Round {tick.round_index}, step {tick.tick_index}</b><br>"
            f"chose <b>{html.escape(chosen)}</b>"
            + (f" on <b>{html.escape(target)}</b>" if target != "\u2014" else "")
        )
        if tick.target_confidence is not None:
            headline += f" &nbsp;\u00b7&nbsp; confidence {tick.target_confidence:.2f}"
        if tick.benefit_score is not None:
            headline += f"<br>expected gain {tick.benefit_score:.2f} of 4"
            if tick.predicted_growth_risk is not None:
                risk = f"{tick.predicted_growth_risk:.0%}"
                headline += f" &nbsp;\u00b7&nbsp; it put the growth risk at {risk}"
        headline += f"<br><i>{html.escape(tick.reason)}</i>"
        self.jev_thinking.setText(headline)

        rows: list[tuple[str, str, float, bool]] = [
            ("which target", name, weight, name == tick.target)
            for name, weight in tick.target_ranking[:_WEIGHT_ROWS]
            if weight > 0.004
        ]
        rows += [
            ("what to do", name, weight, name == tick.action)
            for name, weight in tick.action_ranking[:_WEIGHT_ROWS]
            if weight > 0.004
        ]

        self.jev_weights.setRowCount(len(rows))
        for row, (stage, name, weight, is_choice) in enumerate(rows):
            label = QTableWidgetItem(f"{name}   ({stage})")
            label.setFlags(label.flags() & ~Qt.ItemIsEditable)
            if is_choice:
                font = label.font()
                font.setBold(True)
                label.setFont(font)
            self.jev_weights.setItem(row, 0, label)

            bar = QProgressBar()
            bar.setRange(0, 100)
            bar.setValue(int(round(weight * 100)))
            bar.setFormat(f"{weight:.2f}")
            bar.setTextVisible(True)
            bar.setMaximumHeight(18)
            if is_choice:
                bar.setStyleSheet(
                    "QProgressBar::chunk { background: #4c7a34; } QProgressBar { border: "
                    "1px solid #c8d2dc; border-radius: 3px; text-align: center; }"
                )
            self.jev_weights.setCellWidget(row, 1, bar)
        self.jev_weights.resizeRowsToContents()

    def _draw_jev_map(self, tick, fluxes) -> None:
        """Redraw the flux map for this move, on the curated map when the model has one."""

        from cmm.visualization import escher_flux_map, network_flux_map

        # "end_round on end_round" is not a move description. Name the reaction only when
        # the move is about one, which also keeps the title short enough not to be clipped.
        move = tick.action or "no move"
        if tick.target != tick.action:
            move = f"{move} on {tick.target}"
        title = f"R{tick.round_index}S{tick.tick_index}  {move}  —  product {tick.product_flux:.4g}"
        try:
            if self._map_path:
                # Metabolite labels off. At full size the curated map reads well, but
                # this panel scales it to roughly a third of that and the labels become
                # overlapping smudges. Reaction names are what a move is about and they
                # survive; the metabolite names are recoverable with the toolbar's zoom.
                figure = escher_flux_map(
                    self._map_path,
                    dict(fluxes),
                    title="",
                    width=9.0,
                    label_metabolites=False,
                )
            else:
                figure = network_flux_map(self.model, dict(fluxes), title="")
            # Anchored to the figure, left-aligned, rather than centred on the axes. The
            # renderer centres its title on the map's own extent, which reaches past the
            # panel once the canvas scales it down, and the first character is lost.
            figure.suptitle(title, x=0.015, ha="left", fontsize=11, fontweight="bold")
        except Exception as exc:  # a map that cannot be drawn must not stop the game
            self.status_label.setText(f"{tick.headline()} (map not drawn: {exc})")
            return
        self._set_figure(self.jev_canvas_holder, "jev", figure)
        self.jev_map_caption.setText(
            "Curated Escher map of this model, coloured and widened by flux. Metabolite "
            "labels are hidden at this size \u2014 use the magnifier to zoom, or the Flux "
            "Map tab for the full-size figure."
            if self._map_path
            else "Schematic of the highest-flux reactions; no curated Escher map fits this "
            "model. Arrow colour and width both scale with the flux magnitude."
        )

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
