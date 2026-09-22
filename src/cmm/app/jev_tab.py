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

from datetime import datetime
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
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

#: The tab's name. "Agent" rather than "JEV": the decision model behind it is JEV today and
#: the provenance of every run records that, but the tab is the place where an agent plays
#: this model and the name should not have to change when the model does.
JEV_TAB_NAME = "Agent"

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
    "The agent needs an OpenRouter API key. Use the “Set key…” button beside the run "
    "button, or the Agent menu, or set OPENROUTER_API_KEY in the environment. Every other "
    "tab works without one."
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
        self._jev_cost = 0.0
        self._jev_stop_requested = False
        self._jev_running = False
        self._jev_run_dir = None

        controls = QGroupBox(
            "Agent — a decision model plays this model, one move at a time"
        )
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
        self.jev_ticks_spin.setValue(40)
        budget_row.addWidget(QLabel("rounds"))
        budget_row.addWidget(self.jev_rounds_spin)
        budget_row.addWidget(QLabel("steps per round"))
        budget_row.addWidget(self.jev_ticks_spin)
        budget_row.addStretch(1)
        form.addRow("Budget:", budget_row)

        # Counted separately because they are different things to build, and because one
        # shared cap starved the run: the seeded OptKnock design takes three deletions on its
        # own, so a total of four left the agent a single edit and every round ended a step
        # or two after adopting it.
        targets_row = QHBoxLayout()
        self.jev_knockouts_spin = QSpinBox()
        self.jev_knockouts_spin.setRange(0, 30)
        self.jev_knockouts_spin.setValue(6)
        self.jev_knockouts_spin.setToolTip(
            "How many genes the design may delete. The seeded OptKnock design uses three of "
            "these on its own, so leave room above that or the agent inherits a full design."
        )
        self.jev_knockdowns_spin = QSpinBox()
        self.jev_knockdowns_spin.setRange(0, 30)
        self.jev_knockdowns_spin.setValue(3)
        self.jev_knockdowns_spin.setToolTip(
            "How many genes the design may weaken to half their wild-type activity — a "
            "promoter or RBS change rather than a deletion."
        )
        targets_row.addWidget(QLabel("gene knock-outs"))
        targets_row.addWidget(self.jev_knockouts_spin)
        targets_row.addWidget(QLabel("gene knock-downs (50%)"))
        targets_row.addWidget(self.jev_knockdowns_spin)
        targets_row.addStretch(1)
        form.addRow("Regulation targets:", targets_row)

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
        self.jev_distinct_check = QCheckBox("each round must find a different design")
        self.jev_distinct_check.setChecked(True)
        self.jev_distinct_check.setToolTip(
            "Withholds one reaction from each design already found, so a later round has to "
            "reach somewhere else — the same integer cut OptKnock uses to enumerate "
            "alternatives. Without it, six rounds produced two distinct designs and four "
            "exact repeats: every round starts from the same wild type and plays the same "
            "game. The best design is still whichever round found it."
        )
        limits_row.addWidget(self.jev_distinct_check)
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
            "- NADPH supply is the cofactor I expect to be limiting."
        )
        self.jev_brief.setToolTip(
            "Guidance, not a rule. The agent is shown this on every step and weighs it, and "
            "it may still disagree with you — a brief cannot widen or narrow what CMM allows. "
            "For something you will not build, use the box below, which CMM enforces."
        )
        self.jev_brief.setMaximumHeight(78)
        form.addRow("Brief for the agent (guidance):", self.jev_brief)

        # The other half of what people write into a brief: not "here is what I know" but
        # "here is what I will not build". That is not a matter of opinion and should not be
        # left to the agent's judgement, so it is a separate box and CMM enforces it.
        self.jev_off_limits = QLineEdit()
        self.jev_off_limits.setPlaceholderText(
            "e.g.  ldhA, PFL, Pentose Phosphate Pathway   (comma separated)"
        )
        self.jev_off_limits.setToolTip(
            "Genes, reactions or subsystems this run may not touch. They are taken off the "
            "board before the agent sees them, and no proven design containing one is "
            "offered, so the constraint holds whatever the agent would have preferred. "
            "Matches a reaction id, a gene id, a gene name or a subsystem name. A name that "
            "matches nothing in the model stops the run rather than being ignored."
        )
        form.addRow("Off limits (enforced):", self.jev_off_limits)

        # A real web search, through OpenRouter, run once per reaction the agent is about to
        # act on: what has been published about editing it for this product, what it cost,
        # and what side effects a flux model cannot predict. The answer is pasted into that
        # reaction's record as data; it cannot widen the move vocabulary.
        web_row = QHBoxLayout()
        self.jev_web_check = QCheckBox(
            "Look up published evidence on the web (once, ~30 s)"
        )
        self.jev_web_check.setToolTip(
            "One web search through OpenRouter, before the first move: which genes have been "
            "deleted or down-regulated to raise this product in this organism, what it cost "
            "in growth, what side effects were reported, and which popular targets did not "
            "work. The answer is shown to the agent on every step.\n\n"
            "It used to be a lookup per candidate reaction, run inside a step — correct, and "
            "unusable: a search takes tens of seconds and eight of them turned a one-minute "
            "run into a ten-minute one."
        )
        self.jev_organism = QLineEdit()
        self.jev_organism.setPlaceholderText("Escherichia coli")
        self.jev_organism.setToolTip(
            "Required for the literature lookup. There is no default: asking the published "
            "record about the wrong species returns an answer that is confident and wrong."
        )
        self.jev_organism.setEnabled(False)
        self.jev_web_check.toggled.connect(self.jev_organism.setEnabled)
        self.jev_web_check.toggled.connect(self._on_jev_web_toggled)
        self.jev_organism.textChanged.connect(lambda _: self._update_jev_run_state())
        web_row.addWidget(self.jev_web_check)
        web_row.addWidget(QLabel("organism:"))
        web_row.addWidget(self.jev_organism, 1)
        form.addRow("Web research:", web_row)

        layout.addWidget(controls)

        run_row = QHBoxLayout()
        self.jev_run_btn = QPushButton("Let the agent play")
        self.jev_run_btn.clicked.connect(self.run_jev_agent)
        # These live OUTSIDE the controls group. The group is disabled wholesale while a run
        # is in flight, and Qt keeps a disabled widget's children disabled whatever you then
        # say about them — so a Stop button parented to it was grey and unclickable for
        # exactly as long as there was something to stop.
        self.jev_stop_btn = QPushButton("Stop")
        self.jev_stop_btn.setEnabled(False)
        self.jev_stop_btn.setToolTip(
            "Finish the current step and stop. Everything played so far is kept and "
            "summarised; only the baseline comparison is skipped."
        )
        self.jev_stop_btn.clicked.connect(self.stop_jev_agent)
        # The credential lives in a menu, which is where nobody looked for it. It is also a
        # prerequisite for this button working at all, so it belongs beside the button.
        self.jev_key_btn = QPushButton("Set key…")
        self.jev_key_btn.clicked.connect(self.set_jev_api_key)
        self.jev_key_label = QLabel("")
        self.jev_key_label.setStyleSheet("color: #5a6b7c; font-size: 11px;")
        # The run already writes a full bundle to a scratch directory; this copies it
        # somewhere the user chose and opens the reading copy. Before this the desktop run
        # wrote nothing at all, so a run watched on screen left no record of itself.
        self.jev_save_btn = QPushButton("Save full report…")
        self.jev_save_btn.setEnabled(False)
        self.jev_save_btn.setToolTip(
            "Write the whole run — report.html plus every table, the agent transcript and "
            "the provenance — to a folder you choose."
        )
        self.jev_save_btn.clicked.connect(self.save_jev_report)
        run_row.addWidget(self.jev_run_btn)
        run_row.addWidget(self.jev_stop_btn)
        run_row.addWidget(self.jev_save_btn)
        run_row.addSpacing(16)
        run_row.addWidget(self.jev_key_btn)
        run_row.addWidget(self.jev_key_label)
        run_row.addStretch(1)
        run_bar = QWidget()
        run_bar.setLayout(run_row)
        layout.addWidget(run_bar)

        # The reason the run button is refusing, beside the button, rather than in a tooltip
        # nobody hovers over. A disabled button that does not say why reads as a bug.
        self.jev_hint = QLabel("")
        self.jev_hint.setWordWrap(True)
        self.jev_hint.setStyleSheet("color: #a0342a; font-size: 11px;")
        self.jev_hint.setVisible(False)
        layout.addWidget(self.jev_hint)

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
        # What the run has spent, live. A decision is about $0.00016, so the number is small
        # — which is exactly why it is worth showing rather than leaving people to guess at
        # an order of magnitude they would reasonably assume was larger.
        self.jev_cost = QLabel("$0.0000")
        self.jev_cost.setMinimumWidth(150)
        self.jev_cost.setAlignment(Qt.AlignCenter)
        self.jev_cost.setToolTip(
            "OpenRouter spend on this run, summed from what each decision actually cost."
        )
        self.jev_cost.setStyleSheet(
            "QLabel { border: 1px solid #c2ccd6; border-radius: 4px; background: #f2f5f8; "
            "color: #23313f; font-weight: bold; padding: 4px; }"
        )
        self.jev_cost.setMinimumHeight(26)
        bars.addWidget(self.jev_cost, 1)
        layout.addLayout(bars)

        self.jev_summary = QLabel(
            "The agent is shown the metabolic state and picks one move at a time from a fixed "
            "set \u2014 delete a gene, weaken it to half, or run an analysis first. CMM "
            "executes each move, re-solves, and redraws the flux map, so you watch the "
            "design happen."
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
        summary_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        # Dragged, not capped. A fixed 120 pixels kept the summary from eating the flux map,
        # and cut the design list off mid-line — the design being the one thing a reader is
        # here for. A splitter gives the map its floor and still lets the summary be opened
        # up to read.
        summary_area.setMinimumHeight(96)

        results = QSplitter(Qt.Horizontal)
        map_box = QVBoxLayout()
        # The move's name as a widget, not as text drawn into the figure. A suptitle is laid
        # out in points against a drawing measured in inches, so it grew relative to the map
        # whenever the canvas rescaled the figure — and it cost a band of blank figure across
        # the top that the map could have used. A label is crisp at any panel size and free.
        self.jev_map_title = QLabel("")
        self.jev_map_title.setStyleSheet(
            "font-weight: bold; font-size: 13px; color: #23313f; padding: 2px 4px;"
        )
        map_box.addWidget(self.jev_map_title)
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
        self.jev_table.setSelectionBehavior(QTableWidget.SelectRows)
        # Clicking a move puts that move's flux distribution back on the map. This is also
        # how you see that a move changed nothing: step through the rows and the picture
        # holds still, which is the honest answer and not a broken redraw.
        self.jev_table.itemSelectionChanged.connect(self._show_selected_jev_move)

        # What the agent weighed, above what it did. A ``choice`` answer carries a
        # probability for every option it was offered, so each move records a complete
        # ranking of the board — showing it turns "the agent picked ICL" into "it picked ICL
        # over FUM by 0.40 to 0.16", which is the difference between an assertion and
        # evidence, and is often where the interesting runner-up is.
        self.jev_thinking = QLabel("The agent has not been asked anything yet.")
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

        # The move log and the target report answer different questions — "what did it do?"
        # and "what did it learn?" — and a reader wants the second one after the run, not
        # while it is playing. Tabs rather than a third pane: the map needs the width.
        self.jev_targets = QTableWidget(0, 5)
        self.jev_targets.setHorizontalHeaderLabels(
            ["Target", "Edit", "Best gain", "For", "Against"]
        )
        target_header = self.jev_targets.horizontalHeader()
        for column, mode in enumerate(
            (
                QHeaderView.ResizeToContents,
                QHeaderView.ResizeToContents,
                QHeaderView.ResizeToContents,
                QHeaderView.Stretch,
                QHeaderView.Stretch,
            )
        ):
            target_header.setSectionResizeMode(column, mode)
        self.jev_targets.verticalHeader().setVisible(False)
        self.jev_targets.setAlternatingRowColors(True)
        self.jev_targets.setWordWrap(True)

        # One row per round: what it was allowed to use, what it built, and what it left.
        # A multi-round run is a portfolio, and a portfolio nobody can see is one design.
        self.jev_rounds = QTableWidget(0, 6)
        self.jev_rounds.setHorizontalHeaderLabels(
            [
                "Round",
                "Question it answered",
                "Engineering",
                "Product",
                "Growth",
                "Left undone",
            ]
        )
        rounds_header = self.jev_rounds.horizontalHeader()
        for column, mode in enumerate(
            (
                QHeaderView.ResizeToContents,
                QHeaderView.ResizeToContents,
                QHeaderView.Stretch,
                QHeaderView.ResizeToContents,
                QHeaderView.ResizeToContents,
                QHeaderView.Stretch,
            )
        ):
            rounds_header.setSectionResizeMode(column, mode)
        self.jev_rounds.verticalHeader().setVisible(False)
        self.jev_rounds.setAlternatingRowColors(True)
        self.jev_rounds.setWordWrap(True)

        self.jev_lower_tabs = QTabWidget()
        self.jev_lower_tabs.addTab(self.jev_table, "Moves")
        self.jev_lower_tabs.addTab(self.jev_rounds, "Rounds")
        self.jev_lower_tabs.addTab(self.jev_targets, "Targets: for and against")

        right = QSplitter(Qt.Vertical)
        right.addWidget(panel)
        right.addWidget(self.jev_lower_tabs)
        right.setSizes([260, 240])
        right.setMinimumWidth(500)
        results.addWidget(right)
        results.setStretchFactor(0, 3)
        results.setStretchFactor(1, 2)
        results.setSizes([780, 520])

        body = QSplitter(Qt.Vertical)
        body.addWidget(summary_area)
        body.addWidget(results)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([150, 620])
        layout.addWidget(body, 1)

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
        self._refresh_jev_key_label()
        self._update_jev_run_state()
        if not exchanges:
            self.jev_summary.setText(
                "This model has no exchange reactions, so there is no product flux to raise."
            )
        elif not has_key:
            self.jev_summary.setText(_NO_KEY_MESSAGE)

        self.jev_table.setRowCount(0)
        self.jev_rounds.setRowCount(0)
        self.jev_targets.setRowCount(0)
        self._jev_frames = []
        self._jev_result = None
        self._jev_cost = 0.0
        self._jev_run_dir = None
        self.jev_save_btn.setEnabled(False)
        self.jev_cost.setText("$0.0000")

    def _on_jev_web_toggled(self, checked: bool) -> None:
        """Turning the lookup on puts the cursor where the one missing thing goes."""

        if checked and not self.jev_organism.text().strip():
            self.jev_organism.setFocus()
        self._update_jev_run_state()

    def _update_jev_run_state(self) -> None:
        """Enable the run only when everything it needs is present, and say what is missing.

        The saying is the part that was wrong. The requirement was real — a literature lookup
        aimed at the wrong species returns an answer that is confident and wrong — but it
        lived in a tooltip, so ticking the web-research box greyed the run button out with no
        visible reason, which reads as a broken button rather than as a question.
        """

        if not hasattr(self, "jev_run_btn"):
            return
        needs_organism = (
            self.jev_web_check.isChecked() and not self.jev_organism.text().strip()
        )
        has_key = getattr(self, "_jev_has_key", False)
        has_exchanges = getattr(self, "_jev_has_exchanges", False)
        ready = has_exchanges and has_key and not needs_organism
        self.jev_run_btn.setEnabled(ready)

        if not has_exchanges:
            hint = (
                "This model has no exchange reactions, so there is no product to raise."
            )
        elif not has_key:
            hint = (
                "No OpenRouter API key. Press “Set key…” \u2014 it is the only credential "
                "CMM uses, and only this tab needs it."
            )
        elif needs_organism:
            hint = (
                "Name the organism for the web lookup, for example “Escherichia coli”. "
                "There is no default on purpose: the published record answered about the "
                "wrong species is confident and wrong. Untick the box to run without it."
            )
        else:
            hint = ""
        if hasattr(self, "jev_hint"):
            self.jev_hint.setText(hint)
            self.jev_hint.setVisible(bool(hint))
        self.jev_run_btn.setToolTip(hint)
        self.jev_organism.setStyleSheet(
            "border: 1px solid #a0342a;" if needs_organism else ""
        )

    def _refresh_jev_key_label(self) -> None:
        """Say where the key in force came from, beside the button that sets it."""

        from cmm.jev import credentials

        if not hasattr(self, "jev_key_label"):
            return
        source = credentials.key_source()
        if source == "environment":
            import os

            text = f"key: {credentials.ENV_VAR} ({credentials.masked(os.environ[credentials.ENV_VAR])})"
        elif source == "saved":
            text = f"key: saved ({credentials.masked(credentials.stored_key())})"
        else:
            text = "no key set"
        self.jev_key_label.setText(text)
        self.jev_key_btn.setText("Change key…" if source != "none" else "Set key…")

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
            "The agent is the only part of CMM that needs a key.\n\n"
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
            f"Saved to {saved}.\nThe Agent tab is ready.",
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
                "No key is set. The Agent tab is disabled until one is.\n\n"
                f"Either set {credentials.ENV_VAR} in the environment, or use "
                "Agent \u25b8 Set OpenRouter API Key."
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
        self._jev_cost = 0.0
        self._jev_stop_requested = False
        self.jev_cost.setText("$0.0000")
        self.jev_round_progress.setRange(0, self.jev_rounds_spin.value())
        self.jev_round_progress.setValue(0)
        self.jev_round_progress.setFormat(f"round 0 of {self.jev_rounds_spin.value()}")
        self.jev_progress.setRange(0, self.jev_ticks_spin.value())
        self.jev_progress.setValue(0)
        self.jev_progress.setFormat(
            f"step 0 of up to {self.jev_ticks_spin.value()} this round"
        )
        self.jev_weights.setRowCount(0)
        self.jev_map_title.setText("")
        self.jev_thinking.setText("Asking the agent for its first move\u2026")
        self.jev_summary.setText(f"The agent is playing for {html.escape(product)}…")

        # Serialize on the UI thread: the worker must not touch this model's solver object.
        from cobra.io import write_sbml_model

        scratch = Path(tempfile.mkdtemp(prefix="cmm-agent-"))
        model_path = scratch / "model.xml"
        write_sbml_model(self.model, str(model_path))
        self._jev_run_dir = scratch / "run"

        config = JevConfig(
            model_path=model_path,
            # Always write the bundle. It costs a directory of CSV and it is the difference
            # between a run someone watched and a run someone can check afterwards.
            output_dir=self._jev_run_dir,
            overwrite=True,
            product=product,
            rounds=self.jev_rounds_spin.value(),
            steps_per_round=self.jev_ticks_spin.value(),
            brief=self.jev_brief.toPlainText(),
            max_knockouts=self.jev_knockouts_spin.value(),
            max_knockdowns=self.jev_knockdowns_spin.value(),
            growth_floor=self.jev_growth_spin.value(),
            candidate_limit=self.jev_board_spin.value(),
            require_distinct_rounds=self.jev_distinct_check.isChecked(),
            enable_web_research=self.jev_web_check.isChecked(),
            organism=self.jev_organism.text().strip(),
            off_limits=tuple(
                name.strip()
                for name in self.jev_off_limits.text().split(",")
                if name.strip()
            ),
        )
        bridge = self._jev_bridge

        def _compute():
            return run_jev_design(
                config,
                on_tick=lambda tick, fluxes: bridge.played.emit(tick, dict(fluxes)),
                # Read from the worker thread, written from the UI thread. A plain flag is
                # enough: it is set once, never cleared mid-run, and a step's delay in
                # observing it costs at most one more move.
                should_stop=lambda: self._jev_stop_requested,
            )

        try:
            result = self._run_jev_visibly(_compute)
        except Exception as exc:
            self.jev_summary.setText(f"The agent run failed: {html.escape(str(exc))}")
            self.status_label.setText("Agent run failed.")
            return

        self._jev_result = result
        self._show_jev_summary(result)
        self._fill_jev_rounds(result)
        self._fill_jev_targets(result)
        self.jev_save_btn.setEnabled(result.run_directory is not None)

    def save_jev_report(self) -> None:
        """Copy the whole run somewhere the user chose, and open the reading copy.

        The bundle is written during the run to a scratch directory, so this is a copy rather
        than a re-render: what gets saved is exactly what was scored, not a second pass over
        the result that could disagree with it.
        """

        from qtpy.QtWidgets import QFileDialog

        source = getattr(self, "_jev_run_dir", None)
        if self._jev_result is None or source is None or not Path(source).is_dir():
            QMessageBox.information(
                self,
                "Save report",
                "Run the agent first; there is no run directory to save.",
            )
            return

        chosen = QFileDialog.getExistingDirectory(
            self, "Where should the run be saved?", str(Path.home())
        )
        if not chosen:
            return

        import shutil

        product = str(self._jev_result.summary()["product"])
        stamp = datetime.now().strftime("%Y%m%dT%H%M%S")
        target = Path(chosen) / f"agent-{product}-{stamp}"
        try:
            shutil.copytree(source, target)
        except OSError as error:
            QMessageBox.warning(
                self, "Save report", f"The run could not be saved: {error}"
            )
            return

        report = target / "report.html"
        QMessageBox.information(
            self,
            "Save report",
            f"Saved to:\n{target}\n\nOpen report.html for the whole run in one page; the "
            "CSV tables, the agent transcript and the provenance are beside it.",
        )
        if report.exists():
            from qtpy.QtCore import QUrl
            from qtpy.QtGui import QDesktopServices

            QDesktopServices.openUrl(QUrl.fromLocalFile(str(report)))
        self.status_label.setText(f"Agent run saved to {target}.")

    def stop_jev_agent(self) -> None:
        """Ask the running game to stop after the step it is on.

        Not a kill: the worker is inside a solver call for most of its life, and interrupting
        that would lose the run. The engine checks this flag once per step, so the wait is one
        move — about a second — and everything played is kept and summarised.
        """

        if not getattr(self, "_jev_stop_requested", False):
            self._jev_stop_requested = True
            self.jev_stop_btn.setEnabled(False)
            self.jev_stop_btn.setText("Stopping…")
            self.jev_progress.setFormat("stopping after this step…")
            self.status_label.setText("Stopping the agent run after the current step.")

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

        self._jev_running = True
        self.jev_run_btn.setEnabled(False)
        self.jev_controls.setEnabled(False)
        self.jev_stop_btn.setEnabled(True)
        self.jev_stop_btn.setText("Stop")
        QTimer.singleShot(0, thread.start)
        try:
            loop.exec_()
        finally:
            thread.quit()
            thread.wait()
            worker.deleteLater()
            self._jev_running = False
            self.jev_controls.setEnabled(True)
            self.jev_run_btn.setEnabled(True)
            self.jev_stop_btn.setEnabled(False)
            self.jev_stop_btn.setText("Stop")
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

        previous = self._jev_frames[-1][1] if self._jev_frames else None
        self._jev_frames.append((tick, fluxes))
        self._jev_cost = getattr(self, "_jev_cost", 0.0) + float(
            tick.decision_cost_usd or 0.0
        )
        self.jev_cost.setText(
            f"${self._jev_cost:.4f}  ·  {len(self._jev_frames)} steps"
        )
        self._append_jev_row(tick)
        self._draw_jev_map(tick, fluxes, previous)
        self._show_jev_thinking(tick)
        self._advance_jev_progress(tick)
        self.status_label.setText(tick.headline())

    def _show_selected_jev_move(self) -> None:
        """Put the selected move's flux distribution back on the map.

        Only between runs. While one is playing, the map already follows the newest move, and
        redrawing from a selection signal meant tearing down and rebuilding the canvas from
        inside the row insert that raised the signal — which aborted the process, worker
        thread and all, rather than failing in any way a reader could connect to its cause.
        """

        if getattr(self, "_jev_running", False):
            return
        rows = {index.row() for index in self.jev_table.selectedIndexes()}
        if len(rows) != 1:
            return
        row = rows.pop()
        if not 0 <= row < len(self._jev_frames):
            return
        tick, fluxes = self._jev_frames[row]
        previous = self._jev_frames[row - 1][1] if row else None
        self._draw_jev_map(tick, fluxes, previous)
        self._show_jev_thinking(tick)
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
        # Silenced while the row goes in: inserting and scrolling both move the selection,
        # and the selection handler redraws the map.
        self.jev_table.blockSignals(True)
        try:
            self._append_jev_row_unguarded(tick)
        finally:
            self.jev_table.blockSignals(False)

    def _append_jev_row_unguarded(self, tick) -> None:
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

    @staticmethod
    def _flux_change(fluxes, previous) -> str:
        """One sentence on what moved since the previous step, or that nothing did.

        This exists because of a fair question: the map is redrawn after every move, so why
        does it so often look identical? Because it *is* identical. A scan changes nothing by
        definition, a move that breached the growth floor was reverted before this frame was
        taken, and a deletion of a reaction already carrying no flux is a real move with no
        immediate consequence. Only a move that stuck and mattered redraws differently, and
        on a typical round that is two or three steps out of twenty.

        Saying so is the fix. A picture that has not changed and does not admit it reads as a
        broken redraw; the same picture with "no flux changed on this step" under it reads as
        the result it is.
        """

        if previous is None:
            return "Wild-type flux distribution, before any move."
        moved = sorted(
            (
                (abs(value - float(previous.get(rid, 0.0))), rid, value)
                for rid, value in fluxes.items()
                if abs(value - float(previous.get(rid, 0.0))) > 1e-6
            ),
            reverse=True,
        )
        if not moved:
            return (
                "No flux changed on this step — the move was a scan, was reverted by the "
                "rules, or touched a reaction that was already carrying nothing. The map is "
                "identical to the previous step on purpose."
            )
        biggest = ", ".join(
            f"{rid} {float(previous.get(rid, 0.0)):.3g}\u2192{value:.3g}"
            for _, rid, value in moved[:3]
        )
        return (
            f"{len(moved)} reaction{'' if len(moved) == 1 else 's'} changed flux on this "
            f"step. Largest: {biggest}."
        )

    def _draw_jev_map(self, tick, fluxes, previous=None) -> None:
        """Redraw the flux map for this move, on the curated map when the model has one."""

        from cmm.visualization import escher_flux_map, network_flux_map

        # "end_round on end_round" is not a move description. Name the reaction only when
        # the move is about one, which also keeps the title short enough not to be clipped.
        move = tick.action or "no move"
        if tick.target != tick.action:
            move = f"{move} on {tick.target}"
        title = (
            f"Round {tick.round_index}, step {tick.tick_index}  ·  {move}  ·  "
            f"product {tick.product_flux:.4g}, growth {tick.growth:.4g}"
        )
        # Author the figure at the width it will be shown at. Type is measured in points and
        # the drawing in inches, so a figure drawn at 9 inches and stretched into a 6.5-inch
        # panel comes out with its text 1.4x too large against the network — which is exactly
        # what it looked like. The canvas is the authority on that width.
        canvas = getattr(self, "_jev_canvas", None)
        panel_px = canvas.width() if canvas is not None else 0
        if panel_px < 200:
            panel_px = max(self.jev_canvas_holder.geometry().width(), 620)
        width = max(5.0, min(panel_px / 100.0, 11.0))
        try:
            if self._map_path:
                # Metabolite labels off, and reactions at rest left unnamed. At full size the
                # curated map reads well; in a panel a third of that, naming all ninety-five
                # reactions spends the space on the half that are doing nothing and the names
                # collide into smudges. Zoom or the Flux Map tab gives the full figure.
                figure = escher_flux_map(
                    self._map_path,
                    dict(fluxes),
                    title="",
                    width=width,
                    label_metabolites=False,
                    label_min_fraction=0.01,
                    font_scale=0.85,
                )
            else:
                figure = network_flux_map(self.model, dict(fluxes), title="")
            self.jev_map_title.setText(title)
        except Exception as exc:  # a map that cannot be drawn must not stop the game
            self.status_label.setText(f"{tick.headline()} (map not drawn: {exc})")
            return
        self._set_figure(self.jev_canvas_holder, "jev", figure)
        what = (
            "Curated Escher map, coloured and widened by flux. Only reactions carrying flux "
            "are named at this size, and metabolites are not \u2014 zoom with the magnifier, "
            "or use the Flux Map tab for the full figure."
            if self._map_path
            else "Schematic of the highest-flux reactions; no curated Escher map fits this "
            "model. Arrow colour and width both scale with the flux magnitude."
        )
        self.jev_map_caption.setText(f"{self._flux_change(fluxes, previous)}  {what}")

    def _fill_jev_rounds(self, result) -> None:
        """What each round was allowed to do, what it engineered, and what it left behind.

        The engineering column is the point of the tab: a round's design, named by its genes,
        is the thing a reader takes away. The question column says why the rounds differ —
        each one is barred from a reaction an earlier design used, so a later row reads "the
        best design without X", which is the design a laboratory that cannot edit X needs.
        """

        records = result.rounds
        self.jev_rounds.setRowCount(len(records))
        best = max((r.product_flux for r in records), default=0.0)
        for row, record in enumerate(records):
            question = (
                "best available"
                if not record.withheld
                else "best without " + ", ".join(record.withheld)
            )
            design = "\n".join(record.interventions) or "nothing was applied"
            if record.repeated is not None:
                design += f"\n(the same design round {record.repeated} found)"
            values = [
                f"R{record.round_index}",
                question,
                design,
                f"{record.product_flux:.4g}",
                f"{record.growth:.4g}",
                "\n".join(f"\u2022 {line}" for line in record.shortfall)
                or record.stopped_because,
            ]
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                item.setTextAlignment(Qt.AlignLeft | Qt.AlignTop)
                if record.product_flux >= best - 1e-9 and best > 0:
                    item.setBackground(QColor("#dce9d6"))
                self.jev_rounds.setItem(row, column, item)
        self.jev_rounds.resizeRowsToContents()
        distinct = len({record.signature for record in records})
        self.jev_lower_tabs.setTabText(
            1, f"Rounds ({distinct} distinct design{'' if distinct == 1 else 's'})"
        )

    def _fill_jev_targets(self, result) -> None:
        """Every target the run weighed, with the case for and against editing it.

        Not scored against each other on purpose. "Raises the product by 0.035" and "needs
        three isozymes deleted" are not the same kind of quantity, and collapsing them into a
        rank would hide the trade from the only person who can make it.
        """

        reports = result.targets()
        self.jev_targets.setRowCount(len(reports))
        for row, report in enumerate(reports):
            gains = [
                g
                for g in (report.deletion_gain, report.knockdown_gain)
                if g is not None
            ]
            best = f"{max(gains):+.4g}" if gains else "not measured"
            values = [
                report.reaction_id,
                ", ".join(report.genes) or "—",
                best,
                "\n".join(f"+ {line}" for line in report.pros) or "—",
                "\n".join(f"\u2212 {line}" for line in report.cons) or "—",
            ]
            for column, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                item.setTextAlignment(Qt.AlignLeft | Qt.AlignTop)
                if report.in_best_design:
                    item.setBackground(QColor("#dce9d6"))
                self.jev_targets.setItem(row, column, item)
        self.jev_targets.resizeRowsToContents()
        if reports:
            self.jev_lower_tabs.setTabText(
                2, f"Targets: for and against ({len(reports)})"
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
        # What each round left undone, which is the half of a round's result that tells the
        # next one what to do. A run that only reports its best score hides its own agenda.
        rounds = "".join(
            (
                "<li><b>round {index}</b> — {question}: "
                "<b>{product:.4g}</b> at growth {growth:.4g}. {stopped}{gap}</li>"
            ).format(
                index=record.round_index,
                question=html.escape(
                    "best design available"
                    if not record.withheld
                    else "best design without " + ", ".join(record.withheld)
                ),
                product=record.product_flux,
                growth=record.growth,
                stopped=html.escape(record.stopped_because or "ended"),
                gap=(
                    "; left undone: " + html.escape("; ".join(record.shortfall))
                    if record.shortfall
                    else ""
                ),
            )
            for record in result.rounds
        )
        distinct = len({record.signature for record in result.rounds})
        round_block = (
            f"<br>Rounds — {distinct} distinct design"
            f"{'' if distinct == 1 else 's'} over {len(result.rounds)}:<ul>{rounds}</ul>"
            if rounds
            else ""
        )

        notes = "".join(
            f"<li>{html.escape(str(note))}</li>" for note in summary["notes"]
        )
        note_block = f"<br>Notes:<ul>{notes}</ul>" if notes else ""

        self.jev_summary.setText(
            f"{headline}<br>Design:<ul>{design}</ul>"
            f"{summary['n_ticks']} moves over {summary['n_rounds']} rounds; "
            f"{usage['calls']} agent decisions costing ${usage['cost_usd']:.4f}."
            f"{round_block}{note_block}"
            "<br><i>This is a computational hypothesis. The agent's choices are not "
            "guaranteed to repeat on a re-run; the CMM solves behind them are.</i>"
        )
        self.jev_cost.setText(
            f"${usage['cost_usd']:.4f}  \u00b7  {usage['calls']} decisions"
        )
        self.status_label.setText(
            f"Agent run complete: {summary['n_ticks']} moves, best {product} "
            f"{best:.4g} mmol gDW-1 h-1, ${usage['cost_usd']:.4f} spent."
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
