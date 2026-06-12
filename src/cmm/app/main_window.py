"""CMM desktop main window.

A thin Qt view over the solver-neutral services. All numerics live in cmm.core / cmm.features
/ cmm.omics; this module only renders results and wires controls, so the GUI inherits none of
the analysis logic.
"""

from __future__ import annotations

import html
from collections.abc import Mapping

from cobra import Model
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor, QFont
from qtpy.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSlider,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from cmm.core import (
    PRESET_MEDIA,
    active_solver,
    apply_medium,
    fba,
    fva,
    pfba,
    reference_state_pfba,
    solver_status,
)
from cmm.features.comparison import moma, reference_flux, room
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.revert import revert_targets
from cmm.omics.differential import differential_expression
from cmm.omics.expression import gene_to_reaction_weights, integrate_expression
from cmm.visualization import (
    escher_flux_map,
    fseof_figure,
    fvseof_figure,
    production_envelope_figure,
)

_STYLE = """
QMainWindow, QWidget { background: #f4f6f8; color: #1c2733; font-size: 12px; }
QLabel { background: transparent; }
QLabel#title { color: white; font-size: 19px; font-weight: 600; }
QLabel#subtitle { color: #d7e3ef; font-size: 12px; }
QWidget#header { background: #1f6feb; }
QGroupBox { font-weight: 600; border: 1px solid #d0d7de; border-radius: 6px;
            margin-top: 10px; padding: 8px; background: white; }
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 4px; }
QPushButton { background: #1f6feb; color: white; border: none; border-radius: 5px;
              padding: 6px 14px; font-weight: 600; }
QPushButton:hover { background: #1a5fd0; }
QTableWidget { background: white; gridline-color: #e4e8ec;
               selection-background-color: #cfe3ff; selection-color: #11243b; }
QHeaderView::section { background: #eaeef2; padding: 5px; border: none;
                       border-right: 1px solid #d8dde2; font-weight: 600; }
QTabBar::tab { padding: 7px 16px; background: #e7ebef; }
QTabBar::tab:selected { background: white; border-top: 2px solid #1f6feb; }
"""


def _read_expression_vector(path: str) -> dict[str, float]:
    """Read a two-column gene-expression CSV/TSV file into ``gene -> value``."""

    import pandas as pd

    frame = pd.read_csv(path, sep=None, engine="python")
    if frame.shape[1] < 2:
        raise ValueError("expression file must have at least two columns: gene, expression")
    gene_col, value_col = frame.columns[:2]
    expression: dict[str, float] = {}
    for gene, value in zip(frame[gene_col], frame[value_col], strict=False):
        if pd.isna(gene) or pd.isna(value):
            continue
        expression[str(gene)] = float(value)
    if not expression:
        raise ValueError("expression file contains no numeric expression values")
    return expression


class CmmMainWindow(QMainWindow):
    """Main workbench window over a single cobra model."""

    def __init__(self, model: Model, parent=None, *, default_product: str | None = None,
                 map_path: str | None = None):
        super().__init__(parent)
        self.model = model
        self._fluxes: dict[str, float] = {}
        self._loading = False
        self._fluxes_stale = False
        self._default_product = default_product
        self._map_path = map_path
        self._revert_source_expression: dict[str, float] | None = None
        self._revert_target_expression: dict[str, float] | None = None
        self.setWindowTitle("CMM — Cellular Metabolic Modeling Workbench")
        self.resize(1180, 740)
        self.setStyleSheet(_STYLE)
        self._build_ui()
        self.load_model(model)

    # -- construction -------------------------------------------------------

    def _build_ui(self) -> None:
        self._build_menus()
        root = QWidget()
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._build_header())

        body = QWidget()
        body_layout = QHBoxLayout(body)
        body_layout.setContentsMargins(12, 12, 12, 12)
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._build_model_panel())
        splitter.addWidget(self._build_tabs())
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 5)
        body_layout.addWidget(splitter)
        outer.addWidget(body, 1)

        self.status_label = QLabel("Ready.")
        self.status_label.setContentsMargins(12, 4, 12, 4)
        outer.addWidget(self.status_label)

    def _build_menus(self) -> None:
        """A desktop-style menu bar: Analysis / Model / Config drive the workbench."""

        bar = self.menuBar()

        analysis = bar.addMenu("&Analysis")
        analysis.addAction("Run FBA", self.run_fba)
        analysis.addAction("Run FVA", self.run_fva)
        analysis.addSeparator()
        analysis.addAction("Theoretical Yield",
                           lambda: self._in_tab("Production", self.run_theoretical_yield))
        analysis.addAction("Production Envelope",
                           lambda: self._in_tab("Production", self.run_production_envelope_plot))
        analysis.addAction("FSEOF Targets",
                           lambda: self._in_tab("Production", self.run_fseof_plot))
        analysis.addSeparator()
        analysis.addAction("Omics Integration (E-Flux2 / LAD)…",
                           lambda: self._goto_tab("Omics"))
        analysis.addAction("Revert Metabolism…", lambda: self._goto_tab("Revert Metabolism"))
        if self._map_path:
            analysis.addSeparator()
            analysis.addAction("Render Flux Map",
                               lambda: self._in_tab("Flux Map", self.render_flux_map))

        model_menu = bar.addMenu("&Model")
        model_menu.addAction("Model Info…", self._show_model_info)

        config = bar.addMenu("&Config")
        status = solver_status(self.model)
        solver_action = config.addAction(f"Solver: {status.summary()}")
        solver_action.setEnabled(False)
        config.addAction("Solver status…", self._show_solver_status)

    def _tab_index(self, name: str) -> int | None:
        return next(
            (i for i in range(self.tabs.count()) if self.tabs.tabText(i) == name), None
        )

    def _goto_tab(self, name: str) -> None:
        index = self._tab_index(name)
        if index is not None:
            self.tabs.setCurrentIndex(index)

    def _in_tab(self, name: str, action) -> None:
        self._goto_tab(name)
        action()

    def _show_solver_status(self) -> None:
        status = solver_status(self.model)
        lines = [
            f"Active solver: {status.name}",
            f"Capabilities: {', '.join(status.capabilities)}",
            f"Recommended (full QP+MILP+MIQP): {'yes' if status.recommended else 'no'}",
            f"Available solvers: {', '.join(status.available)}",
        ]
        if status.warning:
            lines.append("")
            lines.append("⚠ " + status.warning)
            QMessageBox.warning(self, "Solver status", "\n".join(lines))
        else:
            QMessageBox.information(self, "Solver status", "\n".join(lines))

    def _show_model_info(self) -> None:
        objective = ", ".join(
            r.id for r in self.model.reactions if r.objective_coefficient != 0
        )
        QMessageBox.information(
            self, "Model Info",
            f"{self.model.id}\n\n"
            f"Reactions: {len(self.model.reactions)}\n"
            f"Metabolites: {len(self.model.metabolites)}\n"
            f"Genes: {len(self.model.genes)}\n"
            f"Objective: {objective or '—'} ({self.model.objective_direction})\n"
            f"Solver: {active_solver(self.model)}",
        )

    def _build_header(self) -> QWidget:
        header = QWidget()
        header.setObjectName("header")
        layout = QVBoxLayout(header)
        layout.setContentsMargins(16, 12, 16, 12)
        title = QLabel("CMM — Cellular Metabolic Modeling Workbench")
        title.setObjectName("title")
        self.subtitle = QLabel("")
        self.subtitle.setObjectName("subtitle")
        layout.addWidget(title)
        layout.addWidget(self.subtitle)
        return header

    def _build_model_panel(self) -> QWidget:
        box = QGroupBox("Model")
        layout = QVBoxLayout(box)
        self.model_info = QLabel("")
        self.model_info.setWordWrap(True)
        layout.addWidget(self.model_info)

        slider_row = QHBoxLayout()
        slider_row.addWidget(QLabel("Flux range filter:"))
        self.flux_slider = QSlider(Qt.Horizontal)
        self.flux_slider.setRange(0, 100)
        self.flux_slider.setToolTip(
            "Highlight reactions whose |flux| is at least this fraction of the maximum "
            "flux (run FBA first)."
        )
        self.flux_slider.valueChanged.connect(self._apply_flux_threshold)
        self.flux_threshold_label = QLabel("0%")
        slider_row.addWidget(self.flux_slider, 1)
        slider_row.addWidget(self.flux_threshold_label)
        layout.addLayout(slider_row)

        self.reaction_table = QTableWidget(0, 4)
        self.reaction_table.setHorizontalHeaderLabels(["Reaction", "Lower", "Upper", "Flux"])
        # Reaction IDs vary in length (short on toy models, long like "Biomass_Ecoli_core" on
        # genome-scale ones); let the id column take the slack while numeric columns fit content.
        _rxn_header = self.reaction_table.horizontalHeader()
        _rxn_header.setSectionResizeMode(0, QHeaderView.Stretch)
        for _col in (1, 2, 3):
            _rxn_header.setSectionResizeMode(_col, QHeaderView.ResizeToContents)
        self.reaction_table.verticalHeader().setVisible(False)
        self.reaction_table.setAlternatingRowColors(True)
        self.reaction_table.itemChanged.connect(self._on_bound_edited)
        layout.addWidget(self.reaction_table, 1)
        hint = QLabel("Double-click a Lower/Upper cell to edit a bound, then re-run FBA.")
        hint.setStyleSheet("color: #5a6b7b; font-style: italic;")
        layout.addWidget(hint)
        return box

    def _build_tabs(self) -> QWidget:
        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_simulation_tab(), "Simulation")
        self.tabs.addTab(self._build_comparison_tab(), "Comparison")
        self.tabs.addTab(self._build_production_tab(), "Production")
        self.tabs.addTab(self._build_omics_tab(), "Omics")
        if self._map_path:
            self.tabs.addTab(self._build_fluxmap_tab(), "Flux Map")
        self.tabs.addTab(self._build_revert_tab(), "Revert Metabolism")
        return self.tabs

    def _build_omics_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Omics integration (expression → flux)")
        row = QHBoxLayout(controls)
        row.addWidget(QLabel("Method:"))
        self.omics_method_combo = QComboBox()
        self.omics_method_combo.addItems(["eflux2", "lad"])
        self.omics_method_combo.setToolTip(
            "E-Flux2: scale bounds by expression, minimize L2 flux. LAD: fit fluxes to expression."
        )
        row.addWidget(self.omics_method_combo)
        load_btn = QPushButton("Load expression CSV…")
        load_btn.clicked.connect(self.load_expression_csv)
        demo_btn = QPushButton("Run on demo expression")
        demo_btn.clicked.connect(self.run_omics_demo)
        row.addWidget(load_btn)
        row.addWidget(demo_btn)
        row.addStretch(1)
        layout.addWidget(controls)

        self.omics_summary = QLabel(
            "Load a gene-expression CSV (columns: gene, expression) or run on demo data."
        )
        self.omics_summary.setWordWrap(True)
        layout.addWidget(self.omics_summary)

        self.omics_table = QTableWidget(0, 2)
        self.omics_table.setHorizontalHeaderLabels(["Reaction", "Predicted flux"])
        self.omics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.omics_table.verticalHeader().setVisible(False)
        self.omics_table.setAlternatingRowColors(True)
        layout.addWidget(self.omics_table, 1)
        return tab

    def _build_fluxmap_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        controls = QHBoxLayout()
        render_btn = QPushButton("Render flux map")
        render_btn.clicked.connect(self.render_flux_map)
        controls.addWidget(render_btn)
        controls.addWidget(QLabel("Escher layout coloured by the current FBA flux."))
        controls.addStretch(1)
        layout.addLayout(controls)
        self.map_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.map_canvas_holder)
        layout.addWidget(holder, 1)
        self._map_canvas = None
        return tab

    def _build_production_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Production design")
        form = QHBoxLayout(controls)
        form.addWidget(QLabel("Target product:"))
        self.product_combo = QComboBox()
        self.product_combo.setMinimumWidth(160)
        self.product_combo.setToolTip("Exchange reaction of the target product")
        self.product_combo.currentTextChanged.connect(self._on_product_changed)
        form.addWidget(self.product_combo, 1)
        form.addWidget(QLabel("substrate:"))
        self.substrate_combo = QComboBox()
        self.substrate_combo.setMinimumWidth(140)
        self.substrate_combo.setToolTip("Carbon substrate (auto = detect from the medium)")
        form.addWidget(self.substrate_combo, 1)
        self.anaerobic_combo = QComboBox()
        self.anaerobic_combo.addItems(["aerobic", "anaerobic"])
        form.addWidget(self.anaerobic_combo)
        prod_fba_btn = QPushButton("Run FBA")
        prod_fba_btn.clicked.connect(self.run_fba)
        yield_btn = QPushButton("Theoretical yield")
        yield_btn.clicked.connect(self.run_theoretical_yield)
        envelope_btn = QPushButton("Production envelope")
        envelope_btn.clicked.connect(self.run_production_envelope_plot)
        fseof_btn = QPushButton("FSEOF targets")
        fseof_btn.clicked.connect(self.run_fseof_plot)
        fvseof_btn = QPushButton("FVSEOF (robust)")
        fvseof_btn.clicked.connect(self.run_fvseof_plot)
        form.addWidget(prod_fba_btn)
        form.addWidget(yield_btn)
        form.addWidget(envelope_btn)
        form.addWidget(fseof_btn)
        form.addWidget(fvseof_btn)
        self._production_buttons = [yield_btn, envelope_btn, fseof_btn, fvseof_btn]
        layout.addWidget(controls)

        self.yield_label = QLabel("Select a product and compute its theoretical yield.")
        self.yield_label.setFont(QFont("", 13, QFont.Bold))
        self.yield_label.setWordWrap(True)
        self.yield_label.setTextFormat(Qt.RichText)
        layout.addWidget(self.yield_label)

        self.production_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.production_canvas_holder)
        layout.addWidget(holder, 1)
        self._production_canvas = None
        return tab

    def _build_comparison_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        controls = QGroupBox("Perturbation response (MOMA / ROOM vs a reference template)")
        row = QHBoxLayout(controls)
        row.addWidget(QLabel("Method:"))
        self.comparison_method_combo = QComboBox()
        self.comparison_method_combo.addItems(["MOMA (L2)", "MOMA (L1)", "ROOM"])
        row.addWidget(self.comparison_method_combo)
        row.addWidget(QLabel("Reference template:"))
        self.template_combo = QComboBox()
        self.template_combo.addItems(["fba", "pfba", "lad", "eflux2"])
        self.template_combo.setToolTip("Wild-type flux template the perturbed state is compared to")
        row.addWidget(self.template_combo)
        row.addWidget(QLabel("Knock out:"))
        self.ko_combo = QComboBox()
        self.ko_combo.setMinimumWidth(140)
        row.addWidget(self.ko_combo)
        run_btn = QPushButton("Run")
        run_btn.clicked.connect(self.run_comparison)
        row.addWidget(run_btn)
        row.addStretch(1)
        layout.addWidget(controls)

        self.comparison_summary = QLabel(
            "Pick a reference template (FBA/pFBA/LAD/E-Flux2), a reaction to knock out, "
            "and a method to predict the perturbed flux state."
        )
        self.comparison_summary.setWordWrap(True)
        layout.addWidget(self.comparison_summary)

        self.comparison_table = QTableWidget(0, 3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )
        self.comparison_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.comparison_table.verticalHeader().setVisible(False)
        self.comparison_table.setAlternatingRowColors(True)
        layout.addWidget(self.comparison_table, 1)
        return tab

    def run_comparison(self) -> None:
        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        ko = self.ko_combo.currentText()
        try:
            expression = None
            if template in ("lad", "eflux2"):
                import numpy as np

                rng = np.random.default_rng(0)
                expression = {g.id: float(rng.uniform(1, 100)) for g in self.model.genes}
            reference = reference_flux(self.model, template, gene_expression=expression)
            with self.model:
                self.model.reactions.get_by_id(ko).bounds = (0.0, 0.0)
                if method_label == "ROOM":
                    result = room(self.model, reference)
                else:
                    result = moma(self.model, reference, linear=method_label == "MOMA (L1)")
        except Exception as exc:
            self.comparison_summary.setText(f"Comparison failed: {html.escape(str(exc))}")
            return
        if result.status != "optimal":
            self.comparison_table.setRowCount(0)
            self.comparison_summary.setText(
                f"<b>{method_label}</b> vs <b>{template}</b>: knocking out "
                f"{html.escape(ko)} makes the model <b>infeasible</b> "
                "(e.g. it cannot meet maintenance demand)."
            )
            self.status_label.setText("Comparison: infeasible perturbation.")
            return
        rows = sorted(result.fluxes.items(), key=lambda kv: -abs(kv[1] - reference.get(kv[0])))
        changed = [(r, v) for r, v in rows if abs(v - reference.get(r)) > 1e-6]
        self.comparison_table.setRowCount(len(changed))
        for i, (rid, flux) in enumerate(changed):
            self.comparison_table.setItem(i, 0, QTableWidgetItem(rid))
            self.comparison_table.setItem(i, 1, QTableWidgetItem(f"{reference.get(rid):.4g}"))
            self.comparison_table.setItem(i, 2, QTableWidgetItem(f"{flux:.4g}"))
        self.comparison_summary.setText(
            f"<b>{method_label}</b> vs <b>{template}</b> template after knocking out "
            f"{html.escape(ko)}: status {result.status}, distance {result.distance:.4g}, "
            f"{len(changed)} reactions changed."
        )
        self.status_label.setText(f"Comparison complete ({method_label}, {template}).")

    def _build_simulation_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        media_row = QHBoxLayout()
        media_row.addWidget(QLabel("Medium:"))
        self.medium_combo = QComboBox()
        self.medium_combo.addItems(sorted(PRESET_MEDIA))
        self.medium_combo.setToolTip("Preset growth medium (sets exchange uptake bounds)")
        media_row.addWidget(self.medium_combo)
        apply_medium_btn = QPushButton("Apply medium")
        apply_medium_btn.clicked.connect(self.apply_selected_medium)
        media_row.addWidget(apply_medium_btn)
        media_row.addStretch(1)
        layout.addLayout(media_row)

        controls = QHBoxLayout()
        fba_btn = QPushButton("Run FBA")
        fba_btn.clicked.connect(self.run_fba)
        pfba_btn = QPushButton("Run pFBA")
        pfba_btn.clicked.connect(self.run_pfba)
        fva_btn = QPushButton("Run FVA")
        fva_btn.clicked.connect(self.run_fva)
        controls.addWidget(fba_btn)
        controls.addWidget(pfba_btn)
        controls.addWidget(fva_btn)
        controls.addWidget(QLabel("fraction:"))
        self.fva_fraction_spin = QDoubleSpinBox()
        self.fva_fraction_spin.setRange(0.0, 1.0)
        self.fva_fraction_spin.setSingleStep(0.1)
        self.fva_fraction_spin.setValue(0.9)
        self.fva_fraction_spin.setToolTip("FVA: fraction of the optimum objective to hold")
        controls.addWidget(self.fva_fraction_spin)
        self.objective_label = QLabel("Objective: —")
        self.objective_label.setFont(QFont("", 13, QFont.Bold))
        controls.addStretch(1)
        controls.addWidget(self.objective_label)
        layout.addLayout(controls)

        self.sim_table = QTableWidget(0, 3)
        self.sim_table.setHorizontalHeaderLabels(["Reaction", "Flux", "FVA range"])
        _sim_header = self.sim_table.horizontalHeader()
        _sim_header.setSectionResizeMode(0, QHeaderView.Stretch)
        _sim_header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        _sim_header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.sim_table.verticalHeader().setVisible(False)
        self.sim_table.setAlternatingRowColors(True)
        layout.addWidget(self.sim_table, 1)
        return tab

    def _build_revert_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Normalization-target prediction (revert metabolism)")
        form = QFormLayout(controls)
        self.method_combo = QComboBox()
        self.method_combo.addItems(["rmta", "mta", "mta_miqp"])
        self.perturbation_combo = QComboBox()
        self.perturbation_combo.addItems(["gene", "reaction"])
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.0, 1.0)
        self.alpha_spin.setSingleStep(0.1)
        self.alpha_spin.setValue(0.9)
        form.addRow("Method:", self.method_combo)
        form.addRow("Knockout level:", self.perturbation_combo)
        form.addRow("Transformation weight α:", self.alpha_spin)

        source_row = QHBoxLayout()
        source_btn = QPushButton("Load source CSV/TSV…")
        source_btn.clicked.connect(lambda: self.load_revert_expression("source"))
        self.revert_source_label = QLabel("not loaded")
        source_row.addWidget(source_btn)
        source_row.addWidget(self.revert_source_label, 1)
        form.addRow("Source expression:", source_row)

        target_row = QHBoxLayout()
        target_btn = QPushButton("Load target CSV/TSV…")
        target_btn.clicked.connect(lambda: self.load_revert_expression("target"))
        self.revert_target_label = QLabel("not loaded")
        target_row.addWidget(target_btn)
        target_row.addWidget(self.revert_target_label, 1)
        form.addRow("Target expression:", target_row)

        self.revert_run_btn = QPushButton("Run Revert")
        self.revert_run_btn.setEnabled(False)
        self.revert_run_btn.clicked.connect(self.run_loaded_revert)
        self.revert_run_btn.setToolTip("Load both source and target expression files to enable.")
        form.addRow("", self.revert_run_btn)
        layout.addWidget(controls)

        self.revert_summary = QLabel(
            "Load source and target expression files to rank normalization targets."
        )
        self.revert_summary.setWordWrap(True)
        layout.addWidget(self.revert_summary)

        self.revert_table = QTableWidget(0, 3)
        self.revert_table.setHorizontalHeaderLabels(["Rank", "Target", "Score"])
        self.revert_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.revert_table.verticalHeader().setVisible(False)
        self.revert_table.setAlternatingRowColors(True)
        layout.addWidget(self.revert_table, 1)
        return tab

    # -- behavior -----------------------------------------------------------

    def load_model(self, model: Model) -> None:
        self.model = model
        solver = active_solver(model)
        self.subtitle.setText(
            f"{model.id}   ·   {len(model.reactions)} reactions   ·   "
            f"{len(model.metabolites)} metabolites   ·   {len(model.genes)} genes   ·   "
            f"solver: {solver}"
        )
        objective = ", ".join(
            r.id for r in model.reactions if r.objective_coefficient != 0
        )
        self.model_info.setText(
            f"<b>{model.id}</b><br>Objective: {objective or '—'} "
            f"({model.objective_direction})"
        )
        self._loading = True
        self.reaction_table.setRowCount(len(model.reactions))
        for row, rxn in enumerate(model.reactions):
            id_item = QTableWidgetItem(rxn.id)
            id_item.setFlags(id_item.flags() & ~Qt.ItemIsEditable)
            flux_item = QTableWidgetItem("—")
            flux_item.setFlags(flux_item.flags() & ~Qt.ItemIsEditable)
            lower_item = QTableWidgetItem(f"{rxn.lower_bound:g}")
            upper_item = QTableWidgetItem(f"{rxn.upper_bound:g}")
            self.reaction_table.setItem(row, 0, id_item)
            self.reaction_table.setItem(row, 1, lower_item)
            self.reaction_table.setItem(row, 2, upper_item)
            self.reaction_table.setItem(row, 3, flux_item)
        self._loading = False

        # Populate the production target selector with exchange reactions.
        self._loading = True
        self.product_combo.clear()
        exchanges = sorted(r.id for r in model.exchanges)
        self.product_combo.addItems(exchanges)
        self.substrate_combo.clear()
        self.substrate_combo.addItems(["auto", *exchanges])
        self.ko_combo.clear()
        self.ko_combo.addItems([r.id for r in model.reactions])
        self._loading = False
        if self._default_product and self._default_product in exchanges:
            self.product_combo.setCurrentText(self._default_product)

        # Production analysis needs exchange reactions; disable it otherwise.
        has_exchanges = bool(exchanges)
        for btn in getattr(self, "_production_buttons", []):
            btn.setEnabled(has_exchanges)
        if not has_exchanges:
            self.yield_label.setText("No exchange reactions in this model — production design unavailable.")

        self._fluxes = {}
        self._fluxes_stale = False
        self.status_label.setText(f"Loaded model '{model.id}'.")

    def set_reaction_bounds(
        self, reaction_id: str, lower: float | None = None, upper: float | None = None
    ) -> None:
        """Programmatically update a reaction bound and reflect it in the table."""

        rxn = self.model.reactions.get_by_id(reaction_id)
        # Assign both bounds at once so a value that crosses the current opposite bound does
        # not trip cobra's per-assignment validation (lower <= upper).
        new_lower = float(lower) if lower is not None else rxn.lower_bound
        new_upper = float(upper) if upper is not None else rxn.upper_bound
        if new_lower > new_upper:
            new_lower, new_upper = new_upper, new_lower
        rxn.bounds = (new_lower, new_upper)
        self._mark_fluxes_stale()
        for row in range(self.reaction_table.rowCount()):
            if self.reaction_table.item(row, 0).text() == reaction_id:
                self._loading = True
                self.reaction_table.item(row, 1).setText(f"{rxn.lower_bound:g}")
                self.reaction_table.item(row, 2).setText(f"{rxn.upper_bound:g}")
                self._loading = False
                break

    def _on_bound_edited(self, item) -> None:
        """Apply an in-table edit of a Lower/Upper cell to the model."""

        if self._loading or item.column() not in (1, 2):
            return
        row = item.row()
        reaction_id = self.reaction_table.item(row, 0).text()
        rxn = self.model.reactions.get_by_id(reaction_id)
        try:
            value = float(item.text())
        except ValueError:
            self._loading = True
            current = rxn.lower_bound if item.column() == 1 else rxn.upper_bound
            item.setText(f"{current:g}")
            self._loading = False
            self.status_label.setText(f"Invalid bound for {reaction_id}; reverted.")
            return
        if item.column() == 1:
            clamped = min(value, rxn.upper_bound)
            rxn.bounds = (clamped, rxn.upper_bound)
        else:
            clamped = max(value, rxn.lower_bound)
            rxn.bounds = (rxn.lower_bound, clamped)
        # Reflect the actually-stored (possibly clamped) value back into the cell.
        self._loading = True
        item.setText(f"{clamped:g}")
        self._loading = False
        note = "" if clamped == value else "  (clamped to keep lower <= upper)"
        self.status_label.setText(
            f"{reaction_id} bounds set to ({rxn.lower_bound:g}, {rxn.upper_bound:g}).{note}"
        )
        self._mark_fluxes_stale()

    def _mark_fluxes_stale(self) -> None:
        """Flag that bounds changed so the displayed flux column is out of date."""

        self._fluxes_stale = True
        if self._fluxes:
            self.status_label.setText("Bounds changed — fluxes are stale; click Run FBA.")

    def run_fba(self) -> None:
        try:
            solution = fba(self.model)
        except Exception as exc:  # surface solver/model errors instead of freezing
            self.objective_label.setText("Objective: error")
            self.status_label.setText(f"FBA failed: {exc}")
            return
        self._fluxes = dict(solution.fluxes)
        self._fluxes_stale = False
        obj = solution.objective_value
        obj_text = f"{obj:.4g}" if obj is not None else "infeasible"
        self.objective_label.setText(f"Objective: {obj_text} ({solution.status})")
        if obj is None or solution.status != "optimal":
            self.status_label.setText(f"FBA status: {solution.status} — model may be infeasible.")
        else:
            self.status_label.setText("FBA complete.")
        self._populate_flux_column()
        self._fill_sim_table()

    def run_pfba(self) -> None:
        try:
            solution = pfba(self.model)
        except Exception as exc:
            self.objective_label.setText("Objective: error")
            self.status_label.setText(f"pFBA failed: {exc}")
            return
        self._fluxes = dict(solution.fluxes)
        self._fluxes_stale = False
        # pFBA's objective value is the minimal total flux; show the growth (objective rxn).
        growth = next(
            (self._fluxes[r.id] for r in self.model.reactions if r.objective_coefficient != 0),
            None,
        )
        growth_text = f"{growth:.4g}" if growth is not None else "—"
        total = sum(abs(v) for v in self._fluxes.values())
        self.objective_label.setText(f"Objective: {growth_text} (pFBA)")
        self.status_label.setText(f"pFBA complete (minimal total flux = {total:.1f}).")
        self._populate_flux_column()
        self._fill_sim_table()

    def apply_selected_medium(self) -> None:
        name = self.medium_combo.currentText()
        try:
            applied = apply_medium(self.model, name)
        except Exception as exc:
            self.status_label.setText(f"Could not apply medium: {exc}")
            return
        # Reflect the new exchange bounds in the reaction table and mark fluxes stale.
        self.load_model(self.model)
        self.status_label.setText(
            f"Applied medium '{name}' ({len(applied)} open exchanges). Re-run FBA."
        )

    def run_fva(self) -> None:
        if not self._fluxes or self._fluxes_stale:
            self.run_fba()
        fraction = self.fva_fraction_spin.value()
        try:
            ranges = fva(self.model, fraction_of_optimum=fraction)
        except Exception as exc:
            self.status_label.setText(f"FVA failed: {exc}")
            return
        self.sim_table.setHorizontalHeaderLabels(["Reaction", "Flux", f"FVA range (f={fraction:g})"])
        self.sim_table.setRowCount(len(self._fluxes))
        for row, (rid, flux) in enumerate(sorted(self._fluxes.items())):
            self.sim_table.setItem(row, 0, QTableWidgetItem(rid))
            self.sim_table.setItem(row, 1, QTableWidgetItem(f"{flux:.3g}"))
            rng = ranges.get(rid)
            text = f"[{rng.minimum:.3g}, {rng.maximum:.3g}]" if rng else "—"
            self.sim_table.setItem(row, 2, QTableWidgetItem(text))
        self.status_label.setText(f"FVA complete (fraction {fraction:g}).")

    def _fill_sim_table(self) -> None:
        self.sim_table.setRowCount(len(self._fluxes))
        for row, (rid, flux) in enumerate(sorted(self._fluxes.items())):
            self.sim_table.setItem(row, 0, QTableWidgetItem(rid))
            self.sim_table.setItem(row, 1, QTableWidgetItem(f"{flux:.3g}"))
            self.sim_table.setItem(row, 2, QTableWidgetItem("—"))

    def _populate_flux_column(self) -> None:
        for row in range(self.reaction_table.rowCount()):
            rid = self.reaction_table.item(row, 0).text()
            flux = self._fluxes.get(rid, 0.0)
            self.reaction_table.setItem(row, 3, QTableWidgetItem(f"{flux:.3g}"))
        self._apply_flux_threshold(self.flux_slider.value())

    def _apply_flux_threshold(self, value: int) -> None:
        self.flux_threshold_label.setText(f"{value}%")
        if not self._fluxes:
            return
        max_abs = max((abs(v) for v in self._fluxes.values()), default=0.0) or 1.0
        threshold = (value / 100.0) * max_abs
        # Recolouring emits itemChanged; guard so it is not mistaken for a bound edit.
        self._loading = True
        for row in range(self.reaction_table.rowCount()):
            rid = self.reaction_table.item(row, 0).text()
            active = abs(self._fluxes.get(rid, 0.0)) >= threshold and threshold >= 0
            color = QColor("#d7f0dd") if active and value > 0 else QColor("white")
            for col in range(4):
                item = self.reaction_table.item(row, col)
                if item is not None:
                    item.setBackground(color)
        self._loading = False

    # -- production design --------------------------------------------------

    def _current_product(self) -> str:
        return self.product_combo.currentText()

    def _current_substrate(self) -> str | None:
        text = self.substrate_combo.currentText()
        return None if text in ("", "auto") else text

    def _is_aerobic(self) -> bool:
        return self.anaerobic_combo.currentText() == "aerobic"

    def render_flux_map(self) -> None:
        """Render the Escher-layout flux map coloured by the current FBA flux."""

        if self._map_path is None:
            return
        if not self._fluxes or self._fluxes_stale:
            self.run_fba()
        try:
            fig = escher_flux_map(self._map_path, self._fluxes,
                                  title=f"{self.model.id} — flux map")
        except Exception as exc:
            self.status_label.setText(f"Flux map failed: {exc}")
            return
        if self._map_canvas is not None:
            old = self._map_canvas.figure
            self.map_canvas_holder.removeWidget(self._map_canvas)
            self._map_canvas.setParent(None)
            old.clf()
        fig.set_dpi(110)
        self._map_canvas = FigureCanvas(fig)
        self._map_canvas.draw()
        self.map_canvas_holder.addWidget(self._map_canvas)
        self.status_label.setText("Flux map rendered (Escher layout).")

    def _on_product_changed(self, product: str) -> None:
        if not self._loading and product:
            self.yield_label.setText(
                f"Compute the theoretical yield, envelope, or FSEOF targets for "
                f"{html.escape(product)}."
            )

    def _set_production_figure(self, fig) -> None:
        if self._production_canvas is not None:
            old_fig = self._production_canvas.figure
            self.production_canvas_holder.removeWidget(self._production_canvas)
            self._production_canvas.setParent(None)
            old_fig.clf()  # release the previous figure so it is not leaked
        # Figures are authored at 300 DPI for file export; render on screen at a lower DPI
        # so the canvas fits the panel instead of clipping.
        fig.set_dpi(100)
        self._production_canvas = FigureCanvas(fig)
        self._production_canvas.draw()
        self.production_canvas_holder.addWidget(self._production_canvas)

    def _run_production(self, action) -> None:
        """Run a production analysis, surfacing any error instead of crashing the UI."""

        try:
            action()
        except Exception as exc:
            self.yield_label.setText(f"Production analysis failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Production analysis failed: {exc}")

    def run_theoretical_yield(self) -> None:
        def _do():
            product = self._current_product()
            result = theoretical_yield(self.model, product,
                                       substrate=self._current_substrate(),
                                       aerobic=self._is_aerobic())
            ceiling = (
                f", carbon ceiling {result.carbon_ceiling:.2f}"
                if result.carbon_ceiling is not None else ""
            )
            co2 = (
                f"; needs net CO₂ fixation ({result.co2_exchange:.1f})"
                if result.co2_fixed and result.exceeds_carbon_ceiling else ""
            )
            self.yield_label.setText(
                f"Theoretical yield of {html.escape(product)}: "
                f"<b>{result.molar_yield:.3f}</b> mol/mol {html.escape(result.substrate)} "
                f"({'aerobic' if result.aerobic else 'anaerobic'}; "
                f"max flux {result.product_flux:.2f}, uptake {result.substrate_uptake:.1f}{ceiling}{co2})."
            )
            self.status_label.setText("Theoretical yield computed.")

        self._run_production(_do)

    def run_production_envelope_plot(self) -> None:
        def _do():
            product = self._current_product()
            aerobic = self._is_aerobic()
            condition = "aerobic" if aerobic else "anaerobic"
            envelope = production_envelope(self.model, product,
                                           substrate=self._current_substrate(),
                                           aerobic=aerobic, points=20)
            fig = production_envelope_figure(
                envelope, title=f"Production envelope — {product} ({condition})"
            )
            self._set_production_figure(fig)
            self.status_label.setText(f"Production envelope computed ({condition}).")

        self._run_production(_do)

    def run_fseof_plot(self) -> None:
        def _do():
            product = self._current_product()
            result = fseof(self.model, product, n_steps=10, aerobic=self._is_aerobic())
            if result.metadata.get("max_product", 0.0) <= 1e-9:
                self.yield_label.setText(
                    f"FSEOF not meaningful: theoretical yield of {html.escape(product)} is zero "
                    "— check the exchange is open and the product is reachable."
                )
                self.status_label.setText("FSEOF: zero theoretical yield.")
                return
            fig = fseof_figure(result, top_n=6, title=f"FSEOF targets — {product}")
            self._set_production_figure(fig)
            # Rank the listed targets the same way the figure does (by flux increase).
            levels = list(result.enforced_levels)
            ranked = sorted(
                result.amplification_targets(),
                key=lambda rid: abs(result.trends.loc[rid, levels[-1]])
                - abs(result.trends.loc[rid, levels[0]]),
                reverse=True,
            )[:6]
            self.yield_label.setText(
                "FSEOF amplification targets: " + html.escape(", ".join(ranked))
            )
            self.status_label.setText("FSEOF complete.")

        self._run_production(_do)

    def run_fvseof_plot(self) -> None:
        def _do():
            product = self._current_product()
            result = fvseof(self.model, product, n_steps=8, aerobic=self._is_aerobic())
            if result.metadata.get("max_product", 0.0) <= 1e-9:
                self.yield_label.setText(
                    f"FVSEOF not meaningful: theoretical yield of {html.escape(product)} is zero."
                )
                self.status_label.setText("FVSEOF: zero theoretical yield.")
                return
            fig = fvseof_figure(result, top_n=5, title=f"FVSEOF robust targets — {product}")
            self._set_production_figure(fig)
            robust = result.robust_targets()[:8]
            self.yield_label.setText(
                "FVSEOF robust amplification targets (forced up): "
                + html.escape(", ".join(robust) if robust else "none")
            )
            self.status_label.setText("FVSEOF complete.")

        self._run_production(_do)

    # -- omics integration --------------------------------------------------

    def run_omics(self, gene_expression: Mapping[str, float]) -> None:
        """Predict a flux distribution from gene expression (E-Flux2 or LAD)."""

        method = self.omics_method_combo.currentText()
        try:
            weights = gene_to_reaction_weights(self.model, dict(gene_expression))
            result = integrate_expression(self.model, dict(gene_expression), method=method)
        except Exception as exc:
            self.omics_summary.setText(f"Omics integration failed: {html.escape(str(exc))}")
            return
        if result.status != "optimal":
            self.omics_summary.setText(
                f"{method.upper()} returned status '{result.status}' "
                f"({len(weights)} reactions mapped from expression)."
            )
            self.omics_table.setRowCount(0)
            return
        ranked = sorted(result.fluxes.items(), key=lambda kv: -abs(kv[1]))
        active = [(rid, v) for rid, v in ranked if abs(v) > 1e-6]
        self.omics_table.setRowCount(len(active))
        for i, (rid, flux) in enumerate(active):
            self.omics_table.setItem(i, 0, QTableWidgetItem(rid))
            self.omics_table.setItem(i, 1, QTableWidgetItem(f"{flux:.4g}"))
        detail = f" [{result.detail}]" if result.detail else ""
        self.omics_summary.setText(
            f"<b>{method.upper()}{detail}</b>: {len(weights)} reactions mapped from "
            f"{len(gene_expression)} genes; {len(active)} active fluxes. "
            f"Objective (FBA optimum) = {result.objective_value:.4g}."
        )
        self.status_label.setText(f"Omics integration complete ({method}).")

    def run_omics_demo(self) -> None:
        """Run omics integration on deterministic demo expression over the model's genes."""

        import numpy as np

        rng = np.random.default_rng(0)
        expression = {g.id: float(rng.uniform(1.0, 100.0)) for g in self.model.genes}
        if not expression:
            self.omics_summary.setText("Model has no genes — omics integration is not applicable.")
            return
        self.run_omics(expression)

    def load_expression_csv(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load gene-expression CSV", "", "CSV files (*.csv *.tsv);;All files (*)"
        )
        if not path:
            return
        try:
            expression = _read_expression_vector(path)
        except Exception as exc:
            self.omics_summary.setText(f"Could not read expression CSV: {html.escape(str(exc))}")
            return
        self.run_omics(expression)

    def _update_revert_run_state(self) -> None:
        ready = (
            self._revert_source_expression is not None
            and self._revert_target_expression is not None
        )
        self.revert_run_btn.setEnabled(ready)

    def load_revert_expression(self, role: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, f"Load {role} expression", "", "CSV/TSV files (*.csv *.tsv);;All files (*)"
        )
        if not path:
            return
        try:
            expression = _read_expression_vector(path)
        except Exception as exc:
            self.revert_summary.setText(
                f"Could not read {html.escape(role)} expression: {html.escape(str(exc))}"
            )
            return
        if role == "source":
            self._revert_source_expression = expression
            self.revert_source_label.setText(f"{len(expression)} genes")
        else:
            self._revert_target_expression = expression
            self.revert_target_label.setText(f"{len(expression)} genes")
        self._update_revert_run_state()
        self.status_label.setText(f"Loaded {role} expression ({len(expression)} genes).")

    def run_loaded_revert(self) -> None:
        if self._revert_source_expression is None or self._revert_target_expression is None:
            self.revert_summary.setText("Load both source and target expression files first.")
            self._update_revert_run_state()
            return
        self.run_revert(self._revert_source_expression, self._revert_target_expression)

    def run_revert(
        self,
        source_expression: Mapping[str, float],
        target_expression: Mapping[str, float],
    ) -> None:
        method = self.method_combo.currentText()
        perturbation = self.perturbation_combo.currentText()
        alpha = self.alpha_spin.value()
        try:
            reference = reference_state_pfba(self.model, name="source")
            direction = differential_expression(
                self.model, source_expression, target_expression, reference=reference
            )
            ranking = revert_targets(
                self.model,
                None,
                reference,
                direction,
                method=method,
                alpha=alpha,
                perturbation=perturbation,
            )
        except Exception as exc:  # surface solver-capability/model errors instead of freezing
            self.revert_summary.setText(f"Revert prediction failed: {exc}")
            self.status_label.setText(f"Revert prediction failed ({method}).")
            return
        rows = ranking.to_records()
        self.revert_table.setRowCount(len(rows))
        for i, record in enumerate(rows):
            self.revert_table.setItem(i, 0, QTableWidgetItem(str(record["rank"])))
            self.revert_table.setItem(i, 1, QTableWidgetItem(record["target_id"]))
            self.revert_table.setItem(i, 2, QTableWidgetItem(f"{record['score']:.4g}"))
            if i == 0:
                for col in range(3):
                    self.revert_table.item(i, col).setBackground(QColor("#d7f0dd"))
        best = ranking.best()
        if best is not None:
            self.revert_summary.setText(
                f"<b>Top normalization target: {best.target_id}</b> "
                f"(score {best.score:.4g}, method {method}, {perturbation} knockout). "
                f"Knocking it out moves the source state most toward the target state."
            )
        self.status_label.setText(f"Revert prediction complete ({method}).")
