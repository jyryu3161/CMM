"""CMM desktop main window.

A thin Qt view over the solver-neutral services. All numerics live in cmm.core / cmm.features
/ cmm.omics; this module only renders results and wires controls, so the GUI inherits none of
the analysis logic.
"""

from __future__ import annotations

import html
import math
import sys
import warnings
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING
from pathlib import Path

from cobra import Model
from cobra.io import from_json, to_json
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from qtpy.QtCore import QEvent, QEventLoop, QObject, Qt, QThread, QTimer, Signal
from qtpy.QtGui import QColor, QFont, QFontMetrics
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QSlider,
    QFrame,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QStyledItemDelegate,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from cmm.core import (
    PRESET_MEDIA,
    Condition,
    ReactionBound,
    active_solver,
    apply_medium,
    fba,
    fva,
    pfba,
    solver_status,
)
from cmm.features._perturbation import (
    blocked_reactions_for_genes,
    gene_perturbations,
    reaction_perturbations,
)
from cmm.features.comparison import (
    ROOM_TOLERANCES,
    batch_comparison,
    knockout_comparison,
    reference_flux,
)
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.response import flux_response
from cmm.features.revert import DEFAULT_EPSILON as REVERT_DEFAULT_EPSILON
from cmm.features.revert import revert_targets
from cmm.features.sampling import (
    random_flux_sampling,
    reference_constrained_sampling,
)
from cmm.features.strain_design import optknock, robustknock
from cmm.features.transformation import transformation_targets
from cmm.omics.conditions import (
    predict_condition_fluxes,
    read_expression_table,
)
from cmm.omics.differential import (
    gene_directions_by_fold_change,
    gene_directions_from_replicates,
    reaction_directions,
    restrict_to_top_changed,
)
from cmm.omics.expression import integrate_expression
from cmm.app.svg_background import svg_background
from cmm.app.jev_tab import JEV_TAB_NAME, JevTabMixin
from cmm.resources import bundled_map_for, map_reaction_ids

if TYPE_CHECKING:  # pragma: no cover - typing only
    import numpy as np
    import pandas as pd

from cmm.visualization import (
    escher_flux_map,
    flux_response_figure,
    fseof_figure,
    fvseof_figure,
    network_flux_map,
    production_envelope_figure,
    sampling_figure,
)

_ASSETS_DIR = Path(__file__).resolve().parent / "assets"
_APP_TITLE = "CMM — Constraint-based Metabolic Modeling"

# "Oxford Slate" — an academic / scientific-software palette: a deep-navy header with a
# restrained steel-blue accent over a cool-paper canvas. Every interactive widget (menus,
# combos, spin boxes, the flux slider, scrollbars, the status strip, tooltips) is styled so
# the platform reads as one cohesive professional tool rather than default Qt chrome.
#
#   navy  #1b2b44   accent #2f5e8f (hover #274f79 / pressed #1f3f61)
#   canvas #eef1f5  panels #ffffff  text #1a2433  muted #5c6b7e  border #d2dae4
_STYLE = """
QMainWindow, QWidget {
    background: #eef1f5; color: #1a2433; font-size: 12px;
    font-family: "Helvetica Neue", "Segoe UI", "Arial", sans-serif;
}

/* Title band */
QWidget#header { background: #1b2b44; border-bottom: 2px solid #2f5e8f; }
QLabel { background: transparent; color: #1a2433; }
QLabel#title { color: #ffffff; font-size: 18px; font-weight: 700; }
QLabel#subtitle { color: #9fb3cc; font-size: 12px; }
QLabel#statusbar { background: #e3e9f1; color: #4a5a70;
                   border-top: 1px solid #cdd6e1; padding: 5px 12px; font-size: 11px; }

/* Menu bar */
QMenuBar { background: #e8edf3; color: #27374b; border-bottom: 1px solid #d4dce6; padding: 2px 6px; }
QMenuBar::item { background: transparent; padding: 5px 12px; border-radius: 4px; }
QMenuBar::item:selected { background: #d6e3f2; color: #16324f; }
QMenu { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1; padding: 4px; }
QMenu::item { padding: 6px 26px 6px 18px; border-radius: 4px; }
QMenu::item:selected { background: #2f5e8f; color: #ffffff; }
QMenu::item:disabled { color: #8a98a8; }
QMenu::separator { height: 1px; background: #dde3ea; margin: 4px 8px; }

/* Cards */
QGroupBox { font-weight: 600; color: #27374b; border: 1px solid #d2dae4; border-radius: 8px;
            margin-top: 12px; padding: 12px 10px 10px 10px; background: #ffffff; }
QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 6px; color: #2f5e8f; }

/* Buttons */
QPushButton { background: #2f5e8f; color: #ffffff; border: 1px solid #2a557f;
              border-radius: 6px; padding: 6px 14px; font-weight: 600; }
QPushButton:hover { background: #274f79; }
QPushButton:pressed { background: #1f3f61; }
QPushButton:disabled { background: #c2cedd; color: #eef1f5; border-color: #c2cedd; }
QPushButton:focus { border: 1px solid #7da3c8; }

/* Combo boxes */
/* A combo and a spin box sit in the same form columns, so their total heights are pinned to
   the same 26px: min-height + both paddings + both 1px borders. Left to their own padding
   they came out 32px and 30px, which put every stacked group of stage boxes on a slightly
   different vertical grid. Change one of these and change the other to match. */
QComboBox { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
            border-radius: 6px; padding: 3px 8px; min-height: 18px; }
QComboBox:hover { border-color: #9bb4d0; }
QComboBox:focus { border-color: #2f5e8f; }
QComboBox::drop-down { border: none; width: 22px; }
QComboBox::down-arrow { image: url(__ASSETS__/down_arrow.svg);
    width: 10px; height: 7px; margin-right: 8px; }
QComboBox QAbstractItemView { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
    selection-background-color: #2f5e8f; selection-color: #ffffff; outline: none; }
QComboBox:disabled { background: #eef1f5; color: #9aa7b6; border-color: #dde3ea; }

/* Spin boxes */
QDoubleSpinBox, QSpinBox { background: #ffffff; color: #1a2433; border: 1px solid #cdd6e1;
    border-radius: 6px; padding: 3px 6px; min-height: 18px; }
QDoubleSpinBox:focus, QSpinBox:focus { border-color: #2f5e8f; }
/* Styling a spin box in QSS drops the native step arrows, so draw them explicitly. */
QDoubleSpinBox::up-button, QSpinBox::up-button { subcontrol-origin: border;
    subcontrol-position: top right; width: 16px; border-left: 1px solid #cdd6e1; }
QDoubleSpinBox::down-button, QSpinBox::down-button { subcontrol-origin: border;
    subcontrol-position: bottom right; width: 16px; border-left: 1px solid #cdd6e1; }
QDoubleSpinBox::up-button:hover, QSpinBox::up-button:hover,
QDoubleSpinBox::down-button:hover, QSpinBox::down-button:hover { background: #eef2f7; }
QDoubleSpinBox::up-arrow, QSpinBox::up-arrow {
    image: url(__ASSETS__/up_arrow.svg); width: 9px; height: 6px; }
QDoubleSpinBox::down-arrow, QSpinBox::down-arrow {
    image: url(__ASSETS__/down_arrow.svg); width: 9px; height: 6px; }
/* An inert control must read as inert: without these a disabled combo or spin box is almost
   indistinguishable from an active one, so users try to edit it and wonder why it does nothing. */
QDoubleSpinBox:disabled, QSpinBox:disabled { background: #eef1f5; color: #9aa7b6;
    border-color: #dde3ea; }
QDoubleSpinBox::up-button:disabled, QSpinBox::up-button:disabled,
QDoubleSpinBox::down-button:disabled, QSpinBox::down-button:disabled {
    border-left-color: #dde3ea; }
QLabel:disabled { color: #9aa7b6; }
QCheckBox:disabled { color: #9aa7b6; }

/* Flux-range slider */
QSlider::groove:horizontal { height: 5px; background: #d2dae4; border-radius: 3px; }
QSlider::sub-page:horizontal { background: #2f5e8f; border-radius: 3px; }
QSlider::handle:horizontal { background: #ffffff; border: 2px solid #2f5e8f;
    width: 14px; height: 14px; margin: -6px 0; border-radius: 8px; }
QSlider::handle:horizontal:hover { border-color: #274f79; }

/* Tables */
QTableWidget { background: #ffffff; alternate-background-color: #f5f8fb; color: #1a2433;
    gridline-color: #e4e9f0; border: 1px solid #d2dae4; border-radius: 6px;
    selection-background-color: #d6e3f2; selection-color: #16324f; }
QTableWidget::item { padding: 3px 6px; }
QHeaderView::section { background: #e8edf3; color: #34465c; padding: 6px 8px; border: none;
    border-right: 1px solid #d8dfe8; border-bottom: 1px solid #d2dae4; font-weight: 600; }

/* Tabs */
QTabWidget::pane { border: 1px solid #d2dae4; border-radius: 6px; top: -1px; background: #ffffff; }
/* No font-weight here: Qt's tabSizeHint ignores a stylesheet font-weight, so a bold label is
   drawn wider than the tab Qt reserved for it and clips. The tab-bar font (incl. its DemiBold
   weight) is set on the widget in code so the size hint measures the real text; the selected
   tab is distinguished by colour + underline. */
QTabBar::tab { background: #e6ebf2; color: #54637a; padding: 7px 13px; margin-right: 2px;
    border: 1px solid #d8dfe8; border-bottom: none;
    border-top-left-radius: 6px; border-top-right-radius: 6px; }
QTabBar::tab:hover { background: #eef2f7; color: #2f5e8f; }
QTabBar::tab:selected { background: #ffffff; color: #1b2b44;
    border-bottom: 2px solid #2f5e8f; margin-bottom: -1px; }

/* Scrollbars */
QScrollBar:vertical { background: transparent; width: 11px; margin: 0; }
QScrollBar::handle:vertical { background: #c2cedd; border-radius: 5px; min-height: 28px; }
QScrollBar::handle:vertical:hover { background: #a9b9cc; }
QScrollBar:horizontal { background: transparent; height: 11px; margin: 0; }
QScrollBar::handle:horizontal { background: #c2cedd; border-radius: 5px; min-width: 28px; }
QScrollBar::handle:horizontal:hover { background: #a9b9cc; }
QScrollBar::add-line, QScrollBar::sub-line { height: 0; width: 0; }
QScrollBar::add-page, QScrollBar::sub-page { background: transparent; }

/* Splitter + tooltips */
QSplitter::handle { background: #dde3ea; }
QSplitter::handle:horizontal { width: 6px; }
QSplitter::handle:hover { background: #c2cedd; }
QToolTip { background: #1b2b44; color: #ffffff; border: 1px solid #12203a;
    padding: 5px 8px; border-radius: 4px; }
"""


def _read_expression_vector(path: str) -> dict[str, float]:
    """Read a two-column gene-expression CSV/TSV file into ``gene -> value``."""

    import pandas as pd

    frame = pd.read_csv(path, sep=None, engine="python")
    if frame.shape[1] < 2:
        raise ValueError(
            "expression file must have at least two columns: gene, expression"
        )
    gene_col, value_col = frame.columns[:2]
    expression: dict[str, float] = {}
    for gene, value in zip(frame[gene_col], frame[value_col], strict=False):
        if pd.isna(gene) or pd.isna(value):
            continue
        gene_id = str(gene)
        if gene_id in expression:
            raise ValueError(f"expression file contains duplicate gene {gene_id!r}")
        numeric = float(value)
        if not math.isfinite(numeric) or numeric < 0:
            raise ValueError(
                f"expression for gene {gene_id!r} must be finite and non-negative"
            )
        expression[gene_id] = numeric
    if not expression:
        raise ValueError("expression file contains no numeric expression values")
    return expression


#: Added before the log2 of a linear expression value so a measured zero stays finite. Matches
#: :func:`cmm.omics.differential.gene_log2_fold_change`'s default, which is what keeps a
#: one-column linear file scoring identically through the old and the new direction routes.
_LOG2_PSEUDOCOUNT = 1.0

#: Widest a stage box's input is allowed to grow. The fields are told to fill their column so
#: that every stage starts at the same x, but the column is as wide as the left panel, and a
#: spin box holding "0.66" stretched across 600px reads as a text field, not a number entry.
#: Below this the fields still fill the column, so they stay equal at any window size.
_FIELD_MAX_WIDTH = 260

#: Slack added to a stage label's measured text width before the shared column is pinned to it.
#: A QLabel draws its text a pixel or two inside its own rect, and the column is fixed rather
#: than merely floored, so an exact fit would clip the widest label on some platforms.
_LABEL_PAD = 4


def _read_expression_matrix(path: str) -> "pd.DataFrame":
    """Read a gene-expression file into a gene-indexed frame of one column per measurement.

    One numeric column is a single measurement per gene and supports a fold-change cut only.
    Two or more are replicates, which is what Yizhak et al.'s (2013) Student's t-test needs.
    The shape is a fact about the file, so it is read and reported rather than asked about.

    Values are linear expression — counts, TPM, FPKM — as everywhere else in the app, so a
    negative one is rejected here rather than left to surface as a solver error later.
    """

    import numpy as np
    import pandas as pd

    frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    numeric = frame.select_dtypes("number")
    if numeric.empty:
        raise ValueError(
            "expression file must have a gene id column followed by at least one numeric "
            "expression column"
        )
    if numeric.isna().any().any():
        raise ValueError("expression file carries missing values")
    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("expression values must be finite and non-negative")
    numeric.index = [str(value) for value in numeric.index]
    if numeric.index.has_duplicates:
        raise ValueError("expression file contains duplicate gene ids")
    return numeric


def _as_expression_matrix(expression) -> "pd.DataFrame":
    """Accept either a gene -> value mapping or an already-read replicate frame.

    ``run_revert`` is called directly by the screenshot script and the smoke tests with plain
    mappings, so widening the tab to replicates must not narrow that entry point.
    """

    import pandas as pd

    if isinstance(expression, pd.DataFrame):
        return expression
    series = pd.Series(
        {str(gene): float(value) for gene, value in expression.items()}, dtype=float
    )
    return series.to_frame("value")


def _expression_in_log2(matrix: "pd.DataFrame") -> "pd.DataFrame":
    """The replicate matrix in log2 space, which is where both direction tests operate.

    Fold change is multiplicative, so only in log space is a doubling and a halving the same
    distance either side of zero. The transform is applied per replicate rather than after
    averaging; for a single column the two are the same number, which is what makes this route
    agree exactly with the fold-change cut this tab applied before it accepted replicates.
    """

    import numpy as np

    return np.log2(matrix + _LOG2_PSEUDOCOUNT)


class _AnalysisWorker(QObject):
    """Runs one compute callable on a worker thread and reports back to the UI thread.

    The callable must be pure computation over the cobra model / services (no Qt or matplotlib
    access) — those stay on the main thread. The result (or the raised exception) is captured
    so the caller can render it (or re-raise) once ``finished`` fires.
    """

    finished = Signal()

    def __init__(self, fn, cleanup=None):
        super().__init__()
        self._fn = fn
        self._cleanup = cleanup
        self.result = None
        self.error: Exception | None = None

    def run(self) -> None:
        try:
            self.result = self._fn()
        except (
            Exception
        ) as exc:  # captured and re-raised on the UI thread by the caller
            self.error = exc
        finally:
            if self._cleanup is not None:
                try:
                    self._cleanup()
                except Exception as exc:
                    if self.error is None:
                        self.error = exc
            self.finished.emit()


class _ComboItemDelegate(QStyledItemDelegate):
    """Left-align combo popup entries and paint disabled ones as clearly inert.

    QSS ``::item:disabled`` styling is unreliable for combo popups (notably on macOS, where it
    does not render), so disabled rows are painted here: a muted neutral fill with faint text
    and no hover/selection highlight. Enabled rows defer to normal styled rendering.
    """

    _DISABLED_TEXT = QColor("#b3bcc7")
    _DISABLED_FILL = QColor("#edeef1")
    _ALIGN = int(Qt.AlignLeft | Qt.AlignVCenter)

    def initStyleOption(self, option, index):  # noqa: N802 (Qt override)
        super().initStyleOption(option, index)
        option.displayAlignment = Qt.AlignLeft | Qt.AlignVCenter

    def paint(self, painter, option, index):
        if not (index.flags() & Qt.ItemIsEnabled):
            painter.save()
            painter.fillRect(option.rect, self._DISABLED_FILL)
            painter.setPen(self._DISABLED_TEXT)
            text = index.data(Qt.DisplayRole)
            # Indent to match the styled left padding of enabled rows.
            rect = option.rect.adjusted(6, 0, -6, 0)
            painter.drawText(rect, self._ALIGN, "" if text is None else str(text))
            painter.restore()
            return
        super().paint(painter, option, index)


class _DisabledHoverCursor(QObject):
    """Show a not-allowed cursor over disabled entries in a combo-box popup."""

    def __init__(self, view):
        super().__init__(view)
        self._view = view

    def eventFilter(self, viewport, event):  # noqa: N802 (Qt override)
        if event.type() == QEvent.MouseMove:
            try:
                pos = event.position().toPoint()  # Qt6
            except AttributeError:
                pos = event.pos()  # Qt5
            index = self._view.indexAt(pos)
            enabled = (
                bool(self._view.model().flags(index) & Qt.ItemIsEnabled)
                if index.isValid()
                else True
            )
            viewport.setCursor(Qt.ArrowCursor if enabled else Qt.ForbiddenCursor)
        return False


_MAP_LAYOUT_ESCHER = "Escher map (curated)"
_MAP_LAYOUT_SCHEMATIC = "Schematic (top reactions)"


class CmmMainWindow(JevTabMixin, QMainWindow):
    """Main platform window over a single cobra model."""

    def __init__(
        self,
        model: Model,
        parent=None,
        *,
        default_product: str | None = None,
        map_path: str | None = None,
    ):
        super().__init__(parent)
        self.model = model
        self._fluxes: dict[
            str, float
        ] = {}  # the last-run distribution (drives left panel/FVA/map)
        # What produced `_fluxes`. Every view of a flux distribution has to be able to name
        # the method behind it — a pFBA map read as an FBA map is a wrong figure, not a
        # cosmetic one — so the label travels with the numbers.
        self._flux_source: str = ""
        self._fba_fluxes: dict[str, float] = {}  # FBA column of the simulation table
        self._pfba_fluxes: dict[str, float] = {}  # pFBA column of the simulation table
        self._loading = False
        self._fluxes_stale = False
        self._default_product = default_product
        # An explicit map wins; otherwise offer a bundled curated map if one suits this
        # model. Either way the Flux Map tab exists — without a map it draws the schematic.
        self._map_path = map_path or bundled_map_for(model)
        # An optional drawing of the same map, laid under the flux strokes. CMM's renderer
        # reuses the layout but not Escher's drawing conventions, so a picture exported from
        # the map itself is the way to get those without embedding a browser.
        self._map_background: "np.ndarray | None" = None
        self._revert_source_expression: "pd.DataFrame | None" = None
        self._revert_target_expression: "pd.DataFrame | None" = None
        # The tab falls back to a fold-change cut on data that cannot carry a t-test; these
        # separate that fallback from a fold-change run the user asked for.
        self._revert_significance_forced = False
        self._revert_significance_syncing = False
        self._stage_labels_aligned = False
        # Real gene expression loaded from a CSV on the Omics tab, reused as the LAD/E-Flux2
        # reference template on the Comparison tab so that path is not left to synthetic data.
        self._omics_expression: dict[str, float] | None = None
        # Which integration methods ("lad"/"eflux2") were actually computed on the Omics tab —
        # only those are offered as a Comparison reference template.
        self._omics_computed_methods: set[str] = set()
        # Last single-run MOMA/ROOM result, kept so the change-threshold can re-filter the
        # table without re-solving.
        self._comparison_cache: dict | None = None
        # Last production result (kind, result) so the FSEOF/FVSEOF "show all" toggle can
        # re-render the table without recomputing.
        self._production_result: tuple[str, object] | None = None
        # Keep hover-cursor event filters alive for combos that grey out entries.
        self._cursor_filters: list[_DisabledHoverCursor] = []
        # A→B transformation finder: two expression vectors → two predicted flux states.
        self._transform_source_expression: dict[str, float] | None = None
        self._transform_target_expression: dict[str, float] | None = None
        # Omics integration source: a gene × condition expression table (1+ conditions).
        # Loading only stores it; Compute predicts a flux column per selected condition.
        self._omics_table_df = None
        # Last computed per-condition fluxes (condition -> {reaction: flux}), kept so the
        # active-only / show-all toggle can re-render without recomputing.
        self._omics_fluxes_by_condition: dict[str, dict[str, float]] | None = None
        self._omics_method_label = ""
        self._omics_conditions_order: list[str] = []
        # Guards re-entrant background runs (see _run_in_background).
        self._busy = False
        self.setWindowTitle(_APP_TITLE)
        self.resize(1160, 760)
        # Qt QSS ``image: url(...)`` needs a real path; forward slashes work cross-platform.
        self.setStyleSheet(_STYLE.replace("__ASSETS__", _ASSETS_DIR.as_posix()))
        self._build_ui()
        self._style_all_combo_popups()
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
        self.status_label.setObjectName("statusbar")
        outer.addWidget(self.status_label)

        self._pin_form_layout_policies()

    def _build_menus(self) -> None:
        """A desktop-style menu bar: Analysis / Model / Config drive the platform."""

        bar = self.menuBar()

        file_menu = bar.addMenu("&File")
        file_menu.addAction("Open Model…", self.open_model_dialog)
        file_menu.addAction("Open Escher Map (JSON)…", self.open_flux_map_dialog)
        file_menu.addAction("Open Map Image…", self.open_map_background_dialog)
        file_menu.addSeparator()
        file_menu.addAction("Export Table to CSV…", self.export_table_csv)
        file_menu.addAction("Save Figure…", self.save_active_figure)
        file_menu.addSeparator()
        file_menu.addAction("Quit", self.close)

        analysis = bar.addMenu("&Analysis")
        analysis.addAction("Run FBA", self.run_fba)
        analysis.addAction("Run FVA", self.run_fva)
        analysis.addSeparator()
        analysis.addAction(
            "Theoretical Yield",
            lambda: self._in_tab("Production", self.run_theoretical_yield),
        )
        analysis.addAction(
            "Production Envelope",
            lambda: self._in_tab("Production", self.run_production_envelope_plot),
        )
        analysis.addAction(
            "FSEOF Targets", lambda: self._in_tab("Production", self.run_fseof_plot)
        )
        analysis.addAction(
            "Strain Design (OptKnock / RobustKnock)…",
            lambda: self._goto_tab("Strain Design"),
        )
        analysis.addSeparator()
        analysis.addAction(
            "Omics Integration (E-Flux2 / LAD)…", lambda: self._goto_tab("Omics")
        )
        analysis.addAction(
            "Revert Metabolism…", lambda: self._goto_tab("Revert Metabolism")
        )
        analysis.addAction(
            "Transformation Targets (A→B)…", lambda: self._goto_tab("Transform (A→B)")
        )
        analysis.addSeparator()
        analysis.addAction(
            "Render Flux Map",
            lambda: self._in_tab("Flux Map", self.render_flux_map),
        )

        jev = bar.addMenu("&Agent")
        jev.addAction("Let the Agent Play…", lambda: self._goto_tab(JEV_TAB_NAME))
        jev.addSeparator()
        jev.addAction("Show Run Progress", self.show_jev_progress)
        jev.addAction("Show Last Decision", self.show_jev_last_decision)
        jev.addSeparator()
        jev.addAction("Set OpenRouter API Key…", self.set_jev_api_key)
        jev.addAction("Clear OpenRouter API Key", self.clear_jev_api_key)
        jev.addAction("Where is my key?", self.show_jev_api_key_status)

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
            self,
            "Model Info",
            f"{self.model.id}\n\n"
            f"Reactions: {len(self.model.reactions)}\n"
            f"Metabolites: {len(self.model.metabolites)}\n"
            f"Genes: {len(self.model.genes)}\n"
            f"Objective: {objective or '—'} ({self.model.objective_direction})\n"
            f"Solver: {active_solver(self.model)}",
        )

    # -- background execution -----------------------------------------------

    def _run_in_background(self, compute, *, label: str = "Working…", cleanup=None):
        """Run ``compute()`` off the UI thread and return its result synchronously.

        The heavy solve runs on a worker thread while a nested event loop keeps the window
        repainting; a modal busy dialog blocks input so the model is not mutated mid-solve.
        Because the nested loop only exits when the worker finishes, callers still see a plain
        synchronous call (tests and the scenario harnesses read results immediately after).

        Re-entrant calls (a backgrounded handler invoking another) run inline so busy dialogs
        and worker threads never stack.
        """

        if self._busy:  # already inside a background run — don't nest threads/dialogs
            return compute()

        self._busy = True
        worker = _AnalysisWorker(compute, cleanup)
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)

        loop = QEventLoop()
        worker.finished.connect(loop.quit)

        # A modal, no-cancel busy dialog: blocks input (menus, table edits, other analyses) so
        # only the worker touches the model, and auto-shows only if the work exceeds
        # minimumDuration so fast solves do not flash a dialog.
        dialog = QProgressDialog(label, None, 0, 0, self)
        dialog.setWindowTitle("CMM")
        dialog.setWindowModality(Qt.ApplicationModal)
        dialog.setMinimumDuration(400)
        dialog.setAutoClose(False)
        dialog.setAutoReset(False)

        # Start the thread only once the event loop is running, so the worker cannot finish and
        # queue loop.quit before exec_() begins (which would otherwise hang the loop).
        QTimer.singleShot(0, thread.start)
        try:
            loop.exec_()
        finally:
            dialog.reset()
            dialog.deleteLater()
            thread.quit()
            thread.wait()
            worker.deleteLater()
            self._busy = False

        if worker.error is not None:
            raise worker.error
        return worker.result

    def _run_model_in_background(self, compute, *, label: str = "Working…"):
        """Run a solve against a solver object created inside the worker thread.

        GLPK's re-entrant build permits independent problem objects in different threads, but
        explicitly forbids a thread from accessing a problem object created by another thread.
        A ``cobra.Model`` owns such a solver object, so passing ``self.model`` directly to a
        worker is unsupported and can make GLPK abort the whole process (observed on Windows
        during FVSEOF's repeated FVA solves).

        Serialize the standard COBRA model state on the UI thread and reconstruct it in the
        worker. Preserve the active solver, objective direction, model tolerance, and public
        solver configuration so the detached solve has the same scientific inputs. The worker
        model is intentionally disposable: GUI edits remain on ``self.model`` and service
        context-manager mutations cannot leak back into it.
        """

        model_json = to_json(self.model)
        solver_name = active_solver(self.model)
        solver_interface = self.model.solver.interface
        solver_json = self.model.solver.to_json()
        default_solver_name = active_solver()
        objective_direction = self.model.objective_direction
        model_tolerance = self.model.tolerance
        configuration = self.model.solver.configuration
        solver_parameters = {
            name: getattr(configuration, name)
            for name in ("presolve", "timeout", "verbosity", "lp_method", "qp_method")
            if hasattr(configuration, name)
        }
        solver_tolerances = configuration.tolerances.to_dict()

        def _compute_with_worker_model():
            worker_model = from_json(model_json)
            glpk_problems = []
            if active_solver(worker_model) in {"glpk", "glpk_exact"}:
                glpk_problems.append(worker_model.solver.problem)
            # COBRA JSON carries the biological network but not custom optlang constraints.
            # Rebuild the complete solver formulation in this thread, then attach it to the
            # detached COBRA objects so programmatic models keep those scientific constraints.
            worker_model._solver = solver_interface.Model.from_json(solver_json)
            if solver_name in {"glpk", "glpk_exact"}:
                glpk_problems.append(worker_model.solver.problem)
            worker_model.objective_direction = objective_direction
            worker_model.tolerance = model_tolerance
            worker_configuration = worker_model.solver.configuration
            for name, value in solver_parameters.items():
                setattr(worker_configuration, name, value)
            for name, value in solver_tolerances.items():
                setattr(worker_configuration.tolerances, name, value)
            try:
                return compute(worker_model)
            finally:
                if sys.platform != "win32" and glpk_problems:
                    # The swiglpk Linux/macOS wheels explicitly disable GLPK re-entrancy, so
                    # glp_free_env would also destroy the UI model's process-wide environment.
                    # Release only the detached problems that this helper created.
                    from swiglpk import glp_delete_prob

                    for problem in glpk_problems:
                        glp_delete_prob(problem)

        cleanup = None
        if sys.platform == "win32" and {solver_name, default_solver_name} & {
            "glpk",
            "glpk_exact",
        }:
            # swiglpk's Windows wheel builds GLPK with TLS. Its re-entrancy contract requires
            # every worker to release that thread-local environment before termination; this
            # does not touch the UI thread's separate environment or original model.
            from swiglpk import glp_free_env

            cleanup = glp_free_env
        return self._run_in_background(
            _compute_with_worker_model, label=label, cleanup=cleanup
        )

    def _pin_form_layout_policies(self) -> None:
        """Lay every form out the same way regardless of the host platform's style.

        ``QFormLayout`` reads its alignment and growth policy from the style, and the styles
        disagree: most give ``AlignLeft|AlignTop`` with ``AllNonFixedFieldsGrow``, while macOS
        gives ``AlignHCenter|AlignTop`` with ``FieldsStayAtSizeHint``. Under the macOS pair
        each box centres its own label-plus-field block at that block's natural width, so a
        column of stage boxes lands at as many different x positions as it has boxes and no
        two inputs share a width - while the same code renders perfectly aligned under any
        other style. Pinning both makes the layout the design's, not the host's.
        """

        for form in self.findChildren(QFormLayout):
            form.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
            form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)

    def showEvent(self, event) -> None:
        """Square the revert tab's four label columns into one, once, on first show.

        Deferred to here rather than done while building: a label's font is only finally
        resolved when the window is realized, and a column pinned to a width measured before
        that is wrong on exactly the machines whose fonts differ from the developer's. The
        same pass caps how wide the inputs may grow.
        """

        super().showEvent(event)
        if not self._stage_labels_aligned:
            self._stage_labels_aligned = True
            self._align_stage_forms(self._revert_stage_forms)

    # -- File menu (open / export) ------------------------------------------

    def open_model_dialog(self) -> None:
        """Load a different SBML model into the running window."""

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open SBML model",
            "",
            "SBML models (*.xml *.sbml *.xml.gz *.json);;All files (*)",
        )
        if not path:
            return
        try:
            import cobra

            model = (
                cobra.io.load_json_model(path)
                if path.endswith(".json")
                else cobra.io.read_sbml_model(path)
            )
        except Exception as exc:
            self.status_label.setText(f"Could not load model: {html.escape(str(exc))}")
            return
        # A new model invalidates any loaded expression / revert state.
        self._omics_expression = None
        self._omics_computed_methods = set()
        self._refresh_comparison_templates()  # disable LAD/E-Flux2 until recomputed
        self._revert_source_expression = None
        self._revert_target_expression = None
        self._transform_source_expression = None
        self._transform_target_expression = None
        self._omics_table_df = None
        self._omics_fluxes_by_condition = None
        self.omics_map_btn.setEnabled(False)
        self.omics_map_combo.clear()
        if hasattr(self, "omics_map_label"):
            self._refresh_omics_map_label()
        self._omics_conditions_order = []
        self.revert_source_label.setText("not loaded")
        self.revert_target_label.setText("not loaded")
        self.transform_source_label.setText("not loaded")
        self.transform_target_label.setText("not loaded")
        self.omics_file_label.setText("no file loaded")
        self.omics_compute_btn.setEnabled(False)
        self.omics_cond_list.clear()
        self.omics_table.setRowCount(0)
        self.transform_table.setRowCount(0)
        self._update_revert_run_state()
        self._update_transform_run_state()
        self.load_model(model)
        self.status_label.setText(f"Loaded model '{model.id}' from {path}.")

    def _active_table(self) -> QTableWidget | None:
        """The primary result table of the current tab, or None for figure-only tabs."""

        name = self.tabs.tabText(self.tabs.currentIndex())
        return {
            "Simulation": self.sim_table,
            "Comparison": self.comparison_table,
            "Production": self.fseof_table,
            "Strain Design": self.sd_table,
            "Flux Response": self.fr_table,
            "Sampling": self.sample_table,
            "Omics": self.omics_table,
            "Revert Metabolism": self.revert_table,
            "Transform (A→B)": self.transform_table,
            # Whichever of the JEV tab's two tables is in front, so Export always exports
            # what the user is looking at.
            JEV_TAB_NAME: self.jev_lower_tabs.currentWidget(),
        }.get(name)

    def _active_figure(self):
        """The matplotlib figure of the current tab, or None if it has no figure."""

        name = self.tabs.tabText(self.tabs.currentIndex())
        if name == "Production" and self._production_canvas is not None:
            return self._production_canvas.figure
        if name == "Flux Response" and self._flux_response_canvas is not None:
            return self._flux_response_canvas.figure
        if name == "Sampling" and self._sampling_canvas is not None:
            return self._sampling_canvas.figure
        if name == "Flux Map" and self._map_canvas is not None:
            return self._map_canvas.figure
        if name == JEV_TAB_NAME and self._jev_canvas is not None:
            return self._jev_canvas.figure
        return None

    def export_table_csv(self) -> None:
        """Write the current tab's result table to a CSV file."""

        import csv

        table = self._active_table()
        if table is None or table.rowCount() == 0:
            self.status_label.setText(
                "No result table to export on this tab — run an analysis first."
            )
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export table to CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if not path:
            return
        headers = [
            table.horizontalHeaderItem(c).text()
            if table.horizontalHeaderItem(c)
            else f"col{c}"
            for c in range(table.columnCount())
        ]
        try:
            with open(path, "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(headers)
                for row in range(table.rowCount()):
                    writer.writerow(
                        [
                            table.item(row, col).text() if table.item(row, col) else ""
                            for col in range(table.columnCount())
                        ]
                    )
        except OSError as exc:
            self.status_label.setText(f"Could not write CSV: {exc}")
            return
        self.status_label.setText(f"Exported {table.rowCount()} rows to {path}.")

    def save_active_figure(self) -> None:
        """Save the current tab's figure (Production or Flux Map) at 300 DPI."""

        fig = self._active_figure()
        if fig is None:
            self.status_label.setText(
                "No figure to save — open the Production or Flux Map tab and run an analysis."
            )
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save figure", "", "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)"
        )
        if not path:
            return
        try:
            fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        except Exception as exc:
            self.status_label.setText(f"Could not save figure: {html.escape(str(exc))}")
            return
        self.status_label.setText(f"Saved figure to {path}.")

    def _build_header(self) -> QWidget:
        header = QWidget()
        header.setObjectName("header")
        layout = QVBoxLayout(header)
        layout.setContentsMargins(16, 12, 16, 12)
        title = QLabel(_APP_TITLE)
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
        self.reaction_table.setHorizontalHeaderLabels(
            ["Reaction", "Lower", "Upper", "Flux"]
        )
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
        hint = QLabel(
            "Double-click a Lower/Upper cell to edit a bound, then re-run FBA."
        )
        hint.setStyleSheet("color: #5c6b7e; font-style: italic;")
        layout.addWidget(hint)
        return box

    def _build_tabs(self) -> QWidget:
        """Tabs ordered so each analysis sits next to the one asking a similar question.

        Simulation/Sampling describe what the model can do with no intervention;
        Comparison/Flux Response predict the consequence of one (discrete and continuous
        respectively — a knockout is a response scan pinned at zero); Production/Strain Design
        propose interventions; the rest are omics-driven.
        """

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_simulation_tab(), "Simulation")
        self.tabs.addTab(self._build_sampling_tab(), "Sampling")
        self.tabs.addTab(self._build_comparison_tab(), "Comparison")
        self.tabs.addTab(self._build_flux_response_tab(), "Flux Response")
        self.tabs.addTab(self._build_production_tab(), "Production")
        self.tabs.addTab(self._build_strain_design_tab(), "Strain Design")
        self.tabs.addTab(self._build_omics_tab(), "Omics")
        self.tabs.addTab(self._build_fluxmap_tab(), "Flux Map")
        self.tabs.addTab(self._build_revert_tab(), "Revert Metabolism")
        self.tabs.addTab(self._build_transformation_tab(), "Transform (A→B)")
        # The agent tab is last: it drives the same services the tabs before it expose,
        # and it is the only one that needs a network credential.
        self.tabs.addTab(self._build_jev_tab(), JEV_TAB_NAME)
        return self.tabs

    def _build_omics_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Omics integration (expression → per-condition flux)")
        cbox = QVBoxLayout(controls)

        # Row 1: method + load (loading only stores the file; Compute runs the integration).
        row = QHBoxLayout()
        row.addWidget(QLabel("Method:"))
        self.omics_method_combo = QComboBox()
        self.omics_method_combo.addItems(["eflux2", "lad"])
        self.omics_method_combo.setToolTip(
            "E-Flux2: scale bounds by expression, minimize L2 flux. LAD: fit fluxes to expression."
        )
        row.addWidget(self.omics_method_combo)
        load_btn = QPushButton("Load expression CSV…")
        load_btn.clicked.connect(self.load_expression_csv)
        row.addWidget(load_btn)
        demo_btn = QPushButton("Use demo expression")
        demo_btn.clicked.connect(self.load_omics_demo)
        row.addWidget(demo_btn)
        row.addStretch(1)
        cbox.addLayout(row)

        # Row 2: the loaded filename and Compute (separate so the method can change and re-run
        # on the same loaded data).
        run_row = QHBoxLayout()
        self.omics_file_label = QLabel("no file loaded")
        self.omics_file_label.setStyleSheet("color: #5c6b7e; font-style: italic;")
        run_row.addWidget(self.omics_file_label, 1)
        self.omics_compute_btn = QPushButton("Compute")
        self.omics_compute_btn.setEnabled(False)
        self.omics_compute_btn.setToolTip(
            "Load an expression source first, then compute the checked conditions."
        )
        self.omics_compute_btn.clicked.connect(self.compute_omics)
        run_row.addWidget(self.omics_compute_btn)
        cbox.addLayout(run_row)

        # Row 2: the conditions detected in the file, each check-selectable.
        cond_row = QHBoxLayout()
        cond_row.addWidget(QLabel("Conditions:"))
        self.omics_cond_list = QListWidget()
        self.omics_cond_list.setSelectionMode(QAbstractItemView.NoSelection)
        self.omics_cond_list.setFlow(QListWidget.LeftToRight)
        self.omics_cond_list.setWrapping(True)
        self.omics_cond_list.setMaximumHeight(64)
        self.omics_cond_list.setToolTip(
            "Columns detected in the file. Check the conditions to predict a flux column for."
        )
        cond_row.addWidget(self.omics_cond_list, 1)
        cbox.addLayout(cond_row)

        # Row 3: send one computed condition to the flux map. An expression-derived flux state
        # is a distribution over the whole network; a 95-row table is the wrong way to read it.
        map_row = QHBoxLayout()
        # The label carries the count. A closed combo shows one item, and directly above it a
        # "Conditions:" checklist shows every condition — so a two-condition run reads as though
        # one went missing. Saying how many are in the box answers that before it is asked, and
        # "Draw" separates this row from the checklist that decides what gets computed.
        self.omics_map_label = QLabel("Draw:")
        map_row.addWidget(self.omics_map_label)
        self.omics_map_combo = QComboBox()
        self.omics_map_combo.setMinimumWidth(160)
        self.omics_map_combo.setToolTip("Which computed condition to draw")
        self.omics_map_combo.currentIndexChanged.connect(self._refresh_omics_map_label)
        map_row.addWidget(self.omics_map_combo)
        self.omics_map_btn = QPushButton("Show on flux map")
        self.omics_map_btn.setEnabled(False)
        self.omics_map_btn.setToolTip(
            "Compute a condition first, then draw it on the map."
        )
        self.omics_map_btn.clicked.connect(self.show_omics_on_flux_map)
        map_row.addWidget(self.omics_map_btn)
        map_row.addStretch(1)
        cbox.addLayout(map_row)
        layout.addWidget(controls)

        self.omics_summary = QLabel(
            "Load a gene × condition CSV (first column = gene id, remaining columns = "
            "conditions) or use demo data, then check conditions and click Compute."
        )
        self.omics_summary.setWordWrap(True)
        layout.addWidget(self.omics_summary)

        # Result-view controls sit with the table so the toggle reads as adjusting the results.
        table_controls = QHBoxLayout()
        table_controls.addWidget(QLabel("Predicted fluxes:"))
        table_controls.addStretch(1)
        self.omics_show_all_check = QCheckBox("Show all reactions")
        self.omics_show_all_check.setToolTip(
            "Off: only reactions with non-zero flux in some condition. On: every reaction."
        )
        self.omics_show_all_check.toggled.connect(self._render_omics_table)
        table_controls.addWidget(self.omics_show_all_check)
        layout.addLayout(table_controls)

        # Columns are set per run: Reaction + one predicted-flux column per computed condition.
        self.omics_table = QTableWidget(0, 1)
        self.omics_table.setHorizontalHeaderLabels(["Reaction"])
        self.omics_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.omics_table.verticalHeader().setVisible(False)
        self.omics_table.setAlternatingRowColors(True)
        layout.addWidget(self.omics_table, 1)
        return tab

    def _build_fluxmap_tab(self) -> QWidget:
        """The flux map: a curated Escher layout when one is available, else a schematic.

        CMM never invents a map layout. A readable metabolic map is hand-drawn, so the Escher
        mode reuses a published layout, and the schematic mode is the honest fallback for a
        model no bundled map fits — it draws only the highest-flux reactions and says so.
        """

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QHBoxLayout()

        # Solve and draw in one click. Without these the tab silently drew FBA and there was
        # nothing to suggest pFBA was equally available — it had to be run from another tab
        # first, which no one would guess. There is no separate "render" button: changing the
        # layout or the reaction count redraws by itself, which is all that button ever did
        # that these two do not, and it invited pressing "FBA" to redraw — quietly replacing
        # an omics or knockout flux state with a fresh FBA solve.
        controls.addWidget(QLabel("draw:"))
        fba_btn = QPushButton("FBA")
        fba_btn.setToolTip("Run FBA and draw it")
        fba_btn.clicked.connect(lambda: self._draw_simulation_on_map(self.run_fba))
        controls.addWidget(fba_btn)
        pfba_btn = QPushButton("pFBA")
        pfba_btn.setToolTip("Run pFBA (minimal total flux) and draw it")
        pfba_btn.clicked.connect(lambda: self._draw_simulation_on_map(self.run_pfba))
        controls.addWidget(pfba_btn)

        controls.addWidget(QLabel("layout:"))
        self.map_layout_combo = QComboBox()
        self.map_layout_combo.addItems([_MAP_LAYOUT_ESCHER, _MAP_LAYOUT_SCHEMATIC])
        self.map_layout_combo.setToolTip(
            "Escher: a curated, publication-style layout covering the reactions it draws.\n"
            "Schematic: no layout file needed — the highest-flux reactions only."
        )
        self.map_layout_combo.currentTextChanged.connect(self._on_map_layout_changed)
        controls.addWidget(self.map_layout_combo)

        self.map_topn_label = QLabel("reactions:")
        controls.addWidget(self.map_topn_label)
        self.map_topn_spin = QSpinBox()
        # Capped low on purpose. The schematic folds the carbon backbone into rows, and a
        # metabolic network branches — past roughly this many reactions the cross-row arrows
        # dominate and the figure stops being readable. That is a curated map's job.
        self.map_topn_spin.setRange(4, 25)
        self.map_topn_spin.setValue(12)
        self.map_topn_spin.setToolTip(
            "How many of the highest-|flux| reactions to draw.\n"
            "Capped at 25 — beyond that the schematic cannot stay readable; "
            "use a curated Escher map instead."
        )
        self.map_topn_spin.valueChanged.connect(self._schedule_map_redraw)
        controls.addWidget(self.map_topn_spin)

        # Named for the file each one takes. "Load map" and "Background" were both true of
        # both buttons — one loads the map's geometry, the other a picture of it — and neither
        # name said which was which.
        load_btn = QPushButton("Map JSON…")
        load_btn.setToolTip(
            "The Escher map file: node coordinates and reaction segments. This is what the "
            "flux is drawn on."
        )
        load_btn.clicked.connect(self.open_flux_map_dialog)
        controls.addWidget(load_btn)

        background_btn = QPushButton("Map image…")
        background_btn.setToolTip(
            "An optional picture of the same map, laid underneath — SVG or PNG exported from "
            "Escher, which draws the network the way Escher does. It changes nothing that is "
            "computed."
        )
        background_btn.clicked.connect(self.open_map_background_dialog)
        controls.addWidget(background_btn)

        self.map_width_check = QCheckBox("width ∝ flux")
        self.map_width_check.setChecked(True)
        self.map_width_check.setToolTip(
            "Scale stroke width with |flux| as well as colour.\n"
            "Turn off for one constant width — colour still carries the magnitude, and a "
            "loaded map image stays visible underneath."
        )
        self.map_width_check.toggled.connect(self._schedule_map_redraw)
        controls.addWidget(self.map_width_check)

        self.map_background_check = QCheckBox("show")
        self.map_background_check.setChecked(True)
        self.map_background_check.setVisible(False)
        self.map_background_check.setToolTip("Hide the background without unloading it")
        self.map_background_check.toggled.connect(self._schedule_map_redraw)
        controls.addWidget(self.map_background_check)
        controls.addStretch(1)
        layout.addLayout(controls)

        self.map_source_label = QLabel()
        self.map_source_label.setWordWrap(True)
        layout.addWidget(self.map_source_label)

        self.map_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.map_canvas_holder)
        layout.addWidget(holder, 1)
        self._map_canvas = None
        self._map_toolbar = None
        # Coalesces a burst of changes — holding the spin arrow steps the value many times a
        # second — into one redraw once the user stops.
        self._map_redraw_timer = QTimer(self)
        self._map_redraw_timer.setSingleShot(True)
        self._map_redraw_timer.setInterval(150)
        self._map_redraw_timer.timeout.connect(self.render_flux_map)
        self._refresh_map_source()
        return tab

    def _schedule_map_redraw(self) -> None:
        """Redraw the loaded flux state after a display change — never solve for one."""

        if self._loading or self._map_canvas is None:
            return  # nothing drawn yet; the draw buttons are the way in
        self._map_redraw_timer.start()

    def _on_map_layout_changed(self, _text: str = "") -> None:
        """Only the schematic has a reaction count to choose, and it says something else."""

        schematic = self.map_layout_combo.currentText() == _MAP_LAYOUT_SCHEMATIC
        self.map_topn_label.setVisible(schematic)
        self.map_topn_spin.setVisible(schematic)
        self._refresh_map_label()
        self._schedule_map_redraw()

    def _refresh_map_label(self) -> None:
        """Describe the layout that is actually on screen, not the one merely available."""

        model_id = html.escape(self.model.id)
        if self.map_layout_combo.currentText() == _MAP_LAYOUT_SCHEMATIC:
            offer = ""
            if self._map_path is not None:
                offer = (
                    " A curated Escher map is available for this model: "
                    f"<b>{html.escape(Path(self._map_path).stem)}</b>."
                )
            self.map_source_label.setText(
                "Schematic layout — the highest-|flux| reactions, currency metabolites "
                f"dropped. No layout file needed.{offer}"
            )
        elif self._map_path is not None:
            drawn = map_reaction_ids(self._map_path)
            present = len(drawn & {rxn.id for rxn in self.model.reactions})
            self.map_source_label.setText(
                f"Escher map: <b>{html.escape(Path(self._map_path).stem)}</b> — "
                f"{present} of its {len(drawn)} reactions are in {model_id}."
            )
        else:  # pragma: no cover - the Escher entry is disabled without a map
            self.map_source_label.setText(f"No Escher map loaded for {model_id}.")

    def _refresh_map_source(self) -> None:
        """Reconcile the layout selector with whether a map suits the current model."""

        combo = self.map_layout_combo
        escher_item = combo.model().item(combo.findText(_MAP_LAYOUT_ESCHER))
        has_map = self._map_path is not None
        escher_item.setEnabled(has_map)
        combo.setCurrentText(_MAP_LAYOUT_ESCHER if has_map else _MAP_LAYOUT_SCHEMATIC)
        self._on_map_layout_changed()
        if not has_map:
            self.map_source_label.setText(
                f"No bundled Escher map fits <b>{html.escape(self.model.id)}</b>, so the "
                "schematic is drawing the highest-|flux| reactions. Load a map for this "
                "model with <i>Map JSON…</i> (from escher.github.io or BiGG)."
            )

    def open_flux_map_dialog(self) -> None:
        """Load a user-supplied Escher map JSON and switch the tab to it."""

        path, _ = QFileDialog.getOpenFileName(
            self, "Open Escher map", "", "Escher maps (*.json);;All files (*)"
        )
        if not path:
            return
        try:
            drawn = map_reaction_ids(path)
        except Exception as exc:
            self.status_label.setText(
                f"Not a readable Escher map: {html.escape(str(exc))}"
            )
            return
        if not drawn:
            self.status_label.setText(
                "That file draws no reactions — not an Escher map."
            )
            return
        present = len(drawn & {rxn.id for rxn in self.model.reactions})
        if present == 0:
            # Rendering would produce an all-grey map; say why rather than draw nothing.
            self.status_label.setText(
                f"None of that map's {len(drawn)} reactions are in {self.model.id} — "
                "it is a map for a different model."
            )
            return
        self._map_path = path
        # The drawing belonged to the previous map, so it no longer describes what is drawn.
        self._clear_map_background()
        self._refresh_map_source()
        self._goto_tab("Flux Map")
        self.render_flux_map()

    def _build_production_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Production design")
        controls_layout = QVBoxLayout(controls)

        # Row 1: the product / substrate / oxygen selectors.
        selectors = QHBoxLayout()
        selectors.addWidget(QLabel("Target product:"))
        self.product_combo = QComboBox()
        self.product_combo.setMinimumWidth(160)
        self.product_combo.setToolTip("Exchange reaction of the target product")
        self.product_combo.currentTextChanged.connect(self._on_product_changed)
        selectors.addWidget(self.product_combo, 1)
        selectors.addWidget(QLabel("substrate:"))
        self.substrate_combo = QComboBox()
        self.substrate_combo.setMinimumWidth(140)
        self.substrate_combo.setToolTip(
            "Carbon substrate (auto = detect from the medium)"
        )
        selectors.addWidget(self.substrate_combo, 1)
        self.anaerobic_combo = QComboBox()
        self.anaerobic_combo.addItems(["aerobic", "anaerobic"])
        selectors.addWidget(self.anaerobic_combo)
        controls_layout.addLayout(selectors)

        # Row 2: the analysis actions, kept on their own line so the controls box does not
        # force the whole window wide.
        actions = QHBoxLayout()
        yield_btn = QPushButton("Theoretical yield")
        yield_btn.clicked.connect(self.run_theoretical_yield)
        envelope_btn = QPushButton("Production envelope")
        envelope_btn.clicked.connect(self.run_production_envelope_plot)
        fseof_btn = QPushButton("FSEOF targets")
        fseof_btn.clicked.connect(self.run_fseof_plot)
        fvseof_btn = QPushButton("FVSEOF (robust)")
        fvseof_btn.clicked.connect(self.run_fvseof_plot)
        actions.addWidget(yield_btn)
        actions.addWidget(envelope_btn)
        actions.addWidget(fseof_btn)
        actions.addWidget(fvseof_btn)
        actions.addStretch(1)
        self._production_buttons = [yield_btn, envelope_btn, fseof_btn, fvseof_btn]
        controls_layout.addLayout(actions)
        layout.addWidget(controls)

        self.yield_label = QLabel("Select a product and compute its theoretical yield.")
        self.yield_label.setFont(QFont("", 13, QFont.Bold))
        self.yield_label.setWordWrap(True)
        self.yield_label.setTextFormat(Qt.RichText)
        layout.addWidget(self.yield_label)

        # Plot on the left, its underlying data table on the right — side by side so the
        # figure and the numbers behind it are visible together (equal 50/50 split).
        results = QSplitter(Qt.Horizontal)
        self.production_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.production_canvas_holder)
        results.addWidget(holder)
        self._production_canvas = None
        self._production_toolbar = None

        # Right side: a header row with the FSEOF/FVSEOF "show all" toggle, above the table.
        right_box = QWidget()
        right_layout = QVBoxLayout(right_box)
        right_layout.setContentsMargins(0, 0, 0, 0)
        table_header = QHBoxLayout()
        table_header.addWidget(QLabel("Result table:"))
        table_header.addStretch(1)
        self.production_show_all_check = QCheckBox("Show all reactions")
        self.production_show_all_check.setToolTip(
            "FSEOF/FVSEOF only. Off: amplify/knockdown targets. On: every reaction "
            "(unchanged ones included)."
        )
        self.production_show_all_check.toggled.connect(self._rerender_production_table)
        table_header.addWidget(self.production_show_all_check)
        right_layout.addLayout(table_header)

        self.fseof_table = QTableWidget(0, 4)
        self._set_production_table_headers(
            ["Reaction", "Flux @ low enforce", "Flux @ high enforce", "Trend"]
        )
        self.fseof_table.verticalHeader().setVisible(False)
        self.fseof_table.setAlternatingRowColors(True)
        right_layout.addWidget(self.fseof_table)
        results.addWidget(right_box)
        results.setStretchFactor(0, 1)
        results.setStretchFactor(1, 1)
        results.setSizes([10000, 10000])
        layout.addWidget(results, 1)
        return tab

    def _build_strain_design_tab(self) -> QWidget:
        """Growth-coupled strain design: OptKnock / RobustKnock over the library service."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Growth-coupled strain design (OptKnock / RobustKnock)")
        form = QVBoxLayout(controls)

        # Row 1: target product + method.
        pick_row = QHBoxLayout()
        pick_row.addWidget(QLabel("Target product:"))
        self.sd_product_combo = QComboBox()
        self.sd_product_combo.setMinimumWidth(150)
        self.sd_product_combo.setToolTip(
            "Exchange reaction of the target product to couple to growth"
        )
        pick_row.addWidget(self.sd_product_combo, 1)
        pick_row.addWidget(QLabel("Method:"))
        self.sd_method_combo = QComboBox()
        self.sd_method_combo.addItems(["optknock", "robustknock"])
        self.sd_method_combo.setToolTip(
            "OptKnock maximizes product at max growth (optimistic); RobustKnock keeps only "
            "designs that guarantee product at max growth (worst case)."
        )
        pick_row.addWidget(self.sd_method_combo)
        # Aeration is a property of the design search, not of whichever tab was last used.
        # optknock/robustknock accepted no condition before 0.4.0, so a run intended to be
        # anaerobic silently returned an aerobic design — the defect that put an aerobic
        # succinate design into the documentation as the anaerobic answer.
        self.sd_anaerobic_combo = QComboBox()
        self.sd_anaerobic_combo.addItems(["aerobic", "anaerobic"])
        self.sd_anaerobic_combo.setToolTip(
            "Aeration the design is searched under, passed as a Condition to "
            "optknock/robustknock and recorded in the result's provenance."
        )
        pick_row.addWidget(self.sd_anaerobic_combo)
        form.addLayout(pick_row)

        # Row 2: search limits + run.
        limit_row = QHBoxLayout()
        limit_row.addWidget(QLabel("max KOs:"))
        self.sd_max_ko_spin = QSpinBox()
        self.sd_max_ko_spin.setRange(1, 6)
        self.sd_max_ko_spin.setValue(3)
        self.sd_max_ko_spin.setToolTip(
            "Maximum number of reaction knockouts per design"
        )
        limit_row.addWidget(self.sd_max_ko_spin)
        limit_row.addWidget(QLabel("solutions:"))
        self.sd_max_sol_spin = QSpinBox()
        self.sd_max_sol_spin.setRange(1, 20)
        self.sd_max_sol_spin.setValue(5)
        self.sd_max_sol_spin.setToolTip("Maximum number of designs to enumerate")
        limit_row.addWidget(self.sd_max_sol_spin)
        self.sd_run_btn = QPushButton("Run design")
        self.sd_run_btn.clicked.connect(self.run_strain_design)
        limit_row.addWidget(self.sd_run_btn)
        limit_row.addStretch(1)
        form.addLayout(limit_row)
        layout.addWidget(controls)

        self.sd_summary = QLabel(
            "OptKnock/RobustKnock search knockout sets that couple product to growth. "
            "Needs a MILP solver (Gurobi/CPLEX)."
        )
        self.sd_summary.setWordWrap(True)
        layout.addWidget(self.sd_summary)

        self.sd_table = QTableWidget(0, 5)
        self.sd_table.setHorizontalHeaderLabels(
            ["Rank", "Knockouts", "Growth", "Max product", "Guaranteed product"]
        )
        _sd_header = self.sd_table.horizontalHeader()
        _sd_header.setSectionResizeMode(1, QHeaderView.Stretch)
        for _col in (0, 2, 3, 4):
            _sd_header.setSectionResizeMode(_col, QHeaderView.ResizeToContents)
        self.sd_table.verticalHeader().setVisible(False)
        self.sd_table.setAlternatingRowColors(True)
        layout.addWidget(self.sd_table, 1)
        return tab

    def run_strain_design(self) -> None:
        product = self.sd_product_combo.currentText()
        if not product:
            self.sd_summary.setText("Select a target product exchange reaction first.")
            return
        method = self.sd_method_combo.currentText()
        max_ko = self.sd_max_ko_spin.value()
        max_sol = self.sd_max_sol_spin.value()
        condition = self._strain_design_condition()
        aeration = self.sd_anaerobic_combo.currentText()
        self.sd_summary.setText(
            f"Searching {method} designs for {html.escape(product)} ({aeration})…"
        )
        self.status_label.setText(f"Running {method} (this can take a while)…")
        self.sd_run_btn.setEnabled(False)
        try:
            search = optknock if method == "optknock" else robustknock
            result = self._run_model_in_background(
                lambda worker_model: search(
                    worker_model,
                    product,
                    max_knockouts=max_ko,
                    max_solutions=max_sol,
                    condition=condition,
                ),
                label=f"Running {method}…",
            )
        except (
            Exception
        ) as exc:  # surface solver-capability/model errors instead of freezing
            self.sd_table.setRowCount(0)
            self.sd_summary.setText(f"Strain design failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Strain design failed ({method}).")
            return
        finally:
            self.sd_run_btn.setEnabled(True)

        designs = result.designs
        self.sd_table.setRowCount(len(designs))
        for i, design in enumerate(designs):
            self.sd_table.setItem(i, 0, QTableWidgetItem(str(i + 1)))
            self.sd_table.setItem(
                i, 1, QTableWidgetItem(", ".join(design.knockouts) or "—")
            )
            self.sd_table.setItem(i, 2, QTableWidgetItem(f"{design.growth:.4g}"))
            self.sd_table.setItem(i, 3, QTableWidgetItem(f"{design.max_product:.4g}"))
            self.sd_table.setItem(
                i, 4, QTableWidgetItem(f"{design.guaranteed_product:.4g}")
            )
            if design.growth_coupled:
                self.sd_table.item(i, 4).setBackground(QColor("#dce9d6"))
        if designs:
            best = result.best()
            self.sd_summary.setText(
                f"<b>{method}</b> for {html.escape(product)} ({aeration}): "
                f"{len(designs)} design(s). "
                f"Best knocks out <b>{html.escape(', '.join(best.knockouts) or '—')}</b> "
                f"(growth {best.growth:.3g} h⁻¹, guaranteed product {best.guaranteed_product:.3g}, "
                f"max product {best.max_product:.3g}). Candidates are restricted to "
                "gene-associated internal reactions; exchange knockouts are not realisable "
                "as gene deletions."
            )
        else:
            self.sd_summary.setText(
                f"<b>{method}</b> found no growth-coupled design for "
                f"{html.escape(product)} ({aeration}) — try a different product, more "
                "knockouts, or check the product can carry flux in the current medium."
            )
        self.status_label.setText(f"Strain design complete ({method}, {aeration}).")

    def _build_flux_response_tab(self) -> QWidget:
        """Scan one reaction's enforced flux and watch a response reaction follow it."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox(
            "Flux response (enforced target flux → maximized response)"
        )
        controls_layout = QVBoxLayout(controls)

        selectors = QHBoxLayout()
        selectors.addWidget(QLabel("Target (scanned):"))
        self.fr_target_combo = QComboBox()
        self.fr_target_combo.setEditable(True)
        self.fr_target_combo.setInsertPolicy(QComboBox.NoInsert)
        self.fr_target_combo.setMinimumWidth(160)
        self.fr_target_combo.setToolTip(
            "Reaction whose flux is fixed at each step of the scan"
        )
        self.fr_target_combo.currentTextChanged.connect(self._on_fr_target_changed)
        selectors.addWidget(self.fr_target_combo, 1)
        selectors.addWidget(QLabel("Response (maximized):"))
        self.fr_response_combo = QComboBox()
        self.fr_response_combo.setEditable(True)
        self.fr_response_combo.setInsertPolicy(QComboBox.NoInsert)
        self.fr_response_combo.setMinimumWidth(160)
        self.fr_response_combo.setToolTip(
            "Reaction maximized at each step. '(objective)' gives the robustness reading; "
            "a product exchange gives the production reading."
        )
        selectors.addWidget(self.fr_response_combo, 1)
        controls_layout.addLayout(selectors)

        options = QHBoxLayout()
        options.addWidget(QLabel("Steps:"))
        self.fr_steps_spin = QSpinBox()
        self.fr_steps_spin.setRange(2, 200)
        self.fr_steps_spin.setValue(20)
        options.addWidget(self.fr_steps_spin)

        options.addWidget(QLabel("Min growth (% of wild type):"))
        self.fr_growth_spin = QDoubleSpinBox()
        self.fr_growth_spin.setRange(0.0, 100.0)
        self.fr_growth_spin.setValue(0.0)
        self.fr_growth_spin.setSingleStep(5.0)
        self.fr_growth_spin.setToolTip(
            "Hold biomass at this fraction of the wild-type optimum across the scan. "
            "0 = unconstrained, which for a product response returns the theoretical "
            "ceiling of a non-growing cell."
        )
        options.addWidget(self.fr_growth_spin)

        # The range is always shown and always editable. There is no auto/manual toggle: the
        # boxes are filled with the target's detected feasible interval whenever the target
        # changes, so "auto" and "manual" would hold the same numbers until you edit them.
        options.addWidget(QLabel("Range min:"))
        self.fr_min_spin = QDoubleSpinBox()
        self.fr_min_spin.setRange(-1e6, 1e6)
        self.fr_min_spin.setDecimals(3)
        self.fr_min_spin.setToolTip(
            "Lowest enforced target flux. Filled with the target's feasible minimum; edit to "
            "scan a narrower window, or past its bounds as a what-if."
        )
        options.addWidget(self.fr_min_spin)
        options.addWidget(QLabel("max:"))
        self.fr_max_spin = QDoubleSpinBox()
        self.fr_max_spin.setRange(-1e6, 1e6)
        self.fr_max_spin.setDecimals(3)
        options.addWidget(self.fr_max_spin)

        self.fr_detect_btn = QPushButton("Detect range")
        self.fr_detect_btn.setToolTip(
            "Reset min/max to the full flux interval this reaction can carry in the current "
            "medium (FVA with the objective unconstrained)"
        )
        self.fr_detect_btn.clicked.connect(self._prefill_fr_range)
        options.addWidget(self.fr_detect_btn)

        self.fr_run_btn = QPushButton("Run scan")
        self.fr_run_btn.clicked.connect(self.run_flux_response)
        options.addWidget(self.fr_run_btn)
        options.addStretch(1)
        controls_layout.addLayout(options)
        layout.addWidget(controls)

        self.fr_summary = QLabel(
            "Pick a target reaction and run the scan to see how the response follows it."
        )
        self.fr_summary.setFont(QFont("", 13, QFont.Bold))
        self.fr_summary.setWordWrap(True)
        self.fr_summary.setTextFormat(Qt.RichText)
        layout.addWidget(self.fr_summary)

        results = QSplitter(Qt.Horizontal)
        self.fr_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.fr_canvas_holder)
        results.addWidget(holder)
        self._flux_response_canvas = None
        self._flux_response_toolbar = None

        right_box = QWidget()
        right_layout = QVBoxLayout(right_box)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(QLabel("Scan points:"))
        self.fr_table = QTableWidget(0, 4)
        self.fr_table.setHorizontalHeaderLabels(
            ["Target flux", "Response flux", "Biomass flux", "Status"]
        )
        self.fr_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.fr_table.verticalHeader().setVisible(False)
        self.fr_table.setAlternatingRowColors(True)
        right_layout.addWidget(self.fr_table)

        # The phase table is the replacement for the old "bottleneck" line: the response
        # curve is concave piecewise linear in the enforced flux, so the exact answer is the
        # list of intervals of constant shadow price and the boundary where it turns
        # negative — not an argmin over finite differences of whatever grid was scanned.
        phase_label = QLabel("Response phases (intervals of constant shadow price):")
        phase_label.setToolTip(
            "Each phase is one linear piece of the response curve. The shadow price is the "
            "exact d(response)/d(target) on that piece, read from the LP dual, so it does "
            "not move with the number of scan steps."
        )
        right_layout.addWidget(phase_label)
        self.fr_phase_table = QTableWidget(0, 4)
        self.fr_phase_table.setHorizontalHeaderLabels(
            ["Target from", "to", "Response at end", "Shadow price"]
        )
        self.fr_phase_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.fr_phase_table.verticalHeader().setVisible(False)
        self.fr_phase_table.setAlternatingRowColors(True)
        right_layout.addWidget(self.fr_phase_table)
        results.addWidget(right_box)
        results.setStretchFactor(0, 1)
        results.setStretchFactor(1, 1)
        results.setSizes([10000, 10000])
        layout.addWidget(results, 1)
        return tab

    def _build_sampling_tab(self) -> QWidget:
        """Random flux sampling: what the feasible space looks like, not one optimum."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Random flux sampling")
        controls_layout = QVBoxLayout(controls)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Mode:"))
        self.sample_mode_combo = QComboBox()
        self.sample_mode_combo.addItems(["uniform", "around a reference"])
        self.sample_mode_combo.setToolTip(
            "'uniform' samples the whole feasible space; 'around a reference' narrows every "
            "reaction to a window around a predicted flux state first."
        )
        self.sample_mode_combo.currentTextChanged.connect(self._on_sample_mode_changed)
        row1.addWidget(self.sample_mode_combo)

        # These two only mean anything in "around a reference" mode, so they are disabled —
        # labels included — in uniform mode rather than sitting there inviting input.
        self.sample_reference_label = QLabel("Reference:")
        row1.addWidget(self.sample_reference_label)
        self.sample_reference_combo = QComboBox()
        self.sample_reference_combo.addItems(["pfba", "fba"])
        self.sample_reference_combo.setToolTip(
            "Which predicted flux state the sampling windows are centred on"
        )
        row1.addWidget(self.sample_reference_combo)

        self.sample_window_label = QLabel("window ±%:")
        row1.addWidget(self.sample_window_label)
        self.sample_window_spin = QDoubleSpinBox()
        self.sample_window_spin.setRange(1.0, 100.0)
        self.sample_window_spin.setValue(20.0)
        self.sample_window_spin.setToolTip(
            "Half-width of the window around each reference flux, as a percentage"
        )
        row1.addWidget(self.sample_window_spin)
        row1.addStretch(1)
        controls_layout.addLayout(row1)
        self._on_sample_mode_changed(self.sample_mode_combo.currentText())

        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Samples:"))
        self.sample_n_spin = QSpinBox()
        self.sample_n_spin.setRange(10, 100000)
        self.sample_n_spin.setValue(1000)
        row2.addWidget(self.sample_n_spin)

        row2.addWidget(QLabel("Sampler:"))
        self.sample_method_combo = QComboBox()
        self.sample_method_combo.addItems(["optgp", "achr"])
        self.sample_method_combo.setToolTip(
            "optgp needs a large sample count to mix; achr converges better for small runs"
        )
        row2.addWidget(self.sample_method_combo)

        row2.addWidget(QLabel("Thinning:"))
        self.sample_thinning_spin = QSpinBox()
        self.sample_thinning_spin.setRange(1, 10000)
        self.sample_thinning_spin.setValue(100)
        self.sample_thinning_spin.setToolTip(
            "Keep every nth iterate to reduce autocorrelation between samples"
        )
        row2.addWidget(self.sample_thinning_spin)

        row2.addWidget(QLabel("Seed:"))
        self.sample_seed_spin = QSpinBox()
        self.sample_seed_spin.setRange(0, 2**31 - 1)
        self.sample_seed_spin.setValue(0)
        self.sample_seed_spin.setToolTip(
            "The same seed reproduces the same samples exactly"
        )
        row2.addWidget(self.sample_seed_spin)

        row2.addWidget(QLabel("Show top:"))
        self.sample_top_spin = QSpinBox()
        self.sample_top_spin.setRange(1, 40)
        self.sample_top_spin.setValue(8)
        self.sample_top_spin.setToolTip(
            "How many of the most variable reactions to draw in the figure"
        )
        row2.addWidget(self.sample_top_spin)

        self.sample_run_btn = QPushButton("Run sampling")
        self.sample_run_btn.clicked.connect(self.run_sampling)
        row2.addWidget(self.sample_run_btn)

        # The table summarizes; the ensemble itself is the primary data, so it gets its own
        # export rather than being reachable only through File > Export Table.
        self.sample_export_btn = QPushButton("Export samples…")
        self.sample_export_btn.setToolTip(
            "Write every drawn sample to CSV: one row per sample, one column per reaction"
        )
        self.sample_export_btn.setEnabled(False)
        self.sample_export_btn.clicked.connect(self.export_sampling_samples)
        row2.addWidget(self.sample_export_btn)
        row2.addStretch(1)
        controls_layout.addLayout(row2)

        # Optional knockouts: sampling a deletion strain answers "what does this mutant's
        # flux space look like", which a wild-type ensemble cannot. Empty = wild type.
        row3 = QHBoxLayout()
        row3.addWidget(QLabel("Knockouts (optional) — level:"))
        self.sample_ko_level_combo = QComboBox()
        self.sample_ko_level_combo.addItems(["reaction", "gene"])
        self.sample_ko_level_combo.setToolTip(
            "Knock out reactions directly, or genes resolved to reactions through the GPR"
        )
        self.sample_ko_level_combo.currentTextChanged.connect(
            lambda level: self._ko_populate("sample_ko", level)
        )
        row3.addWidget(self.sample_ko_level_combo)
        row3.addStretch(1)
        controls_layout.addLayout(row3)
        controls_layout.addLayout(self._build_ko_picker("sample_ko", max_height=110))

        layout.addWidget(controls)
        self._sampling_result = None

        self.sample_summary = QLabel(
            "Sampling shows how much a flux can vary across equally valid solutions — "
            "run it to see which predictions are forced and which are one of many."
        )
        self.sample_summary.setFont(QFont("", 13, QFont.Bold))
        self.sample_summary.setWordWrap(True)
        self.sample_summary.setTextFormat(Qt.RichText)
        layout.addWidget(self.sample_summary)

        results = QSplitter(Qt.Horizontal)
        self.sample_canvas_holder = QVBoxLayout()
        holder = QWidget()
        holder.setLayout(self.sample_canvas_holder)
        results.addWidget(holder)
        self._sampling_canvas = None
        self._sampling_toolbar = None

        right_box = QWidget()
        right_layout = QVBoxLayout(right_box)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(QLabel("Per-reaction statistics:"))
        self.sample_table = QTableWidget(0, 6)
        self.sample_table.setHorizontalHeaderLabels(
            ["Reaction", "Mean", "Std", "Min", "Median", "Max"]
        )
        self.sample_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.sample_table.verticalHeader().setVisible(False)
        self.sample_table.setAlternatingRowColors(True)
        right_layout.addWidget(self.sample_table)
        results.addWidget(right_box)
        results.setStretchFactor(0, 1)
        results.setStretchFactor(1, 1)
        results.setSizes([10000, 10000])
        layout.addWidget(results, 1)
        return tab

    def _on_fr_target_changed(self, target: str) -> None:
        # A range belongs to the reaction it was detected for; re-seed on a switch so a scan
        # never runs another reaction's numbers.
        if self._loading:
            return
        self._prefill_fr_range()

    def _prefill_fr_range(self) -> None:
        """Fill the range boxes with the target's feasible interval in the current medium.

        Runs on target change and on demand from "Detect range", so the boxes never hold a
        placeholder or a previous reaction's numbers.
        """

        target = self.fr_target_combo.currentText().strip()
        if not target or target not in self.model.reactions:
            return
        try:
            detected = self._run_model_in_background(
                lambda worker_model: fva(
                    worker_model, reactions=[target], fraction_of_optimum=0.0
                )[target],
                label="Detecting the target's flux range…",
            )
        except Exception as exc:
            self.status_label.setText(f"Could not detect the flux range: {exc}")
            return
        self.fr_min_spin.setValue(float(detected.minimum))
        self.fr_max_spin.setValue(float(detected.maximum))
        self.status_label.setText(
            f"{target} can carry {detected.minimum:.4g} to {detected.maximum:.4g} "
            "mmol gDW⁻¹ h⁻¹ in this medium."
        )

    def _on_sample_mode_changed(self, mode: str) -> None:
        around_reference = mode == "around a reference"
        for widget in (
            self.sample_reference_label,
            self.sample_reference_combo,
            self.sample_window_label,
            self.sample_window_spin,
        ):
            widget.setEnabled(around_reference)

    def run_flux_response(self) -> None:
        """Run the flux-response scan and render its curve, table, and summary."""

        try:
            target = self.fr_target_combo.currentText().strip()
            if not target:
                self.fr_summary.setText("Select a target reaction first.")
                return
            response_text = self.fr_response_combo.currentText().strip()
            response = None if response_text in ("", "(objective)") else response_text
            growth_percent = self.fr_growth_spin.value()
            biomass_fraction = None if growth_percent <= 0 else growth_percent / 100.0
            steps = self.fr_steps_spin.value()
            target_min = self.fr_min_spin.value()
            target_max = self.fr_max_spin.value()

            result = self._run_model_in_background(
                lambda worker_model: flux_response(
                    worker_model,
                    target,
                    response,
                    n_steps=steps,
                    biomass_fraction=biomass_fraction,
                    target_min=target_min,
                    target_max=target_max,
                ),
                label="Running flux response scan…",
            )
        except Exception as exc:
            self.fr_summary.setText(f"Flux response failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Flux response failed: {exc}")
            return

        self._set_figure(
            self.fr_canvas_holder, "flux_response", flux_response_figure(result)
        )
        self._fill_flux_response_table(result)

        optimum = result.optimum()
        window = result.feasible_range()
        parts = [
            f"<b>{html.escape(result.response)}</b> vs enforced "
            f"<b>{html.escape(result.target)}</b>."
        ]
        if optimum is not None:
            parts.append(
                f"Best response {optimum.response_flux:.4g} at target "
                f"{optimum.target_flux:.4g} (growth {optimum.biomass_flux:.4g})."
            )
        if window is not None:
            parts.append(
                f"Feasible from {window[0]:.4g} to {window[1]:.4g} (FVA on the scanned "
                "reaction, not the grid)."
            )
        # The response curve is the optimal-value function of an LP in one right-hand side,
        # so it is concave piecewise linear: each phase is an interval of constant shadow
        # price, and the limit is the exact boundary where the price crosses the threshold.
        limit = result.limit
        if limit.found and limit.target_flux is not None:
            parts.append(
                f"Response falls beyond <b>{limit.target_flux:.4g}</b> — shadow price "
                f"{limit.shadow_price_before:.4g} → {limit.shadow_price_after:.4g} per unit "
                f"of {html.escape(result.target)}."
            )
        else:
            parts.append(f"No response limit: {html.escape(limit.message)}.")
        if result.phases:
            parts.append(
                f"{len(result.phases)} phase(s) of constant shadow price; at the wild-type "
                f"flux the price is {self._wild_type_shadow_price(result)}."
            )
        self.fr_summary.setText(" ".join(parts))
        self._fill_flux_response_phase_table(result)
        self.status_label.setText("Flux response scan complete.")

    @staticmethod
    def _wild_type_shadow_price(result) -> str:
        """The exact d(response)/d(target) at the wild-type target flux, or why not."""

        wild_type = result.wild_type.get("target_flux")
        if wild_type is None:
            return "—"
        try:
            return f"{result.shadow_price_at(float(wild_type)):.4g}"
        except ValueError:
            # The wild-type flux can sit outside the scanned window; say so rather than
            # printing a number from a phase that does not contain it.
            return "— (wild type is outside the scanned window)"

    def _fill_flux_response_phase_table(self, result) -> None:
        """List the response phases: one row per interval of constant shadow price."""

        frame = result.phases_frame()
        self.fr_phase_table.setRowCount(len(frame))
        for row, (_, record) in enumerate(frame.iterrows()):
            values = [
                f"{record['target_low']:.4g}",
                f"{record['target_high']:.4g}",
                f"{record['response_high']:.4g}",
                f"{record['shadow_price']:.4g}",
            ]
            for col, text in enumerate(values):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                # A negative shadow price means the response is paying for the enforced
                # flux; flag it the same way an infeasible scan point is flagged.
                if col == 3 and record["shadow_price"] < 0:
                    item.setForeground(QColor("#B00020"))
                self.fr_phase_table.setItem(row, col, item)

    def _fill_flux_response_table(self, result) -> None:
        frame = result.to_frame()
        self.fr_table.setRowCount(len(frame))
        for row, (_, record) in enumerate(frame.iterrows()):
            values = [
                f"{record['target_flux']:.4g}",
                "—"
                if math.isnan(record["response_flux"])
                else f"{record['response_flux']:.4g}",
                "—"
                if math.isnan(record["biomass_flux"])
                else f"{record['biomass_flux']:.4g}",
                str(record["status"]),
            ]
            for col, text in enumerate(values):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                if record["status"] != "optimal":
                    item.setForeground(QColor("#B00020"))
                self.fr_table.setItem(row, col, item)

    def _sampling_knockout_condition(self) -> Condition | None:
        """The Sampling tab's chosen knockouts as a Condition, or None for wild type."""

        level = self.sample_ko_level_combo.currentText()
        blocked = self._ko_blocked_reactions("sample_ko", level)
        if not blocked:
            return None
        return Condition(
            name=f"knockout:{','.join(self._ko_targets('sample_ko'))}",
            bounds=tuple(
                ReactionBound(reaction_id=rid, lower_bound=0.0, upper_bound=0.0)
                for rid in blocked
            ),
        )

    @staticmethod
    def _reference_under(model: Model, condition: Condition | None, method: str):
        """Build a reference flux state with the condition applied, then restore the model."""

        with model:
            if condition is not None:
                condition.apply_to(model)
            return reference_flux(model, method)

    def run_sampling(self) -> None:
        """Draw a flux ensemble and summarize how much each reaction actually varies."""

        try:
            n = self.sample_n_spin.value()
            method = self.sample_method_combo.currentText()
            thinning = self.sample_thinning_spin.value()
            seed = self.sample_seed_spin.value()
            around_reference = (
                self.sample_mode_combo.currentText() == "around a reference"
            )
            condition = self._sampling_knockout_condition()

            if around_reference:
                reference_method = self.sample_reference_combo.currentText()
                half_width = self.sample_window_spin.value() / 100.0
                result = self._run_model_in_background(
                    lambda worker_model: reference_constrained_sampling(
                        worker_model,
                        # The reference must come from the knocked-out model: a wild-type
                        # reference would put deleted reactions outside their own windows.
                        self._reference_under(
                            worker_model, condition, reference_method
                        ),
                        n=n,
                        condition=condition,
                        method=method,
                        thinning=thinning,
                        seed=seed,
                        min_fraction=max(0.0, 1.0 - half_width),
                        max_fraction=1.0 + half_width,
                    ),
                    label="Sampling around the reference…",
                )
            else:
                result = self._run_model_in_background(
                    lambda worker_model: random_flux_sampling(
                        worker_model,
                        n=n,
                        condition=condition,
                        method=method,
                        thinning=thinning,
                        seed=seed,
                    ),
                    label="Sampling the feasible space…",
                )
        except Exception as exc:
            # Drop the previous ensemble: leaving it on screen (and exportable) next to a
            # failure message would attribute one run's samples to another run's settings.
            self._sampling_result = None
            self.sample_export_btn.setEnabled(False)
            self.sample_table.setRowCount(0)
            if condition is not None:
                # A deletion that leaves no feasible space is a lethal knockout, which is a
                # result about the strain rather than a failure of the run.
                self.sample_summary.setText(
                    f"No flux distribution could be sampled with these knockouts applied "
                    f"({html.escape(', '.join(self._ko_targets('sample_ko')))}). The set is "
                    f"probably lethal, or leaves a feasible space of zero volume. "
                    f"<br><small>{html.escape(str(exc))}</small>"
                )
                self.status_label.setText("Sampling: knockout set appears lethal.")
                return
            self.sample_summary.setText(f"Sampling failed: {html.escape(str(exc))}")
            self.status_label.setText(f"Sampling failed: {exc}")
            return

        top_n = self.sample_top_spin.value()
        self._set_figure(
            self.sample_canvas_holder,
            "sampling",
            sampling_figure(result, top_n=top_n),
        )
        self._fill_sampling_table(result)
        self._sampling_result = result
        self.sample_export_btn.setEnabled(True)

        statistics = result.statistics()
        varying = int((statistics["std"] > 1e-6).sum())
        if condition is None:
            strain = "wild type"
        else:
            targets = self._ko_targets("sample_ko")
            level = self.sample_ko_level_combo.currentText()
            strain = (
                f"<b>{len(targets)}</b> {level} knockout(s): "
                f"{html.escape(', '.join(targets))}"
            )
        self.sample_summary.setText(
            f"Drew <b>{result.n_samples}</b> samples ({html.escape(result.method)}, "
            f"seed {result.seed}, thinning {thinning}) on {strain}. "
            f"<b>{varying}</b> of {len(statistics)} reactions vary across the ensemble; "
            "the rest are fixed by the constraints."
        )
        self.status_label.setText("Flux sampling complete.")

    def export_sampling_samples(self) -> None:
        """Write the full sampled ensemble to CSV — one row per sample, one column per reaction."""

        result = self._sampling_result
        if result is None:
            self.status_label.setText("Run sampling first, then export the ensemble.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export samples to CSV", "", "CSV files (*.csv);;All files (*)"
        )
        if not path:
            return
        frame = result.to_frame()
        try:
            frame.to_csv(path, index_label="sample")
        except OSError as exc:
            self.status_label.setText(f"Could not write CSV: {exc}")
            return
        self.status_label.setText(
            f"Exported {len(frame)} samples x {len(frame.columns)} reactions to {path}."
        )

    def _fill_sampling_table(self, result) -> None:
        statistics = result.statistics().sort_values("std", ascending=False)
        self.sample_table.setRowCount(len(statistics))
        for row, (reaction_id, record) in enumerate(statistics.iterrows()):
            values = [
                str(reaction_id),
                f"{record['mean']:.4g}",
                f"{record['std']:.4g}",
                f"{record['minimum']:.4g}",
                f"{record['median']:.4g}",
                f"{record['maximum']:.4g}",
            ]
            for col, text in enumerate(values):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.sample_table.setItem(row, col, item)

    def _build_comparison_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        controls = QGroupBox(
            "Perturbation response (MOMA / ROOM vs a reference template)"
        )
        controls_layout = QVBoxLayout(controls)

        row = QHBoxLayout()
        row.addWidget(QLabel("Method:"))
        self.comparison_method_combo = QComboBox()
        self.comparison_method_combo.addItems(["MOMA (L2)", "MOMA (L1)", "ROOM"])
        self.comparison_method_combo.currentTextChanged.connect(
            self._on_comparison_method_changed
        )
        row.addWidget(self.comparison_method_combo)
        # Shlomi et al. publish two tolerance pairs for two different questions, and the
        # choice is not cosmetic: 531 vs 401 switches over the same 35 knockouts.
        self.room_use_case_combo = QComboBox()
        self.room_use_case_combo.addItems(sorted(ROOM_TOLERANCES))
        self.room_use_case_combo.setCurrentText("flux_prediction")
        self.room_use_case_combo.setToolTip(
            "ROOM tolerance pair (Shlomi et al. 2005): flux_prediction δ=0.03 ε=1e-3, "
            "lethality δ=0.1 ε=1e-2. Ignored by MOMA."
        )
        row.addWidget(self.room_use_case_combo)
        row.addWidget(QLabel("Reference template:"))
        self.template_combo = QComboBox()
        # pFBA first, and therefore the default: `reference_flux` defaults to "pfba" too, so
        # the GUI and the Python API no longer disagree about the reference state. FBA is
        # Segre et al.'s own protocol and stays available; it is non-unique (8 identical
        # calls on iJO1366 returned 3 distinct vectors), which is why pFBA is the default.
        self.template_combo.addItems(["pfba", "fba", "lad", "eflux2"])
        self.template_combo.setToolTip(
            "Wild-type flux template the perturbed state is compared to. pFBA is the default "
            "in both the GUI and the Python API; FBA is Segrè et al.'s protocol but has "
            "alternate optima. LAD/E-Flux2 reuse the expression computed on the Omics tab, "
            "and stay disabled until that is available."
        )
        # LAD/E-Flux2 need an omics result; grey them out (delegate + not-allowed cursor) until
        # one exists so they cannot be picked.
        self._style_disabled_items(self.template_combo)
        self._refresh_comparison_templates()
        row.addWidget(self.template_combo)
        row.addWidget(QLabel("Knockout level:"))
        self.ko_level_combo = QComboBox()
        self.ko_level_combo.addItems(["reaction", "gene"])
        self.ko_level_combo.setToolTip(
            "Knock out reactions, or genes (resolved to reactions via GPR)"
        )
        self.ko_level_combo.currentTextChanged.connect(self._populate_ko_list)
        row.addWidget(self.ko_level_combo)
        row.addStretch(1)
        controls_layout.addLayout(row)

        # Row 2: the optional batch product and the single-run significance threshold.
        opts = QHBoxLayout()
        opts.addWidget(QLabel("Batch product (optional):"))
        self.batch_product_combo = QComboBox()
        self.batch_product_combo.setMinimumWidth(140)
        self.batch_product_combo.setToolTip(
            "Optional. When set, the Batch table adds the product's wild-type and post-knockout "
            "flux columns."
        )
        opts.addWidget(self.batch_product_combo, 1)
        # Significant-change threshold for the single-run table: a reaction is shown only when
        # |perturbed - reference| exceeds this fraction of |reference| (plus a small floor).
        # This matches ROOM's own significance criterion and hides alternate-optimum drift.
        opts.addWidget(QLabel("Significant change ≥"))
        self.comparison_threshold_spin = QDoubleSpinBox()
        self.comparison_threshold_spin.setRange(0.0, 100.0)
        self.comparison_threshold_spin.setSingleStep(1.0)
        self.comparison_threshold_spin.setValue(3.0)
        self.comparison_threshold_spin.setSuffix(" % of ref")
        self.comparison_threshold_spin.setToolTip(
            "Show a reaction only when |perturbed − reference| > this % of |reference| + 0.001. "
            "0% shows every reaction that moved at all."
        )
        self.comparison_threshold_spin.valueChanged.connect(
            self._render_comparison_table
        )
        opts.addWidget(self.comparison_threshold_spin)
        controls_layout.addLayout(opts)

        # Row 3: the knock-out actions.
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Knock out:"))
        run_btn = QPushButton("Run (selected as one KO)")
        run_btn.setToolTip("MOMA/ROOM with all selected targets knocked out together")
        run_btn.clicked.connect(self.run_comparison)
        batch_btn = QPushButton("Batch (each separately)")
        batch_btn.setToolTip(
            "Run MOMA/ROOM once per target as a separate single knockout — the selected "
            "targets, or all targets of this level when none are selected."
        )
        batch_btn.clicked.connect(self.run_batch_comparison)
        self.cmp_map_btn = QPushButton("Show on flux map")
        self.cmp_map_btn.setEnabled(False)
        self.cmp_map_btn.setToolTip(
            "Draw the perturbed flux distribution from the last single run on the map."
        )
        self.cmp_map_btn.clicked.connect(self.show_comparison_on_flux_map)
        row2.addWidget(run_btn)
        row2.addWidget(batch_btn)
        row2.addWidget(self.cmp_map_btn)
        row2.addStretch(1)
        controls_layout.addLayout(row2)

        controls_layout.addLayout(self._build_ko_picker("ko"))
        layout.addWidget(controls)

        self.comparison_summary = QLabel(
            "Pick a reference template (FBA/pFBA/LAD/E-Flux2), a knockout level (reaction/gene), "
            "one or more targets, and a method to predict the perturbed flux state."
        )
        self.comparison_summary.setWordWrap(True)
        layout.addWidget(self.comparison_summary)

        self.comparison_table = QTableWidget(0, 3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )
        self.comparison_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        self.comparison_table.verticalHeader().setVisible(False)
        self.comparison_table.setAlternatingRowColors(True)
        layout.addWidget(self.comparison_table, 1)
        return tab

    # -- Knockout picker (shared by the Comparison and Sampling tabs) --------------------
    #
    # Two panels: left = searchable catalogue of all targets, right = the chosen knockout
    # set. Add/remove moves ids between them, so the selection is always plainly visible.
    # Every widget is stored under ``<prefix>_*`` so one implementation serves both tabs.

    def _build_ko_picker(self, prefix: str, *, max_height: int = 150) -> QHBoxLayout:
        picker = QHBoxLayout()

        left = QVBoxLayout()
        search = QLineEdit()
        search.setPlaceholderText("Search targets…")
        search.setClearButtonEnabled(True)
        search.textChanged.connect(lambda text, p=prefix: self._ko_filter(p, text))
        left.addWidget(search)
        available = QListWidget()
        available.setSelectionMode(QAbstractItemView.ExtendedSelection)
        available.setMaximumHeight(max_height)
        available.setToolTip("Double-click (or select and Add →) to add a target.")
        available.itemDoubleClicked.connect(
            lambda item, p=prefix: self._ko_add(p, [item.text()])
        )
        left.addWidget(available)
        picker.addLayout(left, 1)

        mid = QVBoxLayout()
        mid.addStretch(1)
        add_btn = QPushButton("Add →")
        add_btn.setToolTip("Add the highlighted targets to the knockout set")
        add_btn.clicked.connect(
            lambda _=False, p=prefix: self._ko_add(
                p, [i.text() for i in getattr(self, f"{p}_available").selectedItems()]
            )
        )
        remove_btn = QPushButton("← Remove")
        remove_btn.setToolTip("Remove the highlighted targets from the knockout set")
        remove_btn.clicked.connect(
            lambda _=False, p=prefix: self._ko_remove(
                p, [i.text() for i in getattr(self, f"{p}_selected").selectedItems()]
            )
        )
        clear_btn = QPushButton("Clear")
        clear_btn.setToolTip("Empty the knockout set")
        clear_btn.clicked.connect(lambda _=False, p=prefix: self._ko_clear(p))
        mid.addWidget(add_btn)
        mid.addWidget(remove_btn)
        mid.addWidget(clear_btn)
        mid.addStretch(1)
        picker.addLayout(mid)

        right = QVBoxLayout()
        selected_label = QLabel("Selected (0):")
        right.addWidget(selected_label)
        selected = QListWidget()
        selected.setSelectionMode(QAbstractItemView.ExtendedSelection)
        selected.setMaximumHeight(max_height)
        selected.setToolTip(
            "The knockout set. Double-click (or select and ← Remove) to drop one."
        )
        selected.itemDoubleClicked.connect(
            lambda item, p=prefix: self._ko_remove(p, [item.text()])
        )
        right.addWidget(selected)
        picker.addLayout(right, 1)

        setattr(self, f"{prefix}_filter", search)
        setattr(self, f"{prefix}_available", available)
        setattr(self, f"{prefix}_selected", selected)
        setattr(self, f"{prefix}_selected_label", selected_label)
        return picker

    def _ko_populate(self, prefix: str, level: str) -> None:
        """Fill the catalogue for the given level; clear the chosen set."""

        if not hasattr(self, f"{prefix}_available"):
            return
        ids = (
            [g.id for g in self.model.genes]
            if level == "gene"
            else [r.id for r in self.model.reactions]
        )
        available = getattr(self, f"{prefix}_available")
        available.clear()
        available.addItems(ids)
        # Switching level invalidates the previous selection: gene ids are not reaction ids.
        getattr(self, f"{prefix}_selected").clear()
        self._ko_filter(prefix, getattr(self, f"{prefix}_filter").text())
        self._ko_update_label(prefix)

    def _ko_filter(self, prefix: str, text: str) -> None:
        """Hide catalogue targets that do not match the search box (case-insensitive)."""

        query = text.strip().lower()
        available = getattr(self, f"{prefix}_available")
        for i in range(available.count()):
            item = available.item(i)
            item.setHidden(bool(query) and query not in item.text().lower())

    def _ko_add(self, prefix: str, ids) -> None:
        """Add target ids to the chosen set, skipping ones already present."""

        selected = getattr(self, f"{prefix}_selected")
        existing = {selected.item(i).text() for i in range(selected.count())}
        for target_id in ids:
            if target_id not in existing:
                selected.addItem(target_id)
                existing.add(target_id)
        self._ko_update_label(prefix)

    def _ko_remove(self, prefix: str, ids) -> None:
        selected = getattr(self, f"{prefix}_selected")
        drop = set(ids)
        for i in range(selected.count() - 1, -1, -1):
            if selected.item(i).text() in drop:
                selected.takeItem(i)
        self._ko_update_label(prefix)

    def _ko_clear(self, prefix: str) -> None:
        getattr(self, f"{prefix}_selected").clear()
        self._ko_update_label(prefix)

    def _ko_update_label(self, prefix: str) -> None:
        label = getattr(self, f"{prefix}_selected_label", None)
        if label is not None:
            count = getattr(self, f"{prefix}_selected").count()
            label.setText(f"Selected ({count}):")

    def _ko_targets(self, prefix: str) -> list[str]:
        selected = getattr(self, f"{prefix}_selected")
        return [selected.item(i).text() for i in range(selected.count())]

    def _ko_blocked_reactions(self, prefix: str, level: str) -> list[str]:
        """Chosen targets resolved to reaction ids, through the GPR when the level is gene."""

        targets = self._ko_targets(prefix)
        if not targets:
            return []
        if level == "gene":
            return list(blocked_reactions_for_genes(self.model, targets))
        return targets

    # Comparison-tab aliases: the tab, its menu entries, and the tests address the picker by
    # these names, so they stay as thin delegates onto the shared implementation.

    def _populate_ko_list(self) -> None:
        if not hasattr(self, "ko_available"):
            return
        self._ko_populate("ko", self.ko_level_combo.currentText())

    def _filter_ko_list(self, text: str) -> None:
        self._ko_filter("ko", text)

    def _add_ko_targets(self, ids) -> None:
        self._ko_add("ko", ids)

    def _remove_ko_targets(self, ids) -> None:
        self._ko_remove("ko", ids)

    def _clear_ko_selected(self) -> None:
        self._ko_clear("ko")

    def _update_ko_selected_label(self) -> None:
        self._ko_update_label("ko")

    def _selected_ko_targets(self) -> list[str]:
        return self._ko_targets("ko")

    def _style_all_combo_popups(self) -> None:
        """Centre every combo's popup entries (Qt defaults to left-aligned)."""

        for combo in self.findChildren(QComboBox):
            combo.setItemDelegate(_ComboItemDelegate(combo))

    def _style_disabled_items(self, combo: QComboBox) -> None:
        """Make a combo's disabled entries clearly inert: faint painting + not-allowed cursor."""

        combo.setItemDelegate(_ComboItemDelegate(combo))
        view = combo.view()
        viewport = view.viewport()
        viewport.setMouseTracking(True)
        filt = _DisabledHoverCursor(view)
        viewport.installEventFilter(filt)
        self._cursor_filters.append(filt)

    def _refresh_comparison_templates(self) -> None:
        """Enable LAD/E-Flux2 reference templates only once an omics result exists."""

        combo = getattr(self, "template_combo", None)
        if combo is None:
            return
        for method in ("lad", "eflux2"):
            index = combo.findText(method)
            if index < 0:
                continue
            item = combo.model().item(index)
            if item is None:
                continue
            ready = method in self._omics_computed_methods
            item.setEnabled(ready)
            item.setToolTip(
                ""
                if ready
                else f"Compute {method.upper()} on the Omics tab to use it as a "
                "reference template."
            )
        # If the current pick just became unavailable, fall back to a template that always
        # runs. pFBA, not FBA: it is the library default and is deterministic.
        current = combo.currentText()
        if current in ("lad", "eflux2") and current not in self._omics_computed_methods:
            combo.setCurrentText("pfba")

    def _comparison_method_key(self) -> str:
        label = self.comparison_method_combo.currentText()
        if label == "ROOM":
            return "room"
        return "moma_l1" if label == "MOMA (L1)" else "moma_l2"

    def _on_comparison_method_changed(self, _label: str) -> None:
        """The ROOM tolerance pair only means anything when ROOM is the method."""

        combo = getattr(self, "room_use_case_combo", None)
        if combo is not None:
            combo.setEnabled(self._comparison_method_key() == "room")

    def _comparison_expression(self):
        """Expression for a LAD/E-Flux2 template: real if loaded, else synthetic (flagged)."""

        template = self.template_combo.currentText()
        if template not in ("lad", "eflux2"):
            return None, False
        if self._omics_expression is not None:
            return self._omics_expression, False
        import numpy as np

        rng = np.random.default_rng(0)
        return {g.id: float(rng.uniform(1, 100)) for g in self.model.genes}, True

    @staticmethod
    def _ko_label(targets: list[str]) -> str:
        shown = ", ".join(targets[:4])
        return shown if len(targets) <= 4 else f"{shown}, +{len(targets) - 4} more"

    def run_comparison(self) -> None:
        self._comparison_cache = None  # drop any prior result until this run succeeds
        self.cmp_map_btn.setEnabled(False)
        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        level = self.ko_level_combo.currentText()
        targets = self._selected_ko_targets()
        if not targets:
            self.comparison_summary.setText(
                f"Add one or more {level}s to the knockout set (search on the left, then "
                "Add →) before running."
            )
            return
        method_key = self._comparison_method_key()
        room_use_case = self.room_use_case_combo.currentText()
        expression, synthetic_expression = self._comparison_expression()
        ko_label = self._ko_label(targets)

        # Restore the single/multi-KO table shape (a prior batch run reshapes it).
        self.comparison_table.setColumnCount(3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )

        def _compute(worker_model):
            reference = reference_flux(
                worker_model, template, gene_expression=expression
            )
            if level == "gene":
                reaction_ids = blocked_reactions_for_genes(worker_model, targets)
            else:
                reaction_ids = tuple(targets)
            result = knockout_comparison(
                worker_model,
                reference,
                reaction_ids,
                method=method_key,
                room_use_case=room_use_case,
            )
            return reference, result, reaction_ids

        try:
            reference, result, reaction_ids = self._run_model_in_background(
                _compute, label=f"Running {method_label}…"
            )
        except Exception as exc:
            self.comparison_summary.setText(
                f"Comparison failed: {html.escape(str(exc))}"
            )
            return

        synthetic_note = (
            " <span style='color:#b45309'>⚠ synthetic demo expression — load a CSV on the "
            "Omics tab for a data-driven template.</span>"
            if synthetic_expression
            else ""
        )
        blocked_note = (
            f" ({len(reaction_ids)} reactions blocked)" if level == "gene" else ""
        )
        if level == "gene" and not reaction_ids:
            self.comparison_table.setRowCount(0)
            self.comparison_summary.setText(
                f"Knocking out gene(s) {html.escape(ko_label)} blocks no reactions "
                "(no GPR effect) — nothing to compare." + synthetic_note
            )
            self.status_label.setText("Comparison: knockout has no effect.")
            return
        if result.status != "optimal":
            self.comparison_table.setRowCount(0)
            self.comparison_summary.setText(
                f"<b>{method_label}</b> vs <b>{template}</b>: knocking out "
                f"{html.escape(ko_label)}{blocked_note} makes the model <b>infeasible</b> "
                "(e.g. it cannot meet maintenance demand)." + synthetic_note
            )
            self.status_label.setText("Comparison: infeasible perturbation.")
            return
        # Cache the raw solve so the change-threshold can re-filter without re-solving.
        self._comparison_cache = {
            "reference": reference,
            "fluxes": dict(result.fluxes),
            "method_label": method_label,
            "template": template,
            "level": level,
            "ko_label": ko_label,
            "blocked_note": blocked_note,
            "status": result.status,
            "quantity": self._comparison_quantity(result),
            "synthetic_note": synthetic_note,
        }
        self.cmp_map_btn.setEnabled(True)
        self._render_comparison_table()

    @staticmethod
    def _comparison_quantity(result) -> str:
        """The one number this method actually reports, named for what it is.

        0.4.0 split ``ComparisonResult.distance`` from the raw solver objective, because the
        three methods optimise three different quantities: MOMA-L2's objective is Sum d^2
        while Segre et al.'s Eq. (4) distance is its square root, and ROOM's objective is a
        *count of switched reactions*, not a distance at all — so ``distance`` is ``None``
        there. Formatting that ``None`` with ``:.4g`` is what used to raise a ``TypeError``
        on the GUI's own ROOM path.
        """

        if result.distance_kind == "none" or result.distance is None:
            switches = result.n_changed_reactions
            if switches is None:
                return f"objective {result.objective_value:.4g}"
            tolerances = (
                f" at δ={result.metadata.get('delta')} ε={result.metadata.get('epsilon')}, "
                f"the {result.metadata.get('room_use_case')} pair"
                if "delta" in result.metadata
                else ""
            )
            return (
                f"{switches:.0f} reactions switched (ROOM's objective — a count, not a "
                f"distance){tolerances}"
            )
        label = {
            "euclidean_l2": "Euclidean (L2) distance",
            "l1": "L1 distance",
        }.get(result.distance_kind, result.distance_kind)
        if result.distance_kind == "euclidean_l2":
            return (
                f"{label} {result.distance:.4g} (QP objective Σd² "
                f"{result.objective_value:.4g})"
            )
        return f"{label} {result.distance:.4g}"

    def _render_comparison_table(self) -> None:
        """Fill the single-run table with reactions whose change exceeds the threshold."""

        cache = self._comparison_cache
        if not cache:
            return
        reference = cache["reference"]
        fluxes = cache["fluxes"]
        delta = self.comparison_threshold_spin.value() / 100.0
        epsilon = (
            1e-3  # absolute floor so reactions near zero reference aren't over-reported
        )
        rows = sorted(fluxes.items(), key=lambda kv: -abs(kv[1] - reference.get(kv[0])))
        changed = [
            (r, v)
            for r, v in rows
            if abs(v - reference.get(r)) > delta * abs(reference.get(r)) + epsilon
        ]
        self.comparison_table.setColumnCount(3)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Reaction", "Reference flux", "Perturbed flux"]
        )
        self.comparison_table.setRowCount(len(changed))
        for i, (rid, flux) in enumerate(changed):
            self.comparison_table.setItem(i, 0, QTableWidgetItem(rid))
            self.comparison_table.setItem(
                i, 1, QTableWidgetItem(f"{reference.get(rid):.4g}")
            )
            self.comparison_table.setItem(i, 2, QTableWidgetItem(f"{flux:.4g}"))
        self.comparison_summary.setText(
            f"<b>{cache['method_label']}</b> vs <b>{cache['template']}</b> template after "
            f"knocking out {cache['level']} {html.escape(cache['ko_label'])}"
            f"{cache['blocked_note']}: status {cache['status']}, "
            f"{cache['quantity']}, {len(changed)} reactions changed "
            f"(≥ {self.comparison_threshold_spin.value():g}% of reference)."
            + cache["synthetic_note"]
        )
        self.status_label.setText(
            f"Comparison complete ({cache['method_label']}, {cache['template']})."
        )

    def run_batch_comparison(self) -> None:
        """Run MOMA/ROOM once per target (each a separate single knockout) into a batch table."""

        import math

        # The batch reshapes the table; drop the single-run cache so the change-threshold
        # spinbox can't re-render a stale single-run over the batch table.
        self._comparison_cache = None
        self.cmp_map_btn.setEnabled(False)

        method_label = self.comparison_method_combo.currentText()
        template = self.template_combo.currentText()
        level = self.ko_level_combo.currentText()
        selected = self._selected_ko_targets()  # empty -> all targets of this level
        method_key = self._comparison_method_key()
        # A screen asks which deletions are lethal, so the library defaults this path to the
        # lethality pair; the GUI honours whatever the selector says and reports it.
        room_use_case = self.room_use_case_combo.currentText()
        expression, synthetic_expression = self._comparison_expression()
        product = self.batch_product_combo.currentText()
        product = None if product in ("", "(none)") else product
        objective_id = next(
            (r.id for r in self.model.reactions if r.objective_coefficient != 0), None
        )

        def _compute(worker_model):
            reference = reference_flux(
                worker_model, template, gene_expression=expression
            )
            if level == "gene":
                perts = gene_perturbations(worker_model, selected or None)
            else:
                perts = reaction_perturbations(worker_model, selected or None)
            rows = batch_comparison(
                worker_model,
                reference,
                perts,
                method=method_key,
                room_use_case=room_use_case,
                product_reaction=product,
            )
            return reference, rows

        self.status_label.setText(f"Running batch {method_label} ({level})…")
        try:
            reference, rows = self._run_model_in_background(
                _compute, label=f"Batch {method_label} over {level}s…"
            )
        except Exception as exc:
            self.comparison_summary.setText(
                f"Batch comparison failed: {html.escape(str(exc))}"
            )
            self.status_label.setText("Batch comparison failed.")
            return

        # Wild-type reference values (constant across knockouts).
        wt_biomass = (
            float(reference.get(objective_id)) if objective_id else float("nan")
        )
        wt_target = float(reference.get(product)) if product else float("nan")

        def _essential(row) -> bool:
            """A knockout is essential if it leaves the cell unable to grow (or infeasible)."""
            if row.status != "optimal" or not math.isfinite(row.objective):
                return True
            if math.isfinite(wt_biomass) and wt_biomass > 1e-9:
                return row.objective < 0.01 * wt_biomass  # under 1% of wild-type growth
            return row.objective < 1e-6

        # Most impactful first: lowest post-knockout growth (essential/lethal) at the top.
        rows = sorted(
            rows,
            key=lambda r: r.objective if math.isfinite(r.objective) else float("-inf"),
        )

        headers = ["Target", "WT Biomass", "KO Biomass", "Essential"]
        if product:
            headers += ["WT target", "KO target"]
        self.comparison_table.setColumnCount(len(headers))
        self.comparison_table.setHorizontalHeaderLabels(headers)
        self.comparison_table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            ko_biomass = "—" if not math.isfinite(r.objective) else f"{r.objective:.4g}"
            cells = [
                r.target_id,
                f"{wt_biomass:.4g}" if math.isfinite(wt_biomass) else "—",
                ko_biomass,
                "yes" if _essential(r) else "no",
            ]
            if product:
                ko_target = (
                    "—"
                    if not math.isfinite(r.product_flux)
                    else f"{r.product_flux:.4g}"
                )
                cells += [
                    f"{wt_target:.4g}" if math.isfinite(wt_target) else "—",
                    ko_target,
                ]
            for col, text in enumerate(cells):
                self.comparison_table.setItem(i, col, QTableWidgetItem(text))

        synthetic_note = (
            " <span style='color:#b45309'>⚠ synthetic demo expression.</span>"
            if synthetic_expression
            else ""
        )
        essential_n = sum(1 for r in rows if _essential(r))
        product_note = f" Target product: {html.escape(product)}." if product else ""
        self.comparison_summary.setText(
            f"<b>Batch {method_label}</b> vs <b>{template}</b> over {len(rows)} {level} "
            f"knockout(s): {essential_n} essential (growth lost). Sorted by post-knockout "
            f"biomass (most impactful first).{product_note}" + synthetic_note
        )
        self.status_label.setText(
            f"Batch comparison complete ({method_label}, {level})."
        )

    def _build_simulation_tab(self) -> QWidget:
        tab = QWidget()
        layout = QVBoxLayout(tab)
        media_row = QHBoxLayout()
        media_row.addWidget(QLabel("Medium:"))
        self.medium_combo = QComboBox()
        self.medium_combo.addItems(sorted(PRESET_MEDIA))
        self.medium_combo.setToolTip(
            "Preset growth medium (sets exchange uptake bounds)"
        )
        # The preset keys are stable, but the display names changed in 0.4.0 — they now state
        # the mineral set's provenance (iJO1366/BiGG, not M9) and that CO2 uptake is closed.
        # Show each one on its own item rather than letting the key imply a standard medium.
        for _index, _key in enumerate(sorted(PRESET_MEDIA)):
            self.medium_combo.setItemData(
                _index, PRESET_MEDIA[_key].name, Qt.ToolTipRole
            )
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
        self.fva_fraction_spin.setToolTip(
            "FVA: fraction of the optimum objective to hold"
        )
        controls.addWidget(self.fva_fraction_spin)
        # Objective on top, pFBA's minimal-total-flux directly beneath it (previously only in
        # the easy-to-miss status bar).
        obj_box = QVBoxLayout()
        obj_box.setSpacing(0)
        self.objective_label = QLabel("Objective: —")
        self.objective_label.setFont(QFont("", 13, QFont.Bold))
        self.pfba_total_label = QLabel("")
        self.pfba_total_label.setStyleSheet("color: #5c6b7e; font-size: 11px;")
        obj_box.addWidget(self.objective_label)
        obj_box.addWidget(self.pfba_total_label)
        controls.addStretch(1)
        controls.addLayout(obj_box)
        layout.addLayout(controls)

        # FBA and pFBA fluxes live in separate columns so running pFBA adds to — rather than
        # overwrites — the FBA result and the two can be compared side by side.
        self.sim_table = QTableWidget(0, 4)
        self.sim_table.setHorizontalHeaderLabels(
            ["Reaction", "FBA flux", "pFBA flux", "FVA range"]
        )
        _sim_header = self.sim_table.horizontalHeader()
        _sim_header.setSectionResizeMode(0, QHeaderView.Stretch)
        for _c in (1, 2, 3):
            _sim_header.setSectionResizeMode(_c, QHeaderView.ResizeToContents)
        self.sim_table.verticalHeader().setVisible(False)
        self.sim_table.setAlternatingRowColors(True)
        layout.addWidget(self.sim_table, 1)
        return tab

    @staticmethod
    def _stage_box(title: str, explanation: str) -> tuple[QGroupBox, QFormLayout]:
        """A parameter group that says what stage it configures before listing its fields.

        The revert tab drives a three-stage pipeline and its parameters mean nothing out of
        that context — 'epsilon' and 'alpha' belong to different stages and answer different
        questions. Each box therefore opens with a sentence naming its stage's job.
        """

        box = QGroupBox(title)
        outer = QVBoxLayout(box)
        # Four of these stack in one column, so the padding is set against the height budget
        # rather than left at the default: generous enough that the grouping reads as grouping,
        # tight enough that the column still fits a laptop window without scrolling.
        outer.setContentsMargins(12, 6, 12, 6)
        outer.setSpacing(5)
        caption = QLabel(explanation)
        caption.setWordWrap(True)
        caption.setStyleSheet("color: #5a6b80; font-size: 11px;")
        outer.addWidget(caption)
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.setVerticalSpacing(5)
        outer.addLayout(form)
        return box, form

    @staticmethod
    def _align_stage_forms(forms: Sequence[QFormLayout]) -> None:
        """Give a column of stage boxes one shared label column, right-aligned to its fields.

        Each ``QFormLayout`` sizes its label column to its own longest label, so four stacked
        boxes put their inputs at four different x positions, and inside one box a short label
        is left stranded against the box edge with its field far to the right. Widening every
        label to the widest in the column lines the inputs up down the whole column, and
        right-aligning keeps each colon beside the field it names.
        """

        labels = []
        for form in forms:
            for row in range(form.rowCount()):
                item = form.itemAt(row, QFormLayout.LabelRole)
                if item is None or item.widget() is None:
                    # A spanning row - the replicate note - reports through FieldRole too, and
                    # it is a wrapped paragraph that wants the whole box, not an input.
                    continue
                labels.append(item.widget())
                field = form.itemAt(row, QFormLayout.FieldRole)
                if field is not None and field.widget() is not None:
                    field.widget().setMaximumWidth(_FIELD_MAX_WIDTH)
        if not labels:
            return
        widest = 0
        for label in labels:
            # A label reports the *unstyled* font's width until it has been polished against
            # the window it lives in - 190px rather than 134px for the longest of these - so
            # measuring any earlier reserves a column half again wider than the text needs.
            label.ensurePolished()
            metrics = QFontMetrics(label.font())
            widest = max(widest, metrics.horizontalAdvance(label.text()) + _LABEL_PAD)
        for label in labels:
            # setFixedWidth, not setMinimumWidth. A minimum is only a floor: a label that
            # renders wider than this still widens its own QFormLayout's label column, and
            # that stage's inputs then sit right of every other stage's. Measuring from the
            # resolved font rather than from sizeHint is what makes the pin safe to apply.
            label.setFixedWidth(widest)
            # The alignment has to be set on the label's own text, not through the layout's
            # setLabelAlignment: that positions the label *widget* inside the column, and
            # every widget is now exactly the column's width.
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

    def _build_revert_tab(self) -> QWidget:
        tab = QWidget()
        # Parameters left, results right: the three stages are read once while setting a run
        # up, whereas the ranking is what the user returns to, and stacking them vertically
        # pushed the table off the bottom of the window on a laptop screen.
        columns = QHBoxLayout(tab)
        left = QWidget()
        layout = QVBoxLayout(left)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(20)

        intro = QLabel(
            "Ranks knockouts by how well they move a <b>source</b> state toward a <b>target</b> one — "
            "a disease state back toward a healthy one, say. Hover any field for its published "
            "value and what it controls."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        # --- inputs -----------------------------------------------------------------------
        inputs_box, inputs_form = self._stage_box(
            "Input — two expression states",
            "CSV or TSV: a gene id column, then one column per replicate. Values are linear "
            "expression (counts, TPM, FPKM), not log-transformed. Two or more replicates on "
            "each side unlock the published t-test. Source = the state to move away from.",
        )
        source_row = QHBoxLayout()
        source_btn = QPushButton("Load source CSV/TSV…")
        source_btn.clicked.connect(lambda: self.load_revert_expression("source"))
        self.revert_source_label = QLabel("not loaded")
        source_row.addWidget(source_btn)
        source_row.addWidget(self.revert_source_label, 1)
        inputs_form.addRow("Source expression:", source_row)

        target_row = QHBoxLayout()
        target_btn = QPushButton("Load target CSV/TSV…")
        target_btn.clicked.connect(lambda: self.load_revert_expression("target"))
        self.revert_target_label = QLabel("not loaded")
        target_row.addWidget(target_btn)
        target_row.addWidget(self.revert_target_label, 1)
        inputs_form.addRow("Target expression:", target_row)
        layout.addWidget(inputs_box)

        # --- stage 1: source reference state ----------------------------------------------
        vref_box, vref_form = self._stage_box(
            "Stage 1 — source flux state (v_ref)",
            "Turns the source expression into fluxes. Everything below is measured against it.",
        )
        self.revert_reference_combo = QComboBox()
        self.revert_reference_combo.addItems(["eflux2", "lad"])
        self.revert_reference_combo.setToolTip(
            "eflux2 (Kim et al. 2016): scales each reaction's bounds by its expression, then "
            "minimises total flux.\n"
            "lad (Lee et al. 2012): fits |flux| to expression-derived targets by least "
            "absolute deviation; needs only an LP solver.\n\n"
            "Yizhak et al. use iMAT, which CMM does not implement. iMAT places no objective "
            "on growth, so the ranking is conditioned on whichever estimator is used. Supply "
            "an external state through the Python API to reproduce the published pipeline."
        )
        vref_form.addRow("Estimator:", self.revert_reference_combo)
        layout.addWidget(vref_box)

        # --- stage 2: direction map -------------------------------------------------------
        direction_box, direction_form = self._stage_box(
            "Stage 2 — which reactions should change",
            "Compares the two states through each reaction's GPR. Reactions left out must stay put.",
        )
        self.revert_fold_change_spin = QDoubleSpinBox()
        self.revert_fold_change_spin.setDecimals(2)
        self.revert_fold_change_spin.setRange(0.0, 20.0)
        self.revert_fold_change_spin.setSingleStep(0.1)
        self.revert_fold_change_spin.setValue(1.0)
        self.revert_fold_change_spin.setToolTip(
            "A gene counts as changed when |log2 fold change| reaches this. 1.0 means a "
            "two-fold change.\n\n"
            "Used only when the significance test above is 'fold change'. Yizhak et al. use a "
            "Student's t-test at P < 0.05 instead, which this tab runs whenever both files "
            "carry at least two replicate columns."
        )
        self.revert_top_changed_spin = QSpinBox()
        self.revert_top_changed_spin.setRange(0, 100000)
        self.revert_top_changed_spin.setSpecialValueText("all (no cut)")
        self.revert_top_changed_spin.setValue(0)
        self.revert_top_changed_spin.setToolTip(
            "Keep only this many of the most differentially expressed changed reactions; the "
            "rest are treated as steady. Yizhak et al. use the top 100-200. Each changed "
            "reaction adds one binary variable to the search, so on a genome-scale model this "
            "cut is also what decides whether the run finishes. 0 keeps every reaction that "
            "clears the fold-change threshold."
        )
        self.revert_significance_combo = QComboBox()
        self.revert_significance_combo.addItem("t-test (P value)", "ttest")
        self.revert_significance_combo.addItem("fold change", "fold_change")
        self.revert_significance_combo.setToolTip(
            "t-test: Yizhak et al.'s (2013) step 2 — Student's t-test per gene across the "
            "replicate columns, split three ways at the P cutoff. This is the published "
            "route and needs at least two replicates in each file.\n"
            "fold change: cut on |log2 fold change| of the column means. All a single "
            "measurement per state supports; say so when reporting a result.\n\n"
            "Whichever is used also orders the 'top changed reactions' cut below: smallest "
            "gene P value for the t-test, largest |log2FC| for fold change."
        )
        self.revert_p_value_spin = QDoubleSpinBox()
        self.revert_p_value_spin.setDecimals(4)
        self.revert_p_value_spin.setRange(0.0001, 1.0)
        self.revert_p_value_spin.setSingleStep(0.01)
        self.revert_p_value_spin.setValue(0.05)
        self.revert_p_value_spin.setToolTip(
            "A gene counts as changed when its two-sided Student's t-test P value falls "
            "below this. Yizhak et al. use P < 0.05.\n\n"
            "The test is run per gene with no multiple-testing correction, as the paper "
            "specifies; the 'top changed reactions' cut is what limits the changed set."
        )
        self.revert_significance_combo.currentIndexChanged.connect(
            self._on_revert_significance_changed
        )
        direction_form.addRow("Significance test:", self.revert_significance_combo)
        direction_form.addRow("t-test P cutoff:", self.revert_p_value_spin)
        direction_form.addRow("Fold-change threshold:", self.revert_fold_change_spin)
        direction_form.addRow("Top changed reactions:", self.revert_top_changed_spin)
        self.revert_replicate_note = QLabel("")
        self.revert_replicate_note.setWordWrap(True)
        # A wrapped label spanning a form row reports one line's height unless its size policy
        # says the height depends on the width. Without this the note is clipped mid-sentence
        # and spills past the box border - and only once the row is wide enough to wrap.
        note_policy = self.revert_replicate_note.sizePolicy()
        note_policy.setHeightForWidth(True)
        self.revert_replicate_note.setSizePolicy(note_policy)
        self.revert_replicate_note.setStyleSheet("color: #5a6b80; font-size: 11px;")
        direction_form.addRow(self.revert_replicate_note)
        layout.addWidget(direction_box)

        # --- stage 3: the search ----------------------------------------------------------
        search_box, search_form = self._stage_box(
            "Stage 3 — knockout search",
            "One optimisation per candidate, ranked by how much of the required change it achieved.",
        )
        self.method_combo = QComboBox()
        self.method_combo.addItems(["rmta", "mta", "rmta_continuous"])
        self.method_combo.setToolTip(
            "rmta (Valcárcel et al. 2019): three solves per candidate — best case, MOMA, and "
            "the direction reversed for the worst case — combined so that a candidate scoring "
            "well whichever way it is pushed is demoted. MIQP; about 3x the cost of mta.\n"
            "mta (Yizhak et al. 2013): the single published MIQP.\n"
            "rmta_continuous: a historical QP heuristic. It is not published rMTA and must "
            "never be reported as such."
        )
        self.perturbation_combo = QComboBox()
        self.perturbation_combo.addItems(["gene", "reaction"])
        self.perturbation_combo.setToolTip(
            "gene: knock out each gene and block whatever reactions its GPR disables — closer "
            "to an experiment.\n"
            "reaction: knock out each reaction directly."
        )
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0.0, 1.0)
        self.alpha_spin.setSingleStep(0.01)
        self.alpha_spin.setValue(0.66)
        self.alpha_spin.setToolTip(
            "Balances the two halves of the objective: higher a rewards making the required "
            "changes, lower a rewards leaving the rest undisturbed.\n\n"
            "Yizhak et al. report a = 0.66 in the main text and find the ranking robust over "
            "0.1-0.9."
        )
        self.revert_epsilon_spin = QDoubleSpinBox()
        self.revert_epsilon_spin.setDecimals(4)
        self.revert_epsilon_spin.setRange(0.0, 1000.0)
        self.revert_epsilon_spin.setSingleStep(0.001)
        self.revert_epsilon_spin.setValue(REVERT_DEFAULT_EPSILON)
        self.revert_epsilon_spin.setToolTip(
            "How far a reaction must move from v_ref before the move counts: v >= v_ref + e to "
            "increase, v <= v_ref - e to decrease. Also floors the score's denominator.\n\n"
            "This is measured in the model's own flux units, so there is no safe default — a "
            "value suited to a core model is not suited to a genome-scale one. Yizhak et al. "
            "derive it per data set from a sampled reference distribution; CMM uses the fixed "
            "value set here, so state it when reporting a result."
        )
        search_form.addRow("Method:", self.method_combo)
        search_form.addRow("Knockout level:", self.perturbation_combo)
        search_form.addRow("Transformation weight α:", self.alpha_spin)
        search_form.addRow("Significant flux change ε:", self.revert_epsilon_spin)

        layout.addWidget(search_box)
        layout.addStretch(1)

        # Four stacked groups do not fit a laptop window, so the stages scroll — but the Run
        # button sits outside that scroll area, pinned under it. Inside, it would be the last
        # thing in the tallest column and therefore the one control never visible on arrival.
        left_scroll = QScrollArea()
        left_scroll.setWidget(left)
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QFrame.NoFrame)
        left_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        left_column = QWidget()
        left_layout = QVBoxLayout(left_column)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(left_scroll, 1)

        self.revert_run_btn = QPushButton("Run Revert")
        self.revert_run_btn.setEnabled(False)
        self.revert_run_btn.clicked.connect(self.run_loaded_revert)
        self.revert_run_btn.setToolTip(
            "Load both source and target expression files to enable."
        )
        left_layout.addWidget(self.revert_run_btn)

        results = QWidget()
        right = QVBoxLayout(results)
        right.setContentsMargins(0, 0, 0, 0)
        self.revert_summary = QLabel(
            "Load source and target expression files to rank normalization targets."
        )
        self.revert_summary.setWordWrap(True)
        right.addWidget(self.revert_summary)

        self.revert_table = QTableWidget(0, 3)
        self.revert_table.setHorizontalHeaderLabels(["Rank", "Target", "Score"])
        header = self.revert_table.horizontalHeader()
        # Rank is at most four digits; giving it an equal third of a wide table wastes the
        # room the identifiers and scores can use.
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        self.revert_table.verticalHeader().setVisible(False)
        self.revert_table.setAlternatingRowColors(True)
        right.addWidget(self.revert_table, 1)

        columns.addWidget(left_column, 1)
        columns.addWidget(results, 1)
        # Aligned at the end of _build_ui rather than here: the labels have to be inside the
        # window before their widths mean anything.
        self._revert_stage_forms = [inputs_form, vref_form, direction_form, search_form]
        # Nothing is loaded yet, so this settles the direction controls into their
        # can't-t-test state and prints the reason rather than leaving them looking usable.
        self._update_revert_direction_controls()
        return tab

    # -- behavior -----------------------------------------------------------

    def load_model(self, model: Model) -> None:
        self.model = model
        # A map belongs to a model, so re-choose one; a user-loaded map that still fits is
        # kept, and one that no longer draws anything is dropped rather than left stale.
        if self._map_path is None or not (
            map_reaction_ids(self._map_path) & {rxn.id for rxn in model.reactions}
        ):
            self._map_path = bundled_map_for(model)
            self._clear_map_background()
        self._refresh_map_source()
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
        self.sd_product_combo.clear()
        self.sd_product_combo.addItems(exchanges)
        self.batch_product_combo.clear()
        self.batch_product_combo.addItems(["(none)", *exchanges])  # optional in batch
        self.substrate_combo.clear()
        self.substrate_combo.addItems(["auto", *exchanges])

        # Flux response scans any reaction, not just exchanges, so it gets the full list.
        reaction_ids = sorted(r.id for r in model.reactions)
        self.fr_target_combo.clear()
        self.fr_target_combo.addItems(reaction_ids)
        self.fr_response_combo.clear()
        self.fr_response_combo.addItems(["(objective)", *reaction_ids])
        if exchanges:
            self.fr_target_combo.setCurrentText(exchanges[0])

        self._populate_ko_list()
        self._ko_populate("sample_ko", self.sample_ko_level_combo.currentText())
        self._loading = False
        # Seed the flux-response range for the initial target; the combo was populated while
        # _loading suppressed the change handler.
        self._prefill_fr_range()
        if self._default_product and self._default_product in exchanges:
            self.product_combo.setCurrentText(self._default_product)
            self.sd_product_combo.setCurrentText(self._default_product)

        # Production and strain design need exchange reactions; disable them otherwise.
        has_exchanges = bool(exchanges)
        for btn in getattr(self, "_production_buttons", []):
            btn.setEnabled(has_exchanges)
        if hasattr(self, "sd_run_btn"):
            self.sd_run_btn.setEnabled(has_exchanges)
        if not has_exchanges:
            self.yield_label.setText(
                "No exchange reactions in this model — production design unavailable."
            )

        # The agent tab reads its exchange lists and its run gate from the loaded model.
        # Guarded like the lines above: load_model also runs before the UI is fully built.
        if hasattr(self, "jev_product_combo"):
            self._refresh_jev_inputs()

        # A sampled ensemble belongs to the model it was drawn from; drop it on a reload so a
        # stale export cannot be attributed to the new model.
        self._sampling_result = None
        if hasattr(self, "sample_export_btn"):
            self.sample_export_btn.setEnabled(False)
            self.sample_table.setRowCount(0)
            self.fr_table.setRowCount(0)

        self._fluxes = {}
        self._fba_fluxes = {}
        self._pfba_fluxes = {}
        self._flux_source = ""
        self._fluxes_stale = False
        if getattr(self, "pfba_total_label", None) is not None:
            self.pfba_total_label.setText("")
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
            self.status_label.setText(
                "Bounds changed — fluxes are stale; click Run FBA."
            )

    def run_fba(self) -> None:
        try:
            solution = fba(self.model)
        except Exception as exc:  # surface solver/model errors instead of freezing
            self.objective_label.setText("Objective: error")
            self.status_label.setText(f"FBA failed: {exc}")
            return
        self._fluxes = dict(solution.fluxes)
        self._fba_fluxes = dict(solution.fluxes)
        self._flux_source = "FBA"
        self._fluxes_stale = False
        obj = solution.objective_value
        obj_text = f"{obj:.4g}" if obj is not None else "infeasible"
        self.objective_label.setText(f"Objective: {obj_text} ({solution.status})")
        self.pfba_total_label.setText("")  # FBA does not minimise total flux
        if obj is None or solution.status != "optimal":
            self.status_label.setText(
                f"FBA status: {solution.status} — model may be infeasible."
            )
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
        self._pfba_fluxes = dict(solution.fluxes)
        self._flux_source = "pFBA"
        self._fluxes_stale = False
        # pFBA's objective value is the minimal total flux; show the growth (objective rxn).
        growth = next(
            (
                self._fluxes[r.id]
                for r in self.model.reactions
                if r.objective_coefficient != 0
            ),
            None,
        )
        growth_text = f"{growth:.4g}" if growth is not None else "—"
        total = sum(abs(v) for v in self._fluxes.values())
        self.objective_label.setText(f"Objective: {growth_text} (pFBA)")
        self.pfba_total_label.setText(f"Minimal total flux: {total:.1f}")
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
        # `apply_to` returns a MediumApplication as of 0.4.0: it still behaves as the mapping
        # of applied uptakes, but it also names the components this model has no exchange for.
        # A preset that silently drops 18 of 24 components is not the medium it claims to be.
        dropped = (
            f" {len(applied.dropped)} declared component(s) have no exchange in this model "
            f"and were dropped: {', '.join(sorted(applied.dropped))}."
            if applied.dropped
            else ""
        )
        self.status_label.setText(
            f"Applied '{applied.medium}' ({len(applied)} open exchanges).{dropped} "
            "Re-run FBA."
        )
        self.medium_combo.setToolTip(
            "Preset growth medium (sets exchange uptake bounds). Currently applied: "
            f"{applied.medium}"
        )

    def run_fva(self) -> None:
        if not self._fluxes or self._fluxes_stale:
            self.run_fba()
        fraction = self.fva_fraction_spin.value()
        try:
            ranges = self._run_model_in_background(
                lambda worker_model: fva(worker_model, fraction_of_optimum=fraction),
                label="Running FVA…",
            )
        except Exception as exc:
            self.status_label.setText(f"FVA failed: {exc}")
            return
        self._fill_sim_table(fva_ranges=ranges, fraction=fraction)
        self.status_label.setText(f"FVA complete (fraction {fraction:g}).")

    def _fill_sim_table(self, *, fva_ranges=None, fraction=None) -> None:
        """Fill the simulation table with FBA/pFBA flux columns (and an FVA range if given)."""

        fva_header = "FVA range" if fraction is None else f"FVA range (f={fraction:g})"
        self.sim_table.setHorizontalHeaderLabels(
            ["Reaction", "FBA flux", "pFBA flux", fva_header]
        )
        ids = sorted(set(self._fba_fluxes) | set(self._pfba_fluxes))
        self.sim_table.setRowCount(len(ids))
        for row, rid in enumerate(ids):
            self.sim_table.setItem(row, 0, QTableWidgetItem(rid))
            fba = self._fba_fluxes.get(rid)
            pfba = self._pfba_fluxes.get(rid)
            self.sim_table.setItem(
                row, 1, QTableWidgetItem(f"{fba:.3g}" if fba is not None else "—")
            )
            self.sim_table.setItem(
                row, 2, QTableWidgetItem(f"{pfba:.3g}" if pfba is not None else "—")
            )
            if fva_ranges is not None:
                rng = fva_ranges.get(rid)
                text = f"[{rng.minimum:.3g}, {rng.maximum:.3g}]" if rng else "—"
            else:
                text = "—"
            self.sim_table.setItem(row, 3, QTableWidgetItem(text))

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
            color = QColor("#dce9d6") if active and value > 0 else QColor("white")
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

    def _aeration_condition(self, aerobic: bool) -> Condition:
        """The selected aeration as a :class:`Condition` — 0.4.0's only way to state it.

        ``aerobic=True|False`` was removed from ``cmm.features.production`` in 0.4.0, so the
        GUI now says what every other caller says. The anaerobic condition closes CO2 uptake
        alongside oxygen, which is what the removed flag did internally and what the
        ``glucose_anaerobic`` preset does through the media layer; without it an anaerobic
        succinate yield reads 1.3906 mol/mol instead of 1.2000, inflated by CO2 a closed
        fermentation cannot supply. Bounds name only exchanges this model actually has.
        """

        if aerobic:
            # Aerobic is the model as conditioned by the medium: no bound overrides, but a
            # named condition so the run records what it was asked for.
            return Condition(
                name="aerobic (model medium as applied)",
                notes="oxygen and CO2 exchanges left as the medium set them",
            )
        bounds = tuple(
            ReactionBound(reaction_id=rid, lower_bound=0.0)
            for rid in ("EX_o2_e", "EX_co2_e")
            if rid in self.model.reactions
        )
        return Condition(
            name="anaerobic (O2 and CO2 uptake closed)",
            bounds=bounds,
            notes="oxygen uptake closed; CO2 uptake closed so the yield is not CO2-inflated",
        )

    def _production_condition(self) -> Condition:
        """The Production tab's aeration selector as a Condition."""

        return self._aeration_condition(self._is_aerobic())

    def _strain_design_condition(self) -> Condition:
        """The Strain Design tab's own aeration selector as a Condition.

        Strain design used to accept no condition at all, so a caller working anaerobically
        received an aerobic design in silence — the defect that put an aerobic succinate
        design into the scenario documentation as the anaerobic answer.
        """

        return self._aeration_condition(
            self.sd_anaerobic_combo.currentText() == "aerobic"
        )

    def _clear_map_background(self) -> None:
        """Drop the loaded drawing and hide its toggle."""

        self._map_background = None
        if hasattr(self, "map_background_check"):
            self.map_background_check.setVisible(False)

    def open_map_background_dialog(self) -> None:
        """Load a drawing of the current map to lay under the flux colouring."""

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open a drawing of this map",
            "",
            "Map drawings (*.svg *.png *.jpg *.jpeg);;All files (*)",
        )
        if not path:
            return
        try:
            if path.lower().endswith(".svg"):
                background = svg_background(path)
            else:
                import matplotlib.image as mpimg

                background = mpimg.imread(path)
        except Exception as exc:
            self.status_label.setText(
                f"Could not read that drawing: {html.escape(str(exc))}"
            )
            return
        self._map_background = background
        self.map_background_check.setVisible(True)
        self.map_background_check.setChecked(True)
        self.status_label.setText(
            f"Background loaded from {html.escape(Path(path).name)}."
        )
        if self._map_canvas is not None:
            self.render_flux_map()

    def _map_line_width(self) -> float | None:
        """``None`` to scale stroke width with flux, or the one width to draw everything at."""

        return None if self.map_width_check.isChecked() else 2.0

    def _draw_simulation_on_map(self, solve) -> None:
        """Run a simulation and draw it, so the map is one click from either method."""

        solve()
        if self._fluxes:  # a failed solve has already reported why
            self.render_flux_map()

    def _refresh_omics_map_label(self) -> None:
        """Say which of the computed conditions the box is showing."""

        total = self.omics_map_combo.count()
        if total > 1:
            self.omics_map_label.setText(
                f"Draw ({self.omics_map_combo.currentIndex() + 1} of {total}):"
            )
        else:
            self.omics_map_label.setText("Draw:")

    def show_omics_on_flux_map(self) -> None:
        """Draw one computed condition's expression-derived flux state."""

        if not self._omics_fluxes_by_condition:
            self.status_label.setText("Compute a condition on the Omics tab first.")
            return
        condition = self.omics_map_combo.currentText()
        fluxes = self._omics_fluxes_by_condition.get(condition)
        if not fluxes:
            self.status_label.setText(f"No computed flux state for '{condition}'.")
            return
        self.show_on_flux_map(fluxes, f"{self._omics_method_label} · {condition}")

    def show_comparison_on_flux_map(self) -> None:
        """Draw the perturbed flux distribution from the last MOMA/ROOM run."""

        cache = self._comparison_cache
        if not cache:
            self.status_label.setText("Run a single MOMA/ROOM comparison first.")
            return
        self.show_on_flux_map(
            cache["fluxes"], f"{cache['method_label']} · {cache['ko_label']}"
        )

    def show_on_flux_map(self, fluxes: dict[str, float], source: str) -> None:
        """Make ``fluxes`` the window's current distribution and draw it on the map.

        The map used to be reachable only from FBA and pFBA, which left out the two results it
        is most useful for — an expression-derived flux state and a knockout's redistribution.
        Those are exactly the cases where a picture beats a table of 95 numbers.
        """

        if not fluxes:
            self.status_label.setText("Nothing to show on the flux map.")
            return
        self._fluxes = dict(fluxes)
        self._flux_source = source
        self._fluxes_stale = False
        self._populate_flux_column()  # the reaction table shows this distribution too
        self._goto_tab("Flux Map")
        self.render_flux_map()

    def render_flux_map(self) -> None:
        """Draw the current flux distribution on the selected layout."""

        if not self._fluxes or self._fluxes_stale:
            self.run_fba()
        if not self._fluxes:  # FBA itself failed; run_fba already reported why
            return
        escher = (
            self.map_layout_combo.currentText() == _MAP_LAYOUT_ESCHER
            and self._map_path is not None
        )
        try:
            source = self._flux_source or "flux"
            if escher:
                show_background = (
                    self._map_background is not None
                    and self.map_background_check.isChecked()
                )
                fig = escher_flux_map(
                    self._map_path,
                    self._fluxes,
                    title=f"{self.model.id} — {source}",
                    background=self._map_background if show_background else None,
                    line_width=self._map_line_width(),
                )
                note = f"Escher layout, {Path(self._map_path).stem}"
                if show_background:
                    note += ", over the loaded drawing"
            else:
                top_n = self.map_topn_spin.value()
                fig = network_flux_map(
                    self.model,
                    self._fluxes,
                    top_n=top_n,
                    title=f"{self.model.id} — {source}, top {top_n} reactions by |flux|",
                    line_width=self._map_line_width(),
                )
                note = f"schematic, top {top_n} reactions"
        except Exception as exc:
            self.status_label.setText(f"Flux map failed: {html.escape(str(exc))}")
            return
        self._set_figure(self.map_canvas_holder, "map", fig)
        self.status_label.setText(f"Flux map rendered — {source} ({note}).")

    def _on_product_changed(self, product: str) -> None:
        if not self._loading and product:
            self.yield_label.setText(
                f"Compute the theoretical yield, envelope, or FSEOF targets for "
                f"{html.escape(product)}."
            )

    def _set_figure(self, holder, key: str, fig) -> None:
        """Swap the figure shown in a tab's canvas holder, releasing the previous one.

        ``key`` names the ``_<key>_canvas`` / ``_<key>_toolbar`` attribute pair so every
        plotting tab shares one implementation instead of repeating the teardown dance.
        """

        canvas = getattr(self, f"_{key}_canvas", None)
        toolbar = getattr(self, f"_{key}_toolbar", None)
        if canvas is not None:
            old_fig = canvas.figure
            holder.removeWidget(canvas)
            canvas.setParent(None)
            if toolbar is not None:
                holder.removeWidget(toolbar)
                toolbar.setParent(None)
            old_fig.clf()  # release the previous figure so it is not leaked
        # Figures are authored at 300 DPI for file export; render on screen at a lower DPI
        # so the canvas fits the panel instead of clipping.
        fig.set_dpi(100)
        canvas = FigureCanvas(fig)
        # A matplotlib toolbar gives interactive zoom/pan and a "Save the figure" button.
        toolbar = NavigationToolbar(canvas, self)
        holder.addWidget(toolbar)
        holder.addWidget(canvas)
        canvas.draw()
        setattr(self, f"_{key}_canvas", canvas)
        setattr(self, f"_{key}_toolbar", toolbar)

    def _set_production_figure(self, fig) -> None:
        self._set_figure(self.production_canvas_holder, "production", fig)

    def _set_production_table_headers(self, labels: list[str], tooltips=None) -> None:
        """Re-shape the shared production results table for whichever analysis just ran."""

        self.fseof_table.setColumnCount(len(labels))
        self.fseof_table.setHorizontalHeaderLabels(labels)
        self.fseof_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        if tooltips:
            for i, tip in enumerate(tooltips):
                header_item = self.fseof_table.horizontalHeaderItem(i)
                if header_item is not None and tip:
                    header_item.setToolTip(tip)

    def _rerender_production_table(self) -> None:
        """Re-render the cached FSEOF/FVSEOF result when the show-all toggle changes."""

        cached = self._production_result
        if not cached:
            return
        kind, result = cached
        if kind == "fseof":
            self._fill_fseof_table(result)
        elif kind == "fvseof":
            self._fill_fvseof_table(result)

    def _fill_fseof_table(self, result) -> None:
        """Show FSEOF targets ranked by change in |flux|; all reactions when 'show all' is on."""

        self._production_result = ("fseof", result)
        show_all = self.production_show_all_check.isChecked()
        levels = list(result.enforced_levels)
        lo, hi = levels[0], levels[-1]
        rows = []
        for rid in result.trends.index:
            cls = str(result.trends.loc[rid, "classification"])
            if not show_all and cls not in ("amplify", "knockdown"):
                continue
            flo = float(result.trends.loc[rid, lo])
            fhi = float(result.trends.loc[rid, hi])
            rows.append((rid, flo, fhi, cls, abs(abs(fhi) - abs(flo))))
        rows.sort(key=lambda r: -r[4])
        self._set_production_table_headers(
            ["Reaction", "Flux @ low enforce", "Flux @ high enforce", "Trend"],
            [
                "Reaction ID",
                f"Reaction's flux when the product exchange is fixed at the lowest scanned "
                f"level ({lo:.3g}, ~10% of the theoretical maximum) and biomass is maximised",
                f"Reaction's flux when the product exchange is fixed at the highest scanned "
                f"level ({hi:.3g}, ~90% of the theoretical maximum) and biomass is maximised",
                "amplify = |flux| rises as the cell is pushed toward product (over-express "
                "target); knockdown = |flux| falls (down-regulate/delete target)",
            ],
        )
        self.fseof_table.setRowCount(len(rows))
        for i, (rid, flo, fhi, cls, _delta) in enumerate(rows):
            self.fseof_table.setItem(i, 0, QTableWidgetItem(rid))
            self.fseof_table.setItem(i, 1, QTableWidgetItem(f"{flo:.4g}"))
            self.fseof_table.setItem(i, 2, QTableWidgetItem(f"{fhi:.4g}"))
            self.fseof_table.setItem(i, 3, QTableWidgetItem(cls))

    def _fill_fvseof_table(self, result) -> None:
        """Show FVSEOF targets with Park's type; all reactions when 'show all' is on."""

        self._production_result = ("fvseof", result)
        show_all = self.production_show_all_check.isChecked()
        levels = list(result.enforced_levels)
        lo, hi = levels[0], levels[-1]
        robust_ids = set(result.robust_targets())
        rows = []
        for rid in result.classification.index:
            cls = str(result.classification.loc[rid])
            if not show_all and cls not in ("amplify", "knockdown"):
                continue
            mlo = float(result.mean.loc[rid, lo])
            mhi = float(result.mean.loc[rid, hi])
            is_robust = rid in robust_ids
            tag = f"{cls} · robust" if is_robust else cls
            park = f"{int(result.park_type.loc[rid])}"
            rows.append((rid, mlo, mhi, park, tag, is_robust, abs(abs(mhi) - abs(mlo))))
        rows.sort(key=lambda r: (-int(r[5]), -r[6]))
        self._set_production_table_headers(
            ["Reaction", "Mean flux @ low", "Mean flux @ high", "Park type", "Trend"],
            [
                "Reaction ID",
                f"V_avg — midpoint of the reaction's FVA flux range at the lowest enforced "
                f"product level ({lo:.3g})",
                f"V_avg — midpoint of the reaction's FVA flux range at the highest enforced "
                f"product level ({hi:.3g})",
                "Park et al. (2012) nine-type index from the joint sign of ΔV_avg and "
                "Δl_sol: 1–3 amplification, 4–6 knockdown, 7–9 neither",
                "amplify (types 1–3) / knockdown (types 4–6); 'robust' is CMM's own flag "
                "(the FVA minimum |flux| also rises monotonically), not Park's criterion",
            ],
        )
        self.fseof_table.setRowCount(len(rows))
        for i, (rid, mlo, mhi, park, tag, _robust, _delta) in enumerate(rows):
            self.fseof_table.setItem(i, 0, QTableWidgetItem(rid))
            self.fseof_table.setItem(i, 1, QTableWidgetItem(f"{mlo:.4g}"))
            self.fseof_table.setItem(i, 2, QTableWidgetItem(f"{mhi:.4g}"))
            self.fseof_table.setItem(i, 3, QTableWidgetItem(park))
            self.fseof_table.setItem(i, 4, QTableWidgetItem(tag))

    def _fill_envelope_table(self, envelope) -> None:
        """List the production-envelope points: growth range attainable at each product flux."""

        self._production_result = (
            "envelope",
            envelope,
        )  # show-all toggle does not apply here
        frame = envelope.to_frame()
        self._set_production_table_headers(
            ["Product flux", "Growth (min)", "Growth (max)"],
            [
                "Enforced product secretion rate along the envelope",
                "Minimum growth rate feasible at this product flux (lower boundary)",
                "Maximum growth rate feasible at this product flux (upper boundary)",
            ],
        )
        self.fseof_table.setRowCount(len(frame))
        for i, (_, record) in enumerate(frame.iterrows()):
            self.fseof_table.setItem(
                i, 0, QTableWidgetItem(f"{record['product_flux']:.4g}")
            )
            self.fseof_table.setItem(
                i, 1, QTableWidgetItem(f"{record['growth_min']:.4g}")
            )
            self.fseof_table.setItem(
                i, 2, QTableWidgetItem(f"{record['growth_max']:.4g}")
            )

    def _run_production(self, action) -> None:
        """Run a production analysis, surfacing any error instead of crashing the UI."""

        try:
            action()
        except Exception as exc:
            self.yield_label.setText(
                f"Production analysis failed: {html.escape(str(exc))}"
            )
            self.status_label.setText(f"Production analysis failed: {exc}")

    def run_theoretical_yield(self) -> None:
        def _do():
            self.fseof_table.setRowCount(0)  # yield is a scalar, not a target scan
            self._production_result = (
                None  # nothing for the show-all toggle to re-render
            )
            product = self._current_product()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = theoretical_yield(
                    self.model,
                    product,
                    substrate=self._current_substrate(),
                    condition=self._production_condition(),
                )
            ceiling = (
                f", carbon ceiling {result.carbon_ceiling:.2f}"
                if result.carbon_ceiling is not None
                else ""
            )
            # 0.4.0 reports how much of the product carbon CO2 uptake supplies, because
            # `co2_fixed` on its own let a 15.9%-inflated yield pass every guard by sitting
            # below the ceiling. Surface the fraction, and the imbalance if one is left over.
            co2 = (
                f"; {result.co2_carbon_fraction:.1%} of the product carbon comes from net "
                f"CO₂ uptake ({result.co2_uptake:.2f})"
                if result.co2_carbon_fraction > 0.0
                else ""
            )
            if result.carbon_imbalance:
                co2 += (
                    "; <span style='color:#b45309'>⚠ carbon imbalance: the excess is not "
                    "explained by CO₂ uptake</span>"
                )
            for warning in caught:
                self.status_label.setText(str(warning.message))
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
            condition = self._production_condition()
            label = self.anaerobic_combo.currentText()
            substrate = self._current_substrate()
            envelope = self._run_model_in_background(
                lambda worker_model: production_envelope(
                    worker_model,
                    product,
                    substrate=substrate,
                    condition=condition,
                    points=20,
                ),
                label="Computing production envelope…",
            )
            fig = production_envelope_figure(
                envelope, title=f"Production envelope — {product} ({label})"
            )
            self._set_production_figure(fig)
            self._fill_envelope_table(envelope)
            self.status_label.setText(f"Production envelope computed ({label}).")

        self._run_production(_do)

    def run_fseof_plot(self) -> None:
        def _do():
            product = self._current_product()
            condition = self._production_condition()
            result = self._run_model_in_background(
                lambda worker_model: fseof(
                    worker_model, product, n_steps=10, condition=condition
                ),
                label="Running FSEOF…",
            )
            if result.metadata.get("max_product", 0.0) <= 1e-9:
                self.yield_label.setText(
                    f"FSEOF not meaningful: theoretical yield of {html.escape(product)} is zero "
                    "— check the exchange is open and the product is reachable."
                )
                self.status_label.setText("FSEOF: zero theoretical yield.")
                return
            fig = fseof_figure(result, top_n=6, title=f"FSEOF targets — {product}")
            self._set_production_figure(fig)
            self._fill_fseof_table(result)
            # Rank the listed targets the same way the figure does (by flux increase).
            levels = list(result.enforced_levels)
            ranked = sorted(
                result.amplification_targets(),
                key=lambda rid: (
                    abs(result.trends.loc[rid, levels[-1]])
                    - abs(result.trends.loc[rid, levels[0]])
                ),
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
            condition = self._production_condition()
            # n_steps is left at the library default (10), which is Park et al.'s stated
            # minimum; the GUI used to force 8, undercutting it.
            result = self._run_model_in_background(
                lambda worker_model: fvseof(worker_model, product, condition=condition),
                label="Running FVSEOF…",
            )
            if result.metadata.get("max_product", 0.0) <= 1e-9:
                self.yield_label.setText(
                    f"FVSEOF not meaningful: theoretical yield of {html.escape(product)} is zero."
                )
                self.status_label.setText("FVSEOF: zero theoretical yield.")
                return
            fig = fvseof_figure(
                result, top_n=5, title=f"FVSEOF amplification targets — {product}"
            )
            self._set_production_figure(fig)
            self._fill_fvseof_table(result)
            # Park's types 1-3, in his priority order (ascending l_sol), are the published
            # selection; `robust` is CMM's own construct and is labelled as such.
            amplify = result.amplification_targets()[:8]
            robust = result.robust_targets()[:8]
            self.yield_label.setText(
                "FVSEOF amplification targets (Park types 1–3, by ascending l_sol): "
                + html.escape(", ".join(amplify) if amplify else "none")
                + "<br>Also rising monotonically under enforced flux "
                "(CMM's own “robust” flag, not Park et al.'s criterion): "
                + html.escape(", ".join(robust) if robust else "none")
            )
            self.status_label.setText("FVSEOF complete.")

        self._run_production(_do)

    # -- omics integration --------------------------------------------------

    def _set_omics_source(self, table, label: str) -> None:
        """Store a loaded gene × condition table and list its conditions (does not compute)."""

        self._omics_table_df = table
        # Discard any previous result so the show-all toggle can't re-render stale data, and
        # so a new source's methods must be recomputed before they back a Comparison template.
        self._omics_fluxes_by_condition = None
        self.omics_map_btn.setEnabled(False)
        self.omics_map_combo.clear()
        if hasattr(self, "omics_map_label"):
            self._refresh_omics_map_label()
        self._omics_conditions_order = []
        self._omics_expression = None
        self._omics_computed_methods = set()
        self._refresh_comparison_templates()
        self.omics_table.setRowCount(0)
        conditions = [str(c) for c in table.columns]
        self.omics_file_label.setText(
            f"{label} ({table.shape[0]} genes × {len(conditions)} condition(s))"
        )
        self.omics_cond_list.blockSignals(True)
        self.omics_cond_list.clear()
        for cond in conditions:
            item = QListWidgetItem(cond)
            item.setFlags(
                (item.flags() | Qt.ItemIsUserCheckable) & ~Qt.ItemIsSelectable
            )
            item.setCheckState(Qt.Checked)  # compute all conditions by default
            self.omics_cond_list.addItem(item)
        self.omics_cond_list.blockSignals(False)
        self.omics_compute_btn.setEnabled(True)
        self.omics_summary.setText(
            f"Loaded <b>{html.escape(label)}</b>: {len(conditions)} condition(s) "
            f"({html.escape(', '.join(conditions))}). Check conditions and click Compute."
        )
        self.status_label.setText(f"Loaded expression source: {label}.")

    def _checked_omics_conditions(self) -> list[str]:
        return [
            self.omics_cond_list.item(i).text()
            for i in range(self.omics_cond_list.count())
            if self.omics_cond_list.item(i).checkState() == Qt.Checked
        ]

    def load_omics_demo(self) -> None:
        """Load a deterministic single-condition demo table (does not compute)."""

        import numpy as np
        import pandas as pd

        genes = [g.id for g in self.model.genes]
        if not genes:
            self.omics_summary.setText(
                "Model has no genes — omics integration is not applicable."
            )
            return
        rng = np.random.default_rng(0)
        values = {g: float(rng.uniform(1.0, 100.0)) for g in genes}
        table = pd.DataFrame({"demo": pd.Series(values)})
        self._set_omics_source(table, "demo expression")

    def load_expression_csv(self) -> None:
        """Load a gene × condition expression CSV/TSV as the source (compute is separate)."""

        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load gene-expression CSV",
            "",
            "CSV/TSV files (*.csv *.tsv);;All files (*)",
        )
        if not path:
            return
        try:
            table = read_expression_table(path)
        except Exception as exc:
            self.omics_summary.setText(
                f"Could not read expression CSV: {html.escape(str(exc))}"
            )
            return
        import os

        self._set_omics_source(table, os.path.basename(path))

    def run_omics_demo(self) -> None:
        """Load the demo expression and compute immediately (used by scenarios/tests)."""

        self.load_omics_demo()
        if self._omics_table_df is not None:
            self.compute_omics()

    def compute_omics(self) -> None:
        """Predict a flux column for each checked condition (E-Flux2 or LAD)."""

        import pandas as pd

        if self._omics_table_df is None:
            self.omics_summary.setText("Load an expression source (CSV or demo) first.")
            return
        conditions = self._checked_omics_conditions()
        if not conditions:
            self.omics_summary.setText("Check at least one condition to compute.")
            return
        method = self.omics_method_combo.currentText()
        table = self._omics_table_df

        def _compute(worker_model):
            return predict_condition_fluxes(
                worker_model, table, method=method, conditions=conditions
            )

        try:
            predictions = self._run_model_in_background(
                _compute,
                label=f"Computing {method} for {len(conditions)} condition(s)…",
            )
        except Exception as exc:
            self.omics_summary.setText(
                f"Omics integration failed: {html.escape(str(exc))}"
            )
            self.status_label.setText("Omics integration failed.")
            return

        # Keep only conditions that produced an optimal flux state; note the rest.
        per_condition: dict[str, dict[str, float]] = {}
        failed: list[str] = []
        for cond in conditions:
            try:
                per_condition[cond] = predictions.fluxes(cond)
            except Exception:
                failed.append(cond)
        ok = [c for c in conditions if c in per_condition]
        if not ok:
            self._omics_fluxes_by_condition = None
            self.omics_map_btn.setEnabled(False)
            self.omics_map_combo.clear()
            self._omics_conditions_order = []
            self.omics_table.setColumnCount(1)
            self.omics_table.setHorizontalHeaderLabels(["Reaction"])
            self.omics_table.setRowCount(0)
            self.omics_summary.setText(
                f"{method.upper()} produced no optimal flux state for the selected conditions."
            )
            return

        # Cache the result so the show-all toggle can re-render without recomputing.
        self._omics_fluxes_by_condition = per_condition
        self._omics_conditions_order = ok
        self._omics_method_label = method.upper()
        self._loading = True
        self.omics_map_combo.clear()
        self.omics_map_combo.addItems(ok)
        self._loading = False
        self._refresh_omics_map_label()
        self.omics_map_btn.setEnabled(True)
        self.omics_map_btn.setToolTip(
            f"Draw the {method.upper()} flux state of the chosen condition on the flux map."
        )
        n_active = len(
            {
                rid
                for cond in ok
                for rid, v in per_condition[cond].items()
                if abs(v) > 1e-6
            }
        )
        self._render_omics_table()

        # Reuse the first computed condition as the Comparison tab's LAD/E-Flux2 reference,
        # and record this method so only it is offered as a reference template.
        first = ok[0]
        self._omics_expression = {
            str(gene): float(value)
            for gene, value in table[first].items()
            if pd.notna(value)
        }
        self._omics_computed_methods.add(method)
        self._refresh_comparison_templates()  # this method is now available as a reference
        note = (
            f" ({len(failed)} infeasible: {html.escape(', '.join(failed))})"
            if failed
            else ""
        )
        self.omics_summary.setText(
            f"<b>{method.upper()}</b>: {len(ok)} condition(s) × {n_active} active "
            f"reactions{note}. Reference for Comparison = '{html.escape(first)}'."
        )
        self.status_label.setText(f"Omics integration complete ({method}).")

    def _render_omics_table(self) -> None:
        """Render the cached per-condition fluxes, honouring the active-only / show-all toggle."""

        per = self._omics_fluxes_by_condition
        ok = self._omics_conditions_order
        if not per or not ok:
            return
        if self.omics_show_all_check.isChecked():
            rxn_set = {rid for cond in ok for rid in per[cond]}
        else:
            rxn_set = {
                rid for cond in ok for rid, v in per[cond].items() if abs(v) > 1e-6
            }
        rxns = sorted(
            rxn_set, key=lambda rid: -max(abs(per[c].get(rid, 0.0)) for c in ok)
        )
        self.omics_table.setColumnCount(1 + len(ok))
        self.omics_table.setHorizontalHeaderLabels(["Reaction", *ok])
        self.omics_table.setRowCount(len(rxns))
        for i, rid in enumerate(rxns):
            self.omics_table.setItem(i, 0, QTableWidgetItem(rid))
            for j, cond in enumerate(ok, start=1):
                self.omics_table.setItem(
                    i, j, QTableWidgetItem(f"{per[cond].get(rid, 0.0):.4g}")
                )

    def _update_revert_run_state(self) -> None:
        ready = (
            self._revert_source_expression is not None
            and self._revert_target_expression is not None
        )
        self.revert_run_btn.setEnabled(ready)
        self._update_revert_direction_controls()

    def _revert_replicate_counts(self) -> tuple[int, int]:
        """Replicate columns currently loaded on each side; 0 where no file is loaded.

        A caller may assign a plain gene -> value mapping instead of a read frame, which is
        one measurement per gene however it arrived.
        """

        counts = [
            0 if state is None else int(getattr(state, "shape", (0, 1))[1])
            for state in (
                self._revert_source_expression,
                self._revert_target_expression,
            )
        ]
        return counts[0], counts[1]

    def _on_revert_significance_changed(self) -> None:
        """A user-driven pick of the test overrides the tab's own fallback for good."""

        if not self._revert_significance_syncing:
            self._revert_significance_forced = False
        self._update_revert_direction_controls()

    def _update_revert_direction_controls(self) -> None:
        """Offer the t-test only when the loaded files can carry it, and say why not.

        A t-test needs at least two replicates on each side. Leaving the option selectable on
        a one-column file would trade a decision the user can make before running for a solver
        error after it, so the option is disabled and the reason is stated next to it.
        """

        combo = self.revert_significance_combo
        n_source, n_target = self._revert_replicate_counts()
        loaded = n_source and n_target
        can_ttest = n_source >= 2 and n_target >= 2
        ttest_index = combo.findData("ttest")
        item = combo.model().item(ttest_index)
        if item is not None:
            item.setEnabled(bool(can_ttest))
        # Dropping to fold change on a file that cannot be tested is the tab's decision, not
        # the user's, so it is remembered and undone once replicates arrive. A fold-change run
        # the user chose while the t-test was available is left alone.
        self._revert_significance_syncing = True
        try:
            if not can_ttest and combo.currentData() == "ttest":
                combo.setCurrentIndex(combo.findData("fold_change"))
                self._revert_significance_forced = True
            elif can_ttest and self._revert_significance_forced:
                combo.setCurrentIndex(ttest_index)
                self._revert_significance_forced = False
        finally:
            self._revert_significance_syncing = False
        significance = combo.currentData()
        self.revert_p_value_spin.setEnabled(significance == "ttest")
        self.revert_fold_change_spin.setEnabled(significance != "ttest")

        if not loaded:
            note = (
                "Load both files to see how many replicates they carry. The t-test needs at "
                "least two columns on each side."
            )
        elif can_ttest:
            note = (
                f"{n_source} source and {n_target} target replicates — Yizhak et al.'s "
                f"t-test is available."
            )
        else:
            note = (
                f"{n_source} source and {n_target} target replicate(s): too few for a "
                f"t-test, so the direction map is cut on fold change. State that the "
                f"published test was not applied when reporting this run."
            )
        self.revert_replicate_note.setText(note)

    def load_revert_expression(self, role: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Load {role} expression",
            "",
            "CSV/TSV files (*.csv *.tsv);;All files (*)",
        )
        if not path:
            return
        try:
            expression = _read_expression_matrix(path)
        except Exception as exc:
            self.revert_summary.setText(
                f"Could not read {html.escape(role)} expression: {html.escape(str(exc))}"
            )
            return
        import os

        name = os.path.basename(path)
        n_genes, n_replicates = expression.shape
        plural = "" if n_replicates == 1 else "s"
        caption = f"{name} ({n_genes} genes x {n_replicates} column{plural})"
        if role == "source":
            self._revert_source_expression = expression
            self.revert_source_label.setText(caption)
        else:
            self._revert_target_expression = expression
            self.revert_target_label.setText(caption)
        self._update_revert_run_state()
        self.status_label.setText(
            f"Loaded {role} expression ({n_genes} genes, {n_replicates} column{plural})."
        )

    def run_loaded_revert(self) -> None:
        if (
            self._revert_source_expression is None
            or self._revert_target_expression is None
        ):
            self.revert_summary.setText(
                "Load both source and target expression files first."
            )
            self._update_revert_run_state()
            return
        self.run_revert(self._revert_source_expression, self._revert_target_expression)

    def run_revert(
        self,
        source_expression: "Mapping[str, float] | pd.DataFrame",
        target_expression: "Mapping[str, float] | pd.DataFrame",
    ) -> None:
        method = self.method_combo.currentText()
        perturbation = self.perturbation_combo.currentText()
        alpha = self.alpha_spin.value()
        epsilon = self.revert_epsilon_spin.value()
        top_changed = self.revert_top_changed_spin.value() or None
        reference_method = self.revert_reference_combo.currentText()
        fold_change = self.revert_fold_change_spin.value()
        p_value_cutoff = self.revert_p_value_spin.value()
        source_matrix = _as_expression_matrix(source_expression)
        target_matrix = _as_expression_matrix(target_expression)
        # A caller may drive this directly with one value per gene, so what is actually
        # runnable is decided from the matrices in hand rather than from the combo's state.
        significance = self.revert_significance_combo.currentData() or "fold_change"
        if significance == "ttest" and (
            source_matrix.shape[1] < 2 or target_matrix.shape[1] < 2
        ):
            significance = "fold_change"

        # E-Flux2 and LAD scale reaction bounds by a linear expression level, so the flux
        # state gets the replicate mean as measured; only the comparison moves to log space.
        source_linear = {
            str(gene): float(value)
            for gene, value in source_matrix.mean(axis=1).items()
        }
        source_log2 = _expression_in_log2(source_matrix)
        target_log2 = _expression_in_log2(target_matrix)

        def _compute(worker_model):
            # LAD takes no objective fraction: it fits |flux| to expression-derived targets
            # rather than holding an optimum, so passing one would be meaningless.
            extra = {"objective_fraction": 1.0} if reference_method == "eflux2" else {}
            source_result = integrate_expression(
                worker_model,
                source_linear,
                method=reference_method,
                **extra,
            )
            if source_result.status != "optimal" or not source_result.fluxes:
                raise ValueError(
                    f"source expression could not produce a valid source-state flux with "
                    f"{reference_method} ({source_result.status})"
                )
            reference = source_result.to_flux_state(
                f"source_expression_state_{reference_method}"
            )
            if significance == "ttest":
                genes = gene_directions_from_replicates(
                    source_log2, target_log2, p_value_cutoff=p_value_cutoff
                )
            else:
                genes = gene_directions_by_fold_change(
                    source_log2,
                    target_log2,
                    up_threshold=fold_change,
                    down_threshold=fold_change,
                )
            gene_dirs = {
                str(gene): int(value) for gene, value in genes["direction"].items()
            }
            direction = reaction_directions(
                worker_model, gene_dirs, reference=reference
            )
            if top_changed is not None:
                # With a t-test run, the smallest gene P value is the closer analogue of the
                # paper's gene-level test than |log2FC|; without one there are no P values.
                direction = restrict_to_top_changed(
                    worker_model,
                    direction,
                    genes["log2_fold_change"].to_dict(),
                    top_changed,
                    gene_p_values=(
                        genes["p_value"].to_dict() if significance == "ttest" else None
                    ),
                )
            if not direction.nonsteady():
                raise ValueError(
                    "no reaction was labelled as changed; loosen the significance threshold "
                    "or check that the two expression states differ"
                )
            ranking = revert_targets(
                worker_model,
                None,
                reference,
                direction,
                method=method,
                alpha=alpha,
                epsilon=epsilon,
                perturbation=perturbation,
            )
            return ranking, genes, direction

        try:
            ranking, genes, direction = self._run_model_in_background(
                _compute, label=f"Running revert ({method})…"
            )
        except (
            Exception
        ) as exc:  # surface solver-capability/model errors instead of freezing
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
                    self.revert_table.item(i, col).setBackground(QColor("#dce9d6"))
        # The direction map is the run's largest hidden assumption — the same model and
        # method rank differently under a t-test and a fold-change cut — so how it was built
        # is stated with the result rather than left in the controls.
        if significance == "ttest":
            test = (
                f"Student's t-test at P &lt; {p_value_cutoff:g} over "
                f"{source_matrix.shape[1]}v{target_matrix.shape[1]} replicates"
            )
        else:
            test = f"|log2 fold change| &ge; {fold_change:g}"
        n_significant = int(genes["significant"].sum())
        evidence = (
            f"Direction map: {test} &mdash; "
            f"{n_significant} of {len(genes)} shared genes changed, "
            f"{len(direction.nonsteady())} reactions labelled."
        )
        best = ranking.best()
        if best is not None:
            self.revert_summary.setText(
                f"<b>Top normalization target: {best.target_id}</b> "
                f"(score {best.score:.4g}, method {method}, {perturbation} knockout). "
                f"Knocking it out moves the source state most toward the target state.<br>"
                f"{evidence}"
            )
        else:
            self.revert_summary.setText(f"No target was scored. {evidence}")
        self.status_label.setText(f"Revert prediction complete ({method}).")

    # -- transformation finder (A -> B) -------------------------------------

    def _build_transformation_tab(self) -> QWidget:
        """Rank knockouts that move a source flux state toward a target flux state (A→B)."""

        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QGroupBox("Transformation targets (move source → target flux state)")
        form = QFormLayout(controls)
        self.transform_method_combo = QComboBox()
        self.transform_method_combo.addItems(["moma", "mta"])
        self.transform_method_combo.setToolTip(
            "moma: rank by how much the minimal-adjustment state moves toward the target. "
            "mta: published single MTA MIQP on the source→target direction."
        )
        self.transform_omics_combo = QComboBox()
        self.transform_omics_combo.addItems(["eflux2", "lad"])
        self.transform_omics_combo.setToolTip(
            "How each expression vector is turned into a flux state"
        )
        self.transform_perturbation_combo = QComboBox()
        self.transform_perturbation_combo.addItems(["gene", "reaction"])
        self.transform_alpha_spin = QDoubleSpinBox()
        self.transform_alpha_spin.setRange(0.0, 1.0)
        self.transform_alpha_spin.setSingleStep(0.01)
        self.transform_alpha_spin.setValue(0.66)
        form.addRow("Method:", self.transform_method_combo)
        form.addRow("Predict states with:", self.transform_omics_combo)
        form.addRow("Knockout level:", self.transform_perturbation_combo)
        form.addRow("Transformation weight α (mta):", self.transform_alpha_spin)

        source_row = QHBoxLayout()
        source_btn = QPushButton("Load source (A) CSV/TSV…")
        source_btn.clicked.connect(lambda: self.load_transform_expression("source"))
        self.transform_source_label = QLabel("not loaded")
        source_row.addWidget(source_btn)
        source_row.addWidget(self.transform_source_label, 1)
        form.addRow("Source expression:", source_row)

        target_row = QHBoxLayout()
        target_btn = QPushButton("Load target (B) CSV/TSV…")
        target_btn.clicked.connect(lambda: self.load_transform_expression("target"))
        self.transform_target_label = QLabel("not loaded")
        target_row.addWidget(target_btn)
        target_row.addWidget(self.transform_target_label, 1)
        form.addRow("Target expression:", target_row)

        self.transform_run_btn = QPushButton("Run transformation")
        self.transform_run_btn.setEnabled(False)
        self.transform_run_btn.clicked.connect(self.run_transformation)
        self.transform_run_btn.setToolTip(
            "Load both source and target expression files to enable."
        )
        form.addRow("", self.transform_run_btn)
        layout.addWidget(controls)

        self.transform_summary = QLabel(
            "Load source (A) and target (B) expression files; each is predicted to a flux "
            "state, then knockouts are ranked by how well they move A toward B."
        )
        self.transform_summary.setWordWrap(True)
        layout.addWidget(self.transform_summary)

        self.transform_table = QTableWidget(0, 3)
        self.transform_table.setHorizontalHeaderLabels(["Rank", "Target", "Score"])
        self.transform_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.Stretch
        )
        self.transform_table.verticalHeader().setVisible(False)
        self.transform_table.setAlternatingRowColors(True)
        layout.addWidget(self.transform_table, 1)
        return tab

    def _update_transform_run_state(self) -> None:
        ready = (
            self._transform_source_expression is not None
            and self._transform_target_expression is not None
        )
        self.transform_run_btn.setEnabled(ready)

    def load_transform_expression(self, role: str) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Load {role} expression",
            "",
            "CSV/TSV files (*.csv *.tsv);;All files (*)",
        )
        if not path:
            return
        try:
            expression = _read_expression_vector(path)
        except Exception as exc:
            self.transform_summary.setText(
                f"Could not read {html.escape(role)} expression: {html.escape(str(exc))}"
            )
            return
        import os

        name = os.path.basename(path)
        if role == "source":
            self._transform_source_expression = expression
            self.transform_source_label.setText(f"{name} ({len(expression)} genes)")
        else:
            self._transform_target_expression = expression
            self.transform_target_label.setText(f"{name} ({len(expression)} genes)")
        self._update_transform_run_state()
        self.status_label.setText(
            f"Loaded {role} (A→B) expression ({len(expression)} genes)."
        )

    def run_transformation(self) -> None:
        if (
            self._transform_source_expression is None
            or self._transform_target_expression is None
        ):
            self.transform_summary.setText(
                "Load both source (A) and target (B) expression files first."
            )
            return
        method = self.transform_method_combo.currentText()
        omics_method = self.transform_omics_combo.currentText()
        perturbation = self.transform_perturbation_combo.currentText()
        alpha = self.transform_alpha_spin.value()
        source_expr = self._transform_source_expression
        target_expr = self._transform_target_expression

        def _compute(worker_model):
            source_state = integrate_expression(
                worker_model, source_expr, method=omics_method
            ).to_flux_state("source")
            target_state = integrate_expression(
                worker_model, target_expr, method=omics_method
            ).to_flux_state("target")
            return transformation_targets(
                worker_model,
                source_state,
                target_state,
                method=method,
                perturbation=perturbation,
                alpha=alpha,
            )

        try:
            ranking = self._run_model_in_background(
                _compute, label=f"Running transformation ({method})…"
            )
        except (
            Exception
        ) as exc:  # surface solver-capability/model errors instead of freezing
            self.transform_summary.setText(
                f"Transformation failed: {html.escape(str(exc))}"
            )
            self.status_label.setText(f"Transformation failed ({method}).")
            return
        rows = ranking.to_records()
        self.transform_table.setRowCount(len(rows))
        for i, record in enumerate(rows):
            self.transform_table.setItem(i, 0, QTableWidgetItem(str(record["rank"])))
            self.transform_table.setItem(i, 1, QTableWidgetItem(record["target_id"]))
            self.transform_table.setItem(
                i, 2, QTableWidgetItem(f"{record['score']:.4g}")
            )
            if i == 0:
                for col in range(3):
                    self.transform_table.item(i, col).setBackground(QColor("#dce9d6"))
        best = ranking.best()
        if best is not None:
            self.transform_summary.setText(
                f"<b>Top A→B target: {html.escape(best.target_id)}</b> "
                f"(score {best.score:.4g}, method {method}, {omics_method} states, "
                f"{perturbation} knockout). Higher score = moves the source state more toward "
                "the target."
            )
        else:
            self.transform_summary.setText("No candidate knockouts scored.")
        self.status_label.setText(f"Transformation complete ({method}).")
