"""Headless smoke test for the Qt shell. Runs the real widgets under the offscreen platform."""

from __future__ import annotations

import json
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402

pytest.importorskip("qtpy")

from cmm.app.main_window import (  # noqa: E402
    CmmMainWindow,
    _expression_in_log2,
    _read_expression_matrix,
    _read_expression_vector,
)
from cmm.features.comparison import ROOM_TOLERANCES, reference_flux  # noqa: E402
from cmm.app.screenshots import (  # noqa: E402
    SOURCE_EXPRESSION,
    TARGET_EXPRESSION,
    build_demo_model,
)
from qtpy.QtCore import Qt  # noqa: E402
from qtpy.QtWidgets import QApplication, QLabel  # noqa: E402


@pytest.fixture(scope="module")
def app():
    try:
        instance = QApplication.instance() or QApplication([])
    except Exception as exc:  # pragma: no cover - no usable Qt platform
        pytest.skip(f"no Qt platform available: {exc}")
    return instance


def _clear_ko(window):
    """Empty the Comparison tab's chosen-knockout set."""

    window._clear_ko_selected()


def _select_ko(window, ids):
    """Set the Comparison tab's chosen-knockout set to exactly the given ids."""

    window._clear_ko_selected()
    window._add_ko_targets(list(ids))


@pytest.mark.requires_miqp
def test_window_builds_and_drives_services(app):
    window = CmmMainWindow(build_demo_model())

    window.run_fba()
    assert "optimal" in window.objective_label.text()
    assert window.reaction_table.rowCount() == 6

    window.run_fva()
    assert window.sim_table.rowCount() == 6

    window.run_revert(SOURCE_EXPRESSION, TARGET_EXPRESSION)
    # The disease-branch gene must be the top-ranked normalization target in the GUI table.
    assert window.revert_table.item(0, 1).text() == "g2"
    assert "g2" in window.revert_summary.text()


def test_bound_editing_updates_model_and_fba(app):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    assert "10" in window.objective_label.text()

    # Programmatic bound edit reflects in the model and the table.
    window.set_reaction_bounds("SUP_A", upper=4.0)
    assert window.model.reactions.SUP_A.upper_bound == 4.0
    window.run_fba()
    assert "4" in window.objective_label.text()  # growth now substrate-limited to 4

    # In-table edit path: changing a cell triggers _on_bound_edited.
    sup_row = next(
        r
        for r in range(window.reaction_table.rowCount())
        if window.reaction_table.item(r, 0).text() == "SUP_A"
    )
    window.reaction_table.item(sup_row, 2).setText("7")
    assert window.model.reactions.SUP_A.upper_bound == 7.0
    window.run_fba()
    assert "7" in window.objective_label.text()

    # Invalid input is rejected and the cell reverts.
    window.reaction_table.item(sup_row, 2).setText("not a number")
    assert window.model.reactions.SUP_A.upper_bound == 7.0


def test_infeasible_fba_does_not_crash(app):
    window = CmmMainWindow(build_demo_model())
    # Force infeasibility: demand uptake the network cannot satisfy.
    window.model.reactions.SUP_A.bounds = (5.0, 5.0)
    window.model.reactions.BIOMASS.lower_bound = 1000.0
    with pytest.warns(UserWarning, match="infeasible"):
        window.run_fba()  # must not raise
    assert (
        "infeasible" in window.objective_label.text()
        or "error" in window.objective_label.text()
    )


def test_crossing_bound_edit_does_not_crash(app):
    window = CmmMainWindow(build_demo_model())
    # Setting lower above the current upper must not raise (atomic assignment + clamp).
    window.set_reaction_bounds("SUP_A", lower=2000.0)
    assert (
        window.model.reactions.SUP_A.lower_bound
        <= window.model.reactions.SUP_A.upper_bound
    )
    # In-table crossing edit clamps and reflects the stored value in the cell.
    sup_row = next(
        r
        for r in range(window.reaction_table.rowCount())
        if window.reaction_table.item(r, 0).text() == "SUP_A"
    )
    window.reaction_table.item(sup_row, 1).setText("99999")
    lo = window.model.reactions.SUP_A.lower_bound
    assert lo <= window.model.reactions.SUP_A.upper_bound
    assert window.reaction_table.item(sup_row, 1).text() == f"{lo:g}"


def test_production_buttons_disabled_without_exchanges(app):
    from cobra import Metabolite, Model, Reaction

    model = Model("no_exchange")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    r = Reaction("R")  # internal (2 metabolites) -> not a boundary/exchange reaction
    r.add_metabolites({a: -1, b: 1})
    r.gene_reaction_rule = "g"
    model.add_reactions([r])
    model.objective = r
    assert len(model.exchanges) == 0
    window = CmmMainWindow(model)
    assert all(not btn.isEnabled() for btn in window._production_buttons)
    # Even if invoked directly, the guarded runner must not raise.
    window.run_theoretical_yield()
    assert (
        "failed" in window.yield_label.text()
        or "unavailable" in window.yield_label.text()
    )


def test_menu_bar_has_expected_menus(app):
    window = CmmMainWindow(build_demo_model())
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    for expected in ("Analysis", "Model", "Config"):
        assert expected in titles


def test_window_uses_constraint_based_branding(app):
    window = CmmMainWindow(build_demo_model())
    expected = "CMM — Constraint-based Metabolic Modeling"
    assert window.windowTitle() == expected
    title = window.findChild(QLabel, "title")
    assert title is not None
    assert title.text() == expected


@pytest.mark.requires_qp
def test_omics_integration_in_gui(app):
    window = CmmMainWindow(build_demo_model())
    window.omics_method_combo.setCurrentText("eflux2")
    window.run_omics_demo()
    assert window.omics_table.rowCount() > 0
    assert "EFLUX2" in window.omics_summary.text()

    window.omics_method_combo.setCurrentText("lad")
    window.run_omics_demo()
    assert "LAD" in window.omics_summary.text()


def test_expression_vector_reads_tsv(tmp_path):
    path = tmp_path / "expression.tsv"
    path.write_text("gene\texpression\ng1\t2.5\ng2\t7\n")

    assert _read_expression_vector(str(path)) == {"g1": 2.5, "g2": 7.0}


def test_expression_vector_rejects_duplicates_and_negative_values(tmp_path):
    duplicate = tmp_path / "duplicate.csv"
    duplicate.write_text("gene,expression\ng1,1\ng1,2\n")
    with pytest.raises(ValueError, match="duplicate"):
        _read_expression_vector(str(duplicate))

    negative = tmp_path / "negative.csv"
    negative.write_text("gene,expression\ng1,-1\n")
    with pytest.raises(ValueError, match="non-negative"):
        _read_expression_vector(str(negative))


def _write_replicates(path, expression, offsets):
    """One CSV per state: a gene column, then one column per replicate.

    ``offsets`` shifts each replicate off the shared mean so the t-test sees within-group
    variance; without any the variance is zero and the P value is undefined.
    """

    header = "gene," + ",".join(f"rep{i + 1}" for i in range(len(offsets)))
    rows = [
        ",".join([gene] + [f"{value * (1 + off):.6g}" for off in offsets])
        for gene, value in expression.items()
    ]
    path.write_text("\n".join([header, *rows]) + "\n")
    return path


def test_expression_matrix_reads_replicate_columns(tmp_path):
    path = tmp_path / "replicates.tsv"
    path.write_text("gene\tr1\tr2\tr3\ng1\t1\t2\t3\ng2\t4\t5\t6\n")

    frame = _read_expression_matrix(str(path))
    assert list(frame.index) == ["g1", "g2"]
    assert frame.shape == (2, 3)


def test_expression_matrix_rejects_duplicate_genes_and_missing_values(tmp_path):
    duplicate = tmp_path / "duplicate.csv"
    duplicate.write_text("gene,r1,r2\ng1,1,2\ng1,3,4\n")
    with pytest.raises(ValueError, match="duplicate"):
        _read_expression_matrix(str(duplicate))

    missing = tmp_path / "missing.csv"
    missing.write_text("gene,r1,r2\ng1,1,\n")
    with pytest.raises(ValueError, match="missing"):
        _read_expression_matrix(str(missing))

    text_only = tmp_path / "text.csv"
    text_only.write_text("gene,label\ng1,up\n")
    with pytest.raises(ValueError, match="numeric"):
        _read_expression_matrix(str(text_only))

    # Expression is linear throughout the app, so a log-transformed table is refused at the
    # file dialog rather than surfacing later as an E-Flux2 solver error.
    negative = tmp_path / "log2.csv"
    negative.write_text("gene,r1,r2\ng1,-1.5,-0.5\n")
    with pytest.raises(ValueError, match="non-negative"):
        _read_expression_matrix(str(negative))


def test_the_comparison_moves_to_log_space_per_replicate():
    import pandas as pd

    linear = pd.DataFrame({"r1": [1.0, 7.0], "r2": [3.0, 15.0]}, index=["g1", "g2"])

    # log2(x + 1) per replicate, not after averaging, so a one-column file reproduces the
    # log-ratio the tab computed before it accepted replicates.
    assert _expression_in_log2(linear)["r1"].tolist() == [1.0, 3.0]
    assert _expression_in_log2(linear)["r2"].tolist() == [2.0, 4.0]


def test_stage_forms_do_not_inherit_the_platform_form_layout_style(app):
    """The stage boxes must lay out identically under macOS's QFormLayout style hints.

    QFormLayout reads its alignment and growth policy from the style. Most styles give
    AlignLeft + AllNonFixedFieldsGrow; macOS gives AlignHCenter + FieldsStayAtSizeHint, and
    under that pair every box centres its own label-plus-field block at that block's natural
    width - so the stages land at different x positions and no two inputs share a width. CI
    renders under a non-macOS style, so the hints are applied here directly rather than by
    swapping the application style, which would crash (QProxyStyle takes ownership).
    """

    from qtpy.QtWidgets import QAbstractSpinBox, QComboBox, QFormLayout, QGroupBox

    def placements(window):
        page = window.tabs.widget(window._tab_index("Revert Metabolism"))
        return {
            (control.mapTo(window, control.rect().topLeft()).x(), control.width())
            for box in page.findChildren(QGroupBox)
            for control in box.findChildren((QComboBox, QAbstractSpinBox))
        }

    window = CmmMainWindow(build_demo_model())
    window.resize(1400, 900)
    window._goto_tab(
        "Revert Metabolism"
    )  # a hidden page has no real geometry to measure
    for form in window._revert_stage_forms:
        form.setFormAlignment(Qt.AlignHCenter | Qt.AlignTop)
        form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)
    window.show()
    for _ in range(5):
        app.processEvents()

    # The test is only worth anything if these hints really do break the layout.
    assert len(placements(window)) > 1, (
        "macOS's hints no longer scatter the stage inputs"
    )

    window._pin_form_layout_policies()
    for _ in range(5):
        app.processEvents()
    fixed = placements(window)
    assert len(fixed) == 1, f"stage inputs disagree on x/width: {sorted(fixed)}"


def test_stage_inputs_do_not_stretch_across_a_wide_panel(app):
    """Inputs fill their column so the stages align, but only up to a readable width.

    They are told to grow because that is what puts every stage's inputs at one x; left
    uncapped, a spin box holding "0.66" runs the full width of the left panel on a large
    display. The wrapped note under the fields is not an input and must still span the box.
    """

    from qtpy.QtWidgets import QAbstractSpinBox, QComboBox, QGroupBox

    window = CmmMainWindow(build_demo_model())
    window.resize(2000, 1000)  # a display wide enough for the cap to bite
    window._goto_tab("Revert Metabolism")
    window.show()
    for _ in range(5):
        app.processEvents()

    page = window.tabs.widget(window._tab_index("Revert Metabolism"))
    controls = [
        control
        for box in page.findChildren(QGroupBox)
        for control in box.findChildren((QComboBox, QAbstractSpinBox))
    ]
    assert controls
    assert {control.width() for control in controls} == {260}

    note = window.revert_replicate_note
    assert note.width() > 400, (
        "the note is an explanation, not a field; it spans the box"
    )
    # A wrapped label reports one line's height unless its policy says otherwise, and is then
    # clipped mid-sentence. Check the text the tab actually shows, at its actual width.
    for text in (
        note.text(),
        "1 source and 1 target replicate(s): too few for a t-test, so the direction map is "
        "cut on fold change. State that the published test was not applied when reporting "
        "this run.",
    ):
        note.setText(text)
        for _ in range(5):
            app.processEvents()
        assert note.height() >= note.heightForWidth(note.width()), text


def test_revert_stage_boxes_share_one_label_column(app):
    """Four stacked stage boxes must put their inputs at one x, not four.

    Each QFormLayout otherwise sizes its label column to its own longest label, which strands
    a short label like 'Method:' against the box edge with its field far to the right.
    """

    from qtpy.QtGui import QFontMetrics
    from qtpy.QtWidgets import QFormLayout

    window = CmmMainWindow(build_demo_model())
    window.show()  # the column is squared on first show, once the fonts are resolved

    labels = []
    for form in window._revert_stage_forms:
        for row in range(form.rowCount()):
            item = form.itemAt(row, QFormLayout.LabelRole)
            if item is not None and item.widget() is not None:
                labels.append(item.widget())
    widths = {label.width() for label in labels}
    assert len(widths) == 1, (
        f"stage labels span several column widths: {sorted(widths)}"
    )
    assert {int(label.alignment()) & int(Qt.AlignRight) for label in labels} == {
        int(Qt.AlignRight)
    }
    # Pinned, not floored: a label allowed to grow past the column would widen its own stage
    # box alone and pull that stage's inputs out of line with every other stage.
    for label in labels:
        assert label.minimumWidth() == label.maximumWidth()
        assert (
            QFontMetrics(label.font()).horizontalAdvance(label.text()) <= label.width()
        )
    # Bound unused space against the actual styled text. A fixed 160px ceiling assumes the
    # host's font metrics instead of distinguishing wider glyphs from excess column padding.
    # Measuring before fonts settle would still fail by leaving an oversized label column.
    widest_text = max(
        QFontMetrics(label.font()).horizontalAdvance(label.text()) for label in labels
    )
    shared_width = widths.pop()
    assert 0 < widest_text <= shared_width <= widest_text + 8, (
        f"shared column {shared_width}px, widest styled text {widest_text}px"
    )


def test_stage_controls_share_one_height(app):
    """Combos and spin boxes in one form column must be the same height.

    They are styled separately, so their paddings drifted apart (32px against 30px) and each
    stage box ended up on its own vertical grid — the row pitch alternated with whichever
    control the row happened to hold.
    """

    from qtpy.QtWidgets import QAbstractSpinBox, QComboBox

    window = CmmMainWindow(build_demo_model())
    page = window.tabs.widget(window._tab_index("Revert Metabolism"))
    heights = {
        control.sizeHint().height()
        for kind in (QComboBox, QAbstractSpinBox)
        for control in page.findChildren(kind)
    }
    assert len(heights) == 1, f"stage controls disagree on height: {sorted(heights)}"
    assert heights.pop() <= 28


def test_revert_ttest_offered_only_with_replicates(app, tmp_path):
    window = CmmMainWindow(build_demo_model())
    combo = window.revert_significance_combo

    # Nothing loaded: the t-test entry is not selectable and the note says why.
    assert combo.model().item(combo.findData("ttest")).isEnabled() is False
    assert "two columns" in window.revert_replicate_note.text()

    single = _write_replicates(tmp_path / "one.csv", SOURCE_EXPRESSION, [0.0])
    window._revert_source_expression = _read_expression_matrix(str(single))
    window._revert_target_expression = _read_expression_matrix(str(single))
    window._update_revert_run_state()
    assert combo.model().item(combo.findData("ttest")).isEnabled() is False
    assert combo.currentData() == "fold_change"
    assert window.revert_fold_change_spin.isEnabled()
    assert not window.revert_p_value_spin.isEnabled()

    triplicate = _write_replicates(
        tmp_path / "three.csv", SOURCE_EXPRESSION, [-0.05, 0.0, 0.05]
    )
    window._revert_source_expression = _read_expression_matrix(str(triplicate))
    window._revert_target_expression = _read_expression_matrix(str(triplicate))
    window._update_revert_run_state()
    assert combo.model().item(combo.findData("ttest")).isEnabled() is True
    # The fold-change fallback was the tab's own, so replicates arriving restore the
    # published test without the user having to notice the combo.
    assert combo.currentData() == "ttest"
    assert window.revert_p_value_spin.isEnabled()
    assert not window.revert_fold_change_spin.isEnabled()
    assert "t-test is available" in window.revert_replicate_note.text()

    # A fold-change run the user chose while the t-test was available survives a reload.
    combo.setCurrentIndex(combo.findData("fold_change"))
    window._update_revert_run_state()
    assert combo.currentData() == "fold_change"


@pytest.mark.requires_miqp
def test_revert_tab_runs_replicates_through_the_ttest(app, tmp_path):
    window = CmmMainWindow(build_demo_model())
    source = _write_replicates(
        tmp_path / "source.csv", SOURCE_EXPRESSION, [-0.02, 0.0, 0.02]
    )
    target = _write_replicates(
        tmp_path / "target.csv", TARGET_EXPRESSION, [-0.02, 0.0, 0.02]
    )
    window._revert_source_expression = _read_expression_matrix(str(source))
    window._revert_target_expression = _read_expression_matrix(str(target))
    window._update_revert_run_state()
    combo = window.revert_significance_combo
    combo.setCurrentIndex(combo.findData("ttest"))

    window.run_loaded_revert()
    # Same disease-branch gene as the single-measurement run: the published t-test route must
    # not change which target wins on data that differ this strongly.
    assert window.revert_table.item(0, 1).text() == "g2"
    assert "t-test" in window.revert_summary.text()
    assert "3v3 replicates" in window.revert_summary.text()


@pytest.mark.requires_miqp
def test_revert_tab_runs_loaded_expressions(app):
    window = CmmMainWindow(build_demo_model())
    assert not window.revert_run_btn.isEnabled()

    window._revert_source_expression = dict(SOURCE_EXPRESSION)
    window._revert_target_expression = dict(TARGET_EXPRESSION)
    window._update_revert_run_state()
    assert window.revert_run_btn.isEnabled()

    window.run_loaded_revert()
    assert window.revert_table.item(0, 1).text() == "g2"
    assert "g2" in window.revert_summary.text()


def test_media_pfba_comparison_in_gui(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)

    # Apply a preset medium and confirm it changed the model's exchange bounds.
    window.medium_combo.setCurrentText("glucose_anaerobic")
    window.apply_selected_medium()
    assert window.model.reactions.EX_o2_e.lower_bound == 0.0

    # pFBA runs and reports growth.
    window.tabs.setCurrentIndex(0)
    window.run_pfba()
    assert "pFBA" in window.objective_label.text()

    # Comparison: MOMA against a pFBA template after a (feasible) knockout.
    window.medium_combo.setCurrentText("glucose_aerobic")
    window.apply_selected_medium()
    window.comparison_method_combo.setCurrentText("MOMA (L1)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["PFK"])  # reroutable aerobically
    window.run_comparison()
    assert "MOMA" in window.comparison_summary.text()
    assert window.comparison_table.rowCount() > 0

    # A lethal knockout is reported as infeasible, not a crash.
    window.medium_combo.setCurrentText("glucose_anaerobic")
    window.apply_selected_medium()
    _select_ko(window, ["PFK"])
    window.run_comparison()
    assert "infeasible" in window.comparison_summary.text()


def test_fvseof_in_gui(app, ecoli_core, monkeypatch):
    import threading

    import cmm.app.main_window as main_window_module

    # GLPK problem objects must be created and used by the same thread. Track both model
    # reconstruction and the real FVSEOF call so this test cannot regress to handing the GUI's
    # main-thread solver object to QThread (which aborts swiglpk on Windows rather than raising).
    ecoli_core.solver = "glpk"
    ecoli_core.tolerance = 1e-8
    ecoli_core.solver.configuration.presolve = True
    ecoli_core.solver.configuration.timeout = 17
    ecoli_core.solver.configuration.tolerances.integrality = 2e-8
    custom_constraint = ecoli_core.problem.Constraint(
        ecoli_core.reactions.PGI.flux_expression,
        ub=999.0,
        name="_gui_worker_constraint",
    )
    ecoli_core.add_cons_vars([custom_constraint])
    ui_model = ecoli_core
    main_tid = threading.get_ident()
    state = {}
    real_from_json = main_window_module.from_json
    real_fvseof = main_window_module.fvseof

    def tracked_from_json(document):
        state["rebuild_tid"] = threading.get_ident()
        return real_from_json(document)

    def tracked_fvseof(model, *args, **kwargs):
        state["solve_tid"] = threading.get_ident()
        state["detached"] = model is not ui_model
        state["solver"] = main_window_module.active_solver(model)
        state["tolerance"] = model.tolerance
        state["integrality"] = model.solver.configuration.tolerances.integrality
        state["presolve"] = model.solver.configuration.presolve
        state["timeout"] = model.solver.configuration.timeout
        state["objective_direction"] = model.objective_direction
        state["glucose_bounds"] = model.reactions.EX_glc__D_e.bounds
        state["custom_constraint_ub"] = model.solver.constraints[
            "_gui_worker_constraint"
        ].ub
        return real_fvseof(model, *args, **kwargs)

    monkeypatch.setattr(main_window_module, "from_json", tracked_from_json)
    monkeypatch.setattr(main_window_module, "fvseof", tracked_fvseof)

    window = CmmMainWindow(ecoli_core)
    window.tabs.setCurrentIndex(window._tab_index("Production"))
    window.product_combo.setCurrentText("EX_succ_e")
    window.anaerobic_combo.setCurrentText("anaerobic")
    window.run_fvseof_plot()
    assert (
        "FVSEOF" in window.yield_label.text() or "robust" in window.yield_label.text()
    )
    assert state["rebuild_tid"] == state["solve_tid"] != main_tid
    assert state["detached"] is True
    assert state["solver"] == "glpk"
    assert state["tolerance"] == pytest.approx(1e-8)
    assert state["integrality"] == pytest.approx(2e-8)
    assert state["presolve"] is True
    assert state["timeout"] == 17
    assert state["objective_direction"] == "max"
    assert state["glucose_bounds"] == (-10.0, 1000.0)
    assert state["custom_constraint_ub"] == 999.0
    # Worker cleanup releases only its thread-local GLPK environment; the UI model remains
    # usable on the main thread after the background solve.
    assert window.model.slim_optimize(error_value=None) is not None


def test_menu_bar_has_file_menu(app):
    window = CmmMainWindow(build_demo_model())
    titles = [a.text().replace("&", "") for a in window.menuBar().actions()]
    assert "File" in titles


def test_strain_design_tab_present_and_disabled_without_exchanges(app):
    from cobra import Metabolite, Model, Reaction

    # Strain Design tab exists on any model.
    window = CmmMainWindow(build_demo_model())
    assert window._tab_index("Strain Design") is not None

    # ...but Run design is disabled when the model has no exchange reactions.
    model = Model("no_exchange")
    a = Metabolite("a_c", compartment="c")
    b = Metabolite("b_c", compartment="c")
    r = Reaction("R")
    r.add_metabolites({a: -1, b: 1})
    model.add_reactions([r])
    model.objective = r
    window2 = CmmMainWindow(model)
    assert not window2.sd_run_btn.isEnabled()


@pytest.mark.requires_qp
def test_transformation_tab_runs(app):
    window = CmmMainWindow(build_demo_model())
    assert window._tab_index("Transform (A→B)") is not None
    assert not window.transform_run_btn.isEnabled()

    window._transform_source_expression = dict(SOURCE_EXPRESSION)
    window._transform_target_expression = dict(TARGET_EXPRESSION)
    window._update_transform_run_state()
    assert window.transform_run_btn.isEnabled()

    window.transform_method_combo.setCurrentText("moma")
    window.transform_omics_combo.setCurrentText("eflux2")
    window.run_transformation()
    # Every gene is a candidate knockout; the ranking table is populated and summarized.
    assert window.transform_table.rowCount() > 0
    assert "A→B target" in window.transform_summary.text()


def test_omics_multi_condition_in_gui(app, tmp_path):
    from cmm.omics.conditions import read_expression_table

    path = tmp_path / "conditions.csv"
    path.write_text(
        "gene,condA,condB\ng1,50,50\ng2,100,1\ng3,1,100\ng5,1,100\ngb,50,50\n"
    )

    window = CmmMainWindow(build_demo_model())
    # Multi-condition is now folded into the Omics tab.
    assert window._tab_index("Omics") is not None
    assert window._tab_index("Multi-condition") is None

    # Loading detects both conditions (checked by default) but does not compute.
    window._set_omics_source(read_expression_table(str(path)), "conditions.csv")
    assert window._checked_omics_conditions() == ["condA", "condB"]
    assert window.omics_table.rowCount() == 0

    # Compute predicts one flux column per checked condition (column name = condition).
    # LAD is an LP so it runs without a QP solver in CI.
    window.omics_method_combo.setCurrentText("lad")
    window.compute_omics()  # must not raise
    headers = [
        window.omics_table.horizontalHeaderItem(i).text()
        for i in range(window.omics_table.columnCount())
    ]
    assert headers == ["Reaction", "condA", "condB"]
    assert window.omics_table.rowCount() > 0

    # Changing the checked set and recomputing reuses the same data (no reload) and drops
    # the unchecked condition's column.
    for i in range(window.omics_cond_list.count()):
        if window.omics_cond_list.item(i).text() == "condB":
            window.omics_cond_list.item(i).setCheckState(Qt.Unchecked)
    window.compute_omics()
    headers = [
        window.omics_table.horizontalHeaderItem(i).text()
        for i in range(window.omics_table.columnCount())
    ]
    assert headers == ["Reaction", "condA"]


def test_export_table_csv(app, tmp_path):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    window.tabs.setCurrentIndex(window._tab_index("Simulation"))

    out = tmp_path / "fluxes.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (str(out), "CSV")
    )
    window.export_table_csv()
    assert out.exists()
    header = out.read_text().splitlines()[0]
    assert header.startswith("Reaction")


@pytest.mark.requires_qp
def test_comparison_gene_knockout(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("MOMA (L2)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("gene")
    _select_ko(
        window, ["g2"]
    )  # disease-branch gene -> blocks R2 -> reroute through R3/R5
    window.run_comparison()
    assert "gene g2" in window.comparison_summary.text()
    assert "1 reactions blocked" in window.comparison_summary.text()
    assert window.comparison_table.rowCount() > 0


@pytest.mark.requires_qp
def test_comparison_multi_knockout(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["R1", "R2"])
    window.run_comparison()
    # Both reactions named in the summary; knocking out R1 is lethal (no substrate route).
    text = window.comparison_summary.text()
    assert "R1" in text and "R2" in text


def test_comparison_requires_a_selection(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    _clear_ko(window)
    window.run_comparison()
    assert "Add one or more" in window.comparison_summary.text()


@pytest.mark.requires_qp
def test_comparison_batch_over_genes(app):
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("MOMA (L2)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("gene")
    _clear_ko(window)  # nothing selected -> batch over all genes
    window.run_batch_comparison()
    assert window.comparison_table.columnCount() == 4  # no target product selected
    assert window.comparison_table.rowCount() == len(window.model.genes)
    assert "Batch" in window.comparison_summary.text()
    # The batch table reports wild-type vs post-knockout biomass and essentiality per target.
    headers = [
        window.comparison_table.horizontalHeaderItem(c).text()
        for c in range(window.comparison_table.columnCount())
    ]
    assert headers == ["Target", "WT Biomass", "KO Biomass", "Essential"]


def test_background_execution_runs_off_thread_and_stays_responsive(app):
    import threading
    import time

    from qtpy.QtCore import QTimer

    window = CmmMainWindow(build_demo_model())
    main_tid = threading.get_ident()
    state = {"pumped": False, "worker_tid": None, "running": False}

    # A main-thread timer that should fire WHILE the worker computes — proving the UI event
    # loop keeps running (the window would repaint) instead of freezing on the solve.
    def on_timer():
        if state["running"]:
            state["pumped"] = True

    QTimer.singleShot(40, on_timer)

    def compute():
        state["worker_tid"] = threading.get_ident()
        state["running"] = True
        time.sleep(0.3)
        state["running"] = False
        return "done"

    result = window._run_in_background(compute, label="test")
    assert result == "done"  # synchronous: result ready on return
    assert state["worker_tid"] != main_tid  # ran off the UI thread
    assert state["pumped"]  # the UI event loop pumped during the compute


def test_background_execution_reraises_errors(app):
    window = CmmMainWindow(build_demo_model())

    def boom():
        raise ValueError("kaboom")

    with pytest.raises(ValueError, match="kaboom"):
        window._run_in_background(boom, label="test")

    # The busy flag is cleared even when the compute raised, so the window stays usable.
    assert window._busy is False


def test_window_renders_non_blank(app):
    window = CmmMainWindow(build_demo_model())
    window.run_fba()
    window.show()
    app.processEvents()

    pixmap = window.grab()
    assert pixmap.width() > 400
    assert pixmap.height() > 300

    image = pixmap.toImage()
    assert not image.isNull()
    # Sample pixels across the canvas; a blank render would be a single uniform color.
    colors = {
        image.pixel(x, y)
        for x in range(0, image.width(), max(1, image.width() // 20))
        for y in range(0, image.height(), max(1, image.height() // 20))
    }
    assert len(colors) > 5


def test_flux_response_tab_runs(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    window.fr_response_combo.setCurrentText("EX_succ_e")
    window.fr_growth_spin.setValue(30.0)
    window.fr_steps_spin.setValue(10)
    window.run_flux_response()

    assert window.fr_table.rowCount() == 10
    assert "EX_succ_e" in window.fr_summary.text()
    assert "PGI" in window.fr_summary.text()
    # The scan table is the tab's export target and the curve is its figure.
    assert window._active_table() is window.fr_table
    assert window._active_figure() is not None


def test_flux_response_tab_reports_infeasible_points(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    window.fr_response_combo.setCurrentText("EX_succ_e")
    window.fr_growth_spin.setValue(30.0)
    window.fr_steps_spin.setValue(12)
    window.run_flux_response()

    statuses = [
        window.fr_table.item(r, 3).text() for r in range(window.fr_table.rowCount())
    ]
    # A growth floor makes the extremes of the range unsolvable; they are shown, not hidden.
    assert "infeasible" in statuses
    infeasible_row = statuses.index("infeasible")
    assert window.fr_table.item(infeasible_row, 1).text() == "—"


def test_flux_response_tab_surfaces_errors(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("NOT_A_REACTION")
    window.run_flux_response()
    assert "failed" in window.fr_summary.text().lower()


def test_flux_response_manual_range(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    # The range boxes are always editable; narrowing them narrows the scan.
    assert window.fr_min_spin.isEnabled()
    window.fr_min_spin.setValue(-5.0)
    window.fr_max_spin.setValue(5.0)
    window.fr_steps_spin.setValue(6)
    window.run_flux_response()
    assert window.fr_table.rowCount() == 6
    assert float(window.fr_table.item(0, 0).text()) == pytest.approx(-5.0, abs=1e-6)


def test_sampling_tab_runs(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.sample_seed_spin.setValue(4)
    window.run_sampling()

    assert window.sample_table.rowCount() == len(ecoli_core.reactions)
    assert "120" in window.sample_summary.text()
    assert "seed 4" in window.sample_summary.text()
    assert window._active_table() is window.sample_table
    assert window._active_figure() is not None


def test_sampling_tab_reference_mode(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    # The reference controls stay disabled until the mode needs them.
    assert not window.sample_reference_combo.isEnabled()
    window.sample_mode_combo.setCurrentText("around a reference")
    assert window.sample_reference_combo.isEnabled()
    assert window.sample_window_spin.isEnabled()

    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()
    assert window.sample_table.rowCount() == len(ecoli_core.reactions)


def test_sampling_table_exports_to_csv(app, ecoli_core, tmp_path):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()

    out = tmp_path / "samples.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (str(out), "CSV")
    )
    window.export_table_csv()
    assert out.exists()
    assert out.read_text().splitlines()[0].startswith("Reaction,Mean,Std")


def test_sampling_exports_raw_ensemble(app, ecoli_core, tmp_path):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    # Nothing to export until a run has happened.
    assert not window.sample_export_btn.isEnabled()

    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()
    assert window.sample_export_btn.isEnabled()

    out = tmp_path / "ensemble.csv"
    from qtpy import QtWidgets

    QtWidgets.QFileDialog.getSaveFileName = staticmethod(
        lambda *a, **k: (str(out), "CSV")
    )
    window.export_sampling_samples()
    assert out.exists()

    import pandas as pd

    frame = pd.read_csv(out, index_col="sample")
    # One row per drawn sample, one column per reaction — the raw ensemble, not the summary.
    assert frame.shape == (120, len(ecoli_core.reactions))


def test_sampling_table_covers_every_reaction(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.sample_top_spin.setValue(5)
    window.run_sampling()

    # The statistics table is complete; only the figure is restricted to the top-N.
    assert window.sample_table.rowCount() == len(ecoli_core.reactions)
    listed = {
        window.sample_table.item(r, 0).text()
        for r in range(window.sample_table.rowCount())
    }
    assert listed == {r.id for r in ecoli_core.reactions}
    assert len(window._sampling_canvas.figure.axes[0].get_yticklabels()) == 5


def test_reloading_a_model_clears_the_ensemble(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_n_spin.setValue(120)
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()
    assert window._sampling_result is not None

    window.load_model(build_demo_model())
    assert window._sampling_result is None
    assert not window.sample_export_btn.isEnabled()
    assert window.sample_table.rowCount() == 0


def test_range_is_seeded_with_the_theoretical_interval(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")

    # Not a placeholder: PGI can carry -50..10 in this medium.
    assert window.fr_min_spin.value() == pytest.approx(-50.0, abs=1e-3)
    assert window.fr_max_spin.value() == pytest.approx(10.0, abs=1e-3)

    # Switching target re-seeds rather than keeping the previous reaction's numbers.
    window.fr_target_combo.setCurrentText("EX_o2_e")
    assert window.fr_min_spin.value() == pytest.approx(-60.0, abs=1e-3)
    assert window.fr_max_spin.value() == pytest.approx(0.0, abs=1e-3)

    # A manual edit survives until the user asks for a reset.
    window.fr_min_spin.setValue(-5.0)
    window._prefill_fr_range()
    assert window.fr_min_spin.value() == pytest.approx(-60.0, abs=1e-3)


def test_reference_controls_are_disabled_in_uniform_mode(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    reference_widgets = (
        window.sample_reference_label,
        window.sample_reference_combo,
        window.sample_window_label,
        window.sample_window_spin,
    )

    assert window.sample_mode_combo.currentText() == "uniform"
    # Uniform sampling has no reference, so neither the combo nor its label invites input.
    assert not any(widget.isEnabled() for widget in reference_widgets)

    window.sample_mode_combo.setCurrentText("around a reference")
    assert all(widget.isEnabled() for widget in reference_widgets)

    window.sample_mode_combo.setCurrentText("uniform")
    assert not any(widget.isEnabled() for widget in reference_widgets)


def test_tab_order_groups_related_analyses(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    names = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert names == [
        "Simulation",
        "Sampling",
        "Comparison",
        "Flux Response",
        "Production",
        "Strain Design",
        "Omics",
        "Flux Map",
        "Revert Metabolism",
        "Transform (A→B)",
        "Agent",
    ]
    # Each new tab sits beside the analysis asking a similar question.
    assert names.index("Sampling") == names.index("Simulation") + 1
    assert names.index("Flux Response") == names.index("Comparison") + 1


def _run_sampling(window, **overrides):
    window.sample_n_spin.setValue(overrides.get("n", 200))
    window.sample_method_combo.setCurrentText("achr")
    window.sample_thinning_spin.setValue(10)
    window.run_sampling()


def test_sampling_applies_reaction_knockouts(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")

    _run_sampling(window)
    wild_type = window._sampling_result.samples["PGI"].mean()
    assert wild_type > 1e-6
    assert "wild type" in window.sample_summary.text()

    window._ko_add("sample_ko", ["PGI"])
    _run_sampling(window)
    # The deleted reaction carries no flux anywhere in the ensemble.
    assert window._sampling_result.samples["PGI"].abs().max() < 1e-6
    assert "PGI" in window.sample_summary.text()
    assert window._sampling_result.metadata["parameters"]["condition"] == "knockout:PGI"
    # The knockout is a scoped condition, not an edit to the loaded model.
    assert ecoli_core.reactions.PGI.bounds == (-1000.0, 1000.0)


def test_sampling_resolves_gene_knockouts_through_the_gpr(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_ko_level_combo.setCurrentText("gene")
    # b4025 is pgi; the GPR must resolve it to the PGI reaction.
    window._ko_add("sample_ko", ["b4025"])
    _run_sampling(window)
    assert window._sampling_result.samples["PGI"].abs().max() < 1e-6
    assert "gene knockout" in window.sample_summary.text()


def test_sampling_level_switch_clears_the_selection(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window._ko_add("sample_ko", ["PGI"])
    assert window._ko_targets("sample_ko") == ["PGI"]
    # Reaction ids are not gene ids, so switching level must not carry them over.
    window.sample_ko_level_combo.setCurrentText("gene")
    assert window._ko_targets("sample_ko") == []
    assert window.sample_ko_selected_label.text() == "Selected (0):"


def test_reference_mode_uses_a_knocked_out_reference(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    window.sample_mode_combo.setCurrentText("around a reference")
    window._ko_add("sample_ko", ["PGI"])
    # A wild-type reference would put PGI outside its own window and raise; the reference
    # has to be computed under the same knockouts.
    _run_sampling(window)
    assert window._sampling_result is not None
    assert window._sampling_result.samples["PGI"].abs().max() < 1e-6


def test_lethal_knockout_is_reported_not_crashed(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Sampling")
    _run_sampling(window)
    assert window._sampling_result is not None

    window.sample_ko_level_combo.setCurrentText("gene")
    window._ko_add("sample_ko", ["b2415"])  # essential in e_coli_core
    _run_sampling(window)

    assert "lethal" in window.sample_summary.text()
    # The previous run's ensemble must not survive as an exportable result.
    assert window._sampling_result is None
    assert not window.sample_export_btn.isEnabled()
    assert window.sample_table.rowCount() == 0


def test_knockout_picker_is_independent_per_tab(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._add_ko_targets(["PGI"])  # Comparison tab picker
    window._ko_add("sample_ko", ["ENO", "FBA"])  # Sampling tab picker
    assert window._selected_ko_targets() == ["PGI"]
    assert window._ko_targets("sample_ko") == ["ENO", "FBA"]
    window._clear_ko_selected()
    assert window._selected_ko_targets() == []
    assert window._ko_targets("sample_ko") == ["ENO", "FBA"]


def test_room_comparison_renders_a_switch_count_not_a_distance(app):
    """Regression: ROOM's ``distance`` is None, and formatting it with ``:.4g`` raised.

    ``comparison.py`` split ``distance`` from the raw solver objective in 0.4.0, and ROOM's
    objective is a *count of switched reactions*, not a distance — so ``distance`` is None
    there. The GUI offers ROOM, and its summary formatted ``result.distance`` unconditionally,
    which raised ``TypeError: unsupported format string passed to NoneType.__format__`` on
    that path. The suite did not catch it because it only ever exercised MOMA.
    """

    # This is a GUI-labelling regression, so use the six-reaction demo instead of turning the
    # check into a slow genome-scale MILP on the open-source CI solver.
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("ROOM")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["R2"])
    window.run_comparison()  # must not raise

    summary = window.comparison_summary.text()
    assert "ROOM" in summary
    assert "reactions switched" in summary
    assert "not a distance" in summary
    # The tolerance pair that produced the count is stated, because the two published pairs
    # differ by 24% on the same screen.
    assert "flux_prediction" in summary
    assert window.comparison_table.rowCount() > 0


@pytest.mark.requires_qp
def test_moma_comparison_names_the_quantity_it_reports(app):
    """MOMA-L2 reports Segrè Eq. (4)'s Euclidean distance, with the QP objective alongside."""

    # Labelling is independent of network scale.  The six-reaction branched model still
    # exercises a nonzero L2 reroute while fitting Gurobi's bundled restricted license.
    window = CmmMainWindow(build_demo_model())
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("MOMA (L2)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["R2"])
    window.run_comparison()

    summary = window.comparison_summary.text()
    assert "Euclidean (L2) distance" in summary
    assert "QP objective" in summary


def test_room_tolerance_selector_is_only_live_for_room(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Comparison")
    window.comparison_method_combo.setCurrentText("MOMA (L1)")
    assert not window.room_use_case_combo.isEnabled()
    window.comparison_method_combo.setCurrentText("ROOM")
    assert window.room_use_case_combo.isEnabled()
    assert set(
        window.room_use_case_combo.itemText(i)
        for i in range(window.room_use_case_combo.count())
    ) == set(ROOM_TOLERANCES)


def test_gui_reference_template_default_matches_the_library(app, ecoli_core):
    """Plan item 2.12: the GUI defaulted to FBA while ``reference_flux`` defaults to pFBA."""

    import inspect

    window = CmmMainWindow(ecoli_core)
    assert window.template_combo.currentText() == "pfba"
    assert inspect.signature(reference_flux).parameters["method"].default == "pfba"


def test_flux_response_tab_shows_phases_and_the_response_limit(app, ecoli_core):
    """``ResponseBottleneck`` is deleted; the tab presents phases and the shadow-price limit."""

    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Flux Response")
    window.fr_target_combo.setCurrentText("PGI")
    window.fr_response_combo.setCurrentText("EX_succ_e")
    window.fr_growth_spin.setValue(30.0)
    window.fr_steps_spin.setValue(10)
    window.run_flux_response()

    summary = window.fr_summary.text()
    assert "bottleneck" not in summary.lower()
    assert "shadow price" in summary
    assert "phase" in summary
    # One row per phase, and the shadow price shown is the exact one from the phase table.
    assert window.fr_phase_table.rowCount() > 0


def test_flux_response_limit_does_not_move_with_the_step_count(app, ecoli_core):
    """The GUI must show a location that is a property of the network, not of the grid."""

    rendered = set()
    for steps in (6, 20, 60):
        window = CmmMainWindow(ecoli_core)
        window._goto_tab("Flux Response")
        window.fr_target_combo.setCurrentText("PGI")
        window.fr_response_combo.setCurrentText("EX_succ_e")
        window.fr_growth_spin.setValue(30.0)
        window.fr_steps_spin.setValue(steps)
        window.run_flux_response()
        # Everything from the limit sentence onwards, plus the whole phase table. The
        # "best response at target X" clause before it is a property of the grid and is
        # allowed to move; the limit and the phases are not.
        summary = window.fr_summary.text()
        limit_clause = summary[summary.index("Response falls beyond") :]
        table = tuple(
            window.fr_phase_table.item(row, col).text()
            for row in range(window.fr_phase_table.rowCount())
            for col in range(window.fr_phase_table.columnCount())
        )
        rendered.add((limit_clause, table))
    # The old ``bottleneck`` field moved by up to 29.5 flux units over this same sweep and
    # its ``found`` flag inverted; the shadow-price limit renders identically at every n.
    assert len(rendered) == 1


def test_production_tab_states_its_aeration_as_a_condition(app, ecoli_core):
    """``aerobic=`` was removed in 0.4.0; the GUI passes a Condition like every other caller."""

    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Production")
    window.product_combo.setCurrentText("EX_succ_e")
    window.anaerobic_combo.setCurrentText("anaerobic")
    window.run_theoretical_yield()

    text = window.yield_label.text()
    assert "anaerobic" in text
    # 1.200 exactly - the CO2-closed value. 1.391 would mean the GUI lost the CO2 closure
    # the removed flag used to apply for it.
    assert "1.200" in text

    condition = window._production_condition()
    assert {b.reaction_id for b in condition.bounds} == {"EX_o2_e", "EX_co2_e"}
    assert all(b.lower_bound == 0.0 for b in condition.bounds)

    window.anaerobic_combo.setCurrentText("aerobic")
    assert window._production_condition().bounds == ()


def test_strain_design_tab_passes_its_aeration_to_the_search(app, ecoli_core):
    """optknock/robustknock accepted no condition before 0.4.0 — the SC-01 defect class."""

    window = CmmMainWindow(ecoli_core)
    window._goto_tab("Strain Design")
    window.sd_anaerobic_combo.setCurrentText("anaerobic")
    condition = window._strain_design_condition()
    assert condition.name.startswith("anaerobic")
    assert {b.reaction_id for b in condition.bounds} == {"EX_o2_e", "EX_co2_e"}


def test_applying_a_medium_reports_what_was_dropped(app, ecoli_core):
    """``apply_to`` returns a MediumApplication; a preset that drops 18 of 24 must say so."""

    window = CmmMainWindow(ecoli_core)
    window.medium_combo.setCurrentText("glucose_aerobic")
    with pytest.warns(UserWarning, match="were not applied"):
        window.apply_selected_medium()
    status = window.status_label.text()
    assert "Glucose minimal, aerobic" in status  # the 0.4.0 display name, not the key
    assert "dropped" in status
    assert "EX_fe2_e" in status


# --- Flux Map tab ---------------------------------------------------------------------
# The map is the one analysis that needs an external layout file, so the tab has to behave
# when that file is absent, wrong, or for another model.


def _flux_map_window(model):
    from cmm.app.main_window import _MAP_LAYOUT_ESCHER, _MAP_LAYOUT_SCHEMATIC

    window = CmmMainWindow(model)
    window._goto_tab("Flux Map")
    return window, _MAP_LAYOUT_ESCHER, _MAP_LAYOUT_SCHEMATIC


def test_flux_map_tab_opens_on_the_bundled_map_without_any_setup(app, ecoli_core):
    """The default model must show a curated map with no file to find and no flag to pass."""

    window, escher, _ = _flux_map_window(ecoli_core)
    assert window._map_path is not None
    assert window.map_layout_combo.currentText() == escher
    assert "95 reactions" in window.map_source_label.text()

    window.render_flux_map()
    assert window._map_canvas is not None
    assert "Escher layout" in window.status_label.text()


def test_flux_map_falls_back_to_the_schematic_when_no_map_fits(app):
    """A model no bundled map describes still gets a figure — and is told what it is."""

    import cobra

    # Two compartments, because cobra locates the external one before listing exchanges.
    model = cobra.Model("unmapped")
    a_e = cobra.Metabolite("a_e", compartment="e")
    a_c = cobra.Metabolite("a_c", compartment="c")
    b_c = cobra.Metabolite("b_c", compartment="c")
    b_e = cobra.Metabolite("b_e", compartment="e")

    def _reaction(rid, stoichiometry, lower=0.0):
        reaction = cobra.Reaction(rid, lower_bound=lower, upper_bound=1000.0)
        reaction.add_metabolites(stoichiometry)
        return reaction

    model.add_reactions(
        [
            _reaction("EX_a_e", {a_e: -1}, lower=-10.0),
            _reaction("A_TRANSPORT", {a_e: -1, a_c: 1}),
            _reaction("GROW", {a_c: -1, b_c: 1}),
            _reaction("B_TRANSPORT", {b_c: -1, b_e: 1}),
            _reaction("EX_b_e", {b_e: -1}),
        ]
    )
    model.objective = "GROW"

    window, _, schematic = _flux_map_window(model)
    assert window._map_path is None
    assert window.map_layout_combo.currentText() == schematic
    assert "No bundled Escher map" in window.map_source_label.text()

    window.render_flux_map()
    assert window._map_canvas is not None
    assert "schematic" in window.status_label.text()


def test_flux_map_switches_between_layouts(app, ecoli_core):
    window, escher, schematic = _flux_map_window(ecoli_core)

    window.map_layout_combo.setCurrentText(schematic)
    window.map_topn_spin.setValue(8)
    window.render_flux_map()
    assert "top 8 reactions" in window.status_label.text()

    window.map_layout_combo.setCurrentText(escher)
    window.render_flux_map()
    assert "Escher layout" in window.status_label.text()


def test_loading_a_map_for_a_different_model_is_refused(app, ecoli_core, tmp_path):
    """An all-grey map is worse than no map: say the ids do not match, draw nothing."""

    from cmm.app import main_window as mw

    foreign = tmp_path / "foreign.json"
    foreign.write_text(
        json.dumps(
            [
                {"map_name": "foreign", "schema": "…1-0-0#"},
                {
                    "reactions": {"1": {"bigg_id": "NOT_IN_ANY_ECOLI_MODEL"}},
                    "nodes": {},
                },
            ]
        )
    )
    window, _, _ = _flux_map_window(ecoli_core)
    before = window._map_path

    monkey = lambda *a, **k: (str(foreign), "")  # noqa: E731 - one-line dialog stub
    original = mw.QFileDialog.getOpenFileName
    mw.QFileDialog.getOpenFileName = staticmethod(monkey)
    try:
        window.open_flux_map_dialog()
    finally:
        mw.QFileDialog.getOpenFileName = original

    assert window._map_path == before  # the good map was not replaced by the bad one
    assert "different model" in window.status_label.text()


def test_loading_a_non_escher_file_is_refused(app, ecoli_core, tmp_path):
    from cmm.app import main_window as mw

    junk = tmp_path / "junk.json"
    junk.write_text("{ not json at all")
    window, _, _ = _flux_map_window(ecoli_core)
    before = window._map_path

    original = mw.QFileDialog.getOpenFileName
    mw.QFileDialog.getOpenFileName = staticmethod(lambda *a, **k: (str(junk), ""))
    try:
        window.open_flux_map_dialog()
    finally:
        mw.QFileDialog.getOpenFileName = original

    assert window._map_path == before
    assert "not a readable escher map" in window.status_label.text().lower()


def test_map_label_describes_the_layout_actually_shown(app, ecoli_core):
    """In schematic mode the caption must not still be describing the Escher map."""

    window, escher, schematic = _flux_map_window(ecoli_core)
    assert window.map_source_label.text().startswith("Escher map:")

    window.map_layout_combo.setCurrentText(schematic)
    caption = window.map_source_label.text()
    assert caption.startswith("Schematic layout")
    assert "No layout file needed" in caption
    assert "available" in caption  # the Escher map is still offered, not hidden

    window.map_layout_combo.setCurrentText(escher)
    assert window.map_source_label.text().startswith("Escher map:")


def test_flux_map_title_names_the_method_behind_the_numbers(app, ecoli_core):
    """A pFBA map read as an FBA map is a wrong figure, so the title has to say which."""

    window, _, _ = _flux_map_window(ecoli_core)

    window.run_fba()
    window.render_flux_map()
    assert window._flux_source == "FBA"
    assert "FBA" in window._map_canvas.figure.axes[0].get_title()

    window.run_pfba()
    window.render_flux_map()
    assert window._flux_source == "pFBA"
    title = window._map_canvas.figure.axes[0].get_title()
    assert "pFBA" in title and "FBA," not in title


def test_omics_flux_state_can_be_drawn_on_the_map(app, tmp_path):
    """The expression-derived distribution is what the map is most useful for."""

    from cmm.omics.conditions import read_expression_table

    path = tmp_path / "conditions.csv"
    path.write_text(
        "gene,condA,condB\ng1,50,50\ng2,100,1\ng3,1,100\ng5,1,100\ngb,50,50\n"
    )

    window = CmmMainWindow(build_demo_model())
    assert not window.omics_map_btn.isEnabled()  # nothing computed yet

    window._set_omics_source(read_expression_table(str(path)), "conditions.csv")
    window.omics_method_combo.setCurrentText("lad")  # an LP, so no QP solver needed
    window.compute_omics()

    assert window.omics_map_btn.isEnabled()
    assert [
        window.omics_map_combo.itemText(i)
        for i in range(window.omics_map_combo.count())
    ] == ["condA", "condB"]

    window.omics_map_combo.setCurrentText("condB")
    window.show_omics_on_flux_map()

    assert window.tabs.tabText(window.tabs.currentIndex()) == "Flux Map"
    assert window._flux_source == "LAD · condB"
    assert "LAD" in window._map_canvas.figure.axes[0].get_title()
    # The distribution drawn is the omics one, not a quietly re-run FBA.
    assert window._fluxes == dict(window._omics_fluxes_by_condition["condB"])


def test_knockout_redistribution_can_be_drawn_on_the_map(app, ecoli_core):
    window = CmmMainWindow(ecoli_core)
    assert not window.cmp_map_btn.isEnabled()

    window.medium_combo.setCurrentText("glucose_aerobic")
    window.apply_selected_medium()
    window.comparison_method_combo.setCurrentText("MOMA (L1)")
    window.template_combo.setCurrentText("pfba")
    window.ko_level_combo.setCurrentText("reaction")
    _select_ko(window, ["PFK"])
    window.run_comparison()

    assert window.cmp_map_btn.isEnabled()
    window.show_comparison_on_flux_map()

    assert window.tabs.tabText(window.tabs.currentIndex()) == "Flux Map"
    assert "MOMA" in window._flux_source and "PFK" in window._flux_source
    assert window._fluxes == dict(window._comparison_cache["fluxes"])
    assert "MOMA" in window._map_canvas.figure.axes[0].get_title()


def test_a_new_model_forgets_which_flux_state_was_drawn(app, ecoli_core):
    window, _, _ = _flux_map_window(ecoli_core)
    window.run_fba()
    assert window._flux_source == "FBA"

    window.load_model(build_demo_model())
    assert window._flux_source == ""
    assert not window._fluxes


def test_flux_map_tab_can_run_either_simulation_itself(app, ecoli_core):
    """pFBA was unreachable from the map: it had to be run on another tab first."""

    window, _, _ = _flux_map_window(ecoli_core)
    assert not window._fluxes  # nothing solved yet

    window._draw_simulation_on_map(window.run_pfba)
    assert window._flux_source == "pFBA"
    assert "pFBA" in window._map_canvas.figure.axes[0].get_title()
    pfba = dict(window._fluxes)

    window._draw_simulation_on_map(window.run_fba)
    assert window._flux_source == "FBA"
    assert "FBA" in window._map_canvas.figure.axes[0].get_title()
    # pFBA minimises total flux subject to the same optimum, so it can never exceed FBA's.
    # It can equal it — on this medium the FBA vertex is already parsimonious — which is
    # exactly why the title has to name the method rather than leave it to the numbers.
    assert (
        sum(abs(v) for v in pfba.values())
        <= sum(abs(v) for v in window._fluxes.values()) + 1e-6
    )


def test_both_tabs_offer_the_map_under_the_same_name(app, ecoli_core):
    """One action, one label — 'Draw' on one tab and 'Show on flux map' on another is two."""

    window = CmmMainWindow(ecoli_core)
    assert (
        window.omics_map_btn.text() == window.cmp_map_btn.text() == "Show on flux map"
    )


def test_changing_the_layout_redraws_without_re_solving(app, ecoli_core):
    """The reason there is no render button: a display change must not touch the numbers.

    Removing it would be wrong if switching layout meant pressing "FBA" to redraw — that
    quietly replaces an omics or knockout flux state with a fresh FBA solve. The display
    controls redraw what is loaded and never solve.
    """

    window, escher, schematic = _flux_map_window(ecoli_core)
    window._draw_simulation_on_map(window.run_pfba)
    loaded = dict(window._fluxes)
    assert window._flux_source == "pFBA"

    window.map_layout_combo.setCurrentText(schematic)
    assert window._map_redraw_timer.isActive(), "layout change scheduled no redraw"
    window._map_redraw_timer.stop()
    window.render_flux_map()
    assert window._flux_source == "pFBA"  # not silently re-solved as FBA
    assert window._fluxes == loaded

    window.map_topn_spin.setValue(18)
    assert window._map_redraw_timer.isActive(), "reaction count changed nothing"
    window._map_redraw_timer.stop()
    window.render_flux_map()
    assert window._fluxes == loaded
    assert "top 18" in window._map_canvas.figure.axes[0].get_title()


def test_display_controls_do_not_solve_before_anything_is_drawn(app, ecoli_core):
    """With no figure yet, a layout change must not quietly run FBA to have something to draw."""

    window, _, schematic = _flux_map_window(ecoli_core)
    assert window._map_canvas is None

    window.map_layout_combo.setCurrentText(schematic)
    assert not window._map_redraw_timer.isActive()
    assert not window._fluxes
    assert window._flux_source == ""


def test_map_background_loads_toggles_and_is_dropped_with_its_map(
    app, ecoli_core, tmp_path
):
    """A drawing describes one map, so it must not survive that map being replaced."""

    from cmm.app.svg_background import svg_background

    drawing = tmp_path / "drawing.svg"
    drawing.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100" '
        'viewBox="0 0 100 100"><rect width="100" height="100" fill="#8ab"/></svg>'
    )

    window, _, _ = _flux_map_window(ecoli_core)
    window._draw_simulation_on_map(window.run_fba)
    assert window._map_background is None
    assert not window.map_background_check.isVisible()

    window._map_background = svg_background(drawing)
    window.map_background_check.setVisible(True)
    window.map_background_check.setChecked(True)
    window.render_flux_map()
    assert len(window._map_canvas.figure.axes[0].images) == 1
    assert "over the loaded drawing" in window.status_label.text()

    # Unchecking hides it without unloading — the array is still there.
    window.map_background_check.setChecked(False)
    window.render_flux_map()
    assert not window._map_canvas.figure.axes[0].images
    assert window._map_background is not None

    # A different model drops the drawing along with the map it described.
    window.load_model(build_demo_model())
    assert window._map_background is None


def test_omics_condition_chooser_says_how_many_it_holds(app, tmp_path):
    """A closed combo shows one item while the checklist above shows every condition.

    Reported as a bug — two conditions computed, one visible. The box was right; nothing said
    the other was one click away. The label now counts them and tracks the selection.
    """

    from cmm.omics.conditions import read_expression_table

    path = tmp_path / "conditions.csv"
    path.write_text(
        "gene,condA,condB\ng1,50,50\ng2,100,1\ng3,1,100\ng5,1,100\ngb,50,50\n"
    )

    window = CmmMainWindow(build_demo_model())
    window._set_omics_source(read_expression_table(str(path)), "conditions.csv")
    window.omics_method_combo.setCurrentText("lad")
    window.compute_omics()

    assert window.omics_map_combo.count() == 2
    assert window.omics_map_label.text() == "Draw (1 of 2):"

    window.omics_map_combo.setCurrentIndex(1)
    assert window.omics_map_label.text() == "Draw (2 of 2):"
    assert window.omics_map_combo.currentText() == "condB"

    # One condition needs no count, and reloading the source resets the label with the box.
    window._set_omics_source(read_expression_table(str(path)), "conditions.csv")
    assert window.omics_map_combo.count() == 0
    assert window.omics_map_label.text() == "Draw:"
