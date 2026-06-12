"""Execute the three GUI scenario harnesses under pytest so they cannot break silently.

success-criteria G5/R2-G4 make the screenshot scenarios a PASS/FAIL gate ("any scenario
errors or produces a blank capture" = FAIL), but nothing ran them in CI: a regression in
``run_scenarios`` / ``run_gui_scenarios`` or a blank capture would not have failed pytest.
These tests run each scenario offscreen, assert the captures are non-trivial, and pin the
headline succinate claim (the anaerobic + byproduct-KO bound edit raises succinate ~0 -> ~10).
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest  # noqa: E402

pytest.importorskip("qtpy")

from cobra.io import load_model  # noqa: E402
from qtpy.QtWidgets import QApplication  # noqa: E402

from cmm.app import genome_scale_scenario, screenshots, succinate_scenario  # noqa: E402
from cmm.core import fba  # noqa: E402

# A rendered workbench grab is tens of KB; a blank/solid-colour PNG compresses to ~1-2 KB.
# This threshold is a non-blank proxy: well above blank, well below the real ~60 KB+ captures.
_MIN_PNG_BYTES = 5000


@pytest.fixture(scope="module")
def app():
    try:
        return QApplication.instance() or QApplication([])
    except Exception as exc:  # pragma: no cover - no usable Qt platform
        pytest.skip(f"no Qt platform available: {exc}")


def _assert_non_blank(paths) -> None:
    assert paths, "scenario produced no captures"
    for path in paths:
        assert path.exists(), f"missing capture: {path}"
        size = path.stat().st_size
        assert size > _MIN_PNG_BYTES, f"capture looks blank ({size} bytes): {path}"


def test_screenshots_scenario_runs(app, tmp_path, monkeypatch):
    """Branched-demo scenario: load, FBA, FVA, flux slider, rMTA, MTA-MIQP (6 captures)."""

    monkeypatch.setattr(screenshots, "OUTPUT_DIR", tmp_path)
    saved = screenshots.run_scenarios()
    assert len(saved) == 6
    _assert_non_blank(saved)


def test_succinate_gui_scenario_runs(app, tmp_path, monkeypatch):
    """e_coli_core succinate scenario: yield, envelope, FSEOF, bound edit, Escher flux map."""

    monkeypatch.setattr(succinate_scenario, "OUTPUT_DIR", tmp_path)
    saved = succinate_scenario.run_gui_scenarios()
    _assert_non_blank(saved)


def test_succinate_design_increases_succinate():
    """Pin the headline claim: the scenario's engineering bounds force succinate 0 -> ~10.

    The scenario only *prints* this; assert it here using the scenario's own ENGINEERING_BOUNDS
    so a regression (set-bounds no-op, or a design that stops excreting succinate) fails CI.
    """

    model = load_model("textbook")
    wt = fba(model).fluxes.get(succinate_scenario.SUCC, 0.0)
    with model:
        for rid, bounds in succinate_scenario.ENGINEERING_BOUNDS.items():
            model.reactions.get_by_id(rid).bounds = bounds
        design = fba(model).fluxes.get(succinate_scenario.SUCC, 0.0)
    assert wt == pytest.approx(0.0, abs=1e-6)
    assert design > 5.0, f"engineered succinate not produced: {design}"


def test_genome_scale_scenario_runs(app, tmp_path, monkeypatch):
    """Genome-scale GUI scenario on the default textbook model: FBA, glucose edit, yield,
    envelope, E-Flux2 omics."""

    monkeypatch.setattr(genome_scale_scenario, "OUTPUT_DIR", tmp_path)
    saved = genome_scale_scenario.run_scenarios(model_path=None, prefix="test_genome")
    _assert_non_blank(saved)
