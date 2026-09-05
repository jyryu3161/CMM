"""Rendering and validating a transformation run, without needing a solver.

The end-to-end test in ``test_transformation_workflow.py`` needs MIQP and therefore skips
wherever Gurobi or CPLEX is absent — which is every CI runner. That would leave the checked-in
R renderer and the completion gate untested in exactly the environment that is supposed to
prove they work, so these build the bundle by hand instead. The bundle is the renderer's only
input, so a hand-written one exercises the same path a real run does.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

import pytest

from cmm.reporting import (
    render_transformation_report,
    validate_transformation_run,
)
from cmm.reporting.transformation import (
    TransformationReportError,
    renderer_script_path,
)

RANKING = """target_id,score,rank,bTS,mTS,wTS
PGI,3.5934,1,0.1369,0.1312,-0.1369
PFK,2.0970,2,0.1041,0.1007,-0.1040
TPI,1.9430,3,0.0975,0.0942,-0.0974
ENO,0.4210,4,0.0210,0.0203,-0.0210
"""
BASELINE = """target_id,moma_score,rank
ENO,-3.5072,1
PGI,-14.8486,2
PFK,-50.0562,3
TPI,-61.2200,4
"""
SWEEP = """epsilon,target_id,score,rank,rank_at_configured_epsilon
0.005,PGI,3.5934,1,1
0.005,PFK,2.0972,2,2
0.005,TPI,1.9431,3,3
0.02,PGI,3.5930,1,1
0.02,TPI,2.0100,2,3
0.02,PFK,1.9800,3,2
"""


def _r_is_ready() -> bool:
    """Whether this machine can run the checked-in renderer at all."""

    import subprocess

    if shutil.which("Rscript") is None:
        return False
    probe = (
        'q(status = as.integer(!all(vapply(c("jsonlite", "ggplot2", "ggrepel", '
        '"svglite", "ragg"), requireNamespace, logical(1), quietly = TRUE))))'
    )
    return subprocess.run(["Rscript", "--vanilla", "-e", probe]).returncode == 0


def _write(root: Path, relative: str, content: str) -> dict[str, object]:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    raw = path.read_bytes()
    return {
        "path": relative,
        "stage": relative.split("/")[0] if "/" in relative else "root",
        "media_type": "text/csv" if relative.endswith(".csv") else "application/json",
        "status": "complete",
        "sha256": hashlib.sha256(raw).hexdigest(),
        "size_bytes": len(raw),
    }


def _rewrite_artifact(root: Path, relative: str, content: str) -> None:
    """Replace an artifact and re-record its hash, as a real re-run would.

    Editing a CSV without updating the manifest is exactly the corruption the completion gate
    exists to catch, so a fixture that did it would be testing the gate, not the renderer.
    """

    (root / relative).write_text(content, encoding="utf-8")
    manifest_path = root / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    raw = (root / relative).read_bytes()
    for entry in manifest["artifacts"].values():
        if entry["path"] == relative:
            entry["sha256"] = hashlib.sha256(raw).hexdigest()
            entry["size_bytes"] = len(raw)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )


@pytest.fixture()
def bundle(tmp_path: Path) -> Path:
    """A minimal but complete transformation run, written directly to disk."""

    root = tmp_path / "run"
    root.mkdir()
    artifacts: dict[str, dict[str, object]] = {}

    def add(role: str, relative: str, content: str) -> None:
        entry = _write(root, relative, content)
        entry["role"] = role
        artifacts[role] = entry

    add("model", "model/model.xml", "<sbml/>\n")
    add("source_expression", "inputs/source.csv", "gene,expression\nb4025,8\n")
    add("target_expression", "inputs/target.csv", "gene,expression\nb4025,1\n")
    add(
        "preflight",
        "01_preflight/preflight.csv",
        "check,status,value,message\nsolver,pass,miqp,\n",
    )
    add(
        "source_reference_fluxes",
        "02_reference/source_reference_fluxes.csv",
        "reaction_id,flux\nPGI,4.86\nPFK,7.48\n",
    )
    add(
        "gene_differential_expression",
        "03_direction/gene_differential_expression.csv",
        "gene_id,log2_fold_change,t_statistic,p_value,significant,direction\nb4025,-3.0,-51.4,1e-6,True,-1\n",
    )
    add(
        "reaction_direction_map",
        "03_direction/reaction_direction_map.csv",
        "reaction_id,direction\nPGI,-1\nPFK,0\n",
    )
    add(
        "transformation_candidates",
        "04_candidates/candidates.csv",
        "target_id\nPGI\nPFK\nTPI\nENO\n",
    )
    add(
        "transformation_ranking",
        "05_transformation/transformation_ranking.csv",
        RANKING,
    )
    add("moma_baseline", "06_validation/moma_baseline.csv", BASELINE)
    add("epsilon_sensitivity", "06_validation/epsilon_sensitivity.csv", SWEEP)

    summary = {
        "method": "rmta",
        "n_candidates": 4,
        "top_target": "PGI",
        "top_score": 3.5934,
        "reference_method": "eflux2",
        "candidate_construction": {"n_open": 4, "n_blocked_removed": 0},
        "direction_construction": {"significance": "ttest", "n_reactions_changed": 20},
    }
    provenance = {
        "model_path": "model/model.xml",
        "model_id": "fixture",
        "solver": "gurobi",
        "method": "rmta",
        "perturbation": "reaction",
        "reference_method": "eflux2",
        "source_expression": "source.csv",
        "target_expression": "target.csv",
        "alpha": 0.66,
        "epsilon": 0.01,
        "reference_state_deviation": (
            "v_ref is a deterministic eflux2 solve; Yizhak et al. (2013) use iMAT, "
            "which CMM does not implement"
        ),
    }
    add("summary", "00_summary.json", json.dumps(summary, indent=2))
    add("provenance", "00_provenance.json", json.dumps(provenance, indent=2))
    add(
        "workflow_configuration",
        "00_config.json",
        json.dumps(
            {
                "method": "rmta",
                "model_path": "model/model.xml",
                "source_expression_path": "inputs/source.csv",
                "target_expression_path": "inputs/target.csv",
            },
            indent=2,
        ),
    )
    add(
        "reproduction_config",
        "scripts/transformation_config.json",
        json.dumps(
            {
                "model_path": "../model/model.xml",
                "source_expression_path": "../inputs/source.csv",
                "target_expression_path": "../inputs/target.csv",
            }
        ),
    )
    add("reproduce_script", "scripts/reproduce.py", "# renderer fixture\n")

    (root / "00_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "workflow": "transformation_target_discovery",
                "status": "complete",
                "artifacts": artifacts,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return root


def test_a_bundle_validates_before_it_has_been_rendered(bundle: Path) -> None:
    """An analysis-only run is valid; it just says it has not been rendered."""

    report = validate_transformation_run(bundle)
    assert report.valid, report.issues
    assert report.phase == "pre-render"
    assert any("has not been rendered" in warning for warning in report.warnings)


@pytest.mark.parametrize("role", ["source_expression", "target_expression"])
def test_bundle_validation_requires_archived_expression(
    bundle: Path, role: str
) -> None:
    manifest_path = bundle / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    entry = manifest["artifacts"][role]
    path = bundle / entry["path"]
    original = path.read_bytes()
    path.write_bytes(original + b"changed,123\n")
    assert any(
        "sha256" in issue for issue in validate_transformation_run(bundle).issues
    )
    path.unlink()
    assert not validate_transformation_run(bundle).valid
    del manifest["artifacts"][role]
    manifest_path.write_text(json.dumps(manifest))
    assert any(role in issue for issue in validate_transformation_run(bundle).issues)


def test_bundle_validation_refuses_external_reproduction_inputs(bundle: Path) -> None:
    relative = "scripts/transformation_config.json"
    config = json.loads((bundle / relative).read_text())
    config["source_expression_path"] = str((bundle / "inputs/source.csv").resolve())
    _rewrite_artifact(bundle, relative, json.dumps(config))
    validation = validate_transformation_run(bundle)
    assert not validation.valid
    assert any(
        "reproduction_config.source_expression_path" in issue
        for issue in validation.issues
    )


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_r_renderer_writes_vector_and_raster_for_every_panel(bundle: Path) -> None:
    report = render_transformation_report(bundle, highlight="PGI")

    manifest = json.loads(report.figure_manifest.read_text(encoding="utf-8"))
    rendered = {
        figure["id"] for figure in manifest["figures"] if figure["status"] == "rendered"
    }
    assert rendered == {
        "fig01_component_scores",
        "fig02_transformation_ranking",
        "fig03_ranking_vs_moma",
        "fig04_epsilon_sensitivity",
    }
    assert manifest["renderer"]["engine"] == "R/ggplot2"
    # The renderer records what produced the artwork, so a figure can be traced to a script.
    assert (
        manifest["renderer"]["script_sha256"]
        == hashlib.sha256(renderer_script_path().read_bytes()).hexdigest()
    )

    for figure in manifest["figures"]:
        for suffix in ("png", "pdf", "svg"):
            output = bundle / figure["outputs"][suffix]
            assert output.is_file() and output.stat().st_size > 0
        assert figure["dpi"] == 300

    linked = report.report_html.read_text(encoding="utf-8")
    standalone = report.report_standalone_html.read_text(encoding="utf-8")
    # The two variants exist precisely because one of them survives being moved.
    assert "src='figures/fig01_component_scores.png'" in linked
    assert "src='figures/" not in standalone
    assert standalone.count("data:image/png;base64,") == 4


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
@pytest.mark.parametrize("failed_count", [1, 4])
def test_moma_comparison_preserves_but_does_not_plot_failed_solves(
    bundle: Path, failed_count: int
) -> None:
    relative = "06_validation/moma_baseline.csv"
    rows = BASELINE.splitlines()[1:]
    content = "target_id,moma_score,rank,status\n"
    for index, row in enumerate(rows):
        target_id, score, rank = row.split(",")
        failed = index >= len(rows) - failed_count
        content += (
            f"{target_id},{'-inf' if failed else score},{rank},"
            f"{'infeasible' if failed else 'optimal'}\n"
        )
    _rewrite_artifact(bundle, relative, content)

    report = render_transformation_report(bundle)
    manifest = json.loads(report.figure_manifest.read_text())
    panel = next(
        figure
        for figure in manifest["figures"]
        if figure["id"] == "fig03_ranking_vs_moma"
    )
    if failed_count == len(rows):
        assert panel["status"] == "skipped"
        assert "no candidate has a successful MOMA solve" in panel["reason"]
    else:
        assert panel["status"] == "rendered"
        assert "Ranks for 3 candidates" in panel["caption"]
        assert "1 unsuccessful MOMA solves" in panel["caption"]
    html = report.report_html.read_text()
    assert "Unsuccessful MOMA solves remain in the CSV" in html
    assert "infeasible" in html
    assert (bundle / relative).read_text() == content
    assert validate_transformation_run(bundle).valid


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_component_score_panel_is_skipped_with_a_reason_under_mta(bundle: Path) -> None:
    # MTA computes one transformation score and reports it as all three components, so the
    # best-versus-worst contrast this panel exists to show does not exist for that method. It
    # must record why rather than fail the run, and rather than draw a diagonal of identical
    # values that a reader would take for a result.
    relative = "05_transformation/transformation_ranking.csv"
    rows = (bundle / relative).read_text(encoding="utf-8").splitlines()
    rewritten = [rows[0]]
    for row in rows[1:]:
        target, score, rank, b_ts, m_ts, _w_ts = row.split(",")
        rewritten.append(",".join([target, score, rank, b_ts, m_ts, b_ts]))
    _rewrite_artifact(bundle, relative, "\n".join(rewritten) + "\n")

    report = render_transformation_report(bundle)
    manifest = json.loads(report.figure_manifest.read_text(encoding="utf-8"))
    panels = {figure["id"]: figure for figure in manifest["figures"]}
    assert panels["fig01_component_scores"]["status"] == "skipped"
    assert "rMTA" in panels["fig01_component_scores"]["reason"]
    # The run is still complete: the required panel is the one every method can draw.
    assert panels["fig02_transformation_ranking"]["status"] == "rendered"
    assert validate_transformation_run(bundle).valid

    assert validate_transformation_run(bundle).valid


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_optional_panels_are_absent_with_a_stated_reason(
    bundle: Path, tmp_path: Path
) -> None:
    """A stage that was switched off must be legible as a choice, not read as a failure."""

    manifest_path = bundle / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for role in ("moma_baseline", "epsilon_sensitivity"):
        (bundle / str(manifest["artifacts"][role]["path"])).write_text(
            "", encoding="utf-8"
        )
        manifest["artifacts"][role].update(
            {
                "status": "skipped",
                "reason": "the stage was disabled",
                "sha256": None,
                "size_bytes": 0,
            }
        )
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    report = render_transformation_report(bundle)
    figures = json.loads(report.figure_manifest.read_text(encoding="utf-8"))["figures"]
    unavailable = {f["id"]: f for f in figures if f["status"] != "rendered"}
    assert set(unavailable) == {"fig03_ranking_vs_moma", "fig04_epsilon_sensitivity"}
    for figure in unavailable.values():
        assert figure["reason"].strip()

    page = report.report_html.read_text(encoding="utf-8")
    assert "Not run" in page
    validation = validate_transformation_run(bundle)
    assert validation.valid, validation.issues


@pytest.mark.skipif(
    not _r_is_ready(), reason="Rscript renderer packages are not installed"
)
def test_a_required_panel_that_cannot_be_drawn_fails_the_render(bundle: Path) -> None:
    """The ranking panel is not optional, so an unusable ranking is an error, not a gap."""

    manifest_path = bundle / "00_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    path = bundle / str(manifest["artifacts"]["transformation_ranking"]["path"])
    path.write_text("target_id,score,rank\n", encoding="utf-8")
    entry = manifest["artifacts"]["transformation_ranking"]
    entry["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    entry["size_bytes"] = path.stat().st_size
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )

    with pytest.raises(TransformationReportError):
        render_transformation_report(bundle)


def test_render_refuses_a_directory_that_is_not_a_transformation_run(
    tmp_path: Path,
) -> None:
    with pytest.raises(TransformationReportError, match="no 00_manifest.json"):
        render_transformation_report(tmp_path)
