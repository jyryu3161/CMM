"""Write a JEV run as a self-contained, inspectable directory.

Same shape as CMM's other workflows, and built with the same writer, so a JEV run can be
opened, diffed and archived the way an SC-01 or SC-02 run can. What is specific to this
workflow is the transcript: because the agent's decisions are not guaranteed to repeat, the
exact request and response behind every move is saved. That is the only thing that makes a
single run auditable, and it is why ``04_agent/transcript.jsonl`` is a required artifact
rather than a debugging convenience.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, replace
import json
from pathlib import Path
import shutil

import pandas as pd
from cobra import Model

from cmm.core.flux_state import FluxState
from cmm.jev.engine import JevConfig, JevResult, JevWorkflowError
from cmm.workflows._bundle import ArtifactRecord, _ArtifactWriter, _jsonable

SCHEMA_VERSION = 2
WORKFLOW_ID = "jev_target_design"

_STAGE_DIRECTORIES = (
    "01_wild_type",
    "02_game",
    "03_design",
    "04_agent",
    "05_baseline",
    "06_targets",
    "model",
    "figures",
)

_OWNED_ROOT_FILES = (
    "report.html",
    "00_config.json",
    "00_provenance.json",
    "00_summary.json",
    "00_manifest.json",
)


def export_run(result: JevResult, *, model: Model, reference: FluxState) -> JevResult:
    """Write ``result`` to ``result.config.output_dir`` and return it with the path attached."""

    config = result.config
    if config.output_dir is None:  # pragma: no cover - guarded by the caller
        raise JevWorkflowError("export_run needs config.output_dir")

    root = Path(config.output_dir).expanduser().resolve()
    _prepare_directory(root, overwrite=config.overwrite)
    for name in _STAGE_DIRECTORIES:
        (root / name).mkdir(parents=True, exist_ok=True)

    writer = _ArtifactWriter(root, error_type=JevWorkflowError)

    # -- the model that was actually solved ---------------------------------
    source = Path(config.model_path)
    archived = f"model/{source.stem or 'model'}.xml"
    try:
        shutil.copyfile(source, root / archived)
    except OSError:
        from cobra.io import write_sbml_model

        write_sbml_model(model, str(root / archived))
    writer.existing(
        archived, stage="model", role="model", media_type="application/sbml+xml"
    )

    # -- 01 wild type -------------------------------------------------------
    writer.csv(
        "01_wild_type/reference_fluxes.csv",
        pd.DataFrame(sorted(reference.fluxes.items()), columns=["reaction_id", "flux"]),
        stage="01_wild_type",
        role="wild_type_reference_fluxes",
        method="pfba",
    )
    writer.json(
        "01_wild_type/wild_type_summary.json",
        {
            "product": config.product,
            "product_flux": result.wild_type_product_flux,
            "growth": result.wild_type_growth,
            "theoretical_max_yield": result.theoretical_max_yield,
            "growth_floor": config.growth_floor,
        },
        stage="01_wild_type",
        role="wild_type_summary",
    )

    # -- 02 the game --------------------------------------------------------
    writer.csv(
        "02_game/ticks.csv",
        result.ticks_frame(),
        stage="02_game",
        role="ticks",
        method="jev_decisions",
        status="complete" if result.ticks else "skipped",
        reason=None if result.ticks else "the run produced no ticks",
    )
    writer.csv(
        "02_game/rounds.csv",
        result.rounds_frame(),
        stage="02_game",
        role="rounds",
        status="complete" if result.rounds else "skipped",
        reason=None if result.rounds else "the run produced no rounds",
    )
    writer.csv(
        "02_game/candidate_rankings.csv",
        _ranking_frame(result),
        stage="02_game",
        role="candidate_rankings",
        method="jev_choice_probabilities",
        status="complete" if result.ticks else "skipped",
        reason=None if result.ticks else "no decisions were made",
    )

    # -- 03 the design ------------------------------------------------------
    writer.csv(
        "03_design/best_design.csv",
        result.interventions_frame(),
        stage="03_design",
        role="best_design",
        status="complete" if result.best_interventions else "skipped",
        reason=(
            None
            if result.best_interventions
            else "no intervention beat the wild type while holding the growth floor"
        ),
    )
    writer.csv(
        "03_design/final_design.csv",
        pd.DataFrame([i.to_record() for i in result.final_interventions]),
        stage="03_design",
        role="final_design",
        status="complete" if result.final_interventions else "skipped",
        reason=(
            None
            if result.final_interventions
            else "the run ended with no interventions applied"
        ),
    )
    writer.csv(
        "03_design/flux_trajectory.csv",
        _trajectory_frame(result),
        stage="03_design",
        role="flux_trajectory",
        method="pfba",
        status="complete" if result.flux_frames else "skipped",
        reason=None if result.flux_frames else "no flux frames were recorded",
    )

    # -- 04 the agent -------------------------------------------------------
    transcript_path = root / "04_agent/transcript.jsonl"
    transcript_path.write_text(
        "".join(
            json.dumps(_jsonable(entry), ensure_ascii=False, sort_keys=True) + "\n"
            for entry in result.transcript
        ),
        encoding="utf-8",
    )
    writer.existing(
        "04_agent/transcript.jsonl",
        stage="04_agent",
        role="agent_transcript",
        media_type="application/x-ndjson",
    )
    writer.csv(
        "04_agent/literature.csv",
        result.literature_frame(),
        stage="04_agent",
        role="literature_evidence",
        method="openrouter_web_plugin",
        status="complete" if result.literature_brief else "skipped",
        reason=(
            None
            if result.literature_brief
            else "web research was not enabled for this run"
        ),
    )
    writer.json(
        "04_agent/usage.json",
        {
            **dict(result.usage),
            "note": (
                "Token and cost totals as reported by OpenRouter for this run, covering "
                "both decision calls and any web research calls."
            ),
        },
        stage="04_agent",
        role="agent_usage",
    )

    # -- 05 what the deterministic methods give on the same problem ---------
    writer.csv(
        "05_baseline/comparison.csv",
        result.baselines_frame(),
        stage="05_baseline",
        role="baseline_comparison",
        method="optknock;robustknock;moma_l2;fseof",
        status="complete" if result.baselines else "skipped",
        reason=(
            None
            if result.baselines
            else "run_baseline_comparison was disabled, so the agent's design has nothing "
            "to be measured against"
        ),
    )

    # -- 06 every target the run weighed, and the case both ways -------------
    # The headline is one design; this is the rest of what the run learned. A reader whose
    # strain has to hold a higher growth rate, or who cannot delete three isozymes, needs the
    # second-best target and the reason it came second.
    reports = result.targets()
    writer.csv(
        "06_targets/targets.csv",
        result.targets_frame(),
        stage="06_targets",
        role="target_report",
        method="measured_deletion_and_knockdown_screen",
        status="complete" if reports else "skipped",
        reason=(
            None
            if reports
            else "no candidate reaction was measured, so there is nothing to report on"
        ),
    )
    writer.json(
        "06_targets/summary.json",
        _jsonable(result.targets_summary()),
        stage="06_targets",
        role="target_report_summary",
    )
    # -- root ---------------------------------------------------------------
    writer.json(
        "00_config.json",
        _config_payload(config, archived),
        stage="root",
        role="workflow_configuration",
    )
    writer.json(
        "00_provenance.json",
        {**dict(result.provenance), "model_path": archived},
        stage="root",
        role="provenance",
    )
    writer.json("00_summary.json", result.summary(), stage="root", role="summary")

    # One file a reader can open. The bundle is the record; this is the reading copy, and it
    # is self-contained on purpose so it can be sent to someone who does not have CMM.
    from cmm.jev.report import render_agent_report

    writer.text(
        "report.html",
        render_agent_report(result),
        stage="root",
        role="agent_report",
        media_type="text/html",
    )

    manifest_record = ArtifactRecord(
        path="00_manifest.json",
        stage="root",
        role="authoritative_artifact_manifest",
        media_type="application/json",
    )
    records = (*writer.records, manifest_record)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "workflow": WORKFLOW_ID,
        "status": (
            "complete"
            if all(record.status == "complete" for record in writer.records)
            else "partial"
        ),
        "artifacts": {
            record.role: {
                key: value for key, value in asdict(record).items() if value is not None
            }
            for record in records
        },
    }
    (root / "00_manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, sort_keys=True, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    return replace(result, run_directory=root)


def _prepare_directory(root: Path, *, overwrite: bool) -> None:
    """Refuse to write over someone else's analysis; clear only what this workflow owns."""

    if root.exists() and not root.is_dir():
        raise FileExistsError(f"output path is not a directory: {root}")
    if root.exists() and any(root.iterdir()):
        if not overwrite:
            raise FileExistsError(
                f"output directory is not empty: {root}; choose a new directory or set "
                "overwrite=True"
            )
        for name in _STAGE_DIRECTORIES:
            shutil.rmtree(root / name, ignore_errors=True)
        for name in _OWNED_ROOT_FILES:
            (root / name).unlink(missing_ok=True)
    root.mkdir(parents=True, exist_ok=True)


def _ranking_frame(result: JevResult) -> pd.DataFrame:
    """Every candidate's probability at every tick: the agent's own ranking, preserved.

    A ``choice`` answer carries a probability for every option, so each target decision is a
    complete ranking of the board rather than a single pick. Keeping it makes it possible to
    ask afterwards what the agent nearly chose — the runner-up is often the more interesting
    row — without re-running anything.
    """

    rows: list[dict[str, object]] = []
    for tick in result.ticks:
        for rank, (option, probability) in enumerate(tick.target_ranking, start=1):
            rows.append(
                {
                    "round": tick.round_index,
                    "tick": tick.tick_index,
                    "rank": rank,
                    "option": option,
                    "probability": probability,
                    "chosen": option == tick.target,
                }
            )
    return pd.DataFrame(
        rows, columns=["round", "tick", "rank", "option", "probability", "chosen"]
    )


def _trajectory_frame(result: JevResult) -> pd.DataFrame:
    """Product and growth flux at every frame, plus the frame index for the figures."""

    product = result.config.product
    rows = []
    for index, fluxes in enumerate(result.flux_frames):
        rows.append(
            {
                "frame": index,
                "product_flux": float(fluxes.get(product, 0.0)),
                "n_reactions": len(fluxes),
            }
        )
    return pd.DataFrame(rows, columns=["frame", "product_flux", "n_reactions"])


def _config_payload(config: JevConfig, archived_model: str) -> Mapping[str, object]:
    """The config as written, with paths rewritten so the bundle is self-describing.

    ``_jsonable`` is declared to return ``object`` because it walks arbitrary values; over a
    dataclass it always returns a mapping, so the narrowing is asserted rather than assumed.
    """

    encoded = _jsonable(asdict(config))
    if not isinstance(
        encoded, Mapping
    ):  # pragma: no cover - asdict always gives a mapping
        raise JevWorkflowError("the workflow config did not encode to a JSON object")
    payload = dict(encoded)
    payload["model_path"] = archived_model
    payload["output_dir"] = "."
    return payload
