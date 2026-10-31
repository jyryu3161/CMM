from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import shutil
from types import SimpleNamespace
from typing import Literal
import warnings

import pandas as pd
import pytest
from cobra import Metabolite

from cmm.core import (
    Condition,
    FluxRange,
    FluxSolution,
    FluxState,
    FvaResult,
    Medium,
    SolverCapabilityError,
)
from cmm.features import (
    BatchComparisonResult,
    BatchComparisonRow,
    FluxResponseResult,
    FseofResult,
    FvseofResult,
    Perturbation,
    ProductionEnvelope,
    ProductionYield,
    ResponseLimit,
    ResponsePoint,
    SamplingResult,
    StrainDesign,
    StrainDesignResult,
)
from cmm.features.production import EnvelopePoint
from cmm.reporting import validate_run
from cmm.workflows import production
from cmm.workflows.production import (
    AmplificationLoopDiagnosticRecord,
    AmplificationLoopDiagnosticResult,
    FluxResponseValidation,
    ProductionWorkflowConfig,
    SamplingConfig,
    SamplingValidation,
    SingleKnockoutRecord,
    ValidationConfig,
    ValidationTarget,
    run_production_target_discovery,
)


def _screen_row(
    target_id: str,
    *,
    product: float,
    growth: float = 0.8,
) -> BatchComparisonRow:
    return BatchComparisonRow(
        target_id=target_id,
        kind="gene",
        status="optimal",
        objective_value=1.0,
        distance=1.0,
        distance_kind="euclidean_l2",
        n_changed_reactions=None,
        objective=growth,
        n_reactions=1,
        product_flux=product,
    )


def _candidate_record(
    target_id: str,
    signature: str,
    *,
    method: Literal["moma_l2", "room"],
    display_rank: int | None,
    selected: bool = False,
) -> SingleKnockoutRecord:
    return SingleKnockoutRecord(
        method=method,
        target_id=target_id,
        blocked_reactions=(signature,),
        blocked_reaction_signature=signature,
        status="optimal",
        growth_rate=0.8,
        growth_fraction=0.8,
        target_production=1.0,
        product_delta=0.0,
        product_fold_change=None,
        objective_value=0.8,
        distance=1.0,
        distance_kind="euclidean_l2",
        n_changed_reactions=None,
        selected=selected,
        method_rank=1 if selected else None,
        display_rank=display_rank,
        improves_product=selected,
    )


def test_config_from_mapping_parses_nested_scientific_condition() -> None:
    config = ProductionWorkflowConfig.from_mapping(
        {
            "model_path": "model.xml",
            "product": "EX_product_e",
            "solver": "gurobi",
            "strain_design_seed": 37,
            "medium": {
                "mode": "explicit",
                "name": "defined",
                "uptake": {"EX_substrate_e": 8.0},
                "required": ["EX_substrate_e"],
            },
            "condition": {
                "name": "anaerobic",
                "bounds": [
                    {
                        "reaction_id": "EX_o2_e",
                        "lower_bound": 0.0,
                        "upper_bound": 1000.0,
                    }
                ],
                "objective": {
                    "coefficients": {"BIOMASS": 1.0},
                    "direction": "max",
                },
            },
            "validation": {
                "flux_response_biomass_fraction": 0.3,
                "sampling_growth_fraction": 0.1,
                "sampling": {"enabled": True, "n": 20, "method": "achr"},
            },
        }
    )

    assert config.solver == "gurobi"
    assert config.strain_design_seed == 37
    assert isinstance(config.medium, Medium)
    assert config.medium.required == frozenset({"EX_substrate_e"})
    assert isinstance(config.condition, Condition)
    assert config.condition.bounds[0].reaction_id == "EX_o2_e"
    assert config.condition.objective is not None
    assert config.validation.sampling.method == "achr"


def test_config_from_json_resolves_paths_from_config_directory(tmp_path) -> None:
    config_dir = tmp_path / "nested"
    config_dir.mkdir()
    config_path = config_dir / "workflow.json"
    config_path.write_text(
        json.dumps(
            {
                "model_path": "models/source.xml",
                "output_dir": "results/run-1",
                "product": "EX_product_e",
                "strain_design_seed": 41,
            }
        ),
        encoding="utf-8",
    )

    config = ProductionWorkflowConfig.from_json(config_path)

    assert config.model_path == (config_dir / "models/source.xml").resolve()
    assert config.output_dir == (config_dir / "results/run-1").resolve()
    assert config.strain_design_seed == 41


def test_default_shortlist_capacity_covers_independent_methods_and_kos() -> None:
    config = ProductionWorkflowConfig(model_path="model.xml", product="PRODUCT")

    assert config.top_amplification_targets_per_method == 10
    assert config.amplification_loop_diagnostic_top_n == 20
    assert config.validation.max_flux_response_targets == 30
    assert config.strain_design_seed == 0

    with pytest.raises(ValueError, match="strain_design_seed"):
        ProductionWorkflowConfig(
            model_path="model.xml",
            product="PRODUCT",
            strain_design_seed=2_000_000_001,
        )

    with pytest.raises(ValueError, match="independent FSEOF and FVSEOF"):
        ProductionWorkflowConfig(
            model_path="model.xml",
            product="PRODUCT",
            amplification_loop_diagnostic_top_n=19,
        )


def test_config_bounds_the_strain_design_search_and_round_trips_the_bound(
    tmp_path,
) -> None:
    """OptKnock/RobustKnock can stall forever inside straindesign; the default must end them."""

    config = ProductionWorkflowConfig(model_path="model.xml", product="PRODUCT")
    assert config.strain_design_max_seconds == 3600.0

    unbounded = ProductionWorkflowConfig(
        model_path="model.xml", product="PRODUCT", strain_design_max_seconds=None
    )
    assert unbounded.strain_design_max_seconds is None

    for bad in (0, -5.0, float("nan"), float("inf"), True, "600"):
        with pytest.raises(ValueError, match="strain_design_max_seconds"):
            ProductionWorkflowConfig(
                model_path="model.xml",
                product="PRODUCT",
                strain_design_max_seconds=bad,  # type: ignore[arg-type]
            )

    config_path = tmp_path / "workflow.json"
    for written in (90.0, None):
        config_path.write_text(
            json.dumps(
                {
                    "model_path": "model.xml",
                    "product": "PRODUCT",
                    "strain_design_max_seconds": written,
                }
            ),
            encoding="utf-8",
        )
        assert (
            ProductionWorkflowConfig.from_json(config_path).strain_design_max_seconds
            == written
        )
    with pytest.raises(ValueError, match="display-ranked single-KO candidate"):
        ProductionWorkflowConfig(
            model_path="model.xml",
            product="PRODUCT",
            validation=ValidationConfig(max_flux_response_targets=29),
        )

    # Loop eligibility controls recommendations, not whether top-ranked amplification
    # candidates receive a response scan, so disabling the diagnostic does not reduce the
    # required response capacity.
    with pytest.raises(ValueError, match="display-ranked single-KO candidate"):
        ProductionWorkflowConfig(
            model_path="model.xml",
            product="PRODUCT",
            run_amplification_loop_diagnostic=False,
            validation=ValidationConfig(max_flux_response_targets=10),
        )


def test_single_knockout_ranking_is_beneficial_unique_and_deterministic() -> None:
    perturbations = (
        Perturbation("g_a", "gene", ("R_same",)),
        Perturbation("g_b", "gene", ("R_same",)),
        Perturbation("g_c", "gene", ("R_other",)),
        Perturbation("g_epsilon", "gene", ("R_epsilon",)),
        Perturbation("g_worse", "gene", ("R_worse",)),
    )
    screen = BatchComparisonResult(
        [
            _screen_row("g_a", product=2.0),
            _screen_row("g_b", product=3.0),
            _screen_row("g_c", product=2.5),
            _screen_row("g_epsilon", product=1.001),
            _screen_row("g_worse", product=0.9),
        ]
    )

    ranked = production._rank_single_knockout_screen(
        screen,
        perturbations,
        method="room",
        wild_type_growth=1.0,
        wild_type_product=1.0,
        limit=5,
        viability_fraction=0.1,
        improvement_threshold=1e-3,
    )

    selected = sorted(
        (item for item in ranked if item.selected),
        key=lambda item: item.method_rank or 0,
    )
    assert [(item.target_id, item.method_rank) for item in selected] == [
        ("g_b", 1),
        ("g_c", 2),
    ]
    assert not next(item for item in ranked if item.target_id == "g_a").selected
    assert not next(item for item in ranked if item.target_id == "g_epsilon").selected
    assert all(item.improves_product for item in selected)
    displayed = sorted(
        (item for item in ranked if item.display_rank is not None),
        key=lambda item: item.display_rank or 0,
    )
    assert [(item.target_id, item.display_rank) for item in displayed] == [
        ("g_b", 1),
        ("g_c", 2),
        ("g_epsilon", 3),
        ("g_worse", 4),
    ]


def test_room_ranking_uses_epsilon_plus_numeric_margin_strictly() -> None:
    room_epsilon = 1e-3
    numeric_margin = 1e-6
    threshold = room_epsilon + numeric_margin
    perturbations = (
        Perturbation("g_at_threshold", "gene", ("R_at_threshold",)),
        Perturbation("g_above_threshold", "gene", ("R_above_threshold",)),
    )
    screen = BatchComparisonResult(
        [
            _screen_row("g_at_threshold", product=threshold),
            _screen_row(
                "g_above_threshold",
                product=math.nextafter(threshold, math.inf),
            ),
        ],
        metadata={"comparison_method": "room", "epsilon": room_epsilon},
    )

    ranked = production._rank_single_knockout_screen(
        screen,
        perturbations,
        method="room",
        wild_type_growth=1.0,
        wild_type_product=0.0,
        limit=5,
        viability_fraction=0.1,
        improvement_threshold=threshold,
    )
    by_target = {item.target_id: item for item in ranked}

    assert by_target["g_at_threshold"].product_delta == threshold
    assert by_target["g_at_threshold"].display_rank == 2
    assert not by_target["g_at_threshold"].improves_product
    assert not by_target["g_at_threshold"].selected
    assert by_target["g_above_threshold"].product_delta > threshold
    assert by_target["g_above_threshold"].display_rank == 1
    assert by_target["g_above_threshold"].improves_product
    assert by_target["g_above_threshold"].selected


def test_preflight_refuses_l2_moma_downgrade(toy_model) -> None:
    toy_model.solver = "glpk"
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_strain_design=False,
        run_amplification=False,
        validation=ValidationConfig(enabled=False),
    )

    with pytest.raises(SolverCapabilityError, match="L2 MOMA") as raised:
        production._preflight_model(toy_model, config)

    assert raised.value.capability == "QP"


def test_flux_response_maps_candidate_reaction_to_product_for_both_candidate_types(
    toy_model, monkeypatch
) -> None:
    calls: list[
        tuple[
            str,
            str | None,
            tuple[float, float],
            float | None,
            float | None,
            float | None,
        ]
    ] = []

    def fake_flux_response(model, target, response=None, **kwargs):
        calls.append(
            (
                target,
                response,
                tuple(model.reactions.get_by_id("SOURCE_A").bounds),
                kwargs.get("biomass_fraction"),
                kwargs.get("target_min"),
                kwargs.get("target_max"),
            )
        )
        return object()

    monkeypatch.setattr(production, "flux_response", fake_flux_response)
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_single_knockout=False,
        run_strain_design=False,
        run_amplification=False,
        validation=ValidationConfig(
            sampling=SamplingConfig(enabled=False),
            flux_response_biomass_fraction=0.3,
        ),
    )
    targets = (
        ValidationTarget(
            target_id="SOURCE_A",
            scan_reaction="SOURCE_A",
            response_reaction="PRODUCT",
            background="wild_type",
            actions=("amplify",),
            source_methods=("fseof", "fvseof"),
        ),
        ValidationTarget(
            target_id="gene_x",
            scan_reaction="SOURCE_A",
            response_reaction="PRODUCT",
            background="wild_type",
            actions=("knockout",),
            source_methods=("moma_l2", "room"),
            blocked_reactions=("SOURCE_A",),
            candidate_scope="all_display_ranked_candidates",
        ),
    )

    results = production._run_flux_response_validation(
        toy_model,
        config,
        "BIOMASS",
        targets,
        FluxState({"SOURCE_A": 4.0}, name="reference"),
    )

    assert all(item.status == "complete" for item in results)
    assert calls == [
        ("SOURCE_A", "PRODUCT", (0.0, 10.0), 0.3, None, None),
        ("SOURCE_A", "PRODUCT", (0.0, 10.0), 0.3, 0.0, 4.0),
    ]
    assert results[1].target.scan_reference_flux == pytest.approx(4.0)
    assert toy_model.reactions.get_by_id("SOURCE_A").bounds == (0.0, 10.0)


def test_flux_response_keeps_multi_reaction_knockout_as_explicit_skip(
    toy_model, monkeypatch
) -> None:
    monkeypatch.setattr(
        production,
        "flux_response",
        lambda *args, **kwargs: pytest.fail(
            "multi-reaction target must not be scanned"
        ),
    )
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_single_knockout=False,
        run_strain_design=False,
        run_amplification=False,
        validation=ValidationConfig(sampling=SamplingConfig(enabled=False)),
    )
    target = ValidationTarget(
        target_id="gene_complex",
        scan_reaction="R1;R2",
        response_reaction="PRODUCT",
        background="wild_type",
        actions=("knockout",),
        source_methods=("moma_l2",),
        blocked_reactions=("R1", "R2"),
        candidate_scope="all_display_ranked_candidates",
        blocked_reaction_signature="R1;R2",
    )

    result = production._run_flux_response_validation(
        toy_model,
        config,
        "BIOMASS",
        (target,),
        FluxState({"R1": 1.0, "R2": 1.0}, name="reference"),
    )[0]

    assert result.status == "skipped"
    assert result.result is None
    assert result.reason is not None and "multiple reactions" in result.reason


def test_zero_reference_knockout_still_receives_full_domain_response(
    toy_model, monkeypatch
) -> None:
    calls: list[dict[str, object]] = []

    def fake_response(model, target, response=None, **kwargs):
        del model, target, response
        calls.append(kwargs)
        return object()

    monkeypatch.setattr(production, "flux_response", fake_response)
    config = ProductionWorkflowConfig(model_path="model.xml", product="PRODUCT")
    target = ValidationTarget(
        target_id="gene_zero",
        scan_reaction="SOURCE_A",
        response_reaction="PRODUCT",
        background="wild_type",
        actions=("knockout",),
        source_methods=("room",),
        blocked_reactions=("SOURCE_A",),
        candidate_scope="all_display_ranked_candidates",
    )

    result = production._run_flux_response_validation(
        toy_model,
        config,
        "BIOMASS",
        (target,),
        FluxState({"SOURCE_A": 0.0}, name="reference"),
    )[0]

    assert result.status == "complete"
    assert result.target.scan_reference_flux == pytest.approx(0.0)
    assert calls[0]["target_min"] is None
    assert calls[0]["target_max"] is None


def test_sampling_comparison_reports_wasserstein_distance() -> None:
    wild_type = SamplingResult(
        pd.DataFrame({"R1": [0.0, 1.0, 2.0], "R2": [2.0, 2.0, 2.0]}),
        method="achr",
        seed=0,
    )
    knockout = SamplingResult(
        pd.DataFrame({"R1": [2.0, 3.0, 4.0], "R2": [2.0, 2.0, 2.0]}),
        method="achr",
        seed=0,
    )

    comparison = production._sampling_comparison(wild_type, knockout).set_index(
        "reaction_id"
    )

    assert comparison.loc["R1", "wasserstein_distance"] == pytest.approx(2.0)
    assert comparison.loc["R2", "wasserstein_distance"] == pytest.approx(0.0)


def test_gene_knockout_mapping_preserves_all_genes_names_and_gprs(
    branched_model,
) -> None:
    branched_model.reactions.R1.name = "entry reaction"
    branched_model.genes.g1.name = "entry enzyme"
    branched_model.reactions.R2.gene_reaction_rule = "g2 or g2_iso"
    branched_model.genes.g2.name = "short branch enzyme"
    branched_model.genes.g2_iso.name = "short branch isozyme"

    mapping = production._gene_knockout_mapping(branched_model).to_frame()

    assert set(mapping["gene_id"]) == {gene.id for gene in branched_model.genes}
    g1 = mapping[mapping["gene_id"] == "g1"].iloc[0]
    assert not bool(g1["inert"])
    assert g1["gene_name"] == "entry enzyme"
    assert g1["blocked_reaction"] == "R1"
    assert g1["reaction_name"] == "entry reaction"
    assert "A" in str(g1["reaction_equation"])
    assert "B" in str(g1["reaction_equation"])
    assert g1["gpr"] == "g1"
    isozyme = mapping[mapping["gene_id"] == "g2_iso"].iloc[0]
    assert bool(isozyme["inert"])
    assert pd.isna(isozyme["blocked_reaction"])
    assert "R2: g2 or g2_iso" in str(isozyme["gpr"])


def test_amplification_loop_diagnostic_flags_declared_capacity_ratio(
    toy_model, monkeypatch
) -> None:
    toy_model.reactions.PRODUCT.lower_bound = 1.5
    calls: list[tuple[object, float]] = []

    def fake_fva(model, reactions=None, **kwargs):
        calls.append((kwargs["loopless"], model.reactions.PRODUCT.lower_bound))
        target = tuple(reactions or ())[0]
        flux_range = (
            FluxRange(-50.0, 50.0)
            if kwargs["loopless"] is False
            else FluxRange(-2.5, 2.5)
        )
        return FvaResult(
            {target: flux_range}, metadata={"loopless": kwargs["loopless"]}
        )

    monkeypatch.setattr(production, "fva", fake_fva)
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_single_knockout=False,
        run_strain_design=False,
        run_amplification=False,
        loopless_capacity_ratio_threshold=0.1,
        validation=ValidationConfig(enabled=False),
    )

    diagnostic = production._run_amplification_loop_diagnostic(
        toy_model,
        config,
        "BIOMASS",
        0.8,
        (("SOURCE_A", ("fseof", "fvseof")),),
        enforced_product_floor=0.0,
    )

    assert calls == [(False, 1.5), ("fastSNP", 1.5)]
    assert len(diagnostic.records) == 1
    record = diagnostic.records[0]
    assert record.loopless_to_standard_capacity_ratio == pytest.approx(0.05)
    assert record.capacity_ratio_threshold == pytest.approx(0.1)
    assert record.loop_artifact_flag is True
    assert record.enforced_product_floor == pytest.approx(1.5)
    assert diagnostic.metadata["parameters"][
        "requested_enforced_product_floor"
    ] == pytest.approx(0.0)
    assert diagnostic.metadata["enforced_product_floor"] == pytest.approx(1.5)


@pytest.mark.parametrize(
    ("standard", "loopless", "reason_fragment"),
    [
        (FluxRange(0.0, 10.0), FluxRange(float("nan"), float("nan")), "loopless"),
        (FluxRange(0.0, 10.0), FluxRange(float("-inf"), float("inf")), "loopless"),
        (FluxRange(float("-inf"), float("inf")), FluxRange(0.0, 1.0), "standard"),
    ],
)
def test_loop_diagnostic_never_clears_non_finite_fva_bounds(
    standard: FluxRange,
    loopless: FluxRange,
    reason_fragment: str,
) -> None:
    record = production._amplification_loop_record(
        1,
        "TARGET",
        ("fseof",),
        standard,
        loopless,
        threshold=0.1,
        product_floor=1.0,
        biomass_floor=0.3,
    )

    assert record.diagnostic_status == "inconclusive"
    assert record.loop_artifact_flag is None
    assert record.loopless_to_standard_capacity_ratio is None
    assert record.reason is not None and reason_fragment in record.reason
    assert all(
        value is None or math.isfinite(value)
        for value in (
            record.standard_minimum,
            record.standard_maximum,
            record.standard_capacity,
            record.loopless_minimum,
            record.loopless_maximum,
            record.loopless_capacity,
        )
    )


def test_artifact_slugs_are_bounded_and_collision_resistant() -> None:
    values = ["A/B", "A?B", "A B", "A_B", "abc", "ABC", "x" * 300]
    slugs = [production._slug(value) for value in values]

    assert len(slugs) == len(set(slugs))
    assert all(len(slug) <= 96 for slug in slugs)
    assert all("/" not in slug and "\\" not in slug for slug in slugs)
    assert production._slug("A/B") == production._slug("A/B")


def test_recommendations_require_validation_and_exclude_loop_artifacts() -> None:
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_single_knockout=False,
        run_strain_design=False,
        run_amplification=True,
        validation=ValidationConfig(
            sampling=SamplingConfig(enabled=True, n=5, method="achr")
        ),
    )
    response = FluxResponseResult(
        target="AMP",
        response="PRODUCT",
        biomass="BIOMASS",
        points=(
            ResponsePoint(1.0, 0.0, 1.0, "optimal"),
            ResponsePoint(2.0, 3.0, 0.5, "optimal"),
        ),
        phases=(),
        limit=ResponseLimit(False, "no response limit"),
        wild_type={
            "target_flux": 1.0,
            "response_flux": 0.0,
            "biomass_flux": 1.0,
        },
    )
    validation_target = ValidationTarget(
        target_id="AMP",
        scan_reaction="AMP",
        response_reaction="PRODUCT",
        background="wild_type",
        actions=("amplify",),
        source_methods=("fseof", "fvseof"),
    )

    class CandidateResult:
        def amplification_targets(self):
            return ["AMP"]

    def workflow_result(flag: bool):
        return SimpleNamespace(
            config=config,
            reference=FluxState(
                {"BIOMASS": 1.0, "PRODUCT": 0.0},
                name="reference",
            ),
            biomass="BIOMASS",
            selected_single_knockouts=(),
            flux_responses=(
                FluxResponseValidation(
                    target=validation_target,
                    status="complete",
                    result=response,
                ),
            ),
            sampling=(),
            fseof_result=CandidateResult(),
            fvseof_result=CandidateResult(),
            amplification_loop_diagnostic=AmplificationLoopDiagnosticResult(
                (
                    AmplificationLoopDiagnosticRecord(
                        rank=1,
                        target="AMP",
                        source_methods=("fseof", "fvseof"),
                        standard_minimum=0.0,
                        standard_maximum=10.0,
                        standard_capacity=10.0,
                        loopless_minimum=0.0,
                        loopless_maximum=5.0,
                        loopless_capacity=5.0,
                        loopless_to_standard_capacity_ratio=0.5,
                        capacity_ratio_threshold=0.1,
                        loop_artifact_flag=flag,
                        diagnostic_status="complete",
                        reason=None,
                        enforced_product_floor=1.0,
                        biomass_floor=0.3,
                    ),
                ),
                {},
            ),
            robustknock_result=None,
        )

    flagged_result = workflow_result(True)
    assert production._recommendations_frame(flagged_result).empty
    publication_targets, status, reason = production._publication_amplification_targets(
        flagged_result, ("AMP",)
    )
    assert publication_targets == ("AMP",)
    assert status == "complete"
    assert reason is not None and "diagnostic-only" in reason
    annotated = production._annotate_loop_diagnostic(
        pd.DataFrame({"target": ["AMP"]}), flagged_result
    )
    assert bool(annotated.iloc[0]["loop_artifact_flag"])
    assert annotated.iloc[0]["loop_diagnostic_status"] == "complete"
    assert "loop_diagnostic_reason" in annotated.columns

    supported_result = workflow_result(False)
    supported_targets, status, reason = production._publication_amplification_targets(
        supported_result, ("AMP",)
    )
    assert supported_targets == ("AMP",)
    assert status == "complete"
    assert reason is None
    supported = production._recommendations_frame(supported_result)
    assert supported[["target", "type", "verdict"]].to_dict("records") == [
        {"target": "AMP", "type": "amplification", "verdict": "support"}
    ]
    assert not bool(supported.iloc[0]["artifact_flag"])


def test_amplification_recommendations_use_independent_method_shortlists() -> None:
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_single_knockout=False,
        run_strain_design=False,
    )

    class Candidates:
        def __init__(self, targets: list[str]) -> None:
            self._targets = targets

        def amplification_targets(self) -> list[str]:
            return self._targets

    def response(target: str, method: str) -> FluxResponseValidation:
        return FluxResponseValidation(
            target=ValidationTarget(
                target_id=target,
                scan_reaction=target,
                response_reaction="PRODUCT",
                background="wild_type",
                actions=("amplify",),
                source_methods=(method,),
            ),
            status="complete",
            result=FluxResponseResult(
                target,
                "PRODUCT",
                "BIOMASS",
                (
                    ResponsePoint(1.0, 0.0, 1.0, "optimal"),
                    ResponsePoint(2.0, 3.0, 0.5, "optimal"),
                ),
                (),
                ResponseLimit(False, "no limit"),
                {
                    "target_flux": 1.0,
                    "response_flux": 0.0,
                    "biomass_flux": 1.0,
                },
            ),
        )

    def diagnostic(target: str, method: str) -> AmplificationLoopDiagnosticRecord:
        return AmplificationLoopDiagnosticRecord(
            rank=1,
            target=target,
            source_methods=(method,),
            standard_minimum=0.0,
            standard_maximum=10.0,
            standard_capacity=10.0,
            loopless_minimum=0.0,
            loopless_maximum=5.0,
            loopless_capacity=5.0,
            loopless_to_standard_capacity_ratio=0.5,
            capacity_ratio_threshold=0.1,
            loop_artifact_flag=False,
            diagnostic_status="complete",
            reason=None,
            enforced_product_floor=1.0,
            biomass_floor=0.3,
        )

    result = SimpleNamespace(
        config=config,
        reference=FluxState(
            {"BIOMASS": 1.0, "PRODUCT": 0.0},
            name="reference",
        ),
        biomass="BIOMASS",
        selected_single_knockouts=(),
        sampling=(),
        fseof_result=Candidates(["F_ONLY"]),
        fvseof_result=Candidates(["V_ONLY"]),
        flux_responses=(
            response("F_ONLY", "fseof"),
            response("V_ONLY", "fvseof"),
        ),
        amplification_loop_diagnostic=AmplificationLoopDiagnosticResult(
            (
                diagnostic("F_ONLY", "fseof"),
                diagnostic("V_ONLY", "fvseof"),
            ),
            {},
        ),
        robustknock_result=None,
    )

    recommendations = production._recommendations_frame(result).set_index("target")
    assert set(recommendations.index) == {"F_ONLY", "V_ONLY"}
    assert recommendations.at["F_ONLY", "proposal_methods"] == "fseof"
    assert recommendations.at["V_ONLY", "proposal_methods"] == "fvseof"
    assert recommendations.at["F_ONLY", "evidence"].startswith("FSEOF candidate")
    assert recommendations.at["V_ONLY", "evidence"].startswith("FVSEOF candidate")
    assert (
        recommendations["reason"]
        .str.contains("cross-method agreement was not required")
        .all()
    )


def test_knockout_recommendation_requires_concordant_sampling_shifts() -> None:
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        validation=ValidationConfig(
            sampling=SamplingConfig(enabled=True, n=5, method="achr")
        ),
    )
    candidate = production.SingleKnockoutRecord(
        method="moma_l2",
        target_id="g1",
        blocked_reactions=("KO",),
        blocked_reaction_signature="KO",
        status="optimal",
        growth_rate=0.8,
        growth_fraction=0.8,
        target_production=1.0,
        product_delta=1.0,
        product_fold_change=None,
        objective_value=1.0,
        distance=1.0,
        distance_kind="euclidean_l2",
        n_changed_reactions=None,
        selected=True,
        method_rank=1,
        improves_product=True,
    )
    response_target = ValidationTarget(
        target_id="g1",
        scan_reaction="KO",
        response_reaction="PRODUCT",
        background="wild_type",
        actions=("knockout",),
        source_methods=("moma_l2",),
        blocked_reactions=("KO",),
        candidate_scope="all_display_ranked_candidates",
        scan_reference_flux=2.0,
    )
    response = FluxResponseResult(
        target="KO",
        response="PRODUCT",
        biomass="BIOMASS",
        points=(
            ResponsePoint(0.0, 1.0, 0.5, "optimal"),
            ResponsePoint(2.0, 0.2, 0.5, "optimal"),
        ),
        phases=(),
        limit=ResponseLimit(False, "no response limit"),
        wild_type={
            "target_flux": 2.0,
            "response_flux": 0.2,
            "biomass_flux": 1.0,
        },
    )

    def result_with_shifts(mean_delta: float, median_delta: float):
        comparison = pd.DataFrame(
            [
                {
                    "reaction_id": "PRODUCT",
                    "mean_delta": mean_delta,
                    "median_delta": median_delta,
                    "knockout_median": 0.0,
                },
                {
                    "reaction_id": "BIOMASS",
                    "mean_delta": 0.0,
                    "median_delta": 0.0,
                    "knockout_median": 0.6,
                },
            ]
        )
        return SimpleNamespace(
            config=config,
            reference=FluxState(
                {"BIOMASS": 1.0, "PRODUCT": 0.0},
                name="reference",
            ),
            biomass="BIOMASS",
            selected_single_knockouts=(candidate,),
            flux_responses=(
                FluxResponseValidation(
                    target=response_target,
                    status="complete",
                    result=response,
                ),
            ),
            sampling=(
                SamplingValidation(
                    target_id="g1",
                    blocked_reactions=("KO",),
                    source_methods=("moma_l2",),
                    status="complete",
                    comparison=comparison,
                ),
            ),
            fseof_result=None,
            fvseof_result=None,
            amplification_loop_diagnostic=None,
            robustknock_result=None,
        )

    # A positive median cannot promote a target when the sampled mean moves down.
    assert production._recommendations_frame(result_with_shifts(-0.1, 0.2)).empty

    supported = production._recommendations_frame(result_with_shifts(0.3, 0.2))
    assert supported["target"].tolist() == ["g1"]
    assert supported.iloc[0]["product_effect"] == pytest.approx(0.2)
    assert "zero-flux titration product delta=0.8" in supported.iloc[0]["evidence"]
    assert "sampling mean product delta=0.3" in supported.iloc[0]["evidence"]
    assert "sampling median product delta=0.2" in supported.iloc[0]["evidence"]
    assert "response product delta" not in supported.iloc[0]["evidence"]


def test_recommendations_rank_effect_within_intervention_class() -> None:
    config = ProductionWorkflowConfig(model_path="model.xml", product="PRODUCT")
    result = SimpleNamespace(
        config=config,
        reference=FluxState(
            {"BIOMASS": 1.0, "PRODUCT": 0.0},
            name="reference",
        ),
        biomass="BIOMASS",
        selected_single_knockouts=(),
        flux_responses=(),
        sampling=(),
        fseof_result=None,
        fvseof_result=None,
        amplification_loop_diagnostic=None,
        optknock_result=None,
        robustknock_result=StrainDesignResult(
            "robustknock",
            "PRODUCT",
            (
                StrainDesign(("A_LOW",), 0.8, 1.0, 1.0),
                StrainDesign(("Z_HIGH",), 0.5, 5.0, 5.0),
                StrainDesign(("B_TIE",), 0.6, 1.0, 1.0),
            ),
        ),
    )

    recommendations = production._recommendations_frame(result)
    assert recommendations["target"].tolist() == ["Z_HIGH", "A_LOW", "B_TIE"]
    assert recommendations["product_effect"].tolist() == [5.0, 1.0, 1.0]
    assert (
        recommendations["evidence"].str.startswith("Reaction-level RobustKnock").all()
    )
    assert recommendations["reason"].str.contains("GPR-resolved").all()

    result.robustknock_result = StrainDesignResult(
        "robustknock",
        "PRODUCT",
        (StrainDesign(("INCUMBENT",), 0.8, 6.0, 5.0),),
        {
            "parameters": {
                "straindesign_search_status": "time_limit",
                "straindesign_search_complete": False,
            }
        },
    )
    assert production._strain_design_search_status(result.robustknock_result) == (
        "time_limit"
    )
    assert not production._strain_design_search_complete(result.robustknock_result)
    assert production._recommendations_frame(result).empty
    status, reason = production._recommendation_artifact_status(result)
    assert status == "partial"
    assert reason is not None and "strain-design search was incomplete" in reason


def test_amplification_validation_runs_flagged_and_unresolved_with_eligibility() -> (
    None
):
    class Candidates:
        def amplification_targets(self):
            return ["CLEAR", "FLAGGED", "UNRESOLVED"]

    def diagnostic(
        rank: int,
        target: str,
        *,
        status: str,
        flag: bool | None,
    ) -> AmplificationLoopDiagnosticRecord:
        return AmplificationLoopDiagnosticRecord(
            rank=rank,
            target=target,
            source_methods=("fseof", "fvseof"),
            standard_minimum=0.0,
            standard_maximum=10.0,
            standard_capacity=10.0,
            loopless_minimum=0.0,
            loopless_maximum=5.0,
            loopless_capacity=5.0,
            loopless_to_standard_capacity_ratio=0.5,
            capacity_ratio_threshold=0.1,
            loop_artifact_flag=flag,
            diagnostic_status=status,
            reason=None,
            enforced_product_floor=1.0,
            biomass_floor=0.3,
        )

    loop_diagnostic = AmplificationLoopDiagnosticResult(
        (
            diagnostic(1, "CLEAR", status="complete", flag=False),
            diagnostic(2, "FLAGGED", status="complete", flag=True),
            diagnostic(3, "UNRESOLVED", status="inconclusive", flag=None),
        ),
        {},
    )
    candidates = Candidates()
    targets = production._validation_targets(
        (),
        candidates,
        candidates,
        loop_diagnostic,
        product="PRODUCT",
        biomass="BIOMASS",
        amplification_limit=3,
        limit=10,
    )

    assert [target.target_id for target in targets] == [
        "CLEAR",
        "FLAGGED",
        "UNRESOLVED",
    ]
    assert all(target.source_methods == ("fseof", "fvseof") for target in targets)
    assert all(
        target.candidate_scope == "all_report_selected_candidates" for target in targets
    )
    by_target = {target.target_id: target for target in targets}
    assert by_target["CLEAR"].loop_diagnostic_eligible is True
    assert by_target["FLAGGED"].loop_diagnostic_eligible is False
    assert by_target["FLAGGED"].loop_artifact_flag is True
    assert by_target["UNRESOLVED"].loop_diagnostic_eligible is False
    assert by_target["UNRESOLVED"].loop_diagnostic_status == "inconclusive"


def test_validation_universe_covers_all_display_kos_and_reported_amplifications() -> (
    None
):
    """The T124700-shaped 9 KO + 14 amplification universe must yield 23 scans."""

    single_records = tuple(
        [
            _candidate_record(
                f"gm{rank}",
                f"R{rank}",
                method="moma_l2",
                display_rank=rank + 1,
            )
            for rank in range(5)
        ]
        + [
            _candidate_record(
                f"gr{rank}",
                "R0" if rank == 0 else f"R{rank + 4}",
                method="room",
                display_rank=rank + 1,
            )
            for rank in range(5)
        ]
    )

    class Candidates:
        def __init__(self, targets: list[str]) -> None:
            self.targets = targets

        def amplification_targets(self) -> list[str]:
            return self.targets

    fseof_candidates = Candidates([f"A{rank}" for rank in range(10)])
    fvseof_candidates = Candidates(
        [f"A{rank}" for rank in range(6)] + [f"V{rank}" for rank in range(6, 10)]
    )
    targets = production._validation_targets(
        single_records,
        fseof_candidates,
        fvseof_candidates,
        None,
        product="PRODUCT",
        biomass="BIOMASS",
        amplification_limit=10,
        limit=30,
    )

    knockout = [
        target
        for target in targets
        if target.candidate_scope == "all_display_ranked_candidates"
    ]
    amplification = [
        target
        for target in targets
        if target.candidate_scope == "all_report_selected_candidates"
    ]
    assert len(knockout) == 9
    assert len(amplification) == 14
    assert len(targets) == 23
    shared_signature = next(
        target for target in knockout if target.blocked_reaction_signature == "R0"
    )
    assert shared_signature.target_id == "gm0"
    assert shared_signature.candidate_target_ids == ("gm0", "gr0")
    assert shared_signature.source_methods == ("moma_l2", "room")
    assert all(
        target.candidate_scope == "all_display_ranked_candidates" for target in knockout
    )
    assert all(
        target.candidate_scope == "all_report_selected_candidates"
        for target in amplification
    )
    assert all(target.background == "wild_type" for target in targets)
    assert all(target.response_reaction == "PRODUCT" for target in targets)
    assert all(
        target.scan_reaction == target.blocked_reactions[0] for target in knockout
    )
    assert all(
        target.loop_diagnostic_status == "unavailable" for target in amplification
    )

    with pytest.raises(
        production.ProductionWorkflowError, match="instead of truncating"
    ):
        production._validation_targets(
            single_records,
            fseof_candidates,
            fvseof_candidates,
            None,
            product="PRODUCT",
            biomass="BIOMASS",
            amplification_limit=10,
            limit=22,
        )


def test_single_knockout_alias_provenance_is_complete_and_order_stable(
    toy_model,
) -> None:
    """A SUCDi-like phenotype retains every gene alias without adding candidate rows."""

    records = tuple(
        _candidate_record(
            target_id,
            "SUCDi",
            method=method,
            display_rank=1 if target_id == "b0721" else None,
        )
        for method in ("room", "moma_l2")
        for target_id in ("b0724", "b0722", "b0721", "b0723")
    )

    def validation_target(
        ordered: tuple[SingleKnockoutRecord, ...],
    ) -> ValidationTarget:
        targets = production._validation_targets(
            ordered,
            None,
            None,
            None,
            product="PRODUCT",
            biomass="BIOMASS",
            amplification_limit=10,
            limit=10,
        )
        assert len(targets) == 1
        return targets[0]

    forward = validation_target(records)
    reverse = validation_target(tuple(reversed(records)))
    expected_aliases = ("b0721", "b0722", "b0723", "b0724")
    assert forward == reverse
    assert forward.target_id == "b0721"
    assert forward.candidate_target_ids == expected_aliases
    assert forward.source_methods == ("moma_l2", "room")

    candidate_rows = production._single_knockout_candidate_export_rows(
        records,
        toy_model,
    )
    assert [(row["method"], row["target_id"]) for row in candidate_rows] == [
        ("moma_l2", "b0721"),
        ("room", "b0721"),
    ]
    assert {row["candidate_target_ids"] for row in candidate_rows} == {
        ";".join(expected_aliases)
    }
    assert sum(bool(row["validation_representative"]) for row in candidate_rows) == 1


def test_sampling_runs_wild_type_plus_every_unique_display_candidate(
    toy_model, monkeypatch
) -> None:
    records = tuple(
        [
            _candidate_record(
                f"gm{rank}",
                f"R{rank}",
                method="moma_l2",
                display_rank=rank + 1,
            )
            for rank in range(5)
        ]
        + [
            _candidate_record(
                f"gr{rank}",
                "R0" if rank == 0 else f"R{rank + 4}",
                method="room",
                display_rank=rank + 1,
            )
            for rank in range(5)
        ]
        + [
            _candidate_record(
                "aa_moma_alias",
                "R0",
                method="moma_l2",
                display_rank=None,
            ),
            _candidate_record(
                "aa_room_alias",
                "R0",
                method="room",
                display_rank=None,
            ),
        ]
    )
    perturbations = tuple(
        Perturbation(record.target_id, "gene", ("SOURCE_A",)) for record in records
    )
    calls = 0

    def fake_sampling(model, **kwargs):
        nonlocal calls
        del model
        calls += 1
        return SamplingResult(
            pd.DataFrame(
                {
                    "SOURCE_A": [1.0, 1.0],
                    "BIOMASS": [0.8, 0.8],
                    "PRODUCT": [0.2, 0.3],
                }
            ),
            method=str(kwargs["method"]),
            seed=int(kwargs["seed"]),
        )

    monkeypatch.setattr(production, "random_flux_sampling", fake_sampling)
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_strain_design=False,
        run_amplification=False,
        validation=ValidationConfig(
            max_flux_response_targets=10,
            sampling=SamplingConfig(n=2, thinning=1),
        ),
    )

    sampling = production._run_sampling_validation(
        toy_model,
        config,
        "BIOMASS",
        1.0,
        records,
        perturbations,
    )

    assert calls == 10
    assert len(sampling) == 10
    assert sampling[0].target_id == "wild_type"
    candidates = sampling[1:]
    assert len(candidates) == 9
    assert all(item.status == "complete" for item in candidates)
    assert all(
        item.candidate_scope == "all_display_ranked_candidates" for item in candidates
    )
    shared_signature = next(
        item for item in candidates if item.blocked_reaction_signature == "R0"
    )
    assert shared_signature.target_id == "gm0"
    assert shared_signature.candidate_target_ids == (
        "aa_moma_alias",
        "aa_room_alias",
        "gm0",
        "gr0",
    )
    coverage = production._validation_coverage(
        tuple(
            ValidationTarget(
                target_id=item.target_id,
                scan_reaction=item.blocked_reactions[0],
                response_reaction="PRODUCT",
                background="wild_type",
                actions=("knockout",),
                source_methods=item.source_methods,
                candidate_scope="all_display_ranked_candidates",
            )
            for item in candidates
        ),
        (),
        sampling,
        sampling_enabled=True,
    )
    assert coverage["sampling_expected"] == 10
    assert coverage["sampling_attempted"] == 10
    assert coverage["sampling_completed"] == 10


def test_sampling_unavailable_keeps_explicit_candidate_skip_rows(
    toy_model, monkeypatch
) -> None:
    record = _candidate_record(
        "g_display",
        "SOURCE_A",
        method="room",
        display_rank=1,
    )
    monkeypatch.setattr(
        production,
        "random_flux_sampling",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("no warmup")),
    )
    config = ProductionWorkflowConfig(
        model_path="model.xml",
        product="PRODUCT",
        run_strain_design=False,
        run_amplification=False,
        validation=ValidationConfig(
            max_flux_response_targets=10,
            sampling=SamplingConfig(n=2, thinning=1),
        ),
    )

    sampling = production._run_sampling_validation(
        toy_model,
        config,
        "BIOMASS",
        1.0,
        (record,),
        (Perturbation("g_display", "gene", ("SOURCE_A",)),),
    )

    assert [item.status for item in sampling] == ["failed", "skipped"]
    assert sampling[1].candidate_scope == "all_display_ranked_candidates"
    assert sampling[1].candidate_target_ids == ("g_display",)
    assert (
        sampling[1].reason is not None
        and "wild-type sampling failed" in sampling[1].reason
    )


def test_amplification_candidate_union_is_rank_first_not_intersection_first() -> None:
    class Candidates:
        def __init__(self, targets: list[str]) -> None:
            self.targets = targets

        def amplification_targets(self) -> list[str]:
            return self.targets

    fseof_candidates = Candidates(["F_ONLY", "COMMON"])
    fvseof_candidates = Candidates(["V_ONLY", "COMMON"])

    capped = production._amplification_candidate_sources(
        fseof_candidates,
        fvseof_candidates,
        per_method_limit=2,
        total_limit=2,
    )
    assert capped == (("F_ONLY", ("fseof",)), ("V_ONLY", ("fvseof",)))

    complete = production._amplification_candidate_sources(
        fseof_candidates,
        fvseof_candidates,
        per_method_limit=2,
        total_limit=4,
    )
    assert complete == (
        ("F_ONLY", ("fseof",)),
        ("V_ONLY", ("fvseof",)),
        ("COMMON", ("fseof", "fvseof")),
    )


def test_public_runner_exports_partial_schema_v2_with_sidecars_and_scripts(
    toy_model, tmp_path, monkeypatch
) -> None:
    substrate = toy_model.reactions.SOURCE_A
    source_metabolite = next(iter(substrate.metabolites))
    substrate.add_metabolites({source_metabolite: -2.0})
    substrate.bounds = (-10.0, 1000.0)
    source_metabolite.formula = "C"
    oxygen = Metabolite("oxygen_e", formula="O2", compartment="c")
    toy_model.add_boundary(
        oxygen,
        type="exchange",
        reaction_id="EX_oxygen_e",
        lb=-1000.0,
        ub=1000.0,
    )
    source = tmp_path / "source.xml"
    source.write_bytes(b"<sbml>source fixture pinned byte-for-byte</sbml>\n")

    def fake_write_sbml(model, path):
        Path(path).write_text(f"<sbml model='{model.id}'/>\n", encoding="utf-8")

    monkeypatch.setattr(production, "read_sbml_model", lambda path: toy_model.copy())
    monkeypatch.setattr(production, "write_sbml_model", fake_write_sbml)
    output = tmp_path / "run"
    config = ProductionWorkflowConfig(
        model_path=source,
        product="PRODUCT",
        substrate="SOURCE_A",
        output_dir=output,
        run_single_knockout=False,
        run_strain_design=False,
        run_amplification=False,
        validation=ValidationConfig(enabled=False),
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = run_production_target_discovery(config)

    manifest = json.loads((output / "00_manifest.json").read_text())
    assert "amplification_target_ranking" in manifest["artifacts"]
    assert "variability_supported_amplification_targets" in manifest["artifacts"]
    pinned = output / manifest["artifacts"]["model"]["path"]
    assert result.run_directory == output.resolve()
    assert manifest["schema_version"] == 2
    assert manifest["status"] == "partial"
    assert manifest["report"]["language"] == "en"
    assert isinstance(manifest["artifacts"], dict)
    assert pinned.read_bytes() == source.read_bytes()
    conditioned_record = next(
        item
        for item in manifest["supplementary_artifacts"]
        if item["role"] == "conditioned_model"
    )
    assert (output / conditioned_record["path"]).is_file()
    assert {
        "01_preflight",
        "02_yield",
        "03_reference",
        "04_single_knockout",
        "05_strain_design",
        "06_amplification",
        "07_validation",
        "figures",
        "scripts",
    } <= {path.name for path in output.iterdir() if path.is_dir()}
    assert (output / "04_single_knockout/single_knockout_moma.csv").is_file()
    assert (output / "06_amplification/fvseof_tidy.csv").is_file()
    assert manifest["artifacts"]["single_knockout_moma"]["status"] == "skipped"
    assert manifest["artifacts"]["optknock"]["status"] == "skipped"
    assert manifest["artifacts"]["fseof_tidy"]["status"] == "skipped"
    assert manifest["artifacts"]["gene_knockout_mapping"]["status"] == "complete"
    summary = json.loads((output / "00_summary.json").read_text())
    assert summary["status"] == "partial"
    assert summary["medium_as_loaded"] is True
    assert any(
        "model medium was used as loaded" in item for item in summary["warnings"]
    )
    assert any("EX_oxygen_e" in item for item in summary["warnings"])
    assert summary["oxygen_uptake"][0]["uptake_limit"] == pytest.approx(1000.0)

    exported_config_path = output / "00_config.json"
    exported_provenance_path = output / "00_provenance.json"
    preflight_metadata_path = output / "01_preflight/preflight.metadata.json"
    exported_config = json.loads(exported_config_path.read_text())
    exported_provenance = json.loads(exported_provenance_path.read_text())
    assert exported_config["model_path"] == manifest["artifacts"]["model"]["path"]
    assert exported_config["output_dir"] == "."
    assert exported_provenance["parameters"]["source_model_path"] == str(
        manifest["artifacts"]["model"]["path"]
    )
    for portable_json in (
        exported_config_path,
        exported_provenance_path,
        preflight_metadata_path,
    ):
        assert str(tmp_path) not in portable_json.read_text()

    reloaded_export = ProductionWorkflowConfig.from_json(exported_config_path)
    assert reloaded_export.model_path == pinned.resolve()
    assert reloaded_export.output_dir == output.resolve()

    declared = [*manifest["artifacts"].values(), *manifest["supplementary_artifacts"]]
    for artifact in declared:
        path = str(artifact.get("path", ""))
        if artifact.get("status") == "complete" and (
            path.endswith(".csv") or path.endswith(".csv.gz")
        ):
            metadata_path = artifact.get("metadata_path")
            assert metadata_path, path
            assert (output / str(metadata_path)).is_file(), path

    reproduction_path = output / "scripts/production_config.json"
    reproduction = json.loads(reproduction_path.read_text())
    assert not Path(reproduction["model_path"]).is_absolute()
    assert not Path(reproduction["output_dir"]).is_absolute()
    reproduced_config = ProductionWorkflowConfig.from_json(reproduction_path)
    assert reproduced_config.model_path == pinned.resolve()
    assert reproduced_config.output_dir == (tmp_path / "run__reproduced").resolve()
    for script_name in ("reproduce.py", "render.py", "validate.py"):
        script = output / "scripts" / script_name
        compile(script.read_text(), str(script), "exec")

    # An explicit rerun must not retain prior report/figure/target artifacts, and a stale
    # stage symlink must never redirect CMM writes outside the exact run root.
    (output / "report.html").write_text("stale report", encoding="utf-8")
    (output / "report_standalone.html").write_text("stale standalone", encoding="utf-8")
    (output / "report_validation.json").write_text("{}", encoding="utf-8")
    (output / "figures/stale.png").write_bytes(b"stale")
    outside = tmp_path / "outside"
    outside.mkdir()
    shutil.rmtree(output / "07_validation")
    (output / "07_validation").symlink_to(outside, target_is_directory=True)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rerun = run_production_target_discovery(replace(config, overwrite=True))

    assert rerun.run_directory == output.resolve()
    assert not (output / "report.html").exists()
    assert not (output / "report_standalone.html").exists()
    assert not (output / "report_validation.json").exists()
    assert not (output / "figures/stale.png").exists()
    assert not (output / "07_validation").is_symlink()
    assert not any(outside.iterdir())
    assert json.loads((output / "00_manifest.json").read_text())["status"] == "partial"

    rerun_manifest = json.loads((output / "00_manifest.json").read_text())
    bundled_model = output / rerun_manifest["artifacts"]["model"]["path"]
    authoritative_bytes = bundled_model.read_bytes()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run_production_target_discovery(
            replace(config, model_path=bundled_model, overwrite=True)
        )
    refreshed_manifest = json.loads((output / "00_manifest.json").read_text())
    refreshed_model = output / refreshed_manifest["artifacts"]["model"]["path"]
    assert refreshed_model.read_bytes() == authoritative_bytes


def test_full_mocked_workflow_is_publication_valid_and_exports_candidate_metadata(
    toy_model, tmp_path, monkeypatch
) -> None:
    substrate = toy_model.reactions.SOURCE_A
    source_metabolite = next(iter(substrate.metabolites))
    substrate.add_metabolites({source_metabolite: -2.0})
    substrate.bounds = (-10.0, 1000.0)
    substrate.gene_reaction_rule = "g1"
    source_metabolite.formula = "C"
    source = tmp_path / "source.xml"
    source.write_bytes(b"<sbml>full mocked workflow source</sbml>\n")

    monkeypatch.setattr(production, "read_sbml_model", lambda path: toy_model.copy())
    monkeypatch.setattr(
        production,
        "write_sbml_model",
        lambda model, path: Path(path).write_text(
            f"<sbml model='{model.id}'/>\n", encoding="utf-8"
        ),
    )
    monkeypatch.setattr(production, "require", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        production,
        "fba",
        lambda model: FluxSolution(
            "optimal",
            1.0,
            {"SOURCE_A": -10.0, "BIOMASS": 1.0, "PRODUCT": 0.0},
            {"method": "fba"},
        ),
    )
    monkeypatch.setattr(
        production,
        "reference_flux",
        lambda model, method: FluxState(
            {"SOURCE_A": -10.0, "BIOMASS": 1.0, "PRODUCT": 0.0},
            name="reference",
            provenance="pfba",
            metadata={"method": method},
        ),
    )
    monkeypatch.setattr(
        production,
        "theoretical_yield",
        lambda model, product, substrate: ProductionYield(
            product=product,
            substrate=substrate or "SOURCE_A",
            product_flux=10.0,
            substrate_uptake=10.0,
            molar_yield=1.0,
            status="optimal",
            aerobic=False,
            carbon_ceiling=1.0,
            product_carbon=1.0,
            metadata={"method": "theoretical_yield"},
        ),
    )
    monkeypatch.setattr(
        production,
        "production_envelope",
        lambda model, product, objective, substrate, points: ProductionEnvelope(
            product,
            objective,
            (EnvelopePoint(0.0, 0.0, 1.0), EnvelopePoint(10.0, 0.0, 0.0)),
            {"method": "production_envelope", "points": points},
        ),
    )
    perturbation = Perturbation("g1", "gene", ("SOURCE_A",))
    monkeypatch.setattr(
        production, "gene_perturbations", lambda model, genes=None: [perturbation]
    )

    def fake_batch(model, reference, perturbations, *, method, **kwargs):
        del model, reference, perturbations, kwargs
        return BatchComparisonResult(
            [_screen_row("g1", product=1.5, growth=0.6)],
            metadata={"comparison_method": method, "epsilon": 1e-3},
        )

    monkeypatch.setattr(production, "batch_comparison", fake_batch)
    strain_design_calls: list[tuple[str, int]] = []
    strain_design_limits: list[float | None] = []

    def fake_optknock(*args, **kwargs):
        del args
        strain_design_calls.append(("optknock", int(kwargs["seed"])))
        strain_design_limits.append(kwargs["max_seconds"])
        return StrainDesignResult(
            "optknock",
            "PRODUCT",
            (StrainDesign(("SOURCE_A",), 0.7, 2.0, 0.0),),
            {"method": "optknock"},
        )

    def fake_robustknock(*args, **kwargs):
        del args
        strain_design_calls.append(("robustknock", int(kwargs["seed"])))
        strain_design_limits.append(kwargs["max_seconds"])
        return StrainDesignResult(
            "robustknock",
            "PRODUCT",
            (StrainDesign(("SOURCE_A",), 0.7, 2.5, 2.0),),
            {"method": "robustknock"},
        )

    monkeypatch.setattr(production, "optknock", fake_optknock)
    monkeypatch.setattr(production, "robustknock", fake_robustknock)
    levels = (0.0, 1.0)
    fseof_result = FseofResult(
        "PRODUCT",
        "BIOMASS",
        levels,
        pd.DataFrame(
            {
                0.0: [1.0],
                1.0: [2.0],
                "slope": [1.0],
                "classification": ["amplify"],
                "actionable": [True],
            },
            index=["SOURCE_A"],
        ),
        {"method": "fseof"},
    )
    mean = pd.DataFrame({0.0: [1.0], 1.0: [2.0]}, index=["SOURCE_A"])
    fvseof_result = FvseofResult(
        "PRODUCT",
        "BIOMASS",
        levels,
        mean=mean,
        forced=mean.copy(),
        capacity=pd.DataFrame({0.0: [1.0], 1.0: [0.5]}, index=["SOURCE_A"]),
        classification=pd.Series({"SOURCE_A": "amplify"}),
        park_type=pd.Series({"SOURCE_A": 1}),
        robust=pd.Series({"SOURCE_A": True}),
        slope=pd.Series({"SOURCE_A": 1.0}),
        capacity_slope=pd.Series({"SOURCE_A": -0.5}),
        actionable=pd.Series({"SOURCE_A": True}),
        metadata={"method": "fvseof"},
    )
    monkeypatch.setattr(production, "fseof", lambda *args, **kwargs: fseof_result)
    monkeypatch.setattr(production, "fvseof", lambda *args, **kwargs: fvseof_result)
    loop_record = AmplificationLoopDiagnosticRecord(
        rank=1,
        target="SOURCE_A",
        source_methods=("fseof", "fvseof"),
        standard_minimum=-10.0,
        standard_maximum=10.0,
        standard_capacity=20.0,
        loopless_minimum=-5.0,
        loopless_maximum=5.0,
        loopless_capacity=10.0,
        loopless_to_standard_capacity_ratio=0.5,
        capacity_ratio_threshold=0.1,
        loop_artifact_flag=False,
        diagnostic_status="complete",
        reason=None,
        enforced_product_floor=1.0,
        biomass_floor=0.3,
    )
    monkeypatch.setattr(
        production,
        "_run_amplification_loop_diagnostic",
        lambda *args, **kwargs: AmplificationLoopDiagnosticResult(
            (loop_record,), {"method": "loopless_fva"}
        ),
    )

    def fake_response(model, target, response=None, *, biomass, **kwargs):
        del model
        if kwargs.get("target_min") is not None:
            points = (
                ResponsePoint(-10.0, 0.2, 0.5, "optimal"),
                ResponsePoint(0.0, 1.0, 0.5, "optimal"),
            )
            wild_type = {
                "target_flux": -10.0,
                "response_flux": 0.2,
                "biomass_flux": 1.0,
            }
        elif target == "SOURCE_A":
            points = (
                ResponsePoint(1.0, 0.0, 1.0, "optimal"),
                ResponsePoint(2.0, 2.0, 0.5, "optimal"),
            )
            wild_type = {
                "target_flux": 1.0,
                "response_flux": 0.0,
                "biomass_flux": 1.0,
            }
        else:
            points = (
                ResponsePoint(0.0, 0.6, 0.6, "optimal"),
                ResponsePoint(1.5, 0.5, 0.5, "optimal"),
            )
            wild_type = {
                "target_flux": 0.0,
                "response_flux": 0.6,
                "biomass_flux": 0.6,
            }
        return FluxResponseResult(
            target,
            response or biomass,
            biomass,
            points,
            (),
            ResponseLimit(False, "no limit"),
            wild_type,
            metadata={"method": "flux_response"},
        )

    monkeypatch.setattr(production, "flux_response", fake_response)

    def fake_sampling(model, **kwargs):
        knocked_out = model.reactions.SOURCE_A.bounds == (0.0, 0.0)
        product_values = [1.0, 1.1, 1.2] if knocked_out else [0.0, 0.1, 0.2]
        biomass_values = [0.6, 0.6, 0.6] if knocked_out else [1.0, 1.0, 1.0]
        return SamplingResult(
            pd.DataFrame(
                {
                    "SOURCE_A": [-5.0, -5.0, -5.0],
                    "BIOMASS": biomass_values,
                    "PRODUCT": product_values,
                }
            ),
            method=str(kwargs["method"]),
            seed=int(kwargs["seed"]),
            metadata={"method": "random_flux_sampling", "knockout": knocked_out},
        )

    monkeypatch.setattr(production, "random_flux_sampling", fake_sampling)
    output = tmp_path / "full-run"
    config = ProductionWorkflowConfig(
        model_path=source,
        product="PRODUCT",
        substrate="SOURCE_A",
        output_dir=output,
        strain_design_seed=53,
        strain_design_max_seconds=90.0,
        top_single_knockouts_per_method=1,
        top_amplification_targets_per_method=1,
        validation=ValidationConfig(
            max_flux_response_targets=4,
            flux_response_steps=3,
            sampling=SamplingConfig(n=3, method="achr", thinning=1),
        ),
    )

    result = run_production_target_discovery(config)
    validated = validate_run(output)

    assert strain_design_calls == [("optknock", 53), ("robustknock", 53)]
    assert strain_design_limits == [90.0, 90.0]
    assert result.summary()["status"] == "complete"
    assert validated.artifact("gene_knockout_mapping", required=False) is not None
    recommendations = pd.read_csv(output / "07_validation/recommendations.csv")
    assert set(recommendations["type"]) == {
        "single_gene_knockout",
        "amplification",
        "multi_knockout",
    }
    manifest = json.loads((output / "00_manifest.json").read_text())
    exported_config = json.loads((output / "00_config.json").read_text())
    exported_provenance = json.loads((output / "00_provenance.json").read_text())
    assert exported_config["strain_design_seed"] == 53
    assert exported_provenance["parameters"]["strain_design_seed"] == 53
    assert exported_config["strain_design_max_seconds"] == 90.0
    assert exported_provenance["parameters"]["strain_design_max_seconds"] == 90.0
    assert "amplification_target_ranking" in manifest["artifacts"]
    assert "variability_supported_amplification_targets" in manifest["artifacts"]
    for filename in ("fseof.csv", "fvseof.csv"):
        ranking = pd.read_csv(output / "06_amplification" / filename)
        assert {
            "amplification_rank",
            "proposal_method",
            "method_rank",
            "report_selected",
            "loop_diagnostic_status",
            "loop_artifact_flag",
            "loopless_to_standard_capacity_ratio",
            "loop_diagnostic_reason",
        } <= set(ranking.columns)
        assert int(ranking["report_selected"].sum()) == 1
    for filename, method in (
        ("fseof_tidy.csv", "fseof"),
        ("fvseof_tidy.csv", "fvseof"),
    ):
        tidy = pd.read_csv(output / "06_amplification" / filename)
        assert set(tidy["proposal_method"]) == {method}
        assert set(tidy["method_rank"]) == {1}
        assert tidy["report_selected"].all()
        assert tidy["loop_diagnostic_status"].eq("complete").all()
    room_metadata = json.loads(
        (output / "04_single_knockout/single_knockout_room.metadata.json").read_text()
    )["analysis_metadata"]
    assert room_metadata["selection_room_epsilon"] == pytest.approx(1e-3)
    assert room_metadata["selection_numeric_margin"] == pytest.approx(
        config.product_improvement_tolerance
    )
    assert room_metadata["selection_product_delta_threshold"] == pytest.approx(
        1e-3 + config.product_improvement_tolerance
    )
    assert (
        room_metadata["selection_product_delta_criterion"]
        == "strictly_greater_than_threshold"
    )
    assert room_metadata["display_rank_limit"] == 1
    assert room_metadata["display_rank_viability_fraction"] == pytest.approx(
        config.viability_fraction
    )
    moma_rows = pd.read_csv(output / "04_single_knockout/single_knockout_moma.csv")
    room_rows = pd.read_csv(output / "04_single_knockout/single_knockout_room.csv")
    knockout_candidates = pd.read_csv(
        output / "04_single_knockout/single_knockout_candidates.csv"
    )
    mapping_rows = pd.read_csv(output / "04_single_knockout/gene_knockout_mapping.csv")
    reaction_columns = {
        "display_rank",
        "target_name",
        "blocked_reactions",
        "blocked_reaction_names",
        "blocked_reaction_equations",
        "blocked_reaction_gprs",
    }
    assert reaction_columns <= set(moma_rows.columns)
    assert reaction_columns <= set(room_rows.columns)
    assert "reaction_equation" in mapping_rows.columns
    assert set(knockout_candidates["candidate_scope"]) == {
        "all_display_ranked_candidates"
    }
    assert {
        "validation_target_id",
        "candidate_target_ids",
        "candidate_source_methods",
        "validation_representative",
    } <= set(knockout_candidates.columns)
    assert moma_rows.loc[moma_rows["target_id"] == "g1", "display_rank"].notna().any()
    assert (
        moma_rows.loc[moma_rows["target_id"] == "g1", "blocked_reaction_equations"]
        .str.len()
        .gt(0)
        .all()
    )
    response_index = pd.read_csv(output / "07_validation/flux_response_index.csv")
    assert {
        "candidate_scope",
        "blocked_reaction_signature",
        "candidate_target_ids",
        "loop_diagnostic_status",
        "loop_artifact_flag",
        "loop_diagnostic_eligible",
        "loop_diagnostic_reason",
        "scan_reference_flux",
        "reason",
    } <= set(response_index.columns)
    assert set(response_index["candidate_scope"]) == {
        "all_display_ranked_candidates",
        "all_report_selected_candidates",
    }
    assert set(response_index["background"]) == {"wild_type"}
    knockout_response = response_index.loc[
        response_index["candidate_scope"] == "all_display_ranked_candidates"
    ].iloc[0]
    assert knockout_response["scan_reaction"] == "SOURCE_A"
    assert knockout_response["response_reaction"] == "PRODUCT"
    assert knockout_response["scan_reference_flux"] == pytest.approx(-10.0)
    response_tidy = pd.read_csv(output / "07_validation/flux_response_tidy.csv")
    assert set(response_tidy["candidate_scope"]) == {
        "all_display_ranked_candidates",
        "all_report_selected_candidates",
    }
    sampling_index = pd.read_csv(output / "07_validation/random_sampling_index.csv")
    assert {
        "candidate_scope",
        "blocked_reaction_signature",
        "candidate_target_ids",
        "reason",
    } <= set(sampling_index.columns)
    summary = json.loads((output / "00_summary.json").read_text())
    assert summary["n_selected_single_knockout_method_rows"] == 2
    assert summary["n_beneficial_single_knockout_candidates"] == 1
    assert summary["validation_candidate_policy"] == {
        "single_knockout": "all_display_ranked_candidates",
        "amplification": "all_report_selected_candidates",
        "gpr_deduplication": "blocked_reaction_signature_representative",
        "flux_response_axes": "candidate_reaction_flux_to_target_product_flux",
        "single_knockout_flux_response": (
            "pre_deletion_reference_to_zero_or_full_domain_when_reference_is_zero"
        ),
        "multi_reaction_knockout_flux_response": "explicitly_unavailable",
    }
    assert summary["validation_coverage"] == {
        "single_knockout_candidates_expected": 1,
        "amplification_candidates_expected": 1,
        "flux_response_expected": 2,
        "flux_response_attempted": 2,
        "flux_response_completed": 2,
        "flux_response_failed": 0,
        "flux_response_skipped": 0,
        "sampling_expected": 2,
        "sampling_attempted": 2,
        "sampling_completed": 2,
        "sampling_failed": 0,
        "sampling_skipped": 0,
    }
    candidate_records = [
        item
        for item in manifest["supplementary_artifacts"]
        if item["role"]
        in {
            "candidate_flux_response",
            "candidate_flux_response_phases",
            "raw_flux_samples",
            "flux_sampling_statistics",
            "wild_type_knockout_distribution_shift",
        }
    ]
    assert candidate_records
    assert all(item["metadata_path"] for item in candidate_records)
    assert all((output / item["metadata_path"]).is_file() for item in candidate_records)
