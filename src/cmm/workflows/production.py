"""Canonical production-target-discovery workflow.

This module is the goal-level orchestration layer for SC-01.  Numerical work stays in the
existing services; the workflow fixes their order, gates solver capabilities before a long
run, ranks single-gene deletions reproducibly, validates predictions, and writes one
machine-readable schema-v2 run directory.

The public entry point deliberately accepts a configuration containing ``model_path``.
That makes the exact input model part of the invocation rather than hidden caller state.  A
loaded-model helper exists only for focused tests and integrations that already own model
loading.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from numbers import Integral
from typing import Literal, cast

from cobra import Model
from cobra.io import read_sbml_model, write_sbml_model
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

from cmm.core import (
    Condition,
    FluxRange,
    FluxSolution,
    FluxState,
    Medium,
    MediumApplication,
    ObjectiveSpec,
    ReactionBound,
    apply_medium,
    fba,
    fva,
    model_fingerprint,
    require,
    run_provenance,
    solver_status,
)
from cmm.features import (
    BatchComparisonResult,
    FluxResponseResult,
    FseofResult,
    FvseofResult,
    Perturbation,
    ProductionEnvelope,
    ProductionYield,
    SamplingResult,
    StrainDesignResult,
    batch_comparison,
    blocked_reactions_for_genes,
    flux_response,
    fseof,
    fvseof,
    gene_perturbations,
    optknock,
    production_envelope,
    random_flux_sampling,
    reference_flux,
    robustknock,
    theoretical_yield,
)
from cmm.features._perturbation import apply_perturbation
from cmm.features.strain_design import STRAIN_DESIGN_MAX_SECONDS, validate_max_seconds

SCHEMA_ID = "cmm.production-target-discovery"
SCHEMA_VERSION = 2

_STAGE_DIRECTORIES = (
    "01_preflight",
    "02_yield",
    "03_reference",
    "04_single_knockout",
    "05_strain_design",
    "06_amplification",
    "07_validation",
    "model",
    "figures",
    "scripts",
)

ReferenceMethod = Literal["fba", "pfba"]
SamplingMethod = Literal["optgp", "achr"]
SingleKnockoutMethod = Literal["moma_l2", "room"]
ArtifactStatus = Literal["complete", "partial", "skipped", "failed"]

# The bundle writer lives in a shared module so this workflow and transformation-target
# discovery emit one format and one validator can serve both. Imported under the original
# names: every reference below is unchanged.
from cmm.workflows._bundle import (  # noqa: E402
    ArtifactRecord,
    _ArtifactWriter,
    _jsonable,
)

LooplessAlgorithm = Literal["cycleFreeFlux", "fastSNP"]

_STRAIN_DESIGN_SEED_MAX = 2_000_000_000

VALIDATION_CANDIDATE_POLICY: Mapping[str, str] = {
    "single_knockout": "all_display_ranked_candidates",
    "amplification": "all_report_selected_candidates",
    "gpr_deduplication": "blocked_reaction_signature_representative",
    "flux_response_axes": "candidate_reaction_flux_to_target_product_flux",
    "single_knockout_flux_response": (
        "pre_deletion_reference_to_zero_or_full_domain_when_reference_is_zero"
    ),
    "multi_reaction_knockout_flux_response": "explicitly_unavailable",
}


class ProductionWorkflowError(RuntimeError):
    """A preflight or orchestration failure that invalidates the scientific run."""


@dataclass(frozen=True)
class SamplingConfig:
    """Random-sampling settings for wild type and all display-ranked KO candidates."""

    enabled: bool = True
    n: int = 1000
    method: SamplingMethod = "achr"
    thinning: int = 100
    processes: int = 1
    seed: int = 0
    store_raw_samples: bool = True

    def validate(self) -> None:
        if self.n < 1:
            raise ValueError("sampling n must be at least 1")
        if self.thinning < 1:
            raise ValueError("sampling thinning must be at least 1")
        if self.processes < 1:
            raise ValueError("sampling processes must be at least 1")
        if self.method not in {"optgp", "achr"}:
            raise ValueError("sampling method must be 'optgp' or 'achr'")
        if self.method == "achr" and self.processes != 1:
            raise ValueError("the ACHR sampler requires processes=1")


@dataclass(frozen=True)
class ValidationConfig:
    """Forward-validation settings for inverse and single-knockout predictions."""

    enabled: bool = True
    max_flux_response_targets: int = 30
    flux_response_steps: int = 20
    flux_response_biomass_fraction: float = 0.3
    sampling_growth_fraction: float = 0.1
    sampling: SamplingConfig = field(default_factory=SamplingConfig)

    def validate(self) -> None:
        if self.max_flux_response_targets < 1:
            raise ValueError("max_flux_response_targets must be at least 1")
        if self.flux_response_steps < 2:
            raise ValueError("flux_response_steps must be at least 2")
        if not 0.0 < self.flux_response_biomass_fraction <= 1.0:
            raise ValueError("flux_response_biomass_fraction must be in (0, 1]")
        if not 0.0 < self.sampling_growth_fraction <= 1.0:
            raise ValueError("sampling_growth_fraction must be in (0, 1]")
        self.sampling.validate()


@dataclass(frozen=True)
class ProductionWorkflowConfig:
    """Complete, serializable invocation of production-target discovery.

    ``medium=None`` means the bounds stored in the model are retained.  It is recorded as
    such; it is not interpreted as an unnamed preset.  Disabling a stage is explicit and is
    represented as ``skipped`` in the manifest, so a missing solver never causes a silent
    method substitution or omission.
    """

    model_path: str | Path
    product: str
    output_dir: str | Path | None = None
    solver: str | None = None
    substrate: str | None = None
    biomass: str | None = None
    medium: Medium | str | None = None
    condition: Condition | None = None
    reference_method: ReferenceMethod = "pfba"
    envelope_points: int = 20
    run_single_knockout: bool = True
    single_knockout_genes: tuple[str, ...] | None = None
    top_single_knockouts_per_method: int = 5
    viability_fraction: float = 0.1
    product_improvement_tolerance: float = 1e-6
    run_strain_design: bool = True
    max_knockouts: int = 3
    optknock_max_solutions: int = 5
    robustknock_max_solutions: int = 8
    strain_design_seed: int = 0
    strain_design_max_seconds: float | None = STRAIN_DESIGN_MAX_SECONDS
    design_min_growth: float = 0.05
    actionable_designs_only: bool = True
    run_amplification: bool = True
    fseof_steps: int = 10
    fvseof_steps: int = 10
    top_amplification_targets_per_method: int = 10
    scan_fraction_min: float = 0.1
    scan_fraction_max: float = 0.9
    fvseof_biomass_fraction: float = 0.95
    run_amplification_loop_diagnostic: bool = True
    amplification_loop_diagnostic_top_n: int = 20
    loopless_capacity_ratio_threshold: float = 0.1
    loopless_algorithm: LooplessAlgorithm = "fastSNP"
    validation: ValidationConfig = field(default_factory=ValidationConfig)
    overwrite: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "model_path", Path(self.model_path))
        if self.output_dir is not None:
            object.__setattr__(self, "output_dir", Path(self.output_dir))
        if self.single_knockout_genes is not None:
            object.__setattr__(
                self,
                "single_knockout_genes",
                tuple(str(gene) for gene in self.single_knockout_genes),
            )
        self.validate()

    def validate(self) -> None:
        if not str(self.model_path):
            raise ValueError("model_path must not be empty")
        if not self.product:
            raise ValueError("product must not be empty")
        if self.reference_method not in {"fba", "pfba"}:
            raise ValueError("reference_method must be 'fba' or 'pfba'")
        if self.envelope_points < 2:
            raise ValueError("envelope_points must be at least 2")
        if self.top_single_knockouts_per_method < 1:
            raise ValueError("top_single_knockouts_per_method must be at least 1")
        if not 0.0 <= self.viability_fraction <= 1.0:
            raise ValueError("viability_fraction must be in [0, 1]")
        if self.product_improvement_tolerance < 0:
            raise ValueError("product_improvement_tolerance must be non-negative")
        if self.max_knockouts < 1:
            raise ValueError("max_knockouts must be at least 1")
        if self.optknock_max_solutions < 1 or self.robustknock_max_solutions < 1:
            raise ValueError("strain-design max_solutions values must be at least 1")
        if (
            isinstance(self.strain_design_seed, bool)
            or not isinstance(self.strain_design_seed, Integral)
            or not 0 <= self.strain_design_seed <= _STRAIN_DESIGN_SEED_MAX
        ):
            raise ValueError(
                "strain_design_seed must be an integer in the Gurobi-compatible range "
                f"[0, {_STRAIN_DESIGN_SEED_MAX}]"
            )
        validate_max_seconds(
            self.strain_design_max_seconds, "strain_design_max_seconds"
        )
        if self.design_min_growth < 0:
            raise ValueError("design_min_growth must be non-negative")
        if self.fseof_steps < 2 or self.fvseof_steps < 2:
            raise ValueError("FSEOF and FVSEOF steps must each be at least 2")
        if self.top_amplification_targets_per_method < 1:
            raise ValueError("top_amplification_targets_per_method must be at least 1")
        if not 0.0 <= self.scan_fraction_min < self.scan_fraction_max <= 1.0:
            raise ValueError(
                "scan fractions must satisfy 0 <= fraction_min < fraction_max <= 1"
            )
        if not 0.0 < self.fvseof_biomass_fraction <= 1.0:
            raise ValueError("fvseof_biomass_fraction must be in (0, 1]")
        if self.amplification_loop_diagnostic_top_n < 1:
            raise ValueError("amplification_loop_diagnostic_top_n must be at least 1")
        required_loop_targets = 2 * self.top_amplification_targets_per_method
        if (
            self.run_amplification
            and self.run_amplification_loop_diagnostic
            and self.amplification_loop_diagnostic_top_n < required_loop_targets
        ):
            raise ValueError(
                "amplification_loop_diagnostic_top_n must accommodate the independent "
                "FSEOF and FVSEOF shortlists (at least "
                f"{required_loop_targets})"
            )
        if not 0.0 <= self.loopless_capacity_ratio_threshold <= 1.0:
            raise ValueError("loopless_capacity_ratio_threshold must be in [0, 1]")
        if self.loopless_algorithm not in {"cycleFreeFlux", "fastSNP"}:
            raise ValueError("loopless_algorithm must be 'cycleFreeFlux' or 'fastSNP'")
        self.validation.validate()
        required_validation_targets = (
            2 * self.top_amplification_targets_per_method
            if self.run_amplification
            else 0
        ) + (
            2 * self.top_single_knockouts_per_method if self.run_single_knockout else 0
        )
        if (
            self.validation.enabled
            and self.validation.max_flux_response_targets < required_validation_targets
        ):
            raise ValueError(
                "validation.max_flux_response_targets must accommodate both independent "
                "report-selected amplification candidate lists and every per-method "
                "display-ranked single-KO candidate "
                f"(at least {required_validation_targets})"
            )

    @classmethod
    def from_json(cls, path: str | Path) -> ProductionWorkflowConfig:
        """Load a workflow configuration from a UTF-8 JSON object."""

        config_path = Path(path).expanduser().resolve()
        payload = json.loads(config_path.read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise ValueError("production workflow JSON must contain an object")
        values = dict(payload)
        for field_name in ("model_path", "output_dir"):
            raw_path = values.get(field_name)
            if raw_path is None:
                continue
            candidate = Path(str(raw_path)).expanduser()
            if not candidate.is_absolute():
                candidate = config_path.parent / candidate
            values[field_name] = candidate.resolve()
        return cls.from_mapping(values)

    @classmethod
    def from_mapping(cls, payload: Mapping[str, object]) -> ProductionWorkflowConfig:
        """Build a config from CLI-friendly mappings, including nested conditions."""

        values = dict(payload)
        values["medium"] = _medium_from_payload(values.get("medium"))
        values["condition"] = _condition_from_payload(values.get("condition"))
        validation = values.get("validation")
        if isinstance(validation, Mapping):
            validation_values = dict(validation)
            sampling = validation_values.get("sampling")
            if isinstance(sampling, Mapping):
                validation_values["sampling"] = SamplingConfig(**dict(sampling))
            values["validation"] = ValidationConfig(**validation_values)
        genes = values.get("single_knockout_genes")
        if isinstance(genes, Sequence) and not isinstance(genes, (str, bytes)):
            values["single_knockout_genes"] = tuple(str(gene) for gene in genes)
        return cls(**values)  # type: ignore[arg-type]


@dataclass(frozen=True)
class PreflightRecord:
    check: str
    status: Literal["pass", "warning"]
    value: object
    message: str


@dataclass(frozen=True)
class SingleKnockoutRecord:
    """One full-screen row plus its deterministic candidate-selection annotation."""

    method: SingleKnockoutMethod
    target_id: str
    blocked_reactions: tuple[str, ...]
    blocked_reaction_signature: str
    status: str
    growth_rate: float
    growth_fraction: float
    target_production: float
    product_delta: float
    product_fold_change: float | None
    objective_value: float
    distance: float | None
    distance_kind: str
    n_changed_reactions: float | None
    selected: bool = False
    method_rank: int | None = None
    display_rank: int | None = None
    improves_product: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "method": self.method,
            "method_rank": self.method_rank,
            "display_rank": self.display_rank,
            "selected": self.selected,
            "target_id": self.target_id,
            "kind": "gene",
            "blocked_reactions": ";".join(self.blocked_reactions),
            "blocked_reaction_signature": self.blocked_reaction_signature,
            "status": self.status,
            "growth_rate": self.growth_rate,
            "objective": self.growth_rate,
            "growth_fraction": self.growth_fraction,
            "target_production": self.target_production,
            "product_flux": self.target_production,
            "product_delta": self.product_delta,
            "product_fold_change": self.product_fold_change,
            "improves_product": self.improves_product,
            "objective_value": self.objective_value,
            "distance": self.distance,
            "distance_kind": self.distance_kind,
            "n_changed_reactions": self.n_changed_reactions,
            "n_blocked_reactions": len(self.blocked_reactions),
        }


@dataclass(frozen=True)
class GeneKnockoutMappingRecord:
    """One gene-to-blocked-reaction row, retaining human names and the full GPR."""

    gene_id: str
    gene_name: str
    inert: bool
    blocked_reaction: str | None
    reaction_name: str | None
    reaction_equation: str | None
    gpr: str | None


@dataclass(frozen=True)
class GeneKnockoutMappingResult:
    records: tuple[GeneKnockoutMappingRecord, ...]
    metadata: Mapping[str, object]

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame(
            [asdict(record) for record in self.records],
            columns=[
                "gene_id",
                "gene_name",
                "inert",
                "blocked_reaction",
                "reaction_name",
                "reaction_equation",
                "gpr",
            ],
        )


@dataclass(frozen=True)
class AmplificationLoopDiagnosticRecord:
    """Loopless-vs-standard FVA evidence for one amplification candidate."""

    rank: int
    target: str
    source_methods: tuple[str, ...]
    standard_minimum: float | None
    standard_maximum: float | None
    standard_capacity: float | None
    loopless_minimum: float | None
    loopless_maximum: float | None
    loopless_capacity: float | None
    loopless_to_standard_capacity_ratio: float | None
    capacity_ratio_threshold: float
    loop_artifact_flag: bool | None
    diagnostic_status: Literal["complete", "inconclusive", "failed"]
    reason: str | None
    enforced_product_floor: float
    biomass_floor: float


@dataclass(frozen=True)
class AmplificationLoopDiagnosticResult:
    records: tuple[AmplificationLoopDiagnosticRecord, ...]
    metadata: Mapping[str, object]

    def to_frame(self) -> pd.DataFrame:
        rows = []
        for record in self.records:
            row = asdict(record)
            row["source_methods"] = ";".join(record.source_methods)
            rows.append(row)
        return pd.DataFrame(
            rows,
            columns=[
                "rank",
                "target",
                "source_methods",
                "standard_minimum",
                "standard_maximum",
                "standard_capacity",
                "loopless_minimum",
                "loopless_maximum",
                "loopless_capacity",
                "loopless_to_standard_capacity_ratio",
                "capacity_ratio_threshold",
                "loop_artifact_flag",
                "diagnostic_status",
                "reason",
                "enforced_product_floor",
                "biomass_floor",
            ],
        )


@dataclass(frozen=True)
class ValidationTarget:
    target_id: str
    scan_reaction: str
    response_reaction: str
    background: Literal["wild_type", "gene_knockout"]
    actions: tuple[str, ...]
    source_methods: tuple[str, ...]
    blocked_reactions: tuple[str, ...] = ()
    candidate_scope: str = ""
    blocked_reaction_signature: str | None = None
    candidate_target_ids: tuple[str, ...] = ()
    loop_diagnostic_status: str | None = None
    loop_artifact_flag: bool | None = None
    loop_diagnostic_eligible: bool | None = None
    loop_diagnostic_reason: str | None = None
    scan_reference_flux: float | None = None


@dataclass(frozen=True)
class FluxResponseValidation:
    target: ValidationTarget
    status: Literal["complete", "failed", "skipped"]
    result: FluxResponseResult | None = None
    error: str | None = None
    reason: str | None = None


@dataclass(frozen=True)
class SamplingValidation:
    target_id: str
    blocked_reactions: tuple[str, ...]
    source_methods: tuple[str, ...]
    status: Literal["complete", "failed", "skipped"]
    result: SamplingResult | None = None
    comparison: pd.DataFrame | None = None
    error: str | None = None
    reason: str | None = None
    candidate_scope: str = ""
    blocked_reaction_signature: str | None = None
    candidate_target_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.comparison is not None:
            object.__setattr__(self, "comparison", self.comparison.copy())


def _validation_coverage(
    targets: Sequence[ValidationTarget],
    responses: Sequence[FluxResponseValidation],
    sampling_results: Sequence[SamplingValidation],
    *,
    sampling_enabled: bool,
) -> dict[str, int]:
    """Return explicit expected/attempted/completed validation counts."""

    knockout_expected = sum(
        target.candidate_scope == "all_display_ranked_candidates"
        or (not target.candidate_scope and target.background == "gene_knockout")
        for target in targets
    )
    amplification_expected = sum(
        target.candidate_scope == "all_report_selected_candidates"
        or (not target.candidate_scope and target.background == "wild_type")
        for target in targets
    )
    sampling_expected = (
        knockout_expected + 1 if sampling_enabled and knockout_expected else 0
    )
    return {
        "single_knockout_candidates_expected": knockout_expected,
        "amplification_candidates_expected": amplification_expected,
        "flux_response_expected": len(targets),
        "flux_response_attempted": sum(
            item.status in {"complete", "failed"} for item in responses
        ),
        "flux_response_completed": sum(item.status == "complete" for item in responses),
        "flux_response_failed": sum(item.status == "failed" for item in responses),
        "flux_response_skipped": sum(item.status == "skipped" for item in responses),
        "sampling_expected": sampling_expected,
        "sampling_attempted": sum(
            item.status in {"complete", "failed"} for item in sampling_results
        ),
        "sampling_completed": sum(
            item.status == "complete" for item in sampling_results
        ),
        "sampling_failed": sum(item.status == "failed" for item in sampling_results),
        "sampling_skipped": sum(item.status == "skipped" for item in sampling_results),
    }


@dataclass(frozen=True)
class ProductionWorkflowResult:
    """All in-memory results and the authoritative exported-artifact index."""

    config: ProductionWorkflowConfig
    preflight: tuple[PreflightRecord, ...]
    medium_application: MediumApplication
    biomass: str
    wild_type: FluxSolution
    reference: FluxState
    theoretical_yield: ProductionYield
    envelope: ProductionEnvelope
    gene_knockout_mapping: GeneKnockoutMappingResult
    single_knockouts: tuple[SingleKnockoutRecord, ...]
    optknock_result: StrainDesignResult | None
    robustknock_result: StrainDesignResult | None
    fseof_result: FseofResult | None
    fvseof_result: FvseofResult | None
    amplification_loop_diagnostic: AmplificationLoopDiagnosticResult | None
    validation_targets: tuple[ValidationTarget, ...]
    flux_responses: tuple[FluxResponseValidation, ...]
    sampling: tuple[SamplingValidation, ...]
    provenance: Mapping[str, object]
    run_directory: Path | None = None
    artifacts: tuple[ArtifactRecord, ...] = ()

    @property
    def selected_single_knockouts(self) -> tuple[SingleKnockoutRecord, ...]:
        return tuple(record for record in self.single_knockouts if record.selected)

    @property
    def single_knockout_validation_candidates(
        self,
    ) -> tuple[SingleKnockoutRecord, ...]:
        """Return every method-specific D1-D5 candidate before GPR deduplication."""

        return tuple(
            record
            for record in self.single_knockouts
            if record.display_rank is not None
        )

    def summary(self) -> dict[str, object]:
        selected = self.selected_single_knockouts
        beneficial = tuple(record for record in selected if record.improves_product)
        validation_coverage = _validation_coverage(
            self.validation_targets,
            self.flux_responses,
            self.sampling,
            sampling_enabled=(
                self.config.validation.enabled
                and self.config.run_single_knockout
                and self.config.validation.sampling.enabled
            ),
        )
        warnings = _summary_warnings(self)
        recommendations = _recommendations_frame(self)
        recommendation_counts = (
            recommendations.groupby("type").size().astype(int).to_dict()
            if not recommendations.empty
            else {}
        )
        return {
            "schema_id": SCHEMA_ID,
            "schema_version": SCHEMA_VERSION,
            "status": _workflow_analysis_status(self),
            "model_id": self.provenance.get("model_id"),
            "model_sha256": self.provenance.get("model_sha256"),
            "source_model_sha256": self.provenance.get("source_model_sha256"),
            "conditioned_model_sha256": self.provenance.get("conditioned_model_sha256"),
            "product": self.config.product,
            "validation_candidate_policy": dict(VALIDATION_CANDIDATE_POLICY),
            "validation_coverage": validation_coverage,
            "substrate": self.theoretical_yield.substrate,
            "biomass": self.biomass,
            "wild_type_growth": self.reference.get(self.biomass),
            "wild_type_product": self.reference.get(self.config.product),
            "theoretical_product_flux": self.theoretical_yield.product_flux,
            "theoretical_molar_yield": self.theoretical_yield.molar_yield,
            "medium": self.medium_application.medium,
            "medium_as_loaded": self.medium_application.medium == "model_as_loaded",
            "co2_carbon_fraction": self.theoretical_yield.co2_carbon_fraction,
            "carbon_imbalance": self.theoretical_yield.carbon_imbalance,
            "oxygen_uptake": self.provenance.get("oxygen_uptake", []),
            "n_single_knockouts_screened": len(
                {record.target_id for record in self.single_knockouts}
            ),
            "n_single_knockout_method_rows": len(self.single_knockouts),
            "n_genes_mapped": len(
                {record.gene_id for record in self.gene_knockout_mapping.records}
            ),
            "n_inert_genes": len(
                {
                    record.gene_id
                    for record in self.gene_knockout_mapping.records
                    if record.inert
                }
            ),
            "n_single_knockout_candidates": len(
                {
                    record.blocked_reaction_signature
                    for record in self.single_knockout_validation_candidates
                }
            ),
            "n_single_knockout_candidate_method_rows": len(
                self.single_knockout_validation_candidates
            ),
            "n_selected_single_knockout_method_rows": len(selected),
            "n_beneficial_single_knockout_candidates": len(
                {record.blocked_reaction_signature for record in beneficial}
            ),
            "n_optknock_designs": (
                len(self.optknock_result.designs) if self.optknock_result else 0
            ),
            "n_robustknock_designs": (
                len(self.robustknock_result.designs) if self.robustknock_result else 0
            ),
            "n_fseof_amplification_targets": (
                len(self.fseof_result.amplification_targets())
                if self.fseof_result
                else 0
            ),
            "n_fvseof_amplification_targets": (
                len(self.fvseof_result.amplification_targets())
                if self.fvseof_result
                else 0
            ),
            "n_flux_response_complete": sum(
                item.status == "complete" for item in self.flux_responses
            ),
            "n_flux_response_failed": sum(
                item.status == "failed" for item in self.flux_responses
            ),
            "n_flux_response_skipped": sum(
                item.status == "skipped" for item in self.flux_responses
            ),
            "n_sampling_complete": sum(
                item.status == "complete" for item in self.sampling
            ),
            "n_sampling_failed": sum(item.status == "failed" for item in self.sampling),
            "n_sampling_skipped": sum(
                item.status == "skipped" for item in self.sampling
            ),
            "n_amplification_loop_artifact_flags": (
                sum(
                    record.loop_artifact_flag is True
                    for record in self.amplification_loop_diagnostic.records
                )
                if self.amplification_loop_diagnostic is not None
                else 0
            ),
            "recommendation_counts": recommendation_counts,
            "warnings": warnings,
            "n_warnings": len(warnings),
        }


def run_production_target_discovery(
    config: ProductionWorkflowConfig,
) -> ProductionWorkflowResult:
    """Load ``config.model_path`` and run canonical production-target discovery.

    The model must be SBML/XML.  Solver capabilities are checked before the knockout or
    strain-design screens.  In particular, L2 MOMA is never replaced by L1 MOMA and exact
    ROOM is never replaced by its LP relaxation.
    """

    model_path = cast(Path, config.model_path).expanduser().resolve()
    if not model_path.is_file():
        raise FileNotFoundError(f"model_path is not a file: {model_path}")
    if model_path.suffix.lower() not in {".xml", ".sbml"}:
        raise ValueError("production workflow currently accepts SBML .xml/.sbml models")
    model = read_sbml_model(str(model_path))
    if config.solver is not None:
        model.solver = config.solver
    return _run_production_target_discovery(model, config)


def _run_production_target_discovery(
    model: Model,
    config: ProductionWorkflowConfig,
) -> ProductionWorkflowResult:
    """Loaded-model implementation used by the public file-backed runner and tests."""

    config.validate()
    original_fingerprint = model_fingerprint(model)
    working = model.copy()
    if config.solver is not None:
        working.solver = config.solver
    medium_application = _condition_model(working, config)
    biomass = _preflight_model(working, config)
    gene_mapping = _gene_knockout_mapping(working)
    oxygen_uptake = _oxygen_exchange_uptake(working)

    wild_type = fba(working)
    if wild_type.status != "optimal" or not wild_type.fluxes:
        raise ProductionWorkflowError(
            f"wild-type FBA is {wild_type.status}; production discovery requires growth"
        )

    reference = reference_flux(working, method=config.reference_method)
    reference_growth = reference.get(biomass)
    if not math.isfinite(reference_growth) or reference_growth <= 1e-9:
        raise ProductionWorkflowError(
            "wild-type reference has no positive biomass flux under the applied condition"
        )

    yield_result = theoretical_yield(
        working,
        config.product,
        config.substrate,
    )
    if yield_result.status != "optimal" or yield_result.product_flux <= 1e-9:
        raise ProductionWorkflowError(
            f"theoretical product flux for {config.product!r} is zero; change the medium, "
            "substrate, or aeration before searching targets"
        )
    envelope = production_envelope(
        working,
        config.product,
        objective=biomass,
        substrate=yield_result.substrate,
        points=config.envelope_points,
    )

    perturbations: tuple[Perturbation, ...] = ()
    single_records: tuple[SingleKnockoutRecord, ...] = ()
    single_metadata: dict[str, object] = {}
    if config.run_single_knockout:
        enumerated = gene_perturbations(
            working,
            genes=config.single_knockout_genes,
        )
        if not enumerated:
            raise ProductionWorkflowError(
                "single-gene screening found no gene knockout that blocks a reaction"
            )
        perturbations = tuple(enumerated)
        moma_screen = batch_comparison(
            working,
            reference,
            perturbations,
            method="moma_l2",
            objective_reaction=biomass,
            product_reaction=config.product,
        )
        room_screen = batch_comparison(
            working,
            reference,
            perturbations,
            method="room",
            room_use_case="flux_prediction",
            objective_reaction=biomass,
            product_reaction=config.product,
        )
        room_epsilon = float(cast(float, room_screen.metadata["epsilon"]))
        room_improvement_threshold = room_epsilon + config.product_improvement_tolerance
        single_records = (
            *_rank_single_knockout_screen(
                moma_screen,
                perturbations,
                method="moma_l2",
                wild_type_growth=reference_growth,
                wild_type_product=reference.get(config.product),
                limit=config.top_single_knockouts_per_method,
                viability_fraction=config.viability_fraction,
                improvement_threshold=config.product_improvement_tolerance,
            ),
            *_rank_single_knockout_screen(
                room_screen,
                perturbations,
                method="room",
                wild_type_growth=reference_growth,
                wild_type_product=reference.get(config.product),
                limit=config.top_single_knockouts_per_method,
                viability_fraction=config.viability_fraction,
                improvement_threshold=room_improvement_threshold,
            ),
        )
        single_metadata = {
            "moma_l2": {
                **moma_screen.metadata,
                "display_rank_limit": config.top_single_knockouts_per_method,
                "display_rank_viability_fraction": config.viability_fraction,
                "display_rank_criterion": (
                    "optimal and viability-qualified; descending product flux, descending "
                    "growth, blocked-reaction signature, target id; unique signatures"
                ),
                "selection_product_delta_threshold": (
                    config.product_improvement_tolerance
                ),
                "selection_product_delta_criterion": "strictly_greater_than_threshold",
            },
            "room": {
                **room_screen.metadata,
                "display_rank_limit": config.top_single_knockouts_per_method,
                "display_rank_viability_fraction": config.viability_fraction,
                "display_rank_criterion": (
                    "optimal and viability-qualified; descending product flux, descending "
                    "growth, blocked-reaction signature, target id; unique signatures"
                ),
                "selection_room_epsilon": room_epsilon,
                "selection_numeric_margin": config.product_improvement_tolerance,
                "selection_product_delta_threshold": room_improvement_threshold,
                "selection_product_delta_criterion": "strictly_greater_than_threshold",
            },
        }

    optknock_result = None
    robustknock_result = None
    if config.run_strain_design:
        optknock_result = optknock(
            working,
            config.product,
            biomass=biomass,
            max_knockouts=config.max_knockouts,
            max_solutions=config.optknock_max_solutions,
            min_growth=config.design_min_growth,
            actionable_only=config.actionable_designs_only,
            seed=config.strain_design_seed,
            max_seconds=config.strain_design_max_seconds,
        )
        robustknock_result = robustknock(
            working,
            config.product,
            biomass=biomass,
            max_knockouts=config.max_knockouts,
            max_solutions=config.robustknock_max_solutions,
            min_growth=config.design_min_growth,
            actionable_only=config.actionable_designs_only,
            seed=config.strain_design_seed,
            max_seconds=config.strain_design_max_seconds,
        )

    fseof_result = None
    fvseof_result = None
    if config.run_amplification:
        fseof_result = fseof(
            working,
            config.product,
            biomass,
            n_steps=config.fseof_steps,
            fraction_min=config.scan_fraction_min,
            fraction_max=config.scan_fraction_max,
        )
        fvseof_result = fvseof(
            working,
            config.product,
            biomass,
            n_steps=config.fvseof_steps,
            fraction_min=config.scan_fraction_min,
            fraction_max=config.scan_fraction_max,
            biomass_fraction=config.fvseof_biomass_fraction,
        )

    loop_diagnostic = None
    if (
        config.run_amplification
        and config.run_amplification_loop_diagnostic
        and (fseof_result is not None or fvseof_result is not None)
    ):
        loop_diagnostic = _run_amplification_loop_diagnostic(
            working,
            config,
            biomass,
            reference_growth,
            _amplification_candidate_sources(
                fseof_result,
                fvseof_result,
                per_method_limit=config.top_amplification_targets_per_method,
                total_limit=config.amplification_loop_diagnostic_top_n,
            ),
            enforced_product_floor=_diagnostic_product_floor(
                fseof_result, fvseof_result
            ),
        )

    validation_targets: tuple[ValidationTarget, ...] = ()
    response_results: tuple[FluxResponseValidation, ...] = ()
    sampling_results: tuple[SamplingValidation, ...] = ()
    if config.validation.enabled:
        validation_targets = _validation_targets(
            single_records,
            fseof_result,
            fvseof_result,
            loop_diagnostic,
            product=config.product,
            biomass=biomass,
            amplification_limit=config.top_amplification_targets_per_method,
            limit=config.validation.max_flux_response_targets,
        )
        response_results = _run_flux_response_validation(
            working,
            config,
            biomass,
            validation_targets,
            reference,
        )
        if config.validation.sampling.enabled and config.run_single_knockout:
            sampling_results = _run_sampling_validation(
                working,
                config,
                biomass,
                reference_growth,
                single_records,
                perturbations,
            )

    provenance: dict[str, object] = {
        **run_provenance(
            working,
            method="production_target_discovery",
            schema_id=SCHEMA_ID,
            schema_version=SCHEMA_VERSION,
            source_model_path=str(cast(Path, config.model_path)),
            source_model_sha256=original_fingerprint,
            product=config.product,
            substrate=yield_result.substrate,
            biomass=biomass,
            medium=medium_application.to_provenance(),
            condition=config.condition.name if config.condition else None,
            reference_method=config.reference_method,
            requested_solver=config.solver,
            strain_design_seed=config.strain_design_seed,
            strain_design_max_seconds=config.strain_design_max_seconds,
            validation_candidate_policy=dict(VALIDATION_CANDIDATE_POLICY),
        ),
        "source_model_sha256": original_fingerprint,
        "conditioned_model_sha256": model_fingerprint(working),
        "validation_candidate_policy": dict(VALIDATION_CANDIDATE_POLICY),
        "validation_coverage": _validation_coverage(
            validation_targets,
            response_results,
            sampling_results,
            sampling_enabled=(
                config.validation.enabled
                and config.run_single_knockout
                and config.validation.sampling.enabled
            ),
        ),
        "medium_application": medium_application.to_provenance(),
        "oxygen_uptake": oxygen_uptake,
        "analysis_metadata": {
            "wild_type": wild_type.metadata,
            "reference": dict(reference.metadata),
            "theoretical_yield": yield_result.metadata,
            "production_envelope": envelope.metadata,
            "single_knockout": single_metadata,
            "gene_knockout_mapping": gene_mapping.metadata,
            "optknock": optknock_result.metadata if optknock_result else None,
            "robustknock": robustknock_result.metadata if robustknock_result else None,
            "fseof": fseof_result.metadata if fseof_result else None,
            "fvseof": fvseof_result.metadata if fvseof_result else None,
            "amplification_loop_diagnostic": (
                loop_diagnostic.metadata if loop_diagnostic else None
            ),
            "flux_response": [
                item.result.metadata
                for item in response_results
                if item.result is not None
            ],
            "sampling": [
                item.result.metadata
                for item in sampling_results
                if item.result is not None
            ],
        },
    }
    preflight = _preflight_records(
        working,
        config,
        biomass,
        medium_application,
        reference,
        yield_result,
    )
    result = ProductionWorkflowResult(
        config=config,
        preflight=preflight,
        medium_application=medium_application,
        biomass=biomass,
        wild_type=wild_type,
        reference=reference,
        theoretical_yield=yield_result,
        envelope=envelope,
        gene_knockout_mapping=gene_mapping,
        single_knockouts=single_records,
        optknock_result=optknock_result,
        robustknock_result=robustknock_result,
        fseof_result=fseof_result,
        fvseof_result=fvseof_result,
        amplification_loop_diagnostic=loop_diagnostic,
        validation_targets=validation_targets,
        flux_responses=response_results,
        sampling=sampling_results,
        provenance=provenance,
    )
    if config.output_dir is not None:
        run_directory, artifacts = _export_result(
            result,
            working,
            model,
            cast(Path, config.output_dir),
        )
        result = replace(
            result,
            run_directory=run_directory,
            artifacts=artifacts,
        )
    return result


def _condition_model(
    model: Model, config: ProductionWorkflowConfig
) -> MediumApplication:
    if config.medium is None:
        application = MediumApplication(
            medium="model_as_loaded",
            applied=dict(model.medium),
            dropped=(),
        )
    else:
        application = apply_medium(model, config.medium)
    if config.condition is not None:
        config.condition.apply_to(model)
    return application


def _oxygen_exchange_uptake(model: Model) -> list[dict[str, object]]:
    """Describe enabled molecular-oxygen exchanges without relying on model-specific ids."""

    medium = dict(model.medium)
    records: list[dict[str, object]] = []
    for reaction in sorted(model.exchanges, key=lambda item: item.id):
        oxygen_metabolites = [
            metabolite
            for metabolite in reaction.metabolites
            if str(metabolite.formula or "").replace(" ", "").upper() == "O2"
            or str(metabolite.name or "").strip().lower()
            in {"o2", "oxygen", "molecular oxygen"}
        ]
        uptake_limit = float(medium.get(reaction.id, 0.0))
        if oxygen_metabolites and uptake_limit > 0.0:
            records.append(
                {
                    "reaction_id": reaction.id,
                    "metabolites": tuple(
                        sorted(metabolite.id for metabolite in oxygen_metabolites)
                    ),
                    "uptake_limit": uptake_limit,
                    "lower_bound": float(reaction.lower_bound),
                    "upper_bound": float(reaction.upper_bound),
                }
            )
    return records


def _preflight_model(model: Model, config: ProductionWorkflowConfig) -> str:
    exchanges = {reaction.id for reaction in model.exchanges}
    if not exchanges:
        raise ProductionWorkflowError(
            "model.exchanges is empty; production design requires explicit exchanges"
        )
    if config.product not in model.reactions:
        raise ProductionWorkflowError(
            f"product reaction {config.product!r} is not in the model"
        )
    if config.product not in exchanges:
        raise ProductionWorkflowError(
            f"product {config.product!r} is not an exchange reaction; do not substitute an "
            "internal reaction for production design"
        )
    if config.substrate is not None and config.substrate not in exchanges:
        raise ProductionWorkflowError(
            f"substrate {config.substrate!r} is not an exchange reaction"
        )
    biomass = _biomass_reaction(model, config.biomass)
    if config.run_single_knockout:
        require("QP", model.solver.interface, feature="single-gene L2 MOMA screen")
        require("MILP", model.solver.interface, feature="single-gene exact ROOM screen")
    if config.run_strain_design:
        require("MILP", model.solver.interface, feature="OptKnock/RobustKnock")
        try:
            import straindesign  # noqa: F401
        except ImportError as error:
            raise ProductionWorkflowError(
                "OptKnock/RobustKnock is enabled but 'straindesign' is not installed"
            ) from error
    return biomass


def _biomass_reaction(model: Model, explicit: str | None) -> str:
    if explicit is not None:
        if explicit not in model.reactions:
            raise ProductionWorkflowError(
                f"biomass reaction {explicit!r} is not in the model"
            )
        return explicit
    objective_reactions = sorted(
        reaction.id
        for reaction in model.reactions
        if abs(float(reaction.objective_coefficient)) > 0.0
    )
    if not objective_reactions:
        raise ProductionWorkflowError(
            "model has no objective reaction; set config.biomass explicitly"
        )
    if len(objective_reactions) > 1:
        raise ProductionWorkflowError(
            "model has a multi-reaction objective; set config.biomass explicitly so the "
            "growth-rate axis is unambiguous"
        )
    return objective_reactions[0]


def _preflight_records(
    model: Model,
    config: ProductionWorkflowConfig,
    biomass: str,
    medium_application: MediumApplication,
    reference: FluxState,
    yield_result: ProductionYield,
) -> tuple[PreflightRecord, ...]:
    status = solver_status(model)
    records = [
        PreflightRecord(
            "model_structure",
            "pass",
            {
                "reactions": len(model.reactions),
                "metabolites": len(model.metabolites),
                "genes": len(model.genes),
            },
            "model contains reactions, metabolites, genes, and an explicit objective",
        ),
        PreflightRecord(
            "exchange_reactions",
            "pass",
            len(model.exchanges),
            "model exposes exchange reactions",
        ),
        PreflightRecord(
            "product_exchange",
            "pass",
            config.product,
            "requested product is an explicit exchange reaction",
        ),
        PreflightRecord(
            "growth",
            "pass",
            reference.get(biomass),
            "reference state has positive biomass flux",
        ),
        PreflightRecord(
            "product_reachability",
            "pass",
            yield_result.product_flux,
            "the requested product has positive theoretical production",
        ),
        PreflightRecord(
            "solver",
            "pass",
            {"name": status.name, "capabilities": status.capabilities},
            "active solver satisfies every enabled workflow stage",
        ),
        PreflightRecord(
            "medium",
            "warning" if medium_application.medium == "model_as_loaded" else "pass",
            medium_application.to_provenance(),
            (
                "model bounds were retained; interpret results under the model-as-loaded "
                "medium"
                if medium_application.medium == "model_as_loaded"
                else "declared medium was applied and recorded"
            ),
        ),
    ]
    return tuple(records)


def _rank_single_knockout_screen(
    screen: BatchComparisonResult,
    perturbations: Sequence[Perturbation],
    *,
    method: SingleKnockoutMethod,
    wild_type_growth: float,
    wild_type_product: float,
    limit: int,
    viability_fraction: float,
    improvement_threshold: float,
) -> tuple[SingleKnockoutRecord, ...]:
    """Assign independent display and strict-benefit ranks to viable phenotype rows."""

    perturbation_by_id = {item.target_id: item for item in perturbations}
    base: list[SingleKnockoutRecord] = []
    for row in screen:
        perturbation = perturbation_by_id[row.target_id]
        growth_fraction = (
            row.objective / wild_type_growth
            if math.isfinite(row.objective) and wild_type_growth > 0
            else float("nan")
        )
        product_delta = (
            row.product_flux - wild_type_product
            if math.isfinite(row.product_flux)
            else float("nan")
        )
        fold_change = (
            row.product_flux / wild_type_product
            if math.isfinite(row.product_flux) and abs(wild_type_product) > 1e-12
            else None
        )
        signature = ";".join(perturbation.reaction_ids)
        base.append(
            SingleKnockoutRecord(
                method=method,
                target_id=row.target_id,
                blocked_reactions=perturbation.reaction_ids,
                blocked_reaction_signature=signature,
                status=row.status,
                growth_rate=row.objective,
                growth_fraction=growth_fraction,
                target_production=row.product_flux,
                product_delta=product_delta,
                product_fold_change=fold_change,
                objective_value=row.objective_value,
                distance=row.distance,
                distance_kind=row.distance_kind,
                n_changed_reactions=row.n_changed_reactions,
                improves_product=(
                    math.isfinite(product_delta)
                    and product_delta > improvement_threshold
                ),
            )
        )

    display_eligible = [
        record
        for record in base
        if record.status == "optimal"
        and math.isfinite(record.growth_rate)
        and math.isfinite(record.target_production)
        and record.growth_fraction >= viability_fraction
    ]
    display_eligible.sort(
        key=lambda record: (
            -record.target_production,
            -record.growth_rate,
            record.blocked_reaction_signature,
            record.target_id,
        )
    )
    display_rank_by_target: dict[str, int] = {}
    seen_signatures: set[str] = set()
    for record in display_eligible:
        if record.blocked_reaction_signature in seen_signatures:
            continue
        seen_signatures.add(record.blocked_reaction_signature)
        display_rank_by_target[record.target_id] = len(display_rank_by_target) + 1
        if len(display_rank_by_target) == limit:
            break

    rank_by_target: dict[str, int] = {}
    seen_signatures = set()
    for record in display_eligible:
        if record.product_delta <= improvement_threshold:
            continue
        if record.blocked_reaction_signature in seen_signatures:
            continue
        seen_signatures.add(record.blocked_reaction_signature)
        rank_by_target[record.target_id] = len(rank_by_target) + 1
        if len(rank_by_target) == limit:
            break

    return tuple(
        replace(
            record,
            selected=record.target_id in rank_by_target,
            method_rank=rank_by_target.get(record.target_id),
            display_rank=display_rank_by_target.get(record.target_id),
        )
        for record in sorted(base, key=lambda item: item.target_id)
    )


def _gene_knockout_mapping(model: Model) -> GeneKnockoutMappingResult:
    records: list[GeneKnockoutMappingRecord] = []
    inert_genes = 0
    for gene in sorted(model.genes, key=lambda item: item.id):
        blocked = blocked_reactions_for_genes(model, (gene.id,))
        if not blocked:
            inert_genes += 1
            associated = sorted(gene.reactions, key=lambda item: item.id)
            records.append(
                GeneKnockoutMappingRecord(
                    gene_id=gene.id,
                    gene_name=str(gene.name or ""),
                    inert=True,
                    blocked_reaction=None,
                    reaction_name=(
                        "; ".join(
                            f"{reaction.id}: {reaction.name or ''}"
                            for reaction in associated
                        )
                        or None
                    ),
                    reaction_equation=(
                        "; ".join(
                            f"{reaction.id}: "
                            f"{reaction.build_reaction_string(use_metabolite_names=False)}"
                            for reaction in associated
                        )
                        or None
                    ),
                    gpr=(
                        "; ".join(
                            f"{reaction.id}: {reaction.gene_reaction_rule}"
                            for reaction in associated
                        )
                        or None
                    ),
                )
            )
            continue
        for reaction_id in blocked:
            reaction = model.reactions.get_by_id(reaction_id)
            records.append(
                GeneKnockoutMappingRecord(
                    gene_id=gene.id,
                    gene_name=str(gene.name or ""),
                    inert=False,
                    blocked_reaction=reaction_id,
                    reaction_name=str(reaction.name or ""),
                    reaction_equation=reaction.build_reaction_string(
                        use_metabolite_names=False
                    ),
                    gpr=str(reaction.gene_reaction_rule),
                )
            )
    metadata = run_provenance(
        model,
        method="gene_knockout_mapping",
        n_genes=len(model.genes),
        n_inert_genes=inert_genes,
        n_mapping_rows=len(records),
        inert_definition="single_gene_deletion_blocks_no_reaction_under_model_gpr",
    )
    return GeneKnockoutMappingResult(tuple(records), metadata)


def _amplification_candidate_sources(
    fseof_result: FseofResult | None,
    fvseof_result: FvseofResult | None,
    *,
    per_method_limit: int,
    total_limit: int,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    sources: dict[str, set[str]] = {}
    method_rank: dict[str, int] = {}
    for method, targets in (
        (
            "fseof",
            fseof_result.amplification_targets()[:per_method_limit]
            if fseof_result is not None
            else [],
        ),
        (
            "fvseof",
            fvseof_result.amplification_targets()[:per_method_limit]
            if fvseof_result is not None
            else [],
        ),
    ):
        for rank, reaction_id in enumerate(targets, start=1):
            sources.setdefault(reaction_id, set()).add(method)
            method_rank[reaction_id] = min(method_rank.get(reaction_id, rank), rank)
    ordered = sorted(
        sources,
        key=lambda reaction_id: (
            method_rank[reaction_id],
            reaction_id,
        ),
    )
    return tuple(
        (reaction_id, tuple(sorted(sources[reaction_id])))
        for reaction_id in ordered[:total_limit]
    )


def _run_amplification_loop_diagnostic(
    model: Model,
    config: ProductionWorkflowConfig,
    biomass: str,
    wild_type_growth: float,
    candidates: Sequence[tuple[str, tuple[str, ...]]],
    *,
    enforced_product_floor: float,
) -> AmplificationLoopDiagnosticResult:
    metadata = run_provenance(
        model,
        method="amplification_loop_diagnostic",
        targets=tuple(reaction_id for reaction_id, _ in candidates),
        standard_loopless=False,
        loopless_algorithm=config.loopless_algorithm,
        capacity_ratio_threshold=config.loopless_capacity_ratio_threshold,
        ratio_definition="loopless_capacity_divided_by_standard_capacity",
        artifact_rule="ratio_below_threshold",
        requested_enforced_product_floor=enforced_product_floor,
        requested_biomass_floor=(
            config.validation.flux_response_biomass_fraction * wild_type_growth
        ),
    )
    if not candidates:
        return AmplificationLoopDiagnosticResult((), metadata)

    target_ids = tuple(reaction_id for reaction_id, _ in candidates)
    actual_product_floor = enforced_product_floor
    actual_biomass_floor = 0.0
    try:
        with model:
            product_reaction = model.reactions.get_by_id(config.product)
            actual_product_floor = max(
                float(product_reaction.lower_bound), enforced_product_floor
            )
            product_reaction.bounds = (
                actual_product_floor,
                float(product_reaction.upper_bound),
            )
            biomass_reaction = model.reactions.get_by_id(biomass)
            model.objective = biomass_reaction
            model.objective_direction = "max"
            maximum_growth = model.slim_optimize(error_value=float("nan"))
            if maximum_growth is None or not math.isfinite(float(maximum_growth)):
                raise ProductionWorkflowError(
                    "product-enforced model is infeasible for loopless FVA diagnostic"
                )
            requested_floor = (
                config.validation.flux_response_biomass_fraction * wild_type_growth
            )
            actual_biomass_floor = max(
                float(biomass_reaction.lower_bound), requested_floor
            )
            if float(maximum_growth) + 1e-9 < actual_biomass_floor:
                raise ProductionWorkflowError(
                    "the configured biomass floor is infeasible at the diagnostic product "
                    f"floor (maximum growth={float(maximum_growth):.6g}, requested "
                    f"floor={actual_biomass_floor:.6g})"
                )
            biomass_reaction.bounds = (
                actual_biomass_floor,
                float(biomass_reaction.upper_bound),
            )
            standard = fva(
                model,
                reactions=target_ids,
                fraction_of_optimum=0.0,
                loopless=False,
                processes=1,
            )
            loopless = fva(
                model,
                reactions=target_ids,
                fraction_of_optimum=0.0,
                loopless=config.loopless_algorithm,
                processes=1,
            )
    except Exception as error:
        reason = f"{type(error).__name__}: {error}"
        records = tuple(
            AmplificationLoopDiagnosticRecord(
                rank=rank,
                target=reaction_id,
                source_methods=methods,
                standard_minimum=None,
                standard_maximum=None,
                standard_capacity=None,
                loopless_minimum=None,
                loopless_maximum=None,
                loopless_capacity=None,
                loopless_to_standard_capacity_ratio=None,
                capacity_ratio_threshold=config.loopless_capacity_ratio_threshold,
                loop_artifact_flag=None,
                diagnostic_status="failed",
                reason=reason,
                enforced_product_floor=actual_product_floor,
                biomass_floor=actual_biomass_floor,
            )
            for rank, (reaction_id, methods) in enumerate(candidates, start=1)
        )
        return AmplificationLoopDiagnosticResult(records, metadata)

    records = tuple(
        _amplification_loop_record(
            rank,
            reaction_id,
            methods,
            standard[reaction_id],
            loopless[reaction_id],
            threshold=config.loopless_capacity_ratio_threshold,
            product_floor=actual_product_floor,
            biomass_floor=actual_biomass_floor,
        )
        for rank, (reaction_id, methods) in enumerate(candidates, start=1)
    )
    metadata = {
        **metadata,
        "standard_fva": standard.metadata,
        "loopless_fva": loopless.metadata,
        "enforced_product_floor": actual_product_floor,
        "actual_biomass_floor": actual_biomass_floor,
        "n_flagged": sum(record.loop_artifact_flag is True for record in records),
    }
    return AmplificationLoopDiagnosticResult(records, metadata)


def _diagnostic_product_floor(
    fseof_result: FseofResult | None,
    fvseof_result: FvseofResult | None,
) -> float:
    minimum_levels = [
        min(float(level) for level in result.enforced_levels)
        for result in (fseof_result, fvseof_result)
        if result is not None and result.enforced_levels
    ]
    # The highest method-specific minimum is present in both scan domains when their grids
    # differ.  It retains the declared 30% biomass floor without silently relaxing either
    # constraint at the high-product end of a scan.
    return max(minimum_levels, default=0.0)


def _amplification_loop_record(
    rank: int,
    reaction_id: str,
    methods: tuple[str, ...],
    standard: FluxRange,
    loopless: FluxRange,
    *,
    threshold: float,
    product_floor: float,
    biomass_floor: float,
) -> AmplificationLoopDiagnosticRecord:
    def finite_capacity(flux_range: FluxRange) -> float | None:
        if not (
            math.isfinite(flux_range.minimum) and math.isfinite(flux_range.maximum)
        ):
            return None
        capacity = max(0.0, flux_range.maximum - flux_range.minimum)
        return capacity if math.isfinite(capacity) else None

    standard_capacity = finite_capacity(standard)
    loopless_capacity = finite_capacity(loopless)
    ratio: float | None
    flag: bool | None
    reason: str | None
    if standard_capacity is None or loopless_capacity is None:
        ratio = None
        flag = None
        status: Literal["complete", "inconclusive", "failed"] = "inconclusive"
        non_finite = []
        if standard_capacity is None:
            non_finite.append("standard")
        if loopless_capacity is None:
            non_finite.append("loopless")
        reason = f"non-finite {' and '.join(non_finite)} FVA bounds or capacity"
    elif standard_capacity <= 1e-9:
        ratio = None
        flag = None
        status = "inconclusive"
        reason = "standard FVA capacity is zero; a capacity ratio is undefined"
    else:
        ratio = loopless_capacity / standard_capacity
        flag = ratio < threshold
        status = "complete"
        reason = (
            "loopless capacity ratio is below the declared artifact threshold"
            if flag
            else None
        )
    return AmplificationLoopDiagnosticRecord(
        rank=rank,
        target=reaction_id,
        source_methods=methods,
        standard_minimum=(
            standard.minimum if math.isfinite(standard.minimum) else None
        ),
        standard_maximum=(
            standard.maximum if math.isfinite(standard.maximum) else None
        ),
        standard_capacity=standard_capacity,
        loopless_minimum=(
            loopless.minimum if math.isfinite(loopless.minimum) else None
        ),
        loopless_maximum=(
            loopless.maximum if math.isfinite(loopless.maximum) else None
        ),
        loopless_capacity=loopless_capacity,
        loopless_to_standard_capacity_ratio=ratio,
        capacity_ratio_threshold=threshold,
        loop_artifact_flag=flag,
        diagnostic_status=status,
        reason=reason,
        enforced_product_floor=product_floor,
        biomass_floor=biomass_floor,
    )


@dataclass(frozen=True)
class _SingleKnockoutValidationGroup:
    """One display-selected phenotype plus every equivalent screen-row alias."""

    display_records: tuple[SingleKnockoutRecord, ...]
    alias_records: tuple[SingleKnockoutRecord, ...]

    @property
    def representative(self) -> SingleKnockoutRecord:
        """Return a deterministic representative drawn only from D1-D5 rows."""

        return min(
            self.display_records,
            key=lambda item: (
                item.target_id,
                item.method,
                item.display_rank if item.display_rank is not None else 0,
            ),
        )

    @property
    def candidate_target_ids(self) -> tuple[str, ...]:
        """Return all gene ids whose knockout has this model phenotype."""

        return tuple(sorted({item.target_id for item in self.alias_records}))

    @property
    def source_methods(self) -> tuple[str, ...]:
        """Return only methods that placed the phenotype in their D1-D5 list."""

        return tuple(sorted({item.method for item in self.display_records}))


def _single_knockout_validation_groups(
    records: Sequence[SingleKnockoutRecord],
) -> tuple[_SingleKnockoutValidationGroup, ...]:
    """Return the complete D1-D5 validation universe grouped by GPR phenotype.

    ``display_rank`` defines a method-specific candidate independently of the stricter
    beneficial/recommendation rank. Candidate signatures come only from those D1-D5 rows,
    while alias provenance comes from every MOMA/ROOM screen row with the same blocked-
    reaction signature. Validation therefore runs once per model phenotype without losing
    non-representative gene ids from its index.
    """

    display_by_signature: dict[str, list[SingleKnockoutRecord]] = {}
    aliases_by_signature: dict[str, list[SingleKnockoutRecord]] = {}
    for record in records:
        aliases_by_signature.setdefault(record.blocked_reaction_signature, []).append(
            record
        )
        if record.display_rank is not None:
            display_by_signature.setdefault(
                record.blocked_reaction_signature, []
            ).append(record)

    groups: list[_SingleKnockoutValidationGroup] = []
    for signature, display_records in display_by_signature.items():
        groups.append(
            _SingleKnockoutValidationGroup(
                display_records=tuple(
                    sorted(
                        display_records,
                        key=lambda item: (
                            item.method,
                            item.display_rank if item.display_rank is not None else 0,
                            item.target_id,
                        ),
                    )
                ),
                alias_records=tuple(
                    sorted(
                        aliases_by_signature[signature],
                        key=lambda item: (
                            item.target_id,
                            item.method,
                            item.display_rank if item.display_rank is not None else 0,
                        ),
                    )
                ),
            )
        )
    groups.sort(
        key=lambda group: (
            min(item.display_rank or 0 for item in group.display_records),
            group.representative.target_id,
            group.representative.blocked_reaction_signature,
        )
    )
    return tuple(groups)


def _validation_targets(
    single_records: Sequence[SingleKnockoutRecord],
    fseof_result: FseofResult | None,
    fvseof_result: FvseofResult | None,
    loop_diagnostic: AmplificationLoopDiagnosticResult | None,
    *,
    product: str,
    biomass: str,
    amplification_limit: int,
    limit: int,
) -> tuple[ValidationTarget, ...]:
    """Build the complete canonical forward-validation universe.

    Amplification candidates are the independent report-selected FSEOF and FVSEOF lists.
    Every candidate receives a response scan even when its loop diagnostic is flagged or
    unresolved; the diagnostic fields prevent those rows from becoming recommendations.
    Single-knockout candidates are every MOMA/ROOM ``display_rank`` row, deduplicated only
    when multiple genes have the same blocked-reaction signature.
    """

    collected: dict[str, tuple[set[str], set[str]]] = {}
    amplification_order: dict[str, int] = {}
    diagnostic_by_target = {
        record.target: record
        for record in (loop_diagnostic.records if loop_diagnostic is not None else ())
    }

    def add(reaction_id: str, action: str, method: str) -> None:
        actions, methods = collected.setdefault(reaction_id, (set(), set()))
        actions.add(action)
        methods.add(method)

    if fseof_result is not None:
        for rank, reaction_id in enumerate(
            fseof_result.amplification_targets()[:amplification_limit], start=1
        ):
            add(reaction_id, "amplify", "fseof")
            amplification_order[reaction_id] = min(
                amplification_order.get(reaction_id, rank), rank
            )
    if fvseof_result is not None:
        for rank, reaction_id in enumerate(
            fvseof_result.amplification_targets()[:amplification_limit], start=1
        ):
            add(reaction_id, "amplify", "fvseof")
            amplification_order[reaction_id] = min(
                amplification_order.get(reaction_id, rank), rank
            )
    amplification_targets: list[ValidationTarget] = []
    for reaction_id, (actions, methods) in collected.items():
        diagnostic = diagnostic_by_target.get(reaction_id)
        diagnostic_status: str
        artifact_flag: bool | None
        diagnostic_eligible: bool
        diagnostic_reason: str | None
        if diagnostic is None:
            diagnostic_status = "unavailable"
            artifact_flag = None
            diagnostic_eligible = False
            diagnostic_reason = (
                "loop diagnostic unavailable; flux response is retained as exploratory "
                "evidence but cannot support a recommendation"
            )
        else:
            diagnostic_status = diagnostic.diagnostic_status
            artifact_flag = diagnostic.loop_artifact_flag
            diagnostic_eligible = (
                diagnostic.diagnostic_status == "complete"
                and diagnostic.loop_artifact_flag is False
            )
            diagnostic_reason = diagnostic.reason
            if diagnostic_reason is None and diagnostic.loop_artifact_flag is True:
                diagnostic_reason = "loop diagnostic flagged a capacity artifact"
            elif diagnostic_reason is None and not diagnostic_eligible:
                diagnostic_reason = "loop diagnostic was not conclusively eligible"
        amplification_targets.append(
            ValidationTarget(
                target_id=reaction_id,
                scan_reaction=reaction_id,
                response_reaction=product,
                background="wild_type",
                actions=tuple(sorted(actions)),
                source_methods=tuple(sorted(methods)),
                candidate_scope="all_report_selected_candidates",
                candidate_target_ids=(reaction_id,),
                loop_diagnostic_status=diagnostic_status,
                loop_artifact_flag=artifact_flag,
                loop_diagnostic_eligible=diagnostic_eligible,
                loop_diagnostic_reason=diagnostic_reason,
            )
        )
    amplification_targets.sort(
        key=lambda item: (amplification_order[item.target_id], item.target_id)
    )

    knockout_targets: list[ValidationTarget] = []
    for group in _single_knockout_validation_groups(single_records):
        representative = group.representative
        blocked_reactions = representative.blocked_reactions
        scan_reaction = (
            blocked_reactions[0]
            if len(blocked_reactions) == 1
            else representative.blocked_reaction_signature
        )
        knockout_targets.append(
            ValidationTarget(
                target_id=representative.target_id,
                scan_reaction=scan_reaction,
                response_reaction=product,
                background="wild_type",
                actions=("knockout",),
                source_methods=group.source_methods,
                blocked_reactions=blocked_reactions,
                candidate_scope="all_display_ranked_candidates",
                blocked_reaction_signature=representative.blocked_reaction_signature,
                candidate_target_ids=group.candidate_target_ids,
            )
        )
    targets = amplification_targets + knockout_targets
    if len(targets) > limit:
        raise ProductionWorkflowError(
            "validation candidate universe exceeds max_flux_response_targets: "
            f"{len(targets)} candidates require capacity but the configured limit is "
            f"{limit}; increase the explicit capacity instead of truncating candidates"
        )
    return tuple(targets)


def _run_flux_response_validation(
    model: Model,
    config: ProductionWorkflowConfig,
    biomass: str,
    targets: Sequence[ValidationTarget],
    reference: FluxState,
) -> tuple[FluxResponseValidation, ...]:
    """Scan candidate-reaction flux against product production under a growth floor.

    Amplification candidates use their full feasible target-reaction domain.  A
    single-reaction knockout candidate is evaluated before deletion by titrating the
    corresponding reaction between its reproducible reference flux and zero.  When that
    reference is already zero, the full feasible reaction domain is scanned so the candidate
    is still covered, but the curve cannot support a causal deletion claim.  Applying the
    deletion first would fix the reaction at zero and make a response curve impossible; the
    complete-deletion phenotype is instead covered by MOMA/ROOM and matched sampling.
    """

    results: list[FluxResponseValidation] = []
    for target in targets:
        is_knockout_candidate = (
            target.candidate_scope == "all_display_ranked_candidates"
            or (not target.candidate_scope and target.background == "gene_knockout")
        )
        executed_target = target
        target_min: float | None = None
        target_max: float | None = None
        if is_knockout_candidate:
            if len(target.blocked_reactions) != 1:
                results.append(
                    FluxResponseValidation(
                        target=target,
                        status="skipped",
                        reason=(
                            "single-axis flux response cannot represent a knockout that "
                            "simultaneously blocks multiple reactions; no reaction was "
                            "chosen silently"
                        ),
                    )
                )
                continue
            reaction_id = target.blocked_reactions[0]
            try:
                reference_flux_value = float(reference.get(reaction_id))
            except (KeyError, TypeError, ValueError) as error:
                results.append(
                    FluxResponseValidation(
                        target=target,
                        status="skipped",
                        reason=(
                            f"reference flux for {reaction_id!r} is unavailable "
                            f"({type(error).__name__}: {error})"
                        ),
                    )
                )
                continue
            executed_target = replace(
                target,
                scan_reaction=reaction_id,
                response_reaction=config.product,
                background="wild_type",
                scan_reference_flux=reference_flux_value,
            )
            if not math.isfinite(reference_flux_value):
                results.append(
                    FluxResponseValidation(
                        target=executed_target,
                        status="skipped",
                        reason=(
                            f"reference flux for {reaction_id!r} is not finite; knockout "
                            "titration is unavailable"
                        ),
                    )
                )
                continue
            if abs(reference_flux_value) > 1e-9:
                target_min = min(0.0, reference_flux_value)
                target_max = max(0.0, reference_flux_value)
        try:
            with model:
                response = flux_response(
                    model,
                    executed_target.scan_reaction,
                    executed_target.response_reaction,
                    biomass=biomass,
                    target_min=target_min,
                    target_max=target_max,
                    n_steps=config.validation.flux_response_steps,
                    biomass_fraction=config.validation.flux_response_biomass_fraction,
                )
            results.append(
                FluxResponseValidation(
                    target=executed_target,
                    status="complete",
                    result=response,
                )
            )
        except (
            Exception
        ) as error:  # each failed target remains an explicit validation row
            results.append(
                FluxResponseValidation(
                    target=executed_target,
                    status="failed",
                    error=f"{type(error).__name__}: {error}",
                )
            )
    return tuple(results)


def _run_sampling_validation(
    model: Model,
    config: ProductionWorkflowConfig,
    biomass: str,
    wild_type_growth: float,
    records: Sequence[SingleKnockoutRecord],
    perturbations: Sequence[Perturbation],
) -> tuple[SamplingValidation, ...]:
    candidate_groups = _single_knockout_validation_groups(records)
    if not candidate_groups:
        return ()

    sampling = config.validation.sampling
    floor = wild_type_growth * config.validation.sampling_growth_fraction
    try:
        with model:
            biomass_reaction = model.reactions.get_by_id(biomass)
            biomass_reaction.bounds = (
                max(float(biomass_reaction.lower_bound), floor),
                float(biomass_reaction.upper_bound),
            )
            wild_type = random_flux_sampling(
                model,
                n=sampling.n,
                method=sampling.method,
                thinning=sampling.thinning,
                processes=sampling.processes,
                seed=sampling.seed,
            )
    except Exception as error:
        error_text = f"{type(error).__name__}: {error}"
        unavailable: list[SamplingValidation] = [
            SamplingValidation(
                target_id="wild_type",
                blocked_reactions=(),
                source_methods=(),
                status="failed",
                error=error_text,
                candidate_scope="reference",
            )
        ]
        for group in candidate_groups:
            representative = group.representative
            unavailable.append(
                SamplingValidation(
                    target_id=representative.target_id,
                    blocked_reactions=representative.blocked_reactions,
                    source_methods=group.source_methods,
                    status="skipped",
                    reason=(
                        "wild-type sampling failed, so the paired knockout distribution "
                        f"was unavailable: {error_text}"
                    ),
                    candidate_scope="all_display_ranked_candidates",
                    blocked_reaction_signature=(
                        representative.blocked_reaction_signature
                    ),
                    candidate_target_ids=group.candidate_target_ids,
                )
            )
        return tuple(unavailable)

    output: list[SamplingValidation] = [
        SamplingValidation(
            target_id="wild_type",
            blocked_reactions=(),
            source_methods=(),
            status="complete",
            result=wild_type,
            candidate_scope="reference",
        )
    ]
    perturbation_by_id = {item.target_id: item for item in perturbations}
    for group in candidate_groups:
        representative = group.representative
        methods = group.source_methods
        candidate_target_ids = group.candidate_target_ids
        perturbation = perturbation_by_id.get(representative.target_id)
        if perturbation is None:
            output.append(
                SamplingValidation(
                    target_id=representative.target_id,
                    blocked_reactions=representative.blocked_reactions,
                    source_methods=methods,
                    status="skipped",
                    reason="canonical knockout perturbation was unavailable",
                    candidate_scope="all_display_ranked_candidates",
                    blocked_reaction_signature=(
                        representative.blocked_reaction_signature
                    ),
                    candidate_target_ids=candidate_target_ids,
                )
            )
            continue
        try:
            with apply_perturbation(model, perturbation):
                biomass_reaction = model.reactions.get_by_id(biomass)
                biomass_reaction.bounds = (
                    max(float(biomass_reaction.lower_bound), floor),
                    float(biomass_reaction.upper_bound),
                )
                mutant = random_flux_sampling(
                    model,
                    n=sampling.n,
                    method=sampling.method,
                    thinning=sampling.thinning,
                    processes=sampling.processes,
                    seed=sampling.seed,
                )
            output.append(
                SamplingValidation(
                    target_id=representative.target_id,
                    blocked_reactions=perturbation.reaction_ids,
                    source_methods=methods,
                    status="complete",
                    result=mutant,
                    comparison=_sampling_comparison(wild_type, mutant),
                    candidate_scope="all_display_ranked_candidates",
                    blocked_reaction_signature=(
                        representative.blocked_reaction_signature
                    ),
                    candidate_target_ids=candidate_target_ids,
                )
            )
        except Exception as error:
            output.append(
                SamplingValidation(
                    target_id=representative.target_id,
                    blocked_reactions=perturbation.reaction_ids,
                    source_methods=methods,
                    status="failed",
                    error=f"{type(error).__name__}: {error}",
                    candidate_scope="all_display_ranked_candidates",
                    blocked_reaction_signature=(
                        representative.blocked_reaction_signature
                    ),
                    candidate_target_ids=candidate_target_ids,
                )
            )
    return tuple(output)


def _sampling_comparison(
    wild_type: SamplingResult,
    mutant: SamplingResult,
) -> pd.DataFrame:
    wild = wild_type.statistics().add_prefix("wild_type_")
    changed = mutant.statistics().add_prefix("knockout_")
    comparison = wild.join(changed, how="outer")
    comparison["mean_delta"] = (
        comparison["knockout_mean"] - comparison["wild_type_mean"]
    )
    comparison["median_delta"] = (
        comparison["knockout_median"] - comparison["wild_type_median"]
    )
    pooled = np.sqrt(
        (comparison["wild_type_std"] ** 2 + comparison["knockout_std"] ** 2) / 2.0
    )
    comparison["standardized_mean_shift"] = comparison["mean_delta"] / pooled.replace(
        0.0, np.nan
    )
    shared_reactions = comparison.index.intersection(
        wild_type.samples.columns
    ).intersection(mutant.samples.columns)
    comparison["wasserstein_distance"] = np.nan
    for reaction_id in shared_reactions:
        comparison.loc[reaction_id, "wasserstein_distance"] = wasserstein_distance(
            wild_type.samples[reaction_id].to_numpy(dtype=float),
            mutant.samples[reaction_id].to_numpy(dtype=float),
        )
    comparison.index.name = "reaction_id"
    return comparison.reset_index()


def _publication_amplification_targets(
    result: ProductionWorkflowResult,
    raw_targets: Sequence[str],
) -> tuple[tuple[str, ...], ArtifactStatus, str | None]:
    """Return a method's independent report shortlist with diagnostic status.

    The report retains explicitly flagged targets so it can disclose and visually segregate
    them.  Main flux-response validation and recommendations remain stricter and admit only
    complete, non-flagged diagnostics.
    """

    if not result.config.run_amplification:
        return (), "skipped", "amplification stage was disabled in the workflow config"
    if not raw_targets:
        return (), "complete", "the method returned no actionable amplification targets"
    diagnostic = result.amplification_loop_diagnostic
    if diagnostic is None:
        return (
            tuple(raw_targets),
            "partial",
            "loopless FVA diagnostic was unavailable; report targets retain an unknown "
            "artifact status",
        )

    by_target = {record.target: record for record in diagnostic.records}
    flagged: list[str] = []
    unresolved: list[str] = []
    for target in raw_targets:
        record = by_target.get(target)
        if record is None or record.diagnostic_status != "complete":
            unresolved.append(target)
        elif record.loop_artifact_flag is True:
            flagged.append(target)
        elif record.loop_artifact_flag is not False:
            unresolved.append(target)

    reasons: list[str] = []
    if flagged:
        reasons.append(
            f"retained {len(flagged)} flagged candidate(s) as diagnostic-only: "
            + ";".join(flagged)
        )
    if unresolved:
        reasons.append(
            f"retained {len(unresolved)} candidate(s) with unresolved loopless status: "
            + ";".join(unresolved)
        )
    return (
        tuple(raw_targets),
        "partial" if unresolved else "complete",
        "; ".join(reasons) or None,
    )


def _annotate_loop_diagnostic(
    frame: pd.DataFrame,
    result: ProductionWorkflowResult,
    *,
    target_column: str = "target",
) -> pd.DataFrame:
    """Attach loop evidence without conflating display and recommendation eligibility."""

    output = frame.copy()
    by_target = {
        record.target: record
        for record in (
            result.amplification_loop_diagnostic.records
            if result.amplification_loop_diagnostic is not None
            else ()
        )
    }

    def value(target: object, field: str) -> object:
        record = by_target.get(str(target))
        return getattr(record, field) if record is not None else None

    output["loop_diagnostic_status"] = output[target_column].map(
        lambda target: value(target, "diagnostic_status")
    )
    output["loop_artifact_flag"] = output[target_column].map(
        lambda target: value(target, "loop_artifact_flag")
    )
    output["loopless_to_standard_capacity_ratio"] = output[target_column].map(
        lambda target: value(target, "loopless_to_standard_capacity_ratio")
    )
    output["loop_diagnostic_reason"] = output[target_column].map(
        lambda target: value(target, "reason")
    )
    return output


def _export_result(
    result: ProductionWorkflowResult,
    conditioned_model: Model,
    source_model: Model,
    output_dir: Path,
) -> tuple[Path, tuple[ArtifactRecord, ...]]:
    root = output_dir.expanduser().resolve()
    configured_source = cast(Path, result.config.model_path).expanduser().resolve()
    # Capture the authoritative bytes before overwrite cleanup: a legitimate in-place rerun
    # may point model_path at this bundle's own model/ copy.
    configured_source_bytes = (
        configured_source.read_bytes() if configured_source.is_file() else None
    )
    if root.exists() and not root.is_dir():
        raise FileExistsError(f"output path exists and is not a directory: {root}")
    if root.exists() and any(root.iterdir()) and not result.config.overwrite:
        raise FileExistsError(
            f"output directory is not empty: {root}; choose a new directory or set "
            "overwrite=True"
        )
    if result.config.overwrite:
        # Replace only CMM-owned paths.  This removes stale target files and publication
        # products without following a malicious/stale stage-directory symlink or deleting
        # unrelated user files that happen to share the requested root.
        for relative in (
            *_STAGE_DIRECTORIES,
            "00_config.json",
            "00_provenance.json",
            "00_summary.json",
            "00_manifest.json",
            "report.html",
            "report_standalone.html",
            "report_validation.json",
        ):
            generated = root / relative
            if generated.is_symlink() or generated.is_file():
                generated.unlink()
            elif generated.is_dir():
                shutil.rmtree(generated)
    root.mkdir(parents=True, exist_ok=True)
    for directory in _STAGE_DIRECTORIES:
        (root / directory).mkdir(parents=True, exist_ok=True)

    writer = _ArtifactWriter(root)
    writer.csv(
        "01_preflight/preflight.csv",
        pd.DataFrame([asdict(record) for record in result.preflight]),
        stage="01_preflight",
        role="preflight_checks",
    )
    writer.csv(
        "02_yield/theoretical_yield.csv",
        result.theoretical_yield.to_frame(),
        stage="02_yield",
        role="theoretical_yield",
        method="theoretical_yield",
    )
    writer.csv(
        "02_yield/production_envelope.csv",
        result.envelope.to_frame(),
        stage="02_yield",
        role="production_envelope",
        method="production_envelope",
    )
    writer.csv(
        "03_reference/wild_type_fluxes.csv",
        result.wild_type.to_frame(),
        stage="03_reference",
        role="wild_type_fluxes",
        method="fba",
    )
    writer.csv(
        "03_reference/reference_fluxes.csv",
        _flux_frame(result.reference.fluxes),
        stage="03_reference",
        role="reference_fluxes",
        method=result.config.reference_method,
    )
    writer.csv(
        "04_single_knockout/gene_knockout_mapping.csv",
        result.gene_knockout_mapping.to_frame(),
        stage="04_single_knockout",
        role="gene_knockout_mapping",
        method="gpr_mapping",
    )

    single_columns = list(
        _single_knockout_export_row(
            SingleKnockoutRecord(
                method="moma_l2",
                target_id="",
                blocked_reactions=(),
                blocked_reaction_signature="",
                status="",
                growth_rate=float("nan"),
                growth_fraction=float("nan"),
                target_production=float("nan"),
                product_delta=float("nan"),
                product_fold_change=None,
                objective_value=float("nan"),
                distance=None,
                distance_kind="",
                n_changed_reactions=None,
            ),
            conditioned_model,
        )
    )
    for method, filename in (
        ("moma_l2", "single_knockout_moma.csv"),
        ("room", "single_knockout_room.csv"),
    ):
        rows = [
            _single_knockout_export_row(record, conditioned_model)
            for record in result.single_knockouts
            if record.method == method
        ]
        writer.csv(
            f"04_single_knockout/{filename}",
            pd.DataFrame(rows, columns=single_columns),
            stage="04_single_knockout",
            role=f"single_knockout_{'moma' if method == 'moma_l2' else 'room'}",
            method=method,
            status="complete" if result.config.run_single_knockout else "skipped",
            reason=(
                None
                if result.config.run_single_knockout
                else "single-knockout stage was disabled in the workflow config"
            ),
        )
    candidate_rows = _single_knockout_candidate_export_rows(
        result.single_knockouts,
        conditioned_model,
    )
    candidate_columns = [
        *single_columns,
        "candidate_scope",
        "validation_target_id",
        "candidate_target_ids",
        "candidate_source_methods",
        "validation_representative",
    ]
    writer.csv(
        "04_single_knockout/single_knockout_candidates.csv",
        pd.DataFrame(
            candidate_rows,
            columns=candidate_columns,
        ),
        stage="04_single_knockout",
        role="ranked_unique_single_gene_candidates",
        method="moma_l2+room",
        status="complete" if result.config.run_single_knockout else "skipped",
        reason=(
            None
            if result.config.run_single_knockout
            else "single-knockout stage was disabled in the workflow config"
        ),
    )
    by_target: dict[str, list[SingleKnockoutRecord]] = {}
    for record in result.single_knockouts:
        by_target.setdefault(record.target_id, []).append(record)
    consensus = pd.DataFrame(
        [
            {
                "target_id": target_id,
                "recommended": any(record.selected for record in records),
                "selected_methods": ";".join(
                    sorted(record.method for record in records if record.selected)
                ),
                "n_methods_selected": sum(record.selected for record in records),
            }
            for target_id, records in sorted(by_target.items())
        ],
        columns=[
            "target_id",
            "recommended",
            "selected_methods",
            "n_methods_selected",
        ],
    )
    writer.csv(
        "04_single_knockout/single_knockout_consensus.csv",
        consensus,
        stage="04_single_knockout",
        role="single_knockout_consensus",
        method="moma_l2+room",
        status="complete" if result.config.run_single_knockout else "skipped",
        reason=(
            None
            if result.config.run_single_knockout
            else "single-knockout stage was disabled in the workflow config"
        ),
    )

    design_columns = [
        "rank",
        "method",
        "product",
        "knockouts",
        "n_knockouts",
        "growth",
        "max_product",
        "guaranteed_product",
        "growth_coupled",
    ]
    for method, design in (
        ("optknock", result.optknock_result),
        ("robustknock", result.robustknock_result),
    ):
        frame = (
            design.to_frame()
            if design is not None
            else pd.DataFrame(columns=design_columns)
        )
        search_complete = _strain_design_search_complete(design)
        design_status: ArtifactStatus = (
            "skipped"
            if not result.config.run_strain_design
            else "complete"
            if search_complete
            else "partial"
        )
        search_status = _strain_design_search_status(design)
        writer.csv(
            f"05_strain_design/{method}.csv",
            frame,
            stage="05_strain_design",
            role=method,
            method=method,
            status=design_status,
            reason=(
                "strain-design stage was disabled in the workflow config"
                if design_status == "skipped"
                else (
                    f"straindesign search ended with status {search_status!r}; returned "
                    "rows are feasible incumbents but their rank/search is incomplete"
                )
                if design_status == "partial"
                else None
            ),
        )

    fseof_tidy_reason: str | None
    if result.fseof_result is None:
        fseof_frame = pd.DataFrame(
            columns=[
                "reaction_id",
                "amplification_rank",
                "proposal_method",
                "method_rank",
                "report_selected",
                "slope",
                "classification",
                "actionable",
                "loop_diagnostic_status",
                "loop_artifact_flag",
                "loopless_to_standard_capacity_ratio",
                "loop_diagnostic_reason",
            ]
        )
        fseof_trajectory = pd.DataFrame(
            columns=[
                "target",
                "enforced_product_flux",
                "reaction_flux",
                "slope",
                "classification",
                "actionable",
                "proposal_method",
                "method_rank",
                "report_selected",
                "loop_diagnostic_status",
                "loop_artifact_flag",
                "loopless_to_standard_capacity_ratio",
                "loop_diagnostic_reason",
            ]
        )
        fseof_tidy_status: ArtifactStatus = (
            "skipped" if not result.config.run_amplification else "failed"
        )
        fseof_tidy_reason = (
            "amplification stage was disabled in the workflow config"
            if not result.config.run_amplification
            else "FSEOF did not return a result"
        )
    else:
        fseof_frame = result.fseof_result.to_frame()
        raw_fseof_targets = tuple(
            result.fseof_result.amplification_targets()[
                : result.config.top_amplification_targets_per_method
            ]
        )
        fseof_targets, fseof_tidy_status, fseof_tidy_reason = (
            _publication_amplification_targets(result, raw_fseof_targets)
        )
        fseof_frame["proposal_method"] = "fseof"
        fseof_frame["method_rank"] = fseof_frame["amplification_rank"]
        fseof_frame["report_selected"] = fseof_frame["amplification_rank"].notna() & (
            fseof_frame["amplification_rank"]
            <= result.config.top_amplification_targets_per_method
        )
        fseof_frame = _annotate_loop_diagnostic(
            fseof_frame,
            result,
            target_column="reaction_id",
        )
        fseof_trajectory = result.fseof_result.trajectory_frame()
        fseof_trajectory = fseof_trajectory[
            fseof_trajectory["reaction_id"].isin(fseof_targets)
        ].rename(columns={"reaction_id": "target", "flux": "reaction_flux"})
        fseof_ranks = {
            target: rank for rank, target in enumerate(raw_fseof_targets, start=1)
        }
        fseof_trajectory["proposal_method"] = "fseof"
        fseof_trajectory["method_rank"] = fseof_trajectory["target"].map(fseof_ranks)
        fseof_trajectory["report_selected"] = True
        fseof_trajectory = _annotate_loop_diagnostic(fseof_trajectory, result)
        fseof_trajectory = fseof_trajectory.sort_values(
            ["method_rank", "enforced_product_flux"], kind="stable"
        )
    writer.csv(
        "06_amplification/fseof.csv",
        fseof_frame,
        stage="06_amplification",
        role="amplification_target_ranking",
        method="fseof",
        status="complete" if result.config.run_amplification else "skipped",
        reason=(
            None
            if result.config.run_amplification
            else "amplification stage was disabled in the workflow config"
        ),
    )
    writer.csv(
        "06_amplification/fseof_tidy.csv",
        fseof_trajectory,
        stage="06_amplification",
        role="fseof_tidy",
        method="fseof",
        status=fseof_tidy_status,
        reason=fseof_tidy_reason,
    )
    fvseof_tidy_reason: str | None
    if result.fvseof_result is None:
        fvseof_frame = pd.DataFrame(
            columns=[
                "reaction_id",
                "amplification_rank",
                "proposal_method",
                "method_rank",
                "report_selected",
                "classification",
                "park_type",
                "robust",
                "slope",
                "capacity_slope",
                "mean_capacity",
                "actionable",
                "loop_diagnostic_status",
                "loop_artifact_flag",
                "loopless_to_standard_capacity_ratio",
                "loop_diagnostic_reason",
            ]
        )
        fvseof_trajectory = pd.DataFrame(
            columns=[
                "target",
                "enforced_product_flux",
                "mean_flux",
                "forced_min_flux",
                "capacity",
                "classification",
                "park_type",
                "robust",
                "actionable",
                "proposal_method",
                "method_rank",
                "report_selected",
                "loop_diagnostic_status",
                "loop_artifact_flag",
                "loopless_to_standard_capacity_ratio",
                "loop_diagnostic_reason",
            ]
        )
        fvseof_tidy_status: ArtifactStatus = (
            "skipped" if not result.config.run_amplification else "failed"
        )
        fvseof_tidy_reason = (
            "amplification stage was disabled in the workflow config"
            if not result.config.run_amplification
            else "FVSEOF did not return a result"
        )
    else:
        fvseof_frame = result.fvseof_result.to_frame()
        raw_fvseof_targets = tuple(
            result.fvseof_result.amplification_targets()[
                : result.config.top_amplification_targets_per_method
            ]
        )
        fvseof_targets, fvseof_tidy_status, fvseof_tidy_reason = (
            _publication_amplification_targets(result, raw_fvseof_targets)
        )
        fvseof_frame["proposal_method"] = "fvseof"
        fvseof_frame["method_rank"] = fvseof_frame["amplification_rank"]
        fvseof_frame["report_selected"] = fvseof_frame["amplification_rank"].notna() & (
            fvseof_frame["amplification_rank"]
            <= result.config.top_amplification_targets_per_method
        )
        fvseof_frame = _annotate_loop_diagnostic(
            fvseof_frame,
            result,
            target_column="reaction_id",
        )
        fvseof_trajectory = result.fvseof_result.trajectory_frame()
        fvseof_trajectory = fvseof_trajectory[
            fvseof_trajectory["reaction_id"].isin(fvseof_targets)
        ].rename(
            columns={
                "reaction_id": "target",
                "forced_min_abs_flux": "forced_min_flux",
            }
        )
        fvseof_ranks = {
            target: rank for rank, target in enumerate(raw_fvseof_targets, start=1)
        }
        fvseof_trajectory["proposal_method"] = "fvseof"
        fvseof_trajectory["method_rank"] = fvseof_trajectory["target"].map(fvseof_ranks)
        fvseof_trajectory["report_selected"] = True
        fvseof_trajectory = _annotate_loop_diagnostic(fvseof_trajectory, result)
        fvseof_trajectory = fvseof_trajectory.sort_values(
            ["method_rank", "enforced_product_flux"], kind="stable"
        )
    writer.csv(
        "06_amplification/fvseof.csv",
        fvseof_frame,
        stage="06_amplification",
        role="variability_supported_amplification_targets",
        method="fvseof",
        status="complete" if result.config.run_amplification else "skipped",
        reason=(
            None
            if result.config.run_amplification
            else "amplification stage was disabled in the workflow config"
        ),
    )
    writer.csv(
        "06_amplification/fvseof_tidy.csv",
        fvseof_trajectory,
        stage="06_amplification",
        role="fvseof_tidy",
        method="fvseof",
        status=fvseof_tidy_status,
        reason=fvseof_tidy_reason,
    )

    _export_validation(writer, result)

    source_suffix = (
        configured_source.suffix
        if configured_source.suffix.lower() in {".xml", ".sbml"}
        else ".xml"
    )
    model_name = f"{_slug(str(source_model.id) or 'model')}{source_suffix}"
    model_relative = f"model/{model_name}"
    if configured_source_bytes is not None:
        (root / model_relative).write_bytes(configured_source_bytes)
    else:
        write_sbml_model(source_model, str(root / model_relative))
    writer.existing(
        model_relative,
        stage="model",
        role="model",
        media_type="application/sbml+xml",
    )
    conditioned_relative = (
        f"model/{_slug(str(conditioned_model.id) or 'model')}__conditioned.xml"
    )
    write_sbml_model(conditioned_model, str(root / conditioned_relative))
    writer.existing(
        conditioned_relative,
        stage="model",
        role="conditioned_model",
        media_type="application/sbml+xml",
    )
    exported_config = _config_payload(result.config)
    exported_config["model_path"] = model_relative
    exported_config["output_dir"] = "."
    exported_provenance = _portable_run_provenance(
        result.provenance,
        model_relative=model_relative,
    )
    writer.json(
        "00_config.json",
        exported_config,
        stage="root",
        role="workflow_configuration",
    )
    writer.json(
        "00_provenance.json",
        exported_provenance,
        stage="root",
        role="provenance",
    )
    writer.json(
        "00_summary.json",
        result.summary(),
        stage="root",
        role="summary",
    )
    _write_analysis_metadata_sidecars(
        writer,
        result,
        exported_provenance=exported_provenance,
    )
    _write_reproduction_scripts(
        writer,
        result,
        root=root,
        model_relative=model_relative,
    )

    manifest_record = ArtifactRecord(
        path="00_manifest.json",
        stage="root",
        role="authoritative_artifact_manifest",
        media_type="application/json",
    )
    all_records = (*writer.records, manifest_record)
    contract_roles = {
        "provenance",
        "summary",
        "model",
        "wild_type_fluxes",
        "theoretical_yield",
        "production_envelope",
        "reference_fluxes",
        "single_knockout_moma",
        "single_knockout_room",
        "optknock",
        "robustknock",
        "amplification_target_ranking",
        "variability_supported_amplification_targets",
        "fseof_tidy",
        "fvseof_tidy",
        "single_knockout_consensus",
        "flux_response_tidy",
        "sampling_tidy",
        "recommendations",
        "gene_knockout_mapping",
        "amplification_loop_diagnostic",
        "reproduction_config",
        "reproduce_script",
        "render_script",
        "validate_script",
    }
    by_role = {record.role: record for record in writer.records}
    artifact_mapping: dict[str, object] = {}
    for role in sorted(contract_roles):
        artifact_record = by_role.get(role)
        if artifact_record is None:
            raise ProductionWorkflowError(
                f"schema-v2 export did not create required artifact role {role!r}"
            )
        entry: dict[str, object] = {
            "path": artifact_record.path,
            "status": artifact_record.status,
            "sha256": artifact_record.sha256,
            "size_bytes": artifact_record.size_bytes,
        }
        if artifact_record.reason is not None:
            entry["reason"] = artifact_record.reason
        metadata_path = artifact_record.metadata_path or writer.metadata_links.get(role)
        if metadata_path is not None:
            entry["metadata_path"] = metadata_path
        artifact_mapping[role] = entry
    primary_roles = {*contract_roles, "authoritative_artifact_manifest"}
    stage_roles = {
        "single_knockout_moma",
        "single_knockout_room",
        "optknock",
        "robustknock",
        "fseof_tidy",
        "fvseof_tidy",
        "flux_response_tidy",
        "sampling_tidy",
        "amplification_loop_diagnostic",
        "recommendations",
    }
    manifest_status = (
        "complete"
        if all(by_role[role].status == "complete" for role in stage_roles)
        else "partial"
    )
    manifest = {
        "schema_id": SCHEMA_ID,
        "schema_version": SCHEMA_VERSION,
        "authoritative": True,
        "status": manifest_status,
        "report": {
            "title": f"CMM production-target discovery: {result.config.product}",
            "product_label": result.config.product,
            "language": "en",
        },
        "directories": list(_STAGE_DIRECTORIES),
        "artifacts": artifact_mapping,
        "supplementary_artifacts": [
            asdict(record) for record in all_records if record.role not in primary_roles
        ],
    }
    (root / "00_manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, sort_keys=True, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    return root, tuple(all_records)


def _write_analysis_metadata_sidecars(
    writer: _ArtifactWriter,
    result: ProductionWorkflowResult,
    *,
    exported_provenance: Mapping[str, object],
) -> None:
    analysis_metadata = exported_provenance.get("analysis_metadata", {})
    metadata = analysis_metadata if isinstance(analysis_metadata, Mapping) else {}
    single_metadata_value = metadata.get("single_knockout", {})
    single_metadata = (
        single_metadata_value if isinstance(single_metadata_value, Mapping) else {}
    )
    room_metadata_value = single_metadata.get("room", {})
    room_metadata = (
        room_metadata_value if isinstance(room_metadata_value, Mapping) else {}
    )
    flux_response_file_metadata = [
        {
            "target": item.target.target_id,
            "background": item.target.background,
            "candidate_scope": item.target.candidate_scope,
            "blocked_reaction_signature": item.target.blocked_reaction_signature,
            "candidate_target_ids": item.target.candidate_target_ids,
            "loop_diagnostic_status": item.target.loop_diagnostic_status,
            "loop_artifact_flag": item.target.loop_artifact_flag,
            "loop_diagnostic_eligible": item.target.loop_diagnostic_eligible,
            "loop_diagnostic_reason": item.target.loop_diagnostic_reason,
            "scan_reference_flux": item.target.scan_reference_flux,
            "status": item.status,
            "error": item.error,
            "reason": item.reason,
            "data_file": (
                f"flux_response__{_slug(item.target.target_id)}.csv"
                if item.result is not None
                else None
            ),
            "phases_file": (
                f"flux_response_phases__{_slug(item.target.target_id)}.csv"
                if item.result is not None
                else None
            ),
            "metadata_file": (
                f"flux_response__{_slug(item.target.target_id)}.metadata.json"
                if item.result is not None
                else None
            ),
            "analysis_metadata": item.result.metadata
            if item.result is not None
            else None,
        }
        for item in result.flux_responses
    ]
    sampling_file_metadata = [
        {
            "target_id": item.target_id,
            "status": item.status,
            "error": item.error,
            "reason": item.reason,
            "candidate_scope": item.candidate_scope,
            "blocked_reaction_signature": item.blocked_reaction_signature,
            "candidate_target_ids": item.candidate_target_ids,
            "samples_file": (
                f"random_sampling__{_slug(item.target_id)}.csv.gz"
                if item.result is not None
                and result.config.validation.sampling.store_raw_samples
                else None
            ),
            "statistics_file": (
                f"random_sampling_statistics__{_slug(item.target_id)}.csv"
                if item.result is not None
                else None
            ),
            "comparison_file": (
                f"random_sampling_comparison__{_slug(item.target_id)}.csv"
                if item.comparison is not None
                else None
            ),
            "metadata_file": (
                f"random_sampling__{_slug(item.target_id)}.metadata.json"
                if item.result is not None
                else None
            ),
            "analysis_metadata": item.result.metadata
            if item.result is not None
            else None,
        }
        for item in result.sampling
    ]
    policy = {
        "single_gene_knockout": (
            "requires a selected beneficial MOMA-L2 or ROOM prediction, a supporting "
            "pre-deletion target-reaction-to-product titration toward zero under the "
            "declared growth floor, concordant positive mean and median paired-sampling "
            "product shifts after complete deletion, and retained growth"
        ),
        "amplification": (
            "requires an independently ranked FSEOF or FVSEOF candidate, a supporting "
            "wild-type target-to-product response, and a completed loopless FVA diagnostic "
            "with no artifact flag; cross-method agreement is recorded but not required"
        ),
        "multi_knockout": (
            "requires a reaction-level RobustKnock design with positive guaranteed product "
            "and retained growth; experimental interpretation requires GPR resolution and "
            "is not exported as a gene recipe"
        ),
        "combined_knockout_amplification": (
            "never recommended because this workflow does not validate combined interventions"
        ),
    }
    specs: tuple[tuple[str, str, str, object], ...] = (
        (
            "preflight_checks",
            "01_preflight/preflight.metadata.json",
            "01_preflight",
            exported_provenance,
        ),
        (
            "theoretical_yield",
            "02_yield/theoretical_yield.metadata.json",
            "02_yield",
            result.theoretical_yield.metadata,
        ),
        (
            "production_envelope",
            "02_yield/production_envelope.metadata.json",
            "02_yield",
            result.envelope.metadata,
        ),
        (
            "wild_type_fluxes",
            "03_reference/wild_type_fluxes.metadata.json",
            "03_reference",
            result.wild_type.metadata,
        ),
        (
            "reference_fluxes",
            "03_reference/reference_fluxes.metadata.json",
            "03_reference",
            dict(result.reference.metadata),
        ),
        (
            "gene_knockout_mapping",
            "04_single_knockout/gene_knockout_mapping.metadata.json",
            "04_single_knockout",
            result.gene_knockout_mapping.metadata,
        ),
        (
            "single_knockout_moma",
            "04_single_knockout/single_knockout_moma.metadata.json",
            "04_single_knockout",
            single_metadata.get("moma_l2"),
        ),
        (
            "single_knockout_room",
            "04_single_knockout/single_knockout_room.metadata.json",
            "04_single_knockout",
            single_metadata.get("room"),
        ),
        (
            "ranked_unique_single_gene_candidates",
            "04_single_knockout/single_knockout_candidates.metadata.json",
            "04_single_knockout",
            {
                "ranking": (
                    "descending product flux, descending growth, blocked-reaction signature, "
                    "target id; product improvement and viability are required"
                ),
                "display_ranking": (
                    "display_rank is independent of recommendation: optimal and "
                    "viability-qualified rows, descending product flux then growth, one "
                    "representative per blocked-reaction signature"
                ),
                "validation_candidate_policy": (
                    "all MOMA-L2 and ROOM rows with display_rank are validation candidates; "
                    "equivalent gene ids sharing a blocked-reaction signature are executed "
                    "once using the declared validation_target_id"
                ),
                "candidate_scope": "all_display_ranked_candidates",
                "top_per_method": result.config.top_single_knockouts_per_method,
                "viability_fraction": result.config.viability_fraction,
                "moma_product_delta_threshold": (
                    result.config.product_improvement_tolerance
                ),
                "room_epsilon": room_metadata.get("selection_room_epsilon"),
                "room_numeric_margin": room_metadata.get("selection_numeric_margin"),
                "room_product_delta_threshold": room_metadata.get(
                    "selection_product_delta_threshold"
                ),
                "product_delta_criterion": "strictly_greater_than_threshold",
            },
        ),
        (
            "single_knockout_consensus",
            "04_single_knockout/single_knockout_consensus.metadata.json",
            "04_single_knockout",
            {"source_roles": ["single_knockout_moma", "single_knockout_room"]},
        ),
        (
            "optknock",
            "05_strain_design/optknock.metadata.json",
            "05_strain_design",
            metadata.get("optknock"),
        ),
        (
            "robustknock",
            "05_strain_design/robustknock.metadata.json",
            "05_strain_design",
            metadata.get("robustknock"),
        ),
        (
            "amplification_target_ranking",
            "06_amplification/fseof.metadata.json",
            "06_amplification",
            {
                "analysis_metadata": metadata.get("fseof"),
                "publication_filter": (
                    "fseof.csv carries the independent FSEOF rank and report_selected "
                    "flag; top report targets retain loop flags for disclosure"
                ),
            },
        ),
        (
            "fseof_tidy",
            "06_amplification/fseof_tidy.metadata.json",
            "06_amplification",
            {
                "analysis_metadata": metadata.get("fseof"),
                "publication_filter": (
                    "independent top-ranked FSEOF trajectories, including flagged targets "
                    "as diagnostic-only; flags never imply recommendation eligibility"
                ),
            },
        ),
        (
            "variability_supported_amplification_targets",
            "06_amplification/fvseof.metadata.json",
            "06_amplification",
            {
                "analysis_metadata": metadata.get("fvseof"),
                "publication_filter": (
                    "fvseof.csv carries the independent FVSEOF rank and report_selected "
                    "flag; top report targets retain loop flags for disclosure"
                ),
            },
        ),
        (
            "fvseof_tidy",
            "06_amplification/fvseof_tidy.metadata.json",
            "06_amplification",
            {
                "analysis_metadata": metadata.get("fvseof"),
                "publication_filter": (
                    "independent top-ranked FVSEOF trajectories, including flagged targets "
                    "as diagnostic-only; flags never imply recommendation eligibility"
                ),
            },
        ),
        (
            "amplification_loop_diagnostic",
            "07_validation/amplification_loop_diagnostic.metadata.json",
            "07_validation",
            (
                result.amplification_loop_diagnostic.metadata
                if result.amplification_loop_diagnostic is not None
                else None
            ),
        ),
        (
            "flux_response_validation_index",
            "07_validation/flux_response_index.metadata.json",
            "07_validation",
            {"files": flux_response_file_metadata},
        ),
        (
            "flux_response_tidy",
            "07_validation/flux_response_tidy.metadata.json",
            "07_validation",
            {"files": flux_response_file_metadata},
        ),
        (
            "single_knockout_sampling_validation_index",
            "07_validation/random_sampling_index.metadata.json",
            "07_validation",
            {"files": sampling_file_metadata},
        ),
        (
            "sampling_tidy",
            "07_validation/sampling_tidy.metadata.json",
            "07_validation",
            {"files": sampling_file_metadata},
        ),
        (
            "recommendations",
            "07_validation/recommendations.metadata.json",
            "07_validation",
            {"evidence_policy": policy},
        ),
    )
    records_by_role = {record.role: record for record in writer.records}
    for role, relative, stage, payload in specs:
        artifact = records_by_role.get(role)
        if artifact is None:
            continue
        writer.metadata(
            role,
            relative,
            {
                "artifact_role": role,
                "artifact_status": artifact.status,
                "artifact_reason": artifact.reason,
                "analysis_metadata": payload,
            },
            stage=stage,
        )


def _write_reproduction_scripts(
    writer: _ArtifactWriter,
    result: ProductionWorkflowResult,
    *,
    root: Path,
    model_relative: str,
) -> None:
    config_payload = _config_payload(result.config)
    config_payload["model_path"] = f"../{model_relative}"
    config_payload["output_dir"] = f"../../{root.name}__reproduced"
    config_payload["overwrite"] = False
    writer.json(
        "scripts/production_config.json",
        config_payload,
        stage="scripts",
        role="reproduction_config",
    )
    writer.text(
        "scripts/reproduce.py",
        """#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

from cmm.workflows.production import ProductionWorkflowConfig, run_production_target_discovery

HERE = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce this CMM production run")
    parser.add_argument("--output-dir", type=Path, help="Override the sibling output directory")
    args = parser.parse_args()
    config = ProductionWorkflowConfig.from_json(HERE / "production_config.json")
    if args.output_dir is not None:
        config = replace(config, output_dir=args.output_dir.resolve(), overwrite=False)
    result = run_production_target_discovery(config)
    print(result.run_directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
""",
        stage="scripts",
        role="reproduce_script",
        executable=True,
    )
    writer.text(
        "scripts/render.py",
        """#!/usr/bin/env python3
from pathlib import Path

from cmm.reporting import render_production_report

RUN_DIR = Path(__file__).resolve().parents[1]
bundle = render_production_report(RUN_DIR, renderer="nature-r")
print(bundle.report.report_standalone_html)
""",
        stage="scripts",
        role="render_script",
        executable=True,
    )
    writer.text(
        "scripts/validate.py",
        """#!/usr/bin/env python3
import json
from pathlib import Path

from cmm.reporting import validate_production_run

RUN_DIR = Path(__file__).resolve().parents[1]
report = validate_production_run(RUN_DIR)
print(json.dumps({"valid": report.valid, "issues": report.issues, "warnings": report.warnings}, indent=2))
raise SystemExit(0 if report.valid else 1)
""",
        stage="scripts",
        role="validate_script",
        executable=True,
    )


def _export_validation(
    writer: _ArtifactWriter,
    result: ProductionWorkflowResult,
) -> None:
    diagnostic = result.amplification_loop_diagnostic
    diagnostic_reason: str | None
    if diagnostic is None:
        diagnostic_frame = AmplificationLoopDiagnosticResult((), {}).to_frame()
        diagnostic_status: ArtifactStatus = "skipped"
        diagnostic_reason = (
            "amplification stage was disabled in the workflow config"
            if not result.config.run_amplification
            else "loopless amplification diagnostic was disabled in the workflow config"
        )
    else:
        diagnostic_frame = diagnostic.to_frame()
        diagnostic_failures = sum(
            record.diagnostic_status == "failed" for record in diagnostic.records
        )
        diagnostic_unresolved = sum(
            record.diagnostic_status != "complete" for record in diagnostic.records
        )
        diagnostic_status = (
            "failed"
            if diagnostic.records and diagnostic_failures == len(diagnostic.records)
            else "partial"
            if diagnostic_unresolved
            else "complete"
        )
        diagnostic_reason = (
            f"{diagnostic_unresolved} candidate diagnostic(s) were not conclusive, "
            f"including {diagnostic_failures} failure(s)"
            if diagnostic_unresolved
            else None
        )
    writer.csv(
        "07_validation/amplification_loop_diagnostic.csv",
        diagnostic_frame,
        stage="07_validation",
        role="amplification_loop_diagnostic",
        method="standard_fva+loopless_fva",
        status=diagnostic_status,
        reason=diagnostic_reason,
    )

    response_index: list[dict[str, object]] = []
    response_tidy: list[pd.DataFrame] = []
    for response_item in result.flux_responses:
        slug = _slug(response_item.target.target_id)
        data_file = (
            f"flux_response__{slug}.csv" if response_item.result is not None else None
        )
        phases_file = (
            f"flux_response_phases__{slug}.csv"
            if response_item.result is not None
            else None
        )
        metadata_file = (
            f"flux_response__{slug}.metadata.json"
            if response_item.result is not None
            else None
        )
        response_index.append(
            {
                "target": response_item.target.target_id,
                "scan_reaction": response_item.target.scan_reaction,
                "response_reaction": response_item.target.response_reaction,
                "background": response_item.target.background,
                "blocked_reactions": ";".join(response_item.target.blocked_reactions),
                "actions": ";".join(response_item.target.actions),
                "source_methods": ";".join(response_item.target.source_methods),
                "candidate_scope": response_item.target.candidate_scope,
                "blocked_reaction_signature": (
                    response_item.target.blocked_reaction_signature
                ),
                "candidate_target_ids": ";".join(
                    response_item.target.candidate_target_ids
                ),
                "loop_diagnostic_status": (response_item.target.loop_diagnostic_status),
                "loop_artifact_flag": response_item.target.loop_artifact_flag,
                "loop_diagnostic_eligible": (
                    response_item.target.loop_diagnostic_eligible
                ),
                "loop_diagnostic_reason": (response_item.target.loop_diagnostic_reason),
                "scan_reference_flux": response_item.target.scan_reference_flux,
                "status": response_item.status,
                "error": response_item.error,
                "reason": response_item.reason,
                "data_file": data_file,
                "phases_file": phases_file,
                "metadata_file": metadata_file,
            }
        )
        if response_item.result is None:
            continue
        assert data_file is not None
        assert phases_file is not None
        assert metadata_file is not None
        tidy = response_item.result.to_frame()
        tidy.insert(0, "target", response_item.target.target_id)
        tidy["scan_reaction"] = response_item.target.scan_reaction
        tidy["response_reaction"] = response_item.target.response_reaction
        tidy["background"] = response_item.target.background
        tidy["candidate_scope"] = response_item.target.candidate_scope
        response_tidy.append(tidy)
        writer.csv(
            f"07_validation/flux_response__{slug}.csv",
            response_item.result.to_frame(),
            stage="07_validation",
            role="candidate_flux_response",
            method="flux_response",
            metadata_path=f"07_validation/{metadata_file}",
        )
        writer.csv(
            f"07_validation/flux_response_phases__{slug}.csv",
            response_item.result.phases_frame(),
            stage="07_validation",
            role="candidate_flux_response_phases",
            method="flux_response",
            metadata_path=f"07_validation/{metadata_file}",
        )
        writer.json(
            f"07_validation/{metadata_file}",
            {
                "target": asdict(response_item.target),
                "status": response_item.status,
                "error": response_item.error,
                "reason": response_item.reason,
                "data_file": data_file,
                "phases_file": phases_file,
                "analysis_metadata": response_item.result.metadata,
            },
            stage="07_validation",
            role="candidate_flux_response_metadata",
        )
    response_failures = sum(item.status == "failed" for item in result.flux_responses)
    response_skips = sum(item.status == "skipped" for item in result.flux_responses)
    response_artifact_status: ArtifactStatus = (
        "skipped"
        if not result.config.validation.enabled
        else "failed"
        if result.flux_responses and response_failures == len(result.flux_responses)
        else "partial"
        if response_failures or response_skips
        else "complete"
    )
    response_reason = (
        "validation stage was disabled in the workflow config"
        if response_artifact_status == "skipped"
        else (
            f"{response_failures} flux-response validation(s) failed; "
            f"{response_skips} were unavailable"
        )
        if response_failures or response_skips
        else None
    )
    writer.csv(
        "07_validation/flux_response_index.csv",
        pd.DataFrame(
            response_index,
            columns=[
                "target",
                "scan_reaction",
                "response_reaction",
                "background",
                "blocked_reactions",
                "actions",
                "source_methods",
                "candidate_scope",
                "blocked_reaction_signature",
                "candidate_target_ids",
                "loop_diagnostic_status",
                "loop_artifact_flag",
                "loop_diagnostic_eligible",
                "loop_diagnostic_reason",
                "scan_reference_flux",
                "status",
                "error",
                "reason",
                "data_file",
                "phases_file",
                "metadata_file",
            ],
        ),
        stage="07_validation",
        role="flux_response_validation_index",
        method="flux_response",
        status=response_artifact_status,
        reason=response_reason,
    )
    response_columns = [
        "target",
        "target_flux",
        "response_flux",
        "biomass_flux",
        "status",
        "scan_reaction",
        "response_reaction",
        "background",
        "candidate_scope",
    ]
    writer.csv(
        "07_validation/flux_response_tidy.csv",
        (
            pd.concat(response_tidy, ignore_index=True)[response_columns]
            if response_tidy
            else pd.DataFrame(columns=response_columns)
        ),
        stage="07_validation",
        role="flux_response_tidy",
        method="flux_response",
        status=response_artifact_status,
        reason=response_reason,
    )

    sampling_index: list[dict[str, object]] = []
    sampling_tidy: list[pd.DataFrame] = []
    store_raw = result.config.validation.sampling.store_raw_samples
    wild_type_item = next(
        (
            item
            for item in result.sampling
            if item.target_id == "wild_type" and item.result is not None
        ),
        None,
    )
    for sampling_item in result.sampling:
        slug = _slug(sampling_item.target_id)
        samples_file = (
            f"random_sampling__{slug}.csv.gz"
            if sampling_item.result is not None and store_raw
            else None
        )
        statistics_file = (
            f"random_sampling_statistics__{slug}.csv"
            if sampling_item.result is not None
            else None
        )
        comparison_file = (
            f"random_sampling_comparison__{slug}.csv"
            if sampling_item.comparison is not None
            else None
        )
        metadata_file = (
            f"random_sampling__{slug}.metadata.json"
            if sampling_item.result is not None
            else None
        )
        sampling_index.append(
            {
                "target_id": sampling_item.target_id,
                "blocked_reactions": ";".join(sampling_item.blocked_reactions),
                "source_methods": ";".join(sampling_item.source_methods),
                "candidate_scope": sampling_item.candidate_scope,
                "blocked_reaction_signature": (
                    sampling_item.blocked_reaction_signature
                ),
                "candidate_target_ids": ";".join(sampling_item.candidate_target_ids),
                "status": sampling_item.status,
                "error": sampling_item.error,
                "reason": sampling_item.reason,
                "samples_file": samples_file,
                "statistics_file": statistics_file,
                "comparison_file": comparison_file,
                "metadata_file": metadata_file,
            }
        )
        if sampling_item.result is None:
            continue
        assert statistics_file is not None
        assert metadata_file is not None
        if store_raw:
            writer.csv(
                f"07_validation/random_sampling__{slug}.csv.gz",
                sampling_item.result.to_frame(),
                stage="07_validation",
                role="raw_flux_samples",
                method="random_flux_sampling",
                metadata_path=f"07_validation/{metadata_file}",
            )
        writer.csv(
            f"07_validation/random_sampling_statistics__{slug}.csv",
            sampling_item.result.statistics().reset_index(),
            stage="07_validation",
            role="flux_sampling_statistics",
            method="random_flux_sampling",
            metadata_path=f"07_validation/{metadata_file}",
        )
        if sampling_item.comparison is not None:
            writer.csv(
                f"07_validation/random_sampling_comparison__{slug}.csv",
                sampling_item.comparison,
                stage="07_validation",
                role="wild_type_knockout_distribution_shift",
                method="random_flux_sampling",
                metadata_path=f"07_validation/{metadata_file}",
            )
        writer.json(
            f"07_validation/{metadata_file}",
            {
                "target_id": sampling_item.target_id,
                "blocked_reactions": sampling_item.blocked_reactions,
                "source_methods": sampling_item.source_methods,
                "status": sampling_item.status,
                "candidate_scope": sampling_item.candidate_scope,
                "blocked_reaction_signature": (
                    sampling_item.blocked_reaction_signature
                ),
                "candidate_target_ids": sampling_item.candidate_target_ids,
                "samples_file": samples_file,
                "statistics_file": statistics_file,
                "comparison_file": comparison_file,
                "analysis_metadata": sampling_item.result.metadata,
                "comparison_definition": (
                    "paired wild-type versus knockout marginal flux summaries, including "
                    "median delta and one-dimensional Wasserstein distance"
                ),
            },
            stage="07_validation",
            role="candidate_sampling_metadata",
        )
        if (
            sampling_item.target_id != "wild_type"
            and wild_type_item is not None
            and wild_type_item.result is not None
        ):
            reactions = {
                result.config.product,
                result.biomass,
                *sampling_item.blocked_reactions,
            }
            if sampling_item.comparison is not None:
                shifted = sampling_item.comparison.nlargest(3, "wasserstein_distance")[
                    "reaction_id"
                ]
                reactions.update(str(value) for value in shifted)
            for condition, sampled in (
                ("wild_type", wild_type_item.result),
                ("knockout", sampling_item.result),
            ):
                tidy = sampled.long_frame()
                tidy = tidy[tidy["reaction_id"].isin(reactions)]
                tidy.insert(0, "condition", condition)
                tidy.insert(0, "target", sampling_item.target_id)
                sampling_tidy.append(tidy)
    sampling_enabled = (
        result.config.validation.enabled
        and result.config.validation.sampling.enabled
        and result.config.run_single_knockout
    )
    sampling_failures = sum(item.status == "failed" for item in result.sampling)
    sampling_skipped = sum(item.status == "skipped" for item in result.sampling)
    wild_type_sampling_failed = any(
        item.target_id == "wild_type" and item.status == "failed"
        for item in result.sampling
    )
    sampling_status: ArtifactStatus = (
        "skipped"
        if not sampling_enabled
        else "failed"
        if wild_type_sampling_failed
        else "partial"
        if sampling_failures or sampling_skipped
        else "complete"
    )
    sampling_reason = (
        "sampling validation was disabled or its single-knockout stage was skipped"
        if sampling_status == "skipped"
        else (
            f"{sampling_failures} sampling validation(s) failed and "
            f"{sampling_skipped} were skipped"
        )
        if sampling_failures or sampling_skipped
        else None
    )
    writer.csv(
        "07_validation/random_sampling_index.csv",
        pd.DataFrame(
            sampling_index,
            columns=[
                "target_id",
                "blocked_reactions",
                "source_methods",
                "candidate_scope",
                "blocked_reaction_signature",
                "candidate_target_ids",
                "status",
                "error",
                "reason",
                "samples_file",
                "statistics_file",
                "comparison_file",
                "metadata_file",
            ],
        ),
        stage="07_validation",
        role="single_knockout_sampling_validation_index",
        method="random_flux_sampling",
        status=sampling_status,
        reason=sampling_reason,
    )
    sampling_columns = ["target", "condition", "reaction_id", "flux", "sample_id"]
    writer.csv(
        "07_validation/sampling_tidy.csv",
        (
            pd.concat(sampling_tidy, ignore_index=True)[sampling_columns]
            if sampling_tidy
            else pd.DataFrame(columns=sampling_columns)
        ),
        stage="07_validation",
        role="sampling_tidy",
        method="random_flux_sampling",
        status=sampling_status,
        reason=sampling_reason,
    )

    recommendations = _recommendations_frame(result)
    recommendation_status, recommendation_reason = _recommendation_artifact_status(
        result
    )
    writer.csv(
        "07_validation/recommendations.csv",
        recommendations,
        stage="07_validation",
        role="recommendations",
        status=recommendation_status,
        reason=recommendation_reason,
    )


_RECOMMENDATION_COLUMNS = [
    "target",
    "type",
    "evidence",
    "verdict",
    "proposal_methods",
    "validation_methods",
    "growth_retained",
    "product_effect",
    "artifact_flag",
    "reason",
]


def _recommendation_artifact_status(
    result: ProductionWorkflowResult,
) -> tuple[ArtifactStatus, str | None]:
    proposal_stages = (
        result.config.run_single_knockout,
        result.config.run_strain_design,
        result.config.run_amplification,
    )
    if not any(proposal_stages):
        return (
            "skipped",
            "all intervention-proposal stages were disabled in the workflow config",
        )
    limitations: list[str] = []
    if result.config.run_single_knockout and (
        not result.config.validation.enabled
        or not result.config.validation.sampling.enabled
    ):
        limitations.append(
            "single-gene recommendations lack required forward validation"
        )
    if result.config.run_amplification and (
        not result.config.validation.enabled
        or not result.config.run_amplification_loop_diagnostic
    ):
        limitations.append("amplification recommendations lack required validation")
    if result.config.run_strain_design and any(
        not _strain_design_search_complete(design)
        for design in (result.optknock_result, result.robustknock_result)
    ):
        limitations.append(
            "strain-design search was incomplete; feasible incumbents are not promoted "
            "as ranked recommendations"
        )
    if any(item.status == "failed" for item in result.flux_responses):
        limitations.append("one or more flux-response validations failed")
    if any(item.status == "skipped" for item in result.flux_responses):
        limitations.append(
            "one or more candidate flux-response validations were unavailable"
        )
    if any(item.status == "failed" for item in result.sampling):
        limitations.append("one or more sampling validations failed")
    if any(item.status == "skipped" for item in result.sampling):
        limitations.append(
            "one or more candidate sampling validations were unavailable"
        )
    if result.amplification_loop_diagnostic is not None and any(
        record.diagnostic_status != "complete"
        for record in result.amplification_loop_diagnostic.records
    ):
        limitations.append("one or more loopless FVA diagnostics were inconclusive")
    if limitations:
        return "partial", "; ".join(limitations)
    return "complete", None


def _recommendations_frame(result: ProductionWorkflowResult) -> pd.DataFrame:
    """Return recommendations only when the workflow produced the declared evidence.

    This table deliberately does not combine a knockout with an amplification target.  Such
    a pair is a new intervention whose phenotype has not been simulated by this workflow.
    Rejected and unvalidated hypotheses remain available in their full source tables.
    """

    rows: list[dict[str, object]] = []
    wild_type_growth = result.reference.get(result.biomass)
    wild_type_product = result.reference.get(result.config.product)
    tolerance = result.config.product_improvement_tolerance

    selected_by_gene: dict[str, list[SingleKnockoutRecord]] = {}
    for record in result.selected_single_knockouts:
        selected_by_gene.setdefault(record.target_id, []).append(record)
    for target_id, records in sorted(selected_by_gene.items()):
        feasibility_support = _knockout_response_support(
            result,
            target_id,
            wild_type_growth=wild_type_growth,
            tolerance=tolerance,
        )
        sampling_support = _knockout_sampling_support(
            result,
            target_id,
            wild_type_growth=wild_type_growth,
            tolerance=tolerance,
        )
        if feasibility_support is None or sampling_support is None:
            continue
        titration_product_effect, feasibility_growth = feasibility_support
        (
            sampling_effect,
            sampling_mean_effect,
            sampling_median_effect,
            sampling_growth,
        ) = sampling_support
        knockout_proposal_methods = tuple(sorted({record.method for record in records}))
        rows.append(
            {
                "target": target_id,
                "type": "single_gene_knockout",
                "evidence": (
                    "pre-deletion zero-flux titration product delta="
                    f"{titration_product_effect:.6g}; sampling mean product delta="
                    f"{sampling_mean_effect:.6g}; sampling median product delta="
                    f"{sampling_median_effect:.6g}"
                ),
                "verdict": "support",
                "proposal_methods": ";".join(knockout_proposal_methods),
                "validation_methods": "flux_response;random_flux_sampling",
                "growth_retained": min(feasibility_growth, sampling_growth),
                "product_effect": sampling_effect,
                "artifact_flag": False,
                "reason": (
                    "suppressing the single blocked reaction from its reference flux to "
                    "zero increased product under the growth floor, and complete deletion "
                    "had concordant positive mean and median product shifts in sampled "
                    "flux distributions"
                ),
            }
        )

    fseof_targets = set(
        (
            result.fseof_result.amplification_targets()[
                : result.config.top_amplification_targets_per_method
            ]
        )
        if result.fseof_result is not None
        else ()
    )
    fvseof_targets = set(
        (
            result.fvseof_result.amplification_targets()[
                : result.config.top_amplification_targets_per_method
            ]
        )
        if result.fvseof_result is not None
        else ()
    )
    loop_by_target = {
        record.target: record
        for record in (
            result.amplification_loop_diagnostic.records
            if result.amplification_loop_diagnostic is not None
            else ()
        )
    }
    amplification_sources = {
        target: tuple(
            method
            for method, targets in (
                ("fseof", fseof_targets),
                ("fvseof", fvseof_targets),
            )
            if target in targets
        )
        for target in fseof_targets | fvseof_targets
    }
    for target in sorted(amplification_sources):
        diagnostic = loop_by_target.get(target)
        if (
            diagnostic is None
            or diagnostic.diagnostic_status != "complete"
            or diagnostic.loop_artifact_flag is not False
        ):
            continue
        response_support = _amplification_response_support(
            result,
            target,
            wild_type_growth=wild_type_growth,
            tolerance=tolerance,
        )
        if response_support is None:
            continue
        product_effect, growth_retained = response_support
        ratio = diagnostic.loopless_to_standard_capacity_ratio
        amplification_proposal_methods = amplification_sources[target]
        method_evidence = "/".join(
            method.upper() for method in amplification_proposal_methods
        )
        rows.append(
            {
                "target": target,
                "type": "amplification",
                "evidence": (
                    f"{method_evidence} candidate; response product delta="
                    f"{product_effect:.6g}; "
                    f"loopless/standard capacity ratio={ratio:.6g}"
                    if ratio is not None
                    else f"{method_evidence} candidate; response product delta="
                    f"{product_effect:.6g}"
                ),
                "verdict": "support",
                "proposal_methods": ";".join(amplification_proposal_methods),
                "validation_methods": "flux_response;standard_fva;loopless_fva",
                "growth_retained": growth_retained,
                "product_effect": product_effect,
                "artifact_flag": False,
                "reason": (
                    "method-ranked amplification candidate increased product under forced "
                    "target flux and passed the declared loop-artifact threshold; "
                    "cross-method agreement was not required"
                ),
            }
        )

    if result.robustknock_result is not None and _strain_design_search_complete(
        result.robustknock_result
    ):
        for design in result.robustknock_result.designs:
            growth_retained = (
                design.growth / wild_type_growth if wild_type_growth > 0.0 else 0.0
            )
            product_effect = design.guaranteed_product - wild_type_product
            if (
                not design.growth_coupled
                or growth_retained < result.config.viability_fraction
                or product_effect <= tolerance
            ):
                continue
            rows.append(
                {
                    "target": ";".join(design.knockouts),
                    "type": "multi_knockout",
                    "evidence": (
                        "Reaction-level RobustKnock guaranteed product="
                        f"{design.guaranteed_product:.6g}; GPR resolution required"
                    ),
                    "verdict": "coupled",
                    "proposal_methods": "robustknock",
                    "validation_methods": "guaranteed_product_at_maximum_growth",
                    "growth_retained": growth_retained,
                    "product_effect": product_effect,
                    "artifact_flag": False,
                    "reason": (
                        "reaction-level growth-coupled design with positive worst-case "
                        "product; requires GPR-resolved gene implementation before "
                        "experiment; rank by guaranteed product, not optimistic maximum"
                    ),
                }
            )

    type_order = {
        "single_gene_knockout": 0,
        "multi_knockout": 1,
        "amplification": 2,
    }
    return pd.DataFrame(
        sorted(
            rows,
            key=lambda row: (
                type_order.get(str(row["type"]), len(type_order)),
                -cast(float, row["product_effect"]),
                str(row["target"]),
            ),
        ),
        columns=_RECOMMENDATION_COLUMNS,
    )


def _knockout_response_support(
    result: ProductionWorkflowResult,
    target_id: str,
    *,
    wild_type_growth: float,
    tolerance: float,
) -> tuple[float, float] | None:
    """Return the product gain at the zero-flux endpoint and retained growth.

    The curve is evaluated before deletion: product is maximized while the candidate's one
    blocked reaction is moved from its reproducible reference flux to zero.  The comparison
    therefore isolates the direction of suppressing that reaction under one fixed growth
    floor; complete-deletion adaptation is evaluated separately by MOMA/ROOM and sampling.
    """

    match = next(
        (
            item
            for item in result.flux_responses
            if (
                item.target.target_id == target_id
                or target_id in item.target.candidate_target_ids
            )
            and (
                item.target.candidate_scope == "all_display_ranked_candidates"
                or (
                    not item.target.candidate_scope
                    and item.target.background == "gene_knockout"
                )
            )
            and item.status == "complete"
            and item.result is not None
        ),
        None,
    )
    if (
        match is None
        or match.result is None
        or match.target.scan_reference_flux is None
        or abs(match.target.scan_reference_flux) <= 1e-9
        or wild_type_growth <= 0.0
    ):
        return None
    feasible = match.result.feasible_points()
    if not feasible:
        return None
    zero_point = min(feasible, key=lambda point: abs(point.target_flux))
    reference_point = min(
        feasible,
        key=lambda point: abs(
            point.target_flux - cast(float, match.target.scan_reference_flux)
        ),
    )
    product_effect = zero_point.response_flux - reference_point.response_flux
    viable_growth = result.config.viability_fraction * wild_type_growth
    if (
        abs(zero_point.target_flux) > 1e-8
        or product_effect <= tolerance
        or zero_point.biomass_flux < viable_growth
    ):
        return None
    return (
        product_effect,
        zero_point.biomass_flux / wild_type_growth,
    )


def _knockout_sampling_support(
    result: ProductionWorkflowResult,
    target_id: str,
    *,
    wild_type_growth: float,
    tolerance: float,
) -> tuple[float, float, float, float] | None:
    match = next(
        (
            item
            for item in result.sampling
            if item.target_id == target_id
            or target_id in item.candidate_target_ids
            and item.status == "complete"
            and item.comparison is not None
        ),
        None,
    )
    if match is None or match.comparison is None or wild_type_growth <= 0.0:
        return None
    comparison = match.comparison.set_index("reaction_id")
    if (
        result.config.product not in comparison.index
        or result.biomass not in comparison.index
    ):
        return None
    product_mean_delta = float(comparison.at[result.config.product, "mean_delta"])
    product_median_delta = float(comparison.at[result.config.product, "median_delta"])
    biomass_median = float(comparison.at[result.biomass, "knockout_median"])
    growth_retained = biomass_median / wild_type_growth
    if (
        not math.isfinite(product_mean_delta)
        or product_mean_delta <= tolerance
        or not math.isfinite(product_median_delta)
        or product_median_delta <= tolerance
        or not math.isfinite(growth_retained)
        or growth_retained < result.config.validation.sampling_growth_fraction
    ):
        return None
    conservative_product_delta = min(product_mean_delta, product_median_delta)
    return (
        conservative_product_delta,
        product_mean_delta,
        product_median_delta,
        growth_retained,
    )


def _amplification_response_support(
    result: ProductionWorkflowResult,
    target_id: str,
    *,
    wild_type_growth: float,
    tolerance: float,
) -> tuple[float, float] | None:
    match = next(
        (
            item
            for item in result.flux_responses
            if item.target.target_id == target_id
            and (
                item.target.candidate_scope == "all_report_selected_candidates"
                or (
                    not item.target.candidate_scope
                    and item.target.background == "wild_type"
                    and "amplify" in item.target.actions
                )
            )
            and item.target.background == "wild_type"
            and item.status == "complete"
            and item.result is not None
        ),
        None,
    )
    if match is None or match.result is None or wild_type_growth <= 0.0:
        return None
    baseline_target = float(match.result.wild_type.get("target_flux", 0.0))
    baseline_product = float(match.result.wild_type.get("response_flux", 0.0))
    viable_growth = (
        result.config.validation.flux_response_biomass_fraction * wild_type_growth
    )
    supported = [
        point
        for point in match.result.feasible_points()
        if point.target_flux > baseline_target + tolerance
        and point.response_flux > baseline_product + tolerance
        and point.biomass_flux >= viable_growth
    ]
    if not supported:
        return None
    best = max(supported, key=lambda point: (point.response_flux, point.biomass_flux))
    return (
        best.response_flux - baseline_product,
        best.biomass_flux / wild_type_growth,
    )


def _workflow_analysis_status(result: ProductionWorkflowResult) -> str:
    required_enabled = (
        result.config.run_single_knockout
        and result.config.run_strain_design
        and result.config.run_amplification
        and result.config.run_amplification_loop_diagnostic
        and result.config.validation.enabled
        and result.config.validation.sampling.enabled
    )
    validation_complete = not any(
        item.status == "failed" for item in result.flux_responses
    ) and not any(item.status in {"failed", "skipped"} for item in result.sampling)
    diagnostic_complete = result.amplification_loop_diagnostic is not None and all(
        record.diagnostic_status == "complete"
        for record in result.amplification_loop_diagnostic.records
    )
    design_searches_complete = all(
        _strain_design_search_complete(design)
        for design in (result.optknock_result, result.robustknock_result)
    )
    return (
        "complete"
        if (
            required_enabled
            and validation_complete
            and diagnostic_complete
            and design_searches_complete
        )
        else "partial"
    )


def _summary_warnings(result: ProductionWorkflowResult) -> list[str]:
    warnings: list[str] = []
    if result.medium_application.medium == "model_as_loaded":
        warnings.append(
            "The model medium was used as loaded; exchange bounds may encode an undocumented "
            "condition and should be reviewed before publication."
        )
    if result.medium_application.dropped:
        warnings.append(
            "The requested medium omitted unavailable component(s): "
            + ";".join(result.medium_application.dropped)
            + "."
        )
    oxygen_value = result.provenance.get("oxygen_uptake", [])
    oxygen_records = (
        oxygen_value
        if isinstance(oxygen_value, Sequence)
        and not isinstance(oxygen_value, (str, bytes))
        else ()
    )
    for raw_record in oxygen_records:
        if not isinstance(raw_record, Mapping):
            continue
        uptake_limit = raw_record.get("uptake_limit")
        if not isinstance(uptake_limit, (int, float)) or uptake_limit < 100.0:
            continue
        warnings.append(
            "The oxygen exchange "
            f"{raw_record.get('reaction_id', '<unknown>')} permits uptake up to "
            f"{float(uptake_limit):g} mmol gDW^-1 h^-1 "
            f"(bounds {raw_record.get('lower_bound')}, {raw_record.get('upper_bound')}); "
            "this effectively unlimited aerobic condition must be stated and should be "
            "replaced by an experimentally justified aeration bound for publication."
        )
    if result.theoretical_yield.co2_carbon_fraction > 1e-9:
        warnings.append(
            "Net CO2 uptake supplies "
            f"{result.theoretical_yield.co2_carbon_fraction:.1%} of predicted product "
            "carbon; interpret yield against the organic-carbon ceiling with this "
            "contribution explicit."
        )
    if result.theoretical_yield.carbon_imbalance:
        warnings.append(
            "The predicted yield exceeds the carbon ceiling and net CO2 uptake does not "
            "explain the excess; inspect model exchanges before interpreting the result."
        )
    disabled = []
    if not result.config.run_single_knockout:
        disabled.append("single-knockout screening")
    if not result.config.run_strain_design:
        disabled.append("strain design")
    if not result.config.run_amplification:
        disabled.append("amplification discovery")
    elif not result.config.run_amplification_loop_diagnostic:
        disabled.append("amplification loopless diagnostic")
    if not result.config.validation.enabled:
        disabled.append("forward validation")
    elif not result.config.validation.sampling.enabled:
        disabled.append("random-sampling validation")
    if disabled:
        warnings.append(
            "Partial analysis: publication-required stage(s) were disabled: "
            + "; ".join(disabled)
            + "."
        )
    if any(item.status == "failed" for item in result.flux_responses):
        warnings.append(
            "At least one flux-response validation failed; see its index row."
        )
    if any(item.status == "failed" for item in result.sampling):
        warnings.append("At least one sampling validation failed; see its index row.")
    if any(item.status == "skipped" for item in result.sampling):
        warnings.append(
            "At least one candidate sampling validation was unavailable and skipped; "
            "see its index row."
        )
    incomplete_designs = [
        method
        for method, design in (
            ("OptKnock", result.optknock_result),
            ("RobustKnock", result.robustknock_result),
        )
        if not _strain_design_search_complete(design)
    ]
    if incomplete_designs:
        warnings.append(
            "Strain-design search did not prove completion for "
            + "; ".join(incomplete_designs)
            + "; exported rows are feasible incumbents, not a complete ranking."
        )
    return warnings


def _strain_design_search_status(result: StrainDesignResult | None) -> str | None:
    if result is None:
        return None
    parameters_value = result.metadata.get("parameters", {})
    parameters = parameters_value if isinstance(parameters_value, Mapping) else {}
    status = parameters.get("straindesign_search_status")
    return str(status) if status is not None else None


def _strain_design_search_complete(result: StrainDesignResult | None) -> bool:
    """Whether the backend proved search completion.

    Older/in-memory result fixtures have no completion field and retain their historical
    treatment.  Every newly computed public result records the field explicitly.
    """

    if result is None:
        return True
    parameters_value = result.metadata.get("parameters", {})
    parameters = parameters_value if isinstance(parameters_value, Mapping) else {}
    return parameters.get("straindesign_search_complete") is not False


def _flux_frame(fluxes: Mapping[str, float]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"reaction_id": reaction_id, "flux": float(flux)}
            for reaction_id, flux in sorted(fluxes.items())
        ],
        columns=["reaction_id", "flux"],
    )


def _single_knockout_export_row(
    record: SingleKnockoutRecord,
    model: Model,
) -> dict[str, object]:
    """Enrich a MOMA/ROOM row with reporter-ready, reaction-level identity fields."""

    row = record.to_dict()
    gene = (
        model.genes.get_by_id(record.target_id)
        if record.target_id in model.genes
        else None
    )
    reactions = [
        model.reactions.get_by_id(reaction_id)
        for reaction_id in record.blocked_reactions
        if reaction_id in model.reactions
    ]
    row.update(
        {
            "target_name": str(gene.name or "") if gene is not None else "",
            "blocked_reaction_names": ";".join(
                str(reaction.name or "") for reaction in reactions
            ),
            "blocked_reaction_equations": ";".join(
                reaction.build_reaction_string(use_metabolite_names=False)
                for reaction in reactions
            ),
            "blocked_reaction_gprs": ";".join(
                str(reaction.gene_reaction_rule) for reaction in reactions
            ),
        }
    )
    return row


def _single_knockout_candidate_export_rows(
    records: Sequence[SingleKnockoutRecord],
    model: Model,
) -> list[dict[str, object]]:
    """Export every method-specific display candidate with its deduplicated run target."""

    rows: list[dict[str, object]] = []
    for group in _single_knockout_validation_groups(records):
        representative_record = group.representative
        for record in group.display_records:
            row = _single_knockout_export_row(record, model)
            row.update(
                {
                    "candidate_scope": "all_display_ranked_candidates",
                    "validation_target_id": representative_record.target_id,
                    "candidate_target_ids": ";".join(group.candidate_target_ids),
                    "candidate_source_methods": ";".join(group.source_methods),
                    "validation_representative": record is representative_record,
                }
            )
            rows.append(row)
    return rows


def _portable_run_provenance(
    provenance: Mapping[str, object],
    *,
    model_relative: str,
) -> dict[str, object]:
    """Return an export-safe provenance copy with run-relative source-model paths.

    The exact source bytes and fingerprints are already pinned in ``model/`` and the hash
    fields.  Persisting the original workstation path adds no reproducibility and leaks a
    user-specific directory into otherwise relocatable bundles.  Rewrite only fields whose
    semantic key is ``source_model_path``; arbitrary user notes remain untouched.
    """

    jsonable = _jsonable(provenance)
    if not isinstance(jsonable, dict):  # pragma: no cover - Mapping always becomes dict
        raise TypeError("workflow provenance must serialize to a JSON object")

    def rewrite(value: object) -> object:
        if isinstance(value, dict):
            return {
                str(key): (
                    model_relative if str(key) == "source_model_path" else rewrite(item)
                )
                for key, item in value.items()
            }
        if isinstance(value, list):
            return [rewrite(item) for item in value]
        return value

    portable = rewrite(jsonable)
    if not isinstance(portable, dict):  # pragma: no cover - rewrite preserves the root
        raise TypeError("workflow provenance must remain a JSON object")
    return cast(dict[str, object], portable)


def _config_payload(config: ProductionWorkflowConfig) -> dict[str, object]:
    payload = asdict(config)
    payload["model_path"] = str(cast(Path, config.model_path).expanduser().resolve())
    payload["output_dir"] = (
        str(cast(Path, config.output_dir).expanduser().resolve())
        if config.output_dir is not None
        else None
    )
    payload["medium"] = _medium_payload(config.medium)
    return cast(dict[str, object], _jsonable(payload))


def _medium_payload(medium: Medium | str | None) -> object:
    if medium is None:
        return {"mode": "model_as_loaded"}
    if isinstance(medium, str):
        return {"mode": "preset", "name": medium}
    return {
        "mode": "explicit",
        "name": medium.name,
        "uptake": dict(medium.uptake),
        "required": sorted(medium.required or ()),
    }


def _medium_from_payload(value: object) -> Medium | str | None:
    if value is None or isinstance(value, (Medium, str)):
        return value
    if not isinstance(value, Mapping):
        raise ValueError("medium must be null, a preset name, or an object")
    mode = value.get("mode")
    if mode == "model_as_loaded":
        return None
    if mode == "preset":
        name = value.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError("preset medium requires a non-empty name")
        return name
    if mode not in {None, "explicit"}:
        raise ValueError(f"unknown medium mode {mode!r}")
    name = value.get("name")
    uptake = value.get("uptake")
    required = value.get("required")
    if not isinstance(name, str) or not name or not isinstance(uptake, Mapping):
        raise ValueError("explicit medium requires name and uptake mapping")
    required_ids = (
        frozenset(str(item) for item in required)
        if isinstance(required, Sequence) and not isinstance(required, (str, bytes))
        else None
    )
    return Medium(
        name=name,
        uptake={str(key): float(rate) for key, rate in uptake.items()},
        required=required_ids,
    )


def _condition_from_payload(value: object) -> Condition | None:
    if value is None or isinstance(value, Condition):
        return value
    if not isinstance(value, Mapping):
        raise ValueError("condition must be null or an object")
    raw_bounds = value.get("bounds", ())
    if not isinstance(raw_bounds, Sequence) or isinstance(raw_bounds, (str, bytes)):
        raise ValueError("condition bounds must be a list")
    bounds: list[ReactionBound] = []
    for raw_bound in raw_bounds:
        if not isinstance(raw_bound, Mapping):
            raise ValueError("each condition bound must be an object")
        bounds.append(
            ReactionBound(
                reaction_id=str(raw_bound.get("reaction_id", "")),
                lower_bound=(
                    None
                    if raw_bound.get("lower_bound") is None
                    else float(cast(float, raw_bound["lower_bound"]))
                ),
                upper_bound=(
                    None
                    if raw_bound.get("upper_bound") is None
                    else float(cast(float, raw_bound["upper_bound"]))
                ),
            )
        )
    raw_objective = value.get("objective")
    objective = None
    if raw_objective is not None:
        if not isinstance(raw_objective, Mapping):
            raise ValueError("condition objective must be an object")
        coefficients = raw_objective.get("coefficients")
        if not isinstance(coefficients, Mapping):
            raise ValueError("condition objective requires a coefficients mapping")
        objective = ObjectiveSpec(
            coefficients={
                str(key): float(coefficient)
                for key, coefficient in coefficients.items()
            },
            direction=cast(
                Literal["max", "min"], raw_objective.get("direction", "max")
            ),
        )
    return Condition(
        name=str(value.get("name", "default")),
        bounds=tuple(bounds),
        objective=objective,
        notes=str(value.get("notes", "")),
    )


def _slug(value: str) -> str:
    """Return a bounded, filesystem-safe and collision-resistant display slug."""

    digest = hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]
    display = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-") or "item"
    # Always append the digest: case-only ids collide on common macOS/Windows filesystems,
    # while punctuation-normalized ids (A/B, A?B, A B) otherwise share one basename.
    max_length = 96
    suffix = f"__{digest}"
    prefix = display[: max_length - len(suffix)].rstrip("._-") or "item"
    return f"{prefix}{suffix}"
