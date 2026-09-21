"""Feature module boundaries and service entry points for the CMM platform."""

from cmm.features._perturbation import (
    Perturbation,
    PerturbationList,
    blocked_reactions_for_genes,
    gene_perturbations,
    grouped_gene_perturbations,
    perturbation_provenance,
    reaction_perturbations,
)
from cmm.features.coupling import (
    PARALLEL_TOLERANCE,
    CoupledSets,
    coupled_reaction_sets,
)
from cmm.features.comparison import (
    ROOM_TOLERANCES,
    BatchComparisonResult,
    BatchComparisonRow,
    ComparisonResult,
    batch_comparison,
    knockout_comparison,
    moma,
    reference_flux,
    room,
)
from cmm.features.production import (
    CarbonUptake,
    FseofResult,
    FvseofResult,
    ProductionEnvelope,
    ProductionYield,
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.response import (
    FluxResponseResult,
    ResponseLimit,
    ResponsePhase,
    ResponsePoint,
    flux_response,
)
from cmm.features.revert import revert_targets, tie_structure
from cmm.features.sampling import (
    SamplingResult,
    random_flux_sampling,
    reference_constrained_sampling,
)
from cmm.features.strain_design import (
    StrainDesign,
    StrainDesignResult,
    optknock,
    robustknock,
)
from cmm.features.transformation import direction_from_states, transformation_targets

# Features with a shipped, tested service in ``src/cmm``. Keep this list in sync with the
# actual modules — it is the public contract for what the platform provides.
INCLUDED_FEATURES = (
    "flux_visualization_slider",
    "fba",
    "fva",
    "pfba",
    "moma",
    "room",
    "batch_moma_room",
    "theoretical_yield",
    "production_envelope",
    "fseof",
    "fvseof",
    "omics_integration",
    "transformation_targets",
    "optknock",
    "robustknock",
    "revert_metabolism",
    "flux_response_analysis",
    "random_flux_sampling",
    "production_target_workflow",
    "publication_reporting",
    "jev_target_design",
)

# Roadmap items described in docs/feature-roadmap.md (Phase 3 / Phase 5) that do NOT yet have
# an implemented service. Listed separately so the manifest never overstates the shipped
# surface. Move an entry into INCLUDED_FEATURES only when its service and tests land.
PLANNED_FEATURES = (
    "dynamic_fba",
    "enzyme_constrained_modeling",
    "fvseof_grouping_constraints",
)

EXCLUDED_FEATURES = (
    # Generic scenario-template/file-format engines remain outside the shipped surface. The
    # concrete SC-01 workflow and publication reporter above do not make these generic.
    "scenario_templates",
    "media_management",
    "scenario_file_formats",
    "elementary_flux_modes",
    "thermodynamics",
    "legacy_minimal_cut_sets",
)

__all__ = [
    "PARALLEL_TOLERANCE",
    "CoupledSets",
    "coupled_reaction_sets",
    "INCLUDED_FEATURES",
    "PLANNED_FEATURES",
    "EXCLUDED_FEATURES",
    "ROOM_TOLERANCES",
    "BatchComparisonResult",
    "BatchComparisonRow",
    "CarbonUptake",
    "ComparisonResult",
    "FluxResponseResult",
    "FseofResult",
    "FvseofResult",
    "Perturbation",
    "PerturbationList",
    "ProductionEnvelope",
    "ProductionYield",
    "ResponseLimit",
    "ResponsePhase",
    "ResponsePoint",
    "SamplingResult",
    "StrainDesign",
    "StrainDesignResult",
    "batch_comparison",
    "blocked_reactions_for_genes",
    "direction_from_states",
    "flux_response",
    "fseof",
    "fvseof",
    "gene_perturbations",
    "grouped_gene_perturbations",
    "knockout_comparison",
    "moma",
    "perturbation_provenance",
    "production_envelope",
    "random_flux_sampling",
    "reaction_perturbations",
    "reference_constrained_sampling",
    "reference_flux",
    "optknock",
    "robustknock",
    "room",
    "revert_targets",
    "theoretical_yield",
    "tie_structure",
    "transformation_targets",
]
