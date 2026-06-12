"""Feature module boundaries and service entry points for the CMM workbench."""

from cmm.features.comparison import ComparisonResult, moma, reference_flux, room
from cmm.features.production import (
    fseof,
    fvseof,
    production_envelope,
    theoretical_yield,
)
from cmm.features.revert import revert_targets
from cmm.features.strain_design import (
    StrainDesign,
    StrainDesignResult,
    optknock,
    robustknock,
)
from cmm.features.transformation import direction_from_states, transformation_targets

# Features with a shipped, tested service in ``src/cmm``. Keep this list in sync with the
# actual modules — it is the public contract for what the workbench provides.
INCLUDED_FEATURES = (
    "flux_visualization_slider",
    "fba",
    "fva",
    "pfba",
    "moma",
    "room",
    "theoretical_yield",
    "production_envelope",
    "fseof",
    "fvseof",
    "omics_integration",
    "transformation_targets",
    "optknock",
    "robustknock",
    "revert_metabolism",
)

# Roadmap items described in docs/feature-roadmap.md (Phase 3 / Phase 5) that do NOT yet have
# an implemented service. Listed separately so the manifest never overstates the shipped
# surface. Move an entry into INCLUDED_FEATURES only when its service and tests land.
PLANNED_FEATURES = (
    "dynamic_fba",
    "batch_moma_room",
    "random_flux_sampling",
    "flux_response_analysis",
    "enzyme_constrained_modeling",
)

EXCLUDED_FEATURES = (
    "scenario_templates",
    "media_management",
    "scenario_file_formats",
    "elementary_flux_modes",
    "thermodynamics",
    "legacy_minimal_cut_sets",
)

__all__ = [
    "INCLUDED_FEATURES",
    "PLANNED_FEATURES",
    "EXCLUDED_FEATURES",
    "ComparisonResult",
    "StrainDesign",
    "StrainDesignResult",
    "direction_from_states",
    "fseof",
    "fvseof",
    "moma",
    "production_envelope",
    "reference_flux",
    "optknock",
    "robustknock",
    "room",
    "revert_targets",
    "theoretical_yield",
    "transformation_targets",
]
