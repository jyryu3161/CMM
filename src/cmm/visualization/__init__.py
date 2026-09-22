"""Publication-quality figures for CMM results (matplotlib, headless-safe)."""

from cmm.visualization.jev import (
    jev_decision_figure,
    jev_design_space_figure,
    jev_progress_figure,
)
from cmm.visualization.figures import (
    escher_flux_map,
    flux_comparison_figure,
    flux_log_change_figure,
    flux_response_figure,
    fseof_figure,
    fvseof_figure,
    network_flux_map,
    production_envelope_figure,
    sampling_figure,
    save_figure,
    yield_figure,
)

__all__ = [
    "escher_flux_map",
    "flux_comparison_figure",
    "flux_log_change_figure",
    "flux_response_figure",
    "fseof_figure",
    "fvseof_figure",
    "jev_decision_figure",
    "jev_design_space_figure",
    "jev_progress_figure",
    "network_flux_map",
    "production_envelope_figure",
    "sampling_figure",
    "save_figure",
    "yield_figure",
]
