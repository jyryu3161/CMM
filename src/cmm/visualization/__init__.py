"""Publication-quality figures for CMM results (matplotlib, headless-safe)."""

try:
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

    # Inside the same guard: the agent figures import matplotlib directly and build on
    # ``figures``, so a base install fails on them for exactly the same reason.
    from cmm.visualization.jev import (
        jev_decision_figure,
        jev_design_space_figure,
        jev_progress_figure,
    )
except ModuleNotFoundError as error:  # pragma: no cover - depends on the install extras
    # Name the extra rather than leaving a bare "No module named 'matplotlib'": a base
    # install carries the solver-backed services but not the figure layer.
    if error.name != "matplotlib":
        raise
    raise ModuleNotFoundError(
        "cmm.visualization requires matplotlib, which ships in the 'desktop' extra: "
        "install it with `pip install 'cmm[desktop]'`. The numerical services in "
        "cmm.core, cmm.features, cmm.omics and cmm.workflows do not need it."
    ) from error

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
