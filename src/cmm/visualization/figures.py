"""Publication-quality figures for production-design results.

Figures are built with `matplotlib.figure.Figure` directly (no pyplot global state), so the
module is import-safe and headless. Styling targets print/paper use: configurable single- or
double-column font presets, clean spines, a colour-blind-safe palette, and readable type.

`network_flux_map` is a schematic flux network (force-directed, edge width ∝ |flux|). It is
intentionally NOT a curated Escher map: hand-laid Escher metabolic maps encode biochemical
layout that this lightweight, dependency-free renderer does not attempt. Use it for a quick
flux overview; use a curated Escher map for a publication network figure.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from matplotlib import __version__ as matplotlib_version
from matplotlib import cm, colormaps
from matplotlib import image as plt_image
from matplotlib.colors import ListedColormap, Normalize, TwoSlopeNorm
from matplotlib.figure import Figure
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

from cmm.features.production import (
    FseofResult,
    FvseofResult,
    ProductionEnvelope,
    ProductionYield,
)
from cmm.features.response import FluxResponseResult
from cmm.features.sampling import SamplingResult

# matplotlib 3.10 replaced violinplot's boolean `vert` with `orientation`; the project
# supports 3.8+, so pick the spelling the installed version accepts.
_HORIZONTAL: dict[str, object] = (
    {"orientation": "horizontal"}
    if tuple(int(part) for part in matplotlib_version.split(".")[:2]) >= (3, 10)
    else {"vert": False}
)

# Okabe-Ito colour-blind-safe palette.
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999"]
_LINESTYLES = ["-", "--", "-.", ":"]
_FONT_PRESETS = {
    1: {"title": 10, "label": 9, "tick": 8, "legend": 8},  # single column (~3.3 in)
    2: {"title": 13, "label": 12, "tick": 10, "legend": 9},  # double column (~6.5 in)
}


def _new_figure(width: float = 6.0, height: float = 4.0, column_width: int = 2):
    fig = Figure(figsize=(width, height), dpi=300)
    fig.set_facecolor("white")
    ax = fig.subplots()
    return fig, ax, _FONT_PRESETS.get(column_width, _FONT_PRESETS[2])


def _style(ax, font, *, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel, fontsize=font["label"])
    ax.set_ylabel(ylabel, fontsize=font["label"])
    ax.set_title(title, fontsize=font["title"], fontweight="bold", pad=10)
    ax.tick_params(labelsize=font["tick"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)


def save_figure(fig: Figure, path: str | Path, dpi: int = 300) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    return path


def production_envelope_figure(
    envelope: ProductionEnvelope,
    *,
    product_label: str | None = None,
    title: str = "Production envelope",
    column_width: int = 2,
) -> Figure:
    """Growth-vs-product phenotypic phase plane with the feasible region shaded."""

    fig, ax, font = _new_figure(column_width=column_width)
    frame = envelope.to_frame().sort_values("product_flux")
    x = frame["product_flux"].to_numpy()
    gmax = frame["growth_max"].to_numpy()
    gmin = frame["growth_min"].to_numpy()
    ax.fill_between(
        x, gmin, gmax, color=PALETTE[0], alpha=0.18, label="feasible region"
    )
    ax.plot(x, gmax, color=PALETTE[0], linewidth=2.0, label="max growth")
    ax.plot(
        x,
        gmin,
        color=PALETTE[0],
        linewidth=1.0,
        linestyle="--",
        alpha=0.7,
        label="min growth",
    )
    label = product_label or envelope.product
    suffix = "  (min growth = 0 across range)" if float(np.max(gmin)) <= 1e-9 else ""
    _style(
        ax,
        font,
        xlabel=f"{label} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        ylabel="growth rate (h$^{-1}$)",
        title=title + suffix,
    )
    ax.legend(fontsize=font["legend"], frameon=False)
    fig.tight_layout()
    return fig


def fseof_figure(
    result: FseofResult,
    *,
    top_n: int = 6,
    title: str = "FSEOF amplification targets",
    column_width: int = 2,
) -> Figure:
    """Flux of the top amplification-target reactions vs enforced product flux."""

    fig, ax, font = _new_figure(width=7.5, height=4.2, column_width=column_width)
    levels = list(result.enforced_levels)
    targets = result.amplification_targets()
    ranked = sorted(
        targets,
        key=lambda rid: (
            abs(result.trends.loc[rid, levels[-1]])
            - abs(result.trends.loc[rid, levels[0]])
        ),
        reverse=True,
    )[:top_n]
    for i, rid in enumerate(ranked):
        ys = [result.trends.loc[rid, lv] for lv in levels]
        # Vary colour AND linestyle/marker so near-coincident lines stay distinguishable.
        ax.plot(
            levels,
            ys,
            marker=["o", "s", "^", "D", "v", "P"][i % 6],
            markersize=6,
            linewidth=2.4 if i < 2 else 1.6,
            linestyle=_LINESTYLES[i % len(_LINESTYLES)],
            color=PALETTE[i % len(PALETTE)],
            label=rid,
        )
    _style(
        ax,
        font,
        xlabel=f"enforced {result.product} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        ylabel="reaction flux (mmol gDW$^{-1}$ h$^{-1}$)",
        title=title,
    )
    ax.legend(
        fontsize=font["legend"],
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
    )
    fig.tight_layout()
    return fig


def fvseof_figure(
    result: FvseofResult,
    *,
    top_n: int = 5,
    title: str = "FVSEOF robust amplification targets",
    column_width: int = 2,
) -> Figure:
    """Mean flux (solid) and forced-minimum |flux| (dashed) of robust targets vs product.

    A target whose forced-minimum rises with enforced product is *robustly* amplified — the
    reaction cannot avoid carrying more flux. That dashed line rising is the FVSEOF signal.
    """

    fig, ax, font = _new_figure(width=7.5, height=4.2, column_width=column_width)
    levels = list(result.enforced_levels)
    level_cols = list(result.mean.columns)
    robust = result.robust_targets() or result.amplification_targets()
    ranked = sorted(
        robust,
        key=lambda rid: (
            abs(result.mean.loc[rid, level_cols[-1]])
            - abs(result.mean.loc[rid, level_cols[0]])
        ),
        reverse=True,
    )[:top_n]
    for i, rid in enumerate(ranked):
        colour = PALETTE[i % len(PALETTE)]
        mean = [result.mean.loc[rid, c] for c in level_cols]
        forced = [result.forced.loc[rid, c] for c in level_cols]
        ax.plot(
            levels,
            mean,
            marker="o",
            markersize=5,
            linewidth=2.0,
            color=colour,
            label=f"{rid} (mean)",
        )
        ax.plot(levels, forced, linestyle="--", linewidth=1.4, color=colour, alpha=0.8)
    _style(
        ax,
        font,
        xlabel=f"enforced {result.product} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        ylabel="flux (mmol gDW$^{-1}$ h$^{-1}$)",
        title=title,
    )
    ax.legend(
        fontsize=font["legend"],
        frameon=False,
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
    )
    ax.text(
        0.0,
        -0.16,
        "solid = mean flux, dashed = forced minimum |flux|",
        transform=ax.transAxes,
        fontsize=font["legend"],
        color="#555555",
    )
    fig.tight_layout()
    return fig


def _infeasible_spans(frame) -> list[tuple[float, float]]:
    """Contiguous runs of infeasible scan points, as x-ranges bounded by their neighbours."""

    xs = frame["target_flux"].to_numpy(dtype=float)
    bad = (frame["status"] != "optimal").to_numpy()
    spans: list[tuple[float, float]] = []
    index = 0
    while index < len(bad):
        if not bad[index]:
            index += 1
            continue
        start = index
        while index + 1 < len(bad) and bad[index + 1]:
            index += 1
        left = xs[start] if start == 0 else (xs[start - 1] + xs[start]) / 2.0
        right = (
            xs[index] if index == len(bad) - 1 else (xs[index] + xs[index + 1]) / 2.0
        )
        spans.append((float(left), float(right)))
        index += 1
    return spans


def flux_response_figure(
    result: FluxResponseResult,
    *,
    title: str | None = None,
    show_biomass: bool = True,
    column_width: int = 2,
) -> Figure:
    """Response curve vs enforced target flux, with the response limit and phases marked.

    The curve is the optimal-value function of an LP parameterised in one right-hand side,
    so it is concave piecewise linear: the light dividers are the phase boundaries where the
    shadow price changes, and the solid marker is ``result.limit`` — the exact target flux
    beyond which the response starts to fall. Neither moves with ``n_steps``, which the
    finite-difference "steepest decline" span they replace did, by up to 29.5 flux units.

    When the response is a product (not biomass itself) the biomass trace is drawn on a
    secondary axis, because the growth cost of each target flux is what decides whether a
    point on the curve is reachable by a living strain.
    """

    # Aspect close to the square-ish panel the GUI gives this plot, so the curve keeps the
    # width it needs instead of being squeezed by a legend reserved outside the axes.
    fig, ax, font = _new_figure(width=6.4, height=4.6, column_width=column_width)
    frame = result.to_frame()
    feasible = frame[frame["status"] == "optimal"]

    # Shade contiguous stretches with no solution at all — the edge of viability. Spans run
    # to the midpoint of the neighbouring feasible point so the boundary is not overstated.
    for start, end in _infeasible_spans(frame):
        ax.axvspan(start, end, color="#BBBBBB", alpha=0.35, linewidth=0)

    ax.plot(
        feasible["target_flux"],
        feasible["response_flux"],
        marker="o",
        markersize=5,
        linewidth=2.2,
        color=PALETTE[0],
        label=f"{result.response} (maximized)",
    )

    wild_type_flux = result.wild_type["target_flux"]
    ax.axvline(
        wild_type_flux,
        color=PALETTE[6],
        linewidth=1.4,
        linestyle=":",
        label=f"wild type ({wild_type_flux:.3g})",
    )

    optimum = result.optimum()
    if optimum is not None:
        ax.plot(
            [optimum.target_flux],
            [optimum.response_flux],
            marker="*",
            markersize=15,
            color=PALETTE[1],
            linestyle="none",
            label=f"optimum ({optimum.response_flux:.3g})",
        )

    # Phase boundaries: where the LP dual (and therefore the slope) changes. Drawn as light
    # dividers, and deliberately unlabelled when there are many, so they read as structure
    # rather than as a claim about any single one of them.
    for phase in result.phases[1:]:
        ax.axvline(
            phase.target_low,
            color="#BBBBBB",
            linewidth=0.8,
            linestyle="-",
            zorder=0,
        )

    limit = result.limit
    if (
        limit.found
        and limit.target_flux is not None
        and limit.response_flux is not None
    ):
        ax.axvline(
            limit.target_flux,
            color=PALETTE[1],
            linewidth=1.6,
            linestyle="--",
            label=(
                f"response limit ({limit.target_flux:.3g}; shadow price "
                f"{limit.shadow_price_after:.3g})"
            ),
        )

    handles, labels = ax.get_legend_handles_labels()
    if show_biomass and result.response != result.biomass and not feasible.empty:
        twin = ax.twinx()
        twin.plot(
            feasible["target_flux"],
            feasible["biomass_flux"],
            linewidth=1.6,
            linestyle="--",
            color=PALETTE[2],
            label=f"{result.biomass}",
        )
        twin.set_ylabel("growth rate (h$^{-1}$)", fontsize=font["label"])
        twin.tick_params(labelsize=font["tick"])
        twin.spines["top"].set_visible(False)
        # Anchor at zero: with a growth floor the biomass trace is constant, and an
        # auto-scaled axis would magnify solver noise into an apparent trend.
        twin.set_ylim(0.0, max(float(feasible["biomass_flux"].max()), 1e-9) * 1.15)
        twin_handles, twin_labels = twin.get_legend_handles_labels()
        handles += twin_handles
        labels += twin_labels
        # A twin is drawn over its host, so the host's legend would sit under the biomass
        # trace. Lift the host above it (hiding its own opaque background so the twin's line
        # still shows) and keep the legend on the host, where "best" can see the response
        # curve it needs to avoid.
        ax.set_zorder(twin.get_zorder() + 1)
        ax.patch.set_visible(False)

    _style(
        ax,
        font,
        xlabel=f"enforced {result.target} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        ylabel=f"{result.response} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        title=title or f"{result.response} response to {result.target}",
    )
    # Inside the axes, not reserved beside them: an outside legend costs this plot roughly a
    # third of its width and pushes the title past the figure edge on a narrow canvas.
    legend = ax.legend(
        handles,
        labels,
        fontsize=font["legend"],
        loc="best",
        framealpha=1.0,
        facecolor="white",
        edgecolor="#CCCCCC",
    )
    legend.get_frame().set_linewidth(0.6)
    fig.set_layout_engine("tight")
    return fig


def sampling_figure(
    result: SamplingResult,
    *,
    reactions: list[str] | None = None,
    top_n: int = 8,
    reference: dict[str, float] | None = None,
    reference_label: str = "reference",
    title: str = "Sampled flux distributions",
    column_width: int = 2,
) -> Figure:
    """Violin plot of sampled flux per reaction, optionally against a reference solution.

    Reactions default to the most variable ones, since a reaction with a near-constant
    sampled flux is exactly the one whose value a single FBA solution already settled. A
    reference marker far out in a violin's tail is the visual form of "this predicted flux
    is one of many alternate optima".
    """

    samples = result.samples
    if reactions is None:
        spread = samples.std().sort_values(ascending=False)
        selected = [str(rid) for rid in spread.index[:top_n]]
    else:
        missing = [rid for rid in reactions if rid not in samples.columns]
        if missing:
            raise KeyError(f"reactions not present in the sample table: {missing}")
        selected = list(reactions)

    if not selected:
        raise ValueError("no reactions to plot")

    fig, ax, font = _new_figure(
        width=7.0, height=max(3.2, 0.5 * len(selected) + 1.6), column_width=column_width
    )
    data = [samples[rid].to_numpy() for rid in selected]
    positions = np.arange(1, len(selected) + 1)
    parts = ax.violinplot(
        data, positions=positions, showmedians=True, widths=0.8, **_HORIZONTAL
    )
    for body in parts["bodies"]:
        body.set_facecolor(PALETTE[0])
        body.set_alpha(0.45)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        if key in parts:
            parts[key].set_color(PALETTE[0])
            parts[key].set_linewidth(1.2)

    if reference is not None:
        values = [reference.get(rid, float("nan")) for rid in selected]
        ax.plot(
            values,
            positions,
            marker="D",
            markersize=7,
            linestyle="none",
            color=PALETTE[1],
            label=reference_label,
            zorder=5,
        )
        ax.legend(fontsize=font["legend"], frameon=False, loc="best")

    ax.set_yticks(positions)
    ax.set_yticklabels(selected, fontsize=font["tick"])
    ax.axvline(0, color="black", linewidth=0.8, alpha=0.6)
    _style(
        ax,
        font,
        xlabel="flux (mmol gDW$^{-1}$ h$^{-1}$)",
        ylabel="reaction",
        title="",
    )
    # Run parameters sit on their own line under the heading. Stacking rather than placing
    # them beside the title keeps them legible on a narrow panel, where a side-by-side
    # secondary title collides, and a caption pinned below the axes hits the x label.
    ax.set_title(
        f"n = {result.n_samples}, {result.method}, seed {result.seed}",
        fontsize=font["legend"],
        color="#555555",
        pad=6,
    )
    fig.suptitle(title, fontsize=font["title"], fontweight="bold")
    # A layout engine re-runs on every draw, so the labels stay clear when the GUI resizes
    # the canvas; a one-shot tight_layout() only fits the size the figure was authored at.
    fig.set_layout_engine("tight")
    return fig


def flux_comparison_figure(
    reference: dict[str, float],
    comparison: dict[str, float],
    reactions: list[str],
    *,
    reference_label: str = "wild type",
    comparison_label: str = "engineered",
    title: str = "Flux comparison",
    column_width: int = 2,
) -> Figure:
    """Grouped bar chart comparing two flux distributions over selected reactions."""

    fig, ax, font = _new_figure(
        width=max(6.0, 0.7 * len(reactions) + 2), height=4.0, column_width=column_width
    )
    idx = np.arange(len(reactions))
    width = 0.38
    ref = [reference.get(r, 0.0) for r in reactions]
    cmp = [comparison.get(r, 0.0) for r in reactions]
    ax.bar(idx - width / 2, ref, width, color=PALETTE[6], label=reference_label)
    ax.bar(idx + width / 2, cmp, width, color=PALETTE[1], label=comparison_label)
    ax.set_xticks(idx)
    ax.set_xticklabels(reactions, rotation=45, ha="right", fontsize=font["tick"])
    ax.axhline(0, color="black", linewidth=0.8)
    _style(
        ax,
        font,
        xlabel="reaction",
        ylabel="flux (mmol gDW$^{-1}$ h$^{-1}$)",
        title=title,
    )
    ax.legend(fontsize=font["legend"], frameon=False)
    fig.tight_layout()
    return fig


def yield_figure(
    yields: list[ProductionYield],
    *,
    substrate_label: str | None = None,
    title: str = "Theoretical molar yield",
    column_width: int = 2,
) -> Figure:
    """Bar chart of theoretical molar yields (e.g. aerobic vs anaerobic).

    A dashed line marks the substrate carbon ceiling; bars above it are annotated with
    ``+CO2`` to flag that the yield is achievable only with net CO2 fixation.
    """

    fig, ax, font = _new_figure(width=5.2, height=4.2, column_width=column_width)
    labels = [f"{'aerobic' if y.aerobic else 'anaerobic'}" for y in yields]
    values = [y.molar_yield for y in yields]
    colors = [PALETTE[0] if y.aerobic else PALETTE[2] for y in yields]
    bars = ax.bar(labels, values, color=colors, width=0.55)
    for bar, y in zip(bars, yields, strict=True):
        text = f"{y.molar_yield:.2f}" + (" +CO$_2$" if y.exceeds_carbon_ceiling else "")
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            y.molar_yield,
            text,
            ha="center",
            va="bottom",
            fontsize=font["tick"],
        )
    ceiling = next(
        (y.carbon_ceiling for y in yields if y.carbon_ceiling is not None), None
    )
    if ceiling is not None:
        ax.axhline(ceiling, color="#999999", linestyle="--", linewidth=1.0)
        ax.text(
            ax.get_xlim()[1],
            ceiling,
            f" carbon ceiling {ceiling:.2f}",
            va="bottom",
            ha="right",
            fontsize=font["legend"],
            color="#555555",
        )
    substrate = substrate_label or (yields[0].substrate if yields else "substrate")
    product = yields[0].product if yields else "product"
    _style(
        ax,
        font,
        xlabel="condition",
        ylabel=f"mol {product} / mol {substrate}",
        title=title,
    )
    fig.tight_layout()
    return fig


# Currency / cofactor metabolites that connect to nearly everything and turn a metabolic
# network drawing into a hairball. Excluded so the map shows the carbon backbone.
_CURRENCY = frozenset(
    {
        "h2o",
        "h",
        "atp",
        "adp",
        "amp",
        "pi",
        "ppi",
        "co2",
        "o2",
        "nh4",
        "nad",
        "nadh",
        "nadp",
        "nadph",
        "fad",
        "fadh2",
        "coa",
        "q8",
        "q8h2",
        "so4",
        "h2",
        "pyr_h",
    }
)


def _base_id(mid: str) -> str:
    return mid.rsplit("_", 1)[0]


def flux_log_change_figure(
    log_changes: dict[str, float],
    *,
    top_n: int = 20,
    title: str = "Flux log2 fold-change",
    source_label: str = "A",
    target_label: str = "B",
    column_width: int = 2,
) -> Figure:
    """Horizontal diverging bar chart of the largest flux log2 fold-changes between conditions."""

    ranked = sorted(log_changes.items(), key=lambda kv: -abs(kv[1]))[:top_n]
    ranked = ranked[::-1]  # largest at the top after barh
    labels = [rid for rid, _ in ranked]
    values = [v for _, v in ranked]

    fig, ax, font = _new_figure(
        width=6.5, height=max(3.5, 0.32 * len(ranked) + 1), column_width=column_width
    )
    colors = [PALETTE[1] if v > 0 else PALETTE[0] for v in values]
    ax.barh(range(len(values)), values, color=colors)
    ax.set_yticks(range(len(values)))
    ax.set_yticklabels(labels, fontsize=font["tick"])
    ax.axvline(0, color="black", linewidth=0.8)
    _style(
        ax,
        font,
        xlabel=f"log$_2$( |flux {target_label}| / |flux {source_label}| )",
        ylabel="reaction",
        title=title,
    )
    ax.grid(True, axis="y", alpha=0.0)
    fig.tight_layout()
    return fig


def network_flux_map(
    model,
    fluxes: dict[str, float],
    *,
    top_n: int = 12,
    title: str = "Flux network (top reactions)",
    seed: int = 0,
    line_width: float | None = None,
) -> Figure:
    """Schematic carbon-backbone flux network: top-|flux| reactions as edges, width ∝ |flux|.

    Currency metabolites (ATP, H2O, CO2, NAD(P)H, ...) are excluded so the layout follows the
    carbon skeleton. A deterministic force-directed layout; this is a quick flux overview, not
    a curated Escher map (see the module docstring).
    """

    fig, ax, font = _new_figure(width=8.0, height=4.8)
    ax.set_axis_off()
    # The network is drawn edge-to-edge with no axes frame, so a title set tight against it
    # reads as a label stuck on the drawing rather than a caption above it.
    ax.set_title(title, fontsize=font["title"], fontweight="bold", pad=14)

    met_ids: list[str] = []
    edges: list[tuple[int, int, float, str]] = []
    index: dict[str, int] = {}

    def _node(mid: str) -> int:
        if mid not in index:
            index[mid] = len(met_ids)
            met_ids.append(mid)
        return index[mid]

    def _carbon_mets(mets):
        return [m for m in mets if _base_id(m.id) not in _CURRENCY]

    ranked = sorted(
        ((rid, flux) for rid, flux in fluxes.items() if abs(flux) > 1e-6),
        key=lambda kv: -abs(kv[1]),
    )
    for rid, flux in ranked:
        rxn = model.reactions.get_by_id(rid)
        if rxn.boundary:
            continue
        reactants = _carbon_mets([m for m, c in rxn.metabolites.items() if c < 0])
        products = _carbon_mets([m for m, c in rxn.metabolites.items() if c > 0])
        if not reactants or not products:
            continue  # pure cofactor reaction — nothing on the carbon backbone to draw
        if flux < 0:
            reactants, products = products, reactants
        edges.append((_node(reactants[0].id), _node(products[0].id), abs(flux), rid))
        if len(edges) >= top_n:
            break

    n = len(met_ids)
    if n == 0:
        ax.text(
            0.5, 0.5, "no internal carbon flux to display", ha="center", va="center"
        )
        return fig

    _ = seed  # retained for API compatibility; the row layout is deterministic
    pos = _component_row_layout(n, [(s, d) for s, d, _, _ in edges])
    max_w = max((w for _, _, w, _ in edges), default=1.0)
    cmap = colormaps["viridis"]
    reaction_font = max(6, font["legend"] - 2)
    for s, d, w, rid in edges:
        x0, y0 = pos[s]
        x1, y1 = pos[d]
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>",
                color=cmap(0.15 + 0.85 * w / max_w),
                lw=line_width if line_width is not None else 1.2 + 5.0 * w / max_w,
                alpha=0.9,
                shrinkA=14,
                shrinkB=14,
            ),
        )
        ax.text(
            (x0 + x1) / 2,
            (y0 + y1) / 2 - 0.055,
            rid,
            fontsize=reaction_font,
            ha="center",
            va="center",
            color="#222222",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#cccccc", alpha=0.85),
        )
    ax.scatter(
        pos[:, 0], pos[:, 1], s=130, color=PALETTE[0], zorder=3, edgecolors="white"
    )
    metabolite_font = max(7, font["tick"] - 2)
    for i, (mid, (x, y)) in enumerate(zip(met_ids, pos, strict=True)):
        above = i % 2 == 0
        y_offset = 0.055 if above else -0.09
        ax.text(
            x,
            y + y_offset,
            _display_metabolite_id(mid),
            fontsize=metabolite_font,
            fontweight="bold",
            ha="center",
            va="bottom" if above else "top",
            color="#11243b",
        )
    # A formula in the margin is not a legend: the reader has to hold "∝ |flux|" in their head
    # and guess which end of viridis is large. Draw the scale the arrows were actually coloured
    # from — the same truncated span passed to `cmap` above, so the bar cannot lie about them.
    bar_cmap = ListedColormap(cmap(np.linspace(0.15, 1.0, 256)))
    mappable = cm.ScalarMappable(norm=Normalize(vmin=0.0, vmax=max_w), cmap=bar_cmap)
    cbar = fig.colorbar(
        mappable, ax=ax, fraction=0.035, pad=0.02, shrink=0.55, aspect=14
    )
    cbar.set_label("|flux| (mmol gDW$^{-1}$ h$^{-1}$)", fontsize=font["legend"])
    cbar.ax.tick_params(labelsize=font["tick"])
    ax.text(
        0.5,
        -0.10,
        (
            "arrow colour scales with |flux|; direction is the net direction"
            if line_width is not None
            else "arrow colour and width both scale with |flux|; "
            "direction is the net direction"
        ),
        transform=ax.transAxes,
        fontsize=font["legend"],
        color="#555555",
        ha="center",
    )
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.14, 1.02)
    # `colorbar` re-anchors its parent right; centre the network in the space that is left.
    ax.set_anchor("C")
    fig.set_layout_engine("constrained")
    return fig


#: Longest chain drawn on one row before it wraps. A single connected component can be far
#: longer than the panel is wide — at 25 reactions the core carbon backbone is one ~30-node
#: chain — and spreading it across one row collapses the spacing until nodes and labels
#: overlap into a smear. Wrapping trades a fold for legible spacing.
_MAX_NODES_PER_ROW = 10


def _component_row_layout(
    n: int, edges: list[tuple[int, int]], per_row: int = _MAX_NODES_PER_ROW
) -> np.ndarray:
    """Pack the carbon-backbone components into compact rows, folding long chains.

    Rows are laid at one fixed node spacing and centred, so an arrow is the same length
    everywhere and a two-node component no longer stretches across the whole panel. Long
    chains fold boustrophedon — each fold reverses direction — so the node ending one row is
    directly above the node starting the next instead of a full panel width away.
    """

    adjacency = {i: set() for i in range(n)}
    for source, dest in edges:
        adjacency[source].add(dest)
        adjacency[dest].add(source)

    seen: set[int] = set()
    components: list[list[int]] = []
    for node in range(n):
        if node in seen:
            continue
        stack = [node]
        seen.add(node)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))

    components.sort(key=lambda nodes: (-len(nodes), nodes[0]))

    rows: list[list[int]] = []
    for component in components:
        ordered = _ordered_component_nodes(component, edges)
        for start in range(0, len(ordered), per_row):
            chunk = ordered[start : start + per_row]
            if (start // per_row) % 2:
                chunk = chunk[::-1]  # fold back, keeping the join between rows short
            rows.append(chunk)

    y_values = np.linspace(0.88, 0.14, len(rows)) if len(rows) > 1 else [0.55]
    step = 0.84 / max(per_row - 1, 1)
    pos = np.zeros((n, 2))
    for chunk, y in zip(rows, y_values, strict=True):
        x0 = 0.5 - step * (len(chunk) - 1) / 2  # centre each row on the panel
        for offset, node in enumerate(chunk):
            pos[node] = (x0 + offset * step, float(y))
    return pos


def _ordered_component_nodes(
    component: list[int], edges: list[tuple[int, int]]
) -> list[int]:
    """Order one component from likely upstream sources to downstream products."""

    nodes = set(component)
    outgoing = {node: [] for node in component}
    indegree = {node: 0 for node in component}
    for source, dest in edges:
        if source in nodes and dest in nodes:
            outgoing[source].append(dest)
            indegree[dest] += 1

    sources = [node for node in component if indegree[node] == 0 and outgoing[node]]
    if not sources:
        sources = [component[0]]

    depths: dict[int, int] = {}
    queue = list(sorted(sources))
    for source in queue:
        depths[source] = 0
    while queue:
        current = queue.pop(0)
        for dest in sorted(outgoing[current]):
            if dest not in depths:
                depths[dest] = depths[current] + 1
                queue.append(dest)

    fallback = max(depths.values(), default=0) + 1
    return sorted(component, key=lambda node: (depths.get(node, fallback), node))


def _display_metabolite_id(mid: str) -> str:
    return _base_id(mid).replace("__", "-")


def _map_canvas_extent(body: dict) -> tuple[float, float, float, float]:
    """Where a drawing of this map belongs, in the map's coordinates.

    An Escher map records the canvas its layout was drawn on, so a picture exported from that
    same map lines up with the node coordinates without the caller nudging it. Escher's y grows
    downward, which is also how the axes are set below, so the extent is given top-edge-last.
    """

    canvas = body.get("canvas")
    if canvas:
        x, y = float(canvas["x"]), float(canvas["y"])
        return (x, x + float(canvas["width"]), y + float(canvas["height"]), y)
    xs = [float(n["x"]) for n in body["nodes"].values()]
    ys = [float(n["y"]) for n in body["nodes"].values()]
    return (min(xs), max(xs), max(ys), min(ys))


def _fit_extent(
    box: tuple[float, float, float, float], width: int, height: int
) -> tuple[float, float, float, float]:
    """Largest sub-box of ``box`` with the image's own proportions, centred in it.

    A picture placed on an extent is stretched to fill it, so a drawing whose proportions do
    not match the map's canvas would be distorted — circles into ellipses, a map skewed out of
    true. Fitting inside and centring leaves margin instead, which is honest about the
    mismatch; an export of this very map matches exactly and is unaffected.
    """

    left, right, bottom, top = box
    box_w, box_h = abs(right - left), abs(bottom - top)
    if not (box_w and box_h and width and height):
        return box
    scale = min(box_w / width, box_h / height)
    half_w, half_h = width * scale / 2, height * scale / 2
    cx, cy = (left + right) / 2, (top + bottom) / 2
    # Preserve the caller's y direction: Escher's y grows downward, so bottom > top here.
    down = bottom > top
    return (
        cx - half_w,
        cx + half_w,
        cy + half_h if down else cy - half_h,
        cy - half_h if down else cy + half_h,
    )


def escher_flux_map(
    map_path: str | Path,
    fluxes: dict[str, float],
    *,
    title: str | None = None,
    abs_max: float | None = None,
    label_metabolites: bool | None = None,
    label_reactions: bool | None = None,
    width: float = 12.0,
    background: "np.ndarray | str | Path | None" = None,
    background_extent: tuple[float, float, float, float] | None = None,
    background_alpha: float = 1.0,
    line_width: float | None = None,
    font_scale: float = 1.0,
    label_min_fraction: float | None = None,
) -> Figure:
    """Render an Escher map (curated node/segment layout) coloured by flux.

    This reuses a standard Escher map JSON's hand-laid coordinates and bezier segments — the
    same layout a curated Escher map uses — so the figure matches that biochemical layout
    instead of an invented force layout. Reaction edge width and colour encode flux
    (diverging: blue = negative/reverse, red = positive/forward). ``map_path`` is supplied by
    the caller (CMM bundles no maps).
        ``line_width`` overrides the flux-proportional stroke width with one constant width, so
    colour alone encodes magnitude. Useful over a ``background``, where a wide stroke covers
    the drawing underneath it.

    ``font_scale`` multiplies every text size. Type is measured in points and the drawing in
    inches, so a caller that authors the figure at one width and displays it at another
    changes the ratio between them: rendered into a 6.5-inch panel, a figure authored at 9
    inches has text 1.4x larger relative to the network than it was drawn to have. Authoring
    at the display width is the real fix and this is the adjustment on top of it.

    ``label_min_fraction`` labels only the reactions carrying at least that fraction of the
    largest flux. On a map where half the network is at zero, naming all of it spends the
    space that the reactions actually doing something need, and the names collide.
    """

    # A background drawing of this map already carries the labels and the metabolite circles,
    # so drawing CMM's on top of them doubles every one and the text turns to mush. Unless the
    # caller asks for them explicitly, the background is left to say what things are called and
    # this function contributes only the flux.
    if label_metabolites is None:
        label_metabolites = background is None
    if label_reactions is None:
        label_reactions = background is None
    draw_nodes = background is None

    with open(map_path) as handle:
        data = json.load(handle)
    body = data[1] if isinstance(data, list) else data
    nodes = body["nodes"]
    reactions = body["reactions"]
    node_xy = {nid: (float(n["x"]), float(n["y"])) for nid, n in nodes.items()}

    fig = Figure(figsize=(width, width * 0.82), dpi=200)
    fig.set_facecolor("white")
    ax = fig.subplots()
    ax.set_axis_off()
    ax.set_aspect("equal")

    if background is not None:
        # A drawing of the same map, laid underneath in the map's own coordinates. An Escher
        # export renders the network the way Escher draws it — arrowheads, node sizes, the
        # lot — which this function does not attempt; the flux strokes then go on top, so the
        # picture is Escher's and the colouring is CMM's. Rasterisation belongs to the caller
        # (SVG needs a renderer, and this module stays free of Qt), so an array is accepted
        # directly and a path only for formats matplotlib itself reads.
        image = background
        if isinstance(image, (str, Path)):
            if str(image).lower().endswith(".svg"):
                raise ValueError(
                    "escher_flux_map cannot rasterise SVG. Render it to an RGBA array first "
                    "and pass that (cmm.app.svg_background does this with Qt)."
                )
            image = plt_image.imread(str(image))
        extent = background_extent or _fit_extent(
            _map_canvas_extent(body), image.shape[1], image.shape[0]
        )
        ax.imshow(
            image,
            extent=extent,
            alpha=background_alpha,
            interpolation="antialiased",
            zorder=0,
            # Not aspect="auto": imshow writes its aspect onto the axes, and "auto" would
            # discard the equal aspect set above, stretching the whole map to the panel.
            aspect="equal",
        )

    signed = [fluxes.get(r["bigg_id"], 0.0) for r in reactions.values()]
    amax = abs_max or max((abs(v) for v in signed), default=1.0) or 1.0
    norm = TwoSlopeNorm(vmin=-amax, vcenter=0.0, vmax=amax)
    cmap = colormaps["coolwarm"]

    for r in reactions.values():
        flux = fluxes.get(r["bigg_id"], 0.0)
        mag = abs(flux)
        if mag <= 1e-9:
            color = "#d7dbe0"
            lw = line_width if line_width is not None else 0.7
            alpha = 0.7
        else:
            color = cmap(norm(flux))
            # Width and colour carry the same number. Doubling it up makes a large flux
            # unmissable, which is the point on a bare layout; over a drawing of the network it
            # is a thick band hiding the very lines it is drawn on, so the caller can ask for
            # one constant width and let colour alone carry the magnitude.
            lw = line_width if line_width is not None else 0.8 + 5.5 * mag / amax
            alpha = 0.95
        for seg in r["segments"].values():
            a = node_xy.get(seg["from_node_id"])
            b = node_xy.get(seg["to_node_id"])
            if a is None or b is None:
                continue
            b1, b2 = seg.get("b1"), seg.get("b2")
            if b1 and b2:
                verts = [a, (b1["x"], b1["y"]), (b2["x"], b2["y"]), b]
                codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
            else:
                verts = [a, b]
                codes = [MplPath.MOVETO, MplPath.LINETO]
            ax.add_patch(
                PathPatch(
                    MplPath(verts, codes),
                    fill=False,
                    edgecolor=color,
                    lw=lw,
                    alpha=alpha,
                    capstyle="round",
                    joinstyle="round",
                )
            )

    xs, ys = [], []
    for n in nodes.values():
        if n.get("node_type") == "metabolite":
            if draw_nodes:
                ax.plot(n["x"], n["y"], "o", ms=2.2, color="#2b3a4a", zorder=3)
            xs.append(n["x"])
            ys.append(n["y"])
            if label_metabolites:
                ax.text(
                    n.get("label_x", n["x"]),
                    n.get("label_y", n["y"]),
                    str(n.get("bigg_id", "")).rsplit("_", 1)[0],
                    fontsize=4.5 * font_scale,
                    ha="left",
                    va="center",
                    color="#11243b",
                    zorder=4,
                )
    if label_reactions:
        floor = (label_min_fraction or 0.0) * amax
        for r in reactions.values():
            if floor and abs(fluxes.get(r["bigg_id"], 0.0)) < floor:
                # A reaction at rest is drawn, in grey, and not named. Its name would cost the
                # same space as a name that means something, and on a map half of which is at
                # zero that is most of the space.
                continue
            ax.text(
                r.get("label_x", 0),
                r.get("label_y", 0),
                r["bigg_id"],
                fontsize=4.5 * font_scale,
                ha="left",
                va="center",
                color="#7a3b00",
                zorder=4,
                fontstyle="italic",
            )

    if xs:
        mx = 0.04 * (max(xs) - min(xs) + 1)
        my = 0.04 * (max(ys) - min(ys) + 1)
        ax.set_xlim(min(xs) - mx, max(xs) + mx)
        ax.set_ylim(max(ys) + my, min(ys) - my)  # Escher y grows downward
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    # shrink/aspect keep this a legible bar. Left to its defaults against a tall map axes it
    # is drawn as a ~1:23 hairline, which reads as a rule rather than a scale.
    cbar = fig.colorbar(
        mappable, ax=ax, fraction=0.035, pad=0.02, shrink=0.55, aspect=14
    )
    cbar.set_label("flux (mmol gDW$^{-1}$ h$^{-1}$)", fontsize=9 * font_scale)
    cbar.ax.tick_params(labelsize=8 * font_scale)
    # `colorbar` re-anchors its parent to (1.0, 0.5) so the axes hugs the bar. With
    # `set_aspect("equal")` the axes box shrinks whenever the canvas is wider than the map,
    # and a right anchor takes every bit of that shrinkage off the left edge — the map and its
    # title drift right inside the panel. The map should sit in the middle of the space it has.
    ax.set_anchor("C")
    if title:
        # Smaller and further off the drawing than the default: an Escher map fills its axes
        # completely, so a large title sitting on the top edge crowds the network it names.
        ax.set_title(title, fontsize=12 * font_scale, fontweight="bold", pad=16)
    # `tight_layout` solves the margins once, at the size the figure was authored. The GUI
    # canvas then stretches the figure to whatever the panel is, the solved margins no longer
    # fit, and the title is clipped against the top edge. A constrained layout re-solves on
    # every draw, so the figure stays correct at any panel size and at 300 DPI on export.
    fig.set_layout_engine("constrained")
    return fig


def transformation_ranking_figure(
    ranking,
    *,
    highlight: str | None = None,
    top_n: int = 10,
    column_width: int = 2,
) -> Figure:
    """Transformation score against rank, with the leading candidates named.

    A ranking's shape is the first thing to read: a score that falls away sharply over the top
    few candidates says something different from one that decays smoothly across hundreds. The
    band marks the top *N* the report discusses so a reader can see whether the cut sits on a
    cliff or in the middle of a plateau.
    """

    import pandas as pd

    frame = pd.DataFrame(ranking)
    if frame.empty:
        raise ValueError("ranking is empty")
    frame = frame.sort_values("score", ascending=False).reset_index(drop=True)
    frame["rank"] = frame.index + 1

    fig, ax, font = _new_figure(width=7.0, height=4.2, column_width=column_width)
    ax.plot(frame["rank"], frame["score"], lw=1.0, color="#4a5568", zorder=2)
    cut = min(top_n, len(frame))
    ax.axvspan(0.5, cut + 0.5, color="#c6f6d5", alpha=0.55, zorder=0)
    ax.scatter(
        frame["rank"].head(cut),
        frame["score"].head(cut),
        s=26,
        color="#2f855a",
        zorder=3,
        label=f"top {cut}",
    )
    if highlight is not None and (frame["target_id"] == highlight).any():
        row = frame[frame["target_id"] == highlight].iloc[0]
        ax.scatter(
            [row["rank"]],
            [row["score"]],
            s=80,
            color="#e53e3e",
            zorder=5,
            label=f"{highlight} (rank {int(row['rank'])})",
        )
    ax.legend(frameon=False, fontsize=font["tick"])
    if len(frame) > 50:
        ax.set_xscale("log")
    _style(
        ax,
        font,
        xlabel=f"rank among {len(frame)} candidates",
        ylabel="transformation score",
        title="Transformation score by rank",
    )
    fig.tight_layout()
    return fig


def transformation_vs_moma_figure(
    ranking,
    moma_baseline,
    *,
    highlight: str | None = None,
    column_width: int = 2,
) -> Figure:
    """Each candidate's transformation rank against its rank under the MOMA baseline.

    Yizhak et al. report MOMA as markedly inferior for this task, so a transformation ranking
    that merely reproduces MOMA's ordering has not shown its signal comes from the method. Points
    far off the diagonal are where the two disagree; a cloud along it is the warning sign.
    """

    import numpy as np
    import pandas as pd

    left = pd.DataFrame(ranking)
    right = pd.DataFrame(moma_baseline)
    if left.empty or right.empty:
        raise ValueError("both the ranking and the MOMA baseline must be non-empty")
    merged = left[["target_id", "rank"]].merge(
        right[["target_id", "rank"]],
        on="target_id",
        suffixes=("_transformation", "_moma"),
    )
    if merged.empty:
        raise ValueError("the ranking and the MOMA baseline share no candidates")

    fig, ax, font = _new_figure(width=5.6, height=5.0, column_width=column_width)
    ax.scatter(
        merged["rank_transformation"],
        merged["rank_moma"],
        s=10,
        color="#a0aec0",
        alpha=0.6,
        zorder=2,
    )
    limit = max(merged["rank_transformation"].max(), merged["rank_moma"].max())
    ax.plot([1, limit], [1, limit], ls="--", lw=1.0, color="#cbd5e0", zorder=1)
    if highlight is not None and (merged["target_id"] == highlight).any():
        row = merged[merged["target_id"] == highlight].iloc[0]
        ax.scatter(
            [row["rank_transformation"]],
            [row["rank_moma"]],
            s=90,
            color="#e53e3e",
            zorder=5,
        )
        ax.annotate(
            f"{highlight}\n{int(row['rank_transformation'])} vs {int(row['rank_moma'])}",
            xy=(row["rank_transformation"], row["rank_moma"]),
            xytext=(0.08, 0.86),
            textcoords="axes fraction",
            fontsize=font["tick"],
            color="#e53e3e",
            arrowprops=dict(arrowstyle="->", color="#e53e3e", lw=1.0),
        )
    # Pearson over ranks is Spearman's rho by definition.
    rho = float(np.corrcoef(merged["rank_transformation"], merged["rank_moma"])[0, 1])
    ax.text(
        0.97,
        0.04,
        f"Spearman \u03c1 = {rho:.3f}",
        transform=ax.transAxes,
        ha="right",
        fontsize=font["tick"],
        color="#4a5568",
    )
    _style(
        ax,
        font,
        xlabel="transformation rank",
        ylabel="MOMA baseline rank",
        title="Transformation rank vs MOMA rank",
    )
    fig.tight_layout()
    return fig


def epsilon_sensitivity_figure(
    epsilon_sweep,
    *,
    top_n: int = 8,
    highlight: str | None = None,
    column_width: int = 2,
) -> Figure:
    """How each leading candidate's rank moves as epsilon changes.

    Epsilon is a flux magnitude that CMM cannot derive the way the source papers do, so it is
    chosen. This figure is what makes that choice auditable: a candidate whose rank is stable
    across the sweep was not produced by the value picked.
    """

    import pandas as pd

    frame = pd.DataFrame(epsilon_sweep)
    if frame.empty:
        raise ValueError("epsilon_sweep is empty")
    baseline = (
        frame.sort_values("epsilon")
        .groupby("target_id")["rank"]
        .first()
        .sort_values()
        .head(top_n)
        .index.tolist()
    )
    if highlight is not None and highlight not in baseline:
        baseline.append(highlight)

    fig, ax, font = _new_figure(width=6.6, height=4.4, column_width=column_width)
    for target in baseline:
        series = frame[frame["target_id"] == target].sort_values("epsilon")
        if series.empty:
            continue
        emphasised = target == highlight
        ax.plot(
            series["epsilon"],
            series["rank"],
            marker="o",
            markersize=4,
            lw=2.0 if emphasised else 1.0,
            color="#e53e3e" if emphasised else None,
            zorder=5 if emphasised else 2,
            label=target,
        )
    ax.set_xscale("log")
    ax.invert_yaxis()  # rank 1 at the top, where a reader expects the best candidate
    ax.legend(frameon=False, fontsize=font["tick"], ncol=2)
    _style(
        ax,
        font,
        xlabel="epsilon (model flux units)",
        ylabel="rank",
        title="Rank sensitivity to the chosen epsilon",
    )
    fig.tight_layout()
    return fig
