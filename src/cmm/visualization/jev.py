"""Figures for a JEV run: the game's progress, and why each move was chosen.

Separate from :mod:`cmm.visualization.figures` because a JEV run is a different kind of
object — a sequence of moves rather than a scan or a screen — and because keeping it apart
means the existing figure module is untouched by this feature.

The flux *map* is not here. A JEV run's flux distribution is an ordinary flux distribution,
so the desktop app redraws it with the existing
:func:`~cmm.visualization.figures.escher_flux_map` and
:func:`~cmm.visualization.figures.network_flux_map` once per tick. Adding a third renderer
would be a second way to draw the same thing.

Both figures follow the conventions of the main figure module: built on
``matplotlib.figure.Figure`` with no pyplot global state, headless-safe, 300 DPI, and the
Okabe-Ito colour-blind-safe palette.
"""

from __future__ import annotations

from collections.abc import Sequence

from matplotlib.figure import Figure

from cmm.visualization.figures import PALETTE, _new_figure, _style

#: How each tick outcome is marked on the progress figure. Colours are chosen so the three
#: kinds of event a reader looks for — a change that stuck, one the rules rejected, and a
#: move that only looked — are distinguishable without relying on hue alone.
_OUTCOME_STYLE: dict[str, tuple[str, str, str]] = {
    "applied": (PALETTE[2], "o", "intervention applied"),
    "undone_by_agent": (PALETTE[4], "v", "withdrawn by the agent"),
    "reverted_infeasible": (PALETTE[1], "X", "reverted: model infeasible"),
    "reverted_growth_floor": (PALETTE[1], "s", "reverted: growth floor"),
    "not_applicable": (PALETTE[6], "P", "move unavailable"),
    "scan": (PALETTE[5], ".", "analysis run, model unchanged"),
    "end_round": (PALETTE[6], "*", "round ended"),
}


def jev_progress_figure(
    result,
    *,
    width: float = 7.5,
    height: float = 4.5,
    column_width: int = 2,
) -> Figure:
    """Product flux and growth rate over the run, one point per tick.

    The two quantities share an x-axis and use twin y-axes, because the whole question a
    reader has about a run is whether the product went up *without* growth going down, and
    that is a question about the two curves' shapes against each other.

    Wild-type levels are drawn as reference lines rather than as the first data point: the
    wild type is not a move, and plotting it as one would imply the agent chose it.
    """

    ticks: Sequence = tuple(result.ticks)
    fig, ax, font = _new_figure(width, height, column_width)

    if not ticks:
        _style(
            ax,
            font,
            xlabel="tick",
            ylabel="product flux (mmol gDW$^{-1}$ h$^{-1}$)",
            title="JEV run produced no ticks",
        )
        ax.text(
            0.5,
            0.5,
            "no moves were played",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=font["label"],
            color="#666666",
        )
        return fig

    x = list(range(1, len(ticks) + 1))
    product = [float(tick.product_flux) for tick in ticks]
    growth = [float(tick.growth) for tick in ticks]

    ax.plot(
        x,
        product,
        color=PALETTE[0],
        linewidth=1.8,
        zorder=2,
        label=f"{result.config.product} flux",
    )
    ax.axhline(
        result.wild_type_product_flux,
        color=PALETTE[0],
        linestyle=":",
        linewidth=1.2,
        alpha=0.7,
        zorder=1,
    )

    seen: set[str] = set()
    for index, tick in zip(x, ticks):
        colour, marker, label = _OUTCOME_STYLE.get(
            tick.outcome, (PALETTE[6], "o", tick.outcome)
        )
        ax.scatter(
            index,
            tick.product_flux,
            color=colour,
            marker=marker,
            s=46 if tick.outcome != "scan" else 18,
            zorder=3,
            edgecolors="white",
            linewidths=0.5,
            label=label if label not in seen else None,
        )
        seen.add(label)

    twin = ax.twinx()
    twin.plot(x, growth, color=PALETTE[1], linewidth=1.4, linestyle="--", zorder=2)
    twin.axhline(
        result.config.growth_floor,
        color=PALETTE[1],
        linestyle=":",
        linewidth=1.2,
        alpha=0.8,
    )
    twin.set_ylabel(
        "growth rate (h$^{-1}$), dashed", fontsize=font["label"], color=PALETTE[1]
    )
    twin.tick_params(labelsize=font["tick"], colors=PALETTE[1])
    twin.spines["top"].set_visible(False)

    # Round boundaries, so a reader can see where checkpoints fell.
    boundary = 0
    for record in result.rounds:
        boundary += record.n_ticks
        if 0 < boundary < len(ticks):
            ax.axvline(boundary + 0.5, color="#cccccc", linewidth=0.8, zorder=0)

    _style(
        ax,
        font,
        xlabel="tick (move number)",
        ylabel=f"{result.config.product} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        title="JEV run: product and growth after every move",
    )
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, fontsize=font["legend"], loc="best", framealpha=0.9)
    fig.tight_layout()
    return fig


def jev_decision_figure(
    tick,
    *,
    top_n: int = 12,
    width: float = 6.5,
    height: float = 4.5,
    column_width: int = 2,
) -> Figure:
    """What the agent nearly chose: its probability over the whole board for one move.

    A ``choice`` answer carries a probability for every option, not just the winner, so a
    single move records a complete ranking of the board. Drawing it turns "the agent picked
    FRD7" into "the agent picked FRD7 over FUM by 0.63 to 0.14", which is the difference
    between an assertion and evidence — and it is often the runner-up that is worth a look.
    """

    ranking = list(tick.target_ranking)[:top_n]
    fig, ax, font = _new_figure(width, height, column_width)

    if not ranking:
        _style(ax, font, xlabel="probability", ylabel="", title="no ranking recorded")
        ax.text(
            0.5,
            0.5,
            "this move carried no probability distribution",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=font["label"],
            color="#666666",
        )
        return fig

    names = [name for name, _ in ranking][::-1]
    values = [value for _, value in ranking][::-1]
    colours = [PALETTE[2] if name == tick.target else PALETTE[0] for name in names]

    positions = range(len(names))
    ax.barh(list(positions), values, color=colours, height=0.7)
    ax.set_yticks(list(positions))
    ax.set_yticklabels(names, fontsize=font["tick"])
    for position, value in zip(positions, values):
        ax.text(
            value + 0.01,
            position,
            f"{value:.2f}",
            va="center",
            fontsize=font["tick"],
            color="#444444",
        )

    chosen = tick.action or "no move"
    _style(
        ax,
        font,
        xlabel="probability the agent assigned to acting here",
        ylabel="",
        title=(
            f"Round {tick.round_index}, tick {tick.tick_index}: "
            f"chose {tick.target} → {chosen}"
        ),
    )
    ax.set_xlim(0, max(values) * 1.18 + 0.02)
    ax.grid(True, axis="x", alpha=0.25, linewidth=0.6)
    ax.grid(False, axis="y")
    fig.tight_layout()
    return fig


__all__ = ["jev_decision_figure", "jev_progress_figure"]
