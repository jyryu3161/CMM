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


#: Plot labels for the methods whose names are sentences. A name that reads well in a table
#: column runs across half the axes here, and two of them overlap into one unreadable line.
_PLOT_LABELS = {
    "best deterministic design + one knockdown (exhaustive)": "proven + 1 knockdown",
    "best amplification (outside the vocabulary)": "best amplification",
    "best single gene deletion": "best single deletion",
}


def _short_method(name: str) -> str:
    """A baseline's name shortened for a plot, where a table column's worth of words will not fit."""

    if name in _PLOT_LABELS:
        return _PLOT_LABELS[name]
    return name.split(" (")[0].strip()


def _merge_labels(
    points: Sequence[tuple[float, float, str]], *, span_x: float, span_y: float
) -> list[tuple[float, float, str]]:
    """One label per place on the plane, not one per thing plotted.

    Designs land on the same point all the time and it is not a coincidence: OptKnock and
    RobustKnock return the same deletions here, and two rounds forced apart by the cut can
    still meet at the same phenotype. Drawn naively the names print on top of one another and
    the figure claims one illegible thing where it should say two legible ones.

    Points are the same place when they are within a thousandth of the plotted range on both
    axes — a distance chosen to be smaller than a marker, so nothing visually distinct is
    ever merged.
    """

    tol_x = max(span_x, 1e-9) / 1000.0
    tol_y = max(span_y, 1e-9) / 1000.0
    merged: list[tuple[float, float, list[str]]] = []
    for x, y, label in points:
        for index, (mx, my, names) in enumerate(merged):
            if abs(mx - x) <= tol_x and abs(my - y) <= tol_y:
                if label not in names:
                    names.append(label)
                merged[index] = (mx, my, names)
                break
        else:
            merged.append((x, y, [label]))
    return [(x, y, ", ".join(names)) for x, y, names in merged]


def jev_design_space_figure(
    result,
    *,
    width: float = 7.5,
    height: float = 5.0,
    column_width: int = 2,
) -> Figure:
    """Every design the run produced, placed on the growth-versus-product plane.

    This is the figure that makes a multi-round run legible as a portfolio rather than as one
    answer. Each round ends on a design; each design is a point; and the question a reader
    actually has — *what does this cost me in growth, and is there a cheaper one?* — is a
    question about where those points sit relative to each other.

    The envelope is the backdrop and is the reason the points mean anything. It is the
    projection of the flux cone onto this plane (Burgard 2003), so it bounds what any design
    whatsoever could reach: a point near the frontier has little left to win, and a point well
    inside it does. Without it, "9.95 at growth 0.055" is a pair of numbers with nothing to be
    measured against.

    The deterministic methods are plotted in the same axes, as squares. A reader comparing the
    agent to OptKnock should be able to do it by looking, not by holding two tables side by
    side — and where the agent's point sits *below and left* of a deterministic one, that is
    worth seeing plainly rather than discovering in a footnote.

    Every design is drawn twice: a filled marker at what it must make, and a hollow one at what
    it could, joined by a line. The gap between them is the design's own uncertainty, and it is
    the thing this plane exists to show — two designs whose best cases coincide can have
    nothing in common once you ask what each one *guarantees*. A design drawn as a single
    point is a design quoted at its best case, which is not what a strain does.
    """

    fig, ax, font = _new_figure(width, height, column_width)
    product = str(result.config.product)

    envelope = tuple(getattr(result, "envelope", ()) or ())
    if envelope:
        # The frontier: the most product reachable at each growth rate, and the least. Drawn
        # as a filled band because the feasible region is what it delimits, and a reader
        # should see a design sitting *inside* a region rather than near a line.
        flux = [point[0] for point in envelope]
        low = [point[1] for point in envelope]
        high = [point[2] for point in envelope]
        ax.fill_betweenx(
            flux, low, high, color=PALETTE[0], alpha=0.10, zorder=0, linewidth=0
        )
        ax.plot(high, flux, color=PALETTE[0], linewidth=1.3, alpha=0.65, zorder=1)
        ax.plot(low, flux, color=PALETTE[0], linewidth=1.0, alpha=0.4, zorder=1)
        ax.plot(
            [],
            [],
            color=PALETTE[0],
            linewidth=1.3,
            alpha=0.65,
            label="feasible envelope",
        )

    ax.scatter(
        [float(result.wild_type_growth)],
        [float(result.wild_type_product_flux)],
        marker="*",
        s=170,
        color="#555555",
        zorder=5,
        label="wild type",
    )

    # The deterministic methods, so the comparison is something a reader can see rather than
    # something they have to assemble from two tables.
    baseline_labels: list[tuple[float, float, str]] = []
    for row in getattr(result, "baselines", ()) or ():
        if row.method == "JEV agent" or not row.design:
            continue
        if row.product_flux != row.product_flux or row.growth != row.growth:
            continue  # a method that could not run has no point on this plane
        guaranteed = getattr(row, "guaranteed_product", None)
        if guaranteed is not None and abs(guaranteed - row.product_flux) > 1e-6:
            ax.plot(
                [float(row.growth), float(row.growth)],
                [float(guaranteed), float(row.product_flux)],
                color=PALETTE[6],
                linewidth=1.0,
                alpha=0.55,
                zorder=3,
            )
        ax.scatter(
            [float(row.growth)],
            [float(row.product_flux)],
            marker="s",
            s=46,
            facecolor="white",
            edgecolor=PALETTE[6],
            linewidth=1.2,
            zorder=4,
        )
        if guaranteed is not None:
            ax.scatter(
                [float(row.growth)],
                [float(guaranteed)],
                marker="s",
                s=46,
                color=PALETTE[6],
                edgecolor=PALETTE[6],
                linewidth=1.2,
                zorder=4,
            )
        baseline_labels.append(
            (
                float(row.growth),
                float(guaranteed if guaranteed is not None else row.product_flux),
                _short_method(row.method),
            )
        )

    rounds = tuple(result.rounds)

    def _score(record) -> float:
        """What a round is worth — the same quantity the run ranked it on."""

        guaranteed = getattr(record, "guaranteed_product", None)
        return float(guaranteed if guaranteed is not None else record.product_flux)

    best = max((_score(record) for record in rounds), default=float("-inf"))
    round_labels: list[tuple[float, float, str]] = []
    for record in rounds:
        score = _score(record)
        is_best = score >= best - 1e-9
        if abs(score - record.product_flux) > 1e-6:
            ax.plot(
                [float(record.growth), float(record.growth)],
                [score, float(record.product_flux)],
                color=PALETTE[2] if is_best else PALETTE[4],
                linewidth=1.0,
                alpha=0.6,
                zorder=5,
            )
            ax.scatter(
                [float(record.growth)],
                [float(record.product_flux)],
                marker="o",
                s=60 if is_best else 40,
                facecolor="white",
                edgecolor=PALETTE[2] if is_best else PALETTE[4],
                linewidth=1.0,
                zorder=6,
            )
        ax.scatter(
            [float(record.growth)],
            [score],
            marker="o",
            s=110 if is_best else 70,
            color=PALETTE[2] if is_best else PALETTE[4],
            edgecolor="#23313f",
            linewidth=1.0 if is_best else 0.6,
            zorder=6,
        )
        round_labels.append((float(record.growth), score, f"R{record.round_index}"))
    if rounds:
        ax.scatter([], [], marker="o", s=70, color=PALETTE[4], label="a round's design")
        ax.scatter(
            [], [], marker="o", s=110, color=PALETTE[2], label="best design found"
        )
    ax.scatter(
        [],
        [],
        marker="s",
        s=46,
        facecolor="white",
        edgecolor=PALETTE[6],
        linewidth=1.2,
        label="deterministic method",
    )
    # Filled is what a design must make and hollow what it could, so a reader never has to be
    # told which of the two numbers a point is.
    ax.scatter([], [], marker="o", s=70, color="#6b7a88", label="guaranteed product")
    ax.scatter(
        [],
        [],
        marker="o",
        s=50,
        facecolor="white",
        edgecolor="#6b7a88",
        linewidth=1.0,
        label="best case (pFBA)",
    )

    floor = float(result.config.growth_floor)
    if floor > 0:
        # Not decoration: every point left of this line was refused by CMM, so the line is
        # the edge of what the run was allowed to keep.
        ax.axvline(
            floor,
            color=PALETTE[1],
            linestyle="--",
            linewidth=1.1,
            alpha=0.8,
            zorder=2,
            label=f"growth floor ({floor:g})",
        )

    _style(
        ax,
        font,
        xlabel="growth rate (h$^{-1}$)",
        ylabel=f"{product} flux (mmol gDW$^{{-1}}$ h$^{{-1}}$)",
        title="What each design guarantees, and what it costs in growth",
    )
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0)

    # Labels last, once the axes know their range, and merged so two designs meeting at one
    # phenotype are named once instead of printing over each other.
    span_x = max(ax.get_xlim()[1] - ax.get_xlim()[0], 1e-9)
    span_y = max(ax.get_ylim()[1] - ax.get_ylim()[0], 1e-9)
    for x, y, label in _merge_labels(round_labels, span_x=span_x, span_y=span_y):
        ax.annotate(
            label,
            (x, y),
            textcoords="offset points",
            xytext=(9, 6),
            fontsize=font["tick"],
            fontweight="bold",
            color="#23313f",
            zorder=7,
        )
    # Baseline labels are stacked rather than merged when two designs land near — but not on —
    # the same point. Merging them would claim they are one design; printing both at the same
    # offset makes two legible names into one unreadable line, which is how the exhaustive
    # sweep and the amplification probe collided at growth 0.053.
    placed: list[tuple[float, float, int]] = []
    for x, y, label in _merge_labels(baseline_labels, span_x=span_x, span_y=span_y):
        row = 0
        while any(
            # Generous on x, because what collides is the *text*, which runs to the right of
            # its point and is far wider than the marker; tight on y, because two labels at
            # different heights do not overlap however close their points are in growth.
            abs(px - x) <= span_x / 4.0 and abs(py - y) <= span_y / 12.0 and prow == row
            for px, py, prow in placed
        ):
            row += 1
        placed.append((x, y, row))
        ax.annotate(
            label,
            (x, y),
            textcoords="offset points",
            xytext=(9, -11 - 12 * row),
            fontsize=font["tick"] - 1,
            color="#5a6b7c",
            zorder=4,
        )

    # Pinned rather than "best". Matplotlib places a best-fit legend around the plotted
    # artists and knows nothing about the annotations, so it kept landing on the method names
    # in the upper right. The lower left is the one corner a design never occupies: a point
    # there makes little product *and* grows slowly, which is dominated by the wild type.
    ax.legend(fontsize=font["tick"], frameon=False, loc="lower left")
    fig.set_layout_engine("constrained")
    return fig


__all__ = [
    "jev_decision_figure",
    "jev_design_space_figure",
    "jev_progress_figure",
]
