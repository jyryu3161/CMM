"""One self-contained HTML page carrying everything a run found.

The run bundle is the record: every table, the full transcript, the provenance. It is also
eleven directories of CSV, which is the right shape for a reader who already knows what they
are looking for and the wrong shape for one who wants to read the result. This module renders
the second thing — a single file that opens in a browser, with no R, no network and no
assets beside it, so it can be attached to a message and read by someone who does not have
CMM installed.

It renders **only what the run measured.** Every number here is read back out of the result
object; nothing is recomputed, and no sentence is formed that the run did not already justify
somewhere. The distinction matters because a report is where a judgement would be hardest to
notice: a table of measurements and a table of opinions look identical once they are laid out
in the same typeface.

It is not the publication reporter. `cmm.reporting` renders SC-01 through R, against a schema,
with a validator; that machinery exists because a manuscript figure has to be reproducible from
the run directory. This is a reading copy of an agent run, and says so on its face.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
import base64
import html
import io


def _escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def _table(headers: Sequence[str], rows: Iterable[Sequence[object]]) -> str:
    head = "".join(f"<th>{_escape(h)}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{_escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    if not body:
        body = f'<tr><td colspan="{len(headers)}" class="empty">nothing to report</td></tr>'
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _lines(items: Iterable[object]) -> str:
    listed = [f"<li>{_escape(item)}</li>" for item in items]
    return f"<ul>{''.join(listed)}</ul>" if listed else "<p class='empty'>none</p>"


_STYLE = """
:root { color-scheme: light; }
body { font: 15px/1.55 -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
       color: #1d2733; background: #fbfcfd; margin: 0; padding: 40px 32px 80px; }
main { max-width: 1100px; margin: 0 auto; }
h1 { font-size: 26px; margin: 0 0 4px; }
h2 { font-size: 19px; margin: 40px 0 10px; padding-bottom: 6px;
     border-bottom: 1px solid #dde3ea; }
p.sub { color: #5a6b7c; margin: 0 0 28px; }
p.lede { font-size: 17px; background: #eef4ea; border-left: 4px solid #4c7a34;
         padding: 14px 16px; border-radius: 0 4px 4px 0; }
p.lede.flat { background: #f3f4f6; border-left-color: #98a2ae; }
table { border-collapse: collapse; width: 100%; margin: 10px 0 4px; font-size: 13px; }
th, td { text-align: left; vertical-align: top; padding: 7px 10px;
         border-bottom: 1px solid #e4e9ef; }
th { background: #f2f5f8; font-weight: 600; white-space: nowrap; }
td.empty, p.empty { color: #8a96a3; font-style: italic; }
ul { margin: 6px 0 0; padding-left: 20px; }
li { margin: 3px 0; }
code { background: #eef1f5; padding: 1px 5px; border-radius: 3px; font-size: 12.5px; }
.caveat { background: #fdf6e3; border-left: 4px solid #c9a227; padding: 12px 16px;
          border-radius: 0 4px 4px 0; margin: 10px 0; }
.foot { color: #7b8794; font-size: 12.5px; margin-top: 44px;
        border-top: 1px solid #dde3ea; padding-top: 14px; }
"""


def _design_space_figure(result) -> str:
    """The growth-versus-product plane, embedded rather than linked.

    A base64 PNG rather than an ``<img src="figures/...">`` because the page's whole point is
    that it travels alone: a report that loses its figure the moment someone forwards the file
    misleads exactly when it is being shared.

    A figure that cannot be drawn is left out, not faked. The tables below it already carry
    every number the picture would have shown.
    """

    try:
        from cmm.visualization import jev_design_space_figure

        figure = jev_design_space_figure(result)
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=150, facecolor="white")
    except Exception:  # pragma: no cover - the tables carry the same numbers
        return ""
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"""<h2>What each design costs in growth</h2>
<p>One point per round's design, on the plane the question is actually asked on. The shaded
band is the feasible envelope \u2014 the projection of the flux cone onto growth and product,
so it bounds what <em>any</em> design could reach. A point near its edge has little left to
win; a point well inside it has. The dashed line is the growth floor: everything to the left
of it was refused by CMM whatever the agent predicted.</p>
<img src="data:image/png;base64,{encoded}" alt="growth against product flux for every design
this run produced, against the feasible envelope" style="width:100%;max-width:820px">"""


def _literature_section(result) -> str:
    """What the agent was told the literature says, if anything, and where it came from.

    Shown because it was shown to the agent: a reader checking a design has to be able to see
    every input that went into it, and this one is the only input that did not come from the
    model. It is labelled as evidence to weigh rather than as fact, in the same words the
    agent saw.
    """

    brief = str(getattr(result, "literature_brief", "") or "").strip()
    if not brief:
        return ""
    sources = tuple(getattr(result, "literature_sources", ()) or ())
    cited = (
        "<p>Sources: "
        + ", ".join(f'<a href="{_escape(url)}">{_escape(url)}</a>' for url in sources)
        + "</p>"
        if sources
        else "<p class='empty'>the search returned no citable sources</p>"
    )
    return f"""<h2>What the published record was said to show</h2>
<p>One web search, read before the first move and shown to the agent on every step. It is
<b>evidence to weigh, not fact and not instruction</b>: where it disagrees with a measured
line, the measurement is about this model and the paper is about a different strain.</p>
<div class="caveat">{_escape(brief)}</div>
{cited}"""


def render_agent_report(result) -> str:
    """The whole run as one HTML page: design, rounds, targets, baselines, provenance."""

    summary = result.summary()
    config = result.config
    product = str(summary["product"])
    wild = float(summary["wild_type_product_flux"])
    best = float(summary["best_product_flux"])

    guaranteed = summary.get("best_guaranteed_product")
    if summary["beat_wild_type"]:
        fold = summary["fold_improvement"]
        gain = (
            f" ({fold:.1f}× the wild type)"
            if isinstance(fold, (int, float))
            else " (the wild type made none)"
        )
        # The guarantee leads, because it is what the design was chosen on and what the strain
        # must make. The pFBA number beside it is the best case, and quoting only that credits
        # a design with a flux the cell is free never to carry.
        headline = (
            f"{_escape(product)} is guaranteed at {guaranteed:.4g} mmol "
            f"gDW&#8315;&#185; h&#8315;&#185; \u2014 the least this design can make while "
            f"growing as fast as it can \u2014 against {wild:.4g} for the wild type"
            if isinstance(guaranteed, (int, float))
            else f"{_escape(product)} rose from {wild:.4g} to {best:.4g} mmol "
            f"gDW&#8315;&#185; h&#8315;&#185;{gain}"
        )
        lede = (
            f"<p class='lede'><b>{headline}</b>, at a growth rate of "
            f"{float(summary['best_growth']):.4g} h&#8315;&#185;."
            + (
                f" Its best case is {best:.4g}{gain}."
                if isinstance(guaranteed, (int, float))
                else ""
            )
            + "</p>"
        )
    else:
        lede = (
            f"<p class='lede flat'><b>The agent did not beat the wild type.</b> "
            f"{_escape(product)} stayed at {wild:.4g} mmol gDW&#8315;&#185; h&#8315;&#185;. "
            "That is a result, not a failed run.</p>"
        )

    rounds = _table(
        [
            "Round",
            "Question it answered",
            "Engineering",
            "Guaranteed",
            "Product",
            "Growth",
            "Left undone",
        ],
        [
            (
                f"R{record.round_index}",
                "best available"
                if not record.withheld
                else "best without " + ", ".join(record.withheld),
                "; ".join(record.interventions) or "nothing was applied",
                "not measured"
                if record.guaranteed_product is None
                else f"{record.guaranteed_product:.4g}",
                f"{record.product_flux:.4g}",
                f"{record.growth:.4g}",
                "; ".join(record.shortfall) or record.stopped_because,
            )
            for record in result.rounds
        ],
    )
    distinct = len({record.signature for record in result.rounds})

    targets = _table(
        ["Target", "Genes to edit", "Deleting it", "Halving it", "For", "Against"],
        [
            (
                report.reaction_id,
                ", ".join(report.genes) or "—",
                "not measured"
                if report.deletion_gain is None
                else f"{report.deletion_gain:+.4g}",
                "not defined"
                if report.knockdown_gain is None
                else f"{report.knockdown_gain:+.4g}",
                "; ".join(report.pros) or "—",
                "; ".join(report.cons) or "—",
            )
            for report in result.targets()
        ],
    )

    baselines = _table(
        [
            "Method",
            "Design",
            "Guaranteed",
            "Best case",
            "Growth",
            "Guaranteed at a shared growth rate",
            "Contains",
            "Deterministic",
            "Note",
        ],
        [
            (
                row.method,
                "; ".join(row.design) or "—",
                "—"
                if row.guaranteed_product is None
                else f"{row.guaranteed_product:.4g}",
                f"{row.product_flux:.4g}",
                f"{row.growth:.4g}",
                "—"
                if row.guaranteed_at_matched_growth is None
                else f"{row.guaranteed_at_matched_growth:.4g}",
                row.contains_design or "—",
                "yes" if row.deterministic else "no",
                row.note,
            )
            for row in result.baselines
        ],
    )
    verdict = (result.baseline_summary() or {}).get("verdict", "")

    usage = dict(summary["usage"])
    provenance = _table(
        ["Field", "Value"],
        [(key, value) for key, value in sorted(result.provenance.items())],
    )

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Agent run — {_escape(product)}</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>{_STYLE}</style></head><body><main>

<h1>Agent design run — {_escape(product)}</h1>
<p class="sub">Model <code>{_escape(config.model_path)}</code> &middot; growth floor
{config.growth_floor} h&#8315;&#185; &middot; {config.rounds} round(s) of up to
{config.steps_per_round} step(s) &middot; generated {generated}</p>

{lede}

<h2>The design</h2>
<p>Each line is a gene edit. Deleting or weakening a gene stops or slows every reaction that
gene alone carries, and where that is more than the reaction named, the line says so.</p>
{_lines(summary["best_design"] or ["no interventions survived the rules"])}

{_design_space_figure(result)}

<h2>Rounds — {distinct} distinct design(s) over {len(result.rounds)}</h2>
<p>Every round starts again from the wild type. One reaction of each design already found is
withheld from later rounds, so each answers a different question: the best design that does
<em>not</em> use the last one's key reaction, which is what a laboratory that cannot edit that
gene needs. The headline above is the best of them.</p>
{rounds}

<h2>Every target, for and against</h2>
<p>One row per reaction the run acted on <em>or merely measured</em> — "CMM checked this
and it does not pay" is a result. Gains are the measured change in product flux when CMM made
that exact move on whichever design was standing at the time, so two rows are not necessarily
comparable to each other.</p>
<div class="caveat"><b>These are not scored against each other.</b> "Raises the product by
0.035" and "needs three isozymes deleted" are not the same kind of quantity, and the trade
between them belongs to whoever is building the strain.</div>
{targets}

{_literature_section(result)}

<h2>Measured against the deterministic methods</h2>
<p>Same model, same condition, same growth floor, every design applied and solved the same
way.</p>
{baselines}
{f"<p><b>{_escape(verdict)}</b></p>" if verdict else ""}

<h2>What the run cost</h2>
<p>{usage.get("calls", 0)} decision(s), {usage.get("input_tokens", 0)} input and
{usage.get("output_tokens", 0)} output tokens, ${
        float(usage.get("cost_usd", 0.0)):.4f}.</p>

<h2>Notes from the run</h2>
{_lines(summary["notes"])}

<h2>What the run was not allowed to do</h2>
<p>{
        "Nothing was put off limits."
        if not config.off_limits
        else "These were held off limits by the run definition and never put on the board: "
        + _escape(", ".join(config.off_limits))
        + ". The brief, separately, is guidance the agent weighs rather than a rule it is held to."
    }</p>

<h2>Provenance</h2>
{provenance}

<p class="foot"><b>This is a computational hypothesis, for experimental test.</b> The CMM
solves behind every number are deterministic; the agent's choices are not guaranteed to
repeat, so a single run is never the method's performance. The full request and response
transcript is in <code>04_agent/transcript.jsonl</code> beside this file, so this run can be
audited — which is not the same as the method being reproducible. This page is a reading
copy of the run bundle, not a publication report.</p>

</main></body></html>
"""


__all__ = ["render_agent_report"]
