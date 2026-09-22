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
import html


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


def render_agent_report(result) -> str:
    """The whole run as one HTML page: design, rounds, targets, baselines, provenance."""

    summary = result.summary()
    config = result.config
    product = str(summary["product"])
    wild = float(summary["wild_type_product_flux"])
    best = float(summary["best_product_flux"])

    if summary["beat_wild_type"]:
        fold = summary["fold_improvement"]
        gain = (
            f" ({fold:.1f}× the wild type)"
            if isinstance(fold, (int, float))
            else " (the wild type made none)"
        )
        lede = (
            f"<p class='lede'><b>{_escape(product)} rose from {wild:.4g} to {best:.4g} "
            f"mmol gDW&#8315;&#185; h&#8315;&#185;{gain}</b>, at a growth rate of "
            f"{float(summary['best_growth']):.4g} h&#8315;&#185;.</p>"
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
        ["Method", "Design", "Product", "Growth", "Deterministic", "Note"],
        [
            (
                row.method,
                "; ".join(row.design) or "—",
                f"{row.product_flux:.4g}",
                f"{row.growth:.4g}",
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

<h2>Measured against the deterministic methods</h2>
<p>Same model, same condition, same growth floor, every design applied and solved the same
way.</p>
{baselines}
{f"<p><b>{_escape(verdict)}</b></p>" if verdict else ""}

<h2>What the run cost</h2>
<p>{usage.get("calls", 0)} decision(s), {usage.get("input_tokens", 0)} input and
{usage.get("output_tokens", 0)} output tokens, ${float(usage.get("cost_usd", 0.0)):.4f}.</p>

<h2>Notes from the run</h2>
{_lines(summary["notes"])}

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
