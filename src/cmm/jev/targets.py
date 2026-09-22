"""Every target the run touched, with what recommends it and what argues against it.

A run's headline is one design. That is not the whole of what it learned: over several rounds
the agent measures dozens of reactions, tries a good few, and has some of them taken away by
the rules — and all of that is evidence about targets a reader may want to weigh differently
than the agent did. A reader whose strain has to hold a higher growth rate, or who cannot
delete three isozymes, wants the second-best target and the reason it came second.

So this module assembles one row per reaction the run ever acted on or measured, and states
its case both ways. Nothing here is a new computation: every number was already produced by
the game and is being read back out. That matters for the same reason the rest of CMM records
provenance — a pro or con invented at report time is an opinion, and an opinion in a results
table is indistinguishable from a measurement.

**Pros and cons are evidence, not a score.** They are deliberately not collapsed into a
ranking: "raises the product by 0.035" and "needs all three isozymes deleted" are not
commensurable, and the trade between them belongs to whoever is building the strain.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

#: Below this a measured change in product flux is reported as no change rather than as a
#: gain, matching the tolerance the engine itself uses when it judges a move.
_NEGLIGIBLE = 1e-6


@dataclass(frozen=True)
class TargetReport:
    """One reaction the run considered, and the case for and against editing it."""

    reaction_id: str
    name: str
    subsystem: str
    genes: tuple[str, ...]
    gene_edit: str
    side_effects: tuple[str, ...]
    #: Measured product change when CMM deleted it, and when it halved it, on whatever design
    #: was standing at the time. ``None`` means the screen never reached it.
    deletion_gain: float | None
    knockdown_gain: float | None
    essential: bool | None
    fseof_slope: float | None
    design_note: str
    #: What happened when the agent actually played it, if it did.
    attempts: int
    applied: int
    reverted: tuple[str, ...]
    #: True when the reaction is in the best design the run found.
    in_best_design: bool
    best_design_mode: str
    pros: tuple[str, ...]
    cons: tuple[str, ...]

    def to_row(self) -> dict[str, object]:
        return {
            "reaction_id": self.reaction_id,
            "name": self.name,
            "subsystem": self.subsystem,
            "genes": ";".join(self.genes),
            "gene_edit": self.gene_edit,
            "n_side_effects": len(self.side_effects),
            "side_effects": ";".join(self.side_effects),
            "deletion_gain": self.deletion_gain,
            "knockdown_gain": self.knockdown_gain,
            "essential": self.essential,
            "fseof_slope": self.fseof_slope,
            "attempts": self.attempts,
            "applied": self.applied,
            "times_reverted": len(self.reverted),
            "in_best_design": self.in_best_design,
            "best_design_mode": self.best_design_mode,
            "pros": " | ".join(self.pros),
            "cons": " | ".join(self.cons),
        }


@dataclass
class _Accumulator:
    """Mutable working state for one reaction while the run's records are walked."""

    evidence: object | None = None
    attempts: int = 0
    applied: int = 0
    reverted: list[str] = field(default_factory=list)


def build_target_reports(result) -> tuple[TargetReport, ...]:
    """Assemble the per-target report from a finished :class:`~cmm.jev.engine.JevResult`.

    Targets are ordered by how much they are worth: the best design's members first, then by
    the largest measured gain. A reaction the run only measured and never played still gets a
    row — "CMM checked this and it does not pay" is a result, and leaving it out would make
    the table a record of the agent's attention rather than of the evidence.
    """

    seen: dict[str, _Accumulator] = {}

    for evidence in getattr(result, "candidates_seen", ()) or ():
        entry = seen.setdefault(evidence.reaction_id, _Accumulator())
        # Later sightings carry later measurements, which are the ones against the design
        # that was actually standing when the round ended.
        entry.evidence = evidence

    for tick in result.ticks:
        intervention = tick.intervention
        if intervention is None:
            continue
        entry = seen.setdefault(intervention.reaction_id, _Accumulator())
        entry.attempts += 1
        if tick.outcome == "applied":
            entry.applied += 1
        elif tick.outcome.startswith("reverted"):
            entry.reverted.append(tick.reason)

    best = {i.reaction_id: i for i in result.best_interventions}
    reports = [
        _report(reaction_id, entry, best.get(reaction_id))
        for reaction_id, entry in seen.items()
    ]
    reports.sort(
        key=lambda report: (
            not report.in_best_design,
            -max(
                report.deletion_gain if report.deletion_gain is not None else 0.0,
                report.knockdown_gain if report.knockdown_gain is not None else 0.0,
            ),
            report.reaction_id,
        )
    )
    return tuple(reports)


def _report(reaction_id: str, entry: _Accumulator, chosen) -> TargetReport:
    evidence = entry.evidence
    deletion = getattr(evidence, "deletion_gain", None)
    knockdown = getattr(evidence, "knockdown_gain", None)
    essential = getattr(evidence, "essential", None)
    slope = getattr(evidence, "fseof_slope", None)
    side_effects = tuple(getattr(evidence, "side_effects", ()) or ())
    if chosen is not None and not side_effects:
        side_effects = chosen.side_effects
    genes = tuple(getattr(evidence, "genes", ()) or ())
    if chosen is not None and not genes:
        genes = chosen.gene_names or chosen.genes

    pros, cons = _weigh(
        deletion=deletion,
        knockdown=knockdown,
        essential=essential,
        slope=slope,
        side_effects=side_effects,
        genes=genes,
        design_note=str(getattr(evidence, "design_note", "") or ""),
        literature=str(getattr(evidence, "literature", "") or ""),
        entry=entry,
        chosen=chosen,
    )
    return TargetReport(
        reaction_id=reaction_id,
        name=str(getattr(evidence, "name", "") or reaction_id),
        subsystem=str(getattr(evidence, "subsystem", "") or ""),
        genes=genes,
        gene_edit=str(getattr(evidence, "gene_edit", "") or ""),
        side_effects=side_effects,
        deletion_gain=deletion,
        knockdown_gain=knockdown,
        essential=essential,
        fseof_slope=slope,
        design_note=str(getattr(evidence, "design_note", "") or ""),
        attempts=entry.attempts,
        applied=entry.applied,
        reverted=tuple(entry.reverted),
        in_best_design=chosen is not None,
        best_design_mode=chosen.mode if chosen is not None else "",
        pros=pros,
        cons=cons,
    )


def _weigh(
    *,
    deletion: float | None,
    knockdown: float | None,
    essential: bool | None,
    slope: float | None,
    side_effects: Sequence[str],
    genes: Sequence[str],
    design_note: str,
    literature: str,
    entry: _Accumulator,
    chosen,
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """The case for and against one target, assembled only from what the run measured."""

    pros: list[str] = []
    cons: list[str] = []

    if chosen is not None:
        pros.append(f"in the best design this run found, as a {chosen.mode}")

    for label, gain in (("deleting it", deletion), ("halving it", knockdown)):
        if gain is None:
            continue
        if gain > _NEGLIGIBLE:
            pros.append(f"CMM measured {label} raising the product by {gain:+.4g}")
        elif gain < -_NEGLIGIBLE:
            cons.append(f"CMM measured {label} lowering the product by {gain:+.4g}")
        else:
            cons.append(f"CMM measured {label} changing the product by 0")

    if essential is True:
        cons.append("deleting it stops growth: only the knockdown is available")
    elif essential is False:
        pros.append("the cell grows without it, so the deletion is available")

    if slope is not None:
        if slope < -_NEGLIGIBLE:
            pros.append(
                f"FSEOF: its flux falls ({slope:+.3g}) as the product is forced up, so the "
                "product pathway does not need it"
            )
        elif slope > _NEGLIGIBLE:
            cons.append(
                f"FSEOF: its flux rises ({slope:+.3g}) with the product, so cutting it works "
                "against the pathway"
            )

    if design_note:
        pros.append(design_note)

    if len(genes) > 1:
        cons.append(f"the edit needs {len(genes)} genes ({', '.join(genes)}), not one")
    elif len(genes) == 1:
        pros.append(f"a single-gene edit ({genes[0]})")

    if side_effects:
        cons.append(
            f"the same gene(s) also run {', '.join(side_effects)}, which the edit stops too"
        )

    unmodelled = tuple(getattr(chosen, "unmodelled", ()) or ()) if chosen else ()
    if unmodelled:
        cons.append(
            f"the edit also touches {', '.join(unmodelled)}, carrying no wild-type flux, so "
            "its effect there is outside what this model can say"
        )

    for reason in dict.fromkeys(entry.reverted):
        cons.append(f"CMM refused it: {reason}")

    if entry.attempts and not entry.applied and not entry.reverted:
        cons.append("the agent proposed it but it could not be applied")

    if literature:
        pros.append(
            "the literature lookup returned published evidence; see literature.csv"
        )

    return tuple(pros), tuple(cons)


def targets_frame(reports: Sequence[TargetReport]) -> pd.DataFrame:
    """The report as a table, in reporting order."""

    return pd.DataFrame(
        [report.to_row() for report in reports],
        columns=[
            "reaction_id",
            "name",
            "subsystem",
            "genes",
            "gene_edit",
            "n_side_effects",
            "side_effects",
            "deletion_gain",
            "knockdown_gain",
            "essential",
            "fseof_slope",
            "attempts",
            "applied",
            "times_reverted",
            "in_best_design",
            "best_design_mode",
            "pros",
            "cons",
        ],
    )


def targets_summary(reports: Sequence[TargetReport]) -> Mapping[str, object]:
    """A short, honest description of what the target table contains."""

    measured = [r for r in reports if r.deletion_gain is not None]
    paying = [
        r
        for r in measured
        if max(
            r.deletion_gain or 0.0,
            r.knockdown_gain if r.knockdown_gain is not None else 0.0,
        )
        > _NEGLIGIBLE
    ]
    return {
        "n_targets": len(reports),
        "n_measured": len(measured),
        "n_that_pay": len(paying),
        "n_in_best_design": sum(1 for r in reports if r.in_best_design),
        "n_with_side_effects": sum(1 for r in reports if r.side_effects),
        "note": (
            "Pros and cons are read back from what the run measured; they are not scored "
            "against each other, because a product gain and a three-gene edit are not the "
            "same kind of quantity. Every gain was measured on whichever design was standing "
            "when the screen ran, so two rows are not necessarily comparable to each other."
        ),
    }


__all__ = [
    "TargetReport",
    "build_target_reports",
    "targets_frame",
    "targets_summary",
]
