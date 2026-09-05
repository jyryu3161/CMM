"""Shared result types for ranked-intervention services.

FSEOF, OptKnock, and revert-metabolism all produce an ordered list of intervention
targets with a score. They share :class:`TargetRanking` so downstream code (export,
GUI tables, comparison) handles one type.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
import math

import pandas as pd


@dataclass(frozen=True)
class TargetScore:
    """One ranked target."""

    target_id: str
    score: float
    detail: Mapping[str, float] = field(default_factory=dict)
    #: Human-readable label for ``target_id`` — a gene symbol or a reaction name, taken from
    #: the model rather than from any imported table. Empty when the model names nothing, which
    #: is why it never replaces ``target_id``: names are not unique and not always present,
    #: whereas the id is the key that provenance and reproduction are written against.
    target_name: str = ""
    #: Optional solve status; absent for services that only return scored targets.
    status: str | None = None

    def __post_init__(self) -> None:
        if not self.target_id:
            raise ValueError("target_id must not be empty")
        object.__setattr__(self, "target_name", str(self.target_name or ""))
        score = float(self.score)
        if math.isnan(score):
            raise ValueError("target score must not be NaN")
        detail = {str(k): float(v) for k, v in self.detail.items()}
        if any(math.isnan(value) for value in detail.values()):
            raise ValueError("target score detail must not contain NaN")
        object.__setattr__(self, "score", score)
        object.__setattr__(self, "detail", detail)


@dataclass(frozen=True)
class TargetRanking:
    """An ordered collection of scored targets with run metadata."""

    method: str
    targets: tuple[TargetScore, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "targets", tuple(self.targets))
        object.__setattr__(self, "metadata", dict(self.metadata))

    def __iter__(self):
        return iter(self.targets)

    def __len__(self) -> int:
        return len(self.targets)

    def sorted(self, descending: bool = True) -> TargetRanking:
        """Order by score then target_id, retaining failed solves after scored targets."""

        ordered = sorted(
            self.targets,
            key=lambda t: (
                t.status is not None and t.status != "optimal",
                -t.score if descending else t.score,
                t.target_id,
            ),
        )
        return TargetRanking(
            method=self.method, targets=tuple(ordered), metadata=self.metadata
        )

    def top(self, n: int) -> tuple[TargetScore, ...]:
        return self.sorted().targets[:n]

    def best(self) -> TargetScore | None:
        ordered = self.sorted().targets
        return next(
            (target for target in ordered if target.status in (None, "optimal")),
            None,
        )

    def to_records(self) -> list[dict]:
        records: list[dict] = []
        # The column only appears when something was actually named, so a ranking over a model
        # that names nothing exports exactly the table it did before.
        named = any(t.target_name for t in self.targets)
        has_status = any(t.status is not None for t in self.targets)
        for rank, t in enumerate(self.sorted().targets, start=1):
            row: dict = {"rank": rank, "target_id": t.target_id}
            if named:
                row["target_name"] = t.target_name
            row["score"] = t.score
            row.update(t.detail)
            if has_status:
                row["status"] = t.status
            records.append(row)
        return records

    def to_frame(self) -> pd.DataFrame:
        """Deterministic export table (rank, target_id, score, then sorted detail cols)."""

        records = self.to_records()
        base_cols = ["rank", "target_id", "score"]
        if any("target_name" in r for r in records):
            base_cols.insert(2, "target_name")
        detail_cols = sorted({k for r in records for k in r} - set(base_cols))
        frame = pd.DataFrame(records, columns=base_cols + detail_cols)
        return frame

    @classmethod
    def from_scores(
        cls,
        method: str,
        scores: Iterable[tuple[str, float]] | Mapping[str, float],
        metadata: Mapping[str, object] | None = None,
    ) -> TargetRanking:
        items = scores.items() if isinstance(scores, Mapping) else scores
        targets = tuple(TargetScore(target_id=str(k), score=float(v)) for k, v in items)
        return cls(
            method=method, targets=targets, metadata=dict(metadata or {})
        ).sorted()
