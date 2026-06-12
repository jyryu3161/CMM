"""Multi-condition omics: predict per-condition flux and compare as log-change.

Given an expression table with several condition columns (one column per condition), predict
a flux distribution per condition with E-Flux2 or LAD, then compare conditions as log2
fold-change of flux magnitude.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass

import pandas as pd
from cobra import Model

from cmm.omics.expression import OmicsFluxResult, integrate_expression


@dataclass(frozen=True)
class ConditionFluxes:
    """Per-condition predicted flux distributions from one expression table."""

    method: str
    results: dict[str, OmicsFluxResult]

    def conditions(self) -> tuple[str, ...]:
        return tuple(self.results.keys())

    def fluxes(self, condition: str) -> dict[str, float]:
        return self.results[condition].fluxes


def read_expression_table(path: str, gene_column: str | None = None) -> pd.DataFrame:
    """Read a gene-by-condition expression table; index = gene id, columns = conditions."""

    frame = pd.read_csv(path, sep=None, engine="python")
    gene_column = gene_column or frame.columns[0]
    frame = frame.set_index(gene_column)
    return frame.apply(pd.to_numeric, errors="coerce").dropna(how="all")


def predict_condition_fluxes(
    model: Model,
    expression: pd.DataFrame,
    *,
    method: str = "eflux2",
    conditions: Iterable[str] | None = None,
    **kwargs,
) -> ConditionFluxes:
    """Predict a flux distribution for each condition column with E-Flux2 or LAD."""

    columns = list(conditions) if conditions is not None else list(expression.columns)
    results: dict[str, OmicsFluxResult] = {}
    for condition in columns:
        gene_expression = {
            str(gene): float(value)
            for gene, value in expression[condition].items()
            if pd.notna(value)
        }
        results[condition] = integrate_expression(model, gene_expression, method=method, **kwargs)
    return ConditionFluxes(method=method, results=results)


def flux_log_change(
    source: dict[str, float],
    target: dict[str, float],
    *,
    reactions: Iterable[str] | None = None,
    pseudocount: float = 1e-3,
) -> dict[str, float]:
    """log2 fold-change of flux *magnitude* between two conditions (target vs source).

    Uses ``log2((|v_target| + pseudo) / (|v_source| + pseudo))`` so zero/near-zero fluxes are
    handled gracefully; the pseudocount bounds the change for reactions that switch on/off.
    """

    keys = set(reactions) if reactions is not None else set(source) | set(target)
    out: dict[str, float] = {}
    for rid in keys:
        a = abs(source.get(rid, 0.0)) + pseudocount
        b = abs(target.get(rid, 0.0)) + pseudocount
        # The zero branches are only reachable when pseudocount=0 (the default 1e-3 bounds
        # every change). Guard both off-states symmetrically: log2(0/0)=0, log2(b/0)=+inf
        # (switch on), log2(0/a)=-inf (switch off) — the last would raise a math domain error.
        if a == 0.0:
            out[rid] = 0.0 if b == 0.0 else math.inf
        elif b == 0.0:
            out[rid] = -math.inf
        else:
            out[rid] = math.log2(b / a)
    return out


def sign_flips(
    source: dict[str, float],
    target: dict[str, float],
    *,
    reactions: Iterable[str] | None = None,
    tol: float = 1e-6,
) -> list[str]:
    """Reactions whose flux direction reverses between the two conditions."""

    keys = set(reactions) if reactions is not None else set(source) | set(target)
    flipped = []
    for rid in keys:
        a, b = source.get(rid, 0.0), target.get(rid, 0.0)
        if abs(a) > tol and abs(b) > tol and (a > 0) != (b > 0):
            flipped.append(rid)
    return sorted(flipped)
