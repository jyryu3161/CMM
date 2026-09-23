#!/usr/bin/env python3
"""Repeat one JEV configuration n times and report the distribution of what it found.

A single JEV run is evidence about that run. The CMM solves in it repeat exactly; the decision
model's choices are not guaranteed to, so a design quoted from one run says nothing about the
method. This runs the same config n times and reports the spread — which is the minimum any
claim about the agent has to rest on.

It reports the **guaranteed product**, because that is what the run ranks designs on. Two runs
whose pFBA numbers look alike can differ completely in what the strain must actually make.

Usage:
    python evals/jev_replicates.py CONFIG --runs 10
    python evals/jev_replicates.py CONFIG --runs 10 --target 17.5858 --json out.json

``--target`` is what an exhaustive search of the same vocabulary reached, when that is known;
the report then also counts how often the agent got there. Without it the spread is still
reported, just with nothing to be a fraction of.

Needs ``OPENROUTER_API_KEY``, and costs whatever n runs cost — about $0.015 a run on
``e_coli_core`` and $0.025 on ``iJO1366`` as measured. Run directories go under ``--out``,
which defaults to a temporary directory, because n run bundles are not something to leave in a
repository by accident.
"""

from __future__ import annotations

import argparse
import json
import statistics
import tempfile
import time
from pathlib import Path
from typing import Any

# Settings that describe the run rather than where it is written. Everything here is copied
# from the config; ``model_path``, ``condition``, ``jev_model`` and ``output_dir`` are handled
# separately because they are paths, objects, or per-replicate.
_CARRIED = frozenset(
    {
        "product",
        "biomass",
        "organism",
        "brief",
        "rounds",
        "steps_per_round",
        "max_knockouts",
        "max_knockdowns",
        "growth_floor",
        "candidate_limit",
        "allow_look_actions",
        "design_max_knockouts",
        "design_max_solutions",
        "seed_with_strain_design",
        "require_distinct_rounds",
        "run_baseline_comparison",
        "measure_cofactor_limits",
        "measure_guaranteed_product",
        "screen_interventions",
        "run_moma",
        "question_set",
        "enable_web_research",
        "max_decisions",
        "max_cost_usd",
        "seed",
    }
)


def _one(base, out: Path) -> dict[str, Any]:
    from cmm.jev import JevConfig, run_jev_design

    provenance = base.to_provenance()
    config = JevConfig.from_mapping(
        {
            **{k: v for k, v in provenance.items() if k in _CARRIED},
            "model_path": str(base.model_path),
            "condition": base.condition,
            "jev_model": base.jev_model,
            "output_dir": str(out),
            "overwrite": True,
        }
    )
    started = time.perf_counter()
    result = run_jev_design(config)
    return {
        "guaranteed_product": result.best_guaranteed_product,
        "product_flux": result.best_product_flux,
        "growth": result.best_growth,
        "design": [i.describe() for i in result.best_interventions],
        "n_edits": len(result.best_interventions),
        "n_ticks": len(result.ticks),
        "cost_usd": float(result.usage.get("cost_usd", 0.0)),
        "seconds": round(time.perf_counter() - started, 1),
        "run_directory": str(result.run_directory) if result.run_directory else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("config", type=Path, help="a JevConfig JSON file")
    parser.add_argument("--runs", type=int, default=10)
    parser.add_argument(
        "--target",
        type=float,
        default=None,
        help="what an exhaustive search of the same vocabulary reached, if known",
    )
    parser.add_argument("--out", type=Path, default=None, help="where run bundles go")
    parser.add_argument("--json", dest="json_out", type=Path, default=None)
    args = parser.parse_args()

    from cmm.jev import JevConfig

    base = JevConfig.from_json(args.config)
    root = args.out or Path(tempfile.mkdtemp(prefix="jev-replicates-"))
    root.mkdir(parents=True, exist_ok=True)

    runs: list[dict[str, Any]] = []
    for index in range(1, args.runs + 1):
        try:
            record = _one(base, root / f"run_{index:02d}")
        except (
            Exception
        ) as error:  # a failed replicate is a data point, not a lost study
            record = {"error": f"{type(error).__name__}: {error}"}
        runs.append(record)
        value = record.get("guaranteed_product")
        shown = f"{value:10.4f}" if isinstance(value, (int, float)) else str(value)
        print(
            f"run {index:3d}: guaranteed {shown}  edits {record.get('n_edits', '-')}"
            f"  ticks {record.get('n_ticks', '-')}  ${record.get('cost_usd', 0.0):.4f}"
            f"  {record.get('seconds', '-')}s",
            flush=True,
        )

    scored = [
        r["guaranteed_product"]
        for r in runs
        if isinstance(r.get("guaranteed_product"), (int, float))
    ]
    report: dict[str, Any] = {
        "config": str(args.config),
        "product": base.product,
        "n_requested": args.runs,
        "n_scored": len(scored),
        "n_failed": sum(1 for r in runs if "error" in r),
        "median_guaranteed_product": statistics.median(scored) if scored else None,
        "min_guaranteed_product": min(scored) if scored else None,
        "max_guaranteed_product": max(scored) if scored else None,
        "n_returned_nothing": sum(1 for v in scored if abs(v) <= 1e-9),
        "total_cost_usd": round(sum(r.get("cost_usd", 0.0) for r in runs), 4),
        "target": args.target,
        "n_reached_target": (
            sum(1 for v in scored if v >= args.target - 1e-3)
            if args.target is not None
            else None
        ),
        "distinct_designs": len({tuple(r.get("design", ())) for r in runs}),
        "runs": runs,
        "note": (
            "The guaranteed product is the quantity the run ranks designs on. A spread here is "
            "the method's spread; a single run is not the method's performance."
        ),
    }
    print(
        f"\nmedian {report['median_guaranteed_product']}, "
        f"range [{report['min_guaranteed_product']}, {report['max_guaranteed_product']}], "
        f"{report['distinct_designs']} distinct design(s), "
        f"{report['n_returned_nothing']}/{len(scored)} returned nothing"
        + (
            f", {report['n_reached_target']}/{len(scored)} reached {args.target}"
            if args.target is not None
            else ""
        )
        + f", ${report['total_cost_usd']} total"
    )
    if args.json_out:
        args.json_out.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"written to {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
