"""Check that the same example produced the same scientific result on every OS.

    python .github/scripts/compare_runs.py {sc01,sc02} RUN_DIR_A RUN_DIR_B [RUN_DIR_C ...]

Each RUN_DIR holds that workflow's COMPARED files, copied from one platform's run directory.
The first directory is the reference. Identifiers, ranks, statuses and row order must match exactly;
numbers must agree within REL_TOL/ABS_TOL, which absorbs floating-point and solver-platform
noise without hiding a different optimum or a re-ordered ranking.

Only results with a unique optimum are compared. ROOM, FSEOF/FVSEOF and sampling outputs are
left out: they may legitimately land on different equally optimal or stochastic flux vectors.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import sys

COMPARED = {
    "sc01": (
        "00_summary.json",
        "02_yield/theoretical_yield.csv",
        "02_yield/production_envelope.csv",
        "04_single_knockout/single_knockout_moma.csv",
        "05_strain_design/optknock.csv",
        "05_strain_design/robustknock.csv",
    ),
    "sc02": (
        "00_summary.json",
        "02_reference/source_reference_fluxes.csv",
        "05_transformation/transformation_ranking.csv",
        "06_validation/epsilon_sensitivity.csv",
        "06_validation/moma_baseline.csv",
    ),
}
# Order-only columns left out of the comparison. The single-knockout display order sorts the
# product/growth values exactly, so knockouts that change nothing (deltas of 0 vs -1.8e-15)
# are ordered by rounding noise and differ by platform. The values that define the order are
# still compared within tolerance, so a real change in them still fails.
ORDER_ONLY_COLUMNS = {
    "04_single_knockout/single_knockout_moma.csv": {"display_rank", "method_rank"},
}
REL_TOL = 1e-6
ABS_TOL = 1e-6


def as_number(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def compare(reference: object, other: object, where: str, problems: list[str]) -> float:
    """Append mismatches to ``problems``; return the largest absolute numeric difference."""

    if isinstance(reference, dict) and isinstance(other, dict):
        if reference.keys() != other.keys():
            problems.append(
                f"{where}: keys differ {sorted(reference)} vs {sorted(other)}"
            )
            return 0.0
        return max(
            (
                compare(reference[k], other[k], f"{where}.{k}", problems)
                for k in reference
            ),
            default=0.0,
        )
    if isinstance(reference, list) and isinstance(other, list):
        if len(reference) != len(other):
            problems.append(f"{where}: length {len(reference)} vs {len(other)}")
            return 0.0
        return max(
            (
                compare(a, b, f"{where}[{i}]", problems)
                for i, (a, b) in enumerate(zip(reference, other))
            ),
            default=0.0,
        )
    a, b = as_number(reference), as_number(other)
    if a is not None and b is not None:
        if not math.isclose(a, b, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            problems.append(f"{where}: {reference!r} vs {other!r}")
        return abs(a - b)
    if reference != other:
        problems.append(f"{where}: {reference!r} vs {other!r}")
    return 0.0


def load(path: Path, ignored: set[str]) -> object:
    if path.suffix == ".json":
        return json.loads(path.read_text(encoding="utf-8"))
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            {key: value for key, value in row.items() if key not in ignored}
            for row in csv.DictReader(handle)
        ]


def main(workflow: str, directories: list[Path]) -> None:
    if workflow not in COMPARED:
        raise SystemExit(
            f"workflow must be one of {sorted(COMPARED)}, got {workflow!r}"
        )
    if len(directories) < 2:
        raise SystemExit("need at least two run directories to compare")
    reference, *others = directories
    problems: list[str] = []
    for name in COMPARED[workflow]:
        ignored = ORDER_ONLY_COLUMNS.get(name, set())
        if ignored:
            print(f"{name}: order-only columns not compared: {sorted(ignored)}")
        for other in others:
            largest = compare(
                load(reference / name, ignored),
                load(other / name, ignored),
                f"{other.name}:{name}",
                problems,
            )
            print(
                f"{name}: {reference.name} vs {other.name}, max |difference| = {largest:.3g}"
            )
    if problems:
        raise SystemExit("cross-platform results differ:\n  " + "\n  ".join(problems))
    print(
        f"all {len(COMPARED[workflow])} {workflow} result files agree across "
        f"{len(directories)} platforms"
    )


if __name__ == "__main__":
    main(sys.argv[1], [Path(argument) for argument in sys.argv[2:]])
