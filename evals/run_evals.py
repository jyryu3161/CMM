#!/usr/bin/env python3
"""Check that agent-produced CMM runs honour the contracts in AGENTS.md.

Every check reads artifacts that already exist in a run directory. Nothing here solves a
metabolic model, so a failing check means the run was produced wrongly, not that the science is
wrong. See evals/README.md for what this deliberately does not measure.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parent.parent
TASK_DIR = Path(__file__).resolve().parent / "tasks"


@dataclass
class Check:
    task: str
    name: str
    passed: bool
    detail: str


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _validate(run: Path, kind: str) -> tuple[bool, str]:
    if kind == "production":
        from cmm.reporting import validate_production_run as validator
    else:
        from cmm.reporting.transformation import (
            validate_transformation_run as validator,
        )
    report = validator(run)
    issues = "; ".join(report.issues) or "none"
    return bool(report.valid), f"valid={report.valid} issues={issues}"


def _provenance_fields(run: Path, fields: list[str]) -> tuple[bool, str]:
    provenance = _load_json(run / "00_provenance.json")
    missing = [f for f in fields if provenance.get(f) in (None, "", [], {})]
    return not missing, "all present" if not missing else f"missing {missing}"


def _condition_bounds(run: Path, required: list[str]) -> tuple[bool, str]:
    if not required:
        return True, "no bound requirement declared"
    condition = _load_json(run / "00_config.json").get("condition") or {}
    declared = {b.get("reaction_id") for b in condition.get("bounds", [])}
    missing = [r for r in required if r not in declared]
    return (
        not missing,
        "explicit" if not missing else f"condition does not pin {missing}",
    )


def _tables_retain_status(run: Path, tables: list[str]) -> tuple[bool, str]:
    notes = []
    for rel in tables:
        path = run / rel
        if not path.exists():
            return False, f"missing table {rel}"
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if not rows or "status" not in rows[0]:
            return False, f"{rel} has no status column"
        infeasible = sum(
            1 for r in rows if (r.get("status") or "").strip() == "infeasible"
        )
        notes.append(
            f"{Path(rel).name}: {len(rows)} rows, {infeasible} infeasible retained"
        )
    return True, "; ".join(notes) or "no table requirement declared"


def _planned_not_shipped() -> tuple[bool, str]:
    """A planned feature must never appear as shipped (AGENTS.md rule 6)."""
    from cmm.features import INCLUDED_FEATURES, PLANNED_FEATURES

    leaked = sorted(set(PLANNED_FEATURES) & set(INCLUDED_FEATURES))
    return (
        not leaked,
        "disjoint" if not leaked else f"planned reported as shipped: {leaked}",
    )


def run_task(task: dict[str, Any], run_dir: Path | None) -> list[Check]:
    name = task["id"]
    run = (run_dir or REPO / task["run_dir"]).resolve()
    if not run.exists():
        return [Check(name, "run_directory_exists", False, f"not found: {run}")]

    expect = task.get("expect", {})
    checks = [Check(name, "run_directory_exists", True, str(run))]
    for label, fn in (
        ("validates", lambda: _validate(run, task["kind"])),
        (
            "provenance_fields",
            lambda: _provenance_fields(run, expect.get("provenance_fields", [])),
        ),
        (
            "condition_explicit",
            lambda: _condition_bounds(run, expect.get("condition_bounds_include", [])),
        ),
        (
            "infeasible_retained",
            lambda: _tables_retain_status(run, expect.get("tables_retain_status", [])),
        ),
    ):
        try:
            passed, detail = fn()
        except Exception as exc:  # a check that cannot run is a failure, not a skip
            passed, detail = False, f"{type(exc).__name__}: {exc}"
        checks.append(Check(name, label, passed, detail))
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--task", type=Path, action="append", help="task file (repeatable)"
    )
    parser.add_argument(
        "--run-dir", type=Path, help="override the task's run directory"
    )
    parser.add_argument("--json", dest="json_out", type=Path, help="write results here")
    args = parser.parse_args()

    task_files = args.task or sorted(TASK_DIR.glob("*.json"))
    if args.run_dir and len(task_files) != 1:
        parser.error("--run-dir applies to exactly one --task")

    checks: list[Check] = []
    for path in task_files:
        checks.extend(run_task(_load_json(path), args.run_dir))
    passed, detail = _planned_not_shipped()
    checks.append(Check("repository", "planned_not_shipped", passed, detail))

    for check in checks:
        print(
            f"[{'PASS' if check.passed else 'FAIL'}] {check.task}/{check.name}: {check.detail}"
        )
    n_pass = sum(1 for c in checks if c.passed)
    print(f"\n{n_pass}/{len(checks)} checks passed")

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "checks": [asdict(c) for c in checks],
                    "n_passed": n_pass,
                    "n_total": len(checks),
                    "pass_rate": round(n_pass / len(checks), 4) if checks else 0.0,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    return 0 if n_pass == len(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
