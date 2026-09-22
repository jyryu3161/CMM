#!/usr/bin/env python3
"""Flag analysis code that bypasses CMM's services (AGENTS.md rule 12).

An agent asked for an analysis CMM does not ship is tempted to write the method inline. The
numbers that come back look ordinary but carry no `run_provenance`, no method contract in
`docs/VALIDATION.md`, and no test, so nobody can check them later.

This is a static reader, not a sandbox: it inspects scripts an agent wrote and reports where
they solve, sample, or mutate a model outside CMM. It cannot prevent anything — an agent that
never writes a file is invisible to it — so treat a clean report as "no bypass found in the
files given", never as proof that CMM produced every number.

Point it at agent-authored scripts only. CMM's own source legitimately calls the solver.

Usage:
    python evals/check_service_use.py analysis.py [more.py ...]
    python evals/check_service_use.py --json report.json analysis.py
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

# Running a solve. Whatever it returns did not pass through a CMM service.
SOLVE_METHODS = frozenset({"optimize", "slim_optimize"})

# cobra's own analysis entry points, each with the CMM function that owns it (AGENTS.md §1).
COBRA_ANALYSIS = {
    "flux_variability_analysis": "cmm.core.fva",
    "pfba": "cmm.core.pfba",
    "moma": "cmm.features.knockout_comparison",
    "room": "cmm.features.knockout_comparison",
    "sample": "cmm.features.random_flux_sampling",
    "production_envelope": "cmm.features.production_envelope",
    "single_gene_deletion": "cmm.features.batch_comparison",
    "single_reaction_deletion": "cmm.features.batch_comparison",
    "double_gene_deletion": "cmm.features.batch_comparison",
    "double_reaction_deletion": "cmm.features.batch_comparison",
    "find_essential_genes": "cmm.features.batch_comparison",
    "find_essential_reactions": "cmm.features.batch_comparison",
}

# Driving a solver directly leaves CMM out of the loop entirely.
SOLVER_MODULES = frozenset({"gurobipy", "cplex", "optlang", "straindesign"})

# Mutating model state by hand loses the condition record that provenance is built from.
# `cmm.core` owns conditions and media; these are reported for review, not as a bypass.
STATE_ATTRS = frozenset(
    {
        "objective",
        "objective_direction",
        "medium",
        "bounds",
        "lower_bound",
        "upper_bound",
    }
)


@dataclass
class Finding:
    path: str
    line: int
    severity: str  # "bypass" (ran the method itself) | "review" (changed state by hand)
    symbol: str
    message: str


class _Visitor(ast.NodeVisitor):
    def __init__(self, path: str) -> None:
        self.path = path
        self.findings: list[Finding] = []
        self._cobra_analysis_names: dict[str, str] = {}

    def _add(self, node: ast.AST, severity: str, symbol: str, message: str) -> None:
        self.findings.append(
            Finding(self.path, getattr(node, "lineno", 0), severity, symbol, message)
        )

    def visit_Import(self, node: ast.Import) -> None:
        for alias in node.names:
            root = alias.name.split(".")[0]
            if root in SOLVER_MODULES:
                self._add(
                    node,
                    "bypass",
                    alias.name,
                    f"imports the solver package '{alias.name}' directly; CMM owns solver "
                    "selection and capability checks (AGENTS.md §2)",
                )
        self.generic_visit(node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> None:
        module = node.module or ""
        root = module.split(".")[0]
        if root in SOLVER_MODULES:
            self._add(
                node,
                "bypass",
                module,
                f"imports from the solver package '{module}' directly; CMM owns solver "
                "selection and capability checks (AGENTS.md §2)",
            )
        # `from cobra.flux_analysis import pfba` binds a bare name we must follow.
        elif module.startswith("cobra"):
            for alias in node.names:
                if alias.name in COBRA_ANALYSIS:
                    self._cobra_analysis_names[alias.asname or alias.name] = alias.name
        self.generic_visit(node)

    def visit_Call(self, node: ast.Call) -> None:
        func = node.func
        if isinstance(func, ast.Attribute):
            if func.attr in SOLVE_METHODS:
                self._add(
                    node,
                    "bypass",
                    f".{func.attr}()",
                    "solves the model directly; use a CMM service so the result carries "
                    "run_provenance (cmm.core.fba / pfba / fva)",
                )
            elif func.attr in COBRA_ANALYSIS and _mentions_cobra(func.value):
                self._add(
                    node,
                    "bypass",
                    f"cobra…{func.attr}()",
                    f"runs cobra's own {func.attr}; CMM ships this as "
                    f"{COBRA_ANALYSIS[func.attr]}",
                )
        elif isinstance(func, ast.Name) and func.id in self._cobra_analysis_names:
            original = self._cobra_analysis_names[func.id]
            self._add(
                node,
                "bypass",
                f"{func.id}()",
                f"runs cobra's own {original}; CMM ships this as {COBRA_ANALYSIS[original]}",
            )
        self.generic_visit(node)

    def visit_Assign(self, node: ast.Assign) -> None:
        for target in node.targets:
            self._check_state_target(node, target)
        self.generic_visit(node)

    def visit_AugAssign(self, node: ast.AugAssign) -> None:
        self._check_state_target(node, node.target)
        self.generic_visit(node)

    def _check_state_target(self, node: ast.AST, target: ast.AST) -> None:
        if isinstance(target, ast.Attribute) and target.attr in STATE_ATTRS:
            self._add(
                node,
                "review",
                f".{target.attr}",
                f"sets '{target.attr}' by hand; cmm.core owns conditions and media, and a "
                "bound changed here is absent from the run's condition record",
            )


def _mentions_cobra(node: ast.AST) -> bool:
    """True when an attribute chain is rooted at something named cobra."""
    while isinstance(node, ast.Attribute):
        node = node.value
    return isinstance(node, ast.Name) and node.id == "cobra"


def check_file(path: Path) -> list[Finding]:
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as exc:
        return [
            Finding(
                str(path), exc.lineno or 0, "bypass", "syntax", f"cannot parse: {exc}"
            )
        ]
    visitor = _Visitor(str(path))
    visitor.visit(tree)
    return visitor.findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "paths", type=Path, nargs="+", help="agent-authored scripts to read"
    )
    parser.add_argument(
        "--json", dest="json_out", type=Path, help="write findings here"
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="fail on 'review' findings too, not only on a bypass",
    )
    args = parser.parse_args()

    findings: list[Finding] = []
    for path in args.paths:
        if not path.exists():
            print(f"error: no such file: {path}", file=sys.stderr)
            return 2
        findings.extend(check_file(path))

    for f in sorted(findings, key=lambda f: (f.path, f.line)):
        print(f"[{f.severity}] {f.path}:{f.line} {f.symbol} — {f.message}")

    bypasses = [f for f in findings if f.severity == "bypass"]
    reviews = [f for f in findings if f.severity == "review"]
    scanned = len(args.paths)
    if not findings:
        print(
            f"no bypass found in {scanned} file(s) — every analysis call goes through CMM"
        )
    else:
        print(
            f"\n{len(bypasses)} bypass, {len(reviews)} to review, across {scanned} file(s)"
        )

    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(
                {
                    "findings": [asdict(f) for f in findings],
                    "n_bypass": len(bypasses),
                    "n_review": len(reviews),
                    "n_files": scanned,
                },
                indent=2,
            ),
            encoding="utf-8",
        )
    return 1 if bypasses or (args.strict and reviews) else 0


if __name__ == "__main__":
    sys.exit(main())
