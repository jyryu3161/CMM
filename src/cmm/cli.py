"""Command line entry point for CMM."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
import json
from pathlib import Path
import sys

from cmm import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cmm")
    parser.add_argument("--version", action="store_true", help="Print the CMM version.")
    commands = parser.add_subparsers(dest="command")

    production = commands.add_parser(
        "production-targets",
        help="Run the canonical production-target-discovery workflow from JSON config.",
    )
    production.add_argument(
        "--config",
        required=True,
        type=Path,
        help="UTF-8 JSON ProductionWorkflowConfig file.",
    )
    production.add_argument(
        "--analysis-only",
        action="store_true",
        help="Write scientific artifacts without invoking the R publication renderer.",
    )
    production.add_argument(
        "--renderer",
        default="nature-r",
        choices=("nature-r",),
        help="Publication renderer used after a successful analysis.",
    )

    transformation = commands.add_parser(
        "transformation-targets",
        help="Rank knockouts that move a source metabolic state toward a target state.",
    )
    transformation.add_argument(
        "--config",
        required=True,
        type=Path,
        help="UTF-8 JSON TransformationWorkflowConfig file.",
    )
    transformation.add_argument(
        "--analysis-only",
        action="store_true",
        help="Write scientific artifacts without rendering figures or the HTML report.",
    )
    transformation.add_argument(
        "--highlight",
        default=None,
        help="Candidate to mark throughout the report, for example the knockout under test.",
    )

    jev = commands.add_parser(
        "jev-design",
        help="Let the JEV decision model play the model to raise a product flux.",
    )
    jev.add_argument(
        "--config",
        required=True,
        type=Path,
        help="UTF-8 JSON JevConfig file.",
    )
    jev.add_argument(
        "--quiet",
        action="store_true",
        help="Do not print each move as it is played.",
    )

    report = commands.add_parser("report", help="Render or validate a schema-v2 run.")
    report_commands = report.add_subparsers(dest="report_command", required=True)
    render = report_commands.add_parser(
        "render", help="Render R figures and linked/standalone HTML reports."
    )
    render.add_argument("run_dir", type=Path)
    render.add_argument("--renderer", default="nature-r", choices=("nature-r",))
    validate = report_commands.add_parser(
        "validate", help="Validate the run manifest and publication source artifacts."
    )
    validate.add_argument("run_dir", type=Path)
    validate.add_argument(
        "--json",
        action="store_true",
        dest="as_json",
        help="Print the validation result as JSON.",
    )
    return parser


def _validation_payload(report) -> dict[str, object]:
    return {
        "valid": report.valid,
        "issues": list(report.issues),
        "warnings": list(report.warnings),
        "run_directory": str(report.run.root) if report.run is not None else None,
    }


def _run_production(args: argparse.Namespace) -> int:
    from cmm.reporting import render_production_report, validate_production_run
    from cmm.workflows import ProductionWorkflowConfig, run_production_target_discovery

    config = ProductionWorkflowConfig.from_json(args.config)
    if config.output_dir is None:
        raise ValueError(
            "production-targets requires config.output_dir so the analysis has a "
            "self-contained run directory"
        )
    result = run_production_target_discovery(config)
    if result.run_directory is None:  # guarded above, retained as an invariant check
        raise RuntimeError(
            "production workflow completed without exporting a run directory"
        )
    if args.analysis_only:
        print(result.run_directory)
        return 0

    bundle = render_production_report(result.run_directory, renderer=args.renderer)
    validation = validate_production_run(result.run_directory)
    validation.raise_for_errors()
    print(
        json.dumps(
            {
                "run_directory": str(result.run_directory),
                "report_html": str(bundle.report.report_html),
                "report_standalone_html": str(bundle.report.report_standalone_html),
                "figure_manifest": str(bundle.figures.path),
                "valid": validation.valid,
            },
            indent=2,
        )
    )
    return 0


def _run_transformation(args: argparse.Namespace) -> int:
    from cmm.workflows.transformation import (
        TransformationWorkflowConfig,
        run_transformation_target_discovery,
    )

    config = TransformationWorkflowConfig.from_json(args.config)
    if config.output_dir is None:
        raise ValueError(
            "transformation-targets requires config.output_dir so the analysis has a "
            "self-contained run directory"
        )
    result = run_transformation_target_discovery(config)
    if result.run_directory is None:  # guarded above, retained as an invariant check
        raise RuntimeError(
            "transformation workflow completed without exporting a run directory"
        )
    summary = result.summary()
    payload: dict[str, object] = {"run_directory": str(result.run_directory)}
    if not args.analysis_only:
        from cmm.reporting import (
            render_transformation_report,
            validate_transformation_run,
        )

        report = render_transformation_report(
            result.run_directory, highlight=args.highlight
        )
        payload["report_html"] = str(report.report_html)
        # The copy to send someone: the linked page loses every figure once it is moved.
        payload["report_standalone_html"] = str(report.report_standalone_html)
        payload["figures"] = [str(path) for path in report.figures]
        # A rendered page is not a finished run. Same gate the CLI's report subcommand applies.
        validation = validate_transformation_run(result.run_directory)
        payload["validation"] = _validation_payload(validation)
    print(
        json.dumps(
            {
                **payload,
                "method": summary["method"],
                "n_candidates": summary["n_candidates"],
                "top_target": summary["top_target"],
                "top_score": summary["top_score"],
                # Stated on every run: the reference state is not the published iMAT one.
                "reference_method": summary["reference_method"],
            },
            indent=2,
        )
    )
    return 0


def _run_jev(args: argparse.Namespace) -> int:
    from cmm.jev import JevConfig, run_jev_design

    config = JevConfig.from_json(args.config)
    if config.output_dir is None:
        raise ValueError(
            "jev-design requires config.output_dir so the run has a self-contained "
            "directory; the agent transcript is the only record of why each move was made"
        )

    # Moves are printed as they happen: a run is a sequence of decisions, and watching it is
    # most of the point. --quiet is there for scripting.
    def announce(tick, _fluxes) -> None:
        print(tick.headline(), flush=True)

    result = run_jev_design(config, on_tick=None if args.quiet else announce)
    summary = result.summary()
    if not args.quiet and result.baselines:
        # What the established methods give on the same problem, printed next to the run
        # rather than left in a CSV, because a design with nothing to compare it to is not a
        # result a reader can weigh.
        print("\nSame problem, same growth floor, every design scored the same way:")
        frame = result.baselines_frame()
        # Width from the data, not a guess: the amplification-headroom row's label is 60
        # characters and a fixed column wrapped it onto the numbers.
        width = max(len(str(row["method"])) for _, row in frame.iterrows())
        for _, row in frame.iterrows():
            flag = "" if row["deterministic"] else "  (not deterministic)"
            print(
                f"  {str(row['method']):{width}s}  product {row['product_flux']:9.4f}  "
                f"growth {row['growth']:7.4f}{flag}"
            )
        print()
    if not args.quiet:
        # The design is the headline; the per-target evidence is the rest of what the run
        # learned, and a reader weighing a different trade-off needs it.
        rounds = result.rounds_frame()
        if not rounds.empty and "shortfall" in rounds:
            print("What each round asked, reached, and left undone:")
            for _, row in rounds.iterrows():
                gap = str(row["shortfall"] or "nothing measurable")
                print(
                    f"  round {int(row['round'])} \u2014 {row['question_answered']}: "
                    f"{row['product_flux']:.4g} at growth {row['growth']:.4g}; "
                    f"{row['stopped_because']}; {gap}"
                )
            distinct = rounds["design_signature"].nunique()
            print(
                f"  {distinct} distinct design(s) over {len(rounds)} round(s). The headline "
                "is the best of them; the rest are what a laboratory that cannot build it "
                "would use."
            )
            print()
    print(
        json.dumps(
            {
                "run_directory": str(result.run_directory),
                "product": summary["product"],
                "wild_type_product_flux": summary["wild_type_product_flux"],
                "best_product_flux": summary["best_product_flux"],
                "best_growth": summary["best_growth"],
                "beat_wild_type": summary["beat_wild_type"],
                "best_design": summary["best_design"],
                "n_ticks": summary["n_ticks"],
                "usage": summary["usage"],
                "baseline_comparison": summary["baseline_comparison"],
                "targets": summary["targets"],
                # Stated on every run: the CMM solves repeat, the agent's choices need not.
                "notes": summary["notes"],
            },
            indent=2,
        )
    )
    return 0


def _workflow_of(run_dir: str | Path) -> str:
    """Read the run's own workflow id, so the caller never has to name it.

    A run directory already says what it is. Asking for a flag would make it possible to
    validate a transformation run against the production gate, which reports a long list of
    missing production artifacts instead of the one useful fact.
    """

    manifest = Path(run_dir).expanduser().resolve() / "00_manifest.json"
    if not manifest.is_file():
        return "unknown"
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return "unknown"
    return (
        str(payload.get("workflow", "unknown"))
        if isinstance(payload, dict)
        else "unknown"
    )


def _run_report(args: argparse.Namespace) -> int:
    from cmm.reporting import (
        render_production_report,
        render_transformation_report,
        validate_production_run,
        validate_transformation_run,
    )

    workflow = _workflow_of(args.run_dir)
    if workflow == "jev_target_design":
        # A JEV run writes a schema-v2 bundle but has no publication renderer or completion
        # gate of its own. Falling through would validate it against the production contract
        # and report a long list of artifacts it was never meant to contain.
        print(
            "this is an agent run; it has no publication renderer or completion gate. "
            "Read report.html for the whole run on one page; the tables behind it are in "
            "02_game/, 03_design/, 05_baseline/ and 06_targets/, the design-space figure "
            "is in figures/, and every request and response is in "
            "04_agent/transcript.jsonl.",
            file=sys.stderr,
        )
        return 1
    transformation = workflow == "transformation_target_discovery"

    if args.report_command == "render":
        if transformation:
            report = render_transformation_report(args.run_dir)
            print(report.report_html)
            print(report.report_standalone_html)
            return 0
        bundle = render_production_report(args.run_dir, renderer=args.renderer)
        print(bundle.report.report_html)
        print(bundle.report.report_standalone_html)
        return 0

    validation = (
        validate_transformation_run(args.run_dir)
        if transformation
        else validate_production_run(args.run_dir)
    )
    payload = _validation_payload(validation)
    if args.as_json:
        print(json.dumps(payload, indent=2))
    elif validation.valid:
        kind = "transformation" if transformation else "production"
        print(f"valid schema-v2 {kind} run: {Path(args.run_dir).resolve()}")
        for warning in validation.warnings:
            print(f"warning: {warning}")
    else:
        for issue in validation.issues:
            print(f"error: {issue}", file=sys.stderr)
    return 0 if validation.valid else 1


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.version:
        print(__version__)
        return 0

    try:
        if args.command == "production-targets":
            return _run_production(args)
        if args.command == "transformation-targets":
            return _run_transformation(args)
        if args.command == "jev-design":
            return _run_jev(args)
        if args.command == "report":
            return _run_report(args)
    except Exception as error:
        print(f"cmm: {type(error).__name__}: {error}", file=sys.stderr)
        return 1

    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
