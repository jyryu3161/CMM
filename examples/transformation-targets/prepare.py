"""Export CMM's existing synthetic source/target fixture; run analyses through the CLI."""

from pathlib import Path

import pandas as pd
from cobra.io import write_sbml_model

from cmm.app.screenshots import (
    SOURCE_EXPRESSION,
    TARGET_EXPRESSION,
    build_demo_model,
)


def main() -> None:
    directory = Path(__file__).resolve().parent / "data"
    if directory.exists():
        raise SystemExit(
            f"Input directory already exists; preserve or move it: {directory}"
        )
    directory.mkdir(parents=True)
    write_sbml_model(build_demo_model(), str(directory / "demo_disease.xml"))
    for name, values in (
        ("source_disease", SOURCE_EXPRESSION),
        ("target_healthy", TARGET_EXPRESSION),
    ):
        pd.DataFrame.from_dict(values, orient="index", columns=["measurement"]).to_csv(
            directory / f"{name}.csv", index_label="gene"
        )
    print(directory)


if __name__ == "__main__":
    main()
