"""Save COBRApy's existing textbook model for the canonical succinate example."""

from pathlib import Path

from cobra.io import load_model, write_sbml_model


def main() -> None:
    destination = Path(__file__).resolve().parent / "data" / "e_coli_core.xml"
    if destination.exists():
        raise SystemExit(
            f"Input already exists; preserve or move it first: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_sbml_model(load_model("textbook"), str(destination))
    print(destination)


if __name__ == "__main__":
    main()
