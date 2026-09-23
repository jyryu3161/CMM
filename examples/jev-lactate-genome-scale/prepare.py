#!/usr/bin/env python3
"""Write the iJO1366 input for the genome-scale D-lactate example.

COBRApy ships iJO1366, so this unpacks it rather than committing a second copy — the same
arrangement the SC-01 and SC-02 examples use. It refuses to overwrite an existing input so a
curated study model is never silently replaced.
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path

import cobra

DEST = Path(__file__).parent / "data" / "iJO1366.xml"


def main() -> int:
    if DEST.exists():
        print(f"Input already exists; preserve or move it first: {DEST}")
        return 1
    source = Path(cobra.__file__).parent / "data" / "iJO1366.xml.gz"
    if not source.is_file():  # pragma: no cover - depends on the cobra install
        print(f"COBRApy does not ship {source.name} in this install")
        return 1
    DEST.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rb") as packed, open(DEST, "wb") as out:
        shutil.copyfileobj(packed, out)
    print(f"wrote {DEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
