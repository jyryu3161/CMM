#!/usr/bin/env bash
#
# CMM source installer — clone the repo, then run ./install.sh
#
# Creates an isolated virtual environment and installs the workbench from this checkout
# (editable), so `python -m cmm.app` and the `cmm` command work immediately. Works on
# macOS, Linux, and Windows (Git Bash / WSL); for native PowerShell use install.ps1.
#
# Usage:
#   ./install.sh                      # core + desktop GUI + strain design + gurobipy, into ./.venv
#   ./install.sh --dev                # also install the test/lint tooling (pytest, ruff)
#   ./install.sh --no-gurobi          # skip gurobipy (open GLPK runs LP only: FBA/pFBA/FVA/FSEOF)
#   ./install.sh --core-only          # core library only (no GUI / strain-design extras)
#   ./install.sh --python python3.12  # use a specific interpreter
#   ./install.sh --venv /path/to/venv # install into a chosen venv directory (default: ./.venv)
#
set -euo pipefail

# Run from the repository root (the directory holding this script) regardless of CWD.
cd "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

EXTRAS="desktop,design"
WITH_GUROBI=1
PYTHON=""
VENV_DIR=".venv"

while [ $# -gt 0 ]; do
  case "$1" in
    --dev)        EXTRAS="${EXTRAS:+$EXTRAS,}dev" ;;
    --no-gurobi)  WITH_GUROBI=0 ;;
    --core-only)  EXTRAS="" ;;
    --python)     PYTHON="${2:?--python needs an interpreter}"; shift ;;
    --venv)       VENV_DIR="${2:?--venv needs a path}"; shift ;;
    -h|--help)    sed -n '2,20p' "$0"; exit 0 ;;
    *)            echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

# --- locate a Python >= 3.10 ------------------------------------------------------------
if [ -z "$PYTHON" ]; then
  for cand in python3 python; do
    if command -v "$cand" >/dev/null 2>&1; then PYTHON="$cand"; break; fi
  done
fi
if [ -z "$PYTHON" ] || ! command -v "$PYTHON" >/dev/null 2>&1; then
  echo "error: no Python interpreter found. Install Python >= 3.10 first." >&2
  exit 1
fi
if ! "$PYTHON" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)'; then
  echo "error: CMM needs Python >= 3.10, but '$PYTHON' is $("$PYTHON" -V 2>&1)." >&2
  exit 1
fi
echo "==> Using $("$PYTHON" -V 2>&1) ($PYTHON)"

# --- create / reuse the virtual environment --------------------------------------------
if [ ! -d "$VENV_DIR" ]; then
  echo "==> Creating virtual environment at $VENV_DIR"
  "$PYTHON" -m venv "$VENV_DIR"
else
  echo "==> Reusing existing virtual environment at $VENV_DIR"
fi
# venv layout is bin/ on Unix, Scripts/ on Windows (Git Bash).
if [ -d "$VENV_DIR/bin" ]; then VPY="$VENV_DIR/bin/python"; else VPY="$VENV_DIR/Scripts/python.exe"; fi

# --- install ---------------------------------------------------------------------------
echo "==> Upgrading pip"
"$VPY" -m pip install --upgrade pip >/dev/null

if [ -n "$EXTRAS" ]; then
  echo "==> Installing CMM (editable) with extras: $EXTRAS"
  "$VPY" -m pip install -e ".[$EXTRAS]"
else
  echo "==> Installing CMM (editable, core only)"
  "$VPY" -m pip install -e .
fi

if [ "$WITH_GUROBI" -eq 1 ]; then
  echo "==> Installing gurobipy (free size-limited license; enables MOMA/ROOM/OptKnock/MTA on small models)"
  "$VPY" -m pip install gurobipy
else
  echo "==> Skipping gurobipy (open GLPK runs FBA/pFBA/FVA/FSEOF; QP/MILP/MIQP features need Gurobi or CPLEX)"
fi

# --- done ------------------------------------------------------------------------------
SOLVER=$("$VPY" -c 'from cmm.core.solvers import solver_status; s = solver_status(); print(s.name, "(" + ", ".join(s.capabilities) + ")")' 2>/dev/null || echo "unknown")
cat <<EOF

==> CMM installed. Active solver: $SOLVER

   Launch the desktop workbench:
     $VPY -m cmm.app
   Or activate the environment first:
     source $VENV_DIR/bin/activate   # (Windows Git Bash: source $VENV_DIR/Scripts/activate)
     python -m cmm.app
EOF
