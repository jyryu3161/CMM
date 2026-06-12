<#
.SYNOPSIS
  CMM source installer for Windows PowerShell. Clone the repo, then run .\install.ps1

.DESCRIPTION
  Creates an isolated virtual environment and installs the workbench from this checkout
  (editable), so `python -m cmm.app` and the `cmm` command work immediately. The macOS/Linux
  equivalent is install.sh.

.EXAMPLE
  .\install.ps1                       # core + desktop GUI + strain design + gurobipy, into .\.venv
  .\install.ps1 -Dev                  # also install the test/lint tooling (pytest, ruff)
  .\install.ps1 -NoGurobi            # skip gurobipy (open GLPK runs LP only)
  .\install.ps1 -CoreOnly            # core library only (no GUI / strain-design extras)
  .\install.ps1 -Python python3.12   # use a specific interpreter
  .\install.ps1 -VenvDir C:\envs\cmm # install into a chosen venv directory (default: .\.venv)
#>
[CmdletBinding()]
param(
  [switch]$Dev,
  [switch]$NoGurobi,
  [switch]$CoreOnly,
  [string]$Python = "python",
  [string]$VenvDir = ".venv"
)

$ErrorActionPreference = "Stop"

# Run from the repository root (the directory holding this script) regardless of CWD.
Set-Location -Path $PSScriptRoot

$extras = if ($CoreOnly) { "" } else { "desktop,design" }
if ($Dev -and -not $CoreOnly) { $extras = "$extras,dev" }
if ($Dev -and $CoreOnly)      { $extras = "dev" }

# --- locate a Python >= 3.10 ----------------------------------------------------------
if (-not (Get-Command $Python -ErrorAction SilentlyContinue)) {
  Write-Error "no Python interpreter '$Python' found. Install Python >= 3.10 first."
}
& $Python -c "import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)"
if ($LASTEXITCODE -ne 0) {
  Write-Error "CMM needs Python >= 3.10, but '$Python' is older."
}
Write-Host "==> Using $(& $Python -V) ($Python)"

# --- create / reuse the virtual environment -------------------------------------------
if (-not (Test-Path $VenvDir)) {
  Write-Host "==> Creating virtual environment at $VenvDir"
  & $Python -m venv $VenvDir
} else {
  Write-Host "==> Reusing existing virtual environment at $VenvDir"
}
$vpy = Join-Path $VenvDir "Scripts\python.exe"
if (-not (Test-Path $vpy)) { $vpy = Join-Path $VenvDir "bin/python" }  # tolerate cross-shell venvs

# --- install --------------------------------------------------------------------------
Write-Host "==> Upgrading pip"
& $vpy -m pip install --upgrade pip | Out-Null

if ($extras -ne "") {
  Write-Host "==> Installing CMM (editable) with extras: $extras"
  & $vpy -m pip install -e ".[$extras]"
} else {
  Write-Host "==> Installing CMM (editable, core only)"
  & $vpy -m pip install -e .
}

if (-not $NoGurobi) {
  Write-Host "==> Installing gurobipy (free size-limited license; enables MOMA/ROOM/OptKnock/MTA on small models)"
  & $vpy -m pip install gurobipy
} else {
  Write-Host "==> Skipping gurobipy (open GLPK runs FBA/pFBA/FVA/FSEOF; QP/MILP/MIQP features need Gurobi or CPLEX)"
}

# --- done -----------------------------------------------------------------------------
$solver = (& $vpy -c "from cmm.core.solvers import solver_status; s = solver_status(); print(s.name, '(' + ', '.join(s.capabilities) + ')')" 2>$null)
if (-not $solver) { $solver = "unknown" }
Write-Host ""
Write-Host "==> CMM installed. Active solver: $solver"
Write-Host ""
Write-Host "   Launch the desktop workbench:"
Write-Host "     $vpy -m cmm.app"
Write-Host "   Or activate the environment first:"
Write-Host "     $VenvDir\Scripts\Activate.ps1"
Write-Host "     python -m cmm.app"
