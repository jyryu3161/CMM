# E. coli succinate production targets

This example executes the existing SC-01 workflow on COBRApy's `textbook` model. Its explicit
condition is glucose uptake of 10 mmol gDW⁻¹ h⁻¹, no oxygen or CO₂ uptake, and the
`glucose_anaerobic` mineral-medium preset. It is a computational demonstration.

Run from the repository root after installing the frozen dependencies and R renderer:

```bash
uv run --frozen --all-extras python examples/production-targets/prepare.py
uv run --frozen --all-extras cmm production-targets --config examples/production-targets/config.json
uv run --frozen --all-extras cmm report validate results/example-production-succinate --json
```

`prepare.py` only exports the model supplied by COBRApy. All numerical stages run through
`cmm production-targets`: MOMA-L2, ROOM, OptKnock, RobustKnock, FSEOF, FVSEOF, loop diagnostics,
flux-response scans, and matched wild-type/knockout sampling. The seeds and candidate-capacity
limits are explicit in [config.json](config.json). Gurobi must have a license large enough for
the model; the restricted license is insufficient for L2 MOMA on `e_coli_core`.

The medium preset contains mineral exchanges that this small model does not represent. The
workflow records those missing components as a warning. Multi-reaction knockouts have explicit
unavailable single-axis response rows; their MOMA/ROOM and paired sampling results remain.

Open `results/example-production-succinate/report_standalone.html` after validation succeeds.
The run retains source CSVs, provenance, 300-DPI PNG and editable PDF/SVG figures. Re-run with
a new `output_dir`; input preparation and the workflow both refuse to overwrite existing work.

For another product or organism, copy this directory, replace the model, and update the exact
exchange and biomass ids, medium and condition before running the same CLI. Model-relative
input paths and the output path resolve from the config's directory. Use the
[customization guide](../../docs/building-custom-workflows.md) for parameter meanings and
the [production skill](../../.agents/skills/cmm-production-engineering/SKILL.md) when asking an
agent to resolve the run.
