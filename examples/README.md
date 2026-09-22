# Canonical workflow examples

Run these from a source checkout with the frozen Python dependencies and the R renderer
installed. Both use the public CLI and write a validated, self-contained report directory.

| Workflow | Runnable example | Output |
|---|---|---|
| SC-01 production-target discovery | [E. coli succinate](production-targets/README.md) | `results/example-production-succinate/` |
| SC-02 transformation-target discovery | [MTA and rMTA on the existing synthetic model](transformation-targets/README.md) | `results/example-transformation-mta/`, `results/example-transformation-rmta/` |
| SC-03 agent design *(not a canonical workflow)* | [JEV on anaerobic succinate](jev-design/README.md) | `results/example-jev-succinate/` |

SC-03 is listed for completeness and is **not** a canonical workflow: it has no publication
renderer and no completion gate, `cmm report` refuses it, and it needs `OPENROUTER_API_KEY`.
Its decisions do not repeat, so one run is never the method's performance.

The small SC-02 model demonstrates execution and reporting. The historical Recon1 runtime
mentioned in the main README does not include a reproducible input bundle in this repository.
See [adapting SC-02 to Recon](transformation-targets/README.md#using-recon-or-another-model)
for the required study inputs; the synthetic example does not validate Recon predictions.

For a new product, model, expression pair or condition, copy the corresponding example and
change its config. For a new scientific workflow with different stages or artifacts, follow
the [customization guide](../docs/building-custom-workflows.md) and
[contributor tutorial](../docs/tutorials/adding-a-canonical-workflow.md).

## Non-obvious patterns

- **Gotcha: `prepare.py` refuses to overwrite an existing input.** It exits with
  `Input already exists; preserve or move it first` so a curated study input is never silently
  replaced. To re-run an example from scratch, move or delete that example's `data/` directory
  first:

  ```bash
  rm -rf examples/production-targets/data      # or move it somewhere safe
  uv run --frozen --all-extras python examples/production-targets/prepare.py
  ```

- **Note: a run directory is not overwritten either.** Each config names an explicit
  `output_dir`; remove the previous run or point the config elsewhere before repeating it.
- **Why: the example inputs are generated, not committed.** `prepare.py` writes them from
  COBRApy's bundled models or CMM's own fixtures, and `MANIFEST.in` prunes both `data/`
  directories, so the distribution stays free of duplicated model files.

## See also

- [`docs/scenarios/`](../docs/scenarios/README.md) — what each workflow's numbers mean
- [`evals/`](../evals/README.md) — contract checks these example runs are expected to pass
- [`AGENTS.md`](../AGENTS.md) — routing, solver gate, and the operating rules a run must honour
