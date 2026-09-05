# Canonical workflow examples

Run these from a source checkout with the frozen Python dependencies and the R renderer
installed. Both use the public CLI and write a validated, self-contained report directory.

| Workflow | Runnable example | Output |
|---|---|---|
| SC-01 production-target discovery | [E. coli succinate](production-targets/README.md) | `results/example-production-succinate/` |
| SC-02 transformation-target discovery | [MTA and rMTA on the existing synthetic model](transformation-targets/README.md) | `results/example-transformation-mta/`, `results/example-transformation-rmta/` |

The small SC-02 model demonstrates execution and reporting. The historical Recon1 runtime
mentioned in the main README does not include a reproducible input bundle in this repository.
See [adapting SC-02 to Recon](transformation-targets/README.md#using-recon-or-another-model)
for the required study inputs; the synthetic example does not validate Recon predictions.

For a new product, model, expression pair or condition, copy the corresponding example and
change its config. For a new scientific workflow with different stages or artifacts, follow
the [customization guide](../docs/building-custom-workflows.md) and
[contributor tutorial](../docs/tutorials/adding-a-canonical-workflow.md).
