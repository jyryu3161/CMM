# Small SC-01 example for the restricted Gurobi licence

This example runs the complete SC-01 production workflow on a synthetic 26-reaction network
(`mixed_acid_demo`): lumped anaerobic mixed-acid fermentation of glucose to lactate, acetate,
ethanol, formate, and succinate. It exists so the whole workflow can run on every operating
system in CI. MOMA-L2 on `e_coli_core` needs 286 QP variables, above the 200-variable cap of
the restricted licence bundled with `pip install gurobipy`; this network stays below it.

It is an execution check, **not a biological model**. Its numbers say nothing about
*E. coli*; use the [succinate example](../production-targets/README.md) with a full licence for
that.

Run from the repository root after installing the frozen dependencies and R renderer:

```bash
uv run --frozen --all-extras python examples/production-targets-small/prepare.py
uv run --frozen --all-extras cmm production-targets --config examples/production-targets-small/config.json
uv run --frozen --all-extras cmm report validate results/example-production-small --json
```

`prepare.py` only writes the SBML model. Every numerical stage runs through
`cmm production-targets`, with the same stages and parameters as the succinate example:
MOMA-L2, ROOM, OptKnock, RobustKnock, FSEOF, FVSEOF, loop diagnostics, flux-response scans,
and matched sampling.

The network is built so that blocking ethanol formation (`ADHE`) forces the surplus NADH
through the reductive succinate branch without lowering growth. RobustKnock reports that
growth-coupled design, while OptKnock returns no design. Its empty table is a result, kept
as such, and the validator warns about it.
