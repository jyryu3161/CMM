# CMM scenarios

What each installed workflow's numbers mean. **These documents do not run anything** — each
workflow is executed from its skill, which carries the command, the entry point and the
capability gates. Read the section here when a result needs interpreting.

| Workflow | Run it from | Read this for |
|---|---|---|
| Production targets | `.agents/skills/cmm-production-engineering/` | [SC-01](SC-01-production-target-discovery.md) |
| Transformation targets | `.agents/skills/cmm-transformation-engineering/` | [SC-02](SC-02-transformation-target-discovery.md) |

Read `docs/agent-reference.md` for signatures when writing narrow API calls. To change workflow
parameters or compose a downstream study, follow
[Building or customizing a CMM workflow](../building-custom-workflows.md). Contributors adding
another installed workflow should use the
[canonical-workflow tutorial](../tutorials/adding-a-canonical-workflow.md).

## Index

| ID | Goal | Requires | Minimum solver | Key outputs |
|---|---|---|---|---|
| [SC-01](SC-01-production-target-discovery.md) | Increase production of a target metabolite, and design a strain where it is guaranteed | confirmed model, product exchange, condition | **QP + MILP**; strain design also needs importable `straindesign` | MOMA/ROOM single-deletion candidates, OptKnock/RobustKnock designs, FSEOF/FVSEOF targets, forward validation |
| [SC-02](SC-02-transformation-target-discovery.md) | Rank knockouts that move a source metabolic state toward a target state | confirmed model, source and target expression, condition | **MIQP** — no LP or QP substitute exists | ranked candidates with their transformation scores, the MOMA baseline, optional epsilon sensitivity |

Each is scientifically complete on its own and neither is a prerequisite for the other. Both
supply their report and artifact contract automatically.

Shared, used by all of them:

- [`_preflight.md`](_preflight.md) — checks that run before any scenario. Skipping these is how
  a run produces a confident wrong answer instead of an error.
- [`_reporting.md`](_reporting.md) — general reporting principles plus the concrete SC-01
  artifact contract. A different canonical workflow needs its own schema and validator.

## How they fit together

Both workflows share the same preflight and reporting contracts:

```
   _preflight ─┬─►  production skill  ─►  SC-01 numbers  ─┬─► _reporting
               └─►  transformation skill ─► SC-02 numbers ─┘
```

- **SC-01 is the spine of a production goal.** Its single-knockout step evaluates every
  nominated deletion independently with MOMA and ROOM; its separate strain-design step searches
  knockout sets with OptKnock/RobustKnock and proves coupling through `guaranteed_product`.
  These outputs answer different questions and remain separate tables and figure panels.
- **A coupling-only request enters SC-01 at the strain-design step** and skips single-deletion
  and amplification analyses unless the user asks for them.
- **SC-02 asks the inverse question.** Instead of pushing flux toward a product it ranks the
  knockout that best moves one measured state toward another. It is complete on its own and is
  not a step of SC-01.
- **Neither workflow silently downgrades.** A missing QP, MILP, MIQP, importable
  `straindesign`, R, or required R package is a surfaced capability result. A user may
  explicitly approve a narrower run, but its report must state which claim is no longer made.

## Scenario document structure

Every `SC-*.md` uses the same shape so a pipeline can be followed without reading prose:

- **YAML front matter** — id, goal, when to use (including common natural-language phrasings), required and
  optional inputs, solver requirements, step list, expected runtime.
- **Objective** — what the run produces, plus measurable success criteria.
- **Pipeline at a glance** — one row per step: question, method, output.
- **Steps** — each with *Goal*, *Preconditions*, *Call* (runnable code), *Outputs*,
  *Artifacts*, *Decision rule*, *Branch*, *Failure → action*, and *Solver*.
- **Cross-checks** — internal consistency to confirm before writing the report.
- **Do not** — the misuses that produce plausible, wrong conclusions.

The *Decision rule* fields carry numeric thresholds. They are documented starting points, not
community standards: state whichever you used in the report.
