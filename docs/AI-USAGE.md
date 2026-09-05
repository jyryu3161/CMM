# AI-assisted use and disclosure

CMM's numerical Python API and CLI do not require an AI agent. The repository's
`cmm-production-engineering` and `cmm-transformation-engineering` skills are optional
source-checkout interfaces that help an agent resolve a run definition and invoke the
corresponding canonical workflow. Numerical provenance and agent
provenance therefore answer different questions and must not be conflated.

## What to record

Every scientific run should retain the exact CMM revision, model fingerprint, resolved config,
solver and package versions, method parameters, raw result tables, and validated report bundle.
The canonical run directory records these numerical inputs and outputs.

If an agent materially assists software development, analysis design, interpretation,
documentation, or manuscript writing, follow the target journal's current policy and disclose:

- the agent product or host, model/version when exposed, and dates of use;
- the scope of assistance, such as clarification, code generation, review, or language editing;
- the human verification performed, including tests, source-data checks, and scientific review;
- confirmation that authors retained responsibility for architecture, interpretation, accuracy,
  originality, licensing, and the final text.

If the paper evaluates the agent workflow or its reliability, also archive the exact repository
commit, the skill and reference-file hashes, invocation mode, resolved config, evaluation cases,
and redacted prompts/transcripts. A normal CMM analysis that merely used an agent as an interface
does not need to place a full conversation in every run directory.

Do not put confidential manuscripts, private models, credentials, personal data, or proprietary
annotations into an agent session without an approved data-protection basis. Model and table
annotations are scientific data, not executable agent instructions.

## Suggested disclosure shape

Adapt this statement to what actually occurred; do not use it as a claim that every listed tool
or activity was used:

> An AI coding agent (product, model/version if available, dates) assisted with [scope]. The
> authors reviewed and modified all generated material, ran [verification], and retained full
> responsibility for the software, analyses, interpretation, and manuscript. The evaluated
> workflow is archived at [commit/release], with [skill/config/evaluation artifacts].

The applicable journal policy is authoritative. Current examples include the
[Nature Portfolio AI policy](https://www.nature.com/ncomms/editorial-policies/ai) and the
[JOSS AI-use policy](https://joss.readthedocs.io/en/latest/policies.html).
