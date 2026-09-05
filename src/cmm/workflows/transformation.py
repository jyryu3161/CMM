"""Transformation-target discovery: rank knockouts that move one metabolic state to another.

Where :mod:`cmm.workflows.production` pushes flux toward a product, this workflow asks the
inverse question — given expression for a *source* state and a *target* state, which single
knockout best moves the source toward the target? The methods are Yizhak et al.'s MTA (2013)
and Valcárcel et al.'s rMTA (2019).

Six stages, mirroring the production workflow's numbering so a reader of one recognises the
other:

``01_preflight``    inputs, gene coverage, solver capability
``02_reference``    the source flux state v_ref
``03_direction``    which reactions must change, and which way
``04_candidates``   the knockout universe
``05_transformation`` one optimisation per candidate, ranked
``06_validation``   the MOMA baseline the source paper compares against

**Two things this workflow does differently from the papers, by construction.** They are
recorded in the run's provenance rather than hidden, and the report is required to state them.

1. v_ref comes from E-Flux2 or LAD. The papers use iMAT, which CMM does not implement. iMAT
   places no objective on growth, whereas E-Flux2 at ``objective_fraction=1.0`` forces a
   growth-maximal state, so the ranking is conditioned on whichever estimator is used.
2. ``epsilon`` is a fixed scalar. The papers derive it per data set from the sampled reference
   distribution. It is measured in the model's own flux units, so there is no safe default;
   :meth:`TransformationWorkflowConfig.suggest_epsilon` computes candidates from v_ref instead.
A reaction-level run additionally collapses coupled reactions to one representative each, using
full flux coupling; see :mod:`cmm.features.coupling`. A gene-level run, the default, instead
deduplicates genes sharing a blocked-reaction signature.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence, cast

from cmm.core.condition import Condition
from cmm.core.media import Medium

if TYPE_CHECKING:  # heavy imports stay out of module import time
    import pandas as pd

ReferenceMethod = Literal["eflux2", "lad"]
SignificanceTest = Literal["ttest", "fold_change"]
TransformationMethod = Literal["mta", "rmta"]
PerturbationLevel = Literal["gene", "reaction"]
ChangedRanking = Literal["p_value", "fold_change"]

#: Yizhak et al. keep "the top 100-200 most differentially expressed reactions" as the changed
#: set. Each one adds a binary variable to the MIQP, so the cut is also what decides whether a
#: genome-scale run finishes.
PUBLISHED_CHANGED_SET_RANGE = (100, 200)

# Same per-replicate transform used by the desktop expression comparison.
_EXPRESSION_PSEUDOCOUNT = 1.0


@dataclass(frozen=True)
class CandidateConfig:
    """How the knockout universe is built.

    The candidate count is the denominator of any "ranked in the top *N*%" statement, so the
    construction is configuration rather than an implementation detail, and it is reported.

    Coupled-set reduction is defined on reactions: Yizhak et al. take "a member from each
    partially coupled set". It therefore applies only when ``perturbation="reaction"``. A
    gene-level run instead deduplicates genes that block the same reaction signature, which is
    the convention the production workflow already uses for its single-knockout screen.
    """

    exclude_blocked: bool = True
    exclude_essential: bool = True
    essential_growth_fraction: float = 0.20  # paper: growth reduced by more than 80%
    #: ``None`` follows the perturbation level, which is what makes the two defaults coherent:
    #: on for a reaction-level run, off for a gene-level one. Set it explicitly to override.
    collapse_coupled_sets: bool | None = None
    explicit: tuple[str, ...] | None = None  # skip construction; use these ids verbatim

    def collapse_for(self, perturbation: PerturbationLevel) -> bool:
        if self.collapse_coupled_sets is None:
            return perturbation == "reaction"
        return self.collapse_coupled_sets

    def validate(self, perturbation: PerturbationLevel) -> None:
        if not 0.0 <= self.essential_growth_fraction < 1.0:
            raise ValueError("essential_growth_fraction must be in [0, 1)")
        if self.collapse_coupled_sets is True and perturbation == "gene":
            raise ValueError(
                "collapse_coupled_sets applies to reaction-level knockouts; a gene-level run "
                "deduplicates on blocked-reaction signature instead. Set perturbation="
                "'reaction' to reproduce the published candidate set, or leave "
                "collapse_coupled_sets unset to follow the perturbation level."
            )
        if self.explicit is not None and not self.explicit:
            raise ValueError("candidates.explicit must be None or a non-empty list")


@dataclass(frozen=True)
class DirectionConfig:
    """How the two expression states become a per-reaction desired direction of change."""

    significance: SignificanceTest = "ttest"
    p_value_cutoff: float = 0.05  # paper: P < 0.05
    up_threshold: float = 1.0  # used only when significance="fold_change"
    down_threshold: float = 1.0
    top_n_changed: int | None = 200
    ranking: ChangedRanking = "p_value"

    def validate(self) -> None:
        if not 0.0 < self.p_value_cutoff <= 1.0:
            raise ValueError("p_value_cutoff must be in (0, 1]")
        if self.up_threshold < 0 or self.down_threshold < 0:
            raise ValueError("fold-change thresholds must be non-negative")
        if self.top_n_changed is not None and self.top_n_changed < 1:
            raise ValueError("top_n_changed must be None or at least 1")
        if self.ranking == "p_value" and self.significance != "ttest":
            raise ValueError(
                "ranking='p_value' needs significance='ttest'; a fold-change run has no "
                "P values to rank on"
            )


@dataclass(frozen=True)
class TransformationValidationConfig:
    """Cross-checks run after the ranking, not further predictions."""

    enabled: bool = True
    #: Yizhak et al. compare MTA against a MOMA baseline and report it as markedly inferior.
    #: Reproducing that contrast is what shows the ranking's signal comes from the method.
    run_moma_baseline: bool = True
    #: Empty means no sweep. The papers report epsilon sensitivity (Yizhak Supp. Fig. S6);
    #: with no way to derive epsilon, reporting the sweep is the honest substitute.
    epsilon_sweep: tuple[float, ...] = ()

    def validate(self) -> None:
        if any(value < 0 for value in self.epsilon_sweep):
            raise ValueError("epsilon_sweep values must be non-negative")


@dataclass(frozen=True)
class TransformationWorkflowConfig:
    """Complete, serializable invocation of transformation-target discovery.

    ``source`` is the state to move *away* from and ``target`` the state to move *toward*.
    Swapping them is a different scientific question with a different answer, and nothing in
    the model can detect the mistake — the agent skill is required to confirm the direction
    rather than infer it from file names.

    Both expression files contain non-negative linear measurements. Reference fluxes use
    their linear replicate means; direction tests use log2(value + 1) per replicate.
    """

    model_path: str | Path
    source_expression_path: str | Path
    target_expression_path: str | Path
    output_dir: str | Path | None = None
    solver: str | None = None
    medium: Medium | str | None = None
    condition: Condition | None = None

    reference_method: ReferenceMethod = "eflux2"
    reference_objective_fraction: float = 1.0

    direction: DirectionConfig = field(default_factory=DirectionConfig)
    candidates: CandidateConfig = field(default_factory=CandidateConfig)

    method: TransformationMethod = "mta"
    perturbation: PerturbationLevel = "gene"
    alpha: float = 0.66  # paper's main-text value, robust over 0.1-0.9
    epsilon: float = 1e-3
    parameter_k: float = 100.0  # rMTA Equation 9

    validation: TransformationValidationConfig = field(
        default_factory=TransformationValidationConfig
    )
    overwrite: bool = False

    def __post_init__(self) -> None:
        for name in ("model_path", "source_expression_path", "target_expression_path"):
            object.__setattr__(self, name, Path(getattr(self, name)))
        if self.output_dir is not None:
            object.__setattr__(self, "output_dir", Path(self.output_dir))
        self.validate()

    def validate(self) -> None:
        for name in ("model_path", "source_expression_path", "target_expression_path"):
            if not str(getattr(self, name)):
                raise ValueError(f"{name} must not be empty")
        if self.source_expression_path == self.target_expression_path:
            raise ValueError(
                "source and target expression must be different files; the transformation "
                "from a state to itself is not defined"
            )
        if self.reference_method not in ("eflux2", "lad"):
            raise ValueError("reference_method must be 'eflux2' or 'lad'")
        if not 0.0 < self.reference_objective_fraction <= 1.0:
            raise ValueError("reference_objective_fraction must be in (0, 1]")
        if self.method not in ("mta", "rmta"):
            raise ValueError("method must be 'mta' or 'rmta'")
        if self.perturbation not in ("gene", "reaction"):
            raise ValueError("perturbation must be 'gene' or 'reaction'")
        if not 0.0 <= self.alpha <= 1.0:
            raise ValueError("alpha must be between 0 and 1")
        if self.epsilon < 0:
            raise ValueError("epsilon must be non-negative")
        if self.parameter_k <= 0:
            raise ValueError("parameter_k must be positive")
        self.direction.validate()
        self.candidates.validate(self.perturbation)
        self.validation.validate()

    @property
    def follows_published_changed_set_size(self) -> bool:
        """Whether the changed-set cut sits in the range Yizhak et al. report using."""

        low, high = PUBLISHED_CHANGED_SET_RANGE
        return self.top_n_changed is not None and low <= self.top_n_changed <= high

    @property
    def top_n_changed(self) -> int | None:
        return self.direction.top_n_changed

    @staticmethod
    def suggest_epsilon(reference_fluxes: Mapping[str, float]) -> dict[str, float]:
        """Epsilon candidates derived from the reference state's own flux magnitudes.

        Yizhak et al. choose epsilon from a sampled distribution of v_ref, which needs the
        sampling CMM avoids for determinism. What can still be honoured is the property that
        makes their choice work: epsilon is a *flux magnitude*, so it must be read off the
        model at hand rather than carried over from another one.

        Returns percentiles of the non-zero |v_ref| so a caller — or an agent interview — can
        offer a number with a basis instead of a default nobody chose.
        """

        import numpy as np

        magnitudes = np.array(
            [abs(v) for v in reference_fluxes.values() if abs(v) > 1e-12]
        )
        if not magnitudes.size:
            return {}
        return {
            "p10": float(np.percentile(magnitudes, 10)),
            "p25": float(np.percentile(magnitudes, 25)),
            "median": float(np.median(magnitudes)),
            "p75": float(np.percentile(magnitudes, 75)),
        }

    def to_provenance(self) -> dict[str, object]:
        """Everything a reader needs to tell this run apart from a differently-configured one."""

        return {
            "model_path": str(self.model_path),
            "source_expression": str(self.source_expression_path),
            "target_expression": str(self.target_expression_path),
            "reference_method": self.reference_method,
            "reference_objective_fraction": self.reference_objective_fraction,
            "expression_scale": "linear",
            "direction_transform": "log2(value + pseudocount)",
            "direction_pseudocount": _EXPRESSION_PSEUDOCOUNT,
            "significance": self.direction.significance,
            "p_value_cutoff": self.direction.p_value_cutoff,
            "up_threshold": self.direction.up_threshold,
            "down_threshold": self.direction.down_threshold,
            "top_n_changed": self.direction.top_n_changed,
            "changed_ranking": self.direction.ranking,
            "method": self.method,
            "perturbation": self.perturbation,
            "alpha": self.alpha,
            "epsilon": self.epsilon,
            "parameter_k": self.parameter_k,
            "candidate_exclude_blocked": self.candidates.exclude_blocked,
            "candidate_exclude_essential": self.candidates.exclude_essential,
            "candidate_essential_growth_fraction": self.candidates.essential_growth_fraction,
            "candidate_collapse_coupled_sets": self.candidates.collapse_for(
                self.perturbation
            ),
            "candidate_explicit": list(self.candidates.explicit)
            if self.candidates.explicit
            else None,
            "epsilon_sweep": list(self.validation.epsilon_sweep),
            "run_moma_baseline": self.validation.run_moma_baseline,
            "medium": self.medium if isinstance(self.medium, str) else None,
            "condition": self.condition.name if self.condition else None,
            # Stated on every run so a reader never has to infer it from the method name.
            "reference_state_deviation": (
                "v_ref is a deterministic "
                f"{self.reference_method} solve; Yizhak et al. (2013) use iMAT, which CMM "
                "does not implement"
            ),
        }

    @classmethod
    def from_json(cls, path: str | Path) -> TransformationWorkflowConfig:
        """Load from a UTF-8 JSON object, resolving relative paths against the config file."""

        config_path = Path(path).expanduser().resolve()
        payload = json.loads(config_path.read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            raise ValueError("transformation workflow JSON must contain an object")
        values = dict(payload)
        for name in (
            "model_path",
            "source_expression_path",
            "target_expression_path",
            "output_dir",
        ):
            raw = values.get(name)
            if raw is None:
                continue
            candidate = Path(str(raw)).expanduser()
            if not candidate.is_absolute():
                candidate = config_path.parent / candidate
            values[name] = candidate.resolve()
        return cls.from_mapping(values)

    @classmethod
    def from_mapping(
        cls, payload: Mapping[str, object]
    ) -> TransformationWorkflowConfig:
        """Build a config from CLI-friendly mappings, including the nested sections."""

        from cmm.workflows.production import (
            _condition_from_payload,
            _medium_from_payload,
        )

        values = dict(payload)
        values["medium"] = _medium_from_payload(values.get("medium"))
        values["condition"] = _condition_from_payload(values.get("condition"))
        for key, factory in (
            ("direction", DirectionConfig),
            ("candidates", CandidateConfig),
            ("validation", TransformationValidationConfig),
        ):
            section = values.get(key)
            if isinstance(section, Mapping):
                section_values = dict(section)
                for tuple_field in ("explicit", "epsilon_sweep"):
                    raw = section_values.get(tuple_field)
                    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
                        section_values[tuple_field] = tuple(raw)
                values[key] = factory(**section_values)  # type: ignore[arg-type]
        return cls(**values)  # type: ignore[arg-type]


# --- execution ---------------------------------------------------------------------------

_STAGE_DIRECTORIES = (
    "01_preflight",
    "02_reference",
    "03_direction",
    "04_candidates",
    "05_transformation",
    "06_validation",
    "model",
    "inputs",
    "scripts",
)


class TransformationWorkflowError(RuntimeError):
    """The transformation workflow could not run as configured."""


@dataclass(frozen=True)
class PreflightRecord:
    check: str
    status: Literal["pass", "warning"]
    value: object
    message: str


@dataclass(frozen=True)
class TransformationWorkflowResult:
    """Everything one transformation run produced, plus where it was written."""

    config: TransformationWorkflowConfig
    provenance: Mapping[str, object]
    preflight: tuple[PreflightRecord, ...]
    reference_fluxes: Mapping[str, float]
    direction: Mapping[str, int]
    candidates: tuple[str, ...]
    ranking: tuple[Mapping[str, object], ...]
    moma_baseline: tuple[Mapping[str, object], ...] = ()
    epsilon_sweep: tuple[Mapping[str, object], ...] = ()
    candidate_filtering: Mapping[str, object] = field(default_factory=dict)
    direction_summary: Mapping[str, object] = field(default_factory=dict)
    run_directory: Path | None = None

    def summary(self) -> dict[str, object]:
        changed = sum(1 for value in self.direction.values() if value != 0)
        best = self.ranking[0] if self.ranking else None
        return {
            "method": self.config.method,
            "perturbation": self.config.perturbation,
            "n_candidates": len(self.candidates),
            "n_scored": len(self.ranking),
            "n_changed_reactions": changed,
            "reference_method": self.config.reference_method,
            "top_target": best["target_id"] if best else None,
            "top_score": best["score"] if best else None,
            # The count is the denominator of any "ranked in the top N%" reading, so it
            # travels with the result rather than being recomputed by a reader.
            "candidate_construction": dict(self.candidate_filtering),
            "direction_construction": dict(self.direction_summary),
        }


def _read_expression(path: Path) -> "pd.DataFrame":
    """Read a gene-expression table: first column gene ids, the rest measurements.

    Two columns is one measurement per gene, which supports a fold-change cut only. Three or
    more are replicates, which is what the published Student's t-test needs. The shape is a
    fact about the file, so it is inspected and reported rather than asked about.
    """

    import numpy as np
    import pandas as pd

    frame = pd.read_csv(path, sep=None, engine="python", index_col=0)
    numeric = frame.select_dtypes("number")
    if numeric.empty:
        raise TransformationWorkflowError(
            f"{path.name} has no numeric expression columns after the gene id column"
        )
    if numeric.isna().any().any():
        raise TransformationWorkflowError(f"{path.name} carries missing values")
    if not np.isfinite(numeric.to_numpy()).all() or (numeric < 0).any().any():
        raise TransformationWorkflowError(
            f"{path.name} must contain finite, non-negative linear expression values"
        )
    numeric.index = [str(value) for value in numeric.index]
    if numeric.index.has_duplicates:
        raise TransformationWorkflowError(f"{path.name} repeats a gene id")
    return numeric


def _gene_directions(
    source: "pd.DataFrame",
    target: "pd.DataFrame",
    config: DirectionConfig,
) -> "pd.DataFrame":
    """Compare linear expression in log2 space, preserving linear reference inputs."""

    import numpy as np

    from cmm.omics.differential import (
        gene_directions_by_fold_change,
        gene_directions_from_replicates,
    )

    source_log2 = np.log2(source + _EXPRESSION_PSEUDOCOUNT)
    target_log2 = np.log2(target + _EXPRESSION_PSEUDOCOUNT)
    try:
        if config.significance == "ttest":
            return gene_directions_from_replicates(
                source_log2, target_log2, p_value_cutoff=config.p_value_cutoff
            )
        return gene_directions_by_fold_change(
            source_log2,
            target_log2,
            up_threshold=config.up_threshold,
            down_threshold=config.down_threshold,
        )
    except ValueError as error:
        raise TransformationWorkflowError(str(error)) from error


def _build_candidates(
    model,
    config: TransformationWorkflowConfig,
) -> tuple[tuple[str, ...], dict[str, object]]:
    """The knockout universe, and a record of how many reactions each filter removed.

    Yizhak et al. drop dead-end reactions, drop essential ones, and then keep one member per
    coupled set. That reduction is what makes each candidate a distinct intervention; without
    it a linear pathway contributes one intervention several times and inflates the
    denominator of any "top *N*%" reading. Every step's count is returned because the report
    is required to state how the candidate set was built.
    """

    from cobra.flux_analysis import find_blocked_reactions, single_reaction_deletion

    from cmm.features import coupled_reaction_sets
    from cmm.features._perturbation import gene_perturbations

    if config.candidates.explicit is not None:
        chosen = tuple(config.candidates.explicit)
        return chosen, {
            "source": "explicit",
            "n_candidates": len(chosen),
            "perturbation": config.perturbation,
        }

    record: dict[str, object] = {
        "source": "constructed",
        "perturbation": config.perturbation,
        "n_reactions": len(model.reactions),
    }
    wild_type = model.slim_optimize()
    record["wild_type_growth"] = float(wild_type) if wild_type == wild_type else None

    open_reactions = [reaction.id for reaction in model.reactions]
    if config.candidates.exclude_blocked:
        # Both this and the deletion screen below are pinned to one process. COBRApy's
        # defaults fork a pool, and under spawn-based multiprocessing that re-imports the
        # caller — any script without an ``if __name__ == '__main__'`` guard then recurses
        # instead of running. A workflow has to be callable from anywhere.
        blocked = set(find_blocked_reactions(model, processes=1))
        open_reactions = [rid for rid in open_reactions if rid not in blocked]
        record["n_blocked_removed"] = len(blocked)

    active_reactions = set(open_reactions)
    threshold = None
    if config.candidates.exclude_essential and wild_type == wild_type and wild_type > 0:
        deletion = single_reaction_deletion(
            model, reaction_list=open_reactions, processes=1
        )
        growth = {
            next(iter(ids)): (0.0 if value != value else float(value))
            for ids, value in zip(deletion["ids"], deletion["growth"])
        }
        threshold = config.candidates.essential_growth_fraction * wild_type
        essential = {rid for rid, value in growth.items() if value < threshold}
        open_reactions = [rid for rid in open_reactions if rid not in essential]
        record["n_essential_removed"] = len(essential)
        record["essential_growth_threshold"] = float(threshold)
    record["n_reactions_allowed"] = len(open_reactions)

    if config.perturbation == "reaction":
        if config.candidates.collapse_for("reaction"):
            sets = coupled_reaction_sets(model, open_reactions)
            record["coupling"] = sets.to_provenance()
            return tuple(sets.representatives), record
        return tuple(sorted(open_reactions)), record

    # A candidate filter does not change the phenotype of a gene deletion. Keep its full
    # signature, then test each distinct deletion for growth: jointly blocking individually
    # nonessential reactions can also be lethal.
    signatures: dict[tuple[str, ...], list[str]] = {}
    inert = 0
    for perturbation in gene_perturbations(model):
        signature = tuple(sorted(perturbation.reaction_ids))
        if not active_reactions.intersection(signature):
            inert += 1
            continue
        signatures.setdefault(signature, []).append(perturbation.target_id)
    gene_groups: dict[str, list[str]] = {}
    essential_genes = 0
    for signature, genes in signatures.items():
        if threshold is not None:
            with model:
                for rid in signature:
                    model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                growth = model.slim_optimize()
                viable = (
                    model.solver.status == "optimal"
                    and math.isfinite(growth)
                    and growth > 0.0
                    and growth >= threshold
                )
            if not viable:
                essential_genes += len(genes)
                continue
        represented = sorted(genes)
        gene_groups[represented[0]] = represented
    record["n_genes_inert"] = inert
    record["n_genes_essential_removed"] = essential_genes
    record["n_distinct_blocked_signatures"] = len(gene_groups)
    record["represented_genes"] = gene_groups
    return tuple(sorted(gene_groups)), record


def run_transformation_target_discovery(
    config: TransformationWorkflowConfig,
) -> TransformationWorkflowResult:
    """Run the six stages and, when ``output_dir`` is set, write a schema-v2 run bundle."""

    from cobra.io import read_sbml_model

    from cmm.core import solvers
    from cmm.core.provenance import run_provenance
    from cmm.features.revert import revert_targets
    from cmm.omics.differential import (
        reaction_directions,
        restrict_to_top_changed,
    )
    from cmm.omics.expression import integrate_expression

    model = read_sbml_model(str(config.model_path))
    if config.solver:
        model.solver = config.solver

    preflight: list[PreflightRecord] = []

    def check(name: str, ok: bool, value: object, message: str) -> None:
        preflight.append(
            PreflightRecord(name, "pass" if ok else "warning", value, message)
        )

    # --- 01_preflight -------------------------------------------------------------------
    # The MIQP gate is a stop, not a warning: rmta_continuous is a different method and must
    # never stand in for a published one that the solver cannot run.
    solvers.require(
        "MIQP", model.solver.interface, feature=f"published {config.method.upper()}"
    )

    source_expression = _read_expression(Path(config.source_expression_path))
    target_expression = _read_expression(Path(config.target_expression_path))
    check(
        "replicates",
        min(source_expression.shape[1], target_expression.shape[1]) >= 2,
        [source_expression.shape[1], target_expression.shape[1]],
        "at least two replicates per state are needed for the published t-test",
    )

    model_genes = {gene.id for gene in model.genes}
    measured = set(source_expression.index) & set(target_expression.index)
    covered = model_genes & measured
    check(
        "gene_coverage",
        bool(model_genes) and len(covered) / len(model_genes) >= 0.5,
        {"model_genes": len(model_genes), "measured_in_model": len(covered)},
        "expression must cover the model's genes; a small overlap usually means the two use "
        "different identifier systems",
    )
    if not covered:
        raise TransformationWorkflowError(
            "no model gene appears in both expression files; check the identifier system"
        )

    with model:
        if config.medium is not None:
            from cmm.core.media import apply_medium

            apply_medium(model, config.medium)
        if config.condition is not None:
            config.condition.apply_to(model)

        growth = model.slim_optimize()
        check(
            "growth",
            growth == growth and growth > 1e-9,
            float(growth) if growth == growth else None,
            "the model must grow under the requested condition",
        )

        # --- 02_reference ---------------------------------------------------------------
        extra: dict[str, Any] = (
            {"objective_fraction": config.reference_objective_fraction}
            if config.reference_method == "eflux2"
            else {}
        )
        omics = integrate_expression(
            model,
            source_expression.mean(axis=1).to_dict(),
            method=config.reference_method,
            **extra,
        )
        if omics.status != "optimal" or not omics.fluxes:
            raise TransformationWorkflowError(
                f"the source expression produced no valid reference state with "
                f"{config.reference_method} ({omics.status})"
            )
        reference = omics.to_flux_state(f"source_state_{config.reference_method}")

        # --- 03_direction ---------------------------------------------------------------
        genes = _gene_directions(source_expression, target_expression, config.direction)
        gene_dirs = {gene: int(value) for gene, value in genes["direction"].items()}
        direction = reaction_directions(model, gene_dirs, reference=reference)
        n_labelled = len(direction.nonsteady())
        if config.direction.top_n_changed is not None:
            direction = restrict_to_top_changed(
                model,
                direction,
                genes["log2_fold_change"].to_dict(),
                config.direction.top_n_changed,
                gene_p_values=genes["p_value"].to_dict()
                if config.direction.ranking == "p_value"
                else None,
            )
        if not direction.nonsteady():
            raise TransformationWorkflowError(
                "no reaction was labelled as changed; loosen the significance threshold or "
                "check that the two expression states differ"
            )
        direction_summary = {
            "significance": config.direction.significance,
            "n_genes_compared": int(len(genes)),
            "n_genes_significant": int(genes["significant"].sum()),
            "n_reactions_labelled": n_labelled,
            "n_reactions_changed": len(direction.nonsteady()),
            "ranking": config.direction.ranking,
            "top_n_changed": config.direction.top_n_changed,
            **dict(direction.metadata),
        }

        # --- 04_candidates --------------------------------------------------------------
        candidates, filtering = _build_candidates(model, config)
        if not candidates:
            raise TransformationWorkflowError(
                "the candidate filters removed every knockout; relax exclude_blocked or "
                "exclude_essential"
            )

        # --- 05_transformation ----------------------------------------------------------
        ranking = revert_targets(
            model,
            None,
            reference,
            direction,
            targets=candidates,
            method=config.method,
            alpha=config.alpha,
            epsilon=config.epsilon,
            parameter_k=config.parameter_k,
            perturbation=config.perturbation,
        )
        # Only rMTA computes bTS, mTS and wTS as separate quantities. For MTA the same score
        # is returned in all four slots, so emitting them as columns would present one number
        # as three independent measurements that happen to agree.
        components = ("bTS", "mTS", "wTS") if config.method == "rmta" else ()
        # The name column appears only when the model named something, so a reconstruction
        # that names nothing writes exactly the table it wrote before.
        named = any(target.target_name for target in ranking.targets)
        rows = tuple(
            {
                "target_id": target.target_id,
                **({"target_name": target.target_name} if named else {}),
                "score": float(target.score),
                "rank": index + 1,
                **{
                    key: float((target.detail or {})[key])
                    for key in components
                    if isinstance((target.detail or {}).get(key), (int, float))
                },
            }
            for index, target in enumerate(ranking.sorted().targets)
        )

        # --- 06_validation --------------------------------------------------------------
        baseline: tuple[Mapping[str, object], ...] = ()
        sweep: tuple[Mapping[str, object], ...] = ()
        if config.validation.enabled:
            if config.validation.run_moma_baseline:
                baseline = _moma_baseline(
                    model,
                    reference,
                    target_expression,
                    candidates,
                    config,
                )
            sweep = _epsilon_sweep(
                model, reference, direction, candidates, config, rows
            )

    provenance = {
        **run_provenance(model, method=f"transformation_{config.method}"),
        **config.to_provenance(),
        "reference_provenance": reference.provenance,
        "direction_provenance": dict(direction.metadata),
        "ranking_provenance": dict(ranking.metadata),
    }
    result = TransformationWorkflowResult(
        config=config,
        provenance=provenance,
        preflight=tuple(preflight),
        reference_fluxes=dict(reference.fluxes),
        direction=dict(direction.directions),
        candidates=candidates,
        ranking=rows,
        moma_baseline=baseline,
        epsilon_sweep=sweep,
        candidate_filtering=filtering,
        direction_summary=direction_summary,
        run_directory=None,
    )
    if config.output_dir is None:
        return result
    return _export(result, genes, model)


def _write_reproduction_scripts(
    writer: Any,
    result: TransformationWorkflowResult,
    *,
    root: Path,
    model_relative: str,
    expression_relatives: Mapping[str, str],
) -> None:
    """Write the three entry points that make the directory replayable on its own.

    Months later the run directory is all that survives; a bundle that cannot say how to
    re-run, re-render or re-check itself leaves that to whoever finds it.
    """

    from dataclasses import asdict

    from cmm.workflows._bundle import _jsonable

    config_payload = cast(dict[str, object], _jsonable(asdict(result.config)))
    config_payload["model_path"] = f"../{model_relative}"
    config_payload.update(
        {field: f"../{relative}" for field, relative in expression_relatives.items()}
    )
    config_payload["output_dir"] = f"../../{root.name}__reproduced"
    config_payload["overwrite"] = False
    writer.json(
        "scripts/transformation_config.json",
        config_payload,
        stage="scripts",
        role="reproduction_config",
    )
    writer.text(
        "scripts/reproduce.py",
        """#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

from cmm.workflows.transformation import (
    TransformationWorkflowConfig,
    run_transformation_target_discovery,
)

HERE = Path(__file__).resolve().parent


def main() -> int:
    parser = argparse.ArgumentParser(description="Reproduce this CMM transformation run")
    parser.add_argument("--output-dir", type=Path, help="Override the sibling output directory")
    args = parser.parse_args()
    config = TransformationWorkflowConfig.from_json(HERE / "transformation_config.json")
    if args.output_dir is not None:
        config = replace(config, output_dir=args.output_dir.resolve(), overwrite=False)
    result = run_transformation_target_discovery(config)
    print(result.run_directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
""",
        stage="scripts",
        role="reproduce_script",
        executable=True,
    )
    writer.text(
        "scripts/render.py",
        """#!/usr/bin/env python3
from pathlib import Path

from cmm.reporting import render_transformation_report

RUN_DIR = Path(__file__).resolve().parents[1]
report = render_transformation_report(RUN_DIR)
print(report.report_standalone_html)
""",
        stage="scripts",
        role="render_script",
        executable=True,
    )
    writer.text(
        "scripts/validate.py",
        """#!/usr/bin/env python3
import json
from pathlib import Path

from cmm.reporting import validate_transformation_run

RUN_DIR = Path(__file__).resolve().parents[1]
report = validate_transformation_run(RUN_DIR)
print(json.dumps({"valid": report.valid, "issues": report.issues, "warnings": report.warnings}, indent=2))
raise SystemExit(0 if report.valid else 1)
""",
        stage="scripts",
        role="validate_script",
        executable=True,
    )


def _moma_baseline(model, reference, target_expression, candidates, config):
    """The MOMA comparison Yizhak et al. report as markedly inferior to MTA.

    It needs a target *flux* vector rather than a direction map, so the target expression is
    put through the same estimator that produced v_ref. Reproducing the contrast is what shows
    a ranking's signal comes from the method rather than from the inputs.
    """

    from cmm.features.transformation import transformation_targets
    from cmm.omics.expression import integrate_expression

    extra = (
        {"objective_fraction": config.reference_objective_fraction}
        if config.reference_method == "eflux2"
        else {}
    )
    omics = integrate_expression(
        model,
        target_expression.mean(axis=1).to_dict(),
        method=config.reference_method,
        **extra,
    )
    if omics.status != "optimal" or not omics.fluxes:
        return ()
    target_state = omics.to_flux_state(f"target_state_{config.reference_method}")
    ranking = transformation_targets(
        model,
        reference,
        target_state,
        method="moma",
        perturbation=config.perturbation,
        targets=list(candidates),
        order=2,
    )
    named = any(t.target_name for t in ranking.targets)
    return tuple(
        {
            "target_id": t.target_id,
            **({"target_name": t.target_name} if named else {}),
            "moma_score": float(t.score),
            "status": t.status,
            "rank": i + 1,
        }
        for i, t in enumerate(ranking.sorted().targets)
    )


def _epsilon_sweep(model, reference, direction, candidates, config, base_rows):
    """Re-rank at each swept epsilon and report how far each top candidate moves.

    Epsilon is a flux magnitude with no derivable value here, so the sensitivity is the honest
    substitute for the paper's derivation — both source papers report one.
    """

    from cmm.features.revert import revert_targets

    if not config.validation.epsilon_sweep:
        return ()
    base_rank = {row["target_id"]: row["rank"] for row in base_rows}
    out: list[Mapping[str, object]] = []
    for epsilon in config.validation.epsilon_sweep:
        ranking = revert_targets(
            model,
            None,
            reference,
            direction,
            targets=list(candidates),
            method=config.method,
            alpha=config.alpha,
            epsilon=epsilon,
            parameter_k=config.parameter_k,
            perturbation=config.perturbation,
        )
        for index, target in enumerate(ranking.sorted().targets):
            out.append(
                {
                    "epsilon": epsilon,
                    "target_id": target.target_id,
                    "score": float(target.score),
                    "rank": index + 1,
                    "rank_at_configured_epsilon": base_rank.get(target.target_id),
                }
            )
    return tuple(out)


def _export(
    result: TransformationWorkflowResult,
    genes: "pd.DataFrame",
    model,
) -> TransformationWorkflowResult:
    """Write the schema-v2 run bundle: one manifest, one role per artifact."""

    import json as _json
    import shutil
    from dataclasses import asdict, replace

    import pandas as pd
    from cobra.io import write_sbml_model

    from cmm.workflows._bundle import ArtifactRecord, _ArtifactWriter, _jsonable

    config = result.config
    if config.output_dir is None:  # guarded by the caller; kept as an invariant check
        raise TransformationWorkflowError("cannot export a run with no output_dir")
    root = Path(config.output_dir).resolve()
    if root.exists() and not root.is_dir():
        raise FileExistsError(f"output path exists and is not a directory: {root}")
    if root.exists() and any(root.iterdir()) and not config.overwrite:
        raise FileExistsError(
            f"output directory is not empty: {root}; choose a new directory or set "
            "overwrite=True"
        )
    # Capture input bytes before cleanup, including when an in-place rerun reads its own
    # archived model and expression tables.
    model_path = Path(config.model_path)
    model_bytes = model_path.read_bytes() if model_path.is_file() else None
    expression_inputs = {
        field: (Path(path).suffix or ".csv", Path(path).read_bytes())
        for field, path in (
            ("source_expression_path", config.source_expression_path),
            ("target_expression_path", config.target_expression_path),
        )
    }
    if config.overwrite:
        # Remove only workflow-owned paths, so permission to rerun is never permission to
        # delete a user's unrelated files that happen to share the directory.
        for relative in (
            *_STAGE_DIRECTORIES,
            "00_config.json",
            "00_provenance.json",
            "00_summary.json",
            "00_manifest.json",
        ):
            owned = root / relative
            if owned.is_symlink() or owned.is_file():
                owned.unlink()
            elif owned.is_dir():
                shutil.rmtree(owned)
    root.mkdir(parents=True, exist_ok=True)
    for directory in _STAGE_DIRECTORIES:
        (root / directory).mkdir(parents=True, exist_ok=True)

    writer = _ArtifactWriter(root, error_type=TransformationWorkflowError)

    model_relative = f"model/{model_path.name}"
    if model_bytes is not None:
        (root / model_relative).write_bytes(model_bytes)
    else:
        write_sbml_model(model, str(root / model_relative))
    writer.existing(
        model_relative,
        stage="model",
        role="model",
        media_type="application/sbml+xml",
    )

    expression_relatives: dict[str, str] = {}
    for input_field, (suffix, payload) in expression_inputs.items():
        role = input_field.removesuffix("_path")
        relative = f"inputs/{role}{suffix}"
        (root / relative).write_bytes(payload)
        expression_relatives[input_field] = relative
        writer.existing(
            relative,
            stage="inputs",
            role=role,
            media_type="text/tab-separated-values" if suffix == ".tsv" else "text/csv",
        )

    writer.csv(
        "01_preflight/preflight.csv",
        pd.DataFrame([asdict(record) for record in result.preflight]),
        stage="01_preflight",
        role="preflight",
    )
    writer.csv(
        "02_reference/source_reference_fluxes.csv",
        pd.DataFrame(
            {
                "reaction_id": list(result.reference_fluxes),
                "flux": list(result.reference_fluxes.values()),
            }
        ),
        stage="02_reference",
        role="source_reference_fluxes",
        method=config.reference_method,
    )
    writer.csv(
        "03_direction/gene_differential_expression.csv",
        genes.reset_index().rename(columns={"index": "gene_id"}),
        stage="03_direction",
        role="gene_differential_expression",
        method=config.direction.significance,
    )
    writer.csv(
        "03_direction/reaction_direction_map.csv",
        pd.DataFrame(
            {
                "reaction_id": list(result.direction),
                "direction": list(result.direction.values()),
            }
        ),
        stage="03_direction",
        role="reaction_direction_map",
    )
    writer.csv(
        "04_candidates/candidates.csv",
        pd.DataFrame({"target_id": list(result.candidates)}),
        stage="04_candidates",
        role="transformation_candidates",
    )
    writer.csv(
        "05_transformation/transformation_ranking.csv",
        pd.DataFrame(list(result.ranking)),
        stage="05_transformation",
        role="transformation_ranking",
        method=config.method,
    )
    writer.csv(
        "06_validation/moma_baseline.csv",
        pd.DataFrame(list(result.moma_baseline)),
        stage="06_validation",
        role="moma_baseline",
        method="moma",
        status="complete" if result.moma_baseline else "skipped",
        reason=None if result.moma_baseline else "run_moma_baseline was disabled",
    )
    writer.csv(
        "06_validation/epsilon_sensitivity.csv",
        pd.DataFrame(list(result.epsilon_sweep)),
        stage="06_validation",
        role="epsilon_sensitivity",
        status="complete" if result.epsilon_sweep else "skipped",
        reason=None
        if result.epsilon_sweep
        else "no epsilon_sweep values were configured",
    )

    exported_config = cast(dict[str, object], _jsonable(asdict(config)))
    exported_config["model_path"] = model_relative
    exported_config.update(expression_relatives)
    exported_config["output_dir"] = "."
    writer.json(
        "00_config.json", exported_config, stage="root", role="workflow_configuration"
    )
    writer.json(
        "00_provenance.json",
        {
            **result.provenance,
            "model_path": model_relative,
            **{
                field.removesuffix("_path"): relative
                for field, relative in expression_relatives.items()
            },
        },
        stage="root",
        role="provenance",
    )
    writer.json("00_summary.json", result.summary(), stage="root", role="summary")

    _write_reproduction_scripts(
        writer,
        result,
        root=root,
        model_relative=model_relative,
        expression_relatives=expression_relatives,
    )

    manifest_record = ArtifactRecord(
        path="00_manifest.json",
        stage="root",
        role="authoritative_artifact_manifest",
        media_type="application/json",
    )
    records = (*writer.records, manifest_record)
    manifest = {
        "schema_version": 2,
        "workflow": "transformation_target_discovery",
        "status": "complete"
        if all(record.status == "complete" for record in writer.records)
        else "partial",
        "artifacts": {
            record.role: {
                key: value for key, value in asdict(record).items() if value is not None
            }
            for record in records
        },
    }
    (root / "00_manifest.json").write_text(
        _json.dumps(_jsonable(manifest), indent=2, sort_keys=True, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    return replace(result, run_directory=root)
