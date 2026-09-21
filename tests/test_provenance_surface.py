"""Every public service that returns numbers must carry the full run-provenance block.

``AGENTS.md`` rule 4 — "every reported number carries its provenance" — is a property of the
*whole* service surface, not of the services that happen to have been checked. The 0.4.0
release added ``timestamp_utc``, ``seed``, ``solver_version`` and ``platform`` to
``run_provenance`` but left the MOMA/ROOM family and ``predict_condition_fluxes`` without any
block at all, which no test caught because no test asked the question across the surface.

This module asks it. :data:`SERVICES` names every public service that returns numbers; the
parametrised test below fails for any of them whose result metadata is missing a field of
:data:`REQUIRED_FIELDS`. A new service is added to the registry when it ships, and a service
that loses its provenance fails here rather than in a benchmark two releases later.

Not in the registry, and why:

- ``gene_log2_fold_change``, ``flux_log_change``, ``sign_flips``, ``tie_structure`` — pure
  functions over dictionaries that never see a model, so there is no model to fingerprint and
  no solver to name.
- ``gene_directions`` / ``reaction_directions`` / ``differential_expression`` — these return
  direction *codes* (+1 / -1 / 0), not measured quantities, and run no solve. Their
  ``DirectionMap.metadata`` records the GPR rule that produced them, and ``revert_targets``
  embeds it under ``direction_provenance`` inside its own full block, which is where the
  numbers derived from them are reported.
- perturbation enumeration (``gene_perturbations`` and friends) — a list of knockouts, not
  numbers; it carries its own coverage counts via ``PerturbationList.provenance()``.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pytest
from cobra import Metabolite, Model, Reaction
from cobra.io import load_model

from cmm.core import fba, fva, pfba, reference_state_pfba
from cmm.features import (
    batch_comparison,
    flux_response,
    fseof,
    fvseof,
    knockout_comparison,
    moma,
    production_envelope,
    random_flux_sampling,
    reaction_perturbations,
    reference_constrained_sampling,
    reference_flux,
    revert_targets,
    room,
    theoretical_yield,
    transformation_targets,
)
from cmm.omics import (
    DirectionMap,
    eflux2,
    gene_to_reaction_weights,
    lad,
    predict_condition_fluxes,
)

#: The fields the 0.4.0 provenance block promises. ``seed`` must be *present*; its value is
#: ``None`` for an unseeded deterministic method, while stochastic/search methods report the
#: seed that was actually forwarded to their backend.
REQUIRED_FIELDS = (
    "model_id",
    "model_sha256",
    "timestamp_utc",
    "seed",
    "solver",
    "solver_version",
    "platform",
    "machine",
    "processor",
    "python",
    "cmm",
    "cobra",
    "numpy",
    "pandas",
    "scipy",
    "parameters",
)

ISO_UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{6}Z$")

#: Every public service that returns numbers. Keep in sync with ``INCLUDED_FEATURES``.
SERVICES = (
    "fba",
    "pfba",
    "fva",
    "theoretical_yield",
    "production_envelope",
    "fseof",
    "fvseof",
    "flux_response",
    "eflux2",
    "lad",
    "predict_condition_fluxes",
    "random_flux_sampling",
    "reference_constrained_sampling",
    "reference_flux",
    "moma",
    "room",
    "knockout_comparison_moma",
    "knockout_comparison_room",
    "batch_comparison",
    "transformation_targets",
    "revert_targets",
    "optknock",
    "robustknock",
)

_QP_SERVICES = {
    "eflux2",
    "moma",
    "knockout_comparison_moma",
    "batch_comparison",
    "transformation_targets",
}
_MIQP_SERVICES = {"revert_targets"}
SERVICE_CASES = tuple(
    pytest.param(service, marks=pytest.mark.requires_qp)
    if service in _QP_SERVICES
    else pytest.param(service, marks=pytest.mark.requires_miqp)
    if service in _MIQP_SERVICES
    else service
    for service in SERVICES
)

SUCC = "EX_succ_e"


@pytest.fixture(scope="module")
def core():
    """One ``e_coli_core`` shared by the sweep; every service restores the model it uses."""

    return load_model("textbook")


def _expression(model) -> dict[str, float]:
    rng = np.random.default_rng(12345)
    return {gene.id: float(rng.uniform(1.0, 100.0)) for gene in model.genes}


def _direction_map(model, target, source) -> DirectionMap:
    directions = {}
    for reaction in model.reactions:
        delta = target.get(reaction.id) - source.get(reaction.id)
        directions[reaction.id] = 1 if delta > 1e-3 else (-1 if delta < -1e-3 else 0)
    return DirectionMap(directions=directions, metadata={"origin": "pfba delta"})


def _comparison_model() -> Model:
    """Small reroutable model for comparison provenance, including ROOM's MILP path."""

    model = Model("provenance_comparison")
    a = Metabolite("A_c", compartment="c")
    b = Metabolite("B_c", compartment="c")
    d = Metabolite("D_c", compartment="c")
    p = Metabolite("P_c", compartment="c")

    def reaction(reaction_id, stoichiometry, gene_rule="", upper_bound=1000.0):
        item = Reaction(reaction_id, lower_bound=0.0, upper_bound=upper_bound)
        item.add_metabolites(stoichiometry)
        item.gene_reaction_rule = gene_rule
        return item

    model.add_reactions(
        [
            reaction("SUP_A", {a: 1.0}, upper_bound=10.0),
            reaction("R1", {a: -1.0, b: 1.0}, "g1"),
            reaction("R2", {b: -1.0, p: 1.0}, "g2"),
            reaction("R3", {b: -1.0, d: 1.0}, "g3"),
            reaction("R5", {d: -1.0, p: 1.0}, "g5"),
            reaction("BIOMASS", {p: -1.0}, "gb"),
        ]
    )
    model.objective = "BIOMASS"
    return model


def _strain_design_model() -> Model:
    """Tiny growth-coupling model for deterministic strain-design provenance checks."""

    model = Model("provenance_strain_design")
    substrate = Metabolite("substrate_c", compartment="c")
    precursor = Metabolite("precursor_c", compartment="c")
    product = Metabolite("product_c", compartment="c")

    def reaction(reaction_id, stoichiometry, gene_rule="", upper_bound=1000.0):
        item = Reaction(reaction_id, lower_bound=0.0, upper_bound=upper_bound)
        item.add_metabolites(stoichiometry)
        item.gene_reaction_rule = gene_rule
        return item

    model.add_reactions(
        [
            reaction("SUPPLY", {substrate: 1.0}, upper_bound=10.0),
            reaction(
                "DIRECT",
                {substrate: -1.0, precursor: 2.0},
                "g_direct",
            ),
            reaction(
                "COUPLED",
                {substrate: -1.0, precursor: 1.0, product: 1.0},
                "g_coupled",
            ),
            reaction("BIOMASS", {precursor: -1.0}, "g_biomass"),
            reaction("PRODUCT", {product: -1.0}),
        ]
    )
    model.objective = "BIOMASS"
    return model


def _directional_metadata(service: str) -> dict[str, object]:
    """Exercise transform/revert provenance on a real reroute small enough for restricted MIQP."""

    model = _comparison_model()
    source = reference_state_pfba(model, name="source")  # pFBA uses the short R2 branch
    with model:
        model.reactions.R2.upper_bound = 0.0
        target = reference_state_pfba(model, name="target")  # forced R3 -> R5 reroute
    targets = [gene.id for gene in model.genes]

    if service == "transformation_targets":
        return dict(
            transformation_targets(
                model, source, target, method="moma", targets=targets
            ).metadata
        )
    if service == "revert_targets":
        return dict(
            revert_targets(
                model,
                None,
                source,
                _direction_map(model, target, source),
                targets=targets,
                method="rmta",
            ).metadata
        )
    raise AssertionError(f"unsupported directional provenance service: {service}")


def _metadata(service: str, model) -> dict[str, object]:
    """Run one service with the smallest arguments that still exercise it."""

    if service == "fba":
        return dict(fba(model).metadata)
    if service == "pfba":
        return dict(pfba(model).metadata)
    if service == "fva":
        return dict(fva(model, processes=1).metadata)
    if service == "theoretical_yield":
        return dict(theoretical_yield(model, SUCC).metadata)
    if service == "production_envelope":
        return dict(production_envelope(model, SUCC, points=5).metadata)
    if service == "fseof":
        return dict(fseof(model, SUCC, n_steps=3).metadata)
    if service == "fvseof":
        return dict(fvseof(model, SUCC, n_steps=3).metadata)
    if service == "flux_response":
        return dict(flux_response(model, "EX_o2_e", n_steps=5).metadata)

    weights = gene_to_reaction_weights(model, _expression(model))
    if service == "eflux2":
        return dict(eflux2(model, weights).metadata)
    if service == "lad":
        return dict(lad(model, weights).metadata)
    if service == "predict_condition_fluxes":
        expression = _expression(model)
        table = pd.DataFrame(
            {
                "cond_a": pd.Series(expression),
                "cond_b": pd.Series({g: v * 2 for g, v in expression.items()}),
            }
        )
        return dict(predict_condition_fluxes(model, table, method="lad").metadata)

    if service in {"transformation_targets", "revert_targets"}:
        return _directional_metadata(service)

    if service in {
        "moma",
        "room",
        "knockout_comparison_moma",
        "knockout_comparison_room",
        "batch_comparison",
    }:
        model = _comparison_model()
    reference = reference_state_pfba(model)
    if service == "random_flux_sampling":
        return dict(
            random_flux_sampling(model, n=50, thinning=10, processes=1, seed=7).metadata
        )
    if service == "reference_constrained_sampling":
        return dict(
            reference_constrained_sampling(
                model, reference, n=50, thinning=10, processes=1, seed=7
            ).metadata
        )
    if service == "reference_flux":
        return dict(reference_flux(model, "pfba").metadata)

    flux_reference = reference_flux(model, "pfba")
    if service == "moma":
        return dict(moma(model, flux_reference).metadata)
    if service == "room":
        return dict(room(model, flux_reference).metadata)
    if service == "knockout_comparison_moma":
        return dict(
            knockout_comparison(
                model, flux_reference, ["R2"], method="moma_l2"
            ).metadata
        )
    if service == "knockout_comparison_room":
        return dict(
            knockout_comparison(model, flux_reference, ["R2"], method="room").metadata
        )
    if service == "batch_comparison":
        perturbations = reaction_perturbations(
            model, [r.id for r in model.reactions][:3]
        )
        return dict(batch_comparison(model, flux_reference, perturbations).metadata)

    # OptKnock/RobustKnock delegate the bilevel MILP to the optional ``straindesign``
    # package, as in ``tests/test_strain_design.py``.
    pytest.importorskip("straindesign")
    from cmm.features import optknock, robustknock

    search = optknock if service == "optknock" else robustknock
    model = _strain_design_model()
    return dict(search(model, "PRODUCT", max_knockouts=1, max_solutions=1).metadata)


@pytest.mark.parametrize("service", SERVICE_CASES)
def test_every_public_service_carries_the_provenance_block(service, core):
    metadata = _metadata(service, core)

    absent = [field for field in REQUIRED_FIELDS if field not in metadata]
    assert not absent, f"{service} is missing provenance fields {absent}"

    # Populated with real values, not placeholders standing in for one.
    assert ISO_UTC.match(str(metadata["timestamp_utc"])), metadata["timestamp_utc"]
    assert metadata["solver"], service
    assert metadata["solver_version"] not in (None, "", "unknown", "not-installed")
    assert metadata["platform"], service
    assert len(str(metadata["model_sha256"])) == 64
    assert isinstance(metadata["parameters"], dict)

    # ``seed`` is recorded, never invented: sampling reports the requested seed and the two
    # strain-design searches report their explicit reproducible backend default.  An unseeded
    # deterministic method reports null rather than a made-up value.
    if service in ("random_flux_sampling", "reference_constrained_sampling"):
        assert metadata["seed"] == 7
    elif service in ("optknock", "robustknock"):
        assert metadata["seed"] == 0
        assert metadata["parameters"]["seed"] == 0
        assert metadata["parameters"]["strain_design_seed"] == 0
    else:
        assert metadata["seed"] is None, f"{service} invented a seed"


def test_the_registry_covers_every_shipped_analysis_feature():
    """A shipped feature with no entry in :data:`SERVICES` is an untested provenance surface."""

    from cmm.features import INCLUDED_FEATURES

    # Features that are not a single numeric service: the GUI slider, and the two family
    # names whose services are registered individually above.
    not_a_numeric_service = {
        "flux_visualization_slider",
        "omics_integration",  # -> eflux2 / lad / predict_condition_fluxes
        "batch_moma_room",  # -> batch_comparison
        "revert_metabolism",  # -> revert_targets
        "flux_response_analysis",  # -> flux_response
        "production_target_workflow",  # composed, tested orchestration
        # Also composed orchestration: every number it reports comes from a service
        # registered above. Its own run-level block is asserted against REQUIRED_FIELDS
        # by test_jev.py::test_the_run_provenance_carries_every_required_field.
        "jev_target_design",
        "publication_reporting",  # artifact validation and rendering, not a numeric solve
    }
    registered = set(SERVICES) | {
        "knockout_comparison",
        "eflux2",
        "lad",
        "predict_condition_fluxes",
        "batch_comparison",
        "revert_targets",
        "flux_response",
    }
    uncovered = {
        feature
        for feature in INCLUDED_FEATURES
        if feature not in not_a_numeric_service and feature not in registered
    }
    assert not uncovered, (
        f"shipped features with no provenance test: {sorted(uncovered)}"
    )


def test_feature_manifest_separates_the_documented_planned_services():
    """The public roadmap and machine-readable feature manifest must name the same backlog."""

    from cmm.features import INCLUDED_FEATURES, PLANNED_FEATURES

    assert PLANNED_FEATURES == (
        "dynamic_fba",
        "enzyme_constrained_modeling",
        "fvseof_grouping_constraints",
    )
    assert set(INCLUDED_FEATURES).isdisjoint(PLANNED_FEATURES)
