from __future__ import annotations

import pytest
from cmm.core.flux_state import FluxState, reference_state_pfba
from cmm.core.solvers import SolverCapabilityError
from cmm.features.transformation import direction_from_states, transformation_targets


def _states(branched_model):
    source = reference_state_pfba(branched_model, name="A")  # routes through R2
    target = FluxState(
        {"SUP_A": 10, "R1": 10, "R2": 0, "R3": 10, "R5": 10, "BIOMASS": 10}, name="B"
    )
    return source, target


def test_direction_from_states():
    source = FluxState({"R1": 10.0, "R2": 0.0, "R3": 5.0})
    target = FluxState({"R1": 0.0, "R2": 8.0, "R3": 5.0})
    direction = direction_from_states(source, target)
    assert direction["R1"] == -1  # decreases
    assert direction["R2"] == 1  # increases
    assert direction["R3"] == 0  # unchanged


@pytest.mark.requires_qp
def test_moma_transformation_ranks_disease_branch_first(branched_model):
    source, target = _states(branched_model)
    ranking = transformation_targets(branched_model, source, target, method="moma")
    assert ranking.method == "transform_moma"
    assert ranking.best().target_id == "g2"
    assert ranking.best().score > 0  # genuinely moves flux toward the target
    scores = {t.target_id: t.score for t in ranking}
    assert scores["g3"] == pytest.approx(0, abs=1e-9)


@pytest.mark.requires_miqp
def test_mta_transformation_ranks_disease_branch_first(branched_model):
    source, target = _states(branched_model)
    ranking = transformation_targets(branched_model, source, target, method="mta")
    assert ranking.method == "transform_mta"
    assert ranking.best().target_id == "g2"


@pytest.mark.requires_qp
def test_transformation_does_not_mutate_model(branched_model):
    source, target = _states(branched_model)
    growth = branched_model.slim_optimize()
    transformation_targets(branched_model, source, target, method="moma")
    assert branched_model.slim_optimize() == pytest.approx(growth, abs=1e-6)


def test_moma_transformation_requires_qp(branched_model):
    branched_model.solver = "glpk"
    source, target = _states(branched_model)
    with pytest.raises(SolverCapabilityError) as exc:
        transformation_targets(branched_model, source, target, method="moma")
    assert exc.value.capability == "QP"


def test_unknown_method_raises(branched_model):
    source, target = _states(branched_model)
    with pytest.raises(ValueError, match="unknown method"):
        transformation_targets(branched_model, source, target, method="bogus")


def test_unknown_perturbation_raises(branched_model):
    source, target = _states(branched_model)
    with pytest.raises(ValueError, match="perturbation"):
        transformation_targets(
            branched_model, source, target, method="moma", perturbation="bogus"
        )


# --- provenance and labelling (round 2) ------------------------------------


@pytest.mark.requires_qp
@pytest.mark.requires_miqp
def test_both_paths_record_the_two_states_they_ran_between(branched_model):
    source, target = _states(branched_model)
    for method, formulation in (
        ("moma", "yizhak_2013_moma_baseline"),
        ("mta", "yizhak_2013_mta_with_flux_state_directions"),
    ):
        ranking = transformation_targets(branched_model, source, target, method=method)
        metadata = ranking.metadata
        # Without the target-state identity an `mta` run is indistinguishable from a
        # revert_targets run, and the `moma` path used to carry no provenance at all.
        assert metadata["source"] == "A"
        assert metadata["target"] == "B"
        assert metadata["transformation_method"] == method
        # The optimisation is Yizhak's; the direction set is not derived the published way,
        # so the tag inherited from revert_targets would overstate the correspondence.
        assert metadata["formulation"] == formulation


@pytest.mark.requires_qp
def test_moma_path_carries_run_provenance_and_tie_structure(branched_model):
    source, target = _states(branched_model)
    ranking = transformation_targets(branched_model, source, target, method="moma")
    metadata = ranking.metadata
    assert len(metadata["model_sha256"]) == 64
    assert metadata["parameters"]["method"] == "transformation_targets"
    assert metadata["n_perturbations"] == len(ranking)
    assert metadata["n_inert_dropped"] == 0
    assert metadata["n_distinct_scores"] >= 1
    assert metadata["largest_tie_block"] >= 1


def test_direction_from_states_labels_its_own_rule():
    source = FluxState({"R1": 10.0}, name="A")
    target = FluxState({"R1": 0.0}, name="B")
    direction = direction_from_states(source, target)
    assert direction.metadata["direction_rule"] == "flux_state_difference"
    assert direction.metadata["from"] == "A" and direction.metadata["to"] == "B"


@pytest.mark.requires_qp
def test_failed_moma_candidate_is_last_and_retains_status(parallel_pathway_model):
    model = parallel_pathway_model
    source = FluxState({"SUPPLY": 10, "R1": 10, "R2": 10, "R3": 0, "BIOMASS": 10})
    target = FluxState({"SUPPLY": 10, "R1": 10, "R2": 9, "R3": 1, "BIOMASS": 10})
    ranking = transformation_targets(model, source, target, targets=["g1", "g2"])
    frame = ranking.to_frame().set_index("target_id")
    assert ranking.best().target_id == "g2"
    assert frame.loc["g2", "score"] < 0
    assert frame.loc["g2", "status"] == "optimal"
    assert frame.loc["g1", "score"] == float("-inf")
    assert frame.loc["g1", "status"] == "infeasible"
    assert list(frame.index) == ["g2", "g1"]
    assert ranking.sorted(descending=False).targets[-1].target_id == "g1"
    assert ranking.metadata["n_nonoptimal"] == 1
    assert model.slim_optimize() == pytest.approx(10)

    failed = transformation_targets(model, source, target, targets=["g1"])
    assert failed.best() is None
    assert len(failed.to_frame()) == 1
