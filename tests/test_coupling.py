from __future__ import annotations

import pytest
from cobra import Metabolite, Model, Reaction

from cmm.features import CoupledSets, coupled_reaction_sets


def _linear_chain() -> Model:
    """A -> B -> C -> D with an inlet and an outlet.

    Every internal reaction carries the same flux in every steady state, so the whole chain is
    one fully coupled set: cutting any one of them stops all of them.
    """

    model = Model("chain")
    mets = {name: Metabolite(f"{name}_c", compartment="c") for name in "abcd"}
    specs = [
        ("IN", {mets["a"]: 1}),
        ("R1", {mets["a"]: -1, mets["b"]: 1}),
        ("R2", {mets["b"]: -1, mets["c"]: 1}),
        ("R3", {mets["c"]: -1, mets["d"]: 1}),
        ("OUT", {mets["d"]: -1}),
    ]
    for rid, stoich in specs:
        rxn = Reaction(rid)
        rxn.bounds = (0.0, 100.0)
        rxn.add_metabolites(stoich)
        model.add_reactions([rxn])
    return model


def _branched() -> Model:
    """A -> B by two parallel routes, so the branches are not coupled to each other."""

    model = Model("branched")
    mets = {name: Metabolite(f"{name}_c", compartment="c") for name in "ab"}
    specs = [
        ("IN", {mets["a"]: 1}),
        ("LEFT", {mets["a"]: -1, mets["b"]: 1}),
        ("RIGHT", {mets["a"]: -1, mets["b"]: 1}),
        ("OUT", {mets["b"]: -1}),
    ]
    for rid, stoich in specs:
        rxn = Reaction(rid)
        rxn.bounds = (0.0, 100.0)
        rxn.add_metabolites(stoich)
        model.add_reactions([rxn])
    return model


def test_linear_chain_is_one_set():
    sets = coupled_reaction_sets(_linear_chain())
    assert len(sets) == 1
    assert set(sets.members("R2")) == {"IN", "R1", "R2", "R3", "OUT"}
    # Alphabetically first, so a rerun picks the same stand-in.
    assert sets.representatives == ("IN",)
    assert sets.representative_for("R3") == "IN"


def test_branches_are_not_coupled_to_each_other():
    sets = coupled_reaction_sets(_branched())
    # IN and OUT still carry the total and are coupled; the two branches split it freely.
    assert sets.groups["LEFT"] != sets.groups["RIGHT"]
    assert sets.groups["IN"] == sets.groups["OUT"]
    assert len(sets) == 3


def test_candidate_selection_preserves_the_unselected_network():
    model = _branched()
    # RIGHT is not a candidate, but still carries flux: LEFT is not forced to equal IN.
    sets = coupled_reaction_sets(model, ["IN", "LEFT", "OUT"])
    assert len(sets) == 2
    assert sets.members("LEFT") == ("LEFT",)
    assert sets.members("IN") == ("IN", "OUT")


def test_parallel_candidates_are_not_merged_after_essential_reactions_are_filtered():
    model = _branched()
    model.objective = "OUT"
    for left, right in ((10, 0), (0, 10)):
        with model:
            model.reactions.LEFT.bounds = (left, left)
            model.reactions.RIGHT.bounds = (right, right)
            assert model.slim_optimize() == pytest.approx(10)
    sets = coupled_reaction_sets(model, ["LEFT", "RIGHT"])
    assert sets.representatives == ("LEFT", "RIGHT")


def test_opposite_sign_ratio_groups_together():
    """v_i = -c * v_j is the same intervention as v_i = c * v_j."""

    model = Model("signed")
    a = Metabolite("a_c", compartment="c")
    forward = Reaction("FWD")
    forward.bounds = (0.0, 100.0)
    forward.add_metabolites({a: 1})
    consume = Reaction("USE")
    consume.bounds = (-100.0, 0.0)  # written in the reverse direction
    consume.add_metabolites({a: 1})
    model.add_reactions([forward, consume])
    sets = coupled_reaction_sets(model)
    assert len(sets) == 1


def test_reactions_with_no_steady_state_flux_are_separate_singletons():
    model = _linear_chain()
    dead = Reaction("DEAD")
    dead.bounds = (0.0, 100.0)
    dead.add_metabolites({Metabolite("x_c", compartment="c"): -1})
    other = Reaction("DEAD2")
    other.bounds = (0.0, 100.0)
    other.add_metabolites({Metabolite("y_c", compartment="c"): -1})
    model.add_reactions([dead, other])
    sets = coupled_reaction_sets(model)
    # They carry no flux, so they are coupled to nothing — including to each other.
    assert sets.groups["DEAD"] != sets.groups["DEAD2"]
    assert {"DEAD", "DEAD2"} <= set(sets.representatives)
    assert sets.metadata["n_reactions_without_steady_state_flux"] == 2


def test_result_is_reproducible():
    model = _linear_chain()
    assert coupled_reaction_sets(model) == coupled_reaction_sets(model)


def test_provenance_reports_the_definition_used():
    sets = coupled_reaction_sets(_linear_chain())
    provenance = sets.to_provenance()
    # A ranking's denominator depends on this, so the definition rides with the numbers.
    assert provenance["coupling"] == "full"
    assert provenance["n_sets"] == 1
    assert provenance["n_reactions"] == 5
    assert provenance["largest_set"] == 5


def test_empty_selection_is_an_empty_result_not_an_error():
    sets = coupled_reaction_sets(_linear_chain(), [])
    assert len(sets) == 0
    assert isinstance(sets, CoupledSets)


def test_unknown_reaction_is_rejected():
    with pytest.raises(KeyError, match="absent from the model"):
        coupled_reaction_sets(_linear_chain(), ["R1", "NOPE"])


def test_non_positive_tolerance_is_rejected():
    with pytest.raises(ValueError, match="tolerance"):
        coupled_reaction_sets(_linear_chain(), tolerance=0.0)


@pytest.mark.genome_scale
def test_matches_the_validation_run_on_e_coli_core(ecoli_core):
    """Sanity check on a real model: sets exist, and every reaction lands in exactly one."""

    sets = coupled_reaction_sets(ecoli_core)
    assert len(sets.groups) == len(ecoli_core.reactions)
    assert 0 < len(sets) <= len(ecoli_core.reactions)
    for rid in list(sets.groups)[:20]:
        assert rid in sets.members(rid)
        assert sets.representative_for(rid) in sets.members(rid)
