from __future__ import annotations

import pytest
from cobra import Metabolite, Model, Reaction


@pytest.fixture
def toy_model() -> Model:
    model = Model("toy")
    metabolite = Metabolite("a_c", compartment="c")

    source = Reaction("SOURCE_A")
    source.lower_bound = 0
    source.upper_bound = 10
    source.add_metabolites({metabolite: 1})

    biomass = Reaction("BIOMASS")
    biomass.lower_bound = 0
    biomass.upper_bound = 1000
    biomass.add_metabolites({metabolite: -1})

    product = Reaction("PRODUCT")
    product.lower_bound = 0
    product.upper_bound = 1000
    product.add_metabolites({metabolite: -1})

    model.add_reactions([source, biomass, product])
    model.objective = biomass
    model.objective_direction = "max"
    return model


@pytest.fixture
def branched_model() -> Model:
    """A two-branch model for perturbation / revert-metabolism tests.

    Topology (all metabolites in compartment c):

        SUP_A: -> A            (substrate supply, <= 10)
        R1:    A -> B          gene g1
        R2:    B -> P          gene g2   [short "disease" branch]
        R3:    B -> D          gene g3   [long "healthy" branch, step 1]
        R5:    D -> P          gene g5   [long "healthy" branch, step 2]
        BIOMASS: P ->          gene gb   (objective)

    Both branches yield the same biomass, but pFBA (minimal total flux) prefers
    the shorter R2 branch, so the source/disease reference routes flux through R2.
    Knocking out R2 forces the R3->R5 branch, which is the deterministic
    normalization target the revert tests expect to rank first.
    """

    model = Model("branched")
    a = Metabolite("A_c", compartment="c")
    b = Metabolite("B_c", compartment="c")
    d = Metabolite("D_c", compartment="c")
    p = Metabolite("P_c", compartment="c")

    def make(rid, stoich, gene, lb=0.0, ub=1000.0):
        r = Reaction(rid)
        r.lower_bound = lb
        r.upper_bound = ub
        r.add_metabolites(stoich)
        if gene:
            r.gene_reaction_rule = gene
        return r

    sup = make("SUP_A", {a: 1}, "", ub=10.0)
    r1 = make("R1", {a: -1, b: 1}, "g1")
    r2 = make("R2", {b: -1, p: 1}, "g2")
    r3 = make("R3", {b: -1, d: 1}, "g3")
    r5 = make("R5", {d: -1, p: 1}, "g5")
    biomass = make("BIOMASS", {p: -1}, "gb")

    model.add_reactions([sup, r1, r2, r3, r5, biomass])
    model.objective = biomass
    model.objective_direction = "max"
    return model


@pytest.fixture
def parallel_pathway_model() -> Model:
    """Two independent routes; g1 also controls their common essential precursor supply."""

    model = Model("parallel_pathways")
    a, b, c = [Metabolite(f"{name}_c", compartment="c") for name in "abc"]
    for rid, stoich, rule, lower in (
        ("SUPPLY", {a: 1}, "", 0),
        ("R1", {a: -1, b: 1}, "g1", 0),
        ("R2", {b: -1, c: 1}, "g1 and g2", 0),
        ("R3", {b: -1, c: 1}, "g3", 0),
        ("BIOMASS", {c: -1}, "", 1),
    ):
        reaction = Reaction(rid)
        reaction.add_metabolites(stoich)
        reaction.bounds = (lower, 10)
        reaction.gene_reaction_rule = rule
        model.add_reactions([reaction])
    model.objective = "BIOMASS"
    return model


@pytest.fixture
def published_mta_model() -> Model:
    """Independent Python construction of the official COBRA Toolbox MTA test network."""

    model = Model("published_mta_toy")
    a = Metabolite("A_c", compartment="c")
    b = Metabolite("B_c", compartment="c")
    c = Metabolite("C_c", compartment="c")
    d = Metabolite("D_c", compartment="c")

    def reaction(rid, metabolites, rule="", lower=0.0, upper=20.0):
        item = Reaction(rid)
        item.bounds = (lower, upper)
        item.add_metabolites(metabolites)
        item.gene_reaction_rule = rule
        return item

    model.add_reactions(
        [
            reaction("r1", {a: 1}, "g1.1 and g6.1"),
            reaction("r2", {a: -1, b: 1}, "g2.1 or g2.2"),
            reaction("r3", {a: -1, c: 1}),
            reaction("r4", {c: -1, d: 1}, "g4.1", lower=-20.0),
            reaction("r5", {a: -1, d: 1}, "g5.1"),
            reaction("r6", {d: -1, b: 1}, "g6.1"),
            reaction("r7", {b: -1}),
        ]
    )
    model.objective = model.reactions.r7
    return model


@pytest.fixture
def ecoli_core():
    """The e_coli_core textbook model (95 reactions) bundled with cobra."""

    from cobra.io import load_model

    return load_model("textbook")
