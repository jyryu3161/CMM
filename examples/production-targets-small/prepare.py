"""Write the small mixed-acid fermentation model used by the restricted-licence SC-01 example."""

from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import write_sbml_model

FORMULAS = {
    "glc__D_e": "C6H12O6",
    "g6p_c": "C6H11O9P",
    "pyr_c": "C3H3O3",
    "accoa_c": "C23H34N7O17P3S",
    "coa_c": "C21H32N7O16P3S",
    "for_c": "CHO2",
    "for_e": "CHO2",
    "co2_c": "CO2",
    "co2_e": "CO2",
    "oaa_c": "C4H2O5",
    "mal_c": "C4H4O5",
    "succ_c": "C4H4O4",
    "succ_e": "C4H4O4",
    "lac__D_c": "C3H5O3",
    "lac__D_e": "C3H5O3",
    "ac_c": "C2H3O2",
    "ac_e": "C2H3O2",
    "etoh_c": "C2H6O",
    "etoh_e": "C2H6O",
    "o2_e": "O2",
    "atp_c": "C10H12N5O13P3",
    "adp_c": "C10H12N5O10P2",
    "nadh_c": "C21H27N7O14P2",
    "nad_c": "C21H26N7O14P2",
}

# (id, stoichiometry, gene rule, lower bound). Lumped anaerobic mixed-acid fermentation:
# glycolysis yields pyruvate and NADH; lactate, ethanol, and the reductive succinate branch
# reoxidise NADH; acetate yields ATP. Biomass draws ATP and NADH, so blocking ethanol (ADHE)
# forces the NADH surplus through succinate without lowering growth.
REACTIONS = (
    ("EX_glc__D_e", {"glc__D_e": -1}, "", -10.0),
    ("EX_o2_e", {"o2_e": -1}, "", 0.0),
    ("EX_co2_e", {"co2_e": -1}, "", 0.0),
    ("EX_succ_e", {"succ_e": -1}, "", 0.0),
    ("EX_lac__D_e", {"lac__D_e": -1}, "", 0.0),
    ("EX_ac_e", {"ac_e": -1}, "", 0.0),
    ("EX_etoh_e", {"etoh_e": -1}, "", 0.0),
    ("EX_for_e", {"for_e": -1}, "", 0.0),
    ("GLCpts", {"glc__D_e": -1, "atp_c": -1, "g6p_c": 1, "adp_c": 1}, "b0001", 0.0),
    (
        "GLYC",
        {"g6p_c": -1, "adp_c": -3, "nad_c": -2, "pyr_c": 2, "atp_c": 3, "nadh_c": 2},
        "b0002",
        0.0,
    ),
    ("LDH_D", {"pyr_c": -1, "nadh_c": -1, "lac__D_c": 1, "nad_c": 1}, "b0003", 0.0),
    (
        "PFL",
        {"pyr_c": -1, "coa_c": -1, "accoa_c": 1, "for_c": 1},
        "b0004 or b0005",
        0.0,
    ),
    ("FHL", {"for_c": -1, "co2_c": 1}, "b0006", 0.0),
    (
        "ACKr",
        {"accoa_c": -1, "adp_c": -1, "ac_c": 1, "atp_c": 1, "coa_c": 1},
        "b0007",
        0.0,
    ),
    (
        "ADHE",
        {"accoa_c": -1, "nadh_c": -2, "etoh_c": 1, "coa_c": 1, "nad_c": 2},
        "b0008",
        0.0,
    ),
    ("PPC", {"pyr_c": -1, "co2_c": -1, "oaa_c": 1}, "b0009", 0.0),
    ("MDH", {"oaa_c": -1, "nadh_c": -1, "mal_c": 1, "nad_c": 1}, "b0010", 0.0),
    (
        "FRD",
        {"mal_c": -1, "nadh_c": -1, "succ_c": 1, "nad_c": 1},
        "b0011 and b0012",
        0.0,
    ),
    ("LACt", {"lac__D_c": -1, "lac__D_e": 1}, "b0014", 0.0),
    ("ACt", {"ac_c": -1, "ac_e": 1}, "b0015", 0.0),
    ("ETOHt", {"etoh_c": -1, "etoh_e": 1}, "", 0.0),
    ("FORt", {"for_c": -1, "for_e": 1}, "b0016", 0.0),
    ("SUCCt", {"succ_c": -1, "succ_e": 1}, "b0017", 0.0),
    ("CO2t", {"co2_c": -1, "co2_e": 1}, "", -1000.0),
    ("ATPM", {"atp_c": -1, "adp_c": 1}, "", 1.0),
    (
        "BIOMASS",
        {
            "pyr_c": -1,
            "accoa_c": -1,
            "oaa_c": -1,
            "atp_c": -20,
            "nadh_c": -3,
            "coa_c": 1,
            "adp_c": 20,
            "nad_c": 3,
        },
        "",
        0.0,
    ),
)


def build_model() -> Model:
    model = Model("mixed_acid_demo")
    metabolites = {
        mid: Metabolite(mid, formula=formula, compartment=mid[-1])
        for mid, formula in FORMULAS.items()
    }
    for rid, stoichiometry, rule, lower in REACTIONS:
        reaction = Reaction(rid, lower_bound=lower, upper_bound=1000.0)
        reaction.add_metabolites(
            {
                metabolites[mid]: coefficient
                for mid, coefficient in stoichiometry.items()
            }
        )
        reaction.gene_reaction_rule = rule
        model.add_reactions([reaction])
    model.objective = "BIOMASS"
    return model


def main() -> None:
    destination = Path(__file__).resolve().parent / "data" / "mixed_acid_demo.xml"
    if destination.exists():
        raise SystemExit(
            f"Input already exists; preserve or move it first: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    write_sbml_model(build_model(), str(destination))
    print(destination)


if __name__ == "__main__":
    main()
