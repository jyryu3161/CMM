from __future__ import annotations

import math

import pandas as pd
import pytest
from cmm.omics.conditions import (
    flux_log_change,
    predict_condition_fluxes,
    read_expression_table,
    sign_flips,
)


def _two_condition_table(ecoli_core):
    genes = [g.id for g in ecoli_core.genes]
    return pd.DataFrame(
        {
            "condA": [20.0] * len(genes),
            "condB": [50.0 if i % 2 == 0 else 5.0 for i in range(len(genes))],
        },
        index=genes,
    )


def test_predict_condition_fluxes_eflux2(ecoli_core):
    table = _two_condition_table(ecoli_core)
    cf = predict_condition_fluxes(ecoli_core, table, method="eflux2")
    assert cf.conditions() == ("condA", "condB")
    assert all(r.status == "optimal" for r in cf.results.values())
    assert len(cf.fluxes("condA")) == len(ecoli_core.reactions)


def test_predict_condition_fluxes_lad(ecoli_core):
    table = _two_condition_table(ecoli_core)
    cf = predict_condition_fluxes(ecoli_core, table, method="lad")
    assert all(r.status == "optimal" for r in cf.results.values())


def test_flux_log_change_magnitude():
    a = {"R1": 1.0, "R2": 4.0, "R3": 0.0}
    b = {"R1": 2.0, "R2": 1.0, "R3": 3.0}
    # pseudocount=0 is exact for nonzero-flux reactions.
    lc = flux_log_change(a, b, reactions=["R1", "R2"], pseudocount=0.0)
    assert lc["R1"] == pytest.approx(1.0)  # log2(2/1)
    assert lc["R2"] == pytest.approx(-2.0)  # log2(1/4)
    # R3 turns on; the default pseudocount makes it a large finite positive change.
    lc_eps = flux_log_change(a, b, pseudocount=1e-3)
    assert lc_eps["R3"] > 5


def test_flux_log_change_handles_zero_without_pseudocount_divzero():
    # default pseudocount avoids division by zero on an off->off reaction.
    lc = flux_log_change({"R": 0.0}, {"R": 0.0})
    assert lc["R"] == pytest.approx(0.0)


def test_flux_log_change_zero_pseudocount_is_guarded():
    # Even with pseudocount=0 (no smoothing), an all-zero or off->on reaction must not raise.
    lc = flux_log_change({"R1": 0.0, "R2": 0.0}, {"R1": 0.0, "R2": 5.0}, pseudocount=0.0)
    assert lc["R1"] == 0.0
    assert lc["R2"] == math.inf


def test_flux_log_change_switch_off_is_negative_inf():
    # A reaction that switches OFF (nonzero source -> zero target) with pseudocount=0 would
    # hit log2(0/a); it must return -inf symmetrically with the switch-on (+inf) case.
    lc = flux_log_change({"R": 7.0}, {"R": 0.0}, pseudocount=0.0)
    assert lc["R"] == -math.inf


def test_sign_flips_detects_reversal():
    a = {"R1": 5.0, "R2": -3.0, "R3": 2.0}
    b = {"R1": -5.0, "R2": -3.0, "R3": 4.0}
    assert sign_flips(a, b) == ["R1"]


def test_read_expression_table(tmp_path):
    path = tmp_path / "expr.csv"
    path.write_text("Gene ID,c1,c2\ng1,10,20\ng2,5,1\n")
    table = read_expression_table(str(path))
    assert list(table.columns) == ["c1", "c2"]
    assert table.loc["g1", "c2"] == 20.0
    assert math.isclose(table.loc["g2", "c1"], 5.0)


def test_read_expression_table_detects_tsv(tmp_path):
    path = tmp_path / "expr.tsv"
    path.write_text("gene\tc1\tc2\ng1\t10\t20\ng2\t5\t1\n")

    table = read_expression_table(str(path))

    assert list(table.columns) == ["c1", "c2"]
    assert table.loc["g1", "c2"] == 20.0
