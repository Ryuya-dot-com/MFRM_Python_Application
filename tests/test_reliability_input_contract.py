"""Reliability input and saved-output contracts; not sampling validation."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def measures(*, eap=True, estimates=(-1., 0., 1.)):
    return pd.DataFrame({
        "Facet": ["Person" if eap else "Rater"] * len(estimates),
        "Estimate": estimates, "SE": [.2] * len(estimates),
        "Infit": [1.] * len(estimates), "Outfit": [1.] * len(estimates),
        "SE_Method": ["MML latent-ability posterior SD from quadrature" if eap else
                      "conditional information approximation: other fitted parameters held fixed"] * len(estimates),
        "SE_Status": ["posterior_eap_sd" if eap else "conditional_approximation"] * len(estimates),
    })


@pytest.mark.parametrize("column,value", [
    ("SE", np.nan), ("SE", np.inf), ("SE", -.2), ("SE", 0.),
    ("Estimate", np.nan), ("Estimate", np.inf),
])
def test_invalid_pair_never_inflates_reliability(column, value):
    source = measures(estimates=(-1., 0., 1., 100.))
    source.loc[3, column] = value
    before = source.copy(deep=True)
    row = app.calc_reliability(source).iloc[0]
    assert pd.isna(row.Reliability) and pd.isna(row.Separation) and pd.isna(row.Strata)
    assert row.ReliabilityStatus == "unassessable"
    assert row.ReliabilityReasonCode == "stat.reliability.invalid_measure_pairs"
    assert row.ValidLevels == 3 and row.Levels == 4
    if column == "SE":
        assert row.SD == pytest.approx(source.Estimate.std(ddof=1))
    pd.testing.assert_frame_equal(source, before)


def test_zero_eap_variance_is_zero_empirical_reliability():
    row = app.calc_reliability(measures(estimates=(0., 0., 0.))).iloc[0]
    assert row.Reliability == 0. and row.Separation == 0.
    assert row.Strata == pytest.approx(1 / 3)
    assert row.ReliabilityStatus == "diagnostic"
    assert row.RMSE == pytest.approx(.2)


def test_zero_se_based_variance_remains_unassessable():
    row = app.calc_reliability(measures(eap=False, estimates=(0., 0., 0.))).iloc[0]
    assert pd.isna(row.Reliability)
    assert row.ReliabilityReasonCode == "stat.reliability.zero_observed_variance"


def test_single_valid_row_does_not_establish_reliability():
    row = app.calc_reliability(measures(estimates=(1.,))).iloc[0]
    assert pd.isna(row.Reliability)
    assert row.ReliabilityReasonCode == "stat.reliability.insufficient_levels"


@pytest.mark.parametrize("eap,expected", [(True, 1 / 1.04), (False, .96)])
def test_positive_control_and_calibration_hold_keep_formula(eap, expected):
    source = measures(eap=eap)
    before = app.calc_reliability(source).iloc[0]
    source["SE_Status"] = "calibration_inference_withheld"
    after = app.calc_reliability(source).iloc[0]
    assert before.Reliability == pytest.approx(expected)
    assert after.Reliability == pytest.approx(expected)
    assert after.ReliabilityStatus == "withheld"


def test_mixed_uncertainty_bases_do_not_select_eap_by_one_row():
    source = measures()
    source.loc[2, "SE_Method"] = "conditional information approximation"
    source.loc[2, "SE_Status"] = "conditional_approximation"
    row = app.calc_reliability(source).iloc[0]
    assert pd.isna(row.Reliability)
    assert row.ReliabilityReasonCode == "stat.reliability.mixed_uncertainty_basis"


def test_legacy_status_can_preserve_technical_eap_value_but_cannot_authorize_public_use():
    source = measures().drop(columns="SE_Method")
    calculated = app.calc_reliability(source)
    assert calculated.Reliability.iloc[0] == pytest.approx(1 / 1.04)
    table, _ = app._publication_reliability_table({"opt": SimpleNamespace(success=True)}, {"measures": source, "reliability": calculated})
    assert table.Use.iloc[0] == "Withheld" and table.Reliability.iloc[0] == ""


@pytest.mark.parametrize("case", ["missing_se", "stale_numbers", "legacy_no_raw"])
def test_saved_reliability_is_requalified_without_changing_saved_numbers(case):
    raw = measures(eap=False)
    saved = app.calc_reliability(raw)
    if case == "missing_se":
        raw.loc[2, "SE"] = np.nan
    elif case == "stale_numbers":
        saved.loc[0, "Reliability"] = .9999
    else:
        saved = saved[["Facet", "Levels", "SD", "RMSE", "Separation", "Strata", "Reliability", "SE_StatusSummary"]]
    before = saved.copy(deep=True)
    diagnostics = {"reliability": saved}
    if case != "legacy_no_raw":
        diagnostics["measures"] = raw
    table, note = app._publication_reliability_table({"opt": SimpleNamespace(success=True)}, diagnostics)
    assert table.Use.iloc[0] == "Withheld"
    assert table[["RMSE", "Separation", "Strata", "Reliability"]].eq("").all().all()
    assert table.SD.iloc[0] == round(saved.SD.iloc[0], 3)
    assert "withheld" in note.lower()
    audit = app.build_statistical_assumption_audit({"opt": SimpleNamespace(success=True)}, diagnostics).set_index("Area")
    assert audit.loc["Reliability / separation uncertainty", "Status"] == "WITHHELD"
    assert "stat.reliability." in audit.loc["Reliability / separation uncertainty", "Evidence"]
    pd.testing.assert_frame_equal(saved, before)


def test_valid_known_basis_keeps_public_values_and_review_audit():
    raw = measures(eap=False)
    saved = app.calc_reliability(raw)
    diagnostics = {"measures": raw, "reliability": saved}
    table, _ = app._publication_reliability_table({"opt": SimpleNamespace(success=True)}, diagnostics)
    assert table.Use.iloc[0] == "Diagnostic" and table.Reliability.iloc[0] == .96
    audit = app.build_statistical_assumption_audit({"opt": SimpleNamespace(success=True)}, diagnostics).set_index("Area")
    assert audit.loc["Reliability / separation uncertainty", "Status"] == "Review"


def test_missing_facet_in_current_measures_does_not_fall_back_to_saved_ready():
    saved = app.calc_reliability(measures(eap=False))
    diagnostics = {"measures": measures(), "reliability": saved}
    table, _ = app._publication_reliability_table({"opt": SimpleNamespace(success=True)}, diagnostics)
    assert table.Use.iloc[0] == "Withheld"


def test_per_facet_qualification_preserves_valid_neighbor_and_raw_saved_values():
    raw = pd.concat([measures(), measures(eap=False)], ignore_index=True)
    saved = app.calc_reliability(raw)
    raw.loc[0, "SE"] = np.nan
    qualified = app.qualified_reliability_table({"opt": SimpleNamespace(success=True)}, {"measures": raw, "reliability": saved})
    pd.testing.assert_frame_equal(qualified[saved.columns.drop(["ReliabilityStatus", "ReliabilityReasonCode"])],
                                  saved.drop(columns=["ReliabilityStatus", "ReliabilityReasonCode"]))
    table, _ = app._publication_reliability_table({"opt": SimpleNamespace(success=True)}, {"measures": raw, "reliability": saved})
    table = table.set_index("Facet")
    assert table.loc["Person", "Use"] == "Withheld"
    assert table.loc["Person", "Reliability"] == ""
    assert table.loc["Rater", "Use"] == "Diagnostic"
    assert table.loc["Rater", "Reliability"] == .96


def test_saved_current_provenance_can_be_checked_without_original_measures():
    saved = app.calc_reliability(measures(eap=False))
    table, _ = app._publication_reliability_table({"opt": SimpleNamespace(success=True)}, {"reliability": saved})
    assert table.Use.iloc[0] == "Diagnostic" and table.Reliability.iloc[0] == .96
