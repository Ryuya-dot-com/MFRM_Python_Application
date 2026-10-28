"""Closed-form controls for local-scale direction and shared-effect cancellation."""

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def _cells(target="Rater", sign=-1, reverse=False):
    rows = []
    for context, bias, se in (("T1", .4, .3), ("T2", -.2, .4)):
        row = {"FacetA": target, "FacetA_Level": "X1", "FacetA_Measure": .8,
               "FacetA_SE": 2., "FacetA_Index": 1, "FacetA_Sign": sign,
               "FacetB": "Task", "FacetB_Level": context, "FacetB_Measure": 0.,
               "FacetB_SE": .2, "FacetB_Index": len(rows)+1, "FacetB_Sign": -1,
               "Bias Size": bias, "S.E.": se, "Obs-Exp Average": bias/2,
               "Observd Count": 10., "ObsN": 10}
        if reverse:
            for suffix in ("", "_Level", "_Measure", "_SE", "_Index", "_Sign"):
                a, b = "FacetA"+suffix, "FacetB"+suffix
                row[a], row[b] = row[b], row[a]
        rows.append(row)
    return pd.DataFrame(rows)


@pytest.mark.parametrize("target,sign", [("Person", 1), ("Rater", -1), ("Rater", 1)])
@pytest.mark.parametrize("reverse", [False, True])
def test_local_measure_preserves_eta_shift_and_shared_effect_cancels(target, sign, reverse):
    row = app.calc_bias_pairwise(_cells(target, sign, reverse), target, "Task").iloc[0]
    assert sign*(row["Local Measure1"]-.8) == pytest.approx(.4)
    assert sign*(row["Local Measure2"]-.8) == pytest.approx(-.2)
    assert row["Contrast"] == pytest.approx(.6/sign)
    # Known nuisance parameters: Var(b1)=.09, Var(b2)=.16, Cov(b1,b2)=0.
    # The target's SE=2 does not enter either conditional local variance.
    assert row["SE1"] == pytest.approx(.3)
    assert row["SE2"] == pytest.approx(.4)
    assert row["SE"] == pytest.approx(.5)
    assert row["Standardized Contrast"] == pytest.approx(1.2/sign)
    assert pd.isna(row["t"]) and pd.isna(row["d.f."]) and pd.isna(row["Prob."])
    assert "conditional" in row["SE Method"].lower()
    assert row["QualificationStatus"] == "WITHHELD"
    assert not bool(row["PublicConclusionAllowed"])
    swapped = app.calc_bias_pairwise(
        _cells(target, sign, reverse).iloc[::-1], target, "Task").iloc[0]
    assert swapped.Contrast == pytest.approx(-row.Contrast)
    assert swapped.SE == pytest.approx(row.SE)


@pytest.mark.parametrize("measure", [.8, 1e15, np.nan])
@pytest.mark.parametrize("target_se", [.01, 100., np.nan])
def test_contrast_does_not_depend_on_shared_main_measure_or_se(measure, target_se):
    table = _cells().assign(FacetA_Measure=measure, FacetA_SE=target_se)
    row = app.calc_bias_pairwise(table, "Rater", "Task").iloc[0]
    assert row.Contrast == pytest.approx(-.6, abs=1e-14)
    assert row.SE == pytest.approx(.5)


def test_missing_or_conflicting_sign_is_not_silently_assumed():
    legacy = _cells().drop(columns="FacetA_Sign")
    with pytest.raises(ValueError, match="target_sign"):
        app.calc_bias_pairwise(legacy, "Rater", "Task")
    explicit = app.calc_bias_pairwise(legacy, "Rater", "Task", target_sign=-1)
    assert explicit.iloc[0].Contrast == pytest.approx(-.6)
    with pytest.raises(ValueError, match="sign"):
        app.calc_bias_pairwise(_cells(), "Rater", "Task", target_sign=1)
    for sign in (0, 2, np.nan):
        with pytest.raises(ValueError, match="sign"):
            app.calc_bias_pairwise(legacy, "Rater", "Task", target_sign=sign)


def test_requested_pair_is_selected_without_fallback_to_another_context():
    assert app.calc_bias_pairwise(_cells(), "Rater", "Unknown").empty
    unrelated = _cells().assign(FacetA="Other", FacetB="Else")
    mixed = pd.concat([unrelated, _cells()], ignore_index=True)
    result = app.calc_bias_pairwise(mixed, "Rater", "Task")
    assert len(result) == 1
    assert result.iloc[0].Contrast == pytest.approx(-.6)


@pytest.mark.parametrize("invalid", ["duplicate", "measure", "sign"])
def test_ambiguous_cells_or_shared_main_effect_are_rejected(invalid):
    table = _cells()
    if invalid == "duplicate":
        table = pd.concat([table, table.iloc[[0]]], ignore_index=True)
    else:
        table.loc[1, "FacetA_Measure" if invalid == "measure" else "FacetA_Sign"] = 1.
    with pytest.raises(ValueError):
        app.calc_bias_pairwise(table, "Rater", "Task")


@pytest.mark.parametrize("cell_se", [-.3, 0., np.nan, np.inf])
def test_unavailable_cell_information_is_not_turned_into_contrast_precision(cell_se):
    table = _cells()
    table.loc[0, "S.E."] = cell_se
    row = app.calc_bias_pairwise(table, "Rater", "Task").iloc[0]
    assert pd.isna(row.SE)
    assert pd.isna(row["Standardized Contrast"])
    assert row["SE Status"] == "not_available"


@pytest.mark.parametrize("model", ["RSM", "PCM", "GPCM"])
@pytest.mark.parametrize("positive", [False, True])
def test_fitted_bias_cells_record_the_actual_model_direction(model, positive):
    result = app.mfrm_estimate(
        app._make_self_test_gradient_data(), "Person", ["Rater", "Task"], "Score",
        model=model, method="JMLE", maxit=2, positive_facets=["Rater"] if positive else [])
    diag = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    bias = app.estimate_bias_interaction(
        result, diag, "Rater", "Task", omit_extreme=False, max_iter=1)["table"]
    assert bias.FacetA_Sign.eq(1 if positive else -1).all()
    assert bias.FacetB_Sign.eq(-1).all()
    pairwise = app.calc_bias_pairwise(bias, "Rater", "Task")
    assert not pairwise.empty
    assert pairwise.QualificationStatus.eq("WITHHELD").all()
    assert pairwise["Prob."].isna().all()


def test_fixed_calibration_bias_fits_each_cell_once(monkeypatch):
    result = app.mfrm_estimate(app._make_self_test_gradient_data(), "Person", ["Rater", "Task"],
                               "Score", model="RSM", method="JMLE", maxit=2)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    fit_cell = app._fit_conditional_bias_cell
    calls = []
    def counted(*args, **kwargs):
        calls.append(1)
        return fit_cell(*args, **kwargs)
    monkeypatch.setattr(app, "_fit_conditional_bias_cell", counted)
    output = app.estimate_bias_interaction(result, diagnostics, "Rater", "Task",
                                           omit_extreme=False, max_iter=100, tol=0.)
    assert len(calls) == len(output["table"])
    assert output["calculation_mode"] == "single_pass_fixed_calibration"
    record = output["iteration"]
    assert len(record) == 1
    assert "Iteration" not in record.columns
    assert record.iloc[0]["AttemptedCells"] == len(calls)
    assert not record.iloc[0]["CalibrationUpdated"]
