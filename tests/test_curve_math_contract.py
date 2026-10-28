"""Curve interpretation, explicit ZSTD formulas and the visible saved result."""
import json
import re

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import engine, visual_maps
from tests.test_figure_workflow import _view, _saved_plot
from tests.test_publication_figures import synthetic_result, synthetic_diagnostics
from tests.test_facets_df_zstd_alignment import small_rsm_jmle_diagnostics


@pytest.mark.parametrize("model", ["PCM", "GPCM"])
def test_mean_parameter_curve_is_not_labeled_as_average_probability(model):
    result = {"params": {"steps_mat": np.array([[-2., 2.], [-.5, .5]]), "slopes": np.array([.5, 2.])},
              "config": {"model": model, "n_cat": 3, "step_facet": "Task", "slope_facet": "Task"},
              "prep": {"rating_min": 0, "levels": {"Task": ["T1", "T2"]}}}
    curves = [app.build_category_probability_curve_data(result, i, theta_grid=np.array([0.])) for i in (None, 0, 1)]
    probabilities = [c["probability_wide"].drop(columns="Theta").to_numpy()[0] for c in curves]
    assert abs(probabilities[0][1] - np.mean(probabilities[1:], axis=0)[1]) > .01
    assert curves[0]["metadata"]["CurveBasis"] == "mean_parameters"
    assert "Hypothetical" in curves[0]["metadata"]["Scope"]
    exported = app.category_probability_curve_export_table(result)
    assert exported.loc[exported.CurveBasis.eq("mean_parameters"), "Scope"].str.contains("Hypothetical").all()
    information = app.information_curve_export_table(result)
    assert information.InformationBasis.eq("single_response_conditional").all()


@pytest.mark.parametrize("value,df", [(.975, 20.), (1.8, .5), (.2, 400.)])
def test_named_transforms_match_facets_formulas_without_changing_legacy(value, df):
    q = np.sqrt(2/df)
    no_offset = (value**(1/3)-1)*3/q
    expected = {"wilson_hilferty": no_offset+q/3, "wilson_hilferty_no_offset": no_offset,
                "linear_normal": (value-1)/q}
    for name, ref in expected.items():
        actual = engine.zstd_from_mnsq_facets(value, df, zstd_transform=name, cap=None)
        assert actual == pytest.approx(ref, abs=1e-13)
    for legacy, ref in ((False, "wilson_hilferty"), (True, "linear_normal")):
        assert engine.zstd_from_mnsq_facets(value, df, whexact=legacy, cap=None) == pytest.approx(expected[ref])
    with pytest.raises(ValueError):
        engine.zstd_from_mnsq_facets(value, df, whexact=True, zstd_transform="wilson_hilferty")
    with pytest.raises(ValueError):
        engine.zstd_from_mnsq_facets(value, df, zstd_transform="unrecognized")


def test_transform_reaches_diagnostics_recalculation_and_saved_metadata(small_rsm_jmle_diagnostics):
    result = small_rsm_jmle_diagnostics
    diag = engine.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[],
                                   zstd_transform="wilson_hilferty_no_offset", facets_zstd_cap=3.)
    assert diag["fit_zstd_transform_id"] == "wilson_hilferty_no_offset"
    assert diag["fit"].FitZSTDTransform.eq("Wilson-Hilferty without offset").all()
    for method in ("facets", "engine", "both"):
        overall, fit = app._recompute_fit_with_method(diag, method, facets_zstd_cap=3.)
        assert overall.FitZSTDTransform.eq("Wilson-Hilferty without offset").all()
        assert fit.FitZSTDTransform.eq("Wilson-Hilferty without offset").all()
        note = visual_maps.fit_convention_note({"fit": fit}, lang="en")
        assert "without offset" in note
        assert ("capped at ±3" in note) == (method == "facets")
        if method != "facets":
            assert "uncapped" in note


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_level_first_curves_and_human_readable_saved_figures(synthetic_result, synthetic_diagnostics, lang):
    result = synthetic_result
    result["config"].update(model="GPCM", method="MML", step_facet="Task", slope_facet="Task", n_cat=3)
    result["params"].update(steps_mat=np.array([[-2., 2.], [-.5, .5]]), slopes=np.array([.5, 2.]))
    result["prep"]["levels"] = {"Task": ["T1", "T2"]}
    at = AppTest.from_function(_view, args=(result, synthetic_diagnostics, lang)).run()
    at.selectbox(key="guided_figures_panel").select("category_probability").run()
    assert not at.exception
    assert at.selectbox(key="category_curve_scope_GPCM_Task").value == 1
    context, layout, _, html = _saved_plot(at)
    assert "Task = T1" in layout["title"]["text"]
    visible = re.sub(r"<script\b.*?</script>", "", html, flags=re.S)
    assert not any(term in visible for term in ("TECHNICAL_ONLY", "SourceOrder", "unqualified", "terminal-point", "FACETS-style"))
    assert context["calibration_hold_reason"]  # Machine-readable evidence remains.
    at.selectbox(key="category_curve_scope_GPCM_Task").select(0).run()
    _, layout, *_ = _saved_plot(at)
    assert ("平均パラメータ" if lang == "ja" else "mean parameters") in layout["title"]["text"]
    assert any(("平均ではありません" if lang == "ja" else "not the mean") in c.value for c in at.caption)
    at.selectbox(key="guided_figures_panel").select("information_curves").run()
    assert not at.exception
    assert at.selectbox(key="information_curve_scope_GPCM_Task").value == 1
    _, layout, *_ = _saved_plot(at)
    assert ("1評定" if lang == "ja" else "Single-response") in layout["title"]["text"]
    assert ("1評定" if lang == "ja" else "Single-response") in layout["yaxis2"]["title"]["text"]


def test_unknown_saved_zstd_settings_are_not_invented():
    note = visual_maps.fit_convention_note({}, lang="en")
    assert "not recorded" in note
    assert "±9" not in note


def test_category_review_uses_levels_even_when_mean_parameter_curve_hides_problem(monkeypatch):
    result = {"params": {"steps_mat": np.array([[2., -2.], [-6., 6.]])},
              "config": {"model": "PCM", "n_cat": 3, "step_facet": "Task"},
              "prep": {"rating_min": 0, "levels": {"Task": ["T1", "T2"]}}}
    # Keep observed evidence neutral to isolate the curve's contribution to the decision.
    monkeypatch.setattr(app, "calc_category_stats", lambda *a, **kw: pd.DataFrame({
        "Category": [0, 1, 2], "Count": [20, 20, 20], "InfitFlag": False,
        "OutfitFlag": False, "FitUnavailable": False, "ZSTDFlag": False}))
    diag = {"obs": pd.DataFrame({"Observed": [0, 1, 2]})}
    reference = app.category_probability_curve_diagnostic_table(app.build_category_probability_curve_data(result))
    assert reference.loc[reference.CategoryValue.eq(1), "IsEverMostProbable"].all()
    evidence = app.rating_scale_category_evidence_table(result, diag).set_index("Category")
    assert evidence.loc[1, "CategoryEvidenceStatus"] == "Review"
    assert evidence.loc[1, "CurveNoDominanceScopes"] == "Task = T1"
    assert evidence.loc[1, "CurveMinimumPeakProbability"] < .35
    assert "Hypothetical" not in evidence.to_string()
    assert "CurvePeakTheta" not in evidence.columns  # No invented pooled peak coordinate.
    dashboard = app.rating_scale_functioning_dashboard(result, diag).set_index("Check")
    assert dashboard.loc["Category characteristic curves", "Status"] == "Review"
    assert "Hypothetical" not in dashboard.to_string()
    candidates = app.rating_scale_recode_candidate_table(result, diag)
    assert candidates.CandidateID.str.startswith("C").any()
    assert "Task = T1" in candidates.to_string()
    assert "Hypothetical" not in candidates.to_string()
