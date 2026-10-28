"""Readable comparison exports retain source values, scope and selected conditions."""
from copy import deepcopy
from io import BytesIO

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_model_choice_facets import result, capture_refits
from tests.test_help_contextual_app_integration import _locale_value


def view():
    import streamlit as st
    import streamlit_app as app
    if st.session_state.get("show_comparison", True):
        app.render_model_choice_guidance(st.session_state["result"])
    else:
        st.write("Other panel")


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_readable_csv_keeps_precision_conditions_holds_and_no_internal_codes(monkeypatch, lang):
    calls = capture_refits(monkeypatch)
    refit = app._refit_candidate_model
    def controlled(original, candidate, **kwargs):
        out = refit(original, candidate, **kwargs)
        out["summary"]["LogLik"] = -120.123456789
        if candidate == "PCM":
            out["opt"].success = False
        return out
    monkeypatch.setattr(app, "_refit_candidate_model", controlled)
    fit = result()
    bundle = app.compute_model_choice_comparison(fit, step_facet="Task", slope_facet="Rater")
    original = bundle["comparison"].copy(deep=True)
    original_result = deepcopy(fit["summary"])
    calls.clear()
    exports = {}
    download = app.st.download_button
    def capture(label, data, **kwargs):
        exports[kwargs["file_name"]] = (data, kwargs)
        return download(label, data, **kwargs)
    monkeypatch.setattr(app.st, "download_button", capture)
    at = AppTest.from_function(view)
    at.session_state["result"], at.session_state["lang"] = fit, lang
    source = app._model_choice_cache_key(fit)
    at.session_state[f"model_choice_selection::{source}"] = {"step_facet": "Task", "slopes": {"Task": "Rater"}}
    cache = app._model_comparison_cache_key(fit, "Task", "Rater")
    at.session_state[f"model_choice_comparison::{cache}"] = bundle
    at.run()
    assert not at.exception and not calls
    payload, options = exports["mfrm_model_comparison.csv"]
    assert payload.startswith(b"\xef\xbb\xbf") and options["on_click"] == "ignore"
    report = pd.read_csv(BytesIO(payload), float_precision="round_trip")
    label = lambda key: _locale_value(lang, "report_tables." + key)
    column = lambda key: label("model_choice_columns." + key)
    assert len(report) == 3 and len(report.columns) == 15
    assert not {"ReasonCode", "RoadmapIssueID", "QualificationStatus", "OutputID", "_section", "p", "Decision"}.intersection(report.columns)
    assert report[column("StepFacet")].tolist() == [label("model_choice_common_steps"), "Task", "Task"]
    assert report[column("SlopeFacet")].tolist() == [label("model_choice_fixed_slopes"), label("model_choice_fixed_slopes"), "Rater"]
    assert report[label("model_choice_export_source")].eq(app.build_result_analysis_identity(fit).analysis_id).all()
    assert report[label("model_choice_export_method")].eq(fit["config"]["method"]).all()
    assert report[label("model_choice_export_limits")].eq(label("model_choice_recommendation_withheld_template")).all()
    for metric in ["LogLik", "AIC", "DeltaAIC", "BIC", "DeltaBIC"]:
        np.testing.assert_allclose(report[column(metric)], original[metric], rtol=0, atol=0, equal_nan=True)
    pcm = report.loc[report[column("Model")].eq("PCM")].iloc[0]
    assert pcm[column("FitStatus")] == label("model_choice_status.refit_termination_unconfirmed")
    assert np.isnan(pcm[column("DeltaAIC")])
    assert report.loc[1, column("LogLik")] != float(at.table[1].value.loc[1, column("LogLik")])
    assert at.table[0].value.loc[1, column("DeltaAIC")] == "—"
    assert at.table[0].value.loc[0, column("AIC")] == "256"
    technical, options = exports["mfrm_model_choice_comparison.csv"]
    technical = pd.read_csv(BytesIO(technical))
    assert options["on_click"] == "ignore"
    assert technical.ComparisonStepFacet.eq("Task").all() and technical.ComparisonSlopeFacet.eq("Rater").all()
    assert not technical.PublicConclusionAllowed.fillna(False).any()
    assert "ReasonCode" in technical and "_section" in technical
    pd.testing.assert_frame_equal(bundle["comparison"], original)
    pd.testing.assert_frame_equal(fit["summary"], original_result)


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_navigation_and_source_restore_keep_choices_and_results_without_refit(monkeypatch, lang):
    calls = capture_refits(monkeypatch)
    fit = result()
    at = AppTest.from_function(view)
    at.session_state["result"], at.session_state["lang"] = fit, lang
    at.run()
    at.selectbox[0].select("Task").run()
    at.selectbox[1].select("Rater").run()
    at.button[0].click().run()
    assert not at.exception and len(calls) == 2
    table = at.table[0].value.copy()
    at.session_state["show_comparison"] = False
    at.run()
    assert not at.selectbox and not at.table
    at.session_state["show_comparison"] = True
    at.run()
    assert not at.exception and len(calls) == 2
    assert [x.value for x in at.selectbox] == ["Task", "Rater"]
    pd.testing.assert_frame_equal(at.table[0].value, table)
    other = deepcopy(fit); other["summary"]["LogLik"] -= 1
    at.session_state["result"] = other
    at.run()
    assert not at.table and at.selectbox[0].value is None
    at.session_state["result"] = fit
    at.run()
    assert len(calls) == 2 and not at.exception
    assert [x.value for x in at.selectbox] == ["Task", "Rater"]
    pd.testing.assert_frame_equal(at.table[0].value, table)
