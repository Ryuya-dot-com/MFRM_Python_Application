"""The comparison question reaches refits, cache ownership and downloaded records."""
from copy import deepcopy

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_model_choice_completion import fits
from tests.test_help_contextual_app_integration import _locale_value


def result(model="RSM"):
    fit = fits()[model]
    fit["config"]["facet_names"] = ["Rater", "Task", "Singleton"]
    fit["prep"]["data"] = pd.DataFrame({
        "Person": ["P1", "P2", "P3", "P4"], "Score": [0, 1, 2, 1],
        "Rater": ["R1", "R1", "R2", "R2"], "Task": ["T1", "T2", "T1", "T2"],
        "Singleton": ["only"] * 4,
    })
    if model != "RSM":
        fit["config"]["step_facet"] = "Task"
    if model == "GPCM":
        fit["config"]["slope_facet"] = "Rater"
    return fit


def capture_refits(monkeypatch):
    calls = []
    def refit(original, candidate, *, step_facet, slope_facet):
        calls.append((candidate, step_facet, slope_facet))
        fit = deepcopy(original)
        fit["config"].update(model=candidate,
            step_facet=None if candidate == "RSM" else step_facet,
            slope_facet=slope_facet if candidate == "GPCM" else None)
        fit["summary"]["Model"] = candidate
        return fit
    monkeypatch.setattr(app, "_refit_candidate_model", refit)
    return calls


@pytest.mark.parametrize("step,slope", [(None, None), ("Person", None), ("Singleton", None),
                                        ("missing", None), ("Task", "missing")])
def test_ambiguous_or_invalid_facets_stop_before_any_refit(monkeypatch, step, slope):
    calls = capture_refits(monkeypatch)
    fit = result()
    original = deepcopy(fit)
    bundle = app.compute_model_choice_comparison(fit, step_facet=step, slope_facet=slope)
    assert not bundle["available"] and "Select observed" in bundle["reason"]
    assert not calls
    assert fit["config"] == original["config"]
    pd.testing.assert_frame_equal(fit["summary"], original["summary"])


@pytest.mark.parametrize("model,step,slope,reused", [
    ("RSM", "Task", "Rater", True), ("PCM", "Task", "Rater", True),
    ("PCM", "Rater", "Rater", False), ("GPCM", "Task", "Rater", True),
    ("GPCM", "Task", "Task", False), ("GPCM", "Rater", "Rater", False),
])
def test_current_fit_is_reused_only_for_matching_structure(monkeypatch, model, step, slope, reused):
    calls = capture_refits(monkeypatch)
    fit = result(model)
    bundle = app.compute_model_choice_comparison(fit, step_facet=step, slope_facet=slope)
    assert bundle["available"]
    assert {c[0] for c in calls} == set(app._MODEL_CHOICE_CANDIDATES) - ({model} if reused else set())
    assert all(c[1:] == (step, slope) for c in calls)
    row = bundle["comparison"].set_index("Model").loc[model]
    assert row.FitStatus == ("current" if reused else "refit_ok")
    assert bundle["comparison_spec"] == {"step_facet": step, "slope_facet": slope}
    assert bundle["source_result_key"] == app._model_choice_cache_key(fit)
    assert bundle["comparison"].StepFacet.tolist() == ["", step, step]
    assert bundle["comparison"].SlopeFacet.tolist() == ["", "", slope]


def test_omitted_facets_preserve_existing_gpcm_specification(monkeypatch):
    calls = capture_refits(monkeypatch)
    bundle = app.compute_model_choice_comparison(result("GPCM"))
    assert bundle["comparison_spec"] == {"step_facet": "Task", "slope_facet": "Rater"}
    assert len(calls) == 2


@pytest.mark.parametrize("config", [None, [], "old"])
def test_malformed_configuration_is_refused_before_reading_facet_settings(monkeypatch, config):
    calls = capture_refits(monkeypatch)
    fit = result(); fit["config"] = config
    bundle = app.compute_model_choice_comparison(fit, step_facet="Task")
    assert not bundle["available"] and not calls


def test_selected_facets_reach_the_real_estimator(monkeypatch):
    fit = app.mfrm_estimate(app._make_self_test_gradient_data(), person_col="Person",
        facet_cols=["Rater", "Task"], score_col="Score", model="RSM", method="JMLE", maxit=40)
    original = deepcopy(fit["summary"])
    estimate = app.mfrm_estimate
    returned = []
    def capture(**kwargs):
        candidate = estimate(**kwargs)
        returned.append(candidate)
        return candidate
    monkeypatch.setattr(app, "mfrm_estimate", capture)
    for facet in ["Task", "Rater"]:
        bundle = app.compute_model_choice_comparison(fit, step_facet=facet)
        assert bundle["available"]
        assert [r["config"]["step_facet"] for r in returned[-2:]] == [facet, facet]
        assert returned[-1]["config"]["slope_facet"] == facet
        assert list(bundle["comparison"].FitStatus)[0].startswith("current")
    assert len(returned) == 4
    pd.testing.assert_frame_equal(fit["summary"], original)


def view():
    import streamlit as st
    import streamlit_app as app
    app.render_model_choice_guidance(st.session_state["result"])


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_selection_run_cache_and_export_form_one_workflow(monkeypatch, lang):
    calls = capture_refits(monkeypatch)
    exported = []
    to_csv = app.to_csv_bytes
    monkeypatch.setattr(app, "to_csv_bytes", lambda frame: exported.append(frame.copy()) or to_csv(frame))
    fit = result()
    at = AppTest.from_function(view)
    at.session_state["result"], at.session_state["lang"] = fit, lang
    at.run()
    assert not at.exception and not calls and not at.table
    assert at.selectbox[0].value is None and at.button[0].disabled
    assert at.selectbox[0].options == ["Rater", "Task"]
    at.selectbox[0].select("Task").run()
    assert not calls and not at.table and not at.button[0].disabled
    assert at.selectbox[1].value == "Task"
    at.button[0].click().run()
    assert not at.exception and len(calls) == 2 and len(at.table) == 2
    assert len(at.warning) == 1 and all(not e.proto.expanded for e in at.expander)
    assert not at.success
    saved = exported[-1]
    assert saved.ComparisonStepFacet.eq("Task").all() and saved.ComparisonSlopeFacet.eq("Task").all()
    assert saved.SourceAnalysisID.eq(app.build_result_analysis_identity(fit).analysis_id).all()
    assert not saved.PublicConclusionAllowed.fillna(False).any()
    first = at.table[0].value.copy()
    at.selectbox[0].select("Rater").run()
    assert len(calls) == 2 and not at.table
    at.button[0].click().run()
    assert len(calls) == 4
    at.selectbox[0].select("Task").run()
    assert len(calls) == 4
    pd.testing.assert_frame_equal(at.table[0].value, first)
    at.selectbox[1].select("Rater").run()
    assert len(calls) == 4 and not at.table
    at.button[0].click().run()
    assert len(calls) == 6 and exported[-1].ComparisonSlopeFacet.eq("Rater").all()
    at.selectbox[0].select("Rater").run()
    at.selectbox[0].select("Task").run()
    assert at.selectbox[1].value == "Rater" and len(calls) == 6 and len(at.table) == 2
    other_lang = "en" if lang == "ja" else "ja"
    at.session_state["lang"] = other_lang
    at.run()
    assert len(calls) == 6 and not at.exception
    assert _locale_value(other_lang, "report_tables.model_choice_columns.StepFacet") in at.table[1].value
    # A different source fit must not display a cached comparison from this one.
    other = deepcopy(fit); other["summary"]["LogLik"] -= 1
    at.session_state["result"] = other
    at.run()
    assert not calls[6:] and not at.table and at.selectbox[0].value is None


@pytest.mark.parametrize("wrong_field", ["source_result_key", "comparison_spec", "comparison_policy"])
def test_mismatched_cached_bundle_is_not_displayed_or_automatically_refitted(monkeypatch, wrong_field):
    calls = capture_refits(monkeypatch)
    fit = result("GPCM")
    bundle = app.compute_model_choice_comparison(fit)
    bundle[wrong_field] = "foreign-or-old"
    calls.clear()
    at = AppTest.from_function(view)
    at.session_state["result"], at.session_state["lang"] = fit, "en"
    key = app._model_comparison_cache_key(fit, "Task", "Rater")
    at.session_state[f"model_choice_comparison::{key}"] = bundle
    at.run()
    assert not at.exception and not at.table and not calls
    assert len(at.button) == 1 and not at.button[0].disabled


@pytest.mark.parametrize("lang", ["ja", "en"])
@pytest.mark.parametrize("kind", ["weights", "regularization", "population", "anchors", "categories", "data", "model"])
def test_unsupported_settings_have_localized_reasons_without_refits(monkeypatch, lang, kind):
    calls = capture_refits(monkeypatch)
    fit = result()
    if kind == "weights":
        fit["prep"]["data"]["Weight"] = 2.
    elif kind == "regularization":
        fit["config"]["facet_regularization_enabled"] = True
    elif kind == "population":
        fit["config"]["population_model"] = {"enabled": True}
    elif kind == "anchors":
        fit["config"]["anchors_active"] = True
    elif kind == "categories":
        fit["config"]["n_cat"] = 1
    elif kind == "data":
        fit["prep"].pop("data")
    elif kind == "model":
        fit["config"].pop("model")
    at = AppTest.from_function(view)
    at.session_state["result"], at.session_state["lang"] = fit, lang
    at.run()
    assert not at.exception and not calls and not at.table and not at.button
    assert _locale_value(lang, f"report_tables.model_choice_unavailable_reasons.{kind}") in at.info[0].value
