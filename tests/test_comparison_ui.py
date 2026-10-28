"""Comparison describes finite matched values and preserves inference holds."""
from copy import deepcopy
import json

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value
from tests.test_measure_report_consistency import fitted_case


def snapshot(values, levels=None, facets=None):
    return {"diagnostics": {"measures": pd.DataFrame({
        "Facet": facets or ["Rater"] * len(values),
        "Level": levels or [f"R{i}" for i in range(len(values))], "Estimate": values})}}


def recorded(output):
    return {**output, "person_col": "Person", "score_col": "Score", "facet_cols": ["Rater", "Task"],
            "_input_data_fingerprint": "same-ratings", "_anchor_source_fingerprint": None,
            "_group_anchor_source_fingerprint": None, "_population_source_fingerprint": None}


def test_perfect_correlation_does_not_erase_location_difference_or_pool_facets():
    a = snapshot([0., 1., 2., 10., 12.], facets=["Rater"]*3+["Task"]*2)
    b = snapshot([5., 6., 7., 12., 10.], facets=["Rater"]*3+["Task"]*2)
    before = deepcopy(a)
    result = app._comparison_measures_diff(a, b)
    stats = result["statistics"].set_index("Facet")
    assert stats.loc["Rater", "Correlation"] == pytest.approx(1.)
    assert stats.loc["Rater", "RMSE"] == pytest.approx(5.)
    assert stats.loc["Task", "Correlation"] == pytest.approx(-1.)
    assert stats.loc["Task", "RMSE"] == pytest.approx(2.)
    assert stats.N_Overlap.tolist() == [3, 2]
    pd.testing.assert_frame_equal(a["diagnostics"]["measures"], before["diagnostics"]["measures"])


def test_nonfinite_and_unmatched_values_do_not_enter_metrics_or_chart():
    a = snapshot([0., 1., np.inf, 4., np.nan, 8.])
    b = snapshot([5., 6., 7., -np.inf, 9., 10.], levels=["R0", "R1", "R2", "R3", "R4", "other"])
    result = app._comparison_measures_diff(a, b)
    row = result["statistics"].iloc[0]
    assert (row.N_A, row.N_B, row.N_Overlap) == (6, 6, 2)
    assert row.RMSE == pytest.approx(5.)
    assert result["merged"].Level.tolist() == ["R0", "R1"]
    assert np.isfinite(result["merged"][["Estimate_A", "Estimate_B", "diff"]]).all().all()


@pytest.mark.parametrize("a,b", [([1., 1.], [2., 3.]), ([1.], [2.])])
def test_undefined_correlation_is_blank_without_divide_warnings(a, b):
    with np.errstate(all="raise"):
        row = app._comparison_measures_diff(snapshot(a), snapshot(b))["statistics"].iloc[0]
    assert np.isnan(row.Correlation) and np.isfinite(row.RMSE)


def test_no_overlap_is_explicit_and_duplicate_or_missing_identifiers_are_rejected():
    a, b = snapshot([0., 1.]), snapshot([0., 1.], levels=["x", "y"])
    result = app._comparison_measures_diff(a, b)
    assert result["n_overlap"] == 0 and result["statistics"].N_Overlap.eq(0).all()
    for levels in (["R0", "R0"], ["R0", None]):
        assert app._comparison_measures_diff(a, snapshot([0., 1.], levels=levels)) is None


def test_input_alignment_does_not_treat_missing_provenance_or_added_anchor_as_match():
    a, b = recorded({}), recorded({})
    assert app._comparison_data_alignment(a, b) == (True, [])
    for key, value, reason in [
        ("_input_data_fingerprint", None, "unknown_data"),
        ("_input_data_fingerprint", "changed", "different_data"),
        ("_anchor_source_fingerprint", "new-anchor", "anchor_inputs"),
        ("_group_anchor_source_fingerprint", "new-group", "group_anchor_inputs"),
        ("_population_source_fingerprint", "new-covariates", "population_inputs"),
        ("person_col", None, "person_column"), ("score_col", "OtherScore", "score_column"),
        ("facet_cols", None, "facet_set"),
    ]:
        aligned, reasons = app._comparison_data_alignment(a, {**b, key: value})
        assert not aligned and reason in reasons
    b.pop("_anchor_source_fingerprint")
    assert app._comparison_data_alignment(a, b) == (False, ["unknown_sources"])


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_real_comparison_holds_and_translations_do_not_modify_or_refit_saved_results(fitted_case, lang, monkeypatch):
    case, result, diagnostics = fitted_case
    held = case in {"facet_boundary", "fixed_sd", "free_sd"}
    a = recorded({"result": deepcopy(result), "diagnostics": deepcopy(diagnostics)})
    b = deepcopy(a)
    # Old saved positive labels must not override the current hold policy.
    if held:
        for snap in (a, b):
            snap["diagnostics"]["measures"]["SE_Status"] = "ok"
            snap["diagnostics"]["measures"]["CI_Status"] = "ok"
    original = deepcopy(a)
    calls = _block_analysis_calls(monkeypatch)
    figures = []
    plot = app.st.plotly_chart

    def capture(fig, **kwargs):
        figures.append(fig.to_plotly_json())
        return plot(fig, **kwargs)

    monkeypatch.setattr(app.st, "plotly_chart", capture)
    at = AppTest.from_string('''
import streamlit_app as app
app.render_comparison_panel(app.st.session_state["a"], app.st.session_state["b"])
''', default_timeout=120)
    at.session_state["lang"], at.session_state["a"], at.session_state["b"] = lang, a, b
    at.run()
    assert not at.exception
    assert _locale_value(lang, "history.criteria_hold") in [w.value for w in at.info]
    hold_text = _locale_value(lang, "comparison.calibration_hold").format(run="A")
    assert (hold_text in [w.value for w in at.warning]) == held
    assert _locale_value(lang, "comparison.measures_help") in [w.value for w in at.caption]
    assert not at.success
    visible = " ".join(w.value for widgets in (at.caption, at.markdown, at.info, at.warning, at.error) for w in widgets)
    for forbidden in ("G2", "TECHNICAL_ONLY", "effectively identical", "either method is defensible"):
        assert forbidden not in visible
    criteria = next(w for w in at.dataframe if "PublicConclusionAllowed" in w.value)
    assert not criteria.value.PublicConclusionAllowed.any()
    assert not {"QualificationStatus", "ReasonCode"}.intersection(criteria.value.columns)
    summary = next(w for w in at.dataframe if "N_Overlap" in w.value)
    assert summary.value.Facet.tolist() == diagnostics["measures"].Facet.drop_duplicates().tolist()
    assert json.loads(summary.proto.columns)["Correlation"]["label"] == _locale_value(lang, "comparison.columns.correlation")
    table = next(w.value for w in at.dataframe if "Use_A" in w.value)
    for run in ("A", "B"):
        label = _locale_value(lang, "comparison.withheld" if held else "comparison.diagnostic")
        assert table[f"Use_{run}"].eq(label).all()
        if held:
            assert table[[f"{c}_{run}" for c in ("Separation", "Reliability", "Strata")]].eq("").all().all()
    assert figures and figures[0]["layout"]["title"]["text"] == _locale_value(lang, "comparison.plot_title")
    for name in ("measures", "reliability"):
        pd.testing.assert_frame_equal(a["diagnostics"][name], original["diagnostics"][name], check_exact=True)
    pd.testing.assert_frame_equal(a["result"]["summary"], original["result"]["summary"], check_exact=True)
    assert calls == {"estimate": 0, "refit": 0}


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_selector_uses_restorable_entries_and_rejects_same_run(lang, monkeypatch):
    def output(values):
        return recorded({**snapshot(values), "result": {"config": {"method": "JMLE"}, "summary": pd.DataFrame({
            "Model": ["RSM"], "Method": ["JMLE"], "Converged": [True], "Iterations": [5]})}})
    a, b = output([0., 1., 2.]), output([5., 6., 7.])
    b.pop("_input_data_fingerprint")
    at = AppTest.from_string('''
import streamlit_app as app
if not app.st.session_state.get("seeded"):
    app.st.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY] = True
    for output in app.st.session_state["seed_outputs"]:
        app.record_run_in_history(output=output, elapsed_sec=1.)
    app.st.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY] = False
    app.record_run_in_history(output=app.st.session_state["seed_outputs"][0], elapsed_sec=1.)
    app.st.session_state["seeded"] = True
app.render_comparison_selector()
''', default_timeout=120)
    at.session_state["lang"], at.session_state["seed_outputs"] = lang, [a, b]
    calls = _block_analysis_calls(monkeypatch)
    at.run()
    assert not at.exception
    assert at.expander[0].label == _locale_value(lang, "comparison.selector_title")
    assert len(at.selectbox(key="compare_run_a").options) == 2
    assert len(set(at.selectbox(key="compare_run_a").options)) == 2
    assert at.selectbox(key="compare_run_a").label == _locale_value(lang, "comparison.run_a")
    at.selectbox(key="compare_run_b").set_value(0)
    at.button(key="compare_trigger").click().run()
    assert not at.exception
    assert _locale_value(lang, "comparison.choose_different") in [w.value for w in at.info]
    at.selectbox(key="compare_run_b").set_value(1)
    at.button(key="compare_trigger").click().run()
    assert not at.exception
    assert _locale_value(lang, "comparison.reasons.unknown_data") in ' '.join(w.value for w in at.warning)
    summary = next(w.value for w in at.dataframe if "N_Overlap" in w.value)
    assert summary.RMSE.iloc[0] == pytest.approx(5.)
    assert summary.Correlation.iloc[0] == pytest.approx(1.)
    assert not at.success and calls == {"estimate": 0, "refit": 0}
