"""Visual uncertainty controls, not interval-coverage or fit validation."""
from contextlib import nullcontext
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


def measures():
    return pd.DataFrame({
        "Facet": ["Rater"] * 3, "Level": ["A", "B", "C"],
        "Estimate": [-.3, 0., .4], "SE": [.15, .14, .16],
        "SE_Status": ["conditional_approximation"] * 3,
        "CI_Status": ["conditional_approximation"] * 3,
        "SE_Method": ["conditional information approximation: other fitted parameters held fixed"] * 3,
        "CI_Method": ["Wald normal CI from the reported SE"] * 3,
    })


@pytest.fixture
def recorded_ui(monkeypatch):
    output = {"figures": [], "captions": [], "warnings": []}
    monkeypatch.setattr(app.st, "columns", lambda *a, **k: [nullcontext(), nullcontext()])
    monkeypatch.setattr(app.st, "subheader", lambda *a, **k: None)
    monkeypatch.setattr(app.st, "caption", lambda text, **k: output["captions"].append(text))
    monkeypatch.setattr(app.st, "warning", lambda text, **k: output["warnings"].append(text))
    monkeypatch.setattr(app.st, "info", lambda *a, **k: None)
    monkeypatch.setattr(app.st, "selectbox", lambda label, options, index=0, **k: options[index])
    monkeypatch.setattr(app.st, "plotly_chart", lambda fig, **k: output["figures"].append(fig))
    monkeypatch.setattr(app, "render_help_popover", lambda *a, **k: None)
    monkeypatch.setattr(app, "get_viz_ci_level", lambda: .95)
    return output


def interval_traces(output):
    return [trace for trace in output["figures"][0].data if trace.mode == "lines"]


def test_already_held_rows_cannot_draw_intervals_without_result_context(recorded_ui):
    frame = measures()
    frame[["SE_Status", "CI_Status"]] = "calibration_inference_withheld"
    app._draw_measures_forest_plotly({"measures": frame})
    assert not interval_traces(recorded_ui)
    assert recorded_ui["warnings"]


@pytest.mark.parametrize("result", [
    {"config": {"method": "MML", "estimate_population_sd": False}},
    {"config": {"method": "MML", "estimate_population_sd": True}},
    {"config": {"method": "JMLE", "jmle_facet_boundary_review": True}},
    {"config": {"method": "JMLE"}, "identifiability": {"status": "eta_rank_deficient"}},
    {"config": {"method": "JMLE"}, "identifiability": {"status": "audit_size_limit"}},
])
@pytest.mark.parametrize("saved_status", ["conditional_approximation", "calibration_inference_withheld"])
def test_current_hold_overrides_forest_saved_status(result, saved_status, recorded_ui):
    frame = measures()
    frame[["SE_Status", "CI_Status"]] = saved_status
    before = frame.copy(deep=True)
    app._draw_measures_forest_plotly({"measures": frame}, result=result)
    assert not interval_traces(recorded_ui)
    marker = [trace for trace in recorded_ui["figures"][0].data if trace.mode == "markers"][0]
    np.testing.assert_allclose(marker.x, frame.Estimate)
    assert recorded_ui["warnings"]
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("status", ["ok", "regularized", "rank_deficient", "not_available", "",
                                   "unsupported_covariance_scope"])
def test_unknown_or_unqualified_basis_cannot_fabricate_forest_intervals(status, recorded_ui):
    frame = measures()
    frame[["SE_Status", "CI_Status"]] = status
    app._draw_measures_forest_plotly({"measures": frame}, result={"config": {"method": "JMLE"}, "opt": SimpleNamespace(success=True)})
    assert not interval_traces(recorded_ui)


def test_valid_conditional_intervals_keep_their_actual_basis(recorded_ui):
    frame = measures()
    app._draw_measures_forest_plotly({"measures": frame}, result={"config": {"method": "JMLE"}, "opt": SimpleNamespace(success=True)})
    lines = interval_traces(recorded_ui)
    assert len(lines) == 6
    np.testing.assert_allclose(lines[0].x, [-.3 - app._ci_z(.95)*.15, -.3 + app._ci_z(.95)*.15])
    text = " ".join(recorded_ui["captions"])
    assert "Conditional intervals" in text
    assert "does not establish a difference" in text
    assert "differ at the matching level" not in text
    assert "conditional" in recorded_ui["figures"][0].layout.title.text.lower()


def test_missing_or_invalid_se_keeps_points_without_inventing_intervals(recorded_ui):
    frame = measures()
    frame["SE"] = [np.nan, 0., -.2]
    app._draw_measures_forest_plotly({"measures": frame}, result={"config": {"method": "JMLE"}, "opt": SimpleNamespace(success=True)})
    assert not interval_traces(recorded_ui)
    marker = [trace for trace in recorded_ui["figures"][0].data if trace.mode == "markers"][0]
    np.testing.assert_allclose(marker.x, frame.Estimate)


def test_missing_uncertainty_columns_keeps_finite_points(recorded_ui):
    frame = measures().drop(columns=["SE", "SE_Status", "CI_Status", "SE_Method", "CI_Method"])
    app._draw_measures_forest_plotly({"measures": frame}, result={"config": {"method": "JMLE"}, "opt": SimpleNamespace(success=True)})
    assert not interval_traces(recorded_ui)
    assert len(recorded_ui["figures"][0].data[0].x) == 3


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("method,free_sd,boundary", [
    ("MML", False, False), ("MML", True, False), ("JMLE", False, True),
])
def test_full_visuals_route_passes_current_hold_to_forest(lang, method, free_sd, boundary):
    def view(lang, method, free_sd, boundary):
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        from tests.test_visual_uncertainty_contract import measures
        st.session_state["lang"] = lang
        st.session_state["visuals_panel"] = "forest_plot"
        result = {"config": {"method": method, "estimate_population_sd": free_sd,
                             "jmle_facet_boundary_review": boundary}}
        with patch.object(app, "render_help_popover"), patch.object(app.st, "plotly_chart") as plot:
            app.show_visuals_section(result, {"measures": measures()}, force_full=True)
        st.session_state["interval_line_count"] = sum(trace.mode == "lines" for trace in plot.call_args.args[0].data)
    at = AppTest.from_function(view, args=(lang, method, free_sd, boundary)).run(timeout=30)
    assert not at.exception
    assert at.session_state["interval_line_count"] == 0
    assert at.warning
    assert all(not item.value.startswith("visuals_top.") for item in [*at.caption, *at.warning])


def test_correct_binary_model_residuals_need_not_be_normal(recorded_ui):
    responses = np.tile([0., 1.], 20)
    probability = .5
    residuals = (responses - probability) / np.sqrt(probability * (1 - probability))
    np.testing.assert_array_equal(np.unique(residuals), [-1., 1.])
    app._draw_residual_qq_plotly({"obs": pd.DataFrame({"StdResidual": residuals})})
    plot = recorded_ui["figures"][0]
    np.testing.assert_array_equal(np.unique(plot.data[0].y), [-1., 1.])
    text = " ".join(recorded_ui["captions"])
    assert "need not follow N(0, 1)" in text
    assert "descriptive" in text and "not a fit test" in text
    assert "Under a well-fitting Rasch model" not in text
    assert "target ≈ 0, 1" not in text
    assert "descriptive" in plot.layout.title.text.lower()
