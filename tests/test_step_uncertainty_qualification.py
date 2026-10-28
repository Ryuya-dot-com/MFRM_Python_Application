"""Current uncertainty holds survive step exports and restored result displays."""

from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


@pytest.fixture(scope="module")
def fitted_steps():
    result = app.mfrm_estimate(
        app._make_self_test_gradient_data(), "Person", ["Rater", "Task"], "Score",
        model="RSM", method="MML", mml_engine="EM", quad_points=5, maxit=2)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    return result, diagnostics


def test_fresh_fixed_sd_step_values_remain_technical(fitted_steps):
    result, diagnostics = fitted_steps
    steps = result["steps"]
    assert steps.SE.gt(0).all()
    assert steps.SE_Status.eq("calibration_inference_withheld").all()
    assert steps.CI_Status.eq("calibration_inference_withheld").all()
    assert steps.CI_Method.str.contains("Normal approximation", regex=False).all()
    assert steps.SE_CovarianceStatus.eq("ok").all()
    assert steps.UncertaintyCaution.str.contains("inference is withheld").all()
    assert diagnostics["measures"].SE_Status.eq("calibration_inference_withheld").all()


@pytest.mark.parametrize("free_sd", [False, True])
def test_saved_step_exports_requalify_without_changing_numbers(fitted_steps, free_sd):
    result, diagnostics = deepcopy(fitted_steps)
    result["config"]["estimate_population_sd"] = free_sd
    saved = result["steps"].assign(SE_Status="ok", CI_Status="ok")
    result["steps"] = saved
    before = saved.copy(deep=True)
    frames = app.build_result_bundle_frames(result, diagnostics)
    downloads, _ = app.collect_download_frames(
        result, diagnostics, {}, pd.DataFrame(), pd.DataFrame(), public_export_mode=False)
    apa = app._collect_apa_exportable_tables(result, diagnostics)
    for table in (frames["step_thresholds"], downloads["steps"], apa["Step thresholds"]):
        assert table.SE_Status.eq("calibration_inference_withheld").all()
        assert table.CI_Status.eq("calibration_inference_withheld").all()
        np.testing.assert_array_equal(table[["Estimate", "SE", "CI_Lower", "CI_Upper"]],
                                      before[["Estimate", "SE", "CI_Lower", "CI_Upper"]])
    pd.testing.assert_frame_equal(result["steps"], before)


def test_saved_facets_tables_get_current_hold_on_download(fitted_steps):
    result, diagnostics = fitted_steps
    source = pd.DataFrame({"Level": ["R1"], "Measure": [.3], "ModelSE": [.2],
                           "SE_Status": ["ok"], "CI_Status": ["ok"]})
    before = source.copy(deep=True)
    downloads, _ = app.collect_download_frames(
        result, diagnostics, {"Rater": source}, pd.DataFrame(), pd.DataFrame(), public_export_mode=False)
    table = downloads["facets_Rater"]
    assert table.SE_Status.eq("calibration_inference_withheld").all()
    assert table.CI_Status.eq("calibration_inference_withheld").all()
    np.testing.assert_array_equal(table.ModelSE, before.ModelSE)
    pd.testing.assert_frame_equal(source, before)


@pytest.mark.parametrize("config", [
    {"method": "MML"},
    {"method": "MML", "estimate_population_sd": True},
    {"method": "JMLE", "jmle_facet_boundary_review": True},
])
def test_saved_fair_average_statuses_cannot_override_current_hold(config):
    source = pd.DataFrame({"Level": ["R1"], "ModelSE": [.2]})
    for prefix in ("FairM", "FairZ"):
        source[f"{prefix}SE"] = .1
        source[f"{prefix}_CI_Lower"] = .2
        source[f"{prefix}_CI_Upper"] = .6
        source[f"{prefix}_SE_Status"] = "ok"
        source[f"{prefix}_CI_Status"] = "ok"
        source[f"{prefix}_SE_QualificationStatus"] = "READY"
        source[f"{prefix}_SE_PublicConclusionAllowed"] = True
    before = source.copy(deep=True)
    table = app.qualify_measure_uncertainty(source, {"config": config})
    for prefix in ("FairM", "FairZ"):
        assert table[f"{prefix}_SE_Status"].eq("calibration_inference_withheld").all()
        assert table[f"{prefix}_CI_Status"].eq("calibration_inference_withheld").all()
        assert table[f"{prefix}_SE_QualificationStatus"].eq("WITHHELD").all()
        assert not table[f"{prefix}_SE_PublicConclusionAllowed"].any()
        for column in ("ModelSE", f"{prefix}SE", f"{prefix}_CI_Lower", f"{prefix}_CI_Upper"):
            np.testing.assert_array_equal(table[column], before[column])
    pd.testing.assert_frame_equal(app.qualify_measure_uncertainty(table, {"config": config}), table)
    pd.testing.assert_frame_equal(source, before)


@pytest.mark.parametrize("with_bounds", [False, True])
def test_guided_saved_result_display_requalifies_steps_and_measures(with_bounds):
    from streamlit.testing.v1 import AppTest

    def view(with_bounds):
        from contextlib import ExitStack
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state["lang"] = "en"
        steps = pd.DataFrame({"Step": ["Step_1"], "Estimate": [.1], "SE": [.2],
                              "SE_Status": ["ok"], "CI_Status": ["ok"]})
        measures = pd.DataFrame({"Facet": ["Person"], "Level": ["P1"],
                                 "Estimate": [.4], "SE": [.2], "SE_Status": ["ok"],
                                 "CI_Status": ["ok"], "Infit": [1.], "Outfit": [1.]})
        if with_bounds:
            for frame in (steps, measures):
                frame["CI_Lower"] = frame.Estimate-.392
                frame["CI_Upper"] = frame.Estimate+.392
        result = {"config": {"method": "MML"}, "steps": steps}
        facets = measures.rename(columns={"SE": "ModelSE"}).copy()
        with ExitStack() as stack:
            for name in ("show_convergence_section", "render_guided_measures_reading_guide_tables",
                         "show_person_fit_section", "show_posterior_scoring_section",
                         "render_eb_shrinkage_section", "render_metric_reading_notes"):
                stack.enter_context(patch.object(app, name, return_value=None))
            app._render_guided_results_section(result, {"measures": measures}, {"Rater": facets})
        st.session_state["saved_statuses"] = [frame.SE_Status.tolist() for frame in (steps, measures, facets)]

    at = AppTest.from_function(view, args=(with_bounds,)).run(timeout=30)
    assert not at.exception
    tables = [frame.value for frame in at.dataframe if "SE_Status" in frame.value]
    assert len(tables) == 3
    for table in tables:
        assert table.SE_Status.eq("calibration_inference_withheld").all()
        assert table.CI_Status.eq("calibration_inference_withheld").all()
        assert ("CI_Lower" in table) == with_bounds
    assert at.session_state["saved_statuses"] == [["ok"], ["ok"], ["ok"]]


def test_full_results_requalify_saved_uncertainty_without_refitting(monkeypatch):
    from tests.test_analysis_journeys import _fit, _prepare, _run
    from tests.test_help_contextual_app_integration import _block_analysis_calls

    at = _prepare("fixed_sd")
    output = _fit(at)
    for frame in (output["result"]["steps"], output["diagnostics"]["measures"]):
        frame["SE_Status"] = "ok"
        frame["CI_Status"] = "ok"
        frame.drop(columns=["CI_Lower", "CI_Upper"], inplace=True)
    for frame in output["report_tables"].values():
        frame["SE_Status"] = "ok"
        frame["CI_Status"] = "ok"
    calls = _block_analysis_calls(monkeypatch)
    at.radio(key="app_view_density").set_value("Full")
    _run(at)
    at.selectbox(key="main_results_panel").set_value("measures")
    _run(at)
    displayed = [frame.value for frame in at.dataframe if "SE_Status" in frame.value]
    # Person and facet tables now carry the same current hold as steps/combined measures.
    assert len(displayed) == 4
    for table in displayed:
        assert table.SE_Status.eq("calibration_inference_withheld").all()
        assert table.CI_Status.eq("calibration_inference_withheld").all()
    for table in displayed[-2:]:  # saved steps and combined measures lacked interval bounds
        assert "CI_Lower" not in table
    for frame in (output["result"]["steps"], output["diagnostics"]["measures"]):
        assert frame.SE_Status.eq("ok").all()
    at.selectbox(key="main_results_panel").set_value("facets_style_tables")
    _run(at)
    displayed_facets = [frame.value for frame in at.dataframe if "ModelSE" in frame.value]
    assert len(displayed_facets) == len(output["report_tables"])
    for table, source in zip(displayed_facets, output["report_tables"].values()):
        assert table.SE_Status.eq("calibration_inference_withheld").all()
        assert table.CI_Status.eq("calibration_inference_withheld").all()
        np.testing.assert_array_equal(table.ModelSE, source.ModelSE)
        assert source.SE_Status.eq("ok").all()
    assert calls == {"estimate": 0, "refit": 0}
