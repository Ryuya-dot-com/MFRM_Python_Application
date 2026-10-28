"""Output qualification controls, not interval-coverage validation."""
import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def eb_measures():
    return pd.DataFrame({
        "Facet": ["Rater"] * 3, "Level": ["R0", "R1", "R2"],
        "Estimate": [-.1, 0., .1], "SE": [.2] * 3,
        "SE_Status": ["calibration_inference_withheld"] * 3,
    })


def test_zero_variance_eb_boundary_is_not_a_zero_width_confidence_interval():
    measures = eb_measures()
    before = measures.copy(deep=True)
    shr = app.compute_empirical_bayes_shrinkage(measures)
    table = shr["measures"]
    np.testing.assert_array_equal(table.ShrinkageTau2, 0.)
    np.testing.assert_array_equal(table.ShrunkEstimate, 0.)
    np.testing.assert_array_equal(table.ShrunkSE, 0.)
    fig = app._make_eb_shrinkage_figure(shr)
    assert fig is not None
    assert all(not trace.error_x.visible and trace.error_x.array is None for trace in fig.data)
    assert "CI" not in fig.layout.title.text
    assert "uncertainty unqualified" in fig.layout.title.text.lower()
    assert table.ShrinkageCI_Status.eq("unqualified_plugin_uncertainty").all()
    assert table.ShrinkageReasonCode.eq("stat.eb_shrinkage.uncertainty_unqualified").all()
    assert table.ShrinkageCaution.str.contains("zero", case=False).all()
    assert table.ShrinkagePriorMean.eq(0).all()
    assert table.ShrinkageCaution.str.contains("logit origin", regex=False).all()
    pd.testing.assert_frame_equal(measures, before)


@pytest.mark.parametrize("config", [
    {"method": "MML", "estimate_population_sd": False},
    {"method": "MML", "estimate_population_sd": True},
    {"method": "JMLE", "jmle_facet_boundary_review": True},
])
def test_assumption_audit_does_not_restore_stale_ready_status(config):
    diagnostics = {
        "uncertainty": {"summary": pd.DataFrame([
            {"Area": "Measure SE/CI", "Status": "ok", "Rows": 3, "Method": "legacy"},
            {"Area": "MML observed-information covariance", "Status": "ok", "Rows": 1,
             "ClaimStatus": "Ready", "Method": "legacy curvature"},
        ])},
        "reliability": pd.DataFrame([{"Facet": "Rater", "SE_StatusSummary": "ok=3"}]),
    }
    original = diagnostics["uncertainty"]["summary"].copy(deep=True)
    audit = app.build_statistical_assumption_audit({"config": config}, diagnostics).set_index("Area")
    for area in ("Measure SE/CI", "MML observed-information covariance", "Reliability / separation uncertainty"):
        assert audit.loc[area, "Status"] == "WITHHELD"
    if config["method"] == "MML" and not config["estimate_population_sd"]:
        assert audit.loc["MML fixed population prior SD", "Status"] == "WITHHELD"
        assert app.FIXED_SD_MML_INFERENCE_HOLD_REASON in audit.loc["MML fixed population prior SD", "Evidence"]
    pd.testing.assert_frame_equal(diagnostics["uncertainty"]["summary"], original)


def test_eap_reliability_formula_survives_inference_hold():
    measures = pd.DataFrame({
        "Facet": ["Person"] * 3, "Estimate": [-2., 0., 2.], "SE": [1.] * 3,
        "Infit": [1.] * 3, "Outfit": [1.] * 3,
        "SE_Status": ["posterior_eap_sd"] * 3,
        "SE_Method": ["MML latent-ability posterior SD from quadrature"] * 3,
    })
    before = app.calc_reliability(measures).iloc[0]
    held = app.qualify_measure_uncertainty(measures, {"config": {"method": "MML"}})
    after = app.calc_reliability(held).iloc[0]
    assert before.Reliability == pytest.approx(.8)
    assert after.Reliability == pytest.approx(before.Reliability)
    assert after.Separation == pytest.approx(before.Separation)
    assert after.Strata == pytest.approx(before.Strata)
    assert "calibration_inference_withheld" in after.SE_StatusSummary
    assert "Var(EAP)" in after.ReliabilityBasis


@pytest.mark.parametrize("free_sd", [False, True])
def test_saved_eb_tables_requalify_at_both_export_boundaries(free_sd):
    result = app.mfrm_estimate(
        app._make_self_test_gradient_data(), "Person", ["Rater", "Task"], "Score",
        model="RSM", method="MML", mml_engine="EM", quad_points=5, maxit=2)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])
    result["config"]["estimate_population_sd"] = free_sd
    shr = app.compute_empirical_bayes_shrinkage(eb_measures())
    for key in ("measures", "report"):
        shr[key] = shr[key].drop(columns=["ShrinkageCI_Status", "ShrinkageReasonCode", "ShrinkageCaution"])
    shr["measures"]["SE_Status"] = "ok"
    original = shr["measures"].copy(deep=True)
    diagnostics["eb_shrinkage"] = shr
    saved_reliability = diagnostics["reliability"].drop(columns=["ReliabilityStatus", "ReliabilityReasonCode"])
    diagnostics["reliability"] = saved_reliability
    bundles = [app.build_result_bundle_frames(result, diagnostics),
               app.collect_download_frames(result, diagnostics, {}, pd.DataFrame(), pd.DataFrame(),
                                           public_export_mode=False)[0]]
    for frames in bundles:
        assert frames["reliability"].ReliabilityStatus.eq("withheld").all()
        pd.testing.assert_frame_equal(frames["reliability"][saved_reliability.columns], saved_reliability)
        for name in ("eb_shrinkage_measures", "eb_shrinkage_report"):
            assert frames[name].ShrinkageCI_Status.eq("unqualified_plugin_uncertainty").all()
        exported = frames["eb_shrinkage_measures"]
        assert exported.SE_Status.eq("calibration_inference_withheld").all()
        pd.testing.assert_frame_equal(exported[original.columns.drop("SE_Status")],
                                      original.drop(columns="SE_Status"))
    pd.testing.assert_frame_equal(shr["measures"], original)


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_eb_live_panel_explains_held_point_only_output(lang):
    from streamlit.testing.v1 import AppTest

    def view(lang):
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state["lang"] = lang
        measures = pd.DataFrame({"Facet": ["Rater"]*3, "Level": ["R0", "R1", "R2"],
                                 "Estimate": [-.1, 0., .1], "SE": [.2]*3, "SE_Status": ["ok"]*3})
        shr = app.compute_empirical_bayes_shrinkage(measures)
        with patch.object(app, "_offer_fig_download"):
            app.render_eb_shrinkage_section({"config": {"method": "MML"}},
                                          {"eb_shrinkage": shr, "eb_shrinkage_enabled": True})
    at = AppTest.from_function(view, args=(lang,)).run(timeout=30)
    assert not at.exception
    assert len(at.warning) == 1
    assert "ShrunkSE" in at.warning[0].value
    assert "report_tables." not in at.warning[0].value
    measures = [frame.value for frame in at.dataframe if "ShrunkEstimate" in frame.value][0]
    assert measures.SE_Status.eq("calibration_inference_withheld").all()
