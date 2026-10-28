"""Failed or unrecorded termination cannot authorize saved model-dependent output."""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from mfrm_app import output_qualification as qualification
from mfrm_app.reporting_facts import current_manuscript_hold
from tests.test_scoring_calibration_hold import ready_evidence
from tests.test_visual_uncertainty_contract import recorded_ui, interval_traces, measures


REASON = "stat.estimation.optimizer_termination_unverified"


@pytest.mark.parametrize("success", [False, None, "False", "True", 1, np.nan])
def test_unverified_termination_requalifies_stale_rows_and_summary(success, recorded_ui):
    result = {"config": {"method": "JMLE"}, "opt": SimpleNamespace(),
              "summary": pd.DataFrame([{"Method": "JMLE", "N": 12, "InferenceReady": True}]),
              "convergence": pd.DataFrame([{"Converged": True, "GradientNorm": .01}])}
    if success is not None:
        result["opt"].success = success
    raw = measures().assign(PersonInferenceReady=True, ReportableEstimate=lambda f: f.Estimate,
                            UncertaintyCaution="Use for screening; intervals can support decisions.")
    before = raw.copy(deep=True)
    saved = app.calc_reliability(raw)
    diagnostics = {"measures": raw, "reliability": saved,
                   "uncertainty": {"summary": pd.DataFrame([{
                       "Area": "Measure SE/CI", "Status": "conditional_approximation", "Rows": 3,
                       "Caution": "Use for screening; intervals can support decisions."}])}}
    qualified = app.qualify_measure_uncertainty(raw, result)
    assert qualified.SE_Status.eq("calibration_inference_withheld").all()
    assert qualified.UncertaintyReasonCode.eq(REASON).all()
    assert not qualified.UncertaintyCaution.str.contains("Use for screening|support decisions").any()
    assert not qualified.PersonInferenceReady.any()
    assert qualified.ReportableEstimate.isna().all()
    np.testing.assert_array_equal(qualified.Estimate, raw.Estimate)
    np.testing.assert_array_equal(qualified.SE, raw.SE)
    reliability = app.qualified_reliability_table(result, diagnostics)
    assert reliability.ReliabilityStatus.eq("withheld").all()
    assert reliability.ReliabilityReasonCode.eq(REASON).all()
    np.testing.assert_array_equal(reliability.Reliability, saved.Reliability)
    summary = app.qualified_estimation_summary(result)
    assert not summary.InferenceReady.any()
    assert summary.InferenceReadinessReason.eq(REASON).all()
    assert result["summary"].InferenceReady.all()
    convergence = app.qualified_convergence_table(result)
    assert not convergence.Converged.any()
    assert convergence.SavedConverged.all()
    assert convergence.InferenceReadinessReason.eq(REASON).all()
    assert convergence.GradientNorm.iloc[0] == .01
    assert result["convergence"].Converged.all()
    assert app.build_first_read_guide_rows(result, diagnostics)[0]["Status"] == "Do not interpret yet"
    assert not app._run_history_extract_summary({"result": result})[0]
    assert not app._comparison_extract_conv({"result": result})[0]
    audit = app.build_statistical_assumption_audit(result, diagnostics).set_index("Area")
    assert audit.loc["Measure SE/CI", "Status"] == "WITHHELD"
    assert audit.loc["Measure SE/CI", "ReasonCode"] == REASON
    assert "support decisions" not in audit.loc["Measure SE/CI", "Implication"]
    app._draw_measures_forest_plotly(diagnostics, result=result)
    assert not interval_traces(recorded_ui) and recorded_ui["warnings"]
    pd.testing.assert_frame_equal(raw, before)
    pd.testing.assert_frame_equal(app.qualify_measure_uncertainty(qualified, result), qualified)


@pytest.mark.parametrize("success", [True, np.bool_(True)])
def test_explicit_success_preserves_conditional_positive_control(success):
    result = {"config": {"method": "JMLE"}, "opt": SimpleNamespace(success=success)}
    assert qualification.optimizer_termination_established(result)
    assert not qualification.calibration_inference_hold_reason(result)
    assert not current_manuscript_hold(result)
    assert app._run_history_extract_summary({"result": result})[0]
    assert app._comparison_extract_conv({"result": result})[0]
    raw = measures()
    pd.testing.assert_frame_equal(app.qualify_measure_uncertainty(raw, result), raw)


@pytest.mark.parametrize("success", [False, None, "False"])
def test_current_hold_overrides_stale_ready_scoring_tables(ready_evidence, success):
    result, diagnostics = ready_evidence
    result["opt"] = SimpleNamespace(success=True)
    old_decision = app.build_scoring_consistency_decision(result, diagnostics)
    old_qc = app.build_quality_control_recommendations(result, diagnostics)
    assert old_decision.iloc[0].Status == "Ready"
    result["opt"] = SimpleNamespace() if success is None else SimpleNamespace(success=success)
    decision = app.build_scoring_consistency_decision(result, diagnostics).set_index("DecisionArea")
    assert decision.loc["Overall scoring consistency", "ReasonCode"] == REASON
    assert decision.loc["Overall scoring consistency", "Status"] == "Do not claim yet"
    assert decision.loc[["Exact/adjacent agreement", "Krippendorff alpha"], "Status"].eq("Ready").all()
    assert app.build_scoring_quality_first_read_summary(result, diagnostics).iloc[0].CanUseApaDraft == "No"
    draft = app.generate_apa_scoring_quality_draft(result, diagnostics, old_decision, old_qc)
    assert "supported publication-oriented interpretation" not in draft
    assert "**Do not claim yet**" in draft
    assert old_decision.iloc[0].Status == "Ready"


@pytest.fixture(scope="module")
def failed_actual_fit():
    result = app.mfrm_estimate(app._make_self_test_gradient_data(), "Person", ["Rater", "Task"],
                               "Score", method="JMLE", maxit=1)
    assert result["opt"].success is False
    assert result["identifiability"]["status"] == "eta_identified"
    assert not result["config"].get("jmle_facet_boundary_review", False)
    return result, app.mfrm_diagnostics(result, compute_pca=False, interaction_pairs=[])


def test_failed_actual_fit_and_restored_snapshot_cannot_export_ready_uncertainty(failed_actual_fit, monkeypatch):
    result, diagnostics = deepcopy(failed_actual_fit)
    # An older saved run can carry conditional status despite failed termination.
    diagnostics["measures"][["SE_Status", "CI_Status"]] = "conditional_approximation"
    result["summary"]["InferenceReady"] = True
    result["person_boundary"]["summary"]["PersonMeasureInferenceReady"] = True
    result["person_boundary"]["persons"]["PersonInferenceReady"] = True
    diagnostics["uncertainty"]["summary"]["Status"] = "conditional_approximation"
    diagnostics["uncertainty"]["summary"]["Caution"] = "Use for screening"
    state = {app._RUN_HISTORY_KEY: [{"run_id": "failed", "output_snapshot": {
        "result": result, "diagnostics": diagnostics}}]}
    monkeypatch.setattr(app.st, "session_state", state)
    assert app.restore_run_from_history("failed")
    restored = state["facets_mode_output"]
    assert restored["result"] is not result
    framesets = [app.build_result_bundle_frames(restored["result"], restored["diagnostics"]),
        app.collect_download_frames(restored["result"], restored["diagnostics"], {},
                                    pd.DataFrame(), pd.DataFrame(), public_export_mode=False)[0]]
    for frames in framesets:
        assert frames["measures"].SE_Status.eq("calibration_inference_withheld").all()
        assert frames["measures"].UncertaintyReasonCode.eq(REASON).all()
        assert frames["reliability"].ReliabilityStatus.eq("withheld").all()
        assert frames["reliability"].ReliabilityReasonCode.eq(REASON).all()
        assert frames["statistical_uncertainty_summary"].Status.eq("WITHHELD").all()
        assert frames["statistical_uncertainty_summary"].ReasonCode.eq(REASON).all()
        assert not frames["statistical_uncertainty_summary"].Caution.str.contains("Use for screening").any()
        assert not frames["jmle_person_boundary_summary"].PersonMeasureInferenceReady.any()
        assert not frames["jmle_person_boundary_persons"].PersonInferenceReady.any()
        np.testing.assert_array_equal(frames["measures"].SE, diagnostics["measures"].SE)
    apa = app._collect_apa_exportable_tables(restored["result"], restored["diagnostics"])
    assert apa["Reliability / separation"].ReliabilityReasonCode.eq(REASON).all()
    assert not apa["Person measures"].PersonInferenceReady.any()
    assert not framesets[0]["person_measures"].PersonInferenceReady.any()
    first_read = app.build_first_read_guide_rows(restored["result"], restored["diagnostics"])
    row = next(row for row in first_read if row["Check"] == "3. Reliability / separation")
    assert row["Status"] == "Do not interpret yet"
    assert "withheld" in row["What it means"]
    assert "Use person ordering claims" not in row["Next action"]
    guidance = app.build_case_interpretation_guidance(restored["result"], restored["diagnostics"])
    assert not guidance.Case.eq("High human-rater facet reliability").any()
    readiness = app.build_final_report_readiness(restored["result"], restored["diagnostics"]).set_index("Check")
    assert readiness.loc["Convergence", "Status"] != "Ready"
    assert readiness.loc["Reliability / separation", "Status"] != "Ready"
    assert diagnostics["measures"].SE_Status.eq("conditional_approximation").all()


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("success", [False, None, "False"])
def test_breadcrumb_and_restored_history_requalify_stale_success(lang, success):
    from streamlit.testing.v1 import AppTest

    def view(lang, success):
        from types import SimpleNamespace
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        st.session_state["lang"] = lang
        result = {"config": {"method": "JMLE", "model": "PCM"},
                  "opt": SimpleNamespace() if success is None else SimpleNamespace(success=success),
                  "summary": pd.DataFrame([{"Converged": True, "Method": "JMLE", "Iterations": 1}])}
        st.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY] = True
        for _ in range(2):
            entry = app.record_run_in_history(output={"result": result, "model_type": "PCM"}, elapsed_sec=.1)
            entry["converged"] = True
        app._render_run_breadcrumb(result, [])
        app.render_run_history_panel()
        app.render_comparison_selector()
    at = AppTest.from_function(view, args=(lang, success)).run(timeout=30)
    assert not at.exception
    table = next(table.value for table in at.dataframe if "converged" in table.value)
    assert not table.converged.any()
    text = " ".join(row.value for row in at.markdown)
    locale = app._load_locale(lang)
    assert locale["run_progress"]["interpretation_hold"] in text
    assert locale["history"]["not_converged"] in text
    assert all(locale["history"]["not_converged"] in option for option in at.selectbox(key="compare_run_a").options)
