import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from mfrm_app.jmle_boundary import facet_score_boundary_audit, HOLD_REASON


def separated_data():
    return pd.DataFrame([
        dict(Person=f"P{p}", Rater=r, Task=t, Score=s)
        for p in range(8) for r, s in [("R1", 1), ("R2", 0)] for t in ["T1", "T2"]
    ])


def test_facet_separation_holds_inference_even_with_interior_person_scores():
    result = app.mfrm_estimate(separated_data(), person_col="Person", facet_cols=["Rater", "Task"],
        score_col="Score", rating_min=0, rating_max=1, keep_original=True)
    summary = result["summary"].iloc[0]
    assert summary.Converged and summary.EtaStructurallyIdentified
    assert summary.ExtremeJMLEPersons == 0
    assert not summary.InferenceReady and not summary.PersonMeasureInferenceReady
    assert summary.InferenceReadinessReason == HOLD_REASON
    assert np.isfinite(result["facets"]["person"].Estimate).all()
    assert result["facets"]["person"].ReportableEstimate.isna().all()
    assert result["facet_boundary"].RequiresReview.sum() == 2
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False)
    assert app.build_first_read_guide_rows(result, diagnostics)[0]["Status"] == "Do not interpret yet"
    readiness = app.build_final_report_readiness(result, diagnostics).set_index("Check")
    assert readiness.loc["Convergence", "Status"] != "Ready"
    gate = app.build_publication_gate_summary(result, diagnostics).set_index("GateArea")
    assert gate.loc["Overall manuscript gate", "GateStatus"] == "Not ready"
    guide = app.build_manuscript_claim_guide(result, diagnostics).set_index("ManuscriptArea")
    claim = guide.loc["Convergence and final interpretability"]
    assert claim.ClaimStatus == "Do not claim"
    assert "finite unpenalized estimate" in claim.SafeManuscriptWording
    sentences = app.build_apa_report_sentence_audit(result, diagnostics).set_index("SentenceRole")
    sentence = sentences.loc["Convergence"]
    assert sentence.DraftSentence == claim.SafeManuscriptWording
    assert sentence.CopyDecision.startswith("Do not copy")
    panel = app.build_report_ready_summary_panel(result, diagnostics).set_index("ReportArea")
    assert panel.loc["Overall manuscript conclusion", "ReportStatus"] == "Do not report yet"
    assert "withheld" in panel.loc["Overall manuscript conclusion", "SafeWording"]
    measures = diagnostics["measures"]
    assert measures.SE_Status.eq("calibration_inference_withheld").all()
    assert measures.CI_Status.eq("calibration_inference_withheld").all()
    assert measures.UncertaintyCaution.str.contains("Technical terminal-point").all()
    assumptions = app.build_statistical_assumption_audit(result, diagnostics).set_index("Area")
    assert assumptions.loc["Reliability / separation uncertainty", "Status"] == "WITHHELD"
    assert "_skip_reason" in app.estimate_bias_interaction(result, diagnostics, "Rater", "Task")
    assert "JMLE facet score-boundary screen" in app._collect_apa_exportable_tables(result, diagnostics)
    frames = app.build_result_bundle_frames(result, diagnostics)
    assert "jmle_facet_boundary" in frames
    pd.testing.assert_frame_equal(frames["report_ready_summary_panel"].set_index("ReportArea"), panel)
    assert frames["measures"].SE_Status.eq("calibration_inference_withheld").all()


@pytest.mark.parametrize("case, expected", [("free", 2), ("interior", 0), ("fixed", 0), ("group_fixed", 0), ("mml", 0)])
def test_screen_respects_scores_and_effect_constraints(case, expected):
    data = separated_data()
    if case == "interior":
        data.loc[data.Person.eq("P0"), "Score"] = 1 - data.loc[data.Person.eq("P0"), "Score"]
    spec = app.build_facet_constraint(["R1", "R2"], centered=True)
    if case == "fixed":
        spec["anchors"] = np.zeros(2)
    elif case == "group_fixed":
        spec["groups"] = ["g1", "g2"]
    prep = dict(data=data, rating_min=0, rating_max=1)
    config = dict(method="MML" if case == "mml" else "JMLE", facet_names=["Rater"], facet_specs={"Rater": spec})
    assert facet_score_boundary_audit(prep, config).RequiresReview.sum() == expected
