"""Current calibration holds must reach scoring decisions and reporting checks."""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


HOLDS = [
    ({"method": "MML", "estimate_population_sd": False}, {}, "stat.mml_fixed_sd.numerical_qualification_pending"),
    ({"method": "MML", "estimate_population_sd": True}, {}, "stat.mml_free_sd.numerical_qualification_pending"),
    ({"method": "JMLE", "jmle_facet_boundary_review": True}, {}, "stat.jmle.facet_score_boundary_review"),
    ({"method": "JMLE"}, {"status": "eta_rank_deficient"}, "stat.jmle.identification_unqualified"),
]


@pytest.fixture
def ready_evidence(monkeypatch):
    monkeypatch.setattr(app, "_agreement_summary_for_facet", lambda *args: pd.DataFrame([
        {"ExactAgreement": .95, "AdjacentAgreement": 1., "Contexts": 20}]))
    monkeypatch.setattr(app, "_ordinal_alpha_summary_for_facet", lambda *args: pd.DataFrame([
        {"Status": "Ready", "KrippendorffAlpha": .95, "PairableUnits": 20, "Elements": 2}]))
    monkeypatch.setattr(app, "_facet_fit_review_summary", lambda *args: ("Ready", "no review flags", "Report fit."))
    monkeypatch.setattr(app, "rating_scale_decision_support_table", lambda *args: pd.DataFrame([
        {"DecisionStatus": "Ready", "PrimaryEvidence": "ordered thresholds", "RecommendedAction": "Retain categories."}]))
    return {"opt": SimpleNamespace(success=True), "config": {"method": "JMLE", "model": "PCM", "facet_names": ["Rater", "Criterion"]}}, {}


@pytest.mark.parametrize("config,audit,reason", HOLDS)
def test_scoring_hold_overrides_ready_diagnostics_and_propagates_to_draft(ready_evidence, config, audit, reason):
    result, diagnostics = ready_evidence
    legacy_decision = app.build_scoring_consistency_decision(result, diagnostics)
    legacy_qc = app.build_quality_control_recommendations(result, diagnostics)
    assert legacy_decision.iloc[0].Status == "Ready"
    result["config"].update(config)
    result["identifiability"] = audit
    saved_result = deepcopy(result)

    table = app.build_scoring_consistency_decision(result, diagnostics).set_index("DecisionArea")
    model_areas = ["Overall scoring consistency", "Rater fit/severity", "Rubric/category functioning", "Category collapse readiness"]
    assert table.loc[model_areas, "Status"].eq("Do not claim yet").all()
    assert table.loc[model_areas, "ReportBlocking"].eq("Yes").all()
    assert table.loc[model_areas, "ReasonCode"].eq(reason).all()
    assert table.loc[["Exact/adjacent agreement", "Krippendorff alpha"], "Status"].eq("Ready").all()
    assert not table.loc[model_areas].astype(str).apply(lambda s: s.str.contains("Proceed with|Copy only|training/moderation", regex=True)).any().any()
    summary = app.build_scoring_quality_first_read_summary(result, diagnostics).iloc[0]
    assert summary.CanUseApaDraft == "No"
    qc = app.build_quality_control_recommendations(result, diagnostics)
    assert qc.Recommendation.eq("Withheld").all()
    assert qc.Action.str.contains("calibration hold").all()
    draft = app.generate_apa_scoring_quality_draft(result, diagnostics, legacy_decision, legacy_qc)
    assert "**Do not claim yet**" in draft
    assert "supported publication-oriented interpretation" not in draft
    assert result == saved_result
    assert legacy_decision.iloc[0].Status == "Ready"


def test_unheld_scoring_keeps_descriptive_positive_control(ready_evidence):
    result, diagnostics = ready_evidence
    table = app.build_scoring_consistency_decision(result, diagnostics)
    assert table.iloc[0].Status == "Ready"
    assert table.iloc[0].ReportBlocking == "No"


def _checklist_view(lang, config, audit, uncertainty):
    from types import SimpleNamespace
    import numpy as np
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    st.session_state["lang"] = lang
    result = {"opt": SimpleNamespace(success=True), "config": {"model": "PCM", **config}, "identifiability": audit}
    measures = pd.DataFrame({"Facet": ["Person", "Rater"], "Estimate": [.1, -.1], "SE": [.2, .2],
        "SE_Status": ["conditional_approximation"] * 2, "SE_Method": ["conditional information"] * 2,
        "CI_Status": ["conditional_approximation"] * 2, "CI_Method": ["Wald normal"] * 2,
        "CI_Level": [.95, .95], "CI_Lower": [-.3, -.5], "CI_Upper": [.5, .3], "Infit": [1., 1.]})
    if uncertainty == "missing_bounds":
        measures = measures.drop(columns=["CI_Lower", "CI_Upper"])
    elif uncertainty == "one_unavailable":
        measures.loc[1, "SE"] = np.nan
    elif uncertainty == "missing_status":
        measures = measures.drop(columns=["SE_Status", "CI_Status"])
    elif uncertainty == "regularized":
        measures.loc[1, "SE_Status"] = "regularized"
        measures.loc[1, "CI_Status"] = "regularized"
    diagnostics = {"measures": measures, "reliability": pd.DataFrame({"Reliability": [.9], "SE_StatusSummary": ["ok"]})}
    app._render_reporting_checklist(result, diagnostics, {"Rater:Criterion": {"saved_status": "Ready"}})


@pytest.mark.parametrize("lang", ["en", "ja"])
@pytest.mark.parametrize("config,audit,reason", HOLDS)
def test_checklist_cannot_check_held_uncertainty_or_bias(lang, config, audit, reason):
    at = AppTest.from_function(_checklist_view, args=(lang, config, audit, "complete")).run(timeout=45)
    assert not at.exception
    text = "\n".join(item.value for item in at.markdown)
    assert at.warning
    for key in ["checklist_uncertainty", "checklist_interval", "checklist_bias_conclusions"]:
        label = app._load_locale(lang)["report_top"][key]
        assert f"- [ ] {label}" in text
    assert "- [x] **Separation" not in text
    assert "- [x] **Fixed chi-square" not in text
    assert "- [x] **RMSE" not in text
    assert "Significant interactions" not in text
    assert "Estimate ± 1.96" not in text


@pytest.mark.parametrize("uncertainty,se_checked,ci_checked", [
    ("complete", True, True), ("missing_bounds", True, False),
    ("one_unavailable", False, False), ("missing_status", False, False), ("regularized", False, False),
])
def test_checklist_requires_each_rows_uncertainty_metadata(uncertainty, se_checked, ci_checked):
    at = AppTest.from_function(_checklist_view, args=("en", {"method": "JMLE"}, {}, uncertainty)).run(timeout=45)
    assert not at.exception
    text = "\n".join(item.value for item in at.markdown)
    for key, checked in [("checklist_uncertainty", se_checked), ("checklist_interval", ci_checked)]:
        label = app._load_locale("en")["report_top"][key]
        assert f"- [{'x' if checked else ' '}] {label}" in text
