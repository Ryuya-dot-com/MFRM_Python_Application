"""The APA display uses qualified prose; a second narrative cannot bypass holds."""

from unittest.mock import MagicMock

import pandas as pd
import pytest

import streamlit_app as app


@pytest.mark.parametrize("method,free_sd", [("MML", False), ("MML", True), ("JMLE", False)])
def test_apa_display_keeps_qualified_draft_and_references_without_second_narrative(monkeypatch, method, free_sd):
    ui = MagicMock()
    ui.columns.side_effect = lambda count: [MagicMock() for _ in range(count)]
    monkeypatch.setattr(app, "st", ui)
    monkeypatch.setattr(app, "t", lambda key, **kwargs: key)
    monkeypatch.setattr(app, "render_metric_reading_notes", lambda *args: None)
    monkeypatch.setattr(app, "_render_compact_dataframe", lambda *args, **kwargs: None)
    empty_builders = (
        "build_case_interpretation_guidance", "build_submission_action_plan",
        "build_apa_report_sentence_audit", "build_role_based_action_memos",
        "build_role_based_action_checklist", "build_critical_final_review_panel",
        "build_reporting_action_bridge", "rating_scale_decision_support_table",
    )
    for name in empty_builders:
        monkeypatch.setattr(app, name, lambda *args, **kwargs: pd.DataFrame())
    for name in (
        "generate_report_ready_decision_brief", "generate_role_based_action_memos_markdown",
        "generate_role_based_action_checklist_markdown", "generate_critical_final_review_markdown",
        "generate_reporting_action_bridge_markdown",
    ):
        monkeypatch.setattr(app, name, lambda *args, **kwargs: "")
    monkeypatch.setattr(app, "current_custom_simulation_sparse_export_frames", lambda: {})
    monkeypatch.setattr(app, "build_current_run_method_reference_audit", lambda *args: pd.DataFrame({
        "PrimaryReferenceKeys": ["Louis_1982"],
    }))
    held = method == "MML"
    action = "Withhold final conclusions." if held else "Review the study context."
    monkeypatch.setattr(app, "build_publication_gate_summary", lambda *args: pd.DataFrame([{
        "GateArea": "Overall manuscript gate", "GateStatus": "Not ready" if held else "Ready",
        "Evidence": "Calibration qualification pending." if held else "Evidence reviewed.",
        "ManuscriptAction": action,
    }]))
    draft = "# APA Results Paragraph Draft\n\n" + (
        "No Results sentence is currently copyable from the APA sentence audit."
        if held else "The eligible descriptive sentence selected by the evidence audit."
    )
    monkeypatch.setattr(app, "generate_report_ready_apa_results_draft", lambda *args, **kwargs: draft)
    result = {
        "config": {"method": method, "estimate_population_sd": free_sd, "model": "PCM",
                   "facet_names": ["Rater"], "n_cat": 3},
        "prep": {"n_obs": 3, "n_person": 1, "rating_min": 0, "rating_max": 2},
        "steps": pd.DataFrame({"Estimate": [-1., 1.]}),
    }
    diagnostics = {
        "obs": pd.DataFrame({"StdResidual": [0., .1, -.1], "Score": [0, 1, 2]}),
        "measures": pd.DataFrame({"Facet": ["Rater"], "Level": ["R1"],
                                  "Estimate": [0.], "SE": [.2], "Infit": [1.]}),
    }
    app._render_apa_report(result, diagnostics)

    displayed = "\n".join(
        str(call.args[0]) for method_name in ("markdown", "code", "caption", "success", "warning", "error")
        for call in getattr(ui, method_name).call_args_list if call.args
    )
    assert draft in displayed
    assert action in displayed
    assert "Louis" in displayed
    assert "Reported 95% intervals" not in displayed
    assert "Global model fit was" not in displayed
    assert "The model fits the data adequately" not in displayed
    assert "categories function as intended" not in displayed
    assert "**Method.**" not in displayed
    assert "Copy-ready plain text" not in [call.args[0] for call in ui.expander.call_args_list]
    downloads = {call.kwargs.get("file_name"): call.kwargs.get("data") for call in ui.download_button.call_args_list}
    assert downloads["apa_results_paragraph_draft.md"] == draft.encode("utf-8")


@pytest.mark.parametrize("method", ["MML", "JMLE"])
@pytest.mark.parametrize("cells,flagged", [(0, 0), (2, 0), (2, 1)])
def test_case_guidance_never_allows_no_bias_claims_for_screened_cells(monkeypatch, method, cells, flagged):
    monkeypatch.setattr(app, "build_final_report_readiness", lambda *args: pd.DataFrame())
    monkeypatch.setattr(app, "build_publication_gate_summary", lambda *args: pd.DataFrame())
    monkeypatch.setattr(app, "build_bias_inference_audit", lambda *args: pd.DataFrame([{
        "CellsScreened": cells, "FlaggedCells": flagged, "SparseCells": 0,
    }]))
    bundle = {"Rater x Criterion": {"table": pd.DataFrame()}} if cells else {}
    frame = app.build_case_interpretation_guidance({"config": {"method": method}}, {}, bundle)
    bias = frame.loc[frame["Case"].eq("Bias/local interaction screen")]
    assert len(bias) == 1
    row = bias.iloc[0]
    assert row["Status"] == "Boundary"
    assert "withheld" in row["SaferWording"]
    assert "no-bias claims were limited" not in row["SaferWording"]
    assert "validated" in row["NextAction"]
    assert "including the cells examined" in row["ManuscriptGuardrail"]
