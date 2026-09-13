from __future__ import annotations

import inspect
import pandas as pd

import streamlit_app as app


def test_result_router_renders_snapshot_before_next_check():
    source = inspect.getsource(app._render_guided_goal_router)

    assert source.index("_render_guided_run_snapshot(action_plan)") < source.index(
        "guided.open_next_button"
    )


def test_bias_panel_surfaces_raw_threshold_evidence_and_download():
    source = inspect.getsource(app.show_bias_section)

    assert "build_bias_decision_stability_audit(" in source
    assert "Bias threshold / rounding evidence" in source
    assert "Download bias threshold audit CSV" in source


def test_fit_styling_uses_raw_decision_frame_not_rounded_display(monkeypatch):
    displayed = pd.DataFrame({"Facet": ["Rater"], "Infit": [1.500]})
    raw = pd.DataFrame({"Facet": ["Rater"], "Infit": [1.5004]})
    seen: list[float] = []
    original = app._decision_stability.classify_fit_mnsq

    def spy(value):
        seen.append(float(value))
        return original(value)

    monkeypatch.setattr(app._decision_stability, "classify_fit_mnsq", spy)
    styler = app.style_fit_columns(displayed, decision_df=raw)
    styler._compute()

    assert 1.5004 in seen
    assert original(raw.loc[0, "Infit"]) == "noisy"
    assert displayed.loc[0, "Infit"] == 1.5


def test_anchor_audit_reports_unique_descriptive_coverage_without_share_gate():
    prep = {
        "facet_names": ["Rater"],
        "levels": {
            "Person": ["P1", "P2"],
            "Rater": ["R1", "R2", "R3", "R4"],
        },
        "data": pd.DataFrame({
            "Person": ["P1", "P1", "P2", "P2"],
            "Rater": ["R1", "R2", "R3", "R4"],
            "Score": [1, 2, 3, 4],
        }),
    }
    fixed = pd.DataFrame({"Facet": ["Rater"], "Level": ["R1"], "Anchor": [0.0]})
    grouped = pd.DataFrame({
        "Facet": ["Rater", "Rater"],
        "Level": ["R1", "R2"],
        "Group": ["Common", "Common"],
        "GroupValue": [0.0, 0.0],
    })

    audit = app.audit_mfrm_anchors(
        prep,
        anchor_df=fixed,
        group_anchor_df=grouped,
        min_common_anchors=2,
        min_obs_per_element=1,
    )
    row = audit["summary"].loc[audit["summary"]["Facet"] == "Rater"].iloc[0]

    assert row["AnchoredLevelsTotal"] == 2
    assert row["UnanchoredLevels"] == 2
    assert row["AnchorShare"] == 0.5
    assert row["Status"] == "Linked"
    assert "no universal anchor percentage" in row["CoverageInterpretation"]


def test_anchor_audit_rejects_nonfinite_group_values_without_implicit_zero():
    prep = {
        "facet_names": ["Rater"],
        "levels": {"Person": ["P1"], "Rater": ["R1", "R2"]},
        "data": pd.DataFrame({
            "Person": ["P1", "P1"],
            "Rater": ["R1", "R2"],
            "Score": [1, 2],
        }),
    }
    grouped = pd.DataFrame({
        "Facet": ["Rater"],
        "Level": ["R1"],
        "Group": ["Common"],
        "GroupValue": [float("nan")],
    })

    audit = app.audit_mfrm_anchors(
        prep,
        group_anchor_df=grouped,
        min_obs_per_element=1,
    )

    assert audit["valid_group_anchors"].empty
    assert "invalid_group_anchor_value" in set(audit["issues"]["Type"])
    rater = audit["summary"].loc[audit["summary"]["Facet"] == "Rater"].iloc[0]
    assert rater["AnchoredLevelsTotal"] == 0


def test_scoring_fit_summary_does_not_convert_missing_fit_to_ready():
    missing = {
        "fit": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Infit": [float("nan")],
            "Outfit": [float("nan")],
        })
    }
    status, evidence, _ = app._facet_fit_review_summary(missing, "Rater")

    assert status == "Missing evidence"
    assert "no finite" in evidence

    partial = {
        "fit": pd.DataFrame({
            "Facet": ["Rater", "Rater"],
            "Level": ["R1", "R2"],
            "Infit": [1.0, float("nan")],
            "Outfit": [1.1, float("nan")],
        })
    }
    status, evidence, _ = app._facet_fit_review_summary(partial, "Rater")

    assert status == "Review"
    assert "R2" in evidence


def _boundary_bias_bundle() -> dict:
    return {
        "Rater x Task": {
            "facet_a": "Rater",
            "facet_b": "Task",
            "table": pd.DataFrame([{
                "FacetA": "Rater",
                "FacetA_Level": "R1",
                "FacetB": "Task",
                "FacetB_Level": "T1",
                "Bias Size": 0.49996,
                "S.E.": 0.20,
                "t": 2.4998,
                "d.f.": 19,
                "Prob.": 0.04996,
                "ObsN": 20,
            }]),
        }
    }


def test_bias_audit_surfaces_threshold_sensitive_cells():
    diagnostics = {
        "obs": pd.DataFrame({
            "Person": ["P1", "P1", "P2", "P2"],
            "Rater": ["R1", "R2", "R1", "R2"],
            "Task": ["T1", "T1", "T2", "T2"],
            "Observed": [1, 2, 2, 3],
        })
    }
    result = {"config": {"facet_names": ["Rater", "Task"]}}

    audit = app.build_bias_inference_audit(
        _boundary_bias_bundle(),
        result,
        diagnostics,
    )
    row = audit.iloc[0]

    assert row["Status"] == "Review"
    assert row["BoundarySensitiveCells"] == 1
    assert row["BoundarySensitiveDecisions"] >= 1
    assert "raw/display" in row["BoundarySensitivitySummary"]


def test_one_click_bundle_contains_fit_and_bias_decision_stability_evidence():
    result = {
        "config": {"facet_names": ["Rater", "Task"]},
        "summary": pd.DataFrame({"Metric": ["Model"], "Value": ["RSM"]}),
        "facets": {},
    }
    diagnostics = {
        "measures": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Estimate": [0.1],
            "Infit": [1.5004],
            "Outfit": [1.0],
        }),
        "fit": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Infit": [1.5004],
            "Outfit": [1.0],
        }),
    }

    frames = app.build_result_bundle_frames(
        result,
        diagnostics,
        all_bias_results=_boundary_bias_bundle(),
    )

    assert "fit_decision_stability_audit" in frames
    assert "fit_decision_stability_summary" in frames
    assert "bias_decision_stability_audit" in frames
    assert frames["fit_decision_stability_summary"].iloc[0]["Status"] == "Review"


def test_full_download_contains_the_same_decision_stability_evidence():
    result = {
        "config": {"facet_names": ["Rater", "Task"]},
        "summary": pd.DataFrame({"Metric": ["Model"], "Value": ["RSM"]}),
        "prep": {},
        "facets": {},
    }
    diagnostics = {
        "measures": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Infit": [1.5004],
            "Outfit": [1.0],
        }),
        "fit": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R1"],
            "Infit": [1.5004],
            "Outfit": [1.0],
        }),
    }

    frames, _ = app.collect_download_frames(
        result,
        diagnostics,
        {},
        pd.DataFrame(),
        pd.DataFrame(),
        public_export_mode=True,
        all_bias_results=_boundary_bias_bundle(),
    )

    assert "fit_decision_stability_audit" in frames
    assert "fit_decision_stability_summary" in frames
    assert "bias_decision_stability_audit" in frames
    assert frames["fit_decision_stability_summary"].iloc[0]["Status"] == "Review"
