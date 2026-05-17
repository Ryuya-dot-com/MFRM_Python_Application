from __future__ import annotations

import pandas as pd

import streamlit_app as app


def _bias_table(*, strong: bool = True, sparse: bool = True) -> pd.DataFrame:
    rows = [
        {
            "FacetA": "Rater",
            "FacetA_Level": "R1",
            "FacetB": "Task",
            "FacetB_Level": "T1",
            "Bias Size": 0.75 if strong else 0.10,
            "S.E.": 0.12,
            "t": 6.25 if strong else 0.83,
            "d.f.": 40.0,
            "Prob.": 0.001 if strong else 0.41,
            "ObsN": 24,
            "LR ChiSq": 11.2 if strong else 0.4,
            "LR Prob.": 0.0008 if strong else 0.53,
            "Profile CI Lower": 0.38 if strong else -0.14,
            "Profile CI Upper": 1.10 if strong else 0.35,
            "Profile CI Level": 0.95,
            "Profile CI Status": "ok",
            "Likelihood Basis": "conditional profile likelihood",
            "SE Basis": "conditional information",
            "Multiplicity Basis": "cellwise screening statistics",
            "InferenceTier": "screening",
        },
        {
            "FacetA": "Rater",
            "FacetA_Level": "R2",
            "FacetB": "Task",
            "FacetB_Level": "T2",
            "Bias Size": -0.62 if sparse else -0.08,
            "S.E.": 0.44,
            "t": -1.41 if sparse else -0.18,
            "d.f.": 4.0,
            "Prob.": 0.19 if sparse else 0.86,
            "ObsN": 3 if sparse else 18,
            "LR ChiSq": 1.9 if sparse else 0.1,
            "LR Prob.": 0.17 if sparse else 0.75,
            "Profile CI Lower": -1.45 if sparse else -0.40,
            "Profile CI Upper": 0.18 if sparse else 0.24,
            "Profile CI Level": 0.95,
            "Profile CI Status": "ok",
            "Likelihood Basis": "conditional profile likelihood",
            "SE Basis": "conditional information",
            "Multiplicity Basis": "cellwise screening statistics",
            "InferenceTier": "screening",
        },
    ]
    return pd.DataFrame(rows)


def _result() -> dict:
    return {
        "config": {
            "model": "GPCM",
            "method": "MML",
            "facet_names": ["Rater", "Task"],
        },
        "prep": {},
    }


def _connected_diagnostics() -> dict:
    return {
        "obs": pd.DataFrame(
            {
                "Person": ["P1", "P1", "P2", "P3"],
                "Rater": ["R1", "R2", "R2", "R1"],
                "Task": ["T1", "T1", "T2", "T2"],
                "Observed": [3, 4, 2, 5],
            }
        )
    }


def _disconnected_diagnostics() -> dict:
    return {
        "obs": pd.DataFrame(
            {
                "Person": ["P1", "P2"],
                "Rater": ["R1", "R2"],
                "Task": ["T1", "T2"],
                "Observed": [3, 4],
            }
        )
    }


def test_dff_bias_screening_table_carries_conditional_inference_metadata():
    dff = app.build_dff_bias_screening_table({"table": _bias_table()})

    expected = {
        "SparseCell",
        "ClaimStatus",
        "LRChiSq",
        "LRProb",
        "ProfileCILower",
        "ProfileCIUpper",
        "ProfileCIStatus",
        "LikelihoodBasis",
        "SEBasis",
        "InferenceTier",
        "MultiplicityBasis",
        "MultiplicityFamily",
        "ReportingCaution",
    }
    assert expected.issubset(dff.columns)
    assert "Strong review" in set(dff["EvidenceLevel"])
    assert "Sparse cell" in set(dff["EvidenceLevel"])
    assert dff["ClaimStatus"].astype(str).str.contains("caveat|Do not claim", regex=True).any()
    assert dff["LikelihoodBasis"].astype(str).str.contains("conditional").all()
    assert dff["MultiplicityFamily"].astype(str).str.contains("alpha=0.05").all()


def test_bias_inference_audit_marks_disconnected_design_as_caveated():
    audit = app.build_bias_inference_audit(
        {"Rater x Task": {"table": _bias_table(), "facet_a": "Rater", "facet_b": "Task"}},
        _result(),
        _disconnected_diagnostics(),
    )

    assert len(audit) == 1
    row = audit.iloc[0]
    assert row["Status"] == "Review"
    assert row["ClaimStatus"] == "Report with caveat"
    assert row["FlaggedCells"] >= 1
    assert row["SparseCells"] >= 1
    assert "disconnected" in row["CommonScaleStatus"]
    assert "Conditional screening" in row["InferenceScope"]
    assert "Holm and BH" in row["MultiplicityScope"]


def test_bias_inference_audit_allows_only_scoped_ready_when_no_flags():
    audit = app.build_bias_inference_audit(
        {"Rater x Task": {"table": _bias_table(strong=False, sparse=False), "facet_a": "Rater", "facet_b": "Task"}},
        _result(),
        _connected_diagnostics(),
    )

    row = audit.iloc[0]
    assert row["Status"] == "Ready"
    assert row["ClaimStatus"] == "Ready"
    assert row["FlaggedCells"] == 0
    assert row["SparseCells"] == 0
    assert "conditional screen" in row["RecommendedAction"]
    assert "connected" in row["CommonScaleStatus"]


def test_statistical_assumption_audit_uses_bias_inference_audit_counts():
    audit = app.build_statistical_assumption_audit(
        _result(),
        _connected_diagnostics(),
        {"Rater x Task": {"table": _bias_table(), "facet_a": "Rater", "facet_b": "Task"}},
    )

    row = audit.loc[audit["Area"] == "Bias / local interaction inference"].iloc[0]
    assert row["Status"] == "Review"
    assert "flagged" in row["Evidence"]
    assert "strong-review" in row["Evidence"]
    assert "bias_inference_audit" in row["RecommendedAction"]


def test_result_bundle_exports_bias_inference_audit():
    frames = app.build_result_bundle_frames(
        _result(),
        _connected_diagnostics(),
        all_bias_results={"Rater x Task": {"table": _bias_table(), "facet_a": "Rater", "facet_b": "Task"}},
    )

    assert "bias_inference_audit" in frames
    assert not frames["bias_inference_audit"].empty
