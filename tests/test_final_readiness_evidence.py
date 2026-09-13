from __future__ import annotations

import io
import json
import zipfile

import pandas as pd
import pytest

import streamlit_app as app
from mfrm_app import evidence
from mfrm_app import readiness


class _Opt:
    def __init__(self, success: bool):
        self.success = success
        self.message = "converged" if success else "iteration limit"


def _result(*, converged: bool) -> dict:
    return {
        "opt": _Opt(converged),
        "config": {
            "model": "RSM",
            "method": "JMLE",
            "facet_names": ["Rater", "Task"],
            "input_data_fingerprint": "fixture_data_sha256",
            "analysis_config_fingerprint": "fixture_config_sha256",
        },
        "prep": {},
        "facets": {"person": pd.DataFrame(), "others": pd.DataFrame()},
    }


def _diagnostics(*, review_residuals: bool) -> dict:
    flagged = 6 if review_residuals else 5
    return {
        "obs": pd.DataFrame(
            {"StdResidual": [2.1] * flagged + [0.0] * (100 - flagged)}
        ),
        "reliability": pd.DataFrame(),
        "pca_enabled": False,
    }


def test_final_readiness_preserves_legacy_columns_and_adds_valid_evidence_ledger():
    result = _result(converged=True)
    first = app.build_final_report_readiness(
        result,
        _diagnostics(review_residuals=False),
        all_bias_results={},
    )
    second = app.build_final_report_readiness(
        result,
        _diagnostics(review_residuals=False),
        all_bias_results={},
    )

    assert tuple(first.columns[:5]) == readiness.LEGACY_READINESS_COLUMNS
    assert {
        "SchemaVersion",
        "EvidenceID",
        "EvidenceKey",
        "AnalysisID",
        "ComputationState",
        "StabilityState",
        "ContractRequired",
        "ReasonCode",
    }.issubset(first.columns)
    assert first["EvidenceID"].is_unique
    assert first["AnalysisID"].nunique() == 1
    assert first["StabilityState"].eq("NOT_ASSESSED").all()
    pd.testing.assert_series_equal(first["EvidenceID"], second["EvidenceID"])
    pd.testing.assert_series_equal(first["AnalysisID"], second["AnalysisID"])
    assert (first["Required"].eq("Yes") == first["ContractRequired"]).all()

    identity = evidence.AnalysisIdentity.from_payload(first.attrs["analysis_identity"])
    records = tuple(
        evidence.EvidenceRecord.from_payload(payload)
        for payload in first.attrs["evidence_records"]
    )
    evidence.validate_evidence_records(records, identity=identity)


def test_final_readiness_separates_hold_caution_and_not_assessable_reviews():
    frame = app.build_final_report_readiness(
        _result(converged=False),
        _diagnostics(review_residuals=True),
        all_bias_results={},
    ).set_index("Check")

    assert frame.loc["Convergence", "Status"] == "Review"
    assert frame.loc["Convergence", "ComputationState"] == "HOLD"
    assert frame.loc["Convergence", "ReasonCode"] == "evidence.run_failed"

    assert frame.loc["Global residual fit", "Status"] == "Review"
    assert frame.loc["Global residual fit", "ComputationState"] == "CAUTION"
    assert frame.loc["Global residual fit", "ReasonCode"] == "evidence.limited"

    assert frame.loc["Dimensionality screen", "Status"] == "Review"
    assert frame.loc["Dimensionality screen", "ComputationState"] == "NOT_ASSESSABLE"
    assert frame.loc["Bias / local interaction", "ComputationState"] == "NOT_ASSESSABLE"
    assert frame.loc["Anchor / linking audit", "ComputationState"] == "NOT_ASSESSABLE"
    assert frame.loc["Anchor / linking audit", "SourceArtifactsJSON"] == "[]"


@pytest.mark.parametrize("converged", [True, False])
def test_free_sd_em_termination_does_not_qualify_final_interpretation(converged):
    result = _result(converged=converged)
    result["config"].update(method="MML", estimate_population_sd=True,
                            estimated_population_sd=1.3, population_prior_sd=1.3)
    diagnostics = _diagnostics(review_residuals=False)
    row = app.build_final_report_readiness(result, diagnostics).set_index("Check").loc["Convergence"]
    assert row["Status"] == "Review"
    assert row["ComputationState"] == "HOLD"
    assert row["Required"] == "Yes"
    assert row["ReasonCode"] == (
        app.FREE_SD_MML_INFERENCE_HOLD_REASON if converged else "evidence.run_failed"
    )
    first_read = app.build_first_read_guide_rows(result, diagnostics)
    assert first_read[0]["Status"] == "Do not interpret yet"
    assert app.first_read_guide_should_expand(first_read)
    action = app.build_guided_action_plan(result, diagnostics).iloc[0]
    assert action["Check"] == "1. Convergence"
    assert action["Status"] == "Do not interpret yet"
    guide = app.build_manuscript_claim_guide(result, diagnostics).set_index("ManuscriptArea")
    assert guide.loc["Convergence and final interpretability", "ClaimStatus"] == "Do not claim"
    audit = app.build_statistical_assumption_audit(result, diagnostics).set_index("Area")
    assert "MML fixed population prior SD" not in audit.index
    assert audit.loc["MML freely estimated population SD", "Status"] == "WITHHELD"


def test_final_readiness_contract_survives_standard_zip_as_frames_and_json():
    result = _result(converged=True)
    readiness_frame = app.build_final_report_readiness(
        result,
        _diagnostics(review_residuals=False),
        all_bias_results={},
    )
    frames = {"final_report_readiness": readiness_frame}
    app._add_final_readiness_contract_frames(frames, readiness_frame)
    assets = app.build_evidence_contract_text_assets(result, readiness_frame)
    payload = app._exports.build_tables_zip(frames, text_assets=assets)

    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        names = set(archive.namelist())
        contract_payload = json.loads(
            archive.read("final_report_readiness_evidence_contract.json")
        )
        identity_payload = json.loads(
            archive.read("final_report_readiness_analysis_identity.json")
        )

    assert "final_report_readiness_analysis_identity.csv" in names
    assert "final_report_readiness_evidence_records.csv" in names
    identity = evidence.AnalysisIdentity.from_payload(identity_payload)
    records = tuple(
        evidence.EvidenceRecord.from_payload(payload)
        for payload in contract_payload["evidence_records"]
    )
    evidence.validate_evidence_records(records, identity=identity)


def test_person_anchor_chain_uses_public_safe_aggregate_artifact():
    result = _result(converged=True)
    result["equating_chain"] = {
        "available": True,
        "edges": pd.DataFrame([{
            "From": "Current run",
            "To": "Anchor baseline via Person",
            "Facet": "Person",
            "CommonLinkedLevels": 2,
            "Strength": "Linked",
        }]),
    }
    readiness_frame = app.build_final_report_readiness(
        result,
        _diagnostics(review_residuals=False),
        all_bias_results={},
    ).set_index("Check")
    assert json.loads(
        readiness_frame.loc["Equating-chain summary", "SourceArtifactsJSON"]
    ) == ["equating_chain_readiness_summary.csv"]

    frames = {
        "equating_chain_edges": result["equating_chain"]["edges"],
        "equating_chain_readiness_summary": (
            app.build_equating_chain_readiness_summary(result)
        ),
    }
    public_frames = app.prepare_download_frames_for_privacy(
        frames,
        public_export_mode=True,
    )
    assert "equating_chain_edges" not in public_frames
    assert "equating_chain_readiness_summary" in public_frames
