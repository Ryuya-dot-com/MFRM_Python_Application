"""S0 fail-closed contracts for statistically sensitive outputs."""

from __future__ import annotations

import numpy as np
import pandas as pd

from mfrm_app import output_qualification as qualification
import streamlit_app as app


def _records_by_output(method: str) -> dict[str, dict[str, object]]:
    records = qualification.records(
        qualification.model_choice_qualifications(method)
    )
    return {str(record["OutputID"]): record for record in records}


def test_jmle_model_choice_is_fail_closed_with_stable_reason_codes():
    records = _records_by_output("JMLE")

    assert records["model_choice.information_criteria"]["QualificationStatus"] == "TECHNICAL_ONLY"
    assert records["model_choice.automatic_recommendation"]["QualificationStatus"] == "WITHHELD"
    assert records["model_choice.automatic_recommendation"]["ReasonCode"] == (
        "model.model_001.jmle_auto_recommendation_withheld"
    )
    assert all(
        record["PublicConclusionAllowed"] is False
        for record in records.values()
    )


def test_every_gpcm_lrt_pair_is_withheld_from_ordinary_chi_square_use():
    records = _records_by_output("MML")
    for output_id in (
        "model_choice.lrt.pcm_gpcm",
        "model_choice.lrt.rsm_gpcm",
    ):
        assert records[output_id]["QualificationStatus"] == "WITHHELD"
        assert records[output_id]["ReasonCode"] == (
            "model.model_002.gpcm_lrt_uncalibrated"
        )
        assert records[output_id]["RoadmapIssueID"] == "MODEL-002"


def test_covariance_rank_deficiency_cannot_be_reported_with_a_caveat():
    item = qualification.covariance_qualification(
        status="regularized",
        rank=3,
        param_count=4,
    )

    assert item.status is qualification.OutputUse.WITHHELD
    assert item.reason_code is qualification.QualificationReason.COVARIANCE_RANK_DEFICIENT
    assert item.public_conclusion_allowed is False
    assert item.raw_technical_export_allowed is True


def test_full_rank_covariance_remains_technical_until_coverage_gate():
    item = qualification.covariance_qualification(
        status="ok",
        rank=4,
        param_count=4,
    )

    assert item.status is qualification.OutputUse.TECHNICAL_ONLY
    assert item.reason_code is qualification.QualificationReason.COVARIANCE_COVERAGE_UNVALIDATED
    assert item.next_gate == "G2"


def test_covariance_audit_exports_the_same_fail_closed_contract():
    eigvals = np.array([8.0, 3.0, 1.0, 1e-18])
    rotation, _ = np.linalg.qr(np.random.default_rng(20260804).normal(size=(4, 4)))
    hessian = rotation @ np.diag(eigvals) @ rotation.T
    covariance, regularized, rank = app._invert_information_matrix(hessian)

    audit = app.build_mml_covariance_audit(
        {
            "cov": covariance,
            "hessian": hessian,
            "status": "regularized",
            "detail": "rank-deficient fixture",
            "regularized": regularized,
            "rank": rank,
        },
        {"config": {"method": "MML", "model": "GPCM"}},
    )
    row = audit.iloc[0]

    assert row["ClaimStatus"] == "Do not claim"
    assert row["QualificationStatus"] == "WITHHELD"
    assert row["ReasonCode"] == "stat.stat_003.covariance_rank_deficient"
    assert bool(row["PublicConclusionAllowed"]) is False


def test_result_bundle_attaches_qualification_to_likelihood_values():
    result = {
        "config": {"model": "RSM", "method": "JMLE", "n_cat": 3},
        "summary": pd.DataFrame({
            "Metric": ["N", "KParams", "LogLik", "AIC", "BIC"],
            "Value": [100.0, 12.0, -50.0, 124.0, 155.0],
        }),
    }

    frames = app.build_result_bundle_frames(result, {})

    assert "output_qualification" in frames
    matrix = frames["output_qualification"]
    assert set(matrix["SchemaVersion"]) == {"mfrm_output_qualification_v1"}
    assert not matrix["PublicConclusionAllowed"].astype(bool).any()

    likelihood = frames["likelihood_information_criteria"]
    assert set(likelihood["QualificationStatus"]) == {"TECHNICAL_ONLY"}
    assert set(likelihood["ReasonCode"]) == {
        "model.model_001.selection_scope_unvalidated"
    }


def test_full_download_and_apa_tables_carry_the_same_qualification():
    result = {
        "config": {
            "model": "RSM",
            "method": "JMLE",
            "facet_names": ["Rater"],
            "n_cat": 3,
        },
        "prep": {"n_obs": 100, "n_person": 20},
        "summary": pd.DataFrame({
            "Metric": ["N", "KParams", "LogLik", "AIC", "BIC"],
            "Value": [100.0, 12.0, -50.0, 124.0, 155.0],
        }),
    }
    diagnostics: dict[str, object] = {}

    frames, _ = app.collect_download_frames(
        result,
        diagnostics,
        {},
        pd.DataFrame(),
        pd.DataFrame(),
        public_export_mode=True,
    )
    assert "output_qualification" in frames
    assert set(frames["likelihood_information_criteria"]["QualificationStatus"]) == {
        "TECHNICAL_ONLY"
    }

    apa_tables = app._collect_apa_exportable_tables(result, diagnostics)
    assert "Output qualification" in apa_tables
    assert set(apa_tables["Estimation summary"]["PublicConclusionAllowed"]) == {
        False
    }


def test_pairwise_bias_contract_is_withheld_until_g3():
    item = qualification.bias_pairwise_qualification()

    assert item.status is qualification.OutputUse.WITHHELD
    assert item.roadmap_issue_id == "BIAS-001;BIAS-002"
    assert item.next_gate == "G3"
