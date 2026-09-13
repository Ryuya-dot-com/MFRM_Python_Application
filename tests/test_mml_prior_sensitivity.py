"""Tests for fixed population-prior-SD sensitivity diagnostics in MML."""

from __future__ import annotations

import copy
import io
import json
import zipfile

import pandas as pd

import streamlit_app as app
from mfrm_app import evidence


def _small_mml_response_data() -> pd.DataFrame:
    rows = []
    for pi, person in enumerate([f"P{i}" for i in range(1, 7)]):
        for ri, rater in enumerate(["R1", "R2"]):
            for ti, task in enumerate(["T1", "T2"]):
                rows.append({
                    "Person": person,
                    "Rater": rater,
                    "Task": task,
                    "Score": (pi + ri + ti) % 3,
                })
    return pd.DataFrame(rows)


def _small_mml_fit() -> dict:
    return app.mfrm_estimate(
        _small_mml_response_data(),
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="RSM",
        method="MML",
        quad_points=5,
        population_prior_sd=1.0,
        maxit=20,
        reltol=1e-3,
        mml_engine="EM",
    )


def test_mml_prior_sensitivity_plan_names_comparisons_and_decision_rule():
    result = {"config": {"method": "MML", "population_prior_sd": 1.0}}
    plan = app.build_mml_prior_sensitivity_plan(result, multipliers=(0.75, 1.0, 1.25))
    assert not plan.empty
    assert {"PrimaryComparisons", "DecisionRule", "RunStatus"}.issubset(plan.columns)
    assert "latent variance is not estimated" in plan["Boundary"].iloc[0]
    assert plan.loc[plan["Multiplier"] == 1.0, "RunStatus"].iloc[0] == "Already fit"


def test_evaluate_mml_prior_sd_sensitivity_refits_and_reports_deltas():
    result = _small_mml_fit()
    bundle = app.evaluate_mml_prior_sd_sensitivity(
        result,
        multipliers=(1.0, 1.25),
        maxit=10,
        reltol=1e-3,
    )
    assert bundle["available"] is True
    summary = bundle["summary"]
    assert {1.0, 1.25}.issubset(set(summary["PopulationPriorSD"].round(8)))
    assert summary["RunOK"].all()
    assert {"LogLik", "AIC", "BIC", "PersonEstimateSD"}.issubset(summary.columns)

    deltas = bundle["measure_deltas"]
    assert not deltas.empty
    assert {"Facet", "RMSE", "MaxAbsDifference", "RankCorrelation", "Status"}.issubset(deltas.columns)
    assert {"Person", "Rater", "Task"}.issubset(set(deltas["Facet"]))

    report = bundle["report"]
    assert not report.empty
    assert report["Area"].iloc[0] == "MML fixed population prior SD sensitivity"
    assert report["Status"].iloc[0] in {"Stable screen", "Review"}
    assert "fixed-population-prior-SD" in report["Interpretation"].iloc[0]

    contract = bundle["evidence_contract"]
    assert contract["schema_version"] == "mfrm_mml_prior_sd_contract_bundle_v1"
    plan = evidence.SensitivityPlan.from_payload(contract["sensitivity_plan"])
    assert "baseline" not in plan.required_variant_ids
    assert len(plan.required_variant_ids) == 1
    records = tuple(
        evidence.SensitivityRecord.from_payload(payload)
        for payload in contract["sensitivity_records"]
    )
    decision = evidence.SensitivityDecision.from_payload(
        contract["sensitivity_decision"],
        plan=plan,
        records=records,
    )
    readiness_identity = app.build_result_analysis_identity(result)
    assert plan.baseline_analysis_id == readiness_identity.analysis_id
    assert decision.baseline_analysis_id == readiness_identity.analysis_id
    assert decision.decision_id == bundle["settings"]["sensitivity_decision_id"]
    assert report["StabilityState"].iloc[0] == decision.stability_state.value
    assert set(bundle["plan"]["SensitivityPlanID"]) == {plan.plan_id}

    readiness = app.build_final_report_readiness(
        result,
        {
            "obs": pd.DataFrame({"StdResidual": [0.0, 0.1]}),
            "reliability": pd.DataFrame(),
            "pca_enabled": False,
        },
        all_bias_results={},
    )
    prior_row = readiness.loc[
        readiness["Check"] == "MML fixed population prior SD"
    ].iloc[0]
    assert prior_row["StabilityState"] == decision.stability_state.value
    assert prior_row["Status"] == (
        "Ready" if decision.stability_state is evidence.StabilityState.STABLE else "Review"
    )
    assert decision.decision_id in prior_row["ObservedJSON"]
    evidence_assets = app.build_evidence_contract_text_assets(result, readiness)
    assert "mml_prior_sd_sensitivity_evidence_contract.json" in evidence_assets
    assert "mml_prior_sd_sensitivity_settings.json" in evidence_assets
    export_frames: dict[str, pd.DataFrame] = {}
    app._add_mml_prior_sensitivity_export_frames(export_frames, result)
    archive_payload = app._exports.build_tables_zip(
        export_frames,
        text_assets=evidence_assets,
    )
    with zipfile.ZipFile(io.BytesIO(archive_payload)) as archive:
        archive_names = set(archive.namelist())
        archived_contract = json.loads(
            archive.read("mml_prior_sd_sensitivity_evidence_contract.json")
        )
    linked_artifacts = {
        artifact
        for record in (
            archived_contract["evidence_records"]
            + archived_contract["sensitivity_records"]
        )
        for artifact in record["source_artifacts"]
    }
    assert linked_artifacts.issubset(archive_names)

    quick_diagnostics = app.mfrm_diagnostics(
        result,
        compute_pca=False,
        compute_marginal=False,
    )
    quick_frames = app.build_result_bundle_frames(
        result,
        quick_diagnostics,
        all_bias_results={},
    )
    quick_readiness = quick_frames["final_report_readiness"]
    quick_assets = app.build_evidence_contract_text_assets(
        result,
        quick_readiness,
    )
    quick_archive_payload = app._exports.build_tables_zip(
        quick_frames,
        text_assets=quick_assets,
    )
    with zipfile.ZipFile(io.BytesIO(quick_archive_payload)) as archive:
        quick_archive_names = set(archive.namelist())
    readiness_artifacts = {
        artifact
        for _, row in quick_readiness.loc[
            quick_readiness["ComputationState"].isin(["AVAILABLE", "CAUTION"])
        ].iterrows()
        for artifact in json.loads(row["SourceArtifactsJSON"])
    }
    assert readiness_artifacts.issubset(quick_archive_names)
    public_frames = app.prepare_download_frames_for_privacy(
        quick_frames,
        public_export_mode=True,
    )
    public_archive_payload = app._exports.build_tables_zip(
        public_frames,
        text_assets=quick_assets,
    )
    with zipfile.ZipFile(io.BytesIO(public_archive_payload)) as archive:
        public_archive_names = set(archive.namelist())
    assert "residuals.csv" not in public_archive_names
    assert "global_residual_fit_summary.csv" in public_archive_names
    assert readiness_artifacts.issubset(public_archive_names)

    valid_contract = bundle["evidence_contract"]
    corrupted_contract = copy.deepcopy(valid_contract)
    corrupted_contract["sensitivity_decision"]["summary"] = "tampered summary"
    result["evidence_bundles"]["mml_prior_sd"][
        "evidence_contract"
    ] = corrupted_contract
    corrupted_readiness = app.build_final_report_readiness(
        result,
        {
            "obs": pd.DataFrame({"StdResidual": [0.0, 0.1]}),
            "reliability": pd.DataFrame(),
            "pca_enabled": False,
        },
        all_bias_results={},
    )
    corrupted_row = corrupted_readiness.loc[
        corrupted_readiness["Check"] == "MML fixed population prior SD"
    ].iloc[0]
    assert corrupted_row["StabilityState"] == "NOT_ASSESSED"
    assert "failed contract validation" in corrupted_row["Evidence"]
    corrupted_assets = app.build_evidence_contract_text_assets(
        result,
        corrupted_readiness,
    )
    assert "mml_prior_sd_sensitivity_evidence_contract.json" not in corrupted_assets
    result["evidence_bundles"]["mml_prior_sd"][
        "evidence_contract"
    ] = valid_contract

    baseline_only = app.evaluate_mml_prior_sd_sensitivity(
        result,
        prior_sds=[result["config"]["population_prior_sd"]],
    )
    assert baseline_only["available"] is False
    assert "non-baseline" in baseline_only["reason"]


def test_mml_prior_sd_sensitivity_refuses_non_mml_result():
    out = app.evaluate_mml_prior_sd_sensitivity({"config": {"method": "JMLE"}})
    assert out["available"] is False
    assert "MML" in out["reason"]
    assert out["summary"].empty
    assert out["evidence_contract"] is None


def test_mml_prior_sd_sensitivity_refuses_free_population_sd_result():
    out = app.evaluate_mml_prior_sd_sensitivity({
        "config": {
            "method": "MML",
            "estimate_population_sd": True,
            "population_prior_sd": 1.0,
        }
    })
    assert out["available"] is False
    assert "not applicable" in out["reason"]
    assert out["evidence_contract"] is None


def test_mml_prior_plan_fails_closed_for_invalid_population_sd():
    plan = app.build_mml_prior_sensitivity_plan({
        "config": {"method": "MML", "population_prior_sd": "bad"}
    })

    assert plan.empty
    assert "no valid" in plan.attrs["contract_reason"]


def test_mml_prior_comparison_coverage_rejects_missing_and_duplicate_levels():
    baseline = {
        "facets": {
            "person": pd.DataFrame(
                {"Person": ["P1", "P2", "P3"], "Estimate": [0.0, 0.2, 0.4]}
            ),
            "others": pd.DataFrame(
                {"Facet": ["Rater", "Rater"], "Level": ["R1", "R2"], "Estimate": [0.1, -0.1]}
            ),
        }
    }
    missing = {
        "facets": {
            "person": pd.DataFrame(
                {"Person": ["P1", "P2"], "Estimate": [0.0, 0.2]}
            ),
            "others": baseline["facets"]["others"].copy(),
        }
    }
    duplicate = {
        "facets": {
            "person": pd.DataFrame(
                {
                    "Person": ["P1", "P2", "P3", "P3"],
                    "Estimate": [0.0, 0.2, 0.4, 0.4],
                }
            ),
            "others": baseline["facets"]["others"].copy(),
        }
    }

    missing_audit = app._mml_prior_comparison_coverage(
        baseline, missing, latent_regression=False
    )
    duplicate_audit = app._mml_prior_comparison_coverage(
        baseline, duplicate, latent_regression=False
    )

    assert missing_audit["complete"] is False
    assert missing_audit["metrics"]["missing_measure_keys"] == 1
    assert duplicate_audit["complete"] is False
    assert duplicate_audit["metrics"]["variant_measure_duplicate_rows"] == 2


def test_mml_prior_evaluation_cannot_be_stable_with_incomplete_refit(monkeypatch):
    result = _small_mml_fit()
    original_estimate = app.mfrm_estimate

    def incomplete_refit(**kwargs):
        refit = original_estimate(**kwargs)
        refit["facets"]["person"] = refit["facets"]["person"].iloc[:-1].copy()
        return refit

    monkeypatch.setattr(app, "mfrm_estimate", incomplete_refit)
    bundle = app.evaluate_mml_prior_sd_sensitivity(
        result,
        prior_sds=(1.0, 1.25),
        maxit=10,
        reltol=1e-3,
    )

    decision = bundle["evidence_contract"]["sensitivity_decision"]
    record = bundle["evidence_contract"]["sensitivity_records"][0]
    assert decision["stability_state"] == "NOT_ASSESSED"
    assert record["comparable"] is False
    assert record["reason_code"] == "evidence.noncomparable"
    assert record["metrics"]["missing_measure_keys"] == 1


def test_mml_prior_refit_attestation_requires_exact_false_execution_flags():
    result = _small_mml_fit()
    controls = app._mml_prior_resolved_refit_controls(result)

    for config_override in (
        {"estimate_population_sd": "true"},
        {"compute_plausible_values": True},
        {"positive_facets": ["Rater"]},
    ):
        run = {
            **result,
            "config": {**result["config"], **config_override},
        }
        attestation = app._mml_prior_refit_attestation(
            result,
            run,
            population_prior_sd=result["config"]["population_prior_sd"],
            refit_controls=controls,
        )
        assert attestation["complete"] is False
        assert attestation["metrics"]["execution_attestation_problem_count"] >= 1
    assert attestation["metrics"]["semantic_invariants_match"] is False


def test_mml_prior_refit_attestation_compares_normalized_not_raw_upload_fingerprint():
    result = _small_mml_fit()
    result["config"]["input_data_fingerprint"] = "raw_upload_sha256:fixture"

    bundle = app.evaluate_mml_prior_sd_sensitivity(
        result,
        prior_sds=(1.0, 1.25),
        maxit=10,
        reltol=1e-3,
    )

    record = bundle["evidence_contract"]["sensitivity_records"][0]
    assert record["metrics"]["execution_attested"] is True
    assert record["metrics"]["execution_attestation_problem_count"] == 0
