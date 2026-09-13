"""Contracts for the Streamlit-free Monte Carlo orchestration layer."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mfrm_app import operating_characteristics as oc


PRECISION_PLAN = (
    Path(__file__).resolve().parents[1]
    / "validation"
    / "operating_characteristics_precision_plan_20260809.json"
)


def _conditions() -> list[dict[str, object]]:
    return [
        {
            "ConditionId": "balanced-null",
            "SeedGroup": "balanced-pair",
            "Design": "balanced",
            "TruthBias": 0.0,
        },
        {
            "ConditionId": "balanced-alt",
            "SeedGroup": "balanced-pair",
            "Design": "balanced",
            "TruthBias": 0.6,
        },
    ]


def test_manifest_is_deterministic_unique_and_records_common_random_numbers():
    first = oc.build_replicate_manifest(
        _conditions(),
        replicates=3,
        base_seed=20260809,
    )
    second = oc.build_replicate_manifest(
        _conditions(),
        replicates=3,
        base_seed=20260809,
    )

    pd.testing.assert_frame_equal(first, second)
    assert not first.duplicated(["ConditionId", "Replicate"]).any()
    assert first["RunId"].is_unique
    assert set(first["SeedCoupling"]) == {"common_random_numbers"}
    paired = first.pivot(index="Replicate", columns="ConditionId", values="Seed")
    assert (paired["balanced-null"] == paired["balanced-alt"]).all()
    assert oc.frame_fingerprint(first) == oc.frame_fingerprint(second)


def test_manifest_refuses_duplicate_conditions_and_implicit_budget_expansion():
    with pytest.raises(ValueError, match="unique"):
        oc.build_replicate_manifest(
            [_conditions()[0], _conditions()[0]],
            replicates=1,
            base_seed=1,
        )
    with pytest.raises(ValueError, match="max_runs"):
        oc.build_replicate_manifest(
            _conditions(),
            replicates=6,
            base_seed=1,
            max_runs=10,
        )


def test_manifest_extension_preserves_prior_rows_and_appends_replicate_indices():
    prior = oc.build_replicate_manifest(
        _conditions(),
        replicates=2,
        base_seed=20260809,
    )
    expanded = oc.build_replicate_manifest(
        _conditions(),
        replicates=4,
        base_seed=20260809,
    )

    audit = oc.audit_manifest_extension(prior, expanded)

    assert audit["Passed"] is True
    assert audit["PriorRuns"] == 4
    assert audit["ExpandedRuns"] == 8
    assert audit["AddedRuns"] == 4
    assert audit["PriorManifestSHA256"] == audit["RetainedRowsSHA256"]


def test_manifest_extension_rejects_changed_prior_seed():
    prior = oc.build_replicate_manifest(
        _conditions(),
        replicates=2,
        base_seed=20260809,
    )
    expanded = oc.build_replicate_manifest(
        _conditions(),
        replicates=4,
        base_seed=20260809,
    )
    prior_run = str(prior.iloc[0]["RunId"])
    expanded.loc[expanded["RunId"].eq(prior_run), "Seed"] += 1

    with pytest.raises(ValueError, match="changed one or more prior row values"):
        oc.audit_manifest_extension(prior, expanded)


def test_canonical_condition_id_is_order_invariant():
    left = oc.canonical_condition_id({"n": 30, "missing": 0.2})
    right = oc.canonical_condition_id({"missing": 0.2, "n": 30})

    assert left == right


def test_registered_precision_plan_attains_declared_binomial_targets():
    plan = oc.load_precision_plan(PRECISION_PLAN)

    assert plan["schema_version"] == oc.PRECISION_PLAN_SCHEMA_VERSION
    assert plan["promotion_gates"]["public_surface_status"] == "Withheld"
    assert plan["budget_rule"]["pilot_excluded_from_confirmatory_performance"] is True
    for entry in plan["confirmatory_precision"].values():
        n = int(entry["minimum_eligible_per_condition"])
        p = float(entry["reference_rate"])
        assert np.sqrt(p * (1 - p) / n) <= float(entry["target_mcse"])


def test_precision_plan_rejects_rounded_primary_decisions(tmp_path):
    plan = json.loads(PRECISION_PLAN.read_text(encoding="utf-8"))
    plan["authoritative_decisions"]["classification_values"] = "display_rounded"
    altered = tmp_path / "altered-plan.json"
    altered.write_text(json.dumps(plan), encoding="utf-8")

    with pytest.raises(ValueError, match="finite unrounded"):
        oc.load_precision_plan(altered)


def test_pilot_eligibility_planning_uses_lower_bound_and_fails_closed():
    complete = oc.plan_confirmatory_attempts(
        20,
        20,
        target_eligible=500,
        maximum_attempts=2500,
    )
    rare = oc.plan_confirmatory_attempts(
        1,
        20,
        target_eligible=500,
        maximum_attempts=2500,
    )
    absent = oc.plan_confirmatory_attempts(
        0,
        20,
        target_eligible=500,
        maximum_attempts=2500,
    )

    assert complete["RequiredFixedAttempts"] == 597
    assert complete["Status"] == "fixed_attempt_budget_available"
    assert rare["Status"] == "blocked_for_redesign_attempt_cap"
    assert absent["RequiredFixedAttempts"] is None
    assert absent["Status"] == "blocked_for_redesign_zero_lower_bound"


def test_binary_summary_keeps_failures_out_of_the_negative_decision_count():
    rows = pd.DataFrame({
        "ConditionId": ["null"] * 5,
        "Engine": ["Python"] * 5,
        "Estimator": ["JMLE"] * 5,
        "TruthPositive": [False] * 5,
        "AnalysisEligible": [True, True, True, False, True],
        "DecisionStrong": [True, False, False, False, pd.NA],
    })

    summary = oc.summarize_binary_operating_characteristics(
        rows,
        decision_columns=["DecisionStrong"],
    ).iloc[0]

    assert summary["Attempts"] == 5
    assert summary["EligibleDecisions"] == 3
    assert summary["UnavailableOrIneligible"] == 2
    assert summary["FalsePositive"] == 1
    assert summary["TrueNegative"] == 2
    assert summary["FalsePositiveRate"] == pytest.approx(1 / 3)
    assert np.isnan(summary["Power"])
    assert str(summary["EvidenceTier"]).startswith("pilot only")


def test_binary_summary_reports_power_separately_from_false_positive_rate():
    rows = pd.DataFrame({
        "ConditionId": ["alternative"] * 4,
        "Engine": ["Python"] * 4,
        "Estimator": ["JMLE"] * 4,
        "TruthPositive": [True] * 4,
        "AnalysisEligible": [True] * 4,
        "DecisionStrong": [True, True, False, True],
    })

    summary = oc.summarize_binary_operating_characteristics(
        rows,
        decision_columns=["DecisionStrong"],
    ).iloc[0]

    assert summary["Power"] == pytest.approx(0.75)
    assert summary["TruePositive"] == 3
    assert summary["FalseNegative"] == 1
    assert np.isnan(summary["FalsePositiveRate"])
    assert 0 <= summary["PowerWilsonLower95"] <= summary["PowerWilsonUpper95"] <= 1


def test_estimation_summary_matches_closed_form_and_does_not_impute_missing_se():
    rows = pd.DataFrame({
        "ConditionId": ["c1"] * 4,
        "Engine": ["Python"] * 4,
        "Estimator": ["JMLE"] * 4,
        "ParameterType": ["Rater"] * 4,
        "ErrorAligned": [-1.0, 0.0, 1.0, 100.0],
        "SE": [1.0, 1.0, np.nan, 1.0],
        "IncludedInSummary": [True, True, True, False],
    })

    summary = oc.summarize_estimation_operating_characteristics(rows).iloc[0]

    assert summary["RowsIncluded"] == 3
    assert summary["RowsUnavailable"] == 1
    assert summary["Bias"] == pytest.approx(0.0)
    assert summary["RMSE"] == pytest.approx(np.sqrt(2 / 3))
    assert summary["MAE"] == pytest.approx(2 / 3)
    assert summary["SEAvailableRate"] == pytest.approx(2 / 3)
    assert summary["CoverageN"] == 2
    assert summary["Coverage"] == pytest.approx(1.0)


def test_run_accounting_retains_failure_stage_and_reason():
    runs = pd.DataFrame({
        "ConditionId": ["c1", "c1", "c1"],
        "Engine": ["Python"] * 3,
        "Estimator": ["JMLE"] * 3,
        "FitReturned": [True, False, True],
        "Converged": [True, False, False],
        "InferenceReady": [True, False, False],
        "AnalysisEligible": [True, False, False],
        "FailureStage": ["", "fit", "bias"],
        "FailureReason": ["", "optimizer exception", "focal cell sparse"],
    })

    summary, reasons = oc.summarize_run_accounting(runs)
    row = summary.iloc[0]

    assert row["Attempts"] == 3
    assert row["FitsReturned"] == 2
    assert row["Converged"] == 1
    assert row["AnalysisEligible"] == 1
    assert set(reasons["FailureStage"]) == {"fit", "bias"}
    assert reasons["Count"].sum() == 2


def test_run_accounting_empty_failure_table_keeps_export_schema():
    runs = pd.DataFrame({
        "ConditionId": ["c1"],
        "Engine": ["Python"],
        "Estimator": ["JMLE"],
        "FitReturned": [True],
        "Converged": [True],
        "InferenceReady": [True],
        "AnalysisEligible": [True],
    })

    _, reasons = oc.summarize_run_accounting(runs)

    assert reasons.empty
    assert {
        "ConditionId", "Engine", "Estimator", "FailureStage",
        "FailureReason", "Count", "ShareOfAttempts",
    }.issubset(reasons.columns)


def test_conclusion_sensitivity_counts_only_comparable_rows():
    rows = pd.DataFrame({
        "ConditionId": ["c1"] * 4,
        "Engine": ["Python"] * 4,
        "Estimator": ["JMLE"] * 4,
        "RawDecision": [True, False, True, False],
        "RoundedDecision": [False, False, True, pd.NA],
    })

    audit = oc.audit_conclusion_sensitivity(
        rows,
        raw_column="RawDecision",
        comparison_columns=["RoundedDecision"],
    ).iloc[0]

    assert audit["ComparableRows"] == 3
    assert audit["UnavailableRows"] == 1
    assert audit["ChangedConclusions"] == 1
    assert audit["ChangedConclusionRate"] == pytest.approx(1 / 3)


def test_first_read_withholds_public_surface_and_labels_smoke_as_pilot():
    accounting = pd.DataFrame({
        "ConditionId": ["balanced", "sparse"],
        "Attempts": [2, 2],
        "Converged": [2, 2],
        "AnalysisEligible": [2, 0],
    })
    binary = pd.DataFrame({
        "EvidenceTier": ["pilot only (<100 eligible replicates)", "no eligible decisions"],
    })
    estimation = pd.DataFrame({"RowsIncluded": [12], "CoverageN": [10]})

    summary = oc.build_operating_characteristics_first_read(
        accounting,
        binary,
        estimation,
        profile="smoke",
    )

    evidence = summary.loc[summary["Check"].eq("Evidence scope")].iloc[0]
    public = summary.loc[summary["Check"].eq("Public application surface")].iloc[0]
    assert evidence["Status"] == "Pilot only"
    assert public["Status"] == "Withheld"
    assert "mfrmr/TAM/immer/sirt" in public["NextAction"]


def test_first_read_handles_missing_accounting_without_making_claims():
    summary = oc.build_operating_characteristics_first_read(
        pd.DataFrame(),
        pd.DataFrame(),
        pd.DataFrame(),
        profile="smoke",
    )

    assert summary.iloc[0]["Status"] == "Missing"
    assert "Do not infer" in summary.iloc[0]["DoNotClaim"]
