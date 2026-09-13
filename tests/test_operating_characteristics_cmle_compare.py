"""Contracts for native Python CMLE–immer CMLE normalization."""

from __future__ import annotations

import pandas as pd
import pytest

from validation import operating_characteristics_cmle_compare as compare


def _runs():
    python = pd.DataFrame(
        {
            "RunId": ["r1", "r2"],
            "ConditionId": ["c1", "c2"],
            "Design": ["balanced", "balanced"],
            "TruthBias": [0.0, 0.6],
            "Replicate": [1, 1],
            "ConditionalDesignEligible": [True, True],
            "FitReturned": [True, True],
            "Converged": [True, True],
            "InferenceReady": [True, True],
            "ParityEligible": [True, True],
            "RequestedConditionEligible": [True, True],
            "ConditionalLogLik": [-10.0, -20.0],
        }
    )
    immer = pd.DataFrame(
        {
            "RunId": ["r1", "r2"],
            "SharedConditionalDesignEligible": [True, True],
            "FitAttempted": [True, True],
            "FitReturned": [True, True],
            "Converged": [True, True],
            "InferenceReady": [True, False],
            "ParityEligible": [True, False],
            "RequestedConditionEligible": [True, False],
            "ConditionalLogLik": [-10.0 - 1e-12, -20.0],
            "GradientSupNorm": [2e-6, 2e-5],
        }
    )
    return python, immer


def test_pairs_keep_returned_numerical_identity_separate_from_readiness():
    python_runs, immer_runs = _runs()
    python_coefficients = pd.DataFrame(
        {
            "RunId": ["r1", "r2"],
            "Parameter": ["p1", "p1"],
            "Estimate": [0.2, 0.3],
            "SE": [0.1, 0.2],
        }
    )
    immer_coefficients = pd.DataFrame(
        {
            "RunId": ["r1", "r2"],
            "Parameter": ["p1", "p1"],
            "Estimate": [0.20000001, 0.30000002],
            "SE": [0.10000001, 0.20000002],
        }
    )

    pairs = compare.build_cmle_pairs(
        python_coefficients, immer_coefficients, python_runs, immer_runs
    )

    assert pairs["BothPresent"].tolist() == [True, True]
    assert pairs["JointParityEligible"].tolist() == [True, False]
    assert pairs.loc[1, "AbsEstimateDifference"] == pytest.approx(2e-8)
    assert pairs.loc[0, "ConditionalLogLikDifferencePythonMinusImmer"] == pytest.approx(1e-12)


def test_agreement_summary_does_not_discard_returned_threshold_failures():
    python_runs, immer_runs = _runs()
    coefficients = pd.DataFrame(
        {"RunId": ["r1", "r2"], "Parameter": ["p1", "p1"], "Estimate": [0.2, 0.3], "SE": [0.1, 0.2]}
    )
    pairs = compare.build_cmle_pairs(
        coefficients,
        coefficients.assign(Estimate=[0.200001, 0.300002]),
        python_runs,
        immer_runs,
    )

    summary = compare.build_agreement_summary(pairs)
    all_returned = summary.query("Scope == 'all_jointly_returned' and Design == 'All'").iloc[0]
    ready = summary.query("Scope == 'joint_inference_ready' and Design == 'All'").iloc[0]

    assert all_returned["RunIds"] == 2
    assert all_returned["FreeCoordinatePairs"] == 2
    assert ready["RunIds"] == 1
    assert ready["FreeCoordinatePairs"] == 1


def test_readiness_summary_reports_runid_and_unique_data_denominators():
    python_runs, immer_runs = _runs()
    # r1 and r2 represent clean/drift anchor conditions with identical rating
    # bytes; RunId counts must not masquerade as independent datasets.
    identities = pd.DataFrame({"RunId": ["r1", "r2"], "DataId": ["d1", "d1"]})

    summary = compare.build_readiness_summary(python_runs, immer_runs, identities)
    row = summary.query("Engine == 'immer' and Metric == 'InferenceReady'").iloc[0]

    assert row["RunIdNumerator"] == 1
    assert row["RunIdDenominator"] == 2
    assert row["UniqueDataNumerator"] == 1
    assert row["UniqueDataDenominator"] == 1
    assert row["UniqueDataRate"] == pytest.approx(1.0)


def test_omitted_bias_sensitivity_pairs_common_random_number_runs_without_detection_claim():
    python_runs = pd.DataFrame(
        {
            "RunId": ["null", "alt"],
            "Design": ["balanced", "balanced"],
            "TruthBias": [0.0, 0.6],
            "Replicate": [1, 1],
            "Seed": [99, 99],
            "ParityEligible": [True, True],
        }
    )
    coefficients = pd.DataFrame(
        {
            "RunId": ["null", "null", "alt", "alt"],
            "Parameter": [
                "facet:Rater:free:R01",
                "facet:Task:free:T01",
                "facet:Rater:free:R01",
                "facet:Task:free:T01",
            ],
            "Estimate": [0.2, 0.1, 0.05, 0.02],
            "SE": [0.1, 0.1, 0.1, 0.1],
        }
    )
    identities = pd.DataFrame(
        {"RunId": ["null", "alt"], "DataId": ["d-null", "d-alt"]}
    )

    pairs, summary = compare.build_omitted_bias_sensitivity(
        coefficients, python_runs, identities
    )

    assert pairs["PairEligible"].all()
    assert set(pairs["InterpretationBoundary"].str.contains("not a CMLE bias estimate")) == {True}
    row = summary.iloc[0]
    assert row["ReplicatePairs"] == 1
    assert row["MeanFocalRaterCoordinateChange"] == pytest.approx(-0.15)
    assert row["MeanFocalTaskCoordinateChange"] == pytest.approx(-0.08)
