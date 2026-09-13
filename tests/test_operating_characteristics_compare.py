"""Contracts for Python–mfrmr operating-characteristics normalization."""

from __future__ import annotations

import pandas as pd
import pytest

from validation import operating_characteristics_compare as compare


def test_parameter_agreement_separates_modes_and_readiness_inclusion():
    python = pd.DataFrame({
        "RunId": ["r1", "r1"],
        "ConditionId": ["c1", "c1"],
        "Design": ["balanced", "balanced"],
        "TruthBias": [0.0, 0.0],
        "Facet": ["Rater", "Rater"],
        "Level": ["R01", "R02"],
        "ComparisonScale": ["mean_aligned_location"] * 2,
        "EstimateAligned": [0.2, -0.2],
        "SE": [0.1, 0.1],
        "IncludedInSummary": [True, True],
    })
    mfrmr = pd.concat([
        pd.DataFrame({
            **{column: python[column] for column in [
                "RunId", "ConditionId", "Design", "TruthBias", "Facet", "Level", "ComparisonScale",
            ]},
            "Mode": "strict",
            "EstimateAligned": [0.201, -0.201],
            "SE": [0.101, 0.101],
            "IncludedInSummary": [True, True],
        }),
        pd.DataFrame({
            **{column: python[column] for column in [
                "RunId", "ConditionId", "Design", "TruthBias", "Facet", "Level", "ComparisonScale",
            ]},
            "Mode": "matched",
            "EstimateAligned": [0.21, -0.21],
            "SE": [0.11, 0.11],
            "IncludedInSummary": [False, False],
        }),
    ], ignore_index=True)

    pairs, summary = compare.build_parameter_agreement(python, mfrmr)

    assert len(pairs) == 4
    strict = summary.loc[summary["Mode"].eq("strict")].iloc[0]
    matched = summary.loc[summary["Mode"].eq("matched")].iloc[0]
    assert strict["BothIncludedPairs"] == 2
    assert strict["MaxAbsDifference"] == pytest.approx(0.001)
    assert matched["BothIncludedPairs"] == 0
    assert matched["MAEDifference"] == pytest.approx(0.01)


def test_bias_agreement_counts_decisions_only_when_both_engines_are_eligible():
    keys = {
        "RunId": ["r1", "r2"],
        "ConditionId": ["c1", "c2"],
        "Design": ["balanced", "sparse"],
        "TruthBias": [0.6, 0.6],
        "TruthPositive": [True, True],
        "Replicate": [1, 1],
        "Seed": [11, 12],
    }
    python = pd.DataFrame({
        **keys,
        "AnalysisEligible": [True, False],
        "BiasEstimate": [0.6, 0.8],
        "BiasSE": [0.1, 0.5],
        "p_holm": [0.001, 1.0],
        "AbsBias": [0.6, 0.8],
        "SparseCell": [False, True],
        "DecisionStrongRaw": [True, pd.NA],
        "DecisionStrongDisplayed": [True, pd.NA],
    })
    mfrmr = pd.DataFrame({
        **keys,
        "Mode": ["strict", "strict"],
        "AnalysisEligible": [True, False],
        "BiasEstimate": [0.6001, pd.NA],
        "BiasSE": [0.1001, pd.NA],
        "p_holm": [0.001, pd.NA],
        "AbsBias": [0.6001, pd.NA],
        "SparseCell": [False, True],
        "DecisionStrongRaw": [True, pd.NA],
        "DecisionStrongDisplayed": [True, pd.NA],
    })

    _, summary = compare.build_bias_agreement(python, mfrmr)
    row = summary.iloc[0]

    assert row["Pairs"] == 2
    assert row["BothEligiblePairs"] == 1
    assert row["DecisionComparablePairs"] == 1
    assert row["StrongDecisionAgreements"] == 1
    assert row["StrongDecisionDisagreements"] == 0


def test_readiness_contrast_retains_structural_rejection_and_gradient_review():
    python = pd.DataFrame({
        "RunId": ["r1", "r2"],
        "FitReturned": [True, True],
        "Converged": [True, True],
        "InferenceReady": [True, True],
        "AnalysisEligible": [True, False],
    })
    mfrmr = pd.DataFrame({
        "RunId": ["r1", "r2", "r1", "r2"],
        "Mode": ["strict", "strict", "matched", "matched"],
        "FitReturned": [True, False, True, False],
        "Converged": [True, False, True, False],
        "InferenceReady": [True, False, False, False],
        "AnalysisEligible": [True, False, False, False],
        "FailureStage": ["", "fit", "readiness", "fit"],
        "FailureReason": ["", "rank deficient", "gradient review", "rank deficient"],
        "GradientNorm": [0.00004, pd.NA, 0.08, pd.NA],
    })

    summary = compare.build_readiness_contrast(python, mfrmr)

    strict = summary.loc[summary["Mode"].eq("strict")].set_index("Contrast")["Runs"]
    matched = summary.loc[summary["Mode"].eq("matched")].set_index("Contrast")["Runs"]
    assert strict["both inference ready"] == 1
    assert strict["Python ready; mfrmr structurally rejected/failed"] == 1
    assert matched["Python ready; mfrmr returned but readiness withheld"] == 1
    assert matched["Python ready; mfrmr structurally rejected/failed"] == 1
