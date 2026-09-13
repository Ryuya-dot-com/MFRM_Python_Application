"""Contracts for the sirt MML operating-characteristics sensitivity adapter."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from validation import operating_characteristics_sirt_compare as compare


def _sirt_runs() -> pd.DataFrame:
    rows = []
    for mode, shift in [(compare.Q30, 0.0), (compare.Q61, 0.01)]:
        for run_id, eligible in [("r1", True), ("r2", False)]:
            rows.append(
                {
                    "RunId": run_id,
                    "Mode": mode,
                    "ConditionId": f"c-{run_id}",
                    "Design": "sparse_missing",
                    "TruthBias": 0.0,
                    "Replicate": int(run_id[-1]),
                    "Seed": int(run_id[-1]),
                    "AnalysisEligible": eligible,
                    "PopulationSD": 1.0 + shift,
                    "LogLik": -10.0 + shift,
                    "EAPReliability": 0.5 + shift,
                }
            )
    return pd.DataFrame(rows)


def _quadrature_parameters() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    common = []
    for mode, shift in [(compare.Q30, 0.0), (compare.Q61, 0.001)]:
        for run_id, eligible in [("r1", True), ("r2", False)]:
            common.append(
                {
                    "RunId": run_id,
                    "Mode": mode,
                    "IncludedInSummary": eligible,
                    "EstimateAligned": 0.2 + shift,
                }
            )
    rater = pd.DataFrame(common).assign(Level="R01")
    item = pd.DataFrame(common).assign(
        Item="T01__C01",
        Estimate=lambda frame: frame["EstimateAligned"]
        + frame["Mode"].eq(compare.Q61).astype(float) * 0.049,
    )
    person = pd.DataFrame(common).rename(columns={"EstimateAligned": "Estimate"}).assign(
        Person="P01"
    )
    return rater, item, person


def test_quadrature_summary_separates_returned_and_analysis_eligible_scopes():
    rater, item, person = _quadrature_parameters()

    pairs, summary = compare.build_quadrature_pairs(
        rater, item, person, _sirt_runs()
    )

    rater_rows = summary.loc[
        summary["ParameterType"].eq("Rater severity")
        & summary["Design"].eq("sparse_missing")
    ].set_index("Scope")
    assert rater_rows.loc["all_jointly_returned", "RunIds"] == 2
    assert rater_rows.loc["jointly_analysis_eligible", "RunIds"] == 1
    assert rater_rows.loc["jointly_analysis_eligible", "ThresholdDenominator"] == 1

    item_rows = pairs.loc[pairs["Parameter"].eq("T01__C01")].set_index(
        "ParameterType"
    )
    assert item_rows.loc["Virtual-item centered location", "AbsDifference"].iloc[0] == pytest.approx(0.001)
    assert item_rows.loc["Virtual-item raw location", "AbsDifference"].iloc[0] == pytest.approx(0.05)


def test_quadrature_threshold_uses_unrounded_floating_point_value():
    rater, item, person = _quadrature_parameters()
    q61 = rater["Mode"].eq(compare.Q61)
    rater.loc[q61, "EstimateAligned"] = (
        rater.loc[q61, "EstimateAligned"] - 0.001
        + np.nextafter(0.001, np.inf)
    )

    _, summary = compare.build_quadrature_pairs(
        rater, item, person, _sirt_runs()
    )

    row = summary.loc[
        summary["ParameterType"].eq("Rater severity")
        & summary["Scope"].eq("jointly_analysis_eligible")
    ].iloc[0]
    assert row["MaxAbsDifference"] > 0.001
    assert row["RunsMaxDifferenceLE1e3"] == 0
    assert row["ThresholdDenominator"] == 1


def test_python_sirt_agreement_excludes_failed_rows_but_retains_them_as_evidence():
    python = pd.DataFrame(
        {
            "RunId": ["r1", "r2"],
            "Level": ["R01", "R01"],
            "ComparisonScale": ["mean_aligned_location"] * 2,
            "ConditionId": ["c1", "c2"],
            "Design": ["sparse_missing"] * 2,
            "TruthBias": [0.0, 0.0],
            "Replicate": [1, 2],
            "Seed": [1, 2],
            "Facet": ["Rater", "Rater"],
            "TruthAligned": [0.0, 0.0],
            "EstimateAligned": [0.1, 9.0],
            "SE": [0.1, 0.1],
            "ErrorAligned": [0.1, 9.0],
            "IncludedInSummary": [True, True],
        }
    )
    sirt = pd.DataFrame(
        {
            "RunId": ["r1", "r2"],
            "Level": ["R01", "R01"],
            "ComparisonScale": ["mean_aligned_location"] * 2,
            "Mode": [compare.Q61] * 2,
            "EstimateAligned": [0.2, -9.0],
            "SE": [0.1, 0.1],
            "ErrorAligned": [0.2, -9.0],
            "IncludedInSummary": [True, False],
            "Anchored": [False, False],
        }
    )

    pairs, summary = compare.build_python_sirt_rater_agreement(python, sirt)

    assert pairs["BothPresent"].sum() == 2
    assert pairs["BothIncluded"].sum() == 1
    row = summary.iloc[0]
    assert row["Pairs"] == 2
    assert row["BothIncludedPairs"] == 1
    assert row["MeanAbsDifference"] == pytest.approx(0.1)
    assert "not parity" in row["InterpretationBoundary"]


def test_anchor_sensitivity_recovers_injected_location_shift_without_robustness_claim():
    recovery_rows = []
    run_rows = []
    for design, estimate, deviance, sigma in [
        ("balanced_large_anchors", -0.2, 100.0, 1.0),
        ("anchor_drift", 0.05, 100.00001, 1.000001),
    ]:
        recovery_rows.append(
            {
                "RunId": design,
                "ConditionId": design,
                "Engine": "sirt",
                "Estimator": "MML",
                "Mode": compare.Q61,
                "Design": design,
                "TruthBias": 0.0,
                "Replicate": 1,
                "Seed": 99,
                "Level": "R01",
                "Truth": -0.2,
                "Estimate": estimate,
                "SE": 1e-5,
                "EstimateAligned": estimate,
                "TruthAligned": -0.2,
                "ErrorAligned": estimate + 0.2,
                "ComparisonScale": "anchor_identified_absolute",
                "Anchored": True,
                "IncludedInSummary": True,
                "CoverageEligible": False,
            }
        )
        run_rows.append(
            {
                "Mode": compare.Q61,
                "Design": design,
                "TruthBias": 0.0,
                "Replicate": 1,
                "Seed": 99,
                "Deviance": deviance,
                "PopulationSD": sigma,
            }
        )

    result = compare.build_anchor_sensitivity(
        pd.DataFrame(recovery_rows), pd.DataFrame(run_rows)
    )

    assert result.loc[0, "EstimateShiftDriftMinusClean"] == pytest.approx(0.25)
    assert result.loc[0, "ShiftErrorFromInjected0p25"] == pytest.approx(0.0)
    assert result.loc[0, "PopulationSDShiftDriftMinusClean"] == pytest.approx(1e-6)


def test_omitted_bias_pairs_require_both_runs_to_pass_the_analysis_gate():
    rows = []
    for truth_bias, included, estimate in [(0.0, True, 0.1), (0.6, False, 0.4)]:
        rows.append(
            {
                "RunId": f"r-{truth_bias}",
                "Mode": compare.Q61,
                "Design": "sparse_missing",
                "TruthBias": truth_bias,
                "Replicate": 1,
                "Seed": 7,
                "Level": "R01",
                "EstimateAligned": estimate,
                "IncludedInSummary": included,
            }
        )

    result = compare.build_omitted_bias_sensitivity(pd.DataFrame(rows))

    assert len(result) == 1
    assert not bool(result.loc[0, "PairEligible"])
    assert result.loc[0, "AbsEstimateChange"] == pytest.approx(0.3)
    assert "not bias detection" in result.loc[0, "InterpretationBoundary"]


def test_sparse_identification_keeps_assumption_based_mml_distinct_from_rank_evidence():
    sirt = pd.DataFrame(
        {"RunId": ["r1"], "Mode": [compare.Q61], "Design": ["sparse_missing"]}
    )
    python = pd.DataFrame(
        {
            "RunId": ["r1"], "Design": ["sparse_missing"],
            "FitReturned": [True], "Converged": [True], "InferenceReady": [True],
        }
    )
    cmle = pd.DataFrame(
        {
            "RunId": ["r1"], "Design": ["sparse_missing"],
            "ConditionalDesignEligible": [False], "ConditionalRank": [5],
            "ConditionalNullity": [7],
        }
    )
    mfrmr = pd.DataFrame(
        {
            "RunId": ["r1"], "Design": ["sparse_missing"],
            "Mode": ["MFRMR_JML_STRICT"], "FitReturned": [False],
            "Converged": [False], "InferenceReady": [False],
            "FailureReason": ["rank deficient"],
        }
    )

    result = compare.build_sparse_identification(sirt, python, cmle, mfrmr)

    assert bool(result.loc[0, "PythonJMLEFitReturned"])
    assert not bool(result.loc[0, "PythonCMLEConditionalDesignEligible"])
    assert not bool(result.loc[0, "MfrmrFitReturned"])
    assert "common Person distribution" in result.loc[0, "IdentificationBoundary"]


def test_first_read_projection_stays_withheld_and_orders_user_actions():
    runs = pd.DataFrame(
        {
            "Mode": [compare.Q30, compare.Q61],
            "FitReturned": [True, True],
            "AnalysisEligible": [True, True],
        }
    )
    quadrature = pd.DataFrame(
        {
            "ParameterType": [
                "Rater severity", "Virtual-item centered location",
                "Virtual-item raw location", "Person EAP",
            ],
            "AbsDifference": [0.001, 0.002, 0.05, 0.08],
            "BothIncluded": [True] * 4,
        }
    )
    agreement = pd.DataFrame(
        {
            "Design": ["balanced_small"],
            "MeanAbsDifference": [0.01],
            "MaxAbsDifference": [0.02],
        }
    )
    anchor = pd.DataFrame({"ShiftErrorFromInjected0p25": [1e-6]})
    sparse = pd.DataFrame({"AnalysisEligible": [True, False]})

    result = compare.build_first_read_summary(
        runs, quadrature, agreement, anchor, sparse
    )

    assert result["Order"].tolist() == list(range(1, 8))
    public = result.loc[result["Section"].eq("Public application surface")].iloc[0]
    assert public["Status"] == "Withheld"
    assert not bool(public["PublicSurfaceEnabled"])
    assert result["UserAction"].str.len().gt(0).all()
