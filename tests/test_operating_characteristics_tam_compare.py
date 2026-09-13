"""Contracts for TAM MML operating-characteristics normalization."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from validation import operating_characteristics_tam_compare as compare


def _runs() -> pd.DataFrame:
    rows = []
    for mode, q_shift in [(compare.Q21, 0.0), (compare.Q61, 0.02)]:
        for run_id, design, eligible in [
            ("r1", "balanced_small", True),
            ("r2", "sparse_missing", mode == compare.Q61),
        ]:
            boundary = run_id == "r2" and mode == compare.Q21
            rows.append(
                {
                    "RunId": run_id,
                    "ConditionId": f"c-{run_id}",
                    "Design": design,
                    "TruthBias": 0.0,
                    "Replicate": int(run_id[-1]),
                    "Seed": int(run_id[-1]),
                    "Mode": mode,
                    "FitReturned": True,
                    "TerminalProgressParsed": True,
                    "AnalysisEligible": eligible,
                    "TerminalMaxParameterChange": 1e-5,
                    "PopulationVariance": 0.0010000001 if boundary else 1.0,
                    "MinimumPopulationVariance": 0.001,
                    "PopulationVarianceBoundaryTolerance": 1e-8,
                    "PopulationVarianceAtLowerBound": boundary,
                    "PopulationSD": np.sqrt(0.0010000001) if boundary else 1.0 + q_shift,
                    "LogLik": -10.0 + q_shift,
                    "EAPReliability": 0.5 + q_shift,
                }
            )
    return pd.DataFrame(rows)


def _quadrature_frames():
    rows = []
    for mode, shift in [(compare.Q21, 0.0), (compare.Q61, 0.002)]:
        for run_id, eligible in [("r1", True), ("r2", mode == compare.Q61)]:
            rows.append(
                {
                    "RunId": run_id,
                    "Mode": mode,
                    "IncludedInSummary": eligible,
                    "EstimateAligned": 0.2 + shift,
                }
            )
    facets = pd.DataFrame(rows).assign(Facet="Rater", Level="R01")
    steps = pd.DataFrame(rows).rename(columns={"EstimateAligned": "Estimate"}).assign(Step=1)
    persons = pd.DataFrame(rows).rename(columns={"EstimateAligned": "Estimate"}).assign(Person="P01")
    surfaces = pd.DataFrame(rows).rename(
        columns={"EstimateAligned": "CumulativeDifficulty"}
    ).assign(GeneralizedItem="C01-raterR01-taskT01", Category=1)
    return facets, steps, persons, surfaces


def test_quadrature_evidence_separates_returned_and_eligible_scopes():
    facets, steps, persons, surfaces = _quadrature_frames()

    pairs, summary, run_pairs = compare.build_quadrature_evidence(
        facets, steps, persons, surfaces, _runs()
    )

    rows = summary.loc[
        summary["ParameterType"].eq("Facet: Rater")
        & summary["Design"].eq("sparse_missing")
    ].set_index("Scope")
    assert rows.loc["all_jointly_returned", "RunIds"] == 1
    assert rows.loc["jointly_analysis_eligible", "RunIds"] == 0
    assert pairs["BothPresent"].all()
    assert len(run_pairs) == 2


def test_quadrature_threshold_uses_unrounded_values():
    facets, steps, persons, surfaces = _quadrature_frames()
    mask = facets["Mode"].eq(compare.Q61) & facets["RunId"].eq("r1")
    facets.loc[mask, "EstimateAligned"] = 0.2 + np.nextafter(0.001, np.inf)

    _, summary, _ = compare.build_quadrature_evidence(
        facets, steps, persons, surfaces, _runs()
    )

    row = summary.loc[
        summary["ParameterType"].eq("Facet: Rater")
        & summary["Design"].eq("balanced_small")
        & summary["Scope"].eq("jointly_analysis_eligible")
    ].iloc[0]
    assert row["MaxAbsDifference"] > 0.001
    assert row["RunsMaxDifferenceLE1e3"] == 0
    assert row["ThresholdBasis"] == "unrounded retained floating-point values"


def test_variance_boundary_audit_exposes_naive_float_miss():
    result = compare.build_variance_boundary_audit(_runs())
    row = result.loc[result["RunId"].eq("r2")].iloc[0]

    assert row["StoredOffsetFromLowerBoundQ21"] == pytest.approx(1e-10)
    assert not bool(row["NaiveExactBoundaryDecisionQ21"])
    assert bool(row["TolerantBoundaryDecisionQ21"])
    assert bool(row["NaiveVsTolerantMismatchQ21"])


def test_python_tam_agreement_keeps_jmle_mml_boundary_explicit():
    python = pd.DataFrame(
        {
            "RunId": ["r1", "r2"], "Facet": ["Rater", "Rater"],
            "Level": ["R01", "R01"], "ComparisonScale": ["mean_aligned_location"] * 2,
            "ConditionId": ["c1", "c2"], "Design": ["balanced_small", "sparse_missing"],
            "TruthBias": [0.0, 0.0], "Replicate": [1, 2], "Seed": [1, 2],
            "TruthAligned": [0.0, 0.0], "EstimateAligned": [0.1, 9.0],
            "SE": [0.1, 0.1], "ErrorAligned": [0.1, 9.0],
            "IncludedInSummary": [True, True],
        }
    )
    tam = pd.DataFrame(
        {
            "RunId": ["r1", "r2"], "Facet": ["Rater", "Rater"],
            "Level": ["R01", "R01"], "ComparisonScale": ["mean_aligned_location"] * 2,
            "Mode": [compare.Q61] * 2, "EstimateAligned": [0.2, -9.0],
            "SE": [0.1, 0.1], "ErrorAligned": [0.2, -9.0],
            "IncludedInSummary": [True, False], "Anchored": [False, False],
            "DerivedConstraint": [False, False], "CoverageEligible": [True, False],
        }
    )

    pairs, summary = compare.build_python_tam_agreement(python, tam)

    assert pairs["BothPresent"].sum() == 2
    assert pairs["BothIncluded"].sum() == 1
    balanced = summary.loc[summary["Design"].eq("balanced_small")].iloc[0]
    assert balanced["MeanAbsDifference"] == pytest.approx(0.1)
    assert "not parity" in balanced["InterpretationBoundary"]


def test_anchor_sensitivity_records_sum_zero_compensation():
    recovery_rows = []
    run_rows = []
    for design, values in [
        ("balanced_large_anchors", {"R01": -0.2, "R02": -0.1, "R03": 0.1, "R04": 0.2}),
        ("anchor_drift", {"R01": 0.05, "R02": 0.15, "R03": -0.15, "R04": -0.05}),
    ]:
        for level, estimate in values.items():
            recovery_rows.append(
                {
                    "RunId": design, "Mode": compare.Q61, "Design": design,
                    "TruthBias": 0.0, "Replicate": 1, "Seed": 99,
                    "Facet": "Rater", "Level": level, "Estimate": estimate,
                    "Anchored": level in {"R01", "R02"},
                }
            )
        run_rows.append(
            {
                "Mode": compare.Q61, "Design": design, "TruthBias": 0.0,
                "Replicate": 1, "Seed": 99, "Deviance": 100.0,
                "PopulationSD": 1.0,
            }
        )

    result = compare.build_anchor_sensitivity(
        pd.DataFrame(recovery_rows), pd.DataFrame(run_rows)
    )

    anchored = result.loc[result["AnchoredClean"]]
    free = result.loc[~result["AnchoredClean"]]
    assert anchored["EstimateShiftDriftMinusClean"].tolist() == pytest.approx([0.25, 0.25])
    assert free["EstimateShiftDriftMinusClean"].tolist() == pytest.approx([-0.25, -0.25])
    assert result["EstimateShiftDriftMinusClean"].sum() == pytest.approx(0.0)


def test_progress_threshold_summary_does_not_upgrade_formatted_text_to_gradient():
    result = compare.build_progress_threshold_summary(_runs())
    primary = result.loc[result["Mode"].eq(compare.Q61)].set_index("Threshold")

    assert primary.loc[1e-5, "RunNumerator"] == 2
    assert primary.loc[1e-6, "RunNumerator"] == 0
    assert "formatted text" in primary.loc[1e-6, "PrecisionBoundary"]


def test_first_read_projection_remains_withheld():
    q = pd.DataFrame(
        {
            "ParameterType": [
                "Facet: Rater", "Person EAP", "Cumulative response-surface difficulty",
            ],
            "AbsDifference": [0.02, 0.3, 0.6],
            "BothIncluded": [True, True, True],
        }
    )
    python_tam = pd.DataFrame(
        {
            "Design": ["balanced_small"], "Facet": ["Rater"],
            "MeanAbsDifference": [0.01], "MaxAbsDifference": [0.02],
        }
    )
    tam_sirt = pd.DataFrame(
        {
            "Design": ["balanced_small"], "MeanAbsDifference": [0.001],
        }
    )
    anchor = pd.DataFrame(
        {
            "ExpectedConstraintRole": [
                "fixed anchor receives +0.25",
                "unanchored Rater compensates under sum-zero constraint",
            ],
            "EstimateShiftDriftMinusClean": [0.25, -0.25],
        }
    )
    variance = pd.DataFrame(
        {
            "NaiveVsTolerantMismatchQ21": [True],
            "NaiveVsTolerantMismatchQ61": [False],
        }
    )
    sparse = pd.DataFrame({"AnalysisEligible": [True, True]})

    result = compare.build_first_read_summary(
        _runs(), q, python_tam, tam_sirt, anchor, variance, sparse
    )

    assert result["Order"].tolist() == list(range(1, 10))
    public = result.loc[result["Section"].eq("Public application surface")].iloc[0]
    assert public["Status"] == "Withheld"
    assert not bool(public["PublicSurfaceEnabled"])
