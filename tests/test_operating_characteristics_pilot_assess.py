"""Contracts for the frozen Python pilot assessment layer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mfrm_app import operating_characteristics as oc
from validation import operating_characteristics_pilot_assess as assess


PLAN = (
    Path(__file__).resolve().parents[1]
    / "validation"
    / "operating_characteristics_precision_plan_20260809.json"
)


def _runs(eligible: int, *, condition: str, design: str) -> pd.DataFrame:
    rows = []
    for replicate in range(1, 21):
        is_eligible = replicate <= eligible
        rows.append({
            "ConditionId": condition,
            "Design": design,
            "TruthBias": 0.0,
            "RunId": f"{condition}::rep-{replicate:05d}",
            "FitReturned": True,
            "Converged": True,
            "InferenceReady": True,
            "AnalysisEligible": is_eligible,
            "StrictAnalysisEligible1e4": False,
            "TerminalGradientSupNorm": 2e-4,
            "ElapsedSeconds": 0.5,
        })
    return pd.DataFrame(rows)


def test_budget_assessment_inflates_from_eligibility_lower_bound_and_caps_sparse():
    runs = pd.concat([
        _runs(20, condition="balanced", design="balanced_small"),
        _runs(4, condition="sparse", design="sparse_missing"),
    ], ignore_index=True)

    result = assess.build_budget_assessment(runs, oc.load_precision_plan(PLAN))
    balanced = result.loc[result["ConditionId"].eq("balanced")].iloc[0]
    sparse = result.loc[result["ConditionId"].eq("sparse")].iloc[0]

    assert balanced["RequiredFixedAttemptsFromRegisteredEligibility"] == 597
    assert balanced["RegisteredBudgetStatus"] == "fixed_attempt_budget_available"
    assert sparse["RequiredFixedAttemptsFromRegisteredEligibility"] == 6200
    assert sparse["RegisteredBudgetStatus"] == "blocked_for_redesign_attempt_cap"
    assert set(result["PlanningDisposition"]) == {
        "blocked_pending_strict_jmle_numerical_qualification"
    }


def test_gradient_readiness_keeps_all_three_thresholds_visible():
    runs = _runs(20, condition="balanced", design="balanced_small")
    runs.loc[0, "TerminalGradientSupNorm"] = 5e-7

    result = assess.build_gradient_readiness(runs)
    condition = result.loc[result["ConditionId"].eq("balanced")]

    assert set(condition["Threshold"]) == {1e-4, 1e-5, 1e-6}
    assert set(condition["GradientReady"]) == {1}
    assert len(result.loc[result["ConditionId"].eq("ALL")]) == 3


def test_anchor_sensitivity_preserves_fixed_and_free_level_contrasts():
    rows = []
    for design, anchored_estimate, free_estimate in (
        ("balanced_large_anchors", 0.1, -0.2),
        ("anchor_drift", 0.35, -0.19),
    ):
        rows.extend([
            {
                "Design": design,
                "TruthBias": 0.0,
                "Replicate": 1,
                "Seed": 1,
                "Facet": "Rater",
                "Level": "R01",
                "Truth": 0.1,
                "Estimate": anchored_estimate,
                "Anchored": True,
            },
            {
                "Design": design,
                "TruthBias": 0.0,
                "Replicate": 1,
                "Seed": 1,
                "Facet": "Rater",
                "Level": "R03",
                "Truth": -0.2,
                "Estimate": free_estimate,
                "Anchored": False,
            },
        ])

    result = assess.build_anchor_sensitivity(pd.DataFrame(rows))

    anchored = result.loc[result["Anchored"]].iloc[0]
    free = result.loc[~result["Anchored"]].iloc[0]
    assert np.isclose(anchored["EstimateShiftDriftMinusClean"], 0.25)
    assert np.isclose(free["EstimateShiftDriftMinusClean"], 0.01)
    assert np.isclose(anchored["ShiftErrorAgainstInput"], 0.0)


def test_first_read_never_enables_public_surface_and_names_sparse_denominators():
    runs = pd.concat([
        _runs(20, condition="balanced", design="balanced_small"),
        _runs(4, condition="sparse-null", design="sparse_missing"),
        _runs(4, condition="sparse-alt", design="sparse_missing").assign(TruthBias=0.6),
    ], ignore_index=True)
    budget = assess.build_budget_assessment(runs, oc.load_precision_plan(PLAN))
    checks = pd.DataFrame({"Check": ["identity"], "Passed": [True], "Evidence": ["ok"]})
    fit_audit = pd.DataFrame({"DisplayDecisionConsistent": [True, False]})
    anchor = pd.DataFrame({
        "Anchored": [True, False],
        "EstimateShiftDriftMinusClean": [0.25, 0.0],
    })

    result = assess.build_first_read(checks, runs, budget, fit_audit, anchor)
    sparse = result.loc[result["Check"].eq("Sparse design")].iloc[0]

    assert "4/20 and 4/20" in sparse["Evidence"]
    assert not result["PublicSurfaceEnabled"].any()
    assert result.iloc[-1]["Status"] == "Withheld"
