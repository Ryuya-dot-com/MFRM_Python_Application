"""Tests for fixed population-prior-SD sensitivity diagnostics in MML."""

from __future__ import annotations

import pandas as pd

import streamlit_app as app


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


def test_mml_prior_sd_sensitivity_refuses_non_mml_result():
    out = app.evaluate_mml_prior_sd_sensitivity({"config": {"method": "JMLE"}})
    assert out["available"] is False
    assert "MML" in out["reason"]
    assert out["summary"].empty
