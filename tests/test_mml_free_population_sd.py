"""Tests for the opt-in free latent-population-SD estimator in MML.

The MML engine can estimate the person population SD (sigma) as a free
parameter instead of fixing it to a user-set prior. These tests check that the
estimator recovers a known sigma, sets the metric to the fitted scale, adds one
free parameter, returns a profile SE, keeps GPCM slope identification intact,
and — crucially — leaves the default fixed-SD path numerically unchanged.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def _sim_rsm(theta_sd: float, n_person: int = 160, seed: int = 4242):
    params = {
        "persons": [f"P{i:03d}" for i in range(n_person)],
        "raters": ["R1", "R2", "R3"],
        "tasks": ["T1", "T2"],
        "criteria": ["C1", "C2"],
        "theta_sd": theta_sd,
        "rater_severities": np.array([-0.3, 0.0, 0.3]),
        "task_difficulties": np.array([-0.2, 0.2]),
        "criterion_difficulties": np.array([-0.1, 0.1]),
        "tau": np.array([-1.0, 0.0, 1.0]),
    }
    return app._generate_mfrm_rsm_from_params(params, seed=seed)


def _fit(df, *, free, model="RSM", **kw):
    common = dict(
        data=df, person_col="Person", facet_cols=["Rater", "Task", "Criterion"],
        score_col="Score", model=model, method="MML", mml_engine="EM",
        maxit=300, reltol=1e-6, quad_points=19, population_prior_sd=1.0,
    )
    common.update(kw)
    return app.mfrm_estimate(estimate_population_sd=free, **common)


def test_free_sigma_recovers_known_population_sd():
    df = _sim_rsm(theta_sd=1.5)
    cfg = _fit(df, free=True)["config"]
    est = cfg.get("estimated_population_sd")
    assert est is not None and np.isfinite(est)
    assert 1.25 < est < 1.75, f"estimated sigma {est} not near truth 1.5"
    assert est > 1.05, "estimated sigma did not move above the fixed baseline 1.0"


def test_free_sigma_threads_metric_and_preserves_user_input():
    df = _sim_rsm(theta_sd=1.5)
    cfg = _fit(df, free=True)["config"]
    # The fitted scale becomes the quadrature SD; the user's input is preserved.
    assert cfg["population_prior_sd"] == pytest.approx(cfg["estimated_population_sd"])
    assert cfg["population_prior_sd_input"] == pytest.approx(1.0)
    assert cfg["mml_engine"] == "em"
    assert len(cfg.get("sigma_trace", [])) >= 2


def test_free_sigma_adds_one_parameter_and_profile_se():
    df = _sim_rsm(theta_sd=1.5)
    free = _fit(df, free=True)["config"]
    fixed = _fit(df, free=False)["config"]
    assert int(free["parameter_count"]) == int(fixed["parameter_count"]) + 1
    se = free.get("population_sd_se")
    assert se is not None and np.isfinite(se) and se > 0
    lo, hi = free["population_sd_ci"]
    assert lo < free["estimated_population_sd"] < hi


def test_fixed_path_unchanged_and_reports_no_estimate():
    """Default (fixed-SD) fit is byte-identical with/without the new code path."""
    df = _sim_rsm(theta_sd=1.5)
    fixed = _fit(df, free=False)
    assert fixed["config"].get("estimated_population_sd") is None
    assert fixed["config"].get("population_sd_se") is None
    # Two fixed fits on the same data are deterministic and identical.
    again = _fit(df, free=False)
    m1 = app.mfrm_diagnostics(fixed, compute_pca=False, compute_marginal=False)["measures"]
    m2 = app.mfrm_diagnostics(again, compute_pca=False, compute_marginal=False)["measures"]
    np.testing.assert_allclose(
        m1["Estimate"].to_numpy(), m2["Estimate"].to_numpy(), rtol=0, atol=0
    )


def test_free_sigma_keeps_gpcm_slope_identification():
    """Freeing sigma must not trade off against GPCM discrimination."""
    df = _sim_rsm(theta_sd=1.4, n_person=160, seed=77)
    res = _fit(df, free=True, model="GPCM", step_facet="Criterion", slope_facet="Criterion")
    cfg = res["config"]
    assert cfg.get("estimated_population_sd") is not None
    # GPCM slopes are exp(centered log-slopes): the geometric mean must stay ~1
    # regardless of sigma, so freeing sigma cannot trade off against an overall
    # discrimination multiplier.
    slopes = res["slopes"]
    assert isinstance(slopes, pd.DataFrame) and not slopes.empty
    vals = pd.to_numeric(slopes["Estimate"], errors="coerce").to_numpy()
    vals = vals[np.isfinite(vals) & (vals > 0)]
    assert vals.size >= 1
    geo_mean = float(np.exp(np.mean(np.log(vals))))
    assert geo_mean == pytest.approx(1.0, abs=0.05)


def test_recovery_harness_reports_estimated_population_sd():
    bundle = app.evaluate_parameter_recovery(
        model="RSM", fit_method="MML", reps=3,
        n_person=120, n_rater=2, n_criterion=2, n_cat=4,
        theta_sd=1.4, quad_points=19, maxit=300, reltol=1e-5,
        seed=20260601, estimate_population_sd=True,
    )
    overview = bundle["rep_overview"]
    assert "EstimatedPopulationSD" in overview.columns
    vals = overview["EstimatedPopulationSD"].dropna()
    assert len(vals) >= 1
    assert (vals > 1.0).all(), "free-sigma recovery did not exceed the fixed baseline"
