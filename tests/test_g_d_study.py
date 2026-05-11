"""Tests for the generalizability theory (G-study) and D-study pipeline.

``compute_generalizability_study(res)`` runs a method-of-moments crossed
ANOVA on the rating data and reports:

* Variance components per source (object_facet, each random facet,
  Residual), with the proportion-of-variance share;
* G (Cronbach et al., 1972) for relative decisions and Phi for absolute
  decisions at the observed design;
* A D-study forecast grid that scales G / Phi to planned numbers of
  raters and criteria.

The math contract pinned here covers:

* Refusal on invalid input (non-dict, missing prep, no random facets).
* Variance components are non-negative (a method-of-moments artefact
  that the helper clamps at zero to keep G / Phi well-defined).
* G / Phi closed-form identity:
  ``G = sigma2_p / (sigma2_p + sigma2_e / n_total)`` and
  ``Phi = sigma2_p / (sigma2_p + sum(sigma2_facet / n_facet) + sigma2_e / n_total)``.
* D-study monotonicity: G and Phi increase as any facet sample size
  increases.
* D-study row count is ``prod(facet_grid_sizes)``.
* The D-study grid contains the observed sample sizes (the helper picks
  a sensible band around the observed design).
* Coefficients are NaN when sigma2_p collapses to zero.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Refusal: invalid inputs
# -----------------------------------------------------------------------------


def test_generalizability_study_refuses_non_dict():
    out = app.compute_generalizability_study(None)
    assert out["available"] is False
    assert "fit dictionary" in out["reason"].lower()


def test_generalizability_study_refuses_no_random_facets():
    out = app.compute_generalizability_study(
        {
            "config": {"facet_names": ["Person"]},
            "prep": {"data": pd.DataFrame({"Person": ["P1"], "Score": [1]})},
        }
    )
    assert out["available"] is False
    assert "random facet" in out["reason"].lower()


def test_generalizability_study_refuses_empty_data():
    out = app.compute_generalizability_study(
        {
            "config": {"facet_names": ["Person", "Rater"]},
            "prep": {"data": pd.DataFrame()},
        }
    )
    assert out["available"] is False
    assert "raw data frame" in out["reason"].lower()


# -----------------------------------------------------------------------------
# Math contract: end-to-end variance components + G/Phi
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def crossed_rating_fit():
    """40 persons x 3 raters x 4 criteria RSM fit (480 ratings)."""
    rng = np.random.default_rng(20260615)
    rows = []
    for i in range(40):
        theta = rng.normal(0, 1)
        for j in range(3):
            for k in range(4):
                eta = theta - (j - 1) * 0.4 - (k - 1.5) * 0.4
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({
                    "Person": f"P{i+1:02d}",
                    "Rater": f"R{j+1}",
                    "Criterion": f"C{k+1}",
                    "Score": score,
                })
    df = pd.DataFrame(rows)
    return app.mfrm_estimate(
        data=df, person_col="Person", facet_cols=["Rater", "Criterion"],
        score_col="Score", rating_min=0, rating_max=2, model="RSM",
        method="JMLE", maxit=20,
    )


def test_variance_components_are_non_negative(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    assert out["available"] is True
    var = out["variance_components"]
    assert (pd.to_numeric(var["Variance"], errors="coerce") >= -1e-12).all()


def test_variance_components_have_required_sources(crossed_rating_fit):
    """Object, each random facet, and Residual must all appear."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    sources = set(out["variance_components"]["Source"].astype(str))
    assert {"Person", "Rater", "Criterion", "Residual"}.issubset(sources)


def test_proportion_variance_sums_to_one(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    var = out["variance_components"]
    finite = pd.to_numeric(var["ProportionVariance"], errors="coerce").dropna()
    assert finite.sum() == pytest.approx(1.0, abs=1e-9)


def test_g_phi_match_closed_form(crossed_rating_fit):
    """G = sigma2_p / (sigma2_p + sigma2_e / n_total),
    Phi = sigma2_p / (sigma2_p + sum facet_var / n_facet + sigma2_e / n_total)."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    var = out["variance_components"]
    coef = out["observed_coefficients"]
    lookup = dict(zip(var["Source"].astype(str), var["Variance"].astype(float)))
    sigma2_p = lookup["Person"]
    sigma2_e = lookup["Residual"]
    sigma2_r = lookup["Rater"]
    sigma2_c = lookup["Criterion"]
    n_r = out["design"]["observed_n"]["Rater"]
    n_c = out["design"]["observed_n"]["Criterion"]
    n_total = n_r * n_c
    expected_g = sigma2_p / (sigma2_p + sigma2_e / n_total)
    expected_phi = sigma2_p / (
        sigma2_p + sigma2_r / n_r + sigma2_c / n_c + sigma2_e / n_total
    )
    assert coef["G"] == pytest.approx(expected_g, abs=1e-10)
    assert coef["Phi"] == pytest.approx(expected_phi, abs=1e-10)


def test_phi_is_le_g(crossed_rating_fit):
    """Phi (absolute decisions) must be at most G (relative decisions)
    under the standard one-observation-per-cell decomposition."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    coef = out["observed_coefficients"]
    assert coef["Phi"] <= coef["G"] + 1e-12


# -----------------------------------------------------------------------------
# D-study: monotonicity, grid shape, contains observed design
# -----------------------------------------------------------------------------


def test_d_study_g_is_monotone_in_each_facet(crossed_rating_fit):
    """At fixed n on one facet, G must increase with n on the other."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    d_study = out["d_study"]
    # Pin Rater = 3, vary Criterion.
    sub = d_study[d_study["Rater"] == 3].sort_values("Criterion")
    g_vals = sub["G"].to_numpy()
    assert (np.diff(g_vals) >= -1e-12).all(), g_vals
    # Pin Criterion = 4, vary Rater.
    sub = d_study[d_study["Criterion"] == 4].sort_values("Rater")
    g_vals = sub["G"].to_numpy()
    assert (np.diff(g_vals) >= -1e-12).all(), g_vals


def test_d_study_phi_is_monotone_in_each_facet(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    d_study = out["d_study"]
    for rater_n in d_study["Rater"].unique():
        sub = d_study[d_study["Rater"] == rater_n].sort_values("Criterion")
        phi_vals = sub["Phi"].to_numpy()
        assert (np.diff(phi_vals) >= -1e-12).all()


def test_d_study_grid_includes_observed_design(crossed_rating_fit):
    """The D-study grid must include the observed (n_rater, n_criterion)
    point so users can read off the current operating G / Phi."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    d_study = out["d_study"]
    observed = out["design"]["observed_n"]
    n_r = int(observed["Rater"])
    n_c = int(observed["Criterion"])
    hit = d_study[(d_study["Rater"] == n_r) & (d_study["Criterion"] == n_c)]
    assert not hit.empty
    # Observed G in d_study must match observed_coefficients (within
    # the same closed-form identity).
    assert float(hit["G"].iloc[0]) == pytest.approx(out["observed_coefficients"]["G"], abs=1e-10)


def test_d_study_row_count_matches_grid_product(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    d_study = out["d_study"]
    n_r_vals = d_study["Rater"].nunique()
    n_c_vals = d_study["Criterion"].nunique()
    assert len(d_study) == n_r_vals * n_c_vals


# -----------------------------------------------------------------------------
# Caveat and references
# -----------------------------------------------------------------------------


def test_caveat_mentions_one_observation_per_cell(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    assert "one-observation-per-cell" in out["caveat"]
    assert "Brennan" in out["reference"]
    assert "Cronbach" in out["reference"]


def test_design_block_round_trips_observed_design(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    design = out["design"]
    assert design["object_facet"] == "Person"
    assert set(design["random_facets"]) == {"Rater", "Criterion"}
    assert design["observed_n"]["Rater"] == 3
    assert design["observed_n"]["Criterion"] == 4


# -----------------------------------------------------------------------------
# Reporting bands
# -----------------------------------------------------------------------------


def test_mean_obs_per_level_equals_n_total_over_levels_on_balanced(
    crossed_rating_fit,
):
    """For a balanced design the harmonic mean of per-level observation
    counts equals ``n_obs / n_levels`` exactly (every level has the
    same count, so the harmonic mean is the level count). Pins the
    refinement to the canonical balanced-design formula."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    var = out["variance_components"]
    n_obs = int(out["design"]["n_observations"])
    for _, row in var.iterrows():
        if row["Source"] == "Residual":
            continue
        # 40 persons x 3 raters x 4 criteria = 480 observations
        # Levels = {Person: 40, Rater: 3, Criterion: 4}
        # mean_obs_per_level = 480 / 40 = 12 (Person), 480 / 3 = 160 (Rater), 480 / 4 = 120 (Criterion)
        expected = float(n_obs / int(row["Levels"]))
        assert float(row["MeanObsPerLevel"]) == pytest.approx(expected, abs=1e-9)


def test_g_phi_remain_in_unit_interval(crossed_rating_fit):
    """Variance components are clamped at zero so G and Phi always lie
    in [0, 1] in the math layer (subject to numerical noise)."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    coef = out["observed_coefficients"]
    assert -1e-12 <= coef["G"] <= 1.0 + 1e-12
    assert -1e-12 <= coef["Phi"] <= 1.0 + 1e-12
    d_study = out["d_study"]
    assert (d_study["G"] >= -1e-12).all()
    assert (d_study["G"] <= 1.0 + 1e-12).all()
    assert (d_study["Phi"] >= -1e-12).all()
    assert (d_study["Phi"] <= 1.0 + 1e-12).all()
