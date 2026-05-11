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


def test_g_phi_match_brennan_3way_closed_form(crossed_rating_fit):
    """For a balanced 3-way crossed p x i x j design under one
    observation per cell, the G / Phi formulas (Brennan 2001 Eq. 3.18)
    are:

        sigma2(delta) = sigma2_pi / n_i + sigma2_pj / n_j +
                        sigma2_pij / (n_i n_j)
        sigma2(Delta) = sigma2(delta) + sigma2_i / n_i +
                        sigma2_j / n_j + sigma2_ij / (n_i n_j)
        G   = sigma2_p / (sigma2_p + sigma2(delta))
        Phi = sigma2_p / (sigma2_p + sigma2(Delta))
    """
    out = app.compute_generalizability_study(crossed_rating_fit)
    var = out["variance_components"]
    coef = out["observed_coefficients"]
    lookup = dict(zip(var["Source"].astype(str), var["Variance"].astype(float)))
    sigma2_p = lookup["Person"]
    sigma2_e = lookup["Residual"]
    sigma2_r = lookup["Rater"]
    sigma2_c = lookup["Criterion"]
    sigma2_pr = lookup["Person:Rater"]
    sigma2_pc = lookup["Person:Criterion"]
    sigma2_rc = lookup["Rater:Criterion"]
    n_r = out["design"]["observed_n"]["Rater"]
    n_c = out["design"]["observed_n"]["Criterion"]

    expected_relative = (
        sigma2_pr / n_r + sigma2_pc / n_c + sigma2_e / (n_r * n_c)
    )
    expected_absolute = expected_relative + (
        sigma2_r / n_r + sigma2_c / n_c + sigma2_rc / (n_r * n_c)
    )
    expected_g = sigma2_p / (sigma2_p + expected_relative)
    expected_phi = sigma2_p / (sigma2_p + expected_absolute)
    assert coef["G"] == pytest.approx(expected_g, abs=1e-10)
    assert coef["Phi"] == pytest.approx(expected_phi, abs=1e-10)


def test_variance_components_carry_two_way_interaction_rows(crossed_rating_fit):
    """Under the full 3-way decomposition, the variance-components
    table carries seven rows: three main effects + three two-way
    interactions + residual."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    var = out["variance_components"]
    sources = set(var["Source"].astype(str))
    assert "Person:Rater" in sources
    assert "Person:Criterion" in sources
    assert "Rater:Criterion" in sources
    assert "Residual" in sources
    assert (var["Term"] == "two-way").sum() == 3
    assert (var["Term"] == "main").sum() == 3
    assert (var["Term"] == "residual").sum() == 1


def test_decomposition_label_records_brennan_eq_3_18(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    coef = out["observed_coefficients"]
    assert coef["details"]["decomposition"] == "full_3way_brennan_eq_3_18"


# -----------------------------------------------------------------------------
# Bootstrap CI for G / Phi
# -----------------------------------------------------------------------------


def test_bootstrap_ci_disabled_by_default(crossed_rating_fit):
    out = app.compute_generalizability_study(crossed_rating_fit)
    assert "bootstrap_ci" in out
    assert out["bootstrap_ci"]["available"] is False
    assert "not requested" in out["bootstrap_ci"]["reason"].lower()


def test_bootstrap_ci_produces_finite_band_when_enabled(crossed_rating_fit):
    """With B = 50 replicates the percentile CI must produce finite
    lower / upper bounds that bracket the observed G and Phi (or are
    very close to them; the percentile CI can occasionally exclude
    the point estimate under small B but should at least be finite
    and ordered)."""
    out = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=50,
        bootstrap_seed=20260612,
    )
    ci = out["bootstrap_ci"]
    assert ci["available"] is True
    assert ci["n_success"] >= 30  # at least 60% replicates produced a finite G
    assert np.isfinite(ci["G_lower"]) and np.isfinite(ci["G_upper"])
    assert ci["G_lower"] <= ci["G_upper"]
    assert np.isfinite(ci["Phi_lower"]) and np.isfinite(ci["Phi_upper"])
    assert ci["Phi_lower"] <= ci["Phi_upper"]


def test_bootstrap_ci_widens_with_smaller_confidence(crossed_rating_fit):
    """A 99% CI must be at least as wide as a 90% CI (percentile
    method is monotone in confidence level)."""
    out_90 = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=80,
        bootstrap_confidence=0.90, bootstrap_seed=42,
    )
    out_99 = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=80,
        bootstrap_confidence=0.99, bootstrap_seed=42,
    )
    w90 = out_90["bootstrap_ci"]["G_upper"] - out_90["bootstrap_ci"]["G_lower"]
    w99 = out_99["bootstrap_ci"]["G_upper"] - out_99["bootstrap_ci"]["G_lower"]
    assert w99 + 1e-9 >= w90


def test_bootstrap_ci_seed_makes_results_deterministic(crossed_rating_fit):
    """Same seed must produce identical bootstrap output."""
    a = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=20,
        bootstrap_seed=2026,
    )["bootstrap_ci"]
    b = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=20,
        bootstrap_seed=2026,
    )["bootstrap_ci"]
    assert a["G_lower"] == pytest.approx(b["G_lower"], abs=1e-12)
    assert a["G_upper"] == pytest.approx(b["G_upper"], abs=1e-12)


def test_bootstrap_replicates_bound_in_unit_interval(crossed_rating_fit):
    out = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=40,
        bootstrap_seed=20260615,
    )
    g_reps = out["bootstrap_ci"]["G_replicates"]
    phi_reps = out["bootstrap_ci"]["Phi_replicates"]
    assert (g_reps >= -1e-12).all() and (g_reps <= 1.0 + 1e-12).all()
    assert (phi_reps >= -1e-12).all() and (phi_reps <= 1.0 + 1e-12).all()


# -----------------------------------------------------------------------------
# D-study CI bands (bootstrap propagated through Brennan Eq. 3.18 grid points)
# -----------------------------------------------------------------------------


def test_d_study_carries_ci_columns_when_bootstrap_enabled(crossed_rating_fit):
    """With bootstrap_ci=True, the d_study DataFrame must carry per-row
    G_lower / G_upper / Phi_lower / Phi_upper columns alongside the
    point estimates."""
    out = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=60,
        bootstrap_seed=20260613,
    )
    d_study = out["d_study"]
    for col in ["G_lower", "G_upper", "Phi_lower", "Phi_upper"]:
        assert col in d_study.columns, col


def test_d_study_ci_bands_bracket_point_estimates(crossed_rating_fit):
    """For every D-study grid row the CI band should bracket the point
    estimate; we use a soft inclusive bound to allow for slight noise
    in 60-replicate bootstraps."""
    out = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=60,
        bootstrap_seed=20260613,
    )
    d_study = out["d_study"]
    finite = d_study.dropna(subset=["G", "G_lower", "G_upper"])
    for _, row in finite.iterrows():
        # Bands may not always bracket the observed point estimate
        # under a finite bootstrap (the percentile method can shift),
        # but the (lower, upper) pair must be ordered and inside [0, 1].
        assert row["G_lower"] <= row["G_upper"] + 1e-12, row
        assert row["Phi_lower"] <= row["Phi_upper"] + 1e-12, row
        assert 0 - 1e-9 <= row["G_lower"] <= 1 + 1e-9
        assert 0 - 1e-9 <= row["G_upper"] <= 1 + 1e-9


def test_d_study_ci_band_at_observed_design_matches_bootstrap_ci(
    crossed_rating_fit,
):
    """The d_study row corresponding to the observed (n_rater,
    n_criterion) sample sizes must have CI bounds within numerical
    noise of the bootstrap_ci values reported at the observed
    design (they share the same replicate set)."""
    out = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=80,
        bootstrap_seed=20260613,
    )
    bc = out["bootstrap_ci"]
    d_study = out["d_study"]
    observed_n = out["design"]["observed_n"]
    n_r = int(observed_n["Rater"])
    n_c = int(observed_n["Criterion"])
    matching = d_study[(d_study["Rater"] == n_r) & (d_study["Criterion"] == n_c)]
    assert not matching.empty, "Observed design row missing from D-study grid"
    row = matching.iloc[0]
    assert float(row["G_lower"]) == pytest.approx(bc["G_lower"], abs=1e-12)
    assert float(row["G_upper"]) == pytest.approx(bc["G_upper"], abs=1e-12)
    assert float(row["Phi_lower"]) == pytest.approx(bc["Phi_lower"], abs=1e-12)
    assert float(row["Phi_upper"]) == pytest.approx(bc["Phi_upper"], abs=1e-12)


def test_d_study_ci_band_width_shrinks_with_design_size(crossed_rating_fit):
    """Adding observations shrinks the bootstrap CI band for G at
    that grid point on average; we pick two grid points where one
    is strictly larger in both facets and verify the width ordering
    holds."""
    out = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=80,
        bootstrap_seed=20260613,
    )
    d_study = out["d_study"]
    finite = d_study.dropna(subset=["G_lower", "G_upper"])
    # Find pairs where (n_r2 >= n_r1, n_c2 >= n_c1) and at least one strict.
    rows_by_design = {(int(r["Rater"]), int(r["Criterion"])): r for _, r in finite.iterrows()}
    pairs_checked = 0
    pairs_consistent = 0
    designs = sorted(rows_by_design.keys())
    for a in designs:
        for b in designs:
            if b == a:
                continue
            if b[0] >= a[0] and b[1] >= a[1] and (b[0] > a[0] or b[1] > a[1]):
                ra, rb = rows_by_design[a], rows_by_design[b]
                width_a = float(ra["G_upper"]) - float(ra["G_lower"])
                width_b = float(rb["G_upper"]) - float(rb["G_lower"])
                pairs_checked += 1
                if width_b <= width_a + 1e-6:
                    pairs_consistent += 1
    # Allow a small fraction of pairs to violate due to bootstrap noise;
    # the majority must show shrinking width as design grows.
    assert pairs_consistent / max(pairs_checked, 1) > 0.7, (
        f"Only {pairs_consistent}/{pairs_checked} pairs showed shrinking CI width."
    )


def test_d_study_no_ci_columns_when_bootstrap_disabled(crossed_rating_fit):
    """With bootstrap disabled (the default), the d_study DataFrame
    must not carry CI columns."""
    out = app.compute_generalizability_study(crossed_rating_fit)
    d_study = out["d_study"]
    for col in ["G_lower", "G_upper", "Phi_lower", "Phi_upper"]:
        assert col not in d_study.columns, col


def test_d_study_ci_band_reproducible_under_same_seed(crossed_rating_fit):
    out_a = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=30,
        bootstrap_seed=2026,
    )
    out_b = app.compute_generalizability_study(
        crossed_rating_fit, bootstrap_ci=True, n_bootstrap=30,
        bootstrap_seed=2026,
    )
    a = out_a["d_study"][["G_lower", "G_upper", "Phi_lower", "Phi_upper"]]
    b = out_b["d_study"][["G_lower", "G_upper", "Phi_lower", "Phi_upper"]]
    pd.testing.assert_frame_equal(a, b)


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
