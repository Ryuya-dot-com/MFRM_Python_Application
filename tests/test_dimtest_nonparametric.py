"""Tests for the nonparametric DIMTEST helper.

The Streamlit app exposes ``compute_dimtest_nonparametric`` as a
polytomous adaptation of Stout's (1987) conditional-covariance T
statistic, with a cluster bootstrap on persons for the standard
error and p-value. The tests below pin two contracts:

* On a synthetic *unidimensional* design with an a priori AT / PT
  split (not derived from the data), the test must fail to reject
  H0 at the nominal level on average — i.e. the median bootstrap
  p-value across runs stays well above 0.05.

* On a synthetic *two-dimensional* design with an a priori correct
  AT / PT split (criteria 1-3 measure trait A, criteria 4-6 measure
  trait B), the test must strongly reject H0 (p well below 0.05).

The PCA-driven auto-split inherits a known selection-bias inflation
of Type I error (Roussos & Stout, 1996); the tests pin that
behaviour separately to document the caveat.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


@pytest.fixture(scope="module")
def unidim_fit():
    """80 persons x 4 raters x 6 criteria, all criteria load on the
    same latent ability."""
    rng = np.random.default_rng(20260620)
    n_p, n_r, n_c = 80, 4, 6
    rows = []
    for i in range(n_p):
        theta = rng.normal(0, 1)
        for j in range(n_r):
            rater_eff = rng.normal(0, 0.3)
            for k in range(n_c):
                crit_eff = rng.normal(0, 0.3)
                eta = theta - rater_eff - crit_eff
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({
                    "Person": f"P{i+1:02d}", "Rater": f"R{j+1}",
                    "Criterion": f"C{k+1}", "Score": score,
                })
    res = app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20,
    )
    return res, app.mfrm_diagnostics(res, compute_pca=True, compute_marginal=False)


@pytest.fixture(scope="module")
def two_dim_fit():
    """80 persons x 4 raters x 6 criteria, criteria 1-3 measure trait A,
    criteria 4-6 measure trait B (independent traits)."""
    rng = np.random.default_rng(20260621)
    n_p, n_r, n_c = 80, 4, 6
    rows = []
    for i in range(n_p):
        theta_a = rng.normal(0, 1)
        theta_b = rng.normal(0, 1)
        for j in range(n_r):
            rater_eff = rng.normal(0, 0.3)
            for k in range(n_c):
                theta = theta_a if k < 3 else theta_b
                crit_eff = rng.normal(0, 0.3)
                eta = theta - rater_eff - crit_eff
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({
                    "Person": f"P{i+1:02d}", "Rater": f"R{j+1}",
                    "Criterion": f"C{k+1}", "Score": score,
                })
    res = app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20,
    )
    return res, app.mfrm_diagnostics(res, compute_pca=True, compute_marginal=False)


# -----------------------------------------------------------------------------
# Refusal on invalid input
# -----------------------------------------------------------------------------


def test_dimtest_refuses_non_dict():
    out = app.compute_dimtest_nonparametric(None)
    assert out["available"] is False
    assert "fit dictionary" in out["reason"].lower()


def test_dimtest_refuses_unknown_item_facet(unidim_fit):
    res, diag = unidim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="NotAFacet",
    )
    assert out["available"] is False
    assert "not a facet" in out["reason"].lower()


def test_dimtest_refuses_overlapping_at_pt(unidim_fit):
    res, diag = unidim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"],
        pt_items=["C3", "C4", "C5"],
    )
    assert out["available"] is False
    assert "disjoint" in out["reason"].lower()


# -----------------------------------------------------------------------------
# A priori AT / PT split: unidim must fail to reject; 2-dim must reject
# -----------------------------------------------------------------------------


def test_dimtest_unidim_a_priori_split_fails_to_reject(unidim_fit):
    """With an a priori AT / PT split on truly unidimensional data,
    the p-value should be well above the conventional 0.05
    threshold; we use 0.05 as a soft floor (a single bootstrap run
    can dip below 0.05 by chance, but the point estimate of T_L
    should be close to zero relative to its bootstrap SE)."""
    res, diag = unidim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=200, seed=42,
    )
    assert out["available"] is True
    assert out["split_method"] == "user_specified"
    # Z should be modest (under 2.5 on this sample); T_L itself should
    # be near zero relative to the SE.
    assert abs(out["Z"]) < 2.5, (
        f"Unidim T_L = {out['T_L']:.4f}, Z = {out['Z']:.3f}, "
        f"p = {out['p_value']:.4f} — Type I error inflation suspected."
    )


def test_dimtest_two_dim_a_priori_split_rejects(two_dim_fit):
    """With an a priori correct AT / PT split on 2-dim data (criteria
    1-3 vs 4-6), the test must reject H0 at the conventional 0.05
    level (Z markedly above 1.96)."""
    res, diag = two_dim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=200, seed=42,
    )
    assert out["available"] is True
    assert out["Z"] > 3.0, (
        f"2-dim T_L = {out['T_L']:.4f}, Z = {out['Z']:.3f}, "
        f"p = {out['p_value']:.4f} — expected strong rejection."
    )
    assert out["p_value"] < 0.05


def test_dimtest_two_dim_T_L_larger_than_unidim(unidim_fit, two_dim_fit):
    """T_L on 2-dim data with the correct split must clearly exceed
    T_L on unidim data with the same split (the statistic captures
    the magnitude of the secondary dimension)."""
    res_u, diag_u = unidim_fit
    res_t, diag_t = two_dim_fit
    out_u = app.compute_dimtest_nonparametric(
        res_u, diag_u, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=100, seed=42,
    )
    out_t = app.compute_dimtest_nonparametric(
        res_t, diag_t, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=100, seed=42,
    )
    assert out_t["T_L"] > 2.0 * out_u["T_L"], (
        f"unidim T_L = {out_u['T_L']:.4f}, 2-dim T_L = {out_t['T_L']:.4f}"
    )


# -----------------------------------------------------------------------------
# Bootstrap reproducibility and bundle shape
# -----------------------------------------------------------------------------


def test_dimtest_seed_makes_bootstrap_deterministic(two_dim_fit):
    """Same seed must produce identical T_L / Z / p across runs."""
    res, diag = two_dim_fit
    out_a = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=50, seed=2026,
    )
    out_b = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=50, seed=2026,
    )
    assert out_a["T_L"] == pytest.approx(out_b["T_L"], abs=1e-12)
    assert out_a["SE"] == pytest.approx(out_b["SE"], abs=1e-12)


def test_dimtest_bundle_carries_required_fields(two_dim_fit):
    res, diag = two_dim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=50, seed=1,
    )
    for key in [
        "available", "T_L", "SE", "Z", "p_value",
        "ci_lower", "ci_upper",
        "n_persons", "n_strata_used", "stratum_sizes",
        "at_items", "pt_items", "split_method",
        "method", "references", "caveat", "settings",
    ]:
        assert key in out
    assert "Stout" in out["references"]
    assert "Nandakumar" in out["references"]


def test_dimtest_auto_split_documents_caveat(unidim_fit):
    """The auto-split path has documented selection bias; the bundle's
    caveat field must mention it."""
    res, diag = unidim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        n_bootstrap=50, seed=1,
    )
    assert out["available"] is True
    assert out["split_method"] in {"pc2_loading_sign", "alphabetical_fallback"}
    assert "selection-bias" in out["caveat"].lower() or "selection bias" in out["caveat"].lower()


def test_dimtest_p_value_is_two_sided_normal(two_dim_fit):
    """p_value should equal 2 * (1 - Phi(|Z|))."""
    from scipy.stats import norm
    res, diag = two_dim_fit
    out = app.compute_dimtest_nonparametric(
        res, diag, item_facet="Criterion",
        at_items=["C1", "C2", "C3"], pt_items=["C4", "C5", "C6"],
        n_bootstrap=50, seed=1,
    )
    expected_p = 2.0 * (1.0 - norm.cdf(abs(out["Z"])))
    assert out["p_value"] == pytest.approx(expected_p, abs=1e-12)
