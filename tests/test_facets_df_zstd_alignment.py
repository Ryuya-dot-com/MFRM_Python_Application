"""Tests for the FACETS / Wright-Masters d.f. and ZSTD reporting alignment.

Two d.f. conventions ship side-by-side:

* **Engine** (default, backwards compatible): ``DF_Infit = sum(Var * w)``,
  ``DF_Outfit = sum(w)``. The Wilson-Hilferty ZSTD assumes the per-cell
  variances are homogeneous enough that the residual chi-square behaves
  with that d.f.

* **FACETS** (Wright & Masters, 1982, Eqs. 4.20 / 4.27): a Welch-
  Satterthwaite d.f. that captures the variance heterogeneity by way of
  the fourth central moment of each per-observation polytomous
  distribution:

      DF_Infit_FACETS  = 2 * (sum sigma^2 w)^2 / sum w * (M4 - sigma^4)
      DF_Outfit_FACETS = 2 * (sum w)^2 / sum w * (M4 / sigma^4 - 1)

  where ``M4 = sum_k P_k * (k - E[X])^4``. The FACETS ZSTD function
  caps the output at ``+/- cap`` (default 9.0) to match FACETS
  software behaviour on extreme cells.

The math contract here pins the closed-form fourth moment, the
Wilson-Hilferty mapping, the cap, the engine/FACETS dispatch, and the
column shape under each dispatch mode.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Unit: FourthCentralMoment closed form on a tiny obs_df derived from probs
# -----------------------------------------------------------------------------


def _fourth_moment_from_probs(probs: np.ndarray, expected_k: np.ndarray) -> np.ndarray:
    k_vals = np.arange(probs.shape[1])
    diff = k_vals[np.newaxis, :] - expected_k[:, np.newaxis]
    return np.sum(probs * (diff ** 4), axis=1)


def test_fourth_central_moment_from_compute_obs_table_matches_closed_form():
    """For a small RSM JMLE fit the FourthCentralMoment column must equal the
    closed-form ``sum_k P_k * (k - E[X])^4`` evaluated from
    ``compute_prob_matrix``."""
    rng = np.random.default_rng(20260515)
    n_person, n_rater = 10, 2
    rows = []
    theta = rng.normal(0, 1, n_person)
    rater_eff = np.array([-0.3, 0.3])
    for i, theta_i in enumerate(theta):
        for j in range(n_rater):
            eta = theta_i - rater_eff[j]
            p1 = 1.0 / (1.0 + np.exp(-eta))
            p2 = 1.0 / (1.0 + np.exp(-(eta - 0.2)))
            score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
            rows.append({"Person": f"P{i+1:02d}", "Rater": f"R{j+1}", "Score": score})
    df = pd.DataFrame(rows)
    res = app.mfrm_estimate(
        data=df, person_col="Person", facet_cols=["Rater"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20, reltol=1e-4,
    )
    obs = app.compute_obs_table(res)
    probs = app.compute_prob_matrix(res)
    expected_k = pd.to_numeric(obs["Expected"], errors="coerce").to_numpy() - res["prep"]["rating_min"]
    expected_m4 = _fourth_moment_from_probs(np.asarray(probs, dtype=float), expected_k)
    actual_m4 = pd.to_numeric(obs["FourthCentralMoment"], errors="coerce").to_numpy()
    assert np.allclose(actual_m4, expected_m4, rtol=0, atol=1e-12)


def test_fourth_central_moment_is_nonneg_and_bounded_below_by_var_squared():
    """For any real random variable ``Var^2 <= M4`` (variance-fourth-moment
    Cauchy-Schwarz). The compute_obs_table column must satisfy that."""
    rng = np.random.default_rng(20260516)
    n_person = 12
    rows = []
    for i in range(n_person):
        theta = rng.normal(0, 1)
        for j in range(3):
            eta = theta - (j - 1) * 0.4
            p1 = 1.0 / (1.0 + np.exp(-eta))
            score = int(rng.uniform() < p1)
            rows.append({"Person": f"P{i+1:02d}", "Task": f"T{j+1}", "Score": score})
    res = app.mfrm_estimate(
        data=pd.DataFrame(rows), person_col="Person", facet_cols=["Task"], score_col="Score",
        rating_min=0, rating_max=1, model="RSM", method="JMLE", maxit=20, reltol=1e-4,
    )
    obs = app.compute_obs_table(res)
    var = pd.to_numeric(obs["Var"], errors="coerce").to_numpy()
    m4 = pd.to_numeric(obs["FourthCentralMoment"], errors="coerce").to_numpy()
    ok = np.isfinite(var) & np.isfinite(m4)
    assert (m4[ok] >= -1e-12).all()
    assert (m4[ok] + 1e-12 >= var[ok] ** 2).all()


# -----------------------------------------------------------------------------
# Unit: zstd_from_mnsq_facets — Wilson-Hilferty + cap
# -----------------------------------------------------------------------------


def test_zstd_from_mnsq_facets_matches_wilson_hilferty_when_below_cap():
    """For a moderate cell the FACETS ZSTD is the same Wilson-Hilferty value
    that the engine convention reports (the difference between the two
    only kicks in via the d.f. input and the cap)."""
    for mnsq, df in [(1.05, 50.0), (0.85, 20.0), (1.30, 100.0), (0.50, 8.0)]:
        z_engine = app.zstd_from_mnsq(mnsq, df, whexact=False)
        z_facets = app.zstd_from_mnsq_facets(mnsq, df, whexact=False, cap=9.0)
        assert z_facets == pytest.approx(z_engine, abs=1e-12)


def test_zstd_from_mnsq_facets_caps_extreme_values():
    """An extreme cell can blow up the WH-z; the cap must clamp it to
    ``+/- cap``. Uses the linear-normal-approximation branch
    (``whexact=True``) so the pre-cap value has a clean closed form
    ``(mnsq - 1) * sqrt(df / 2)``."""
    # whexact=True: z = (10 - 1) * sqrt(200 / 2) = 9 * 10 = 90 → cap to 9.
    z = app.zstd_from_mnsq_facets(10.0, 200.0, whexact=True, cap=9.0)
    assert z == pytest.approx(9.0, abs=1e-12)
    # z = (0.05 - 1) * sqrt(200 / 2) = -0.95 * 10 = -9.5 → cap to -9.
    z = app.zstd_from_mnsq_facets(0.05, 200.0, whexact=True, cap=9.0)
    assert z == pytest.approx(-9.0, abs=1e-12)
    # Wilson-Hilferty branch hits the cap on a tiny-df pathological cell.
    z = app.zstd_from_mnsq_facets(10.0, 1e-4, whexact=False, cap=9.0)
    assert z == pytest.approx(9.0, abs=1e-12)


def test_zstd_from_mnsq_facets_returns_nan_for_invalid_inputs():
    assert np.isnan(app.zstd_from_mnsq_facets(np.nan, 10.0))
    assert np.isnan(app.zstd_from_mnsq_facets(1.0, np.nan))
    assert np.isnan(app.zstd_from_mnsq_facets(1.0, 0.0))
    assert np.isnan(app.zstd_from_mnsq_facets(-1.0, 10.0))


def test_zstd_from_mnsq_facets_disables_cap_when_none_or_zero():
    """Passing cap=None or cap <= 0 must skip the clamp."""
    # whexact=True: |z| = (10 - 1) * sqrt(200/2) = 90, well outside cap.
    z = app.zstd_from_mnsq_facets(10.0, 200.0, whexact=True, cap=None)
    assert np.isfinite(z) and abs(z) == pytest.approx(90.0, abs=1e-12)
    z = app.zstd_from_mnsq_facets(10.0, 200.0, whexact=True, cap=0.0)
    assert np.isfinite(z) and abs(z) == pytest.approx(90.0, abs=1e-12)


# -----------------------------------------------------------------------------
# Unit: FACETS d.f. Welch-Satterthwaite formulas
# -----------------------------------------------------------------------------


def test_facets_fit_df_terms_match_wright_masters_formula():
    """Hand-rolled fixture: 3 observations with known var/M4/weight; the
    denominator sums must match the Wright-Masters definitions."""
    var = np.array([0.20, 0.50, 0.10], dtype=float)
    fourth = np.array([0.08, 0.32, 0.05], dtype=float)
    weight = np.array([1.0, 2.0, 1.0], dtype=float)
    expected_infit = float(np.sum(weight * (fourth - var ** 2)))
    expected_outfit = float(np.sum(weight * (fourth / (var ** 2) - 1.0)))
    actual_infit, actual_outfit = app._facets_fit_df_terms(var, fourth, weight)
    assert actual_infit == pytest.approx(expected_infit, abs=1e-12)
    assert actual_outfit == pytest.approx(expected_outfit, abs=1e-12)


def test_facets_fit_df_matches_welch_satterthwaite_definition():
    """DF = 2 * (cell info)^2 / denominator, by the Welch-Satterthwaite
    derivation. Pin both Infit and Outfit branches."""
    sum_var_w = 5.0
    sum_w = 30.0
    denom_infit = 1.5
    denom_outfit = 4.0
    df_infit, df_outfit = app._facets_fit_df(sum_var_w, sum_w, denom_infit, denom_outfit)
    assert df_infit == pytest.approx(2.0 * 5.0 ** 2 / 1.5, abs=1e-12)
    assert df_outfit == pytest.approx(2.0 * 30.0 ** 2 / 4.0, abs=1e-12)


def test_facets_fit_df_returns_nan_when_denominators_collapse():
    """A near-zero denominator must surface NaN rather than an inf."""
    df_infit, df_outfit = app._facets_fit_df(5.0, 30.0, 0.0, 0.0)
    assert np.isnan(df_infit)
    assert np.isnan(df_outfit)


# -----------------------------------------------------------------------------
# Dispatch: engine / facets / both
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def small_rsm_jmle_diagnostics():
    """Reused across the dispatch tests so the fit cost is paid once."""
    rng = np.random.default_rng(20260517)
    n_person, n_rater, n_criterion = 25, 2, 3
    rows = []
    for i in range(n_person):
        theta = rng.normal(0, 1)
        for j in range(n_rater):
            for k in range(n_criterion):
                eta = theta - (j - 0.5) * 0.4 - (k - 1) * 0.4
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append({"Person": f"P{i+1:02d}", "Rater": f"R{j+1}",
                             "Criterion": f"C{k+1}", "Score": score})
    df = pd.DataFrame(rows)
    res = app.mfrm_estimate(
        data=df, person_col="Person", facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE", maxit=20, reltol=1e-4,
    )
    return res


def test_calc_facet_fit_engine_mode_is_backward_compatible(small_rsm_jmle_diagnostics):
    """In default ``"engine"`` mode the table must carry the original
    columns and *not* expose any FACETS suffix columns; this guarantees
    backward compatibility with downstream consumers built before the
    alignment shipped."""
    obs = app.compute_obs_table(small_rsm_jmle_diagnostics)
    tbl = app.calc_facet_fit(obs, ["Person", "Rater", "Criterion"])
    assert {"Infit", "Outfit", "InfitZSTD", "OutfitZSTD",
            "DF_Infit", "DF_Outfit"}.issubset(tbl.columns)
    for col in ["DF_Infit_FACETS", "DF_Outfit_FACETS",
                "InfitZSTD_FACETS", "OutfitZSTD_FACETS",
                "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE",
                "FitDfMethod", "FitZSTDTransform", "FitZSTDCap"]:
        assert col not in tbl.columns


def test_calc_facet_fit_facets_mode_swaps_primary_zstd(small_rsm_jmle_diagnostics):
    """In ``"facets"`` mode the primary ZSTD columns carry the FACETS values;
    the original engine values must still be available under
    ``*_ENGINE`` suffix columns so callers can compare conventions."""
    obs = app.compute_obs_table(small_rsm_jmle_diagnostics)
    engine_tbl = app.calc_facet_fit(obs, ["Rater", "Criterion"])
    facets_tbl = app.calc_facet_fit(obs, ["Rater", "Criterion"], fit_df_method="facets")
    # Engine values preserved as *_ENGINE.
    assert "InfitZSTD_ENGINE" in facets_tbl.columns
    assert "OutfitZSTD_ENGINE" in facets_tbl.columns
    # Engine *_ENGINE column matches the engine-mode primary value.
    pd.testing.assert_series_equal(
        facets_tbl["InfitZSTD_ENGINE"].reset_index(drop=True),
        engine_tbl["InfitZSTD"].reset_index(drop=True),
        check_names=False,
    )
    # Primary InfitZSTD now reflects the FACETS d.f.
    assert "InfitZSTD_FACETS" in facets_tbl.columns
    pd.testing.assert_series_equal(
        facets_tbl["InfitZSTD"].reset_index(drop=True),
        facets_tbl["InfitZSTD_FACETS"].reset_index(drop=True),
        check_names=False,
    )
    # Metadata columns set.
    assert (facets_tbl["FitDfMethod"] == "facets_wright_masters").all()
    assert (facets_tbl["FitZSTDTransform"] == "Wilson-Hilferty").all()
    assert (facets_tbl["FitZSTDCap"] == 9.0).all()


def test_calc_facet_fit_both_mode_keeps_engine_primary_with_facets_suffix(
    small_rsm_jmle_diagnostics,
):
    """In ``"both"`` mode the primary ZSTD stays on the engine convention
    (preserving backward-compatible interpretations) but the FACETS
    suffix columns are available for side-by-side comparison."""
    obs = app.compute_obs_table(small_rsm_jmle_diagnostics)
    engine_tbl = app.calc_facet_fit(obs, ["Rater", "Criterion"])
    both_tbl = app.calc_facet_fit(obs, ["Rater", "Criterion"], fit_df_method="both")
    pd.testing.assert_series_equal(
        both_tbl["InfitZSTD"].reset_index(drop=True),
        engine_tbl["InfitZSTD"].reset_index(drop=True),
        check_names=False,
    )
    for col in ["InfitZSTD_FACETS", "OutfitZSTD_FACETS",
                "DF_Infit_FACETS", "DF_Outfit_FACETS",
                "InfitZSTD_ENGINE", "OutfitZSTD_ENGINE"]:
        assert col in both_tbl.columns
    assert (both_tbl["FitDfMethod"] == "engine_primary_facets_available").all()


def test_calc_overall_fit_dispatch_matches_per_facet_dispatch(small_rsm_jmle_diagnostics):
    """The overall-fit table must respect the same dispatch logic as the
    per-facet table — the two should never disagree on convention."""
    obs = app.compute_obs_table(small_rsm_jmle_diagnostics)
    overall_engine = app.calc_overall_fit(obs)
    overall_facets = app.calc_overall_fit(obs, fit_df_method="facets")
    overall_both = app.calc_overall_fit(obs, fit_df_method="both")
    assert "FitDfMethod" not in overall_engine.columns
    assert overall_facets.loc[0, "FitDfMethod"] == "facets_wright_masters"
    assert overall_both.loc[0, "FitDfMethod"] == "engine_primary_facets_available"


def test_mfrm_diagnostics_propagates_fit_df_method(small_rsm_jmle_diagnostics):
    """The bundle returned by ``mfrm_diagnostics`` must surface the
    requested method and pin its row-level value on every row of the
    fit table; downstream UI code reads ``diagnostics['fit_df_method']``
    as the canonical method label."""
    diag_engine = app.mfrm_diagnostics(small_rsm_jmle_diagnostics, compute_pca=False)
    diag_facets = app.mfrm_diagnostics(
        small_rsm_jmle_diagnostics, compute_pca=False, fit_df_method="facets"
    )
    diag_both = app.mfrm_diagnostics(
        small_rsm_jmle_diagnostics, compute_pca=False, fit_df_method="both"
    )
    assert diag_engine["fit_df_method"] == "engine"
    assert diag_facets["fit_df_method"] == "facets"
    assert diag_both["fit_df_method"] == "both"
    # The fit table's status column matches the bundle-level method.
    assert "FitDfMethod" not in diag_engine["fit"].columns
    assert (diag_facets["fit"]["FitDfMethod"] == "facets_wright_masters").all()
    assert (diag_both["fit"]["FitDfMethod"] == "engine_primary_facets_available").all()


# -----------------------------------------------------------------------------
# Wilson-Hilferty transform label
# -----------------------------------------------------------------------------


def test_fit_zstd_transform_label_distinguishes_wh_from_linear_normal():
    assert app.fit_zstd_transform_label(whexact=False) == "Wilson-Hilferty"
    assert app.fit_zstd_transform_label(whexact=True) == "linear normal approximation"
