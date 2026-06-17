"""Tests for the polytomous person-fit indices: lz and Snijders lz*.

The Streamlit app ports the mfrmr 0.2.0 person-fit pipeline:

* Drasgow, Levine, and Williams (1985) define the standardised
  polytomous log-likelihood

      lz = (l - E[l]) / sqrt(Var[l])

  where l, E[l], Var[l] are aggregated per person from the per-item
  log-probability of the observed category.

* Snijders (2001, Eq. 16) corrects lz for the fact that theta_n was
  estimated. For JML estimates the correction substitutes the
  conditional information at the JML estimate:

      c_n           = Cov[l, S] / I(theta)
      corrected_var = Var[l]    - Cov[l, S]^2 / I(theta)
      lz*           = (l - E[l] - c_n * S) / sqrt(corrected_var)

  Under JML the score S sums to zero by construction, but its
  variance contribution is what restores the N(0, 1) null.

The tests below pin the math contract directly (hand-built obs
frames with known closed-form lz/lz*), and then verify that the
end-to-end pipeline emits the correct status fields on both JMLE
and MML fits.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

import streamlit_app as app


# -----------------------------------------------------------------------------
# Unit: person_fit_threshold_table
# -----------------------------------------------------------------------------


def test_person_fit_threshold_table_returns_5pct_and_1pct():
    tbl = app.person_fit_threshold_table()
    assert list(tbl["Threshold"]) == ["5pct", "1pct"]
    assert tbl.loc[0, "TwoSidedAlpha"] == pytest.approx(0.05)
    assert tbl.loc[1, "TwoSidedAlpha"] == pytest.approx(0.01)
    assert tbl.loc[0, "AbsZ"] == pytest.approx(float(norm.ppf(0.975)))
    assert tbl.loc[1, "AbsZ"] == pytest.approx(float(norm.ppf(0.995)))
    assert tbl.loc[0, "Rule"] == "|z| > 1.96"
    assert tbl.loc[1, "Rule"] == "|z| > 2.58"


# -----------------------------------------------------------------------------
# Math contract: closed-form lz on a hand-built obs frame
# -----------------------------------------------------------------------------


def _build_synthetic_obs(
    *,
    person_id: str = "P1",
    pr_observed=(0.6, 0.4, 0.7),
    item_entropy=(-0.9, -0.7, -0.8),
    item_var_logp=(0.20, 0.15, 0.25),
    extra_cols: dict | None = None,
) -> pd.DataFrame:
    """Build a per-observation diagnostics frame with hand-set inputs.

    Returns a DataFrame with the minimum columns required by
    ``compute_person_fit_indices`` (``Person``, ``PrObserved``,
    ``ItemEntropy``, ``ItemVarLogP``), optionally augmented with
    additional columns (e.g. the lz* covariance / information terms).
    """
    n = len(pr_observed)
    base = {
        "Person": [person_id] * n,
        "PrObserved": list(pr_observed),
        "ItemEntropy": list(item_entropy),
        "ItemVarLogP": list(item_var_logp),
    }
    if extra_cols:
        for k, v in extra_cols.items():
            base[k] = list(v)
    return pd.DataFrame(base)


def test_lz_matches_closed_form_on_hand_built_obs():
    """lz = (sum log P - sum E[log P]) / sqrt(sum Var[log P]) by definition."""
    pr_observed = (0.6, 0.4, 0.7)
    item_entropy = (-0.9, -0.7, -0.8)
    item_var_logp = (0.20, 0.15, 0.25)
    obs = _build_synthetic_obs(
        pr_observed=pr_observed,
        item_entropy=item_entropy,
        item_var_logp=item_var_logp,
    )
    # Compute the contract by hand.
    log_p = np.log(np.array(pr_observed))
    ll = float(np.sum(log_p))
    e_ll = float(np.sum(item_entropy))
    var_ll = float(np.sum(item_var_logp))
    expected_lz = (ll - e_ll) / np.sqrt(var_ll)

    # Run through the public pipeline; res = None means lz is still
    # computed but lz* is left as ``fit_required``.
    out = app.compute_person_fit_indices({"config": {}, "prep": {"levels": {"Person": ["P1"]}}}, obs=obs)
    row = out.iloc[0]
    assert row["N"] == 3
    assert row["LogLik"] == pytest.approx(ll)
    assert row["lz"] == pytest.approx(expected_lz, abs=1e-12)
    # No fit => lz* must be NaN with a non-success status.
    assert pd.isna(row["lz_star"])
    assert row["lz_star_status"] != "computed_jml_conditional_calibration"


# -----------------------------------------------------------------------------
# Math contract: closed-form lz* on a hand-built obs frame
# -----------------------------------------------------------------------------


def test_lz_star_matches_closed_form_on_hand_built_obs():
    """lz* = (l - E[l] - c_n * S) / sqrt(Var[l] - Cov^2 / I) by Snijders (2001)."""
    pr_observed = (0.55, 0.62, 0.48)
    item_entropy = (-0.95, -0.80, -1.02)
    item_var_logp = (0.18, 0.22, 0.27)
    item_cov = (0.05, -0.04, 0.08)
    score_info = (0.40, 0.55, 0.45)
    obs_score_deriv = (0.10, -0.08, 0.12)

    obs = _build_synthetic_obs(
        pr_observed=pr_observed,
        item_entropy=item_entropy,
        item_var_logp=item_var_logp,
        extra_cols={
            "ItemLogPScoreCov": item_cov,
            "ScoreInformation": score_info,
            "ObservedScoreDerivative": obs_score_deriv,
        },
    )

    log_p = np.log(np.array(pr_observed))
    ll = float(np.sum(log_p))
    e_ll = float(np.sum(item_entropy))
    var_ll = float(np.sum(item_var_logp))
    info_total = float(np.sum(score_info))
    cov_total = float(np.sum(item_cov))
    score_sum = float(np.sum(obs_score_deriv))
    c_n = cov_total / info_total
    corrected_var = var_ll - (cov_total ** 2) / info_total
    expected_lz_star = (ll - e_ll - c_n * score_sum) / np.sqrt(corrected_var)

    res = {
        "config": {"method": "JMLE"},
        "prep": {"levels": {"Person": ["P1"]}},
    }
    out = app.compute_person_fit_indices(res, obs=obs)
    row = out.iloc[0]
    assert row["lz_star_status"] == "computed_jml_conditional_calibration"
    assert row["lz_star"] == pytest.approx(expected_lz_star, abs=1e-12)
    assert row["lz_star_c"] == pytest.approx(c_n, abs=1e-12)
    assert row["lz_star_variance"] == pytest.approx(corrected_var, abs=1e-12)


def test_lz_star_eap_population_corrected_closed_form_for_mml():
    """MML/EAP lz* uses the population-prior correction when terms exist."""
    pr_observed = (0.55, 0.62, 0.48)
    item_entropy = (-0.95, -0.80, -1.02)
    item_var_logp = (0.18, 0.22, 0.27)
    item_cov = (0.05, -0.04, 0.08)
    score_info = (0.40, 0.55, 0.45)
    obs_score_deriv = (0.10, -0.08, 0.12)
    sigma = 1.25
    obs = _build_synthetic_obs(
        pr_observed=pr_observed,
        item_entropy=item_entropy,
        item_var_logp=item_var_logp,
        extra_cols={
            "ItemLogPScoreCov": item_cov,
            "ScoreInformation": score_info,
            "ObservedScoreDerivative": obs_score_deriv,
        },
    )

    log_p = np.log(np.array(pr_observed))
    ll = float(np.sum(log_p))
    e_ll = float(np.sum(item_entropy))
    var_ll = float(np.sum(item_var_logp))
    info_total = float(np.sum(score_info))
    cov_total = float(np.sum(item_cov))
    score_sum = float(np.sum(obs_score_deriv))
    p = 1.0 / (sigma ** 2)
    denom = info_total + p
    c_n = cov_total / denom
    corrected_var = var_ll - (cov_total ** 2) * (info_total + 2.0 * p) / (denom ** 2)
    expected_lz_star = (ll - e_ll - c_n * score_sum) / np.sqrt(corrected_var)

    res = {
        "config": {"method": "MML", "population_prior_sd": sigma},
        "prep": {"levels": {"Person": ["P1"]}},
    }
    out = app.compute_person_fit_indices(res, obs=obs)
    row = out.iloc[0]
    assert row["lz_star_status"] == "computed_eap_population_corrected"
    assert row["lz_star"] == pytest.approx(expected_lz_star, abs=1e-12)
    assert row["lz_star_c"] == pytest.approx(c_n, abs=1e-12)
    assert row["lz_star_variance"] == pytest.approx(corrected_var, abs=1e-12)
    assert row["ReportIndex"] == "lz_star"
    assert row["ReportValue"] == pytest.approx(expected_lz_star, abs=1e-12)
    assert "Sinharay" in row["ReportCaveat"]


def test_lz_star_eap_reduces_to_jml_as_sigma_grows():
    """As sigma -> inf (p -> 0) the EAP correction collapses onto the JML formula."""
    extra = {
        "ItemLogPScoreCov": (0.05, -0.04, 0.08),
        "ScoreInformation": (0.40, 0.55, 0.45),
        "ObservedScoreDerivative": (0.10, -0.08, 0.12),
    }
    obs = _build_synthetic_obs(
        pr_observed=(0.55, 0.62, 0.48),
        item_entropy=(-0.95, -0.80, -1.02),
        item_var_logp=(0.18, 0.22, 0.27),
        extra_cols=extra,
    )
    jml = app.compute_person_fit_indices(
        {"config": {"method": "JMLE"}, "prep": {"levels": {"Person": ["P1"]}}}, obs=obs
    ).iloc[0]
    eap = app.compute_person_fit_indices(
        {"config": {"method": "MML", "population_prior_sd": 1.0e6},
         "prep": {"levels": {"Person": ["P1"]}}}, obs=obs
    ).iloc[0]
    assert eap["lz_star_status"] == "computed_eap_population_corrected"
    assert jml["lz_star_status"] == "computed_jml_conditional_calibration"
    assert eap["lz_star"] == pytest.approx(jml["lz_star"], abs=1e-6)
    assert eap["lz_star_c"] == pytest.approx(jml["lz_star_c"], abs=1e-6)
    assert eap["lz_star_variance"] == pytest.approx(jml["lz_star_variance"], abs=1e-6)


def test_lz_star_eap_unavailable_when_population_sd_missing():
    """MML without a usable population SD falls back to unadjusted lz."""
    obs = _build_synthetic_obs(
        extra_cols={
            "ItemLogPScoreCov": (0.05, -0.04, 0.08),
            "ScoreInformation": (0.40, 0.55, 0.45),
            "ObservedScoreDerivative": (0.10, -0.08, 0.12),
        }
    )
    res = {"config": {"method": "MML", "population_prior_sd": 0.0},
           "prep": {"levels": {"Person": ["P1"]}}}
    row = app.compute_person_fit_indices(res, obs=obs).iloc[0]
    assert row["lz_star_status"] == "eap_population_sd_unavailable"
    assert pd.isna(row["lz_star"])
    assert np.isfinite(row["lz"])
    assert row["ReportIndex"] == "lz"
    assert row["ReportValue"] == pytest.approx(row["lz"])


def test_lz_star_status_diagnostics_missing_when_obs_lacks_columns():
    """Without the lz*-specific columns the pipeline reports the status
    without raising."""
    obs = _build_synthetic_obs()  # no extra_cols => no lz* terms
    res = {"config": {"method": "JMLE"}, "prep": {"levels": {"Person": ["P1"]}}}
    out = app.compute_person_fit_indices(res, obs=obs)
    row = out.iloc[0]
    assert row["lz_star_status"] == "diagnostics_missing_snijders_terms"
    assert pd.isna(row["lz_star"])
    assert row["ReportIndex"] == "lz"


def test_lz_star_status_insufficient_information_on_degenerate_obs():
    """Degenerate inputs (zero info, zero variance) must report
    ``insufficient_information`` and leave the statistic NaN."""
    obs = _build_synthetic_obs(
        pr_observed=(0.5,),
        item_entropy=(-0.7,),
        item_var_logp=(0.0,),  # zero variance => can't standardise
        extra_cols={
            "ItemLogPScoreCov": [0.0],
            "ScoreInformation": [0.0],  # zero info => corrected_var undefined
            "ObservedScoreDerivative": [0.0],
        },
    )
    res = {"config": {"method": "JMLE"}, "prep": {"levels": {"Person": ["P1"]}}}
    out = app.compute_person_fit_indices(res, obs=obs)
    row = out.iloc[0]
    assert row["lz_star_status"] == "insufficient_information"
    assert pd.isna(row["lz_star"])


# -----------------------------------------------------------------------------
# Threshold flag logic
# -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "lz_value,expected_5,expected_1,expected_level",
    [
        (0.5, False, False, "none"),
        (-1.5, False, False, "none"),
        (2.0, True, False, "5pct"),    # |z| ∈ (1.96, 2.58]
        (-2.3, True, False, "5pct"),
        (3.0, True, True, "1pct"),     # |z| > 2.58
        (-3.5, True, True, "1pct"),
    ],
)
def test_lz_flag_levels_follow_normal_quantile_rule(
    lz_value, expected_5, expected_1, expected_level
):
    """The flag columns must use |z| > qnorm(0.975) and qnorm(0.995) on
    the lz value (falling back when lz* is unavailable)."""
    log_p = np.array([lz_value / np.sqrt(1.0)])
    obs = pd.DataFrame(
        {
            "Person": ["P1"],
            # Choose PrObserved/Entropy/Var such that lz = lz_value:
            # log P = lz * sqrt(Var) + E[log P], with Var = 1, E[log P] = 0
            # so log P = lz_value, PrObserved = exp(lz_value).
            "PrObserved": [float(np.exp(lz_value))],
            "ItemEntropy": [0.0],
            "ItemVarLogP": [1.0],
        }
    )
    res = {"config": {"method": "MML"}, "prep": {"levels": {"Person": ["P1"]}}}
    out = app.compute_person_fit_indices(res, obs=obs)
    row = out.iloc[0]
    assert row["lz"] == pytest.approx(lz_value, abs=1e-12)
    assert bool(row["lz_flag_5pct"]) is expected_5
    assert bool(row["lz_flag_1pct"]) is expected_1
    assert row["ReportFlagLevel"] == expected_level
    if expected_level == "1pct":
        assert row["ReviewStatus"] == "review_1pct"
    elif expected_level == "5pct":
        assert row["ReviewStatus"] == "review_5pct"
    elif expected_level == "none":
        assert row["ReviewStatus"] == "not_flagged"


# -----------------------------------------------------------------------------
# ReportIndex selection: prefers lz_star when available
# -----------------------------------------------------------------------------


def test_report_index_prefers_lz_star_when_snijders_succeeded():
    obs = _build_synthetic_obs(
        extra_cols={
            "ItemLogPScoreCov": [0.05, -0.04, 0.08],
            "ScoreInformation": [0.40, 0.55, 0.45],
            "ObservedScoreDerivative": [0.10, -0.08, 0.12],
        }
    )
    res = {"config": {"method": "JMLE"}, "prep": {"levels": {"Person": ["P1"]}}}
    out = app.compute_person_fit_indices(res, obs=obs)
    row = out.iloc[0]
    assert row["ReportIndex"] == "lz_star"
    assert row["ReportValue"] == pytest.approx(row["lz_star"])
    # Caveat must mention the conditional-calibration scope.
    assert "Snijders" in row["ReportCaveat"]


# -----------------------------------------------------------------------------
# Missing required columns => clear error
# -----------------------------------------------------------------------------


def test_compute_person_fit_indices_raises_on_missing_required_columns():
    obs = pd.DataFrame(
        {
            "Person": ["P1"],
            "PrObserved": [0.5],
            "ItemEntropy": [-0.7],
            # ItemVarLogP missing
        }
    )
    res = {"config": {"method": "JMLE"}, "prep": {"levels": {"Person": ["P1"]}}}
    with pytest.raises(ValueError, match="missing required person-fit columns"):
        app.compute_person_fit_indices(res, obs=obs)


# -----------------------------------------------------------------------------
# End-to-end: small JMLE fit produces both lz and lz* with correct status
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def small_gpcm_jmle_fit():
    """30 persons x 2 raters x 3 criteria GPCM JMLE fit.

    Uses the same data shape as the bias-inference fixture but a JMLE
    estimator, so that lz* receives the JML-method path rather than
    the MML fall-back. ``maxit`` is kept small to keep the fixture
    cost low; the test only needs a converged enough fit to populate
    all six per-observation diagnostic columns.
    """
    rng = np.random.default_rng(20260514)
    n_person, n_rater, n_criterion = 30, 2, 3
    theta = rng.normal(0, 1, n_person)
    rater_eff = np.array([-0.4, 0.4])
    crit_eff = np.array([-0.5, 0.0, 0.5])
    rows = []
    for i, theta_i in enumerate(theta):
        for j in range(n_rater):
            for k in range(n_criterion):
                eta = theta_i - rater_eff[j] - crit_eff[k]
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                score = int(rng.uniform() < p1) + int(rng.uniform() < p2)
                rows.append(
                    {
                        "Person": f"P{i+1:02d}",
                        "Rater": f"R{j+1}",
                        "Criterion": f"C{k+1}",
                        "Score": score,
                    }
                )
    df = pd.DataFrame(rows)
    res = app.mfrm_estimate(
        data=df,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="GPCM",
        method="JMLE",
        step_facet="Criterion",
        maxit=25,
        reltol=1e-4,
    )
    return res


def test_compute_person_fit_indices_on_small_jmle_fit(small_gpcm_jmle_fit):
    out = app.compute_person_fit_indices(small_gpcm_jmle_fit)
    assert len(out) == 30
    # Every person must have a finite lz; lz* succeeds on the JMLE path
    # for the majority of persons (a few may trip the corrected_var
    # guard rail on synthetic data, which is expected).
    assert out["lz"].notna().all()
    success_mask = out["lz_star_status"].eq("computed_jml_conditional_calibration")
    assert success_mask.mean() >= 0.5  # most persons converge
    # On the success path lz and lz* must both be finite.
    finite_subset = out[success_mask]
    assert finite_subset["lz"].notna().all()
    assert finite_subset["lz_star"].notna().all()
    # ReportIndex tracks the available statistic.
    assert (out.loc[success_mask, "ReportIndex"] == "lz_star").all()
    assert (out.loc[~success_mask, "ReportIndex"] == "lz").all()


def test_compute_person_fit_indices_columns_are_stable(small_gpcm_jmle_fit):
    """The returned column order is part of the public contract because
    downstream UI and export code reads by column name."""
    out = app.compute_person_fit_indices(small_gpcm_jmle_fit)
    expected_cols = [
        "Person", "N", "LogLik", "lz", "lz_star", "lz_star_status",
        "lz_star_c", "lz_star_variance",
        "lz_flag_5pct", "lz_flag_1pct",
        "lz_star_flag_5pct", "lz_star_flag_1pct",
        "ReportIndex", "ReportValue", "ReportFlagLevel", "ReportFlag",
        "ReviewStatus", "ReviewReason", "ReportCaveat",
    ]
    assert list(out.columns) == expected_cols


def test_compute_person_fit_indices_caveats_describe_pathway(small_gpcm_jmle_fit):
    out = app.compute_person_fit_indices(small_gpcm_jmle_fit)
    success_mask = out["lz_star_status"].eq("computed_jml_conditional_calibration")
    assert out.loc[success_mask, "ReportCaveat"].str.contains("Snijders").all()
    fallback_caveats = out.loc[~success_mask, "ReportCaveat"]
    assert fallback_caveats.str.contains("lz_star_status").all()
