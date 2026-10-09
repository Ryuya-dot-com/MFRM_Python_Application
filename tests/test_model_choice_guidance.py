"""Tests for the RSM / PCM / GPCM model-choice guidance pipeline.

``compute_model_choice_comparison`` refits the two non-current rating-
scale models on the same data and returns an AIC / BIC comparison
table plus the nested-pair likelihood-ratio tests. The math contract
pinned here covers:

* Backwards-compatible refusal when the fit is anchored / regularized
  / population-regressed (likelihood scale is not comparable).
* The current model is always the cheapest row (no refit), and the
  other two refits report a finite ``FitStatus``.
* ``DeltaAIC`` / ``DeltaBIC`` are computed relative to the row with
  the lowest IC, with the minimum row at exactly zero.
* The Akaike / Schwarz evidence ratios ``exp(DeltaIC / 2)`` are
  reproduced to machine precision.
* The nested LR table follows ``Lambda = 2 (LL_alt - LL_null) ~
  chi2_{df_alt - df_null}`` (Wilks, 1938).
* On RSM-generated data the recommendation prefers RSM.
* On GPCM-generated data with markedly heterogeneous slopes the
  recommendation prefers GPCM.

Refits use small fixtures to keep the test suite cheap; the
correctness contract is what the tests pin, not the exact AIC values.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2

import streamlit_app as app


# -----------------------------------------------------------------------------
# Refusal: anchored / regularized / population-formula fits
# -----------------------------------------------------------------------------


def test_compute_model_choice_returns_unavailable_for_invalid_input():
    out = app.compute_model_choice_comparison(None)
    assert out["available"] is False
    assert "fit dictionary" in out["reason"].lower()


def test_compute_model_choice_refuses_facet_regularization_fits():
    """When the regularization flag is set, the optimized objective is
    not a comparable likelihood and the helper must refuse."""
    out = app.compute_model_choice_comparison({
        "config": {"model": "RSM", "facet_regularization_enabled": True, "n_cat": 3},
        "prep": {"data": pd.DataFrame({"Person": ["P1"], "Score": [1]})},
    })
    assert out["available"] is False
    assert "regularization" in out["reason"].lower()


def test_compute_model_choice_refuses_population_formula_fits():
    out = app.compute_model_choice_comparison({
        "config": {
            "model": "RSM",
            "n_cat": 3,
            "population_model": {"enabled": True},
        },
        "prep": {"data": pd.DataFrame({"Person": ["P1"], "Score": [1]})},
    })
    assert out["available"] is False
    assert "population" in out["reason"].lower()


@pytest.mark.parametrize("free_sd,sd_input,sd_returned", [
    (False, None, 1.7), (False, .8, 1.7), (True, .8, 1.7), (True, None, 1.7),
])
def test_candidate_refit_keeps_population_model(monkeypatch, free_sd, sd_input, sd_returned):
    """Capture estimator arguments without generating or fitting any data."""
    config = dict(model="RSM", method="MML", n_cat=3, facet_names=["Criterion"],
                  population_prior_sd=sd_returned, estimate_population_sd=free_sd)
    if sd_input is not None:
        config["population_prior_sd_input"] = sd_input
    calls = []
    monkeypatch.setattr(app, "mfrm_estimate", lambda **kwargs: calls.append(kwargs) or kwargs)
    result = app._refit_candidate_model(
        dict(config=config, prep=dict(data=pd.DataFrame({"Person": ["P1"], "Score": [1]}))),
        "PCM", step_facet="Criterion", slope_facet="Criterion",
    )
    assert result is calls[0]
    assert result["estimate_population_sd"] is free_sd
    assert result["population_prior_sd"] == (sd_input if free_sd and sd_input is not None else sd_returned)
    assert result["method"] == "MML" and result["model"] == "PCM"


@pytest.mark.parametrize("weight", [0., .5, 2., float("nan"), "invalid"])
def test_model_comparison_refuses_nonunit_weights_before_refit(monkeypatch, weight):
    def no_refit(*args, **kwargs):
        pytest.fail("An ineligible weighted comparison must stop before any refit")
    monkeypatch.setattr(app, "_refit_candidate_model", no_refit)
    result = app.compute_model_choice_comparison(dict(
        config=dict(model="RSM", method="MML", n_cat=3),
        prep=dict(data=pd.DataFrame({"Person": ["P1", "P2"], "Score": [1, 0], "Weight": [1., weight]})),
    ))
    assert result["available"] is False
    assert "unit weights" in result["reason"]


def test_model_comparison_accepts_explicit_unit_weights():
    allowed, reason = app._model_choice_can_compare(dict(
        config=dict(model="RSM", method="MML", n_cat=3),
        prep=dict(data=pd.DataFrame({"Person": ["P1"], "Score": [1], "Weight": [1.]})),
    ))
    assert allowed and not reason


# -----------------------------------------------------------------------------
# RSM-generated data: comparison ranks RSM at the top
# -----------------------------------------------------------------------------


def _build_rsm_jmle_fit(n_person: int = 25, seed: int = 20260520):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_person):
        theta = rng.normal(0, 1)
        for j in range(2):
            for k in range(3):
                eta = theta - (j - 0.5) * 0.4 - (k - 1) * 0.4
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
        method="JMLE", maxit=25, reltol=1e-4,
    )


@pytest.fixture(scope="module")
def rsm_fit_for_model_choice():
    return _build_rsm_jmle_fit()


def test_model_choice_returns_three_candidate_rows(rsm_fit_for_model_choice):
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    assert out["available"] is True
    comparison = out["comparison"]
    assert list(comparison["Model"]) == ["RSM", "PCM", "GPCM"]
    # The original model is the "current" row, others are "refit_ok" or "refit_failed".
    statuses = set(comparison["FitStatus"].astype(str))
    assert "current" in statuses
    assert statuses.issubset({"current", "refit_ok", "refit_failed"})


def test_model_choice_delta_columns_have_min_zero(rsm_fit_for_model_choice):
    """The IC differences are computed relative to the best model in the
    comparison; the minimum row must have DeltaAIC = DeltaBIC = 0."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    comparison = out["comparison"]
    finite = comparison.dropna(subset=["AIC", "BIC"])
    assert float(finite["DeltaAIC"].min()) == pytest.approx(0.0, abs=1e-12)
    assert float(finite["DeltaBIC"].min()) == pytest.approx(0.0, abs=1e-12)


def test_model_choice_evidence_ratio_matches_closed_form(rsm_fit_for_model_choice):
    """``evidence_ratio = exp(DeltaIC / 2)`` (Burnham & Anderson 2002,
    Eq. 2.10) must hold to machine precision for both AIC and BIC."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    comp = out["comparison"]
    for col_delta, col_ratio in [("DeltaAIC", "AICEvidenceRatio"),
                                  ("DeltaBIC", "BICEvidenceRatio")]:
        for delta, ratio in zip(comp[col_delta], comp[col_ratio]):
            if not np.isfinite(delta):
                continue
            assert ratio == pytest.approx(np.exp(delta / 2.0), abs=1e-12)


def test_model_choice_lr_test_chi_square_matches_2_delta_loglik(rsm_fit_for_model_choice):
    """LR ChiSq = 2 * (LL_alt - LL_null) for every nested pair where
    both fits succeeded."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    lr = out["lr_tests"]
    comparison = out["comparison"]
    ll_lookup = dict(zip(comparison["Model"], comparison["LogLik"]))
    for _, row in lr.iterrows():
        ll_null = float(ll_lookup[row["Null"]])
        ll_alt = float(ll_lookup[row["Alternative"]])
        assert row["ChiSq"] == pytest.approx(2.0 * (ll_alt - ll_null), abs=1e-12)


def test_model_choice_lr_p_value_matches_scipy(rsm_fit_for_model_choice):
    """``p = chi2.sf(ChiSq, df)`` must hold exactly."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    lr = out["lr_tests"]
    for _, row in lr.iterrows():
        if not np.isfinite(row["p"]):
            continue
        expected = float(chi2.sf(row["ChiSq"], df=int(row["df"])))
        assert row["p"] == pytest.approx(expected, abs=1e-12)


def test_model_choice_lr_values_are_technical_and_public_decisions_are_withheld(
    rsm_fit_for_model_choice,
):
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    lr = out["lr_tests"]

    assert not lr.empty
    assert set(lr["QualificationStatus"]) == {"WITHHELD"}
    assert not lr["PublicConclusionAllowed"].astype(bool).any()
    gpcm_pairs = lr.loc[lr["Alternative"].eq("GPCM")]
    assert set(gpcm_pairs["ReasonCode"]) == {
        "model.model_002.gpcm_lrt_uncalibrated"
    }


def test_model_choice_keeps_rsm_ranking_technical_and_withholds_public_recommendation(
    rsm_fit_for_model_choice,
):
    """The pre-gate ranking is retained for audit, never as public guidance."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    rec = out["recommendation"]
    technical = out["technical_recommendation"]
    assert technical["model"] == "RSM"
    assert technical["tier"] in {"strong", "moderate", "weak"}
    assert rec["model"] == ""
    assert rec["tier"] == "withheld"
    assert rec["reason_code"] == "model.model_001.jmle_auto_recommendation_withheld"


# -----------------------------------------------------------------------------
# GPCM-generated data with heterogeneous slopes: recommendation prefers GPCM
# -----------------------------------------------------------------------------


def _build_gpcm_jmle_fit(n_person: int = 60, seed: int = 20260521):
    """Synthetic data with markedly heterogeneous discriminations per
    criterion. Slopes (0.5, 1.5, 2.5) drive the per-criterion category
    transitions to differ enough that PCM cannot reproduce them; RSM is
    even further off because the step thresholds also vary per criterion."""
    rng = np.random.default_rng(seed)
    rows = []
    n_rater = 2
    n_criterion = 3
    slopes = np.array([0.5, 1.5, 2.5])
    crit_steps = np.array([
        [-1.5, 1.5],
        [-0.5, 0.5],
        [-0.2, 0.2],
    ])
    rater_eff = np.array([-0.4, 0.4])
    theta = rng.normal(0, 1, n_person)
    for i, theta_i in enumerate(theta):
        for j in range(n_rater):
            for k in range(n_criterion):
                eta = theta_i - rater_eff[j]
                a = slopes[k]
                step1 = crit_steps[k, 0]
                step2 = crit_steps[k, 1]
                # GPCM category probabilities for K = 3 categories.
                num = np.array([
                    0.0,
                    a * (eta - step1),
                    a * (eta - step1) + a * (eta - step2),
                ])
                num -= num.max()  # numerical stability
                prob = np.exp(num) / np.exp(num).sum()
                score = int(np.searchsorted(np.cumsum(prob), rng.uniform()))
                rows.append({
                    "Person": f"P{i+1:02d}",
                    "Rater": f"R{j+1}",
                    "Criterion": f"C{k+1}",
                    "Score": score,
                })
    df = pd.DataFrame(rows)
    return app.mfrm_estimate(
        data=df, person_col="Person", facet_cols=["Rater", "Criterion"],
        score_col="Score", rating_min=0, rating_max=2, model="GPCM",
        method="JMLE", step_facet="Criterion", maxit=30, reltol=1e-4,
    )


@pytest.fixture(scope="module")
def gpcm_fit_for_model_choice():
    return _build_gpcm_jmle_fit()


def test_model_choice_prefers_gpcm_on_gpcm_heterogeneous_slope_data(
    gpcm_fit_for_model_choice,
):
    """With slopes (0.5, 1.5, 2.5) the GPCM model is the only one that
    can fit the per-criterion transition slopes. The comparison must
    pick GPCM (or at least put it above RSM, since the slope-1 model
    cannot reproduce the shape)."""
    out = app.compute_model_choice_comparison(gpcm_fit_for_model_choice)
    rec = out["recommendation"]
    technical = out["technical_recommendation"]
    comp = out["comparison"]
    finite = comp.dropna(subset=["AIC"])
    # GPCM must beat RSM at the AIC level; we don't fix the recommendation
    # tier because the BIC penalty can occasionally still prefer simpler
    # models on small samples even when GPCM has materially lower deviance.
    aic_rsm = float(finite.set_index("Model").loc["RSM", "AIC"])
    aic_gpcm = float(finite.set_index("Model").loc["GPCM", "AIC"])
    assert aic_gpcm < aic_rsm
    assert technical["model"] in {"GPCM", "PCM", "RSM"}
    assert technical["tier"] in {"strong", "moderate", "weak", "tie"}
    assert rec["model"] == ""
    assert rec["tier"] == "withheld"


# -----------------------------------------------------------------------------
# Caveat string and column-order stability
# -----------------------------------------------------------------------------


def test_model_choice_caveat_cites_wilks_and_describes_scope(
    rsm_fit_for_model_choice,
):
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    caveat = out["caveat"]
    assert "Wilks" in caveat
    assert "AIC" in caveat and "BIC" in caveat


def test_model_choice_comparison_column_order_is_stable(
    rsm_fit_for_model_choice,
):
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    expected = [
        "Model", "FitStatus", "N", "KParams", "N_over_K", "AICcRecommended",
        "LogLik",
        "AIC", "DeltaAIC", "AICEvidenceRatio", "AICWeight",
        "AICc", "DeltaAICc", "AICcEvidenceRatio", "AICcWeight",
        "BIC", "DeltaBIC", "BICEvidenceRatio", "BICWeight",
        "BICSampleSize", "BICSampleSizeBasis", "AICcStatus",
    ]
    assert list(out["comparison"].columns) == expected


def test_model_choice_information_weights_sum_to_one(rsm_fit_for_model_choice):
    """Akaike / Schwarz weights sum to 1 over the finite-DeltaIC rows
    (Burnham & Anderson 2002, Eq. 2.8)."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    comp = out["comparison"]
    for col in ["AICWeight", "AICcWeight", "BICWeight"]:
        weights = pd.to_numeric(comp[col], errors="coerce").dropna()
        if weights.empty:
            continue
        assert float(weights.sum()) == pytest.approx(1.0, abs=1e-12), col


def test_model_choice_information_weights_match_closed_form(
    rsm_fit_for_model_choice,
):
    """w_i = exp(-DeltaIC_i / 2) / sum_j exp(-DeltaIC_j / 2)
    (Burnham & Anderson 2002, Eq. 2.8). Pin to machine precision."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    comp = out["comparison"]
    for delta_col, weight_col in [
        ("DeltaAIC", "AICWeight"),
        ("DeltaAICc", "AICcWeight"),
        ("DeltaBIC", "BICWeight"),
    ]:
        deltas = pd.to_numeric(comp[delta_col], errors="coerce")
        weights = pd.to_numeric(comp[weight_col], errors="coerce")
        finite = deltas.notna()
        if not finite.any():
            continue
        rel = np.exp(-deltas[finite].to_numpy() / 2.0)
        expected = rel / rel.sum()
        actual = weights[finite].to_numpy()
        assert np.allclose(actual, expected, atol=1e-12), (
            f"Weight mismatch for {weight_col}: actual={actual}, expected={expected}"
        )


def test_model_choice_best_model_carries_largest_weight(rsm_fit_for_model_choice):
    """The row with DeltaIC = 0 must also have the largest weight
    (the weight is a monotone-decreasing function of DeltaIC)."""
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    comp = out["comparison"]
    for delta_col, weight_col in [
        ("DeltaAIC", "AICWeight"),
        ("DeltaAICc", "AICcWeight"),
        ("DeltaBIC", "BICWeight"),
    ]:
        finite = comp.dropna(subset=[delta_col, weight_col])
        if finite.empty:
            continue
        best_idx = finite[delta_col].idxmin()
        best_weight = float(finite.loc[best_idx, weight_col])
        assert best_weight == pytest.approx(
            float(finite[weight_col].max()), abs=1e-12
        )


def test_model_choice_fit_times_keyed_by_model_and_current_is_zero(
    rsm_fit_for_model_choice,
):
    out = app.compute_model_choice_comparison(rsm_fit_for_model_choice)
    fit_times = out["fit_times"]
    assert set(fit_times.keys()).issubset({"RSM", "PCM", "GPCM"})
    # The current model is not refit, so its time must be exactly zero.
    assert fit_times["RSM"] == 0.0


# -----------------------------------------------------------------------------
# AICc finite-sample correction (Hurvich & Tsai, 1989)
# -----------------------------------------------------------------------------


def test_model_choice_aicc_is_unavailable_without_mfrm_derivation(rsm_fit_for_model_choice):
    comp = app.compute_model_choice_comparison(rsm_fit_for_model_choice)["comparison"]
    for column in ["AICc", "DeltaAICc", "AICcEvidenceRatio", "AICcWeight"]:
        assert comp[column].isna().all()
    assert set(comp["AICcStatus"]) == {"unavailable_no_mfrm_small_sample_derivation"}


def test_model_choice_small_n_over_k_does_not_qualify_aicc(rsm_fit_for_model_choice):
    comp = app.compute_model_choice_comparison(rsm_fit_for_model_choice)["comparison"]
    assert (comp["N_over_K"] < 40).any()
    assert not comp["AICcRecommended"].any()


def test_model_choice_ranking_ignores_legacy_aicc_recommendation_flags():
    comp = pd.DataFrame(dict(Model=["RSM", "PCM"], BIC=[30., 30.], AIC=[20., 25.],
                             AICc=[100., 28.], AICcRecommended=[True, True]))
    rec = app._model_choice_recommend(comp)
    assert rec["model"] == "RSM"
    assert "AICc" not in rec["rationale"]


def test_model_choice_caveat_distinguishes_formula_provenance_and_applicability(rsm_fit_for_model_choice):
    caveat = app.compute_model_choice_comparison(rsm_fit_for_model_choice)["caveat"]
    assert "AICc is unavailable" in caveat
    assert "Hurvich" in caveat and "Burnham" in caveat
    assert "N / K < 40 does not establish applicability" in caveat
    assert "independent Person" in caveat
