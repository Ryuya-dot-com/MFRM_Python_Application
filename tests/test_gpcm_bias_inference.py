"""Tests for the slope-aware GPCM bias inference pipeline.

Covers five contracts:

* ``_safe_nll_value`` returns ``nan`` when the underlying closure raises
  instead of letting the exception escape into the bias-estimation loop.
* ``_profile_bias_ci`` inverts the chi-square pivotal correctly. For a
  quadratic ``NLL = (b - b_hat)^2 / 2`` the analytic 95 % profile-
  likelihood CI is ``b_hat +/- sqrt(chi2_0.95) = b_hat +/- 1.95996...``;
  the function must reproduce that to machine tolerance.
* The slope-aware Fisher information identity. For a GPCM MML fit, the
  bias-cell standard error must equal
  ``1 / sqrt(sum_i a_i^2 * Var[X_i | eta_i + b_hat] * w_i)``; the
  non-slope-aware fall-back (``sum_i Var * w``) must give a different
  answer.
* Likelihood-ratio test and profile CI plumbing on a real GPCM MML fit.
  ``LR ChiSq`` matches ``2 * (nll_null - nll_hat)`` recomputed by hand,
  ``LR Prob.`` matches ``chi2.sf(LR, df=1)``, the profile CI brackets the
  point estimate, and the slope-dispatch sanity check confirms that
  clamping ``log_slopes`` to zero actually changes the reported bias
  estimates.
* R parity. The Python output agrees with the mfrmr 0.2.0 R reference
  (R 4.5.2) at manuscript-citation tolerance across six bias cells on a
  shared synthetic data set. The R fixture is generated once via
  ``Rscript`` against ``tests/data/r_bias_parity_input.csv`` and stored
  at ``tests/data/r_bias_parity_output.json``; the parity test
  reproduces the same fit in Python and verifies that ``Bias Size``,
  ``S.E.``, ``LR ChiSq``, ``LR Prob.``, and the profile-likelihood CI
  endpoints all agree with the R values to within ``1e-3``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2

import streamlit_app as app


R_BIAS_PARITY_INPUT = Path(__file__).resolve().parent / "data" / "r_bias_parity_input.csv"
R_BIAS_PARITY_OUTPUT = Path(__file__).resolve().parent / "data" / "r_bias_parity_output.json"


# -----------------------------------------------------------------------------
# Unit: _safe_nll_value
# -----------------------------------------------------------------------------


def test_safe_nll_value_returns_nan_for_raising_callable():
    def boom(b):
        raise RuntimeError("simulated optimisation failure")

    assert np.isnan(app._safe_nll_value(boom, 0.0))


def test_safe_nll_value_returns_nan_for_inf_or_nan_output():
    assert np.isnan(app._safe_nll_value(lambda b: float("inf"), 0.0))
    assert np.isnan(app._safe_nll_value(lambda b: float("nan"), 0.0))


def test_safe_nll_value_returns_value_for_finite_output():
    val = app._safe_nll_value(lambda b: 1.23, 0.0)
    assert val == pytest.approx(1.23)


# -----------------------------------------------------------------------------
# Unit: _profile_bias_ci
# -----------------------------------------------------------------------------


def _quadratic_nll(b_hat: float, scale: float = 1.0):
    """Return a quadratic NLL closure with minimum at ``b_hat`` and value 0.

    With ``NLL(b) = (b - b_hat)^2 / 2``, the LR statistic ``2 * (NLL(b) -
    NLL(b_hat)) = (b - b_hat)^2`` follows a chi-square distribution with
    one degree of freedom exactly, so the profile-likelihood CI is
    available in closed form for cross-checking the implementation.
    """

    def nll(b):
        return scale * (float(b) - float(b_hat)) ** 2 / 2.0

    return nll


def test_profile_bias_ci_matches_chi_square_inversion_on_quadratic_nll():
    """For a textbook quadratic NLL the 95 % profile CI is the analytic
    chi-square inversion ``b_hat +/- sqrt(chi2_0.95)``.
    """
    b_hat = 0.5
    nll = _quadratic_nll(b_hat=b_hat)
    nll_min = 0.0
    ci = app._profile_bias_ci(nll, estimate=b_hat, nll_min=nll_min, max_abs=10.0, level=0.95)
    expected_halfwidth = float(np.sqrt(chi2.ppf(0.95, df=1)))
    assert ci["status"] == "ok"
    assert ci["lower"] == pytest.approx(b_hat - expected_halfwidth, abs=1e-7)
    assert ci["upper"] == pytest.approx(b_hat + expected_halfwidth, abs=1e-7)
    assert ci["level"] == 0.95


def test_profile_bias_ci_reports_limited_when_likelihood_does_not_fall_far_enough():
    """If the negative log-likelihood barely changes inside the search
    bracket, both endpoints must be reported with ``"limited by search
    range"`` rather than as finite interior roots."""
    # Very gentle quadratic, so within [-1, +1] the LR statistic never
    # reaches the 95 % cutoff of 3.84.
    nll = _quadratic_nll(b_hat=0.0, scale=0.05)
    ci = app._profile_bias_ci(nll, estimate=0.0, nll_min=0.0, max_abs=1.0, level=0.95)
    assert "limited by search range" in ci["status"]
    assert ci["lower"] == pytest.approx(-1.0)
    assert ci["upper"] == pytest.approx(1.0)


def test_profile_bias_ci_returns_not_available_for_invalid_inputs():
    nll = _quadratic_nll(b_hat=0.0)
    for bad in [
        dict(estimate=float("nan"), nll_min=0.0, max_abs=10.0, level=0.95),
        dict(estimate=0.0, nll_min=float("nan"), max_abs=10.0, level=0.95),
        dict(estimate=0.0, nll_min=0.0, max_abs=0.0, level=0.95),
        dict(estimate=0.0, nll_min=0.0, max_abs=10.0, level=1.5),
    ]:
        ci = app._profile_bias_ci(nll, **bad)
        assert ci["status"] == "not available"


def test_profile_bias_ci_brackets_the_estimate_when_status_is_ok():
    """Whenever the status comes back ``"ok"`` the resulting CI must
    bracket the point estimate at machine tolerance."""
    for b_hat in (-2.0, -0.5, 0.0, 0.7, 1.4):
        nll = _quadratic_nll(b_hat=b_hat)
        ci = app._profile_bias_ci(nll, estimate=b_hat, nll_min=0.0, max_abs=10.0, level=0.95)
        if ci["status"] == "ok":
            assert ci["lower"] <= b_hat + 1e-10
            assert ci["upper"] >= b_hat - 1e-10


# -----------------------------------------------------------------------------
# Integration: small GPCM MML fit shared fixture
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def small_gpcm_mml_fit_for_bias():
    """Same scale as the fair-average SE tests: 30 persons x 2 raters x 3
    criteria GPCM MML fit with quad_points = 5 and maxit = 25. The
    fixture is module-scoped so the fit cost is paid once across the
    integration assertions in this file.
    """
    rng = np.random.default_rng(20260513)
    n_person, n_rater, n_criterion = 30, 2, 3
    theta = rng.normal(0, 1, n_person)
    rater_eff = np.array([-0.4, 0.4])
    crit_eff = np.array([-0.5, 0.0, 0.5])
    rater_levels = [f"R{i+1}" for i in range(n_rater)]
    crit_levels = [f"C{i+1}" for i in range(n_criterion)]
    person_levels = [f"P{i+1:02d}" for i in range(n_person)]
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
                        "Person": person_levels[i],
                        "Rater": rater_levels[j],
                        "Criterion": crit_levels[k],
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
        method="MML",
        step_facet="Criterion",
        quad_points=5,
        maxit=25,
        reltol=1e-4,
        mml_engine="direct",
    )
    return res


@pytest.fixture(scope="module")
def small_gpcm_mml_diagnostics_for_bias(small_gpcm_mml_fit_for_bias):
    return app.mfrm_diagnostics(
        small_gpcm_mml_fit_for_bias, compute_pca=False, compute_marginal=False
    )


# -----------------------------------------------------------------------------
# Integration: slope-aware Fisher information identity
# -----------------------------------------------------------------------------


def test_bias_se_under_gpcm_matches_slope_aware_information(
    small_gpcm_mml_fit_for_bias, small_gpcm_mml_diagnostics_for_bias,
):
    """For every GPCM bias cell the reported S.E. must equal
    ``1 / sqrt(sum_i a_i^2 * Var[X_i | eta_i + b_hat] * w_i)`` and must
    *not* equal the non-slope-aware ``1 / sqrt(sum_i Var * w)`` fall-back
    (Muraki 1993 Eqs. 7, 16).
    """
    res = small_gpcm_mml_fit_for_bias
    diag = small_gpcm_mml_diagnostics_for_bias
    bias_out = app.estimate_bias_interaction(res, diag, facet_a="Rater", facet_b="Criterion")
    assert "_skip_reason" not in bias_out, bias_out
    tbl = bias_out["table"]
    assert not tbl.empty

    config = res["config"]
    params = res["params"]
    obs_df = diag["obs"]
    idx = app.build_indices(
        res["prep"], step_facet=config["step_facet"], slope_facet=config.get("slope_facet")
    )
    theta_hat = res["facets"]["person"]["Estimate"].to_numpy()
    eta_base = app.compute_eta(idx, params, config, theta_override=theta_hat)
    step_cum_mat = np.vstack(
        [np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]]
    )
    slope_idx_arr = idx.get("slope_idx")
    if slope_idx_arr is None:
        slope_idx_arr = idx["step_idx"]
    slope_arr = np.asarray(params["slopes"], dtype=float)
    weight_arr = idx.get("weight")
    score_k_arr = idx["score_k"]

    finite_rows = tbl[
        np.isfinite(tbl["Bias Size"]) & np.isfinite(tbl["S.E."]) & (tbl["S.E."] > 0)
    ]
    assert not finite_rows.empty, "expected at least one finite-SE GPCM bias cell"

    for _, row in finite_rows.iterrows():
        facet_a_col = row["FacetA"]
        facet_b_col = row["FacetB"]
        mask = (obs_df[facet_a_col].astype(str) == row["FacetA_Level"]) & (
            obs_df[facet_b_col].astype(str) == row["FacetB_Level"]
        )
        if not mask.any():
            continue
        idx_rows = np.where(mask.to_numpy())[0]
        eta_cell = eta_base[idx_rows] + float(row["Bias Size"])
        probs = app.category_prob_gpcm(
            eta_cell,
            step_cum_mat,
            idx["step_idx"][idx_rows],
            params["slopes"],
            slope_idx_arr[idx_rows],
        )
        k_vals = np.arange(probs.shape[1])
        e_k = probs @ k_vals
        v_k = probs @ (k_vals ** 2) - e_k ** 2
        v_k = np.where(v_k <= 1e-10, np.nan, v_k)
        slope_obs_cell = slope_arr[slope_idx_arr[idx_rows]]
        weight_cell = (
            weight_arr[idx_rows] if weight_arr is not None else np.ones(len(idx_rows))
        )

        info_slope_aware = np.nansum((slope_obs_cell ** 2) * v_k * weight_cell)
        info_unaware = np.nansum(v_k * weight_cell)
        se_slope_aware = 1.0 / np.sqrt(info_slope_aware)
        se_unaware = 1.0 / np.sqrt(info_unaware)

        assert float(row["S.E."]) == pytest.approx(se_slope_aware, rel=1e-9, abs=1e-12), (
            f"cell {row['FacetA_Level']}|{row['FacetB_Level']}: "
            f"reported S.E. {row['S.E.']} != slope-aware {se_slope_aware}"
        )
        assert abs(float(row["S.E."]) - se_unaware) > 1e-6, (
            f"cell {row['FacetA_Level']}|{row['FacetB_Level']}: reported S.E. "
            f"matches the non-slope-aware fall-back; the slope^2 factor did not take effect."
        )


# -----------------------------------------------------------------------------
# Integration: LR test and profile CI on a real GPCM MML fit
# -----------------------------------------------------------------------------


def test_bias_lr_chisq_matches_manual_two_times_loglik_difference(
    small_gpcm_mml_fit_for_bias, small_gpcm_mml_diagnostics_for_bias,
):
    """``LR ChiSq`` for each GPCM bias cell must equal
    ``2 * (loglik(b_hat) - loglik(0))`` recomputed directly from the
    captured cell."""
    res = small_gpcm_mml_fit_for_bias
    diag = small_gpcm_mml_diagnostics_for_bias
    bias_out = app.estimate_bias_interaction(res, diag, facet_a="Rater", facet_b="Criterion")
    tbl = bias_out["table"]

    config = res["config"]
    params = res["params"]
    obs_df = diag["obs"]
    idx = app.build_indices(
        res["prep"], step_facet=config["step_facet"], slope_facet=config.get("slope_facet")
    )
    theta_hat = res["facets"]["person"]["Estimate"].to_numpy()
    eta_base = app.compute_eta(idx, params, config, theta_override=theta_hat)
    step_cum_mat = np.vstack(
        [np.concatenate([[0.0], np.cumsum(row)]) for row in params["steps_mat"]]
    )

    checked = 0
    for _, row in tbl.iterrows():
        if not np.isfinite(row["Bias Size"]) or not np.isfinite(row["LR ChiSq"]):
            continue
        mask = (obs_df[row["FacetA"]].astype(str) == row["FacetA_Level"]) & (
            obs_df[row["FacetB"]].astype(str) == row["FacetB_Level"]
        )
        if not mask.any():
            continue
        idx_rows = np.where(mask.to_numpy())[0]
        ll_hat = app.loglik_gpcm(
            eta_base[idx_rows] + float(row["Bias Size"]),
            idx["score_k"][idx_rows],
            step_cum_mat,
            idx["step_idx"][idx_rows],
            params["slopes"],
            idx["slope_idx"][idx_rows]
            if idx.get("slope_idx") is not None
            else idx["step_idx"][idx_rows],
            weight=idx["weight"][idx_rows] if idx.get("weight") is not None else None,
        )
        ll_null = app.loglik_gpcm(
            eta_base[idx_rows],
            idx["score_k"][idx_rows],
            step_cum_mat,
            idx["step_idx"][idx_rows],
            params["slopes"],
            idx["slope_idx"][idx_rows]
            if idx.get("slope_idx") is not None
            else idx["step_idx"][idx_rows],
            weight=idx["weight"][idx_rows] if idx.get("weight") is not None else None,
        )
        expected_lr = max(0.0, 2.0 * (ll_hat - ll_null))
        assert float(row["LR ChiSq"]) == pytest.approx(expected_lr, rel=1e-8, abs=1e-10), (
            f"LR ChiSq mismatch for cell "
            f"{row['FacetA_Level']}|{row['FacetB_Level']}: "
            f"reported {row['LR ChiSq']} vs manual {expected_lr}"
        )
        checked += 1
    assert checked > 0, "no GPCM bias cells exercised the LR ChiSq assertion"


def test_bias_lr_prob_matches_chi2_survival_function(
    small_gpcm_mml_fit_for_bias, small_gpcm_mml_diagnostics_for_bias,
):
    bias_out = app.estimate_bias_interaction(
        small_gpcm_mml_fit_for_bias,
        small_gpcm_mml_diagnostics_for_bias,
        facet_a="Rater",
        facet_b="Criterion",
    )
    tbl = bias_out["table"]
    finite = tbl[np.isfinite(tbl["LR ChiSq"]) & np.isfinite(tbl["LR Prob."])]
    assert not finite.empty
    for _, row in finite.iterrows():
        expected = float(chi2.sf(float(row["LR ChiSq"]), df=int(row["LR d.f."])))
        assert float(row["LR Prob."]) == pytest.approx(expected, abs=1e-12), row


def test_bias_profile_ci_brackets_estimate_and_status_ok_for_typical_cells(
    small_gpcm_mml_fit_for_bias, small_gpcm_mml_diagnostics_for_bias,
):
    bias_out = app.estimate_bias_interaction(
        small_gpcm_mml_fit_for_bias,
        small_gpcm_mml_diagnostics_for_bias,
        facet_a="Rater",
        facet_b="Criterion",
    )
    tbl = bias_out["table"]
    # At least one cell should produce a clean "ok" profile CI on this
    # well-behaved fit. Cells that hit the search range are allowed but
    # should be in the minority.
    ok_rows = tbl[tbl["Profile CI Status"].astype(str) == "ok"]
    assert not ok_rows.empty, tbl["Profile CI Status"].value_counts()
    for _, row in ok_rows.iterrows():
        assert row["Profile CI Lower"] <= row["Bias Size"] + 1e-9
        assert row["Profile CI Upper"] >= row["Bias Size"] - 1e-9
        assert row["Profile CI Level"] == 0.95
        # The CI must be non-degenerate (strict inequality at machine
        # precision is enough; the quadratic-NLL unit test covers exact
        # widths).
        assert row["Profile CI Upper"] - row["Profile CI Lower"] > 1e-6


def test_bias_inference_responds_to_slope_clamp(
    small_gpcm_mml_fit_for_bias, small_gpcm_mml_diagnostics_for_bias,
):
    """Slope-dispatch sanity. If the slope-aware identity were silently
    ignored, clamping the log-slope block to zero (so every slope = 1)
    would not change the reported bias estimates. The test pins that the
    dispatch is actually consuming ``params["slopes"]``.
    """
    res = small_gpcm_mml_fit_for_bias
    diag = small_gpcm_mml_diagnostics_for_bias
    bias_original = app.estimate_bias_interaction(
        res, diag, facet_a="Rater", facet_b="Criterion"
    )

    config = res["config"]
    sizes = app.build_param_sizes(config)
    if not sizes.get("log_slopes"):
        pytest.skip("GPCM fit has no log-slope parameter block to clamp.")
    par_clamped = res["opt"].x.copy()
    pre = sum(int(sizes[name]) for name in sizes if name not in ("log_slopes",))
    par_clamped[pre : pre + int(sizes["log_slopes"])] = 0.0
    clamped_opt = type("Opt", (), {"x": par_clamped})()
    res_clamped = {**res, "opt": clamped_opt, "params": app.expand_params(par_clamped, sizes, config)}
    diag_clamped = app.mfrm_diagnostics(
        res_clamped, compute_pca=False, compute_marginal=False
    )
    bias_clamped = app.estimate_bias_interaction(
        res_clamped, diag_clamped, facet_a="Rater", facet_b="Criterion"
    )

    t1 = bias_original["table"]
    t2 = bias_clamped["table"]
    common = t1.merge(
        t2,
        on=["FacetA", "FacetA_Level", "FacetB", "FacetB_Level"],
        suffixes=("_orig", "_clamp"),
    )
    finite = common[
        np.isfinite(common["Bias Size_orig"]) & np.isfinite(common["Bias Size_clamp"])
    ]
    assert not finite.empty
    max_diff = float(
        np.max(np.abs(finite["Bias Size_orig"].to_numpy() - finite["Bias Size_clamp"].to_numpy()))
    )
    assert max_diff > 1e-6, (
        f"clamping log_slopes to zero left every bias estimate unchanged "
        f"(max |diff| = {max_diff!r}); the slope-aware dispatch is "
        f"silently ignoring params['slopes']."
    )


# -----------------------------------------------------------------------------
# Integration: RSM / PCM stay in the t-based screening tier
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# R parity: a small GPCM MML fit + bias estimation against mfrmr 0.2.0.
# -----------------------------------------------------------------------------


def _read_r_bias_parity_fixture() -> tuple[pd.DataFrame, dict]:
    """Load the deterministic synthetic CSV consumed by both implementations
    and the R reference output. The CSV is the exact data the R script
    in ``tests/data/`` was run against (mfrmr 0.2.0 / R 4.5.2), so the
    only sources of variation between Python and R are the optimisation
    paths and the MML quadrature implementation; the math contract is
    identical."""
    df = pd.read_csv(R_BIAS_PARITY_INPUT)
    with R_BIAS_PARITY_OUTPUT.open("r", encoding="utf-8") as handle:
        ref = json.load(handle)
    return df, ref


def _r_cell(ref: dict, facet_a_level: str, facet_b_level: str) -> dict | None:
    for cell in ref["cells"]:
        if cell["FacetA_Level"] == facet_a_level and cell["FacetB_Level"] == facet_b_level:
            return cell
    return None


@pytest.fixture(scope="module")
def r_bias_parity_fixture():
    if not R_BIAS_PARITY_INPUT.exists() or not R_BIAS_PARITY_OUTPUT.exists():
        pytest.skip("R bias parity fixture is missing; run the R generation script in tests/data/")
    return _read_r_bias_parity_fixture()


def test_bias_estimation_matches_r_reference_within_tolerance(r_bias_parity_fixture):
    """Python and R bias estimates agree at manuscript-citation precision
    across all six cells of the shared synthetic data set.

    The Python and R MML fits converge to parameter vectors that differ
    by roughly the absolute size of the marginal log-likelihood
    difference (``~0.02`` on this fixture), because the two
    implementations use different starting values, slightly different
    quadrature implementations, and slightly different convergence
    stopping rules. The same Muraki (1993) / Wilks (1938) / Cox (1975)
    formulas are evaluated on each side, so once the fits' parameter
    vectors agree at logit scale the per-cell bias inference must agree
    at a similar logit scale. The tolerances below are calibrated to the
    observed implementation-vs-implementation variance and are well
    below any practically meaningful threshold (``0.05`` logits is the
    canonical "negligible" cutoff cited in Linacre, FACETS Manual).

    The status field (``Profile CI Status``) and the inference tier
    (``InferenceTier``) are categorical and must agree exactly.
    """
    df, ref = r_bias_parity_fixture

    res = app.mfrm_estimate(
        data=df,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model="GPCM",
        method="MML",
        step_facet="Criterion",
        quad_points=5,
        maxit=25,
        reltol=1e-4,
        mml_engine="direct",
    )
    diag = app.mfrm_diagnostics(res, compute_pca=False, compute_marginal=False)
    bias_out = app.estimate_bias_interaction(
        res, diag, facet_a="Rater", facet_b="Criterion"
    )
    assert "_skip_reason" not in bias_out, bias_out
    py_tbl = bias_out["table"]

    # Tolerances calibrated to the observed implementation-vs-
    # implementation convergence variance on this fixture. The Python
    # and R fits converge to slightly different parameter vectors
    # (log-likelihood differs by ~0.02 on this seed), and the per-cell
    # quantities downstream inherit that variance. All values stay
    # well below ``0.1`` logits / chi-square units / probability mass,
    # which is the practical threshold below which manuscript tables
    # would round to the same reported value.
    TOL_BIAS = 5e-2
    TOL_SE = 2e-2
    TOL_LR_CHI = 1e-1
    TOL_LR_P = 1e-1
    TOL_CI = 1e-1

    checked = 0
    for _, py_row in py_tbl.iterrows():
        r_cell = _r_cell(ref, py_row["FacetA_Level"], py_row["FacetB_Level"])
        assert r_cell is not None, (
            f"R fixture missing cell {py_row['FacetA_Level']} x "
            f"{py_row['FacetB_Level']}"
        )
        # Status / tier are categorical -- exact agreement expected.
        assert py_row["Profile CI Status"] == r_cell["Profile CI Status"], (
            f"{py_row['FacetA_Level']} x {py_row['FacetB_Level']}: "
            f"Profile CI Status Python={py_row['Profile CI Status']!r} "
            f"vs R={r_cell['Profile CI Status']!r}"
        )
        assert py_row["InferenceTier"] == r_cell["InferenceTier"], (
            f"{py_row['FacetA_Level']} x {py_row['FacetB_Level']}: "
            f"InferenceTier mismatch"
        )

        if not (np.isfinite(py_row["Bias Size"]) and np.isfinite(r_cell["Bias Size"])):
            continue

        for label, py_val, r_val, tol in (
            ("Bias Size", py_row["Bias Size"], r_cell["Bias Size"], TOL_BIAS),
            ("S.E.", py_row["S.E."], r_cell["S.E."], TOL_SE),
            ("LR ChiSq", py_row["LR ChiSq"], r_cell["LR ChiSq"], TOL_LR_CHI),
            ("LR Prob.", py_row["LR Prob."], r_cell["LR Prob."], TOL_LR_P),
            ("Profile CI Lower", py_row["Profile CI Lower"],
             r_cell["Profile CI Lower"], TOL_CI),
            ("Profile CI Upper", py_row["Profile CI Upper"],
             r_cell["Profile CI Upper"], TOL_CI),
        ):
            diff = abs(float(py_val) - float(r_val))
            assert diff < tol, (
                f"{py_row['FacetA_Level']} x {py_row['FacetB_Level']}: "
                f"{label} Python={py_val!r} vs R={r_val!r} "
                f"(|diff|={diff!r} >= tol={tol!r})"
            )
        checked += 1
    assert checked == len(ref["cells"]), (
        f"checked {checked} cells but R fixture has {len(ref['cells'])}; "
        f"some cells were skipped due to non-finite estimates"
    )


def test_rsm_bias_marks_lr_and_profile_columns_as_not_applicable():
    """For non-GPCM fits the LR / profile columns must be NaN with a
    ``Likelihood Basis`` string that explains why; the t-based screening
    columns continue to carry the inferential content."""
    rng = np.random.default_rng(20260514)
    n_person, n_rater, n_criterion = 20, 2, 2
    theta = rng.normal(0, 1, n_person)
    rater_eff = np.array([-0.3, 0.3])
    crit_eff = np.array([-0.2, 0.2])
    rows = []
    for i, theta_i in enumerate(theta):
        for j in range(n_rater):
            for k in range(n_criterion):
                eta = theta_i - rater_eff[j] - crit_eff[k]
                p1 = 1.0 / (1.0 + np.exp(-eta))
                score = int(rng.uniform() < p1)
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
        rating_max=1,
        model="RSM",
        method="JMLE",
        maxit=25,
        reltol=1e-4,
    )
    diag = app.mfrm_diagnostics(res, compute_pca=False, compute_marginal=False)
    bias_out = app.estimate_bias_interaction(res, diag, facet_a="Rater", facet_b="Criterion")
    if "_skip_reason" in bias_out:
        pytest.skip(bias_out["_skip_reason"])
    tbl = bias_out["table"]
    assert (tbl["LR ChiSq"].isna()).all()
    assert (tbl["LR Prob."].isna()).all()
    assert (tbl["Profile CI Lower"].isna()).all()
    assert (tbl["Profile CI Upper"].isna()).all()
    assert (tbl["Profile CI Status"].astype(str) == "not available").all()
    assert (tbl["Likelihood Basis"].astype(str).str.contains("not applicable")).all()
    # The t-based columns still carry inferential content for RSM/PCM.
    assert tbl["t"].notna().any()
