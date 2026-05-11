"""Tests for the MML observed-information covariance and the structural
delta-method standard error / confidence interval of the slope-aware GPCM
fair-average.

Unit contracts:
    * ``_invert_information_matrix`` recovers the inverse of a positive-definite
      input, flags near-singular spectra as regularized, and rejects
      non-symmetric or non-finite input.
    * ``_finite_difference_gradient`` agrees with the analytical gradient of a
      simple quadratic at machine tolerance and falls back to one-sided
      differences when the symmetric perturbation is non-finite.
    * ``_build_param_slices`` returns disjoint, contiguous slices that cover
      the full parameter vector exactly.

Integration contracts (require a small GPCM MML fit):
    * ``compute_mml_parameter_covariance`` returns ``status = "not_applicable"``
      for a JMLE fit (no marginal likelihood) and finite, symmetric covariance
      with ``status in {"ok", "regularized"}`` for a small GPCM MML fit.
    * ``add_gpcm_fair_average_delta_se`` annotates each facet's table with the
      delta-method SE / CI columns; Person rows are marked ``"not available"``
      because MML EAPs are not part of the structural Hessian.
    * ``fair_average_table(fit_se=False)`` is identical to ``calc_facets_report_tbls``;
      ``fair_average_table(fit_se=True)`` adds the SE columns.

The integration tests share one MML fit fixture to keep total runtime
reasonable -- the fit itself plus the Hessian and per-row delta-method
gradient runs comfortably under a minute on a typical laptop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Unit: _invert_information_matrix
# -----------------------------------------------------------------------------


def test_invert_information_matrix_recovers_inverse_of_psd_input():
    """For an identity Hessian the covariance is the identity itself."""
    info = np.eye(4)
    cov, regularized, rank = app._invert_information_matrix(info)
    assert cov is not None
    assert np.allclose(cov, np.eye(4), atol=1e-12)
    assert regularized is False
    assert rank == 4


def test_invert_information_matrix_recovers_inverse_of_generic_psd():
    """For a generic positive-definite matrix the inverse round-trips at 1e-10."""
    rng = np.random.default_rng(0)
    A = rng.standard_normal((5, 5))
    info = A @ A.T + np.eye(5)  # PSD with positive diagonal
    cov, regularized, rank = app._invert_information_matrix(info)
    assert cov is not None
    assert regularized is False
    assert rank == 5
    assert np.allclose(info @ cov, np.eye(5), atol=1e-10)
    # Symmetric output
    assert np.allclose(cov, cov.T, atol=1e-14)


def test_invert_information_matrix_flags_near_singular_spectrum():
    """A rank-deficient information matrix triggers regularization."""
    # Build a matrix with one near-zero eigenvalue.
    eigvals = np.array([10.0, 5.0, 1.0, 1e-16])
    Q, _ = np.linalg.qr(np.random.default_rng(1).standard_normal((4, 4)))
    info = Q @ np.diag(eigvals) @ Q.T
    cov, regularized, rank = app._invert_information_matrix(info)
    assert cov is not None
    assert regularized is True
    assert rank == 3


def test_invert_information_matrix_rejects_non_finite():
    info = np.array([[1.0, np.nan], [np.nan, 1.0]])
    cov, regularized, rank = app._invert_information_matrix(info)
    assert cov is None
    assert regularized is False
    assert rank == 0


def test_invert_information_matrix_rejects_empty_or_misshaped():
    assert app._invert_information_matrix(None) == (None, False, 0)
    assert app._invert_information_matrix(np.array([])) == (None, False, 0)
    cov, regularized, rank = app._invert_information_matrix(np.zeros((2, 3)))
    assert cov is None and regularized is False and rank == 0


# -----------------------------------------------------------------------------
# Unit: _finite_difference_gradient
# -----------------------------------------------------------------------------


def test_finite_difference_gradient_matches_closed_form_quadratic():
    """For ``f(x) = x_0^2 + 2 * x_1 + 3 * x_2`` the central-difference
    gradient matches the analytic gradient ``(2*x_0, 2, 3)`` at 1e-7."""

    def f(p):
        return p[0] ** 2 + 2.0 * p[1] + 3.0 * p[2]

    par = np.array([1.5, -0.5, 0.7])
    grad = app._finite_difference_gradient(f, par)
    expected = np.array([2.0 * par[0], 2.0, 3.0])
    assert np.allclose(grad, expected, atol=1e-7)


def test_finite_difference_gradient_falls_back_to_one_sided_for_nan_branch():
    """If the symmetric perturbation hits a NaN branch, the routine falls
    back to a one-sided difference using the base evaluation."""

    def f(p):
        if p[0] < 0:
            return float("nan")
        return p[0] ** 2

    par = np.array([1e-12])
    grad = app._finite_difference_gradient(f, par)
    # Forward difference: (f(par + h) - f(par)) / h with f(par) = 1e-24 and
    # f(par + h) approx h^2 + 2*par*h. h = 1e-5 here, so grad approx 2 * par + h
    # which is dominated by h itself (1e-5). The point of the test is finiteness,
    # not the exact value -- assert finite output and a sane sign.
    assert np.isfinite(grad[0])
    assert grad[0] >= 0


def test_finite_difference_gradient_empty_par():
    grad = app._finite_difference_gradient(lambda p: 0.0, np.array([]))
    assert grad.size == 0


# -----------------------------------------------------------------------------
# Unit: _build_param_slices
# -----------------------------------------------------------------------------


def test_build_param_slices_covers_full_par_vector_disjointly():
    from collections import OrderedDict

    sizes = OrderedDict([("theta", 0), ("a", 3), ("b", 2), ("steps", 4)])
    slices = app._build_param_slices(sizes)
    par = np.arange(9)
    seen = np.zeros(9, dtype=bool)
    for name, sl in slices.items():
        block = par[sl]
        # Disjoint coverage: every covered index must be hit at most once.
        idx_range = range(sl.start, sl.stop)
        for k in idx_range:
            assert not seen[k], f"slice {name} re-covers index {k}"
            seen[k] = True
        assert block.size == sizes[name]
    # All indices that map to non-empty sections must be covered.
    assert seen.sum() == sum(int(v) for v in sizes.values())


# -----------------------------------------------------------------------------
# Integration: fit fixture + SE machinery on a small GPCM MML fit.
# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def small_gpcm_mml_fit():
    """Run a small GPCM MML fit once and reuse across the integration tests.

    The data is a 30-person x 2-rater x 3-criterion GPCM design with a
    rating scale of 0..2. ``quad_points = 5`` and ``maxit = 25`` keep the
    fit and the downstream Hessian / delta-method evaluations under a
    minute even on a modest laptop while still exercising every code
    path in the SE pipeline.
    """
    rng = np.random.default_rng(20260512)
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
def small_gpcm_mml_diagnostics(small_gpcm_mml_fit):
    return app.mfrm_diagnostics(small_gpcm_mml_fit, compute_pca=False, compute_marginal=False)


# -----------------------------------------------------------------------------
# Integration: compute_mml_parameter_covariance
# -----------------------------------------------------------------------------


def test_compute_mml_parameter_covariance_is_not_applicable_for_jmle():
    """Without a marginal likelihood the structural covariance is undefined."""
    res = {"config": {"method": "JMLE", "model": "RSM"}}
    out = app.compute_mml_parameter_covariance(res)
    assert out["status"] == "not_applicable"
    assert out["cov"] is None
    assert "MML" in out["detail"]


def test_compute_mml_parameter_covariance_returns_ok_for_small_gpcm_mml_fit(
    small_gpcm_mml_fit,
):
    cov_info = app.compute_mml_parameter_covariance(small_gpcm_mml_fit)
    assert cov_info["status"] in {"ok", "regularized"}, cov_info
    cov = cov_info["cov"]
    assert cov is not None
    assert np.all(np.isfinite(cov))
    # Symmetric.
    assert np.allclose(cov, cov.T, atol=1e-10)
    # Positive semi-definite (eigenvalues >= -tol after regularization).
    eigvals = np.linalg.eigvalsh((cov + cov.T) / 2)
    assert eigvals.min() >= -1e-10, eigvals.min()
    # Shape matches the optimisation parameter vector.
    par = app._get_opt_par(small_gpcm_mml_fit)
    assert cov.shape == (par.size, par.size)
    # Hessian round-trips: I_obs * Cov approx I, modulo any regularization.
    hess = cov_info["hessian"]
    assert hess is not None
    product = hess @ cov
    assert np.allclose(product, np.eye(par.size), atol=1e-6) or cov_info["regularized"]


# -----------------------------------------------------------------------------
# Integration: add_gpcm_fair_average_delta_se
# -----------------------------------------------------------------------------


_SE_COLUMNS = (
    "FairMSE",
    "FairM_CI_Lower",
    "FairM_CI_Upper",
    "FairM_CI_Level",
    "FairM_SE_Method",
    "FairM_SE_Status",
    "FairM_SE_Detail",
    "FairZSE",
    "FairZ_CI_Lower",
    "FairZ_CI_Upper",
    "FairZ_CI_Level",
    "FairZ_SE_Method",
    "FairZ_SE_Status",
    "FairZ_SE_Detail",
)


def test_delta_method_se_marks_person_rows_unavailable_and_others_finite(
    small_gpcm_mml_fit, small_gpcm_mml_diagnostics,
):
    raw_tbls = app.calc_facets_report_tbls(small_gpcm_mml_fit, small_gpcm_mml_diagnostics)
    annotated = app.add_gpcm_fair_average_delta_se(raw_tbls, small_gpcm_mml_fit)

    assert set(annotated.keys()) == set(raw_tbls.keys())

    for facet, tbl in annotated.items():
        for col in _SE_COLUMNS:
            assert col in tbl.columns, f"{facet} table missing column {col}"

    # Person rows must carry status "not available" and NaN SE.
    person_tbl = annotated["Person"]
    assert (person_tbl["FairM_SE_Status"] == "not available").all()
    assert person_tbl["FairMSE"].isna().all()
    assert person_tbl["FairZ_SE_Status"].notna().all()
    assert person_tbl["FairZSE"].isna().all()

    # Non-Person rows under GPCM should be either ok or regularized, with at
    # least some finite SEs.
    for facet, tbl in annotated.items():
        if facet == "Person":
            continue
        non_person_status = tbl["FairM_SE_Status"]
        assert non_person_status.isin({"ok", "regularized"}).any(), (
            f"{facet} expected at least one ok / regularized SE row, got {non_person_status.tolist()!r}"
        )
        finite_se = tbl.loc[tbl["FairM_SE_Status"].isin({"ok", "regularized"}), "FairMSE"]
        assert finite_se.notna().all() and (finite_se > 0).all()


def test_delta_method_ci_lies_within_rating_bounds(
    small_gpcm_mml_fit, small_gpcm_mml_diagnostics,
):
    raw_tbls = app.calc_facets_report_tbls(small_gpcm_mml_fit, small_gpcm_mml_diagnostics)
    annotated = app.add_gpcm_fair_average_delta_se(raw_tbls, small_gpcm_mml_fit, ci_level=0.95)
    prep = small_gpcm_mml_fit["prep"]
    rating_min = float(prep["rating_min"])
    rating_max = float(prep["rating_max"])

    for facet, tbl in annotated.items():
        if facet == "Person":
            continue
        finite_mask = tbl["FairM_SE_Status"].isin({"ok", "regularized"}) & tbl["FairMSE"].notna()
        if not finite_mask.any():
            continue
        sub = tbl.loc[finite_mask]
        # CI must be within [rating_min, rating_max] and bracket the estimate.
        assert (sub["FairM_CI_Lower"] >= rating_min - 1e-9).all()
        assert (sub["FairM_CI_Upper"] <= rating_max + 1e-9).all()
        assert (sub["FairM_CI_Lower"] <= sub["FairM"] + 1e-9).all()
        assert (sub["FairM_CI_Upper"] >= sub["FairM"] - 1e-9).all()


def test_delta_method_is_not_applicable_for_non_gpcm_fits():
    """An RSM-fit-shaped result must report ``status = "not_applicable"`` rather
    than emitting a number."""
    res = {
        "config": {"method": "MML", "model": "RSM"},
        "prep": {"rating_min": 0, "rating_max": 2},
        "opt": type("Opt", (), {"x": np.zeros(3)})(),
    }
    raw_tbls = {"Person": pd.DataFrame({"Level": ["P1"], "FairM": [1.0], "FairZ": [1.0]})}
    annotated = app.add_gpcm_fair_average_delta_se(raw_tbls, res)
    assert annotated["Person"]["FairM_SE_Status"].iloc[0] == "not_applicable"
    assert np.isnan(annotated["Person"]["FairMSE"].iloc[0])


# -----------------------------------------------------------------------------
# Integration: fair_average_table public wrapper
# -----------------------------------------------------------------------------


def test_fair_average_table_passes_through_when_fair_se_false(
    small_gpcm_mml_fit, small_gpcm_mml_diagnostics,
):
    """``fair_se=False`` returns ``calc_facets_report_tbls`` output verbatim."""
    plain = app.calc_facets_report_tbls(small_gpcm_mml_fit, small_gpcm_mml_diagnostics)
    wrapped = app.fair_average_table(
        small_gpcm_mml_fit, small_gpcm_mml_diagnostics, fair_se=False
    )
    assert set(plain.keys()) == set(wrapped.keys())
    for facet in plain:
        pd.testing.assert_frame_equal(plain[facet], wrapped[facet])


def test_fair_average_table_adds_se_columns_when_fair_se_true(
    small_gpcm_mml_fit, small_gpcm_mml_diagnostics,
):
    annotated = app.fair_average_table(
        small_gpcm_mml_fit, small_gpcm_mml_diagnostics, fair_se=True, ci_level=0.95
    )
    for facet, tbl in annotated.items():
        for col in _SE_COLUMNS:
            assert col in tbl.columns, f"{facet} missing {col} when fair_se=True"
        # CI level column reflects the request.
        assert (tbl["FairM_CI_Level"] == 0.95).all()
        assert (tbl["FairZ_CI_Level"] == 0.95).all()
