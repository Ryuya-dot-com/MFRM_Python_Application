"""Mathematical-contract tests for the slope-aware GPCM fair-average kernel.

These tests pin the kernel of ``expected_score_from_eta`` and the wiring in
``calc_facets_report_tbls`` that selects the per-row slope and per-row
cumulative thresholds for the Linacre fair-average.

Kernel contract:

* At ``slope == 1`` the kernel reduces byte-for-byte to the PCM/RSM Linacre
  expected score (Masters 1982; Andrich 1978).
* A worked example matches a hand-computed value to better than ``1e-6``.
* The numerical derivative of ``E[X|eta]`` with respect to ``eta`` equals
  ``a * Var[X|eta]`` (Muraki 1993, Eq. 16).
* Invariance under the GPCM identification rescaling
  ``(eta, delta_cum, a) -> (eta / c, delta_cum / c, a * c)``.
* Degenerate slopes (zero, negative, non-finite) and non-finite ``eta``
  both return ``np.nan`` instead of silently falling back to a number.
* Boundary fixtures: singleton category, binary category, extreme slopes.
* 27 cases generated directly from the mfrmr R kernel agree to ``1e-10``.

Integration contract for ``calc_facets_report_tbls``:

* Non slope-facet rows are evaluated at slope = 1 under the GPCM model
  (the geometric-mean-one identification), not at the arithmetic mean of
  the estimated slopes.
* Slope-facet rows are evaluated at that row's own discrimination and
  that row's own step cumulative thresholds.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, softmax

import streamlit_app as app


R_FIXTURE_PATH = Path(__file__).resolve().parent / "data" / "r_kernel_fixture.json"


def _pcm_reference_expected_score(eta: float, step_cum, rating_min: int) -> float:
    """Independent textbook PCM expected score (no slope multiplication)."""
    step_cum = np.asarray(step_cum, dtype=float)
    k_vals = np.arange(len(step_cum), dtype=float)
    log_num = k_vals * eta - step_cum
    probs = softmax(log_num)
    return float(rating_min + np.sum(probs * k_vals))


def test_kernel_at_slope_one_matches_pcm_reference():
    cases = [
        (-1.5, [0.0, -0.5, 0.0, 0.5], 0),
        (-0.3, [0.0, -0.5, 0.0, 0.5], 0),
        (0.0, [0.0, -0.5, 0.0, 0.5], 0),
        (0.5, [0.0, -0.5, 0.0, 0.5], 1),
        (1.2, [0.0, -1.0, 0.5], 1),
        (-0.7, [0.0, -1.0, 0.5], 0),
        (2.0, [0.0, 0.2], 0),
    ]
    for eta, step_cum, r_min in cases:
        ref = _pcm_reference_expected_score(eta, step_cum, r_min)
        out = app.expected_score_from_eta(eta, step_cum, r_min, slope=1.0)
        assert abs(out - ref) < 1e-15, (
            f"slope=1 path must match PCM reference at eta={eta}, "
            f"step_cum={step_cum}: ref={ref!r} vs out={out!r}"
        )


def test_kernel_worked_example_matches_hand_calculation():
    # Worked example: K = 4 ordered categories, per-step delta = (-0.5, 0, 0.5)
    # so the cumulative threshold vector is (0, -0.5, -0.5, 0). At eta = 0.3
    # and slope = 1.2 the category probabilities are
    #     (0.097089, 0.253566, 0.363447, 0.285898)
    # and the expected score (rating_min = 0) is 1.8381511.
    fa = app.expected_score_from_eta(
        eta=0.3, step_cum=[0.0, -0.5, -0.5, 0.0], rating_min=0, slope=1.2
    )
    assert abs(fa - 1.8381511) < 1e-6, fa


def test_numerical_derivative_matches_slope_times_variance():
    # Muraki (1993, Eq. 16): d E[X|eta] / d eta = a * Var[X|eta]. For the
    # worked-example fixture above, Var[K] = 0.901626 and slope = 1.2 give
    # the analytic value 1.0819512.
    delta_cum = [0.0, -0.5, -0.5, 0.0]
    h = 1e-6
    fa_plus = app.expected_score_from_eta(0.3 + h, delta_cum, 0, slope=1.2)
    fa_minus = app.expected_score_from_eta(0.3 - h, delta_cum, 0, slope=1.2)
    deriv_num = (fa_plus - fa_minus) / (2 * h)
    assert abs(deriv_num - 1.0819512) < 1e-4, deriv_num


def test_kernel_invariant_under_slope_rescaling():
    # The GPCM identification permits rescaling (eta, delta_cum, slope) by
    # any positive constant c via (eta / c, delta_cum / c, slope * c) without
    # changing the probability vector. The fair-average must therefore be
    # invariant.
    delta_cum = np.array([0.0, -0.5, -0.5, 0.0], dtype=float)
    c = 2.0
    fa_orig = app.expected_score_from_eta(0.3, delta_cum, 0, slope=1.2)
    fa_resc = app.expected_score_from_eta(0.3 / c, delta_cum / c, 0, slope=1.2 * c)
    assert abs(fa_orig - fa_resc) < 1e-12, (fa_orig, fa_resc)


def test_degenerate_slope_returns_nan_not_silent_fallback():
    step_cum = [0.0, -0.5, 0.5]
    fa_one = app.expected_score_from_eta(0.3, step_cum, 0, slope=1.0)
    fa_zero = app.expected_score_from_eta(0.3, step_cum, 0, slope=0.0)
    fa_neg = app.expected_score_from_eta(0.3, step_cum, 0, slope=-0.5)
    fa_nan = app.expected_score_from_eta(0.3, step_cum, 0, slope=float("nan"))
    assert np.isfinite(fa_one)
    assert np.isnan(fa_zero), fa_zero
    assert np.isnan(fa_neg), fa_neg
    assert np.isnan(fa_nan), fa_nan


def test_kernel_actually_consumes_slope_argument():
    # Sanity dispatch: if the kernel silently ignored its slope argument
    # (for example by clamping to 1.0), two distinct positive slopes would
    # produce identical expected scores. Pinning a measurable response
    # protects against such silent regressions in future refactors.
    eta = 0.3
    step_cum = [0.0, -0.5, -0.5, 0.0]
    fa_low = app.expected_score_from_eta(eta, step_cum, 0, slope=0.5)
    fa_high = app.expected_score_from_eta(eta, step_cum, 0, slope=1.5)
    assert abs(fa_low - fa_high) > 1e-3, (fa_low, fa_high)


# -----------------------------------------------------------------------------
# Boundary fixtures: singleton category, binary, non-finite eta, extreme slope.
# -----------------------------------------------------------------------------


def test_kernel_single_category_returns_rating_min():
    """With K = 1 the only legal score is rating_min itself."""
    for eta in [-1.0, 0.0, 1.0]:
        for slope in [0.5, 1.0, 2.0]:
            for r_min in [0, 1, 5]:
                fa = app.expected_score_from_eta(eta, [0.0], r_min, slope=slope)
                assert abs(fa - r_min) < 1e-15, (eta, slope, r_min, fa)


def test_kernel_binary_matches_logistic_sigmoid():
    """For K = 2 the expected score reduces to ``rating_min + sigmoid(slope * (eta - sc_1))``."""
    for eta in [-0.5, 0.0, 0.3]:
        for slope in [0.5, 1.0, 1.7]:
            for sc_1 in [-0.5, 0.0, 0.5]:
                fa = app.expected_score_from_eta(eta, [0.0, sc_1], rating_min=0, slope=slope)
                expected = float(expit(slope * (eta - sc_1)))
                assert abs(fa - expected) < 1e-12, (eta, slope, sc_1, fa, expected)


def test_kernel_returns_nan_for_non_finite_eta():
    step_cum = [0.0, -0.5, 0.5]
    for eta in [float("inf"), float("-inf"), float("nan")]:
        fa = app.expected_score_from_eta(eta, step_cum, rating_min=0, slope=1.0)
        assert np.isnan(fa), (eta, fa)


def test_kernel_handles_extreme_slope_without_overflow():
    """Tiny and large positive slopes both keep the kernel numerically sane."""
    step_cum = [0.0, -0.5, 0.5]

    # Tiny slope flattens the distribution toward a uniform over K = 3
    # categories, whose mean is exactly (K - 1) / 2 = 1.
    fa_flat = app.expected_score_from_eta(0.0, step_cum, rating_min=0, slope=1e-6)
    assert abs(fa_flat - 1.0) < 1e-3, fa_flat

    # Large slope with positive eta concentrates mass on the top category.
    fa_sharp_top = app.expected_score_from_eta(2.0, step_cum, rating_min=0, slope=10.0)
    assert np.isfinite(fa_sharp_top), fa_sharp_top
    assert fa_sharp_top > 1.95, fa_sharp_top

    # Large slope with negative eta concentrates mass on the bottom category.
    fa_sharp_bot = app.expected_score_from_eta(-2.0, step_cum, rating_min=0, slope=10.0)
    assert np.isfinite(fa_sharp_bot), fa_sharp_bot
    assert fa_sharp_bot < 0.05, fa_sharp_bot


# -----------------------------------------------------------------------------
# Parity against the mfrmr R kernel: 27 cases at machine tolerance.
# -----------------------------------------------------------------------------


def test_kernel_matches_r_reference_fixture():
    """Each Python value must agree with the mfrmr R kernel within 1e-10.

    The fixture was generated by ``mfrmr:::expected_score_from_eta_gpcm()``
    (mfrmr 0.2.0, R 4.5.2) over 3 etas x 3 slopes x 3 step-cumulative
    profiles. The agreement at machine tolerance is the strongest evidence
    we can offer that the Python kernel and the R reference encode the
    same Muraki (1992) Eqs. 2-3, 10 form.
    """
    with R_FIXTURE_PATH.open("r", encoding="utf-8") as handle:
        cases = json.load(handle)
    assert len(cases) == 27, f"expected 27 R fixture cases, got {len(cases)}"
    for case in cases:
        py = app.expected_score_from_eta(
            case["eta"],
            case["step_cum"],
            case["rating_min"],
            slope=case["slope"],
        )
        r_value = case["fair_average"]
        diff = abs(py - r_value)
        assert diff < 1e-10, (
            f"profile={case['profile']!r} eta={case['eta']} slope={case['slope']}: "
            f"Python={py!r} R={r_value!r} diff={diff!r}"
        )


# -----------------------------------------------------------------------------
# Integration contract for calc_facets_report_tbls: slope_list / step_cum_list
# wiring under GPCM.
# -----------------------------------------------------------------------------


def _gpcm_two_criterion_fixture():
    """Two-Criterion GPCM fit fixture with slopes whose geometric mean is one
    and whose arithmetic mean is 1.25, so the two contracts disagree and the
    test can distinguish them."""
    rating_min, rating_max = 0, 2
    step_mat = np.array([[-0.5, 0.5], [-0.3, 0.3]], dtype=float)
    slopes = np.array([0.5, 2.0], dtype=float)
    res = {
        "prep": {
            "rating_min": rating_min,
            "rating_max": rating_max,
            "levels": {
                "Person": ["P1", "P2"],
                "Rater": ["R1"],
                "Criterion": ["C1", "C2"],
            },
        },
        "config": {
            "method": "MML",
            "model": "GPCM",
            "step_facet": "Criterion",
            "facet_names": ["Rater", "Criterion"],
            "facet_signs": {"Rater": -1, "Criterion": -1},
            "theta_spec": None,
            "facet_specs": {},
        },
        "params": {
            "theta": np.array([0.0, 0.5], dtype=float),
            "steps_mat": step_mat,
            "slopes": slopes,
            "facets": {
                "Rater": np.array([0.0], dtype=float),
                "Criterion": np.array([-0.2, 0.2], dtype=float),
            },
        },
        "facets": {
            "person": pd.DataFrame({"Person": ["P1", "P2"], "Estimate": [0.0, 0.5]}),
        },
    }
    diagnostics = {
        "obs": pd.DataFrame(
            {
                "Person": ["P1", "P1", "P2", "P2"],
                "Rater": ["R1", "R1", "R1", "R1"],
                "Criterion": ["C1", "C2", "C1", "C2"],
                "Observed": [1, 1, 2, 1],
                "Weight": [1.0, 1.0, 1.0, 1.0],
            }
        ),
        "measures": pd.DataFrame(
            {
                "Facet": ["Person", "Person", "Rater", "Criterion", "Criterion"],
                "Level": ["P1", "P2", "R1", "C1", "C2"],
                "Estimate": [0.0, 0.5, 0.0, -0.2, 0.2],
                "SE": [0.1] * 5,
                "Infit": [1.0] * 5,
                "Outfit": [1.0] * 5,
                "InfitZSTD": [0.0] * 5,
                "OutfitZSTD": [0.0] * 5,
                "PTMEA": [0.5] * 5,
            }
        ),
    }
    return res, diagnostics


def test_facets_report_uses_unit_slope_on_non_slope_facet_rows():
    """Person rows under GPCM must use slope = 1, not the arithmetic mean.

    With slopes (0.5, 2.0) the geometric mean is 1 and the arithmetic mean
    is 1.25. Under the sum-to-zero log-slope identification, slope = 1 is
    the neutral discrimination on non slope-facet rows; the historic
    Python fallback used the arithmetic mean and so reported a systematically
    over-discriminated fair average there.
    """
    res, diagnostics = _gpcm_two_criterion_fixture()
    rating_min = res["prep"]["rating_min"]

    # Hand-derived eta for Person P2 under the function's own conventions:
    #   theta_mean = mean([0.0, 0.5]) = 0.25 (unused for Person row)
    #   facet means cancel because they are all zero in this fixture
    #   eta_m = Estimate_P2 + 0 = 0.5
    # step_cum on a non step-facet row collapses to column-mean cumsum:
    #   step_cum_mean = (0, cumsum(colMeans(step_mat))) = (0, -0.4, 0)
    step_cum_mean = [0.0, -0.4, 0.0]
    expected_at_slope_one = app.expected_score_from_eta(
        0.5, step_cum_mean, rating_min, slope=1.0
    )
    arithmetic_mean_slope = float(np.mean(res["params"]["slopes"]))
    fair_under_old_contract = app.expected_score_from_eta(
        0.5, step_cum_mean, rating_min, slope=arithmetic_mean_slope
    )
    # Premise of the test: the two contracts disagree on this fixture.
    assert abs(expected_at_slope_one - fair_under_old_contract) > 1e-3

    tbls = app.calc_facets_report_tbls(res, diagnostics)
    person_tbl = tbls["Person"]
    p2_fair_m = float(person_tbl.loc[person_tbl["Level"] == "P2", "FairM"].iloc[0])

    assert abs(p2_fair_m - expected_at_slope_one) < 1e-12, (
        f"Person row FairM must use slope = 1: expected {expected_at_slope_one}, "
        f"got {p2_fair_m}."
    )
    assert abs(p2_fair_m - fair_under_old_contract) > 1e-3, (
        "Person row FairM still equals the arithmetic-mean-slope path; "
        "the bug fix did not take effect."
    )


def test_facets_report_uses_per_level_slope_on_slope_facet_rows():
    """Each Criterion row must use its own slope and its own step thresholds."""
    res, diagnostics = _gpcm_two_criterion_fixture()
    rating_min = res["prep"]["rating_min"]
    slopes = res["params"]["slopes"]

    tbls = app.calc_facets_report_tbls(res, diagnostics)
    criterion_tbl = tbls["Criterion"]

    # Per-level expectations: eta on a slope-facet row is
    #   eta_m = theta_mean + sum(sign[F'] * mean[F']) over F' != Criterion
    #             + sign[Criterion] * Estimate_Criterion
    # All non-Criterion facet means are zero, so:
    #   eta_m(C1) = 0.25 + 0 + (-1)*(-0.2) = 0.45
    #   eta_m(C2) = 0.25 + 0 + (-1)*( 0.2) = 0.05
    # step_cum follows the row's own level:
    #   step_cum(C1) = (0, cumsum([-0.5, 0.5])) = (0, -0.5, 0)
    #   step_cum(C2) = (0, cumsum([-0.3, 0.3])) = (0, -0.3, 0)
    expected_c1 = app.expected_score_from_eta(
        0.45, [0.0, -0.5, 0.0], rating_min, slope=float(slopes[0])
    )
    expected_c2 = app.expected_score_from_eta(
        0.05, [0.0, -0.3, 0.0], rating_min, slope=float(slopes[1])
    )

    c1_fair_m = float(criterion_tbl.loc[criterion_tbl["Level"] == "C1", "FairM"].iloc[0])
    c2_fair_m = float(criterion_tbl.loc[criterion_tbl["Level"] == "C2", "FairM"].iloc[0])

    assert abs(c1_fair_m - expected_c1) < 1e-12, (c1_fair_m, expected_c1)
    assert abs(c2_fair_m - expected_c2) < 1e-12, (c2_fair_m, expected_c2)

    # And the two rows must in fact differ -- otherwise the per-row slope is
    # collapsing to a shared value somewhere upstream.
    assert abs(c1_fair_m - c2_fair_m) > 1e-3, (c1_fair_m, c2_fair_m)
