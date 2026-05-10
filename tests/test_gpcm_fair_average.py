"""Mathematical-contract tests for the slope-aware GPCM fair-average kernel.

These tests pin the kernel of ``expected_score_from_eta``:

1. At ``slope == 1`` the kernel reduces byte-for-byte to the PCM/RSM
   Linacre expected score (Masters 1982; Andrich 1978).
2. A worked example matches a hand-computed value to better than
   ``1e-6`` absolute (cross-checks with the mfrmr R reference fixture).
3. The numerical derivative of ``E[X|eta]`` with respect to ``eta``
   equals ``a * Var[X|eta]`` (Muraki 1993, Eq. 16).
4. The kernel is invariant under the GPCM identification rescaling
   ``(eta, delta_cum, a) -> (eta / c, delta_cum / c, a * c)``.
5. Degenerate slopes (zero, negative, non-finite) return ``np.nan``
   instead of silently falling back to slope = 1.
6. The kernel actually consumes its ``slope`` argument: distinct slopes
   produce distinct expected scores.
"""

from __future__ import annotations

import numpy as np
from scipy.special import softmax

import streamlit_app as app


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
