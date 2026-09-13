"""Numerical contracts for fixed-calibration Warm Person scoring."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mfrm_app.person_scoring import (
    evaluate_fixed_calibration_score,
    score_fixed_calibration_persons,
)


def _binary_rows(scores: list[list[int]]):
    n_person = len(scores)
    n_items = len(scores[0])
    return (
        np.repeat([f"P{index + 1}" for index in range(n_person)], n_items),
        np.asarray(scores, dtype=int).reshape(-1),
        np.zeros((n_person * n_items, 2), dtype=float),
        np.tile(np.array([[0.0, 1.0]]), (n_person * n_items, 1)),
    )


def test_binary_rasch_wle_matches_closed_form_at_both_extremes_and_interior():
    scores = [
        [0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 1, 1, 1, 1],
    ]
    person, observed, intercepts, slopes = _binary_rows(scores)
    result = score_fixed_calibration_persons(person, observed, intercepts, slopes)
    raw_scores = np.array([0.0, 2.0, 5.0])
    probabilities = (raw_scores + 0.5) / 6.0
    expected = np.log(probabilities / (1.0 - probabilities))
    expected_se = np.sqrt(1.0 / (5.0 * probabilities * (1.0 - probabilities)))

    assert np.allclose(result["Estimate"], expected, rtol=0.0, atol=1e-10)
    assert np.allclose(result["StandardError"], expected_se, rtol=0.0, atol=1e-10)
    assert result["ExtremeScoreDirection"].tolist() == ["all_minimum", "interior", "all_maximum"]
    assert result["Status"].eq("ok").all()
    assert result["AdjustedScoreResidual"].max() <= 1e-10
    assert result["EstimateRole"].eq("fixed_calibration_wle_not_jmle_mle").all()


def test_unadjusted_mle_fails_closed_for_extremes_instead_of_clamping():
    person, observed, intercepts, slopes = _binary_rows([
        [0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 1, 1, 1, 1],
    ])
    result = score_fixed_calibration_persons(
        person,
        observed,
        intercepts,
        slopes,
        method="MLE",
    ).set_index("Person")

    assert result.loc["P1", "Status"] == "no_finite_likelihood_score_root"
    assert pd.isna(result.loc["P1", "Estimate"])
    assert result.loc["P3", "Status"] == "no_finite_likelihood_score_root"
    assert pd.isna(result.loc["P3", "Estimate"])
    assert result.loc["P2", "Status"] == "ok"
    assert result.loc["P2", "Estimate"] == pytest.approx(np.log(2 / 3), abs=1e-10)


def test_integer_row_weight_matches_literal_row_duplication():
    observed = np.array([0, 2, 1])
    intercepts = np.array([
        [0.0, -0.4, -1.2],
        [0.0, 0.2, -0.7],
        [0.0, -0.1, -0.9],
    ])
    slopes = np.tile(np.array([[0.0, 1.0, 2.0]]), (3, 1))
    weights = np.array([2.0, 3.0, 1.0])
    weighted = score_fixed_calibration_persons(
        ["P1"] * 3,
        observed,
        intercepts,
        slopes,
        row_weights=weights,
    ).iloc[0]
    repeated_index = np.repeat(np.arange(3), weights.astype(int))
    duplicated = score_fixed_calibration_persons(
        ["P1"] * len(repeated_index),
        observed[repeated_index],
        intercepts[repeated_index],
        slopes[repeated_index],
    ).iloc[0]

    assert weighted["Estimate"] == pytest.approx(duplicated["Estimate"], abs=1e-10)
    assert weighted["StandardError"] == pytest.approx(duplicated["StandardError"], abs=1e-10)
    assert weighted["Information"] == pytest.approx(duplicated["Information"], abs=1e-10)


def test_information_derivative_and_weighted_score_match_finite_differences():
    theta = 0.37
    observed = np.array([0, 2, 1, 3])
    intercepts = np.array([
        [0.0, -0.2, -1.0, -2.0],
        [0.0, 0.1, -0.8, -1.9],
        [0.0, -0.5, -0.9, -1.1],
        [0.0, 0.3, -0.2, -1.0],
    ])
    slopes = np.array([
        [0.0, 0.8, 1.6, 2.4],
        [0.0, 1.1, 2.2, 3.3],
        [0.0, 0.6, 1.2, 1.8],
        [0.0, 1.4, 2.8, 4.2],
    ])
    weights = np.array([1.0, 0.75, 2.0, 1.25])
    h = 1e-5
    center = evaluate_fixed_calibration_score(
        theta, observed, intercepts, slopes, row_weights=weights
    )
    lower = evaluate_fixed_calibration_score(
        theta - h, observed, intercepts, slopes, row_weights=weights
    )
    upper = evaluate_fixed_calibration_score(
        theta + h, observed, intercepts, slopes, row_weights=weights
    )
    information_derivative_fd = (upper["Information"] - lower["Information"]) / (2 * h)
    weighted_objective_derivative_fd = (
        upper["WeightedLogLikelihood"] - lower["WeightedLogLikelihood"]
    ) / (2 * h)

    assert center["InformationDerivative"] == pytest.approx(information_derivative_fd, abs=1e-6)
    assert center["AdjustedScore"] == pytest.approx(weighted_objective_derivative_fd, abs=1e-6)


def test_missing_mask_omits_invalid_placeholder_rows_before_validation():
    person = np.array(["P1", "P1", "P1"])
    observed = np.array([0, -999, 1])
    intercepts = np.zeros((3, 2))
    slopes = np.tile(np.array([[0.0, 1.0]]), (3, 1))
    weights = np.array([1.0, np.nan, 1.0])
    result = score_fixed_calibration_persons(
        person,
        observed,
        intercepts,
        slopes,
        row_weights=weights,
        observed_mask=[True, False, True],
    ).iloc[0]

    assert result["Status"] == "ok"
    assert int(result["ObservedRows"]) == 2
    assert float(result["EffectiveWeight"]) == 2.0


def test_invalid_retained_weight_fails_closed():
    person, observed, intercepts, slopes = _binary_rows([[0, 1]])
    with pytest.raises(ValueError, match="strictly positive"):
        score_fixed_calibration_persons(
            person,
            observed,
            intercepts,
            slopes,
            row_weights=[1.0, 0.0],
        )


def test_root_search_retains_one_sided_finite_points_when_opposite_tail_saturates():
    # A large negative intercept places the WLE root far to the right. At the
    # symmetric negative endpoint the binary information has already rounded
    # below the configured minimum; the finite positive-side points must still
    # bracket the root rather than returning an unavailable status.
    person = ["P1"] * 5
    observed = [1, 1, 0, 0, 0]
    intercepts = np.tile(np.array([[0.0, -25.0]]), (5, 1))
    slopes = np.tile(np.array([[0.0, 1.0]]), (5, 1))
    result = score_fixed_calibration_persons(
        person,
        observed,
        intercepts,
        slopes,
    ).iloc[0]

    expected = 25.0 + np.log((2.0 + 0.5) / (5.0 - 2.0 + 0.5))
    assert result["Status"] == "ok"
    assert result["Estimate"] == pytest.approx(expected, abs=1e-10)
    assert result["BracketLower"] >= 0.0


def test_direct_score_evaluator_rejects_invalid_support_and_noninteger_response():
    intercepts = np.zeros((1, 3))
    slopes = np.array([[0.0, 1.0, 2.0]])
    with pytest.raises(ValueError, match="integer category"):
        evaluate_fixed_calibration_score(0.0, [1.5], intercepts, slopes)
    with pytest.raises(ValueError, match="at least two response categories"):
        evaluate_fixed_calibration_score(
            0.0,
            [0],
            intercepts,
            slopes,
            category_available=np.array([[True, False, False]]),
        )
    with pytest.raises(ValueError, match="finite intercept and slope"):
        evaluate_fixed_calibration_score(
            0.0,
            [0],
            np.array([[0.0, np.nan, 0.0]]),
            slopes,
            category_available=np.array([[True, True, False]]),
        )
