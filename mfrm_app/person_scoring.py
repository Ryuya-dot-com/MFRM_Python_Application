"""Fixed-calibration Person scoring for one-dimensional polytomous models.

The module is deliberately independent of Streamlit and joint calibration.
Warm weighted-likelihood estimates (WLEs) are adjusted-score Person summaries,
not finite JMLE maximum-likelihood estimates and not MML/EAP scores.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import brentq


ESTIMATOR_LABEL_WLE = "WLE_fixed_calibration"
ESTIMATOR_LABEL_MLE = "MLE_fixed_calibration"


def _validated_inputs(
    person_index: Sequence[object],
    observed_category: Sequence[int],
    category_intercepts: np.ndarray,
    category_slopes: np.ndarray,
    *,
    row_weights: Sequence[float] | None,
    category_available: np.ndarray | None,
    observed_mask: Sequence[bool] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    persons = np.asarray(person_index, dtype=object).reshape(-1)
    observed = np.asarray(observed_category).reshape(-1)
    intercepts = np.asarray(category_intercepts, dtype=float)
    slopes = np.asarray(category_slopes, dtype=float)
    if intercepts.ndim != 2 or slopes.ndim != 2 or intercepts.shape != slopes.shape:
        raise ValueError("category_intercepts and category_slopes must be equal-shaped 2D arrays.")
    n_rows, n_categories = intercepts.shape
    if n_categories < 2:
        raise ValueError("Fixed-calibration Person scoring requires at least two categories.")
    if len(persons) != n_rows or len(observed) != n_rows:
        raise ValueError("Person, observed-category, and category-array row counts must match.")
    if row_weights is None:
        weights = np.ones(n_rows, dtype=float)
    else:
        weights = np.asarray(row_weights, dtype=float).reshape(-1)
        if len(weights) != n_rows:
            raise ValueError("row_weights length must match category-array rows.")
    available = (
        np.isfinite(intercepts) & np.isfinite(slopes)
        if category_available is None
        else np.asarray(category_available, dtype=bool)
    )
    if available.shape != intercepts.shape:
        raise ValueError("category_available must match category-array shape.")
    mask = (
        np.ones(n_rows, dtype=bool)
        if observed_mask is None
        else np.asarray(observed_mask, dtype=bool).reshape(-1)
    )
    if len(mask) != n_rows:
        raise ValueError("observed_mask length must match category-array rows.")
    persons = persons[mask]
    observed = observed[mask]
    intercepts = intercepts[mask]
    slopes = slopes[mask]
    weights = weights[mask]
    available = available[mask]
    if not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError("Every retained row weight must be finite and strictly positive.")
    if not np.all(available.sum(axis=1) >= 2):
        raise ValueError("Every retained row must support at least two response categories.")
    if not np.all(np.equal(observed, np.floor(observed))):
        raise ValueError("observed_category must contain integer category indices.")
    observed = observed.astype(int)
    if np.any((observed < 0) | (observed >= n_categories)):
        raise ValueError("observed_category contains an index outside the category arrays.")
    if np.any(~available[np.arange(len(observed)), observed]):
        raise ValueError("An observed category is unavailable under its row-specific support.")
    intercepts = np.where(available, intercepts, -np.inf)
    slopes = np.where(available, slopes, 0.0)
    return persons, observed, intercepts, slopes, weights, available


def _probabilities(
    theta: float,
    intercepts: np.ndarray,
    slopes: np.ndarray,
    available: np.ndarray,
) -> np.ndarray:
    logits = intercepts + float(theta) * slopes
    logits = np.where(available, logits, -np.inf)
    row_max = np.max(logits, axis=1, keepdims=True)
    exponentiated = np.where(available, np.exp(logits - row_max), 0.0)
    denominator = exponentiated.sum(axis=1, keepdims=True)
    if np.any(~np.isfinite(denominator) | (denominator <= 0)):
        raise FloatingPointError("Category probabilities could not be normalized.")
    return exponentiated / denominator


def evaluate_fixed_calibration_score(
    theta: float,
    observed_category: Sequence[int],
    category_intercepts: np.ndarray,
    category_slopes: np.ndarray,
    *,
    row_weights: Sequence[float] | None = None,
    category_available: np.ndarray | None = None,
    weighted_likelihood: bool = True,
    information_minimum: float = 1e-12,
) -> dict[str, float]:
    """Evaluate likelihood and Warm adjusted-score moments at ``theta``."""
    if not np.isfinite(theta):
        raise ValueError("theta must be finite.")
    if not np.isfinite(information_minimum) or float(information_minimum) < 0:
        raise ValueError("information_minimum must be finite and non-negative.")
    observed_raw = np.asarray(observed_category).reshape(-1)
    intercepts = np.asarray(category_intercepts, dtype=float)
    slopes = np.asarray(category_slopes, dtype=float)
    if intercepts.shape != slopes.shape or intercepts.ndim != 2:
        raise ValueError("category arrays must be equal-shaped 2D arrays.")
    if intercepts.shape[1] < 2:
        raise ValueError("Fixed-calibration Person scoring requires at least two categories.")
    if len(observed_raw) != intercepts.shape[0]:
        raise ValueError("observed_category length must match category-array rows.")
    if not np.all(np.equal(observed_raw, np.floor(observed_raw))):
        raise ValueError("observed_category must contain integer category indices.")
    observed = observed_raw.astype(int)
    weights = (
        np.ones(len(observed), dtype=float)
        if row_weights is None else np.asarray(row_weights, dtype=float).reshape(-1)
    )
    if len(weights) != len(observed) or not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError("row_weights must be finite, positive, and row-aligned.")
    available = (
        np.isfinite(intercepts) & np.isfinite(slopes)
        if category_available is None else np.asarray(category_available, dtype=bool)
    )
    if available.shape != intercepts.shape:
        raise ValueError("category_available must match category-array shape.")
    if not np.all(available.sum(axis=1) >= 2):
        raise ValueError("Every row must support at least two response categories.")
    if np.any(available & (~np.isfinite(intercepts) | ~np.isfinite(slopes))):
        raise ValueError("Every available category must have finite intercept and slope values.")
    if np.any((observed < 0) | (observed >= intercepts.shape[1])):
        raise ValueError("observed_category contains an index outside the category arrays.")
    if np.any(~available[np.arange(len(observed)), observed]):
        raise ValueError("An observed category is unavailable under its row-specific support.")
    intercepts = np.where(available, intercepts, -np.inf)
    slopes = np.where(available, slopes, 0.0)
    probabilities = _probabilities(theta, intercepts, slopes, available)
    mean = np.sum(probabilities * slopes, axis=1)
    centered = slopes - mean[:, None]
    variance = np.sum(probabilities * centered**2, axis=1)
    third_central = np.sum(probabilities * centered**3, axis=1)
    observed_sufficient = slopes[np.arange(len(observed)), observed]
    likelihood_score = float(np.sum(weights * (observed_sufficient - mean)))
    information = float(np.sum(weights * variance))
    information_derivative = float(np.sum(weights * third_central))
    if not np.isfinite(information) or information <= float(information_minimum):
        warm_adjustment = np.nan
        adjusted_score = np.nan if weighted_likelihood else likelihood_score
    else:
        warm_adjustment = 0.5 * information_derivative / information
        adjusted_score = likelihood_score + warm_adjustment if weighted_likelihood else likelihood_score
    observed_probability = probabilities[np.arange(len(observed)), observed]
    with np.errstate(divide="ignore"):
        log_likelihood = float(np.sum(weights * np.log(observed_probability)))
    weighted_log_likelihood = (
        log_likelihood + 0.5 * np.log(information)
        if weighted_likelihood and information > float(information_minimum)
        else log_likelihood
    )
    return {
        "LikelihoodScore": likelihood_score,
        "Information": information,
        "InformationDerivative": information_derivative,
        "WarmAdjustment": float(warm_adjustment),
        "AdjustedScore": float(adjusted_score),
        "LogLikelihood": log_likelihood,
        "WeightedLogLikelihood": float(weighted_log_likelihood),
    }


def _score_one_person(
    observed: np.ndarray,
    intercepts: np.ndarray,
    slopes: np.ndarray,
    weights: np.ndarray,
    available: np.ndarray,
    *,
    weighted_likelihood: bool,
    maximum_absolute_theta: float,
    root_absolute_tolerance: float,
    root_relative_tolerance: float,
    score_residual_tolerance: float,
    information_minimum: float,
) -> dict[str, object]:
    method_label = ESTIMATOR_LABEL_WLE if weighted_likelihood else ESTIMATOR_LABEL_MLE

    def evaluate(theta: float) -> dict[str, float]:
        return evaluate_fixed_calibration_score(
            theta,
            observed,
            intercepts,
            slopes,
            row_weights=weights,
            category_available=available,
            weighted_likelihood=weighted_likelihood,
            information_minimum=information_minimum,
        )

    minimum_categories = np.argmax(available, axis=1)
    maximum_categories = available.shape[1] - 1 - np.argmax(available[:, ::-1], axis=1)
    all_minimum = bool(np.all(observed == minimum_categories))
    all_maximum = bool(np.all(observed == maximum_categories))
    direction = "all_minimum" if all_minimum else "all_maximum" if all_maximum else "interior"

    # Exact extreme patterns have no finite unadjusted MLE.  Detect them from
    # the response pattern rather than accepting a floating-point zero after
    # the probabilities saturate at a large finite theta.
    if not weighted_likelihood and (all_minimum or all_maximum):
        return {
            "Estimator": method_label,
            "Status": "no_finite_likelihood_score_root",
            "Estimate": np.nan,
            "StandardError": np.nan,
            "ExtremeScorePattern": True,
            "ExtremeScoreDirection": direction,
            "AdjustedScoreResidual": np.nan,
            "LikelihoodScore": np.nan,
            "Information": np.nan,
            "InformationDerivative": np.nan,
            "WarmAdjustment": np.nan,
            "BracketLower": np.nan,
            "BracketUpper": np.nan,
            "Iterations": 0,
            "EstimateRole": "fixed_calibration_mle",
        }

    bracket_lower = np.nan
    bracket_upper = np.nan
    radius = min(4.0, float(maximum_absolute_theta))
    root_result = None
    finite_scores: dict[float, float] = {}

    def register(theta: float) -> None:
        score = evaluate(theta)["AdjustedScore"]
        if np.isfinite(score):
            finite_scores[float(theta)] = float(score)

    register(0.0)
    while radius <= float(maximum_absolute_theta) + np.finfo(float).eps:
        register(-float(radius))
        register(float(radius))
        ordered = sorted(finite_scores.items())
        exact = next(((theta, score) for theta, score in ordered if score == 0.0), None)
        if exact is not None:
            estimate = exact[0]
            bracket_lower = bracket_upper = exact[0]
            break
        sign_change = next(
            (
                (left_theta, right_theta)
                for (left_theta, left_score), (right_theta, right_score)
                in zip(ordered[:-1], ordered[1:])
                if left_score * right_score < 0
            ),
            None,
        )
        if sign_change is not None:
            bracket_lower, bracket_upper = sign_change
            estimate, root_result = brentq(
                lambda value: evaluate(value)["AdjustedScore"],
                bracket_lower,
                bracket_upper,
                xtol=float(root_absolute_tolerance),
                rtol=max(float(root_relative_tolerance), 4 * np.finfo(float).eps),
                full_output=True,
                disp=False,
            )
            break
        if radius >= float(maximum_absolute_theta):
            estimate = np.nan
            break
        radius = min(float(maximum_absolute_theta), radius * 2.0)
    else:
        estimate = np.nan

    if not np.isfinite(estimate):
        return {
            "Estimator": method_label,
            "Status": "no_finite_adjusted_score_root" if weighted_likelihood else "no_finite_likelihood_score_root",
            "Estimate": np.nan,
            "StandardError": np.nan,
            "ExtremeScorePattern": bool(all_minimum or all_maximum),
            "ExtremeScoreDirection": direction,
            "AdjustedScoreResidual": np.nan,
            "LikelihoodScore": np.nan,
            "Information": np.nan,
            "InformationDerivative": np.nan,
            "WarmAdjustment": np.nan,
            "BracketLower": bracket_lower,
            "BracketUpper": bracket_upper,
            "Iterations": 0,
            "EstimateRole": (
                "fixed_calibration_wle_not_jmle_mle"
                if weighted_likelihood else "fixed_calibration_mle"
            ),
        }
    final = evaluate(float(estimate))
    residual = abs(float(final["AdjustedScore"]))
    information = float(final["Information"])
    standard_error = (
        float(np.sqrt(1.0 / information))
        if np.isfinite(information) and information > float(information_minimum)
        else np.nan
    )
    status = "ok" if residual <= float(score_residual_tolerance) else "adjusted_score_residual_above_tolerance"
    return {
        "Estimator": method_label,
        "Status": status,
        "Estimate": float(estimate),
        "StandardError": standard_error,
        "ExtremeScorePattern": bool(all_minimum or all_maximum),
        "ExtremeScoreDirection": direction,
        "AdjustedScoreResidual": residual,
        "LikelihoodScore": float(final["LikelihoodScore"]),
        "Information": information,
        "InformationDerivative": float(final["InformationDerivative"]),
        "WarmAdjustment": float(final["WarmAdjustment"]),
        "BracketLower": bracket_lower,
        "BracketUpper": bracket_upper,
        "Iterations": int(getattr(root_result, "iterations", 0) or 0),
        "EstimateRole": (
            "fixed_calibration_wle_not_jmle_mle"
            if weighted_likelihood else "fixed_calibration_mle"
        ),
    }


def score_fixed_calibration_persons(
    person_index: Sequence[object],
    observed_category: Sequence[int],
    category_intercepts: np.ndarray,
    category_slopes: np.ndarray,
    *,
    row_weights: Sequence[float] | None = None,
    category_available: np.ndarray | None = None,
    observed_mask: Sequence[bool] | None = None,
    person_levels: Iterable[object] | None = None,
    method: str = "WLE",
    maximum_absolute_theta: float = 60.0,
    root_absolute_tolerance: float = 1e-10,
    root_relative_tolerance: float = 1e-12,
    score_residual_tolerance: float = 1e-8,
    information_minimum: float = 1e-12,
) -> pd.DataFrame:
    """Score Persons with fixed category intercepts and theta coefficients."""
    if str(method).upper() not in {"WLE", "MLE"}:
        raise ValueError("method must be 'WLE' or 'MLE'.")
    if not np.isfinite(maximum_absolute_theta) or float(maximum_absolute_theta) <= 0:
        raise ValueError("maximum_absolute_theta must be finite and positive.")
    persons, observed, intercepts, slopes, weights, available = _validated_inputs(
        person_index,
        observed_category,
        category_intercepts,
        category_slopes,
        row_weights=row_weights,
        category_available=category_available,
        observed_mask=observed_mask,
    )
    levels = (
        list(dict.fromkeys(persons.tolist()))
        if person_levels is None else list(person_levels)
    )
    rows = []
    for person in levels:
        selected = np.fromiter((value == person for value in persons), dtype=bool, count=len(persons))
        if not np.any(selected):
            rows.append({
                "Person": person,
                "Estimator": ESTIMATOR_LABEL_WLE if str(method).upper() == "WLE" else ESTIMATOR_LABEL_MLE,
                "Status": "no_observed_rows",
                "Estimate": np.nan,
                "StandardError": np.nan,
                "ObservedRows": 0,
                "EffectiveWeight": 0.0,
                "ExtremeScorePattern": pd.NA,
                "ExtremeScoreDirection": "unavailable",
                "AdjustedScoreResidual": np.nan,
                "LikelihoodScore": np.nan,
                "Information": np.nan,
                "InformationDerivative": np.nan,
                "WarmAdjustment": np.nan,
                "BracketLower": np.nan,
                "BracketUpper": np.nan,
                "Iterations": 0,
                "EstimateRole": "unavailable",
            })
            continue
        scored = _score_one_person(
            observed[selected],
            intercepts[selected],
            slopes[selected],
            weights[selected],
            available[selected],
            weighted_likelihood=str(method).upper() == "WLE",
            maximum_absolute_theta=maximum_absolute_theta,
            root_absolute_tolerance=root_absolute_tolerance,
            root_relative_tolerance=root_relative_tolerance,
            score_residual_tolerance=score_residual_tolerance,
            information_minimum=information_minimum,
        )
        rows.append({
            "Person": person,
            **scored,
            "ObservedRows": int(selected.sum()),
            "EffectiveWeight": float(weights[selected].sum()),
        })
    return pd.DataFrame(rows)


__all__ = [
    "ESTIMATOR_LABEL_MLE",
    "ESTIMATOR_LABEL_WLE",
    "evaluate_fixed_calibration_score",
    "score_fixed_calibration_persons",
]
