"""Research-only sensitivity surfaces for Person-fit MnSq thresholds.

This module does not alter the application's canonical 0.50/1.50/2.00
decision contract.  It evaluates alternative, explicitly supplied ordered
threshold triplets on retained unrounded values so that display precision and
substantive threshold dependence can be audited separately.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import pandas as pd


CANONICAL_FIT_THRESHOLDS = (0.5, 1.5, 2.0)


def validate_fit_thresholds(
    overfit_upper: Any,
    acceptable_upper: Any,
    noisy_upper: Any,
) -> tuple[float, float, float]:
    """Return a finite, strictly ordered MnSq threshold triplet."""

    try:
        values = tuple(
            float(value)
            for value in (overfit_upper, acceptable_upper, noisy_upper)
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("Fit thresholds must be finite numeric values.") from exc
    if not np.isfinite(values).all():
        raise ValueError("Fit thresholds must be finite numeric values.")
    if values[0] <= 0.0 or not values[0] < values[1] < values[2]:
        raise ValueError(
            "Fit thresholds must be positive and strictly ordered: "
            "overfit_upper < acceptable_upper < noisy_upper."
        )
    return values


def classify_fit_mnsq_at_thresholds(
    value: Any,
    *,
    overfit_upper: float,
    acceptable_upper: float,
    noisy_upper: float,
) -> str:
    """Classify one raw MnSq with explicit endpoint semantics."""

    thresholds = validate_fit_thresholds(
        overfit_upper, acceptable_upper, noisy_upper
    )
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "unavailable"
    if not np.isfinite(numeric):
        return "unavailable"
    if numeric < thresholds[0]:
        return "overfit"
    if numeric <= thresholds[1]:
        return "acceptable"
    if numeric <= thresholds[2]:
        return "noisy"
    return "distorting"


def classify_fit_mnsq_array(
    values: Sequence[float],
    thresholds: Sequence[float],
) -> np.ndarray:
    """Vectorized raw-value classifier for a validated threshold triplet."""

    if len(thresholds) != 3:
        raise ValueError("thresholds must contain exactly three values.")
    lower, acceptable, noisy = validate_fit_thresholds(*thresholds)
    numeric = np.asarray(values, dtype=float).reshape(-1)
    out = np.full(len(numeric), "unavailable", dtype=object)
    finite = np.isfinite(numeric)
    out[finite & (numeric < lower)] = "overfit"
    out[finite & (numeric >= lower) & (numeric <= acceptable)] = "acceptable"
    out[finite & (numeric > acceptable) & (numeric <= noisy)] = "noisy"
    out[finite & (numeric > noisy)] = "distorting"
    return out


def _threshold_triplets(
    triplets: Iterable[Sequence[float]],
) -> list[tuple[float, float, float]]:
    retained: list[tuple[float, float, float]] = []
    for triplet in triplets:
        if len(triplet) != 3:
            raise ValueError("Every threshold triplet must contain three values.")
        retained.append(validate_fit_thresholds(*triplet))
    if not retained:
        raise ValueError("At least one threshold triplet is required.")
    if len(retained) != len(set(retained)):
        raise ValueError("Threshold triplets must be unique.")
    return retained


def _validated_transition_frame(
    frame: pd.DataFrame,
    statistics: Sequence[str],
    group_columns: Sequence[str],
) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise ValueError("A non-empty Person-draw DataFrame is required.")
    required = set(group_columns)
    for statistic in statistics:
        required.update({statistic, f"Baseline{statistic}"})
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError("Person-draw frame is missing columns: " + ", ".join(missing))
    work = frame.copy()
    for statistic in statistics:
        work[statistic] = pd.to_numeric(work[statistic], errors="coerce")
        work[f"Baseline{statistic}"] = pd.to_numeric(
            work[f"Baseline{statistic}"], errors="coerce"
        )
        if not np.isfinite(
            work[[statistic, f"Baseline{statistic}"]].to_numpy(dtype=float)
        ).all():
            raise ValueError("Ready Person-draw MnSq inputs must be finite.")
    return work


def evaluate_fit_threshold_surface(
    frame: pd.DataFrame,
    *,
    statistics: Sequence[str] = ("Infit", "Outfit"),
    threshold_triplets: Iterable[Sequence[float]],
    group_columns: Sequence[str] = ("Model", "Lane"),
    canonical_thresholds: Sequence[float] = CANONICAL_FIT_THRESHOLDS,
) -> pd.DataFrame:
    """Evaluate per-statistic and any-statistic transition surfaces.

    Every classification is recomputed from the unrounded statistic and its
    unrounded baseline.  The returned shares are descriptive proportions of
    retained Person-replicates, not independent-binomial estimates.
    """

    statistic_names = tuple(str(value) for value in statistics)
    if not statistic_names or len(statistic_names) != len(set(statistic_names)):
        raise ValueError("statistics must contain unique column names.")
    groups = tuple(str(value) for value in group_columns)
    triplets = _threshold_triplets(threshold_triplets)
    canonical = validate_fit_thresholds(*canonical_thresholds)
    work = _validated_transition_frame(frame, statistic_names, groups)
    rows: list[dict[str, object]] = []
    group_key: str | list[str] = groups[0] if len(groups) == 1 else list(groups)
    for group_values, group in work.groupby(group_key, sort=False, dropna=False):
        group_tuple = group_values if isinstance(group_values, tuple) else (group_values,)
        identity = dict(zip(groups, group_tuple))
        for lower, acceptable, noisy in triplets:
            statistic_transition_masks: list[np.ndarray] = []
            for statistic in statistic_names:
                replicate_class = classify_fit_mnsq_array(
                    group[statistic].to_numpy(dtype=float),
                    (lower, acceptable, noisy),
                )
                baseline_class = classify_fit_mnsq_array(
                    group[f"Baseline{statistic}"].to_numpy(dtype=float),
                    (lower, acceptable, noisy),
                )
                changed = replicate_class != baseline_class
                statistic_transition_masks.append(changed)
                rows.append(
                    {
                        **identity,
                        "Statistic": statistic,
                        "OverfitUpper": lower,
                        "AcceptableUpper": acceptable,
                        "NoisyUpper": noisy,
                        "CanonicalThresholds": (lower, acceptable, noisy)
                        == canonical,
                        "PersonReplicates": int(len(group)),
                        "ClassTransitions": int(changed.sum()),
                        "TransitionShare": float(changed.mean()),
                        "ClassificationInput": "finite_unrounded_mnsq",
                        "DependentPersonReplicates": True,
                    }
                )
            any_changed = np.logical_or.reduce(statistic_transition_masks)
            rows.append(
                {
                    **identity,
                    "Statistic": "Any",
                    "OverfitUpper": lower,
                    "AcceptableUpper": acceptable,
                    "NoisyUpper": noisy,
                    "CanonicalThresholds": (lower, acceptable, noisy)
                    == canonical,
                    "PersonReplicates": int(len(group)),
                    "ClassTransitions": int(any_changed.sum()),
                    "TransitionShare": float(any_changed.mean()),
                    "ClassificationInput": "finite_unrounded_mnsq",
                    "DependentPersonReplicates": True,
                }
            )
    return pd.DataFrame(rows)


def evaluate_display_precision_sensitivity(
    frame: pd.DataFrame,
    *,
    decimals: Iterable[int],
    statistics: Sequence[str] = ("Infit", "Outfit"),
    group_columns: Sequence[str] = ("Model", "Lane"),
    thresholds: Sequence[float] = CANONICAL_FIT_THRESHOLDS,
) -> pd.DataFrame:
    """Compare raw and rounded class/transition decisions by precision."""

    threshold_values = validate_fit_thresholds(*thresholds)
    precisions: list[int] = []
    for value in decimals:
        if isinstance(value, (bool, np.bool_)):
            raise ValueError("display decimals must be non-negative integers.")
        numeric = int(value)
        if float(value) != numeric or numeric < 0:
            raise ValueError("display decimals must be non-negative integers.")
        precisions.append(numeric)
    if not precisions or len(precisions) != len(set(precisions)):
        raise ValueError("display decimals must be a non-empty unique sequence.")
    statistic_names = tuple(str(value) for value in statistics)
    groups = tuple(str(value) for value in group_columns)
    work = _validated_transition_frame(frame, statistic_names, groups)
    rows: list[dict[str, object]] = []
    group_key: str | list[str] = groups[0] if len(groups) == 1 else list(groups)
    for group_values, group in work.groupby(group_key, sort=False, dropna=False):
        group_tuple = group_values if isinstance(group_values, tuple) else (group_values,)
        identity = dict(zip(groups, group_tuple))
        for precision in precisions:
            raw_any = np.zeros(len(group), dtype=bool)
            rounded_any = np.zeros(len(group), dtype=bool)
            for statistic in statistic_names:
                raw = group[statistic].to_numpy(dtype=float)
                baseline_raw = group[f"Baseline{statistic}"].to_numpy(dtype=float)
                raw_class = classify_fit_mnsq_array(raw, threshold_values)
                baseline_raw_class = classify_fit_mnsq_array(
                    baseline_raw, threshold_values
                )
                rounded_class = classify_fit_mnsq_array(
                    np.round(raw, decimals=precision), threshold_values
                )
                baseline_rounded_class = classify_fit_mnsq_array(
                    np.round(baseline_raw, decimals=precision), threshold_values
                )
                raw_transition = raw_class != baseline_raw_class
                rounded_transition = rounded_class != baseline_rounded_class
                raw_any |= raw_transition
                rounded_any |= rounded_transition
                rows.append(
                    {
                        **identity,
                        "DisplayDecimals": precision,
                        "Statistic": statistic,
                        "PersonReplicates": int(len(group)),
                        "RawClassTransitions": int(raw_transition.sum()),
                        "RoundedClassTransitions": int(rounded_transition.sum()),
                        "RoundedMinusRawTransitions": int(
                            rounded_transition.sum() - raw_transition.sum()
                        ),
                        "TransitionIndicatorDisagreements": int(
                            np.sum(raw_transition != rounded_transition)
                        ),
                        "ReplicateClassMismatches": int(
                            np.sum(raw_class != rounded_class)
                        ),
                        "BaselineClassMismatches": int(
                            np.sum(baseline_raw_class != baseline_rounded_class)
                        ),
                        "RawClassificationRetained": True,
                    }
                )
            rows.append(
                {
                    **identity,
                    "DisplayDecimals": precision,
                    "Statistic": "Any",
                    "PersonReplicates": int(len(group)),
                    "RawClassTransitions": int(raw_any.sum()),
                    "RoundedClassTransitions": int(rounded_any.sum()),
                    "RoundedMinusRawTransitions": int(
                        rounded_any.sum() - raw_any.sum()
                    ),
                    "TransitionIndicatorDisagreements": int(
                        np.sum(raw_any != rounded_any)
                    ),
                    "ReplicateClassMismatches": np.nan,
                    "BaselineClassMismatches": np.nan,
                    "RawClassificationRetained": True,
                }
            )
    return pd.DataFrame(rows)


__all__ = [
    "CANONICAL_FIT_THRESHOLDS",
    "classify_fit_mnsq_array",
    "classify_fit_mnsq_at_thresholds",
    "evaluate_display_precision_sensitivity",
    "evaluate_fit_threshold_surface",
    "validate_fit_thresholds",
]
