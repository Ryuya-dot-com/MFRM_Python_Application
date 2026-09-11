"""Post-result nominal-class decomposition for MnSq threshold surfaces.

This addendum is isolated from :mod:`fit_threshold_sensitivity` so the source
identity recorded by the immutable parent threshold-surface evidence remains
unchanged.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd

from mfrm_app.fit_threshold_sensitivity import (
    CANONICAL_FIT_THRESHOLDS,
    classify_fit_mnsq_array,
    validate_fit_thresholds,
)


FIT_CLASSES = ("overfit", "acceptable", "noisy", "distorting")


def _validated_triplets(
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


def _validated_frame(
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


def evaluate_fit_class_transition_surface(
    frame: pd.DataFrame,
    *,
    statistics: Sequence[str] = ("Infit", "Outfit"),
    threshold_triplets: Iterable[Sequence[float]],
    group_columns: Sequence[str] = ("Model", "Lane"),
    canonical_thresholds: Sequence[float] = CANONICAL_FIT_THRESHOLDS,
) -> pd.DataFrame:
    """Return complete baseline-by-replicate class matrices across thresholds.

    All 16 cells are emitted for every group, statistic, and threshold triplet,
    including zero-count cells. Nominal labels are retained instead of being
    reduced to changed/not-changed or an unsupported ordinal severity score.
    """

    statistic_names = tuple(str(value) for value in statistics)
    if not statistic_names or len(statistic_names) != len(set(statistic_names)):
        raise ValueError("statistics must contain unique column names.")
    groups = tuple(str(value) for value in group_columns)
    triplets = _validated_triplets(threshold_triplets)
    canonical = validate_fit_thresholds(*canonical_thresholds)
    work = _validated_frame(frame, statistic_names, groups)
    rows: list[dict[str, object]] = []
    group_key: str | list[str] = groups[0] if len(groups) == 1 else list(groups)
    for group_values, group in work.groupby(group_key, sort=False, dropna=False):
        group_tuple = group_values if isinstance(group_values, tuple) else (group_values,)
        identity = dict(zip(groups, group_tuple))
        denominator = int(len(group))
        for lower, acceptable, noisy in triplets:
            for statistic in statistic_names:
                replicate_class = classify_fit_mnsq_array(
                    group[statistic].to_numpy(dtype=float),
                    (lower, acceptable, noisy),
                )
                baseline_class = classify_fit_mnsq_array(
                    group[f"Baseline{statistic}"].to_numpy(dtype=float),
                    (lower, acceptable, noisy),
                )
                for baseline_label in FIT_CLASSES:
                    for replicate_label in FIT_CLASSES:
                        count = int(
                            np.sum(
                                (baseline_class == baseline_label)
                                & (replicate_class == replicate_label)
                            )
                        )
                        rows.append(
                            {
                                **identity,
                                "Statistic": statistic,
                                "OverfitUpper": lower,
                                "AcceptableUpper": acceptable,
                                "NoisyUpper": noisy,
                                "CanonicalThresholds": (
                                    lower,
                                    acceptable,
                                    noisy,
                                )
                                == canonical,
                                "BaselineClass": baseline_label,
                                "ReplicateClass": replicate_label,
                                "ClassChanged": baseline_label != replicate_label,
                                "PersonReplicates": denominator,
                                "CellCount": count,
                                "CellShare": float(count / denominator),
                                "ClassificationInput": "finite_unrounded_mnsq",
                                "DependentPersonReplicates": True,
                            }
                        )
    return pd.DataFrame(rows)


__all__ = ["FIT_CLASSES", "evaluate_fit_class_transition_surface"]
