"""Contracts for research-only MnSq threshold and display sensitivity."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mfrm_app.decision_stability import classify_fit_mnsq
from mfrm_app.fit_threshold_sensitivity import (
    CANONICAL_FIT_THRESHOLDS,
    classify_fit_mnsq_array,
    classify_fit_mnsq_at_thresholds,
    evaluate_display_precision_sensitivity,
    evaluate_fit_threshold_surface,
    validate_fit_thresholds,
)
from mfrm_app.fit_threshold_severity import (
    FIT_CLASSES,
    evaluate_fit_class_transition_surface,
)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Model": ["RSM"] * 3,
            "Lane": ["fixed"] * 3,
            "Person": ["P1", "P2", "P3"],
            "Infit": [0.49, 0.60, 1.51],
            "BaselineInfit": [0.60, 0.60, 1.40],
            "Outfit": [0.60, 1.00, 1.40],
            "BaselineOutfit": [0.60, 1.00, 1.40],
        }
    )


def test_explicit_threshold_endpoints_and_nextafter_semantics():
    lower, acceptable, noisy = CANONICAL_FIT_THRESHOLDS
    classify = lambda value: classify_fit_mnsq_at_thresholds(
        value,
        overfit_upper=lower,
        acceptable_upper=acceptable,
        noisy_upper=noisy,
    )

    assert classify(np.nextafter(lower, -np.inf)) == "overfit"
    assert classify(lower) == "acceptable"
    assert classify(acceptable) == "acceptable"
    assert classify(np.nextafter(acceptable, np.inf)) == "noisy"
    assert classify(noisy) == "noisy"
    assert classify(np.nextafter(noisy, np.inf)) == "distorting"
    assert classify(np.nan) == "unavailable"


def test_canonical_vector_classifier_matches_application_contract():
    values = np.array([0.2, 0.5, 1.2, 1.5, 1.8, 2.0, 2.1, np.nan])
    expected = np.array([classify_fit_mnsq(value) for value in values], dtype=object)
    observed = classify_fit_mnsq_array(values, CANONICAL_FIT_THRESHOLDS)
    assert np.array_equal(observed, expected)


@pytest.mark.parametrize(
    "thresholds",
    [
        (0.0, 1.5, 2.0),
        (0.5, 0.5, 2.0),
        (0.5, 2.0, 2.0),
        (0.5, np.nan, 2.0),
    ],
)
def test_threshold_triplets_fail_closed_when_not_positive_finite_and_ordered(thresholds):
    with pytest.raises(ValueError, match="threshold"):
        validate_fit_thresholds(*thresholds)


def test_surface_recomputes_raw_transitions_and_any_union():
    surface = evaluate_fit_threshold_surface(
        _frame(),
        threshold_triplets=[(0.5, 1.5, 2.0), (0.45, 1.5, 2.0)],
    )
    canonical = surface.loc[surface["CanonicalThresholds"]].set_index("Statistic")
    alternative = surface.loc[
        surface["OverfitUpper"].eq(0.45)
    ].set_index("Statistic")

    assert int(canonical.loc["Infit", "ClassTransitions"]) == 2
    assert int(canonical.loc["Outfit", "ClassTransitions"]) == 0
    assert int(canonical.loc["Any", "ClassTransitions"]) == 2
    assert int(alternative.loc["Infit", "ClassTransitions"]) == 1
    assert int(alternative.loc["Any", "ClassTransitions"]) == 1
    assert surface["ClassificationInput"].eq("finite_unrounded_mnsq").all()
    assert surface["DependentPersonReplicates"].all()


def test_threshold_surface_is_row_order_invariant():
    triplets = [(0.5, 1.5, 2.0), (0.55, 1.4, 1.9)]
    first = evaluate_fit_threshold_surface(_frame(), threshold_triplets=triplets)
    second = evaluate_fit_threshold_surface(
        _frame().sample(frac=1.0, random_state=19).reset_index(drop=True),
        threshold_triplets=triplets,
    )
    columns = [
        "Model",
        "Lane",
        "Statistic",
        "OverfitUpper",
        "AcceptableUpper",
        "NoisyUpper",
    ]
    pd.testing.assert_frame_equal(
        first.sort_values(columns).reset_index(drop=True),
        second.sort_values(columns).reset_index(drop=True),
    )


def test_display_precision_reports_but_does_not_replace_raw_transition():
    frame = pd.DataFrame(
        {
            "Model": ["RSM"],
            "Lane": ["joint"],
            "Infit": [1.5004],
            "BaselineInfit": [1.4],
            "Outfit": [1.0],
            "BaselineOutfit": [1.0],
        }
    )
    result = evaluate_display_precision_sensitivity(
        frame, decimals=[3, 4]
    )
    infit = result.loc[result["Statistic"].eq("Infit")].set_index(
        "DisplayDecimals"
    )

    assert int(infit.loc[3, "RawClassTransitions"]) == 1
    assert int(infit.loc[3, "RoundedClassTransitions"]) == 0
    assert int(infit.loc[3, "TransitionIndicatorDisagreements"]) == 1
    assert int(infit.loc[3, "ReplicateClassMismatches"]) == 1
    assert int(infit.loc[4, "RoundedClassTransitions"]) == 1
    assert int(infit.loc[4, "TransitionIndicatorDisagreements"]) == 0
    assert result["RawClassificationRetained"].all()


def test_class_transition_surface_emits_zeros_and_matches_binary_surface():
    triplets = [(0.5, 1.5, 2.0), (0.45, 1.5, 2.0)]
    matrices = evaluate_fit_class_transition_surface(
        _frame(), threshold_triplets=triplets
    )
    binary = evaluate_fit_threshold_surface(_frame(), threshold_triplets=triplets)

    assert len(matrices) == 2 * 2 * len(FIT_CLASSES) ** 2
    assert set(matrices["BaselineClass"]) == set(FIT_CLASSES)
    assert set(matrices["ReplicateClass"]) == set(FIT_CLASSES)
    keys = [
        "Model",
        "Lane",
        "Statistic",
        "OverfitUpper",
        "AcceptableUpper",
        "NoisyUpper",
    ]
    totals = matrices.groupby(keys, as_index=False)["CellCount"].sum()
    assert totals["CellCount"].eq(3).all()
    off_diagonal = (
        matrices.loc[matrices["ClassChanged"]]
        .groupby(keys, as_index=False)["CellCount"]
        .sum()
        .rename(columns={"CellCount": "MatrixTransitions"})
    )
    expected = binary.loc[binary["Statistic"].ne("Any"), keys + ["ClassTransitions"]]
    compared = expected.merge(off_diagonal, on=keys, validate="one_to_one")
    assert compared["ClassTransitions"].eq(compared["MatrixTransitions"]).all()
    assert matrices["CellCount"].eq(0).any()


def test_noisy_boundary_can_relabel_severity_while_binary_transition_is_flat():
    frame = pd.DataFrame(
        {
            "Model": ["RSM"] * 3,
            "Lane": ["fixed"] * 3,
            "Infit": [1.92, 2.02, 2.08],
            "BaselineInfit": [1.0, 1.0, 1.0],
        }
    )
    triplets = [(0.5, 1.5, 1.9), (0.5, 1.5, 2.1)]
    binary = evaluate_fit_threshold_surface(
        frame,
        statistics=("Infit",),
        threshold_triplets=triplets,
    )
    matrices = evaluate_fit_class_transition_surface(
        frame,
        statistics=("Infit",),
        threshold_triplets=triplets,
    )

    assert binary.loc[binary["Statistic"].eq("Infit"), "ClassTransitions"].tolist() == [3, 3]
    replicate_counts = (
        matrices.groupby(["NoisyUpper", "ReplicateClass"])["CellCount"]
        .sum()
        .unstack(fill_value=0)
    )
    assert int(replicate_counts.loc[1.9, "distorting"]) == 3
    assert int(replicate_counts.loc[1.9, "noisy"]) == 0
    assert int(replicate_counts.loc[2.1, "distorting"]) == 0
    assert int(replicate_counts.loc[2.1, "noisy"]) == 3
