"""Contracts for fixed-calibration CMLE-WLE Person MnSq."""

from __future__ import annotations

import copy

import numpy as np
import pandas as pd

from mfrm_app.cmle import fit_cmle
from mfrm_app.cmle_wle_fit import (
    CMLE_WLE_PERSON_FIT_LABEL,
    PERSON_FIT_REFERENCE,
    compute_cmle_wle_person_fit,
    compute_fixed_calibration_person_fit,
)


def _frame(*, include_extremes: bool = True, missing: bool = False) -> pd.DataFrame:
    scores = {
        "P1": [0, 1, 1, 2],
        "P2": [1, 0, 2, 1],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 1],
        "P5": [0, 2, 1, 1],
        "P6": [2, 0, 1, 1],
    }
    if include_extremes:
        scores["P7"] = [0, 0, 0, 0]
        scores["P8"] = [2, 2, 2, 2]
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    rows = [
        (person, rater, criterion, score)
        for person, values in scores.items()
        for (rater, criterion), score in zip(units, values)
    ]
    frame = pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])
    if missing:
        remove = ((frame["Person"] == "P2") & (frame["Criterion"] == "C2")) | (
            (frame["Person"] == "P5")
            & (frame["Rater"] == "R2")
            & (frame["Criterion"] == "C1")
        )
        frame = frame.loc[~remove].reset_index(drop=True)
    return frame


def _fit(model: str, *, missing: bool = False) -> dict:
    return fit_cmle(
        _frame(missing=missing),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def test_fixed_fit_matches_direct_probability_moments_and_aggregates():
    intercepts = np.array([[0.0, -0.4, -1.1], [0.0, 0.2, -0.3]])
    theta = np.array([0.3, 0.3])
    observed = np.array([1, 2])
    result = compute_fixed_calibration_person_fit(
        ["P1", "P1"],
        observed,
        intercepts,
        theta,
        row_weights=[1.0, 2.0],
    )

    scores = np.arange(3, dtype=float)
    logits = intercepts + theta[:, None] * scores
    probabilities = np.exp(logits - logits.max(axis=1, keepdims=True))
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    expected = probabilities @ scores
    centered = scores[None, :] - expected[:, None]
    variance = np.sum(probabilities * centered**2, axis=1)
    fourth = np.sum(probabilities * centered**4, axis=1)
    residual = observed - expected
    weights = np.array([1.0, 2.0])
    infit = np.sum(weights * residual**2) / np.sum(weights * variance)
    outfit = np.sum(weights * residual**2 / variance) / np.sum(weights)

    observation = result["observations"]
    person = result["persons"].iloc[0]
    assert np.max(np.abs(observation["ExpectedInternalCategory"] - expected)) < 1e-15
    assert np.max(np.abs(observation["Variance"] - variance)) < 1e-15
    assert np.max(np.abs(observation["FourthCentralMoment"] - fourth)) < 1e-15
    assert abs(person["Infit"] - infit) < 1e-15
    assert abs(person["Outfit"] - outfit) < 1e-15
    assert bool(person["PersonFitReady"])
    assert not bool(person["OutlierTrimmingApplied"])


def test_cmle_wle_fit_is_ready_for_rsm_pcm_missingness_and_exact_extremes():
    for model in ("RSM", "PCM"):
        result = compute_cmle_wle_person_fit(_fit(model, missing=True))
        persons = result["persons"].set_index("Person")
        observations = result["observations"]

        assert persons["PersonFitReady"].all()
        assert np.isfinite(persons[["Infit", "Outfit"]]).all().all()
        assert persons.loc["P7", "ExtremeScoreDirection"] == "all_minimum"
        assert persons.loc["P8", "ExtremeScoreDirection"] == "all_maximum"
        assert bool(persons.loc["P7", "ExtremeScorePattern"])
        assert bool(persons.loc["P8", "ExtremeScorePattern"])
        assert observations["FitContributionAvailable"].all()
        assert len(observations) == len(_frame(missing=True))
        assert result["contract"]["EstimatorLabel"] == CMLE_WLE_PERSON_FIT_LABEL
        assert result["contract"]["ReferenceImplementation"] == PERSON_FIT_REFERENCE
        assert not result["contract"]["ZSTDReady"]


def test_cmle_wle_fit_is_row_and_coefficient_order_invariant():
    fit = _fit("PCM", missing=True)
    original = compute_cmle_wle_person_fit(fit)["persons"].sort_values("Person")

    permuted = copy.deepcopy(fit)
    order = np.random.default_rng(29).permutation(len(permuted["design"].data))
    permuted["design"].data = permuted["design"].data.iloc[order].reset_index(drop=True)
    permuted["design"].row_design = permuted["design"].row_design[order]
    permuted["design"].row_offset = permuted["design"].row_offset[order]
    permuted["design"].unit_codes = [permuted["design"].unit_codes[index] for index in order]
    permuted["coefficients"] = permuted["coefficients"].sample(
        frac=1.0, random_state=31
    ).reset_index(drop=True)
    reordered = compute_cmle_wle_person_fit(permuted)["persons"].sort_values("Person")

    for column in ("WLEEstimate", "Infit", "Outfit"):
        assert np.max(np.abs(original[column].to_numpy() - reordered[column].to_numpy())) < 1e-12
    assert original["InfitClass"].tolist() == reordered["InfitClass"].tolist()
    assert original["OutfitClass"].tolist() == reordered["OutfitClass"].tolist()


def test_nonpositive_numerical_variance_is_unavailable_not_zero():
    result = compute_fixed_calibration_person_fit(
        ["P1"],
        [0],
        np.array([[0.0, 1000.0]]),
        [0.0],
    )
    observation = result["observations"].iloc[0]
    person = result["persons"].iloc[0]

    assert observation["Variance"] == 0.0
    assert not bool(observation["FitContributionAvailable"])
    assert not bool(person["PersonFitReady"])
    assert np.isnan(person["Infit"])
    assert np.isnan(person["Outfit"])
    assert person["InfitClass"] == "unavailable"
    assert person["PersonFitStatus"] == "unavailable_nonfinite_or_nonpositive_variance"


def test_unrounded_mnsq_drives_decision_when_display_rounding_crosses_boundary():
    target = 1.5004
    result = compute_fixed_calibration_person_fit(
        ["P1"],
        [0],
        np.array([[0.0, np.log(target)]]),
        [0.0],
        display_decimals=3,
    )
    person = result["persons"].iloc[0]

    assert abs(person["Infit"] - target) < 1e-12
    assert abs(person["Outfit"] - target) < 1e-12
    assert person["InfitClass"] == "noisy"
    assert person["InfitDisplay"] == 1.5
    assert person["InfitBoundaryStatus"] == "display_rounding_boundary"
    assert not bool(person["InfitDisplayDecisionConsistent"])


def test_zstd_and_p_values_remain_explicitly_withheld():
    result = compute_cmle_wle_person_fit(_fit("RSM"))
    persons = result["persons"]

    assert persons["InfitZSTD"].isna().all()
    assert persons["OutfitZSTD"].isna().all()
    assert persons["FitPValue"].isna().all()
    assert persons["ZSTDStatus"].str.startswith("withheld:").all()
    assert persons["PValueStatus"].str.startswith("withheld:").all()
    assert not result["contract"]["PValueReady"]
    assert not result["contract"]["ConfidenceIntervalReady"]
    assert not result["contract"]["TotalStandardErrorReady"]
