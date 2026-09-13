"""Calibration-draw sensitivity contracts for exact CMLE plus Warm WLE."""

from __future__ import annotations

import copy

import numpy as np
import pandas as pd
import pytest

from mfrm_app.cmle import fit_cmle
from mfrm_app.cmle_wle_uncertainty import cmle_wle_calibration_sensitivity


def _fit(model: str = "RSM", include_extremes: bool = False) -> dict:
    scores = {
        "P1": [0, 0, 1, 2],
        "P2": [1, 0, 2, 2],
        "P3": [2, 1, 1, 0],
        "P4": [1, 2, 0, 2],
        "P5": [0, 2, 1, 0],
        "P6": [2, 0, 2, 1],
    }
    if include_extremes:
        scores.update({"P7": [0, 0, 0, 0], "P8": [2, 2, 2, 2]})
    units = [("R1", "C1"), ("R1", "C2"), ("R2", "C1"), ("R2", "C2")]
    frame = pd.DataFrame(
        [
            (person, rater, criterion, score)
            for person, values in scores.items()
            for (rater, criterion), score in zip(units, values)
        ],
        columns=["Person", "Rater", "Criterion", "Score"],
    )
    return fit_cmle(
        frame,
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def test_zero_covariance_scale_returns_zero_calibration_dispersion():
    result = cmle_wle_calibration_sensitivity(
        _fit(), n_draws=40, seed=11, covariance_scale=0.0
    )
    persons = result["persons"]
    assert np.max(persons["CalibrationDrawSD"].abs()) < 1e-12
    assert np.max(persons["CalibrationDrawMeanMinusBaseline"].abs()) < 1e-12
    assert np.max(
        np.abs(persons["QuadratureSensitivitySE"] - persons["StandardError"])
    ) < 1e-12
    assert not persons["CalibrationDrawIntervalIsConfidenceInterval"].any()
    assert not persons["QuadratureSensitivityIsInferenceQualified"].any()


def test_same_seed_is_exactly_reproducible_and_draw_retention_is_optional():
    fit = _fit()
    first = cmle_wle_calibration_sensitivity(
        fit, n_draws=80, seed=20260810, retain_draws=True
    )
    second = cmle_wle_calibration_sensitivity(
        fit, n_draws=80, seed=20260810, retain_draws=False
    )
    columns = [
        "CalibrationDrawMean",
        "CalibrationDrawSD",
        "CalibrationDrawQ025",
        "CalibrationDrawQ500",
        "CalibrationDrawQ975",
    ]
    assert np.array_equal(
        first["persons"][columns].to_numpy(),
        second["persons"][columns].to_numpy(),
    )
    assert len(first["draws"]) == 80 * len(first["persons"])
    assert second["draws"].empty


def test_covariance_scale_increases_median_draw_sensitivity():
    fit = _fit()
    medians = []
    for scale in (0.5, 1.0, 2.0):
        result = cmle_wle_calibration_sensitivity(
            fit, n_draws=500, seed=20260810, covariance_scale=scale
        )
        medians.append(float(result["summary"].iloc[0]["MedianCalibrationDrawSD"]))
    assert medians[0] < medians[1] < medians[2]


def test_pcm_exact_extremes_remain_finite_across_calibration_draws():
    result = cmle_wle_calibration_sensitivity(
        _fit("PCM", include_extremes=True),
        n_draws=100,
        seed=19,
        covariance_scale=1.0,
    )
    persons = result["persons"].set_index("Person")
    assert bool(result["summary"].iloc[0]["AllDrawPersonScoresSuccessful"])
    assert persons.loc[["P7", "P8"], "ExtremeScorePattern"].all()
    assert np.isfinite(persons.loc[["P7", "P8"], "CalibrationDrawSD"]).all()
    assert persons.loc[["P7", "P8"], "CalibrationDrawsSuccessful"].eq(100).all()


def test_resource_and_covariance_fail_closed_checks():
    fit = _fit()
    with pytest.raises(ValueError, match="row-draw work"):
        cmle_wle_calibration_sensitivity(
            fit,
            n_draws=100,
            seed=1,
            maximum_row_draw_work=100,
        )
    invalid = copy.deepcopy(fit)
    invalid["covariance"] = np.array(invalid["covariance"], copy=True)
    invalid["covariance"][0, 0] = -10.0
    with pytest.raises(ValueError, match="non-positive-semidefinite"):
        cmle_wle_calibration_sensitivity(invalid, n_draws=10, seed=1)

    asymmetric = copy.deepcopy(fit)
    asymmetric["covariance"] = np.array(asymmetric["covariance"], copy=True)
    asymmetric["covariance"][0, 1] += 0.01
    with pytest.raises(ValueError, match="non-symmetric"):
        cmle_wle_calibration_sensitivity(asymmetric, n_draws=10, seed=1)
