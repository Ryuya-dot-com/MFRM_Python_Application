"""Contracts for the repository-only exact-CMLE plus Warm-WLE bridge."""

from __future__ import annotations

import copy

import numpy as np
import pandas as pd
import pytest

from mfrm_app.cmle import fit_cmle
from mfrm_app.cmle_person_scoring import (
    CMLE_WLE_ESTIMATOR_LABEL,
    score_cmle_persons_wle,
)
from mfrm_app.person_scoring import score_fixed_calibration_persons


def _frame(include_extremes: bool = False) -> pd.DataFrame:
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
    return pd.DataFrame(rows, columns=["Person", "Rater", "Criterion", "Score"])


def _fit(model: str = "RSM", include_extremes: bool = False) -> dict:
    return fit_cmle(
        _frame(include_extremes),
        person_col="Person",
        facet_cols=["Rater", "Criterion"],
        score_col="Score",
        rating_min=0,
        rating_max=2,
        model=model,
        step_facet="Criterion" if model == "PCM" else None,
        gtol=1e-8,
    )


def test_cmle_wle_bridge_equals_direct_generic_scoring_and_fitted_surfaces():
    fit = _fit("RSM")
    bridge = score_cmle_persons_wle(fit)
    design = fit["design"]
    coefficient_map = fit["coefficients"].set_index("Parameter")["Estimate"]
    parameters = coefficient_map.reindex(design.parameter_names).to_numpy(dtype=float)
    intercepts = np.einsum("rkp,p->rk", design.row_design, parameters)
    slopes = np.tile(np.arange(design.n_categories, dtype=float), (len(design.data), 1))
    direct = score_fixed_calibration_persons(
        design.data[design.person_col].astype(str),
        design.data[design.score_col].to_numpy(dtype=int) - design.rating_min,
        intercepts,
        slopes,
        person_levels=fit["person_status"]["Person"].astype(str).tolist(),
    )

    assert np.max(np.abs(bridge["Estimate"] - direct["Estimate"])) < 1e-12
    assert np.max(np.abs(bridge["StandardError"] - direct["StandardError"])) < 1e-12
    assert bridge["CombinedEstimator"].eq(CMLE_WLE_ESTIMATOR_LABEL).all()
    assert bridge["CalibrationInferenceReady"].all()
    assert not bridge["CalibrationUncertaintyPropagated"].any()
    assert bridge["WLEPersonScoringReady"].all()

    unit_columns = list(design.facet_cols)
    surface = fit["surfaces"].sort_values(unit_columns + ["InternalCategory"])
    first_rows = design.data.reset_index().drop_duplicates(unit_columns)
    reconstructed = []
    for row in first_rows.itertuples(index=False):
        row_index = int(row.index)
        for category in range(design.n_categories):
            reconstructed.append(
                {
                    **{column: str(getattr(row, column)) for column in unit_columns},
                    "InternalCategory": category,
                    "Reconstructed": intercepts[row_index, category],
                }
            )
    check = surface.merge(
        pd.DataFrame(reconstructed),
        on=[*unit_columns, "InternalCategory"],
        validate="one_to_one",
    )
    assert np.max(
        np.abs(check["LogKernelWithoutPerson"] - check["Reconstructed"])
    ) < 1e-12


def test_pcm_bridge_returns_finite_wle_for_conditional_extreme_persons():
    fit = _fit("PCM", include_extremes=True)
    result = score_cmle_persons_wle(fit).set_index("Person")

    assert bool(fit["summary"].iloc[0]["InferenceReady"])
    assert result.loc["P7", "ConditionalStatus"] == "extreme_low"
    assert result.loc["P8", "ConditionalStatus"] == "extreme_high"
    assert result.loc["P7", "ExtremeScoreDirection"] == "all_minimum"
    assert result.loc["P8", "ExtremeScoreDirection"] == "all_maximum"
    assert np.isfinite(result.loc[["P7", "P8"], "Estimate"]).all()
    assert result.loc[["P7", "P8"], "WLEPersonScoringReady"].all()
    assert result["AdjustedScoreResidual"].max() <= 1e-8


def test_bridge_is_coefficient_order_invariant_but_rejects_identity_loss():
    fit = _fit("RSM")
    original = score_cmle_persons_wle(fit)
    permuted = copy.deepcopy(fit)
    permuted["coefficients"] = permuted["coefficients"].sample(
        frac=1.0, random_state=17
    ).reset_index(drop=True)
    reordered = score_cmle_persons_wle(permuted)
    assert np.max(np.abs(original["Estimate"] - reordered["Estimate"])) < 1e-12

    missing = copy.deepcopy(fit)
    missing["coefficients"] = missing["coefficients"].iloc[:-1].copy()
    with pytest.raises(ValueError, match="coefficient identity"):
        score_cmle_persons_wle(missing)


def test_bridge_rejects_calibration_not_marked_inference_ready():
    fit = _fit("RSM")
    fit["summary"] = fit["summary"].copy()
    fit["summary"].loc[0, "InferenceReady"] = False
    with pytest.raises(ValueError, match="inference-ready"):
        score_cmle_persons_wle(fit)
