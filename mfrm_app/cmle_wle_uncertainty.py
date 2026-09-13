"""Calibration-draw sensitivity for the repository-only CMLE-WLE bridge.

This module samples the asymptotic exact-CMLE free-coordinate covariance and
rescored fit-sample Persons with fixed-calibration Warm WLE.  The resulting
dispersion is a sensitivity diagnostic, not an inference-qualified total
standard error or confidence interval: calibration and Person scores reuse the
same responses, and their dependence is not represented by independent normal
coefficient draws.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from mfrm_app.cmle_person_scoring import (
    _ready_cmle_inputs,
    score_cmle_persons_wle,
)
from mfrm_app.person_scoring import score_fixed_calibration_persons


CALIBRATION_SENSITIVITY_SCOPE = (
    "asymptotic CMLE free-coordinate draws; same-sample dependence omitted; "
    "not an inference-qualified total SE or interval"
)


def _positive_integer(value: object, name: str) -> int:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive integer.") from exc
    if not np.isfinite(numeric) or numeric < 1 or abs(numeric - round(numeric)) > 1e-12:
        raise ValueError(f"{name} must be a positive integer.")
    return int(round(numeric))


def _validated_quantiles(values: Sequence[float]) -> np.ndarray:
    quantiles = np.asarray(values, dtype=float).reshape(-1)
    if (
        len(quantiles) < 2
        or not np.isfinite(quantiles).all()
        or np.any((quantiles <= 0.0) | (quantiles >= 1.0))
        or np.any(np.diff(quantiles) <= 0.0)
    ):
        raise ValueError("quantiles must be finite, strictly increasing values inside (0, 1).")
    return quantiles


def _covariance_square_root(covariance: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    covariance = np.asarray(covariance, dtype=float)
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("CMLE covariance must be square.")
    if not np.isfinite(covariance).all():
        raise ValueError("CMLE covariance must be finite.")
    asymmetry = float(np.max(np.abs(covariance - covariance.T))) if covariance.size else 0.0
    entry_scale = max(float(np.max(np.abs(covariance))), 1.0) if covariance.size else 1.0
    symmetry_tolerance = max(1e-12, entry_scale * 1e-10)
    if asymmetry > symmetry_tolerance:
        raise ValueError(
            "CMLE covariance is materially non-symmetric; "
            f"maximum asymmetry={asymmetry:.6g}, tolerance={symmetry_tolerance:.6g}."
        )
    symmetric = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    tolerance = max(1e-12, scale * 1e-10)
    minimum = float(np.min(eigenvalues))
    if minimum < -tolerance:
        raise ValueError(
            "CMLE covariance is materially non-positive-semidefinite; "
            f"minimum eigenvalue={minimum:.6g}, tolerance={tolerance:.6g}."
        )
    clipped = np.clip(eigenvalues, 0.0, None)
    square_root = eigenvectors @ np.diag(np.sqrt(clipped))
    return square_root, {
        "CovarianceAsymmetryMaxAbs": asymmetry,
        "CovarianceSymmetryTolerance": symmetry_tolerance,
        "CovarianceMinimumEigenvalue": minimum,
        "CovarianceEigenvalueTolerance": tolerance,
        "CovarianceClippedEigenvalues": int(np.sum(eigenvalues < 0.0)),
    }


def cmle_wle_calibration_sensitivity(
    cmle_result: dict[str, object],
    *,
    n_draws: int = 1000,
    seed: int,
    covariance_scale: float = 1.0,
    quantiles: Sequence[float] = (0.025, 0.5, 0.975),
    maximum_draws: int = 10_000,
    maximum_row_draw_work: int = 5_000_000,
    retain_draws: bool = False,
) -> dict[str, object]:
    """Rescore fit-sample Persons over asymptotic CMLE coefficient draws."""
    draws_requested = _positive_integer(n_draws, "n_draws")
    maximum_draws = _positive_integer(maximum_draws, "maximum_draws")
    maximum_row_draw_work = _positive_integer(
        maximum_row_draw_work, "maximum_row_draw_work"
    )
    if draws_requested > maximum_draws:
        raise ValueError(
            f"n_draws={draws_requested} exceeds maximum_draws={maximum_draws}."
        )
    try:
        seed_value = int(seed)
    except (TypeError, ValueError) as exc:
        raise ValueError("seed must be an explicit integer.") from exc
    if isinstance(seed, float) and (not np.isfinite(seed) or seed != seed_value):
        raise ValueError("seed must be an explicit integer.")
    covariance_scale = float(covariance_scale)
    if not np.isfinite(covariance_scale) or covariance_scale < 0.0:
        raise ValueError("covariance_scale must be finite and non-negative.")
    quantile_values = _validated_quantiles(quantiles)

    design, parameters = _ready_cmle_inputs(cmle_result)
    work = int(len(design.data) * draws_requested)
    if work > maximum_row_draw_work:
        raise ValueError(
            f"Requested row-draw work {work} exceeds maximum_row_draw_work="
            f"{maximum_row_draw_work}."
        )
    covariance = cmle_result.get("covariance")
    if covariance is None:
        raise ValueError("Inference-ready CMLE covariance is unavailable.")
    covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (design.n_parameters, design.n_parameters):
        raise ValueError("CMLE covariance dimensions differ from the fitted design.")
    square_root, covariance_audit = _covariance_square_root(covariance)

    baseline = score_cmle_persons_wle(cmle_result).copy()
    person_levels = baseline["Person"].astype(str).tolist()
    persons_by_row = design.data[design.person_col].astype(str).to_numpy(dtype=object)
    observed = (
        pd.to_numeric(design.data[design.score_col], errors="raise").to_numpy(dtype=int)
        - int(design.rating_min)
    )
    category_slopes = np.tile(
        np.arange(design.n_categories, dtype=float),
        (len(design.data), 1),
    )
    weights = (
        pd.to_numeric(design.data["__cmle_weight__"], errors="raise").to_numpy(dtype=float)
        if "__cmle_weight__" in design.data.columns
        else np.ones(len(design.data), dtype=float)
    )

    rng = np.random.default_rng(seed_value)
    standard_normal = rng.standard_normal((draws_requested, design.n_parameters))
    parameter_draws = (
        parameters[None, :]
        + covariance_scale * (standard_normal @ square_root.T)
    )
    estimates = np.full((draws_requested, len(person_levels)), np.nan, dtype=float)
    conditional_ses = np.full_like(estimates, np.nan)
    statuses = np.empty((draws_requested, len(person_levels)), dtype=object)
    draw_parts: list[pd.DataFrame] = []
    for draw_number, draw in enumerate(parameter_draws):
        intercepts = design.row_offset + np.einsum(
            "rkp,p->rk", design.row_design, draw, optimize=True
        )
        scored = score_fixed_calibration_persons(
            persons_by_row,
            observed,
            intercepts,
            category_slopes,
            row_weights=weights,
            person_levels=person_levels,
            method="WLE",
        )
        estimates[draw_number, :] = scored["Estimate"].to_numpy(dtype=float)
        conditional_ses[draw_number, :] = scored["StandardError"].to_numpy(dtype=float)
        statuses[draw_number, :] = scored["Status"].astype(str).to_numpy()
        if retain_draws:
            draw_part = scored[
                ["Person", "Status", "Estimate", "StandardError", "AdjustedScoreResidual"]
            ].copy()
            draw_part.insert(0, "Draw", draw_number + 1)
            draw_parts.append(draw_part)

    successful = statuses == "ok"
    finite = np.isfinite(estimates) & np.isfinite(conditional_ses) & successful
    calibration_sd = np.array(
        [
            float(np.std(estimates[finite[:, column], column], ddof=1))
            if int(finite[:, column].sum()) >= 2 else np.nan
            for column in range(len(person_levels))
        ]
    )
    draw_mean = np.array(
        [
            float(np.mean(estimates[finite[:, column], column]))
            if int(finite[:, column].sum()) else np.nan
            for column in range(len(person_levels))
        ]
    )
    draw_quantiles = np.column_stack(
        [
            np.array(
                [
                    float(np.quantile(estimates[finite[:, column], column], quantile))
                    if int(finite[:, column].sum()) else np.nan
                    for column in range(len(person_levels))
                ]
            )
            for quantile in quantile_values
        ]
    )
    persons = baseline.copy()
    persons["CalibrationDrawsRequested"] = draws_requested
    persons["CalibrationDrawsSuccessful"] = finite.sum(axis=0).astype(int)
    persons["CalibrationDrawMean"] = draw_mean
    persons["CalibrationDrawMeanMinusBaseline"] = draw_mean - persons["Estimate"].to_numpy(dtype=float)
    persons["CalibrationDrawSD"] = calibration_sd
    for index, quantile in enumerate(quantile_values):
        label = f"CalibrationDrawQ{int(round(1000 * quantile)):03d}"
        persons[label] = draw_quantiles[:, index]
    baseline_se = persons["StandardError"].to_numpy(dtype=float)
    persons["QuadratureSensitivitySE"] = np.sqrt(
        np.square(baseline_se) + np.square(calibration_sd)
    )
    persons["CalibrationSDToConditionalSERatio"] = np.divide(
        calibration_sd,
        baseline_se,
        out=np.full_like(calibration_sd, np.nan),
        where=np.isfinite(baseline_se) & (baseline_se > 0),
    )
    baseline_estimate = persons["Estimate"].to_numpy(dtype=float)
    same_sign = np.sign(estimates) == np.sign(baseline_estimate)[None, :]
    persons["CalibrationDrawSameSignShare"] = np.array(
        [
            float(np.mean(same_sign[finite[:, column], column]))
            if int(finite[:, column].sum()) else np.nan
            for column in range(len(person_levels))
        ]
    )
    persons["CalibrationSensitivityScope"] = CALIBRATION_SENSITIVITY_SCOPE
    persons["CalibrationDrawIntervalIsConfidenceInterval"] = False
    persons["QuadratureSensitivityIsInferenceQualified"] = False

    all_successful = bool(finite.all())
    summary = pd.DataFrame(
        [
            {
                "Method": "CMLE_WLE_asymptotic_calibration_draw_sensitivity",
                "Model": design.model,
                "Persons": len(persons),
                "Rows": len(design.data),
                "FreeCalibrationParameters": design.n_parameters,
                "DrawsRequested": draws_requested,
                "Seed": seed_value,
                "CovarianceScale": covariance_scale,
                "RowDrawWork": work,
                "AllDrawPersonScoresSuccessful": all_successful,
                "FailedDrawPersonScores": int((~finite).sum()),
                "MedianCalibrationDrawSD": float(np.nanmedian(calibration_sd)),
                "MaximumCalibrationDrawSD": float(np.nanmax(calibration_sd)),
                "MedianCalibrationSDToConditionalSERatio": float(
                    np.nanmedian(persons["CalibrationSDToConditionalSERatio"])
                ),
                "MaximumCalibrationSDToConditionalSERatio": float(
                    np.nanmax(persons["CalibrationSDToConditionalSERatio"])
                ),
                "MaximumAbsDrawMeanMinusBaseline": float(
                    np.nanmax(np.abs(persons["CalibrationDrawMeanMinusBaseline"]))
                ),
                **covariance_audit,
                "CalibrationDrawIntervalIsConfidenceInterval": False,
                "QuadratureSensitivityIsInferenceQualified": False,
                "Scope": CALIBRATION_SENSITIVITY_SCOPE,
            }
        ]
    )
    return {
        "summary": summary,
        "persons": persons,
        "draws": (
            pd.concat(draw_parts, ignore_index=True)
            if retain_draws else pd.DataFrame()
        ),
    }


__all__ = [
    "CALIBRATION_SENSITIVITY_SCOPE",
    "cmle_wle_calibration_sensitivity",
]
