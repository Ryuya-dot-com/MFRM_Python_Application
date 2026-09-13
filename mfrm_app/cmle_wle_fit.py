"""Fixed-calibration Person MnSq for the exact-CMLE plus Warm-WLE bridge.

The module is intentionally independent of Streamlit and of bootstrap
orchestration.  It evaluates category moments at a fixed Person measure and
an exact-CMLE structural calibration, then reports untrimmed Person Infit and
Outfit mean squares.  The primary parity target is ``sirt::pcm.fit`` with
identical fixed ``theta`` and cumulative category difficulties ``b_k=-q_k``.

ZSTD, p-values, confidence intervals, and total standard errors are withheld:
their reference distributions and uncertainty composition are outside this
validated phase.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from mfrm_app.cmle_person_scoring import (
    CMLE_WLE_ESTIMATOR_LABEL,
    _ready_cmle_inputs,
    score_cmle_persons_wle,
)
from mfrm_app.decision_stability import (
    FIT_DISPLAY_DECIMALS,
    audit_fit_decision_stability,
    evaluate_fit_mnsq,
    summarize_fit_decision_stability,
    worst_fit_classification,
)


CMLE_WLE_PERSON_FIT_LABEL = "CMLE_WLE_fixed_calibration_person_mnsq"
PERSON_FIT_REFERENCE = "sirt::pcm.fit_untrimmed_fixed_theta"
WITHHELD_REFERENCE_DISTRIBUTION = (
    "withheld: CMLE-WLE finite-sample reference distribution and effective "
    "degrees of freedom are not validated"
)


def _as_1d(values: Sequence[Any], name: str, n_rows: int | None = None) -> np.ndarray:
    out = np.asarray(values).reshape(-1)
    if n_rows is not None and len(out) != int(n_rows):
        raise ValueError(f"{name} length must match the category-kernel rows.")
    return out


def _validate_fixed_fit_inputs(
    person_index: Sequence[object],
    observed_category: Sequence[int],
    category_intercepts: np.ndarray,
    person_estimate: Sequence[float],
    *,
    row_weights: Sequence[float] | None,
    category_available: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    intercepts = np.asarray(category_intercepts, dtype=float)
    if intercepts.ndim != 2 or intercepts.shape[1] < 2:
        raise ValueError("category_intercepts must be a 2D array with at least two categories.")
    n_rows, n_categories = intercepts.shape
    persons = _as_1d(person_index, "person_index", n_rows).astype(object)
    if pd.isna(persons).any():
        raise ValueError("person_index must not contain missing values.")
    persons = np.asarray([str(value) for value in persons], dtype=object)
    observed_raw = _as_1d(observed_category, "observed_category", n_rows)
    try:
        observed_numeric = np.asarray(observed_raw, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("observed_category must contain finite integer indices.") from exc
    if not np.all(np.isfinite(observed_numeric)) or not np.all(
        observed_numeric == np.floor(observed_numeric)
    ):
        raise ValueError("observed_category must contain finite integer indices.")
    observed = observed_numeric.astype(int)
    if np.any((observed < 0) | (observed >= n_categories)):
        raise ValueError("observed_category contains an index outside the category arrays.")
    theta = np.asarray(
        _as_1d(person_estimate, "person_estimate", n_rows), dtype=float
    )
    if not np.all(np.isfinite(theta)):
        raise ValueError("Every row-aligned person_estimate must be finite.")
    weights = (
        np.ones(n_rows, dtype=float)
        if row_weights is None
        else np.asarray(_as_1d(row_weights, "row_weights", n_rows), dtype=float)
    )
    if not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError("Every row weight must be finite and strictly positive.")
    available = (
        np.isfinite(intercepts)
        if category_available is None
        else np.asarray(category_available, dtype=bool)
    )
    if available.shape != intercepts.shape:
        raise ValueError("category_available must match category_intercepts.")
    if np.any(available & ~np.isfinite(intercepts)):
        raise ValueError("Every available category intercept must be finite.")
    if np.any(available.sum(axis=1) < 2):
        raise ValueError("Every row must have at least two available categories.")
    if np.any(~available[np.arange(n_rows), observed]):
        raise ValueError("An observed category is unavailable under its row support.")
    return persons, observed, intercepts, theta, weights, available


def _fixed_category_moments(
    intercepts: np.ndarray,
    theta: np.ndarray,
    available: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    scores = np.arange(intercepts.shape[1], dtype=float)
    logits = intercepts + theta[:, None] * scores[None, :]
    logits = np.where(available, logits, -np.inf)
    row_max = np.max(logits, axis=1, keepdims=True)
    with np.errstate(over="ignore", under="ignore", invalid="ignore"):
        exp_logits = np.where(available, np.exp(logits - row_max), 0.0)
    denominators = exp_logits.sum(axis=1, keepdims=True)
    if np.any(~np.isfinite(denominators) | (denominators <= 0)):
        raise FloatingPointError("Category probabilities could not be normalized.")
    probabilities = exp_logits / denominators
    expected = probabilities @ scores
    centered = scores[None, :] - expected[:, None]
    variance = np.sum(probabilities * centered**2, axis=1)
    fourth = np.sum(probabilities * centered**4, axis=1)
    return probabilities, expected, variance, fourth


def compute_fixed_calibration_person_fit(
    person_index: Sequence[object],
    observed_category: Sequence[int],
    category_intercepts: np.ndarray,
    person_estimate: Sequence[float],
    *,
    row_weights: Sequence[float] | None = None,
    category_available: np.ndarray | None = None,
    person_levels: Sequence[object] | None = None,
    display_decimals: int = FIT_DISPLAY_DECIMALS,
) -> dict[str, object]:
    """Compute untrimmed fixed-theta Person Infit and Outfit MnSq.

    ``category_intercepts`` contains the category log kernels without the
    Person term.  ``person_estimate`` is row-aligned.  Missing responses must
    be removed by the caller; they are never converted to category zero.
    """

    decimals = int(display_decimals)
    if decimals < 0:
        raise ValueError("display_decimals must be non-negative.")
    persons, observed, intercepts, theta, weights, available = (
        _validate_fixed_fit_inputs(
            person_index,
            observed_category,
            category_intercepts,
            person_estimate,
            row_weights=row_weights,
            category_available=category_available,
        )
    )
    probabilities, expected, variance, fourth = _fixed_category_moments(
        intercepts, theta, available
    )
    row_id = np.arange(len(persons))
    residual = observed.astype(float) - expected
    observed_probability = probabilities[row_id, observed]
    usable_variance = np.isfinite(variance) & (variance > 0.0)
    standardized_residual = np.full(len(persons), np.nan, dtype=float)
    standardized_squared = np.full(len(persons), np.nan, dtype=float)
    standardized_residual[usable_variance] = (
        residual[usable_variance] / np.sqrt(variance[usable_variance])
    )
    standardized_squared[usable_variance] = (
        residual[usable_variance] ** 2 / variance[usable_variance]
    )
    observation = pd.DataFrame(
        {
            "ObservationIndex": row_id,
            "Person": persons,
            "ObservedInternalCategory": observed,
            "PersonEstimate": theta,
            "ExpectedInternalCategory": expected,
            "Variance": variance,
            "FourthCentralMoment": fourth,
            "Residual": residual,
            "StandardizedResidual": standardized_residual,
            "StandardizedSquaredResidual": standardized_squared,
            "ObservedProbability": observed_probability,
            "RowWeight": weights,
            "FitContributionAvailable": usable_variance
            & np.isfinite(residual)
            & np.isfinite(fourth)
            & np.isfinite(observed_probability),
        }
    )

    observed_levels = list(dict.fromkeys(persons.tolist()))
    if person_levels is None:
        levels = observed_levels
    else:
        levels = [str(value) for value in person_levels]
        if len(levels) != len(set(levels)) or set(levels) != set(observed_levels):
            raise ValueError("person_levels must identify every retained Person exactly once.")

    person_rows: list[dict[str, object]] = []
    for person in levels:
        group = observation.loc[observation["Person"].eq(person)]
        available_fit = bool(group["FitContributionAvailable"].all()) and len(group) > 0
        squared_residual_weighted = float(
            np.sum(group["RowWeight"] * group["Residual"] ** 2)
        )
        variance_weighted = float(np.sum(group["RowWeight"] * group["Variance"]))
        outfit_numerator = float(
            np.sum(group["RowWeight"] * group["StandardizedSquaredResidual"])
        )
        outfit_denominator = float(group["RowWeight"].sum())
        available_fit = bool(
            available_fit
            and np.isfinite(squared_residual_weighted)
            and np.isfinite(variance_weighted)
            and variance_weighted > 0.0
            and np.isfinite(outfit_numerator)
            and np.isfinite(outfit_denominator)
            and outfit_denominator > 0.0
        )
        infit = squared_residual_weighted / variance_weighted if available_fit else np.nan
        outfit = outfit_numerator / outfit_denominator if available_fit else np.nan
        infit_eval = evaluate_fit_mnsq(infit, display_decimals=decimals)
        outfit_eval = evaluate_fit_mnsq(outfit, display_decimals=decimals)
        person_rows.append(
            {
                "Person": person,
                "NObservations": int(len(group)),
                "ObservationWeight": outfit_denominator,
                "InfitNumerator": squared_residual_weighted,
                "InfitDenominator": variance_weighted,
                "OutfitNumerator": outfit_numerator,
                "OutfitDenominator": outfit_denominator,
                "Infit": float(infit),
                "Outfit": float(outfit),
                "InfitClass": infit_eval["RawDecision"],
                "OutfitClass": outfit_eval["RawDecision"],
                "WorstFitClass": worst_fit_classification(infit, outfit),
                "InfitDisplay": infit_eval["DisplayValue"],
                "OutfitDisplay": outfit_eval["DisplayValue"],
                "InfitBoundaryStatus": infit_eval["BoundaryStatus"],
                "OutfitBoundaryStatus": outfit_eval["BoundaryStatus"],
                "InfitDecisionStable": bool(infit_eval["DecisionStable"]),
                "OutfitDecisionStable": bool(outfit_eval["DecisionStable"]),
                "InfitDisplayDecisionConsistent": bool(
                    infit_eval["DisplayDecisionConsistent"]
                ),
                "OutfitDisplayDecisionConsistent": bool(
                    outfit_eval["DisplayDecisionConsistent"]
                ),
                "PersonFitReady": available_fit,
                "PersonFitStatus": (
                    "ready_untrimmed_fixed_theta"
                    if available_fit
                    else "unavailable_nonfinite_or_nonpositive_variance"
                ),
                "DecisionInput": "finite_unrounded_mnsq",
                "DisplayDecimals": decimals,
                "OutlierTrimmingApplied": False,
                "InfitZSTD": np.nan,
                "OutfitZSTD": np.nan,
                "ZSTDStatus": WITHHELD_REFERENCE_DISTRIBUTION,
                "FitPValue": np.nan,
                "PValueStatus": WITHHELD_REFERENCE_DISTRIBUTION,
            }
        )
    person_frame = pd.DataFrame(person_rows)
    audit = audit_fit_decision_stability(person_frame, display_decimals=decimals)
    audit_summary = summarize_fit_decision_stability(audit)
    return {
        "observations": observation,
        "persons": person_frame,
        "decision_audit": audit,
        "decision_summary": audit_summary,
        "contract": {
            "EstimatorLabel": CMLE_WLE_PERSON_FIT_LABEL,
            "ReferenceImplementation": PERSON_FIT_REFERENCE,
            "PersonMeasureRole": "fixed_calibration_wle",
            "OutlierTrimmingApplied": False,
            "DecisionsUseUnroundedValues": True,
            "ZSTDReady": False,
            "PValueReady": False,
            "ConfidenceIntervalReady": False,
            "TotalStandardErrorReady": False,
        },
    }


def compute_cmle_wle_person_fit(
    cmle_result: dict[str, object],
    *,
    display_decimals: int = FIT_DISPLAY_DECIMALS,
) -> dict[str, object]:
    """Compute fit-sample Person MnSq from an inference-ready CMLE-WLE fit."""

    design, parameters = _ready_cmle_inputs(cmle_result)
    wle = score_cmle_persons_wle(cmle_result).copy()
    if not bool(wle["WLEPersonScoringReady"].all()):
        raise ValueError("Every retained Person must have a ready fixed-calibration WLE.")
    estimates = wle.assign(Person=wle["Person"].astype(str)).set_index("Person")[
        "Estimate"
    ]
    persons = design.data[design.person_col].astype(str).to_numpy(dtype=object)
    row_theta = pd.Series(persons).map(estimates).to_numpy(dtype=float)
    if not np.isfinite(row_theta).all():
        raise ValueError("CMLE design rows could not be matched to finite WLE Person measures.")
    category_intercepts = design.row_offset + np.einsum(
        "rkp,p->rk", design.row_design, parameters, optimize=True
    )
    observed = (
        pd.to_numeric(design.data[design.score_col], errors="raise").to_numpy(dtype=int)
        - int(design.rating_min)
    )
    weights = (
        pd.to_numeric(design.data["__cmle_weight__"], errors="raise").to_numpy(
            dtype=float
        )
        if "__cmle_weight__" in design.data.columns
        else np.ones(len(design.data), dtype=float)
    )
    levels = cmle_result["person_status"]["Person"].astype(str).tolist()
    result = compute_fixed_calibration_person_fit(
        persons,
        observed,
        category_intercepts,
        row_theta,
        row_weights=weights,
        person_levels=levels,
        display_decimals=display_decimals,
    )
    observations = result["observations"].copy()
    observations["ObservedCategory"] = observed + int(design.rating_min)
    observations["ExpectedCategory"] = (
        observations["ExpectedInternalCategory"] + int(design.rating_min)
    )
    observations["VirtualUnit"] = [
        "|".join(
            f"{facet}={str(design.data.iloc[row_index][facet])}"
            for facet in design.facet_cols
        )
        for row_index in range(len(design.data))
    ]
    for facet in design.facet_cols:
        observations[facet] = design.data[facet].astype(str).to_numpy()

    persons_out = result["persons"].merge(
        wle[
            [
                "Person",
                "Estimate",
                "StandardError",
                "ExtremeScorePattern",
                "ExtremeScoreDirection",
                "ConditionalStatus",
            ]
        ].rename(
            columns={
                "Estimate": "WLEEstimate",
                "StandardError": "ConditionalWLEStandardError",
            }
        ),
        on="Person",
        how="left",
        validate="one_to_one",
    )
    persons_out["CombinedEstimator"] = CMLE_WLE_ESTIMATOR_LABEL
    persons_out["PersonFitEstimator"] = CMLE_WLE_PERSON_FIT_LABEL
    persons_out["ReferenceImplementation"] = PERSON_FIT_REFERENCE
    persons_out["ScoringPopulation"] = "cmle_fit_sample"
    persons_out["CalibrationUncertaintyPropagated"] = False
    result["observations"] = observations
    result["persons"] = persons_out
    result["contract"] = {
        **result["contract"],
        "StructuralEstimator": "CMLE_exact_conditional",
        "CombinedEstimator": CMLE_WLE_ESTIMATOR_LABEL,
        "Models": ["RSM", "PCM"],
        "ScoringPopulation": "cmle_fit_sample",
        "CalibrationUncertaintyPropagated": False,
    }
    return result


__all__ = [
    "CMLE_WLE_PERSON_FIT_LABEL",
    "PERSON_FIT_REFERENCE",
    "WITHHELD_REFERENCE_DISTRIBUTION",
    "compute_cmle_wle_person_fit",
    "compute_fixed_calibration_person_fit",
]
