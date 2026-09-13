"""Repository-only CMLE structural calibration to Warm-WLE Person bridge.

The two estimator stages remain explicit: exact conditional maximum
likelihood estimates the supported RSM/PCM structural calibration, then Warm
weighted likelihood scores Persons while treating that calibration as fixed.
The reported Person standard error therefore does not propagate CMLE
calibration uncertainty.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mfrm_app.cmle import CMLEDesign
from mfrm_app.person_scoring import score_fixed_calibration_persons


CMLE_WLE_ESTIMATOR_LABEL = (
    "CMLE_exact_structural_plus_WLE_fixed_calibration_person"
)


def _ready_cmle_inputs(cmle_result: dict[str, object]) -> tuple[CMLEDesign, np.ndarray]:
    if not isinstance(cmle_result, dict):
        raise ValueError("cmle_result must be a fit_cmle result dictionary.")
    summary = cmle_result.get("summary")
    config = cmle_result.get("config")
    design = cmle_result.get("design")
    coefficients = cmle_result.get("coefficients")
    if not isinstance(summary, pd.DataFrame) or len(summary) != 1:
        raise ValueError("CMLE summary must contain exactly one fit row.")
    if not isinstance(config, dict) or str(config.get("method")) != "CMLE":
        raise ValueError("Person bridge requires an explicitly labelled CMLE fit.")
    if not bool(summary.iloc[0].get("InferenceReady", False)):
        raise ValueError("CMLE calibration must be inference-ready before Person scoring.")
    if not isinstance(design, CMLEDesign):
        raise ValueError("CMLE design is unavailable or has an unexpected type.")
    if design.model not in {"RSM", "PCM"}:
        raise ValueError("CMLE-to-WLE bridge supports RSM and PCM calibration only.")
    if not isinstance(coefficients, pd.DataFrame):
        raise ValueError("CMLE free-coefficient table is unavailable.")
    required = {"Parameter", "Estimate"}
    if not required.issubset(coefficients.columns):
        raise ValueError("CMLE coefficient table lacks Parameter or Estimate.")
    names = coefficients["Parameter"].astype(str)
    if names.duplicated().any() or set(names) != set(design.parameter_names):
        raise ValueError("CMLE coefficient identity differs from the fitted design.")
    lookup = coefficients.assign(Parameter=names).set_index("Parameter")["Estimate"]
    parameters = pd.to_numeric(
        lookup.reindex(design.parameter_names), errors="coerce"
    ).to_numpy(dtype=float)
    if parameters.shape != (design.n_parameters,) or not np.isfinite(parameters).all():
        raise ValueError("CMLE coefficients are incomplete or non-finite.")
    if design.row_design.shape != (
        len(design.data),
        design.n_categories,
        design.n_parameters,
    ):
        raise ValueError("CMLE row design is stale or dimensionally inconsistent.")
    if design.row_offset.shape != (len(design.data), design.n_categories):
        raise ValueError("CMLE row offset is stale or dimensionally inconsistent.")
    if not np.isfinite(design.row_offset).all():
        raise ValueError("CMLE row offset contains non-finite values.")
    return design, parameters


def score_cmle_persons_wle(cmle_result: dict[str, object]) -> pd.DataFrame:
    """Score fit-sample Persons from an inference-ready exact CMLE calibration.

    This bridge is intentionally limited to rows retained inside the fitted
    ``CMLEDesign``. New Persons and unseen virtual response units require a
    separate prediction-design contract.
    """
    design, parameters = _ready_cmle_inputs(cmle_result)
    category_intercepts = design.row_offset + np.einsum(
        "rkp,p->rk", design.row_design, parameters, optimize=True
    )
    category_slopes = np.tile(
        np.arange(design.n_categories, dtype=float),
        (len(design.data), 1),
    )
    observed = (
        pd.to_numeric(design.data[design.score_col], errors="raise").to_numpy(dtype=int)
        - int(design.rating_min)
    )
    persons = design.data[design.person_col].astype(str).to_numpy(dtype=object)
    person_status = cmle_result.get("person_status")
    if not isinstance(person_status, pd.DataFrame) or "Person" not in person_status.columns:
        raise ValueError("CMLE Person-status table is unavailable.")
    levels = person_status["Person"].astype(str).tolist()
    if len(levels) != len(set(levels)) or set(levels) != set(persons):
        raise ValueError("CMLE Person identity differs between design and status table.")
    weights = (
        pd.to_numeric(design.data["__cmle_weight__"], errors="raise").to_numpy(dtype=float)
        if "__cmle_weight__" in design.data.columns
        else np.ones(len(design.data), dtype=float)
    )
    scored = score_fixed_calibration_persons(
        persons,
        observed,
        category_intercepts,
        category_slopes,
        row_weights=weights,
        person_levels=levels,
        method="WLE",
    )
    status = person_status.copy()
    status["Person"] = status["Person"].astype(str)
    scored["Person"] = scored["Person"].astype(str)
    out = scored.merge(status, on="Person", how="left", validate="one_to_one")
    out["CombinedEstimator"] = CMLE_WLE_ESTIMATOR_LABEL
    out["StructuralEstimator"] = "CMLE_exact_conditional"
    out["PersonEstimator"] = "WLE_fixed_calibration"
    out["ScoringPopulation"] = "cmle_fit_sample"
    out["CalibrationInferenceReady"] = True
    out["CalibrationUncertaintyPropagated"] = False
    out["StandardErrorScope"] = (
        "conditional Person information with fitted CMLE calibration treated as fixed"
    )
    out["WLEPersonScoringReady"] = out["Status"].eq("ok")
    out["ReportableWLEEstimate"] = out["Estimate"].where(
        out["WLEPersonScoringReady"]
    )
    out["EstimateRole"] = "fixed_calibration_wle_not_jmle_mle"
    return out


__all__ = ["CMLE_WLE_ESTIMATOR_LABEL", "score_cmle_persons_wle"]
