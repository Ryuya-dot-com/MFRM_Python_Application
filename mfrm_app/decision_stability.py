"""Floating-point and display-rounding audits for threshold decisions.

The statistical decision contract in this module is deliberately small:

* classifications always use the finite, unrounded value;
* display rounding never changes the stored decision;
* values close enough to a threshold to be hidden by display rounding are
  labelled for review instead of being silently presented as unambiguous;
* exact endpoint semantics are explicit and tested.

This is not an uncertainty interval for a fitted statistic.  It is an audit
of numerical and presentation sensitivity around a configured decision rule.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any

import numpy as np
import pandas as pd


FIT_MNSQ_THRESHOLDS: tuple[float, ...] = (0.5, 1.5, 2.0)
FIT_DISPLAY_DECIMALS = 3
BIAS_DISPLAY_DECIMALS = 4
_IDENTIFIER_COLUMNS = (
    "Facet",
    "Level",
    "Person",
    "Element",
    "FacetPair",
    "FacetA",
    "FacetA_Level",
    "FacetB",
    "FacetB_Level",
)


def classify_fit_mnsq(value: Any) -> str:
    """Classify a mean-square fit value using explicit endpoint semantics.

    Exact 0.50 and 1.50 are acceptable, and exact 2.00 is noisy.  Thus the
    bands are ``<0.50``, ``[0.50, 1.50]``, ``(1.50, 2.00]``, and ``>2.00``.
    """

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "unavailable"
    if not np.isfinite(numeric):
        return "unavailable"
    if numeric < 0.5:
        return "overfit"
    if numeric <= 1.5:
        return "acceptable"
    if numeric <= 2.0:
        return "noisy"
    return "distorting"


def fit_classification_rank(classification: Any) -> int:
    """Return a severity rank suitable for selecting a row's worst result."""

    return {
        "unavailable": -1,
        "acceptable": 0,
        "overfit": 1,
        "noisy": 2,
        "distorting": 3,
    }.get(str(classification), -1)


def worst_fit_classification(*values: Any) -> str:
    """Return the worst finite classification among one or more MNSQ values."""

    classes = [classify_fit_mnsq(value) for value in values]
    available = [value for value in classes if value != "unavailable"]
    if not available:
        return "unavailable"
    return max(available, key=fit_classification_rank)


def fit_mnsq_review_mask(values: Any) -> pd.Series:
    """Return ``True`` for finite MNSQ values outside the acceptable band."""

    if isinstance(values, pd.Series):
        series = pd.to_numeric(values, errors="coerce")
    else:
        series = pd.to_numeric(pd.Series(values), errors="coerce")
    classes = series.map(classify_fit_mnsq)
    return classes.isin({"overfit", "noisy", "distorting"})


def numerical_tolerance(value: float, threshold: float, *, ulps: int = 64) -> float:
    """Return a scale-aware tolerance used only to label boundary proximity."""

    scale = max(1.0, abs(float(value)), abs(float(threshold)))
    return float(max(1, int(ulps)) * np.finfo(float).eps * scale)


def _boundary_evaluation(
    value: Any,
    *,
    thresholds: Iterable[float],
    classifier: Callable[[float], str],
    display_decimals: int,
    transform: Callable[[float], float] | None = None,
) -> dict[str, Any]:
    """Evaluate numerical/display proximity without changing the raw decision."""

    decimals = max(0, int(display_decimals))
    try:
        raw_value = float(value)
    except (TypeError, ValueError):
        raw_value = np.nan
    if not np.isfinite(raw_value):
        return {
            "RawValue": raw_value,
            "DecisionValue": np.nan,
            "DisplayValue": np.nan,
            "NearestThreshold": np.nan,
            "DistanceToThreshold": np.nan,
            "NumericalTolerance": np.nan,
            "DisplayHalfUnit": 0.5 * 10.0 ** (-decimals),
            "BoundaryStatus": "unavailable",
            "DecisionStable": False,
            "RawDecision": "unavailable",
            "DisplayDecision": "unavailable",
            "DisplayDecisionConsistent": True,
        }

    decision_value = transform(raw_value) if transform else raw_value
    display_raw = float(np.round(raw_value, decimals=decimals))
    display_decision_value = transform(display_raw) if transform else display_raw
    threshold_values = np.asarray(tuple(float(v) for v in thresholds), dtype=float)
    threshold_values = threshold_values[np.isfinite(threshold_values)]
    if threshold_values.size:
        distances = np.abs(threshold_values - decision_value)
        nearest = float(threshold_values[int(np.argmin(distances))])
        distance = float(np.min(distances))
        num_tol = numerical_tolerance(decision_value, nearest)
    else:
        nearest = distance = num_tol = np.nan
    display_half_unit = 0.5 * 10.0 ** (-decimals)
    if np.isfinite(distance) and distance <= num_tol:
        boundary_status = "numerical_boundary"
    elif np.isfinite(distance) and distance <= display_half_unit + num_tol:
        boundary_status = "display_rounding_boundary"
    else:
        boundary_status = "stable"
    raw_decision = classifier(decision_value)
    display_decision = classifier(display_decision_value)
    return {
        "RawValue": raw_value,
        "DecisionValue": decision_value,
        "DisplayValue": display_raw,
        "NearestThreshold": nearest,
        "DistanceToThreshold": distance,
        "NumericalTolerance": num_tol,
        "DisplayHalfUnit": display_half_unit,
        "BoundaryStatus": boundary_status,
        "DecisionStable": boundary_status == "stable",
        "RawDecision": raw_decision,
        "DisplayDecision": display_decision,
        "DisplayDecisionConsistent": raw_decision == display_decision,
    }


def evaluate_fit_mnsq(value: Any, *, display_decimals: int = FIT_DISPLAY_DECIMALS) -> dict[str, Any]:
    """Audit one MNSQ value against all canonical fit thresholds."""

    result = _boundary_evaluation(
        value,
        thresholds=FIT_MNSQ_THRESHOLDS,
        classifier=classify_fit_mnsq,
        display_decimals=display_decimals,
    )
    result["DecisionRule"] = (
        "unrounded MNSQ: <0.50 overfit; 0.50-1.50 acceptable; "
        ">1.50-2.00 noisy; >2.00 distorting"
    )
    return result


def _zstd_classifier(value: float) -> str:
    return "review" if float(value) >= 2.0 else "acceptable"


def _probability_classifier(alpha: float) -> Callable[[float], str]:
    return lambda value: "significant" if float(value) < float(alpha) else "not_significant"


def _practical_classifier(threshold: float) -> Callable[[float], str]:
    return lambda value: "practical_review" if float(value) >= float(threshold) else "below_practical_threshold"


def _audit_rows(
    frame: pd.DataFrame,
    specifications: Iterable[dict[str, Any]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return pd.DataFrame()
    work = frame.reset_index(drop=False).rename(columns={"index": "SourceIndex"})
    identifiers = [column for column in _IDENTIFIER_COLUMNS if column in work.columns]
    for source_row, row in work.iterrows():
        identity = {column: row.get(column) for column in identifiers}
        for spec in specifications:
            column = str(spec["column"])
            if column not in work.columns:
                continue
            evaluation = _boundary_evaluation(
                row.get(column),
                thresholds=spec["thresholds"],
                classifier=spec["classifier"],
                display_decimals=int(spec["display_decimals"]),
                transform=spec.get("transform"),
            )
            rows.append({
                "SourceRow": int(source_row),
                "SourceIndex": row.get("SourceIndex"),
                **identity,
                "Statistic": column,
                "DecisionRule": str(spec["decision_rule"]),
                **evaluation,
            })
    return pd.DataFrame(rows)


def audit_fit_decision_stability(
    frame: pd.DataFrame,
    *,
    display_decimals: int = FIT_DISPLAY_DECIMALS,
) -> pd.DataFrame:
    """Return row-level fit threshold and display-rounding evidence."""

    mnsq_columns = (
        "Infit",
        "Outfit",
        "InfitMnSq",
        "OutfitMnSq",
        "Infit_MNSQ",
        "Outfit_MNSQ",
    )
    zstd_columns = ("InfitZSTD", "OutfitZSTD", "InfitZStd", "OutfitZStd", "Infit_ZSTD", "Outfit_ZSTD")
    specifications: list[dict[str, Any]] = []
    for column in mnsq_columns:
        specifications.append({
            "column": column,
            "thresholds": FIT_MNSQ_THRESHOLDS,
            "classifier": classify_fit_mnsq,
            "display_decimals": display_decimals,
            "decision_rule": (
                "unrounded MNSQ: <0.50 overfit; 0.50-1.50 acceptable; "
                ">1.50-2.00 noisy; >2.00 distorting"
            ),
        })
    for column in zstd_columns:
        specifications.append({
            "column": column,
            "thresholds": (2.0,),
            "classifier": _zstd_classifier,
            "display_decimals": display_decimals,
            "transform": abs,
            "decision_rule": "unrounded |ZSTD|: >=2.00 review; <2.00 acceptable",
        })
    return _audit_rows(frame, specifications)


def summarize_fit_decision_stability(audit: pd.DataFrame) -> pd.DataFrame:
    """Return a one-row user-facing summary of a fit stability audit."""

    columns = [
        "Status",
        "StatisticsAudited",
        "AvailableStatistics",
        "UnavailableStatistics",
        "StableStatistics",
        "NumericalBoundaryStatistics",
        "DisplayBoundaryStatistics",
        "DisplayDecisionMismatches",
        "DecisionContract",
        "RecommendedAction",
    ]
    if not isinstance(audit, pd.DataFrame) or audit.empty:
        return pd.DataFrame([{
            "Status": "Missing",
            "StatisticsAudited": 0,
            "AvailableStatistics": 0,
            "UnavailableStatistics": 0,
            "StableStatistics": 0,
            "NumericalBoundaryStatistics": 0,
            "DisplayBoundaryStatistics": 0,
            "DisplayDecisionMismatches": 0,
            "DecisionContract": "Decisions require finite unrounded fit statistics.",
            "RecommendedAction": "Compute element fit before applying threshold-based labels.",
        }], columns=columns)
    statuses = audit.get("BoundaryStatus", pd.Series(dtype=str)).astype(str)
    consistent = audit.get("DisplayDecisionConsistent", pd.Series(True, index=audit.index)).fillna(True).astype(bool)
    numerical_n = int(statuses.eq("numerical_boundary").sum())
    display_n = int(statuses.eq("display_rounding_boundary").sum())
    unavailable_n = int(statuses.eq("unavailable").sum())
    available_n = int(len(audit) - unavailable_n)
    mismatch_n = int((~consistent).sum())
    if available_n == 0:
        status = "Missing"
    elif numerical_n + display_n + mismatch_n + unavailable_n > 0:
        status = "Review"
    else:
        status = "Ready"
    return pd.DataFrame([{
        "Status": status,
        "StatisticsAudited": int(len(audit)),
        "AvailableStatistics": available_n,
        "UnavailableStatistics": unavailable_n,
        "StableStatistics": int(statuses.eq("stable").sum()),
        "NumericalBoundaryStatistics": numerical_n,
        "DisplayBoundaryStatistics": display_n,
        "DisplayDecisionMismatches": mismatch_n,
        "DecisionContract": "All classifications use finite unrounded values; displayed values are presentation only.",
        "RecommendedAction": (
            "No finite fit statistics were available; compute fit before applying threshold-based labels."
            if status == "Missing" else
            "Inspect boundary rows and report the raw value, rule, and sensitivity; do not treat a rounded display as the decision input."
            if status == "Review" and unavailable_n == 0 else
            "Inspect boundary rows and unavailable statistics; report raw values and do not treat missing or rounded output as stable evidence."
            if status == "Review" else
            "No fit statistic is close enough to a configured threshold to be hidden by the selected display precision."
        ),
    }], columns=columns)


def audit_bias_decision_stability(
    frame: pd.DataFrame,
    *,
    alpha: float = 0.05,
    practical_logit: float = 0.50,
    display_decimals: int = BIAS_DISPLAY_DECIMALS,
) -> pd.DataFrame:
    """Audit p-value and practical-effect thresholds in a DFF/bias table."""

    specifications = [
        {
            "column": column,
            "thresholds": (float(alpha),),
            "classifier": _probability_classifier(alpha),
            "display_decimals": display_decimals,
            "decision_rule": f"unrounded {column} < {float(alpha):.12g}",
        }
        for column in ("p_holm", "p_bh")
    ]
    specifications.append({
        "column": "AbsBias",
        "thresholds": (float(practical_logit),),
        "classifier": _practical_classifier(practical_logit),
        "display_decimals": display_decimals,
        "decision_rule": f"unrounded |bias| >= {float(practical_logit):.12g}",
    })
    return _audit_rows(frame, specifications)


def summarize_boundary_audit(audit: pd.DataFrame) -> dict[str, int]:
    """Return compact counts shared by bias and other boundary audits."""

    if not isinstance(audit, pd.DataFrame) or audit.empty:
        return {
            "audited": 0,
            "numerical_boundary": 0,
            "display_rounding_boundary": 0,
            "display_decision_mismatch": 0,
        }
    status = audit.get("BoundaryStatus", pd.Series(dtype=str)).astype(str)
    consistent = audit.get("DisplayDecisionConsistent", pd.Series(True, index=audit.index)).fillna(True).astype(bool)
    return {
        "audited": int(len(audit)),
        "numerical_boundary": int(status.eq("numerical_boundary").sum()),
        "display_rounding_boundary": int(status.eq("display_rounding_boundary").sum()),
        "display_decision_mismatch": int((~consistent).sum()),
    }
