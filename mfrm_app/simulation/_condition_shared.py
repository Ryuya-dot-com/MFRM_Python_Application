"""Shared validation and identity helpers for scientific conditions.

This private module is deliberately standard-library only.  Leaf condition
modules import it directly so the public ``conditions`` compatibility facade
never becomes an internal dependency.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields
import hashlib
import json
import math


MODEL_RSM = "RSM"

FORMAL_MONTE_CARLO_MIN_REPLICATES = 20
MAX_MONTE_CARLO_REPLICATES = 1_000_000
MAX_ESTIMATOR_ITERATIONS = 10_000
MAX_ABS_ADJACENT_THRESHOLD = 100.0
MAX_THRESHOLD_SPAN = 200.0
MAX_GENERATING_STANDARD_DEVIATION = 10.0
MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN = 10.0
MAX_RATING_MIN_ABS = 2**31 - 1
MAX_SIMULATION_PERSON_INDEX = 999_999
MAX_SIMULATION_RATER_INDEX = 999_999
MAX_SIMULATION_ARTIFACT_INDEX = 999
MAX_SIMULATION_CRITERION_INDEX = 999


class SimulationConditionValidationError(ValueError):
    """Raised when a scientific simulation condition is ambiguous."""


def _strict_mapping(
    value: Mapping,
    cls: type,
    *,
    label: str,
) -> dict:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a {cls.__name__} or mapping")
    expected = tuple(field.name for field in fields(cls))
    expected_set = frozenset(expected)
    supplied = set(value)
    if any(not isinstance(name, str) for name in supplied):
        raise SimulationConditionValidationError(
            f"Saved {label} field names must be strings"
        )
    missing = expected_set - supplied
    unknown = supplied - expected_set
    if missing:
        raise SimulationConditionValidationError(
            f"Saved {label} is missing fields: {sorted(missing)}"
        )
    if unknown:
        raise SimulationConditionValidationError(
            f"Saved {label} contains unknown fields: {sorted(unknown)}"
        )
    return {name: value[name] for name in expected}


def _tuple_payload_field(payload: Mapping, name: str, *, label: str) -> tuple:
    """Normalize one saved JSON array without leaking container TypeErrors."""
    value = payload[name]
    if not isinstance(value, (list, tuple)):
        raise SimulationConditionValidationError(
            f"Saved {label} field {name!r} must be an array"
        )
    return tuple(value)


def _require_int(
    name: str,
    value: object,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise SimulationConditionValidationError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise SimulationConditionValidationError(f"{name} must be >= {minimum}")
    if maximum is not None and value > maximum:
        raise SimulationConditionValidationError(f"{name} must be <= {maximum}")
    return value


def _require_real(
    name: str,
    value: object,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    strictly_positive: bool = False,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SimulationConditionValidationError(f"{name} must be numeric")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise SimulationConditionValidationError(f"{name} must be finite") from exc
    if not math.isfinite(result):
        raise SimulationConditionValidationError(f"{name} must be finite")
    if strictly_positive and result <= 0:
        raise SimulationConditionValidationError(f"{name} must be positive")
    if minimum is not None and result < minimum:
        raise SimulationConditionValidationError(f"{name} must be >= {minimum}")
    if maximum is not None and result > maximum:
        raise SimulationConditionValidationError(f"{name} must be <= {maximum}")
    return result


def _require_literal(name: str, value: object, expected: str) -> None:
    if value != expected:
        raise SimulationConditionValidationError(
            f"{name} must remain {expected!r} in the v1 contract"
        )


def _canonical_json_value(value: object) -> object:
    """Canonicalize JSON primitives, including semantically equal signed zero."""
    if isinstance(value, float):
        return 0.0 if value == 0.0 else value
    if isinstance(value, dict):
        return {
            key: _canonical_json_value(child)
            for key, child in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_canonical_json_value(child) for child in value]
    return value


def _json_safe(payload: dict) -> dict:
    normalized = _canonical_json_value(payload)
    if not isinstance(normalized, dict):  # Internal contract: payload is a dict.
        raise TypeError("canonical JSON payload must remain a dictionary")
    json.dumps(normalized, allow_nan=False, ensure_ascii=False)
    return normalized


def _fingerprint(payload: Mapping, *, length: int = 16) -> str:
    if (
        isinstance(length, bool)
        or not isinstance(length, int)
        or not 8 <= length <= 64
    ):
        raise ValueError("fingerprint length must be an integer between 8 and 64")
    encoded = json.dumps(
        dict(payload),
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:length]


__all__ = [
    "FORMAL_MONTE_CARLO_MIN_REPLICATES",
    "MAX_ABS_ADJACENT_THRESHOLD",
    "MAX_CENTRAL_TENDENCY_LOG_SCALE_MEAN",
    "MAX_ESTIMATOR_ITERATIONS",
    "MAX_GENERATING_STANDARD_DEVIATION",
    "MAX_MONTE_CARLO_REPLICATES",
    "MAX_RATING_MIN_ABS",
    "MAX_SIMULATION_ARTIFACT_INDEX",
    "MAX_SIMULATION_CRITERION_INDEX",
    "MAX_SIMULATION_PERSON_INDEX",
    "MAX_SIMULATION_RATER_INDEX",
    "MAX_THRESHOLD_SPAN",
    "MODEL_RSM",
    "SimulationConditionValidationError",
]
