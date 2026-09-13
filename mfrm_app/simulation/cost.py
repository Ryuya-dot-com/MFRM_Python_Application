"""Rating-session, score-record, and explicit weighted-cost contracts."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
import string

from .schema import DesignSpecV1, design_spec_fingerprint, normalize_design_spec


DESIGN_WORKLOAD_VERSION = "mfrm_design_workload_v1"
COST_POLICY_VERSION = "mfrm_cost_policy_v1"
WEIGHTED_DESIGN_COST_VERSION = "mfrm_weighted_design_cost_v1"


def _require_nonnegative_int(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _require_finite_nonnegative(
    name: str,
    value: object,
    *,
    strictly_positive: bool = False,
) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite numeric scalar")
    try:
        numeric = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite numeric scalar") from exc
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be finite")
    if strictly_positive and numeric <= 0:
        raise ValueError(f"{name} must be > 0")
    if not strictly_positive and numeric < 0:
        raise ValueError(f"{name} must be >= 0")
    return 0.0 if numeric == 0.0 else numeric


def _require_fingerprint(name: str, value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 16
        or any(character not in string.hexdigits.lower() for character in value)
        or value != value.lower()
    ):
        raise ValueError(f"{name} must be a 16-character lowercase hex fingerprint")
    return value


def _json_safe_dict(value: object) -> dict:
    payload = asdict(value)
    json.dumps(payload, allow_nan=False)
    return payload


@dataclass(frozen=True, slots=True, kw_only=True)
class DesignWorkloadV1:
    schema_version: str
    design_spec_fingerprint: str
    score_records_per_session: int
    primary_rating_sessions: int
    bridge_rating_sessions: int
    common_anchor_rating_sessions: int
    total_rating_sessions: int
    primary_score_records: int
    bridge_score_records: int
    common_anchor_score_records: int
    total_score_records: int
    calibration_events: int

    def __post_init__(self) -> None:
        if self.schema_version != DESIGN_WORKLOAD_VERSION:
            raise ValueError(
                f"Unsupported design workload version: {self.schema_version!r}"
            )
        _require_fingerprint(
            "design_spec_fingerprint",
            self.design_spec_fingerprint,
        )
        if (
            isinstance(self.score_records_per_session, bool)
            or not isinstance(self.score_records_per_session, int)
            or self.score_records_per_session < 1
        ):
            raise ValueError("score_records_per_session must be a positive integer")
        count_fields = (
            "primary_rating_sessions",
            "bridge_rating_sessions",
            "common_anchor_rating_sessions",
            "total_rating_sessions",
            "primary_score_records",
            "bridge_score_records",
            "common_anchor_score_records",
            "total_score_records",
            "calibration_events",
        )
        for name in count_fields:
            _require_nonnegative_int(name, getattr(self, name))
        if self.primary_rating_sessions == 0:
            raise ValueError("primary_rating_sessions must be positive")
        if self.total_rating_sessions != (
            self.primary_rating_sessions
            + self.bridge_rating_sessions
            + self.common_anchor_rating_sessions
        ):
            raise ValueError("total_rating_sessions must equal its component sum")
        for prefix in ("primary", "bridge", "common_anchor"):
            if getattr(self, f"{prefix}_score_records") != (
                getattr(self, f"{prefix}_rating_sessions")
                * self.score_records_per_session
            ):
                raise ValueError(
                    f"{prefix}_score_records must equal rating sessions times "
                    "score_records_per_session"
                )
        if self.total_score_records != (
            self.primary_score_records
            + self.bridge_score_records
            + self.common_anchor_score_records
        ):
            raise ValueError("total_score_records must equal its component sum")

    def to_dict(self) -> dict:
        return _json_safe_dict(self)


@dataclass(frozen=True, slots=True, kw_only=True)
class CostPolicyV1:
    """Explicit stacked unit prices for the three workload currencies.

    ``rating_session_unit_cost`` prices the artifact-by-rater act.
    ``score_record_unit_cost`` separately prices every stored/scored long row,
    including the first row in a session. The latter is zero by default, so
    long-format expansion never increases the default comparison currency.
    """

    rating_session_unit_cost: float = 1.0
    score_record_unit_cost: float = 0.0
    calibration_event_unit_cost: float = 0.0
    cost_unit: str = "rating-session-equivalent"
    schema_version: str = COST_POLICY_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != COST_POLICY_VERSION:
            raise ValueError(f"Unsupported cost policy version: {self.schema_version!r}")
        object.__setattr__(
            self,
            "rating_session_unit_cost",
            _require_finite_nonnegative(
                "rating_session_unit_cost",
                self.rating_session_unit_cost,
            ),
        )
        object.__setattr__(
            self,
            "score_record_unit_cost",
            _require_finite_nonnegative(
                "score_record_unit_cost",
                self.score_record_unit_cost,
            ),
        )
        object.__setattr__(
            self,
            "calibration_event_unit_cost",
            _require_finite_nonnegative(
                "calibration_event_unit_cost",
                self.calibration_event_unit_cost,
            ),
        )
        if not isinstance(self.cost_unit, str) or not self.cost_unit.strip():
            raise ValueError("cost_unit must be a non-empty string")
        object.__setattr__(self, "cost_unit", self.cost_unit.strip())
        if not any(
            (
                self.rating_session_unit_cost,
                self.score_record_unit_cost,
                self.calibration_event_unit_cost,
            )
        ):
            raise ValueError("at least one workload unit cost must be > 0")

    def to_dict(self) -> dict:
        return _json_safe_dict(self)


_COST_POLICY_FIELDS = tuple(field.name for field in fields(CostPolicyV1))
_COST_POLICY_FIELD_SET = frozenset(_COST_POLICY_FIELDS)


def normalize_cost_policy(value: CostPolicyV1 | Mapping) -> CostPolicyV1:
    """Strictly parse a saved policy without filling omitted identity fields."""
    if isinstance(value, CostPolicyV1):
        return value
    if not isinstance(value, Mapping):
        raise TypeError("cost policy must be a CostPolicyV1 or mapping")
    supplied = set(value)
    if any(not isinstance(name, str) for name in supplied):
        raise ValueError("Saved cost policy field names must be strings")
    missing = _COST_POLICY_FIELD_SET - supplied
    unknown = supplied - _COST_POLICY_FIELD_SET
    if missing:
        raise ValueError(f"Saved cost policy is missing fields: {sorted(missing)}")
    if unknown:
        raise ValueError(f"Saved cost policy contains unknown fields: {sorted(unknown)}")
    return CostPolicyV1(**{name: value[name] for name in _COST_POLICY_FIELDS})


def cost_policy_to_dict(value: CostPolicyV1 | Mapping) -> dict:
    """Return a canonical JSON-primitive policy representation."""
    return normalize_cost_policy(value).to_dict()


def cost_policy_fingerprint(
    value: CostPolicyV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    """Return a stable SHA-256 prefix for a canonical cost policy."""
    if (
        isinstance(length, bool)
        or not isinstance(length, int)
        or not 8 <= length <= 64
    ):
        raise ValueError("fingerprint length must be an integer between 8 and 64")
    encoded = json.dumps(
        cost_policy_to_dict(value),
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:length]


@dataclass(frozen=True, slots=True, kw_only=True)
class WeightedDesignCostV1:
    schema_version: str
    design_spec_fingerprint: str
    cost_policy_schema_version: str
    cost_policy_fingerprint: str
    cost_unit: str
    rating_session_unit_cost: float
    score_record_unit_cost: float
    calibration_event_unit_cost: float
    total_rating_sessions: int
    total_score_records: int
    calibration_events: int
    rating_session_cost: float
    score_record_cost: float
    calibration_event_cost: float
    total_cost: float

    def __post_init__(self) -> None:
        if self.schema_version != WEIGHTED_DESIGN_COST_VERSION:
            raise ValueError(
                f"Unsupported weighted design cost version: {self.schema_version!r}"
            )
        if self.cost_policy_schema_version != COST_POLICY_VERSION:
            raise ValueError(
                "Unsupported embedded cost policy version: "
                f"{self.cost_policy_schema_version!r}"
            )
        _require_fingerprint(
            "design_spec_fingerprint",
            self.design_spec_fingerprint,
        )
        _require_fingerprint(
            "cost_policy_fingerprint",
            self.cost_policy_fingerprint,
        )
        if not isinstance(self.cost_unit, str) or not self.cost_unit.strip():
            raise ValueError("cost_unit must be a non-empty string")
        if self.cost_unit != self.cost_unit.strip():
            raise ValueError("cost_unit must be canonical without outer whitespace")
        for name, strictly_positive in (
            ("rating_session_unit_cost", False),
            ("score_record_unit_cost", False),
            ("calibration_event_unit_cost", False),
            ("rating_session_cost", False),
            ("score_record_cost", False),
            ("calibration_event_cost", False),
            ("total_cost", False),
        ):
            numeric = _require_finite_nonnegative(
                name,
                getattr(self, name),
                strictly_positive=strictly_positive,
            )
            object.__setattr__(self, name, numeric)
        for name in (
            "total_rating_sessions",
            "total_score_records",
            "calibration_events",
        ):
            _require_nonnegative_int(name, getattr(self, name))
        embedded_policy = CostPolicyV1(
            rating_session_unit_cost=self.rating_session_unit_cost,
            score_record_unit_cost=self.score_record_unit_cost,
            calibration_event_unit_cost=self.calibration_event_unit_cost,
            cost_unit=self.cost_unit,
            schema_version=self.cost_policy_schema_version,
        )
        if cost_policy_fingerprint(embedded_policy) != self.cost_policy_fingerprint:
            raise ValueError(
                "cost_policy_fingerprint does not match the embedded unit prices"
            )
        expected_components = (
            (
                "rating_session_cost",
                _finite_cost_product(
                    self.total_rating_sessions,
                    self.rating_session_unit_cost,
                    name="rating_session_cost",
                ),
            ),
            (
                "score_record_cost",
                _finite_cost_product(
                    self.total_score_records,
                    self.score_record_unit_cost,
                    name="score_record_cost",
                ),
            ),
            (
                "calibration_event_cost",
                _finite_cost_product(
                    self.calibration_events,
                    self.calibration_event_unit_cost,
                    name="calibration_event_cost",
                ),
            ),
        )
        for name, expected in expected_components:
            if getattr(self, name) != expected:
                raise ValueError(f"{name} does not match count times unit cost")
        expected_total = (
            self.rating_session_cost
            + self.score_record_cost
            + self.calibration_event_cost
        )
        if self.total_cost != expected_total:
            raise ValueError("total_cost must equal its component sum")

    def to_dict(self) -> dict:
        return _json_safe_dict(self)


def compute_design_workload(
    value: DesignSpecV1 | Mapping,
) -> DesignWorkloadV1:
    """Compute scheduled workload without materializing assignment rows."""
    spec = normalize_design_spec(value)
    primary_sessions = (
        spec.n_persons
        * spec.n_artifacts_per_person
        * spec.raters_per_artifact
    )
    bridge_sessions = (
        spec.n_bridge_artifacts
        * spec.bridge_extra_raters_per_artifact
    )
    common_anchor_sessions = (
        spec.n_common_anchor_artifacts
        * spec.n_groups
        * spec.common_anchor_raters_per_group
    )
    total_sessions = primary_sessions + bridge_sessions + common_anchor_sessions
    score_multiplier = spec.score_records_per_session
    primary_records = primary_sessions * score_multiplier
    bridge_records = bridge_sessions * score_multiplier
    common_anchor_records = common_anchor_sessions * score_multiplier
    return DesignWorkloadV1(
        schema_version=DESIGN_WORKLOAD_VERSION,
        design_spec_fingerprint=design_spec_fingerprint(spec),
        score_records_per_session=score_multiplier,
        primary_rating_sessions=primary_sessions,
        bridge_rating_sessions=bridge_sessions,
        common_anchor_rating_sessions=common_anchor_sessions,
        total_rating_sessions=total_sessions,
        primary_score_records=primary_records,
        bridge_score_records=bridge_records,
        common_anchor_score_records=common_anchor_records,
        total_score_records=(
            primary_records + bridge_records + common_anchor_records
        ),
        calibration_events=spec.calibration_events,
    )


def _finite_cost_product(count: int, unit_cost: float, *, name: str) -> float:
    try:
        value = count * unit_cost
    except OverflowError as exc:
        raise ValueError(f"{name} exceeds the finite weighted-cost range") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} exceeds the finite weighted-cost range")
    return float(value)


def apply_cost_policy(
    workload: DesignWorkloadV1,
    policy: CostPolicyV1 | None = None,
) -> WeightedDesignCostV1:
    """Apply explicit unit prices without silently combining currencies."""
    if not isinstance(workload, DesignWorkloadV1):
        raise TypeError("workload must be a DesignWorkloadV1")
    resolved = policy if policy is not None else CostPolicyV1()
    if not isinstance(resolved, CostPolicyV1):
        raise TypeError("policy must be a CostPolicyV1 or None")
    session_cost = _finite_cost_product(
        workload.total_rating_sessions,
        resolved.rating_session_unit_cost,
        name="rating_session_cost",
    )
    record_cost = _finite_cost_product(
        workload.total_score_records,
        resolved.score_record_unit_cost,
        name="score_record_cost",
    )
    calibration_cost = _finite_cost_product(
        workload.calibration_events,
        resolved.calibration_event_unit_cost,
        name="calibration_event_cost",
    )
    total_cost = session_cost + record_cost + calibration_cost
    if not math.isfinite(total_cost):
        raise ValueError("total_cost exceeds the finite weighted-cost range")
    return WeightedDesignCostV1(
        schema_version=WEIGHTED_DESIGN_COST_VERSION,
        design_spec_fingerprint=workload.design_spec_fingerprint,
        cost_policy_schema_version=resolved.schema_version,
        cost_policy_fingerprint=cost_policy_fingerprint(resolved),
        cost_unit=resolved.cost_unit,
        rating_session_unit_cost=resolved.rating_session_unit_cost,
        score_record_unit_cost=resolved.score_record_unit_cost,
        calibration_event_unit_cost=resolved.calibration_event_unit_cost,
        total_rating_sessions=workload.total_rating_sessions,
        total_score_records=workload.total_score_records,
        calibration_events=workload.calibration_events,
        rating_session_cost=session_cost,
        score_record_cost=record_cost,
        calibration_event_cost=calibration_cost,
        total_cost=float(total_cost),
    )


__all__ = [
    "COST_POLICY_VERSION",
    "DESIGN_WORKLOAD_VERSION",
    "WEIGHTED_DESIGN_COST_VERSION",
    "CostPolicyV1",
    "DesignWorkloadV1",
    "WeightedDesignCostV1",
    "apply_cost_policy",
    "compute_design_workload",
    "cost_policy_fingerprint",
    "cost_policy_to_dict",
    "normalize_cost_policy",
]
