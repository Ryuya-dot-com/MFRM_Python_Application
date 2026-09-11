"""Versioned, JSON-safe contracts for prospective rating designs.

The design contract describes scheduled field work. It does not contain
rater-effect distributions, an estimator, simulation replicates, missingness,
or a fitted-data reference; those belong to later simulation-condition and
study contracts. Keeping those layers separate lets the application compare
rating effort without starting JMLE, MML, or Bayesian computation.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, fields
import hashlib
import json


DESIGN_SPEC_VERSION = "mfrm_design_spec_v1"
ASSIGNMENT_ALGORITHM_V1 = "balanced_round_robin_v1"

FULLY_CROSSED = "fully_crossed"
SPIRAL = "spiral"
NESTED_BRIDGE = "nested_bridge"
DISCONNECTED_GROUPS = "disconnected_groups"
DESIGN_PLAN_CHOICES = (
    FULLY_CROSSED,
    SPIRAL,
    NESTED_BRIDGE,
    DISCONNECTED_GROUPS,
)


class DesignValidationError(ValueError):
    """Raised when a prospective design is ambiguous or internally invalid."""


def _require_int(name: str, value: object, *, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise DesignValidationError(f"{name} must be an integer")
    if value < minimum:
        comparator = "positive" if minimum == 1 else f">= {minimum}"
        raise DesignValidationError(f"{name} must be {comparator}")
    return value


@dataclass(frozen=True, slots=True, kw_only=True)
class DesignSpecV1:
    """Scheduled assignment and workload identity for one rating design.

    ``n_raters`` is the number of active raters, not the size of an available
    pool. Every compiled v1 schedule must therefore assign at least one session
    to every declared rater.

    A ``RatingSession`` is one unique artifact-by-rater scoring act. One
    session may create several long-format ``ScoreRecord`` rows, such as one
    row per analytic-rubric criterion. Common anchor artifacts are additional,
    externally held artifacts shared across groups and are not included in
    ``n_persons``. Calibration events are tracked separately and never create
    observational graph edges.

    The assignment algorithm is part of the saved identity. Its exact v1 rules
    are frozen in the design-contract documentation and compiler golden tests.
    """

    plan: str
    n_persons: int
    n_raters: int
    n_artifacts_per_person: int = 1
    score_records_per_session: int = 1
    raters_per_artifact: int = 1
    n_groups: int = 1
    n_bridge_artifacts: int = 0
    bridge_extra_raters_per_artifact: int = 0
    n_common_anchor_artifacts: int = 0
    common_anchor_raters_per_group: int = 0
    calibration_events: int = 0
    assignment_seed: int = 20260721
    assignment_algorithm: str = ASSIGNMENT_ALGORITHM_V1
    schema_version: str = DESIGN_SPEC_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != DESIGN_SPEC_VERSION:
            raise DesignValidationError(
                f"Unsupported design spec version: {self.schema_version!r}"
            )
        if self.assignment_algorithm != ASSIGNMENT_ALGORITHM_V1:
            raise DesignValidationError(
                "Unsupported assignment algorithm: "
                f"{self.assignment_algorithm!r}"
            )
        if not isinstance(self.plan, str) or self.plan not in DESIGN_PLAN_CHOICES:
            choices = ", ".join(DESIGN_PLAN_CHOICES)
            raise DesignValidationError(
                f"Unknown design plan {self.plan!r}; expected one of: {choices}"
            )

        positive_fields = (
            "n_persons",
            "n_raters",
            "n_artifacts_per_person",
            "score_records_per_session",
            "raters_per_artifact",
            "n_groups",
        )
        nonnegative_fields = (
            "n_bridge_artifacts",
            "bridge_extra_raters_per_artifact",
            "n_common_anchor_artifacts",
            "common_anchor_raters_per_group",
            "calibration_events",
            "assignment_seed",
        )
        for name in positive_fields:
            _require_int(name, getattr(self, name), minimum=1)
        for name in nonnegative_fields:
            _require_int(name, getattr(self, name), minimum=0)
        if self.assignment_seed > 2**63 - 1:
            raise DesignValidationError("assignment_seed must be <= 2**63 - 1")
        if self.raters_per_artifact > self.n_raters:
            raise DesignValidationError(
                "raters_per_artifact cannot exceed n_raters"
            )

        if self.plan == FULLY_CROSSED:
            self._validate_fully_crossed()
        elif self.plan == SPIRAL:
            self._validate_spiral()
        elif self.plan == NESTED_BRIDGE:
            self._validate_nested_bridge()
        else:
            self._validate_disconnected_groups()

    def _validate_fully_crossed(self) -> None:
        if self.raters_per_artifact != self.n_raters:
            raise DesignValidationError(
                "fully_crossed requires raters_per_artifact == n_raters"
            )
        self._require_neutral_structure()

    def _validate_spiral(self) -> None:
        if not 1 <= self.raters_per_artifact < self.n_raters:
            raise DesignValidationError(
                "spiral requires 1 <= raters_per_artifact < n_raters"
            )
        available_slots = (
            self.n_persons
            * self.n_artifacts_per_person
            * self.raters_per_artifact
        )
        if available_slots < self.n_raters:
            raise DesignValidationError(
                "spiral has too few primary rating slots to activate every rater"
            )
        self._require_neutral_structure()

    def _validate_nested_bridge(self) -> None:
        if self.raters_per_artifact != 1:
            raise DesignValidationError(
                "nested_bridge requires one primary rater per artifact"
            )
        if self.n_groups != 1:
            raise DesignValidationError("nested_bridge requires n_groups == 1")
        total_artifacts = self.n_persons * self.n_artifacts_per_person
        if self.n_bridge_artifacts > total_artifacts:
            raise DesignValidationError(
                "n_bridge_artifacts cannot exceed the primary artifact count"
            )
        if self.n_bridge_artifacts == 0:
            if self.bridge_extra_raters_per_artifact != 0:
                raise DesignValidationError(
                    "zero bridge artifacts require zero extra bridge raters"
                )
        elif not 1 <= self.bridge_extra_raters_per_artifact < self.n_raters:
            raise DesignValidationError(
                "bridge artifacts require between 1 and n_raters - 1 extra raters"
            )
        # Primary nesting is at Person level, so repeated artifacts do not add
        # another primary-rater activation opportunity.
        activation_slots = self.n_persons + (
            self.n_bridge_artifacts
            * self.bridge_extra_raters_per_artifact
        )
        if activation_slots < self.n_raters:
            raise DesignValidationError(
                "nested_bridge has too few person/bridge slots to activate every rater"
            )
        self._require_neutral_common_anchor()

    def _require_neutral_structure(self) -> None:
        if self.n_groups != 1:
            raise DesignValidationError(f"{self.plan} requires n_groups == 1")
        if self.n_bridge_artifacts or self.bridge_extra_raters_per_artifact:
            raise DesignValidationError(
                f"{self.plan} does not accept bridge-specific fields"
            )
        self._require_neutral_common_anchor()

    def _require_neutral_common_anchor(self) -> None:
        if (
            self.n_common_anchor_artifacts
            or self.common_anchor_raters_per_group
        ):
            raise DesignValidationError(
                f"{self.plan} does not accept common-anchor-specific fields"
            )

    def _validate_disconnected_groups(self) -> None:
        if self.raters_per_artifact != 1:
            raise DesignValidationError(
                "disconnected_groups requires one primary rater per artifact"
            )
        if not 2 <= self.n_groups <= min(self.n_persons, self.n_raters):
            raise DesignValidationError(
                "disconnected_groups requires "
                "2 <= n_groups <= min(n_persons, n_raters)"
            )
        if self.n_bridge_artifacts or self.bridge_extra_raters_per_artifact:
            raise DesignValidationError(
                "disconnected_groups cannot contain observational bridge fields"
            )
        if self.n_common_anchor_artifacts == 0:
            if self.common_anchor_raters_per_group != 0:
                raise DesignValidationError(
                    "zero common anchor artifacts require zero anchor raters per group"
                )
        else:
            raters_in_smallest_group = self.n_raters // self.n_groups
            if not (
                1
                <= self.common_anchor_raters_per_group
                <= raters_in_smallest_group
            ):
                raise DesignValidationError(
                    "common_anchor_raters_per_group exceeds the smallest rater group"
                )

        person_quotient, person_remainder = divmod(
            self.n_persons,
            self.n_groups,
        )
        rater_quotient, rater_remainder = divmod(
            self.n_raters,
            self.n_groups,
        )
        common_anchor_slots = (
            self.n_common_anchor_artifacts
            * self.common_anchor_raters_per_group
        )
        # Balanced group sizes are piecewise constant and change only at the
        # two remainder boundaries. Check those boundaries in O(1), even when
        # a formula-only prospective design declares billions of groups.
        boundary_indices = {
            0,
            person_remainder,
            rater_remainder,
            self.n_groups - 1,
        }
        for group_index in boundary_indices:
            if not 0 <= group_index < self.n_groups:
                continue
            person_count = person_quotient + (
                1 if group_index < person_remainder else 0
            )
            rater_count = rater_quotient + (
                1 if group_index < rater_remainder else 0
            )
            if person_count + common_anchor_slots < rater_count:
                raise DesignValidationError(
                    "disconnected_groups has too few person/common-anchor slots "
                    "to activate every rater in each group"
                )


_DESIGN_SPEC_FIELDS = tuple(field.name for field in fields(DesignSpecV1))
_DESIGN_SPEC_FIELD_SET = frozenset(_DESIGN_SPEC_FIELDS)


def normalize_design_spec(value: DesignSpecV1 | Mapping) -> DesignSpecV1:
    """Strictly parse a saved design without filling omitted identity fields."""
    if isinstance(value, DesignSpecV1):
        return value
    if not isinstance(value, Mapping):
        raise TypeError("design spec must be a DesignSpecV1 or mapping")
    supplied = set(value)
    if any(not isinstance(name, str) for name in supplied):
        raise DesignValidationError("Saved design spec field names must be strings")
    missing = _DESIGN_SPEC_FIELD_SET - supplied
    unknown = supplied - _DESIGN_SPEC_FIELD_SET
    if missing:
        raise DesignValidationError(
            f"Saved design spec is missing fields: {sorted(missing)}"
        )
    if unknown:
        raise DesignValidationError(
            f"Saved design spec contains unknown fields: {sorted(unknown)}"
        )
    return DesignSpecV1(**{name: value[name] for name in _DESIGN_SPEC_FIELDS})


def design_spec_to_dict(value: DesignSpecV1 | Mapping) -> dict:
    """Return a canonical JSON-primitive representation."""
    spec = normalize_design_spec(value)
    payload = asdict(spec)
    json.dumps(payload, allow_nan=False)
    return payload


def design_spec_fingerprint(
    value: DesignSpecV1 | Mapping,
    *,
    length: int = 16,
) -> str:
    """Return a stable SHA-256 prefix for a canonical design spec."""
    if (
        isinstance(length, bool)
        or not isinstance(length, int)
        or not 8 <= length <= 64
    ):
        raise ValueError("fingerprint length must be an integer between 8 and 64")
    payload = design_spec_to_dict(value)
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:length]


__all__ = [
    "ASSIGNMENT_ALGORITHM_V1",
    "DESIGN_PLAN_CHOICES",
    "DESIGN_SPEC_VERSION",
    "DISCONNECTED_GROUPS",
    "FULLY_CROSSED",
    "NESTED_BRIDGE",
    "SPIRAL",
    "DesignSpecV1",
    "DesignValidationError",
    "design_spec_fingerprint",
    "design_spec_to_dict",
    "normalize_design_spec",
]
