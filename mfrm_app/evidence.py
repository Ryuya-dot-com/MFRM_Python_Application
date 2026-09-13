"""Versioned evidence and decision contracts for the standalone Python app.

The Streamlit entrypoint historically used several human-facing status
vocabularies.  This module provides a small machine-facing layer that keeps
three questions separate:

* was an analysis computable;
* was its conclusion stable under required sensitivity checks; and
* how should the application route the user next.

The contracts intentionally do not encode a pass/fail verdict about a person,
rater, institution, or assessment.  They are dependency-free apart from the
optional pandas conversion helper and never invoke an external engine.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from datetime import date, datetime
from enum import Enum
import hashlib
import json
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence, TypeVar


ANALYSIS_IDENTITY_SCHEMA_VERSION = "mfrm_analysis_identity_v1"
EVIDENCE_RECORD_SCHEMA_VERSION = "mfrm_evidence_record_v1"
SENSITIVITY_RECORD_SCHEMA_VERSION = "mfrm_sensitivity_record_v1"
SENSITIVITY_RULE_SCHEMA_VERSION = "mfrm_sensitivity_rule_v1"
SENSITIVITY_PLAN_SCHEMA_VERSION = "mfrm_sensitivity_plan_v1"
SENSITIVITY_DECISION_SCHEMA_VERSION = "mfrm_sensitivity_decision_v1"
DECISION_RECORD_SCHEMA_VERSION = "mfrm_decision_record_v1"
REASON_REGISTRY_VERSION = "mfrm_reason_registry_v1"
MAX_CONTRACT_INTEGER_BITS = 4096


class ContractValidationError(ValueError):
    """Raised when a versioned evidence payload violates its contract."""


class ContractVersionError(ContractValidationError):
    """Raised when a payload uses a schema version this code cannot read."""


class ComputationState(str, Enum):
    """Whether evidence was computed and is interpretable for its purpose."""

    AVAILABLE = "AVAILABLE"
    CAUTION = "CAUTION"
    HOLD = "HOLD"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"


class StabilityState(str, Enum):
    """Whether a conclusion survives its required sensitivity checks."""

    STABLE = "STABLE"
    CONDITIONALLY_STABLE = "CONDITIONALLY_STABLE"
    SENSITIVE = "SENSITIVE"
    NOT_ASSESSED = "NOT_ASSESSED"


class DecisionDisposition(str, Enum):
    """Workflow routing, not a substantive accept/reject classification."""

    USE = "USE"
    USE_WITH_CAVEAT = "USE_WITH_CAVEAT"
    WITHHOLD = "WITHHOLD"
    BOUNDARY_ONLY = "BOUNDARY_ONLY"
    NOT_EVALUATED = "NOT_EVALUATED"


class ReasonCode(str, Enum):
    """Core reason constants; contracts also preserve valid extension codes."""

    PREREQUISITE_MISSING = "evidence.prerequisite_missing"
    INSUFFICIENT_OBSERVATIONS = "data.insufficient_observations"
    INSUFFICIENT_LEVELS = "data.insufficient_levels"
    INSUFFICIENT_OVERLAP = "data.insufficient_overlap"
    DISCONNECTED_DESIGN = "design.disconnected"
    UNSUPPORTED_MODEL = "model.unsupported"
    RUN_FAILED = "evidence.run_failed"
    NONCOMPARABLE = "evidence.noncomparable"
    NOT_APPLICABLE = "evidence.not_applicable"
    LIMITED_EVIDENCE = "evidence.limited"
    COMPATIBILITY_ONLY = "scope.compatibility_only"
    UNKNOWN_LEGACY_STATUS = "contract.unknown_legacy_status"
    SENSITIVITY_NOT_RUN = "sensitivity.not_run"
    SENSITIVITY_STABLE = "sensitivity.stable"
    SENSITIVITY_THRESHOLD_NEAR = "sensitivity.threshold_near"
    SENSITIVITY_CONCLUSION_CHANGED = "sensitivity.conclusion_changed"
    SENSITIVITY_ESTIMATE_SHIFT = "sensitivity.estimate_shift"
    SENSITIVITY_RANK_SHIFT = "sensitivity.rank_shift"
    SENSITIVITY_MISSING_REQUIRED = "sensitivity.missing_required"
    SENSITIVITY_RULE_NOT_EVALUABLE = "sensitivity.rule_not_evaluable"


class LegacyVocabulary(str, Enum):
    """Named legacy vocabularies whose identical labels have different meaning."""

    INPUT_READINESS = "input_readiness"
    RESOURCE_PREFLIGHT = "resource_preflight"
    FINAL_READINESS = "final_readiness"
    GUIDED = "guided"


class LegacyDispositionVocabulary(str, Enum):
    """Legacy labels that route claims rather than describe computation."""

    CLAIM_GATE = "claim_gate"
    REPORT_DECISION = "report_decision"


_TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]*$")
_REASON_RE = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+$")
_EnumT = TypeVar("_EnumT", bound=Enum)


def _require_text(value: object, name: str, *, token: bool = False) -> str:
    if not isinstance(value, str):
        raise ContractValidationError(f"{name} must be a string")
    text = value.strip()
    if not text:
        raise ContractValidationError(f"{name} must be a non-empty string")
    if token and not _TOKEN_RE.fullmatch(text):
        raise ContractValidationError(
            f"{name} must contain only letters, numbers, '.', '_', ':', or '-'"
        )
    return text


def _validate_payload_shape(
    payload: Mapping[str, object],
    *,
    name: str,
    keys: set[str],
    reason_registry: bool = False,
) -> None:
    """Reject missing, unknown, or version-registry fields before coercion."""
    if not isinstance(payload, Mapping):
        raise ContractValidationError(f"{name} payload must be a mapping")
    actual = set(payload)
    if actual != keys:
        missing = sorted(keys.difference(actual))
        unknown = sorted((repr(key) for key in actual.difference(keys)))
        raise ContractValidationError(
            f"{name} payload shape mismatch; missing={missing!r}, unknown={unknown!r}"
        )
    if reason_registry and payload.get("reason_registry_version") != REASON_REGISTRY_VERSION:
        raise ContractVersionError(
            f"Unsupported reason registry version: {payload.get('reason_registry_version')!r}"
        )


def _payload_sequence(payload: Mapping[str, object], key: str) -> tuple[object, ...]:
    value = payload[key]
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or isinstance(value, Mapping)
    ):
        raise ContractValidationError(f"{key} must be a JSON array")
    return tuple(value)


def _coerce_enum(enum_type: type[_EnumT], value: _EnumT | str, name: str) -> _EnumT:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(str(value))
    except ValueError as exc:
        allowed = ", ".join(member.value for member in enum_type)
        raise ContractValidationError(f"Unknown {name} {value!r}; expected one of: {allowed}") from exc


def _coerce_reason_code(value: ReasonCode | str, name: str = "reason code") -> str:
    text = value.value if isinstance(value, ReasonCode) else str(value).strip()
    if not _REASON_RE.fullmatch(text):
        raise ContractValidationError(
            f"{name} must be a namespaced lowercase code such as 'evidence.run_failed'"
        )
    return text


def _normalize_json(value: object) -> object:
    """Normalize supported values before deterministic JSON serialization."""
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int):
        if value.bit_length() > MAX_CONTRACT_INTEGER_BITS:
            raise ContractValidationError(
                "Contract integers cannot exceed 4096 bits"
            )
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ContractValidationError("Contract payloads cannot contain NaN or infinity")
        return value
    if isinstance(value, Enum):
        return _normalize_json(value.value)
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return {
            item.name: _normalize_json(getattr(value, item.name))
            for item in fields(value)
        }
    if isinstance(value, Mapping):
        normalized: dict[str, object] = {}
        for raw_key, raw_value in value.items():
            if not isinstance(raw_key, str):
                raise ContractValidationError("Contract mapping keys must be strings")
            normalized[raw_key] = _normalize_json(raw_value)
        return {key: normalized[key] for key in sorted(normalized)}
    if isinstance(value, (list, tuple)):
        return [_normalize_json(item) for item in value]
    if isinstance(value, (set, frozenset)):
        items = [_normalize_json(item) for item in value]
        return sorted(items, key=canonical_json)
    item_method = getattr(value, "item", None)
    if callable(item_method):
        try:
            return _normalize_json(item_method())
        except (TypeError, ValueError):
            pass
    raise ContractValidationError(
        f"Unsupported contract value of type {type(value).__name__}; convert it to a JSON value first"
    )


def _freeze_json(value: object) -> object:
    """Normalize and recursively freeze JSON-compatible contract content."""
    normalized = _normalize_json(value)

    def freeze(item: object) -> object:
        if isinstance(item, dict):
            return MappingProxyType({key: freeze(val) for key, val in item.items()})
        if isinstance(item, list):
            return tuple(freeze(val) for val in item)
        return item

    return freeze(normalized)


def _thaw_json(value: object) -> object:
    """Return a defensive plain JSON copy of frozen contract content."""
    if isinstance(value, Mapping):
        return {str(key): _thaw_json(val) for key, val in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(val) for val in value]
    return value


def _freeze_mapping(value: Mapping[str, object], name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ContractValidationError(f"{name} must be a mapping")
    frozen = _freeze_json(value)
    if not isinstance(frozen, Mapping):
        raise ContractValidationError(f"{name} must normalize to a mapping")
    return frozen


def canonical_json(payload: object) -> str:
    """Return deterministic, strict JSON for identifiers and archive payloads."""
    try:
        return json.dumps(
            _normalize_json(payload),
            sort_keys=True,
            ensure_ascii=False,
            allow_nan=False,
            separators=(",", ":"),
        )
    except ContractValidationError:
        raise
    except (OverflowError, ValueError) as exc:
        raise ContractValidationError(
            "Contract payload cannot be represented as canonical JSON"
        ) from exc


def payload_fingerprint(payload: object, *, length: int = 64) -> str:
    """Return a deterministic SHA-256 prefix for a JSON-compatible payload."""
    if not 8 <= int(length) <= 64:
        raise ContractValidationError("fingerprint length must be between 8 and 64")
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()[: int(length)]


def _analysis_core_payload(
    *,
    input_data_fingerprint: str,
    config_fingerprint: str,
    analysis_type: str,
    variant_id: str,
    baseline_reference: str,
    engine_version: str,
    seed: int | None,
) -> dict[str, object]:
    return {
        "schema_version": ANALYSIS_IDENTITY_SCHEMA_VERSION,
        "input_data_fingerprint": input_data_fingerprint,
        "config_fingerprint": config_fingerprint,
        "analysis_type": analysis_type,
        "variant_id": variant_id,
        "baseline_reference": baseline_reference,
        "engine_version": engine_version,
        "seed": seed,
    }


@dataclass(frozen=True)
class AnalysisIdentity:
    """Deterministic identity and lineage for one baseline or variant analysis."""

    analysis_id: str
    input_data_fingerprint: str
    config_fingerprint: str
    analysis_type: str
    variant_id: str
    baseline_analysis_id: str
    engine_version: str
    seed: int | None = None
    schema_version: str = ANALYSIS_IDENTITY_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != ANALYSIS_IDENTITY_SCHEMA_VERSION:
            raise ContractVersionError(
                f"Unsupported analysis identity version: {self.schema_version!r}"
            )
        for name in (
            "analysis_id",
            "input_data_fingerprint",
            "config_fingerprint",
            "analysis_type",
            "variant_id",
            "baseline_analysis_id",
            "engine_version",
        ):
            value = _require_text(
                getattr(self, name),
                name,
                token=name in {"analysis_id", "analysis_type", "variant_id", "baseline_analysis_id"},
            )
            object.__setattr__(self, name, value)
        if isinstance(self.seed, bool) or (self.seed is not None and not isinstance(self.seed, int)):
            raise ContractValidationError("seed must be an integer or None")
        baseline_reference = "" if self.variant_id == "baseline" else self.baseline_analysis_id
        if self.variant_id == "baseline" and self.baseline_analysis_id != self.analysis_id:
            raise ContractValidationError("A baseline analysis must point to its own AnalysisID")
        if self.variant_id != "baseline" and self.baseline_analysis_id == self.analysis_id:
            raise ContractValidationError("A variant must point to a distinct baseline AnalysisID")
        expected = "ana_" + payload_fingerprint(
            _analysis_core_payload(
                input_data_fingerprint=self.input_data_fingerprint,
                config_fingerprint=self.config_fingerprint,
                analysis_type=self.analysis_type,
                variant_id=self.variant_id,
                baseline_reference=baseline_reference,
                engine_version=self.engine_version,
                seed=self.seed,
            ),
            length=24,
        )
        if self.analysis_id != expected:
            raise ContractValidationError(
                f"AnalysisID does not match the contract payload; expected {expected!r}"
            )

    @property
    def is_baseline(self) -> bool:
        return self.variant_id == "baseline"

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "analysis_id": self.analysis_id,
            "input_data_fingerprint": self.input_data_fingerprint,
            "config_fingerprint": self.config_fingerprint,
            "analysis_type": self.analysis_type,
            "variant_id": self.variant_id,
            "baseline_analysis_id": self.baseline_analysis_id,
            "engine_version": self.engine_version,
            "seed": self.seed,
        }

    def to_row(self) -> dict[str, object]:
        return {
            "SchemaVersion": self.schema_version,
            "AnalysisID": self.analysis_id,
            "InputDataFingerprint": self.input_data_fingerprint,
            "ConfigFingerprint": self.config_fingerprint,
            "AnalysisType": self.analysis_type,
            "VariantID": self.variant_id,
            "BaselineAnalysisID": self.baseline_analysis_id,
            "EngineVersion": self.engine_version,
            "Seed": self.seed,
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> "AnalysisIdentity":
        _validate_payload_shape(
            payload,
            name="AnalysisIdentity",
            keys={
                "schema_version",
                "analysis_id",
                "input_data_fingerprint",
                "config_fingerprint",
                "analysis_type",
                "variant_id",
                "baseline_analysis_id",
                "engine_version",
                "seed",
            },
        )
        return cls(
            schema_version=payload["schema_version"],
            analysis_id=payload["analysis_id"],
            input_data_fingerprint=payload["input_data_fingerprint"],
            config_fingerprint=payload["config_fingerprint"],
            analysis_type=payload["analysis_type"],
            variant_id=payload["variant_id"],
            baseline_analysis_id=payload["baseline_analysis_id"],
            engine_version=payload["engine_version"],
            seed=payload["seed"],
        )


def build_analysis_identity(
    *,
    input_data_fingerprint: str,
    resolved_settings: Mapping[str, object],
    analysis_type: str,
    engine_version: str,
    variant_id: str = "baseline",
    baseline: AnalysisIdentity | str | None = None,
    seed: int | None = None,
) -> AnalysisIdentity:
    """Build a deterministic analysis identity from data, settings, and lineage."""
    data_fp = _require_text(input_data_fingerprint, "input_data_fingerprint")
    analysis_type = _require_text(analysis_type, "analysis_type", token=True)
    variant_id = _require_text(variant_id, "variant_id", token=True)
    engine_version = _require_text(engine_version, "engine_version")
    if isinstance(seed, bool) or (seed is not None and not isinstance(seed, int)):
        raise ContractValidationError("seed must be an integer or None")
    config_fp = payload_fingerprint(resolved_settings)

    if variant_id == "baseline":
        if baseline is not None:
            raise ContractValidationError("A baseline analysis cannot name another baseline")
        baseline_reference = ""
    else:
        if baseline is None:
            raise ContractValidationError("A variant analysis must name its baseline AnalysisID")
        baseline_reference = baseline.analysis_id if isinstance(baseline, AnalysisIdentity) else str(baseline)
        baseline_reference = _require_text(
            baseline_reference, "baseline_analysis_id", token=True
        )

    analysis_id = "ana_" + payload_fingerprint(
        _analysis_core_payload(
            input_data_fingerprint=data_fp,
            config_fingerprint=config_fp,
            analysis_type=analysis_type,
            variant_id=variant_id,
            baseline_reference=baseline_reference,
            engine_version=engine_version,
            seed=seed,
        ),
        length=24,
    )
    return AnalysisIdentity(
        analysis_id=analysis_id,
        input_data_fingerprint=data_fp,
        config_fingerprint=config_fp,
        analysis_type=analysis_type,
        variant_id=variant_id,
        baseline_analysis_id=analysis_id if variant_id == "baseline" else baseline_reference,
        engine_version=engine_version,
        seed=seed,
    )


def _text_tuple(values: Iterable[object], name: str) -> tuple[str, ...]:
    normalized = tuple(_require_text(value, name) for value in values)
    if len(set(normalized)) != len(normalized):
        raise ContractValidationError(f"{name} values must be unique")
    return normalized


_SENSITIVITY_RULE_OPERATORS = {"lt", "lte", "gt", "gte", "eq", "ne"}
_MISSING_RULE_VALUE = object()


def _normalize_sensitivity_rule(rule: Mapping[str, object]) -> Mapping[str, object]:
    """Validate and freeze the executable sensitivity-rule schema."""
    if not isinstance(rule, Mapping):
        raise ContractValidationError("decision_rule must be a mapping")
    expected_keys = {"schema_version", "combine", "criteria"}
    if set(rule) != expected_keys:
        raise ContractValidationError(
            "decision_rule must contain exactly schema_version, combine, and criteria"
        )
    if rule.get("schema_version") != SENSITIVITY_RULE_SCHEMA_VERSION:
        raise ContractVersionError(
            f"Unsupported sensitivity rule version: {rule.get('schema_version')!r}"
        )
    combine = str(rule.get("combine", "")).strip()
    if combine not in {"any", "all"}:
        raise ContractValidationError("decision_rule combine must be 'any' or 'all'")
    raw_criteria = rule.get("criteria")
    if (
        not isinstance(raw_criteria, Sequence)
        or isinstance(raw_criteria, (str, bytes, bytearray))
        or not raw_criteria
    ):
        raise ContractValidationError("decision_rule criteria must be a non-empty sequence")
    criteria: list[dict[str, object]] = []
    for index, raw in enumerate(raw_criteria):
        if not isinstance(raw, Mapping) or set(raw) != {"metric", "operator", "value"}:
            raise ContractValidationError(
                f"decision_rule criterion {index} must contain exactly metric, operator, and value"
            )
        metric = _require_text(raw.get("metric"), f"criterion {index} metric", token=True)
        operator = str(raw.get("operator", "")).strip()
        if operator not in _SENSITIVITY_RULE_OPERATORS:
            raise ContractValidationError(
                f"Unknown decision_rule operator {operator!r}; expected one of "
                f"{sorted(_SENSITIVITY_RULE_OPERATORS)!r}"
            )
        value = _normalize_json(raw.get("value"))
        if isinstance(value, (dict, list)):
            raise ContractValidationError("decision_rule criterion values must be JSON scalars")
        if operator in {"lt", "lte", "gt", "gte"} and (
            isinstance(value, bool) or not isinstance(value, (int, float))
        ):
            raise ContractValidationError(
                f"decision_rule operator {operator!r} requires a finite numeric value"
            )
        criteria.append({"metric": metric, "operator": operator, "value": value})
    criteria.sort(key=canonical_json)
    if len({canonical_json(item) for item in criteria}) != len(criteria):
        raise ContractValidationError("decision_rule criteria must be unique")
    return _freeze_mapping(
        {
            "schema_version": SENSITIVITY_RULE_SCHEMA_VERSION,
            "combine": combine,
            "criteria": criteria,
        },
        "decision_rule",
    )


def _evaluate_sensitivity_rule(
    rule: Mapping[str, object],
    metrics: Mapping[str, object],
) -> bool | None:
    """Evaluate a validated rule using three-valued, fail-closed logic."""
    normalized_rule = _normalize_sensitivity_rule(rule)
    outcomes: list[bool | None] = []
    for criterion in normalized_rule["criteria"]:
        metric_name = str(criterion["metric"])
        observed = metrics.get(metric_name, _MISSING_RULE_VALUE)
        if observed is _MISSING_RULE_VALUE:
            outcomes.append(None)
            continue
        operator = str(criterion["operator"])
        expected = criterion["value"]
        if operator in {"lt", "lte", "gt", "gte"}:
            if (
                isinstance(observed, bool)
                or not isinstance(observed, (int, float))
                or (isinstance(observed, float) and not math.isfinite(observed))
            ):
                outcomes.append(None)
                continue
            outcomes.append(
                {
                    "lt": observed < expected,
                    "lte": observed <= expected,
                    "gt": observed > expected,
                    "gte": observed >= expected,
                }[operator]
            )
        else:
            try:
                comparable = _normalize_json(observed)
            except ContractValidationError:
                outcomes.append(None)
                continue
            if (
                comparable is None
                or expected is None
                or isinstance(comparable, (dict, list))
                or isinstance(expected, (dict, list))
            ):
                outcomes.append(None)
                continue
            if isinstance(expected, bool):
                same_type = isinstance(comparable, bool)
            elif isinstance(expected, (int, float)) and not isinstance(expected, bool):
                same_type = isinstance(comparable, (int, float)) and not isinstance(
                    comparable, bool
                )
            elif isinstance(expected, str):
                same_type = isinstance(comparable, str)
            else:
                same_type = False
            if not same_type:
                outcomes.append(None)
                continue
            equal = comparable == expected
            outcomes.append(equal if operator == "eq" else not equal)

    if normalized_rule["combine"] == "any":
        if True in outcomes:
            return True
        return None if None in outcomes else False
    if False in outcomes:
        return False
    return None if None in outcomes else True


@dataclass(frozen=True)
class EvidenceRecord:
    """One traceable piece of evidence for a named analysis question."""

    evidence_id: str
    evidence_key: str
    analysis_id: str
    domain: str
    question: str
    computation_state: ComputationState
    stability_state: StabilityState
    summary: str
    interpretation_boundary: str
    recommended_action: str
    required: bool = True
    reason_code: str | None = None
    reason_args: Mapping[str, object] = field(default_factory=dict)
    scope: Mapping[str, object] = field(default_factory=dict)
    prerequisites: tuple[str, ...] = ()
    observed: Mapping[str, object] = field(default_factory=dict)
    uncertainty: Mapping[str, object] = field(default_factory=dict)
    next_inspection: str = ""
    source_artifacts: tuple[str, ...] = ()
    sensitivity_ids: tuple[str, ...] = ()
    schema_version: str = EVIDENCE_RECORD_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != EVIDENCE_RECORD_SCHEMA_VERSION:
            raise ContractVersionError(
                f"Unsupported evidence record version: {self.schema_version!r}"
            )
        for name in (
            "evidence_id",
            "evidence_key",
            "analysis_id",
            "domain",
            "question",
            "summary",
            "interpretation_boundary",
            "recommended_action",
        ):
            value = _require_text(
                getattr(self, name),
                name,
                token=name in {"evidence_id", "evidence_key", "analysis_id"},
            )
            object.__setattr__(self, name, value)
        object.__setattr__(
            self,
            "computation_state",
            _coerce_enum(ComputationState, self.computation_state, "computation state"),
        )
        object.__setattr__(
            self,
            "stability_state",
            _coerce_enum(StabilityState, self.stability_state, "stability state"),
        )
        if not isinstance(self.required, bool):
            raise ContractValidationError("required must be a boolean")
        if self.reason_code is not None:
            object.__setattr__(self, "reason_code", _coerce_reason_code(self.reason_code))
        object.__setattr__(self, "prerequisites", _text_tuple(self.prerequisites, "prerequisite"))
        object.__setattr__(self, "source_artifacts", _text_tuple(self.source_artifacts, "source artifact"))
        object.__setattr__(self, "sensitivity_ids", _text_tuple(self.sensitivity_ids, "sensitivity ID"))
        for name in ("reason_args", "scope", "observed", "uncertainty"):
            object.__setattr__(self, name, _freeze_mapping(getattr(self, name), name))
        if not isinstance(self.next_inspection, str):
            raise ContractValidationError("next_inspection must be a string")

        if self.computation_state is ComputationState.AVAILABLE and self.reason_code is not None:
            raise ContractValidationError("AVAILABLE evidence cannot carry a warning reason code")
        if self.computation_state is not ComputationState.AVAILABLE and self.reason_code is None:
            raise ContractValidationError(
                f"{self.computation_state.value} evidence requires a stable ReasonCode"
            )
        if self.computation_state in {ComputationState.AVAILABLE, ComputationState.CAUTION}:
            if not self.source_artifacts:
                raise ContractValidationError("Computed evidence requires at least one source artifact")
        elif self.stability_state is not StabilityState.NOT_ASSESSED:
            raise ContractValidationError(
                "HOLD or NOT_ASSESSABLE evidence cannot claim conclusion stability"
            )
        if (
            self.stability_state is not StabilityState.NOT_ASSESSED
            and not self.sensitivity_ids
        ):
            raise ContractValidationError(
                "Assessed conclusion stability requires at least one SensitivityID"
            )
        expected_id = "evi_" + payload_fingerprint(
            {
                "schema_version": self.schema_version,
                "analysis_id": self.analysis_id,
                "evidence_key": self.evidence_key,
            },
            length=24,
        )
        if self.evidence_id != expected_id:
            raise ContractValidationError(
                f"EvidenceID does not match the contract payload; expected {expected_id!r}"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "evidence_id": self.evidence_id,
            "evidence_key": self.evidence_key,
            "analysis_id": self.analysis_id,
            "domain": self.domain,
            "question": self.question,
            "computation_state": self.computation_state.value,
            "stability_state": self.stability_state.value,
            "required": self.required,
            "reason_code": self.reason_code,
            "reason_registry_version": REASON_REGISTRY_VERSION,
            "reason_args": _thaw_json(self.reason_args),
            "summary": self.summary,
            "scope": _thaw_json(self.scope),
            "prerequisites": list(self.prerequisites),
            "observed": _thaw_json(self.observed),
            "uncertainty": _thaw_json(self.uncertainty),
            "interpretation_boundary": self.interpretation_boundary,
            "recommended_action": self.recommended_action,
            "next_inspection": self.next_inspection,
            "source_artifacts": list(self.source_artifacts),
            "sensitivity_ids": list(self.sensitivity_ids),
        }

    def to_row(self) -> dict[str, object]:
        payload = self.to_payload()
        return {
            "SchemaVersion": self.schema_version,
            "EvidenceID": self.evidence_id,
            "EvidenceKey": self.evidence_key,
            "AnalysisID": self.analysis_id,
            "Domain": self.domain,
            "Question": self.question,
            "ComputationState": self.computation_state.value,
            "StabilityState": self.stability_state.value,
            "Required": self.required,
            "ReasonCode": self.reason_code or "",
            "ReasonRegistryVersion": REASON_REGISTRY_VERSION,
            "ReasonArgsJSON": canonical_json(payload["reason_args"]),
            "Summary": self.summary,
            "ScopeJSON": canonical_json(payload["scope"]),
            "PrerequisitesJSON": canonical_json(payload["prerequisites"]),
            "ObservedJSON": canonical_json(payload["observed"]),
            "UncertaintyJSON": canonical_json(payload["uncertainty"]),
            "InterpretationBoundary": self.interpretation_boundary,
            "RecommendedAction": self.recommended_action,
            "NextInspection": self.next_inspection,
            "SourceArtifactsJSON": canonical_json(payload["source_artifacts"]),
            "SensitivityIDsJSON": canonical_json(payload["sensitivity_ids"]),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> "EvidenceRecord":
        _validate_payload_shape(
            payload,
            name="EvidenceRecord",
            keys={
                "schema_version",
                "evidence_id",
                "evidence_key",
                "analysis_id",
                "domain",
                "question",
                "computation_state",
                "stability_state",
                "required",
                "reason_code",
                "reason_registry_version",
                "reason_args",
                "summary",
                "scope",
                "prerequisites",
                "observed",
                "uncertainty",
                "interpretation_boundary",
                "recommended_action",
                "next_inspection",
                "source_artifacts",
                "sensitivity_ids",
            },
            reason_registry=True,
        )
        return cls(
            schema_version=payload["schema_version"],
            evidence_id=payload["evidence_id"],
            evidence_key=payload["evidence_key"],
            analysis_id=payload["analysis_id"],
            domain=payload["domain"],
            question=payload["question"],
            computation_state=payload["computation_state"],
            stability_state=payload["stability_state"],
            required=payload["required"],
            reason_code=payload["reason_code"],
            reason_args=payload["reason_args"],
            summary=payload["summary"],
            scope=payload["scope"],
            prerequisites=_payload_sequence(payload, "prerequisites"),
            observed=payload["observed"],
            uncertainty=payload["uncertainty"],
            interpretation_boundary=payload["interpretation_boundary"],
            recommended_action=payload["recommended_action"],
            next_inspection=payload["next_inspection"],
            source_artifacts=_payload_sequence(payload, "source_artifacts"),
            sensitivity_ids=_payload_sequence(payload, "sensitivity_ids"),
        )


def make_evidence_record(
    identity: AnalysisIdentity | str,
    *,
    evidence_key: str,
    domain: str,
    question: str,
    computation_state: ComputationState | str,
    summary: str,
    interpretation_boundary: str,
    recommended_action: str,
    stability_state: StabilityState | str = StabilityState.NOT_ASSESSED,
    reason_code: ReasonCode | str | None = None,
    required: bool = True,
    reason_args: Mapping[str, object] | None = None,
    scope: Mapping[str, object] | None = None,
    prerequisites: Iterable[str] = (),
    observed: Mapping[str, object] | None = None,
    uncertainty: Mapping[str, object] | None = None,
    next_inspection: str = "",
    source_artifacts: Iterable[str] = (),
    sensitivity_ids: Iterable[str] = (),
) -> EvidenceRecord:
    """Create a validated evidence record with a deterministic EvidenceID."""
    analysis_id = identity.analysis_id if isinstance(identity, AnalysisIdentity) else str(identity)
    evidence_key = _require_text(evidence_key, "evidence_key", token=True)
    evidence_id = "evi_" + payload_fingerprint(
        {
            "schema_version": EVIDENCE_RECORD_SCHEMA_VERSION,
            "analysis_id": analysis_id,
            "evidence_key": evidence_key,
        },
        length=24,
    )
    return EvidenceRecord(
        evidence_id=evidence_id,
        evidence_key=evidence_key,
        analysis_id=analysis_id,
        domain=domain,
        question=question,
        computation_state=computation_state,
        stability_state=stability_state,
        required=required,
        reason_code=reason_code,
        reason_args=reason_args or {},
        summary=summary,
        scope=scope or {},
        prerequisites=tuple(prerequisites),
        observed=observed or {},
        uncertainty=uncertainty or {},
        interpretation_boundary=interpretation_boundary,
        recommended_action=recommended_action,
        next_inspection=next_inspection,
        source_artifacts=tuple(source_artifacts),
        sensitivity_ids=tuple(sensitivity_ids),
    )


EVIDENCE_ROW_COLUMNS = (
    "SchemaVersion",
    "EvidenceID",
    "EvidenceKey",
    "AnalysisID",
    "Domain",
    "Question",
    "ComputationState",
    "StabilityState",
    "Required",
    "ReasonCode",
    "ReasonRegistryVersion",
    "ReasonArgsJSON",
    "Summary",
    "ScopeJSON",
    "PrerequisitesJSON",
    "ObservedJSON",
    "UncertaintyJSON",
    "InterpretationBoundary",
    "RecommendedAction",
    "NextInspection",
    "SourceArtifactsJSON",
    "SensitivityIDsJSON",
)


def assert_matching_identity(
    records: Iterable[EvidenceRecord],
    identity: AnalysisIdentity | str,
) -> None:
    """Reject evidence that was produced by a different analysis identity."""
    expected = identity.analysis_id if isinstance(identity, AnalysisIdentity) else str(identity)
    mismatched = sorted({record.analysis_id for record in records if record.analysis_id != expected})
    if mismatched:
        raise ContractValidationError(
            f"Evidence identity mismatch: expected {expected!r}, found {mismatched!r}"
        )


def validate_evidence_records(
    records: Iterable[EvidenceRecord],
    *,
    identity: AnalysisIdentity | str | None = None,
) -> tuple[EvidenceRecord, ...]:
    """Validate uniqueness and optional analysis identity for an evidence ledger."""
    normalized = tuple(records)
    if any(not isinstance(record, EvidenceRecord) for record in normalized):
        raise ContractValidationError("All evidence ledger entries must be EvidenceRecord instances")
    ids = [record.evidence_id for record in normalized]
    if len(ids) != len(set(ids)):
        raise ContractValidationError("EvidenceID values must be unique within a ledger")
    if identity is not None:
        assert_matching_identity(normalized, identity)
    return normalized


def evidence_records_to_frame(
    records: Iterable[EvidenceRecord],
    *,
    identity: AnalysisIdentity | str | None = None,
):
    """Return a DataFrame with a fixed CSV/export column order."""
    import pandas as pd

    normalized = validate_evidence_records(records, identity=identity)
    return pd.DataFrame(
        [record.to_row() for record in normalized],
        columns=EVIDENCE_ROW_COLUMNS,
    )


_LEGACY_STATE_MAP: dict[LegacyVocabulary, dict[str, ComputationState | None]] = {
    LegacyVocabulary.INPUT_READINESS: {
        "ok": ComputationState.AVAILABLE,
        "warning": ComputationState.CAUTION,
        "issue": ComputationState.HOLD,
    },
    LegacyVocabulary.RESOURCE_PREFLIGHT: {
        "OK": ComputationState.AVAILABLE,
        "Review": ComputationState.CAUTION,
        "Block": ComputationState.HOLD,
    },
    LegacyVocabulary.FINAL_READINESS: {
        "Ready": ComputationState.AVAILABLE,
        "OK": ComputationState.AVAILABLE,
        "Review": None,
        "Missing": ComputationState.NOT_ASSESSABLE,
        "Not ready": ComputationState.HOLD,
    },
    LegacyVocabulary.GUIDED: {
        "OK": ComputationState.AVAILABLE,
        "Caution": ComputationState.CAUTION,
        "Review": None,
        "Skipped": ComputationState.NOT_ASSESSABLE,
        "Do not interpret yet": ComputationState.HOLD,
    },
}

_LEGACY_DISPOSITION_MAP: dict[
    LegacyDispositionVocabulary,
    dict[str, DecisionDisposition],
] = {
    LegacyDispositionVocabulary.CLAIM_GATE: {
        "Ready": DecisionDisposition.USE,
        "Report with caveat": DecisionDisposition.USE_WITH_CAVEAT,
        "Boundary": DecisionDisposition.BOUNDARY_ONLY,
        "Do not claim": DecisionDisposition.WITHHOLD,
        "Not ready": DecisionDisposition.WITHHOLD,
    },
    LegacyDispositionVocabulary.REPORT_DECISION: {
        "Ready": DecisionDisposition.USE,
        "Needs caveat": DecisionDisposition.USE_WITH_CAVEAT,
        "Needs rerun": DecisionDisposition.WITHHOLD,
        "Do not report yet": DecisionDisposition.WITHHOLD,
    },
}


def normalize_legacy_status(
    status: object,
    *,
    vocabulary: LegacyVocabulary | str,
    explicit_state: ComputationState | str | None = None,
) -> ComputationState:
    """Map computation labels, requiring context for ambiguous legacy values."""
    vocab = _coerce_enum(LegacyVocabulary, vocabulary, "legacy vocabulary")
    label = str(status).strip()
    try:
        mapped = _LEGACY_STATE_MAP[vocab][label]
    except KeyError as exc:
        raise ContractValidationError(
            f"{ReasonCode.UNKNOWN_LEGACY_STATUS.value}: {label!r} is not valid for {vocab.value}"
        ) from exc
    if explicit_state is not None:
        explicit = _coerce_enum(
            ComputationState,
            explicit_state,
            "explicit computation state",
        )
        if mapped is not None and explicit is not mapped:
            raise ContractValidationError(
                f"Explicit state {explicit.value} conflicts with unambiguous {vocab.value} label {label!r}"
            )
        return explicit
    if mapped is None:
        raise ContractValidationError(
            f"{ReasonCode.UNKNOWN_LEGACY_STATUS.value}: {label!r} is ambiguous for "
            f"{vocab.value}; supply explicit_state"
        )
    return mapped


def normalize_legacy_disposition(
    status: object,
    *,
    vocabulary: LegacyDispositionVocabulary | str,
) -> DecisionDisposition:
    """Map legacy claim/report routing without altering computation state."""
    vocab = _coerce_enum(
        LegacyDispositionVocabulary,
        vocabulary,
        "legacy disposition vocabulary",
    )
    label = str(status).strip()
    try:
        return _LEGACY_DISPOSITION_MAP[vocab][label]
    except KeyError as exc:
        raise ContractValidationError(
            f"{ReasonCode.UNKNOWN_LEGACY_STATUS.value}: {label!r} is not valid for {vocab.value}"
        ) from exc


@dataclass(frozen=True)
class SensitivityRecord:
    """Outcome for one required or optional baseline-to-variant comparison."""

    sensitivity_id: str
    baseline_analysis_id: str
    variant_analysis_id: str
    variant_id: str
    conclusion_id: str
    required: bool
    credible: bool
    computation_state: ComputationState
    comparable: bool
    summary: str
    evidence_ids: tuple[str, ...]
    source_artifacts: tuple[str, ...]
    reason_code: str | None = None
    metrics: Mapping[str, object] = field(default_factory=dict)
    schema_version: str = SENSITIVITY_RECORD_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != SENSITIVITY_RECORD_SCHEMA_VERSION:
            raise ContractVersionError(
                f"Unsupported sensitivity record version: {self.schema_version!r}"
            )
        for name in (
            "sensitivity_id",
            "baseline_analysis_id",
            "variant_analysis_id",
            "variant_id",
            "conclusion_id",
            "summary",
        ):
            value = _require_text(
                getattr(self, name),
                name,
                token=name != "summary",
            )
            object.__setattr__(self, name, value)
        if self.baseline_analysis_id == self.variant_analysis_id:
            raise ContractValidationError("Sensitivity variants must differ from their baseline")
        if not isinstance(self.required, bool) or not isinstance(self.credible, bool):
            raise ContractValidationError("required and credible must be booleans")
        object.__setattr__(
            self,
            "computation_state",
            _coerce_enum(ComputationState, self.computation_state, "computation state"),
        )
        if self.reason_code is not None:
            object.__setattr__(self, "reason_code", _coerce_reason_code(self.reason_code))
        if not isinstance(self.comparable, bool):
            raise ContractValidationError("comparable must be a boolean")
        object.__setattr__(self, "evidence_ids", _text_tuple(self.evidence_ids, "evidence ID"))
        object.__setattr__(self, "source_artifacts", _text_tuple(self.source_artifacts, "source artifact"))
        object.__setattr__(self, "metrics", _freeze_mapping(self.metrics, "metrics"))

        computed = self.computation_state in {
            ComputationState.AVAILABLE,
            ComputationState.CAUTION,
        }
        if computed:
            if not self.comparable:
                raise ContractValidationError(
                    "A computed sensitivity comparison must be comparable"
                )
            if not self.evidence_ids or not self.source_artifacts:
                raise ContractValidationError(
                    "A computed sensitivity comparison requires EvidenceID and source artifact links"
                )
        elif self.comparable:
            raise ContractValidationError(
                "An uncomputed sensitivity comparison cannot claim comparability"
            )
        if self.computation_state is ComputationState.AVAILABLE and self.reason_code is not None:
            raise ContractValidationError("AVAILABLE sensitivity evidence cannot carry a warning reason")
        if self.computation_state is not ComputationState.AVAILABLE and self.reason_code is None:
            raise ContractValidationError(
                f"{self.computation_state.value} sensitivity evidence requires a stable ReasonCode"
            )
        expected_id = "sen_" + payload_fingerprint(
            {
                "schema_version": self.schema_version,
                "baseline_analysis_id": self.baseline_analysis_id,
                "variant_analysis_id": self.variant_analysis_id,
                "conclusion_id": self.conclusion_id,
            },
            length=24,
        )
        if self.sensitivity_id != expected_id:
            raise ContractValidationError(
                f"SensitivityID does not match the contract payload; expected {expected_id!r}"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "sensitivity_id": self.sensitivity_id,
            "baseline_analysis_id": self.baseline_analysis_id,
            "variant_analysis_id": self.variant_analysis_id,
            "variant_id": self.variant_id,
            "conclusion_id": self.conclusion_id,
            "required": self.required,
            "credible": self.credible,
            "computation_state": self.computation_state.value,
            "comparable": self.comparable,
            "reason_code": self.reason_code,
            "reason_registry_version": REASON_REGISTRY_VERSION,
            "summary": self.summary,
            "metrics": _thaw_json(self.metrics),
            "evidence_ids": list(self.evidence_ids),
            "source_artifacts": list(self.source_artifacts),
        }

    def to_row(self) -> dict[str, object]:
        return {
            "SchemaVersion": self.schema_version,
            "SensitivityID": self.sensitivity_id,
            "BaselineAnalysisID": self.baseline_analysis_id,
            "VariantAnalysisID": self.variant_analysis_id,
            "VariantID": self.variant_id,
            "ConclusionID": self.conclusion_id,
            "Required": self.required,
            "Credible": self.credible,
            "ComputationState": self.computation_state.value,
            "Comparable": self.comparable,
            "ReasonCode": self.reason_code or "",
            "ReasonRegistryVersion": REASON_REGISTRY_VERSION,
            "Summary": self.summary,
            "MetricsJSON": canonical_json(self.metrics),
            "EvidenceIDsJSON": canonical_json(self.evidence_ids),
            "SourceArtifactsJSON": canonical_json(self.source_artifacts),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> "SensitivityRecord":
        _validate_payload_shape(
            payload,
            name="SensitivityRecord",
            keys={
                "schema_version",
                "sensitivity_id",
                "baseline_analysis_id",
                "variant_analysis_id",
                "variant_id",
                "conclusion_id",
                "required",
                "credible",
                "computation_state",
                "comparable",
                "reason_code",
                "reason_registry_version",
                "summary",
                "metrics",
                "evidence_ids",
                "source_artifacts",
            },
            reason_registry=True,
        )
        return cls(
            schema_version=payload["schema_version"],
            sensitivity_id=payload["sensitivity_id"],
            baseline_analysis_id=payload["baseline_analysis_id"],
            variant_analysis_id=payload["variant_analysis_id"],
            variant_id=payload["variant_id"],
            conclusion_id=payload["conclusion_id"],
            required=payload["required"],
            credible=payload["credible"],
            computation_state=payload["computation_state"],
            comparable=payload["comparable"],
            reason_code=payload["reason_code"],
            summary=payload["summary"],
            metrics=payload["metrics"],
            evidence_ids=_payload_sequence(payload, "evidence_ids"),
            source_artifacts=_payload_sequence(payload, "source_artifacts"),
        )


def make_sensitivity_record(
    baseline: AnalysisIdentity,
    variant: AnalysisIdentity,
    *,
    conclusion_id: str,
    required: bool,
    credible: bool,
    computation_state: ComputationState | str,
    comparable: bool,
    summary: str,
    evidence_ids: Iterable[str] = (),
    source_artifacts: Iterable[str] = (),
    reason_code: ReasonCode | str | None = None,
    metrics: Mapping[str, object] | None = None,
) -> SensitivityRecord:
    """Create a validated sensitivity comparison linked to its baseline."""
    if not baseline.is_baseline:
        raise ContractValidationError("baseline must be a baseline AnalysisIdentity")
    if variant.is_baseline or variant.baseline_analysis_id != baseline.analysis_id:
        raise ContractValidationError("variant must point to the supplied baseline AnalysisID")
    if variant.input_data_fingerprint != baseline.input_data_fingerprint:
        raise ContractValidationError(
            "Sensitivity comparisons require the baseline input-data fingerprint"
        )
    if variant.analysis_type != baseline.analysis_type:
        raise ContractValidationError(
            "Sensitivity comparisons require the baseline analysis type"
        )
    conclusion_id = _require_text(conclusion_id, "conclusion_id", token=True)
    sensitivity_id = "sen_" + payload_fingerprint(
        {
            "schema_version": SENSITIVITY_RECORD_SCHEMA_VERSION,
            "baseline_analysis_id": baseline.analysis_id,
            "variant_analysis_id": variant.analysis_id,
            "conclusion_id": conclusion_id,
        },
        length=24,
    )
    return SensitivityRecord(
        sensitivity_id=sensitivity_id,
        baseline_analysis_id=baseline.analysis_id,
        variant_analysis_id=variant.analysis_id,
        variant_id=variant.variant_id,
        conclusion_id=conclusion_id,
        required=required,
        credible=credible,
        computation_state=computation_state,
        comparable=comparable,
        reason_code=reason_code,
        summary=summary,
        metrics=metrics or {},
        evidence_ids=tuple(evidence_ids),
        source_artifacts=tuple(source_artifacts),
    )


@dataclass(frozen=True)
class SensitivityPlan:
    """Prespecified variant identities and executable rule for one conclusion."""

    plan_id: str
    baseline_analysis_id: str
    conclusion_id: str
    baseline_identity_payload: Mapping[str, object]
    required_variant_analysis_ids: Mapping[str, object]
    required_variant_identity_payloads: Mapping[str, object]
    decision_rule: Mapping[str, object]
    schema_version: str = SENSITIVITY_PLAN_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != SENSITIVITY_PLAN_SCHEMA_VERSION:
            raise ContractVersionError(
                f"Unsupported sensitivity plan version: {self.schema_version!r}"
            )
        for name in ("plan_id", "baseline_analysis_id", "conclusion_id"):
            object.__setattr__(
                self,
                name,
                _require_text(getattr(self, name), name, token=True),
            )
        baseline_payload = _freeze_mapping(
            self.baseline_identity_payload,
            "baseline_identity_payload",
        )
        baseline_identity = AnalysisIdentity.from_payload(_thaw_json(baseline_payload))
        if (
            not baseline_identity.is_baseline
            or baseline_identity.analysis_id != self.baseline_analysis_id
        ):
            raise ContractValidationError(
                "baseline_identity_payload must resolve the plan baseline AnalysisID"
            )
        object.__setattr__(self, "baseline_identity_payload", baseline_payload)
        required_variants = _freeze_mapping(
            self.required_variant_analysis_ids,
            "required_variant_analysis_ids",
        )
        if not required_variants:
            raise ContractValidationError("A sensitivity plan requires at least one variant")
        for variant_id, analysis_id in required_variants.items():
            _require_text(variant_id, "required variant ID", token=True)
            if not isinstance(analysis_id, str):
                raise ContractValidationError(
                    "required variant AnalysisIDs must be strings"
                )
            _require_text(analysis_id, "required variant AnalysisID", token=True)
            if variant_id == "baseline":
                raise ContractValidationError("A sensitivity plan cannot use baseline as a variant")
            if analysis_id == self.baseline_analysis_id:
                raise ContractValidationError(
                    "A sensitivity variant AnalysisID must differ from its baseline"
                )
        if len(set(required_variants.values())) != len(required_variants):
            raise ContractValidationError(
                "Required variants must map one-to-one to distinct AnalysisIDs"
            )
        object.__setattr__(self, "required_variant_analysis_ids", required_variants)
        variant_payloads = _freeze_mapping(
            self.required_variant_identity_payloads,
            "required_variant_identity_payloads",
        )
        if set(variant_payloads) != set(required_variants):
            raise ContractValidationError(
                "required_variant_identity_payloads must contain exactly the planned variants"
            )
        for variant_id, raw_payload in variant_payloads.items():
            if not isinstance(raw_payload, Mapping):
                raise ContractValidationError(
                    "Each required variant identity payload must be a mapping"
                )
            variant = AnalysisIdentity.from_payload(_thaw_json(raw_payload))
            if variant.variant_id != variant_id:
                raise ContractValidationError(
                    "A variant identity payload does not match its planned VariantID"
                )
            if variant.analysis_id != required_variants[variant_id]:
                raise ContractValidationError(
                    "A variant identity payload does not match its planned AnalysisID"
                )
            if variant.is_baseline or variant.baseline_analysis_id != baseline_identity.analysis_id:
                raise ContractValidationError(
                    "Every planned variant identity must resolve to the plan baseline"
                )
            if variant.input_data_fingerprint != baseline_identity.input_data_fingerprint:
                raise ContractValidationError(
                    "Planned variant identities must use the baseline input data"
                )
            if variant.analysis_type != baseline_identity.analysis_type:
                raise ContractValidationError(
                    "Planned variant identities must use the baseline analysis type"
                )
        object.__setattr__(
            self,
            "required_variant_identity_payloads",
            variant_payloads,
        )
        object.__setattr__(
            self,
            "decision_rule",
            _normalize_sensitivity_rule(self.decision_rule),
        )
        expected_id = "spl_" + payload_fingerprint(
            {
                "schema_version": self.schema_version,
                "baseline_analysis_id": self.baseline_analysis_id,
                "conclusion_id": self.conclusion_id,
                "baseline_identity_payload": self.baseline_identity_payload,
                "required_variant_analysis_ids": self.required_variant_analysis_ids,
                "required_variant_identity_payloads": self.required_variant_identity_payloads,
                "decision_rule": self.decision_rule,
            },
            length=24,
        )
        if self.plan_id != expected_id:
            raise ContractValidationError(
                f"SensitivityPlanID does not match the contract payload; expected {expected_id!r}"
            )

    @property
    def required_variant_ids(self) -> tuple[str, ...]:
        return tuple(self.required_variant_analysis_ids)

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "plan_id": self.plan_id,
            "baseline_analysis_id": self.baseline_analysis_id,
            "conclusion_id": self.conclusion_id,
            "baseline_identity_payload": _thaw_json(self.baseline_identity_payload),
            "required_variant_analysis_ids": _thaw_json(
                self.required_variant_analysis_ids
            ),
            "required_variant_identity_payloads": _thaw_json(
                self.required_variant_identity_payloads
            ),
            "decision_rule": _thaw_json(self.decision_rule),
        }

    def to_row(self) -> dict[str, object]:
        return {
            "SchemaVersion": self.schema_version,
            "SensitivityPlanID": self.plan_id,
            "BaselineAnalysisID": self.baseline_analysis_id,
            "ConclusionID": self.conclusion_id,
            "BaselineIdentityJSON": canonical_json(self.baseline_identity_payload),
            "RequiredVariantIDsJSON": canonical_json(self.required_variant_ids),
            "RequiredVariantAnalysisIDsJSON": canonical_json(
                self.required_variant_analysis_ids
            ),
            "RequiredVariantIdentitiesJSON": canonical_json(
                self.required_variant_identity_payloads
            ),
            "DecisionRuleJSON": canonical_json(self.decision_rule),
        }

    @classmethod
    def from_payload(cls, payload: Mapping[str, object]) -> "SensitivityPlan":
        _validate_payload_shape(
            payload,
            name="SensitivityPlan",
            keys={
                "schema_version",
                "plan_id",
                "baseline_analysis_id",
                "conclusion_id",
                "baseline_identity_payload",
                "required_variant_analysis_ids",
                "required_variant_identity_payloads",
                "decision_rule",
            },
        )
        return cls(
            schema_version=payload["schema_version"],
            plan_id=payload["plan_id"],
            baseline_analysis_id=payload["baseline_analysis_id"],
            conclusion_id=payload["conclusion_id"],
            baseline_identity_payload=payload["baseline_identity_payload"],
            required_variant_analysis_ids=payload["required_variant_analysis_ids"],
            required_variant_identity_payloads=payload[
                "required_variant_identity_payloads"
            ],
            decision_rule=payload["decision_rule"],
        )


def make_sensitivity_rule(
    *,
    criteria: Iterable[Mapping[str, object]],
    combine: str = "any",
) -> Mapping[str, object]:
    """Build an immutable, order-independent executable sensitivity rule."""
    return _normalize_sensitivity_rule(
        {
            "schema_version": SENSITIVITY_RULE_SCHEMA_VERSION,
            "combine": combine,
            "criteria": tuple(criteria),
        }
    )


def make_sensitivity_plan(
    baseline: AnalysisIdentity,
    *,
    conclusion_id: str,
    required_variants: Iterable[AnalysisIdentity],
    decision_rule: Mapping[str, object],
) -> SensitivityPlan:
    """Create a fingerprinted plan that fixes exact variants before aggregation."""
    if not isinstance(baseline, AnalysisIdentity) or not baseline.is_baseline:
        raise ContractValidationError("A sensitivity plan must reference a baseline analysis")
    conclusion_id = _require_text(conclusion_id, "conclusion_id", token=True)
    variants = tuple(required_variants)
    if not variants or any(not isinstance(variant, AnalysisIdentity) for variant in variants):
        raise ContractValidationError(
            "required_variants must contain at least one AnalysisIdentity"
        )
    variant_map: dict[str, str] = {}
    for variant in variants:
        if variant.is_baseline or variant.baseline_analysis_id != baseline.analysis_id:
            raise ContractValidationError(
                "Every planned variant must point to the supplied baseline AnalysisID"
            )
        if variant.input_data_fingerprint != baseline.input_data_fingerprint:
            raise ContractValidationError(
                "Sensitivity variants must use the baseline input-data fingerprint"
            )
        if variant.analysis_type != baseline.analysis_type:
            raise ContractValidationError(
                "Sensitivity variants must use the baseline analysis type"
            )
        if variant.variant_id in variant_map:
            raise ContractValidationError("Required variant IDs must be unique")
        variant_map[variant.variant_id] = variant.analysis_id
    variant_map = {key: variant_map[key] for key in sorted(variant_map)}
    variant_payloads = {
        variant.variant_id: variant.to_payload()
        for variant in sorted(variants, key=lambda item: item.variant_id)
    }
    frozen_rule = _normalize_sensitivity_rule(decision_rule)
    plan_id = "spl_" + payload_fingerprint(
        {
            "schema_version": SENSITIVITY_PLAN_SCHEMA_VERSION,
            "baseline_analysis_id": baseline.analysis_id,
            "conclusion_id": conclusion_id,
            "baseline_identity_payload": baseline.to_payload(),
            "required_variant_analysis_ids": variant_map,
            "required_variant_identity_payloads": variant_payloads,
            "decision_rule": frozen_rule,
        },
        length=24,
    )
    return SensitivityPlan(
        plan_id=plan_id,
        baseline_analysis_id=baseline.analysis_id,
        conclusion_id=conclusion_id,
        baseline_identity_payload=baseline.to_payload(),
        required_variant_analysis_ids=variant_map,
        required_variant_identity_payloads=variant_payloads,
        decision_rule=frozen_rule,
    )


@dataclass(frozen=True)
class SensitivityDecision:
    """Deterministic aggregate of the required credible sensitivity variants."""

    decision_id: str
    plan_id: str
    baseline_analysis_id: str
    conclusion_id: str
    stability_state: StabilityState
    reason_code: str
    summary: str
    required_variant_analysis_ids: Mapping[str, object]
    sensitivity_ids: tuple[str, ...]
    variant_sensitivity_ids: Mapping[str, object]
    record_fingerprints: Mapping[str, object]
    changed_variant_ids: tuple[str, ...] = ()
    missing_variant_ids: tuple[str, ...] = ()
    failed_variant_ids: tuple[str, ...] = ()
    noncredible_variant_ids: tuple[str, ...] = ()
    caution_variant_ids: tuple[str, ...] = ()
    schema_version: str = SENSITIVITY_DECISION_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != SENSITIVITY_DECISION_SCHEMA_VERSION:
            raise ContractVersionError(
                f"Unsupported sensitivity decision version: {self.schema_version!r}"
            )
        for name in (
            "decision_id",
            "plan_id",
            "baseline_analysis_id",
            "conclusion_id",
            "summary",
        ):
            value = _require_text(getattr(self, name), name, token=name != "summary")
            object.__setattr__(self, name, value)
        object.__setattr__(
            self,
            "stability_state",
            _coerce_enum(StabilityState, self.stability_state, "stability state"),
        )
        object.__setattr__(self, "reason_code", _coerce_reason_code(self.reason_code))
        required_variants = _freeze_mapping(
            self.required_variant_analysis_ids,
            "required_variant_analysis_ids",
        )
        if not required_variants:
            raise ContractValidationError("A sensitivity decision requires planned variants")
        for variant_id, analysis_id in required_variants.items():
            _require_text(variant_id, "required variant ID", token=True)
            if not isinstance(analysis_id, str):
                raise ContractValidationError(
                    "required variant AnalysisIDs must be strings"
                )
            _require_text(analysis_id, "required variant AnalysisID", token=True)
        if len(set(required_variants.values())) != len(required_variants):
            raise ContractValidationError(
                "Required variants must map one-to-one to distinct AnalysisIDs"
            )
        object.__setattr__(self, "required_variant_analysis_ids", required_variants)
        object.__setattr__(self, "sensitivity_ids", tuple(sorted(_text_tuple(self.sensitivity_ids, "sensitivity ID"))))
        object.__setattr__(self, "changed_variant_ids", tuple(sorted(_text_tuple(self.changed_variant_ids, "variant ID"))))
        object.__setattr__(self, "missing_variant_ids", tuple(sorted(_text_tuple(self.missing_variant_ids, "variant ID"))))
        object.__setattr__(self, "failed_variant_ids", tuple(sorted(_text_tuple(self.failed_variant_ids, "variant ID"))))
        object.__setattr__(self, "noncredible_variant_ids", tuple(sorted(_text_tuple(self.noncredible_variant_ids, "variant ID"))))
        object.__setattr__(self, "caution_variant_ids", tuple(sorted(_text_tuple(self.caution_variant_ids, "variant ID"))))
        object.__setattr__(
            self,
            "variant_sensitivity_ids",
            _freeze_mapping(self.variant_sensitivity_ids, "variant_sensitivity_ids"),
        )
        object.__setattr__(
            self,
            "record_fingerprints",
            _freeze_mapping(self.record_fingerprints, "record_fingerprints"),
        )
        for sensitivity_id, fingerprint in self.record_fingerprints.items():
            _require_text(sensitivity_id, "record fingerprint SensitivityID", token=True)
            if not isinstance(fingerprint, str) or not re.fullmatch(r"[0-9a-f]{64}", fingerprint):
                raise ContractValidationError(
                    "record_fingerprints values must be full lowercase SHA-256 digests"
                )
        if set(self.record_fingerprints) != set(self.sensitivity_ids):
            raise ContractValidationError(
                "record_fingerprints must contain exactly the linked SensitivityIDs"
            )
        planned = set(self.required_variant_analysis_ids)
        variant_links = set(self.variant_sensitivity_ids)
        if not variant_links.issubset(planned):
            raise ContractValidationError(
                "variant_sensitivity_ids keys must be required variants"
            )
        linked_ids = tuple(self.variant_sensitivity_ids.values())
        if any(not isinstance(value, str) for value in linked_ids):
            raise ContractValidationError("variant_sensitivity_ids values must be SensitivityIDs")
        if len(linked_ids) != len(set(linked_ids)):
            raise ContractValidationError(
                "Each planned variant must link to a distinct SensitivityID"
            )
        if set(linked_ids) != set(self.sensitivity_ids):
            raise ContractValidationError(
                "variant_sensitivity_ids values must match sensitivity_ids"
            )
        if variant_links.intersection(self.missing_variant_ids):
            raise ContractValidationError("A variant cannot be both linked and missing")
        if variant_links.union(self.missing_variant_ids) != planned:
            raise ContractValidationError(
                "Linked and missing variants must partition the sensitivity plan"
            )
        for name in (
            "changed_variant_ids",
            "missing_variant_ids",
            "failed_variant_ids",
            "noncredible_variant_ids",
            "caution_variant_ids",
        ):
            if not set(getattr(self, name)).issubset(planned):
                raise ContractValidationError(f"{name} must be a subset of required variants")
        for name in (
            "changed_variant_ids",
            "failed_variant_ids",
            "noncredible_variant_ids",
            "caution_variant_ids",
        ):
            if not set(getattr(self, name)).issubset(variant_links):
                raise ContractValidationError(f"{name} must refer to linked variants")
        unusable = set(self.failed_variant_ids).union(self.noncredible_variant_ids)
        changed = set(self.changed_variant_ids)
        if changed.intersection(unusable):
            raise ContractValidationError(
                "Changed variants must be successfully evaluated and credible"
            )
        usable = variant_links.difference(unusable)
        if self.stability_state is StabilityState.STABLE:
            if self.reason_code != ReasonCode.SENSITIVITY_STABLE.value:
                raise ContractValidationError("STABLE requires sensitivity.stable")
            if changed or any((self.missing_variant_ids, self.failed_variant_ids, self.noncredible_variant_ids, self.caution_variant_ids)):
                raise ContractValidationError("STABLE cannot contain incomplete or caution variants")
        elif self.stability_state is StabilityState.SENSITIVE:
            if not changed or self.reason_code != ReasonCode.SENSITIVITY_CONCLUSION_CHANGED.value:
                raise ContractValidationError(
                    "SENSITIVE requires changed variants and sensitivity.conclusion_changed"
                )
        elif self.stability_state is StabilityState.CONDITIONALLY_STABLE:
            allowed = {
                ReasonCode.SENSITIVITY_MISSING_REQUIRED.value,
                ReasonCode.SENSITIVITY_THRESHOLD_NEAR.value,
                ReasonCode.SENSITIVITY_RULE_NOT_EVALUABLE.value,
                ReasonCode.LIMITED_EVIDENCE.value,
            }
            if changed or not usable or self.reason_code not in allowed:
                raise ContractValidationError(
                    "CONDITIONALLY_STABLE requires an unchanged usable result plus a limitation"
                )
            if not any((self.missing_variant_ids, self.failed_variant_ids, self.noncredible_variant_ids, self.caution_variant_ids)):
                raise ContractValidationError(
                    "CONDITIONALLY_STABLE requires at least one recorded limitation"
                )
        else:
            allowed = {
                ReasonCode.SENSITIVITY_NOT_RUN.value,
                ReasonCode.SENSITIVITY_MISSING_REQUIRED.value,
                ReasonCode.SENSITIVITY_RULE_NOT_EVALUABLE.value,
            }
            if changed or usable or self.reason_code not in allowed:
                raise ContractValidationError(
                    "NOT_ASSESSED requires no usable evaluated sensitivity result"
                )
        if self.stability_state is not StabilityState.NOT_ASSESSED and not self.sensitivity_ids:
            raise ContractValidationError("Assessed stability requires SensitivityIDs")
        expected_id = "sdc_" + payload_fingerprint(
            {
                "schema_version": self.schema_version,
                "reason_registry_version": REASON_REGISTRY_VERSION,
                "plan_id": self.plan_id,
                "baseline_analysis_id": self.baseline_analysis_id,
                "conclusion_id": self.conclusion_id,
                "stability_state": self.stability_state,
                "reason_code": self.reason_code,
                "summary": self.summary,
                "required_variant_analysis_ids": self.required_variant_analysis_ids,
                "sensitivity_ids": self.sensitivity_ids,
                "variant_sensitivity_ids": self.variant_sensitivity_ids,
                "record_fingerprints": self.record_fingerprints,
                "changed_variant_ids": self.changed_variant_ids,
                "missing_variant_ids": self.missing_variant_ids,
                "failed_variant_ids": self.failed_variant_ids,
                "noncredible_variant_ids": self.noncredible_variant_ids,
                "caution_variant_ids": self.caution_variant_ids,
            },
            length=24,
        )
        if self.decision_id != expected_id:
            raise ContractValidationError(
                f"SensitivityDecisionID does not match the contract payload; expected {expected_id!r}"
            )

    @property
    def required_variant_ids(self) -> tuple[str, ...]:
        return tuple(self.required_variant_analysis_ids)

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "decision_id": self.decision_id,
            "plan_id": self.plan_id,
            "baseline_analysis_id": self.baseline_analysis_id,
            "conclusion_id": self.conclusion_id,
            "stability_state": self.stability_state.value,
            "reason_code": self.reason_code,
            "reason_registry_version": REASON_REGISTRY_VERSION,
            "summary": self.summary,
            "required_variant_analysis_ids": _thaw_json(
                self.required_variant_analysis_ids
            ),
            "sensitivity_ids": list(self.sensitivity_ids),
            "variant_sensitivity_ids": _thaw_json(self.variant_sensitivity_ids),
            "record_fingerprints": _thaw_json(self.record_fingerprints),
            "changed_variant_ids": list(self.changed_variant_ids),
            "missing_variant_ids": list(self.missing_variant_ids),
            "failed_variant_ids": list(self.failed_variant_ids),
            "noncredible_variant_ids": list(self.noncredible_variant_ids),
            "caution_variant_ids": list(self.caution_variant_ids),
        }

    def to_row(self) -> dict[str, object]:
        payload = self.to_payload()
        return {
            "SchemaVersion": self.schema_version,
            "SensitivityDecisionID": self.decision_id,
            "SensitivityPlanID": self.plan_id,
            "BaselineAnalysisID": self.baseline_analysis_id,
            "ConclusionID": self.conclusion_id,
            "StabilityState": self.stability_state.value,
            "ReasonCode": self.reason_code,
            "ReasonRegistryVersion": REASON_REGISTRY_VERSION,
            "Summary": self.summary,
            "RequiredVariantIDsJSON": canonical_json(self.required_variant_ids),
            "RequiredVariantAnalysisIDsJSON": canonical_json(
                payload["required_variant_analysis_ids"]
            ),
            "SensitivityIDsJSON": canonical_json(payload["sensitivity_ids"]),
            "VariantSensitivityIDsJSON": canonical_json(payload["variant_sensitivity_ids"]),
            "RecordFingerprintsJSON": canonical_json(payload["record_fingerprints"]),
            "ChangedVariantIDsJSON": canonical_json(payload["changed_variant_ids"]),
            "MissingVariantIDsJSON": canonical_json(payload["missing_variant_ids"]),
            "FailedVariantIDsJSON": canonical_json(payload["failed_variant_ids"]),
            "NoncredibleVariantIDsJSON": canonical_json(payload["noncredible_variant_ids"]),
            "CautionVariantIDsJSON": canonical_json(payload["caution_variant_ids"]),
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
        *,
        plan: SensitivityPlan | None = None,
        records: Iterable[SensitivityRecord] | object = _MISSING_RULE_VALUE,
    ) -> "SensitivityDecision":
        """Restore only when the linked plan and records reproduce the payload."""
        if plan is None or records is _MISSING_RULE_VALUE:
            raise ContractValidationError(
                "SensitivityDecision.from_payload requires linked plan and records"
            )
        _validate_payload_shape(
            payload,
            name="SensitivityDecision",
            keys={
                "schema_version",
                "decision_id",
                "plan_id",
                "baseline_analysis_id",
                "conclusion_id",
                "stability_state",
                "reason_code",
                "reason_registry_version",
                "summary",
                "required_variant_analysis_ids",
                "sensitivity_ids",
                "variant_sensitivity_ids",
                "record_fingerprints",
                "changed_variant_ids",
                "missing_variant_ids",
                "failed_variant_ids",
                "noncredible_variant_ids",
                "caution_variant_ids",
            },
            reason_registry=True,
        )
        candidate = cls(
            schema_version=payload["schema_version"],
            decision_id=payload["decision_id"],
            plan_id=payload["plan_id"],
            baseline_analysis_id=payload["baseline_analysis_id"],
            conclusion_id=payload["conclusion_id"],
            stability_state=payload["stability_state"],
            reason_code=payload["reason_code"],
            summary=payload["summary"],
            required_variant_analysis_ids=payload["required_variant_analysis_ids"],
            sensitivity_ids=_payload_sequence(payload, "sensitivity_ids"),
            variant_sensitivity_ids=payload["variant_sensitivity_ids"],
            record_fingerprints=payload["record_fingerprints"],
            changed_variant_ids=_payload_sequence(payload, "changed_variant_ids"),
            missing_variant_ids=_payload_sequence(payload, "missing_variant_ids"),
            failed_variant_ids=_payload_sequence(payload, "failed_variant_ids"),
            noncredible_variant_ids=_payload_sequence(
                payload, "noncredible_variant_ids"
            ),
            caution_variant_ids=_payload_sequence(payload, "caution_variant_ids"),
        )
        expected = synthesize_sensitivity_decision(tuple(records), plan=plan)
        if candidate != expected:
            raise ContractValidationError(
                "Sensitivity decision payload does not reproduce from its linked plan and records"
            )
        return candidate


def synthesize_sensitivity_decision(
    records: Iterable[SensitivityRecord],
    *,
    plan: SensitivityPlan,
) -> SensitivityDecision:
    """Aggregate sensitivity records using an order-independent fail-closed rule."""
    if not isinstance(plan, SensitivityPlan):
        raise ContractValidationError("plan must be a SensitivityPlan")
    all_records = tuple(records)
    if any(not isinstance(record, SensitivityRecord) for record in all_records):
        raise ContractValidationError("All sensitivity entries must be SensitivityRecord instances")
    if any(record.conclusion_id != plan.conclusion_id for record in all_records):
        raise ContractValidationError("Sensitivity records do not match the plan ConclusionID")
    if any(record.baseline_analysis_id != plan.baseline_analysis_id for record in all_records):
        raise ContractValidationError("Sensitivity records do not match the plan baseline AnalysisID")
    by_variant = {record.variant_id: record for record in all_records}
    if len(by_variant) != len(all_records):
        raise ContractValidationError("A conclusion may have only one sensitivity record per VariantID")
    expected = set(plan.required_variant_analysis_ids)
    unplanned_required = sorted(
        record.variant_id
        for record in all_records
        if record.variant_id not in expected and record.required
    )
    if unplanned_required:
        raise ContractValidationError(
            f"Required sensitivity records are absent from the plan: {unplanned_required!r}"
        )
    required_candidates = tuple(by_variant[variant] for variant in plan.required_variant_ids if variant in by_variant)
    for record in required_candidates:
        if not record.required:
            raise ContractValidationError(
                f"Planned variant {record.variant_id!r} must be marked required"
            )
        expected_analysis_id = plan.required_variant_analysis_ids[record.variant_id]
        if record.variant_analysis_id != expected_analysis_id:
            raise ContractValidationError(
                f"Variant {record.variant_id!r} does not match its prespecified AnalysisID"
            )
    missing = tuple(sorted(expected.difference(by_variant)))
    sensitivity_ids = tuple(sorted(record.sensitivity_id for record in required_candidates))
    record_fingerprints = {
        record.sensitivity_id: payload_fingerprint(record.to_payload())
        for record in required_candidates
    }
    variant_sensitivity_ids = {
        record.variant_id: record.sensitivity_id
        for record in required_candidates
    }
    noncredible = tuple(sorted(record.variant_id for record in required_candidates if not record.credible))
    rule_outcomes = {
        record.variant_id: (
            _evaluate_sensitivity_rule(plan.decision_rule, record.metrics)
            if record.computation_state in {ComputationState.AVAILABLE, ComputationState.CAUTION}
            and record.comparable
            else None
        )
        for record in required_candidates
    }
    rule_unevaluable = tuple(sorted(
        record.variant_id
        for record in required_candidates
        if record.computation_state in {ComputationState.AVAILABLE, ComputationState.CAUTION}
        and record.comparable
        and rule_outcomes[record.variant_id] is None
    ))
    failed = tuple(sorted(
        record.variant_id
        for record in required_candidates
        if record.computation_state not in {ComputationState.AVAILABLE, ComputationState.CAUTION}
        or not record.comparable
        or rule_outcomes[record.variant_id] is None
    ))
    cautions = tuple(sorted(
        record.variant_id
        for record in required_candidates
        if record.computation_state is ComputationState.CAUTION
    ))
    completed = tuple(
        record
        for record in required_candidates
        if record.credible
        if record.computation_state in {ComputationState.AVAILABLE, ComputationState.CAUTION}
        and record.comparable
        and rule_outcomes[record.variant_id] is not None
    )
    changed = tuple(sorted(
        record.variant_id
        for record in completed
        if rule_outcomes[record.variant_id] is True
    ))

    def decision(
        state: StabilityState,
        reason_code: ReasonCode | str,
        summary: str,
    ) -> SensitivityDecision:
        decision_id = "sdc_" + payload_fingerprint(
            {
                "schema_version": SENSITIVITY_DECISION_SCHEMA_VERSION,
                "reason_registry_version": REASON_REGISTRY_VERSION,
                "plan_id": plan.plan_id,
                "baseline_analysis_id": plan.baseline_analysis_id,
                "conclusion_id": plan.conclusion_id,
                "stability_state": state,
                "reason_code": _coerce_reason_code(reason_code),
                "summary": summary,
                "required_variant_analysis_ids": plan.required_variant_analysis_ids,
                "sensitivity_ids": sensitivity_ids,
                "variant_sensitivity_ids": variant_sensitivity_ids,
                "record_fingerprints": record_fingerprints,
                "changed_variant_ids": changed,
                "missing_variant_ids": missing,
                "failed_variant_ids": failed,
                "noncredible_variant_ids": noncredible,
                "caution_variant_ids": cautions,
            },
            length=24,
        )
        return SensitivityDecision(
            decision_id=decision_id,
            plan_id=plan.plan_id,
            baseline_analysis_id=plan.baseline_analysis_id,
            conclusion_id=plan.conclusion_id,
            stability_state=state,
            reason_code=_coerce_reason_code(reason_code),
            summary=summary,
            required_variant_analysis_ids=plan.required_variant_analysis_ids,
            sensitivity_ids=sensitivity_ids,
            variant_sensitivity_ids=variant_sensitivity_ids,
            record_fingerprints=record_fingerprints,
            changed_variant_ids=changed,
            missing_variant_ids=missing,
            failed_variant_ids=failed,
            noncredible_variant_ids=noncredible,
            caution_variant_ids=cautions,
        )

    if changed:
        return decision(
            StabilityState.SENSITIVE,
            ReasonCode.SENSITIVITY_CONCLUSION_CHANGED,
            f"{len(changed)} required credible variant(s) changed the conclusion.",
        )

    incomplete = bool(missing or failed or noncredible)
    if not completed:
        reason = ReasonCode.SENSITIVITY_MISSING_REQUIRED
        if rule_unevaluable and not missing and set(failed) == set(rule_unevaluable):
            reason = ReasonCode.SENSITIVITY_RULE_NOT_EVALUABLE
        return decision(
            StabilityState.NOT_ASSESSED,
            reason,
            "No required credible sensitivity variant was successfully compared.",
        )

    if incomplete or cautions:
        caution_reasons = {
            record.reason_code
            for record in completed
            if record.computation_state is ComputationState.CAUTION
        }
        return decision(
            StabilityState.CONDITIONALLY_STABLE,
            (
                ReasonCode.SENSITIVITY_RULE_NOT_EVALUABLE
                if incomplete and rule_unevaluable and not missing
                else ReasonCode.SENSITIVITY_MISSING_REQUIRED
                if incomplete
                else (
                    ReasonCode.SENSITIVITY_THRESHOLD_NEAR
                    if ReasonCode.SENSITIVITY_THRESHOLD_NEAR.value in caution_reasons
                    else ReasonCode.LIMITED_EVIDENCE
                )
            ),
            (
                f"The conclusion was unchanged in {len(completed)} completed variant(s), "
                "but the required set was incomplete or limited."
            ),
        )

    return decision(
        StabilityState.STABLE,
        ReasonCode.SENSITIVITY_STABLE,
        f"The conclusion was unchanged in all {len(completed)} required credible variant(s).",
    )


@dataclass(frozen=True)
class DecisionRecord:
    """Evidence-linked workflow disposition for one user-facing decision."""

    decision_id: str
    decision_key: str
    analysis_id: str
    conclusion_id: str
    title: str
    computation_state: ComputationState
    stability_state: StabilityState
    disposition: DecisionDisposition
    rationale: str
    interpretation_boundary: str
    recommended_action: str
    evidence_ids: tuple[str, ...]
    required_evidence_ids: tuple[str, ...]
    evidence_fingerprints: Mapping[str, object]
    boundary_only: bool = False
    sensitivity_ids: tuple[str, ...] = ()
    sensitivity_decision_id: str | None = None
    sensitivity_plan_id: str | None = None
    sensitivity_decision_fingerprint: str | None = None
    reason_codes: tuple[str, ...] = ()
    schema_version: str = DECISION_RECORD_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.schema_version != DECISION_RECORD_SCHEMA_VERSION:
            raise ContractVersionError(
                f"Unsupported decision record version: {self.schema_version!r}"
            )
        for name in (
            "decision_id",
            "decision_key",
            "analysis_id",
            "conclusion_id",
            "title",
            "rationale",
            "interpretation_boundary",
            "recommended_action",
        ):
            value = _require_text(
                getattr(self, name),
                name,
                token=name in {"decision_id", "decision_key", "analysis_id", "conclusion_id"},
            )
            object.__setattr__(self, name, value)
        object.__setattr__(
            self,
            "computation_state",
            _coerce_enum(ComputationState, self.computation_state, "computation state"),
        )
        object.__setattr__(
            self,
            "stability_state",
            _coerce_enum(StabilityState, self.stability_state, "stability state"),
        )
        object.__setattr__(
            self,
            "disposition",
            _coerce_enum(DecisionDisposition, self.disposition, "decision disposition"),
        )
        object.__setattr__(self, "evidence_ids", _text_tuple(self.evidence_ids, "evidence ID"))
        object.__setattr__(
            self,
            "required_evidence_ids",
            _text_tuple(self.required_evidence_ids, "required evidence ID"),
        )
        object.__setattr__(
            self,
            "evidence_fingerprints",
            _freeze_mapping(self.evidence_fingerprints, "evidence_fingerprints"),
        )
        if not isinstance(self.boundary_only, bool):
            raise ContractValidationError("boundary_only must be a boolean")
        object.__setattr__(self, "sensitivity_ids", _text_tuple(self.sensitivity_ids, "sensitivity ID"))
        for name in ("sensitivity_decision_id", "sensitivity_plan_id"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _require_text(value, name, token=True))
        object.__setattr__(
            self,
            "reason_codes",
            tuple(_coerce_reason_code(code) for code in self.reason_codes),
        )
        if not self.evidence_ids:
            raise ContractValidationError("A decision record requires at least one EvidenceID")
        if not self.required_evidence_ids:
            raise ContractValidationError("A decision record requires at least one required EvidenceID")
        if not set(self.required_evidence_ids).issubset(self.evidence_ids):
            raise ContractValidationError("required_evidence_ids must be a subset of evidence_ids")
        if set(self.evidence_fingerprints) != set(self.evidence_ids):
            raise ContractValidationError(
                "evidence_fingerprints must contain exactly the linked EvidenceIDs"
            )
        for evidence_id, fingerprint in self.evidence_fingerprints.items():
            _require_text(evidence_id, "evidence fingerprint EvidenceID", token=True)
            if not isinstance(fingerprint, str) or not re.fullmatch(r"[0-9a-f]{64}", fingerprint):
                raise ContractValidationError(
                    "evidence_fingerprints values must be full lowercase SHA-256 digests"
                )
        if len(self.reason_codes) != len(set(self.reason_codes)):
            raise ContractValidationError("reason_codes must be unique")
        has_sensitivity_decision = self.sensitivity_decision_id is not None
        has_sensitivity_plan = self.sensitivity_plan_id is not None
        if has_sensitivity_decision != has_sensitivity_plan:
            raise ContractValidationError(
                "SensitivityDecisionID and SensitivityPlanID must be provided together"
            )
        if has_sensitivity_decision:
            if (
                not isinstance(self.sensitivity_decision_fingerprint, str)
                or not re.fullmatch(
                    r"[0-9a-f]{64}", self.sensitivity_decision_fingerprint
                )
            ):
                raise ContractValidationError(
                    "A linked sensitivity decision requires its full payload fingerprint"
                )
        elif self.sensitivity_decision_fingerprint is not None:
            raise ContractValidationError(
                "A sensitivity fingerprint requires a SensitivityDecisionID"
            )
        if self.stability_state is not StabilityState.NOT_ASSESSED and not has_sensitivity_decision:
            raise ContractValidationError(
                "Assessed stability requires SensitivityDecisionID and SensitivityPlanID"
            )
        if self.stability_state is not StabilityState.NOT_ASSESSED and not self.sensitivity_ids:
            raise ContractValidationError("Assessed stability requires SensitivityIDs")
        if not has_sensitivity_decision and self.sensitivity_ids:
            raise ContractValidationError(
                "SensitivityIDs require a SensitivityDecisionID and SensitivityPlanID"
            )
        if not any(code.startswith("sensitivity.") for code in self.reason_codes):
            raise ContractValidationError(
                "A decision record requires an explicit sensitivity reason code"
            )
        if (
            not has_sensitivity_decision
            and ReasonCode.SENSITIVITY_NOT_RUN.value not in self.reason_codes
        ):
            raise ContractValidationError(
                "A decision without a sensitivity plan requires sensitivity.not_run"
            )
        expected_disposition = derive_decision_disposition(
            self.computation_state,
            self.stability_state,
            boundary_only=self.boundary_only,
        )
        if self.disposition is not expected_disposition:
            raise ContractValidationError(
                f"Decision disposition is inconsistent; expected {expected_disposition.value}"
            )
        expected_id = "dec_" + payload_fingerprint(
            {
                "schema_version": self.schema_version,
                "reason_registry_version": REASON_REGISTRY_VERSION,
                "analysis_id": self.analysis_id,
                "decision_key": self.decision_key,
                "conclusion_id": self.conclusion_id,
                "title": self.title,
                "computation_state": self.computation_state,
                "stability_state": self.stability_state,
                "disposition": self.disposition,
                "reason_codes": self.reason_codes,
                "rationale": self.rationale,
                "interpretation_boundary": self.interpretation_boundary,
                "recommended_action": self.recommended_action,
                "evidence_ids": self.evidence_ids,
                "required_evidence_ids": self.required_evidence_ids,
                "evidence_fingerprints": self.evidence_fingerprints,
                "boundary_only": self.boundary_only,
                "sensitivity_ids": self.sensitivity_ids,
                "sensitivity_decision_id": self.sensitivity_decision_id,
                "sensitivity_plan_id": self.sensitivity_plan_id,
                "sensitivity_decision_fingerprint": self.sensitivity_decision_fingerprint,
            },
            length=24,
        )
        if self.decision_id != expected_id:
            raise ContractValidationError(
                f"DecisionID does not match the contract payload; expected {expected_id!r}"
            )

    def to_payload(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "decision_id": self.decision_id,
            "decision_key": self.decision_key,
            "analysis_id": self.analysis_id,
            "conclusion_id": self.conclusion_id,
            "title": self.title,
            "computation_state": self.computation_state.value,
            "stability_state": self.stability_state.value,
            "disposition": self.disposition.value,
            "reason_codes": list(self.reason_codes),
            "reason_registry_version": REASON_REGISTRY_VERSION,
            "rationale": self.rationale,
            "interpretation_boundary": self.interpretation_boundary,
            "recommended_action": self.recommended_action,
            "evidence_ids": list(self.evidence_ids),
            "required_evidence_ids": list(self.required_evidence_ids),
            "evidence_fingerprints": _thaw_json(self.evidence_fingerprints),
            "boundary_only": self.boundary_only,
            "sensitivity_ids": list(self.sensitivity_ids),
            "sensitivity_decision_id": self.sensitivity_decision_id,
            "sensitivity_plan_id": self.sensitivity_plan_id,
            "sensitivity_decision_fingerprint": self.sensitivity_decision_fingerprint,
        }

    def to_row(self) -> dict[str, object]:
        payload = self.to_payload()
        return {
            "SchemaVersion": self.schema_version,
            "DecisionID": self.decision_id,
            "DecisionKey": self.decision_key,
            "AnalysisID": self.analysis_id,
            "ConclusionID": self.conclusion_id,
            "Title": self.title,
            "ComputationState": self.computation_state.value,
            "StabilityState": self.stability_state.value,
            "Disposition": self.disposition.value,
            "ReasonCodesJSON": canonical_json(payload["reason_codes"]),
            "ReasonRegistryVersion": REASON_REGISTRY_VERSION,
            "Rationale": self.rationale,
            "InterpretationBoundary": self.interpretation_boundary,
            "RecommendedAction": self.recommended_action,
            "EvidenceIDsJSON": canonical_json(payload["evidence_ids"]),
            "RequiredEvidenceIDsJSON": canonical_json(payload["required_evidence_ids"]),
            "EvidenceFingerprintsJSON": canonical_json(payload["evidence_fingerprints"]),
            "BoundaryOnly": self.boundary_only,
            "SensitivityIDsJSON": canonical_json(payload["sensitivity_ids"]),
            "SensitivityDecisionID": self.sensitivity_decision_id or "",
            "SensitivityPlanID": self.sensitivity_plan_id or "",
            "SensitivityDecisionFingerprint": self.sensitivity_decision_fingerprint or "",
        }

    @classmethod
    def from_payload(
        cls,
        payload: Mapping[str, object],
        *,
        identity: AnalysisIdentity | None = None,
        evidence: Iterable[EvidenceRecord] | object = _MISSING_RULE_VALUE,
        sensitivity_plan: SensitivityPlan | None = None,
        sensitivity_records: Iterable[SensitivityRecord] | object = _MISSING_RULE_VALUE,
    ) -> "DecisionRecord":
        """Restore only when all linked evidence and sensitivity inputs reproduce it."""
        if identity is None or evidence is _MISSING_RULE_VALUE:
            raise ContractValidationError(
                "DecisionRecord.from_payload requires identity and linked evidence"
            )
        _validate_payload_shape(
            payload,
            name="DecisionRecord",
            keys={
                "schema_version",
                "decision_id",
                "decision_key",
                "analysis_id",
                "conclusion_id",
                "title",
                "computation_state",
                "stability_state",
                "disposition",
                "reason_codes",
                "reason_registry_version",
                "rationale",
                "interpretation_boundary",
                "recommended_action",
                "evidence_ids",
                "required_evidence_ids",
                "evidence_fingerprints",
                "boundary_only",
                "sensitivity_ids",
                "sensitivity_decision_id",
                "sensitivity_plan_id",
                "sensitivity_decision_fingerprint",
            },
            reason_registry=True,
        )
        candidate = cls(
            schema_version=payload["schema_version"],
            decision_id=payload["decision_id"],
            decision_key=payload["decision_key"],
            analysis_id=payload["analysis_id"],
            conclusion_id=payload["conclusion_id"],
            title=payload["title"],
            computation_state=payload["computation_state"],
            stability_state=payload["stability_state"],
            disposition=payload["disposition"],
            reason_codes=_payload_sequence(payload, "reason_codes"),
            rationale=payload["rationale"],
            interpretation_boundary=payload["interpretation_boundary"],
            recommended_action=payload["recommended_action"],
            evidence_ids=_payload_sequence(payload, "evidence_ids"),
            required_evidence_ids=_payload_sequence(payload, "required_evidence_ids"),
            evidence_fingerprints=payload["evidence_fingerprints"],
            boundary_only=payload["boundary_only"],
            sensitivity_ids=_payload_sequence(payload, "sensitivity_ids"),
            sensitivity_decision_id=payload["sensitivity_decision_id"],
            sensitivity_plan_id=payload["sensitivity_plan_id"],
            sensitivity_decision_fingerprint=payload[
                "sensitivity_decision_fingerprint"
            ],
        )
        records = (
            ()
            if sensitivity_records is _MISSING_RULE_VALUE
            else tuple(sensitivity_records)
        )
        validate_decision_bundle(
            candidate,
            identity=identity,
            evidence=tuple(evidence),
            sensitivity_plan=sensitivity_plan,
            sensitivity_records=records,
        )
        return candidate


def _aggregate_computation_state(records: Sequence[EvidenceRecord]) -> ComputationState:
    required = tuple(record for record in records if record.required)
    if not required:
        raise ContractValidationError("A decision requires at least one required evidence record")
    states = {record.computation_state for record in required}
    if ComputationState.HOLD in states:
        return ComputationState.HOLD
    if ComputationState.NOT_ASSESSABLE in states:
        return ComputationState.NOT_ASSESSABLE
    if ComputationState.CAUTION in states:
        return ComputationState.CAUTION
    return ComputationState.AVAILABLE


def derive_decision_disposition(
    computation_state: ComputationState | str,
    stability_state: StabilityState | str,
    *,
    boundary_only: bool = False,
) -> DecisionDisposition:
    """Route a decision without turning the route into a validity verdict."""
    computation = _coerce_enum(ComputationState, computation_state, "computation state")
    stability = _coerce_enum(StabilityState, stability_state, "stability state")
    if computation is ComputationState.HOLD:
        return DecisionDisposition.WITHHOLD
    if computation is ComputationState.NOT_ASSESSABLE:
        return DecisionDisposition.NOT_EVALUATED
    if boundary_only:
        return DecisionDisposition.BOUNDARY_ONLY
    if computation is ComputationState.AVAILABLE and stability is StabilityState.STABLE:
        return DecisionDisposition.USE
    return DecisionDisposition.USE_WITH_CAVEAT


def make_decision_record(
    identity: AnalysisIdentity,
    *,
    decision_key: str,
    conclusion_id: str,
    title: str,
    evidence: Iterable[EvidenceRecord],
    rationale: str,
    interpretation_boundary: str,
    recommended_action: str,
    sensitivity_plan: SensitivityPlan | None = None,
    sensitivity_records: Iterable[SensitivityRecord] = (),
    boundary_only: bool = False,
) -> DecisionRecord:
    """Create a decision linked to validated evidence and optional sensitivity."""
    if not isinstance(identity, AnalysisIdentity) or not identity.is_baseline:
        raise ContractValidationError("A decision must reference a baseline AnalysisIdentity")
    analysis_id = identity.analysis_id
    records = validate_evidence_records(evidence, identity=analysis_id)
    if not records:
        raise ContractValidationError("A decision requires at least one evidence record")
    conclusion_id = _require_text(conclusion_id, "conclusion_id", token=True)
    decision_key = _require_text(decision_key, "decision_key", token=True)
    title = _require_text(title, "title")
    rationale = _require_text(rationale, "rationale")
    interpretation_boundary = _require_text(
        interpretation_boundary, "interpretation_boundary"
    )
    recommended_action = _require_text(recommended_action, "recommended_action")
    if not isinstance(boundary_only, bool):
        raise ContractValidationError("boundary_only must be a boolean")
    sensitivity_inputs = tuple(sensitivity_records)
    linked_evidence_ids = {record.evidence_id for record in records}
    unresolved_sensitivity_evidence = sorted(
        {
            evidence_id
            for record in sensitivity_inputs
            for evidence_id in record.evidence_ids
            if evidence_id not in linked_evidence_ids
        }
    )
    if unresolved_sensitivity_evidence:
        raise ContractValidationError(
            "Sensitivity records reference EvidenceIDs outside the decision bundle: "
            f"{unresolved_sensitivity_evidence!r}"
        )
    if sensitivity_plan is not None:
        sensitivity_decision = synthesize_sensitivity_decision(
            sensitivity_inputs,
            plan=sensitivity_plan,
        )
        if sensitivity_decision.baseline_analysis_id != analysis_id:
            raise ContractValidationError("Sensitivity decision does not match the decision AnalysisID")
        if sensitivity_decision.conclusion_id != conclusion_id:
            raise ContractValidationError("Sensitivity decision does not match the decision ConclusionID")
        resolved_stability = sensitivity_decision.stability_state
        sensitivity_ids = sensitivity_decision.sensitivity_ids
        sensitivity_reason = (sensitivity_decision.reason_code,)
        sensitivity_decision_id = sensitivity_decision.decision_id
        sensitivity_plan_id = sensitivity_decision.plan_id
        sensitivity_decision_fingerprint = payload_fingerprint(
            sensitivity_decision.to_payload()
        )
    else:
        if sensitivity_inputs:
            raise ContractValidationError(
                "Sensitivity records require a matching SensitivityPlan"
            )
        resolved_stability = StabilityState.NOT_ASSESSED
        sensitivity_ids = ()
        sensitivity_reason = (ReasonCode.SENSITIVITY_NOT_RUN.value,)
        sensitivity_decision_id = None
        sensitivity_plan_id = None
        sensitivity_decision_fingerprint = None
    computation = _aggregate_computation_state(records)
    disposition = derive_decision_disposition(
        computation,
        resolved_stability,
        boundary_only=boundary_only,
    )
    reasons = tuple(
        sorted(
            {
                *(record.reason_code for record in records if record.reason_code is not None),
                *sensitivity_reason,
            },
            key=str,
        )
    )
    evidence_ids = tuple(sorted(record.evidence_id for record in records))
    required_evidence_ids = tuple(
        sorted(record.evidence_id for record in records if record.required)
    )
    evidence_fingerprints = {
        record.evidence_id: payload_fingerprint(record.to_payload())
        for record in sorted(records, key=lambda item: item.evidence_id)
    }
    identity_payload = {
        "schema_version": DECISION_RECORD_SCHEMA_VERSION,
        "reason_registry_version": REASON_REGISTRY_VERSION,
        "analysis_id": analysis_id,
        "decision_key": decision_key,
        "conclusion_id": conclusion_id,
        "title": title,
        "computation_state": computation,
        "stability_state": resolved_stability,
        "disposition": disposition,
        "reason_codes": reasons,
        "rationale": rationale,
        "interpretation_boundary": interpretation_boundary,
        "recommended_action": recommended_action,
        "evidence_ids": evidence_ids,
        "required_evidence_ids": required_evidence_ids,
        "evidence_fingerprints": evidence_fingerprints,
        "boundary_only": boundary_only,
        "sensitivity_ids": tuple(sorted(sensitivity_ids)),
        "sensitivity_decision_id": sensitivity_decision_id,
        "sensitivity_plan_id": sensitivity_plan_id,
        "sensitivity_decision_fingerprint": sensitivity_decision_fingerprint,
    }
    decision_id = "dec_" + payload_fingerprint(identity_payload, length=24)
    return DecisionRecord(
        decision_id=decision_id,
        decision_key=decision_key,
        analysis_id=analysis_id,
        conclusion_id=conclusion_id,
        title=title,
        computation_state=computation,
        stability_state=resolved_stability,
        disposition=disposition,
        reason_codes=reasons,
        rationale=rationale,
        interpretation_boundary=interpretation_boundary,
        recommended_action=recommended_action,
        evidence_ids=evidence_ids,
        required_evidence_ids=required_evidence_ids,
        evidence_fingerprints=evidence_fingerprints,
        boundary_only=boundary_only,
        sensitivity_ids=tuple(sorted(sensitivity_ids)),
        sensitivity_decision_id=sensitivity_decision_id,
        sensitivity_plan_id=sensitivity_plan_id,
        sensitivity_decision_fingerprint=sensitivity_decision_fingerprint,
    )


def validate_decision_bundle(
    decision: DecisionRecord,
    *,
    identity: AnalysisIdentity,
    evidence: Iterable[EvidenceRecord],
    sensitivity_plan: SensitivityPlan | None = None,
    sensitivity_records: Iterable[SensitivityRecord] = (),
) -> DecisionRecord:
    """Recompute a decision from its complete linked bundle and require equality."""
    if not isinstance(decision, DecisionRecord):
        raise ContractValidationError("decision must be a DecisionRecord")
    expected = make_decision_record(
        identity,
        decision_key=decision.decision_key,
        conclusion_id=decision.conclusion_id,
        title=decision.title,
        evidence=tuple(evidence),
        rationale=decision.rationale,
        interpretation_boundary=decision.interpretation_boundary,
        recommended_action=decision.recommended_action,
        sensitivity_plan=sensitivity_plan,
        sensitivity_records=tuple(sensitivity_records),
        boundary_only=decision.boundary_only,
    )
    if decision != expected:
        raise ContractValidationError(
            "Decision payload does not reproduce from its linked evidence and sensitivity bundle"
        )
    return decision
