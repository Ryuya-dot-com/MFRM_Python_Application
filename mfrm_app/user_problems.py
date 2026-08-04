"""Privacy-safe user-problem contracts for the standalone Python app.

The public application needs to explain a failed step without copying an
exception, a path, an uploaded value, or an internal identifier into the UI.
This module keeps that boundary structural:

* a :class:`UserProblemSpec` contains reviewed, language-independent keys;
* a :class:`UserProblemNotice` stores only a registered problem code, a
  bounded occurrence phase, and a random support reference; and
* exception text is inspected transiently for coarse classification and is
  never retained by either contract or its serialized payload.

The module intentionally has no Streamlit, pandas, estimator, filesystem, or
environment dependency.  Rendering and confidential exception logging remain
adapter responsibilities.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import json
import re
import secrets
from types import MappingProxyType
from typing import Mapping


USER_PROBLEM_SCHEMA_VERSION = "mfrm_user_problem_v1"
USER_PROBLEM_SUPPORT_REFERENCE_LABEL_KEY = "problems.support_reference_label"
_SUPPORT_TOKEN_BYTES = 12
_SUPPORT_REFERENCE_RE = re.compile(
    r"^MFRM-[0-9A-F]{8}-[0-9A-F]{8}-[0-9A-F]{8}$"
)
_STABLE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:[._:-][a-z0-9][a-z0-9_-]*)*$")
_LOCALE_KEY_RE = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z0-9_]+)+$")


class UserProblemContractError(ValueError):
    """Raised when a user-problem contract violates its safe schema."""


class UserProblemSeverity(str, Enum):
    """Presentation priority without exposing an exception category."""

    INFORMATION = "information"
    CAUTION = "caution"
    BLOCKED = "blocked"


class UserProblemPhase(str, Enum):
    """Coarse application phase used to select a bounded classifier."""

    PARSE = "parse"
    ESTIMATION = "estimation"
    RESIDUAL_PCA = "residual_pca"
    RESOURCE = "resource"
    APPLICATION = "application"


class UserProblemCode(str, Enum):
    """Stable public problem codes; values are not exception class names."""

    INPUT_PARSE_FAILED = "input.parse_failed"
    INPUT_MAPPING_INVALID = "input.mapping_invalid"
    INPUT_NO_USABLE_ROWS = "input.no_usable_rows"
    ESTIMATION_IDENTIFICATION_FAILED = "estimation.identification_failed"
    ESTIMATION_NONCONVERGENCE = "estimation.nonconvergence"
    ESTIMATION_RATING_SCALE_INVALID = "estimation.rating_scale_invalid"
    ESTIMATION_ANCHOR_INVALID = "estimation.anchor_invalid"
    RESOURCE_HOSTED_LIMIT = "resource.hosted_limit"
    DIAGNOSTIC_RESIDUAL_PCA_UNAVAILABLE = "diagnostic.residual_pca_unavailable"
    DIAGNOSTIC_RESIDUAL_PCA_FAILED = "diagnostic.residual_pca_failed"
    PROBLEM_UNEXPECTED = "problem.unexpected"


# Concise aliases make the enums convenient at adapter call sites while the
# explicit names above remain the canonical public API.
ProblemSeverity = UserProblemSeverity
ProblemPhase = UserProblemPhase
ProblemCode = UserProblemCode


def _coerce_enum(enum_type: type[Enum], value: object, field_name: str) -> Enum:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(value)
    except (TypeError, ValueError) as exc:
        # Do not interpolate the rejected value: callers may have received it
        # from an unsafe exception or request payload.
        raise UserProblemContractError(f"Unsupported {field_name}") from exc


def _require_stable_id(value: object, field_name: str, *, prefix: str | None = None) -> str:
    if not isinstance(value, str) or value != value.strip() or not _STABLE_ID_RE.fullmatch(value):
        raise UserProblemContractError(
            f"{field_name} must be a lowercase language-independent stable ID"
        )
    if prefix is not None and not value.startswith(prefix):
        raise UserProblemContractError(f"{field_name} must use the {prefix} namespace")
    return value


def _require_locale_key(value: object, field_name: str, *, prefix: str) -> str:
    if not isinstance(value, str) or value != value.strip() or not _LOCALE_KEY_RE.fullmatch(value):
        raise UserProblemContractError(
            f"{field_name} must be a namespaced lowercase locale key"
        )
    if not value.startswith(prefix):
        raise UserProblemContractError(f"{field_name} must use the {prefix} namespace")
    return value


def _require_unique(values: tuple[str, ...], field_name: str) -> None:
    if len(values) != len(set(values)):
        raise UserProblemContractError(f"{field_name} must not contain duplicates")


@dataclass(frozen=True, slots=True)
class UserProblemSpec:
    """Reviewed metadata for one coarse, reversible user problem."""

    problem_code: UserProblemCode | str
    severity: UserProblemSeverity | str
    phase: UserProblemPhase | str
    title_key: str
    body_key: str
    action_keys: tuple[str, ...]
    action_target_ids: tuple[str, ...]
    help_topic_id: str

    def __post_init__(self) -> None:
        code = _coerce_enum(UserProblemCode, self.problem_code, "problem code")
        severity = _coerce_enum(UserProblemSeverity, self.severity, "problem severity")
        phase = _coerce_enum(UserProblemPhase, self.phase, "problem phase")
        object.__setattr__(self, "problem_code", code.value)
        object.__setattr__(self, "severity", severity)
        object.__setattr__(self, "phase", phase)

        expected_prefix = f"problems.{code.value}."
        title_key = _require_locale_key(self.title_key, "title_key", prefix=expected_prefix)
        body_key = _require_locale_key(self.body_key, "body_key", prefix=expected_prefix)
        if not title_key.endswith(".title") or not body_key.endswith(".body"):
            raise UserProblemContractError(
                "title_key and body_key must end in .title and .body"
            )
        object.__setattr__(self, "title_key", title_key)
        object.__setattr__(self, "body_key", body_key)

        if isinstance(self.action_keys, (str, bytes, bytearray)):
            raise UserProblemContractError("action_keys must be a tuple of locale keys")
        if isinstance(self.action_target_ids, (str, bytes, bytearray)):
            raise UserProblemContractError(
                "action_target_ids must be a tuple of stable IDs"
            )
        action_keys = tuple(
            _require_locale_key(key, f"action_keys[{index}]", prefix="problems.actions.")
            for index, key in enumerate(self.action_keys)
        )
        action_targets = tuple(
            _require_stable_id(
                target,
                f"action_target_ids[{index}]",
                prefix="target.",
            )
            for index, target in enumerate(self.action_target_ids)
        )
        if not action_keys or len(action_keys) != len(action_targets):
            raise UserProblemContractError(
                "action_keys and action_target_ids must be non-empty paired tuples"
            )
        _require_unique(action_keys, "action_keys")
        _require_unique(action_targets, "action_target_ids")
        object.__setattr__(self, "action_keys", action_keys)
        object.__setattr__(self, "action_target_ids", action_targets)
        object.__setattr__(
            self,
            "help_topic_id",
            _require_stable_id(self.help_topic_id, "help_topic_id", prefix="help."),
        )

    def to_payload(self) -> dict[str, object]:
        """Return reviewed metadata only; no runtime exception is accepted."""

        return {
            "problem_code": self.problem_code,
            "severity": self.severity.value,
            "phase": self.phase.value,
            "title_key": self.title_key,
            "body_key": self.body_key,
            "action_keys": list(self.action_keys),
            "action_target_ids": list(self.action_target_ids),
            "help_topic_id": self.help_topic_id,
        }


def _spec(
    code: UserProblemCode,
    *,
    severity: UserProblemSeverity,
    phase: UserProblemPhase,
    actions: tuple[tuple[str, str], ...],
    help_topic_id: str,
) -> UserProblemSpec:
    """Build one registry entry from its stable code."""

    locale_base = f"problems.{code.value}"
    return UserProblemSpec(
        problem_code=code,
        severity=severity,
        phase=phase,
        title_key=f"{locale_base}.title",
        body_key=f"{locale_base}.body",
        action_keys=tuple(action for action, _ in actions),
        action_target_ids=tuple(target for _, target in actions),
        help_topic_id=help_topic_id,
    )


_USER_PROBLEM_SPEC_ITEMS = (
    _spec(
        UserProblemCode.INPUT_PARSE_FAILED,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.PARSE,
        actions=(("problems.actions.review_input_format", "target.input.format_delimiter"),),
        help_topic_id="help.data.long_format",
    ),
    _spec(
        UserProblemCode.INPUT_MAPPING_INVALID,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.PARSE,
        actions=(("problems.actions.review_column_mapping", "target.input.column_mapping"),),
        help_topic_id="help.data.mapping",
    ),
    _spec(
        UserProblemCode.INPUT_NO_USABLE_ROWS,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.PARSE,
        actions=(("problems.actions.review_data_audit", "target.input.data_audit"),),
        help_topic_id="help.data.missingness",
    ),
    _spec(
        UserProblemCode.ESTIMATION_IDENTIFICATION_FAILED,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.ESTIMATION,
        actions=(("problems.actions.review_design_constraints", "target.settings.design_constraints"),),
        help_topic_id="help.design.coverage",
    ),
    _spec(
        UserProblemCode.ESTIMATION_NONCONVERGENCE,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.ESTIMATION,
        actions=(
            ("problems.actions.review_estimation_settings", "target.settings.estimation"),
            ("problems.actions.review_data_audit", "target.input.data_audit"),
        ),
        help_topic_id="help.run.nonconvergence",
    ),
    _spec(
        UserProblemCode.ESTIMATION_RATING_SCALE_INVALID,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.ESTIMATION,
        actions=(("problems.actions.review_score_categories", "target.input.score_categories"),),
        help_topic_id="help.results.categories",
    ),
    _spec(
        UserProblemCode.ESTIMATION_ANCHOR_INVALID,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.ESTIMATION,
        actions=(("problems.actions.review_anchor_settings", "target.settings.anchors"),),
        help_topic_id="help.run.estimation_failed",
    ),
    _spec(
        UserProblemCode.RESOURCE_HOSTED_LIMIT,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.RESOURCE,
        actions=(("problems.actions.reduce_analysis_scope", "target.settings.analysis_scope"),),
        help_topic_id="help.run.large_analysis",
    ),
    _spec(
        UserProblemCode.DIAGNOSTIC_RESIDUAL_PCA_UNAVAILABLE,
        severity=UserProblemSeverity.INFORMATION,
        phase=UserProblemPhase.RESIDUAL_PCA,
        actions=(("problems.actions.review_diagnostic_request", "target.settings.diagnostics"),),
        help_topic_id="help.results.residual_structure",
    ),
    _spec(
        UserProblemCode.DIAGNOSTIC_RESIDUAL_PCA_FAILED,
        severity=UserProblemSeverity.CAUTION,
        phase=UserProblemPhase.RESIDUAL_PCA,
        actions=(("problems.actions.review_observation_coverage", "target.input.data_audit"),),
        help_topic_id="help.results.residual_structure",
    ),
    _spec(
        UserProblemCode.PROBLEM_UNEXPECTED,
        severity=UserProblemSeverity.BLOCKED,
        phase=UserProblemPhase.APPLICATION,
        actions=(("problems.actions.open_problem_help", "target.help.problem_resolution"),),
        help_topic_id="help.problem.unexpected",
    ),
)


USER_PROBLEM_SPECS: Mapping[str, UserProblemSpec] = MappingProxyType(
    {spec.problem_code: spec for spec in _USER_PROBLEM_SPEC_ITEMS}
)
USER_PROBLEM_CODES = tuple(code.value for code in UserProblemCode)


def validate_user_problem_registry(
    registry: Mapping[str, UserProblemSpec] = USER_PROBLEM_SPECS,
) -> None:
    """Validate exact code coverage and key-to-record identity."""

    if not isinstance(registry, Mapping):
        raise UserProblemContractError("user-problem registry must be a mapping")
    actual = tuple(registry)
    if set(actual) != set(USER_PROBLEM_CODES) or len(actual) != len(USER_PROBLEM_CODES):
        raise UserProblemContractError(
            "user-problem registry must contain every supported problem code exactly once"
        )
    for code, spec in registry.items():
        if not isinstance(spec, UserProblemSpec) or code != spec.problem_code:
            raise UserProblemContractError(
                "user-problem registry keys must match their UserProblemSpec code"
            )


validate_user_problem_registry()


def get_user_problem_spec(problem_code: UserProblemCode | str) -> UserProblemSpec:
    """Resolve one registered specification without echoing rejected input."""

    code = _coerce_enum(UserProblemCode, problem_code, "problem code")
    return USER_PROBLEM_SPECS[code.value]


def generate_support_reference() -> str:
    """Return a random 96-bit support reference.

    There is deliberately no exception, input-data, path, filename,
    identifier, digest, token, or prebuilt-reference argument from which a
    caller could derive or inject the reference.
    """

    token = secrets.token_bytes(_SUPPORT_TOKEN_BYTES)
    token_hex = token.hex().upper()
    return "MFRM-" + "-".join(
        token_hex[index : index + 8] for index in range(0, len(token_hex), 8)
    )


def validate_support_reference(support_reference: object) -> str:
    """Reject a support value unless it has the random-reference format."""

    if not isinstance(support_reference, str) or not _SUPPORT_REFERENCE_RE.fullmatch(
        support_reference
    ):
        raise UserProblemContractError(
            "support_reference must use the MFRM random-reference format"
        )
    return support_reference


@dataclass(frozen=True, slots=True, init=False)
class UserProblemNotice:
    """Serializable public notice containing no exception-derived field.

    Only the three bounded values below can enter the instance.  All
    presentation metadata is projected from the immutable registry, so
    arbitrary prose cannot be attached to a notice by a renderer or exception
    adapter.
    """

    problem_code: UserProblemCode | str
    occurrence_phase: UserProblemPhase | str
    support_reference: str

    def __init__(self) -> None:
        """Prevent callers from attaching a chosen value to a public notice."""

        raise UserProblemContractError(
            "Use build_user_problem_notice() to create a public notice"
        )

    @property
    def spec(self) -> UserProblemSpec:
        return USER_PROBLEM_SPECS[self.problem_code]

    @property
    def severity(self) -> UserProblemSeverity:
        return self.spec.severity

    @property
    def phase(self) -> UserProblemPhase:
        return UserProblemPhase(self.occurrence_phase)

    @property
    def title_key(self) -> str:
        return self.spec.title_key

    @property
    def body_key(self) -> str:
        return self.spec.body_key

    @property
    def action_keys(self) -> tuple[str, ...]:
        return self.spec.action_keys

    @property
    def action_target_ids(self) -> tuple[str, ...]:
        return self.spec.action_target_ids

    @property
    def help_topic_id(self) -> str:
        return self.spec.help_topic_id

    def to_payload(self) -> dict[str, object]:
        """Serialize only stable registry projections and the random reference."""

        return {
            "schema_version": USER_PROBLEM_SCHEMA_VERSION,
            **self.spec.to_payload(),
            "phase": self.phase.value,
            "support_reference": self.support_reference,
        }

    def to_json(self) -> str:
        """Return deterministic JSON suitable for a UI adapter or safe test."""

        return json.dumps(
            self.to_payload(),
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )


def build_user_problem_notice(
    problem_code: UserProblemCode | str,
    *,
    occurrence_phase: UserProblemPhase | str | None = None,
) -> UserProblemNotice:
    """Create a notice for a known non-exception or classified condition."""

    spec = get_user_problem_spec(problem_code)
    phase = (
        spec.phase
        if occurrence_phase is None
        else _coerce_enum(UserProblemPhase, occurrence_phase, "problem phase")
    )
    if (
        spec.problem_code != UserProblemCode.PROBLEM_UNEXPECTED.value
        and phase is not spec.phase
    ):
        raise UserProblemContractError(
            "occurrence phase must match the registered problem phase"
        )
    return _construct_user_problem_notice(spec, phase)


def _construct_user_problem_notice(
    spec: UserProblemSpec,
    occurrence_phase: UserProblemPhase,
) -> UserProblemNotice:
    """Construct from registry-validated values after bounded classification."""

    notice = object.__new__(UserProblemNotice)
    object.__setattr__(notice, "problem_code", spec.problem_code)
    object.__setattr__(notice, "occurrence_phase", occurrence_phase.value)
    object.__setattr__(
        notice,
        "support_reference",
        validate_support_reference(generate_support_reference()),
    )
    return notice


def _exception_match_text(exc: BaseException) -> str:
    """Read bounded classifier input without returning or storing it."""

    if not isinstance(exc, BaseException):
        raise TypeError("exc must be an exception instance")
    try:
        message = str(exc)
    except Exception:
        message = ""
    # Class names can identify common dependency exceptions, but they are also
    # transient and never become notice or log-context fields.
    return f"{type(exc).__name__} {message}".casefold()[:8192]


def _contains_any(text: str, markers: tuple[str, ...]) -> bool:
    return any(marker in text for marker in markers)


_MAPPING_MARKERS = (
    "column mapping",
    "required column",
    "missing column",
    "column not found",
    "unknown column",
    "person_col",
    "score_col",
    "facet_cols",
)
_NO_ROWS_MARKERS = (
    "empty dataframe",
    "empty data frame",
    "emptydataerror",
    "no rows",
    "zero rows",
    "0 rows",
    "no usable",
    "after filtering",
    "no observations",
)
_IDENTIFICATION_MARKERS = (
    "singular",
    "rank deficient",
    "rank-deficient",
    "not identifiable",
    "identifiability",
    "identification",
    "positive semi-definite",
    "positive semidefinite",
    "information matrix",
)
_NONCONVERGENCE_MARKERS = (
    "did not converge",
    "does not converge",
    "nonconvergence",
    "non-convergence",
    "maximum iterations",
    "max iterations",
    "iteration limit",
    "stopping criterion",
)
_RATING_SCALE_MARKERS = (
    "rating scale",
    "rating-scale",
    "score category",
    "score categories",
    "category structure",
    "non-contiguous score",
    "noncontiguous score",
    "all observations are extreme",
)
_ANCHOR_MARKERS = (
    "anchor constraint",
    "anchor setting",
    "anchor table",
    "anchor row",
    "unknown anchor",
    "invalid anchor",
)
_RESOURCE_MARKERS = (
    "memoryerror",
    "out of memory",
    "memory limit",
    "hosted limit",
    "resource limit",
    "safety envelope",
    "analysis too large",
    "request too large",
    "budget exceeded",
)
_PCA_UNAVAILABLE_MARKERS = (
    "residual pca not requested",
    "residual pca not enabled",
    "too few comparable residual profiles",
    "too few complete residual profiles",
    "insufficient observations for residual pca",
    "insufficient comparable residual profiles",
    "no comparable residual profiles",
    "residual comparison graph is too sparse",
    "residual pca prerequisite is missing",
    "residual pca prerequisite is absent",
    "no residual profiles",
)


def classify_user_problem(
    exc: BaseException,
    *,
    phase: UserProblemPhase | str,
) -> UserProblemSpec:
    """Classify an exception into a stable, privacy-safe problem spec.

    The phase is supplied by the calling adapter rather than inferred from a
    traceback.  Matching text determines only a coarse registered code; it is
    immediately discarded and never copied into the returned object.
    """

    normalized_phase = _coerce_enum(UserProblemPhase, phase, "problem phase")
    text = _exception_match_text(exc)

    if normalized_phase is UserProblemPhase.PARSE:
        if _contains_any(text, _NO_ROWS_MARKERS):
            code = UserProblemCode.INPUT_NO_USABLE_ROWS
        elif isinstance(exc, KeyError) or _contains_any(text, _MAPPING_MARKERS):
            code = UserProblemCode.INPUT_MAPPING_INVALID
        else:
            code = UserProblemCode.INPUT_PARSE_FAILED
    elif normalized_phase is UserProblemPhase.ESTIMATION:
        if _contains_any(text, _RESOURCE_MARKERS):
            code = UserProblemCode.RESOURCE_HOSTED_LIMIT
        elif _contains_any(text, _ANCHOR_MARKERS):
            code = UserProblemCode.ESTIMATION_ANCHOR_INVALID
        elif _contains_any(text, _RATING_SCALE_MARKERS):
            code = UserProblemCode.ESTIMATION_RATING_SCALE_INVALID
        elif _contains_any(text, _NO_ROWS_MARKERS):
            code = UserProblemCode.INPUT_NO_USABLE_ROWS
        elif _contains_any(text, _MAPPING_MARKERS):
            code = UserProblemCode.INPUT_MAPPING_INVALID
        elif _contains_any(text, _NONCONVERGENCE_MARKERS):
            code = UserProblemCode.ESTIMATION_NONCONVERGENCE
        elif _contains_any(text, _IDENTIFICATION_MARKERS):
            code = UserProblemCode.ESTIMATION_IDENTIFICATION_FAILED
        else:
            code = UserProblemCode.PROBLEM_UNEXPECTED
    elif normalized_phase is UserProblemPhase.RESIDUAL_PCA:
        code = (
            UserProblemCode.DIAGNOSTIC_RESIDUAL_PCA_UNAVAILABLE
            if _contains_any(text, _PCA_UNAVAILABLE_MARKERS)
            else UserProblemCode.DIAGNOSTIC_RESIDUAL_PCA_FAILED
        )
    elif normalized_phase is UserProblemPhase.RESOURCE:
        code = (
            UserProblemCode.RESOURCE_HOSTED_LIMIT
            if isinstance(exc, MemoryError) or _contains_any(text, _RESOURCE_MARKERS)
            else UserProblemCode.PROBLEM_UNEXPECTED
        )
    else:
        code = UserProblemCode.PROBLEM_UNEXPECTED

    return USER_PROBLEM_SPECS[code.value]


def classify_parse_exception(exc: BaseException) -> UserProblemSpec:
    """Classify an exception raised while reading or validating a table."""

    return classify_user_problem(exc, phase=UserProblemPhase.PARSE)


def classify_estimation_exception(exc: BaseException) -> UserProblemSpec:
    """Classify an exception raised while fitting the requested model."""

    return classify_user_problem(exc, phase=UserProblemPhase.ESTIMATION)


def classify_residual_pca_exception(exc: BaseException) -> UserProblemSpec:
    """Classify an exception raised by the residual-structure screen."""

    return classify_user_problem(exc, phase=UserProblemPhase.RESIDUAL_PCA)


def classify_resource_exception(exc: BaseException) -> UserProblemSpec:
    """Classify an exception raised by a resource limit or preflight step."""

    return classify_user_problem(exc, phase=UserProblemPhase.RESOURCE)


def notice_from_exception(
    exc: BaseException,
    *,
    phase: UserProblemPhase | str,
) -> UserProblemNotice:
    """Classify an exception and return a notice that retains none of it."""

    occurrence_phase = _coerce_enum(UserProblemPhase, phase, "problem phase")
    spec = classify_user_problem(exc, phase=occurrence_phase)
    return _construct_user_problem_notice(spec, occurrence_phase)


def safe_log_context(notice: UserProblemNotice) -> Mapping[str, str]:
    """Return allowlisted structured metadata for confidential logging.

    A logger may attach the original exception separately under its own access
    controls.  This context is safe to retain or aggregate because it contains
    no exception text, traceback, data value, path, environment setting, or
    analysis identity.
    """

    if not isinstance(notice, UserProblemNotice):
        raise TypeError("notice must be a UserProblemNotice")
    return MappingProxyType(
        {
            "event": "mfrm.user_problem",
            "schema_version": USER_PROBLEM_SCHEMA_VERSION,
            "problem_code": notice.problem_code,
            "severity": notice.severity.value,
            "phase": notice.phase.value,
            "support_reference": notice.support_reference,
        }
    )


__all__ = [
    "ProblemCode",
    "ProblemPhase",
    "ProblemSeverity",
    "USER_PROBLEM_CODES",
    "USER_PROBLEM_SCHEMA_VERSION",
    "USER_PROBLEM_SUPPORT_REFERENCE_LABEL_KEY",
    "USER_PROBLEM_SPECS",
    "UserProblemCode",
    "UserProblemContractError",
    "UserProblemNotice",
    "UserProblemPhase",
    "UserProblemSeverity",
    "UserProblemSpec",
    "build_user_problem_notice",
    "classify_estimation_exception",
    "classify_parse_exception",
    "classify_residual_pca_exception",
    "classify_resource_exception",
    "classify_user_problem",
    "generate_support_reference",
    "get_user_problem_spec",
    "notice_from_exception",
    "safe_log_context",
    "validate_support_reference",
    "validate_user_problem_registry",
]
