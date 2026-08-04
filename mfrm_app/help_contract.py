"""Immutable contracts for the application's contextual Help graph.

This module deliberately contains no rendering, translation, data-frame, or
analysis code.  It gives the Streamlit adapter stable, language-independent
identifiers to navigate between a source surface, a review-tracked Help topic, and
an exact return/focus target.

The concrete built-in topics live outside this module.  Keeping the schemas
and graph validation generic lets tests validate a complete registry without
loading Streamlit or the statistical core.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from enum import Enum
import re
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence, TypeVar


HELP_CONTRACT_SCHEMA_VERSION = "mfrm_help_contract_v1"


class HelpContractValidationError(ValueError):
    """Raised when a Help schema or registry violates its contract."""


class HelpContextPolicy(str, Enum):
    """Whether a Help link may display current-analysis content."""

    STATIC_ONLY = "static_only"
    OPTIONAL_CURRENT = "optional_current"
    REQUIRED_CURRENT = "required_current"


class HelpDataContext(str, Enum):
    """Privacy-safe provenance needed to isolate sample and real-data Help."""

    SAMPLE = "sample"
    REAL = "real"


class HelpAnalysisPhase(str, Enum):
    """Identity phase used only to prevent draft/fitted context confusion."""

    DRAFT = "draft"
    FITTED = "fitted"


_EnumT = TypeVar("_EnumT", bound=Enum)
_STABLE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:[._:-][a-z0-9][a-z0-9_-]*)*$")
_LOCALE_KEY_RE = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z0-9_]+)+$")
_OPAQUE_REFERENCE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,255}$")
_CONTENT_VERSION_RE = re.compile(r"^[1-9][0-9]*\.[0-9]+\.[0-9]+$")
_AUDIENCE_LAYERS = frozenset({"standard", "method", "technical"})
_REVIEW_STATES = frozenset({"draft", "in_review", "reviewed"})


def _coerce_enum(enum_type: type[_EnumT], value: _EnumT | str, name: str) -> _EnumT:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(str(value))
    except ValueError as exc:
        allowed = ", ".join(item.value for item in enum_type)
        raise HelpContractValidationError(
            f"Unknown {name} {value!r}; expected one of: {allowed}"
        ) from exc


def _require_stable_id(value: object, name: str) -> str:
    if not isinstance(value, str):
        raise HelpContractValidationError(f"{name} must be a string")
    if value != value.strip() or not _STABLE_ID_RE.fullmatch(value):
        raise HelpContractValidationError(
            f"{name} must be a lowercase, language-independent stable ID"
        )
    return value


def _require_locale_key(value: object, name: str) -> str:
    if not isinstance(value, str):
        raise HelpContractValidationError(f"{name} must be a string")
    if value != value.strip() or not _LOCALE_KEY_RE.fullmatch(value):
        raise HelpContractValidationError(
            f"{name} must be a namespaced lowercase locale key"
        )
    return value


def _require_opaque_reference(value: object, name: str) -> str:
    """Validate an in-session identity without treating it as a route ID."""

    if not isinstance(value, str):
        raise HelpContractValidationError(f"{name} must be a string")
    if value != value.strip() or not _OPAQUE_REFERENCE_RE.fullmatch(value):
        raise HelpContractValidationError(
            f"{name} must be an opaque in-session identity, not free text or a path"
        )
    return value


def _stable_ids(values: Iterable[object], name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise HelpContractValidationError(f"{name} must be a sequence of stable IDs")
    normalized = tuple(
        _require_stable_id(value, f"{name}[{index}]")
        for index, value in enumerate(values)
    )
    _require_unique(normalized, name)
    return normalized


def _locale_keys(values: Iterable[object], name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise HelpContractValidationError(f"{name} must be a sequence of locale keys")
    normalized = tuple(
        _require_locale_key(value, f"{name}[{index}]")
        for index, value in enumerate(values)
    )
    _require_unique(normalized, name)
    return normalized


def _opaque_references(values: Iterable[object], name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise HelpContractValidationError(f"{name} must be a sequence of references")
    normalized = tuple(
        _require_opaque_reference(value, f"{name}[{index}]")
        for index, value in enumerate(values)
    )
    _require_unique(normalized, name)
    return normalized


def _require_unique(values: Sequence[str], name: str) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    if duplicates:
        raise HelpContractValidationError(
            f"{name} must be unique; duplicates={sorted(duplicates)!r}"
        )


def _freeze_json(value: object, name: str) -> object:
    """Recursively freeze the small JSON-like presentation-state payload."""

    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        frozen: dict[str, object] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise HelpContractValidationError(f"{name} keys must be strings")
            frozen[key] = _freeze_json(item, f"{name}.{key}")
        return MappingProxyType(dict(sorted(frozen.items())))
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return tuple(_freeze_json(item, name) for item in value)
    raise HelpContractValidationError(
        f"{name} contains unsupported value type {type(value).__name__}"
    )


def _freeze_mapping(value: Mapping[str, object], name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise HelpContractValidationError(f"{name} must be a mapping")
    frozen = _freeze_json(value, name)
    if not isinstance(frozen, Mapping):  # defensive; mappings normalize to mappings
        raise HelpContractValidationError(f"{name} must normalize to a mapping")
    return frozen


def _freeze_review_metadata(value: Mapping[str, object]) -> Mapping[str, object]:
    """Validate accountable content review metadata before freezing it."""

    if not isinstance(value, Mapping):
        raise HelpContractValidationError("version_and_review must be a mapping")
    required_fields = {
        "content_version",
        "owner_role",
        "review_state",
        "last_reviewed",
        "reviewer_roles",
    }
    if set(value) != required_fields:
        raise HelpContractValidationError(
            "version_and_review must contain exactly content_version, owner_role, "
            "review_state, last_reviewed, and reviewer_roles"
        )

    content_version = value["content_version"]
    if not isinstance(content_version, str) or not _CONTENT_VERSION_RE.fullmatch(
        content_version
    ):
        raise HelpContractValidationError(
            "content_version must use MAJOR.MINOR.PATCH with a positive major"
        )
    owner_role = _require_stable_id(value["owner_role"], "owner_role")
    review_state = value["review_state"]
    if review_state not in _REVIEW_STATES:
        raise HelpContractValidationError(
            "review_state must be draft, in_review, or reviewed"
        )
    last_reviewed = value["last_reviewed"]
    if not isinstance(last_reviewed, str):
        raise HelpContractValidationError("last_reviewed must be an ISO date")
    try:
        date.fromisoformat(last_reviewed)
    except ValueError as exc:
        raise HelpContractValidationError("last_reviewed must be an ISO date") from exc

    reviewer_roles = _stable_ids(value["reviewer_roles"], "reviewer_roles")
    if review_state in {"in_review", "reviewed"} and not reviewer_roles:
        raise HelpContractValidationError(
            "in-review and reviewed Help must name at least one reviewer role"
        )
    if owner_role in reviewer_roles:
        raise HelpContractValidationError(
            "owner_role must not also be listed as a reviewer"
        )
    return MappingProxyType(
        {
            "content_version": content_version,
            "last_reviewed": last_reviewed,
            "owner_role": owner_role,
            "review_state": review_state,
            "reviewer_roles": reviewer_roles,
        }
    )


@dataclass(frozen=True)
class HelpSection:
    """One stable subsection within a Help topic."""

    section_id: str
    title_key: str
    content_keys: tuple[str, ...]
    audience_layers: tuple[str, ...] = ("standard", "method", "technical")

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "section_id", _require_stable_id(self.section_id, "section_id")
        )
        object.__setattr__(self, "title_key", _require_locale_key(self.title_key, "title_key"))
        content_keys = _locale_keys(self.content_keys, "content_keys")
        if not content_keys:
            raise HelpContractValidationError("content_keys must not be empty")
        object.__setattr__(self, "content_keys", content_keys)
        layers = tuple(str(layer) for layer in self.audience_layers)
        if not layers or any(layer not in _AUDIENCE_LAYERS for layer in layers):
            raise HelpContractValidationError(
                "audience_layers must contain only standard, method, or technical"
            )
        _require_unique(layers, "audience_layers")
        object.__setattr__(self, "audience_layers", layers)


@dataclass(frozen=True)
class HelpTopic:
    """Review-tracked, locale-independent metadata for a Help topic.

    ``lifecycle_states`` is an applicability declaration for a rendering
    adapter, not a permission computed by this graph.  An adapter must call
    :meth:`supports_lifecycle` before presenting restricted topic content.
    """

    help_topic_id: str
    title_key: str
    summary_key: str
    sections: tuple[HelpSection, ...]
    version_and_review: Mapping[str, object]
    concept_ids: tuple[str, ...] = ()
    audience_layers: tuple[str, ...] = ("standard", "method", "technical")
    lifecycle_states: tuple[str, ...] = ("all",)
    applicability: Mapping[str, object] = field(default_factory=dict)
    prerequisite_keys: tuple[str, ...] = ()
    computed_key: str | None = None
    can_show_key: str | None = None
    cannot_show_key: str | None = None
    next_check_key: str | None = None
    next_action_key: str | None = None
    safe_report_guard_key: str | None = None
    safe_report_key: str | None = None
    avoid_report_key: str | None = None
    claim_boundary_ids: tuple[str, ...] = ()
    reference_ids: tuple[str, ...] = ()
    related_target_ids: tuple[str, ...] = ()
    search_alias_keys: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "help_topic_id",
            _require_stable_id(self.help_topic_id, "help_topic_id"),
        )
        object.__setattr__(self, "title_key", _require_locale_key(self.title_key, "title_key"))
        object.__setattr__(
            self, "summary_key", _require_locale_key(self.summary_key, "summary_key")
        )

        if isinstance(self.sections, (str, bytes, bytearray)):
            raise HelpContractValidationError("sections must be a sequence")
        sections = tuple(self.sections)
        if not sections or any(not isinstance(section, HelpSection) for section in sections):
            raise HelpContractValidationError(
                "sections must contain at least one HelpSection"
            )
        _require_unique(tuple(section.section_id for section in sections), "section IDs")
        object.__setattr__(self, "sections", sections)

        object.__setattr__(self, "concept_ids", _stable_ids(self.concept_ids, "concept_ids"))
        layers = tuple(str(layer) for layer in self.audience_layers)
        if not layers or any(layer not in _AUDIENCE_LAYERS for layer in layers):
            raise HelpContractValidationError(
                "audience_layers must contain only standard, method, or technical"
            )
        _require_unique(layers, "audience_layers")
        object.__setattr__(self, "audience_layers", layers)

        lifecycle_states = _stable_ids(self.lifecycle_states, "lifecycle_states")
        if not lifecycle_states:
            raise HelpContractValidationError("lifecycle_states must not be empty")
        if "all" in lifecycle_states and len(lifecycle_states) != 1:
            raise HelpContractValidationError(
                "the all lifecycle state cannot be combined with restricted states"
            )
        object.__setattr__(self, "lifecycle_states", lifecycle_states)
        object.__setattr__(
            self, "applicability", _freeze_mapping(self.applicability, "applicability")
        )
        object.__setattr__(
            self,
            "prerequisite_keys",
            _locale_keys(self.prerequisite_keys, "prerequisite_keys"),
        )

        for field_name in (
            "computed_key",
            "can_show_key",
            "cannot_show_key",
            "next_check_key",
            "next_action_key",
            "safe_report_guard_key",
            "safe_report_key",
            "avoid_report_key",
        ):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self, field_name, _require_locale_key(value, field_name)
                )

        if (self.safe_report_key is None) != (self.safe_report_guard_key is None):
            raise HelpContractValidationError(
                "safe_report_key and safe_report_guard_key must be declared together"
            )

        object.__setattr__(
            self,
            "claim_boundary_ids",
            _stable_ids(self.claim_boundary_ids, "claim_boundary_ids"),
        )
        object.__setattr__(
            self, "reference_ids", _stable_ids(self.reference_ids, "reference_ids")
        )
        object.__setattr__(
            self,
            "related_target_ids",
            _stable_ids(self.related_target_ids, "related_target_ids"),
        )
        object.__setattr__(
            self,
            "search_alias_keys",
            _locale_keys(self.search_alias_keys, "search_alias_keys"),
        )
        object.__setattr__(
            self,
            "version_and_review",
            _freeze_review_metadata(self.version_and_review),
        )

        section_layers = {layer for section in sections for layer in section.audience_layers}
        unknown_section_layers = section_layers.difference(layers)
        if unknown_section_layers:
            raise HelpContractValidationError(
                "section audience layers must be declared by the parent topic; "
                f"unknown={sorted(unknown_section_layers)!r}"
            )

    @property
    def section_ids(self) -> tuple[str, ...]:
        """Return stable section IDs in declared reading order."""

        return tuple(section.section_id for section in self.sections)

    def has_section(self, section_id: str) -> bool:
        return section_id in self.section_ids

    def supports_lifecycle(self, lifecycle_state: str) -> bool:
        """Return whether an adapter may render this topic in one lifecycle."""

        normalized = _require_stable_id(lifecycle_state, "lifecycle_state")
        return "all" in self.lifecycle_states or normalized in self.lifecycle_states


@dataclass(frozen=True)
class HelpTarget:
    """A stable application location and the focus markers it owns."""

    target_id: str
    surface_id: str
    focus_ids: tuple[str, ...]
    presentation_state: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "target_id", _require_stable_id(self.target_id, "target_id"))
        object.__setattr__(
            self, "surface_id", _require_stable_id(self.surface_id, "surface_id")
        )
        focus_ids = _stable_ids(self.focus_ids, "focus_ids")
        if not focus_ids:
            raise HelpContractValidationError("focus_ids must not be empty")
        object.__setattr__(self, "focus_ids", focus_ids)
        object.__setattr__(
            self,
            "presentation_state",
            _freeze_mapping(self.presentation_state, "presentation_state"),
        )


@dataclass(frozen=True)
class HelpLink:
    """A stable edge from one application target to one Help subsection."""

    help_link_id: str
    source_target_id: str
    help_topic_id: str
    section_id: str | None
    return_target_id: str
    return_focus_id: str
    context_policy: HelpContextPolicy | str = HelpContextPolicy.STATIC_ONLY

    def __post_init__(self) -> None:
        for field_name in (
            "help_link_id",
            "source_target_id",
            "help_topic_id",
            "return_target_id",
            "return_focus_id",
        ):
            object.__setattr__(
                self,
                field_name,
                _require_stable_id(getattr(self, field_name), field_name),
            )
        if self.section_id is not None:
            object.__setattr__(
                self,
                "section_id",
                _require_stable_id(self.section_id, "section_id"),
            )
        object.__setattr__(
            self,
            "context_policy",
            _coerce_enum(HelpContextPolicy, self.context_policy, "context policy"),
        )


@dataclass(frozen=True)
class HelpContextBinding:
    """Private, in-session references used to verify dynamic Help context.

    The binding contains identity references only.  It intentionally has no
    estimates, diagnostic values, readiness state, decision state, free text,
    or tutorial/learning completion field.
    """

    data_context: HelpDataContext | str
    analysis_id: str | None = None
    study_context_id: str | None = None
    evidence_ids: tuple[str, ...] = ()
    evidence_issue_id: str | None = None
    claim_boundary_ids: tuple[str, ...] = ()
    analysis_phase: HelpAnalysisPhase | str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "data_context",
            _coerce_enum(HelpDataContext, self.data_context, "data context"),
        )
        for field_name in ("analysis_id", "study_context_id", "evidence_issue_id"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self,
                    field_name,
                    _require_opaque_reference(value, field_name),
                )
        object.__setattr__(
            self,
            "evidence_ids",
            _opaque_references(self.evidence_ids, "evidence_ids"),
        )
        object.__setattr__(
            self,
            "claim_boundary_ids",
            _stable_ids(self.claim_boundary_ids, "claim_boundary_ids"),
        )
        if self.analysis_phase is not None:
            object.__setattr__(
                self,
                "analysis_phase",
                _coerce_enum(HelpAnalysisPhase, self.analysis_phase, "analysis phase"),
            )

    @property
    def has_current_analysis(self) -> bool:
        """Whether this binding can identify a particular analysis."""

        return self.analysis_id is not None


@dataclass(frozen=True)
class HelpRegistry:
    """A validated, immutable Help graph used by pure navigation code."""

    topics: tuple[HelpTopic, ...]
    links: tuple[HelpLink, ...]
    targets: tuple[HelpTarget, ...]
    entry_topic_ids: tuple[str, ...] = ()
    system_topic_ids: tuple[str, ...] = ()
    require_zero_orphans: bool = True

    def __post_init__(self) -> None:
        if any(not isinstance(topic, HelpTopic) for topic in self.topics):
            raise HelpContractValidationError("topics must contain only HelpTopic values")
        if any(not isinstance(link, HelpLink) for link in self.links):
            raise HelpContractValidationError("links must contain only HelpLink values")
        if any(not isinstance(target, HelpTarget) for target in self.targets):
            raise HelpContractValidationError("targets must contain only HelpTarget values")
        object.__setattr__(self, "topics", tuple(self.topics))
        object.__setattr__(self, "links", tuple(self.links))
        object.__setattr__(self, "targets", tuple(self.targets))
        object.__setattr__(
            self,
            "entry_topic_ids",
            _stable_ids(self.entry_topic_ids, "entry_topic_ids"),
        )
        object.__setattr__(
            self,
            "system_topic_ids",
            _stable_ids(self.system_topic_ids, "system_topic_ids"),
        )
        if not isinstance(self.require_zero_orphans, bool):
            raise HelpContractValidationError("require_zero_orphans must be a boolean")
        validate_help_registry(
            self.topics,
            self.links,
            self.targets,
            entry_topic_ids=self.entry_topic_ids,
            system_topic_ids=self.system_topic_ids,
            require_zero_orphans=self.require_zero_orphans,
        )

    @property
    def topic_map(self) -> Mapping[str, HelpTopic]:
        return MappingProxyType({topic.help_topic_id: topic for topic in self.topics})

    @property
    def link_map(self) -> Mapping[str, HelpLink]:
        return MappingProxyType({link.help_link_id: link for link in self.links})

    @property
    def target_map(self) -> Mapping[str, HelpTarget]:
        return MappingProxyType({target.target_id: target for target in self.targets})

    def topic(self, help_topic_id: str) -> HelpTopic | None:
        return self.topic_map.get(help_topic_id)

    def link(self, help_link_id: str) -> HelpLink | None:
        return self.link_map.get(help_link_id)

    def target(self, target_id: str) -> HelpTarget | None:
        return self.target_map.get(target_id)


def _duplicates(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def validate_help_registry(
    topics: Sequence[HelpTopic],
    links: Sequence[HelpLink],
    targets: Sequence[HelpTarget],
    *,
    entry_topic_ids: Sequence[str] = (),
    system_topic_ids: Sequence[str] = (),
    known_concept_ids: Iterable[str] | None = None,
    known_claim_boundary_ids: Iterable[str] | None = None,
    known_reference_ids: Iterable[str] | None = None,
    require_zero_orphans: bool = True,
) -> None:
    """Validate references, sections, return focus, and graph reachability.

    ``entry_topic_ids`` are topics directly reachable from Help home.
    ``system_topic_ids`` are reserved destinations such as a fail-closed
    fallback and must not also be advertised as Help-home entries.  Every
    other topic must have an incoming :class:`HelpLink` when zero-orphan
    validation is enabled.  Targets are reachable when a link or a topic
    references them.
    """

    if any(not isinstance(topic, HelpTopic) for topic in topics):
        raise HelpContractValidationError("topics must contain only HelpTopic values")
    if any(not isinstance(link, HelpLink) for link in links):
        raise HelpContractValidationError("links must contain only HelpLink values")
    if any(not isinstance(target, HelpTarget) for target in targets):
        raise HelpContractValidationError("targets must contain only HelpTarget values")
    if not isinstance(require_zero_orphans, bool):
        raise HelpContractValidationError("require_zero_orphans must be a boolean")

    topic_ids = tuple(topic.help_topic_id for topic in topics)
    link_ids = tuple(link.help_link_id for link in links)
    target_ids = tuple(target.target_id for target in targets)
    for name, identifiers in (
        ("help_topic_id", topic_ids),
        ("help_link_id", link_ids),
        ("target_id", target_ids),
    ):
        duplicates = _duplicates(identifiers)
        if duplicates:
            raise HelpContractValidationError(
                f"Duplicate {name} values: {duplicates!r}"
            )

    focus_owners: dict[str, str] = {}
    ambiguous_focus_ids: set[str] = set()
    for target in targets:
        for focus_id in target.focus_ids:
            previous_owner = focus_owners.setdefault(focus_id, target.target_id)
            if previous_owner != target.target_id:
                ambiguous_focus_ids.add(focus_id)
    if ambiguous_focus_ids:
        raise HelpContractValidationError(
            "focus_ids must have one global target owner; duplicates="
            f"{sorted(ambiguous_focus_ids)!r}"
        )

    topic_map = {topic.help_topic_id: topic for topic in topics}
    target_map = {target.target_id: target for target in targets}
    entries = _stable_ids(entry_topic_ids, "entry_topic_ids")
    system_topics = _stable_ids(system_topic_ids, "system_topic_ids")
    overlap = sorted(set(entries).intersection(system_topics))
    if overlap:
        raise HelpContractValidationError(
            f"Help topics cannot be both entry and system topics: {overlap!r}"
        )
    unknown_entries = sorted(set((*entries, *system_topics)).difference(topic_map))
    if unknown_entries:
        raise HelpContractValidationError(
            f"Unknown Help entry/system topic IDs: {unknown_entries!r}"
        )

    for topic in topics:
        unknown_targets = sorted(set(topic.related_target_ids).difference(target_map))
        if unknown_targets:
            raise HelpContractValidationError(
                f"Help topic {topic.help_topic_id!r} references unknown targets: "
                f"{unknown_targets!r}"
            )

    for link in links:
        topic = topic_map.get(link.help_topic_id)
        if topic is None:
            raise HelpContractValidationError(
                f"Help link {link.help_link_id!r} references unknown topic "
                f"{link.help_topic_id!r}"
            )
        if link.section_id is not None and not topic.has_section(link.section_id):
            raise HelpContractValidationError(
                f"Help link {link.help_link_id!r} references unknown section "
                f"{link.section_id!r} in topic {link.help_topic_id!r}"
            )
        if link.source_target_id not in target_map:
            raise HelpContractValidationError(
                f"Help link {link.help_link_id!r} references unknown source target "
                f"{link.source_target_id!r}"
            )
        return_target = target_map.get(link.return_target_id)
        if return_target is None:
            raise HelpContractValidationError(
                f"Help link {link.help_link_id!r} references unknown return target "
                f"{link.return_target_id!r}"
            )
        if link.return_focus_id not in return_target.focus_ids:
            raise HelpContractValidationError(
                f"Help link {link.help_link_id!r} return focus "
                f"{link.return_focus_id!r} is not owned by target "
                f"{link.return_target_id!r}"
            )

    _validate_known_ids(
        topics,
        attribute="concept_ids",
        known_ids=known_concept_ids,
        label="concept",
    )
    _validate_known_ids(
        topics,
        attribute="claim_boundary_ids",
        known_ids=known_claim_boundary_ids,
        label="ClaimBoundary",
    )
    _validate_known_ids(
        topics,
        attribute="reference_ids",
        known_ids=known_reference_ids,
        label="reference",
    )

    if require_zero_orphans:
        reachable_topics = set((*entries, *system_topics))
        reachable_topics.update(link.help_topic_id for link in links)
        orphan_topics = sorted(set(topic_ids).difference(reachable_topics))

        reachable_targets: set[str] = set()
        for link in links:
            reachable_targets.add(link.source_target_id)
            reachable_targets.add(link.return_target_id)
        for topic in topics:
            reachable_targets.update(topic.related_target_ids)
        orphan_targets = sorted(set(target_ids).difference(reachable_targets))
        if orphan_topics or orphan_targets:
            raise HelpContractValidationError(
                "Help graph contains orphan nodes; "
                f"topics={orphan_topics!r}, targets={orphan_targets!r}"
            )


def _validate_known_ids(
    topics: Sequence[HelpTopic],
    *,
    attribute: str,
    known_ids: Iterable[str] | None,
    label: str,
) -> None:
    if known_ids is None:
        return
    known = set(_stable_ids(known_ids, f"known_{attribute}"))
    for topic in topics:
        values = set(getattr(topic, attribute))
        unknown = sorted(values.difference(known))
        if unknown:
            raise HelpContractValidationError(
                f"Help topic {topic.help_topic_id!r} references unknown {label} IDs: "
                f"{unknown!r}"
            )


def build_help_registry(
    topics: Iterable[HelpTopic],
    links: Iterable[HelpLink],
    targets: Iterable[HelpTarget],
    *,
    entry_topic_ids: Iterable[str] = (),
    system_topic_ids: Iterable[str] = (),
    require_zero_orphans: bool = True,
) -> HelpRegistry:
    """Build and validate an immutable Help registry from arbitrary iterables."""

    return HelpRegistry(
        topics=tuple(topics),
        links=tuple(links),
        targets=tuple(targets),
        entry_topic_ids=tuple(entry_topic_ids),
        system_topic_ids=tuple(system_topic_ids),
        require_zero_orphans=require_zero_orphans,
    )


__all__ = [
    "HELP_CONTRACT_SCHEMA_VERSION",
    "HelpAnalysisPhase",
    "HelpContextBinding",
    "HelpContextPolicy",
    "HelpContractValidationError",
    "HelpDataContext",
    "HelpLink",
    "HelpRegistry",
    "HelpSection",
    "HelpTarget",
    "HelpTopic",
    "build_help_registry",
    "validate_help_registry",
]
