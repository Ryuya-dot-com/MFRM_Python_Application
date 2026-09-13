"""Immutable audience-layer terminology contracts for the Python app.

The registry separates a concept's stable identity from the words exposed on
three audiences:

``standard``
    Task, meaning, limitation, and next-action language used by the default UI.
``method``
    Reviewed SLA and measurement terminology with an interpretation boundary.
``technical``
    Exact support, reproducibility, and machine-contract labels.

Only locale keys are returned for localized layers.  Exact implementation
labels are retained, but a projection never carries labels or aliases from a
different layer.  This makes accidental exposure of ``AnalysisID`` and similar
terms less likely while preserving them for technical exports and support.

This module deliberately has no Streamlit, pandas, estimator, or locale-file
dependency.  A later rendering adapter is responsible for resolving locale
keys to versioned English and Japanese copy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from enum import Enum
import re
from types import MappingProxyType
from typing import Iterable, Iterator, Mapping


TERMINOLOGY_REGISTRY_VERSION = "mfrm_terminology_v1"


class TerminologyValidationError(ValueError):
    """Raised when a term or terminology registry violates its contract."""


class UnknownConceptError(KeyError):
    """Raised when a stable concept ID is not registered."""


class LayerNotAllowedError(ValueError):
    """Raised when a concept has no registered projection for an audience layer."""


class AudienceLayer(str, Enum):
    """Reviewed terminology layers; these are not analysis modes."""

    STANDARD = "standard"
    METHOD = "method"
    TECHNICAL = "technical"


class ReviewState(str, Enum):
    """Review state of terminology content and its interpretation boundary."""

    DRAFT = "draft"
    IN_REVIEW = "in_review"
    REVIEWED = "reviewed"


_STABLE_ID_RE = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+$")
_LOCALE_KEY_RE = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+$")
_VERSION_RE = re.compile(r"^[1-9][0-9]*\.[0-9]+\.[0-9]+$")
_REGISTRY_VERSION_RE = re.compile(r"^[a-z][a-z0-9_]*_v[1-9][0-9]*$")
_ROLE_RE = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$")


def _require_text(value: object, name: str) -> str:
    if not isinstance(value, str):
        raise TerminologyValidationError(f"{name} must be a string")
    text = value.strip()
    if not text:
        raise TerminologyValidationError(f"{name} must be a non-empty string")
    return text


def _require_stable_id(value: object, name: str) -> str:
    text = _require_text(value, name)
    if not _STABLE_ID_RE.fullmatch(text):
        raise TerminologyValidationError(
            f"{name} must be a lowercase namespaced ID such as "
            "'evidence.availability'"
        )
    return text


def _require_locale_key(value: object, name: str) -> str:
    text = _require_text(value, name)
    if not _LOCALE_KEY_RE.fullmatch(text) or not text.startswith("terms."):
        raise TerminologyValidationError(
            f"{name} must be a lowercase locale key below 'terms.'"
        )
    return text


def _optional_locale_key(value: object, name: str) -> str | None:
    if value is None:
        return None
    return _require_locale_key(value, name)


def _require_tuple(value: object, name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes, bytearray)):
        raise TerminologyValidationError(f"{name} must be a sequence, not text")
    try:
        return tuple(value)  # type: ignore[arg-type]
    except TypeError as exc:
        raise TerminologyValidationError(f"{name} must be a sequence") from exc


def _locale_key_tuple(value: object, name: str, *, nonempty: bool = True) -> tuple[str, ...]:
    raw_items = _require_tuple(value, name)
    items = tuple(
        _require_locale_key(item, f"{name}[{index}]")
        for index, item in enumerate(raw_items)
    )
    if nonempty and not items:
        raise TerminologyValidationError(f"{name} must contain at least one key")
    if len(items) != len(set(items)):
        raise TerminologyValidationError(f"{name} must not contain duplicate keys")
    return items


def _stable_id_tuple(value: object, name: str, *, nonempty: bool = True) -> tuple[str, ...]:
    raw_items = _require_tuple(value, name)
    items = tuple(
        _require_stable_id(item, f"{name}[{index}]")
        for index, item in enumerate(raw_items)
    )
    if nonempty and not items:
        raise TerminologyValidationError(f"{name} must contain at least one ID")
    if len(items) != len(set(items)):
        raise TerminologyValidationError(f"{name} must not contain duplicate IDs")
    return items


def _coerce_layer(value: AudienceLayer | str, name: str = "audience layer") -> AudienceLayer:
    if isinstance(value, AudienceLayer):
        return value
    try:
        return AudienceLayer(str(value))
    except ValueError as exc:
        allowed = ", ".join(layer.value for layer in AudienceLayer)
        raise TerminologyValidationError(
            f"Unknown {name} {value!r}; expected one of: {allowed}"
        ) from exc


@dataclass(frozen=True)
class VersionReviewMetadata:
    """Version and accountable review metadata for a term definition."""

    content_version: str
    owner_role: str
    last_reviewed: date
    review_state: ReviewState
    reviewer_roles: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        version = _require_text(self.content_version, "content_version")
        if not _VERSION_RE.fullmatch(version):
            raise TerminologyValidationError(
                "content_version must use MAJOR.MINOR.PATCH with a positive major version"
            )

        owner = _require_text(self.owner_role, "owner_role")
        if not _ROLE_RE.fullmatch(owner):
            raise TerminologyValidationError(
                "owner_role must be a stable lowercase role ID"
            )
        if type(self.last_reviewed) is not date:
            raise TerminologyValidationError("last_reviewed must be a date")

        try:
            state = (
                self.review_state
                if isinstance(self.review_state, ReviewState)
                else ReviewState(str(self.review_state))
            )
        except ValueError as exc:
            raise TerminologyValidationError(
                f"Unknown review_state {self.review_state!r}"
            ) from exc

        raw_reviewers = _require_tuple(self.reviewer_roles, "reviewer_roles")
        reviewers = tuple(
            _require_text(role, f"reviewer_roles[{index}]")
            for index, role in enumerate(raw_reviewers)
        )
        if any(not _ROLE_RE.fullmatch(role) for role in reviewers):
            raise TerminologyValidationError(
                "reviewer_roles must contain stable lowercase role IDs"
            )
        if len(reviewers) != len(set(reviewers)):
            raise TerminologyValidationError("reviewer_roles must not contain duplicates")
        if owner in reviewers:
            raise TerminologyValidationError(
                "owner_role must not also be listed as a reviewer"
            )
        if state is ReviewState.REVIEWED and not reviewers:
            raise TerminologyValidationError(
                "reviewed terminology must name at least one reviewer role"
            )

        object.__setattr__(self, "content_version", version)
        object.__setattr__(self, "owner_role", owner)
        object.__setattr__(self, "review_state", state)
        object.__setattr__(self, "reviewer_roles", reviewers)


@dataclass(frozen=True)
class TermProjection:
    """The minimum safe term payload for one audience layer.

    Localized layers expose ``label_key``.  The technical layer exposes
    ``literal_label``.  A projection intentionally carries neither aliases nor
    labels belonging to another audience layer.
    """

    concept_id: str
    layer: AudienceLayer
    label_key: str | None
    literal_label: str | None
    definition_key: str
    prohibited_claim_keys: tuple[str, ...]
    related_help_topic_ids: tuple[str, ...]
    content_version: str
    review_state: ReviewState

    def __post_init__(self) -> None:
        concept_id = _require_stable_id(self.concept_id, "concept_id")
        layer = _coerce_layer(self.layer)
        label_key = _optional_locale_key(self.label_key, "label_key")
        literal_label = (
            None
            if self.literal_label is None
            else _require_text(self.literal_label, "literal_label")
        )
        if (label_key is None) == (literal_label is None):
            raise TerminologyValidationError(
                "A projection must expose exactly one of label_key or literal_label"
            )
        if layer is AudienceLayer.TECHNICAL and literal_label is None:
            raise TerminologyValidationError(
                "A technical projection must expose only a literal label"
            )
        if layer is not AudienceLayer.TECHNICAL and label_key is None:
            raise TerminologyValidationError(
                "A localized projection must expose only a locale label key"
            )
        definition_key = _require_locale_key(self.definition_key, "definition_key")
        prohibited_claim_keys = _locale_key_tuple(
            self.prohibited_claim_keys,
            "prohibited_claim_keys",
        )
        related_help_topic_ids = _stable_id_tuple(
            self.related_help_topic_ids,
            "related_help_topic_ids",
        )
        content_version = _require_text(self.content_version, "content_version")
        if not _VERSION_RE.fullmatch(content_version):
            raise TerminologyValidationError(
                "content_version must use MAJOR.MINOR.PATCH with a positive major version"
            )
        try:
            review_state = (
                self.review_state
                if isinstance(self.review_state, ReviewState)
                else ReviewState(str(self.review_state))
            )
        except ValueError as exc:
            raise TerminologyValidationError(
                f"Unknown review_state {self.review_state!r}"
            ) from exc

        object.__setattr__(self, "concept_id", concept_id)
        object.__setattr__(self, "layer", layer)
        object.__setattr__(self, "label_key", label_key)
        object.__setattr__(self, "literal_label", literal_label)
        object.__setattr__(self, "definition_key", definition_key)
        object.__setattr__(self, "prohibited_claim_keys", prohibited_claim_keys)
        object.__setattr__(self, "related_help_topic_ids", related_help_topic_ids)
        object.__setattr__(self, "content_version", content_version)
        object.__setattr__(self, "review_state", review_state)

    @property
    def label_reference(self) -> str:
        """Return the locale key or exact technical label for this projection."""

        return self.label_key if self.label_key is not None else str(self.literal_label)


@dataclass(frozen=True)
class TermSpec:
    """One stable concept and its review-tracked audience terminology."""

    concept_id: str
    standard_label_key: str | None
    method_label_key: str | None
    technical_label: str | None
    definition_key: str
    aliases_keys: tuple[str, ...]
    allowed_layers: frozenset[AudienceLayer]
    prohibited_claim_keys: tuple[str, ...]
    related_help_topic_ids: tuple[str, ...]
    version_and_review: VersionReviewMetadata

    def __post_init__(self) -> None:
        concept_id = _require_stable_id(self.concept_id, "concept_id")
        standard_key = _optional_locale_key(
            self.standard_label_key, "standard_label_key"
        )
        method_key = _optional_locale_key(self.method_label_key, "method_label_key")
        technical_label = (
            None
            if self.technical_label is None
            else _require_text(self.technical_label, "technical_label")
        )
        definition_key = _require_locale_key(self.definition_key, "definition_key")
        aliases = _locale_key_tuple(self.aliases_keys, "aliases_keys")
        prohibited = _locale_key_tuple(
            self.prohibited_claim_keys, "prohibited_claim_keys"
        )
        help_topics = _stable_id_tuple(
            self.related_help_topic_ids, "related_help_topic_ids"
        )

        raw_layers = _require_tuple(self.allowed_layers, "allowed_layers")
        layers = frozenset(
            _coerce_layer(layer, f"allowed_layers[{index}]")
            for index, layer in enumerate(raw_layers)
        )
        if not layers:
            raise TerminologyValidationError(
                "allowed_layers must contain at least one audience layer"
            )

        labels: Mapping[AudienceLayer, str | None] = {
            AudienceLayer.STANDARD: standard_key,
            AudienceLayer.METHOD: method_key,
            AudienceLayer.TECHNICAL: technical_label,
        }
        for layer, label in labels.items():
            if layer in layers and label is None:
                raise TerminologyValidationError(
                    f"{concept_id!r} allows {layer.value!r} but has no label for it"
                )
            if layer not in layers and label is not None:
                raise TerminologyValidationError(
                    f"{concept_id!r} defines a {layer.value!r} label but does not "
                    "allow that layer"
                )

        localized_keys = tuple(
            key
            for key in (standard_key, method_key, definition_key, *aliases, *prohibited)
            if key is not None
        )
        if len(localized_keys) != len(set(localized_keys)):
            raise TerminologyValidationError(
                f"{concept_id!r} reuses a locale key across different term fields"
            )
        if not isinstance(self.version_and_review, VersionReviewMetadata):
            raise TerminologyValidationError(
                "version_and_review must be VersionReviewMetadata"
            )

        object.__setattr__(self, "concept_id", concept_id)
        object.__setattr__(self, "standard_label_key", standard_key)
        object.__setattr__(self, "method_label_key", method_key)
        object.__setattr__(self, "technical_label", technical_label)
        object.__setattr__(self, "definition_key", definition_key)
        object.__setattr__(self, "aliases_keys", aliases)
        object.__setattr__(self, "allowed_layers", layers)
        object.__setattr__(self, "prohibited_claim_keys", prohibited)
        object.__setattr__(self, "related_help_topic_ids", help_topics)

    @property
    def aliases(self) -> tuple[str, ...]:
        """Stable locale keys for review-tracked search aliases."""

        return self.aliases_keys

    @property
    def locale_keys(self) -> frozenset[str]:
        """All locale keys required to render and search this term."""

        return frozenset(
            key
            for key in (
                self.standard_label_key,
                self.method_label_key,
                self.definition_key,
                *self.aliases_keys,
                *self.prohibited_claim_keys,
            )
            if key is not None
        )

    def allows(self, layer: AudienceLayer | str) -> bool:
        """Return whether this concept has registered content for ``layer``."""

        return _coerce_layer(layer) in self.allowed_layers

    def project(self, layer: AudienceLayer | str) -> TermProjection:
        """Return a data-minimized projection for one registered layer."""

        normalized_layer = _coerce_layer(layer)
        if normalized_layer not in self.allowed_layers:
            raise LayerNotAllowedError(
                f"Concept {self.concept_id!r} is not available in the "
                f"{normalized_layer.value!r} layer"
            )

        label_key: str | None = None
        literal_label: str | None = None
        if normalized_layer is AudienceLayer.STANDARD:
            label_key = self.standard_label_key
        elif normalized_layer is AudienceLayer.METHOD:
            label_key = self.method_label_key
        else:
            literal_label = self.technical_label

        return TermProjection(
            concept_id=self.concept_id,
            layer=normalized_layer,
            label_key=label_key,
            literal_label=literal_label,
            definition_key=self.definition_key,
            prohibited_claim_keys=self.prohibited_claim_keys,
            related_help_topic_ids=self.related_help_topic_ids,
            content_version=self.version_and_review.content_version,
            review_state=self.version_and_review.review_state,
        )


def validate_term_registry(terms: Iterable[TermSpec]) -> tuple[TermSpec, ...]:
    """Validate cross-term identity and search-key uniqueness.

    The returned tuple is sorted by stable concept ID, making downstream
    glossary and audit generation deterministic.
    """

    try:
        items = tuple(terms)
    except TypeError as exc:
        raise TerminologyValidationError("terms must be iterable") from exc
    if not items:
        raise TerminologyValidationError("A terminology registry cannot be empty")
    if any(not isinstance(term, TermSpec) for term in items):
        raise TerminologyValidationError(
            "A terminology registry may contain only TermSpec instances"
        )

    by_id: dict[str, TermSpec] = {}
    aliases: dict[str, str] = {}
    technical_labels: dict[str, str] = {}
    locale_key_owner: dict[str, str] = {}
    for term in items:
        if term.concept_id in by_id:
            raise TerminologyValidationError(
                f"Duplicate concept_id {term.concept_id!r}"
            )
        by_id[term.concept_id] = term

        for alias_key in term.aliases_keys:
            previous = aliases.get(alias_key)
            if previous is not None:
                raise TerminologyValidationError(
                    f"Alias key {alias_key!r} is shared by {previous!r} and "
                    f"{term.concept_id!r}"
                )
            aliases[alias_key] = term.concept_id

        if term.technical_label is not None:
            folded = term.technical_label.casefold()
            previous = technical_labels.get(folded)
            if previous is not None:
                raise TerminologyValidationError(
                    f"Technical label {term.technical_label!r} is shared by "
                    f"{previous!r} and {term.concept_id!r}"
                )
            technical_labels[folded] = term.concept_id

        for locale_key in term.locale_keys:
            previous = locale_key_owner.get(locale_key)
            if previous is not None:
                raise TerminologyValidationError(
                    f"Locale key {locale_key!r} is shared by {previous!r} and "
                    f"{term.concept_id!r}"
                )
            locale_key_owner[locale_key] = term.concept_id

    return tuple(sorted(items, key=lambda term: term.concept_id))


@dataclass(frozen=True)
class TerminologyRegistry:
    """Validated, deterministic collection of immutable term specifications."""

    version: str
    terms: tuple[TermSpec, ...]
    _by_id: Mapping[str, TermSpec] = field(init=False, repr=False, compare=False)
    _by_alias_key: Mapping[str, TermSpec] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        version = _require_text(self.version, "registry version")
        if not _REGISTRY_VERSION_RE.fullmatch(version):
            raise TerminologyValidationError(
                "registry version must look like 'mfrm_terminology_v1'"
            )
        terms = validate_term_registry(self.terms)
        by_id = MappingProxyType({term.concept_id: term for term in terms})
        by_alias = MappingProxyType(
            {
                alias_key: term
                for term in terms
                for alias_key in term.aliases_keys
            }
        )
        object.__setattr__(self, "version", version)
        object.__setattr__(self, "terms", terms)
        object.__setattr__(self, "_by_id", by_id)
        object.__setattr__(self, "_by_alias_key", by_alias)

    def __iter__(self) -> Iterator[TermSpec]:
        return iter(self.terms)

    def __len__(self) -> int:
        return len(self.terms)

    def get(self, concept_id: str) -> TermSpec | None:
        """Return a term by exact stable ID, or ``None`` when absent."""

        return self._by_id.get(concept_id)

    def require(self, concept_id: str) -> TermSpec:
        """Return a term by exact stable ID or raise ``UnknownConceptError``."""

        term = self.get(concept_id)
        if term is None:
            raise UnknownConceptError(concept_id)
        return term

    def resolve_alias_key(self, alias_key: str) -> TermSpec | None:
        """Resolve one exact localized alias key without substring matching."""

        return self._by_alias_key.get(alias_key)

    def project(
        self, concept_id: str, layer: AudienceLayer | str
    ) -> TermProjection:
        """Project a registered concept into one audience layer."""

        return self.require(concept_id).project(layer)

    def project_all(self, layer: AudienceLayer | str) -> tuple[TermProjection, ...]:
        """Project every term registered for ``layer`` in deterministic order."""

        normalized_layer = _coerce_layer(layer)
        return tuple(
            term.project(normalized_layer)
            for term in self.terms
            if term.allows(normalized_layer)
        )

    @property
    def locale_keys(self) -> frozenset[str]:
        """All locale keys required by the registry."""

        return frozenset(key for term in self.terms for key in term.locale_keys)


_INTERNAL_REVIEW = VersionReviewMetadata(
    content_version="1.0.0",
    owner_role="product_content_owner",
    last_reviewed=date(2026, 7, 24),
    review_state=ReviewState.IN_REVIEW,
    reviewer_roles=("ux_reviewer",),
)

_METHOD_REVIEW = VersionReviewMetadata(
    content_version="1.0.0",
    owner_role="measurement_method_owner",
    last_reviewed=date(2026, 7, 24),
    review_state=ReviewState.IN_REVIEW,
    reviewer_roles=("sla_researcher", "psychometric_reviewer"),
)


def _term(
    concept_id: str,
    *,
    standard: str | None = None,
    method: str | None = None,
    technical: str | None = None,
    definition: str,
    aliases: tuple[str, ...],
    prohibited_claims: tuple[str, ...],
    help_topics: tuple[str, ...],
    review: VersionReviewMetadata,
) -> TermSpec:
    layers = frozenset(
        layer
        for layer, label in (
            (AudienceLayer.STANDARD, standard),
            (AudienceLayer.METHOD, method),
            (AudienceLayer.TECHNICAL, technical),
        )
        if label is not None
    )
    return TermSpec(
        concept_id=concept_id,
        standard_label_key=standard,
        method_label_key=method,
        technical_label=technical,
        definition_key=definition,
        aliases_keys=aliases,
        allowed_layers=layers,
        prohibited_claim_keys=prohibited_claims,
        related_help_topic_ids=help_topics,
        version_and_review=review,
    )


TERM_SPECS: tuple[TermSpec, ...] = (
    _term(
        "download.contents_record",
        standard="terms.download_contents.standard_label",
        technical="manifest",
        definition="terms.download_contents.definition",
        aliases=(
            "terms.download_contents.alias.contents_check",
            "terms.download_contents.alias.manifest",
        ),
        prohibited_claims=("terms.download_contents.claim.not_privacy_guarantee",),
        help_topics=("help.downloads.privacy",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "estimation.jmle",
        method="terms.jmle.method_label",
        technical="JMLE",
        definition="terms.jmle.definition",
        aliases=("terms.jmle.alias.acronym", "terms.jmle.alias.full_name"),
        prohibited_claims=("terms.jmle.claim.not_automatically_unbiased",),
        help_topics=("help.methods.jmle_mml",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "estimation.mml",
        method="terms.mml.method_label",
        technical="MML",
        definition="terms.mml.definition",
        aliases=("terms.mml.alias.acronym", "terms.mml.alias.full_name"),
        prohibited_claims=("terms.mml.claim.not_automatically_preferred",),
        help_topics=("help.methods.jmle_mml",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "evidence.availability",
        standard="terms.evidence_availability.standard_label",
        technical="ComputationState",
        definition="terms.evidence_availability.definition",
        aliases=(
            "terms.evidence_availability.alias.status",
            "terms.evidence_availability.alias.computation_state",
        ),
        prohibited_claims=("terms.evidence_availability.claim.not_validity_proof",),
        help_topics=("help.report.claim_limits",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "evidence.limit_reason",
        standard="terms.limit_reason.standard_label",
        technical="ReasonCode",
        definition="terms.limit_reason.definition",
        aliases=(
            "terms.limit_reason.alias.reason",
            "terms.limit_reason.alias.reason_code",
        ),
        prohibited_claims=("terms.limit_reason.claim.not_person_judgment",),
        help_topics=("help.report.claim_limits",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "evidence.planned_check_stability",
        standard="terms.planned_check_stability.standard_label",
        technical="StabilityState",
        definition="terms.planned_check_stability.definition",
        aliases=(
            "terms.planned_check_stability.alias.stability",
            "terms.planned_check_stability.alias.stability_state",
        ),
        prohibited_claims=("terms.planned_check_stability.claim.not_general_robustness",),
        help_topics=("help.report.claim_limits",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "evidence.summary",
        standard="terms.evidence_summary.standard_label",
        technical="EvidenceRecord",
        definition="terms.evidence_summary.definition",
        aliases=(
            "terms.evidence_summary.alias.summary",
            "terms.evidence_summary.alias.evidence_record",
        ),
        prohibited_claims=("terms.evidence_summary.claim.not_sufficient_by_itself",),
        help_topics=("help.results.first_read",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "fit.infit",
        method="terms.infit.method_label",
        technical="Infit",
        definition="terms.infit.definition",
        aliases=("terms.infit.alias.information_weighted", "terms.infit.alias.mnsq"),
        prohibited_claims=("terms.infit.claim.not_automatic_exclusion",),
        help_topics=("help.results.fit",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "fit.outfit",
        method="terms.outfit.method_label",
        technical="Outfit",
        definition="terms.outfit.definition",
        aliases=("terms.outfit.alias.outlier_sensitive", "terms.outfit.alias.mnsq"),
        prohibited_claims=("terms.outfit.claim.not_automatic_exclusion",),
        help_topics=("help.results.fit",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "measurement.separation",
        method="terms.separation.method_label",
        technical="separation",
        definition="terms.separation.definition",
        aliases=(
            "terms.separation.alias.index",
            "terms.separation.alias.reliability",
        ),
        prohibited_claims=("terms.separation.claim.not_validity_proof",),
        help_topics=("help.results.measures_targeting",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "model.gpcm",
        method="terms.gpcm.method_label",
        technical="GPCM",
        definition="terms.gpcm.definition",
        aliases=("terms.gpcm.alias.acronym", "terms.gpcm.alias.full_name"),
        prohibited_claims=("terms.gpcm.claim.not_interchangeable",),
        help_topics=("help.methods.rsm_pcm_gpcm",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "model.mfrm",
        method="terms.mfrm.method_label",
        technical="MFRM",
        definition="terms.mfrm.definition",
        aliases=("terms.mfrm.alias.acronym", "terms.mfrm.alias.full_name"),
        prohibited_claims=("terms.mfrm.claim.not_validity_proof",),
        help_topics=("help.methods.rsm_pcm_gpcm",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "model.pcm",
        method="terms.pcm.method_label",
        technical="PCM",
        definition="terms.pcm.definition",
        aliases=("terms.pcm.alias.acronym", "terms.pcm.alias.full_name"),
        prohibited_claims=("terms.pcm.claim.not_interchangeable",),
        help_topics=("help.methods.rsm_pcm_gpcm",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "model.rsm",
        method="terms.rsm.method_label",
        technical="RSM",
        definition="terms.rsm.definition",
        aliases=("terms.rsm.alias.acronym", "terms.rsm.alias.full_name"),
        prohibited_claims=("terms.rsm.claim.not_interchangeable",),
        help_topics=("help.methods.rsm_pcm_gpcm",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "reproducibility.input_match",
        standard="terms.input_match.standard_label",
        technical="fingerprint",
        definition="terms.input_match.definition",
        aliases=(
            "terms.input_match.alias.match_value",
            "terms.input_match.alias.fingerprint",
        ),
        prohibited_claims=("terms.input_match.claim.not_data_quality_evidence",),
        help_topics=("help.downloads.repeat_analysis",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "reproducibility.script",
        standard="terms.reproducibility_script.standard_label",
        technical="reproducibility script",
        definition="terms.reproducibility_script.definition",
        aliases=(
            "terms.reproducibility_script.alias.rerun",
            "terms.reproducibility_script.alias.runner",
        ),
        prohibited_claims=("terms.reproducibility_script.claim.not_result_validation",),
        help_topics=("help.downloads.repeat_analysis",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "residual_structure.pca",
        method="terms.residual_pca.method_label",
        technical="residual PCA",
        definition="terms.residual_pca.definition",
        aliases=(
            "terms.residual_pca.alias.residual_pca",
            "terms.residual_pca.alias.pca_residuals",
        ),
        prohibited_claims=("terms.residual_pca.claim.not_unidimensionality_proof",),
        help_topics=("help.results.residual_structure",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "runtime.browser_retention",
        standard="terms.browser_retention.standard_label",
        technical="app session state",
        definition="terms.browser_retention.definition",
        aliases=(
            "terms.browser_retention.alias.browser_session",
            "terms.browser_retention.alias.session_storage",
        ),
        prohibited_claims=("terms.browser_retention.claim.not_durable_backup",),
        help_topics=("help.downloads.repeat_analysis",),
        review=_INTERNAL_REVIEW,
    ),
    _term(
        "scale.logit",
        method="terms.logit.method_label",
        technical="logit",
        definition="terms.logit.definition",
        aliases=("terms.logit.alias.scale", "terms.logit.alias.log_odds"),
        prohibited_claims=("terms.logit.claim.not_raw_score_unit",),
        help_topics=("help.results.measures_targeting",),
        review=_METHOD_REVIEW,
    ),
    _term(
        "workflow.analysis_reference",
        standard="terms.analysis_reference.standard_label",
        technical="AnalysisID",
        definition="terms.analysis_reference.definition",
        aliases=(
            "terms.analysis_reference.alias.reference",
            "terms.analysis_reference.alias.analysis_id",
        ),
        prohibited_claims=("terms.analysis_reference.claim.not_quality_marker",),
        help_topics=("help.downloads.repeat_analysis",),
        review=_INTERNAL_REVIEW,
    ),
)


TERMINOLOGY_REGISTRY = TerminologyRegistry(
    version=TERMINOLOGY_REGISTRY_VERSION,
    terms=TERM_SPECS,
)


INTERNAL_TECHNICAL_LABELS: Mapping[str, str] = MappingProxyType(
    {
        "workflow.analysis_reference": "AnalysisID",
        "evidence.summary": "EvidenceRecord",
        "evidence.limit_reason": "ReasonCode",
        "evidence.availability": "ComputationState",
        "evidence.planned_check_stability": "StabilityState",
        "reproducibility.script": "reproducibility script",
        "reproducibility.input_match": "fingerprint",
        "download.contents_record": "manifest",
        "runtime.browser_retention": "app session state",
    }
)


SCIENTIFIC_CONCEPT_IDS = frozenset(
    {
        "model.mfrm",
        "model.rsm",
        "model.pcm",
        "model.gpcm",
        "estimation.jmle",
        "estimation.mml",
        "scale.logit",
        "fit.infit",
        "fit.outfit",
        "residual_structure.pca",
        "measurement.separation",
    }
)


__all__ = [
    "AudienceLayer",
    "INTERNAL_TECHNICAL_LABELS",
    "LayerNotAllowedError",
    "ReviewState",
    "SCIENTIFIC_CONCEPT_IDS",
    "TERMINOLOGY_REGISTRY",
    "TERMINOLOGY_REGISTRY_VERSION",
    "TERM_SPECS",
    "TermProjection",
    "TermSpec",
    "TerminologyRegistry",
    "TerminologyValidationError",
    "UnknownConceptError",
    "VersionReviewMetadata",
    "validate_term_registry",
]
