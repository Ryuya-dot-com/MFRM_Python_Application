"""Pure, locale-independent state machine for the optional sample guide.

The guide is a presentation and learning projection over the production
workflow.  This module intentionally imports neither Streamlit nor pandas and
contains no estimator callback.  A guidance event may retain references to a
data set, fit, or evidence record, but it can never create or change scientific
state.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
import re


GUIDANCE_SCHEMA_VERSION = "mfrm_guidance_v1"
DEFAULT_SAMPLE_CONTEXT_ID = "sample.mfrm.writing.v1"


class GuidanceContractError(ValueError):
    """Raised when guide state or an event violates the guidance contract."""


class GuideLifecycle(str, Enum):
    NOT_STARTED = "NOT_STARTED"
    ACTIVE = "ACTIVE"
    SKIPPED = "SKIPPED"
    COMPLETED = "COMPLETED"


class GuideRoute(str, Enum):
    SAMPLE = "sample"
    OWN_DATA = "own_data"
    WITHOUT_GUIDE = "without_guide"


class GuideNode(str, Enum):
    WELCOME = "welcome"
    DATA_CHECK = "data_check"
    ESTIMATE = "estimate"
    EVIDENCE_REVIEW = "evidence_review"
    ARCHIVE = "archive"


class LearningState(str, Enum):
    LOCKED = "LOCKED"
    AVAILABLE = "AVAILABLE"
    ACTIVE = "ACTIVE"
    COMPLETE = "COMPLETE"


class GuidanceEventType(str, Enum):
    START_SAMPLE = "START_SAMPLE"
    START_OWN_DATA = "START_OWN_DATA"
    SKIP = "SKIP"
    EXIT = "EXIT"
    RESUME = "RESUME"
    RESTART = "RESTART"
    NODE_REVIEWED = "NODE_REVIEWED"
    FORMATIVE_ANSWERED = "FORMATIVE_ANSWERED"
    DATA_OR_MAPPING_CHANGED = "DATA_OR_MAPPING_CHANGED"
    DRAFT_SPEC_CHANGED = "DRAFT_SPEC_CHANGED"
    FIT_SUCCEEDED = "FIT_SUCCEEDED"
    FIT_FAILED = "FIT_FAILED"
    EVIDENCE_UPDATED = "EVIDENCE_UPDATED"
    ARCHIVE_CREATED = "ARCHIVE_CREATED"
    VIEW_CHANGED = "VIEW_CHANGED"
    LANGUAGE_CHANGED = "LANGUAGE_CHANGED"
    HELP_OPENED = "HELP_OPENED"
    HELP_RETURNED = "HELP_RETURNED"


@dataclass(frozen=True, slots=True)
class GuideNodeDefinition:
    """Stable catalog entry rendered through locale keys at the UI edge."""

    node_id: GuideNode
    ordinal: int
    title_key: str
    body_key: str
    stage_ids: tuple[str, ...]


GUIDE_NODES = (
    GuideNodeDefinition(
        GuideNode.WELCOME,
        1,
        "guide.welcome_title",
        "guide.welcome_body",
        ("plan",),
    ),
    GuideNodeDefinition(
        GuideNode.DATA_CHECK,
        2,
        "guide.data_check_title",
        "guide.data_check_body",
        ("plan", "check"),
    ),
    GuideNodeDefinition(
        GuideNode.ESTIMATE,
        3,
        "guide.estimate_title",
        "guide.estimate_body",
        ("estimate",),
    ),
    GuideNodeDefinition(
        GuideNode.EVIDENCE_REVIEW,
        4,
        "guide.evidence_review_title",
        "guide.evidence_review_body",
        ("diagnose", "stress", "interpret"),
    ),
    GuideNodeDefinition(
        GuideNode.ARCHIVE,
        5,
        "guide.archive_title",
        "guide.archive_body",
        ("decide", "archive"),
    ),
)
GUIDE_NODE_ORDER = tuple(item.node_id for item in GUIDE_NODES)
_NODE_INDEX = {node_id: index for index, node_id in enumerate(GUIDE_NODE_ORDER)}
_OPAQUE_REFERENCE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,255}$")
_PASSIVE_EVENTS = frozenset(
    {
        GuidanceEventType.VIEW_CHANGED,
        GuidanceEventType.LANGUAGE_CHANGED,
        GuidanceEventType.HELP_OPENED,
        GuidanceEventType.HELP_RETURNED,
    }
)


@dataclass(frozen=True, slots=True)
class GuideNodeProgress:
    node_id: GuideNode
    state: LearningState
    learning_evidence_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "node_id", _coerce_enum(GuideNode, self.node_id, "node"))
        object.__setattr__(self, "state", _coerce_enum(LearningState, self.state, "learning state"))
        evidence_ids = _references(self.learning_evidence_ids, "learning_evidence_ids")
        object.__setattr__(self, "learning_evidence_ids", evidence_ids)


@dataclass(frozen=True, slots=True)
class GuidanceState:
    """Session-scoped guide state with no scientific decision fields."""

    lifecycle: GuideLifecycle
    route_id: GuideRoute | None
    active_node_id: GuideNode | None
    last_valid_node_id: GuideNode
    node_states: tuple[GuideNodeProgress, ...]
    bound_study_context_id: str | None = None
    bound_data_fingerprint: str | None = None
    bound_analysis_id: str | None = None
    reviewed_evidence_ids: tuple[str, ...] = ()
    sample_context_id: str | None = None
    skip_origin_node_id: GuideNode | None = None
    workflow_complete: bool = False
    learning_complete: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "lifecycle", _coerce_enum(GuideLifecycle, self.lifecycle, "lifecycle"))
        if self.route_id is not None:
            object.__setattr__(self, "route_id", _coerce_enum(GuideRoute, self.route_id, "route"))
        if self.active_node_id is not None:
            object.__setattr__(self, "active_node_id", _coerce_enum(GuideNode, self.active_node_id, "active node"))
        object.__setattr__(self, "last_valid_node_id", _coerce_enum(GuideNode, self.last_valid_node_id, "last valid node"))
        if self.skip_origin_node_id is not None:
            object.__setattr__(self, "skip_origin_node_id", _coerce_enum(GuideNode, self.skip_origin_node_id, "skip origin"))

        states = tuple(self.node_states)
        if tuple(item.node_id for item in states) != GUIDE_NODE_ORDER:
            raise GuidanceContractError("node_states must contain every guide node exactly once in canonical order")
        object.__setattr__(self, "node_states", states)
        object.__setattr__(self, "reviewed_evidence_ids", _references(self.reviewed_evidence_ids, "reviewed_evidence_ids"))
        for name in (
            "bound_study_context_id",
            "bound_data_fingerprint",
            "bound_analysis_id",
            "sample_context_id",
        ):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _reference(value, name))

        if not isinstance(self.workflow_complete, bool) or not isinstance(self.learning_complete, bool):
            raise GuidanceContractError("completion fields must be booleans")
        active_nodes = tuple(item.node_id for item in states if item.state is LearningState.ACTIVE)
        if self.lifecycle in {GuideLifecycle.NOT_STARTED, GuideLifecycle.ACTIVE}:
            if active_nodes != (self.active_node_id,):
                raise GuidanceContractError("an open guide state requires exactly one active node")
        elif active_nodes or self.active_node_id is not None:
            raise GuidanceContractError("a skipped or completed guide cannot retain an active node")

        if self.lifecycle is GuideLifecycle.NOT_STARTED:
            if self.route_id is not None or self.active_node_id is not GuideNode.WELCOME:
                raise GuidanceContractError("the initial guide must be unbound at welcome")
        elif self.route_id is None:
            raise GuidanceContractError("a non-initial guide requires a route")

        if self.lifecycle is GuideLifecycle.SKIPPED and self.skip_origin_node_id is None:
            raise GuidanceContractError("a skipped guide requires its exit origin")
        if self.lifecycle is GuideLifecycle.COMPLETED:
            if not self.workflow_complete or not self.learning_complete:
                raise GuidanceContractError("a completed guide requires independent workflow and learning completion")
            if any(item.state is not LearningState.COMPLETE for item in states):
                raise GuidanceContractError("all nodes must be complete when the guide is completed")

        if self.route_id is GuideRoute.SAMPLE and self.sample_context_id is None:
            raise GuidanceContractError("the sample route requires an isolated sample context")
        if self.route_id is not GuideRoute.SAMPLE and self.sample_context_id is not None:
            raise GuidanceContractError("only the sample route may retain a sample context")

    @classmethod
    def initial(cls) -> "GuidanceState":
        return cls(
            lifecycle=GuideLifecycle.NOT_STARTED,
            route_id=None,
            active_node_id=GuideNode.WELCOME,
            last_valid_node_id=GuideNode.WELCOME,
            node_states=_initial_node_states(),
        )

    def node_progress(self, node_id: GuideNode | str) -> GuideNodeProgress:
        normalized = _coerce_enum(GuideNode, node_id, "node")
        return self.node_states[_NODE_INDEX[normalized]]

    def to_payload(self) -> dict[str, object]:
        """Return a deterministic, privacy-safe session payload."""

        return {
            "schema_version": GUIDANCE_SCHEMA_VERSION,
            "lifecycle": self.lifecycle.value,
            "route_id": self.route_id.value if self.route_id else None,
            "active_node_id": self.active_node_id.value if self.active_node_id else None,
            "last_valid_node_id": self.last_valid_node_id.value,
            "node_states": [
                {
                    "node_id": item.node_id.value,
                    "state": item.state.value,
                    "learning_evidence_ids": list(item.learning_evidence_ids),
                }
                for item in self.node_states
            ],
            "bound_study_context_id": self.bound_study_context_id,
            "bound_data_fingerprint": self.bound_data_fingerprint,
            "bound_analysis_id": self.bound_analysis_id,
            "reviewed_evidence_ids": list(self.reviewed_evidence_ids),
            "sample_context_id": self.sample_context_id,
            "skip_origin_node_id": self.skip_origin_node_id.value if self.skip_origin_node_id else None,
            "workflow_complete": self.workflow_complete,
            "learning_complete": self.learning_complete,
        }


@dataclass(frozen=True, slots=True)
class GuidanceEvent:
    event_type: GuidanceEventType
    data_fingerprint: str | None = None
    analysis_id: str | None = None
    study_context_id: str | None = None
    evidence_ids: tuple[str, ...] = ()
    learning_evidence_id: str | None = None
    sample_context_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "event_type", _coerce_enum(GuidanceEventType, self.event_type, "event"))
        for name in (
            "data_fingerprint",
            "analysis_id",
            "study_context_id",
            "learning_evidence_id",
            "sample_context_id",
        ):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _reference(value, name))
        object.__setattr__(self, "evidence_ids", _references(self.evidence_ids, "evidence_ids"))


def reduce_guidance(state: GuidanceState, event: GuidanceEvent) -> GuidanceState:
    """Apply one exhaustive guide event without touching scientific state."""

    if not isinstance(state, GuidanceState):
        raise TypeError("state must be a GuidanceState")
    if not isinstance(event, GuidanceEvent):
        raise TypeError("event must be a GuidanceEvent")
    kind = event.event_type
    if kind in _PASSIVE_EVENTS:
        return state
    if kind is GuidanceEventType.START_SAMPLE:
        return _start_route(
            GuideRoute.SAMPLE,
            sample_context_id=event.sample_context_id or DEFAULT_SAMPLE_CONTEXT_ID,
        )
    if kind is GuidanceEventType.START_OWN_DATA:
        return _start_route(GuideRoute.OWN_DATA)
    if kind in {GuidanceEventType.SKIP, GuidanceEventType.EXIT}:
        return _skip(state)
    if kind is GuidanceEventType.RESUME:
        return _resume(state)
    if kind is GuidanceEventType.RESTART:
        return _restart(state)
    if kind is GuidanceEventType.NODE_REVIEWED:
        return _complete_data_check(state, event)
    if kind is GuidanceEventType.FIT_SUCCEEDED:
        return _fit_succeeded(state, event)
    if kind is GuidanceEventType.FIT_FAILED:
        return _fit_failed(state, event)
    if kind is GuidanceEventType.FORMATIVE_ANSWERED:
        return _formative_answered(state, event)
    if kind is GuidanceEventType.ARCHIVE_CREATED:
        return _archive_created(state, event)
    if kind is GuidanceEventType.DATA_OR_MAPPING_CHANGED:
        return _invalidate_from(state, GuideNode.DATA_CHECK, data_fingerprint=event.data_fingerprint)
    if kind is GuidanceEventType.DRAFT_SPEC_CHANGED:
        return _invalidate_from(state, GuideNode.ESTIMATE)
    if kind is GuidanceEventType.EVIDENCE_UPDATED:
        return _invalidate_from(state, GuideNode.EVIDENCE_REVIEW, evidence_ids=event.evidence_ids)
    raise GuidanceContractError(f"unhandled guidance event: {kind.value}")


def guide_node_definition(node_id: GuideNode | str) -> GuideNodeDefinition:
    normalized = _coerce_enum(GuideNode, node_id, "node")
    return GUIDE_NODES[_NODE_INDEX[normalized]]


def _initial_node_states() -> tuple[GuideNodeProgress, ...]:
    return tuple(
        GuideNodeProgress(
            node_id=node_id,
            state=LearningState.ACTIVE if node_id is GuideNode.WELCOME else LearningState.LOCKED,
        )
        for node_id in GUIDE_NODE_ORDER
    )


def _start_route(route: GuideRoute, *, sample_context_id: str | None = None) -> GuidanceState:
    states = _states_through(GuideNode.DATA_CHECK, completed_before=True)
    return GuidanceState(
        lifecycle=GuideLifecycle.ACTIVE,
        route_id=route,
        active_node_id=GuideNode.DATA_CHECK,
        last_valid_node_id=GuideNode.WELCOME,
        node_states=states,
        sample_context_id=sample_context_id,
    )


def _skip(state: GuidanceState) -> GuidanceState:
    if state.lifecycle in {GuideLifecycle.SKIPPED, GuideLifecycle.COMPLETED}:
        return state
    origin = state.active_node_id or GuideNode.WELCOME
    route = state.route_id or GuideRoute.WITHOUT_GUIDE
    states = tuple(
        replace(item, state=LearningState.AVAILABLE)
        if item.state is LearningState.ACTIVE
        else item
        for item in state.node_states
    )
    return replace(
        state,
        lifecycle=GuideLifecycle.SKIPPED,
        route_id=route,
        active_node_id=None,
        node_states=states,
        skip_origin_node_id=origin,
    )


def _resume(state: GuidanceState) -> GuidanceState:
    if state.lifecycle is not GuideLifecycle.SKIPPED:
        return state
    if state.route_id is GuideRoute.WITHOUT_GUIDE:
        return state
    candidates = [
        item.node_id
        for item in state.node_states
        if item.state in {LearningState.AVAILABLE, LearningState.ACTIVE}
    ]
    node_id = candidates[0] if candidates else _first_incomplete_node(state)
    return replace(
        state,
        lifecycle=GuideLifecycle.ACTIVE,
        active_node_id=node_id,
        node_states=_activate(state.node_states, node_id),
        skip_origin_node_id=None,
    )


def _restart(state: GuidanceState) -> GuidanceState:
    if state.lifecycle is GuideLifecycle.NOT_STARTED:
        return state
    if state.route_id is GuideRoute.WITHOUT_GUIDE:
        return GuidanceState.initial()
    return GuidanceState(
        lifecycle=GuideLifecycle.ACTIVE,
        route_id=state.route_id,
        active_node_id=GuideNode.DATA_CHECK,
        last_valid_node_id=GuideNode.WELCOME,
        node_states=_states_through(GuideNode.DATA_CHECK, completed_before=True),
        bound_study_context_id=state.bound_study_context_id,
        bound_data_fingerprint=state.bound_data_fingerprint,
        bound_analysis_id=state.bound_analysis_id,
        reviewed_evidence_ids=(),
        sample_context_id=state.sample_context_id,
    )


def _complete_data_check(state: GuidanceState, event: GuidanceEvent) -> GuidanceState:
    _require_active(state, GuideNode.DATA_CHECK)
    evidence_id = event.learning_evidence_id or "learning.data_roles_reviewed"
    states = _complete_and_advance(state.node_states, GuideNode.DATA_CHECK, GuideNode.ESTIMATE, evidence_id)
    return replace(
        state,
        active_node_id=GuideNode.ESTIMATE,
        last_valid_node_id=GuideNode.DATA_CHECK,
        node_states=states,
        bound_data_fingerprint=event.data_fingerprint or state.bound_data_fingerprint,
        bound_study_context_id=event.study_context_id or state.bound_study_context_id,
    )


def _fit_succeeded(state: GuidanceState, event: GuidanceEvent) -> GuidanceState:
    _require_active(state, GuideNode.ESTIMATE)
    if event.analysis_id is None:
        raise GuidanceContractError("FIT_SUCCEEDED requires analysis_id")
    states = _complete_and_advance(
        state.node_states,
        GuideNode.ESTIMATE,
        GuideNode.EVIDENCE_REVIEW,
        event.learning_evidence_id or "learning.fit_outcome_reviewed",
    )
    return replace(
        state,
        active_node_id=GuideNode.EVIDENCE_REVIEW,
        last_valid_node_id=GuideNode.ESTIMATE,
        node_states=states,
        bound_analysis_id=event.analysis_id,
        bound_data_fingerprint=event.data_fingerprint or state.bound_data_fingerprint,
    )


def _fit_failed(state: GuidanceState, event: GuidanceEvent) -> GuidanceState:
    _require_active(state, GuideNode.ESTIMATE)
    evidence_id = event.learning_evidence_id or "learning.fit_failure_reviewed"
    current = state.node_progress(GuideNode.ESTIMATE)
    updated = replace(
        current,
        learning_evidence_ids=_merge_references(current.learning_evidence_ids, (evidence_id,)),
    )
    states = _replace_progress(state.node_states, updated)
    return replace(state, node_states=states, bound_analysis_id=None)


def _formative_answered(state: GuidanceState, event: GuidanceEvent) -> GuidanceState:
    _require_active(state, GuideNode.EVIDENCE_REVIEW)
    if event.learning_evidence_id is None:
        raise GuidanceContractError("FORMATIVE_ANSWERED requires learning_evidence_id")
    states = _complete_and_advance(
        state.node_states,
        GuideNode.EVIDENCE_REVIEW,
        GuideNode.ARCHIVE,
        event.learning_evidence_id,
    )
    return replace(
        state,
        active_node_id=GuideNode.ARCHIVE,
        last_valid_node_id=GuideNode.EVIDENCE_REVIEW,
        node_states=states,
        reviewed_evidence_ids=_merge_references(state.reviewed_evidence_ids, event.evidence_ids),
    )


def _archive_created(state: GuidanceState, event: GuidanceEvent) -> GuidanceState:
    _require_active(state, GuideNode.ARCHIVE)
    evidence_id = event.learning_evidence_id or "learning.sample_record_reviewed"
    archive = state.node_progress(GuideNode.ARCHIVE)
    completed_archive = replace(
        archive,
        state=LearningState.COMPLETE,
        learning_evidence_ids=_merge_references(archive.learning_evidence_ids, (evidence_id,)),
    )
    return replace(
        state,
        lifecycle=GuideLifecycle.COMPLETED,
        active_node_id=None,
        last_valid_node_id=GuideNode.ARCHIVE,
        node_states=_replace_progress(state.node_states, completed_archive),
        workflow_complete=True,
        learning_complete=True,
    )


def _invalidate_from(
    state: GuidanceState,
    node_id: GuideNode,
    *,
    data_fingerprint: str | None = None,
    evidence_ids: tuple[str, ...] = (),
) -> GuidanceState:
    if state.lifecycle is GuideLifecycle.NOT_STARTED:
        return state
    if state.route_id is GuideRoute.WITHOUT_GUIDE:
        return state
    states: list[GuideNodeProgress] = []
    start = _NODE_INDEX[node_id]
    for index, item in enumerate(state.node_states):
        if index < start:
            states.append(item)
        elif index == start:
            states.append(GuideNodeProgress(item.node_id, LearningState.ACTIVE))
        else:
            states.append(GuideNodeProgress(item.node_id, LearningState.LOCKED))
    return replace(
        state,
        lifecycle=GuideLifecycle.ACTIVE,
        active_node_id=node_id,
        last_valid_node_id=GUIDE_NODE_ORDER[max(0, start - 1)],
        node_states=tuple(states),
        bound_data_fingerprint=data_fingerprint or state.bound_data_fingerprint,
        bound_analysis_id=None if start <= _NODE_INDEX[GuideNode.ESTIMATE] else state.bound_analysis_id,
        reviewed_evidence_ids=tuple(evidence_ids),
        skip_origin_node_id=None,
        workflow_complete=False,
        learning_complete=False,
    )


def _states_through(active_node: GuideNode, *, completed_before: bool) -> tuple[GuideNodeProgress, ...]:
    active_index = _NODE_INDEX[active_node]
    return tuple(
        GuideNodeProgress(
            node_id=node_id,
            state=(
                LearningState.COMPLETE
                if completed_before and index < active_index
                else LearningState.ACTIVE
                if index == active_index
                else LearningState.LOCKED
            ),
        )
        for index, node_id in enumerate(GUIDE_NODE_ORDER)
    )


def _complete_and_advance(
    states: tuple[GuideNodeProgress, ...],
    current_node: GuideNode,
    next_node: GuideNode,
    learning_evidence_id: str,
) -> tuple[GuideNodeProgress, ...]:
    current = states[_NODE_INDEX[current_node]]
    completed = replace(
        current,
        state=LearningState.COMPLETE,
        learning_evidence_ids=_merge_references(current.learning_evidence_ids, (learning_evidence_id,)),
    )
    output = _replace_progress(states, completed)
    return _activate(output, next_node)


def _activate(states: tuple[GuideNodeProgress, ...], node_id: GuideNode) -> tuple[GuideNodeProgress, ...]:
    return tuple(
        replace(item, state=LearningState.ACTIVE)
        if item.node_id is node_id
        else replace(item, state=LearningState.AVAILABLE)
        if item.state is LearningState.ACTIVE
        else item
        for item in states
    )


def _replace_progress(
    states: tuple[GuideNodeProgress, ...], updated: GuideNodeProgress
) -> tuple[GuideNodeProgress, ...]:
    return tuple(updated if item.node_id is updated.node_id else item for item in states)


def _first_incomplete_node(state: GuidanceState) -> GuideNode:
    for item in state.node_states:
        if item.state is not LearningState.COMPLETE:
            return item.node_id
    return GuideNode.ARCHIVE


def _require_active(state: GuidanceState, node_id: GuideNode) -> None:
    if state.lifecycle is not GuideLifecycle.ACTIVE or state.active_node_id is not node_id:
        raise GuidanceContractError(f"{node_id.value} is not the active guide node")


def _merge_references(existing: tuple[str, ...], added: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(dict.fromkeys((*existing, *added)))


def _reference(value: object, name: str) -> str:
    if not isinstance(value, str) or value != value.strip() or not _OPAQUE_REFERENCE_RE.fullmatch(value):
        raise GuidanceContractError(f"{name} must be an opaque, privacy-safe reference")
    return value


def _references(values: tuple[str, ...], name: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes, bytearray)):
        raise GuidanceContractError(f"{name} must be a sequence")
    normalized = tuple(_reference(value, f"{name}[{index}]") for index, value in enumerate(values))
    if len(normalized) != len(set(normalized)):
        raise GuidanceContractError(f"{name} must be unique")
    return normalized


def _coerce_enum(enum_type: type[Enum], value: object, name: str):
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(str(value))
    except ValueError as exc:
        raise GuidanceContractError(f"unknown {name}: {value!r}") from exc
