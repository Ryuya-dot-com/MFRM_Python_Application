"""Pure reducer for locale-independent Help navigation.

Only presentation-route identities and optional context references enter this
module.  There is no Streamlit session access and no field or argument for an
estimate, evidence conclusion, scientific readiness, decision, or learning
completion.  Adapters may therefore open, browse, invalidate, and return from
Help without changing those states.
"""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from enum import Enum

from .help_contract import (
    HelpContextBinding,
    HelpContextPolicy,
    HelpRegistry,
    _coerce_enum,
    _require_stable_id,
)


SAFE_FALLBACK_TOPIC_ID = "help.fallback.unavailable"
SAFE_FALLBACK_SECTION_ID = "overview"


class _RouteConstructionPermit:
    """One-use permit preventing public construction or dataclass replacement."""

    __slots__ = ("active",)

    def __init__(self) -> None:
        self.active = True


class HelpContextStatus(str, Enum):
    """Whether a dynamic current-analysis Help block may be rendered."""

    CURRENT = "CURRENT"
    STATIC_ONLY = "STATIC_ONLY"
    STALE = "STALE"
    MISSING_REQUIRED = "MISSING_REQUIRED"


class HelpFallbackReason(str, Enum):
    """Coarse, localizable reasons for a safe generic Help route."""

    LINK_UNAVAILABLE = "help.link_unavailable"
    TOPIC_UNAVAILABLE = "help.topic_unavailable"
    SECTION_UNAVAILABLE = "help.section_unavailable"
    RETURN_TARGET_UNAVAILABLE = "help.return_target_unavailable"


@dataclass(frozen=True)
class HelpRouteLocation:
    """One locale-independent topic/section pair in route history."""

    help_topic_id: str
    section_id: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "help_topic_id",
            _require_stable_id(self.help_topic_id, "help_topic_id"),
        )
        object.__setattr__(
            self, "section_id", _require_stable_id(self.section_id, "section_id")
        )


@dataclass(frozen=True)
class HelpRouteState:
    """Session-scoped Help presentation state.

    A closed state may retain the exact return target/focus so the rendering
    adapter can restore it.  Topic identities and open-route history remain
    stable across locale changes because no translated title or locale code is
    stored here; closing Help clears that route history.
    """

    is_open: bool = False
    help_topic_id: str | None = None
    section_id: str | None = None
    help_link_id: str | None = None
    source_target_id: str | None = None
    return_target_id: str | None = None
    return_focus_id: str | None = None
    context_policy: HelpContextPolicy = HelpContextPolicy.STATIC_ONLY
    context_binding: HelpContextBinding | None = None
    context_status: HelpContextStatus = HelpContextStatus.STATIC_ONLY
    history: tuple[HelpRouteLocation, ...] = ()
    is_fallback: bool = False
    fallback_reason: HelpFallbackReason | None = None
    _construction_permit: InitVar[object] = None

    def __post_init__(self, _construction_permit: object) -> None:
        permit = _construction_permit
        if not isinstance(permit, _RouteConstructionPermit) or not permit.active:
            raise ValueError(
                "HelpRouteState must be created by the Help navigation reducers"
            )
        permit.active = False
        if not isinstance(self.is_open, bool):
            raise TypeError("is_open must be a boolean")
        for field_name in (
            "help_topic_id",
            "section_id",
            "help_link_id",
            "source_target_id",
            "return_target_id",
            "return_focus_id",
        ):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(
                    self,
                    field_name,
                    _require_stable_id(value, field_name),
                )
        object.__setattr__(
            self,
            "context_policy",
            _coerce_enum(HelpContextPolicy, self.context_policy, "context policy"),
        )
        object.__setattr__(
            self,
            "context_status",
            _coerce_enum(HelpContextStatus, self.context_status, "context status"),
        )
        if self.context_binding is not None and not isinstance(
            self.context_binding, HelpContextBinding
        ):
            raise TypeError("context_binding must be a HelpContextBinding or None")
        if isinstance(self.history, (str, bytes, bytearray)):
            raise TypeError("history must be a sequence of HelpRouteLocation values")
        history = tuple(self.history)
        if any(not isinstance(item, HelpRouteLocation) for item in history):
            raise TypeError("history must contain only HelpRouteLocation values")
        object.__setattr__(self, "history", history)
        if not isinstance(self.is_fallback, bool):
            raise TypeError("is_fallback must be a boolean")
        if self.fallback_reason is not None:
            object.__setattr__(
                self,
                "fallback_reason",
                _coerce_enum(
                    HelpFallbackReason, self.fallback_reason, "fallback reason"
                ),
            )

        if self.is_open and (self.help_topic_id is None or self.section_id is None):
            raise ValueError("An open Help route requires a topic and section")
        if not self.is_open and (
            self.help_topic_id is not None
            or self.section_id is not None
            or self.help_link_id is not None
            or self.context_binding is not None
            or self.history
            or self.is_fallback
            or self.fallback_reason is not None
        ):
            raise ValueError("A closed Help route cannot retain active Help content")
        if not self.is_open and (
            self.context_policy is not HelpContextPolicy.STATIC_ONLY
            or self.context_status is not HelpContextStatus.STATIC_ONLY
        ):
            raise ValueError("A closed Help route must keep neutral static context")
        if (self.return_target_id is None) != (self.return_focus_id is None):
            raise ValueError("return_target_id and return_focus_id must appear together")
        if self.context_policy is HelpContextPolicy.STATIC_ONLY:
            if self.context_binding is not None:
                raise ValueError("STATIC_ONLY Help cannot retain a context binding")
            if self.context_status is not HelpContextStatus.STATIC_ONLY:
                raise ValueError("STATIC_ONLY policy requires STATIC_ONLY status")
        elif self.context_binding is None:
            expected_without_binding = (
                HelpContextStatus.STATIC_ONLY
                if self.context_policy is HelpContextPolicy.OPTIONAL_CURRENT
                else HelpContextStatus.MISSING_REQUIRED
            )
            if self.context_status is not expected_without_binding:
                raise ValueError(
                    "Help context status is inconsistent with a missing binding"
                )
        elif self.context_status is HelpContextStatus.STATIC_ONLY:
            raise ValueError(
                "A non-static Help route with a binding cannot use STATIC_ONLY status"
            )
        if self.is_fallback != (self.fallback_reason is not None):
            raise ValueError("Fallback routes require exactly one fallback reason")
        if self.is_fallback and (
            self.help_topic_id != SAFE_FALLBACK_TOPIC_ID
            or self.section_id != SAFE_FALLBACK_SECTION_ID
        ):
            raise ValueError("Fallback routes must use the reserved safe destination")
        if (
            self.is_open
            and self.help_topic_id == SAFE_FALLBACK_TOPIC_ID
            and not self.is_fallback
        ):
            raise ValueError(
                "The reserved safe destination must be a fallback route"
            )
        if self.is_fallback and (
            self.context_policy is not HelpContextPolicy.STATIC_ONLY
            or self.context_binding is not None
            or self.context_status is not HelpContextStatus.STATIC_ONLY
        ):
            raise ValueError("Fallback routes must keep neutral static context")

    @classmethod
    def closed(cls) -> "HelpRouteState":
        """Return the initial route state."""

        return _new_route_state()

    @property
    def return_destination(self) -> tuple[str, str] | None:
        if self.return_target_id is None or self.return_focus_id is None:
            return None
        return self.return_target_id, self.return_focus_id


def evaluate_context_binding(
    context_policy: HelpContextPolicy | str,
    binding: HelpContextBinding | None,
    current_context: HelpContextBinding | None,
) -> HelpContextStatus:
    """Classify whether a current-analysis Help projection is safe to show."""

    policy = _coerce_enum(HelpContextPolicy, context_policy, "context policy")
    if policy is HelpContextPolicy.STATIC_ONLY:
        return HelpContextStatus.STATIC_ONLY
    if binding is None:
        if policy is HelpContextPolicy.OPTIONAL_CURRENT:
            return HelpContextStatus.STATIC_ONLY
        return HelpContextStatus.MISSING_REQUIRED
    if not isinstance(binding, HelpContextBinding):
        raise TypeError("binding must be a HelpContextBinding or None")
    if current_context is not None and not isinstance(
        current_context, HelpContextBinding
    ):
        raise TypeError("current_context must be a HelpContextBinding or None")

    # Any dynamic statement must be tied to a particular analysis.  A binding
    # with only sample/real provenance is useful for isolation but insufficient
    # for presenting current-analysis findings.
    if not binding.has_current_analysis:
        return HelpContextStatus.MISSING_REQUIRED
    if current_context is None or not current_context.has_current_analysis:
        return HelpContextStatus.MISSING_REQUIRED

    if binding.data_context is not current_context.data_context:
        return HelpContextStatus.STALE
    if binding.analysis_id != current_context.analysis_id:
        return HelpContextStatus.STALE
    if not _optional_reference_matches(
        binding.study_context_id, current_context.study_context_id
    ):
        return HelpContextStatus.STALE
    if not _optional_reference_matches(
        binding.evidence_issue_id, current_context.evidence_issue_id
    ):
        return HelpContextStatus.STALE
    if (
        binding.analysis_phase is not None
        and binding.analysis_phase is not current_context.analysis_phase
    ):
        return HelpContextStatus.STALE
    if not set(binding.evidence_ids).issubset(current_context.evidence_ids):
        return HelpContextStatus.STALE
    if not set(binding.claim_boundary_ids).issubset(
        current_context.claim_boundary_ids
    ):
        return HelpContextStatus.STALE
    return HelpContextStatus.CURRENT


def _optional_reference_matches(bound: str | None, current: str | None) -> bool:
    return bound is None or bound == current


def open_help_link(
    state: HelpRouteState,
    registry: HelpRegistry,
    help_link_id: str,
    *,
    context_binding: HelpContextBinding | None = None,
    current_context: HelpContextBinding | None = None,
    fallback_return_target_id: str | None = None,
    fallback_return_focus_id: str | None = None,
) -> HelpRouteState:
    """Open one registered link or a visible safe fallback.

    Reopening the same link with the same verified context is idempotent and
    does not grow route history.
    """

    _require_state_and_registry(state, registry)
    try:
        normalized_link_id = _require_stable_id(help_link_id, "help_link_id")
    except Exception:
        return safe_fallback_route(
            state,
            HelpFallbackReason.LINK_UNAVAILABLE,
            registry=registry,
            return_target_id=fallback_return_target_id,
            return_focus_id=fallback_return_focus_id,
        )
    link = registry.link(normalized_link_id)
    if link is None:
        return safe_fallback_route(
            state,
            HelpFallbackReason.LINK_UNAVAILABLE,
            registry=registry,
            return_target_id=fallback_return_target_id,
            return_focus_id=fallback_return_focus_id,
        )
    if link.help_topic_id == SAFE_FALLBACK_TOPIC_ID:
        return safe_fallback_route(
            state,
            HelpFallbackReason.LINK_UNAVAILABLE,
            registry=registry,
            return_target_id=fallback_return_target_id,
            return_focus_id=fallback_return_focus_id,
        )

    topic = registry.topic(link.help_topic_id)
    if topic is None:  # validated registries prevent this; keep runtime fail-closed
        return safe_fallback_route(
            state, HelpFallbackReason.TOPIC_UNAVAILABLE, registry=registry
        )
    section_id = link.section_id or topic.sections[0].section_id
    if not topic.has_section(section_id):
        return safe_fallback_route(
            state, HelpFallbackReason.SECTION_UNAVAILABLE, registry=registry
        )

    effective_binding = (
        None
        if link.context_policy is HelpContextPolicy.STATIC_ONLY
        else context_binding
    )
    context_status = evaluate_context_binding(
        link.context_policy, effective_binding, current_context
    )
    history = _history_before_transition(state, link.help_topic_id, section_id)
    candidate = _new_route_state(
        is_open=True,
        help_topic_id=link.help_topic_id,
        section_id=section_id,
        help_link_id=link.help_link_id,
        source_target_id=link.source_target_id,
        return_target_id=link.return_target_id,
        return_focus_id=link.return_focus_id,
        context_policy=link.context_policy,
        context_binding=effective_binding,
        context_status=context_status,
        history=history,
    )
    return state if candidate == state else candidate


def change_help_topic(
    state: HelpRouteState,
    registry: HelpRegistry,
    help_topic_id: str,
    *,
    section_id: str | None = None,
    context_policy: HelpContextPolicy | str = HelpContextPolicy.STATIC_ONLY,
    context_binding: HelpContextBinding | None = None,
    current_context: HelpContextBinding | None = None,
) -> HelpRouteState:
    """Change an active Help topic while preserving its origin and return.

    A delayed topic-change event after Help closes is ignored; only an
    explicit open-link event may start a new Help session.
    """

    _require_state_and_registry(state, registry)
    if not state.is_open:
        return state
    try:
        normalized_topic_id = _require_stable_id(help_topic_id, "help_topic_id")
    except Exception:
        return safe_fallback_route(
            state, HelpFallbackReason.TOPIC_UNAVAILABLE, registry=registry
        )
    topic = registry.topic(normalized_topic_id)
    if topic is None:
        return safe_fallback_route(
            state, HelpFallbackReason.TOPIC_UNAVAILABLE, registry=registry
        )
    if normalized_topic_id == SAFE_FALLBACK_TOPIC_ID:
        return safe_fallback_route(
            state, HelpFallbackReason.TOPIC_UNAVAILABLE, registry=registry
        )
    if section_id is None:
        normalized_section_id = topic.sections[0].section_id
    else:
        try:
            normalized_section_id = _require_stable_id(section_id, "section_id")
        except Exception:
            return safe_fallback_route(
                state, HelpFallbackReason.SECTION_UNAVAILABLE, registry=registry
            )
    if not topic.has_section(normalized_section_id):
        return safe_fallback_route(
            state, HelpFallbackReason.SECTION_UNAVAILABLE, registry=registry
        )

    policy = _coerce_enum(HelpContextPolicy, context_policy, "context policy")
    effective_binding = (
        None if policy is HelpContextPolicy.STATIC_ONLY else context_binding
    )
    status = evaluate_context_binding(policy, effective_binding, current_context)
    history = _history_before_transition(
        state, normalized_topic_id, normalized_section_id
    )
    candidate = _new_route_state(
        is_open=True,
        help_topic_id=normalized_topic_id,
        section_id=normalized_section_id,
        help_link_id=None,
        source_target_id=state.source_target_id,
        return_target_id=state.return_target_id,
        return_focus_id=state.return_focus_id,
        context_policy=policy,
        context_binding=effective_binding,
        context_status=status,
        history=history,
    )
    return state if candidate == state else candidate


def return_from_help(state: HelpRouteState) -> HelpRouteState:
    """Close Help while retaining only the presentation return destination."""

    if not isinstance(state, HelpRouteState):
        raise TypeError("state must be a HelpRouteState")
    if not state.is_open:
        return state
    return _new_route_state(
        is_open=False,
        source_target_id=state.source_target_id,
        return_target_id=state.return_target_id,
        return_focus_id=state.return_focus_id,
        history=(),
    )


def invalidate_help_route(
    state: HelpRouteState,
    registry: HelpRegistry,
    *,
    current_context: HelpContextBinding | None = None,
) -> HelpRouteState:
    """Revalidate an active route after registry or analysis identity changes."""

    _require_state_and_registry(state, registry)
    if not state.is_open:
        return state
    if state.is_fallback:
        return safe_fallback_route(
            state,
            state.fallback_reason or HelpFallbackReason.TOPIC_UNAVAILABLE,
            registry=registry,
        )
    topic = registry.topic(state.help_topic_id or "")
    if topic is None:
        return safe_fallback_route(
            state, HelpFallbackReason.TOPIC_UNAVAILABLE, registry=registry
        )
    if state.section_id is None or not topic.has_section(state.section_id):
        return safe_fallback_route(
            state, HelpFallbackReason.SECTION_UNAVAILABLE, registry=registry
        )
    if state.help_link_id is not None:
        link = registry.link(state.help_link_id)
        linked_topic = None if link is None else registry.topic(link.help_topic_id)
        resolved_link_section = (
            None
            if link is None or linked_topic is None
            else link.section_id or linked_topic.sections[0].section_id
        )
        if link is None or linked_topic is None or (
            link.help_topic_id != state.help_topic_id
            or resolved_link_section != state.section_id
            or link.source_target_id != state.source_target_id
            or link.return_target_id != state.return_target_id
            or link.return_focus_id != state.return_focus_id
            or link.context_policy is not state.context_policy
        ):
            return safe_fallback_route(
                state, HelpFallbackReason.LINK_UNAVAILABLE, registry=registry
            )
    if state.return_target_id is not None:
        target = registry.target(state.return_target_id)
        if (
            target is None
            or state.return_focus_id is None
            or state.return_focus_id not in target.focus_ids
        ):
            return safe_fallback_route(
                state,
                HelpFallbackReason.RETURN_TARGET_UNAVAILABLE,
                registry=registry,
            )
    status = evaluate_context_binding(
        state.context_policy, state.context_binding, current_context
    )
    if status is state.context_status:
        return state
    return _new_route_state(
        is_open=state.is_open,
        help_topic_id=state.help_topic_id,
        section_id=state.section_id,
        help_link_id=state.help_link_id,
        source_target_id=state.source_target_id,
        return_target_id=state.return_target_id,
        return_focus_id=state.return_focus_id,
        context_policy=state.context_policy,
        context_binding=state.context_binding,
        context_status=status,
        history=state.history,
        is_fallback=state.is_fallback,
        fallback_reason=state.fallback_reason,
    )


def retain_route_for_locale_change(
    state: HelpRouteState, *, locale: str | None = None
) -> HelpRouteState:
    """Retain the exact route when presentation locale changes.

    ``locale`` is accepted so adapters can call this at their locale boundary,
    but it is intentionally neither interpreted nor persisted as route
    identity.
    """

    if not isinstance(state, HelpRouteState):
        raise TypeError("state must be a HelpRouteState")
    if locale is not None and not isinstance(locale, str):
        raise TypeError("locale must be a string or None")
    return state


def safe_fallback_route(
    state: HelpRouteState,
    reason: HelpFallbackReason | str,
    *,
    registry: HelpRegistry,
    return_target_id: str | None = None,
    return_focus_id: str | None = None,
) -> HelpRouteState:
    """Render-safe route for a missing link/topic/section; never a silent no-op."""

    if not isinstance(state, HelpRouteState):
        raise TypeError("state must be a HelpRouteState")
    if not isinstance(registry, HelpRegistry):
        raise TypeError("registry must be a HelpRegistry")
    fallback_reason = _coerce_enum(HelpFallbackReason, reason, "fallback reason")

    fallback_topic = registry.topic(SAFE_FALLBACK_TOPIC_ID)
    if (
        fallback_topic is None
        or not fallback_topic.has_section(SAFE_FALLBACK_SECTION_ID)
        or SAFE_FALLBACK_TOPIC_ID not in registry.system_topic_ids
    ):
        raise ValueError(
            "registry must reserve help.fallback.unavailable#overview as a "
            "system topic before a safe fallback route can be created"
        )

    selected_source = state.source_target_id if state.is_open else None
    selected_target = state.return_target_id if state.is_open else None
    selected_focus = state.return_focus_id if state.is_open else None
    if return_target_id is not None or return_focus_id is not None:
        if return_target_id is None or return_focus_id is None:
            raise ValueError("fallback return target and focus must appear together")
        selected_target = _require_stable_id(
            return_target_id, "fallback return_target_id"
        )
        selected_focus = _require_stable_id(
            return_focus_id, "fallback return_focus_id"
        )
        selected_source = selected_target

    if selected_target is not None:
        target = registry.target(selected_target)
        if target is None or selected_focus not in target.focus_ids:
            selected_source = None
            selected_target = None
            selected_focus = None

    history = _history_before_transition(
        state, SAFE_FALLBACK_TOPIC_ID, SAFE_FALLBACK_SECTION_ID
    )
    candidate = _new_route_state(
        is_open=True,
        help_topic_id=SAFE_FALLBACK_TOPIC_ID,
        section_id=SAFE_FALLBACK_SECTION_ID,
        source_target_id=selected_source,
        return_target_id=selected_target,
        return_focus_id=selected_focus,
        context_policy=HelpContextPolicy.STATIC_ONLY,
        context_status=HelpContextStatus.STATIC_ONLY,
        history=history,
        is_fallback=True,
        fallback_reason=fallback_reason,
    )
    return state if candidate == state else candidate


def _history_before_transition(
    state: HelpRouteState, help_topic_id: str, section_id: str
) -> tuple[HelpRouteLocation, ...]:
    if not state.is_open or state.help_topic_id is None or state.section_id is None:
        return ()
    current = HelpRouteLocation(state.help_topic_id, state.section_id)
    destination = HelpRouteLocation(help_topic_id, section_id)
    if current == destination:
        return state.history
    if state.history and state.history[-1] == current:
        return state.history
    return (*state.history, current)


def _new_route_state(**values: object) -> HelpRouteState:
    """Construct one validated state with a one-use reducer permit."""

    return HelpRouteState(
        **values,
        _construction_permit=_RouteConstructionPermit(),
    )


def _require_state_and_registry(
    state: HelpRouteState, registry: HelpRegistry
) -> None:
    if not isinstance(state, HelpRouteState):
        raise TypeError("state must be a HelpRouteState")
    if not isinstance(registry, HelpRegistry):
        raise TypeError("registry must be a HelpRegistry")


# Short reducer aliases for adapters that prefer event-style names.
open_help = open_help_link
change_help = change_help_topic
return_help = return_from_help
invalidate_help = invalidate_help_route


__all__ = [
    "SAFE_FALLBACK_SECTION_ID",
    "SAFE_FALLBACK_TOPIC_ID",
    "HelpContextStatus",
    "HelpFallbackReason",
    "HelpRouteLocation",
    "HelpRouteState",
    "change_help",
    "change_help_topic",
    "evaluate_context_binding",
    "invalidate_help",
    "invalidate_help_route",
    "open_help",
    "open_help_link",
    "retain_route_for_locale_change",
    "return_from_help",
    "return_help",
    "safe_fallback_route",
]
