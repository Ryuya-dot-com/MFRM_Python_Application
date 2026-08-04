from __future__ import annotations

from dataclasses import fields, replace
import inspect

import pytest

from mfrm_app import help_navigation as navigation
from mfrm_app.help_contract import (
    HelpContextBinding,
    HelpLink,
    HelpRegistry,
    HelpSection,
    HelpTarget,
    HelpTopic,
)
from mfrm_app.help_navigation import (
    SAFE_FALLBACK_SECTION_ID,
    SAFE_FALLBACK_TOPIC_ID,
    HelpContextStatus,
    HelpFallbackReason,
    HelpRouteState,
    change_help_topic,
    evaluate_context_binding,
    invalidate_help_route,
    open_help_link,
    retain_route_for_locale_change,
    return_from_help,
)


def _topic(help_topic_id: str, target_id: str) -> HelpTopic:
    key_root = help_topic_id.replace("help.", "help_topics.")
    return HelpTopic(
        help_topic_id=help_topic_id,
        title_key=f"{key_root}.title",
        summary_key=f"{key_root}.summary",
        sections=(
            HelpSection(
                section_id="overview",
                title_key=f"{key_root}.overview_title",
                content_keys=(f"{key_root}.overview_body",),
            ),
            HelpSection(
                section_id="limits",
                title_key=f"{key_root}.limits_title",
                content_keys=(f"{key_root}.limits_body",),
            ),
        ),
        version_and_review={
            "content_version": "1.0.0",
            "owner_role": "help_content_owner",
            "review_state": "in_review",
            "last_reviewed": "2026-07-24",
            "reviewer_roles": ("sla_researcher", "ux_reviewer"),
        },
        claim_boundary_ids=("claim.interpretation_bounded",),
        reference_ids=("reference.help_method",),
        related_target_ids=(target_id,),
    )


def _registry() -> HelpRegistry:
    fit_target = HelpTarget(
        target_id="screen.results.fit",
        surface_id="surface.results",
        focus_ids=("focus.results.fit_chart",),
    )
    report_target = HelpTarget(
        target_id="screen.report.claims",
        surface_id="surface.report",
        focus_ids=("focus.report.claims",),
    )
    fit_topic = _topic("help.results.fit", fit_target.target_id)
    report_topic = _topic("help.report.claim_limits", report_target.target_id)
    fallback_topic = _topic(SAFE_FALLBACK_TOPIC_ID, fit_target.target_id)
    return HelpRegistry(
        topics=(fit_topic, report_topic, fallback_topic),
        links=(
            HelpLink(
                help_link_id="link.results.fit_chart",
                source_target_id=fit_target.target_id,
                help_topic_id=fit_topic.help_topic_id,
                section_id="limits",
                return_target_id=fit_target.target_id,
                return_focus_id="focus.results.fit_chart",
                context_policy="required_current",
            ),
            HelpLink(
                help_link_id="link.report.claims",
                source_target_id=report_target.target_id,
                help_topic_id=report_topic.help_topic_id,
                section_id="overview",
                return_target_id=report_target.target_id,
                return_focus_id="focus.report.claims",
                context_policy="static_only",
            ),
        ),
        targets=(fit_target, report_target),
        system_topic_ids=(SAFE_FALLBACK_TOPIC_ID,),
    )


def _context(
    *,
    data_context: str = "real",
    analysis_id: str | None = "analysis-42",
    study_context_id: str | None = "study-7",
    evidence_ids: tuple[str, ...] = ("evidence-3",),
    evidence_issue_id: str | None = "issue-2",
    claim_boundary_ids: tuple[str, ...] = ("claim.interpretation_bounded",),
    analysis_phase: str = "fitted",
) -> HelpContextBinding:
    return HelpContextBinding(
        data_context=data_context,
        analysis_id=analysis_id,
        study_context_id=study_context_id,
        evidence_ids=evidence_ids,
        evidence_issue_id=evidence_issue_id,
        claim_boundary_ids=claim_boundary_ids,
        analysis_phase=analysis_phase,
    )


def _validated_route(**kwargs) -> HelpRouteState:
    """Exercise constructor invariants through the module-private factory."""

    return navigation._new_route_state(**kwargs)


def test_open_and_return_preserve_exact_target_and_focus():
    registry = _registry()
    context = _context()

    opened = open_help_link(
        HelpRouteState.closed(),
        registry,
        "link.results.fit_chart",
        context_binding=context,
        current_context=context,
    )

    assert opened.is_open
    assert opened.help_topic_id == "help.results.fit"
    assert opened.section_id == "limits"
    assert opened.context_status is HelpContextStatus.CURRENT
    assert opened.return_destination == (
        "screen.results.fit",
        "focus.results.fit_chart",
    )

    returned = return_from_help(opened)
    assert not returned.is_open
    assert returned.help_topic_id is None
    assert returned.context_binding is None
    assert returned.return_destination == opened.return_destination


def test_double_open_is_idempotent_and_does_not_grow_history():
    registry = _registry()
    context = _context()
    first = open_help_link(
        HelpRouteState.closed(),
        registry,
        "link.results.fit_chart",
        context_binding=context,
        current_context=context,
    )
    second = open_help_link(
        first,
        registry,
        "link.results.fit_chart",
        context_binding=context,
        current_context=context,
    )

    assert second is first
    assert second.history == ()


def test_locale_change_retains_stable_route_identity_exactly():
    registry = _registry()
    opened = open_help_link(
        HelpRouteState.closed(), registry, "link.report.claims"
    )

    assert retain_route_for_locale_change(opened, locale="ja") is opened
    assert retain_route_for_locale_change(opened, locale="en") is opened


def test_topic_change_preserves_origin_and_return_but_defaults_to_static_content():
    registry = _registry()
    context = _context()
    opened = open_help_link(
        HelpRouteState.closed(),
        registry,
        "link.results.fit_chart",
        context_binding=context,
        current_context=context,
    )

    changed = change_help_topic(opened, registry, "help.report.claim_limits")

    assert changed.help_topic_id == "help.report.claim_limits"
    assert changed.section_id == "overview"
    assert changed.source_target_id == opened.source_target_id
    assert changed.return_destination == opened.return_destination
    assert changed.context_status is HelpContextStatus.STATIC_ONLY
    assert changed.context_binding is None
    assert changed.history[-1].help_topic_id == "help.results.fit"

    closed = return_from_help(changed)
    reopened = open_help_link(closed, registry, "link.report.claims")
    assert closed.history == ()
    assert reopened.history == ()


def test_delayed_topic_change_cannot_reopen_a_closed_help_session():
    registry = _registry()
    opened = open_help_link(
        HelpRouteState.closed(), registry, "link.results.fit_chart"
    )
    closed = return_from_help(opened)

    delayed = change_help_topic(closed, registry, "help.report.claim_limits")

    assert delayed is closed
    assert not delayed.is_open
    assert delayed.return_destination == opened.return_destination


def test_new_fallback_open_does_not_inherit_a_closed_sessions_return_target():
    registry = _registry()
    opened = open_help_link(
        HelpRouteState.closed(), registry, "link.results.fit_chart"
    )
    closed = return_from_help(opened)

    fallback = open_help_link(closed, registry, "link.unknown")

    assert fallback.is_fallback
    assert fallback.source_target_id is None
    assert fallback.return_destination is None


def test_reserved_fallback_topic_cannot_be_opened_as_ordinary_help():
    registry = _registry()
    opened = open_help_link(
        HelpRouteState.closed(), registry, "link.results.fit_chart"
    )

    changed = change_help_topic(opened, registry, SAFE_FALLBACK_TOPIC_ID)

    assert changed.is_fallback
    assert changed.fallback_reason is HelpFallbackReason.TOPIC_UNAVAILABLE
    with pytest.raises(ValueError, match="reserved safe destination"):
        _validated_route(
            is_open=True,
            help_topic_id=SAFE_FALLBACK_TOPIC_ID,
            section_id=SAFE_FALLBACK_SECTION_ID,
        )


def test_unknown_link_topic_or_section_routes_to_visible_safe_fallback():
    registry = _registry()
    unknown_link = open_help_link(
        HelpRouteState.closed(),
        registry,
        "link.unknown",
        fallback_return_target_id="screen.results.fit",
        fallback_return_focus_id="focus.results.fit_chart",
    )
    assert unknown_link.is_open
    assert unknown_link.is_fallback
    assert unknown_link.help_topic_id == SAFE_FALLBACK_TOPIC_ID
    assert unknown_link.fallback_reason is HelpFallbackReason.LINK_UNAVAILABLE
    assert unknown_link.return_destination == (
        "screen.results.fit",
        "focus.results.fit_chart",
    )

    unknown_topic = change_help_topic(
        unknown_link, registry, "help.unknown.topic"
    )
    assert unknown_topic.fallback_reason is HelpFallbackReason.TOPIC_UNAVAILABLE

    opened = open_help_link(
        HelpRouteState.closed(), registry, "link.report.claims"
    )
    unknown_section = change_help_topic(
        opened,
        registry,
        "help.report.claim_limits",
        section_id="missing",
    )
    assert unknown_section.fallback_reason is HelpFallbackReason.SECTION_UNAVAILABLE


def test_fallback_fails_closed_when_registry_omits_the_reserved_destination():
    registry = _registry()
    registry_without_fallback = HelpRegistry(
        topics=tuple(
            topic
            for topic in registry.topics
            if topic.help_topic_id != SAFE_FALLBACK_TOPIC_ID
        ),
        links=registry.links,
        targets=registry.targets,
        system_topic_ids=(),
    )

    with pytest.raises(ValueError, match="safe fallback route"):
        open_help_link(
            HelpRouteState.closed(), registry_without_fallback, "link.unknown"
        )

    registry_with_public_fallback = HelpRegistry(
        topics=registry.topics,
        links=registry.links,
        targets=registry.targets,
        entry_topic_ids=(SAFE_FALLBACK_TOPIC_ID,),
    )
    with pytest.raises(ValueError, match="system topic"):
        open_help_link(
            HelpRouteState.closed(), registry_with_public_fallback, "link.unknown"
        )


def test_context_policy_distinguishes_static_current_stale_and_missing():
    bound = _context()
    assert (
        evaluate_context_binding("static_only", bound, bound)
        is HelpContextStatus.STATIC_ONLY
    )
    assert (
        evaluate_context_binding("optional_current", None, None)
        is HelpContextStatus.STATIC_ONLY
    )
    assert (
        evaluate_context_binding("required_current", None, bound)
        is HelpContextStatus.MISSING_REQUIRED
    )
    assert (
        evaluate_context_binding("required_current", _context(analysis_id=None), bound)
        is HelpContextStatus.MISSING_REQUIRED
    )
    assert (
        evaluate_context_binding("required_current", bound, None)
        is HelpContextStatus.MISSING_REQUIRED
    )
    assert (
        evaluate_context_binding("required_current", bound, bound)
        is HelpContextStatus.CURRENT
    )
    assert (
        evaluate_context_binding(
            "required_current", bound, _context(analysis_id="analysis-99")
        )
        is HelpContextStatus.STALE
    )


def test_sample_and_real_contexts_never_bind_to_each_other():
    bound_sample = _context(data_context="sample")
    current_real = _context(data_context="real")

    assert (
        evaluate_context_binding("required_current", bound_sample, current_real)
        is HelpContextStatus.STALE
    )


@pytest.mark.parametrize(
    "current_context",
    (
        _context(study_context_id="study-8"),
        _context(evidence_issue_id="issue-9"),
        _context(analysis_phase="draft"),
        _context(evidence_ids=("evidence-4",)),
        _context(claim_boundary_ids=("claim.other_boundary",)),
    ),
)
def test_every_dynamic_identity_dimension_fails_closed_on_mismatch(current_context):
    assert (
        evaluate_context_binding("required_current", _context(), current_context)
        is HelpContextStatus.STALE
    )


def test_current_context_may_be_a_verified_superset_of_bound_evidence():
    current = _context(
        evidence_ids=("evidence-3", "evidence-4"),
        claim_boundary_ids=(
            "claim.interpretation_bounded",
            "claim.additional_boundary",
        ),
    )

    assert (
        evaluate_context_binding("required_current", _context(), current)
        is HelpContextStatus.CURRENT
    )


def test_invalidation_suppresses_dynamic_content_without_closing_help():
    registry = _registry()
    bound = _context()
    opened = open_help_link(
        HelpRouteState.closed(),
        registry,
        "link.results.fit_chart",
        context_binding=bound,
        current_context=bound,
    )

    invalidated = invalidate_help_route(
        opened, registry, current_context=_context(data_context="sample")
    )

    assert invalidated.is_open
    assert invalidated.help_topic_id == opened.help_topic_id
    assert invalidated.context_status is HelpContextStatus.STALE
    assert invalidated.return_destination == opened.return_destination


def test_invalidation_rejects_reused_link_id_with_changed_semantics():
    registry = _registry()
    context = _context()
    opened = open_help_link(
        HelpRouteState.closed(),
        registry,
        "link.results.fit_chart",
        context_binding=context,
        current_context=context,
    )
    changed_link = replace(
        registry.link("link.results.fit_chart"),
        help_topic_id="help.report.claim_limits",
        section_id="overview",
    )
    changed_registry = HelpRegistry(
        topics=registry.topics,
        links=(changed_link, registry.link("link.report.claims")),
        targets=registry.targets,
        entry_topic_ids=("help.results.fit",),
        system_topic_ids=registry.system_topic_ids,
    )

    invalidated = invalidate_help_route(
        opened, changed_registry, current_context=context
    )

    assert invalidated.is_fallback
    assert invalidated.fallback_reason is HelpFallbackReason.LINK_UNAVAILABLE
    assert invalidated.help_topic_id == SAFE_FALLBACK_TOPIC_ID
    assert invalidated.section_id == SAFE_FALLBACK_SECTION_ID
    assert invalidated.return_destination == opened.return_destination


def test_route_contract_has_no_scientific_or_learning_state_channel():
    route_fields = {item.name for item in fields(HelpRouteState)}
    forbidden = {
        "computation_state",
        "stability_state",
        "decision_state",
        "learning_complete",
        "diagnostic_result",
        "estimate",
    }
    assert route_fields.isdisjoint(forbidden)
    assert "_construction_permit" not in route_fields

    reducer_parameters = set(inspect.signature(open_help_link).parameters)
    assert reducer_parameters.isdisjoint(forbidden)


@pytest.mark.parametrize(
    "overrides",
    (
        {
            "context_policy": "static_only",
            "context_binding": _context(),
            "context_status": "CURRENT",
        },
        {
            "context_policy": "static_only",
            "context_status": "MISSING_REQUIRED",
        },
        {
            "context_policy": "required_current",
            "context_status": "STATIC_ONLY",
        },
        {
            "context_policy": "optional_current",
            "context_binding": _context(),
            "context_status": "STATIC_ONLY",
        },
    ),
)
def test_route_state_rejects_forged_context_policy_status_combinations(overrides):
    with pytest.raises(ValueError, match="context|STATIC_ONLY"):
        _validated_route(
            is_open=True,
            help_topic_id="help.results.fit",
            section_id="overview",
            **overrides,
        )


def test_route_state_rejects_dynamic_context_on_fallback_and_closed_routes():
    context = _context()
    with pytest.raises(ValueError, match="Fallback routes|neutral static"):
        _validated_route(
            is_open=True,
            help_topic_id=SAFE_FALLBACK_TOPIC_ID,
            section_id=SAFE_FALLBACK_SECTION_ID,
            context_policy="optional_current",
            context_binding=context,
            context_status="CURRENT",
            is_fallback=True,
            fallback_reason="help.topic_unavailable",
        )

    with pytest.raises(ValueError, match="neutral static"):
        _validated_route(
            context_policy="optional_current",
            context_status="STATIC_ONLY",
        )


def test_route_state_public_constructor_is_sealed_and_closed_history_is_empty():
    with pytest.raises(ValueError, match="navigation reducers"):
        HelpRouteState()
    with pytest.raises(ValueError, match="navigation reducers"):
        replace(
            HelpRouteState.closed(),
            context_policy="optional_current",
            context_binding=_context(),
            context_status="CURRENT",
        )
    stale = open_help_link(
        HelpRouteState.closed(),
        _registry(),
        "link.results.fit_chart",
        context_binding=_context(),
        current_context=_context(analysis_id="analysis-99"),
    )
    assert stale.context_status is HelpContextStatus.STALE
    with pytest.raises(ValueError, match="navigation reducers"):
        replace(stale, context_status="CURRENT")
    with pytest.raises(ValueError, match="closed Help route"):
        _validated_route(
            history=(
                navigation.HelpRouteLocation(
                    "help.results.fit",
                    "overview",
                ),
            )
        )
