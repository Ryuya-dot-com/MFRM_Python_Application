from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from mfrm_app.help_contract import (
    HelpContextBinding,
    HelpContractValidationError,
    HelpLink,
    HelpRegistry,
    HelpSection,
    HelpTarget,
    HelpTopic,
    validate_help_registry,
)


def _review() -> dict[str, object]:
    return {
        "content_version": "1.0.0",
        "owner_role": "help_content_owner",
        "review_state": "in_review",
        "last_reviewed": "2026-07-24",
        "reviewer_roles": ("sla_researcher", "ux_reviewer"),
    }


def _topic(
    help_topic_id: str = "help.results.fit",
    *,
    related_target_ids: tuple[str, ...] = ("screen.results.fit",),
) -> HelpTopic:
    key_root = help_topic_id.replace("help.", "help_topics.")
    return HelpTopic(
        help_topic_id=help_topic_id,
        title_key=f"{key_root}.title",
        summary_key=f"{key_root}.summary",
        sections=(
            HelpSection(
                section_id="overview",
                title_key=f"{key_root}.sections.overview_title",
                content_keys=(f"{key_root}.sections.overview_body",),
            ),
            HelpSection(
                section_id="limits",
                title_key=f"{key_root}.sections.limits_title",
                content_keys=(f"{key_root}.sections.limits_body",),
            ),
        ),
        version_and_review=_review(),
        concept_ids=("concept.model_fit",),
        lifecycle_states=("no_data", "fitted"),
        applicability={"models": ["rsm", "pcm"]},
        prerequisite_keys=(f"{key_root}.prerequisite",),
        computed_key=f"{key_root}.computed",
        can_show_key=f"{key_root}.can_show",
        cannot_show_key=f"{key_root}.cannot_show",
        next_check_key=f"{key_root}.next_check",
        next_action_key=f"{key_root}.next_action",
        safe_report_guard_key="help_nav.safe_report_guard",
        safe_report_key=f"{key_root}.safe_report",
        avoid_report_key=f"{key_root}.avoid_report",
        claim_boundary_ids=("claim.fit_not_validity",),
        reference_ids=("reference.wright1994",),
        related_target_ids=related_target_ids,
        search_alias_keys=(f"{key_root}.aliases.fit",),
    )


def _target(
    target_id: str = "screen.results.fit",
    *,
    focus_ids: tuple[str, ...] = ("focus.results.fit_chart",),
) -> HelpTarget:
    return HelpTarget(
        target_id=target_id,
        surface_id="surface.results",
        focus_ids=focus_ids,
        presentation_state={"panel": "fit", "expanded": True},
    )


def _link(
    *,
    help_topic_id: str = "help.results.fit",
    section_id: str = "limits",
    source_target_id: str = "screen.results.fit",
    return_target_id: str = "screen.results.fit",
    return_focus_id: str = "focus.results.fit_chart",
) -> HelpLink:
    return HelpLink(
        help_link_id="link.results.fit_chart",
        source_target_id=source_target_id,
        help_topic_id=help_topic_id,
        section_id=section_id,
        return_target_id=return_target_id,
        return_focus_id=return_focus_id,
        context_policy="optional_current",
    )


def test_contract_schemas_are_frozen_and_nested_metadata_is_immutable():
    topic = _topic()
    target = _target()
    binding = HelpContextBinding(
        data_context="real",
        analysis_id="analysis-42",
        evidence_ids=("evidence-3",),
        claim_boundary_ids=("claim.fit_not_validity",),
        analysis_phase="fitted",
    )

    with pytest.raises(FrozenInstanceError):
        topic.help_topic_id = "help.changed"  # type: ignore[misc]
    with pytest.raises(TypeError):
        topic.applicability["models"] = ("gpcm",)  # type: ignore[index]
    with pytest.raises(TypeError):
        target.presentation_state["panel"] = "other"  # type: ignore[index]
    with pytest.raises(FrozenInstanceError):
        binding.analysis_id = "other"  # type: ignore[misc]

    for unsafe_reference in (
        "/Users/private/ratings.csv",
        "student@example.org",
        "a participant name",
    ):
        with pytest.raises(HelpContractValidationError, match="opaque"):
            HelpContextBinding(data_context="real", analysis_id=unsafe_reference)


@pytest.mark.parametrize(
    "invalid_id",
    ["Analysis Workflow", "Help.Results.Fit", "help..fit", "結果.fit", " help.fit"],
)
def test_route_contracts_reject_display_text_and_invalid_stable_ids(invalid_id: str):
    with pytest.raises(HelpContractValidationError, match="stable ID"):
        _target(invalid_id)


def test_topic_rejects_duplicate_or_unknown_section_layers():
    section = HelpSection(
        section_id="overview",
        title_key="help_topics.fit.overview_title",
        content_keys=("help_topics.fit.overview_body",),
    )
    with pytest.raises(HelpContractValidationError, match="section IDs"):
        HelpTopic(
            help_topic_id="help.results.fit",
            title_key="help_topics.fit.title",
            summary_key="help_topics.fit.summary",
            sections=(section, section),
            version_and_review=_review(),
        )
    with pytest.raises(HelpContractValidationError, match="parent topic"):
        HelpTopic(
            help_topic_id="help.results.fit",
            title_key="help_topics.fit.title",
            summary_key="help_topics.fit.summary",
            sections=(section,),
            version_and_review=_review(),
            audience_layers=("standard",),
        )
    with pytest.raises(HelpContractValidationError, match="cannot be combined"):
        HelpTopic(
            help_topic_id="help.results.fit",
            title_key="help_topics.fit.title",
            summary_key="help_topics.fit.summary",
            sections=(section,),
            version_and_review=_review(),
            lifecycle_states=("all", "fitted"),
        )


def test_topic_lifecycle_applicability_is_explicit_adapter_metadata():
    topic = _topic()

    assert topic.supports_lifecycle("no_data")
    assert topic.supports_lifecycle("fitted")
    assert not topic.supports_lifecycle("draft")


def test_topic_requires_a_guard_for_every_report_wording_example():
    topic = _topic()
    with pytest.raises(HelpContractValidationError, match="declared together"):
        HelpTopic(
            help_topic_id=topic.help_topic_id,
            title_key=topic.title_key,
            summary_key=topic.summary_key,
            sections=topic.sections,
            version_and_review=_review(),
            safe_report_key="help_topics.results.fit.safe_report",
        )


def test_topic_rejects_unaccountable_review_metadata():
    topic = _topic()
    invalid_review = {
        "content_version": "1.0.0",
        "owner_role": "help_content_owner",
        "review_state": "reviewed",
        "last_reviewed": "2026-07-24",
        "reviewer_roles": (),
    }

    with pytest.raises(HelpContractValidationError, match="reviewer role"):
        HelpTopic(
            help_topic_id=topic.help_topic_id,
            title_key=topic.title_key,
            summary_key=topic.summary_key,
            sections=topic.sections,
            version_and_review=invalid_review,
        )


def test_registry_accepts_complete_zero_orphan_graph():
    topic = _topic()
    target = _target()
    link = _link()

    registry = HelpRegistry(topics=(topic,), links=(link,), targets=(target,))

    assert registry.topic(topic.help_topic_id) is topic
    assert registry.link(link.help_link_id) is link
    assert registry.target(target.target_id) is target


def test_registry_rejects_missing_topic_target_and_related_target_references():
    topic = _topic()
    target = _target()
    with pytest.raises(HelpContractValidationError, match="unknown topic"):
        HelpRegistry(
            topics=(topic,),
            links=(_link(help_topic_id="help.missing.topic"),),
            targets=(target,),
        )

    with pytest.raises(HelpContractValidationError, match="unknown source target"):
        HelpRegistry(
            topics=(topic,),
            links=(_link(source_target_id="screen.missing"),),
            targets=(target,),
        )

    with pytest.raises(HelpContractValidationError, match="unknown targets"):
        HelpRegistry(
            topics=(_topic(related_target_ids=("screen.missing",)),),
            links=(_link(),),
            targets=(target,),
        )


def test_registry_rejects_unknown_section_and_inconsistent_return_focus():
    topic = _topic()
    target = _target()
    with pytest.raises(HelpContractValidationError, match="unknown section"):
        HelpRegistry(
            topics=(topic,),
            links=(_link(section_id="not_registered"),),
            targets=(target,),
        )

    with pytest.raises(HelpContractValidationError, match="is not owned"):
        HelpRegistry(
            topics=(topic,),
            links=(_link(return_focus_id="focus.results.other"),),
            targets=(target,),
        )


def test_registry_rejects_duplicate_ids_and_orphan_nodes():
    topic = _topic()
    target = _target()
    with pytest.raises(HelpContractValidationError, match="Duplicate help_topic_id"):
        HelpRegistry(
            topics=(topic, topic),
            links=(_link(),),
            targets=(target,),
        )

    unlinked_topic = _topic(
        "help.report.claim_limits", related_target_ids=("screen.report.claims",)
    )
    report_target = _target(
        "screen.report.claims", focus_ids=("focus.report.claims",)
    )
    with pytest.raises(HelpContractValidationError, match="orphan nodes"):
        HelpRegistry(
            topics=(topic, unlinked_topic),
            links=(_link(),),
            targets=(target, report_target),
        )

    # Help-home entry topics and their related targets are explicitly reachable.
    registry = HelpRegistry(
        topics=(topic, unlinked_topic),
        links=(_link(),),
        targets=(target, report_target),
        entry_topic_ids=("help.report.claim_limits",),
    )
    assert registry.topic("help.report.claim_limits") is unlinked_topic


def test_registry_keeps_system_topics_out_of_help_home_entries():
    topic = _topic()
    target = _target()

    registry = HelpRegistry(
        topics=(topic,),
        links=(),
        targets=(target,),
        system_topic_ids=(topic.help_topic_id,),
    )

    assert registry.entry_topic_ids == ()
    assert registry.system_topic_ids == (topic.help_topic_id,)
    with pytest.raises(HelpContractValidationError, match="both entry and system"):
        HelpRegistry(
            topics=(topic,),
            links=(),
            targets=(target,),
            entry_topic_ids=(topic.help_topic_id,),
            system_topic_ids=(topic.help_topic_id,),
        )


def test_registry_rejects_a_focus_id_owned_by_multiple_targets():
    topic = _topic()
    target = _target()
    second_target = _target(
        "screen.results.other",
        focus_ids=("focus.results.fit_chart",),
    )

    with pytest.raises(HelpContractValidationError, match="global target owner"):
        HelpRegistry(
            topics=(topic,),
            links=(_link(),),
            targets=(target, second_target),
            entry_topic_ids=(topic.help_topic_id,),
            require_zero_orphans=False,
        )


def test_registry_checks_external_concept_boundary_and_reference_catalogs():
    topic = _topic()
    target = _target()
    link = _link()

    validate_help_registry(
        (topic,),
        (link,),
        (target,),
        known_concept_ids=("concept.model_fit",),
        known_claim_boundary_ids=("claim.fit_not_validity",),
        known_reference_ids=("reference.wright1994",),
    )
    with pytest.raises(HelpContractValidationError, match="unknown ClaimBoundary"):
        validate_help_registry(
            (topic,),
            (link,),
            (target,),
            known_claim_boundary_ids=("claim.other",),
        )
