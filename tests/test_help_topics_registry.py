from __future__ import annotations

from dataclasses import fields, is_dataclass
from enum import Enum
import re
from typing import Mapping

from mfrm_app.help_contract import HelpContextBinding, HelpContextPolicy
from mfrm_app.help_navigation import (
    SAFE_FALLBACK_SECTION_ID,
    SAFE_FALLBACK_TOPIC_ID,
    HelpContextStatus,
    HelpRouteState,
    open_help_link,
    return_from_help,
)
from mfrm_app.help_topics import (
    CANONICAL_HELP_TOPIC_IDS,
    DECLARED_HELP_CLAIM_BOUNDARY_IDS,
    DECLARED_HELP_REFERENCE_IDS,
    FALLBACK_HELP_SECTION_ID,
    FALLBACK_HELP_TOPIC_ID,
    HELP_ENTRY_TOPIC_IDS,
    HELP_ALL_LOCALE_KEYS,
    HELP_LINKS,
    HELP_NAV_LOCALE_KEYS,
    HELP_POPOVER_LINK_IDS,
    HELP_REGISTRY,
    HELP_REQUIRED_LOCALE_KEYS,
    HELP_SURFACE_CATALOG_REVIEW_STATE,
    HELP_SURFACE_TARGET_CATALOG,
    HELP_SYSTEM_TOPIC_IDS,
    HELP_TARGETS,
    HELP_TOPICS,
    HELP_SUPPORT_CATALOG_REVIEW_STATE,
    required_help_locale_keys,
    validate_builtin_help_registry,
)
from mfrm_app.terminology import TERMINOLOGY_REGISTRY
from mfrm_app.user_problems import USER_PROBLEM_SPECS


EXPECTED_CANONICAL_TOPIC_IDS = {
    "help.get_started.overview",
    "help.data.long_format",
    "help.data.mapping",
    "help.data.missingness",
    "help.design.coverage",
    "help.run.large_analysis",
    "help.run.estimation_failed",
    "help.run.nonconvergence",
    "help.problem.unexpected",
    "help.results.first_read",
    "help.results.measures_targeting",
    "help.results.fit",
    "help.results.categories",
    "help.results.residual_structure",
    "help.results.rater_evidence",
    "help.results.differential_interaction",
    "help.report.claim_limits",
    "help.downloads.privacy",
    "help.downloads.repeat_analysis",
    "help.methods.rsm_pcm_gpcm",
    "help.methods.jmle_mml",
    "help.glossary",
}


def test_registry_uses_the_complete_task_centered_topic_catalog():
    assert set(CANONICAL_HELP_TOPIC_IDS) == EXPECTED_CANONICAL_TOPIC_IDS
    assert len(CANONICAL_HELP_TOPIC_IDS) == 22
    assert len(HELP_TOPICS) == 23
    assert len(HELP_LINKS) == 74
    assert EXPECTED_CANONICAL_TOPIC_IDS.issubset(HELP_REGISTRY.topic_map)
    assert FALLBACK_HELP_TOPIC_ID in HELP_REGISTRY.topic_map
    assert set(HELP_ENTRY_TOPIC_IDS) == EXPECTED_CANONICAL_TOPIC_IDS
    assert set(HELP_SYSTEM_TOPIC_IDS) == {FALLBACK_HELP_TOPIC_ID}
    assert set(HELP_ENTRY_TOPIC_IDS).isdisjoint(HELP_SYSTEM_TOPIC_IDS)
    assert all(topic.claim_boundary_ids for topic in HELP_TOPICS)
    assert all(topic.reference_ids for topic in HELP_TOPICS)
    assert all(topic.lifecycle_states == ("all",) for topic in HELP_TOPICS)
    assert all(
        topic.version_and_review.get("review_state") == "in_review"
        for topic in HELP_TOPICS
    )
    assert all(
        topic.version_and_review.get("reviewer_roles")
        == ("sla_researcher", "ux_reviewer")
        for topic in HELP_TOPICS
    )
    assert {
        boundary_id
        for topic in HELP_TOPICS
        for boundary_id in topic.claim_boundary_ids
    } == DECLARED_HELP_CLAIM_BOUNDARY_IDS
    assert {
        reference_id
        for topic in HELP_TOPICS
        for reference_id in topic.reference_ids
    } == DECLARED_HELP_REFERENCE_IDS
    assert HELP_SUPPORT_CATALOG_REVIEW_STATE == "provisional"
    assert HELP_SURFACE_CATALOG_REVIEW_STATE == "provisional"


def test_required_locale_key_catalog_is_complete_and_value_independent():
    assert HELP_REQUIRED_LOCALE_KEYS == required_help_locale_keys(HELP_TOPICS)
    assert HELP_ALL_LOCALE_KEYS == HELP_REQUIRED_LOCALE_KEYS | HELP_NAV_LOCALE_KEYS
    assert len(
        {
            key
            for key in HELP_REQUIRED_LOCALE_KEYS
            if key.startswith("help_topics.")
        }
    ) == 437
    assert len(HELP_NAV_LOCALE_KEYS) == 24
    for topic in HELP_TOPICS:
        assert topic.title_key in HELP_REQUIRED_LOCALE_KEYS
        assert topic.summary_key in HELP_REQUIRED_LOCALE_KEYS
        assert set(topic.search_alias_keys).issubset(HELP_REQUIRED_LOCALE_KEYS)
        assert topic.safe_report_guard_key == "help_nav.safe_report_guard"
        assert topic.safe_report_guard_key in HELP_REQUIRED_LOCALE_KEYS
        for section in topic.sections:
            assert section.title_key in HELP_REQUIRED_LOCALE_KEYS
            assert set(section.content_keys).issubset(HELP_REQUIRED_LOCALE_KEYS)


def test_every_user_problem_topic_and_action_target_is_registered():
    topic_ids = set(HELP_REGISTRY.topic_map)
    target_ids = set(HELP_REGISTRY.target_map)

    assert {
        spec.help_topic_id for spec in USER_PROBLEM_SPECS.values()
    }.issubset(topic_ids)
    assert {
        target_id
        for spec in USER_PROBLEM_SPECS.values()
        for target_id in spec.action_target_ids
    }.issubset(target_ids)

    linked_problem_targets = {
        link.source_target_id
        for link in HELP_LINKS
        if link.help_link_id.startswith("link.problem.")
    }
    assert {
        target_id
        for spec in USER_PROBLEM_SPECS.values()
        for target_id in spec.action_target_ids
    }.issubset(linked_problem_targets)


def test_surface_catalog_is_independent_and_exactly_consumed():
    assert len(HELP_SURFACE_TARGET_CATALOG) == 42
    assert set(HELP_SURFACE_TARGET_CATALOG) == set(HELP_REGISTRY.target_map)
    assert {
        target.target_id: (
            target.surface_id,
            target.focus_ids[0],
            target.presentation_state["panel_id"],
        )
        for target in HELP_TARGETS
    } == dict(HELP_SURFACE_TARGET_CATALOG)
    assert len(
        {
            focus_id
            for _, focus_id, _ in HELP_SURFACE_TARGET_CATALOG.values()
        }
    ) == len(HELP_SURFACE_TARGET_CATALOG)


def test_every_terminology_help_reference_resolves_exactly():
    referenced_ids = {
        topic_id
        for term in TERMINOLOGY_REGISTRY.terms
        for topic_id in term.related_help_topic_ids
    }

    assert referenced_ids.issubset(HELP_REGISTRY.topic_map)
    assert referenced_ids.issubset(EXPECTED_CANONICAL_TOPIC_IDS)


def test_existing_popover_keys_resolve_only_through_exact_registered_links():
    assert len(HELP_POPOVER_LINK_IDS) == 18
    assert set(HELP_POPOVER_LINK_IDS) == {
        "scree",
        "fit_scatter",
        "pathway_map",
        "misfit_ranking",
        "zstd_distribution",
        "qq_residuals",
        "category_probability",
        "category_usage",
        "threshold_map",
        "coverage_heatmap",
        "wright_map",
        "ecdf_measures",
        "forest_measures",
        "facet_distribution",
        "bias_heatmap",
        "classical_dif",
        "rater_agreement",
        "mml_person_sd",
    }
    for popover_key, link_id in HELP_POPOVER_LINK_IDS.items():
        assert HELP_REGISTRY.link(link_id) is not None
        assert HELP_REGISTRY.link(link_id.replace(popover_key, popover_key.upper())) is None


def test_fallback_topic_and_section_match_the_navigation_reservation():
    assert FALLBACK_HELP_TOPIC_ID == SAFE_FALLBACK_TOPIC_ID
    assert FALLBACK_HELP_SECTION_ID == SAFE_FALLBACK_SECTION_ID
    fallback = HELP_REGISTRY.topic(FALLBACK_HELP_TOPIC_ID)
    assert fallback is not None
    assert fallback.has_section(FALLBACK_HELP_SECTION_ID)
    assert "target.help.home" in fallback.related_target_ids

    routed = open_help_link(
        HelpRouteState.closed(), HELP_REGISTRY, "link.not_registered"
    )
    assert routed.is_open
    assert routed.is_fallback
    assert routed.help_topic_id == FALLBACK_HELP_TOPIC_ID
    assert routed.section_id == FALLBACK_HELP_SECTION_ID


def test_registry_is_zero_orphan_and_all_return_focus_pairs_are_exact():
    validate_builtin_help_registry()
    assert HELP_REGISTRY.require_zero_orphans
    assert len(HELP_REGISTRY.topic_map) == len(HELP_TOPICS)
    assert len(HELP_REGISTRY.target_map) == len(HELP_TARGETS)
    assert len(HELP_REGISTRY.link_map) == len(HELP_LINKS)

    for link in HELP_LINKS:
        topic = HELP_REGISTRY.topic(link.help_topic_id)
        target = HELP_REGISTRY.target(link.return_target_id)
        assert topic is not None
        assert link.section_id is None or topic.has_section(link.section_id)
        assert target is not None
        assert link.return_focus_id in target.focus_ids


def test_context_policies_cover_static_optional_and_verified_current_routes():
    policies = {link.context_policy for link in HELP_LINKS}
    assert policies == {
        HelpContextPolicy.STATIC_ONLY,
        HelpContextPolicy.OPTIONAL_CURRENT,
        HelpContextPolicy.REQUIRED_CURRENT,
    }
    assert all(
        HELP_REGISTRY.link(link_id).context_policy
        is HelpContextPolicy.REQUIRED_CURRENT
        for link_id in HELP_POPOVER_LINK_IDS.values()
    )


def test_open_and_return_smoke_uses_stable_ids_without_locale_values():
    context = HelpContextBinding(
        data_context="real",
        study_context_id="study-7",
        analysis_id="analysis-42",
        evidence_ids=("evidence-3",),
        claim_boundary_ids=("claim.fit_not_automatic_exclusion",),
        analysis_phase="fitted",
    )
    opened = open_help_link(
        HelpRouteState.closed(),
        HELP_REGISTRY,
        HELP_POPOVER_LINK_IDS["fit_scatter"],
        context_binding=context,
        current_context=context,
    )

    assert opened.help_topic_id == "help.results.fit"
    assert opened.section_id == "fit-scatter"
    assert opened.context_status is HelpContextStatus.CURRENT
    assert opened.return_destination == (
        "target.results.figure.fit_scatter",
        "focus.results.figure.fit_scatter",
    )

    returned = return_from_help(opened)
    assert not returned.is_open
    assert returned.return_destination == opened.return_destination


def _metadata_strings(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    if isinstance(value, Enum):
        return (str(value.value),)
    if is_dataclass(value) and not isinstance(value, type):
        return tuple(
            item
            for field_info in fields(value)
            for item in _metadata_strings(getattr(value, field_info.name))
        )
    if isinstance(value, Mapping):
        return tuple(
            item
            for key, nested in value.items()
            for item in (*_metadata_strings(key), *_metadata_strings(nested))
        )
    if isinstance(value, (tuple, list, set, frozenset)):
        return tuple(item for nested in value for item in _metadata_strings(nested))
    return (str(value),)


def test_active_help_metadata_contains_no_external_execution_call_to_action():
    metadata = "\n".join(
        _metadata_strings((HELP_TOPICS, HELP_TARGETS, HELP_LINKS))
    ).casefold()
    tokens = set(re.findall(r"[a-z0-9_]+", metadata))
    prohibited_tokens = {
        "tam",
        "conquest",
        "mirt",
        "rscript",
        "julia",
        "cmdstanpy",
        "posterior_viewer",
    }
    prohibited_phrases = {
        "cross-engine",
        "cross engine",
        "external runtime",
        "external engine",
    }

    assert tokens.isdisjoint(prohibited_tokens)
    assert all(phrase not in metadata for phrase in prohibited_phrases)
