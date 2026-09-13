"""Reviewed stable IDs for the first standalone Help registry.

The registry contains metadata and locale keys, not localized prose.  It is a
pure projection over :mod:`mfrm_app.help_contract`: importing it neither loads
the application UI nor invokes an analysis.  Long task-centred topic IDs are
the canonical routes.  Existing popover keys map explicitly to registered
links; invented topic aliases and translated-title routing are not accepted.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Iterable, Mapping

from .help_contract import (
    HelpContextPolicy,
    HelpLink,
    HelpRegistry,
    HelpSection,
    HelpTarget,
    HelpTopic,
    validate_help_registry,
)
from .terminology import TERMINOLOGY_REGISTRY
from .user_problems import USER_PROBLEM_SPECS


HELP_TOPIC_REGISTRY_VERSION = "mfrm_help_topics_v1"


HELP_NAV_LOCALE_KEYS = frozenset(
    {
        "help_nav.open_help",
        "help_nav.open_full_guide",
        "help_nav.help_home",
        "help_nav.back_to_source",
        "help_nav.returned_status",
        "help_nav.topic_label",
        "help_nav.prerequisites_heading",
        "help_nav.can_show_heading",
        "help_nav.cannot_show_heading",
        "help_nav.next_check_heading",
        "help_nav.next_action_heading",
        "help_nav.details_heading",
        "help_nav.computed_heading",
        "help_nav.reporting_example_heading",
        "help_nav.avoid_report_heading",
        "help_nav.search_label",
        "help_nav.search_placeholder",
        "help_nav.static_only_notice",
        "help_nav.stale_context_notice",
        "help_nav.missing_context_notice",
        "help_nav.safe_report_guard",
        "help_nav.unavailable_title",
        "help_nav.unavailable_body",
        "help_nav.unavailable_home_action",
    }
)


CANONICAL_HELP_TOPIC_IDS = (
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
)

FALLBACK_HELP_TOPIC_ID = "help.fallback.unavailable"
FALLBACK_HELP_SECTION_ID = "overview"


_TOPIC_SECTIONS: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "help.get_started.overview": (
            "overview",
            "shortest-path",
            "claim-boundaries",
        ),
        "help.data.long_format": (
            "overview",
            "row-meaning",
            "required-columns",
            "missingness",
        ),
        "help.data.mapping": ("overview", "roles", "check-mapping"),
        "help.data.missingness": (
            "overview",
            "planned-missing",
            "unexpected-missing",
            "limits",
        ),
        "help.design.coverage": (
            "overview",
            "observed-coverage",
            "connectedness",
            "limits",
        ),
        "help.run.large_analysis": (
            "overview",
            "retained-state",
            "safe-scope",
        ),
        "help.run.estimation_failed": (
            "overview",
            "preserved-state",
            "next-check",
        ),
        "help.run.nonconvergence": (
            "overview",
            "meaning",
            "does-not-mean",
            "next-check",
        ),
        "help.problem.unexpected": (
            "overview",
            "first-check",
            "get-support",
        ),
        "help.results.first_read": (
            "overview",
            "reading-order",
            "claim-boundaries",
        ),
        "help.results.measures_targeting": (
            "overview",
            "wright-map",
            "ecdf-measures",
            "forest-measures",
            "facet-distribution",
            "limits",
        ),
        "help.results.fit": (
            "overview",
            "fit-scatter",
            "pathway-map",
            "misfit-ranking",
            "zstd-distribution",
            "qq-residuals",
            "limits",
        ),
        "help.results.categories": (
            "overview",
            "category-probability",
            "category-usage",
            "threshold-map",
            "limits",
        ),
        "help.results.residual_structure": (
            "overview",
            "residual-pca",
            "limits",
            "next-check",
        ),
        "help.results.rater_evidence": (
            "overview",
            "agreement",
            "severity",
            "differentiation",
            "limits",
        ),
        "help.results.differential_interaction": (
            "overview",
            "bias-heatmap",
            "classical-dif",
            "fairness-boundary",
        ),
        "help.report.claim_limits": (
            "overview",
            "supported",
            "limited",
            "unavailable",
            "wording",
        ),
        "help.downloads.privacy": (
            "overview",
            "included",
            "excluded",
            "sharing-check",
        ),
        "help.downloads.repeat_analysis": (
            "overview",
            "settings",
            "evidence",
            "limitations",
        ),
        "help.methods.rsm_pcm_gpcm": (
            "overview",
            "rsm",
            "pcm",
            "gpcm",
            "choice-boundary",
        ),
        "help.methods.jmle_mml": (
            "overview",
            "jmle",
            "mml",
            "person-distribution-sd",
            "choice-boundary",
        ),
        "help.glossary": ("overview", "search", "concept-boundaries"),
        FALLBACK_HELP_TOPIC_ID: (FALLBACK_HELP_SECTION_ID,),
    }
)


_PRIMARY_TARGET_BY_TOPIC: Mapping[str, str] = MappingProxyType(
    {
        "help.get_started.overview": "target.input.start",
        "help.data.long_format": "target.input.format_delimiter",
        "help.data.mapping": "target.input.column_mapping",
        "help.data.missingness": "target.input.data_audit",
        "help.design.coverage": "target.settings.design_constraints",
        "help.run.large_analysis": "target.settings.analysis_scope",
        "help.run.estimation_failed": "target.settings.estimation",
        "help.run.nonconvergence": "target.settings.estimation",
        "help.problem.unexpected": "target.help.problem_resolution",
        "help.results.first_read": "target.results.first_read",
        "help.results.measures_targeting": "target.results.measures_targeting",
        "help.results.fit": "target.results.fit",
        "help.results.categories": "target.results.categories",
        "help.results.residual_structure": "target.results.residual_structure",
        "help.results.rater_evidence": "target.results.rater_evidence",
        "help.results.differential_interaction": (
            "target.results.differential_interaction"
        ),
        "help.report.claim_limits": "target.report.claim_limits",
        "help.downloads.privacy": "target.downloads.privacy",
        "help.downloads.repeat_analysis": "target.downloads.repeat_analysis",
        "help.methods.rsm_pcm_gpcm": "target.settings.models",
        "help.methods.jmle_mml": "target.settings.estimation",
        "help.glossary": "target.help.glossary",
        FALLBACK_HELP_TOPIC_ID: "target.help.home",
    }
)


_CONCEPT_IDS_BY_TOPIC: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "help.get_started.overview": ("evidence.summary",),
        "help.results.first_read": (
            "evidence.summary",
            "evidence.availability",
            "evidence.limit_reason",
            "evidence.planned_check_stability",
        ),
        "help.results.measures_targeting": (
            "scale.logit",
            "measurement.separation",
        ),
        "help.results.fit": ("fit.infit", "fit.outfit"),
        "help.results.categories": ("model.rsm", "model.pcm", "model.gpcm"),
        "help.results.residual_structure": ("residual_structure.pca",),
        "help.results.rater_evidence": (
            "model.mfrm",
            "measurement.separation",
        ),
        "help.results.differential_interaction": ("model.mfrm",),
        "help.report.claim_limits": (
            "evidence.availability",
            "evidence.limit_reason",
            "evidence.planned_check_stability",
        ),
        "help.downloads.privacy": (
            "download.contents_record",
            "runtime.browser_retention",
        ),
        "help.downloads.repeat_analysis": (
            "workflow.analysis_reference",
            "reproducibility.input_match",
            "reproducibility.script",
        ),
        "help.methods.rsm_pcm_gpcm": (
            "model.mfrm",
            "model.rsm",
            "model.pcm",
            "model.gpcm",
        ),
        "help.methods.jmle_mml": ("estimation.jmle", "estimation.mml"),
        "help.glossary": tuple(
            term.concept_id for term in TERMINOLOGY_REGISTRY.terms
        ),
    }
)


_CLAIM_BOUNDARY_BY_TOPIC: Mapping[str, str] = MappingProxyType(
    {
        "help.data.long_format": "claim.data_shape_not_quality",
        "help.data.mapping": "claim.mapping_not_design_validity",
        "help.data.missingness": "claim.missingness_requires_context",
        "help.design.coverage": "claim.coverage_not_representativeness",
        "help.run.large_analysis": "claim.resource_limit_not_result",
        "help.run.estimation_failed": "claim.failure_not_substantive_result",
        "help.run.nonconvergence": "claim.nonconvergence_limits_interpretation",
        "help.problem.unexpected": "claim.unexpected_failure_not_result",
        "help.results.first_read": "claim.single_screen_not_sufficient",
        "help.results.measures_targeting": "claim.measure_not_person_judgment",
        "help.results.fit": "claim.fit_not_automatic_exclusion",
        "help.results.categories": "claim.category_pattern_not_quality_verdict",
        "help.results.residual_structure": "claim.residual_screen_not_dimension_proof",
        "help.results.rater_evidence": "claim.rater_indicators_remain_distinct",
        "help.results.differential_interaction": "claim.interaction_not_fairness_verdict",
        "help.report.claim_limits": "claim.reporting_follows_available_evidence",
        "help.downloads.privacy": "claim.package_check_not_privacy_guarantee",
        "help.downloads.repeat_analysis": "claim.repeatability_not_validation",
        "help.methods.rsm_pcm_gpcm": "claim.model_choice_requires_design",
        "help.methods.jmle_mml": "claim.estimator_choice_requires_purpose",
        "help.glossary": "claim.term_label_not_conclusion",
        "help.get_started.overview": "claim.workflow_completion_not_validity",
        FALLBACK_HELP_TOPIC_ID: "claim.help_unavailable_preserves_source_limit",
    }
)


_REFERENCE_ID_BY_TOPIC: Mapping[str, str] = MappingProxyType(
    {
        "help.get_started.overview": "reference.product.safe_workflow",
        "help.data.long_format": "reference.product.input_contract",
        "help.data.mapping": "reference.product.input_contract",
        "help.data.missingness": "reference.measurement.missingness",
        "help.design.coverage": "reference.measurement.connected_design",
        "help.run.large_analysis": "reference.product.resource_boundary",
        "help.run.estimation_failed": "reference.product.problem_recovery",
        "help.run.nonconvergence": "reference.measurement.estimation_checks",
        "help.problem.unexpected": "reference.product.problem_recovery",
        "help.results.first_read": "reference.measurement.evidence_sequence",
        "help.results.measures_targeting": "reference.measurement.scale_interpretation",
        "help.results.fit": "reference.measurement.fit_interpretation",
        "help.results.categories": "reference.measurement.category_interpretation",
        "help.results.residual_structure": "reference.measurement.residual_structure",
        "help.results.rater_evidence": "reference.measurement.rater_interpretation",
        "help.results.differential_interaction": (
            "reference.measurement.interaction_interpretation"
        ),
        "help.report.claim_limits": "reference.product.reporting_boundary",
        "help.downloads.privacy": "reference.product.sharing_boundary",
        "help.downloads.repeat_analysis": "reference.product.repeat_analysis",
        "help.methods.rsm_pcm_gpcm": "reference.measurement.response_models",
        "help.methods.jmle_mml": "reference.measurement.estimation_methods",
        "help.glossary": "reference.measurement.glossary",
        FALLBACK_HELP_TOPIC_ID: "reference.product.problem_recovery",
    }
)


# These declarations are independent of the topic-to-ID maps above so a typo
# cannot certify itself during graph validation.  They are provisional R1
# identifiers, not yet a bibliographic or ClaimBoundary record catalog; visible
# method Help remains gated on the review described in the roadmap.
HELP_SUPPORT_CATALOG_REVIEW_STATE = "provisional"
DECLARED_HELP_CLAIM_BOUNDARY_IDS = frozenset(
    {
        "claim.category_pattern_not_quality_verdict",
        "claim.coverage_not_representativeness",
        "claim.data_shape_not_quality",
        "claim.estimator_choice_requires_purpose",
        "claim.failure_not_substantive_result",
        "claim.fit_not_automatic_exclusion",
        "claim.help_unavailable_preserves_source_limit",
        "claim.interaction_not_fairness_verdict",
        "claim.mapping_not_design_validity",
        "claim.measure_not_person_judgment",
        "claim.missingness_requires_context",
        "claim.model_choice_requires_design",
        "claim.nonconvergence_limits_interpretation",
        "claim.package_check_not_privacy_guarantee",
        "claim.rater_indicators_remain_distinct",
        "claim.repeatability_not_validation",
        "claim.reporting_follows_available_evidence",
        "claim.residual_screen_not_dimension_proof",
        "claim.resource_limit_not_result",
        "claim.single_screen_not_sufficient",
        "claim.term_label_not_conclusion",
        "claim.unexpected_failure_not_result",
        "claim.workflow_completion_not_validity",
    }
)
DECLARED_HELP_REFERENCE_IDS = frozenset(
    {
        "reference.measurement.category_interpretation",
        "reference.measurement.connected_design",
        "reference.measurement.estimation_checks",
        "reference.measurement.estimation_methods",
        "reference.measurement.evidence_sequence",
        "reference.measurement.fit_interpretation",
        "reference.measurement.glossary",
        "reference.measurement.interaction_interpretation",
        "reference.measurement.missingness",
        "reference.measurement.rater_interpretation",
        "reference.measurement.residual_structure",
        "reference.measurement.response_models",
        "reference.measurement.scale_interpretation",
        "reference.product.input_contract",
        "reference.product.problem_recovery",
        "reference.product.repeat_analysis",
        "reference.product.reporting_boundary",
        "reference.product.resource_boundary",
        "reference.product.safe_workflow",
        "reference.product.sharing_boundary",
    }
)


_POPOVER_DESTINATIONS = (
    ("scree", "help.results.residual_structure", "residual-pca"),
    ("fit_scatter", "help.results.fit", "fit-scatter"),
    ("pathway_map", "help.results.fit", "pathway-map"),
    ("misfit_ranking", "help.results.fit", "misfit-ranking"),
    ("zstd_distribution", "help.results.fit", "zstd-distribution"),
    ("qq_residuals", "help.results.fit", "qq-residuals"),
    ("category_probability", "help.results.categories", "category-probability"),
    ("category_usage", "help.results.categories", "category-usage"),
    ("threshold_map", "help.results.categories", "threshold-map"),
    ("coverage_heatmap", "help.design.coverage", "observed-coverage"),
    ("wright_map", "help.results.measures_targeting", "wright-map"),
    ("ecdf_measures", "help.results.measures_targeting", "ecdf-measures"),
    ("forest_measures", "help.results.measures_targeting", "forest-measures"),
    ("facet_distribution", "help.results.measures_targeting", "facet-distribution"),
    (
        "bias_heatmap",
        "help.results.differential_interaction",
        "bias-heatmap",
    ),
    (
        "classical_dif",
        "help.results.differential_interaction",
        "classical-dif",
    ),
    ("rater_agreement", "help.results.rater_evidence", "agreement"),
    ("mml_person_sd", "help.methods.jmle_mml", "person-distribution-sd"),
)


def _ordered_unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(values))


def _locale_root(help_topic_id: str) -> str:
    return "help_topics." + help_topic_id.removeprefix("help.")


def _locale_section(section_id: str) -> str:
    return section_id.replace("-", "_")


def _sections(help_topic_id: str) -> tuple[HelpSection, ...]:
    root = _locale_root(help_topic_id)
    return tuple(
        HelpSection(
            section_id=section_id,
            title_key=f"{root}.sections.{_locale_section(section_id)}.title",
            content_keys=(
                f"{root}.sections.{_locale_section(section_id)}.body",
            ),
        )
        for section_id in _TOPIC_SECTIONS[help_topic_id]
    )


def _figure_target_id(popover_key: str) -> str:
    return f"target.results.figure.{popover_key}"


def _related_targets(help_topic_id: str) -> tuple[str, ...]:
    figure_targets = (
        _figure_target_id(popover_key)
        for popover_key, topic_id, _ in _POPOVER_DESTINATIONS
        if topic_id == help_topic_id
    )
    return _ordered_unique(
        (_PRIMARY_TARGET_BY_TOPIC[help_topic_id], *figure_targets)
    )


def _topic(help_topic_id: str) -> HelpTopic:
    root = _locale_root(help_topic_id)
    return HelpTopic(
        help_topic_id=help_topic_id,
        title_key=f"{root}.title",
        summary_key=f"{root}.summary",
        sections=_sections(help_topic_id),
        concept_ids=_CONCEPT_IDS_BY_TOPIC.get(help_topic_id, ()),
        lifecycle_states=("all",),
        applicability={"product_scope": "python", "content_kind": "help"},
        prerequisite_keys=(f"{root}.prerequisites",),
        computed_key=f"{root}.computed",
        can_show_key=f"{root}.can_show",
        cannot_show_key=f"{root}.cannot_show",
        next_check_key=f"{root}.next_check",
        next_action_key=f"{root}.next_action",
        safe_report_guard_key="help_nav.safe_report_guard",
        safe_report_key=f"{root}.safe_report",
        avoid_report_key=f"{root}.avoid_report",
        claim_boundary_ids=(_CLAIM_BOUNDARY_BY_TOPIC[help_topic_id],),
        reference_ids=(_REFERENCE_ID_BY_TOPIC[help_topic_id],),
        related_target_ids=_related_targets(help_topic_id),
        search_alias_keys=(f"{root}.search_alias",),
        version_and_review={
            "content_version": "1.0.0",
            "owner_role": "help_content_owner",
            "review_state": "in_review",
            "last_reviewed": "2026-07-24",
            "reviewer_roles": ("sla_researcher", "ux_reviewer"),
        },
    )


_CANONICAL_TOPICS = tuple(
    _topic(help_topic_id)
    for help_topic_id in (*CANONICAL_HELP_TOPIC_IDS, FALLBACK_HELP_TOPIC_ID)
)
HELP_TOPICS = _CANONICAL_TOPICS


# This catalog is deliberately independent of topics, links, popovers, and
# problem actions.  A typo in one of those consumers therefore cannot create a
# target that certifies itself.  R2 must bind each provisional panel/focus ID
# to the real Streamlit surface and verify the adapter before visible rollout.
HELP_SURFACE_CATALOG_REVIEW_STATE = "provisional"
HELP_SURFACE_TARGET_CATALOG: Mapping[str, tuple[str, str, str]] = MappingProxyType(
    {
        "target.downloads.privacy": (
            "surface.downloads",
            "focus.downloads.privacy",
            "panel.downloads.privacy",
        ),
        "target.downloads.repeat_analysis": (
            "surface.downloads",
            "focus.downloads.repeat_analysis",
            "panel.downloads.repeat_analysis",
        ),
        "target.help.glossary": (
            "surface.help",
            "focus.help.glossary",
            "panel.help.glossary",
        ),
        "target.help.home": (
            "surface.help",
            "focus.help.home",
            "panel.help.home",
        ),
        "target.help.problem_resolution": (
            "surface.help",
            "focus.help.problem_resolution",
            "panel.help.problem_resolution",
        ),
        "target.input.column_mapping": (
            "surface.input",
            "focus.input.column_mapping",
            "panel.input.column_mapping",
        ),
        "target.input.data_audit": (
            "surface.input",
            "focus.input.data_audit",
            "panel.input.data_audit",
        ),
        "target.input.format_delimiter": (
            "surface.input",
            "focus.input.format_delimiter",
            "panel.input.format_delimiter",
        ),
        "target.input.score_categories": (
            "surface.input",
            "focus.input.score_categories",
            "panel.input.score_categories",
        ),
        "target.input.start": (
            "surface.input",
            "focus.input.start",
            "panel.input.start",
        ),
        "target.report.claim_limits": (
            "surface.report",
            "focus.report.claim_limits",
            "panel.report.claim_limits",
        ),
        "target.results.categories": (
            "surface.results",
            "focus.results.categories",
            "panel.results.categories",
        ),
        "target.results.differential_interaction": (
            "surface.results",
            "focus.results.differential_interaction",
            "panel.results.differential_interaction",
        ),
        "target.results.figure.bias_heatmap": (
            "surface.results",
            "focus.results.figure.bias_heatmap",
            "panel.results.figure.bias_heatmap",
        ),
        "target.results.figure.category_probability": (
            "surface.results",
            "focus.results.figure.category_probability",
            "panel.results.figure.category_probability",
        ),
        "target.results.figure.category_usage": (
            "surface.results",
            "focus.results.figure.category_usage",
            "panel.results.figure.category_usage",
        ),
        "target.results.figure.classical_dif": (
            "surface.results",
            "focus.results.figure.classical_dif",
            "panel.results.figure.classical_dif",
        ),
        "target.results.figure.coverage_heatmap": (
            "surface.results",
            "focus.results.figure.coverage_heatmap",
            "panel.results.figure.coverage_heatmap",
        ),
        "target.results.figure.ecdf_measures": (
            "surface.results",
            "focus.results.figure.ecdf_measures",
            "panel.results.figure.ecdf_measures",
        ),
        "target.results.figure.facet_distribution": (
            "surface.results",
            "focus.results.figure.facet_distribution",
            "panel.results.figure.facet_distribution",
        ),
        "target.results.figure.fit_scatter": (
            "surface.results",
            "focus.results.figure.fit_scatter",
            "panel.results.figure.fit_scatter",
        ),
        "target.results.figure.forest_measures": (
            "surface.results",
            "focus.results.figure.forest_measures",
            "panel.results.figure.forest_measures",
        ),
        "target.results.figure.misfit_ranking": (
            "surface.results",
            "focus.results.figure.misfit_ranking",
            "panel.results.figure.misfit_ranking",
        ),
        "target.results.figure.mml_person_sd": (
            "surface.results",
            "focus.results.figure.mml_person_sd",
            "panel.results.figure.mml_person_sd",
        ),
        "target.results.figure.pathway_map": (
            "surface.results",
            "focus.results.figure.pathway_map",
            "panel.results.figure.pathway_map",
        ),
        "target.results.figure.qq_residuals": (
            "surface.results",
            "focus.results.figure.qq_residuals",
            "panel.results.figure.qq_residuals",
        ),
        "target.results.figure.rater_agreement": (
            "surface.results",
            "focus.results.figure.rater_agreement",
            "panel.results.figure.rater_agreement",
        ),
        "target.results.figure.scree": (
            "surface.results",
            "focus.results.figure.scree",
            "panel.results.figure.scree",
        ),
        "target.results.figure.threshold_map": (
            "surface.results",
            "focus.results.figure.threshold_map",
            "panel.results.figure.threshold_map",
        ),
        "target.results.figure.wright_map": (
            "surface.results",
            "focus.results.figure.wright_map",
            "panel.results.figure.wright_map",
        ),
        "target.results.figure.zstd_distribution": (
            "surface.results",
            "focus.results.figure.zstd_distribution",
            "panel.results.figure.zstd_distribution",
        ),
        "target.results.first_read": (
            "surface.results",
            "focus.results.first_read",
            "panel.results.first_read",
        ),
        "target.results.fit": (
            "surface.results",
            "focus.results.fit",
            "panel.results.fit",
        ),
        "target.results.measures_targeting": (
            "surface.results",
            "focus.results.measures_targeting",
            "panel.results.measures_targeting",
        ),
        "target.results.rater_evidence": (
            "surface.results",
            "focus.results.rater_evidence",
            "panel.results.rater_evidence",
        ),
        "target.results.residual_structure": (
            "surface.results",
            "focus.results.residual_structure",
            "panel.results.residual_structure",
        ),
        "target.settings.analysis_scope": (
            "surface.settings",
            "focus.settings.analysis_scope",
            "panel.settings.analysis_scope",
        ),
        "target.settings.anchors": (
            "surface.settings",
            "focus.settings.anchors",
            "panel.settings.anchors",
        ),
        "target.settings.design_constraints": (
            "surface.settings",
            "focus.settings.design_constraints",
            "panel.settings.design_constraints",
        ),
        "target.settings.diagnostics": (
            "surface.settings",
            "focus.settings.diagnostics",
            "panel.settings.diagnostics",
        ),
        "target.settings.estimation": (
            "surface.settings",
            "focus.settings.estimation",
            "panel.settings.estimation",
        ),
        "target.settings.models": (
            "surface.settings",
            "focus.settings.models",
            "panel.settings.models",
        ),
    }
)


def _focus_id(target_id: str) -> str:
    try:
        return HELP_SURFACE_TARGET_CATALOG[target_id][1]
    except KeyError as exc:
        raise ValueError(f"Target is not declared by the Help surface catalog: {target_id}") from exc


HELP_TARGETS = tuple(
    HelpTarget(
        target_id=target_id,
        surface_id=surface_id,
        focus_ids=(focus_id,),
        presentation_state={"panel_id": panel_id, "surface_id": surface_id},
    )
    for target_id, (surface_id, focus_id, panel_id) in HELP_SURFACE_TARGET_CATALOG.items()
)


def _topic_suffix(help_topic_id: str) -> str:
    return help_topic_id.removeprefix("help.")


def _topic_context_policy(help_topic_id: str) -> HelpContextPolicy:
    if help_topic_id.startswith("help.results.") or help_topic_id in {
        "help.design.coverage",
        "help.run.estimation_failed",
        "help.run.nonconvergence",
        "help.report.claim_limits",
    }:
        return HelpContextPolicy.OPTIONAL_CURRENT
    return HelpContextPolicy.STATIC_ONLY


_HOME_LINKS = tuple(
    HelpLink(
        help_link_id=f"link.home.{_topic_suffix(help_topic_id)}",
        source_target_id="target.help.home",
        help_topic_id=help_topic_id,
        section_id="overview",
        return_target_id="target.help.home",
        return_focus_id=_focus_id("target.help.home"),
        context_policy=HelpContextPolicy.STATIC_ONLY,
    )
    for help_topic_id in CANONICAL_HELP_TOPIC_IDS
)

_CONTEXT_LINKS = tuple(
    HelpLink(
        help_link_id=f"link.context.{_topic_suffix(help_topic_id)}",
        source_target_id=_PRIMARY_TARGET_BY_TOPIC[help_topic_id],
        help_topic_id=help_topic_id,
        section_id="overview",
        return_target_id=_PRIMARY_TARGET_BY_TOPIC[help_topic_id],
        return_focus_id=_focus_id(_PRIMARY_TARGET_BY_TOPIC[help_topic_id]),
        context_policy=_topic_context_policy(help_topic_id),
    )
    for help_topic_id in CANONICAL_HELP_TOPIC_IDS
)

_PROBLEM_LINKS = tuple(
    HelpLink(
        help_link_id=f"link.problem.{problem_code}.{index}",
        source_target_id=target_id,
        help_topic_id=spec.help_topic_id,
        section_id="overview",
        return_target_id=target_id,
        return_focus_id=_focus_id(target_id),
        context_policy=(
            HelpContextPolicy.STATIC_ONLY
            if problem_code.startswith("input.") or problem_code == "problem.unexpected"
            else HelpContextPolicy.OPTIONAL_CURRENT
        ),
    )
    for problem_code, spec in USER_PROBLEM_SPECS.items()
    for index, target_id in enumerate(spec.action_target_ids, start=1)
)

_POPOVER_LINKS = tuple(
    HelpLink(
        help_link_id=f"link.popover.{popover_key}",
        source_target_id=_figure_target_id(popover_key),
        help_topic_id=help_topic_id,
        section_id=section_id,
        return_target_id=_figure_target_id(popover_key),
        return_focus_id=_focus_id(_figure_target_id(popover_key)),
        context_policy=HelpContextPolicy.REQUIRED_CURRENT,
    )
    for popover_key, help_topic_id, section_id in _POPOVER_DESTINATIONS
)

HELP_POPOVER_LINK_IDS: Mapping[str, str] = MappingProxyType(
    {
        popover_key: f"link.popover.{popover_key}"
        for popover_key, _, _ in _POPOVER_DESTINATIONS
    }
)

HELP_LINKS = (*_HOME_LINKS, *_CONTEXT_LINKS, *_PROBLEM_LINKS, *_POPOVER_LINKS)
HELP_ENTRY_TOPIC_IDS = CANONICAL_HELP_TOPIC_IDS
HELP_SYSTEM_TOPIC_IDS = (FALLBACK_HELP_TOPIC_ID,)

HELP_REGISTRY = HelpRegistry(
    topics=HELP_TOPICS,
    links=HELP_LINKS,
    targets=HELP_TARGETS,
    entry_topic_ids=HELP_ENTRY_TOPIC_IDS,
    system_topic_ids=HELP_SYSTEM_TOPIC_IDS,
    require_zero_orphans=True,
)


def required_help_locale_keys(
    topics: Iterable[HelpTopic] = HELP_TOPICS,
) -> frozenset[str]:
    """Return every locale key required to render the supplied topic metadata."""

    keys: set[str] = set()
    for topic in topics:
        keys.update((topic.title_key, topic.summary_key))
        keys.update(topic.prerequisite_keys)
        keys.update(topic.search_alias_keys)
        keys.update(
            key
            for key in (
                topic.computed_key,
                topic.can_show_key,
                topic.cannot_show_key,
                topic.next_check_key,
                topic.next_action_key,
                topic.safe_report_guard_key,
                topic.safe_report_key,
                topic.avoid_report_key,
            )
            if key is not None
        )
        for section in topic.sections:
            keys.add(section.title_key)
            keys.update(section.content_keys)
    return frozenset(keys)


HELP_REQUIRED_LOCALE_KEYS = required_help_locale_keys()
HELP_ALL_LOCALE_KEYS = HELP_REQUIRED_LOCALE_KEYS | HELP_NAV_LOCALE_KEYS


def validate_builtin_help_registry() -> None:
    """Validate the Help graph and its terminology/problem references."""

    concept_ids = tuple(term.concept_id for term in TERMINOLOGY_REGISTRY.terms)
    validate_help_registry(
        HELP_TOPICS,
        HELP_LINKS,
        HELP_TARGETS,
        entry_topic_ids=HELP_ENTRY_TOPIC_IDS,
        system_topic_ids=HELP_SYSTEM_TOPIC_IDS,
        known_concept_ids=concept_ids,
        known_claim_boundary_ids=DECLARED_HELP_CLAIM_BOUNDARY_IDS,
        known_reference_ids=DECLARED_HELP_REFERENCE_IDS,
        require_zero_orphans=True,
    )

    topic_ids = set(HELP_REGISTRY.topic_map)
    target_ids = set(HELP_REGISTRY.target_map)
    problem_topic_ids = {
        spec.help_topic_id for spec in USER_PROBLEM_SPECS.values()
    }
    problem_target_ids = {
        target_id
        for spec in USER_PROBLEM_SPECS.values()
        for target_id in spec.action_target_ids
    }
    terminology_topic_ids = {
        help_topic_id
        for term in TERMINOLOGY_REGISTRY.terms
        for help_topic_id in term.related_help_topic_ids
    }
    missing_problem_topics = sorted(problem_topic_ids.difference(topic_ids))
    missing_problem_targets = sorted(problem_target_ids.difference(target_ids))
    missing_terminology_topics = sorted(terminology_topic_ids.difference(topic_ids))
    referenced_target_ids = {
        "target.help.home",
        *problem_target_ids,
        *(_PRIMARY_TARGET_BY_TOPIC.values()),
        *(
            _figure_target_id(popover_key)
            for popover_key, _, _ in _POPOVER_DESTINATIONS
        ),
    }
    catalog_target_ids = set(HELP_SURFACE_TARGET_CATALOG)
    missing_catalog_targets = sorted(referenced_target_ids.difference(catalog_target_ids))
    unused_catalog_targets = sorted(catalog_target_ids.difference(referenced_target_ids))
    if (
        missing_problem_topics
        or missing_problem_targets
        or missing_terminology_topics
        or missing_catalog_targets
        or unused_catalog_targets
    ):
        raise ValueError(
            "Built-in Help registry has unresolved cross-contract references; "
            f"problem_topics={missing_problem_topics!r}, "
            f"problem_targets={missing_problem_targets!r}, "
            f"terminology_topics={missing_terminology_topics!r}, "
            f"missing_catalog_targets={missing_catalog_targets!r}, "
            f"unused_catalog_targets={unused_catalog_targets!r}"
        )

    fallback = HELP_REGISTRY.topic(FALLBACK_HELP_TOPIC_ID)
    if fallback is None or not fallback.has_section(FALLBACK_HELP_SECTION_ID):
        raise ValueError("Built-in Help registry must include its safe fallback topic")


validate_builtin_help_registry()


__all__ = [
    "CANONICAL_HELP_TOPIC_IDS",
    "DECLARED_HELP_CLAIM_BOUNDARY_IDS",
    "DECLARED_HELP_REFERENCE_IDS",
    "FALLBACK_HELP_SECTION_ID",
    "FALLBACK_HELP_TOPIC_ID",
    "HELP_ENTRY_TOPIC_IDS",
    "HELP_ALL_LOCALE_KEYS",
    "HELP_LINKS",
    "HELP_NAV_LOCALE_KEYS",
    "HELP_POPOVER_LINK_IDS",
    "HELP_REGISTRY",
    "HELP_REQUIRED_LOCALE_KEYS",
    "HELP_SURFACE_CATALOG_REVIEW_STATE",
    "HELP_SURFACE_TARGET_CATALOG",
    "HELP_SYSTEM_TOPIC_IDS",
    "HELP_TARGETS",
    "HELP_TOPICS",
    "HELP_TOPIC_REGISTRY_VERSION",
    "HELP_SUPPORT_CATALOG_REVIEW_STATE",
    "required_help_locale_keys",
    "validate_builtin_help_registry",
]
