"""Contract tests for audience-layer terminology and seed concepts."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import pytest

from mfrm_app.terminology import (
    AudienceLayer,
    INTERNAL_TECHNICAL_LABELS,
    LayerNotAllowedError,
    ReviewState,
    SCIENTIFIC_CONCEPT_IDS,
    TERMINOLOGY_REGISTRY,
    TERMINOLOGY_REGISTRY_VERSION,
    TermProjection,
    TermSpec,
    TerminologyRegistry,
    TerminologyValidationError,
    UnknownConceptError,
    VersionReviewMetadata,
    validate_term_registry,
)


MODULE_PATH = Path(__file__).resolve().parents[1] / "mfrm_app" / "terminology.py"

PLANNED_HELP_TOPIC_IDS = {
    "help.downloads.privacy",
    "help.downloads.repeat_analysis",
    "help.methods.jmle_mml",
    "help.methods.rsm_pcm_gpcm",
    "help.report.claim_limits",
    "help.results.first_read",
    "help.results.fit",
    "help.results.measures_targeting",
    "help.results.residual_structure",
}


def _term(concept_id: str) -> TermSpec:
    return TERMINOLOGY_REGISTRY.require(concept_id)


def test_default_registry_is_valid_unique_and_deterministic():
    validated = validate_term_registry(tuple(TERMINOLOGY_REGISTRY))
    concept_ids = tuple(term.concept_id for term in validated)

    assert TERMINOLOGY_REGISTRY.version == TERMINOLOGY_REGISTRY_VERSION
    assert len(concept_ids) == len(set(concept_ids))
    assert concept_ids == tuple(sorted(concept_ids))
    assert len(concept_ids) >= 20


def test_seed_registry_contains_internal_and_scientific_concepts():
    registered = {term.concept_id for term in TERMINOLOGY_REGISTRY}

    assert set(INTERNAL_TECHNICAL_LABELS) <= registered
    assert SCIENTIFIC_CONCEPT_IDS <= registered
    assert {
        "estimation.jmle",
        "estimation.mml",
        "fit.infit",
        "fit.outfit",
        "model.mfrm",
        "residual_structure.pca",
    } <= registered


def test_seed_terms_reference_only_planned_help_topic_ids():
    referenced = {
        help_topic_id
        for term in TERMINOLOGY_REGISTRY
        for help_topic_id in term.related_help_topic_ids
    }

    assert referenced == PLANNED_HELP_TOPIC_IDS


def test_every_term_has_locale_alias_boundary_help_and_review_metadata():
    for term in TERMINOLOGY_REGISTRY:
        assert term.definition_key.startswith("terms.")
        assert term.aliases_keys
        assert term.aliases == term.aliases_keys
        assert term.prohibited_claim_keys
        assert term.related_help_topic_ids
        assert term.locale_keys
        assert term.version_and_review.content_version == "1.0.0"
        assert term.version_and_review.review_state is ReviewState.IN_REVIEW
        assert term.version_and_review.reviewer_roles


def test_internal_exact_labels_are_retained_only_by_technical_projection():
    for concept_id, exact_label in INTERNAL_TECHNICAL_LABELS.items():
        standard = TERMINOLOGY_REGISTRY.project(concept_id, AudienceLayer.STANDARD)
        technical = TERMINOLOGY_REGISTRY.project(
            concept_id, AudienceLayer.TECHNICAL
        )

        assert standard.label_key is not None
        assert standard.literal_label is None
        assert standard.label_reference != exact_label
        assert exact_label.casefold() not in standard.label_reference.casefold()
        assert technical.label_key is None
        assert technical.literal_label == exact_label
        assert technical.label_reference == exact_label

        # A renderer receiving a Standard projection cannot inspect labels or
        # aliases from another layer by following fields on the projection.
        assert not hasattr(standard, "technical_label")
        assert not hasattr(standard, "aliases_keys")


def test_scientific_concepts_are_reviewed_for_method_projection():
    expected_labels = {
        "model.mfrm": "MFRM",
        "estimation.jmle": "JMLE",
        "estimation.mml": "MML",
        "residual_structure.pca": "residual PCA",
        "fit.infit": "Infit",
        "fit.outfit": "Outfit",
    }

    for concept_id in SCIENTIFIC_CONCEPT_IDS:
        term = _term(concept_id)
        projection = term.project(AudienceLayer.METHOD)
        assert AudienceLayer.METHOD in term.allowed_layers
        assert projection.label_key == term.method_label_key
        assert projection.literal_label is None
        assert projection.prohibited_claim_keys

    for concept_id, exact_label in expected_labels.items():
        assert _term(concept_id).project("technical").literal_label == exact_label


def test_residual_pca_boundary_does_not_encode_unidimensionality_proof():
    projection = TERMINOLOGY_REGISTRY.project(
        "residual_structure.pca", AudienceLayer.METHOD
    )

    assert any(
        "not_unidimensionality_proof" in key
        for key in projection.prohibited_claim_keys
    )


def test_unreviewed_layer_projection_fails_closed():
    residual_pca = _term("residual_structure.pca")

    assert not residual_pca.allows(AudienceLayer.STANDARD)
    with pytest.raises(LayerNotAllowedError, match="standard"):
        residual_pca.project(AudienceLayer.STANDARD)


def test_registry_rejects_duplicate_concept_ids():
    first = _term("estimation.jmle")
    second = replace(_term("estimation.mml"), concept_id=first.concept_id)

    with pytest.raises(TerminologyValidationError, match="Duplicate concept_id"):
        TerminologyRegistry(
            version="test_terminology_v1",
            terms=(first, second),
        )


def test_registry_rejects_alias_keys_shared_by_two_concepts():
    first = _term("fit.infit")
    second_source = _term("fit.outfit")
    second = replace(
        second_source,
        aliases_keys=(first.aliases_keys[0], second_source.aliases_keys[1]),
    )

    with pytest.raises(TerminologyValidationError, match="Alias key .* is shared"):
        validate_term_registry((first, second))


def test_registry_rejects_duplicate_technical_labels_case_insensitively():
    first = _term("estimation.jmle")
    second = replace(_term("estimation.mml"), technical_label="jmle")

    with pytest.raises(TerminologyValidationError, match="Technical label"):
        validate_term_registry((first, second))


def test_term_rejects_allowed_layer_without_a_corresponding_label():
    source = _term("fit.infit")

    with pytest.raises(TerminologyValidationError, match="allows 'method'"):
        replace(source, method_label_key=None)


def test_term_rejects_label_for_a_layer_not_declared_allowed():
    source = _term("fit.infit")

    with pytest.raises(TerminologyValidationError, match="does not allow"):
        replace(
            source,
            standard_label_key="terms.infit.standard_label",
        )


def test_term_rejects_invalid_ids_and_locale_keys():
    source = _term("fit.infit")

    with pytest.raises(TerminologyValidationError, match="namespaced ID"):
        replace(source, concept_id="Infit")
    with pytest.raises(TerminologyValidationError, match="below 'terms.'"):
        replace(source, definition_key="help.infit.definition")


def test_term_rejects_duplicate_keys_inside_one_field():
    source = _term("fit.infit")

    with pytest.raises(TerminologyValidationError, match="duplicate keys"):
        replace(source, aliases_keys=(source.aliases_keys[0],) * 2)


def test_reviewed_metadata_requires_an_independent_reviewer_role():
    with pytest.raises(TerminologyValidationError, match="at least one reviewer"):
        VersionReviewMetadata(
            content_version="1.0.0",
            owner_role="measurement_method_owner",
            last_reviewed=_term("fit.infit").version_and_review.last_reviewed,
            review_state=ReviewState.REVIEWED,
        )

    with pytest.raises(TerminologyValidationError, match="must not also"):
        VersionReviewMetadata(
            content_version="1.0.0",
            owner_role="measurement_method_owner",
            last_reviewed=_term("fit.infit").version_and_review.last_reviewed,
            review_state=ReviewState.REVIEWED,
            reviewer_roles=("measurement_method_owner",),
        )


def test_registry_lookup_and_alias_resolution_are_exact():
    analysis = _term("workflow.analysis_reference")
    alias_key = analysis.aliases_keys[0]

    assert TERMINOLOGY_REGISTRY.get(analysis.concept_id) is analysis
    assert TERMINOLOGY_REGISTRY.resolve_alias_key(alias_key) is analysis
    assert TERMINOLOGY_REGISTRY.resolve_alias_key(alias_key.upper()) is None
    assert TERMINOLOGY_REGISTRY.resolve_alias_key("analysis") is None
    with pytest.raises(UnknownConceptError):
        TERMINOLOGY_REGISTRY.require("workflow.not_registered")


def test_project_all_returns_only_terms_allowed_for_that_layer():
    standard = TERMINOLOGY_REGISTRY.project_all(AudienceLayer.STANDARD)
    method = TERMINOLOGY_REGISTRY.project_all(AudienceLayer.METHOD)
    technical = TERMINOLOGY_REGISTRY.project_all(AudienceLayer.TECHNICAL)

    assert standard
    assert method
    assert len(technical) == len(TERMINOLOGY_REGISTRY)
    assert {item.concept_id for item in method} == SCIENTIFIC_CONCEPT_IDS
    assert all(item.layer is AudienceLayer.STANDARD for item in standard)
    assert all(item.label_key is not None for item in standard)
    assert all(item.literal_label is None for item in standard)


def test_contract_objects_are_frozen():
    term = _term("workflow.analysis_reference")
    projection = term.project(AudienceLayer.STANDARD)

    with pytest.raises(FrozenInstanceError):
        term.concept_id = "workflow.changed"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        projection.label_key = "terms.changed"  # type: ignore[misc]
    with pytest.raises(FrozenInstanceError):
        TERMINOLOGY_REGISTRY.version = "changed_v1"  # type: ignore[misc]


def test_projection_constructor_enforces_layer_and_contract_fields():
    source = _term("workflow.analysis_reference").project(AudienceLayer.STANDARD)

    with pytest.raises(TerminologyValidationError, match="localized projection"):
        TermProjection(
            concept_id=source.concept_id,
            layer=AudienceLayer.STANDARD,
            label_key=None,
            literal_label="AnalysisID",
            definition_key=source.definition_key,
            prohibited_claim_keys=source.prohibited_claim_keys,
            related_help_topic_ids=source.related_help_topic_ids,
            content_version=source.content_version,
            review_state=source.review_state,
        )
    with pytest.raises(TerminologyValidationError, match="below 'terms.'"):
        replace(source, definition_key="help.invalid.definition")
    with pytest.raises(TerminologyValidationError, match="at least one key"):
        replace(source, prohibited_claim_keys=())
    with pytest.raises(TerminologyValidationError, match="MAJOR.MINOR.PATCH"):
        replace(source, content_version="draft")


def test_terminology_core_has_no_streamlit_or_pandas_import():
    tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"), filename=str(MODULE_PATH))
    imported_roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported_roots.add(node.module.split(".", 1)[0])

    assert imported_roots.isdisjoint({"streamlit", "pandas"})
