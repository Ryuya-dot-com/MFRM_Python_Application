from __future__ import annotations

from dataclasses import fields
import inspect

import pytest

import mfrm_app.guidance as guidance_module
from mfrm_app.guidance import (
    DEFAULT_SAMPLE_CONTEXT_ID,
    GUIDE_NODES,
    GUIDE_NODE_ORDER,
    GuidanceContractError,
    GuidanceEvent,
    GuidanceEventType,
    GuidanceState,
    GuideLifecycle,
    GuideNode,
    GuideRoute,
    LearningState,
    reduce_guidance,
)


def _event(event_type: GuidanceEventType, **kwargs) -> GuidanceEvent:
    return GuidanceEvent(event_type=event_type, **kwargs)


def _sample_at_estimate() -> GuidanceState:
    started = reduce_guidance(GuidanceState.initial(), _event(GuidanceEventType.START_SAMPLE))
    return reduce_guidance(
        started,
        _event(
            GuidanceEventType.NODE_REVIEWED,
            data_fingerprint="data.sample.v1",
        ),
    )


def _sample_at_evidence() -> GuidanceState:
    return reduce_guidance(
        _sample_at_estimate(),
        _event(
            GuidanceEventType.FIT_SUCCEEDED,
            analysis_id="analysis.sample.1",
            data_fingerprint="data.sample.v1",
        ),
    )


def test_catalog_uses_five_unique_stable_nodes_and_locale_keys() -> None:
    assert GUIDE_NODE_ORDER == (
        GuideNode.WELCOME,
        GuideNode.DATA_CHECK,
        GuideNode.ESTIMATE,
        GuideNode.EVIDENCE_REVIEW,
        GuideNode.ARCHIVE,
    )
    assert [item.ordinal for item in GUIDE_NODES] == [1, 2, 3, 4, 5]
    assert len({item.node_id for item in GUIDE_NODES}) == 5
    assert all(item.title_key.startswith("guide.") for item in GUIDE_NODES)
    assert all(item.body_key.startswith("guide.") for item in GUIDE_NODES)


def test_guidance_core_has_no_streamlit_pandas_or_estimator_dependency() -> None:
    source = inspect.getsource(guidance_module)

    assert "import streamlit" not in source
    assert "import pandas" not in source
    assert "mfrm_estimate" not in source
    assert "claim_ready" not in source


def test_initial_state_is_unbound_and_focuses_welcome() -> None:
    state = GuidanceState.initial()

    assert state.lifecycle is GuideLifecycle.NOT_STARTED
    assert state.route_id is None
    assert state.active_node_id is GuideNode.WELCOME
    assert state.node_progress(GuideNode.WELCOME).state is LearningState.ACTIVE
    assert all(
        state.node_progress(node).state is LearningState.LOCKED
        for node in GUIDE_NODE_ORDER[1:]
    )
    assert "claim_ready" not in {field.name for field in fields(GuidanceState)}


def test_sample_start_binds_isolated_context_without_analysis() -> None:
    state = reduce_guidance(GuidanceState.initial(), _event(GuidanceEventType.START_SAMPLE))

    assert state.lifecycle is GuideLifecycle.ACTIVE
    assert state.route_id is GuideRoute.SAMPLE
    assert state.active_node_id is GuideNode.DATA_CHECK
    assert state.sample_context_id == DEFAULT_SAMPLE_CONTEXT_ID
    assert state.bound_analysis_id is None
    assert state.node_progress(GuideNode.WELCOME).state is LearningState.COMPLETE


def test_own_data_start_has_no_sample_context() -> None:
    state = reduce_guidance(GuidanceState.initial(), _event(GuidanceEventType.START_OWN_DATA))

    assert state.route_id is GuideRoute.OWN_DATA
    assert state.active_node_id is GuideNode.DATA_CHECK
    assert state.sample_context_id is None


def test_data_review_and_fit_advance_without_promoting_scientific_state() -> None:
    estimate = _sample_at_estimate()
    evidence = _sample_at_evidence()

    assert estimate.active_node_id is GuideNode.ESTIMATE
    assert estimate.bound_data_fingerprint == "data.sample.v1"
    assert evidence.active_node_id is GuideNode.EVIDENCE_REVIEW
    assert evidence.bound_analysis_id == "analysis.sample.1"
    assert not evidence.workflow_complete
    assert not evidence.learning_complete


def test_formative_answer_preserves_analysis_and_cannot_create_claim_readiness() -> None:
    evidence = _sample_at_evidence()
    archive = reduce_guidance(
        evidence,
        _event(
            GuidanceEventType.FORMATIVE_ANSWERED,
            learning_evidence_id="learning.answer.bounded_claim",
            evidence_ids=("evidence.sample.convergence",),
        ),
    )

    assert archive.active_node_id is GuideNode.ARCHIVE
    assert archive.bound_analysis_id == evidence.bound_analysis_id
    assert archive.reviewed_evidence_ids == ("evidence.sample.convergence",)
    assert "claim_ready" not in archive.to_payload()


def test_archive_completion_is_learning_and_workflow_completion_only() -> None:
    archive = reduce_guidance(
        _sample_at_evidence(),
        _event(
            GuidanceEventType.FORMATIVE_ANSWERED,
            learning_evidence_id="learning.answer.bounded_claim",
        ),
    )
    completed = reduce_guidance(
        archive,
        _event(GuidanceEventType.ARCHIVE_CREATED),
    )

    assert completed.lifecycle is GuideLifecycle.COMPLETED
    assert completed.active_node_id is None
    assert completed.workflow_complete
    assert completed.learning_complete
    assert all(item.state is LearningState.COMPLETE for item in completed.node_states)


@pytest.mark.parametrize(
    "event_type",
    [
        GuidanceEventType.VIEW_CHANGED,
        GuidanceEventType.LANGUAGE_CHANGED,
        GuidanceEventType.HELP_OPENED,
        GuidanceEventType.HELP_RETURNED,
    ],
)
def test_presentation_events_are_identity_preserving_noops(event_type) -> None:
    state = _sample_at_evidence()

    assert reduce_guidance(state, _event(event_type)) is state


def test_exit_resume_and_restart_preserve_fit_binding() -> None:
    evidence = _sample_at_evidence()
    skipped = reduce_guidance(evidence, _event(GuidanceEventType.EXIT))
    resumed = reduce_guidance(skipped, _event(GuidanceEventType.RESUME))
    restarted = reduce_guidance(skipped, _event(GuidanceEventType.RESTART))

    assert skipped.lifecycle is GuideLifecycle.SKIPPED
    assert skipped.skip_origin_node_id is GuideNode.EVIDENCE_REVIEW
    assert resumed.active_node_id is GuideNode.EVIDENCE_REVIEW
    assert resumed.bound_analysis_id == evidence.bound_analysis_id
    assert restarted.active_node_id is GuideNode.DATA_CHECK
    assert restarted.bound_analysis_id == evidence.bound_analysis_id
    assert restarted.node_progress(GuideNode.ESTIMATE).state is LearningState.LOCKED


def test_skip_from_welcome_opens_normal_route_without_sample_context() -> None:
    skipped = reduce_guidance(GuidanceState.initial(), _event(GuidanceEventType.SKIP))

    assert skipped.lifecycle is GuideLifecycle.SKIPPED
    assert skipped.route_id is GuideRoute.WITHOUT_GUIDE
    assert skipped.skip_origin_node_id is GuideNode.WELCOME
    assert skipped.sample_context_id is None
    assert reduce_guidance(skipped, _event(GuidanceEventType.RESUME)) is skipped


def test_data_and_spec_changes_reopen_first_affected_node() -> None:
    evidence = _sample_at_evidence()
    data_changed = reduce_guidance(
        evidence,
        _event(
            GuidanceEventType.DATA_OR_MAPPING_CHANGED,
            data_fingerprint="data.sample.v2",
        ),
    )
    spec_changed = reduce_guidance(evidence, _event(GuidanceEventType.DRAFT_SPEC_CHANGED))

    assert data_changed.active_node_id is GuideNode.DATA_CHECK
    assert data_changed.bound_data_fingerprint == "data.sample.v2"
    assert data_changed.bound_analysis_id is None
    assert spec_changed.active_node_id is GuideNode.ESTIMATE
    assert spec_changed.bound_analysis_id is None
    assert spec_changed.node_progress(GuideNode.DATA_CHECK).state is LearningState.COMPLETE


def test_fit_success_requires_an_analysis_reference() -> None:
    with pytest.raises(GuidanceContractError, match="requires analysis_id"):
        reduce_guidance(
            _sample_at_estimate(),
            _event(GuidanceEventType.FIT_SUCCEEDED),
        )


def test_fit_failure_records_learning_but_stays_at_estimate_without_binding() -> None:
    estimate = _sample_at_estimate()
    failed = reduce_guidance(
        estimate,
        _event(
            GuidanceEventType.FIT_FAILED,
            learning_evidence_id="learning.fit_nonconvergence_reviewed",
        ),
    )

    assert failed.active_node_id is GuideNode.ESTIMATE
    assert failed.bound_analysis_id is None
    assert failed.node_progress(GuideNode.ESTIMATE).state is LearningState.ACTIVE
    assert failed.node_progress(GuideNode.ESTIMATE).learning_evidence_ids == (
        "learning.fit_nonconvergence_reviewed",
    )


def test_state_payload_is_deterministic_and_contains_only_stable_values() -> None:
    state = _sample_at_evidence()

    assert state.to_payload() == state.to_payload()
    assert state.to_payload()["schema_version"] == "mfrm_guidance_v1"
    assert state.to_payload()["active_node_id"] == "evidence_review"
    assert state.to_payload()["route_id"] == "sample"
