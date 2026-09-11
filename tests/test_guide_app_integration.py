from __future__ import annotations

from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app.guidance import GuideLifecycle, GuideNode, GuideRoute


APPTEST_TIMEOUT = 180


def _new_app() -> AppTest:
    return AppTest.from_file("streamlit_app.py").run(timeout=45)


def test_landing_owns_the_page_and_offers_three_optional_routes() -> None:
    at = _new_app()

    assert not at.exception
    assert [(button.key, button.label) for button in at.button] == [
        ("onboarding_quickstart", "Learn with a sample"),
        ("onboarding_dismiss", "Start with my data"),
        ("onboarding_skip_guide", "Continue without the guide"),
        ("mfrm_help_open_global", "Open Help"),
    ]
    assert {radio.key for radio in at.radio} == {"lang", "app_view_density"}
    assert not any(radio.key == "data_source_flat" for radio in at.radio)
    assert not any(radio.key == "data_source_class" for radio in at.radio)
    assert not any(
        selectbox.key and str(selectbox.key).startswith("facets_mode_")
        for selectbox in at.selectbox
    )
    state = at.session_state[app._GUIDANCE_STATE_KEY]
    assert state.lifecycle is GuideLifecycle.NOT_STARTED
    assert state.active_node_id is GuideNode.WELCOME


def test_own_data_route_opens_paste_preflight_without_starting_sample() -> None:
    at = _new_app()
    at.button(key="onboarding_dismiss").click()
    at.run(timeout=45)

    assert not at.exception
    state = at.session_state[app._GUIDANCE_STATE_KEY]
    assert state.lifecycle is GuideLifecycle.SKIPPED
    assert state.route_id is GuideRoute.OWN_DATA
    assert state.sample_context_id is None
    assert at.session_state["data_source_flat"] == "paste"
    assert any("pasted and uploaded rating files" in item.value for item in at.warning)


def test_language_switch_preserves_the_active_sample_node_and_bindings() -> None:
    at = _new_app()
    at.button(key="onboarding_quickstart").click()
    at.run(timeout=45)
    before = at.session_state[app._GUIDANCE_STATE_KEY]

    at.radio(key="lang").set_value("ja")
    at.run(timeout=45)

    assert not at.exception
    assert at.session_state[app._GUIDANCE_STATE_KEY] == before
    assert at.button(key="sample_guide_data_check_continue").label == "対応づけを確認して次へ"
    assert [(metric.label, metric.value) for metric in at.metric] == [
        ("評価数", "960"),
        ("対象者数", "30"),
        ("Facet数", "3"),
    ]


def test_sample_data_check_can_exit_and_resume_without_estimation() -> None:
    at = _new_app()
    at.button(key="onboarding_quickstart").click()
    at.run(timeout=45)

    state = at.session_state[app._GUIDANCE_STATE_KEY]
    assert state.active_node_id is GuideNode.DATA_CHECK
    assert [(metric.label, metric.value) for metric in at.metric] == [
        ("Ratings", "960"),
        ("People", "30"),
        ("Facets", "3"),
    ]
    assert "facets_mode_output" not in at.session_state

    at.button(key="sample_guide_exit_data_check").click()
    at.run(timeout=45)
    exited = at.session_state[app._GUIDANCE_STATE_KEY]
    assert exited.lifecycle is GuideLifecycle.SKIPPED
    assert exited.skip_origin_node_id is GuideNode.DATA_CHECK
    assert "facets_mode_output" not in at.session_state

    at.button(key="sample_guide_resume").click()
    at.run(timeout=45)
    resumed = at.session_state[app._GUIDANCE_STATE_KEY]
    assert resumed.lifecycle is GuideLifecycle.ACTIVE
    assert resumed.active_node_id is GuideNode.DATA_CHECK
    assert "facets_mode_output" not in at.session_state


def test_sample_context_restores_an_existing_workspace_on_exit() -> None:
    at = AppTest.from_string(
        """
import streamlit as st
import streamlit_app as app
from mfrm_app import guidance

app._apply_pending_guide_workspace_restore()
if not st.session_state.get("_guide_harness_seeded"):
    st.session_state["_guide_harness_seeded"] = True
    st.session_state["facets_mode_workflow_mode"] = "Advanced controls"
    st.session_state["facets_mode_output"] = {"sentinel": "real-result"}

if st.button("Start sample", key="harness_start"):
    app._start_sample_guide()
    st.rerun()
if st.button("Exit sample", key="harness_exit"):
    app._request_sample_workspace_restore()
    app._dispatch_guidance_event(guidance.GuidanceEventType.EXIT)
    st.rerun()
if st.button("Resume sample", key="harness_resume"):
    app._enter_sample_guide_workspace()
    app._dispatch_guidance_event(guidance.GuidanceEventType.RESUME)
    st.rerun()
""",
        default_timeout=45,
    ).run()

    at.button(key="harness_start").click()
    at.run()
    assert "facets_mode_output" not in at.session_state
    assert "facets_mode_workflow_mode" not in at.session_state
    assert app._GUIDE_SAVED_WORKSPACE_KEY in at.session_state

    at.button(key="harness_exit").click()
    at.run()
    assert at.session_state["facets_mode_output"] == {"sentinel": "real-result"}
    assert at.session_state["facets_mode_workflow_mode"] == "Advanced controls"
    assert app._GUIDE_SAVED_WORKSPACE_KEY not in at.session_state
    assert at.session_state[app._GUIDANCE_STATE_KEY].lifecycle is GuideLifecycle.SKIPPED

    at.button(key="harness_resume").click()
    at.run()
    assert "facets_mode_output" not in at.session_state
    assert "facets_mode_workflow_mode" not in at.session_state
    assert app._GUIDE_SAVED_WORKSPACE_KEY in at.session_state
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is GuideNode.DATA_CHECK

    at.button(key="harness_exit").click()
    at.run()
    assert at.session_state["facets_mode_output"] == {"sentinel": "real-result"}
    assert at.session_state["facets_mode_workflow_mode"] == "Advanced controls"


def test_complete_sample_guide_uses_real_fit_and_keeps_analysis_identity() -> None:
    at = _new_app()
    at.button(key="onboarding_quickstart").click()
    at.run(timeout=45)
    at.button(key="sample_guide_data_check_continue").click()
    at.run(timeout=45)
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is GuideNode.ESTIMATE

    at.button(key="sample_guide_estimate_run").click()
    at.run(timeout=APPTEST_TIMEOUT)
    assert not at.exception
    evidence_state = at.session_state[app._GUIDANCE_STATE_KEY]
    assert evidence_state.active_node_id is GuideNode.EVIDENCE_REVIEW
    assert evidence_state.bound_analysis_id == app.build_result_analysis_identity(
        at.session_state["facets_mode_output"]["result"]
    ).analysis_id

    at.radio(key="sample_guide_formative_answer").set_value("overclaim")
    at.button(key="sample_guide_formative_check").click()
    at.run(timeout=60)
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is GuideNode.EVIDENCE_REVIEW
    assert any("goes beyond the evidence" in item.value for item in at.error)

    at.radio(key="sample_guide_formative_answer").set_value("bounded")
    at.button(key="sample_guide_formative_check").click()
    at.run(timeout=60)
    assert at.session_state[app._GUIDANCE_STATE_KEY].active_node_id is GuideNode.ARCHIVE

    analysis_id = at.session_state[app._GUIDANCE_STATE_KEY].bound_analysis_id
    at.button(key="sample_guide_finish").click()
    at.run(timeout=60)
    completed = at.session_state[app._GUIDANCE_STATE_KEY]
    assert not at.exception
    assert completed.lifecycle is GuideLifecycle.COMPLETED
    assert completed.workflow_complete
    assert completed.learning_complete
    assert completed.bound_analysis_id == analysis_id
    assert app.build_result_analysis_identity(
        at.session_state["facets_mode_output"]["result"]
    ).analysis_id == analysis_id
    assert any("Sample guide completed" in item.value for item in at.success)
