from __future__ import annotations

from mfrm_app import ux


def test_data_source_classification_uses_stable_widget_ids() -> None:
    assert ux.classify_data_source("scenario:writing_essay") is ux.DataSourceKind.SAMPLE
    assert ux.classify_data_source("simulate") is ux.DataSourceKind.SIMULATION
    assert ux.classify_data_source("paste") is ux.DataSourceKind.PASTE
    assert ux.classify_data_source("upload") is ux.DataSourceKind.UPLOAD
    assert ux.classify_data_source("future-source") is ux.DataSourceKind.UNKNOWN


def test_privacy_policy_fails_closed_for_unregistered_sources() -> None:
    assert not ux.data_source_may_contain_user_data(ux.DataSourceKind.SAMPLE)
    assert not ux.data_source_may_contain_user_data(ux.DataSourceKind.SIMULATION)
    assert ux.data_source_may_contain_user_data(ux.DataSourceKind.PASTE)
    assert ux.data_source_may_contain_user_data(ux.DataSourceKind.UPLOAD)
    assert ux.data_source_may_contain_user_data("future-source")


def test_workspace_presentation_progressively_discloses_setup() -> None:
    choose = ux.workspace_presentation(has_data=False, has_result=False)
    review = ux.workspace_presentation(has_data=True, has_result=False)
    results = ux.workspace_presentation(has_data=True, has_result=True)

    assert choose.phase is ux.WorkspacePhase.CHOOSE_DATA
    assert not choose.show_pre_run_checks
    assert review.phase is ux.WorkspacePhase.REVIEW_AND_RUN
    assert review.show_pre_run_checks
    assert not review.collapse_input
    assert results.phase is ux.WorkspacePhase.READ_RESULTS
    assert results.collapse_input
    assert not results.show_pre_run_checks


def test_workflow_shell_has_one_primary_action_owner_for_every_phase() -> None:
    source = ux.resolve_workflow_shell(has_data=False, has_result=False)
    setup = ux.resolve_workflow_shell(has_data=True, has_result=False)
    running = ux.resolve_workflow_shell(
        has_data=True,
        has_result=False,
        estimation_running=True,
    )
    first_read = ux.resolve_workflow_shell(
        has_data=True,
        has_result=True,
        selected_result_section="first_read",
    )
    detail = ux.resolve_workflow_shell(
        has_data=True,
        has_result=True,
        selected_result_section="diagnostics",
    )
    archive = ux.resolve_workflow_shell(
        has_data=True,
        has_result=True,
        selected_result_section="report_export",
    )

    assert [shell.phase for shell in (source, setup, running, first_read, detail, archive)] == [
        ux.WorkflowPhase.SOURCE,
        ux.WorkflowPhase.SETUP,
        ux.WorkflowPhase.RUN,
        ux.WorkflowPhase.FIRST_READ,
        ux.WorkflowPhase.EVIDENCE_DETAIL,
        ux.WorkflowPhase.ARCHIVE,
    ]
    assert source.primary_action_surface is ux.WorkflowSurface.SOURCE_WORKSPACE
    assert setup.primary_action_surface is ux.WorkflowSurface.SETUP_WORKSPACE
    assert running.primary_action_surface is ux.WorkflowSurface.RUN_STATUS
    for shell in (first_read, detail, archive):
        assert shell.primary_action_surface is ux.WorkflowSurface.RESULT_GOAL_ROUTER
        assert shell.show_result_router
        assert shell.supporting_surfaces == (ux.WorkflowSurface.RESULT_SECTION,)
        assert shell.collapse_setup


def test_workflow_shell_fails_stale_result_section_back_to_first_read() -> None:
    shell = ux.resolve_workflow_shell(
        has_data=True,
        has_result=True,
        selected_result_section="future-section",
    )

    assert shell.phase is ux.WorkflowPhase.FIRST_READ
    assert shell.selected_result_section == "start"


def test_result_section_registry_has_stable_order_and_unique_ids() -> None:
    assert ux.RESULT_SECTION_IDS == (
        "start",
        "first_read",
        "results",
        "diagnostics",
        "figures",
        "report_export",
        "learn",
    )
    assert len(ux.RESULT_SECTION_IDS) == len(set(ux.RESULT_SECTION_IDS))
