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
