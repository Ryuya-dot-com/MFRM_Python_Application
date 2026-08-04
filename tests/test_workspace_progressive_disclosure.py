from __future__ import annotations

import inspect

import streamlit_app as app


def test_setup_workspace_collapses_raw_rows_and_pre_run_checks_after_results() -> None:
    source = inspect.getsource(app.render_analysis_setup_workspace)

    assert "resolve_workflow_shell" in source
    assert "shell.collapse_setup" in source
    assert "app.setup_after_run_expander" in source
    assert "app.input_preview_rows_expander_template" in source
    assert "shell.show_pre_run_checks" in source


def test_successful_run_clears_setup_surface_before_results() -> None:
    source = inspect.getsource(app.run_facets_mode)

    assert "setup_surface = render_analysis_setup_workspace" in source
    assert "setup_surface.empty()" in source


def test_onboarding_is_action_first_and_supporting_steps_are_collapsed() -> None:
    source = inspect.getsource(app.render_onboarding_banner)

    assert "onboarding.title" in source
    assert 'type="primary"' in source
    assert 'expanded=False' in source
    assert "guided_first_run_route_table" not in source


def test_result_router_keeps_only_one_primary_action_above_supporting_detail() -> None:
    source = inspect.getsource(app._render_guided_goal_router)

    primary = "_render_action_card(primary_row, primary=True)"
    supporting = 'with st.expander(t("guided.goal_supporting_detail_expander")'
    secondary = "_render_action_card(hub_row, primary=False)"
    assert primary in source
    assert supporting in source
    assert secondary in source
    assert source.index(primary) < source.index(supporting) < source.index(secondary)


def test_default_result_route_has_no_duplicate_first_read_overview() -> None:
    source = inspect.getsource(app.run_facets_mode)

    assert source.count("_render_guided_goal_router(") == 1
    assert "workflow_shell.show_result_router" in source
    assert 'st.subheader(t("guided.overview_subheader"))' not in source
    assert "_render_guided_action_plan(" not in source

    section_source = inspect.getsource(app._render_guided_essential_tabs)
    assert "_render_guided_action_plan(" in section_source


def test_result_section_labels_cover_the_workflow_shell_registry() -> None:
    assert tuple(app.GUIDED_SECTION_I18N_KEYS) == app.GUIDED_SECTION_IDS
