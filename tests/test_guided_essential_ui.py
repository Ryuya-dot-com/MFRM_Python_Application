from __future__ import annotations

import inspect
from types import SimpleNamespace

import pandas as pd

import streamlit_app as app
from mfrm_app import help_popovers


def test_guided_action_plan_preserves_detail_and_routes_priority():
    result = {
        "opt": SimpleNamespace(success=False, message="maximum iterations reached"),
        "config": {"model": "PCM"},
    }
    diagnostics = {
        "obs": pd.DataFrame({"StdResidual": [0.1, 2.4, -0.3, 0.0]}),
        "reliability": pd.DataFrame({"Facet": ["person"], "Reliability": [0.63]}),
        "pca_enabled": False,
        "marginal_fit": {},
    }

    plan = app.build_guided_action_plan(result, diagnostics, all_bias_results={})

    required = {
        "Priority",
        "SourceOrder",
        "Area",
        "Check",
        "Status",
        "WhatItSays",
        "NextAction",
        "DetailLocation",
    }
    assert required.issubset(plan.columns)
    assert not plan.empty
    assert plan.iloc[0]["Status"] == "Do not interpret yet"
    assert plan.iloc[0]["DetailLocation"] == "Report"
    assert plan["WhatItSays"].astype(str).str.len().gt(0).all()
    assert plan["NextAction"].astype(str).str.len().gt(0).all()


def test_guided_detail_navigator_keeps_six_tab_path_visible():
    navigator = app._guided_detail_navigator_table()

    assert list(navigator.columns) == ["Question", "StartHere", "ThenOpen", "Why"]
    assert set(navigator["StartHere"]) >= {
        app.t("guided.tab_first_read"),
        app.t("guided.tab_results"),
        app.t("guided.tab_diagnostics"),
        app.t("guided.tab_figures"),
        app.t("guided.tab_report_export"),
        app.t("guided.tab_learn"),
    }


def test_guided_measures_reading_guide_separates_inspection_from_claims():
    guide = app.guided_measures_reading_guide_table()

    assert list(guide.columns) == [
        app.t("guided.measures_col_surface"),
        app.t("guided.measures_col_read_first"),
        app.t("guided.measures_col_answers"),
        app.t("guided.measures_col_do_not"),
        app.t("guided.measures_col_safe_use"),
        app.t("guided.measures_col_open_detail"),
    ]
    assert len(guide) == 6
    joined = " ".join(guide.astype(str).to_numpy().ravel().tolist())
    assert app.t("guided.measures_combined_surface") in set(guide[app.t("guided.measures_col_surface")])
    assert app.t("guided.measures_person_surface") in set(guide[app.t("guided.measures_col_surface")])
    assert "Do not" in joined
    assert "First Read" in joined
    assert "Diagnostics" in joined or "diagnostics" in joined
    assert guide.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_results_section_renders_measures_reading_guide():
    source = inspect.getsource(app._render_guided_results_section)

    assert "render_guided_measures_reading_guide_tables()" in source
    assert "guided.measures_reading_heading" in source


def test_guided_measures_reading_guide_renders_in_readable_slices():
    source = inspect.getsource(app.render_guided_measures_reading_guide_tables)

    assert "guided_measures_reading_guide_table()" in source
    assert "guided.measures_inspection_heading" in source
    assert "guided.measures_boundary_heading" in source
    assert "guided.measures_route_heading" in source
    assert "guided.measures_col_do_not" in source
    assert "guided.measures_col_safe_use" in source


def test_help_section_renders_measures_reading_guide_reference():
    source = inspect.getsource(app.show_help_section)

    assert "render_guided_measures_reading_guide_tables()" in source
    assert "help.measures_reading_expander" in source


def test_help_popover_japanese_overlay_covers_every_topic():
    missing = help_popovers.missing_translation_fields(
        app._HELP_POPOVER_LIBRARY.keys(),
        lang="ja",
    )

    assert missing == []
    fallback = app._HELP_POPOVER_LIBRARY["category_usage"]
    localized = help_popovers.localized_help_popover_topic(
        "category_usage",
        lang="ja",
        fallback=fallback,
    )
    assert localized["title"] != fallback["title"]
    assert "カテゴリ" in localized["title"]
    assert localized["what"]
    assert localized["how"]
    assert localized["watch"]


def test_render_help_popover_uses_localized_topic_overlay():
    source = inspect.getsource(app.render_help_popover)

    assert "_help_popovers.localized_help_popover_topic" in source
    assert "st.session_state.get(\"lang\", DEFAULT_LANG)" in source


def test_language_switch_uses_lightweight_result_fast_path():
    source = inspect.getsource(app.main)

    assert "on_change=_mark_language_switch" in source
    assert "_consume_language_switch_fast_path(app_mode)" in source
    assert "render_language_switch_fast_path()" in source


def test_language_fast_path_consumes_pending_state_only_for_results(monkeypatch):
    state = {}

    class _FakeSt:
        session_state = state

    monkeypatch.setattr(app, "st", _FakeSt())
    state["_language_switch_pending"] = True
    assert app._consume_language_switch_fast_path("FACETS-mode estimation") is False
    assert state["_language_switch_pending"] is False

    state["_language_switch_pending"] = True
    state["facets_mode_output"] = {"result": {}}
    assert app._consume_language_switch_fast_path("FACETS-mode estimation") is True
    assert state["_language_switch_pending"] is False


def test_guided_section_selector_keeps_all_detail_sections_in_order():
    options = app.guided_section_selector_options()

    assert list(options["SectionId"]) == [
        "start",
        "first_read",
        "results",
        "diagnostics",
        "figures",
        "report_export",
        "learn",
    ]
    assert list(options.columns) == ["SectionId", app.t("guided.section_col_section")]
    assert list(options[app.t("guided.section_col_section")]) == [
        app.t("guided.tab_start"),
        app.t("guided.tab_first_read"),
        app.t("guided.tab_results"),
        app.t("guided.tab_diagnostics"),
        app.t("guided.tab_figures"),
        app.t("guided.tab_report_export"),
        app.t("guided.tab_learn"),
    ]


def test_guided_section_reading_order_covers_every_essential_section():
    expected_columns = [
        app.t("guided.section_col_section"),
        app.t("guided.goal_focus_col_open"),
        app.t("guided.next_click_col_do"),
        app.t("guided.interpret_col_do_not_conclude"),
    ]
    for section_id in app.GUIDED_SECTION_IDS:
        cue = app.guided_section_reading_order_table(section_id)

        assert cue.shape == (1, 4)
        assert list(cue.columns) == expected_columns
        assert cue.iloc[0][expected_columns[0]] == app._guided_section_label(section_id)
        assert cue.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()
        assert "Do not" in cue.iloc[0][expected_columns[3]]

    fallback = app.guided_section_reading_order_table("not_a_section")
    start = app.guided_section_reading_order_table("start")
    assert fallback.equals(start)


def test_guided_essential_tabs_show_reading_order_before_section_body():
    source = inspect.getsource(app._render_guided_essential_tabs)

    cue_call = "render_guided_section_reading_order(selected_section)"
    section_body = "if selected_section == \"start\":"
    assert cue_call in source
    assert source.index(cue_call) < source.index(section_body)


def test_guided_diagnostics_selector_lazy_renders_one_panel():
    options = app.guided_diagnostic_panel_options()

    assert list(options["PanelId"]) == [
        "fit_details",
        "dimensionality",
        "wright_map",
        "visuals",
        "bias_interaction",
        "categories_steps",
        "agreement",
        "facet_dashboard",
        "prediction_simulation",
    ]
    assert list(options.columns) == [
        "PanelId",
        app.t("guided.diagnostics_panel_col_panel"),
    ]

    source = inspect.getsource(app._render_guided_diagnostics_section)
    assert "st.segmented_control" in source
    assert "guided_diagnostics_panel" in source
    assert "st.tabs" not in source


def test_guided_figures_selector_lazy_renders_one_panel():
    options = app.guided_figure_panel_options()

    assert list(options["PanelId"]) == [
        "wright_map",
        "visuals",
        "export_status",
    ]
    assert list(options.columns) == [
        "PanelId",
        app.t("guided.figures_panel_col_panel"),
    ]

    source = inspect.getsource(app._render_guided_figures_section)
    assert "st.segmented_control" in source
    assert "guided_figures_panel" in source
    assert "st.tabs" not in source


def test_guided_first_run_route_table_links_sample_run_to_current_focus():
    table = app.guided_first_run_route_table()

    assert list(table.columns) == [
        app.t("guided.first_run_route_col_step"),
        app.t("guided.first_run_route_col_action"),
        app.t("guided.first_run_route_col_where"),
        app.t("guided.first_run_route_col_why"),
        app.t("guided.first_run_route_col_done"),
    ]
    assert len(table) == 5
    assert list(table[app.t("guided.first_run_route_col_step")]) == [
        app.t(f"guided.first_run_route_{step}_step") for step in range(1, 6)
    ]
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())
    assert app.t("guided.goal_focus_heading").split()[0] in joined
    assert app.t("guided.first_run_route_2_where") in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_action_plan_display_uses_reader_facing_columns():
    result = {
        "opt": SimpleNamespace(success=False, message="maximum iterations reached"),
        "config": {"model": "PCM"},
    }
    diagnostics = {
        "obs": pd.DataFrame({"StdResidual": [0.1, 2.4, -0.3, 0.0]}),
        "reliability": pd.DataFrame({"Facet": ["person"], "Reliability": [0.63]}),
        "pca_enabled": False,
        "marginal_fit": {},
    }
    plan = app.build_guided_action_plan(result, diagnostics, all_bias_results={})

    compact = app.guided_action_plan_display_frame(plan, include_detail=False)
    detailed = app.guided_action_plan_display_frame(plan, include_detail=True)

    assert app.t("guided.col_priority") in compact.columns
    assert app.t("guided.col_raw_status") not in compact.columns
    assert app.t("guided.col_raw_status") in detailed.columns
    assert compact.iloc[0][app.t("guided.col_status")] == app.t("guided.status_do_not_interpret")
    assert compact.iloc[0][app.t("guided.col_detail_location")] == app.t("guided.detail_location_report")


def test_guided_status_legend_explains_every_status_level():
    legend = app.guided_status_legend_table()

    assert list(legend.columns) == [
        app.t("guided.col_status"),
        app.t("guided.legend_meaning_col"),
        app.t("guided.legend_action_col"),
    ]
    assert set(legend[app.t("guided.col_status")]) >= {
        app.t("guided.status_do_not_interpret"),
        app.t("guided.status_caution"),
        app.t("guided.status_review"),
        app.t("guided.status_skipped"),
        app.t("guided.status_ok"),
    }
    assert legend[app.t("guided.legend_meaning_col")].astype(str).str.len().gt(0).all()
    assert legend[app.t("guided.legend_action_col")].astype(str).str.len().gt(0).all()


def test_guided_interpretation_readiness_summary_pauses_for_hard_stop():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Estimation",
        "Check": "Convergence",
        "Status": "Do not interpret yet",
        "WhatItSays": "Optimizer did not converge.",
        "NextAction": "Fix estimation first.",
        "DetailLocation": "Report",
    }])

    summary = app.guided_interpretation_readiness_summary_table(plan)
    row = summary.iloc[0]

    assert list(summary.columns) == [
        app.t("guided.interpret_col_status"),
        app.t("guided.interpret_col_main_item"),
        app.t("guided.interpret_col_open_next"),
        app.t("guided.interpret_col_do_not_conclude"),
        app.t("guided.interpret_col_safe_output"),
    ]
    assert row[app.t("guided.interpret_col_status")] == app.t("guided.interpret_status_pause")
    assert row[app.t("guided.interpret_col_open_next")] == app.t("guided.detail_location_report")
    assert "Optimizer did not converge." in row[app.t("guided.interpret_col_main_item")]
    assert row[app.t("guided.interpret_col_do_not_conclude")] == app.t("guided.interpret_pause_do_not")


def test_guided_interpretation_readiness_summary_ok_keeps_diagnostics_visible():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Reporting",
        "Check": "Report readiness",
        "Status": "OK",
        "WhatItSays": "No high-priority blocker was detected.",
        "NextAction": "Continue to report guidance.",
        "DetailLocation": "Report",
    }])

    summary = app.guided_interpretation_readiness_summary_table(plan)
    row = summary.iloc[0]

    assert row[app.t("guided.interpret_col_status")] == app.t("guided.interpret_status_ok")
    assert row[app.t("guided.interpret_col_open_next")] == app.t("guided.interpret_ok_open_next")
    assert row[app.t("guided.interpret_col_do_not_conclude")] == app.t("guided.interpret_ok_do_not")
    assert "diagnostics" in row[app.t("guided.interpret_col_safe_output")].lower()
    assert "Ready to interpret" not in " ".join(summary.astype(str).to_numpy().ravel())


def test_guided_interpretation_readiness_help_table_documents_states():
    table = app.guided_interpretation_readiness_help_table()

    assert list(table.columns) == [
        app.t("help.guided_interpret_col_status"),
        app.t("help.guided_interpret_col_when"),
        app.t("guided.interpret_col_do_not_conclude"),
        app.t("guided.interpret_col_safe_output"),
    ]
    assert len(table) == 5
    assert app.t("guided.interpret_status_pause") in set(table[app.t("help.guided_interpret_col_status")])
    assert app.t("guided.interpret_status_ok") in set(table[app.t("help.guided_interpret_col_status")])
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_action_plan_renders_interpretation_readiness_summary():
    source = inspect.getsource(app._render_guided_action_plan)

    assert "guided_interpretation_readiness_summary_table(plan)" in source
    assert "guided.interpret_heading" in source


def test_guided_goal_routes_cover_common_user_intents():
    routes = app.guided_goal_route_table()

    assert set(routes["GoalId"]) == {
        "check_readiness",
        "inspect_measures",
        "diagnose_problem",
        "prepare_report",
        "export_share",
        "learn_terms",
    }
    expected_columns = {
        "GoalId",
        app.t("guided.goal_col_goal"),
        app.t("guided.goal_col_when"),
        app.t("guided.goal_col_start"),
        app.t("guided.goal_col_next"),
        app.t("guided.goal_col_detail"),
    }
    assert expected_columns.issubset(routes.columns)
    assert routes[app.t("guided.goal_col_goal")].astype(str).str.len().gt(0).all()
    assert routes[app.t("guided.goal_col_next")].astype(str).str.len().gt(0).all()
    assert set(routes[app.t("guided.goal_col_start")]) >= {
        app.t("guided.tab_first_read"),
        app.t("guided.tab_results"),
        app.t("guided.tab_diagnostics"),
        app.t("guided.tab_report_export"),
        app.t("guided.tab_learn"),
    }


def test_guided_figures_section_is_a_direct_visual_surface():
    source = inspect.getsource(app._render_guided_figures_section)

    assert "show_wright_map_section(result, diagnostics)" in source
    assert "show_visuals_section(result, diagnostics" in source
    assert "downloads.figures_skipped_info" in source
    assert "guided.figures_location_markdown" in source
    assert "guided.figures_panel_select_caption" in source


def test_visuals_section_lazy_renders_one_plot_panel():
    source = inspect.getsource(app.show_visuals_section)

    assert "visuals_panel" in source
    assert "st.selectbox" in source
    assert "st.tabs" not in source
    assert 'if selected_visual_panel == "category_probability"' in source
    assert 'if selected_visual_panel == "ecdf"' in source


def test_downloads_privacy_mode_is_defined_before_export_builders():
    source = inspect.getsource(app._render_downloads)

    assert "st.segmented_control" in source
    assert "downloads_panel" in source
    assert "st.tabs" not in source
    checkbox_idx = source.index("public_export_mode = st.checkbox")
    stan_manifest_idx = source.index("stan_reproducibility_package_assets(")
    assert checkbox_idx < stan_manifest_idx


def test_full_mode_main_results_panel_is_lazy_rendered():
    source = inspect.getsource(app.run_facets_mode)

    assert "main_results_panel" in source
    assert "st.selectbox" in source
    assert "main_tabs.panel_select_caption" in source
    assert "tabs = st.tabs([" not in source
    assert 'default="data"' not in source
    assert 'if selected_main_panel == "downloads"' in source


def test_help_section_is_lazy_rendered_by_topic_selector():
    source = inspect.getsource(app.show_help_section)

    assert "help_panel" in source
    assert "st.selectbox" in source
    assert "st.tabs" not in source
    assert "default=visible_labels[0]" not in source
    assert 'if selected_help_label == "Analysis Workflow"' in source
    assert 'if selected_help_label == "Public Beta"' in source


def test_publication_figure_payloads_use_plotly_helpers_with_expected_inputs():
    source = inspect.getsource(app._publication_figure_payloads)

    assert "_make_wright_map_export_figure(person_tbl, facet_tbl, step_tbl)" in source
    assert "build_category_probability_curve_data(result)" in source
    assert "_make_category_probability_curve_figures(curve)" in source
    assert "_make_facet_distribution_export_figure(diagnostics.get(\"measures\", pd.DataFrame()))" in source


def test_run_staleness_tracks_report_scaling_and_visual_outputs():
    source = inspect.getsource(app.run_facets_mode)

    for token in [
        '"_totalscore"',
        '"_omit_unobserved"',
        '"_xtreme"',
        '"_umean"',
        '"_uscale"',
        '"_udecimals"',
        '"_render_interactive_plots"',
        '"_generate_figure_exports"',
        '"_visualization_preferences"',
        "mfrm_new_design_prediction_bundle",
        "mfrm_refit_simulation_bundle",
    ]:
        assert token in source
    assert "FACETS total-score reporting option" in source
    assert "figure export bundle option" in source


def test_guided_goal_step_tables_have_three_actionable_steps():
    expected_goals = {
        "check_readiness",
        "inspect_measures",
        "diagnose_problem",
        "prepare_report",
        "export_share",
        "learn_terms",
    }
    for goal_id in expected_goals:
        steps = app.guided_goal_step_table(goal_id)
        assert len(steps) == 3
        assert list(steps[app.t("guided.goal_step_col_step")]) == [1, 2, 3]
        assert steps[app.t("guided.goal_step_col_open")].astype(str).str.len().gt(0).all()
        assert steps[app.t("guided.goal_step_col_action")].astype(str).str.len().gt(0).all()
        assert steps[app.t("guided.goal_step_col_done")].astype(str).str.len().gt(0).all()


def test_guided_goal_step_table_falls_back_to_readiness_goal():
    fallback = app.guided_goal_step_table("not_a_known_goal")
    readiness = app.guided_goal_step_table("check_readiness")

    assert fallback.equals(readiness)


def test_guided_goal_choice_help_table_routes_uncertain_users_to_goals():
    table = app.guided_goal_choice_help_table()
    goal_col = app.t("help.guided_choice_col_choose_goal")

    assert list(table.columns) == [
        app.t("help.guided_choice_col_question"),
        app.t("help.guided_choice_col_choose_goal"),
        app.t("help.guided_choice_col_first_open"),
        app.t("help.guided_choice_col_hold_off"),
        app.t("help.guided_choice_col_ready_when"),
    ]
    assert len(table) == len(app.GUIDED_GOAL_IDS)
    assert list(table[goal_col]) == [
        app.t(f"guided.goal_{goal_id}") for goal_id in app.GUIDED_GOAL_IDS
    ]
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_reader_profile_help_table_keeps_detail_by_reader_scope():
    table = app.guided_reader_profile_help_table()
    profile_col = app.t("help.guided_reader_col_profile")

    assert list(table.columns) == [
        app.t("help.guided_reader_col_profile"),
        app.t("help.guided_reader_col_read_first"),
        app.t("help.guided_reader_col_keep_detail"),
        app.t("help.guided_reader_col_stop_condition"),
        app.t("help.guided_reader_col_output"),
    ]
    assert list(table[profile_col]) == [
        app.t("help.guided_reader_first_pass_profile"),
        app.t("help.guided_reader_applied_user_profile"),
        app.t("help.guided_reader_analyst_profile"),
        app.t("help.guided_reader_manuscript_author_profile"),
        app.t("help.guided_reader_reviewer_auditor_profile"),
    ]
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_focus_class_help_table_documents_routing_classes():
    table = app.guided_focus_class_help_table()
    class_col = app.t("help.guided_focus_class_col_class")

    assert list(table.columns) == [
        app.t("help.guided_focus_class_col_class"),
        app.t("help.guided_focus_class_col_when"),
        app.t("help.guided_focus_class_col_user_action"),
        app.t("help.guided_focus_class_col_examples"),
    ]
    assert list(table[class_col]) == [
        app.t("help.guided_focus_class_hard_label"),
        app.t("help.guided_focus_class_claim_caveat_label"),
        app.t("help.guided_focus_class_goal_review_label"),
        app.t("help.guided_focus_class_clean_label"),
    ]
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_next_click_help_table_documents_summary_labels():
    table = app.guided_next_click_help_table()
    label_col = app.t("help.guided_next_click_col_label")

    assert list(table.columns) == [
        app.t("help.guided_next_click_col_label"),
        app.t("help.guided_next_click_col_when"),
        app.t("help.guided_next_click_col_open"),
        app.t("help.guided_next_click_col_keep_detail"),
    ]
    assert list(table[label_col]) == [
        app.t("guided.next_click_run_first"),
        app.t("guided.next_click_hard_stop"),
        app.t("guided.next_click_selected_goal"),
        app.t("guided.next_click_claim_caveat"),
        app.t("guided.next_click_goal_route"),
        app.t("guided.next_click_review"),
    ]
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_screen_reading_order_help_table_keeps_main_path_first():
    table = app.guided_screen_reading_order_help_table()
    surface_col = app.t("help.guided_screen_order_col_surface")

    assert list(table.columns) == [
        app.t("help.guided_screen_order_col_order"),
        app.t("help.guided_screen_order_col_surface"),
        app.t("help.guided_screen_order_col_use"),
        app.t("help.guided_screen_order_col_detail_policy"),
    ]
    assert list(table[surface_col]) == [
        app.t("help.guided_screen_order_current_focus_surface"),
        app.t("help.guided_screen_order_next_click_surface"),
        app.t("help.guided_screen_order_goal_locator_surface"),
        app.t("help.guided_screen_order_guardrail_surface"),
        app.t("help.guided_screen_order_supporting_detail_surface"),
        app.t("help.guided_screen_order_claim_boundary_surface"),
    ]
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()
    assert app.t("guided.goal_supporting_detail_expander") in set(table[surface_col])


def test_guided_goal_locator_help_table_lists_every_goal_detail_place():
    table = app.guided_goal_locator_help_table()
    goal_col = app.t("help.guided_goal_locator_col_goal")

    assert list(table.columns) == [
        app.t("help.guided_goal_locator_col_goal"),
        app.t("guided.goal_locator_col_open"),
        app.t("guided.goal_locator_col_place"),
        app.t("guided.goal_locator_col_check"),
        app.t("guided.goal_locator_col_use"),
    ]
    assert len(table) == len(app.GUIDED_GOAL_IDS) * 3
    assert set(table[goal_col]) == {
        app.t(f"guided.goal_{goal_id}") for goal_id in app.GUIDED_GOAL_IDS
    }
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()
    assert "Combined measures" in " ".join(table.astype(str).to_numpy().ravel().tolist())
    assert "evidence matrix" in " ".join(table.astype(str).to_numpy().ravel().tolist())


def test_guided_goal_behavior_help_table_covers_all_goals_and_caveat_modes():
    table = app.guided_goal_behavior_help_table()
    goal_col = app.t("help.guided_goal_behavior_col_goal")
    caveat_col = app.t("help.guided_goal_behavior_col_if_caveat")

    assert len(table) == len(app.GUIDED_GOAL_IDS)
    assert set(table[goal_col]) == {
        app.t(f"guided.goal_{goal_id}") for goal_id in app.GUIDED_GOAL_IDS
    }
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()

    inspect_row = table.loc[table[goal_col].eq(app.t("guided.goal_inspect_measures"))].iloc[0]
    report_row = table.loc[table[goal_col].eq(app.t("guided.goal_prepare_report"))].iloc[0]
    assert inspect_row[caveat_col] == app.t("help.guided_goal_behavior_continue_with_caveat")
    assert report_row[caveat_col] == app.t("help.guided_goal_behavior_main_focus")


def test_guided_claim_boundary_help_table_preserves_reporting_scope():
    table = app.guided_claim_boundary_help_table()
    focus_col = app.t("help.guided_claim_boundary_col_focus_type")

    assert list(table.columns) == [
        app.t("help.guided_claim_boundary_col_focus_type"),
        app.t("help.guided_claim_boundary_col_safe_now"),
        app.t("help.guided_claim_boundary_col_not_yet"),
        app.t("help.guided_claim_boundary_col_confirm"),
    ]
    assert list(table[focus_col]) == [
        app.t("help.guided_focus_class_hard_label"),
        app.t("help.guided_focus_class_claim_caveat_label"),
        app.t("help.guided_focus_class_goal_review_label"),
        app.t("help.guided_focus_class_clean_label"),
    ]
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_reproducibility_guardrail_preserves_professional_checks():
    table = app.guided_reproducibility_guardrail_table()
    check_col = app.t("guided.repro_col_check")
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())

    assert list(table.columns) == [
        app.t("guided.repro_col_check"),
        app.t("guided.repro_col_verify"),
        app.t("guided.repro_col_artifact"),
        app.t("guided.repro_col_do_not"),
    ]
    assert len(table) == 6
    assert app.t("guided.repro_model_setup_check") in set(table[check_col])
    assert app.t("guided.repro_claim_trace_check") in set(table[check_col])
    assert app.t("guided.repro_rerun_check") in set(table[check_col])
    assert "Claim-to-evidence" in joined
    assert "external" in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_report_claim_trace_table_badges_generated_claims():
    result = {
        "opt": SimpleNamespace(success=True, message="converged"),
        "config": {
            "model": "PCM",
            "method": "JMLE",
            "facet_names": ["Rater", "Task"],
            "app_version": "test",
            "n_cat": 4,
        },
        "prep": {
            "n_obs": 20,
            "n_person": 5,
            "rating_min": 0,
            "rating_max": 3,
            "unused_score_categories": [],
        },
    }
    diagnostics = {
        "obs": pd.DataFrame({
            "StdResidual": [0.1, 0.2, -0.1],
            "Expected": [1.0, 2.0, 3.0],
            "Observed": [1.0, 2.0, 3.0],
        }),
        "reliability": pd.DataFrame({"Facet": ["person"], "Reliability": [0.70]}),
        "pca_enabled": False,
        "marginal_fit": {},
    }

    table = app.guided_report_claim_trace_table(result, diagnostics, {})
    item_col = app.t("guided.claim_trace_col_item")
    status_col = app.t("guided.claim_trace_col_status")
    evidence_col = app.t("guided.claim_trace_col_evidence")
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())

    assert list(table.columns) == [
        app.t("guided.claim_trace_col_item"),
        app.t("guided.claim_trace_col_status"),
        app.t("guided.claim_trace_col_evidence"),
        app.t("guided.claim_trace_col_before_use"),
        app.t("guided.claim_trace_col_open"),
    ]
    assert len(table) == 5
    assert app.t("guided.claim_trace_overall_item") in set(table[item_col])
    assert app.t("guided.claim_trace_apa_item") in set(table[item_col])
    assert app.t("guided.claim_trace_matrix_item") in set(table[item_col])
    assert table[status_col].astype(str).str.contains(r"\[CAVEAT\]|\[READY\]|\[REVIEW\]").any()
    assert "apa_report_sentence_audit.csv" in joined
    assert "claim_to_evidence_matrix.csv" in joined
    assert "publication_gate_summary.csv" in " ".join(table[evidence_col].astype(str).tolist())
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_report_claim_trace_help_table_documents_evidence_files():
    table = app.guided_report_claim_trace_help_table()
    file_col = app.t("help.claim_trace_col_file")
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())

    assert list(table.columns) == [
        app.t("help.claim_trace_col_file"),
        app.t("help.claim_trace_col_use"),
        app.t("help.claim_trace_col_checks"),
        app.t("help.claim_trace_col_boundary"),
    ]
    assert len(table) == 5
    assert app.t("help.claim_trace_publication_gate_file") in set(table[file_col])
    assert app.t("help.claim_trace_apa_audit_file") in set(table[file_col])
    assert app.t("help.claim_trace_claim_matrix_file") in set(table[file_col])
    assert "CopyDecision" in joined
    assert "ReviewerQuestion" in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_export_share_preflight_table_maps_share_targets_to_artifacts():
    result = {
        "opt": SimpleNamespace(success=True, message="converged"),
        "config": {
            "model": "PCM",
            "method": "JMLE",
            "facet_names": ["Rater", "Task"],
            "app_version": "test",
            "n_cat": 4,
        },
        "prep": {
            "n_obs": 20,
            "n_person": 5,
            "rating_min": 0,
            "rating_max": 3,
            "unused_score_categories": [],
        },
    }
    diagnostics = {
        "obs": pd.DataFrame({
            "StdResidual": [0.1, 0.2, -0.1],
            "Expected": [1.0, 2.0, 3.0],
            "Observed": [1.0, 2.0, 3.0],
        }),
        "reliability": pd.DataFrame({"Facet": ["person"], "Reliability": [0.70]}),
        "pca_enabled": False,
        "marginal_fit": {},
    }

    table = app.guided_export_share_preflight_table(result, diagnostics, {})
    audience_col = app.t("guided.export_preflight_col_audience")
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())

    assert list(table.columns) == [
        app.t("guided.export_preflight_col_audience"),
        app.t("guided.export_preflight_col_use"),
        app.t("guided.export_preflight_col_required_checks"),
        app.t("guided.export_preflight_col_keep"),
        app.t("guided.export_preflight_col_stop"),
    ]
    assert len(table) == 6
    assert app.t("guided.export_preflight_public_package_audience") in set(table[audience_col])
    assert app.t("guided.export_preflight_manuscript_review_audience") in set(table[audience_col])
    assert app.t("guided.export_preflight_external_validation_audience") in set(table[audience_col])
    assert "export_privacy_manifest.csv" in joined
    assert "claim_to_evidence_matrix.csv" in joined
    assert "figure_manifest.csv" in joined
    assert "external_simulation_reference_inventory.csv" in joined
    assert "Current gate:" in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_export_share_preflight_help_table_documents_share_boundaries():
    table = app.guided_export_share_preflight_help_table()
    audience_col = app.t("help.export_preflight_col_audience")
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())

    assert list(table.columns) == [
        app.t("help.export_preflight_col_audience"),
        app.t("help.export_preflight_col_use"),
        app.t("help.export_preflight_col_checks"),
        app.t("help.export_preflight_col_boundary"),
    ]
    assert len(table) == 6
    assert app.t("help.export_preflight_public_package_audience") in set(table[audience_col])
    assert app.t("help.export_preflight_reviewer_repo_audience") in set(table[audience_col])
    assert "export_privacy_manifest.csv" in joined
    assert "visual_qa_preflight.csv" in joined
    assert "formal de-identification" in joined
    assert "prospective-validity" in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_report_export_renders_reproducibility_guardrail():
    source = inspect.getsource(app._render_guided_report_export_section)

    assert "build_report_ready_summary_panel(" in source
    assert "guided.report_ready_summary_heading" in source
    assert "report_ready_summary_panel.csv" in source
    assert "generate_report_ready_decision_brief(" in source
    assert "generate_report_ready_apa_results_draft(" in source
    apa_draft_source = inspect.getsource(app.generate_report_ready_apa_results_draft)
    assert "build_apa_bridge_revision_plan(" in apa_draft_source
    assert "Bridge-Constrained Use Conditions" in apa_draft_source
    assert "Bridge Revision Plan" in apa_draft_source
    assert "build_role_based_action_memos(" in source
    assert "generate_role_based_action_memos_markdown(" in source
    assert "build_role_based_action_checklist(" in source
    assert "apply_role_based_action_checklist_done_state(" in source
    assert "role_based_action_checklist_progress_summary(" in source
    assert "generate_role_based_action_checklist_markdown(" in source
    assert "build_report_ready_reanalysis_checklist(" in source
    assert "generate_report_ready_reanalysis_checklist_markdown(" in source
    assert "report_ready_decision_brief.md" in source
    assert "apa_results_paragraph_draft.md" in source
    assert "role_based_action_memos.csv" in source
    assert "role_based_action_memos.md" in source
    assert "role_based_action_checklist.csv" in source
    assert "role_based_action_checklist.md" in source
    assert "report_ready_reanalysis_checklist.csv" in source
    assert "report_ready_reanalysis_checklist.md" in source
    assert "build_critical_final_review_panel(" in source
    assert "generate_critical_final_review_markdown(" in source
    assert "critical_final_review_panel.csv" in source
    assert "critical_final_review.md" in source
    assert "guided.critical_review_heading" in source
    assert "build_status_rationale_drilldown(" in source
    assert "build_result_reading_path(" in source
    assert "build_operational_decision_board(" in source
    assert "build_reporting_action_bridge(" in source
    assert "generate_reporting_action_bridge_markdown(" in source
    assert "status_rationale_drilldown.csv" in source
    assert "result_reading_path.csv" in source
    assert "operational_decision_board.csv" in source
    assert "reporting_action_bridge.csv" in source
    assert "reporting_action_bridge.md" in source
    assert "guided.audience_boards_heading" in source
    assert "guided.reporting_bridge_heading" in source
    assert "st.data_editor(" in source
    assert "st.progress(" in source
    assert "guided.role_checklist_open_only_label" in source
    assert "guided_reproducibility_guardrail_table()" in source
    assert "guided.repro_heading" in source
    assert "guided_report_claim_trace_table(" in source
    assert "guided.claim_trace_heading" in source
    assert "guided_export_share_preflight_table(" in source
    assert "guided.export_preflight_heading" in source
    assert "show_report_section(" in source


def test_help_section_renders_reproducibility_guardrail_reference():
    source = inspect.getsource(app.show_help_section)

    assert "guided_reproducibility_guardrail_table()" in source
    assert "help.repro_guardrail_expander" in source
    assert "guided_stan_posterior_reproducibility_help_table()" in source
    assert "help.stan_repro_heading" in source
    assert "guided_uto_claim_boundary_help_table()" in source
    assert "guided_uto_design_audit_field_help_table()" in source
    assert "help.uto_audit_heading" in source
    assert "stan_posterior_handoff_checklist()" in source
    assert "stan_run_manifest_template()" in source
    assert "stan_posterior_reproducibility_handoff_markdown()" in source
    assert "guided_report_claim_trace_help_table()" in source
    assert "help.claim_trace_expander" in source
    assert "guided_export_share_preflight_help_table()" in source
    assert "help.export_preflight_expander" in source


def test_guided_stan_posterior_reproducibility_help_table_keeps_manifest_detail():
    table = app.guided_stan_posterior_reproducibility_help_table()
    checkpoint_col = app.t("help.stan_repro_col_checkpoint")

    assert list(table.columns) == [
        app.t("help.stan_repro_col_checkpoint"),
        app.t("help.stan_repro_col_evidence"),
        app.t("help.stan_repro_col_check"),
        app.t("help.stan_repro_col_claim_boundary"),
    ]
    assert len(table) == 9
    joined = " ".join(table.astype(str).to_numpy().ravel().tolist())
    assert app.t("help.stan_repro_prior_decision_checkpoint") in set(table[checkpoint_col])
    assert app.t("help.stan_repro_uto_design_audit_checkpoint") in set(table[checkpoint_col])
    assert "MFRM_Complete_Stan_Reproducibility_Package.zip" in joined
    assert "mfrm_uto_bayesian_mfrm_design_audit.csv" in joined
    assert "mfrm_stan_prior_decision_log_template.csv" in joined
    assert "stan_run_manifest.json" in joined
    assert "stan_file_sha256" in joined
    assert "sigma_theta_prior_scale" in joined
    assert "adapt_delta" in joined
    assert "Rhat" in joined
    assert "ESS" in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_uto_design_audit_help_tables_define_claim_boundaries():
    claim_table = app.guided_uto_claim_boundary_help_table()
    field_table = app.guided_uto_design_audit_field_help_table()

    assert list(claim_table.columns) == [
        app.t("help.uto_claim_col_audit_row"),
        app.t("help.uto_claim_col_ready_when"),
        app.t("help.uto_claim_col_review_means"),
        app.t("help.uto_claim_col_report_boundary"),
        app.t("help.uto_claim_col_verify"),
    ]
    assert len(claim_table) == 6
    joined_claims = " ".join(claim_table.astype(str).to_numpy().ravel().tolist())
    assert "C_dim" in joined_claims
    assert "S=1" in joined_claims
    assert "arbitrary row order" in joined_claims
    assert "mfrm_uto_bayesian_mfrm_design_audit.csv" not in joined_claims
    assert claim_table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()

    assert list(field_table.columns) == [
        app.t("help.uto_field_col_column"),
        app.t("help.uto_field_col_meaning"),
        app.t("help.uto_field_col_use"),
    ]
    assert len(field_table) == 9
    assert {
        "DesignRole",
        "MappedFacet",
        "Status",
        "ClaimReadiness",
        "AllowedClaimBoundary",
        "UtoReferenceLayer",
    }.issubset(set(field_table[app.t("help.uto_field_col_column")]))
    assert field_table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_stan_posterior_reproducibility_handoff_markdown_is_portable():
    handoff = app.stan_posterior_reproducibility_handoff_markdown()

    assert "# Stan Posterior Reproducibility Handoff" in handoff
    assert "mfrm_model.stan" in handoff
    assert "mfrm_stan_data.json" in handoff
    assert "stan_run_manifest.json" in handoff
    assert "stan_run_manifest.csv" in handoff
    assert "Posterior Viewer" in handoff
    assert "Prior Setting And Sensitivity" in handoff
    assert "sigma_theta_prior_scale" in handoff
    assert "facet_prior_scale" in handoff
    assert "step_prior_scale" in handoff
    assert "Uto-family" in handoff
    assert "/Users/" not in handoff
    assert "C:/Users/" not in handoff


def test_guided_situation_walkthrough_help_table_covers_common_run_states():
    table = app.guided_situation_walkthrough_help_table()
    situation_col = app.t("help.guided_situation_col_situation")

    assert list(table.columns) == [
        app.t("help.guided_situation_col_situation"),
        app.t("help.guided_situation_col_first_open"),
        app.t("help.guided_situation_col_inspect"),
        app.t("help.guided_situation_col_record"),
        app.t("help.guided_situation_col_avoid"),
    ]
    assert len(table) == 8
    assert app.t("help.guided_situation_nonconvergence_situation") in set(table[situation_col])
    assert app.t("help.guided_situation_fit_pca_situation") in set(table[situation_col])
    assert app.t("help.guided_situation_manuscript_export_situation") in set(table[situation_col])
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_guided_simulation_help_table_covers_edge_case_expectations():
    table = app.guided_simulation_help_table()
    focus_col = app.t("help.guided_sim_col_expected_focus")

    assert list(table.columns) == [
        app.t("help.guided_sim_col_situation"),
        app.t("help.guided_sim_col_expected_focus"),
        app.t("help.guided_sim_col_test_reason"),
    ]
    assert len(table) == 6
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()
    assert app.t("help.guided_sim_nonconvergence_expected_focus") in set(table[focus_col])
    assert app.t("help.guided_sim_claim_caveat_expected_focus") in set(table[focus_col])
    assert app.t("help.guided_sim_goal_plus_caveat_expected_focus") in set(table[focus_col])


def test_guided_goal_current_focus_keeps_claim_caveat_while_allowing_goal():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 2,
        "Area": "Fit",
        "Check": "Global residual screen",
        "Status": "Caution",
        "WhatItSays": "Residual evidence needs review.",
        "NextAction": "Open fit details before interpreting measures.",
        "DetailLocation": "Fit Details",
    }])

    focus = app.guided_goal_current_focus_table("inspect_measures", plan)
    goal_row = focus.iloc[0]
    caveat_row = focus.iloc[1]

    assert list(focus.columns) == [
        app.t("guided.goal_focus_col_focus"),
        app.t("guided.goal_focus_col_why"),
        app.t("guided.goal_focus_col_open"),
        app.t("guided.goal_focus_col_action"),
    ]
    assert len(focus) == 2
    assert app.t("guided.goal_focus_selected_goal_prefix") in goal_row[app.t("guided.goal_focus_col_focus")]
    assert goal_row[app.t("guided.goal_focus_col_open")] == app.t("guided.tab_results")
    assert app.t("guided.goal_focus_claim_caveat_prefix") in caveat_row[app.t("guided.goal_focus_col_focus")]
    assert "P1:" in caveat_row[app.t("guided.goal_focus_col_focus")]
    assert app.t("guided.status_caution") in caveat_row[app.t("guided.goal_focus_col_focus")]
    assert app.t("guided.area_fit") in caveat_row[app.t("guided.goal_focus_col_focus")]
    assert "Residual evidence needs review." in caveat_row[app.t("guided.goal_focus_col_why")]
    assert caveat_row[app.t("guided.goal_focus_col_open")] == app.t("guided.detail_location_fit_details")
    assert caveat_row[app.t("guided.goal_focus_col_action")] == "Open fit details before interpreting measures."


def test_guided_goal_current_focus_keeps_hard_flags_above_goal_relevance():
    plan = pd.DataFrame([
        {
            "Priority": "P1",
            "SourceOrder": 1,
            "Area": "Reporting",
            "Check": "Convergence status",
            "Status": "Do not interpret yet",
            "WhatItSays": "The optimizer did not finish cleanly.",
            "NextAction": "Fix estimation before reading measures.",
            "DetailLocation": "Report",
        },
        {
            "Priority": "P2",
            "SourceOrder": 2,
            "Area": "Measures",
            "Check": "Reliability summary",
            "Status": "Review",
            "WhatItSays": "Reliability needs context.",
            "NextAction": "Inspect measures.",
            "DetailLocation": "Report",
        },
    ])

    focus = app.guided_goal_current_focus_table("inspect_measures", plan)
    row = focus.iloc[0]

    assert "Convergence status" in row[app.t("guided.goal_focus_col_focus")]
    assert app.t("guided.status_do_not_interpret") in row[app.t("guided.goal_focus_col_focus")]
    assert row[app.t("guided.goal_focus_col_action")] == "Fix estimation before reading measures."


def test_guided_goal_current_focus_prefers_goal_relevant_review_rows():
    plan = pd.DataFrame([
        {
            "Priority": "P1",
            "SourceOrder": 1,
            "Area": "Reporting",
            "Check": "Report wording review",
            "Status": "Review",
            "WhatItSays": "Report text needs a caveat.",
            "NextAction": "Review claim guidance.",
            "DetailLocation": "Report",
        },
        {
            "Priority": "P2",
            "SourceOrder": 2,
            "Area": "Measures",
            "Check": "Reliability summary",
            "Status": "Review",
            "WhatItSays": "Reliability needs context.",
            "NextAction": "Inspect combined measures.",
            "DetailLocation": "Report",
        },
    ])

    focus = app.guided_goal_current_focus_table("inspect_measures", plan)
    row = focus.iloc[0]

    assert "Reliability summary" in row[app.t("guided.goal_focus_col_focus")]
    assert app.t("guided.area_measures") in row[app.t("guided.goal_focus_col_focus")]
    assert app.t("guided.goal_focus_goal_relevant_why") in row[app.t("guided.goal_focus_col_why")]
    assert row[app.t("guided.goal_focus_col_action")] == "Inspect combined measures."


def test_guided_goal_current_focus_keeps_pca_caveat_below_measure_route():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Assumptions",
        "Check": "Dimensionality",
        "Status": "Caution",
        "WhatItSays": "Residual PCA needs content review.",
        "NextAction": "Open dimensionality before final claims.",
        "DetailLocation": "Dimensionality",
    }])

    focus = app.guided_goal_current_focus_table("inspect_measures", plan)
    goal_row = focus.iloc[0]
    caveat_row = focus.iloc[1]

    assert len(focus) == 2
    assert app.t("guided.goal_focus_selected_goal_prefix") in goal_row[app.t("guided.goal_focus_col_focus")]
    assert goal_row[app.t("guided.goal_focus_col_open")] == app.t("guided.tab_results")
    assert app.t("guided.goal_focus_claim_caveat_prefix") in caveat_row[app.t("guided.goal_focus_col_focus")]
    assert "Dimensionality" in caveat_row[app.t("guided.goal_focus_col_focus")]
    assert caveat_row[app.t("guided.goal_focus_col_action")] == "Open dimensionality before final claims."


def test_guided_goal_current_focus_uses_selected_route_when_clean():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Reporting",
        "Check": "Report readiness",
        "Status": "OK",
        "WhatItSays": "No high-priority blocker was detected.",
        "NextAction": "Continue to report guidance.",
        "DetailLocation": "Report",
    }])

    focus = app.guided_goal_current_focus_table("prepare_report", plan)
    route = app._guided_goal_route_row("prepare_report")
    row = focus.iloc[0]

    assert row[app.t("guided.goal_focus_col_focus")] == app.t("guided.goal_focus_clean_focus")
    assert row[app.t("guided.goal_focus_col_why")] == app.t("guided.goal_focus_prepare_report_clean_why")
    assert row[app.t("guided.goal_focus_col_open")] == route[app.t("guided.goal_col_start")]
    assert row[app.t("guided.goal_focus_col_action")] == route[app.t("guided.goal_col_next")]


def test_guided_goal_current_focus_falls_back_for_unknown_goal():
    focus = app.guided_goal_current_focus_table("not_a_known_goal", pd.DataFrame())
    readiness = app._guided_goal_route_row("check_readiness")
    row = focus.iloc[0]

    assert row[app.t("guided.goal_focus_col_focus")] == app.t("guided.goal_focus_unavailable_focus")
    assert row[app.t("guided.goal_focus_col_open")] == readiness[app.t("guided.goal_col_start")]


def test_guided_next_click_summary_runs_first_when_no_action_plan():
    summary = app.guided_next_click_summary_table("inspect_measures", pd.DataFrame())
    row = summary.iloc[0]

    assert list(summary.columns) == [
        app.t("guided.next_click_col_next"),
        app.t("guided.next_click_col_open"),
        app.t("guided.next_click_col_do"),
        app.t("guided.next_click_col_keep_visible"),
    ]
    assert row[app.t("guided.next_click_col_next")] == app.t("guided.next_click_run_first")
    assert row[app.t("guided.next_click_col_open")] == app.t("guided.tab_start")
    assert row[app.t("guided.next_click_col_do")] == app.t("guided.goal_focus_unavailable_why")
    assert row[app.t("guided.next_click_col_keep_visible")] == app.t("guided.next_click_keep_run_first")


def test_guided_next_click_summary_prioritizes_hard_stop():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Reporting",
        "Check": "Convergence status",
        "Status": "Do not interpret yet",
        "WhatItSays": "The optimizer did not finish cleanly.",
        "NextAction": "Fix estimation before reading measures.",
        "DetailLocation": "Report",
    }])

    row = app.guided_next_click_summary_table("inspect_measures", plan).iloc[0]

    assert row[app.t("guided.next_click_col_next")] == app.t("guided.next_click_hard_stop")
    assert row[app.t("guided.next_click_col_open")] == app.t("guided.detail_location_report")
    assert row[app.t("guided.next_click_col_do")] == "Fix estimation before reading measures."
    assert row[app.t("guided.next_click_col_keep_visible")] == app.t("guided.next_click_keep_hard_stop")


def test_guided_next_click_summary_keeps_caveat_for_measure_route():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 2,
        "Area": "Fit",
        "Check": "Global residual screen",
        "Status": "Caution",
        "WhatItSays": "Residual evidence needs review.",
        "NextAction": "Open fit details before interpreting measures.",
        "DetailLocation": "Fit Details",
    }])

    row = app.guided_next_click_summary_table("inspect_measures", plan).iloc[0]
    keep_visible = row[app.t("guided.next_click_col_keep_visible")]

    assert row[app.t("guided.next_click_col_next")] == app.t("guided.next_click_selected_goal")
    assert row[app.t("guided.next_click_col_open")] == app.t("guided.tab_results")
    assert row[app.t("guided.next_click_col_do")] == app.t("guided.goal_inspect_measures_next")
    assert app.t("guided.next_click_keep_prefix") in keep_visible
    assert app.t("guided.goal_focus_claim_caveat_prefix") in keep_visible
    assert "Global residual screen" in keep_visible


def test_guided_next_click_summary_opens_claim_caveat_before_report():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Assumptions",
        "Check": "Dimensionality",
        "Status": "Caution",
        "WhatItSays": "Residual PCA needs content review.",
        "NextAction": "Open dimensionality before final claims.",
        "DetailLocation": "Dimensionality",
    }])

    row = app.guided_next_click_summary_table("prepare_report", plan).iloc[0]

    assert row[app.t("guided.next_click_col_next")] == app.t("guided.next_click_claim_caveat")
    assert row[app.t("guided.next_click_col_open")] == app.t("guided.detail_location_dimensionality")
    assert row[app.t("guided.next_click_col_do")] == "Open dimensionality before final claims."
    assert row[app.t("guided.next_click_col_keep_visible")] == app.t("guided.next_click_keep_claim_caveat")


def test_guided_next_click_summary_uses_goal_route_when_clean():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Reporting",
        "Check": "Report readiness",
        "Status": "OK",
        "WhatItSays": "No high-priority blocker was detected.",
        "NextAction": "Continue to report guidance.",
        "DetailLocation": "Report",
    }])
    route = app._guided_goal_route_row("prepare_report")

    row = app.guided_next_click_summary_table("prepare_report", plan).iloc[0]

    assert row[app.t("guided.next_click_col_next")] == app.t("guided.next_click_goal_route")
    assert row[app.t("guided.next_click_col_open")] == route[app.t("guided.goal_col_start")]
    assert row[app.t("guided.next_click_col_do")] == route[app.t("guided.goal_col_next")]
    assert row[app.t("guided.next_click_col_keep_visible")] == app.t("guided.next_click_keep_evidence")


def test_guided_goal_guardrail_summary_keeps_not_yet_decision_visible():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Reporting",
        "Check": "Convergence status",
        "Status": "Do not interpret yet",
        "WhatItSays": "The optimizer did not finish cleanly.",
        "NextAction": "Fix estimation before reading measures.",
        "DetailLocation": "Report",
    }])

    guardrail = app.guided_goal_guardrail_summary_table("inspect_measures", plan)
    row = guardrail.iloc[0]

    assert list(guardrail.columns) == [
        app.t("guided.goal_guardrail_col_do_not_decide"),
        app.t("guided.goal_guardrail_col_keep_visible"),
        app.t("guided.goal_guardrail_col_expected_output"),
    ]
    assert row[app.t("guided.goal_guardrail_col_do_not_decide")] == app.t(
        "guided.goal_brief_inspect_measures_do_not_decide"
    )
    assert row[app.t("guided.goal_guardrail_col_keep_visible")] == app.t(
        "guided.next_click_keep_hard_stop"
    )
    assert row[app.t("guided.goal_guardrail_col_expected_output")] == app.t(
        "guided.goal_brief_inspect_measures_output"
    )


def test_guided_goal_guardrail_summary_falls_back_to_readiness_goal():
    fallback = app.guided_goal_guardrail_summary_table("not_a_known_goal", pd.DataFrame())
    readiness = app.guided_goal_guardrail_summary_table("check_readiness", pd.DataFrame())

    assert fallback.equals(readiness)


def test_guided_section_id_for_target_maps_common_detail_labels():
    assert app.guided_section_id_for_target(app.t("guided.tab_results")) == "results"
    assert app.guided_section_id_for_target("Fit Details / Dimensionality") == "diagnostics"
    assert app.guided_section_id_for_target("Downloads and manuscript binder") == "report_export"
    assert app.guided_section_id_for_target("Glossary and interpretation guide") == "learn"


def test_guided_action_hub_cards_keep_primary_detail_and_boundary_actions():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Fit",
        "Check": "Global residual screen",
        "Status": "Caution",
        "WhatItSays": "Residual evidence needs review.",
        "NextAction": "Open fit details before interpreting measures.",
        "DetailLocation": "Fit Details",
    }])

    cards = app.guided_action_hub_cards("inspect_measures", plan)

    assert len(cards) == 3
    assert set(cards["ActionId"]) == {"primary", "detail", "boundary"}
    assert {
        app.t("guided.action_hub_col_role"),
        app.t("guided.action_hub_col_open"),
        app.t("guided.action_hub_col_action"),
        app.t("guided.action_hub_col_detail"),
        "SectionId",
        app.t("guided.action_hub_col_button"),
    }.issubset(cards.columns)
    assert cards.loc[cards["ActionId"].eq("primary"), "SectionId"].iloc[0] == "results"
    assert cards.loc[cards["ActionId"].eq("detail"), "SectionId"].iloc[0] == "results"
    assert cards.loc[cards["ActionId"].eq("boundary"), "SectionId"].iloc[0] == "report_export"
    joined = " ".join(cards.astype(str).to_numpy().ravel())
    assert app.t("guided.action_hub_primary_role") in joined
    assert app.t("guided.action_hub_detail_role") in joined
    assert app.t("guided.action_hub_boundary_role") in joined


def test_guided_goal_router_renders_action_hub_with_section_buttons():
    source = inspect.getsource(app._render_guided_goal_router)

    assert "guided_action_hub_cards" in source
    assert "guided_essential_section" in source
    assert "st.button" in source
    assert "st.rerun()" in source
    assert "action_hub_detail_expander" in source


def test_guided_progress_checklist_is_checkbox_ready_and_caveat_sensitive():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Fit",
        "Check": "Global residual screen",
        "Status": "Caution",
        "WhatItSays": "Residual evidence needs review.",
        "NextAction": "Open fit details before interpreting measures.",
        "DetailLocation": "Fit Details",
    }])

    checklist = app.guided_progress_checklist("inspect_measures", plan)
    done_col = app.t("guided.progress_col_done")
    checkpoint_col = app.t("guided.progress_col_checkpoint")
    state_col = app.t("guided.progress_col_state")

    assert list(checklist.columns) == [
        done_col,
        app.t("guided.progress_col_checkpoint"),
        app.t("guided.progress_col_state"),
        app.t("guided.progress_col_open"),
        app.t("guided.progress_col_action"),
        app.t("guided.progress_col_completion"),
        app.t("guided.progress_col_boundary"),
    ]
    assert len(checklist) == 6
    assert checklist[done_col].map(type).eq(bool).all()
    assert app.t("guided.progress_state_user_check") in checklist[state_col].tolist()
    assert bool(
        checklist.loc[
            checklist[checkpoint_col].eq(app.t("guided.progress_goal_checkpoint")),
            done_col,
        ].iloc[0]
    )
    assert not bool(
        checklist.loc[
            checklist[checkpoint_col].eq(app.t("guided.progress_claim_checkpoint")),
            done_col,
        ].iloc[0]
    )
    joined = " ".join(checklist.astype(str).to_numpy().ravel())
    assert app.t("guided.progress_claim_boundary") in joined
    assert "presumed ability" not in joined.lower()


def test_guided_progress_checklist_marks_clean_claim_boundary_ready():
    plan = pd.DataFrame([{
        "Priority": "P1",
        "SourceOrder": 1,
        "Area": "Reporting",
        "Check": "Report readiness",
        "Status": "OK",
        "WhatItSays": "No high-priority blocker was detected.",
        "NextAction": "Continue to report guidance.",
        "DetailLocation": "Report",
    }])

    checklist = app.guided_progress_checklist("prepare_report", plan)
    done_col = app.t("guided.progress_col_done")
    checkpoint_col = app.t("guided.progress_col_checkpoint")
    state_col = app.t("guided.progress_col_state")

    claim_row = checklist.loc[
        checklist[checkpoint_col].eq(app.t("guided.progress_claim_checkpoint"))
    ].iloc[0]
    focus_row = checklist.loc[
        checklist[checkpoint_col].eq(app.t("guided.progress_focus_checkpoint"))
    ].iloc[0]

    assert bool(claim_row[done_col])
    assert bool(focus_row[done_col])
    assert claim_row[state_col] == app.t("guided.progress_state_ready")


def test_guided_goal_router_renders_progress_checklist_with_checkbox_column():
    source = inspect.getsource(app._render_guided_goal_router)

    assert "guided_progress_checklist" in source
    assert "st.data_editor" in source
    assert "CheckboxColumn" in source
    assert "mfrm_guided_progress_checklist.csv" in source


def test_guided_goal_decision_brief_has_prompt_scope_output_and_caveat():
    expected_goals = {
        "check_readiness",
        "inspect_measures",
        "diagnose_problem",
        "prepare_report",
        "export_share",
        "learn_terms",
    }
    expected_columns = [
        app.t("guided.goal_brief_col_prompt"),
        app.t("guided.goal_brief_col_decide"),
        app.t("guided.goal_brief_col_do_not_decide"),
        app.t("guided.goal_brief_col_output"),
        app.t("guided.goal_brief_col_caveat"),
    ]
    for goal_id in expected_goals:
        brief = app.guided_goal_decision_brief_table(goal_id)
        assert list(brief.columns) == expected_columns
        assert len(brief) == 1
        assert brief.iloc[0].astype(str).str.len().gt(0).all()


def test_guided_goal_decision_brief_falls_back_to_readiness_goal():
    fallback = app.guided_goal_decision_brief_table("not_a_known_goal")
    readiness = app.guided_goal_decision_brief_table("check_readiness")

    assert fallback.equals(readiness)


def test_guided_goal_evidence_checklists_have_guardrails():
    expected_goals = {
        "check_readiness",
        "inspect_measures",
        "diagnose_problem",
        "prepare_report",
        "export_share",
        "learn_terms",
    }
    expected_columns = [
        app.t("guided.goal_evidence_col_evidence"),
        app.t("guided.goal_evidence_col_reason"),
        app.t("guided.goal_evidence_col_caveat"),
    ]
    for goal_id in expected_goals:
        checklist = app.guided_goal_evidence_checklist(goal_id)
        assert list(checklist.columns) == expected_columns
        assert len(checklist) == 3
        assert checklist.iloc[0].astype(str).str.len().gt(0).all()
        assert checklist[app.t("guided.goal_evidence_col_caveat")].astype(str).str.len().gt(0).all()


def test_guided_goal_evidence_checklist_falls_back_to_readiness_goal():
    fallback = app.guided_goal_evidence_checklist("not_a_known_goal")
    readiness = app.guided_goal_evidence_checklist("check_readiness")

    assert fallback.equals(readiness)


def test_guided_goal_detail_locator_tables_point_to_specific_places():
    expected_goals = {
        "check_readiness",
        "inspect_measures",
        "diagnose_problem",
        "prepare_report",
        "export_share",
        "learn_terms",
    }
    expected_columns = [
        app.t("guided.goal_locator_col_open"),
        app.t("guided.goal_locator_col_place"),
        app.t("guided.goal_locator_col_check"),
        app.t("guided.goal_locator_col_use"),
    ]
    for goal_id in expected_goals:
        locator = app.guided_goal_detail_locator_table(goal_id)
        assert list(locator.columns) == expected_columns
        assert len(locator) == 3
        assert locator.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()

    measures = app.guided_goal_detail_locator_table("inspect_measures")
    report = app.guided_goal_detail_locator_table("prepare_report")
    assert "Combined measures" in " ".join(measures.iloc[0].astype(str).tolist())
    assert "evidence matrix" in " ".join(report.astype(str).to_numpy().ravel().tolist())


def test_guided_goal_detail_locator_falls_back_to_readiness_goal():
    fallback = app.guided_goal_detail_locator_table("not_a_known_goal")
    readiness = app.guided_goal_detail_locator_table("check_readiness")

    assert fallback.equals(readiness)
