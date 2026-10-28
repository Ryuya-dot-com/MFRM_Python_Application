"""Tests for `build_readiness_report` — the pre-estimation data-quality gate.

The readiness panel is the first structured feedback users get on
their data. If it mis-classifies a bad upload as "[OK] Ready" or a
clean upload as "[ISSUE] Issues", users waste time. These tests pin the
severity logic on representative fixtures.
"""

from __future__ import annotations

import inspect

import pandas as pd
import numpy as np
import pytest

import streamlit_app as app


def _default_scenario_df() -> pd.DataFrame:
    return app.sample_mfrm_data_by_key("writing_essay", seed=0)


def _find_check(report: dict, name: str) -> dict | None:
    for c in report.get("checks", []):
        if c["name"] == name:
            return c
    return None


def test_empty_data_triggers_issue():
    report = app.build_readiness_report(
        data=pd.DataFrame(),
        person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert report["overall"] == "issue"
    assert report["n_issues"] >= 1


def test_input_readiness_next_action_loads_data_when_empty():
    report = app.build_readiness_report(
        data=pd.DataFrame(),
        person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )

    summary = app.build_input_readiness_next_action_table(report)
    row = summary.iloc[0]

    assert list(summary.columns) == [
        app.t("input_readiness.col_status"),
        app.t("input_readiness.col_next_action"),
        app.t("input_readiness.col_open"),
        app.t("input_readiness.col_reason"),
    ]
    assert row[app.t("input_readiness.col_status")] == app.t("input_readiness.status_issue")
    assert row[app.t("input_readiness.col_next_action")] == app.t("input_readiness.action_load_data")
    assert row[app.t("input_readiness.col_open")] == app.t("input_readiness.open_data_source")


def test_input_readiness_panel_copy_uses_i18n_keys():
    source = inspect.getsource(app.render_readiness_panel)

    assert "input_readiness.panel_ok_template" in source
    assert "input_readiness.panel_warning_template" in source
    assert "input_readiness.panel_issue_template" in source
    assert "input_readiness.detail_expander_ok" in source
    assert "input_readiness.detail_expander_action_needed" in source
    assert "Ready to run." not in source
    assert "Proceed with caution." not in source
    assert "Not ready to run." not in source
    assert "Show data quality detail" not in source
    assert "Show all checks" not in source


def test_input_readiness_panel_templates_are_formattable():
    ok_text = app.t(
        "input_readiness.panel_ok_template",
        icon="🟢",
        status="OK",
        n_checks=3,
    )
    warning_text = app.t(
        "input_readiness.panel_warning_template",
        icon="🟡",
        status="CAUTION",
        n_warnings=2,
    )
    issue_text = app.t(
        "input_readiness.panel_issue_template",
        icon="🔴",
        status="ISSUE",
        n_issues=1,
    )

    assert "input" in ok_text.lower()
    assert "Ready to interpret" not in ok_text
    assert "{n_checks}" not in ok_text
    assert "input" in warning_text.lower()
    assert "2" in warning_text
    assert "{n_warnings}" not in warning_text
    assert "input" in issue_text.lower()
    assert "1" in issue_text
    assert "{n_issues}" not in issue_text


def test_clean_sample_scenario_is_ok_or_warning():
    """The built-in Writing Essay scenario is clean — must not produce [ISSUE]."""
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    assert report["overall"] in {"ok", "warning"}
    assert report["n_issues"] == 0


def test_input_readiness_next_action_points_to_facet_mapping():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater"],
    )

    row = app.build_input_readiness_next_action_table(report).iloc[0]

    assert row[app.t("input_readiness.col_status")] == app.t("input_readiness.status_issue")
    assert row[app.t("input_readiness.col_next_action")] == app.t("input_readiness.action_select_two_facets")
    assert row[app.t("input_readiness.col_open")] == app.t("input_readiness.open_column_mapping")


def test_input_readiness_next_action_points_to_score_column_for_no_variation():
    persons = [f"P{i:02d}" for i in range(1, 31)]
    df = pd.DataFrame({
        "Person": persons * 4,
        "Rater": ["R1", "R2", "R3", "R4"] * 30,
        "Task": ["T1", "T2"] * 60,
        "Score": [3] * 120,
    })
    report = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater", "Task"],
    )

    row = app.build_input_readiness_next_action_table(report).iloc[0]

    assert row[app.t("input_readiness.col_status")] == app.t("input_readiness.status_issue")
    assert row[app.t("input_readiness.col_next_action")] == app.t("input_readiness.action_fix_score_column")
    assert row[app.t("input_readiness.col_open")] == app.t("input_readiness.open_column_mapping")


def test_input_readiness_detail_table_preserves_all_checks_and_detail():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )

    table = app.build_input_readiness_detail_table(report["checks"])

    assert list(table.columns) == [
        app.t("input_readiness.col_status"),
        app.t("input_readiness.col_check"),
        app.t("input_readiness.col_meaning"),
        app.t("input_readiness.col_impact"),
        app.t("input_readiness.col_next_action"),
        app.t("input_readiness.col_open"),
        app.t("input_readiness.col_technical_detail"),
    ]
    assert len(table) == len(report["checks"])
    assert table[app.t("input_readiness.col_check")].astype(str).str.len().gt(0).all()
    assert table[app.t("input_readiness.col_meaning")].astype(str).str.len().gt(0).all()
    assert table[app.t("input_readiness.col_impact")].astype(str).str.len().gt(0).all()
    assert table[app.t("input_readiness.col_technical_detail")].astype(str).str.len().gt(0).all()


def test_input_readiness_detail_table_orders_actionable_checks_first():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater"],
    )

    table = app.build_input_readiness_detail_table(report["checks"])
    first = table.iloc[0]

    assert first[app.t("input_readiness.col_status")].startswith("🔴 [ISSUE]")
    assert first[app.t("input_readiness.col_next_action")] == app.t(
        "input_readiness.action_select_two_facets"
    )
    assert first[app.t("input_readiness.col_open")] == app.t("input_readiness.open_column_mapping")


def test_input_readiness_detail_table_explains_score_check():
    persons = [f"P{i:02d}" for i in range(1, 31)]
    df = pd.DataFrame({
        "Person": persons * 4,
        "Rater": ["R1", "R2", "R3", "R4"] * 30,
        "Task": ["T1", "T2"] * 60,
        "Score": [3] * 120,
    })
    report = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater", "Task"],
    )

    table = app.build_input_readiness_detail_table(report["checks"])
    score_row = table[
        table[app.t("input_readiness.col_next_action")]
        == app.t("input_readiness.action_fix_score_column")
    ].iloc[0]

    assert score_row[app.t("input_readiness.col_meaning")] == app.t(
        "input_readiness.guidance_score_meaning"
    )
    assert score_row[app.t("input_readiness.col_impact")] == app.t(
        "input_readiness.guidance_score_impact"
    )
    assert score_row[app.t("input_readiness.col_technical_detail")]


def test_input_readiness_detail_table_has_fallback_guidance_for_unknown_check():
    checks = [{
        "name": "custom_quality_gate",
        "severity": "warning",
        "headline": "Custom quality gate triggered",
        "detail": "A project-specific readiness check needs review.",
    }]

    table = app.build_input_readiness_detail_table(checks)
    row = table.iloc[0]

    assert row[app.t("input_readiness.col_meaning")] == app.t(
        "input_readiness.guidance_flagged_check_meaning"
    )
    assert row[app.t("input_readiness.col_impact")] == app.t(
        "input_readiness.guidance_flagged_check_impact"
    )
    assert row[app.t("input_readiness.col_next_action")] == app.t(
        "input_readiness.action_review_flagged_check"
    )


def test_missing_person_column_is_flagged():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df, person_col="NotThere", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    check = _find_check(report, "n_persons")
    assert check is not None
    assert check["severity"] == "issue"


def test_single_facet_is_flagged():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater"],
    )
    check = _find_check(report, "n_facets")
    assert check is not None
    assert check["severity"] == "issue"


def test_column_role_overlap_is_flagged():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Person",  # overlap
        facet_cols=["Rater", "Task"],
    )
    check = _find_check(report, "column_overlap")
    assert check is not None
    assert check["severity"] == "issue"


def test_small_dataset_shows_warning():
    # 25 observations, below the 100-obs warning threshold
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(1, 6)] * 5,
        "Rater": (["R1"] * 5 + ["R2"] * 5 + ["R3"] * 5 + ["R1"] * 5 + ["R2"] * 5),
        "Score": [1, 2, 3, 4, 5] * 5,
    })
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater"],
    )
    # The single-facet check will raise an issue, but n_observations
    # should flag the sample size.
    check = _find_check(report, "n_observations")
    assert check is not None
    assert check["severity"] == "warning"


def _small_crossed_data():
    return pd.DataFrame([
        {"Person": f"P{p}", "Rater": f"R{r}", "Task": f"T{t}",
         "Score": (p + r + t) % 4}
        for p in range(4) for r in range(2) for t in range(2)
    ])


@pytest.mark.parametrize("case", ["small", "mostly_missing", "one_combination", "constant_persons"])
def test_usable_input_cautions_do_not_block(case):
    data = _small_crossed_data()
    if case == "mostly_missing":
        data = pd.concat([data, data.assign(Score=np.nan), data.assign(Score=np.nan)], ignore_index=True)
    elif case == "one_combination":
        data = data.groupby("Person", sort=False).head(1)
    elif case == "constant_persons":
        data = pd.concat([data.assign(Person=lambda d: d.Person + str(i)) for i in range(3)], ignore_index=True)
        data["Score"] = data.Person.str[1].astype(int) % 3
    report = app.build_readiness_report(data=data, person_col="Person", score_col="Score", facet_cols=["Rater", "Task"])
    assert report["overall"] == "warning"
    assert report["n_issues"] == 0
    # The engine's existing preprocessing accepts the same input; no fit needed.
    prep = app.prepare_mfrm_data(data, "Person", ["Rater", "Task"], "Score")
    assert len(prep["data"]) > 0


@pytest.mark.parametrize("case", ["fractional", "infinite_score", "infinite_weight", "no_rows", "retained_constant", "near_integer_constant", "overlap", "duplicate_names", "missing_column"])
def test_input_blocks_agree_with_existing_preprocessor(case):
    data = _small_crossed_data().assign(Weight=1.0)
    score_col = "Score"
    if case == "fractional":
        data["Score"] = data.Score.astype(float)
        data.loc[0, "Score"] = .5
    elif case == "infinite_score":
        data["Score"] = data.Score.astype(float)
        data.loc[0, "Score"] = np.inf
    elif case == "infinite_weight":
        data.loc[0, "Weight"] = np.inf
    elif case == "no_rows":
        data["Weight"] = 0
    elif case == "retained_constant":
        data.loc[data.Score.ne(1), "Weight"] = 0
    elif case == "near_integer_constant":
        data["Score"] = 1.0
        data.loc[0, "Score"] = 1.000000001
    elif case == "overlap":
        score_col = "Weight"
    elif case == "duplicate_names":
        data = pd.concat([data, data[["Score"]]], axis=1)
    elif case == "missing_column":
        data = data.drop(columns="Task")
    report = app.build_readiness_report(data=data, person_col="Person", score_col=score_col,
                                        facet_cols=["Rater", "Task"], weight_col="Weight")
    assert report["overall"] == "issue"
    with pytest.raises(ValueError):
        app.prepare_mfrm_data(data, "Person", ["Rater", "Task"], score_col, weight_col="Weight")


def test_person_count_review_does_not_ask_to_remap_a_valid_column():
    report = app.build_readiness_report(data=_small_crossed_data(), person_col="Person", score_col="Score",
                                        facet_cols=["Rater", "Task"])
    check = _find_check(report, "n_persons")
    assert check["severity"] == "warning"
    summary = app.build_input_readiness_next_action_table({"overall": "warning", "checks": [check]})
    assert summary.iloc[0][app.t("input_readiness.col_next_action")] == app.t("input_readiness.action_review_sample_size")


def test_report_includes_outlier_findings():
    """v0.2.9-beta: readiness report must integrate outlier detection."""
    df = _default_scenario_df()
    # Inject a zero-variance person
    inj = pd.DataFrame({
        "Person": ["P_CONST"] * 12,
        "Rater": ["R1", "R2", "R3", "R4"] * 3,
        "Task": ["Task1", "Task2"] * 6,
        "Criterion": ["Content"] * 12,
        "Score": [3] * 12,
    })
    df = pd.concat([df, inj], ignore_index=True)
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    # Zero-variance check should appear somewhere in the report.
    names = {c["name"] for c in report["checks"]}
    assert "zero_variance_persons" in names


def test_report_overall_is_max_severity():
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    # Overall must be the highest severity of any individual check
    severities = [c["severity"] for c in report["checks"]]
    priority = {"ok": 0, "warning": 1, "issue": 2}
    expected = max(priority[s] for s in severities) if severities else 0
    actual = priority.get(report["overall"], -1)
    assert actual == expected


def test_report_counts_match_checks():
    df = _default_scenario_df()
    # Inject enough anomalies to get at least one warning
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    n_issues = sum(1 for c in report["checks"] if c["severity"] == "issue")
    n_warnings = sum(1 for c in report["checks"] if c["severity"] == "warning")
    assert report["n_issues"] == n_issues
    assert report["n_warnings"] == n_warnings


def test_input_readiness_before_run_help_table_documents_pre_run_checks():
    table = app.input_readiness_before_run_help_table()

    assert list(table.columns) == [
        app.t("help.input_readiness_col_check"),
        app.t("help.input_readiness_col_why"),
        app.t("help.input_readiness_col_fix"),
        app.t("help.input_readiness_col_blocking"),
    ]
    assert len(table) == 6
    assert app.t("help.input_readiness_score_column_check") in set(
        table[app.t("help.input_readiness_col_check")]
    )
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_input_readiness_check_reference_help_table_reuses_panel_guidance():
    table = app.input_readiness_check_reference_help_table()

    assert list(table.columns) == [
        app.t("help.input_readiness_reference_col_check"),
        app.t("input_readiness.col_meaning"),
        app.t("input_readiness.col_impact"),
        app.t("input_readiness.col_next_action"),
        app.t("input_readiness.col_open"),
    ]
    assert len(table) >= 12
    assert app.t("input_readiness.guidance_score_meaning") in set(
        table[app.t("input_readiness.col_meaning")]
    )
    assert app.t("input_readiness.guidance_score_saturation_impact") in set(
        table[app.t("input_readiness.col_impact")]
    )
    assert app.t("input_readiness.action_review_flagged_check") in set(
        table[app.t("input_readiness.col_next_action")]
    )
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_input_readiness_ok_boundary_help_table_documents_interpretation_limit():
    table = app.input_readiness_ok_boundary_help_table()

    assert list(table.columns) == [
        app.t("help.input_readiness_ok_col_stage"),
        app.t("help.input_readiness_ok_col_means"),
        app.t("help.input_readiness_ok_col_not_guaranteed"),
        app.t("help.input_readiness_ok_col_next_check"),
    ]
    assert len(table) == 5
    joined = "\n".join(table.astype(str).to_numpy().ravel())
    assert "First Read" in joined
    assert "Convergence" in joined
    assert "model fit" in joined
    assert "Ready to interpret" not in joined
    assert table.astype(str).apply(lambda column: column.str.len().gt(0).all()).all()


def test_help_section_renders_input_readiness_reference_table():
    source = inspect.getsource(app.show_help_section)

    assert "input_readiness_ok_boundary_help_table()" in source
    assert "help.input_readiness_ok_expander" in source
    assert "input_readiness_check_reference_help_table()" in source
    assert "help.input_readiness_reference_expander" in source
