"""Tests for `build_readiness_report` — the pre-estimation data-quality gate.

The readiness panel is the first structured feedback users get on
their data. If it mis-classifies a bad upload as "🟢 Ready" or a
clean upload as "🔴 Issues", users waste time. These tests pin the
severity logic on representative fixtures.
"""

from __future__ import annotations

import pandas as pd

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


def test_clean_sample_scenario_is_ok_or_warning():
    """The built-in Writing Essay scenario is clean — must not produce 🔴."""
    df = _default_scenario_df()
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
    )
    assert report["overall"] in {"ok", "warning"}
    assert report["n_issues"] == 0


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
    assert check["severity"] in {"warning", "issue"}


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
