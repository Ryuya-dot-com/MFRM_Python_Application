from __future__ import annotations

import numpy as np
import pandas as pd

import streamlit_app as app


def test_response_data_audit_marks_likelihood_rows_and_reasons():
    df = pd.DataFrame({
        "Person": ["P1", "P2", "P3", "P4", "P5"],
        "Rater": ["R1", "R1", "R2", None, "R2"],
        "Task": ["T1", "T1", "T2", "T2", "T2"],
        "Score": [1, np.nan, 3, 4, 5],
        "Weight": [1, 1, 0, 1, 1],
    })

    audit = app.build_response_data_audit(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
        weight_col="Weight",
    )

    rows = audit["rows"]
    assert rows["IncludedInLikelihood"].tolist() == [True, False, False, False, True]
    assert "missing score" in rows.loc[1, "ExclusionReason"]
    assert "non-positive weight" in rows.loc[2, "ExclusionReason"]
    assert "missing Rater" in rows.loc[3, "ExclusionReason"]

    summary = audit["summary"]
    included = summary.loc[summary["Item"] == "Rows included in likelihood", "Count"].iloc[0]
    excluded = summary.loc[summary["Item"] == "Rows excluded before likelihood", "Count"].iloc[0]
    assert included == 2
    assert excluded == 3


def test_prepare_mfrm_data_keeps_audit_and_excludes_blank_facet_rows():
    df = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2"],
        "Rater": ["R1", " ", "R1", "R2"],
        "Task": ["T1", "T1", "T2", "T2"],
        "Score": [0, 1, 1, 2],
    })

    prep = app.prepare_mfrm_data(
        df,
        person_col="Person",
        facet_cols=["Rater", "Task"],
        score_col="Score",
    )

    assert len(prep["data"]) == 3
    assert "row_audit" in prep
    excluded = prep["excluded_rows"]
    assert len(excluded) == 1
    assert int(excluded.iloc[0]["InputRow"]) == 2
    assert "missing Rater" in excluded.iloc[0]["ExclusionReason"]


def test_readiness_report_warns_when_rows_will_be_excluded():
    df = pd.DataFrame({
        "Person": ["P1", "P2", "P3", "P4"] * 30,
        "Rater": ["R1", "R2", "R1", "R2"] * 30,
        "Task": ["T1", "T1", "T2", "T2"] * 30,
        "Score": ([0, 1, 2, np.nan] * 30),
    })

    report = app.build_readiness_report(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater", "Task"],
    )

    likelihood = next(c for c in report["checks"] if c["name"] == "likelihood_rows")
    assert likelihood["severity"] == "warning"
    assert "excluded before likelihood" in likelihood["headline"]

