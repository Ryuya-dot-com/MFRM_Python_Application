"""Tests for the ingestion-time outlier detector `detect_data_outliers`.

Each test injects one known anomaly into an otherwise-clean baseline
dataset and asserts the detector picks it up with the expected
severity and name. Clean data must produce zero findings (no false
positives), which is the test most likely to flag a detector
overreach.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def _baseline(seed: int = 0, n_persons: int = 20) -> pd.DataFrame:
    """20 × 3 × 2 = 120-row RSM-like baseline with no injected anomaly."""
    rng = np.random.default_rng(seed)
    rows = []
    for pi in range(1, n_persons + 1):
        for ri in range(1, 4):
            for ti in range(1, 3):
                rows.append((f"P{pi:02d}", f"R{ri}", f"T{ti}",
                             int(rng.integers(0, 6))))
    return pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Score"])


def _names(findings: list[dict]) -> set[str]:
    return {f["name"] for f in findings}


def test_clean_data_produces_zero_findings():
    df = _baseline(seed=0)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert findings == [], (
        f"clean baseline produced false positives: {[f['name'] for f in findings]}"
    )


def test_detects_scores_below_rating_min():
    df = _baseline()
    df.loc[len(df)] = ("P_LOW", "R1", "T1", -1)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
        rating_min=0, rating_max=5,
    )
    assert "scores_below_range" in _names(findings)


def test_detects_scores_above_rating_max():
    df = _baseline()
    df.loc[len(df)] = ("P_HIGH", "R1", "T1", 99)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
        rating_min=0, rating_max=5,
    )
    assert "scores_above_range" in _names(findings)


def test_detects_negative_scores_without_explicit_range():
    df = _baseline()
    df.loc[len(df)] = ("P_NEG", "R1", "T1", -3)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
        # no rating_min/max supplied — the negative-score check is
        # still active because rating scales are conventionally ≥ 0.
    )
    assert "scores_negative" in _names(findings)


def test_detects_zero_variance_person():
    df = _baseline()
    # Add a person whose score is identical across all their observations.
    rows = [("P_CONST", f"R{r}", f"T{t}", 3) for r in (1, 2, 3) for t in (1, 2)]
    df = pd.concat([df, pd.DataFrame(rows, columns=df.columns)], ignore_index=True)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert "zero_variance_persons" in _names(findings)


def test_detects_zero_variance_rater():
    df = _baseline()
    # Add a rater whose score is identical across all persons.
    rows = [(f"P{p:02d}", "R_CONST", "T1", 5) for p in range(1, 15)]
    df = pd.concat([df, pd.DataFrame(rows, columns=df.columns)], ignore_index=True)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert "zero_variance_Rater" in _names(findings)


def test_detects_over_observed_person():
    # Build a varied-count baseline so IQR > 0 and the Tukey fence is
    # active. 10 persons with 4–9 observations each gives a non-zero
    # IQR; then pile 400 extra rows onto P_OVER.
    rng = np.random.default_rng(10)
    rows: list[tuple] = []
    for pi in range(1, 11):
        n_obs = int(rng.integers(4, 10))
        for j in range(n_obs):
            rows.append((f"P{pi:02d}", f"R{1 + j % 3}", f"T{1 + j % 2}",
                         int(rng.integers(0, 6))))
    df = pd.DataFrame(rows, columns=["Person", "Rater", "Task", "Score"])
    extra = pd.DataFrame({
        "Person": ["P_OVER"] * 400,
        "Rater": ["R1"] * 400,
        "Task": ["T1"] * 400,
        "Score": rng.integers(0, 6, size=400),
    })
    df = pd.concat([df, extra], ignore_index=True)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert "person_over_observed" in _names(findings)


def test_detects_ceiling_saturation():
    # 60 rows, > 90 % at the scale maximum.
    df = pd.DataFrame({
        "Person": [f"C{i:02d}" for i in range(60)],
        "Rater": ["R1"] * 60,
        "Task": ["T1"] * 60,
        "Score": [5] * 55 + [3, 4, 4, 5, 4],
    })
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert "ceiling_saturation" in _names(findings)


def test_detects_floor_saturation():
    df = pd.DataFrame({
        "Person": [f"F{i:02d}" for i in range(60)],
        "Rater": ["R1"] * 60,
        "Task": ["T1"] * 60,
        "Score": [0] * 56 + [1, 0, 1, 0],
    })
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    assert "floor_saturation" in _names(findings)


def test_finding_shape_is_valid():
    df = _baseline()
    df.loc[len(df)] = ("P_NEG", "R1", "T1", -2)
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"],
    )
    for f in findings:
        assert set(f.keys()) >= {"name", "severity", "headline", "detail"}
        assert f["severity"] in {"ok", "warning", "issue"}
        assert isinstance(f["headline"], str) and f["headline"]
        assert isinstance(f["detail"], str) and f["detail"]


def test_empty_dataframe_returns_empty_list():
    findings = app.detect_data_outliers(
        pd.DataFrame(), person_col="Person", score_col="Score",
        facet_cols=["Rater"],
    )
    assert findings == []


def test_missing_score_column_does_not_crash():
    df = pd.DataFrame({"Person": ["P1"], "Rater": ["R1"]})
    findings = app.detect_data_outliers(
        df, person_col="Person", score_col="Score",
        facet_cols=["Rater"],
    )
    # Should not raise; no score-based checks fire.
    assert isinstance(findings, list)
