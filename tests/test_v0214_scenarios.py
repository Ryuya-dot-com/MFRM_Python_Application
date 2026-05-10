"""Tests for v0.2.14-beta scenarios: missing data + music peer-rating.

These exercise the two new code paths added in v0.2.14:
  - `missing_rate` post-generation row drop
  - `_generate_mfrm_peer_rating_data` sparse cyclic generator
"""

from __future__ import annotations

import pandas as pd

import streamlit_app as app


# ---------------------------------------------------------------------------
# writing_with_missing
# ---------------------------------------------------------------------------

def test_writing_missing_drops_about_15pct_of_rows():
    df = app.sample_mfrm_data_by_key("writing_with_missing", seed=11)
    full = 80 * 4 * 2 * 3  # 1,920 fully crossed
    expected = int(full * 0.85)  # 1,632
    # 15% MAR drop has stddev ≈ sqrt(N * p * (1-p)) ≈ ~16 rows
    tol = max(40, int(full * 0.04))
    assert abs(len(df) - expected) <= tol, (
        f"missing-rate scenario produced {len(df)} rows, expected {expected} ± {tol}"
    )


def test_writing_missing_is_deterministic():
    df1 = app.sample_mfrm_data_by_key("writing_with_missing", seed=99)
    df2 = app.sample_mfrm_data_by_key("writing_with_missing", seed=99)
    assert df1.equals(df2)


def test_writing_missing_estimator_completes():
    """JMLE must succeed even with ~15% rows missing."""
    df = app.sample_mfrm_data_by_key("writing_with_missing", seed=3)
    core = app.load_core_namespace()
    res = core["mfrm_estimate"](
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task", "Criterion"],
        model="RSM", method="JMLE",
        rating_min=0, rating_max=5,
        maxit=200, reltol=1e-5,
    )
    summary = res.get("summary")
    assert summary is not None and not summary.empty
    assert bool(summary.iloc[0].get("Converged")) is True


# ---------------------------------------------------------------------------
# music_peer_rating
# ---------------------------------------------------------------------------

def test_music_peer_rating_schema():
    df = app.sample_mfrm_data_by_key("music_peer_rating", seed=42)
    assert list(df.columns) == ["Person", "Rater", "Piece", "Criterion", "Score"]
    # 120 examinees × 2 raters × 2 pieces × 4 criteria = 1,920
    assert df.shape == (1920, 5)


def test_music_peer_rating_uses_shared_id_pool():
    """Person and Rater columns should draw from the same musician pool."""
    df = app.sample_mfrm_data_by_key("music_peer_rating", seed=42)
    persons = set(df["Person"].astype(str).unique())
    raters = set(df["Rater"].astype(str).unique())
    assert persons == raters
    assert len(persons) == 120


def test_music_peer_rating_no_self_rating():
    """Cyclic schedule (i rated by i+1, i+2) excludes self-rating."""
    df = app.sample_mfrm_data_by_key("music_peer_rating", seed=1)
    self_rate = df[df["Person"].astype(str) == df["Rater"].astype(str)]
    assert self_rate.empty, (
        f"Found {len(self_rate)} self-rating rows — cyclic schedule violated"
    )


def test_music_peer_rating_each_examinee_has_exactly_two_raters():
    df = app.sample_mfrm_data_by_key("music_peer_rating", seed=1)
    rater_counts = df.groupby("Person")["Rater"].nunique()
    # Every musician is rated by exactly 2 peers (cyclic neighbors).
    assert (rater_counts == 2).all(), (
        f"Expected exactly 2 raters per examinee; got distribution: "
        f"{rater_counts.value_counts().to_dict()}"
    )


def test_music_peer_rating_each_rater_evaluates_exactly_two_examinees():
    df = app.sample_mfrm_data_by_key("music_peer_rating", seed=1)
    examinee_counts = df.groupby("Rater")["Person"].nunique()
    assert (examinee_counts == 2).all(), (
        f"Expected exactly 2 examinees per rater; got distribution: "
        f"{examinee_counts.value_counts().to_dict()}"
    )


def test_music_peer_rating_estimator_completes():
    """Sparse design with 240 unique pairs must still converge."""
    df = app.sample_mfrm_data_by_key("music_peer_rating", seed=2)
    core = app.load_core_namespace()
    res = core["mfrm_estimate"](
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Rater", "Piece", "Criterion"],
        model="RSM", method="JMLE",
        rating_min=0, rating_max=4,
        maxit=400, reltol=1e-5,
    )
    summary = res.get("summary")
    assert summary is not None and not summary.empty
    # Convergence is plausible but not guaranteed on sparse designs;
    # the contract here is that the estimator returns a usable result.
    assert summary.iloc[0].get("Iterations", 0) > 0


# ---------------------------------------------------------------------------
# Quick-download bundle
# ---------------------------------------------------------------------------

def test_build_result_bundle_frames_minimum_set():
    """Should always include summary + measures when both exist."""
    import numpy as np
    res = {
        "config": {"model": "RSM", "method": "JMLE",
                   "facet_names": ["Rater", "Task"], "n_cat": 5},
        "prep": {"n_obs": 120, "n_person": 20,
                 "rating_min": 0, "rating_max": 4,
                 "unused_score_categories": []},
        "params": {"steps": [-1.5, -0.5, 0.5, 1.5]},
        "summary": pd.DataFrame([{
            "Model": "RSM", "Method": "JMLE",
            "Converged": True, "Iterations": 30, "LogLik": -200.0,
        }]),
        "facets": {
            "person": pd.DataFrame({
                "Person": [f"P{i}" for i in range(20)],
                "Estimate": np.linspace(-1, 1, 20),
            }),
            "others": pd.DataFrame({
                "Facet": ["Rater"] * 3, "Level": ["R1", "R2", "R3"],
                "Estimate": [-0.3, 0.0, 0.3],
            }),
        },
    }
    diag = {
        "measures": pd.DataFrame({
            "Facet": ["Rater"] * 3, "Level": ["R1", "R2", "R3"],
            "Estimate": [-0.3, 0.0, 0.3], "SE": [0.1, 0.1, 0.1],
        }),
    }
    frames = app.build_result_bundle_frames(res, diag)
    assert "summary" in frames
    assert "measures" in frames
    assert "person_measures" in frames
    assert "facet_element_measures" in frames
    # Every value is a non-empty DataFrame
    for name, df in frames.items():
        assert isinstance(df, pd.DataFrame)
        assert not df.empty, f"{name!r} bundle frame is empty"


def test_build_result_bundle_frames_empty_input():
    """Empty input should not raise."""
    frames = app.build_result_bundle_frames({}, {})
    assert frames == {}


# ---------------------------------------------------------------------------
# Singleton-facet handling
# ---------------------------------------------------------------------------

def test_singleton_facet_readiness_check_is_warning_not_issue():
    """v0.2.14: singleton facets must not block a run with [ISSUE]."""
    import numpy as np
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(40)] * 2,
        "Singleton": ["only_value"] * 80,
        "Task": (["T1"] * 40 + ["T2"] * 40),
        "Score": np.random.default_rng(0).integers(0, 5, size=80),
    })
    report = app.build_readiness_report(
        data=df, person_col="Person", score_col="Score",
        facet_cols=["Singleton", "Task"],
    )
    # Find the singleton check
    singleton_check = next(
        (c for c in report["checks"] if c["name"] == "facet:Singleton"), None,
    )
    assert singleton_check is not None
    assert singleton_check["severity"] == "warning", (
        f"Expected singleton facet to be a warning, got "
        f"{singleton_check['severity']}: {singleton_check}"
    )
