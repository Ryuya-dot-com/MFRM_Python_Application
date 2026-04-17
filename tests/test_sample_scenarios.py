"""Contract tests for the four built-in sample-data scenarios.

The scenario registry is the user's first touch with the app, so
drift here (shape change, citation unresolved, non-determinism) is
the kind of thing that breaks screenshots and onboarding docs. These
tests sit in pytest so they run independently of the inline
self-test suite.
"""

from __future__ import annotations

import pandas as pd
import pytest

import streamlit_app as app


ALL_KEYS = list(app.SAMPLE_DATA_SCENARIOS.keys())


def test_registry_has_four_scenarios():
    assert len(app.SAMPLE_DATA_SCENARIOS) >= 4


def test_default_key_is_in_registry():
    assert app.DEFAULT_SAMPLE_SCENARIO_KEY in app.SAMPLE_DATA_SCENARIOS


@pytest.mark.parametrize("key", ALL_KEYS)
def test_scenario_metadata_shape(key: str):
    meta = app.SAMPLE_DATA_SCENARIOS[key]
    for field in ("label", "short", "description", "dimensions",
                  "n_obs", "build_params", "citations"):
        assert field in meta, f"{key!r} missing {field!r}"
    dims = meta["dimensions"]
    for dkey in ("persons", "raters", "tasks", "criteria", "n_cat"):
        assert dkey in dims, f"{key!r} dimensions missing {dkey!r}"
    assert meta["n_obs"] == (
        dims["persons"] * dims["raters"] * dims["tasks"] * dims["criteria"]
    ), f"{key!r}: n_obs does not match dimension product"


@pytest.mark.parametrize("key", ALL_KEYS)
def test_scenario_generates_expected_shape(key: str):
    meta = app.SAMPLE_DATA_SCENARIOS[key]
    df = app.sample_mfrm_data_by_key(key, seed=42)
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (meta["n_obs"], 5)


@pytest.mark.parametrize("key", ALL_KEYS)
def test_scenario_score_bounds(key: str):
    meta = app.SAMPLE_DATA_SCENARIOS[key]
    df = app.sample_mfrm_data_by_key(key, seed=42)
    score_col = df.columns[-1]
    scores = pd.to_numeric(df[score_col], errors="coerce")
    assert scores.notna().all(), f"{key!r}: non-numeric scores found"
    assert scores.min() >= 0, f"{key!r}: negative score"
    assert scores.max() <= meta["dimensions"]["n_cat"] - 1, (
        f"{key!r}: score exceeds n_cat - 1"
    )


@pytest.mark.parametrize("key", ALL_KEYS)
def test_scenario_facets_have_multiple_levels(key: str):
    df = app.sample_mfrm_data_by_key(key, seed=42)
    # Facet columns should have ≥ 2 unique levels — except pro-forma
    # singleton Rater / Scorer facets used by single-scorer binary
    # testlet designs. MFRM simply centers those to 0; they contribute
    # zero variance and are not a bug.
    for col in df.columns[:-1]:
        min_levels = 1 if col in ("Rater", "Scorer") else 2
        assert df[col].nunique() >= min_levels, (
            f"{key!r}: column {col!r} has < {min_levels} level(s)"
        )
        assert df[col].notna().all(), f"{key!r}: column {col!r} has NaN"


@pytest.mark.parametrize("key", ALL_KEYS)
def test_scenario_citations_resolve(key: str):
    meta = app.SAMPLE_DATA_SCENARIOS[key]
    for cite in meta["citations"]:
        mapped = app._CITATION_TO_KEY.get(cite)
        assert mapped is not None, (
            f"{key!r}: citation {cite!r} not in _CITATION_TO_KEY"
        )
        assert mapped in app._APA_REFERENCE_LIBRARY, (
            f"{key!r}: citation {cite!r} → {mapped!r} not in library"
        )


@pytest.mark.parametrize("key", ALL_KEYS)
def test_scenario_is_deterministic(key: str):
    df1 = app.sample_mfrm_data_by_key(key, seed=99)
    df2 = app.sample_mfrm_data_by_key(key, seed=99)
    assert df1.equals(df2), f"{key!r} non-deterministic for fixed seed"


def test_sample_mfrm_data_legacy_alias_preserved():
    """Backward-compat: the v0.1 alias still returns the 960-row writing-essay."""
    legacy = app.sample_mfrm_data(seed=20240101)
    assert legacy.shape == (960, 5)
    assert list(legacy.columns) == ["Person", "Rater", "Task", "Criterion", "Score"]


def test_unknown_scenario_falls_back_to_default():
    """Robustness: an unknown key should not raise; it should return the default."""
    df_default = app.sample_mfrm_data_by_key(app.DEFAULT_SAMPLE_SCENARIO_KEY, seed=7)
    df_unknown = app.sample_mfrm_data_by_key("__does_not_exist__", seed=7)
    assert df_unknown.equals(df_default)


def test_large_writing_has_enough_persons_for_pca():
    """The PCA-ready scenario must have ≥ 100 persons for stable residual PCA."""
    df = app.sample_mfrm_data_by_key("large_writing_pca", seed=1)
    assert df["Person"].nunique() >= 100
