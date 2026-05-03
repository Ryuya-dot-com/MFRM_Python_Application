"""Tests for v0.2.13-beta binary / testlet support.

Verifies that the Reading-testlet binary scenario:
  1. Generates a 0/1 long-format DataFrame with the expected schema.
  2. Runs through the JMLE estimator to a converged solution.
  3. All 5 citations resolve to entries in `_APA_REFERENCE_LIBRARY`.
  4. The Help tab's "MFRM vs testlet models" section is present.

These are the kinds of tests I should have written BEFORE shipping
the scenario — they are the cheapest guard against the class of
regression that keeps getting flagged.
"""

from __future__ import annotations

import pandas as pd

import streamlit_app as app


def test_reading_testlet_binary_schema():
    df = app.sample_mfrm_data_by_key("reading_testlet_binary", seed=7)
    assert list(df.columns) == ["Person", "Scorer", "Text", "Item", "Score"]
    assert df.shape == (2400, 5)
    # Score must be strictly binary
    uniq = sorted(pd.to_numeric(df["Score"]).unique())
    assert uniq == [0, 1], f"expected [0, 1] scores, got {uniq}"
    # Design dimensions match the registry spec
    meta = app.SAMPLE_DATA_SCENARIOS["reading_testlet_binary"]
    dims = meta["dimensions"]
    assert df["Person"].nunique() == dims["persons"] == 100
    assert df["Scorer"].nunique() == dims["raters"] == 1  # single scorer OK
    assert df["Text"].nunique() == dims["tasks"] == 6
    assert df["Item"].nunique() == dims["criteria"] == 4


def test_reading_testlet_binary_is_deterministic():
    df1 = app.sample_mfrm_data_by_key("reading_testlet_binary", seed=99)
    df2 = app.sample_mfrm_data_by_key("reading_testlet_binary", seed=99)
    assert df1.equals(df2)


def test_reading_testlet_binary_passes_mfrm_estimator():
    """End-to-end: binary 2,400-row dataset must converge in JMLE."""
    df = app.sample_mfrm_data_by_key("reading_testlet_binary", seed=0)
    core = app.load_core_namespace()
    res = core["mfrm_estimate"](
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Scorer", "Text", "Item"],
        model="RSM", method="JMLE",
        rating_min=0, rating_max=1,
        maxit=150, reltol=1e-5,
    )
    summary = res.get("summary")
    assert summary is not None and not summary.empty
    assert bool(summary.iloc[0].get("Converged")) is True
    # n_cat = 2 means Rasch (1PL) with a single Andrich threshold
    assert int(res["config"].get("n_cat")) == 2
    steps = res.get("params", {}).get("steps")
    assert steps is not None and len(steps) == 1


def test_all_testlet_citations_resolve():
    meta = app.SAMPLE_DATA_SCENARIOS["reading_testlet_binary"]
    for cite in meta["citations"]:
        key = app._CITATION_TO_KEY.get(cite)
        assert key is not None, f"{cite!r} not in _CITATION_TO_KEY"
        assert key in app._APA_REFERENCE_LIBRARY, (
            f"{cite!r} → {key!r} not in _APA_REFERENCE_LIBRARY"
        )
        entry = app._APA_REFERENCE_LIBRARY[key]
        # Every reference should include a 4-digit year
        import re
        assert re.search(r"\b(19|20)\d{2}\b", entry), (
            f"{key!r} reference body lacks a year: {entry!r}"
        )


def test_new_v0_2_13_references_exist():
    """Pin the 5 testlet references added in v0.2.13-beta."""
    for key in [
        "Wainer_Kiely_1987",
        "Bradlow_Wainer_Wang_1999",
        "Wang_Bradlow_Wainer_2002",
        "Rijmen_2010",
        "DeMars_2006",
    ]:
        assert key in app._APA_REFERENCE_LIBRARY, f"missing {key!r}"


def test_binary_score_range_works_in_generator():
    """Generator must produce 0/1 when tau has a single threshold."""
    params = app._params_reading_testlet_binary()
    assert len(params["tau"]) == 1, "binary should have 1 threshold (n_cat=2)"
    df = app._generate_mfrm_rsm_from_params(params, seed=1)
    assert set(df["Score"].unique()) == {0, 1}


def test_readme_help_section_contains_testlet_content():
    """The Help tab → Model Capability section must explain MFRM vs testlet.

    The content lives in `locales/en.json` (key
    ``help.model_capability_body``) since the i18n refactor moved help
    markdown out of the Python source. This test pins a handful of
    keywords there so a future refactor cannot silently drop the
    comparison table from the canonical English source.
    """
    import json
    from pathlib import Path
    en_path = Path("locales/en.json")
    src = en_path.read_text(encoding="utf-8")
    body = json.loads(src).get("help", {}).get("model_capability_body", "")
    must_contain = [
        "Testlet-format Data",
        "Bradlow, Wainer, and Wang",
        "local-independence",
        "TESTLET_RI",
        "TESTLET_BIFACTOR",
        "Wainer & Kiely",
    ]
    for phrase in must_contain:
        assert phrase in body, f"Help doc missing {phrase!r}"
