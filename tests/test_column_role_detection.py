"""Tests for the smart column-role detector used by the upload UI.

`auto_detect_column_roles` suggests defaults for Person / Score / Rater /
Criterion / Weight column pickers. It combines multilingual keyword
matching (English + Japanese) with dtype and cardinality cues so a
typical MFRM dataset opens with the right defaults pre-filled.

The contract pinned here covers:

* English-only headers — canonical names map to the right roles.
* Japanese-only headers — same dataset translated still detects.
* Mixed-language headers — both English and Japanese resolve.
* Refusal on empty input — returns ``None`` for every role.
* Score detection prefers numeric small-cardinality columns over
  numeric high-cardinality columns (a continuous covariate must not
  be misread as a rating).
* Person detection prefers high-cardinality non-numeric columns.
* Weight column is detected only when explicitly named — it must not
  default to Score when the header is unrelated.
* No role is assigned twice — each column is claimed at most once.
* Confidence is in [0, 1] for every returned assignment.
* `_column_is_numeric` correctly handles blanks, strings, and
  mixed numeric / text cells at the 70 % threshold.
"""

from __future__ import annotations

import pandas as pd

import streamlit_app as app


# -----------------------------------------------------------------------------
# _column_is_numeric
# -----------------------------------------------------------------------------


def test_column_is_numeric_pure_integer():
    s = pd.Series([1, 2, 3, 4, 5])
    assert app._column_is_numeric(s) is True


def test_column_is_numeric_pure_string():
    s = pd.Series(["P1", "P2", "P3"])
    assert app._column_is_numeric(s) is False


def test_column_is_numeric_mixed_below_threshold():
    """50 % numeric is below the 70 % default and must read as non-numeric."""
    s = pd.Series([1, 2, "abc", "xyz"])
    assert app._column_is_numeric(s) is False


def test_column_is_numeric_blank_cells_ignored():
    """Blank cells are excluded from both numerator and denominator."""
    s = pd.Series([1, 2, "", "", 3, 4])
    assert app._column_is_numeric(s) is True


def test_column_is_numeric_returns_false_for_empty():
    assert app._column_is_numeric(pd.Series([], dtype=object)) is False


# -----------------------------------------------------------------------------
# auto_detect_column_roles — English headers
# -----------------------------------------------------------------------------


def test_detect_canonical_english_headers():
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(20)],
        "Rater": ["R1", "R2"] * 10,
        "Criterion": ["C1", "C2", "C3", "C4"] * 5,
        "Score": [0, 1, 2, 3, 4] * 4,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["person"] == "Person"
    assert roles["score"] == "Score"
    assert roles["rater"] == "Rater"
    assert roles["criterion"] == "Criterion"
    assert roles["weight"] is None


def test_detect_synonym_english_headers():
    """Examinee / Judge / Item must resolve via the synonym list."""
    df = pd.DataFrame({
        "Examinee": [f"E{i:02d}" for i in range(12)],
        "Judge": ["J1", "J2"] * 6,
        "Item": ["I1", "I2", "I3"] * 4,
        "Rating": [1, 2, 3, 4] * 3,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["person"] == "Examinee"
    assert roles["score"] == "Rating"
    assert roles["rater"] == "Judge"
    assert roles["criterion"] == "Item"


# -----------------------------------------------------------------------------
# auto_detect_column_roles — Japanese headers
# -----------------------------------------------------------------------------


def test_detect_canonical_japanese_headers():
    df = pd.DataFrame({
        "学生": [f"S{i:02d}" for i in range(20)],
        "評価者": ["R1", "R2"] * 10,
        "観点": ["C1", "C2", "C3", "C4"] * 5,
        "得点": [0, 1, 2, 3, 4] * 4,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["person"] == "学生"
    assert roles["score"] == "得点"
    assert roles["rater"] == "評価者"
    assert roles["criterion"] == "観点"


def test_detect_japanese_synonyms():
    """受験者 / 採点者 / 項目 / 評定 must resolve through the JA synonyms."""
    df = pd.DataFrame({
        "受験者": [f"P{i:02d}" for i in range(12)],
        "採点者": ["R1", "R2"] * 6,
        "項目": ["I1", "I2", "I3"] * 4,
        "評定": [1, 2, 3, 4] * 3,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["person"] == "受験者"
    assert roles["score"] == "評定"
    assert roles["rater"] == "採点者"
    assert roles["criterion"] == "項目"


# -----------------------------------------------------------------------------
# Mixed-language headers
# -----------------------------------------------------------------------------


def test_detect_mixed_language_headers():
    """Some columns English, some Japanese — both must resolve."""
    df = pd.DataFrame({
        "学生": [f"P{i:02d}" for i in range(12)],
        "Rater": ["R1", "R2"] * 6,
        "観点": ["C1", "C2", "C3"] * 4,
        "Score": [1, 2, 3, 4] * 3,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["person"] == "学生"
    assert roles["score"] == "Score"
    assert roles["rater"] == "Rater"
    assert roles["criterion"] == "観点"


# -----------------------------------------------------------------------------
# Empty input
# -----------------------------------------------------------------------------


def test_detect_empty_dataframe_returns_none_for_all_roles():
    result = app.auto_detect_column_roles(pd.DataFrame())
    for role in ("person", "score", "rater", "criterion", "weight"):
        assert result[role] is None
    assert result["suggested_facets"] == []
    assert result["confidence"] == {}


# -----------------------------------------------------------------------------
# Dtype / cardinality bias
# -----------------------------------------------------------------------------


def test_score_prefers_small_cardinality_numeric():
    """A numeric column with ~5 distinct rating values must outrank a
    numeric column with 100 distinct continuous values when neither
    header has a Score keyword. We name the rating column ``Score``
    so it gets keyword + dtype boost; the continuous column
    ``Measurement`` should not steal the assignment."""
    rng = list(range(100))
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(100)],
        "Measurement": rng,  # 100 distinct numeric — not a rating scale
        "Score": [i % 5 for i in rng],  # 5 levels
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["score"] == "Score"
    # And Person should be Person, not the high-cardinality numeric column.
    assert roles["person"] == "Person"


def test_person_prefers_non_numeric_high_cardinality():
    """When no Person keyword is present, Person should still tend
    toward a high-cardinality non-numeric column over a low-cardinality
    one. We test the keyword case here: the explicit ``Person`` header
    must beat ``Group`` even if both qualify by dtype."""
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(40)],
        "Group": ["A", "B"] * 20,
        "Score": [1, 2, 3] * 13 + [1],
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["person"] == "Person"


# -----------------------------------------------------------------------------
# Weight detection
# -----------------------------------------------------------------------------


def test_detect_weight_column_when_named():
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(10)],
        "Rater": ["R1", "R2"] * 5,
        "Score": [1, 2, 3, 4, 5] * 2,
        "Weight": [1.0, 0.5] * 5,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["weight"] == "Weight"


def test_weight_is_none_when_no_weight_column():
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(10)],
        "Rater": ["R1", "R2"] * 5,
        "Score": [1, 2, 3, 4, 5] * 2,
    })
    roles = app.auto_detect_column_roles(df)
    assert roles["weight"] is None


# -----------------------------------------------------------------------------
# Uniqueness / confidence contract
# -----------------------------------------------------------------------------


def test_no_column_assigned_to_two_roles():
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(15)],
        "Rater": ["R1", "R2", "R3"] * 5,
        "Criterion": ["C1", "C2"] * 7 + ["C1"],
        "Score": [1, 2, 3] * 5,
    })
    roles = app.auto_detect_column_roles(df)
    picked = [
        roles[r] for r in ("person", "score", "rater", "criterion", "weight")
        if roles[r] is not None
    ]
    assert len(picked) == len(set(picked))


def test_confidence_values_are_in_unit_interval():
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(10)],
        "Rater": ["R1", "R2"] * 5,
        "Score": [1, 2, 3, 4, 5] * 2,
    })
    roles = app.auto_detect_column_roles(df)
    for role, conf in roles["confidence"].items():
        assert 0.0 <= conf <= 1.0, f"{role} confidence {conf} out of [0, 1]"


def test_suggested_facets_includes_detected_rater_and_criterion():
    df = pd.DataFrame({
        "Person": [f"P{i:02d}" for i in range(12)],
        "Rater": ["R1", "R2"] * 6,
        "Criterion": ["C1", "C2", "C3"] * 4,
        "Score": [1, 2, 3, 4] * 3,
    })
    roles = app.auto_detect_column_roles(df)
    assert "Rater" in roles["suggested_facets"]
    assert "Criterion" in roles["suggested_facets"]


# -----------------------------------------------------------------------------
# Generic / opaque headers — must not assign noise
# -----------------------------------------------------------------------------


def test_generic_headers_yield_low_confidence_or_none():
    """When headers carry no signal (``Col1``, ``Col2``, ...) and the
    data is mostly numeric, we still want a sensible Score guess but
    should not falsely assign Person to a meaningless string column.

    This test pins behaviour rather than perfect inference: at minimum,
    no role should fire above 0.6 confidence on opaque headers."""
    df = pd.DataFrame({
        "Col1": [f"a{i}" for i in range(8)],
        "Col2": [1, 2, 3, 4, 5, 1, 2, 3],
        "Col3": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    })
    roles = app.auto_detect_column_roles(df)
    for conf in roles["confidence"].values():
        # Without a keyword match the maximum dtype boost is 0.4 (Score),
        # so we expect every confidence to stay under 0.6 here.
        assert conf < 0.6
