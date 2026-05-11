"""Tests for the wide-format Excel upload helper.

`detect_wide_format_columns` heuristically classifies an uploaded table
as wide or long; `apply_wide_to_long_pivot` melts a wide-format
DataFrame into the canonical long format the MFRM pipeline expects.

The contract pinned here covers:

* Detection on canonical wide / long fixtures (one numeric vs many).
* Pivot output shape: ``n_id_rows * n_score_cols`` long rows.
* Pivot preserves id columns verbatim.
* Score cells are coerced to numeric; non-numeric / blank cells are
  dropped so the likelihood pipeline never sees zeros for missing.
* Refusal on overlapping id / score columns, empty score_cols, or
  name collisions with the new facet column.
* End-to-end: a wide-format Excel-style fixture pivoted to long must
  fit successfully via ``mfrm_estimate``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Detection
# -----------------------------------------------------------------------------


def test_detect_classifies_long_format_table_as_long():
    """A typical long-format table (Person, Rater, Criterion, Score) has
    exactly one numeric column and should be flagged as long."""
    df = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2"],
        "Rater": ["R1", "R2", "R1", "R2"],
        "Criterion": ["C1", "C1", "C1", "C1"],
        "Score": [3, 2, 4, 3],
    })
    result = app.detect_wide_format_columns(df)
    assert result["looks_wide"] is False


def test_detect_classifies_wide_format_table_as_wide():
    """A wide-format table (Person, Rater, C1, C2, C3) with three numeric
    score columns must be flagged as wide."""
    df = pd.DataFrame({
        "Person": ["P1", "P1", "P2", "P2"],
        "Rater": ["R1", "R2", "R1", "R2"],
        "C1": [3, 2, 4, 3],
        "C2": [4, 3, 3, 2],
        "C3": [2, 4, 3, 3],
    })
    result = app.detect_wide_format_columns(df)
    assert result["looks_wide"] is True
    assert set(result["probable_id_cols"]) >= {"Person", "Rater"}
    assert set(result["probable_score_cols"]) == {"C1", "C2", "C3"}


def test_detect_returns_empty_for_empty_dataframe():
    result = app.detect_wide_format_columns(pd.DataFrame())
    assert result["looks_wide"] is False
    assert "empty" in result["reason"].lower()


# -----------------------------------------------------------------------------
# Pivot output shape
# -----------------------------------------------------------------------------


def test_pivot_produces_expected_row_count():
    """For an n_id_rows x n_score_cols wide table with no missing
    cells, the long output has exactly ``n_id_rows * n_score_cols``
    rows."""
    df = pd.DataFrame({
        "Person": ["P1", "P2", "P3"],
        "Rater": ["R1", "R1", "R1"],
        "C1": [3, 2, 4],
        "C2": [4, 3, 3],
        "C3": [2, 4, 3],
        "C4": [3, 3, 2],
    })
    pivoted = app.apply_wide_to_long_pivot(
        df, id_cols=["Person", "Rater"],
        score_cols=["C1", "C2", "C3", "C4"],
        new_facet_name="Criterion", score_col_name="Score",
    )
    assert len(pivoted) == 3 * 4
    assert set(pivoted.columns) == {"Person", "Rater", "Criterion", "Score"}


def test_pivot_preserves_id_values_per_score_column():
    """Each ID row must appear once per score column in the output."""
    df = pd.DataFrame({
        "Person": ["P1", "P2"],
        "Rater": ["R1", "R2"],
        "C1": [3, 2],
        "C2": [4, 3],
    })
    pivoted = app.apply_wide_to_long_pivot(
        df, id_cols=["Person", "Rater"], score_cols=["C1", "C2"],
        new_facet_name="Criterion", score_col_name="Score",
    )
    # P1 + R1 appears twice (once for C1, once for C2).
    p1_rows = pivoted[(pivoted["Person"] == "P1") & (pivoted["Rater"] == "R1")]
    assert len(p1_rows) == 2
    assert set(p1_rows["Criterion"]) == {"C1", "C2"}
    # The scores match the original wide cells.
    assert int(p1_rows[p1_rows["Criterion"] == "C1"]["Score"].iloc[0]) == 3
    assert int(p1_rows[p1_rows["Criterion"] == "C2"]["Score"].iloc[0]) == 4


def test_pivot_drops_blank_and_non_numeric_score_cells():
    """Blank or non-numeric cells in the wide table must be dropped
    from the long output (treated as missing, not zero)."""
    df = pd.DataFrame({
        "Person": ["P1", "P2"],
        "Rater": ["R1", "R2"],
        "C1": [3, ""],
        "C2": ["x", 4],
        "C3": [2, 3],
    })
    pivoted = app.apply_wide_to_long_pivot(
        df, id_cols=["Person", "Rater"],
        score_cols=["C1", "C2", "C3"],
        new_facet_name="Criterion", score_col_name="Score",
    )
    # 6 wide cells minus 2 blanks/non-numeric = 4 long rows.
    assert len(pivoted) == 4
    # No NaN in the Score column.
    assert pivoted["Score"].notna().all()


# -----------------------------------------------------------------------------
# Refusal on bad input
# -----------------------------------------------------------------------------


def test_pivot_refuses_overlapping_id_and_score_columns():
    df = pd.DataFrame({"Person": ["P1"], "C1": [3]})
    with pytest.raises(ValueError, match="overlap"):
        app.apply_wide_to_long_pivot(
            df, id_cols=["Person", "C1"], score_cols=["C1"],
        )


def test_pivot_refuses_empty_score_cols():
    df = pd.DataFrame({"Person": ["P1"], "C1": [3]})
    with pytest.raises(ValueError, match="score_cols is empty"):
        app.apply_wide_to_long_pivot(
            df, id_cols=["Person"], score_cols=[],
        )


def test_pivot_refuses_missing_columns():
    df = pd.DataFrame({"Person": ["P1"], "C1": [3]})
    with pytest.raises(ValueError, match="not in input"):
        app.apply_wide_to_long_pivot(
            df, id_cols=["Person"], score_cols=["NonexistentCol"],
        )


def test_pivot_refuses_name_collision_with_id_column():
    """The new facet column and the score column name cannot collide
    with an existing id column name."""
    df = pd.DataFrame({"Criterion": ["P1"], "C1": [3]})
    with pytest.raises(ValueError, match="collides with an id column"):
        app.apply_wide_to_long_pivot(
            df, id_cols=["Criterion"], score_cols=["C1"],
            new_facet_name="Criterion",
        )


# -----------------------------------------------------------------------------
# End-to-end: pivot then fit
# -----------------------------------------------------------------------------


def test_pivoted_long_format_fits_via_mfrm_estimate():
    """A wide-format Excel-style fixture pivoted via the helper must
    parse and fit successfully through the standard pipeline."""
    rng = np.random.default_rng(20260612)
    n_p, n_r, n_c = 20, 2, 4
    rows = []
    for i in range(n_p):
        theta = rng.normal(0, 1)
        for j in range(n_r):
            row = {"Person": f"P{i+1:02d}", "Rater": f"R{j+1}"}
            for k in range(n_c):
                eta = theta - (j - 0.5) * 0.3 - (k - 1.5) * 0.3
                p1 = 1.0 / (1.0 + np.exp(-eta))
                p2 = 1.0 / (1.0 + np.exp(-(eta - 0.3)))
                row[f"C{k+1}"] = int(rng.uniform() < p1) + int(rng.uniform() < p2)
            rows.append(row)
    wide_df = pd.DataFrame(rows)

    # Confirm the heuristic flags it as wide.
    assert app.detect_wide_format_columns(wide_df)["looks_wide"]

    long_df = app.apply_wide_to_long_pivot(
        wide_df, id_cols=["Person", "Rater"],
        score_cols=["C1", "C2", "C3", "C4"],
        new_facet_name="Criterion", score_col_name="Score",
    )
    assert len(long_df) == n_p * n_r * n_c

    res = app.mfrm_estimate(
        data=long_df, person_col="Person",
        facet_cols=["Rater", "Criterion"], score_col="Score",
        rating_min=0, rating_max=2, model="RSM", method="JMLE",
        maxit=15, reltol=1e-3,
    )
    assert bool(res["summary"].iloc[0]["Converged"])
    # All expected (Rater, Criterion) levels present.
    assert int(res["facets"]["others"][res["facets"]["others"]["Facet"] == "Rater"].shape[0]) == n_r
    assert int(res["facets"]["others"][res["facets"]["others"]["Facet"] == "Criterion"].shape[0]) == n_c
