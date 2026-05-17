"""Tests for the APA 7 table preset helpers and the monochrome palette.

Two presets ship side-by-side with the existing Plotly / publication
infrastructure:

* ``apa_table_to_markdown`` and ``apa_table_to_html`` re-emit any
  DataFrame as a manuscript-ready APA 7 table with the table number,
  caption (italic), header / body grid, and a ``Note.`` block. The
  output is fully self-contained: the Markdown output is GitHub-
  flavoured pipe Markdown that pastes into Quarto / RMarkdown source,
  and the HTML output uses inline CSS so it can be embedded in any
  manuscript HTML template.

* ``get_monochrome_palette`` returns a grayscale-only colour palette
  for print-friendly APA figure output; ``get_monochrome_dash_patterns``
  pairs Plotly dash patterns so distinct lines under the monochrome
  preset can be told apart by both shade and line style.

The math contract pinned here covers:

* Empty input returns the empty string for both formatters.
* The Markdown output is a valid GitHub pipe table (header row,
  separator, one body row per DataFrame row).
* Float values respect the ``digits`` argument (default 3).
* ``Note.`` and the caption appear only when supplied; the table number
  appears only when ``table_number`` is passed.
* HTML cells are HTML-escaped so ``<script>`` payloads cannot break
  out of the cell.
* Boolean values render as ``"Yes" / "No"`` and ``NaN`` renders as
  the empty string.
* The monochrome palette is a list of 8 grayscale hex strings; the
  ``n`` argument cycles / truncates correctly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


# -----------------------------------------------------------------------------
# Markdown output
# -----------------------------------------------------------------------------


def test_apa_markdown_returns_empty_string_on_empty_dataframe():
    assert app.apa_table_to_markdown(pd.DataFrame()) == ""
    assert app.apa_table_to_markdown(None) == ""


def test_apa_markdown_contains_pipe_table_structure():
    df = pd.DataFrame({"A": [1.0, 2.0], "B": [3.5, 4.5]})
    out = app.apa_table_to_markdown(df)
    # Header row, separator row, two body rows.
    lines = [line for line in out.split("\n") if line.strip()]
    pipe_lines = [line for line in lines if "|" in line]
    assert len(pipe_lines) == 4
    # Separator row uses '---'.
    assert any(set(line.replace("|", "").strip().split()) == {"---"} for line in pipe_lines)


def test_apa_markdown_includes_caption_and_note_when_supplied():
    df = pd.DataFrame({"A": [1.0]})
    out = app.apa_table_to_markdown(
        df, caption="G coefficients", note="G = relative.", table_number=2
    )
    assert "**Table 2**" in out
    assert "*G coefficients*" in out
    assert "*Note.* G = relative." in out


def test_apa_markdown_omits_caption_when_not_supplied():
    df = pd.DataFrame({"A": [1.0]})
    out = app.apa_table_to_markdown(df)
    assert "Table" not in out
    assert "*Note.*" not in out


def test_apa_markdown_respects_digits_argument():
    df = pd.DataFrame({"A": [1.234567]})
    out_default = app.apa_table_to_markdown(df)
    out_two = app.apa_table_to_markdown(df, digits=2)
    out_five = app.apa_table_to_markdown(df, digits=5)
    assert "1.235" in out_default
    assert "1.23" in out_two and "1.235" not in out_two
    assert "1.23457" in out_five


def test_apa_markdown_handles_nan_and_booleans():
    df = pd.DataFrame({
        "X": [1.0, float("nan"), 3.0],
        "Flag": [True, False, True],
    })
    out = app.apa_table_to_markdown(df)
    # NaN renders as the empty string.
    lines = [line for line in out.split("\n") if "|" in line]
    # Second body row: "1 |  |" style (empty between separators).
    nan_row = lines[3]  # header, sep, P1, P2(nan)
    assert "|  |" in nan_row.replace("| ", "|").replace(" |", "|") or "  " in nan_row
    # Booleans render as Yes/No.
    assert "Yes" in out
    assert "No" in out


# -----------------------------------------------------------------------------
# HTML output
# -----------------------------------------------------------------------------


def test_apa_html_returns_empty_string_on_empty_dataframe():
    assert app.apa_table_to_html(pd.DataFrame()) == ""
    assert app.apa_table_to_html(None) == ""


def test_apa_html_contains_table_thead_tbody_and_borders():
    df = pd.DataFrame({"A": [1.0]})
    out = app.apa_table_to_html(df)
    assert "<table" in out
    assert "<thead>" in out and "<tbody>" in out
    # APA 7 horizontal rules: top border on table, bottom border on header,
    # bottom border on table.
    assert "border-top: 2px solid" in out
    assert "border-bottom: 2px solid" in out
    assert "border-bottom: 1px solid" in out  # header row


def test_apa_html_includes_table_number_caption_and_note():
    df = pd.DataFrame({"A": [1.0]})
    out = app.apa_table_to_html(
        df, caption="Caption text", note="Note text.", table_number=3
    )
    assert "Table 3" in out
    assert "Caption text" in out
    assert "Note text." in out


def test_apa_html_escapes_html_special_characters_in_cells():
    """The cell content must be HTML-escaped so payloads like
    ``<script>`` cannot break out of the cell into a script tag."""
    df = pd.DataFrame({"A": ["<script>alert(1)</script>"]})
    out = app.apa_table_to_html(df)
    assert "<script>" not in out
    assert "&lt;script&gt;" in out


def test_apa_html_escapes_caption_and_note_special_characters():
    df = pd.DataFrame({"A": [1.0]})
    out = app.apa_table_to_html(
        df, caption="A & B <evil>", note="Note <em>tag</em> & more",
    )
    assert "<evil>" not in out
    assert "&amp;" in out


# -----------------------------------------------------------------------------
# Monochrome palette
# -----------------------------------------------------------------------------


def test_monochrome_palette_returns_eight_grayscale_hex_strings():
    palette = app.get_monochrome_palette()
    assert len(palette) == 8
    for c in palette:
        assert isinstance(c, str) and c.startswith("#") and len(c) == 7
        # Hex digits only.
        assert all(ch in "0123456789abcdefABCDEF" for ch in c[1:])


def test_monochrome_palette_truncates_when_n_below_base_length():
    palette = app.get_monochrome_palette(3)
    assert len(palette) == 3
    assert palette == list(app.get_monochrome_palette())[:3]


def test_monochrome_palette_cycles_when_n_above_base_length():
    palette = app.get_monochrome_palette(12)
    assert len(palette) == 12
    base = app.get_monochrome_palette()
    assert palette[: len(base)] == base
    assert palette[len(base):] == base[: 12 - len(base)]


def test_monochrome_palette_returns_empty_for_n_zero_or_negative():
    assert app.get_monochrome_palette(0) == []
    assert app.get_monochrome_palette(-5) == []


def test_monochrome_dash_patterns_are_valid_plotly_dashes():
    valid_dashes = {"solid", "dash", "dot", "dashdot", "longdash", "longdashdot"}
    for dash in app.get_monochrome_dash_patterns():
        assert dash in valid_dashes


def test_monochrome_dash_patterns_cycle_under_n_argument():
    out = app.get_monochrome_dash_patterns(10)
    assert len(out) == 10
    base = app.get_monochrome_dash_patterns()
    assert out[: len(base)] == base


# -----------------------------------------------------------------------------
# Round-trip: a representative report table renders cleanly
# -----------------------------------------------------------------------------


def test_apa_markdown_and_html_round_trip_a_realistic_table():
    """A realistic g/phi-style report table must produce clean output
    in both formats."""
    df = pd.DataFrame({
        "Facet": ["Person", "Rater", "Criterion"],
        "Variance": [0.171, 0.013, 0.052],
        "Proportion": [0.276, 0.021, 0.084],
        "Flag": [True, False, False],
    })
    md = app.apa_table_to_markdown(
        df, caption="G-study variance components", note="See Brennan (2001).",
        digits=3, table_number=1,
    )
    assert "**Table 1**" in md
    assert "0.171" in md and "0.276" in md
    html = app.apa_table_to_html(
        df, caption="G-study variance components", note="See Brennan (2001).",
        digits=3, table_number=1,
    )
    assert "Table 1" in html
    assert "0.171" in html and "Yes" in html and "No" in html
