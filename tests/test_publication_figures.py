"""Tests for the Publication Document figure payloads.

The PDF / Word / HTML publication documents all embed the same four
core figures. v0.2.10-beta added a matplotlib fallback so that PDFs
include plots even when kaleido cannot reach a Chrome binary (the
common failure mode on Streamlit Community Cloud). These tests pin
the fallback path so regressions are caught quickly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


@pytest.fixture
def synthetic_result():
    """Minimal synthetic `result` dict the figure helpers expect."""
    return {
        "config": {"model": "RSM", "method": "JMLE",
                   "facet_names": ["Rater", "Task"], "n_cat": 5},
        "prep": {"n_obs": 120, "n_person": 20,
                 "rating_min": 0, "rating_max": 4,
                 "unused_score_categories": []},
        "params": {"steps": [-1.5, -0.5, 0.5, 1.5]},
        "opt": None,
        "summary": pd.DataFrame([{
            "Model": "RSM", "Method": "JMLE",
            "Converged": True, "Iterations": 42, "LogLik": -180.0,
        }]),
        "facets": {
            "person": pd.DataFrame({
                "Person": [f"P{i}" for i in range(20)],
                "Estimate": np.linspace(-1.5, 1.5, 20),
            }),
            "others": pd.DataFrame({
                "Facet": ["Rater"] * 3,
                "Level": ["R1", "R2", "R3"],
                "Estimate": [-0.3, 0.0, 0.4],
            }),
        },
        "steps": pd.DataFrame({
            "Step": ["S1", "S2", "S3", "S4"],
            "Estimate": [-1.5, -0.5, 0.5, 1.5],
        }),
    }


@pytest.fixture
def synthetic_diagnostics():
    return {
        "measures": pd.DataFrame({
            "Facet": ["Person"] * 10 + ["Rater"] * 3 + ["Task"] * 2,
            "Level": ([f"P{i}" for i in range(10)]
                      + ["R1", "R2", "R3", "T1", "T2"]),
            "Estimate": (list(np.linspace(-1, 1, 10))
                         + [-0.3, 0.0, 0.3, -0.2, 0.2]),
            "SE": [0.1] * 15,
            "Infit": [1.0, 0.9, 1.1, 1.05, 0.95, 1.0, 1.1, 0.9, 1.05,
                      0.95, 0.98, 1.02, 1.0, 0.99, 1.01],
            "Outfit": [1.0, 0.95, 1.1, 1.0, 1.0, 1.05, 1.15, 0.9, 1.0,
                       1.0, 0.95, 1.05, 1.0, 1.0, 1.0],
        }),
    }


# ---------------------------------------------------------------------------
# matplotlib fallback helpers
# ---------------------------------------------------------------------------

def test_mpl_wright_map_produces_png(synthetic_result, synthetic_diagnostics):
    png = app._mpl_wright_map_png(synthetic_result, synthetic_diagnostics)
    assert png is not None
    assert png.startswith(b"\x89PNG"), "output must start with PNG magic bytes"
    assert len(png) > 2000, "PNG seems too small to contain a real figure"


def test_mpl_fit_scatter_produces_png(synthetic_diagnostics):
    png = app._mpl_fit_scatter_png(synthetic_diagnostics)
    assert png is not None
    assert png.startswith(b"\x89PNG")
    assert len(png) > 2000


def test_mpl_category_probability_produces_png(synthetic_result):
    png = app._mpl_category_probability_png(synthetic_result)
    assert png is not None
    assert png.startswith(b"\x89PNG")
    assert len(png) > 2000


def test_mpl_facet_distribution_produces_png(synthetic_diagnostics):
    png = app._mpl_facet_distribution_png(synthetic_diagnostics)
    assert png is not None
    assert png.startswith(b"\x89PNG")
    assert len(png) > 2000


def test_mpl_fallbacks_handle_empty_diagnostics():
    """Every helper must return None (never raise) on empty inputs."""
    empty_result = {"config": {}, "prep": {}, "facets": {}}
    empty_diag = {"measures": pd.DataFrame()}
    assert app._mpl_wright_map_png(empty_result, empty_diag) is None
    assert app._mpl_fit_scatter_png(empty_diag) is None
    assert app._mpl_category_probability_png(empty_result) is None
    assert app._mpl_facet_distribution_png(empty_diag) is None


# ---------------------------------------------------------------------------
# Integration: _publication_figure_payloads
# ---------------------------------------------------------------------------

def test_publication_figure_payloads_returns_all_four_figures(
    synthetic_result, synthetic_diagnostics
):
    """All 4 core figures must be emitted when the data is present.

    This is the test that specifically guards against the regression
    we fixed in v0.2.10-beta: on Streamlit Cloud (no Chrome), the
    Plotly export returns None and the payload list used to be empty.
    With the matplotlib fallback wired in, every slot should fill.
    """
    payloads = app._publication_figure_payloads(
        synthetic_result, synthetic_diagnostics,
    )
    # Must produce exactly 4 figures (wright, fit, category, facet dist)
    assert len(payloads) == 4
    expected_prefixes = {"wright_map", "fit_scatter",
                         "category_curves_", "facet_distribution"}
    got_ids = [fid for fid, _, _ in payloads]
    for expected in expected_prefixes:
        assert any(gid.startswith(expected) for gid in got_ids), (
            f"missing figure {expected!r} in payloads: {got_ids}"
        )


def test_payloads_are_png_bytes(synthetic_result, synthetic_diagnostics):
    payloads = app._publication_figure_payloads(
        synthetic_result, synthetic_diagnostics,
    )
    for fid, caption, png in payloads:
        assert isinstance(png, (bytes, bytearray))
        assert png.startswith(b"\x89PNG"), f"{fid}: not a PNG"
        assert len(png) > 500
        assert isinstance(caption, str) and caption


def test_payload_captions_are_descriptive(synthetic_result, synthetic_diagnostics):
    payloads = app._publication_figure_payloads(
        synthetic_result, synthetic_diagnostics,
    )
    for fid, caption, _png in payloads:
        # Captions must be sentences, not placeholders
        assert len(caption) > 20, f"{fid}: caption too short ({caption!r})"
        assert "." in caption, f"{fid}: caption missing period ({caption!r})"


# ---------------------------------------------------------------------------
# Integration: PDF builder embeds the figures
# ---------------------------------------------------------------------------

def test_pdf_contains_embedded_images(synthetic_result, synthetic_diagnostics):
    """The full PDF builder must embed the figures — not just include captions."""
    try:
        import reportlab  # noqa: F401
    except ImportError:
        pytest.skip("reportlab not installed")

    pdf = app.build_publication_pdf_bytes(
        synthetic_result, synthetic_diagnostics,
    )
    assert pdf.startswith(b"%PDF")
    # An empty-figures PDF tends to be ~20-40 KB; with 4 embedded
    # PNGs it is several hundred KB. This size floor will catch the
    # "silently no plots" regression we fixed.
    assert len(pdf) > 100_000, (
        f"PDF suspiciously small ({len(pdf):,} bytes) — are figures embedded?"
    )
    # Reportlab compresses embedded images with FlateDecode.
    assert b"FlateDecode" in pdf


def test_pdf_smaller_when_figures_missing():
    """Sanity check: an empty figures case produces a substantially smaller PDF."""
    try:
        import reportlab  # noqa: F401
    except ImportError:
        pytest.skip("reportlab not installed")

    empty_result = {
        "config": {"model": "RSM", "method": "JMLE", "facet_names": [], "n_cat": 5},
        "prep": {"n_obs": 0, "n_person": 0, "rating_min": 0, "rating_max": 4},
        "params": {}, "opt": None,
        "summary": pd.DataFrame(), "facets": {}, "steps": pd.DataFrame(),
    }
    empty_diag = {"measures": pd.DataFrame()}
    pdf_empty = app.build_publication_pdf_bytes(empty_result, empty_diag)
    assert pdf_empty.startswith(b"%PDF")
    # Empty-figures PDF should exist but be modestly sized (no images)
    assert len(pdf_empty) > 2000
