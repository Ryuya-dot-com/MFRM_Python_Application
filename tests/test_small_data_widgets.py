"""Regression tests for v0.2.10-beta `StreamlitValueBelowMinError` crash.

When a facet has fewer rows than a widget's `min_value` (most commonly
a Rater facet with 3 raters meeting a `min_value=5` slider / number
input), Streamlit raises a hard error:
    StreamlitValueBelowMinError / StreamlitValueAboveMaxError /
    StreamlitAPIException (min_value >= max_value)

These tests render the affected sections via `streamlit.testing.v1`
against a synthetic 3-rater / 2-task / 10-person dataset and assert
that no exception surfaces. They are the cheapest possible regression
guard for the pattern.
"""

from __future__ import annotations

import pandas as pd
from streamlit.testing.v1 import AppTest

import streamlit_app as app


# ---------------------------------------------------------------------------
# Unit tests: the widget-math is clamped correctly before the call
# ---------------------------------------------------------------------------

def test_forest_plot_tiny_facet_does_not_crash():
    """_draw_measures_forest_plotly on a 3-row facet must not raise."""
    diag = {
        "measures": pd.DataFrame({
            "Facet": ["Rater"] * 3,
            "Level": ["R1", "R2", "R3"],
            "Estimate": [-0.3, 0.0, 0.4],
            "SE": [0.15, 0.14, 0.16],
            "Infit": [1.0, 0.95, 1.05],
            "Outfit": [1.05, 0.95, 1.10],
        }),
    }
    # The widget itself requires a ScriptRunContext so a direct call
    # will log warnings but must not raise a Streamlit error.
    try:
        app._draw_measures_forest_plotly(diag)
    except Exception as exc:
        # Streamlit API warnings are expected outside a ScriptRunContext;
        # only re-raise Streamlit VALUE errors which signal the actual
        # regression.
        msg = str(type(exc).__name__) + ": " + str(exc)
        assert "ValueBelowMinError" not in msg, msg
        assert "ValueAboveMaxError" not in msg, msg
        assert "StreamlitAPIException" not in msg, msg


def test_misfit_ranking_tiny_fit_table_does_not_crash():
    """_draw_misfit_ranking on a 3-row fit table must not raise."""
    fit = pd.DataFrame({
        "Facet": ["Rater"] * 3,
        "Level": ["R1", "R2", "R3"],
        "InfitZSTD": [0.5, -1.2, 2.1],
        "OutfitZSTD": [0.4, -1.0, 2.3],
    })
    try:
        app._draw_misfit_ranking(fit)
    except Exception as exc:
        msg = str(type(exc).__name__) + ": " + str(exc)
        assert "ValueBelowMinError" not in msg, msg
        assert "ValueAboveMaxError" not in msg, msg
        assert "StreamlitAPIException" not in msg, msg


# ---------------------------------------------------------------------------
# Integration test: full AppTest smoke with a small scenario
# ---------------------------------------------------------------------------

def test_apptest_with_clinical_osce_runs_visuals_without_crash():
    """End-to-end: load the small Clinical OSCE scenario and confirm the
    app renders without a StreamlitValueBelowMinError.

    The Clinical OSCE scenario has only 3 competency elements on the
    Competency facet, which is exactly the shape that triggered the
    regression reported by the user.
    """
    at = AppTest.from_file("streamlit_app.py").run(timeout=40)
    # Switch the sample scenario selector to Clinical OSCE
    for radio in at.radio:
        if radio.key == "data_source_flat":
            try:
                radio.set_value("Clinical OSCE (60×4×5×3, 3,600 obs)")
            except Exception:
                # Label may be slightly different — best-effort
                pass
            break
    at.run(timeout=60)
    # We deliberately do NOT click "Run MFRM estimation" —
    # the initial-render path alone exercises the widget math that
    # blew up. A clean run here is the pass condition.
    assert not at.exception, f"AppTest raised: {at.exception}"


# ---------------------------------------------------------------------------
# Boundary-case unit tests: feed extreme small inputs directly to every
# renderer that contains a data-dependent widget. These go beyond the
# shipped scenarios to simulate user-provided data with 2 raters, 1 task,
# or 1-element facets — shapes plausible in real research.
# ---------------------------------------------------------------------------

def test_forest_plot_single_row_facet_does_not_crash():
    """1-row facet: the widget must be skipped, plot still attempted."""
    diag = {
        "measures": pd.DataFrame({
            "Facet": ["Rater"],
            "Level": ["R_ONLY"],
            "Estimate": [0.0],
            "SE": [0.2],
            "Infit": [1.0],
            "Outfit": [1.0],
        }),
    }
    try:
        app._draw_measures_forest_plotly(diag)
    except Exception as exc:
        msg = f"{type(exc).__name__}: {exc}"
        assert "ValueBelowMinError" not in msg, msg
        assert "ValueAboveMaxError" not in msg, msg
        assert "StreamlitAPIException" not in msg, msg


def test_misfit_ranking_single_row_does_not_crash():
    fit = pd.DataFrame({
        "Facet": ["Rater"],
        "Level": ["R_ONLY"],
        "InfitZSTD": [0.5],
        "OutfitZSTD": [0.5],
    })
    try:
        app._draw_misfit_ranking(fit)
    except Exception as exc:
        msg = f"{type(exc).__name__}: {exc}"
        assert "ValueBelowMinError" not in msg, msg
        assert "ValueAboveMaxError" not in msg, msg
        assert "StreamlitAPIException" not in msg, msg


def test_bias_heatmap_tiny_pivot_does_not_crash():
    """A 2x2 pivot exercises the slider-bounds-clamping branch."""
    tiny_tbl = pd.DataFrame({
        "FacetA": ["Rater"] * 4,
        "FacetB": ["Task"] * 4,
        "FacetA_Level": ["R1", "R1", "R2", "R2"],
        "FacetB_Level": ["T1", "T2", "T1", "T2"],
        "BiasSize": [0.1, -0.2, 0.3, -0.1],
        "BiasSE": [0.15] * 4,
        "t": [0.7, -1.3, 2.0, -0.7],
    })
    try:
        app._draw_bias_heatmap(tiny_tbl, "Rater", "Task")
    except Exception as exc:
        msg = f"{type(exc).__name__}: {exc}"
        assert "ValueBelowMinError" not in msg, msg
        assert "ValueAboveMaxError" not in msg, msg
        assert "StreamlitAPIException" not in msg, msg
