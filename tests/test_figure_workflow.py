"""Figure navigation and actual saved payloads must describe the same result."""
import json
import re
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_publication_figures import synthetic_result, synthetic_diagnostics


def _view(result, diagnostics, lang):
    from unittest.mock import patch
    import streamlit as st
    import streamlit_app as app
    st.session_state["lang"] = lang
    st.session_state["downloads"] = {}
    original = st.download_button

    def capture(label, data, *args, **kwargs):
        st.session_state["downloads"][kwargs.get("key")] = data
        return original(label, data, *args, **kwargs)

    with patch.object(app, "mfrm_estimate", side_effect=AssertionError("Figure viewing must not fit")), \
         patch.object(app, "mfrm_diagnostics", side_effect=AssertionError("Figure viewing must not recompute diagnostics")), \
         patch.object(st, "download_button", capture):
        app._render_guided_figures_section(result, diagnostics,
            result_render_plots=False, result_generate_figures=False,
            action_plan=diagnostics.get("review_notes"))


def _saved_plot(at, index=0):
    data = at.session_state["downloads"]
    html = [value for key, value in data.items() if key.startswith("figure_download_")][index].decode()
    start = re.search(r'Plotly.newPlot\(\s*"mfrm-figure",\s*', html).end()
    traces, end = json.JSONDecoder().raw_decode(html[start:])
    layout, _ = json.JSONDecoder().raw_decode(html[start+end:].lstrip().lstrip(",").lstrip())
    context = json.loads(re.search(r'<script id="mfrm-figure-context" type="application/json">(.*?)</script>', html, re.S).group(1))
    rendered = json.loads(at.get("plotly_chart")[index].proto.spec)
    assert traces == rendered["data"]
    assert layout == rendered["layout"]
    assert context == layout["meta"]["mfrm"]
    assert not re.search(r'<script[^>]+src=', html)
    return context, layout, traces, html


@pytest.mark.parametrize("lang", ["ja", "en"])
@pytest.mark.parametrize("held", [False, True])
def test_core_maps_save_current_controls_without_fit(synthetic_result, synthetic_diagnostics, lang, held):
    result, diagnostics = synthetic_result, synthetic_diagnostics
    result["opt"] = SimpleNamespace(success=True)
    if held:
        result["config"]["method"] = "MML"
    diagnostics["measures"]["InfitZSTD"] = np.linspace(-3., 4., len(diagnostics["measures"]))
    original = diagnostics["measures"].copy(deep=True)
    analysis_id = app.build_result_analysis_identity(result).analysis_id
    at = AppTest.from_function(_view, args=(result, diagnostics, lang)).run()
    assert not at.exception
    assert at.selectbox(key="guided_figures_panel").value == "wright_map"
    assert not any(w.key == "wright_map_panel" for w in at.selectbox)
    context, *_ = _saved_plot(at)
    assert context["analysis_id"] == analysis_id
    assert bool(context["calibration_hold_reason"]) == held
    assert bool(at.warning) == held

    at.selectbox(key="guided_figures_panel").select("pathway_map").run()
    assert not at.exception
    context, layout, traces, _ = _saved_plot(at)
    assert len(layout["shapes"]) == 1 and layout["shapes"][0]["x0"] == 1
    at.selectbox(key="pathway_fit_statistic").select("Approx. Infit ZSTD").run()
    at.selectbox(key="pathway_display_facet").select("Rater").run()
    assert not at.exception
    context, layout, traces, _ = _saved_plot(at)
    assert sorted(shape["x0"] for shape in layout["shapes"]) == [-2, 0, 2]
    assert [trace["name"] for trace in traces] == ["Rater"]
    assert context["analysis_id"] == analysis_id
    assert "Müller" in context["reading_guide"][1]

    at.selectbox(key="guided_figures_panel").select("category_probability").run()
    assert not at.exception
    context, *_ = _saved_plot(at)
    assert context["analysis_id"] == analysis_id
    assert len(at.get("plotly_chart")) == 2
    _saved_plot(at, index=1)
    pd.testing.assert_frame_equal(diagnostics["measures"], original, check_exact=True)


def test_all_existing_figures_remain_reachable(synthetic_result, synthetic_diagnostics):
    synthetic_result["opt"] = SimpleNamespace(success=True)
    synthetic_diagnostics["obs"] = pd.DataFrame({
        "StdResidual": np.linspace(-2, 2, 20), "Expected": np.linspace(.1, 3.9, 20),
        "Observed": np.tile([0, 1, 2, 3, 4], 4),
    })
    at = AppTest.from_function(_view, args=(synthetic_result, synthetic_diagnostics, "ja")).run()
    for panel in app.GUIDED_FIGURE_PANEL_IDS:
        at.selectbox(key="guided_figures_panel").select(panel).run()
        assert not at.exception, panel
        assert at.get("plotly_chart"), panel
        _saved_plot(at)


def test_png_receives_the_exact_annotated_plot(synthetic_result, synthetic_diagnostics, monkeypatch):
    captured = []
    def png(fig, **kwargs):
        assert kwargs == {"publication": False}
        captured.append(json.loads(fig.to_json()))
        return b"\x89PNG\r\n\x1a\nserializer-test-double"
    monkeypatch.setattr(app, "fig_to_png_bytes", png)
    at = AppTest.from_function(_view, args=(synthetic_result, synthetic_diagnostics, "ja")).run()
    at.radio(key="figure_format_wright_map_0").set_value("PNG").run()
    assert not at.exception
    assert captured == [json.loads(at.get("plotly_chart")[0].proto.spec)]
    assert captured[0]["layout"]["meta"]["mfrm"]["calibration_hold_reason"]
    monkeypatch.setattr(app, "fig_to_png_bytes", lambda *a, **k: None)
    at.run()
    assert not at.exception and at.info
    assert not any(key.startswith("figure_download_") for key in at.session_state["downloads"])
    at.radio(key="figure_format_wright_map_0").set_value("HTML").run()
    _saved_plot(at)


def test_diagnostic_review_notes_survive_html_save(synthetic_result, synthetic_diagnostics):
    plan = pd.DataFrame([{"Priority": "P1", "SourceOrder": 1, "Area": "Assumptions",
        "Check": "Dimensionality", "Status": "Caution", "WhatItSays": "Review residual stability.",
        "NextAction": "Read the residual PCA.", "DetailLocation": "Dimensionality"}])
    synthetic_diagnostics["review_notes"] = plan
    at = AppTest.from_function(_view, args=(synthetic_result, synthetic_diagnostics, "ja")).run()
    assert not at.exception
    context, layout, _, html = _saved_plot(at)
    assert context["review_notes"] == plan.to_dict(orient="records")
    assert "Review residual stability." in html and "<table" in html
    assert "確認事項はHTMLに記録" in layout["annotations"][-1]["text"]


def test_nearly_identical_measures_are_not_magnified(synthetic_result, synthetic_diagnostics):
    synthetic_diagnostics["measures"]["Estimate"] = np.linspace(-1e-11, 1e-11, 15)
    at = AppTest.from_function(_view, args=(synthetic_result, synthetic_diagnostics, "ja")).run()
    at.selectbox(key="guided_figures_panel").select("pathway_map").run()
    assert not at.exception
    _, layout, _, _ = _saved_plot(at)
    assert layout["yaxis"]["range"] == [-.5, .5]


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_reference_band_counts_do_not_assign_overall_fit_verdict(lang):
    def view(lang, count):
        import streamlit as st
        import pandas as pd
        import streamlit_app as app
        st.session_state["lang"] = lang
        app._fit_summary_callout(pd.DataFrame({"Infit": [1.51] * count + [1.] * (100-count)}))
    for count in (0, 10, 25, 26, 100):
        at = AppTest.from_function(view, args=(lang, count)).run()
        assert not at.exception
        assert not at.error and not at.warning and not at.success
        assert len(at.caption) == 1 and str(count) in at.caption[0].value
