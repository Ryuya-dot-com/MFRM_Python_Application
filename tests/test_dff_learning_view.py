"""DFF reading views preserve fitted values, missingness and inference holds."""
from io import BytesIO
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_analysis_journeys import _run
from tests.test_help_contextual_app_integration import _block_analysis_calls, _locale_value


@pytest.fixture(scope="module", params=["JMLE", "MML"])
def fitted(request):
    method = request.param
    result = app.mfrm_estimate(
        app._make_self_test_gradient_data(), person_col="Person", score_col="Score",
        facet_cols=["Rater", "Task"], method=method, model="RSM",
        noncenter_facet="Person" if method == "JMLE" else "Task",
        mml_engine="direct", quad_points=15, maxit=200, reltol=1e-6,
    )
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False, compute_marginal=False)
    bias = app.estimate_bias_interaction(result, diagnostics, "Rater", "Task")
    assert bool(result["summary"].iloc[0].Converged)
    return result, diagnostics, bias


@pytest.mark.parametrize("language", ["ja", "en"])
def test_first_view_puts_plot_and_five_columns_before_collapsed_details(fitted, language, monkeypatch):
    result, diagnostics, bias = fitted
    original_table = bias["table"].copy(deep=True)
    calls = _block_analysis_calls(monkeypatch)
    figures, downloads = [], {}
    original_plot, original_download = app.st.plotly_chart, app.st.download_button

    def capture_plot(fig, **kwargs):
        figures.append(fig)
        return original_plot(fig, **kwargs)

    def capture_download(label, data, *args, **kwargs):
        downloads[kwargs.get("file_name")] = data
        return original_download(label, data, *args, **kwargs)

    monkeypatch.setattr(app.st, "plotly_chart", capture_plot)
    monkeypatch.setattr(app.st, "download_button", capture_download)
    at = AppTest.from_string(
        'import streamlit as st\nimport streamlit_app as app\n'
        'app.show_bias_section(st.session_state["bias"], {}, '
        'result=st.session_state["result"], diagnostics=st.session_state["diagnostics"])',
        default_timeout=120,
    )
    for key, value in dict(lang=language, result=result, diagnostics=diagnostics, bias=bias).items():
        at.session_state[key] = value
    _run(at)
    children = list(at.main.children.values())
    if len(children) == 1 and children[0].type == "flex_container" and hasattr(children[0], "children"):
        children = list(children[0].children.values())
    details = next(e for e in at.expander if e.label == _locale_value(language, "bias_interaction.technical_details_label"))
    assert not details.proto.expanded
    assert len([e for e in children if e.type == "dataframe"]) == 1
    assert next(i for i, e in enumerate(children) if e.type == "plotly_chart") < children.index(details)
    expected = original_table[["FacetA_Level", "FacetB_Level", "PositiveWeightObsN", "Obs-Exp Average", "Bias Size"]].rename(columns={
        "FacetA_Level": "Rater", "FacetB_Level": "Task",
        "PositiveWeightObsN": _locale_value(language, "bias_interaction.cell_count_label"),
        "Obs-Exp Average": _locale_value(language, "bias_interaction.observed_difference_label"),
        "Bias Size": _locale_value(language, "bias_interaction.heatmap_metric_bias_size_display"),
    })
    pd.testing.assert_frame_equal(at.dataframe[0].value, expected)
    assert _locale_value(language, "bias_interaction.interpretation_hold") in [e.value for e in at.info]
    assert _locale_value(language, "bias_interaction.reading_guide") in [e.value for e in at.markdown]
    heatmap = next(fig for fig in figures if fig.data[0].type == "heatmap")
    assert _locale_value(language, "bias_interaction.heatmap_metric_bias_size_display") in heatmap.layout.title.text
    assert heatmap.data[0].text is None
    assert _locale_value(language, "bias_interaction.cell_count_label") in str(heatmap.data[0].customdata)

    figures.clear()
    at.radio(key="bias_heatmap_metric").set_value("Standardized Bias")
    _run(at)
    heatmap = next(fig for fig in figures if fig.data[0].type == "heatmap")
    expected_z = original_table.pivot_table(index="FacetA_Level", columns="FacetB_Level", values="Standardized Bias", dropna=False)
    np.testing.assert_allclose(heatmap.data[0].z, expected_z, equal_nan=True)
    pd.testing.assert_frame_equal(at.dataframe[0].value, expected)
    pd.testing.assert_frame_equal(bias["table"], original_table, check_exact=True)
    assert downloads["mfrm_bias_table.csv"] == app.to_csv_bytes(original_table)
    screening = pd.read_csv(BytesIO(downloads["mfrm_dff_bias_screening.csv"]))
    assert screening.ClaimStatus.eq("Do not claim").all()
    assert calls == {"estimate": 0, "refit": 0}


@pytest.mark.parametrize("language", ["ja", "en"])
def test_skipped_bias_calculation_explains_how_to_run_it(language):
    at = AppTest.from_string("import streamlit_app as app\napp.show_bias_section(None, {})")
    at.session_state["lang"] = language
    _run(at)
    assert _locale_value(language, "bias_interaction.no_results_info") in [e.value for e in at.info]
    assert not at.dataframe


@pytest.mark.parametrize("count_source", ["PositiveWeightObsN", "ObsN", "Observd Count"])
def test_reading_table_and_hover_distinguish_counts_from_weight_sum(fitted, monkeypatch, count_source):
    """Saved results with different count metadata must never rename weights as N."""
    result, diagnostics, original_bias = fitted
    bias = deepcopy(original_bias)
    table = bias["table"]
    table["PositiveWeightObsN"], table["ObsN"], table["Observd Count"] = 2, 3, 999.5
    if count_source != "PositiveWeightObsN":
        table.drop(columns="PositiveWeightObsN", inplace=True)
    if count_source == "Observd Count":
        table.drop(columns="ObsN", inplace=True)
    figures = []
    original_plot = app.st.plotly_chart
    def capture_plot(fig, **kwargs):
        figures.append(fig)
        return original_plot(fig, **kwargs)
    monkeypatch.setattr(app.st, "plotly_chart", capture_plot)
    calls = _block_analysis_calls(monkeypatch)
    at = AppTest.from_string(
        'import streamlit as st\nimport streamlit_app as app\n'
        'app.show_bias_section(st.session_state["bias"], {}, '
        'result=st.session_state["result"], diagnostics=st.session_state["diagnostics"])',
        default_timeout=120,
    )
    for key, value in dict(lang="ja", result=result, diagnostics=diagnostics, bias=bias).items():
        at.session_state[key] = value
    _run(at)
    key = "weight_sum_label" if count_source == "Observd Count" else "cell_count_label"
    label = _locale_value("ja", f"bias_interaction.{key}")
    expected = {"PositiveWeightObsN": 2, "ObsN": 3, "Observd Count": 999.5}[count_source]
    assert at.dataframe[0].value[label].eq(expected).all()
    hover = str(next(fig for fig in figures if fig.data[0].type == "heatmap").data[0].customdata)
    displayed = f"{expected:.3f}" if count_source == "Observd Count" else str(expected)
    assert f"{label}: {displayed}<br>" in hover
    if count_source != "Observd Count":
        assert "999.5" not in hover
    assert calls == {"estimate": 0, "refit": 0}
