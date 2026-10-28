"""Input guidance and real Web/native figure parity without extra model fits."""
from contextlib import nullcontext
from types import SimpleNamespace
import re

import numpy as np
import pandas as pd
import pytest

from desktop_app import demonstration_data
from mfrm_app import basic_analysis as basic, desktop_guidance as guide, visual_maps as maps


@pytest.fixture(scope="module")
def analysis():
    return basic.run_basic_analysis(demonstration_data(), person_col="Person", facet_cols=["Rater", "Task"],
                                    score_col="Score", rating_min=0, rating_max=3)


def test_csv_preserves_japanese_and_literal_identifiers(tmp_path):
    path = tmp_path / "評定.csv"
    path.write_text("対象者,評価者,得点\n001,NA,0\n002,NA,1\n", encoding="utf-8-sig")
    frame = guide.read_desktop_csv(path)
    assert frame["対象者"].tolist() == ["001", "002"]
    assert frame["評価者"].tolist() == ["NA", "NA"]
    assert guide.input_review(frame, dict(person_col="対象者", facet_cols=["評価者"], score_col="得点",
                                         rating_min=0, rating_max=5))["included"] == 2


@pytest.mark.parametrize("content,match", [
    ("Person,Rater,Rater\nP1,R1,2", "重複"),
    ("Person,Rater,Score\nP1,R1,2,3", "列数"),
    ("Person,Rater,Score\n", "評定データ"),
    ("Person,,Score\nP1,R1,0", "列名"),
    ("Person,Rater,Score\nP1,\"R1,0", "引用符"),
])
def test_malformed_csv_has_actionable_error(tmp_path, content, match):
    path = tmp_path / "bad.csv"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError, match=match):
        guide.read_desktop_csv(path)


def test_csv_encoding_is_explained(tmp_path):
    path = tmp_path / "shift_jis.csv"
    path.write_bytes("対象者,評価者,得点\n001,01,0".encode("cp932"))
    with pytest.raises(ValueError, match="CSV UTF-8"):
        guide.read_desktop_csv(path)


def test_input_review_matches_shared_exclusions_and_detects_duplicate_cells():
    frame = pd.DataFrame({"P": ["001", "002", "", "003", "001"],
                          "R": ["R1"] * 5, "S": ["0", "2", "1", "NA", "1"]})
    reviewed = guide.input_review(frame, dict(person_col="P", facet_cols=["R"], score_col="S", rating_min=0, rating_max=5))
    assert reviewed["included"] == 3 and reviewed["excluded"] == 2
    assert reviewed["duplicate_rows"] == 2 and len(reviewed["messages"]) == 2


@pytest.mark.parametrize("value,match", [(1.5,"小数"), (np.inf,"無限大"), (6,"指定範囲")])
def test_score_errors_precede_worker(value, match):
    frame = pd.DataFrame({"P": ["A", "B"], "R": ["R1", "R1"], "S": [0, value]})
    with pytest.raises(ValueError, match=match):
        guide.input_review(frame, dict(person_col="P", facet_cols=["R"], score_col="S", rating_min=0, rating_max=5))


def test_bundled_example_is_the_web_example():
    import streamlit_app as web
    expected = web.sample_mfrm_data_by_key("writing_essay", seed=20240101).astype(str)
    pd.testing.assert_frame_equal(guide.writing_example(), expected)


def test_wright_actual_web_and_native_coordinates_agree(analysis):
    import streamlit_app as web
    facets = analysis["result"]["facets"]
    browser = web._make_wright_map_export_figure(facets["person"], facets["others"])
    native = basic.basic_figure(analysis, "wright")
    np.testing.assert_array_equal(browser.data[0].y, facets["person"]["Estimate"])
    for artist, rows in native._mfrm_hover_points:
        trace = next(trace for trace in browser.data if trace.name == rows["Facet"].iloc[0])
        np.testing.assert_array_equal(trace.y, np.asarray(artist.get_offsets())[:, 1])
        assert list(trace.customdata[:, 0]) == rows["Level"].astype(str).tolist()


@pytest.mark.parametrize("statistic,label", [("Infit", "Infit MnSq"), ("InfitZSTD", "Approx. Infit ZSTD")])
def test_pathway_actual_web_and_native_coordinates_agree(analysis, monkeypatch, statistic, label):
    import streamlit_app as web
    captured = []
    noop = lambda *args, **kwargs: None
    ui = SimpleNamespace(session_state={"lang": "en"},columns=lambda *a: [nullcontext(), nullcontext()], expander=lambda *a, **k: nullcontext(),
                         subheader=noop, caption=noop, markdown=noop, info=noop, selectbox=lambda *a, **k: label,
                         plotly_chart=lambda fig, **k: captured.append(fig))
    monkeypatch.setattr(web, "st", ui)
    monkeypatch.setattr(web, "render_help_popover", noop)
    monkeypatch.setattr(web, "render_chart_guide", noop)
    web._draw_pathway_map_plotly(analysis["diagnostics"])
    native = basic.basic_figure(analysis, "pathway", statistic=statistic)
    for artist, rows in native._mfrm_hover_points:
        trace = next(trace for trace in captured[0].data if trace.name == rows["Facet"].iloc[0])
        np.testing.assert_array_equal(np.column_stack([trace.x, trace.y]), artist.get_offsets())
    text = str(web._CHART_GUIDE_LIBRARY["pathway_map"])
    assert "Infit MnSq (default)" in text and "1.4" not in text
    references = sorted(float(line.get_xdata()[0]) for line in native.axes[0].lines)
    assert references == ([1.] if statistic == "Infit" else [-2., 0., 2.])
    np.testing.assert_array_equal(captured[0].layout.xaxis.range, native.axes[0].get_xlim())
    assert not native.axes[0].patches  # No shaded pass region.


def test_nonfinite_coordinates_and_unconfigured_facet_are_handled():
    measures = pd.DataFrame({"Facet": ["評価者"] * 4, "Level": ["a", "b", "c", "d"],
                             "Estimate": [0, np.inf, 2, 3], "InfitZSTD": [1, 2, np.nan, 4]})
    fig = maps.pathway_figure({}, {"measures": measures}, statistic="InfitZSTD")
    assert len(fig._mfrm_hover_points) == 1
    np.testing.assert_array_equal(fig._mfrm_hover_points[0][0].get_offsets(), [[1, 0], [4, 3]])
    empty = maps.pathway_figure({}, {})
    assert any("No finite" in text.get_text() for text in empty.axes[0].texts)


def test_maps_exports_and_html_never_refit_and_svg_ids_are_unique(analysis, monkeypatch):
    def forbidden(*a, **k):
        raise AssertionError("Rendering must not fit or recompute diagnostics")
    monkeypatch.setattr(basic.engine, "mfrm_estimate", forbidden)
    monkeypatch.setattr(basic.engine, "mfrm_diagnostics", forbidden)
    for kind in maps.MAP_GUIDES:
        assert basic.figure_bytes(analysis, kind, format="png").startswith(b"\x89PNG")
    document = basic.basic_report_html(analysis)
    assert document.count("<svg") == 3
    ids = re.findall(r'\bid="([^"]+)"', document)
    assert len(ids) == len(set(ids))
    for target in re.findall(r'(?:href="#|url\(#)([^"\)]+)', document):
        assert target in ids
    assert f'href="{guide.WEB_APP_URL}"' in document
    assert "±2" in document


def test_map_hold_is_visible_and_no_uncertainty_is_implied(analysis, monkeypatch):
    monkeypatch.setattr(maps.output_qualification, "calibration_inference_hold_reason", lambda result: "held")
    for kind in ("wright", "pathway", "category"):
        fig = basic.basic_figure(analysis, kind)
        assert any("Successful estimation has not been confirmed" in text.get_text() for text in fig.texts)
        assert all(not axes.containers or all(type(c).__name__ != "ErrorbarContainer" for c in axes.containers)
                   for axes in fig.axes)


def test_report_uses_the_selected_fit_statistic_and_facet(analysis):
    default = basic.basic_report_html(analysis)
    selected = basic.basic_report_html(analysis, pathway_statistic="InfitZSTD", pathway_facet="Rater")
    assert "<!-- Infit MnSq -->" in default
    assert "<!-- Approximate Infit ZSTD -->" not in default
    assert "<!-- Approximate Infit ZSTD -->" in selected
    assert "<!-- Pathway map · Rater -->" in selected


def test_zstd_depends_on_degrees_of_freedom_at_the_same_mnsq():
    from mfrm_app.engine import zstd_from_mnsq_facets
    assert abs(zstd_from_mnsq_facets(1.2, 20)) < 2
    assert zstd_from_mnsq_facets(1.2, 500) > 2

@pytest.mark.parametrize("lang", ["ja", "en"])
@pytest.mark.parametrize("chart,kind,reference", [
    ("wright_map", "wright", "Boone_2008"),
    ("pathway_map", "pathway", "Muller_2020"),
    ("category_probability", "category", "Linacre_RatingScale_2002"),
])
def test_web_reading_guides_share_native_text_and_resolvable_sources(lang, chart, kind, reference):
    from streamlit.testing.v1 import AppTest
    import streamlit_app as web

    def view(lang, chart):
        import streamlit as st
        import streamlit_app as web
        st.session_state["lang"] = lang
        web.render_chart_guide(chart)

    at = AppTest.from_function(view, args=(lang, chart)).run()
    assert not at.exception
    guide = (maps.MAP_GUIDES if lang == "ja" else maps.MAP_GUIDES_EN)[kind]
    assert any(text.value == guide[1] for text in at.markdown)
    assert not at.expander[0].proto.expanded
    assert web._APA_REFERENCE_LIBRARY[reference] in web.collect_cited_references(guide[1])
