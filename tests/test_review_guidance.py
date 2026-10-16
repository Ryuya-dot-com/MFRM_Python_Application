"""Regression checks for the Japanese first-read path and sample guidance."""
from pathlib import Path

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


def test_japanese_sample_run_has_translated_first_read_and_current_testlet_scope():
    at = AppTest.from_file(Path(app.__file__)).run(timeout=60)
    at.radio(key="lang").set_value("ja").run(timeout=60)
    at.button(key="onboarding_skip_guide").click().run(timeout=60)
    at.button(key="facets_mode_run_primary").click().run(timeout=120)
    assert not at.exception
    text = "\n".join(str(x.value) for kind in ("info", "warning", "error", "markdown") for x in getattr(at, kind))
    assert "残差の第1主成分の固有値" in text
    for stale in ("Current run", "First residual PCA eigenvalue", "Stability audit:", "estimand契約", "G3適格性"):
        assert stale not in text
    details = next(e for e in at.expander if e.label == "再現のための実行情報")
    assert not details.proto.expanded
    at.selectbox(key="data_source_scenario").set_value("reading_testlet_binary").run(timeout=60)
    assert not at.exception
    descriptions = "\n".join(str(x.value) for x in at.markdown)
    assert "項目群のランダム効果の推定やStanコード生成の操作画面はありません" in descriptions
    assert "TESTLET_RI" not in descriptions


def test_display_localization_keeps_numerical_statuses_and_values(monkeypatch):
    from types import SimpleNamespace
    result = {"opt": SimpleNamespace(success=True), "config": {"method": "JMLE"}}
    diagnostics = {"pca_enabled": True, "pca": {"eigenvalues": [3.74],
        "stability_table": pd.DataFrame([{"PCAStabilityStatus": "Review",
            "SensitivityFlagSummary": "bootstrap loading correlation 0.06"}])}}
    monkeypatch.setitem(app.st.session_state, "lang", "ja")
    raw = app.build_first_read_guide_rows(result, diagnostics)
    localized = app.build_first_read_guide_rows(result, diagnostics, localize=True)
    assert [(r["Check"], r["Status"]) for r in raw] == [(r["Check"], r["Status"]) for r in localized]
    row = next(r for r in localized if r["Check"] == "4. Dimensionality")
    assert "3.74" in row["What it means"] and "要確認" in row["What it means"]
    assert "First residual PCA" not in row["What it means"]


@pytest.mark.parametrize("lang", ["en", "ja"])
def test_real_facet_boundary_hold_reaches_report_screen(lang):
    def view(lang):
        import pandas as pd
        import streamlit as st
        import streamlit_app as app
        from tests.test_jmle_facet_boundary import separated_data

        st.session_state["lang"] = lang
        if "fit" not in st.session_state:
            st.session_state["fit"] = app.mfrm_estimate(
                separated_data(), person_col="Person", facet_cols=["Rater", "Task"],
                score_col="Score", rating_min=0, rating_max=1, keep_original=True,
            )
            st.session_state["diagnostics"] = app.mfrm_diagnostics(st.session_state["fit"], compute_pca=False)
        st.session_state["expected_hold"] = app.t("guided.export_state_hold")
        st.session_state["draft_ready"] = app.t("guided.export_state_draft")
        app._render_guided_report_export_section(
            st.session_state["fit"], st.session_state["diagnostics"], {},
            pd.DataFrame(), pd.DataFrame(), None, {}, generate_figures=False,
        )

    at = AppTest.from_function(view, args=(lang,)).run(timeout=60)
    assert not at.exception
    assert at.session_state["expected_hold"] in [warning.value for warning in at.warning]
    assert at.session_state["draft_ready"] not in [caption.value for caption in at.caption]
    at.segmented_control(key="guided_export_task").set_value("review").run(timeout=60)
    assert not at.exception
    assert at.session_state["expected_hold"] in [warning.value for warning in at.warning]
