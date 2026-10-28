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
        st.session_state["jmle_hold"] = app.t("guided.review_hold_jmle")
        gate = app.build_publication_gate_summary(st.session_state["fit"], st.session_state["diagnostics"], {})
        st.session_state["source_findings"] = gate.loc[
            gate["GateArea"].ne("Overall manuscript gate"), "Evidence"].tolist()
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
    assert at.session_state["jmle_hold"] in [item.value for item in at.markdown]
    detail_text = '\n'.join(item.value for expander in at.expander for item in expander.markdown)
    assert all(text in detail_text for text in at.session_state["source_findings"])


def test_review_grouping_preserves_distinct_findings_statuses_and_unknown_checks():
    source = pd.DataFrame([
        ('Overall manuscript gate', 'Not ready', 'overall', 'overall action'),
        ('Reliability / separation', 'Report with caveat', 'rater separation', 'do not call it agreement'),
        ('Bias / local interaction', 'Report with caveat', 'screen not run', 'record it as unexamined'),
        ('Reliability and separation', 'Boundary', 'person uncertainty withheld', 'keep uncertainty hold'),
        ('Bias / local interaction', 'Boundary', 'inference unqualified', 'no bias or no-bias conclusion'),
        ('Convergence', 'Not ready', 'termination unknown', 'inspect termination'),
        ('Convergence and final interpretability', 'Not ready', 'identification unqualified', 'inspect design'),
        ('Future diagnostic A', 'Future status', 'new evidence A', 'new action A'),
        ('Future diagnostic B', 'Ready', 'new evidence B', 'new action B'),
    ], columns=['GateArea', 'GateStatus', 'Evidence', 'ManuscriptAction'])
    # A duplicated index must not cause either duplication or loss during grouping.
    source.index = [0] * len(source)
    before = source.copy(deep=True)
    grouped = app._group_guided_report_checks(source)
    assert [key for key, _ in grouped] == ['calibration', 'reliability', 'bias', 'other', 'other']
    restored = pd.concat([frame for _, frame in grouped])
    assert sorted(restored.itertuples(index=False, name=None)) == sorted(source.iloc[1:].itertuples(index=False, name=None))
    pd.testing.assert_frame_equal(source, before)
    assert set(grouped[1][1]['GateStatus']) == {'Report with caveat', 'Boundary'}
    assert set(grouped[2][1]['Evidence']) == {'screen not run', 'inference unqualified'}
    grouped[0][1].iloc[0, 2] = 'display copy changed'
    pd.testing.assert_frame_equal(source, before)
    assert app._group_guided_report_checks(pd.DataFrame()) == []
    assert app._group_guided_report_checks(source.drop(columns='GateStatus')) == []


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_review_questions_show_all_statuses_and_keep_source_details_optional(lang):
    def view(lang):
        from types import SimpleNamespace
        from unittest.mock import patch
        import pandas as pd
        import streamlit as st
        import streamlit_app as app

        st.session_state.setdefault('lang', lang)
        st.radio('Language', ['en', 'ja'], key='lang')
        st.session_state.setdefault('result', {'analysis_id': 'same-run',
            'config': {'method': 'MML', 'estimate_population_sd': True},
            'opt': SimpleNamespace(success=False)})
        gate = pd.DataFrame([
            ('Overall manuscript gate', 'Not ready', '', ''),
            ('Convergence', 'Not ready', 'termination failed', 'inspect termination'),
            ('Convergence and final interpretability', 'Not ready', 'inference hold', 'inspect qualification'),
            ('Reliability / separation', 'Report with caveat', 'rater evidence', 'rater action'),
            ('Reliability and separation', 'Boundary', 'person evidence', 'person action'),
            ('Bias / local interaction', 'Report with caveat', 'bias screen', 'bias screen action'),
            ('Bias / local interaction', 'Boundary', 'bias hold', 'bias hold action'),
            ('Future diagnostic', 'Future status', 'unknown evidence', 'unknown action'),
        ], columns=['GateArea', 'GateStatus', 'Evidence', 'ManuscriptAction'])
        with patch.object(app, 'mfrm_estimate', side_effect=AssertionError('UI must not refit')):
            app._render_guided_report_checks(gate, st.session_state['result'])
        st.session_state['expected'] = {
            'questions': [app.t(f'guided.review_topics.{key}.question')
                          for key in ('calibration', 'reliability', 'bias', 'other')],
            'statuses': [app.t(f'guided.review_status.{key}')
                         for key in ('hold', 'caveat', 'boundary', 'unknown')],
            'holds': [app.t(f'guided.review_hold_{key}') for key in ('termination', 'mml')],
        }

    at = AppTest.from_function(view, args=(lang,)).run(timeout=45)
    for current_lang in (lang, 'ja' if lang == 'en' else 'en'):
        if current_lang != lang:
            at.radio(key='lang').set_value(current_lang).run()
        assert not at.exception
        main_text = '\n'.join(node.value for node in at.main.children.values()
                              if node.type in ('markdown', 'caption'))
        expected = at.session_state['expected']
        for question in expected['questions']:
            assert main_text.count(question) == 1
        for message in expected['statuses'] + expected['holds']:
            assert message in main_text
        assert len(at.expander) == 4
        assert len({expander.label for expander in at.expander}) == 4
        assert all('{number}' not in expander.label and '{count}' not in expander.label
                   for expander in at.expander)
        assert all(not expander.proto.expanded for expander in at.expander)
        raw_text = '\n'.join(node.value for expander in at.expander for node in expander.markdown)
        for text in ('termination failed', 'inference hold', 'rater evidence', 'rater action',
                     'person evidence', 'person action', 'bias screen', 'bias screen action',
                     'bias hold', 'bias hold action', 'unknown evidence', 'unknown action'):
            assert text in raw_text
            assert text not in main_text
        assert 'guided.review_' not in main_text
        assert at.session_state['result']['analysis_id'] == 'same-run'
        if current_lang == 'ja':
            assert 'Reliability / separation' not in main_text
            assert 'Report with caveat' not in main_text
            assert 'Boundary' not in main_text
            assert '正常に終了したことを確認できていません' in main_text
            assert 'バイアスあり／なしの推論を保留' in main_text


def test_review_translation_covers_current_gate_topics_without_changing_gate():
    import json
    from types import SimpleNamespace

    result = {'opt': SimpleNamespace(success=True), 'config': {'method': 'JMLE'}}
    gate = app.build_publication_gate_summary(result, {}, {})
    before = gate.copy(deep=True)
    grouped = app._group_guided_report_checks(gate)
    assert grouped and all(topic != 'other' for topic, _ in grouped)
    assert sum(len(rows) for _, rows in grouped) == len(gate) - 1
    pd.testing.assert_frame_equal(gate, before)
    for lang in ('en', 'ja'):
        translations = json.loads((Path(app.__file__).parent / 'locales' / f'{lang}.json').read_text())['guided']
        assert set(translations['review_topics']) == set(app._GUIDED_REVIEW_TOPICS) | {'other'}
        assert all(value for topic in translations['review_topics'].values() for value in topic.values())
