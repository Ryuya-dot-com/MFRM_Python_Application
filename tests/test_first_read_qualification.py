"""First-read and report claims follow retained qualification, not large values."""
from copy import deepcopy
from io import BytesIO
from types import SimpleNamespace
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
import pytest
from docx import Document
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_reliability_input_contract import measures
from tests.test_report_run_facts import writing_run


def result_record():
    return {'opt': SimpleNamespace(success=True),
            'config': {'method': 'JMLE', 'model': 'RSM', 'facet_names': ['Rater', 'Task']}}


def reliability_diagnostics():
    raw = measures(eap=False)
    raw['Facet'] = 'Person'
    return {'measures': raw, 'reliability': app.calc_reliability(raw), 'pca_enabled': False}


def unavailable_pca():
    obs = pd.DataFrame({'Person': ['P1', 'P2', 'P3'], 'Rater': ['R1']*3,
                        'Task': ['T1']*3, 'StdResidual': [-1., 0., 1.]})
    assert app.compute_pca_overall(obs, ['Rater', 'Task']) is None
    return {'pca_enabled': True, 'pca': None,
            'pca_reason': app.diagnose_pca_skip_reason(obs, ['Rater', 'Task'])}


def technical_bias(*, skip=False, unavailable=False):
    # Presentation fixture with actual current audit metadata, not a fitted study.
    table = pd.DataFrame({'FacetA': ['Person'], 'FacetB': ['Rater'],
                          'FacetA_Level': ['P1'], 'FacetB_Level': ['R1'],
                          'Bias Size': [.1], 'S.E.': [.2], 't': [.5], 'ObsN': [20],
                          'LR Prob.': [.6], 'Cell Status': ['interior'],
                          'Primary P Method': ['unit_weight_conditional_lr_v1'],
                          'Reference Status': ['unit_weight_conditional_lr_candidate']})
    if unavailable:
        table['Reference Status'] = 'not_available'
    bias = {'Person x Rater': {'table': table}}
    if skip:
        bias['Rater x Task'] = {'_skip_reason': 'No usable cells'}
    return bias


def check_row(result, diagnostics, check, bias=None, **kwargs):
    return next(row for row in app.build_first_read_guide_rows(result, diagnostics, bias, **kwargs)
                if row['Check'] == check)


@pytest.mark.parametrize('case,reason', [
    ('stale', 'Saved values do not match'),
    ('legacy', 'saved calculation could not be verified'),
    ('missing_se', 'estimates or standard errors are missing'),
    ('basis', 'standard-error method or its qualification'),
])
def test_unqualified_reliability_never_authorizes_ranking(case, reason):
    diagnostics = reliability_diagnostics()
    if case == 'stale':
        diagnostics['reliability'].loc[0, 'Reliability'] = .9999
    elif case == 'legacy':
        diagnostics.pop('measures')
        diagnostics['reliability'].drop(columns='ReliabilityCalculationVersion', inplace=True)
    elif case == 'missing_se':
        diagnostics['measures'].loc[0, 'SE'] = np.nan
    else:
        diagnostics['measures'].drop(columns='SE_Method', inplace=True)
    original = deepcopy(diagnostics)
    row = check_row(result_record(), diagnostics, '3. Reliability / separation')
    assert row['Status'] == 'Review' and reason in row['What it means']
    assert 'Person reliability is 1.00' not in row['What it means']
    assert 'Person reliability is 0.96' not in row['What it means']
    assert 'Use person ordering claims' not in row['Next action']
    assert '3. Reliability / separation (Review)' in app._publication_review_notes(result_record(), diagnostics)
    for name, value in original.items():
        if isinstance(value, pd.DataFrame):
            pd.testing.assert_frame_equal(diagnostics[name], value)


@pytest.mark.parametrize('held_facet', ['Person', 'Rater'])
def test_one_held_facet_preserves_valid_neighbor(held_facet):
    diagnostics = reliability_diagnostics()
    diagnostics['measures'] = pd.concat([diagnostics['measures'], measures(eap=False)], ignore_index=True)
    diagnostics['reliability'] = app.calc_reliability(diagnostics['measures'])
    saved = diagnostics['reliability']
    saved.loc[saved.Facet.eq(held_facet), 'Reliability'] = .9999
    row = check_row(result_record(), diagnostics, '3. Reliability / separation')
    valid_facet = 'Rater' if held_facet == 'Person' else 'Person'
    assert row['Status'] == 'Review'
    assert f'cannot be interpreted for {held_facet}' in row['What it means']
    assert f'Descriptive interpretation remains available for: {valid_facet}' in row['What it means']
    table, _ = app._publication_reliability_table(result_record(), diagnostics)
    assert table.set_index('Facet').loc[valid_facet, 'Use'] == 'Diagnostic'
    assert table.set_index('Facet').loc[held_facet, 'Use'] == 'Withheld'


def test_valid_reliability_retains_value_and_success():
    diagnostics = reliability_diagnostics()
    row = check_row(result_record(), diagnostics, '3. Reliability / separation')
    assert row['Status'] == 'OK' and '0.96' in row['What it means']
    assert diagnostics['reliability'].Reliability.iloc[0] == pytest.approx(.96)


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_pca_unavailable_keeps_actual_localizable_reason(monkeypatch, lang):
    monkeypatch.setitem(app.st.session_state, 'lang', lang)
    diagnostics = unavailable_pca()
    row = check_row(result_record(), diagnostics, '4. Dimensionality', localize=True)
    assert row['Status'] == 'Review'
    assert app._render_reason(diagnostics['pca_reason']) in row['What it means']
    assert ('could not be assessed' if lang == 'en' else '結果を確認できません') in row['What it means']
    report_row = check_row(result_record(), diagnostics, '4. Dimensionality')
    assert 'item-combination' in report_row['What it means']
    assert '項目' not in report_row['What it means']


@pytest.mark.parametrize('stability', [None, '', 'Unknown', 'Review', 'Ready', 'Stable screen'])
def test_pca_stability_requires_affirmative_evidence(stability):
    diagnostics = {'pca_enabled': True, 'pca': {'eigenvalues': [1.5]}}
    if stability is not None:
        diagnostics['pca']['stability_table'] = pd.DataFrame([{'PCAStabilityStatus': stability}])
    row = check_row(result_record(), diagnostics, '4. Dimensionality')
    ready = stability in {'Ready', 'Stable screen'}
    assert (row['Status'] == 'OK') == ready
    final = app.build_final_report_readiness(result_record(), diagnostics)
    final_row = final.loc[final.Check.eq('Dimensionality screen')].iloc[0]
    assert (final_row.Status == 'Ready') == ready
    if not ready and stability != 'Review':
        assert 'Stability has not been established' in row['What it means']
        assert 'stable enough' not in row['Next action']


@pytest.mark.parametrize('eigenvalues', [None, [], [np.nan], [np.inf]])
def test_unusable_pca_values_do_not_become_a_screen(eigenvalues):
    row = check_row(result_record(), {'pca_enabled': True, 'pca': {'eigenvalues': eigenvalues}}, '4. Dimensionality')
    assert row['Status'] == 'Review' and 'could not be assessed' in row['What it means']
    assert 'nan' not in row['What it means'] and 'inf' not in row['What it means']


def test_disabled_pca_stays_distinct_from_unavailable():
    row = check_row(result_record(), {'pca_enabled': False}, '4. Dimensionality')
    assert row['Status'] == 'Skipped' and 'turned off' in row['What it means']


@pytest.mark.parametrize('skip', [False, True])
def test_zero_bias_flags_keep_inference_hold_and_unassessed_pairs(skip):
    result, bias = result_record(), technical_bias(skip=skip)
    audit = app.build_bias_inference_audit(bias, result, {})
    assert audit.FlaggedCells.sum() == 0 and audit.CellsScreened.sum() == 1
    row = check_row(result, {}, '5. Bias interactions', bias)
    assert row['Status'] == 'Review'
    assert '0 of 1' in row['What it means'] and 'Bias inference remains withheld' in row['What it means']
    if skip:
        assert 'Rater x Task' in row['What it means'] and 'No usable cells' in row['What it means']
    assert '5. Bias interactions (Review)' in app._publication_review_notes(result, {}, bias)
    plan = app.build_guided_action_plan(result, {}, bias)
    bias_row = plan.loc[plan.Check.eq('5. Bias interactions')].iloc[0]
    assert app._guided_focus_class(bias_row) != 'hard_stop'


def test_unavailable_bias_cells_are_not_disguised_as_no_flags():
    row = check_row(result_record(), {}, '5. Bias interactions', technical_bias(unavailable=True))
    assert row['Status'] != 'OK'
    assert 'Inference was unavailable for 1 cell(s)' in row['What it means']


@pytest.mark.parametrize('include_details', [False, True])
@pytest.mark.parametrize('kind', ['pdf', 'word', 'html'])
def test_all_report_formats_retain_the_three_qualification_failures(writing_run, monkeypatch, kind, include_details):
    from reportlab.platypus import Paragraph, SimpleDocTemplate
    result, _ = writing_run
    diagnostics = reliability_diagnostics()
    diagnostics['reliability'].loc[0, 'Reliability'] = .9999
    diagnostics.update(unavailable_pca())
    monkeypatch.setitem(app.st.session_state, 'lang', 'ja')
    monkeypatch.setattr(app, '_publication_figure_payloads', lambda *a: [])
    monkeypatch.setattr(app, 'generate_method_appendix_text', lambda *a: 'Method appendix')
    monkeypatch.setattr(app, 'generate_manuscript_reporting_template', lambda *a: 'Writing appendix')
    pdf_text = []
    original_build = SimpleDocTemplate.build
    def capture(self, story, *args, **kwargs):
        pdf_text.extend(part.getPlainText() for part in story if isinstance(part, Paragraph))
        return original_build(self, story, *args, **kwargs)
    monkeypatch.setattr(SimpleDocTemplate, 'build', capture)
    payload = getattr(app, f'build_publication_{kind}_bytes')(
        result, diagnostics, technical_bias(skip=True), include_details=include_details)
    if kind == 'word':
        text = '\n'.join(p.text for p in Document(BytesIO(payload)).paragraphs)
    elif kind == 'html':
        text = '\n'.join(ET.fromstring(payload.decode().split('<!DOCTYPE html>\n', 1)[1]).itertext())
    else:
        assert payload.startswith(b'%PDF')
        text = '\n'.join(pdf_text)
    notes = text.split('Computation and interpretation checks', 1)[1].split('Results tables', 1)[0]
    for phrase in ('Saved values do not match', 'could not be assessed', 'item-combination',
                   'Bias inference remains withheld', 'Rater x Task', 'No usable cells'):
        assert phrase in notes
    assert 'Person reliability is 1.00' not in notes and '項目' not in notes


def test_rendered_first_read_switches_language_without_estimation():
    def view():
        from unittest.mock import patch
        import streamlit as st
        import streamlit_app as app
        from tests.test_first_read_qualification import result_record, reliability_diagnostics, unavailable_pca, technical_bias
        st.radio('Language', ['en', 'ja'], key='lang')
        diagnostics = reliability_diagnostics()
        diagnostics['reliability'].loc[0, 'Reliability'] = .9999
        diagnostics.update(unavailable_pca())
        with patch.object(app, 'mfrm_estimate', side_effect=AssertionError('unexpected refit')):
            app._show_first_read_guide(result_record(), diagnostics, technical_bias(skip=True))
    at = AppTest.from_function(view).run(timeout=45)
    assert not at.exception
    for lang in ('en', 'ja'):
        if lang == 'ja':
            at.radio(key='lang').set_value(lang).run(timeout=45)
        assert not at.exception
        text = '\n'.join(str(item.value) for kind in ('markdown', 'info', 'caption', 'warning') for item in getattr(at, kind))
        text += '\n' + '\n'.join(item.value.to_string(index=False) for item in at.dataframe)
        for phrase in (('Saved values do not match', 'could not be assessed', 'Bias inference remains withheld')
                       if lang == 'en' else ('保存値が', '結果を確認できません', 'バイアスに関する推論は保留中')):
            assert phrase in text
