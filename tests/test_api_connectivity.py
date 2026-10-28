"""Retained analyses own optional results and the settings used to report them."""
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
from tests.test_first_read_qualification import technical_bias
from tests.test_report_run_facts import writing_run


@pytest.mark.parametrize('mutation', ['data', 'parameters', 'optimizer', 'anchor'])
def test_optional_cache_distinguishes_same_shape_analyses(writing_run, mutation):
    result, _ = writing_run
    result['params'] = {'theta': np.array([-.3, .3])}
    result['opt'].x = np.array([-.3, .3])
    changed = deepcopy(result)
    if mutation == 'data':
        changed['prep']['data'].loc[0, 'Score'] = 1
    elif mutation == 'parameters':
        changed['params']['theta'][0] = -.4
    elif mutation == 'optimizer':
        changed['opt'].x[0] = -.4
    else:
        changed['config']['theta_spec'] = {'anchors': {'P1': -.4}}
    key = app._model_choice_cache_key(result)
    assert key == app._model_choice_cache_key(deepcopy(result))
    assert key != app._model_choice_cache_key(changed)
    with np.printoptions(precision=1), pd.option_context('display.max_rows', 1):
        assert key == app._model_choice_cache_key(result)


@pytest.mark.parametrize('ownership', ['current', 'foreign', 'legacy', 'missing'])
def test_shared_export_rejects_foreign_optional_tables(writing_run, monkeypatch, ownership):
    result, diagnostics = writing_run
    other = deepcopy(result)
    other['prep']['data'].loc[0, 'Score'] = 1
    state = {}
    for name, table_key in [('new_design_prediction', 'table'), ('refit_simulation', 'summary')]:
        state[f'mfrm_{name}_bundle'] = {table_key: pd.DataFrame({'Source': ['owned']}),
                                      'computation_settings': {'replicates': 2}}
        if ownership != 'missing':
            state[f'mfrm_{name}_meta'] = ((app._model_choice_cache_key(result if ownership == 'current' else other),)
                                         if ownership != 'legacy' else ('JMLE', 'PCM', 6, -10))
    monkeypatch.setattr(app.st, 'session_state', state)
    frames, _ = app.collect_download_frames(result, diagnostics, {}, pd.DataFrame(), pd.DataFrame(), public_export_mode=False)
    for key in ('new_design_predictions', 'refit_simulation_summary'):
        assert (key in frames) == (ownership == 'current')
    assert bool(app._current_prediction_bundle(result, 'mfrm_refit_simulation_bundle')) == (ownership == 'current')


def dif_view():
    from unittest.mock import patch
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    from tests.test_classical_dif import _two_group_data
    st.radio('Language', ['en', 'ja'], key='lang')
    show = st.checkbox('DIF panel', value=True, key='show_dif')
    if 'fit' not in st.session_state:
        data, groups = _two_group_data(n_person=30)
        st.session_state['fit'] = {'prep': {'data': data, 'rating_min': 0, 'rating_max': 4},
            'config': {'facet_names': ['Rater', 'Task', 'Criterion']},
            'population': {'person_data': groups.rename('Group').rename_axis('Person').reset_index()}}
        st.session_state['calls'] = 0
    compute = app.analyze_classical_dif
    def counted(*args, **kwargs):
        st.session_state['calls'] += 1
        return compute(*args, **kwargs)
    if show:
        with patch.object(app, 'analyze_classical_dif', counted), \
             patch.object(app, 'mfrm_estimate', side_effect=AssertionError('Unexpected MFRM fit')):
            app.show_classical_dif_section(st.session_state['fit'], {})


def has_dif_download(at):
    return any(button.proto.id.endswith('-dl_dif_screen') for button in at.get('download_button'))


def test_dif_survives_language_navigation_and_invalidates_changed_inputs():
    at = AppTest.from_function(dif_view).run(timeout=45)
    assert not at.exception
    at.number_input(key='dif_min_n').set_value(5)
    at.button(key='dif_run').click().run(timeout=45)
    assert not at.exception and has_dif_download(at)
    saved = deepcopy(at.session_state['mfrm_classical_dif']['result'])
    for target in ('language', 'away', 'back', 'rerun'):
        if target == 'language':
            at.radio(key='lang').set_value('ja')
        elif target in ('away', 'back'):
            at.checkbox(key='show_dif').set_value(target == 'back')
        at.run(timeout=45)
        assert not at.exception and at.session_state['calls'] == 1
        if target != 'away':
            assert has_dif_download(at) and at.number_input(key='dif_min_n').value == 5
            retained = at.session_state['mfrm_classical_dif']['result']
            for key in ('table', 'audit'):
                pd.testing.assert_frame_equal(retained[key], saved[key])
            assert retained['caveats'] == saved['caveats'] and retained['settings'] == saved['settings']
    at.number_input(key='dif_alpha').set_value(.1).run(timeout=45)
    assert not at.exception and not has_dif_download(at) and at.session_state['calls'] == 1
    at.button(key='dif_run').click().run(timeout=45)
    assert not at.exception and has_dif_download(at) and at.session_state['calls'] == 2
    # Same values with different Person associations must not reuse the old screen.
    groups = at.session_state['fit']['population']['person_data']
    groups.loc[[0, 1], 'Person'] = groups.loc[[1, 0], 'Person'].to_numpy()
    at.run(timeout=45)
    assert not at.exception and not has_dif_download(at) and at.session_state['calls'] == 2
    at.button(key='dif_run').click().run(timeout=45)
    assert not at.exception and has_dif_download(at)
    at.session_state['fit']['prep']['data'].loc[0, 'Score'] += 1
    at.run(timeout=45)
    assert not at.exception and not has_dif_download(at) and at.session_state['calls'] == 3


def configured_bias():
    bias = technical_bias()
    bias['Person x Rater']['screening_settings'] = {'alpha': .1, 'min_n': 100, 'practical_logit': 1.}
    bias['Rater x Task'] = deepcopy(technical_bias()['Person x Rater'])
    bias['Rater x Task']['table']['FacetA'] = 'Rater'
    bias['Rater x Task']['table']['FacetB'] = 'Task'
    return bias


def test_bias_pair_settings_reach_summary_audits_and_common_csv(writing_run, monkeypatch):
    monkeypatch.setattr(app.st, 'session_state', {})
    result, diagnostics = writing_run
    bias = configured_bias()
    raw = deepcopy(bias['Person x Rater']['table'])
    custom = app.build_dff_bias_screening_table(bias['Person x Rater'])
    default = app.build_dff_bias_screening_table(bias['Person x Rater'], alpha=.05, min_n=5, practical_logit=.5)
    assert custom.Flag.sum() == 1 and default.Flag.sum() == 0
    summary = app.summarize_dff_bias_screening(custom).set_index('Metric')
    assert 'Holm p < .1' in summary.index and 'BH q < .1' in summary.index
    audit = app.build_bias_inference_audit(bias, result, diagnostics).set_index('FacetPair')
    assert audit.loc['Person x Rater', 'ScreenMinN'] == 100
    assert audit.loc['Rater x Task', 'ScreenMinN'] == 5
    assert not audit.PublicConclusionAllowed.any()
    frames, _ = app.collect_download_frames(result, diagnostics, {}, pd.DataFrame(), pd.DataFrame(),
                                           public_export_mode=False, all_bias_results=bias)
    exported = frames[app._frame_bundle.safe_frame_key('dff_bias', 'Person x Rater', max_label_length=40)]
    pd.testing.assert_frame_equal(exported, custom)
    notes = app._publication_review_notes(result, diagnostics, bias)
    assert '1 of 2' in notes and 'minimum ratings per cell=100' in notes and 'minimum ratings per cell=5' in notes
    assert 'Bias inference remains withheld' in notes
    pd.testing.assert_frame_equal(raw, bias['Person x Rater']['table'])


@pytest.mark.parametrize('kind', ['pdf', 'word', 'html'])
@pytest.mark.parametrize('details', [False, True])
def test_bias_settings_reach_each_report_format(writing_run, monkeypatch, kind, details):
    from reportlab.platypus import Paragraph, SimpleDocTemplate
    result, diagnostics = writing_run
    monkeypatch.setattr(app, '_publication_figure_payloads', lambda *args: [])
    monkeypatch.setattr(app, 'generate_method_appendix_text', lambda *args: 'Method appendix')
    monkeypatch.setattr(app, 'generate_manuscript_reporting_template', lambda *args: 'Writing appendix')
    paragraphs = []
    original = SimpleDocTemplate.build
    def capture(self, story, *args, **kwargs):
        paragraphs.extend(p.getPlainText() for p in story if isinstance(p, Paragraph))
        return original(self, story, *args, **kwargs)
    monkeypatch.setattr(SimpleDocTemplate, 'build', capture)
    payload = getattr(app, f'build_publication_{kind}_bytes')(result, diagnostics, configured_bias(), include_details=details)
    if kind == 'word':
        text = '\n'.join(p.text for p in Document(BytesIO(payload)).paragraphs)
    elif kind == 'html':
        text = '\n'.join(ET.fromstring(payload.decode().split('<!DOCTYPE html>\n', 1)[1]).itertext())
    else:
        assert payload.startswith(b'%PDF')
        text = '\n'.join(paragraphs)
    assert 'minimum ratings per cell=100' in text and 'minimum ratings per cell=5' in text and 'alpha=0.1' in text
    assert 'Bias inference remains withheld' in text


def bias_view():
    from io import BytesIO
    from unittest.mock import patch
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    from tests.test_api_connectivity import configured_bias
    st.radio('Language', ['en', 'ja'], key='lang')
    visible = st.checkbox('Bias panel', value=True, key='show_bias')
    st.session_state.setdefault('bias', configured_bias())
    download = st.download_button
    def capture(*args, **kwargs):
        if kwargs.get('file_name') == 'mfrm_dff_bias_screening.csv':
            st.session_state['screen_csv'] = pd.read_csv(BytesIO(kwargs['data']))
        return download(*args, **kwargs)
    if visible:
        with patch.object(app, '_draw_bias_heatmap'), patch.object(app, '_offer_fig_download'), \
             patch.object(app.st, 'download_button', capture), \
             patch.object(app, 'mfrm_estimate', side_effect=AssertionError('Unexpected fit')):
            app.show_bias_section(None, {}, st.session_state['bias'])
    st.session_state['notes'] = app._publication_review_notes({'config': {}}, {}, st.session_state['bias'])


def test_bias_widget_changes_persist_per_pair_and_feed_report():
    at = AppTest.from_function(bias_view).run(timeout=45)
    assert not at.exception
    assert at.number_input(key='dff_min_n_person_x_rater').value == 100
    at.number_input(key='dff_alpha_person_x_rater').set_value(.15)
    at.number_input(key='dff_min_n_person_x_rater').set_value(25)
    at.number_input(key='dff_practical_person_x_rater').set_value(.8).run(timeout=45)
    assert not at.exception
    for action in ('language', 'away', 'back', 'pair', 'restore_pair'):
        if action == 'language':
            at.radio(key='lang').set_value('ja')
        elif action in ('away', 'back'):
            at.checkbox(key='show_bias').set_value(action == 'back')
        else:
            at.selectbox(key='bias_pair_selector').set_value('Rater x Task' if action == 'pair' else 'Person x Rater')
        at.run(timeout=45)
        assert not at.exception
        settings = at.session_state['bias']['Person x Rater']['screening_settings']
        assert settings == {'alpha': .15, 'min_n': 25, 'practical_logit': .8}
        assert 'minimum ratings per cell=25' in at.session_state['notes'] and 'alpha=0.15' in at.session_state['notes']
        if action not in ('away', 'pair'):
            assert at.number_input(key='dff_min_n_person_x_rater').value == 25
            assert at.session_state['screen_csv']['ScreenMinN'].eq(25).all()
        if action == 'pair':
            assert at.number_input(key='dff_min_n_rater_x_task').value == 5


def test_runner_carries_bias_settings_and_ignores_foreign_simulation_settings(writing_run, monkeypatch):
    import ast
    result, diagnostics = writing_run
    state = {'prediction_sim_replicates': 987, 'refit_sim_replicates': 49,
             'mfrm_prediction_simulation_bundle': {'computation_settings': {'replicates': 456}},
             'mfrm_prediction_simulation_meta': ('foreign',)}
    monkeypatch.setattr(app.st, 'session_state', state)
    bias = configured_bias()
    assets = app.python_native_reproducibility_assets(result, diagnostics, bias)
    source = assets['mfrm_app_engine_runner.py']
    namespace = {'__name__': 'test_runner', '__file__': '/tmp/mfrm_app_engine_runner.py'}
    exec(compile(source, '<generated runner>', 'exec'), namespace)
    assert namespace['SIMULATION_REPLICATES'] == 100
    assert namespace['REFIT_SIMULATION_REPLICATES'] == 3
    assert namespace['BIAS_SCREENING_SETTINGS']['Person x Rater']['min_n'] == 100
    # Execute the generated bias handoff with recorded cells, without refitting.
    module = ast.parse(source)
    main = next(node for node in module.body if isinstance(node, ast.FunctionDef) and node.name == 'main')
    bias_loop = next(node for node in main.body if isinstance(node, ast.If) and 'bias_scan_pairs and' in ast.unparse(node.test))
    selection = next(node for node in main.body if isinstance(node, ast.If) and 'BIAS_PAIRS is not None' in ast.unparse(node.test))
    frames = {}
    namespace.update(result=result, diagnostics=diagnostics,
                     app=SimpleNamespace(estimate_bias_interaction=lambda *a, **k: deepcopy(bias['Person x Rater']),
                                         build_dff_bias_screening_table=app.build_dff_bias_screening_table),
                     write_csv=lambda name, value: frames.update({name: value}))
    exec(compile(ast.Module(body=[selection, bias_loop], type_ignores=[]), '<bias export>', 'exec'), namespace)
    assert frames['dff_bias_Person_Rater'].ScreenMinN.eq(100).all()
    assert frames['dff_bias_Person_Rater'].Flag.all()
