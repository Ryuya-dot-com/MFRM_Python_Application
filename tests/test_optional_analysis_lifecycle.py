"""Optional inputs and explicit API arguments survive their consumer boundary."""
import ast
from copy import deepcopy

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from tests.test_report_run_facts import writing_run


def prediction_view():
    from unittest.mock import patch
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    st.radio('Language', ['en', 'ja'], key='lang')
    show = st.checkbox('Prediction panel', value=True, key='panel')
    st.session_state.setdefault('result', {
        'config': {'method': 'JMLE', 'model': 'RSM', 'facet_names': ['Rater']},
        'prep': {'data': pd.DataFrame({'Person': ['001'], 'Rater': ['R1'], 'Score': [1]})}})
    st.session_state.setdefault('calls', {'predict': 0, 'simulate': 0, 'refit': 0})
    def prediction(result, data, **kwargs):
        st.session_state['calls']['predict'] += 1
        return {'available': True, 'table': data.copy()}
    def simulation(result, **kwargs):
        st.session_state['calls']['simulate'] += 1
        return {'available': True, 'summary': pd.DataFrame({'MeanScore': [1.], 'MeanMinusObserved': [0.]})}
    def refit(result, **kwargs):
        st.session_state['calls']['refit'] += 1
        return {'available': True, 'summary': pd.DataFrame({'Converged': [True], 'PersonReliability': [.8]})}
    if show:
        with patch.object(app, 'render_fixed_density_assignment_sensitivity'), \
             patch.object(app, 'mfrm_estimate', side_effect=AssertionError('Unexpected fit')):
            app.show_prediction_simulation_section(st.session_state['result'], {}, core={
                'compute_fitted_prediction_table': lambda *a: pd.DataFrame(),
                'predict_mfrm_design': prediction, 'simulate_from_fitted_model': simulation,
                'simulate_refit_design': refit, 'evaluate_design_from_fitted': lambda *a, **k: {}})


def test_prediction_inputs_and_results_survive_navigation():
    at = AppTest.from_function(prediction_view).run(timeout=45)
    assert not at.exception
    at.text_area(key='new_design_prediction_text').set_value('Person,Rater\n001,R1')
    at.text_input(key='new_design_prediction_score_col').set_value('Rating')
    at.number_input(key='prediction_sim_replicates').set_value(2)
    at.number_input(key='prediction_sim_seed').set_value(17)
    at.checkbox(key='prediction_sim_row_draws').check()
    at.number_input(key='refit_sim_replicates').set_value(2)
    at.number_input(key='refit_sim_seed').set_value(19)
    at.slider(key='refit_sim_missing_rate').set_value(.15)
    at.button(key='run_new_design_prediction').click().run(timeout=45)
    at.button(key='run_fitted_simulation').click().run(timeout=45)
    at.button(key='run_refit_simulation').click().run(timeout=45)
    assert not at.exception and at.session_state['calls'] == {'predict': 1, 'simulate': 1, 'refit': 1}
    before = {name: deepcopy(at.session_state[f'mfrm_{name}_bundle'])
              for name in ('new_design_prediction', 'prediction_simulation', 'refit_simulation')}
    at.checkbox(key='panel').uncheck().run(timeout=45)
    at.radio(key='lang').set_value('ja').run(timeout=45)
    at.checkbox(key='panel').check().run(timeout=45)
    assert not at.exception
    assert at.text_area(key='new_design_prediction_text').value == 'Person,Rater\n001,R1'
    assert at.number_input(key='prediction_sim_replicates').value == 2
    assert at.number_input(key='refit_sim_seed').value == 19
    for name, bundle in before.items():
        current = at.session_state[f'mfrm_{name}_bundle']
        assert current['computation_settings'] == bundle['computation_settings']
        field = 'table' if name == 'new_design_prediction' else 'summary'
        pd.testing.assert_frame_equal(current[field], bundle[field])
    assert at.session_state['calls'] == {'predict': 1, 'simulate': 1, 'refit': 1}
    # Only changed scientific inputs invalidate the matching optional result.
    at.number_input(key='refit_sim_seed').set_value(23).run(timeout=45)
    assert not at.exception and 'mfrm_refit_simulation_bundle' not in at.session_state
    assert 'mfrm_new_design_prediction_bundle' in at.session_state
    at.session_state['result']['prep']['data'].loc[0, 'Score'] = 2
    at.run(timeout=45)
    assert not at.exception
    assert all(f'mfrm_{name}_bundle' not in at.session_state for name in before)
    assert at.session_state['calls'] == {'predict': 1, 'simulate': 1, 'refit': 1}


@pytest.mark.parametrize('pairs,mode,selected,expected', [
    ([('Person', 'Rater')], 'None', None, [('Person', 'Rater')]),
    ([], 'All facet pairs', None, []),
    (None, 'Selected pair', ['Person', 'Rater'], [('Person', 'Rater')]),
    (None, 'All facet pairs', None, [('Person', 'Rater')]),
])
def test_runner_consumes_explicit_pair_argument(writing_run, monkeypatch, pairs, mode, selected, expected):
    result, diagnostics = writing_run
    result['config'].update(bias_mode=mode, selected_bias_pair=selected)
    monkeypatch.setattr(app.st, 'session_state', {})
    source = app._generate_app_engine_runner_script(result, diagnostics, bias_pairs=pairs)
    namespace = {'__file__': '/tmp/mfrm_runner.py', '__name__': 'test_runner'}
    exec(compile(source, '<generated runner>', 'exec'), namespace)
    tree = ast.parse(source)
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'main')
    selection = next(n for n in main.body if isinstance(n, ast.If)
                     and any(isinstance(t, ast.Name) and t.id == 'bias_scan_pairs'
                             for child in ast.walk(n) if isinstance(child, ast.Assign) for t in child.targets))
    exec(compile(ast.Module(body=[selection], type_ignores=[]), '<pair selection>', 'exec'), namespace)
    assert namespace['bias_scan_pairs'] == expected


def dif_upload_view():
    from io import BytesIO
    from unittest.mock import patch
    import streamlit as st
    import streamlit_app as app
    from tests.test_classical_dif import _two_group_data
    st.radio('Language', ['en', 'ja'], key='lang')
    event = st.selectbox('Upload event', ['selected', 'removed'], key='upload_event')
    if 'result' not in st.session_state:
        data, groups = _two_group_data(n_person=30)
        st.session_state['result'] = {'prep': {'data': data, 'rating_min': 0, 'rating_max': 4},
                                     'config': {'facet_names': ['Rater', 'Criterion']}}
        st.session_state['groups_csv'] = groups.rename('Group').rename_axis('Person').to_csv().encode()
    def uploader(*args, **kwargs):
        file = None
        if event == 'selected':
            file = BytesIO(st.session_state['groups_csv']); file.name = 'groups.csv'
        if st.session_state.get('previous_event') != event:
            st.session_state[kwargs['key']] = file
            if kwargs.get('on_change'):
                kwargs['on_change']()
        st.session_state['previous_event'] = event
        return file
    with patch.object(app.st, 'file_uploader', uploader):
        app.show_classical_dif_section(st.session_state['result'], {})


def test_removing_dif_upload_does_not_silently_reuse_removed_groups():
    at = AppTest.from_function(dif_upload_view).run(timeout=45)
    assert not at.exception
    at.number_input(key='dif_min_n').set_value(5)
    at.button(key='dif_run').click().run(timeout=45)
    assert not at.exception and 'result' in at.session_state['mfrm_classical_dif']
    at.selectbox(key='upload_event').set_value('removed').run(timeout=45)
    assert not at.exception
    assert not any(button.key == 'dif_run' for button in at.button)
    assert 'result' not in at.session_state['mfrm_classical_dif']


def retained_upload_view():
    from io import BytesIO
    from unittest.mock import patch
    import streamlit as st
    import streamlit_app as app
    st.radio('Language', ['en', 'ja'], key='lang')
    source = st.selectbox('Analysis', ['A', 'B'], key='source')
    show = st.checkbox('Upload panel', value=True, key='panel')
    if st.session_state.get('saved_source') != source:
        st.session_state['saved'] = {}
        st.session_state['saved_source'] = source
    st.session_state.setdefault('visible_uploads', {})
    st.session_state.setdefault('change_count', 0)
    def changed():
        st.session_state['change_count'] += 1
    def uploader(*args, **kwargs):
        key = kwargs['key']
        if 'upload_event' in st.session_state:
            content = st.session_state.pop('upload_event')
            file = None
            if content is not None:
                file = BytesIO(content); file.name = 'same-name.csv'
            st.session_state['visible_uploads'][key] = file
            st.session_state[key] = file
            kwargs['on_change']()
        return st.session_state['visible_uploads'].get(key)
    if show:
        with patch.object(app.st, 'file_uploader', uploader):
            file = app._retained_analysis_upload('Groups', key='upload', source_key=source,
                saved=st.session_state['saved'], on_change=changed, type=['csv'])
        st.session_state['content'] = file.read() if file is not None else None
    else:
        # AppTest cannot drive the native upload widget; model its documented
        # cleanup without invoking the explicit file-change callback.
        st.session_state['visible_uploads'].clear()


@pytest.mark.parametrize('lang', ['en', 'ja'])
def test_retained_upload_reuses_bytes_and_allows_explicit_removal(lang):
    at = AppTest.from_function(retained_upload_view)
    at.session_state['lang'] = lang
    at.session_state['upload_event'] = b'Person,Group\n001,A\n002,B\n'
    at.run(timeout=45)
    assert not at.exception
    before = at.session_state['content']
    assert at.session_state['change_count'] == 1
    at.run(timeout=45)
    assert at.session_state['content'] == before  # Reading again is not EOF.
    at.checkbox(key='panel').uncheck().run(timeout=45)
    at.checkbox(key='panel').check().run(timeout=45)
    assert not at.exception and at.session_state['content'] == before
    assert at.session_state['change_count'] == 1
    assert 'same-name.csv' in at.caption[0].value
    at.button(key='upload::A::clear').click().run(timeout=45)
    assert not at.exception and at.session_state['content'] is None
    assert at.session_state['change_count'] == 2
    at.session_state['upload_event'] = b'Person,Group\n001,B\n002,A\n'
    at.run(timeout=45)
    assert at.session_state['content'] != before
    at.session_state['upload_event'] = None  # Native clear invokes on_change.
    at.run(timeout=45)
    assert not at.exception and at.session_state['content'] is None
    at.session_state['upload_event'] = before
    at.run(timeout=45)
    at.selectbox(key='source').set_value('B').run(timeout=45)
    assert not at.exception and at.session_state['content'] is None


@pytest.mark.parametrize('generator', ['engine', 'portable_python', 'portable_r'])
@pytest.mark.parametrize('ownership', ['foreign', 'current', 'missing'])
def test_runner_uses_only_the_requested_analyses_input_mapping(writing_run, monkeypatch, ownership, generator):
    result, diagnostics = writing_run
    stored = deepcopy(result)
    if ownership == 'foreign':
        stored['prep']['data'].loc[0, 'Score'] = 1
    output = {'person_col': 'OtherPerson', 'facet_cols': ['OtherRater'], 'score_col': 'OtherScore',
              '_maxit': 987, '_keep_original': True}
    if ownership != 'missing':
        output['result'] = stored
    monkeypatch.setattr(app.st, 'session_state', {'facets_mode_output': output})
    current = ownership == 'current'
    if generator == 'portable_r':
        source = app._generate_repro_r_script(result)
        person = 'OtherPerson' if current else 'Person'
        maximum = 987 if current else result['config']['maxit']
        assert f'PERSON_COL      <- "{person}"' in source
        assert f'MAXIT           <- {maximum}L' in source
        assert ('OtherRater' in source) == current
        return
    source = (app._generate_app_engine_runner_script(result, diagnostics) if generator == 'engine'
              else app._generate_repro_python_script(result))
    compile(source, '<generated runner>', 'exec')
    # Inspect the generated constants, without loading data or fitting a model.
    namespace = {}
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in {'PERSON_COL', 'FACET_COLS', 'SCORE_COL', 'MAXIT'}:
                namespace[node.targets[0].id] = ast.literal_eval(node.value)
    assert namespace['PERSON_COL'] == ('OtherPerson' if current else 'Person')
    assert namespace['FACET_COLS'] == (['OtherRater'] if current else ['Rater'])
    assert namespace['SCORE_COL'] == ('OtherScore' if current else 'Score')
    assert namespace['MAXIT'] == (987 if current else result['config']['maxit'])
