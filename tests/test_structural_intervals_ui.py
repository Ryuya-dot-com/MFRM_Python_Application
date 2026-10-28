"""Real fitted result -> selected difference -> common table/figure/export/replay."""
import json

import pytest
from streamlit.testing.v1 import AppTest


def view(lang, model):
    from unittest.mock import patch
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    from mfrm_app import engine
    if 'fit' not in st.session_state:
        data = pd.DataFrame([dict(Person=f'{p:03d}', Rater=f'R{r}', Task=f'T{t}', Score=(p+2*r+t)%3)
                            for p in range(30) for r in range(3) for t in range(2)])
        st.session_state['fit'] = engine.mfrm_estimate(data,'Person',['Rater','Task'],'Score',
            method='MML', model=model, step_facet='Task', quad_points=61, maxit=100)
        st.session_state['lang'] = lang
        st.session_state['evaluations'] = 0
    original = engine.mfrm_loglik_mml_value_grad
    def evaluate(*args, **kwargs):
        st.session_state['evaluations'] += 1
        return original(*args, **kwargs)
    def no_refit(*args, **kwargs):
        raise AssertionError('Displaying/computing intervals must not refit')
    with patch.object(engine,'mfrm_loglik_mml_value_grad',evaluate), patch.object(engine,'minimize',no_refit):
        app._render_structural_intervals(st.session_state['fit'], key_prefix='test')


@pytest.mark.parametrize('lang', ['en','ja'])
@pytest.mark.parametrize('model', ['RSM','PCM'])
def test_selected_intervals_survive_language_change_and_saved_record_restore(lang, model):
    import streamlit_app as app
    from mfrm_app import structural_intervals as intervals
    at = AppTest.from_function(view, args=(lang,model)).run(timeout=60)
    assert not at.exception and at.session_state['evaluations'] == 0
    prefix = 'test_structural_intervals'
    at.checkbox(key=prefix+'_difference').check().run()
    at.selectbox(key=prefix+'_second').set_value('R2').run()
    at.button(key=prefix+'_calculate').click().run()
    assert not at.exception and at.session_state['evaluations'] == 2
    result = at.session_state['fit']
    record = result['_structural_intervals']
    assert record['table'][-1]['Target'] == 'R0 − R2'
    assert all(row['CIEligible'] and not row['CoverageQualified'] for row in record['table'])
    assert len(at.get('plotly_chart')) == 1 and len(at.get('download_button')) == 3
    fields = app._load_locale(lang)['structural_intervals']
    assert fields['CIEligible'] in at.dataframe[0].value
    assert any(fields['available'] == x.value for x in at.info)
    # JSON restoration uses the same validator as the uploader; no optimizer or likelihood call.
    restored = json.loads(intervals.interval_json(record))
    result['_structural_intervals'] = restored
    at.session_state['lang'] = 'ja' if lang == 'en' else 'en'
    at.run()
    assert not at.exception and at.session_state['evaluations'] == 2
    assert result['_structural_intervals'] == record
    assert len(at.get('download_button')) == 3
    # An actual score edit invalidates the saved intervals without an automatic calculation.
    result['prep']['data'].iloc[0, result['prep']['data'].columns.get_loc('Score')] = 2
    at.run()
    assert not at.exception and at.warning and at.session_state['evaluations'] == 2
    assert not at.get('download_button') and not at.get('plotly_chart')


@pytest.mark.parametrize('store_snapshot', [True, False])
def test_interval_work_is_retained_only_in_the_matching_opted_in_history(store_snapshot):
    import streamlit_app as app
    from copy import deepcopy
    from unittest.mock import patch
    from mfrm_app import engine

    at = AppTest.from_function(view, args=('en', 'RSM')).run(timeout=60)
    fit = at.session_state['fit']
    output = dict(result=fit)
    # Use the real history/storage functions against this session, including a
    # newer run with the same fitted values: equality must not select the run.
    with patch.object(app.st, 'session_state', at.session_state):
        at.session_state['facets_mode_output'] = output
        at.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY] = store_snapshot
        first = app.record_run_in_history(output=output, elapsed_sec=1.)
        second = app.record_run_in_history(output=deepcopy(output), elapsed_sec=1.)
    at.button(key='test_structural_intervals_calculate').click().run()
    assert not at.exception
    record = fit['_structural_intervals']
    if not store_snapshot:
        assert first['output_snapshot'] == second['output_snapshot'] == {}
        return
    assert first['output_snapshot']['result']['_structural_intervals'] == record
    assert first['output_snapshot']['result']['_structural_intervals'] is not record
    assert '_structural_intervals' not in second['output_snapshot']['result']
    with patch.object(app.st, 'session_state', at.session_state), \
         patch.object(engine, 'minimize', side_effect=AssertionError('Unexpected refit')), \
         patch.object(engine, 'mfrm_loglik_mml_value_grad', side_effect=AssertionError('Unexpected evaluation')):
        assert app.restore_run_from_history(first['run_id'])
        restored = at.session_state['facets_mode_output']['result']
        saved = json.loads(json.dumps(record))
        app._structural_intervals.read_structural_intervals(restored, saved)
        app._store_structural_interval_record(restored, saved)
        assert app.restore_run_from_history(first['run_id'])
        assert at.session_state['facets_mode_output']['result']['_structural_intervals'] == record
        assert '_structural_intervals' not in second['output_snapshot']['result']
