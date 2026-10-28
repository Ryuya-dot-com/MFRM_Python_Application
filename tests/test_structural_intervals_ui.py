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
