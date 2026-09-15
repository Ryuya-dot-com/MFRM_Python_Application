"""Returned-point checking is opt-in, current-fit-bound and never a refit."""
import json

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app


def view(lang):
    from unittest.mock import patch
    import pandas as pd
    import streamlit as st
    import streamlit_app as app
    if 'fit' not in st.session_state:
        data=pd.DataFrame([dict(Person=f'{p:03d}',Rater=f'R{r}',Task=f'T{t}',Criterion=f'C{c}',Score=(p+2*r+t+c)%4)
            for p in range(8) for r in range(2) for t in range(2) for c in range(2)])
        st.session_state['fit']=app.mfrm_estimate(data=data,person_col='Person',facet_cols=['Rater','Task','Criterion'],
            score_col='Score',model='PCM',method='MML',step_facet='Criterion',noncenter_facet='Criterion',
            population_prior_sd=1.1,quad_points=31,maxit=40)
        st.session_state['lang']=lang
        st.session_state['evaluations']=0
    fit=st.session_state['fit']
    original=app._analysis_input_assets(fit)
    actual=app.mfrm_loglik_mml_value_grad
    def evaluate(*args,**kwargs):
        st.session_state['evaluations']+=1
        value,gradient=actual(*args,**kwargs)
        if st.session_state.get('flat'):
            gradient=gradient*0. # UI-only seam; mathematical checks use the actual gradient.
        return value,gradient
    def no_refit(*args,**kwargs):
        raise AssertionError('A returned-point check must not optimize or launch an engine')
    with patch.object(app,'mfrm_loglik_mml_value_grad',evaluate), patch.object(app,'minimize',no_refit), \
         patch.object(app._pcm_native,'run_process',no_refit):
        app.show_convergence_section(fit,key_prefix='test')
        if st.session_state.get('two_contexts'):
            app.show_convergence_section(fit,key_prefix='report')
    assert app._analysis_input_assets(fit)==original


@pytest.mark.parametrize('lang', ['en','ja'])
def test_compact_ui_waits_for_click_and_keeps_integration_unassessed(lang):
    at=AppTest.from_function(view,args=(lang,)).run(timeout=60)
    assert not at.exception and not at.success and at.session_state['evaluations']==0
    assert len(at.expander)==1 and not at.expander[0].proto.expanded
    assert len(at.button)==1 and not at.get('download_button')
    # All technical tables are inside the initially closed details panel.
    assert len(at.expander[0].dataframe)==len(at.dataframe)
    at.button(key='test_returned_point_check').click().run()
    assert not at.exception and at.session_state['evaluations']==1 and not at.success
    record=at.session_state['_fixed_sd_returned_point_check']
    assert record['stationarity_pass'] is None and record['quadrature_sensitivity_pass'] is None
    assert record['inference_ready'] is False
    text='\n'.join(m.value for m in at.markdown)
    assert app._load_locale(lang)['estimation_subsections']['returned_point']['integration'] in text
    assert ('未確認' if lang=='ja' else 'Not checked') in text
    fit=at.session_state['fit'];prepared=app.prepare_fixed_sd_pcm_check(fit)
    actual=app.evaluate_fixed_sd_pcm_check(prepared)
    np.testing.assert_array_equal(record['gradient'],actual['gradient'])
    assert app._analysis_input_assets(fit)['analysis_snapshot.json']==app.build_cross_engine_validation_bundle(fit)['analysis_snapshot.json']
    at.session_state['lang']='ja' if lang=='en' else 'en'
    at.run()
    assert not at.exception and at.session_state['evaluations']==1
    at.session_state['flat']=True
    at.button(key='test_returned_point_check').click().run()
    assert not at.exception and not at.success and at.session_state['evaluations']==2
    assert at.session_state['_fixed_sd_returned_point_check']['gradient_supnorm']==0.
    assert at.info and not at.warning
    # Alter the same result object while leaving its recorded run fingerprints untouched.
    fit=at.session_state['fit'];fit['config']['quad_points']=61
    at.run()
    assert not at.exception and at.session_state['evaluations']==2
    assert '_fixed_sd_returned_point_check' not in at.session_state
    at.button(key='test_returned_point_check').click().run()
    assert not at.exception and at.warning
    assert 'error' in at.session_state['_fixed_sd_returned_point_check']


def test_unsupported_or_inconsistent_inputs_clear_cached_check_and_multiple_contexts_work():
    at=AppTest.from_function(view,args=('ja',)).run(timeout=60)
    at.button(key='test_returned_point_check').click().run()
    calls=at.session_state['evaluations']
    fit=at.session_state['fit']
    fit['prep']['data']['score_k']=3-fit['prep']['data']['score_k']
    at.run()
    assert not at.exception and not at.button and at.session_state['evaluations']==calls
    assert '_fixed_sd_returned_point_check' not in at.session_state
    assert any(app._load_locale('ja')['estimation_subsections']['returned_point']['unavailable'] in c.value for c in at.caption)
    fit['prep']['data']['score_k']=fit['prep']['data']['Score'].to_numpy()
    at.session_state['two_contexts']=True
    at.run()
    assert not at.exception and len(at.button)==2
    assert not at.success


def test_input_identity_includes_current_parameters_and_settings_without_truncated_arrays():
    at=AppTest.from_function(view,args=('en',)).run(timeout=60)
    fit=at.session_state['fit'];original=app.prepare_fixed_sd_pcm_check(fit)
    assets=app._analysis_input_assets(fit)
    snapshot=json.loads(assets['analysis_snapshot.json'])
    np.testing.assert_array_equal(snapshot['returned_free_coordinates'],fit['opt'].x)
    fit['opt'].message+=' changed'
    assert app.prepare_fixed_sd_pcm_check(fit)['input_id']!=original['input_id']
    fit['config']['population_prior_sd']=1.6
    assert app.prepare_fixed_sd_pcm_check(fit)['input_id']!=original['input_id']
    with pytest.raises(ValueError,match='Reported NLL'):
        app.evaluate_fixed_sd_pcm_check(app.prepare_fixed_sd_pcm_check(fit))
