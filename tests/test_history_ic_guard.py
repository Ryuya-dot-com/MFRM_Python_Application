"""History differences require matching recorded likelihood conditions and termination."""
from copy import deepcopy
from io import BytesIO
import json

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import engine
from tests.test_ic_sensitivity import retained_fits
from tests.test_help_contextual_app_integration import _locale_value


def outputs():
    fit=retained_fits()['RSM']
    fit['config'].update(population_model={'enabled':False},facet_regularization_enabled=False,population_sd_bounds=(.05,10.))
    fit['summary']['Deviance']=-2*fit['summary'].LogLik
    return [{'result':deepcopy(fit),'model_type':'RSM','_est_method':'MML','facet_cols':['Rater','Task','Criterion']} for _ in range(2)]


def history(items):
    return [{'timestamp':str(i),'output_snapshot':item} for i,item in enumerate(items)]


@pytest.mark.parametrize('output',[None,{}, {'result':None}])
def test_missing_legacy_result_has_no_comparability_key(output):
    assert app.extract_run_likelihood_row(output)['ComparabilityKey']==''


def test_matching_completed_history_differences_preserve_values_and_holds():
    items=outputs();items[1]['result']['summary']['LogLik']+=2
    items[1]['result']['summary']['Deviance']-=4;items[1]['result']['summary']['AIC']-=4
    original=deepcopy(items)
    table=app.build_run_history_likelihood_table(history(items))
    assert table.ComparableIC.all() and table.AlgorithmTerminated.all()
    assert table.ICComparisonStatus.eq('matched_completed').all()
    for metric in ['AIC','BIC','Deviance']:np.testing.assert_allclose(table['Delta'+metric],[4,0])
    qualified=app.qualify_likelihood_information_table(table,app.build_output_qualification_table({'config':{'method':'MML'}}))
    assert not qualified.PublicConclusionAllowed.any()
    assert app.build_likelihood_descriptives(table).set_index('Metric').loc['AIC','Mean']==pytest.approx(table.AIC.mean())
    for a,b in zip(items,original):pd.testing.assert_frame_equal(a['result']['summary'],b['result']['summary'])


@pytest.mark.parametrize('change',['data','weights','fixed_sd','free_sd','points','population','anchors','signs','model','method','input_fingerprint'])
def test_different_conditions_never_receive_differences(change):
    items=outputs();fit=items[1]['result'];cfg=fit['config']
    if change=='data':fit['prep']['data'].loc[0,'score_k']+=1
    if change=='weights':fit['prep']['data'].loc[0,'Weight']=2.
    if change=='fixed_sd':cfg['population_prior_sd']=1.7
    if change=='free_sd':cfg['estimate_population_sd']=True
    if change=='points':cfg['quad_points']=61
    if change=='population':cfg['population_model']={'enabled':True,'X':np.arange(6).reshape(-1,1)}
    if change=='anchors':cfg['facet_specs']['Rater']['anchors'][0]=.5
    if change=='signs':cfg['facet_signs']['Rater']=1
    if change=='model':cfg['model']='PCM';fit['summary']['Model']='PCM';cfg['step_facet']='Criterion'
    if change=='method':cfg['method']='JMLE';cfg['theta_spec']={'n_params':6};fit['summary']['Method']='JMLE'
    if change=='input_fingerprint':cfg['input_data_fingerprint']='another-input'
    table=app.build_run_history_likelihood_table(history(items))
    assert not table.ComparableIC.any()
    assert table.ICComparisonStatus.eq('different_conditions').all()
    assert table[['DeltaAIC','DeltaBIC','DeltaDeviance']].isna().all().all()
    assert table.AIC.notna().all()
    assert app.build_likelihood_descriptives(table).empty


@pytest.mark.parametrize('change',['data','points','sd','mode','population','facets','constraints','regularization'])
def test_identically_missing_or_unsupported_conditions_are_not_a_match(change):
    items=outputs()
    for item in items:
        fit=item['result'];cfg=fit['config']
        if change=='data':fit['prep']={}
        if change=='points':cfg['quad_points']=None
        if change=='sd':cfg.pop('population_prior_sd')
        if change=='mode':cfg.pop('estimate_population_sd')
        if change=='population':cfg.pop('population_model')
        if change=='facets':cfg['facet_names']=None
        if change=='constraints':cfg.pop('facet_specs')
        if change=='regularization':cfg['facet_regularization_enabled']=True
    table=app.build_run_history_likelihood_table(history(items))
    status='unsupported_scope' if change=='regularization' else 'unknown_conditions'
    assert not table.ComparableIC.any() and table.ICComparisonStatus.eq(status).all()
    assert table[['DeltaAIC','DeltaBIC','DeltaDeviance']].isna().all().all()


@pytest.mark.parametrize('success',[False,None,'true'])
def test_unfinished_lowest_value_cannot_be_the_history_reference(success):
    items=outputs();fit=items[1]['result'];fit['opt'].success=success
    fit['summary']['Converged']=True;fit['summary']['AIC']-=20
    table=app.build_run_history_likelihood_table(history(items))
    assert table.ComparableIC.all()  # Conditions match; completion is a separate condition.
    assert table.AlgorithmTerminated.tolist()==[True,False]
    assert table.ICComparisonStatus.eq('termination_unconfirmed').all()
    assert table[['DeltaAIC','DeltaBIC','DeltaDeviance']].isna().all().all()
    assert table.AIC.iloc[1]<table.AIC.iloc[0]


@pytest.mark.parametrize('value',[np.nan,np.inf,-np.inf])
def test_nonfinite_criterion_or_single_run_has_no_delta(value):
    items=outputs();items[1]['result']['summary']['AIC']=value
    table=app.build_run_history_likelihood_table(history(items))
    assert table.ICComparisonStatus.eq('nonfinite').all() and table.DeltaAIC.isna().all()
    table=app.build_run_history_likelihood_table(history(items[:1]))
    assert table.ICComparisonStatus.eq('insufficient_runs').all() and table.DeltaAIC.isna().all()


def test_free_sd_estimates_and_optimizer_settings_are_not_fixed_model_conditions():
    items=outputs()
    for item in items:item['result']['config']['estimate_population_sd']=True
    items[1]['result']['config'].update(population_prior_sd=1.8,estimated_population_sd=1.8,reltol=1e-10,maxit=1000,mml_engine='em')
    table=app.build_run_history_likelihood_table(history(items))
    assert table.ComparableIC.all() and table.ICComparisonStatus.eq('matched_completed').all()
    items[1]['result']['config']['population_sd_bounds']=(.1,5.)
    assert not app.build_run_history_likelihood_table(history(items)).ComparableIC.any()


def test_compact_mml_bic_requires_explicit_recorded_person_basis_and_preserves_old_value():
    entry=dict(method='MML',model='RSM',n_obs=24,k_params=5,loglik=-10.,aic=30.,bic=99.,converged=True)
    old=app.build_run_history_likelihood_table([entry]).iloc[0]
    assert pd.isna(old.BIC) and old.SavedBIC==99
    modern={**entry,'bic_sample_size':6,'bic_sample_size_basis':'independent_person_marginals'}
    row=app.build_run_history_likelihood_table([modern]).iloc[0]
    assert row.BIC==pytest.approx(20+5*np.log(6)) and row.BICSampleSize==6 and row.SavedBIC==99
    assert not row.ComparableIC and pd.isna(row.DeltaBIC)


@pytest.mark.parametrize('lang',['ja','en'])
def test_history_and_pair_view_explain_omitted_differences_and_export_consistently(monkeypatch,lang):
    items=outputs();items[1]['result']['opt'].success=False
    # Stored raw BIC is stale: compact and detailed views must both rebuild it.
    for item in items:item['result']['summary']['BIC']=999.
    at=AppTest.from_string('''
import streamlit_app as app
if not app.st.session_state.get('seeded'):
 app.st.session_state[app._RUN_HISTORY_STORE_SNAPSHOTS_KEY]=True
 for output in app.st.session_state['outputs']:app.record_run_in_history(output=output,elapsed_sec=1.)
 for entry in app.st.session_state[app._RUN_HISTORY_KEY]:entry['bic']=999.
 app.st.session_state['seeded']=True
app.render_run_history_panel()
app.render_comparison_panel(*app.st.session_state['outputs'])
''')
    at.session_state['lang']=lang;at.session_state['outputs']=items
    exports={};download=app.st.download_button
    def capture(label,data,**kwargs):exports[kwargs['file_name']]=(data,kwargs);return download(label,data,**kwargs)
    monkeypatch.setattr(app.st,'download_button',capture)
    monkeypatch.setattr(engine,'mfrm_estimate',lambda *a,**k:pytest.fail('No fitting while reading'))
    at.run();assert not at.exception
    note=_locale_value(lang,'history.criteria_status.termination_unconfirmed')
    assert [c.value for c in at.caption].count(note)==2
    payload,options=exports['mfrm_run_history_likelihood_criteria.csv'];table=pd.read_csv(BytesIO(payload))
    assert options['on_click']=='ignore' and table.ICComparisonNote.eq(note).all()
    assert table[['DeltaAIC','DeltaBIC','DeltaDeviance']].isna().all().all()
    assert not table.PublicConclusionAllowed.any() and table.ComparableIC.all()
    assert table.SavedBIC.eq(999.).all()
    assert at.dataframe[0].value.BIC.ne(999.).all()
    np.testing.assert_allclose(at.dataframe[0].value.BIC,table.BIC.iloc[::-1].round(3))
    for frame in at.dataframe:
        if 'DeltaAIC' in frame.value:assert frame.value.DeltaAIC.isna().all()
        assert 'ICComparisonStatus' not in frame.value
    entries=at.session_state[app._RUN_HISTORY_KEY]
    assert all(e['bic_sample_size']==6 for e in entries)
    assert all(e['bic']==999. for e in entries)  # Saved original retained; only the view is rebuilt.
    at.run();assert not at.exception
    assert exports['mfrm_run_history_likelihood_criteria.csv'][0]==payload
