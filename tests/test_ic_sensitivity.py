"""Checks use unchanged points, keep provenance, and cannot qualify inference."""
from copy import deepcopy
from io import BytesIO
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import ic_sensitivity as check
from tests.test_mml_grouped_likelihood import problem
from validation.mml_grouped_likelihood import GroupedMmlLikelihood
from tests.test_help_contextual_app_integration import _locale_value
from tests.test_model_choice_exports import view


def retained_fits(q=15):
    data=problem('RSM','regular')[0]
    prep=app.prepare_mfrm_data(data,'Person',['Rater','Task','Criterion'],'Score',rating_min=0,rating_max=3)
    fits={}
    for model in check.MODELS:
        _,cfg,sizes,_,_=problem(model,'regular')
        cfg.update(quad_points=q,population_prior_sd=1.2,estimate_population_sd=False)
        point=np.linspace(-.2,.3,sum(sizes.values()))
        rule=app.gauss_hermite_normal(q,sd=1.2)
        independent=GroupedMmlLikelihood(data,cfg)
        nll=independent.evaluate(np.r_[point,np.log(1.2)],rule['nodes']/1.2,rule['weights'])['nll']
        k=len(point);n=cfg['n_person']
        fits[model]=dict(config=cfg,prep=prep,opt=SimpleNamespace(x=point,success=True,message='fixture'),
            summary=pd.DataFrame([dict(Model=model,Method='MML',N=len(data),Persons=n,KParams=k,
                LogLik=-nll,AIC=2*nll+2*k,BIC=2*nll+k*np.log(n))]))
    return fits


def test_all_rules_match_independent_formula_and_pair_changes_without_mutation():
    fits=retained_fits();original=deepcopy(fits)
    inputs=check.capture_inputs(fits);identity=check.input_id(inputs)
    result=check.evaluate(inputs)
    assert len(result['evaluations'])==9 and len(result['pairs'])==6
    assert result['input_id']==identity and not result['inference_qualified'] and not result['refitted']
    for model in check.MODELS:
        f=fits[model];independent=GroupedMmlLikelihood(f['prep']['data'],f['config'])
        for row in result['evaluations'].query('Model == @model').itertuples():
            rule=app.gauss_hermite_normal(row.Points,sd=1.2)
            ref=independent.evaluate(np.r_[f['opt'].x,np.log(1.2)],rule['nodes']/1.2,rule['weights'])
            assert row.LogLik==pytest.approx(-ref['nll'],abs=1e-10)
            assert row.GradientSupnorm==pytest.approx(max(abs(ref['gradient'][:-1])),abs=1e-10)
            assert row.ICChange==pytest.approx(row.AIC-f['summary'].AIC.iloc[0],abs=1e-10)
            assert not row.InferenceQualified and not row.Refitted and not row.ContinuousReferenceChecked
        np.testing.assert_array_equal(f['opt'].x,original[model]['opt'].x)
        pd.testing.assert_frame_equal(f['summary'],original[model]['summary'])
    pairs=result['pairs'];np.testing.assert_allclose(pairs.FinerRuleDifference-pairs.OriginalDifference,pairs.DifferenceChange,atol=1e-12)
    fits['RSM']['opt'].x[0]+=1
    assert check.input_id(inputs)==identity  # retained snapshot is independent of later result edits


@pytest.mark.parametrize('mutation', ['free_sd','method','weights','large','bad_k','bad_levels','fractional_q'])
def test_unsupported_or_inconsistent_inputs_are_refused(mutation):
    fits=retained_fits();f=fits['RSM']
    if mutation=='free_sd':f['config']['estimate_population_sd']=True
    if mutation=='method':f['config']['method']='JMLE'
    if mutation=='weights':f['prep']['data']['Weight']=2.
    if mutation=='large':f['prep']['data']=pd.concat([f['prep']['data']]*100,ignore_index=True)
    if mutation=='bad_k':f['summary']['KParams']+=1
    if mutation=='bad_levels':f['prep']['levels']['Rater']=list(reversed(f['prep']['levels']['Rater']))
    if mutation=='fractional_q':f['config']['quad_points']=15.5
    with pytest.raises((ValueError,TypeError)):
        check.capture_inputs(fits)


def test_changed_points_and_mismatched_likelihoods_are_not_silently_used():
    inputs=check.capture_inputs(retained_fits())
    inputs['RSM']['point'][0]+=.2
    with pytest.raises(ValueError,match='changed_inputs'):check.evaluate(inputs)
    fits=retained_fits();fits['PCM']['config']['population_prior_sd']=1.3
    with pytest.raises(ValueError,match='different_likelihoods'):check.evaluate(check.capture_inputs(fits))
    fits=retained_fits();fits['RSM']['summary']['LogLik']-=.5
    with pytest.raises(ValueError,match='reported_likelihood_mismatch'):check.evaluate(check.capture_inputs(fits))


def test_unknown_termination_is_not_coerced_to_success():
    fits=retained_fits();fits['RSM']['opt'].success='true'
    assert check.capture_inputs(fits)['RSM']['terminated'] is False


@pytest.mark.parametrize('lang',['ja','en'])
def test_explicit_check_retains_results_across_navigation_language_and_csv(monkeypatch,lang):
    fits=retained_fits();calls=[]
    monkeypatch.setattr(app,'_refit_candidate_model',lambda res,model,**kw:fits[model])
    bundle=app.compute_model_choice_comparison(fits['RSM'],step_facet='Criterion',slope_facet='Criterion')
    real=check.evaluate
    def counted(inputs):
        calls.append(1);return real(inputs)
    monkeypatch.setattr(check,'evaluate',counted)
    monkeypatch.setattr(app,'_refit_candidate_model',lambda *a,**k:pytest.fail('No refitting during check'))
    exports={};download=app.st.download_button
    def capture(label,data,**kwargs):
        exports[kwargs['file_name']]=(data,kwargs);return download(label,data,**kwargs)
    monkeypatch.setattr(app.st,'download_button',capture)
    at=AppTest.from_function(view)
    at.session_state['result']=fits['RSM'];at.session_state['lang']=lang
    source=app._model_choice_cache_key(fits['RSM']);cache=app._model_comparison_cache_key(fits['RSM'],'Criterion','Criterion')
    at.session_state[f'model_choice_selection::{source}']={'step_facet':'Criterion','slopes':{'Criterion':'Criterion'}}
    at.session_state[f'model_choice_comparison::{cache}']=bundle
    at.run();assert not calls and not at.exception
    before=at.table[0].value.copy()
    label=_locale_value(lang,'report_tables.model_choice_numerics.button')
    next(button for button in at.button if button.label==label).click().run()
    assert not at.exception and len(calls)==1
    payload,options=exports['mfrm_ic_numerical_check.csv'];table=pd.read_csv(BytesIO(payload))
    assert len(table)==15 and options['on_click']=='ignore'
    assert not table.InferenceQualified.any() and not table.Refitted.any()
    assert table.SourceAnalysisID.eq(source.split('|',1)[-1]).all()
    at.session_state['show_comparison']=False;at.run()
    at.session_state['show_comparison']=True;at.run()
    assert not at.exception and len(calls)==1
    pd.testing.assert_frame_equal(at.table[0].value,before)
    at.session_state['lang']='en' if lang=='ja' else 'ja';at.run()
    assert not at.exception and len(calls)==1
    assert len(at.table)==3
    exports.clear()
    at.session_state[f'model_choice_comparison::{cache}']['numerical_check']['inference_qualified']=True
    at.run()
    assert not at.exception and len(calls)==1
    assert 'mfrm_ic_numerical_check.csv' not in exports and len(at.table)==2
