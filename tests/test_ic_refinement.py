"""Optional precision refits retain originals, failed attempts and inference holds."""
from copy import deepcopy
from io import BytesIO
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import engine, ic_sensitivity as check
from tests.test_ic_sensitivity import retained_fits
from tests.test_model_choice_exports import view
from tests.test_help_contextual_app_integration import _locale_value
from validation.mml_grouped_likelihood import GroupedMmlLikelihood


def unchanged(start, *args, **kwargs):
    return SimpleNamespace(x=start.copy(), success=True, message='controlled stop', nit=1, nfev=1)


def test_checked_coordinates_match_independent_formula_and_preserve_originals(monkeypatch):
    fits=retained_fits(); inputs=check.capture_inputs(fits); original=deepcopy(inputs)
    monkeypatch.setattr(engine,'mfrm_direct_mml',unchanged)
    result=check.refine(inputs)
    assert result['input_id']==check.input_id(original)
    assert check.input_id(inputs)==check.input_id(original)
    assert result['pairs'].Included.all()
    for row in result['candidates'].itertuples():
        assert row.CalculationStatus=='completed' and row.Feasible and not row.Boundary
        assert row.CalculationReason==''
        assert row.FitPoints==121 and row.CheckPoints==243
        assert not row.InferenceQualified and not row.ContinuousReferenceChecked and not row.OriginalResultReplaced
        f=fits[row.Model];point=np.asarray(result['points'][row.Model])
        reference=GroupedMmlLikelihood(f['prep']['data'],f['config'])
        rule=engine.gauss_hermite_normal(row.CheckPoints,sd=1.2)
        value=reference.evaluate(np.r_[point,np.log(1.2)],rule['nodes']/1.2,rule['weights'])
        assert row.CheckLogLik==pytest.approx(-value['nll'],abs=1e-10)
        assert row.CheckGradientSupnorm==pytest.approx(max(abs(value['gradient'][:-1])),abs=1e-10)
        assert row.AIC+row.CheckICChange==pytest.approx(2*value['nll']+2*row.KParams,abs=1e-10)
        assert row.BIC-row.AIC==pytest.approx((np.log(6)-2)*row.KParams,abs=1e-10)
        assert json.loads(row.OptimizerOptions)['ftol']==1e-12
    for row in result['pairs'].itertuples():
        a=result['candidates'].set_index('Model').loc[row.First]
        b=result['candidates'].set_index('Model').loc[row.Second]
        assert row.CheckedDifference==pytest.approx(b[row.Criterion]+b.CheckICChange-a[row.Criterion]-a.CheckICChange)


@pytest.mark.parametrize('problem',['changed','likelihood','different_sd'])
def test_all_originals_checked_before_any_optimization(monkeypatch,problem):
    fits=retained_fits()
    if problem=='likelihood': fits['GPCM']['summary']['LogLik']-=.1
    if problem=='different_sd': fits['GPCM']['config']['population_prior_sd']=1.3
    inputs=check.capture_inputs(fits)
    if problem=='changed': inputs['GPCM']['point'][0]+=.1
    monkeypatch.setattr(engine,'mfrm_direct_mml',lambda *a,**k:pytest.fail('Optimization must not start'))
    with pytest.raises(ValueError): check.refine(inputs)


@pytest.mark.parametrize('failure',['unfinished','string_success','nonfinite','wrong_shape','derived_bound','worse','exception','bad_gradient'])
def test_bad_candidate_retained_but_excluded_from_pairs(monkeypatch,failure):
    inputs=check.capture_inputs(retained_fits());calls=[]
    def optimizer(start,idx,cfg,sizes,rule,**kwargs):
        calls.append(cfg['model']);opt=unchanged(start)
        if cfg['model']!='GPCM': return opt
        if failure=='unfinished': opt.success=False
        if failure=='string_success': opt.success='true'
        if failure=='nonfinite': opt.x[0]=np.nan
        if failure=='wrong_shape': opt.x=opt.x[:-1]
        if failure=='derived_bound':
            assert sizes['log_slopes']==2
            opt.x[-2:]=2.0  # free coordinates legal, derived coordinate -4 outside [-3,3]
        if failure=='worse': opt.x[0]+=50
        if failure=='exception': raise RuntimeError('controlled failure')
        return opt
    monkeypatch.setattr(engine,'mfrm_direct_mml',optimizer)
    if failure=='bad_gradient':
        real=engine.mfrm_loglik_mml_value_grad
        def evaluate(point,idx,cfg,sizes,rule):
            nll,grad=real(point,idx,cfg,sizes,rule)
            if cfg['model']=='GPCM' and cfg['quad_points']==121: grad=grad[:-1]
            return nll,grad
        monkeypatch.setattr(engine,'mfrm_loglik_mml_value_grad',evaluate)
    result=check.refine(inputs);row=result['candidates'].set_index('Model').loc['GPCM']
    assert calls==list(check.MODELS)
    assert row.CalculationStatus in ('failed','unfinished')
    bad=result['pairs'].query("First=='GPCM' or Second=='GPCM'")
    assert not bad.Included.any() and bad.Difference.isna().all() and bad.CheckedDifference.isna().all()
    assert result['pairs'].query("Second=='PCM'").Included.all()
    assert not result['inference_qualified'] and not result['original_result_replaced']
    expected_reason={'unfinished':'termination_unconfirmed','string_success':'termination_unconfirmed',
        'nonfinite':'invalid_candidate','wrong_shape':'invalid_candidate','derived_bound':'infeasible_candidate',
        'worse':'objective_worsened','exception':'calculation_failed','bad_gradient':'nonfinite_evaluation'}
    assert row.CalculationReason==expected_reason[failure]
    if failure=='derived_bound': assert row.Failure=='infeasible_candidate' and not row.Feasible
    if failure=='exception': assert row.Failure=='controlled failure'


def test_boundary_is_flagged_without_treating_raw_gradient_as_certificate(monkeypatch):
    inputs=check.capture_inputs(retained_fits())
    def optimizer(start,idx,cfg,*args,**kwargs):
        opt=unchanged(start)
        if cfg['model']=='GPCM':opt.x[-2:]=[1.5,1.5]
        return opt
    monkeypatch.setattr(engine,'mfrm_direct_mml',optimizer)
    row=check.refine(inputs)['candidates'].set_index('Model').loc['GPCM']
    assert row.Feasible and row.Boundary and not row.InferenceQualified


@pytest.mark.parametrize('lang',['ja','en'])
def test_refit_requires_click_and_preserves_main_export_across_navigation(monkeypatch,lang):
    fits=retained_fits();calls=[]
    monkeypatch.setattr(app,'_refit_candidate_model',lambda res,model,**kw:fits[model])
    bundle=app.compute_model_choice_comparison(fits['RSM'],step_facet='Criterion',slope_facet='Criterion')
    def optimizer(*args,**kwargs): calls.append(1);return unchanged(*args,**kwargs)
    monkeypatch.setattr(engine,'mfrm_direct_mml',optimizer)
    exports={};download=app.st.download_button
    def capture(label,data,**kwargs):
        exports[kwargs['file_name']]=(data,kwargs);return download(label,data,**kwargs)
    monkeypatch.setattr(app.st,'download_button',capture)
    at=AppTest.from_function(view);at.session_state['result']=fits['RSM'];at.session_state['lang']=lang
    source=app._model_choice_cache_key(fits['RSM']);cache=app._model_comparison_cache_key(fits['RSM'],'Criterion','Criterion')
    at.session_state[f'model_choice_selection::{source}']={'step_facet':'Criterion','slopes':{'Criterion':'Criterion'}}
    state=f'model_choice_comparison::{cache}';at.session_state[state]=bundle
    at.run();assert not at.exception and not calls
    original_export=exports['mfrm_model_comparison.csv'][0];original_table=at.table[0].value.copy()
    label=lambda key:_locale_value(lang,'report_tables.model_choice_numerics.'+key)
    next(b for b in at.button if b.label==label('button')).click().run()
    assert not at.exception and not calls and 'mfrm_ic_refinement.csv' not in exports
    next(b for b in at.button if b.label==label('refinement.button')).click().run()
    assert not at.exception and len(calls)==3
    payload,options=exports['mfrm_ic_refinement.csv'];table=pd.read_csv(BytesIO(payload))
    assert len(table)==9 and options['on_click']=='ignore'
    assert table.Refitted.all() and not table.InferenceQualified.any() and not table.OriginalResultReplaced.any()
    assert table.SourceAnalysisID.eq(source.split('|',1)[-1]).all()
    assert table.ComparisonStepFacet.eq('Criterion').all() and table.ComparisonSlopeFacet.eq('Criterion').all()
    assert table.RefinementInputID.eq(check.input_id(bundle['numerical_inputs'])).all()
    for row in table.query("_section=='refined_candidates'").itertuples():
        assert json.loads(row.CandidateCoordinates)==bundle['numerical_inputs'][row.Model]['point'].tolist()
    assert exports['mfrm_model_comparison.csv'][0]==original_export
    pd.testing.assert_frame_equal(original_table,at.table[0].value)
    at.session_state['show_comparison']=False;at.run();at.session_state['show_comparison']=True;at.run()
    at.session_state['lang']='en' if lang=='ja' else 'ja';at.run()
    assert not at.exception and len(calls)==3 and len(at.table)==4
    heading=_locale_value('en' if lang=='ja' else 'ja','report_tables.model_choice_numerics.refinement.heading')
    assert any(element.value=='**'+heading+'**' for element in at.markdown)
    exports.clear();at.session_state[state]['numerical_refinement']['inference_qualified']=True;at.run()
    assert not at.exception and len(calls)==3 and 'mfrm_ic_refinement.csv' not in exports


def _refinement_view():
    import streamlit as st
    import streamlit_app as app
    bundle=st.session_state['refinement_bundle']
    app._render_model_comparison_refinement(bundle,bundle['numerical_inputs'],'fixture')


@pytest.mark.parametrize('lang',['ja','en'])
@pytest.mark.parametrize('scenario',['mixed','all_failed','legacy','boundary'])
def test_readable_outcomes_match_csv_without_raw_errors_or_new_calculation(monkeypatch,lang,scenario):
    inputs=check.capture_inputs(retained_fits())
    monkeypatch.setattr(engine,'mfrm_direct_mml',unchanged)
    record=check.refine(inputs);table=record['candidates']
    reasons={}
    if scenario=='mixed':
        for model,reason in [('PCM','termination_unconfirmed'),('GPCM','objective_worsened')]:
            table.loc[table.Model.eq(model),['CalculationStatus','CalculationReason']]=['unfinished',reason]
            reasons[model]=reason
    if scenario=='all_failed':
        table['CalculationStatus']='failed'
        table['CalculationReason']='calculation_failed'
        table['Failure']='INTERNAL_TEST_TRACE'
        table[['AIC','BIC','CheckICChange','CheckGradientSupnorm']]=np.nan
        reasons=dict.fromkeys(check.MODELS,'calculation_failed')
    if scenario=='legacy':
        table.loc[table.Model.eq('PCM'),'CalculationStatus']='unfinished'
        table.drop(columns='CalculationReason',inplace=True)
        reasons['PCM']='unknown'
    if scenario=='boundary':
        table.loc[table.Model.eq('GPCM'),['Boundary','Feasible']]=[True,True]
        table.loc[table.Model.eq('PCM'),['Boundary','Feasible','CalculationStatus','CalculationReason']]=[True,False,'failed','infeasible_candidate']
        reasons['PCM']='infeasible_candidate'
    # Table rendering uses saved outcomes; numerical and pair rules are tested above.
    monkeypatch.setattr(check,'refine',lambda *a,**k:pytest.fail('Reading must not refit'))
    exports={};download=app.st.download_button
    def capture(label,data,**kwargs):
        exports[kwargs['file_name']]=(data,kwargs);return download(label,data,**kwargs)
    monkeypatch.setattr(app.st,'download_button',capture)
    bundle=dict(numerical_inputs=inputs,numerical_refinement=record,source_result_key='saved|source',
        comparison_spec=dict(step_facet='Criterion',slope_facet='Criterion'))
    at=AppTest.from_function(_refinement_view);at.session_state['lang']=lang
    at.session_state['refinement_bundle']=bundle;at.run()
    assert not at.exception and not at.button and not at.success
    tr=lambda key:_locale_value(lang,'report_tables.model_choice_numerics.refinement.'+key)
    assert any(element.value=='**'+tr('heading')+'**' for element in at.markdown)
    captions=[element.value for element in at.caption]
    payload,options=exports['mfrm_ic_refinement.csv'];csv=pd.read_csv(BytesIO(payload)).query("_section=='refined_candidates'").set_index('Model')
    assert options['on_click']=='ignore' and not csv.InferenceQualified.any()
    for model in check.MODELS:
        if model in reasons:
            note=tr('excluded').format(reason=tr('reasons.'+reasons[model]))
            assert f'{model}: {note}' in captions
            assert csv.loc[model,'CalculationNote']==note
        else: assert pd.isna(csv.loc[model,'CalculationNote'])
        status=table.set_index('Model').loc[model,'CalculationStatus']
        assert csv.loc[model,'CalculationStatusLabel']==tr('status_'+status)
    if scenario=='boundary':
        assert tr('boundary').format(models='GPCM') in captions
        assert tr('boundary').format(models='PCM') not in captions
    assert 'INTERNAL_TEST_TRACE' not in '\n'.join(captions)
    assert not {'CalculationReason','Failure','Feasible'}.intersection(at.table[0].value.columns)
