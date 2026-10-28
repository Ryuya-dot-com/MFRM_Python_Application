"""Transport the sealed joint kernel through the real app, retaining G0/G1 holds."""
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from mfrm_app import mml_joint
from tests.test_mfrm_recovery_design import design
from tests.test_mml_grouped_likelihood import problem
from validation.mfrm_recovery_design import generate_population_sample
from validation.mml_grouped_fit import fit_grouped_joint, normal_rule, OPTIONS, LIMITS
from validation.mml_population_study import prepare_retained_problem, json_bytes


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('assignment',['regular','sparse','unbalanced'])
def test_app_kernel_is_exactly_the_sealed_finite_rule_with_full_information(model,assignment):
    data,config,sizes,idx,reference=problem(model,assignment,duplicate=True)
    actual=mml_joint.GroupedMmlLikelihood(data,config)
    assert mml_joint.OPTIONS==OPTIONS and mml_joint.LIMITS==LIMITS
    for name in ['counts','design','slope_design','cell_totals','k']:
        np.testing.assert_array_equal(getattr(actual,name),getattr(reference,name))
    assert actual.person_labels==reference.person_labels
    assert actual.structural_recession_audit()==reference.structural_recession_audit()
    point=np.r_[np.linspace(-.4,.5,sum(sizes.values())),np.log(1.2)]
    for q in [31,1281]:
        z,w=normal_rule(q)
        expected=reference.evaluate(point,z,w,information=True)
        got=actual.evaluate(point,z,w,information=True)
        for key in expected:np.testing.assert_array_equal(got[key],expected[key])


def sample(model='RSM',assignment='regular',persons=24):
    return generate_population_sample(replace(design(model,assignment,persons=persons),rating_min=0),91537)


def arguments(value,model='RSM',**kw):
    return dict(data=value['data'],person_col='Person',facet_cols=['Rater','Task','Criterion'],
        score_col='Score',rating_min=0,rating_max=3,keep_original=True,model=model,
        method='MML',noncenter_facet='Criterion',step_facet=None if model=='RSM' else 'Criterion',
        mml_engine='joint',estimate_population_sd=True,quad_points=31,**kw)


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('assignment',['regular','sparse','unbalanced'])
def test_real_app_joint_fit_matches_registered_estimator_and_retains_inference_hold(model,assignment):
    import hashlib
    value=sample(model,assignment)
    csv=value['data'].to_csv(index=False)
    retained=dict(value,csv_utf8=csv,csv_sha256=hashlib.sha256(csv.encode()).hexdigest())
    data,p,invalid=prepare_retained_problem(retained,model,31)
    assert invalid is None
    reference=fit_grouped_joint(data,p,q=31)
    result=app.mfrm_estimate(**arguments(value,model,compute_plausible_values=True,n_plausible_values=3))
    actual=result['joint_mml']
    for key in ['point','covariance','evaluated','checks','numerical_pass','status','recession']:
        left,right=actual[key],reference[key]
        if key=='evaluated':
            left=dict(left);right={k:v for k,v in right.items() if k!='posterior'}
        assert json_bytes(left)==json_bytes(right),key
    for actual_run,expected_run in zip(actual['runs'],reference['runs'],strict=True):
        assert json_bytes({k:v for k,v in actual_run.items() if k!='evaluations'})==json_bytes(expected_run)
    assert result['config']['mml_engine']=='joint'
    assert 'mml_engine_requested_before_sd_override' not in result['config']
    np.testing.assert_array_equal(result['opt'].x,reference['point'][:-1])
    sigma=np.exp(reference['point'][-1])
    assert result['config']['estimated_population_sd']==sigma
    assert result['config']['population_prior_sd']==sigma
    assert result['summary'].LogLik.iloc[0]==pytest.approx(-reference['evaluated']['nll'],abs=1e-10)
    np.testing.assert_allclose(result['facets']['person'].Estimate,reference['evaluated']['eap'],atol=1e-12,rtol=0)
    np.testing.assert_allclose(result['posterior']['scores'].PosteriorSD,reference['evaluated']['sd'],atol=1e-12,rtol=0)
    assert not result['summary'].InferenceReady.iloc[0]
    assert not actual['scientific_inference_ready'] and not actual['quadrature_accuracy_qualified']
    assert result['config']['population_sd_se'] is None
    assert result['config']['population_sd_ci'] is None
    assert result['convergence'].GradientScope.iloc[0]=='structural coordinates and log(population SD)'
    assert result['convergence'].JointStationarityPass.iloc[0]==reference['numerical_pass']
    assert result['summary'].KParams.iloc[0]==len(reference['point'])
    assert result['config']['parameterization']['optimizer_dimension']==len(reference['point'])
    assert app._normalize_mml_engine_contract('joint')=='joint'


@pytest.mark.parametrize('case',['fixed_sd','jmle','weights','anchors','dummy','penalty','noncenter_person','population_regression','q1','single_category','structural_ray','rank_deficient'])
def test_joint_scope_and_prefit_guards_run_before_optimizer(monkeypatch,case):
    value=sample(persons=24);kwargs=arguments(value)
    if case=='fixed_sd':kwargs['estimate_population_sd']=False
    elif case=='jmle':kwargs['method']='JMLE'
    elif case=='weights':
        value['data']['Weight']=1.5;kwargs['weight_col']='Weight'
    elif case=='anchors':kwargs['anchor_df']=pd.DataFrame([dict(Facet='Rater',Level='R0',Value=0.)])
    elif case=='dummy':kwargs['dummy_facets']=['Task']
    elif case=='penalty':kwargs['facet_regularization']=dict(enabled=True,facets=['Rater'],sd=1.)
    elif case=='noncenter_person':kwargs['noncenter_facet']='Person'
    elif case=='population_regression':
        people=sorted(value['data'].Person.unique())
        kwargs.update(person_data=pd.DataFrame(dict(Person=people,Age=np.arange(len(people)))),population_formula='~ Age')
    elif case=='q1':kwargs['quad_points']=1
    elif case=='single_category':value['data']['Score']=0
    elif case=='structural_ray':
        value['data']['Score']=np.where(value['data'].Rater=='R0',0,3)
    elif case=='rank_deficient':value['data']['Task']=value['data']['Rater']
    def no_fit(*a,**kw):pytest.fail('An unsupported/unidentified joint fit reached the optimizer')
    monkeypatch.setattr(mml_joint,'minimize',no_fit)
    with pytest.raises(ValueError):app.mfrm_estimate(**kwargs)


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
def test_joint_diagnostics_and_both_download_paths_preserve_technical_scope(model):
    result=app.mfrm_estimate(**arguments(sample(model),model))
    diag=app.mfrm_diagnostics(result,compute_pca=False)
    uncertainty=diag['uncertainty']
    assert uncertainty['covariance']['status']=='unsupported_covariance_scope'
    assert uncertainty['structural_se'][['StructuralSE','StructuralCI_Lower','StructuralCI_Upper']].isna().all().all()
    assert not uncertainty['covariance_audit'].PublicConclusionAllowed.any()
    # Keep the technical uncertainty basis while carrying the calibration hold.
    persons=diag['measures'].Facet=='Person'
    assert diag['measures'].SE_Status.eq('calibration_inference_withheld').all()
    assert diag['measures'].CI_Status.eq('calibration_inference_withheld').all()
    assert diag['measures'].loc[persons,'SE_Method'].str.contains('posterior SD').all()
    assert diag['measures'].loc[~persons,'SE_Method'].str.contains('conditional information').all()
    assert result['steps'].SE.isna().all()
    expected=app.build_joint_mml_audit_frames(result)
    quick=app.build_result_bundle_frames(result,diag)
    bundles=[quick]
    # The shared full-download collector also assembles many unrelated report
    # panels. Exercise its joint transport once; check all three models in
    # the core bundle and in separately executed generated runners below.
    if model=='RSM':
        full,_=app.collect_download_frames(result,diag,{},app.compute_scorefile(result),
            app.compute_residual_file(result),public_export_mode=False)
        bundles.append(full)
    for frames in bundles:
        for name,frame in expected.items():
            pd.testing.assert_frame_equal(frames[name],frame)
            assert not frame.ScientificInferenceReady.any()
            assert not frame.PublicConclusionAllowed.any()
        qualification=frames['output_qualification'].set_index('OutputID')
        assert qualification.loc['mml.joint_numerical_audit','QualificationStatus']=='TECHNICAL_ONLY'
        assert qualification.loc['mml.joint_numerical_audit','RawTechnicalExportAllowed']
        for key in ['mml.structural_covariance','mml.population_sd_uncertainty']:
            assert qualification.loc[key,'QualificationStatus']=='WITHHELD'
            assert not qualification.loc[key,'PublicConclusionAllowed']


def test_stale_joint_evidence_is_not_exported():
    from copy import deepcopy
    result=app.mfrm_estimate(**arguments(sample()))
    for kind in ['point','sigma','quadrature','input']:
        changed=deepcopy(result)
        if kind=='point':changed['opt'].x[0]+=.01
        elif kind=='sigma':changed['config']['population_prior_sd']+=.01
        elif kind=='quadrature':changed['config']['quad_points']+=2
        else:changed['prep']['data'].loc[0,'Score']=(changed['prep']['data'].loc[0,'Score']+1)%4
        with pytest.raises(ValueError,match='evidence differs'):
            app.build_joint_mml_audit_frames(changed)


@pytest.mark.parametrize('model,engine',[('RSM','joint'),('PCM','joint'),('GPCM','joint'),('RSM','EM')])
def test_generated_runner_reproduces_free_sd_from_original_start_in_another_process(tmp_path,model,engine):
    import subprocess
    import sys
    from pathlib import Path
    value=sample(model)
    kwargs=arguments(value,model,population_prior_sd=1.2)
    kwargs['mml_engine']=engine
    result=app.mfrm_estimate(**kwargs)
    diag=app.mfrm_diagnostics(result,compute_pca=False)
    source=app._generate_app_engine_runner_script(result,diag)
    assert 'ESTIMATE_POPULATION_SD = True' in source
    assert 'POPULATION_PRIOR_SD = 1.2' in source
    runner=tmp_path/'runner.py';runner.write_text(source)
    datafile=tmp_path/'input.csv';value['data'].to_csv(datafile,index=False)
    output=tmp_path/'output'
    with (tmp_path/'runner.log').open('w') as log:
        run=subprocess.run([sys.executable,'-B',str(runner),'--input',str(datafile),
            '--app-path',str(Path(app.__file__).resolve()),'--output-dir',str(output)],
            cwd=tmp_path,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    assert run.returncode==0,(tmp_path/'runner.log').read_text()[-6000:]
    summary=pd.read_csv(output/'mfrm_summary.csv')
    for key in ['LogLik','KParams']:
        assert summary[key].iloc[0]==pytest.approx(result['summary'][key].iloc[0],abs=1e-10,rel=0)
    assert not summary.InferenceReady.any()
    assert (output/'mfrm_person_measures.csv').read_text()==result['facets']['person'].to_csv(index=False)
    assert (output/'mfrm_facet_measures.csv').read_text()==result['facets']['others'].to_csv(index=False)
    for name,frame in app.build_joint_mml_audit_frames(result).items():
        assert (output/f'mfrm_{name}.csv').read_text()==frame.to_csv(index=False)
    qualification=pd.read_csv(output/'mfrm_output_qualification.csv').set_index('OutputID')
    assert qualification.loc['mml.population_sd_uncertainty','QualificationStatus']=='WITHHELD'
    assert not qualification.PublicConclusionAllowed.any()
