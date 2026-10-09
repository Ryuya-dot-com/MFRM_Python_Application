"""Standalone fair-average uncertainty must carry the active G0 contract."""
from contextlib import nullcontext
from copy import deepcopy
from io import StringIO
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from tests.test_fair_se_cache_identity import fixture as identity_fixture, ui_cache_block


@pytest.fixture
def inputs(monkeypatch):
    result=dict(config=dict(model='GPCM',method='MML',estimate_population_sd=False),
        opt=SimpleNamespace(x=np.array([.1,.2])),prep=dict(rating_min=0.,rating_max=3.))
    tables={facet:pd.DataFrame(dict(Level=['C0','C1'] if facet=='Criterion' else ['P0','P1'],
        FairM=[.8,.8],FairZ=[.9,.9])) for facet in ['Criterion','Person']}
    covariance=dict(status='ok',rank=2,regularized=False,cov=np.array([[.04,.01],[.01,.09]]))
    monkeypatch.setattr(app,'_compute_fair_average_value_from_par',
        lambda res,point,facet,level,metric: .3+point[0]+2*point[1] if metric=='FairM' else 1+point[0]-point[1])
    return result,tables,covariance


@pytest.mark.parametrize('status,rank,qualified,reason',[
    ('ok',2,'TECHNICAL_ONLY','stat.stat_003.covariance_coverage_unvalidated'),
    ('regularized',1,'WITHHELD','stat.stat_003.covariance_rank_deficient'),
    ('ok',1,'WITHHELD','stat.stat_003.covariance_rank_deficient'),
])
def test_computed_values_and_active_qualification_are_both_retained(inputs,status,rank,qualified,reason):
    result,tables,covariance=inputs;before=deepcopy(tables)
    covariance.update(status=status,rank=rank,regularized=status=='regularized',
        audit=pd.DataFrame([dict(QualificationStatus='READY',PublicConclusionAllowed=True)]))
    output=app.add_gpcm_fair_average_delta_se(tables,result,covariance)
    for metric,gradient in [('FairM',np.array([1.,2.])),('FairZ',np.array([1.,-1.]))]:
        frame=output['Criterion'];se=np.sqrt(gradient@covariance['cov']@gradient)
        np.testing.assert_allclose(frame[metric+'SE'],se,atol=1e-10)
        assert set(frame[metric+'_SE_QualificationStatus'])=={qualified}
        assert set(frame[metric+'_SE_ReasonCode'])=={reason}
        assert frame[metric+'_SE_RawTechnicalExportAllowed'].all()
        assert not frame[metric+'_SE_PublicConclusionAllowed'].any()
        assert frame[metric+'_SE_Detail'].str.contains('coverage|confirmatory',case=False).all()
        persons=output['Person']
        assert persons[metric+'SE'].isna().all()
        assert set(persons[metric+'_SE_QualificationStatus'])=={'WITHHELD'}
        assert not persons[metric+'_SE_PublicConclusionAllowed'].any()
        assert not persons[metric+'_SE_RawTechnicalExportAllowed'].any()
    for facet in tables:pd.testing.assert_frame_equal(tables[facet],before[facet])


@pytest.mark.parametrize('case',['free_sd','penalty','missing_covariance','dimension','JMLE','PCM'])
def test_unsupported_scope_cannot_be_qualified_by_supplied_covariance(inputs,case):
    result,tables,covariance=inputs
    if case=='free_sd':result['config']['estimate_population_sd']=True
    elif case=='penalty':result['config']['facet_regularization_enabled']=True
    elif case=='missing_covariance':covariance.update(cov=None,status='fallback')
    elif case=='dimension':covariance['cov']=np.eye(3)
    elif case=='JMLE':result['config']['method']='JMLE'
    else:result['config']['model']='PCM'
    output=app.add_gpcm_fair_average_delta_se(tables,result,covariance)
    expected='NOT_APPLICABLE' if case in {'JMLE','PCM'} else 'WITHHELD'
    for frame in output.values():
        for metric in ['FairM','FairZ']:
            assert frame[[metric+'SE',metric+'_CI_Lower',metric+'_CI_Upper']].isna().all().all()
            assert set(frame[metric+'_SE_QualificationStatus'])=={expected}
            assert not frame[metric+'_SE_PublicConclusionAllowed'].any()
            assert not frame[metric+'_SE_RawTechnicalExportAllowed'].any()


def test_failed_delta_derivative_does_not_inherit_full_rank_qualification(inputs,monkeypatch):
    result,tables,covariance=inputs;original=app._compute_fair_average_value_from_par
    monkeypatch.setattr(app,'_compute_fair_average_value_from_par',
        lambda res,point,facet,level,metric: np.nan if level=='C1' else original(res,point,facet,level,metric))
    frame=app.add_gpcm_fair_average_delta_se(tables,result,covariance)['Criterion']
    assert list(frame.FairM_SE_QualificationStatus)==['TECHNICAL_ONLY','WITHHELD']
    assert list(frame.FairM_SE_RawTechnicalExportAllowed)==[True,False]
    assert frame.FairMSE.iloc[0]>0 and np.isnan(frame.FairMSE.iloc[1])


def test_standalone_csv_retains_status_reason_and_no_public_conclusion(inputs):
    result,tables,covariance=inputs;output=app.add_gpcm_fair_average_delta_se(tables,result,covariance)
    combined=pd.concat([frame.assign(Facet=facet) for facet,frame in output.items()],ignore_index=True)
    restored=pd.read_csv(StringIO(app.to_csv_bytes(combined).decode('utf-8-sig')))
    for metric in ['FairM','FairZ']:
        assert set(restored[metric+'_SE_QualificationStatus'])=={'TECHNICAL_ONLY','WITHHELD'}
        assert restored[metric+'_SE_ReasonCode'].notna().all()
        assert not restored[metric+'_SE_PublicConclusionAllowed'].any()


def test_old_v2_annotation_cache_is_not_reused_after_qualification_change():
    result,tables=identity_fixture();current=app._facets_fair_se_cache_key(result,.95,tables)
    assert current is not None
    old_v2=('_fair_se_facets_tables_cache_v2',*current[1:]);calls=[]
    ui=SimpleNamespace(session_state={old_v2:{'Criterion':pd.DataFrame({'legacy':[True]})}},
        spinner=lambda *a:nullcontext(),caption=lambda *a:None)
    def annotate(*args,**kwargs):calls.append(True);return {'Criterion':pd.DataFrame({'current':[True]})}
    scope=dict(fair_se_enabled=True,fit_supports_fair_se=True,result=result,report_tables=tables,display_tables=tables,
        ci_level=.95,st=ui,t=lambda s:s,_standalone_ui_text=lambda **text:text['en'],
        _facets_fair_se_cache_key=app._facets_fair_se_cache_key,add_gpcm_fair_average_delta_se=annotate)
    code=ui_cache_block();exec(code,scope);exec(code,scope)
    assert calls==[True] and 'current' in scope['display_tables']['Criterion']


@pytest.mark.parametrize('lang',['en','ja'])
def test_current_caption_describes_observed_information_and_unqualified_coverage(lang,monkeypatch):
    result,tables=identity_fixture();captions=[]
    ui=SimpleNamespace(session_state={'lang':lang},spinner=lambda *a:nullcontext(),caption=captions.append)
    monkeypatch.setattr(app,'st',ui)
    scope=dict(fair_se_enabled=True,fit_supports_fair_se=True,result=result,report_tables=tables,display_tables=tables,
        ci_level=.95,st=ui,t=lambda s:s,_standalone_ui_text=app._standalone_ui_text,
        _facets_fair_se_cache_key=app._facets_fair_se_cache_key,add_gpcm_fair_average_delta_se=lambda *a,**k:tables)
    exec(ui_cache_block(),scope)
    assert ('observed information' if lang=='en' else '観測情報') in captions[0]
    assert ('unqualified' if lang=='en' else '未適格') in captions[0]
    assert ('regularization' if lang=='en' else '正則化') in captions[0]
    assert 'Fisher' not in captions[0] and 'MLE' not in captions[0]
