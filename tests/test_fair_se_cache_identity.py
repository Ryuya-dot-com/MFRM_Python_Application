"""A parameter vector alone cannot identify likelihood curvature or a report."""
import ast
from contextlib import nullcontext
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


def fixture():
    result=dict(opt=SimpleNamespace(x=np.array([.1,.2]),success=True),
        config=dict(model='GPCM',method='MML',quad_points=31,population_prior_sd=1.,
                    estimate_population_sd=False,facet_names=['Criterion'],
                    input_data_fingerprint='same-upload',analysis_config_fingerprint='same-stale-config'),
        prep=dict(data=pd.DataFrame(dict(Person=['P0','P1'],Criterion=['C0','C0'],Score=[0,1],Weight=[1.,1.]))),
        params=dict(slopes=np.array([1.,1.])),
        facets=dict(person=pd.DataFrame(dict(Level=['P0','P1'],Estimate=[-.3,.4]))))
    tables=dict(Criterion=pd.DataFrame(dict(Level=['C0','C1'],FairM=[.8,1.2],FairZ=[1.,1.],ModelSE=[.2,.3])))
    return result,tables


@pytest.mark.parametrize('mutation',['response','weight','person','prior','free_sd','quadrature','constraint','params','EAP','point'])
def test_same_optimum_does_not_authorize_different_input_or_calibration(mutation):
    result,tables=fixture();changed=deepcopy(result)
    if mutation=='response':changed['prep']['data'].loc[0,'Score']=1
    elif mutation=='weight':changed['prep']['data'].loc[0,'Weight']=2.
    elif mutation=='person':changed['prep']['data'].loc[1,'Person']='P0'
    elif mutation=='prior':changed['config']['population_prior_sd']=1.2
    elif mutation=='free_sd':changed['config']['estimate_population_sd']=True
    elif mutation=='quadrature':changed['config']['quad_points']=61
    elif mutation=='constraint':changed['config']['noncenter_facet']='Criterion'
    elif mutation=='params':changed['params']['slopes'][0]=1.1
    elif mutation=='EAP':changed['facets']['person'].loc[0,'Estimate']+=.1
    else:changed['opt'].x[0]+=.1
    assert app._facets_fair_se_cache_key(result,.95,tables)!=app._facets_fair_se_cache_key(changed,.95,tables)


@pytest.mark.parametrize('mutation',['value','diagnostic','label','order','column_type'])
def test_all_displayed_and_downloaded_table_values_own_the_cache(mutation):
    result,tables=fixture();changed=deepcopy(tables)
    if mutation=='value':changed['Criterion'].loc[0,'FairM']+=.1
    elif mutation=='diagnostic':changed['Criterion'].loc[0,'ModelSE']+=.1
    elif mutation=='label':changed['Criterion'].loc[0,'Level']='C2'
    elif mutation=='order':changed['Criterion']=changed['Criterion'].iloc[::-1]
    else:changed['Criterion']=changed['Criterion'].rename(columns={'FairM':1})
    assert app._facets_fair_se_cache_key(result,.95,tables)!=app._facets_fair_se_cache_key(result,.95,changed)


def test_identical_copy_reuses_key_and_ci_change_does_not():
    result,tables=fixture();key=app._facets_fair_se_cache_key(result,.95,tables)
    assert key is not None
    assert key==app._facets_fair_se_cache_key(deepcopy(result),.95,deepcopy(tables))
    assert key!=app._facets_fair_se_cache_key(result,.9,tables)
    assert key!=('_fair_se_facets_tables_cache',hash(result['opt'].x.tobytes()),.95)


def test_legacy_dictionary_optimizer_coordinates_are_bound_explicitly():
    result,tables=fixture();result['opt']={'x':np.array([.1,.2])}
    before=app._facets_fair_se_cache_key(result,.95,tables)
    result['opt']['x'][0]+=.1
    assert before is not None and before!=app._facets_fair_se_cache_key(result,.95,tables)


@pytest.mark.parametrize('mutation',['no_data','empty_data','no_point','nonfinite_point','unknown_type','no_tables','invalid_ci'])
def test_incomplete_or_unsupported_identity_disables_reuse(mutation):
    result,tables=fixture();ci=.95
    if mutation=='no_data':del result['prep']['data']
    elif mutation=='empty_data':result['prep']['data']=pd.DataFrame()
    elif mutation=='no_point':del result['opt']
    elif mutation=='nonfinite_point':result['opt'].x[0]=np.nan
    elif mutation=='unknown_type':result['config']['unsupported']=object()
    elif mutation=='no_tables':tables=None
    else:ci=np.inf
    assert app._facets_fair_se_cache_key(result,ci,tables) is None


def ui_cache_block():
    tree=ast.parse(Path(app.__file__).read_text())
    node=next(n for n in ast.walk(tree) if isinstance(n,ast.If)
        and ast.unparse(n.test)=='fair_se_enabled and fit_supports_fair_se and report_tables')
    return compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),app.__file__,'exec')


def test_actual_ui_block_reuses_only_current_tables_and_computes_without_a_cache_key():
    result,tables=fixture();calls=[]
    def annotate(raw,res,ci_level):
        calls.append((raw,res,ci_level));out=deepcopy(raw)
        out['Criterion']['FairMSE']=float(res['prep'].get('data',pd.DataFrame()).shape[0])
        return out
    ui=SimpleNamespace(session_state={},spinner=lambda *a:nullcontext(),caption=lambda *a:None)
    old_key=('_fair_se_facets_tables_cache',hash(result['opt'].x.tobytes()),.95)
    ui.session_state[old_key]={'Criterion':pd.DataFrame({'wrong':[True]})}
    scope=dict(fair_se_enabled=True,fit_supports_fair_se=True,result=result,report_tables=tables,
        display_tables=tables,ci_level=.95,st=ui,t=lambda s:s,
        _standalone_ui_text=lambda **text:text['en'],
        _facets_fair_se_cache_key=app._facets_fair_se_cache_key,add_gpcm_fair_average_delta_se=annotate)
    code=ui_cache_block()
    exec(code,scope);exec(code,scope)
    assert len(calls)==1 and list(scope['display_tables']['Criterion'].FairMSE)==[2.,2.]
    scope['result']=deepcopy(result)
    scope['result']['prep']['data']=pd.concat([result['prep']['data']]*2,ignore_index=True)
    exec(code,scope)
    assert len(calls)==2 and list(scope['display_tables']['Criterion'].FairMSE)==[4.,4.]
    scope['report_tables']=deepcopy(tables);scope['report_tables']['Criterion'].loc[0,'FairM']=.9
    exec(code,scope)
    assert len(calls)==3 and scope['display_tables']['Criterion'].FairM.iloc[0]==.9
    del scope['result']['prep']['data']
    exec(code,scope);exec(code,scope)
    assert len(calls)==5 and list(scope['display_tables']['Criterion'].FairMSE)==[0.,0.]
    assert 'wrong' not in scope['display_tables']['Criterion']
