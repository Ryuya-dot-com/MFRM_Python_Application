"""Report-local reuse must never make changed inputs or caller edits stale."""
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from tests.test_app_joint_mml import sample, arguments


def test_complete_inputs_copy_isolation_and_scope_expiry():
    calls=[]
    @app._reuse_report_frame
    def summarize(array, table, settings):
        calls.append(1)
        return pd.DataFrame({'Sum':[array.sum()+table['x'].sum()],
            'Type':[str(table['x'].dtype)],'Index':[str(table.index.tolist())],
            'Setting':[settings['model']]})
    a=np.arange(2000.,dtype=float);df=pd.DataFrame({'x':[1,2]});cfg={'model':'PCM'}
    @app._with_report_calculation_cache
    def report():
        original=summarize(a,df,cfg)
        original.loc[0,'Sum']=-999
        assert summarize(a,df,cfg).Sum.iloc[0]==a.sum()+3
        assert len(calls)==1
        a[1234]+=.125
        assert summarize(a,df,cfg).Sum.iloc[0]==a.sum()+3
        df['x']=df.x.astype(float)
        assert summarize(a,df,cfg).Type.iloc[0]=='float64'
        df.index=[8,3]
        assert summarize(a,df,cfg).Index.iloc[0]=='[8, 3]'
        cfg['model']='GPCM'
        assert summarize(a,df,cfg).Setting.iloc[0]=='GPCM'
        assert len(calls)==5
    report()
    assert app._report_calculation_cache.get() is None
    summarize(a,df,cfg);summarize(a,df,cfg)
    assert len(calls)==7


def test_failure_fallback_nested_scope_and_bounded_retention():
    calls=[]
    @app._reuse_report_frame
    def frame(value, opaque=None):
        calls.append(value)
        if value<0:raise ValueError('invalid calculation')
        return pd.DataFrame({'Value':[value]})
    @app._with_report_calculation_cache
    def inner():
        frame(0);frame(0)
    @app._with_report_calculation_cache
    def report():
        inner();frame(0)
        assert calls==[0]
        for i in range(1,12):frame(i)
        assert len(app._report_calculation_cache.get())==8
        frame(0)
        assert calls.count(0)==2
        opaque=lambda:None
        frame(12,opaque);frame(12,opaque)
        assert calls.count(12)==2
        for _ in range(2):
            with pytest.raises(ValueError,match='invalid calculation'):frame(-1)
        assert calls.count(-1)==2
        raise RuntimeError('report failed')
    with pytest.raises(RuntimeError,match='report failed'):report()
    assert app._report_calculation_cache.get() is None
    inner()
    assert calls.count(0)==3


def test_parallel_exports_keep_separate_scopes():
    barrier=Barrier(2)
    @app._reuse_report_frame
    def frame(value):return pd.DataFrame({'Value':[value]})
    @app._with_report_calculation_cache
    def export(value):
        cache=app._report_calculation_cache.get()
        barrier.wait(timeout=10)
        one=frame(value);one.loc[0,'Value']=-100
        two=frame(value)
        barrier.wait(timeout=10)
        return id(cache),two.Value.iloc[0]
    with ThreadPoolExecutor(max_workers=2) as pool:
        one,two=list(pool.map(export,[1,2]))
    assert one[0]!=two[0] and [one[1],two[1]]==[1,2]
    assert app._report_calculation_cache.get() is None


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
def test_actual_category_tables_match_fresh_computation_after_input_changes(model):
    result=app.mfrm_estimate(**arguments(sample(model),model))
    obs=app.compute_obs_table(result)
    @app._with_report_calculation_cache
    def check():
        for step in range(4):
            if step==1:
                key='steps' if model=='RSM' else 'steps_mat'
                result['params'][key][0]+=.02
            elif step==2:
                obs.loc[obs.index[0],'StdSq']+=.5
                result['config']['population_prior_sd']+=.1
            elif step==3:
                result['prep']['rating_min']=1
                result['prep']['rating_max']=4
            for fn,args in [(app.calc_category_stats,(obs,result)),
                            (app.category_curve_diagnostic_scope_table,(result,))]:
                expected=fn.__wrapped__(*args)
                actual=fn(*args)
                pd.testing.assert_frame_equal(actual,expected,check_exact=True)
                if not actual.empty:
                    actual.iloc[0,0]=-999 if pd.api.types.is_numeric_dtype(actual.dtypes.iloc[0]) else 'caller edit'
                pd.testing.assert_frame_equal(fn(*args),expected,check_exact=True)
    check()
