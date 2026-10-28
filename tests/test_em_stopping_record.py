"""Persist the EM condition actually tested, without implying numerical qualification."""
from io import BytesIO
import pickle

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import OptimizeResult
from streamlit.testing.v1 import AppTest

import streamlit_app as app
from mfrm_app import engine
from tests.test_mml_termination_qualification import small_problem


@pytest.mark.parametrize('free_sd',[False,True])
def test_em_records_checked_change_not_the_duplicate_final_trace_entry(free_sd):
    point,idx,cfg,sizes,quad=small_problem()
    cfg['estimate_population_sd']=free_sd
    fit=engine.mfrm_em_mml(point,idx,cfg,sizes,quad,maxit=200,reltol=1e-6)
    assert fit.success
    assert fit.ll_trace[-1]==fit.ll_trace[-2]
    expected=abs(fit.ll_trace[-2]-fit.ll_trace[-3])/(abs(fit.ll_trace[-3])+1e-10)
    assert 0 < expected < 1e-6
    assert fit.em_stopping_relative_loglik_change==expected
    assert fit.em_m_step_options==dict(method='L-BFGS-B',maxiter=50,gtol=1e-5,ftol=1e-9)
    table=engine.build_convergence_summary(fit,cfg,-fit.raw_fun)
    result=dict(opt=fit,config=cfg,convergence=table)
    restored=pickle.loads(pickle.dumps(result))
    qualified=app.qualified_convergence_table(restored)
    saved=pd.read_csv(BytesIO(app.to_csv_bytes(qualified)),float_precision='round_trip')
    assert saved.EMRelativeLogLikTolerance.iloc[0]==1e-6
    assert saved.EMStoppingRelativeLogLikChange.iloc[0]==expected
    assert saved.EMMStepFunctionTolerance.iloc[0]==1e-9
    assert saved.EMMStepGradientTolerance.iloc[0]==1e-5
    assert not saved.QuadratureAccuracyQualified.iloc[0]
    assert not saved.ScientificInferenceReady.iloc[0]


def test_unfinished_em_does_not_invent_a_successful_stopping_change():
    point,idx,cfg,sizes,quad=small_problem()
    fit=engine.mfrm_em_mml(point,idx,cfg,sizes,quad,maxit=1)
    assert not fit.success and np.isnan(fit.em_stopping_relative_loglik_change)
    table=engine.build_convergence_summary(fit,cfg,-fit.raw_fun)
    assert np.isnan(table.EMStoppingRelativeLogLikChange.iloc[0])
    assert table.EMMStepMaxIterations.iloc[0]==50


def test_old_or_non_em_results_do_not_acquire_inferred_controls():
    old=OptimizeResult(success=True,fun=2.,message='retained',nit=3)
    table=engine.build_convergence_summary(old,dict(method='MML',reltol=1e-6),-2.)
    assert not any(column.startswith('EMMStep') or column.startswith('EMRelative')
                   or column.startswith('EMStopping') for column in table.columns)


@pytest.mark.parametrize('language',['ja','en'])
def test_existing_detail_view_displays_saved_controls_without_estimation(monkeypatch,language):
    point,idx,cfg,sizes,quad=small_problem()
    fit=engine.mfrm_em_mml(point,idx,cfg,sizes,quad,maxit=200)
    table=engine.build_convergence_summary(fit,cfg,-fit.raw_fun)
    def view():
        import streamlit as st
        import streamlit_app as app
        app.show_convergence_section(st.session_state['result'],allow_check=False)
    at=AppTest.from_function(view)
    at.session_state['lang']=language
    at.session_state['result']=dict(config=cfg,opt=fit,convergence=table)
    monkeypatch.setattr(engine,'mfrm_em_mml',lambda *a,**k:pytest.fail('No estimation on display'))
    at.run()
    assert not at.exception
    detail=next(d.value for d in at.dataframe if d.value.iloc[:,0].eq('EMRelativeLogLikTolerance').any())
    assert detail.iloc[:,0].eq('EMStoppingRelativeLogLikChange').any()
    expected=str(fit.em_stopping_relative_loglik_change)
    assert detail.loc[detail.iloc[:,0].eq('EMStoppingRelativeLogLikChange')].iloc[0,1]==expected
