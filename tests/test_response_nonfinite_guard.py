"""Infinite input is a domain error; ordinary missing/zero rows keep their policy."""
import numpy as np
import pandas as pd
import pytest

import streamlit_app as app


@pytest.mark.parametrize('column',['Score','Weight'])
@pytest.mark.parametrize('value',[np.inf,-np.inf,'inf','-Infinity'])
def test_infinite_input_is_reported_and_rejected_before_numerical_work(column,value):
    data=pd.DataFrame(dict(Person=['p0','p0','p1','p1'],Rater=['r0','r1','r0','r1'],
                           Score=[0,1,2,1],Weight=[1,1,1,1]),dtype=object)
    data.loc[0,column]=value
    kwargs=dict(person_col='Person',facet_cols=['Rater'],score_col='Score',weight_col='Weight')
    audit=app.build_response_data_audit(data,**kwargs)
    row=audit['rows'].iloc[0]
    assert row.LikelihoodAction=='error_before_fit'
    assert not row.IncludedInLikelihood
    assert row[column+'NonFinite']
    assert row.ExclusionReason=='non-finite '+column.lower()
    assert int(audit['summary'].set_index('Item').loc['Non-finite '+column.lower()+' error','Count'])==1
    with pytest.raises(ValueError,match='must (contain |be )finite'):
        app.prepare_mfrm_data(data,**kwargs,rating_min=0,rating_max=2)


def test_finite_missing_and_zero_weight_rows_preserve_likelihood_membership():
    data=pd.DataFrame(dict(Person=['p0','p0','p1','p1','p2','p3'],Rater=['r0','r1','r0','r1','r0','r1'],
                           Score=[0,1,2,None,1,2],Weight=[1,2,1,1,0,None]))
    kwargs=dict(person_col='Person',facet_cols=['Rater'],score_col='Score',weight_col='Weight')
    audit=app.build_response_data_audit(data,**kwargs)
    assert audit['rows'].IncludedInLikelihood.tolist()==[True,True,True,False,False,False]
    assert not audit['rows'][['ScoreNonFinite','WeightNonFinite']].to_numpy().any()
    prep=app.prepare_mfrm_data(data,**kwargs,rating_min=0,rating_max=2)
    assert prep['data'].Weight.tolist()==[1.,2.,1.]
    assert prep['data'].score_k.tolist()==[0,1,2]


@pytest.mark.parametrize('column',['Score','Weight'])
@pytest.mark.parametrize('value',[np.inf,-np.inf])
@pytest.mark.parametrize('otherwise_excluded',[False,True])
def test_exported_python_rejects_infinite_values_before_optimization(monkeypatch,column,value,otherwise_excluded):
    import scipy.optimize
    data=pd.DataFrame(dict(Person=['p0','p0','p1','p1'],Rater=['r0','r1','r0','r1'],
                           Score=[0.,1.,2.,1.],Weight=[1.,1.,1.,1.]))
    data.loc[0,column]=value
    if otherwise_excluded:
        data.loc[0,'Weight' if column=='Score' else 'Score']=0 if column=='Score' else np.nan
    monkeypatch.setattr(pd,'read_csv',lambda *a,**kw:data.copy())
    def no_fit(*a,**kw):
        pytest.fail('Optimization started with an infinite input value')
    monkeypatch.setattr(scipy.optimize,'minimize',no_fit)
    script=app._generate_repro_python_script(dict(config=dict(model='RSM',method='JMLE',
        facet_names=['Rater'],weight_col='Weight',n_cat=3)))
    with pytest.raises(ValueError,match='must (contain |be )finite'):
        exec(compile(script,'exported_mfrm.py','exec'),{'__name__':'__main__'})
