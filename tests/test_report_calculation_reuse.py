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


@pytest.mark.parametrize("builder_name", ["build_final_report_readiness", "build_report_ready_summary_panel",
                                        "build_manuscript_claim_guide", "build_claim_to_evidence_matrix"])
@pytest.mark.parametrize("case", ["facet_boundary", "fixed_sd", "free_sd"])
def test_readiness_reuse_preserves_evidence_and_rechecks_changed_fit(monkeypatch, builder_name, case):
    from tests.test_jmle_facet_boundary import separated_data
    data = separated_data() if case == "facet_boundary" else app._make_self_test_gradient_data()
    options = ({"method": "MML", "mml_engine": "EM", "quad_points": 15,
                "estimate_population_sd": case == "free_sd"} if case != "facet_boundary" else {})
    result = app.mfrm_estimate(data, person_col="Person", facet_cols=["Rater", "Task"],
        score_col="Score", maxit=100, **options)
    diagnostics = app.mfrm_diagnostics(result, compute_pca=False)
    builder = getattr(app, builder_name)
    fresh = builder.__wrapped__
    expected = fresh(result, diagnostics)
    identity = app.build_result_analysis_identity
    calls = []

    def count_identity(value):
        calls.append(1)
        return identity(value)

    monkeypatch.setattr(app, "build_result_analysis_identity", count_identity)

    @app._with_report_calculation_cache
    def report():
        first = builder(result, diagnostics)
        pd.testing.assert_frame_equal(first, expected, check_exact=True)
        count_before_reuse = len(calls)
        status_column = next(column for column in ("Status", "ReportStatus", "ClaimStatus") if column in first)
        first.loc[0, status_column] = "caller edit"
        if "evidence_records" in first.attrs:
            first.attrs["evidence_records"][0]["caller_edit"] = True
        again = builder(result, diagnostics)
        pd.testing.assert_frame_equal(again, expected, check_exact=True)
        assert again.attrs == expected.attrs
        assert count_before_reuse > 0 and len(calls) == count_before_reuse
        result["opt"].success = False
        actual = builder(result, diagnostics)
        assert len(calls) > count_before_reuse
        pd.testing.assert_frame_equal(actual, fresh(result, diagnostics), check_exact=True)

    report()
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


def test_agreement_reuse_rechecks_ratings_weights_model_and_selected_facet(monkeypatch):
    result = app.mfrm_estimate(app._make_self_test_gradient_data(), "Person", ["Rater", "Task"], "Score", maxit=100)
    diagnostics = {"obs": app.compute_obs_table(result)}
    builder = app._agreement_summary_for_facet
    original = app.calc_interrater_agreement
    calls = []

    def count(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(app, "calc_interrater_agreement", count)

    @app._with_report_calculation_cache
    def report():
        for step in range(5):
            target = "Task" if step == 1 else "Rater"
            if step == 2:
                diagnostics["obs"].loc[0, "Observed"] += 1
            elif step == 3:
                diagnostics["obs"]["Weight"] = 1.0
                diagnostics["obs"].loc[0, "Weight"] = 2.0
            elif step == 4:
                result["params"]["facets"]["Rater"][0] += .1
            before = len(calls)
            actual = builder(result, diagnostics, ["Rater", "Task"], target)
            assert len(calls) == before + 1
            expected = builder.__wrapped__(result, diagnostics, ["Rater", "Task"], target)
            pd.testing.assert_frame_equal(actual, expected, check_exact=True)
            assert not actual.empty
            actual.loc[actual.index[0], "ExactAgreement"] = -999
            before = len(calls)
            pd.testing.assert_frame_equal(builder(result, diagnostics, ["Rater", "Task"], target), expected, check_exact=True)
            assert len(calls) == before

    report()
    assert app._report_calculation_cache.get() is None
