"""Reference transport controls using explicitly excluded development seeds."""
import hashlib
import json

import numpy as np
import pytest

from validation.mml_population_protocol import design, make_plan, planned_runs
from validation.mml_population_references import export_reference_bundle
from validation.mml_population_study import AttemptLedger, execute_registered_run
from tests.test_mml_population_study import fixture_plan
from dataclasses import asdict


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('categories',[2,6])
def test_reference_export_preserves_same_points_nodes_and_declared_support(tmp_path,model,categories):
    plan,scenario,_=fixture_plan()
    scenario.update(design=asdict(design(model,'regular',persons=12,categories=categories)),
        primary_model=model,fit_models=[model],paired_rsm_pcm_comparison=False)
    run=next(planned_runs(plan));ledger=AttemptLedger(tmp_path/'a.sqlite',plan_sha256='a'*64)
    execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    out=tmp_path/'export';result=export_reference_bundle(out,ledger,scenario,1,plan,'a'*64)
    expected=5 if model!='GPCM' else 4
    assert result['required_cases']==expected
    for name,sha in json.loads((out/'manifest.json').read_bytes()).items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest()==sha
    cases=json.loads((out/'finite_input.json').read_bytes())['cases']
    native=json.loads((out/'native_rowwise.json').read_bytes())
    outputs=json.loads((out/'python_finite.json').read_bytes())
    assert len(cases)==len(native)==len(outputs)==2
    for case,reference,fit in zip(cases,native,outputs,strict=True):
        assert case['config']['n_cat']==categories
        assert min(case['weights'])>0 and sum(case['weights'])==pytest.approx(1.)
        assert reference['nll']==pytest.approx(fit['nll'],abs=1e-9)
        np.testing.assert_allclose(reference['gradient'],fit['gradient'],rtol=0,atol=1e-9)
        assert len(fit['expanded'])==len(fit['expanded_jacobian'])
    if model!='GPCM':
        raw=ledger.get(run['run_id'],'input_retained')['csv_utf8'].encode()
        csv_files=list((out/'tam_input').glob('*.csv'));assert len(csv_files)==1
        assert csv_files[0].read_bytes()==raw
    ledger.close()


def test_missing_results_remain_in_reference_denominator_and_selection_is_fixed(tmp_path):
    plan,scenario,_=fixture_plan();ledger=AttemptLedger(tmp_path/'a.sqlite',plan_sha256='a'*64)
    result=export_reference_bundle(tmp_path/'missing',ledger,scenario,1,plan,'a'*64)
    assert result['required_cases']==9 and result['available_cases']==0
    reference=json.loads((tmp_path/'missing/reference_plan.json').read_bytes())
    assert reference['missing_required_cases']==9 and reference['seed'] is None
    with pytest.raises(ValueError,match='not selected'):
        export_reference_bundle(tmp_path/'wrong',ledger,scenario,2,plan,'a'*64)
    ledger.close()


def test_reference_interruption_and_failure_never_trigger_numerical_replay(tmp_path):
    from validation.mml_population_reference_execution import _once
    ledger=AttemptLedger(tmp_path/'a.sqlite',plan_sha256='a'*64);calls=[]
    def interrupt():
        calls.append('started');raise KeyboardInterrupt('injected interruption')
    with pytest.raises(KeyboardInterrupt):_once(ledger,'case1','R',interrupt,lambda:None)
    assert _once(ledger,'case1','R',interrupt,lambda:None)['status']=='interrupted_reference_retained'
    assert calls==['started']
    def fail():
        calls.append('failed');raise ValueError('injected failure')
    assert _once(ledger,'case2','R',fail,lambda:None)['status']=='reference_failed'
    assert _once(ledger,'case2','R',fail,lambda:None)['status']=='reference_failed'
    assert calls==['started','failed']
    assert ledger.verify()['events']==4
    ledger.close()


def test_reference_distance_rejects_silent_broadcasting_and_nonfinite_values():
    from validation.mml_population_reference_review import distance
    with pytest.raises(ValueError):distance([1,2],[[1,2]])
    with pytest.raises(ValueError):distance([np.nan],[np.nan])
    assert distance([1,2],[1,3])==1
