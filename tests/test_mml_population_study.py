"""Pre-registration execution fault controls; all generated data are fixtures."""
import copy
import sqlite3

import numpy as np
import pytest

from validation.mfrm_recovery_metrics import calibration_se_summary
from validation.mml_population_protocol import make_plan, planned_runs
from validation.mml_population_study import AttemptLedger, execute_registered_run, prepare_retained_problem


def fixture_plan(*, negative=None, model='PCM'):
    plan = make_plan()
    if negative:
        scenario = next(s for s in plan['scenarios'] if s['id']=='NEG_'+model+'_'+negative)
    else:
        scenario = next(s for s in plan['scenarios'] if s['id'].startswith('PCM_regular_N48_K4'))
        scenario['design']['persons'] = 12
    scenario['id'] = 'DEVELOPMENT_FIXTURE_'+scenario['id']; scenario['replicates'] = 1
    plan.update(study_id='DEVELOPMENT_TEST_ONLY_population', master_seed=412711,
        classification='DEVELOPMENT_ONLY', scenarios=[scenario])
    plan['estimator'].update(quadrature=[11,21], primary_q=11, check_q=21)
    return plan, scenario, next(planned_runs(plan))


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
@pytest.mark.parametrize('mutation',['all_minimum','all_maximum','unused_internal_category','aliased_rater_task'])
def test_negative_control_status_is_mathematical_and_never_calls_optimizer(tmp_path,monkeypatch,model,mutation):
    import validation.mml_grouped_fit as module
    plan,scenario,run=fixture_plan(negative=mutation,model=model)
    ledger=AttemptLedger(tmp_path/'attempts.sqlite',plan_sha256='a'*64)
    def forbidden(*args,**kwargs): raise AssertionError('Certified no-fit control reached numerical optimization')
    monkeypatch.setattr(module,'minimize',forbidden)
    execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    for q in [11,21]:
        assert ledger.get(run['run_id'],f'{model}_Q{q}_finished')['status']==scenario['expected_status']
    ledger.close()


def test_plan_has_the_full_registered_denominator_and_distinct_128_bit_seeds_without_generating_data():
    plan = make_plan(); runs = list(planned_runs(plan))
    assert len(runs) == 102018
    assert len({r['run_id'] for r in runs}) == len(runs)
    assert len({r['seed'] for r in runs}) == len(runs)
    assert sum(r['cohort']=='population_recovery' for r in runs) == 102000
    assert max(r['seed'] for r in runs) > 2**100
    assert len([s for s in plan['scenarios'] if s['paired_rsm_pcm_comparison']]) == 24
    assert plan['monte_carlo']['expected_mcse_at_95'] == pytest.approx(np.sqrt(.95*.05/2000))
    assert not plan['scientific_inference_ready']


def test_sqlite_events_are_immutable_bound_and_detect_payload_changes(tmp_path):
    path = tmp_path/'attempts.sqlite'; ledger = AttemptLedger(path, plan_sha256='a'*64)
    ledger.append('run1','started',dict(seed=2**110+3, value=np.array([1.,np.nan])))
    assert ledger.get('run1','started')['seed'] == 2**110+3
    assert ledger.get('run1','started')['value'][1] == {'nonfinite':'nan'}
    assert ledger.verify()['events'] == 1
    with pytest.raises(sqlite3.IntegrityError): ledger.append('run1','started',{})
    with pytest.raises(sqlite3.IntegrityError): ledger.connection.execute('DELETE FROM event')
    ledger.connection.rollback()
    with pytest.raises(ValueError,match='another frozen plan'): AttemptLedger(path,plan_sha256='b'*64)
    ledger.connection.execute('DROP TRIGGER no_event_update')
    ledger.connection.execute("UPDATE event SET raw_sha256=?", ('0'*64,));ledger.connection.commit()
    with pytest.raises(ValueError,match='payload'): ledger.verify()
    ledger.close()


@pytest.mark.parametrize('negative,expected', [('missing_declared_rater','input_missing_declared_level'),
    ('fractional_weight','input_nonunit_weight')])
def test_invalid_retained_inputs_are_no_fit_outcomes_in_every_planned_lane(tmp_path,monkeypatch,negative,expected):
    import validation.mml_population_study as module
    plan,scenario,run=fixture_plan(negative=negative);ledger=AttemptLedger(tmp_path/'attempts.sqlite',plan_sha256='a'*64)
    def forbidden(*args,**kwargs): raise AssertionError('Invalid input must not reach fitting')
    monkeypatch.setattr(module,'fit_grouped_joint',forbidden)
    execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    for q in [11,21]: assert ledger.get(run['run_id'],f'PCM_Q{q}_finished')['status']==expected
    assert ledger.get(run['run_id'],'terminal')['truth_available']
    assert ledger.get(run['run_id'],'PCM_assessment')['primary_fit_status']==expected
    assert ledger.verify()['events']==8
    ledger.close()


def test_execution_errors_preserve_both_attempts_and_known_truth_without_retries(tmp_path,monkeypatch):
    import validation.mml_population_study as module
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'attempts.sqlite',plan_sha256='a'*64)
    called=[]
    def fail(*args,**kwargs):
        called.append(kwargs['q']);raise RuntimeError('injected fit failure')
    monkeypatch.setattr(module,'fit_grouped_joint',fail)
    execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert called==[11,21,11,21]
    for model in ['RSM','PCM']:
        assert ledger.get(run['run_id'],model+'_assessment')['status']=='ineligible_candidate'
        assert ledger.get(run['run_id'],model+'_Q11_failed')['error']['message']=='injected fit failure'
    assert not ledger.get(run['run_id'],'terminal')['comparison']['available']
    endpoints=ledger.get(run['run_id'],'PCM_assessment')['endpoints']
    assert endpoints['person']['truth'] is not None
    assert all(not row['eligible'] for row in endpoints['calibration'])
    with pytest.raises(ValueError,match='already started'):
        execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert called==[11,21,11,21];ledger.close()


def test_interrupted_attempt_is_not_replayed(tmp_path,monkeypatch):
    import validation.mml_population_study as module
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'attempts.sqlite',plan_sha256='a'*64)
    def interrupt(*args,**kwargs): raise KeyboardInterrupt('injected interruption')
    monkeypatch.setattr(module,'fit_grouped_joint',interrupt)
    with pytest.raises(KeyboardInterrupt): execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert ledger.get(run['run_id'],'RSM_Q11_started') is not None
    assert ledger.get(run['run_id'],'terminal') is None
    with pytest.raises(ValueError,match='already started'):
        execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    ledger.close()


@pytest.mark.parametrize('alteration', ['seed','cohort','scenario'])
def test_workers_cannot_relabel_registered_cases(tmp_path,alteration):
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'attempts.sqlite',plan_sha256='a'*64)
    if alteration=='seed': run['seed']+=1
    elif alteration=='cohort': run['cohort']='negative_control'
    else:
        scenario=copy.deepcopy(scenario);scenario['design']['population_sd']=2.
    with pytest.raises(ValueError,match='frozen plan'):
        execute_registered_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert ledger.verify()['events']==0;ledger.close()


def test_se_ratio_mcse_keeps_error_variance_and_model_variance_dependence():
    error=np.array([-2.,-1.,1.,2.]);se=np.abs(error)
    result=calibration_se_summary(error,se)
    assert result['ratio']==pytest.approx(np.sqrt(3/4))
    assert result['mcse_ratio']<1e-8  # First-order ratio variation cancels here.
    shifted=calibration_se_summary(error+100,se)
    scaled=calibration_se_summary(error*7,se*7)
    assert shifted['ratio']==pytest.approx(result['ratio'])
    assert scaled['ratio']==pytest.approx(result['ratio'])
    assert shifted['mcse_ratio']==pytest.approx(result['mcse_ratio'],abs=1e-8)
    assert scaled['mcse_ratio']==pytest.approx(result['mcse_ratio'],abs=1e-8)
    independent=calibration_se_summary(error,np.ones(4))
    assert independent['mcse_ratio']>0
    missing=calibration_se_summary([1,np.nan],[.1,np.nan])
    assert missing['planned_replicates']==2 and missing['eligible_replicates']==1
    assert missing['ratio'] is None
