"""No qualification seeds; retained single attempts and independent lanes."""
from dataclasses import asdict
import hashlib

import pytest

from validation.bias_population_protocol import make_bias_plan
from validation.bias_population_analysis import summarize_bias_scenario
from validation import bias_population_study as study
from validation.mfrm_recovery_design import replicate_seed
from validation.mml_population_protocol import design
from validation.mml_population_study import AttemptLedger


def fixture_plan(model='PCM'):
    plan=make_bias_plan();plan.update(classification='DEVELOPMENT_ONLY',study_id='DEVELOPMENT_ONLY_G3_EXECUTOR_TEST',
        master_seed=91443853,planned_population_datasets=1)
    scenario=plan['scenarios'][0];scenario.update(id=model+'_test',design=asdict(design(model,'regular',persons=12,categories=2)),
        replicates=1,independent_reference_replicates=[1])
    plan['scenarios']=[scenario];plan['estimator'].update(quadrature=[31,61],primary_q=31,check_q=61)
    run=dict(run_id=scenario['id']+'__0001',scenario_id=scenario['id'],replicate=1,
        seed=replicate_seed(plan['master_seed'],plan['study_id']+':'+scenario['id'],1))
    return plan,scenario,run


def test_draft_refuses_generation_even_with_a_supplied_verifier(tmp_path):
    plan=make_bias_plan();scenario=plan['scenarios'][0]
    ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    with pytest.raises(ValueError,match='does not authorize'):
        study.execute_bias_run(ledger,{},scenario,plan,verify_environment=lambda:None)
    assert ledger.verify()['events']==0;ledger.close()


@pytest.mark.parametrize('model',['RSM','PCM','GPCM'])
def test_calibration_failure_keeps_known_lane_input_and_every_attempt(tmp_path,monkeypatch,model):
    plan,scenario,run=fixture_plan(model);ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    calls=[]
    def broken_fit(data,problem,*,q):
        assert ledger.get(run['run_id'],'input_retained') is not None
        assert ledger.get(run['run_id'],f'Q{q}_started') is not None
        calls.append(q);raise ArithmeticError('Injected development fit failure')
    monkeypatch.setattr(study,'fit_grouped_joint',broken_fit)
    study.execute_bias_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    retained=ledger.get(run['run_id'],'input_retained')
    assert hashlib.sha256(retained['csv_utf8'].encode()).hexdigest()==retained['csv_sha256']
    assert ledger.get(run['run_id'],'known_calibration_known_person')['declared_cells']==6
    assert ledger.get(run['run_id'],'fitted_calibration_EAP')['status']=='calibration_unavailable'
    assert ledger.get(run['run_id'],'terminal')['status']=='completed'
    assert calls==[31,61]
    summary=summarize_bias_scenario(ledger,scenario,plan)
    assert summary['lanes'][0]['lane_status_counts']=={'computed':1}
    assert summary['lanes'][3]['cells'][0]['rules']['raw']['all_planned_identification_interval']==[0.,1.]
    before=ledger.verify()
    with pytest.raises(ValueError,match='replay is prohibited'):
        study.execute_bias_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert ledger.verify()==before and calls==[31,61];ledger.close()


def test_environment_interruption_preserves_started_input_without_replay(tmp_path):
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    calls=[]
    def verify():
        calls.append(None)
        if len(calls)>1:raise RuntimeError('Injected environment change')
    with pytest.raises(RuntimeError,match='environment change'):
        study.execute_bias_run(ledger,run,scenario,plan,verify_environment=verify)
    assert ledger.get(run['run_id'],'input_retained') is not None
    assert ledger.get(run['run_id'],'terminal') is None
    assert summarize_bias_scenario(ledger,scenario,plan)['status']=='incomplete'
    with pytest.raises(ValueError,match='replay is prohibited'):
        study.execute_bias_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    ledger.close()


def test_generation_failure_is_terminal_and_not_replaced(tmp_path,monkeypatch):
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    def fail(*args):raise RuntimeError('Injected generation failure')
    monkeypatch.setattr(study,'generate_shifted_sample',fail)
    study.execute_bias_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert ledger.get(run['run_id'],'input_retained') is None
    assert ledger.get(run['run_id'],'terminal')['status']=='generation_failed'
    result=summarize_bias_scenario(ledger,scenario,plan)
    assert result['execution_status_counts']=={'generation_failed':1}
    assert all(lane['incomplete_families']==1 for lane in result['lanes'])
    ledger.close()


def test_wrong_seed_is_rejected_before_attempt_start(tmp_path):
    plan,scenario,run=fixture_plan();ledger=AttemptLedger(tmp_path/'events.sqlite',plan_sha256='fixture')
    run['seed']+=1
    with pytest.raises(ValueError,match='Unregistered RunId or seed'):
        study.execute_bias_run(ledger,run,scenario,plan,verify_environment=lambda:None)
    assert ledger.verify()['events']==0;ledger.close()
