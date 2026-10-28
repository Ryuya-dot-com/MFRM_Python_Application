"""Run against each prepared source bundle, without generating any study data."""
import importlib.util

import pytest

from validation import mml_population_study as mml
from validation.mml_population_protocol import make_plan, planned_runs


KINDS=['mml']
if importlib.util.find_spec('validation.bias_population_study') is not None:
    KINDS.append('g3')


def executor(kind):
    if kind=='mml':
        plan=make_plan()
        return mml,mml.execute_registered_run,'generate_population_sample',plan,next(planned_runs(plan))
    from validation import bias_population_study as g3
    from validation.bias_population_protocol import make_registered_bias_plan,planned_bias_runs
    plan=make_registered_bias_plan()
    return g3,g3.execute_bias_run,'generate_shifted_sample',plan,next(planned_bias_runs(plan))


@pytest.mark.parametrize('kind',KINDS)
def test_binding_precedes_input_and_no_replay(kind,tmp_path,monkeypatch):
    mod,run_once,generator,plan,run=executor(kind)
    scenario=plan['scenarios'][0]
    ledger=mml.AttemptLedger(tmp_path/'attempt.sqlite',plan_sha256='retained-original-study')
    binding={'schema':'mfrm-run-implementation-v1','study_T0_sha256':'original','amendment_T0_sha256':'candidate'}
    def fail_before_data(*args):
        assert ledger.get(run['run_id'],'execution_implementation')==binding
        raise RuntimeError('injected before generating data')
    monkeypatch.setattr(mod,generator,fail_before_data)
    try:
        run_once(ledger,run,scenario,plan,verify_environment=lambda:binding)
        phases=[row[0] for row in ledger.connection.execute('SELECT phase FROM event ORDER BY sequence')]
        assert phases==['generation_started','execution_implementation','terminal']
        assert ledger.get(run['run_id'],'generation_started')==run
        assert ledger.get(run['run_id'],'terminal')['status']=='generation_failed'
        chain=ledger.verify()
        with pytest.raises(ValueError,match='already started'):
            run_once(ledger,run,scenario,plan,verify_environment=lambda:binding)
        assert ledger.verify()==chain
    finally:ledger.close()


@pytest.mark.parametrize('kind',KINDS)
def test_failed_verification_starts_no_attempt(kind,tmp_path):
    _,run_once,_,plan,run=executor(kind)
    ledger=mml.AttemptLedger(tmp_path/'attempt.sqlite',plan_sha256='retained-original-study')
    def refused():raise ValueError('source changed')
    try:
        with pytest.raises(ValueError,match='source changed'):
            run_once(ledger,run,plan['scenarios'][0],plan,verify_environment=refused)
        assert not ledger.started(run['run_id'])
        assert ledger.verify()['events']==0
    finally:ledger.close()
